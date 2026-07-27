#!/usr/bin/env python
"""Run a small, deliberately non-singular axisymmetric control experiment.

The physical field is swirl-only,

    u = u^theta(r,t) e_theta,  u^theta = r u1,

and is independent of z on a periodic cylinder.  Meridional velocity and
omega1 = omega^theta/r are zero.  Radial pressure balances the centripetal
term, while u1 obeys

    partial_t u1 = nu (partial_rr + 3/r partial_r) u1.

The Gaussian analytic solution is evolved with an independent tridiagonal
Crank--Nicolson implementation.  This is a control for the numerical
infrastructure, not a singularity search.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ns_certificate_lab._integrity import strict_json_loads, write_with_digest
from ns_certificate_lab.provenance import collect_runtime_provenance


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Evolution:
    r: np.ndarray
    times: np.ndarray
    profiles: np.ndarray
    peak_u_theta: np.ndarray
    peak_vorticity: np.ndarray
    peak_gradient: np.ndarray
    energy_per_z_length: np.ndarray
    all_finite: bool


def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    """Integrate on NumPy 1.24--1.x and NumPy 2.x without deprecation noise."""

    modern = getattr(np, "trapezoid", None)
    if modern is not None:
        return float(modern(y, x))
    return float(np.trapz(y, x))


def _inside_repository(path: Path) -> bool:
    try:
        path.resolve().relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def analytic_u1(
    r: np.ndarray,
    t: float,
    *,
    viscosity: float,
    amplitude: float,
    sigma: float,
) -> np.ndarray:
    """Four-radial-dimensional Gaussian heat solution."""

    width2 = sigma * sigma + 4.0 * viscosity * t
    scale = (sigma * sigma / width2) ** 2
    return amplitude * scale * np.exp(-(r * r) / width2)


def _operator_coefficients(r: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tridiagonal coefficients for d_rr + 3/r d_r.

    The last radial point is a prescribed Dirichlet boundary and is therefore
    excluded from the unknown vector.  At the axis, even reflection gives
    L f(0) = 4 f_rr(0) = 8 (f_1-f_0)/dr^2 + O(dr^2).
    """

    h = float(r[1] - r[0])
    count = r.size - 1
    lower = np.zeros(count, dtype=np.float64)
    diagonal = np.full(count, -2.0 / (h * h), dtype=np.float64)
    upper = np.zeros(count, dtype=np.float64)
    diagonal[0] = -8.0 / (h * h)
    upper[0] = 8.0 / (h * h)
    indices = np.arange(1, count)
    radii = r[indices]
    lower[indices] = 1.0 / (h * h) - 3.0 / (2.0 * radii * h)
    upper[indices] = 1.0 / (h * h) + 3.0 / (2.0 * radii * h)
    return lower, diagonal, upper


def _apply_operator(
    values: np.ndarray,
    boundary_value: float,
    lower: np.ndarray,
    diagonal: np.ndarray,
    upper: np.ndarray,
) -> np.ndarray:
    result = diagonal * values
    result[1:] += lower[1:] * values[:-1]
    result[:-1] += upper[:-1] * values[1:]
    result[-1] += upper[-1] * boundary_value
    return result


def _solve_tridiagonal(
    lower: np.ndarray,
    diagonal: np.ndarray,
    upper: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    """Thomas solve, kept independent of the package differentiation code."""

    count = diagonal.size
    c_prime = np.zeros(count - 1, dtype=np.float64)
    d_prime = np.empty(count, dtype=np.float64)
    pivot = float(diagonal[0])
    if not np.isfinite(pivot) or abs(pivot) < 1e-15:
        raise ArithmeticError("invalid first tridiagonal pivot")
    c_prime[0] = upper[0] / pivot
    d_prime[0] = rhs[0] / pivot
    for index in range(1, count):
        pivot = diagonal[index] - lower[index] * c_prime[index - 1]
        if not np.isfinite(pivot) or abs(pivot) < 1e-15:
            raise ArithmeticError(f"invalid tridiagonal pivot at {index}")
        if index < count - 1:
            c_prime[index] = upper[index] / pivot
        d_prime[index] = (
            rhs[index] - lower[index] * d_prime[index - 1]
        ) / pivot
    solution = np.empty(count, dtype=np.float64)
    solution[-1] = d_prime[-1]
    for index in range(count - 2, -1, -1):
        solution[index] = d_prime[index] - c_prime[index] * solution[index + 1]
    return solution


def _diagnostics(r: np.ndarray, u1: np.ndarray) -> tuple[float, float, float, float]:
    u_theta = r * u1
    edge_order = 2 if r.size >= 3 else 1
    d_u_theta = np.gradient(u_theta, r, edge_order=edge_order)
    d_u1 = np.gradient(u1, r, edge_order=edge_order)
    omega_z = 2.0 * u1 + r * d_u1
    energy = math.pi * _trapezoid((u1 * u1) * r**3, r)
    return (
        float(np.max(np.abs(u_theta))),
        float(np.max(np.abs(omega_z))),
        float(np.max(np.abs(d_u_theta))),
        energy,
    )


def evolve(
    *,
    nr: int,
    r_max: float,
    final_time: float,
    viscosity: float,
    amplitude: float,
    sigma: float,
    dt_over_dr: float,
    sample_count: int,
    exact_outer_boundary: bool,
) -> Evolution:
    if nr < 9:
        raise ValueError("nr must be at least 9")
    r = np.linspace(0.0, r_max, nr, dtype=np.float64)
    dr = float(r[1] - r[0])
    requested_dt = dt_over_dr * dr
    step_count = max(1, int(math.ceil(final_time / requested_dt)))
    dt = final_time / step_count
    times_all = np.linspace(0.0, final_time, step_count + 1)
    sample_indices = np.unique(
        np.rint(np.linspace(0, step_count, sample_count)).astype(int)
    )

    lower_l, diagonal_l, upper_l = _operator_coefficients(r)
    theta = 0.5
    matrix_lower = -theta * dt * viscosity * lower_l
    matrix_diagonal = 1.0 - theta * dt * viscosity * diagonal_l
    matrix_upper = -theta * dt * viscosity * upper_l[:-1]

    current = analytic_u1(
        r,
        0.0,
        viscosity=viscosity,
        amplitude=amplitude,
        sigma=sigma,
    )
    if not exact_outer_boundary:
        current[-1] = 0.0

    saved_times: list[float] = []
    saved_profiles: list[np.ndarray] = []
    peak_u_theta: list[float] = []
    peak_vorticity: list[float] = []
    peak_gradient: list[float] = []
    energies: list[float] = []
    sample_set = set(int(value) for value in sample_indices)

    def save(time: float, profile: np.ndarray) -> None:
        peak_u, peak_omega, peak_grad, energy = _diagnostics(r, profile)
        saved_times.append(time)
        saved_profiles.append(profile.copy())
        peak_u_theta.append(peak_u)
        peak_vorticity.append(peak_omega)
        peak_gradient.append(peak_grad)
        energies.append(energy)

    save(0.0, current)
    for step in range(step_count):
        old_time = float(times_all[step])
        new_time = float(times_all[step + 1])
        if exact_outer_boundary:
            old_boundary = float(
                analytic_u1(
                    r[-1:],
                    old_time,
                    viscosity=viscosity,
                    amplitude=amplitude,
                    sigma=sigma,
                )[0]
            )
            new_boundary = float(
                analytic_u1(
                    r[-1:],
                    new_time,
                    viscosity=viscosity,
                    amplitude=amplitude,
                    sigma=sigma,
                )[0]
            )
        else:
            old_boundary = 0.0
            new_boundary = 0.0
        unknown = current[:-1]
        rhs = unknown + (1.0 - theta) * dt * viscosity * _apply_operator(
            unknown,
            old_boundary,
            lower_l,
            diagonal_l,
            upper_l,
        )
        rhs[-1] += theta * dt * viscosity * upper_l[-1] * new_boundary
        next_unknown = _solve_tridiagonal(
            matrix_lower,
            matrix_diagonal,
            matrix_upper,
            rhs,
        )
        current = np.concatenate((next_unknown, np.array([new_boundary])))
        if step + 1 in sample_set:
            save(new_time, current)

    arrays = [
        np.asarray(saved_times),
        np.asarray(saved_profiles),
        np.asarray(peak_u_theta),
        np.asarray(peak_vorticity),
        np.asarray(peak_gradient),
        np.asarray(energies),
    ]
    return Evolution(
        r=r,
        times=arrays[0],
        profiles=arrays[1],
        peak_u_theta=arrays[2],
        peak_vorticity=arrays[3],
        peak_gradient=arrays[4],
        energy_per_z_length=arrays[5],
        all_finite=bool(all(np.all(np.isfinite(array)) for array in arrays)),
    )


def _observed_orders(errors: list[float], resolutions: list[int]) -> list[float]:
    orders: list[float] = []
    for coarse, fine, n_coarse, n_fine in zip(
        errors[:-1], errors[1:], resolutions[:-1], resolutions[1:]
    ):
        h_ratio = (n_fine - 1) / (n_coarse - 1)
        orders.append(float(math.log(coarse / fine) / math.log(h_ratio)))
    return orders


def _relative_l2(r: np.ndarray, error: np.ndarray, reference: np.ndarray) -> float:
    numerator = _trapezoid(error * error * r**3, r)
    denominator = _trapezoid(reference * reference * r**3, r)
    return math.sqrt(numerator / denominator)


def _blowup_guard(evolution: Evolution) -> dict[str, Any]:
    """Reject a blow-up interpretation unless sustained resolved growth exists."""

    values = evolution.peak_vorticity
    midpoint = max(1, values.size // 2)
    tail = values[midpoint:]
    monotone_growth = bool(
        tail.size >= 4
        and np.all(np.diff(tail) >= 0.0)
        and tail[-1] > 1.05 * tail[0]
    )
    return {
        "classification": (
            "candidate_requires_multiresolution_fit"
            if monotone_growth
            else "bounded_or_decaying_control"
        ),
        "sustained_tail_growth": monotone_growth,
        "reciprocal_blowup_fit_attempted": False,
        "reason": (
            "growth gate passed; a separate preregistered fit would be required"
            if monotone_growth
            else "peak physical vorticity lacks sustained monotone growth"
        ),
        "important_limitation": (
            "This heuristic can reject an obvious false positive; it cannot "
            "prove regularity or singularity."
        ),
    }


def _write_csv(path: Path, header: list[str], rows: list[list[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _polyline(
    x: np.ndarray,
    y: np.ndarray,
    *,
    left: float,
    top: float,
    width: float,
    height: float,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
) -> str:
    x_den = max(x_range[1] - x_range[0], np.finfo(float).eps)
    y_den = max(y_range[1] - y_range[0], np.finfo(float).eps)
    points = []
    for x_value, y_value in zip(x, y):
        px = left + width * (float(x_value) - x_range[0]) / x_den
        py = top + height * (1.0 - (float(y_value) - y_range[0]) / y_den)
        points.append(f"{px:.3f},{py:.3f}")
    return " ".join(points)


def _write_svg(
    path: Path,
    *,
    title: str,
    x_label: str,
    y_label: str,
    series: list[tuple[str, np.ndarray, np.ndarray, str]],
) -> None:
    width, height = 800.0, 500.0
    left, top, plot_width, plot_height = 90.0, 55.0, 660.0, 365.0
    all_x = np.concatenate([entry[1] for entry in series])
    all_y = np.concatenate([entry[2] for entry in series])
    x_range = (float(np.min(all_x)), float(np.max(all_x)))
    y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
    padding = 0.05 * max(y_max - y_min, 1e-12)
    y_range = (y_min - padding, y_max + padding)
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" '
            f'height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">'
        ),
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2:.1f}" y="28" text-anchor="middle" '
        f'font-family="sans-serif" font-size="18">{title}</text>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" '
        f'y2="{top + plot_height}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" '
        'stroke="black"/>',
        f'<text x="{left + plot_width/2}" y="475" text-anchor="middle" '
        f'font-family="sans-serif" font-size="14">{x_label}</text>',
        f'<text x="20" y="{top + plot_height/2}" text-anchor="middle" '
        f'transform="rotate(-90 20 {top + plot_height/2})" '
        f'font-family="sans-serif" font-size="14">{y_label}</text>',
    ]
    for index, (label, x_values, y_values, colour) in enumerate(series):
        points = _polyline(
            x_values,
            y_values,
            left=left,
            top=top,
            width=plot_width,
            height=plot_height,
            x_range=x_range,
            y_range=y_range,
        )
        lines.append(
            f'<polyline points="{points}" fill="none" stroke="{colour}" '
            'stroke-width="2"/>'
        )
        legend_y = 70 + 22 * index
        lines.append(
            f'<line x1="590" y1="{legend_y}" x2="620" y2="{legend_y}" '
            f'stroke="{colour}" stroke-width="2"/>'
        )
        lines.append(
            f'<text x="628" y="{legend_y + 5}" font-family="sans-serif" '
            f'font-size="12">{label}</text>'
        )
    lines.append(
        f'<text x="{left}" y="445" font-family="monospace" font-size="11">'
        f'x=[{x_range[0]:.4g}, {x_range[1]:.4g}], '
        f'y=[{y_range[0]:.4g}, {y_range[1]:.4g}]</text>'
    )
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    # Capture once before generated files can change the checkout's dirty bit.
    run_provenance = collect_runtime_provenance()
    np.random.seed(int(config["seed"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite non-empty evidence directory: {output_dir}"
        )

    viscosity = float(config["viscosity"])
    amplitude = float(config["amplitude"])
    sigma = float(config["sigma"])
    final_time = float(config["final_time"])
    r_max = float(config["r_max"])
    resolutions = [int(value) for value in config["resolutions"]]
    dt_over_dr = float(config["dt_over_dr"])
    sample_count = int(config["sample_count"])
    normalization = config["normalization"]
    length_scale = float(normalization["length_scale_L0"])
    velocity_scale = float(normalization["velocity_scale_U0"])
    time_scale = length_scale / velocity_scale
    reynolds_number = velocity_scale * length_scale / viscosity

    evolutions: list[Evolution] = []
    relative_errors: list[float] = []
    max_errors: list[float] = []
    for resolution in resolutions:
        evolution = evolve(
            nr=resolution,
            r_max=r_max,
            final_time=final_time,
            viscosity=viscosity,
            amplitude=amplitude,
            sigma=sigma,
            dt_over_dr=dt_over_dr,
            sample_count=sample_count,
            exact_outer_boundary=True,
        )
        exact = analytic_u1(
            evolution.r,
            final_time,
            viscosity=viscosity,
            amplitude=amplitude,
            sigma=sigma,
        )
        error = evolution.profiles[-1] - exact
        relative_errors.append(_relative_l2(evolution.r, error, exact))
        max_errors.append(float(np.max(np.abs(error))))
        evolutions.append(evolution)

    orders = _observed_orders(relative_errors, resolutions)
    boundary_config = config["boundary_sensitivity"]
    boundary_evolutions: list[Evolution] = []
    for radius in boundary_config["radii"]:
        nr = int(round(float(radius) * int(boundary_config["points_per_unit"]))) + 1
        boundary_evolutions.append(
            evolve(
                nr=nr,
                r_max=float(radius),
                final_time=final_time,
                viscosity=viscosity,
                amplitude=amplitude,
                sigma=sigma,
                dt_over_dr=dt_over_dr,
                sample_count=sample_count,
                exact_outer_boundary=False,
            )
        )
    small, large = boundary_evolutions
    mask = small.r <= float(boundary_config["compare_radius"]) + 1e-14
    large_on_small = np.interp(
        small.r[mask],
        large.r,
        large.profiles[-1],
    )
    boundary_difference = float(
        np.max(np.abs(small.profiles[-1, mask] - large_on_small))
    )

    finest = evolutions[-1]
    energy_differences = np.diff(finest.energy_per_z_length)
    energy_scale = max(float(finest.energy_per_z_length[0]), 1e-300)
    maximum_relative_energy_increase = float(
        max(0.0, np.max(energy_differences, initial=0.0) / energy_scale)
    )
    initial_peak = max(float(finest.peak_vorticity[0]), 1e-300)
    peak_factor = float(np.max(finest.peak_vorticity) / initial_peak)
    acceptance = config["acceptance"]
    checks = {
        "all_values_finite": bool(all(item.all_finite for item in evolutions)),
        "errors_decrease_under_refinement": bool(
            all(fine < coarse for coarse, fine in zip(relative_errors, relative_errors[1:]))
        ),
        "observed_order_passes": bool(
            orders and min(orders) >= float(acceptance["minimum_observed_order"])
        ),
        "finest_error_passes": bool(
            relative_errors[-1]
            <= float(acceptance["maximum_finest_relative_l2_error"])
        ),
        "energy_nonincreasing_with_tolerance": bool(
            maximum_relative_energy_increase
            <= float(acceptance["energy_relative_increase_tolerance"])
        ),
        "stable_peak": bool(
            peak_factor <= float(acceptance["stability_peak_factor"])
        ),
        "outer_boundary_not_dominant": bool(
            boundary_difference
            <= float(boundary_config["maximum_allowed_difference"])
        ),
    }

    config_snapshot = output_dir / "config.snapshot.json"
    config_snapshot.write_text(
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _write_csv(
        output_dir / "convergence.csv",
        ["nr", "dr", "relative_l2_error", "max_error", "order_from_previous"],
        [
            [
                resolution,
                r_max / (resolution - 1),
                relative_errors[index],
                max_errors[index],
                "" if index == 0 else orders[index - 1],
            ]
            for index, resolution in enumerate(resolutions)
        ],
    )
    _write_csv(
        output_dir / "diagnostics.csv",
        [
            "time",
            "peak_u_theta",
            "peak_vorticity",
            "peak_radial_gradient_u_theta",
            "energy_per_unit_z",
        ],
        [
            [
                finest.times[index],
                finest.peak_u_theta[index],
                finest.peak_vorticity[index],
                finest.peak_gradient[index],
                finest.energy_per_z_length[index],
            ]
            for index in range(finest.times.size)
        ],
    )
    np.savez_compressed(
        output_dir / "profiles.npz",
        r=finest.r,
        times=finest.times,
        u1=finest.profiles,
    )
    _write_svg(
        output_dir / "energy.svg",
        title="Non-singular control: energy decay",
        x_label="physical time",
        y_label="kinetic energy per unit z length",
        series=[
            (
                "Crank-Nicolson",
                finest.times,
                finest.energy_per_z_length,
                "#1f77b4",
            )
        ],
    )
    profile_colours = ["#444444", "#2ca02c", "#d62728"]
    profile_indices = [0, finest.times.size // 2, finest.times.size - 1]
    _write_svg(
        output_dir / "profiles.svg",
        title="Smooth swirl profile u1(r,t)",
        x_label="r",
        y_label="u1",
        series=[
            (
                f"t={finest.times[index]:.3f}",
                finest.r,
                finest.profiles[index],
                profile_colours[position],
            )
            for position, index in enumerate(profile_indices)
        ],
    )

    summary: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Non-singular numerical control only. This output is neither a "
            "singularity candidate nor a mathematical proof."
        ),
        "method": {
            "equation": "u1_t = nu * (u1_rr + 3*u1_r/r)",
            "domain": (
                f"0<=r<={r_max}, z periodic with period "
                f"{float(config['z_period'])}; fields are z-independent"
            ),
            "time_integrator": "Crank-Nicolson, theta=1/2",
            "spatial_operator": "independent second-order centered radial stencil",
            "axis_limit": "L u1(0) = 4*u1_rr(0)",
            "outer_boundary_for_convergence": "analytic Dirichlet trace",
        },
        "physical_interpretation": {
            "coordinate_system": config["coordinate_system"],
            "units": config["units"],
            "normalization": {
                "length_scale_L0": length_scale,
                "velocity_scale_U0": velocity_scale,
                "time_scale_L0_over_U0": time_scale,
                "reynolds_number_U0_L0_over_nu": reynolds_number,
            },
        },
        "reproducibility": {
            "seed": int(config["seed"]),
            "config_sha256": _sha256(config_snapshot),
            "runtime_provenance": run_provenance,
        },
        "convergence": {
            "resolutions": resolutions,
            "relative_l2_errors": relative_errors,
            "max_errors": max_errors,
            "observed_orders": orders,
        },
        "stability_and_physics": {
            "all_finite": finest.all_finite,
            "initial_energy_per_unit_z": float(finest.energy_per_z_length[0]),
            "final_energy_per_unit_z": float(finest.energy_per_z_length[-1]),
            "maximum_relative_energy_increase": maximum_relative_energy_increase,
            "peak_vorticity_factor": peak_factor,
            "boundary_inner_max_difference": boundary_difference,
            "boundary_radii": boundary_config["radii"],
            "boundary_compare_radius": boundary_config["compare_radius"],
        },
        "false_positive_guard": _blowup_guard(finest),
        "acceptance_checks": checks,
        "accepted_as_control": bool(all(checks.values())),
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    artifact_names = [
        "config.snapshot.json",
        "convergence.csv",
        "diagnostics.csv",
        "energy.svg",
        "profiles.npz",
        "profiles.svg",
        "summary.json",
    ]
    manifest = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "files": {
            name: {
                "sha256": _sha256(output_dir / name),
                "bytes": (output_dir / name).stat().st_size,
            }
            for name in artifact_names
        },
    }
    manifest_bytes = (
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    write_with_digest(
        output_dir / "manifest.json",
        manifest_bytes,
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "baseline.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "baseline",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="baseline config",
    )
    if not isinstance(config, dict):
        parser.error("baseline config must be a JSON object")
    summary = run(config, output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["accepted_as_control"] else 2


if __name__ == "__main__":
    sys.exit(main())
