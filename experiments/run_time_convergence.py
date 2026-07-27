#!/usr/bin/env python
"""Run a fixed-spatial-grid time-step convergence control.

This experiment reuses only the independent radial Crank--Nicolson benchmark
integrator from :mod:`experiments.run_baseline`.  The physical radial grid is
held fixed while the requested time steps are exactly

    dt, dt/2, dt/4.

Every run is compared with the analytic Gaussian swirl diffusion solution.
The experiment is a consistency check for a smooth, non-singular control; it
is not evidence for or against a Navier--Stokes singularity.
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

from experiments.run_baseline import Evolution, analytic_u1, evolve
from ns_certificate_lab._integrity import strict_json_loads, write_with_digest
from ns_certificate_lab.provenance import collect_runtime_provenance


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_TIME_STEP_DIVISORS = (1, 2, 4)


@dataclass(frozen=True)
class TimeRun:
    """One member of the fixed-grid time refinement study."""

    requested_dt: float
    actual_dt: float
    step_count: int
    evolution: Evolution
    relative_l2_error: float
    maximum_absolute_error: float
    auxiliary_boundary_all_time_inner_max_difference: float
    auxiliary_boundary_final_inner_max_difference: float


def _inside_repository(path: Path) -> bool:
    try:
        path.resolve().relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _trapezoid(values: np.ndarray, coordinates: np.ndarray) -> float:
    modern = getattr(np, "trapezoid", None)
    if modern is not None:
        return float(modern(values, coordinates))
    return float(np.trapz(values, coordinates))


def _relative_l2(
    r: np.ndarray,
    error: np.ndarray,
    reference: np.ndarray,
) -> float:
    numerator = _trapezoid(error * error * r**3, r)
    denominator = _trapezoid(reference * reference * r**3, r)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("analytic reference has invalid weighted L2 norm")
    return math.sqrt(numerator / denominator)


def _finite_positive(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite positive number")
    converted = float(value)
    if not math.isfinite(converted) or converted <= 0.0:
        raise ValueError(f"{name} must be a finite positive number")
    return converted


def _integer(value: Any, *, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _validate_config(config: dict[str, Any]) -> None:
    required = {
        "schema_version",
        "experiment_id",
        "coordinate_system",
        "normalization",
        "units",
        "seed",
        "viscosity",
        "amplitude",
        "sigma",
        "final_time",
        "r_max",
        "z_period",
        "nr",
        "base_time_step",
        "time_step_divisors",
        "sample_count",
        "boundary_sensitivity",
        "acceptance",
    }
    missing = sorted(required - config.keys())
    if missing:
        raise ValueError(f"time convergence config is missing: {missing}")
    if config["schema_version"] != 1:
        raise ValueError("unsupported time convergence config schema")
    if not isinstance(config["experiment_id"], str) or not config["experiment_id"]:
        raise ValueError("experiment_id must be a nonempty string")
    _integer(config["seed"], name="seed", minimum=0)
    _integer(config["nr"], name="nr", minimum=9)
    _integer(config["sample_count"], name="sample_count", minimum=2)
    for name in (
        "viscosity",
        "amplitude",
        "sigma",
        "final_time",
        "r_max",
        "z_period",
        "base_time_step",
    ):
        _finite_positive(config[name], name=name)
    divisors = config["time_step_divisors"]
    if (
        not isinstance(divisors, list)
        or any(isinstance(value, bool) or not isinstance(value, int) for value in divisors)
        or tuple(divisors) != EXPECTED_TIME_STEP_DIVISORS
    ):
        raise ValueError("time_step_divisors must be exactly [1, 2, 4]")

    final_time = float(config["final_time"])
    base_dt = float(config["base_time_step"])
    for divisor in EXPECTED_TIME_STEP_DIVISORS:
        ratio = final_time / (base_dt / divisor)
        if not math.isclose(ratio, round(ratio), rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("final_time must be an integer multiple of every time step")

    boundary = config["boundary_sensitivity"]
    if not isinstance(boundary, dict):
        raise ValueError("boundary_sensitivity must be an object")
    boundary_role = boundary.get("role")
    if not isinstance(boundary_role, str) or not boundary_role:
        raise ValueError("boundary_sensitivity role must be a nonempty string")
    radii = boundary.get("radii")
    if (
        not isinstance(radii, list)
        or len(radii) != 2
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or float(value) <= 0.0
            for value in radii
        )
        or float(radii[0]) >= float(radii[1])
    ):
        raise ValueError("boundary radii must be two increasing positive values")
    points_per_unit = _integer(
        boundary.get("points_per_unit"),
        name="boundary points_per_unit",
        minimum=4,
    )
    compare_radius = _finite_positive(
        boundary.get("compare_radius"),
        name="boundary compare_radius",
    )
    if compare_radius >= float(radii[0]):
        raise ValueError("boundary compare_radius must be inside the smaller cylinder")
    _finite_positive(
        boundary.get("maximum_allowed_difference"),
        name="maximum_allowed_difference",
    )
    if int(round(float(radii[0]) * points_per_unit)) + 1 < 9:
        raise ValueError("boundary grid must contain at least nine radial points")

    acceptance = config["acceptance"]
    if not isinstance(acceptance, dict):
        raise ValueError("acceptance must be an object")
    for name in (
        "minimum_observed_temporal_order",
        "maximum_finest_relative_l2_error",
        "energy_relative_increase_tolerance",
        "stability_peak_factor",
    ):
        _finite_positive(acceptance.get(name), name=f"acceptance {name}")


def temporal_orders(errors: list[float], time_steps: list[float]) -> list[float]:
    """Return observed orders using the actual time-step ratios."""

    if len(errors) != len(time_steps) or len(errors) < 2:
        raise ValueError("errors and time_steps must have equal length >= 2")
    orders: list[float] = []
    for coarse_error, fine_error, coarse_dt, fine_dt in zip(
        errors[:-1],
        errors[1:],
        time_steps[:-1],
        time_steps[1:],
    ):
        values = (coarse_error, fine_error, coarse_dt, fine_dt)
        if any(not math.isfinite(value) or value <= 0.0 for value in values):
            raise ValueError("errors and time steps must be finite and positive")
        if fine_dt >= coarse_dt:
            raise ValueError("time steps must be strictly decreasing")
        orders.append(
            float(
                math.log(coarse_error / fine_error)
                / math.log(coarse_dt / fine_dt)
            )
        )
    return orders


def step_doubling_diagnostics(
    runs: list[TimeRun],
    analytic_reference: np.ndarray,
) -> tuple[list[float], list[float]]:
    """Measure temporal order from consecutive numerical-solution differences.

    Because every member uses the identical spatial grid, the leading spatial
    discretization error is common to the three runs and largely cancels in
    ``u_dt - u_dt/2``.  The returned differences are normalized by the fixed
    analytic-reference norm; this normalization cancels from their ratio.
    """

    if len(runs) < 3:
        raise ValueError("step doubling requires at least three runs")
    if any(
        not np.array_equal(runs[0].evolution.r, item.evolution.r)
        for item in runs[1:]
    ):
        raise ValueError("step-doubling runs must share one spatial grid")
    differences = [
        _relative_l2(
            coarse.evolution.r,
            coarse.evolution.profiles[-1] - fine.evolution.profiles[-1],
            analytic_reference,
        )
        for coarse, fine in zip(runs[:-1], runs[1:])
    ]
    orders = temporal_orders(
        differences,
        [item.actual_dt for item in runs[:-1]],
    )
    return differences, orders


def _evolve_at_requested_dt(
    *,
    nr: int,
    r_max: float,
    final_time: float,
    viscosity: float,
    amplitude: float,
    sigma: float,
    requested_dt: float,
    sample_count: int,
    exact_outer_boundary: bool,
) -> tuple[Evolution, float, int]:
    dr = r_max / (nr - 1)
    expected_steps = int(round(final_time / requested_dt))
    actual_dt = final_time / expected_steps
    if not math.isclose(actual_dt, requested_dt, rel_tol=0.0, abs_tol=1e-14):
        raise ValueError("requested time step does not exactly partition final_time")
    evolution = evolve(
        nr=nr,
        r_max=r_max,
        final_time=final_time,
        viscosity=viscosity,
        amplitude=amplitude,
        sigma=sigma,
        dt_over_dr=requested_dt / dr,
        sample_count=max(sample_count, expected_steps + 1),
        exact_outer_boundary=exact_outer_boundary,
    )
    observed_steps = evolution.times.size - 1
    if observed_steps != expected_steps:
        raise RuntimeError(
            "baseline integrator did not preserve every requested time sample"
        )
    if not np.allclose(
        np.diff(evolution.times),
        requested_dt,
        rtol=0.0,
        atol=2e-15,
    ):
        raise RuntimeError("baseline integrator changed the requested time step")
    return evolution, actual_dt, expected_steps


def _auxiliary_boundary_radius_sensitivity(
    config: dict[str, Any],
    *,
    requested_dt: float,
) -> tuple[float, float]:
    """Return all-time and final inner differences on auxiliary R=3,4 grids.

    This intentionally does not use the main ``R=5`` convergence grid.  It is
    a small boundary-radius sensitivity test, not a bound on main-domain
    truncation error.
    """

    boundary = config["boundary_sensitivity"]
    evolutions: list[Evolution] = []
    for radius_value in boundary["radii"]:
        radius = float(radius_value)
        nr = int(round(radius * int(boundary["points_per_unit"]))) + 1
        evolution, _, _ = _evolve_at_requested_dt(
            nr=nr,
            r_max=radius,
            final_time=float(config["final_time"]),
            viscosity=float(config["viscosity"]),
            amplitude=float(config["amplitude"]),
            sigma=float(config["sigma"]),
            requested_dt=requested_dt,
            sample_count=int(config["sample_count"]),
            exact_outer_boundary=False,
        )
        evolutions.append(evolution)
    small, large = evolutions
    if not np.array_equal(small.times, large.times):
        raise RuntimeError("auxiliary boundary-radius runs have different times")
    compare_radius = float(boundary["compare_radius"])
    mask = small.r <= compare_radius + 1e-14
    differences: list[float] = []
    for index in range(small.times.size):
        large_on_small = np.interp(
            small.r[mask],
            large.r,
            large.profiles[index],
        )
        differences.append(
            float(
                np.max(
                    np.abs(small.profiles[index, mask] - large_on_small)
                )
            )
        )
    return max(differences), differences[-1]


def evaluate(config: dict[str, Any]) -> tuple[list[TimeRun], list[float]]:
    """Evaluate the temporal study without writing artifacts."""

    _validate_config(config)
    np.random.seed(int(config["seed"]))
    nr = int(config["nr"])
    r_max = float(config["r_max"])
    final_time = float(config["final_time"])
    viscosity = float(config["viscosity"])
    amplitude = float(config["amplitude"])
    sigma = float(config["sigma"])
    sample_count = int(config["sample_count"])
    time_steps = [
        float(config["base_time_step"]) / divisor
        for divisor in EXPECTED_TIME_STEP_DIVISORS
    ]

    runs: list[TimeRun] = []
    reference_r: np.ndarray | None = None
    for requested_dt in time_steps:
        evolution, actual_dt, step_count = _evolve_at_requested_dt(
            nr=nr,
            r_max=r_max,
            final_time=final_time,
            viscosity=viscosity,
            amplitude=amplitude,
            sigma=sigma,
            requested_dt=requested_dt,
            sample_count=sample_count,
            exact_outer_boundary=True,
        )
        if reference_r is None:
            reference_r = evolution.r.copy()
        elif not np.array_equal(reference_r, evolution.r):
            raise RuntimeError("spatial grid changed during time-step refinement")
        exact = analytic_u1(
            evolution.r,
            final_time,
            viscosity=viscosity,
            amplitude=amplitude,
            sigma=sigma,
        )
        error = evolution.profiles[-1] - exact
        boundary_all_time, boundary_final = (
            _auxiliary_boundary_radius_sensitivity(
                config,
                requested_dt=requested_dt,
            )
        )
        runs.append(
            TimeRun(
                requested_dt=requested_dt,
                actual_dt=actual_dt,
                step_count=step_count,
                evolution=evolution,
                relative_l2_error=_relative_l2(evolution.r, error, exact),
                maximum_absolute_error=float(np.max(np.abs(error))),
                auxiliary_boundary_all_time_inner_max_difference=(
                    boundary_all_time
                ),
                auxiliary_boundary_final_inner_max_difference=boundary_final,
            )
        )
    errors = [item.relative_l2_error for item in runs]
    actual_time_steps = [item.actual_dt for item in runs]
    return runs, temporal_orders(errors, actual_time_steps)


def _write_csv(path: Path, header: list[str], rows: list[list[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Run the study and write a checksummed evidence directory."""

    _validate_config(config)
    run_provenance = collect_runtime_provenance()
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite non-empty evidence directory: {output_dir}"
        )
    runs, orders = evaluate(config)
    errors = [item.relative_l2_error for item in runs]
    time_steps = [item.actual_dt for item in runs]
    analytic_reference = analytic_u1(
        runs[0].evolution.r,
        float(config["final_time"]),
        viscosity=float(config["viscosity"]),
        amplitude=float(config["amplitude"]),
        sigma=float(config["sigma"]),
    )
    step_doubling_differences, step_doubling_orders = (
        step_doubling_diagnostics(runs, analytic_reference)
    )

    energy_increase_fractions: list[float] = []
    peak_vorticity_factors: list[float] = []
    for item in runs:
        energies = item.evolution.energy_per_z_length
        initial_energy = max(float(energies[0]), 1e-300)
        energy_increase_fractions.append(
            float(max(0.0, np.max(np.diff(energies), initial=0.0) / initial_energy))
        )
        initial_vorticity = max(float(item.evolution.peak_vorticity[0]), 1e-300)
        peak_vorticity_factors.append(
            float(np.max(item.evolution.peak_vorticity) / initial_vorticity)
        )

    acceptance = config["acceptance"]
    boundary_limit = float(
        config["boundary_sensitivity"]["maximum_allowed_difference"]
    )
    finite_diagnostics = (
        errors
        + orders
        + step_doubling_differences
        + step_doubling_orders
        + energy_increase_fractions
        + peak_vorticity_factors
        + [item.maximum_absolute_error for item in runs]
        + [
            item.auxiliary_boundary_all_time_inner_max_difference
            for item in runs
        ]
        + [
            item.auxiliary_boundary_final_inner_max_difference
            for item in runs
        ]
        + [
            float(value)
            for item in runs
            for values in (
                item.evolution.energy_per_z_length,
                item.evolution.peak_vorticity,
            )
            for value in values
        ]
    )
    checks = {
        "time_steps_are_dt_dt_over_2_dt_over_4": bool(
            np.allclose(
                time_steps,
                [
                    float(config["base_time_step"]),
                    float(config["base_time_step"]) / 2.0,
                    float(config["base_time_step"]) / 4.0,
                ],
                rtol=0.0,
                atol=1e-14,
            )
        ),
        "spatial_grid_identical": bool(
            all(
                np.array_equal(runs[0].evolution.r, item.evolution.r)
                for item in runs[1:]
            )
        ),
        "all_values_finite": bool(
            all(item.evolution.all_finite for item in runs)
            and all(
                math.isfinite(value)
                for value in finite_diagnostics
            )
        ),
        "analytic_errors_decrease": bool(
            all(fine < coarse for coarse, fine in zip(errors, errors[1:]))
        ),
        "observed_temporal_order_passes": bool(
            orders
            and min(orders)
            >= float(acceptance["minimum_observed_temporal_order"])
        ),
        "step_doubling_order_passes": bool(
            step_doubling_orders
            and min(step_doubling_orders)
            >= float(acceptance["minimum_observed_temporal_order"])
        ),
        "finest_error_passes": bool(
            errors[-1]
            <= float(acceptance["maximum_finest_relative_l2_error"])
        ),
        "energy_nonincreasing_with_tolerance": bool(
            max(energy_increase_fractions)
            <= float(acceptance["energy_relative_increase_tolerance"])
        ),
        "stable_peak_vorticity": bool(
            max(peak_vorticity_factors)
            <= float(acceptance["stability_peak_factor"])
        ),
        "auxiliary_boundary_radius_all_time_test_passes": bool(
            max(
                item.auxiliary_boundary_all_time_inner_max_difference
                for item in runs
            )
            <= boundary_limit
        ),
        "auxiliary_boundary_radius_final_time_test_passes": bool(
            max(
                item.auxiliary_boundary_final_inner_max_difference
                for item in runs
            )
            <= boundary_limit
        ),
    }

    config_snapshot = output_dir / "config.snapshot.json"
    config_snapshot.write_text(
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _write_csv(
        output_dir / "time_convergence.csv",
        [
            "requested_dt",
            "actual_dt",
            "step_count",
            "nr",
            "dr",
            "relative_l2_error",
            "maximum_absolute_error",
            "analytic_error_order_from_previous",
            "step_doubling_difference_to_next",
            "step_doubling_order_from_previous_difference",
            "initial_energy_per_unit_z",
            "final_energy_per_unit_z",
            "maximum_vorticity_over_time",
            "final_vorticity",
            "auxiliary_boundary_all_time_inner_max_difference",
            "auxiliary_boundary_final_inner_max_difference",
        ],
        [
            [
                item.requested_dt,
                item.actual_dt,
                item.step_count,
                int(config["nr"]),
                float(config["r_max"]) / (int(config["nr"]) - 1),
                item.relative_l2_error,
                item.maximum_absolute_error,
                "" if index == 0 else orders[index - 1],
                (
                    ""
                    if index == len(runs) - 1
                    else step_doubling_differences[index]
                ),
                (
                    ""
                    if index == 0 or index == len(runs) - 1
                    else step_doubling_orders[index - 1]
                ),
                item.evolution.energy_per_z_length[0],
                item.evolution.energy_per_z_length[-1],
                np.max(item.evolution.peak_vorticity),
                item.evolution.peak_vorticity[-1],
                item.auxiliary_boundary_all_time_inner_max_difference,
                item.auxiliary_boundary_final_inner_max_difference,
            ]
            for index, item in enumerate(runs)
        ],
    )
    _write_csv(
        output_dir / "time_diagnostics.csv",
        ["dt", "time", "energy_per_unit_z", "peak_vorticity"],
        [
            [
                item.actual_dt,
                item.evolution.times[index],
                item.evolution.energy_per_z_length[index],
                item.evolution.peak_vorticity[index],
            ]
            for item in runs
            for index in range(item.evolution.times.size)
        ],
    )
    np.savez_compressed(
        output_dir / "final_profiles.npz",
        r=runs[0].evolution.r,
        dt=np.asarray(time_steps, dtype=np.float64),
        u1_final=np.stack([item.evolution.profiles[-1] for item in runs]),
        u1_exact=analytic_reference,
    )

    normalization = config["normalization"]
    length_scale = float(normalization["length_scale_L0"])
    velocity_scale = float(normalization["velocity_scale_U0"])
    summary: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Fixed-grid temporal convergence of a smooth non-singular control. "
            "This is neither a singularity candidate nor a mathematical proof."
        ),
        "method": {
            "equation": "u1_t = nu * (u1_rr + 3*u1_r/r)",
            "time_integrator": "Crank-Nicolson, theta=1/2",
            "spatial_operator": "independent second-order centered radial stencil",
            "fixed_spatial_grid": {
                "nr": int(config["nr"]),
                "r_max": float(config["r_max"]),
                "dr": float(config["r_max"]) / (int(config["nr"]) - 1),
            },
            "time_step_sequence": "dt, dt/2, dt/4",
            "outer_boundary_for_analytic_error": "analytic Dirichlet trace",
            "auxiliary_boundary_radius_test": (
                "homogeneous outer Dirichlet traces on auxiliary R=3 and R=4 "
                "cylinders; this is not a direct boundary-error proof for the "
                "main R=5 convergence domain"
            ),
        },
        "physical_interpretation": {
            "coordinate_system": config["coordinate_system"],
            "units": config["units"],
            "normalization": {
                "length_scale_L0": length_scale,
                "velocity_scale_U0": velocity_scale,
                "time_scale_L0_over_U0": length_scale / velocity_scale,
                "reynolds_number_U0_L0_over_nu": (
                    velocity_scale * length_scale / float(config["viscosity"])
                ),
            },
        },
        "reproducibility": {
            "seed": int(config["seed"]),
            "config_sha256": _sha256(config_snapshot),
            "runtime_provenance": run_provenance,
        },
        "temporal_convergence": {
            "requested_time_steps": [item.requested_dt for item in runs],
            "actual_time_steps": time_steps,
            "step_counts": [item.step_count for item in runs],
            "relative_l2_errors": errors,
            "maximum_absolute_errors": [
                item.maximum_absolute_error for item in runs
            ],
            "analytic_error_observed_orders": orders,
            "analytic_error_order_definition": (
                "log(E_dt/E_dt_over_2) / log(dt/(dt/2)), where E is the "
                "raw weighted-L2 error against the analytic solution"
            ),
            "step_doubling_relative_differences": step_doubling_differences,
            "step_doubling_observed_orders": step_doubling_orders,
            "step_doubling_order_definition": (
                "log(||u_dt-u_dt_over_2||/||u_dt_over_2-u_dt_over_4||) "
                "/ log(2), with all profiles on the identical radial grid"
            ),
            "spatial_error_note": (
                "Raw analytic errors include the fixed spatial error and can "
                "eventually plateau. Step-doubling differences are also "
                "reported because their common spatial error largely cancels."
            ),
        },
        "energy_and_vorticity": [
            {
                "dt": item.actual_dt,
                "initial_energy_per_unit_z": float(
                    item.evolution.energy_per_z_length[0]
                ),
                "final_energy_per_unit_z": float(
                    item.evolution.energy_per_z_length[-1]
                ),
                "maximum_relative_energy_increase": energy_increase_fractions[index],
                "maximum_vorticity_over_time": float(
                    np.max(item.evolution.peak_vorticity)
                ),
                "final_vorticity": float(item.evolution.peak_vorticity[-1]),
                "peak_vorticity_factor": peak_vorticity_factors[index],
            }
            for index, item in enumerate(runs)
        ],
        "auxiliary_boundary_radius_test": {
            "scope": config["boundary_sensitivity"]["role"],
            "main_domain_r_max": float(config["r_max"]),
            "limitation": (
                "The auxiliary R=3 versus R=4 comparison does not directly "
                "measure or prove the truncation error on the main R=5 domain."
            ),
            "radii": config["boundary_sensitivity"]["radii"],
            "points_per_unit": int(
                config["boundary_sensitivity"]["points_per_unit"]
            ),
            "compare_radius": float(
                config["boundary_sensitivity"]["compare_radius"]
            ),
            "inner_differences_by_dt": [
                {
                    "dt": item.actual_dt,
                    "all_time_max_difference": (
                        item.auxiliary_boundary_all_time_inner_max_difference
                    ),
                    "final_time_difference": (
                        item.auxiliary_boundary_final_inner_max_difference
                    ),
                }
                for item in runs
            ],
        },
        "acceptance_checks": checks,
        "accepted_as_temporal_control": bool(all(checks.values())),
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    artifact_names = [
        "config.snapshot.json",
        "final_profiles.npz",
        "summary.json",
        "time_convergence.csv",
        "time_diagnostics.csv",
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
    write_with_digest(
        output_dir / "manifest.json",
        (
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8"),
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "baseline_time_convergence.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "baseline_time_convergence",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="time convergence config",
    )
    if not isinstance(config, dict):
        parser.error("time convergence config must be a JSON object")
    try:
        summary = run(config, output_dir)
    except (ValueError, FileExistsError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["accepted_as_temporal_control"] else 2


if __name__ == "__main__":
    sys.exit(main())
