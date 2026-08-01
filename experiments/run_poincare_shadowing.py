r"""Continuum-to-lattice shadowing scaffolding for the doubling Poincare map.

Verification-sprint workstream E.  Pilot A (``run_renormalized_cascade.py``)
found **no** periodic orbit in its scanned box, so no orbit-based Floquet
radius exists and no orbit-shadowing theorem can be attempted.  What can still
be measured is the *scaffolding* a shadowing argument would need:

1. the restriction ``R_N`` and interpolation ``I_N`` operators that connect a
   continuum cell profile ``U_a(xi)`` to lattice Fourier coefficients;
2. the **consistency error** of the lattice Leray convolution against the
   continuum quadratic form ``Q``, measured by comparing scale ``N`` and scale
   ``2N`` restrictions of the *same* continuum profile through the doubling
   pairing ``2k <-> k``;
3. the same consistency at the level of one full evolve+pullback+renormalize
   stage (the one-stage Poincare map ``P``);
4. the **local Lipschitz factor** ``L`` of ``P`` along the measured
   trajectory, from which the contraction margin ``1 - L`` follows.

Pre-registered acceptance rule (fixed before the numbers were produced): the
orbit-shadowing proof candidate survives only if

    contraction margin (1 - L)  >  map-level consistency error.

Everything here is binary64 (``numpy`` only, no ``scipy``); the fitted exponent
``sigma`` is reported as a **diagnostic exponent, not a proven rate**, per the
repository's TM-22 discipline.  No continuum enclosure, no PDE statement.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    write_with_digest,
)
from ns_certificate_lab.leray_response_relay import mean_energy
from ns_certificate_lab.mesoscopic_cloud_scaling import MesoscopicCloudConfig
import ns_certificate_lab.mesoscopic_local_fft as local_fft
import ns_certificate_lab.renormalized_cascade as cascade
from ns_certificate_lab.mesoscopic_galerkin import build_angle_box_parent

OUTPUT_SCHEMA = "ns-certificate-lab/poincare-shadowing-scaffolding/v1"
STATUS = "BINARY64 VERIFICATION DIAGNOSTIC / NOT A PROOF"

# eta = 3/16 fixed-relative two-box family.  ``floor(3 N / 16)`` realises the
# relative half-width exactly when 16 | N; the other N are kept as recorded
# geometry-mismatched controls, never deleted.
ETA_NUMERATOR = 3
ETA_DENOMINATOR = 16
CHILD_DIRECTION = (2, 1, 1)  # p + q with p=(1,1,0), q=(1,0,1)

# The continuum phase slope acts on the cell variable xi = q / N.  It is chosen
# so that the N = 16 member reproduces the repository's default integer-offset
# slope exactly.
CONTINUUM_PHASE_SLOPE = (
    0.173 * ETA_DENOMINATOR,
    -0.119 * ETA_DENOMINATOR,
    0.071 * ETA_DENOMINATOR,
)
LATTICE_PHASE_SLOPE = (0.173, -0.119, 0.071)

PERTURBATION_RELATIVE_SIZE = 1.0e-6
PERTURBATION_COUNT = 3
RANDOM_SEED = 20260802


# --------------------------------------------------------------------------
# Fits
# --------------------------------------------------------------------------


def log_log_fit(
    abscissae: Sequence[float], values: Sequence[float]
) -> dict[str, float | int | None]:
    """Return ``sigma`` and ``R^2`` for ``value ~ C * x ** (-sigma)``."""

    x = np.log(np.asarray(abscissae, dtype=np.float64))
    y = np.log(np.asarray(values, dtype=np.float64))
    if x.size < 3:
        return {
            "points": int(x.size),
            "sigma": None,
            "log_prefactor": None,
            "r_squared": None,
        }
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    residual = float(np.sum((y - predicted) ** 2))
    total = float(np.sum((y - float(np.mean(y))) ** 2))
    r_squared = 1.0 - residual / total if total > 0.0 else None
    return {
        "points": int(x.size),
        "sigma": float(-slope),
        "log_prefactor": float(intercept),
        "r_squared": r_squared,
    }


# --------------------------------------------------------------------------
# Part A -- convolution consistency
# --------------------------------------------------------------------------


def cloud_width(base_scale: int) -> int:
    return max(1, (ETA_NUMERATOR * base_scale) // ETA_DENOMINATOR)


def relative_width(base_scale: int) -> float:
    return cloud_width(base_scale) / float(base_scale)


def child_band_coefficients(
    base_scale: int, half_width: int, *, continuum_profile: bool
) -> np.ndarray:
    """Return the Leray convolution on the child band of one cloud.

    With ``u_hat(k) = N^-2 Psi(k/N)`` the lattice convolution coefficient at
    ``k = N xi`` approximates ``-Q(Psi,Psi)(xi)`` with an ``N^0`` prefactor, so
    the returned blocks at two scales are directly comparable.
    """

    width = cloud_width(base_scale)
    slope = (
        tuple(value / base_scale for value in CONTINUUM_PHASE_SLOPE)
        if continuum_profile
        else LATTICE_PHASE_SLOPE
    )
    config = MesoscopicCloudConfig(
        base_scale=base_scale,
        gamma=1.0,
        width_override=width,
        phase_slope=slope,
    )
    _parents, _parent_blocks, blocks, _bytes = local_fft._run_with_blocks(
        config, maximum_working_bytes=1_500_000_000
    )
    center = tuple(base_scale * value for value in CHILD_DIRECTION)
    lower = tuple(value - half_width for value in center)
    upper = tuple(value + half_width for value in center)
    field = local_fft._field_on_region(blocks, lower, upper)
    del blocks, _parents, _parent_blocks
    return np.asarray(field, dtype=np.complex128)


def convolution_consistency(
    base_scale: int, *, continuum_profile: bool
) -> dict[str, Any]:
    """Compare ``Q_N`` and ``Q_{2N}`` through the doubling pairing."""

    coarse_half = cloud_width(base_scale) - 1
    fine_half = cloud_width(2 * base_scale) - 1
    common = min(coarse_half, fine_half // 2)
    coarse = child_band_coefficients(
        base_scale, common, continuum_profile=continuum_profile
    )
    fine = child_band_coefficients(
        2 * base_scale, 2 * common, continuum_profile=continuum_profile
    )
    paired = fine[:, ::2, ::2, ::2]
    if paired.shape != coarse.shape:
        raise AssertionError("doubling pairing produced mismatched blocks")
    numerator = float(np.linalg.norm(coarse - paired))
    denominator = float(np.linalg.norm(paired))
    if not denominator > 0.0:
        raise ValueError("the fine child band vanished")
    geometry_matched = math.isclose(
        relative_width(base_scale),
        relative_width(2 * base_scale),
        rel_tol=0.0,
        abs_tol=0.0,
    )
    return {
        "base_scale": base_scale,
        "partner_scale": 2 * base_scale,
        "width_coarse": cloud_width(base_scale),
        "width_fine": cloud_width(2 * base_scale),
        "relative_width_coarse": relative_width(base_scale),
        "relative_width_fine": relative_width(2 * base_scale),
        "geometry_matched": bool(geometry_matched),
        "common_half_width": common,
        "common_mode_count": int((2 * common + 1) ** 3),
        "fine_band_norm": denominator,
        "relative_error": numerator / denominator,
    }


# --------------------------------------------------------------------------
# Parts B and C -- one-stage Poincare map
# --------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class StageOperator:
    """One evolve + doubling pullback + optional sea drop + renormalize."""

    config: cascade.CascadeConfig
    wave_squared: np.ndarray
    galerkin_mask: np.ndarray
    linf: np.ndarray
    steps: int

    @property
    def pullback_cutoff(self) -> int:
        return self.config.cutoff // 2

    def apply(self, state: np.ndarray) -> np.ndarray:
        fixed = dataclasses.replace(
            self.config, base_steps=self.steps, max_steps=self.steps
        )
        evolved, used = cascade._evolve(
            state,
            config=fixed,
            wave_squared=self.wave_squared,
            galerkin_mask=self.galerkin_mask,
        )
        if used != self.steps:
            raise AssertionError("the stage integrator changed its step count")
        pulled = cascade.doubling_pullback(evolved)
        if self.config.drop_below is not None:
            keep = ~(self.linf < self.config.drop_below)
            pulled = np.asarray(
                pulled * keep[None, ...], dtype=np.complex128
            )
        energy = mean_energy(pulled)
        if not energy > 0.0:
            raise ValueError("the pullback emptied the state")
        target = self.config.energy_constant / float(self.config.scale)
        return np.asarray(
            pulled * math.sqrt(target / energy), dtype=np.complex128
        )


def build_stage(
    *,
    scale: int,
    width: int,
    grid_size: int,
    energy_constant: float,
    drop_below: int | None,
) -> tuple[StageOperator, np.ndarray]:
    config = cascade.CascadeConfig(
        scale=scale,
        width=width,
        grid_size=grid_size,
        energy_constant=energy_constant,
        stages=1,
        drop_below=drop_below,
    )
    config.validate()
    if 2 * config.cutoff >= grid_size // 2:
        raise ValueError("grid is too small for exact padded Galerkin products")
    parent, _child = build_angle_box_parent(
        grid_size,
        scale=scale,
        width=width,
        energy_constant=energy_constant,
    )
    kx, ky, kz = cascade._frequency_mesh(grid_size)
    wave_squared = kx * kx + ky * ky + kz * kz
    linf = np.maximum(np.abs(kx), np.maximum(np.abs(ky), np.abs(kz)))
    galerkin_mask = cascade._galerkin_mask(grid_size, config.cutoff)
    advective_rate = config.cutoff * cascade._sup_velocity(parent)
    steps = config.base_steps
    if advective_rate > 0.0:
        needed = int(
            math.ceil(config.stage_time * advective_rate * config.cfl_safety)
        )
        steps = min(max(steps, needed), config.max_steps)
    operator = StageOperator(
        config=config,
        wave_squared=wave_squared,
        galerkin_mask=galerkin_mask,
        linf=linf,
        steps=steps,
    )
    return operator, parent


def profile_block(state: np.ndarray, scale: int, half: int) -> np.ndarray:
    """Return ``Psi(k/scale) = scale^2 u_hat(k)`` on ``|k|_inf <= half``."""

    grid = state.shape[1]
    index = np.arange(-half, half + 1, dtype=np.int64) % grid
    block = state[np.ix_(range(3), index, index, index)]
    return np.asarray(float(scale) ** 2 * block, dtype=np.complex128)


def map_consistency(
    coarse: StageOperator,
    coarse_state: np.ndarray,
    fine: StageOperator,
    fine_state: np.ndarray,
) -> dict[str, Any]:
    """Compare one stage at ``N_0`` and at ``2 N_0`` through ``2k <-> k``."""

    if fine.config.scale != 2 * coarse.config.scale:
        raise ValueError("the fine stage must sit at twice the coarse scale")
    coarse_image = coarse.apply(coarse_state)
    fine_image = fine.apply(fine_state)
    half = min(coarse.pullback_cutoff, fine.pullback_cutoff // 2)
    coarse_profile = profile_block(coarse_image, coarse.config.scale, half)
    fine_profile = profile_block(fine_image, fine.config.scale, 2 * half)
    paired = fine_profile[:, ::2, ::2, ::2]
    if paired.shape != coarse_profile.shape:
        raise AssertionError("doubling pairing produced mismatched blocks")
    numerator = float(np.linalg.norm(coarse_profile - paired))
    denominator = float(np.linalg.norm(paired))
    overlap = abs(complex(np.vdot(coarse_profile, paired))) / max(
        float(np.linalg.norm(coarse_profile)) * denominator,
        np.finfo(float).tiny,
    )
    return {
        "coarse_scale": coarse.config.scale,
        "coarse_width": coarse.config.width,
        "coarse_grid": coarse.config.grid_size,
        "coarse_steps": coarse.steps,
        "fine_scale": fine.config.scale,
        "fine_width": fine.config.width,
        "fine_grid": fine.config.grid_size,
        "fine_steps": fine.steps,
        "drop_below": coarse.config.drop_below,
        "profiles_are_relatively_matched": bool(
            (coarse.config.width - 1) * fine.config.scale
            == (fine.config.width - 1) * coarse.config.scale
        ),
        "common_half_width": half,
        "common_mode_count": int((2 * half + 1) ** 3),
        "coarse_profile_norm": float(np.linalg.norm(coarse_profile)),
        "fine_profile_norm": denominator,
        "shape_overlap": overlap,
        "relative_error": numerator / denominator,
    }


def divergence_free_perturbation(
    generator: np.random.Generator,
    reference: np.ndarray,
    galerkin_mask: np.ndarray,
    grid_size: int,
) -> np.ndarray:
    """Return a real, mean-zero, divergence-free Fourier perturbation."""

    shape = (3, grid_size, grid_size, grid_size)
    raw = generator.standard_normal(shape) + 1.0j * generator.standard_normal(
        shape
    )
    raw = np.asarray(raw * galerkin_mask[None, ...], dtype=np.complex128)
    index = (-np.arange(grid_size, dtype=np.int64)) % grid_size
    reflected = np.take(raw, index, axis=1)
    reflected = np.take(reflected, index, axis=2)
    reflected = np.take(reflected, index, axis=3)
    raw = 0.5 * (raw + np.conjugate(reflected))
    kx, ky, kz = cascade._frequency_mesh(grid_size)
    wave = np.stack((kx, ky, kz), axis=0)
    wave_squared = np.sum(wave * wave, axis=0)
    projection = np.sum(wave * raw, axis=0)
    safe = np.where(wave_squared > 0.0, wave_squared, 1.0)
    raw = raw - wave * (projection / safe)[None, ...]
    raw[:, 0, 0, 0] = 0.0
    norm = float(np.linalg.norm(raw))
    if not norm > 0.0:
        raise ValueError("the random perturbation vanished")
    target = PERTURBATION_RELATIVE_SIZE * float(np.linalg.norm(reference))
    return np.asarray(raw * (target / norm), dtype=np.complex128)


def contraction_margin(
    operator: StageOperator, state: np.ndarray, generator: np.random.Generator
) -> dict[str, Any]:
    """Return the local Lipschitz factors of ``P`` at ``state``."""

    image = operator.apply(state)
    factors: list[dict[str, float]] = []
    for trial in range(PERTURBATION_COUNT):
        delta = divergence_free_perturbation(
            generator,
            state,
            operator.galerkin_mask,
            operator.config.grid_size,
        )
        perturbed_image = operator.apply(
            np.asarray(state + delta, dtype=np.complex128)
        )
        delta_norm = float(np.linalg.norm(delta))
        response = float(np.linalg.norm(perturbed_image - image))
        factors.append(
            {
                "trial": trial,
                "delta_norm": delta_norm,
                "response_norm": response,
                "lipschitz": response / delta_norm,
                "divergence_defect": cascade._divergence_defect(delta),
                "reality_defect": cascade._reality_defect(delta),
            }
        )
    values = [entry["lipschitz"] for entry in factors]
    return {
        "scale": operator.config.scale,
        "width": operator.config.width,
        "grid_size": operator.config.grid_size,
        "energy_constant": operator.config.energy_constant,
        "drop_below": operator.config.drop_below,
        "steps": operator.steps,
        "state_norm": float(np.linalg.norm(state)),
        "trials": factors,
        "lipschitz_min": min(values),
        "lipschitz_mean": float(np.mean(values)),
        "lipschitz_max": max(values),
        "contraction_margin_from_max": 1.0 - max(values),
    }


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------


def run(output_dir: Path) -> dict[str, Any]:
    scales = (8, 16, 24, 32, 48, 64)

    convolution_rows: list[dict[str, Any]] = []
    for base_scale in scales:
        row = convolution_consistency(base_scale, continuum_profile=True)
        row["profile"] = "continuum-fixed"
        convolution_rows.append(row)
        print(
            f"[A] N={base_scale:3d}->{2 * base_scale:3d} "
            f"matched={row['geometry_matched']} "
            f"err={row['relative_error']:.6e}",
            flush=True,
        )
    control_rows: list[dict[str, Any]] = []
    for base_scale in scales:
        row = convolution_consistency(base_scale, continuum_profile=False)
        row["profile"] = "lattice-fixed-phase-slope-null-control"
        control_rows.append(row)
        print(
            f"[A-null] N={base_scale:3d} err={row['relative_error']:.6e}",
            flush=True,
        )

    clean = [row for row in convolution_rows if row["geometry_matched"]]
    convolution_fit = log_log_fit(
        [row["base_scale"] for row in clean],
        [row["relative_error"] for row in clean],
    )
    control_fit = log_log_fit(
        [row["base_scale"] for row in control_rows if row["geometry_matched"]],
        [
            row["relative_error"]
            for row in control_rows
            if row["geometry_matched"]
        ],
    )

    map_rows: list[dict[str, Any]] = []
    for energy_constant in (1.0, 100.0):
        coarse, coarse_state = build_stage(
            scale=4,
            width=2,
            grid_size=64,
            energy_constant=energy_constant,
            drop_below=None,
        )
        for width, grid_size, label in ((2, 80, "prescribed"), (3, 96, "profile-matched")):
            fine, fine_state = build_stage(
                scale=8,
                width=width,
                grid_size=grid_size,
                energy_constant=energy_constant,
                drop_below=None,
            )
            row = map_consistency(coarse, coarse_state, fine, fine_state)
            row["energy_constant"] = energy_constant
            row["pairing"] = label
            map_rows.append(row)
            print(
                f"[B] c_E={energy_constant} {label} "
                f"err={row['relative_error']:.6e} "
                f"overlap={row['shape_overlap']:.6f}",
                flush=True,
            )

    generator = np.random.default_rng(RANDOM_SEED)
    contraction_rows: list[dict[str, Any]] = []
    for energy_constant in (1.0, 100.0):
        operator, state = build_stage(
            scale=4,
            width=2,
            grid_size=64,
            energy_constant=energy_constant,
            drop_below=3,
        )
        row = contraction_margin(operator, state, generator)
        contraction_rows.append(row)
        print(
            f"[C] c_E={energy_constant} L in "
            f"[{row['lipschitz_min']:.6e}, {row['lipschitz_max']:.6e}] "
            f"margin={row['contraction_margin_from_max']:.6e}",
            flush=True,
        )

    matched_map_errors = [
        row["relative_error"]
        for row in map_rows
        if row["pairing"] == "profile-matched"
    ]
    worst_map_error = max(matched_map_errors)
    worst_lipschitz = max(row["lipschitz_max"] for row in contraction_rows)
    margin = 1.0 - worst_lipschitz
    verdict = {
        "pre_registered_rule": (
            "orbit-shadowing proof candidate only if the contraction margin "
            "(1 - L) strictly exceeds the map-level consistency error"
        ),
        "worst_map_consistency_error": worst_map_error,
        "worst_local_lipschitz": worst_lipschitz,
        "contraction_margin": margin,
        "rule_satisfied": bool(margin > worst_map_error),
        "orbit_exists_in_scanned_box": False,
        "conclusion": (
            "KILLED: pilot A found no periodic orbit, so no Floquet radius "
            "exists, and the pre-registered margin rule is not met either"
            if not margin > worst_map_error
            else "CONDITIONAL: the margin rule is met but no orbit exists to "
            "shadow, so only the scaffolding survives"
        ),
    }

    summary = {
        "schema": OUTPUT_SCHEMA,
        "status": STATUS,
        "arithmetic": "binary64 (numpy); no interval enclosure, no exact lane",
        "operators": {
            "restriction": (
                "(R_N U)_{N a + q} = N^-2 U_a(q/N) for q in the cell lattice"
            ),
            "interpolation": (
                "(I_N u)_a(xi) = trilinear interpolation of N^2 u_{N a + q} "
                "at the nodes xi = q/N"
            ),
            "left_inverse": "R_N I_N = id on lattice data (nodal exactness)",
            "interpolation_defect_bound": (
                "|| (I_N R_N - id) U ||_inf <= sqrt(3) || grad U ||_inf / N "
                "on C^1 cell profiles"
            ),
        },
        "convolution_consistency": {
            "rows": convolution_rows,
            "null_control_rows": control_rows,
            "fit_on_geometry_matched_rows": convolution_fit,
            "null_control_fit": control_fit,
            "fit_discipline": (
                "diagnostic exponent per TM-22; not a proven convergence rate"
            ),
        },
        "map_consistency": {"rows": map_rows},
        "contraction": {"rows": contraction_rows},
        "verdict": verdict,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    write_with_digest(
        output_dir / "summary.json", canonical_json_bytes(summary)
    )
    return summary


CSV_HEADER = (
    "kind",
    "n_or_pair",
    "c_E",
    "value",
    "sigma_fit",
    "r_squared",
    "notes",
)


def write_csv(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fit = summary["convolution_consistency"]["fit_on_geometry_matched_rows"]
    control_fit = summary["convolution_consistency"]["null_control_fit"]
    rows: list[tuple[Any, ...]] = []
    for row in summary["convolution_consistency"]["rows"]:
        matched = row["geometry_matched"]
        rows.append(
            (
                "convolution_consistency",
                f"{row['base_scale']}->{row['partner_scale']}",
                "",
                f"{row['relative_error']:.9e}",
                f"{fit['sigma']:.6f}" if matched and fit["sigma"] else "",
                f"{fit['r_squared']:.6f}"
                if matched and fit["r_squared"] is not None
                else "",
                (
                    "float; continuum-fixed profile; eta=3/16 exact on both "
                    f"members; {row['common_mode_count']} paired modes"
                    if matched
                    else "float; RECORDED BUT EXCLUDED FROM FIT: realized "
                    f"relative width {row['relative_width_coarse']:.4f} vs "
                    f"{row['relative_width_fine']:.4f} (floor(3N/16) is exact "
                    "only for 16|N)"
                ),
            )
        )
    for row in summary["convolution_consistency"]["null_control_rows"]:
        if not row["geometry_matched"]:
            continue
        rows.append(
            (
                "convolution_consistency",
                f"{row['base_scale']}->{row['partner_scale']}",
                "",
                f"{row['relative_error']:.9e}",
                f"{control_fit['sigma']:.6f}"
                if control_fit["sigma"] is not None
                else "",
                f"{control_fit['r_squared']:.6f}"
                if control_fit["r_squared"] is not None
                else "",
                "float; NULL CONTROL: N-independent integer phase slope, so "
                "the profile is not fixed in xi and the error must not decay",
            )
        )
    for row in summary["map_consistency"]["rows"]:
        rows.append(
            (
                "map_consistency",
                f"{row['coarse_scale']}->{row['fine_scale']}",
                f"{row['energy_constant']:g}",
                f"{row['relative_error']:.9e}",
                "",
                "",
                (
                    f"float; {row['pairing']} pair "
                    f"(W={row['coarse_width']}/{row['fine_width']}, "
                    f"grid={row['coarse_grid']}/{row['fine_grid']}, "
                    f"steps={row['coarse_steps']}/{row['fine_steps']}); "
                    f"drop_below={row['drop_below']}; "
                    f"shape overlap {row['shape_overlap']:.6f}"
                ),
            )
        )
    for row in summary["contraction"]["rows"]:
        for trial in row["trials"]:
            rows.append(
                (
                    "contraction",
                    f"{row['scale']}",
                    f"{row['energy_constant']:g}",
                    f"{trial['lipschitz']:.9e}",
                    "",
                    "",
                    (
                        f"float; trial {trial['trial']}; "
                        f"|delta|/|u|={PERTURBATION_RELATIVE_SIZE:g}; "
                        f"front-only drop_below={row['drop_below']}; "
                        f"steps={row['steps']}; margin "
                        f"1-L={1.0 - trial['lipschitz']:.6e}"
                    ),
                )
            )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(CSV_HEADER)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/verification_sprint_v1/shadowing"),
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(
            "docs/research_notes/verification_sprint_v1/"
            "poincare_shadowing_scaling.csv"
        ),
    )
    arguments = parser.parse_args()
    summary = run(arguments.output_dir)
    write_csv(arguments.csv, summary)
    print(json.dumps(summary["verdict"], indent=2))


if __name__ == "__main__":
    main()
