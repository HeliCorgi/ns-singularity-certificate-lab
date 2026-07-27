#!/usr/bin/env python
"""Run the early-time Hou finite-cylinder observation on uniform fixed grids.

The initial datum is E-29,

    u1(0,r,z) = 12000 (1-r^2)^18 sin(2 pi z) / (1 + 12.5 sin^2(pi z)),
    omega1(0,r,z) = 0,

integrated with the audited two-stage viscosity protocol of E-30 by the
production integrator :mod:`ns_certificate_lab.nonlinear_cylinder`.  Every
snapshot is checkpointed and re-examined with the independent finite-cylinder
elliptic solver B.

This experiment is a *numerical observation*.  A uniform fixed grid cannot
resolve the adaptive mesh scales of the source calculation, so nothing here
reproduces, confirms or refutes any published amplification figure, and none
of it is evidence for singularity formation.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field as dataclass_field
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

from ns_certificate_lab._integrity import strict_json_loads, write_with_digest
from ns_certificate_lab.finite_cylinder_poisson import solve_finite_cylinder_poisson
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    DIAGNOSTIC_FIELDS,
    RELATIVE_DIAGNOSTIC_FIELDS,
    ConstrainedState,
    cartesian_vorticity,
    constrain_state,
    hou_initial_swirl,
    integrate,
    normalize_viscosity_schedule,
    save_checkpoint,
    viscosity_at,
)
from ns_certificate_lab.provenance import collect_runtime_provenance


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

# E-29: the audited amplitude is part of the equation, not a tunable input.
E29_AMPLITUDE = 12000.0
# E-29b: derived (not published) initial norms at amplitude_scale = 1.
E29B_MAX_ABS_U1 = 3265.9863
E29B_MAX_CARTESIAN_VORTICITY = 7569.62

# --------------------------------------------------------------- CFL policy
#
# ``adaptive_time_step`` sizes a step from the velocity maxima of the state at
# the START of that step; ``compute_diagnostics`` reports the effective CFL
# from the state the step PRODUCED.  With the advective bound active and
# frozen velocities the two coincide exactly, because
# ``max(|u^r|/dr, |u^z|/dz)`` is the reciprocal of ``min(dr/|u^r|, dz/|u^z|)``.
# The recorded effective CFL therefore exceeds ``cfl_coefficient`` by exactly
# the intra-step relative growth of the dominant directional velocity ratio.
# Clipping to snapshot / viscosity-switch / ``t_final`` targets and the
# ``max_time_step`` cap only ever shrink ``dt``, so they can only lower the
# effective CFL.  The delegated HH21 guidance recorded in
# ``docs/hou_setup_audit.md`` bounds the one-step relative growth of the
# maximum by one percent; the shipped ``outputs/hou_early_time_v1`` run
# measured 0.23 percent at its worst step.  The acceptance check therefore
# compares against ``C (1 + cfl_excess_tolerance)`` instead of ``C``.
DEFAULT_CFL_EXCESS_TOLERANCE = 0.05
V1_CFL_COEFFICIENT = 0.1
V1_MAXIMUM_ADVECTIVE_CFL = 0.10022676304450114
V1_ADVECTIVE_CFL_EXCESS_RATIO = (
    V1_MAXIMUM_ADVECTIVE_CFL / V1_CFL_COEFFICIENT - 1.0
)

CFL_RULE = (
    "dt <= C * min(dr / max|u^r|, dz / max|u^z|) evaluated at the START of the "
    "step, intersected with the viscous bound C * min(dr,dz)^2 / (4 nu), the "
    "max_time_step cap and the clipping to the next snapshot, viscosity "
    "switch or t_final target"
)
EFFECTIVE_CFL_DEFINITION = (
    "dt * max(max|u^r| / dr, max|u^z| / dz) evaluated on the state AFTER the "
    "step (diagnostics key advective_cfl)"
)
CFL_EXCESS_MECHANISM = (
    "The step is sized before the velocities grow. With the advective bound "
    "active and frozen velocities the post-step effective CFL equals C "
    "exactly, because max(|u^r|/dr, |u^z|/dz) = 1 / min(dr/|u^r|, dz/|u^z|); "
    "the measured excess is therefore exactly the one-step relative growth of "
    "the dominant directional velocity ratio. Step clipping and the "
    "max_time_step cap only shrink dt and so only lower the effective CFL."
)


def advective_cfl_within_tolerance(
    max_effective_cfl: float,
    *,
    cfl_coefficient: float,
    cfl_excess_tolerance: float = DEFAULT_CFL_EXCESS_TOLERANCE,
) -> bool:
    """Return whether a post-step effective CFL is inside the allowed excess.

    The rule is ``max_effective_cfl <= C * (1 + cfl_excess_tolerance)``.  A
    non-finite measurement never passes.  ``cfl_excess_tolerance = 0`` demands
    the post-step CFL never exceed ``C`` at all, which the start-of-step
    selection rule cannot guarantee (the v1 run recorded
    ``0.10022676304450114`` against ``C = 0.1``).
    """

    coefficient = _number(cfl_coefficient, name="cfl_coefficient", positive=True)
    tolerance = _number(cfl_excess_tolerance, name="cfl_excess_tolerance")
    if tolerance < 0.0:
        raise ValueError("cfl_excess_tolerance must be nonnegative")
    if isinstance(max_effective_cfl, bool) or not isinstance(
        max_effective_cfl, (int, float)
    ):
        return False
    measured = float(max_effective_cfl)
    if not math.isfinite(measured):
        return False
    return measured <= coefficient * (1.0 + tolerance)


LIMITATIONS = (
    "A uniform fixed grid cannot resolve the adaptive mesh scales (minimum "
    "spacing of order 1e-8) used by the source calculation; every quantity "
    "reported here is resolution limited.",
    "This is a numerical observation, not a reproduction claim and not "
    "evidence of singularity formation.",
    "The full z period is integrated and the odd symmetry of E-28 is "
    "monitored as a diagnostic rather than imposed.",
    "Solver A and solver B share the same second-order accuracy class, so "
    "their agreement bounds discretization consistency, not continuum error.",
    "The wall vorticity uses the repository's E-31 discretization, which the "
    "source calculation does not print; no filter of any kind is applied.",
    "Floating-point arithmetic only: there is no interval enclosure and no "
    "certified bound anywhere in this pipeline.",
)


@dataclass(frozen=True)
class SnapshotRecord:
    """One checkpointed snapshot and its independent elliptic cross-check.

    ``psi_cross_solver_relative_difference`` is
    ``psi_cross_solver_max_abs_difference / psi_cross_solver_relative_denominator``
    and that denominator is ``max |psi1|`` of solver A *at this snapshot*, i.e.
    the same number as ``max_abs_psi1``.
    """

    time: float
    checkpoint: str
    psi_cross_solver_max_abs_difference: float
    psi_cross_solver_relative_difference: float
    psi_cross_solver_relative_denominator: float
    psi_cross_solver_argmax_r: float
    psi_cross_solver_argmax_z: float
    solver_b_boundary_error_max: float
    solver_b_physical_defect_max: float
    max_abs_psi1: float
    max_abs_u1: float
    argmax_u1_r: float
    argmax_u1_z: float
    max_cartesian_vorticity: float
    argmax_cartesian_vorticity_r: float
    argmax_cartesian_vorticity_z: float


@dataclass(frozen=True)
class ResolutionResult:
    """Outcome of one resolution, including an explicit failure slot."""

    grid: AxisymmetricGrid
    completed: bool
    failure: str | None
    step_count: int
    final_time: float
    history: tuple[dict[str, float], ...]
    snapshots: tuple[SnapshotRecord, ...]
    viscosity_switch_times: tuple[float, ...]
    initial_max_abs_u1: float
    initial_max_cartesian_vorticity: float
    initial_norm_relative_errors: dict[str, float]
    state: ConstrainedState | None = dataclass_field(default=None)


def _sha256(path_or_bytes: Path | bytes) -> str:
    if isinstance(path_or_bytes, bytes):
        return hashlib.sha256(path_or_bytes).hexdigest()
    digest = hashlib.sha256()
    with path_or_bytes.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _prepare_output(output_dir: Path) -> None:
    if output_dir.exists():
        if not output_dir.is_dir():
            raise FileExistsError(f"output path is not a directory: {output_dir}")
        if any(output_dir.iterdir()):
            raise FileExistsError(
                f"refusing to overwrite nonempty output directory: {output_dir}"
            )
    else:
        output_dir.mkdir(parents=True)


def _inside_repository(path: Path) -> bool:
    try:
        path.relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def _number(value: object, *, name: str, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "positive and " if positive else ""
        raise ValueError(f"{name} must be {qualifier}finite")
    return result


def _positive_integer(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def validate_config(config: dict[str, Any]) -> None:
    """Reject any config with missing, unknown or out-of-range entries."""

    required = {
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "coordinate_system",
        "normalization_description",
        "units",
        "seed",
        "r_max",
        "z_period",
        "resolutions",
        "amplitude_scale",
        "viscosity_schedule",
        "t_final",
        "snapshot_times",
        "cfl_coefficient",
        "max_time_step",
        "max_steps",
        "diagnostic_stride",
        "acceptance",
    }
    # Optional keys keep every shipped config valid without edits; each one has
    # a documented default.
    optional = {"cfl_excess_tolerance"}
    if (
        not isinstance(config, dict)
        or not required <= set(config)
        or not set(config) <= required | optional
    ):
        raise ValueError("hou early time config has missing or unknown keys")
    if config["schema_version"] != 1:
        raise ValueError("schema_version must be 1")
    if config["experiment_id"] != "hou_early_time_v1":
        raise ValueError("experiment_id must be hou_early_time_v1")
    for key in (
        "description",
        "interpretation",
        "coordinate_system",
        "normalization_description",
    ):
        if not isinstance(config[key], str) or not config[key].strip():
            raise ValueError(f"{key} must be a nonempty string")
    if not isinstance(config["units"], dict) or not config["units"]:
        raise ValueError("units must be a nonempty object")
    if isinstance(config["seed"], bool) or not isinstance(config["seed"], int):
        raise ValueError("seed must be an integer")
    _number(config["r_max"], name="r_max", positive=True)
    _number(config["z_period"], name="z_period", positive=True)

    resolutions = config["resolutions"]
    if not isinstance(resolutions, list) or len(resolutions) < 2:
        raise ValueError("at least two resolutions are required")
    for item in resolutions:
        if (
            not isinstance(item, list)
            or len(item) != 2
            or any(isinstance(value, bool) or not isinstance(value, int) for value in item)
        ):
            raise ValueError("each resolution must be an [nr, nz] integer pair")
        if item[0] < 4 or item[1] < 5:
            raise ValueError("resolutions are below AxisymmetricGrid minima")
    if any(a[0] >= b[0] for a, b in zip(resolutions, resolutions[1:])):
        raise ValueError("radial resolutions must be strictly increasing")

    _number(config["amplitude_scale"], name="amplitude_scale", positive=True)
    normalize_viscosity_schedule(config["viscosity_schedule"])
    t_final = _number(config["t_final"], name="t_final", positive=True)

    snapshots = config["snapshot_times"]
    if not isinstance(snapshots, list) or not snapshots:
        raise ValueError("snapshot_times must be a nonempty list")
    previous = -math.inf
    for value in snapshots:
        moment = _number(value, name="snapshot_times entry")
        if moment < 0.0 or moment > t_final:
            raise ValueError("snapshot times must lie in [0, t_final]")
        if moment <= previous:
            raise ValueError("snapshot times must be strictly increasing")
        previous = moment
    if abs(snapshots[-1] - t_final) > 1.0e-12:
        raise ValueError("the last snapshot time must be t_final")

    _number(config["cfl_coefficient"], name="cfl_coefficient", positive=True)
    if "cfl_excess_tolerance" in config:
        tolerance = _number(
            config["cfl_excess_tolerance"], name="cfl_excess_tolerance"
        )
        if tolerance < 0.0:
            raise ValueError("cfl_excess_tolerance must be nonnegative")
    _number(config["max_time_step"], name="max_time_step", positive=True)
    _positive_integer(config["max_steps"], name="max_steps")
    _positive_integer(config["diagnostic_stride"], name="diagnostic_stride")

    acceptance = config["acceptance"]
    expected_acceptance = {
        "maximum_circulation_growth_ratio",
        "maximum_energy_growth_ratio",
        "maximum_initial_norm_relative_error",
        "maximum_odd_symmetry_defect_ratio",
    }
    if not isinstance(acceptance, dict) or set(acceptance) != expected_acceptance:
        raise ValueError("acceptance object has missing or unknown keys")
    for key in expected_acceptance:
        _number(acceptance[key], name=f"acceptance.{key}", positive=True)


def build_grid(config: dict[str, Any], *, nr: int, nz: int) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=float(config["r_max"]),
        z_min=0.0,
        z_max=float(config["z_period"]),
        periodic_z=True,
    )


def initial_norms(
    grid: AxisymmetricGrid,
    *,
    amplitude_scale: float,
) -> dict[str, float]:
    """Measure the E-29b initial norms on ``grid`` and their relative errors."""

    u1 = hou_initial_swirl(
        grid,
        amplitude=E29_AMPLITUDE,
        amplitude_scale=amplitude_scale,
    )
    state = constrain_state(grid, u1, np.zeros(grid.shape))
    magnitude = np.sqrt(
        sum(
            component * component
            for component in cartesian_vorticity(grid, state.u1, state.omega1)
        )
    )
    measured_u1 = float(np.max(np.abs(u1)))
    measured_vorticity = float(np.max(magnitude))
    expected_u1 = E29B_MAX_ABS_U1 * float(amplitude_scale)
    expected_vorticity = E29B_MAX_CARTESIAN_VORTICITY * float(amplitude_scale)
    return {
        "max_abs_u1": measured_u1,
        "max_cartesian_vorticity": measured_vorticity,
        "expected_max_abs_u1": expected_u1,
        "expected_max_cartesian_vorticity": expected_vorticity,
        "max_abs_u1_relative_error": abs(measured_u1 - expected_u1) / expected_u1,
        "max_cartesian_vorticity_relative_error": (
            abs(measured_vorticity - expected_vorticity) / expected_vorticity
        ),
    }


def _cross_check_snapshot(
    grid: AxisymmetricGrid,
    state: ConstrainedState,
) -> dict[str, float]:
    """Re-solve the snapshot elliptic problem with the independent solver B.

    The relative difference is divided by ``max |psi1|`` of solver A at this
    snapshot; that denominator and the ``(r, z)`` location of the largest
    absolute difference are returned alongside it so the report never shows a
    relative number without the scale it was taken against.
    """

    solution = solve_finite_cylinder_poisson(
        state.omega1,
        grid,
        outer_boundary=0.0,
        condition_mode_indices=(),
    )
    deviation = np.abs(state.psi1 - solution.psi)
    row, column = divmod(int(np.argmax(deviation)), grid.nz)
    difference = float(deviation[row, column])
    denominator = float(np.max(np.abs(state.psi1)))
    scale = max(1.0e-300, denominator)
    return {
        "max_abs_difference": difference,
        "relative_difference": difference / scale,
        "relative_denominator": denominator,
        "argmax_r": float(grid.r[row]),
        "argmax_z": float(grid.z[column]),
        "boundary_error_max": float(solution.diagnostics.boundary_error_max),
        "physical_defect_max": float(
            solution.diagnostics.physical_cross_stencil_defect_max
        ),
    }


def evolve_resolution(
    config: dict[str, Any],
    *,
    nr: int,
    nz: int,
    checkpoint_dir: Path | None = None,
    provenance: dict[str, Any] | None = None,
) -> ResolutionResult:
    """Integrate one resolution, checkpointing and cross-checking snapshots."""

    validate_config(config)
    grid = build_grid(config, nr=nr, nz=nz)
    amplitude_scale = float(config["amplitude_scale"])
    schedule = normalize_viscosity_schedule(config["viscosity_schedule"])
    u1_initial = hou_initial_swirl(
        grid,
        amplitude=E29_AMPLITUDE,
        amplitude_scale=amplitude_scale,
    )
    norms = initial_norms(grid, amplitude_scale=amplitude_scale)
    snapshots: list[SnapshotRecord] = []

    def on_snapshot(
        time: float,
        state: ConstrainedState,
        diagnostics: dict[str, float],
    ) -> None:
        cross = _cross_check_snapshot(grid, state)
        difference = cross["max_abs_difference"]
        name = ""
        if checkpoint_dir is not None:
            name = f"checkpoint_nr{nr}_nz{nz}_t{len(snapshots):03d}.npz"
            save_checkpoint(
                checkpoint_dir / name,
                grid=grid,
                state=state,
                time=time,
                viscosity=viscosity_at(schedule, time),
                seed=int(config["seed"]),
                config=config,
                provenance=provenance,
                metadata={
                    "experiment_id": config["experiment_id"],
                    "nr": nr,
                    "nz": nz,
                    "amplitude_scale": amplitude_scale,
                    "cross_solver_psi_max_abs_difference": difference,
                },
            )
        snapshots.append(
            SnapshotRecord(
                time=float(time),
                checkpoint=name,
                psi_cross_solver_max_abs_difference=difference,
                psi_cross_solver_relative_difference=cross["relative_difference"],
                psi_cross_solver_relative_denominator=cross["relative_denominator"],
                psi_cross_solver_argmax_r=cross["argmax_r"],
                psi_cross_solver_argmax_z=cross["argmax_z"],
                solver_b_boundary_error_max=cross["boundary_error_max"],
                solver_b_physical_defect_max=cross["physical_defect_max"],
                max_abs_psi1=float(np.max(np.abs(state.psi1))),
                max_abs_u1=float(diagnostics["max_abs_u1"]),
                argmax_u1_r=float(diagnostics["argmax_u1_r"]),
                argmax_u1_z=float(diagnostics["argmax_u1_z"]),
                max_cartesian_vorticity=float(
                    diagnostics["max_cartesian_vorticity"]
                ),
                argmax_cartesian_vorticity_r=float(
                    diagnostics["argmax_cartesian_vorticity_r"]
                ),
                argmax_cartesian_vorticity_z=float(
                    diagnostics["argmax_cartesian_vorticity_z"]
                ),
            )
        )

    failure: str | None = None
    try:
        result = integrate(
            grid,
            u1=u1_initial,
            omega1=np.zeros(grid.shape),
            t_final=float(config["t_final"]),
            viscosity_schedule=schedule,
            cfl_coefficient=float(config["cfl_coefficient"]),
            max_time_step=float(config["max_time_step"]),
            max_steps=int(config["max_steps"]),
            snapshot_times=[float(value) for value in config["snapshot_times"]],
            diagnostic_stride=int(config["diagnostic_stride"]),
            on_snapshot=on_snapshot,
        )
    except (FloatingPointError, ArithmeticError, ValueError) as exc:
        return ResolutionResult(
            grid=grid,
            completed=False,
            failure=f"{type(exc).__name__}: {exc}",
            step_count=0,
            final_time=0.0,
            history=(),
            snapshots=tuple(snapshots),
            viscosity_switch_times=(),
            initial_max_abs_u1=norms["max_abs_u1"],
            initial_max_cartesian_vorticity=norms["max_cartesian_vorticity"],
            initial_norm_relative_errors={
                "max_abs_u1": norms["max_abs_u1_relative_error"],
                "max_cartesian_vorticity": norms[
                    "max_cartesian_vorticity_relative_error"
                ],
            },
        )
    if not result.completed:
        failure = (
            f"integration stopped after {result.step_count} steps at "
            f"t={result.time!r} before t_final"
        )
    return ResolutionResult(
        grid=grid,
        completed=bool(result.completed),
        failure=failure,
        step_count=result.step_count,
        final_time=float(result.time),
        history=result.history,
        snapshots=tuple(snapshots),
        viscosity_switch_times=result.viscosity_switch_times,
        initial_max_abs_u1=norms["max_abs_u1"],
        initial_max_cartesian_vorticity=norms["max_cartesian_vorticity"],
        initial_norm_relative_errors={
            "max_abs_u1": norms["max_abs_u1_relative_error"],
            "max_cartesian_vorticity": norms[
                "max_cartesian_vorticity_relative_error"
            ],
        },
        state=result.state,
    )


def _series(result: ResolutionResult, name: str) -> list[float]:
    return [record[name] for record in result.history]


# The record-only relativized constraint diagnostics, persisted per resolution
# as a time series next to the absolute quantities they relativize.
CONSTRAINT_SERIES_FIELDS: tuple[str, ...] = (
    "time",
    "divergence_residual_max",
    *RELATIVE_DIAGNOSTIC_FIELDS,
)

DENOMINATOR_DEFINITIONS = {
    "divergence_residual_relative": (
        "max over the grid of |d_r u^r| + |u^r/r| + |d_z u^z| (key "
        "divergence_relative_denominator), built from the same three terms "
        "whose sum is the E-02 residual, with the same regular axis limit "
        "u^r/r -> d_r u^r at r = 0"
    ),
    "divergence_pointwise_ratio_max": (
        "the pointwise term sum |d_r u^r| + |u^r/r| + |d_z u^z| at each grid "
        "point; the reported value is the largest pointwise ratio, i.e. the "
        "reciprocal of the smallest TM-09 cancellation index kappa_res"
    ),
    "axis_parity_relative_u1": (
        "max |d_r u1| over the grid (key axis_parity_relative_denominator_u1)"
    ),
    "axis_parity_relative_omega1": (
        "max |d_r omega1| over the grid (key "
        "axis_parity_relative_denominator_omega1)"
    ),
    "snapshot_cross_solver_psi_relative_difference": (
        "max |psi1| of solver A at the same snapshot (snapshots.csv column "
        "psi_cross_solver_relative_denominator, equal to max_abs_psi1)"
    ),
}


def _relativization_row(record: dict[str, float]) -> dict[str, float]:
    """Extract the relativized constraint entries of one diagnostic record."""

    return {name: float(record[name]) for name in CONSTRAINT_SERIES_FIELDS}


def _cross_solver_relativization(
    snapshots: tuple[SnapshotRecord, ...],
) -> dict[str, Any] | None:
    """Report the worst relative cross-solver snapshot with its denominator."""

    if not snapshots:
        return None
    worst = max(snapshots, key=lambda item: item.psi_cross_solver_relative_difference)
    return {
        "relative_difference_denominator": DENOMINATOR_DEFINITIONS[
            "snapshot_cross_solver_psi_relative_difference"
        ],
        "worst_relative_snapshot_time": worst.time,
        "worst_relative_difference": worst.psi_cross_solver_relative_difference,
        "worst_relative_max_abs_difference": (
            worst.psi_cross_solver_max_abs_difference
        ),
        "worst_relative_denominator": worst.psi_cross_solver_relative_denominator,
        "worst_relative_argmax_r": worst.psi_cross_solver_argmax_r,
        "worst_relative_argmax_z": worst.psi_cross_solver_argmax_z,
    }


def resolution_metrics(
    result: ResolutionResult,
    *,
    cfl_coefficient: float | None = None,
) -> dict[str, Any]:
    """Reduce one resolution's history to the recorded acceptance quantities."""

    identity = {
        "completed": result.completed,
        "failure": result.failure,
        "step_count": result.step_count,
        "final_time": result.final_time,
        "initial_max_abs_u1_measured": result.initial_max_abs_u1,
        "initial_max_cartesian_vorticity_measured": (
            result.initial_max_cartesian_vorticity
        ),
        "initial_norm_relative_errors": result.initial_norm_relative_errors,
    }
    if not result.history:
        return identity
    energy = _series(result, "energy")
    circulation = _series(result, "circulation_max")
    swirl = _series(result, "max_abs_u1")
    amplification = _series(result, "amplification")
    odd_defect = _series(result, "odd_symmetry_defect")
    finite = all(
        math.isfinite(value)
        for record in result.history
        for value in record.values()
    )
    maximum_cfl = max(_series(result, "advective_cfl"))
    coefficient = (
        None
        if cfl_coefficient is None
        else _number(cfl_coefficient, name="cfl_coefficient", positive=True)
    )
    return {
        **identity,
        "all_diagnostics_finite": finite,
        "initial_energy": energy[0],
        "final_energy": energy[-1],
        "maximum_energy_growth_ratio": max(energy) / energy[0] - 1.0,
        "initial_circulation_max": circulation[0],
        "final_circulation_max": circulation[-1],
        "maximum_circulation_growth_ratio": max(circulation) / circulation[0] - 1.0,
        "initial_max_abs_u1": swirl[0],
        "final_max_abs_u1": swirl[-1],
        "minimum_max_abs_u1": min(swirl),
        "early_swirl_maximum_decreases": bool(len(swirl) > 1 and swirl[1] < swirl[0]),
        "final_amplification": amplification[-1],
        "maximum_amplification": max(amplification),
        "final_argmax_u1_r": result.history[-1]["argmax_u1_r"],
        "final_argmax_u1_z": result.history[-1]["argmax_u1_z"],
        "final_argmax_cartesian_vorticity_r": result.history[-1][
            "argmax_cartesian_vorticity_r"
        ],
        "final_argmax_cartesian_vorticity_z": result.history[-1][
            "argmax_cartesian_vorticity_z"
        ],
        "maximum_odd_symmetry_defect": max(odd_defect),
        "maximum_odd_symmetry_defect_ratio": max(odd_defect) / max(swirl),
        "maximum_axis_parity_defect": max(_series(result, "axis_parity_defect")),
        "maximum_divergence_residual": max(
            _series(result, "divergence_residual_max")
        ),
        # ---- record-only relativized constraint diagnostics (TM-09) --------
        "maximum_divergence_residual_relative": max(
            _series(result, "divergence_residual_relative")
        ),
        "maximum_divergence_pointwise_ratio": max(
            _series(result, "divergence_pointwise_ratio_max")
        ),
        "maximum_axis_parity_relative_u1": max(
            _series(result, "axis_parity_relative_u1")
        ),
        "maximum_axis_parity_relative_omega1": max(
            _series(result, "axis_parity_relative_omega1")
        ),
        "constraint_relativization": {
            "gated": False,
            "denominator_definitions": dict(DENOMINATOR_DEFINITIONS),
            "initial": _relativization_row(result.history[0]),
            "final": _relativization_row(result.history[-1]),
        },
        "constraint_relativization_series": {
            name: _series(result, name) for name in CONSTRAINT_SERIES_FIELDS
        },
        "cross_solver_relativization": _cross_solver_relativization(
            result.snapshots
        ),
        "maximum_wall_u1_abs": max(_series(result, "wall_u1_max_abs")),
        "minimum_dt": min(record["dt"] for record in result.history[1:])
        if len(result.history) > 1
        else 0.0,
        "maximum_advective_cfl": maximum_cfl,
        # Post-step measurement against the start-of-step selection rule; see
        # CFL_EXCESS_MECHANISM and the summary cfl_policy block.
        "maximum_advective_cfl_excess_ratio": (
            None if coefficient is None else maximum_cfl / coefficient - 1.0
        ),
        "viscosity_switch_times": list(result.viscosity_switch_times),
        "snapshot_cross_solver_psi_max_abs_difference": (
            max(
                item.psi_cross_solver_max_abs_difference for item in result.snapshots
            )
            if result.snapshots
            else None
        ),
        "snapshot_cross_solver_psi_relative_difference": (
            max(
                item.psi_cross_solver_relative_difference for item in result.snapshots
            )
            if result.snapshots
            else None
        ),
    }


def evaluate(
    config: dict[str, Any],
    *,
    checkpoint_dir: Path | None = None,
    provenance: dict[str, Any] | None = None,
) -> tuple[list[ResolutionResult], list[dict[str, Any]], dict[str, bool]]:
    validate_config(config)
    acceptance = config["acceptance"]
    cfl_coefficient = float(config["cfl_coefficient"])
    cfl_excess_tolerance = cfl_excess_tolerance_of(config)
    results = [
        evolve_resolution(
            config,
            nr=int(nr),
            nz=int(nz),
            checkpoint_dir=checkpoint_dir,
            provenance=provenance,
        )
        for nr, nz in config["resolutions"]
    ]
    metrics = [
        resolution_metrics(result, cfl_coefficient=cfl_coefficient)
        for result in results
    ]
    completed = [item for item in metrics if item.get("all_diagnostics_finite") is not None]
    finest = metrics[-1]
    checks = {
        "all_resolutions_completed": all(item["completed"] for item in metrics),
        "all_diagnostics_finite": bool(completed)
        and all(item["all_diagnostics_finite"] for item in completed),
        "energy_non_increasing": bool(completed)
        and all(
            item["maximum_energy_growth_ratio"]
            <= float(acceptance["maximum_energy_growth_ratio"])
            for item in completed
        ),
        "circulation_max_principle": bool(completed)
        and all(
            item["maximum_circulation_growth_ratio"]
            <= float(acceptance["maximum_circulation_growth_ratio"])
            for item in completed
        ),
        "early_swirl_maximum_decreases": bool(completed)
        and all(item["early_swirl_maximum_decreases"] for item in completed),
        "odd_symmetry_preserved": bool(completed)
        and all(
            item["maximum_odd_symmetry_defect_ratio"]
            <= float(acceptance["maximum_odd_symmetry_defect_ratio"])
            for item in completed
        ),
        "initial_norms_match_e29b": max(
            finest["initial_norm_relative_errors"]["max_abs_u1"],
            finest["initial_norm_relative_errors"]["max_cartesian_vorticity"],
        )
        <= float(acceptance["maximum_initial_norm_relative_error"]),
        "cross_solver_elliptic_agreement_recorded": all(
            item.get("snapshot_cross_solver_psi_max_abs_difference") is not None
            and math.isfinite(
                float(item["snapshot_cross_solver_psi_max_abs_difference"])
            )
            for item in completed
        ),
        # The recorded CFL is a post-step measurement of a step sized at the
        # start of the step, so it is accepted against C (1 + tolerance).
        "advective_cfl_within_tolerance": bool(completed)
        and all(
            advective_cfl_within_tolerance(
                item["maximum_advective_cfl"],
                cfl_coefficient=cfl_coefficient,
                cfl_excess_tolerance=cfl_excess_tolerance,
            )
            for item in completed
        ),
    }
    return results, metrics, checks


def cfl_excess_tolerance_of(config: dict[str, Any]) -> float:
    """Return the configured CFL excess tolerance or its documented default."""

    if "cfl_excess_tolerance" not in config:
        return DEFAULT_CFL_EXCESS_TOLERANCE
    return _number(config["cfl_excess_tolerance"], name="cfl_excess_tolerance")


def cfl_policy(
    config: dict[str, Any],
    metrics: list[dict[str, Any]],
) -> dict[str, Any]:
    """Return the summary block that fixes the CFL acceptance semantics."""

    coefficient = float(config["cfl_coefficient"])
    tolerance = cfl_excess_tolerance_of(config)
    measured = [
        float(item["maximum_advective_cfl"])
        for item in metrics
        if item.get("maximum_advective_cfl") is not None
    ]
    worst = max(measured) if measured else None
    return {
        "rule": CFL_RULE,
        "effective_cfl_definition": EFFECTIVE_CFL_DEFINITION,
        "why_the_effective_cfl_can_exceed_the_coefficient": CFL_EXCESS_MECHANISM,
        "cfl_coefficient": coefficient,
        "cfl_excess_tolerance": tolerance,
        "cfl_excess_tolerance_source": (
            "config key cfl_excess_tolerance"
            if "cfl_excess_tolerance" in config
            else f"default {DEFAULT_CFL_EXCESS_TOLERANCE}"
        ),
        "accepted_effective_cfl_bound": coefficient * (1.0 + tolerance),
        "maximum_effective_cfl": worst,
        "maximum_effective_cfl_excess_ratio": (
            None if worst is None else worst / coefficient - 1.0
        ),
        "intra_step_growth_guidance": (
            "docs/hou_setup_audit.md records the delegated HH21 "
            "(arXiv:2102.06663) prescription: advective and viscous CFL "
            "numbers 0.1 with a one-step relative growth of the maximum below "
            "one percent, which bounds the excess above"
        ),
        "v1_reference_measurement": {
            "run": "outputs/hou_early_time_v1 (193x384 resolution)",
            "cfl_coefficient": V1_CFL_COEFFICIENT,
            "maximum_advective_cfl": V1_MAXIMUM_ADVECTIVE_CFL,
            "excess_ratio": V1_ADVECTIVE_CFL_EXCESS_RATIO,
            "excess_percent": 100.0 * V1_ADVECTIVE_CFL_EXCESS_RATIO,
        },
    }


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Run every resolution and write a checksummed evidence bundle."""

    validate_config(config)
    provenance = collect_runtime_provenance()
    config_bytes = (
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    _prepare_output(output_dir)
    (output_dir / "config.snapshot.json").write_bytes(config_bytes)
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir()
    results, metrics, checks = evaluate(
        config,
        checkpoint_dir=checkpoint_dir,
        provenance=provenance,
    )

    diagnostic_rows = [
        {"nr": result.grid.nr, "nz": result.grid.nz, **record}
        for result in results
        for record in result.history
    ]
    if not diagnostic_rows:
        raise RuntimeError("no diagnostics were recorded")
    _write_csv(
        output_dir / "diagnostics.csv",
        ["nr", "nz", *DIAGNOSTIC_FIELDS],
        diagnostic_rows,
    )

    snapshot_rows = [
        {
            "nr": result.grid.nr,
            "nz": result.grid.nz,
            "time": item.time,
            "checkpoint": item.checkpoint,
            "max_abs_psi1": item.max_abs_psi1,
            "psi_cross_solver_max_abs_difference": (
                item.psi_cross_solver_max_abs_difference
            ),
            "psi_cross_solver_relative_difference": (
                item.psi_cross_solver_relative_difference
            ),
            "psi_cross_solver_relative_denominator": (
                item.psi_cross_solver_relative_denominator
            ),
            "psi_cross_solver_argmax_r": item.psi_cross_solver_argmax_r,
            "psi_cross_solver_argmax_z": item.psi_cross_solver_argmax_z,
            "solver_b_boundary_error_max": item.solver_b_boundary_error_max,
            "solver_b_physical_defect_max": item.solver_b_physical_defect_max,
            "max_abs_u1": item.max_abs_u1,
            "argmax_u1_r": item.argmax_u1_r,
            "argmax_u1_z": item.argmax_u1_z,
            "max_cartesian_vorticity": item.max_cartesian_vorticity,
            "argmax_cartesian_vorticity_r": item.argmax_cartesian_vorticity_r,
            "argmax_cartesian_vorticity_z": item.argmax_cartesian_vorticity_z,
        }
        for result in results
        for item in result.snapshots
    ]
    if snapshot_rows:
        _write_csv(
            output_dir / "snapshots.csv",
            list(snapshot_rows[0]),
            snapshot_rows,
        )

    arrays: dict[str, np.ndarray] = {}
    for result in results:
        suffix = f"nr{result.grid.nr}_nz{result.grid.nz}"
        arrays[f"r_{suffix}"] = result.grid.r
        arrays[f"z_{suffix}"] = result.grid.z
        if result.state is not None:
            arrays[f"u1_final_{suffix}"] = result.state.u1
            arrays[f"omega1_final_{suffix}"] = result.state.omega1
            arrays[f"psi1_final_{suffix}"] = result.state.psi1
        for name in ("time", "max_abs_u1", "max_cartesian_vorticity", "amplification"):
            arrays[f"{name}_{suffix}"] = np.asarray(
                _series(result, name) if result.history else [],
                dtype=np.float64,
            )
    np.savez_compressed(output_dir / "trajectories.npz", **arrays)

    summary: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Early-time numerical observation of the E-29 initial datum on "
            "uniform fixed grids. It is not a reproduction of any published "
            "amplification value, not a singularity candidate, and not a proof."
        ),
        "method": {
            "equations": "E-11, E-12, E-13, E-14",
            "time_integrator": "explicit Heun/RK2",
            "time_step_rule": (
                "adaptive min(C dr/max|u^r|, C dz/max|u^z|, "
                "C min(dr,dz)^2/(4 nu)) capped by max_time_step"
            ),
            "wall_conditions": "E-27 with the second-order Thom form E-31",
            "poisson_solver": (
                "solver A: Fourier-z / r^3-flux finite volume with homogeneous "
                "outer Dirichlet data"
            ),
            "independent_cross_check": (
                "solver B: finite_cylinder_poisson non-divergence radial "
                "stencil, re-solved at every snapshot"
            ),
            "viscosity_protocol": "E-30 two-stage piecewise-constant schedule",
            "initial_data": "E-29 with amplitude 12000 times amplitude_scale",
            "symmetry": "full z period, odd symmetry monitored not imposed",
            "constraint_relativization": (
                "the E-02 divergence residual and the E-16c axis parity defect "
                "are reported both absolutely and divided by their own "
                "cancellation scales (TM-09); the relative numbers are "
                "recorded and no acceptance check reads them"
            ),
        },
        "cfl_policy": cfl_policy(config, metrics),
        "reproducibility": {
            "seed": int(config["seed"]),
            "config_sha256": _sha256(config_bytes),
            "runtime_provenance": provenance,
            "command": (
                "python experiments/run_hou_early_time.py --config "
                "configs/hou_early_time.json --output-dir outputs/hou_early_time"
            ),
        },
        "derived_reference_norms": {
            "source": "E-29b (derived by this repository, absent from the paper)",
            "max_abs_u1": E29B_MAX_ABS_U1,
            "max_cartesian_vorticity": E29B_MAX_CARTESIAN_VORTICITY,
        },
        "resolutions": [
            {"nr": result.grid.nr, "nz": result.grid.nz, "dr": result.grid.dr,
             "dz": result.grid.dz, **metric}
            for result, metric in zip(results, metrics)
        ],
        "acceptance_checks": checks,
        "accepted_as_early_time_observation": bool(all(checks.values())),
        "limitations": list(LIMITATIONS),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    artifact_names = [
        "config.snapshot.json",
        "diagnostics.csv",
        "summary.json",
        "trajectories.npz",
    ]
    if snapshot_rows:
        artifact_names.append("snapshots.csv")
    artifact_names.extend(
        sorted(
            path.relative_to(output_dir).as_posix()
            for path in checkpoint_dir.rglob("*")
            if path.is_file()
        )
    )
    manifest = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "files": {
            name: {
                "sha256": _sha256(output_dir / name),
                "bytes": (output_dir / name).stat().st_size,
            }
            for name in sorted(artifact_names)
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
        default=REPOSITORY_ROOT / "configs" / "hou_early_time.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "hou_early_time",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="hou early time config",
    )
    if not isinstance(config, dict):
        parser.error("config must be a JSON object")
    summary = run(config, output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["accepted_as_early_time_observation"] else 2


if __name__ == "__main__":
    sys.exit(main())
