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
    """One checkpointed snapshot and its independent elliptic cross-check."""

    time: float
    checkpoint: str
    psi_cross_solver_max_abs_difference: float
    psi_cross_solver_relative_difference: float
    solver_b_boundary_error_max: float
    solver_b_physical_defect_max: float
    max_abs_psi1: float


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
    if not isinstance(config, dict) or set(config) != required:
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
) -> tuple[float, float, float, float]:
    """Re-solve the snapshot elliptic problem with the independent solver B."""

    solution = solve_finite_cylinder_poisson(
        state.omega1,
        grid,
        outer_boundary=0.0,
        condition_mode_indices=(),
    )
    difference = float(np.max(np.abs(state.psi1 - solution.psi)))
    scale = max(1.0e-300, float(np.max(np.abs(state.psi1))))
    return (
        difference,
        difference / scale,
        float(solution.diagnostics.boundary_error_max),
        float(solution.diagnostics.physical_cross_stencil_defect_max),
    )


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
        _diagnostics: dict[str, float],
    ) -> None:
        difference, relative, boundary, defect = _cross_check_snapshot(grid, state)
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
                psi_cross_solver_relative_difference=relative,
                solver_b_boundary_error_max=boundary,
                solver_b_physical_defect_max=defect,
                max_abs_psi1=float(np.max(np.abs(state.psi1))),
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


def resolution_metrics(result: ResolutionResult) -> dict[str, Any]:
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
        "maximum_odd_symmetry_defect": max(odd_defect),
        "maximum_odd_symmetry_defect_ratio": max(odd_defect) / max(swirl),
        "maximum_axis_parity_defect": max(_series(result, "axis_parity_defect")),
        "maximum_divergence_residual": max(
            _series(result, "divergence_residual_max")
        ),
        "maximum_wall_u1_abs": max(_series(result, "wall_u1_max_abs")),
        "minimum_dt": min(record["dt"] for record in result.history[1:])
        if len(result.history) > 1
        else 0.0,
        "maximum_advective_cfl": max(_series(result, "advective_cfl")),
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
    metrics = [resolution_metrics(result) for result in results]
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
    }
    return results, metrics, checks


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
            "solver_b_boundary_error_max": item.solver_b_boundary_error_max,
            "solver_b_physical_defect_max": item.solver_b_physical_defect_max,
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
        },
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
