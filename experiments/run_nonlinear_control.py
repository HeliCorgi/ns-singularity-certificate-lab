#!/usr/bin/env python
"""Run a forced smooth nonlinear Poisson/time-integration negative control.

At every Heun right-hand-side evaluation this experiment solves

    -(d_rr + 3 d_r/r + d_zz) psi1 = omega1

with the independent Fourier/finite-volume solver, recovers meridional
velocity, and evaluates both transformed nonlinear equations.  The analytic
``manufactured_fields`` trajectory supplies a time-dependent body force and
all outer traces.  Consequently this is deliberately *not* an unforced
general Navier--Stokes calculation, a singularity candidate, or a proof.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Literal

import numpy as np

from ns_certificate_lab._integrity import strict_json_loads, write_with_digest
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.manufactured import manufactured_fields
from ns_certificate_lab.operators import (
    derivative_r,
    derivative_z,
    laplacian_5d_formal,
    recover_velocity,
)
from ns_certificate_lab.poisson import PoissonSolution, solve_streamfunction_poisson
from ns_certificate_lab.provenance import collect_runtime_provenance


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
Fault = Literal[
    "poisson_sign",
    "omega_swirl_sign",
    "velocity_z_sign",
    "outer_trace",
]
VALID_FAULTS: tuple[Fault, ...] = (
    "poisson_sign",
    "omega_swirl_sign",
    "velocity_z_sign",
    "outer_trace",
)


@dataclass(frozen=True)
class FlowSnapshot:
    time: float
    energy: float
    max_vorticity: float
    exact_energy: float
    exact_max_vorticity: float
    poisson_discrete_max: float
    poisson_pde_rms: float
    poisson_pde_max: float
    solver_boundary_defect: float
    analytic_boundary_trace_error: float


@dataclass(frozen=True)
class ResolutionResult:
    grid: AxisymmetricGrid
    step_count: int
    actual_dt: float
    u1: np.ndarray
    omega1: np.ndarray
    psi1: np.ndarray
    exact_u1: np.ndarray
    exact_omega1: np.ndarray
    exact_psi1: np.ndarray
    relative_errors: dict[str, float]
    snapshots: tuple[FlowSnapshot, ...]
    maximum_rhs_discrete_poisson_residual: float
    maximum_rhs_pde_poisson_residual: float
    maximum_rhs_solver_boundary_defect: float
    outer_u1_trace_error: float
    outer_omega1_trace_error: float
    poisson_metadata: dict[str, Any]


@dataclass
class _ResidualAccumulator:
    discrete_max: float = 0.0
    pde_max: float = 0.0
    boundary_max: float = 0.0

    def add(self, solution: PoissonSolution) -> None:
        self.discrete_max = max(
            self.discrete_max,
            float(np.max(np.abs(solution.discrete_residual[:-1]))),
        )
        self.pde_max = max(
            self.pde_max,
            float(np.max(np.abs(solution.pde_residual[:-1]))),
        )
        self.boundary_max = max(
            self.boundary_max,
            float(solution.metadata["outer_boundary_max_abs_defect"]),
        )


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


def validate_config(config: dict[str, Any]) -> None:
    required = {
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "coordinate_system",
        "normalization_description",
        "units",
        "seed",
        "viscosity",
        "final_time",
        "r_max",
        "z_min",
        "z_max",
        "radial_resolutions",
        "axial_resolutions",
        "time_step_coefficient",
        "diagnostic_stride",
        "fault_amplitude",
        "acceptance",
    }
    if set(config) != required:
        raise ValueError("nonlinear control config has missing or unknown keys")
    if config["schema_version"] != 1:
        raise ValueError("schema_version must be 1")
    for key in (
        "experiment_id",
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
    viscosity = _number(config["viscosity"], name="viscosity")
    if viscosity < 0.0:
        raise ValueError("viscosity must be nonnegative")
    _number(config["final_time"], name="final_time", positive=True)
    _number(config["r_max"], name="r_max", positive=True)
    z_min = _number(config["z_min"], name="z_min")
    z_max = _number(config["z_max"], name="z_max")
    if z_max <= z_min:
        raise ValueError("z_max must exceed z_min")
    radial = config["radial_resolutions"]
    axial = config["axial_resolutions"]
    if (
        not isinstance(radial, list)
        or not isinstance(axial, list)
        or len(radial) != 3
        or len(axial) != 3
        or any(isinstance(value, bool) or not isinstance(value, int) for value in radial + axial)
    ):
        raise ValueError("exactly three integer radial and axial resolutions are required")
    if any(value < 4 for value in radial) or any(value < 5 for value in axial):
        raise ValueError("resolutions are below AxisymmetricGrid minima")
    if any(a >= b for a, b in zip(radial, radial[1:])) or any(
        a >= b for a, b in zip(axial, axial[1:])
    ):
        raise ValueError("resolutions must be strictly increasing")
    _number(
        config["time_step_coefficient"],
        name="time_step_coefficient",
        positive=True,
    )
    if (
        isinstance(config["diagnostic_stride"], bool)
        or not isinstance(config["diagnostic_stride"], int)
        or config["diagnostic_stride"] < 1
    ):
        raise ValueError("diagnostic_stride must be a positive integer")
    _number(config["fault_amplitude"], name="fault_amplitude", positive=True)
    acceptance = config["acceptance"]
    expected_acceptance = {
        "minimum_observed_order",
        "maximum_finest_relative_error",
        "maximum_discrete_poisson_residual",
        "maximum_outer_trace_error",
        "maximum_energy_relative_error",
        "maximum_vorticity_relative_error",
    }
    if not isinstance(acceptance, dict) or set(acceptance) != expected_acceptance:
        raise ValueError("acceptance object has missing or unknown keys")
    for key in expected_acceptance:
        _number(acceptance[key], name=f"acceptance.{key}", positive=True)


def _analytic_outer_traces(
    grid: AxisymmetricGrid,
    time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fields = manufactured_fields(grid, time=time)
    return fields.u1[-1].copy(), fields.omega1[-1].copy(), fields.psi1[-1].copy()


def _poisson(
    grid: AxisymmetricGrid,
    omega1: np.ndarray,
    time: float,
    *,
    fault: Fault | None,
    fault_amplitude: float,
    estimate_condition: bool,
) -> tuple[PoissonSolution, np.ndarray, np.ndarray]:
    _, _, exact_trace = _analytic_outer_traces(grid, time)
    supplied_trace = exact_trace.copy()
    if fault == "outer_trace":
        supplied_trace += fault_amplitude * np.cos(grid.z)
    solution = solve_streamfunction_poisson(
        grid,
        omega1,
        supplied_trace,
        estimate_condition=estimate_condition,
    )
    psi1 = solution.psi1.copy()
    if fault == "poisson_sign":
        psi1 *= -1.0
    return solution, psi1, exact_trace


def _right_hand_side(
    grid: AxisymmetricGrid,
    *,
    u1: np.ndarray,
    omega1: np.ndarray,
    time: float,
    viscosity: float,
    fault: Fault | None,
    fault_amplitude: float,
    residuals: _ResidualAccumulator,
) -> tuple[np.ndarray, np.ndarray]:
    solution, psi1, _ = _poisson(
        grid,
        omega1,
        time,
        fault=fault,
        fault_amplitude=fault_amplitude,
        estimate_condition=False,
    )
    residuals.add(solution)
    u_r, u_z = recover_velocity(grid, psi1)
    if fault == "velocity_z_sign":
        u_z = -u_z
    u1_r = derivative_r(grid, u1, even_at_axis=True)
    u1_z = derivative_z(grid, u1)
    omega1_r = derivative_r(grid, omega1, even_at_axis=True)
    omega1_z = derivative_z(grid, omega1)
    exact = manufactured_fields(grid, time=time)
    swirl_source = derivative_z(grid, u1 * u1)
    if fault == "omega_swirl_sign":
        swirl_source = -swirl_source
    rhs_u1 = (
        -u_r * u1_r
        - u_z * u1_z
        + 2.0 * u1 * derivative_z(grid, psi1)
        + viscosity * laplacian_5d_formal(grid, u1)
        + exact.forcing_u1(viscosity)
    )
    rhs_omega1 = (
        -u_r * omega1_r
        - u_z * omega1_z
        + swirl_source
        + viscosity * laplacian_5d_formal(grid, omega1)
        + exact.forcing_omega1(viscosity)
    )
    # These are prescribed time-dependent outer state traces, not numerical
    # PDE rows.  Their explicit use is part of the finite-cylinder problem.
    rhs_u1[-1] = exact.du1_dt[-1]
    rhs_omega1[-1] = exact.domega1_dt[-1]
    if not (np.all(np.isfinite(rhs_u1)) and np.all(np.isfinite(rhs_omega1))):
        raise FloatingPointError("nonlinear control RHS became non-finite")
    return rhs_u1, rhs_omega1


def _weighted_squared_norm(grid: AxisymmetricGrid, values: np.ndarray) -> float:
    array = np.asarray(values, dtype=np.float64)
    if array.shape[:2] != grid.shape:
        raise ValueError("weighted norm field has incompatible leading shape")
    radial_weights = grid.r.copy()
    radial_weights[0] *= 0.5
    radial_weights[-1] *= 0.5
    squared = array * array
    if array.ndim == 3:
        squared = np.sum(squared, axis=2)
    return float(grid.dr * grid.dz * np.sum(radial_weights[:, None] * squared))


def _relative_error(
    grid: AxisymmetricGrid,
    numerical: np.ndarray,
    exact: np.ndarray,
) -> float:
    denominator = _weighted_squared_norm(grid, exact)
    if denominator <= 0.0:
        raise ValueError("relative-error reference norm must be positive")
    return math.sqrt(
        _weighted_squared_norm(grid, numerical - exact) / denominator
    )


def _velocity_and_vorticity(
    grid: AxisymmetricGrid,
    *,
    u1: np.ndarray,
    omega1: np.ndarray,
    psi1: np.ndarray,
    velocity_z_fault: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    u_r, u_z = recover_velocity(grid, psi1)
    if velocity_z_fault:
        u_z = -u_z
    radius = grid.r[:, None]
    velocity = np.stack((u_r, radius * u1, u_z), axis=2)
    vorticity = np.stack(
        (
            -radius * derivative_z(grid, u1),
            radius * omega1,
            2.0 * u1 + radius * derivative_r(grid, u1, even_at_axis=True),
        ),
        axis=2,
    )
    return velocity, vorticity


def _flow_snapshot(
    grid: AxisymmetricGrid,
    *,
    u1: np.ndarray,
    omega1: np.ndarray,
    time: float,
    fault: Fault | None,
    fault_amplitude: float,
) -> tuple[FlowSnapshot, PoissonSolution, np.ndarray]:
    solution, psi1, exact_trace = _poisson(
        grid,
        omega1,
        time,
        fault=fault,
        fault_amplitude=fault_amplitude,
        estimate_condition=False,
    )
    exact = manufactured_fields(grid, time=time)
    velocity, vorticity = _velocity_and_vorticity(
        grid,
        u1=u1,
        omega1=omega1,
        psi1=psi1,
        velocity_z_fault=fault == "velocity_z_sign",
    )
    exact_velocity, exact_vorticity = _velocity_and_vorticity(
        grid,
        u1=exact.u1,
        omega1=exact.omega1,
        psi1=exact.psi1,
    )
    energy = math.pi * _weighted_squared_norm(grid, velocity)
    exact_energy = math.pi * _weighted_squared_norm(grid, exact_velocity)
    return (
        FlowSnapshot(
            time=float(time),
            energy=energy,
            max_vorticity=float(np.max(np.linalg.norm(vorticity, axis=2))),
            exact_energy=exact_energy,
            exact_max_vorticity=float(
                np.max(np.linalg.norm(exact_vorticity, axis=2))
            ),
            poisson_discrete_max=float(
                np.max(np.abs(solution.discrete_residual[:-1]))
            ),
            poisson_pde_rms=float(
                np.sqrt(np.mean(solution.pde_residual[:-1] ** 2))
            ),
            poisson_pde_max=float(np.max(np.abs(solution.pde_residual[:-1]))),
            solver_boundary_defect=float(
                solution.metadata["outer_boundary_max_abs_defect"]
            ),
            analytic_boundary_trace_error=float(
                np.max(np.abs(psi1[-1] - exact_trace))
            ),
        ),
        solution,
        psi1,
    )


def evolve_resolution(
    config: dict[str, Any],
    *,
    nr: int,
    nz: int,
    fault: Fault | None = None,
) -> ResolutionResult:
    """Evolve one grid.  ``fault`` exists only for explicit rejection tests."""

    validate_config(config)
    if fault is not None and fault not in VALID_FAULTS:
        raise ValueError(f"unknown injected fault: {fault}")
    grid = AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=float(config["r_max"]),
        z_min=float(config["z_min"]),
        z_max=float(config["z_max"]),
        periodic_z=True,
    )
    final_time = float(config["final_time"])
    target_dt = float(config["time_step_coefficient"]) * grid.dr**2
    step_count = max(1, math.ceil(final_time / target_dt))
    dt = final_time / step_count
    viscosity = float(config["viscosity"])
    fault_amplitude = float(config["fault_amplitude"])
    diagnostic_stride = int(config["diagnostic_stride"])
    initial = manufactured_fields(grid, time=0.0)
    u1 = initial.u1.copy()
    omega1 = initial.omega1.copy()
    residuals = _ResidualAccumulator()
    snapshots: list[FlowSnapshot] = []
    snapshot, _, _ = _flow_snapshot(
        grid,
        u1=u1,
        omega1=omega1,
        time=0.0,
        fault=fault,
        fault_amplitude=fault_amplitude,
    )
    snapshots.append(snapshot)

    for step in range(step_count):
        time = step * dt
        k1_u1, k1_omega1 = _right_hand_side(
            grid,
            u1=u1,
            omega1=omega1,
            time=time,
            viscosity=viscosity,
            fault=fault,
            fault_amplitude=fault_amplitude,
            residuals=residuals,
        )
        stage_u1 = u1 + dt * k1_u1
        stage_omega1 = omega1 + dt * k1_omega1
        next_time = (step + 1) * dt
        trace_u1, trace_omega1, _ = _analytic_outer_traces(grid, next_time)
        stage_u1[-1] = trace_u1
        stage_omega1[-1] = trace_omega1
        k2_u1, k2_omega1 = _right_hand_side(
            grid,
            u1=stage_u1,
            omega1=stage_omega1,
            time=next_time,
            viscosity=viscosity,
            fault=fault,
            fault_amplitude=fault_amplitude,
            residuals=residuals,
        )
        u1 = u1 + 0.5 * dt * (k1_u1 + k2_u1)
        omega1 = omega1 + 0.5 * dt * (k1_omega1 + k2_omega1)
        u1[-1] = trace_u1
        omega1[-1] = trace_omega1
        if not (np.all(np.isfinite(u1)) and np.all(np.isfinite(omega1))):
            raise FloatingPointError("nonlinear control state became non-finite")
        if (step + 1) % diagnostic_stride == 0 or step + 1 == step_count:
            snapshot, _, _ = _flow_snapshot(
                grid,
                u1=u1,
                omega1=omega1,
                time=next_time,
                fault=fault,
                fault_amplitude=fault_amplitude,
            )
            snapshots.append(snapshot)

    final_exact = manufactured_fields(grid, time=final_time)
    final_solution, psi1, exact_trace = _poisson(
        grid,
        omega1,
        final_time,
        fault=fault,
        fault_amplitude=fault_amplitude,
        estimate_condition=True,
    )
    numerical_velocity, _ = _velocity_and_vorticity(
        grid,
        u1=u1,
        omega1=omega1,
        psi1=psi1,
        velocity_z_fault=fault == "velocity_z_sign",
    )
    exact_velocity, _ = _velocity_and_vorticity(
        grid,
        u1=final_exact.u1,
        omega1=final_exact.omega1,
        psi1=final_exact.psi1,
    )
    relative_errors = {
        "u1": _relative_error(grid, u1, final_exact.u1),
        "omega1": _relative_error(grid, omega1, final_exact.omega1),
        "psi1": _relative_error(grid, psi1, final_exact.psi1),
        "velocity": _relative_error(grid, numerical_velocity, exact_velocity),
    }
    trace_u1, trace_omega1, _ = _analytic_outer_traces(grid, final_time)
    if not all(math.isfinite(value) for value in relative_errors.values()):
        raise FloatingPointError("nonlinear control produced non-finite errors")
    # The final solve is also part of the recorded elliptic audit.
    residuals.add(final_solution)
    metadata = dict(final_solution.metadata)
    metadata["analytic_outer_trace_max_abs_error"] = float(
        np.max(np.abs(psi1[-1] - exact_trace))
    )
    return ResolutionResult(
        grid=grid,
        step_count=step_count,
        actual_dt=dt,
        u1=u1,
        omega1=omega1,
        psi1=psi1,
        exact_u1=final_exact.u1,
        exact_omega1=final_exact.omega1,
        exact_psi1=final_exact.psi1,
        relative_errors=relative_errors,
        snapshots=tuple(snapshots),
        maximum_rhs_discrete_poisson_residual=residuals.discrete_max,
        maximum_rhs_pde_poisson_residual=residuals.pde_max,
        maximum_rhs_solver_boundary_defect=residuals.boundary_max,
        outer_u1_trace_error=float(np.max(np.abs(u1[-1] - trace_u1))),
        outer_omega1_trace_error=float(
            np.max(np.abs(omega1[-1] - trace_omega1))
        ),
        poisson_metadata=metadata,
    )


def observed_orders(
    errors: list[float],
    spacings: list[float],
) -> list[float]:
    if len(errors) != len(spacings) or len(errors) < 2:
        raise ValueError("errors and spacings must have matching length >= 2")
    if any(
        not math.isfinite(value) or value <= 0.0
        for value in errors + spacings
    ):
        raise ValueError("errors and spacings must be positive and finite")
    if any(a <= b for a, b in zip(spacings, spacings[1:])):
        raise ValueError("spacings must be strictly decreasing")
    return [
        math.log(errors[index] / errors[index + 1])
        / math.log(spacings[index] / spacings[index + 1])
        for index in range(len(errors) - 1)
    ]


def evaluate(
    config: dict[str, Any],
    *,
    fault: Fault | None = None,
) -> tuple[list[ResolutionResult], dict[str, list[float]], dict[str, bool]]:
    validate_config(config)
    results = [
        evolve_resolution(config, nr=nr, nz=nz, fault=fault)
        for nr, nz in zip(
            config["radial_resolutions"],
            config["axial_resolutions"],
        )
    ]
    spacings = [result.grid.dr for result in results]
    orders = {
        name: observed_orders(
            [result.relative_errors[name] for result in results],
            spacings,
        )
        for name in ("u1", "omega1", "psi1", "velocity")
    }
    acceptance = config["acceptance"]
    maximum_error = float(acceptance["maximum_finest_relative_error"])
    final_snapshot = results[-1].snapshots[-1]
    energy_relative_error = abs(
        final_snapshot.energy - final_snapshot.exact_energy
    ) / final_snapshot.exact_energy
    vorticity_relative_error = abs(
        final_snapshot.max_vorticity - final_snapshot.exact_max_vorticity
    ) / final_snapshot.exact_max_vorticity
    checks = {
        "all_field_errors_decrease": all(
            all(fine < coarse for coarse, fine in zip(
                [result.relative_errors[name] for result in results],
                [result.relative_errors[name] for result in results][1:],
            ))
            for name in orders
        ),
        "minimum_observed_order": min(
            order for values in orders.values() for order in values
        )
        >= float(acceptance["minimum_observed_order"]),
        "finest_field_errors": max(results[-1].relative_errors.values())
        <= maximum_error,
        "discrete_poisson_residual": max(
            result.maximum_rhs_discrete_poisson_residual for result in results
        )
        <= float(acceptance["maximum_discrete_poisson_residual"]),
        "outer_trace": max(
            max(
                result.maximum_rhs_solver_boundary_defect,
                result.outer_u1_trace_error,
                result.outer_omega1_trace_error,
                float(
                    result.poisson_metadata[
                        "analytic_outer_trace_max_abs_error"
                    ]
                ),
            )
            for result in results
        )
        <= float(acceptance["maximum_outer_trace_error"]),
        "energy_matches_oracle": energy_relative_error
        <= float(acceptance["maximum_energy_relative_error"]),
        "maximum_vorticity_matches_oracle": vorticity_relative_error
        <= float(acceptance["maximum_vorticity_relative_error"]),
    }
    return results, orders, checks


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Run the preregistered control and write checksummed evidence."""

    validate_config(config)
    provenance = collect_runtime_provenance()
    config_bytes = (
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    _prepare_output(output_dir)
    (output_dir / "config.snapshot.json").write_bytes(config_bytes)
    results, orders, checks = evaluate(config)

    convergence_rows: list[dict[str, Any]] = []
    for index, result in enumerate(results):
        row: dict[str, Any] = {
            "nr": result.grid.nr,
            "nz": result.grid.nz,
            "dr": result.grid.dr,
            "dz": result.grid.dz,
            "step_count": result.step_count,
            "actual_dt": result.actual_dt,
            "u1_relative_l2_error": result.relative_errors["u1"],
            "omega1_relative_l2_error": result.relative_errors["omega1"],
            "psi1_relative_l2_error": result.relative_errors["psi1"],
            "velocity_relative_l2_error": result.relative_errors["velocity"],
            "maximum_rhs_discrete_poisson_residual": (
                result.maximum_rhs_discrete_poisson_residual
            ),
            "maximum_rhs_independent_pde_poisson_residual": (
                result.maximum_rhs_pde_poisson_residual
            ),
            "solver_boundary_defect": result.maximum_rhs_solver_boundary_defect,
            "analytic_psi_outer_trace_error": result.poisson_metadata[
                "analytic_outer_trace_max_abs_error"
            ],
            "u1_outer_trace_error": result.outer_u1_trace_error,
            "omega1_outer_trace_error": result.outer_omega1_trace_error,
            "final_energy": result.snapshots[-1].energy,
            "final_exact_energy": result.snapshots[-1].exact_energy,
            "maximum_vorticity_over_samples": max(
                item.max_vorticity for item in result.snapshots
            ),
            "maximum_exact_vorticity_over_samples": max(
                item.exact_max_vorticity for item in result.snapshots
            ),
        }
        for name in ("u1", "omega1", "psi1", "velocity"):
            row[f"{name}_order_from_previous"] = (
                "" if index == 0 else orders[name][index - 1]
            )
        convergence_rows.append(row)
    convergence_fields = list(convergence_rows[0])
    _write_csv(output_dir / "convergence.csv", convergence_fields, convergence_rows)

    diagnostic_rows = [
        {
            "nr": result.grid.nr,
            "nz": result.grid.nz,
            "time": snapshot.time,
            "energy": snapshot.energy,
            "exact_energy": snapshot.exact_energy,
            "max_vorticity": snapshot.max_vorticity,
            "exact_max_vorticity": snapshot.exact_max_vorticity,
            "poisson_discrete_max": snapshot.poisson_discrete_max,
            "poisson_independent_pde_rms": snapshot.poisson_pde_rms,
            "poisson_independent_pde_max": snapshot.poisson_pde_max,
            "solver_boundary_defect": snapshot.solver_boundary_defect,
            "analytic_boundary_trace_error": snapshot.analytic_boundary_trace_error,
        }
        for result in results
        for snapshot in result.snapshots
    ]
    _write_csv(
        output_dir / "diagnostics.csv",
        list(diagnostic_rows[0]),
        diagnostic_rows,
    )
    arrays: dict[str, np.ndarray] = {}
    for result in results:
        suffix = f"nr{result.grid.nr}_nz{result.grid.nz}"
        arrays[f"r_{suffix}"] = result.grid.r
        arrays[f"z_{suffix}"] = result.grid.z
        arrays[f"u1_final_{suffix}"] = result.u1
        arrays[f"omega1_final_{suffix}"] = result.omega1
        arrays[f"psi1_final_{suffix}"] = result.psi1
        arrays[f"u1_exact_{suffix}"] = result.exact_u1
        arrays[f"omega1_exact_{suffix}"] = result.exact_omega1
        arrays[f"psi1_exact_{suffix}"] = result.exact_psi1
        arrays[f"times_{suffix}"] = np.asarray(
            [item.time for item in result.snapshots],
            dtype=np.float64,
        )
        arrays[f"energy_{suffix}"] = np.asarray(
            [item.energy for item in result.snapshots],
            dtype=np.float64,
        )
        arrays[f"max_vorticity_{suffix}"] = np.asarray(
            [item.max_vorticity for item in result.snapshots],
            dtype=np.float64,
        )
    np.savez_compressed(output_dir / "trajectories.npz", **arrays)

    final_snapshot = results[-1].snapshots[-1]
    summary: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Forced smooth manufactured negative control. It is not an "
            "unforced general Navier-Stokes solution, a singularity candidate, "
            "or a proof."
        ),
        "method": {
            "time_integrator": "explicit Heun/RK2",
            "time_step_rule": "dt <= coefficient * dr^2",
            "poisson_equation": "-L5 psi1 = omega1",
            "poisson_solver": (
                "independent Fourier-z / r^3-flux finite-volume solver"
            ),
            "poisson_called_at_every_rhs": True,
            "outer_boundary": (
                "explicit analytic time-dependent Dirichlet traces for "
                "psi1, u1, and omega1"
            ),
            "forcing": (
                "analytic manufactured forcing evaluated at every RK stage"
            ),
        },
        "reproducibility": {
            "seed": int(config["seed"]),
            "config_sha256": _sha256(config_bytes),
            "runtime_provenance": provenance,
        },
        "resolutions": [
            {
                "nr": result.grid.nr,
                "nz": result.grid.nz,
                "dr": result.grid.dr,
                "dz": result.grid.dz,
                "step_count": result.step_count,
                "actual_dt": result.actual_dt,
                "relative_errors": result.relative_errors,
                "maximum_rhs_discrete_poisson_residual": (
                    result.maximum_rhs_discrete_poisson_residual
                ),
                "maximum_rhs_independent_pde_poisson_residual": (
                    result.maximum_rhs_pde_poisson_residual
                ),
                "poisson_metadata": result.poisson_metadata,
                "final_energy": result.snapshots[-1].energy,
                "final_exact_energy": result.snapshots[-1].exact_energy,
                "maximum_vorticity_over_samples": max(
                    item.max_vorticity for item in result.snapshots
                ),
                "maximum_exact_vorticity_over_samples": max(
                    item.exact_max_vorticity for item in result.snapshots
                ),
            }
            for result in results
        ],
        "observed_orders": orders,
        "finest_energy_relative_error": abs(
            final_snapshot.energy - final_snapshot.exact_energy
        )
        / final_snapshot.exact_energy,
        "finest_max_vorticity_relative_error": abs(
            final_snapshot.max_vorticity - final_snapshot.exact_max_vorticity
        )
        / final_snapshot.exact_max_vorticity,
        "acceptance_checks": checks,
        "accepted_as_forced_smooth_control": bool(all(checks.values())),
        "limitations": [
            "The forcing is chosen from the analytic oracle.",
            "Exact finite-cylinder outer traces are prescribed.",
            "The Poisson solve is finite-dimensional and has no interval bound.",
            "This does not test an unforced general solution or singularity.",
        ],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    artifact_names = [
        "config.snapshot.json",
        "convergence.csv",
        "diagnostics.csv",
        "summary.json",
        "trajectories.npz",
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
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False)
            + "\n"
        ).encode("utf-8"),
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "nonlinear_manufactured_control.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "nonlinear_manufactured_control",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="nonlinear manufactured control config",
    )
    if not isinstance(config, dict):
        parser.error("config must be a JSON object")
    summary = run(config, output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["accepted_as_forced_smooth_control"] else 2


if __name__ == "__main__":
    sys.exit(main())
