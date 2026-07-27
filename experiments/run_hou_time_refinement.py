#!/usr/bin/env python
r"""Fixed-grid time-step refinement study for the Hou early-time run.

One spatial grid, one final time, three *fixed* time steps ``dt``, ``dt/2``,
``dt/4``.  Everything that :mod:`experiments.run_hou_early_time` varies with
resolution is held frozen here, so the only quantity that changes between the
three members is the temporal discretization.  Comparing the members at the
common final time therefore separates the temporal error from the spatial
error that the resolution ladder of ``outputs/hou_early_time_v1`` and
``outputs/hou_early_time_v2_hires`` already measured.

Compared at ``t_final`` (user requirement, task 3)
--------------------------------------------------

* ``||omega||_inf`` amplification (and the raw ``max |omega|`` it comes from),
* ``||u_1||_inf`` and the grid location of its argmax,
* the kinetic energy E-20a,
* the circulation max-principle violation E-23,
* the *relative* E-02 divergence residual and the *relative* E-16c axis parity
  defect (TM-09 relativization, both recorded not gated),
* the independent solver-B elliptic cross-check difference at ``t_final``.

Smooth scalars are analysed by step doubling: ``|Q(dt) - Q(dt/2)|``,
``|Q(dt/2) - Q(dt/4)|`` and the observed order ``log2`` of their ratio.  Heun
is second order, so an order near two is the expected outcome.  When the two
differences fall below the explicit measurement floor
``acceptance.temporal_error_floor_relative * |Q(dt/4)|`` the order is *not*
reported as a number and the quantity is recorded as
``temporal_error_below_measurement_floor``; that outcome is an accepted
result with an explanation, not a failure, because a temporal error that
small is exactly the statement this experiment is trying to establish.

The constraint metrics (relative divergence, relative axis parity, circulation
violation, cross-solver difference) are expected to be *spatially* dominated
and therefore roughly independent of ``dt``.  Confirming that expectation is
itself half of the temporal/spatial separation statement, so their three
values and their spread across the ladder are recorded explicitly in
``temporal_vs_spatial_separation`` rather than being gated.

Time-step policy
----------------

``nonlinear_cylinder.integrate`` runs in its uniform ``fixed_time_step`` mode.
It re-derives the step count from the requested step with a ceiling and then
uses the *uniform* step ``(t_final - t_0) / step_count``, so the run lands
exactly on ``t_final`` without a short final step; the tiny clipping the
integrator still applies to the last step (accumulated floating-point drift of
``t + dt``) is measured and recorded per level.  When ``t_final`` is not an
integer multiple of the requested ``dt`` the effective step is therefore
slightly *smaller* than requested, which is recorded as
``requested_dt_divides_t_final = false`` together with the relative
adjustment.  The refinement ladder stays exact as long as the three step
counts are ``n``, ``2n``, ``4n``; that is validated and gated.

No adaptive step is taken.  The CFL numbers are *monitored*: because the
adaptive advective bound is ``C / max(|u^r|/dr, |u^z|/dz)`` and the recorded
``advective_cfl`` diagnostic is ``dt * max(|u^r|/dr, |u^z|/dz)``, the fixed
step satisfies the adaptive advective bound exactly when the recorded
effective CFL stays at or below ``C``; the same equivalence holds for the
viscous bound and ``viscous_cfl``.  Both are accepted against
``C * (1 + cfl_excess_tolerance)`` with the same documented 5 percent
tolerance as ``experiments/run_hou_early_time.py``.  An exceedance is
recorded as a failed acceptance check, never raised.

Nothing here is a proof, an interval enclosure, or evidence for singularity
formation.  A second-order-in-time convergence rate on a fixed uniform grid
says nothing about the continuum limit of the spatial discretization.
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

from experiments.run_hou_early_time import (
    DEFAULT_CFL_EXCESS_TOLERANCE,
    E29_AMPLITUDE,
    E29B_MAX_ABS_U1,
    E29B_MAX_CARTESIAN_VORTICITY,
    advective_cfl_within_tolerance,
    build_grid,
    initial_norms,
)
from ns_certificate_lab._integrity import strict_json_loads, write_with_digest
from ns_certificate_lab.finite_cylinder_poisson import solve_finite_cylinder_poisson
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    DIAGNOSTIC_FIELDS,
    RELATIVE_DIAGNOSTIC_FIELDS,
    ConstrainedState,
    hou_initial_swirl,
    integrate,
    normalize_viscosity_schedule,
    save_checkpoint,
    viscosity_at,
)
from ns_certificate_lab.provenance import collect_runtime_provenance


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPERIMENT_ID = "hou_time_refinement_v1"

# Exactly three members: dt, dt/2, dt/4.
EXPECTED_TIME_STEP_COUNT = 3
# Relative tolerance for "each entry is exactly half of the previous one".
HALVING_RELATIVE_TOLERANCE = 1.0e-12
# Relative tolerance for "t_final is an integer multiple of this dt".
DIVISIBILITY_RELATIVE_TOLERANCE = 1.0e-9
# ``integrate`` re-derives its own step count as ceil((t_final - t_0) / dt).
# Handing it exactly ``t_final / n`` can round to ``n + 1`` (for example
# ``2e-5 / (2e-5 / 10) == 10.000000000000002``), which would silently break the
# halving ladder.  Nudging the requested step up by one part in 1e9 makes that
# ceiling reproduce the step count validated here for every n < 1e9.  The
# uniform step the integrator then uses is ``(t_final - t_0) / step_count``,
# i.e. the recorded effective dt, and is unaffected by the nudge.
STEP_COUNT_NUDGE = 1.0 + 1.0e-9
# Heun/RK2: the leading temporal error term is O(dt^2).
NOMINAL_TEMPORAL_ORDER = 2.0

# The smooth scalars whose temporal order is measured by step doubling.
# ``final_amplification`` and ``final_max_cartesian_vorticity`` are exactly
# proportional on a fixed grid (identical reference norm), so only the former
# is analysed; ``final_max_abs_omega1`` is the independent vorticity scalar.
SMOOTH_QUANTITIES: tuple[tuple[str, str], ...] = (
    ("final_amplification", "amplification"),
    ("final_max_abs_u1", "max_abs_u1"),
    ("final_energy", "energy"),
    ("final_max_abs_omega1", "max_abs_omega1"),
)

# The metrics that are expected to be spatially dominated, i.e. essentially
# independent of dt.  Recorded with their spread across the ladder; not gated
# on a trend (the circulation growth ratio is gated separately against the
# preregistered E-23 bound).
CONSTRAINT_QUANTITIES: tuple[str, ...] = (
    "final_divergence_residual_relative",
    "final_divergence_pointwise_ratio_max",
    "final_axis_parity_relative_u1",
    "final_axis_parity_relative_omega1",
    "maximum_circulation_growth_ratio",
    "final_cross_solver_psi_relative_difference",
)

# The argmax locations compared across the ladder.  They live on the shared
# discrete grid, so exact equality is the likely outcome; agreement is
# reported to within one grid cell and never forced.
ARGMAX_QUANTITIES: tuple[tuple[str, str, str], ...] = (
    ("argmax_u1", "final_argmax_u1_r", "final_argmax_u1_z"),
    (
        "argmax_cartesian_vorticity",
        "final_argmax_cartesian_vorticity_r",
        "final_argmax_cartesian_vorticity_z",
    ),
)

# ---------------------------------------------------------------- spatial ladder
#
# The already-shipped spatial refinement of the identical E-29/E-30 early-time
# run at the identical t_final = T_1.  These four amplification values are read
# back from the two summaries at run time (see
# :func:`read_spatial_amplification_ladder`); the constants below are the
# audited copies used to detect drift and to keep this module usable when the
# outputs directory is absent.
SPATIAL_AMPLIFICATION_LADDER: tuple[dict[str, Any], ...] = (
    {
        "nr": 65,
        "nz": 128,
        "final_amplification": 6.114768458177527,
        "summary": "outputs/hou_early_time_v1/summary.json",
    },
    {
        "nr": 129,
        "nz": 256,
        "final_amplification": 12.695695243738214,
        "summary": "outputs/hou_early_time_v1/summary.json",
    },
    {
        "nr": 193,
        "nz": 384,
        "final_amplification": 15.627954940635405,
        "summary": "outputs/hou_early_time_v1/summary.json",
    },
    {
        "nr": 257,
        "nz": 512,
        "final_amplification": 17.258778340188396,
        "summary": "outputs/hou_early_time_v2_hires/summary.json",
    },
)

SPATIAL_LADDER_SUMMARIES: tuple[str, ...] = (
    "outputs/hou_early_time_v1/summary.json",
    "outputs/hou_early_time_v2_hires/summary.json",
)

LIMITATIONS = (
    "This measures the temporal discretization error of the shipped Heun "
    "integrator on one fixed uniform grid. It says nothing about the "
    "continuum limit of the spatial discretization.",
    "A uniform fixed grid cannot resolve the adaptive mesh scales (minimum "
    "spacing of order 1e-8) used by the source calculation; every quantity "
    "reported here is resolution limited.",
    "This is a numerical observation, not a reproduction claim and not "
    "evidence of singularity formation.",
    "Step doubling estimates the temporal error of the coarser member; it is "
    "not a bound on the total error, which is dominated by the fixed spatial "
    "discretization at this resolution.",
    "Solver A and solver B share the same second-order accuracy class, so "
    "their agreement bounds discretization consistency, not continuum error.",
    "Floating-point arithmetic only: there is no interval enclosure and no "
    "certified bound anywhere in this pipeline.",
)

# The cross-solver helper mirrors ``_cross_check_snapshot`` in
# ``experiments/run_hou_early_time.py``.  That function is private there and
# was not promoted (the module is finalized), so it is reimplemented here with
# identical semantics; see :func:`cross_check_state`.
KNOWN_GAPS = (
    "The solver-B cross-check helper duplicates the private "
    "_cross_check_snapshot of experiments/run_hou_early_time.py because that "
    "helper is not part of the module's public surface; the two must be kept "
    "in step by hand.",
)


# ------------------------------------------------------------------- helpers


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


def cfl_within_tolerance(
    value: float,
    *,
    cfl_coefficient: float,
    cfl_excess_tolerance: float,
) -> bool:
    """Return whether a recorded CFL number is inside ``C (1 + tolerance)``.

    This is the directional generalization of
    :func:`experiments.run_hou_early_time.advective_cfl_within_tolerance` to
    the viscous channel; the advective channel keeps calling the shipped
    function so that both experiments share one policy implementation.
    """

    coefficient = _number(cfl_coefficient, name="cfl_coefficient", positive=True)
    tolerance = _number(cfl_excess_tolerance, name="cfl_excess_tolerance")
    if tolerance < 0.0:
        raise ValueError("cfl_excess_tolerance must be nonnegative")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    measured = float(value)
    if not math.isfinite(measured):
        return False
    return measured <= coefficient * (1.0 + tolerance)


# ------------------------------------------------------------ configuration


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
        "resolution",
        "amplitude_scale",
        "viscosity_schedule",
        "t_final",
        "time_steps",
        "cfl_coefficient",
        "max_steps",
        "diagnostic_stride",
        "acceptance",
    }
    optional = {"cfl_excess_tolerance"}
    if (
        not isinstance(config, dict)
        or not required <= set(config)
        or not set(config) <= required | optional
    ):
        raise ValueError("hou time refinement config has missing or unknown keys")
    if config["schema_version"] != 1:
        raise ValueError("schema_version must be 1")
    if config["experiment_id"] != EXPERIMENT_ID:
        raise ValueError(f"experiment_id must be {EXPERIMENT_ID}")
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

    resolution = config["resolution"]
    if (
        not isinstance(resolution, list)
        or len(resolution) != 2
        or any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in resolution
        )
    ):
        raise ValueError("resolution must be a single [nr, nz] integer pair")
    if resolution[0] < 4 or resolution[1] < 5:
        raise ValueError("resolution is below AxisymmetricGrid minima")

    _number(config["amplitude_scale"], name="amplitude_scale", positive=True)
    schedule = normalize_viscosity_schedule(config["viscosity_schedule"])
    t_final = _number(config["t_final"], name="t_final", positive=True)
    if len(schedule) > 1 and schedule[1][0] < t_final:
        raise ValueError(
            "the time refinement study requires a single viscosity stage over "
            "[0, t_final]: no switch time may fall below t_final"
        )

    validate_time_steps(config["time_steps"])

    _number(config["cfl_coefficient"], name="cfl_coefficient", positive=True)
    if "cfl_excess_tolerance" in config:
        tolerance = _number(
            config["cfl_excess_tolerance"], name="cfl_excess_tolerance"
        )
        if tolerance < 0.0:
            raise ValueError("cfl_excess_tolerance must be nonnegative")
    _positive_integer(config["max_steps"], name="max_steps")
    _positive_integer(config["diagnostic_stride"], name="diagnostic_stride")

    acceptance = config["acceptance"]
    expected_acceptance = {
        "maximum_circulation_growth_ratio",
        "maximum_energy_growth_ratio",
        "temporal_order_lower_bound",
        "temporal_order_upper_bound",
        "minimum_quantities_at_second_order",
        "temporal_error_floor_relative",
    }
    if not isinstance(acceptance, dict) or set(acceptance) != expected_acceptance:
        raise ValueError("acceptance object has missing or unknown keys")
    for key in (
        "maximum_circulation_growth_ratio",
        "maximum_energy_growth_ratio",
        "temporal_order_lower_bound",
        "temporal_order_upper_bound",
        "temporal_error_floor_relative",
    ):
        _number(acceptance[key], name=f"acceptance.{key}", positive=True)
    lower = float(acceptance["temporal_order_lower_bound"])
    upper = float(acceptance["temporal_order_upper_bound"])
    if not lower < upper:
        raise ValueError(
            "acceptance.temporal_order_lower_bound must be below "
            "acceptance.temporal_order_upper_bound"
        )
    minimum = _positive_integer(
        acceptance["minimum_quantities_at_second_order"],
        name="acceptance.minimum_quantities_at_second_order",
    )
    if minimum > len(SMOOTH_QUANTITIES):
        raise ValueError(
            "acceptance.minimum_quantities_at_second_order exceeds the number "
            "of analysed smooth quantities"
        )


def validate_time_steps(time_steps: object) -> tuple[float, ...]:
    """Return the validated ``(dt, dt/2, dt/4)`` ladder.

    Exactly three entries, strictly decreasing, each one *exactly* half of the
    previous entry to within :data:`HALVING_RELATIVE_TOLERANCE`.
    """

    if (
        not isinstance(time_steps, list)
        or len(time_steps) != EXPECTED_TIME_STEP_COUNT
    ):
        raise ValueError(
            f"time_steps must be a list of {EXPECTED_TIME_STEP_COUNT} entries"
        )
    values = [
        _number(item, name="time_steps entry", positive=True) for item in time_steps
    ]
    for coarse, fine in zip(values, values[1:]):
        if fine >= coarse:
            raise ValueError("time_steps must be strictly decreasing")
        if abs(fine - 0.5 * coarse) > HALVING_RELATIVE_TOLERANCE * coarse:
            raise ValueError(
                "time_steps must halve exactly: every entry must be one half "
                "of the previous entry"
            )
    return tuple(values)


# ------------------------------------------------------------ step planning


@dataclass(frozen=True)
class StepPlan:
    """The uniform step actually taken for one requested ``dt``."""

    requested_dt: float
    step_count: int
    effective_dt: float
    requested_dt_divides_t_final: bool
    effective_dt_relative_adjustment: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "requested_dt": self.requested_dt,
            "step_count": self.step_count,
            "effective_dt": self.effective_dt,
            "requested_dt_divides_t_final": self.requested_dt_divides_t_final,
            "effective_dt_relative_adjustment": (
                self.effective_dt_relative_adjustment
            ),
        }


def plan_step(t_final: float, requested_dt: float) -> StepPlan:
    """Return the uniform step that lands exactly on ``t_final``.

    If ``t_final`` is an integer multiple of ``requested_dt`` (to within
    :data:`DIVISIBILITY_RELATIVE_TOLERANCE`) the requested step is used
    unchanged.  Otherwise the step count is rounded *up*, which shortens every
    step slightly instead of leaving a ragged final step; the shortening is
    recorded as ``effective_dt_relative_adjustment``.
    """

    stop = _number(t_final, name="t_final", positive=True)
    step = _number(requested_dt, name="requested_dt", positive=True)
    raw = stop / step
    nearest = float(round(raw))
    divides = nearest >= 1.0 and abs(raw - nearest) <= (
        DIVISIBILITY_RELATIVE_TOLERANCE * max(1.0, raw)
    )
    step_count = int(nearest) if divides else int(math.ceil(raw))
    if step_count < 1:
        raise ValueError("t_final is smaller than one requested time step")
    effective = stop / step_count
    return StepPlan(
        requested_dt=step,
        step_count=step_count,
        effective_dt=effective,
        requested_dt_divides_t_final=bool(divides),
        effective_dt_relative_adjustment=abs(effective - step) / step,
    )


def plan_ladder(t_final: float, time_steps: list[float]) -> tuple[StepPlan, ...]:
    """Return the three step plans for a validated ``dt, dt/2, dt/4`` ladder."""

    validate_time_steps(list(time_steps))
    return tuple(plan_step(t_final, value) for value in time_steps)


def ladder_halves_exactly(plans: tuple[StepPlan, ...]) -> bool:
    """Return whether the *effective* steps still halve exactly.

    Step counts ``n, 2n, 4n`` give effective steps ``T/n, T/(2n), T/(4n)``;
    division by two is exact in binary floating point, so the ladder is exact
    whenever the step counts double.
    """

    if len(plans) != EXPECTED_TIME_STEP_COUNT:
        return False
    base = plans[0].step_count
    if any(
        item.step_count != base * 2**index for index, item in enumerate(plans)
    ):
        return False
    return all(
        fine.effective_dt == 0.5 * coarse.effective_dt
        for coarse, fine in zip(plans, plans[1:])
    )


# -------------------------------------------------------- solver-B crosscheck


def cross_check_state(
    grid: AxisymmetricGrid,
    state: ConstrainedState,
) -> dict[str, float]:
    """Re-solve the elliptic problem of ``state`` with the independent solver B.

    Identical in semantics to the private ``_cross_check_snapshot`` of
    ``experiments/run_hou_early_time.py`` (see :data:`KNOWN_GAPS`): the
    relative difference is divided by ``max |psi1|`` of solver A at the same
    state, and both that denominator and the ``(r, z)`` location of the
    largest absolute difference are returned so no relative number is ever
    reported without the scale it was taken against.
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


# ------------------------------------------------------------------ one level


@dataclass(frozen=True)
class TimeStepRun:
    """Outcome of one fixed-``dt`` member, including an explicit failure slot."""

    index: int
    plan: StepPlan
    grid: AxisymmetricGrid
    completed: bool
    failure: str | None
    observed_step_count: int
    final_time: float
    history: tuple[dict[str, float], ...]
    cross_check: dict[str, float] | None
    checkpoint: str
    state: ConstrainedState | None = dataclass_field(default=None)


def _series(history: tuple[dict[str, float], ...], name: str) -> list[float]:
    return [record[name] for record in history]


def evolve_time_step(
    config: dict[str, Any],
    *,
    index: int,
    plan: StepPlan,
    checkpoint_dir: Path | None = None,
    provenance: dict[str, Any] | None = None,
) -> TimeStepRun:
    """Integrate one member of the ladder with a fixed uniform time step."""

    validate_config(config)
    nr, nz = (int(value) for value in config["resolution"])
    grid = build_grid(config, nr=nr, nz=nz)
    amplitude_scale = float(config["amplitude_scale"])
    schedule = normalize_viscosity_schedule(config["viscosity_schedule"])
    t_final = float(config["t_final"])
    u1_initial = hou_initial_swirl(
        grid,
        amplitude=E29_AMPLITUDE,
        amplitude_scale=amplitude_scale,
    )
    try:
        result = integrate(
            grid,
            u1=u1_initial,
            omega1=np.zeros(grid.shape),
            t_final=t_final,
            viscosity_schedule=schedule,
            fixed_time_step=plan.effective_dt * STEP_COUNT_NUDGE,
            max_steps=int(config["max_steps"]),
            diagnostic_stride=int(config["diagnostic_stride"]),
        )
    except (FloatingPointError, ArithmeticError, ValueError) as exc:
        return TimeStepRun(
            index=index,
            plan=plan,
            grid=grid,
            completed=False,
            failure=f"{type(exc).__name__}: {exc}",
            observed_step_count=0,
            final_time=0.0,
            history=(),
            cross_check=None,
            checkpoint="",
        )

    failure: str | None = None
    if not result.completed:
        failure = (
            f"integration stopped after {result.step_count} steps at "
            f"t={result.time!r} before t_final"
        )
    elif result.step_count != plan.step_count:
        failure = (
            f"integrator took {result.step_count} steps where the validated "
            f"plan requires {plan.step_count}"
        )

    cross = cross_check_state(grid, result.state)
    name = ""
    if checkpoint_dir is not None:
        name = f"checkpoint_dt{index:02d}_n{plan.step_count}.npz"
        save_checkpoint(
            checkpoint_dir / name,
            grid=grid,
            state=result.state,
            time=float(result.time),
            viscosity=viscosity_at(schedule, float(result.time)),
            seed=int(config["seed"]),
            config=config,
            provenance=provenance,
            metadata={
                "experiment_id": config["experiment_id"],
                "nr": nr,
                "nz": nz,
                "amplitude_scale": amplitude_scale,
                "requested_dt": plan.requested_dt,
                "effective_dt": plan.effective_dt,
                "step_count": plan.step_count,
                "cross_solver_psi_max_abs_difference": cross["max_abs_difference"],
            },
        )
    return TimeStepRun(
        index=index,
        plan=plan,
        grid=grid,
        completed=bool(result.completed) and failure is None,
        failure=failure,
        observed_step_count=int(result.step_count),
        final_time=float(result.time),
        history=result.history,
        cross_check=cross,
        checkpoint=name,
        state=result.state,
    )


def level_metrics(
    run_result: TimeStepRun,
    *,
    cfl_coefficient: float,
) -> dict[str, Any]:
    """Reduce one member to its recorded final-time comparison quantities."""

    identity: dict[str, Any] = {
        **run_result.plan.as_dict(),
        "level_index": run_result.index,
        "nr": run_result.grid.nr,
        "nz": run_result.grid.nz,
        "dr": run_result.grid.dr,
        "dz": run_result.grid.dz,
        "completed": run_result.completed,
        "failure": run_result.failure,
        "observed_step_count": run_result.observed_step_count,
        "final_time": run_result.final_time,
        "checkpoint": run_result.checkpoint,
    }
    if not run_result.history:
        return identity
    history = run_result.history
    final = history[-1]
    energy = _series(history, "energy")
    circulation = _series(history, "circulation_max")
    advective = _series(history, "advective_cfl")
    viscous = _series(history, "viscous_cfl")
    maximum_advective = max(advective)
    maximum_viscous = max(viscous)
    recorded_final_dt = float(final["dt"])
    cross = run_result.cross_check or {}
    return {
        **identity,
        "diagnostic_record_count": len(history),
        "all_diagnostics_finite": all(
            math.isfinite(value) for record in history for value in record.values()
        ),
        # The integrator clips the very last step so that the run lands on
        # t_final despite the accumulated drift of t + dt.
        "final_recorded_dt": recorded_final_dt,
        "final_step_relative_deviation": abs(
            recorded_final_dt - run_result.plan.effective_dt
        )
        / run_result.plan.effective_dt,
        # ---- the compared final-time quantities -----------------------------
        "final_amplification": float(final["amplification"]),
        "final_max_cartesian_vorticity": float(final["max_cartesian_vorticity"]),
        "final_max_abs_u1": float(final["max_abs_u1"]),
        "final_max_abs_omega1": float(final["max_abs_omega1"]),
        "final_max_abs_psi1": float(final["max_abs_psi1"]),
        "final_energy": float(final["energy"]),
        "final_enstrophy": float(final["enstrophy"]),
        "final_argmax_u1_r": float(final["argmax_u1_r"]),
        "final_argmax_u1_z": float(final["argmax_u1_z"]),
        "final_argmax_cartesian_vorticity_r": float(
            final["argmax_cartesian_vorticity_r"]
        ),
        "final_argmax_cartesian_vorticity_z": float(
            final["argmax_cartesian_vorticity_z"]
        ),
        "initial_energy": energy[0],
        "maximum_energy_growth_ratio": max(energy) / energy[0] - 1.0,
        "initial_circulation_max": circulation[0],
        "final_circulation_max": circulation[-1],
        "maximum_circulation_growth_ratio": max(circulation) / circulation[0] - 1.0,
        # ---- relativized constraint diagnostics (TM-09), recorded not gated -
        "final_divergence_residual_max": float(final["divergence_residual_max"]),
        "final_divergence_relative_denominator": float(
            final["divergence_relative_denominator"]
        ),
        "final_divergence_residual_relative": float(
            final["divergence_residual_relative"]
        ),
        "final_divergence_residual_argmax_r": float(
            final["divergence_residual_argmax_r"]
        ),
        "final_divergence_residual_argmax_z": float(
            final["divergence_residual_argmax_z"]
        ),
        "final_divergence_pointwise_ratio_max": float(
            final["divergence_pointwise_ratio_max"]
        ),
        "final_axis_parity_defect_u1": float(final["axis_parity_defect_u1"]),
        "final_axis_parity_defect_omega1": float(final["axis_parity_defect_omega1"]),
        "final_axis_parity_relative_u1": float(final["axis_parity_relative_u1"]),
        "final_axis_parity_relative_omega1": float(
            final["axis_parity_relative_omega1"]
        ),
        "final_axis_parity_relative_denominator_u1": float(
            final["axis_parity_relative_denominator_u1"]
        ),
        "final_axis_parity_relative_denominator_omega1": float(
            final["axis_parity_relative_denominator_omega1"]
        ),
        "final_odd_symmetry_defect": float(final["odd_symmetry_defect"]),
        "maximum_wall_u1_abs": max(_series(history, "wall_u1_max_abs")),
        # ---- independent solver-B elliptic cross-check at t_final -----------
        "final_cross_solver_psi_max_abs_difference": cross.get("max_abs_difference"),
        "final_cross_solver_psi_relative_difference": cross.get(
            "relative_difference"
        ),
        "final_cross_solver_psi_relative_denominator": cross.get(
            "relative_denominator"
        ),
        "final_cross_solver_psi_argmax_r": cross.get("argmax_r"),
        "final_cross_solver_psi_argmax_z": cross.get("argmax_z"),
        "final_solver_b_boundary_error_max": cross.get("boundary_error_max"),
        "final_solver_b_physical_defect_max": cross.get("physical_defect_max"),
        # ---- CFL monitoring (fixed dt, no adaptive selection) ---------------
        "maximum_advective_cfl": maximum_advective,
        "maximum_viscous_cfl": maximum_viscous,
        "maximum_advective_cfl_excess_ratio": maximum_advective / cfl_coefficient
        - 1.0,
        "maximum_viscous_cfl_excess_ratio": maximum_viscous / cfl_coefficient - 1.0,
    }


# --------------------------------------------------------- temporal analysis


def step_doubling_order(
    coarse: float,
    middle: float,
    fine: float,
    *,
    floor_relative: float,
) -> dict[str, Any]:
    """Return the step-doubling differences and the observed order of one scalar.

    ``coarse``/``middle``/``fine`` are ``Q(dt)``, ``Q(dt/2)`` and ``Q(dt/4)``.
    The order is ``log2(|Q(dt)-Q(dt/2)| / |Q(dt/2)-Q(dt/4)|)``.  When either
    difference falls at or below ``floor_relative * |Q(dt/4)|`` the order is
    not reported as a number; the status is then
    ``temporal_error_below_measurement_floor``, which is an accepted outcome
    with an explanation rather than a failure.
    """

    scale = abs(float(fine))
    floor = float(floor_relative) * scale
    values = (float(coarse), float(middle), float(fine))
    if not all(math.isfinite(value) for value in values):
        return {
            "coarse_value": values[0],
            "middle_value": values[1],
            "fine_value": values[2],
            "difference_coarse_middle": None,
            "difference_middle_fine": None,
            "measurement_floor": floor,
            "observed_order": None,
            "status": "non_finite",
            "differences_decrease": False,
        }
    first = abs(values[0] - values[1])
    second = abs(values[1] - values[2])
    common = {
        "coarse_value": values[0],
        "middle_value": values[1],
        "fine_value": values[2],
        "difference_coarse_middle": first,
        "difference_middle_fine": second,
        "measurement_floor": floor,
        "differences_decrease": bool(second < first),
    }
    if first <= floor or second <= floor:
        return {
            **common,
            "observed_order": None,
            "status": "temporal_error_below_measurement_floor",
        }
    return {
        **common,
        "observed_order": math.log2(first / second),
        "status": "measured",
    }


ANALYSABLE_STATUSES: frozenset[str] = frozenset(
    {"measured", "temporal_error_below_measurement_floor"}
)


def members_share_final_time(metrics: list[dict[str, Any]]) -> bool:
    """Return whether every member actually reached the same final time.

    Step doubling compares the members *at one instant*; if a member stopped
    early (an exhausted step budget, say) its value is not comparable and no
    order may be reported.
    """

    times = [item.get("final_time") for item in metrics]
    if any(
        value is None
        or isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        for value in times
    ):
        return False
    numbers = [float(value) for value in times]
    scale = max(abs(value) for value in numbers)
    if scale == 0.0:
        return True
    return (max(numbers) - min(numbers)) <= 1.0e-12 * scale


def temporal_order_analysis(
    metrics: list[dict[str, Any]],
    *,
    floor_relative: float,
) -> dict[str, dict[str, Any]]:
    """Run :func:`step_doubling_order` on every analysed smooth quantity."""

    if len(metrics) != EXPECTED_TIME_STEP_COUNT:
        raise ValueError("temporal order analysis needs exactly three members")
    comparable = members_share_final_time(metrics)
    analysis: dict[str, dict[str, Any]] = {}
    for name, source in SMOOTH_QUANTITIES:
        if not comparable:
            analysis[name] = {
                "diagnostic_key": source,
                "status": "members_disagree_on_final_time",
                "observed_order": None,
                "differences_decrease": False,
            }
            continue
        if any(name not in item for item in metrics):
            analysis[name] = {
                "diagnostic_key": source,
                "status": "unavailable",
                "observed_order": None,
                "differences_decrease": False,
            }
            continue
        analysis[name] = {
            "diagnostic_key": source,
            **step_doubling_order(
                metrics[0][name],
                metrics[1][name],
                metrics[2][name],
                floor_relative=floor_relative,
            ),
        }
    return analysis


def _circular_gap(first: float, second: float, period: float) -> float:
    raw = abs(float(first) - float(second)) % period
    return min(raw, period - raw)


def argmax_agreement(
    metrics: list[dict[str, Any]],
    *,
    z_period: float,
) -> dict[str, dict[str, Any]]:
    """Report the three argmax locations and whether they agree within a cell.

    The locations live on the single shared grid, so exact equality is the
    likely outcome.  This reports it; it never forces it.
    """

    report: dict[str, dict[str, Any]] = {}
    for name, radial_key, axial_key in ARGMAX_QUANTITIES:
        if any(radial_key not in item or axial_key not in item for item in metrics):
            report[name] = {"status": "unavailable"}
            continue
        radial = [float(item[radial_key]) for item in metrics]
        axial = [float(item[axial_key]) for item in metrics]
        dr = float(metrics[0]["dr"])
        dz = float(metrics[0]["dz"])
        radial_gap = max(abs(a - b) for a in radial for b in radial)
        axial_gap = max(
            _circular_gap(a, b, z_period) for a in axial for b in axial
        )
        report[name] = {
            "status": "recorded",
            "r": radial,
            "z": axial,
            "dr": dr,
            "dz": dz,
            "maximum_radial_separation": radial_gap,
            "maximum_axial_separation": axial_gap,
            "identical_across_the_ladder": bool(
                radial_gap == 0.0 and axial_gap == 0.0
            ),
            "agree_within_one_grid_cell": bool(
                radial_gap <= dr * (1.0 + 1.0e-12)
                and axial_gap <= dz * (1.0 + 1.0e-12)
            ),
        }
    return report


def constraint_trend(metrics: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Report the three values of every constraint metric and their spread.

    These quantities are expected to be dominated by the *spatial*
    discretization and therefore to be essentially independent of ``dt``.  The
    spread across the ladder is recorded so that expectation can be checked;
    nothing here is gated.
    """

    report: dict[str, dict[str, Any]] = {}
    for name in CONSTRAINT_QUANTITIES:
        values = [item.get(name) for item in metrics]
        if any(
            value is None
            or isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            for value in values
        ):
            report[name] = {"status": "unavailable", "values": values}
            continue
        numbers = [float(value) for value in values]
        largest = max(abs(value) for value in numbers)
        spread = max(numbers) - min(numbers)
        report[name] = {
            "status": "recorded",
            "values": numbers,
            "absolute_spread": spread,
            "relative_spread": (spread / largest) if largest > 0.0 else 0.0,
            "monotone_in_dt": bool(
                all(a <= b for a, b in zip(numbers, numbers[1:]))
                or all(a >= b for a, b in zip(numbers, numbers[1:]))
            ),
        }
    return report


def read_spatial_amplification_ladder(
    repository_root: Path = REPOSITORY_ROOT,
) -> dict[str, Any]:
    """Read the shipped spatial refinement ladder from the two summaries.

    ``outputs/hou_early_time_v1/summary.json`` and
    ``outputs/hou_early_time_v2_hires/summary.json`` record the identical
    E-29/E-30 early-time run at the identical ``t_final = T_1`` on four
    resolutions.  Their ``final_amplification`` values are the spatial half of
    the temporal/spatial separation statement.
    """

    found: dict[tuple[int, int], dict[str, Any]] = {}
    read_files: list[str] = []
    missing: list[str] = []
    for relative in SPATIAL_LADDER_SUMMARIES:
        path = repository_root / relative
        if not path.is_file():
            missing.append(relative)
            continue
        payload = strict_json_loads(
            path.read_text(encoding="utf-8"),
            label=f"hou early time summary {relative}",
        )
        read_files.append(relative)
        if not isinstance(payload, dict):
            continue
        for item in payload.get("resolutions", []):
            key = (int(item["nr"]), int(item["nz"]))
            found.setdefault(
                key,
                {
                    "nr": key[0],
                    "nz": key[1],
                    "final_amplification": float(item["final_amplification"]),
                    "summary": relative,
                },
            )
    rungs = [dict(item) for item in SPATIAL_AMPLIFICATION_LADDER]
    matches = True
    for rung in rungs:
        observed = found.get((int(rung["nr"]), int(rung["nz"])))
        if observed is None:
            rung["read_from_summary"] = False
            matches = False
            continue
        rung["read_from_summary"] = True
        rung["summary"] = observed["summary"]
        rung["final_amplification_read"] = observed["final_amplification"]
        if observed["final_amplification"] != rung["final_amplification"]:
            matches = False
    amplifications = [float(item["final_amplification"]) for item in rungs]
    differences = [
        abs(fine - coarse)
        for coarse, fine in zip(amplifications, amplifications[1:])
    ]
    return {
        "source_summaries": list(SPATIAL_LADDER_SUMMARIES),
        "summaries_read": read_files,
        "summaries_missing": missing,
        "matches_recorded_reference": bool(matches),
        "rungs": rungs,
        "final_amplifications": amplifications,
        "consecutive_differences": differences,
        "smallest_consecutive_difference": min(differences) if differences else None,
        "largest_consecutive_difference": max(differences) if differences else None,
    }


def temporal_vs_spatial_separation(
    metrics: list[dict[str, Any]],
    analysis: dict[str, dict[str, Any]],
    constraints: dict[str, dict[str, Any]],
    *,
    repository_root: Path = REPOSITORY_ROOT,
) -> dict[str, Any]:
    """Contrast the measured temporal error with the recorded spatial ladder."""

    spatial = read_spatial_amplification_ladder(repository_root)
    amplification = analysis.get("final_amplification", {})
    temporal_difference = amplification.get("difference_middle_fine")
    richardson = (
        None
        if temporal_difference is None
        else float(temporal_difference) / (2.0**NOMINAL_TEMPORAL_ORDER - 1.0)
    )
    smallest_spatial = spatial["smallest_consecutive_difference"]
    ratio: float | None = None
    dominant = "undetermined"
    if richardson is not None and smallest_spatial is not None:
        if richardson > 0.0:
            ratio = float(smallest_spatial) / richardson
            dominant = "spatial" if ratio > 1.0 else "temporal"
        elif float(smallest_spatial) > 0.0:
            # The temporal difference vanished at the recorded precision, so
            # the ratio is not a finite number; the spatial term is the only
            # one left.
            dominant = "spatial"
    resolution = (
        f"{int(metrics[0]['nr'])}x{int(metrics[0]['nz'])}" if metrics else "unknown"
    )
    return {
        "statement": (
            "The three members share one spatial grid and one final time, so "
            "their mutual differences contain no spatial error at all: they "
            "are a pure measurement of the temporal discretization. The "
            "spatial error is measured separately by the resolution ladder of "
            "the early-time runs at the identical t_final."
        ),
        "fixed_resolution": resolution,
        "temporal": {
            "quantity": "final_amplification",
            "step_doubling_difference_dt_over_2_to_dt_over_4": temporal_difference,
            "richardson_error_estimate_at_finest_dt": richardson,
            "richardson_definition": (
                "|Q(dt/2) - Q(dt/4)| / (2^p - 1) with p = 2 (Heun); the "
                "standard estimate of the temporal error remaining in Q(dt/4)"
            ),
            "status": amplification.get("status"),
            "measurement_floor": amplification.get("measurement_floor"),
        },
        "spatial": spatial,
        "spatial_to_temporal_ratio": ratio,
        "spatial_to_temporal_ratio_definition": (
            "smallest consecutive spatial-ladder amplification difference "
            "divided by the Richardson temporal error estimate at the finest dt"
        ),
        "dominant_error_source": dominant,
        "constraint_metrics_expected_spatially_dominated": {
            "expectation": (
                "The relative divergence residual, the relative axis parity "
                "defect, the circulation max-principle violation and the "
                "solver-B cross-check difference are all set by the spatial "
                "stencils, so they should be essentially independent of dt. "
                "Their spread across the ladder is recorded here; a small "
                "spread confirms the separation and is itself the "
                "temporal/spatial separation statement for those metrics."
            ),
            "gated": False,
            "metrics": constraints,
        },
    }


# ------------------------------------------------------------------ evaluate


def cfl_excess_tolerance_of(config: dict[str, Any]) -> float:
    """Return the configured CFL excess tolerance or its documented default."""

    if "cfl_excess_tolerance" not in config:
        return DEFAULT_CFL_EXCESS_TOLERANCE
    return _number(config["cfl_excess_tolerance"], name="cfl_excess_tolerance")


def cfl_policy(
    config: dict[str, Any],
    metrics: list[dict[str, Any]],
) -> dict[str, Any]:
    """Return the summary block that fixes the CFL monitoring semantics."""

    coefficient = float(config["cfl_coefficient"])
    tolerance = cfl_excess_tolerance_of(config)
    advective = [
        float(item["maximum_advective_cfl"])
        for item in metrics
        if item.get("maximum_advective_cfl") is not None
    ]
    viscous = [
        float(item["maximum_viscous_cfl"])
        for item in metrics
        if item.get("maximum_viscous_cfl") is not None
    ]
    worst_advective = max(advective) if advective else None
    worst_viscous = max(viscous) if viscous else None
    return {
        "mode": "fixed uniform time step; no adaptive selection is performed",
        "monitoring_only": True,
        "rule": (
            "the fixed dt is accepted when the recorded effective CFL numbers "
            "stay at or below cfl_coefficient * (1 + cfl_excess_tolerance)"
        ),
        "equivalence": (
            "the adaptive advective bound is C / max(|u^r|/dr, |u^z|/dz) and "
            "the recorded advective_cfl is dt * max(|u^r|/dr, |u^z|/dz), so "
            "dt <= bound is exactly advective_cfl <= C; the viscous bound "
            "C min(dr,dz)^2 / (4 nu) and viscous_cfl = 4 nu dt / min(dr,dz)^2 "
            "are related the same way"
        ),
        "effective_cfl_definition": (
            "advective_cfl = dt * max(max|u^r| / dr, max|u^z| / dz) and "
            "viscous_cfl = 4 nu dt / min(dr, dz)^2, both evaluated on the "
            "state AFTER the step"
        ),
        "cfl_coefficient": coefficient,
        "cfl_excess_tolerance": tolerance,
        "cfl_excess_tolerance_source": (
            "config key cfl_excess_tolerance"
            if "cfl_excess_tolerance" in config
            else f"default {DEFAULT_CFL_EXCESS_TOLERANCE}"
        ),
        "accepted_effective_cfl_bound": coefficient * (1.0 + tolerance),
        "maximum_effective_advective_cfl": worst_advective,
        "maximum_effective_viscous_cfl": worst_viscous,
        "maximum_effective_advective_cfl_excess_ratio": (
            None if worst_advective is None else worst_advective / coefficient - 1.0
        ),
        "maximum_effective_viscous_cfl_excess_ratio": (
            None if worst_viscous is None else worst_viscous / coefficient - 1.0
        ),
        "exceedance_policy": (
            "an exceedance is recorded as a failed acceptance check "
            "(fixed_step_within_adaptive_cfl_bound), never raised"
        ),
    }


def evaluate(
    config: dict[str, Any],
    *,
    checkpoint_dir: Path | None = None,
    provenance: dict[str, Any] | None = None,
    repository_root: Path = REPOSITORY_ROOT,
) -> tuple[
    list[TimeStepRun],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, bool],
]:
    """Run the three members and reduce them to metrics, analysis and checks."""

    validate_config(config)
    acceptance = config["acceptance"]
    coefficient = float(config["cfl_coefficient"])
    tolerance = cfl_excess_tolerance_of(config)
    floor_relative = float(acceptance["temporal_error_floor_relative"])
    plans = plan_ladder(float(config["t_final"]), list(config["time_steps"]))
    halves = ladder_halves_exactly(plans)

    runs = [
        evolve_time_step(
            config,
            index=index,
            plan=plan,
            checkpoint_dir=checkpoint_dir,
            provenance=provenance,
        )
        for index, plan in enumerate(plans)
    ]
    metrics = [
        level_metrics(item, cfl_coefficient=coefficient) for item in runs
    ]
    usable = [item for item in metrics if "final_amplification" in item]

    analysis = temporal_order_analysis(metrics, floor_relative=floor_relative)
    argmax = argmax_agreement(metrics, z_period=float(config["z_period"]))
    constraints = constraint_trend(metrics)
    separation = temporal_vs_spatial_separation(
        metrics,
        analysis,
        constraints,
        repository_root=repository_root,
    )

    measured = [
        item for item in analysis.values() if item.get("status") == "measured"
    ]
    lower = float(acceptance["temporal_order_lower_bound"])
    upper = float(acceptance["temporal_order_upper_bound"])
    in_band = [
        item
        for item in measured
        if lower <= float(item["observed_order"]) <= upper
    ]
    required = min(
        int(acceptance["minimum_quantities_at_second_order"]), len(measured)
    )
    below_floor = [
        item
        for item in analysis.values()
        if item.get("status") == "temporal_error_below_measurement_floor"
    ]

    checks = {
        "all_time_steps_completed": bool(usable)
        and len(usable) == EXPECTED_TIME_STEP_COUNT
        and all(item["completed"] for item in metrics),
        "all_diagnostics_finite": bool(usable)
        and all(item["all_diagnostics_finite"] for item in usable),
        "effective_time_steps_halve_exactly": bool(halves),
        "energy_non_increasing": bool(usable)
        and all(
            item["maximum_energy_growth_ratio"]
            <= float(acceptance["maximum_energy_growth_ratio"])
            for item in usable
        ),
        "circulation_max_principle": bool(usable)
        and all(
            item["maximum_circulation_growth_ratio"]
            <= float(acceptance["maximum_circulation_growth_ratio"])
            for item in usable
        ),
        # Heun is second order.  Quantities whose step-doubling differences sit
        # at or below the measurement floor are an accepted outcome with an
        # explanation and are excluded from the count instead of failing.
        "temporal_order_second_order_or_below_floor": bool(usable)
        and len(analysis) == len(SMOOTH_QUANTITIES)
        and all(
            item.get("status") in ANALYSABLE_STATUSES for item in analysis.values()
        )
        and len(in_band) >= required,
        "step_doubling_differences_decrease": bool(usable)
        and all(item["differences_decrease"] for item in measured),
        "fixed_step_within_adaptive_cfl_bound": bool(usable)
        and all(
            advective_cfl_within_tolerance(
                item["maximum_advective_cfl"],
                cfl_coefficient=coefficient,
                cfl_excess_tolerance=tolerance,
            )
            and cfl_within_tolerance(
                item["maximum_viscous_cfl"],
                cfl_coefficient=coefficient,
                cfl_excess_tolerance=tolerance,
            )
            for item in usable
        ),
        "cross_solver_elliptic_agreement_recorded": bool(usable)
        and all(
            item.get("final_cross_solver_psi_max_abs_difference") is not None
            and math.isfinite(
                float(item["final_cross_solver_psi_max_abs_difference"])
            )
            for item in usable
        ),
    }
    report = {
        "temporal_order": analysis,
        "temporal_order_definition": (
            "log2(|Q(dt) - Q(dt/2)| / |Q(dt/2) - Q(dt/4)|) with all three "
            "members on the identical spatial grid and identical t_final"
        ),
        "temporal_order_band": [lower, upper],
        "members_share_final_time": members_share_final_time(metrics),
        "measured_quantity_count": len(measured),
        "quantities_within_the_order_band": [
            name
            for name, item in analysis.items()
            if item.get("status") == "measured"
            and lower <= float(item["observed_order"]) <= upper
        ],
        "quantities_below_measurement_floor": [
            name
            for name, item in analysis.items()
            if item.get("status") == "temporal_error_below_measurement_floor"
        ],
        "quantities_required_within_the_band": required,
        "all_quantities_below_measurement_floor": bool(
            len(below_floor) == len(SMOOTH_QUANTITIES)
        ),
        "measurement_floor_definition": (
            "acceptance.temporal_error_floor_relative * |Q(dt/4)|; a quantity "
            "whose step-doubling differences fall at or below it is recorded "
            "as temporal_error_below_measurement_floor and is an accepted "
            "outcome with an explanation, not a failure"
        ),
        "argmax_locations": argmax,
        "constraint_metric_trends": constraints,
        "step_plans": [item.as_dict() for item in plans],
        "effective_time_steps_halve_exactly": bool(halves),
        "temporal_vs_spatial_separation": separation,
    }
    return runs, metrics, report, checks


# -------------------------------------------------------------------- output


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


LEVEL_CSV_FIELDS: tuple[str, ...] = (
    "level_index",
    "requested_dt",
    "effective_dt",
    "step_count",
    "observed_step_count",
    "requested_dt_divides_t_final",
    "effective_dt_relative_adjustment",
    "final_recorded_dt",
    "final_step_relative_deviation",
    "nr",
    "nz",
    "dr",
    "dz",
    "completed",
    "final_time",
    "final_amplification",
    "final_max_cartesian_vorticity",
    "final_max_abs_u1",
    "final_max_abs_omega1",
    "final_energy",
    "final_argmax_u1_r",
    "final_argmax_u1_z",
    "final_argmax_cartesian_vorticity_r",
    "final_argmax_cartesian_vorticity_z",
    "maximum_circulation_growth_ratio",
    "maximum_energy_growth_ratio",
    "final_divergence_residual_relative",
    "final_divergence_pointwise_ratio_max",
    "final_axis_parity_relative_u1",
    "final_axis_parity_relative_omega1",
    "final_cross_solver_psi_max_abs_difference",
    "final_cross_solver_psi_relative_difference",
    "final_cross_solver_psi_relative_denominator",
    "final_solver_b_boundary_error_max",
    "final_solver_b_physical_defect_max",
    "maximum_advective_cfl",
    "maximum_viscous_cfl",
    "checkpoint",
)


def run(
    config: dict[str, Any],
    output_dir: Path,
    *,
    repository_root: Path = REPOSITORY_ROOT,
) -> dict[str, Any]:
    """Run the three members and write a checksummed evidence bundle."""

    validate_config(config)
    provenance = collect_runtime_provenance()
    config_bytes = (
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    _prepare_output(output_dir)
    (output_dir / "config.snapshot.json").write_bytes(config_bytes)
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir()
    runs, metrics, report, checks = evaluate(
        config,
        checkpoint_dir=checkpoint_dir,
        provenance=provenance,
        repository_root=repository_root,
    )

    diagnostic_rows = [
        {
            "level_index": item.index,
            "requested_dt": item.plan.requested_dt,
            "effective_dt": item.plan.effective_dt,
            **record,
        }
        for item in runs
        for record in item.history
    ]
    if not diagnostic_rows:
        raise RuntimeError("no diagnostics were recorded")
    _write_csv(
        output_dir / "diagnostics.csv",
        ["level_index", "requested_dt", "effective_dt", *DIAGNOSTIC_FIELDS],
        diagnostic_rows,
    )
    _write_csv(
        output_dir / "time_refinement.csv",
        list(LEVEL_CSV_FIELDS),
        metrics,
    )

    arrays: dict[str, np.ndarray] = {
        "r": runs[0].grid.r,
        "z": runs[0].grid.z,
        "requested_dt": np.asarray(
            [item.plan.requested_dt for item in runs], dtype=np.float64
        ),
        "effective_dt": np.asarray(
            [item.plan.effective_dt for item in runs], dtype=np.float64
        ),
    }
    for item in runs:
        suffix = f"dt{item.index:02d}"
        if item.state is not None:
            arrays[f"u1_final_{suffix}"] = item.state.u1
            arrays[f"omega1_final_{suffix}"] = item.state.omega1
            arrays[f"psi1_final_{suffix}"] = item.state.psi1
        for name in (
            "time",
            "energy",
            "max_abs_u1",
            "max_abs_omega1",
            "max_cartesian_vorticity",
            "amplification",
        ):
            arrays[f"{name}_{suffix}"] = np.asarray(
                _series(item.history, name) if item.history else [],
                dtype=np.float64,
            )
    np.savez_compressed(output_dir / "trajectories.npz", **arrays)

    nr, nz = (int(value) for value in config["resolution"])
    grid = build_grid(config, nr=nr, nz=nz)
    norms = initial_norms(grid, amplitude_scale=float(config["amplitude_scale"]))

    summary: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Fixed-grid, fixed-step temporal refinement of the E-29 initial "
            "datum at the audited early time T_1. It measures the temporal "
            "discretization error of the shipped Heun integrator and nothing "
            "else: it is not a reproduction of any published amplification "
            "value, not a singularity candidate, and not a proof."
        ),
        "method": {
            "equations": "E-11, E-12, E-13, E-14",
            "time_integrator": "explicit Heun/RK2 (nominal order 2)",
            "time_step_rule": (
                "fixed uniform dt per member; the ladder is dt, dt/2, dt/4 "
                "with step counts n, 2n, 4n landing exactly on t_final"
            ),
            "fixed_spatial_grid": {
                "nr": nr,
                "nz": nz,
                "dr": grid.dr,
                "dz": grid.dz,
                "r_max": float(config["r_max"]),
                "z_period": float(config["z_period"]),
            },
            "wall_conditions": "E-27 with the second-order Thom form E-31",
            "poisson_solver": (
                "solver A: Fourier-z / r^3-flux finite volume with homogeneous "
                "outer Dirichlet data"
            ),
            "independent_cross_check": (
                "solver B: finite_cylinder_poisson non-divergence radial "
                "stencil, re-solved on the final state of every member"
            ),
            "viscosity_protocol": (
                "single E-30 stage: t_final = T_1 lies below the switch time "
                "t_0 = 0.00227375, so no viscosity switch occurs"
            ),
            "initial_data": "E-29 with amplitude 12000 times amplitude_scale",
            "symmetry": "full z period, odd symmetry monitored not imposed",
            "constraint_relativization": (
                "the E-02 divergence residual and the E-16c axis parity defect "
                "are reported both absolutely and divided by their own "
                "cancellation scales (TM-09); the relative numbers are "
                "recorded and no acceptance check reads them"
            ),
            "relativized_diagnostic_fields": list(RELATIVE_DIAGNOSTIC_FIELDS),
        },
        "cfl_policy": cfl_policy(config, metrics),
        "reproducibility": {
            "seed": int(config["seed"]),
            "config_sha256": _sha256(config_bytes),
            "runtime_provenance": provenance,
            "command": (
                "python experiments/run_hou_time_refinement.py --config "
                "configs/hou_time_refinement.json --output-dir "
                "outputs/hou_time_refinement"
            ),
        },
        "derived_reference_norms": {
            "source": "E-29b (derived by this repository, absent from the paper)",
            "max_abs_u1": E29B_MAX_ABS_U1,
            "max_cartesian_vorticity": E29B_MAX_CARTESIAN_VORTICITY,
            "measured_max_abs_u1": norms["max_abs_u1"],
            "measured_max_cartesian_vorticity": norms["max_cartesian_vorticity"],
            "max_abs_u1_relative_error": norms["max_abs_u1_relative_error"],
            "max_cartesian_vorticity_relative_error": norms[
                "max_cartesian_vorticity_relative_error"
            ],
        },
        "time_steps": [dict(item) for item in metrics],
        # The separation block is promoted to the top level of the summary and
        # is therefore removed from the nested analysis rather than duplicated.
        "temporal_analysis": {
            name: value
            for name, value in report.items()
            if name != "temporal_vs_spatial_separation"
        },
        "temporal_vs_spatial_separation": report["temporal_vs_spatial_separation"],
        "acceptance_checks": checks,
        "accepted_as_temporal_refinement_study": bool(all(checks.values())),
        "known_gaps": list(KNOWN_GAPS),
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
        "time_refinement.csv",
        "trajectories.npz",
    ]
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
        default=REPOSITORY_ROOT / "configs" / "hou_time_refinement.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "hou_time_refinement",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="hou time refinement config",
    )
    if not isinstance(config, dict):
        parser.error("config must be a JSON object")
    summary = run(config, output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["accepted_as_temporal_refinement_study"] else 2


if __name__ == "__main__":
    sys.exit(main())
