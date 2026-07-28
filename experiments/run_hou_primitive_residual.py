#!/usr/bin/env python
r"""Pressure-free primitive Navier--Stokes residual on real Hou checkpoints.

``outputs/hou_snapshot_cartesian_audit_v1/summary.json`` records an explicit
gap (``primitive_navier_stokes_residual_gap``): the independent Cartesian audit
cannot evaluate the primitive momentum residual on the shipped Hou checkpoints
because (i) a checkpoint stores no pressure and (ii) the shipped snapshot
spacing ``5e-4`` is about five hundred integrator steps, so a two-point time
difference over it is not a controlled ``u_t``.

This experiment closes that gap without inventing a pressure.

Primary result: the pressure-free residuals
----------------------------------------------------------------------------

For the exact incompressible equations the momentum residual

``R(u) = u_t + (u.grad)u - nu*Laplacian(u)``

equals ``-grad p``, so ``curl R = 0`` pointwise.  Two different discrete forms
of that statement are evaluated:

* ``curl R``, the curl of the assembled momentum residual, and
* ``S = omega_t + (u.grad)omega - (omega.grad)u - nu*Laplacian(omega)`` with
  ``omega = curl u``, the Cartesian vorticity-transport residual.

They agree in the continuum for a divergence-free field but are genuinely
different discretizations of the nonlinear term (TM-14), so both are reported.
All stencils come from :mod:`ns_certificate_lab.cartesian_momentum_audit`,
which uses the public surface of :mod:`ns_certificate_lab.cartesian_validation`
and imports no cylindrical operator, PDE or Poisson module.

The controlled time derivative
----------------------------------------------------------------------------

Both quantities need ``u_t``.  This experiment therefore *produces* the states
it audits: it integrates the E-29 datum with the shipped fixed-step mode and
saves a symmetric fan of schema-v2 candidates at the offsets
:data:`STATE_OFFSETS` around one target time, each offset being an exact whole
number of integrator steps.  The centered differences at spacings ``dt``,
``2 dt`` and ``4 dt`` then give three independent approximations of the same
``u_t``, whose mutual differences measure the temporal error and whose ratio
measures its order (nominally two).  Nothing here relies on the ``5e-4``
snapshot spacing of the shipped run.

Non-circularity
----------------------------------------------------------------------------

The evaluation never reads the in-memory integrator state.  Every state is
written to disk with :func:`nonlinear_cylinder.save_checkpoint`, reloaded with
``artifacts.load_candidate`` (which verifies the recorded array digests) and
mapped to the Cartesian box by the existing adapter.  The digests of the
reloaded archives are recorded in the summary and in the manifest.

Secondary result: a scoped pressure recovery
----------------------------------------------------------------------------

``R = -grad p`` is inverted in the discrete least-squares sense on the same
Cartesian box (see :func:`cartesian_momentum_audit.recover_pressure`).  The
box is periodic in ``z`` but its ``x`` and ``y`` faces are an artificial cut
through the source cylinder, where the true pressure is not determined by data
inside the box.  The recovered field is therefore a **scoped diagnostic and
never a validated pressure**; only the interior consistency
``max|R + grad p| / (term sum)`` is reported and it is **record-only**.  The
whole path is verified separately on an analytic field with a closed-form
pressure, where the recovered pressure gradient converges at second order.

Limits
----------------------------------------------------------------------------

Nothing here is an interval enclosure, a discretization error bound, a
reproduction claim, or evidence for or against singularity formation.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np

from experiments.run_hou_early_time import (
    DEFAULT_CFL_EXCESS_TOLERANCE,
    E29_AMPLITUDE,
    advective_cfl_within_tolerance,
    build_grid,
)
from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    write_with_digest,
)
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.cartesian_candidate_adapter import (
    reconstruct_loaded_candidate_on_cartesian,
)
from ns_certificate_lab.cartesian_momentum_audit import (
    DENOMINATORS,
    curl_momentum_residual,
    defect_report,
    divergence_report,
    finite_values,
    observed_order,
    ratio,
    recover_pressure,
    region_offset,
    region_window,
    rms,
    velocity_scales,
    vorticity_transport_residual,
)
from ns_certificate_lab.cartesian_validation import (
    UniformCartesianGrid,
    cartesian_curl,
)
from ns_certificate_lab.nonlinear_cylinder import (
    ConstrainedState,
    hou_initial_swirl,
    integrate,
    normalize_viscosity_schedule,
    save_checkpoint,
    viscosity_at,
)
from ns_certificate_lab.provenance import collect_runtime_provenance

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPECTED_EXPERIMENT_ID = "hou_primitive_residual_v1"
EXPECTED_COORDINATE_SYSTEM = (
    "uniform Cartesian (x,y,z) box [-a,a]^2 x [z_min,z_max) inscribed in the "
    "finite cylinder of the saved checkpoints"
)
SUMMARY_SCHEMA_VERSION = 1
MANIFEST_SCHEMA = "ns-certificate-lab/hou-primitive-residual-manifest/v1"

#: Saved states, in whole integrator steps relative to the target time.  The
#: symmetric offsets 1, 2 and 4 give centered time differences at spacings
#: ``dt``, ``2 dt`` and ``4 dt``; three spacings are the minimum from which an
#: *order* can be measured rather than merely an error estimate.
STATE_OFFSETS: tuple[int, ...] = (-4, -2, -1, 0, 1, 2, 4)

#: Half-spacings, in integrator steps, of the three centered differences.
TIME_DIFFERENCE_HALF_WIDTHS: tuple[int, ...] = (1, 2, 4)

#: ``integrate`` re-derives its own step count as ``ceil(t_final / dt)``.
#: Handing it exactly ``t_final / N`` can round up to ``N + 1``; nudging the
#: requested step up by one part in 1e9 makes that ceiling reproduce ``N``.
#: Identical in purpose to ``run_hou_time_refinement.STEP_COUNT_NUDGE``.
STEP_COUNT_NUDGE = 1.0 + 1.0e-9

#: Relative tolerance for "the target time is an integer number of steps".
DIVISIBILITY_RELATIVE_TOLERANCE = 1.0e-9

#: Heun/RK2 produces the states; the centered time difference itself is the
#: quantity whose order is measured here, and it is nominally second order.
NOMINAL_TIME_DIFFERENCE_ORDER = 2.0

LIMITATIONS: tuple[str, ...] = (
    "The Cartesian reconstruction samples the cylindrical nodal arrays with "
    "BILINEAR interpolation, and the residuals evaluated here take up to three "
    "nested derivatives of that interpolant.  Every number below therefore "
    "contains an interpolation error that this experiment does not separate "
    "from a genuine inconsistency of the states.  The audit-grid refinement "
    "block is the only handle on that separation and it is recorded, not "
    "assumed.",
    "The uniform Cartesian audit box is generally coarser than the source "
    "cylindrical grid, so the reported defects are dominated by audit "
    "truncation rather than by the solver.  A defect that shrinks under audit "
    "refinement is audit truncation.",
    "The x and y faces of the audit box are an artificial cut through the "
    "source cylinder.  The interior_margin norms exist because the one-sided "
    "closures there are not a physical boundary condition of the source "
    "problem.",
    "The recovered pressure is a scoped diagnostic and NOT a validated "
    "pressure: the x,y boundary treatment is a least-squares artefact, not the "
    "boundary behaviour of the true pressure of the source problem.  Its "
    "consistency metric is record-only.",
    "The time difference is centered on the saved states, so it measures the "
    "time derivative of the DISCRETE trajectory, not of any continuum "
    "solution.  A small pressure-free residual therefore certifies that the "
    "discrete trajectory is consistent with the Cartesian momentum equation to "
    "the audit's own truncation level; it certifies nothing about the "
    "continuum limit.",
    "Floating-point binary64 arithmetic with tolerances only.  There is no "
    "interval enclosure, no discretization error bound, and nothing here is "
    "evidence for or against singularity formation.",
)

KNOWN_GAPS: tuple[str, ...] = (
    "The two pressure-free forms are numerically indistinguishable on a "
    "decaying two-dimensional Taylor-Green field because that field's "
    "advection term is an exact gradient, so BOTH nonlinear contributions "
    "vanish identically.  The manufactured cross-check that actually "
    "exercises the nonlinear discretizations therefore uses a separate "
    "nontrivial solenoidal field (see tests/test_hou_primitive_residual.py).",
    "This experiment audits the states IT produces.  It does not evaluate the "
    "pressure-free residual on the shipped outputs/hou_early_time_v1 "
    "checkpoints, whose 5e-4 spacing still admits no controlled u_t.  The "
    "shipped snapshots remain un-audited in this respect.",
    "The recovered pressure is never compared with an independently computed "
    "pressure of the source problem, because none exists in this repository.",
    "At the E-30 viscosity nu=5e-4 and this audit resolution the viscous term "
    "of both pressure-free residuals is roughly 1e-4 of the time-derivative "
    "and advection terms (see term_magnitudes in every region block).  A "
    "sign error in the viscous term is therefore NOT detectable on these "
    "states: it moves the residual by far less than the audit truncation "
    "already does.  That fault injection is exercised on a manufactured field "
    "with a viscosity large enough to make the term visible, and the negative "
    "result on the real states is recorded rather than hidden.",
)

_SQRT_TWO = math.sqrt(2.0)


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------


def _nonempty_text(value: Any, *, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _finite_positive(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite positive number")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a finite positive number")
    return result


def _finite_number(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _positive_integer(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _shape(value: Any, *, name: str) -> tuple[int, int, int]:
    if (
        not isinstance(value, list)
        or len(value) != 3
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
    ):
        raise ValueError(f"{name} must be a list of three integers")
    if any(item < 4 for item in value):
        raise ValueError(f"{name} entries must each be at least 4")
    return (int(value[0]), int(value[1]), int(value[2]))


def validate_config(config: Mapping[str, Any]) -> None:
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
        "time_step",
        "target_time",
        "cfl_coefficient",
        "max_steps",
        "diagnostic_stride",
        "cartesian_grid",
        "pressure_recovery",
        "hard_gates",
        "record_only_thresholds",
    }
    optional = {"cfl_excess_tolerance"}
    if (
        not isinstance(config, dict)
        or not required <= set(config)
        or not set(config) <= required | optional
    ):
        raise ValueError("hou primitive residual config has missing or unknown keys")
    if config["schema_version"] != SUMMARY_SCHEMA_VERSION:
        raise ValueError("schema_version must be 1")
    if config["experiment_id"] != EXPECTED_EXPERIMENT_ID:
        raise ValueError(f"experiment_id must be {EXPECTED_EXPERIMENT_ID}")
    for key in ("description", "interpretation", "normalization_description"):
        _nonempty_text(config[key], name=key)
    if config["coordinate_system"] != EXPECTED_COORDINATE_SYSTEM:
        raise ValueError("coordinate_system must equal the audited canonical value")
    if not isinstance(config["units"], dict) or not config["units"]:
        raise ValueError("units must be a nonempty object")
    if isinstance(config["seed"], bool) or not isinstance(config["seed"], int):
        raise ValueError("seed must be an integer")

    _finite_positive(config["r_max"], name="r_max")
    _finite_positive(config["z_period"], name="z_period")
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
    _finite_positive(config["amplitude_scale"], name="amplitude_scale")

    schedule = normalize_viscosity_schedule(config["viscosity_schedule"])
    time_step = _finite_positive(config["time_step"], name="time_step")
    target_time = _finite_positive(config["target_time"], name="target_time")
    plan = plan_states(target_time=target_time, time_step=time_step)
    if len(schedule) > 1 and schedule[1][0] < plan.final_time:
        raise ValueError(
            "the primitive residual study requires a single viscosity stage "
            "over [0, target_time + max_offset*dt]: no switch time may fall "
            "below the final integration time"
        )

    _finite_positive(config["cfl_coefficient"], name="cfl_coefficient")
    if "cfl_excess_tolerance" in config:
        tolerance = _finite_number(
            config["cfl_excess_tolerance"], name="cfl_excess_tolerance"
        )
        if tolerance < 0.0:
            raise ValueError("cfl_excess_tolerance must be nonnegative")
    max_steps = _positive_integer(config["max_steps"], name="max_steps")
    if max_steps < plan.step_count:
        raise ValueError(
            "max_steps is below the step count the target time requires"
        )
    _positive_integer(config["diagnostic_stride"], name="diagnostic_stride")

    grid_spec = config["cartesian_grid"]
    grid_keys = {
        "half_width",
        "z_min",
        "z_max",
        "z_periodic",
        "shape",
        "coarse_shape",
        "interior_margin",
    }
    if not isinstance(grid_spec, dict) or set(grid_spec) != grid_keys:
        raise ValueError("cartesian_grid has missing or unknown keys")
    half_width = _finite_positive(grid_spec["half_width"], name="half_width")
    if half_width * _SQRT_TWO > float(config["r_max"]) * (1.0 + 1.0e-12):
        raise ValueError(
            "Cartesian audit box corners leave the source cylinder: "
            "half_width*sqrt(2) must not exceed r_max"
        )
    z_min = _finite_number(grid_spec["z_min"], name="z_min")
    z_max = _finite_number(grid_spec["z_max"], name="z_max")
    if z_max <= z_min:
        raise ValueError("z_max must exceed z_min")
    if not isinstance(grid_spec["z_periodic"], bool):
        raise ValueError("z_periodic must be a boolean")
    fine = _shape(grid_spec["shape"], name="shape")
    coarse = _shape(grid_spec["coarse_shape"], name="coarse_shape")
    if any(c >= f for c, f in zip(coarse, fine)):
        raise ValueError(
            "coarse_shape must be strictly coarser than shape on every axis"
        )
    margin = grid_spec["interior_margin"]
    if isinstance(margin, bool) or not isinstance(margin, int) or margin < 0:
        raise ValueError("interior_margin must be a nonnegative integer")
    if any(size <= 2 * margin for size in (fine[0], fine[1], coarse[0], coarse[1])):
        raise ValueError("interior_margin leaves an empty norm region")

    recovery = config["pressure_recovery"]
    recovery_keys = {"enabled", "tolerance", "max_iterations"}
    if not isinstance(recovery, dict) or set(recovery) != recovery_keys:
        raise ValueError("pressure_recovery has missing or unknown keys")
    if not isinstance(recovery["enabled"], bool):
        raise ValueError("pressure_recovery.enabled must be a boolean")
    _finite_positive(recovery["tolerance"], name="pressure_recovery.tolerance")
    _positive_integer(
        recovery["max_iterations"], name="pressure_recovery.max_iterations"
    )

    gates = config["hard_gates"]
    expected_gates = {
        "maximum_relative_curl_momentum_residual_rms",
        "maximum_relative_vorticity_transport_residual_rms",
        "time_difference_order_lower_bound",
        "time_difference_order_upper_bound",
    }
    if not isinstance(gates, dict) or set(gates) != expected_gates:
        raise ValueError("hard_gates has missing or unknown keys")
    for key in expected_gates:
        _finite_positive(gates[key], name=f"hard_gates.{key}")
    if not (
        float(gates["time_difference_order_lower_bound"])
        < float(gates["time_difference_order_upper_bound"])
    ):
        raise ValueError(
            "hard_gates.time_difference_order_lower_bound must be below "
            "hard_gates.time_difference_order_upper_bound"
        )

    thresholds = config["record_only_thresholds"]
    expected_thresholds = {
        "reference_relative_curl_momentum_residual_max",
        "reference_relative_vorticity_transport_residual_max",
        "reference_relative_pressure_consistency_max",
        "reference_relative_divergence_rms",
    }
    if not isinstance(thresholds, dict) or set(thresholds) != expected_thresholds:
        raise ValueError("record_only_thresholds has missing or unknown keys")
    for key in expected_thresholds:
        _finite_positive(thresholds[key], name=f"record_only_thresholds.{key}")


# --------------------------------------------------------------------------
# step planning
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class StatePlan:
    """The exact step grid on which the saved fan of states lives."""

    requested_time_step: float
    target_time: float
    steps_to_target: int
    step_count: int
    effective_time_step: float
    final_time: float
    effective_target_time: float
    requested_step_divides_target: bool
    effective_step_relative_adjustment: float

    def state_time(self, offset: int) -> float:
        """Physical time of the state ``offset`` whole steps from the target."""

        return float(self.steps_to_target + int(offset)) * self.effective_time_step

    def snapshot_times(self) -> tuple[float, ...]:
        return tuple(self.state_time(offset) for offset in STATE_OFFSETS)

    def as_dict(self) -> dict[str, Any]:
        return {
            "requested_time_step": self.requested_time_step,
            "target_time": self.target_time,
            "steps_to_target": self.steps_to_target,
            "step_count": self.step_count,
            "effective_time_step": self.effective_time_step,
            "final_time": self.final_time,
            "effective_target_time": self.effective_target_time,
            "requested_step_divides_target": self.requested_step_divides_target,
            "effective_step_relative_adjustment": (
                self.effective_step_relative_adjustment
            ),
            "state_offsets": list(STATE_OFFSETS),
            "state_times": list(self.snapshot_times()),
        }


def plan_states(*, target_time: float, time_step: float) -> StatePlan:
    """Return the step grid that places every saved state on a step boundary.

    ``target_time`` must be an integer number of ``time_step`` values (to
    within :data:`DIVISIBILITY_RELATIVE_TOLERANCE`) and must be far enough from
    zero that the earliest saved state, ``max(STATE_OFFSETS)`` steps before it,
    still has a strictly positive time.  Anything else is rejected: a fan that
    does not sit on exact step boundaries would silently turn the centered
    differences into uncontrolled approximations.
    """

    stop = _finite_positive(target_time, name="target_time")
    step = _finite_positive(time_step, name="time_step")
    raw = stop / step
    nearest = float(round(raw))
    divides = nearest >= 1.0 and abs(raw - nearest) <= (
        DIVISIBILITY_RELATIVE_TOLERANCE * max(1.0, raw)
    )
    if not divides:
        raise ValueError(
            "target_time must be an integer multiple of time_step so that "
            "every saved state lands on an exact integrator step boundary"
        )
    steps_to_target = int(nearest)
    reach = max(STATE_OFFSETS)
    if steps_to_target <= reach:
        raise ValueError(
            "target_time is not reachable: it must be more than "
            f"{reach} time steps from the initial time so that the earliest "
            "saved state has a positive time"
        )
    step_count = steps_to_target + reach
    final_time = float(step_count) * step
    effective = final_time / float(step_count)
    return StatePlan(
        requested_time_step=step,
        target_time=stop,
        steps_to_target=steps_to_target,
        step_count=step_count,
        effective_time_step=effective,
        final_time=final_time,
        effective_target_time=float(steps_to_target) * effective,
        requested_step_divides_target=True,
        effective_step_relative_adjustment=abs(effective - step) / step,
    )


# --------------------------------------------------------------------------
# production of the audited states
# --------------------------------------------------------------------------


def state_filename(offset: int) -> str:
    """Deterministic checkpoint name for one saved offset."""

    sign = "m" if int(offset) < 0 else "p"
    return f"state_{sign}{abs(int(offset)):02d}dt.npz"


@dataclass(frozen=True)
class ProductionResult:
    """Outcome of the integration that produced the saved fan of states."""

    plan: StatePlan
    completed: bool
    failure: str | None
    observed_step_count: int
    final_time: float
    history: tuple[dict[str, float], ...]
    saved: dict[int, str]
    grid_nr: int
    grid_nz: int
    grid_dr: float
    grid_dz: float
    viscosity: float


def produce_states(
    config: Mapping[str, Any],
    checkpoint_dir: Path,
    *,
    provenance: Mapping[str, Any] | None = None,
) -> ProductionResult:
    """Integrate the E-29 datum and save the fan of schema-v2 candidates."""

    validate_config(config)
    plan = plan_states(
        target_time=float(config["target_time"]),
        time_step=float(config["time_step"]),
    )
    nr, nz = (int(value) for value in config["resolution"])
    grid = build_grid(config, nr=nr, nz=nz)
    schedule = normalize_viscosity_schedule(config["viscosity_schedule"])
    amplitude_scale = float(config["amplitude_scale"])
    u1_initial = hou_initial_swirl(
        grid,
        amplitude=E29_AMPLITUDE,
        amplitude_scale=amplitude_scale,
    )
    wanted = {
        plan.state_time(offset): int(offset) for offset in STATE_OFFSETS
    }
    saved: dict[int, str] = {}
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def on_snapshot(
        moment: float,
        state: ConstrainedState,
        _diagnostics: dict[str, float],
    ) -> None:
        offset = wanted.get(moment)
        if offset is None:
            # Fall back to the nearest requested time; ``integrate`` compares
            # snapshot times with an absolute tolerance, so an exact dict hit
            # is expected but not guaranteed.
            best = min(wanted, key=lambda value: abs(value - moment))
            if abs(best - moment) > 1.0e-12:
                return
            offset = wanted[best]
        name = state_filename(offset)
        save_checkpoint(
            checkpoint_dir / name,
            grid=grid,
            state=state,
            time=float(moment),
            viscosity=viscosity_at(schedule, float(moment)),
            seed=int(config["seed"]),
            config=config,
            provenance=provenance,
            metadata={
                "experiment_id": config["experiment_id"],
                "nr": nr,
                "nz": nz,
                "amplitude_scale": amplitude_scale,
                "step_offset_from_target": int(offset),
                "effective_time_step": plan.effective_time_step,
                "steps_to_target": plan.steps_to_target,
            },
        )
        saved[int(offset)] = name

    failure: str | None = None
    try:
        result = integrate(
            grid,
            u1=u1_initial,
            omega1=np.zeros(grid.shape),
            t_final=plan.final_time,
            viscosity_schedule=schedule,
            fixed_time_step=plan.effective_time_step * STEP_COUNT_NUDGE,
            max_steps=int(config["max_steps"]),
            diagnostic_stride=int(config["diagnostic_stride"]),
            snapshot_times=plan.snapshot_times(),
            on_snapshot=on_snapshot,
        )
    except (FloatingPointError, ArithmeticError, ValueError) as exc:
        return ProductionResult(
            plan=plan,
            completed=False,
            failure=f"{type(exc).__name__}: {exc}",
            observed_step_count=0,
            final_time=0.0,
            history=(),
            saved=dict(saved),
            grid_nr=nr,
            grid_nz=nz,
            grid_dr=float(grid.dr),
            grid_dz=float(grid.dz),
            viscosity=viscosity_at(schedule, 0.0),
        )

    if not result.completed:
        failure = (
            f"integration stopped after {result.step_count} steps at "
            f"t={result.time!r} before the planned final time"
        )
    elif result.step_count != plan.step_count:
        failure = (
            f"integrator took {result.step_count} steps where the validated "
            f"plan requires {plan.step_count}"
        )
    elif set(saved) != set(STATE_OFFSETS):
        failure = (
            "the integrator did not deliver every requested state: missing "
            f"offsets {sorted(set(STATE_OFFSETS) - set(saved))}"
        )
    return ProductionResult(
        plan=plan,
        completed=failure is None,
        failure=failure,
        observed_step_count=int(result.step_count),
        final_time=float(result.time),
        history=result.history,
        saved=dict(saved),
        grid_nr=nr,
        grid_nz=nz,
        grid_dr=float(grid.dr),
        grid_dz=float(grid.dz),
        viscosity=viscosity_at(schedule, float(plan.effective_target_time)),
    )


# --------------------------------------------------------------------------
# reloading through the saved artifacts
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ReloadedState:
    """One saved candidate, reloaded and reconstructed on a Cartesian box."""

    offset: int
    path: Path
    archive_sha256: str
    physical_time: float
    viscosity: float
    velocity: np.ndarray
    vorticity_e18b: np.ndarray


def build_cartesian_grid(
    grid_spec: Mapping[str, Any],
    shape: tuple[int, int, int],
) -> UniformCartesianGrid:
    """Build the audit box ``[-a,a]^2 x [z_min, z_max)`` at ``shape``."""

    half_width = float(grid_spec["half_width"])
    return UniformCartesianGrid.uniform(
        shape=shape,
        x_bounds=(-half_width, half_width),
        y_bounds=(-half_width, half_width),
        z_bounds=(float(grid_spec["z_min"]), float(grid_spec["z_max"])),
        periodic=(False, False, bool(grid_spec["z_periodic"])),
    )


def reload_states(
    checkpoint_dir: Path,
    grid: UniformCartesianGrid,
    *,
    offsets: Sequence[int] = STATE_OFFSETS,
) -> dict[int, ReloadedState]:
    """Reload every saved state from disk and reconstruct it on ``grid``.

    The evaluation path deliberately starts here: nothing that follows can see
    the integrator's in-memory arrays, so the reported residuals are computed
    from the same bytes an external reader would get (TM-14).
    """

    states: dict[int, ReloadedState] = {}
    for offset in offsets:
        path = checkpoint_dir / state_filename(int(offset))
        candidate = load_candidate(path)
        if candidate.description is None:
            raise ValueError(f"saved state lacks a v2 description: {path}")
        arrays = reconstruct_loaded_candidate_on_cartesian(candidate, grid)
        states[int(offset)] = ReloadedState(
            offset=int(offset),
            path=path,
            archive_sha256=sha256_file(path),
            physical_time=float(candidate.description.physical_time),
            viscosity=float(candidate.description.viscosity),
            velocity=arrays.velocity,
            vorticity_e18b=arrays.vorticity_e18b,
        )
    return states


# --------------------------------------------------------------------------
# time differences
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class CenteredDifference:
    """One centered time difference and the interval it was taken over."""

    half_width: int
    spacing: float
    derivative: np.ndarray


def centered_differences(
    states: Mapping[int, ReloadedState],
    *,
    field: str = "velocity",
) -> dict[int, CenteredDifference]:
    """Return the centered time differences at every half width.

    The spacing is taken from the *recorded* physical times of the reloaded
    checkpoints, never from the nominal step, so a mismatch between the plan
    and what the integrator actually did cannot hide here.
    """

    result: dict[int, CenteredDifference] = {}
    for half_width in TIME_DIFFERENCE_HALF_WIDTHS:
        early = states[-half_width]
        late = states[half_width]
        spacing = late.physical_time - early.physical_time
        if not math.isfinite(spacing) or spacing <= 0.0:
            raise ValueError(
                "reloaded states are not ordered in time: the centered "
                "difference has a non-positive spacing"
            )
        result[half_width] = CenteredDifference(
            half_width=int(half_width),
            spacing=float(spacing),
            derivative=(getattr(late, field) - getattr(early, field)) / spacing,
        )
    return result


def richardson_time_difference_report(
    differences: Mapping[int, CenteredDifference],
    *,
    window: tuple[slice, slice, slice],
) -> dict[str, Any]:
    r"""Measure the order of the centered time difference itself.

    ``D_h = u_t + C h^2 + O(h^4)`` gives ``D_{2h} - D_h = 3 C h^2`` and
    ``D_{4h} - D_{2h} = 12 C h^2``, so the ratio of those two differences is
    ``4`` and its base-two logarithm is the observed order.  The remaining
    error of the tightest difference is estimated by ``|D_{2h} - D_h| / 3``.
    """

    selector = (slice(None), *window)
    one, two, four = (differences[width] for width in TIME_DIFFERENCE_HALF_WIDTHS)
    first = one.derivative[selector]
    second = two.derivative[selector]
    third = four.derivative[selector]
    scale = rms(first)
    gap_low = rms(second - first)
    gap_high = rms(third - second)
    order: float | None = None
    if gap_low > 0.0 and gap_high > 0.0:
        value = math.log2(gap_high / gap_low)
        order = float(value) if math.isfinite(value) else None
    estimate = gap_low / 3.0
    return {
        "definition": (
            "D_h = (u(t+h) - u(t-h)) / 2h for h = dt, 2dt, 4dt; the observed "
            "order is log2(||D_4h - D_2h|| / ||D_2h - D_h||), which is 2 when "
            "the centered difference is second order"
        ),
        "spacings": {
            str(width): differences[width].spacing
            for width in TIME_DIFFERENCE_HALF_WIDTHS
        },
        "derivative_rms": {
            str(width): rms(differences[width].derivative[selector])
            for width in TIME_DIFFERENCE_HALF_WIDTHS
        },
        "difference_rms_2h_minus_h": gap_low,
        "difference_rms_4h_minus_2h": gap_high,
        "differences_decrease": bool(gap_low < gap_high),
        "observed_order": order,
        "richardson_error_estimate_at_dt": estimate,
        "richardson_error_estimate_relative": ratio(estimate, scale),
        "derivative_scale_rms": scale,
    }


# --------------------------------------------------------------------------
# the pressure-free evaluation
# --------------------------------------------------------------------------


def _term_fractions(
    terms: Mapping[str, np.ndarray],
    window: tuple[Any, ...],
) -> dict[str, Any]:
    """RMS of every term over ``window`` plus its share of the total."""

    magnitudes = {name: rms(field[window]) for name, field in terms.items()}
    total = sum(magnitudes.values())
    return {
        "rms": magnitudes,
        "fraction_of_term_rms_sum": {
            name: (ratio(value, total) if total > 0.0 else None)
            for name, value in magnitudes.items()
        },
        "term_rms_sum": float(total),
    }


def evaluate_on_grid(
    grid: UniformCartesianGrid,
    states: Mapping[int, ReloadedState],
    *,
    viscosity: float,
    interior_margin: int,
    pressure_recovery: Mapping[str, Any],
    viscous_sign: float = 1.0,
    velocity_t_override: np.ndarray | None = None,
) -> dict[str, Any]:
    """Evaluate both pressure-free residuals on one Cartesian audit grid.

    ``viscous_sign`` and ``velocity_t_override`` exist only for the fault
    injections: the physical values are ``1.0`` and ``None``, and any other
    value deliberately breaks the evaluated equation.
    """

    centre = states[0]
    window = region_window(grid, int(interior_margin))
    offset = region_offset(grid, int(interior_margin))
    full_window = region_window(grid, 0)
    full_offset = region_offset(grid, 0)
    mesh = grid.mesh()

    velocity_differences = centered_differences(states, field="velocity")
    if velocity_t_override is None:
        velocity_t = velocity_differences[1].derivative
        velocity_t_source = "centered difference at half width 1 (spacing dt)"
    else:
        velocity_t = grid.validate_vector(
            velocity_t_override, name="velocity_t_override"
        )
        velocity_t_source = "caller supplied override (fault injection)"

    vorticity = {
        key: cartesian_curl(grid, state.velocity) for key, state in states.items()
    }
    if velocity_t_override is None:
        omega_early = vorticity[-1]
        omega_late = vorticity[1]
        spacing = velocity_differences[1].spacing
        vorticity_t = (omega_late - omega_early) / spacing
    else:
        vorticity_t = cartesian_curl(grid, velocity_t)

    curl_form = curl_momentum_residual(
        grid,
        velocity=centre.velocity,
        velocity_t=velocity_t,
        viscosity=float(viscosity),
        viscous_sign=float(viscous_sign),
    )
    transport_form = vorticity_transport_residual(
        grid,
        velocity=centre.velocity,
        vorticity=vorticity[0],
        vorticity_t=vorticity_t,
        viscosity=float(viscosity),
        viscous_sign=float(viscous_sign),
    )

    regions: dict[str, Any] = {}
    for label, active_window, active_offset, margin in (
        ("full", full_window, full_offset, 0),
        ("interior", window, offset, int(interior_margin)),
    ):
        vector_window = (slice(None), *active_window)
        block: dict[str, Any] = {
            "margin": int(margin),
            "term_magnitudes": {
                "definition": (
                    "RMS of each term over this region, and that RMS divided "
                    "by the sum of the three (respectively four) term RMS "
                    "values.  A term whose fraction is far below the reported "
                    "relative residual cannot be probed by this measurement at "
                    "all: corrupting it would move the residual by less than "
                    "the audit truncation already does."
                ),
                "curl_momentum": _term_fractions(
                    {
                        "curl_time_derivative": curl_form.curl_time_derivative,
                        "curl_advection": curl_form.curl_advection,
                        "curl_viscous": curl_form.curl_viscous,
                    },
                    vector_window,
                ),
                "vorticity_transport": _term_fractions(
                    {
                        "time_derivative": transport_form.time_derivative,
                        "advection": transport_form.advection,
                        "stretching": transport_form.stretching,
                        "viscous": transport_form.viscous,
                    },
                    vector_window,
                ),
            },
            "scales": velocity_scales(
                grid,
                centre.velocity,
                vorticity[0],
                window=active_window,
            ),
            "curl_momentum_residual": defect_report(
                grid,
                defect=curl_form.residual,
                term_sum=curl_form.term_sum,
                window=active_window,
                offset=active_offset,
                mesh=mesh,
            ),
            "vorticity_transport_residual": defect_report(
                grid,
                defect=transport_form.residual,
                term_sum=transport_form.term_sum,
                window=active_window,
                offset=active_offset,
                mesh=mesh,
            ),
            "form_difference": defect_report(
                grid,
                defect=curl_form.residual - transport_form.residual,
                term_sum=np.abs(curl_form.residual)
                + np.abs(transport_form.residual),
                window=active_window,
                offset=active_offset,
                mesh=mesh,
            ),
            "divergence": divergence_report(
                grid,
                centre.velocity,
                window=active_window,
                offset=active_offset,
                mesh=mesh,
            ),
        }
        regions[label] = block

    richardson = richardson_time_difference_report(
        velocity_differences, window=window
    )

    # Record-only: the same residual assembled from the 2dt-spaced difference.
    wide = velocity_differences[2]
    wide_curl = curl_momentum_residual(
        grid,
        velocity=centre.velocity,
        velocity_t=wide.derivative,
        viscosity=float(viscosity),
        viscous_sign=float(viscous_sign),
    )
    selector = (slice(None), *window)
    spacing_sensitivity = {
        "definition": (
            "the curl-of-momentum residual re-assembled with the 2dt-spaced "
            "centered difference instead of the dt-spaced one; a residual "
            "dominated by the spatial audit truncation barely moves"
        ),
        "relative_rms_at_dt": regions["interior"]["curl_momentum_residual"][
            "relative_rms"
        ],
        "relative_rms_at_2dt": ratio(
            rms(wide_curl.residual[selector]),
            float(np.max(np.abs(wide_curl.term_sum[selector]))),
        ),
        "absolute_rms_at_dt": regions["interior"]["curl_momentum_residual"]["rms"],
        "absolute_rms_at_2dt": rms(wide_curl.residual[selector]),
    }

    recovery: dict[str, Any] = {
        "enabled": bool(pressure_recovery["enabled"]),
        "status": "disabled",
    }
    if bool(pressure_recovery["enabled"]):
        recovered = recover_pressure(
            grid,
            curl_form.momentum.residual,
            momentum_term_sum=curl_form.momentum.term_sum,
            tolerance=float(pressure_recovery["tolerance"]),
            max_iterations=int(pressure_recovery["max_iterations"]),
        )
        recovery = {
            "enabled": True,
            "status": "recorded",
            "scope": (
                "SCOPED DIAGNOSTIC, NOT A VALIDATED PRESSURE.  The x and y "
                "faces of the audit box are an artificial cut through the "
                "source cylinder, so the least-squares problem supplies its "
                "own boundary behaviour there rather than the boundary "
                "behaviour of the true pressure.  Only the interior "
                "consistency is reported and it is record-only."
            ),
            "method": (
                "minimize ||grad_h p + R||^2 over the whole box; the normal "
                "equations are div grad p = -div R written with the audit's "
                "own gradient and its exact algebraic transpose"
            ),
            "iterations": recovered.iterations,
            "converged": recovered.converged,
            "normal_equation_relative_residual": (
                recovered.normal_equation_relative_residual
            ),
            "null_space_dimension": recovered.null_space_dimension,
            "pressure_max_abs": float(np.max(np.abs(recovered.pressure))),
            "consistency": {
                label: defect_report(
                    grid,
                    defect=recovered.consistency,
                    term_sum=recovered.term_sum,
                    window=active_window,
                    offset=active_offset,
                    mesh=mesh,
                )
                for label, active_window, active_offset in (
                    ("full", full_window, full_offset),
                    ("interior", window, offset),
                )
            },
        }

    return {
        "shape": list(grid.shape),
        "spacings": list(grid.spacings),
        "periodic": [bool(value) for value in grid.periodic],
        "interior_margin": int(interior_margin),
        "viscosity": float(viscosity),
        "velocity_t_source": velocity_t_source,
        "centre_physical_time": centre.physical_time,
        "regions": regions,
        "time_difference_richardson": richardson,
        "time_difference_spacing_sensitivity": spacing_sensitivity,
        "pressure_recovery": recovery,
    }


# --------------------------------------------------------------------------
# whole evaluation
# --------------------------------------------------------------------------


def evaluate(
    config: Mapping[str, Any],
    checkpoint_dir: Path,
    *,
    viscosity: float,
) -> dict[str, Any]:
    """Evaluate the saved fan on the fine and the coarse audit grid."""

    validate_config(config)
    grid_spec = config["cartesian_grid"]
    margin = int(grid_spec["interior_margin"])
    evaluations: dict[str, Any] = {}
    reloaded_digests: dict[str, dict[str, Any]] = {}
    for label, shape_key in (("fine", "shape"), ("coarse", "coarse_shape")):
        shape = _shape(grid_spec[shape_key], name=shape_key)
        grid = build_cartesian_grid(grid_spec, shape)
        states = reload_states(checkpoint_dir, grid)
        if not reloaded_digests:
            reloaded_digests = {
                str(offset): {
                    "file": state.path.name,
                    "archive_sha256": state.archive_sha256,
                    "physical_time": state.physical_time,
                    "viscosity": state.viscosity,
                }
                for offset, state in sorted(states.items())
            }
        evaluations[label] = evaluate_on_grid(
            grid,
            states,
            viscosity=float(viscosity),
            interior_margin=margin,
            pressure_recovery=config["pressure_recovery"],
        )

    fine = evaluations["fine"]["regions"]["interior"]
    coarse = evaluations["coarse"]["regions"]["interior"]
    spacing_ratio = max(evaluations["coarse"]["spacings"]) / max(
        evaluations["fine"]["spacings"]
    )
    refinement: dict[str, Any] = {"spacing_ratio": float(spacing_ratio)}
    for name in (
        "curl_momentum_residual",
        "vorticity_transport_residual",
        "form_difference",
        "divergence",
    ):
        coarse_value = coarse[name]["relative_rms"]
        fine_value = fine[name]["relative_rms"]
        refinement[name] = {
            "relative_rms_coarse": coarse_value,
            "relative_rms_fine": fine_value,
            "observed_order": observed_order(
                coarse_value, fine_value, float(spacing_ratio)
            ),
            "decreases_under_refinement": bool(
                coarse_value is not None
                and fine_value is not None
                and fine_value < coarse_value
            ),
        }
    return {
        "evaluations": evaluations,
        "audit_grid_refinement": refinement,
        "reloaded_states": reloaded_digests,
    }


def cfl_excess_tolerance_of(config: Mapping[str, Any]) -> float:
    """Return the configured CFL excess tolerance or its documented default."""

    if "cfl_excess_tolerance" not in config:
        return DEFAULT_CFL_EXCESS_TOLERANCE
    return _finite_number(
        config["cfl_excess_tolerance"], name="cfl_excess_tolerance"
    )


def acceptance(
    config: Mapping[str, Any],
    production: ProductionResult,
    analysis: Mapping[str, Any],
) -> tuple[dict[str, bool], dict[str, bool]]:
    """Assemble the hard acceptance checks and the record-only checks."""

    gates = config["hard_gates"]
    thresholds = config["record_only_thresholds"]
    fine = analysis["evaluations"]["fine"]
    interior = fine["regions"]["interior"]
    richardson = fine["time_difference_richardson"]
    order = richardson["observed_order"]
    coefficient = float(config["cfl_coefficient"])
    tolerance = cfl_excess_tolerance_of(config)
    advective = [float(record["advective_cfl"]) for record in production.history]
    viscous = [float(record["viscous_cfl"]) for record in production.history]

    def within(name: str, gate_key: str) -> bool:
        value = interior[name]["relative_rms"]
        return value is not None and float(value) <= float(gates[gate_key])

    checks = {
        "integration_completed": bool(production.completed),
        "every_requested_state_was_saved": set(production.saved)
        == set(STATE_OFFSETS),
        "all_states_finite": all(
            math.isfinite(value)
            for record in production.history
            for value in record.values()
        ),
        "all_metrics_finite": all(
            math.isfinite(value) for value in finite_values(analysis)
        ),
        "time_difference_order_within_band": (
            order is not None
            and float(gates["time_difference_order_lower_bound"])
            <= float(order)
            <= float(gates["time_difference_order_upper_bound"])
        ),
        "time_difference_errors_decrease_with_spacing": bool(
            richardson["differences_decrease"]
        ),
        "relative_curl_momentum_residual_rms_within_gate": within(
            "curl_momentum_residual",
            "maximum_relative_curl_momentum_residual_rms",
        ),
        "relative_vorticity_transport_residual_rms_within_gate": within(
            "vorticity_transport_residual",
            "maximum_relative_vorticity_transport_residual_rms",
        ),
        "residuals_decrease_under_audit_grid_refinement": bool(
            analysis["audit_grid_refinement"]["curl_momentum_residual"][
                "decreases_under_refinement"
            ]
            and analysis["audit_grid_refinement"]["vorticity_transport_residual"][
                "decreases_under_refinement"
            ]
        ),
        "fixed_step_within_adaptive_cfl_bound": bool(advective)
        and advective_cfl_within_tolerance(
            max(advective),
            cfl_coefficient=coefficient,
            cfl_excess_tolerance=tolerance,
        )
        and max(viscous) <= coefficient * (1.0 + tolerance),
    }

    def below(name: str, threshold_key: str, key: str = "relative_max") -> bool:
        value = interior[name][key]
        return value is not None and float(value) <= float(thresholds[threshold_key])

    consistency = fine["pressure_recovery"].get("consistency")
    record_only = {
        "relative_curl_momentum_residual_max_below_reference": below(
            "curl_momentum_residual",
            "reference_relative_curl_momentum_residual_max",
        ),
        "relative_vorticity_transport_residual_max_below_reference": below(
            "vorticity_transport_residual",
            "reference_relative_vorticity_transport_residual_max",
        ),
        "relative_divergence_rms_below_reference": below(
            "divergence",
            "reference_relative_divergence_rms",
            key="relative_rms",
        ),
        "pressure_recovery_converged": bool(
            fine["pressure_recovery"].get("converged", False)
        ),
        "relative_pressure_consistency_max_below_reference": bool(
            consistency is not None
            and consistency["interior"]["relative_max"] is not None
            and float(consistency["interior"]["relative_max"])
            <= float(thresholds["reference_relative_pressure_consistency_max"])
        ),
    }
    return checks, record_only


# --------------------------------------------------------------------------
# evidence bundle
# --------------------------------------------------------------------------


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
        path.resolve().relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def metric_rows(analysis: Mapping[str, Any]) -> list[dict[str, Any]]:
    """One CSV row per (audit grid, region, quantity)."""

    rows: list[dict[str, Any]] = []
    for label, evaluation in sorted(analysis["evaluations"].items()):
        for region, block in sorted(evaluation["regions"].items()):
            for quantity in (
                "curl_momentum_residual",
                "vorticity_transport_residual",
                "form_difference",
                "divergence",
            ):
                report = block[quantity]
                rows.append(
                    {
                        "audit_grid": label,
                        "audit_shape": "x".join(
                            str(value) for value in evaluation["shape"]
                        ),
                        "region": region,
                        "margin": block["margin"],
                        "quantity": quantity,
                        "rms": report["rms"],
                        "max": report["max"],
                        "term_sum_max": report["term_sum_max"],
                        "relative_rms": (
                            "" if report["relative_rms"] is None
                            else report["relative_rms"]
                        ),
                        "relative_max": (
                            "" if report["relative_max"] is None
                            else report["relative_max"]
                        ),
                        "argmax_x": report["argmax"]["x"],
                        "argmax_y": report["argmax"]["y"],
                        "argmax_z": report["argmax"]["z"],
                        "argmax_radius": report["argmax"]["radius"],
                        "argmax_component": (
                            "" if report["argmax"]["component"] is None
                            else report["argmax"]["component"]
                        ),
                        **{
                            f"scale_{key}": value
                            for key, value in sorted(block["scales"].items())
                        },
                    }
                )
    return rows


def run(config: Mapping[str, Any], output_dir: Path) -> dict[str, Any]:
    """Produce the states, evaluate them from disk and write the bundle."""

    validate_config(config)
    if not _inside_repository(output_dir):
        raise ValueError("output directory must remain inside this repository")
    provenance = collect_runtime_provenance()
    config_bytes = canonical_json_bytes(dict(config))
    _prepare_output(output_dir)

    checkpoint_dir = output_dir / "checkpoints"
    production = produce_states(config, checkpoint_dir, provenance=provenance)
    if not production.saved:
        raise RuntimeError(
            "the integration produced no saved states: "
            f"{production.failure or 'unknown failure'}"
        )
    analysis = evaluate(config, checkpoint_dir, viscosity=production.viscosity)
    checks, record_only = acceptance(config, production, analysis)

    summary: dict[str, Any] = {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "experiment_id": config["experiment_id"],
        "interpretation": config["interpretation"],
        "closes_recorded_gap": {
            "gap": (
                "outputs/hou_snapshot_cartesian_audit_v1/summary.json key "
                "primitive_navier_stokes_residual_gap"
            ),
            "closed": (
                "the PRESSURE-FREE half: curl of the momentum residual and the "
                "Cartesian vorticity-transport residual, both evaluated on "
                "real saved states with a controlled centered u_t"
            ),
            "still_open": (
                "the pressure-carrying half: no independent pressure of the "
                "source problem exists in this repository, so the recovered "
                "pressure is a scoped least-squares diagnostic and is "
                "record-only"
            ),
        },
        "method": {
            "production": (
                "ns_certificate_lab.nonlinear_cylinder.integrate in its fixed "
                "uniform time-step mode from the E-29 datum; the audited "
                "states are written as schema-v2 candidates at exact whole "
                "step offsets around the target time"
            ),
            "reconstruction": (
                "ns_certificate_lab.cartesian_candidate_adapter: private "
                "second-order (r,z) stencils and bilinear sampling of the "
                "RELOADED nodal arrays"
            ),
            "checker": (
                "ns_certificate_lab.cartesian_momentum_audit on top of the "
                "public surface of ns_certificate_lab.cartesian_validation; "
                "neither imports operators.py, pde.py, poisson.py, "
                "finite_cylinder_poisson.py or nonlinear_cylinder.py"
            ),
            "non_circularity": (
                "every evaluated array is reloaded from disk with "
                "artifacts.load_candidate, which verifies the recorded array "
                "digests; the integrator's in-memory state is never read by "
                "the evaluation"
            ),
            "curl_momentum_residual": (
                "curl(u_t) + curl((u.grad)u) - nu*curl(Laplacian u), which "
                "vanishes for an exact solution because R = -grad p"
            ),
            "vorticity_transport_residual": (
                "omega_t + (u.grad)omega - (omega.grad)u - nu*Laplacian(omega) "
                "with omega = curl u, a different discretization of the same "
                "continuum statement"
            ),
            "time_derivative": (
                "centered difference of the reloaded states at half widths "
                "dt, 2dt and 4dt; the order of the difference itself is "
                "measured from the three of them"
            ),
            "norm_region": (
                "interior_margin trims the non-periodic x and y axes only; "
                "the periodic z seam is always retained"
            ),
        },
        "denominators": dict(DENOMINATORS),
        "state_plan": production.plan.as_dict(),
        "production": {
            "completed": production.completed,
            "failure": production.failure,
            "observed_step_count": production.observed_step_count,
            "final_time": production.final_time,
            "diagnostic_record_count": len(production.history),
            "resolution": [production.grid_nr, production.grid_nz],
            "dr": production.grid_dr,
            "dz": production.grid_dz,
            "viscosity_at_target": production.viscosity,
            "maximum_advective_cfl": (
                max(float(record["advective_cfl"]) for record in production.history)
                if production.history
                else None
            ),
            "maximum_viscous_cfl": (
                max(float(record["viscous_cfl"]) for record in production.history)
                if production.history
                else None
            ),
            "cfl_coefficient": float(config["cfl_coefficient"]),
            "cfl_excess_tolerance": cfl_excess_tolerance_of(config),
            "saved_states": {
                str(offset): name for offset, name in sorted(production.saved.items())
            },
        },
        "cartesian_grid": dict(config["cartesian_grid"]),
        "pressure_recovery_policy": {
            "record_only": True,
            "scope": (
                "SCOPED DIAGNOSTIC, NOT A VALIDATED PRESSURE.  The x,y faces "
                "of the audit box are an artificial cut through the source "
                "cylinder; the least-squares problem supplies its own "
                "boundary behaviour there.  Only the interior consistency "
                "max|R + grad p| / (term sum) is meaningful, and no "
                "acceptance gate reads it."
            ),
            "verified_on": (
                "an analytic decaying Taylor-Green field with a closed-form "
                "pressure, where the recovered pressure gradient converges at "
                "second order (tests/test_hou_primitive_residual.py)"
            ),
        },
        "hard_gates": dict(config["hard_gates"]),
        "record_only_thresholds": dict(config["record_only_thresholds"]),
        "reloaded_states": analysis["reloaded_states"],
        "evaluations": analysis["evaluations"],
        "audit_grid_refinement": analysis["audit_grid_refinement"],
        "acceptance_checks": checks,
        "record_only_checks": record_only,
        "accepted_as_pressure_free_primitive_residual_study": bool(
            all(checks.values())
        ),
        "known_gaps": list(KNOWN_GAPS),
        "limitations": list(LIMITATIONS),
        "reproducibility": {
            "seed": int(config["seed"]),
            "config_sha256": sha256_bytes(config_bytes),
            "runtime_provenance": provenance,
            "command": (
                "python -m experiments.run_hou_primitive_residual "
                "--config configs/hou_primitive_residual.json "
                "--output-dir outputs/hou_primitive_residual_v1"
            ),
        },
    }

    write_with_digest(output_dir / "config.snapshot.json", config_bytes)
    write_with_digest(output_dir / "summary.json", canonical_json_bytes(summary))
    rows = metric_rows(analysis)
    with (output_dir / "metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    artifact_names = [
        "config.snapshot.json",
        "config.snapshot.json.sha256",
        "metrics.csv",
        "summary.json",
        "summary.json.sha256",
    ]
    artifact_names.extend(
        sorted(
            path.relative_to(output_dir).as_posix()
            for path in checkpoint_dir.rglob("*")
            if path.is_file()
        )
    )
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "experiment_id": config["experiment_id"],
        "files": {
            name: {
                "sha256": sha256_file(output_dir / name),
                "bytes": (output_dir / name).stat().st_size,
            }
            for name in sorted(artifact_names)
        },
    }
    write_with_digest(output_dir / "manifest.json", canonical_json_bytes(manifest))
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "hou_primitive_residual.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "hou_primitive_residual_v1",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    try:
        config = strict_json_loads(
            config_path.read_text(encoding="utf-8"),
            label="hou primitive residual config",
        )
        if not isinstance(config, dict):
            raise ValueError("config must be a JSON object")
        summary = run(config, output_dir)
    except (OSError, ValueError, FileExistsError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary["acceptance_checks"], indent=2, sort_keys=True))
    return 0 if summary["accepted_as_pressure_free_primitive_residual_study"] else 2


if __name__ == "__main__":
    sys.exit(main())
