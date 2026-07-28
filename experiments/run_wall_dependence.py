#!/usr/bin/env python
r"""Preregistered wall-dependence experiment for the Hou finite-cylinder run.

``docs/wall_dependence_prereg.md`` asks one question: does the Hou growth
mechanism (E-29 datum, E-27 wall conditions) depend on the wall at ``r = 1``?
The design fixed there is a nested family of cylinders
``R_wall in {1, 1.5, 2, 3}`` sharing an *identical core discretization*
``dr = 1/192``, ``dz = 1/384`` (so ``nr = 192 R_wall + 1`` and ``nz = 384``),
plus a resolution-consistency pair ``R_wall in {1, 2}`` at the coarser
``dr = 1/128``.

Because the plain E-29 datum cannot be extended past ``r = 1``, every member
starts from the E-32 ``C^infinity`` compact-support family
(:mod:`ns_certificate_lab.wall_dependence`): E-29 multiplied by a smooth cutoff
of ``r^2`` that is exactly one for ``r <= 0.9`` and exactly zero for
``r >= 0.95``.  Multiplication by the literal ``1.0`` is exact, so *the core
initial datum is bit-identical across every wall radius*, and the datum
vanishes strictly inside every wall.

What is measured (preregistration section 3), on ``Omega_c = {r <= 0.9}`` at
every snapshot
--------------------------------------------------------------------------

1. the amplification ``A_R(t) = max_{Omega_c} |omega|(t) / ||omega(0)||_inf``
   against the *shared* reference norm,
2. the ``argmax`` location of ``u_1``,
3. ``max_{Omega_c} |psi_{1,z}|``, the meridional flow strength,
4. the core kinetic energy (E-20a restricted to ``Omega_c``),
5. the elliptic nonlocal contribution: one core-restricted ``omega_1`` solved on
   two different wall radii and compared nodewise on ``Omega_c``,
6. the tail amplitudes ``max_{r>1.1} |u_1|, |omega_1|`` for ``R_wall > 1``.

Decision rule (preregistration section 4) -- RECORD ONLY
--------------------------------------------------------

``S(R,R') = |A_R(T_1) - A_{R'}(T_1)| / A_{R'}(T_1)`` for adjacent radii.  The
run is classified ``wall_dependent`` when ``S`` on the largest pair exceeds
``0.20`` or the ``argmax`` location moves by more than ten percent in norm;
``wall_effect_small`` when ``S`` decreases with radius and the largest pair is
at or below ``0.05``; ``undecided`` otherwise.  **The classification is a
scientific conclusion, not a gate.**  The acceptance checks gate run health
only: finiteness, energy, circulation, CFL and the core-identity requirements.

Nothing here is a proof, an interval enclosure, or evidence for singularity
formation, and a finite nested family on a uniform grid is not the limit
"wall to infinity".
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
from typing import Any, Sequence

import numpy as np

from experiments.run_hou_early_time import (
    CFL_EXCESS_MECHANISM,
    CFL_RULE,
    DEFAULT_CFL_EXCESS_TOLERANCE,
    EFFECTIVE_CFL_DEFINITION,
    E29B_MAX_ABS_U1,
    E29B_MAX_CARTESIAN_VORTICITY,
    V1_ADVECTIVE_CFL_EXCESS_RATIO,
    V1_CFL_COEFFICIENT,
    V1_MAXIMUM_ADVECTIVE_CFL,
    advective_cfl_within_tolerance,
)
from experiments.run_hou_time_refinement import cross_check_state
from ns_certificate_lab._integrity import strict_json_loads, write_with_digest
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    DIAGNOSTIC_FIELDS,
    RELATIVE_DIAGNOSTIC_FIELDS,
    ConstrainedState,
    cartesian_vorticity,
    constrain_state,
    integrate,
    kinetic_energy,
    normalize_viscosity_schedule,
    save_checkpoint,
    solve_poisson,
    viscosity_at,
)
from ns_certificate_lab.operators import derivative_z
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.wall_dependence import (
    E32_SUP_DEVIATION_BOUND,
    core_index_count,
    envelope_initial_swirl,
    initial_data_acceptance,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPERIMENT_ID = "wall_dependence_v1"

#: Preregistration section 4 thresholds.  These are fixed by
#: ``docs/wall_dependence_prereg.md`` and must not be tuned.
PREREGISTERED_DECISION_KEYS: tuple[str, ...] = (
    "wall_dependent_amplification_ratio",
    "wall_dependent_argmax_relative_displacement",
    "wall_effect_small_amplification_ratio",
)
#: Thresholds the preregistration leaves unfixed; chosen here and disclosed.
IMPLEMENTATION_DECISION_KEYS: tuple[str, ...] = (
    "resolution_consistency_relative_tolerance",
)

CLASSIFICATIONS: tuple[str, ...] = (
    "wall_dependent",
    "wall_effect_small",
    "undecided",
)

#: Preregistration section 2 check (iv): the E-32 envelope must leave the E-29b
#: derived norms unchanged to this relative tolerance.  Preregistered, not a
#: config knob.  The comparison partner is the plain E-29 datum *on the same
#: grid*, not the continuum E-29b constant: the constants are attained at
#: ``r = 0`` and ``r = 1/sqrt(37)``, which a uniform grid of this size resolves
#: only to about ``1e-4``, whereas the envelope perturbs nothing there at all.
PREREG_DERIVED_NORM_RELATIVE_TOLERANCE = 1.0e-12

# ``docs/wall_dependence_prereg.md`` section 6, reproduced verbatim.  The
# preregistration is written in Japanese; the English renderings that follow in
# :data:`LIMITATIONS` are translations of these four bullets and add nothing.
PREREG_SECTION_6_VERBATIM: tuple[str, ...] = (
    "z 周期 1 固定のまま半径だけ拡大するため、領域の aspect 比が変わる。"
    "完全な「壁を無限遠へ」の極限ではない。",
    "一様格子の core 解像度は Hou の適応格子に遠く及ばず、本実験が測るのは"
    "「この離散化・この時間区間における壁半径感度」である。",
    "T_1 までの早期区間に限定する。後期(粘性切替後)の壁依存性は別実験。",
    "外側境界 psi_1(R_wall)=0 は依然として人工的な有限領域条件であり、"
    "全空間 Green 関数処理(handoff 8.3)の代替ではない。",
)

LIMITATIONS: tuple[str, ...] = (
    "The z period stays fixed at one while only the radius grows, so the "
    "aspect ratio of the domain changes across the family. This is not the "
    "complete 'wall to infinity' limit.",
    "The uniform-grid core resolution is far below the adaptive mesh of the "
    "source calculation, so what this experiment measures is the wall-radius "
    "sensitivity of THIS discretization over THIS time interval.",
    "Only the early interval up to T_1 is integrated. Wall dependence after "
    "the viscosity switch is a separate experiment.",
    "The outer boundary condition psi_1(R_wall) = 0 is still an artificial "
    "finite-domain condition and is not a substitute for a whole-space Green "
    "function treatment.",
    "The E-32 family is not the E-29 datum: it agrees with E-29 bitwise only "
    "on r <= 0.9 and deviates by up to 3.4e-10 in the transition band, so "
    "growth measured here is not a reproduction of any published value.",
    "Floating-point arithmetic only: there is no interval enclosure and no "
    "certified bound anywhere in this pipeline.",
)

KNOWN_GAPS: tuple[str, ...] = (
    "The solver-B snapshot cross-check reuses the public cross_check_state of "
    "experiments/run_hou_time_refinement.py, which itself duplicates the "
    "private _cross_check_snapshot of experiments/run_hou_early_time.py; the "
    "three must be kept in step by hand.",
    "experiments/run_hou_early_time.initial_norms measures the plain E-29 "
    "datum, so the E-32 envelope norms are measured by a local reimplementation "
    "(measure_initial_state) with the same semantics.",
    "The trapezoidal radial weights of the E-20a energy measure are private to "
    "ns_certificate_lab.nonlinear_cylinder (_radial_weights, "
    "_weighted_squared_norm), so the Omega_c-restricted core energy "
    "reimplements them locally in core_energy.",
    "experiments/run_hou_early_time.build_grid reads a single top-level r_max, "
    "while this experiment varies the wall radius per member, so the grid "
    "construction is reimplemented in build_member_grid.",
)

# The core-restricted snapshot columns, in the order of preregistration
# section 3.
CORE_SNAPSHOT_FIELDS: tuple[str, ...] = (
    "group",
    "wall_radius",
    "nr",
    "nz",
    "snapshot_index",
    "time",
    "checkpoint",
    "core_max_cartesian_vorticity",
    "core_amplification",
    "core_amplification_denominator",
    "core_argmax_u1_r",
    "core_argmax_u1_z",
    "core_argmax_abs_u1_r",
    "core_argmax_abs_u1_z",
    "core_max_abs_u1",
    "core_max_abs_psi1",
    "core_max_abs_psi1_z",
    "core_energy",
    "tail_max_abs_u1",
    "tail_max_abs_omega1",
    "tail_row_count",
    "global_max_cartesian_vorticity",
    "global_max_abs_u1",
    "global_max_abs_psi1",
    "psi_cross_solver_max_abs_difference",
    "psi_cross_solver_relative_difference",
    "psi_cross_solver_relative_denominator",
    "psi_cross_solver_argmax_r",
    "psi_cross_solver_argmax_z",
    "solver_b_boundary_error_max",
    "solver_b_physical_defect_max",
)

DENOMINATOR_DEFINITIONS: dict[str, str] = {
    "core_amplification": (
        "the shared initial norm ||omega(0)||_inf of the member group, "
        "recorded as core_amplification_denominator; every member of a group "
        "starts from the bit-identical core datum, so one denominator serves "
        "the whole group and the amplifications are directly comparable"
    ),
    "elliptic_nonlocal_relative_difference": (
        "max |psi1| over Omega_c of the LARGER-radius member for the same "
        "core-restricted source (key relative_denominator)"
    ),
    "amplification_separation": (
        "A_{R'}(T_1), the amplification of the LARGER radius of the pair (key "
        "denominator)"
    ),
    "argmax_relative_displacement": (
        "the Euclidean norm of the LARGER-radius member's argmax location "
        "(r, z) with z measured on the periodic circle (key denominator)"
    ),
    "snapshot_cross_solver_psi_relative_difference": (
        "max |psi1| of solver A at the same snapshot (key "
        "psi_cross_solver_relative_denominator)"
    ),
}


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


def radius_slug(radius: float) -> str:
    """Return a filename-safe token for a wall radius (``1.5 -> R1p5``)."""

    text = f"{float(radius):g}".replace(".", "p").replace("-", "m")
    return f"R{text}"


def circular_gap(first: float, second: float, period: float) -> float:
    """Return the shortest distance between two ``z`` values on the circle."""

    raw = abs(float(first) - float(second)) % float(period)
    return min(raw, float(period) - raw)


# ------------------------------------------------------------ member specification


@dataclass(frozen=True)
class MemberSpec:
    """One wall radius at one core spacing."""

    group: str
    wall_radius: float
    points_per_unit_radius: int
    nr: int
    nz: int
    primary: bool

    @property
    def dr(self) -> float:
        return self.wall_radius / (self.nr - 1)

    @property
    def label(self) -> str:
        return f"{self.group}/{radius_slug(self.wall_radius)}"

    def as_dict(self) -> dict[str, Any]:
        return {
            "group": self.group,
            "wall_radius": self.wall_radius,
            "points_per_unit_radius": self.points_per_unit_radius,
            "nr": self.nr,
            "nz": self.nz,
            "primary": self.primary,
        }


def radial_point_count(points_per_unit_radius: int, wall_radius: float) -> int:
    """Return ``nr = points_per_unit_radius * R_wall + 1``, refusing non-integers.

    Fixing the core spacing across the family is the whole point of the design:
    a radius whose ``points_per_unit_radius * R_wall`` is not an integer cannot
    place ``r = R_wall`` on the shared lattice and is rejected rather than
    rounded.
    """

    count = _positive_integer(
        points_per_unit_radius, name="points_per_unit_radius"
    )
    radius = _number(wall_radius, name="wall_radius", positive=True)
    exact = count * radius
    nearest = round(exact)
    if abs(exact - nearest) > 1.0e-9 * max(1.0, abs(exact)):
        raise ValueError(
            f"wall radius {radius!r} does not yield an integer radial interval "
            f"count at {count} points per unit radius (got {exact!r})"
        )
    return int(nearest) + 1


def build_members(config: dict[str, Any]) -> tuple[MemberSpec, ...]:
    """Expand the configured member groups into individual member specs."""

    members: list[MemberSpec] = []
    for group in config["member_groups"]:
        label = str(group["label"])
        per_unit = _positive_integer(
            group["points_per_unit_radius"], name="points_per_unit_radius"
        )
        nz = _positive_integer(group["nz"], name="nz")
        primary = bool(group["primary"])
        for radius in group["wall_radii"]:
            value = _number(radius, name="wall_radii entry", positive=True)
            members.append(
                MemberSpec(
                    group=label,
                    wall_radius=value,
                    points_per_unit_radius=per_unit,
                    nr=radial_point_count(per_unit, value),
                    nz=nz,
                    primary=primary,
                )
            )
    return tuple(members)


def build_member_grid(config: dict[str, Any], spec: MemberSpec) -> AxisymmetricGrid:
    """Build one member's grid.

    ``run_hou_early_time.build_grid`` reads a single top-level ``r_max``; this
    experiment varies the wall radius per member, so the construction is local
    (see :data:`KNOWN_GAPS`).
    """

    return AxisymmetricGrid.uniform(
        nr=spec.nr,
        nz=spec.nz,
        r_max=spec.wall_radius,
        z_min=0.0,
        z_max=float(config["z_period"]),
        periodic_z=True,
    )


def unit_reference_grid(config: dict[str, Any], spec: MemberSpec) -> AxisymmetricGrid:
    """Build the unit-radius grid at a group's core spacing.

    The E-32 initial-data acceptance measures compare against E-29, which is
    defined only on ``r <= 1``, so they are evaluated on this grid rather than
    on a wide member.
    """

    return AxisymmetricGrid.uniform(
        nr=spec.points_per_unit_radius + 1,
        nz=spec.nz,
        r_max=1.0,
        z_min=0.0,
        z_max=float(config["z_period"]),
        periodic_z=True,
    )


# ------------------------------------------------------------- configuration


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
        "z_period",
        "member_groups",
        "core_radius",
        "tail_radius",
        "envelope",
        "amplitude_scale",
        "viscosity_schedule",
        "t_final",
        "snapshot_times",
        "cfl_coefficient",
        "max_time_step",
        "max_steps",
        "diagnostic_stride",
        "decision_rule",
        "acceptance",
    }
    optional = {"cfl_excess_tolerance"}
    if (
        not isinstance(config, dict)
        or not required <= set(config)
        or not set(config) <= required | optional
    ):
        raise ValueError("wall dependence config has missing or unknown keys")
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
    _number(config["z_period"], name="z_period", positive=True)

    validate_member_groups(config["member_groups"])

    core_radius = _number(config["core_radius"], name="core_radius", positive=True)
    tail_radius = _number(config["tail_radius"], name="tail_radius", positive=True)
    if tail_radius <= core_radius:
        raise ValueError("tail_radius must exceed core_radius")

    envelope = config["envelope"]
    if not isinstance(envelope, dict) or set(envelope) != {"rho1", "rho2"}:
        raise ValueError("envelope must be an object with exactly rho1 and rho2")
    rho1 = _number(envelope["rho1"], name="envelope.rho1")
    rho2 = _number(envelope["rho2"], name="envelope.rho2")
    if rho1 < 0.0:
        raise ValueError("envelope.rho1 must be nonnegative")
    if rho1 >= rho2:
        raise ValueError("envelope.rho1 must be strictly below envelope.rho2")
    if rho1 < core_radius * core_radius:
        raise ValueError(
            "envelope.rho1 must be at least core_radius^2 so that the E-32 "
            "cutoff is exactly one on the whole core"
        )
    smallest = min(
        _number(radius, name="wall_radii entry", positive=True)
        for group in config["member_groups"]
        for radius in group["wall_radii"]
    )
    if rho2 > smallest * smallest:
        raise ValueError(
            "envelope.rho2 must not exceed the smallest wall radius squared: "
            "the E-32 datum must vanish strictly inside every wall"
        )

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

    decision = config["decision_rule"]
    expected_decision = set(PREREGISTERED_DECISION_KEYS) | set(
        IMPLEMENTATION_DECISION_KEYS
    )
    if not isinstance(decision, dict) or set(decision) != expected_decision:
        raise ValueError("decision_rule object has missing or unknown keys")
    for key in expected_decision:
        _number(decision[key], name=f"decision_rule.{key}", positive=True)
    if (
        float(decision["wall_effect_small_amplification_ratio"])
        >= float(decision["wall_dependent_amplification_ratio"])
    ):
        raise ValueError(
            "decision_rule.wall_effect_small_amplification_ratio must be below "
            "decision_rule.wall_dependent_amplification_ratio"
        )

    acceptance = config["acceptance"]
    expected_acceptance = {
        "maximum_circulation_growth_ratio",
        "maximum_energy_growth_ratio",
        "maximum_initial_norm_relative_error",
        "maximum_odd_symmetry_defect_ratio",
        "maximum_core_reference_relative_spread",
        "maximum_initial_axis_parity_relative",
    }
    if not isinstance(acceptance, dict) or set(acceptance) != expected_acceptance:
        raise ValueError("acceptance object has missing or unknown keys")
    for key in expected_acceptance:
        _number(acceptance[key], name=f"acceptance.{key}", positive=True)

    # Expanding the members validates the integer radial counts of every group.
    build_members(config)


def validate_member_groups(groups: object) -> None:
    """Validate the member-group block, including the integer ``nr`` rule."""

    if not isinstance(groups, list) or not groups:
        raise ValueError("member_groups must be a nonempty list")
    labels: set[str] = set()
    primary_count = 0
    for group in groups:
        if not isinstance(group, dict) or set(group) != {
            "label",
            "points_per_unit_radius",
            "nz",
            "wall_radii",
            "primary",
        }:
            raise ValueError("each member group has missing or unknown keys")
        label = group["label"]
        if not isinstance(label, str) or not label.strip():
            raise ValueError("member group label must be a nonempty string")
        if label in labels:
            raise ValueError("member group labels must be unique")
        labels.add(label)
        _positive_integer(
            group["points_per_unit_radius"], name="points_per_unit_radius"
        )
        nz = _positive_integer(group["nz"], name="nz")
        if nz < 5:
            raise ValueError("nz is below the AxisymmetricGrid minimum")
        if not isinstance(group["primary"], bool):
            raise ValueError("member group primary must be a boolean")
        primary_count += 1 if group["primary"] else 0
        radii = group["wall_radii"]
        if not isinstance(radii, list) or len(radii) < 2:
            raise ValueError("each member group needs at least two wall radii")
        values = [
            _number(item, name="wall_radii entry", positive=True) for item in radii
        ]
        if any(a >= b for a, b in zip(values, values[1:])):
            raise ValueError("wall radii must be strictly increasing")
        for value in values:
            nr = radial_point_count(int(group["points_per_unit_radius"]), value)
            if nr < 4:
                raise ValueError("nr is below the AxisymmetricGrid minimum")
    if primary_count != 1:
        raise ValueError("exactly one member group must be marked primary")


def cfl_excess_tolerance_of(config: dict[str, Any]) -> float:
    """Return the configured CFL excess tolerance or its documented default."""

    if "cfl_excess_tolerance" not in config:
        return DEFAULT_CFL_EXCESS_TOLERANCE
    return _number(config["cfl_excess_tolerance"], name="cfl_excess_tolerance")


# ------------------------------------------------------- core-restricted measures


def core_energy(
    grid: AxisymmetricGrid,
    state: ConstrainedState,
    *,
    core_rows: int,
) -> float:
    r"""Return the E-20a kinetic energy restricted to ``Omega_c``.

    ``E = pi \int ((u^r)^2 + (u^theta)^2 + (u^z)^2) r dr dz`` with the same
    trapezoidal radial weighting as
    :func:`ns_certificate_lab.nonlinear_cylinder.kinetic_energy`, halving the
    two end rows of the integration range.  Those weights are module private
    there, so they are reimplemented here (see :data:`KNOWN_GAPS`).
    """

    if core_rows < 2 or core_rows > grid.nr:
        raise ValueError("core_rows must select at least two radial rows")
    radius = grid.r[:core_rows, None]
    total = (
        state.u_r[:core_rows] ** 2
        + (radius * state.u1[:core_rows]) ** 2
        + state.u_z[:core_rows] ** 2
    )
    weights = grid.r[:core_rows].copy()
    weights[0] *= 0.5
    weights[-1] *= 0.5
    return float(
        math.pi * grid.dr * grid.dz * np.sum(weights[:, None] * total)
    )


def measure_initial_state(
    grid: AxisymmetricGrid,
    *,
    amplitude_scale: float,
    rho1: float,
    rho2: float,
    core_rows: int,
) -> dict[str, float]:
    """Measure the E-32 initial datum's norms on one member grid.

    ``run_hou_early_time.initial_norms`` measures the plain E-29 datum, so the
    envelope needs its own measurement (see :data:`KNOWN_GAPS`).
    """

    u1 = envelope_initial_swirl(
        grid,
        amplitude_scale=amplitude_scale,
        rho1=rho1,
        rho2=rho2,
    )
    state = constrain_state(grid, u1, np.zeros(grid.shape))
    magnitude = np.sqrt(
        sum(
            component * component
            for component in cartesian_vorticity(grid, state.u1, state.omega1)
        )
    )
    full = float(np.max(magnitude))
    core = float(np.max(magnitude[:core_rows]))
    measured_u1 = float(np.max(np.abs(u1)))
    expected_u1 = E29B_MAX_ABS_U1 * float(amplitude_scale)
    expected_vorticity = E29B_MAX_CARTESIAN_VORTICITY * float(amplitude_scale)
    return {
        "max_cartesian_vorticity_full": full,
        "max_cartesian_vorticity_core": core,
        "core_equals_full": float(full == core),
        "max_abs_u1_full": measured_u1,
        "max_abs_u1_core": float(np.max(np.abs(u1[:core_rows]))),
        "energy": kinetic_energy(grid, state),
        "core_energy": core_energy(grid, state, core_rows=core_rows),
        "e29b_max_abs_u1_relative_error": abs(measured_u1 - expected_u1)
        / expected_u1,
        "e29b_max_cartesian_vorticity_relative_error": abs(
            full - expected_vorticity
        )
        / expected_vorticity,
    }


@dataclass(frozen=True)
class CoreSnapshot:
    """One checkpointed snapshot reduced to the preregistered core measures."""

    index: int
    time: float
    checkpoint: str
    core_max_cartesian_vorticity: float
    core_amplification: float
    core_amplification_denominator: float
    core_argmax_u1_r: float
    core_argmax_u1_z: float
    core_argmax_abs_u1_r: float
    core_argmax_abs_u1_z: float
    core_max_abs_u1: float
    core_max_abs_psi1: float
    core_max_abs_psi1_z: float
    core_energy: float
    tail_max_abs_u1: float | None
    tail_max_abs_omega1: float | None
    tail_row_count: int
    global_max_cartesian_vorticity: float
    global_max_abs_u1: float
    global_max_abs_psi1: float
    cross: dict[str, float]

    def as_row(self, spec: MemberSpec) -> dict[str, Any]:
        return {
            "group": spec.group,
            "wall_radius": spec.wall_radius,
            "nr": spec.nr,
            "nz": spec.nz,
            "snapshot_index": self.index,
            "time": self.time,
            "checkpoint": self.checkpoint,
            "core_max_cartesian_vorticity": self.core_max_cartesian_vorticity,
            "core_amplification": self.core_amplification,
            "core_amplification_denominator": self.core_amplification_denominator,
            "core_argmax_u1_r": self.core_argmax_u1_r,
            "core_argmax_u1_z": self.core_argmax_u1_z,
            "core_argmax_abs_u1_r": self.core_argmax_abs_u1_r,
            "core_argmax_abs_u1_z": self.core_argmax_abs_u1_z,
            "core_max_abs_u1": self.core_max_abs_u1,
            "core_max_abs_psi1": self.core_max_abs_psi1,
            "core_max_abs_psi1_z": self.core_max_abs_psi1_z,
            "core_energy": self.core_energy,
            "tail_max_abs_u1": self.tail_max_abs_u1,
            "tail_max_abs_omega1": self.tail_max_abs_omega1,
            "tail_row_count": self.tail_row_count,
            "global_max_cartesian_vorticity": self.global_max_cartesian_vorticity,
            "global_max_abs_u1": self.global_max_abs_u1,
            "global_max_abs_psi1": self.global_max_abs_psi1,
            "psi_cross_solver_max_abs_difference": self.cross["max_abs_difference"],
            "psi_cross_solver_relative_difference": self.cross[
                "relative_difference"
            ],
            "psi_cross_solver_relative_denominator": self.cross[
                "relative_denominator"
            ],
            "psi_cross_solver_argmax_r": self.cross["argmax_r"],
            "psi_cross_solver_argmax_z": self.cross["argmax_z"],
            "solver_b_boundary_error_max": self.cross["boundary_error_max"],
            "solver_b_physical_defect_max": self.cross["physical_defect_max"],
        }


def core_snapshot_measures(
    grid: AxisymmetricGrid,
    state: ConstrainedState,
    *,
    index: int,
    time: float,
    checkpoint: str,
    core_rows: int,
    tail_radius: float,
    reference_vorticity: float,
    wall_radius: float,
) -> CoreSnapshot:
    """Reduce one snapshot to the preregistered section 3 quantities."""

    magnitude = np.sqrt(
        sum(
            component * component
            for component in cartesian_vorticity(grid, state.u1, state.omega1)
        )
    )
    core_magnitude = magnitude[:core_rows]
    core_u1 = state.u1[:core_rows]
    core_psi1 = state.psi1[:core_rows]
    psi1_z = derivative_z(grid, state.psi1)

    row, column = divmod(int(np.argmax(core_u1)), grid.nz)
    abs_row, abs_column = divmod(int(np.argmax(np.abs(core_u1))), grid.nz)

    tail = grid.r > tail_radius
    tail_rows = int(np.count_nonzero(tail))
    tail_u1: float | None = None
    tail_omega1: float | None = None
    if wall_radius > 1.0 and tail_rows > 0:
        tail_u1 = float(np.max(np.abs(state.u1[tail])))
        tail_omega1 = float(np.max(np.abs(state.omega1[tail])))

    denominator = float(reference_vorticity)
    core_max = float(np.max(core_magnitude))
    return CoreSnapshot(
        index=index,
        time=float(time),
        checkpoint=checkpoint,
        core_max_cartesian_vorticity=core_max,
        core_amplification=core_max / max(denominator, 1.0e-300),
        core_amplification_denominator=denominator,
        core_argmax_u1_r=float(grid.r[row]),
        core_argmax_u1_z=float(grid.z[column]),
        core_argmax_abs_u1_r=float(grid.r[abs_row]),
        core_argmax_abs_u1_z=float(grid.z[abs_column]),
        core_max_abs_u1=float(np.max(np.abs(core_u1))),
        core_max_abs_psi1=float(np.max(np.abs(core_psi1))),
        core_max_abs_psi1_z=float(np.max(np.abs(psi1_z[:core_rows]))),
        core_energy=core_energy(grid, state, core_rows=core_rows),
        tail_max_abs_u1=tail_u1,
        tail_max_abs_omega1=tail_omega1,
        tail_row_count=tail_rows,
        global_max_cartesian_vorticity=float(np.max(magnitude)),
        global_max_abs_u1=float(np.max(np.abs(state.u1))),
        global_max_abs_psi1=float(np.max(np.abs(state.psi1))),
        cross=cross_check_state(grid, state),
    )


# ------------------------------------------------------------------ one member


@dataclass(frozen=True)
class MemberResult:
    """Outcome of one member, including an explicit failure slot."""

    spec: MemberSpec
    grid: AxisymmetricGrid
    core_rows: int
    completed: bool
    failure: str | None
    step_count: int
    final_time: float
    history: tuple[dict[str, float], ...]
    snapshots: tuple[CoreSnapshot, ...]
    core_omega1: tuple[np.ndarray, ...]
    initial: dict[str, float]
    reference_vorticity: float
    state: ConstrainedState | None = dataclass_field(default=None)


def _series(history: Sequence[dict[str, float]], name: str) -> list[float]:
    return [record[name] for record in history]


def evolve_member(
    config: dict[str, Any],
    spec: MemberSpec,
    *,
    reference_vorticity: float,
    checkpoint_dir: Path | None = None,
    provenance: dict[str, Any] | None = None,
) -> MemberResult:
    """Integrate one member, checkpointing and reducing every snapshot."""

    validate_config(config)
    grid = build_member_grid(config, spec)
    core_radius = float(config["core_radius"])
    tail_radius = float(config["tail_radius"])
    core_rows = core_index_count(grid, core_radius=core_radius)
    amplitude_scale = float(config["amplitude_scale"])
    rho1 = float(config["envelope"]["rho1"])
    rho2 = float(config["envelope"]["rho2"])
    schedule = normalize_viscosity_schedule(config["viscosity_schedule"])
    u1_initial = envelope_initial_swirl(
        grid,
        amplitude_scale=amplitude_scale,
        rho1=rho1,
        rho2=rho2,
    )
    initial = measure_initial_state(
        grid,
        amplitude_scale=amplitude_scale,
        rho1=rho1,
        rho2=rho2,
        core_rows=core_rows,
    )

    snapshots: list[CoreSnapshot] = []
    core_sources: list[np.ndarray] = []

    def on_snapshot(
        time: float,
        state: ConstrainedState,
        diagnostics: dict[str, float],
    ) -> None:
        index = len(snapshots)
        name = ""
        if checkpoint_dir is not None:
            name = (
                f"checkpoint_{spec.group}_{radius_slug(spec.wall_radius)}"
                f"_t{index:03d}.npz"
            )
        record = core_snapshot_measures(
            grid,
            state,
            index=index,
            time=time,
            checkpoint=name,
            core_rows=core_rows,
            tail_radius=tail_radius,
            reference_vorticity=reference_vorticity,
            wall_radius=spec.wall_radius,
        )
        if checkpoint_dir is not None:
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
                    "group": spec.group,
                    "wall_radius": spec.wall_radius,
                    "nr": spec.nr,
                    "nz": spec.nz,
                    "amplitude_scale": amplitude_scale,
                    "core_radius": core_radius,
                    "core_rows": core_rows,
                    "core_amplification": record.core_amplification,
                    "cross_solver_psi_max_abs_difference": record.cross[
                        "max_abs_difference"
                    ],
                },
            )
        snapshots.append(record)
        core_sources.append(np.array(state.omega1[:core_rows], dtype=np.float64))

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
            reference_vorticity_max=reference_vorticity,
        )
    except (FloatingPointError, ArithmeticError, ValueError) as exc:
        return MemberResult(
            spec=spec,
            grid=grid,
            core_rows=core_rows,
            completed=False,
            failure=f"{type(exc).__name__}: {exc}",
            step_count=0,
            final_time=0.0,
            history=(),
            snapshots=tuple(snapshots),
            core_omega1=tuple(core_sources),
            initial=initial,
            reference_vorticity=float(reference_vorticity),
        )
    if not result.completed:
        failure = (
            f"integration stopped after {result.step_count} steps at "
            f"t={result.time!r} before t_final"
        )
    return MemberResult(
        spec=spec,
        grid=grid,
        core_rows=core_rows,
        completed=bool(result.completed),
        failure=failure,
        step_count=int(result.step_count),
        final_time=float(result.time),
        history=result.history,
        snapshots=tuple(snapshots),
        core_omega1=tuple(core_sources),
        initial=initial,
        reference_vorticity=float(reference_vorticity),
        state=result.state,
    )


def member_metrics(
    result: MemberResult,
    *,
    cfl_coefficient: float,
) -> dict[str, Any]:
    """Reduce one member's history to the recorded acceptance quantities."""

    identity: dict[str, Any] = {
        **result.spec.as_dict(),
        "dr": result.grid.dr,
        "dz": result.grid.dz,
        "core_rows": result.core_rows,
        "completed": result.completed,
        "failure": result.failure,
        "step_count": result.step_count,
        "final_time": result.final_time,
        "reference_vorticity": result.reference_vorticity,
        "initial": dict(result.initial),
        "snapshot_count": len(result.snapshots),
    }
    if not result.history:
        return identity
    history = result.history
    energy = _series(history, "energy")
    circulation = _series(history, "circulation_max")
    swirl = _series(history, "max_abs_u1")
    odd_defect = _series(history, "odd_symmetry_defect")
    maximum_cfl = max(_series(history, "advective_cfl"))
    final_snapshot = result.snapshots[-1] if result.snapshots else None
    return {
        **identity,
        "diagnostic_record_count": len(history),
        "all_diagnostics_finite": all(
            math.isfinite(value) for record in history for value in record.values()
        ),
        "initial_energy": energy[0],
        "final_energy": energy[-1],
        "maximum_energy_growth_ratio": max(energy) / energy[0] - 1.0,
        "initial_circulation_max": circulation[0],
        "final_circulation_max": circulation[-1],
        "maximum_circulation_growth_ratio": max(circulation) / circulation[0] - 1.0,
        "initial_max_abs_u1": swirl[0],
        "final_max_abs_u1": swirl[-1],
        "final_global_amplification": history[-1]["amplification"],
        "maximum_global_amplification": max(_series(history, "amplification")),
        "maximum_odd_symmetry_defect": max(odd_defect),
        "maximum_odd_symmetry_defect_ratio": max(odd_defect) / max(swirl),
        "maximum_axis_parity_defect": max(_series(history, "axis_parity_defect")),
        "maximum_divergence_residual": max(
            _series(history, "divergence_residual_max")
        ),
        "maximum_divergence_residual_relative": max(
            _series(history, "divergence_residual_relative")
        ),
        "maximum_divergence_pointwise_ratio": max(
            _series(history, "divergence_pointwise_ratio_max")
        ),
        "maximum_axis_parity_relative_u1": max(
            _series(history, "axis_parity_relative_u1")
        ),
        "maximum_axis_parity_relative_omega1": max(
            _series(history, "axis_parity_relative_omega1")
        ),
        "maximum_wall_u1_abs": max(_series(history, "wall_u1_max_abs")),
        "minimum_dt": min(record["dt"] for record in history[1:])
        if len(history) > 1
        else 0.0,
        "maximum_advective_cfl": maximum_cfl,
        "maximum_advective_cfl_excess_ratio": maximum_cfl / cfl_coefficient - 1.0,
        "maximum_viscous_cfl": max(_series(history, "viscous_cfl")),
        "final_core_amplification": (
            None if final_snapshot is None else final_snapshot.core_amplification
        ),
        "final_core_max_cartesian_vorticity": (
            None
            if final_snapshot is None
            else final_snapshot.core_max_cartesian_vorticity
        ),
        "final_core_argmax_u1_r": (
            None if final_snapshot is None else final_snapshot.core_argmax_u1_r
        ),
        "final_core_argmax_u1_z": (
            None if final_snapshot is None else final_snapshot.core_argmax_u1_z
        ),
        "final_core_max_abs_psi1_z": (
            None if final_snapshot is None else final_snapshot.core_max_abs_psi1_z
        ),
        "final_core_energy": (
            None if final_snapshot is None else final_snapshot.core_energy
        ),
        "final_tail_max_abs_u1": (
            None if final_snapshot is None else final_snapshot.tail_max_abs_u1
        ),
        "final_tail_max_abs_omega1": (
            None if final_snapshot is None else final_snapshot.tail_max_abs_omega1
        ),
        "snapshot_cross_solver_psi_max_abs_difference": (
            max(item.cross["max_abs_difference"] for item in result.snapshots)
            if result.snapshots
            else None
        ),
        "snapshot_cross_solver_psi_relative_difference": (
            max(item.cross["relative_difference"] for item in result.snapshots)
            if result.snapshots
            else None
        ),
    }


# ------------------------------------------------------------- core identity


def core_identity_report(
    grids: Sequence[AxisymmetricGrid],
    specs: Sequence[MemberSpec],
    *,
    core_rows: int,
) -> dict[str, Any]:
    """Verify that a group's members really do share one core discretization.

    The requirement is exact, not approximate: ``np.linspace(0, R, nr)``
    evaluates ``r_j = j * (R / (nr - 1))`` and ``R / (nr - 1)`` is the identical
    float for every member (``1/192 == 1.5/288 == 2/384 == 3/576`` bitwise), so
    the shared core nodes must compare bitwise equal.  Nodewise comparison of
    the elliptic solutions on ``Omega_c`` needs no interpolation *because* of
    this, and the check is what licenses that.
    """

    if not grids:
        raise ValueError("a member group must contain at least one grid")
    spacings = [item.dr for item in grids]
    reference = grids[0]
    nodes_identical = all(
        np.array_equal(item.r[:core_rows], reference.r[:core_rows])
        for item in grids
    )
    axial_identical = all(np.array_equal(item.z, reference.z) for item in grids)
    spread = max(spacings) - min(spacings)
    return {
        "members": [item.label for item in specs],
        "radial_spacings": spacings,
        "radial_spacing_identical": bool(len(set(spacings)) == 1),
        "radial_spacing_absolute_spread": float(spread),
        "radial_spacing_relative_spread": float(spread / spacings[0]),
        "axial_spacings": [item.dz for item in grids],
        "axial_grid_identical": bool(axial_identical),
        "core_rows": int(core_rows),
        "core_radial_nodes_bitwise_identical": bool(nodes_identical),
        "core_radius_of_last_core_node": float(reference.r[core_rows - 1]),
        "radial_point_counts": [item.nr for item in specs],
        "radial_point_counts_are_integers": True,
        "interpretation": (
            "the shared core nodes are bitwise identical, so the Omega_c "
            "comparisons in this experiment are exact nodewise comparisons "
            "with no interpolation anywhere"
        ),
    }


def core_reference_report(
    results_or_initials: Sequence[dict[str, float]],
    specs: Sequence[MemberSpec],
) -> dict[str, Any]:
    """Report the shared ``||omega(0)||_inf`` and its spread across a group."""

    core_values = [
        float(item["max_cartesian_vorticity_core"]) for item in results_or_initials
    ]
    full_values = [
        float(item["max_cartesian_vorticity_full"]) for item in results_or_initials
    ]
    core_energies = [float(item["core_energy"]) for item in results_or_initials]
    reference = core_values[0]
    spread = max(core_values) - min(core_values)
    energy_spread = max(core_energies) - min(core_energies)
    return {
        "members": [item.label for item in specs],
        "reference": reference,
        "reference_definition": (
            "||omega(0)||_inf of the first member of the group; the initial "
            "datum is compactly supported inside r < 0.95 and its maximum sits "
            "at r = 1/sqrt(37), so the core and whole-domain norms coincide"
        ),
        "initial_core_maxima": core_values,
        "initial_full_maxima": full_values,
        "core_equals_full_for_every_member": bool(
            all(a == b for a, b in zip(core_values, full_values))
        ),
        "absolute_spread": float(spread),
        "relative_spread": float(spread / reference) if reference > 0.0 else 0.0,
        "bitwise_identical": bool(len(set(core_values)) == 1),
        "initial_core_energies": core_energies,
        "initial_core_energy_relative_spread": (
            float(energy_spread / core_energies[0]) if core_energies[0] > 0.0 else 0.0
        ),
    }


# ------------------------------------------------- elliptic nonlocal contribution


def elliptic_nonlocal_pair(
    core_source: np.ndarray,
    *,
    small_grid: AxisymmetricGrid,
    large_grid: AxisymmetricGrid,
    core_rows: int,
) -> dict[str, float]:
    r"""Solve one core-restricted ``omega_1`` on two wall radii and compare.

    Preregistration item 5.  The *same* core-restricted source (the
    larger-radius member's ``omega_1`` zeroed outside ``Omega_c``) is embedded
    in both domains and solved with the E-27 homogeneous outer Dirichlet
    condition.  Because the two grids share the core nodes bitwise the
    comparison on ``Omega_c`` is exact and needs no interpolation; that sharing
    is asserted here rather than assumed.
    """

    source = np.asarray(core_source, dtype=np.float64)
    if source.shape != (core_rows, small_grid.nz):
        raise ValueError("the core source has the wrong shape")
    if small_grid.nz != large_grid.nz:
        raise ValueError("the two members must share the axial grid")
    if not np.array_equal(
        small_grid.r[:core_rows], large_grid.r[:core_rows]
    ) or not np.array_equal(small_grid.z, large_grid.z):
        raise ValueError(
            "the two members do not share the core nodes bitwise, so a "
            "nodewise Omega_c comparison would be meaningless"
        )
    small_rhs = np.zeros(small_grid.shape, dtype=np.float64)
    small_rhs[:core_rows] = source
    large_rhs = np.zeros(large_grid.shape, dtype=np.float64)
    large_rhs[:core_rows] = source
    psi_small = solve_poisson(small_grid, small_rhs)
    psi_large = solve_poisson(large_grid, large_rhs)
    deviation = np.abs(psi_small[:core_rows] - psi_large[:core_rows])
    row, column = divmod(int(np.argmax(deviation)), small_grid.nz)
    difference = float(deviation[row, column])
    denominator = float(np.max(np.abs(psi_large[:core_rows])))
    return {
        "max_abs_difference": difference,
        "relative_difference": difference / max(denominator, 1.0e-300),
        "relative_denominator": denominator,
        "small_core_max_abs_psi1": float(np.max(np.abs(psi_small[:core_rows]))),
        "large_core_max_abs_psi1": denominator,
        "argmax_r": float(small_grid.r[row]),
        "argmax_z": float(small_grid.z[column]),
        "source_core_max_abs_omega1": float(np.max(np.abs(source))),
    }


def elliptic_nonlocal_report(
    results: Sequence[MemberResult],
    *,
    core_rows: int,
) -> list[dict[str, Any]]:
    """Run :func:`elliptic_nonlocal_pair` on every radius pair and snapshot."""

    report: list[dict[str, Any]] = []
    ordered = sorted(results, key=lambda item: item.spec.wall_radius)
    for index, small in enumerate(ordered):
        for large in ordered[index + 1 :]:
            count = min(len(small.snapshots), len(large.core_omega1))
            for snapshot in range(count):
                identity = {
                    "group": small.spec.group,
                    "smaller_wall_radius": small.spec.wall_radius,
                    "larger_wall_radius": large.spec.wall_radius,
                    "snapshot_index": snapshot,
                    "time": large.snapshots[snapshot].time,
                    "source": (
                        "the larger-radius member's omega1 restricted to "
                        "Omega_c and embedded in both domains"
                    ),
                }
                try:
                    measures = elliptic_nonlocal_pair(
                        large.core_omega1[snapshot],
                        small_grid=small.grid,
                        large_grid=large.grid,
                        core_rows=core_rows,
                    )
                except ValueError as exc:
                    # The nodewise Omega_c comparison is only meaningful when
                    # the two members share the core lattice bitwise.  A pair
                    # that does not is recorded as an explicit failure so the
                    # acceptance check falls over instead of the run crashing.
                    report.append(
                        {
                            **identity,
                            "status": "core_nodes_not_shared",
                            "error": f"{type(exc).__name__}: {exc}",
                        }
                    )
                    continue
                report.append({**identity, "status": "measured", **measures})
    return report


# ------------------------------------------------------------- decision rule


def amplification_separations(
    amplifications: Sequence[float],
    radii: Sequence[float],
) -> list[dict[str, Any]]:
    """Return ``S(R,R') = |A_R - A_R'| / A_R'`` for every adjacent radius pair."""

    report: list[dict[str, Any]] = []
    for index in range(len(radii) - 1):
        smaller = float(amplifications[index])
        larger = float(amplifications[index + 1])
        finite = math.isfinite(smaller) and math.isfinite(larger) and larger > 0.0
        report.append(
            {
                "smaller_wall_radius": float(radii[index]),
                "larger_wall_radius": float(radii[index + 1]),
                "amplification_smaller": smaller,
                "amplification_larger": larger,
                "separation": (
                    abs(smaller - larger) / larger if finite else None
                ),
                "denominator": larger,
                "finite": bool(finite),
            }
        )
    return report


def argmax_displacement(
    first: tuple[float, float],
    second: tuple[float, float],
    *,
    z_period: float,
) -> dict[str, Any]:
    """Return the relative displacement between two ``argmax`` locations.

    ``second`` is the reference (the larger radius): the displacement norm is
    divided by ``||second||``.  The axial separation is measured on the
    periodic circle, so the wrap-around of ``z`` never inflates it.
    """

    radial = float(first[0]) - float(second[0])
    axial = circular_gap(first[1], second[1], z_period)
    displacement = math.hypot(radial, axial)
    denominator = math.hypot(float(second[0]), float(second[1]))
    return {
        "radial_separation": radial,
        "axial_separation": axial,
        "displacement": displacement,
        "denominator": denominator,
        "relative_displacement": (
            displacement / denominator if denominator > 0.0 else None
        ),
        "reference_r": float(second[0]),
        "reference_z": float(second[1]),
    }


def resolution_consistency_report(
    separations_by_group: dict[str, list[dict[str, Any]]],
    amplifications_by_group: dict[str, dict[float, float]],
    *,
    primary_group: str,
    tolerance: float,
) -> dict[str, Any]:
    """Compare the wall-effect measurement across the two core spacings.

    The preregistration's resolution-consistency check asks whether the
    *wall-effect* measurement is stable under refinement, not whether the
    amplification itself is: the shipped resolution ladder of
    ``outputs/hou_early_time_v1`` already shows that ``A`` itself moves a long
    way between ``dr = 1/128`` and ``dr = 1/192``, so demanding stability of
    ``A`` would fail by construction and would say nothing about the wall.  The
    compared statistic is therefore ``S`` on the radius pair the two groups
    share.
    """

    others = [name for name in separations_by_group if name != primary_group]
    shared: list[dict[str, Any]] = []
    stable = True
    for name in others:
        primary_radii = set(amplifications_by_group[primary_group])
        secondary_radii = set(amplifications_by_group[name])
        common = sorted(primary_radii & secondary_radii)
        if len(common) < 2:
            stable = False
            shared.append(
                {
                    "group": name,
                    "status": "insufficient_shared_radii",
                    "shared_radii": common,
                }
            )
            continue
        low, high = common[0], common[-1]
        entry: dict[str, Any] = {"group": name, "shared_radii": [low, high]}
        values: dict[str, float | None] = {}
        for label, source in (
            (primary_group, amplifications_by_group[primary_group]),
            (name, amplifications_by_group[name]),
        ):
            a_low = source.get(low)
            a_high = source.get(high)
            if (
                a_low is None
                or a_high is None
                or not math.isfinite(a_low)
                or not math.isfinite(a_high)
                or a_high <= 0.0
            ):
                values[label] = None
            else:
                values[label] = abs(a_low - a_high) / a_high
        entry["separation_by_group"] = values
        primary_value = values[primary_group]
        secondary_value = values[name]
        if (
            primary_value is None
            or secondary_value is None
            or primary_value <= 0.0
        ):
            entry["relative_change"] = None
            entry["stable"] = False
            stable = False
        else:
            change = abs(secondary_value - primary_value) / primary_value
            entry["relative_change"] = change
            entry["stable"] = bool(change <= tolerance)
            stable = stable and entry["stable"]
        entry["amplification_by_group"] = {
            primary_group: {
                str(radius): amplifications_by_group[primary_group].get(radius)
                for radius in (low, high)
            },
            name: {
                str(radius): amplifications_by_group[name].get(radius)
                for radius in (low, high)
            },
        }
        shared.append(entry)
    return {
        "statistic": (
            "S(R_low, R_high) = |A_low - A_high| / A_high on the radius pair "
            "shared by the two core spacings"
        ),
        "tolerance": float(tolerance),
        "tolerance_source": (
            "not fixed by docs/wall_dependence_prereg.md; chosen by this "
            "implementation and recorded here"
        ),
        "primary_group": primary_group,
        "comparisons": shared,
        "stable_under_refinement": bool(stable and bool(others)),
        "why_not_the_amplification_itself": (
            "the shipped hou_early_time resolution ladder moves the "
            "amplification from 12.70 at dr = 1/128 to 15.63 at dr = 1/192 on "
            "the plain E-29 datum, so the amplification is strongly resolution "
            "limited; only the wall-effect ratio S can be asked to be stable"
        ),
    }


def classify_wall_dependence(
    separations: Sequence[dict[str, Any]],
    displacement: dict[str, Any] | None,
    *,
    thresholds: dict[str, Any],
    resolution_stable: bool,
) -> dict[str, Any]:
    """Apply preregistration section 4 exactly.  RECORD ONLY; never a gate.

    ``wall_dependent`` when ``S`` on the largest adjacent pair exceeds
    ``wall_dependent_amplification_ratio`` or the ``argmax`` location moves by
    more than ``wall_dependent_argmax_relative_displacement`` in norm;
    ``wall_effect_small`` when ``S`` strictly decreases with radius and the
    largest pair is at or below ``wall_effect_small_amplification_ratio``;
    ``undecided`` otherwise.

    The preregistration additionally conditions the ``wall_dependent`` verdict
    on the resolution-consistency check ("both measurements pass the section 2
    resolution-consistency check", and section 2 says to hold the verdict when
    the measurement is unstable).  That qualifier is applied last and the
    unqualified verdict is preserved in ``classification_before_qualifier``.
    """

    values = [item["separation"] for item in separations]
    finite = [
        value for value in values if value is not None and math.isfinite(value)
    ]
    relative_displacement = (
        None if displacement is None else displacement.get("relative_displacement")
    )
    wall_dependent_ratio = float(thresholds["wall_dependent_amplification_ratio"])
    wall_dependent_shift = float(
        thresholds["wall_dependent_argmax_relative_displacement"]
    )
    wall_effect_small_ratio = float(
        thresholds["wall_effect_small_amplification_ratio"]
    )

    reasons: list[str] = []
    if len(finite) != len(values) or not values:
        classification = "undecided"
        reasons.append("at least one adjacent separation is not a finite number")
        strictly_decreasing = False
        non_increasing = False
        largest = None
    else:
        largest = float(values[-1])
        strictly_decreasing = all(a > b for a, b in zip(values, values[1:]))
        non_increasing = all(a >= b for a, b in zip(values, values[1:]))
        exceeds_ratio = largest > wall_dependent_ratio
        exceeds_shift = (
            relative_displacement is not None
            and math.isfinite(relative_displacement)
            and relative_displacement > wall_dependent_shift
        )
        if exceeds_ratio or exceeds_shift:
            classification = "wall_dependent"
            if exceeds_ratio:
                reasons.append(
                    f"S on the largest pair is {largest!r} > {wall_dependent_ratio!r}"
                )
            if exceeds_shift:
                reasons.append(
                    "the argmax location moved by "
                    f"{relative_displacement!r} > {wall_dependent_shift!r} in norm"
                )
        elif strictly_decreasing and largest <= wall_effect_small_ratio:
            classification = "wall_effect_small"
            reasons.append(
                "S decreases strictly with radius and the largest pair is "
                f"{largest!r} <= {wall_effect_small_ratio!r}"
            )
        else:
            classification = "undecided"
            reasons.append(
                "neither the wall-dependent nor the wall-effect-small "
                "condition holds"
            )

    unqualified = classification
    if classification == "wall_dependent" and not resolution_stable:
        classification = "undecided"
        reasons.append(
            "preregistration section 4: the wall_dependent verdict additionally "
            "requires both measurements to pass the section 2 "
            "resolution-consistency check, and they do not, so the verdict is "
            "held"
        )
    # Section 4 attaches the resolution-consistency clause to the
    # wall_dependent branch only, while section 2 says more generally to hold
    # the verdict ("判定保留") when the measurement is unstable.  Both readings
    # are recorded so neither is lost and neither threshold is touched.
    strict_hold = (
        "undecided"
        if (classification != "undecided" and not resolution_stable)
        else classification
    )
    return {
        "classification": classification,
        "classification_before_qualifier": unqualified,
        "classification_with_section_2_hold": strict_hold,
        "qualifier_readings": (
            "classification_before_qualifier applies the section 4 thresholds "
            "alone; classification additionally applies the section 4 clause "
            "that a wall_dependent verdict needs the resolution-consistency "
            "check; classification_with_section_2_hold applies the more "
            "general section 2 instruction to hold ANY decisive verdict when "
            "the wall-effect measurement is not stable under refinement"
        ),
        "allowed_classifications": list(CLASSIFICATIONS),
        "separations": [item["separation"] for item in separations],
        "largest_pair_separation": largest,
        "separations_strictly_decreasing": bool(strictly_decreasing),
        "separations_non_increasing": bool(non_increasing),
        "monotonicity_convention": (
            "the wall_effect_small branch requires STRICT decrease; the "
            "non-strict boolean is recorded alongside so the weaker reading is "
            "auditable"
        ),
        "argmax_relative_displacement": relative_displacement,
        "argmax_displacement": displacement,
        "thresholds": {
            "wall_dependent_amplification_ratio": wall_dependent_ratio,
            "wall_dependent_argmax_relative_displacement": wall_dependent_shift,
            "wall_effect_small_amplification_ratio": wall_effect_small_ratio,
            "source": "docs/wall_dependence_prereg.md section 4 (preregistered)",
        },
        "resolution_consistency_stable": bool(resolution_stable),
        "reasons": reasons,
        "gated": False,
        "interpretation": (
            "This classification is a recorded scientific conclusion about one "
            "finite family on one uniform discretization over one early time "
            "interval. It is not an acceptance gate, it is not a proof of wall "
            "independence in the whole space, and a wall_effect_small outcome "
            "is explicitly NOT evidence of a Clay candidate."
        ),
        "denominator_definitions": {
            "separation": DENOMINATOR_DEFINITIONS["amplification_separation"],
            "argmax_relative_displacement": DENOMINATOR_DEFINITIONS[
                "argmax_relative_displacement"
            ],
        },
    }


# ---------------------------------------------------------------- evaluation


def cfl_policy(
    config: dict[str, Any],
    metrics: Sequence[dict[str, Any]],
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
        "v1_reference_measurement": {
            "run": "outputs/hou_early_time_v1 (193x384 resolution)",
            "cfl_coefficient": V1_CFL_COEFFICIENT,
            "maximum_advective_cfl": V1_MAXIMUM_ADVECTIVE_CFL,
            "excess_ratio": V1_ADVECTIVE_CFL_EXCESS_RATIO,
        },
    }


def evaluate(
    config: dict[str, Any],
    *,
    checkpoint_dir: Path | None = None,
    provenance: dict[str, Any] | None = None,
) -> tuple[
    list[MemberResult],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, bool],
]:
    """Run every member and reduce the family to metrics, analysis and checks."""

    validate_config(config)
    acceptance = config["acceptance"]
    coefficient = float(config["cfl_coefficient"])
    tolerance = cfl_excess_tolerance_of(config)
    core_radius = float(config["core_radius"])
    amplitude_scale = float(config["amplitude_scale"])
    rho1 = float(config["envelope"]["rho1"])
    rho2 = float(config["envelope"]["rho2"])
    specs = build_members(config)

    groups: dict[str, list[MemberSpec]] = {}
    for spec in specs:
        groups.setdefault(spec.group, []).append(spec)
    primary_group = next(
        spec.group for spec in specs if spec.primary
    )

    results: list[MemberResult] = []
    core_identity: dict[str, Any] = {}
    core_reference: dict[str, Any] = {}
    initial_data: dict[str, Any] = {}
    core_rows_by_group: dict[str, int] = {}

    for label, members in groups.items():
        grids = [build_member_grid(config, spec) for spec in members]
        rows = core_index_count(grids[0], core_radius=core_radius)
        core_rows_by_group[label] = rows
        core_identity[label] = core_identity_report(grids, members, core_rows=rows)
        initials = [
            measure_initial_state(
                grid,
                amplitude_scale=amplitude_scale,
                rho1=rho1,
                rho2=rho2,
                core_rows=rows,
            )
            for grid in grids
        ]
        core_reference[label] = core_reference_report(initials, members)
        initial_data[label] = initial_data_acceptance(
            unit_reference_grid(config, members[0]),
            amplitude_scale=amplitude_scale,
            rho1=rho1,
            rho2=rho2,
            core_radius=core_radius,
        )
        reference = float(core_reference[label]["reference"])
        for spec in members:
            results.append(
                evolve_member(
                    config,
                    spec,
                    reference_vorticity=reference,
                    checkpoint_dir=checkpoint_dir,
                    provenance=provenance,
                )
            )

    metrics = [member_metrics(item, cfl_coefficient=coefficient) for item in results]
    usable = [item for item in metrics if "all_diagnostics_finite" in item]

    elliptic: list[dict[str, Any]] = []
    for label, members in groups.items():
        group_results = [item for item in results if item.spec.group == label]
        elliptic.extend(
            elliptic_nonlocal_report(
                group_results,
                core_rows=core_rows_by_group[label],
            )
        )

    separations_by_group: dict[str, list[dict[str, Any]]] = {}
    amplifications_by_group: dict[str, dict[float, float]] = {}
    displacement_by_group: dict[str, dict[str, Any] | None] = {}
    for label, members in groups.items():
        ordered = sorted(
            (item for item in results if item.spec.group == label),
            key=lambda item: item.spec.wall_radius,
        )
        radii = [item.spec.wall_radius for item in ordered]
        finals = [
            item.snapshots[-1].core_amplification if item.snapshots else float("nan")
            for item in ordered
        ]
        amplifications_by_group[label] = {
            radius: value for radius, value in zip(radii, finals)
        }
        separations_by_group[label] = amplification_separations(finals, radii)
        if len(ordered) >= 2 and ordered[-1].snapshots and ordered[-2].snapshots:
            displacement_by_group[label] = argmax_displacement(
                (
                    ordered[-2].snapshots[-1].core_argmax_u1_r,
                    ordered[-2].snapshots[-1].core_argmax_u1_z,
                ),
                (
                    ordered[-1].snapshots[-1].core_argmax_u1_r,
                    ordered[-1].snapshots[-1].core_argmax_u1_z,
                ),
                z_period=float(config["z_period"]),
            )
        else:
            displacement_by_group[label] = None

    consistency = resolution_consistency_report(
        separations_by_group,
        amplifications_by_group,
        primary_group=primary_group,
        tolerance=float(
            config["decision_rule"]["resolution_consistency_relative_tolerance"]
        ),
    )
    decision = classify_wall_dependence(
        separations_by_group[primary_group],
        displacement_by_group[primary_group],
        thresholds=config["decision_rule"],
        resolution_stable=bool(consistency["stable_under_refinement"]),
    )

    checks = {
        "all_members_completed": bool(metrics)
        and all(item["completed"] for item in metrics),
        "all_diagnostics_finite": bool(usable)
        and all(item["all_diagnostics_finite"] for item in usable),
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
        "odd_symmetry_preserved": bool(usable)
        and all(
            item["maximum_odd_symmetry_defect_ratio"]
            <= float(acceptance["maximum_odd_symmetry_defect_ratio"])
            for item in usable
        ),
        "advective_cfl_within_tolerance": bool(usable)
        and all(
            advective_cfl_within_tolerance(
                item["maximum_advective_cfl"],
                cfl_coefficient=coefficient,
                cfl_excess_tolerance=tolerance,
            )
            for item in usable
        ),
        # ---- the core-identity requirements of the preregistered design -----
        "core_discretization_identical": bool(core_identity)
        and all(
            item["radial_spacing_identical"]
            and item["core_radial_nodes_bitwise_identical"]
            and item["axial_grid_identical"]
            for item in core_identity.values()
        ),
        "initial_core_reference_shared": bool(core_reference)
        and all(
            item["relative_spread"]
            <= float(acceptance["maximum_core_reference_relative_spread"])
            for item in core_reference.values()
        ),
        # ---- the preregistered E-32 initial-data checks (i)-(vi) ------------
        "initial_data_family_accepted": bool(initial_data)
        and all(
            item["core_bit_identical"]
            and item["sup_deviation_within_acceptance"]
            and item["sup_deviation_within_analytic_bound"]
            and item["exact_zero_outside_support"]
            and item["cutoff_minimum"] >= 0.0
            and item["cutoff_maximum"] <= 1.0
            and item["cutoff_monotone_non_increasing"]
            and item["derived_norms"]["max_abs_u1_relative_change"]
            <= PREREG_DERIVED_NORM_RELATIVE_TOLERANCE
            and item["derived_norms"]["max_cartesian_vorticity_relative_change"]
            <= PREREG_DERIVED_NORM_RELATIVE_TOLERANCE
            and item["fourth_difference"]["bounded"]
            and item["fourth_difference"]["finite"]
            and item["axis_parity"]["envelope_relative"]
            <= float(acceptance["maximum_initial_axis_parity_relative"])
            for item in initial_data.values()
        ),
        # The separate grid-convergence statement: the measured envelope norms
        # approach the derived E-29b constants at this resolution.
        "initial_norms_match_e29b": bool(metrics)
        and all(
            max(
                float(item["initial"]["e29b_max_abs_u1_relative_error"]),
                float(
                    item["initial"][
                        "e29b_max_cartesian_vorticity_relative_error"
                    ]
                ),
            )
            <= float(acceptance["maximum_initial_norm_relative_error"])
            for item in metrics
        ),
        "elliptic_nonlocal_contribution_recorded": bool(elliptic)
        and all(
            item.get("status") == "measured"
            and item.get("max_abs_difference") is not None
            and math.isfinite(float(item["max_abs_difference"]))
            and item.get("relative_denominator") is not None
            and math.isfinite(float(item["relative_denominator"]))
            for item in elliptic
        ),
        "cross_solver_elliptic_agreement_recorded": bool(usable)
        and all(
            item.get("snapshot_cross_solver_psi_max_abs_difference") is not None
            and math.isfinite(
                float(item["snapshot_cross_solver_psi_max_abs_difference"])
            )
            for item in usable
        ),
    }

    report = {
        "member_groups": {
            label: {
                "primary": label == primary_group,
                "members": [spec.as_dict() for spec in members],
                "core_identity": core_identity[label],
                "core_reference": core_reference[label],
                "initial_data_acceptance": initial_data[label],
                "final_core_amplification_by_radius": {
                    str(radius): value
                    for radius, value in amplifications_by_group[label].items()
                },
                "amplification_separations": separations_by_group[label],
                "argmax_displacement": displacement_by_group[label],
            }
            for label, members in groups.items()
        },
        "primary_group": primary_group,
        "elliptic_nonlocal_contribution": elliptic,
        "elliptic_nonlocal_definition": (
            "preregistration item 5: the larger-radius member's omega1 is "
            "restricted to Omega_c, embedded in BOTH members' domains, solved "
            "with the E-27 homogeneous outer Dirichlet condition, and the two "
            "psi1 fields are compared nodewise on Omega_c; the grids share the "
            "core nodes bitwise so no interpolation is involved"
        ),
        "resolution_consistency": consistency,
        "denominator_definitions": dict(DENOMINATOR_DEFINITIONS),
        "measured_quantities": (
            "preregistration section 3, all restricted to Omega_c = {r <= "
            f"{core_radius}}}: amplification against the shared "
            "||omega(0)||_inf, the argmax location of u1, max |psi1_z|, the "
            "E-20a core energy, the pairwise elliptic nonlocal contribution "
            "and the tail amplitudes beyond r = "
            f"{float(config['tail_radius'])}"
        ),
    }
    return results, metrics, {**report, "decision": decision}, checks


# -------------------------------------------------------------------- output


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


MEMBER_CSV_FIELDS: tuple[str, ...] = (
    "group",
    "wall_radius",
    "nr",
    "nz",
    "dr",
    "dz",
    "core_rows",
    "primary",
    "completed",
    "step_count",
    "final_time",
    "reference_vorticity",
    "initial_energy",
    "final_energy",
    "maximum_energy_growth_ratio",
    "maximum_circulation_growth_ratio",
    "final_core_amplification",
    "final_core_max_cartesian_vorticity",
    "final_core_argmax_u1_r",
    "final_core_argmax_u1_z",
    "final_core_max_abs_psi1_z",
    "final_core_energy",
    "final_tail_max_abs_u1",
    "final_tail_max_abs_omega1",
    "final_global_amplification",
    "maximum_odd_symmetry_defect_ratio",
    "maximum_divergence_residual_relative",
    "maximum_axis_parity_relative_u1",
    "maximum_axis_parity_relative_omega1",
    "maximum_wall_u1_abs",
    "maximum_advective_cfl",
    "maximum_viscous_cfl",
    "snapshot_cross_solver_psi_max_abs_difference",
    "snapshot_cross_solver_psi_relative_difference",
)

ELLIPTIC_CSV_FIELDS: tuple[str, ...] = (
    "group",
    "smaller_wall_radius",
    "larger_wall_radius",
    "snapshot_index",
    "time",
    "status",
    "max_abs_difference",
    "relative_difference",
    "relative_denominator",
    "small_core_max_abs_psi1",
    "large_core_max_abs_psi1",
    "argmax_r",
    "argmax_z",
    "source_core_max_abs_omega1",
)


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Run every member and write a checksummed evidence bundle."""

    validate_config(config)
    provenance = collect_runtime_provenance()
    config_bytes = (
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    _prepare_output(output_dir)
    (output_dir / "config.snapshot.json").write_bytes(config_bytes)
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir()
    results, metrics, report, checks = evaluate(
        config,
        checkpoint_dir=checkpoint_dir,
        provenance=provenance,
    )

    diagnostic_rows = [
        {
            "group": item.spec.group,
            "wall_radius": item.spec.wall_radius,
            "nr": item.spec.nr,
            "nz": item.spec.nz,
            **record,
        }
        for item in results
        for record in item.history
    ]
    if not diagnostic_rows:
        raise RuntimeError("no diagnostics were recorded")
    _write_csv(
        output_dir / "diagnostics.csv",
        ["group", "wall_radius", "nr", "nz", *DIAGNOSTIC_FIELDS],
        diagnostic_rows,
    )
    _write_csv(output_dir / "members.csv", list(MEMBER_CSV_FIELDS), metrics)
    snapshot_rows = [
        snapshot.as_row(item.spec) for item in results for snapshot in item.snapshots
    ]
    if snapshot_rows:
        _write_csv(
            output_dir / "core_snapshots.csv",
            list(CORE_SNAPSHOT_FIELDS),
            snapshot_rows,
        )
    if report["elliptic_nonlocal_contribution"]:
        _write_csv(
            output_dir / "elliptic_nonlocal.csv",
            list(ELLIPTIC_CSV_FIELDS),
            report["elliptic_nonlocal_contribution"],
        )

    arrays: dict[str, np.ndarray] = {}
    for item in results:
        suffix = f"{item.spec.group}_{radius_slug(item.spec.wall_radius)}"
        arrays[f"r_{suffix}"] = item.grid.r
        arrays[f"z_{suffix}"] = item.grid.z
        if item.state is not None:
            arrays[f"u1_final_{suffix}"] = item.state.u1
            arrays[f"omega1_final_{suffix}"] = item.state.omega1
            arrays[f"psi1_final_{suffix}"] = item.state.psi1
        for name in ("time", "max_abs_u1", "max_cartesian_vorticity", "amplification"):
            arrays[f"{name}_{suffix}"] = np.asarray(
                _series(item.history, name) if item.history else [],
                dtype=np.float64,
            )
        arrays[f"core_amplification_{suffix}"] = np.asarray(
            [snapshot.core_amplification for snapshot in item.snapshots],
            dtype=np.float64,
        )
        arrays[f"snapshot_time_{suffix}"] = np.asarray(
            [snapshot.time for snapshot in item.snapshots], dtype=np.float64
        )
    np.savez_compressed(output_dir / "trajectories.npz", **arrays)

    decision = report["decision"]
    summary: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Preregistered wall-dependence observation of the E-32 "
            "compact-support initial-data family on nested finite cylinders "
            "sharing one core discretization. It is not a reproduction of any "
            "published amplification value, not a singularity candidate, and "
            "not a proof."
        ),
        "preregistration": "docs/wall_dependence_prereg.md",
        "method": {
            "equations": "E-11, E-12, E-13, E-14",
            "initial_data": (
                "E-32a: the E-29 datum times the E-32b smooth cutoff of r^2 "
                f"with rho1 = {float(config['envelope']['rho1'])!r} and rho2 = "
                f"{float(config['envelope']['rho2'])!r}; omega1(0) = 0"
            ),
            "initial_data_core_identity": (
                "the cutoff is the literal 1.0 for r <= 0.9, and multiplication "
                "by 1.0 is exact, so every wall radius starts from a "
                "bit-identical core datum"
            ),
            "time_integrator": "explicit Heun/RK2",
            "time_step_rule": (
                "adaptive min(C dr/max|u^r|, C dz/max|u^z|, "
                "C min(dr,dz)^2/(4 nu)) capped by max_time_step"
            ),
            "wall_conditions": (
                "E-27 applied at r = R_wall with the second-order Thom form "
                "E-31; the integrator applies both at the outermost radial row, "
                "so no solver change is needed for a larger wall radius"
            ),
            "poisson_solver": (
                "solver A: Fourier-z / r^3-flux finite volume with homogeneous "
                "outer Dirichlet data"
            ),
            "independent_cross_check": (
                "solver B: finite_cylinder_poisson non-divergence radial "
                "stencil, re-solved at every snapshot"
            ),
            "viscosity_protocol": (
                "E-30 stage one only: t_final = T_1 lies below the switch time "
                "t_0 = 0.00227375"
            ),
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
            # This module imports experiments.run_hou_early_time, so the
            # repository root has to be importable.  Running it as a module
            # from the repository root is the invocation that satisfies that
            # with PYTHONPATH=src alone.
            "command": (
                "PYTHONPATH=src python -m experiments.run_wall_dependence "
                "--config configs/wall_dependence.json --output-dir "
                "outputs/wall_dependence"
            ),
        },
        "derived_reference_norms": {
            "source": "E-29b (derived by this repository, absent from the paper)",
            "max_abs_u1": E29B_MAX_ABS_U1,
            "max_cartesian_vorticity": E29B_MAX_CARTESIAN_VORTICITY,
            "e32_sup_deviation_bound": E32_SUP_DEVIATION_BOUND,
        },
        "members": metrics,
        "analysis": {
            name: value for name, value in report.items() if name != "decision"
        },
        "wall_dependence_classification": decision["classification"],
        "wall_dependence_decision": decision,
        "acceptance_checks": checks,
        "accepted_as_wall_dependence_observation": bool(all(checks.values())),
        "known_gaps": list(KNOWN_GAPS),
        "limitations": list(LIMITATIONS),
        "limitations_preregistration_section_6_verbatim": list(
            PREREG_SECTION_6_VERBATIM
        ),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    artifact_names = [
        "config.snapshot.json",
        "diagnostics.csv",
        "members.csv",
        "summary.json",
        "trajectories.npz",
    ]
    if snapshot_rows:
        artifact_names.append("core_snapshots.csv")
    if report["elliptic_nonlocal_contribution"]:
        artifact_names.append("elliptic_nonlocal.csv")
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
        default=REPOSITORY_ROOT / "configs" / "wall_dependence.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "wall_dependence",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="wall dependence config",
    )
    if not isinstance(config, dict):
        parser.error("config must be a JSON object")
    summary = run(config, output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["accepted_as_wall_dependence_observation"] else 2


if __name__ == "__main__":
    sys.exit(main())
