"""Gate 6: mid-amplitude calibration, amplitude continuation, and certificates.

Gate 5 established the differentiated entry point at an amplitude where nothing
happened.  This experiment does four things:

1. **Calibration.**  Measure, as a function of box tightness and resolution, the
   core difference between outer boundary treatments against the Richardson
   estimate of the discretisation error.  This is a *pilot*: its numbers are not
   candidate evidence, they only select the verification setting.
2. **Verification.**  At the selected setting, separate every discretisation
   factor one at a time -- ``dr``, ``dz``, joint, ``dt``, two integrators,
   ``R_max``, ``Z_max`` -- and evaluate the preregistered Gate-6 criteria.
3. **Continuation.**  Sweep amplitude and shape over the fixed initial-data
   family, rank with a composite gate, and take the best three to higher
   resolution.  Candidates that fail promotion are rejected explicitly.
4. **Certificates.**  Demonstrate the domain-expansion trigger with its
   before/after invariants, and build and independently verify an
   interval-arithmetic certificate for one snapshot.

No result here is a proof about Navier--Stokes, and a candidate that survives
this gate would still be a numerical object.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from ns_certificate_lab.critical_l3 import analyze_critical_l3
from ns_certificate_lab.domain_expansion import (
    ExpansionEvent,
    ExpansionPolicy,
    embed_in_larger_box,
    expansion_is_required,
)
from ns_certificate_lab.free_space_recovery import (
    recover_free_space_velocity,
    source_moments,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.initial_data import SwirlFamily
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.snapshot_certificate import (
    build_snapshot_certificate,
    verify_snapshot_certificate,
)
from ns_certificate_lab.tail_propagation import (
    FieldBounds,
    PotentialErrorBounds,
    advection_defect_bound,
    gronwall_state_error,
)
from ns_certificate_lab.whole_space_evolution import (
    WholeSpaceEvolution,
    cartesian_audit_of_state,
)

EXPECTED_SCHEMA_VERSION = 1
EXPECTED_EXPERIMENT = "whole_space_gate6_calibration_and_continuation"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("config root must be an object")
    return data


def _validate_config(config: dict[str, Any]) -> None:
    required = {
        "schema_version", "experiment", "viscosity", "family", "calibration",
        "verification", "continuation", "domain_expansion", "snapshot_certificate",
    }
    if set(config) != required:
        raise ValueError("Gate 6 config has missing or unknown keys")
    if config["schema_version"] != EXPECTED_SCHEMA_VERSION:
        raise ValueError("unsupported Gate 6 config schema")
    if config["experiment"] != EXPECTED_EXPERIMENT:
        raise ValueError("experiment must equal the audited canonical value")
    if float(config["viscosity"]) <= 0.0:
        raise ValueError("a Clay candidate needs a fixed positive viscosity")


def _write_json(path: Path, data: object) -> None:
    text = (
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    )
    path.write_text(text, encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _order(coarse: float, fine: float) -> float:
    return math.log(coarse / fine) / math.log(2.0)


def _evolve(
    family: SwirlFamily,
    *,
    viscosity: float,
    nr: int,
    nz: int,
    r_max: float,
    z_max: float,
    steps: int,
    dt: float,
    mode: str = "dipole",
    integrator: str = "rk4",
):
    grid = AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=r_max, z_min=-z_max, z_max=z_max, periodic_z=False
    )
    evolution = WholeSpaceEvolution(grid, viscosity=viscosity, boundary_mode=mode)
    u1, omega1 = family.initial_state(grid)
    for _ in range(steps):
        u1, omega1 = evolution.step(u1, omega1, dt, integrator)
        if not (np.all(np.isfinite(u1)) and np.all(np.isfinite(omega1))):
            raise FloatingPointError("state left the finite range")
    return grid, evolution, u1, omega1


def _core_mask(grid: AxisymmetricGrid, fraction: float) -> np.ndarray:
    r_mesh, z_mesh = grid.mesh()
    radius = fraction * min(float(grid.r[-1]), float(grid.z[-1]))
    return np.hypot(r_mesh, z_mesh) <= radius


# --------------------------------------------------------------------------- #
# 1. calibration                                                               #
# --------------------------------------------------------------------------- #


def _calibration(config: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    block = config["calibration"]
    viscosity = float(config["viscosity"])
    family = SwirlFamily(
        amplitude=float(block["amplitude"]),
        radial_support=float(config["family"]["radial_support"]),
        axial_support=float(config["family"]["axial_support"]),
        axial_concentration=float(config["family"]["axial_concentration"]),
    )
    rows: list[dict[str, Any]] = []
    for r_max, z_max in block["boxes"]:
        for nr in block["radial_points"]:
            nr = int(nr)
            nz = 2 * nr - 1
            states: dict[str, np.ndarray] = {}
            grid = None
            for mode in block["boundary_modes"]:
                grid, _evolution, _u1, omega1 = _evolve(
                    family, viscosity=viscosity, nr=nr, nz=nz,
                    r_max=float(r_max), z_max=float(z_max),
                    steps=int(block["steps"]), dt=float(block["time_step"]), mode=mode,
                )
                states[mode] = omega1
            mask = _core_mask(grid, float(block["core_fraction"]))
            scale = max(float(np.max(np.abs(states["dipole"]))), 1.0e-300)

            def difference(a: str, b: str) -> float:
                return float(np.max(np.abs((states[a] - states[b]) * mask))) / scale

            fine_grid, _e, _u, fine = _evolve(
                family, viscosity=viscosity, nr=2 * nr - 1, nz=2 * nz - 1,
                r_max=float(r_max), z_max=float(z_max),
                steps=int(block["steps"]), dt=float(block["time_step"]), mode="dipole",
            )
            richardson = float(
                np.max(np.abs((states["dipole"] - fine[::2, ::2]) * mask))
            ) / scale
            moments = source_moments(grid, states["dipole"])
            rows.append(
                {
                    "r_max": float(r_max),
                    "z_max": float(z_max),
                    "nr": nr,
                    "nz": nz,
                    "omega_scale": scale,
                    "zero_vs_monopole": difference("zero", "monopole"),
                    "monopole_vs_dipole": difference("monopole", "dipole"),
                    "dipole_vs_quadrupole": difference("dipole", "quadrupole"),
                    "zero_vs_dipole": difference("zero", "dipole"),
                    "richardson": richardson,
                    "boundary_over_richardson": difference("zero", "dipole")
                    / max(richardson, 1.0e-300),
                    "relative_monopole": abs(moments.monopole) / scale,
                    "relative_axial_dipole": abs(moments.axial_dipole) / scale,
                    "relative_axial_quadrupole": abs(moments.axial_quadrupole) / scale,
                }
            )
    best = max(rows, key=lambda row: row["boundary_over_richardson"])
    # The degenerate traces are compared as *ratios* against the one trace step
    # that is not degenerate.  An absolute threshold would be arbitrary; the
    # ratio is what "these two traces are the same object" means numerically.
    live = max(row["monopole_vs_dipole"] for row in rows)
    degenerate = max(
        max(row["zero_vs_monopole"] for row in rows),
        max(row["dipole_vs_quadrupole"] for row in rows),
    )
    resolution_groups: list[float] = []
    for r_max, z_max in {(row["r_max"], row["z_max"]) for row in rows}:
        group = [
            row["zero_vs_dipole"]
            for row in rows
            if row["r_max"] == r_max and row["z_max"] == z_max
        ]
        resolution_groups.append(
            (max(group) - min(group)) / max(max(group), 1.0e-300)
        )
    finding = {
        "relative_monopole_max": max(row["relative_monopole"] for row in rows),
        "monopole_vanishes_by_compact_support": max(
            row["relative_monopole"] for row in rows
        )
        < 1.0e-12,
        "degenerate_over_live_trace_ratio": degenerate / max(live, 1.0e-300),
        "monopole_trace_equals_zero_trace": degenerate / max(live, 1.0e-300) < 1.0e-9,
        "quadrupole_trace_equals_dipole_trace": (
            max(row["dipole_vs_quadrupole"] for row in rows) / max(live, 1.0e-300)
            < 1.0e-9
        ),
        "boundary_difference_resolution_spread": max(resolution_groups),
        "boundary_difference_is_resolution_independent": max(resolution_groups) < 0.15,
        "best_boundary_over_richardson": best["boundary_over_richardson"],
        "best_setting": {
            "r_max": best["r_max"], "z_max": best["z_max"], "nr": best["nr"]
        },
    }
    return rows, finding


# --------------------------------------------------------------------------- #
# 2. verification: one factor at a time                                        #
# --------------------------------------------------------------------------- #


def _verification(config: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    block = config["verification"]
    viscosity = float(config["viscosity"])
    family = SwirlFamily(
        amplitude=float(block["amplitude"]),
        radial_support=float(config["family"]["radial_support"]),
        axial_support=float(config["family"]["axial_support"]),
        axial_concentration=float(config["family"]["axial_concentration"]),
    )
    base = dict(
        viscosity=viscosity, nr=int(block["nr"]), nz=int(block["nz"]),
        r_max=float(block["r_max"]), z_max=float(block["z_max"]),
        steps=int(block["steps"]), dt=float(block["time_step"]),
    )
    rows: list[dict[str, Any]] = []

    def record(label: str, **overrides) -> dict[str, Any]:
        settings = {**base, **overrides}
        grid, evolution, u1, omega1 = _evolve(family, **settings)
        recovery = recover_free_space_velocity(
            evolution.solver, omega1, support_threshold=1.0e-10
        )
        report = analyze_critical_l3(
            grid, recovery.u_r, grid.mesh()[0] * u1, recovery.u_z, shell_count=8
        )
        row = {
            "family": label,
            "nr": grid.nr, "nz": grid.nz,
            "r_max": float(grid.r[-1]), "z_max": float(grid.z[-1]),
            "time_step": settings["dt"], "steps": settings["steps"],
            "integrator": settings.get("integrator", "rk4"),
            "boundary_mode": settings.get("mode", "dipole"),
            "max_omega1": float(np.max(np.abs(omega1))),
            "max_u1": float(np.max(np.abs(u1))),
            "l3_norm": report.l3_norm,
            "radial_scale": report.radial_scale,
            "axial_scale": report.axial_scale,
            "effective_shell_count": report.effective_shell_count,
            "outer_radial_fraction": report.outer_radial_fraction,
            "poisson_residual": recovery.poisson_residual_max,
            "tail_bound_available": recovery.tail_bound_available,
            "gradient_tail_bound": (
                recovery.interior_tail.gradient if recovery.tail_bound_available else None
            ),
            "velocity_tail_bound": (
                recovery.velocity_tail_bound if recovery.tail_bound_available else None
            ),
        }
        rows.append(row)
        return row

    for nr in block["radial_points"]:
        record("radial", nr=int(nr))
    for nz in block["axial_points"]:
        record("axial", nz=int(nz))
    for nr in block["radial_points"]:
        record("joint", nr=int(nr), nz=2 * int(nr) - 1)
    for factor in block["time_factors"]:
        record(
            "time", steps=int(block["steps"]) * int(factor),
            dt=float(block["time_step"]) / int(factor),
        )
    for integrator in block["integrators"]:
        record("integrator", integrator=str(integrator))
    spacing = float(block["r_max"]) / (int(block["nr"]) - 1)
    for r_max in block["radial_extents"]:
        record(
            "radial_domain", r_max=float(r_max),
            nr=int(round(float(r_max) / spacing)) + 1,
        )
    for z_max in block["axial_extents"]:
        record(
            "axial_domain", z_max=float(z_max),
            nz=2 * int(round(float(z_max) / spacing)) + 1,
        )

    def series(name: str, key: str = "max_omega1") -> list[float]:
        return [row[key] for row in rows if row["family"] == name]

    joint = series("joint")
    time_series = series("time")
    integrators = series("integrator")
    radial_domain = series("radial_domain")
    axial_domain = series("axial_domain")
    signal = abs(joint[-1])
    summary = {
        "joint_differences": [
            abs(joint[i] - joint[-1]) / signal for i in range(len(joint) - 1)
        ],
        "radial_differences": [
            abs(v - series("radial")[-1]) / signal for v in series("radial")[:-1]
        ],
        "axial_differences": [
            abs(v - series("axial")[-1]) / signal for v in series("axial")[:-1]
        ],
        "time_differences": [
            abs(v - time_series[-1]) / signal for v in time_series[:-1]
        ],
        "integrator_spread": abs(integrators[0] - integrators[-1]) / signal,
        "radial_domain_spread": (max(radial_domain) - min(radial_domain)) / signal,
        "axial_domain_spread": (max(axial_domain) - min(axial_domain)) / signal,
    }
    return rows, summary


# --------------------------------------------------------------------------- #
# 3. amplitude and shape continuation                                          #
# --------------------------------------------------------------------------- #


def _composite_metrics(
    family: SwirlFamily, config: dict[str, Any], spec: dict[str, Any]
) -> dict[str, Any]:
    viscosity = float(config["viscosity"])
    grid = AxisymmetricGrid.uniform(
        nr=int(spec["nr"]), nz=int(spec["nz"]), r_max=float(spec["r_max"]),
        z_min=-float(spec["z_max"]), z_max=float(spec["z_max"]), periodic_z=False,
    )
    evolution = WholeSpaceEvolution(grid, viscosity=viscosity, boundary_mode="dipole")
    u1, omega1 = family.initial_state(grid)
    r_mesh, _z = grid.mesh()
    initial_recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    initial = analyze_critical_l3(
        grid, initial_recovery.u_r, r_mesh * u1, initial_recovery.u_z, shell_count=8
    )
    diverged = False
    try:
        for _ in range(int(spec["steps"])):
            u1, omega1 = evolution.step(u1, omega1, float(spec["time_step"]), "rk4")
            if not (np.all(np.isfinite(u1)) and np.all(np.isfinite(omega1))):
                diverged = True
                break
    except FloatingPointError:
        diverged = True
    if diverged:
        return {"diverged": True}
    recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    final = analyze_critical_l3(
        grid, recovery.u_r, r_mesh * u1, recovery.u_z, shell_count=8
    )
    return {
        "diverged": False,
        "l3_growth": final.l3_norm / max(initial.l3_norm, 1.0e-300),
        "l3_norm": final.l3_norm,
        "max_omega1": float(np.max(np.abs(omega1))),
        "max_velocity": final.max_velocity,
        "radial_scale_ratio": final.radial_scale / max(initial.radial_scale, 1.0e-300),
        "axial_scale_ratio": final.axial_scale / max(initial.axial_scale, 1.0e-300),
        "effective_shell_count": final.effective_shell_count,
        "shell_count_change": final.effective_shell_count - initial.effective_shell_count,
        "critical_scaling_product": final.critical_scaling_product,
        "outer_radial_fraction": final.outer_radial_fraction,
        "shell_fractions": list(final.shell_fractions),
        "poisson_residual": recovery.poisson_residual_max,
        "tail_bound_available": recovery.tail_bound_available,
        "gradient_tail_bound": (
            recovery.interior_tail.gradient if recovery.tail_bound_available else None
        ),
        "neglected_source_fraction": recovery.moments.neglected_fraction,
    }


def _continuation(config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    block = config["continuation"]
    pilot = block["pilot"]
    rows: list[dict[str, Any]] = []
    for amplitude in block["amplitudes"]:
        for radial in block["radial_supports"]:
            for axial in block["axial_supports"]:
                for concentration in block["concentrations"]:
                    family = SwirlFamily(
                        amplitude=float(amplitude),
                        radial_support=float(radial),
                        axial_support=float(axial),
                        axial_concentration=float(concentration),
                    )
                    metrics = _composite_metrics(family, config, pilot)
                    entry = {**family.as_dict(), **metrics, "label": family.label()}
                    entry["composite_score"] = (
                        0.0
                        if metrics.get("diverged")
                        else (
                            max(metrics["l3_growth"] - 1.0, 0.0)
                            + max(1.0 - metrics["radial_scale_ratio"], 0.0)
                            + max(1.0 - metrics["axial_scale_ratio"], 0.0)
                            + max(metrics["shell_count_change"], 0.0)
                        )
                    )
                    rows.append(entry)
    ranked = sorted(rows, key=lambda row: row["composite_score"], reverse=True)
    promoted: list[dict[str, Any]] = []
    for entry in ranked[: int(block["top_candidates"])]:
        family = SwirlFamily(
            amplitude=entry["amplitude"], radial_support=entry["radial_support"],
            axial_support=entry["axial_support"],
            axial_concentration=entry["axial_concentration"],
        )
        series = []
        for nr in block["promotion_resolutions"]:
            spec = {**pilot, "nr": int(nr), "nz": 2 * int(nr) - 1}
            series.append({"nr": int(nr), **_composite_metrics(family, config, spec)})
        l3 = [row["l3_growth"] for row in series if not row.get("diverged")]
        widths = [row["radial_scale_ratio"] for row in series if not row.get("diverged")]
        promotion = {
            "same_trend_at_three_resolutions": len(l3) == 3
            and (all(v > 1.0 for v in l3) or all(v <= 1.0 for v in l3)),
            "global_l3_increases": len(l3) == 3 and min(l3) > 1.0,
            "critical_width_shrinks": len(widths) == 3 and max(widths) < 1.0,
            "shell_count_increases": all(
                row.get("shell_count_change", -1.0) > 0.0 for row in series
            ),
            "tail_bound_below_signal": all(
                row.get("tail_bound_available", False)
                and row["gradient_tail_bound"]
                < 0.1 * abs(row.get("max_omega1", 0.0)) + 1.0e-300
                for row in series
            ),
        }
        promotion["promoted"] = all(promotion.values())
        promoted.append(
            {
                "label": entry["label"],
                **family.as_dict(),
                "resolution_series": series,
                "promotion": promotion,
                "verdict": "promoted" if promotion["promoted"] else "rejected",
                "rejection_reasons": [
                    key for key, value in promotion.items()
                    if key != "promoted" and not value
                ],
            }
        )
    return ranked, promoted


# --------------------------------------------------------------------------- #
# 4. domain expansion and the snapshot certificate                             #
# --------------------------------------------------------------------------- #


def _domain_expansion(config: dict[str, Any]) -> dict[str, Any]:
    block = config["domain_expansion"]
    policy = ExpansionPolicy(**{k: v for k, v in block["policy"].items()})
    probe = block["trigger_probe"]
    fire, reasons = expansion_is_required(policy, **probe)
    viscosity = float(config["viscosity"])
    family = SwirlFamily(
        amplitude=float(config["calibration"]["amplitude"]),
        radial_support=float(config["family"]["radial_support"]),
        axial_support=float(config["family"]["axial_support"]),
        axial_concentration=float(config["family"]["axial_concentration"]),
    )
    grid, evolution, u1, omega1 = _evolve(
        family, viscosity=viscosity, nr=49, nz=97, r_max=2.0, z_max=2.4,
        steps=40, dt=4.0e-4,
    )
    r_mesh, _ = grid.mesh()
    before_recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    before = analyze_critical_l3(
        grid, before_recovery.u_r, r_mesh * u1, before_recovery.u_z, shell_count=8
    )
    energy_before = evolution._kinetic_energy(u1, before_recovery.u_r, before_recovery.u_z)
    divergence_before, _curl = evolution._consistency_checks(
        omega1, before_recovery.dpsi1_dz, before_recovery.u_r, before_recovery.u_z
    )

    new_grid, (new_u1, new_omega1), defect = embed_in_larger_box(
        grid, (u1, omega1), policy
    )
    new_evolution = WholeSpaceEvolution(
        new_grid, viscosity=viscosity, boundary_mode="dipole"
    )
    new_r, _ = new_grid.mesh()
    after_recovery = recover_free_space_velocity(
        new_evolution.solver, new_omega1, support_threshold=1.0e-10
    )
    after = analyze_critical_l3(
        new_grid, after_recovery.u_r, new_r * new_u1, after_recovery.u_z, shell_count=8
    )
    energy_after = new_evolution._kinetic_energy(
        new_u1, after_recovery.u_r, after_recovery.u_z
    )
    divergence_after, _curl2 = new_evolution._consistency_checks(
        new_omega1, after_recovery.dpsi1_dz, after_recovery.u_r, after_recovery.u_z
    )
    event = ExpansionEvent(
        step=40, time=40 * 4.0e-4, reasons=reasons,
        old_r_max=float(grid.r[-1]), old_z_max=float(grid.z[-1]),
        new_r_max=float(new_grid.r[-1]), new_z_max=float(new_grid.z[-1]),
        energy_before=energy_before, energy_after=energy_after,
        l3_before=before.l3_norm, l3_after=after.l3_norm,
        omega_max_before=float(np.max(np.abs(omega1))),
        omega_max_after=float(np.max(np.abs(new_omega1))),
        shell_mass_before=before.shell_masses, shell_mass_after=after.shell_masses,
        divergence_before=divergence_before, divergence_after=divergence_after,
        interpolation_defect=defect,
    )
    return {
        "policy": policy.as_dict(),
        "probe": probe,
        "trigger_fired": fire,
        "trigger_reasons": list(reasons),
        "event": event.as_dict(),
    }


def _snapshot_certificate(config: dict[str, Any]) -> dict[str, Any]:
    block = config["snapshot_certificate"]
    viscosity = float(config["viscosity"])
    family = SwirlFamily(
        amplitude=float(config["calibration"]["amplitude"]),
        radial_support=float(config["family"]["radial_support"]),
        axial_support=float(config["family"]["axial_support"]),
        axial_concentration=float(config["family"]["axial_concentration"]),
    )
    grid, evolution, u1, omega1 = _evolve(
        family, viscosity=viscosity, nr=int(block["nr"]), nz=int(block["nz"]),
        r_max=float(block["r_max"]), z_max=float(block["z_max"]),
        steps=int(block["steps"]), dt=float(block["time_step"]),
    )
    recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    certificate = build_snapshot_certificate(
        grid, u1, omega1, recovery.psi1, viscosity=viscosity,
        green_tail_bound=recovery.interior_tail.value,
        derivative_tail_bound=recovery.interior_tail.gradient,
        interior_radius=float(block["interior_radius"]),
        precision_bits=int(block["precision_bits"]),
    )
    payload = json.loads(json.dumps(certificate.as_dict()))
    checks = verify_snapshot_certificate(payload)
    return {"certificate": payload, "independent_checks": checks}


def _tail_propagation(config: dict[str, Any]) -> dict[str, Any]:
    viscosity = float(config["viscosity"])
    family = SwirlFamily(
        amplitude=float(config["calibration"]["amplitude"]),
        radial_support=float(config["family"]["radial_support"]),
        axial_support=float(config["family"]["axial_support"]),
        axial_concentration=float(config["family"]["axial_concentration"]),
    )
    grid, evolution, u1, omega1 = _evolve(
        family, viscosity=viscosity, nr=73, nz=145, r_max=3.0, z_max=3.6,
        steps=40, dt=4.0e-4,
    )
    from ns_certificate_lab.operators import derivative_r, derivative_z

    recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    if not recovery.tail_bound_available:
        raise RuntimeError(
            "the tail-propagation demonstration needs a box the source stays "
            "inside; widen it rather than fabricating a bound"
        )
    potential = PotentialErrorBounds(
        value=recovery.interior_tail.value,
        gradient=recovery.interior_tail.gradient,
        hessian=recovery.interior_tail.hessian,
        radial_extent=float(grid.r[-1]),
    )
    fields = FieldBounds(
        u1_max=float(np.max(np.abs(u1))),
        omega1_max=float(np.max(np.abs(omega1))),
        u1_gradient_max=max(
            float(np.max(np.abs(derivative_r(grid, u1, even_at_axis=True)))),
            float(np.max(np.abs(derivative_z(grid, u1)))),
        ),
        omega1_gradient_max=max(
            float(np.max(np.abs(derivative_r(grid, omega1, even_at_axis=True)))),
            float(np.max(np.abs(derivative_z(grid, omega1)))),
        ),
        velocity_radial_max=float(np.max(np.abs(recovery.u_r))),
        velocity_axial_max=float(np.max(np.abs(recovery.u_z))),
        psi_axial_derivative_max=float(np.max(np.abs(recovery.dpsi1_dz))),
    )
    defect = advection_defect_bound(potential, fields)
    horizon = 40 * 4.0e-4
    gronwall = gronwall_state_error(
        initial_error=0.0, constant_defect=defect.total,
        lipschitz=defect.state_lipschitz, horizon=horizon,
    )
    return {
        "potential": potential.as_dict(),
        "fields": fields.as_dict(),
        "defect": defect.as_dict(),
        "gronwall": gronwall.as_dict(),
        "signal": float(np.max(np.abs(omega1))),
        "bound_below_signal": gronwall.bound < float(np.max(np.abs(omega1))),
    }


# --------------------------------------------------------------------------- #


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = _load_json(config_path)
    _validate_config(config)
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {output_dir}")
    provenance = collect_runtime_provenance()
    output_dir.mkdir(parents=True)

    calibration_rows, calibration_finding = _calibration(config)
    verification_rows, verification_summary = _verification(config)
    continuation_rows, promoted = _continuation(config)
    expansion = _domain_expansion(config)
    certificate = _snapshot_certificate(config)
    propagation = _tail_propagation(config)

    thresholds = config["calibration"]["thresholds"]
    # The preregistered Gate-6 criterion below is reported exactly as written and
    # was NOT retuned after the data were seen.  The degeneracy statements are
    # findings, not pass criteria, and live in `calibration_finding`.
    gate6 = {
        "boundary_difference_exceeds_richardson": (
            calibration_finding["best_boundary_over_richardson"]
            >= float(thresholds["boundary_over_richardson"])
        ),
        "spatial_refinement_converges": (
            verification_summary["joint_differences"][0]
            > verification_summary["joint_differences"][-1]
        ),
        "time_refinement_converges": (
            verification_summary["time_differences"][0]
            >= verification_summary["time_differences"][-1]
        ),
        "integrators_agree_within_1e_3": (
            verification_summary["integrator_spread"] < 1.0e-3
        ),
        "radial_domain_independent_within_1e_3": (
            verification_summary["radial_domain_spread"] < 1.0e-3
        ),
        "axial_domain_independent_within_1e_3": (
            verification_summary["axial_domain_spread"] < 1.0e-3
        ),
        "expansion_trigger_fires": bool(expansion["trigger_fired"]),
        "expansion_preserves_invariants": bool(
            expansion["event"]["invariants_preserved"]
        ),
        "snapshot_certificate_verifies": bool(
            certificate["independent_checks"]["all_passed"]
        ),
        "tail_propagation_bound_below_signal": bool(propagation["bound_below_signal"]),
        "every_candidate_resolved": all(
            entry["verdict"] in {"promoted", "rejected"} for entry in promoted
        ),
    }
    # Quadratic-response probe: if max|omega1| scales exactly like A^2 across the
    # continuation, the run never left the weakly nonlinear regime and the
    # ranking cannot see nonlinear behaviour at all.  That is a property of the
    # run, so it is measured rather than assumed.
    quadratic = [
        row for row in continuation_rows
        if row["radial_support"] == 1.2 and row["axial_support"] == 1.5
        and row["axial_concentration"] == 2.0 and not row.get("diverged")
    ]
    quadratic.sort(key=lambda row: row["amplitude"])
    residuals = [
        abs(
            quadratic[i + 1]["max_omega1"] / quadratic[i]["max_omega1"]
            - (quadratic[i + 1]["amplitude"] / quadratic[i]["amplitude"]) ** 2
        )
        / (quadratic[i + 1]["amplitude"] / quadratic[i]["amplitude"]) ** 2
        for i in range(len(quadratic) - 1)
    ] if len(quadratic) > 1 else []
    gate6["continuation_left_the_quadratic_regime"] = bool(
        residuals and max(residuals) > 1.0e-2
    )
    gate6["all_passed"] = all(gate6.values())

    summary = {
        "schema_version": 1,
        "experiment": EXPECTED_EXPERIMENT,
        "statement": (
            "mid-amplitude calibration, one-factor verification, amplitude and "
            "shape continuation, domain expansion and a single-snapshot interval "
            "certificate; no proof about Navier-Stokes is claimed"
        ),
        "viscosity": float(config["viscosity"]),
        "family": config["family"],
        "calibration": calibration_rows,
        "calibration_finding": calibration_finding,
        "verification": verification_rows,
        "verification_summary": verification_summary,
        "continuation": continuation_rows,
        "promoted": promoted,
        "domain_expansion": expansion,
        "snapshot_certificate": certificate,
        "tail_propagation": propagation,
        "quadratic_response_residuals": residuals,
        "gate6": gate6,
        "limitations": [
            "the calibration criterion 'core boundary difference >= 8x Richardson' "
            "is NOT met and the threshold was not retuned after seeing the data; "
            "the reason is recorded in calibration_finding",
            "for compactly supported data the five-dimensional monopole and the "
            "axial quadrupole of omega1 vanish by symmetry, so the zero and "
            "monopole traces coincide, as do the dipole and quadrupole ones",
            "binary64 everywhere except the snapshot certificate",
            "the snapshot certificate encloses discrete quantities only; the "
            "discretisation error itself is not enclosed, so PO-05 stays open",
            "the tail-propagation Gronwall bound carries an assumed solve "
            "operator norm of zero, so it bounds the tail contribution only",
            "the continuation is a short-time fixed-viscosity sweep; it cannot "
            "see a singularity even if one existed",
            "the continuation never left the quadratic-response regime: "
            "max|omega1| scales as A^2 to better than one percent across a "
            "factor of ten in amplitude, so the ranking is blind to nonlinear "
            "behaviour and no candidate could have been promoted on merit",
            "at mid amplitude the time refinement and the integrator comparison "
            "are still uninformative: both agree to roundoff, so dt is not a "
            "limiting error for this problem class",
        ],
        "reproducibility": {"runtime_provenance": provenance},
    }

    written: list[Path] = []
    for name, payload in (("summary.json", summary), ("config.snapshot.json", config)):
        path = output_dir / name
        _write_json(path, payload)
        written.append(path)

    def write_csv(name: str, rows: list[dict[str, Any]], fields: list[str]) -> None:
        path = output_dir / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key) for key in fields})
        written.append(path)

    write_csv(
        "calibration.csv", calibration_rows,
        ["r_max", "z_max", "nr", "nz", "omega_scale", "zero_vs_monopole",
         "monopole_vs_dipole", "dipole_vs_quadrupole", "zero_vs_dipole",
         "richardson", "boundary_over_richardson", "relative_monopole",
         "relative_axial_dipole", "relative_axial_quadrupole"],
    )
    write_csv(
        "verification.csv", verification_rows,
        ["family", "nr", "nz", "r_max", "z_max", "time_step", "steps", "integrator",
         "boundary_mode", "max_omega1", "max_u1", "l3_norm", "radial_scale",
         "axial_scale", "effective_shell_count", "outer_radial_fraction",
         "poisson_residual", "gradient_tail_bound", "velocity_tail_bound"],
    )
    write_csv(
        "continuation.csv", continuation_rows,
        ["label", "amplitude", "radial_support", "axial_support",
         "axial_concentration", "aspect_ratio", "composite_score", "l3_growth",
         "max_omega1", "radial_scale_ratio", "axial_scale_ratio",
         "shell_count_change", "outer_radial_fraction", "diverged"],
    )

    manifest = {
        "schema_version": 1,
        "payloads": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in written
        },
    }
    manifest_path = output_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    (output_dir / "manifest.json.sha256").write_text(
        _sha256(manifest_path) + "\n", encoding="ascii"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.config, args.output_dir)
    print(json.dumps(summary["gate6"], ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
