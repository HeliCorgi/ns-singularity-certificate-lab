"""Gate 5: derivative tail bounds, velocity recovery, and a small-amplitude run.

Gate 4 certified the free-space potential **value**.  This experiment promotes
that to a whole-space nonlinear entry point by measuring, in order:

1. the explicit **derivative** tail bounds of
   ``ns_certificate_lab.free_space_recovery`` against a manufactured source with
   an exact closed-form solution, at every derivative order they cover;
2. the **analytic-Green** evaluation path, which differentiates the kernel and
   quadratures the source, against the same exact solution;
3. the **velocity recovery API** — space convergence, domain convergence, axis
   regularity, independent Cartesian divergence and curl, fault injection;
4. a **small-amplitude nonlinear whole-space run** from smooth compactly
   supported pure-swirl data, refined independently in ``dt``, ``dr``, ``dz``,
   ``R_max``, ``Z_max``, integrator and outer boundary order;
5. the **dyadic cascade toy model**, which decides whether low-mode-only forcing
   can drive high-shell amplitude through nonlinear transfer alone.

Passing does not license a Hou-amplitude whole-space run: it licenses the next
design step, which is resolution.  Nothing here is a proof of anything about
Navier--Stokes.
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

from ns_certificate_lab.cascade_toy import CascadeModel, decompose_forcing_roles
from ns_certificate_lab.free_space_recovery import (
    green_quadrature,
    multipole_boundary_trace,
    multipole_tail_bounds,
    recover_free_space_velocity,
    source_moments,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.whole_space_evolution import (
    SwirlBump,
    WholeSpaceEvolution,
    cartesian_audit_of_state,
    initial_state,
)
from ns_certificate_lab.whole_space_gate import (
    BumpSpec,
    ManufacturedFreeSpaceField,
    WholeSpaceEllipticSolver,
)

EXPECTED_SCHEMA_VERSION = 1
EXPECTED_EXPERIMENT = "whole_space_gate5_derivative_and_nonlinear"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("config root must be an object")
    return data


def _validate_config(config: dict[str, Any]) -> None:
    required = {
        "schema_version",
        "experiment",
        "source",
        "tail_bounds",
        "green_path",
        "recovery",
        "evolution",
        "cascade",
    }
    if set(config) != required:
        raise ValueError("Gate 5 config has missing or unknown keys")
    if config["schema_version"] != EXPECTED_SCHEMA_VERSION:
        raise ValueError("unsupported Gate 5 config schema")
    if config["experiment"] != EXPECTED_EXPERIMENT:
        raise ValueError("experiment must equal the audited canonical value")


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


def _order(coarse: float, fine: float, refinement: float = 2.0) -> float:
    return math.log(coarse / fine) / math.log(refinement)


def _field(block: dict[str, Any]) -> ManufacturedFreeSpaceField:
    return ManufacturedFreeSpaceField(
        bumps=tuple(
            BumpSpec(
                center=float(b["center"]),
                radius=float(b["radius"]),
                amplitude=float(b["amplitude"]),
                power=int(b["power"]),
            )
            for b in block["bumps"]
        )
    )


# --------------------------------------------------------------------------- #
# 1. derivative tail bounds                                                    #
# --------------------------------------------------------------------------- #


def _tail_bound_rows(config: dict[str, Any], field: ManufacturedFreeSpaceField):
    block = config["tail_bounds"]
    spec = block["moment_grid"]
    grid = AxisymmetricGrid.uniform(
        nr=int(spec["nr"]),
        nz=int(spec["nz"]),
        r_max=float(spec["r_max"]),
        z_min=-float(spec["z_max"]),
        z_max=float(spec["z_max"]),
        periodic_z=False,
    )
    r_mesh, z_mesh = grid.mesh()
    moments = source_moments(grid, field.omega1(r_mesh, z_mesh))
    step = float(block["finite_difference_step"])
    samples = int(block["samples_per_sphere"])
    rows: list[dict[str, Any]] = []
    for truncation in ("monopole", "dipole"):
        for radius in block["evaluation_radii"]:
            radius = float(radius)
            bound = multipole_tail_bounds(
                moments, evaluation_radius=radius, truncation=truncation
            )
            angle = np.linspace(0.0, math.pi, samples)
            rr = radius * np.sin(angle)
            zz = radius * np.cos(angle)

            def approx(a, c):
                return multipole_boundary_trace(grid, moments, a, c, truncation=truncation)

            value_error = float(np.max(np.abs(field.psi1(rr, zz) - approx(rr, zz))))
            grad_r = (approx(rr + step, zz) - approx(rr - step, zz)) / (2.0 * step)
            grad_z = (approx(rr, zz + step) - approx(rr, zz - step)) / (2.0 * step)
            gradient_error = float(
                np.max(
                    np.hypot(
                        field.dpsi1_dr(rr, zz) - grad_r, field.dpsi1_dz(rr, zz) - grad_z
                    )
                )
            )
            hess_rr = (approx(rr + step, zz) - 2.0 * approx(rr, zz) + approx(rr - step, zz)) / step**2
            hess_zz = (approx(rr, zz + step) - 2.0 * approx(rr, zz) + approx(rr, zz - step)) / step**2
            hessian_error = max(
                float(np.max(np.abs(field.d2psi1_drr(rr, zz) - hess_rr))),
                float(np.max(np.abs(field.d2psi1_dzz(rr, zz) - hess_zz))),
            )
            rows.append(
                {
                    "truncation": truncation,
                    "evaluation_radius": radius,
                    "distance": bound.distance,
                    "value_error": value_error,
                    "value_bound": bound.value,
                    "value_ratio": value_error / bound.value,
                    "gradient_error": gradient_error,
                    "gradient_bound": bound.gradient,
                    "gradient_ratio": gradient_error / bound.gradient,
                    "hessian_error": hessian_error,
                    "hessian_bound": bound.hessian,
                    "hessian_ratio": hessian_error / bound.hessian,
                    "third_bound": bound.third,
                }
            )
    understated = multipole_tail_bounds(
        source_moments(grid, field.omega1(r_mesh, z_mesh)),
        evaluation_radius=float(block["evaluation_radii"][0]),
        truncation="monopole",
    )
    factor = float(block["understated_moment_factor"])
    fault = {
        "understated_value_bound": understated.value / factor,
        "understated_gradient_bound": understated.gradient / factor,
        "measured_value_error": rows[0]["value_error"],
        "measured_gradient_error": rows[0]["gradient_error"],
        "value_bound_broken": rows[0]["value_error"] > understated.value / factor,
        "gradient_bound_broken": rows[0]["gradient_error"]
        > understated.gradient / factor,
    }
    return moments, rows, fault


# --------------------------------------------------------------------------- #
# 2. the analytic-Green path                                                   #
# --------------------------------------------------------------------------- #


def _green_rows(config: dict[str, Any], field: ManufacturedFreeSpaceField):
    block = config["green_path"]
    targets_r = np.asarray(block["targets_r"], dtype=np.float64)
    targets_z = np.asarray(block["targets_z"], dtype=np.float64)
    exact = {
        "psi": field.psi1(targets_r, targets_z),
        "dr": field.dpsi1_dr(targets_r, targets_z),
        "dz": field.dpsi1_dz(targets_r, targets_z),
        "drr": field.d2psi1_drr(targets_r, targets_z),
        "dzz": field.d2psi1_dzz(targets_r, targets_z),
        "drz": field.d2psi1_drz(targets_r, targets_z),
    }
    rows: list[dict[str, Any]] = []
    previous: dict[str, float] | None = None
    for spec in block["grids"]:
        grid = AxisymmetricGrid.uniform(
            nr=int(spec["nr"]),
            nz=int(spec["nz"]),
            r_max=float(block["r_max"]),
            z_min=-float(block["z_max"]),
            z_max=float(block["z_max"]),
            periodic_z=False,
        )
        r_mesh, z_mesh = grid.mesh()
        out = green_quadrature(
            grid,
            field.omega1(r_mesh, z_mesh),
            targets_r,
            targets_z,
            quadrature_nodes=int(block["quadrature_nodes"]),
        )
        row: dict[str, Any] = {"nr": grid.nr, "nz": grid.nz}
        for key, reference in exact.items():
            scale = max(float(np.max(np.abs(reference))), 1.0e-300)
            error = float(np.max(np.abs(out[key] - reference))) / scale
            row[f"relative_{key}"] = error
            if previous is not None:
                row[f"order_{key}"] = _order(previous[key], error)
        previous = {key: row[f"relative_{key}"] for key in exact}
        rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# 3. the velocity recovery API                                                 #
# --------------------------------------------------------------------------- #


def _recovery_rows(config: dict[str, Any], field: ManufacturedFreeSpaceField):
    block = config["recovery"]
    extent = float(block["extent"])
    interior_radius = float(block["interior_radius"])
    rows: list[dict[str, Any]] = []
    previous: dict[str, float] | None = None
    for nr in block["radial_points"]:
        nr = int(nr)
        grid = AxisymmetricGrid.uniform(
            nr=nr, nz=2 * nr - 1, r_max=extent, z_min=-extent, z_max=extent,
            periodic_z=False,
        )
        r_mesh, z_mesh = grid.mesh()
        solver = WholeSpaceEllipticSolver.build(grid, boundary_mode="monopole")
        recovery = recover_free_space_velocity(solver, field.omega1(r_mesh, z_mesh))
        mask = np.hypot(r_mesh, z_mesh) <= interior_radius
        exact_r = -r_mesh * field.dpsi1_dz(r_mesh, z_mesh)
        exact_z = 2.0 * field.psi1(r_mesh, z_mesh) + r_mesh * field.dpsi1_dr(
            r_mesh, z_mesh
        )
        errors = {
            "psi": float(np.max(np.abs((recovery.psi1 - field.psi1(r_mesh, z_mesh)) * mask))),
            "dr": float(np.max(np.abs((recovery.dpsi1_dr - field.dpsi1_dr(r_mesh, z_mesh)) * mask))),
            "dz": float(np.max(np.abs((recovery.dpsi1_dz - field.dpsi1_dz(r_mesh, z_mesh)) * mask))),
            "u_r": float(np.max(np.abs((recovery.u_r - exact_r) * mask))),
            "u_z": float(np.max(np.abs((recovery.u_z - exact_z) * mask))),
        }
        row: dict[str, Any] = {
            "nr": nr,
            "nz": grid.nz,
            "dr": grid.dr,
            "axis_u_r_max": float(np.max(np.abs(recovery.u_r[0]))),
            "axis_psi_r_max": float(np.max(np.abs(recovery.dpsi1_dr[0]))),
            "poisson_residual_max": recovery.poisson_residual_max,
            "interior_gradient_tail_bound": recovery.interior_tail.gradient,
            "velocity_tail_bound": recovery.velocity_tail_bound,
        }
        for key, value in errors.items():
            row[f"error_{key}"] = value
            if previous is not None:
                row[f"order_{key}"] = _order(previous[key], value)
        previous = errors
        rows.append(row)
    return rows


def _domain_rows(config: dict[str, Any], field: ManufacturedFreeSpaceField):
    block = config["recovery"]
    spacing = float(block["domain_spacing"])
    fixed = float(block["domain_fixed_extent"])
    interior_radius = float(block["interior_radius"])
    rows: list[dict[str, Any]] = []
    for enlarge in ("radial", "axial"):
        reference = None
        collected = []
        sizes = [float(v) for v in block["domain_sizes"]] + [4.0 * float(block["domain_sizes"][-1])]
        for size in sizes:
            r_max = size if enlarge == "radial" else fixed
            z_max = fixed if enlarge == "radial" else size
            grid = AxisymmetricGrid.uniform(
                nr=int(round(r_max / spacing)) + 1,
                nz=2 * int(round(z_max / spacing)) + 1,
                r_max=r_max, z_min=-z_max, z_max=z_max, periodic_z=False,
            )
            r_mesh, z_mesh = grid.mesh()
            solver = WholeSpaceEllipticSolver.build(grid, boundary_mode="monopole")
            recovery = recover_free_space_velocity(solver, field.omega1(r_mesh, z_mesh))
            mask = np.hypot(r_mesh, z_mesh) <= interior_radius
            exact_r = -r_mesh * field.dpsi1_dz(r_mesh, z_mesh)
            collected.append((size, r_max, z_max, grid, recovery, mask, exact_r))
        reference = collected[-1]
        for size, r_max, z_max, grid, recovery, mask, _exact in collected[:-1]:
            ref_grid = reference[3]
            ref = reference[4]
            axial_offset = int(round((grid.z[0] - ref_grid.z[0]) / grid.dz))
            axial_index = np.arange(grid.nz) + axial_offset
            keep = (axial_index >= 0) & (axial_index < ref_grid.nz)
            radial_index = np.arange(grid.nr)
            selected_r = ref.dpsi1_dr[np.ix_(radial_index, axial_index[keep])]
            selected_z = ref.dpsi1_dz[np.ix_(radial_index, axial_index[keep])]
            local_mask = mask[:, keep]
            rows.append(
                {
                    "enlarged": enlarge,
                    "size": size,
                    "r_max": r_max,
                    "z_max": z_max,
                    "gradient_common_diff": max(
                        float(np.max(np.abs((recovery.dpsi1_dr[:, keep] - selected_r) * local_mask))),
                        float(np.max(np.abs((recovery.dpsi1_dz[:, keep] - selected_z) * local_mask))),
                    ),
                    "gradient_tail_bound": recovery.interior_tail.gradient,
                }
            )
    return rows


def _recovery_faults(config: dict[str, Any], field: ManufacturedFreeSpaceField):
    block = config["recovery"]
    extent = float(block["extent"])
    grid = AxisymmetricGrid.uniform(
        nr=129, nz=257, r_max=extent, z_min=-extent, z_max=extent, periodic_z=False
    )
    r_mesh, z_mesh = grid.mesh()
    omega = field.omega1(r_mesh, z_mesh)
    mask = np.hypot(r_mesh, z_mesh) <= float(block["interior_radius"])
    healthy = recover_free_space_velocity(
        WholeSpaceEllipticSolver.build(grid, boundary_mode="monopole"), omega
    )
    faulty_axis = recover_free_space_velocity(
        WholeSpaceEllipticSolver.build(
            grid, boundary_mode="monopole", axis_coefficient=4.0
        ),
        omega,
    )
    faulty_drift = recover_free_space_velocity(
        WholeSpaceEllipticSolver.build(
            grid, boundary_mode="monopole", radial_drift_coefficient=1.0
        ),
        omega,
    )
    exact_dr = field.dpsi1_dr(r_mesh, z_mesh)

    def error(recovery) -> float:
        return float(np.max(np.abs((recovery.dpsi1_dr - exact_dr) * mask)))

    return {
        "healthy_gradient_error": error(healthy),
        "axis_coefficient_gradient_error": error(faulty_axis),
        "radial_drift_gradient_error": error(faulty_drift),
        "axis_fault_detected": error(faulty_axis) > 2.0 * error(healthy),
        "drift_fault_detected": error(faulty_drift) > 10.0 * error(healthy),
    }


def _cartesian_rows(config: dict[str, Any], field: ManufacturedFreeSpaceField):
    from ns_certificate_lab.whole_space_gate import audit_cartesian_recovery

    block = config["recovery"]
    rows = []
    previous = None
    for points in block["cartesian_points"]:
        report = audit_cartesian_recovery(
            field, half_width=float(block["cartesian_half_width"]), points=int(points)
        )
        row = report.as_dict()
        if previous is not None:
            row["divergence_order"] = _order(previous[0], report.divergence_max)
            row["curl_order"] = _order(previous[1], report.curl_theta_error_max)
        previous = (report.divergence_max, report.curl_theta_error_max)
        rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# 4. the small-amplitude nonlinear run                                         #
# --------------------------------------------------------------------------- #


def _evolution_rows(config: dict[str, Any]):
    block = config["evolution"]
    base = block["base"]
    bump = SwirlBump(
        center=float(block["swirl"]["center"]),
        radius=float(block["swirl"]["radius"]),
        amplitude=float(block["swirl"]["amplitude"]),
        power=int(block["swirl"]["power"]),
    )

    def run(
        *, nr, nz, r_max, z_max, steps, dt, integrator="rk4", boundary_mode="monopole"
    ):
        grid = AxisymmetricGrid.uniform(
            nr=nr, nz=nz, r_max=r_max, z_min=-z_max, z_max=z_max, periodic_z=False
        )
        evolution = WholeSpaceEvolution(
            grid, viscosity=float(block["viscosity"]), boundary_mode=boundary_mode
        )
        u1, omega1 = initial_state(grid, bump)
        return grid, evolution, evolution.run(
            u1, omega1, steps=steps, time_step=dt, integrator=integrator,
            snapshot_every=steps,
        )

    rows: list[dict[str, Any]] = []

    def record(label: str, summary, extra: dict[str, Any] | None = None) -> None:
        final = summary.final
        entry = {
            "family": label,
            "integrator": summary.integrator,
            "boundary_mode": summary.boundary_mode,
            "nr": summary.nr,
            "nz": summary.nz,
            "r_max": summary.r_max,
            "z_max": summary.z_max,
            "steps": summary.steps,
            "time_step": summary.time_step,
            "final_time": summary.final_time,
            "final_max_omega1": final.max_omega1,
            "final_energy": final.kinetic_energy,
            "final_l3": final.l3_norm,
            "effective_shell_count": final.effective_shell_count,
            "max_relative_divergence": summary.streaming_max_relative_divergence,
            "max_curl_defect": summary.streaming_max_curl_defect,
            "max_poisson_residual": summary.streaming_max_poisson_residual,
            "max_outer_band_fraction": summary.streaming_max_outer_band_fraction,
            "max_energy_increase": summary.streaming_max_energy_increase,
            "boundary_sensitivity": final.boundary_sensitivity,
            "derivative_tail_bound": final.derivative_tail_bound,
            "velocity_tail_bound": final.velocity_tail_bound,
            "low_frequency_fraction": final.low_frequency_fraction,
            "neglected_source_fraction": final.neglected_source_fraction,
            "shell_fractions": list(final.shell_fractions),
            "diverged": summary.diverged,
        }
        if extra:
            entry.update(extra)
        rows.append(entry)

    # time refinement
    for factor in block["time_refinement"]:
        factor = int(factor)
        _g, _e, summary = run(
            nr=int(base["nr"]), nz=int(base["nz"]), r_max=float(base["r_max"]),
            z_max=float(base["z_max"]), steps=int(base["steps"]) * factor,
            dt=float(base["time_step"]) / factor,
        )
        record("time_refinement", summary, {"refinement_factor": factor})
    # radial refinement at fixed nz
    for nr in block["radial_points"]:
        _g, _e, summary = run(
            nr=int(nr), nz=int(base["nz"]), r_max=float(base["r_max"]),
            z_max=float(base["z_max"]), steps=int(base["steps"]),
            dt=float(base["time_step"]),
        )
        record("radial_refinement", summary)
    # axial refinement at fixed nr
    for nz in block["axial_points"]:
        _g, _e, summary = run(
            nr=int(base["nr"]), nz=int(nz), r_max=float(base["r_max"]),
            z_max=float(base["z_max"]), steps=int(base["steps"]),
            dt=float(base["time_step"]),
        )
        record("axial_refinement", summary)
    # joint refinement, the honest second-order check
    for nr in block["radial_points"]:
        _g, _e, summary = run(
            nr=int(nr), nz=2 * int(nr) - 1, r_max=float(base["r_max"]),
            z_max=float(base["z_max"]), steps=int(base["steps"]),
            dt=float(base["time_step"]),
        )
        record("joint_refinement", summary)
    # domain enlargement at fixed spacing
    spacing = float(base["r_max"]) / (int(base["nr"]) - 1)
    for r_max in block["radial_extents"]:
        r_max = float(r_max)
        _g, _e, summary = run(
            nr=int(round(r_max / spacing)) + 1, nz=int(base["nz"]),
            r_max=r_max, z_max=float(base["z_max"]),
            steps=int(base["steps"]), dt=float(base["time_step"]),
        )
        record("radial_domain", summary)
    for z_max in block["axial_extents"]:
        z_max = float(z_max)
        _g, _e, summary = run(
            nr=int(base["nr"]), nz=2 * int(round(z_max / spacing)) + 1,
            r_max=float(base["r_max"]), z_max=z_max,
            steps=int(base["steps"]), dt=float(base["time_step"]),
        )
        record("axial_domain", summary)
    # integrators
    for integrator in block["integrators"]:
        _g, _e, summary = run(
            nr=int(base["nr"]), nz=int(base["nz"]), r_max=float(base["r_max"]),
            z_max=float(base["z_max"]), steps=int(base["steps"]),
            dt=float(base["time_step"]), integrator=str(integrator),
        )
        record("integrator", summary)
    # boundary order
    for mode in block["boundary_modes"]:
        _g, _e, summary = run(
            nr=int(base["nr"]), nz=int(base["nz"]), r_max=float(base["r_max"]),
            z_max=float(base["z_max"]), steps=int(base["steps"]),
            dt=float(base["time_step"]), boundary_mode=str(mode),
        )
        record("boundary_mode", summary)

    grid, evolution, summary = run(
        nr=int(base["nr"]), nz=int(base["nz"]), r_max=float(base["r_max"]),
        z_max=float(base["z_max"]), steps=int(base["steps"]),
        dt=float(base["time_step"]),
    )
    u1, omega1 = initial_state(grid, bump)
    for _ in range(int(base["steps"])):
        u1, omega1 = evolution.step(u1, omega1, float(base["time_step"]), "rk4")
    cartesian = cartesian_audit_of_state(
        grid, u1, omega1, evolution.solver.solve(omega1),
        half_width=float(block["cartesian_half_width"]),
        points=int(block["cartesian_points"]),
    )
    return rows, cartesian


# --------------------------------------------------------------------------- #
# 5. the cascade toy model                                                     #
# --------------------------------------------------------------------------- #


def _cascade_rows(config: dict[str, Any]):
    block = config["cascade"]
    shells = int(block["shells"])
    initial = np.zeros(shells, dtype=np.float64)
    initial[0] = float(block["seed_amplitude"])
    rows: list[dict[str, Any]] = []
    for case in block["cases"]:
        model = CascadeModel(
            shells=shells,
            viscosity=float(case["viscosity"]),
            forced_shells=int(case["forced_shells"]),
            forcing_amplitude=float(case["forcing"]),
            high_shell_start=int(block["high_shell_start"]),
        )
        run = model.run(
            initial, steps=int(block["steps"]), time_step=float(block["time_step"])
        )
        entry = run.as_dict()
        entry["label"] = case["label"]
        rows.append(entry)
    return rows


# --------------------------------------------------------------------------- #


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = _load_json(config_path)
    _validate_config(config)
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {output_dir}")
    provenance = collect_runtime_provenance()
    output_dir.mkdir(parents=True)

    field = _field(config["source"])
    moments, tail_rows, tail_fault = _tail_bound_rows(config, field)
    green_rows = _green_rows(config, field)
    recovery_rows = _recovery_rows(config, field)
    domain_rows = _domain_rows(config, field)
    recovery_faults = _recovery_faults(config, field)
    cartesian_rows = _cartesian_rows(config, field)
    evolution_rows, evolution_cartesian = _evolution_rows(config)
    cascade_rows = _cascade_rows(config)

    def family(name: str) -> list[dict[str, Any]]:
        return [row for row in evolution_rows if row["family"] == name]

    joint = family("joint_refinement")
    joint_orders = [
        _order(joint[i]["max_curl_defect"], joint[i + 1]["max_curl_defect"])
        for i in range(len(joint) - 1)
    ]
    time_rows = family("time_refinement")
    time_differences = [
        abs(time_rows[i]["final_max_omega1"] - time_rows[-1]["final_max_omega1"])
        for i in range(len(time_rows) - 1)
    ]
    integrator_values = [row["final_max_omega1"] for row in family("integrator")]
    integrator_spread = (
        (max(integrator_values) - min(integrator_values)) / max(abs(max(integrator_values)), 1.0e-300)
    )
    boundary_values = [row["final_max_omega1"] for row in family("boundary_mode")]
    boundary_spread = (
        abs(boundary_values[0] - boundary_values[-1]) / max(abs(boundary_values[0]), 1.0e-300)
    )
    radial_domain = [row["final_max_omega1"] for row in family("radial_domain")]
    axial_domain = [row["final_max_omega1"] for row in family("axial_domain")]

    forced = next(row for row in cascade_rows if row["label"] == "low_mode_forced_viscous")
    unforced = next(row for row in cascade_rows if row["label"] == "unforced_viscous")

    acceptance = {
        "tail_bounds_dominate_every_order": all(
            row["value_ratio"] <= 1.0
            and row["gradient_ratio"] <= 1.0
            and row["hessian_ratio"] <= 1.0
            for row in tail_rows
        ),
        "understated_moment_breaks_the_bound": (
            tail_fault["value_bound_broken"] and tail_fault["gradient_bound_broken"]
        ),
        "green_path_converges_on_every_derivative": all(
            green_rows[-1][f"order_{key}"] > 2.5
            for key in ("psi", "dr", "dz", "drr", "dzz", "drz")
        ),
        "green_path_relative_error_below_1e_7": all(
            green_rows[-1][f"relative_{key}"] < 1.0e-7
            for key in ("psi", "dr", "dz", "drr", "dzz", "drz")
        ),
        "recovery_second_order_in_every_field": all(
            recovery_rows[-1][f"order_{key}"] > 1.8
            for key in ("psi", "dr", "dz", "u_r", "u_z")
        ),
        "axis_regularity_exact": all(
            row["axis_u_r_max"] == 0.0 and row["axis_psi_r_max"] == 0.0
            for row in recovery_rows
        ),
        "gradient_domain_error_decreases": all(
            group[0]["gradient_common_diff"] > group[-1]["gradient_common_diff"]
            for group in (
                [r for r in domain_rows if r["enlarged"] == "radial"],
                [r for r in domain_rows if r["enlarged"] == "axial"],
            )
        ),
        "gradient_tail_bound_dominates_domain_error": all(
            row["gradient_common_diff"] <= row["gradient_tail_bound"] for row in domain_rows
        ),
        "recovery_axis_fault_detected": recovery_faults["axis_fault_detected"],
        "recovery_drift_fault_detected": recovery_faults["drift_fault_detected"],
        "cartesian_divergence_and_curl_second_order": (
            cartesian_rows[-1]["divergence_order"] > 1.8
            and cartesian_rows[-1]["curl_order"] > 1.8
        ),
        "evolution_never_diverged": all(not row["diverged"] for row in evolution_rows),
        "evolution_energy_never_increases": all(
            row["max_energy_increase"] <= 1.0e-12 for row in evolution_rows
        ),
        "evolution_joint_refinement_second_order": min(joint_orders) > 1.8,
        "evolution_time_refinement_converges": (
            time_differences[0] > time_differences[-1]
        ),
        "evolution_integrators_agree_within_1e_3": integrator_spread < 1.0e-3,
        "evolution_boundary_order_agrees_within_1e_6": boundary_spread < 1.0e-6,
        "evolution_radial_domain_converges": abs(
            radial_domain[-1] - radial_domain[-2]
        ) < abs(radial_domain[1] - radial_domain[0]) + 1.0e-30,
        "evolution_axial_domain_converges": abs(
            axial_domain[-1] - axial_domain[-2]
        ) < abs(axial_domain[1] - axial_domain[0]) + 1.0e-30,
        "evolution_outer_band_negligible": all(
            row["max_outer_band_fraction"] < 1.0e-12 for row in evolution_rows
        ),
        "evolution_cartesian_curl_consistent": evolution_cartesian["curl_relative"] < 0.1,
        "cascade_transfer_is_energy_conserving": all(
            row["nonlinear_transfer_defect"] < 1.0e-12 for row in cascade_rows
        ),
        "cascade_low_mode_forcing_drives_high_shells": (
            forced["direct_high_shell_injection"] == 0.0
            and forced["max_high_shell_amplitude"]
            > 1.0e6 * unforced["max_high_shell_amplitude"]
        ),
        "cascade_viscosity_competes_with_the_indirect_route": (
            next(
                row for row in cascade_rows
                if row["label"] == "low_mode_forced_high_viscosity"
            )["max_high_shell_amplitude"]
            < forced["max_high_shell_amplitude"]
        ),
    }
    acceptance["all_passed"] = all(acceptance.values())

    summary = {
        "schema_version": 1,
        "experiment": EXPECTED_EXPERIMENT,
        "statement": (
            "derivative tail bounds, free-space velocity recovery and a "
            "small-amplitude nonlinear whole-space run; a gate on the entry "
            "point, not a statement about any candidate"
        ),
        "source_moments": moments.as_dict(),
        "green_derivative_constants": "A_m with ||D^m G5|| = A_m |Z|^-3-m; see free_space_recovery",
        "tail_bounds": tail_rows,
        "tail_bound_fault_injection": tail_fault,
        "green_path": green_rows,
        "recovery": recovery_rows,
        "recovery_domain": domain_rows,
        "recovery_fault_injection": recovery_faults,
        "cartesian_audit": cartesian_rows,
        "evolution": evolution_rows,
        "evolution_cartesian_audit": evolution_cartesian,
        "evolution_joint_orders": joint_orders,
        "evolution_integrator_spread": integrator_spread,
        "evolution_boundary_spread": boundary_spread,
        "cascade": cascade_rows,
        "forcing_roles": [role.as_dict() for role in decompose_forcing_roles()],
        "acceptance": acceptance,
        "limitations": [
            "binary64 throughout; no interval enclosure anywhere",
            "the derivative tail bounds are continuum statements about "
            "psi_infinity and the continuum box solution; the discretisation "
            "error of the finite-difference solve is measured, not bounded",
            "the analytic-Green path is restricted to targets outside the "
            "source support, where the reduced kernel is smooth",
            "the nonlinear run is small amplitude and short time; it says "
            "nothing about a Hou-amplitude candidate",
            "u1 and omega1 are pinned to zero on the outer boundary, which is "
            "exact only while the support stays inside; outer_band_fraction "
            "measures the violation and it is reported for every run",
            "the multipole tail bound needs a numerical support radius, so a "
            "threshold is applied and the neglected source mass is reported",
            "the cascade model is not Navier-Stokes: it has no geometry and no "
            "phase, so it can only answer design-feasibility questions",
        ],
        "reproducibility": {"runtime_provenance": provenance},
    }

    paths = {
        "summary.json": summary,
        "config.snapshot.json": config,
    }
    written = []
    for name, payload in paths.items():
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
        "tail_bounds.csv", tail_rows,
        ["truncation", "evaluation_radius", "distance", "value_error", "value_bound",
         "value_ratio", "gradient_error", "gradient_bound", "gradient_ratio",
         "hessian_error", "hessian_bound", "hessian_ratio", "third_bound"],
    )
    write_csv(
        "green_path.csv", green_rows,
        ["nr", "nz"] + [f"relative_{k}" for k in ("psi", "dr", "dz", "drr", "dzz", "drz")]
        + [f"order_{k}" for k in ("psi", "dr", "dz", "drr", "dzz", "drz")],
    )
    write_csv(
        "recovery.csv", recovery_rows,
        ["nr", "nz", "dr", "error_psi", "order_psi", "error_dr", "order_dr",
         "error_dz", "order_dz", "error_u_r", "order_u_r", "error_u_z", "order_u_z",
         "axis_u_r_max", "axis_psi_r_max", "poisson_residual_max",
         "interior_gradient_tail_bound", "velocity_tail_bound"],
    )
    write_csv(
        "evolution.csv", evolution_rows,
        ["family", "integrator", "boundary_mode", "nr", "nz", "r_max", "z_max",
         "steps", "time_step", "final_max_omega1", "final_energy", "final_l3",
         "max_relative_divergence", "max_curl_defect", "max_poisson_residual",
         "max_outer_band_fraction", "max_energy_increase", "boundary_sensitivity",
         "derivative_tail_bound", "velocity_tail_bound", "low_frequency_fraction",
         "neglected_source_fraction", "effective_shell_count", "diverged"],
    )
    write_csv(
        "cascade.csv", cascade_rows,
        ["label", "viscosity", "forced_shells", "forcing_amplitude",
         "initial_bandwidth", "max_bandwidth", "max_active_shell",
         "max_high_shell_amplitude", "direct_high_shell_injection",
         "nonlinear_transfer_defect", "energy_final", "bandwidth_grew", "diverged"],
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

    if not acceptance["all_passed"]:
        raise RuntimeError("whole-space Gate 5 failed acceptance checks")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.config, args.output_dir)
    print(json.dumps(summary["acceptance"], ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
