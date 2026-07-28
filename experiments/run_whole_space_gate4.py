"""Run Gate 4: the linear whole-space elliptic gate with non-periodic ``z``.

``PLAN.md`` Phase 2.8 forbids any whole-space *nonlinear* run until a linear
gate has removed both finite-domain artefacts of the shipped evidence — the
radial Dirichlet wall and the periodic axial direction — against a reference
that is not itself a discretisation.  This experiment is that gate.  It writes

* a **grid-refinement table** against an exact closed-form free-space solution,
* an **independent ``R_max`` / ``Z_max`` enlargement table** whose truncation
  component is isolated by common-interior differencing,
* a **tail-error table** comparing the measured boundary-data error with the
  a posteriori monopole bound,
* the **zero-frequency stress table**, which quantifies exactly how badly a
  periodic axial box misrepresents the far field,
* the **period-image / radial-wall separation**, done at the continuum level
  because at every affordable spacing the measured errors are discretisation
  dominated, and
* an **independent Cartesian divergence and curl audit** of the recovered
  velocity.

Passing licenses coupling the solver into a time-stepper.  It licenses nothing
about the resulting evolution.
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

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.whole_space_gate import (
    BumpSpec,
    ManufacturedFreeSpaceField,
    WholeSpaceEllipticSolver,
    audit_cartesian_recovery,
    axial_spectrum,
    domain_truncation_series,
    monopole_tail_bound,
    separate_period_and_wall_errors,
    solve_and_compare,
    zero_frequency_stress,
)

EXPECTED_SCHEMA_VERSION = 1
EXPECTED_EXPERIMENT = "whole_space_gate4_linear"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("config root must be an object")
    return data


def _plain_int(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    return int(value)


def _finite_number(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    converted = float(value)
    if not math.isfinite(converted):
        raise ValueError(f"{name} must be a finite number")
    return converted


def _validate_config(config: dict[str, Any]) -> None:
    required = {
        "schema_version",
        "experiment",
        "source",
        "balanced_source",
        "interior_radius",
        "grid_refinement",
        "domain_enlargement",
        "zero_frequency",
        "separation",
        "cartesian",
        "fault_injection",
    }
    if set(config) != required:
        raise ValueError("Gate 4 config has missing or unknown keys")
    if _plain_int(config["schema_version"], name="schema_version") != (
        EXPECTED_SCHEMA_VERSION
    ):
        raise ValueError("unsupported Gate 4 config schema")
    if config["experiment"] != EXPECTED_EXPERIMENT:
        raise ValueError("experiment must equal the audited canonical value")
    for key in ("source", "balanced_source"):
        block = config[key]
        if not isinstance(block, dict) or set(block) != {"label", "bumps"}:
            raise ValueError(f"{key} must have exactly label and bumps")
        if not isinstance(block["bumps"], list) or not block["bumps"]:
            raise ValueError(f"{key}.bumps must be a non-empty list")
        for index, bump in enumerate(block["bumps"]):
            if set(bump) != {"center", "radius", "amplitude", "power"}:
                raise ValueError(f"{key}.bumps[{index}] has unexpected keys")
            _finite_number(bump["center"], name="center")
            _finite_number(bump["radius"], name="radius")
            _finite_number(bump["amplitude"], name="amplitude")
            _plain_int(bump["power"], name="power")
    if _finite_number(config["interior_radius"], name="interior_radius") <= 0.0:
        raise ValueError("interior_radius must be positive")


def _field(block: dict[str, Any]) -> ManufacturedFreeSpaceField:
    return ManufacturedFreeSpaceField(
        bumps=tuple(
            BumpSpec(
                center=float(bump["center"]),
                radius=float(bump["radius"]),
                amplitude=float(bump["amplitude"]),
                power=int(bump["power"]),
            )
            for bump in block["bumps"]
        )
    )


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


def _order(coarse: float, fine: float, refinement: float) -> float:
    return math.log(coarse / fine) / math.log(refinement)


def _finite_or_none(value: float) -> float | None:
    return value if math.isfinite(value) else None


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = _load_json(config_path)
    _validate_config(config)
    if output_dir.exists():
        raise FileExistsError(
            f"refusing to overwrite existing output directory: {output_dir}"
        )
    provenance = collect_runtime_provenance()
    output_dir.mkdir(parents=True)

    field = _field(config["source"])
    balanced = _field(config["balanced_source"])
    interior_radius = float(config["interior_radius"])

    # -- table 1: grid refinement against the exact reference ---------------- #
    refinement_rows: list[dict[str, Any]] = []
    for mode in ("zero", "monopole"):
        previous: float | None = None
        for nr in config["grid_refinement"]["radial_points"]:
            row = solve_and_compare(
                field,
                r_max=float(config["grid_refinement"]["extent"]),
                z_max=float(config["grid_refinement"]["extent"]),
                nr=int(nr),
                nz=2 * int(nr) - 1,
                boundary_mode=mode,
                interior_radius=interior_radius,
            )
            payload = row.as_dict()
            payload["observed_order_from_previous"] = (
                None if previous is None else _order(previous, row.interior_relative_l2, 2.0)
            )
            previous = row.interior_relative_l2
            refinement_rows.append(payload)

    # -- table 2: independent domain enlargement ----------------------------- #
    enlargement = config["domain_enlargement"]
    truncation_rows: list[dict[str, Any]] = []
    truncation_rates: dict[str, list[float]] = {}
    for mode in ("zero", "monopole"):
        for enlarge in ("radial", "axial"):
            rows = domain_truncation_series(
                field,
                sizes=[float(value) for value in enlargement["sizes"]],
                fixed_extent=float(enlargement["fixed_extent"]),
                spacing=float(enlargement["spacing"]),
                enlarge=enlarge,
                boundary_mode=mode,
                interior_radius=interior_radius,
                reference_size=float(enlargement["reference_size"]),
            )
            values = [row.common_interior_max_diff for row in rows]
            sizes = [float(value) for value in enlargement["sizes"]]
            rates = [
                _order(values[i], values[i + 1], sizes[i + 1] / sizes[i])
                for i in range(len(values) - 1)
            ]
            truncation_rates[f"{mode}_{enlarge}"] = rates
            for row, size in zip(rows, sizes):
                payload = row.as_dict()
                payload["enlarged"] = enlarge
                payload["enlarged_extent"] = size
                payload["tail_bound"] = _finite_or_none(row.tail_bound)
                truncation_rows.append(payload)

    # -- table 3: tail bound versus measured boundary error ------------------ #
    tail_rows: list[dict[str, Any]] = []
    for row in truncation_rows:
        if row["boundary_mode"] != "monopole":
            continue
        bound = row["tail_bound"]
        tail_rows.append(
            {
                "r_max": row["r_max"],
                "z_max": row["z_max"],
                "enlarged": row["enlarged"],
                "boundary_data_error_max": row["boundary_data_error_max"],
                "common_interior_max_diff": row["common_interior_max_diff"],
                "tail_bound": bound,
                "bound_dominates_boundary_error": bound is not None
                and row["boundary_data_error_max"] <= bound,
                "bound_dominates_interior_difference": bound is not None
                and row["common_interior_max_diff"] <= bound,
                "sharpness": None
                if not bound
                else row["boundary_data_error_max"] / bound,
            }
        )

    # -- table 4: zero-frequency stress -------------------------------------- #
    zero_block = config["zero_frequency"]
    zero_rows = [
        row.as_dict()
        for row in zero_frequency_stress(
            field,
            radial_extents=[float(value) for value in zero_block["radial_extents"]],
            axial_period=float(zero_block["axial_period"]),
        )
    ]
    spectrum_rows: list[dict[str, Any]] = []
    for half in zero_block["axial_half_extents"]:
        half = float(half)
        grid = AxisymmetricGrid.uniform(
            nr=33,
            nz=2 * int(half * 8) + 1,
            r_max=8.0,
            z_min=-half,
            z_max=half,
            periodic_z=False,
        )
        payload = axial_spectrum(grid).as_dict()
        payload["continuum_lowest_eigenvalue"] = (math.pi / (2.0 * half)) ** 2
        spectrum_rows.append(payload)

    # -- table 5: period-image versus radial-wall separation ----------------- #
    separation_block = config["separation"]
    separation_rows = [
        separate_period_and_wall_errors(
            field,
            radial_extent=float(separation_block["radial_extent"]),
            axial_half_period=float(half),
            spacing=float(separation_block["spacing"]),
            interior_radius=float(separation_block["interior_radius"]),
        ).as_dict()
        for half in separation_block["axial_half_periods"]
    ]

    # -- table 6: independent Cartesian audit -------------------------------- #
    cartesian_rows: list[dict[str, Any]] = []
    previous_curl: float | None = None
    previous_divergence: float | None = None
    for points in config["cartesian"]["points"]:
        report = audit_cartesian_recovery(
            field,
            half_width=float(config["cartesian"]["half_width"]),
            points=int(points),
        )
        payload = report.as_dict()
        payload["curl_order_from_previous"] = (
            None if previous_curl is None else _order(previous_curl, report.curl_theta_error_max, 2.0)
        )
        payload["divergence_order_from_previous"] = (
            None
            if previous_divergence is None
            else _order(previous_divergence, report.divergence_max, 2.0)
        )
        previous_curl = report.curl_theta_error_max
        previous_divergence = report.divergence_max
        cartesian_rows.append(payload)

    # -- the balanced (zero net mass) control -------------------------------- #
    balanced_row = solve_and_compare(
        balanced,
        r_max=float(config["grid_refinement"]["extent"]),
        z_max=float(config["grid_refinement"]["extent"]),
        nr=int(config["grid_refinement"]["radial_points"][-1]),
        nz=2 * int(config["grid_refinement"]["radial_points"][-1]) - 1,
        boundary_mode="monopole",
        interior_radius=interior_radius,
    ).as_dict()
    balanced_row["net_five_dimensional_mass"] = balanced.mass

    # -- fault injection ------------------------------------------------------ #
    fault_block = config["fault_injection"]
    fault_rows: list[dict[str, Any]] = []
    for label, kwargs in (
        ("healthy", {}),
        ("axis_coefficient", {"axis_coefficient": float(fault_block["axis_coefficient"])}),
        (
            "radial_drift_coefficient",
            {"radial_drift_coefficient": float(fault_block["radial_drift_coefficient"])},
        ),
    ):
        axis_errors: list[float] = []
        for nr in config["grid_refinement"]["radial_points"]:
            grid = AxisymmetricGrid.uniform(
                nr=int(nr),
                nz=2 * int(nr) - 1,
                r_max=float(config["grid_refinement"]["extent"]),
                z_min=-float(config["grid_refinement"]["extent"]),
                z_max=float(config["grid_refinement"]["extent"]),
                periodic_z=False,
            )
            solver = WholeSpaceEllipticSolver.build(
                grid, boundary_mode="monopole", **kwargs
            )
            r_mesh, z_mesh = grid.mesh()
            error = solver.solve(field.omega1(r_mesh, z_mesh)) - field.psi1(
                r_mesh, z_mesh
            )
            axis_errors.append(
                float(np.max(np.abs(error[0, np.abs(grid.z) <= interior_radius])))
            )
        fault_rows.append(
            {
                "label": label,
                "axis_error_by_resolution": axis_errors,
                "observed_orders": [
                    _order(axis_errors[i], axis_errors[i + 1], 2.0)
                    for i in range(len(axis_errors) - 1)
                ],
            }
        )
    healthy_axis = fault_rows[0]["axis_error_by_resolution"]
    understated_bound = monopole_tail_bound(
        source_l1_norm_5d=field.l1_norm_5d / float(fault_block["understated_source_factor"]),
        support_radius=field.support_radius,
        boundary_min_radius=float(enlargement["sizes"][0]),
    )
    first_monopole_tail = next(
        row for row in tail_rows if row["r_max"] == float(enlargement["sizes"][0])
    )

    # -- preregistered acceptance -------------------------------------------- #
    monopole_refinement_orders = [
        row["observed_order_from_previous"]
        for row in refinement_rows
        if row["boundary_mode"] == "monopole"
        and row["observed_order_from_previous"] is not None
    ]
    zero_refinement_orders = [
        row["observed_order_from_previous"]
        for row in refinement_rows
        if row["boundary_mode"] == "zero"
        and row["observed_order_from_previous"] is not None
    ]
    acceptance = {
        "grid_refinement_second_order_with_monopole_trace": min(
            monopole_refinement_orders
        )
        >= 1.85,
        "zero_trace_saturates_as_expected": min(zero_refinement_orders) < 1.0,
        "discrete_residual_below_1e_9": all(
            row["residual_max"] < 1.0e-9 for row in refinement_rows
        ),
        "radial_truncation_rate_zero_trace_at_least_2_7": min(
            truncation_rates["zero_radial"]
        )
        >= 2.7,
        "radial_truncation_rate_monopole_trace_at_least_3_5": min(
            truncation_rates["monopole_radial"]
        )
        >= 3.5,
        "axial_truncation_rate_monopole_trace_at_least_3_5": min(
            truncation_rates["monopole_axial"]
        )
        >= 3.5,
        "tail_bound_dominates_every_monopole_row": all(
            row["bound_dominates_boundary_error"]
            and row["bound_dominates_interior_difference"]
            for row in tail_rows
        ),
        "no_zero_axial_mode_at_any_extent": all(
            not row["zero_mode_present"] and row["smallest_magnitude"] > 0.0
            for row in spectrum_rows
        ),
        "axial_spectrum_matches_continuum_within_1e_3": all(
            abs(row["smallest_magnitude"] / row["continuum_lowest_eigenvalue"] - 1.0)
            < 1.0e-3
            for row in spectrum_rows
        ),
        "periodic_zero_mode_overstates_far_field_for_large_boxes": zero_rows[-1][
            "ratio"
        ]
        > 1.0,
        "period_image_component_falls_with_the_period": separation_rows[0][
            "period_image_exact"
        ]
        > separation_rows[-1]["period_image_exact"],
        "period_image_sum_is_converged": all(
            row["period_image_tail_bound"] < 1.0e-3 * row["period_image_exact"]
            for row in separation_rows
        ),
        "cartesian_divergence_second_order": (
            cartesian_rows[-1]["divergence_order_from_previous"] >= 1.8
        ),
        "cartesian_curl_second_order": (
            cartesian_rows[-1]["curl_order_from_previous"] >= 1.8
        ),
        "cartesian_curl_relative_below_1e_2": (
            cartesian_rows[-1]["curl_theta_relative"] < 1.0e-2
        ),
        "balanced_source_has_zero_net_mass": abs(
            balanced_row["net_five_dimensional_mass"]
        )
        < 1.0e-15,
        "fault_axis_coefficient_inflates_the_axis_error": all(
            bad > 2.0 * good
            for good, bad in zip(healthy_axis, fault_rows[1]["axis_error_by_resolution"])
        ),
        "fault_radial_drift_stops_convergence": max(fault_rows[2]["observed_orders"])
        < 0.5,
        "fault_understated_source_norm_breaks_the_bound": (
            first_monopole_tail["boundary_data_error_max"] > understated_bound
        ),
    }
    acceptance["all_passed"] = all(acceptance.values())

    summary = {
        "schema_version": 1,
        "experiment": EXPECTED_EXPERIMENT,
        "statement": (
            "linear whole-space elliptic gate with non-periodic z, measured "
            "against an exact closed-form free-space reference; a discretisation "
            "and domain-truncation gate, not a statement about any evolution"
        ),
        "operator": "-(d_rr + 3/r d_r + d_zz) psi1 = omega1",
        "axial_path": (
            "dense discrete sine transform of type I built from numpy.sin; no FFT "
            "and no convention shared with the periodic solvers"
        ),
        "reference": (
            "exact: five-dimensional Newton formula for a superposition of "
            "compactly supported radial bumps, closed form, no quadrature"
        ),
        "source_mass": field.mass,
        "source_support_radius": field.support_radius,
        "source_l1_norm_5d": field.l1_norm_5d,
        "grid_refinement": refinement_rows,
        "domain_enlargement": truncation_rows,
        "domain_enlargement_rates": truncation_rates,
        "tail_error": tail_rows,
        "zero_frequency": zero_rows,
        "axial_spectrum": spectrum_rows,
        "period_wall_separation": separation_rows,
        "cartesian_audit": cartesian_rows,
        "balanced_source_control": balanced_row,
        "fault_injection": fault_rows,
        "understated_source_bound": understated_bound,
        "acceptance": acceptance,
        "coupling_contract": WholeSpaceEllipticSolver.coupling_contract(),
        "limitations": [
            "binary64 throughout, no outward rounding and no interval enclosure",
            "the a posteriori bound covers the continuum truncation error only; "
            "the discretisation error is measured, not bounded, because the "
            "radial row i=1 is not an M-matrix",
            "the manufactured source is C^{power-1}, not C^infinity: a C^infinity "
            "bump has no elementary free-space potential and the reference would "
            "stop being independent",
            "the measured period-image and radial-wall errors are discretisation "
            "dominated at every affordable spacing, which is why the separation "
            "is reported from the exact reference instead",
            "a wrong axis coefficient does not break the observed order, only its "
            "amplitude; an order-only gate would miss that fault",
            "this is a linear gate and licenses no claim about the nonlinear "
            "whole-space evolution",
        ],
        "reproducibility": {"runtime_provenance": provenance},
    }

    summary_path = output_dir / "summary.json"
    refinement_path = output_dir / "grid_refinement.csv"
    enlargement_path = output_dir / "domain_enlargement.csv"
    tail_path = output_dir / "tail_error.csv"
    cartesian_path = output_dir / "cartesian_audit.csv"
    config_snapshot_path = output_dir / "config.snapshot.json"

    _write_json(summary_path, summary)
    _write_json(config_snapshot_path, config)

    def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key) for key in fieldnames})

    _write_csv(
        refinement_path,
        refinement_rows,
        [
            "boundary_mode",
            "r_max",
            "z_max",
            "nr",
            "nz",
            "dr",
            "dz",
            "interior_relative_l2",
            "interior_max_abs",
            "observed_order_from_previous",
            "boundary_data_error_max",
            "residual_max",
        ],
    )
    _write_csv(
        enlargement_path,
        truncation_rows,
        [
            "boundary_mode",
            "enlarged",
            "enlarged_extent",
            "r_max",
            "z_max",
            "common_interior_max_diff",
            "boundary_data_error_max",
            "tail_bound",
        ],
    )
    _write_csv(
        tail_path,
        tail_rows,
        [
            "enlarged",
            "r_max",
            "z_max",
            "boundary_data_error_max",
            "common_interior_max_diff",
            "tail_bound",
            "bound_dominates_boundary_error",
            "bound_dominates_interior_difference",
            "sharpness",
        ],
    )
    _write_csv(
        cartesian_path,
        cartesian_rows,
        [
            "points",
            "velocity_max",
            "divergence_max",
            "divergence_order_from_previous",
            "curl_theta_error_max",
            "curl_theta_relative",
            "curl_order_from_previous",
            "curl_axial_max",
        ],
    )

    manifest = {
        "schema_version": 1,
        "payloads": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in (
                summary_path,
                refinement_path,
                enlargement_path,
                tail_path,
                cartesian_path,
                config_snapshot_path,
            )
        },
    }
    manifest_path = output_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    (output_dir / "manifest.json.sha256").write_text(
        _sha256(manifest_path) + "\n", encoding="ascii"
    )

    if not acceptance["all_passed"]:
        raise RuntimeError("whole-space Gate 4 failed acceptance checks")
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
