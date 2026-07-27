"""Run the finite-cylinder Poisson manufactured-solution gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from ns_certificate_lab.finite_cylinder_poisson import solve_finite_cylinder_poisson
from ns_certificate_lab.grid import AxisymmetricGrid


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("config root must be an object")
    return data


def _write_json(path: Path, data: object) -> None:
    text = json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.write_text(text, encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manufactured(
    grid: AxisymmetricGrid,
    *,
    mode: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r, z = grid.mesh()
    q = 2.0 * np.pi * mode / float(grid.z_period)
    a, b = -0.31, 0.07
    c, d, e = 0.22, -0.19, 0.035
    constant = 1.1 + a * r**2 + b * r**4
    oscillatory_radial = c + d * r**2 + e * r**4
    cosine = np.cos(q * z)
    psi = constant + oscillatory_radial * cosine
    omega = -(8.0 * a + 24.0 * b * r**2)
    omega += (q * q * oscillatory_radial - 8.0 * d - 24.0 * e * r**2) * cosine
    return psi, omega, psi[-1, :]


def _weighted_relative_l2(error: np.ndarray, exact: np.ndarray, grid: AxisymmetricGrid) -> float:
    weights = grid.r[:, None] ** 3
    return float(np.sqrt(np.sum(weights * error**2) / np.sum(weights * exact**2)))


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True)

    config = _load_json(config_path)
    mode = int(config["mode"])
    rows: list[dict[str, Any]] = []
    final_arrays: dict[str, np.ndarray] | None = None

    for item in config["resolutions"]:
        nr = int(item["nr"])
        nz = int(item["nz"])
        grid = AxisymmetricGrid.uniform(
            nr=nr,
            nz=nz,
            r_max=float(config["r_max"]),
            z_min=float(config["z_min"]),
            z_max=float(config["z_max"]),
            periodic_z=True,
        )
        exact, omega, boundary = _manufactured(grid, mode=mode)
        result = solve_finite_cylinder_poisson(
            omega,
            grid,
            outer_boundary=boundary,
            condition_mode_indices=tuple(int(v) for v in config["condition_mode_indices"]),
        )
        error = result.psi - exact
        row = {
            "nr": nr,
            "nz": nz,
            "dr": grid.dr,
            "dz": grid.dz,
            "weighted_relative_l2_error": _weighted_relative_l2(error, exact, grid),
            "max_abs_error": float(np.max(np.abs(error))),
            "fourier_matrix_residual_max": result.diagnostics.fourier_matrix_residual_max,
            "physical_cross_stencil_defect_rms": result.diagnostics.physical_cross_stencil_defect_rms,
            "physical_cross_stencil_defect_max": result.diagnostics.physical_cross_stencil_defect_max,
            "boundary_error_max": result.diagnostics.boundary_error_max,
            "condition_numbers": {
                str(index): value for index, value in result.diagnostics.condition_numbers.items()
            },
        }
        rows.append(row)
        final_arrays = {
            "r": grid.r,
            "z": grid.z,
            "psi_exact": exact,
            "omega": omega,
            "outer_boundary": boundary,
            "psi_numerical": result.psi,
            "error": error,
        }

    for index in range(len(rows) - 1):
        coarse = float(rows[index]["weighted_relative_l2_error"])
        fine = float(rows[index + 1]["weighted_relative_l2_error"])
        h_coarse = float(rows[index]["dr"])
        h_fine = float(rows[index + 1]["dr"])
        rows[index + 1]["observed_order_from_previous"] = float(
            np.log(coarse / fine) / np.log(h_coarse / h_fine)
        )
    rows[0]["observed_order_from_previous"] = None

    orders = [float(row["observed_order_from_previous"]) for row in rows[1:]]
    cross_defect_orders = [
        float(
            np.log(
                float(rows[i]["physical_cross_stencil_defect_max"])
                / float(rows[i + 1]["physical_cross_stencil_defect_max"])
            )
            / np.log(float(rows[i]["dr"]) / float(rows[i + 1]["dr"]))
        )
        for i in range(len(rows) - 1)
    ]
    acceptance = {
        "errors_strictly_decrease": all(
            float(rows[i + 1]["weighted_relative_l2_error"])
            < float(rows[i]["weighted_relative_l2_error"])
            for i in range(len(rows) - 1)
        ),
        "minimum_observed_order_at_least_1_85": min(orders) >= 1.85,
        "cross_stencil_defect_strictly_decreases": all(
            float(rows[i + 1]["physical_cross_stencil_defect_max"])
            < float(rows[i]["physical_cross_stencil_defect_max"])
            for i in range(len(rows) - 1)
        ),
        "minimum_cross_stencil_order_at_least_1_85": min(cross_defect_orders) >= 1.85,
        "maximum_fourier_matrix_residual_below_2e_9": max(
            float(row["fourier_matrix_residual_max"]) for row in rows
        ) < 2.0e-9,
        "maximum_boundary_error_below_2e_13": max(
            float(row["boundary_error_max"]) for row in rows
        ) < 2.0e-13,
    }
    acceptance["all_passed"] = all(acceptance.values())

    summary = {
        "schema_version": 1,
        "experiment": "finite_cylinder_poisson_gate",
        "statement": "floating-point finite-cylinder manufactured-solution gate; not a continuum proof",
        "operator": "-(d_rr + 3/r d_r + d_zz)",
        "axis_row": "8*(psi[0]-psi[1])/dr^2 in the negative operator",
        "outer_boundary": "explicit nonzero Dirichlet data at r=R",
        "rows": rows,
        "acceptance": acceptance,
        "limitations": [
            "binary64 arithmetic without outward rounding",
            "finite cylinder with periodic z",
            "no full-space Green-tail or domain-truncation bound",
            "condition numbers are dense estimates for selected Fourier modes only",
        ],
    }

    summary_path = output_dir / "summary.json"
    csv_path = output_dir / "convergence.csv"
    arrays_path = output_dir / "finest_solution.npz"
    config_snapshot_path = output_dir / "config.snapshot.json"

    _write_json(summary_path, summary)
    _write_json(config_snapshot_path, config)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "nr",
            "nz",
            "dr",
            "dz",
            "weighted_relative_l2_error",
            "max_abs_error",
            "observed_order_from_previous",
            "fourier_matrix_residual_max",
            "physical_cross_stencil_defect_rms",
            "physical_cross_stencil_defect_max",
            "boundary_error_max",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})
    assert final_arrays is not None
    np.savez(arrays_path, **final_arrays)

    manifest = {
        "schema_version": 1,
        "payloads": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in (summary_path, csv_path, arrays_path, config_snapshot_path)
        },
    }
    manifest_path = output_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    (output_dir / "manifest.json.sha256").write_text(_sha256(manifest_path) + "\n", encoding="ascii")

    if not acceptance["all_passed"]:
        raise RuntimeError("finite-cylinder Poisson gate failed acceptance checks")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.config, args.output_dir)
    print(json.dumps(summary["acceptance"], ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
