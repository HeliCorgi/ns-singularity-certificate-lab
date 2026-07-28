"""Quantify the periodic-low-frequency transition in radial wall sensitivity."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from ns_certificate_lab.finite_cylinder_poisson import solve_finite_cylinder_poisson
from ns_certificate_lab.free_space_poisson import solve_l5_free_radial_periodized_z
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.wall_sensitivity import predicted_zero_dirichlet_correction


def relative_l2(error: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(error.ravel()) / np.linalg.norm(reference.ravel()))


def run(output: Path) -> dict[str, object]:
    rows = []
    for period in (1.0, 2.0, 4.0, 8.0, 16.0):
        grid = AxisymmetricGrid.uniform(
            nr=257,
            nz=64,
            r_max=4.0,
            z_min=0.0,
            z_max=period,
            periodic_z=True,
        )
        rr, zz = grid.mesh()
        k = 2.0 * math.pi / period
        source = np.exp(-rr**2) * np.cos(k * zz)
        free, _ = solve_l5_free_radial_periodized_z(source, grid.r, grid.dz)
        finite = solve_finite_cylinder_poisson(
            source,
            grid,
            outer_boundary=0.0,
            condition_mode_indices=(),
        ).psi
        radial_prediction = predicted_zero_dirichlet_correction(
            grid.r,
            k,
            float(free[-1, 0]),
        )
        prediction = radial_prediction[:, None] * np.cos(k * grid.z)[None, :]
        interior = grid.r <= 2.0
        wall_error = relative_l2((finite - free)[interior], free[interior])
        prediction_defect = relative_l2(
            (finite - free - prediction)[interior],
            (finite - free)[interior],
        )
        rows.append(
            {
                "period": period,
                "k": k,
                "kR": k * grid.r[-1],
                "free_boundary_amplitude": float(free[-1, 0]),
                "relative_wall_error_interior": wall_error,
                "modal_correction_prediction_defect": prediction_defect,
            }
        )

    zero_grid = AxisymmetricGrid.uniform(
        nr=257,
        nz=16,
        r_max=4.0,
        z_min=0.0,
        z_max=8.0,
        periodic_z=True,
    )
    rr, _ = zero_grid.mesh()
    zero_source = np.exp(-rr**2)
    zero_free, _ = solve_l5_free_radial_periodized_z(
        zero_source,
        zero_grid.r,
        zero_grid.dz,
    )
    zero_finite = solve_finite_cylinder_poisson(
        zero_source,
        zero_grid,
        outer_boundary=0.0,
        condition_mode_indices=(),
    ).psi
    zero_prediction = predicted_zero_dirichlet_correction(
        zero_grid.r,
        0.0,
        float(zero_free[-1, 0]),
    )[:, None]
    interior = zero_grid.r <= 2.0
    zero_row = {
        "period": None,
        "k": 0.0,
        "kR": 0.0,
        "free_boundary_amplitude": float(zero_free[-1, 0]),
        "relative_wall_error_interior": relative_l2(
            (zero_finite - zero_free)[interior], zero_free[interior]
        ),
        "modal_correction_prediction_defect": relative_l2(
            (zero_finite - zero_free - zero_prediction)[interior],
            (zero_finite - zero_free)[interior],
        ),
    }

    result = {
        "status": "pass"
        if (
            rows[-2]["relative_wall_error_interior"]
            > 5.0 * rows[-3]["relative_wall_error_interior"]
            and rows[-1]["relative_wall_error_interior"]
            > 4.0 * rows[-2]["relative_wall_error_interior"]
            and zero_row["relative_wall_error_interior"]
            > rows[-1]["relative_wall_error_interior"]
            and zero_row["modal_correction_prediction_defect"] < 5.0e-4
        )
        else "fail",
        "interpretation": (
            "Fixed nonzero Fourier modes have exponentially weak radial-wall "
            "influence when kR is large.  As the axial period grows, k falls "
            "and the wall error rises.  The zero mode has an algebraic R^-2 tail."
        ),
        "nonzero_modes": rows,
        "zero_mode": zero_row,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "low_frequency_wall_gate.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
