"""Run the free-radial, padded-z Poisson manufactured-solution gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ns_certificate_lab.free_space_poisson import solve_l5_free_radial_periodized_z


def relative_l2(error: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(error.ravel()) / np.linalg.norm(reference.ravel()))


def run(output: Path) -> dict[str, object]:
    radial_points = 181
    z_points = 80
    radial_max = 6.0
    z_half_width = 2.0
    a = 0.9

    r = np.linspace(0.0, radial_max, radial_points)
    dz = 2.0 * z_half_width / z_points
    z = (np.arange(z_points) - z_points // 2) * dz
    rr, zz = np.meshgrid(r, z, indexing="ij")
    exact = np.exp(-a * (rr**2 + zz**2))
    source = (10.0 * a - 4.0 * a * a * (rr**2 + zz**2)) * exact
    interior = (rr <= 3.2) & (np.abs(zz) <= 1.0)

    cases = []
    previous = None
    for pad_factor in (1, 2, 4):
        computed, diagnostics = solve_l5_free_radial_periodized_z(
            source,
            r,
            dz,
            pad_factor=pad_factor,
        )
        error = relative_l2((computed - exact)[interior], exact[interior])
        case = {
            "pad_factor": pad_factor,
            "relative_l2_interior": error,
            "period_image_reduction_from_previous": (
                None if previous is None else previous / error
            ),
            "diagnostics": diagnostics.__dict__,
        }
        cases.append(case)
        previous = error

    result = {
        "status": "pass"
        if (
            cases[1]["relative_l2_interior"] < 0.5 * cases[0]["relative_l2_interior"]
            and cases[2]["relative_l2_interior"] < cases[1]["relative_l2_interior"]
            and cases[2]["relative_l2_interior"] < 4.5e-3
        )
        else "fail",
        "interpretation": (
            "The radial boundary is removed analytically.  z remains periodized; "
            "padding sensitivity measures, but does not rigorously bound, image error."
        ),
        "cases": cases,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "free_space_poisson_gate.json").write_text(
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
