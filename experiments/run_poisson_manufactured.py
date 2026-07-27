#!/usr/bin/env python
"""Run the independent finite-cylinder Poisson manufactured control.

This is a smooth elliptic convergence test.  It is neither a Navier--Stokes
singularity candidate nor a continuum proof.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_file,
    strict_json_loads,
    write_with_digest,
)
from ns_certificate_lab.diagnostics import (
    write_diagnostics_csv,
    write_diagnostics_json,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.poisson import (
    PoissonSolution,
    solve_streamfunction_poisson,
)
from ns_certificate_lab.provenance import collect_runtime_provenance

REPOSITORY = Path(__file__).resolve().parents[1]
EXPECTED_EXPERIMENT_ID = "independent_poisson_manufactured_control"
EXPECTED_COORDINATE_SYSTEM = "axisymmetric finite cylinder (r,z), periodic z"
EXPECTED_EQUATION = "-(d_rr + 3/r d_r + d_zz) psi1 = omega1"


def _inside_repository(path: Path) -> bool:
    try:
        path.resolve().relative_to(REPOSITORY)
    except ValueError:
        return False
    return True


def _finite_positive(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite positive number")
    converted = float(value)
    if not math.isfinite(converted) or converted <= 0.0:
        raise ValueError(f"{name} must be a finite positive number")
    return converted


def _validate_config(config: dict[str, Any]) -> None:
    required = {
        "schema_version",
        "experiment_id",
        "coordinate_system",
        "equation",
        "seed",
        "r_max",
        "z_min",
        "z_max",
        "radial_resolutions",
        "axial_points_per_radial_interval",
        "acceptance",
    }
    if set(config) != required:
        raise ValueError("Poisson config has missing or unknown keys")
    if config["schema_version"] != 1:
        raise ValueError("unsupported Poisson config schema")
    canonical_labels = {
        "experiment_id": EXPECTED_EXPERIMENT_ID,
        "coordinate_system": EXPECTED_COORDINATE_SYSTEM,
        "equation": EXPECTED_EQUATION,
    }
    for name, expected in canonical_labels.items():
        if config[name] != expected:
            raise ValueError(f"{name} must equal the audited canonical value")
    seed = config["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2**32 - 1:
        raise ValueError("seed must be an integer in [0, 2**32 - 1]")
    _finite_positive(config["r_max"], name="r_max")
    z_min = config["z_min"]
    z_max = config["z_max"]
    if (
        isinstance(z_min, bool)
        or isinstance(z_max, bool)
        or not isinstance(z_min, (int, float))
        or not isinstance(z_max, (int, float))
        or not math.isfinite(float(z_min))
        or not math.isfinite(float(z_max))
        or float(z_max) <= float(z_min)
    ):
        raise ValueError("z_min and z_max must define a finite positive period")
    resolutions = config["radial_resolutions"]
    if (
        not isinstance(resolutions, list)
        or len(resolutions) < 3
        or any(isinstance(value, bool) or not isinstance(value, int) for value in resolutions)
        or any(value < 5 for value in resolutions)
        or resolutions != sorted(set(resolutions))
    ):
        raise ValueError(
            "radial_resolutions must contain at least three increasing integers"
        )
    multiplier = config["axial_points_per_radial_interval"]
    if isinstance(multiplier, bool) or not isinstance(multiplier, int) or multiplier < 1:
        raise ValueError("axial_points_per_radial_interval must be a positive integer")
    acceptance = config["acceptance"]
    if not isinstance(acceptance, dict) or set(acceptance) != {
        "minimum_observed_order",
        "minimum_independent_residual_order",
        "maximum_finest_rms_error",
        "maximum_finest_independent_residual_rms",
        "maximum_relative_discrete_residual",
        "maximum_relative_imaginary_leakage",
        "maximum_boundary_defect",
    }:
        raise ValueError("Poisson acceptance object is invalid")
    for name, value in acceptance.items():
        _finite_positive(value, name=f"acceptance {name}")


def _manufactured(
    grid: AxisymmetricGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an analytic profile with Fourier modes 0, 1, and 2."""

    r, z = grid.mesh()
    fundamental_wave_number = 2.0 * np.pi / grid.z_period
    phase = fundamental_wave_number * (z - float(grid.z[0]))
    p1 = 1.0 + 0.2 * r**2 + 0.1 * r**4
    p2 = 0.3 * (1.0 - r**2 + 0.25 * r**4)
    psi = (
        0.15 * r**2
        + p1 * np.cos(phase)
        + p2 * np.sin(2.0 * phase)
    )
    # L5(r**(2q)) = 4*q*(q+1)*r**(2q-2), derived analytically.
    omega = (
        -1.2
        + (
            fundamental_wave_number**2 * p1
            - (1.6 + 2.4 * r**2)
        )
        * np.cos(phase)
        + (
            4.0 * fundamental_wave_number**2 * p2
            - 0.3 * (-8.0 + 6.0 * r**2)
        )
        * np.sin(2.0 * phase)
    )
    return psi, omega, psi[-1].copy()


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=np.float64) ** 2)))


def _orders(errors: list[float], spacings: list[float]) -> list[float]:
    if len(errors) != len(spacings) or len(errors) < 2:
        raise ValueError("at least two errors and spacings are required")
    result: list[float] = []
    for coarse_error, fine_error, coarse_h, fine_h in zip(
        errors[:-1],
        errors[1:],
        spacings[:-1],
        spacings[1:],
    ):
        if min(coarse_error, fine_error, coarse_h, fine_h) <= 0.0:
            raise ValueError("convergence inputs must be positive")
        result.append(
            float(
                math.log(coarse_error / fine_error)
                / math.log(coarse_h / fine_h)
            )
        )
    return result


def evaluate(
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, np.ndarray], list[float]]:
    """Evaluate all resolutions without writing evidence."""

    _validate_config(config)
    np.random.seed(int(config["seed"]))
    records: list[dict[str, Any]] = []
    spacings: list[float] = []
    errors: list[float] = []
    finest: dict[str, np.ndarray] | None = None
    multiplier = int(config["axial_points_per_radial_interval"])
    z_min = float(config["z_min"])
    z_max = float(config["z_max"])

    for nr in config["radial_resolutions"]:
        grid = AxisymmetricGrid.uniform(
            nr=int(nr),
            nz=multiplier * (int(nr) - 1),
            r_max=float(config["r_max"]),
            z_min=z_min,
            z_max=z_max,
            periodic_z=True,
        )
        exact, omega, boundary = _manufactured(grid)
        solution: PoissonSolution = solve_streamfunction_poisson(
            grid,
            omega,
            boundary,
        )
        error = solution.psi1 - exact
        omega_scale = max(_rms(omega[:-1]), np.finfo(np.float64).tiny)
        record = {
            "nr": grid.nr,
            "nz": grid.nz,
            "dr": grid.dr,
            "dz": grid.dz,
            "psi1_rms_error": _rms(error),
            "psi1_max_abs_error": float(np.max(np.abs(error))),
            "relative_discrete_residual_rms": (
                float(solution.metadata["discrete_residual_rms_interior"])
                / omega_scale
            ),
            "discrete_residual_max_abs": float(
                solution.metadata["discrete_residual_max_abs_interior"]
            ),
            "independent_pde_residual_rms": float(
                solution.metadata["independent_pde_residual_rms_interior"]
            ),
            "independent_pde_residual_max_abs": float(
                solution.metadata["independent_pde_residual_max_abs_interior"]
            ),
            "outer_boundary_max_abs_defect": float(
                solution.metadata["outer_boundary_max_abs_defect"]
            ),
            "axis_radial_coefficient_dimensionless": float(
                solution.metadata["axis_radial_coefficient_dimensionless"]
            ),
            "zero_mode_condition_number_inf": float(
                solution.metadata["zero_mode_condition_number_inf"]
            ),
            "solution_relative_imaginary_leakage": float(
                solution.metadata["solution_relative_imaginary_leakage"]
            ),
            "residual_relative_imaginary_leakage": float(
                solution.metadata["residual_relative_imaginary_leakage"]
            ),
        }
        records.append(record)
        spacings.append(grid.dr)
        errors.append(record["psi1_rms_error"])
        finest = {
            "r": grid.r,
            "z": grid.z,
            "psi1_numerical": solution.psi1,
            "psi1_exact": exact,
            "omega1": omega,
            "outer_dirichlet": boundary,
            "discrete_residual": solution.discrete_residual,
            "independent_pde_residual": solution.pde_residual,
        }
    assert finest is not None
    observed_orders = _orders(errors, spacings)
    for index, record in enumerate(records):
        record["observed_order_from_previous"] = (
            None if index == 0 else observed_orders[index - 1]
        )
    return records, finest, observed_orders


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Run the control and write checksummed JSON, CSV and NPZ evidence."""

    _validate_config(config)
    if not _inside_repository(output_dir):
        raise ValueError("output directory must remain inside this repository")
    provenance = collect_runtime_provenance()
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite non-empty evidence directory: {output_dir}"
        )
    records, finest, observed_orders = evaluate(config)
    acceptance = config["acceptance"]
    independent_residuals = [
        float(record["independent_pde_residual_rms"]) for record in records
    ]
    spacings = [float(record["dr"]) for record in records]
    independent_residual_orders = _orders(independent_residuals, spacings)
    finite_values = [
        float(value)
        for record in records
        for value in record.values()
        if isinstance(value, (int, float)) and value is not None
    ] + observed_orders + independent_residual_orders
    checks = {
        "all_diagnostics_finite": all(math.isfinite(value) for value in finite_values),
        "analytic_errors_decrease": all(
            fine["psi1_rms_error"] < coarse["psi1_rms_error"]
            for coarse, fine in zip(records[:-1], records[1:])
        ),
        "observed_order_passes": (
            min(observed_orders)
            >= float(acceptance["minimum_observed_order"])
        ),
        "finest_error_passes": (
            records[-1]["psi1_rms_error"]
            <= float(acceptance["maximum_finest_rms_error"])
        ),
        "independent_pde_residual_decreases": all(
            fine < coarse
            for coarse, fine in zip(
                independent_residuals[:-1],
                independent_residuals[1:],
            )
        ),
        "independent_pde_residual_order_passes": (
            min(independent_residual_orders)
            >= float(acceptance["minimum_independent_residual_order"])
        ),
        "finest_independent_pde_residual_passes": (
            independent_residuals[-1]
            <= float(
                acceptance["maximum_finest_independent_residual_rms"]
            )
        ),
        "discrete_residual_passes": (
            max(record["relative_discrete_residual_rms"] for record in records)
            <= float(acceptance["maximum_relative_discrete_residual"])
        ),
        "fft_imaginary_leakage_passes": (
            max(
                max(
                    record["solution_relative_imaginary_leakage"],
                    record["residual_relative_imaginary_leakage"],
                )
                for record in records
            )
            <= float(acceptance["maximum_relative_imaginary_leakage"])
        ),
        "boundary_defect_passes": (
            max(record["outer_boundary_max_abs_defect"] for record in records)
            <= float(acceptance["maximum_boundary_defect"])
        ),
        "axis_coefficient_is_eight": all(
            record["axis_radial_coefficient_dimensionless"] == 8.0
            for record in records
        ),
    }
    summary = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Smooth manufactured elliptic convergence control; not a "
            "continuum proof and not singularity evidence."
        ),
        "method": {
            "equation": config["equation"],
            "radial": "second-order r^3-flux finite volume",
            "rhs_sampling": (
                "nodal omega1 values stand in for radial control-volume "
                "averages; this is second-order consistent for the smooth "
                "even manufactured field"
            ),
            "axis_limit": "8*(psi1[1]-psi1[0])/dr^2",
            "axial": "NumPy FFT Fourier coefficients",
            "linear_solve": "locally implemented per-mode Thomas algorithm",
            "outer_boundary": "nonzero z-dependent analytic Dirichlet trace",
            "manufactured_modes": [0, 1, 2],
            "condition_estimate": (
                "Unweighted coordinate-basis induced infinity-norm condition "
                "number of the zero-Fourier-mode radial matrix only; not a "
                "weighted coercivity bound"
            ),
        },
        "observed_orders": observed_orders,
        "independent_pde_residual_orders": independent_residual_orders,
        "records": records,
        "checks": checks,
        "accepted_as_smooth_poisson_control": bool(all(checks.values())),
        "reproducibility": {
            "seed": int(config["seed"]),
            "runtime_provenance": provenance,
        },
    }

    write_with_digest(
        output_dir / "config.snapshot.json",
        canonical_json_bytes(config),
    )
    write_diagnostics_json(
        output_dir / "diagnostics.json",
        records=records,
        metadata={
            "experiment": config["experiment_id"],
            "observed_orders": observed_orders,
            "checks": checks,
            "provenance": provenance,
        },
    )
    write_diagnostics_csv(
        output_dir / "diagnostics.csv",
        records=records,
    )
    np.savez_compressed(output_dir / "finest_fields.npz", **finest)
    write_with_digest(
        output_dir / "summary.json",
        canonical_json_bytes(summary),
    )

    artifact_names = [
        "config.snapshot.json",
        "diagnostics.csv",
        "diagnostics.json",
        "finest_fields.npz",
        "summary.json",
    ]
    manifest = {
        "schema": "ns-certificate-lab/poisson-manufactured-manifest/v1",
        "experiment_id": config["experiment_id"],
        "files": {
            name: {
                "sha256": sha256_file(output_dir / name),
                "bytes": (output_dir / name).stat().st_size,
            }
            for name in artifact_names
        },
    }
    write_with_digest(
        output_dir / "manifest.json",
        canonical_json_bytes(manifest),
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY / "configs" / "poisson_manufactured.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY / "outputs" / "poisson_manufactured",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    try:
        config = strict_json_loads(
            config_path.read_text(encoding="utf-8"),
            label="Poisson manufactured config",
        )
        if not isinstance(config, dict):
            raise ValueError("Poisson manufactured config must be a JSON object")
        summary = run(config, output_dir)
    except (OSError, ValueError, FileExistsError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["accepted_as_smooth_poisson_control"] else 2


if __name__ == "__main__":
    sys.exit(main())
