"""Run the inexpensive independent manufactured-solution audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from ns_certificate_lab._integrity import strict_json_loads
from ns_certificate_lab.artifacts import CandidateDescription, save_candidate
from ns_certificate_lab.axis import check_axis_regularity
from ns_certificate_lab.cartesian import (
    cartesian_curl_consistency_defect,
    reconstruct_cartesian_velocity,
)
from ns_certificate_lab.config import save_run_config, seed_everything
from ns_certificate_lab.convergence import assess_convergence, rms
from ns_certificate_lab.diagnostics import (
    write_diagnostics_csv,
    write_diagnostics_json,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.manufactured import manufactured_fields
from ns_certificate_lab.operators import (
    divergence_physical,
    elliptic_defect,
    recover_velocity,
)
from ns_certificate_lab.pde import compute_pde_residuals
from ns_certificate_lab.provenance import collect_runtime_provenance

REPOSITORY = Path(__file__).resolve().parents[1]


def _inside_repository(path: Path) -> bool:
    try:
        path.resolve().relative_to(REPOSITORY)
    except ValueError:
        return False
    return True


def _load_parameters(path: Path) -> dict[str, Any]:
    parameters = strict_json_loads(
        path.read_text(encoding="utf-8"),
        label="manufactured config",
    )
    if not isinstance(parameters, dict):
        raise ValueError("manufactured config must be a JSON object")
    resolutions = parameters.get("radial_resolutions")
    if (
        not isinstance(resolutions, list)
        or len(resolutions) < 3
        or any(isinstance(value, bool) or not isinstance(value, int) for value in resolutions)
        or any(value < 5 for value in resolutions)
        or resolutions != sorted(set(resolutions))
    ):
        raise ValueError("radial_resolutions must be at least three increasing integers")
    return parameters


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        raise ValueError("config and output paths must remain inside this repository")
    # Snapshot provenance once before any generated file can change git_dirty.
    # Every artifact from this logical run receives this same validated object.
    run_provenance = collect_runtime_provenance()
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite non-empty evidence directory: {output_dir}"
        )
    parameters = _load_parameters(config_path)
    seed = int(parameters["seed"])
    seed_everything(seed)
    viscosity = float(parameters["viscosity"])
    time = float(parameters["time"])
    normalization = parameters["normalization"]
    length_scale = float(normalization["length_scale_L0"])
    velocity_scale = float(normalization["velocity_scale_U0"])
    time_scale = length_scale / velocity_scale
    reynolds_number = velocity_scale * length_scale / viscosity
    records: list[dict[str, Any]] = []
    finest: tuple[AxisymmetricGrid, Any] | None = None

    for nr in parameters["radial_resolutions"]:
        grid = AxisymmetricGrid.uniform(
            nr=nr,
            nz=2 * (nr - 1),
            r_max=float(parameters["r_max"]),
            z_min=float(parameters["z_min"]),
            z_max=float(parameters["z_max"]),
            periodic_z=True,
        )
        fields = manufactured_fields(grid, time=time)
        numerical_u_r, numerical_u_z = recover_velocity(grid, fields.psi1)
        cartesian_velocity = reconstruct_cartesian_velocity(
            grid,
            u1=fields.u1,
            psi1=fields.psi1,
            theta=(0.0, 1.0471975511965976, 1.5707963267948966),
        )
        cartesian_curl_defect = cartesian_curl_consistency_defect(
            grid,
            velocity=cartesian_velocity,
            omega1=fields.omega1,
        )
        residuals = compute_pde_residuals(
            grid,
            u1=fields.u1,
            omega1=fields.omega1,
            psi1=fields.psi1,
            du1_dt=fields.du1_dt,
            domega1_dt=fields.domega1_dt,
            viscosity=viscosity,
        )
        axis = check_axis_regularity(
            grid,
            even_fields={
                "u1": fields.u1,
                "omega1": fields.omega1,
                "psi1": fields.psi1,
                "u_z": numerical_u_z,
            },
            odd_fields={
                "u_r": numerical_u_r,
                "u_theta": grid.r[:, None] * fields.u1,
            },
            tolerance=1e-8,
        )
        u1_forcing = fields.forcing_u1(viscosity)
        omega1_forcing = fields.forcing_omega1(viscosity)
        u1_term_rms = {
            name: rms(value)
            for name, value in residuals.u1.as_dict().items()
        }
        u1_term_rms["forced_error"] = rms(
            residuals.u1.residual - u1_forcing
        )
        omega1_term_rms = {
            name: rms(value)
            for name, value in residuals.omega1.as_dict().items()
        }
        omega1_term_rms["forced_error"] = rms(
            residuals.omega1.residual - omega1_forcing
        )
        records.append(
            {
                "nr": nr,
                "nz": grid.nz,
                "spacing": max(grid.dr, grid.dz),
                "velocity_recovery_rms": rms(numerical_u_r - fields.u_r)
                + rms(numerical_u_z - fields.u_z),
                "physical_divergence_rms": rms(
                    divergence_physical(grid, fields.u_r, fields.u_z)
                ),
                "elliptic_defect_rms": rms(
                    elliptic_defect(grid, fields.psi1, fields.omega1)
                ),
                "cartesian_curl_defect_rms": rms(cartesian_curl_defect),
                "u1_forced_residual_rms": u1_term_rms["forced_error"],
                "omega1_forced_residual_rms": omega1_term_rms["forced_error"],
                "axis_regularity_passed": axis.passed,
                "axis_regularity": axis.to_dict(),
                "pde_term_rms": {
                    "u1": u1_term_rms,
                    "omega1": omega1_term_rms,
                },
            }
        )
        finest = (grid, fields)

    spacings = [float(record["spacing"]) for record in records]
    convergence = {}
    for key in (
        "velocity_recovery_rms",
        "physical_divergence_rms",
        "elliptic_defect_rms",
        "cartesian_curl_defect_rms",
        "u1_forced_residual_rms",
        "omega1_forced_residual_rms",
    ):
        convergence[key] = assess_convergence(
            [float(record[key]) for record in records],
            spacings,
            minimum_order=1.8,
        ).to_dict()

    save_run_config(
        output_dir / "run_config.json",
        parameters=parameters,
        seed=seed,
        provenance=run_provenance,
    )
    write_diagnostics_json(
        output_dir / "diagnostics.json",
        records=records,
        metadata={
            "experiment": "manufactured_solution",
            "convergence": convergence,
            "interpretation": "smooth consistency test; not singularity evidence",
            "provenance": run_provenance,
        },
    )
    csv_records = [
        {
            key: value
            for key, value in record.items()
            if isinstance(value, (str, int, float, bool)) or value is None
        }
        for record in records
    ]
    write_diagnostics_csv(output_dir / "diagnostics.csv", records=csv_records)
    assert finest is not None
    grid, fields = finest
    save_candidate(
        output_dir / "manufactured_candidate.npz",
        grid=grid,
        fields={
            "u1": fields.u1,
            "omega1": fields.omega1,
            "psi1": fields.psi1,
        },
        config=parameters,
        seed=seed,
        description=CandidateDescription(
            representation="nodal little-endian float64 arrays on a tensor grid",
            coordinate_system=(
                "axisymmetric cylindrical half-plane coordinates (r,z), "
                "with periodic z and r>=0"
            ),
            units={
                "r": "dimensionless length",
                "z": "dimensionless length",
                "u1": "dimensionless transformed swirl u^theta/r",
                "omega1": "dimensionless transformed azimuthal vorticity omega^theta/r",
                "psi1": "dimensionless transformed streamfunction psi^theta/r",
            },
            normalization=(
                f"dimensionless with L0={length_scale:.17g}, "
                f"U0={velocity_scale:.17g}, time scale L0/U0="
                f"{time_scale:.17g}, Re=U0*L0/nu={reynolds_number:.17g}"
            ),
            physical_time=time,
            viscosity=viscosity,
            basis_convention=(
                "nodal values ordered (n_r,n_z); u^r=-r*d_z psi1; "
                "u^z=2*psi1+r*d_r psi1; -L5 psi1=omega1"
            ),
        ),
        metadata={
            "kind": "manufactured smooth field",
            "is_singularity_claim": False,
            "generator": {
                "schema": "ns-certificate-lab/manufactured-generator/v1",
                "factory": "ns_certificate_lab.manufactured.manufactured_fields",
            },
            "physical_time": time,
            "viscosity": viscosity,
            "normalization": {
                "length_scale_L0": length_scale,
                "velocity_scale_U0": velocity_scale,
                "time_scale_L0_over_U0": time_scale,
                "reynolds_number_U0_L0_over_nu": reynolds_number,
            },
        },
        provenance=run_provenance,
    )
    return {
        "records": records,
        "convergence": convergence,
        "all_convergence_checks_passed": all(
            report["passed"] for report in convergence.values()
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY / "configs" / "manufactured.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY / "outputs" / "manufactured",
    )
    args = parser.parse_args()
    result = run(args.config.resolve(), args.output_dir.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_convergence_checks_passed"] else 2


if __name__ == "__main__":
    sys.exit(main())
