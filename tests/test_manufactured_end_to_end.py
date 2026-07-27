from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from experiments.run_manufactured import run
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.axis import check_axis_regularity
from ns_certificate_lab.cartesian import (
    cartesian_curl_consistency_defect,
    reconstruct_candidate_cartesian_velocity,
)
from ns_certificate_lab.config import load_run_config
from ns_certificate_lab.convergence import rms
from ns_certificate_lab.diagnostics import load_diagnostics_json
from ns_certificate_lab.manufactured import manufactured_fields
from ns_certificate_lab.operators import recover_velocity
from ns_certificate_lab.pde import compute_pde_residuals

REPOSITORY = Path(__file__).resolve().parents[1]


def test_saved_candidate_reconstructs_detailed_manufactured_diagnostics(
    scratch_dir: Path,
) -> None:
    result = run(
        REPOSITORY / "configs" / "manufactured.json",
        scratch_dir,
    )
    assert result["all_convergence_checks_passed"]
    diagnostics = load_diagnostics_json(scratch_dir / "diagnostics.json")
    assert diagnostics["metadata"]["provenance"]["source_fingerprint_sha256"]
    finest_record = diagnostics["records"][-1]

    loaded = load_candidate(scratch_dir / "manufactured_candidate.npz")
    run_config = load_run_config(scratch_dir / "run_config.json")
    assert loaded.description is not None
    assert loaded.provenance is not None
    assert (
        run_config["provenance"]
        == diagnostics["metadata"]["provenance"]
        == loaded.provenance
    )
    assert loaded.metadata["generator"] == {
        "factory": "ns_certificate_lab.manufactured.manufactured_fields",
        "schema": "ns-certificate-lab/manufactured-generator/v1",
    }
    time = loaded.description.physical_time
    viscosity = loaded.description.viscosity
    assert loaded.metadata["physical_time"] == time
    assert loaded.metadata["viscosity"] == viscosity
    assert loaded.metadata["normalization"] == {
        "length_scale_L0": 1.0,
        "velocity_scale_U0": 1.0,
        "time_scale_L0_over_U0": 1.0,
        "reynolds_number_U0_L0_over_nu": 20.0,
    }
    assert "L0=1" in loaded.description.normalization
    assert "Re=U0*L0/nu=20" in loaded.description.normalization

    regenerated = manufactured_fields(loaded.grid, time=time)
    for name in ("u1", "omega1", "psi1"):
        assert np.array_equal(loaded.fields[name], getattr(regenerated, name))

    cartesian_velocity = reconstruct_candidate_cartesian_velocity(
        loaded,
        theta=(0.0, 1.0471975511965976, 1.5707963267948966),
    )
    cartesian_curl_defect = cartesian_curl_consistency_defect(
        loaded.grid,
        velocity=cartesian_velocity,
        omega1=loaded.fields["omega1"],
    )
    assert finest_record["cartesian_curl_defect_rms"] == pytest.approx(
        rms(cartesian_curl_defect)
    )

    u_r, u_z = recover_velocity(loaded.grid, loaded.fields["psi1"])
    axis = check_axis_regularity(
        loaded.grid,
        even_fields={
            "u1": loaded.fields["u1"],
            "omega1": loaded.fields["omega1"],
            "psi1": loaded.fields["psi1"],
            "u_z": u_z,
        },
        odd_fields={
            "u_r": u_r,
            "u_theta": loaded.grid.r[:, None] * loaded.fields["u1"],
        },
        tolerance=1e-8,
    )
    recorded_axis = finest_record["axis_regularity"]
    assert recorded_axis["passed"] == axis.passed
    assert recorded_axis["tolerance"] == axis.tolerance
    assert recorded_axis["violations"] == list(axis.violations)
    assert recorded_axis["metrics"] == pytest.approx(axis.metrics)

    residuals = compute_pde_residuals(
        loaded.grid,
        u1=loaded.fields["u1"],
        omega1=loaded.fields["omega1"],
        psi1=loaded.fields["psi1"],
        du1_dt=regenerated.du1_dt,
        domega1_dt=regenerated.domega1_dt,
        viscosity=viscosity,
    )
    required_terms = {
        "time",
        "advection_r",
        "advection_z",
        "source",
        "diffusion",
        "residual",
        "forced_error",
    }
    for name, equation, forcing in (
        ("u1", residuals.u1, regenerated.forcing_u1(viscosity)),
        (
            "omega1",
            residuals.omega1,
            regenerated.forcing_omega1(viscosity),
        ),
    ):
        recorded = finest_record["pde_term_rms"][name]
        assert set(recorded) == required_terms
        for term_name, values in equation.as_dict().items():
            assert recorded[term_name] == pytest.approx(rms(values))
        assert recorded["forced_error"] == pytest.approx(
            rms(equation.residual - forcing)
        )
