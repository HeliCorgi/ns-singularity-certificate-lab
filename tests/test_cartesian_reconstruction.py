from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ns_certificate_lab.artifacts import (
    CandidateDescription,
    load_candidate,
    save_candidate,
)
from ns_certificate_lab.cartesian import (
    cartesian_curl_consistency_defect,
    curl_consistency_defect,
    reconstruct_candidate_cartesian_velocity,
    reconstruct_cartesian_velocity,
    recover_meridional_velocity_independent,
)
from ns_certificate_lab.convergence import assess_convergence, max_abs, rms
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.manufactured import manufactured_fields


def _description() -> CandidateDescription:
    return CandidateDescription(
        representation="nodal little-endian float64 arrays",
        coordinate_system="axisymmetric cylindrical half-plane (r,z)",
        units={
            "r": "dimensionless length",
            "z": "dimensionless length",
            "u1": "dimensionless transformed swirl",
            "omega1": "dimensionless transformed azimuthal vorticity",
            "psi1": "dimensionless transformed streamfunction",
        },
        normalization="dimensionless manufactured normalization",
        physical_time=0.2,
        viscosity=0.05,
        basis_convention=(
            "array axes (r,z); u^r=-r psi1_z; "
            "u^z=2 psi1+r psi1_r; -L5 psi1=omega1"
        ),
    )


def test_independent_physical_curl_converges_to_r_omega1() -> None:
    errors: list[float] = []
    spacings: list[float] = []
    for nr in (17, 33, 65):
        grid = AxisymmetricGrid.uniform(
            nr=nr,
            nz=2 * (nr - 1),
            r_max=1.0,
            periodic_z=True,
        )
        fields = manufactured_fields(grid, time=0.2)
        u_r, u_z = recover_meridional_velocity_independent(grid, fields.psi1)
        defect = curl_consistency_defect(
            grid,
            u_r=u_r,
            u_z=u_z,
            omega1=fields.omega1,
        )
        errors.append(rms(defect))
        spacings.append(max(grid.dr, grid.dz))

    report = assess_convergence(errors, spacings, minimum_order=1.8)
    assert report.passed, report.to_dict()


def test_cartesian_reconstruction_matches_cylindrical_map_and_axis() -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=32, r_max=1.0)
    fields = manufactured_fields(grid)
    angles = np.array([0.0, 0.5 * np.pi])
    velocity = reconstruct_cartesian_velocity(
        grid,
        u1=fields.u1,
        psi1=fields.psi1,
        theta=angles,
    )
    u_r, u_z = recover_meridional_velocity_independent(grid, fields.psi1)
    u_theta = grid.r[:, None] * fields.u1

    assert velocity.x.shape == (grid.nr, grid.nz, 2)
    assert np.allclose(velocity.u_x[:, :, 0], u_r)
    assert np.allclose(velocity.u_y[:, :, 0], u_theta)
    assert np.allclose(velocity.u_x[:, :, 1], -u_theta)
    assert np.allclose(velocity.u_y[:, :, 1], u_r)
    assert np.allclose(velocity.u_z[:, :, 0], u_z)
    assert max_abs(velocity.u_x[0]) == 0.0
    assert max_abs(velocity.u_y[0]) == 0.0


def test_saved_candidate_can_be_reconstructed_in_cartesian_space(
    scratch_dir: Path,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=32, r_max=1.0)
    fields = manufactured_fields(grid)
    path = scratch_dir / "candidate.npz"
    save_candidate(
        path,
        grid=grid,
        fields={
            "u1": fields.u1,
            "omega1": fields.omega1,
            "psi1": fields.psi1,
        },
        config={"kind": "manufactured"},
        seed=1729,
        description=_description(),
    )
    loaded = load_candidate(path)
    reconstructed = reconstruct_candidate_cartesian_velocity(
        loaded,
        theta=np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False),
    )
    assert reconstructed.u_x.shape == (17, 32, 8)
    assert np.all(np.isfinite(reconstructed.u_x))
    assert np.all(np.isfinite(reconstructed.u_y))
    assert np.all(np.isfinite(reconstructed.u_z))
    defect = cartesian_curl_consistency_defect(
        loaded.grid,
        velocity=reconstructed,
        omega1=loaded.fields["omega1"],
    )
    assert rms(defect) < 0.01


def test_independent_curl_detects_wrong_vorticity_sign() -> None:
    grid = AxisymmetricGrid.uniform(nr=65, nz=128, r_max=1.0)
    fields = manufactured_fields(grid)
    u_r, u_z = recover_meridional_velocity_independent(grid, fields.psi1)
    correct = curl_consistency_defect(
        grid,
        u_r=u_r,
        u_z=u_z,
        omega1=fields.omega1,
    )
    wrong = curl_consistency_defect(
        grid,
        u_r=u_r,
        u_z=u_z,
        omega1=-fields.omega1,
    )
    assert rms(correct) < 0.01
    assert rms(wrong) > 1.0


def test_independent_fourth_order_audit_rejects_too_small_radial_grid() -> None:
    grid = AxisymmetricGrid.uniform(nr=4, nz=8, r_max=1.0)
    field = np.zeros(grid.shape)
    with pytest.raises(ValueError, match="nr >= 5"):
        recover_meridional_velocity_independent(grid, field)
