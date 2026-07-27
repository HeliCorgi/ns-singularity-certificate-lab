from __future__ import annotations

import numpy as np

from ns_certificate_lab.axis import check_axis_regularity
from ns_certificate_lab.convergence import assess_convergence, rms
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.manufactured import manufactured_fields
from ns_certificate_lab.operators import (
    derivative_r,
    derivative_rr,
    derivative_z,
    derivative_zz,
    divergence_physical,
    elliptic_defect,
    recover_velocity,
)


def _grids() -> list[AxisymmetricGrid]:
    return [
        AxisymmetricGrid.uniform(
            nr=nr,
            nz=2 * (nr - 1),
            r_max=1.0,
            periodic_z=True,
        )
        for nr in (17, 33, 65)
    ]


def test_independent_first_and_second_derivatives_converge_at_second_order() -> None:
    errors: list[float] = []
    spacings: list[float] = []
    time = 0.2
    amplitude = np.exp(-time)
    for grid in _grids():
        fields = manufactured_fields(grid, time=time)
        r, _ = grid.mesh()
        exact_rr = amplitude * (-4.0 + 12.0 * r * r) * np.cos(grid.z)[None, :]
        exact_zz = -fields.psi1
        combined_error = sum(
            (
                rms(
                    derivative_r(grid, fields.psi1, even_at_axis=True)
                    - fields.dpsi1_dr
                ),
                rms(derivative_z(grid, fields.psi1) - fields.dpsi1_dz),
                rms(
                    derivative_rr(grid, fields.psi1, even_at_axis=True)
                    - exact_rr
                ),
                rms(derivative_zz(grid, fields.psi1) - exact_zz),
            )
        )
        errors.append(combined_error)
        spacings.append(max(grid.dr, grid.dz))

    report = assess_convergence(errors, spacings, minimum_order=1.8)
    assert report.passed, report.to_dict()


def test_velocity_recovery_converges_to_independent_closed_form() -> None:
    errors: list[float] = []
    spacings: list[float] = []
    for grid in _grids():
        fields = manufactured_fields(grid)
        u_r, u_z = recover_velocity(grid, fields.psi1)
        errors.append(rms(u_r - fields.u_r) + rms(u_z - fields.u_z))
        spacings.append(max(grid.dr, grid.dz))

    report = assess_convergence(errors, spacings, minimum_order=1.8)
    assert report.passed, report.to_dict()


def test_physical_three_dimensional_divergence_converges_to_zero() -> None:
    errors: list[float] = []
    spacings: list[float] = []
    for grid in _grids():
        fields = manufactured_fields(grid)
        defect = divergence_physical(grid, fields.u_r, fields.u_z)
        errors.append(rms(defect))
        spacings.append(max(grid.dr, grid.dz))

    report = assess_convergence(errors, spacings, minimum_order=1.8)
    assert report.passed, report.to_dict()
    assert errors[-1] < 1e-3


def test_elliptic_relation_converges_to_independent_exact_omega1() -> None:
    errors: list[float] = []
    spacings: list[float] = []
    for grid in _grids():
        fields = manufactured_fields(grid)
        errors.append(rms(elliptic_defect(grid, fields.psi1, fields.omega1)))
        spacings.append(max(grid.dr, grid.dz))

    report = assess_convergence(errors, spacings, minimum_order=1.8)
    assert report.passed, report.to_dict()


def test_manufactured_fields_obey_axis_even_odd_conditions() -> None:
    grid = AxisymmetricGrid.uniform(nr=33, nz=64, r_max=1.0)
    fields = manufactured_fields(grid)
    u_r, u_z = recover_velocity(grid, fields.psi1)
    report = check_axis_regularity(
        grid,
        even_fields={
            "psi1": fields.psi1,
            "omega1": fields.omega1,
            "u1": fields.u1,
            "u_z": u_z,
        },
        odd_fields={
            "u_r": u_r,
            "u_theta": grid.r[:, None] * fields.u1,
        },
        tolerance=1e-9,
    )
    assert report.passed, report.to_dict()
