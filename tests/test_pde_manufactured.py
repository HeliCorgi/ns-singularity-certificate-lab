from __future__ import annotations

import numpy as np

from ns_certificate_lab.convergence import assess_convergence, rms
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.manufactured import manufactured_fields
from ns_certificate_lab.pde import compute_pde_residuals


def test_forced_pde_residuals_converge_at_second_order() -> None:
    viscosity = 0.05
    u1_errors: list[float] = []
    omega1_errors: list[float] = []
    spacings: list[float] = []
    for nr in (17, 33, 65):
        grid = AxisymmetricGrid.uniform(
            nr=nr,
            nz=2 * (nr - 1),
            r_max=1.0,
        )
        fields = manufactured_fields(grid, time=0.2)
        residuals = compute_pde_residuals(
            grid,
            u1=fields.u1,
            omega1=fields.omega1,
            psi1=fields.psi1,
            du1_dt=fields.du1_dt,
            domega1_dt=fields.domega1_dt,
            viscosity=viscosity,
        )
        u1_errors.append(
            rms(residuals.u1.residual - fields.forcing_u1(viscosity))
        )
        omega1_errors.append(
            rms(residuals.omega1.residual - fields.forcing_omega1(viscosity))
        )
        spacings.append(max(grid.dr, grid.dz))

    u1_report = assess_convergence(u1_errors, spacings, minimum_order=1.8)
    omega1_report = assess_convergence(omega1_errors, spacings, minimum_order=1.8)
    assert u1_report.passed, u1_report.to_dict()
    assert omega1_report.passed, omega1_report.to_dict()


def test_residual_is_exact_signed_sum_of_individual_terms() -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=32, r_max=1.0)
    fields = manufactured_fields(grid)
    result = compute_pde_residuals(
        grid,
        u1=fields.u1,
        omega1=fields.omega1,
        psi1=fields.psi1,
        du1_dt=fields.du1_dt,
        domega1_dt=fields.domega1_dt,
        viscosity=0.01,
    )
    for equation in (result.u1, result.omega1):
        independently_assembled = (
            equation.time
            + equation.advection_r
            + equation.advection_z
            - equation.source
            - equation.diffusion
        )
        assert np.array_equal(equation.residual, independently_assembled)
