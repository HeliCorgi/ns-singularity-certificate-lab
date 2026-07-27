from __future__ import annotations

import numpy as np
import pytest

from ns_certificate_lab.finite_cylinder_poisson import (
    apply_finite_cylinder_system,
    assemble_radial_mode_matrix,
    axial_wavenumbers,
    solve_finite_cylinder_poisson,
)
from ns_certificate_lab.grid import AxisymmetricGrid


def _manufactured_fields(
    grid: AxisymmetricGrid,
    *,
    mode: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r, z = grid.mesh()
    q = 2.0 * np.pi * mode / float(grid.z_period)

    a, b = -0.31, 0.07
    c, d, e = 0.22, -0.19, 0.035
    radial_constant = 1.1 + a * r**2 + b * r**4
    radial_mode = c + d * r**2 + e * r**4
    cosine = np.cos(q * z)
    psi = radial_constant + radial_mode * cosine

    omega_constant = -(8.0 * a + 24.0 * b * r**2)
    omega_mode = (q * q * radial_mode - 8.0 * d - 24.0 * e * r**2) * cosine
    omega = omega_constant + omega_mode
    return psi, omega, psi[-1, :]


def _weighted_relative_l2(error: np.ndarray, exact: np.ndarray, grid: AxisymmetricGrid) -> float:
    # Formal L5 radial measure.  The axis has zero weight and does not cause a
    # singular quadrature factor.
    weights = grid.r[:, None] ** 3
    numerator = np.sum(weights * error * error)
    denominator = np.sum(weights * exact * exact)
    return float(np.sqrt(numerator / denominator))


def test_axis_row_has_coefficient_eight_and_outer_identity() -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=12, r_max=2.0)
    matrix = assemble_radial_mode_matrix(grid, q_squared=0.0)
    inv_dr2 = 1.0 / grid.dr**2

    assert matrix[0, 0] == pytest.approx(8.0 * inv_dr2)
    assert matrix[0, 1] == pytest.approx(-8.0 * inv_dr2)
    assert np.count_nonzero(matrix[0, 2:]) == 0

    expected_boundary = np.zeros(grid.nr)
    expected_boundary[-1] = 1.0
    np.testing.assert_allclose(matrix[-1, :].real, expected_boundary)
    np.testing.assert_allclose(matrix[-1, :].imag, 0.0)


def test_fourier_wavenumbers_match_period() -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=16, r_max=1.0, z_max=4.0 * np.pi)
    q = axial_wavenumbers(grid)
    assert q[1] == pytest.approx(0.5)
    assert q[-1] == pytest.approx(-0.5)


def test_nonzero_boundary_manufactured_solution_converges_second_order() -> None:
    sizes = ((17, 32), (33, 64), (65, 128))
    errors: list[float] = []
    cross_defects: list[float] = []

    for nr, nz in sizes:
        grid = AxisymmetricGrid.uniform(nr=nr, nz=nz, r_max=1.7)
        exact, omega, boundary = _manufactured_fields(grid)
        result = solve_finite_cylinder_poisson(
            omega,
            grid,
            outer_boundary=boundary,
            condition_mode_indices=(),
        )
        errors.append(_weighted_relative_l2(result.psi - exact, exact, grid))
        cross_defects.append(result.diagnostics.physical_cross_stencil_defect_max)
        assert result.diagnostics.boundary_error_max < 2.0e-13
        assert result.diagnostics.fourier_matrix_residual_max < 2.0e-9

    orders = [np.log(errors[i] / errors[i + 1]) / np.log(2.0) for i in range(2)]
    assert errors[2] < errors[1] < errors[0]
    assert min(orders) > 1.85
    defect_orders = [np.log(cross_defects[i] / cross_defects[i + 1]) / np.log(2.0) for i in range(2)]
    assert cross_defects[2] < cross_defects[1] < cross_defects[0]
    assert min(defect_orders) > 1.85


def test_global_sign_is_not_silently_reversed() -> None:
    grid = AxisymmetricGrid.uniform(nr=33, nz=32, r_max=1.4)
    exact, omega, boundary = _manufactured_fields(grid, mode=2)

    correct = solve_finite_cylinder_poisson(
        omega,
        grid,
        outer_boundary=boundary,
        condition_mode_indices=(),
    ).psi
    wrong = solve_finite_cylinder_poisson(
        -omega,
        grid,
        outer_boundary=boundary,
        condition_mode_indices=(),
    ).psi

    correct_error = _weighted_relative_l2(correct - exact, exact, grid)
    wrong_error = _weighted_relative_l2(wrong - exact, exact, grid)
    assert correct_error < 2.0e-3
    assert wrong_error > 50.0 * correct_error


def test_physical_space_system_residual_uses_dirichlet_row() -> None:
    grid = AxisymmetricGrid.uniform(nr=25, nz=40, r_max=1.2)
    exact, omega, boundary = _manufactured_fields(grid, mode=4)
    lhs = apply_finite_cylinder_system(exact, grid)

    # The outer equation is psi(R,z)=g(z), not the differential equation.
    np.testing.assert_allclose(lhs[-1, :], boundary, rtol=0.0, atol=0.0)
    assert np.max(np.abs(lhs[:-1, :] - omega[:-1, :])) > 0.0  # continuum/discrete defect


def test_condition_numbers_are_recorded_for_requested_modes() -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=24, r_max=1.0)
    exact, omega, boundary = _manufactured_fields(grid)
    result = solve_finite_cylinder_poisson(
        omega,
        grid,
        outer_boundary=boundary,
        condition_mode_indices=(0, 1, grid.nz // 2),
    )
    assert set(result.diagnostics.condition_numbers) == {0, 1, grid.nz // 2}
    assert all(np.isfinite(value) and value > 1.0 for value in result.diagnostics.condition_numbers.values())


@pytest.mark.parametrize(
    "omega,boundary,match",
    [
        (np.zeros((8, 12)), np.zeros(12), "shape"),
        (np.full((9, 12), np.nan), np.zeros(12), "finite"),
        (np.zeros((9, 12)), np.zeros(11), "shape"),
        (np.zeros((9, 12)), np.full(12, np.inf), "finite"),
    ],
)
def test_invalid_inputs_are_rejected(
    omega: np.ndarray,
    boundary: np.ndarray,
    match: str,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=12, r_max=1.0)
    with pytest.raises(ValueError, match=match):
        solve_finite_cylinder_poisson(omega, grid, outer_boundary=boundary)


def test_nonperiodic_grid_is_rejected() -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=12, r_max=1.0, periodic_z=False)
    with pytest.raises(ValueError, match="periodic"):
        solve_finite_cylinder_poisson(
            np.zeros(grid.shape),
            grid,
            outer_boundary=np.zeros(grid.nz),
        )
