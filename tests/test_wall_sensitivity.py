from __future__ import annotations

import math

import numpy as np

from ns_certificate_lab.finite_cylinder_poisson import solve_finite_cylinder_poisson
from ns_certificate_lab.free_space_poisson import solve_l5_free_radial_periodized_z
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.wall_sensitivity import (
    predicted_free_boundary_from_source,
    predicted_zero_dirichlet_correction,
    regular_mode_ratio,
)


def _relative_l2(error: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(error.ravel()) / np.linalg.norm(reference.ravel()))


def test_zero_mode_boundary_value_equals_radial_moment_formula() -> None:
    r = np.linspace(0.0, 4.0, 257)
    source_radial = np.exp(-r**2)
    source = np.repeat(source_radial[:, None], 16, axis=1)
    free, _ = solve_l5_free_radial_periodized_z(source, r, dz=0.5)
    predicted = predicted_free_boundary_from_source(source_radial, r, 0.0)
    np.testing.assert_allclose(
        free[-1, :], predicted, rtol=2.0e-12, atol=2.0e-12
    )



def test_nonzero_mode_boundary_value_matches_source_tail_formula() -> None:
    period = 10.0
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
    source_radial = np.exp(-grid.r**2)
    source = source_radial[:, None] * np.cos(k * zz)
    free, _ = solve_l5_free_radial_periodized_z(source, grid.r, grid.dz)
    predicted = predicted_free_boundary_from_source(source_radial, grid.r, k)
    np.testing.assert_allclose(free[-1, 0], predicted, rtol=4.0e-7, atol=2.0e-12)

def test_zero_mode_wall_error_is_asymptotically_constant() -> None:
    grid = AxisymmetricGrid.uniform(
        nr=257,
        nz=16,
        r_max=4.0,
        z_min=0.0,
        z_max=8.0,
        periodic_z=True,
    )
    rr, _ = grid.mesh()
    source = np.exp(-rr**2)
    free, _ = solve_l5_free_radial_periodized_z(source, grid.r, grid.dz)
    finite = solve_finite_cylinder_poisson(
        source,
        grid,
        outer_boundary=0.0,
        condition_mode_indices=(),
    ).psi
    predicted = predicted_zero_dirichlet_correction(
        grid.r,
        0.0,
        float(free[-1, 0]),
    )[:, None]
    interior = grid.r <= 2.0
    defect = _relative_l2((finite - free - predicted)[interior], (finite - free)[interior])
    assert defect < 5.0e-4


def test_low_nonzero_mode_wall_correction_shape() -> None:
    period = 16.0
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
    radial_correction = predicted_zero_dirichlet_correction(
        grid.r,
        k,
        float(free[-1, 0]),
    )
    predicted = radial_correction[:, None] * np.cos(k * grid.z)[None, :]
    interior = grid.r <= 2.0
    defect = _relative_l2((finite - free - predicted)[interior], (finite - free)[interior])
    assert defect < 2.0e-3


def test_wall_sensitivity_grows_as_lowest_periodic_wavenumber_falls() -> None:
    errors = []
    for period in (2.0, 4.0, 8.0, 16.0):
        grid = AxisymmetricGrid.uniform(
            nr=129,
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
        interior = grid.r <= 2.0
        errors.append(_relative_l2((finite - free)[interior], free[interior]))

    # The L=4 value can be slightly below L=2 because finite-difference error
    # dominates when the true wall correction is extremely small.  Once the
    # low-frequency tail is resolved, sensitivity rises decisively.
    assert errors[2] > 5.0 * errors[1]
    assert errors[3] > 4.0 * errors[2]


def test_regular_mode_ratio_is_normalized_and_axis_finite() -> None:
    r = np.linspace(0.0, 4.0, 65)
    ratio = regular_mode_ratio(r, 0.4)
    assert np.isfinite(ratio[0])
    np.testing.assert_allclose(ratio[-1], 1.0, rtol=1e-12, atol=1e-12)
    assert np.all(np.diff(ratio) > 0.0)
