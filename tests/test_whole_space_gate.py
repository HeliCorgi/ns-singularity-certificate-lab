"""Tests for Gate 4: the linear whole-space elliptic gate with non-periodic z.

The manufactured reference is exact and closed-form, so the convergence tests
here measure a genuine discretisation error and not a self-consistency residual.
Every fault-injection test asserts that the injected defect is detected, as
AGENTS.md requires.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.whole_space_gate import (
    BumpSpec,
    ManufacturedFreeSpaceField,
    WholeSpaceEllipticSolver,
    audit_cartesian_recovery,
    axial_spectrum,
    domain_truncation_series,
    monopole_tail_bound,
    separate_period_and_wall_errors,
    solve_and_compare,
    source_l1_norm_5d,
    source_monopole_moment,
    zero_frequency_stress,
)

MONOPOLE = ManufacturedFreeSpaceField(
    bumps=(
        BumpSpec(center=0.5, radius=1.2, amplitude=1.0),
        BumpSpec(center=-0.8, radius=0.9, amplitude=0.6),
    )
)
#: A balanced pair with zero net five-dimensional mass: the far field decays one
#: order faster, which is the sharpest probe of the truncation error.
DIPOLE = ManufacturedFreeSpaceField(
    bumps=(
        BumpSpec(center=0.9, radius=1.2, amplitude=1.0),
        BumpSpec(center=-0.9, radius=0.9, amplitude=-((1.2 / 0.9) ** 5)),
    )
)


def _observed_order(coarse: float, fine: float, refinement: float = 2.0) -> float:
    return math.log(coarse / fine) / math.log(refinement)


# --------------------------------------------------------------------------- #
# Item 1/2: the manufactured source and its exact free-space solution           #
# --------------------------------------------------------------------------- #


def test_source_is_compactly_supported_in_both_directions():
    field = MONOPOLE
    outside_r = field.omega1(np.array([3.0, 5.0]), np.array([0.0, 0.0]))
    outside_z = field.omega1(np.array([0.0, 0.0]), np.array([4.0, -4.0]))
    assert np.all(outside_r == 0.0)
    assert np.all(outside_z == 0.0)
    assert field.support_radius == pytest.approx(1.7)


def test_reference_solves_the_pde_and_the_residual_is_second_order():
    """The closed-form reference is exact: its PDE residual is pure truncation."""
    residuals = []
    for nr in (129, 257):
        grid = AxisymmetricGrid.uniform(
            nr=nr, nz=2 * nr - 1, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=False
        )
        r_mesh, z_mesh = grid.mesh()
        psi = MONOPOLE.psi1(r_mesh, z_mesh)
        omega = MONOPOLE.omega1(r_mesh, z_mesh)
        laplacian = np.zeros_like(psi)
        laplacian[1:-1, :] = (
            psi[2:, :] - 2.0 * psi[1:-1, :] + psi[:-2, :]
        ) / grid.dr**2 + 3.0 / grid.r[1:-1, None] * (psi[2:, :] - psi[:-2, :]) / (
            2.0 * grid.dr
        )
        laplacian[1:-1, 1:-1] += (
            psi[1:-1, 2:] - 2.0 * psi[1:-1, 1:-1] + psi[1:-1, :-2]
        ) / grid.dz**2
        residuals.append(
            float(np.max(np.abs(-laplacian[1:-1, 1:-1] - omega[1:-1, 1:-1])))
        )
    assert residuals[1] < residuals[0]
    assert _observed_order(*residuals) > 1.85


def test_reference_derivatives_match_finite_differences():
    step = 1.0e-6
    radius = np.array([0.4, 1.1, 2.5])
    axial = np.array([0.2, -0.7, 1.3])
    numeric_r = (
        MONOPOLE.psi1(radius + step, axial) - MONOPOLE.psi1(radius - step, axial)
    ) / (2.0 * step)
    numeric_z = (
        MONOPOLE.psi1(radius, axial + step) - MONOPOLE.psi1(radius, axial - step)
    ) / (2.0 * step)
    assert np.allclose(MONOPOLE.dpsi1_dr(radius, axial), numeric_r, atol=1.0e-9)
    assert np.allclose(MONOPOLE.dpsi1_dz(radius, axial), numeric_z, atol=1.0e-9)


def test_far_field_approaches_the_monopole_and_the_dipole_cancels_it():
    """``3R^3 psi -> M``.  The approach is only asymptotic: the exact field is
    a superposition of *shifted* bumps, so the monopole about the origin is a
    limit and not an identity at finite radius."""
    assert DIPOLE.mass == pytest.approx(0.0, abs=1.0e-15)
    ratios = [
        float(MONOPOLE.psi1(np.array([0.0]), np.array([radius]))[0])
        * 3.0
        * radius**3
        / MONOPOLE.mass
        for radius in (50.0, 200.0, 1000.0)
    ]
    assert ratios[0] > ratios[1] > ratios[2] > 1.0
    assert ratios[-1] == pytest.approx(1.0, abs=2.0e-3)
    # The balanced pair has no monopole at all, so its field decays one order
    # faster: the ratio to a monopole field falls off like 1/R.
    ratio = [
        abs(float(DIPOLE.psi1(np.array([0.0]), np.array([radius]))[0]))
        / abs(float(MONOPOLE.psi1(np.array([0.0]), np.array([radius]))[0]))
        for radius in (100.0, 1000.0)
    ]
    assert ratio[0] > ratio[1]
    assert ratio[1] < 0.01


def test_bump_rejects_invalid_specifications():
    with pytest.raises(ValueError):
        BumpSpec(center=0.0, radius=-1.0, amplitude=1.0)
    with pytest.raises(ValueError):
        BumpSpec(center=0.0, radius=1.0, amplitude=1.0, power=1)
    with pytest.raises(ValueError):
        ManufacturedFreeSpaceField(bumps=())


def test_source_moments_from_grid_data_match_the_closed_form():
    grid = AxisymmetricGrid.uniform(
        nr=513, nz=1025, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=False
    )
    r_mesh, z_mesh = grid.mesh()
    omega = MONOPOLE.omega1(r_mesh, z_mesh)
    assert source_monopole_moment(grid, omega) == pytest.approx(
        MONOPOLE.mass, rel=1.0e-4
    )
    assert source_l1_norm_5d(grid, omega) == pytest.approx(
        MONOPOLE.l1_norm_5d, rel=1.0e-4
    )


# --------------------------------------------------------------------------- #
# The solver                                                                    #
# --------------------------------------------------------------------------- #


def test_solver_rejects_a_periodic_grid_and_bad_modes():
    periodic = AxisymmetricGrid.uniform(
        nr=17, nz=32, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=True
    )
    with pytest.raises(ValueError):
        WholeSpaceEllipticSolver.build(periodic)
    aperiodic = AxisymmetricGrid.uniform(
        nr=17, nz=33, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=False
    )
    with pytest.raises(ValueError):
        WholeSpaceEllipticSolver.build(aperiodic, boundary_mode="reflecting")


def test_discrete_residual_is_at_roundoff():
    grid = AxisymmetricGrid.uniform(
        nr=97, nz=193, r_max=8.0, z_min=-8.0, z_max=8.0, periodic_z=False
    )
    solver = WholeSpaceEllipticSolver.build(grid)
    r_mesh, z_mesh = grid.mesh()
    omega = MONOPOLE.omega1(r_mesh, z_mesh)
    psi = solver.solve(omega)
    scale = float(np.max(np.abs(omega)))
    assert float(np.max(np.abs(solver.residual(psi, omega)))) < 1.0e-9 * scale


def test_solve_is_pure_and_linear():
    grid = AxisymmetricGrid.uniform(
        nr=65, nz=129, r_max=8.0, z_min=-8.0, z_max=8.0, periodic_z=False
    )
    solver = WholeSpaceEllipticSolver.build(grid)
    r_mesh, z_mesh = grid.mesh()
    omega = MONOPOLE.omega1(r_mesh, z_mesh)
    first = solver.solve(omega)
    assert np.array_equal(first, solver.solve(omega))
    assert np.allclose(solver.solve(-omega), -first, rtol=0.0, atol=1.0e-15)


def _refinement_errors(boundary_mode, extent=8.0):
    return [
        solve_and_compare(
            MONOPOLE,
            r_max=extent,
            z_max=extent,
            nr=nr,
            nz=2 * nr - 1,
            boundary_mode=boundary_mode,
            interior_radius=3.0,
        ).interior_relative_l2
        for nr in (65, 129, 257)
    ]


def test_grid_refinement_is_second_order_with_the_monopole_trace():
    errors = _refinement_errors("monopole")
    orders = [_observed_order(errors[i], errors[i + 1]) for i in range(len(errors) - 1)]
    assert all(error > 0.0 for error in errors)
    assert min(orders) > 1.85


def test_zero_trace_saturates_and_that_is_why_the_monopole_trace_exists():
    """A zero outer trace stops converging: truncation dominates immediately.

    This is a *finding*, not a defect.  On the same box the monopole trace
    keeps second order, which is the quantitative argument for using it.
    """
    zero_errors = _refinement_errors("zero")
    monopole_errors = _refinement_errors("monopole")
    zero_orders = [
        _observed_order(zero_errors[i], zero_errors[i + 1])
        for i in range(len(zero_errors) - 1)
    ]
    assert min(zero_orders) < 1.0
    assert zero_errors[-1] > 5.0 * monopole_errors[-1]


def test_axial_spectrum_has_no_zero_mode_and_matches_the_continuum():
    for z_max in (4.0, 8.0):
        grid = AxisymmetricGrid.uniform(
            nr=33,
            nz=2 * int(z_max * 8) + 1,
            r_max=8.0,
            z_min=-z_max,
            z_max=z_max,
            periodic_z=False,
        )
        report = axial_spectrum(grid)
        assert not report.zero_mode_present
        assert report.smallest_magnitude > 0.0
        assert report.smallest_magnitude == pytest.approx(
            (math.pi / (2.0 * z_max)) ** 2, rel=1.0e-3
        )


def test_axial_spectrum_rejects_a_periodic_grid():
    periodic = AxisymmetricGrid.uniform(
        nr=17, nz=32, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=True
    )
    with pytest.raises(ValueError):
        axial_spectrum(periodic)


# --------------------------------------------------------------------------- #
# Item 3: independent domain enlargement                                        #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("boundary_mode", "minimum_rate"),
    [("zero", 2.7), ("monopole", 3.5)],
)
def test_radial_enlargement_reaches_the_predicted_truncation_rate(
    boundary_mode, minimum_rate
):
    """Zero trace truncates at ``R^-3``; the monopole trace at ``R^-4``."""
    rows = domain_truncation_series(
        MONOPOLE,
        sizes=(5.0, 7.5, 11.25),
        fixed_extent=20.0,
        spacing=0.125 if boundary_mode == "monopole" else 0.25,
        enlarge="radial",
        boundary_mode=boundary_mode,
        interior_radius=3.0,
        reference_size=45.0,
    )
    values = [row.common_interior_max_diff for row in rows]
    rates = [
        _observed_order(values[i], values[i + 1], refinement=1.5)
        for i in range(len(values) - 1)
    ]
    assert min(rates) > minimum_rate


def test_axial_enlargement_is_independent_of_the_radial_extent():
    rows = domain_truncation_series(
        MONOPOLE,
        sizes=(5.0, 7.5, 11.25),
        fixed_extent=20.0,
        spacing=0.25,
        enlarge="axial",
        boundary_mode="monopole",
        interior_radius=3.0,
        reference_size=45.0,
    )
    values = [row.common_interior_max_diff for row in rows]
    assert values[0] > values[1] > values[2]
    assert _observed_order(values[0], values[1], refinement=1.5) > 3.5


def test_truncation_series_validates_its_inputs():
    with pytest.raises(ValueError):
        domain_truncation_series(
            MONOPOLE,
            sizes=(5.0,),
            fixed_extent=20.0,
            spacing=0.5,
            enlarge="radial",
            boundary_mode="zero",
            interior_radius=3.0,
        )
    with pytest.raises(ValueError):
        domain_truncation_series(
            MONOPOLE,
            sizes=(5.0, 7.5),
            fixed_extent=20.0,
            spacing=0.5,
            enlarge="diagonal",
            boundary_mode="zero",
            interior_radius=3.0,
        )
    with pytest.raises(ValueError):
        domain_truncation_series(
            MONOPOLE,
            sizes=(5.0, 7.5),
            fixed_extent=20.0,
            spacing=0.5,
            enlarge="radial",
            boundary_mode="zero",
            interior_radius=9.0,
        )


# --------------------------------------------------------------------------- #
# Item 4: low- and zero-frequency stress                                        #
# --------------------------------------------------------------------------- #


def test_periodic_zero_mode_overstates_the_far_field_linearly_in_radius():
    """The exact ratio is ``2R/L``: the periodic artefact grows with the box."""
    period = 12.0
    rows = zero_frequency_stress(
        MONOPOLE, radial_extents=(4.0, 8.0, 16.0, 32.0), axial_period=period
    )
    for row in rows:
        assert row.ratio == pytest.approx(2.0 * row.radial_extent / period, rel=1.0e-9)
    assert rows[-1].ratio > rows[0].ratio
    assert rows[-1].ratio > 4.0


def test_zero_frequency_stress_validates_its_inputs():
    with pytest.raises(ValueError):
        zero_frequency_stress(MONOPOLE, radial_extents=(4.0,), axial_period=0.0)
    with pytest.raises(ValueError):
        zero_frequency_stress(MONOPOLE, radial_extents=(1.0,), axial_period=8.0)


# --------------------------------------------------------------------------- #
# Item 5: period-image error separated from the radial-wall error               #
# --------------------------------------------------------------------------- #


def test_period_image_component_falls_when_the_period_grows():
    reports = [
        separate_period_and_wall_errors(
            MONOPOLE,
            radial_extent=8.0,
            axial_half_period=half,
            spacing=0.5,
            interior_radius=2.5,
        )
        for half in (6.0, 12.0)
    ]
    assert reports[0].period_image_exact > reports[1].period_image_exact
    for report in reports:
        # The truncated image sum is converged far below the value it reports.
        assert report.period_image_tail_bound < 1.0e-3 * report.period_image_exact
        assert report.radial_wall_bound > 0.0
        assert report.nonperiodic_truncation_bound > 0.0
    # Documented gate limitation: at any affordable spacing the measured errors
    # are dominated by discretisation, which is why the separation is exact.
    assert reports[0].measured_finite_cylinder > reports[0].period_image_exact


def test_separation_rejects_an_interior_region_outside_the_box():
    with pytest.raises(ValueError):
        separate_period_and_wall_errors(
            MONOPOLE,
            radial_extent=8.0,
            axial_half_period=3.0,
            spacing=0.5,
            interior_radius=5.0,
        )


# --------------------------------------------------------------------------- #
# Item 6/7: Cartesian recovery and independent divergence/curl audit            #
# --------------------------------------------------------------------------- #


def test_cartesian_divergence_and_curl_converge_at_second_order():
    reports = [
        audit_cartesian_recovery(MONOPOLE, half_width=3.0, points=points)
        for points in (33, 65, 129)
    ]
    divergence = [report.divergence_max for report in reports]
    curl = [report.curl_theta_error_max for report in reports]
    assert _observed_order(divergence[1], divergence[2]) > 1.8
    assert _observed_order(curl[1], curl[2]) > 1.8
    assert reports[-1].curl_theta_relative < 1.0e-2
    assert reports[-1].curl_axial_max < 1.0e-3 * reports[-1].velocity_max


def test_cartesian_audit_validates_its_inputs():
    with pytest.raises(ValueError):
        audit_cartesian_recovery(MONOPOLE, half_width=3.0, points=4)
    with pytest.raises(ValueError):
        audit_cartesian_recovery(MONOPOLE, half_width=3.0, points=32, interior_trim=0)


# --------------------------------------------------------------------------- #
# Item 8: the a posteriori tail bound                                           #
# --------------------------------------------------------------------------- #


def test_tail_bound_dominates_the_measured_truncation_for_the_monopole_trace():
    rows = domain_truncation_series(
        MONOPOLE,
        sizes=(5.0, 7.5, 11.25),
        fixed_extent=20.0,
        spacing=0.25,
        enlarge="radial",
        boundary_mode="monopole",
        interior_radius=3.0,
        reference_size=45.0,
    )
    for row in rows:
        assert row.boundary_data_error_max <= row.tail_bound
        assert row.common_interior_max_diff <= row.tail_bound


def test_tail_bound_requires_the_boundary_outside_the_support():
    with pytest.raises(ValueError):
        monopole_tail_bound(
            source_l1_norm_5d=1.0, support_radius=2.0, boundary_min_radius=2.0
        )
    with pytest.raises(ValueError):
        monopole_tail_bound(
            source_l1_norm_5d=-1.0, support_radius=1.0, boundary_min_radius=4.0
        )


def test_tail_bound_scales_like_the_fourth_inverse_power():
    first = monopole_tail_bound(
        source_l1_norm_5d=1.0, support_radius=1.0, boundary_min_radius=11.0
    )
    second = monopole_tail_bound(
        source_l1_norm_5d=1.0, support_radius=1.0, boundary_min_radius=21.0
    )
    assert _observed_order(first, second) == pytest.approx(4.0, rel=1.0e-9)


# --------------------------------------------------------------------------- #
# Item 9: fault injection                                                       #
# --------------------------------------------------------------------------- #


def _axis_errors(**solver_kwargs):
    """Max error on the axis row itself.

    The weighted interior norm carries the factor ``r^3`` and therefore cannot
    see an axis-row defect at all; the axis has to be probed directly.
    """
    errors = []
    for nr in (65, 129, 257):
        grid = AxisymmetricGrid.uniform(
            nr=nr, nz=2 * nr - 1, r_max=8.0, z_min=-8.0, z_max=8.0, periodic_z=False
        )
        solver = WholeSpaceEllipticSolver.build(
            grid, boundary_mode="monopole", **solver_kwargs
        )
        r_mesh, z_mesh = grid.mesh()
        error = solver.solve(MONOPOLE.omega1(r_mesh, z_mesh)) - MONOPOLE.psi1(
            r_mesh, z_mesh
        )
        errors.append(float(np.max(np.abs(error[0, np.abs(grid.z) <= 3.0]))))
    return errors


def test_fault_injection_wrong_axis_coefficient_is_detected():
    """The audited axis row is ``8(psi1-psi0)/dr^2``.

    Recorded limitation: replacing ``8`` by ``4`` does **not** destroy the
    observed order -- the defective row occupies a control volume of measure
    ``O(dr^4)`` in the five-dimensional weighting, so its contribution still
    converges.  It inflates the axis error by a stable factor instead, so a
    convergence-order gate alone would miss this fault and the amplitude has
    to be compared as well.
    """
    healthy = _axis_errors()
    faulty = _axis_errors(axis_coefficient=4.0)
    healthy_orders = [
        _observed_order(healthy[i], healthy[i + 1]) for i in range(len(healthy) - 1)
    ]
    assert min(healthy_orders) > 1.85
    assert all(bad > 2.0 * good for good, bad in zip(healthy, faulty))


def test_fault_injection_wrong_radial_drift_coefficient_is_detected():
    """The operator is ``d_rr + 3/r d_r + d_zz``; ``1/r`` must break convergence."""
    faulty = _axis_errors(radial_drift_coefficient=1.0)
    orders = [_observed_order(faulty[i], faulty[i + 1]) for i in range(len(faulty) - 1)]
    assert max(orders) < 0.5


def test_fault_injection_dropping_the_outer_coupling_is_detected():
    """Forgetting to move the ``r=R_max`` Dirichlet trace into the RHS."""
    grid = AxisymmetricGrid.uniform(
        nr=161, nz=769, r_max=2.5, z_min=-6.0, z_max=6.0, periodic_z=False
    )
    r_mesh, z_mesh = grid.mesh()
    omega = MONOPOLE.omega1(r_mesh, z_mesh)
    exact = MONOPOLE.psi1(r_mesh, z_mesh)
    mask = np.hypot(r_mesh, z_mesh) <= 2.0
    healthy = WholeSpaceEllipticSolver.build(grid, boundary_mode="monopole")
    faulty = dataclasses.replace(healthy, outer_coupling=0.0)
    healthy_error = float(np.max(np.abs((healthy.solve(omega) - exact) * mask)))
    faulty_error = float(np.max(np.abs((faulty.solve(omega) - exact) * mask)))
    assert faulty_error > 5.0 * healthy_error


def test_fault_injection_understated_source_norm_breaks_the_tail_bound():
    rows = domain_truncation_series(
        MONOPOLE,
        sizes=(5.0, 7.5),
        fixed_extent=20.0,
        spacing=0.5,
        enlarge="radial",
        boundary_mode="monopole",
        interior_radius=3.0,
        reference_size=30.0,
    )
    understated = monopole_tail_bound(
        source_l1_norm_5d=MONOPOLE.l1_norm_5d / 1000.0,
        support_radius=MONOPOLE.support_radius,
        boundary_min_radius=5.0,
    )
    assert rows[0].boundary_data_error_max <= rows[0].tail_bound
    assert rows[0].boundary_data_error_max > understated


def test_fault_injection_collapsed_pivot_is_reported():
    grid = AxisymmetricGrid.uniform(
        nr=33, nz=65, r_max=8.0, z_min=-8.0, z_max=8.0, periodic_z=False
    )
    solver = WholeSpaceEllipticSolver.build(grid, pivot_floor=1.0e6)
    r_mesh, z_mesh = grid.mesh()
    with pytest.raises(ZeroDivisionError):
        solver.solve(MONOPOLE.omega1(r_mesh, z_mesh))


def test_fault_injection_sign_flipped_axial_derivative_breaks_the_curl_identity():
    class SignFlipped:
        omega1 = staticmethod(MONOPOLE.omega1)
        psi1 = staticmethod(MONOPOLE.psi1)
        dpsi1_dr = staticmethod(MONOPOLE.dpsi1_dr)

        @staticmethod
        def dpsi1_dz(r, z):
            return -MONOPOLE.dpsi1_dz(r, z)

    healthy = audit_cartesian_recovery(MONOPOLE, half_width=3.0, points=65)
    faulty = audit_cartesian_recovery(SignFlipped(), half_width=3.0, points=65)
    assert faulty.curl_theta_relative > 0.2
    assert faulty.curl_theta_relative > 10.0 * healthy.curl_theta_relative
    assert faulty.divergence_max > 100.0 * healthy.divergence_max


# --------------------------------------------------------------------------- #
# Item 10: the interface the nonlinear evolution will consume                   #
# --------------------------------------------------------------------------- #


def test_velocity_recovery_matches_the_analytic_profile():
    grid = AxisymmetricGrid.uniform(
        nr=257, nz=513, r_max=8.0, z_min=-8.0, z_max=8.0, periodic_z=False
    )
    solver = WholeSpaceEllipticSolver.build(grid)
    r_mesh, z_mesh = grid.mesh()
    psi = solver.solve(MONOPOLE.omega1(r_mesh, z_mesh))
    u_r, u_z = solver.recover_meridional_velocity(psi)
    exact_r = -r_mesh * MONOPOLE.dpsi1_dz(r_mesh, z_mesh)
    exact_z = 2.0 * MONOPOLE.psi1(r_mesh, z_mesh) + r_mesh * MONOPOLE.dpsi1_dr(
        r_mesh, z_mesh
    )
    interior = np.hypot(r_mesh, z_mesh) <= 3.0
    scale = float(np.max(np.abs(exact_z)))
    assert float(np.max(np.abs((u_r - exact_r) * interior))) < 5.0e-3 * scale
    assert float(np.max(np.abs((u_z - exact_z) * interior))) < 5.0e-3 * scale
    assert np.all(u_r[0] == 0.0)


def test_coupling_contract_states_what_is_not_established():
    contract = WholeSpaceEllipticSolver.coupling_contract()
    assert contract["state_mutated_by_solve"] is False
    assert any("nonlinear" in item for item in contract["not_established"])
    assert any("interval" in item for item in contract["not_established"])
