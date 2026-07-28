"""Tests for the derivative tail bounds and the free-space velocity recovery."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ns_certificate_lab.free_space_recovery import (
    GREEN_CONSTANT,
    GREEN_DERIVATIVE_CONSTANTS,
    green_derivative_constant,
    green_quadrature,
    harmonic_interior_factor,
    interior_derivative_tail_bounds,
    multipole_boundary_trace,
    multipole_tail_bounds,
    recover_free_space_velocity,
    source_moments,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.whole_space_gate import (
    BumpSpec,
    ManufacturedFreeSpaceField,
    WholeSpaceEllipticSolver,
)

FIELD = ManufacturedFreeSpaceField(
    bumps=(
        BumpSpec(center=0.5, radius=1.2, amplitude=1.0),
        BumpSpec(center=-0.8, radius=0.9, amplitude=0.6),
    )
)


def _grid(nr: int, nz: int, extent: float = 8.0) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=extent, z_min=-extent, z_max=extent, periodic_z=False
    )


def _moments(nr: int = 385, nz: int = 769, extent: float = 6.0):
    grid = _grid(nr, nz, extent)
    r_mesh, z_mesh = grid.mesh()
    return grid, source_moments(grid, FIELD.omega1(r_mesh, z_mesh))


def _order(coarse: float, fine: float, refinement: float = 2.0) -> float:
    return math.log(coarse / fine) / math.log(refinement)


# --------------------------------------------------------------------------- #
# The exact second derivatives of the manufactured reference                    #
# --------------------------------------------------------------------------- #


def test_exact_second_derivatives_satisfy_the_pde_at_machine_precision():
    """With closed-form second derivatives the PDE residual is roundoff.

    This is the sharpest available confirmation that the Gate-4 reference is
    exact rather than merely accurate: no discretisation is involved.
    """
    r = np.array([0.3, 0.9, 1.6, 2.4, 4.0])
    z = np.array([0.1, -0.5, 0.8, 1.4, -2.0])
    laplacian = (
        FIELD.d2psi1_drr(r, z) + 3.0 / r * FIELD.dpsi1_dr(r, z) + FIELD.d2psi1_dzz(r, z)
    )
    assert float(np.max(np.abs(-laplacian - FIELD.omega1(r, z)))) < 1.0e-14


def test_exact_second_derivatives_match_finite_differences():
    step = 1.0e-5
    r = np.array([0.4, 1.1, 2.5])
    z = np.array([0.2, -0.7, 1.3])
    numeric_rr = (FIELD.dpsi1_dr(r + step, z) - FIELD.dpsi1_dr(r - step, z)) / (2 * step)
    numeric_zz = (FIELD.dpsi1_dz(r, z + step) - FIELD.dpsi1_dz(r, z - step)) / (2 * step)
    numeric_rz = (FIELD.dpsi1_dr(r, z + step) - FIELD.dpsi1_dr(r, z - step)) / (2 * step)
    assert np.allclose(FIELD.d2psi1_drr(r, z), numeric_rr, atol=1.0e-8)
    assert np.allclose(FIELD.d2psi1_dzz(r, z), numeric_zz, atol=1.0e-8)
    assert np.allclose(FIELD.d2psi1_drz(r, z), numeric_rz, atol=1.0e-8)


# --------------------------------------------------------------------------- #
# Main work A: the derivative tail bounds                                      #
# --------------------------------------------------------------------------- #


def test_green_derivative_constants_are_the_documented_multiples():
    expected = (1.0, 3.0, 12.0, 150.0, 1620.0, 21420.0)
    for order, multiple in enumerate(expected):
        assert green_derivative_constant(order) == pytest.approx(
            GREEN_CONSTANT * multiple
        )
    assert len(GREEN_DERIVATIVE_CONSTANTS) == len(expected)
    with pytest.raises(ValueError):
        green_derivative_constant(-1)
    with pytest.raises(ValueError):
        green_derivative_constant(len(expected))


def test_gradient_and_hessian_constants_are_exact_for_the_kernel():
    """``A_1 = 3 c_5`` and ``A_2 = 12 c_5`` are sharp, so check them numerically."""
    step = 1.0e-6

    def kernel(x: np.ndarray) -> np.ndarray:
        return GREEN_CONSTANT * np.linalg.norm(x, axis=-1) ** -3

    direction = np.array([0.3, -0.5, 0.7, 0.2, 0.4])
    direction = direction / np.linalg.norm(direction)
    point = 1.7 * direction
    gradient = np.array(
        [
            (kernel(point + step * np.eye(5)[i]) - kernel(point - step * np.eye(5)[i]))
            / (2 * step)
            for i in range(5)
        ]
    )
    assert float(np.linalg.norm(gradient)) == pytest.approx(
        green_derivative_constant(1) * 1.7**-4, rel=1.0e-6
    )
    # Radial second derivative equals the largest Hessian eigenvalue, 12 c_5 R^-5.
    radial = (
        kernel(point * (1 + step)) - 2 * kernel(point) + kernel(point * (1 - step))
    ) / (step * 1.7) ** 2
    assert radial == pytest.approx(green_derivative_constant(2) * 1.7**-5, rel=1.0e-4)


@pytest.mark.parametrize("truncation", ["monopole", "dipole"])
@pytest.mark.parametrize("radius", [2.5, 4.0])
def test_tail_bounds_dominate_the_measured_error_at_every_order(truncation, radius):
    grid, moments = _moments()
    bound = multipole_tail_bounds(
        moments, evaluation_radius=radius, truncation=truncation
    )
    step = 2.0e-4
    angle = np.linspace(0.0, math.pi, 121)
    rr = radius * np.sin(angle)
    zz = radius * np.cos(angle)

    def approx(a, c):
        return multipole_boundary_trace(grid, moments, a, c, truncation=truncation)

    value = float(np.max(np.abs(FIELD.psi1(rr, zz) - approx(rr, zz))))
    grad_r = (approx(rr + step, zz) - approx(rr - step, zz)) / (2 * step)
    grad_z = (approx(rr, zz + step) - approx(rr, zz - step)) / (2 * step)
    gradient = float(
        np.max(np.hypot(FIELD.dpsi1_dr(rr, zz) - grad_r, FIELD.dpsi1_dz(rr, zz) - grad_z))
    )
    hess = (approx(rr + step, zz) - 2 * approx(rr, zz) + approx(rr - step, zz)) / step**2
    hessian = float(np.max(np.abs(FIELD.d2psi1_drr(rr, zz) - hess)))
    assert 0.0 < value <= bound.value
    assert 0.0 < gradient <= bound.gradient
    assert 0.0 < hessian <= bound.hessian
    assert bound.third > 0.0
    assert bound[0] == bound.value and bound[3] == bound.third


def test_dipole_truncation_is_sharper_than_monopole_far_away():
    _grid_, moments = _moments()
    mono = multipole_tail_bounds(moments, evaluation_radius=6.0, truncation="monopole")
    dip = multipole_tail_bounds(moments, evaluation_radius=6.0, truncation="dipole")
    assert dip.value < mono.value


def test_tail_bounds_scale_with_the_predicted_powers():
    _grid_, moments = _moments()
    # Anchor on the measured support radius so the distances are exactly 10
    # and 20; otherwise the discrete support radius shifts the observed power.
    near = multipole_tail_bounds(
        moments, evaluation_radius=moments.support_radius + 10.0, truncation="monopole"
    )
    far = multipole_tail_bounds(
        moments, evaluation_radius=moments.support_radius + 20.0, truncation="monopole"
    )
    assert _order(near.value, far.value) == pytest.approx(4.0, rel=1.0e-6)
    assert _order(near.gradient, far.gradient) == pytest.approx(5.0, rel=1.0e-6)
    assert _order(near.hessian, far.hessian) == pytest.approx(6.0, rel=1.0e-6)


def test_tail_bounds_validate_their_inputs():
    _grid_, moments = _moments()
    with pytest.raises(ValueError):
        multipole_tail_bounds(moments, evaluation_radius=1.0, truncation="monopole")
    with pytest.raises(ValueError):
        multipole_tail_bounds(moments, evaluation_radius=5.0, truncation="octupole")
    with pytest.raises(ValueError):
        multipole_boundary_trace(
            _grid_, moments, np.array([1.0]), np.array([1.0]), truncation="octupole"
        )


def test_harmonic_interior_factor_is_the_classical_estimate():
    assert harmonic_interior_factor(0, 2.0) == 1.0
    assert harmonic_interior_factor(1, 2.0) == pytest.approx(2.5)
    assert harmonic_interior_factor(2, 2.0) == pytest.approx(25.0)
    with pytest.raises(ValueError):
        harmonic_interior_factor(1, 0.0)
    with pytest.raises(ValueError):
        harmonic_interior_factor(-1, 1.0)


def test_interior_bounds_grow_with_the_derivative_order():
    bounds = interior_derivative_tail_bounds(1.0e-6, distance_to_boundary=2.0)
    assert bounds.value < bounds.gradient < bounds.hessian < bounds.third
    assert bounds.truncation == "box_truncation"
    with pytest.raises(ValueError):
        interior_derivative_tail_bounds(-1.0, distance_to_boundary=1.0)


def test_source_moments_report_the_neglected_mass_of_a_threshold():
    grid = _grid(257, 513, 6.0)
    r_mesh, z_mesh = grid.mesh()
    omega = FIELD.omega1(r_mesh, z_mesh)
    exact = source_moments(grid, omega)
    assert exact.neglected_fraction == 0.0
    assert exact.support_radius == pytest.approx(FIELD.support_radius, rel=5.0e-3)
    truncated = source_moments(grid, omega, support_threshold=1.0e-2)
    assert truncated.support_radius < exact.support_radius
    assert 0.0 < truncated.neglected_fraction < 1.0e-2
    assert exact.newtonian_mass == pytest.approx(FIELD.mass, rel=1.0e-3)


# --------------------------------------------------------------------------- #
# Main work B: the two independent evaluation paths                            #
# --------------------------------------------------------------------------- #


def test_analytic_green_path_matches_the_exact_solution_and_converges():
    """Kernel differentiated analytically, source quadratured; no solver used."""
    targets_r = np.array([2.2, 2.6, 3.0, 0.1, 2.0])
    targets_z = np.array([0.3, -1.0, 1.5, 2.4, 2.2])
    exact = {
        "psi": FIELD.psi1(targets_r, targets_z),
        "dr": FIELD.dpsi1_dr(targets_r, targets_z),
        "dz": FIELD.dpsi1_dz(targets_r, targets_z),
        "drr": FIELD.d2psi1_drr(targets_r, targets_z),
        "dzz": FIELD.d2psi1_dzz(targets_r, targets_z),
        "drz": FIELD.d2psi1_drz(targets_r, targets_z),
    }
    errors = []
    for nr in (81, 161):
        grid = AxisymmetricGrid.uniform(
            nr=nr, nz=2 * nr - 1, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=False
        )
        r_mesh, z_mesh = grid.mesh()
        out = green_quadrature(
            grid, FIELD.omega1(r_mesh, z_mesh), targets_r, targets_z,
            quadrature_nodes=32,
        )
        errors.append(
            {
                key: float(np.max(np.abs(out[key] - reference)))
                / max(float(np.max(np.abs(reference))), 1.0e-300)
                for key, reference in exact.items()
            }
        )
    for key in exact:
        assert errors[1][key] < errors[0][key]
        assert _order(errors[0][key], errors[1][key]) > 2.5
        assert errors[1][key] < 1.0e-5


def test_green_quadrature_validates_its_inputs():
    grid = _grid(33, 65, 4.0)
    r_mesh, z_mesh = grid.mesh()
    omega = FIELD.omega1(r_mesh, z_mesh)
    with pytest.raises(ValueError):
        green_quadrature(grid, omega, np.array([1.0, 2.0]), np.array([1.0]))
    with pytest.raises(ValueError):
        green_quadrature(
            grid, omega, np.array([3.0]), np.array([0.0]), quadrature_nodes=1
        )


@pytest.mark.parametrize("nr", [65, 129])
def test_recovery_reproduces_the_exact_velocity(nr):
    grid = _grid(nr, 2 * nr - 1, 8.0)
    r_mesh, z_mesh = grid.mesh()
    solver = WholeSpaceEllipticSolver.build(grid, boundary_mode="monopole")
    recovery = recover_free_space_velocity(solver, FIELD.omega1(r_mesh, z_mesh))
    mask = np.hypot(r_mesh, z_mesh) <= 3.0
    exact_r = -r_mesh * FIELD.dpsi1_dz(r_mesh, z_mesh)
    exact_z = 2.0 * FIELD.psi1(r_mesh, z_mesh) + r_mesh * FIELD.dpsi1_dr(r_mesh, z_mesh)
    scale = float(np.max(np.abs(exact_z)))
    assert float(np.max(np.abs((recovery.u_r - exact_r) * mask))) < 0.05 * scale
    assert float(np.max(np.abs((recovery.u_z - exact_z) * mask))) < 0.05 * scale
    assert recovery.poisson_residual_max < 1.0e-9


def test_recovery_is_second_order_in_every_returned_field():
    errors = []
    for nr in (65, 129, 257):
        grid = _grid(nr, 2 * nr - 1, 8.0)
        r_mesh, z_mesh = grid.mesh()
        solver = WholeSpaceEllipticSolver.build(grid, boundary_mode="monopole")
        recovery = recover_free_space_velocity(solver, FIELD.omega1(r_mesh, z_mesh))
        mask = np.hypot(r_mesh, z_mesh) <= 3.0
        errors.append(
            {
                "psi": float(np.max(np.abs((recovery.psi1 - FIELD.psi1(r_mesh, z_mesh)) * mask))),
                "dr": float(np.max(np.abs((recovery.dpsi1_dr - FIELD.dpsi1_dr(r_mesh, z_mesh)) * mask))),
                "dz": float(np.max(np.abs((recovery.dpsi1_dz - FIELD.dpsi1_dz(r_mesh, z_mesh)) * mask))),
            }
        )
    for key in ("psi", "dr", "dz"):
        orders = [
            _order(errors[i][key], errors[i + 1][key]) for i in range(len(errors) - 1)
        ]
        assert min(orders) > 1.8


def test_recovery_enforces_axis_regularity_exactly():
    grid = _grid(65, 129, 8.0)
    r_mesh, z_mesh = grid.mesh()
    solver = WholeSpaceEllipticSolver.build(grid)
    recovery = recover_free_space_velocity(solver, FIELD.omega1(r_mesh, z_mesh))
    assert np.all(recovery.u_r[0] == 0.0)
    assert np.all(recovery.dpsi1_dr[0] == 0.0)


def test_recovery_reports_a_velocity_tail_bound_and_validates_the_interior_radius():
    grid = _grid(65, 129, 8.0)
    r_mesh, z_mesh = grid.mesh()
    solver = WholeSpaceEllipticSolver.build(grid)
    recovery = recover_free_space_velocity(solver, FIELD.omega1(r_mesh, z_mesh))
    assert recovery.velocity_tail_bound > recovery.interior_tail.value
    payload = recovery.as_dict()
    assert payload["boundary_tail"]["truncation"] == "monopole"
    assert payload["interior_tail"]["truncation"] == "box_truncation"
    with pytest.raises(ValueError):
        recover_free_space_velocity(
            solver, FIELD.omega1(r_mesh, z_mesh), interior_radius=99.0
        )


# --------------------------------------------------------------------------- #
# Fault injection                                                              #
# --------------------------------------------------------------------------- #


def test_fault_injection_wrong_axis_coefficient_inflates_the_gradient_error():
    grid = _grid(129, 257, 8.0)
    r_mesh, z_mesh = grid.mesh()
    omega = FIELD.omega1(r_mesh, z_mesh)
    exact = FIELD.dpsi1_dr(r_mesh, z_mesh)
    mask = np.hypot(r_mesh, z_mesh) <= 3.0

    def error(**kwargs) -> float:
        solver = WholeSpaceEllipticSolver.build(grid, boundary_mode="monopole", **kwargs)
        recovery = recover_free_space_velocity(solver, omega)
        return float(np.max(np.abs((recovery.dpsi1_dr - exact) * mask)))

    assert error(axis_coefficient=4.0) > 2.0 * error()
    assert error(radial_drift_coefficient=1.0) > 10.0 * error()


def test_fault_injection_understated_moment_breaks_the_tail_bound():
    grid, moments = _moments()
    honest = multipole_tail_bounds(moments, evaluation_radius=2.5, truncation="monopole")
    understated = honest.value / 1000.0
    angle = np.linspace(0.0, math.pi, 61)
    rr = 2.5 * np.sin(angle)
    zz = 2.5 * np.cos(angle)
    measured = float(
        np.max(
            np.abs(
                FIELD.psi1(rr, zz)
                - multipole_boundary_trace(grid, moments, rr, zz, truncation="monopole")
            )
        )
    )
    assert measured <= honest.value
    assert measured > understated
