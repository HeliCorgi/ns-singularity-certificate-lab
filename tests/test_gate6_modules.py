"""Tests for the Gate-6 modules: initial data, tail propagation, expansion, intervals."""

from __future__ import annotations

from fractions import Fraction
import json
import math

import numpy as np
import pytest

from ns_certificate_lab.cartesian_validation import (
    UniformCartesianGrid,
    cartesian_divergence,
)
from ns_certificate_lab.domain_expansion import (
    ExpansionPolicy,
    embed_in_larger_box,
    expansion_is_required,
)
from ns_certificate_lab.free_space_recovery import (
    multipole_boundary_trace,
    multipole_tail_bounds,
    recover_free_space_velocity,
    source_moments,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.initial_data import (
    SwirlFamily,
    smooth_bump,
    smooth_bump_derivative,
)
from ns_certificate_lab.snapshot_certificate import (
    Interval,
    build_snapshot_certificate,
    verify_snapshot_certificate,
)
from ns_certificate_lab.tail_propagation import (
    FieldBounds,
    PotentialErrorBounds,
    advection_defect_bound,
    gronwall_state_error,
    velocity_error_bounds,
)
from ns_certificate_lab.whole_space_evolution import WholeSpaceEvolution
from ns_certificate_lab.whole_space_gate import (
    BumpSpec,
    ManufacturedFreeSpaceField,
    WholeSpaceEllipticSolver,
)

FAMILY = SwirlFamily(
    amplitude=10.0, radial_support=1.2, axial_support=1.5, axial_concentration=0.5
)


# --------------------------------------------------------------------------- #
# The explicit initial-data family                                             #
# --------------------------------------------------------------------------- #


def test_bump_is_smooth_and_compactly_supported():
    s = np.array([-0.5, 0.0, 0.25, 0.9, 0.99, 1.0, 1.5])
    values = smooth_bump(s)
    assert values[0] == 0.0  # negative argument is outside the support
    assert values[1] == pytest.approx(math.exp(-1.0))
    assert values[-1] == 0.0 and values[-2] == 0.0
    # exp(-1/(1-0.99)) = exp(-100): tiny but nonzero, and every derivative
    # vanishes faster still, which is why the extension by zero is C^infinity.
    assert 0.0 < values[4] < 1.0e-40
    step = 1.0e-7
    inner = np.array([0.2, 0.5, 0.8])
    numeric = (smooth_bump(inner + step) - smooth_bump(inner - step)) / (2.0 * step)
    assert np.allclose(smooth_bump_derivative(inner), numeric, rtol=1.0e-5)


def test_family_rejects_bad_parameters():
    with pytest.raises(ValueError):
        SwirlFamily(amplitude=1.0, radial_support=0.0, axial_support=1.0)
    with pytest.raises(ValueError):
        SwirlFamily(amplitude=1.0, radial_support=1.0, axial_support=1.0,
                    axial_concentration=-1.0)
    with pytest.raises(ValueError):
        SwirlFamily(amplitude=math.nan, radial_support=1.0, axial_support=1.0)


def test_initial_data_is_compactly_supported_and_pure_swirl():
    grid = AxisymmetricGrid.uniform(
        nr=65, nz=129, r_max=2.5, z_min=-2.5, z_max=2.5, periodic_z=False
    )
    u1, omega1 = FAMILY.initial_state(grid)
    assert np.all(omega1 == 0.0)
    r_mesh, z_mesh = grid.mesh()
    outside = (r_mesh >= FAMILY.radial_support) | (
        np.abs(z_mesh) >= FAMILY.axial_support
    )
    assert np.all(u1[outside] == 0.0)
    # Odd in z, which is why the five-dimensional monopole of omega1 vanishes.
    mirrored = FAMILY.u1(r_mesh, -z_mesh)
    assert np.allclose(u1[:, :], -mirrored, atol=1.0e-15)


def test_initial_data_is_cartesian_smooth_across_the_axis():
    """``u = u_1(r^2,z)(-y,x,0)`` is smooth because ``u_1`` is a function of ``r^2``.

    A field built from ``chi(r/R_0)`` instead would have a corner at ``r=0``;
    the second difference across the axis detects that immediately.
    """
    z = 0.4
    h = 1.0e-3
    offsets = np.array([-2.0, -1.0, 0.0, 1.0, 2.0]) * h
    points = np.stack(
        (offsets, np.zeros_like(offsets), np.full_like(offsets, z)), axis=1
    )
    velocity = FAMILY.cartesian_velocity(points)
    # u_y(x,0,z) = x * u_1(|x|, z) must be smooth and odd in x.
    uy = velocity[:, 1]
    assert uy[2] == pytest.approx(0.0, abs=1.0e-15)
    assert uy[1] == pytest.approx(-uy[3], rel=1.0e-9)
    fourth = uy[0] - 4.0 * uy[1] + 6.0 * uy[2] - 4.0 * uy[3] + uy[4]
    assert abs(fourth) < 1.0e-9 * max(abs(uy).max(), 1.0e-30)


def test_initial_data_is_divergence_free_on_a_cartesian_grid():
    """``div u = 0`` exactly in the continuum; the stencil converges to it.

    The analytic cancellation is ``-y u_{1,x} + x u_{1,y} = 0``, which a centred
    difference reproduces only to ``O(h^2)``.  The bump has very large high
    derivatives near the edge of its support, so the *constant* in that ``O(h^2)``
    is large and only the observed order is meaningful.
    """
    errors = []
    for points in (41, 81, 161):
        axis = np.linspace(-1.6, 1.6, points)
        grid = UniformCartesianGrid(x=axis, y=axis, z=axis)
        x, y, z = grid.mesh()
        swirl = FAMILY.u1(np.hypot(x, y), z)
        velocity = np.stack((-y * swirl, x * swirl, np.zeros_like(swirl)), axis=0)
        divergence = cartesian_divergence(grid, velocity)
        core = (slice(2, -2), slice(2, -2), slice(2, -2))
        errors.append(
            float(np.max(np.abs(divergence[core])))
            / float(np.max(np.abs(velocity)))
        )
    orders = [
        math.log(errors[i] / errors[i + 1]) / math.log(2.0)
        for i in range(len(errors) - 1)
    ]
    assert errors[0] > errors[1] > errors[2]
    # The observed order climbs toward two (1.58 then 1.81 at the tested grids);
    # the sharp edge of the bump keeps the coarse pair out of the asymptotic
    # range, so the increase is the meaningful statement.
    assert min(orders) > 1.4
    assert orders[-1] > orders[0]


def test_initial_data_has_finite_energy_and_l3_and_axis_regularity():
    grid = AxisymmetricGrid.uniform(
        nr=129, nz=257, r_max=2.5, z_min=-2.5, z_max=2.5, periodic_z=False
    )
    assert 0.0 < FAMILY.kinetic_energy(grid) < math.inf
    assert 0.0 < FAMILY.l3_norm(grid) < math.inf
    # u^theta = r u_1 vanishes on the axis, and u_1 is even in r.
    assert FAMILY.swirl_velocity(np.array([0.0]), np.array([0.3]))[0] == 0.0
    assert FAMILY.u1(np.array([0.05]), np.array([0.3]))[0] == pytest.approx(
        FAMILY.u1(np.array([-0.05]), np.array([0.3]))[0], rel=1.0e-14
    )


def test_family_labels_and_aspect_ratio():
    assert FAMILY.aspect_ratio == pytest.approx(1.2 / 1.5)
    assert FAMILY.support_radius == pytest.approx(math.hypot(1.2, 1.5))
    assert "A10" in FAMILY.label()


# --------------------------------------------------------------------------- #
# Tail propagation                                                             #
# --------------------------------------------------------------------------- #


def test_velocity_error_bounds_are_the_stated_formulas():
    potential = PotentialErrorBounds(
        value=1.0e-3, gradient=2.0e-3, hessian=5.0e-3, radial_extent=4.0
    )
    radial, axial = velocity_error_bounds(potential)
    assert radial == pytest.approx(4.0 * 2.0e-3)
    assert axial == pytest.approx(2.0 * 1.0e-3 + 4.0 * 2.0e-3)


def test_potential_and_field_bounds_reject_bad_input():
    with pytest.raises(ValueError):
        PotentialErrorBounds(value=-1.0, gradient=0.0, hessian=0.0, radial_extent=1.0)
    with pytest.raises(ValueError):
        PotentialErrorBounds(value=0.0, gradient=0.0, hessian=0.0, radial_extent=0.0)
    with pytest.raises(ValueError):
        FieldBounds(
            u1_max=-1.0, omega1_max=0.0, u1_gradient_max=0.0,
            omega1_gradient_max=0.0, velocity_radial_max=0.0,
            velocity_axial_max=0.0, psi_axial_derivative_max=0.0,
        )


def test_advection_defect_grows_with_the_potential_error():
    fields = FieldBounds(
        u1_max=1.0, omega1_max=2.0, u1_gradient_max=3.0, omega1_gradient_max=4.0,
        velocity_radial_max=0.5, velocity_axial_max=0.6, psi_axial_derivative_max=0.7,
    )
    small = advection_defect_bound(
        PotentialErrorBounds(value=1e-6, gradient=1e-6, hessian=1e-6, radial_extent=2.0),
        fields,
    )
    large = advection_defect_bound(
        PotentialErrorBounds(value=1e-3, gradient=1e-3, hessian=1e-3, radial_extent=2.0),
        fields,
    )
    assert large.total > small.total
    assert large.advection_u1 > 0.0 and large.advection_omega1 > 0.0
    assert large.stretching_source == 0.0  # the tail does not enter it directly
    assert large.state_lipschitz == small.state_lipschitz  # state-only constant
    assert set(large.as_dict()) >= {"total", "u1_equation", "omega1_equation"}
    with pytest.raises(ValueError):
        advection_defect_bound(
            PotentialErrorBounds(value=0.0, gradient=0.0, hessian=0.0, radial_extent=1.0),
            fields, solve_operator_norm=-1.0,
        )


def test_gronwall_bound_matches_the_closed_form_and_is_monotone():
    bound = gronwall_state_error(
        initial_error=0.1, constant_defect=0.2, lipschitz=0.5, horizon=2.0
    )
    assert bound.bound == pytest.approx((0.1 + 0.2 * 2.0) * math.exp(0.5 * 2.0))
    zero = gronwall_state_error(
        initial_error=0.0, constant_defect=0.0, lipschitz=1.0, horizon=1.0
    )
    assert zero.bound == 0.0
    longer = gronwall_state_error(
        initial_error=0.1, constant_defect=0.2, lipschitz=0.5, horizon=3.0
    )
    assert longer.bound > bound.bound
    with pytest.raises(ValueError):
        gronwall_state_error(
            initial_error=-1.0, constant_defect=0.0, lipschitz=0.0, horizon=1.0
        )


# --------------------------------------------------------------------------- #
# Multipole hierarchy                                                          #
# --------------------------------------------------------------------------- #


def test_quadrupole_trace_beats_dipole_which_beats_monopole():
    field = ManufacturedFreeSpaceField(
        bumps=(
            BumpSpec(center=0.5, radius=1.2, amplitude=1.0),
            BumpSpec(center=-0.8, radius=0.9, amplitude=0.6),
        )
    )
    grid = AxisymmetricGrid.uniform(
        nr=257, nz=513, r_max=6.0, z_min=-6.0, z_max=6.0, periodic_z=False
    )
    r_mesh, z_mesh = grid.mesh()
    moments = source_moments(grid, field.omega1(r_mesh, z_mesh))
    angle = np.linspace(0.0, math.pi, 61)
    radius = 5.0
    rr, zz = radius * np.sin(angle), radius * np.cos(angle)
    errors = {}
    for truncation in ("monopole", "dipole", "quadrupole"):
        approx = multipole_boundary_trace(grid, moments, rr, zz, truncation=truncation)
        errors[truncation] = float(np.max(np.abs(field.psi1(rr, zz) - approx)))
        bound = multipole_tail_bounds(
            moments, evaluation_radius=radius, truncation=truncation
        )
        assert errors[truncation] <= bound.value
    assert errors["quadrupole"] < errors["dipole"] < errors["monopole"]


def test_solver_accepts_every_boundary_mode_and_rejects_others():
    grid = AxisymmetricGrid.uniform(
        nr=33, nz=65, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=False
    )
    for mode in ("zero", "monopole", "dipole", "quadrupole"):
        WholeSpaceEllipticSolver.build(grid, boundary_mode=mode)
    with pytest.raises(ValueError):
        WholeSpaceEllipticSolver.build(grid, boundary_mode="octupole")


def test_tail_bound_reports_unavailability_when_the_source_reaches_the_wall():
    grid = AxisymmetricGrid.uniform(
        nr=33, nz=65, r_max=1.4, z_min=-1.7, z_max=1.7, periodic_z=False
    )
    evolution = WholeSpaceEvolution(grid, viscosity=5.0e-3, boundary_mode="dipole")
    u1, omega1 = FAMILY.initial_state(grid)
    for _ in range(5):
        u1, omega1 = evolution.step(u1, omega1, 4.0e-4, "rk4")
    recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    assert not recovery.tail_bound_available
    assert recovery.interior_tail.as_dict()["value"] is None


# --------------------------------------------------------------------------- #
# Domain expansion                                                             #
# --------------------------------------------------------------------------- #


def test_expansion_policy_validates_and_reports_reasons():
    policy = ExpansionPolicy()
    fire, reasons = expansion_is_required(
        policy, outer_band_fraction=0.0, tail_bound_ratio=0.0,
        multipole_remainder_ratio=0.0, boundary_difference_ratio=0.0,
    )
    assert not fire and reasons == ()
    fire, reasons = expansion_is_required(
        policy, outer_band_fraction=1.0, tail_bound_ratio=1.0,
        multipole_remainder_ratio=1.0, boundary_difference_ratio=1.0,
    )
    assert fire and len(reasons) == 4
    with pytest.raises(ValueError):
        ExpansionPolicy(outer_band_fraction=0.0)
    with pytest.raises(ValueError):
        ExpansionPolicy(radial_growth_cells=0)


def test_embedding_is_zero_padding_with_aligned_nodes():
    grid = AxisymmetricGrid.uniform(
        nr=33, nz=65, r_max=2.0, z_min=-2.0, z_max=2.0, periodic_z=False
    )
    u1, omega1 = FAMILY.initial_state(grid)
    policy = ExpansionPolicy(radial_growth_cells=8, axial_growth_cells=8)
    new_grid, (new_u1, new_omega1), defect = embed_in_larger_box(
        grid, (u1, omega1), policy
    )
    assert defect < 1.0e-12
    assert new_grid.nr == grid.nr + 8 and new_grid.nz == grid.nz + 16
    assert new_grid.dr == pytest.approx(grid.dr) and new_grid.dz == pytest.approx(grid.dz)
    assert np.array_equal(new_u1[: grid.nr, 8 : 8 + grid.nz], u1)
    assert float(np.max(np.abs(new_u1))) == pytest.approx(float(np.max(np.abs(u1))))
    padded = new_u1.copy()
    padded[: grid.nr, 8 : 8 + grid.nz] = 0.0
    assert np.all(padded == 0.0)
    with pytest.raises(ValueError):
        periodic = AxisymmetricGrid.uniform(
            nr=17, nz=32, r_max=2.0, z_min=-2.0, z_max=2.0, periodic_z=True
        )
        embed_in_larger_box(periodic, (np.zeros(periodic.shape),), policy)


# --------------------------------------------------------------------------- #
# The interval certificate                                                     #
# --------------------------------------------------------------------------- #


def test_interval_arithmetic_is_exact_and_outward():
    a = Interval(Fraction(1, 3), Fraction(1, 2))
    b = Interval(Fraction(-1, 5), Fraction(1, 7))
    assert (a + b).lower == Fraction(1, 3) - Fraction(1, 5)
    assert (a * b).lower == min(
        Fraction(1, 3) * Fraction(-1, 5), Fraction(1, 2) * Fraction(-1, 5)
    )
    assert (a - a).contains(0)
    assert a.magnitude == Fraction(1, 2)
    assert b.contains_zero
    rounded = a.round_outward(8)
    assert rounded.lower <= a.lower and rounded.upper >= a.upper
    assert Interval.exact(0.5).lower == Fraction(1, 2)
    with pytest.raises(ValueError):
        Interval(Fraction(1), Fraction(0))


def test_snapshot_certificate_encloses_the_float_result_and_verifies():
    grid = AxisymmetricGrid.uniform(
        nr=25, nz=49, r_max=2.0, z_min=-2.4, z_max=2.4, periodic_z=False
    )
    evolution = WholeSpaceEvolution(grid, viscosity=5.0e-3, boundary_mode="dipole")
    u1, omega1 = FAMILY.initial_state(grid)
    for _ in range(5):
        u1, omega1 = evolution.step(u1, omega1, 4.0e-4, "rk4")
    recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    certificate = build_snapshot_certificate(
        grid, u1, omega1, recovery.psi1, viscosity=5.0e-3,
        green_tail_bound=1.0e-6, derivative_tail_bound=1.0e-5,
        interior_radius=1.4, precision_bits=64,
    )
    payload = json.loads(json.dumps(certificate.as_dict()))
    checks = verify_snapshot_certificate(payload)
    assert checks["all_passed"]
    assert checks["energy_encloses_float"]
    assert checks["residual_dominates_float"]
    # The residual bound is a genuine rational upper bound, not a float compare.
    assert Fraction(payload["poisson_residual_upper"]) >= 0
    assert payload["interior_nodes"] > 0


def test_snapshot_certificate_checker_rejects_a_tampered_payload():
    grid = AxisymmetricGrid.uniform(
        nr=17, nz=33, r_max=2.0, z_min=-2.4, z_max=2.4, periodic_z=False
    )
    evolution = WholeSpaceEvolution(grid, viscosity=5.0e-3, boundary_mode="dipole")
    u1, omega1 = FAMILY.initial_state(grid)
    u1, omega1 = evolution.step(u1, omega1, 4.0e-4, "rk4")
    recovery = recover_free_space_velocity(
        evolution.solver, omega1, support_threshold=1.0e-10
    )
    payload = json.loads(
        json.dumps(
            build_snapshot_certificate(
                grid, u1, omega1, recovery.psi1, viscosity=5.0e-3,
                green_tail_bound=1.0e-6, derivative_tail_bound=1.0e-5,
                interior_radius=1.2, precision_bits=64,
            ).as_dict()
        )
    )
    assert verify_snapshot_certificate(payload)["all_passed"]
    tampered = dict(payload)
    tampered["poisson_residual_upper"] = "0"
    tampered["float_witness"] = {**payload["float_witness"], "poisson_residual": 1.0}
    assert not verify_snapshot_certificate(tampered)["residual_dominates_float"]
    shrunk = dict(payload)
    shrunk["energy_interval"] = ["0", "0"]
    assert not verify_snapshot_certificate(shrunk)["energy_encloses_float"]
    with pytest.raises(ValueError):
        verify_snapshot_certificate({"kind": "something else"})


def test_snapshot_certificate_refuses_an_oversized_grid():
    grid = AxisymmetricGrid.uniform(
        nr=257, nz=513, r_max=2.0, z_min=-2.0, z_max=2.0, periodic_z=False
    )
    zeros = np.zeros(grid.shape)
    with pytest.raises(ValueError):
        build_snapshot_certificate(
            grid, zeros, zeros, zeros, viscosity=1.0e-3,
            green_tail_bound=0.0, derivative_tail_bound=0.0, interior_radius=1.0,
        )
