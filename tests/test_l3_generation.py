"""Tests for the critical-norm generation rate, the mixed families and the search."""

from __future__ import annotations

from fractions import Fraction
import json
import math

import numpy as np
import pytest

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.initial_data import family_S
from ns_certificate_lab.l3_certificate import (
    bilinear_cell_range,
    exp_interval,
    sqrt_interval,
)
from ns_certificate_lab.l3_generation import (
    AxisymmetricPressureSolver,
    MixedField,
    compact_poisson_residual,
    l3_generation_rate,
    parity_selection_rule,
    pressure_exterior_tail_bound,
    pure_swirl_no_go,
)
from ns_certificate_lab.l3_optimizer import (
    KILL_CONDITION_CRITICAL_REYNOLDS,
    SEARCH_BASIS,
    critical_reynolds,
    evaluate_shape,
    normalise_amplitudes,
)
from ns_certificate_lab.mixed_initial_data import (
    GeneratorComponent,
    MIXED_FAMILIES,
    MixedFamily,
    family_M1,
    family_M2,
    family_M3,
)
from ns_certificate_lab.snapshot_certificate import Interval


def _grid(nr: int = 41, nz: int = 81) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=3.0, z_min=-3.6, z_max=3.6, periodic_z=False
    )


# --------------------------------------------------------------------------- #
# the three-dimensional pressure solve                                          #
# --------------------------------------------------------------------------- #


def test_pressure_solver_uses_the_three_dimensional_coefficients() -> None:
    """The 3D axisymmetric Laplacian has radial drift 1 and axis row 4, against
    3 and 8 for the five-dimensional operator.  Silently using the 5D values
    solves a different equation and would give a wrong pressure everywhere."""
    assert AxisymmetricPressureSolver.RADIAL_DRIFT == 1.0
    assert AxisymmetricPressureSolver.AXIS_COEFFICIENT == 4.0


def test_pressure_solver_is_second_order_on_a_manufactured_solution() -> None:
    """``p = exp(-(r^2+z^2))`` gives ``Laplacian_3 p = (4(r^2+z^2) - 6) p``.
    A first-order rate would mean the axis row or the drift term is wrong."""
    errors = []
    for n in (33, 65, 129):
        grid = AxisymmetricGrid.uniform(
            nr=n, nz=2 * n - 1, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=False
        )
        r, z = grid.mesh()
        exact = np.exp(-(r**2 + z**2))
        source = -(4.0 * (r**2 + z**2) - 6.0) * exact
        solver = AxisymmetricPressureSolver.build(grid)
        numeric = solver.solve(
            source, trace=(exact[-1, :], exact[:, 0], exact[:, -1])
        )
        errors.append(float(np.max(np.abs(numeric - exact))))
    for coarse, fine in zip(errors, errors[1:]):
        assert math.log(coarse / fine) / math.log(2.0) > 1.8


def test_pressure_solver_fault_injection_is_detectable() -> None:
    """Building the solver with the five-dimensional coefficients must change
    the answer by far more than the discretisation error."""
    grid = AxisymmetricGrid.uniform(
        nr=49, nz=97, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=False
    )
    r, z = grid.mesh()
    exact = np.exp(-(r**2 + z**2))
    source = -(4.0 * (r**2 + z**2) - 6.0) * exact
    trace = (exact[-1, :], exact[:, 0], exact[:, -1])
    good = AxisymmetricPressureSolver.build(grid).solve(source, trace=trace)
    bad = AxisymmetricPressureSolver.build(
        grid, axis_coefficient=8.0, radial_drift_coefficient=3.0
    ).solve(source, trace=trace)
    assert np.max(np.abs(bad - exact)) > 20.0 * np.max(np.abs(good - exact))


def test_pressure_solver_rejects_a_periodic_grid() -> None:
    with pytest.raises(ValueError):
        AxisymmetricPressureSolver.build(
            AxisymmetricGrid.uniform(nr=17, nz=17, r_max=1.0)
        )


def test_pressure_tail_bound_scales_like_the_inverse_cube() -> None:
    """The monopole and dipole of the pressure source vanish for a compactly
    supported divergence-free field, so the pressure decays like the quadrupole."""
    first = pressure_exterior_tail_bound(1.0, 2.0)
    second = pressure_exterior_tail_bound(1.0, 4.0)
    assert first / second == pytest.approx(8.0)
    with pytest.raises(ValueError):
        pressure_exterior_tail_bound(1.0, 0.0)


# --------------------------------------------------------------------------- #
# the mixed families                                                           #
# --------------------------------------------------------------------------- #


def test_mixed_families_are_divergence_free_to_machine_precision() -> None:
    """The construction makes the divergence vanish as an algebraic identity, so
    what is measured here is the difference operator's truncation error and
    nothing else.  A value of order the field size would mean the analytic
    partials are wrong."""
    grid = _grid()
    for builder in MIXED_FAMILIES.values():
        field = builder()
        analytic = field.field(grid)
        divergence = analytic.divergence()
        scale = float(np.max(np.abs(analytic.gradient()["rr"])))
        assert float(np.max(np.abs(divergence))) <= 1.0e-12 * scale


def test_analytic_and_difference_gradients_agree_under_refinement() -> None:
    """The analytic gradient is the certified one; the difference gradient is the
    independent check.  They must converge to each other, or one of them is wrong.

    The refinement starts at ``nr = 97`` deliberately.  Measured on this family
    the difference is 1.33, 1.54, 1.76 at ``nr = 33, 49, 65`` and only then falls
    to 0.87, 0.74, 0.36, 0.22 at 97, 129, 193, 257: below about ``nr = 97`` the
    edge of the ``chi`` bump is not resolved at all and the comparison is not in
    its asymptotic range, so it would be measuring nothing.  The observed order
    over the resolved range is about 1.55, below two, which is the known
    behaviour of these data -- ``chi`` has enormous high derivatives near the
    edge of its support.  This is precisely why the certified gradient is the
    analytic one and the difference operator is only the cross-check.
    """
    errors = []
    for n in (97, 129, 193, 257):
        grid = AxisymmetricGrid.uniform(
            nr=n, nz=2 * n - 1, r_max=3.0, z_min=-3.6, z_max=3.6, periodic_z=False
        )
        field = family_M2().field(grid)
        exact = field.gradient()
        approximate = field.finite_difference_gradient()
        interior = (slice(2, -2), slice(2, -2))
        errors.append(
            max(
                float(np.max(np.abs((exact[k] - approximate[k])[interior])))
                for k in exact
            )
        )
    assert errors == sorted(errors, reverse=True)
    order = math.log(errors[0] / errors[-1]) / math.log(256.0 / 96.0)
    assert order > 1.3


def test_generator_partials_match_central_differences() -> None:
    """Each analytic partial is checked against a difference quotient of the
    generator itself, which catches a wrong chain-rule factor that the
    divergence identity would not see."""
    component = GeneratorComponent(
        amplitude=1.7, radial_support=1.1, axial_support=0.9,
        axial_center=0.2, axial_concentration=0.4, odd_axial=True,
    )
    r = np.array([0.3, 0.6, 0.9])
    z = np.array([-0.2, 0.1, 0.5])
    partials = component.partials(r, z)
    step = 1.0e-6
    ds = (
        component.evaluate(np.sqrt(r**2 + step), z)
        - component.evaluate(np.sqrt(np.maximum(r**2 - step, 0.0)), z)
    ) / (2.0 * step)
    dz = (component.evaluate(r, z + step) - component.evaluate(r, z - step)) / (
        2.0 * step
    )
    assert np.allclose(partials["ds"], ds, rtol=1e-4, atol=1e-8)
    assert np.allclose(partials["dz"], dz, rtol=1e-5, atol=1e-9)


def test_family_M1_recovers_the_pure_swirl_baseline(  # noqa: N802
) -> None:
    """M1 with zero meridional amplitude must reproduce family S exactly, or the
    bridge back to every Gate 7 result is broken."""
    grid = _grid()
    r, z = grid.mesh()
    assert np.array_equal(
        family_M1(swirl_amplitude=7.0, meridional_amplitude=0.0).u1(r, z),
        family_S(amplitude=7.0).u1(r, z),
    )


def test_all_stream_generators_are_odd_in_z() -> None:
    """The parity rule kills the pressure term for an even stream generator, so a
    family that shipped one would be structurally incapable of the thing it is
    for."""
    for builder in MIXED_FAMILIES.values():
        for component in builder().stream:
            assert component.odd_axial is True


def test_families_are_compactly_supported_and_axis_regular() -> None:
    grid = _grid()
    r, z = grid.mesh()
    for builder in MIXED_FAMILIES.values():
        family = builder()
        field = family.field(grid)
        outside = np.hypot(r, z) > family.support_radius
        assert np.all(field.speed[outside] == 0.0)
        # u^r and u^theta vanish on the axis; u^z need not.
        assert np.all(field.u_r[0] == 0.0)
        assert np.all(np.abs(field.u_theta[0]) < 1.0e-14)


def test_family_with_amplitudes_round_trips() -> None:
    family = family_M3()
    values = np.arange(1.0, len(family.amplitudes) + 1.0)
    assert np.allclose(family.with_amplitudes(values).amplitudes, values)
    with pytest.raises(ValueError):
        family.with_amplitudes([1.0])


# --------------------------------------------------------------------------- #
# the generation rate                                                          #
# --------------------------------------------------------------------------- #


def test_transport_contribution_is_negligible_against_the_others() -> None:
    """The identity says it is exactly zero.  A residual comparable to the
    pressure term would mean the quadrature cannot resolve the pressure term
    either, so this is the gate on trusting ``J`` at all."""
    grid = _grid(65, 129)
    solver = AxisymmetricPressureSolver.build(grid)
    field = family_M2().field(grid)
    rate = l3_generation_rate(field, viscosity=5.0e-3, solver=solver)
    assert abs(rate.transport_residual) < 0.05 * abs(rate.viscous)


def test_pure_swirl_has_exactly_zero_pressure_contribution() -> None:
    """This is the no-go, measured: for ``u^r = u^z = 0`` the pressure integrand
    is identically zero pointwise, not merely small."""
    grid = _grid()
    field = family_M1(meridional_amplitude=0.0).field(grid)
    verdict = pure_swirl_no_go(field)
    assert verdict["is_pure_swirl"] is True
    rate = l3_generation_rate(field, viscosity=5.0e-3)
    assert rate.pressure == 0.0
    assert rate.transport_residual == 0.0
    assert rate.total < 0.0


def test_pure_swirl_no_go_records_what_it_does_not_exclude() -> None:
    grid = _grid()
    verdict = pure_swirl_no_go(family_M1(meridional_amplitude=0.0).field(grid))
    assert any("later time" in item for item in verdict["does_not_exclude"])
    assert any("Clay" in item for item in verdict["does_not_exclude"])


def test_viscous_contribution_is_never_positive() -> None:
    grid = _grid()
    solver = AxisymmetricPressureSolver.build(grid)
    for builder in MIXED_FAMILIES.values():
        rate = l3_generation_rate(
            builder().field(grid), viscosity=1.0e-2, solver=solver
        )
        assert rate.viscous <= 0.0


def test_pressure_is_odd_and_viscous_is_even_under_a_sign_flip() -> None:
    """``p`` is quadratic in ``u`` and ``div(|u|u)`` is odd, so ``P`` flips sign
    while ``V`` does not.  This is what lets the search read the good sign off
    rather than search for it, so it is checked to full precision."""
    grid = _grid()
    solver = AxisymmetricPressureSolver.build(grid)
    plus = l3_generation_rate(
        family_M2(10.0, 4.0).field(grid), viscosity=5.0e-3, solver=solver
    )
    minus = l3_generation_rate(
        family_M2(-10.0, -4.0).field(grid), viscosity=5.0e-3, solver=solver
    )
    assert minus.pressure == pytest.approx(-plus.pressure, rel=1e-10)
    assert minus.viscous == pytest.approx(plus.viscous, rel=1e-12)


def test_generation_rate_rejects_a_vanishing_viscosity() -> None:
    grid = _grid(17, 33)
    with pytest.raises(ValueError):
        l3_generation_rate(family_M1().field(grid), viscosity=0.0)


def test_generation_rate_serialises_without_non_finite_values() -> None:
    grid = _grid(25, 49)
    rate = l3_generation_rate(family_M1().field(grid), viscosity=1.0e-2)
    json.dumps(rate.as_dict(), allow_nan=False)


def test_compact_residual_is_at_machine_precision() -> None:
    """The residual uses the solver's own stencil, so anything above rounding
    level means the linear solve is wrong -- as opposed to the discretisation
    being inaccurate, which is a different question."""
    grid = _grid(49, 97)
    solver = AxisymmetricPressureSolver.build(grid)
    field = family_M3().field(grid)
    source = field.pressure_source()
    pressure = solver.solve(source)
    residual = compact_poisson_residual(grid, pressure, source)
    scale = float(np.max(np.abs(source)))
    assert float(np.max(residual[1:-1, 1:-1])) < 1.0e-10 * scale


# --------------------------------------------------------------------------- #
# the parity selection rule                                                    #
# --------------------------------------------------------------------------- #


def test_parity_rule_depends_only_on_the_stream_generator() -> None:
    for swirl in ("odd", "even", "none"):
        assert parity_selection_rule(swirl, "even")[
            "pressure_term_can_be_nonzero"
        ] is False
        assert parity_selection_rule(swirl, "odd")[
            "pressure_term_can_be_nonzero"
        ] is True
    with pytest.raises(ValueError):
        parity_selection_rule("odd", "sideways")


def test_even_stream_generator_cancels_the_pressure_term_numerically() -> None:
    """The rule is derived by a parity argument; this checks the derivation
    against the actual integral, which is the only way to catch a sign slip in
    the argument itself."""
    grid = _grid(49, 97)
    solver = AxisymmetricPressureSolver.build(grid)
    even = MixedFamily(
        name="even-stream",
        swirl=(
            GeneratorComponent(
                amplitude=10.0, radial_support=1.2, axial_support=1.5,
                axial_concentration=0.5, odd_axial=True,
            ),
        ),
        stream=(
            GeneratorComponent(
                amplitude=4.0, radial_support=0.9, axial_support=0.8,
                odd_axial=False,
            ),
        ),
        reference_length=1.2,
    )
    rate = l3_generation_rate(even.field(grid), viscosity=5.0e-3, solver=solver)
    assert abs(rate.pressure) < 1.0e-12 * abs(rate.viscous)


# --------------------------------------------------------------------------- #
# the search                                                                   #
# --------------------------------------------------------------------------- #


def test_critical_reynolds_uses_one_power_of_the_length() -> None:
    """``P ~ A^4 L^2`` and ``|V| ~ nu A^3 L`` with ``A`` the physical velocity
    amplitude, so the ratio carries ``AL/nu``.  Two powers of ``L`` would be the
    swirl-variable convention and would change the objective by a factor ``L``."""
    grid = _grid()
    field = family_M2(-10.0, -4.0).field(grid)
    rate = l3_generation_rate(field, viscosity=5.0e-3)
    assert rate.pressure > 0.0
    doubled = critical_reynolds(rate, 2.0 * 1.2)
    single = critical_reynolds(rate, 1.2)
    assert doubled == pytest.approx(2.0 * single)


def test_critical_reynolds_is_infinite_for_a_useless_pressure_sign() -> None:
    grid = _grid()
    rate = l3_generation_rate(family_M2(10.0, 4.0).field(grid), viscosity=5.0e-3)
    if rate.pressure <= 0.0:
        assert not math.isfinite(critical_reynolds(rate, 1.2))


def test_positive_generation_is_reached_above_the_critical_reynolds() -> None:
    """The whole point of the objective: below ``Re_crit`` the rate is negative
    and above it positive, with the crossing where the formula says."""
    grid = _grid(65, 129)
    solver = AxisymmetricPressureSolver.build(grid)
    field = family_M2(-10.0, -4.0).field(grid)
    amplitude = float(np.max(field.speed))
    length = 1.2
    probe = l3_generation_rate(field, viscosity=1.0e-3, solver=solver)
    assert probe.pressure > 0.0
    reynolds = critical_reynolds(probe, length)
    assert math.isfinite(reynolds)
    below = l3_generation_rate(
        field, viscosity=amplitude * length / (0.5 * reynolds), solver=solver
    )
    above = l3_generation_rate(
        field, viscosity=amplitude * length / (2.0 * reynolds), solver=solver
    )
    assert below.total < 0.0
    assert above.total > 0.0


def test_evaluate_shape_picks_the_sign_that_makes_the_pressure_positive() -> None:
    grid = _grid(33, 65)
    score = evaluate_shape(
        np.array([1.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0]), grid=grid
    )
    assert score.pressure > 0.0
    assert score.sign in (1.0, -1.0)
    assert math.isfinite(score.critical_reynolds)


def test_normalise_amplitudes_removes_the_scale_direction() -> None:
    values = np.array([3.0, -4.0])
    assert np.allclose(normalise_amplitudes(values), [0.6, -0.8])
    assert np.allclose(
        normalise_amplitudes(values), normalise_amplitudes(10.0 * values)
    )
    with pytest.raises(ValueError):
        normalise_amplitudes(np.zeros(3))


def test_search_basis_breaks_the_axial_symmetry() -> None:
    """A basis confined to the symmetric subspace could not leave it, so the
    optimiser would be searching a set on which the objective is structurally
    suppressed."""
    offsets = [c.axial_center for c in (*SEARCH_BASIS.swirl, *SEARCH_BASIS.stream)]
    assert any(offset != 0.0 for offset in offsets)


def test_kill_condition_is_recorded_and_positive() -> None:
    assert KILL_CONDITION_CRITICAL_REYNOLDS > 0.0


# --------------------------------------------------------------------------- #
# interval primitives used by the certificate                                  #
# --------------------------------------------------------------------------- #


def test_exp_interval_encloses_the_true_value() -> None:
    for value in (-6.0, -1.5, -0.3, 0.0, 0.7, 2.5):
        box = exp_interval(Interval(Fraction(value), Fraction(value)))
        assert float(box.lower) <= math.exp(value) <= float(box.upper)
        assert float(box.upper) - float(box.lower) < 1.0e-9 * max(
            1.0, math.exp(value)
        )


def test_exp_interval_is_monotone_and_rejects_short_series() -> None:
    box = exp_interval(Interval(Fraction(-1), Fraction(1)))
    assert float(box.lower) <= math.exp(-1.0)
    assert float(box.upper) >= math.exp(1.0)
    with pytest.raises(ValueError):
        exp_interval(Interval(Fraction(0), Fraction(0)), terms=2)


def test_sqrt_interval_encloses_and_rejects_negatives() -> None:
    box = sqrt_interval(Interval(Fraction(2), Fraction(9)))
    assert float(box.lower) <= math.sqrt(2.0)
    assert float(box.upper) >= 3.0
    with pytest.raises(ValueError):
        sqrt_interval(Interval(Fraction(-1), Fraction(1)))


def test_bilinear_cell_range_is_the_exact_corner_hull() -> None:
    """The bilinear interpolant is a convex combination of its corners, so the
    hull is the range, not merely a bound.  Widening it would be sound but
    wasteful; narrowing it would be wrong."""
    corners = (Fraction(1), Fraction(-2), Fraction(3), Fraction(0))
    box = bilinear_cell_range(*corners)
    assert box.lower == Fraction(-2)
    assert box.upper == Fraction(3)
