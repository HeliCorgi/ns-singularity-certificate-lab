"""Tests for the Gaussian--Hermite mixed family and its rigorous tail bound."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
import json
import math

import numpy as np
import pytest

from ns_certificate_lab.gaussian_hermite import (
    GaussianGenerator,
    GaussianMixedFamily,
    GaussianPolynomial,
    gaussian_search_basis,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.l3_generation import (
    AxisymmetricPressureSolver,
    l3_generation_rate,
)
from ns_certificate_lab.l3_optimizer import evaluate_shape, require_clay_admissible


def _grid(nr: int = 41, nz: int = 81) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=3.0, z_min=-3.6, z_max=3.6, periodic_z=False
    )


#: A generic amplitude vector: every slot nonzero, no accidental symmetry, so a
#: cancellation between generators cannot hide a wrong formula in one of them.
GENERIC_AMPLITUDES = (1.0, -0.7, 0.5, 0.6, 0.8, 0.4, -0.5, 0.9)


# --------------------------------------------------------------------------- #
# the polynomial recursions and the generator partials                         #
# --------------------------------------------------------------------------- #


def test_polynomial_derivative_recursions_are_exact() -> None:
    """The recursions are the whole mechanism behind the analytic gradient, so
    they are pinned on a case small enough to differentiate by hand:
    ``P = s z`` with ``E = exp(-2s - 3z^2)`` gives ``d_s(PE) = (z - 2sz)E`` and
    ``d_z(PE) = (s - 6sz^2)E``.  A wrong coefficient here would poison every
    gradient and every tail bound downstream."""
    poly = GaussianPolynomial({(1, 1): 1})
    assert poly.d_s(Fraction(2)) == GaussianPolynomial({(0, 1): 1, (1, 1): -2})
    assert poly.d_z(Fraction(3)) == GaussianPolynomial({(1, 0): 1, (1, 2): -6})
    # Rational arithmetic must be exact, not floating: 1/3 stays 1/3.
    third = GaussianPolynomial({(0, 0): Fraction(1, 3)})
    assert third.d_s(Fraction(1, 7)) == GaussianPolynomial(
        {(0, 0): Fraction(-1, 21)}
    )


def test_polynomial_rejects_bad_exponents() -> None:
    with pytest.raises(ValueError):
        GaussianPolynomial({(-1, 0): 1})


def test_generator_partials_match_central_differences() -> None:
    """Each analytic partial is checked against a difference quotient of the
    generator, which catches a wrong recursion coefficient that the divergence
    identity would not see (the divergence cancels for *any* stream generator,
    right or wrong).  First derivatives are differenced from ``evaluate``; the
    second derivatives are differenced from the analytic first partials, so
    each recursion step is verified against the one below it."""
    generator = GaussianGenerator(
        polynomial=GaussianPolynomial(
            {(0, 1): Fraction(3, 2), (1, 1): -1, (2, 3): Fraction(1, 4)}
        ),
        alpha=Fraction(3, 2), beta=Fraction(4, 5), amplitude=1.3,
    )
    r = np.array([0.35, 0.7, 1.1])
    z = np.array([-0.4, 0.2, 0.8])
    s = r * r
    partials = generator.partials(r, z)
    step = 1.0e-5

    def at(s_values: np.ndarray, z_values: np.ndarray) -> dict[str, np.ndarray]:
        return generator.partials(np.sqrt(s_values), z_values)

    fd_ds = (
        generator.evaluate(np.sqrt(s + step), z)
        - generator.evaluate(np.sqrt(s - step), z)
    ) / (2.0 * step)
    fd_dz = (
        generator.evaluate(r, z + step) - generator.evaluate(r, z - step)
    ) / (2.0 * step)
    fd_dss = (at(s + step, z)["ds"] - at(s - step, z)["ds"]) / (2.0 * step)
    fd_dzz = (at(s, z + step)["dz"] - at(s, z - step)["dz"]) / (2.0 * step)
    fd_dsz = (at(s, z + step)["ds"] - at(s, z - step)["ds"]) / (2.0 * step)
    assert np.allclose(partials["ds"], fd_ds, rtol=1e-7, atol=1e-10)
    assert np.allclose(partials["dz"], fd_dz, rtol=1e-7, atol=1e-10)
    assert np.allclose(partials["dss"], fd_dss, rtol=1e-6, atol=1e-9)
    assert np.allclose(partials["dzz"], fd_dzz, rtol=1e-6, atol=1e-9)
    assert np.allclose(partials["dsz"], fd_dsz, rtol=1e-6, atol=1e-9)


def test_generator_rejects_nonpositive_widths() -> None:
    with pytest.raises(ValueError):
        GaussianGenerator(
            polynomial=GaussianPolynomial({(0, 1): 1}), alpha=Fraction(0),
            beta=Fraction(1),
        )


# --------------------------------------------------------------------------- #
# the constructed field                                                        #
# --------------------------------------------------------------------------- #


def test_gaussian_field_is_divergence_free_to_machine_precision() -> None:
    """The construction makes the divergence vanish as an algebraic identity,
    so what is measured is floating-point cancellation, not truncation.  A
    value of order the gradient size would mean the analytic partials are
    wrong -- the same gate the chi families pass, and it must hold here with a
    completely different differentiation mechanism behind it."""
    grid = _grid()
    family = gaussian_search_basis().with_amplitudes(GENERIC_AMPLITUDES)
    analytic = family.field(grid)
    divergence = analytic.divergence()
    scale = float(np.max(np.abs(analytic.gradient()["rr"])))
    assert float(np.max(np.abs(divergence))) <= 1.0e-12 * scale


def test_analytic_and_difference_gradients_converge_at_second_order() -> None:
    """The analytic gradient is the certified one; the difference operators are
    the independent check, and they must converge to it at the operators' own
    order or one of the two is wrong.

    The refinement starts at ``nr = 33`` deliberately, in pointed contrast to
    the chi families, which need ``nr = 97`` before the comparison even enters
    its asymptotic range and never do better than order ``1.55``: Gaussians
    have mild derivatives, so second order must appear immediately.  Measured
    orders per step are 1.78, 1.95, 1.99.  Failure to clear 1.8 here would mean
    either a wrong analytic partial or that the family has lost the very
    mildness that justifies its existence."""
    family = gaussian_search_basis().with_amplitudes(GENERIC_AMPLITUDES)
    errors = []
    for n in (33, 65, 129, 257):
        grid = AxisymmetricGrid.uniform(
            nr=n, nz=2 * n - 1, r_max=3.0, z_min=-3.6, z_max=3.6, periodic_z=False
        )
        field = family.field(grid)
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
    order = math.log(errors[0] / errors[-1]) / math.log(256.0 / 32.0)
    assert order > 1.8


# --------------------------------------------------------------------------- #
# parity                                                                       #
# --------------------------------------------------------------------------- #


def test_all_stream_generators_of_the_basis_are_odd_in_z() -> None:
    """The parity rule kills the pressure term for an even stream generator, so
    an even stream slot would be a design variable structurally incapable of
    the one thing the basis exists for.  The parity is read off the exponent
    dictionary, so this is a structural check, not a sampling."""
    basis = gaussian_search_basis()
    for generator in basis.stream:
        assert generator.axial_parity == "odd"


def test_even_stream_gaussian_family_cancels_the_pressure_term_numerically() -> None:
    """The selection rule was derived for the chi families by a parity argument
    that never used compact support; this checks the derivation transfers to
    Gaussian generators by measuring the actual integral.  A surviving pressure
    term would mean the rule's premise (the recovery map's parity bookkeeping)
    fails off the compactly supported class."""
    grid = _grid(49, 97)
    solver = AxisymmetricPressureSolver.build(grid)
    even = GaussianMixedFamily(
        name="even-stream-gaussian",
        swirl=(
            GaussianGenerator(
                polynomial=GaussianPolynomial({(0, 1): 1}),
                alpha=Fraction(3), beta=Fraction(3), amplitude=10.0,
            ),
        ),
        stream=(
            GaussianGenerator(
                polynomial=GaussianPolynomial({(0, 0): 1, (1, 0): Fraction(-1, 2)}),
                alpha=Fraction(5, 2), beta=Fraction(5, 2), amplitude=4.0,
            ),
        ),
        reference_length=1.0,
    )
    assert even.stream[0].axial_parity == "even"
    rate = l3_generation_rate(even.field(grid), viscosity=5.0e-3, solver=solver)
    assert abs(rate.pressure) < 1.0e-12 * abs(rate.viscous)


# --------------------------------------------------------------------------- #
# the tail bound                                                               #
# --------------------------------------------------------------------------- #


def test_tail_bound_is_a_true_bound_on_spheres() -> None:
    """The bound claims to dominate ``sup |u|`` outside every sphere; it is
    sampled against the actual field on the spheres themselves, where the
    exterior supremum is attained or approached.  A single sampled speed above
    the bound would mean the bound is not a certificate and every use of it --
    the support radius, the box-truncation argument -- collapses."""
    family = gaussian_search_basis().with_amplitudes(GENERIC_AMPLITUDES)
    angles = np.linspace(0.0, np.pi, 721)
    previous = math.inf
    for radius in (1.0, 1.5, 2.0, 2.5, 3.0):
        u_r, u_theta, u_z = family.velocity(
            radius * np.sin(angles), radius * np.cos(angles)
        )
        sampled = float(np.max(np.sqrt(u_r**2 + u_theta**2 + u_z**2)))
        bound = family.tail_bound(radius)
        assert sampled <= bound
        # Monotone nonincreasing, as the exterior supremum must be.
        assert bound <= previous
        previous = bound
    with pytest.raises(ValueError):
        family.tail_bound(-1.0)


def test_tail_bound_meets_the_preregistered_truncation_level() -> None:
    """The docstring of :func:`gaussian_search_basis` commits to a truncation
    level below ``1e-7`` relative at the standard box, whose exterior lies
    beyond radius 3.  This pins that commitment to the implementation: if a
    width or a bound coefficient drifts, the preregistration and the code
    disagree and this test says so."""
    family = gaussian_search_basis().with_amplitudes(np.ones(8))
    grid = _grid()
    scale = float(np.max(family.field(grid).speed))
    assert family.tail_bound(3.0) <= 1.0e-7 * scale


def test_support_radius_certifies_its_threshold() -> None:
    """``support_radius`` must return a radius at which the bound itself has
    been evaluated and found below the threshold -- conservative by
    construction, so this failing would mean the bisection returned the wrong
    side."""
    family = gaussian_search_basis().with_amplitudes(GENERIC_AMPLITUDES)
    threshold = 1.0e-9
    radius = family.support_radius(threshold)
    assert family.tail_bound(radius) <= threshold
    # The bound is far above the threshold one length unit closer in, i.e. the
    # returned radius is meaningfully tight rather than a lazy overestimate.
    assert family.tail_bound(max(radius - 1.0, 0.0)) > threshold
    # An identically zero family needs no exclusion sphere at all.
    assert gaussian_search_basis().support_radius(threshold) == 0.0
    with pytest.raises(ValueError):
        family.support_radius(0.0)


# --------------------------------------------------------------------------- #
# the optimiser contract                                                       #
# --------------------------------------------------------------------------- #


def test_evaluate_shape_scores_the_gaussian_seed() -> None:
    """The whole point of duck-typing the family: the optimiser's entry point
    must consume it unchanged.  The seed mirrors the deterministic start of
    ``optimise`` -- first swirl and first stream slot -- and a non-finite
    critical Reynolds number would mean the Gaussian pressure term vanished,
    i.e. the parity bookkeeping or the field construction is broken."""
    grid = _grid(33, 65)
    score = evaluate_shape(
        np.array([1.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0]),
        basis=gaussian_search_basis(), grid=grid,
    )
    assert score.pressure > 0.0
    assert score.sign in (1.0, -1.0)
    assert math.isfinite(score.critical_reynolds)
    assert score.critical_reynolds > 0.0


def test_with_amplitudes_round_trips_and_rejects_wrong_length() -> None:
    basis = gaussian_search_basis()
    values = np.arange(1.0, len(basis.amplitudes) + 1.0)
    assert np.allclose(basis.with_amplitudes(values).amplitudes, values)
    # The original is untouched: the family is a value, not a mutable buffer.
    assert np.allclose(basis.amplitudes, 0.0)
    with pytest.raises(ValueError):
        basis.with_amplitudes([1.0, 2.0])


def test_family_is_clay_admissible_and_the_guard_accepts_it() -> None:
    """Gaussian--Hermite data are Schwartz class, so unlike a spline surrogate
    they are legitimate smooth decaying data and the guard must let them
    through; and the guard must still bite when the flag is withdrawn, or the
    admissibility firewall has silently stopped working for this type."""
    family = gaussian_search_basis().with_amplitudes(GENERIC_AMPLITUDES)
    assert family.clay_admissible is True
    require_clay_admissible(family)
    with pytest.raises(ValueError):
        require_clay_admissible(replace(family, clay_admissible=False))


def test_family_serialises_with_exact_rational_data() -> None:
    """``as_dict`` must survive ``json.dumps`` (results files) and must carry
    the widths as exact fraction strings, because the tail bound's rigour rests
    on those rationals and a float round trip would quietly discard it."""
    family = gaussian_search_basis().with_amplitudes(GENERIC_AMPLITUDES)
    payload = family.as_dict()
    json.dumps(payload, allow_nan=False)
    assert payload["clay_admissible"] is True
    assert payload["stream"][3]["alpha"] == "5/2"
    assert all(entry["axial_parity"] == "odd" for entry in payload["stream"])
