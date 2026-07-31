"""Tests for the Gaussian--Hermite interval enclosures and their certificate."""

from __future__ import annotations

import copy
from fractions import Fraction
import math

import numpy as np
import pytest

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.interval_gaussian import (
    GaussianFamily,
    GaussianTerm,
    as_gaussian_term,
    build_gaussian_generation_certificate,
    cell_enclosure_gaussian,
    clear_gaussian_cache,
    exterior_viscous_tail,
    gaussian_box,
    gaussian_tail_moment,
    interval_poly_eval,
    seed_gaussian_family,
)
from ns_certificate_lab.l3_certificate import verify_generation_certificate
from ns_certificate_lab.snapshot_certificate import Interval


def _box(lower: Fraction, upper: Fraction) -> Interval:
    return Interval(Fraction(lower), Fraction(upper))


# --------------------------------------------------------------------------- #
# floating-point reference evaluation, independent of the module's closure     #
# --------------------------------------------------------------------------- #
#
# The enclosures are checked against finite differences of the CLOSED FORM of
# each term, not against the module's own exact derivative polynomials: a sign
# slip in the derivative closure would poison both sides of a comparison that
# reused it, whereas a difference quotient of P e^{-alpha s - beta z^2} knows
# nothing about how the module differentiates.


def _value(terms, s: float, z: float) -> float:
    total = 0.0
    for term in terms:
        poly = sum(float(c) * s**m * z**n for m, n, c in term.monomials)
        total += poly * math.exp(-(float(term.alpha) * s + float(term.beta) * z * z))
    return total


def _fd_partials(terms, s: float, z: float, step: float = 1.0e-5) -> dict[str, float]:
    centre = _value(terms, s, z)
    return {
        "value": centre,
        "ds": (_value(terms, s + step, z) - _value(terms, s - step, z)) / (2 * step),
        "dss": (_value(terms, s + step, z) - 2 * centre + _value(terms, s - step, z))
        / step**2,
        "dz": (_value(terms, s, z + step) - _value(terms, s, z - step)) / (2 * step),
        "dzz": (_value(terms, s, z + step) - 2 * centre + _value(terms, s, z - step))
        / step**2,
        "dsz": (
            _value(terms, s + step, z + step) - _value(terms, s + step, z - step)
            - _value(terms, s - step, z + step) + _value(terms, s - step, z - step)
        ) / (4 * step**2),
    }


# --------------------------------------------------------------------------- #
# the sharp primitives                                                         #
# --------------------------------------------------------------------------- #


def test_gaussian_box_is_sharp_on_a_small_box() -> None:
    """The Gaussian's range over a box is attained at corners, so the enclosure
    must be the true range up to series and rounding slop -- the whole reason
    this basis replaces the flat bump.  Width ratio measured at 1.0000002."""
    s_box = _box(Fraction(1, 5), Fraction(2, 5))
    z_box = _box(Fraction(1, 10), Fraction(3, 10))
    alpha, beta = Fraction(2, 3), Fraction(3, 2)
    enclosure = gaussian_box(alpha, beta, s_box, z_box, bits=64)
    true_upper = math.exp(-float(alpha) * 0.2 - float(beta) * 0.01)
    true_lower = math.exp(-float(alpha) * 0.4 - float(beta) * 0.09)
    ratio = (float(enclosure.upper) - float(enclosure.lower)) / (true_upper - true_lower)
    assert ratio < 1.02
    # math.exp carries its own ulp of error, so containment of the float
    # reference is asserted only up to that granularity.
    assert float(enclosure.lower) <= true_lower + 1.0e-12
    assert true_upper <= float(enclosure.upper) + 1.0e-12


def test_gaussian_box_encloses_pointwise_samples() -> None:
    s_box = _box(Fraction(0), Fraction(1))
    z_box = _box(Fraction(-1, 2), Fraction(3, 4))
    enclosure = gaussian_box(1, Fraction(1, 3), s_box, z_box, bits=48)
    for s in np.linspace(0.0, 1.0, 9):
        for z in np.linspace(-0.5, 0.75, 9):
            sample = math.exp(-(s + z * z / 3.0))
            assert float(enclosure.lower) - 1e-12 <= sample <= float(enclosure.upper) + 1e-12


def test_gaussian_box_rejects_bad_arguments() -> None:
    with pytest.raises(ValueError):
        gaussian_box(0, 1, _box(0, 1), _box(0, 1))
    with pytest.raises(ValueError):
        gaussian_box(1, 1, _box(-1, 1), _box(0, 1))


def test_interval_poly_eval_encloses_dense_samples() -> None:
    """Mixed signs and a straddling z-box exercise both the parity-sharp even
    powers and the cancellation the Horner intersection can pick up."""
    poly = {
        (0, 0): Fraction(1), (1, 0): Fraction(-2), (0, 2): Fraction(3),
        (2, 1): Fraction(1, 2), (1, 3): Fraction(-3, 4),
    }
    s_box = _box(Fraction(1, 5), Fraction(1, 2))
    z_box = _box(Fraction(-2, 5), Fraction(3, 10))
    enclosure = interval_poly_eval(poly, s_box, z_box, bits=48)
    for s in np.linspace(0.2, 0.5, 13):
        for z in np.linspace(-0.4, 0.3, 13):
            sample = 1.0 - 2.0 * s + 3.0 * z**2 + 0.5 * s**2 * z - 0.75 * s * z**3
            assert float(enclosure.lower) - 1e-12 <= sample <= float(enclosure.upper) + 1e-12
    with pytest.raises(ValueError):
        interval_poly_eval(poly, _box(Fraction(-1), Fraction(1)), z_box)


def test_gaussian_term_rejects_a_basis_without_decay() -> None:
    """A term with a non-positive rate would make every exterior tail bound
    silently meaningless, so it must refuse to exist at all."""
    with pytest.raises(ValueError):
        GaussianTerm.from_triple({(0, 1): 1}, 0, 1)
    with pytest.raises(ValueError):
        GaussianTerm.from_triple({(0, 1): 1}, 1, -1)
    with pytest.raises(ValueError):
        GaussianTerm.from_triple({(-1, 0): 1}, 1, 1)


def test_adapter_accepts_triples_and_attribute_objects() -> None:
    """The wiring contract with the symbolic gaussian_hermite module: a plain
    (poly, alpha, beta) triple and any object with those attributes must both
    land on the same term, so the two modules never need each other's names."""
    triple = ({(1, 2): Fraction(5, 3)}, Fraction(2), Fraction(1, 2))
    from_triple = as_gaussian_term(triple)
    assert from_triple == GaussianTerm.from_triple(*triple)
    assert as_gaussian_term(from_triple) is from_triple

    class Foreign:
        poly = {(1, 2): Fraction(5, 3)}
        alpha = Fraction(2)
        beta = Fraction(1, 2)

    assert as_gaussian_term(Foreign()) == from_triple

    family = seed_gaussian_family()
    as_dict = {"swirl": list(family.swirl), "stream": list(family.stream)}
    r_box, z_box = _box(Fraction(1, 2), Fraction(3, 5)), _box(Fraction(1, 5), Fraction(3, 10))
    direct = cell_enclosure_gaussian(family, r_box, z_box, bits=48)
    adapted = cell_enclosure_gaussian(as_dict, r_box, z_box, bits=48)
    assert direct.divergence.as_pair() == adapted.divergence.as_pair()


# --------------------------------------------------------------------------- #
# the cell enclosure                                                           #
# --------------------------------------------------------------------------- #


def test_divergence_enclosure_contains_zero_and_is_narrow() -> None:
    """The acceptance test that the Gaussian basis kills the flat-bump
    catastrophe.  The divergence is exactly zero by the stream-function
    identity, so its interval must contain zero and its width measures nothing
    but enclosure quality.  Measured on the seed family at bits=48, no cell
    splitting (subdivisions = 1 semantics), on the mid-support cell
    ``[1/2, 1/2 + 2^-12]^2``: divergence = [-1.259e-3, +1.259e-3], width
    2.517e-3, RMS gradient bound sqrt(gradient_squared.upper) = 3.980, ratio
    6.32e-4.  The flat-bump enclosures gave +-1.45e2 around zero -- five
    orders of magnitude wider -- because chi'' carries (1-sigma)^{-4}; the
    Gaussian closure has no such factor and the residual width above is the
    plain dependency loss of a product over a cell, linear in the cell size."""
    clear_gaussian_cache()
    family = seed_gaussian_family()
    h = Fraction(1, 4096)
    r_box = Interval(Fraction(1, 2), Fraction(1, 2) + h)
    z_box = Interval(Fraction(1, 2), Fraction(1, 2) + h)
    enclosure = cell_enclosure_gaussian(family, r_box, z_box, terms=32, bits=48)
    assert enclosure.divergence.contains_zero
    width = float(enclosure.divergence.upper - enclosure.divergence.lower)
    gradient_scale = math.sqrt(float(enclosure.gradient_squared.upper))
    assert width < 1.0e-3 * gradient_scale
    assert width < 3.0e-3  # the absolute headline: the flat bump gave ~2.9e2


def test_cell_enclosure_contains_dense_float_samples() -> None:
    """Speed, viscous integrand and the flux magnitude bound must contain a
    dense sample of the true field values, computed here from difference
    quotients of the closed form rather than from the module's own derivative
    polynomials (see the header note on independence)."""
    family = seed_gaussian_family()
    r_box = _box(Fraction(2, 5), Fraction(1, 2))
    z_box = _box(Fraction(1, 5), Fraction(3, 10))
    enclosure = cell_enclosure_gaussian(family, r_box, z_box, terms=32, bits=48)
    tolerance = 1.0e-4
    for r in np.linspace(0.4, 0.5, 7):
        for z in np.linspace(0.2, 0.3, 7):
            s = r * r
            swirl = _fd_partials(family.swirl, s, z)
            stream = _fd_partials(family.stream, s, z)
            u_t = r * swirl["value"]
            u_r = -r * stream["dz"]
            u_z = 2.0 * stream["value"] + 2.0 * s * stream["ds"]
            gradient = (
                -(stream["dz"] + 2.0 * s * stream["dsz"]),
                swirl["value"] + 2.0 * s * swirl["ds"],
                2.0 * r * (4.0 * stream["ds"] + 2.0 * s * stream["dss"]),
                -swirl["value"],
                -stream["dz"],
                0.0,
                -r * stream["dzz"],
                r * swirl["dz"],
                2.0 * stream["dz"] + 2.0 * s * stream["dsz"],
            )
            speed = math.sqrt(u_r**2 + u_t**2 + u_z**2)
            gradient_sq = sum(g * g for g in gradient)
            assert float(enclosure.speed.lower) - tolerance <= speed
            assert speed <= float(enclosure.speed.upper) + tolerance
            viscous = 2.0 * speed * gradient_sq
            assert float(enclosure.viscous_integrand.lower) - tolerance <= viscous
            assert viscous <= float(enclosure.viscous_integrand.upper) + tolerance
            kato = speed * math.sqrt(gradient_sq)
            assert kato <= float(enclosure.flux_magnitude.upper) + tolerance


# --------------------------------------------------------------------------- #
# the tail bounds                                                              #
# --------------------------------------------------------------------------- #


def test_gaussian_tail_moment_dominates_a_numerical_integral() -> None:
    """The bound must dominate a fine trapezoid quadrature and must not be
    grossly loose -- the m=0 Mills ratio is the worst case and is still within
    19% here.  The m=1 descent base is EXACT, and the trapezoid rule slightly
    over-estimates a convex integrand, so domination is asserted with a 1e-6
    relative allowance for the quadrature's own bias; everything the bound
    could get wrong is far larger than that."""
    for gamma, R, m in (
        (1.0, 1.0, 0), (1.0, 1.0, 1), (0.5, 1.5, 2), (2.0, 1.0, 4), (0.75, 2.0, 3),
    ):
        bound = gaussian_tail_moment(Fraction(gamma), Fraction(R), m, bits=64)
        grid = np.linspace(R, R + 12.0 / math.sqrt(gamma), 200001)
        integrand = grid**m * np.exp(-2.0 * gamma * grid**2)
        numeric = float(
            np.sum((integrand[1:] + integrand[:-1]) * 0.5 * (grid[1] - grid[0]))
        )
        assert numeric <= float(bound) * (1.0 + 1.0e-6)
        assert float(bound) <= 1.2 * numeric + 1.0e-30


def test_gaussian_tail_moment_decays_and_rejects_bad_arguments() -> None:
    small = gaussian_tail_moment(1, 3, 2, bits=64)
    smaller = gaussian_tail_moment(1, 4, 2, bits=64)
    assert Fraction(0) < smaller < small
    with pytest.raises(ValueError):
        gaussian_tail_moment(0, 1, 0)
    with pytest.raises(ValueError):
        gaussian_tail_moment(1, 0, 0)  # the Mills bound genuinely needs R > 0
    with pytest.raises(ValueError):
        gaussian_tail_moment(1, 1, -1)


def test_exterior_tail_is_finite_positive_and_decays_with_the_box() -> None:
    """The tails must be finite (the flat-bump basis got this for free from
    compact support; the Gaussian basis has to earn it) and must fall like
    exp(-3 gamma d^2) as the box grows, or the certificate could never make
    the exterior contribution negligible by enlarging the domain."""
    family = seed_gaussian_family()
    near = exterior_viscous_tail(family, 3, 3, bits=64)
    far = exterior_viscous_tail(family, 4, 4, bits=64)
    for tail in (near, far):
        assert tail.viscous_tail > 0
        assert tail.flux_tail > 0
        assert tail.cubed_flux_tail > 0
        assert tail.speed_coefficient > 0 and tail.gradient_coefficient > 0
    assert far.viscous_tail < near.viscous_tail
    assert far.flux_tail < near.flux_tail
    # e^{-3 gamma (16 - 9)} with gamma = 1/2 is about 2.8e-5; the measured
    # ratio 3.6e-5 confirms the advertised decay rate, not merely decrease.
    assert far.viscous_tail < near.viscous_tail * Fraction(1, 10000)


def test_exterior_tail_envelope_dominates_the_field_on_the_boundary() -> None:
    """|u| <= C_u e^{-gamma rho^2} is claimed on ALL of space; spot-check it
    where it matters most, along the boundary of the certificate box."""
    family = seed_gaussian_family()
    tail = exterior_viscous_tail(family, 3, 3, bits=64)
    gamma = float(tail.gamma)
    for r, z in ((3.0, 0.0), (3.0, 1.5), (0.0, 3.0), (2.0, 3.0), (3.0, 3.0)):
        s = r * r
        swirl = _fd_partials(family.swirl, s, z)
        stream = _fd_partials(family.stream, s, z)
        u_t = r * swirl["value"]
        u_r = -r * stream["dz"]
        u_z = 2.0 * stream["value"] + 2.0 * s * stream["ds"]
        speed = math.sqrt(u_r**2 + u_t**2 + u_z**2)
        envelope = float(tail.speed_coefficient) * math.exp(-gamma * (s + z * z))
        assert speed <= envelope * (1.0 + 1.0e-9)


# --------------------------------------------------------------------------- #
# the end-to-end certificate                                                   #
# --------------------------------------------------------------------------- #


def _grid() -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=9, nz=17, r_max=3.0, z_min=-3.0, z_max=3.0, periodic_z=False
    )


def _certificate():
    return build_gaussian_generation_certificate(
        seed_gaussian_family(),
        _grid(),
        viscosity=5.0e-3,
        precision_bits=48,
        exp_terms=28,
        subdivisions=1,
    )


def test_certificate_payload_verifies_and_is_internally_consistent() -> None:
    """The verifier is the SAME one the flat-bump certificate uses -- the
    checker must not care which basis produced the enclosures -- and the tail
    accounting must reconcile exactly: interior + exterior = total, as
    Fractions, not as floats."""
    payload = _certificate().as_dict()
    report = verify_generation_certificate(payload)
    assert report["verified"] is True, report["failures"]

    rigorous = payload["rigorous"]
    assert rigorous["enclosure_basis"] == "gaussian_hermite"
    margin = Fraction(rigorous["margin"])
    assert margin == (
        Fraction(rigorous["pressure_lower_bound_discrete"])
        + Fraction(rigorous["viscous_lower_bound"])
    )
    total_viscous = Fraction(rigorous["viscous_integral_upper_bound"])
    interior = Fraction(rigorous["interior_viscous_integral_upper_bound"])
    tail = Fraction(rigorous["exterior_viscous_tail"])
    assert total_viscous == interior + tail
    assert tail > 0
    flux_total = Fraction(rigorous["flux_integral_upper_bound"])
    flux_interior = Fraction(rigorous["interior_flux_integral_upper_bound"])
    flux_tail = Fraction(rigorous["exterior_flux_tail"])
    assert flux_total == 3 * (flux_interior + flux_tail)
    # Every cell's divergence interval contained zero, so the recorded largest
    # excursion is exactly the largest half-width the sweep saw.
    assert Fraction(rigorous["divergence_excursion"]) > 0
    assert any("Gaussian-Hermite" in str(claim) for claim in payload["claims"])
    assert any("Clay" in str(claim) for claim in payload["claims"])


def test_certificate_tamper_rejections_still_fire() -> None:
    """The Gaussian payload must fail the shared verifier for the same reasons
    the flat-bump payload would: a certificate whose checker cannot reject a
    forged margin is not a certificate."""
    payload = _certificate().as_dict()

    forged = copy.deepcopy(payload)
    forged["rigorous"]["margin"] = str(Fraction(forged["rigorous"]["margin"]) + 1)
    assert verify_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["hypotheses"]["P1_pressure_discretisation"]["statement"] = "trust me"
    assert verify_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["hypotheses"]["P1_pressure_discretisation"]["proved"] = True
    assert verify_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["claims"] = [c for c in forged["claims"] if "Clay" not in str(c)]
    assert verify_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["pressure_error_threshold"] = str(
        Fraction(forged["rigorous"]["pressure_error_threshold"]) + 1
    )
    assert verify_generation_certificate(forged)["verified"] is False


def test_certificate_builder_rejects_bad_inputs() -> None:
    family = seed_gaussian_family()
    with pytest.raises(ValueError):
        build_gaussian_generation_certificate(family, _grid(), viscosity=0.0)
    with pytest.raises(ValueError):
        build_gaussian_generation_certificate(
            family, _grid(), viscosity=1.0e-2, subdivisions=0
        )
    with pytest.raises(ValueError):
        build_gaussian_generation_certificate(
            family,
            AxisymmetricGrid.uniform(nr=9, nz=16, r_max=3.0, z_min=-3.0, z_max=3.0),
            viscosity=1.0e-2,
        )
    with pytest.raises(ValueError):
        # The exterior tail needs a box surrounding the origin in z.
        build_gaussian_generation_certificate(
            family,
            AxisymmetricGrid.uniform(
                nr=9, nz=17, r_max=3.0, z_min=0.5, z_max=3.0, periodic_z=False
            ),
            viscosity=1.0e-2,
        )
    surrogate = GaussianFamily(
        name="surrogate",
        swirl=seed_gaussian_family().swirl,
        stream=(),
        clay_admissible=False,
    )
    with pytest.raises(ValueError):
        build_gaussian_generation_certificate(surrogate, _grid(), viscosity=1.0e-2)
