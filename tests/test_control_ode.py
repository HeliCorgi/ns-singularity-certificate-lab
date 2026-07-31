"""Tests for the rigorous scalar control ODE: Picard box, Bernstein, Taylor."""

from __future__ import annotations

import json
from fractions import Fraction

import pytest

from ns_certificate_lab.control_ode import (
    BERNSTEIN_HULL_PROPERTY,
    CONTROL_LIMITATIONS,
    ControlCoefficients,
    PICARD_SELF_MAPPING,
    bernstein_coefficients,
    build_control_certificate,
    prove_rough_enclosure,
    taylor_coefficients,
    taylor_tube,
    verify_control_certificate,
)
from ns_certificate_lab.l3_certificate import exp_interval
from ns_certificate_lab.snapshot_certificate import Interval


def _thin(value: Fraction | int) -> Interval:
    return Interval(Fraction(value), Fraction(value))


def _span(lower: Fraction | int, upper: Fraction | int) -> Interval:
    return Interval(Fraction(lower), Fraction(upper))


def _coefficients(linear, quadratic, forcing) -> ControlCoefficients:
    return ControlCoefficients(_thin(linear), _thin(quadratic), _thin(forcing))


#: A decaying problem with a mild quadratic term: the sort of control ODE a
#: viscosity-dominated slab produces.
GENTLE = _coefficients(Fraction(-1, 2), Fraction(1, 10), Fraction(1, 100))

#: The Riccati equation ``R' = R^2`` from ``R_0 = 1``, which blows up at
#: ``t = 1``.  Asked to certify a slab of length two it cannot succeed, and
#: that is the point of using it.
STIFF = _coefficients(0, 1, 0)

#: The whole-space shape of the control ODE: a positive linear term because the
#: ``-nu R`` damping of the torus theory is absent, a quadratic term, and a
#: small residual forcing.
REALISTIC = ControlCoefficients(_thin(2), _thin(1), _span(0, Fraction(1, 10**6)))


# --------------------------------------------------------------------------- #
# the interval extension of the right-hand side                                #
# --------------------------------------------------------------------------- #


def test_the_quadratic_term_uses_a_square_and_not_a_self_product() -> None:
    """``[-1,1] * [-1,1]`` is ``[-1,1]``; the range of ``x^2`` there is ``[0,1]``.

    The difference is not cosmetic.  A control ODE whose state interval
    straddles zero would otherwise be handed a spurious negative quadratic
    contribution, which on the decaying side is an unearned bound and on the
    growing side an unnecessary loss.
    """
    coefficients = _coefficients(0, 1, 0)
    field = coefficients.field(_span(-1, 1))
    assert (field.lower, field.upper) == (Fraction(0), Fraction(1))


def test_the_lipschitz_bound_is_the_polynomial_one() -> None:
    coefficients = _coefficients(Fraction(-1, 2), 3, 7)
    box = _span(-2, 5)
    # |a| + 2|b| max|B| = 1/2 + 2*3*5.
    assert coefficients.lipschitz_bound(box) == Fraction(61, 2)


# --------------------------------------------------------------------------- #
# stage B: the self-mapping proof                                              #
# --------------------------------------------------------------------------- #


def test_self_mapping_succeeds_immediately_on_a_gentle_problem() -> None:
    """A decaying slab should close on the first box the heuristic proposes.

    If it needed doublings here the starting radius ``2h|F(R_0)|`` would be
    mis-scaled, and every later problem would pay for it in retries.
    """
    rough = prove_rough_enclosure(
        GENTLE, _span(Fraction(9, 10), Fraction(11, 10)), Fraction(1, 20)
    )
    assert rough.proved
    assert rough.doublings == 0
    assert rough.halvings == 0
    assert rough.step == Fraction(1, 20)
    # The proved inclusion, re-checked here rather than trusted.
    assert rough.box.lower <= rough.image.lower
    assert rough.image.upper <= rough.box.upper
    assert rough.box.lower <= rough.initial.lower
    assert rough.initial.upper <= rough.box.upper


def test_self_mapping_reports_failure_on_a_stiff_slab_without_inflating() -> None:
    """``R' = R^2`` from ``R_0 = 1`` blows up at ``t = 1``, so a slab of length
    two admits no self-mapping box at all: the inclusion needs ``4h <= 1``.

    The requirement is not merely that the function returns ``proved = False``.
    It must also not hand back a wide box with a quiet flag, because a box that
    fails the inclusion proves nothing however wide it is.  So the returned box
    is checked to genuinely fail the test it was supposed to pass.
    """
    rough = prove_rough_enclosure(STIFF, _thin(1), Fraction(2), max_halvings=0)
    assert rough.proved is False
    assert rough.halvings == 0
    assert rough.step == Fraction(2)
    escapes = (
        rough.image.upper > rough.box.upper or rough.image.lower < rough.box.lower
    )
    assert escapes
    assert "failure" in rough.verdict


def test_halving_the_slab_rescues_a_failed_self_mapping() -> None:
    """The same problem, with the outer budget released.

    The inclusion ``1 + h(1+rho)^2 <= 1 + rho`` is feasible exactly when
    ``4h <= 1``, and at ``h = 1/4`` only for the single radius ``rho = 1``,
    which a doubling search cannot be expected to land on.  The first length
    with an open interval of feasible radii is ``h = 1/8``, so that is where
    the search must stop: two, one, a half and a quarter all fail, four
    halvings in all.
    """
    rough = prove_rough_enclosure(STIFF, _thin(1), Fraction(2), max_halvings=6)
    assert rough.proved
    assert rough.step == Fraction(1, 8)
    assert rough.halvings == 4
    assert rough.requested_step == Fraction(2)


def test_self_mapping_rejects_a_nonpositive_slab_and_a_zero_floor() -> None:
    with pytest.raises(ValueError):
        prove_rough_enclosure(GENTLE, _thin(1), Fraction(0))
    with pytest.raises(ValueError):
        prove_rough_enclosure(
            GENTLE, _thin(1), Fraction(1, 10), inflation_floor=Fraction(0)
        )


def test_a_stationary_datum_still_gets_a_box() -> None:
    """``F([R_0])`` vanishes at ``R_0 = 0`` for a forcing-free problem, so the
    leading term of the starting radius is zero and only the floor keeps the
    search alive.  Without it the box would be degenerate and could never be
    doubled anywhere."""
    rough = prove_rough_enclosure(_coefficients(1, 1, 0), _thin(0), Fraction(1, 10))
    assert rough.proved
    assert rough.box.lower < 0 < rough.box.upper


def test_a_tube_refuses_to_be_built_without_a_proved_box() -> None:
    rough = prove_rough_enclosure(STIFF, _thin(1), Fraction(2), max_halvings=0)
    with pytest.raises(ValueError):
        taylor_tube(STIFF, rough)


# --------------------------------------------------------------------------- #
# the Bernstein form                                                           #
# --------------------------------------------------------------------------- #


def _evaluate(monomial, point: Fraction) -> Fraction:
    return sum(
        (item.lower * point**k for k, item in enumerate(monomial)), Fraction(0)
    )


def test_bernstein_hull_is_exact_when_the_extrema_are_at_the_endpoints() -> None:
    """``s^2`` on ``[0,2]`` is monotone, so its range is the hull of its end
    values and the Bernstein form must reproduce that with nothing to spare.

    A conversion that lost tightness on a monomial this simple would be losing
    it everywhere, and the loss would be invisible in the certificate because
    an over-wide enclosure is still a valid one.
    """
    monomial = [_thin(0), _thin(0), _thin(1)]
    coefficients = bernstein_coefficients(monomial, Fraction(2))
    hull = Interval.hull(coefficients)
    assert (hull.lower, hull.upper) == (Fraction(0), Fraction(4))


def test_bernstein_hull_strictly_contains_a_range_with_an_interior_extremum()\
        -> None:
    """``s^2 - s`` on ``[0,1]`` has its minimum ``-1/4`` at ``s = 1/2``.

    The convex hull property gives a *valid* enclosure, not a tight one: the
    control polygon dips to ``-1/2``, twice as far.  This is exactly the price
    of not using a derivative bound, and it is recorded here so that a later
    tightening (subdivision, degree elevation) can be measured against it.
    """
    monomial = [_thin(0), _thin(-1), _thin(1)]
    coefficients = bernstein_coefficients(monomial, Fraction(1))
    hull = Interval.hull(coefficients)
    assert hull.lower == Fraction(-1, 2)
    assert hull.upper == Fraction(0)
    assert hull.lower < Fraction(-1, 4)
    # Valid on the whole interval, sampled densely enough to catch a slip.
    for index in range(0, 1001):
        point = Fraction(index, 1000)
        assert hull.lower <= _evaluate(monomial, point) <= hull.upper


def test_bernstein_conversion_reproduces_the_polynomial_at_both_endpoints() -> None:
    """``b_0 = c_0`` and ``b_p`` is the polynomial at ``s = h``; those two are
    interpolated exactly by the Bezier form, so they anchor the conversion."""
    monomial = [_thin(2), _thin(-3), _thin(5), _thin(-1)]
    step = Fraction(3, 7)
    coefficients = bernstein_coefficients(monomial, step)
    assert coefficients[0].contains(_evaluate(monomial, Fraction(0)))
    assert coefficients[-1].contains(_evaluate(monomial, step))


# --------------------------------------------------------------------------- #
# stage C: the Taylor tube against exactly known solutions                     #
# --------------------------------------------------------------------------- #


def test_tube_is_exact_on_the_affine_solution() -> None:
    """``a = b = 0`` and ``e`` constant gives ``R(t) = R_0 + e t`` exactly.

    Every Taylor coefficient above the first vanishes, the remainder vanishes
    identically, and the Bernstein hull of a degree-one polynomial is its range.
    The tube must therefore be the exact range ``[R_0, R_0 + e h]``, with no
    slack at all -- an equality, not a containment, because any slack here would
    be pure interval overestimation with nothing to blame it on.
    """
    coefficients = _coefficients(0, 0, 3)
    certificate = build_control_certificate(
        coefficients, _thin(2), Fraction(1, 4), order=5
    )
    block = certificate.as_dict()["tube"]
    tube = Interval.from_pair(block["tube"])
    assert (tube.lower, tube.upper) == (Fraction(2), Fraction(11, 4))
    remainder = Interval.from_pair(block["remainder"])
    assert (remainder.lower, remainder.upper) == (Fraction(0), Fraction(0))
    for index in range(0, 21):
        point = Fraction(index, 20) * Fraction(1, 4)
        assert tube.contains(Fraction(2) + 3 * point)


def test_tube_contains_the_exponential_solution_rigorously() -> None:
    """``b = e = 0`` gives ``R(t) = R_0 e^{a t}``, and the comparison is made
    against a rigorous rational enclosure of the exponential rather than a
    binary64 value, so a failure here cannot be a rounding artefact.

    The tube is also required to be tight: the truncated series at order eight
    over a slab of length one half leaves a remainder around ``10^{-11}``, so a
    tube wider than the true range by more than a part in ``10^6`` would mean
    the interval evaluation, not the truncation, is dominating.
    """
    coefficients = _coefficients(Fraction(1, 2), 0, 0)
    step = Fraction(1, 2)
    certificate = build_control_certificate(
        coefficients, _thin(1), step, order=8
    )
    tube = Interval.from_pair(certificate.as_dict()["tube"]["tube"])
    for index in range(0, 11):
        point = Fraction(index, 10) * step
        exact = exp_interval(Interval.exact(point / 2))
        assert tube.lower <= exact.lower and exact.upper <= tube.upper
    true_range = exp_interval(Interval(Fraction(0), Fraction(1, 4)))
    slack = (tube.upper - tube.lower) - (true_range.upper - true_range.lower)
    assert slack >= 0
    assert slack < Fraction(1, 10**6)


def test_taylor_coefficients_reproduce_the_exponential_series() -> None:
    """``a_k = R_0 a^k / k!`` for the linear problem; the recurrence has to
    generate exactly that, so a slip in the ``1/(k+1)`` factor shows here."""
    coefficients = _coefficients(Fraction(1, 2), 0, 0)
    series = taylor_coefficients(coefficients, _thin(1), 6)
    factorial = 1
    for k, item in enumerate(series):
        if k:
            factorial *= k
        assert item.contains(Fraction(1, 2**k * factorial))
        assert item.upper - item.lower < Fraction(1, 2**60)


def test_taylor_coefficients_reproduce_the_riccati_series() -> None:
    """``R' = R^2`` from ``R_0 = 1`` has ``R(t) = 1/(1-t)``, whose coefficients
    are all one.  This is the case that exercises the convolution, including
    the diagonal square term."""
    series = taylor_coefficients(STIFF, _thin(1), 6)
    for item in series:
        assert item.contains(Fraction(1))
        assert item.upper - item.lower < Fraction(1, 2**60)


def test_the_tube_never_escapes_the_rough_box_it_was_built_on() -> None:
    """Both enclose the same solution, so their intersection does too; the
    refined tube is what a downstream step should quote."""
    certificate = build_control_certificate(
        REALISTIC, _thin(10), Fraction(1, 50), order=8
    )
    payload = certificate.as_dict()
    box = Interval.from_pair(payload["self_mapping"]["box"])
    refined = Interval.from_pair(payload["tube"]["refined_tube"])
    tube = Interval.from_pair(payload["tube"]["tube"])
    assert box.lower <= refined.lower and refined.upper <= box.upper
    assert tube.lower <= refined.lower and refined.upper <= tube.upper
    # The tube is much the tighter of the two on a slab this short.
    assert payload["tube"]["tightness_ratio"] < 0.25


# --------------------------------------------------------------------------- #
# the certificate and its independent checker                                  #
# --------------------------------------------------------------------------- #


def _payload(**kwargs) -> dict:
    certificate = build_control_certificate(
        GENTLE, _span(Fraction(9, 10), Fraction(11, 10)), Fraction(1, 20),
        order=6, **kwargs,
    )
    return json.loads(json.dumps(certificate.as_dict(), allow_nan=False))


def test_certificate_verifies_and_proves_existence() -> None:
    verdict = verify_control_certificate(_payload())
    assert verdict["verified"] is True
    assert verdict["proves_existence"] is True
    assert verdict["failures"] == []


def test_certificate_records_the_audited_theorem_statements() -> None:
    payload = _payload()
    assert payload["self_mapping"]["statement"] == PICARD_SELF_MAPPING
    assert payload["tube"]["bernstein_property"] == BERNSTEIN_HULL_PROPERTY


def test_certificate_records_the_three_open_limitations() -> None:
    """None of the three is closed by this module, and the payload has to say
    so in its own text or a reader could take the tube for a PDE statement."""
    payload = _payload()
    assert set(payload["limitations"]) == set(CONTROL_LIMITATIONS)
    assert all(entry["closed"] is False for entry in payload["limitations"].values())
    absent = payload["limitations"]["R3_has_no_spectral_gap"]["statement"]
    assert "torus" in absent and "-nu R" in absent
    assert "HS-5" in payload["limitations"]["HS5_coefficients_to_pde"]["statement"]


def test_a_failed_slab_still_produces_a_consistent_certificate() -> None:
    """A failure has to be recordable.  The payload is internally consistent --
    it verifies -- while proving nothing, and the two facts are reported
    separately so that neither can be read as the other."""
    certificate = build_control_certificate(
        STIFF, _thin(1), Fraction(2), max_halvings=0
    )
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    assert payload["self_mapping"]["proved"] is False
    assert payload["tube"] is None
    assert payload["picard"] is None
    verdict = verify_control_certificate(payload)
    assert verdict["verified"] is True
    assert verdict["proves_existence"] is False


def test_checker_rejects_an_inverted_enclosure() -> None:
    payload = _payload()
    lower, upper = payload["tube"]["tube"]
    payload["tube"]["tube"] = [upper, lower]
    verdict = verify_control_certificate(payload)
    assert verdict["verified"] is False
    assert verdict["proves_existence"] is False
    assert verdict["failures"]


def test_checker_rejects_a_self_mapping_claim_its_own_field_contradicts() -> None:
    """Widening ``F(B)`` while leaving the flag and the image alone is the
    natural forgery: the payload still says the box is self-mapping, but its
    own recorded field no longer sends the box into itself."""
    payload = _payload()
    lower, _upper = payload["self_mapping"]["field_on_box"]
    payload["self_mapping"]["field_on_box"] = [lower, "100"]
    verdict = verify_control_certificate(payload)
    assert verdict["verified"] is False
    assert any("self-mapping inclusion" in failure for failure in verdict["failures"])


def test_checker_rejects_a_narrowed_bernstein_enclosure() -> None:
    """The convex hull property is an inequality in one direction only.  An
    enclosure tighter than the hull of its own coefficients is not a sharper
    result, it is an unsupported one."""
    payload = _payload()
    coefficients = [
        Interval.from_pair(pair)
        for pair in payload["tube"]["bernstein_coefficients"]
    ]
    hull = Interval.hull(coefficients)
    shrunk = Fraction(hull.lower + hull.upper, 2)
    payload["tube"]["polynomial_enclosure"] = [str(shrunk), str(hull.upper)]
    verdict = verify_control_certificate(payload)
    assert verdict["verified"] is False
    assert any("convex hull" in failure for failure in verdict["failures"])


def test_checker_rejects_a_dropped_clay_disclaimer() -> None:
    payload = _payload()
    payload["claims"] = ["the control ODE closes"]
    verdict = verify_control_certificate(payload)
    assert verdict["verified"] is False
    assert any("Clay" in failure for failure in verdict["failures"])


def test_checker_rejects_a_narrowed_remainder_coefficient() -> None:
    """The remainder is the only place the rough box enters the tube, so a
    payload that shrinks it is quietly discarding the one hypothesis the tube
    depends on."""
    payload = _payload()
    lower, upper = payload["tube"]["remainder_coefficient"]
    middle = Fraction(Fraction(lower) + Fraction(upper), 2)
    payload["tube"]["remainder_coefficient"] = [str(middle), str(middle)]
    assert verify_control_certificate(payload)["verified"] is False


def test_checker_rejects_an_edited_theorem_statement() -> None:
    payload = _payload()
    payload["self_mapping"]["statement"] = "the box is fine"
    assert verify_control_certificate(payload)["verified"] is False


def test_checker_rejects_a_promoted_limitation() -> None:
    payload = _payload()
    payload["limitations"]["HS5_coefficients_to_pde"]["closed"] = True
    verdict = verify_control_certificate(payload)
    assert verdict["verified"] is False
    assert any("closed" in failure for failure in verdict["failures"])


def test_checker_rejects_a_step_that_disagrees_with_the_halvings() -> None:
    payload = _payload()
    payload["slab"]["halvings"] = 3
    assert verify_control_certificate(payload)["verified"] is False


def test_checker_returns_a_verdict_on_garbage_rather_than_raising() -> None:
    for payload in ({}, {"certificate": "scalar_control_ode_slab"}, {"slab": 3},
                    {"certificate": "scalar_control_ode_slab", "coefficients": 7}):
        verdict = verify_control_certificate(payload)
        assert verdict["verified"] is False
        assert verdict["proves_existence"] is False
        assert verdict["failures"]


# --------------------------------------------------------------------------- #
# the whole-space shape of the problem                                         #
# --------------------------------------------------------------------------- #


def test_the_certified_slab_is_a_fraction_of_the_blow_up_time() -> None:
    """``R' = 2R + R^2`` from ``R_0 = 10`` blows up at
    ``t* = (1/2) log(6/5) = 0.0912``, and the self-mapping test closes up to
    about ``h = 0.023``.

    That ratio is the honest cost of the construction and is worth pinning: a
    quarter of the blow-up time per slab, with no linear damping available to
    extend it, because on ``R^3`` the ``-nu R`` term of the torus control ODE
    is not there to be used.
    """
    assert prove_rough_enclosure(
        REALISTIC, _thin(10), Fraction(1, 40), max_halvings=0
    ).proved is False
    closing = prove_rough_enclosure(
        REALISTIC, _thin(10), Fraction(1, 50), max_halvings=0
    )
    assert closing.proved
    assert closing.doublings <= 2
