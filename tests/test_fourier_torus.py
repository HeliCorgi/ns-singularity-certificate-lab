"""Tests for the Track-P exact Fourier machinery and the torus slab certificate."""

from __future__ import annotations

from fractions import Fraction
import json
import math

import pytest

from ns_certificate_lab.fourier_torus import (
    TORUS_FAMILIES,
    TrigVector,
    adot_squared_upper,
    advection,
    family_P1,
    family_P2,
    family_P3,
    galerkin_modes,
    galerkin_rhs,
    gradient_part,
    leray,
    residual_tail,
    sup_derivative_bound,
    taylor_coefficients,
)
from ns_certificate_lab.torus_aposteriori import (
    EXTERNAL_THEOREMS,
    SQRT8_UPPER,
    build_torus_slab_certificate,
    prove_galerkin_box,
    verify_torus_slab_certificate,
)


# --------------------------------------------------------------------------- #
# exact machinery                                                              #
# --------------------------------------------------------------------------- #


def test_families_are_exactly_divergence_free_and_mean_zero() -> None:
    """Divergence freedom is per-mode orthogonality, so a failure here means a
    coefficient vector was typed against the wrong wavevector."""
    for name, builder in TORUS_FAMILIES.items():
        field = builder()
        assert not field.divergence().cleaned().terms, name
        assert (0, 0, 0) not in field.coefficient_table(), name


def test_families_are_inside_the_default_galerkin_band() -> None:
    for builder in TORUS_FAMILIES.values():
        for k in builder().coefficient_table():
            assert 1 <= sum(c * c for c in k) <= 4


def test_product_matches_pointwise_evaluation() -> None:
    """The product-to-sum table is where a sign error would live.  Sampling the
    exact product against the float product of the factors at irrational points
    catches every entry of the table at once."""
    u = family_P3()
    f = u.components[0]
    g = u.components[1]
    product = f * g
    for point in [(0.3, 1.1, -0.7), (2.1, 0.5, 0.9), (-1.3, 2.2, 0.4)]:
        assert product.evaluate(point) == pytest.approx(
            f.evaluate(point) * g.evaluate(point), abs=1e-12
        )


def test_derivative_matches_finite_differences() -> None:
    u = family_P2()
    scalar = u.components[2]
    step = 1.0e-6
    for direction in range(3):
        derivative = scalar.derivative(direction)
        base = [0.4, -0.8, 1.7]
        forward = list(base)
        backward = list(base)
        forward[direction] += step
        backward[direction] -= step
        numeric = (
            scalar.evaluate(tuple(forward)) - scalar.evaluate(tuple(backward))
        ) / (2.0 * step)
        assert derivative.evaluate(tuple(base)) == pytest.approx(numeric, abs=1e-6)


def test_advection_mean_vanishes_exactly() -> None:
    """The mean of ``(u . grad) u`` is zero for divergence-free ``u`` by parts.
    Exact arithmetic must reproduce that as an identity, not approximately —
    anything else means the convolution folded a mode onto ``k = 0``."""
    for builder in TORUS_FAMILIES.values():
        u = builder()
        mean = advection(u, u).coefficient_table().get((0, 0, 0))
        if mean is not None:
            assert all(value == 0 for vector in mean for value in vector)


def test_leray_output_is_divergence_free_and_idempotent() -> None:
    u = family_P3()
    w = advection(u, u)
    projected = leray(w)
    # Divergence in Fourier: k . a_k = 0 per mode.
    for k, (a, b) in projected.coefficient_table().items():
        if k == (0, 0, 0):
            continue
        assert sum(Fraction(k[i]) * a[i] for i in range(3)) == 0
        assert sum(Fraction(k[i]) * b[i] for i in range(3)) == 0
    twice = leray(projected)
    for k, (a, b) in projected.coefficient_table().items():
        again = twice.coefficient_table()[k]
        assert a == again[0] and b == again[1]


def test_leray_orthogonality_is_pythagorean() -> None:
    """``|w|^2 = |Pw|^2 + |(I-P)w|^2`` exactly — the projection is orthogonal,
    and in exact rationals the identity must hold to the last digit."""
    u = family_P2()
    w = advection(u, u)
    total = w.sobolev_sq(0)
    projected = leray(w).sobolev_sq(0)
    gradient = gradient_part(w).sobolev_sq(0)
    assert total == projected + gradient


def test_family_P3_is_pressure_driven() -> None:  # noqa: N802
    """P3's design property, measured: the part the pressure removes is larger
    than half the part that survives."""
    u = family_P3()
    w = advection(u, u)
    assert gradient_part(w).sobolev_sq(0) * 2 > leray(w).sobolev_sq(0)


def test_nonlinearity_is_nonzero_for_every_family() -> None:
    """A family whose nonlinear term vanished (a Beltrami flow, for instance)
    would make the whole a posteriori exercise vacuous: the Galerkin trajectory
    would be exact heat decay and the residual would be zero."""
    for name, builder in TORUS_FAMILIES.items():
        u = builder()
        assert leray(advection(u, u)).sobolev_sq(0) > 0, name


def test_residual_tail_lives_strictly_outside_the_band() -> None:
    u = family_P1()
    tail = residual_tail(u, cutoff_sq=4)
    for k in tail.coefficient_table():
        if k == (0, 0, 0):
            continue
        assert sum(c * c for c in k) > 4
    assert tail.sobolev_sq(4) > 0


def test_sobolev_norms_are_monotone_in_the_order() -> None:
    """``|k| >= 1`` on mean-zero fields makes the homogeneous ladder monotone;
    the control inequality's damping step depends on exactly this."""
    u = family_P2()
    values = [u.sobolev_sq(n) for n in range(5)]
    for lower, higher in zip(values, values[1:]):
        assert lower <= higher


def test_taylor_recursion_matches_a_small_time_step() -> None:
    """The Taylor coefficients must actually solve the Galerkin ODE: comparing
    the degree-3 polynomial against one tiny explicit Euler-refined step checks
    the recursion including the 1/(m+1) factor."""
    nu = Fraction(1, 10)
    u = family_P1()
    coefficients = taylor_coefficients(u, viscosity=nu, cutoff_sq=4, order=3)
    t = 1.0e-3
    # Fine RK4-free reference: 64 explicit Euler micro-steps in float.
    state = {
        (k, i, slot): float(pair[slot])
        for i in range(3)
        for k, pair in u.components[i].terms.items()
        for slot in range(2)
    }
    steps = 64
    dt = t / steps
    current = u
    for _ in range(steps):
        rhs = galerkin_rhs(current, viscosity=nu, cutoff_sq=4)
        # float Euler on top of exact structures, adequate at this tolerance
        new_components = []
        for i in range(3):
            terms = {}
            keys = set(current.components[i].terms) | set(rhs.components[i].terms)
            for k in keys:
                base = current.components[i].terms.get(k, [Fraction(0), Fraction(0)])
                slope = rhs.components[i].terms.get(k, [Fraction(0), Fraction(0)])
                terms[k] = [
                    base[0] + Fraction(dt) * slope[0],
                    base[1] + Fraction(dt) * slope[1],
                ]
            new_components.append(type(current.components[i])(terms, current.ops))
        current = TrigVector(tuple(new_components))
    taylor_value = {}
    for m, cm in enumerate(coefficients):
        for i in range(3):
            for k, pair in cm.components[i].terms.items():
                for slot in range(2):
                    taylor_value[(k, i, slot)] = taylor_value.get(
                        (k, i, slot), 0.0
                    ) + float(pair[slot]) * t**m
    for key, value in taylor_value.items():
        k, i, slot = key
        reference = float(
            current.components[i].terms.get(k, [Fraction(0), Fraction(0)])[slot]
        )
        assert value == pytest.approx(reference, abs=5e-5)


def test_adot_bound_dominates_a_partial_sum_and_tightens() -> None:
    """The bound must sit above every partial sum (it is an upper bound on the
    full lattice sum) and decrease as the exact region grows."""
    coarse = adot_squared_upper(6)
    fine = adot_squared_upper(20)
    partial = Fraction(0)
    for kx in range(-6, 7):
        for ky in range(-6, 7):
            for kz in range(-6, 7):
                n = kx * kx + ky * ky + kz * kz
                if n:
                    partial += Fraction(1, n * n)
    assert fine < coarse
    assert fine > partial
    with pytest.raises(ValueError):
        adot_squared_upper(0)


def test_sup_derivative_bound_dominates_sampled_derivatives() -> None:
    """M_1 must dominate |grad u| pointwise; sampling the actual gradient at
    scattered points is the only check that catches a dropped sqrt."""
    u = family_P1()
    bound = float(sup_derivative_bound(u, 1))
    step = 1.0e-5
    for point in [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0), (0.7, -1.1, 2.9)]:
        frobenius = 0.0
        for i in range(3):
            for j in range(3):
                forward = list(point)
                backward = list(point)
                forward[i] += step
                backward[i] -= step
                derivative = (
                    u.components[j].evaluate(tuple(forward))
                    - u.components[j].evaluate(tuple(backward))
                ) / (2.0 * step)
                frobenius += derivative**2
        assert math.sqrt(frobenius) <= bound + 1e-6


def test_galerkin_modes_counts() -> None:
    # |k|^2 = 1: 6 vectors, = 2: 12, = 3: 8, = 4: 6 -> 32 vectors, 16 canonical.
    assert len(galerkin_modes(4)) == 16
    assert len(galerkin_modes(1)) == 3


def test_from_modes_refuses_the_mean_mode() -> None:
    with pytest.raises(ValueError):
        TrigVector.from_modes([((0, 0, 0), (1, 0, 0), (0, 0, 0))])


# --------------------------------------------------------------------------- #
# the slab certificate                                                         #
# --------------------------------------------------------------------------- #


def _certificate():
    return build_torus_slab_certificate(
        "P1", viscosity=Fraction(1, 10), cutoff_sq=4, step=Fraction(1, 2048),
        bits=64,
    )


def test_slab_certificate_closes_and_verifies() -> None:
    """The end-to-end statement: box proved, control slab closed, checker green.

    A failure here is the headline failure — the periodic lane's first rigorous
    certificate no longer closes."""
    payload = _certificate().as_dict()
    assert payload["galerkin_box"]["proved"] is True
    assert payload["conclusion"]["proved"] is True
    verdict = verify_torus_slab_certificate(payload)
    assert verdict["verified"], verdict["failures"]


def test_slab_certificate_distance_is_meaningfully_small() -> None:
    """The bound must be far below the trajectory's own norm to mean anything:
    a tube wider than the solution would be a true statement about nothing."""
    payload = _certificate().as_dict()
    distance = Fraction(payload["conclusion"]["distance_h4_dot"])
    datum_norm_sq = Fraction(payload["datum"]["h4_norm_sq"])
    assert distance * distance < datum_norm_sq / 100


def test_slab_certificate_is_json_safe_and_carries_the_distinction() -> None:
    payload = _certificate().as_dict()
    json.dumps(payload, allow_nan=False)
    distinction = payload["distinction"]
    assert "REMAINS" in distinction["fixed_band_trajectory"]
    assert "NOT excluded" in distinction["finite_band_datum"]


def test_slab_checker_rejects_a_forged_external_theorem() -> None:
    payload = _certificate().as_dict()
    payload["external_theorems"]["EXT-P1"]["proved"] = True
    assert not verify_torus_slab_certificate(payload)["verified"]


def test_slab_checker_rejects_a_lean_axiomatisation_claim() -> None:
    """The session rule is that the external classical theorems are never
    inserted as Lean axioms; a payload recording otherwise must fail."""
    payload = _certificate().as_dict()
    payload["external_theorems"]["EXT-P2"]["axiomatised_in_lean"] = True
    assert not verify_torus_slab_certificate(payload)["verified"]


def test_slab_checker_rejects_a_tampered_coefficient() -> None:
    payload = _certificate().as_dict()
    payload["constants"]["K1"] = "1"
    assert not verify_torus_slab_certificate(payload)["verified"]


def test_slab_checker_rejects_a_dropped_disclaimer() -> None:
    payload = _certificate().as_dict()
    payload["claims"] = ["everything is fine"]
    assert not verify_torus_slab_certificate(payload)["verified"]


def test_slab_checker_returns_a_verdict_on_garbage() -> None:
    for garbage in ({}, {"certificate": "torus_slab_aposteriori"}, {"datum": 3}):
        verdict = verify_torus_slab_certificate(garbage)
        assert verdict["verified"] is False
        assert verdict["failures"]


def test_certificate_rejects_a_vanishing_viscosity() -> None:
    with pytest.raises(ValueError):
        build_torus_slab_certificate("P1", viscosity=Fraction(0))


def test_certificate_rejects_an_unknown_family() -> None:
    with pytest.raises(ValueError):
        build_torus_slab_certificate("Q9")


def test_galerkin_box_reports_failure_rather_than_inflating() -> None:
    """An absurdly long slab must come back ``proved: False`` (or with the step
    halved), never silently widened into meaninglessness."""
    box = prove_galerkin_box(
        family_P1(), viscosity=Fraction(1, 10), cutoff_sq=4,
        step=Fraction(50), max_doublings=3, max_halvings=1, bits=48,
    )
    assert (not box.proved) or box.step < Fraction(50)


def test_sqrt8_constant_is_an_upper_root() -> None:
    assert SQRT8_UPPER * SQRT8_UPPER >= 8
    assert EXTERNAL_THEOREMS.keys() == {"EXT-P1", "EXT-P2", "EXT-P3"}
