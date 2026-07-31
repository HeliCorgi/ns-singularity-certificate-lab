"""Tests for the ``n = 3`` Kato-constant certificate (Lane K)."""

from __future__ import annotations

from fractions import Fraction
import json
import math

import pytest

from ns_certificate_lab.fourier_torus import (
    TrigScalar,
    TrigVector,
    _IntervalOps,
    adot_squared_upper,
    family_P1,
)
from ns_certificate_lab.kato_constant import (
    CONVENTIONS,
    a4_upper,
    a6_upper,
    build_kato_certificate,
    c_kato_upper,
    c_shift_upper,
    g3_upper,
    verify_kato_certificate,
)
from ns_certificate_lab.snapshot_certificate import Interval


# --------------------------------------------------------------------------- #
# lattice sums                                                                 #
# --------------------------------------------------------------------------- #


def test_tiny_lattice_partial_sums_by_hand() -> None:
    """The ``N = 1`` box has exactly 26 nonzero points, checkable by hand.

    The 6 unit vectors (``|p|^2 = 1``) contribute exactly ``6 * 1^{-4} = 6``
    to ``A_4`` and ``6 * 1^{-6} = 6`` to ``A_6``.  The 12 face diagonals
    (``|p|^2 = 2``) and the 8 corners (``|p|^2 = 3``) complete the box:

    ``A_4`` partial: ``6/1 + 12/4 + 8/9 = 89/9``;
    ``A_6`` partial: ``6/1 + 12/8 + 8/27 = 421/54``.
    """
    unit_vectors_a4 = 6 * Fraction(1, 1**2)
    unit_vectors_a6 = 6 * Fraction(1, 1**3)
    assert unit_vectors_a4 == 6
    assert unit_vectors_a6 == 6
    partial_a4 = unit_vectors_a4 + 12 * Fraction(1, 2**2) + 8 * Fraction(1, 3**2)
    partial_a6 = unit_vectors_a6 + 12 * Fraction(1, 2**3) + 8 * Fraction(1, 3**3)
    assert partial_a4 == Fraction(89, 9)
    assert partial_a6 == Fraction(421, 54)
    # partial + tail, with the tails at their closed forms 26/N and 9/N^3.
    assert a4_upper(1) == Fraction(89, 9) + Fraction(26, 1)
    assert a6_upper(1) == Fraction(421, 54) + Fraction(9, 1)


def test_a4_upper_equals_the_adot_mechanics_exactly() -> None:
    """``a4_upper`` re-exports the ``adot_squared_upper`` mechanics: the
    shell-grouped sum must equal the pointwise sum to the last digit."""
    for cut in (1, 4, 9):
        assert a4_upper(cut) == adot_squared_upper(cut)


def test_lattice_bounds_reject_a_degenerate_cut() -> None:
    for bad in (0, -3):
        with pytest.raises(ValueError):
            a4_upper(bad)
        with pytest.raises(ValueError):
            a6_upper(bad)


def test_bounds_are_monotone_improving_across_cuts() -> None:
    """Doubling the exact region can only tighten an upper bound built from
    partial sum + dominating tail; strict decrease across 10 -> 20 -> 40."""
    a4_values = [a4_upper(cut) for cut in (10, 20, 40)]
    a6_values = [a6_upper(cut) for cut in (10, 20, 40)]
    g3_values = [g3_upper(cut, bits=64) for cut in (10, 20, 40)]
    for coarse, fine in zip(a4_values, a4_values[1:]):
        assert fine < coarse
    for coarse, fine in zip(a6_values, a6_values[1:]):
        assert fine < coarse
    for coarse, fine in zip(g3_values, g3_values[1:]):
        assert fine < coarse


def test_g3_matches_the_documented_table() -> None:
    """The derivation note records G3 <= 49.9451 / 49.3672 / 49.1752 at
    N = 20 / 40 / 60; the certificate must reproduce them within rounding."""
    documented = {20: 49.9451, 40: 49.3672, 60: 49.1752}
    for cut, expected in documented.items():
        assert float(g3_upper(cut, bits=64)) == pytest.approx(expected, abs=1e-3)


# --------------------------------------------------------------------------- #
# per-trajectory band constants                                                #
# --------------------------------------------------------------------------- #


def test_c_kato_single_mode_by_hand() -> None:
    """One mode at ``k = (1,0,0)`` with ``a = (0,3,0)``, ``b = (0,0,4)``:
    ``sqrt(|a|^2 + |b|^2) = sqrt(9 + 16) = 5`` exactly, ``|p| = 1``, so
    ``C_kato = 6 * (1 + 1) * 5 = 60`` — and the dyadic square roots of
    perfect squares are exact, so equality is exact."""
    field = TrigVector.from_modes([((1, 0, 0), (0, 3, 0), (0, 0, 4))])
    assert c_kato_upper(field, bits=64) == 60


def test_c_shift_single_mode_by_hand() -> None:
    """Same field: ``C_shift = |j| (1 + |j|)^3 sqrt(|a|^2+|b|^2)
    = 1 * 2^3 * 5 = 40`` exactly."""
    field = TrigVector.from_modes([((1, 0, 0), (0, 3, 0), (0, 0, 4))])
    assert c_shift_upper(field, bits=64) == 40


def test_helpers_bound_an_irrational_mode_tightly_from_above() -> None:
    """At ``k = (1,1,0)`` with ``a = (1,-1,0)``: coefficient ``sqrt 2``,
    ``|p| = sqrt 2``.  Hand values ``C_kato = 6 * sqrt2 * 3 * sqrt2 = 36``
    and ``C_shift = sqrt2 (1+sqrt2)^3 sqrt2 = 14 + 10 sqrt2``; the computed
    rational bounds must sit above them and within outward rounding."""
    field = TrigVector.from_modes([((1, 1, 0), (1, -1, 0), (0, 0, 0))])
    kato = c_kato_upper(field, bits=64)
    shift = c_shift_upper(field, bits=64)
    assert kato >= 36
    assert kato - 36 < Fraction(1, 10**9)
    exact_shift = 14 + 10 * math.sqrt(2)
    assert float(shift) >= exact_shift - 1e-12
    assert float(shift) == pytest.approx(exact_shift, abs=1e-9)


def test_helpers_dominate_a_float_recomputation_on_family_P1() -> None:  # noqa: N802
    """The helpers on the preregistered P1 datum against a float re-run of
    the section 5a/5b formulas: an upper bound, and a tight one."""
    field = family_P1()
    kato_hand = 0.0
    shift_hand = 0.0
    for k, (a, b) in field.coefficient_table().items():
        n = sum(c * c for c in k)
        norm = math.sqrt(n)
        coefficient = math.sqrt(
            sum(float(x) ** 2 for x in a) + sum(float(x) ** 2 for x in b)
        )
        kato_hand += 6 * (norm + norm**3) * coefficient
        shift_hand += norm * (1 + norm) ** 3 * coefficient
    kato = c_kato_upper(field, bits=64)
    shift = c_shift_upper(field, bits=64)
    assert float(kato) >= kato_hand - 1e-9
    assert float(kato) == pytest.approx(kato_hand, abs=1e-6)
    assert float(shift) >= shift_hand - 1e-9
    assert float(shift) == pytest.approx(shift_hand, abs=1e-6)


def test_helpers_accept_interval_coefficients() -> None:
    """A Galerkin-box-style field with interval coefficients: the magnitudes
    of ``[-3, 3]`` and ``[0, 4]`` are 3 and 4, so the single-mode hand values
    60 and 40 must be reproduced exactly."""
    ops = _IntervalOps(64)
    zero = Interval(Fraction(0), Fraction(0))
    empty = TrigScalar({}, ops)
    cos_slot = TrigScalar(
        {(1, 0, 0): [Interval(Fraction(-3), Fraction(3)), zero]}, ops
    )
    sin_slot = TrigScalar(
        {(1, 0, 0): [zero, Interval(Fraction(0), Fraction(4))]}, ops
    )
    field = TrigVector((empty, cos_slot, sin_slot))
    assert c_kato_upper(field, bits=64) == 60
    assert c_shift_upper(field, bits=64) == 40


# --------------------------------------------------------------------------- #
# the certificate and its checker                                              #
# --------------------------------------------------------------------------- #


def _payload(cut: int = 10) -> dict:
    return build_kato_certificate(lattice_cut=cut, bits=64)


def test_checker_accepts_the_builder_output() -> None:
    payload = _payload()
    json.dumps(payload, allow_nan=False)
    verdict = verify_kato_certificate(payload)
    assert verdict["verified"], verdict["failures"]
    assert verdict["failures"] == []


def test_default_build_matches_the_note_and_verifies() -> None:
    """The default cut is 40; the note's headline numbers must be reproduced
    and the checker (including the doubled-cut recomputation at 80) green."""
    payload = build_kato_certificate()
    assert payload["lattice_cut"] == 40
    assert float(Fraction(payload["a4"]["upper"])) == pytest.approx(16.92444, abs=1e-4)
    assert float(Fraction(payload["g3"]["upper"])) == pytest.approx(49.3672, abs=1e-3)
    mono = payload["monotonicity"]
    assert Fraction(mono["g3_upper_at_double_cut"]) <= Fraction(mono["g3_upper_at_cut"])
    verdict = verify_kato_certificate(payload)
    assert verdict["verified"], verdict["failures"]


def test_builder_rejects_degenerate_inputs() -> None:
    with pytest.raises(ValueError):
        build_kato_certificate(lattice_cut=0)
    with pytest.raises(ValueError):
        build_kato_certificate(bits=4)


def test_checker_rejects_a_forged_g3() -> None:
    payload = _payload()
    payload["g3"]["upper"] = "1"
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_rejects_a_forged_a4_partial_sum() -> None:
    payload = _payload()
    payload["a4"]["partial_sum"] = str(Fraction(payload["a4"]["partial_sum"]) / 2)
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_rejects_a_forged_a4_tail() -> None:
    payload = _payload()
    payload["a4"]["tail"] = "1/1000000"
    payload["a4"]["upper"] = str(
        Fraction(payload["a4"]["partial_sum"]) + Fraction(1, 1000000)
    )
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_rejects_dropped_or_altered_conventions() -> None:
    dropped = _payload()
    del dropped["conventions"]
    assert not verify_kato_certificate(dropped)["verified"]

    altered = _payload()
    altered["conventions"] = list(CONVENTIONS)
    altered["conventions"][1] = "measure | unnormalised Lebesgue"
    assert not verify_kato_certificate(altered)["verified"]

    truncated = _payload()
    truncated["conventions"] = list(CONVENTIONS)[:-1]
    assert not verify_kato_certificate(truncated)["verified"]


def test_checker_rejects_the_comparison_promoted_to_input() -> None:
    flipped = _payload()
    flipped["literature_comparison"]["verification_only"] = False
    assert not verify_kato_certificate(flipped)["verified"]

    promoted = _payload()
    promoted["literature_comparison"]["role"] = (
        "input to the chain linear coefficient"
    )
    assert not verify_kato_certificate(promoted)["verified"]


def test_checker_rejects_a_forged_monotonicity_witness() -> None:
    payload = _payload()
    payload["monotonicity"]["a4_upper_at_double_cut"] = "1"
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_rejects_a_forged_sqrt() -> None:
    payload = _payload()
    payload["g3"]["sqrt_a4_upper"] = "1"
    payload["g3"]["upper"] = "12"
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_rejects_a_tampered_k3_coefficient() -> None:
    payload = _payload()
    payload["k3"]["coeff_v3_w3"] = str(Fraction(payload["k3"]["coeff_v3_w3"]) / 2)
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_rejects_a_dropped_disclaimer() -> None:
    payload = _payload()
    payload["claims"] = ["everything is fine"]
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_rejects_dropped_sharpening_levers() -> None:
    payload = _payload()
    payload["sharpening_levers"] = payload["sharpening_levers"][:2]
    assert not verify_kato_certificate(payload)["verified"]


def test_checker_returns_a_verdict_on_garbage() -> None:
    for garbage in ({}, {"schema": "kato_h3_constant"}, {"a4": 3}, [], "junk"):
        verdict = verify_kato_certificate(garbage)
        assert verdict["verified"] is False
        assert verdict["failures"]
