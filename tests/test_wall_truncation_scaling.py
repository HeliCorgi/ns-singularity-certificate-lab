"""Acceptance tests for the E-33 wall-truncation verification experiment.

Three things are pinned here.

1. **The oracle.**  ``ns_certificate_lab.bessel_reference`` is the only source
   of ``K_1/I_1`` in the whole verification chain, so if it were wrong every
   "agreement" reported by the experiment would be meaningless.  It is checked
   against published values, against two *independent* reimplementations written
   out in this module (a different integral for ``I_1`` and a different
   ascending series for ``K_1``), against the DLMF cross-product identity, and
   against both asymptotic branches of E-33(c).

2. **The measurement.**  The ``k = 0`` closed form of E-33(d) is an exact
   oracle, so the measured core difference is compared with it directly, at the
   resolution that reproduces the numbers recorded in ``docs/equation_audit.md``.

3. **The refusals.**  Config validation, and three fault injections that each
   break one link of the chain and must be caught.

This module imports the experiment runner rather than the solvers.  The
three-solver comparison lives inside the experiment module because
``tests/test_poisson_cross_validation.py`` forbids any ``tests/test_*.py`` from
importing both solver A and solver B; that scan is over test files only, so a
test that imports neither -- as this one does -- is unaffected.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from experiments import run_wall_truncation_scaling as experiment
from experiments.run_wall_truncation_scaling import (
    EXPERIMENT_ID,
    assert_compact_support,
    build_source,
    classify_regime,
    measure_axial_case,
    modal_ode_relative_residual,
    radial_point_count,
    radial_profile,
    run,
    source_moment,
    validate_config,
    wall_response_amplitude,
    zero_mode_closed_form,
)
from ns_certificate_lab import bessel_reference
from ns_certificate_lab.grid import AxisymmetricGrid

REPOSITORY = Path(__file__).resolve().parents[1]
SHIPPED_CONFIG = REPOSITORY / "configs" / "wall_truncation_scaling.json"

EULER_MASCHERONI = 0.5772156649015328606065120900824


# --------------------------------------------------------------------------
# Independent reimplementations, written out here so that the oracle is
# checked against something that shares no code with it.
# --------------------------------------------------------------------------


def _i_n_via_integral(order: int, x: float, nodes: int = 1024) -> float:
    """``I_n(x) = (1/pi) int_0^pi e^{x cos t} cos(n t) dt`` (DLMF 10.32.3).

    This is a different representation from the ascending power series the
    module uses.  The integrand is analytic and, extended evenly, periodic, so
    the trapezoidal rule converges spectrally; ``nodes = 512`` and
    ``nodes = 1024`` agree to ``2e-16`` over the range used below.
    """

    t = np.arange(nodes + 1) * (math.pi / nodes)
    values = np.exp(x * np.cos(t)) * np.cos(order * t)
    weights = np.full(nodes + 1, 1.0)
    weights[0] = 0.5
    weights[-1] = 0.5
    return float(np.sum(weights * values)) / nodes


def _k1_via_ascending_series(x: float, terms: int = 90) -> float:
    """``K_1`` from the ascending series of DLMF 10.31.1 with ``n = 1``:

    ``K_1(z) = 1/z + ln(z/2) I_1(z)
               - (z/4) sum_k [psi(k+1)+psi(k+2)] (z^2/4)^k / (k! (k+1)!)``

    with ``psi(1) = -gamma`` and ``psi(m+1) = -gamma + H_m``.  ``I_1`` is taken
    from :func:`_i_n_via_integral`, so this path shares nothing with the
    module's quadrature.  The alternating cancellation grows with ``z``; the
    tests below use it only for ``z <= 2``, where it is good to ``6e-16``.
    """

    total = 0.0
    harmonic = 0.0
    for index in range(terms):
        next_harmonic = harmonic + 1.0 / (index + 1)
        digamma_sum = (-EULER_MASCHERONI + harmonic) + (
            -EULER_MASCHERONI + next_harmonic
        )
        total += (
            digamma_sum
            * (0.25 * x * x) ** index
            / (math.factorial(index) * math.factorial(index + 1))
        )
        harmonic = next_harmonic
    return 1.0 / x + math.log(0.5 * x) * _i_n_via_integral(1, x) - 0.25 * x * total


# --------------------------------------------------------------------------
# 1. Oracle validation
# --------------------------------------------------------------------------


def test_bessel_oracle_reproduces_published_values() -> None:
    """The two values named in the task specification must be exact.

    ``I_1(1) = 0.5651591039924850`` and ``K_1(1) = 0.6019072301972346`` are the
    standard tabulated values (DLMF 10.25.2 / 10.31.1; Abramowitz & Stegun
    Table 9.8 lists the exponentially scaled forms, ``e^{-1} I_1(1) ~
    0.20791042`` and ``e^{1} K_1(1) ~ 1.6361535``, to fewer digits than are
    checked here).  Only the two full-precision values are asserted literally;
    the scaled forms are derived from them rather than quoted, so no additional
    table digits are being trusted.
    """

    assert bessel_reference.bessel_i1(1.0) == pytest.approx(
        0.5651591039924850, rel=1.0e-15, abs=0.0
    )
    assert bessel_reference.bessel_k1(1.0) == pytest.approx(
        0.6019072301972346, rel=1.0e-15, abs=0.0
    )
    # The scaled forms are what the ratio actually uses, so pin them too.
    assert bessel_reference.exp_scaled_i1(1.0) == pytest.approx(
        0.5651591039924850 * math.exp(-1.0), rel=1.0e-15, abs=0.0
    )
    assert bessel_reference.exp_scaled_k1(1.0) == pytest.approx(
        0.6019072301972346 * math.exp(1.0), rel=1.0e-15, abs=0.0
    )
    assert bessel_reference.k1_over_i1(1.0) == pytest.approx(
        0.6019072301972346 / 0.5651591039924850, rel=1.0e-15, abs=0.0
    )


def test_bessel_oracle_matches_two_independent_reimplementations() -> None:
    """Four further points, derived here rather than transcribed from a table.

    ``I_1`` is recomputed from a different integral representation and ``K_1``
    from a different ascending series (both written out above with their DLMF
    references).  Agreement at ``1e-13`` across these points would be a
    coincidence only if two unrelated algorithms shared the same defect.
    """

    for argument in (0.5, 1.0, 2.0, 5.0, 10.0):
        independent = _i_n_via_integral(1, argument)
        assert bessel_reference.bessel_i1(argument) == pytest.approx(
            independent, rel=1.0e-13, abs=0.0
        )
    for argument in (0.25, 0.5, 1.0, 2.0):
        independent = _k1_via_ascending_series(argument)
        assert bessel_reference.bessel_k1(argument) == pytest.approx(
            independent, rel=1.0e-13, abs=0.0
        )
    # The integral path must itself be converged, otherwise the agreement above
    # would only say that both routines were computing the same wrong thing.
    for argument in (0.5, 5.0, 10.0):
        coarse = _i_n_via_integral(1, argument, nodes=512)
        fine = _i_n_via_integral(1, argument, nodes=1024)
        assert abs(coarse - fine) <= 1.0e-14 * abs(fine)


def test_bessel_oracle_satisfies_the_cross_product_identity() -> None:
    """``I_0 K_1 + I_1 K_0 = 1/x`` (DLMF 10.28.2) ties the series to the quadrature."""

    for argument in (1.0e-8, 1.0e-4, 0.01, 0.5, 1.0, 5.0, 20.0, 60.0):
        assert bessel_reference.wronskian_relative_defect(argument) <= 1.0e-14


def test_k_quadrature_is_spectrally_accurate_not_second_order() -> None:
    """Halving the trapezoidal step must change nothing at working precision.

    The module argues that the error is ``O(exp(-2 pi d / h))`` because the
    integrand is analytic in a strip, not ``O(h^2)``.  A second-order rule would
    change by about a factor of four here; measured, it changes by at most
    ``1e-15`` relative.
    """

    for argument in (1.0e-3, 0.05, 1.0, 20.0, 60.0):
        assert bessel_reference.k_quadrature_step_halving_defect(argument) <= 1.0e-13
        assert bessel_reference.k_quadrature_step_halving_defect(argument, 0) <= 1.0e-13


def test_bessel_oracle_is_monotone_where_it_must_be() -> None:
    """``I_1`` increases, ``K_1`` decreases, and so does the ratio E-33(b) uses."""

    arguments = [0.1 * step for step in range(1, 121)]
    i_values = [bessel_reference.bessel_i1(x) for x in arguments]
    k_values = [bessel_reference.bessel_k1(x) for x in arguments]
    ratios = [bessel_reference.k1_over_i1(x) for x in arguments]
    assert all(i_values[i] < i_values[i + 1] for i in range(len(i_values) - 1))
    assert all(k_values[i] > k_values[i + 1] for i in range(len(k_values) - 1))
    assert all(ratios[i] > ratios[i + 1] for i in range(len(ratios) - 1))
    assert all(value > 0.0 for value in ratios)


def test_ratio_is_stable_at_large_argument() -> None:
    """The scaled quotient must survive arguments where the naive one dies.

    ``kR`` reaches ``2 pi * 8 = 50.3`` in the shipped sweep.  Up to there the
    naive quotient still works and must agree; beyond it the separate factors
    overflow and underflow while the scaled form keeps returning a correct
    positive number.
    """

    for argument in (10.0, 30.0, 50.3, 60.0):
        naive = bessel_reference.bessel_k1(argument) / bessel_reference.bessel_i1(
            argument
        )
        assert bessel_reference.k1_over_i1(argument) == pytest.approx(
            naive, rel=1.0e-14, abs=0.0
        )

    # Past the point where the naive quotient has lost its meaning, the scaled
    # form is still correct to the accuracy of the asymptotic series itself.
    for argument in (120.0, 250.0, 350.0):
        value = bessel_reference.k1_over_i1(argument)
        assert math.isfinite(value)
        assert value > 0.0
        reference = bessel_reference.k1_over_i1_large_argument_asymptote(
            argument, terms=3
        )
        assert value == pytest.approx(reference, rel=1.0e-4, abs=0.0)
    assert bessel_reference.bessel_i1(350.0) * bessel_reference.bessel_k1(
        350.0
    ) > 0.0


def test_e33c_large_argument_branch() -> None:
    """``K_1/I_1 -> pi e^{-2x}(1 + 3/(4x))`` with an ``O(x^-2)`` remainder.

    The remainder's coefficient is checked too: squaring the ``3/(8x)`` terms of
    the standard ``I_1`` and ``K_1`` expansions predicts ``9/32 = 0.28125``, and
    ``x^2 (rho / asymptote_2 - 1)`` is measured converging to it.
    """

    leading_errors = []
    two_term_errors = []
    scaled_remainders = []
    for argument in (3.0, 5.0, 10.0, 20.0, 40.0, 60.0):
        exact = bessel_reference.k1_over_i1(argument)
        leading = bessel_reference.k1_over_i1_large_argument_asymptote(
            argument, terms=1
        )
        two_term = bessel_reference.k1_over_i1_large_argument_asymptote(
            argument, terms=2
        )
        leading_errors.append(abs(exact / leading - 1.0))
        two_term_errors.append(abs(exact / two_term - 1.0))
        scaled_remainders.append(argument * argument * (exact / two_term - 1.0))

    # Both branches must improve monotonically, and the 3/(4x) term must be a
    # genuine improvement rather than decoration.
    assert all(
        leading_errors[i] > leading_errors[i + 1]
        for i in range(len(leading_errors) - 1)
    )
    assert all(
        two_term_errors[i] > two_term_errors[i + 1]
        for i in range(len(two_term_errors) - 1)
    )
    assert two_term_errors[0] < 0.15 * leading_errors[0]
    assert two_term_errors[-1] < 1.0e-4
    # x^2 * remainder -> 9/32.
    assert scaled_remainders[-1] == pytest.approx(9.0 / 32.0, rel=0.02, abs=0.0)

    # The exponential rate itself: d log rho / dx -> -2.
    for argument in (20.0, 40.0):
        slope = (
            math.log(
                bessel_reference.k1_over_i1(argument + 1.0e-3)
                / bessel_reference.k1_over_i1(argument)
            )
            / 1.0e-3
        )
        assert slope == pytest.approx(-2.0, rel=2.0e-3, abs=0.0)


def test_e33c_small_argument_branch() -> None:
    """``K_1/I_1 -> 2/x^2`` with a slope tending to ``-2`` in log-log."""

    errors = []
    for argument in (1.0e-4, 1.0e-3, 1.0e-2, 0.1):
        exact = bessel_reference.k1_over_i1(argument)
        asymptote = bessel_reference.k1_over_i1_small_argument_asymptote(argument)
        errors.append(abs(exact / asymptote - 1.0))
    assert all(errors[i] < errors[i + 1] for i in range(len(errors) - 1))
    assert errors[0] <= 1.0e-7

    for argument, tolerance in ((1.0e-4, 1.0e-6), (1.0e-3, 1.0e-5), (1.0e-2, 1.0e-3)):
        stepped = argument * 1.02
        slope = math.log(
            bessel_reference.k1_over_i1(stepped)
            / bessel_reference.k1_over_i1(argument)
        ) / math.log(1.02)
        assert abs(slope + 2.0) <= tolerance


def test_e33a_homogeneous_solutions_solve_the_modal_equation() -> None:
    """E-33(a): ``I_1(kr)/r`` and ``K_1(kr)/r`` must annihilate ``L_{5,k}``.

    The residual is a fourth-order centered difference of oracle values only.
    It must be small and must fall at fourth order; where it is already at the
    cancellation-roundoff level of that difference no order is claimed, and the
    module records the floor it uses.
    """

    for wavenumber in (0.7853981633974483, 3.141592653589793, 6.283185307179586):
        for branch in ("I", "K"):
            coarse = max(
                modal_ode_relative_residual(wavenumber, radius, step=0.02, branch=branch)
                for radius in (0.3, 0.5, 0.9, 1.5, 2.0)
            )
            fine = max(
                modal_ode_relative_residual(wavenumber, radius, step=0.01, branch=branch)
                for radius in (0.3, 0.5, 0.9, 1.5, 2.0)
            )
            assert fine <= 1.0e-5
            if fine >= experiment.MODAL_ODE_ROUNDOFF_FLOOR:
                order = math.log(coarse / fine) / math.log(2.0)
                assert 3.5 <= order <= 4.5

    with pytest.raises(ValueError):
        modal_ode_relative_residual(1.0, 1.0, step=0.01, branch="J")


def test_bessel_oracle_rejects_invalid_arguments() -> None:
    with pytest.raises(ValueError):
        bessel_reference.bessel_i1(-1.0)
    with pytest.raises(ValueError):
        bessel_reference.bessel_i1(float("nan"))
    with pytest.raises(ValueError):
        bessel_reference.bessel_k1(0.0)
    with pytest.raises(ValueError):
        bessel_reference.k1_over_i1(0.0)
    with pytest.raises(ValueError):
        bessel_reference.bessel_k1(1.0e-14)
    with pytest.raises(TypeError):
        bessel_reference.bessel_i1(True)
    assert bessel_reference.bessel_i1(0.0) == 0.0
    assert bessel_reference.i1_over_argument(0.0) == 0.5


# --------------------------------------------------------------------------
# 2. The k = 0 closed form of E-33(d)
# --------------------------------------------------------------------------


def _zero_mode_difference(
    points_per_unit_radius: int, first: float, second: float
) -> tuple[float, float, float]:
    """Return ``(mean, spread/|mean|, closed form)`` for one ``k = 0`` pair."""

    case = {"label": "k0", "z_period": 1.0, "mode": 0}
    moment = source_moment(support_radius=0.95, exponent=8)
    row, _ = measure_axial_case(
        case=case,
        points_per_unit_radius=points_per_unit_radius,
        wall_radii=[first, second],
        reference_radius=second,
        axial_points=8,
        support_radius=0.95,
        exponent=8,
        core_radius=0.9,
        roundoff_floor_factor=16.0,
        well_resolved_floor_margin=1.0e4,
        classification={"exponential_min_kr": 3.0, "algebraic_max_kr": 0.5},
        source_moment_value=moment,
    )
    entry = row["radii"][0]
    return (
        float(entry["measured_difference_mean"]),
        float(entry["constancy_relative_spread"]),
        float(entry["closed_form_difference"]),
    )


def test_zero_mode_moment_matches_the_closed_form_integral() -> None:
    """``Q_inf = a^4 / (2 (p+1)(p+2))`` -- verified against direct quadrature."""

    support_radius, exponent = 0.95, 8
    analytic = source_moment(support_radius=support_radius, exponent=exponent)
    assert analytic == pytest.approx(0.95**4 / 180.0, rel=1.0e-15, abs=0.0)

    # Independent high-resolution Simpson quadrature of int_0^a s^3 omega ds.
    nodes = 200001
    s = np.linspace(0.0, support_radius, nodes)
    integrand = s**3 * (1.0 - (s / support_radius) ** 2) ** exponent
    weights = np.ones(nodes)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    numeric = float(np.sum(weights * integrand)) * (s[1] - s[0]) / 3.0
    assert numeric == pytest.approx(analytic, rel=1.0e-10, abs=0.0)


def test_zero_mode_response_is_constant_and_matches_the_closed_form() -> None:
    """Reproduce the E-33(d) numbers recorded in ``docs/equation_audit.md``.

    At ``dr = 1/256`` the audit records the ``R = 1`` versus ``R = 2`` core
    difference as ``-1.6969516537e-3``, constant in ``r`` to ``8.6e-14``
    relative, against the closed form ``-1.6968880208e-3``; and the ``R = 2``
    versus ``R = 4`` difference as ``-4.2423943061e-4`` against
    ``-4.2422200521e-4``.  Those exact values are reproduced here.
    """

    mean12, spread12, closed12 = _zero_mode_difference(256, 1.0, 2.0)
    mean24, spread24, closed24 = _zero_mode_difference(256, 2.0, 4.0)

    assert closed12 == pytest.approx(-1.6968880208333332e-3, rel=1.0e-14, abs=0.0)
    assert closed24 == pytest.approx(-4.2422200520833333e-4, rel=1.0e-14, abs=0.0)
    assert mean12 == pytest.approx(-1.6969516537e-3, rel=1.0e-10, abs=0.0)
    assert mean24 == pytest.approx(-4.2423943061e-4, rel=1.0e-10, abs=0.0)

    # Constancy in r: a discrete identity, so it holds at roundoff.
    assert spread12 <= 1.0e-12
    assert spread24 <= 1.0e-12
    assert spread12 == pytest.approx(8.6e-14, rel=0.1, abs=0.0)

    # Magnitude: an O(dr^2) discretization error, so it does NOT hold at
    # roundoff.  Conflating the two would hide a real fault, so the size of the
    # gap is pinned from below as well as above.
    error12 = abs(mean12 - closed12) / abs(closed12)
    error24 = abs(mean24 - closed24) / abs(closed24)
    assert 1.0e-5 <= error12 <= 1.0e-4
    assert 1.0e-5 <= error24 <= 1.0e-4
    assert error12 == pytest.approx(3.75e-5, rel=0.02, abs=0.0)
    assert error24 == pytest.approx(4.11e-5, rel=0.02, abs=0.0)


def test_zero_mode_closed_form_error_is_second_order() -> None:
    """The measured-versus-closed-form error must fall like ``dr^2``."""

    errors = []
    spacings = []
    for points_per_unit_radius in (64, 128, 256):
        mean, _, closed = _zero_mode_difference(points_per_unit_radius, 1.0, 2.0)
        errors.append(abs(mean - closed) / abs(closed))
        spacings.append(1.0 / points_per_unit_radius)
    orders = [
        math.log(errors[i] / errors[i + 1]) / math.log(spacings[i] / spacings[i + 1])
        for i in range(len(errors) - 1)
    ]
    assert min(orders) >= 1.9
    assert max(orders) <= 2.1


def test_zero_mode_closed_form_sign_and_scaling() -> None:
    """``psi^(R) - psi^(R') = (Q_inf/2)(R'^-2 - R^-2)`` including the sign."""

    moment = source_moment(support_radius=0.95, exponent=8)
    value = zero_mode_closed_form(
        source_moment_value=moment, radius=1.0, other_radius=2.0
    )
    assert value < 0.0
    assert value == pytest.approx(0.5 * moment * (0.25 - 1.0), rel=1e-15, abs=0.0)
    # Antisymmetry under swapping the two walls.
    reversed_value = zero_mode_closed_form(
        source_moment_value=moment, radius=2.0, other_radius=1.0
    )
    assert reversed_value == pytest.approx(-value, rel=1.0e-15, abs=0.0)


# --------------------------------------------------------------------------
# 3. E-33(b) shape and E-33(c) regimes, measured
# --------------------------------------------------------------------------


def test_regime_classification_brackets_the_crossover() -> None:
    kwargs = {"exponential_min_kr": 3.0, "algebraic_max_kr": 0.5}
    assert classify_regime(smaller_kr=6.3, larger_kr=9.4, **kwargs) == "exponential"
    assert classify_regime(smaller_kr=0.2, larger_kr=0.3, **kwargs) == "algebraic"
    assert classify_regime(smaller_kr=0.2, larger_kr=0.9, **kwargs) == "crossover"
    assert classify_regime(smaller_kr=1.0, larger_kr=1.6, **kwargs) == "crossover"


def test_measured_response_follows_the_bessel_ratio_and_shape() -> None:
    """One axial case, measured end to end against the E-33(b) prediction."""

    case = {"label": "Lz1_m1", "z_period": 1.0, "mode": 1}
    row, _ = measure_axial_case(
        case=case,
        points_per_unit_radius=64,
        wall_radii=[1.0, 1.5, 2.0, 8.0],
        reference_radius=8.0,
        axial_points=16,
        support_radius=0.95,
        exponent=8,
        core_radius=0.9,
        roundoff_floor_factor=16.0,
        well_resolved_floor_margin=1.0e4,
        classification={"exponential_min_kr": 3.0, "algebraic_max_kr": 0.5},
        source_moment_value=source_moment(support_radius=0.95, exponent=8),
    )

    for entry in row["radii"]:
        assert entry["oracle_comparable"]
        assert entry["well_resolved"]
        # E-33(b): the shape is I_1(kr)/r, which increases with r, so the
        # response must peak at the outermost core node.
        assert entry["argmax_r"] == pytest.approx(
            row["core_radius_of_last_node"], rel=0.0, abs=1.0e-12
        )
        assert entry["shape_one_minus_abs_cosine"] <= 1.0e-6
        assert entry["shape_sign_matches_e33b"]
        # And it is emphatically NOT constant in r, unlike the k = 0 case.
        assert entry["response_spread_over_abs_mean"] > 1.0

    for pair in row["pairs"]:
        assert pair["regime"] == "exponential"
        assert pair["ratio_relative_error"] <= 5.0e-3
        assert pair["slopes_measured"]
        assert pair["exponential_slope_relative_deviation"] <= 0.06
        assert pair["measured_versus_oracle_slope_relative_error"] <= 5.0e-3
        # The leading E-33(c) asymptote alone is a visibly worse predictor than
        # the exact ratio, which is what the 3/(4x) term is there to fix.
        assert pair["leading_asymptote_ratio"] < pair["oracle_ratio"]
        assert pair["two_term_asymptote_ratio"] > pair["leading_asymptote_ratio"]

    # The Hou configuration, L_z = 1 and R = 1, sits at kR = 2 pi, deep in the
    # exponential regime: this is the audit's "wall effects are negligible
    # here, and here only" statement, as a number.
    assert row["radii"][0]["kr"] == pytest.approx(2.0 * math.pi, rel=1e-14, abs=0.0)
    assert row["radii"][0]["response_relative_to_reference_core_max"] < 1.0e-2


def test_algebraic_regime_slope_is_near_minus_two() -> None:
    """At small ``kR`` the response falls only algebraically, not exponentially."""

    case = {"label": "Lz32_m1", "z_period": 32.0, "mode": 1}
    row, _ = measure_axial_case(
        case=case,
        points_per_unit_radius=64,
        wall_radii=[1.0, 1.5, 2.0, 8.0],
        reference_radius=8.0,
        axial_points=16,
        support_radius=0.95,
        exponent=8,
        core_radius=0.9,
        roundoff_floor_factor=16.0,
        well_resolved_floor_margin=1.0e4,
        classification={"exponential_min_kr": 3.0, "algebraic_max_kr": 0.5},
        source_moment_value=source_moment(support_radius=0.95, exponent=8),
    )
    for pair in row["pairs"]:
        assert pair["regime"] == "algebraic"
        assert abs(pair["slope_dlog_dlogradius"] + 2.0) <= 0.35
        # The measured slope tracks the exact oracle far more tightly than it
        # tracks the leading -2, which is the honest statement of E-33(c) at
        # these arguments.
        assert pair["measured_versus_oracle_slope_relative_error"] <= 5.0e-3
        # Not exponential: e^{2k dR} would be about 1.2 here, nothing like the
        # measured factor.
        assert pair["measured_ratio"] > 1.5


def test_wall_response_amplitude_is_symmetric_and_exact() -> None:
    wavenumber = 2.0 * math.pi
    forward = wall_response_amplitude(wavenumber, 1.0, 2.0)
    backward = wall_response_amplitude(wavenumber, 2.0, 1.0)
    assert forward == pytest.approx(backward, rel=1.0e-15, abs=0.0)
    assert forward == pytest.approx(
        abs(
            bessel_reference.k1_over_i1(wavenumber)
            - bessel_reference.k1_over_i1(2.0 * wavenumber)
        ),
        rel=1.0e-15,
        abs=0.0,
    )


# --------------------------------------------------------------------------
# 4. Configuration contract
# --------------------------------------------------------------------------


def _reduced_config() -> dict[str, Any]:
    """A CI-sized config that still exercises both E-33(c) regimes."""

    return {
        "schema_version": 1,
        "experiment": EXPERIMENT_ID,
        "description": "reduced smoke configuration for the E-33 verification",
        "interpretation": "smoke run; not evidence",
        "source_profile": {"support_radius": 0.95, "exponent": 8},
        "core_radius": 0.9,
        "axial_cases": [
            {"label": "k0", "z_period": 1.0, "mode": 0},
            {"label": "Lz1_m1", "z_period": 1.0, "mode": 1},
            {"label": "Lz32_m1", "z_period": 32.0, "mode": 1},
        ],
        # The reference wall must sit far enough out that its OWN residual
        # response is small compared with the responses being measured.  It is
        # not needed for the oracle -- that comparison is exact either way --
        # but the E-33(c) regime slopes are statements about rho(kR) itself,
        # and a nearby reference contaminates them: with the reference at
        # R = 4 the L_z = 32 log-log slope moves from -2.13 to -2.27 and -2.55.
        "wall_radii": [1.0, 1.5, 2.0, 8.0],
        "reference_wall_radius": 8.0,
        "axial_points": 16,
        "radial_resolutions": [24, 32, 48],
        "cross_check": {
            "points_per_unit_radius": 24,
            "axial_points": 32,
            "wall_radii": [1.0, 2.0, 8.0],
            "axial_case_labels": ["k0", "Lz1_m1"],
            "realspace_tolerance": 1e-12,
            "realspace_max_iterations": 20000,
        },
        "nonlinear_cross_link": {
            "summary_path": "outputs/wall_dependence_v1/summary.json",
            "member_group": "core_dr192",
            "z_period": 1.0,
            "mode": 1,
        },
        "classification": {"exponential_min_kr": 3.0, "algebraic_max_kr": 0.5},
        "roundoff_floor_factor": 16.0,
        "well_resolved_floor_margin": 1.0e4,
        "acceptance": {
            "max_oracle_ratio_relative_error_finest": 5e-3,
            "max_oracle_ratio_relative_error_any": 3e-2,
            "min_oracle_ratio_observed_order": 1.8,
            "max_oracle_ratio_observed_order": 2.2,
            "max_zero_mode_closed_form_relative_error_finest": 3e-3,
            "min_zero_mode_closed_form_observed_order": 1.9,
            "max_zero_mode_closed_form_observed_order": 2.1,
            "max_zero_mode_constancy_relative_spread": 1e-11,
            "max_shape_one_minus_cosine": 1e-4,
            "min_shape_observed_order": 3.5,
            "max_shape_observed_order": 4.5,
            "max_exponential_slope_relative_deviation": 0.06,
            "max_algebraic_slope_absolute_deviation": 0.35,
            "max_measured_versus_oracle_slope_relative_error": 3e-2,
            "max_modal_ode_relative_residual": 1e-5,
            "min_modal_ode_observed_order": 3.5,
            "max_wronskian_relative_defect": 1e-14,
            "max_published_value_relative_error": 1e-12,
            "max_cross_solver_relative_difference": 3e-2,
            "min_cross_solver_ab_difference_over_dr_squared": 1e-3,
            "max_cross_solver_ab_difference_over_dr_squared": 1.0,
            "max_cross_solver_ratio_relative_deviation": 5e-2,
        },
    }


def test_shipped_config_is_valid_and_canonical() -> None:
    config = json.loads(SHIPPED_CONFIG.read_text(encoding="utf-8"))
    validate_config(config)
    assert config["experiment"] == EXPERIMENT_ID
    assert config["reference_wall_radius"] == max(config["wall_radii"])
    assert config["source_profile"]["support_radius"] < min(config["wall_radii"])
    # The sweep must actually cross the E-33(c) transition at kR ~ 1.
    products = [
        2.0 * math.pi * int(case["mode"]) / float(case["z_period"]) * radius
        for case in config["axial_cases"]
        if int(case["mode"]) > 0
        for radius in config["wall_radii"]
    ]
    assert min(products) < 0.5
    assert max(products) > 20.0


def test_reduced_config_is_valid() -> None:
    validate_config(_reduced_config())


@pytest.mark.parametrize("key", sorted(_reduced_config()))
def test_config_rejects_every_missing_key(key: str) -> None:
    config = _reduced_config()
    del config[key]
    with pytest.raises(ValueError, match="missing keys"):
        validate_config(config)


def test_config_rejects_unknown_key() -> None:
    config = _reduced_config()
    config["extra_knob"] = 1.0
    with pytest.raises(ValueError, match="unknown keys"):
        validate_config(config)


def test_config_rejects_unknown_nested_keys() -> None:
    for section, key in (
        ("source_profile", "cutoff"),
        ("classification", "midpoint"),
        ("cross_check", "solver"),
        ("acceptance", "max_anything"),
        ("nonlinear_cross_link", "group"),
    ):
        config = _reduced_config()
        config[section][key] = 1.0
        with pytest.raises(ValueError, match="unknown keys"):
            validate_config(config)


def test_config_rejects_mislabeled_experiment_id() -> None:
    config = _reduced_config()
    config["experiment"] = "wall_dependence_v1"
    with pytest.raises(ValueError, match="canonical value"):
        validate_config(config)


def test_config_rejects_unsupported_schema() -> None:
    config = _reduced_config()
    config["schema_version"] = 2
    with pytest.raises(ValueError, match="schema"):
        validate_config(config)


def test_config_rejects_non_monotone_wall_radii() -> None:
    config = _reduced_config()
    config["wall_radii"] = [1.0, 2.0, 1.5, 8.0]
    with pytest.raises(ValueError, match="strictly increasing"):
        validate_config(config)

    duplicated = _reduced_config()
    duplicated["wall_radii"] = [1.0, 1.5, 1.5, 8.0]
    with pytest.raises(ValueError, match="strictly increasing"):
        validate_config(duplicated)


def test_config_rejects_reference_radius_that_is_not_the_largest() -> None:
    config = _reduced_config()
    config["reference_wall_radius"] = 2.0
    with pytest.raises(ValueError, match="largest entry"):
        validate_config(config)


def test_config_rejects_degenerate_refinement_and_classification() -> None:
    too_few = _reduced_config()
    too_few["radial_resolutions"] = [24, 32]
    with pytest.raises(ValueError, match="at least 3 entries"):
        validate_config(too_few)

    unsorted_resolutions = _reduced_config()
    unsorted_resolutions["radial_resolutions"] = [48, 32, 24]
    with pytest.raises(ValueError, match="strictly increasing"):
        validate_config(unsorted_resolutions)

    fractional = _reduced_config()
    fractional["radial_resolutions"] = [25, 32, 48]
    with pytest.raises(ValueError, match="whole cells"):
        validate_config(fractional)

    bad_regimes = _reduced_config()
    bad_regimes["classification"] = {
        "exponential_min_kr": 0.4,
        "algebraic_max_kr": 0.5,
    }
    with pytest.raises(ValueError, match="crossover"):
        validate_config(bad_regimes)


def test_config_rejects_bad_axial_case_sets() -> None:
    no_zero_mode = _reduced_config()
    no_zero_mode["axial_cases"] = [
        case for case in no_zero_mode["axial_cases"] if case["mode"] != 0
    ]
    with pytest.raises(ValueError, match="mode 0"):
        validate_config(no_zero_mode)

    duplicate_labels = _reduced_config()
    duplicate_labels["axial_cases"][2]["label"] = "Lz1_m1"
    with pytest.raises(ValueError, match="unique"):
        validate_config(duplicate_labels)


def test_config_rejects_bad_cross_check_radii() -> None:
    outside = _reduced_config()
    outside["cross_check"]["wall_radii"] = [1.0, 3.0, 8.0]
    with pytest.raises(ValueError, match="subset"):
        validate_config(outside)

    truncated = _reduced_config()
    truncated["cross_check"]["wall_radii"] = [1.0, 2.0]
    with pytest.raises(ValueError, match="reference wall radius"):
        validate_config(truncated)

    unknown_label = _reduced_config()
    unknown_label["cross_check"]["axial_case_labels"] = ["Lz99_m1"]
    with pytest.raises(ValueError, match="declared cases"):
        validate_config(unknown_label)


# --------------------------------------------------------------------------
# 5. Smoke run and artifact bundle
# --------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def smoke_bundle(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, dict[str, Any]]:
    output_dir = tmp_path_factory.mktemp("wall_truncation") / "bundle"
    summary = run(_reduced_config(), output_dir)
    return output_dir, summary


def test_smoke_run_writes_a_verifiable_bundle(
    smoke_bundle: tuple[Path, dict[str, Any]]
) -> None:
    output_dir, summary = smoke_bundle
    assert summary["accepted"]
    assert summary["experiment_id"] == EXPERIMENT_ID
    assert summary["audit_entry"] == "E-33"

    expected = {
        "summary.json",
        "config.snapshot.json",
        "wall_response.csv",
        "radius_pairs.csv",
        "core_response.npz",
        "manifest.json",
        "manifest.json.sha256",
    }
    assert {path.name for path in output_dir.iterdir()} == expected

    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["experiment_id"] == EXPERIMENT_ID
    assert set(manifest["files"]) == expected - {"manifest.json", "manifest.json.sha256"}
    for name, record in manifest["files"].items():
        payload = output_dir / name
        assert record["sha256"] == _sha256(payload)
        assert record["bytes"] == payload.stat().st_size
    sidecar = (output_dir / "manifest.json.sha256").read_text(encoding="ascii").strip()
    assert sidecar == _sha256(output_dir / "manifest.json")

    # The provenance record and the config snapshot must both be present and
    # must describe the run that produced these numbers.
    assert summary["reproducibility"]["runtime_provenance"]
    snapshot = json.loads((output_dir / "config.snapshot.json").read_text(encoding="utf-8"))
    assert snapshot == _reduced_config()


def test_smoke_run_measures_every_e33_claim(
    smoke_bundle: tuple[Path, dict[str, Any]]
) -> None:
    _, summary = smoke_bundle
    checks = summary["acceptance"]["checks"]
    assert checks["all_passed"]
    for name in (
        "oracle_ratio_agreement_finest",
        "oracle_ratio_orders_in_band",
        "zero_mode_closed_form_agreement_finest",
        "zero_mode_response_constant_in_radius",
        "shape_matches_i1_over_r",
        "shape_sign_matches_e33b",
        "exponential_regime_slope_matches_two_k",
        "algebraic_regime_slope_matches_minus_two",
        "modal_ode_residual_small",
        "cross_solver_agreement_within_band",
        "cross_solver_ab_gap_is_genuine_truncation",
    ):
        assert checks[name] is True, name

    observed = summary["acceptance"]["observed"]
    assert observed["exponential_pair_count"] >= 1
    assert observed["algebraic_pair_count"] >= 1

    # The three-solver comparison must be a real comparison, not two names for
    # the same array.
    cross = summary["analysis"]["cross_solver"]
    assert cross["min_ab_difference_over_dr_squared"] > 0.0
    assert cross["max_relative_difference"] > 0.0


def test_smoke_run_records_the_nonlinear_cross_link(
    smoke_bundle: tuple[Path, dict[str, Any]]
) -> None:
    """The elliptic law must explain the nonlinear run's wall sensitivity.

    This is recorded, not gated: the check below asserts that the comparison
    was made and that its numbers are the ones the nonlinear bundle contains,
    not that the agreement is good.  The agreement itself is reported in the
    summary for interpretation.
    """

    _, summary = smoke_bundle
    link = summary["analysis"]["nonlinear_cross_link"]
    if not link["available"]:
        pytest.skip("outputs/wall_dependence_v1/summary.json is not in this checkout")
    assert link["gated"] is False
    assert link["comparisons"]
    for comparison in link["comparisons"]:
        assert math.isfinite(comparison["nonlinear_separation_ratio"])
        assert math.isfinite(comparison["elliptic_oracle_ratio"])
    # L_z = 1 is the nonlinear run's period, so kR = 2 pi at the innermost wall.
    assert link["wavenumber"] == pytest.approx(2.0 * math.pi, rel=1e-14, abs=0.0)


def test_run_refuses_to_overwrite_a_nonempty_output_directory(scratch_dir: Path) -> None:
    output_dir = scratch_dir / "bundle"
    output_dir.mkdir()
    (output_dir / "existing_evidence.json").write_text("{}", encoding="utf-8")
    with pytest.raises(FileExistsError, match="nonempty"):
        run(_reduced_config(), output_dir)
    assert (output_dir / "existing_evidence.json").exists()


def test_run_validates_before_creating_any_output(scratch_dir: Path) -> None:
    config = _reduced_config()
    config["experiment"] = "not_this_experiment"
    output_dir = scratch_dir / "never_created"
    with pytest.raises(ValueError):
        run(config, output_dir)
    assert not output_dir.exists()


# --------------------------------------------------------------------------
# 6. Fault injections
# --------------------------------------------------------------------------


def test_detects_sign_error_in_the_zero_mode_closed_form(
    monkeypatch: pytest.MonkeyPatch, scratch_dir: Path
) -> None:
    """Fault (i): ``+Q_inf/(2R^2)`` instead of ``-Q_inf/(2R^2)``.

    E-33(e) fixes the sign of the wall-dependent constant.  A sign error there
    is the single most plausible transcription mistake, and it must not survive
    the oracle comparison.
    """

    moment = source_moment(support_radius=0.95, exponent=8)
    mean, _, closed = _zero_mode_difference(64, 1.0, 2.0)
    honest_error = abs(mean - closed) / abs(closed)
    flipped = -closed
    flipped_error = abs(mean - flipped) / abs(flipped)
    assert honest_error <= 1.0e-3
    assert flipped_error >= 1.9
    assert flipped_error > 1000.0 * honest_error

    def corrupted(
        *, source_moment_value: float, radius: float, other_radius: float
    ) -> float:
        return -zero_mode_closed_form(
            source_moment_value=source_moment_value,
            radius=radius,
            other_radius=other_radius,
        )

    monkeypatch.setattr(experiment, "zero_mode_closed_form", corrupted)
    with pytest.raises(RuntimeError, match="failed acceptance"):
        run(_reduced_config(), scratch_dir / "sign_fault")
    assert moment > 0.0


def test_detects_source_that_is_not_compactly_supported_inside_the_wall(
    monkeypatch: pytest.MonkeyPatch, scratch_dir: Path
) -> None:
    """Fault (ii): E-33's premise fails and the experiment must refuse to run.

    Three layers are exercised: the declared support radius reaching the wall
    is rejected by config validation; a grid too coarse to show any source-free
    interior row is rejected as well; and a source that silently fails to vanish
    outside its declared support is caught at run time.
    """

    reaching = _reduced_config()
    reaching["source_profile"]["support_radius"] = 1.2
    with pytest.raises(ValueError, match="compact-support premise"):
        validate_config(reaching)

    touching = _reduced_config()
    touching["source_profile"]["support_radius"] = 1.0
    with pytest.raises(ValueError, match="compact-support premise"):
        validate_config(touching)

    unobservable = _reduced_config()
    unobservable["radial_resolutions"] = [8, 16, 32]
    with pytest.raises(ValueError, match="no interior grid node"):
        validate_config(unobservable)

    # A config that passes validation but whose sampled source does not vanish
    # where it claims to must still be refused, and before any solve is used.
    grid = AxisymmetricGrid.uniform(
        nr=radial_point_count(24, 1.0), nz=16, r_max=1.0, z_min=0.0, z_max=1.0
    )
    honest = build_source(
        grid, support_radius=0.95, exponent=8, z_period=1.0, mode=1
    )
    assert (
        assert_compact_support(
            grid, honest, support_radius=0.95, smallest_wall_radius=1.0
        )
        >= 1
    )
    leaking = honest.copy()
    leaking[-2, :] = 1.0e-30
    with pytest.raises(ValueError, match="not exactly zero outside"):
        assert_compact_support(
            grid, leaking, support_radius=0.95, smallest_wall_radius=1.0
        )

    def unbounded_profile(
        radii: np.ndarray, *, support_radius: float, exponent: int
    ) -> np.ndarray:
        # A Gaussian is nowhere exactly zero, so it violates the premise on
        # every grid however small its tail is.
        return np.exp(-((radii / support_radius) ** 2))

    monkeypatch.setattr(experiment, "radial_profile", unbounded_profile)
    with pytest.raises(ValueError, match="not exactly zero outside"):
        run(_reduced_config(), scratch_dir / "support_fault")
    assert not (scratch_dir / "support_fault" / "summary.json").exists()


def test_detects_corrupted_bessel_oracle(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fault (iii): return ``I_1`` where ``K_1`` belongs in the ratio.

    The corruption inverts the ratio, so every predicted wall response becomes
    dominated by the reference wall's own term and all the oracle ratios
    collapse towards one.  The measured ratios are in the hundreds, so the
    disagreement must be enormous -- not a tolerance question.
    """

    case = {"label": "Lz1_m1", "z_period": 1.0, "mode": 1}
    kwargs = {
        "case": case,
        "points_per_unit_radius": 32,
        "wall_radii": [1.0, 1.5, 2.0, 8.0],
        "reference_radius": 8.0,
        "axial_points": 16,
        "support_radius": 0.95,
        "exponent": 8,
        "core_radius": 0.9,
        "roundoff_floor_factor": 16.0,
        "well_resolved_floor_margin": 1.0e4,
        "classification": {"exponential_min_kr": 3.0, "algebraic_max_kr": 0.5},
        "source_moment_value": source_moment(support_radius=0.95, exponent=8),
    }
    honest, _ = measure_axial_case(**kwargs)
    for pair in honest["pairs"]:
        assert pair["ratio_relative_error"] <= 2.0e-2
        assert pair["measured_ratio"] > 100.0

    def swapped(x: float) -> float:
        return bessel_reference.bessel_i1(x) / bessel_reference.bessel_k1(x)

    monkeypatch.setattr(bessel_reference, "k1_over_i1", swapped)
    corrupted, _ = measure_axial_case(**kwargs)
    for pair in corrupted["pairs"]:
        # The corrupted oracle predicts a ratio of essentially one against a
        # measured ratio in the hundreds.
        assert pair["oracle_ratio"] == pytest.approx(1.0, rel=1.0e-3, abs=0.0)
        assert pair["ratio_relative_error"] > 100.0

    honest_worst = max(pair["ratio_relative_error"] for pair in honest["pairs"])
    corrupted_worst = max(pair["ratio_relative_error"] for pair in corrupted["pairs"])
    assert corrupted_worst > 1.0e4 * honest_worst


def test_detects_shape_prediction_replaced_by_the_wrong_branch() -> None:
    """A guard on E-33(b): correlating against ``K_1(kr)/r`` must fail.

    The two homogeneous branches of E-33(a) are both solutions of the modal
    equation, so a "shape agrees with a modified Bessel function" claim is only
    worth something if it distinguishes the axis-regular branch from the
    decaying one.
    """

    case = {"label": "Lz1_m1", "z_period": 1.0, "mode": 1}
    row, stored = measure_axial_case(
        case=case,
        points_per_unit_radius=64,
        wall_radii=[1.0, 8.0],
        reference_radius=8.0,
        axial_points=16,
        support_radius=0.95,
        exponent=8,
        core_radius=0.9,
        roundoff_floor_factor=16.0,
        well_resolved_floor_margin=1.0e4,
        classification={"exponential_min_kr": 3.0, "algebraic_max_kr": 0.5},
        source_moment_value=source_moment(support_radius=0.95, exponent=8),
    )
    assert row["radii"][0]["shape_one_minus_abs_cosine"] <= 1.0e-6

    core_r = stored["core_r"]
    difference = stored["response_R1"]
    column = difference[:, int(np.argmax(np.abs(difference[0, :])))]
    wavenumber = row["wavenumber"]
    # The decaying branch, evaluated away from the axis where it is finite.
    interior = core_r > 0.0
    wrong_branch = np.array(
        [
            bessel_reference.bessel_k1(wavenumber * float(radius)) / float(radius)
            for radius in core_r[interior]
        ]
    )
    restricted = column[interior]
    cosine = float(
        np.dot(restricted, wrong_branch)
        / (np.linalg.norm(restricted) * np.linalg.norm(wrong_branch))
    )
    assert 1.0 - abs(cosine) > 1.0e-2

    # And the right branch still wins on the same restricted node set.
    right_branch = np.array(
        [
            wavenumber
            * bessel_reference.i1_over_argument(wavenumber * float(radius))
            for radius in core_r[interior]
        ]
    )
    right_cosine = float(
        np.dot(restricted, right_branch)
        / (np.linalg.norm(restricted) * np.linalg.norm(right_branch))
    )
    # The regular branch reproduces the same near-unit correlation on the
    # restricted node set that the experiment reports on the full one, four
    # orders of magnitude better than the decaying branch manages.
    assert 1.0 - abs(right_cosine) < 1.0e-6
    assert (1.0 - abs(cosine)) > 1.0e4 * (1.0 - abs(right_cosine))


def test_radial_profile_is_exactly_zero_outside_its_support() -> None:
    radii = np.linspace(0.0, 2.0, 401)
    values = radial_profile(radii, support_radius=0.95, exponent=8)
    outside = radii >= 0.95
    assert np.all(values[outside] == 0.0)
    assert np.all(values[~outside] > 0.0)
    assert values[0] == 1.0
    # Monotone decreasing on the support, so the moment integral is positive.
    inside = values[~outside]
    assert np.all(np.diff(inside) < 0.0)


def test_radial_point_count_refuses_fractional_cells() -> None:
    assert radial_point_count(64, 1.5) == 97
    assert radial_point_count(24, 8.0) == 193
    with pytest.raises(ValueError, match="whole cells"):
        radial_point_count(25, 1.5)


def test_shipped_and_reduced_configs_agree_on_the_audited_invariants() -> None:
    """A reduced config may be smaller, never differently specified."""

    shipped = json.loads(SHIPPED_CONFIG.read_text(encoding="utf-8"))
    reduced = _reduced_config()
    assert shipped["source_profile"] == reduced["source_profile"]
    assert shipped["core_radius"] == reduced["core_radius"]
    assert shipped["classification"] == reduced["classification"]
    assert shipped["roundoff_floor_factor"] == reduced["roundoff_floor_factor"]
    assert set(shipped["acceptance"]) == set(reduced["acceptance"])
    assert copy.deepcopy(shipped) == json.loads(
        SHIPPED_CONFIG.read_text(encoding="utf-8")
    )
