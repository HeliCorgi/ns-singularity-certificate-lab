"""Tests for the diverging-bandwidth exponent conditions."""

from __future__ import annotations

import math

import pytest

from ns_certificate_lab.shell_constraints import (
    MAXIMUM_BANDWIDTH_EXPONENT,
    ShellExponents,
    classify_shell_exponents,
    feasible_sigma_interval,
    feasible_sigma_lower_is_attained,
    scan_shell_exponents,
    smooth_force_shell_bound,
)


def _point(gamma: float, sigma: float, beta: float = 1.0) -> ShellExponents:
    return ShellExponents(bandwidth=gamma, peak_energy=sigma, spectral_slope=beta)


def test_feasible_interval_is_empty_outside_zero_one():
    assert feasible_sigma_interval(0.0) is None
    assert feasible_sigma_interval(-0.3) is None
    assert feasible_sigma_interval(MAXIMUM_BANDWIDTH_EXPONENT) is None
    assert feasible_sigma_interval(1.5) is None


@pytest.mark.parametrize("gamma", [0.1, 0.25, 0.5, 0.75, 0.9])
def test_feasible_interval_endpoints_match_the_derivation(gamma):
    interval = feasible_sigma_interval(gamma)
    assert interval is not None
    lower, upper = interval
    assert lower == pytest.approx(max(0.0, 2.0 * gamma - 1.0))
    assert upper == pytest.approx(gamma)
    assert lower < upper


@pytest.mark.parametrize("gamma", [0.1, 0.3, 0.5, 0.7, 0.95])
def test_interior_points_are_admissible_and_endpoints_are_not(gamma):
    lower, upper = feasible_sigma_interval(gamma)
    middle = 0.5 * (lower + upper)
    assert classify_shell_exponents(_point(gamma, middle)).feasible
    # Both endpoints fail, and they fail for the stated reasons.
    at_upper = classify_shell_exponents(_point(gamma, upper))
    assert "critical_norm_diverges" in at_upper.failed_conditions
    at_lower = classify_shell_exponents(_point(gamma, lower))
    # The lower endpoint is attained exactly below gamma = 1/2, where it is the
    # non-strict energy bound sigma >= 0 rather than the strict dissipation one.
    assert at_lower.feasible == feasible_sigma_lower_is_attained(gamma)
    if not at_lower.feasible:
        assert "dissipation_integrable" in at_lower.failed_conditions


def test_bandwidth_at_or_above_one_is_excluded_for_every_sigma():
    """``gamma >= 1`` is the sharp no-go: no `sigma` survives."""
    for gamma in (1.0, 1.25, 2.0):
        for sigma in (-1.0, 0.0, 0.5, 1.0, 1.5, 3.0):
            verdict = classify_shell_exponents(_point(gamma, sigma))
            assert not verdict.feasible


def test_bounded_bandwidth_is_excluded_by_the_finite_mode_no_go():
    verdict = classify_shell_exponents(_point(0.0, 0.0))
    assert "bandwidth_diverges" in verdict.failed_conditions


def test_negative_energy_exponent_is_excluded():
    verdict = classify_shell_exponents(_point(0.5, -0.1))
    assert "energy_bounded" in verdict.failed_conditions


def test_flat_spectrum_is_outside_the_historical_positive_slope_branch():
    """The corrected finite-floor classifier treats beta=0 piecewise.

    This assertion is retained only for compatibility with the preregistered
    beta-positive scan; see ``test_zeno_packet_relay.py`` for the correction.
    """
    verdict = classify_shell_exponents(_point(0.5, 0.25, beta=0.0))
    assert "spectrum_summable" in verdict.failed_conditions


def test_classical_type_one_self_similar_point_is_excluded_by_dissipation():
    """``gamma = 1/2``, ``sigma = 0`` is the parabolic self-similar scaling.

    Its enstrophy integral is logarithmically divergent, so the point sits
    exactly on the excluded boundary.  This is the exponent-level shadow of the
    known non-existence of self-similar blow-up.
    """
    verdict = classify_shell_exponents(_point(0.5, 0.0))
    assert "dissipation_integrable" in verdict.failed_conditions
    assert verdict.exponents.dissipation_time_exponent == pytest.approx(-1.0)


def test_derived_time_exponents_are_the_stated_formulas():
    point = _point(0.4, 0.3)
    assert point.dissipation_time_exponent == pytest.approx(0.3 - 0.8)
    assert point.critical_norm_time_exponent == pytest.approx(0.5 * (0.3 - 0.4))
    assert point.sobolev_exponent_threshold == pytest.approx(0.3 / 0.8)
    assert math.isinf(_point(0.0, 0.0).sobolev_exponent_threshold)


def test_scan_reports_a_nonempty_but_strictly_bounded_feasible_set():
    gammas = [0.05 * i for i in range(1, 30)]
    sigmas = [0.05 * i for i in range(0, 30)]
    verdicts = scan_shell_exponents(gammas, sigmas, [1.0])
    feasible = [v for v in verdicts if v.feasible]
    assert feasible, "the feasible triangle must not be empty"
    assert len(feasible) < len(verdicts), "the conditions must exclude something"
    for verdict in feasible:
        assert 0.0 < verdict.exponents.bandwidth < MAXIMUM_BANDWIDTH_EXPONENT
        assert verdict.exponents.peak_energy < verdict.exponents.bandwidth


def test_exponents_reject_non_finite_input():
    with pytest.raises(ValueError):
        ShellExponents(bandwidth=math.inf, peak_energy=0.0, spectral_slope=1.0)
    with pytest.raises(ValueError):
        ShellExponents(bandwidth=0.5, peak_energy=math.nan, spectral_slope=1.0)


def test_smooth_force_bound_decays_faster_than_any_power():
    """The force is invisible in the high-shell budget: that is the point."""
    for order in (1, 4, 16):
        assert smooth_force_shell_bound(1024.0, order, 1.0) < smooth_force_shell_bound(
            1024.0, order - 1, 1.0
        )
    assert smooth_force_shell_bound(2.0**20, 8, 1.0) < 1.0e-40
    with pytest.raises(ValueError):
        smooth_force_shell_bound(0.0, 1, 1.0)
    with pytest.raises(ValueError):
        smooth_force_shell_bound(2.0, -1, 1.0)
    with pytest.raises(ValueError):
        smooth_force_shell_bound(2.0, 1, -1.0)


def test_verdict_dict_is_json_friendly():
    payload = classify_shell_exponents(_point(0.4, 0.3)).as_dict()
    assert payload["verdict"] == "conditionally_admissible"
    assert payload["failed_conditions"] == []
    payload = classify_shell_exponents(_point(1.5, 0.3)).as_dict()
    assert payload["verdict"] == "excluded"
    assert payload["failed_conditions"]
