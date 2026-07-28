from __future__ import annotations

import numpy as np
import pytest

from ns_certificate_lab.scaling_constraints import CandidateFamily
from ns_certificate_lab.scaling_fit import (
    fit_joint_scaling,
    fit_log_power_law,
    scan_joint_scaling_fits,
    summarize_fit_scan,
)


def test_log_power_law_recovers_exact_exponent() -> None:
    singular_time = 1.0
    times = np.array([0.1, 0.3, 0.5, 0.7, 0.8])
    values = 2.5 * (singular_time - times) ** 0.45
    fit = fit_log_power_law(times, values, singular_time=singular_time)
    assert fit.exponent == pytest.approx(0.45, abs=1.0e-12)
    assert fit.prefactor == pytest.approx(2.5, abs=1.0e-12)
    assert fit.log_rmse < 1.0e-12
    assert fit.r_squared == pytest.approx(1.0)


def test_joint_fit_recovers_euler_like_type_ii_family() -> None:
    singular_time = 1.0
    times = np.array([0.2, 0.4, 0.6, 0.75, 0.85])
    tau = singular_time - times
    amplitudes = 3.0 * tau ** (-0.55)
    radial = 0.7 * tau**0.45
    axial = 0.9 * tau**0.45
    fit = fit_joint_scaling(
        times,
        amplitudes,
        radial,
        axial,
        singular_time=singular_time,
        balance_tolerance=1.0e-8,
    )
    assert fit.alpha == pytest.approx(0.55, abs=1.0e-12)
    assert fit.beta_r == pytest.approx(0.45, abs=1.0e-12)
    assert fit.beta_z == pytest.approx(0.45, abs=1.0e-12)
    assert fit.candidate_family is CandidateFamily.EULER_LIKE_TYPE_II
    assert fit.passes_conditional_core_gates


def test_expanding_widths_cannot_pass_shrinking_core_gate() -> None:
    singular_time = 1.0
    times = np.array([0.2, 0.4, 0.6, 0.75, 0.85])
    tau = singular_time - times
    amplitudes = tau ** (-0.55)
    radial = tau ** (-0.10)
    axial = tau ** (-0.05)
    fit = fit_joint_scaling(
        times,
        amplitudes,
        radial,
        axial,
        singular_time=singular_time,
    )
    assert fit.beta_r < 0.0
    assert fit.beta_z < 0.0
    assert fit.constraint_report is None
    assert not fit.passes_conditional_core_gates


def test_scan_finds_true_candidate_time_for_exact_synthetic_data() -> None:
    true_time = 1.0
    times = np.array([0.2, 0.4, 0.6, 0.75, 0.85])
    tau = true_time - times
    amplitudes = tau ** (-0.55)
    radial = tau**0.45
    axial = tau**0.45
    fits = scan_joint_scaling_fits(
        times,
        amplitudes,
        radial,
        axial,
        singular_times=[0.95, 1.0, 1.1],
        minimum_points=4,
        balance_tolerance=1.0e-8,
    )
    summary = summarize_fit_scan(fits)
    best = summary["best_unconstrained"]
    assert isinstance(best, dict)
    assert best["singular_time"] == pytest.approx(true_time)
    assert summary["passing_fit_count"] > 0


def test_fit_validation() -> None:
    with pytest.raises(ValueError):
        fit_log_power_law([0.0, 0.5], [1.0, 2.0], singular_time=0.5)
    with pytest.raises(ValueError):
        fit_log_power_law([0.0, 0.5], [1.0, -2.0], singular_time=1.0)
    with pytest.raises(ValueError):
        scan_joint_scaling_fits(
            [0.0, 0.5],
            [1.0, 2.0],
            [1.0, 0.5],
            [1.0, 0.5],
            singular_times=[1.0],
            minimum_points=3,
        )
