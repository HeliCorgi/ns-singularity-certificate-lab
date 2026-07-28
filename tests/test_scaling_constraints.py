from __future__ import annotations

import math

import pytest

from ns_certificate_lab.scaling_constraints import (
    BalanceFamily,
    CandidateFamily,
    PowerLawScaling,
    analyze_power_law_scaling,
    generic_family_bounds,
    isotropic_type_ii_wedge_contains,
    scan_power_law_scalings,
)


def test_lp_and_serrin_exponents_match_change_of_variables() -> None:
    scaling = PowerLawScaling(alpha=0.55, beta_r=0.45, beta_z=0.45)
    assert scaling.volume_exponent == pytest.approx(1.35)
    assert scaling.lp_growth_exponent(3.0) == pytest.approx(0.10)
    assert scaling.lp_growth_exponent(6.0) == pytest.approx(0.325)
    assert scaling.serrin_time_exponent(6.0) == pytest.approx(1.30)


def test_standard_parabolic_isotropic_scaling_has_bounded_l3() -> None:
    report = analyze_power_law_scaling(
        PowerLawScaling(alpha=0.5, beta_r=0.5, beta_z=0.5)
    )
    assert report.balance_family is BalanceFamily.TRIPLE_PARABOLIC
    assert report.energy_margin == pytest.approx(0.5)
    assert report.critical_l3_margin == pytest.approx(0.0)
    assert not report.critical_l3_divergent
    assert not report.passes_conditional_core_gates


def test_isotropic_type_ii_wedge_point_passes() -> None:
    report = analyze_power_law_scaling(
        PowerLawScaling(alpha=0.55, beta_r=0.45, beta_z=0.45)
    )
    assert report.balance_family is BalanceFamily.TIME_INERTIAL
    assert report.candidate_family is CandidateFamily.EULER_LIKE_TYPE_II
    assert report.energy_margin == pytest.approx(0.25)
    assert report.critical_l3_margin == pytest.approx(0.30)
    assert report.dissipation_margin == pytest.approx(0.35)
    assert report.passes_conditional_core_gates
    assert isotropic_type_ii_wedge_contains(0.55, 0.45)


def test_isotropic_type_ii_wedge_boundaries() -> None:
    assert isotropic_type_ii_wedge_contains(0.6, 0.4)
    assert not isotropic_type_ii_wedge_contains(0.5, 0.5)
    assert not isotropic_type_ii_wedge_contains(0.61, 0.39)
    assert not isotropic_type_ii_wedge_contains(0.54, 0.45)


def test_anisotropic_parabolic_family_passes() -> None:
    # gamma=1/2, B=1.2: all three modeled PDE terms balance, while the
    # anisotropic core L3 mass grows and total dissipation remains integrable.
    report = analyze_power_law_scaling(
        PowerLawScaling(alpha=0.5, beta_r=0.5, beta_z=0.2)
    )
    assert report.volume_exponent == pytest.approx(1.2)
    assert report.balance_family is BalanceFamily.TRIPLE_PARABOLIC
    assert report.candidate_family is CandidateFamily.ANISOTROPIC_PARABOLIC
    assert report.passes_conditional_core_gates


def test_quasi_steady_viscous_inertial_family_passes() -> None:
    # gamma=alpha=0.6 and B=1.5 lies in 4 gamma - 1 < B < 3 gamma.
    report = analyze_power_law_scaling(
        PowerLawScaling(alpha=0.6, beta_r=0.6, beta_z=0.3)
    )
    assert report.volume_exponent == pytest.approx(1.5)
    assert report.balance_family is BalanceFamily.INERTIAL_VISCOUS
    assert report.candidate_family is CandidateFamily.QUASI_STEADY_VISCOUS_INERTIAL
    assert report.passes_conditional_core_gates


def test_single_dominant_term_requires_cancellation() -> None:
    report = analyze_power_law_scaling(
        PowerLawScaling(alpha=0.7, beta_r=0.5, beta_z=0.5)
    )
    assert report.balance_family is BalanceFamily.REQUIRES_LEADING_CANCELLATION
    assert report.candidate_family is CandidateFamily.NONE
    assert not report.passes_conditional_core_gates


def test_finite_energy_and_dissipation_reject_bad_candidates() -> None:
    energy_bad = analyze_power_law_scaling(
        PowerLawScaling(alpha=0.8, beta_r=0.4, beta_z=0.4)
    )
    assert not energy_bad.finite_energy_compatible

    dissipation_bad = analyze_power_law_scaling(
        PowerLawScaling(alpha=0.6, beta_r=0.6, beta_z=1.0e-2)
    )
    assert not dissipation_bad.finite_dissipation_compatible


def test_generic_family_bounds_are_exact_at_representative_points() -> None:
    first = generic_family_bounds(0.45)
    assert first["family"] == CandidateFamily.EULER_LIKE_TYPE_II.value
    assert first["alpha"] == pytest.approx(0.55)
    assert first["volume_lower"] == pytest.approx(1.10)
    assert first["volume_upper"] == pytest.approx(1.35)

    middle = generic_family_bounds(0.5)
    assert middle["family"] == CandidateFamily.ANISOTROPIC_PARABOLIC.value
    assert middle["volume_lower"] == pytest.approx(1.0)
    assert middle["volume_upper"] == pytest.approx(1.5)
    assert middle["lower_open"] and middle["upper_open"]

    third = generic_family_bounds(0.6)
    assert third["family"] == CandidateFamily.QUASI_STEADY_VISCOUS_INERTIAL.value
    assert third["alpha"] == pytest.approx(0.6)
    assert third["volume_lower"] == pytest.approx(1.4)
    assert third["volume_upper"] == pytest.approx(1.8)

    assert not generic_family_bounds(0.3)["feasible"]
    assert not generic_family_bounds(1.0)["feasible"]


def test_scan_returns_all_three_generic_families() -> None:
    reports = scan_power_law_scalings(
        alphas=[0.5, 0.55, 0.6],
        beta_rs=[0.2, 0.3, 0.4, 0.45, 0.5, 0.6],
        beta_zs=[0.2, 0.3, 0.4, 0.45, 0.5, 0.6],
    )
    families = {report.candidate_family for report in reports}
    assert CandidateFamily.EULER_LIKE_TYPE_II in families
    assert CandidateFamily.ANISOTROPIC_PARABOLIC in families
    assert CandidateFamily.QUASI_STEADY_VISCOUS_INERTIAL in families


def test_input_validation() -> None:
    with pytest.raises(ValueError):
        PowerLawScaling(alpha=math.nan, beta_r=0.5, beta_z=0.5)
    with pytest.raises(ValueError):
        PowerLawScaling(alpha=0.5, beta_r=0.0, beta_z=0.5)
    with pytest.raises(ValueError):
        PowerLawScaling(alpha=-0.1, beta_r=0.5, beta_z=0.5)
    with pytest.raises(ValueError):
        PowerLawScaling(alpha=0.5, beta_r=0.5, beta_z=0.5).serrin_time_exponent(3.0)
