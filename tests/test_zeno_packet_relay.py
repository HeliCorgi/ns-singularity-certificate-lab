"""Checks for the finite-floor shell correction and pressure-relay pilot."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ns_certificate_lab.zeno_packet_relay import (
    MovingFrontLaw,
    classify_finite_floor_shell_exponents,
    critical_staircase_metrics,
    finite_floor_shell_sums,
    gaussian_vortex_moment,
    normalized_pressure_stretch,
    pressure_hessian_quadrupole,
    zeno_scaling_exponents,
)


def test_critical_staircase_has_bounded_energy_and_growing_critical_mass():
    for level in (4, 8, 16):
        report = critical_staircase_metrics(level)
        assert report.energy == pytest.approx(2.0 - 2.0 ** (-level))
        assert report.enstrophy == pytest.approx(2.0 ** (level + 1) - 1.0)
        assert report.critical_besov == pytest.approx(level + 1.0)
        assert report.critical_l3_cubed == pytest.approx(level + 1.0)
        assert report.vorticity_linf_scale == pytest.approx(report.bandwidth**2)
        assert report.nonlinear_l2_scale == pytest.approx(report.bandwidth**1.5)


def test_finite_floor_beta_minus_one_is_not_excluded_by_energy():
    """The torus has levels 0..J, not an infinite tail below level zero."""

    verdict = classify_finite_floor_shell_exponents(0.5, 0.5, -1.0)
    assert verdict.energy_bounded
    assert verdict.dissipation_integrable
    assert verdict.critical_upper_bound_not_uniform
    assert verdict.logarithmic_critical_boundary
    assert verdict.kinematically_admissible

    for level in (8, 16, 24):
        front = 2.0**level
        peak_energy = 1.0 / front
        energy, enstrophy, critical_besov = finite_floor_shell_sums(
            level,
            spectral_slope=-1.0,
            peak_energy=peak_energy,
        )
        assert energy == pytest.approx(2.0 - 2.0 ** (-level))
        assert enstrophy / front == pytest.approx(2.0 - 2.0 ** (-level))
        assert critical_besov == pytest.approx(level + 1.0)


def test_flat_finite_floor_has_the_expected_logarithm():
    level = 20
    energy, _, _ = finite_floor_shell_sums(
        level,
        spectral_slope=0.0,
        peak_energy=1.0,
    )
    assert energy == pytest.approx(level + 1.0)
    assert not classify_finite_floor_shell_exponents(0.5, 0.0, 0.0).energy_bounded
    assert classify_finite_floor_shell_exponents(0.5, 0.1, 0.0).energy_bounded


def test_far_field_pressure_hessian_has_favourable_transverse_stretch():
    # This is the normalized energy moment of the explicit divergence-free
    # Schwartz field W=(0,-z,y) exp(-|x|^2/2), not an unrealizable rank-one
    # algebraic tensor.
    moment = gaussian_vortex_moment(np.array([1.0, 0.0, 0.0]))
    displacement = np.array([1.0, 0.0, 0.0])
    hessian = pressure_hessian_quadrupole(moment, displacement)
    expected = np.diag([-12.0, 6.0, 6.0]) / (4.0 * math.pi)
    assert np.allclose(hessian, expected, rtol=0.0, atol=1.0e-14)
    assert np.trace(hessian) == pytest.approx(0.0, abs=1.0e-14)
    assert normalized_pressure_stretch(
        moment, displacement, np.array([1.0, 0.0, 0.0])
    ) == pytest.approx(12.0)
    assert normalized_pressure_stretch(
        moment, displacement, np.array([0.0, 1.0, 0.0])
    ) == pytest.approx(-6.0)


def test_gaussian_vortex_moment_is_physically_realizable_rank_two_tensor():
    moment = gaussian_vortex_moment(np.array([0.0, 0.0, 4.0]))
    assert np.allclose(moment, np.diag([0.5, 0.5, 0.0]))
    assert np.trace(moment) == pytest.approx(1.0)
    assert np.linalg.eigvalsh(moment).min() == pytest.approx(0.0)


def test_near_endpoint_slope_is_not_fuzzed_to_logarithmic_boundary():
    verdict = classify_finite_floor_shell_exponents(
        0.5, 0.5, -1.0 + 5.0e-13
    )
    assert not verdict.logarithmic_critical_boundary
    assert not verdict.critical_upper_bound_not_uniform
    assert not verdict.kinematically_admissible


def test_pressure_hessian_obeys_distance_minus_five_scaling():
    moment = np.diag([2.0, 1.0, 0.5])
    near = pressure_hessian_quadrupole(moment, np.array([2.0, 0.0, 0.0]))
    far = pressure_hessian_quadrupole(moment, np.array([4.0, 0.0, 0.0]))
    assert np.allclose(near, 32.0 * far, rtol=1.0e-14, atol=1.0e-14)


def test_moving_front_law_has_parabolic_finite_time_growth():
    law = MovingFrontLaw(initial_bandwidth=2.0, coefficient=0.25, initial_time=1.0)
    assert law.blowup_time == pytest.approx(1.5)
    for time in (1.05, 1.2, 1.4):
        bandwidth = law.bandwidth(time)
        remaining = law.blowup_time - time
        assert bandwidth * math.sqrt(2.0 * law.coefficient * remaining) == pytest.approx(
            1.0
        )
    with pytest.raises(ValueError):
        law.bandwidth(law.blowup_time)

    converted = MovingFrontLaw("2", "0.25", "1")
    assert converted.blowup_time == pytest.approx(1.5)


def test_zeno_scaling_window_meets_dissipation_and_bkm_necessary_conditions():
    report = zeno_scaling_exponents(0.5)
    assert report["energy"] == 0.0
    assert report["enstrophy"] == pytest.approx(-0.5)
    assert report["vorticity_linf"] == pytest.approx(-1.0)
    assert report["global_l3_cubed"] == "log(1/tau)"


@pytest.mark.parametrize(
    "call",
    [
        lambda: critical_staircase_metrics(-1),
        lambda: critical_staircase_metrics(1024),
        lambda: critical_staircase_metrics(2, base=1.0),
        lambda: finite_floor_shell_sums(2, spectral_slope=0.0, peak_energy=-1.0),
        lambda: pressure_hessian_quadrupole(np.eye(2), np.ones(3)),
        lambda: pressure_hessian_quadrupole(np.eye(3), np.zeros(3)),
        lambda: pressure_hessian_quadrupole(
            np.diag([2.0, -1.0, 1.0]), np.ones(3)
        ),
        lambda: normalized_pressure_stretch(np.eye(3), np.ones(3), np.zeros(3)),
        lambda: gaussian_vortex_moment(np.zeros(3)),
        lambda: MovingFrontLaw(1.0, -1.0),
        lambda: zeno_scaling_exponents(0.0),
    ],
)
def test_invalid_inputs_are_rejected(call):
    with pytest.raises(ValueError):
        call()
