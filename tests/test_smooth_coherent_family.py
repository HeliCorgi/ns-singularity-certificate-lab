"""Tests for the smoothly truncated coherent family and its continuum profile.

The tests are grouped by the claim they defend:

* the weight ``chi`` is what it is advertised to be, and is exactly rational;
* the exact moment laws hold with **zero** residual (not "small" residual);
* the exact rational lane and the dealiased FFT lane agree to roundoff;
* the continuum profile ``V`` is solenoidal and has the predicted
  ``pi^2 (v0 + yhat (yhat.v0))/|y|`` tail;
* ``curl(V . grad V)`` is nonzero, computed two independent ways
  (closed-form axisymmetric reduction vs Cartesian finite differences), and
  the discrete lattice convolution converges to the continuum coefficient.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from ns_certificate_lab.smooth_coherent_family import (
    asymptotic_curl_coefficient,
    band_waves,
    capacity_metrics,
    chi_exact,
    chi_float,
    convolution_coefficient,
    curl_coefficient,
    dealias_grid,
    exact_moment_laws,
    exact_nonlinear_mode,
    exact_nonlinear_norm_sq,
    exact_scaled_coefficients,
    meridional_profiles,
    radial_profiles,
    spectral_field,
    vector_profile,
)

V0 = (1, 2, 3)


# --------------------------------------------------------------------------
# the weight
# --------------------------------------------------------------------------


def test_chi_is_exactly_one_on_the_inner_half_and_zero_at_the_rim() -> None:
    band = 12
    for wave in band_waves(band):
        norm_sq = sum(c * c for c in wave)
        weight = chi_exact(norm_sq, band)
        assert isinstance(weight, Fraction)
        assert 0 <= weight <= 1
        if 4 * norm_sq <= band * band:
            assert weight == 1
        if norm_sq == band * band:
            assert weight == 0


def test_chi_is_nonincreasing_in_the_radius() -> None:
    band = 20
    previous = Fraction(1)
    for norm_sq in range(1, band * band + 1):
        weight = chi_exact(norm_sq, band)
        assert weight <= previous
        previous = weight


def test_chi_float_matches_the_exact_weight() -> None:
    band = 16
    for norm_sq in range(1, band * band + 1):
        exact = float(chi_exact(norm_sq, band))
        approximate = float(chi_float(math.sqrt(norm_sq) / band))
        # The float path squares a square root and then raises the result to
        # the ninth power, so it loses a few ulps relative to the rational
        # evaluation; the exact lane is the certificate lane.
        assert approximate == pytest.approx(exact, abs=1.0e-12)


def test_chi_is_c1_across_the_two_break_points() -> None:
    # A C^4 join: the first derivative must match to quadrature accuracy.
    for break_radius in (0.5, 1.0):
        step = 1.0e-5
        left = (
            float(chi_float(break_radius - step)) - float(chi_float(break_radius - 2 * step))
        ) / step
        right = (
            float(chi_float(break_radius + 2 * step)) - float(chi_float(break_radius + step))
        ) / step
        assert left == pytest.approx(0.0, abs=1.0e-6)
        assert right == pytest.approx(0.0, abs=1.0e-6)


# --------------------------------------------------------------------------
# the exact family laws
# --------------------------------------------------------------------------


@pytest.mark.parametrize("band", [4, 6, 8])
def test_moment_laws_hold_with_exactly_zero_residual(band: int) -> None:
    report = exact_moment_laws(band, V0)
    assert report["H0_law_exact"] is True
    assert report["H1_law_exact"] is True
    assert report["u_at_origin_law_exact"] is True
    assert report["H0_residual"] == "0"
    assert report["H1_residual"] == "0"
    assert report["u_at_origin_residual"] == ["0", "0", "0"]


@pytest.mark.parametrize("band", [4, 6, 8])
def test_family_is_real_and_divergence_free(band: int) -> None:
    report = exact_moment_laws(band, V0)
    assert report["real_cosine_field"] is True
    assert report["divergence_free"] is True


def test_the_origin_value_is_parallel_to_the_seed() -> None:
    report = exact_moment_laws(6, V0)
    value = [Fraction(entry) for entry in report["u_at_origin"]]
    # u_N(0) = (2/3) v0 sum chi/|k|^2, hence exactly proportional to v0.
    assert value[1] == 2 * value[0]
    assert value[2] == 3 * value[0]
    assert value[0] > 0


def test_scaled_coefficients_are_integral_and_transverse() -> None:
    scale, table = exact_scaled_coefficients(6, V0)
    assert isinstance(scale, int) and scale > 0
    assert table
    for wave, vector in table.items():
        assert all(isinstance(component, int) for component in vector)
        assert sum(c * v for c, v in zip(wave, vector)) == 0


def test_the_rim_modes_are_annihilated_by_the_weight() -> None:
    band = 5
    _, table = exact_scaled_coefficients(band, V0)
    for wave in band_waves(band):
        if sum(c * c for c in wave) == band * band:
            assert wave not in table


# --------------------------------------------------------------------------
# exact lane vs FFT lane
# --------------------------------------------------------------------------


def test_exact_and_dealiased_fft_capacities_agree() -> None:
    band = 6
    exact = exact_nonlinear_norm_sq(band, V0)
    measured = capacity_metrics(spectral_field(dealias_grid(band), band, V0), band)
    assert measured["nonlinear_sq"] == pytest.approx(
        exact["nonlinear_sq_float"], rel=1.0e-12
    )
    assert measured["K"] == pytest.approx(exact["K_float"], rel=1.0e-12)
    assert measured["N0_sq"] == pytest.approx(exact["N0_sq_float"], rel=1.0e-12)


def test_exact_single_mode_matches_the_full_fft_coefficient() -> None:
    band = 8
    wave = (4, 4, 4)
    exact = exact_nonlinear_mode(band, wave, V0)
    from ns_certificate_lab.leray_response_relay import leray_advection

    grid = dealias_grid(band)
    nonlinear = leray_advection(
        spectral_field(grid, band, V0), spectral_field(grid, band, V0)
    )
    index = tuple(c % grid for c in wave)
    measured = np.array([nonlinear[(c, *index)] for c in range(3)])
    # N_hat(k) = i P_k b_k, so the FFT coefficient is purely imaginary.
    assert np.max(np.abs(measured.real)) < 1.0e-9 * np.max(np.abs(measured))
    assert np.linalg.norm(measured.imag) == pytest.approx(
        exact["leray_b_k_norm"], rel=1.0e-10
    )


def test_dealias_grid_clears_the_repository_guard() -> None:
    for band in (3, 8, 17):
        grid = dealias_grid(band)
        assert 2 * band < grid // 2


# --------------------------------------------------------------------------
# the continuum profile
# --------------------------------------------------------------------------


@pytest.mark.parametrize("radius", [0.3, 1.0, 3.0, 9.0, 27.0])
def test_the_profile_is_solenoidal(radius: float) -> None:
    # r alpha' + 2 alpha - 2A = 0 is equivalent to div V = 0 and is an
    # independent check on g, h and their derivatives simultaneously.
    report = meridional_profiles(radius)
    scale = abs(report["alpha"]) + abs(report["A"])
    assert abs(report["divergence_residual"]) < 1.0e-11 * scale


def test_the_profile_has_the_predicted_inverse_distance_tail() -> None:
    for radius in (64.0, 96.0):
        report = vector_profile(radius)
        assert report["r_times_A"] == pytest.approx(math.pi**2, rel=2.0e-3)
        assert report["r_times_B"] == pytest.approx(math.pi**2, rel=2.0e-3)


def test_the_untruncated_closed_form_is_recovered_by_the_radial_reduction() -> None:
    # With chi == 1 the reduction gives g = 2 pi^2 / r and h = 0.  The
    # truncated g must approach 2 pi^2 / r and h must decay.
    for radius in (64.0, 96.0):
        report = radial_profiles(radius)
        assert report["g"] == pytest.approx(2.0 * math.pi**2 / radius, rel=3.0e-3)
        assert abs(report["h"]) < 5.0e-3 * report["g"]


# --------------------------------------------------------------------------
# (V-NONDEG): the decisive input
# --------------------------------------------------------------------------


@pytest.mark.parametrize("radius", [1.0, 2.0, 4.0, 8.0, 16.0])
def test_the_curl_of_the_continuum_nonlinearity_is_nonzero(radius: float) -> None:
    value = curl_coefficient(radius, V0)
    assert value < 0.0
    assert abs(value) > 1.0e-3


def test_the_curl_amplitude_approaches_its_closed_form_tail() -> None:
    for radius in (64.0, 96.0):
        assert curl_coefficient(radius, V0) == pytest.approx(
            asymptotic_curl_coefficient(radius, V0), rel=0.05
        )


def test_axisymmetric_curl_matches_cartesian_finite_differences() -> None:
    """Second, independent implementation of ``curl(V . grad V)``.

    The reduction claims ``curl(V.grad V) = Z(r) sin(theta) cos(theta) e_phi``.
    Here ``V`` is rebuilt from ``A`` and ``B`` in Cartesian coordinates and
    the curl of the advection term is taken by central differences, with the
    seed pointing along ``e_1`` so no coordinate frame is shared with the
    derivation.
    """

    seed = np.array([1.0, 0.0, 0.0])

    def field(point: np.ndarray) -> np.ndarray:
        radius = float(np.linalg.norm(point))
        profile = vector_profile(radius, order=200)
        direction = point / radius
        return profile["A"] * seed + profile["B"] * direction * float(direction @ seed)

    def advection(point: np.ndarray, step: float) -> np.ndarray:
        value = field(point)
        total = np.zeros(3)
        for axis in range(3):
            shift = np.zeros(3)
            shift[axis] = step
            total += value[axis] * (field(point + shift) - field(point - shift)) / (2 * step)
        return total

    step = 1.0e-4
    for point in (np.array([2.0, 1.0, 0.5]), np.array([4.0, -2.0, 1.0])):
        jacobian = np.zeros((3, 3))
        for axis in range(3):
            shift = np.zeros(3)
            shift[axis] = step
            jacobian[:, axis] = (
                advection(point + shift, step) - advection(point - shift, step)
            ) / (2 * step)
        measured = np.array(
            [
                jacobian[2, 1] - jacobian[1, 2],
                jacobian[0, 2] - jacobian[2, 0],
                jacobian[1, 0] - jacobian[0, 1],
            ]
        )
        radius = float(np.linalg.norm(point))
        direction = point / radius
        cosine = float(direction @ seed)
        sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
        azimuth = np.cross(seed, direction)
        azimuth /= np.linalg.norm(azimuth)
        predicted = meridional_profiles(radius, order=200)["Z"] * sine * cosine * azimuth
        assert np.linalg.norm(measured - predicted) < 1.0e-4 * np.linalg.norm(predicted)


def test_convolution_coefficient_is_converged_and_nonzero() -> None:
    zeta = np.array([0.5, 0.5, 0.5])
    coarse = convolution_coefficient(zeta, V0, 30, 30, 30)
    fine = convolution_coefficient(zeta, V0, 60, 60, 60)
    assert np.linalg.norm(fine - coarse) < 1.0e-5 * np.linalg.norm(fine)
    direction = zeta / np.linalg.norm(zeta)
    projected = fine - direction * float(direction @ fine)
    assert np.linalg.norm(projected) > 1.0


def test_discrete_lattice_convolution_converges_to_the_continuum_coefficient() -> None:
    """The Riemann-sum claim, exact on the discrete side.

    ``b_k`` is computed in exact rational arithmetic at two bands and
    Richardson-extrapolated in ``1/N``; the result must agree with the
    independently quadratured continuum coefficient ``P_zeta C(zeta)``.
    """

    zeta = np.array([0.5, 0.5, 0.5])
    coarse = exact_nonlinear_mode(16, (8, 8, 8), V0)["leray_b_k_norm"]
    fine = exact_nonlinear_mode(32, (16, 16, 16), V0)["leray_b_k_norm"]
    extrapolated = 2.0 * fine - coarse

    continuum = convolution_coefficient(zeta, V0, 60, 60, 60)
    direction = zeta / np.linalg.norm(zeta)
    projected = float(
        np.linalg.norm(continuum - direction * float(direction @ continuum))
    )

    assert coarse < fine < projected
    assert extrapolated == pytest.approx(projected, rel=5.0e-3)
