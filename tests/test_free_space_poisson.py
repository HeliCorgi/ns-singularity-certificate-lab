from __future__ import annotations

import math

import numpy as np
import pytest

from ns_certificate_lab.bessel_reference import exp_scaled_i1, exp_scaled_k1
from ns_certificate_lab.free_space_poisson import (
    _bessel_product_i1k1,
    radial_green_matrix,
    solve_l5_free_radial_periodized_z,
)


def _relative_l2(error: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(error.ravel()) / np.linalg.norm(reference.ravel()))


def test_cephes_scaled_bessel_product_matches_scipy_when_available() -> None:
    scipy_special = pytest.importorskip("scipy.special")
    k = 7.3
    r_less = np.array([0.0, 1.0e-4, 0.03, 0.4, 2.0, 9.0])
    r_greater = np.array([0.2, 0.5, 0.4, 1.2, 5.0, 12.0])
    approximate = _bessel_product_i1k1(k, r_less, r_greater, backend="cephes")
    reference = _bessel_product_i1k1(k, r_less, r_greater, backend="scipy")
    assert np.allclose(approximate, reference, rtol=3.0e-7, atol=1.0e-12)
    assert np.all(np.isfinite(approximate))
    assert scipy_special is not None


def test_cephes_scaled_bessel_product_matches_the_in_repo_oracle() -> None:
    """Same check as the scipy test, against this repository's own oracle.

    SciPy is deliberately absent from this project's dependency set, so the
    test above always skips here and its verification would be lost.
    :mod:`ns_certificate_lab.bessel_reference` is the scipy-free ``I``/``K``
    oracle built for E-33 (series ``I_nu``, integral-representation ``K_nu``,
    Wronskian-validated); it is an independent implementation of the same
    special functions, so it can genuinely falsify the cephes polynomials.

    ``exp_scaled_i1(x) = e^{-x} I_1(x)`` and ``exp_scaled_k1(x) = e^{x}
    K_1(x)``, hence the kernel is ``scaled_i(a) scaled_k(b) e^{a-b}/(r_< r_>)``
    with ``a = k r_<``, ``b = k r_>``; at ``r_< = 0`` the module uses the exact
    limit ``I_1(k r)/r -> k/2``, which is reproduced here independently.
    """

    k = 7.3
    r_less = np.array([0.0, 1.0e-4, 0.03, 0.4, 2.0, 9.0])
    r_greater = np.array([0.2, 0.5, 0.4, 1.2, 5.0, 12.0])
    approximate = _bessel_product_i1k1(k, r_less, r_greater, backend="cephes")

    expected = np.empty_like(approximate)
    for index, (rl, rg) in enumerate(zip(r_less, r_greater)):
        b = k * float(rg)
        scaled_k = exp_scaled_k1(b)
        if rl == 0.0:
            expected[index] = 0.5 * k * scaled_k * math.exp(-b) / rg
        else:
            a = k * float(rl)
            expected[index] = (
                exp_scaled_i1(a) * scaled_k * math.exp(a - b) / (rl * rg)
            )

    assert np.allclose(approximate, expected, rtol=3.0e-7, atol=1.0e-12)
    assert np.all(np.isfinite(approximate))
    # The oracle is an independent path: it must not reproduce the cephes
    # polynomials bit for bit, or it would not be able to falsify them.
    assert not np.array_equal(approximate, expected)


def test_zero_mode_green_kernel_has_exact_formula() -> None:
    r = np.linspace(0.0, 2.0, 9)
    green = radial_green_matrix(r, 0.0)
    for i, ri in enumerate(r):
        for j, rj in enumerate(r):
            if ri == 0.0 and rj == 0.0:
                assert green[i, j] == 0.0
            else:
                assert green[i, j] == pytest.approx(0.5 / max(ri, rj) ** 2)
    assert np.allclose(green, green.T)


def test_nonzero_green_kernel_is_symmetric_and_positive() -> None:
    r = np.linspace(0.0, 3.0, 31)
    green = radial_green_matrix(r, 2.0 * math.pi)
    assert np.allclose(green, green.T, rtol=1.0e-13, atol=1.0e-13)
    assert np.all(green[1:, 1:] > 0.0)
    assert np.all(np.isfinite(green))


def test_periodic_single_mode_gaussian_manufactured_solution() -> None:
    radial_points = 161
    z_points = 96
    radial_max = 6.0
    period = 2.0 * math.pi
    a = 0.8
    mode = 3

    r = np.linspace(0.0, radial_max, radial_points)
    z = np.arange(z_points) * period / z_points
    rr, zz = np.meshgrid(r, z, indexing="ij")
    exact = np.exp(-a * rr**2) * np.cos(mode * zz)
    source = (8.0 * a - 4.0 * a * a * rr**2 + mode**2) * exact

    computed, diagnostics = solve_l5_free_radial_periodized_z(
        source,
        r,
        period / z_points,
        bessel_backend="cephes",
    )
    interior = r <= 3.5
    relative = _relative_l2(computed[interior] - exact[interior], exact[interior])
    assert relative < 1.8e-3
    assert diagnostics.pad_factor == 1
    assert diagnostics.padded_period == pytest.approx(period)


def test_zero_mode_gaussian_manufactured_solution() -> None:
    radial_points = 241
    z_points = 16
    radial_max = 7.0
    a = 0.65
    r = np.linspace(0.0, radial_max, radial_points)
    rr = r[:, None]
    exact_radial = np.exp(-a * rr**2)
    exact = np.repeat(exact_radial, z_points, axis=1)
    source = (8.0 * a - 4.0 * a * a * rr**2) * exact

    computed, _ = solve_l5_free_radial_periodized_z(source, r, dz=0.3)
    interior = r <= 4.0
    relative = _relative_l2(computed[interior] - exact[interior], exact[interior])
    assert relative < 8.0e-4


def test_gaussian_whole_space_padding_reduces_period_image_error() -> None:
    radial_points = 181
    z_points = 80
    radial_max = 6.0
    z_half_width = 2.0
    a = 0.9

    r = np.linspace(0.0, radial_max, radial_points)
    dz = 2.0 * z_half_width / z_points
    z = (np.arange(z_points) - z_points // 2) * dz
    rr, zz = np.meshgrid(r, z, indexing="ij")
    exact = np.exp(-a * (rr**2 + zz**2))
    source = (10.0 * a - 4.0 * a * a * (rr**2 + zz**2)) * exact

    errors = []
    for pad_factor in (1, 2, 4):
        computed, _ = solve_l5_free_radial_periodized_z(
            source,
            r,
            dz,
            pad_factor=pad_factor,
        )
        interior = (rr <= 3.2) & (np.abs(zz) <= 1.0)
        errors.append(_relative_l2((computed - exact)[interior], exact[interior]))

    assert errors[1] < 0.5 * errors[0]
    assert errors[2] < errors[1]
    assert errors[2] < 4.5e-3


def test_input_validation() -> None:
    r = np.linspace(0.0, 1.0, 9)
    source = np.zeros((9, 16))
    with pytest.raises(ValueError):
        solve_l5_free_radial_periodized_z(source, r, dz=-1.0)
    with pytest.raises(ValueError):
        solve_l5_free_radial_periodized_z(source, r, dz=0.1, pad_factor=0)
    with pytest.raises(ValueError):
        solve_l5_free_radial_periodized_z(source[:, :3], r, dz=0.1)
    bad = source.copy()
    bad[1, 1] = np.nan
    with pytest.raises(ValueError):
        solve_l5_free_radial_periodized_z(bad, r, dz=0.1)
