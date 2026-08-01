"""Unit checks for the continuum-to-lattice shadowing scaffolding."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest


_SPEC = importlib.util.spec_from_file_location(
    "run_poincare_shadowing",
    Path(__file__).resolve().parents[1]
    / "experiments"
    / "run_poincare_shadowing.py",
)
assert _SPEC is not None and _SPEC.loader is not None
shadowing = importlib.util.module_from_spec(_SPEC)
# ``dataclasses`` resolves annotations through ``sys.modules``, so the module
# must be registered before it is executed.
sys.modules[_SPEC.name] = shadowing
_SPEC.loader.exec_module(shadowing)


def _hat_kernel(y: np.ndarray) -> np.ndarray:
    return np.prod(np.maximum(0.0, 1.0 - np.abs(y)), axis=-1)


def _interpolate(nodes: np.ndarray, values: np.ndarray, point: np.ndarray, n: int):
    weights = _hat_kernel(n * point[None, :] - nodes)
    return weights @ values


def test_relative_width_is_exact_only_for_multiples_of_sixteen() -> None:
    for scale in (16, 32, 48, 64):
        assert shadowing.relative_width(scale) == pytest.approx(3.0 / 16.0)
    for scale in (8, 24):
        assert shadowing.relative_width(scale) != pytest.approx(3.0 / 16.0)


def test_restriction_after_interpolation_is_the_identity_on_nodes() -> None:
    """``R_N I_N = id``: the hat kernel is nodally exact."""

    n = 5
    axis = np.arange(-2, 3, dtype=np.int64)
    nodes = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    ).reshape(-1, 3)
    generator = np.random.default_rng(11)
    coefficients = generator.standard_normal(nodes.shape[0]) + 1.0j * (
        generator.standard_normal(nodes.shape[0])
    )
    profile_values = float(n) ** 2 * coefficients
    for node in nodes:
        point = node.astype(np.float64) / n
        interpolated = _interpolate(nodes, profile_values, point, n)
        restricted = float(n) ** (-2) * interpolated
        index = int(np.flatnonzero((nodes == node).all(axis=1))[0])
        assert restricted == pytest.approx(coefficients[index], abs=1e-12)


def test_interpolation_defect_on_a_c1_profile_is_first_order() -> None:
    """``|| I_N R_N U - U ||_inf <= sqrt(3) ||grad U||_inf / N`` and decays."""

    def profile(point: np.ndarray) -> np.ndarray:
        return np.sin(2.0 * point[..., 0]) + 0.5 * np.cos(
            1.5 * point[..., 1] + point[..., 2]
        )

    gradient_bound = 2.0 + 0.5 * math.sqrt(2.0) * 1.5
    previous = None
    for n in (8, 16, 32):
        axis = np.arange(-n // 4, n // 4 + 1, dtype=np.int64)
        nodes = np.stack(
            np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
        ).reshape(-1, 3)
        values = profile(nodes.astype(np.float64) / n)
        generator = np.random.default_rng(7)
        samples = generator.uniform(
            -0.2, 0.2, size=(200, 3)
        )
        defect = 0.0
        for point in samples:
            interpolated = _interpolate(nodes, values, point, n)
            defect = max(defect, abs(interpolated - profile(point)))
        assert defect <= math.sqrt(3.0) * gradient_bound / n
        if previous is not None:
            assert defect < previous
        previous = defect


def test_log_log_fit_recovers_a_planted_exponent() -> None:
    abscissae = [16.0, 32.0, 48.0, 64.0]
    values = [3.0 * value ** (-1.75) for value in abscissae]
    fit = shadowing.log_log_fit(abscissae, values)
    assert fit["sigma"] == pytest.approx(1.75, rel=1e-9)
    assert fit["r_squared"] == pytest.approx(1.0, abs=1e-12)


def test_log_log_fit_refuses_fewer_than_three_points() -> None:
    fit = shadowing.log_log_fit([16.0, 32.0], [1.0, 0.5])
    assert fit["sigma"] is None
    assert fit["r_squared"] is None
    assert fit["points"] == 2


def test_perturbation_is_real_divergence_free_and_correctly_sized() -> None:
    grid = 16
    mask = np.zeros((grid, grid, grid), dtype=np.bool_)
    frequency = np.rint(
        np.fft.fftfreq(grid, d=1.0 / grid)
    ).astype(np.int64)
    kx, ky, kz = np.meshgrid(frequency, frequency, frequency, indexing="ij")
    mask[(np.abs(kx) <= 3) & (np.abs(ky) <= 3) & (np.abs(kz) <= 3)] = True
    reference = np.zeros((3, grid, grid, grid), dtype=np.complex128)
    reference[0, 1, 0, 0] = 1.0
    reference[0, grid - 1, 0, 0] = 1.0
    generator = np.random.default_rng(3)
    delta = shadowing.divergence_free_perturbation(
        generator, reference, mask, grid
    )
    ratio = float(np.linalg.norm(delta)) / float(np.linalg.norm(reference))
    assert ratio == pytest.approx(shadowing.PERTURBATION_RELATIVE_SIZE, rel=1e-12)
    import ns_certificate_lab.renormalized_cascade as cascade

    assert cascade._divergence_defect(delta) < 1e-12
    assert cascade._reality_defect(delta) < 1e-12
    assert np.all(delta[:, 0, 0, 0] == 0.0)


def test_doubling_pairing_selects_the_common_child_band() -> None:
    """``fine[::2]`` on ``|o| <= 2m`` is exactly the image of ``|o| <= m``."""

    m = 3
    axis = np.arange(-2 * m, 2 * m + 1)
    block = np.zeros((3, axis.size, axis.size, axis.size), dtype=np.complex128)
    for i, x in enumerate(axis):
        for j, y in enumerate(axis):
            for k, z in enumerate(axis):
                block[:, i, j, k] = complex(x, 10 * y + 100 * z)
    paired = block[:, ::2, ::2, ::2]
    assert paired.shape == (3, 2 * m + 1, 2 * m + 1, 2 * m + 1)
    coarse_axis = np.arange(-m, m + 1)
    for i, x in enumerate(coarse_axis):
        for j, y in enumerate(coarse_axis):
            for k, z in enumerate(coarse_axis):
                assert paired[0, i, j, k] == complex(2 * x, 20 * y + 200 * z)
