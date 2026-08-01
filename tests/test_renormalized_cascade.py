"""Tests for the renormalized Galerkin cascade pilot."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ns_certificate_lab.leray_response_relay import mean_energy
from ns_certificate_lab.mesoscopic_galerkin import build_angle_box_parent
from ns_certificate_lab.renormalized_cascade import (
    CascadeConfig,
    doubling_pullback,
    linf_shell_energies,
    run_renormalized_cascade,
)


def _integer_frequencies(grid_size: int) -> np.ndarray:
    return np.rint(
        np.fft.fftfreq(grid_size, d=1.0 / grid_size)
    ).astype(np.int64)


def _mode_index(k: tuple[int, int, int], grid_size: int):
    return tuple(int(component) % grid_size for component in k)


def test_pullback_moves_even_modes_and_scales_by_four() -> None:
    grid = 32
    field = np.zeros((3, grid, grid, grid), dtype=np.complex128)
    wave = (6, -4, 2)
    coefficient = np.array([2.0 - 1.0j, 1.0j, 3.0], dtype=np.complex128)
    field[(slice(None), *_mode_index(wave, grid))] = coefficient
    field[(slice(None), *_mode_index((-6, 4, -2), grid))] = np.conjugate(
        coefficient
    )
    pulled = doubling_pullback(field)
    target = pulled[(slice(None), *_mode_index((3, -2, 1), grid))]
    np.testing.assert_allclose(target, 4.0 * coefficient)
    reflected = pulled[(slice(None), *_mode_index((-3, 2, -1), grid))]
    np.testing.assert_allclose(reflected, np.conjugate(4.0 * coefficient))
    pulled[(slice(None), *_mode_index((3, -2, 1), grid))] = 0.0
    pulled[(slice(None), *_mode_index((-3, 2, -1), grid))] = 0.0
    assert float(np.linalg.norm(pulled)) == 0.0


def test_pullback_drops_odd_modes_entirely() -> None:
    grid = 32
    field = np.zeros((3, grid, grid, grid), dtype=np.complex128)
    field[(slice(None), *_mode_index((5, 3, -1), grid))] = 1.0 + 0.5j
    pulled = doubling_pullback(field)
    assert float(np.linalg.norm(pulled)) == 0.0


def test_pullback_preserves_divergence_free_structure() -> None:
    grid = 32
    parent, _ = build_angle_box_parent(grid, scale=4, width=2)
    pulled = doubling_pullback(parent)
    frequency = _integer_frequencies(grid).astype(np.float64)
    mesh = np.meshgrid(frequency, frequency, frequency, indexing="ij")
    divergence = 1.0j * sum(
        mesh[axis] * pulled[axis] for axis in range(3)
    )
    assert float(np.linalg.norm(divergence)) <= 1e-12 * max(
        float(np.linalg.norm(pulled)), 1.0
    )


def test_pullback_energy_matches_even_mode_bookkeeping() -> None:
    grid = 32
    parent, _ = build_angle_box_parent(grid, scale=4, width=2)
    frequency = _integer_frequencies(grid)
    even = frequency % 2 == 0
    even_mask = (
        even[:, None, None] & even[None, :, None] & even[None, None, :]
    )
    even_energy = mean_energy(
        np.asarray(parent * even_mask[None, ...], dtype=np.complex128)
    )
    pulled_energy = mean_energy(doubling_pullback(parent))
    assert pulled_energy == pytest.approx(16.0 * even_energy, rel=1e-12)


def test_linf_shell_energies_partition_total() -> None:
    grid = 32
    parent, _ = build_angle_box_parent(grid, scale=4, width=2)
    shells = linf_shell_energies(parent, (4, 8))
    assert len(shells) == 3
    assert sum(shells) == pytest.approx(mean_energy(parent), rel=1e-12)


def test_cascade_smoke_runs_and_records_consistent_diagnostics() -> None:
    config = CascadeConfig(
        scale=4,
        width=2,
        grid_size=64,
        stages=2,
        base_steps=4,
        max_steps=64,
        drop_below=3,
    )
    result = run_renormalized_cascade(config)
    assert result.completed_stages == 2
    assert result.terminated_reason == "completed"
    for record in result.stage_records:
        assert record.pullback_gain > 0.0
        assert 0.0 <= record.dropped_sea_fraction <= 1.0
        assert record.reality_defect <= 1e-10
        assert record.divergence_defect <= 1e-10
        assert record.evolved_energy_ratio <= 1.0 + 1e-9
    payload = result.as_dict()
    assert payload["status"].startswith("BINARY64")
    assert payload["config"]["cutoff"] == 10


def test_cascade_rejects_undersized_grid() -> None:
    config = CascadeConfig(scale=4, width=2, grid_size=32, stages=1)
    with pytest.raises(ValueError):
        run_renormalized_cascade(config)


def test_config_validation_rejects_bad_values() -> None:
    with pytest.raises(ValueError):
        CascadeConfig(tau=-1.0).validate()
    with pytest.raises(ValueError):
        CascadeConfig(stages=0).validate()
    with pytest.raises(ValueError):
        CascadeConfig(drop_below=0).validate()
    with pytest.raises(ValueError):
        CascadeConfig(max_steps=1, base_steps=8).validate()
