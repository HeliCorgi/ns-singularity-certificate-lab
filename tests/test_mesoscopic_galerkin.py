from __future__ import annotations

import math

import numpy as np
import pytest

import ns_certificate_lab.mesoscopic_galerkin as mesoscopic
from ns_certificate_lab.leray_response_relay import (
    gradient_l2_squared,
    mean_energy,
)
from ns_certificate_lab.mesoscopic_galerkin import (
    build_angle_box_parent,
    run_small_mesoscopic_galerkin,
)


def test_angle_box_parent_is_real_divergence_free_and_has_empty_child() -> None:
    parent, child = build_angle_box_parent(64, scale=4, width=2)
    assert mean_energy(parent) == pytest.approx(0.25, rel=1.0e-14)
    assert mean_energy(parent * child[None, ...]) == 0.0
    physical = np.fft.ifftn(parent, axes=(1, 2, 3)) * 64**3
    assert np.max(np.abs(physical.imag)) < 1.0e-12

    frequency = np.fft.fftfreq(64, d=1.0 / 64)
    wave = np.stack(
        np.meshgrid(frequency, frequency, frequency, indexing="ij"), axis=0
    )
    divergence = 1.0j * np.sum(wave * parent, axis=0)
    assert np.linalg.norm(divergence) < 1.0e-12

    parent_support = np.any(np.abs(parent) > 0.0, axis=0)
    assert np.count_nonzero(parent_support) == 4 * 3**3
    assert np.count_nonzero(child) == 2 * 3**3
    outer_sumset_wave = 4 * mesoscopic._C + np.array((2, 0, 0))
    assert not child[tuple(int(value) % 64 for value in outer_sumset_wave)]


def test_full_galerkin_generates_child_and_refines_in_time() -> None:
    coarse = run_small_mesoscopic_galerkin(steps=8)
    fine = run_small_mesoscopic_galerkin(steps=16)
    assert fine.initial_child_energy == 0.0
    assert fine.frozen_child_energy_ratio > 0.0
    assert fine.full_child_energy_ratio > 0.0
    assert fine.full_to_frozen_ratio > 0.0
    assert fine.dealias_margin == 12
    assert fine.child_core_half_width == 1
    assert fine.child_full_sumset_half_width == 2
    assert fine.child_core_mode_count == 2 * 3**3
    assert fine.child_full_sumset_mode_count == 2 * 5**3
    assert 0.0 < fine.child_core_forcing_fraction < 1.0
    assert fine.final_total_energy_ratio < 1.0
    assert fine.maximum_energy_increase < 1.0e-12
    assert fine.final_reality_defect < 1.0e-12
    assert fine.final_divergence_defect < 1.0e-12
    assert fine.full_child_energy_ratio == pytest.approx(
        coarse.full_child_energy_ratio, rel=2.0e-5
    )


def test_single_triad_rhs_sign_energy_law_and_frozen_heat_factor() -> None:
    scale = 2
    grid_size = 32
    viscosity = 1.0 / 40.0
    tau = 0.25
    parent, _ = build_angle_box_parent(
        grid_size, scale=scale, width=1
    )
    kx, ky, kz = mesoscopic._mesh(grid_size)
    wave_squared = kx * kx + ky * ky + kz * kz
    mask = mesoscopic._galerkin_mask(grid_size, 2 * scale)
    rhs = mesoscopic._rhs(
        parent,
        viscosity=viscosity,
        wave_squared=wave_squared,
        galerkin_mask=mask,
    )
    child_wave = tuple(int(value) % grid_size for value in scale * mesoscopic._C)
    assert rhs[(slice(None), *child_wave)] == pytest.approx(
        np.array((1.0, -1.0, -1.0)) / 3.0, abs=2.0e-15
    )
    assert float(np.vdot(parent, rhs).real) == pytest.approx(
        -viscosity * gradient_l2_squared(parent), abs=2.0e-15
    )

    result = run_small_mesoscopic_galerkin(
        scale=scale,
        width=1,
        grid_size=grid_size,
        viscosity=viscosity,
        tau=tau,
        steps=4,
    )
    child_wave_squared = 6.0 * scale * scale
    final_time = tau / (scale * scale)
    heat_factor = -math.expm1(
        -viscosity * child_wave_squared * final_time
    ) / (viscosity * child_wave_squared)
    expected_frozen_ratio = scale * heat_factor * heat_factor / 3.0
    assert result.frozen_child_energy_ratio == pytest.approx(
        expected_frozen_ratio, rel=2.0e-14
    )
    assert result.child_core_forcing_fraction == pytest.approx(1.0)

    tiny_tau = 1.0e-12
    tiny = run_small_mesoscopic_galerkin(
        scale=scale,
        width=1,
        grid_size=grid_size,
        viscosity=viscosity,
        tau=tiny_tau,
        steps=1,
    )
    tiny_time = tiny_tau / (scale * scale)
    tiny_heat_factor = -math.expm1(
        -viscosity * child_wave_squared * tiny_time
    ) / (viscosity * child_wave_squared)
    assert tiny.frozen_child_energy_ratio == pytest.approx(
        scale * tiny_heat_factor * tiny_heat_factor / 3.0,
        rel=2.0e-14,
    )


def test_rhs_projects_ambient_state_to_the_galerkin_cube() -> None:
    grid_size = 32
    parent, _ = build_angle_box_parent(
        grid_size, scale=2, width=1
    )
    kx, ky, kz = mesoscopic._mesh(grid_size)
    wave_squared = kx * kx + ky * ky + kz * kz
    mask = mesoscopic._galerkin_mask(grid_size, 4)
    clean = mesoscopic._rhs(
        parent,
        viscosity=1.0 / 40.0,
        wave_squared=wave_squared,
        galerkin_mask=mask,
    )
    polluted = parent.copy()
    polluted[1, 5, 0, 0] = 0.1
    polluted[1, -5 % grid_size, 0, 0] = 0.1
    projected = mesoscopic._rhs(
        polluted,
        viscosity=1.0 / 40.0,
        wave_squared=wave_squared,
        galerkin_mask=mask,
    )
    assert projected == pytest.approx(clean, abs=2.0e-15)


def test_rk4_time_refinement_shows_high_order_convergence() -> None:
    values = {
        steps: run_small_mesoscopic_galerkin(
            scale=2,
            width=1,
            grid_size=32,
            tau=1.0,
            steps=steps,
        ).full_child_energy_ratio
        for steps in (2, 4, 8, 16)
    }
    coarse_difference = abs(values[2] - values[4])
    middle_difference = abs(values[4] - values[8])
    fine_difference = abs(values[8] - values[16])
    assert coarse_difference / middle_difference > 8.0
    assert middle_difference / fine_difference > 8.0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"scale": 4, "width": 3},
        {"grid_size": 32},
        {"scale": 2, "width": 1, "grid_size": 17},
        {"steps": 0},
        {"tau": 0.0},
    ],
)
def test_small_galerkin_rejects_invalid_geometry_or_resolution(kwargs) -> None:
    with pytest.raises(ValueError):
        run_small_mesoscopic_galerkin(**kwargs)
