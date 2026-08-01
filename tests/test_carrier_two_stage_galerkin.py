from __future__ import annotations

import json
import math

import numpy as np
import pytest

import ns_certificate_lab.carrier_two_stage_galerkin as carrier
from ns_certificate_lab.carrier_two_stage_galerkin import (
    build_partial_carrier_parent,
    run_carrier_two_stage_galerkin,
)
from ns_certificate_lab.leray_response_relay import mean_energy


def test_four_exact_parents_are_real_divergence_free_and_normalized() -> None:
    scale = 2
    grid_size = 32
    field = build_partial_carrier_parent(grid_size, scale=scale)
    assert mean_energy(field) == pytest.approx(0.5, rel=2.0e-15)

    support = np.any(np.abs(field) > 0.0, axis=0)
    assert np.count_nonzero(support) == 8
    physical = np.fft.ifftn(field, axes=(1, 2, 3)) * grid_size**3
    assert np.max(np.abs(physical.imag)) < 2.0e-15

    frequency = np.fft.fftfreq(grid_size, d=1.0 / grid_size)
    wave = np.stack(
        np.meshgrid(frequency, frequency, frequency, indexing="ij"), axis=0
    )
    assert np.linalg.norm(1.0j * np.sum(wave * field, axis=0)) < 2.0e-15

    amplitude = math.sqrt(0.5)
    a1 = tuple(value % grid_size for value in (2, 2, 0))
    b1 = tuple(value % grid_size for value in (2, 0, 2))
    a2 = tuple(value % grid_size for value in (0, 2, -2))
    b2 = tuple(value % grid_size for value in (2, 0, -2))
    assert field[(slice(None), *a1)] == pytest.approx(
        np.array((0.0, 0.0, -0.5j * amplitude)), abs=2.0e-15
    )
    assert field[(slice(None), *b1)] == pytest.approx(
        np.array((0.0, 0.5 * amplitude, 0.0)), abs=2.0e-15
    )
    assert field[(slice(None), *a2)] == pytest.approx(
        np.array((0.5 * amplitude, 0.0, 0.0)), abs=2.0e-15
    )
    assert field[(slice(None), *b2)] == pytest.approx(
        np.array((0.0, 0.5 * amplitude, 0.0)), abs=2.0e-15
    )


def test_two_stage_run_separates_cross_talk_and_has_exact_shell_budget() -> None:
    result = run_carrier_two_stage_galerkin(
        scale=2, grid_size=32, steps=8
    )
    groups = dict(result.mode_groups)
    assert groups["first_child_one"] == ((-4, -2, -2), (4, 2, 2))
    assert groups["first_child_two"] == ((-2, -2, 4), (2, 2, -4))
    assert groups["cross_talk_one"] == ((-4, -2, 2), (4, 2, -2))
    assert groups["cross_talk_two"] == ((-2, 2, -4), (2, -2, 4))
    assert groups["grandchild_sum"] == ((-6, -4, 2), (6, 4, -2))
    assert groups["grandchild_difference"] == ((-2, 0, -6), (2, 0, 6))
    assert len(set().union(*(set(waves) for waves in groups.values()))) == 20

    assert result.initial_parent_energy == pytest.approx(0.5, rel=2.0e-15)
    assert result.initial_named_nonparent_energy == 0.0
    assert result.dealias_margin == 4
    assert result.final_first_child_energy_ratio > 0.0
    assert result.final_cross_talk_energy_ratio > 0.9 * (
        result.final_first_child_energy_ratio
    )
    assert result.final_grandchild_energy_ratio > 0.0
    assert result.history[-1]["grandchild_sum_energy_ratio"] > (
        result.history[-1]["grandchild_difference_energy_ratio"]
    )
    assert result.grandchild_binary64_resolved
    assert result.grandchild_resolution_margin > 1.0e6
    assert result.final_total_energy_ratio < 1.0
    assert result.maximum_energy_increase < 1.0e-14
    assert result.maximum_shell_budget_residual < 3.0e-15
    assert result.initial_energy_identity_residual < 1.0e-14
    assert result.final_energy_identity_residual < 1.0e-14
    assert result.final_reality_defect < 1.0e-14
    assert result.final_divergence_defect < 1.0e-14

    for record in result.history:
        classified = (
            record["parent_energy_ratio"]
            + record["first_child_energy_ratio"]
            + record["cross_talk_energy_ratio"]
            + record["grandchild_energy_ratio"]
            + record["remainder_energy_ratio"]
        )
        assert classified == pytest.approx(
            record["total_energy_ratio"], abs=2.0e-15
        )


def test_rk4_time_refinement_stabilizes_even_the_small_grandchild_signal() -> None:
    results = {
        steps: run_carrier_two_stage_galerkin(
            scale=2, grid_size=32, steps=steps
        )
        for steps in (2, 4, 8)
    }
    for attribute in (
        "final_total_energy_ratio",
        "final_first_child_energy_ratio",
        "final_cross_talk_energy_ratio",
        "final_grandchild_energy_ratio",
    ):
        coarse_difference = abs(
            getattr(results[2], attribute) - getattr(results[4], attribute)
        )
        fine_difference = abs(
            getattr(results[4], attribute) - getattr(results[8], attribute)
        )
        assert coarse_difference / fine_difference > 8.0

    assert results[8].final_grandchild_energy_ratio == pytest.approx(
        results[4].final_grandchild_energy_ratio, rel=2.0e-5
    )


def test_initial_grandchildren_have_no_resolved_direct_rhs() -> None:
    scale = 2
    grid_size = 32
    cutoff = 3 * scale
    parent = build_partial_carrier_parent(grid_size, scale=scale)
    kx, ky, kz = carrier._mesh(grid_size)
    wave_squared = kx * kx + ky * ky + kz * kz
    rhs = carrier._rhs(
        parent,
        viscosity=1.0 / 40.0,
        wave_squared=wave_squared,
        galerkin_mask=carrier._galerkin_mask(grid_size, cutoff),
    )
    groups = carrier._carrier_groups(scale)
    grandchild_mask = carrier._mode_mask(
        grid_size,
        groups["grandchild_sum"] + groups["grandchild_difference"],
    )
    full_rhs_energy = 0.5 * float(np.vdot(rhs, rhs).real)
    noise_energy = carrier._masked_energy(rhs, grandchild_mask)
    assert noise_energy / full_rhs_energy < 1.0e-30


def test_result_payload_is_json_safe_and_explicitly_rejected() -> None:
    result = run_carrier_two_stage_galerkin(
        scale=2, grid_size=32, steps=2
    )
    payload = result.as_dict()
    json.dumps(payload, allow_nan=False)
    assert payload["classification"] == "partial_rejected_cross_talk"
    assert payload["dealias_verified"] is True
    assert "pathway-contaminated" in payload["interpretation"]
    assert "NOT A RELAY" in payload["status"]
    assert len(payload["history"]) == 3


@pytest.mark.parametrize(
    "kwargs",
    [
        {"scale": 0},
        {"scale": 2, "grid_size": 24},
        {"scale": 2, "grid_size": 32, "galerkin_cutoff": 5},
        {"scale": 2, "grid_size": 32, "steps": 0},
        {"scale": 2, "grid_size": 32, "tau": 0.0},
        {"scale": 2, "grid_size": 32, "time_multiples": 0.0},
        {"scale": 2, "grid_size": 32, "viscosity": 0.0},
        {"scale": 2, "grid_size": 32, "energy_constant": 0.0},
    ],
)
def test_invalid_time_or_padding_geometry_is_rejected(kwargs) -> None:
    with pytest.raises(ValueError):
        run_carrier_two_stage_galerkin(**kwargs)
