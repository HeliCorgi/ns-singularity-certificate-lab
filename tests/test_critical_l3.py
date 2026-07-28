"""Tests for critical-L3 scaling and multiscale diagnostics."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ns_certificate_lab.critical_l3 import (
    analyze_critical_l3,
    anisotropic_l3_scaling_factor,
    axisymmetric_l3_cubed,
    critical_growth_series,
    rescaled_l3_cubed,
    state_critical_l3,
)
from ns_certificate_lab.grid import AxisymmetricGrid


def _grid(*, radial_scale: float = 1.0, axial_scale: float = 1.0):
    return AxisymmetricGrid.uniform(
        nr=65,
        nz=128,
        r_max=2.0 * radial_scale,
        z_min=0.0,
        z_max=4.0 * axial_scale,
        periodic_z=True,
    )


def _base_components(grid: AxisymmetricGrid):
    r, z = grid.mesh()
    phase = 2.0 * math.pi * z / float(grid.z_period)
    envelope = np.exp(-2.0 * (r / grid.r[-1]) ** 2)
    return (
        0.20 * r * envelope * np.sin(phase),
        0.75 * r * envelope * (1.0 + 0.15 * np.cos(phase)),
        0.35 * envelope * np.cos(phase),
    )



def test_constant_velocity_matches_exact_cylinder_volume() -> None:
    grid = _grid()
    zeros = np.zeros(grid.shape)
    ones = np.ones(grid.shape)
    measured = axisymmetric_l3_cubed(grid, zeros, zeros, ones)
    exact_volume = math.pi * float(grid.r[-1]) ** 2 * float(grid.z_period)
    assert measured == pytest.approx(exact_volume, rel=2.0e-15)

def test_anisotropic_change_of_variables_is_exact_on_matched_grids() -> None:
    base_grid = _grid()
    base = _base_components(base_grid)
    base_l3 = axisymmetric_l3_cubed(base_grid, *base)

    amplitude = 3.5
    radial_scale = 0.4
    axial_scale = 1.7
    physical_grid = _grid(radial_scale=radial_scale, axial_scale=axial_scale)
    physical = tuple(amplitude * component for component in base)
    measured = axisymmetric_l3_cubed(physical_grid, *physical)
    predicted = rescaled_l3_cubed(
        base_l3,
        amplitude=amplitude,
        radial_scale=radial_scale,
        axial_scale=axial_scale,
    )
    assert measured == pytest.approx(predicted, rel=2.0e-14, abs=1.0e-14)


def test_standard_isotropic_scaling_preserves_l3() -> None:
    for scale in (0.125, 0.5, 2.0, 8.0):
        assert anisotropic_l3_scaling_factor(1.0 / scale, scale, scale) == pytest.approx(
            1.0, rel=2.0e-15
        )


def test_shell_partition_sums_to_total_and_is_translation_invariant() -> None:
    grid = _grid()
    components = _base_components(grid)
    report = analyze_critical_l3(grid, *components, shell_count=9)
    rolled = tuple(np.roll(component, 19, axis=1) for component in components)
    translated = analyze_critical_l3(grid, *rolled, shell_count=9)

    assert sum(report.shell_masses) == pytest.approx(report.l3_cubed, rel=2.0e-15)
    assert sum(report.shell_fractions) == pytest.approx(1.0, rel=2.0e-15)
    assert translated.l3_norm == pytest.approx(report.l3_norm, rel=2.0e-15)
    assert translated.radial_scale == pytest.approx(report.radial_scale, rel=2.0e-15)
    assert translated.axial_scale == pytest.approx(report.axial_scale, rel=2.0e-13)
    assert translated.critical_scaling_product == pytest.approx(
        report.critical_scaling_product, rel=2.0e-13
    )
    assert translated.shell_fractions == pytest.approx(
        report.shell_fractions, rel=2.0e-13, abs=2.0e-15
    )


def test_state_adapter_uses_physical_swirl_r_times_u1() -> None:
    grid = _grid()
    u_r, u_theta, u_z = _base_components(grid)
    u1 = np.zeros_like(u_theta)
    u1[1:] = u_theta[1:] / grid.r[1:, None]
    direct = analyze_critical_l3(grid, u_r, u_theta, u_z)
    adapted = state_critical_l3(grid, u1=u1, u_r=u_r, u_z=u_z)
    assert adapted.l3_cubed == pytest.approx(direct.l3_cubed, rel=2.0e-15)
    assert adapted.critical_scaling_product == pytest.approx(
        direct.critical_scaling_product, rel=2.0e-15
    )


def test_growth_series_reports_known_exponential_rates() -> None:
    grid = _grid()
    components = _base_components(grid)
    times = (0.0, 0.25, 0.75)
    reports = [
        analyze_critical_l3(
            grid,
            *(math.exp(time) * component for component in components),
        )
        for time in times
    ]
    growth = critical_growth_series(times, reports)
    # L3 norm scales linearly with amplitude; Q scales cubically.
    assert growth.log_l3_slopes == pytest.approx((1.0, 1.0), rel=2.0e-13)
    assert growth.log_scaling_product_slopes == pytest.approx(
        (3.0, 3.0), rel=2.0e-13
    )


def test_rejects_invalid_inputs() -> None:
    grid = _grid()
    components = _base_components(grid)
    with pytest.raises(ValueError, match="positive and finite"):
        anisotropic_l3_scaling_factor(0.0, 1.0, 1.0)
    with pytest.raises(ValueError, match="shell_count"):
        analyze_critical_l3(grid, *components, shell_count=1)
    with pytest.raises(ValueError, match="nonzero velocity"):
        analyze_critical_l3(grid, *(np.zeros(grid.shape) for _ in range(3)))
