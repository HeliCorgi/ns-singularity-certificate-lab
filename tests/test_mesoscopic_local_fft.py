import math

import numpy as np
import pytest

from ns_certificate_lab.mesoscopic_cloud_scaling import (
    MesoscopicCloudConfig,
    build_sparse_parent,
    exact_sparse_leray_convolution,
    measure_mesoscopic_cloud,
)
from ns_certificate_lab.mesoscopic_local_fft import (
    local_fft_leray_coefficients,
    measure_local_fft_cloud,
    measure_local_fft_mesoscopic_metrics,
)


def _assert_all_coefficients_match(config: MesoscopicCloudConfig) -> None:
    parent = build_sparse_parent(config)
    sparse = exact_sparse_leray_convolution(parent)
    local_fft = local_fft_leray_coefficients(config)
    assert local_fft.keys() == sparse.keys()
    scale = max(float(np.linalg.norm(value)) for value in sparse.values())
    for wave in sorted(sparse):
        assert local_fft[wave] == pytest.approx(
            sparse[wave], rel=2.0e-13, abs=2.0e-15 * max(scale, 1.0)
        )


def _norm_squared_in_boxes(
    field: dict[tuple[int, int, int], np.ndarray],
    centers: list[tuple[int, int, int]],
    half_width: int,
) -> float:
    return sum(
        float(np.vdot(value, value).real)
        for wave, value in field.items()
        if any(
            all(
                abs(component - center_component) <= half_width
                for component, center_component in zip(wave, center)
            )
            for center in centers
        )
    )


def _negate(wave: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(-value for value in wave)


def test_width_one_matches_every_sparse_coefficient_and_parseval_norm() -> None:
    config = MesoscopicCloudConfig(
        base_scale=8,
        gamma=1.0,
        width_override=1,
        phase_slope=(0.13, -0.07, 0.05),
    )
    _assert_all_coefficients_match(config)
    diagnostics = measure_local_fft_cloud(config)
    sparse = exact_sparse_leray_convolution(build_sparse_parent(config))
    sparse_squared = sum(float(np.vdot(value, value).real) for value in sparse.values())
    assert diagnostics.input_side_length == 1
    assert diagnostics.padded_side_length == 1
    assert diagnostics.zero_padding_is_alias_free
    assert diagnostics.output_block_count == 9
    assert diagnostics.unique_output_mode_count == len(sparse)
    assert diagnostics.full_nonlinear_norm_squared == pytest.approx(
        sparse_squared, rel=2.0e-14
    )
    assert diagnostics.child_forcing_fraction == pytest.approx(1.0, abs=2.0e-15)
    assert diagnostics.parent_divergence_relative < 2.0e-16
    assert diagnostics.parent_reality_relative == 0.0
    assert diagnostics.nonlinear_divergence_relative < 2.0e-16
    assert diagnostics.nonlinear_reality_relative < 2.0e-15
    assert "sum_k" in diagnostics.parseval_convention
    assert "K=2L-1" in diagnostics.convolution_convention


def test_width_two_matches_every_coefficient_norm_tag_and_sparse_backend() -> None:
    config = MesoscopicCloudConfig(
        base_scale=8,
        gamma=1.0,
        width_override=2,
        energy_constant=1.25,
        phase_slope=(0.13, -0.07, 0.05),
        exact_target_limit=1000,
    )
    _assert_all_coefficients_match(config)
    diagnostics = measure_local_fft_cloud(config)
    sparse_metrics = measure_mesoscopic_cloud(config)
    sparse = exact_sparse_leray_convolution(build_sparse_parent(config))
    full_squared = sum(float(np.vdot(value, value).real) for value in sparse.values())

    assert diagnostics.input_side_length == 3
    assert diagnostics.padded_side_length == 5
    assert diagnostics.full_nonlinear_norm_squared == pytest.approx(
        full_squared, rel=8.0e-14
    )
    assert diagnostics.full_nonlinear_norm == pytest.approx(
        sparse_metrics.full_nonlinear_norm, rel=8.0e-14
    )
    assert diagnostics.child_core_norm == pytest.approx(
        sparse_metrics.a_critical, rel=8.0e-14
    )
    assert diagnostics.child_forcing_fraction == pytest.approx(
        sparse_metrics.child_forcing_fraction, rel=8.0e-14
    )
    assert diagnostics.off_child_main_ratio == pytest.approx(
        sparse_metrics.off_chain_main_ratio, rel=8.0e-14
    )
    assert diagnostics.low_side_forcing_ratio == pytest.approx(
        sparse_metrics.low_side_forcing_ratio, rel=8.0e-14
    )
    assert diagnostics.difference_sideband_main_ratio == pytest.approx(
        sparse_metrics.difference_sideband_main_ratio, rel=8.0e-14
    )
    assert diagnostics.self_interaction_cross_interaction_ratio == pytest.approx(
        sparse_metrics.self_interaction_cross_interaction_ratio, rel=8.0e-14
    )
    assert diagnostics.gain_g == pytest.approx(sparse_metrics.gain_g, rel=8.0e-14)
    assert diagnostics.heat_factor_h == pytest.approx(
        sparse_metrics.heat_factor_h, rel=8.0e-14
    )
    assert diagnostics.duhamel_energy_ratio == pytest.approx(
        sparse_metrics.duhamel_energy_ratio, rel=8.0e-14
    )
    assert diagnostics.duhamel_identity_relative_error < 2.0e-14
    assert diagnostics.energy_cancellation_residual < 2.0e-14


def test_overlapping_output_boxes_are_globally_combined_before_norms() -> None:
    config = MesoscopicCloudConfig(
        base_scale=3,
        gamma=1.0,
        width_override=2,
        energy_constant=1.0,
        phase_slope=(0.19, -0.11, 0.03),
        exact_target_limit=1000,
    )
    _assert_all_coefficients_match(config)
    diagnostics = measure_local_fft_cloud(config)
    parent = build_sparse_parent(config)
    sparse = exact_sparse_leray_convolution(parent)
    full_squared = sum(float(np.vdot(value, value).real) for value in sparse.values())
    p = tuple(config.base_scale * value for value in config.p_direction)
    q = tuple(config.base_scale * value for value in config.q_direction)
    child = tuple(a + b for a, b in zip(p, q))
    difference = tuple(a - b for a, b in zip(p, q))
    full_half = 2 * (config.width - 1)
    core_half = config.width - 1
    child_core_squared = _norm_squared_in_boxes(
        sparse, [child, _negate(child)], core_half
    )
    child_full_squared = _norm_squared_in_boxes(
        sparse, [child, _negate(child)], full_half
    )
    difference_squared = _norm_squared_in_boxes(
        sparse, [difference, _negate(difference)], full_half
    )
    low_squared = _norm_squared_in_boxes(sparse, [(0, 0, 0)], full_half)

    assert diagnostics.output_boxes_overlap
    assert not diagnostics.channel_tags_are_orthogonal
    assert diagnostics.unique_output_mode_count == len(sparse)
    assert diagnostics.full_nonlinear_norm_squared == pytest.approx(
        full_squared, rel=1.0e-13
    )
    assert diagnostics.child_core_norm_squared == pytest.approx(
        child_core_squared, rel=1.0e-13
    )
    assert diagnostics.child_full_norm_squared == pytest.approx(
        child_full_squared, rel=1.0e-13
    )
    assert diagnostics.difference_full_norm_squared == pytest.approx(
        difference_squared, rel=1.0e-13
    )
    assert diagnostics.low_full_norm_squared == pytest.approx(
        low_squared, rel=1.0e-13
    )
    assert diagnostics.off_child_norm_squared == pytest.approx(
        full_squared - child_core_squared, rel=1.0e-13
    )
    pairing = sum(
        complex(np.vdot(value, sparse.get(wave, np.zeros(3))))
        for wave, value in parent.items()
    )
    assert diagnostics.energy_pairing_real == pytest.approx(pairing.real, abs=2.0e-15)
    assert diagnostics.energy_pairing_imag == pytest.approx(pairing.imag, abs=2.0e-15)
    assert diagnostics.energy_cancellation_residual < 2.0e-14
    assert diagnostics.nonlinear_divergence_relative < 2.0e-15
    assert diagnostics.nonlinear_reality_relative < 2.0e-15
    assert diagnostics.full_nonlinear_norm_squared == pytest.approx(
        diagnostics.self_interaction_norm_squared
        + diagnostics.cross_interaction_norm_squared
        + 2.0 * diagnostics.self_cross_real_pairing,
        rel=2.0e-14,
    )


def test_in_memory_measurement_is_deterministic() -> None:
    config = MesoscopicCloudConfig(
        base_scale=8,
        gamma=1.0,
        width_override=2,
    )
    first = measure_local_fft_cloud(config)
    second = measure_local_fft_cloud(config)
    for name in (
        "parent_norm_squared",
        "full_nonlinear_norm_squared",
        "child_core_norm_squared",
        "child_full_norm_squared",
        "difference_full_norm_squared",
        "low_full_norm_squared",
        "energy_pairing_real",
        "nonlinear_divergence_relative",
        "nonlinear_reality_relative",
        "duhamel_energy_ratio",
    ):
        assert getattr(first, name) == pytest.approx(
            getattr(second, name), rel=0.0, abs=0.0
        )
    assert first.estimated_scratch_bytes == 0
    assert first.estimated_peak_working_bytes > 0


def test_resource_and_materialization_limits_fail_before_large_allocation() -> None:
    config = MesoscopicCloudConfig(
        base_scale=64,
        gamma=1.0,
        width_override=27,
    )
    with pytest.raises(ValueError, match="working memory"):
        measure_local_fft_cloud(
            config,
            maximum_working_bytes=100_000_000,
        )
    with pytest.raises(ValueError, match="conservative output count"):
        local_fft_leray_coefficients(config, maximum_output_modes=1000)


@pytest.mark.parametrize("bad_limit", [0, -1, True, 1.5])
def test_memory_limit_is_a_positive_integer(bad_limit) -> None:
    config = MesoscopicCloudConfig(base_scale=8, gamma=1.0, width_override=1)
    with pytest.raises(ValueError, match="maximum_working_bytes"):
        measure_local_fft_cloud(config, maximum_working_bytes=bad_limit)


def test_parent_overlap_is_rejected_before_fft_convolution() -> None:
    config = MesoscopicCloudConfig(base_scale=2, gamma=1.0, width_override=2)
    with pytest.raises(ValueError, match="parent boxes overlap"):
        measure_local_fft_cloud(config)


def test_public_metrics_adapter_matches_exact_small_row_and_structure() -> None:
    config = MesoscopicCloudConfig(
        base_scale=8,
        gamma=1.0,
        width_override=2,
        exact_target_limit=1000,
    )
    adapted = measure_local_fft_mesoscopic_metrics(config)
    sparse = measure_mesoscopic_cloud(config)
    assert adapted.measurement_method == (
        "exact-zero-padded-local-fft-global-combination"
    )
    assert adapted.measurement_available
    assert adapted.targets_per_band == 5**3
    for name in (
        "parent_norm_squared",
        "effective_mode_count",
        "parent_min_wavenumber",
        "parent_max_wavenumber",
        "parent_rms_wavenumber",
        "a_critical",
        "a_unit",
        "gain_g",
        "heat_factor_h",
        "duhamel_norm_ratio",
        "duhamel_energy_ratio",
        "child_forcing_fraction",
        "off_chain_main_ratio",
        "low_side_forcing_ratio",
        "difference_sideband_main_ratio",
        "full_nonlinear_norm",
    ):
        assert getattr(adapted, name) == pytest.approx(
            getattr(sparse, name), rel=2.0e-13
        )
    assert adapted.effective_mode_count <= adapted.mode_count
    assert adapted.duhamel_norm_ratio**2 == pytest.approx(
        adapted.duhamel_energy_ratio, rel=2.0e-14
    )
    assert adapted.duhamel_energy_ratio <= adapted.effective_duhamel_upper_bound
    assert adapted.duhamel_upper_bound_method.startswith("phase-independent")
    assert adapted.energy_cancellation_available
    assert adapted.energy_cancellation_residual < 2.0e-14


def test_public_adapter_marks_overlapping_tags_but_keeps_exact_full_metrics() -> None:
    config = MesoscopicCloudConfig(
        base_scale=3,
        gamma=0.8,
        width_override=2,
        exact_target_limit=1000,
    )
    adapted = measure_local_fft_mesoscopic_metrics(config)
    exact = measure_local_fft_cloud(config)
    assert not adapted.full_output_bands_disjoint
    assert not adapted.channel_ratios_are_orthogonal
    assert adapted.forcing_metric_method.startswith("exact-local-fft-full-global")
    assert adapted.full_nonlinear_norm == pytest.approx(
        exact.full_nonlinear_norm, rel=2.0e-14
    )
    assert adapted.child_forcing_fraction == pytest.approx(
        exact.child_forcing_fraction, rel=2.0e-14
    )
    assert adapted.self_interaction_cross_interaction_ratio == pytest.approx(
        exact.self_interaction_cross_interaction_ratio, rel=2.0e-14
    )
    assert not adapted.asymptotic_screen_pass
    assert "generic_duhamel_ratio_decays_for_gamma_below_one" in (
        adapted.asymptotic_screen_reasons
    )
