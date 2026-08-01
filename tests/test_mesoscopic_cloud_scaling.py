import math

import numpy as np
import pytest

import ns_certificate_lab.mesoscopic_cloud_scaling as mesoscopic
from ns_certificate_lab.leray_response_relay import leray_advection
from ns_certificate_lab.mesoscopic_cloud_scaling import (
    MesoscopicCloudConfig,
    build_sparse_parent,
    exact_sparse_leray_convolution,
    measure_mesoscopic_cloud,
)


def _child_mask(config: MesoscopicCloudConfig, grid_size: int) -> np.ndarray:
    frequency = np.rint(
        np.fft.fftfreq(grid_size, d=1.0 / grid_size)
    ).astype(np.int64)
    kx, ky, kz = np.meshgrid(frequency, frequency, frequency, indexing="ij")
    child = config.base_scale * np.asarray(
        tuple(a + b for a, b in zip(config.p_direction, config.q_direction)),
        dtype=np.int64,
    )
    half = config.width - 1
    positive = (
        (np.abs(kx - child[0]) <= half)
        & (np.abs(ky - child[1]) <= half)
        & (np.abs(kz - child[2]) <= half)
    )
    negative = (
        (np.abs(kx + child[0]) <= half)
        & (np.abs(ky + child[1]) <= half)
        & (np.abs(kz + child[2]) <= half)
    )
    return positive | negative


def test_exact_small_measurement_has_the_required_structures_and_identity() -> None:
    config = MesoscopicCloudConfig(
        base_scale=8,
        gamma=0.4,
        energy_constant=1.25,
        exact_target_limit=1000,
    )
    metrics = measure_mesoscopic_cloud(config)
    assert metrics.width == 2
    assert metrics.measurement_method == "exact-target"
    assert metrics.parent_boxes_disjoint
    assert metrics.child_band_isolated
    assert metrics.channel_bands_disjoint
    assert metrics.full_output_bands_disjoint
    assert metrics.measurement_available
    assert metrics.parent_norm_squared == pytest.approx(2.5 / 8.0, rel=2.0e-14)
    assert metrics.normalization_relative_error < 2.0e-14
    assert metrics.divergence_relative < 2.0e-16
    assert metrics.reality_relative == 0.0
    assert metrics.a_critical > 0.0
    assert metrics.a_unit > 0.0
    assert metrics.chi_support > 0.0
    assert metrics.gain_g > 0.0
    assert metrics.duhamel_norm_ratio**2 == pytest.approx(
        metrics.duhamel_energy_ratio, rel=2.0e-14
    )
    assert metrics.duhamel_energy_ratio == pytest.approx(
        metrics.duhamel_identity_rhs, rel=2.0e-14
    )
    assert metrics.duhamel_identity_relative_error < 2.0e-14
    assert 0.0 < metrics.child_forcing_fraction < 1.0
    assert metrics.low_side_forcing_ratio >= 0.0
    assert metrics.off_chain_main_ratio >= 0.0
    assert metrics.difference_sideband_main_ratio >= 0.0
    assert metrics.self_interaction_cross_interaction_ratio >= 0.0
    assert metrics.energy_cancellation_available
    assert metrics.energy_cancellation_method.startswith("exact-sparse-all-pairs")
    assert metrics.energy_cancellation_residual < 2.0e-14
    assert abs(metrics.energy_cancellation_pairing) < 2.0e-14
    assert metrics.full_nonlinear_norm > 0.0
    assert sum(
        metrics.channel_ratios[name]
        for name in (
            "child",
            "child_spill",
            "difference",
            "difference_spill",
            "harmonic",
            "same_carrier_spill",
            "low",
            "low_spill",
        )
    ) == pytest.approx(1.0, rel=2.0e-14)


def test_exact_target_sum_matches_sparse_pair_enumeration_and_fft() -> None:
    config = MesoscopicCloudConfig(
        base_scale=8,
        gamma=0.4,
        energy_constant=1.0,
        phase_slope=(0.13, -0.07, 0.05),
        exact_target_limit=1000,
    )
    metrics = measure_mesoscopic_cloud(config)
    parent = build_sparse_parent(config)
    for wave, coefficient in parent.items():
        reflected = tuple(-component for component in wave)
        assert parent[reflected] == pytest.approx(
            np.conjugate(coefficient), abs=2.0e-15
        )
        assert np.dot(np.asarray(wave, dtype=float), coefficient) == pytest.approx(
            0.0, abs=2.0e-15
        )
    assert sum(float(np.vdot(value, value).real) for value in parent.values()) == (
        pytest.approx(2.0 * config.energy_constant / config.base_scale, rel=2.0e-14)
    )
    nonlinear = exact_sparse_leray_convolution(parent)
    child_center = tuple(
        config.base_scale * (a + b)
        for a, b in zip(config.p_direction, config.q_direction)
    )
    half = config.width - 1

    def in_child(wave: tuple[int, int, int]) -> bool:
        return any(
            all(abs(k - sign * center) <= half for k, center in zip(wave, child_center))
            for sign in (-1, 1)
        )

    sparse_a = math.sqrt(
        sum(
            float(np.vdot(value, value).real)
            for wave, value in nonlinear.items()
            if in_child(wave)
        )
    )
    assert metrics.a_critical == pytest.approx(sparse_a, rel=2.0e-13)
    full_squared = sum(
        float(np.vdot(value, value).real) for value in nonlinear.values()
    )
    assert metrics.child_forcing_fraction == pytest.approx(
        sparse_a**2 / full_squared, rel=3.0e-13
    )
    assert metrics.off_chain_main_ratio == pytest.approx(
        math.sqrt((full_squared - sparse_a**2) / sparse_a**2), rel=3.0e-13
    )
    pairing = sum(
        float(np.vdot(value, nonlinear.get(wave, np.zeros(3))).real)
        for wave, value in parent.items()
    )
    parent_squared = sum(
        float(np.vdot(value, value).real) for value in parent.values()
    )
    support_radius_factor = max(
        math.sqrt(sum(component * component for component in wave))
        for wave in parent
    ) / config.base_scale
    coefficient_l1 = sum(
        math.sqrt(float(np.vdot(value, value).real))
        for value in parent.values()
    )
    effective_modes = coefficient_l1**2 / parent_squared
    universal_bound = (
        2.0
        * support_radius_factor**2
        * config.tau**2
        * config.energy_constant
        * len(parent)
        / config.base_scale**3
    )
    effective_bound = universal_bound * effective_modes / len(parent)
    assert metrics.support_radius_factor_kappa == pytest.approx(
        support_radius_factor, rel=2.0e-14
    )
    assert metrics.universal_duhamel_upper_bound == pytest.approx(
        universal_bound, rel=3.0e-14
    )
    assert metrics.effective_duhamel_upper_bound == pytest.approx(
        effective_bound, rel=3.0e-14
    )
    assert metrics.duhamel_energy_ratio <= effective_bound
    assert metrics.duhamel_to_universal_upper_bound == pytest.approx(
        metrics.duhamel_energy_ratio / universal_bound, rel=3.0e-14
    )
    assert metrics.duhamel_to_effective_upper_bound == pytest.approx(
        metrics.duhamel_energy_ratio / effective_bound, rel=3.0e-14
    )
    cancellation = abs(pairing) / math.sqrt(parent_squared * full_squared)
    assert metrics.energy_cancellation_pairing == pytest.approx(
        pairing, abs=2.0e-15
    )
    assert metrics.energy_cancellation_residual == pytest.approx(
        cancellation, abs=2.0e-15
    )

    def norm_squared_near(
        centers: list[tuple[int, int, int]], selected_half: int
    ) -> float:
        return sum(
            float(np.vdot(value, value).real)
            for wave, value in nonlinear.items()
            if any(
                all(
                    abs(component - center) <= selected_half
                    for component, center in zip(wave, selected_center)
                )
                for selected_center in centers
            )
        )

    p_center = tuple(config.base_scale * value for value in config.p_direction)
    q_center = tuple(config.base_scale * value for value in config.q_direction)
    difference_center = tuple(
        config.base_scale * (a - b)
        for a, b in zip(config.p_direction, config.q_direction)
    )
    full_half = 2 * (config.width - 1)

    def negate(wave):
        return tuple(-value for value in wave)

    difference_squared = norm_squared_near(
        [difference_center, negate(difference_center)], full_half
    )
    low_squared = norm_squared_near([(0, 0, 0)], full_half)
    harmonic_centers = [
        tuple(2 * value for value in p_center),
        tuple(-2 * value for value in p_center),
        tuple(2 * value for value in q_center),
        tuple(-2 * value for value in q_center),
    ]
    harmonic_squared = norm_squared_near(harmonic_centers, full_half)
    child_full_squared = norm_squared_near(
        [child_center, negate(child_center)], full_half
    )
    assert metrics.low_side_forcing_ratio == pytest.approx(
        math.sqrt(low_squared / sparse_a**2), rel=3.0e-13
    )
    assert metrics.difference_sideband_main_ratio == pytest.approx(
        math.sqrt(difference_squared / sparse_a**2), rel=3.0e-13
    )
    assert metrics.self_interaction_cross_interaction_ratio == pytest.approx(
        math.sqrt(
            (harmonic_squared + low_squared)
            / (child_full_squared + difference_squared)
        ),
        rel=3.0e-13,
    )

    sparse_full_a_squared = sum(
        float(np.vdot(value, value).real)
        for wave, value in nonlinear.items()
        if any(
            all(
                abs(k - sign * center) <= full_half
                for k, center in zip(wave, child_center)
            )
            for sign in (-1, 1)
        )
    )
    assert metrics.channel_norms["child_spill"] ** 2 == pytest.approx(
        sparse_full_a_squared - sparse_a**2, rel=3.0e-13, abs=2.0e-14
    )

    grid_size = 48
    dense = np.zeros((3, grid_size, grid_size, grid_size), dtype=np.complex128)
    for wave, coefficient in parent.items():
        index = tuple(component % grid_size for component in wave)
        dense[(slice(None), *index)] = coefficient
    dense_nonlinear = leray_advection(dense, dense)
    mask = _child_mask(config, grid_size)
    fft_a = math.sqrt(float(np.vdot(
        dense_nonlinear * mask[None, ...],
        dense_nonlinear * mask[None, ...],
    ).real))
    assert metrics.a_critical == pytest.approx(fft_a, rel=3.0e-13)


def test_critical_rescaling_changes_a_but_not_unit_shape_or_gain() -> None:
    smaller = measure_mesoscopic_cloud(
        MesoscopicCloudConfig(
            base_scale=8,
            gamma=0.4,
            energy_constant=0.5,
            exact_target_limit=1000,
        )
    )
    larger = measure_mesoscopic_cloud(
        MesoscopicCloudConfig(
            base_scale=8,
            gamma=0.4,
            energy_constant=2.0,
            exact_target_limit=1000,
        )
    )
    assert larger.a_critical / smaller.a_critical == pytest.approx(4.0, rel=2.0e-13)
    assert larger.a_unit == pytest.approx(smaller.a_unit, rel=2.0e-13)
    assert larger.gain_g == pytest.approx(smaller.gain_g, rel=2.0e-13)
    assert larger.duhamel_energy_ratio / smaller.duhamel_energy_ratio == pytest.approx(
        4.0, rel=2.0e-13
    )


def test_g_growth_cannot_pass_the_subcritical_duhamel_screen() -> None:
    config = MesoscopicCloudConfig(
        base_scale=16,
        gamma=0.8,
        exact_target_limit=0,
        strata_per_axis=1,
    )
    metrics = measure_mesoscopic_cloud(config)
    assert metrics.predicted_g_grows
    assert metrics.predicted_g_exponent == pytest.approx(0.2)
    assert metrics.predicted_duhamel_exponent == pytest.approx(-0.6)
    assert metrics.predicted_duhamel_decays
    assert not metrics.asymptotic_screen_pass
    assert "generic_duhamel_ratio_decays_for_gamma_below_one" in (
        metrics.asymptotic_screen_reasons
    )
    assert metrics.box_overlap_counts
    assert not metrics.measurement_available
    assert metrics.support_radius_factor_kappa is None
    assert metrics.universal_duhamel_upper_bound is None
    assert metrics.effective_duhamel_upper_bound is None
    assert metrics.duhamel_upper_bound_method == "unavailable-parent-overlap"


@pytest.mark.parametrize("rho", [0.10, 0.15, 0.20])
def test_fixed_relative_width_override_is_the_gamma_one_control_family(rho) -> None:
    width = math.floor(rho * 32)
    config = MesoscopicCloudConfig(
        base_scale=32,
        gamma=1.0,
        width_override=width,
        exact_target_limit=100,
        strata_per_axis=2,
    )
    metrics = measure_mesoscopic_cloud(config)
    assert metrics.width == width
    assert metrics.width_source == "override"
    assert metrics.relative_width == pytest.approx(width / 32.0)
    assert metrics.predicted_duhamel_exponent == 0.0
    assert not metrics.predicted_duhamel_decays
    assert metrics.parent_boxes_disjoint
    assert metrics.child_band_isolated
    assert metrics.channel_bands_disjoint
    assert metrics.full_output_bands_disjoint
    assert metrics.measurement_available


def test_stratified_measurement_is_deterministic_and_weighted() -> None:
    config = MesoscopicCloudConfig(
        base_scale=32,
        gamma=1.0,
        width_override=5,
        exact_target_limit=100,
        strata_per_axis=2,
    )
    first = measure_mesoscopic_cloud(config)
    second = measure_mesoscopic_cloud(config)
    assert first.measurement_method == "deterministic-stratified-target"
    assert first.targets_per_band == 8
    assert first.a_critical == second.a_critical
    assert first.duhamel_energy_ratio == second.duhamel_energy_ratio
    assert first.channel_ratios == second.channel_ratios


def test_child_core_uses_d_but_relay_clearance_uses_full_two_d_sumset() -> None:
    config = MesoscopicCloudConfig(
        base_scale=8,
        gamma=1.0,
        width_override=4,
        exact_target_limit=0,
        strata_per_axis=1,
    )
    metrics = measure_mesoscopic_cloud(config)
    # Core parent/child boxes have half-width d=3 and are geometrically
    # disjoint since N>2d.  The reachable child source has half-width 2d and
    # still touches the old parent because N<=3d, so the relay screen rejects.
    assert metrics.half_width == 3
    assert metrics.parent_boxes_disjoint
    assert metrics.channel_bands_disjoint
    assert not metrics.child_band_isolated
    assert not metrics.full_output_bands_disjoint
    assert not metrics.channel_ratios_are_orthogonal
    assert any(key.startswith("full_child") for key in metrics.box_overlap_counts)
    assert "child_band_not_isolated" in metrics.finite_screen_reasons
    assert metrics.child_forcing_fraction is None
    assert metrics.off_chain_main_ratio is None
    assert metrics.self_interaction_cross_interaction_ratio is None
    assert metrics.low_side_forcing_ratio is not None
    assert metrics.difference_sideband_main_ratio is not None
    assert metrics.forcing_metric_method.startswith(
        "overlapping-band-projections-only"
    )


def test_large_full_energy_diagnostic_is_explicitly_unavailable_not_zero() -> None:
    metrics = measure_mesoscopic_cloud(
        MesoscopicCloudConfig(
            base_scale=32,
            gamma=1.0,
            width_override=5,
            exact_target_limit=100,
            strata_per_axis=2,
            exact_energy_pair_limit=0,
        )
    )
    assert not metrics.energy_cancellation_available
    assert metrics.energy_cancellation_method.startswith("unavailable-pair-limit")
    assert metrics.energy_cancellation_residual is None
    assert metrics.energy_cancellation_pairing is None
    assert metrics.full_nonlinear_norm is None
    assert metrics.child_forcing_fraction is not None
    assert 0.0 <= metrics.child_forcing_fraction <= 1.0


def test_exact_full_metric_survives_overlap_without_calling_tags_orthogonal() -> None:
    metrics = measure_mesoscopic_cloud(
        MesoscopicCloudConfig(
            base_scale=3,
            gamma=1.0,
            width_override=2,
            exact_target_limit=1000,
        )
    )
    assert metrics.parent_boxes_disjoint
    assert not metrics.child_band_isolated
    assert not metrics.full_output_bands_disjoint
    assert not metrics.channel_ratios_are_orthogonal
    assert metrics.energy_cancellation_available
    assert metrics.energy_cancellation_residual < 2.0e-14
    assert metrics.forcing_metric_method == "exact-sparse-full"
    assert 0.0 < metrics.child_forcing_fraction < 1.0
    assert metrics.off_chain_main_ratio > 0.0
    # This label would require disjoint harmonic/low and child/difference
    # bands.  Exact full energy does not make an overlapping tag orthogonal.
    assert metrics.self_interaction_cross_interaction_ratio is None


@pytest.mark.parametrize(
    ("phase_p", "phase_q", "phase_slope"),
    [
        (0.0, 0.0, (0.0, 0.0, 0.0)),
        (1.7, -0.4, (0.31, -0.23, 0.19)),
        (-2.2, 2.8, (-0.47, 0.29, 0.11)),
    ],
)
def test_phase_independent_bounds_hold_for_arbitrary_affine_phases(
    phase_p, phase_q, phase_slope
) -> None:
    metrics = measure_mesoscopic_cloud(
        MesoscopicCloudConfig(
            base_scale=8,
            gamma=0.4,
            phase_p=phase_p,
            phase_q=phase_q,
            phase_slope=phase_slope,
            exact_target_limit=1000,
        )
    )
    assert metrics.measurement_method == "exact-target"
    assert metrics.duhamel_energy_ratio <= metrics.effective_duhamel_upper_bound
    assert (
        metrics.effective_duhamel_upper_bound
        <= metrics.universal_duhamel_upper_bound
    )
    assert 0.0 <= metrics.duhamel_to_effective_upper_bound <= 1.0


def test_measurement_above_phase_independent_bound_is_rejected(monkeypatch) -> None:
    original = mesoscopic._measure_band

    def inflated_heat(*args, **kwargs):
        norm_squared, heat_squared = original(*args, **kwargs)
        return norm_squared, heat_squared * 1.0e6

    monkeypatch.setattr(mesoscopic, "_measure_band", inflated_heat)
    with pytest.raises(AssertionError, match="exceeds the .* upper bound"):
        measure_mesoscopic_cloud(
            MesoscopicCloudConfig(
                base_scale=8,
                gamma=0.4,
                exact_target_limit=1000,
            )
        )


@pytest.mark.parametrize("bad_width", [0, -1, True])
def test_width_override_is_validated(bad_width) -> None:
    with pytest.raises(ValueError, match="width_override"):
        MesoscopicCloudConfig(
            base_scale=8,
            gamma=1.0,
            width_override=bad_width,
        )
