"""Tests for the exact spectral front-gap ledger."""

from __future__ import annotations

from fractions import Fraction

import pytest

from ns_certificate_lab.exact_leray_relay import build_exact_relay_triad
from ns_certificate_lab.fourier_torus import TrigVector
from ns_certificate_lab.mesoscopic_cloud_scaling import (
    MesoscopicCloudConfig,
    build_sparse_parent,
)
from ns_certificate_lab.spectral_front_monotone import (
    front_gap_identity,
    front_wavenumber,
    inverse_square_lattice_sum,
    lemma_k_certificate,
    sparse_parent_delocalization,
)


def _four_parent_field() -> TrigVector:
    third = Fraction(1, 3)
    return TrigVector.from_modes(
        [
            ((1, 1, 0), (0, 0, 0), (0, 0, 1)),
            ((1, 0, 1), (0, 1, 0), (0, 0, 0)),
            ((0, 1, -1), (third, 2 * third, 2 * third), (0, 0, 0)),
            ((1, 0, -1), (2 * third, third, 2 * third), (0, 0, 0)),
        ]
    )


def test_front_gap_identity_on_relay_triad_all_orders() -> None:
    field = build_exact_relay_triad()
    for order in (0, 1, 2):
        record = front_gap_identity(
            field, order=order, viscosity=Fraction(1, 40)
        )
        assert record.gap_total >= 0
        assert 0 <= record.saturation_deficit <= 1
        assert record.closable_upper >= record.log_bandwidth_derivative
        assert record.h_r > 0 and record.h_r1 > 0
        payload = record.as_dict()
        assert payload["order"] == order


def test_front_gap_identity_matches_exact_variance() -> None:
    field = _four_parent_field()
    record = front_gap_identity(field, order=0, viscosity=Fraction(1, 10))
    reconstructed = record.h_r2 / record.h_r - record.bandwidth_squared**2
    assert record.variance_direct == reconstructed
    assert record.gap_total == (
        record.closable_upper - record.log_bandwidth_derivative
    )


def test_single_mode_field_has_zero_nonlinear_ledger() -> None:
    field = TrigVector.from_modes([((1, 1, 0), (0, 0, 0), (0, 0, 1))])
    record = front_gap_identity(field, order=0, viscosity=1)
    assert record.g_r == 0
    assert record.closable_upper == 0
    assert record.log_bandwidth_derivative == 0
    assert record.saturation_deficit == 0


def test_two_shell_linear_field_narrows_bandwidth() -> None:
    field = TrigVector.from_modes(
        [
            ((1, 0, 0), (0, 0, 0), (0, 1, 0)),
            ((0, 3, 0), (1, 0, 0), (0, 0, 0)),
        ]
    )
    record = front_gap_identity(field, order=0, viscosity=Fraction(1, 4))
    if record.g_r == 0:
        assert record.log_bandwidth_derivative < 0
        assert record.gap_total == -record.log_bandwidth_derivative


def test_inverse_square_lattice_sum_small_values() -> None:
    assert inverse_square_lattice_sum(1) == 6
    assert inverse_square_lattice_sum(2) == 12
    assert inverse_square_lattice_sum(3) == Fraction(44, 3)


def test_lemma_k_certificate_on_relay_triad() -> None:
    field = build_exact_relay_triad()
    record = lemma_k_certificate(field)
    assert Fraction(record["margin"]) >= 0
    assert record["band_norm_squared"] == 6
    assert Fraction(record["front_wavenumber_K"]) == front_wavenumber(field)


def test_front_wavenumber_positive_on_triad() -> None:
    assert front_wavenumber(build_exact_relay_triad()) > 0


def test_sparse_parent_delocalization_equal_amplitudes() -> None:
    parent = {
        (4, 4, 0): (0.0, 0.0, 1.0),
        (4, 0, 4): (0.0, 1.0, 0.0),
    }
    record = sparse_parent_delocalization(parent, scale=4)
    assert record["support_mode_count"] == 4.0
    assert record["effective_mode_count"] == pytest.approx(4.0)
    assert record["mu_delocalization"] == pytest.approx(4.0 / 64.0)


def test_sparse_parent_delocalization_on_mesoscopic_family() -> None:
    config = MesoscopicCloudConfig(
        base_scale=16, gamma=1.0, width_override=3
    )
    parent = build_sparse_parent(config)
    record = sparse_parent_delocalization(parent, scale=16)
    assert 0.0 < record["effective_over_support"] <= 1.0
    assert record["mu_delocalization"] > 0.0


def test_front_gap_identity_rejects_bad_inputs() -> None:
    field = build_exact_relay_triad()
    with pytest.raises(ValueError):
        front_gap_identity(field, order=-1, viscosity=1)
    with pytest.raises(ValueError):
        front_gap_identity(field, order=0, viscosity=0)
    with pytest.raises(ValueError):
        inverse_square_lattice_sum(0)
    with pytest.raises(ValueError):
        sparse_parent_delocalization({}, scale=4)
