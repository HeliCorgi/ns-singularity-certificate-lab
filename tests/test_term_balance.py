from __future__ import annotations

import numpy as np
import pytest

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.pde import EquationTerms
from ns_certificate_lab.term_balance import (
    balance_gate,
    critical_core_mask,
    derivative_sensitivity,
    summarize_terms,
    time_derivatives,
)


def _grid() -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=17, nz=32, r_max=1.0, z_min=0.0, z_max=1.0, periodic_z=True
    )


def test_nonuniform_quadratic_time_derivative_is_exact() -> None:
    times = [0.0, 0.2, 0.7, 1.4]
    base = np.arange(12, dtype=np.float64).reshape(3, 4) + 1.0
    fields = [(time * time + 3.0 * time + 2.0) * base for time in times]
    derivatives = time_derivatives(times, fields)
    for time, derivative in zip(times, derivatives, strict=True):
        np.testing.assert_allclose(derivative, (2.0 * time + 3.0) * base, rtol=2e-14, atol=2e-14)


def test_derivative_sensitivity_vanishes_for_linear_series() -> None:
    times = [0.0, 0.3, 1.0]
    base = np.ones((4, 5), dtype=np.float64)
    fields = [(2.0 * time - 1.0) * base for time in times]
    sensitivity = derivative_sensitivity(times, fields)
    for item in sensitivity:
        np.testing.assert_allclose(item, 0.0, atol=2e-14)


def test_summary_detects_exact_time_inviscid_balance() -> None:
    grid = _grid()
    one = np.ones(grid.shape)
    zero = np.zeros(grid.shape)
    terms = EquationTerms(
        time=-one,
        advection_r=one,
        advection_z=zero,
        source=zero,
        diffusion=zero,
        residual=zero,
    )
    report = summarize_terms(grid, terms, time_sensitivity=zero)
    assert report.normalized_residual == pytest.approx(0.0)
    assert report.dominant_balance == "time_inviscid"


def test_summary_detects_source_advection_cancellation() -> None:
    grid = _grid()
    one = np.ones(grid.shape)
    zero = np.zeros(grid.shape)
    terms = EquationTerms(
        time=zero,
        advection_r=one,
        advection_z=zero,
        source=one,
        diffusion=zero,
        residual=zero,
    )
    report = summarize_terms(grid, terms, time_sensitivity=zero)
    assert report.advection_source_cancellation == pytest.approx(0.0)


def test_large_residual_is_unresolved() -> None:
    grid = _grid()
    one = np.ones(grid.shape)
    zero = np.zeros(grid.shape)
    terms = EquationTerms(
        time=zero,
        advection_r=one,
        advection_z=zero,
        source=zero,
        diffusion=zero,
        residual=one,
    )
    report = summarize_terms(grid, terms, time_sensitivity=zero)
    assert report.dominant_balance == "unresolved_or_other"


def test_periodic_core_mask_wraps_across_seam() -> None:
    grid = _grid()
    mask = critical_core_mask(
        grid,
        center_z=0.99,
        radial_scale=0.2,
        axial_scale=0.05,
        radius_multiplier=1.0,
    )
    assert mask[0, 0]
    assert mask[0, -1]


def test_bad_times_rejected() -> None:
    fields = [np.zeros((2, 2)), np.ones((2, 2))]
    with pytest.raises(ValueError, match="strictly increasing"):
        time_derivatives([0.0, 0.0], fields)


def test_balance_gate_rejects_sparse_derivative() -> None:
    grid = _grid()
    one = np.ones(grid.shape)
    zero = np.zeros(grid.shape)
    terms = EquationTerms(
        time=-one, advection_r=one, advection_z=zero, source=zero,
        diffusion=zero, residual=zero,
    )
    report = summarize_terms(grid, terms, time_sensitivity=0.5 * one)
    decision = balance_gate(report)
    assert not decision.passed
    assert decision.reasons == ("time_derivative_sensitivity",)
