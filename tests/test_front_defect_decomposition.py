"""Exact defect decomposition of the closable front bound (Lambda O-9)."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from ns_certificate_lab.exact_leray_relay import build_exact_relay_triad
from ns_certificate_lab.fourier_torus import family_P1, family_P2, family_P3
from ns_certificate_lab.spectral_front_monotone import (
    front_defect_decomposition,
    front_gap_identity,
    full_nonlinear_power,
)

FIELDS = {
    "relay_s1": build_exact_relay_triad(scale=1),
    "relay_s2": build_exact_relay_triad(scale=2),
    "relay_small_child": build_exact_relay_triad(
        scale=1, child_cosine=Fraction(1, 32)
    ),
    "P1": family_P1(),
    "P2": family_P2(),
    "P3": family_P3(),
}


@pytest.mark.parametrize("name", sorted(FIELDS))
@pytest.mark.parametrize("order", [0, 1, 2])
@pytest.mark.parametrize("viscosity", [Fraction(1, 40), Fraction(1, 10)])
@pytest.mark.parametrize("convention", ["full", "in_support"])
def test_defects_are_nonnegative_and_telescope(
    name, order, viscosity, convention
):
    record = front_defect_decomposition(
        FIELDS[name], order=order, viscosity=viscosity, convention=convention
    )
    # Exactly rational half of the telescoping.
    assert record.delta_sign >= 0
    assert record.rational_remainder >= 0
    assert record.gap_total == record.delta_sign + record.rational_remainder
    # Enclosed half: the three radical defects sum to the rational remainder.
    assert record.split_residual[0] <= 0 <= record.split_residual[1]
    for enclosure in (
        record.delta_cs_modal,
        record.delta_cs_vector,
        record.delta_square_completion,
    ):
        assert enclosure[0] <= enclosure[1]
        assert enclosure[1] >= 0
    assert all(value >= 0 for value in record.per_mode_cs_defect.values())


@pytest.mark.parametrize("name", sorted(FIELDS))
def test_in_support_convention_reproduces_the_published_ledger(name):
    record = front_defect_decomposition(
        FIELDS[name], order=0, viscosity=Fraction(1, 40), convention="in_support"
    )
    legacy = front_gap_identity(FIELDS[name], order=0, viscosity=Fraction(1, 40))
    assert record.gap_total == legacy.gap_total
    assert record.saturation_deficit == legacy.saturation_deficit
    assert record.g_r == legacy.g_r


@pytest.mark.parametrize("name", sorted(FIELDS))
def test_full_nonlinear_power_dominates_the_in_support_moment(name):
    record = front_defect_decomposition(
        FIELDS[name], order=0, viscosity=Fraction(1), convention="full"
    )
    assert record.g_full >= record.g_in
    assert record.g_full == full_nonlinear_power(FIELDS[name], 0)[0]


@pytest.mark.parametrize("name", sorted(FIELDS))
@pytest.mark.parametrize("viscosity", [Fraction(1, 40), Fraction(1, 3)])
def test_optimal_deficit_is_the_minimum_over_viscosity(name, viscosity):
    record = front_defect_decomposition(
        FIELDS[name], order=0, viscosity=viscosity, convention="full"
    )
    assert record.optimal_deficit <= record.saturation_deficit
    if record.optimal_viscosity is not None:
        at_optimum = front_defect_decomposition(
            FIELDS[name],
            order=0,
            viscosity=record.optimal_viscosity,
            convention="full",
        )
        assert at_optimum.saturation_deficit == record.optimal_deficit


def test_relay_family_closed_form_matches_the_exact_ledger():
    """d*_full = 15D^2/(2+18D^2), d*_in = 6D^2/(2+9D^2) at B=C=1."""

    for scale in (1, 2, 3):
        for child in (Fraction(1), Fraction(1, 8), Fraction(1, 128)):
            field = build_exact_relay_triad(scale=scale, child_cosine=child)
            t = child * child
            for convention, expected in (
                ("full", 15 * t / (2 + 18 * t)),
                ("in_support", 6 * t / (2 + 9 * t)),
            ):
                record = front_defect_decomposition(
                    field,
                    order=0,
                    viscosity=Fraction(1),
                    convention=convention,
                )
                assert record.optimal_deficit == expected


def test_relay_family_deficit_has_no_positive_floor():
    """The infimum of the r=0 deficit over divergence-free fields is 0."""

    previous = None
    for child in (Fraction(1, 8), Fraction(1, 32), Fraction(1, 128)):
        record = front_defect_decomposition(
            build_exact_relay_triad(scale=1, child_cosine=child),
            order=0,
            viscosity=Fraction(1),
            convention="full",
        )
        assert record.optimal_deficit > 0
        if previous is not None:
            assert record.optimal_deficit < previous
        previous = record.optimal_deficit
    assert previous < Fraction(1, 1000)


def test_energy_neutrality_is_automatic_at_r0_saturation():
    """a_k = gamma (x_k - mu) e_k forces T_0 = gamma (H_1 - mu H_0) = 0."""

    for field in FIELDS.values():
        record = front_defect_decomposition(
            field, order=0, viscosity=Fraction(1), convention="full"
        )
        assert record.t_r == 0
        assert record.h_r1 == record.bandwidth_squared * record.h_r


def test_analytic_gradient_matches_central_differences():
    """Guard for the float search lane of run_lambda_o9_defect_search."""

    import importlib.util
    import sys
    from pathlib import Path

    path = (
        Path(__file__).resolve().parents[1]
        / "experiments"
        / "run_lambda_o9_defect_search.py"
    )
    spec = importlib.util.spec_from_file_location("_lambda_o9", path)
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolve annotations through sys.modules; register first.
    sys.modules["_lambda_o9"] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop("_lambda_o9", None)

    band = module.Band(2)
    rng = np.random.default_rng(7)
    theta = rng.normal(size=band.parameter_count)
    for full in (True, False):
        value, gradient, _ = module.ratio_and_gradient(band, theta, full=full)
        assert value > 0.0
        step = 1e-6
        for index in rng.choice(band.parameter_count, 5, replace=False):
            plus, minus = theta.copy(), theta.copy()
            plus[index] += step
            minus[index] -= step
            difference = (
                module.ratio_and_gradient(band, plus, full=full)[0]
                - module.ratio_and_gradient(band, minus, full=full)[0]
            ) / (2 * step)
            assert gradient[index] == pytest.approx(difference, rel=1e-4, abs=1e-9)
