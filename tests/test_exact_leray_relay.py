from fractions import Fraction

import pytest

from ns_certificate_lab.exact_leray_relay import (
    build_exact_relay_triad,
    exact_relay_metrics,
    fixed_cardinality_scaling,
)
from ns_certificate_lab.fourier_torus import family_P2
from ns_certificate_lab.modal_front_actions import (
    h3_bandwidth_factorization,
    modal_growth_identity,
)


def test_exact_three_mode_relay_has_positive_post_viscous_child_flux() -> None:
    metrics = exact_relay_metrics(viscosity=Fraction(1, 40))
    assert metrics.parent_energy == Fraction(1, 2)
    assert metrics.child_energy == Fraction(3, 256)
    assert metrics.child_flux == Fraction(1, 16)
    assert metrics.parent_flux == -metrics.child_flux
    assert metrics.child_viscous_loss == Fraction(9, 2560)
    assert metrics.child_net == Fraction(151, 2560)
    assert metrics.total_nonlinear_energy_defect == 0
    assert not metrics.difference_mode_present
    assert metrics.child_coefficient_matches_formula


def test_exact_relay_closed_form_survives_integer_scaling() -> None:
    metrics = exact_relay_metrics(
        viscosity=Fraction(2, 7),
        scale=3,
        parent_sine=2,
        parent_cosine=Fraction(3, 2),
        child_cosine=Fraction(1, 11),
    )
    assert metrics.child_flux == Fraction(9, 22)
    assert metrics.off_chain_nonlinear_l2_squared == Fraction(675, 3872)
    assert metrics.off_chain_nonlinear_l2_squared == metrics.formula_off_chain_l2_squared
    assert metrics.total_nonlinear_energy_defect == 0


def test_fixed_cardinality_relay_is_explicitly_rejected_at_high_scale() -> None:
    scaling = fixed_cardinality_scaling()
    assert scaling["flux_to_viscosity_ratio"] == "N^-3/2"
    assert scaling["verdict"].startswith("REJECTED")


def test_modal_covariance_and_h3_factorization_are_exact() -> None:
    field = build_exact_relay_triad()
    for order in range(3):
        identity = modal_growth_identity(
            field, order=order, viscosity=Fraction(1, 40)
        )
        assert identity.covariance_front_term == identity.direct_front_term
        assert identity.cauchy_schwarz_gap >= 0
        assert identity.square_completion_gap >= 0
    factorization = h3_bandwidth_factorization(field)
    assert factorization["H3"] == factorization["reconstructed_H3"]


def test_modal_growth_sign_and_nontrivial_cauchy_gap_are_anchored() -> None:
    """P2 has three radii, so this is not the two-point equality case."""

    identity = modal_growth_identity(
        family_P2(), order=1, viscosity=Fraction(1, 40)
    )
    assert identity.t_r == 0
    assert identity.t_r1 == Fraction(-1, 72)
    assert identity.cauchy_schwarz_gap == Fraction(5617, 3086550)
    assert identity.log_bandwidth_derivative == Fraction(-29, 1900)
    assert identity.square_completion_gap > 0


@pytest.mark.parametrize(
    "call",
    [
        lambda: build_exact_relay_triad(scale=0),
        lambda: exact_relay_metrics(viscosity=0),
        lambda: modal_growth_identity(
            build_exact_relay_triad(), order=-1, viscosity=1
        ),
    ],
)
def test_exact_relay_invalid_inputs_are_rejected(call) -> None:
    with pytest.raises(ValueError):
        call()
