from __future__ import annotations

from fractions import Fraction

import pytest

from ns_certificate_lab.exact_carrier_search import (
    DEFAULT_WAVEVECTORS,
    ModeSpec,
    search_exact_carrier_gadget,
)


def _coefficient(field, wave):
    cosine, sine = field.coefficient_table()[wave]
    return tuple(cosine), tuple(sine)


def test_default_search_is_deterministic_and_honestly_classified() -> None:
    first = search_exact_carrier_gadget()
    second = search_exact_carrier_gadget()
    assert first.partial_candidate_found
    assert not first.partial_search_exhaustive
    assert first.partial_gadget is not None
    assert first.as_dict() == second.as_dict()
    assert first.wavevectors == DEFAULT_WAVEVECTORS
    assert len(first.wavevectors) == 9
    assert first.mode_spec_count == 24
    assert first.unordered_mode_pair_tests == 276
    assert first.viable_relay_pair_count == 152
    assert first.partial_eligible_relay_pair_tests == 1
    assert first.classification == "partial_rejected_cross_talk"
    assert "REJECTED CROSS-TALK" in first.status
    assert "SIMPLE RECURSION REJECTED" in first.status
    payload = first.as_dict()
    assert "found" not in payload
    assert payload["partial_candidate_found"] is True
    assert payload["strict_search"]["found"] is False
    assert "nonexhaustive" in first.reason
    assert "222/2483" in first.reason
    assert "interact exactly to zero" in first.reason
    assert "REJECTED" in first.reason


def test_first_relay_recovers_known_difference_cancelling_leray_triad() -> None:
    certificate = search_exact_carrier_gadget()
    assert certificate.partial_gadget is not None
    relay = certificate.partial_gadget.relay_one
    assert relay.left_spec == ModeSpec((1, 1, 0), (0, 0, 1), "sin")
    assert relay.right_spec == ModeSpec((1, 0, 1), (0, 1, 0), "cos")
    interaction = relay.interaction
    assert interaction.sum_wave == (2, 1, 1)
    assert interaction.difference_wave == (0, 1, -1)
    assert not interaction.difference_output.coefficient_table()
    assert not interaction.other_output.coefficient_table()
    assert _coefficient(interaction.sum_output, (2, 1, 1)) == (
        (Fraction(-1, 3), Fraction(1, 3), Fraction(1, 3)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert _coefficient(interaction.child, (2, 1, 1)) == (
        (Fraction(1, 3), Fraction(-1, 3), Fraction(-1, 3)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert interaction.signed_parent_to_child_flux == Fraction(1, 6)
    assert interaction.parent_flux_in_populated_triad == Fraction(-1, 6)
    assert interaction.child_flux_in_populated_triad == Fraction(1, 6)
    assert interaction.populated_triad_energy_defect == 0


def test_second_relay_and_child_child_outputs_are_exact() -> None:
    certificate = search_exact_carrier_gadget()
    assert certificate.partial_gadget is not None
    second = certificate.partial_gadget.relay_two
    assert certificate.partial_gadget.common_parent_shell_squared == 2
    assert certificate.partial_gadget.common_child_shell_squared == 6
    assert second.left_spec == ModeSpec((0, 1, -1), (1, 0, 0), "cos")
    assert second.right_spec == ModeSpec((1, 0, -1), (0, 1, 0), "cos")
    interaction = second.interaction
    assert interaction.sum_wave == (1, 1, -2)
    assert interaction.difference_wave == (1, -1, 0)
    assert not interaction.difference_output.coefficient_table()
    assert not interaction.other_output.coefficient_table()
    assert _coefficient(interaction.child, (1, 1, -2)) == (
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)),
    )
    assert interaction.signed_parent_to_child_flux == Fraction(1, 6)
    assert interaction.populated_triad_energy_defect == 0

    next_stage = certificate.partial_gadget.next_interaction
    assert next_stage.sum_wave == (3, 2, -1)
    assert next_stage.difference_wave == (1, 0, 3)
    assert _coefficient(next_stage.sum_output, (3, 2, -1)) == (
        (Fraction(1, 7), Fraction(-5, 21), Fraction(-1, 21)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert _coefficient(next_stage.difference_output, (1, 0, 3)) == (
        (Fraction(-1, 5), Fraction(1, 3), Fraction(1, 15)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert next_stage.signed_parent_to_child_flux == Fraction(5, 126)
    assert Fraction(next_stage.difference_output.sobolev_sq(0)) == Fraction(7, 90)
    assert not next_stage.other_output.coefficient_table()
    assert next_stage.populated_triad_energy_defect == 0


def test_diagonal_cross_talk_rejects_the_partial_candidate_exactly() -> None:
    certificate = search_exact_carrier_gadget()
    assert certificate.partial_gadget is not None
    diagnostic = certificate.partial_gadget.diagonal_cross_talk
    assert diagnostic.target_shell_squared == 6
    assert not diagnostic.passes
    assert _coefficient(
        diagnostic.a1_b2_target_shell_output,
        (2, 1, -1),
    ) == (
        (Fraction(1, 3), Fraction(-1, 3), Fraction(1, 3)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert _coefficient(
        diagnostic.a2_b1_target_shell_output,
        (1, -1, 2),
    ) == (
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(1, 3), Fraction(-1, 3), Fraction(-1, 3)),
    )


def test_strict_known_relay_extension_search_is_an_exhaustive_negative() -> None:
    strict = search_exact_carrier_gadget().strict_search
    assert not strict.found
    assert strict.exhaustive
    assert strict.fixed_known_relay_available
    assert strict.eligible_second_relay_orientations_tested == 16
    assert len(strict.orientation_records) == 16
    assert [record.orientation_index for record in strict.orientation_records] == list(
        range(1, 17)
    )
    assert all(not record.strict_passes for record in strict.orientation_records)
    assert all(record.rejection_reasons for record in strict.orientation_records)
    assert all(
        record.diagonal_cross_talk.as_dict()["a1_b2_output"] is not None
        and record.diagonal_cross_talk.as_dict()["a2_b1_output"] is not None
        for record in strict.orientation_records
    )
    assert sum(
        record.required_next_interaction_passes
        for record in strict.orientation_records
    ) == 8
    assert sum(
        record.diagonal_cross_talk.passes for record in strict.orientation_records
    ) == 4
    assert strict.gadget is None
    assert "NEGATIVE STRICT CROSS-TALK EXTENSION" in strict.status
    assert "fixed known relay" in strict.reason
    assert "not a no-go result outside" in strict.reason


def test_leakage_fluxes_recursion_stop_and_energy_are_exact() -> None:
    certificate = search_exact_carrier_gadget()
    assert certificate.partial_gadget is not None
    gadget = certificate.partial_gadget
    cross_support = set(gadget.relay_cross_output.coefficient_table())
    assert (3, 2, -1) in cross_support
    assert (1, 0, 3) in cross_support
    assert gadget.relay_cross_output_l2_squared == Fraction(2483, 1890)
    assert gadget.intended_next_output_l2_squared == Fraction(37, 315)
    assert gadget.intended_fraction_of_relay_cross_power == Fraction(222, 2483)
    assert gadget.simple_binary_recursion_rejected
    assert gadget.grandchild_mutual_output_l2_squared == 0
    assert not gadget.grandchild_mutual_output.coefficient_table()
    assert dict(gadget.full_field_flux_by_wave) == {
        (0, 1, -1): Fraction(-1, 12),
        (1, 0, -1): Fraction(-1, 12),
        (1, 0, 1): Fraction(-1, 12),
        (1, 1, -2): Fraction(11, 63),
        (1, 1, 0): Fraction(-1, 12),
        (2, 1, 1): Fraction(5, 42),
        (3, 2, -1): Fraction(5, 126),
    }
    assert sum(dict(gadget.full_field_flux_by_wave).values(), Fraction(0)) == 0
    assert gadget.relay_one.interaction.populated_triad_energy_defect == 0
    assert gadget.relay_two.interaction.populated_triad_energy_defect == 0
    assert gadget.next_interaction.populated_triad_energy_defect == 0
    assert gadget.full_field_energy_defect == 0
    for field in (
        gadget.relay_one.interaction.child,
        gadget.relay_two.interaction.child,
        gadget.next_interaction.child,
        gadget.relay_cross_output,
        gadget.grandchild_mutual_output,
    ):
        assert not field.divergence().cleaned().terms


def test_empty_gadget_scope_returns_an_honest_negative_certificate() -> None:
    certificate = search_exact_carrier_gadget(((1, 0, 0),))
    assert not certificate.partial_candidate_found
    assert certificate.partial_search_exhaustive
    assert certificate.partial_gadget is None
    assert certificate.viable_relay_pair_count == 0
    assert certificate.partial_eligible_relay_pair_tests == 0
    assert not certificate.strict_search.found
    assert certificate.strict_search.exhaustive
    assert not certificate.strict_search.fixed_known_relay_available
    assert "STATED FINITE SEARCH ONLY" in certificate.status
    assert "not a no-go result outside" in certificate.reason


@pytest.mark.parametrize(
    "waves",
    [
        ((0, 0, 0),),
        ((-1, 0, 0),),
        ((1, 0, 0), (1, 0, 0)),
        tuple((index + 1, 0, 0) for index in range(17)),
    ],
)
def test_invalid_custom_search_scopes_are_rejected(waves) -> None:
    with pytest.raises(ValueError):
        search_exact_carrier_gadget(waves)
