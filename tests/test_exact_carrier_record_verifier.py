from __future__ import annotations

import copy

import pytest

from ns_certificate_lab.exact_carrier_record_verifier import (
    verify_serialized_expanded_carrier_certificate,
    verify_serialized_strict_orientation_records,
)
from ns_certificate_lab.exact_carrier_search import search_exact_carrier_gadget
from ns_certificate_lab.expanded_carrier_search import search_expanded_known_relay


def _payload() -> dict[str, object]:
    return search_exact_carrier_gadget().as_dict()


def test_all_sixteen_strict_records_are_independently_recomputed() -> None:
    verdict = verify_serialized_strict_orientation_records(_payload())
    assert verdict == {
        "implementation": (
            "independent-complex-Fourier-ordered-pair-Fraction-convolution"
        ),
        "records_verified": 16,
        "strict_pass_count": 0,
        "energy_identities_verified": 25,
        "rejection_histogram": {
            "diagonal_parent_cross_talk_hits_target_child_shell": 12,
            "next_child_interaction_lacks_required_nonzero_sum_or_difference": 8,
        },
        "verified": True,
    }


def test_coefficient_table_fault_injection_is_detected() -> None:
    payload = copy.deepcopy(_payload())
    record = payload["strict_search"]["orientation_records"][0]
    interaction = record["relay_two"]["interaction"]
    first_wave = next(iter(interaction["sum_output"]))
    interaction["sum_output"][first_wave]["cosine"][0] = "999"
    with pytest.raises(ValueError, match="coefficient table"):
        verify_serialized_strict_orientation_records(payload)


def test_rejection_reason_fault_injection_is_detected() -> None:
    payload = copy.deepcopy(_payload())
    payload["strict_search"]["orientation_records"][0][
        "rejection_reasons"
    ] = []
    with pytest.raises(ValueError, match="rejection predicate"):
        verify_serialized_strict_orientation_records(payload)


def test_missing_orientation_record_invalidates_exhaustive_inventory() -> None:
    payload = copy.deepcopy(_payload())
    payload["strict_search"]["orientation_records"].pop()
    with pytest.raises(ValueError, match="inventory"):
        verify_serialized_strict_orientation_records(payload)


def test_expanded_certificate_is_independently_recomputed() -> None:
    verdict = verify_serialized_expanded_carrier_certificate(
        search_expanded_known_relay()
    )
    assert verdict == {
        "implementation": (
            "independent-complex-Fourier-ordered-pair-Fraction-convolution"
        ),
        "structural_records_verified": 256,
        "strict_exact_witnesses_verified": 2,
        "partial_exact_witnesses_verified": 10,
        "full_search_arithmetic_stats_verified": {
            "maximum_denominator": 21,
            "maximum_denominator_bit_length": 5,
            "maximum_numerator_bit_length": 8,
            "fraction_component_count": 1824,
        },
        "all_strict_hits_have_lower_shell_return": True,
        "verified": True,
    }


def test_expanded_scaled_parent_fault_injection_is_detected() -> None:
    payload = copy.deepcopy(search_expanded_known_relay())
    payload["strict_hits"][0]["amplitude_normalization"]["left"][
        "field_scale_applied"
    ] = "1"
    with pytest.raises(ValueError, match="field polarization"):
        verify_serialized_expanded_carrier_certificate(payload)


def test_expanded_structural_failure_fault_injection_is_detected() -> None:
    payload = copy.deepcopy(search_expanded_known_relay())
    payload["all_structural_orientation_records"][0][
        "disposition"
    ] = "strict_hit_orientation"
    with pytest.raises(ValueError, match="structural record"):
        verify_serialized_expanded_carrier_certificate(payload)


def test_expanded_full_energy_fault_injection_is_detected() -> None:
    payload = copy.deepcopy(search_expanded_known_relay())
    payload["strict_hits"][0]["full_gadget"][
        "full_field_energy_defect"
    ] = "1"
    with pytest.raises(ValueError, match="full gadget"):
        verify_serialized_expanded_carrier_certificate(payload)
