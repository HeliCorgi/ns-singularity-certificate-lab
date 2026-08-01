from __future__ import annotations

from fractions import Fraction
import math

import pytest

from ns_certificate_lab.exact_carrier_search import (
    ModeSpec,
    _build_interaction,
    _mode_band,
)
from ns_certificate_lab.expanded_carrier_search import (
    ExpandedCarrierSearchConfig,
    canonical_waves_in_box,
    primitive_polarizations,
    projected_mixed_channels,
    search_expanded_known_relay,
)


def test_bounded_wave_and_primitive_polarization_catalog_is_canonical() -> None:
    waves = canonical_waves_in_box(2)
    assert len(waves) == 62
    assert len(canonical_waves_in_box(3)) == 171
    assert sum(len(primitive_polarizations(wave, 2)) for wave in waves) == 372
    for wave in waves:
        assert max(abs(value) for value in wave) <= 2
        assert next(value for value in wave if value) > 0
        directions = primitive_polarizations(wave, 2)
        assert len(directions) == len(set(directions))
        for vector in directions:
            assert sum(a * b for a, b in zip(wave, vector)) == 0
            assert math.gcd(
                math.gcd(abs(vector[0]), abs(vector[1])), abs(vector[2])
            ) == 1
            assert next(value for value in vector if value) > 0


@pytest.mark.parametrize(
    ("left", "right"),
    (
        (
            ModeSpec((1, 1, 0), (0, 0, 1), "sin"),
            ModeSpec((1, 0, 1), (0, 1, 0), "cos"),
        ),
        (
            ModeSpec((0, 1, -1), (1, 2, 2), "cos"),
            ModeSpec((1, 0, -1), (2, 1, 2), "cos"),
        ),
    ),
)
def test_fast_projected_zero_and_norms_match_real_trigonometric_algebra(
    left: ModeSpec,
    right: ModeSpec,
) -> None:
    fast = projected_mixed_channels(
        left.wave,
        left.polarization,
        right.wave,
        right.polarization,
    )
    exact = _build_interaction(
        left.field(),
        right.field(),
        left_wave=left.wave,
        right_wave=right.wave,
    )
    assert exact is not None
    for prefix in ("sum", "difference"):
        vector = fast[f"{prefix}_vector"]
        wave = fast[f"{prefix}_wave"]
        assert isinstance(vector, tuple) and isinstance(wave, tuple)
        canonical_wave = next(
            candidate
            for candidate in exact.mixed_output.coefficient_table()
            if candidate == exact.__getattribute__(f"{prefix}_wave")
        ) if any(vector) else exact.__getattribute__(f"{prefix}_wave")
        band = _mode_band(exact.mixed_output, canonical_wave)
        fast_l2_squared = sum(value * value for value in vector) / 8
        assert Fraction(band.sobolev_sq(0)) == fast_l2_squared


def test_expanded_bound_two_search_is_exhaustive_and_finds_two_labelled_roles() -> None:
    payload = search_expanded_known_relay()
    assert payload["classification"] == "strict_finite_carrier_hit"
    assert payload["scope"]["exhaustive"] is True
    assert payload["counts"] == {
        "ordered_second_parent_wave_pairs_examined": 3540,
        "scope_eligible_ordered_wave_pairs": 4,
        "polarization_pair_tests": 256,
        "strict_hit_orientation_count": 2,
        "unique_strict_hit_count": 2,
        "unique_cross_talk_partial_count": 10,
        "saved_top_partial_count": 10,
        "strict_hits_without_lower_shell_return": 0,
        "saved_structural_record_count": 256,
    }
    assert payload["rejection_histogram"] == {
        "second_sum_leray_zero": 16,
        "second_difference_leray_nonzero": 216,
        "next_sum_leray_zero": 12,
        "next_difference_leray_zero": 0,
        "diagonal_cross_talk_hits_child_shell": 10,
        "strict_hit_orientation": 2,
    }
    assert len(payload["all_structural_orientation_records"]) == 256
    assert [
        record["record_index"]
        for record in payload["all_structural_orientation_records"]
    ] == list(range(1, 257))
    assert payload["exact_arithmetic_stats_over_full_fast_search"] == {
        "maximum_denominator": 21,
        "maximum_denominator_bit_length": 5,
        "maximum_numerator_bit_length": 8,
        "fraction_component_count": 1824,
    }


def test_unit_normalized_strict_hit_tables_and_low_return_obstruction_are_exact() -> None:
    payload = search_expanded_known_relay()
    assert len(payload["strict_hits"]) == 2
    preferred = payload["strict_hits"][0]
    assert preferred["second_left"]["wave"] == [0, 1, -1]
    assert preferred["second_left"]["polarization"] == [1, 2, 2]
    assert preferred["second_right"]["wave"] == [1, 0, -1]
    assert preferred["second_right"]["polarization"] == [2, 1, 2]
    for side in ("left", "right"):
        normalization = preferred["amplitude_normalization"][side]
        assert normalization["primitive_norm_squared"] == 9
        assert normalization["field_scale_applied"] == "1/3"
        assert normalization["unit_normalization_exactly_rational"] is True
    assert preferred["second_relay"]["interaction"][
        "signed_parent_to_child_flux"
    ] == "25/486"
    assert preferred["next_interaction"][
        "signed_parent_to_child_flux"
    ] == "125/10206"
    assert preferred["next_interaction"]["populated_triad_energy_defect"] == "0"
    assert preferred["full_gadget"]["full_field_energy_defect"] == "0"
    assert preferred["full_gadget"]["grandchild_mutual_output_l2_squared"] == "0"
    assert preferred["full_gadget"]["simple_binary_recursion_rejected"] is True
    assert preferred["diagonal_cross_talk"]["passes"] is True
    audit = preferred["all_unintended_parent_pair_audit"]
    assert audit["all_target_child_shell_outputs_zero"] is True
    assert audit["lower_shell_return_nonzero"] is True
    assert all(pair["target_shell_zero"] for pair in audit["pairs"])
    assert preferred["four_parent_forcing_budget"] == {
        "parent_l2_squared": "2",
        "full_nonlinear_l2_squared": "565/486",
        "target_child_shell_l2_squared": "53/243",
        "intended_C1_C2_mode_l2_squared": "53/243",
        "target_shell_contamination_l2_squared": "0",
        "low_shell_l2_squared": "17/18",
        "high_shell_l2_squared": "0",
        "target_shell_forcing_fraction": "106/565",
        "intended_mode_forcing_fraction": "106/565",
        "low_to_intended_power_ratio": "459/106",
        "energy_pairing": "0",
        "energy_identity_exact": True,
        "classification": "TARGET-SHELL-CLEAN-BUT-LOW-RETURN-DOMINATED",
    }
    assert payload["relay_quality_screen"]["passing_orientation_count"] == 0


def test_parametric_family_and_affine_phase_witness_are_recorded() -> None:
    family = search_expanded_known_relay()["parametric_rational_family"]
    assert family["exact_formulas"] == {
        "relay_two_signed_flux": "25*a^2*d^2/6",
        "next_sum_signed_flux": "125*a^2*d^2/126",
        "next_difference_l2_squared": "35*a^2*d^2/18",
        "next_difference_to_sum_norm_ratio": "7/5",
        "full_energy_defect": "0",
    }
    assert family["affine_phase_witness"] == {
        "units": "quarter turns (pi/2)",
        "x0": [-1, 0, 0],
        "alpha_A": 0,
        "alpha_B": 1,
        "resulting_phases": {
            "A1_at_p": "sin",
            "A2_at_r": "cos",
            "B1_at_q": "cos",
            "B2_at_s": "cos",
        },
    }

    classification = search_expanded_known_relay()[
        "continuous_polarization_classification"
    ]
    assert classification["all_cross_target_shell_zero_equations"][
        "equivalent_factored_system"
    ] == ["a*(e+4*f)=0", "b*(2*e-f)=0"]
    assert [
        branch["constraints"]
        for branch in classification["complete_nonzero_branches"]
    ] == [["a=0", "f=2*e"], ["b=0", "e=-4*f"]]
    assert classification["unit_normalized_equivalence"] == {
        "branches_related_by_cubic_symmetry": True,
        "relay_two_signed_flux": "25/486",
        "next_sum_signed_flux": "125/10206",
        "lower_shell_return_eliminated": False,
    }


def test_coordinate_size_one_search_is_an_exact_negative_control() -> None:
    payload = search_expanded_known_relay(
        ExpandedCarrierSearchConfig(polarization_component_bound=1)
    )
    assert payload["classification"] == "finite_scope_negative"
    assert payload["counts"]["polarization_pair_tests"] == 64
    assert payload["counts"]["strict_hit_orientation_count"] == 0


def test_resource_caps_fail_before_an_incomplete_certificate_is_claimed() -> None:
    with pytest.raises(ValueError, match="ordered wave-pair resource cap"):
        search_expanded_known_relay(
            ExpandedCarrierSearchConfig(maximum_ordered_wave_pair_tests=1)
        )
    with pytest.raises(ValueError, match="polarization-pair resource cap"):
        search_expanded_known_relay(
            ExpandedCarrierSearchConfig(maximum_polarization_pair_tests=1)
        )
    with pytest.raises(ValueError, match="structural-record save cap"):
        search_expanded_known_relay(
            ExpandedCarrierSearchConfig(maximum_saved_structural_records=1)
        )
