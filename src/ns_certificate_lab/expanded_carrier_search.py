r"""Expanded exact search around the known periodic Leray relay.

The search fixes the known relay

``(1,1,0), e3*sin`` + ``(1,0,1), e2*cos`` -> ``(2,1,1)``

and exhausts a bounded alphabet of primitive integer polarization directions
for a second relay.  Only exact integer/:class:`fractions.Fraction` arithmetic
is used for the structural prefilter.  Full saved witnesses are rebuilt with
the real trigonometric Fourier algebra in :mod:`exact_carrier_search`.

This is a finite carrier-level algebra certificate.  It is not an invariant
set, an iteratable cascade, a Galerkin orbit, or a Navier--Stokes singularity.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import math
from typing import Iterable

from .exact_carrier_search import (
    ModeSpec,
    RelayPairWitness,
    _build_gadget,
    _build_interaction,
    _canonical_wave,
    _cross_talk_diagnostic,
    _field_table,
    _is_zero,
    _mixed_output,
    _shell_band,
)
from .fourier_torus import advection, leray
from .torus_chain import l2_inner

__all__ = [
    "ExpandedCarrierSearchConfig",
    "canonical_waves_in_box",
    "primitive_polarizations",
    "projected_mixed_channels",
    "search_expanded_known_relay",
]


Wave = tuple[int, int, int]
IntVector = tuple[int, int, int]
RatVector = tuple[Fraction, Fraction, Fraction]

SCHEMA = "ns-certificate-lab/expanded-exact-carrier-search/v1"
STATUS_FOUND = "STRICT EXACT FINITE CARRIER HIT / NOT A PDE OR CASCADE"
STATUS_NEGATIVE = "NEGATIVE IN THE STATED EXPANDED FINITE SEARCH ONLY"

_P: Wave = (1, 1, 0)
_Q: Wave = (1, 0, 1)
_C1: Wave = (2, 1, 1)
_E3: IntVector = (0, 0, 1)
_E2: IntVector = (0, 1, 0)
_ZERO: Wave = (0, 0, 0)


def _dot(left: Iterable[int | Fraction], right: Iterable[int | Fraction]) -> Fraction:
    return sum(
        (Fraction(a) * Fraction(b) for a, b in zip(left, right)),
        Fraction(0),
    )


def _add(left: Wave, right: Wave) -> Wave:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def _subtract(left: Wave, right: Wave) -> Wave:
    return tuple(a - b for a, b in zip(left, right))  # type: ignore[return-value]


def _norm_squared(vector: Iterable[int | Fraction]) -> Fraction:
    return _dot(vector, vector)


def _project(vector: Iterable[int | Fraction], wave: Wave) -> RatVector:
    denominator = Fraction(sum(component * component for component in wave))
    if denominator == 0:
        raise ValueError("cannot apply the Leray projection at zero wave")
    values = tuple(Fraction(value) for value in vector)
    longitudinal = _dot(wave, values) / denominator
    return tuple(
        value - Fraction(component) * longitudinal
        for value, component in zip(values, wave)
    )  # type: ignore[return-value]


def _zero(vector: RatVector) -> bool:
    return all(value == 0 for value in vector)


def _canonical_direction(vector: IntVector) -> IntVector:
    for component in vector:
        if component > 0:
            return vector
        if component < 0:
            return tuple(-value for value in vector)  # type: ignore[return-value]
    raise ValueError("zero has no polarization direction")


def canonical_waves_in_box(bound: int) -> tuple[Wave, ...]:
    """Return every nonzero half-lattice wave with ``|k|_inf <= bound``."""

    if isinstance(bound, bool) or not isinstance(bound, int) or not 1 <= bound <= 3:
        raise ValueError("wave bound must be an integer in [1, 3]")
    waves: list[Wave] = []
    for raw in product(range(-bound, bound + 1), repeat=3):
        wave: Wave = tuple(raw)  # type: ignore[assignment]
        if wave == _ZERO:
            continue
        canonical, sign = _canonical_wave(wave)
        if sign > 0 and canonical == wave:
            waves.append(wave)
    return tuple(sorted(waves))


def primitive_polarizations(wave: Wave, component_bound: int) -> tuple[IntVector, ...]:
    """Enumerate sign-quotiented primitive integer directions perpendicular to k."""

    if (
        isinstance(component_bound, bool)
        or not isinstance(component_bound, int)
        or not 1 <= component_bound <= 4
    ):
        raise ValueError("polarization bound must be an integer in [1, 4]")
    directions: set[IntVector] = set()
    for raw in product(range(-component_bound, component_bound + 1), repeat=3):
        vector: IntVector = tuple(raw)  # type: ignore[assignment]
        if vector == _ZERO or sum(a * b for a, b in zip(wave, vector)) != 0:
            continue
        divisor = math.gcd(math.gcd(abs(vector[0]), abs(vector[1])), abs(vector[2]))
        if divisor != 1:
            continue
        directions.add(_canonical_direction(vector))
    return tuple(sorted(directions))


def projected_mixed_channels(
    left_wave: Wave,
    left_polarization: Iterable[int | Fraction],
    right_wave: Wave,
    right_polarization: Iterable[int | Fraction],
) -> dict[str, object]:
    """Exact phase-quotiented sum/difference Leray channel vectors.

    Cosine/sine choices only change the trigonometric phase and an overall
    sign of these two vectors.  They do not change whether either channel is
    zero or the squared radius of its support.
    """

    left = tuple(Fraction(value) for value in left_polarization)
    right = tuple(Fraction(value) for value in right_polarization)
    left_on_right = _dot(left, right_wave)
    right_on_left = _dot(right, left_wave)
    sum_raw = tuple(
        left_on_right * right_component + right_on_left * left_component
        for left_component, right_component in zip(left, right)
    )
    difference_raw = tuple(
        left_on_right * right_component - right_on_left * left_component
        for left_component, right_component in zip(left, right)
    )
    sum_wave = _add(left_wave, right_wave)
    difference_wave = _subtract(left_wave, right_wave)
    if sum_wave == _ZERO or difference_wave == _ZERO:
        raise ValueError("mixed channel requires distinct non-opposite waves")
    return {
        "sum_wave": sum_wave,
        "difference_wave": difference_wave,
        "sum_vector": _project(sum_raw, sum_wave),
        "difference_vector": _project(difference_raw, difference_wave),
    }


@dataclass(frozen=True)
class ExpandedCarrierSearchConfig:
    """Finite resource limits and alphabet bounds for the exact search."""

    wave_bound: int = 2
    polarization_component_bound: int = 2
    maximum_ordered_wave_pair_tests: int = 10_000
    maximum_polarization_pair_tests: int = 200_000
    maximum_saved_structural_records: int = 5_000
    maximum_saved_unique_hits: int = 64
    top_partial_limit: int = 16

    def __post_init__(self) -> None:
        canonical_waves_in_box(self.wave_bound)
        primitive_polarizations((1, 0, 0), self.polarization_component_bound)
        for name, value in (
            ("maximum_ordered_wave_pair_tests", self.maximum_ordered_wave_pair_tests),
            ("maximum_polarization_pair_tests", self.maximum_polarization_pair_tests),
            (
                "maximum_saved_structural_records",
                self.maximum_saved_structural_records,
            ),
            ("maximum_saved_unique_hits", self.maximum_saved_unique_hits),
            ("top_partial_limit", self.top_partial_limit),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")


@dataclass(frozen=True)
class _FastCandidate:
    left_wave: Wave
    left_polarization: IntVector
    right_wave: Wave
    right_polarization: IntVector
    diagonal_leakage_squared: Fraction
    next_sum_strength_squared: Fraction

    def unordered_key(self) -> tuple[tuple[Wave, IntVector], tuple[Wave, IntVector]]:
        return tuple(
            sorted(
                (
                    (self.left_wave, self.left_polarization),
                    (self.right_wave, self.right_polarization),
                )
            )
        )  # type: ignore[return-value]

    def ordering_key(self) -> tuple[object, ...]:
        return (
            self.diagonal_leakage_squared,
            -self.next_sum_strength_squared,
            self.left_wave,
            self.left_polarization,
            self.right_wave,
            self.right_polarization,
        )


def _channel_strength(vector: RatVector) -> Fraction:
    return _norm_squared(vector)


def _vector_strings(vector: RatVector) -> list[str]:
    return [str(value) for value in vector]


def _update_fraction_stats(stats: dict[str, int], value: object) -> None:
    values: Iterable[object]
    if isinstance(value, tuple):
        values = value
    else:
        values = (value,)
    for item in values:
        if not isinstance(item, Fraction):
            continue
        stats["maximum_denominator"] = max(
            stats["maximum_denominator"], item.denominator
        )
        stats["maximum_denominator_bit_length"] = max(
            stats["maximum_denominator_bit_length"], item.denominator.bit_length()
        )
        stats["maximum_numerator_bit_length"] = max(
            stats["maximum_numerator_bit_length"], abs(item.numerator).bit_length()
        )
        stats["fraction_component_count"] += 1


def _target_shell_leakage(
    left_wave: Wave,
    left_polarization: Iterable[int | Fraction],
    right_wave: Wave,
    right_polarization: Iterable[int | Fraction],
    *,
    target_shell_squared: int,
    fraction_stats: dict[str, int] | None = None,
) -> Fraction:
    channels = projected_mixed_channels(
        left_wave,
        left_polarization,
        right_wave,
        right_polarization,
    )
    leakage = Fraction(0)
    for prefix in ("sum", "difference"):
        wave = channels[f"{prefix}_wave"]
        vector = channels[f"{prefix}_vector"]
        if not isinstance(wave, tuple) or not isinstance(vector, tuple):
            raise AssertionError("internal channel type mismatch")
        if fraction_stats is not None:
            _update_fraction_stats(fraction_stats, vector)
        if sum(component * component for component in wave) == target_shell_squared:
            leakage += _channel_strength(vector)  # phase-independent scale
    return leakage


def _maximum_denominator(value: object) -> int:
    maximum = 1
    if isinstance(value, str):
        try:
            maximum = Fraction(value).denominator
        except (ValueError, ZeroDivisionError):
            pass
    elif isinstance(value, dict):
        maximum = max((_maximum_denominator(item) for item in value.values()), default=1)
    elif isinstance(value, list):
        maximum = max((_maximum_denominator(item) for item in value), default=1)
    return maximum


def _rational_unit_scale(vector: IntVector) -> Fraction | None:
    norm_squared = sum(value * value for value in vector)
    norm = math.isqrt(norm_squared)
    return Fraction(1, norm) if norm * norm == norm_squared else None


def _normalization_payload(
    vector: IntVector,
    scale: Fraction,
    *,
    requested: bool,
) -> dict[str, object]:
    norm_squared = sum(value * value for value in vector)
    norm = math.isqrt(norm_squared)
    return {
        "primitive_polarization": list(vector),
        "primitive_norm_squared": norm_squared,
        "primitive_norm": str(norm) if norm * norm == norm_squared else f"sqrt({norm_squared})",
        "field_scale_applied": str(scale),
        "unit_normalization_requested": requested,
        "unit_normalization_exactly_rational": norm * norm == norm_squared,
        "field_polarization": [str(scale * value) for value in vector],
    }


def _unintended_pair_payload(
    label: str,
    left,
    right,
    *,
    target_shell_squared: int,
) -> dict[str, object]:
    output = _mixed_output(left, right)
    target = _shell_band(output, target_shell_squared)
    shells = sorted(
        {
            sum(component * component for component in wave)
            for wave in output.coefficient_table()
        }
    )
    return {
        "label": label,
        "full_output": _field_table(output),
        "target_shell_output": _field_table(target),
        "target_shell_zero": _is_zero(target),
        "nonzero_output_shells_squared": shells,
    }


def _parent_forcing_budget(known: RelayPairWitness, second: RelayPairWitness) -> dict[str, object]:
    parent = (
        known.interaction.left
        + known.interaction.right
        + second.interaction.left
        + second.interaction.right
    ).cleaned()
    nonlinear = leray(advection(parent, parent)).cleaned()
    target_shell = _shell_band(nonlinear, 6)
    intended_waves = {
        known.interaction.sum_wave,
        second.interaction.sum_wave,
    }
    intended = nonlinear.restrict(lambda wave: wave in intended_waves).cleaned()
    low = nonlinear.restrict(
        lambda wave: sum(component * component for component in wave) < 6
    ).cleaned()
    high = nonlinear.restrict(
        lambda wave: sum(component * component for component in wave) > 6
    ).cleaned()
    parent_norm = Fraction(parent.sobolev_sq(0))
    full_power = Fraction(nonlinear.sobolev_sq(0))
    target_shell_power = Fraction(target_shell.sobolev_sq(0))
    intended_power = Fraction(intended.sobolev_sq(0))
    low_power = Fraction(low.sobolev_sq(0))
    high_power = Fraction(high.sobolev_sq(0))
    pairing = Fraction(l2_inner(parent, nonlinear))
    return {
        "parent_l2_squared": str(parent_norm),
        "full_nonlinear_l2_squared": str(full_power),
        "target_child_shell_l2_squared": str(target_shell_power),
        "intended_C1_C2_mode_l2_squared": str(intended_power),
        "target_shell_contamination_l2_squared": str(
            target_shell_power - intended_power
        ),
        "low_shell_l2_squared": str(low_power),
        "high_shell_l2_squared": str(high_power),
        "target_shell_forcing_fraction": str(target_shell_power / full_power),
        "intended_mode_forcing_fraction": str(intended_power / full_power),
        "low_to_intended_power_ratio": str(low_power / intended_power),
        "energy_pairing": str(pairing),
        "energy_identity_exact": pairing == 0,
        "classification": (
            "TARGET-SHELL-CLEAN-BUT-LOW-RETURN-DOMINATED"
            if target_shell_power == intended_power and low_power > intended_power
            else "TARGET-SHELL-CONTAMINATED"
        ),
    }


def _exact_candidate_payload(
    candidate: _FastCandidate,
    *,
    rational_unit_normalization: bool,
) -> dict[str, object]:
    known_left = ModeSpec(_P, _E3, "sin")
    known_right = ModeSpec(_Q, _E2, "cos")
    known_interaction = _build_interaction(
        known_left.field(),
        known_right.field(),
        left_wave=_P,
        right_wave=_Q,
    )
    if known_interaction is None:
        raise AssertionError("the fixed known relay unexpectedly vanished")
    known = RelayPairWitness(known_left, known_right, known_interaction)

    second_left = ModeSpec(
        candidate.left_wave,
        candidate.left_polarization,
        "cos",
    )
    second_right = ModeSpec(
        candidate.right_wave,
        candidate.right_polarization,
        "cos",
    )
    left_unit_scale = _rational_unit_scale(candidate.left_polarization)
    right_unit_scale = _rational_unit_scale(candidate.right_polarization)
    can_normalize = left_unit_scale is not None and right_unit_scale is not None
    normalization_applied = rational_unit_normalization and can_normalize
    left_scale = left_unit_scale if normalization_applied else Fraction(1)
    right_scale = right_unit_scale if normalization_applied else Fraction(1)
    if left_scale is None or right_scale is None:
        raise AssertionError("internal normalization selection failed")
    second_left_field = second_left.field().scale(left_scale)
    second_right_field = second_right.field().scale(right_scale)
    second_interaction = _build_interaction(
        second_left_field,
        second_right_field,
        left_wave=second_left.wave,
        right_wave=second_right.wave,
    )
    if second_interaction is None or not _is_zero(second_interaction.difference_output):
        raise AssertionError("fast second-relay predicate disagrees with exact algebra")
    second = RelayPairWitness(second_left, second_right, second_interaction)
    next_interaction = _build_interaction(
        known.interaction.child,
        second.interaction.child,
        left_wave=known.interaction.sum_wave,
        right_wave=second.interaction.sum_wave,
    )
    if next_interaction is None or _is_zero(next_interaction.difference_output):
        raise AssertionError("fast next-stage predicate disagrees with exact algebra")
    cross_talk = _cross_talk_diagnostic(known, second, 6)
    exact_leakage_zero = cross_talk.passes
    if exact_leakage_zero != (candidate.diagonal_leakage_squared == 0):
        raise AssertionError("fast diagonal predicate disagrees with exact algebra")
    gadget = _build_gadget(known, second, next_interaction, (2, 6), cross_talk)
    unintended_pairs = [
        _unintended_pair_payload(
            "A1+A2",
            known.interaction.left,
            second.interaction.left,
            target_shell_squared=6,
        ),
        _unintended_pair_payload(
            "B1+B2",
            known.interaction.right,
            second.interaction.right,
            target_shell_squared=6,
        ),
        _unintended_pair_payload(
            "A1+B2",
            known.interaction.left,
            second.interaction.right,
            target_shell_squared=6,
        ),
        _unintended_pair_payload(
            "A2+B1",
            second.interaction.left,
            known.interaction.right,
            target_shell_squared=6,
        ),
    ]
    all_unintended_target_zero = all(
        bool(record["target_shell_zero"]) for record in unintended_pairs
    )
    has_lower_shell_return = any(
        any(int(shell) < 6 for shell in record["nonzero_output_shells_squared"])
        for record in unintended_pairs
    )
    payload: dict[str, object] = {
        "representative_phase_choice": {
            "fixed_left": "sin",
            "fixed_right": "cos",
            "second_left": "cos",
            "second_right": "cos",
        },
        "second_left": second_left.as_dict(),
        "second_right": second_right.as_dict(),
        "amplitude_normalization": {
            "policy": (
                "rational unit normalization"
                if normalization_applied
                else (
                    "primitive exact representative; requested unit scale is algebraic"
                    if rational_unit_normalization
                    else "primitive integer representative"
                )
            ),
            "reconstruction_rule": (
                "multiply each primitive ModeSpec field by its recorded "
                "field_scale_applied before recomputing interactions"
            ),
            "left": _normalization_payload(
                candidate.left_polarization,
                left_scale,
                requested=rational_unit_normalization,
            ),
            "right": _normalization_payload(
                candidate.right_polarization,
                right_scale,
                requested=rational_unit_normalization,
            ),
        },
        "exact_scaled_second_parent_fields": {
            "left": _field_table(second.interaction.left),
            "right": _field_table(second.interaction.right),
        },
        "diagonal_target_shell_leakage_strength_squared": str(
            candidate.diagonal_leakage_squared
        ),
        "phase_quotiented_next_sum_strength_squared": str(
            candidate.next_sum_strength_squared
        ),
        "second_relay": second.as_dict(),
        "next_interaction": next_interaction.as_dict(),
        "diagonal_cross_talk": cross_talk.as_dict(),
        "all_unintended_parent_pair_audit": {
            "pairs": unintended_pairs,
            "all_target_child_shell_outputs_zero": all_unintended_target_zero,
            "lower_shell_return_nonzero": has_lower_shell_return,
            "classification": (
                "target-shell-clean-but-lower-shell-return-present"
                if all_unintended_target_zero and has_lower_shell_return
                else "target-shell-contaminated-or-no-lower-return"
            ),
        },
        "four_parent_forcing_budget": _parent_forcing_budget(known, second),
        "full_gadget": gadget.as_dict(),
        "exact_predicates": {
            "second_difference_zero": True,
            "second_signed_transfer_positive": (
                second_interaction.signed_parent_to_child_flux > 0
            ),
            "next_sum_nonzero": not _is_zero(next_interaction.sum_output),
            "next_difference_nonzero": not _is_zero(
                next_interaction.difference_output
            ),
            "next_signed_transfer_positive": (
                next_interaction.signed_parent_to_child_flux > 0
            ),
            "diagonal_target_shell_zero": cross_talk.passes,
            "all_unintended_parent_pairs_target_shell_zero": (
                all_unintended_target_zero
            ),
            "relay_two_energy_defect_zero": (
                second_interaction.populated_triad_energy_defect == 0
            ),
            "next_energy_defect_zero": (
                next_interaction.populated_triad_energy_defect == 0
            ),
            "full_field_energy_defect_zero": gadget.full_field_energy_defect == 0,
        },
    }
    payload["maximum_saved_fraction_denominator"] = _maximum_denominator(payload)
    return payload


def _parametric_family_payload() -> dict[str, object]:
    """The exact rational family exposed by the first bounded strict hit."""

    return {
        "waves": {
            "r": [0, 1, -1],
            "s": [1, 0, -1],
            "C2": [1, 1, -2],
        },
        "polarizations": {
            "R": "a*(1,2,2)",
            "S": "d*(2,1,2)",
            "parameters": "a,d are arbitrary nonzero rational numbers",
        },
        "derivation_constraints": {
            "general_perpendicular_forms": "R=(a,b,b), S=(c,d,c)",
            "second_difference_leray_zero": "a*c=b*d",
            "A1_plus_B2_target_shell_zero": "c=2*d",
            "solution_family": "b=2*a, c=2*d",
        },
        "exact_formulas": {
            "relay_two_signed_flux": "25*a^2*d^2/6",
            "next_sum_signed_flux": "125*a^2*d^2/126",
            "next_difference_l2_squared": "35*a^2*d^2/18",
            "next_difference_to_sum_norm_ratio": "7/5",
            "full_energy_defect": "0",
        },
        "unit_normalized_representative": {
            "a": "1/3",
            "d": "1/3",
            "R": ["1/3", "2/3", "2/3"],
            "S": ["2/3", "1/3", "2/3"],
            "relay_two_signed_flux": "25/486",
            "next_sum_signed_flux": "125/10206",
        },
        "affine_phase_witness": {
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
        },
        "obstructions_still_present": {
            "diagonal_outputs_globally_zero": False,
            "lower_shell_return": True,
            "simple_aligned_grandchild_recursion_zero": True,
            "interpretation": (
                "the family clears the intended child shell but is not a closed "
                "or iteratable relay"
            ),
        },
    }


def _continuous_polarization_classification_payload() -> dict[str, object]:
    """Complete algebraic classification on the displayed carrier geometry."""

    return {
        "status": "EXACT SYMBOLIC CLASSIFICATION ON FIXED FOUR-WAVE GEOMETRY",
        "geometry": {
            "p": [1, 1, 0],
            "q": [1, 0, 1],
            "r": [0, 1, -1],
            "s": [1, 0, -1],
        },
        "general_perpendicular_polarizations": {
            "P": "(a,-a,b)",
            "Q": "t*(a,b,-a)",
            "R": "(e,f,f)",
            "S": "u*(f,e,f)",
            "nonzero_parameters": "t,u and both relay child amplitudes are nonzero",
        },
        "relay_identities": {
            "relay_one_difference": "0",
            "relay_one_child_scalar": "2*t*(2*a-b)*(a+b)/3",
            "relay_two_difference": "0",
            "relay_two_child_scalar": "2*u*(e+2*f)*(e-f)/3",
        },
        "all_cross_target_shell_zero_equations": {
            "F_minus": "a*e+4*a*f-2*b*e+b*f=0",
            "F_plus": "a*e+4*a*f+2*b*e-b*f=0",
            "equivalent_factored_system": [
                "a*(e+4*f)=0",
                "b*(2*e-f)=0",
            ],
        },
        "complete_nonzero_branches": [
            {
                "branch": "I",
                "constraints": ["a=0", "f=2*e"],
                "representative": {
                    "P": [0, 0, 1],
                    "Q": [0, 1, 0],
                    "R": [1, 2, 2],
                    "S": [2, 1, 2],
                },
            },
            {
                "branch": "II",
                "constraints": ["b=0", "e=-4*f"],
                "representative": {
                    "P": [1, -1, 0],
                    "Q": [1, 0, -1],
                    "R": [-4, 1, 1],
                    "S": [1, -4, 1],
                },
            },
        ],
        "completeness_argument": {
            "a_and_b_both_nonzero": (
                "e+4*f=0 and 2*e-f=0 have determinant -9, hence e=f=0, "
                "contradicting a nonzero second relay"
            ),
            "remaining_cases": "a=0 gives branch I; b=0 gives branch II",
        },
        "unit_normalized_equivalence": {
            "branches_related_by_cubic_symmetry": True,
            "relay_two_signed_flux": "25/486",
            "next_sum_signed_flux": "125/10206",
            "lower_shell_return_eliminated": False,
        },
        "classification": (
            "TWO CUBICALLY EQUIVALENT STRICT BRANCHES; NO GENERIC THIRD BRANCH"
        ),
        "claim_boundary": (
            "classification fixes the four carrier wavevectors and says nothing "
            "about other geometries or dynamical closure"
        ),
    }


def _fixed_relay_payload() -> dict[str, object]:
    left = ModeSpec(_P, _E3, "sin")
    right = ModeSpec(_Q, _E2, "cos")
    interaction = _build_interaction(
        left.field(), right.field(), left_wave=_P, right_wave=_Q
    )
    if interaction is None:
        raise AssertionError("known relay is missing")
    return RelayPairWitness(left, right, interaction).as_dict()


def search_expanded_known_relay(
    config: ExpandedCarrierSearchConfig | None = None,
) -> dict[str, object]:
    """Exhaust the fully shell-coherent extension of the fixed known relay."""

    config = ExpandedCarrierSearchConfig() if config is None else config
    waves = canonical_waves_in_box(config.wave_bound)
    polarizations = {
        wave: primitive_polarizations(wave, config.polarization_component_bound)
        for wave in waves
    }
    second_waves = tuple(wave for wave in waves if wave not in (_P, _Q))
    all_ordered_wave_pairs = tuple(
        (left, right)
        for left in second_waves
        for right in second_waves
        if left != right
    )
    if len(all_ordered_wave_pairs) > config.maximum_ordered_wave_pair_tests:
        raise ValueError("ordered wave-pair resource cap exceeded before search")

    scope_wave_pairs = tuple(
        (left, right)
        for left, right in all_ordered_wave_pairs
        if sum(value * value for value in left) == 2
        and sum(value * value for value in right) == 2
        and sum(value * value for value in _add(left, right)) == 6
        and _add(left, right) != _C1
    )
    polarization_pair_tests = sum(
        len(polarizations[left]) * len(polarizations[right])
        for left, right in scope_wave_pairs
    )
    if polarization_pair_tests > config.maximum_polarization_pair_tests:
        raise ValueError("polarization-pair resource cap exceeded before search")
    if polarization_pair_tests > config.maximum_saved_structural_records:
        raise ValueError(
            "structural-record save cap exceeded; refusing an incomplete certificate"
        )

    rejection_histogram = {
        "second_sum_leray_zero": 0,
        "second_difference_leray_nonzero": 0,
        "next_sum_leray_zero": 0,
        "next_difference_leray_zero": 0,
        "diagonal_cross_talk_hits_child_shell": 0,
        "strict_hit_orientation": 0,
    }
    fraction_stats = {
        "maximum_denominator": 1,
        "maximum_denominator_bit_length": 1,
        "maximum_numerator_bit_length": 1,
        "fraction_component_count": 0,
    }
    hit_orientations: list[_FastCandidate] = []
    partials: list[_FastCandidate] = []
    structural_records: list[dict[str, object]] = []
    for left_wave, right_wave in scope_wave_pairs:
        for left_polarization in polarizations[left_wave]:
            for right_polarization in polarizations[right_wave]:
                second = projected_mixed_channels(
                    left_wave,
                    left_polarization,
                    right_wave,
                    right_polarization,
                )
                second_sum = second["sum_vector"]
                second_difference = second["difference_vector"]
                second_child_wave = second["sum_wave"]
                if not isinstance(second_sum, tuple) or not isinstance(
                    second_difference, tuple
                ) or not isinstance(second_child_wave, tuple):
                    raise AssertionError("internal channel type mismatch")
                _update_fraction_stats(fraction_stats, second_sum)
                _update_fraction_stats(fraction_stats, second_difference)
                structural_record: dict[str, object] = {
                    "record_index": len(structural_records) + 1,
                    "second_left": {
                        "wave": list(left_wave),
                        "primitive_polarization": list(left_polarization),
                    },
                    "second_right": {
                        "wave": list(right_wave),
                        "primitive_polarization": list(right_polarization),
                    },
                    "second_sum_wave": list(second_child_wave),
                    "second_sum_projected_vector": _vector_strings(second_sum),
                    "second_difference_wave": list(second["difference_wave"]),
                    "second_difference_projected_vector": _vector_strings(
                        second_difference
                    ),
                    "next_sum_wave": None,
                    "next_sum_projected_vector": None,
                    "next_difference_wave": None,
                    "next_difference_projected_vector": None,
                    "diagonal_target_shell_leakage_strength_squared": None,
                    "disposition": None,
                }
                if _zero(second_sum):
                    rejection_histogram["second_sum_leray_zero"] += 1
                    structural_record["disposition"] = "second_sum_leray_zero"
                    structural_records.append(structural_record)
                    continue
                if not _zero(second_difference):
                    rejection_histogram["second_difference_leray_nonzero"] += 1
                    structural_record[
                        "disposition"
                    ] = "second_difference_leray_nonzero"
                    structural_records.append(structural_record)
                    continue
                next_channels = projected_mixed_channels(
                    _C1,
                    projected_mixed_channels(_P, _E3, _Q, _E2)["sum_vector"],
                    second_child_wave,
                    second_sum,
                )
                next_sum = next_channels["sum_vector"]
                next_difference = next_channels["difference_vector"]
                if not isinstance(next_sum, tuple) or not isinstance(
                    next_difference, tuple
                ):
                    raise AssertionError("internal channel type mismatch")
                _update_fraction_stats(fraction_stats, next_sum)
                _update_fraction_stats(fraction_stats, next_difference)
                structural_record["next_sum_wave"] = list(
                    next_channels["sum_wave"]
                )
                structural_record["next_sum_projected_vector"] = _vector_strings(
                    next_sum
                )
                structural_record["next_difference_wave"] = list(
                    next_channels["difference_wave"]
                )
                structural_record[
                    "next_difference_projected_vector"
                ] = _vector_strings(next_difference)
                if _zero(next_sum):
                    rejection_histogram["next_sum_leray_zero"] += 1
                    structural_record["disposition"] = "next_sum_leray_zero"
                    structural_records.append(structural_record)
                    continue
                if _zero(next_difference):
                    rejection_histogram["next_difference_leray_zero"] += 1
                    structural_record[
                        "disposition"
                    ] = "next_difference_leray_zero"
                    structural_records.append(structural_record)
                    continue
                leakage = _target_shell_leakage(
                    _P,
                    _E3,
                    right_wave,
                    right_polarization,
                    target_shell_squared=6,
                    fraction_stats=fraction_stats,
                ) + _target_shell_leakage(
                    left_wave,
                    left_polarization,
                    _Q,
                    _E2,
                    target_shell_squared=6,
                    fraction_stats=fraction_stats,
                )
                candidate = _FastCandidate(
                    left_wave=left_wave,
                    left_polarization=left_polarization,
                    right_wave=right_wave,
                    right_polarization=right_polarization,
                    diagonal_leakage_squared=leakage,
                    next_sum_strength_squared=_channel_strength(next_sum),
                )
                structural_record[
                    "diagonal_target_shell_leakage_strength_squared"
                ] = str(leakage)
                if leakage:
                    rejection_histogram[
                        "diagonal_cross_talk_hits_child_shell"
                    ] += 1
                    partials.append(candidate)
                    structural_record[
                        "disposition"
                    ] = "diagonal_cross_talk_hits_child_shell"
                    structural_records.append(structural_record)
                    continue
                rejection_histogram["strict_hit_orientation"] += 1
                hit_orientations.append(candidate)
                structural_record["disposition"] = "strict_hit_orientation"
                structural_records.append(structural_record)

    if sum(rejection_histogram.values()) != polarization_pair_tests:
        raise AssertionError("rejection histogram does not partition the search")

    unique_hits: dict[
        tuple[Wave, IntVector, Wave, IntVector], _FastCandidate
    ] = {}
    for candidate in sorted(hit_orientations, key=_FastCandidate.ordering_key):
        key = (
            candidate.left_wave,
            candidate.left_polarization,
            candidate.right_wave,
            candidate.right_polarization,
        )
        unique_hits.setdefault(key, candidate)
    if len(unique_hits) > config.maximum_saved_unique_hits:
        raise ValueError("unique strict-hit save cap exceeded; refusing incomplete certificate")

    unique_partials: dict[
        tuple[Wave, IntVector, Wave, IntVector], _FastCandidate
    ] = {}
    for candidate in sorted(partials, key=_FastCandidate.ordering_key):
        key = (
            candidate.left_wave,
            candidate.left_polarization,
            candidate.right_wave,
            candidate.right_polarization,
        )
        unique_partials.setdefault(key, candidate)
    saved_partials = list(unique_partials.values())[: config.top_partial_limit]

    strict_payloads = [
        _exact_candidate_payload(candidate, rational_unit_normalization=True)
        for candidate in unique_hits.values()
    ]
    partial_payloads = [
        _exact_candidate_payload(candidate, rational_unit_normalization=False)
        for candidate in saved_partials
    ]
    found = bool(strict_payloads)
    strict_without_low_return = sum(
        not bool(
            payload["all_unintended_parent_pair_audit"][
                "lower_shell_return_nonzero"
            ]
        )
        for payload in strict_payloads
    )
    status = STATUS_FOUND if found else STATUS_NEGATIVE
    total_polarization_directions = sum(len(values) for values in polarizations.values())
    return {
        "schema": SCHEMA,
        "status": status,
        "classification": (
            "strict_finite_carrier_hit" if found else "finite_scope_negative"
        ),
        "claim_boundary": (
            "Exact finite carrier algebra only. No closure, invariant cone, "
            "scale iteration, PDE trajectory, or singularity is certified."
        ),
        "search_config": {
            "wave_bound_l_infinity": config.wave_bound,
            "primitive_polarization_component_bound": (
                config.polarization_component_bound
            ),
            "input_rational_denominator": 1,
            "maximum_ordered_wave_pair_tests": (
                config.maximum_ordered_wave_pair_tests
            ),
            "maximum_polarization_pair_tests": (
                config.maximum_polarization_pair_tests
            ),
            "maximum_saved_structural_records": (
                config.maximum_saved_structural_records
            ),
            "maximum_saved_unique_hits": config.maximum_saved_unique_hits,
            "top_partial_limit": config.top_partial_limit,
        },
        "phase_quotient": {
            "applied": True,
            "represented_phase_pairs": ["cos/cos", "cos/sin", "sin/cos", "sin/sin"],
            "reason": (
                "each cosine/sine mixture multiplies a given frequency channel "
                "by a nonzero phase/sign scalar; the projected vector's zero "
                "predicate and support shell radius are unchanged. No cancellation "
                "between distinct parent pairs is used"
            ),
            "independent_phase_optimization_used": False,
        },
        "deduplication": {
            "wave_sign": "canonical half lattice",
            "polarization_scale_and_sign": (
                "primitive gcd-one vector with first nonzero component positive"
            ),
            "second_relay_exchange": (
                "not quotiented: A2/B2 exchange changes labelled diagonal pairs"
            ),
        },
        "scope": {
            "name": "fully_shell_coherent_fixed_known_relay",
            "fixed_parent_shell_squared": 2,
            "second_parent_shell_squared": 2,
            "fixed_child_shell_squared": 6,
            "second_child_shell_squared": 6,
            "second_child_distinct_from_fixed_child": True,
            "requires_second_difference_zero": True,
            "requires_next_sum_and_difference_nonzero": True,
            "diagonal_forbidden_shells_squared": [6],
            "exhaustive": True,
        },
        "catalog": {
            "canonical_wave_count": len(waves),
            "wavevectors": [list(wave) for wave in waves],
            "polarizations_by_wave": {
                ",".join(str(value) for value in wave): [
                    list(vector) for vector in polarizations[wave]
                ]
                for wave in waves
            },
            "total_primitive_polarization_direction_count": (
                total_polarization_directions
            ),
        },
        "counts": {
            "ordered_second_parent_wave_pairs_examined": len(all_ordered_wave_pairs),
            "scope_eligible_ordered_wave_pairs": len(scope_wave_pairs),
            "polarization_pair_tests": polarization_pair_tests,
            "strict_hit_orientation_count": len(hit_orientations),
            "unique_strict_hit_count": len(strict_payloads),
            "unique_cross_talk_partial_count": len(unique_partials),
            "saved_top_partial_count": len(partial_payloads),
            "strict_hits_without_lower_shell_return": strict_without_low_return,
            "saved_structural_record_count": len(structural_records),
        },
        "exact_arithmetic_stats_over_full_fast_search": fraction_stats,
        "rejection_histogram": rejection_histogram,
        "all_structural_orientation_records": structural_records,
        "fixed_known_relay": _fixed_relay_payload(),
        "strict_hits": strict_payloads,
        "top_partial_candidates": partial_payloads,
        "parametric_rational_family": (
            _parametric_family_payload() if found else None
        ),
        "continuous_polarization_classification": (
            _continuous_polarization_classification_payload() if found else None
        ),
        "relay_quality_screen": {
            "scope": "all four unintended parent pairs for every strict hit",
            "requires_no_output_below_child_shell_squared_6": True,
            "passing_orientation_count": strict_without_low_return,
            "exhaustive_in_recorded_alphabet": True,
            "classification": (
                "FINITE-SCOPE NEGATIVE: EVERY STRICT HIT RETURNS TO LOWER SHELLS"
                if found and strict_without_low_return == 0
                else "quality screen not uniformly negative"
            ),
        },
        "verdict": (
            "A strict finite carrier extension exists in this bounded alphabet; "
            "the earlier coordinate-axis-polarization negative does not extend "
            "to primitive integer polarizations. The displayed simple aligned "
            "binary recursion is exactly rejected; broader closure remains open."
            if found
            else
            "No strict extension exists in this exact finite alphabet; this has "
            "no force outside the recorded bounds."
        ),
    }
