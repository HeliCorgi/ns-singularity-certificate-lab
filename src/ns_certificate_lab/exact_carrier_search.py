r"""Deterministic exact search for a finite two-carrier relay gadget.

This module searches a deliberately tiny, completely enumerated class of real
Fourier modes on the normalised ``(2*pi)``-periodic three-torus.  It is a
finite-algebra discovery tool, not a Navier--Stokes orbit or a cascade proof.

The default search space is fixed as follows.

* canonical nonzero wavevectors with coordinates in ``{-1, 0, 1}`` and
  squared length in ``{1, 2}`` (nine wavevectors);
* unit coordinate polarizations perpendicular to the wavevector;
* unit-amplitude cosine and sine phases;
* unordered pairs with distinct wavevectors.

For a pair of modes ``A, B`` the exact mixed Leray output is

``P((A . grad) B + (B . grad) A)``.

Its sum-frequency component is populated with the negative of that output,
so the signed parent-to-child transfer is its exact squared ``L2`` norm.  Two
such relays with disjoint parent wavevectors are retained only when all four
parents share one Fourier radius, the two children share a second radius, and
both first-generation difference outputs vanish.  The two children must then
have nonzero sum- and difference-frequency Leray outputs.

Such a pair is only an algebraic *partial candidate*.  It is screened
separately against the diagonal parent cross interactions ``A1+B2`` and
``A2+B1`` on the target child shell.  The default run returns the first partial
candidate in deterministic order, explicitly marking that search
nonexhaustive.  Separately, it fixes the known relay below and exhausts every
compatible orientation of a second relay under the strict cross-talk screen.
Every coefficient and energy identity is evaluated by :class:`TrigVector`
and :class:`fractions.Fraction`.

The known relay ``(1,1,0), e3*sin`` with ``(1,0,1), e2*cos`` is placed first
in the deterministic enumeration.  A strict failure is an exhaustive
*negative certificate only for extensions of that fixed relay in the stated
finite pool*; it says nothing about larger or different carrier families.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from typing import Iterable, Literal

from .fourier_torus import TrigVector, advection, leray
from .torus_chain import l2_inner

__all__ = [
    "DEFAULT_WAVEVECTORS",
    "CarrierGadget",
    "CarrierSearchCertificate",
    "CrossTalkDiagnostic",
    "ExactInteraction",
    "ModeSpec",
    "RelayPairWitness",
    "StrictOrientationRecord",
    "StrictCrossTalkCertificate",
    "search_exact_carrier_gadget",
]


Wavevector = tuple[int, int, int]
IntegerVector = tuple[int, int, int]
Phase = Literal["cos", "sin"]

SEARCH_SCHEMA = "ns-certificate-lab/exact-carrier-search/v3"
PARTIAL_STATUS = (
    "PARTIAL / REJECTED CROSS-TALK / SIMPLE RECURSION REJECTED / "
    "NOT A PDE OR CASCADE"
)
STRICT_FOUND_STATUS = "STRICT FINITE GADGET / NOT A PDE OR CASCADE"
NEGATIVE_STATUS = "NEGATIVE IN THE STATED FINITE SEARCH ONLY"
STRICT_NEGATIVE_STATUS = (
    "NEGATIVE STRICT CROSS-TALK EXTENSION OF KNOWN RELAY / FINITE SEARCH ONLY"
)
_ZERO_VECTOR: IntegerVector = (0, 0, 0)
_COORDINATE_AXES: tuple[IntegerVector, ...] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)
_PHASES: tuple[Phase, ...] = ("cos", "sin")
_MAX_CUSTOM_WAVEVECTORS = 16


def _norm_squared(wave: Wavevector) -> int:
    return sum(component * component for component in wave)


def _canonical_wave(wave: Wavevector) -> tuple[Wavevector, int]:
    """Return the canonical half-lattice representative and sine sign."""

    for component in wave:
        if component > 0:
            return wave, 1
        if component < 0:
            return tuple(-value for value in wave), -1  # type: ignore[return-value]
    raise ValueError("the zero wavevector has no carrier representative")


def _default_wavevectors() -> tuple[Wavevector, ...]:
    waves: list[Wavevector] = []
    for wave in product(range(-1, 2), repeat=3):
        candidate: Wavevector = tuple(wave)  # type: ignore[assignment]
        if _norm_squared(candidate) not in (1, 2):
            continue
        _, sign = _canonical_wave(candidate)
        if sign > 0:
            waves.append(candidate)
    return tuple(sorted(waves))


DEFAULT_WAVEVECTORS: tuple[Wavevector, ...] = _default_wavevectors()


def _field_table(field: TrigVector) -> dict[str, dict[str, list[str]]]:
    output: dict[str, dict[str, list[str]]] = {}
    for wave, (cosine, sine) in sorted(field.cleaned().coefficient_table().items()):
        output[",".join(str(component) for component in wave)] = {
            "cosine": [str(Fraction(value)) for value in cosine],
            "sine": [str(Fraction(value)) for value in sine],
        }
    return output


def _is_zero(field: TrigVector) -> bool:
    return not field.cleaned().coefficient_table()


def _require_divergence_free(field: TrigVector, *, label: str) -> None:
    if field.divergence().cleaned().terms:
        raise AssertionError(f"{label} is not exactly divergence free")


@dataclass(frozen=True)
class ModeSpec:
    """One unit-amplitude real plane wave in the finite search alphabet."""

    wave: Wavevector
    polarization: IntegerVector
    phase: Phase

    def __post_init__(self) -> None:
        if len(self.wave) != 3 or any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in self.wave
        ):
            raise ValueError("wave must be an integer three-vector")
        if self.wave == _ZERO_VECTOR or _canonical_wave(self.wave)[0] != self.wave:
            raise ValueError("wave must be a nonzero canonical wavevector")
        if len(self.polarization) != 3 or any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in self.polarization
        ):
            raise ValueError("polarization must be an integer three-vector")
        if self.polarization == _ZERO_VECTOR:
            raise ValueError("polarization must be nonzero")
        if sum(
            wave_component * polarization_component
            for wave_component, polarization_component in zip(
                self.wave, self.polarization
            )
        ):
            raise ValueError("polarization must be perpendicular to its wavevector")
        if self.phase not in _PHASES:
            raise ValueError("phase must be 'cos' or 'sin'")

    def field(self) -> TrigVector:
        cosine = self.polarization if self.phase == "cos" else _ZERO_VECTOR
        sine = self.polarization if self.phase == "sin" else _ZERO_VECTOR
        field = TrigVector.from_modes([(self.wave, cosine, sine)]).cleaned()
        _require_divergence_free(field, label="search mode")
        return field

    def as_dict(self) -> dict[str, object]:
        return {
            "wave": list(self.wave),
            "polarization": list(self.polarization),
            "phase": self.phase,
        }


@dataclass(frozen=True)
class ExactInteraction:
    """Exact sum/difference diagnostics for two single Fourier modes."""

    left_wave: Wavevector
    right_wave: Wavevector
    left: TrigVector
    right: TrigVector
    mixed_output: TrigVector
    sum_wave: Wavevector
    difference_wave: Wavevector
    sum_output: TrigVector
    difference_output: TrigVector
    other_output: TrigVector
    child: TrigVector
    signed_parent_to_child_flux: Fraction
    parent_flux_in_populated_triad: Fraction
    child_flux_in_populated_triad: Fraction
    populated_triad_energy_defect: Fraction

    def parent_field(self) -> TrigVector:
        return (self.left + self.right).cleaned()

    def populated_triad(self) -> TrigVector:
        return (self.left + self.right + self.child).cleaned()

    def as_dict(self) -> dict[str, object]:
        return {
            "left_wave": list(self.left_wave),
            "right_wave": list(self.right_wave),
            "sum_wave": list(self.sum_wave),
            "difference_wave": list(self.difference_wave),
            "mixed_output": _field_table(self.mixed_output),
            "sum_output": _field_table(self.sum_output),
            "difference_output": _field_table(self.difference_output),
            "other_output": _field_table(self.other_output),
            "child": _field_table(self.child),
            "signed_parent_to_child_flux": str(
                self.signed_parent_to_child_flux
            ),
            "parent_flux_in_populated_triad": str(
                self.parent_flux_in_populated_triad
            ),
            "child_flux_in_populated_triad": str(
                self.child_flux_in_populated_triad
            ),
            "populated_triad_energy_defect": str(
                self.populated_triad_energy_defect
            ),
        }


@dataclass(frozen=True)
class RelayPairWitness:
    """A searched parent pair together with its exact child interaction."""

    left_spec: ModeSpec
    right_spec: ModeSpec
    interaction: ExactInteraction

    def parent_waves(self) -> frozenset[Wavevector]:
        return frozenset((self.left_spec.wave, self.right_spec.wave))

    def carrier_shells(self) -> tuple[int, int] | None:
        """Return common parent and child squared radii, when they exist."""

        left_radius = _norm_squared(self.left_spec.wave)
        right_radius = _norm_squared(self.right_spec.wave)
        if left_radius != right_radius:
            return None
        return left_radius, _norm_squared(self.interaction.sum_wave)

    def as_dict(self) -> dict[str, object]:
        shells = self.carrier_shells()
        return {
            "left": self.left_spec.as_dict(),
            "right": self.right_spec.as_dict(),
            "carrier_shells_squared": (
                None if shells is None else list(shells)
            ),
            "interaction": self.interaction.as_dict(),
        }


@dataclass(frozen=True)
class CrossTalkDiagnostic:
    """Exact diagonal parent cross-talk on the intended child shell."""

    target_shell_squared: int
    a1_b2_output: TrigVector
    a2_b1_output: TrigVector
    a1_b2_target_shell_output: TrigVector
    a2_b1_target_shell_output: TrigVector
    passes: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "condition": (
                "A1+B2 and A2+B1 must vanish on the target child shell"
            ),
            "target_shell_squared": self.target_shell_squared,
            "a1_b2_output": _field_table(self.a1_b2_output),
            "a2_b1_output": _field_table(self.a2_b1_output),
            "a1_b2_target_shell_output": _field_table(
                self.a1_b2_target_shell_output
            ),
            "a2_b1_target_shell_output": _field_table(
                self.a2_b1_target_shell_output
            ),
            "passes": self.passes,
        }


@dataclass(frozen=True)
class CarrierGadget:
    """Two exact relays, with an explicit strict-screen classification."""

    classification: str
    relay_one: RelayPairWitness
    relay_two: RelayPairWitness
    next_interaction: ExactInteraction
    common_parent_shell_squared: int
    common_child_shell_squared: int
    diagonal_cross_talk: CrossTalkDiagnostic
    relay_cross_output: TrigVector
    relay_cross_output_l2_squared: Fraction
    intended_next_output_l2_squared: Fraction
    intended_fraction_of_relay_cross_power: Fraction
    grandchild_mutual_output: TrigVector
    grandchild_mutual_output_l2_squared: Fraction
    simple_binary_recursion_rejected: bool
    full_field_flux_by_wave: tuple[tuple[Wavevector, Fraction], ...]
    full_field_energy_defect: Fraction

    def as_dict(self) -> dict[str, object]:
        return {
            "classification": self.classification,
            "relay_one": self.relay_one.as_dict(),
            "relay_two": self.relay_two.as_dict(),
            "next_interaction": self.next_interaction.as_dict(),
            "common_parent_shell_squared": self.common_parent_shell_squared,
            "common_child_shell_squared": self.common_child_shell_squared,
            "diagonal_cross_talk": self.diagonal_cross_talk.as_dict(),
            "relay_cross_output": _field_table(self.relay_cross_output),
            "relay_cross_output_l2_squared": str(
                self.relay_cross_output_l2_squared
            ),
            "intended_next_output_l2_squared": str(
                self.intended_next_output_l2_squared
            ),
            "intended_fraction_of_relay_cross_power": str(
                self.intended_fraction_of_relay_cross_power
            ),
            "grandchild_mutual_output": _field_table(
                self.grandchild_mutual_output
            ),
            "grandchild_mutual_output_l2_squared": str(
                self.grandchild_mutual_output_l2_squared
            ),
            "simple_binary_recursion_rejected": (
                self.simple_binary_recursion_rejected
            ),
            "full_field_flux_by_wave": {
                ",".join(str(component) for component in wave): str(flux)
                for wave, flux in self.full_field_flux_by_wave
            },
            "full_field_energy_defect": str(self.full_field_energy_defect),
        }


@dataclass(frozen=True)
class StrictOrientationRecord:
    """One fully disclosed orientation in the fixed-relay strict search."""

    orientation_index: int
    relay_two: RelayPairWitness
    next_interaction: ExactInteraction | None
    diagonal_cross_talk: CrossTalkDiagnostic
    required_next_interaction_passes: bool
    strict_passes: bool
    rejection_reasons: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "orientation_index": self.orientation_index,
            "relay_two": self.relay_two.as_dict(),
            "next_interaction": (
                None
                if self.next_interaction is None
                else self.next_interaction.as_dict()
            ),
            "diagonal_cross_talk": self.diagonal_cross_talk.as_dict(),
            "required_next_interaction_passes": (
                self.required_next_interaction_passes
            ),
            "strict_passes": self.strict_passes,
            "rejection_reasons": list(self.rejection_reasons),
        }


@dataclass(frozen=True)
class StrictCrossTalkCertificate:
    """Result of the exhaustive strict search with the known relay fixed."""

    status: str
    found: bool
    exhaustive: bool
    fixed_known_relay_available: bool
    eligible_second_relay_orientations_tested: int
    orientation_records: tuple[StrictOrientationRecord, ...]
    gadget: CarrierGadget | None
    reason: str

    def as_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "found": self.found,
            "exhaustive": self.exhaustive,
            "scope": "extensions of the fixed known relay in this finite pool",
            "fixed_known_relay_available": self.fixed_known_relay_available,
            "eligible_second_relay_orientations_tested": (
                self.eligible_second_relay_orientations_tested
            ),
            "orientation_records": [
                record.as_dict() for record in self.orientation_records
            ],
            "gadget": None if self.gadget is None else self.gadget.as_dict(),
            "reason": self.reason,
        }


@dataclass(frozen=True)
class CarrierSearchCertificate:
    """Separated partial-candidate and strict finite-search results."""

    schema: str
    classification: str
    status: str
    wavevectors: tuple[Wavevector, ...]
    mode_spec_count: int
    unordered_mode_pair_tests: int
    viable_relay_pair_count: int
    partial_eligible_relay_pair_tests: int
    partial_candidate_found: bool
    partial_search_exhaustive: bool
    search_rules: tuple[str, ...]
    partial_gadget: CarrierGadget | None
    strict_search: StrictCrossTalkCertificate
    reason: str

    def as_dict(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "classification": self.classification,
            "status": self.status,
            "partial_candidate_found": self.partial_candidate_found,
            "partial_search_exhaustive": self.partial_search_exhaustive,
            "search_space": {
                "wavevectors": [list(wave) for wave in self.wavevectors],
                "mode_spec_count": self.mode_spec_count,
                "polarizations": "unit coordinate axes perpendicular to k",
                "phases": list(_PHASES),
                "amplitude": "1",
                "rules": list(self.search_rules),
            },
            "counts": {
                "unordered_mode_pair_tests": self.unordered_mode_pair_tests,
                "viable_relay_pair_count": self.viable_relay_pair_count,
                "partial_eligible_relay_pair_tests": (
                    self.partial_eligible_relay_pair_tests
                ),
            },
            "partial_gadget": (
                None
                if self.partial_gadget is None
                else self.partial_gadget.as_dict()
            ),
            "strict_search": self.strict_search.as_dict(),
            "reason": self.reason,
        }


_KNOWN_LEFT = ModeSpec((1, 1, 0), (0, 0, 1), "sin")
_KNOWN_RIGHT = ModeSpec((1, 0, 1), (0, 1, 0), "cos")


def _mode_specs(wavevectors: tuple[Wavevector, ...]) -> tuple[ModeSpec, ...]:
    allowed = set(wavevectors)
    modes: list[ModeSpec] = []
    for known in (_KNOWN_LEFT, _KNOWN_RIGHT):
        if known.wave in allowed:
            modes.append(known)
    for wave in wavevectors:
        for polarization in _COORDINATE_AXES:
            if sum(a * b for a, b in zip(wave, polarization)):
                continue
            for phase in _PHASES:
                candidate = ModeSpec(wave, polarization, phase)
                if candidate not in modes:
                    modes.append(candidate)
    return tuple(modes)


def _mode_band(field: TrigVector, wave: Wavevector) -> TrigVector:
    return field.restrict(lambda candidate: candidate == wave).cleaned()


def _build_interaction(
    left: TrigVector,
    right: TrigVector,
    *,
    left_wave: Wavevector,
    right_wave: Wavevector,
) -> ExactInteraction | None:
    _require_divergence_free(left, label="left interaction mode")
    _require_divergence_free(right, label="right interaction mode")
    sum_wave, _ = _canonical_wave(
        tuple(a + b for a, b in zip(left_wave, right_wave))  # type: ignore[arg-type]
    )
    difference_wave, _ = _canonical_wave(
        tuple(a - b for a, b in zip(left_wave, right_wave))  # type: ignore[arg-type]
    )
    mixed = leray(advection(left, right) + advection(right, left)).cleaned()
    sum_output = _mode_band(mixed, sum_wave)
    if _is_zero(sum_output):
        return None
    difference_output = _mode_band(mixed, difference_wave)
    other_output = mixed.restrict(
        lambda wave: wave not in (sum_wave, difference_wave)
    ).cleaned()
    child = (-sum_output).cleaned()
    _require_divergence_free(mixed, label="mixed Leray output")
    _require_divergence_free(child, label="aligned child")

    signed_injection = -l2_inner(child, mixed)
    if signed_injection <= 0:
        raise AssertionError("an aligned nonzero child did not have positive transfer")

    parent = (left + right).cleaned()
    triad = (parent + child).cleaned()
    nonlinear = leray(advection(triad, triad)).cleaned()
    parent_flux = -l2_inner(parent, nonlinear)
    child_flux = -l2_inner(child, nonlinear)
    energy_defect = l2_inner(triad, nonlinear)
    if energy_defect != 0 or parent_flux + child_flux != 0:
        raise AssertionError("the exact populated-triad energy identity failed")
    if child_flux != signed_injection:
        raise AssertionError("populating the aligned child changed its own flux")

    return ExactInteraction(
        left_wave=left_wave,
        right_wave=right_wave,
        left=left,
        right=right,
        mixed_output=mixed,
        sum_wave=sum_wave,
        difference_wave=difference_wave,
        sum_output=sum_output,
        difference_output=difference_output,
        other_output=other_output,
        child=child,
        signed_parent_to_child_flux=signed_injection,
        parent_flux_in_populated_triad=parent_flux,
        child_flux_in_populated_triad=child_flux,
        populated_triad_energy_defect=energy_defect,
    )


def _mixed_output(left: TrigVector, right: TrigVector) -> TrigVector:
    """Return the exact symmetric mixed Leray interaction."""

    output = leray(advection(left, right) + advection(right, left)).cleaned()
    _require_divergence_free(output, label="mixed cross interaction")
    return output


def _shell_band(field: TrigVector, shell_squared: int) -> TrigVector:
    return field.restrict(
        lambda wave: _norm_squared(wave) == shell_squared
    ).cleaned()


def _cross_talk_diagnostic(
    relay_one: RelayPairWitness,
    relay_two: RelayPairWitness,
    target_shell_squared: int,
) -> CrossTalkDiagnostic:
    """Test the ordered diagonal interactions A1+B2 and A2+B1."""

    a1_b2 = _mixed_output(
        relay_one.interaction.left,
        relay_two.interaction.right,
    )
    a2_b1 = _mixed_output(
        relay_two.interaction.left,
        relay_one.interaction.right,
    )
    a1_b2_target = _shell_band(a1_b2, target_shell_squared)
    a2_b1_target = _shell_band(a2_b1, target_shell_squared)
    return CrossTalkDiagnostic(
        target_shell_squared=target_shell_squared,
        a1_b2_output=a1_b2,
        a2_b1_output=a2_b1,
        a1_b2_target_shell_output=a1_b2_target,
        a2_b1_target_shell_output=a2_b1_target,
        passes=_is_zero(a1_b2_target) and _is_zero(a2_b1_target),
    )


def _reverse_relay(relay: RelayPairWitness) -> RelayPairWitness:
    interaction = _build_interaction(
        relay.interaction.right,
        relay.interaction.left,
        left_wave=relay.right_spec.wave,
        right_wave=relay.left_spec.wave,
    )
    if interaction is None:
        raise AssertionError("reversing a viable relay removed its sum output")
    return RelayPairWitness(
        left_spec=relay.right_spec,
        right_spec=relay.left_spec,
        interaction=interaction,
    )


def _known_relay(
    relay_candidates: Iterable[RelayPairWitness],
) -> RelayPairWitness | None:
    for relay in relay_candidates:
        if (
            relay.left_spec == _KNOWN_LEFT
            and relay.right_spec == _KNOWN_RIGHT
        ):
            return relay
        if (
            relay.left_spec == _KNOWN_RIGHT
            and relay.right_spec == _KNOWN_LEFT
        ):
            return _reverse_relay(relay)
    return None


def _compatible_carrier_shells(
    relay_one: RelayPairWitness,
    relay_two: RelayPairWitness,
) -> tuple[int, int] | None:
    if relay_one.parent_waves() & relay_two.parent_waves():
        return None
    if relay_one.interaction.sum_wave == relay_two.interaction.sum_wave:
        return None
    carrier_shells = relay_one.carrier_shells()
    if carrier_shells is None or relay_two.carrier_shells() != carrier_shells:
        return None
    if not _is_zero(relay_one.interaction.difference_output):
        return None
    if not _is_zero(relay_two.interaction.difference_output):
        return None
    return carrier_shells


def _next_interaction(
    relay_one: RelayPairWitness,
    relay_two: RelayPairWitness,
) -> ExactInteraction | None:
    interaction = _build_interaction(
        relay_one.interaction.child,
        relay_two.interaction.child,
        left_wave=relay_one.interaction.sum_wave,
        right_wave=relay_two.interaction.sum_wave,
    )
    if interaction is None or _is_zero(interaction.difference_output):
        return None
    return interaction


def _build_gadget(
    relay_one: RelayPairWitness,
    relay_two: RelayPairWitness,
    next_interaction: ExactInteraction,
    carrier_shells: tuple[int, int],
    diagonal_cross_talk: CrossTalkDiagnostic,
) -> CarrierGadget:
    first_triad = relay_one.interaction.populated_triad()
    second_triad = relay_two.interaction.populated_triad()
    cross_output = _mixed_output(first_triad, second_triad)
    full_field = (first_triad + second_triad + next_interaction.child).cleaned()
    _require_divergence_free(full_field, label="full carrier gadget")
    full_nonlinear = leray(advection(full_field, full_field)).cleaned()
    full_energy_defect = l2_inner(full_field, full_nonlinear)
    if full_energy_defect != 0:
        raise AssertionError("the full exact carrier energy identity failed")
    full_field_flux_by_wave = tuple(
        (
            wave,
            -l2_inner(_mode_band(full_field, wave), full_nonlinear),
        )
        for wave in sorted(full_field.coefficient_table())
    )
    if sum((flux for _, flux in full_field_flux_by_wave), Fraction(0)) != 0:
        raise AssertionError("the exact modal flux decomposition failed")

    difference_child = (-next_interaction.difference_output).cleaned()
    grandchild_mutual_output = _mixed_output(
        next_interaction.child,
        difference_child,
    )
    grandchild_mutual_norm = Fraction(
        grandchild_mutual_output.sobolev_sq(0)
    )
    relay_cross_norm = Fraction(cross_output.sobolev_sq(0))
    intended_next_norm = Fraction(
        (
            next_interaction.sum_output
            + next_interaction.difference_output
        ).sobolev_sq(0)
    )
    intended_fraction = (
        intended_next_norm / relay_cross_norm
        if relay_cross_norm
        else Fraction(0)
    )
    classification = (
        "strict_cross_talk_pass"
        if diagonal_cross_talk.passes
        else "partial_rejected_cross_talk"
    )
    return CarrierGadget(
        classification=classification,
        relay_one=relay_one,
        relay_two=relay_two,
        next_interaction=next_interaction,
        common_parent_shell_squared=carrier_shells[0],
        common_child_shell_squared=carrier_shells[1],
        diagonal_cross_talk=diagonal_cross_talk,
        relay_cross_output=cross_output,
        relay_cross_output_l2_squared=relay_cross_norm,
        intended_next_output_l2_squared=intended_next_norm,
        intended_fraction_of_relay_cross_power=intended_fraction,
        grandchild_mutual_output=grandchild_mutual_output,
        grandchild_mutual_output_l2_squared=grandchild_mutual_norm,
        simple_binary_recursion_rejected=grandchild_mutual_norm == 0,
        full_field_flux_by_wave=full_field_flux_by_wave,
        full_field_energy_defect=full_energy_defect,
    )


def _strict_known_relay_search(
    relay_candidates: tuple[RelayPairWitness, ...],
) -> StrictCrossTalkCertificate:
    known = _known_relay(relay_candidates)
    if known is None:
        return StrictCrossTalkCertificate(
            status=STRICT_NEGATIVE_STATUS,
            found=False,
            exhaustive=True,
            fixed_known_relay_available=False,
            eligible_second_relay_orientations_tested=0,
            orientation_records=(),
            gadget=None,
            reason=(
                "the fixed known relay is unavailable in the supplied pool; "
                "therefore no extension of it exists in this finite scope"
            ),
        )

    tested = 0
    records: list[StrictOrientationRecord] = []
    first_strict_gadget: CarrierGadget | None = None
    for relay in relay_candidates:
        carrier_shells = _compatible_carrier_shells(known, relay)
        if carrier_shells is None:
            continue
        next_interaction = _next_interaction(known, relay)
        for oriented_relay in (relay, _reverse_relay(relay)):
            tested += 1
            cross_talk = _cross_talk_diagnostic(
                known,
                oriented_relay,
                carrier_shells[1],
            )
            rejection_reasons: list[str] = []
            if next_interaction is None:
                rejection_reasons.append(
                    "next_child_interaction_lacks_required_nonzero_sum_or_difference"
                )
            if not cross_talk.passes:
                rejection_reasons.append(
                    "diagonal_parent_cross_talk_hits_target_child_shell"
                )
            strict_passes = not rejection_reasons
            records.append(
                StrictOrientationRecord(
                    orientation_index=tested,
                    relay_two=oriented_relay,
                    next_interaction=next_interaction,
                    diagonal_cross_talk=cross_talk,
                    required_next_interaction_passes=next_interaction is not None,
                    strict_passes=strict_passes,
                    rejection_reasons=tuple(rejection_reasons),
                )
            )
            if strict_passes and first_strict_gadget is None:
                if next_interaction is None:
                    raise AssertionError("strict pass lacks a next interaction")
                first_strict_gadget = _build_gadget(
                    known,
                    oriented_relay,
                    next_interaction,
                    carrier_shells,
                    cross_talk,
                )

    if first_strict_gadget is not None:
        return StrictCrossTalkCertificate(
            status=STRICT_FOUND_STATUS,
            found=True,
            exhaustive=True,
            fixed_known_relay_available=True,
            eligible_second_relay_orientations_tested=tested,
            orientation_records=tuple(records),
            gadget=first_strict_gadget,
            reason=(
                "at least one strict extension of the fixed known relay was "
                "found; the entire stated finite scope was still exhausted"
            ),
        )
    return StrictCrossTalkCertificate(
        status=STRICT_NEGATIVE_STATUS,
        found=False,
        exhaustive=True,
        fixed_known_relay_available=True,
        eligible_second_relay_orientations_tested=tested,
        orientation_records=tuple(records),
        gadget=None,
        reason=(
            f"all {tested} compatible oriented second relays extending the "
            "fixed known relay were tested; none kept both diagonal cross "
            "interactions off the target child shell while producing nonzero "
            "child sum and difference outputs.  This is not a no-go result "
            "outside that fixed finite scope"
        ),
    )


def _validate_wavevectors(wavevectors: Iterable[Wavevector]) -> tuple[Wavevector, ...]:
    values = tuple(wavevectors)
    if len(values) > _MAX_CUSTOM_WAVEVECTORS:
        raise ValueError(
            f"at most {_MAX_CUSTOM_WAVEVECTORS} wavevectors may be searched"
        )
    if len(set(values)) != len(values):
        raise ValueError("wavevector pool must not contain duplicates")
    for wave in values:
        if len(wave) != 3 or any(
            isinstance(component, bool) or not isinstance(component, int)
            for component in wave
        ):
            raise ValueError("every wavevector must be an integer three-vector")
        if wave == _ZERO_VECTOR or _canonical_wave(wave)[0] != wave:
            raise ValueError("every wavevector must be nonzero and canonical")
    return tuple(sorted(values))


def search_exact_carrier_gadget(
    wavevectors: Iterable[Wavevector] | None = None,
) -> CarrierSearchCertificate:
    """Return separate first-hit and strict results for the finite alphabet.

    Passing ``None`` selects :data:`DEFAULT_WAVEVECTORS`.  A custom pool is
    useful for exact negative controls and is capped at sixteen canonical
    wavevectors so that this function remains a small CPU-only search.
    """

    waves = (
        DEFAULT_WAVEVECTORS
        if wavevectors is None
        else _validate_wavevectors(wavevectors)
    )
    specs = _mode_specs(waves)
    relay_candidates: list[RelayPairWitness] = []
    mode_pair_tests = 0
    for left_spec, right_spec in combinations(specs, 2):
        mode_pair_tests += 1
        if left_spec.wave == right_spec.wave:
            continue
        interaction = _build_interaction(
            left_spec.field(),
            right_spec.field(),
            left_wave=left_spec.wave,
            right_wave=right_spec.wave,
        )
        if interaction is None:
            continue
        relay_candidates.append(
            RelayPairWitness(
                left_spec=left_spec,
                right_spec=right_spec,
                interaction=interaction,
            )
        )

    relay_candidates_tuple = tuple(relay_candidates)
    partial_pair_tests = 0
    partial_gadget: CarrierGadget | None = None
    for relay_one, relay_two in combinations(relay_candidates, 2):
        carrier_shells = _compatible_carrier_shells(relay_one, relay_two)
        if carrier_shells is None:
            continue
        partial_pair_tests += 1
        next_interaction = _next_interaction(relay_one, relay_two)
        if next_interaction is None:
            continue
        cross_talk = _cross_talk_diagnostic(
            relay_one,
            relay_two,
            carrier_shells[1],
        )
        partial_gadget = _build_gadget(
            relay_one,
            relay_two,
            next_interaction,
            carrier_shells,
            cross_talk,
        )
        break

    strict_search = _strict_known_relay_search(relay_candidates_tuple)
    partial_found = partial_gadget is not None
    partial_exhaustive = not partial_found
    if strict_search.found:
        classification = "strict_gadget_found"
        status = STRICT_FOUND_STATUS
        reason = (
            "the strict fixed-known-relay search found a finite exact witness; "
            "this is not a Navier--Stokes orbit, invariant cone, or cascade"
        )
    elif partial_gadget is not None:
        if partial_gadget.diagonal_cross_talk.passes:
            classification = "partial_candidate_passes_diagonal_screen"
            status = (
                "PARTIAL FIRST-HIT PASSES DIAGONAL SCREEN / "
                "FIXED-KNOWN STRICT SEARCH NEGATIVE"
            )
            reason = (
                "the nonexhaustive first-hit search found a candidate passing "
                "the diagonal screen, but it is outside the separately scoped "
                "strict extension result and is not a cascade"
            )
        else:
            classification = "partial_rejected_cross_talk"
            status = PARTIAL_STATUS
            reason = (
                "the first same-shell algebraic candidate fails the exact "
                "diagonal target-shell cross-talk condition; the first-hit "
                "partial search is explicitly nonexhaustive, while the "
                "fixed-known-relay strict search is exhaustively negative.  "
                "Its intended next sum/difference outputs contain only "
                f"{partial_gadget.intended_fraction_of_relay_cross_power} "
                "of the exact relay-cross L2 power, and its two aligned "
                "grandchildren interact exactly to zero, so the simple "
                "binary recursion is REJECTED"
            )
    else:
        classification = "negative_in_finite_search"
        status = NEGATIVE_STATUS
        reason = (
            "no partial carrier pair exists in the explicitly listed finite "
            "search space; this is not a no-go result outside that space"
        )

    return CarrierSearchCertificate(
        schema=SEARCH_SCHEMA,
        classification=classification,
        status=status,
        wavevectors=waves,
        mode_spec_count=len(specs),
        unordered_mode_pair_tests=mode_pair_tests,
        viable_relay_pair_count=len(relay_candidates),
        partial_eligible_relay_pair_tests=partial_pair_tests,
        partial_candidate_found=partial_found,
        partial_search_exhaustive=partial_exhaustive,
        search_rules=(
            "all mode pairs in the supplied finite alphabet were enumerated",
            "the partial search stops at its first candidate, if one exists",
            "partial relays must share parent and child squared radii",
            "both first-generation difference outputs must vanish exactly",
            "a partial candidate requires nonzero child sum and difference outputs",
            "strict cross-talk means A1+B2 and A2+B1 avoid the child shell",
            "the strict search fixes and includes the known exact relay",
            "both orientations of every compatible strict second relay are tested",
            "strict negative scope is only the fixed relay and finite alphabet",
        ),
        partial_gadget=partial_gadget,
        strict_search=strict_search,
        reason=reason,
    )
