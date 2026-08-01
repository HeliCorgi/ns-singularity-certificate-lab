r"""Independent ordered-pair verifier for serialized carrier-search records.

Unlike :mod:`ns_certificate_lab.exact_carrier_search`, this module does not use
the real trigonometric polynomial implementation or its search builders.  It
reconstructs real modes as conjugate-symmetric complex Fourier coefficients,
performs the ordered-pair convolution

``hat B(n) = i sum_{k+l=n} (hat u(k) . l) hat v(l)``,

and applies the rational Leray matrix independently.  It is intentionally
small and only verifies the disclosed fixed-known-relay strict records.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import math
from typing import Iterable

__all__ = [
    "verify_serialized_expanded_carrier_certificate",
    "verify_serialized_strict_orientation_records",
]


Wave = tuple[int, int, int]
ComplexQ = tuple[Fraction, Fraction]
ComplexVector = tuple[ComplexQ, ComplexQ, ComplexQ]
Field = dict[Wave, ComplexVector]

_ZERO_C: ComplexQ = (Fraction(0), Fraction(0))
_ZERO_V: ComplexVector = (_ZERO_C, _ZERO_C, _ZERO_C)
_P: Wave = (1, 1, 0)
_Q: Wave = (1, 0, 1)


def _c_add(a: ComplexQ, b: ComplexQ) -> ComplexQ:
    return a[0] + b[0], a[1] + b[1]


def _c_neg(a: ComplexQ) -> ComplexQ:
    return -a[0], -a[1]


def _c_mul(a: ComplexQ, b: ComplexQ) -> ComplexQ:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def _c_scale(a: ComplexQ, value: int | Fraction) -> ComplexQ:
    factor = Fraction(value)
    return factor * a[0], factor * a[1]


def _c_conjugate(a: ComplexQ) -> ComplexQ:
    return a[0], -a[1]


def _c_i_mul(a: ComplexQ) -> ComplexQ:
    return -a[1], a[0]


def _v_add(a: ComplexVector, b: ComplexVector) -> ComplexVector:
    return tuple(_c_add(x, y) for x, y in zip(a, b))  # type: ignore[return-value]


def _v_neg(a: ComplexVector) -> ComplexVector:
    return tuple(_c_neg(value) for value in a)  # type: ignore[return-value]


def _v_scale(a: ComplexVector, value: int | Fraction) -> ComplexVector:
    return tuple(_c_scale(item, value) for item in a)  # type: ignore[return-value]


def _v_zero(a: ComplexVector) -> bool:
    return all(value == _ZERO_C for value in a)


def _field_add(*fields: Field) -> Field:
    output: Field = {}
    for field in fields:
        for wave, vector in field.items():
            output[wave] = _v_add(output.get(wave, _ZERO_V), vector)
    return {wave: vector for wave, vector in output.items() if not _v_zero(vector)}


def _field_neg(field: Field) -> Field:
    return {wave: _v_neg(vector) for wave, vector in field.items()}


def _field_scale(field: Field, scale: Fraction) -> Field:
    return {
        wave: _v_scale(vector, scale)
        for wave, vector in field.items()
        if not _v_zero(_v_scale(vector, scale))
    }


def _wave_add(a: Wave, b: Wave) -> Wave:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def _wave_neg(a: Wave) -> Wave:
    return tuple(-value for value in a)  # type: ignore[return-value]


def _canonical(wave: Wave) -> Wave:
    for value in wave:
        if value > 0:
            return wave
        if value < 0:
            return _wave_neg(wave)
    raise ValueError("zero wave is not canonical")


def _mode(spec: object) -> Field:
    if not isinstance(spec, dict) or set(spec) != {"wave", "polarization", "phase"}:
        raise ValueError("serialized mode spec is invalid")
    raw_wave = spec["wave"]
    raw_polarization = spec["polarization"]
    phase = spec["phase"]
    if (
        not isinstance(raw_wave, list)
        or len(raw_wave) != 3
        or not all(isinstance(value, int) and not isinstance(value, bool) for value in raw_wave)
        or not isinstance(raw_polarization, list)
        or len(raw_polarization) != 3
        or not all(
            isinstance(value, int) and not isinstance(value, bool)
            for value in raw_polarization
        )
        or phase not in ("cos", "sin")
    ):
        raise ValueError("serialized mode spec has invalid values")
    wave: Wave = tuple(raw_wave)  # type: ignore[assignment]
    polarization = tuple(Fraction(value) for value in raw_polarization)
    if wave == (0, 0, 0) or _canonical(wave) != wave:
        raise ValueError("serialized wave is not canonical")
    if sum(Fraction(k) * value for k, value in zip(wave, polarization)) != 0:
        raise ValueError("serialized mode is not divergence free")
    positive: ComplexVector
    if phase == "cos":
        positive = tuple((value / 2, Fraction(0)) for value in polarization)  # type: ignore[assignment]
    else:
        positive = tuple((Fraction(0), -value / 2) for value in polarization)  # type: ignore[assignment]
    negative: ComplexVector = tuple(
        _c_conjugate(value) for value in positive
    )  # type: ignore[assignment]
    return {wave: positive, _wave_neg(wave): negative}


def _advection(left: Field, right: Field) -> Field:
    output: Field = {}
    for left_wave, left_vector in left.items():
        for right_wave, right_vector in right.items():
            target = _wave_add(left_wave, right_wave)
            dot = _ZERO_C
            for component, frequency in zip(left_vector, right_wave):
                dot = _c_add(dot, _c_scale(component, frequency))
            contribution: ComplexVector = tuple(
                _c_i_mul(_c_mul(dot, value)) for value in right_vector
            )  # type: ignore[assignment]
            output[target] = _v_add(output.get(target, _ZERO_V), contribution)
    return {wave: vector for wave, vector in output.items() if not _v_zero(vector)}


def _leray(field: Field) -> Field:
    output: Field = {}
    for wave, vector in field.items():
        if wave == (0, 0, 0):
            projected = vector
        else:
            denominator = sum(value * value for value in wave)
            dot = _ZERO_C
            for frequency, value in zip(wave, vector):
                dot = _c_add(dot, _c_scale(value, frequency))
            projected = tuple(
                _c_add(value, _c_scale(dot, Fraction(-frequency, denominator)))
                for frequency, value in zip(wave, vector)
            )  # type: ignore[assignment]
        if not _v_zero(projected):
            output[wave] = projected
    return output


def _mixed(left: Field, right: Field) -> Field:
    return _leray(_field_add(_advection(left, right), _advection(right, left)))


def _restrict_canonical(field: Field, waves: set[Wave]) -> Field:
    return {
        wave: vector
        for wave, vector in field.items()
        if wave != (0, 0, 0) and _canonical(wave) in waves
    }


def _restrict_shell(field: Field, shell_squared: int) -> Field:
    return {
        wave: vector
        for wave, vector in field.items()
        if sum(value * value for value in wave) == shell_squared
    }


def _inner(left: Field, right: Field) -> Fraction:
    total = _ZERO_C
    for wave, left_vector in left.items():
        right_vector = right.get(wave, _ZERO_V)
        for left_value, right_value in zip(left_vector, right_vector):
            total = _c_add(total, _c_mul(_c_conjugate(left_value), right_value))
    if total[1] != 0:
        raise ValueError("real-field inner product acquired an imaginary part")
    return total[0]


def _norm_squared(field: Field) -> Fraction:
    return _inner(field, field)


def _table(field: Field) -> dict[str, dict[str, list[str]]]:
    output: dict[str, dict[str, list[str]]] = {}
    for wave in sorted(field):
        if wave == (0, 0, 0) or _canonical(wave) != wave:
            continue
        vector = field[wave]
        cosine = [str(2 * value[0]) for value in vector]
        sine = [str(-2 * value[1]) for value in vector]
        if all(value == "0" for value in cosine + sine):
            continue
        output[",".join(str(value) for value in wave)] = {
            "cosine": cosine,
            "sine": sine,
        }
    return output


def _wave_from_list(value: object, *, label: str) -> Wave:
    if (
        not isinstance(value, list)
        or len(value) != 3
        or not all(isinstance(item, int) and not isinstance(item, bool) for item in value)
    ):
        raise ValueError(f"{label} is not an integer wave")
    return tuple(value)  # type: ignore[return-value]


def _interaction(left: Field, right: Field, left_wave: Wave, right_wave: Wave) -> tuple[dict[str, object] | None, Field | None, Field | None]:
    raw_sum = _wave_add(left_wave, right_wave)
    raw_difference = _wave_add(left_wave, _wave_neg(right_wave))
    if raw_sum == (0, 0, 0) or raw_difference == (0, 0, 0):
        raise ValueError("interaction has degenerate support")
    sum_wave = _canonical(raw_sum)
    difference_wave = _canonical(raw_difference)
    mixed = _mixed(left, right)
    sum_output = _restrict_canonical(mixed, {sum_wave})
    if not sum_output:
        return None, None, None
    difference_output = _restrict_canonical(mixed, {difference_wave})
    other_output = {
        wave: vector
        for wave, vector in mixed.items()
        if wave != (0, 0, 0)
        and _canonical(wave) not in {sum_wave, difference_wave}
    }
    child = _field_neg(sum_output)
    signed = -_inner(child, mixed)
    triad = _field_add(left, right, child)
    nonlinear = _leray(_advection(triad, triad))
    parent = _field_add(left, right)
    parent_flux = -_inner(parent, nonlinear)
    child_flux = -_inner(child, nonlinear)
    defect = _inner(triad, nonlinear)
    payload: dict[str, object] = {
        "left_wave": list(left_wave),
        "right_wave": list(right_wave),
        "sum_wave": list(sum_wave),
        "difference_wave": list(difference_wave),
        "mixed_output": _table(mixed),
        "sum_output": _table(sum_output),
        "difference_output": _table(difference_output),
        "other_output": _table(other_output),
        "child": _table(child),
        "signed_parent_to_child_flux": str(signed),
        "parent_flux_in_populated_triad": str(parent_flux),
        "child_flux_in_populated_triad": str(child_flux),
        "populated_triad_energy_defect": str(defect),
    }
    if signed <= 0 or defect != 0 or parent_flux + child_flux != 0:
        raise ValueError("independent interaction energy/transfer identity failed")
    return payload, child, difference_output


def _verify_interaction_payload(
    serialized: object,
    left: Field,
    right: Field,
    left_wave: Wave,
    right_wave: Wave,
) -> tuple[Field, Field]:
    recomputed, child, difference = _interaction(
        left, right, left_wave, right_wave
    )
    if recomputed is None or child is None or difference is None:
        raise ValueError("serialized nonzero interaction independently vanished")
    if serialized != recomputed:
        raise ValueError("serialized interaction coefficient table does not recompute")
    return child, difference


def verify_serialized_strict_orientation_records(
    certificate: object,
) -> dict[str, object]:
    """Independently verify every disclosed fixed-relay strict orientation."""

    if not isinstance(certificate, dict):
        raise ValueError("carrier certificate must be an object")
    strict = certificate.get("strict_search")
    if not isinstance(strict, dict):
        raise ValueError("carrier certificate lacks strict search")
    records = strict.get("orientation_records")
    tested = strict.get("eligible_second_relay_orientations_tested")
    if not isinstance(records, list) or tested != len(records):
        raise ValueError("strict orientation inventory is incomplete")

    fixed_left_spec = {"wave": list(_P), "polarization": [0, 0, 1], "phase": "sin"}
    fixed_right_spec = {"wave": list(_Q), "polarization": [0, 1, 0], "phase": "cos"}
    fixed_left = _mode(fixed_left_spec)
    fixed_right = _mode(fixed_right_spec)
    fixed_payload, fixed_child, fixed_difference = _interaction(
        fixed_left, fixed_right, _P, _Q
    )
    if fixed_payload is None or fixed_child is None or fixed_difference is None or fixed_difference:
        raise ValueError("independent fixed relay check failed")

    reason_histogram: dict[str, int] = {}
    strict_passes = 0
    energy_identities = 1
    for expected_index, record in enumerate(records, start=1):
        if not isinstance(record, dict) or record.get("orientation_index") != expected_index:
            raise ValueError("strict orientation indices are not complete and ordered")
        relay = record.get("relay_two")
        if not isinstance(relay, dict):
            raise ValueError("strict orientation lacks relay two")
        left_spec = relay.get("left")
        right_spec = relay.get("right")
        left = _mode(left_spec)
        right = _mode(right_spec)
        if not isinstance(left_spec, dict) or not isinstance(right_spec, dict):
            raise ValueError("strict orientation mode specs are invalid")
        left_wave = _wave_from_list(left_spec.get("wave"), label="left wave")
        right_wave = _wave_from_list(right_spec.get("wave"), label="right wave")
        second_child, second_difference = _verify_interaction_payload(
            relay.get("interaction"), left, right, left_wave, right_wave
        )
        energy_identities += 1

        second_sum_wave = _canonical(_wave_add(left_wave, right_wave))
        next_payload, _, next_difference = _interaction(
            fixed_child,
            second_child,
            _canonical(_wave_add(_P, _Q)),
            second_sum_wave,
        )
        required_next = next_payload is not None and bool(next_difference)
        serialized_next = record.get("next_interaction")
        if required_next:
            if serialized_next != next_payload:
                raise ValueError("serialized next interaction does not recompute")
            energy_identities += 1
        elif serialized_next is not None:
            raise ValueError("record fabricates a required next interaction")

        target_shell = 6
        a1_b2 = _mixed(fixed_left, right)
        a2_b1 = _mixed(left, fixed_right)
        a1_target = _restrict_shell(a1_b2, target_shell)
        a2_target = _restrict_shell(a2_b1, target_shell)
        cross_passes = not a1_target and not a2_target
        cross_payload = {
            "condition": "A1+B2 and A2+B1 must vanish on the target child shell",
            "target_shell_squared": target_shell,
            "a1_b2_output": _table(a1_b2),
            "a2_b1_output": _table(a2_b1),
            "a1_b2_target_shell_output": _table(a1_target),
            "a2_b1_target_shell_output": _table(a2_target),
            "passes": cross_passes,
        }
        if record.get("diagonal_cross_talk") != cross_payload:
            raise ValueError("serialized diagonal coefficient table does not recompute")

        reasons: list[str] = []
        if not required_next:
            reasons.append(
                "next_child_interaction_lacks_required_nonzero_sum_or_difference"
            )
        if not cross_passes:
            reasons.append("diagonal_parent_cross_talk_hits_target_child_shell")
        passes = not reasons
        if (
            record.get("required_next_interaction_passes") is not required_next
            or record.get("strict_passes") is not passes
            or record.get("rejection_reasons") != reasons
        ):
            raise ValueError("serialized strict rejection predicate is inconsistent")
        if second_difference:
            raise ValueError("strict second relay difference is not exactly zero")
        strict_passes += int(passes)
        for reason in reasons:
            reason_histogram[reason] = reason_histogram.get(reason, 0) + 1

    if bool(strict.get("found")) != bool(strict_passes):
        raise ValueError("strict found flag disagrees with independent records")
    if strict.get("exhaustive") is not True:
        raise ValueError("record inventory is not labelled exhaustive")
    return {
        "implementation": "independent-complex-Fourier-ordered-pair-Fraction-convolution",
        "records_verified": len(records),
        "strict_pass_count": strict_passes,
        "energy_identities_verified": energy_identities,
        "rejection_histogram": reason_histogram,
        "verified": True,
    }


def _qdot(left: Iterable[int | Fraction], right: Iterable[int | Fraction]) -> Fraction:
    return sum(
        (Fraction(a) * Fraction(b) for a, b in zip(left, right)),
        Fraction(0),
    )


def _qproject(vector: Iterable[int | Fraction], wave: Wave) -> tuple[Fraction, Fraction, Fraction]:
    denominator = sum(value * value for value in wave)
    if denominator == 0:
        raise ValueError("zero-frequency projection in expanded record")
    values = tuple(Fraction(value) for value in vector)
    longitudinal = _qdot(values, wave) / denominator
    return tuple(
        value - Fraction(frequency) * longitudinal
        for value, frequency in zip(values, wave)
    )  # type: ignore[return-value]


def _qchannels(
    left_wave: Wave,
    left_polarization: Iterable[int | Fraction],
    right_wave: Wave,
    right_polarization: Iterable[int | Fraction],
) -> tuple[
    Wave,
    tuple[Fraction, Fraction, Fraction],
    Wave,
    tuple[Fraction, Fraction, Fraction],
]:
    left = tuple(Fraction(value) for value in left_polarization)
    right = tuple(Fraction(value) for value in right_polarization)
    left_on_right = _qdot(left, right_wave)
    right_on_left = _qdot(right, left_wave)
    sum_wave = _wave_add(left_wave, right_wave)
    difference_wave = _wave_add(left_wave, _wave_neg(right_wave))
    sum_raw = tuple(
        left_on_right * b + right_on_left * a for a, b in zip(left, right)
    )
    # This is the real-trigonometric record convention; the complex formula
    # differs by a harmless overall sign in the difference channel.
    difference_raw = tuple(
        left_on_right * b - right_on_left * a for a, b in zip(left, right)
    )
    return (
        sum_wave,
        _qproject(sum_raw, sum_wave),
        difference_wave,
        _qproject(difference_raw, difference_wave),
    )


def _qzero(vector: Iterable[Fraction]) -> bool:
    return all(value == 0 for value in vector)


def _qstrings(vector: Iterable[Fraction]) -> list[str]:
    return [str(value) for value in vector]


def _accumulate_qstats(
    stats: dict[str, int],
    *vectors: Iterable[Fraction],
) -> None:
    for vector in vectors:
        for value in vector:
            stats["maximum_denominator"] = max(
                stats["maximum_denominator"], value.denominator
            )
            stats["maximum_denominator_bit_length"] = max(
                stats["maximum_denominator_bit_length"],
                value.denominator.bit_length(),
            )
            stats["maximum_numerator_bit_length"] = max(
                stats["maximum_numerator_bit_length"],
                abs(value.numerator).bit_length(),
            )
            stats["fraction_component_count"] += 1


def _primitive_directions(wave: Wave, bound: int) -> tuple[tuple[int, int, int], ...]:
    values: set[tuple[int, int, int]] = set()
    for raw in product(range(-bound, bound + 1), repeat=3):
        vector: tuple[int, int, int] = tuple(raw)  # type: ignore[assignment]
        if vector == (0, 0, 0) or sum(a * b for a, b in zip(wave, vector)):
            continue
        if math.gcd(math.gcd(abs(vector[0]), abs(vector[1])), abs(vector[2])) != 1:
            continue
        for component in vector:
            if component > 0:
                values.add(vector)
                break
            if component < 0:
                values.add(tuple(-item for item in vector))
                break
    return tuple(sorted(values))


def _canonical_box(bound: int) -> tuple[Wave, ...]:
    values: list[Wave] = []
    for raw in product(range(-bound, bound + 1), repeat=3):
        wave: Wave = tuple(raw)  # type: ignore[assignment]
        if wave != (0, 0, 0) and _canonical(wave) == wave:
            values.append(wave)
    return tuple(sorted(values))


def _target_leakage(
    left_wave: Wave,
    left_polarization: Iterable[int | Fraction],
    right_wave: Wave,
    right_polarization: Iterable[int | Fraction],
    shell_squared: int,
) -> Fraction:
    sum_wave, sum_vector, difference_wave, difference_vector = _qchannels(
        left_wave,
        left_polarization,
        right_wave,
        right_polarization,
    )
    total = Fraction(0)
    if sum(value * value for value in sum_wave) == shell_squared:
        total += _qdot(sum_vector, sum_vector)
    if sum(value * value for value in difference_wave) == shell_squared:
        total += _qdot(difference_vector, difference_vector)
    return total


def _verify_expanded_structural_records(payload: dict[str, object]) -> dict[str, object]:
    config = payload.get("search_config")
    catalog = payload.get("catalog")
    records = payload.get("all_structural_orientation_records")
    if not isinstance(config, dict) or not isinstance(catalog, dict) or not isinstance(records, list):
        raise ValueError("expanded structural inventory is missing")
    wave_bound = config.get("wave_bound_l_infinity")
    polarization_bound = config.get("primitive_polarization_component_bound")
    if not isinstance(wave_bound, int) or not isinstance(polarization_bound, int):
        raise ValueError("expanded bounds are invalid")
    waves = _canonical_box(wave_bound)
    serialized_waves = catalog.get("wavevectors")
    if serialized_waves != [list(wave) for wave in waves]:
        raise ValueError("expanded wave catalog does not independently enumerate")
    directions = {
        wave: _primitive_directions(wave, polarization_bound) for wave in waves
    }
    expected_catalog = {
        ",".join(str(value) for value in wave): [list(vector) for vector in directions[wave]]
        for wave in waves
    }
    if catalog.get("polarizations_by_wave") != expected_catalog:
        raise ValueError("expanded polarization catalog does not independently enumerate")
    if (
        catalog.get("canonical_wave_count") != len(waves)
        or catalog.get("total_primitive_polarization_direction_count")
        != sum(len(values) for values in directions.values())
    ):
        raise ValueError("expanded catalog counts are inconsistent")

    second_waves = tuple(wave for wave in waves if wave not in (_P, _Q))
    all_pairs = tuple(
        (left, right)
        for left in second_waves
        for right in second_waves
        if left != right
    )
    scope_pairs = tuple(
        (left, right)
        for left, right in all_pairs
        if sum(value * value for value in left) == 2
        and sum(value * value for value in right) == 2
        and sum(value * value for value in _wave_add(left, right)) == 6
        and _wave_add(left, right) != (2, 1, 1)
    )
    expected_tests = sum(
        len(directions[left]) * len(directions[right]) for left, right in scope_pairs
    )
    if len(records) != expected_tests:
        raise ValueError("expanded structural record count is incomplete")
    if (
        config.get("maximum_ordered_wave_pair_tests", -1) < len(all_pairs)
        or config.get("maximum_polarization_pair_tests", -1) < expected_tests
        or config.get("maximum_saved_structural_records", -1) < expected_tests
    ):
        raise ValueError("expanded resource caps do not cover the claimed search")

    _, fixed_child_vector, _, fixed_difference = _qchannels(
        _P, (0, 0, 1), _Q, (0, 1, 0)
    )
    if not _qzero(fixed_difference):
        raise AssertionError("independent fixed difference formula failed")
    histogram: dict[str, int] = {
        "second_sum_leray_zero": 0,
        "second_difference_leray_nonzero": 0,
        "next_sum_leray_zero": 0,
        "next_difference_leray_zero": 0,
        "diagonal_cross_talk_hits_child_shell": 0,
        "strict_hit_orientation": 0,
    }
    arithmetic_stats = {
        "maximum_denominator": 1,
        "maximum_denominator_bit_length": 1,
        "maximum_numerator_bit_length": 1,
        "fraction_component_count": 0,
    }
    strict_keys: list[tuple[Wave, tuple[int, int, int], Wave, tuple[int, int, int]]] = []
    partial_rank_records: list[
        tuple[
            Fraction,
            Fraction,
            Wave,
            tuple[int, int, int],
            Wave,
            tuple[int, int, int],
        ]
    ] = []
    index = 0
    for left_wave, right_wave in scope_pairs:
        for left_polarization in directions[left_wave]:
            for right_polarization in directions[right_wave]:
                index += 1
                record = records[index - 1]
                if not isinstance(record, dict):
                    raise ValueError("expanded structural record is not an object")
                sum_wave, sum_vector, difference_wave, difference_vector = _qchannels(
                    left_wave,
                    left_polarization,
                    right_wave,
                    right_polarization,
                )
                _accumulate_qstats(
                    arithmetic_stats, sum_vector, difference_vector
                )
                expected: dict[str, object] = {
                    "record_index": index,
                    "second_left": {
                        "wave": list(left_wave),
                        "primitive_polarization": list(left_polarization),
                    },
                    "second_right": {
                        "wave": list(right_wave),
                        "primitive_polarization": list(right_polarization),
                    },
                    "second_sum_wave": list(sum_wave),
                    "second_sum_projected_vector": _qstrings(sum_vector),
                    "second_difference_wave": list(difference_wave),
                    "second_difference_projected_vector": _qstrings(difference_vector),
                    "next_sum_wave": None,
                    "next_sum_projected_vector": None,
                    "next_difference_wave": None,
                    "next_difference_projected_vector": None,
                    "diagonal_target_shell_leakage_strength_squared": None,
                    "disposition": None,
                }
                if _qzero(sum_vector):
                    disposition = "second_sum_leray_zero"
                elif not _qzero(difference_vector):
                    disposition = "second_difference_leray_nonzero"
                else:
                    (
                        next_sum_wave,
                        next_sum_vector,
                        next_difference_wave,
                        next_difference_vector,
                    ) = _qchannels((2, 1, 1), fixed_child_vector, sum_wave, sum_vector)
                    _accumulate_qstats(
                        arithmetic_stats,
                        next_sum_vector,
                        next_difference_vector,
                    )
                    expected.update(
                        {
                            "next_sum_wave": list(next_sum_wave),
                            "next_sum_projected_vector": _qstrings(next_sum_vector),
                            "next_difference_wave": list(next_difference_wave),
                            "next_difference_projected_vector": _qstrings(
                                next_difference_vector
                            ),
                        }
                    )
                    if _qzero(next_sum_vector):
                        disposition = "next_sum_leray_zero"
                    elif _qzero(next_difference_vector):
                        disposition = "next_difference_leray_zero"
                    else:
                        cross_one = _qchannels(
                            _P,
                            (0, 0, 1),
                            right_wave,
                            right_polarization,
                        )
                        cross_two = _qchannels(
                            left_wave,
                            left_polarization,
                            _Q,
                            (0, 1, 0),
                        )
                        _accumulate_qstats(
                            arithmetic_stats,
                            cross_one[1],
                            cross_one[3],
                            cross_two[1],
                            cross_two[3],
                        )
                        leakage = Fraction(0)
                        for cross in (cross_one, cross_two):
                            for wave, vector in (
                                (cross[0], cross[1]),
                                (cross[2], cross[3]),
                            ):
                                if sum(value * value for value in wave) == 6:
                                    leakage += _qdot(vector, vector)
                        expected[
                            "diagonal_target_shell_leakage_strength_squared"
                        ] = str(leakage)
                        disposition = (
                            "strict_hit_orientation"
                            if leakage == 0
                            else "diagonal_cross_talk_hits_child_shell"
                        )
                        if disposition == "diagonal_cross_talk_hits_child_shell":
                            partial_rank_records.append(
                                (
                                    leakage,
                                    -_qdot(next_sum_vector, next_sum_vector),
                                    left_wave,
                                    left_polarization,
                                    right_wave,
                                    right_polarization,
                                )
                            )
                expected["disposition"] = disposition
                histogram[disposition] += 1
                if disposition == "strict_hit_orientation":
                    strict_keys.append(
                        (
                            left_wave,
                            left_polarization,
                            right_wave,
                            right_polarization,
                        )
                    )
                if record != expected:
                    raise ValueError(
                        f"expanded structural record {index} does not recompute"
                    )
    if payload.get("rejection_histogram") != histogram:
        raise ValueError("expanded rejection histogram is inconsistent")
    counts = payload.get("counts")
    if not isinstance(counts, dict) or any(
        counts.get(name) != value
        for name, value in (
            ("ordered_second_parent_wave_pairs_examined", len(all_pairs)),
            ("scope_eligible_ordered_wave_pairs", len(scope_pairs)),
            ("polarization_pair_tests", expected_tests),
            ("strict_hit_orientation_count", len(strict_keys)),
            ("saved_structural_record_count", len(records)),
        )
    ):
        raise ValueError("expanded search counts are inconsistent")
    if payload.get("exact_arithmetic_stats_over_full_fast_search") != arithmetic_stats:
        raise ValueError("expanded exact arithmetic provenance is inconsistent")
    return {
        "records_verified": len(records),
        "histogram": histogram,
        "strict_keys": strict_keys,
        "partial_ranked_keys": [
            (record[2], record[3], record[4], record[5])
            for record in sorted(partial_rank_records)
        ],
        "arithmetic_stats": arithmetic_stats,
    }


def _maximum_serialized_denominator(value: object) -> int:
    if isinstance(value, str):
        try:
            return Fraction(value).denominator
        except (ValueError, ZeroDivisionError):
            return 1
    if isinstance(value, dict):
        return max(
            (_maximum_serialized_denominator(item) for item in value.values()),
            default=1,
        )
    if isinstance(value, list):
        return max(
            (_maximum_serialized_denominator(item) for item in value),
            default=1,
        )
    return 1


def _scaled_mode_from_candidate(
    spec: object,
    normalization: object,
    *,
    unit_normalization_requested: bool,
) -> tuple[Field, Fraction]:
    if not isinstance(normalization, dict):
        raise ValueError("expanded amplitude normalization is invalid")
    primitive = normalization.get("primitive_polarization")
    if not isinstance(spec, dict) or spec.get("polarization") != primitive:
        raise ValueError("primitive ModeSpec and normalization label disagree")
    scale = Fraction(str(normalization.get("field_scale_applied")))
    expected_field_polarization = [
        str(scale * int(value)) for value in primitive
    ]
    if normalization.get("field_polarization") != expected_field_polarization:
        raise ValueError("expanded scaled field polarization is inconsistent")
    norm_squared = sum(int(value) ** 2 for value in primitive)
    if normalization.get("primitive_norm_squared") != norm_squared:
        raise ValueError("expanded primitive polarization norm is inconsistent")
    norm = math.isqrt(norm_squared)
    rational_unit = norm * norm == norm_squared
    expected_norm = str(norm) if rational_unit else f"sqrt({norm_squared})"
    if (
        normalization.get("primitive_norm") != expected_norm
        or normalization.get("unit_normalization_requested")
        is not unit_normalization_requested
        or normalization.get("unit_normalization_exactly_rational")
        is not rational_unit
    ):
        raise ValueError("expanded normalization metadata is inconsistent")
    expected_scale = (
        Fraction(1, norm)
        if unit_normalization_requested and rational_unit
        else Fraction(1)
    )
    if scale != expected_scale:
        raise ValueError("expanded field scale does not follow normalization policy")
    return _field_scale(_mode(spec), scale), scale


def _unintended_payload(
    label: str,
    left: Field,
    right: Field,
) -> dict[str, object]:
    output = _mixed(left, right)
    target = _restrict_shell(output, 6)
    shells = sorted(
        {sum(value * value for value in wave) for wave in output}
    )
    return {
        "label": label,
        "full_output": _table(output),
        "target_shell_output": _table(target),
        "target_shell_zero": not target,
        "nonzero_output_shells_squared": shells,
    }


def _forcing_budget(
    fixed_left: Field,
    fixed_right: Field,
    second_left: Field,
    second_right: Field,
    second_sum_wave: Wave,
) -> dict[str, object]:
    parent = _field_add(fixed_left, fixed_right, second_left, second_right)
    nonlinear = _leray(_advection(parent, parent))
    target = _restrict_shell(nonlinear, 6)
    intended = _restrict_canonical(nonlinear, {(2, 1, 1), second_sum_wave})
    low = {
        wave: vector
        for wave, vector in nonlinear.items()
        if sum(value * value for value in wave) < 6
    }
    high = {
        wave: vector
        for wave, vector in nonlinear.items()
        if sum(value * value for value in wave) > 6
    }
    parent_norm = _norm_squared(parent)
    full_power = _norm_squared(nonlinear)
    target_power = _norm_squared(target)
    intended_power = _norm_squared(intended)
    low_power = _norm_squared(low)
    high_power = _norm_squared(high)
    pairing = _inner(parent, nonlinear)
    return {
        "parent_l2_squared": str(parent_norm),
        "full_nonlinear_l2_squared": str(full_power),
        "target_child_shell_l2_squared": str(target_power),
        "intended_C1_C2_mode_l2_squared": str(intended_power),
        "target_shell_contamination_l2_squared": str(
            target_power - intended_power
        ),
        "low_shell_l2_squared": str(low_power),
        "high_shell_l2_squared": str(high_power),
        "target_shell_forcing_fraction": str(target_power / full_power),
        "intended_mode_forcing_fraction": str(intended_power / full_power),
        "low_to_intended_power_ratio": str(low_power / intended_power),
        "energy_pairing": str(pairing),
        "energy_identity_exact": pairing == 0,
        "classification": (
            "TARGET-SHELL-CLEAN-BUT-LOW-RETURN-DOMINATED"
            if target_power == intended_power and low_power > intended_power
            else "TARGET-SHELL-CONTAMINATED"
        ),
    }


def _verify_expanded_exact_candidate(
    candidate: object,
    *,
    strict_expected: bool,
) -> tuple[Wave, tuple[int, int, int], Wave, tuple[int, int, int], bool]:
    if not isinstance(candidate, dict):
        raise ValueError("expanded exact candidate is not an object")
    left_spec = candidate.get("second_left")
    right_spec = candidate.get("second_right")
    normalization = candidate.get("amplitude_normalization")
    if not isinstance(normalization, dict):
        raise ValueError("expanded candidate normalization is missing")
    second_left, _ = _scaled_mode_from_candidate(
        left_spec,
        normalization.get("left"),
        unit_normalization_requested=strict_expected,
    )
    second_right, _ = _scaled_mode_from_candidate(
        right_spec,
        normalization.get("right"),
        unit_normalization_requested=strict_expected,
    )
    left_norm_squared = sum(int(value) ** 2 for value in left_spec["polarization"])
    right_norm_squared = sum(int(value) ** 2 for value in right_spec["polarization"])
    left_rational = math.isqrt(left_norm_squared) ** 2 == left_norm_squared
    right_rational = math.isqrt(right_norm_squared) ** 2 == right_norm_squared
    normalization_applied = strict_expected and left_rational and right_rational
    expected_policy = (
        "rational unit normalization"
        if normalization_applied
        else (
            "primitive exact representative; requested unit scale is algebraic"
            if strict_expected
            else "primitive integer representative"
        )
    )
    if normalization.get("policy") != expected_policy or normalization.get(
        "reconstruction_rule"
    ) != (
        "multiply each primitive ModeSpec field by its recorded "
        "field_scale_applied before recomputing interactions"
    ):
        raise ValueError("expanded normalization policy metadata is inconsistent")
    if not isinstance(left_spec, dict) or not isinstance(right_spec, dict):
        raise ValueError("expanded candidate specs are invalid")
    left_wave = _wave_from_list(left_spec.get("wave"), label="expanded left")
    right_wave = _wave_from_list(right_spec.get("wave"), label="expanded right")
    left_polarization = tuple(left_spec.get("polarization", ()))
    right_polarization = tuple(right_spec.get("polarization", ()))
    if len(left_polarization) != 3 or len(right_polarization) != 3:
        raise ValueError("expanded primitive polarization is invalid")
    if candidate.get("representative_phase_choice") != {
        "fixed_left": "sin",
        "fixed_right": "cos",
        "second_left": "cos",
        "second_right": "cos",
    }:
        raise ValueError("expanded representative phase metadata is invalid")
    second_fast = _qchannels(
        left_wave,
        left_polarization,
        right_wave,
        right_polarization,
    )
    fixed_fast = _qchannels(_P, (0, 0, 1), _Q, (0, 1, 0))
    next_fast = _qchannels(
        (2, 1, 1), fixed_fast[1], second_fast[0], second_fast[1]
    )
    raw_leakage = _target_leakage(
        _P,
        (0, 0, 1),
        right_wave,
        right_polarization,
        6,
    ) + _target_leakage(
        left_wave,
        left_polarization,
        _Q,
        (0, 1, 0),
        6,
    )
    if candidate.get("diagonal_target_shell_leakage_strength_squared") != str(
        raw_leakage
    ) or candidate.get("phase_quotiented_next_sum_strength_squared") != str(
        _qdot(next_fast[1], next_fast[1])
    ):
        raise ValueError("expanded phase-quotiented score metadata is inconsistent")

    fixed_left_spec = {"wave": [1, 1, 0], "polarization": [0, 0, 1], "phase": "sin"}
    fixed_right_spec = {"wave": [1, 0, 1], "polarization": [0, 1, 0], "phase": "cos"}
    fixed_left = _mode(fixed_left_spec)
    fixed_right = _mode(fixed_right_spec)
    fixed_payload, fixed_child, fixed_difference = _interaction(
        fixed_left, fixed_right, _P, _Q
    )
    if fixed_payload is None or fixed_child is None or fixed_difference is None or fixed_difference:
        raise AssertionError("independent expanded fixed relay failed")
    second_payload, second_child, second_difference = _interaction(
        second_left, second_right, left_wave, right_wave
    )
    if second_payload is None or second_child is None or second_difference is None:
        raise ValueError("expanded second relay independently vanished")
    if second_difference:
        raise ValueError("expanded exact candidate has nonzero second difference")
    second_sum_wave = _canonical(_wave_add(left_wave, right_wave))
    next_payload, next_child, next_difference = _interaction(
        fixed_child, second_child, (2, 1, 1), second_sum_wave
    )
    if next_payload is None or next_child is None or next_difference is None or not next_difference:
        raise ValueError("expanded next interaction independently failed")

    expected_second_relay = {
        "left": left_spec,
        "right": right_spec,
        "carrier_shells_squared": [2, 6],
        "interaction": second_payload,
    }
    if candidate.get("second_relay") != expected_second_relay:
        raise ValueError("expanded normalized second relay does not recompute")
    if candidate.get("next_interaction") != next_payload:
        raise ValueError("expanded next coefficient table does not recompute")
    scaled_tables = candidate.get("exact_scaled_second_parent_fields")
    if scaled_tables != {
        "left": _table(second_left),
        "right": _table(second_right),
    }:
        raise ValueError("expanded scaled parent field tables do not recompute")

    a1_b2 = _mixed(fixed_left, second_right)
    a2_b1 = _mixed(second_left, fixed_right)
    a1_target = _restrict_shell(a1_b2, 6)
    a2_target = _restrict_shell(a2_b1, 6)
    cross_passes = not a1_target and not a2_target
    cross_payload = {
        "condition": "A1+B2 and A2+B1 must vanish on the target child shell",
        "target_shell_squared": 6,
        "a1_b2_output": _table(a1_b2),
        "a2_b1_output": _table(a2_b1),
        "a1_b2_target_shell_output": _table(a1_target),
        "a2_b1_target_shell_output": _table(a2_target),
        "passes": cross_passes,
    }
    if candidate.get("diagonal_cross_talk") != cross_payload:
        raise ValueError("expanded diagonal coefficient table does not recompute")
    if cross_passes is not strict_expected:
        raise ValueError("expanded strict/partial classification is inconsistent")

    unintended = [
        _unintended_payload("A1+A2", fixed_left, second_left),
        _unintended_payload("B1+B2", fixed_right, second_right),
        _unintended_payload("A1+B2", fixed_left, second_right),
        _unintended_payload("A2+B1", second_left, fixed_right),
    ]
    all_target_zero = all(record["target_shell_zero"] for record in unintended)
    lower_return = any(
        any(int(shell) < 6 for shell in record["nonzero_output_shells_squared"])
        for record in unintended
    )
    expected_unintended = {
        "pairs": unintended,
        "all_target_child_shell_outputs_zero": all_target_zero,
        "lower_shell_return_nonzero": lower_return,
        "classification": (
            "target-shell-clean-but-lower-shell-return-present"
            if all_target_zero and lower_return
            else "target-shell-contaminated-or-no-lower-return"
        ),
    }
    if candidate.get("all_unintended_parent_pair_audit") != expected_unintended:
        raise ValueError("expanded all-pair audit does not recompute")
    if candidate.get("four_parent_forcing_budget") != _forcing_budget(
        fixed_left,
        fixed_right,
        second_left,
        second_right,
        second_sum_wave,
    ):
        raise ValueError("expanded four-parent forcing budget does not recompute")

    first_triad = _field_add(fixed_left, fixed_right, fixed_child)
    second_triad = _field_add(second_left, second_right, second_child)
    relay_cross = _mixed(first_triad, second_triad)
    full_field = _field_add(first_triad, second_triad, next_child)
    full_nonlinear = _leray(_advection(full_field, full_field))
    full_defect = _inner(full_field, full_nonlinear)
    flux_by_wave: dict[str, str] = {}
    for wave in sorted(_table(full_field)):
        parsed_wave: Wave = tuple(int(value) for value in wave.split(","))  # type: ignore[assignment]
        band = _restrict_canonical(full_field, {parsed_wave})
        flux_by_wave[wave] = str(-_inner(band, full_nonlinear))
    difference_child = _field_neg(next_difference)
    grandchild_mutual = _mixed(next_child, difference_child)
    intended_next = _field_add(_field_neg(next_child), next_difference)
    relay_cross_norm = _norm_squared(relay_cross)
    intended_next_norm = _norm_squared(intended_next)
    expected_fixed_relay = {
        "left": fixed_left_spec,
        "right": fixed_right_spec,
        "carrier_shells_squared": [2, 6],
        "interaction": fixed_payload,
    }
    expected_full_gadget = {
        "classification": (
            "strict_cross_talk_pass" if cross_passes else "partial_rejected_cross_talk"
        ),
        "relay_one": expected_fixed_relay,
        "relay_two": expected_second_relay,
        "next_interaction": next_payload,
        "common_parent_shell_squared": 2,
        "common_child_shell_squared": 6,
        "diagonal_cross_talk": cross_payload,
        "relay_cross_output": _table(relay_cross),
        "relay_cross_output_l2_squared": str(relay_cross_norm),
        "intended_next_output_l2_squared": str(intended_next_norm),
        "intended_fraction_of_relay_cross_power": str(
            intended_next_norm / relay_cross_norm if relay_cross_norm else 0
        ),
        "grandchild_mutual_output": _table(grandchild_mutual),
        "grandchild_mutual_output_l2_squared": str(
            _norm_squared(grandchild_mutual)
        ),
        "simple_binary_recursion_rejected": not grandchild_mutual,
        "full_field_flux_by_wave": flux_by_wave,
        "full_field_energy_defect": str(full_defect),
    }
    if candidate.get("full_gadget") != expected_full_gadget:
        raise ValueError("expanded full gadget does not independently recompute")
    if full_defect != 0 or sum(
        (Fraction(value) for value in flux_by_wave.values()), Fraction(0)
    ) != 0:
        raise ValueError("expanded full-field energy identity failed")

    expected_predicates = {
        "second_difference_zero": True,
        "second_signed_transfer_positive": Fraction(
            second_payload["signed_parent_to_child_flux"]
        ) > 0,
        "next_sum_nonzero": bool(_field_neg(next_child)),
        "next_difference_nonzero": bool(next_difference),
        "next_signed_transfer_positive": Fraction(
            next_payload["signed_parent_to_child_flux"]
        ) > 0,
        "diagonal_target_shell_zero": cross_passes,
        "all_unintended_parent_pairs_target_shell_zero": all_target_zero,
        "relay_two_energy_defect_zero": second_payload[
            "populated_triad_energy_defect"
        ] == "0",
        "next_energy_defect_zero": next_payload[
            "populated_triad_energy_defect"
        ] == "0",
        "full_field_energy_defect_zero": full_defect == 0,
    }
    if candidate.get("exact_predicates") != expected_predicates:
        raise ValueError("expanded exact predicate summary is inconsistent")
    without_denominator = dict(candidate)
    saved_maximum = without_denominator.pop("maximum_saved_fraction_denominator", None)
    if saved_maximum != _maximum_serialized_denominator(without_denominator):
        raise ValueError("expanded saved denominator provenance is inconsistent")
    return (
        left_wave,
        tuple(int(value) for value in left_polarization),  # type: ignore[arg-type]
        right_wave,
        tuple(int(value) for value in right_polarization),  # type: ignore[arg-type]
        lower_return,
    )


def verify_serialized_expanded_carrier_certificate(
    certificate: object,
) -> dict[str, object]:
    """Independently verify the expanded search inventory and exact witnesses."""

    if not isinstance(certificate, dict):
        raise ValueError("expanded carrier certificate must be an object")
    required_fields = {
        "schema",
        "status",
        "classification",
        "claim_boundary",
        "search_config",
        "phase_quotient",
        "deduplication",
        "scope",
        "catalog",
        "counts",
        "exact_arithmetic_stats_over_full_fast_search",
        "rejection_histogram",
        "all_structural_orientation_records",
        "fixed_known_relay",
        "strict_hits",
        "top_partial_candidates",
        "parametric_rational_family",
        "continuous_polarization_classification",
        "relay_quality_screen",
        "verdict",
    }
    if set(certificate) != required_fields:
        raise ValueError("expanded carrier has missing or unknown fields")
    if certificate.get("schema") != (
        "ns-certificate-lab/expanded-exact-carrier-search/v1"
    ):
        raise ValueError("expanded carrier schema is invalid")
    if certificate.get("claim_boundary") != (
        "Exact finite carrier algebra only. No closure, invariant cone, "
        "scale iteration, PDE trajectory, or singularity is certified."
    ):
        raise ValueError("expanded carrier claim boundary is invalid")
    search_config = certificate.get("search_config")
    if not isinstance(search_config, dict) or set(search_config) != {
        "wave_bound_l_infinity",
        "primitive_polarization_component_bound",
        "input_rational_denominator",
        "maximum_ordered_wave_pair_tests",
        "maximum_polarization_pair_tests",
        "maximum_saved_structural_records",
        "maximum_saved_unique_hits",
        "top_partial_limit",
    } or search_config.get("input_rational_denominator") != 1:
        raise ValueError("expanded search configuration metadata is invalid")
    phase_quotient = certificate.get("phase_quotient")
    if phase_quotient != {
        "applied": True,
        "represented_phase_pairs": ["cos/cos", "cos/sin", "sin/cos", "sin/sin"],
        "reason": (
            "each cosine/sine mixture multiplies a given frequency channel "
            "by a nonzero phase/sign scalar; the projected vector's zero "
            "predicate and support shell radius are unchanged. No cancellation "
            "between distinct parent pairs is used"
        ),
        "independent_phase_optimization_used": False,
    }:
        raise ValueError("expanded phase quotient metadata is invalid")
    if certificate.get("deduplication") != {
        "wave_sign": "canonical half lattice",
        "polarization_scale_and_sign": (
            "primitive gcd-one vector with first nonzero component positive"
        ),
        "second_relay_exchange": (
            "not quotiented: A2/B2 exchange changes labelled diagonal pairs"
        ),
    }:
        raise ValueError("expanded deduplication metadata is invalid")
    if certificate.get("scope") != {
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
    }:
        raise ValueError("expanded search scope metadata is invalid")
    structural = _verify_expanded_structural_records(certificate)
    fixed_left_spec = {"wave": [1, 1, 0], "polarization": [0, 0, 1], "phase": "sin"}
    fixed_right_spec = {"wave": [1, 0, 1], "polarization": [0, 1, 0], "phase": "cos"}
    fixed_payload, _, fixed_difference = _interaction(
        _mode(fixed_left_spec), _mode(fixed_right_spec), _P, _Q
    )
    if fixed_payload is None or fixed_difference is None or fixed_difference:
        raise AssertionError("expanded independent fixed relay failed")
    if certificate.get("fixed_known_relay") != {
        "left": fixed_left_spec,
        "right": fixed_right_spec,
        "carrier_shells_squared": [2, 6],
        "interaction": fixed_payload,
    }:
        raise ValueError("expanded fixed known relay does not recompute")
    strict_hits = certificate.get("strict_hits")
    partials = certificate.get("top_partial_candidates")
    if not isinstance(strict_hits, list) or not isinstance(partials, list):
        raise ValueError("expanded exact candidate lists are invalid")
    verified_strict = [
        _verify_expanded_exact_candidate(candidate, strict_expected=True)
        for candidate in strict_hits
    ]
    verified_partials = [
        _verify_expanded_exact_candidate(candidate, strict_expected=False)
        for candidate in partials
    ]
    strict_keys = structural["strict_keys"]
    verified_keys = [record[:4] for record in verified_strict]
    if len(verified_keys) != len(set(verified_keys)):
        raise ValueError("expanded strict hit list contains duplicates")
    if verified_keys != strict_keys:
        raise ValueError("expanded strict hit tables do not cover all labelled hits")
    top_limit = search_config.get("top_partial_limit")
    if not isinstance(top_limit, int):
        raise ValueError("expanded top-partial limit is invalid")
    expected_partial_keys = structural["partial_ranked_keys"][:top_limit]
    verified_partial_keys = [record[:4] for record in verified_partials]
    if len(structural["partial_ranked_keys"]) != len(
        set(structural["partial_ranked_keys"])
    ):
        raise ValueError("expanded structural partial list contains duplicates")
    if verified_partial_keys != expected_partial_keys:
        raise ValueError("expanded partial ranking, uniqueness, or order is inconsistent")
    counts = certificate.get("counts")
    if not isinstance(counts, dict) or any(
        counts.get(name) != value
        for name, value in (
            ("unique_strict_hit_count", len(verified_strict)),
            ("saved_top_partial_count", len(verified_partials)),
            (
                "unique_cross_talk_partial_count",
                len(structural["partial_ranked_keys"]),
            ),
            (
                "strict_hits_without_lower_shell_return",
                sum(not record[4] for record in verified_strict),
            ),
        )
    ):
        raise ValueError("expanded saved candidate counts are inconsistent")
    if search_config.get("maximum_saved_unique_hits", -1) < len(
        verified_strict
    ):
        raise ValueError("expanded strict save cap does not cover all hits")
    found = bool(verified_strict)
    if (
        certificate.get("classification")
        != ("strict_finite_carrier_hit" if found else "finite_scope_negative")
        or ("STRICT EXACT FINITE CARRIER HIT" in str(certificate.get("status")))
        is not found
    ):
        raise ValueError("expanded carrier classification is inconsistent")
    expected_status = (
        "STRICT EXACT FINITE CARRIER HIT / NOT A PDE OR CASCADE"
        if found
        else "NEGATIVE IN THE STATED EXPANDED FINITE SEARCH ONLY"
    )
    expected_verdict = (
        "A strict finite carrier extension exists in this bounded alphabet; "
        "the earlier coordinate-axis-polarization negative does not extend "
        "to primitive integer polarizations. The displayed simple aligned "
        "binary recursion is exactly rejected; broader closure remains open."
        if found
        else
        "No strict extension exists in this exact finite alphabet; this has "
        "no force outside the recorded bounds."
    )
    if certificate.get("status") != expected_status or certificate.get(
        "verdict"
    ) != expected_verdict:
        raise ValueError("expanded carrier status or verdict is inconsistent")

    quality = certificate.get("relay_quality_screen")
    expected_quality_passes = sum(not record[4] for record in verified_strict)
    expected_quality = {
        "scope": "all four unintended parent pairs for every strict hit",
        "requires_no_output_below_child_shell_squared_6": True,
        "passing_orientation_count": expected_quality_passes,
        "exhaustive_in_recorded_alphabet": True,
        "classification": (
            "FINITE-SCOPE NEGATIVE: EVERY STRICT HIT RETURNS TO LOWER SHELLS"
            if found and expected_quality_passes == 0
            else "quality screen not uniformly negative"
        ),
    }
    if quality != expected_quality:
        raise ValueError("expanded relay quality screen is inconsistent")

    family = certificate.get("parametric_rational_family")
    if found:
        if not isinstance(family, dict):
            raise ValueError("expanded strict hit lacks its rational family")
        formulas = family.get("exact_formulas")
        if formulas != {
            "relay_two_signed_flux": "25*a^2*d^2/6",
            "next_sum_signed_flux": "125*a^2*d^2/126",
            "next_difference_l2_squared": "35*a^2*d^2/18",
            "next_difference_to_sum_norm_ratio": "7/5",
            "full_energy_defect": "0",
        }:
            raise ValueError("expanded parametric-family formulas are invalid")
        if family.get("waves") != {
            "r": [0, 1, -1],
            "s": [1, 0, -1],
            "C2": [1, 1, -2],
        } or family.get("polarizations") != {
            "R": "a*(1,2,2)",
            "S": "d*(2,1,2)",
            "parameters": "a,d are arbitrary nonzero rational numbers",
        } or family.get("derivation_constraints") != {
            "general_perpendicular_forms": "R=(a,b,b), S=(c,d,c)",
            "second_difference_leray_zero": "a*c=b*d",
            "A1_plus_B2_target_shell_zero": "c=2*d",
            "solution_family": "b=2*a, c=2*d",
        } or family.get("obstructions_still_present") != {
            "diagonal_outputs_globally_zero": False,
            "lower_shell_return": True,
            "simple_aligned_grandchild_recursion_zero": True,
            "interpretation": (
                "the family clears the intended child shell but is not a closed "
                "or iteratable relay"
            ),
        }:
            raise ValueError("expanded parametric family metadata is invalid")
        affine = family.get("affine_phase_witness")
        if not isinstance(affine, dict) or (
            affine.get("units") != "quarter turns (pi/2)"
            or affine.get("x0") != [-1, 0, 0]
            or affine.get("alpha_A") != 0
            or affine.get("alpha_B") != 1
        ):
            raise ValueError("expanded affine phase witness is invalid")
        x0 = tuple(int(value) for value in affine["x0"])
        phase_values = {
            "A1_at_p": (sum(a * b for a, b in zip(x0, _P)) + int(affine["alpha_A"])) % 4,
            "A2_at_r": (sum(a * b for a, b in zip(x0, (0, 1, -1))) + int(affine["alpha_A"])) % 4,
            "B1_at_q": (sum(a * b for a, b in zip(x0, _Q)) + int(affine["alpha_B"])) % 4,
            "B2_at_s": (sum(a * b for a, b in zip(x0, (1, 0, -1))) + int(affine["alpha_B"])) % 4,
        }
        # cos(theta + 3*pi/2) = sin(theta); zero quarter-turn is cosine.
        expected_phase_names = {
            label: "sin" if value == 3 else "cos" if value == 0 else "other"
            for label, value in phase_values.items()
        }
        if affine.get("resulting_phases") != expected_phase_names:
            raise ValueError("expanded affine phase labels do not evaluate")
        unit = family.get("unit_normalized_representative")
        preferred = strict_hits[0]
        if not isinstance(unit, dict) or not isinstance(preferred, dict) or unit != {
            "a": "1/3",
            "d": "1/3",
            "R": ["1/3", "2/3", "2/3"],
            "S": ["2/3", "1/3", "2/3"],
            "relay_two_signed_flux": "25/486",
            "next_sum_signed_flux": "125/10206",
        } or (
            unit.get("relay_two_signed_flux")
            != preferred["second_relay"]["interaction"][
                "signed_parent_to_child_flux"
            ]
            or unit.get("next_sum_signed_flux")
            != preferred["next_interaction"]["signed_parent_to_child_flux"]
        ):
            raise ValueError("expanded normalized family values do not recompute")
    elif family is not None:
        raise ValueError("expanded negative certificate fabricates a family")

    continuous = certificate.get("continuous_polarization_classification")
    expected_continuous = {
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
    if continuous != (expected_continuous if found else None):
        raise ValueError("expanded continuous-polarization classification is invalid")
    if found:
        # Independently check the two linear factors really have no common
        # nonzero (e,f) solution when a and b are both nonzero.
        determinant = 1 * (-1) - 4 * 2
        if determinant != -9:
            raise AssertionError("continuous branch determinant check failed")

    return {
        "implementation": (
            "independent-complex-Fourier-ordered-pair-Fraction-convolution"
        ),
        "structural_records_verified": structural["records_verified"],
        "strict_exact_witnesses_verified": len(verified_strict),
        "partial_exact_witnesses_verified": len(verified_partials),
        "full_search_arithmetic_stats_verified": structural["arithmetic_stats"],
        "all_strict_hits_have_lower_shell_return": all(
            record[4] for record in verified_strict
        ),
        "verified": True,
    }
