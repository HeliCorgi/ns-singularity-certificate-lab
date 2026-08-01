r"""Exact zero-padded local-FFT convolution for mesoscopic Leray clouds.

This is a finite binary64 convolution backend, not a singularity proof.  It
computes every coefficient of ``P((u . grad)u)`` for the two angled parent
boxes used by :mod:`ns_certificate_lab.mesoscopic_cloud_scaling`.  The FFT is
only an algorithm for a *linear* convolution of offset arrays: an input side
has length ``L = 2 W - 1`` and every transform is padded to
``K = 2 L - 1 = 4 W - 3``.  Hence no periodic wrap-around is present.

The Fourier-series convention is

``u(x) = sum_k u_hat[k] exp(i k.x)`` and
``||u||_2^2 = sum_k |u_hat[k]|^2``

where the torus measure is normalized to one.  NumPy's default FFT pair is
used only for coefficient convolution: ``ifft(fft(a) fft(b))`` is the linear
convolution after the stated zero padding.  There is consequently no extra
``K^3`` Parseval factor in any reported Fourier-series norm.

Different carrier sums can have overlapping output boxes.  They are retained
as separate in-memory blocks while the convolution is built, then summed on a
disjoint partition of global wave-number space before norms, tagged
projections, divergence, reality, and energy pairing are evaluated.  Thus an
overlap is neither discarded nor counted twice.
"""

from __future__ import annotations

import itertools
import math
import time
from dataclasses import asdict, dataclass
from typing import Iterator, Mapping, Sequence

import numpy as np
import numpy.typing as npt

from ns_certificate_lab.mesoscopic_cloud_scaling import (
    MesoscopicCloudConfig,
    MesoscopicCloudMetrics,
)
import ns_certificate_lab.mesoscopic_cloud_scaling as _scaling


Wave = tuple[int, int, int]
ComplexArray = npt.NDArray[np.complex128]

_P_POLARIZATION = np.asarray((0.0, 0.0, 1.0), dtype=np.complex128)
_Q_POLARIZATION = np.asarray((0.0, 1.0, 0.0), dtype=np.complex128)
_SPECIES = ("p+", "p-", "q+", "q-")

__all__ = [
    "LocalFFTCloudDiagnostics",
    "local_fft_leray_coefficients",
    "measure_local_fft_cloud",
    "measure_local_fft_mesoscopic_metrics",
]


def _add(left: Wave, right: Wave) -> Wave:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def _scale(multiplier: int, wave: Wave) -> Wave:
    return tuple(multiplier * value for value in wave)  # type: ignore[return-value]


def _species_direction(config: MesoscopicCloudConfig, species: str) -> Wave:
    p = tuple(config.p_direction)
    q = tuple(config.q_direction)
    if species == "p+":
        return p  # type: ignore[return-value]
    if species == "p-":
        return _scale(-1, p)  # type: ignore[arg-type]
    if species == "q+":
        return q  # type: ignore[return-value]
    if species == "q-":
        return _scale(-1, q)  # type: ignore[arg-type]
    raise ValueError(f"unknown species {species!r}")


def _species_family(species: str) -> str:
    if species.startswith("p"):
        return "p"
    if species.startswith("q"):
        return "q"
    raise ValueError(f"unknown species {species!r}")


def _normalization(config: MesoscopicCloudConfig) -> float:
    width = config.width
    one_dimensional_l2 = (2.0 * width * width + 1.0) / (3.0 * width)
    raw_norm_squared = 4.0 * one_dimensional_l2**3
    target_norm_squared = 2.0 * config.energy_constant / config.base_scale
    return math.sqrt(target_norm_squared / raw_norm_squared)


def _positive_species_array(
    config: MesoscopicCloudConfig,
    family: str,
    normalization: float,
) -> ComplexArray:
    """Return one positive box indexed by offsets ``[-d,d]^3``."""

    width = config.width
    d = width - 1
    axis = np.arange(-d, d + 1, dtype=np.int64)
    ox = axis[:, None, None]
    oy = axis[None, :, None]
    oz = axis[None, None, :]
    if family == "p":
        direction = np.asarray(config.p_direction, dtype=np.int64)
        polarization = _P_POLARIZATION
        phase_offset = config.phase_p
    elif family == "q":
        direction = np.asarray(config.q_direction, dtype=np.int64)
        polarization = _Q_POLARIZATION
        phase_offset = config.phase_q
    else:
        raise ValueError(f"unknown family {family!r}")

    kx = config.base_scale * direction[0] + ox
    ky = config.base_scale * direction[1] + oy
    kz = config.base_scale * direction[2] + oz
    wave_squared = (
        kx.astype(np.float64) ** 2
        + ky.astype(np.float64) ** 2
        + kz.astype(np.float64) ** 2
    )
    if np.any(wave_squared == 0.0):
        raise ValueError("a parent box contains the zero mode")
    dot = (
        kx * polarization[0]
        + ky * polarization[1]
        + kz * polarization[2]
    )
    projected = np.empty((3, axis.size, axis.size, axis.size), dtype=np.complex128)
    projected[0] = polarization[0] - kx * dot / wave_squared
    projected[1] = polarization[1] - ky * dot / wave_squared
    projected[2] = polarization[2] - kz * dot / wave_squared
    projected_norm = np.sqrt(np.sum(np.abs(projected) ** 2, axis=0))
    if np.any(projected_norm <= np.finfo(np.float64).tiny):
        raise ValueError("a projected polarization vanished")
    projected /= projected_norm[None, ...]

    slope = np.asarray(config.phase_slope, dtype=np.float64)
    phase = np.exp(
        1.0j
        * (
            phase_offset
            + slope[0] * ox
            + slope[1] * oy
            + slope[2] * oz
        )
    )
    weight = (
        (1.0 - np.abs(ox) / float(width))
        * (1.0 - np.abs(oy) / float(width))
        * (1.0 - np.abs(oz) / float(width))
    )
    projected *= (normalization * weight * phase)[None, ...]
    return projected


def _parent_arrays(
    config: MesoscopicCloudConfig,
) -> dict[str, ComplexArray]:
    normalization = _normalization(config)
    positive_p = _positive_species_array(config, "p", normalization)
    positive_q = _positive_species_array(config, "q", normalization)
    return {
        "p+": positive_p,
        "p-": np.conjugate(positive_p[:, ::-1, ::-1, ::-1]),
        "q+": positive_q,
        "q-": np.conjugate(positive_q[:, ::-1, ::-1, ::-1]),
    }


@dataclass(frozen=True)
class _TagBox:
    center: Wave
    half_width: int


@dataclass
class _Block:
    center: Wave
    half_width: int
    name: str
    interaction_class: str
    data: ComplexArray


def _carrier_groups(
    config: MesoscopicCloudConfig,
) -> dict[Wave, list[tuple[str, str]]]:
    groups: dict[Wave, list[tuple[str, str]]] = {}
    for left in _SPECIES:
        for right in _SPECIES:
            center_direction = _add(
                _species_direction(config, left),
                _species_direction(config, right),
            )
            center = _scale(config.base_scale, center_direction)
            groups.setdefault(center, []).append((left, right))
    return groups


def _unordered_pairs(
    ordered: Sequence[tuple[str, str]],
) -> list[tuple[str, str]]:
    remaining = set(ordered)
    result: list[tuple[str, str]] = []
    while remaining:
        left, right = min(remaining)
        remaining.remove((left, right))
        if left != right:
            reverse = (right, left)
            if reverse not in remaining:
                raise AssertionError("carrier group is not closed under pair reversal")
            remaining.remove(reverse)
        result.append((left, right))
    return result


def _interaction_class(pairs: Sequence[tuple[str, str]]) -> str:
    classes = {
        "self" if _species_family(left) == _species_family(right) else "cross"
        for left, right in pairs
    }
    if len(classes) != 1:
        raise AssertionError("one carrier block mixes self and cross interactions")
    return next(iter(classes))


def _fft_components(array: ComplexArray, padded_length: int) -> list[ComplexArray]:
    return [
        np.asarray(
            np.fft.fftn(array[component], s=(padded_length,) * 3, axes=(0, 1, 2)),
            dtype=np.complex128,
        )
        for component in range(3)
    ]


def _accumulate_unordered_pair(
    raw: ComplexArray,
    center: Wave,
    left: str,
    right: str,
    parents: Mapping[str, ComplexArray],
    padded_length: int,
    residual_axis: npt.NDArray[np.int64],
) -> None:
    left_fft = _fft_components(parents[left], padded_length)
    right_fft = left_fft if left == right else _fft_components(
        parents[right], padded_length
    )
    k_axes = (
        center[0] + residual_axis[:, None, None],
        center[1] + residual_axis[None, :, None],
        center[2] + residual_axis[None, None, :],
    )
    for first in range(3):
        for second in range(first, 3):
            if left == right:
                spectrum = left_fft[first] * left_fft[second]
            else:
                spectrum = (
                    left_fft[first] * right_fft[second]
                    + right_fft[first] * left_fft[second]
                )
            convolution = np.asarray(
                np.fft.ifftn(spectrum), dtype=np.complex128
            )
            # Since every parent coefficient is transverse to its own wave,
            #   sum_(a+b=k) (u(a).b) u_j(b)
            # = k_m sum_(a+b=k) u_m(a) u_j(b).
            raw[second] += 1.0j * k_axes[first] * convolution
            if first != second:
                raw[first] += 1.0j * k_axes[second] * convolution
            del convolution, spectrum
    del left_fft
    if right != left:
        del right_fft


def _project_leray_in_place(
    raw: ComplexArray,
    center: Wave,
    residual_axis: npt.NDArray[np.int64],
) -> None:
    kx = center[0] + residual_axis[:, None, None]
    ky = center[1] + residual_axis[None, :, None]
    kz = center[2] + residual_axis[None, None, :]
    wave_squared = (
        kx.astype(np.float64) ** 2
        + ky.astype(np.float64) ** 2
        + kz.astype(np.float64) ** 2
    )
    factor = raw[0] * kx + raw[1] * ky + raw[2] * kz
    np.divide(factor, wave_squared, out=factor, where=wave_squared != 0.0)
    factor[wave_squared == 0.0] = 0.0
    raw[0] -= kx * factor
    raw[1] -= ky * factor
    raw[2] -= kz * factor
    if all(abs(value) <= residual_axis[-1] for value in center):
        index = tuple(int(-value + residual_axis[-1]) for value in center)
        if all(0 <= value < residual_axis.size for value in index):
            raw[(slice(None), *index)] = 0.0


def _build_output_blocks(
    config: MesoscopicCloudConfig,
    parents: Mapping[str, ComplexArray],
) -> list[_Block]:
    d = config.width - 1
    padded_length = 4 * config.width - 3
    residual_axis = np.arange(-2 * d, 2 * d + 1, dtype=np.int64)
    blocks: list[_Block] = []
    for number, (center, ordered_pairs) in enumerate(
        sorted(_carrier_groups(config).items())
    ):
        raw = np.zeros((3, padded_length, padded_length, padded_length), dtype=np.complex128)
        for left, right in _unordered_pairs(ordered_pairs):
            _accumulate_unordered_pair(
                raw,
                center,
                left,
                right,
                parents,
                padded_length,
                residual_axis,
            )
        _project_leray_in_place(raw, center, residual_axis)
        blocks.append(
            _Block(
                center=center,
                half_width=2 * d,
                name=f"carrier_{center[0]}_{center[1]}_{center[2]}",
                interaction_class=_interaction_class(ordered_pairs),
                data=raw,
            )
        )
    return blocks


def _parent_blocks(
    config: MesoscopicCloudConfig,
    parents: Mapping[str, ComplexArray],
) -> list[_Block]:
    d = config.width - 1
    return [
        _Block(
            center=_scale(config.base_scale, _species_direction(config, species)),
            half_width=d,
            name=species,
            interaction_class="parent",
            data=array,
        )
        for species, array in parents.items()
    ]


def _block_bounds(block: _Block) -> tuple[Wave, Wave]:
    lower = tuple(value - block.half_width for value in block.center)
    upper = tuple(value + block.half_width for value in block.center)
    return lower, upper  # type: ignore[return-value]


def _covers(block: _Block, lower: Wave, upper: Wave) -> bool:
    block_lower, block_upper = _block_bounds(block)
    return all(
        block_lower[axis] <= lower[axis] and upper[axis] <= block_upper[axis]
        for axis in range(3)
    )


def _read_block_region(block: _Block, lower: Wave, upper: Wave) -> ComplexArray:
    if not _covers(block, lower, upper):
        raise ValueError("requested region is not contained in block")
    block_lower, _ = _block_bounds(block)
    slices = tuple(
        slice(lower[axis] - block_lower[axis], upper[axis] - block_lower[axis] + 1)
        for axis in range(3)
    )
    return np.asarray(block.data[(slice(None), *slices)])


def _field_on_region(
    blocks: Sequence[_Block], lower: Wave, upper: Wave
) -> ComplexArray:
    shape = tuple(upper[axis] - lower[axis] + 1 for axis in range(3))
    field = np.zeros((3, *shape), dtype=np.complex128)
    for block in blocks:
        if _covers(block, lower, upper):
            field += _read_block_region(block, lower, upper)
    return field


def _partition_cells(blocks: Sequence[_Block]) -> Iterator[tuple[Wave, Wave]]:
    edges: list[list[int]] = []
    for axis in range(3):
        values: set[int] = set()
        for block in blocks:
            lower, upper = _block_bounds(block)
            values.add(lower[axis])
            values.add(upper[axis] + 1)
        edges.append(sorted(values))
    intervals = [
        [(values[index], values[index + 1] - 1) for index in range(len(values) - 1)]
        for values in edges
    ]
    for x_interval, y_interval, z_interval in itertools.product(*intervals):
        lower = (x_interval[0], y_interval[0], z_interval[0])
        upper = (x_interval[1], y_interval[1], z_interval[1])
        if any(_covers(block, lower, upper) for block in blocks):
            yield lower, upper


def _mask_for_boxes(lower: Wave, upper: Wave, boxes: Sequence[_TagBox]) -> npt.NDArray[np.bool_]:
    shape = tuple(upper[axis] - lower[axis] + 1 for axis in range(3))
    mask = np.zeros(shape, dtype=bool)
    x = np.arange(lower[0], upper[0] + 1, dtype=np.int64)[:, None, None]
    y = np.arange(lower[1], upper[1] + 1, dtype=np.int64)[None, :, None]
    z = np.arange(lower[2], upper[2] + 1, dtype=np.int64)[None, None, :]
    for box in boxes:
        mask |= (
            (np.abs(x - box.center[0]) <= box.half_width)
            & (np.abs(y - box.center[1]) <= box.half_width)
            & (np.abs(z - box.center[2]) <= box.half_width)
        )
    return mask


def _tag_boxes(config: MesoscopicCloudConfig) -> dict[str, list[_TagBox]]:
    n = config.base_scale
    d = config.width - 1
    p = tuple(config.p_direction)
    q = tuple(config.q_direction)
    child = _add(p, q)  # type: ignore[arg-type]
    difference = tuple(a - b for a, b in zip(p, q))
    return {
        "child_core": [
            _TagBox(_scale(n, child), d),
            _TagBox(_scale(-n, child), d),
        ],
        "child_full": [
            _TagBox(_scale(n, child), 2 * d),
            _TagBox(_scale(-n, child), 2 * d),
        ],
        "difference_full": [
            _TagBox(_scale(n, difference), 2 * d),
            _TagBox(_scale(-n, difference), 2 * d),
        ],
        "difference_core": [
            _TagBox(_scale(n, difference), d),
            _TagBox(_scale(-n, difference), d),
        ],
        "harmonic_full": [
            _TagBox(_scale(2 * n, p), 2 * d),
            _TagBox(_scale(-2 * n, p), 2 * d),
            _TagBox(_scale(2 * n, q), 2 * d),
            _TagBox(_scale(-2 * n, q), 2 * d),
        ],
        "harmonic_core": [
            _TagBox(_scale(2 * n, p), d),
            _TagBox(_scale(-2 * n, p), d),
            _TagBox(_scale(2 * n, q), d),
            _TagBox(_scale(-2 * n, q), d),
        ],
        "low_full": [_TagBox((0, 0, 0), 2 * d)],
        "low_core": [_TagBox((0, 0, 0), d)],
    }


def _box_intersection(left: _TagBox, right: _TagBox) -> bool:
    return all(
        abs(a - b) <= left.half_width + right.half_width
        for a, b in zip(left.center, right.center)
    )


def _channel_tags_are_disjoint(tags: Mapping[str, Sequence[_TagBox]]) -> bool:
    selected = (
        "child_full",
        "difference_full",
        "harmonic_full",
        "low_full",
    )
    for first_index, first_name in enumerate(selected):
        for second_name in selected[first_index + 1 :]:
            if any(
                _box_intersection(left, right)
                for left in tags[first_name]
                for right in tags[second_name]
            ):
                return False
    return True


@dataclass
class _GlobalSums:
    unique_output_modes: int = 0
    full_squared: float = 0.0
    child_core_squared: float = 0.0
    child_full_squared: float = 0.0
    child_heat_squared: float = 0.0
    difference_core_squared: float = 0.0
    difference_full_squared: float = 0.0
    harmonic_core_squared: float = 0.0
    harmonic_full_squared: float = 0.0
    low_core_squared: float = 0.0
    low_full_squared: float = 0.0
    self_squared: float = 0.0
    cross_squared: float = 0.0
    self_cross_real_pairing: float = 0.0
    energy_pairing: complex = 0.0j
    divergence_squared: float = 0.0
    divergence_denominator_squared: float = 0.0
    reality_defect_squared: float = 0.0


def _global_sums(
    config: MesoscopicCloudConfig,
    output_blocks: Sequence[_Block],
    parent_blocks: Sequence[_Block],
) -> _GlobalSums:
    result = _GlobalSums()
    tags = _tag_boxes(config)
    partition_blocks = list(output_blocks) + list(parent_blocks)
    self_blocks = [b for b in output_blocks if b.interaction_class == "self"]
    cross_blocks = [b for b in output_blocks if b.interaction_class == "cross"]
    for lower, upper in _partition_cells(partition_blocks):
        active_output = [b for b in output_blocks if _covers(b, lower, upper)]
        if not active_output:
            continue
        field = _field_on_region(active_output, lower, upper)
        density = np.sum(np.abs(field) ** 2, axis=0)
        result.unique_output_modes += int(density.size)
        result.full_squared += float(np.sum(density))

        masks = {
            name: _mask_for_boxes(lower, upper, boxes)
            for name, boxes in tags.items()
        }
        result.child_core_squared += float(np.sum(density[masks["child_core"]]))
        result.child_full_squared += float(np.sum(density[masks["child_full"]]))
        result.difference_full_squared += float(
            np.sum(density[masks["difference_full"]])
        )
        result.difference_core_squared += float(
            np.sum(density[masks["difference_core"]])
        )
        result.harmonic_full_squared += float(
            np.sum(density[masks["harmonic_full"]])
        )
        result.harmonic_core_squared += float(
            np.sum(density[masks["harmonic_core"]])
        )
        result.low_full_squared += float(np.sum(density[masks["low_full"]]))
        result.low_core_squared += float(np.sum(density[masks["low_core"]]))

        x = np.arange(lower[0], upper[0] + 1, dtype=np.float64)[:, None, None]
        y = np.arange(lower[1], upper[1] + 1, dtype=np.float64)[None, :, None]
        z = np.arange(lower[2], upper[2] + 1, dtype=np.float64)[None, None, :]
        wave_squared = x * x + y * y + z * z
        scaled_squared = wave_squared / float(config.base_scale**2)
        argument = config.viscosity * config.tau * scaled_squared
        heat = np.empty_like(scaled_squared)
        nonzero = scaled_squared != 0.0
        heat[nonzero] = -np.expm1(-argument[nonzero]) / (
            config.viscosity * scaled_squared[nonzero]
        )
        heat[~nonzero] = config.tau
        child_mask = masks["child_core"]
        result.child_heat_squared += float(
            np.sum((heat[child_mask] ** 2) * density[child_mask])
        )

        divergence = x * field[0] + y * field[1] + z * field[2]
        result.divergence_squared += float(np.sum(np.abs(divergence) ** 2))
        result.divergence_denominator_squared += float(
            np.sum(wave_squared * density)
        )

        negative_lower = tuple(-value for value in upper)
        negative_upper = tuple(-value for value in lower)
        negative_active = [
            b for b in output_blocks if _covers(b, negative_lower, negative_upper)
        ]
        negative_field = _field_on_region(
            negative_active, negative_lower, negative_upper
        )[:, ::-1, ::-1, ::-1]
        reality_defect = field - np.conjugate(negative_field)
        result.reality_defect_squared += float(
            np.sum(np.abs(reality_defect) ** 2)
        )

        active_parent = [b for b in parent_blocks if _covers(b, lower, upper)]
        if active_parent:
            parent = _field_on_region(active_parent, lower, upper)
            result.energy_pairing += complex(np.vdot(parent, field))

        active_self = [b for b in self_blocks if _covers(b, lower, upper)]
        active_cross = [b for b in cross_blocks if _covers(b, lower, upper)]
        self_field = _field_on_region(active_self, lower, upper)
        cross_field = _field_on_region(active_cross, lower, upper)
        result.self_squared += float(np.sum(np.abs(self_field) ** 2))
        result.cross_squared += float(np.sum(np.abs(cross_field) ** 2))
        result.self_cross_real_pairing += float(np.vdot(self_field, cross_field).real)
    return result


@dataclass(frozen=True)
class _FieldResiduals:
    norm_squared: float
    divergence_relative: float
    reality_relative: float
    coefficient_l1: float
    minimum_wavenumber: float
    maximum_wavenumber: float
    gradient_norm_squared: float


def _parent_residuals(blocks: Sequence[_Block]) -> _FieldResiduals:
    norm_squared = 0.0
    divergence_squared = 0.0
    divergence_denominator = 0.0
    reality_squared = 0.0
    coefficient_l1 = 0.0
    minimum_wavenumber = math.inf
    maximum_wavenumber = 0.0
    for lower, upper in _partition_cells(blocks):
        active = [b for b in blocks if _covers(b, lower, upper)]
        field = _field_on_region(active, lower, upper)
        density = np.sum(np.abs(field) ** 2, axis=0)
        norm_squared += float(np.sum(density))
        x = np.arange(lower[0], upper[0] + 1, dtype=np.float64)[:, None, None]
        y = np.arange(lower[1], upper[1] + 1, dtype=np.float64)[None, :, None]
        z = np.arange(lower[2], upper[2] + 1, dtype=np.float64)[None, None, :]
        wave_squared = x * x + y * y + z * z
        coefficient_l1 += float(np.sum(np.sqrt(density)))
        radii = np.sqrt(wave_squared)
        minimum_wavenumber = min(minimum_wavenumber, float(np.min(radii)))
        maximum_wavenumber = max(maximum_wavenumber, float(np.max(radii)))
        divergence = x * field[0] + y * field[1] + z * field[2]
        divergence_squared += float(np.sum(np.abs(divergence) ** 2))
        divergence_denominator += float(np.sum(wave_squared * density))
        neg_lower = tuple(-value for value in upper)
        neg_upper = tuple(-value for value in lower)
        neg_active = [b for b in blocks if _covers(b, neg_lower, neg_upper)]
        reflected = _field_on_region(neg_active, neg_lower, neg_upper)[
            :, ::-1, ::-1, ::-1
        ]
        reality_squared += float(
            np.sum(np.abs(field - np.conjugate(reflected)) ** 2)
        )
    return _FieldResiduals(
        norm_squared=norm_squared,
        divergence_relative=math.sqrt(
            divergence_squared / max(divergence_denominator, np.finfo(float).tiny)
        ),
        reality_relative=math.sqrt(
            reality_squared / max(norm_squared, np.finfo(float).tiny)
        ),
        coefficient_l1=coefficient_l1,
        minimum_wavenumber=minimum_wavenumber,
        maximum_wavenumber=maximum_wavenumber,
        gradient_norm_squared=divergence_denominator,
    )


def _parent_boxes_disjoint(blocks: Sequence[_Block]) -> bool:
    for left, right in itertools.combinations(blocks, 2):
        left_lower, left_upper = _block_bounds(left)
        right_lower, right_upper = _block_bounds(right)
        if all(
            max(left_lower[a], right_lower[a]) <= min(left_upper[a], right_upper[a])
            for a in range(3)
        ):
            return False
    return True


def _output_boxes_overlap(blocks: Sequence[_Block]) -> bool:
    for left, right in itertools.combinations(blocks, 2):
        left_lower, left_upper = _block_bounds(left)
        right_lower, right_upper = _block_bounds(right)
        if all(
            max(left_lower[a], right_lower[a]) <= min(left_upper[a], right_upper[a])
            for a in range(3)
        ):
            return True
    return False


@dataclass(frozen=True)
class LocalFFTCloudDiagnostics:
    """Full finite-convolution diagnostics for one cloud."""

    base_scale: int
    gamma: float
    width: int
    half_width: int
    input_side_length: int
    padded_side_length: int
    zero_padding_is_alias_free: bool
    parent_mode_count: int
    output_block_count: int
    unique_output_mode_count: int
    parent_boxes_disjoint: bool
    output_boxes_overlap: bool
    channel_tags_are_orthogonal: bool
    target_parent_norm_squared: float
    parent_norm_squared: float
    parent_effective_mode_count: float
    parent_min_wavenumber: float
    parent_max_wavenumber: float
    parent_rms_wavenumber: float
    normalization_relative_error: float
    parent_divergence_relative: float
    parent_reality_relative: float
    full_nonlinear_norm_squared: float
    full_nonlinear_norm: float
    child_core_norm_squared: float
    child_core_norm: float
    child_full_norm_squared: float
    child_forcing_fraction: float
    off_child_norm_squared: float
    off_child_main_ratio: float
    difference_core_norm_squared: float
    difference_full_norm_squared: float
    difference_sideband_main_ratio: float
    harmonic_core_norm_squared: float
    harmonic_full_norm_squared: float
    low_core_norm_squared: float
    low_full_norm_squared: float
    low_side_forcing_ratio: float
    self_interaction_norm_squared: float
    cross_interaction_norm_squared: float
    self_interaction_cross_interaction_ratio: float
    self_cross_real_pairing: float
    energy_pairing_real: float
    energy_pairing_imag: float
    energy_cancellation_residual: float
    nonlinear_divergence_relative: float
    nonlinear_reality_relative: float
    a_critical: float
    a_unit: float
    gain_g: float
    heat_factor_h: float
    duhamel_energy_ratio: float
    duhamel_identity_rhs: float
    duhamel_identity_relative_error: float
    estimated_peak_working_bytes: int
    estimated_scratch_bytes: int
    elapsed_seconds: float
    parseval_convention: str
    convolution_convention: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _resource_estimate(width: int) -> int:
    padded = 4 * width - 3
    scalar = padded**3 * np.dtype(np.complex128).itemsize
    # Nine retained three-component output blocks cost 27 scalars.  Six cached
    # transforms, the current raw block, parent arrays, and generous FFT/cell
    # temporaries are bounded by another 16 padded scalars.
    return 43 * scalar


def _validate_resources(
    config: MesoscopicCloudConfig,
    maximum_working_bytes: int,
) -> int:
    if (
        isinstance(maximum_working_bytes, bool)
        or not isinstance(maximum_working_bytes, int)
        or maximum_working_bytes <= 0
    ):
        raise ValueError("maximum_working_bytes must be positive")
    working = _resource_estimate(config.width)
    if working > maximum_working_bytes:
        raise ValueError(
            f"estimated working memory {working} exceeds limit {maximum_working_bytes}"
        )
    return working


def _run_with_blocks(
    config: MesoscopicCloudConfig,
    *,
    maximum_working_bytes: int,
) -> tuple[
    dict[str, ComplexArray],
    list[_Block],
    list[_Block],
    int,
]:
    working = _validate_resources(config, maximum_working_bytes)
    parents = _parent_arrays(config)
    parent_blocks = _parent_blocks(config, parents)
    if not _parent_boxes_disjoint(parent_blocks):
        raise ValueError("parent boxes overlap; local species normalization is ambiguous")
    outputs = _build_output_blocks(config, parents)
    return parents, parent_blocks, outputs, working


def measure_local_fft_cloud(
    config: MesoscopicCloudConfig,
    *,
    maximum_working_bytes: int = 1_500_000_000,
) -> LocalFFTCloudDiagnostics:
    """Compute all local convolutions and global diagnostics without aliasing."""

    started = time.perf_counter()
    parents, parent_blocks, outputs, working = _run_with_blocks(
        config,
        maximum_working_bytes=maximum_working_bytes,
    )
    try:
        parent = _parent_residuals(parent_blocks)
        sums = _global_sums(config, outputs, parent_blocks)
        target_parent = 2.0 * config.energy_constant / config.base_scale
        child = sums.child_core_squared
        full = sums.full_squared
        if not child > 0.0 or not full > 0.0:
            raise ValueError("nonlinear child or full norm vanished")
        off_child = max(full - child, 0.0)
        cancellation_denominator = math.sqrt(parent.norm_squared * full)
        a_critical = math.sqrt(child)
        a_unit = a_critical / parent.norm_squared
        gain = a_critical / (config.base_scale**2 * parent.norm_squared)
        heat_factor = sums.child_heat_squared / child
        duhamel = (
            sums.child_heat_squared
            / config.base_scale**4
            / parent.norm_squared
        )
        identity_rhs = heat_factor * gain * gain * parent.norm_squared
        identity_error = abs(duhamel - identity_rhs) / max(
            abs(duhamel), np.finfo(float).tiny
        )
        diagnostics = LocalFFTCloudDiagnostics(
            base_scale=config.base_scale,
            gamma=config.gamma,
            width=config.width,
            half_width=config.width - 1,
            input_side_length=2 * config.width - 1,
            padded_side_length=4 * config.width - 3,
            zero_padding_is_alias_free=True,
            parent_mode_count=4 * (2 * config.width - 1) ** 3,
            output_block_count=len(outputs),
            unique_output_mode_count=sums.unique_output_modes,
            parent_boxes_disjoint=True,
            output_boxes_overlap=_output_boxes_overlap(outputs),
            channel_tags_are_orthogonal=_channel_tags_are_disjoint(_tag_boxes(config)),
            target_parent_norm_squared=target_parent,
            parent_norm_squared=parent.norm_squared,
            parent_effective_mode_count=(
                parent.coefficient_l1**2 / parent.norm_squared
            ),
            parent_min_wavenumber=parent.minimum_wavenumber,
            parent_max_wavenumber=parent.maximum_wavenumber,
            parent_rms_wavenumber=math.sqrt(
                parent.gradient_norm_squared / parent.norm_squared
            ),
            normalization_relative_error=abs(parent.norm_squared - target_parent)
            / target_parent,
            parent_divergence_relative=parent.divergence_relative,
            parent_reality_relative=parent.reality_relative,
            full_nonlinear_norm_squared=full,
            full_nonlinear_norm=math.sqrt(full),
            child_core_norm_squared=child,
            child_core_norm=a_critical,
            child_full_norm_squared=sums.child_full_squared,
            child_forcing_fraction=child / full,
            off_child_norm_squared=off_child,
            off_child_main_ratio=math.sqrt(off_child / child),
            difference_core_norm_squared=sums.difference_core_squared,
            difference_full_norm_squared=sums.difference_full_squared,
            difference_sideband_main_ratio=math.sqrt(
                sums.difference_full_squared / child
            ),
            harmonic_core_norm_squared=sums.harmonic_core_squared,
            harmonic_full_norm_squared=sums.harmonic_full_squared,
            low_core_norm_squared=sums.low_core_squared,
            low_full_norm_squared=sums.low_full_squared,
            low_side_forcing_ratio=math.sqrt(sums.low_full_squared / child),
            self_interaction_norm_squared=sums.self_squared,
            cross_interaction_norm_squared=sums.cross_squared,
            self_interaction_cross_interaction_ratio=math.sqrt(
                sums.self_squared / sums.cross_squared
            ),
            self_cross_real_pairing=sums.self_cross_real_pairing,
            energy_pairing_real=float(sums.energy_pairing.real),
            energy_pairing_imag=float(sums.energy_pairing.imag),
            energy_cancellation_residual=abs(sums.energy_pairing.real)
            / max(cancellation_denominator, np.finfo(float).tiny),
            nonlinear_divergence_relative=math.sqrt(
                sums.divergence_squared
                / max(sums.divergence_denominator_squared, np.finfo(float).tiny)
            ),
            nonlinear_reality_relative=math.sqrt(
                sums.reality_defect_squared / full
            ),
            a_critical=a_critical,
            a_unit=a_unit,
            gain_g=gain,
            heat_factor_h=heat_factor,
            duhamel_energy_ratio=duhamel,
            duhamel_identity_rhs=identity_rhs,
            duhamel_identity_relative_error=identity_error,
            estimated_peak_working_bytes=working,
            estimated_scratch_bytes=0,
            elapsed_seconds=time.perf_counter() - started,
            parseval_convention=(
                "normalized-torus Fourier series: ||f||_2^2=sum_k |f_hat(k)|^2"
            ),
            convolution_convention=(
                "NumPy backward FFT; input side L=2W-1 zero-padded to "
                "K=2L-1, so ifft(fft(a)fft(b)) is linear convolution"
            ),
        )
        return diagnostics
    finally:
        del parents, parent_blocks, outputs


def local_fft_leray_coefficients(
    config: MesoscopicCloudConfig,
    *,
    maximum_output_modes: int = 2_000_000,
    maximum_working_bytes: int = 1_500_000_000,
) -> dict[Wave, ComplexArray]:
    """Materialize the combined global coefficient map for small cross-checks."""

    if (
        isinstance(maximum_output_modes, bool)
        or not isinstance(maximum_output_modes, int)
        or maximum_output_modes < 1
    ):
        raise ValueError("maximum_output_modes must be a positive integer")
    width = config.width
    conservative = 9 * (4 * width - 3) ** 3
    if conservative > maximum_output_modes:
        raise ValueError(
            f"conservative output count {conservative} exceeds maximum_output_modes"
        )
    parents, parent_blocks, outputs, _ = _run_with_blocks(
        config,
        maximum_working_bytes=maximum_working_bytes,
    )
    coefficients: dict[Wave, ComplexArray] = {}
    for lower, upper in _partition_cells(outputs):
        active = [block for block in outputs if _covers(block, lower, upper)]
        field = _field_on_region(active, lower, upper)
        for local in np.ndindex(field.shape[1:]):
            wave = tuple(lower[axis] + local[axis] for axis in range(3))
            coefficients[wave] = np.array(
                field[(slice(None), *local)], copy=True
            )
    del parents, parent_blocks, outputs
    return coefficients


def measure_local_fft_mesoscopic_metrics(
    config: MesoscopicCloudConfig,
    *,
    maximum_working_bytes: int = 1_500_000_000,
) -> MesoscopicCloudMetrics:
    """Adapt the full FFT result to the repository's scaling-row schema.

    Unlike the target quadrature backend, this adapter always uses the actual
    globally combined nonlinear field for the full norm and every tagged
    projection.  Tag ratios are explicitly marked non-orthogonal when their
    boxes overlap.  The self/cross ratio is nevertheless available in that
    case because it is computed from the two algebraic interaction fields,
    not by pretending overlapping geometric tags are disjoint.
    """

    (
        overlaps,
        parent_disjoint,
        child_isolated,
        channels_disjoint,
        full_outputs_disjoint,
    ) = _scaling._overlap_diagnostics(config)
    if not parent_disjoint:
        # The old path returns an honest unavailable row without attempting a
        # convolution or assigning coefficients twice in an overlapping
        # parent support.
        return _scaling.measure_mesoscopic_cloud(config)

    exact = measure_local_fft_cloud(
        config, maximum_working_bytes=maximum_working_bytes
    )
    mode_count = exact.parent_mode_count
    effective_modes = exact.parent_effective_mode_count
    a_unit = exact.a_unit
    chi_support = a_unit / (config.base_scale * math.sqrt(mode_count))
    chi_effective = a_unit / (
        config.base_scale * math.sqrt(effective_modes)
    )

    def spill(full: float, core: float) -> float:
        return max(full - core, 0.0)

    channel_squared = {
        "child": exact.child_core_norm_squared,
        "child_spill": spill(
            exact.child_full_norm_squared, exact.child_core_norm_squared
        ),
        "difference": exact.difference_core_norm_squared,
        "difference_spill": spill(
            exact.difference_full_norm_squared,
            exact.difference_core_norm_squared,
        ),
        "harmonic": exact.harmonic_core_norm_squared,
        "same_carrier_spill": spill(
            exact.harmonic_full_norm_squared,
            exact.harmonic_core_norm_squared,
        ),
        "low": exact.low_core_norm_squared,
        "low_spill": spill(
            exact.low_full_norm_squared, exact.low_core_norm_squared
        ),
    }
    tagged_total = sum(channel_squared.values())
    channel_norms = {
        name: math.sqrt(value) for name, value in channel_squared.items()
    }
    channel_ratios = {
        name: value / tagged_total for name, value in channel_squared.items()
    }
    channel_ratios["off_child"] = 1.0 - channel_ratios["child"]

    support_radius_factor = exact.parent_max_wavenumber / config.base_scale
    bound_prefactor = (
        2.0
        * support_radius_factor**2
        * config.tau**2
        * config.energy_constant
        / config.base_scale**3
    )
    universal_bound = bound_prefactor * mode_count
    effective_bound = bound_prefactor * effective_modes
    if exact.duhamel_energy_ratio > effective_bound * (
        1.0 + 65536.0 * np.finfo(float).eps
    ):
        raise AssertionError("exact Duhamel ratio exceeds effective-count bound")

    finite_reasons: list[str] = []
    if not child_isolated:
        finite_reasons.append("child_band_not_isolated")
    if not channels_disjoint:
        finite_reasons.append("channel_bands_overlap")
    if exact.duhamel_energy_ratio < config.required_duhamel_ratio:
        finite_reasons.append("duhamel_ratio_below_required_target")
    asymptotic_reasons = list(finite_reasons)
    if config.gamma < 1.0:
        asymptotic_reasons.append(
            "generic_duhamel_ratio_decays_for_gamma_below_one"
        )
    p = np.asarray(config.p_direction, dtype=np.float64)
    q = np.asarray(config.q_direction, dtype=np.float64)
    return MesoscopicCloudMetrics(
        base_scale=config.base_scale,
        gamma=config.gamma,
        width=config.width,
        half_width=config.width - 1,
        width_source=(
            "override" if config.width_override is not None else "power_floor"
        ),
        relative_width=config.width / config.base_scale,
        mode_count=mode_count,
        effective_mode_count=effective_modes,
        angle_sine=float(
            np.linalg.norm(np.cross(p, q)) / (np.linalg.norm(p) * np.linalg.norm(q))
        ),
        box_overlap_counts=overlaps,
        parent_boxes_disjoint=True,
        child_band_isolated=child_isolated,
        channel_bands_disjoint=channels_disjoint,
        full_output_bands_disjoint=full_outputs_disjoint,
        measurement_available=True,
        measurement_method="exact-zero-padded-local-fft-global-combination",
        targets_per_band=exact.padded_side_length**3,
        critical_energy=config.energy_constant / config.base_scale,
        parent_norm_squared=exact.parent_norm_squared,
        normalization_relative_error=exact.normalization_relative_error,
        divergence_relative=exact.parent_divergence_relative,
        reality_relative=exact.parent_reality_relative,
        parent_divergence_relative=exact.parent_divergence_relative,
        parent_reality_relative=exact.parent_reality_relative,
        nonlinear_divergence_relative=exact.nonlinear_divergence_relative,
        nonlinear_reality_relative=exact.nonlinear_reality_relative,
        parent_min_wavenumber=exact.parent_min_wavenumber,
        parent_max_wavenumber=exact.parent_max_wavenumber,
        parent_rms_wavenumber=exact.parent_rms_wavenumber,
        duhamel_upper_bound_method=(
            "phase-independent-bernstein-heat-support-and-effective-count"
        ),
        support_radius_factor_kappa=support_radius_factor,
        universal_duhamel_upper_bound=universal_bound,
        effective_duhamel_upper_bound=effective_bound,
        duhamel_to_universal_upper_bound=(
            exact.duhamel_energy_ratio / universal_bound
        ),
        duhamel_to_effective_upper_bound=(
            exact.duhamel_energy_ratio / effective_bound
        ),
        a_critical=exact.a_critical,
        a_unit=exact.a_unit,
        chi_support=chi_support,
        chi_effective=chi_effective,
        gain_g=exact.gain_g,
        heat_factor_h=exact.heat_factor_h,
        duhamel_norm_ratio=math.sqrt(exact.duhamel_energy_ratio),
        duhamel_energy_ratio=exact.duhamel_energy_ratio,
        duhamel_identity_rhs=exact.duhamel_identity_rhs,
        duhamel_identity_relative_error=exact.duhamel_identity_relative_error,
        channel_norms=channel_norms,
        channel_ratios=channel_ratios,
        channel_ratios_are_orthogonal=full_outputs_disjoint,
        forcing_metric_method=(
            "exact-local-fft-full-global-overlap-aware-with-self-cross-decomposition"
        ),
        child_forcing_fraction=exact.child_forcing_fraction,
        low_side_forcing_ratio=exact.low_side_forcing_ratio,
        off_chain_main_ratio=exact.off_child_main_ratio,
        difference_sideband_main_ratio=exact.difference_sideband_main_ratio,
        self_interaction_cross_interaction_ratio=(
            exact.self_interaction_cross_interaction_ratio
        ),
        energy_cancellation_available=True,
        energy_cancellation_method="exact-zero-padded-local-fft-full-pairing",
        energy_cancellation_residual=exact.energy_cancellation_residual,
        energy_cancellation_pairing=exact.energy_pairing_real,
        full_nonlinear_norm=exact.full_nonlinear_norm,
        predicted_a_unit_exponent=1.0 + 1.5 * config.gamma,
        predicted_a_critical_exponent=1.5 * config.gamma,
        predicted_g_exponent=1.5 * config.gamma - 1.0,
        predicted_duhamel_exponent=3.0 * config.gamma - 3.0,
        predicted_g_grows=1.5 * config.gamma - 1.0 > 0.0,
        predicted_duhamel_decays=config.gamma < 1.0,
        finite_screen_pass=not finite_reasons,
        asymptotic_screen_pass=not asymptotic_reasons,
        finite_screen_reasons=tuple(finite_reasons),
        asymptotic_screen_reasons=tuple(asymptotic_reasons),
    )
