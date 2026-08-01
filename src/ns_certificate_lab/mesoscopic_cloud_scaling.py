r"""Sparse scaling screen for an angled mesoscopic Leray cloud.

This module is a falsification aid, not a singularity certificate.  It studies
two Fourier boxes centred at ``N p`` and ``N q`` with the fixed, order-one
angle

``p = (1, 1, 0)``, ``q = (1, 0, 1)``.

The positive-frequency coefficients have Fejer weights, affine phases, and a
mode-by-mode Leray polarization.  Negative frequencies are supplied by exact
Hermitian reflection.  The whole parent is then normalized to the critical
energy law ``E_N = c_E / N``.

Parent and designated child *cores* use the same half-width ``W - 1``.  The
larger half-width ``2(W - 1)`` is only the parent-parent sumset: its outer ring
is reported as child, difference, same-carrier, or low-frequency spill and is
never included in ``A_N`` or ``G_N``.

For small boxes every target coefficient is summed.  Larger boxes use a
deterministic tensor-stratified target quadrature.  A target coefficient is
still obtained by summing *all* of its parent decompositions; only the outer
sum over targets is sampled.  Thus the method does not allocate a dense
three-dimensional FFT grid and does not mistake a sampled pair sum for a
coherent convolution.

The principal identity checked here is, for frozen parent forcing and
``t = tau N^-2``,

.. math::

   {\|v(t)\|_2^2\over\|u_N\|_2^2}
   ={2c_E\over N} H_N G_N^2,
   \qquad
   G_N={A_N\over N^2\|u_N\|_2^2}.

Consequently growth of ``G_N`` alone is never an acceptance condition.  In
the optimistic ``M_N ~ N^(3 gamma)`` saturation model the displayed ratio is
``N^(3 gamma - 3)`` and decays for every ``gamma < 1``.

There is also a phase-independent finite-``N`` ceiling.  If ``M`` is the
parent support count and ``kappa=max{|k|: u_hat(k)!=0}/N``, Bernstein,
Parseval, and heat contraction give

.. math::

   D_N \le {2\kappa^2\tau^2 c_E M\over N^3}.

Replacing ``M`` by ``(||u_hat||_1/||u||_2)^2`` gives the sharper bound for the
constructed weighted field.  Both are checked against every available
measurement and do not use its phases or channel tags.

The forcing ratios use squared norms only for the child fraction.  Quantities
named ``*_ratio`` are L2-norm ratios.  A full sparse all-pairs computation also
reports the dimensionless cancellation residual
``|<u,B(u,u)>|/(||u||_2 ||B(u,u)||_2)``.  When its configured pair limit is
exceeded, the residual is ``None`` and its method records the rejected pair
count; no numerical zero is substituted.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import asdict, dataclass
from typing import Iterable, Mapping

import numpy as np
import numpy.typing as npt


FloatArray = npt.NDArray[np.float64]
ComplexArray = npt.NDArray[np.complex128]
Wave = tuple[int, int, int]


P_DIRECTION: Wave = (1, 1, 0)
Q_DIRECTION: Wave = (1, 0, 1)
P_POLARIZATION = np.asarray((0.0, 0.0, 1.0), dtype=np.complex128)
Q_POLARIZATION = np.asarray((0.0, 1.0, 0.0), dtype=np.complex128)

_SPECIES = ("p+", "p-", "q+", "q-")

__all__ = [
    "MesoscopicCloudConfig",
    "MesoscopicCloudMetrics",
    "P_DIRECTION",
    "Q_DIRECTION",
    "build_sparse_parent",
    "exact_sparse_leray_convolution",
    "measure_mesoscopic_cloud",
]


def _finite(value: float, *, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _integer_wave(value: Iterable[int], *, name: str) -> Wave:
    components = tuple(value)
    if len(components) != 3 or any(
        isinstance(component, bool) or not isinstance(component, int)
        for component in components
    ):
        raise ValueError(f"{name} must be an integer three-vector")
    return components  # type: ignore[return-value]


@dataclass(frozen=True)
class MesoscopicCloudConfig:
    """Input for one sparse mesoscopic-cloud measurement.

    ``width_override`` is useful for the fixed-relative boundary family.  For
    example, use ``gamma=1`` and ``width_override=floor(rho*N)`` to compare
    ``rho in {0.10, 0.15, 0.20}`` with the sublinear ``N**gamma`` family.
    ``gamma`` remains the claimed family exponent used by the asymptotic
    screen; the realized integer width is always recorded separately.
    """

    base_scale: int
    gamma: float
    energy_constant: float = 1.0
    viscosity: float = 1.0 / 40.0
    tau: float = 1.0
    width_override: int | None = None
    p_direction: Wave = P_DIRECTION
    q_direction: Wave = Q_DIRECTION
    phase_p: float = -math.pi / 2.0
    phase_q: float = 0.0
    phase_slope: tuple[float, float, float] = (0.173, -0.119, 0.071)
    exact_target_limit: int = 1000
    strata_per_axis: int = 4
    exact_energy_pair_limit: int = 500_000
    required_duhamel_ratio: float = 0.5

    def __post_init__(self) -> None:
        if (
            isinstance(self.base_scale, bool)
            or not isinstance(self.base_scale, int)
            or self.base_scale < 2
        ):
            raise ValueError("base_scale must be an integer at least two")
        gamma = _finite(self.gamma, name="gamma")
        if gamma <= 0.0:
            raise ValueError("gamma must be positive")
        if _finite(self.energy_constant, name="energy_constant") <= 0.0:
            raise ValueError("energy_constant must be positive")
        if _finite(self.viscosity, name="viscosity") <= 0.0:
            raise ValueError("viscosity must be positive")
        if _finite(self.tau, name="tau") <= 0.0:
            raise ValueError("tau must be positive")
        if self.width_override is not None and (
            isinstance(self.width_override, bool)
            or not isinstance(self.width_override, int)
            or self.width_override < 1
        ):
            raise ValueError("width_override must be a positive integer")
        p = np.asarray(_integer_wave(self.p_direction, name="p_direction"), dtype=float)
        q = np.asarray(_integer_wave(self.q_direction, name="q_direction"), dtype=float)
        if np.linalg.norm(p) == 0.0 or np.linalg.norm(q) == 0.0:
            raise ValueError("carrier directions must be nonzero")
        sine = float(np.linalg.norm(np.cross(p, q)) / (np.linalg.norm(p) * np.linalg.norm(q)))
        if sine < 0.25:
            raise ValueError("p_direction and q_direction must have an order-one angle")
        if np.array_equal(p + q, np.zeros(3)) or np.array_equal(p - q, np.zeros(3)):
            raise ValueError("sum and difference carrier directions must be nonzero")
        slope = tuple(_finite(value, name="phase_slope") for value in self.phase_slope)
        if len(slope) != 3:
            raise ValueError("phase_slope must be a three-vector")
        _finite(self.phase_p, name="phase_p")
        _finite(self.phase_q, name="phase_q")
        if (
            isinstance(self.exact_target_limit, bool)
            or not isinstance(self.exact_target_limit, int)
            or self.exact_target_limit < 0
        ):
            raise ValueError("exact_target_limit must be a nonnegative integer")
        if (
            isinstance(self.strata_per_axis, bool)
            or not isinstance(self.strata_per_axis, int)
            or not 1 <= self.strata_per_axis <= 32
        ):
            raise ValueError("strata_per_axis must lie between one and 32")
        if (
            isinstance(self.exact_energy_pair_limit, bool)
            or not isinstance(self.exact_energy_pair_limit, int)
            or self.exact_energy_pair_limit < 0
        ):
            raise ValueError(
                "exact_energy_pair_limit must be a nonnegative integer"
            )
        ratio = _finite(self.required_duhamel_ratio, name="required_duhamel_ratio")
        if ratio <= 0.0:
            raise ValueError("required_duhamel_ratio must be positive")

    @property
    def width(self) -> int:
        if self.width_override is not None:
            return self.width_override
        return max(1, int(math.floor(self.base_scale**self.gamma)))


@dataclass(frozen=True)
class MesoscopicCloudMetrics:
    """Diagnostics from :func:`measure_mesoscopic_cloud`."""

    base_scale: int
    gamma: float
    width: int
    half_width: int
    width_source: str
    relative_width: float
    mode_count: int
    effective_mode_count: float | None
    angle_sine: float
    box_overlap_counts: dict[str, int]
    parent_boxes_disjoint: bool
    child_band_isolated: bool
    channel_bands_disjoint: bool
    full_output_bands_disjoint: bool
    measurement_available: bool
    measurement_method: str | None
    targets_per_band: int
    critical_energy: float
    parent_norm_squared: float | None
    normalization_relative_error: float | None
    divergence_relative: float | None
    reality_relative: float | None
    parent_divergence_relative: float | None
    parent_reality_relative: float | None
    nonlinear_divergence_relative: float | None
    nonlinear_reality_relative: float | None
    parent_min_wavenumber: float | None
    parent_max_wavenumber: float | None
    parent_rms_wavenumber: float | None
    duhamel_upper_bound_method: str
    support_radius_factor_kappa: float | None
    universal_duhamel_upper_bound: float | None
    effective_duhamel_upper_bound: float | None
    duhamel_to_universal_upper_bound: float | None
    duhamel_to_effective_upper_bound: float | None
    a_critical: float | None
    a_unit: float | None
    chi_support: float | None
    chi_effective: float | None
    gain_g: float | None
    heat_factor_h: float | None
    duhamel_norm_ratio: float | None
    duhamel_energy_ratio: float | None
    duhamel_identity_rhs: float | None
    duhamel_identity_relative_error: float | None
    channel_norms: dict[str, float] | None
    channel_ratios: dict[str, float] | None
    channel_ratios_are_orthogonal: bool
    forcing_metric_method: str
    child_forcing_fraction: float | None
    low_side_forcing_ratio: float | None
    off_chain_main_ratio: float | None
    difference_sideband_main_ratio: float | None
    self_interaction_cross_interaction_ratio: float | None
    energy_cancellation_available: bool
    energy_cancellation_method: str
    energy_cancellation_residual: float | None
    energy_cancellation_pairing: float | None
    full_nonlinear_norm: float | None
    predicted_a_unit_exponent: float
    predicted_a_critical_exponent: float
    predicted_g_exponent: float
    predicted_duhamel_exponent: float
    predicted_g_grows: bool
    predicted_duhamel_decays: bool
    finite_screen_pass: bool
    asymptotic_screen_pass: bool
    finite_screen_reasons: tuple[str, ...]
    asymptotic_screen_reasons: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class _Box:
    name: str
    center: Wave
    half_width: int
    kind: str


def _add(left: Wave, right: Wave) -> Wave:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def _subtract(left: Wave, right: Wave) -> Wave:
    return tuple(a - b for a, b in zip(left, right))  # type: ignore[return-value]


def _scale(multiplier: int, wave: Wave) -> Wave:
    return tuple(multiplier * component for component in wave)  # type: ignore[return-value]


def _box_intersection_count(left: _Box, right: _Box) -> int:
    count = 1
    for left_center, right_center in zip(left.center, right.center):
        lower = max(
            left_center - left.half_width,
            right_center - right.half_width,
        )
        upper = min(
            left_center + left.half_width,
            right_center + right.half_width,
        )
        if lower > upper:
            return 0
        count *= upper - lower + 1
    return count


def _boxes(
    config: MesoscopicCloudConfig,
) -> tuple[list[_Box], list[_Box], list[_Box]]:
    n = config.base_scale
    d = config.width - 1
    p = _integer_wave(config.p_direction, name="p_direction")
    q = _integer_wave(config.q_direction, name="q_direction")
    c = _add(p, q)
    difference = _subtract(p, q)
    parents = [
        _Box("parent_p+", _scale(n, p), d, "parent"),
        _Box("parent_p-", _scale(-n, p), d, "parent"),
        _Box("parent_q+", _scale(n, q), d, "parent"),
        _Box("parent_q-", _scale(-n, q), d, "parent"),
    ]
    core_outputs = [
        _Box("child+", _scale(n, c), d, "child"),
        _Box("child-", _scale(-n, c), d, "child"),
        _Box("difference+", _scale(n, difference), d, "difference"),
        _Box("difference-", _scale(-n, difference), d, "difference"),
        _Box("harmonic_p+", _scale(2 * n, p), d, "harmonic"),
        _Box("harmonic_p-", _scale(-2 * n, p), d, "harmonic"),
        _Box("harmonic_q+", _scale(2 * n, q), d, "harmonic"),
        _Box("harmonic_q-", _scale(-2 * n, q), d, "harmonic"),
        _Box("low", (0, 0, 0), d, "low"),
    ]
    full_outputs = [
        _Box(f"full_{box.name}", box.center, 2 * d, box.kind)
        for box in core_outputs
    ]
    return parents, core_outputs, full_outputs


def _overlap_diagnostics(
    config: MesoscopicCloudConfig,
) -> tuple[dict[str, int], bool, bool, bool, bool]:
    parents, core_outputs, full_outputs = _boxes(config)
    overlaps: dict[str, int] = {}
    parent_disjoint = True
    for left, right in itertools.combinations(parents, 2):
        count = _box_intersection_count(left, right)
        if count:
            overlaps[f"{left.name}|{right.name}"] = count
            parent_disjoint = False

    channel_disjoint = True
    for left, right in itertools.combinations(core_outputs, 2):
        count = _box_intersection_count(left, right)
        if count:
            overlaps[f"{left.name}|{right.name}"] = count
            channel_disjoint = False

    child_isolated = True
    children = [box for box in core_outputs if box.kind == "child"]
    nonchildren = [box for box in core_outputs if box.kind != "child"]
    for child in children:
        # A child core has half-width d, but its parent-parent source reaches
        # the full half-width 2d.  Keeping that reachable source disjoint from
        # the old parent is the conservative N > 3d relay-clearance condition
        # for the default p,q geometry.
        full_child = next(
            box for box in full_outputs if box.name == f"full_{child.name}"
        )
        for parent in parents:
            count = _box_intersection_count(full_child, parent)
            if count:
                overlaps[f"{full_child.name}|{parent.name}"] = count
                child_isolated = False
        for other in nonchildren:
            if _box_intersection_count(child, other):
                child_isolated = False
    if _box_intersection_count(children[0], children[1]):
        child_isolated = False

    full_output_disjoint = True
    for left, right in itertools.combinations(full_outputs, 2):
        count = _box_intersection_count(left, right)
        if count:
            overlaps[f"{left.name}|{right.name}"] = count
            full_output_disjoint = False
    return (
        overlaps,
        parent_disjoint,
        child_isolated,
        channel_disjoint,
        full_output_disjoint,
    )


def _species_center(config: MesoscopicCloudConfig, species: str) -> Wave:
    p = _integer_wave(config.p_direction, name="p_direction")
    q = _integer_wave(config.q_direction, name="q_direction")
    if species == "p+":
        return p
    if species == "p-":
        return _scale(-1, p)
    if species == "q+":
        return q
    if species == "q-":
        return _scale(-1, q)
    raise ValueError(f"unknown species {species!r}")


def _weight(offsets: npt.NDArray[np.int64], width: int) -> FloatArray:
    return np.prod(1.0 - np.abs(offsets) / float(width), axis=1)


def _normalization_scale(config: MesoscopicCloudConfig) -> float:
    width = config.width
    one_dimensional_l2 = (2.0 * width * width + 1.0) / (3.0 * width)
    raw_norm_squared = 4.0 * one_dimensional_l2**3
    target_norm_squared = 2.0 * config.energy_constant / config.base_scale
    return math.sqrt(target_norm_squared / raw_norm_squared)


def _positive_coefficients(
    config: MesoscopicCloudConfig,
    family: str,
    offsets: npt.NDArray[np.int64],
    normalization: float,
) -> ComplexArray:
    if family == "p":
        center = np.asarray(config.p_direction, dtype=np.int64)
        base = P_POLARIZATION
        phase_offset = config.phase_p
    elif family == "q":
        center = np.asarray(config.q_direction, dtype=np.int64)
        base = Q_POLARIZATION
        phase_offset = config.phase_q
    else:
        raise ValueError(f"unknown positive family {family!r}")
    waves = config.base_scale * center[None, :] + offsets
    wave_squared = np.sum(waves * waves, axis=1).astype(np.float64)
    dot = np.sum(waves * base[None, :], axis=1)
    projected = base[None, :] - waves * (dot / wave_squared)[:, None]
    projected_norm = np.sqrt(np.sum(np.abs(projected) ** 2, axis=1))
    if np.any(projected_norm <= np.finfo(np.float64).tiny):
        raise ValueError("a carrier polarization became parallel to its wavevector")
    projected = projected / projected_norm[:, None]
    slope = np.asarray(config.phase_slope, dtype=np.float64)
    phase = np.exp(1.0j * (phase_offset + offsets @ slope))
    amplitude = normalization * _weight(offsets, config.width) * phase
    return np.asarray(amplitude[:, None] * projected, dtype=np.complex128)


def _coefficients(
    config: MesoscopicCloudConfig,
    species: str,
    offsets: npt.NDArray[np.int64],
    normalization: float,
) -> ComplexArray:
    if species == "p+":
        return _positive_coefficients(config, "p", offsets, normalization)
    if species == "q+":
        return _positive_coefficients(config, "q", offsets, normalization)
    if species == "p-":
        return np.conjugate(
            _positive_coefficients(config, "p", -offsets, normalization)
        )
    if species == "q-":
        return np.conjugate(
            _positive_coefficients(config, "q", -offsets, normalization)
        )
    raise ValueError(f"unknown species {species!r}")


def _offset_array(half_width: int) -> npt.NDArray[np.int64]:
    axis = np.arange(-half_width, half_width + 1, dtype=np.int64)
    mesh = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.stack(mesh, axis=-1).reshape(-1, 3)


def _structure_metrics(
    config: MesoscopicCloudConfig,
    normalization: float,
) -> tuple[float, float, float, float, float, float]:
    """Return norm, divergence, l1, min/max wave, and gradient norm squared."""

    d = config.width - 1
    axis = np.arange(-d, d + 1, dtype=np.int64)
    norm_squared = 0.0
    divergence_squared = 0.0
    coefficient_l1 = 0.0
    gradient_squared = 0.0
    minimum = math.inf
    maximum = 0.0
    # Slice in q_x so memory is O(W^2), not O(W^3).
    yz = np.stack(
        np.meshgrid(axis, axis, indexing="ij"), axis=-1
    ).reshape(-1, 2)
    for qx in axis:
        offsets = np.column_stack(
            (
                np.full(yz.shape[0], qx, dtype=np.int64),
                yz[:, 0],
                yz[:, 1],
            )
        )
        for family, center in (
            ("p", np.asarray(config.p_direction, dtype=np.int64)),
            ("q", np.asarray(config.q_direction, dtype=np.int64)),
        ):
            coefficients = _positive_coefficients(
                config, family, offsets, normalization
            )
            waves = config.base_scale * center[None, :] + offsets
            magnitudes_squared = np.sum(np.abs(coefficients) ** 2, axis=1)
            wave_squared = np.sum(waves * waves, axis=1).astype(np.float64)
            divergence = np.sum(waves * coefficients, axis=1)
            # Factor two supplies the Hermitian negative box.
            norm_squared += 2.0 * float(np.sum(magnitudes_squared))
            divergence_squared += 2.0 * float(np.sum(np.abs(divergence) ** 2))
            coefficient_l1 += 2.0 * float(np.sum(np.sqrt(magnitudes_squared)))
            gradient_squared += 2.0 * float(
                np.sum(wave_squared * magnitudes_squared)
            )
            radii = np.sqrt(wave_squared)
            minimum = min(minimum, float(np.min(radii)))
            maximum = max(maximum, float(np.max(radii)))
    return (
        norm_squared,
        divergence_squared,
        coefficient_l1,
        minimum,
        maximum,
        gradient_squared,
    )


def _pair_contribution(
    config: MesoscopicCloudConfig,
    left_species: str,
    right_species: str,
    target: Wave,
    normalization: float,
) -> ComplexArray:
    d = config.width - 1
    left_center = np.asarray(_species_center(config, left_species), dtype=np.int64)
    right_center = np.asarray(_species_center(config, right_species), dtype=np.int64)
    residual = np.asarray(target, dtype=np.int64) - config.base_scale * (
        left_center + right_center
    )
    if np.any(np.abs(residual) > 2 * d):
        return np.zeros(3, dtype=np.complex128)
    lower = np.maximum(-d, residual - d)
    upper = np.minimum(d, residual + d)
    if np.any(lower > upper):
        return np.zeros(3, dtype=np.complex128)
    axes = [
        np.arange(lower[index], upper[index] + 1, dtype=np.int64)
        for index in range(3)
    ]
    mesh = np.meshgrid(*axes, indexing="ij")
    left_offsets = np.stack(mesh, axis=-1).reshape(-1, 3)
    right_offsets = residual[None, :] - left_offsets
    left = _coefficients(config, left_species, left_offsets, normalization)
    right = _coefficients(config, right_species, right_offsets, normalization)
    right_waves = config.base_scale * right_center[None, :] + right_offsets
    directional = np.sum(left * right_waves, axis=1)
    return np.asarray(
        np.sum(1.0j * directional[:, None] * right, axis=0),
        dtype=np.complex128,
    )


def _target_coefficient(
    config: MesoscopicCloudConfig,
    target: Wave,
    normalization: float,
) -> ComplexArray:
    accumulated = np.zeros(3, dtype=np.complex128)
    for left_species in _SPECIES:
        for right_species in _SPECIES:
            accumulated += _pair_contribution(
                config,
                left_species,
                right_species,
                target,
                normalization,
            )
    wave = np.asarray(target, dtype=np.float64)
    wave_squared = float(np.dot(wave, wave))
    if wave_squared == 0.0:
        return np.zeros(3, dtype=np.complex128)
    return np.asarray(
        accumulated - wave * (np.dot(wave, accumulated) / wave_squared),
        dtype=np.complex128,
    )


def _target_quadrature(
    output_half_width: int,
    *,
    exact_target_limit: int,
    strata_per_axis: int,
) -> tuple[list[tuple[Wave, int]], str]:
    axis = np.arange(
        -output_half_width, output_half_width + 1, dtype=np.int64
    )
    total = int(axis.size**3)
    if total <= exact_target_limit:
        return [
            ((int(x), int(y), int(z)), 1)
            for x, y, z in itertools.product(axis, repeat=3)
        ], "exact-target"

    pieces = [piece for piece in np.array_split(axis, strata_per_axis) if piece.size]
    samples: list[tuple[Wave, int]] = []
    for x_piece, y_piece, z_piece in itertools.product(pieces, repeat=3):
        offset = (
            int(x_piece[(x_piece.size - 1) // 2]),
            int(y_piece[(y_piece.size - 1) // 2]),
            int(z_piece[(z_piece.size - 1) // 2]),
        )
        weight = int(x_piece.size * y_piece.size * z_piece.size)
        samples.append((offset, weight))
    if sum(weight for _, weight in samples) != total:
        raise AssertionError("stratified target weights do not cover the target box")
    return samples, "deterministic-stratified-target"


def _heat_multiplier(config: MesoscopicCloudConfig, target: Wave) -> float:
    scaled_squared = sum(component * component for component in target) / float(
        config.base_scale**2
    )
    if scaled_squared == 0.0:
        return config.tau
    argument = config.viscosity * config.tau * scaled_squared
    return -math.expm1(-argument) / (config.viscosity * scaled_squared)


def _measure_band(
    config: MesoscopicCloudConfig,
    center_direction: Wave,
    normalization: float,
    quadrature: list[tuple[Wave, int]],
) -> tuple[float, float]:
    norm_squared = 0.0
    heat_squared = 0.0
    center = _scale(config.base_scale, center_direction)
    for offset, weight in quadrature:
        target = _add(center, offset)
        coefficient = _target_coefficient(config, target, normalization)
        magnitude_squared = float(np.vdot(coefficient, coefficient).real)
        multiplier = _heat_multiplier(config, target)
        norm_squared += weight * magnitude_squared
        heat_squared += weight * multiplier * multiplier * magnitude_squared
    return norm_squared, heat_squared


def build_sparse_parent(
    config: MesoscopicCloudConfig,
    *,
    maximum_modes: int = 10000,
) -> dict[Wave, ComplexArray]:
    """Materialize the normalized sparse parent for an independent small check."""

    overlaps, parent_disjoint, _, _, _ = _overlap_diagnostics(config)
    if not parent_disjoint:
        parent_keys = [key for key in overlaps if key.startswith("parent_")]
        raise ValueError(f"parent boxes overlap: {parent_keys}")
    mode_count = 4 * (2 * config.width - 1) ** 3
    if mode_count > maximum_modes:
        raise ValueError("sparse parent exceeds maximum_modes")
    normalization = _normalization_scale(config)
    offsets = _offset_array(config.width - 1)
    output: dict[Wave, ComplexArray] = {}
    for species in _SPECIES:
        center = np.asarray(_species_center(config, species), dtype=np.int64)
        waves = config.base_scale * center[None, :] + offsets
        coefficients = _coefficients(config, species, offsets, normalization)
        for wave, coefficient in zip(waves, coefficients):
            key = tuple(int(component) for component in wave)
            if key in output:
                raise ValueError("parent boxes overlap after lattice construction")
            output[key] = np.asarray(coefficient, dtype=np.complex128)
    return output


def exact_sparse_leray_convolution(
    parent: Mapping[Wave, npt.ArrayLike],
    *,
    maximum_pairs: int = 4_000_000,
) -> dict[Wave, ComplexArray]:
    """Return the exact finite pair enumeration of ``P((u.grad)u)``.

    "Exact" here means that every finite Fourier pair is enumerated; the
    coefficient arithmetic remains complex binary64.
    """

    items = [
        (wave, np.asarray(coefficient, dtype=np.complex128))
        for wave, coefficient in parent.items()
    ]
    if any(len(wave) != 3 or coefficient.shape != (3,) for wave, coefficient in items):
        raise ValueError("parent must map three-dimensional waves to three-vectors")
    if len(items) ** 2 > maximum_pairs:
        raise ValueError("sparse convolution exceeds maximum_pairs")
    raw: dict[Wave, ComplexArray] = {}
    for left_wave, left in items:
        for right_wave, right in items:
            target = _add(left_wave, right_wave)
            directional = np.dot(left, np.asarray(right_wave, dtype=np.float64))
            term = 1.0j * directional * right
            if target in raw:
                raw[target] += term
            else:
                raw[target] = np.asarray(term, dtype=np.complex128)
    projected: dict[Wave, ComplexArray] = {}
    for target, coefficient in raw.items():
        wave = np.asarray(target, dtype=np.float64)
        wave_squared = float(np.dot(wave, wave))
        if wave_squared == 0.0:
            projected[target] = np.zeros(3, dtype=np.complex128)
        else:
            projected[target] = np.asarray(
                coefficient - wave * (np.dot(wave, coefficient) / wave_squared),
                dtype=np.complex128,
            )
    return projected


@dataclass(frozen=True)
class _ExactFullDiagnostic:
    available: bool
    method: str
    full_norm_squared: float | None = None
    pairing: float | None = None
    cancellation_residual: float | None = None
    child_core_squared: float | None = None
    child_full_squared: float | None = None
    difference_full_squared: float | None = None
    harmonic_full_squared: float | None = None
    low_full_squared: float | None = None


def _wave_in_box(wave: Wave, box: _Box) -> bool:
    return all(
        abs(component - center) <= box.half_width
        for component, center in zip(wave, box.center)
    )


def _norm_squared_in_boxes(
    field: Mapping[Wave, npt.ArrayLike], boxes: Iterable[_Box]
) -> float:
    selected = tuple(boxes)
    return sum(
        float(np.vdot(value, value).real)
        for wave, value in field.items()
        if any(_wave_in_box(wave, box) for box in selected)
    )


def _exact_full_diagnostic(
    config: MesoscopicCloudConfig,
) -> _ExactFullDiagnostic:
    mode_count = 4 * (2 * config.width - 1) ** 3
    pair_count = mode_count * mode_count
    if pair_count > config.exact_energy_pair_limit:
        return _ExactFullDiagnostic(
            available=False,
            method=(
                "unavailable-pair-limit: "
                f"{pair_count}>{config.exact_energy_pair_limit}"
            ),
        )
    parent = build_sparse_parent(config, maximum_modes=mode_count)
    nonlinear = exact_sparse_leray_convolution(
        parent, maximum_pairs=config.exact_energy_pair_limit
    )
    parent_norm_squared = sum(
        float(np.vdot(value, value).real) for value in parent.values()
    )
    full_norm_squared = sum(
        float(np.vdot(value, value).real) for value in nonlinear.values()
    )
    pairing = float(
        sum(
            np.vdot(value, nonlinear.get(wave, np.zeros(3))).real
            for wave, value in parent.items()
        )
    )
    denominator = math.sqrt(parent_norm_squared * full_norm_squared)
    residual = (
        abs(pairing) / denominator
        if denominator > 0.0
        else None
    )
    _, core_outputs, full_outputs = _boxes(config)
    return _ExactFullDiagnostic(
        available=True,
        method=f"exact-sparse-all-pairs:{pair_count}",
        full_norm_squared=full_norm_squared,
        pairing=pairing,
        cancellation_residual=residual,
        child_core_squared=_norm_squared_in_boxes(
            nonlinear, (box for box in core_outputs if box.kind == "child")
        ),
        child_full_squared=_norm_squared_in_boxes(
            nonlinear, (box for box in full_outputs if box.kind == "child")
        ),
        difference_full_squared=_norm_squared_in_boxes(
            nonlinear,
            (box for box in full_outputs if box.kind == "difference"),
        ),
        harmonic_full_squared=_norm_squared_in_boxes(
            nonlinear,
            (box for box in full_outputs if box.kind == "harmonic"),
        ),
        low_full_squared=_norm_squared_in_boxes(
            nonlinear, (box for box in full_outputs if box.kind == "low")
        ),
    )


def _unavailable_metrics(
    config: MesoscopicCloudConfig,
    overlaps: dict[str, int],
    parent_disjoint: bool,
    child_isolated: bool,
    channels_disjoint: bool,
    full_outputs_disjoint: bool,
) -> MesoscopicCloudMetrics:
    gamma = config.gamma
    reasons = ["parent_boxes_overlap"]
    if not child_isolated:
        reasons.append("child_band_not_isolated")
    if not channels_disjoint:
        reasons.append("channel_bands_overlap")
    finite_reasons = tuple(dict.fromkeys(reasons + ["measurement_unavailable"]))
    asymptotic = list(finite_reasons)
    if gamma < 1.0:
        asymptotic.append("generic_duhamel_ratio_decays_for_gamma_below_one")
    width = config.width
    p = np.asarray(config.p_direction, dtype=float)
    q = np.asarray(config.q_direction, dtype=float)
    return MesoscopicCloudMetrics(
        base_scale=config.base_scale,
        gamma=gamma,
        width=width,
        half_width=width - 1,
        width_source="override" if config.width_override is not None else "power_floor",
        relative_width=width / config.base_scale,
        mode_count=4 * (2 * width - 1) ** 3,
        effective_mode_count=None,
        angle_sine=float(np.linalg.norm(np.cross(p, q)) / (np.linalg.norm(p) * np.linalg.norm(q))),
        box_overlap_counts=overlaps,
        parent_boxes_disjoint=parent_disjoint,
        child_band_isolated=child_isolated,
        channel_bands_disjoint=channels_disjoint,
        full_output_bands_disjoint=full_outputs_disjoint,
        measurement_available=False,
        measurement_method=None,
        targets_per_band=0,
        critical_energy=config.energy_constant / config.base_scale,
        parent_norm_squared=None,
        normalization_relative_error=None,
        divergence_relative=None,
        reality_relative=None,
        parent_divergence_relative=None,
        parent_reality_relative=None,
        nonlinear_divergence_relative=None,
        nonlinear_reality_relative=None,
        parent_min_wavenumber=None,
        parent_max_wavenumber=None,
        parent_rms_wavenumber=None,
        duhamel_upper_bound_method="unavailable-parent-overlap",
        support_radius_factor_kappa=None,
        universal_duhamel_upper_bound=None,
        effective_duhamel_upper_bound=None,
        duhamel_to_universal_upper_bound=None,
        duhamel_to_effective_upper_bound=None,
        a_critical=None,
        a_unit=None,
        chi_support=None,
        chi_effective=None,
        gain_g=None,
        heat_factor_h=None,
        duhamel_norm_ratio=None,
        duhamel_energy_ratio=None,
        duhamel_identity_rhs=None,
        duhamel_identity_relative_error=None,
        channel_norms=None,
        channel_ratios=None,
        channel_ratios_are_orthogonal=False,
        forcing_metric_method="unavailable-parent-overlap",
        child_forcing_fraction=None,
        low_side_forcing_ratio=None,
        off_chain_main_ratio=None,
        difference_sideband_main_ratio=None,
        self_interaction_cross_interaction_ratio=None,
        energy_cancellation_available=False,
        energy_cancellation_method="unavailable-parent-overlap",
        energy_cancellation_residual=None,
        energy_cancellation_pairing=None,
        full_nonlinear_norm=None,
        predicted_a_unit_exponent=1.0 + 1.5 * gamma,
        predicted_a_critical_exponent=1.5 * gamma,
        predicted_g_exponent=1.5 * gamma - 1.0,
        predicted_duhamel_exponent=3.0 * gamma - 3.0,
        predicted_g_grows=1.5 * gamma - 1.0 > 0.0,
        predicted_duhamel_decays=gamma < 1.0,
        finite_screen_pass=False,
        asymptotic_screen_pass=False,
        finite_screen_reasons=finite_reasons,
        asymptotic_screen_reasons=tuple(asymptotic),
    )


def measure_mesoscopic_cloud(config: MesoscopicCloudConfig) -> MesoscopicCloudMetrics:
    """Measure one cloud and apply finite and asymptotic rejection screens."""

    (
        overlaps,
        parent_disjoint,
        child_isolated,
        channels_disjoint,
        full_outputs_disjoint,
    ) = _overlap_diagnostics(config)
    if not parent_disjoint:
        return _unavailable_metrics(
            config,
            overlaps,
            parent_disjoint,
            child_isolated,
            channels_disjoint,
            full_outputs_disjoint,
        )

    width = config.width
    d = width - 1
    normalization = _normalization_scale(config)
    (
        parent_norm_squared,
        divergence_squared,
        coefficient_l1,
        minimum_wave,
        maximum_wave,
        gradient_squared,
    ) = _structure_metrics(config, normalization)
    target_norm_squared = 2.0 * config.energy_constant / config.base_scale
    normalization_error = abs(parent_norm_squared - target_norm_squared) / target_norm_squared
    divergence_relative = math.sqrt(divergence_squared / max(gradient_squared, np.finfo(float).tiny))
    effective_modes = coefficient_l1**2 / parent_norm_squared
    rms_wave = math.sqrt(gradient_squared / parent_norm_squared)

    core_quadrature, core_method = _target_quadrature(
        d,
        exact_target_limit=config.exact_target_limit,
        strata_per_axis=config.strata_per_axis,
    )
    full_quadrature, full_method = _target_quadrature(
        2 * d,
        exact_target_limit=config.exact_target_limit,
        strata_per_axis=config.strata_per_axis,
    )
    method = (
        "exact-target"
        if core_method == full_method == "exact-target"
        else "deterministic-stratified-target"
    )
    p = _integer_wave(config.p_direction, name="p_direction")
    q = _integer_wave(config.q_direction, name="q_direction")
    child = _add(p, q)
    difference = _subtract(p, q)
    child_core_positive, child_heat_positive = _measure_band(
        config, child, normalization, core_quadrature
    )
    child_full_positive, _ = _measure_band(
        config, child, normalization, full_quadrature
    )
    difference_core_positive, _ = _measure_band(
        config, difference, normalization, core_quadrature
    )
    difference_full_positive, _ = _measure_band(
        config, difference, normalization, full_quadrature
    )
    harmonic_p_core_positive, _ = _measure_band(
        config, _scale(2, p), normalization, core_quadrature
    )
    harmonic_p_full_positive, _ = _measure_band(
        config, _scale(2, p), normalization, full_quadrature
    )
    harmonic_q_core_positive, _ = _measure_band(
        config, _scale(2, q), normalization, core_quadrature
    )
    harmonic_q_full_positive, _ = _measure_band(
        config, _scale(2, q), normalization, full_quadrature
    )
    low_core, _ = _measure_band(
        config, (0, 0, 0), normalization, core_quadrature
    )
    low_full, _ = _measure_band(
        config, (0, 0, 0), normalization, full_quadrature
    )

    def spill(full: float, core: float) -> float:
        # The two deterministic quadratures are independent approximations in
        # sampled mode.  Clamp a possible tiny negative sampling difference;
        # exact-target mode has literal set inclusion and needs no correction.
        return max(full - core, 0.0)

    channel_squared = {
        "child": 2.0 * child_core_positive,
        "child_spill": 2.0 * spill(child_full_positive, child_core_positive),
        "difference": 2.0 * difference_core_positive,
        "difference_spill": 2.0
        * spill(difference_full_positive, difference_core_positive),
        "harmonic": 2.0 * harmonic_p_core_positive
        + 2.0 * harmonic_q_core_positive,
        "same_carrier_spill": 2.0
        * spill(harmonic_p_full_positive, harmonic_p_core_positive)
        + 2.0
        * spill(harmonic_q_full_positive, harmonic_q_core_positive),
        "low": low_core,
        "low_spill": spill(low_full, low_core),
    }
    total_channel_squared = sum(channel_squared.values())
    if not total_channel_squared > 0.0 or not channel_squared["child"] > 0.0:
        raise ValueError("the measured nonlinear response is zero")
    channel_norms = {
        name: math.sqrt(max(value, 0.0)) for name, value in channel_squared.items()
    }
    channel_ratios = {
        name: value / total_channel_squared for name, value in channel_squared.items()
    }
    channel_ratios["off_child"] = 1.0 - channel_ratios["child"]

    exact_full = _exact_full_diagnostic(config)
    if exact_full.available:
        if (
            exact_full.full_norm_squared is None
            or exact_full.child_core_squared is None
            or exact_full.child_full_squared is None
            or exact_full.difference_full_squared is None
            or exact_full.harmonic_full_squared is None
            or exact_full.low_full_squared is None
        ):
            raise AssertionError("available exact full diagnostic is incomplete")
        metric_child_squared = exact_full.child_core_squared
        metric_full_squared = exact_full.full_norm_squared
        metric_low_squared = exact_full.low_full_squared
        metric_difference_squared = exact_full.difference_full_squared
        child_forcing_fraction = metric_child_squared / metric_full_squared
        off_chain_main_ratio = math.sqrt(
            max(metric_full_squared - metric_child_squared, 0.0)
            / metric_child_squared
        )
        low_side_forcing_ratio = math.sqrt(
            metric_low_squared / metric_child_squared
        )
        difference_sideband_main_ratio = math.sqrt(
            metric_difference_squared / metric_child_squared
        )
        if full_outputs_disjoint:
            self_cross_ratio = math.sqrt(
                (
                    exact_full.harmonic_full_squared
                    + exact_full.low_full_squared
                )
                / (
                    exact_full.child_full_squared
                    + exact_full.difference_full_squared
                )
            )
        else:
            self_cross_ratio = None
        forcing_metric_method = "exact-sparse-full"
    else:
        metric_child_squared = channel_squared["child"]
        sampled_child_full_squared = (
            channel_squared["child"] + channel_squared["child_spill"]
        )
        sampled_difference_full_squared = (
            channel_squared["difference"]
            + channel_squared["difference_spill"]
        )
        sampled_harmonic_full_squared = (
            channel_squared["harmonic"]
            + channel_squared["same_carrier_spill"]
        )
        sampled_low_full_squared = (
            channel_squared["low"] + channel_squared["low_spill"]
        )
        low_side_forcing_ratio = math.sqrt(
            sampled_low_full_squared / metric_child_squared
        )
        difference_sideband_main_ratio = math.sqrt(
            sampled_difference_full_squared / metric_child_squared
        )
        if full_outputs_disjoint:
            disjoint_full_squared = total_channel_squared
            child_forcing_fraction = (
                metric_child_squared / disjoint_full_squared
            )
            off_chain_main_ratio = math.sqrt(
                max(disjoint_full_squared - metric_child_squared, 0.0)
                / metric_child_squared
            )
            self_cross_ratio = math.sqrt(
                (sampled_harmonic_full_squared + sampled_low_full_squared)
                / (
                    sampled_child_full_squared
                    + sampled_difference_full_squared
                )
            )
            forcing_metric_method = f"disjoint-band-{method}"
        else:
            child_forcing_fraction = None
            off_chain_main_ratio = None
            self_cross_ratio = None
            forcing_metric_method = f"overlapping-band-projections-only-{method}"

    a_critical = math.sqrt(channel_squared["child"])
    a_unit = a_critical / parent_norm_squared
    mode_count = 4 * (2 * width - 1) ** 3
    chi_support = a_unit / (config.base_scale * math.sqrt(mode_count))
    chi_effective = a_unit / (config.base_scale * math.sqrt(effective_modes))
    gain = a_critical / (config.base_scale**2 * parent_norm_squared)
    heat_weighted_squared = 2.0 * child_heat_positive
    heat_factor = heat_weighted_squared / channel_squared["child"]
    duhamel = (
        heat_weighted_squared
        / config.base_scale**4
        / parent_norm_squared
    )
    identity_rhs = heat_factor * gain * gain * parent_norm_squared
    identity_error = abs(duhamel - identity_rhs) / max(
        abs(duhamel), np.finfo(float).tiny
    )
    support_radius_factor = maximum_wave / config.base_scale
    universal_prefactor = (
        2.0
        * support_radius_factor**2
        * config.tau**2
        * config.energy_constant
        / config.base_scale**3
    )
    universal_duhamel_bound = universal_prefactor * mode_count
    effective_duhamel_bound = universal_prefactor * effective_modes
    roundoff_factor = 65536.0 * np.finfo(np.float64).eps
    if effective_duhamel_bound > universal_duhamel_bound * (
        1.0 + roundoff_factor
    ):
        raise AssertionError(
            "effective support count exceeds the universal support count"
        )
    for name, bound in (
        ("support-count", universal_duhamel_bound),
        ("effective-count", effective_duhamel_bound),
    ):
        tolerance = roundoff_factor * max(
            bound, np.finfo(np.float64).tiny
        )
        if duhamel > bound + tolerance:
            raise AssertionError(
                f"measured Duhamel ratio exceeds the {name} upper bound: "
                f"{duhamel} > {bound}"
            )
    duhamel_to_universal = duhamel / universal_duhamel_bound
    duhamel_to_effective = duhamel / effective_duhamel_bound

    finite_reasons: list[str] = []
    if not child_isolated:
        finite_reasons.append("child_band_not_isolated")
    if not channels_disjoint:
        finite_reasons.append("channel_bands_overlap")
    if duhamel < config.required_duhamel_ratio:
        finite_reasons.append("duhamel_ratio_below_required_target")
    finite_pass = not finite_reasons
    asymptotic_reasons = list(finite_reasons)
    if config.gamma < 1.0:
        asymptotic_reasons.append(
            "generic_duhamel_ratio_decays_for_gamma_below_one"
        )
    asymptotic_pass = not asymptotic_reasons
    p_array = np.asarray(p, dtype=float)
    q_array = np.asarray(q, dtype=float)
    gamma = config.gamma
    return MesoscopicCloudMetrics(
        base_scale=config.base_scale,
        gamma=gamma,
        width=width,
        half_width=d,
        width_source="override" if config.width_override is not None else "power_floor",
        relative_width=width / config.base_scale,
        mode_count=mode_count,
        effective_mode_count=effective_modes,
        angle_sine=float(
            np.linalg.norm(np.cross(p_array, q_array))
            / (np.linalg.norm(p_array) * np.linalg.norm(q_array))
        ),
        box_overlap_counts=overlaps,
        parent_boxes_disjoint=True,
        child_band_isolated=child_isolated,
        channel_bands_disjoint=channels_disjoint,
        full_output_bands_disjoint=full_outputs_disjoint,
        measurement_available=True,
        measurement_method=method,
        targets_per_band=len(full_quadrature),
        critical_energy=config.energy_constant / config.base_scale,
        parent_norm_squared=parent_norm_squared,
        normalization_relative_error=normalization_error,
        divergence_relative=divergence_relative,
        reality_relative=0.0,
        parent_divergence_relative=divergence_relative,
        parent_reality_relative=0.0,
        nonlinear_divergence_relative=None,
        nonlinear_reality_relative=None,
        parent_min_wavenumber=minimum_wave,
        parent_max_wavenumber=maximum_wave,
        parent_rms_wavenumber=rms_wave,
        duhamel_upper_bound_method=(
            "phase-independent-bernstein-heat-support-and-effective-count"
        ),
        support_radius_factor_kappa=support_radius_factor,
        universal_duhamel_upper_bound=universal_duhamel_bound,
        effective_duhamel_upper_bound=effective_duhamel_bound,
        duhamel_to_universal_upper_bound=duhamel_to_universal,
        duhamel_to_effective_upper_bound=duhamel_to_effective,
        a_critical=a_critical,
        a_unit=a_unit,
        chi_support=chi_support,
        chi_effective=chi_effective,
        gain_g=gain,
        heat_factor_h=heat_factor,
        duhamel_norm_ratio=math.sqrt(duhamel),
        duhamel_energy_ratio=duhamel,
        duhamel_identity_rhs=identity_rhs,
        duhamel_identity_relative_error=identity_error,
        channel_norms=channel_norms,
        channel_ratios=channel_ratios,
        channel_ratios_are_orthogonal=full_outputs_disjoint,
        forcing_metric_method=forcing_metric_method,
        child_forcing_fraction=child_forcing_fraction,
        low_side_forcing_ratio=low_side_forcing_ratio,
        off_chain_main_ratio=off_chain_main_ratio,
        difference_sideband_main_ratio=difference_sideband_main_ratio,
        self_interaction_cross_interaction_ratio=self_cross_ratio,
        energy_cancellation_available=exact_full.available,
        energy_cancellation_method=exact_full.method,
        energy_cancellation_residual=exact_full.cancellation_residual,
        energy_cancellation_pairing=exact_full.pairing,
        full_nonlinear_norm=(
            math.sqrt(exact_full.full_norm_squared)
            if exact_full.full_norm_squared is not None
            else None
        ),
        predicted_a_unit_exponent=1.0 + 1.5 * gamma,
        predicted_a_critical_exponent=1.5 * gamma,
        predicted_g_exponent=1.5 * gamma - 1.0,
        predicted_duhamel_exponent=3.0 * gamma - 3.0,
        predicted_g_grows=1.5 * gamma - 1.0 > 0.0,
        predicted_duhamel_decays=gamma < 1.0,
        finite_screen_pass=finite_pass,
        asymptotic_screen_pass=asymptotic_pass,
        finite_screen_reasons=tuple(finite_reasons),
        asymptotic_screen_reasons=tuple(asymptotic_reasons),
    )
