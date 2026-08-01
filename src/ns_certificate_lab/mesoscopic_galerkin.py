r"""Small full-Galerkin check for an initially empty Leray child band.

This module is intentionally a *small cross-check*, not the large-
``(N, gamma)`` scaling engine.  It constructs two angle-separated parent
boxes around ``N p`` and ``N q`` for

``p=(1,1,0)``, ``q=(1,0,1)``

with affine Fourier phases and modewise Leray-projected polarizations.  The
designated child *core* around ``N(p+q)`` has half-width ``width-1`` and is
empty at time zero.  The full parent-parent geometric sumset has twice that
half-width; its outer ring is deliberately excluded from the child observable
and is reported separately.  We compare

1. the exact frozen-parent heat-Duhamel response in that child box; and
2. an RK4 integration of the unforced Fourier--Galerkin Navier--Stokes ODE.

The Galerkin state is projected to a fixed cubic cutoff after every right-hand
side evaluation.  The physical FFT grid is padded so the product of any two
retained modes lies strictly below Nyquist; hence this calculation does not
confuse Galerkin truncation with FFT aliasing.  It is binary64 and has no
continuum error enclosure.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt

from .leray_response_relay import (
    gradient_l2_squared,
    leray_advection,
    mean_energy,
)


ComplexArray = npt.NDArray[np.complex128]
BoolArray = npt.NDArray[np.bool_]

__all__ = [
    "GalerkinRelayResult",
    "build_angle_box_parent",
    "run_small_mesoscopic_galerkin",
]


_P = np.array((1, 1, 0), dtype=np.int64)
_Q = np.array((1, 0, 1), dtype=np.int64)
_C = _P + _Q


def _frequencies(grid_size: int) -> npt.NDArray[np.float64]:
    return np.asarray(
        np.fft.fftfreq(grid_size, d=1.0 / grid_size), dtype=np.float64
    )


def _mesh(grid_size: int):
    frequency = _frequencies(grid_size)
    return tuple(np.meshgrid(frequency, frequency, frequency, indexing="ij"))


def _projected_unit(
    waves: npt.NDArray[np.float64], polarization: npt.ArrayLike
) -> npt.NDArray[np.float64]:
    vector = np.asarray(polarization, dtype=np.float64)
    norm_squared = np.sum(waves * waves, axis=-1)
    projected = vector - waves * (
        np.sum(waves * vector, axis=-1) / norm_squared
    )[..., None]
    norms = np.linalg.norm(projected, axis=-1)
    if np.any(norms <= np.finfo(np.float64).eps):
        raise ValueError("a carrier polarization became parallel to its mode")
    return np.asarray(projected / norms[..., None], dtype=np.float64)


def _box_offsets(width: int) -> npt.NDArray[np.int64]:
    if isinstance(width, bool) or not isinstance(width, int) or width < 1:
        raise ValueError("width must be a positive integer")
    half_width = width - 1
    coordinates = np.arange(-half_width, half_width + 1, dtype=np.int64)
    mesh = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    return np.stack(mesh, axis=-1).reshape(-1, 3)


def _indices(waves: npt.NDArray[np.int64], grid_size: int):
    reduced = np.mod(waves, grid_size)
    return tuple(reduced[:, axis] for axis in range(3))


def _negative_reflection(field: ComplexArray) -> ComplexArray:
    n = field.shape[1]
    indices = (-np.arange(n, dtype=np.int64)) % n
    reflected = np.take(field, indices, axis=1)
    reflected = np.take(reflected, indices, axis=2)
    return np.take(reflected, indices, axis=3)


def _reality_defect(field: ComplexArray) -> float:
    reflected = _negative_reflection(field)
    denominator = max(float(np.linalg.norm(field)), np.finfo(float).tiny)
    return float(np.linalg.norm(field - np.conjugate(reflected)) / denominator)


def _divergence_defect(field: ComplexArray) -> float:
    n = field.shape[1]
    wave = np.stack(_mesh(n), axis=0)
    divergence = 1.0j * np.sum(wave * field, axis=0)
    numerator = float(np.linalg.norm(divergence))
    denominator = max(
        math.sqrt(gradient_l2_squared(field)), np.finfo(float).tiny
    )
    return numerator / denominator


def _child_mask(grid_size: int, *, scale: int, width: int) -> BoolArray:
    offsets = _box_offsets(width)
    positive = scale * _C[None, :] + offsets
    negative = -positive
    mask = np.zeros((grid_size, grid_size, grid_size), dtype=np.bool_)
    mask[_indices(positive, grid_size)] = True
    mask[_indices(negative, grid_size)] = True
    return mask


def _full_child_sumset_mask(
    grid_size: int, *, scale: int, width: int
) -> BoolArray:
    half_width = width - 1
    offsets = _box_offsets(2 * half_width + 1)
    positive = scale * _C[None, :] + offsets
    negative = -positive
    mask = np.zeros((grid_size, grid_size, grid_size), dtype=np.bool_)
    mask[_indices(positive, grid_size)] = True
    mask[_indices(negative, grid_size)] = True
    return mask


def build_angle_box_parent(
    grid_size: int,
    *,
    scale: int,
    width: int,
    energy_constant: float = 1.0,
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> tuple[ComplexArray, BoolArray]:
    r"""Return a critical-energy real parent and its initially empty child mask.

    The offset set is ``{-(width-1), ..., width-1}^3``.  The conservative
    condition ``scale > 3*(width-1)`` keeps parent, low, and child boxes
    separated.  Positive representatives receive affine phases; negative
    representatives are inserted as their exact complex conjugates.
    """

    if isinstance(grid_size, bool) or not isinstance(grid_size, int):
        raise ValueError("grid_size must be an integer")
    if isinstance(scale, bool) or not isinstance(scale, int) or scale < 2:
        raise ValueError("scale must be an integer at least two")
    if (
        isinstance(energy_constant, bool)
        or not math.isfinite(energy_constant)
        or energy_constant <= 0.0
    ):
        raise ValueError("energy_constant must be finite and positive")
    offsets = _box_offsets(width)
    half_width = width - 1
    if scale <= 3 * half_width:
        raise ValueError("parent, low, and child boxes are not separated")
    maximum_parent_coordinate = scale + half_width
    if 2 * maximum_parent_coordinate >= grid_size // 2:
        raise ValueError("parent self-product reaches the FFT Nyquist mode")

    field = np.zeros((3, grid_size, grid_size, grid_size), dtype=np.complex128)
    translation_vector = np.asarray(translation, dtype=np.float64)
    if translation_vector.shape != (3,) or not np.all(
        np.isfinite(translation_vector)
    ):
        raise ValueError("translation must be a finite three-vector")
    for center, polarization, phase_offset in (
        (scale * _P, (0.0, 0.0, 1.0), -0.5 * math.pi),
        (scale * _Q, (0.0, 1.0, 0.0), 0.0),
    ):
        waves = center[None, :] + offsets
        vectors = _projected_unit(waves.astype(np.float64), polarization)
        phases = np.exp(
            1.0j * (waves @ translation_vector + phase_offset)
        )
        coefficients = vectors * phases[:, None]
        field[(slice(None), *_indices(waves, grid_size))] += coefficients.T
        field[(slice(None), *_indices(-waves, grid_size))] += np.conjugate(
            coefficients.T
        )

    target_energy = energy_constant / float(scale)
    current_energy = mean_energy(field)
    field *= math.sqrt(target_energy / current_energy)
    child = _child_mask(grid_size, scale=scale, width=width)
    if np.any(child & np.any(np.abs(field) > 0.0, axis=0)):
        raise AssertionError("the child mask is not empty at time zero")
    return field, child


def _galerkin_mask(grid_size: int, cutoff: int) -> BoolArray:
    kx, ky, kz = _mesh(grid_size)
    return np.asarray(
        (np.abs(kx) <= cutoff)
        & (np.abs(ky) <= cutoff)
        & (np.abs(kz) <= cutoff),
        dtype=np.bool_,
    )


def _rhs(
    state: ComplexArray,
    *,
    viscosity: float,
    wave_squared: npt.NDArray[np.float64],
    galerkin_mask: BoolArray,
) -> ComplexArray:
    retained = np.asarray(
        state * galerkin_mask[None, ...], dtype=np.complex128
    )
    nonlinear = leray_advection(retained, retained)
    return np.asarray(
        (-nonlinear - viscosity * wave_squared[None, ...] * retained)
        * galerkin_mask[None, ...],
        dtype=np.complex128,
    )


@dataclass(frozen=True)
class GalerkinRelayResult:
    """Binary64 comparison between frozen and evolving-parent child filling."""

    scale: int
    width: int
    grid_size: int
    galerkin_cutoff: int
    dealias_margin: int
    steps: int
    tau: float
    viscosity: float
    parent_energy: float
    initial_child_energy: float
    child_core_half_width: int
    child_full_sumset_half_width: int
    child_core_mode_count: int
    child_full_sumset_mode_count: int
    child_core_forcing_fraction: float
    frozen_child_energy_ratio: float
    full_child_energy_ratio: float
    full_to_frozen_ratio: float
    final_total_energy_ratio: float
    maximum_energy_increase: float
    final_reality_defect: float
    final_divergence_defect: float
    history: tuple[dict[str, float], ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "status": "SMALL BINARY64 GALERKIN CROSS-CHECK / NOT CONTINUUM",
            "scale": self.scale,
            "width": self.width,
            "grid_size": self.grid_size,
            "galerkin_cutoff": self.galerkin_cutoff,
            "dealias_margin": self.dealias_margin,
            "steps": self.steps,
            "tau": self.tau,
            "viscosity": self.viscosity,
            "parent_energy": self.parent_energy,
            "initial_child_energy": self.initial_child_energy,
            "child_core_half_width": self.child_core_half_width,
            "child_full_sumset_half_width": self.child_full_sumset_half_width,
            "child_core_mode_count": self.child_core_mode_count,
            "child_full_sumset_mode_count": self.child_full_sumset_mode_count,
            "child_core_forcing_fraction": self.child_core_forcing_fraction,
            "frozen_child_energy_ratio": self.frozen_child_energy_ratio,
            "full_child_energy_ratio": self.full_child_energy_ratio,
            "full_to_frozen_ratio": self.full_to_frozen_ratio,
            "final_total_energy_ratio": self.final_total_energy_ratio,
            "maximum_energy_increase": self.maximum_energy_increase,
            "final_reality_defect": self.final_reality_defect,
            "final_divergence_defect": self.final_divergence_defect,
            "history": list(self.history),
        }


def run_small_mesoscopic_galerkin(
    *,
    scale: int = 4,
    width: int = 2,
    grid_size: int = 64,
    viscosity: float = 1.0 / 40.0,
    energy_constant: float = 1.0,
    tau: float = 0.25,
    steps: int = 16,
) -> GalerkinRelayResult:
    """Run the registered small full-Galerkin versus frozen-parent check."""

    if (
        isinstance(viscosity, bool)
        or not math.isfinite(viscosity)
        or viscosity <= 0.0
    ):
        raise ValueError("viscosity must be finite and positive")
    if isinstance(tau, bool) or not math.isfinite(tau) or tau <= 0.0:
        raise ValueError("tau must be finite and positive")
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
        raise ValueError("steps must be a positive integer")
    half_width = width - 1
    cutoff = 2 * scale + 2 * half_width
    # leray_advection rejects products whose coordinate support reaches
    # floor(grid_size/2).  This direct form is correct for odd as well as even
    # grids; the tempting 4*cutoff < grid_size shortcut fails for odd n.
    if 2 * cutoff >= grid_size // 2:
        raise ValueError("grid is too small for exact padded Galerkin products")

    parent, child_mask = build_angle_box_parent(
        grid_size,
        scale=scale,
        width=width,
        energy_constant=energy_constant,
    )
    parent_energy = mean_energy(parent)
    initial_child_energy = mean_energy(
        np.asarray(parent * child_mask[None, ...], dtype=np.complex128)
    )
    nonlinear = leray_advection(parent, parent)
    forcing = np.asarray(-nonlinear * child_mask[None, ...], dtype=np.complex128)
    full_child_mask = _full_child_sumset_mask(
        grid_size, scale=scale, width=width
    )
    if np.any(full_child_mask & np.any(np.abs(parent) > 0.0, axis=0)):
        raise AssertionError("the full child sumset is not empty at time zero")
    full_child_forcing = np.asarray(
        -nonlinear * full_child_mask[None, ...], dtype=np.complex128
    )
    core_forcing_squared = float(np.vdot(forcing, forcing).real)
    full_forcing_squared = float(np.vdot(full_child_forcing, full_child_forcing).real)
    if not full_forcing_squared > 0.0:
        raise ValueError("the full child sumset has no resolved forcing")
    child_core_forcing_fraction = core_forcing_squared / full_forcing_squared
    kx, ky, kz = _mesh(grid_size)
    wave_squared = kx * kx + ky * ky + kz * kz
    final_time = tau / float(scale * scale)
    heat_factor = np.full(wave_squared.shape, final_time, dtype=np.float64)
    nonzero_wave = wave_squared > 0.0
    argument = viscosity * wave_squared[nonzero_wave] * final_time
    heat_factor[nonzero_wave] = -np.expm1(-argument) / (
        viscosity * wave_squared[nonzero_wave]
    )
    frozen = forcing * heat_factor[None, ...]
    frozen_ratio = mean_energy(frozen) / parent_energy

    mask = _galerkin_mask(grid_size, cutoff)
    state = parent.copy()
    dt = final_time / steps
    initial_total_energy = mean_energy(state)
    previous_energy = initial_total_energy
    maximum_increase = 0.0
    records: list[dict[str, float]] = [
        {
            "time": 0.0,
            "child_energy_ratio": 0.0,
            "total_energy_ratio": 1.0,
        }
    ]
    for index in range(steps):
        k1 = _rhs(
            state,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=mask,
        )
        k2 = _rhs(
            state + 0.5 * dt * k1,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=mask,
        )
        k3 = _rhs(
            state + 0.5 * dt * k2,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=mask,
        )
        k4 = _rhs(
            state + dt * k3,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=mask,
        )
        state = np.asarray(
            (state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4))
            * mask[None, ...],
            dtype=np.complex128,
        )
        total_energy = mean_energy(state)
        maximum_increase = max(maximum_increase, total_energy - previous_energy)
        previous_energy = total_energy
        child_energy = mean_energy(
            np.asarray(state * child_mask[None, ...], dtype=np.complex128)
        )
        records.append(
            {
                "time": (index + 1) * dt,
                "child_energy_ratio": child_energy / parent_energy,
                "total_energy_ratio": total_energy / initial_total_energy,
            }
        )

    full_ratio = records[-1]["child_energy_ratio"]
    return GalerkinRelayResult(
        scale=scale,
        width=width,
        grid_size=grid_size,
        galerkin_cutoff=cutoff,
        dealias_margin=grid_size // 2 - 2 * cutoff,
        steps=steps,
        tau=tau,
        viscosity=viscosity,
        parent_energy=parent_energy,
        initial_child_energy=initial_child_energy,
        child_core_half_width=half_width,
        child_full_sumset_half_width=2 * half_width,
        child_core_mode_count=int(np.count_nonzero(child_mask)),
        child_full_sumset_mode_count=int(np.count_nonzero(full_child_mask)),
        child_core_forcing_fraction=child_core_forcing_fraction,
        frozen_child_energy_ratio=frozen_ratio,
        full_child_energy_ratio=full_ratio,
        full_to_frozen_ratio=full_ratio / frozen_ratio,
        final_total_energy_ratio=previous_energy / initial_total_energy,
        maximum_energy_increase=maximum_increase,
        final_reality_defect=_reality_defect(state),
        final_divergence_defect=_divergence_defect(state),
        history=tuple(records),
    )
