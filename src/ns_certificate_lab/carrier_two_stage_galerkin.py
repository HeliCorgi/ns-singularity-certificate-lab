r"""Small Fourier--Galerkin failure diagnostic for the partial carrier gadget.

The exact finite search in :mod:`ns_certificate_lab.exact_carrier_search`
finds two first-generation relays, but rejects their union because diagonal
parent interactions populate the target child shell.  This module starts
only the four searched parent modes and integrates the true unforced,
viscous, Leray-projected Fourier--Galerkin Navier--Stokes ODE.  It records
separately

* the two intended first children;
* the two exact diagonal cross-talk children;
* the intended sum and difference grandchildren; and
* all remaining Galerkin modes.

The Fourier grid is padding, not the Galerkin truncation: every product of
two modes in the retained cube lies strictly below FFT Nyquist.  Thus the
quadratic convolution is alias-free.  The time integration is binary64 RK4,
however, and supplies neither a continuum error enclosure nor a relay proof.
In particular, a populated grandchild mode cannot be attributed solely to
the intended child--child interaction because cross-talk paths are present.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np
import numpy.typing as npt

from .exact_carrier_search import ModeSpec, search_exact_carrier_gadget
from .leray_response_relay import (
    gradient_l2_squared,
    leray_advection,
    mean_energy,
)


ComplexArray = npt.NDArray[np.complex128]
BoolArray = npt.NDArray[np.bool_]
Wavevector = tuple[int, int, int]

__all__ = [
    "CarrierTwoStageResult",
    "build_partial_carrier_parent",
    "run_carrier_two_stage_galerkin",
]


_STATUS = (
    "BINARY64 GALERKIN FAILURE DIAGNOSTIC FOR PARTIAL REJECTED GADGET / "
    "NOT A RELAY OR CONTINUUM PROOF"
)
_GROUP_ORDER = (
    "parents",
    "first_child_one",
    "first_child_two",
    "cross_talk_one",
    "cross_talk_two",
    "grandchild_sum",
    "grandchild_difference",
)


def _frequencies(grid_size: int) -> npt.NDArray[np.float64]:
    return np.asarray(
        np.fft.fftfreq(grid_size, d=1.0 / grid_size), dtype=np.float64
    )


def _mesh(grid_size: int):
    frequency = _frequencies(grid_size)
    return tuple(np.meshgrid(frequency, frequency, frequency, indexing="ij"))


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
    wave = np.stack(_mesh(field.shape[1]), axis=0)
    divergence = 1.0j * np.sum(wave * field, axis=0)
    denominator = max(
        math.sqrt(gradient_l2_squared(field)), np.finfo(float).tiny
    )
    return float(np.linalg.norm(divergence) / denominator)


def _canonical_with_negatives(
    waves: tuple[Wavevector, ...], scale: int
) -> tuple[Wavevector, ...]:
    scaled: set[Wavevector] = set()
    for wave in waves:
        positive = tuple(scale * component for component in wave)
        negative = tuple(-component for component in positive)
        scaled.add(positive)  # type: ignore[arg-type]
        scaled.add(negative)  # type: ignore[arg-type]
    return tuple(sorted(scaled))


def _single_support(field) -> tuple[Wavevector, ...]:
    support = tuple(sorted(field.coefficient_table()))
    if len(support) != 1:
        raise AssertionError("an exact carrier channel is not a single mode")
    return support


def _carrier_groups(scale: int) -> dict[str, tuple[Wavevector, ...]]:
    certificate = search_exact_carrier_gadget()
    gadget = certificate.partial_gadget
    if gadget is None or gadget.classification != "partial_rejected_cross_talk":
        raise AssertionError("the registered partial rejected gadget is unavailable")

    relay_one = gadget.relay_one
    relay_two = gadget.relay_two
    parent_waves = (
        relay_one.left_spec.wave,
        relay_one.right_spec.wave,
        relay_two.left_spec.wave,
        relay_two.right_spec.wave,
    )
    base_groups: dict[str, tuple[Wavevector, ...]] = {
        "parents": parent_waves,
        "first_child_one": (relay_one.interaction.sum_wave,),
        "first_child_two": (relay_two.interaction.sum_wave,),
        "cross_talk_one": _single_support(
            gadget.diagonal_cross_talk.a1_b2_target_shell_output
        ),
        "cross_talk_two": _single_support(
            gadget.diagonal_cross_talk.a2_b1_target_shell_output
        ),
        "grandchild_sum": (gadget.next_interaction.sum_wave,),
        "grandchild_difference": (gadget.next_interaction.difference_wave,),
    }
    groups = {
        name: _canonical_with_negatives(waves, scale)
        for name, waves in base_groups.items()
    }
    occupied: set[Wavevector] = set()
    for name in _GROUP_ORDER:
        overlap = occupied.intersection(groups[name])
        if overlap:
            raise AssertionError(f"carrier diagnostic groups overlap: {overlap}")
        occupied.update(groups[name])
    return groups


def _mode_mask(
    grid_size: int, waves: tuple[Wavevector, ...]
) -> BoolArray:
    mask = np.zeros((grid_size, grid_size, grid_size), dtype=np.bool_)
    for wave in waves:
        mask[tuple(component % grid_size for component in wave)] = True
    return mask


def _galerkin_mask(grid_size: int, cutoff: int) -> BoolArray:
    kx, ky, kz = _mesh(grid_size)
    return np.asarray(
        (np.abs(kx) <= cutoff)
        & (np.abs(ky) <= cutoff)
        & (np.abs(kz) <= cutoff),
        dtype=np.bool_,
    )


def _insert_trig_mode(
    field: ComplexArray, spec: ModeSpec, scale: int
) -> None:
    wave = tuple(scale * component for component in spec.wave)
    index = tuple(component % field.shape[1] for component in wave)
    negative_index = tuple((-component) % field.shape[1] for component in wave)
    polarization = np.asarray(spec.polarization, dtype=np.float64)
    if spec.phase == "cos":
        positive = 0.5 * polarization.astype(np.complex128)
    elif spec.phase == "sin":
        positive = -0.5j * polarization
    else:  # pragma: no cover - protected by ModeSpec
        raise AssertionError("unknown exact carrier phase")
    field[(slice(None), *index)] += positive
    field[(slice(None), *negative_index)] += np.conjugate(positive)


def build_partial_carrier_parent(
    grid_size: int,
    *,
    scale: int = 4,
    energy_constant: float = 1.0,
) -> ComplexArray:
    """Build the four exact searched parents with total energy ``c_E/N``."""

    if isinstance(grid_size, bool) or not isinstance(grid_size, int) or grid_size < 8:
        raise ValueError("grid_size must be an integer at least eight")
    if isinstance(scale, bool) or not isinstance(scale, int) or scale < 1:
        raise ValueError("scale must be a positive integer")
    if (
        isinstance(energy_constant, bool)
        or not math.isfinite(energy_constant)
        or energy_constant <= 0.0
    ):
        raise ValueError("energy_constant must be finite and positive")
    if scale >= grid_size // 2:
        raise ValueError("parent carrier reaches the FFT Nyquist mode")

    certificate = search_exact_carrier_gadget()
    gadget = certificate.partial_gadget
    if gadget is None:
        raise AssertionError("the registered partial gadget is unavailable")
    specs = (
        gadget.relay_one.left_spec,
        gadget.relay_one.right_spec,
        gadget.relay_two.left_spec,
        gadget.relay_two.right_spec,
    )
    field = np.zeros((3, grid_size, grid_size, grid_size), dtype=np.complex128)
    for spec in specs:
        _insert_trig_mode(field, spec, scale)

    target_energy = energy_constant / float(scale)
    field *= math.sqrt(target_energy / mean_energy(field))
    if _reality_defect(field) > 64.0 * np.finfo(float).eps:
        raise AssertionError("constructed carrier parent is not Hermitian")
    if _divergence_defect(field) > 64.0 * np.finfo(float).eps:
        raise AssertionError("constructed carrier parent is not divergence free")
    return field


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


def _masked_energy(field: ComplexArray, mask: BoolArray) -> float:
    coefficients = field[:, mask]
    return 0.5 * float(np.vdot(coefficients, coefficients).real)


def _energy_identity_residual(
    state: ComplexArray,
    rhs: ComplexArray,
    *,
    viscosity: float,
) -> float:
    observed = float(np.vdot(state, rhs).real)
    expected = -viscosity * gradient_l2_squared(state)
    scale = max(abs(expected), np.finfo(float).tiny)
    return abs(observed - expected) / scale


def _history_record(
    state: ComplexArray,
    *,
    time: float,
    scale: int,
    parent_energy: float,
    group_masks: Mapping[str, BoolArray],
    remainder_mask: BoolArray,
) -> tuple[dict[str, float], float]:
    energies = {
        name: _masked_energy(state, group_masks[name])
        for name in _GROUP_ORDER
    }
    remainder = _masked_energy(state, remainder_mask)
    total = mean_energy(state)
    classified = sum(energies.values()) + remainder
    record = {
        "time": time,
        "scaled_time": time * scale * scale,
        "parent_energy_ratio": energies["parents"] / parent_energy,
        "first_child_one_energy_ratio": (
            energies["first_child_one"] / parent_energy
        ),
        "first_child_two_energy_ratio": (
            energies["first_child_two"] / parent_energy
        ),
        "first_child_energy_ratio": (
            energies["first_child_one"] + energies["first_child_two"]
        )
        / parent_energy,
        "cross_talk_one_energy_ratio": (
            energies["cross_talk_one"] / parent_energy
        ),
        "cross_talk_two_energy_ratio": (
            energies["cross_talk_two"] / parent_energy
        ),
        "cross_talk_energy_ratio": (
            energies["cross_talk_one"] + energies["cross_talk_two"]
        )
        / parent_energy,
        "grandchild_sum_energy_ratio": (
            energies["grandchild_sum"] / parent_energy
        ),
        "grandchild_difference_energy_ratio": (
            energies["grandchild_difference"] / parent_energy
        ),
        "grandchild_energy_ratio": (
            energies["grandchild_sum"]
            + energies["grandchild_difference"]
        )
        / parent_energy,
        "remainder_energy_ratio": remainder / parent_energy,
        "total_energy_ratio": total / parent_energy,
    }
    budget_residual = abs(classified - total) / max(
        total, np.finfo(float).tiny
    )
    return record, budget_residual


@dataclass(frozen=True)
class CarrierTwoStageResult:
    """One alias-free binary64 Galerkin trajectory and channel budget."""

    scale: int
    grid_size: int
    galerkin_cutoff: int
    dealias_margin: int
    steps: int
    tau: float
    time_multiples: float
    final_time: float
    viscosity: float
    energy_constant: float
    initial_parent_energy: float
    initial_named_nonparent_energy: float
    mode_groups: tuple[tuple[str, tuple[Wavevector, ...]], ...]
    final_parent_energy_ratio: float
    final_first_child_energy_ratio: float
    final_cross_talk_energy_ratio: float
    final_grandchild_energy_ratio: float
    final_remainder_energy_ratio: float
    final_total_energy_ratio: float
    maximum_energy_increase: float
    maximum_shell_budget_residual: float
    initial_energy_identity_residual: float
    final_energy_identity_residual: float
    final_reality_defect: float
    final_divergence_defect: float
    initial_grandchild_rhs_noise_ratio: float
    grandchild_roundoff_floor_ratio: float
    grandchild_resolution_margin: float
    grandchild_binary64_resolved: bool
    history: tuple[dict[str, float], ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "status": _STATUS,
            "classification": "partial_rejected_cross_talk",
            "interpretation": (
                "grandchild energy is pathway-contaminated and is not evidence "
                "for an iterable child-child relay"
            ),
            "scale": self.scale,
            "grid_size": self.grid_size,
            "galerkin_cutoff": self.galerkin_cutoff,
            "dealias_margin": self.dealias_margin,
            "dealias_verified": self.dealias_margin > 0,
            "steps": self.steps,
            "tau": self.tau,
            "time_multiples": self.time_multiples,
            "final_time": self.final_time,
            "viscosity": self.viscosity,
            "energy_constant": self.energy_constant,
            "initial_parent_energy": self.initial_parent_energy,
            "initial_named_nonparent_energy": (
                self.initial_named_nonparent_energy
            ),
            "mode_groups": {
                name: [list(wave) for wave in waves]
                for name, waves in self.mode_groups
            },
            "final_parent_energy_ratio": self.final_parent_energy_ratio,
            "final_first_child_energy_ratio": (
                self.final_first_child_energy_ratio
            ),
            "final_cross_talk_energy_ratio": (
                self.final_cross_talk_energy_ratio
            ),
            "final_grandchild_energy_ratio": (
                self.final_grandchild_energy_ratio
            ),
            "final_remainder_energy_ratio": (
                self.final_remainder_energy_ratio
            ),
            "final_total_energy_ratio": self.final_total_energy_ratio,
            "maximum_energy_increase": self.maximum_energy_increase,
            "maximum_shell_budget_residual": (
                self.maximum_shell_budget_residual
            ),
            "initial_energy_identity_residual": (
                self.initial_energy_identity_residual
            ),
            "final_energy_identity_residual": self.final_energy_identity_residual,
            "final_reality_defect": self.final_reality_defect,
            "final_divergence_defect": self.final_divergence_defect,
            "initial_grandchild_rhs_noise_ratio": (
                self.initial_grandchild_rhs_noise_ratio
            ),
            "grandchild_roundoff_floor_ratio": (
                self.grandchild_roundoff_floor_ratio
            ),
            "grandchild_resolution_margin": self.grandchild_resolution_margin,
            "grandchild_binary64_resolved": self.grandchild_binary64_resolved,
            "history": list(self.history),
        }


def run_carrier_two_stage_galerkin(
    *,
    scale: int = 4,
    grid_size: int = 64,
    galerkin_cutoff: int | None = None,
    viscosity: float = 1.0 / 40.0,
    energy_constant: float = 1.0,
    tau: float = 0.25,
    time_multiples: float = 2.0,
    steps: int = 16,
) -> CarrierTwoStageResult:
    r"""Integrate the partial four-parent gadget to ``m*tau*N^{-2}``.

    The default cutoff is the smallest simple cube containing both intended
    grandchildren, ``|k_j| <= 3N``.  The padding requirement
    ``2*cutoff < floor(grid_size/2)`` is checked before any FFT product.
    """

    if isinstance(scale, bool) or not isinstance(scale, int) or scale < 1:
        raise ValueError("scale must be a positive integer")
    if isinstance(grid_size, bool) or not isinstance(grid_size, int):
        raise ValueError("grid_size must be an integer")
    if (
        isinstance(viscosity, bool)
        or not math.isfinite(viscosity)
        or viscosity <= 0.0
    ):
        raise ValueError("viscosity must be finite and positive")
    if (
        isinstance(energy_constant, bool)
        or not math.isfinite(energy_constant)
        or energy_constant <= 0.0
    ):
        raise ValueError("energy_constant must be finite and positive")
    if isinstance(tau, bool) or not math.isfinite(tau) or tau <= 0.0:
        raise ValueError("tau must be finite and positive")
    if (
        isinstance(time_multiples, bool)
        or not math.isfinite(time_multiples)
        or time_multiples <= 0.0
    ):
        raise ValueError("time_multiples must be finite and positive")
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
        raise ValueError("steps must be a positive integer")

    groups = _carrier_groups(scale)
    required_cutoff = max(
        abs(component)
        for waves in groups.values()
        for wave in waves
        for component in wave
    )
    cutoff = 3 * scale if galerkin_cutoff is None else galerkin_cutoff
    if isinstance(cutoff, bool) or not isinstance(cutoff, int):
        raise ValueError("galerkin_cutoff must be an integer")
    if cutoff < required_cutoff:
        raise ValueError("Galerkin cutoff excludes a named carrier channel")
    if 2 * cutoff >= grid_size // 2:
        raise ValueError("grid is too small for exact padded Galerkin products")

    parent = build_partial_carrier_parent(
        grid_size, scale=scale, energy_constant=energy_constant
    )
    parent_energy = mean_energy(parent)
    group_masks = {
        name: _mode_mask(grid_size, groups[name]) for name in _GROUP_ORDER
    }
    galerkin_mask = _galerkin_mask(grid_size, cutoff)
    named_mask = np.zeros(galerkin_mask.shape, dtype=np.bool_)
    for mask in group_masks.values():
        if np.any(named_mask & mask):
            raise AssertionError("named carrier masks overlap")
        named_mask |= mask
    if np.any(named_mask & ~galerkin_mask):
        raise AssertionError("a named carrier mask lies outside Galerkin cutoff")
    remainder_mask = np.asarray(galerkin_mask & ~named_mask, dtype=np.bool_)

    initial_named_nonparent = sum(
        _masked_energy(parent, group_masks[name])
        for name in _GROUP_ORDER
        if name != "parents"
    )
    if initial_named_nonparent != 0.0:
        raise AssertionError("a named child channel is not initially empty")

    kx, ky, kz = _mesh(grid_size)
    wave_squared = kx * kx + ky * ky + kz * kz
    initial_rhs = _rhs(
        parent,
        viscosity=viscosity,
        wave_squared=wave_squared,
        galerkin_mask=galerkin_mask,
    )
    grandchild_mask = np.asarray(
        group_masks["grandchild_sum"]
        | group_masks["grandchild_difference"],
        dtype=np.bool_,
    )
    grandchild_rhs_noise = _masked_energy(initial_rhs, grandchild_mask)
    initial_rhs_energy = 0.5 * float(np.vdot(initial_rhs, initial_rhs).real)
    initial_grandchild_rhs_noise_ratio = grandchild_rhs_noise / max(
        initial_rhs_energy, np.finfo(float).tiny
    )
    initial_identity_residual = _energy_identity_residual(
        parent, initial_rhs, viscosity=viscosity
    )

    final_time = time_multiples * tau / float(scale * scale)
    dt = final_time / steps
    state = parent.copy()
    previous_energy = parent_energy
    maximum_energy_increase = 0.0
    maximum_budget_residual = 0.0
    first_record, first_budget = _history_record(
        state,
        time=0.0,
        scale=scale,
        parent_energy=parent_energy,
        group_masks=group_masks,
        remainder_mask=remainder_mask,
    )
    records = [first_record]
    maximum_budget_residual = max(maximum_budget_residual, first_budget)

    for index in range(steps):
        k1 = _rhs(
            state,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        k2 = _rhs(
            state + 0.5 * dt * k1,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        k3 = _rhs(
            state + 0.5 * dt * k2,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        k4 = _rhs(
            state + dt * k3,
            viscosity=viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        state = np.asarray(
            (
                state
                + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            )
            * galerkin_mask[None, ...],
            dtype=np.complex128,
        )
        total_energy = mean_energy(state)
        maximum_energy_increase = max(
            maximum_energy_increase, total_energy - previous_energy
        )
        previous_energy = total_energy
        record, budget_residual = _history_record(
            state,
            time=(index + 1) * dt,
            scale=scale,
            parent_energy=parent_energy,
            group_masks=group_masks,
            remainder_mask=remainder_mask,
        )
        records.append(record)
        maximum_budget_residual = max(
            maximum_budget_residual, budget_residual
        )

    final_rhs = _rhs(
        state,
        viscosity=viscosity,
        wave_squared=wave_squared,
        galerkin_mask=galerkin_mask,
    )
    final_identity_residual = _energy_identity_residual(
        state, final_rhs, viscosity=viscosity
    )
    final = records[-1]

    # A deliberately conservative binary64 floor: propagate both the actual
    # forbidden initial-channel FFT residue and a 4096-epsilon fraction of the
    # full initial right-hand-side norm for the whole time interval.
    epsilon_factor = 4096.0 * np.finfo(np.float64).eps
    rhs_roundoff_energy = max(
        grandchild_rhs_noise,
        epsilon_factor * epsilon_factor * initial_rhs_energy,
    )
    roundoff_floor_ratio = (
        final_time * final_time * rhs_roundoff_energy / parent_energy
    )
    grandchild_ratio = final["grandchild_energy_ratio"]
    resolution_margin = grandchild_ratio / max(
        roundoff_floor_ratio, np.finfo(float).tiny
    )
    resolved = bool(resolution_margin >= 100.0)

    return CarrierTwoStageResult(
        scale=scale,
        grid_size=grid_size,
        galerkin_cutoff=cutoff,
        dealias_margin=grid_size // 2 - 2 * cutoff,
        steps=steps,
        tau=tau,
        time_multiples=time_multiples,
        final_time=final_time,
        viscosity=viscosity,
        energy_constant=energy_constant,
        initial_parent_energy=parent_energy,
        initial_named_nonparent_energy=initial_named_nonparent,
        mode_groups=tuple((name, groups[name]) for name in _GROUP_ORDER),
        final_parent_energy_ratio=final["parent_energy_ratio"],
        final_first_child_energy_ratio=final["first_child_energy_ratio"],
        final_cross_talk_energy_ratio=final["cross_talk_energy_ratio"],
        final_grandchild_energy_ratio=grandchild_ratio,
        final_remainder_energy_ratio=final["remainder_energy_ratio"],
        final_total_energy_ratio=final["total_energy_ratio"],
        maximum_energy_increase=maximum_energy_increase,
        maximum_shell_budget_residual=maximum_budget_residual,
        initial_energy_identity_residual=initial_identity_residual,
        final_energy_identity_residual=final_identity_residual,
        final_reality_defect=_reality_defect(state),
        final_divergence_defect=_divergence_defect(state),
        initial_grandchild_rhs_noise_ratio=(
            initial_grandchild_rhs_noise_ratio
        ),
        grandchild_roundoff_floor_ratio=roundoff_floor_ratio,
        grandchild_resolution_margin=resolution_margin,
        grandchild_binary64_resolved=resolved,
        history=tuple(records),
    )
