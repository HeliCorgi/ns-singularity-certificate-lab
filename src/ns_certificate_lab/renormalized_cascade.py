r"""Renormalized Galerkin cascade: forward integration of the doubling map.

Discovery-phase pilot for the fixed-relative-width Leray cloud lane.  The
candidate document poses shape closure as a projective fixed/periodic-point
problem for the continuum doubling response operator ``T_2`` (its eq (6.10)),
which so far has only been attacked as a one-step frozen-response optimization.
This module instead *integrates the renormalization flow forward*:

1. evolve the full Fourier--Galerkin Navier--Stokes ODE for one parabolic
   stage time ``tau / scale**2`` (RK4, exact padded products, no aliasing);
2. apply the exact lattice doubling pullback -- ``T_2`` with ``lambda = 2``
   restricted to the integer lattice is even-mode subsampling,

   ``w(k) = 4 v(2k)``,

   which preserves Hermitian reality and exact divergence-freeness
   (``k . w(k) = 2 (2k) . v(2k) = 0``) and whose energy bookkeeping matches
   the continuum relation ``E(T_2 U) = 2 * (child-band energy)`` through the
   factor ``16 * (1/8)`` of coefficient scaling times even-decimation;
3. optionally discard the sub-front sea (all modes with ``max_i |k_i|`` below
   a floor) to isolate the moving front;
4. record the pre-normalization gain ``A_j = ||w||^2 / ||u_j||^2`` and shape
   overlaps at lags 1..3, then projectively renormalize the energy back to
   the critical target ``c_E / scale`` and repeat.

If the projective doubling map has an *attracting* fixed point or periodic
orbit with positive gain, this iteration converges to it without any
optimization; if every orbit decays or disperses in the scanned parameter
box, the recorded gain and overlap sequences falsify attracting-orbit closure
there.  By the ``c_E``-collapse observation (all quadratic channels scale as
``c_E**1.5`` against viscous ``c_E``), a *stable positive* gain plateau at any
single ``c_E`` is the survival signal; the absolute size of the plateau is
tunable by ``c_E`` and is not itself the pass/fail criterion.

Status: binary64 discovery diagnostic.  No continuum enclosure, no proof, and
stopping at any finite stage says nothing about the PDE.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt

from .leray_response_relay import leray_advection, mean_energy
from .mesoscopic_galerkin import build_angle_box_parent

ComplexArray = npt.NDArray[np.complex128]
BoolArray = npt.NDArray[np.bool_]

__all__ = [
    "CascadeConfig",
    "CascadeResult",
    "CascadeStageRecord",
    "doubling_pullback",
    "linf_shell_energies",
    "run_renormalized_cascade",
]


def _grid_size(field: ComplexArray) -> int:
    if field.ndim != 4 or field.shape[0] != 3:
        raise ValueError("expected a (3, n, n, n) Fourier field")
    n = field.shape[1]
    if field.shape[1:] != (n, n, n):
        raise ValueError("expected a cubic Fourier grid")
    return n


def _integer_frequencies(grid_size: int) -> npt.NDArray[np.int64]:
    return np.rint(
        np.fft.fftfreq(grid_size, d=1.0 / grid_size)
    ).astype(np.int64)


def _frequency_mesh(grid_size: int):
    frequency = _integer_frequencies(grid_size).astype(np.float64)
    return tuple(np.meshgrid(frequency, frequency, frequency, indexing="ij"))


def doubling_pullback(field: ComplexArray) -> ComplexArray:
    """Return ``w`` with ``w(k) = 4 field(2k)`` on the integer lattice.

    Only target modes ``k`` with ``2k`` strictly below the Nyquist row are
    filled; everything else is zero.  This is the ``lambda = 2`` quadratic
    renormalisation map of the candidate document restricted to the lattice.
    """

    n = _grid_size(field)
    frequency = _integer_frequencies(n)
    limit = (n // 2 - 1) // 2
    valid = np.where(np.abs(frequency) <= limit)[0]
    source = (2 * frequency[valid]) % n
    out = np.zeros_like(field)
    out[np.ix_(range(3), valid, valid, valid)] = (
        4.0 * field[np.ix_(range(3), source, source, source)]
    )
    return out


def linf_shell_energies(
    field: ComplexArray, band_edges: tuple[int, ...]
) -> tuple[float, ...]:
    """Return energies in ``max_i |k_i|`` bands split at ``band_edges``.

    ``band_edges = (a, b)`` produces bands ``[0, a)``, ``[a, b)``,
    ``[b, infinity)``.
    """

    n = _grid_size(field)
    kx, ky, kz = _frequency_mesh(n)
    linf = np.maximum(np.abs(kx), np.maximum(np.abs(ky), np.abs(kz)))
    density = 0.5 * np.sum(np.abs(field) ** 2, axis=0)
    edges = (0, *band_edges, None)
    energies = []
    for low, high in zip(edges[:-1], edges[1:]):
        mask = linf >= low if high is None else (linf >= low) & (linf < high)
        energies.append(float(np.sum(density[mask])))
    return tuple(energies)


def _reality_defect(field: ComplexArray) -> float:
    n = field.shape[1]
    indices = (-np.arange(n, dtype=np.int64)) % n
    reflected = np.take(field, indices, axis=1)
    reflected = np.take(reflected, indices, axis=2)
    reflected = np.take(reflected, indices, axis=3)
    denominator = max(float(np.linalg.norm(field)), np.finfo(float).tiny)
    return float(
        np.linalg.norm(field - np.conjugate(reflected)) / denominator
    )


def _divergence_defect(field: ComplexArray) -> float:
    n = field.shape[1]
    wave = np.stack(_frequency_mesh(n), axis=0)
    divergence = 1.0j * np.sum(wave * field, axis=0)
    weight = np.sum(wave * wave, axis=0)
    denominator = math.sqrt(
        float(np.sum(weight[None, ...] * np.abs(field) ** 2).real)
    )
    return float(
        np.linalg.norm(divergence) / max(denominator, np.finfo(float).tiny)
    )


def _sup_velocity(field: ComplexArray) -> float:
    n = field.shape[1]
    physical = np.fft.ifftn(field, axes=(1, 2, 3)) * float(n**3)
    return float(np.max(np.sqrt(np.sum(np.abs(physical) ** 2, axis=0))))


def _shape_overlap(left: ComplexArray, right: ComplexArray) -> float:
    numerator = abs(complex(np.vdot(left, right)))
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator <= 0.0:
        return 0.0
    return numerator / denominator


@dataclass(frozen=True)
class CascadeConfig:
    """Configuration of one renormalized-cascade run."""

    scale: int = 4
    width: int = 2
    grid_size: int = 64
    viscosity: float = 1.0 / 40.0
    energy_constant: float = 1.0
    tau: float = 0.25
    base_steps: int = 16
    max_steps: int = 512
    stages: int = 10
    drop_below: int | None = 3
    cfl_safety: float = 8.0

    def validate(self) -> None:
        if isinstance(self.scale, bool) or not isinstance(self.scale, int):
            raise ValueError("scale must be an integer")
        if isinstance(self.width, bool) or not isinstance(self.width, int):
            raise ValueError("width must be an integer")
        for name in ("viscosity", "energy_constant", "tau", "cfl_safety"):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not math.isfinite(float(value))
                or float(value) <= 0.0
            ):
                raise ValueError(f"{name} must be finite and positive")
        for name in ("base_steps", "max_steps", "stages"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.max_steps < self.base_steps:
            raise ValueError("max_steps must be at least base_steps")
        if self.drop_below is not None and (
            isinstance(self.drop_below, bool)
            or not isinstance(self.drop_below, int)
            or self.drop_below < 1
        ):
            raise ValueError("drop_below must be None or a positive integer")

    @property
    def cutoff(self) -> int:
        return 2 * self.scale + 2 * (self.width - 1)

    @property
    def stage_time(self) -> float:
        return self.tau / float(self.scale * self.scale)


@dataclass(frozen=True)
class CascadeStageRecord:
    """Diagnostics for one evolve-pullback-renormalize stage."""

    stage: int
    steps: int
    evolved_energy_ratio: float
    child_band_energy_fraction: float
    pullback_gain: float
    dropped_sea_fraction: float
    overlap_lag1: float
    overlap_lag2: float
    overlap_lag3: float
    shell_energies: tuple[float, ...]
    sup_velocity: float
    reality_defect: float
    divergence_defect: float

    def as_dict(self) -> dict[str, object]:
        return {
            "stage": self.stage,
            "steps": self.steps,
            "evolved_energy_ratio": self.evolved_energy_ratio,
            "child_band_energy_fraction": self.child_band_energy_fraction,
            "pullback_gain": self.pullback_gain,
            "dropped_sea_fraction": self.dropped_sea_fraction,
            "overlap_lag1": self.overlap_lag1,
            "overlap_lag2": self.overlap_lag2,
            "overlap_lag3": self.overlap_lag3,
            "shell_energies": list(self.shell_energies),
            "sup_velocity": self.sup_velocity,
            "reality_defect": self.reality_defect,
            "divergence_defect": self.divergence_defect,
        }


@dataclass(frozen=True)
class CascadeResult:
    """Full history of one renormalized-cascade run."""

    config: CascadeConfig
    completed_stages: int
    terminated_reason: str
    stage_records: tuple[CascadeStageRecord, ...]

    @property
    def gains(self) -> tuple[float, ...]:
        return tuple(record.pullback_gain for record in self.stage_records)

    def as_dict(self) -> dict[str, object]:
        return {
            "status": "BINARY64 DISCOVERY DIAGNOSTIC / NOT A PROOF",
            "config": {
                "scale": self.config.scale,
                "width": self.config.width,
                "grid_size": self.config.grid_size,
                "viscosity": self.config.viscosity,
                "energy_constant": self.config.energy_constant,
                "tau": self.config.tau,
                "base_steps": self.config.base_steps,
                "max_steps": self.config.max_steps,
                "stages": self.config.stages,
                "drop_below": self.config.drop_below,
                "cfl_safety": self.config.cfl_safety,
                "cutoff": self.config.cutoff,
                "stage_time": self.config.stage_time,
            },
            "completed_stages": self.completed_stages,
            "terminated_reason": self.terminated_reason,
            "stage_records": [
                record.as_dict() for record in self.stage_records
            ],
        }


def _galerkin_mask(grid_size: int, cutoff: int) -> BoolArray:
    kx, ky, kz = _frequency_mesh(grid_size)
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


def _evolve(
    state: ComplexArray,
    *,
    config: CascadeConfig,
    wave_squared: npt.NDArray[np.float64],
    galerkin_mask: BoolArray,
) -> tuple[ComplexArray, int]:
    stage_time = config.stage_time
    advective_rate = config.cutoff * _sup_velocity(state)
    steps = config.base_steps
    if advective_rate > 0.0:
        needed = int(
            math.ceil(stage_time * advective_rate * config.cfl_safety)
        )
        steps = min(max(steps, needed), config.max_steps)
    dt = stage_time / steps
    for _ in range(steps):
        k1 = _rhs(
            state,
            viscosity=config.viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        k2 = _rhs(
            state + 0.5 * dt * k1,
            viscosity=config.viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        k3 = _rhs(
            state + 0.5 * dt * k2,
            viscosity=config.viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        k4 = _rhs(
            state + dt * k3,
            viscosity=config.viscosity,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        state = np.asarray(
            (state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4))
            * galerkin_mask[None, ...],
            dtype=np.complex128,
        )
        # The mean mode is exactly zero for the continuous dynamics; zeroing
        # the accumulated FFT roundoff there is hygiene, not physics.
        state[:, 0, 0, 0] = 0.0
    return state, steps


def run_renormalized_cascade(config: CascadeConfig) -> CascadeResult:
    """Iterate evolve -> doubling pullback -> renormalize and record gains."""

    config.validate()
    if 2 * config.cutoff >= config.grid_size // 2:
        raise ValueError("grid is too small for exact padded Galerkin products")

    parent, _child_mask = build_angle_box_parent(
        config.grid_size,
        scale=config.scale,
        width=config.width,
        energy_constant=config.energy_constant,
    )
    target_energy = config.energy_constant / float(config.scale)

    kx, ky, kz = _frequency_mesh(config.grid_size)
    wave_squared = kx * kx + ky * ky + kz * kz
    linf = np.maximum(np.abs(kx), np.maximum(np.abs(ky), np.abs(kz)))
    galerkin_mask = _galerkin_mask(config.grid_size, config.cutoff)
    # Child band of the *current* stage before pullback: modes whose halved
    # image lies at or above the front floor.  Used only as a diagnostic.
    child_band = np.asarray(
        (linf >= config.scale + config.width)
        & (linf <= 2 * (config.scale + config.width - 1)),
        dtype=np.bool_,
    )

    state = parent
    history: list[ComplexArray] = [state.copy()]
    records: list[CascadeStageRecord] = []
    terminated = "completed"
    # Projective renormalization re-amplifies the binary64 roundoff floor by
    # 1/sqrt(gain) per stage; once the estimated floor reaches a tenth of the
    # renormalized field the iteration is measuring noise and must stop.
    noise_floor_estimate = 1.0e-15
    for stage in range(config.stages):
        if noise_floor_estimate > 0.1:
            terminated = "float_noise_floor_reached"
            break
        initial_energy = mean_energy(state)
        evolved, steps = _evolve(
            state,
            config=config,
            wave_squared=wave_squared,
            galerkin_mask=galerkin_mask,
        )
        evolved_energy = mean_energy(evolved)
        child_energy = mean_energy(
            np.asarray(evolved * child_band[None, ...], dtype=np.complex128)
        )
        pulled = doubling_pullback(evolved)
        pulled_energy_full = mean_energy(pulled)
        if config.drop_below is not None:
            sea_mask = linf < config.drop_below
            pulled = np.asarray(
                pulled * (~sea_mask)[None, ...], dtype=np.complex128
            )
        pulled_energy = mean_energy(pulled)
        dropped_fraction = (
            (pulled_energy_full - pulled_energy) / pulled_energy_full
            if pulled_energy_full > 0.0
            else 0.0
        )
        gain = pulled_energy / initial_energy
        if gain > 0.0:
            noise_floor_estimate /= math.sqrt(gain)
        if not gain > 1.0e-20:
            terminated = "pullback_empty"
            records.append(
                CascadeStageRecord(
                    stage=stage,
                    steps=steps,
                    evolved_energy_ratio=evolved_energy / initial_energy,
                    child_band_energy_fraction=child_energy / initial_energy,
                    pullback_gain=gain if gain > 0.0 else 0.0,
                    dropped_sea_fraction=dropped_fraction,
                    overlap_lag1=0.0,
                    overlap_lag2=0.0,
                    overlap_lag3=0.0,
                    shell_energies=linf_shell_energies(
                        evolved,
                        (config.scale, 2 * config.scale),
                    ),
                    sup_velocity=_sup_velocity(evolved),
                    reality_defect=_reality_defect(evolved),
                    divergence_defect=_divergence_defect(evolved),
                )
            )
            break
        renormalized = np.asarray(
            pulled * math.sqrt(target_energy / pulled_energy),
            dtype=np.complex128,
        )
        overlaps = []
        for lag in (1, 2, 3):
            if len(history) >= lag:
                overlaps.append(
                    _shape_overlap(renormalized, history[-lag])
                )
            else:
                overlaps.append(0.0)
        records.append(
            CascadeStageRecord(
                stage=stage,
                steps=steps,
                evolved_energy_ratio=evolved_energy / initial_energy,
                child_band_energy_fraction=child_energy / initial_energy,
                pullback_gain=gain,
                dropped_sea_fraction=dropped_fraction,
                overlap_lag1=overlaps[0],
                overlap_lag2=overlaps[1],
                overlap_lag3=overlaps[2],
                shell_energies=linf_shell_energies(
                    renormalized,
                    (config.scale, 2 * config.scale),
                ),
                sup_velocity=_sup_velocity(renormalized),
                reality_defect=_reality_defect(renormalized),
                divergence_defect=_divergence_defect(renormalized),
            )
        )
        state = renormalized
        history.append(state.copy())
        if len(history) > 4:
            history.pop(0)

    return CascadeResult(
        config=config,
        completed_stages=len(records),
        terminated_reason=terminated,
        stage_records=tuple(records),
    )
