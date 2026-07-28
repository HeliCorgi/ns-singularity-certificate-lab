r"""A finite dyadic cascade model for the role of a smooth force.

``docs/research_notes/track_f_shell_constraints.md`` §6 shows that a Clay-smooth
force has Fourier coefficients decaying faster than any polynomial, so the
*direct* injection term ``<u_j, f_j>`` is invisible in the high-shell budget.
That statement was previously over-read as "forcing gives Track F no advantage".
It does not: the force can shape the **low** shells, and the nonlinear triad
transfer can carry that shaping upward.  Deciding whether the indirect route can
sustain a growing bandwidth is a dynamical question, not a corollary of a decay
estimate.

This module is the smallest honest instrument for that question: the classical
dyadic (Desnyansky--Novikov / Katz--Pavlovic) shell model

.. math::

   \dot a_j = k_{j-1} a_{j-1}^2 - k_j a_j a_{j+1} - \nu k_j^2 a_j + f_j,
   \qquad k_j = k_0 2^j,\ j = 0..N-1,

with ``a_{-1} = a_N = 0``.  Its nonlinear part is **exactly energy conserving**:
``sum_j a_j (k_{j-1}a_{j-1}^2 - k_j a_j a_{j+1})`` telescopes to zero, which is
the shell-model shadow of the trilinear cancellation formalized in Lean as
``NSSingularity.advectionForm_eq_zero`` (F-12).  Each ``j`` couples ``(j-1,j)``
and ``(j,j+1)``, so ``N >= 5`` gives at least four connected triads.

Scope, stated before any number is produced.
* **This is not Navier--Stokes and not a Clay counterexample.**  A shell model
  discards all geometry, all phase information and all spatial structure.  It
  can only answer design-feasibility questions of the form "is the mechanism
  arithmetically possible", never "does it happen in the PDE".
* Its inviscid version is known to blow up in finite time, which is precisely
  why it is *useless* as evidence for blow-up and *useful* as a test of whether
  low-mode forcing can drive bandwidth growth against viscosity.
* Nothing here is formalized, and nothing here is a bound on the PDE.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "CascadeModel",
    "CascadeRun",
    "ForcingRole",
    "decompose_forcing_roles",
]


@dataclass(frozen=True)
class CascadeRun:
    """Streaming diagnostics of one dyadic-cascade integration."""

    shells: int
    viscosity: float
    forced_shells: int
    forcing_amplitude: float
    steps: int
    time_step: float
    final_time: float
    initial_bandwidth: float
    final_bandwidth: float
    max_bandwidth: float
    final_active_shell: int
    max_active_shell: int
    energy_initial: float
    energy_final: float
    max_energy: float
    nonlinear_transfer_defect: float
    max_high_shell_amplitude: float
    direct_high_shell_injection: float
    diverged: bool

    @property
    def bandwidth_grew(self) -> bool:
        return self.max_active_shell > 0 and self.max_bandwidth > self.initial_bandwidth

    def as_dict(self) -> dict[str, object]:
        return {
            "shells": self.shells,
            "viscosity": self.viscosity,
            "forced_shells": self.forced_shells,
            "forcing_amplitude": self.forcing_amplitude,
            "steps": self.steps,
            "time_step": self.time_step,
            "final_time": self.final_time,
            "initial_bandwidth": self.initial_bandwidth,
            "final_bandwidth": self.final_bandwidth,
            "max_bandwidth": self.max_bandwidth,
            "final_active_shell": self.final_active_shell,
            "max_active_shell": self.max_active_shell,
            "energy_initial": self.energy_initial,
            "energy_final": self.energy_final,
            "max_energy": self.max_energy,
            "nonlinear_transfer_defect": self.nonlinear_transfer_defect,
            "max_high_shell_amplitude": self.max_high_shell_amplitude,
            "direct_high_shell_injection": self.direct_high_shell_injection,
            "bandwidth_grew": self.bandwidth_grew,
            "diverged": self.diverged,
        }


@dataclass(frozen=True)
class CascadeModel:
    """The dyadic shell model with low-mode forcing only."""

    shells: int
    viscosity: float
    base_wavenumber: float = 1.0
    forced_shells: int = 2
    forcing_amplitude: float = 0.0
    high_shell_start: int = 4

    def __post_init__(self) -> None:
        if self.shells < 5:
            raise ValueError("at least five shells are needed for four triads")
        if self.viscosity < 0.0:
            raise ValueError("viscosity must be nonnegative")
        if self.base_wavenumber <= 0.0:
            raise ValueError("base_wavenumber must be positive")
        if not 0 <= self.forced_shells <= self.shells:
            raise ValueError("forced_shells must lie in [0, shells]")
        if not 0 < self.high_shell_start < self.shells:
            raise ValueError("high_shell_start must select a nonempty high band")
        if self.forced_shells > self.high_shell_start:
            raise ValueError(
                "forced_shells must stay strictly below high_shell_start so the "
                "force never injects directly into the band being measured"
            )

    @property
    def wavenumbers(self) -> FloatArray:
        return self.base_wavenumber * 2.0 ** np.arange(self.shells, dtype=np.float64)

    @property
    def forcing(self) -> FloatArray:
        """Force supported strictly in the low shells; zero in the high band."""
        out = np.zeros(self.shells, dtype=np.float64)
        out[: self.forced_shells] = self.forcing_amplitude
        return out

    # -- dynamics ----------------------------------------------------------- #

    def nonlinear(self, a: FloatArray) -> FloatArray:
        """The energy-conserving triad transfer."""
        k = self.wavenumbers
        below = np.zeros_like(a)
        below[1:] = k[:-1] * a[:-1] ** 2
        above = np.zeros_like(a)
        above[:-1] = k[:-1] * a[:-1] * a[1:]
        return below - above

    def transfer_defect(self, a: FloatArray) -> float:
        """``|sum_j a_j N_j(a)|`` relative to the scale: exactly zero in theory."""
        value = float(np.dot(a, self.nonlinear(a)))
        scale = max(float(np.max(np.abs(a))) ** 3 * float(self.wavenumbers[-1]), 1.0e-300)
        return abs(value) / scale

    def rhs(self, a: FloatArray) -> FloatArray:
        k = self.wavenumbers
        return self.nonlinear(a) - self.viscosity * k**2 * a + self.forcing

    # -- diagnostics --------------------------------------------------------- #

    def bandwidth(self, a: FloatArray) -> float:
        """Spectral centroid ``sum_j j a_j^2 / sum_j a_j^2``: a smooth bandwidth."""
        weight = a**2
        total = float(np.sum(weight))
        if total <= 0.0:
            return 0.0
        return float(np.dot(np.arange(self.shells, dtype=np.float64), weight)) / total

    def active_shell(self, a: FloatArray, *, threshold: float = 1.0e-8) -> int:
        """Highest shell whose amplitude exceeds ``threshold`` times the peak."""
        peak = float(np.max(np.abs(a)))
        if peak <= 0.0:
            return 0
        active = np.nonzero(np.abs(a) > threshold * peak)[0]
        return int(active[-1]) if active.size else 0

    def run(
        self,
        initial: FloatArray,
        *,
        steps: int,
        time_step: float,
        threshold: float = 1.0e-8,
    ) -> CascadeRun:
        """Integrate with RK4 and stream the cascade diagnostics."""
        if steps <= 0 or time_step <= 0.0:
            raise ValueError("steps and time_step must be positive")
        a = np.asarray(initial, dtype=np.float64).copy()
        if a.shape != (self.shells,):
            raise ValueError(f"initial must have shape ({self.shells},)")

        energy0 = float(np.sum(a**2))
        bandwidth0 = self.bandwidth(a)
        max_bandwidth = bandwidth0
        max_active = self.active_shell(a, threshold=threshold)
        max_energy = energy0
        max_defect = self.transfer_defect(a)
        max_high = float(np.max(np.abs(a[self.high_shell_start :])))
        diverged = False

        with np.errstate(over="ignore", invalid="ignore"):
            for _ in range(steps):
                k1 = self.rhs(a)
                k2 = self.rhs(a + 0.5 * time_step * k1)
                k3 = self.rhs(a + 0.5 * time_step * k2)
                k4 = self.rhs(a + time_step * k3)
                a = a + time_step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
                if not np.all(np.isfinite(a)):
                    diverged = True
                    break
                max_bandwidth = max(max_bandwidth, self.bandwidth(a))
                max_active = max(max_active, self.active_shell(a, threshold=threshold))
                max_energy = max(max_energy, float(np.sum(a**2)))
                max_defect = max(max_defect, self.transfer_defect(a))
                max_high = max(
                    max_high, float(np.max(np.abs(a[self.high_shell_start :])))
                )

        return CascadeRun(
            shells=self.shells,
            viscosity=self.viscosity,
            forced_shells=self.forced_shells,
            forcing_amplitude=self.forcing_amplitude,
            steps=steps,
            time_step=time_step,
            final_time=steps * time_step,
            initial_bandwidth=bandwidth0,
            final_bandwidth=self.bandwidth(a) if not diverged else math.inf,
            max_bandwidth=max_bandwidth,
            final_active_shell=self.active_shell(a, threshold=threshold)
            if not diverged
            else self.shells - 1,
            max_active_shell=max_active,
            energy_initial=energy0,
            energy_final=float(np.sum(a**2)) if not diverged else math.inf,
            max_energy=max_energy,
            nonlinear_transfer_defect=max_defect,
            max_high_shell_amplitude=max_high,
            direct_high_shell_injection=float(
                np.max(np.abs(self.forcing[self.high_shell_start :]))
            ),
            diverged=diverged,
        )


@dataclass(frozen=True)
class ForcingRole:
    """The five channels a smooth force could in principle use."""

    name: str
    available_to_a_clay_force: bool
    reason: str

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "available_to_a_clay_force": self.available_to_a_clay_force,
            "reason": self.reason,
        }


def decompose_forcing_roles() -> tuple[ForcingRole, ...]:
    """Separate the roles a Clay-admissible force can and cannot play.

    The point of the decomposition is that exactly one channel is closed by the
    decay estimate of ``track_f_shell_constraints.md`` §6, and the other four
    are not.  Reading the closed channel as "forcing does not help" is the error
    this function exists to prevent.
    """
    return (
        ForcingRole(
            name="direct_high_frequency_injection",
            available_to_a_clay_force=False,
            reason=(
                "a C-infinity force has Fourier coefficients decaying faster "
                "than any polynomial, so <u_j, f_j> is O(lambda_j^-m) for every "
                "m and cannot feed a high shell directly"
            ),
        ),
        ForcingRole(
            name="low_frequency_injection",
            available_to_a_clay_force=True,
            reason=(
                "the decay estimate says nothing about O(1) wavenumbers, where "
                "the force is unconstrained apart from smoothness"
            ),
        ),
        ForcingRole(
            name="shell_flux_through_nonlinear_triads",
            available_to_a_clay_force=True,
            reason=(
                "energy placed in low shells is transported upward by the "
                "energy-conserving triad term; the force never appears in that "
                "term but it sets the amplitudes that drive it"
            ),
        ),
        ForcingRole(
            name="phase_and_alignment_control",
            available_to_a_clay_force=True,
            reason=(
                "the force can choose the relative phase and geometric "
                "alignment of the low modes, which the PDE nonlinearity is "
                "sensitive to; a shell model cannot even represent this channel"
            ),
        ),
        ForcingRole(
            name="sustaining_flux_above_viscous_damping",
            available_to_a_clay_force=True,
            reason=(
                "a force acting for all time can hold the low shells at an "
                "amplitude for which the triad flux exceeds nu k^2 damping in "
                "the shells above; whether that suffices in the PDE is open"
            ),
        ),
    )
