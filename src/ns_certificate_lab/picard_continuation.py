r"""Picard ladder and dimensionless-time continuation for the whole-space run.

Gate 6 discovered that thirty-two runs had all stayed inside the first Picard
iterate, diagnosed by the exact ``A^2`` scaling of ``max|omega_1|``.  Detecting
that after the fact is not enough: this module integrates the Picard iterates
**alongside** the full solution, so the distance to them is a measured quantity
at every checkpoint rather than an inference from a scaling law.

The ladder
----------
For the mild formulation, the Picard iteration keeps the viscous term implicit
and evaluates the nonlinear terms at the previous iterate:

.. math::

   \partial_t u^{(n+1)} = -u^{(n)}\!\cdot\!\nabla u^{(n)}_1
       + 2u^{(n)}_1\partial_z\psi^{(n)} + \nu\mathcal L_5 u^{(n+1)}, \\
   \partial_t \omega^{(n+1)} = -u^{(n)}\!\cdot\!\nabla \omega^{(n)}
       + \partial_z\bigl((u^{(n)}_1)^2\bigr) + \nu\mathcal L_5 \omega^{(n+1)} ,

both starting from the same datum.  Level ``0`` is pure diffusion, so its
vorticity stays zero and its velocity vanishes; level ``1`` is therefore the
familiar "vorticity grows linearly from the stretching source" picture, and
level ``2`` is the first iterate that feels advection.

Integrating all levels **and** the full nonlinear solution with the same
integrator and the same steps makes the comparison exact: no interpolation
between trajectories, no snapshot differencing.

Why ``tau = O(1)`` is not the answer
------------------------------------
Pure swirl advects nothing.  The meridional velocity that does the advecting is
generated at order ``tau`` through the elliptic solve, whose kernel carries the
small constant ``1/(8 pi^2)``.  The advection-to-source ratio therefore grows
like ``C tau^2`` with ``C`` of order ``10^{-3}`` for these shapes, so leaving the
first Picard iterate needs ``tau`` of order ``sqrt(0.1/C)``, not ``1``.  That is
measured here rather than assumed; see
:func:`ns_certificate_lab.nondimensional.nonlinear_time_estimate`.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .operators import derivative_r, derivative_z, laplacian_5d_formal
from .whole_space_gate import WholeSpaceEllipticSolver

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "ContinuationCheckpoint",
    "ContinuationResult",
    "PicardLadder",
    "DEPARTURE_THRESHOLDS",
    "PREREGISTERED_TAU",
    "PREREGISTERED_REYNOLDS",
]

#: Preregistered dimensionless-time checkpoints.  Recorded before any run.
PREREGISTERED_TAU: tuple[float, ...] = (0.025, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0)

#: Preregistered Reynolds-number series.  Recorded before any run.
PREREGISTERED_REYNOLDS: tuple[float, ...] = (10.0, 25.0, 50.0, 100.0, 200.0, 400.0)

#: Preregistered thresholds of the first-Picard-departure gate.
DEPARTURE_THRESHOLDS = {
    "relative_distance_to_first_picard": 0.10,
    "resolution_stability": 0.20,
    "front_points": 10.0,
}


@dataclass(frozen=True)
class ContinuationCheckpoint:
    """Every diagnostic recorded at one preregistered ``tau``."""

    tau: float
    physical_time: float
    steps: int
    accepted_time_step: float
    max_u1: float
    max_omega1: float
    max_speed: float
    distance_to_first_picard: float
    distance_to_second_picard: float
    advection_over_source: float
    viscous_over_source: float
    rhs_change_rate: float
    kinetic_energy: float
    enstrophy: float
    l3_norm: float
    core_l3_fraction: float
    shell_fractions: tuple[float, ...]
    effective_shell_count: float
    radial_scale: float
    axial_scale: float
    shell_flux_upward: float
    outer_band_fraction: float
    relative_divergence: float
    curl_defect: float
    poisson_residual: float
    front_points_radial: float
    front_points_axial: float

    def as_dict(self) -> dict[str, object]:
        payload = {
            key: getattr(self, key)
            for key in self.__dataclass_fields__
            if key != "shell_fractions"
        }
        payload["shell_fractions"] = list(self.shell_fractions)
        return payload


@dataclass(frozen=True)
class ContinuationResult:
    """One completed continuation run."""

    label: str
    family: str
    reynolds: float
    amplitude: float
    reference_length: float
    viscosity: float
    nr: int
    nz: int
    r_max: float
    z_max: float
    integrator: str
    checkpoints: tuple[ContinuationCheckpoint, ...]
    reached_tau: float
    stopped_reason: str

    @property
    def final(self) -> ContinuationCheckpoint:
        return self.checkpoints[-1]

    def as_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "family": self.family,
            "reynolds": self.reynolds,
            "amplitude": self.amplitude,
            "reference_length": self.reference_length,
            "viscosity": self.viscosity,
            "nr": self.nr,
            "nz": self.nz,
            "r_max": self.r_max,
            "z_max": self.z_max,
            "integrator": self.integrator,
            "reached_tau": self.reached_tau,
            "stopped_reason": self.stopped_reason,
            "checkpoints": [c.as_dict() for c in self.checkpoints],
        }


class PicardLadder:
    """Full solution plus Picard levels 0, 1 and 2, integrated together."""

    def __init__(
        self,
        grid: AxisymmetricGrid,
        *,
        viscosity: float,
        boundary_mode: str = "dipole",
        cfl: float = 0.25,
        diffusion_safety: float = 0.2,
        outer_band_cells: int = 4,
    ) -> None:
        if grid.periodic_z:
            raise ValueError("the whole-space ladder requires a non-periodic z grid")
        if viscosity <= 0.0:
            raise ValueError("a Clay candidate requires a fixed positive viscosity")
        self.grid = grid
        self.viscosity = float(viscosity)
        self.cfl = float(cfl)
        self.diffusion_safety = float(diffusion_safety)
        self.outer_band_cells = int(outer_band_cells)
        self.solver = WholeSpaceEllipticSolver.build(grid, boundary_mode=boundary_mode)
        r_mesh, z_mesh = grid.mesh()
        self._r = r_mesh
        self._z = z_mesh
        self._volume = 2.0 * math.pi * r_mesh * grid.dr * grid.dz

    # -- pieces --------------------------------------------------------------- #

    def _pin(self, field: FloatArray) -> FloatArray:
        out = field.copy()
        out[-1, :] = 0.0
        out[:, 0] = 0.0
        out[:, -1] = 0.0
        return out

    def velocity(self, omega1: FloatArray):
        psi = self.solver.solve(omega1)
        psi_r = derivative_r(self.grid, psi, even_at_axis=True)
        psi_z = derivative_z(self.grid, psi)
        u_r = -self._r * psi_z
        u_r[0] = 0.0
        u_z = 2.0 * psi + self._r * psi_r
        return psi, psi_z, u_r, u_z

    def _nonlinear(self, u1: FloatArray, omega1: FloatArray):
        """The terms the Picard iteration freezes at the previous level."""
        _psi, psi_z, u_r, u_z = self.velocity(omega1)
        grid = self.grid
        swirl = (
            -u_r * derivative_r(grid, u1, even_at_axis=True)
            - u_z * derivative_z(grid, u1)
            + 2.0 * u1 * psi_z
        )
        vorticity = (
            -u_r * derivative_r(grid, omega1, even_at_axis=True)
            - u_z * derivative_z(grid, omega1)
            + derivative_z(grid, u1 * u1)
        )
        return swirl, vorticity, u_r, u_z, psi_z

    def rhs(self, state: tuple[FloatArray, ...]) -> tuple[FloatArray, ...]:
        """``state = (u0,w0, u1,w1, u2,w2, uF,wF)``; returns the same layout."""
        grid = self.grid
        nu = self.viscosity
        u0, w0, u1, w1, u2, w2, uF, wF = state
        # level 0: pure diffusion
        d_u0 = nu * laplacian_5d_formal(grid, u0)
        d_w0 = nu * laplacian_5d_formal(grid, w0)
        # level 1: nonlinear terms from level 0
        s0, v0, _ur0, _uz0, _pz0 = self._nonlinear(u0, w0)
        d_u1 = s0 + nu * laplacian_5d_formal(grid, u1)
        d_w1 = v0 + nu * laplacian_5d_formal(grid, w1)
        # level 2: nonlinear terms from level 1
        s1, v1, _ur1, _uz1, _pz1 = self._nonlinear(u1, w1)
        d_u2 = s1 + nu * laplacian_5d_formal(grid, u2)
        d_w2 = v1 + nu * laplacian_5d_formal(grid, w2)
        # the full nonlinear solution
        sF, vF, _urF, _uzF, _pzF = self._nonlinear(uF, wF)
        d_uF = sF + nu * laplacian_5d_formal(grid, uF)
        d_wF = vF + nu * laplacian_5d_formal(grid, wF)
        return tuple(
            self._pin(field)
            for field in (d_u0, d_w0, d_u1, d_w1, d_u2, d_w2, d_uF, d_wF)
        )

    # -- time stepping -------------------------------------------------------- #

    def stable_time_step(self, state: tuple[FloatArray, ...]) -> float:
        """Adaptive step from the advective CFL and the diffusion limit."""
        grid = self.grid
        _psi, _pz, u_r, u_z = self.velocity(state[7])
        speed = max(
            float(np.max(np.abs(u_r))) / grid.dr,
            float(np.max(np.abs(u_z))) / grid.dz,
            1.0e-300,
        )
        advective = self.cfl / speed
        # L5 near the axis contributes 8/dr^2; the axial part 2/dz^2.
        diffusive = self.diffusion_safety / (
            self.viscosity * (8.0 / grid.dr**2 + 2.0 / grid.dz**2)
        )
        return min(advective, diffusive)

    def step(
        self, state: tuple[FloatArray, ...], dt: float, integrator: str
    ) -> tuple[FloatArray, ...]:
        def axpy(base, increment, factor):
            return tuple(b + factor * i for b, i in zip(base, increment))

        k1 = self.rhs(state)
        if integrator == "ssprk3":
            a1 = axpy(state, k1, dt)
            k2 = self.rhs(a1)
            a2 = tuple(
                0.75 * s + 0.25 * (a + dt * k)
                for s, a, k in zip(state, a1, k2)
            )
            k3 = self.rhs(a2)
            return tuple(
                self._pin(s / 3.0 + (2.0 / 3.0) * (a + dt * k))
                for s, a, k in zip(state, a2, k3)
            )
        if integrator != "rk4":
            raise ValueError("integrator must be 'ssprk3' or 'rk4'")
        k2 = self.rhs(axpy(state, k1, 0.5 * dt))
        k3 = self.rhs(axpy(state, k2, 0.5 * dt))
        k4 = self.rhs(axpy(state, k3, dt))
        return tuple(
            self._pin(s + dt / 6.0 * (a + 2.0 * b + 2.0 * c + d))
            for s, a, b, c, d in zip(state, k1, k2, k3, k4)
        )

    # -- diagnostics ---------------------------------------------------------- #

    def _front_points(self, field: FloatArray) -> tuple[float, float]:
        """Points across the 10--90 front of ``|field|`` in each direction."""
        magnitude = np.abs(field)
        peak = float(magnitude.max())
        if peak <= 0.0:
            return (0.0, 0.0)
        index = np.unravel_index(int(np.argmax(magnitude)), magnitude.shape)

        def count(profile: np.ndarray) -> float:
            high = profile >= 0.9 * peak
            low = profile >= 0.1 * peak
            return float(low.sum() - high.sum()) / 2.0

        return (
            count(magnitude[:, index[1]]),
            count(magnitude[index[0], :]),
        )

    def diagnostics(
        self,
        state: tuple[FloatArray, ...],
        *,
        tau: float,
        physical_time: float,
        steps: int,
        dt: float,
        previous_rhs: tuple[FloatArray, ...] | None,
        shell_count: int = 8,
    ) -> tuple[ContinuationCheckpoint, tuple[FloatArray, ...]]:
        from .critical_l3 import analyze_critical_l3
        from .free_space_recovery import recover_free_space_velocity

        grid = self.grid
        u0, _w0, _u1, w1, _u2, w2, uF, wF = state
        rhs = self.rhs(state)
        scale = max(float(np.max(np.abs(wF))), 1.0e-300)
        recovery = recover_free_space_velocity(
            self.solver, wF, support_threshold=1.0e-10
        )
        report = analyze_critical_l3(
            grid, recovery.u_r, self._r * uF, recovery.u_z, shell_count=shell_count
        )
        source = derivative_z(grid, uF * uF)
        source_scale = max(float(np.max(np.abs(source))), 1.0e-300)
        advection = -recovery.u_r * derivative_r(
            grid, wF, even_at_axis=True
        ) - recovery.u_z * derivative_z(grid, wF)
        viscous = self.viscosity * laplacian_5d_formal(grid, wF)
        rhs_change = (
            float(np.max(np.abs(rhs[7] - previous_rhs[7])))
            / max(float(np.max(np.abs(rhs[7]))), 1.0e-300)
            if previous_rhs is not None
            else 0.0
        )
        band = self.outer_band_cells
        magnitude = np.abs(uF) + np.abs(wF)
        total = float(np.sum(magnitude * self._volume)) or 1.0
        mask = np.zeros_like(magnitude, dtype=bool)
        mask[-band:, :] = True
        mask[:, :band] = True
        mask[:, -band:] = True
        core = np.hypot(self._r, self._z) <= 0.5 * min(
            float(grid.r[-1]), float(grid.z[-1])
        )
        speed = np.sqrt(recovery.u_r**2 + recovery.u_z**2 + (self._r * uF) ** 2)
        l3_density = speed**3 * self._volume
        core_fraction = float(np.sum(l3_density * core)) / max(
            float(np.sum(l3_density)), 1.0e-300
        )
        shells = report.shell_fractions
        upward = float(sum(shells[len(shells) // 2 :]))
        radial_points, axial_points = self._front_points(wF)
        divergence = derivative_r(grid, recovery.u_r) + derivative_z(grid, recovery.u_z)
        divergence[1:] += recovery.u_r[1:] / grid.r[1:, None]
        divergence[0] += -recovery.dpsi1_dz[0]
        divergence_scale = max(
            float(np.max(np.abs(derivative_r(grid, recovery.u_r)))), 1.0e-300
        )
        curl = derivative_z(grid, recovery.u_r) - derivative_r(grid, recovery.u_z)
        checkpoint = ContinuationCheckpoint(
            tau=tau,
            physical_time=physical_time,
            steps=steps,
            accepted_time_step=dt,
            max_u1=float(np.max(np.abs(uF))),
            max_omega1=float(np.max(np.abs(wF))),
            max_speed=float(np.max(speed)),
            distance_to_first_picard=float(np.max(np.abs(wF - w1))) / scale,
            distance_to_second_picard=float(np.max(np.abs(wF - w2))) / scale,
            advection_over_source=float(np.max(np.abs(advection))) / source_scale,
            viscous_over_source=float(np.max(np.abs(viscous))) / source_scale,
            rhs_change_rate=rhs_change,
            kinetic_energy=0.5
            * float(np.sum((recovery.u_r**2 + recovery.u_z**2 + (self._r * uF) ** 2)
                           * self._volume)),
            enstrophy=float(np.sum((self._r * wF) ** 2 * self._volume)),
            l3_norm=report.l3_norm,
            core_l3_fraction=core_fraction,
            shell_fractions=shells,
            effective_shell_count=report.effective_shell_count,
            radial_scale=report.radial_scale,
            axial_scale=report.axial_scale,
            shell_flux_upward=upward,
            outer_band_fraction=float(np.sum(magnitude * self._volume * mask)) / total,
            relative_divergence=float(
                np.max(np.abs(divergence[1:-1, 1:-1]))
            ) / divergence_scale,
            curl_defect=float(
                np.max(np.abs((curl - self._r * wF)[1:-1, 1:-1]))
            ) / max(float(np.max(np.abs(self._r * wF))), 1.0e-300),
            poisson_residual=recovery.poisson_residual_max,
            front_points_radial=radial_points,
            front_points_axial=axial_points,
        )
        return checkpoint, rhs

    # -- the continuation ------------------------------------------------------ #

    def run(
        self,
        u1_initial: FloatArray,
        *,
        amplitude: float,
        tau_checkpoints: Sequence[float] = PREREGISTERED_TAU,
        integrator: str = "rk4",
        max_steps: int = 200000,
        max_time_step: float | None = None,
    ) -> tuple[list[ContinuationCheckpoint], str, tuple[FloatArray, ...]]:
        """Advance to each preregistered ``tau`` in turn, stopping on any fault."""
        grid = self.grid
        u0 = self._pin(grid.validate_field(u1_initial, name="u1"))
        zero = np.zeros_like(u0)
        state = (u0, zero.copy(), u0.copy(), zero.copy(), u0.copy(), zero.copy(),
                 u0.copy(), zero.copy())
        checkpoints: list[ContinuationCheckpoint] = []
        previous_rhs: tuple[FloatArray, ...] | None = None
        time = 0.0
        steps = 0
        reason = "completed"
        for target in tau_checkpoints:
            target_time = target / amplitude
            while time < target_time:
                dt = self.stable_time_step(state)
                if max_time_step is not None:
                    dt = min(dt, max_time_step)
                dt = min(dt, target_time - time)
                if dt <= 0.0 or not math.isfinite(dt):
                    reason = "time step collapsed"
                    break
                state = self.step(state, dt, integrator)
                time += dt
                steps += 1
                if not all(np.all(np.isfinite(field)) for field in state):
                    reason = "non-finite state"
                    break
                if steps > max_steps:
                    reason = "step budget exhausted"
                    break
            if reason != "completed":
                break
            checkpoint, previous_rhs = self.diagnostics(
                state, tau=target, physical_time=time, steps=steps,
                dt=self.stable_time_step(state), previous_rhs=previous_rhs,
            )
            if checkpoint.max_omega1 <= 0.0:
                reason = "vorticity vanished"
                break
            if checkpoint.front_points_radial < 1.0:
                reason = "front unresolved"
                break
            checkpoints.append(checkpoint)
        return checkpoints, reason, state
