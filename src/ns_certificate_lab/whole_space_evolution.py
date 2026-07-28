r"""Small-amplitude nonlinear time evolution on the whole space.

This is the first nonlinear use of the Gate-4 elliptic path.  It deliberately
does **not** run a Hou-amplitude candidate: the point is to establish that the
non-periodic-``z`` free-space elliptic recovery composes correctly with a
time-stepper, and to measure convergence in space, time, domain and integrator
before any strongly nonlinear run is attempted.

Equations (audited entries E-11--E-14 of ``docs/equation_audit.md``):

.. math::

   \partial_t u_1 + u^r\partial_r u_1 + u^z\partial_z u_1
       = 2u_1\partial_z\psi_1 + \nu\mathcal L_5 u_1, \\
   \partial_t \omega_1 + u^r\partial_r \omega_1 + u^z\partial_z \omega_1
       = \partial_z(u_1^2) + \nu\mathcal L_5 \omega_1, \\
   -\mathcal L_5\psi_1 = \omega_1,\quad
   u^r = -r\partial_z\psi_1,\quad u^z = 2\psi_1 + r\partial_r\psi_1 .

Initial data is **pure swirl**: ``u_1(0)`` is a smooth compactly supported
axis-even bump and ``omega_1(0) = 0``.  The corresponding three-dimensional
field ``u = r u_1 \hat e_\theta`` is exactly divergence free, smooth, and
compactly supported, so it is admissible Clay data restricted to a bounded set;
``omega_1`` is then generated purely by the audited vortex-stretching term
``\partial_z(u_1^2)``.  That is the same mechanism as the Hou setup, run at an
amplitude where nothing is expected to happen.

Boundary treatment.  ``u_1`` and ``omega_1`` are pinned to zero on the outer
boundary.  That is exact for the initial datum and is only an approximation once
the support spreads, so ``outer_band_fraction`` measures how much of the field
has reached the outer band and every run reports it.  The elliptic solve uses
the multipole trace of :mod:`ns_certificate_lab.free_space_recovery`.

Nothing here is a proof.  Every number is binary64 without outward rounding.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
import math
from typing import Callable, Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .operators import (
    derivative_r,
    derivative_z,
    laplacian_5d_formal,
)
from .whole_space_gate import WholeSpaceEllipticSolver

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "EvolutionSnapshot",
    "EvolutionSummary",
    "SwirlBump",
    "WholeSpaceEvolution",
    "INTEGRATORS",
]

#: Butcher-style coefficients are inlined; these are the supported names.
INTEGRATORS = ("heun", "ssprk3", "rk4")


@dataclass(frozen=True)
class SwirlBump:
    """Smooth compactly supported axis-even pure-swirl initial datum.

    ``u_1(r,z) = amplitude (1 - (R/radius)^2)^power`` inside the five-dimensional
    ball ``R = sqrt(r^2 + (z-center)^2) <= radius``, zero outside.  It is
    ``C^{power-1}``, even in ``r``, and compactly supported in both directions.
    """

    center: float = 0.0
    radius: float = 1.0
    amplitude: float = 1.0
    power: int = 6

    def __post_init__(self) -> None:
        if not math.isfinite(self.radius) or self.radius <= 0.0:
            raise ValueError("radius must be positive and finite")
        if self.power < 2:
            raise ValueError("power must be at least 2 so the datum is C^1")

    def evaluate(self, r: FloatArray, z: FloatArray) -> FloatArray:
        radius5 = np.hypot(r, z - self.center)
        s = radius5 / self.radius
        out = np.zeros_like(radius5)
        inside = s < 1.0
        out[inside] = self.amplitude * (1.0 - s[inside] ** 2) ** self.power
        return out


@dataclass(frozen=True)
class EvolutionSnapshot:
    """Diagnostics at one recorded time."""

    time: float
    step: int
    kinetic_energy: float
    enstrophy_proxy: float
    l3_norm: float
    shell_fractions: tuple[float, ...]
    effective_shell_count: float
    max_relative_divergence: float
    max_curl_defect: float
    poisson_residual_max: float
    derivative_tail_bound: float
    velocity_tail_bound: float
    boundary_sensitivity: float
    outer_band_fraction: float
    low_frequency_fraction: float
    neglected_source_fraction: float
    max_u1: float
    max_omega1: float
    max_speed: float

    def as_dict(self) -> dict[str, object]:
        return {
            "time": self.time,
            "step": self.step,
            "kinetic_energy": self.kinetic_energy,
            "enstrophy_proxy": self.enstrophy_proxy,
            "l3_norm": self.l3_norm,
            "shell_fractions": list(self.shell_fractions),
            "effective_shell_count": self.effective_shell_count,
            "max_relative_divergence": self.max_relative_divergence,
            "max_curl_defect": self.max_curl_defect,
            "poisson_residual_max": self.poisson_residual_max,
            "derivative_tail_bound": self.derivative_tail_bound,
            "velocity_tail_bound": self.velocity_tail_bound,
            "boundary_sensitivity": self.boundary_sensitivity,
            "outer_band_fraction": self.outer_band_fraction,
            "low_frequency_fraction": self.low_frequency_fraction,
            "neglected_source_fraction": self.neglected_source_fraction,
            "max_u1": self.max_u1,
            "max_omega1": self.max_omega1,
            "max_speed": self.max_speed,
        }


@dataclass(frozen=True)
class EvolutionSummary:
    """One completed run: snapshots plus every-step streaming maxima."""

    integrator: str
    steps: int
    time_step: float
    final_time: float
    viscosity: float
    r_max: float
    z_max: float
    nr: int
    nz: int
    boundary_mode: str
    snapshots: tuple[EvolutionSnapshot, ...]
    streaming_max_u1: float
    streaming_max_omega1: float
    streaming_max_speed: float
    streaming_max_relative_divergence: float
    streaming_max_curl_defect: float
    streaming_max_poisson_residual: float
    streaming_max_outer_band_fraction: float
    streaming_max_energy_increase: float
    diverged: bool

    @property
    def final(self) -> EvolutionSnapshot:
        return self.snapshots[-1]

    def as_dict(self) -> dict[str, object]:
        return {
            "integrator": self.integrator,
            "steps": self.steps,
            "time_step": self.time_step,
            "final_time": self.final_time,
            "viscosity": self.viscosity,
            "r_max": self.r_max,
            "z_max": self.z_max,
            "nr": self.nr,
            "nz": self.nz,
            "boundary_mode": self.boundary_mode,
            "snapshots": [snapshot.as_dict() for snapshot in self.snapshots],
            "streaming_max_u1": self.streaming_max_u1,
            "streaming_max_omega1": self.streaming_max_omega1,
            "streaming_max_speed": self.streaming_max_speed,
            "streaming_max_relative_divergence": (
                self.streaming_max_relative_divergence
            ),
            "streaming_max_curl_defect": self.streaming_max_curl_defect,
            "streaming_max_poisson_residual": self.streaming_max_poisson_residual,
            "streaming_max_outer_band_fraction": (
                self.streaming_max_outer_band_fraction
            ),
            "streaming_max_energy_increase": self.streaming_max_energy_increase,
            "diverged": self.diverged,
        }


class WholeSpaceEvolution:
    """Fixed-step nonlinear evolution on a non-periodic-``z`` free-space box."""

    def __init__(
        self,
        grid: AxisymmetricGrid,
        *,
        viscosity: float,
        boundary_mode: str = "monopole",
        outer_band_cells: int = 4,
        support_threshold: float = 1.0e-10,
    ) -> None:
        if grid.periodic_z:
            raise ValueError("the whole-space evolution requires a non-periodic z grid")
        if viscosity < 0.0:
            raise ValueError("viscosity must be nonnegative")
        if not 1 <= outer_band_cells < min(grid.nr, grid.nz) // 2:
            raise ValueError("outer_band_cells must select a proper outer band")
        self.grid = grid
        self.viscosity = float(viscosity)
        self.boundary_mode = boundary_mode
        self.outer_band_cells = int(outer_band_cells)
        # Viscous diffusion makes |omega1| positive everywhere at roundoff, so
        # exact support detection would report the whole box and no multipole
        # tail bound would apply.  The threshold is a declared choice and the
        # neglected mass it drops is reported in every snapshot.
        self.support_threshold = float(support_threshold)
        self.solver = WholeSpaceEllipticSolver.build(grid, boundary_mode=boundary_mode)
        self._alt_solver = WholeSpaceEllipticSolver.build(
            grid, boundary_mode="zero" if boundary_mode == "monopole" else "monopole"
        )
        r_mesh, z_mesh = grid.mesh()
        self._r = r_mesh
        self._z = z_mesh
        self._volume = 2.0 * math.pi * grid.r[:, None] * np.ones_like(r_mesh)
        self._volume = self._volume * grid.dr * grid.dz

    # -- boundary bookkeeping ---------------------------------------------- #

    def _pin(self, field: FloatArray) -> FloatArray:
        out = field.copy()
        out[-1, :] = 0.0
        out[:, 0] = 0.0
        out[:, -1] = 0.0
        return out

    # -- the right-hand side ------------------------------------------------ #

    def velocity(self, omega1: FloatArray, *, solver=None):
        """Return ``(psi1, psi1_r, psi1_z, u^r, u^z)`` for the given vorticity."""
        active = self.solver if solver is None else solver
        psi = active.solve(omega1)
        psi_r = derivative_r(self.grid, psi, even_at_axis=True)
        psi_z = derivative_z(self.grid, psi)
        u_r = -self._r * psi_z
        u_r[0] = 0.0
        u_z = 2.0 * psi + self._r * psi_r
        return psi, psi_r, psi_z, u_r, u_z

    def rhs(
        self, u1: FloatArray, omega1: FloatArray, *, solver=None
    ) -> tuple[FloatArray, FloatArray]:
        grid = self.grid
        _psi, _psi_r, psi_z, u_r, u_z = self.velocity(omega1, solver=solver)
        u1_r = derivative_r(grid, u1, even_at_axis=True)
        u1_z = derivative_z(grid, u1)
        w_r = derivative_r(grid, omega1, even_at_axis=True)
        w_z = derivative_z(grid, omega1)
        du1 = (
            -u_r * u1_r
            - u_z * u1_z
            + 2.0 * u1 * psi_z
            + self.viscosity * laplacian_5d_formal(grid, u1)
        )
        dw = (
            -u_r * w_r
            - u_z * w_z
            + derivative_z(grid, u1 * u1)
            + self.viscosity * laplacian_5d_formal(grid, omega1)
        )
        return self._pin(du1), self._pin(dw)

    # -- integrators --------------------------------------------------------- #

    def step(
        self, u1: FloatArray, omega1: FloatArray, dt: float, integrator: str
    ) -> tuple[FloatArray, FloatArray]:
        if integrator not in INTEGRATORS:
            raise ValueError(f"integrator must be one of {INTEGRATORS}")
        k1u, k1w = self.rhs(u1, omega1)
        if integrator == "heun":
            pu, pw = u1 + dt * k1u, omega1 + dt * k1w
            k2u, k2w = self.rhs(pu, pw)
            return (
                self._pin(u1 + 0.5 * dt * (k1u + k2u)),
                self._pin(omega1 + 0.5 * dt * (k1w + k2w)),
            )
        if integrator == "ssprk3":
            a1u, a1w = u1 + dt * k1u, omega1 + dt * k1w
            k2u, k2w = self.rhs(a1u, a1w)
            a2u = 0.75 * u1 + 0.25 * (a1u + dt * k2u)
            a2w = 0.75 * omega1 + 0.25 * (a1w + dt * k2w)
            k3u, k3w = self.rhs(a2u, a2w)
            return (
                self._pin(u1 / 3.0 + (2.0 / 3.0) * (a2u + dt * k3u)),
                self._pin(omega1 / 3.0 + (2.0 / 3.0) * (a2w + dt * k3w)),
            )
        k2u, k2w = self.rhs(u1 + 0.5 * dt * k1u, omega1 + 0.5 * dt * k1w)
        k3u, k3w = self.rhs(u1 + 0.5 * dt * k2u, omega1 + 0.5 * dt * k2w)
        k4u, k4w = self.rhs(u1 + dt * k3u, omega1 + dt * k3w)
        return (
            self._pin(u1 + dt / 6.0 * (k1u + 2.0 * k2u + 2.0 * k3u + k4u)),
            self._pin(omega1 + dt / 6.0 * (k1w + 2.0 * k2w + 2.0 * k3w + k4w)),
        )

    # -- diagnostics ---------------------------------------------------------- #

    def _kinetic_energy(self, u1: FloatArray, u_r: FloatArray, u_z: FloatArray) -> float:
        speed_squared = u_r**2 + u_z**2 + (self._r * u1) ** 2
        return 0.5 * float(np.sum(speed_squared * self._volume))

    def _consistency_checks(
        self,
        omega1: FloatArray,
        psi_z: FloatArray,
        u_r: FloatArray,
        u_z: FloatArray,
    ) -> tuple[float, float]:
        """Relative physical divergence and the curl consistency defect.

        The divergence is the audited three-dimensional form
        ``d_r u^r + u^r/r + d_z u^z`` (E-15).  On the axis ``u^r/r`` is replaced
        by its exact continuous extension ``-psi_z``, so nothing divides by
        ``r``.  The curl check compares ``d_z u^r - d_r u^z`` against the
        audited identity ``omega^theta = r omega_1``; the two sides share no
        code, because the left side differentiates the recovered velocity and
        the right side is the elliptic source itself.
        """
        grid = self.grid
        dr_ur = derivative_r(grid, u_r)
        dz_uz = derivative_z(grid, u_z)
        divergence = dr_ur + dz_uz
        divergence[1:] += u_r[1:] / grid.r[1:, None]
        divergence[0] += -psi_z[0]
        scale = max(
            float(np.max(np.abs(dr_ur))), float(np.max(np.abs(dz_uz))), 1.0e-300
        )
        interior = (slice(1, -1), slice(1, -1))
        relative_divergence = float(np.max(np.abs(divergence[interior])) / scale)
        curl_theta = derivative_z(grid, u_r) - derivative_r(grid, u_z)
        expected = self._r * omega1
        curl_scale = max(float(np.max(np.abs(expected))), 1.0e-300)
        curl_defect = float(
            np.max(np.abs((curl_theta - expected)[interior])) / curl_scale
        )
        return relative_divergence, curl_defect

    def _outer_band_fraction(self, u1: FloatArray, omega1: FloatArray) -> float:
        band = self.outer_band_cells
        magnitude = np.abs(u1) + np.abs(omega1)
        total = float(np.sum(magnitude * self._volume))
        if total <= 0.0:
            return 0.0
        mask = np.zeros_like(magnitude, dtype=bool)
        mask[-band:, :] = True
        mask[:, :band] = True
        mask[:, -band:] = True
        return float(np.sum(magnitude * self._volume * mask)) / total

    def _low_frequency_fraction(self, omega1: FloatArray, *, share: float = 0.1) -> float:
        """Fraction of the axial sine spectrum energy in the lowest modes."""
        interior = omega1[:, 1:-1]
        spectrum = interior @ self.solver.sine.T
        power = np.sum(spectrum**2, axis=0)
        total = float(np.sum(power))
        if total <= 0.0:
            return 0.0
        cut = max(1, int(round(share * power.size)))
        return float(np.sum(power[:cut])) / total

    # -- the run -------------------------------------------------------------- #

    def run(
        self,
        u1_initial: FloatArray,
        omega1_initial: FloatArray,
        *,
        steps: int,
        time_step: float,
        integrator: str = "rk4",
        snapshot_every: int = 0,
        shell_count: int = 8,
    ) -> EvolutionSummary:
        from .critical_l3 import analyze_critical_l3

        if steps <= 0 or time_step <= 0.0:
            raise ValueError("steps and time_step must be positive")
        grid = self.grid
        u1 = self._pin(grid.validate_field(u1_initial, name="u1"))
        omega1 = self._pin(grid.validate_field(omega1_initial, name="omega1"))
        stride = snapshot_every if snapshot_every > 0 else max(1, steps // 4)

        from .free_space_recovery import recover_free_space_velocity

        snapshots: list[EvolutionSnapshot] = []
        streaming = {
            "u1": 0.0,
            "omega1": 0.0,
            "speed": 0.0,
            "divergence": 0.0,
            "curl": 0.0,
            "poisson": 0.0,
            "outer": 0.0,
            "energy_increase": 0.0,
        }
        previous_energy: float | None = None
        diverged = False

        def record(step: int, time: float) -> EvolutionSnapshot:
            nonlocal previous_energy
            recovery = recover_free_space_velocity(
                self.solver, omega1, support_threshold=self.support_threshold
            )
            u_r, u_z = recovery.u_r, recovery.u_z
            energy = self._kinetic_energy(u1, u_r, u_z)
            if previous_energy is not None:
                streaming["energy_increase"] = max(
                    streaming["energy_increase"],
                    (energy - previous_energy) / max(previous_energy, 1.0e-300),
                )
            previous_energy = energy
            relative_divergence, curl_defect = self._consistency_checks(
                omega1, recovery.dpsi1_dz, u_r, u_z
            )
            report = analyze_critical_l3(
                grid, u_r, self._r * u1, u_z, shell_count=shell_count
            )
            alt = self.rhs(u1, omega1, solver=self._alt_solver)
            base = self.rhs(u1, omega1)
            sensitivity = max(
                float(np.max(np.abs(alt[0] - base[0]))),
                float(np.max(np.abs(alt[1] - base[1]))),
            ) / max(
                float(np.max(np.abs(base[0]))), float(np.max(np.abs(base[1]))), 1.0e-300
            )
            outer = self._outer_band_fraction(u1, omega1)
            speed = float(np.max(np.sqrt(u_r**2 + u_z**2 + (self._r * u1) ** 2)))
            streaming["u1"] = max(streaming["u1"], float(np.max(np.abs(u1))))
            streaming["omega1"] = max(streaming["omega1"], float(np.max(np.abs(omega1))))
            streaming["speed"] = max(streaming["speed"], speed)
            streaming["divergence"] = max(streaming["divergence"], relative_divergence)
            streaming["curl"] = max(streaming["curl"], curl_defect)
            streaming["poisson"] = max(
                streaming["poisson"], recovery.poisson_residual_max
            )
            streaming["outer"] = max(streaming["outer"], outer)
            return EvolutionSnapshot(
                time=time,
                step=step,
                kinetic_energy=energy,
                enstrophy_proxy=float(
                    np.sum((self._r * omega1) ** 2 * self._volume)
                ),
                l3_norm=report.l3_norm,
                shell_fractions=report.shell_fractions,
                effective_shell_count=report.effective_shell_count,
                max_relative_divergence=relative_divergence,
                max_curl_defect=curl_defect,
                poisson_residual_max=recovery.poisson_residual_max,
                derivative_tail_bound=recovery.interior_tail.gradient,
                velocity_tail_bound=recovery.velocity_tail_bound,
                boundary_sensitivity=sensitivity,
                outer_band_fraction=outer,
                low_frequency_fraction=self._low_frequency_fraction(omega1),
                neglected_source_fraction=recovery.moments.neglected_fraction,
                max_u1=float(np.max(np.abs(u1))),
                max_omega1=float(np.max(np.abs(omega1))),
                max_speed=speed,
            )

        snapshots.append(record(0, 0.0))
        for index in range(steps):
            u1, omega1 = self.step(u1, omega1, time_step, integrator)
            if not (np.all(np.isfinite(u1)) and np.all(np.isfinite(omega1))):
                diverged = True
                break
            streaming["u1"] = max(streaming["u1"], float(np.max(np.abs(u1))))
            streaming["omega1"] = max(streaming["omega1"], float(np.max(np.abs(omega1))))
            if (index + 1) % stride == 0 or index == steps - 1:
                snapshots.append(record(index + 1, (index + 1) * time_step))

        return EvolutionSummary(
            integrator=integrator,
            steps=steps,
            time_step=time_step,
            final_time=steps * time_step,
            viscosity=self.viscosity,
            r_max=float(grid.r[-1]),
            z_max=float(grid.z[-1]),
            nr=grid.nr,
            nz=grid.nz,
            boundary_mode=self.boundary_mode,
            snapshots=tuple(snapshots),
            streaming_max_u1=streaming["u1"],
            streaming_max_omega1=streaming["omega1"],
            streaming_max_speed=streaming["speed"],
            streaming_max_relative_divergence=streaming["divergence"],
            streaming_max_curl_defect=streaming["curl"],
            streaming_max_poisson_residual=streaming["poisson"],
            streaming_max_outer_band_fraction=streaming["outer"],
            streaming_max_energy_increase=streaming["energy_increase"],
            diverged=diverged,
        )


def initial_state(
    grid: AxisymmetricGrid, bump: SwirlBump
) -> tuple[FloatArray, FloatArray]:
    """Pure-swirl initial data: ``u_1`` from the bump, ``omega_1 = 0``."""
    r_mesh, z_mesh = grid.mesh()
    u1 = bump.evaluate(r_mesh, z_mesh)
    u1[-1, :] = 0.0
    u1[:, 0] = 0.0
    u1[:, -1] = 0.0
    return u1, np.zeros_like(u1)


def _bilinear_sample(
    grid: AxisymmetricGrid, field: FloatArray, radius: FloatArray, axial: FloatArray
) -> FloatArray:
    """Bilinear sampling of a grid field; out-of-range points return zero."""
    r_index = (radius - grid.r[0]) / grid.dr
    z_index = (axial - grid.z[0]) / grid.dz
    inside = (
        (r_index >= 0.0)
        & (r_index <= grid.nr - 1)
        & (z_index >= 0.0)
        & (z_index <= grid.nz - 1)
    )
    i0 = np.clip(np.floor(r_index).astype(int), 0, grid.nr - 2)
    j0 = np.clip(np.floor(z_index).astype(int), 0, grid.nz - 2)
    tr = np.clip(r_index - i0, 0.0, 1.0)
    tz = np.clip(z_index - j0, 0.0, 1.0)
    value = (
        field[i0, j0] * (1.0 - tr) * (1.0 - tz)
        + field[i0 + 1, j0] * tr * (1.0 - tz)
        + field[i0, j0 + 1] * (1.0 - tr) * tz
        + field[i0 + 1, j0 + 1] * tr * tz
    )
    return np.where(inside, value, 0.0)


def cartesian_audit_of_state(
    grid: AxisymmetricGrid,
    u1: FloatArray,
    omega1: FloatArray,
    psi1: FloatArray,
    *,
    half_width: float,
    points: int,
    interior_trim: int = 2,
) -> dict[str, float]:
    r"""Independent Cartesian divergence and curl audit of a numerical state.

    The regularized component map (E-18) is

    .. math::

        u_x = -x\psi_{1,z} - y u_1,\quad
        u_y = -y\psi_{1,z} + x u_1,\quad
        u_z = 2\psi_1 + r\psi_{1,r},

    and the audited vorticity components for an axisymmetric field with swirl
    give the Cartesian curl

    .. math::

        \omega_x = -x\,\partial_z u_1 - y\,\omega_1,\quad
        \omega_y = -y\,\partial_z u_1 + x\,\omega_1,\quad
        \omega_z = 2u_1 + r\,\partial_r u_1 .

    Both are checked with the Cartesian stencils of
    :mod:`ns_certificate_lab.cartesian_validation`, which never call a
    cylindrical operator.  Bilinear sampling from the cylindrical grid is the
    one shared step and it is the dominant error source, so the numbers below
    are consistency evidence, not a bound.
    """
    from .cartesian_validation import (
        UniformCartesianGrid,
        cartesian_curl,
        cartesian_divergence,
    )

    if points < 8:
        raise ValueError("at least eight points per axis are required")
    psi_r = derivative_r(grid, psi1, even_at_axis=True)
    psi_z = derivative_z(grid, psi1)
    u1_r = derivative_r(grid, u1, even_at_axis=True)
    u1_z = derivative_z(grid, u1)

    axis = np.linspace(-half_width, half_width, points, dtype=np.float64)
    cartesian = UniformCartesianGrid(x=axis, y=axis, z=axis)
    x, y, z = cartesian.mesh()
    radius = np.hypot(x, y)

    def sample(field: FloatArray) -> FloatArray:
        return _bilinear_sample(grid, field, radius, z)

    swirl = sample(u1)
    velocity = np.stack(
        (
            -x * sample(psi_z) - y * swirl,
            -y * sample(psi_z) + x * swirl,
            2.0 * sample(psi1) + radius * sample(psi_r),
        ),
        axis=0,
    )
    divergence = cartesian_divergence(cartesian, velocity)
    curl = cartesian_curl(cartesian, velocity)
    expected = np.stack(
        (
            -x * sample(u1_z) - y * sample(omega1),
            -y * sample(u1_z) + x * sample(omega1),
            2.0 * swirl + radius * sample(u1_r),
        ),
        axis=0,
    )
    trim = slice(interior_trim, -interior_trim)
    core = (trim, trim, trim)
    velocity_scale = max(float(np.max(np.abs(velocity[(slice(None), *core)]))), 1.0e-300)
    curl_scale = max(float(np.max(np.abs(expected[(slice(None), *core)]))), 1.0e-300)
    return {
        "points": points,
        "velocity_max": velocity_scale,
        "divergence_max": float(np.max(np.abs(divergence[core]))),
        "divergence_relative": float(np.max(np.abs(divergence[core])))
        * float(axis[1] - axis[0])
        / velocity_scale,
        "curl_error_max": float(np.max(np.abs((curl - expected)[(slice(None), *core)]))),
        "curl_relative": float(
            np.max(np.abs((curl - expected)[(slice(None), *core)]))
        )
        / curl_scale,
    }
