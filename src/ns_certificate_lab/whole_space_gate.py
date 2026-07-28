r"""Gate 4: the linear whole-space elliptic gate with non-periodic ``z``.

Scope
-----
``docs/whole_space_transition.md`` §7 and ``PLAN.md`` Phase 2.8 Gate 4 require,
before any whole-space *nonlinear* run, a **linear** gate that removes both
finite-domain artefacts of the existing evidence at once:

* the radial Dirichlet wall at ``r = R``, and
* the periodic ``z`` direction, whose period images were previously only
  probed by zero padding.

This module supplies that gate.  It solves

.. math::

   -L_5\psi_1 = \omega_1,\qquad
   L_5 = \partial_{rr} + \frac{3}{r}\partial_r + \partial_{zz}

on a finite box ``[0,R_{\max}] × [-Z_{\max},Z_{\max}]`` with **non-periodic**
``z``, against an **exact closed-form free-space reference**.

Why the reference is exact and non-circular
-------------------------------------------
``L_5`` is the Laplacian of :math:`\mathbb R^5` acting on functions of
:math:`(r,z)` with :math:`r=|x'|`, :math:`x'\in\mathbb R^4`.  Writing
:math:`R=\sqrt{r^2+z^2}` for the full five-dimensional radius, a function of
``R`` alone obeys :math:`L_5\psi = \psi'' + 4\psi'/R`.  Newton's theorem then
gives, for a source supported in ``R ≤ a``,

.. math::

   \psi_\infty(R)=\int_R^\infty s^{-4}m(s)\,ds,\qquad
   m(s)=\int_0^s t^4\omega_1(t)\,dt ,

and :math:`\psi_\infty(R)=M/(3R^3)` for ``R ≥ a`` with ``M = m(∞)``.  Choosing
:math:`\omega_1 = c\,(1-(R/a)^2)^p` makes every integral a finite elementary
sum, so the reference costs no quadrature and no discretisation.  Because
``L_5`` is invariant under ``z``-translation, superposing bumps centred at
different ``z_0`` keeps the reference exact while producing genuine
``z``-structure — including a balanced pair with ``M = 0``, whose far field
decays one order faster and which is the sharpest available probe of the
truncation error.

Independence of the axial path
------------------------------
Every existing solver in this repository diagonalises ``z`` with
``numpy.fft``.  ``docs/formalization_map.md`` and ``STATUS.md`` list the
absence of a **non-Fourier independent axial path** as an open gate.  Here the
axial direction carries homogeneous Dirichlet data at ``z = ±Z_max`` after the
boundary trace is moved to the right-hand side, so the interior second
difference is diagonalised by the **discrete sine transform of type I**,
applied as a dense orthogonal matrix built from ``numpy.sin`` — no FFT, no
shared wavenumber array, no shared convention with the periodic solvers.  Each
axial eigenvalue then gives a radial tridiagonal system solved by a locally
implemented Thomas sweep.

What this module is not
-----------------------
* It is a **linear** gate.  It says nothing about the nonlinear evolution and
  does not license a whole-space Hou run on its own.
* Every error it reports is binary64 without outward rounding.  The one
  statement with proof status is :func:`monopole_tail_bound`, and even there
  only the *continuum domain-truncation* part is bounded (by the maximum
  principle for the five-dimensional Laplacian); the discretisation error is
  measured, not bounded, because the radial row ``i=1`` of this stencil is not
  an M-matrix and carries no discrete maximum principle.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Iterable, Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "BumpSpec",
    "CartesianGateReport",
    "DomainConvergenceRow",
    "ManufacturedFreeSpaceField",
    "TailBoundReport",
    "WholeSpaceEllipticSolver",
    "AxialSpectrumReport",
    "SeparationReport",
    "TruncationRow",
    "audit_cartesian_recovery",
    "axial_spectrum",
    "cylindrical_axis_free_space_solid_angle",
    "domain_truncation_series",
    "ZeroFrequencyRow",
    "separate_period_and_wall_errors",
    "zero_frequency_stress",
    "solve_and_compare",
    "monopole_tail_bound",
    "source_monopole_moment",
    "source_l1_norm_5d",
]

#: ``|S^3|``, the surface measure of the unit sphere in ``ℝ^4``.
SPHERE3_MEASURE = 2.0 * math.pi**2
#: ``|S^4|``, the surface measure of the unit sphere in ``ℝ^5``.
SPHERE4_MEASURE = 8.0 * math.pi**2 / 3.0


def cylindrical_axis_free_space_solid_angle() -> tuple[float, float]:
    """Return ``(|S^3|, |S^4|)``, the two measures this module depends on."""
    return SPHERE3_MEASURE, SPHERE4_MEASURE


# --------------------------------------------------------------------------- #
# 1. Smooth compactly supported manufactured source with an exact free-space   #
#    solution on non-periodic ``z``                                            #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class BumpSpec:
    """One five-dimensional radial bump centred on the axis at ``z = center``.

    The source contribution is ``amplitude * (1 - (R/radius)^2)^power`` inside
    ``R ≤ radius`` and exactly zero outside, with
    ``R = sqrt(r^2 + (z-center)^2)``.  It is ``C^{power-1}`` and compactly
    supported in both ``r`` and ``z``.

    ``power`` is deliberately finite: a ``C^∞`` bump such as
    ``exp(-1/(1-s^2))`` has no elementary antiderivative, so its free-space
    solution would itself have to be quadratured and would stop being an
    independent reference.  ``power = 6`` gives ``C^5``, comfortably more than
    a second-order stencil needs.
    """

    center: float
    radius: float
    amplitude: float
    power: int = 6

    def __post_init__(self) -> None:
        if not math.isfinite(self.center):
            raise ValueError("center must be finite")
        if not math.isfinite(self.radius) or self.radius <= 0.0:
            raise ValueError("radius must be positive and finite")
        if not math.isfinite(self.amplitude):
            raise ValueError("amplitude must be finite")
        if self.power < 2:
            raise ValueError("power must be at least 2 so the source is C^1")

    @property
    def _binomials(self) -> list[float]:
        return [float(math.comb(self.power, j)) * (-1.0) ** j for j in range(self.power + 1)]

    @property
    def mass(self) -> float:
        """``m(∞) = ∫_0^∞ t^4 ω(t) dt`` for this bump."""
        total = sum(
            coefficient / (5.0 + 2.0 * j)
            for j, coefficient in enumerate(self._binomials)
        )
        return self.amplitude * self.radius**5 * total

    @property
    def l1_mass(self) -> float:
        """``∫_0^∞ t^4 |ω(t)| dt`` — the absolute moment used by the bound."""
        total = sum(
            coefficient / (5.0 + 2.0 * j)
            for j, coefficient in enumerate(self._binomials)
        )
        # (1-s^2)^p >= 0 on the support, so |ω| = |amplitude| (1-s^2)^p.
        return abs(self.amplitude) * self.radius**5 * total

    def support_radius_from_origin(self) -> float:
        """Radius of a five-dimensional ball about the origin covering the support."""
        return abs(self.center) + self.radius

    def omega(self, r: FloatArray, z: FloatArray) -> FloatArray:
        radius5 = np.hypot(np.asarray(r, dtype=np.float64), np.asarray(z, dtype=np.float64) - self.center)
        s = radius5 / self.radius
        inside = s < 1.0
        out = np.zeros_like(radius5)
        out[inside] = self.amplitude * (1.0 - s[inside] ** 2) ** self.power
        return out

    def _radial_potential(self, radius5: FloatArray) -> FloatArray:
        """``ψ(R)`` for this bump, exact and elementary."""
        s = radius5 / self.radius
        out = np.empty_like(radius5)
        far = s >= 1.0
        near = ~far
        mass = self.mass
        out[far] = mass / (3.0 * np.maximum(radius5[far], 1.0e-300) ** 3)
        if np.any(near):
            sn = s[near]
            inner = np.zeros_like(sn)
            for j, coefficient in enumerate(self._binomials):
                inner += (coefficient / (5.0 + 2.0 * j)) * (
                    1.0 - sn ** (2 + 2 * j)
                ) / (2.0 + 2.0 * j)
            out[near] = self.amplitude * self.radius**2 * inner + mass / (
                3.0 * self.radius**3
            )
        return out

    def _radial_potential_slope_over_radius(self, radius5: FloatArray) -> FloatArray:
        """``ψ'(R)/R = -m(R)/R^5``, regular at ``R = 0``."""
        s = radius5 / self.radius
        out = np.empty_like(radius5)
        far = s >= 1.0
        near = ~far
        out[far] = -self.mass / np.maximum(radius5[far], 1.0e-300) ** 5
        if np.any(near):
            sn = s[near]
            series = np.zeros_like(sn)
            for j, coefficient in enumerate(self._binomials):
                series += (coefficient / (5.0 + 2.0 * j)) * sn ** (2 * j)
            out[near] = -self.amplitude * series
        return out


@dataclass(frozen=True)
class ManufacturedFreeSpaceField:
    """A superposition of :class:`BumpSpec` with its exact free-space solution.

    ``L_5`` commutes with translations in ``z``, so the exact solution of the
    superposed source is the superposition of the shifted radial solutions.
    Nothing here is discretised: ``psi1``, ``dpsi1_dr`` and ``dpsi1_dz`` are
    closed-form and can be evaluated at arbitrary points, which is what makes
    them usable both as the solver reference and as the analytic profile of the
    independent Cartesian audit.
    """

    bumps: tuple[BumpSpec, ...]

    def __post_init__(self) -> None:
        if not self.bumps:
            raise ValueError("at least one bump is required")

    @property
    def mass(self) -> float:
        """``M = m(∞)``.  The far field is ``M/(3R^3)``."""
        return sum(bump.mass for bump in self.bumps)

    @property
    def monopole_moment_5d(self) -> float:
        """``∫ ω dV_5 = |S^4| M``."""
        return SPHERE4_MEASURE * self.mass

    @property
    def l1_norm_5d(self) -> float:
        """An upper bound for ``∫ |ω| dV_5`` (exact for a single bump)."""
        return SPHERE4_MEASURE * sum(bump.l1_mass for bump in self.bumps)

    @property
    def support_radius(self) -> float:
        """Radius of a five-dimensional ball about the origin covering the source."""
        return max(bump.support_radius_from_origin() for bump in self.bumps)

    def omega1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        total = np.zeros(np.broadcast(r_array, z_array).shape, dtype=np.float64)
        for bump in self.bumps:
            total = total + bump.omega(np.broadcast_to(r_array, total.shape),
                                       np.broadcast_to(z_array, total.shape))
        return total

    def psi1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        total = np.zeros(shape, dtype=np.float64)
        rb = np.broadcast_to(r_array, shape)
        zb = np.broadcast_to(z_array, shape)
        for bump in self.bumps:
            total = total + bump._radial_potential(np.hypot(rb, zb - bump.center))
        return total

    def _slope_over_radius(self, rb: FloatArray, zb: FloatArray) -> list[FloatArray]:
        return [
            bump._radial_potential_slope_over_radius(np.hypot(rb, zb - bump.center))
            for bump in self.bumps
        ]

    def dpsi1_dr(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        rb = np.broadcast_to(r_array, shape)
        zb = np.broadcast_to(z_array, shape)
        total = np.zeros(shape, dtype=np.float64)
        for bump, slope in zip(self.bumps, self._slope_over_radius(rb, zb)):
            total = total + slope * rb
        return total

    def dpsi1_dz(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        rb = np.broadcast_to(r_array, shape)
        zb = np.broadcast_to(z_array, shape)
        total = np.zeros(shape, dtype=np.float64)
        for bump, slope in zip(self.bumps, self._slope_over_radius(rb, zb)):
            total = total + slope * (zb - bump.center)
        return total

    def far_field(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        """The monopole approximation ``M/(3R^3)`` about the origin."""
        radius5 = np.hypot(
            np.asarray(r, dtype=np.float64), np.asarray(z, dtype=np.float64)
        )
        return self.mass / (3.0 * np.maximum(radius5, 1.0e-300) ** 3)


# --------------------------------------------------------------------------- #
# 2. Source moments computed from grid data (never from the reference)         #
# --------------------------------------------------------------------------- #


def _trapezoid(values: FloatArray, spacing: float, axis: int) -> FloatArray:
    """Composite trapezoid rule, implemented locally for NumPy-version safety."""
    moved = np.moveaxis(np.asarray(values, dtype=np.float64), axis, 0)
    if moved.shape[0] < 2:
        raise ValueError("at least two samples are required to integrate")
    weights = np.full(moved.shape[0], spacing, dtype=np.float64)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    return np.tensordot(weights, moved, axes=(0, 0))


def source_monopole_moment(grid: AxisymmetricGrid, omega1: npt.ArrayLike) -> float:
    """``M = m(∞) = (1/|S^4|) ∫ ω dV_5`` from grid data alone.

    ``dV_5 = |S^3| r^3 dr dz``, so ``M = (|S^3|/|S^4|) ∫∫ ω r^3 dr dz``.  This
    is what the monopole boundary condition consumes; it never reads the
    analytic reference, so using it keeps the solver non-circular.
    """
    values = grid.validate_field(omega1, name="omega1")
    weighted = values * grid.r[:, None] ** 3
    radial = _trapezoid(_trapezoid(weighted, grid.dz, axis=1), grid.dr, axis=0)
    return float(SPHERE3_MEASURE / SPHERE4_MEASURE * radial)


def source_l1_norm_5d(grid: AxisymmetricGrid, omega1: npt.ArrayLike) -> float:
    """``∫ |ω| dV_5`` from grid data alone."""
    values = np.abs(grid.validate_field(omega1, name="omega1"))
    weighted = values * grid.r[:, None] ** 3
    radial = _trapezoid(_trapezoid(weighted, grid.dz, axis=1), grid.dr, axis=0)
    return float(SPHERE3_MEASURE * radial)


# --------------------------------------------------------------------------- #
# 3. The a posteriori tail bound                                               #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TailBoundReport:
    """A posteriori bound on the far-field error of the monopole truncation."""

    support_radius: float
    boundary_min_radius: float
    source_l1_norm_5d: float
    bound: float
    measured_boundary_error: float

    @property
    def bound_dominates(self) -> bool:
        return self.measured_boundary_error <= self.bound

    @property
    def sharpness(self) -> float:
        """Measured error divided by the bound; ``≤ 1`` when the bound holds."""
        return self.measured_boundary_error / self.bound if self.bound > 0.0 else math.inf

    def as_dict(self) -> dict[str, object]:
        return {
            "support_radius": self.support_radius,
            "boundary_min_radius": self.boundary_min_radius,
            "source_l1_norm_5d": self.source_l1_norm_5d,
            "bound": self.bound,
            "measured_boundary_error": self.measured_boundary_error,
            "bound_dominates": self.bound_dominates,
            "sharpness": self.sharpness,
        }


def monopole_tail_bound(
    *,
    source_l1_norm_5d: float,
    support_radius: float,
    boundary_min_radius: float,
) -> float:
    r"""Rigorous bound on ``|ψ_∞ - M_5 G|`` outside the source support.

    With ``G(X) = 1/(8π^2|X|^3)`` the fundamental solution of ``-Δ`` on
    ``ℝ^5``, a source supported in ``|Y| ≤ a`` gives

    .. math::

       |\psi_\infty(X) - M_5 G(X)|
         = \Bigl|\int (G(X-Y)-G(X))\,\omega(Y)\,dY\Bigr|
         \le a\,\sup_{|Y|\le a}\|\nabla G(X-Y)\|\,\|\omega\|_{L^1(dV_5)}
         \le \frac{3a}{8\pi^2}\,
             \frac{\|\omega\|_{L^1(dV_5)}}{(|X|-a)^4},

    using ``|∇G(Z)| = 3/(8π^2|Z|^4)`` and ``|X-Y| ≥ |X|-a``.

    Evaluated at ``|X| = boundary_min_radius`` this bounds the Dirichlet data
    error of the monopole boundary condition on the whole box boundary, because
    the bound decreases in ``|X|``.  By the maximum principle for ``L_5`` the
    same number bounds the *continuum* domain-truncation error everywhere
    inside the box.  It does **not** bound the discretisation error, and it is
    not an interval-arithmetic result: the inputs are binary64 quadratures.
    """
    if source_l1_norm_5d < 0.0:
        raise ValueError("source_l1_norm_5d must be nonnegative")
    if support_radius <= 0.0:
        raise ValueError("support_radius must be positive")
    if boundary_min_radius <= support_radius:
        raise ValueError(
            "boundary_min_radius must exceed support_radius for the bound to apply"
        )
    gap = boundary_min_radius - support_radius
    return 3.0 * support_radius * source_l1_norm_5d / (8.0 * math.pi**2 * gap**4)


# --------------------------------------------------------------------------- #
# 4. The non-periodic-z solver                                                 #
# --------------------------------------------------------------------------- #


def _sine_transform_matrix(interior: int) -> FloatArray:
    """Dense DST-I matrix ``S[m-1, j-1] = sin(π j m /(N+1))``, no FFT involved.

    ``S`` is symmetric with ``S @ S = ((N+1)/2) I``, so the inverse transform is
    ``(2/(N+1)) S``.  Building it from :func:`numpy.sin` keeps the axial path
    independent of every other solver in this repository, all of which share
    ``numpy.fft`` conventions.
    """
    if interior < 1:
        raise ValueError("at least one interior axial node is required")
    indices = np.arange(1, interior + 1, dtype=np.float64)
    return np.sin(np.pi * np.outer(indices, indices) / (interior + 1.0))


def _thomas_batch(
    lower: FloatArray,
    diagonal: FloatArray,
    upper: FloatArray,
    rhs: FloatArray,
    *,
    pivot_floor: float,
) -> FloatArray:
    """Solve a batch of tridiagonal systems.  ``rhs`` has shape ``(n, batch)``.

    ``diagonal`` has shape ``(n, batch)`` because the axial eigenvalue shifts
    it; ``lower`` and ``upper`` have shape ``(n,)``.  A local implementation is
    used deliberately: no external sparse solver, and no code shared with the
    periodic solvers.
    """
    size = rhs.shape[0]
    c = np.zeros_like(rhs)
    d = np.zeros_like(rhs)
    pivot = diagonal[0].copy()
    if np.any(np.abs(pivot) < pivot_floor):
        raise ZeroDivisionError("tridiagonal pivot collapsed at the axis row")
    c[0] = upper[0] / pivot
    d[0] = rhs[0] / pivot
    smallest = float(np.min(np.abs(pivot)))
    for i in range(1, size):
        pivot = diagonal[i] - lower[i] * c[i - 1]
        smallest = min(smallest, float(np.min(np.abs(pivot))))
        if np.any(np.abs(pivot) < pivot_floor):
            raise ZeroDivisionError(f"tridiagonal pivot collapsed at radial row {i}")
        if i < size - 1:
            c[i] = upper[i] / pivot
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / pivot
    solution = np.empty_like(rhs)
    solution[-1] = d[-1]
    for i in range(size - 2, -1, -1):
        solution[i] = d[i] - c[i] * solution[i + 1]
    return solution


@dataclass(frozen=True)
class WholeSpaceEllipticSolver:
    """``-L_5 ψ = ω`` on a finite box with **non-periodic** ``z``.

    Unknowns are the radial nodes ``r_0 = 0 … r_{nr-2}`` and the *interior*
    axial nodes ``z_1 … z_{nz-2}``.  The Dirichlet traces at ``r = R_max`` and
    ``z = ±Z_max`` are moved to the right-hand side, so the interior axial
    operator is the plain second difference and is diagonalised exactly by the
    dense DST-I of :func:`_sine_transform_matrix`.  Every axial eigenvalue
    ``λ_m = -4 sin²(mπ/(2(N+1)))/Δz²`` is strictly negative: unlike the
    periodic case there is **no zero axial mode**, which is precisely the
    structural difference Gate 4 exists to expose.

    Item 10 of the gate — a reusable API for the nonlinear evolution — is the
    reason this is a stateful object: the transform matrix and the radial
    stencil are built once and :meth:`solve` may then be called every stage of
    every time step.
    """

    grid: AxisymmetricGrid
    boundary_mode: str
    sine: FloatArray
    axial_eigenvalues: FloatArray
    radial_lower: FloatArray
    radial_upper: FloatArray
    radial_diagonal_base: FloatArray
    outer_coupling: float
    axis_coefficient: float
    pivot_floor: float

    @classmethod
    def build(
        cls,
        grid: AxisymmetricGrid,
        *,
        boundary_mode: str = "monopole",
        axis_coefficient: float = 8.0,
        radial_drift_coefficient: float = 3.0,
        pivot_floor: float = 1.0e-30,
    ) -> "WholeSpaceEllipticSolver":
        """Assemble the solver.

        ``axis_coefficient`` and ``radial_drift_coefficient`` are exposed only
        so the fault-injection tests can corrupt the audited values ``8`` and
        ``3``; production callers must leave them alone.
        """
        if grid.periodic_z:
            raise ValueError("Gate 4 requires a non-periodic z grid")
        if boundary_mode not in {"zero", "monopole"}:
            raise ValueError("boundary_mode must be 'zero' or 'monopole'")
        interior = grid.nz - 2
        if interior < 1:
            raise ValueError("at least three axial nodes are required")

        dr = grid.dr
        dz = grid.dz
        radial_unknowns = grid.nr - 1  # r = R_max is Dirichlet data

        modes = np.arange(1, interior + 1, dtype=np.float64)
        eigenvalues = -4.0 / dz**2 * np.sin(np.pi * modes / (2.0 * (interior + 1.0))) ** 2

        lower = np.zeros(radial_unknowns, dtype=np.float64)
        upper = np.zeros(radial_unknowns, dtype=np.float64)
        diagonal = np.zeros(radial_unknowns, dtype=np.float64)

        # Axis row: (L_r ψ)_0 = axis_coefficient (ψ_1 - ψ_0)/dr².
        diagonal[0] = axis_coefficient / dr**2
        upper[0] = -axis_coefficient / dr**2
        for i in range(1, radial_unknowns):
            drift = radial_drift_coefficient / (2.0 * grid.r[i] * dr)
            lower[i] = -(1.0 / dr**2 - drift)
            diagonal[i] = 2.0 / dr**2
            if i < radial_unknowns - 1:
                upper[i] = -(1.0 / dr**2 + drift)
        # Coupling of the last unknown to the Dirichlet node at r = R_max.
        outer_drift = radial_drift_coefficient / (
            2.0 * grid.r[radial_unknowns - 1] * dr
        )
        outer_coupling = 1.0 / dr**2 + outer_drift

        return cls(
            grid=grid,
            boundary_mode=boundary_mode,
            sine=_sine_transform_matrix(interior),
            axial_eigenvalues=eigenvalues,
            radial_lower=lower,
            radial_upper=upper,
            radial_diagonal_base=diagonal,
            outer_coupling=float(outer_coupling),
            axis_coefficient=float(axis_coefficient),
            pivot_floor=float(pivot_floor),
        )

    # -- boundary data ----------------------------------------------------- #

    def boundary_trace(self, omega1: npt.ArrayLike) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Dirichlet data ``(outer_r, low_z, high_z)`` for the chosen mode.

        For ``boundary_mode='monopole'`` the trace is ``M/(3R^3)`` with ``M``
        obtained from the *discrete source* by :func:`source_monopole_moment`.
        The analytic reference is never consulted.
        """
        grid = self.grid
        outer = np.zeros(grid.nz, dtype=np.float64)
        low = np.zeros(grid.nr, dtype=np.float64)
        high = np.zeros(grid.nr, dtype=np.float64)
        if self.boundary_mode == "zero":
            return outer, low, high
        mass = source_monopole_moment(grid, omega1)
        r_max = float(grid.r[-1])
        z_low = float(grid.z[0])
        z_high = float(grid.z[-1])
        outer = mass / (3.0 * np.hypot(r_max, grid.z) ** 3)
        low = mass / (3.0 * np.hypot(grid.r, z_low) ** 3)
        high = mass / (3.0 * np.hypot(grid.r, z_high) ** 3)
        return outer, low, high

    # -- the solve --------------------------------------------------------- #

    def solve(self, omega1: npt.ArrayLike) -> FloatArray:
        """Return ``ψ_1`` on the full grid, including the Dirichlet nodes."""
        grid = self.grid
        source = grid.validate_field(omega1, name="omega1")
        outer, low, high = self.boundary_trace(source)

        interior = grid.nz - 2
        radial_unknowns = grid.nr - 1

        # -L5 ψ = ω, with known traces moved to the right-hand side.
        rhs = source[:radial_unknowns, 1:-1].copy()
        rhs[:, 0] += low[:radial_unknowns] / grid.dz**2
        rhs[:, -1] += high[:radial_unknowns] / grid.dz**2
        rhs[-1, :] += self.outer_coupling * outer[1:-1]

        # DST-I along z, then one tridiagonal solve per axial eigenvalue.
        transformed = rhs @ self.sine.T
        diagonal = (
            self.radial_diagonal_base[:, None]
            - self.axial_eigenvalues[None, :]
        )
        solution_hat = _thomas_batch(
            self.radial_lower,
            diagonal,
            self.radial_upper,
            transformed,
            pivot_floor=self.pivot_floor,
        )
        interior_solution = (solution_hat @ self.sine.T) * (2.0 / (interior + 1.0))

        psi = np.zeros(grid.shape, dtype=np.float64)
        psi[:radial_unknowns, 1:-1] = interior_solution
        psi[-1, :] = outer
        psi[:, 0] = low
        psi[:, -1] = high
        psi[-1, 0] = 0.5 * (outer[0] + low[-1])
        psi[-1, -1] = 0.5 * (outer[-1] + high[-1])
        return psi

    # -- diagnostics ------------------------------------------------------- #

    # -- item 10: the interface the nonlinear evolution will consume -------- #

    def recover_meridional_velocity(
        self, psi1: npt.ArrayLike
    ) -> tuple[FloatArray, FloatArray]:
        """``(u^r, u^z)`` from ``ψ_1`` through the audited E-14 recovery.

        Delegates to :func:`ns_certificate_lab.operators.recover_velocity` so
        that the whole-space path uses exactly the same recovery as every
        periodic result already in the repository; only the elliptic solve
        differs.
        """
        from .operators import recover_velocity

        return recover_velocity(self.grid, psi1)

    @staticmethod
    def coupling_contract() -> dict[str, object]:
        """Machine-readable statement of what the nonlinear solver may assume.

        Gate 4 is a **linear** gate.  Passing it licenses coupling this solver
        into a time-stepper; it does not license any claim about the resulting
        evolution, and every entry under ``not_established`` remains open.
        """
        return {
            "solves": "-L5 psi1 = omega1",
            "axial_boundary": "non-periodic Dirichlet at z = +/- Z_max",
            "radial_boundary": "Dirichlet at r = R_max",
            "outer_trace": "zero, or the five-dimensional monopole M/(3R^3) "
            "with M taken from the discrete source by quadrature",
            "reusable": "build() once, solve() every stage of every time step",
            "state_mutated_by_solve": False,
            "velocity_recovery": "operators.recover_velocity (E-14)",
            "established": [
                "second-order convergence against an exact closed-form free-space "
                "reference at fixed domain",
                "independent R_max and Z_max enlargement with the truncation "
                "component isolated by common-interior differencing",
                "no zero axial mode exists, unlike the periodic operator",
                "an a posteriori bound on the continuum truncation error",
                "independent Cartesian divergence and curl audit of the recovery",
            ],
            "not_established": [
                "any statement about the nonlinear evolution",
                "an interval-arithmetic enclosure of any quantity",
                "a discrete maximum principle (radial row i=1 is not an M-matrix)",
                "resolution of a Hou-type front on a uniform whole-space grid",
            ],
        }

    def residual(self, psi1: npt.ArrayLike, omega1: npt.ArrayLike) -> FloatArray:
        """``-L_5 ψ - ω`` on interior nodes, by a second physical-space path."""
        grid = self.grid
        psi = grid.validate_field(psi1, name="psi1")
        source = grid.validate_field(omega1, name="omega1")
        dr, dz = grid.dr, grid.dz
        laplacian = np.zeros_like(psi)
        laplacian[1:-1, :] = (psi[2:, :] - 2.0 * psi[1:-1, :] + psi[:-2, :]) / dr**2
        laplacian[1:-1, :] += (
            3.0
            / grid.r[1:-1, None]
            * (psi[2:, :] - psi[:-2, :])
            / (2.0 * dr)
        )
        laplacian[0, :] = self.axis_coefficient * (psi[1, :] - psi[0, :]) / dr**2
        axial = np.zeros_like(psi)
        axial[:, 1:-1] = (psi[:, 2:] - 2.0 * psi[:, 1:-1] + psi[:, :-2]) / dz**2
        out = -(laplacian + axial) - source
        out[-1, :] = 0.0
        out[:, 0] = 0.0
        out[:, -1] = 0.0
        return out


# --------------------------------------------------------------------------- #
# 5. Domain convergence and error separation                                   #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DomainConvergenceRow:
    """One ``(R_max, Z_max, n_r, n_z)`` entry of the Gate-4 convergence table."""

    r_max: float
    z_max: float
    nr: int
    nz: int
    dr: float
    dz: float
    boundary_mode: str
    interior_relative_l2: float
    interior_max_abs: float
    boundary_data_error_max: float
    tail_bound: float
    residual_max: float

    def as_dict(self) -> dict[str, object]:
        return {
            "r_max": self.r_max,
            "z_max": self.z_max,
            "nr": self.nr,
            "nz": self.nz,
            "dr": self.dr,
            "dz": self.dz,
            "boundary_mode": self.boundary_mode,
            "interior_relative_l2": self.interior_relative_l2,
            "interior_max_abs": self.interior_max_abs,
            "boundary_data_error_max": self.boundary_data_error_max,
            "tail_bound": self.tail_bound,
            "residual_max": self.residual_max,
        }


def _interior_mask(grid: AxisymmetricGrid, *, radius: float) -> tuple[FloatArray, FloatArray]:
    r_mesh, z_mesh = grid.mesh()
    return np.hypot(r_mesh, z_mesh) <= radius, r_mesh


def solve_and_compare(
    field: ManufacturedFreeSpaceField,
    *,
    r_max: float,
    z_max: float,
    nr: int,
    nz: int,
    boundary_mode: str = "monopole",
    interior_radius: float,
    axis_coefficient: float = 8.0,
    radial_drift_coefficient: float = 3.0,
) -> DomainConvergenceRow:
    """Run one box and measure the interior error against the exact reference."""
    grid = AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=r_max,
        z_min=-z_max,
        z_max=z_max,
        periodic_z=False,
    )
    r_mesh, z_mesh = grid.mesh()
    omega = field.omega1(r_mesh, z_mesh)
    exact = field.psi1(r_mesh, z_mesh)
    solver = WholeSpaceEllipticSolver.build(
        grid,
        boundary_mode=boundary_mode,
        axis_coefficient=axis_coefficient,
        radial_drift_coefficient=radial_drift_coefficient,
    )
    psi = solver.solve(omega)
    error = psi - exact

    mask, _ = _interior_mask(grid, radius=interior_radius)
    if not np.any(mask):
        raise ValueError("interior_radius selects no grid point")
    weights = grid.r[:, None] ** 3
    numerator = float(np.sum(weights * error**2 * mask))
    denominator = float(np.sum(weights * exact**2 * mask))
    relative = math.sqrt(numerator / denominator) if denominator > 0.0 else math.inf

    outer, low, high = solver.boundary_trace(omega)
    boundary_error = max(
        float(np.max(np.abs(outer - exact[-1, :]))),
        float(np.max(np.abs(low - exact[:, 0]))),
        float(np.max(np.abs(high - exact[:, -1]))),
    )
    boundary_min_radius = min(r_max, z_max)
    bound = (
        monopole_tail_bound(
            source_l1_norm_5d=source_l1_norm_5d(grid, omega),
            support_radius=field.support_radius,
            boundary_min_radius=boundary_min_radius,
        )
        if boundary_min_radius > field.support_radius
        else math.inf
    )
    return DomainConvergenceRow(
        r_max=float(r_max),
        z_max=float(z_max),
        nr=int(nr),
        nz=int(nz),
        dr=grid.dr,
        dz=grid.dz,
        boundary_mode=boundary_mode,
        interior_relative_l2=relative,
        interior_max_abs=float(np.max(np.abs(error * mask))),
        boundary_data_error_max=boundary_error,
        tail_bound=bound,
        residual_max=float(np.max(np.abs(solver.residual(psi, omega)))),
    )


# --------------------------------------------------------------------------- #
# 6/7. Cartesian recovery and independent divergence/curl audit                #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CartesianGateReport:
    """Independent Cartesian check of the recovered whole-space velocity."""

    points: int
    velocity_max: float
    divergence_max: float
    divergence_relative: float
    curl_theta_error_max: float
    curl_theta_relative: float
    curl_axial_max: float

    def as_dict(self) -> dict[str, object]:
        return {
            "points": self.points,
            "velocity_max": self.velocity_max,
            "divergence_max": self.divergence_max,
            "divergence_relative": self.divergence_relative,
            "curl_theta_error_max": self.curl_theta_error_max,
            "curl_theta_relative": self.curl_theta_relative,
            "curl_axial_max": self.curl_axial_max,
        }


def audit_cartesian_recovery(
    field: ManufacturedFreeSpaceField,
    *,
    half_width: float,
    points: int,
    interior_trim: int = 2,
) -> CartesianGateReport:
    r"""Recover ``u`` on a Cartesian box and check ``∇·u`` and ``∇×u``.

    With no swirl the regularized recovery (E-14/E-18) is

    .. math::

       u_x = -x\,\psi_{1,z},\quad u_y = -y\,\psi_{1,z},\quad
       u_z = 2\psi_1 + r\,\psi_{1,r},

    which never divides by ``r``.  Two identities are then checked with
    Cartesian stencils that do not call any cylindrical operator:

    * ``∇·u = 0`` (E-15), and
    * ``∇×u = ω_1 \,(-y, x, 0)``, because the only nonzero vorticity component
      is ``ω^θ = r\omega_1`` and ``\hat e_θ = (-y,x,0)/r``.

    The second identity is the one that matters for Gate 4: it ties the
    elliptic solve ``-L_5ψ_1 = ω_1`` to the *physical* curl through a path that
    shares no code with the solver.
    """
    from .cartesian_validation import (
        RegularizedAxisymmetricProfile,
        UniformCartesianGrid,
        cartesian_curl,
        cartesian_divergence,
        reconstruct_axisymmetric_regularized,
    )

    if points < 8:
        raise ValueError("at least eight points per axis are required")
    if interior_trim < 1:
        raise ValueError("interior_trim must be at least one to drop one-sided rows")

    axis = np.linspace(-half_width, half_width, points, dtype=np.float64)
    grid = UniformCartesianGrid(x=axis, y=axis, z=axis)
    profile = RegularizedAxisymmetricProfile(
        u1=lambda r, z: np.zeros_like(r),
        psi1=field.psi1,
        dpsi1_dr=field.dpsi1_dr,
        dpsi1_dz=field.dpsi1_dz,
    )
    velocity = reconstruct_axisymmetric_regularized(grid, profile)
    divergence = cartesian_divergence(grid, velocity)
    curl = cartesian_curl(grid, velocity)

    x, y, z = grid.mesh()
    radius = np.hypot(x, y)
    omega1 = field.omega1(radius, z)
    expected = np.stack((-y * omega1, x * omega1, np.zeros_like(omega1)), axis=0)

    trim = slice(interior_trim, -interior_trim)
    core = (trim, trim, trim)
    velocity_scale = float(np.max(np.abs(velocity[(slice(None), *core)])))
    curl_scale = float(np.max(np.abs(expected[(slice(None), *core)])))
    divergence_max = float(np.max(np.abs(divergence[core])))
    curl_error = np.abs(curl - expected)
    curl_theta_error = float(np.max(curl_error[0][core]))
    curl_theta_error = max(curl_theta_error, float(np.max(curl_error[1][core])))
    return CartesianGateReport(
        points=points,
        velocity_max=velocity_scale,
        divergence_max=divergence_max,
        divergence_relative=divergence_max / max(velocity_scale, 1.0e-300) * float(
            axis[1] - axis[0]
        ),
        curl_theta_error_max=curl_theta_error,
        curl_theta_relative=curl_theta_error / max(curl_scale, 1.0e-300),
        curl_axial_max=float(np.max(np.abs(curl[2][core]))),
    )


# --------------------------------------------------------------------------- #
# 8. Domain truncation isolated from discretisation error                      #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TruncationRow:
    """Common-interior difference between one box and the largest box.

    Enlarging the box at **fixed** ``Δr`` and ``Δz`` leaves the discretisation
    error essentially unchanged, so the difference between two solutions on
    their common interior isolates the *domain-truncation* component.  This is
    the only way to see the truncation rate at all: a direct comparison with
    the exact reference saturates at the discretisation floor, which
    ``outputs/whole_space_gate4_v1`` records explicitly.
    """

    r_max: float
    z_max: float
    boundary_mode: str
    common_interior_max_diff: float
    boundary_data_error_max: float
    tail_bound: float

    def as_dict(self) -> dict[str, object]:
        return {
            "r_max": self.r_max,
            "z_max": self.z_max,
            "boundary_mode": self.boundary_mode,
            "common_interior_max_diff": self.common_interior_max_diff,
            "boundary_data_error_max": self.boundary_data_error_max,
            "tail_bound": self.tail_bound,
        }


def _solve_on_box(
    field: ManufacturedFreeSpaceField,
    *,
    r_max: float,
    z_max: float,
    spacing: float,
    boundary_mode: str,
) -> tuple[AxisymmetricGrid, FloatArray, FloatArray, float]:
    """Solve on ``[0,r_max] x [-z_max,z_max]`` at the requested uniform spacing."""
    nr = int(round(r_max / spacing)) + 1
    nz = 2 * int(round(z_max / spacing)) + 1
    grid = AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=r_max, z_min=-z_max, z_max=z_max, periodic_z=False
    )
    r_mesh, z_mesh = grid.mesh()
    omega = field.omega1(r_mesh, z_mesh)
    solver = WholeSpaceEllipticSolver.build(grid, boundary_mode=boundary_mode)
    psi = solver.solve(omega)
    exact = field.psi1(r_mesh, z_mesh)
    outer, low, high = solver.boundary_trace(omega)
    boundary_error = max(
        float(np.max(np.abs(outer - exact[-1, :]))),
        float(np.max(np.abs(low - exact[:, 0]))),
        float(np.max(np.abs(high - exact[:, -1]))),
    )
    return grid, psi, exact, boundary_error


def domain_truncation_series(
    field: ManufacturedFreeSpaceField,
    *,
    sizes: Sequence[float],
    fixed_extent: float,
    spacing: float,
    enlarge: str,
    boundary_mode: str,
    interior_radius: float,
    reference_size: float | None = None,
) -> list[TruncationRow]:
    """Enlarge ``R_max`` or ``Z_max`` at fixed spacing and isolate truncation.

    ``enlarge`` is ``'radial'`` or ``'axial'``; the other extent is held at
    ``fixed_extent``.  The reference solution is taken on a *separate* box of
    extent ``reference_size`` (default: four times the largest entry of
    ``sizes``), never on the largest member of the series itself — otherwise
    the last row would report an identically zero difference and the observed
    rate would be biased by the reference's own truncation error.
    """
    if enlarge not in {"radial", "axial"}:
        raise ValueError("enlarge must be 'radial' or 'axial'")
    ordered = [float(value) for value in sizes]
    if len(ordered) < 2 or any(b <= a for a, b in zip(ordered, ordered[1:])):
        raise ValueError("sizes must be strictly increasing with at least two entries")
    if interior_radius >= min(min(ordered), fixed_extent):
        raise ValueError("interior_radius must lie inside the smallest box")

    reference_extent = float(reference_size) if reference_size is not None else 4.0 * ordered[-1]
    if reference_extent <= ordered[-1]:
        raise ValueError("reference_size must exceed the largest entry of sizes")

    solutions = []
    for size in ordered:
        r_max = size if enlarge == "radial" else fixed_extent
        z_max = fixed_extent if enlarge == "radial" else size
        grid, psi, _exact, boundary_error = _solve_on_box(
            field,
            r_max=r_max,
            z_max=z_max,
            spacing=spacing,
            boundary_mode=boundary_mode,
        )
        solutions.append((grid, psi, r_max, z_max, boundary_error))

    reference_grid, reference_psi, _, _ = _solve_on_box(
        field,
        r_max=reference_extent if enlarge == "radial" else fixed_extent,
        z_max=fixed_extent if enlarge == "radial" else reference_extent,
        spacing=spacing,
        boundary_mode=boundary_mode,
    )
    rows: list[TruncationRow] = []
    for grid, psi, r_max, z_max, boundary_error in solutions:
        radial_index = np.arange(grid.nr)
        axial_offset = int(round((grid.z[0] - reference_grid.z[0]) / grid.dz))
        axial_index = np.arange(grid.nz) + axial_offset
        keep = (axial_index >= 0) & (axial_index < reference_grid.nz)
        selected = reference_psi[np.ix_(radial_index, axial_index[keep])]
        r_mesh, z_mesh = grid.mesh()
        mask = (np.hypot(r_mesh, z_mesh) <= interior_radius)[:, keep]
        rows.append(
            TruncationRow(
                r_max=r_max,
                z_max=z_max,
                boundary_mode=boundary_mode,
                common_interior_max_diff=float(
                    np.max(np.abs((psi[:, keep] - selected) * mask))
                ),
                boundary_data_error_max=boundary_error,
                tail_bound=(
                    monopole_tail_bound(
                        source_l1_norm_5d=field.l1_norm_5d,
                        support_radius=field.support_radius,
                        boundary_min_radius=min(r_max, z_max),
                    )
                    if min(r_max, z_max) > field.support_radius
                    else math.inf
                ),
            )
        )
    return rows


# --------------------------------------------------------------------------- #
# 9. Period-image error separated from the radial-wall error                   #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SeparationReport:
    """Period-image error and radial-wall error, separated at the continuum level.

    Measuring the two components by *differencing solver outputs* does not work:
    at every spacing this repository can afford, the discretisation error of all
    four treatments (``≈10^{-3}`` at ``h=1/8``) is orders of magnitude larger
    than the period-image contamination it is supposed to expose.  That failure
    is recorded in ``outputs/whole_space_gate4_v1`` rather than hidden.

    The separation below is therefore done exactly, using the closed-form
    free-space reference, and the measured numbers are reported alongside only
    as context.

    * ``period_image_exact`` -- ``max |Σ_{n≠0} ψ_∞(r, z+nL_z)|`` over the
      interior region.  This is precisely what a periodic ``z`` treatment adds
      to the free-space field and what a non-periodic treatment does not: for
      Gate 4 this component is **exactly zero**, not merely small.
    * ``radial_wall_bound`` -- ``max_z |ψ_periodised(R_max, z)|``.  The finite
      cylinder pins ``ψ=0`` at ``r=R_max``, so its difference from the
      periodised free field is ``L_5``-harmonic with exactly that boundary
      trace, and the maximum principle bounds it by this number everywhere
      inside.
    * ``nonperiodic_truncation_bound`` -- the a posteriori monopole tail bound
      of :func:`monopole_tail_bound`, which plays the same role for Gate 4.

    The three bounds are continuum statements about the boundary treatment.
    They are computed in binary64 from an exact formula; they are not interval
    arithmetic.
    """

    interior_radius: float
    axial_half_period: float
    radial_extent: float
    spacing: float
    period_image_exact: float
    period_image_tail_bound: float
    radial_wall_bound: float
    nonperiodic_truncation_bound: float
    measured_finite_cylinder: float
    measured_free_radial_periodic_z: float
    measured_nonperiodic_zero: float
    measured_nonperiodic_monopole: float

    @property
    def period_image_removed(self) -> bool:
        """Gate 4 removes the period-image component exactly, not approximately."""
        return True

    @property
    def image_to_truncation_ratio(self) -> float:
        denominator = max(self.nonperiodic_truncation_bound, 1.0e-300)
        return self.period_image_exact / denominator

    def as_dict(self) -> dict[str, object]:
        return {
            "interior_radius": self.interior_radius,
            "axial_half_period": self.axial_half_period,
            "radial_extent": self.radial_extent,
            "spacing": self.spacing,
            "period_image_exact": self.period_image_exact,
            "period_image_tail_bound": self.period_image_tail_bound,
            "radial_wall_bound": self.radial_wall_bound,
            "nonperiodic_truncation_bound": self.nonperiodic_truncation_bound,
            "measured_finite_cylinder": self.measured_finite_cylinder,
            "measured_free_radial_periodic_z": self.measured_free_radial_periodic_z,
            "measured_nonperiodic_zero": self.measured_nonperiodic_zero,
            "measured_nonperiodic_monopole": self.measured_nonperiodic_monopole,
            "period_image_removed": self.period_image_removed,
            "image_to_truncation_ratio": self.image_to_truncation_ratio,
        }


def _periodic_image_sum(
    field: ManufacturedFreeSpaceField,
    r: FloatArray,
    z: FloatArray,
    *,
    period: float,
    images: int,
) -> FloatArray:
    """``Σ_{n≠0} ψ_∞(r, z + nL_z)`` truncated at ``|n| ≤ images``."""
    total = np.zeros(np.broadcast(r, z).shape, dtype=np.float64)
    for n in range(1, images + 1):
        total = total + field.psi1(r, z + n * period)
        total = total + field.psi1(r, z - n * period)
    return total


def _image_tail_bound(
    field: ManufacturedFreeSpaceField,
    *,
    period: float,
    images: int,
    interior_radius: float,
) -> float:
    """Bound the images beyond ``|n| > images`` by ``Σ |M|/(3 d_n^3)``."""
    reach = interior_radius + field.support_radius
    gap = (images + 1) * period - reach
    if gap <= 0.0:
        return math.inf
    # Σ_{n>images} 2|M|/(3 (nL - reach)^3) ≤ 2|M|/(3) * (1/gap^3 + ∫ …)
    tail = 1.0 / gap**3 + 1.0 / (2.0 * period * gap**2)
    return 2.0 * abs(field.mass) / 3.0 * tail


def separate_period_and_wall_errors(
    field: ManufacturedFreeSpaceField,
    *,
    radial_extent: float,
    axial_half_period: float,
    spacing: float,
    interior_radius: float,
    images: int = 400,
) -> SeparationReport:
    """Separate the two finite-domain artefacts exactly; measure all four solvers."""
    from .finite_cylinder_poisson import solve_finite_cylinder_poisson
    from .free_space_poisson import solve_l5_free_radial_periodized_z

    if interior_radius >= min(radial_extent, axial_half_period):
        raise ValueError("interior_radius must lie inside every box")
    period = 2.0 * axial_half_period

    nr = int(round(radial_extent / spacing)) + 1
    nz_periodic = 2 * int(round(axial_half_period / spacing))
    periodic_grid = AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz_periodic,
        r_max=radial_extent,
        z_min=-axial_half_period,
        z_max=axial_half_period,
        periodic_z=True,
    )
    r_mesh, z_mesh = periodic_grid.mesh()
    omega_periodic = field.omega1(r_mesh, z_mesh)
    exact_periodic = field.psi1(r_mesh, z_mesh)
    mask_periodic = np.hypot(r_mesh, z_mesh) <= interior_radius

    images_interior = _periodic_image_sum(
        field, r_mesh, z_mesh, period=period, images=images
    )
    period_image_exact = float(np.max(np.abs(images_interior * mask_periodic)))
    wall_radius = np.full(periodic_grid.nz, radial_extent)
    periodised_wall = field.psi1(wall_radius, periodic_grid.z) + _periodic_image_sum(
        field, wall_radius, periodic_grid.z, period=period, images=images
    )

    cylinder = solve_finite_cylinder_poisson(
        omega_periodic,
        periodic_grid,
        outer_boundary=0.0,
        condition_mode_indices=(),
    ).psi
    free_radial, _ = solve_l5_free_radial_periodized_z(
        omega_periodic, periodic_grid.r, periodic_grid.dz
    )

    measured: dict[str, float] = {
        "measured_finite_cylinder": float(
            np.max(np.abs((cylinder - exact_periodic) * mask_periodic))
        ),
        "measured_free_radial_periodic_z": float(
            np.max(np.abs((free_radial - exact_periodic) * mask_periodic))
        ),
    }
    for mode, key in (
        ("zero", "measured_nonperiodic_zero"),
        ("monopole", "measured_nonperiodic_monopole"),
    ):
        grid, psi, exact, _ = _solve_on_box(
            field,
            r_max=radial_extent,
            z_max=axial_half_period,
            spacing=spacing,
            boundary_mode=mode,
        )
        r_np, z_np = grid.mesh()
        mask = np.hypot(r_np, z_np) <= interior_radius
        measured[key] = float(np.max(np.abs((psi - exact) * mask)))

    boundary_min_radius = min(radial_extent, axial_half_period)
    truncation_bound = (
        monopole_tail_bound(
            source_l1_norm_5d=field.l1_norm_5d,
            support_radius=field.support_radius,
            boundary_min_radius=boundary_min_radius,
        )
        if boundary_min_radius > field.support_radius
        else math.inf
    )
    return SeparationReport(
        interior_radius=interior_radius,
        axial_half_period=axial_half_period,
        radial_extent=radial_extent,
        spacing=spacing,
        period_image_exact=period_image_exact,
        period_image_tail_bound=_image_tail_bound(
            field, period=period, images=images, interior_radius=interior_radius
        ),
        radial_wall_bound=float(np.max(np.abs(periodised_wall))),
        nonperiodic_truncation_bound=truncation_bound,
        **measured,
    )


# --------------------------------------------------------------------------- #
# 11. Zero-frequency stress test: the periodic R^-2 tail versus R^-3           #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ZeroFrequencyRow:
    """Far-field magnitude of the periodic zero mode against the free field."""

    radial_extent: float
    axial_period: float
    periodic_zero_mode_far_field: float
    free_space_far_field: float

    @property
    def ratio(self) -> float:
        return self.periodic_zero_mode_far_field / max(self.free_space_far_field, 1.0e-300)

    def as_dict(self) -> dict[str, object]:
        return {
            "radial_extent": self.radial_extent,
            "axial_period": self.axial_period,
            "periodic_zero_mode_far_field": self.periodic_zero_mode_far_field,
            "free_space_far_field": self.free_space_far_field,
            "ratio": self.ratio,
        }


def zero_frequency_stress(
    field: ManufacturedFreeSpaceField,
    *,
    radial_extents: Sequence[float],
    axial_period: float,
) -> list[ZeroFrequencyRow]:
    r"""Compare the periodic zero axial mode with the free-space monopole.

    For a periodic ``z`` of period ``L`` the ``k=0`` axial mode obeys the
    four-dimensional radial Poisson equation whose whole-space kernel is
    ``G_0(r,ρ) = 1/(2\max(r,ρ)^2)``.  Outside the source therefore

    .. math::

       \psi_{k=0}(r)=\frac{1}{2r^2}\int \rho^3 f_0(\rho)\,d\rho
                    =\frac{M_5}{4\pi^2 L\,r^{2}} ,

    an **algebraic ``r^{-2}`` tail**, whereas the true free-space field decays
    like ``M_5/(8\pi^2 R^3)``.  Their ratio at ``|X| = R`` is exactly ``2R/L``:
    the periodic treatment overstates the far field by a factor that *grows
    linearly with the radial extent*.  This is the quantitative form of
    ``docs/low_frequency_wall_obstruction.md`` and the reason a periodic ``z``
    box cannot be enlarged into a whole-space computation.

    Gate 4 removes the mode entirely: :func:`axial_spectrum` shows every
    non-periodic Dirichlet eigenvalue is strictly negative.
    """
    if axial_period <= 0.0:
        raise ValueError("axial_period must be positive")
    moment5 = field.monopole_moment_5d
    rows: list[ZeroFrequencyRow] = []
    for extent in radial_extents:
        if extent <= field.support_radius:
            raise ValueError("radial_extents must exceed the source support radius")
        rows.append(
            ZeroFrequencyRow(
                radial_extent=float(extent),
                axial_period=float(axial_period),
                periodic_zero_mode_far_field=abs(
                    moment5 / (4.0 * math.pi**2 * axial_period * extent**2)
                ),
                free_space_far_field=abs(moment5 / (8.0 * math.pi**2 * extent**3)),
            )
        )
    return rows


# --------------------------------------------------------------------------- #
# 10. Low- and zero-frequency stress test                                      #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class AxialSpectrumReport:
    """The axial eigenvalue spectrum of the non-periodic Dirichlet operator.

    The periodic solvers carry an exact zero axial mode whose radial Green
    function has the slow algebraic tail ``R^{-2}``
    (``docs/low_frequency_wall_obstruction.md``).  With non-periodic ``z`` and
    Dirichlet ends every eigenvalue satisfies ``lambda_m < 0`` strictly, so
    **there is no zero mode at all**; the slowest decaying object becomes the
    five-dimensional monopole ``R^{-3}``.  ``smallest_magnitude`` records how
    close the spectrum gets to zero as ``Z_max`` grows, which is the honest
    residual of the old obstruction.
    """

    z_max: float
    axial_points: int
    smallest_magnitude: float
    largest_magnitude: float
    zero_mode_present: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "z_max": self.z_max,
            "axial_points": self.axial_points,
            "smallest_magnitude": self.smallest_magnitude,
            "largest_magnitude": self.largest_magnitude,
            "zero_mode_present": self.zero_mode_present,
        }


def axial_spectrum(grid: AxisymmetricGrid) -> AxialSpectrumReport:
    """Report the axial spectrum used by :class:`WholeSpaceEllipticSolver`."""
    if grid.periodic_z:
        raise ValueError("axial_spectrum describes the non-periodic operator")
    interior = grid.nz - 2
    modes = np.arange(1, interior + 1, dtype=np.float64)
    eigenvalues = (
        -4.0 / grid.dz**2 * np.sin(np.pi * modes / (2.0 * (interior + 1.0))) ** 2
    )
    magnitudes = np.abs(eigenvalues)
    return AxialSpectrumReport(
        z_max=float(grid.z[-1]),
        axial_points=grid.nz,
        smallest_magnitude=float(np.min(magnitudes)),
        largest_magnitude=float(np.max(magnitudes)),
        zero_mode_present=bool(np.any(magnitudes == 0.0)),
    )
