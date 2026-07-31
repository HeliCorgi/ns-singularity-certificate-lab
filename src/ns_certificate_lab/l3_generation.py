r"""The initial generation rate of the critical norm, and the pure-swirl no-go.

Gate 7 found that the critical norm ``F(t) = \int_{\mathbb R^3}|u|^3\,dx``
decayed in all eighteen runs.  This module explains that outcome rather than
merely recording it, and turns the explanation into a search objective.

The differentiation identity
----------------------------
For a smooth, rapidly decaying, divergence-free solution of

.. math::  u_t + (u\cdot\nabla)u = -\nabla p + \nu\Delta u ,\qquad \nabla\cdot u = 0

one has

.. math::

   F'(t) = 3\int p\,\nabla\!\cdot\!(|u|u)\,dx
     - 3\nu\int\Bigl(|u|\,|\nabla u|^2
       + |u|^{-1}\sum_k (u\cdot\partial_k u)^2\Bigr)dx .

The derivation is in :doc:`../../docs/research_notes/l3_generation_rate`.  Three
features of it drive everything here.

1. **The transport term vanishes identically.**  For any divergence-free field,
   ``|u|\,u\cdot(u\cdot\nabla)u = \tfrac13 u\cdot\nabla(|u|^3)``, whose integral
   is ``-\int(\nabla\cdot u)|u|^3 = 0``.  Advection therefore contributes
   *nothing* to the initial rate.  The only way to make ``F'(0)`` positive is
   through the **pressure**.
2. **The viscous term is a sum of two nonnegative integrands**, so it always
   pushes ``F'`` down.  Using ``u\cdot\partial_k u = |u|\,\partial_k|u|`` it
   collapses to ``-3\nu\int|u|\bigl(|\nabla u|^2 + |\nabla|u||^2\bigr)dx``.
3. **A pure swirl kills the pressure term.**  If ``u_0 = u^\theta(r,z)e_\theta``
   then ``|u_0|`` depends only on ``(r,z)`` and ``|u_0|u_0 = g(r,z)e_\theta``,
   whose cylindrical divergence ``r^{-1}\partial_\theta g`` is zero.  Hence
   ``F'(0)\le 0``, strictly unless ``u_0\equiv 0``.

That last point is the **pure-swirl no-go**: no purely azimuthal initial datum
can increase the global critical norm at ``t = 0``.  Every family the repository
has run so far — S, A, H — is a pure swirl.  What the theorem forces is only
the *initial* sign ``F'(0) \le 0``: vortex stretching generates meridional flow
immediately, after which the pressure channel switches on with indefinite sign
and the theorem says nothing.  The sustained decay Gate 7 observed over
``0 < \tau \le 1`` is a numerical observation — consistent with the theorem,
not implied by it.  See :func:`pure_swirl_no_go`.

Numerical form
--------------
In cylindrical coordinates with ``u = (u^r, u^\theta, u^z)(r,z)``:

* ``\nabla\cdot(|u|u) = u\cdot\nabla|u| = u^r\partial_r q + u^z\partial_z q``
  with ``q = |u|``;
* ``\sum_k(u\cdot\partial_k u)^2 = |u|^2|\nabla q|^2``, so the second viscous
  integrand is ``|u|\,|\nabla q|^2``;
* the pressure term is integrated **back** by parts to
  ``-3\int q\,(u^r\partial_r p + u^z\partial_z p)\,dx``, which removes the only
  division by ``q`` from the dominant contribution.

The remaining ``q`` in a denominator, ``|\nabla q|^2 |u| = (G_r^2+G_z^2)/q``
with ``G = \tfrac12\nabla(|u|^2)`` smooth, is bounded by ``|u||\nabla u|^2`` and
is defined to be zero on the zero set of ``u``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import ClassVar

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .operators import derivative_r, derivative_z

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "AxisymmetricPressureSolver",
    "GenerationRate",
    "MixedField",
    "compact_poisson_residual",
    "l3_generation_rate",
    "parity_selection_rule",
    "pressure_exterior_tail_bound",
    "pure_swirl_no_go",
]


# --------------------------------------------------------------------------- #
# the free-space pressure Poisson solve, in three dimensions                    #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class AxisymmetricPressureSolver:
    r"""``-\Delta_3 p = f`` on a box, with a **non-periodic** ``z``.

    The three-dimensional axisymmetric Laplacian is
    ``\Delta_3 = \partial_{rr} + r^{-1}\partial_r + \partial_{zz}`` — **not**
    the five-dimensional ``\mathcal L_5`` that the stream-function solver uses.
    The two differ in the radial drift coefficient (``1`` against ``3``) and
    therefore in the axis row: for a smooth even ``f``, the radial part of
    ``\Delta_3`` at ``r = 0`` is ``2f''(0)``, discretised as
    ``4(f_1-f_0)/\Delta r^2``, against ``8(f_1-f_0)/\Delta r^2`` for
    ``\mathcal L_5``.  Confusing them is the single easiest way to get a wrong
    pressure, so the coefficients are named constants and are pinned by a test.

    The outer Dirichlet trace uses the **exact exterior representation** of the
    pressure of a divergence-free field.  With ``-\Delta p = \partial_i u_j
    \partial_j u_i = \partial_i\partial_j(u_iu_j)``,

    .. math::

       p(x) = \frac1{4\pi}\int
         \frac{3n_in_j-\delta_{ij}}{|x-y|^3}\,u_i(y)u_j(y)\,dy ,
       \qquad n = \frac{x-y}{|x-y|} ,

    valid for every ``x`` outside the support of ``u``.  The monopole and dipole
    both vanish — ``\int\partial_i\partial_j(u_iu_j) = 0`` and
    ``\int x_k\partial_i\partial_j(u_iu_j) = 0`` for compactly supported ``u``,
    being total derivatives — so the pressure decays like ``|x|^{-3}`` with the
    quadrupole coefficient ``\int u_ku_l\,dy``.  That is what
    :func:`pressure_exterior_tail_bound` bounds.
    """

    #: Radial drift coefficient of the three-dimensional axisymmetric Laplacian.
    RADIAL_DRIFT: ClassVar[float] = 1.0
    #: Axis-row coefficient: the radial part of ``\Delta_3`` is ``2f''(0)``.
    AXIS_COEFFICIENT: ClassVar[float] = 4.0

    grid: AxisymmetricGrid
    sine: FloatArray
    axial_eigenvalues: FloatArray
    radial_lower: FloatArray
    radial_upper: FloatArray
    radial_diagonal_base: FloatArray
    outer_coupling: float
    pivot_floor: float

    @classmethod
    def build(
        cls,
        grid: AxisymmetricGrid,
        *,
        axis_coefficient: float | None = None,
        radial_drift_coefficient: float | None = None,
        pivot_floor: float = 1.0e-30,
    ) -> "AxisymmetricPressureSolver":
        """Assemble the solver.

        The two coefficients are exposed only so fault-injection tests can
        corrupt the audited values ``4`` and ``1``; production callers must
        leave them alone.  Passing the five-dimensional values ``8`` and ``3``
        here silently solves the wrong equation, which is exactly what the
        fault-injection test checks is detectable.
        """
        if grid.periodic_z:
            raise ValueError("the free-space pressure solve requires non-periodic z")
        interior = grid.nz - 2
        if interior < 1:
            raise ValueError("at least three axial nodes are required")
        axis = cls.AXIS_COEFFICIENT if axis_coefficient is None else axis_coefficient
        drift = (
            cls.RADIAL_DRIFT
            if radial_drift_coefficient is None
            else radial_drift_coefficient
        )

        from .whole_space_gate import _sine_transform_matrix

        dr, dz = grid.dr, grid.dz
        radial_unknowns = grid.nr - 1
        modes = np.arange(1, interior + 1, dtype=np.float64)
        eigenvalues = -4.0 / dz**2 * np.sin(np.pi * modes / (2.0 * (interior + 1.0))) ** 2

        lower = np.zeros(radial_unknowns, dtype=np.float64)
        upper = np.zeros(radial_unknowns, dtype=np.float64)
        diagonal = np.zeros(radial_unknowns, dtype=np.float64)
        diagonal[0] = axis / dr**2
        upper[0] = -axis / dr**2
        for i in range(1, radial_unknowns):
            shift = drift / (2.0 * grid.r[i] * dr)
            lower[i] = -(1.0 / dr**2 - shift)
            diagonal[i] = 2.0 / dr**2
            if i < radial_unknowns - 1:
                upper[i] = -(1.0 / dr**2 + shift)
        outer_shift = drift / (2.0 * grid.r[radial_unknowns - 1] * dr)
        return cls(
            grid=grid,
            sine=_sine_transform_matrix(interior),
            axial_eigenvalues=eigenvalues,
            radial_lower=lower,
            radial_upper=upper,
            radial_diagonal_base=diagonal,
            outer_coupling=float(1.0 / dr**2 + outer_shift),
            pivot_floor=float(pivot_floor),
        )

    # -- the exterior representation ----------------------------------------- #

    def exterior_trace(
        self,
        u_r: FloatArray,
        u_theta: FloatArray,
        u_z: FloatArray,
        *,
        support_threshold: float = 1.0e-12,
        azimuthal_nodes: int = 24,
    ) -> tuple[FloatArray, FloatArray, FloatArray]:
        r"""Dirichlet data on the three outer faces from the exact formula.

        The integral over the source is done in the meridional plane with the
        physical measure ``r\,dr\,dz`` and a midpoint rule in ``\theta``; the
        ``\theta`` integral is exact for the ``\theta``-independent part and
        converges spectrally for the rest, since the integrand is smooth and
        periodic.  ``u^\theta`` enters only through ``u^\theta u^\theta``, which
        is ``\theta``-independent, so the only ``\theta`` dependence comes from
        the direction vector ``n``.
        """
        grid = self.grid
        r_source, z_source = grid.mesh()
        support = (
            np.abs(u_r) + np.abs(u_theta) + np.abs(u_z)
        ) > support_threshold
        if not np.any(support):
            shape_outer = np.zeros(grid.nz, dtype=np.float64)
            return (shape_outer, np.zeros(grid.nr), np.zeros(grid.nr))

        weight = 2.0 * np.pi * r_source * grid.dr * grid.dz / azimuthal_nodes
        phi = (np.arange(azimuthal_nodes) + 0.5) * (2.0 * np.pi / azimuthal_nodes)
        cos, sin = np.cos(phi)[:, None], np.sin(phi)[:, None]

        # Source rings, built once: shape (azimuthal_nodes, active_nodes).
        rs = r_source[support][None, :]
        zs = z_source[support][None, :]
        ws = weight[support][None, :]
        ur, ut, uz = u_r[support][None, :], u_theta[support][None, :], u_z[support][None, :]
        sx, sy, sz = rs * cos, rs * sin, zs + 0.0 * cos
        vx = ur * cos - ut * sin
        vy = ur * sin + ut * cos
        vz = uz + 0.0 * cos
        speed_sq = vx * vx + vy * vy + vz * vz
        weights = ws + 0.0 * cos

        def trace(r_target: FloatArray, z_target: FloatArray) -> FloatArray:
            values = np.zeros(r_target.shape, dtype=np.float64)
            for index in range(r_target.size):
                # Target on the phi = 0 half-plane; axisymmetry makes that general.
                tx, tz = float(r_target.flat[index]), float(z_target.flat[index])
                dx, dy, dz_ = tx - sx, -sy, tz - sz
                safe = np.maximum(np.sqrt(dx * dx + dy * dy + dz_ * dz_), 1.0e-300)
                dot = (dx * vx + dy * vy + dz_ * vz) / safe
                kernel = (3.0 * dot * dot - speed_sq) / safe**3
                values.flat[index] = float(np.sum(kernel * weights)) / (4.0 * math.pi)
            return values

        outer = trace(np.full(grid.nz, float(grid.r[-1])), grid.z)
        low = trace(grid.r, np.full(grid.nr, float(grid.z[0])))
        high = trace(grid.r, np.full(grid.nr, float(grid.z[-1])))
        return outer, low, high

    # -- the solve ------------------------------------------------------------ #

    def solve(
        self,
        source: npt.ArrayLike,
        *,
        trace: tuple[FloatArray, FloatArray, FloatArray] | None = None,
    ) -> FloatArray:
        r"""Return ``p`` on the full grid, solving ``-\Delta_3 p = `` ``source``."""
        from .whole_space_gate import _thomas_batch

        grid = self.grid
        values = grid.validate_field(source, name="source")
        if trace is None:
            outer = np.zeros(grid.nz, dtype=np.float64)
            low = np.zeros(grid.nr, dtype=np.float64)
            high = np.zeros(grid.nr, dtype=np.float64)
        else:
            outer, low, high = trace

        interior = grid.nz - 2
        radial_unknowns = grid.nr - 1
        rhs = values[:radial_unknowns, 1:-1].copy()
        rhs[:, 0] += low[:radial_unknowns] / grid.dz**2
        rhs[:, -1] += high[:radial_unknowns] / grid.dz**2
        rhs[-1, :] += self.outer_coupling * outer[1:-1]

        transformed = rhs @ self.sine.T
        diagonal = (
            self.radial_diagonal_base[:, None] - self.axial_eigenvalues[None, :]
        )
        solution_hat = _thomas_batch(
            self.radial_lower, diagonal, self.radial_upper, transformed,
            pivot_floor=self.pivot_floor,
        )
        inner = (solution_hat @ self.sine.T) * (2.0 / (interior + 1.0))

        pressure = np.zeros(grid.shape, dtype=np.float64)
        pressure[:radial_unknowns, 1:-1] = inner
        pressure[-1, :] = outer
        pressure[:, 0] = low
        pressure[:, -1] = high
        pressure[-1, 0] = 0.5 * (outer[0] + low[-1])
        pressure[-1, -1] = 0.5 * (outer[-1] + high[-1])
        return pressure


def compact_poisson_residual(
    grid: AxisymmetricGrid, pressure: FloatArray, source: FloatArray
) -> FloatArray:
    r"""``|-\Delta_h p - f|`` with the solver's own three-point stencil."""
    dr, dz = grid.dr, grid.dz
    residual = np.zeros(grid.shape, dtype=np.float64)
    interior = pressure[1:-1, 1:-1]
    radial = (
        pressure[2:, 1:-1] - 2.0 * interior + pressure[:-2, 1:-1]
    ) / dr**2 + AxisymmetricPressureSolver.RADIAL_DRIFT / grid.r[1:-1, None] * (
        pressure[2:, 1:-1] - pressure[:-2, 1:-1]
    ) / (2.0 * dr)
    axial = (pressure[1:-1, 2:] - 2.0 * interior + pressure[1:-1, :-2]) / dz**2
    residual[1:-1, 1:-1] = np.abs(-(radial + axial) - source[1:-1, 1:-1])
    # Axis row: the radial part of Delta_3 is 2 f''(0), discretised as
    # AXIS_COEFFICIENT (p_1 - p_0)/dr^2.
    axis_radial = (
        AxisymmetricPressureSolver.AXIS_COEFFICIENT
        * (pressure[1, 1:-1] - pressure[0, 1:-1])
        / dr**2
    )
    axis_axial = (
        pressure[0, 2:] - 2.0 * pressure[0, 1:-1] + pressure[0, :-2]
    ) / dz**2
    residual[0, 1:-1] = np.abs(-(axis_radial + axis_axial) - source[0, 1:-1])
    return residual


def parity_selection_rule(swirl_parity: str, stream_parity: str) -> dict[str, object]:
    r"""Decide, before any computation, whether the pressure term can survive.

    The derivation is in the module docstring of
    :mod:`ns_certificate_lab.mixed_initial_data`.  The pressure integrand has
    parity ``-\sigma_\psi`` in ``z``, so on a ``z``-symmetric domain it
    integrates to zero unless the stream generator is **odd**.  The swirl parity
    never enters.

    ``'none'`` means the generator has no definite parity -- a bump centred away
    from ``z = 0``, say -- and then the rule gives no obstruction.
    """
    allowed = {"odd", "even", "none"}
    if swirl_parity not in allowed or stream_parity not in allowed:
        raise ValueError("parities must be 'odd', 'even' or 'none'")
    if stream_parity == "even":
        verdict = False
        reason = (
            "the stream generator is even in z, so the pressure integrand is odd "
            "and cancels exactly on a z-symmetric domain"
        )
    elif stream_parity == "odd":
        verdict = True
        reason = (
            "the stream generator is odd in z, so the pressure integrand is even "
            "and does not cancel"
        )
    else:
        verdict = True
        reason = (
            "the stream generator has no definite parity, so no cancellation is "
            "forced"
        )
    return {
        "swirl_parity": swirl_parity,
        "stream_parity": stream_parity,
        "pressure_term_can_be_nonzero": verdict,
        "reason": reason,
        "swirl_parity_is_irrelevant": True,
    }


def pressure_exterior_tail_bound(
    kinetic_l2_squared: float, distance: float
) -> float:
    r"""Rigorous bound on ``|p|`` outside the support of ``u``.

    From the exterior representation, ``|3n_in_j-\delta_{ij}|`` has operator
    norm ``2`` (eigenvalues ``2, -1, -1``), so

    .. math::
       |p(x)| \le \frac{2}{4\pi}\int\frac{|u(y)|^2}{|x-y|^3}dy
              \le \frac{\|u\|_{L^2}^2}{2\pi\,d^3},
       \qquad d = \operatorname{dist}(x,\operatorname{supp}u).

    This is a genuine upper bound, not an asymptotic: it holds at every finite
    distance outside the support.
    """
    if distance <= 0.0:
        raise ValueError("the distance to the support must be positive")
    if kinetic_l2_squared < 0.0:
        raise ValueError("a squared L2 norm cannot be negative")
    return kinetic_l2_squared / (2.0 * math.pi * distance**3)


# --------------------------------------------------------------------------- #
# the mixed field and its gradients                                            #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MixedField:
    r"""An axisymmetric velocity with **both** swirl and meridional flow.

    Built from two scalar generators through the audited recovery

    .. math::
       u^\theta = r u_1,\quad u^r = -r\partial_z\psi_1,\quad
       u^z = 2\psi_1 + r\partial_r\psi_1 ,

    which is divergence free for **any** ``\psi_1``:
    ``\partial_r u^r + u^r/r + \partial_z u^z
    = -r\partial_{rz}\psi_1 - \partial_z\psi_1
      + 2\partial_z\psi_1 + \partial_z\psi_1 + r\partial_{rz}\psi_1 = 0``.

    The nine physical components of ``\nabla u`` in the orthonormal cylindrical
    frame are

    ======  ==========================  ======================
    entry   value                       entry
    ======  ==========================  ======================
    ``rr``  ``\partial_r u^r``          ``\theta r``: ``-u^\theta/r``
    ``r\theta``  ``\partial_r u^\theta``  ``\theta\theta``: ``u^r/r``
    ``rz``  ``\partial_r u^z``          ``\theta z``: ``0``
    ``zr``  ``\partial_z u^r``          ``zz``: ``\partial_z u^z``
    ``z\theta``  ``\partial_z u^\theta``
    ======  ==========================  ======================

    The two ``1/r`` entries are the ones a Cartesian-blind implementation gets
    wrong; they are what makes ``|\nabla u|^2`` differ from the naive sum of
    squared meridional derivatives.
    """

    grid: AxisymmetricGrid
    u_r: FloatArray
    u_theta: FloatArray
    u_z: FloatArray
    #: Supplied when the generators can differentiate themselves exactly.  A
    #: finite-difference gradient of a chi-bump field carries a large truncation
    #: error -- the bump has enormous high derivatives -- and that error would
    #: land directly in ``J``, so an analytic gradient is used whenever one is
    #: available and the difference operators are kept only as a cross-check.
    analytic_gradient: dict[str, FloatArray] | None = None

    @staticmethod
    def from_generators(
        grid: AxisymmetricGrid, u1: npt.ArrayLike, psi1: npt.ArrayLike
    ) -> "MixedField":
        r_mesh, _z = grid.mesh()
        swirl = grid.validate_field(u1, name="u1")
        stream = grid.validate_field(psi1, name="psi1")
        psi_r = derivative_r(grid, stream, even_at_axis=True)
        psi_z = derivative_z(grid, stream)
        u_r = -r_mesh * psi_z
        u_r[0] = 0.0
        u_z = 2.0 * stream + r_mesh * psi_r
        return MixedField(grid=grid, u_r=u_r, u_theta=r_mesh * swirl, u_z=u_z)

    @property
    def speed(self) -> FloatArray:
        return np.sqrt(self.u_r**2 + self.u_theta**2 + self.u_z**2)

    def gradient(self) -> dict[str, FloatArray]:
        """The nine physical components, as a dictionary keyed ``'rr'`` etc.

        Returns the analytic gradient when the constructor was given one; the
        finite-difference fallback below is what a field with no closed form
        gets, and :meth:`finite_difference_gradient` exposes it separately so
        the two can be compared.
        """
        if self.analytic_gradient is not None:
            return self.analytic_gradient
        return self.finite_difference_gradient()

    def finite_difference_gradient(self) -> dict[str, FloatArray]:
        """The same nine components from the audited difference operators."""
        grid = self.grid
        r_mesh, _z = grid.mesh()
        safe_r = r_mesh.copy()
        safe_r[0] = 1.0
        # u^r and u^theta vanish on the axis like r, so the quotients are finite;
        # their limits are d_r u^r and d_r u^theta.
        u_r_over_r = np.where(r_mesh > 0.0, self.u_r / safe_r, 0.0)
        u_t_over_r = np.where(r_mesh > 0.0, self.u_theta / safe_r, 0.0)
        d_r_ur = derivative_r(grid, self.u_r)
        d_r_ut = derivative_r(grid, self.u_theta)
        u_r_over_r[0] = d_r_ur[0]
        u_t_over_r[0] = d_r_ut[0]
        return {
            "rr": d_r_ur,
            "rt": d_r_ut,
            "rz": derivative_r(grid, self.u_z, even_at_axis=True),
            "tr": -u_t_over_r,
            "tt": u_r_over_r,
            "tz": np.zeros_like(self.u_r),
            "zr": derivative_z(grid, self.u_r),
            "zt": derivative_z(grid, self.u_theta),
            "zz": derivative_z(grid, self.u_z),
        }

    def gradient_norm_squared(self, gradient: dict[str, FloatArray] | None = None):
        entries = self.gradient() if gradient is None else gradient
        return sum(value * value for value in entries.values())

    def divergence(self, gradient: dict[str, FloatArray] | None = None) -> FloatArray:
        entries = self.gradient() if gradient is None else gradient
        return entries["rr"] + entries["tt"] + entries["zz"]

    def pressure_source(
        self, gradient: dict[str, FloatArray] | None = None
    ) -> FloatArray:
        r"""``\partial_i u_j\partial_j u_i = \operatorname{tr}((\nabla u)^2)``."""
        g = self.gradient() if gradient is None else gradient
        return (
            g["rr"] ** 2
            + g["tt"] ** 2
            + g["zz"] ** 2
            + 2.0 * g["rt"] * g["tr"]
            + 2.0 * g["rz"] * g["zr"]
            + 2.0 * g["tz"] * g["zt"]
        )

    def volume_element(self) -> FloatArray:
        r"""``2\pi r\,\Delta r\,\Delta z`` — the **physical three-dimensional**
        measure, not the ``r^3`` measure of the five-dimensional lift."""
        r_mesh, _z = self.grid.mesh()
        return 2.0 * math.pi * r_mesh * self.grid.dr * self.grid.dz

    def kinetic_energy(self) -> float:
        return 0.5 * float(np.sum(self.speed**2 * self.volume_element()))

    def l3_cubed(self) -> float:
        return float(np.sum(self.speed**3 * self.volume_element()))

    def support_radius(self, threshold: float = 1.0e-12) -> float:
        r_mesh, z_mesh = self.grid.mesh()
        mask = self.speed > threshold
        if not np.any(mask):
            return 0.0
        return float(np.max(np.hypot(r_mesh[mask], z_mesh[mask])))


# --------------------------------------------------------------------------- #
# the generation rate                                                          #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class GenerationRate:
    r"""``J(u_0) = \frac{d}{dt}\|u\|_{L^3}^3\big|_{t=0}`` and its parts."""

    total: float
    pressure: float
    viscous: float
    transport_residual: float
    l3_cubed: float
    kinetic_energy: float
    max_speed: float
    max_divergence: float
    relative_divergence: float
    poisson_residual: float
    pressure_tail_bound: float
    support_radius: float
    outer_gap: float
    viscosity: float

    @property
    def normalised(self) -> float:
        """``J`` divided by ``\\|u\\|_{L^3}^3`` — a rate with units of inverse time."""
        return self.total / max(self.l3_cubed, 1.0e-300)

    @property
    def positive(self) -> bool:
        return self.total > 0.0

    def as_dict(self) -> dict[str, float | bool]:
        payload: dict[str, float | bool] = {
            name: float(getattr(self, name)) for name in self.__dataclass_fields__
        }
        payload["normalised"] = float(self.normalised)
        payload["positive"] = bool(self.positive)
        return payload


def l3_generation_rate(
    field: MixedField,
    *,
    viscosity: float,
    solver: AxisymmetricPressureSolver | None = None,
    support_threshold: float = 1.0e-12,
    azimuthal_nodes: int = 24,
    interior_margin: int = 2,
) -> GenerationRate:
    r"""Evaluate ``J(u_0)`` with the three contributions kept apart.

    The pressure contribution is evaluated in the integrated-by-parts form

    .. math::  3\int p\,\nabla\!\cdot\!(|u|u)\,dx = -3\int |u|\,u\cdot\nabla p\,dx ,

    which avoids dividing by ``|u|`` in the dominant term.  The transport
    contribution is *known* to be zero; it is computed anyway and returned as
    ``transport_residual``, because a residual that is not small means the
    quadrature or the divergence-free property has failed and the other two
    numbers cannot be trusted either.
    """
    if viscosity <= 0.0:
        raise ValueError("a Clay candidate requires a fixed positive viscosity")
    grid = field.grid
    gradient = field.gradient()
    speed = field.speed
    volume = field.volume_element()

    pressure_solver = (
        AxisymmetricPressureSolver.build(grid) if solver is None else solver
    )
    source = field.pressure_source(gradient)
    trace = pressure_solver.exterior_trace(
        field.u_r, field.u_theta, field.u_z,
        support_threshold=support_threshold, azimuthal_nodes=azimuthal_nodes,
    )
    pressure = pressure_solver.solve(source, trace=trace)

    interior = slice(interior_margin, -interior_margin)
    mask = np.zeros(grid.shape, dtype=bool)
    mask[interior, interior] = True

    p_r = derivative_r(grid, pressure, even_at_axis=True)
    p_z = derivative_z(grid, pressure)
    pressure_density = -3.0 * speed * (field.u_r * p_r + field.u_z * p_z)

    # d_k q where q = |u|:  q d_k q = (1/2) d_k(q^2) is smooth, so the quotient
    # below is the only division by q and it is dominated by |u||grad u|^2.
    square = speed**2
    g_r = 0.5 * derivative_r(grid, square, even_at_axis=True)
    g_z = 0.5 * derivative_z(grid, square)
    gradient_square = field.gradient_norm_squared(gradient)
    safe_speed = np.where(speed > 0.0, speed, 1.0)
    # Kato: |grad |u|| <= |grad u| pointwise a.e., so q |grad q|^2 <= q |grad u|^2.
    # Taking the minimum is therefore the sharp analytic bound rather than a
    # numerical fudge, and it removes the 0/0 on the zero set of u.
    grad_speed_term = np.minimum(
        np.where(speed > 0.0, (g_r**2 + g_z**2) / safe_speed, 0.0),
        speed * gradient_square,
    )
    viscous_density = -3.0 * viscosity * (speed * gradient_square + grad_speed_term)

    # The transport term, which the identity says is exactly zero.
    advection_r = (
        field.u_r * gradient["rr"] + field.u_z * gradient["zr"]
        - field.u_theta**2 / np.where(grid.mesh()[0] > 0.0, grid.mesh()[0], 1.0)
    )
    advection_t = (
        field.u_r * gradient["rt"] + field.u_z * gradient["zt"]
        + field.u_r * field.u_theta / np.where(
            grid.mesh()[0] > 0.0, grid.mesh()[0], 1.0
        )
    )
    advection_z = field.u_r * gradient["rz"] + field.u_z * gradient["zz"]
    transport_density = -3.0 * speed * (
        field.u_r * advection_r + field.u_theta * advection_t + field.u_z * advection_z
    )

    # The residual uses the solver's OWN three-point stencil, so it verifies the
    # linear solve.  A wider stencil would instead re-measure truncation error,
    # which on chi-bump data is dominated by the fourth derivative of chi and
    # says nothing about whether the solve is right.
    residual = compact_poisson_residual(grid, pressure, source)

    divergence = field.divergence(gradient)
    divergence_scale = max(float(np.max(np.abs(gradient["rr"]))), 1.0e-300)
    radius = field.support_radius(support_threshold)
    gap = min(float(grid.r[-1]), float(grid.z[-1])) - radius
    energy = field.kinetic_energy()
    tail = (
        pressure_exterior_tail_bound(2.0 * energy, gap) if gap > 0.0 else math.inf
    )

    return GenerationRate(
        total=float(np.sum((pressure_density + viscous_density) * volume)),
        pressure=float(np.sum(pressure_density * volume)),
        viscous=float(np.sum(viscous_density * volume)),
        transport_residual=float(np.sum(transport_density * volume)),
        l3_cubed=field.l3_cubed(),
        kinetic_energy=energy,
        max_speed=float(np.max(speed)),
        max_divergence=float(np.max(np.abs(divergence[mask]))),
        relative_divergence=float(np.max(np.abs(divergence[mask]))) / divergence_scale,
        poisson_residual=float(np.max(residual[mask])),
        pressure_tail_bound=float(tail) if math.isfinite(tail) else -1.0,
        support_radius=radius,
        outer_gap=gap,
        viscosity=float(viscosity),
    )


def pure_swirl_no_go(field: MixedField, *, tolerance: float = 1.0e-12) -> dict:
    r"""Certify that a datum is a pure swirl, and record what that excludes.

    A pure swirl is ``u^r = u^z = 0``.  For such a datum ``|u| = |u^\theta(r,z)|``
    depends only on ``(r,z)``, so ``|u|u = g(r,z)e_\theta`` and

    .. math::  \nabla\cdot(|u|u) = r^{-1}\partial_\theta g = 0 .

    The pressure contribution to ``J`` therefore vanishes *identically*, leaving
    ``J = -3\nu\int(|u||\nabla u|^2 + |u||\nabla|u||^2)\le 0``, with equality
    only for ``u\equiv 0``.

    **What this does not exclude.**  It is a statement at ``t = 0``.  The
    solution immediately generates meridional flow through vortex stretching, so
    ``\nabla\cdot(|u|u)`` becomes nonzero for ``t>0`` and the pressure term
    switches on.  Nothing here forbids ``\|u(t)\|_{L^3}`` from growing later, and
    nothing here bears on the Clay problem.
    """
    meridional = max(
        float(np.max(np.abs(field.u_r))), float(np.max(np.abs(field.u_z)))
    )
    swirl = float(np.max(np.abs(field.u_theta)))
    is_pure = meridional <= tolerance * max(swirl, 1.0)
    return {
        "is_pure_swirl": bool(is_pure),
        "max_meridional_speed": meridional,
        "max_swirl_speed": swirl,
        "pressure_contribution_vanishes": bool(is_pure),
        "initial_rate_sign": "non-positive" if is_pure else "undetermined",
        "excludes": (
            "a positive initial generation rate of the global critical norm"
            if is_pure
            else None
        ),
        "does_not_exclude": [
            "growth of the critical norm at any later time",
            "meridional flow generated by vortex stretching for t > 0",
            "anything about the Clay problem",
        ],
    }
