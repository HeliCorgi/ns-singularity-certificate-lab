r"""Derivative tail bounds and free-space velocity recovery.

Gate 4 (``whole_space_gate.py``) certified the free-space potential **value**:
its a posteriori bound covers ``|psi - psi_mono|`` and propagates into the box by
the maximum principle.  That is not enough to enter a nonlinear whole-space run,
because the momentum equation consumes *derivatives*:

    u^r = -r psi1_z,   u^z = 2 psi1 + r psi1_r,

and the advection term consumes derivatives of those again.  A maximum principle
says nothing about derivatives.  This module supplies the missing half.

Two independent constructions
-----------------------------
* **Analytic Green differentiation.**  The five-dimensional fundamental solution
  ``G_5(X) = 1/(8 pi^2 |X|^3)`` is ``(-3)``-homogeneous, so ``D^m G_5`` is
  ``(-3-m)``-homogeneous and ``||D^m G_5(Z)|| = A_m |Z|^{-3-m}`` **exactly**,
  with ``A_m`` a pure number computed once from the closed-form derivative
  tensors.  Taylor's theorem with the mean-value remainder then bounds every
  derivative of the multipole truncation error.  No maximum principle is used.
* **Interior elliptic estimate.**  The box-truncation error is ``L_5``-harmonic,
  so the classical interior estimate ``|D^k h(x)| <= (n k / r)^k sup_{B_r} |h|``
  with ``n = 5`` converts the boundary-data bound into interior derivative
  bounds.  This *is* an elliptic estimate rather than a maximum principle, and
  it is the only place where the box enters.

Both are continuum statements about ``psi_infinity`` and the continuum box
solution.  The discretisation error of the finite-difference solve is measured,
not bounded — the radial row ``i=1`` is not an M-matrix and this repository has
no validated numerics.

Recovery paths
--------------
``recover_free_space_velocity`` differentiates the *solver output* with the
audited finite-difference operators; that is the path a time-stepper uses.
``green_quadrature`` instead differentiates the *kernel* analytically and
quadratures the source, sharing no code with the solver or the difference
operators.  It is restricted to target points outside the source support, where
the reduced kernel is smooth; that is exactly the region the tail bounds
describe.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .whole_space_gate import SPHERE3_MEASURE, SPHERE4_MEASURE

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "DerivativeTailBounds",
    "FreeSpaceRecovery",
    "GREEN_DERIVATIVE_CONSTANTS",
    "SourceMoments",
    "green_derivative_constant",
    "green_quadrature",
    "harmonic_interior_factor",
    "interior_derivative_tail_bounds",
    "multipole_boundary_trace",
    "multipole_tail_bounds",
    "recover_free_space_velocity",
    "source_moments",
]

#: ``c_5`` in ``G_5(X) = c_5 |X|^{-3}``.
GREEN_CONSTANT = 1.0 / (8.0 * math.pi**2)

#: ``A_m`` with ``||D^m G_5(Z)|| = A_m |Z|^{-3-m}``, as multiples of ``c_5``.
#:
#: ``A_0 = 1`` and ``A_1 = 3`` are exact.  ``A_2 = 12`` is exact: the Hessian of
#: ``|Z|^{-3}`` is ``-3|Z|^{-5} I + 15 |Z|^{-7} Z (x) Z``, whose eigenvalues are
#: ``12|Z|^{-5}`` along ``Z`` and ``-3|Z|^{-5}`` transverse.  ``A_3 = 150`` and
#: ``A_4 = 1620`` are triangle-inequality upper bounds on the symmetric
#: multilinear forms obtained from the closed-form derivative tensors
#:
#:   D^3 = 15(delta Z + delta Z + delta Z)|Z|^-7 - 105 ZZZ |Z|^-9,
#:   D^4 = 15(3 delta delta)|Z|^-7 - 105(6 delta ZZ)|Z|^-9 + 945 ZZZZ |Z|^-11,
#:   D^5 = -315(delta delta Z)|Z|^-9 + 105(12 delta delta Z)|Z|^-9
#:         + 945(6+4 delta ZZZ)|Z|^-11 - 10395 ZZZZZ |Z|^-13,
#:
#: giving 45+105 = 150, 45+630+945 = 1620 and 315+1260+5670+3780+10395 = 21420.
#: The observed sharp value of ``A_3`` is ``60``, so the high-order bounds are
#: conservative by a known factor rather than by "some constant".
#: ``A_6 <= 335160`` follows from the same rule: differentiating a term
#: ``c delta^a Z^b |Z|^{-p}`` produces terms of total magnitude at most
#: ``c(b+p)`` on unit vectors.  Applied to the ``A_4`` term list that rule
#: reproduces ``A_5 = 21420`` exactly, which is why it is trusted for ``A_6``.
_GREEN_DERIVATIVE_MULTIPLES: tuple[float, ...] = (
    1.0, 3.0, 12.0, 150.0, 1620.0, 21420.0, 335160.0,
)
GREEN_DERIVATIVE_CONSTANTS: tuple[float, ...] = tuple(
    GREEN_CONSTANT * value for value in _GREEN_DERIVATIVE_MULTIPLES
)


def green_derivative_constant(order: int) -> float:
    """``A_m`` such that ``||D^m G_5(Z)|| <= A_m |Z|^{-3-m}``."""
    if order < 0 or order >= len(GREEN_DERIVATIVE_CONSTANTS):
        raise ValueError(
            f"order must lie in 0..{len(GREEN_DERIVATIVE_CONSTANTS) - 1}"
        )
    return GREEN_DERIVATIVE_CONSTANTS[order]


# --------------------------------------------------------------------------- #
# Source moments                                                               #
# --------------------------------------------------------------------------- #


def _trapezoid(values: FloatArray, spacing: float, axis: int) -> FloatArray:
    moved = np.moveaxis(np.asarray(values, dtype=np.float64), axis, 0)
    weights = np.full(moved.shape[0], spacing, dtype=np.float64)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    return np.tensordot(weights, moved, axes=(0, 0))


@dataclass(frozen=True)
class SourceMoments:
    """Moments of a compactly supported axisymmetric source, in ``dV_5``.

    ``dV_5 = |S^3| rho^3 d rho d zeta`` for an axisymmetric integrand.  All
    entries are computed from grid data by the trapezoid rule; none of them
    reads an analytic reference, so a boundary condition built from them stays
    non-circular.
    """

    absolute_mass: float
    first_absolute_moment: float
    second_absolute_moment: float
    monopole: float
    axial_dipole: float
    support_radius: float
    neglected_absolute_mass: float
    third_absolute_moment: float
    axial_quadrupole: float
    transverse_quadrupole: float

    @property
    def newtonian_mass(self) -> float:
        """``M = m(infinity)``: the far field is ``M/(3R^3)``."""
        return self.monopole / SPHERE4_MEASURE

    def as_dict(self) -> dict[str, float]:
        return {
            "absolute_mass": self.absolute_mass,
            "first_absolute_moment": self.first_absolute_moment,
            "second_absolute_moment": self.second_absolute_moment,
            "monopole": self.monopole,
            "axial_dipole": self.axial_dipole,
            "support_radius": self.support_radius,
            "third_absolute_moment": self.third_absolute_moment,
            "axial_quadrupole": self.axial_quadrupole,
            "transverse_quadrupole": self.transverse_quadrupole,
            "neglected_absolute_mass": self.neglected_absolute_mass,
            "neglected_fraction": self.neglected_fraction,
        }

    @property
    def neglected_fraction(self) -> float:
        """Share of ``int |omega| dV_5`` outside the reported support radius.

        A positive threshold makes every bound in this module **conditional**:
        the derivation assumes the source vanishes outside ``support_radius``.
        This number is what the assumption costs, and callers must check it.
        """
        if self.absolute_mass <= 0.0:
            return 0.0
        return self.neglected_absolute_mass / self.absolute_mass


def source_moments(
    grid: AxisymmetricGrid,
    omega1: npt.ArrayLike,
    *,
    support_threshold: float = 0.0,
) -> SourceMoments:
    """Compute the moments Gate-5 needs, including the support radius.

    ``support_radius`` is the smallest five-dimensional radius containing every
    node where ``|omega1|`` exceeds ``support_threshold`` times its maximum.  A
    strictly positive threshold is a *choice* and makes the bounds conditional
    on it, so the default is exact support.
    """
    values = grid.validate_field(omega1, name="omega1")
    r_mesh, z_mesh = grid.mesh()
    radius5 = np.hypot(r_mesh, z_mesh)
    weight = grid.r[:, None] ** 3

    def integrate(field: FloatArray) -> float:
        return float(
            _trapezoid(_trapezoid(field * weight, grid.dz, axis=1), grid.dr, axis=0)
        )

    magnitude = np.abs(values)
    peak = float(magnitude.max())
    active = magnitude > support_threshold * peak if peak > 0.0 else magnitude > 0.0
    support = float(radius5[active].max()) if np.any(active) else 0.0
    outside = radius5 > support
    return SourceMoments(
        neglected_absolute_mass=SPHERE3_MEASURE
        * integrate(np.where(outside, magnitude, 0.0)),
        absolute_mass=SPHERE3_MEASURE * integrate(magnitude),
        first_absolute_moment=SPHERE3_MEASURE * integrate(magnitude * radius5),
        second_absolute_moment=SPHERE3_MEASURE * integrate(magnitude * radius5**2),
        third_absolute_moment=SPHERE3_MEASURE * integrate(magnitude * radius5**3),
        monopole=SPHERE3_MEASURE * integrate(values),
        axial_dipole=SPHERE3_MEASURE * integrate(values * z_mesh),
        axial_quadrupole=SPHERE3_MEASURE * integrate(values * z_mesh**2),
        # Q_ij = q_t delta_ij on the four transverse axes, so
        # 4 q_t = int |y'|^2 omega dV_5.
        transverse_quadrupole=0.25 * SPHERE3_MEASURE * integrate(values * r_mesh**2),
        support_radius=support,
    )


# --------------------------------------------------------------------------- #
# Main work A: derivative tail bounds                                          #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DerivativeTailBounds:
    """Bounds on ``|D^k(psi - psi_truncated)|`` for ``k = 0,1,2,3``."""

    truncation: str
    distance: float
    value: float
    gradient: float
    hessian: float
    third: float

    def __getitem__(self, order: int) -> float:
        return (self.value, self.gradient, self.hessian, self.third)[order]

    @property
    def available(self) -> bool:
        """False when the source reaches the evaluation sphere and no bound applies.

        A box tight enough for the outer boundary treatment to matter is often
        too tight for the multipole bound to exist at all.  That tension is real
        and is reported rather than papered over with a fabricated number.
        """
        return math.isfinite(self.value)

    @staticmethod
    def unavailable(truncation: str, distance: float) -> "DerivativeTailBounds":
        return DerivativeTailBounds(
            truncation=truncation, distance=distance,
            value=math.inf, gradient=math.inf, hessian=math.inf, third=math.inf,
        )

    def as_dict(self) -> dict[str, object]:
        def finite(value: float) -> float | None:
            return value if math.isfinite(value) else None

        return {
            "truncation": self.truncation,
            "distance": self.distance,
            "available": self.available,
            "value": finite(self.value),
            "gradient": finite(self.gradient),
            "hessian": finite(self.hessian),
            "third": finite(self.third),
        }


def multipole_tail_bounds(
    moments: SourceMoments,
    *,
    evaluation_radius: float,
    truncation: str = "monopole",
) -> DerivativeTailBounds:
    r"""Explicit bounds on every derivative of the multipole truncation error.

    Let the source be supported in ``|Y| <= a`` and set ``d = |X| - a > 0``.

    **Monopole truncation.**  ``psi(X) - M_5 G_5(X) = int [G_5(X-Y) - G_5(X)]
    omega(Y) dY``, so differentiating ``k`` times under the integral and applying
    the mean value theorem along the segment ``X - theta Y``,

    .. math::

        |D^k \psi - D^k(M_5 G_5)|(X)
          \le A_{k+1}\, \frac{I_1}{d^{\,4+k}},\qquad
        I_1 = \int |Y|\,|\omega(Y)|\,dV_5 ,

    because ``|X - theta Y| >= |X| - a = d`` and ``||D^{k+1}G_5(Z)||
    = A_{k+1}|Z|^{-4-k}`` exactly, by homogeneity.

    **Dipole truncation.**  Subtracting also ``-P . grad G_5(X)`` with
    ``P = int Y omega dV_5`` and using the second-order Taylor remainder,

    .. math::

        |D^k \psi - D^k(M_5 G_5 - P\cdot\nabla G_5)|(X)
          \le \tfrac12 A_{k+2}\, \frac{I_2}{d^{\,5+k}},\qquad
        I_2 = \int |Y|^2 |\omega|\,dV_5 .

    Neither line uses a maximum principle: both come from differentiating the
    Green kernel and bounding the remainder of a Taylor expansion.
    """
    if truncation not in {"monopole", "dipole", "quadrupole"}:
        raise ValueError("truncation must be 'monopole', 'dipole' or 'quadrupole'")
    distance = float(evaluation_radius) - moments.support_radius
    if distance <= 0.0:
        raise ValueError(
            "evaluation_radius must exceed the source support radius; the "
            "source has spread to the evaluation sphere and no multipole tail "
            "bound applies"
        )
    if truncation == "monopole":
        moment = moments.first_absolute_moment
        shift, factor = 1, 1.0
    elif truncation == "dipole":
        moment = moments.second_absolute_moment
        shift, factor = 2, 0.5
    else:
        moment = moments.third_absolute_moment
        shift, factor = 3, 1.0 / 6.0
    values = [
        factor
        * green_derivative_constant(order + shift)
        * moment
        / distance ** (3 + shift + order)
        for order in range(4)
    ]
    return DerivativeTailBounds(
        truncation=truncation,
        distance=distance,
        value=values[0],
        gradient=values[1],
        hessian=values[2],
        third=values[3],
    )


def harmonic_interior_factor(order: int, distance: float) -> float:
    r"""``(n k / rho)^k`` with ``n = 5``: the interior estimate for harmonic ``h``.

    For ``h`` harmonic on ``B_rho(x) subset R^n``, ``|D^k h(x)| <= (n k /
    rho)^k sup_{B_rho} |h|``.  This is an interior elliptic estimate, not a
    maximum principle: the maximum principle only gives the ``k = 0`` case.
    """
    if order < 0:
        raise ValueError("order must be nonnegative")
    if distance <= 0.0:
        raise ValueError("distance must be positive")
    if order == 0:
        return 1.0
    return (5.0 * order / distance) ** order


def interior_derivative_tail_bounds(
    boundary_value_bound: float,
    *,
    distance_to_boundary: float,
) -> DerivativeTailBounds:
    r"""Propagate a boundary-data bound to interior derivative bounds.

    The box-truncation error ``e = psi_box - psi_infinity`` satisfies
    ``L_5 e = 0`` inside the box, so it is harmonic in ``R^5``.  With
    ``sup_{\partial Box} |e| <= eps_0`` the maximum principle gives ``|e| <=
    eps_0`` and the interior estimate gives ``|D^k e| <= (5k/rho)^k eps_0`` at
    distance ``rho`` from the boundary.

    Continuum statement only.  The discrete solve adds a discretisation error
    that this does not bound.
    """
    if boundary_value_bound < 0.0:
        raise ValueError("boundary_value_bound must be nonnegative")
    return DerivativeTailBounds(
        truncation="box_truncation",
        distance=float(distance_to_boundary),
        value=boundary_value_bound,
        gradient=harmonic_interior_factor(1, distance_to_boundary)
        * boundary_value_bound,
        hessian=harmonic_interior_factor(2, distance_to_boundary)
        * boundary_value_bound,
        third=harmonic_interior_factor(3, distance_to_boundary)
        * boundary_value_bound,
    )


def multipole_boundary_trace(
    grid: AxisymmetricGrid,
    moments: SourceMoments,
    radius: FloatArray,
    axial: FloatArray,
    *,
    truncation: str = "monopole",
) -> FloatArray:
    r"""Evaluate the multipole approximation at the requested points.

    ``psi ~ M_5 G_5(X) - P_z d_z G_5(X)`` with ``G_5 = c_5 |X|^{-3}`` and
    ``d_z G_5 = -3 c_5 z |X|^{-5}``, so the dipole term is
    ``+3 c_5 P_z z |X|^{-5}``.  Transverse dipole components vanish by
    axisymmetry, which is why only ``P_z`` appears.
    """
    if truncation not in {"monopole", "dipole", "quadrupole"}:
        raise ValueError("truncation must be 'monopole', 'dipole' or 'quadrupole'")
    radius5 = np.hypot(
        np.asarray(radius, dtype=np.float64), np.asarray(axial, dtype=np.float64)
    )
    safe = np.maximum(radius5, 1.0e-300)
    axial_values = np.asarray(axial, dtype=np.float64)
    out = moments.monopole * GREEN_CONSTANT / safe**3
    if truncation in {"dipole", "quadrupole"}:
        out = out + 3.0 * GREEN_CONSTANT * moments.axial_dipole * axial_values / safe**5
    if truncation == "quadrupole":
        # (1/2) Q : D^2 G reduces, by axisymmetry and harmonicity of G away from
        # the origin, to (1/2)(q_z - q_t)(15 z^2 |X|^-7 - 3 |X|^-5) c_5.
        gap = moments.axial_quadrupole - moments.transverse_quadrupole
        out = out + 0.5 * GREEN_CONSTANT * gap * (
            15.0 * axial_values**2 / safe**7 - 3.0 / safe**5
        )
    return out


# --------------------------------------------------------------------------- #
# The analytic-Green quadrature path                                           #
# --------------------------------------------------------------------------- #


def _gauss_chebyshev_second_kind(nodes: int) -> tuple[FloatArray, FloatArray]:
    """Nodes and weights for ``int_{-1}^{1} f(t) sqrt(1-t^2) dt``.

    Exact for polynomials of degree ``2*nodes - 1``.  Building the rule here
    keeps the Green path free of any shared quadrature code.
    """
    if nodes < 2:
        raise ValueError("at least two quadrature nodes are required")
    index = np.arange(1, nodes + 1, dtype=np.float64)
    angle = index * math.pi / (nodes + 1.0)
    return np.cos(angle), math.pi / (nodes + 1.0) * np.sin(angle) ** 2


def green_quadrature(
    grid: AxisymmetricGrid,
    omega1: npt.ArrayLike,
    targets_r: npt.ArrayLike,
    targets_z: npt.ArrayLike,
    *,
    quadrature_nodes: int = 48,
) -> dict[str, FloatArray]:
    r"""Evaluate ``psi`` and its derivatives by differentiating the kernel.

    For an axisymmetric source the ``S^3`` average of the kernel reduces to a
    one-dimensional integral,

    .. math::

        K(r,z;\rho,\zeta) = 4\pi c_5 \int_{-1}^{1}
            \bigl(r^2+\rho^2-2r\rho t+(z-\zeta)^2\bigr)^{-3/2}\sqrt{1-t^2}\,dt,

    using ``int_{S^3} f(\hat x\cdot\hat y)\,d\sigma = |S^2| \int_{-1}^1 f(t)
    \sqrt{1-t^2}\,dt``; the weight is matched exactly by Gauss--Chebyshev of the
    second kind.  Every derivative below differentiates the kernel
    **analytically** before quadrature — no finite differences and no solver.

    Valid only for target points outside the source support, where the kernel is
    smooth; that is the region the tail bounds describe.  This routine does not
    check that condition, because "outside the support" depends on the support
    threshold the caller chose.
    """
    values = grid.validate_field(omega1, name="omega1")
    r_target = np.asarray(targets_r, dtype=np.float64).ravel()
    z_target = np.asarray(targets_z, dtype=np.float64).ravel()
    if r_target.shape != z_target.shape:
        raise ValueError("target coordinate arrays must have the same shape")

    nodes, weights = _gauss_chebyshev_second_kind(quadrature_nodes)
    rho = grid.r
    zeta = grid.z
    # Source quadrature weights: trapezoid in (rho, zeta) with the rho^3 measure.
    wr = np.full(rho.size, grid.dr)
    wr[0] *= 0.5
    wr[-1] *= 0.5
    wz = np.full(zeta.size, grid.dz)
    wz[0] *= 0.5
    wz[-1] *= 0.5
    density = values * (rho[:, None] ** 3) * wr[:, None] * wz[None, :]

    prefactor = 4.0 * math.pi * GREEN_CONSTANT
    out = {
        key: np.zeros(r_target.shape, dtype=np.float64)
        for key in ("psi", "dr", "dz", "drr", "dzz", "drz")
    }
    for index in range(r_target.size):
        r0 = r_target[index]
        z0 = z_target[index]
        dz = z0 - zeta[None, :]
        # shape (n_rho, n_zeta, n_nodes)
        cross = r0 * rho[:, None, None] * nodes[None, None, :]
        squared = (
            r0**2
            + rho[:, None, None] ** 2
            - 2.0 * cross
            + dz[:, :, None] ** 2
        )
        distance = np.sqrt(np.maximum(squared, 1.0e-300))
        inv3 = distance**-3
        inv5 = distance**-5
        inv7 = distance**-7
        dr_component = r0 - rho[:, None, None] * nodes[None, None, :]
        dz_component = dz[:, :, None]
        kernels = {
            "psi": inv3,
            "dr": -3.0 * dr_component * inv5,
            "dz": -3.0 * dz_component * inv5,
            "drr": -3.0 * inv5 + 15.0 * dr_component**2 * inv7,
            "dzz": -3.0 * inv5 + 15.0 * dz_component**2 * inv7,
            "drz": 15.0 * dr_component * dz_component * inv7,
        }
        for key, kernel in kernels.items():
            reduced = prefactor * np.tensordot(kernel, weights, axes=(2, 0))
            out[key][index] = float(np.sum(reduced * density))
    return out


# --------------------------------------------------------------------------- #
# Main work B: the free-space velocity recovery API                            #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class FreeSpaceRecovery:
    """Everything a whole-space time-stepper needs from one elliptic solve."""

    psi1: FloatArray
    dpsi1_dr: FloatArray
    dpsi1_dz: FloatArray
    d2psi1_drr: FloatArray
    d2psi1_dzz: FloatArray
    d2psi1_drz: FloatArray
    u_r: FloatArray
    u_z: FloatArray
    moments: SourceMoments
    boundary_tail: DerivativeTailBounds
    interior_tail: DerivativeTailBounds
    poisson_residual_max: float

    @property
    def tail_bound_available(self) -> bool:
        return self.interior_tail.available

    @property
    def velocity_tail_bound(self) -> float:
        """Bound on the recovered velocity error induced by the elliptic tail.

        ``u^r = -r psi_z`` and ``u^z = 2 psi + r psi_r``, so a bound ``E_0`` on
        ``|psi|`` and ``E_1`` on ``|grad psi|`` gives ``|delta u| <= 2 E_0 +
        R_max E_1`` on the box.
        """
        return 2.0 * self.interior_tail.value + self._radius_max * (
            self.interior_tail.gradient
        )

    _radius_max: float = 0.0

    def as_dict(self) -> dict[str, object]:
        return {
            "moments": self.moments.as_dict(),
            "boundary_tail": self.boundary_tail.as_dict(),
            "interior_tail": self.interior_tail.as_dict(),
            "velocity_tail_bound": self.velocity_tail_bound,
            "poisson_residual_max": self.poisson_residual_max,
        }


def recover_free_space_velocity(
    solver,
    omega1: npt.ArrayLike,
    *,
    interior_radius: float | None = None,
    truncation: str = "monopole",
    support_threshold: float = 0.0,
) -> FreeSpaceRecovery:
    """Solve, differentiate, and attach the derivative tail bounds.

    ``solver`` is a :class:`ns_certificate_lab.whole_space_gate.WholeSpaceEllipticSolver`.
    Derivatives use the audited operators of
    :mod:`ns_certificate_lab.operators`, with the even-at-axis regularity
    condition imposed on the radial derivatives of ``psi1`` — ``psi1`` is a
    smooth even function of ``r`` (E-16), so ``psi1_r(0) = 0`` exactly.

    ``interior_radius`` selects where the interior derivative bound is quoted;
    it defaults to the largest ball centred on the origin that stays one grid
    spacing away from the outer boundary.
    """
    from .operators import derivative_r, derivative_rr, derivative_z, derivative_zz

    grid = solver.grid
    source = grid.validate_field(omega1, name="omega1")
    psi = solver.solve(source)

    psi_r = derivative_r(grid, psi, even_at_axis=True)
    psi_z = derivative_z(grid, psi)
    psi_rr = derivative_rr(grid, psi, even_at_axis=True)
    psi_zz = derivative_zz(grid, psi)
    psi_rz = derivative_z(grid, psi_r)

    r_column = grid.r[:, None]
    u_r = -r_column * psi_z
    u_r[0] = 0.0
    u_z = 2.0 * psi + r_column * psi_r

    moments = source_moments(grid, source, support_threshold=support_threshold)
    boundary_radius = min(float(grid.r[-1]), float(abs(grid.z[0])), float(grid.z[-1]))
    if boundary_radius > moments.support_radius:
        boundary = multipole_tail_bounds(
            moments, evaluation_radius=boundary_radius, truncation=truncation
        )
    else:
        # The source reaches the box boundary: no multipole tail bound exists.
        # Reporting that is the honest outcome; inventing one is not.
        boundary = DerivativeTailBounds.unavailable(
            truncation, boundary_radius - moments.support_radius
        )
    if interior_radius is None:
        interior_radius = boundary_radius - max(grid.dr, grid.dz)
    distance = boundary_radius - float(interior_radius)
    if distance <= 0.0:
        raise ValueError("interior_radius must lie strictly inside the box")
    interior = (
        interior_derivative_tail_bounds(boundary.value, distance_to_boundary=distance)
        if boundary.available
        else DerivativeTailBounds.unavailable("box_truncation", distance)
    )
    return FreeSpaceRecovery(
        psi1=psi,
        dpsi1_dr=psi_r,
        dpsi1_dz=psi_z,
        d2psi1_drr=psi_rr,
        d2psi1_dzz=psi_zz,
        d2psi1_drz=psi_rz,
        u_r=u_r,
        u_z=u_z,
        moments=moments,
        boundary_tail=boundary,
        interior_tail=interior,
        poisson_residual_max=float(np.max(np.abs(solver.residual(psi, source)))),
        _radius_max=float(grid.r[-1]),
    )
