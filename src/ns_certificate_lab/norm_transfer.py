r"""Exact transfer of Sobolev norms between a 3D field and its 2D generators.

An axisymmetric velocity field on ``R^3`` is stored in this repository by two
scalars on the meridional half-plane, ``u1`` and ``psi1``, through

.. math::

   u^\theta = r\,u_1,\qquad
   u^r = -r\,\partial_z\psi_1,\qquad
   u^z = 2\psi_1 + r\,\partial_r\psi_1 .

Theorem NT states that this storage is not merely convenient but an *isometry*
with constant one: every Sobolev norm of the three-dimensional field is a
weighted norm of the two generators, with no unknown constants and no
inequality.  Writing ``<f,g>_mu = int_{r>0} f g r^3 dr dz``,
``L5 = d_rr + (3/r) d_r + d_zz``, ``omega1 = -L5 psi1`` and

    E_k(f) = || L5^(k/2) f ||_mu                 for even k,
    E_k(f) = || grad L5^((k-1)/2) f ||_mu        for odd k,

with ``|grad g|^2 = g_r^2 + g_z^2``, the statement is

    ||u||^2_{Hdot^n(R^3)} = 2 pi [ E_n(u1)^2 + E_{n+1}(psi1)^2 ],
    ||u||^2_{H^n(R^3)}    = 2 pi [ ||u1||^2_mu + E_n(u1)^2
                                   + <psi1,omega1>_mu + E_{n-1}(omega1)^2 ].

The second line is the first line added to the ``L^2`` norm, rewritten so that
both meridional contributions are expressed through ``omega1``; the two
rewrites used are ``E_{n-1}(omega1) = E_{n+1}(psi1)``, which is immediate from
``omega1 = -L5 psi1``, and ``<psi1,omega1>_mu = E_1(psi1)^2``, which is one
integration by parts in the ``r^3`` weight.

Why the ladder closes at all is worth stating, because it is the whole content
of the theorem.  The curl of an axisymmetric field is again axisymmetric, and
in generator variables it acts by

    curl : (u1, psi1) |-> (omega1, u1),   omega1 = -L5 psi1,

since ``omega^r = -r d_z u1``, ``omega^z = 2 u1 + r d_r u1`` reproduce the
meridional map with ``psi1`` replaced by ``u1``, while ``omega^theta = r
omega1``.  For a decaying divergence-free field ``||u||_{Hdot^n}`` is the
``L^2`` norm of the ``n``-th curl, so each step of the ladder advances ``u1``
by one level and ``psi1`` by one level, which is exactly the displayed
formula.  This also explains the off-by-one that is unavoidable in every
argument built on these variables: ``u1`` enters at level ``n`` and ``psi1``
at level ``n+1``.

TWO WEIGHTS, AND WHY THEY MUST NEVER BE CONFUSED
------------------------------------------------
The single easiest error in this repository is to integrate a generator
against the physical measure, or a physical field against the ladder measure.
The two weights are

    mu-weight   ``r^3 dr dz``   for the generators u1, psi1, omega1,
    rho-weight  ``r   dr dz``   for physical fields, with ``dx = 2 pi r dr dz``.

The ``r^3`` weight is not a convention.  ``L5`` is the radial part of the
Laplacian on ``R^5 = R^4 x R`` acting on functions of ``(|x'|, z)``; the
measure on ``R^5`` in those coordinates is ``|S^3| r^3 dr dz``, so ``L5`` is
self-adjoint for ``<.,.>_mu`` and for nothing else.  Every integration by
parts in the ladder, including the one that produces
``<psi1,omega1>_mu = E_1(psi1)^2``, fails outright in the ``r`` weight.
Conversely a physical integral such as ``int |u|^3 dx`` is an integral over
``R^3`` and carries ``2 pi r dr dz``; see :mod:`ns_certificate_lab.critical_l3`,
whose ``axisymmetric_volume_weights`` equals ``2 pi`` times :func:`rho_weights`.
The factor ``r^2`` between the two weights is precisely the ``r`` in
``u^theta = r u1`` and the ``r`` in ``u^r = -r d_z psi1``, squared.

The API names are deliberately unmistakable: :func:`mu_inner`, :func:`mu_norm`
for the ladder weight and :func:`rho_inner`, :func:`rho_norm` for the physical
weight.  Neither name is a default; there is no unqualified ``inner``.

WHAT THIS MODULE IS NOT
-----------------------
Everything here is a second-order finite-difference evaluation on a truncated
half-plane.  It measures the isometry, it does not prove it, and it says
nothing about whole-space truncation.  It also gives no PDE stability
information whatsoever; see :func:`five_dimensional_algebra_warning` and the
section it summarises for the five distinct reasons why the tidy ``H^s``
algebra of the five-dimensional lift does not transfer to the physical
three-dimensional equation.

WHY THE 5D LIFT DOES NOT GIVE 3D STABILITY
------------------------------------------
The ladder above is an identity between norms.  It is tempting to conclude
that one may work with ``u1`` and ``psi1`` as if they solved a nice equation
on ``R^5``, apply the usual ``H^s`` product estimates there, and read the
answer back in ``R^3``.  That fails for five separate reasons, any one of
which is fatal on its own.

First, the natural five-dimensional lift is not incompressible.  The lifted
field has ``div_5 = 2 u^r / r``, which does not vanish, so the five-dimensional
transport is compressible and the antisymmetry ``<u.grad U, U> = 0`` that
makes the transport term free in the energy estimate simply does not hold.
What the failed cancellation leaves behind is not a small remainder: it is
exactly the physical vortex-stretching term, the one term the whole problem is
about.

Second, there is no five-dimensional pressure to work with.  The pressure was
eliminated when the three-dimensional curl was taken to reach the
``(u1, psi1)`` variables, so the lifted system has no Leray projection, no
Riesz-transform structure and no elliptic gain from a pressure equation.  Any
estimate that would ordinarily be closed by moving a term into the pressure
has nowhere to move it.

Third, the two generators do not sit at a common level.  In
``||u||^2_{Hdot^n} = 2 pi [E_n(u1)^2 + E_{n+1}(psi1)^2]`` the stream generator
carries one derivative more than the swirl generator.  An argument that treats
both at a single five-dimensional order is therefore off by one derivative in
the meridional half, and the error is always in the unsafe direction: the
quantity one actually needs is the one with the extra derivative.

Fourth, multiplication by ``r`` -- the map that connects the generators to the
physical components -- is not an ``H^s(R^5)`` multiplier for ``s >= 2``.  The
obstruction is at the axis, not at infinity, because ``r = |x'|`` is a distance
function and is not smooth where it vanishes.  Compact support therefore buys
nothing at all: cutting the field off far away leaves the axis untouched.

Fifth, the embedding thresholds do not match.  ``H^s(R^5)`` embeds in
``L^infinity`` only for ``s > 5/2`` while ``H^s(R^3)`` needs only ``s > 3/2``,
and the reconciliation between the two is carried entirely by the factor ``r``
from the fourth point.  So the apparently harmless move of "working in five
dimensions" costs a full derivative in every embedding, and the cost is paid
back only by the multiplier that has just been ruled out.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid, FloatArray
from .operators import (
    derivative_r,
    derivative_z,
    laplacian_5d_formal,
    recover_velocity,
)

__all__ = [
    "CartesianCrossCheck",
    "HomogeneousSobolevTransfer",
    "SobolevTransfer",
    "apply_ladder_operator",
    "cartesian_sobolev_norm",
    "five_dimensional_algebra_warning",
    "homogeneous_sobolev_norm_from_generators",
    "isometry_defect",
    "ladder_norm",
    "mu_inner",
    "mu_norm",
    "mu_weights",
    "rho_inner",
    "rho_norm",
    "rho_weights",
    "sobolev_norm_from_generators",
]


def _tensor_weights(grid: AxisymmetricGrid, radial_power: int) -> FloatArray:
    """Tensor-product quadrature weights for ``r^p dr dz``.

    The radial rule is trapezoidal and includes the outer wall, so a field that
    has not decayed by ``r_max`` is not silently truncated at half weight.  The
    axial rule is the rectangle rule on a periodic grid, whose endpoint is
    omitted by construction, and the trapezoidal rule otherwise.  These are the
    same conventions as :func:`ns_certificate_lab.critical_l3.axisymmetric_volume_weights`
    so that the two modules cannot disagree about quadrature while appearing to
    disagree about the theorem.
    """

    radial = np.full(grid.nr, grid.dr, dtype=np.float64)
    radial[0] *= 0.5
    radial[-1] *= 0.5
    radial *= grid.r**radial_power

    axial = np.full(grid.nz, grid.dz, dtype=np.float64)
    if not grid.periodic_z:
        axial[0] *= 0.5
        axial[-1] *= 0.5

    weights = radial[:, None] * axial[None, :]
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise FloatingPointError("invalid quadrature weights")
    return weights


def mu_weights(grid: AxisymmetricGrid) -> FloatArray:
    """Quadrature weights for the ladder measure ``r^3 dr dz``.

    This is the measure for which ``L5`` is self-adjoint.  Use it for ``u1``,
    ``psi1``, ``omega1`` and anything obtained from them by ``L5`` or by the
    meridional gradient.  Never use it for a physical field.
    """

    return _tensor_weights(grid, 3)


def rho_weights(grid: AxisymmetricGrid) -> FloatArray:
    """Quadrature weights for the physical measure ``r dr dz``.

    A three-dimensional integral of an axisymmetric integrand is ``2 pi`` times
    the sum against these weights, the ``2 pi`` coming from the azimuthal
    direction.  Use it for physical velocity components, never for a generator.
    """

    return _tensor_weights(grid, 1)


def mu_inner(
    grid: AxisymmetricGrid,
    f: npt.ArrayLike,
    g: npt.ArrayLike,
) -> float:
    """Return ``int f g r^3 dr dz``, the ladder inner product."""

    left = grid.validate_field(f, name="f")
    right = grid.validate_field(g, name="g")
    value = float(np.sum(left * right * mu_weights(grid)))
    if not math.isfinite(value):
        raise FloatingPointError("mu inner product is not finite")
    return value


def mu_norm(grid: AxisymmetricGrid, f: npt.ArrayLike) -> float:
    """Return ``(int f^2 r^3 dr dz)^(1/2)``, the ladder norm."""

    return math.sqrt(max(0.0, mu_inner(grid, f, f)))


def rho_inner(
    grid: AxisymmetricGrid,
    f: npt.ArrayLike,
    g: npt.ArrayLike,
) -> float:
    """Return ``int f g r dr dz``, the physical meridional inner product.

    The physical three-dimensional integral of ``f g`` for axisymmetric ``f``
    and ``g`` is ``2 pi`` times this number.  The factor is not included here
    because the theorem carries it explicitly and hiding it in a helper is how
    a stray ``2 pi`` gets lost.
    """

    left = grid.validate_field(f, name="f")
    right = grid.validate_field(g, name="g")
    value = float(np.sum(left * right * rho_weights(grid)))
    if not math.isfinite(value):
        raise FloatingPointError("rho inner product is not finite")
    return value


def rho_norm(grid: AxisymmetricGrid, f: npt.ArrayLike) -> float:
    """Return ``(int f^2 r dr dz)^(1/2)``, the physical meridional norm."""

    return math.sqrt(max(0.0, rho_inner(grid, f, f)))


def _validated_order(value: object, *, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")
    order = int(value)
    if order < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return order


def apply_ladder_operator(
    grid: AxisymmetricGrid,
    f: npt.ArrayLike,
    power: int,
) -> FloatArray:
    """Return ``L5^power f`` by repeated application of the assembled operator.

    The operator is applied whole, never expanded into ``d_rr``, ``3 d_r / r``
    and ``d_zz`` and recombined afterwards.  For ``power >= 2`` that expansion
    would be arithmetically equivalent only in exact arithmetic on an exactly
    axis-even field: the individual terms of ``L5^m`` carry negative powers of
    ``r`` that diverge separately at the axis and cancel only after assembly.
    Applying the assembled operator repeatedly keeps every intermediate finite.
    """

    repeats = _validated_order(power, name="power", minimum=0)
    result = grid.validate_field(f, name="f")
    for _ in range(repeats):
        result = laplacian_5d_formal(grid, result)
    return result


def ladder_norm(
    grid: AxisymmetricGrid,
    f: npt.ArrayLike,
    k: int,
) -> float:
    """Return ``E_k(f)``: the ``k``-th rung of the ladder for one generator.

    For even ``k`` this is ``|| L5^(k/2) f ||_mu`` and for odd ``k`` it is
    ``|| grad L5^((k-1)/2) f ||_mu`` with ``|grad g|^2 = g_r^2 + g_z^2``.  The
    alternation is not a definitional choice; it is what the curl recursion
    produces, an odd rung being the ``L^2`` norm of a meridional field and an
    even rung the ``L^2`` norm of a swirl field.

    The radial derivative is taken with the axis-even closure, which is the
    correct parity for ``L5^m f``.  The axis row carries zero quadrature weight
    in the ``r^3`` measure and so cannot contribute to the returned number
    directly; it matters only through the interior stencils that read it.
    """

    rung = _validated_order(k, name="k", minimum=0)
    base = apply_ladder_operator(grid, f, rung // 2)
    if rung % 2 == 0:
        return mu_norm(grid, base)
    radial = derivative_r(grid, base, even_at_axis=True)
    axial = derivative_z(grid, base)
    weights = mu_weights(grid)
    value = float(np.sum((radial * radial + axial * axial) * weights))
    if not math.isfinite(value) or value < 0.0:
        raise FloatingPointError("ladder gradient norm is not finite")
    return math.sqrt(value)


@dataclass(frozen=True)
class HomogeneousSobolevTransfer:
    """Pieces of ``||u||^2_{Hdot^n(R^3)} = 2 pi [E_n(u1)^2 + E_{n+1}(psi1)^2]``.

    Both pieces are kept separately because which one dominates is the whole
    diagnostic value of the decomposition: a candidate whose homogeneous norm
    is carried by ``swirl_squared`` is a swirl-driven candidate, and one
    carried by ``stream_squared`` is a meridional-cell candidate, and the two
    fail for different reasons.
    """

    n: int
    swirl_squared: float
    stream_squared: float
    bracket: float
    total_squared: float

    @property
    def norm(self) -> float:
        """The homogeneous norm itself rather than its square."""

        return math.sqrt(max(0.0, self.total_squared))

    @property
    def swirl_fraction(self) -> float:
        """Share of the homogeneous norm held by the swirl generator."""

        if self.bracket <= 0.0:
            return 0.0
        return self.swirl_squared / self.bracket

    def as_dict(self) -> dict[str, object]:
        return {
            "n": int(self.n),
            "swirl_squared": float(self.swirl_squared),
            "stream_squared": float(self.stream_squared),
            "bracket": float(self.bracket),
            "total_squared": float(self.total_squared),
            "norm": float(self.norm),
            "swirl_fraction": float(self.swirl_fraction),
        }


@dataclass(frozen=True)
class SobolevTransfer:
    """Pieces of the inhomogeneous ``H^n`` identity of Theorem NT.

    The four bracket terms are, in order, ``||u1||^2_mu`` (the swirl part of
    the ``L^2`` norm), ``E_n(u1)^2`` (the swirl part of the top order),
    ``<psi1,omega1>_mu`` (the meridional part of the ``L^2`` norm, equal to
    ``E_1(psi1)^2`` after one integration by parts) and ``E_{n-1}(omega1)^2``
    (the meridional part of the top order, equal to ``E_{n+1}(psi1)^2``).  They
    are reported separately so that a caller can see which one dominates
    without recomputing anything.
    """

    n: int
    swirl_l2_squared: float
    swirl_top_squared: float
    stream_l2_squared: float
    stream_top_squared: float
    bracket: float
    total_squared: float

    @property
    def norm(self) -> float:
        return math.sqrt(max(0.0, self.total_squared))

    @property
    def dominant_term(self) -> str:
        """Name of the largest of the four bracket terms."""

        terms = {
            "swirl_l2_squared": self.swirl_l2_squared,
            "swirl_top_squared": self.swirl_top_squared,
            "stream_l2_squared": self.stream_l2_squared,
            "stream_top_squared": self.stream_top_squared,
        }
        return max(terms, key=lambda key: terms[key])

    def as_dict(self) -> dict[str, object]:
        return {
            "n": int(self.n),
            "swirl_l2_squared": float(self.swirl_l2_squared),
            "swirl_top_squared": float(self.swirl_top_squared),
            "stream_l2_squared": float(self.stream_l2_squared),
            "stream_top_squared": float(self.stream_top_squared),
            "bracket": float(self.bracket),
            "total_squared": float(self.total_squared),
            "norm": float(self.norm),
            "dominant_term": self.dominant_term,
        }


def homogeneous_sobolev_norm_from_generators(
    grid: AxisymmetricGrid,
    u1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    n: int,
) -> HomogeneousSobolevTransfer:
    """Evaluate ``||u||^2_{Hdot^n(R^3)} = 2 pi [E_n(u1)^2 + E_{n+1}(psi1)^2]``.

    The homogeneous norm is the one that scales cleanly, so it is the form used
    when comparing candidates across the amplitude-length rescaling; the
    inhomogeneous norm of :func:`sobolev_norm_from_generators` adds the ``L^2``
    part, which does not scale the same way.
    """

    order = _validated_order(n, name="n", minimum=0)
    swirl = ladder_norm(grid, u1, order) ** 2
    stream = ladder_norm(grid, psi1, order + 1) ** 2
    bracket = swirl + stream
    return HomogeneousSobolevTransfer(
        n=order,
        swirl_squared=swirl,
        stream_squared=stream,
        bracket=bracket,
        total_squared=2.0 * math.pi * bracket,
    )


def sobolev_norm_from_generators(
    grid: AxisymmetricGrid,
    u1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    n: int,
) -> SobolevTransfer:
    """Evaluate the inhomogeneous identity of Theorem NT for ``n >= 1``.

    The identity reads

        ||u||^2_{H^n} = 2 pi [ ||u1||^2_mu + E_n(u1)^2
                               + <psi1,omega1>_mu + E_{n-1}(omega1)^2 ]

    with the convention ``||u||^2_{H^n} = ||u||^2_{L^2} + ||u||^2_{Hdot^n}``.

    ``n = 0`` is rejected rather than special-cased.  The formula contains
    ``E_{n-1}(omega1)``, which has no meaning at ``n = 0``, and the honest
    ``n = 0`` statement is that ``H^0`` is ``L^2``, which is precisely
    :func:`homogeneous_sobolev_norm_from_generators` at ``n = 0``.  Silently
    collapsing the two would double count the swirl term and hide the mistake
    behind a factor of two, so the caller is asked to say which one is meant.
    """

    order = _validated_order(n, name="n", minimum=1)
    swirl_field = grid.validate_field(u1, name="u1")
    stream_field = grid.validate_field(psi1, name="psi1")
    omega1 = -laplacian_5d_formal(grid, stream_field)

    swirl_l2 = mu_inner(grid, swirl_field, swirl_field)
    swirl_top = ladder_norm(grid, swirl_field, order) ** 2
    stream_l2 = mu_inner(grid, stream_field, omega1)
    stream_top = ladder_norm(grid, omega1, order - 1) ** 2
    bracket = swirl_l2 + swirl_top + stream_l2 + stream_top
    return SobolevTransfer(
        n=order,
        swirl_l2_squared=swirl_l2,
        swirl_top_squared=swirl_top,
        stream_l2_squared=stream_l2,
        stream_top_squared=stream_top,
        bracket=bracket,
        total_squared=2.0 * math.pi * bracket,
    )


@dataclass(frozen=True)
class CartesianCrossCheck:
    """Independent evaluation of ``||u||^2_{Hdot^n}`` in physical variables.

    ``swirl_part`` collects the contribution of ``u^theta`` and its
    derivatives, ``meridional_part`` that of ``u^r`` and ``u^z``.  The split
    matches the two bracket terms of
    :class:`HomogeneousSobolevTransfer` term by term, which makes a
    disagreement immediately attributable to one generator rather than to the
    pair.
    """

    n: int
    swirl_part: float
    meridional_part: float
    total_squared: float

    @property
    def norm(self) -> float:
        return math.sqrt(max(0.0, self.total_squared))

    def as_dict(self) -> dict[str, object]:
        return {
            "n": int(self.n),
            "swirl_part": float(self.swirl_part),
            "meridional_part": float(self.meridional_part),
            "total_squared": float(self.total_squared),
            "norm": float(self.norm),
        }


def cartesian_sobolev_norm(
    grid: AxisymmetricGrid,
    u1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    n: int,
) -> CartesianCrossCheck:
    r"""Recompute ``||u||^2_{Hdot^n(R^3)}`` without using the ladder at all.

    This is the cross-check that makes Theorem NT testable rather than
    asserted.  It shares no line of arithmetic with
    :func:`homogeneous_sobolev_norm_from_generators`: it reconstructs the
    physical components ``u^r``, ``u^theta`` and ``u^z``, differentiates those,
    and integrates against the physical measure ``2 pi r dr dz``.  In
    particular ``L5`` is never applied and the ``r^3`` weight never appears, so
    a wrong power of ``r`` or a missing ``2 pi`` in the ladder shows up here as
    a gross disagreement rather than as a shared error.

    For ``n = 0`` the integrand is ``|u|^2``.  For ``n = 1`` it is the
    Frobenius norm of the full Cartesian velocity gradient, which for an
    axisymmetric field is the exact pointwise identity

    .. math::

       |\nabla u|^2 = (\partial_r u^r)^2 + (\partial_z u^r)^2
                    + (\partial_r u^\theta)^2 + (\partial_z u^\theta)^2
                    + (\partial_r u^z)^2 + (\partial_z u^z)^2
                    + (u^r/r)^2 + (u^\theta/r)^2 .

    The two apparently singular ratios are evaluated from the generators by the
    exact algebraic identities ``u^theta / r = u1`` and ``u^r / r = -d_z
    psi1``, so nothing is divided by zero and no limit is taken numerically.
    The identity itself holds because the Frobenius norm is invariant under the
    rotation taking the cylindrical orthonormal frame to the Cartesian one; the
    two ratio terms are the frame-rotation terms of the cylindrical gradient
    and they are exactly what a naive ``(d_r, d_z)`` count would omit.

    WHAT IS CROSS-CHECKED.  The constant ``2 pi``; the ``r^3`` versus ``r``
    weight; the boundary term that makes the ``4 r psi1^2 + 4 r^2 psi1
    d_r psi1`` cross term in ``|u^z|^2`` integrate to zero; and, at ``n = 1``,
    the identification of ``E_2(psi1)`` with the swirl component of the
    vorticity and of ``E_1(u1)`` with its meridional components.

    WHAT IS NOT CROSS-CHECKED.  Nothing above ``n = 1``: for ``n >= 2`` the
    ladder rests on the curl recursion, and only its internal consistency is
    testable here.  Nor is the whole-space statement: both sides are truncated
    to the same finite half-plane, so a field that has not decayed by ``r_max``
    makes both sides equally wrong.  Nor is smoothness: agreement of two
    quadratures says nothing about whether the stored arrays extend to a smooth
    axisymmetric field, which is the business of
    :func:`ns_certificate_lab.axis.check_axis_regularity`.
    """

    order = _validated_order(n, name="n", minimum=0)
    if order > 1:
        raise ValueError(
            "the Cartesian cross-check is implemented for n = 0 and n = 1 only; "
            "higher rungs would need repeated Cartesian curls and would test "
            "the quadrature rather than the theorem"
        )
    swirl_rate = grid.validate_field(u1, name="u1")
    stream = grid.validate_field(psi1, name="psi1")
    u_r, u_z = recover_velocity(grid, stream)
    u_theta = grid.r[:, None] * swirl_rate
    two_pi = 2.0 * math.pi

    if order == 0:
        swirl_part = two_pi * rho_inner(grid, u_theta, u_theta)
        meridional_part = two_pi * (
            rho_inner(grid, u_r, u_r) + rho_inner(grid, u_z, u_z)
        )
        return CartesianCrossCheck(
            n=order,
            swirl_part=swirl_part,
            meridional_part=meridional_part,
            total_squared=swirl_part + meridional_part,
        )

    # The exact axis-regular forms of the two frame-rotation ratios.
    theta_ratio = swirl_rate
    radial_ratio = -derivative_z(grid, stream)

    swirl_density = (
        derivative_r(grid, u_theta, even_at_axis=False) ** 2
        + derivative_z(grid, u_theta) ** 2
        + theta_ratio**2
    )
    meridional_density = (
        derivative_r(grid, u_r, even_at_axis=False) ** 2
        + derivative_z(grid, u_r) ** 2
        + derivative_r(grid, u_z, even_at_axis=True) ** 2
        + derivative_z(grid, u_z) ** 2
        + radial_ratio**2
    )
    weights = rho_weights(grid)
    swirl_part = two_pi * float(np.sum(swirl_density * weights))
    meridional_part = two_pi * float(np.sum(meridional_density * weights))
    total = swirl_part + meridional_part
    if not math.isfinite(total):
        raise FloatingPointError("Cartesian cross-check is not finite")
    return CartesianCrossCheck(
        n=order,
        swirl_part=swirl_part,
        meridional_part=meridional_part,
        total_squared=total,
    )


def isometry_defect(
    grid: AxisymmetricGrid,
    u1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    n: int,
) -> float:
    """Relative disagreement between the ladder and the physical evaluation.

    Both sides are second-order accurate and use different stencils, so the
    returned number is expected to fall like the square of the mesh spacing.
    It is a discretisation diagnostic; a value that stops falling under
    refinement indicates a genuine error in one of the two routes, whereas a
    value that falls at second order indicates only that finite differences are
    finite differences.
    """

    ladder = homogeneous_sobolev_norm_from_generators(grid, u1, psi1, n)
    physical = cartesian_sobolev_norm(grid, u1, psi1, n)
    scale = max(abs(ladder.total_squared), abs(physical.total_squared))
    if scale <= 0.0:
        return 0.0
    return abs(ladder.total_squared - physical.total_squared) / scale


def five_dimensional_algebra_warning() -> tuple[str, ...]:
    """Return the five reasons the 5D lift gives no 3D PDE stability.

    The reasons are returned as data, not merely written in prose, so that a
    report generator can quote them and a test can assert that all five are
    still present.  Losing one of them silently is how a norm identity gets
    mistaken for an existence theorem.
    """

    return (
        "The natural five-dimensional lift is not divergence free: div_5 = "
        "2 u^r / r.  The five-dimensional transport is therefore compressible "
        "and the antisymmetry <u.grad U, U> = 0 fails.  What the failed "
        "cancellation leaves behind is exactly the physical vortex-stretching "
        "term.",
        "There is no five-dimensional pressure.  It was eliminated by the "
        "three-dimensional curl that produced the generator variables, so the "
        "lifted system has no Leray projection and no elliptic gain to absorb "
        "a term into.",
        "The generators do not sit at a common level: u1 enters at order n "
        "and psi1 at order n+1.  Any argument carrying both at one "
        "five-dimensional order is off by one derivative in the meridional "
        "half, always in the unsafe direction.",
        "Multiplication by r is not an H^s(R^5) multiplier for s >= 2, and the "
        "obstruction is the axis rather than infinity, so compact support buys "
        "nothing.",
        "H^s(R^5) embeds in L^infinity only for s > 5/2 while H^s(R^3) needs "
        "only s > 3/2, and the whole reconciliation is carried by the factor "
        "r, that is, by the multiplier just ruled out.",
    )
