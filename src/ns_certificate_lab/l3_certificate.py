r"""A rigorous interval certificate for a positive initial generation rate.

The object certified is

.. math::  J(u_0) = \underbrace{3\int p\,\nabla\!\cdot\!(|u|u)\,dx}_{P}
     \underbrace{-\,3\nu\int|u|\bigl(|\nabla u|^2+|\nabla|u||^2\bigr)dx}_{V},

with the advective contribution *exactly* zero by the identity
``|u|\,u\cdot(u\cdot\nabla)u = \tfrac13 u\cdot\nabla(|u|^3)``.

What makes a rigorous bound possible here
-----------------------------------------
The datum is not a grid array.  It is a finite linear combination of explicit
generators, so every quantity in ``V`` — the velocity, all nine components of
``\nabla u``, ``|u|`` and ``\nabla|u|`` — is an explicit composition of
``\exp``, powers and rational functions of ``(r^2, z)``.  Interval-evaluating
those compositions on a **cell** ``[r_i,r_{i+1}]\times[z_j,z_{j+1}]`` rather
than at a node gives a rigorous enclosure of the integrand over the whole cell,
so the quadrature error is *enclosed* rather than estimated.  No Lipschitz
hypothesis and no inflation factor are needed anywhere in ``V``.

The one thing that is not explicit is the pressure, which solves a Poisson
equation.  The certificate therefore does **not** assert ``P`` directly.  It
proves the conditional statement

.. math::

   \|p - p_h\|_{L^\infty(\Omega)} \le \tau_\ast
   \quad\Longrightarrow\quad J(u_0) > 0 ,

where ``p_h`` is the computed pressure, ``\Omega`` is the box, and

.. math::  \tau_\ast = \frac{P_h^{\rm lo} + V^{\rm lo}}{3\int|g|\,dx},
   \qquad g = \nabla\!\cdot\!(|u|u) ,

with ``P_h^{\rm lo}`` a rigorous lower bound on ``3\int p_h g\,dx`` and
``V^{\rm lo}`` a rigorous lower bound on ``V``.  Both of those, and
``\int|g|dx``, are rigorous: ``g`` is explicit, and ``p_h`` is enclosed by the
**multilinear interpolant** of its nodal values, whose range over a cell is
*exactly* the corner hull (see :func:`bilinear_cell_range`).

``\tau_\ast`` is an unconditional, checkable number.  Turning the implication
into an unconditional ``J>0`` needs a bound on ``\|p-p_h\|_\infty``, which is
recorded as the single named hypothesis ``P1`` and is **not proved here**.  It
is the pressure-side analogue of ``HS-5``.

Arithmetic
----------
Exact rational intervals throughout
(:class:`~ns_certificate_lab.snapshot_certificate.Interval`), with outward
rounding.  ``\exp`` is the only transcendental involved and it is bounded by
:func:`exp_interval`, which uses a truncated Taylor series with a rigorous
geometric remainder — never a floating-point ``math.exp``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "GenerationCertificate",
    "PRESSURE_HYPOTHESIS",
    "bilinear_cell_range",
    "build_generation_certificate",
    "exp_interval",
    "reciprocal_interval",
    "sqrt_interval",
    "verify_generation_certificate",
]

PRESSURE_HYPOTHESIS = (
    "P1 (pressure discretisation): ||p - p_h||_{L^inf(Omega)} <= tau_star, with "
    "p the free-space solution of -Laplacian p = d_i u_j d_j u_i and p_h the "
    "computed pressure.  Unproved here: it needs an a posteriori bound on a "
    "Poisson solve, which is the pressure-side analogue of HS-5."
)


# --------------------------------------------------------------------------- #
# rational interval elementary functions                                       #
# --------------------------------------------------------------------------- #


def reciprocal_interval(value: Interval) -> Interval:
    """``1/x`` for an interval that does not contain zero."""
    if value.contains_zero:
        raise ZeroDivisionError("cannot invert an interval containing zero")
    return Interval(min(1 / value.lower, 1 / value.upper),
                    max(1 / value.lower, 1 / value.upper))


def exp_interval(value: Interval, *, terms: int = 40) -> Interval:
    r"""A rigorous rational enclosure of ``\exp`` on an interval.

    ``\exp`` is increasing, so the range on ``[a,b]`` is ``[e^a, e^b]`` and only
    a lower bound at ``a`` and an upper bound at ``b`` are needed.

    For ``x \le 0`` write ``e^x = 1/e^{-x}`` with ``-x \ge 0``, so it suffices to
    bound ``e^y`` for ``y \ge 0``.  Truncating ``\sum y^k/k!`` after ``n`` terms
    leaves a remainder bounded, when ``y \le n/2``, by

    .. math::  \sum_{k>n}\frac{y^k}{k!} \le \frac{y^n}{n!}\cdot
               \frac{y/(n+1)}{1-y/(n+1)} ,

    a geometric bound.  Arguments larger than ``n/2`` are reduced by
    ``e^y = (e^{y/2^m})^{2^m}`` until they are small, which keeps the remainder
    tiny and the rationals of bounded size.
    """
    if terms < 4:
        raise ValueError("at least four Taylor terms are required")

    guard = Fraction(2) ** 96

    def bounds_nonneg(y: Fraction) -> tuple[Fraction, Fraction]:
        """Rational bounds on e^y for y >= 0.

        Two things keep the rationals small.  The argument is rounded to a
        dyadic bracket first -- near the edge of a bump support 1/(1-sigma)
        is a fraction with an enormous denominator, and carrying it exactly
        through thirty-two Taylor terms squares that denominator thirty-two
        times.  And a very large argument is cut off outright: ln 2 > 0.6931,
        so y >= 700 gives e^{-y} <= 2^{-1000}, which is a rigorous bound
        and not an approximation.
        """
        if y < 0:
            raise ValueError('bounds_nonneg needs a nonnegative argument')
        if y > 700:
            # e^y >= 2^1000 and, at the reciprocal, e^{-y} <= 2^{-1000}.
            return Fraction(2) ** 1000, Fraction(2) ** 100000
        low_arg = Fraction(math.floor(y * guard), 1) / guard
        high_arg = Fraction(math.ceil(y * guard), 1) / guard

        def series(argument: Fraction, *, up: bool) -> Fraction:
            value = argument
            squarings = 0
            while value > Fraction(1, 2):
                value /= 2
                squarings += 1
            partial = Fraction(0)
            term = Fraction(1)
            for k in range(terms):
                if k:
                    term = term * value / k
                partial += term
            ratio = value / (terms + 1)
            remainder = term * value / terms * ratio / (1 - ratio)
            result = partial + remainder if up else partial
            for _ in range(squarings):
                result = result * result
                result = (
                    Fraction(math.ceil(result * guard), 1) / guard
                    if up
                    else Fraction(math.floor(result * guard), 1) / guard
                )
                if result < 0:
                    result = Fraction(0)
            return result

        return series(low_arg, up=False), series(high_arg, up=True)

    def lower(x: Fraction) -> Fraction:
        if x >= 0:
            return bounds_nonneg(x)[0]
        return 1 / bounds_nonneg(-x)[1]

    def upper(x: Fraction) -> Fraction:
        if x >= 0:
            return bounds_nonneg(x)[1]
        return 1 / bounds_nonneg(-x)[0]

    return Interval(lower(value.lower), upper(value.upper))


def sqrt_interval(value: Interval, *, bits: int = DEFAULT_PRECISION_BITS) -> Interval:
    """A rigorous rational enclosure of ``sqrt`` on a nonnegative interval."""
    if value.lower < 0:
        raise ValueError("cannot take the square root of a negative interval")

    def root(x: Fraction, *, up: bool) -> Fraction:
        if x == 0:
            return Fraction(0)
        scale = Fraction(2) ** bits
        # Round the argument outward to a dyadic first.  Without it the
        # denominator of x can be 2^(k*bits) after a product chain, and the
        # integer handed to isqrt grows without bound -- which is what made an
        # unrounded version thousands of times slower than the rest of the
        # enclosure put together.
        dyadic = (
            Fraction(math.ceil(x * scale), 1) / scale
            if up
            else Fraction(math.floor(x * scale), 1) / scale
        )
        if dyadic < 0:
            dyadic = Fraction(0)
        target = int(dyadic * scale * scale)
        integer = math.isqrt(target)
        if up and integer * integer < target:
            integer += 1
        return Fraction(integer) / scale

    return Interval(root(value.lower, up=False), root(value.upper, up=True))


def bilinear_cell_range(
    corner_00: Fraction, corner_10: Fraction, corner_01: Fraction, corner_11: Fraction
) -> Interval:
    r"""The **exact** range of a bilinear interpolant over its cell.

    On ``[0,1]^2`` the bilinear interpolant is

    .. math::  B(\xi,\eta) = \sum_{a,b\in\{0,1\}}
       \lambda_a(\xi)\lambda_b(\eta)\,v_{ab},\qquad
       \lambda_0(t)=1-t,\ \lambda_1(t)=t ,

    and the four weights are nonnegative and sum to one.  ``B`` is therefore a
    **convex combination** of the corner values at every point of the cell, so
    its range is exactly the hull of those four values — no inflation, no
    Lipschitz constant, no hypothesis.

    This is the reason the certificate interpolates multilinearly rather than
    with anything higher order: a bicubic would need a Bernstein conversion to
    get an enclosure, and a piecewise-constant reconstruction would need a
    variation bound.  The multilinear interpolant is the one whose cell range is
    free.
    """
    values = (corner_00, corner_10, corner_01, corner_11)
    return Interval(min(values), max(values))


# --------------------------------------------------------------------------- #
# interval evaluation of the generators on a cell                              #
# --------------------------------------------------------------------------- #


def _chi_interval(argument: Interval, *, terms: int) -> Interval:
    r"""``\chi(\sigma) = \exp(-1/(1-\sigma))`` on ``[0,1)``, enclosed exactly.

    ``\chi`` is **decreasing**: ``\chi' = -\chi/(1-\sigma)^2 < 0``.  So the range
    on ``[a,b]`` is ``[\chi(b), \chi(a)]`` and monotonicity gives the sharp
    enclosure with no dependency-problem widening.  Outside ``[0,1)`` the value
    is zero by definition, and an interval straddling ``1`` is enclosed by
    ``[0, \chi(a)]``.
    """
    low, high = argument.lower, argument.upper
    if low >= 1:
        return Interval(Fraction(0), Fraction(0))
    left = max(low, Fraction(0))
    if high >= 1:
        upper = exp_interval(
            Interval(-1 / (1 - left), -1 / (1 - left)), terms=terms
        ).upper
        return Interval(Fraction(0), upper)
    lower_value = exp_interval(
        Interval(-1 / (1 - high), -1 / (1 - high)), terms=terms
    ).lower
    upper_value = exp_interval(
        Interval(-1 / (1 - left), -1 / (1 - left)), terms=terms
    ).upper
    return Interval(lower_value, upper_value)


# --------------------------------------------------------------------------- #
# the certificate                                                              #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class GenerationCertificate:
    """A machine-checkable enclosure of the initial generation rate."""

    payload: dict[str, object]

    def as_dict(self) -> dict[str, object]:
        return self.payload


def _cell_hull_from_nodes(values: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Corner-hull min and max of a nodal array, one entry per cell."""
    stack = np.stack(
        (values[:-1, :-1], values[1:, :-1], values[:-1, 1:], values[1:, 1:])
    )
    return stack.min(axis=0), stack.max(axis=0)


def build_generation_certificate(
    family,
    grid: AxisymmetricGrid,
    *,
    viscosity: float,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    exp_terms: int = 40,
    subdivisions: int = 1,
) -> GenerationCertificate:
    r"""Build the certificate for one mixed family on one grid.

    The rigorous pieces are, cell by cell in exact rational arithmetic:

    * an enclosure of ``|u|``, ``|\nabla u|^2`` and ``|\nabla|u||^2`` over the
      **whole cell**, from interval evaluation of the analytic generators;
    * hence an upper bound on ``|V|`` and a lower bound on ``V``;
    * an upper bound on ``\int|g|dx`` with ``g = \nabla\cdot(|u|u)``;
    * a lower bound on ``3\int p_h g\,dx`` using the multilinear enclosure of
      ``p_h`` and the interval enclosure of ``g``;
    * the threshold ``\tau_\ast`` above which the pressure error would destroy
      the conclusion.

    ``subdivisions`` splits each grid cell before the interval evaluation.  The
    dependency problem makes a single interval evaluation over a wide cell
    pessimistic; subdividing narrows it at linear cost and is the standard
    remedy.
    """
    from .l3_generation import AxisymmetricPressureSolver, l3_generation_rate

    if viscosity <= 0.0:
        raise ValueError("a Clay candidate requires a fixed positive viscosity")
    if grid.periodic_z:
        raise ValueError("the generation certificate targets non-periodic z grids")
    if subdivisions < 1:
        raise ValueError("subdivisions must be at least one")

    field = family.field(grid)
    solver = AxisymmetricPressureSolver.build(grid)
    rate = l3_generation_rate(field, viscosity=viscosity, solver=solver)

    gradient = field.gradient()
    speed = field.speed
    gradient_square = sum(v * v for v in gradient.values())
    square = speed**2
    from .operators import derivative_r, derivative_z

    g_r = 0.5 * derivative_r(grid, square, even_at_axis=True)
    g_z = 0.5 * derivative_z(grid, square)
    safe = np.where(speed > 0.0, speed, 1.0)
    grad_speed = np.minimum(
        np.where(speed > 0.0, (g_r**2 + g_z**2) / safe, 0.0), speed * gradient_square
    )
    viscous_integrand = speed * gradient_square + grad_speed
    # g = div(|u| u) = u . grad |u|, evaluated as (u^r G_r + u^z G_z)/|u|, which
    # is bounded by |u| |grad u| and is defined to be zero on the zero set.
    flux = np.minimum(
        np.abs(
            np.where(speed > 0.0, (field.u_r * g_r + field.u_z * g_z) / safe, 0.0)
        ),
        speed * np.sqrt(gradient_square),
    )

    pressure_solver = solver
    trace = pressure_solver.exterior_trace(field.u_r, field.u_theta, field.u_z)
    pressure = pressure_solver.solve(field.pressure_source(gradient), trace=trace)

    # -- rigorous cell enclosures ------------------------------------------- #
    #
    # The analytic quantities are enclosed by INTERVAL EVALUATION of the
    # generators over each cell, never by a hull of nodal values: a nodal hull
    # is not an enclosure of a general function, and inferring one from stored
    # samples is exactly the move the certificate discipline forbids.  The
    # computed pressure is different -- it *is* a grid function, and the
    # multilinear interpolant of a grid function has its cell range exactly
    # equal to the corner hull, so there the hull is a theorem.
    from .interval_generators import cell_enclosure

    r_nodes = [Fraction(float(v)) for v in grid.r]
    z_nodes = [Fraction(float(v)) for v in grid.z]
    two_pi_upper = Fraction(6283185307179587, 10**15)  # > 2 pi
    p_lo, p_hi = _cell_hull_from_nodes(pressure)

    viscous_upper = Fraction(0)
    flux_upper = Fraction(0)
    pressure_lower = Fraction(0)
    widest_flux = Fraction(0)
    divergence_excursion = Fraction(0)
    cells = 0
    step_r = Fraction(float(grid.dr)) / subdivisions
    step_z = Fraction(float(grid.dz)) / subdivisions

    for i in range(grid.nr - 1):
        for j in range(grid.nz - 1):
            pressure_box = Interval(
                Fraction(float(p_lo[i, j])), Fraction(float(p_hi[i, j]))
            )
            for a in range(subdivisions):
                r_box = Interval(
                    r_nodes[i] + a * step_r, r_nodes[i] + (a + 1) * step_r
                )
                # 2 pi r dr dz, bounded above by the largest radius of the box.
                weight = two_pi_upper * r_box.upper * step_r * step_z
                for b in range(subdivisions):
                    z_box = Interval(
                        z_nodes[j] + b * step_z, z_nodes[j] + (b + 1) * step_z
                    )
                    enclosure = cell_enclosure(
                        family, r_box, z_box, terms=exp_terms, bits=precision_bits
                    )
                    cells += 1
                    viscous_upper += weight * max(
                        enclosure.viscous_integrand.upper, Fraction(0)
                    )
                    flux_upper += weight * enclosure.flux_magnitude.upper
                    products = (
                        pressure_box.lower * enclosure.flux.lower,
                        pressure_box.lower * enclosure.flux.upper,
                        pressure_box.upper * enclosure.flux.lower,
                        pressure_box.upper * enclosure.flux.upper,
                    )
                    pressure_lower += weight * min(products)
                    width = enclosure.flux.upper - enclosure.flux.lower
                    widest_flux = max(widest_flux, width)
                    divergence_excursion = max(
                        divergence_excursion, enclosure.divergence.magnitude
                    )
        viscous_upper = Interval(viscous_upper, viscous_upper).round_outward(
            precision_bits
        ).upper
        flux_upper = Interval(flux_upper, flux_upper).round_outward(
            precision_bits
        ).upper
        pressure_lower = Interval(pressure_lower, pressure_lower).round_outward(
            precision_bits
        ).lower

    nu = Fraction(float(viscosity))
    viscous_bound = -3 * nu * viscous_upper          # a lower bound on V
    pressure_bound = 3 * pressure_lower              # a lower bound on 3 int p_h g
    flux_integral = 3 * flux_upper                   # an upper bound on 3 int |g|
    margin = pressure_bound + viscous_bound
    threshold = margin / flux_integral if flux_integral > 0 else Fraction(0)

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate": "l3_generation_rate",
        "family": family.as_dict(),
        "grid": {
            "nr": grid.nr, "nz": grid.nz,
            "r_max": float(grid.r[-1]),
            "z_min": float(grid.z[0]), "z_max": float(grid.z[-1]),
            "subdivisions": subdivisions,
        },
        "viscosity": float(viscosity),
        "rigorous": {
            "viscous_lower_bound": str(viscous_bound),
            "viscous_integral_upper_bound": str(viscous_upper),
            "pressure_lower_bound_discrete": str(pressure_bound),
            "flux_integral_upper_bound": str(flux_integral),
            "margin": str(margin),
            "pressure_error_threshold": str(threshold),
            "margin_is_positive": bool(margin > 0),
            "cells_enclosed": cells,
            "widest_flux_enclosure": str(widest_flux),
            "divergence_excursion": str(divergence_excursion),
        },
        "hypotheses": {
            "P1_pressure_discretisation": {
                "statement": PRESSURE_HYPOTHESIS,
                "threshold": str(threshold),
                "proved": False,
            }
        },
        "floating_point_reference": rate.as_dict(),
        "claims": [
            "Rigorous and unconditional: the stated lower bound on the viscous "
            "contribution, the upper bound on the flux integral, and the lower "
            "bound on the discrete pressure contribution, all in exact rational "
            "arithmetic with cell-interior enclosures.",
            "Conditional on P1: if the computed pressure is within the stated "
            "threshold of the true free-space pressure in L^infinity on the box, "
            "then J(u_0) > 0.",
            "The advective contribution is exactly zero by an algebraic "
            "identity, not by cancellation of computed numbers.",
            "Nothing here is a singularity result and nothing here bears on the "
            "Clay problem.  A positive initial generation rate of the critical "
            "norm is a screening criterion, not evidence of blow-up.",
        ],
    }
    return GenerationCertificate(payload)


def verify_generation_certificate(payload: dict[str, object]) -> dict[str, object]:
    """Independently re-check a generation payload from its own contents."""
    try:
        return _verify(payload)
    except Exception as error:  # noqa: BLE001 -- adversarial input is a failure
        return {"verified": False, "failures": [f"payload rejected: {error}"]}


def _verify(payload: dict[str, object]) -> dict[str, object]:
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(payload.get("certificate") == "l3_generation_rate", "wrong certificate kind")
    require(payload.get("schema_version") == 1, "unsupported schema version")
    require(float(payload["viscosity"]) > 0.0, "viscosity must be positive")

    block = payload.get("rigorous")
    require(isinstance(block, dict), "missing rigorous block")
    if isinstance(block, dict):
        viscous = Fraction(str(block["viscous_lower_bound"]))
        pressure = Fraction(str(block["pressure_lower_bound_discrete"]))
        flux = Fraction(str(block["flux_integral_upper_bound"]))
        margin = Fraction(str(block["margin"]))
        threshold = Fraction(str(block["pressure_error_threshold"]))
        require(viscous <= 0, "the viscous contribution must be non-positive")
        require(flux >= 0, "the flux integral bound must be nonnegative")
        require(margin == pressure + viscous, "the margin does not equal P_h + V")
        require(
            bool(block["margin_is_positive"]) == (margin > 0),
            "the margin flag contradicts the margin",
        )
        if flux > 0:
            require(
                threshold * flux == margin,
                "the threshold does not equal the margin divided by 3 int |g|",
            )
        require(
            (margin > 0) == (threshold > 0),
            "a positive threshold requires a positive margin",
        )

    hypotheses = payload.get("hypotheses")
    require(isinstance(hypotheses, dict), "missing hypotheses")
    if isinstance(hypotheses, dict):
        require(
            set(hypotheses) == {"P1_pressure_discretisation"},
            "the payload must carry exactly the one named hypothesis",
        )
        entry = hypotheses["P1_pressure_discretisation"]
        require(entry.get("proved") is False, "P1 claims to be proved; it is not")
        require(
            entry.get("statement") == PRESSURE_HYPOTHESIS,
            "the P1 statement does not match the audited text",
        )
        if isinstance(block, dict):
            require(
                Fraction(str(entry["threshold"]))
                == Fraction(str(block["pressure_error_threshold"])),
                "the hypothesis threshold disagrees with the rigorous block",
            )

    claims = payload.get("claims")
    require(
        isinstance(claims, list) and any("Clay" in str(c) for c in claims),
        "the payload must record that it says nothing about the Clay problem",
    )
    return {"verified": not failures, "failures": failures}
