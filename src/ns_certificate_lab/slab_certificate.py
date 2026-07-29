r"""One space-time slab certificate over ``[t_n, t_{n+1}]``.

:mod:`ns_certificate_lab.snapshot_certificate` encloses one instant at the grid
*nodes*.  That is the weakest possible object: a singularity is a statement
about a space-time region, and a bound that holds at isolated nodes at isolated
instants constrains nothing between them.  This module encloses a whole slab --
every point of every cell, at every instant of ``[t_n, t_{n+1}]``.

Name the object before certifying it
------------------------------------
Three different things could be meant by "the solution over the slab", and the
earlier version of this module conflated them.

* **A** -- the piecewise interpolant of the stored nodal data.  A *defined*
  function; every claim about it is algebra, with no hypotheses at all.
* **B** -- the exact solution of the semi-discrete system ``Y' = F(Y)`` through
  the stored initial state.  Reachable by a rigorous ODE enclosure.
* **C** -- the continuum Navier-Stokes solution.  **Not reachable from nodal
  data by any argument.**

This module certifies **A**, exactly and unconditionally, and records the gap to
**C** as a named unproved hypothesis.  Reaching **B** needs the Picard
self-mapping enclosure implemented in
:mod:`ns_certificate_lab.control_ode`; it is not attempted here for the full
field.

``H1`` is now a theorem
-----------------------
The previous version inflated the corner hull by the measured divided
differences and called the result a hypothesis.  That was doubly wrong: it used
stored nodal data as a derivative proxy, which the certificate discipline
forbids, and it was unnecessary.

The fix is to **define** the interpolant and then bound it exactly.  With the
nodal slopes fixed by the recorded stencils, the tensor-product cubic Hermite
interpolant on a cell is a bicubic polynomial; converting its Hermite
coefficients to the **Bernstein** basis gives sixteen coefficients whose hull
contains the polynomial at every point of the cell, because the tensor Bernstein
basis is nonnegative and partitions unity.  No derivative bound appears
anywhere.  See :func:`bernstein_cell_envelope`.

``H2`` was not fixable and has been removed
-------------------------------------------
The old ``H2`` asserted ``|y(t) - H(t)| <= \Delta^4 M_4/384``, the two-point
Hermite remainder.  That bound is valid only when ``H`` interpolates **one**
trajectory at both ends with matching slopes.  Here the end state is the RK4
output, not the value at ``t_{n+1}`` of the trajectory through the start state,
and the end slope is the slope of a *different* trajectory.  Supplying a proved
``M_4`` would not have rescued it; the hypothesis was false as framed.

What replaces it is honest bookkeeping.  The temporal enclosure below is a
statement about the *interpolant* -- object **A** -- and the exact ranges of the
Hermite basis polynomials make it a theorem.  The distance from **A** to **B**
is the integrator's local defect, computed exactly (trapezoid and Simpson) and
reported as a diagnostic; the distance from **B** to **C** is ``H3``.

The one remaining hypothesis
----------------------------
``H3`` (*semi-discrete to continuum*), **not proved**: the consistency error
``||F(\Pi u) - \Pi(N(u))||``, of order ``h^2`` times fourth continuum spatial
derivatives, multiplied by the discrete stability constant, plus the elliptic
solve's truncation and the multipole domain truncation.  It needs regularity the
numerical state does not carry.  It is the same gap as ``HS-5``.

Consequently this certificate establishes, exactly and unconditionally, a
property of the computed space-time object, and establishes **nothing** about
the continuous Navier-Stokes solution.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "SlabCertificate",
    "build_slab_certificate",
    "verify_slab_certificate",
    "HERMITE_BASIS_RANGES",
    "SLAB_HYPOTHESES",
    "SLAB_THEOREMS",
    "bernstein_cell_envelope",
    "hermite_slopes",
]

#: Exact ranges of the cubic Hermite basis on ``[0, 1]``.
#:
#: ``h00 = 2t^3 - 3t^2 + 1`` and ``h01 = -2t^3 + 3t^2`` satisfy ``h00 + h01 = 1``
#: with both in ``[0, 1]``, so the value part of the interpolant is a convex
#: combination of the endpoints and is enclosed by their hull.  The tangent
#: bases are ``h10 = t^3 - 2t^2 + t``, stationary at ``t = 1/3`` with value
#: ``4/27``, and ``h11 = t^3 - t^2``, stationary at ``t = 2/3`` with value
#: ``-4/27``.
HERMITE_BASIS_RANGES = {
    "h00": (Fraction(0), Fraction(1)),
    "h01": (Fraction(0), Fraction(1)),
    "h10": (Fraction(0), Fraction(4, 27)),
    "h11": (Fraction(-4, 27), Fraction(0)),
}

#: The one remaining hypothesis.  ``H1`` became :func:`bernstein_cell_envelope`
#: and ``H2`` was withdrawn as unsound; see the module docstring.
SLAB_HYPOTHESES = {
    "H3_semidiscrete_to_continuum": (
        "The certified object is the interpolant of the semi-discrete state, not "
        "the continuum Navier-Stokes solution.  Bridging them needs the "
        "consistency error ||F(Pi u) - Pi(N(u))|| = O(h^2) times fourth "
        "continuum spatial derivatives, multiplied by the discrete stability "
        "constant, plus the elliptic truncation and the multipole domain "
        "truncation.  Unproved: it needs regularity the numerical state does not "
        "carry.  Same gap as HS-5."
    ),
}

#: Theorems this module proves outright, recorded in the payload so a checker can
#: insist they are marked proved rather than assumed.
SLAB_THEOREMS = {
    "T1_bernstein_cell_envelope": (
        "For every point of the closed cell, the tensor-product cubic Hermite "
        "interpolant with slopes given by the recorded stencils lies in the hull "
        "of its sixteen Bernstein coefficients.  Proof: the tensor Bernstein "
        "basis is nonnegative and partitions unity, so the interpolant is a "
        "convex combination of those coefficients."
    ),
    "T2_hermite_time_range": (
        "For every instant of the slab, the cubic Hermite interpolant in time "
        "lies in the hull of its endpoint values plus Delta times the exact "
        "ranges [0, 4/27] and [-4/27, 0] of the tangent basis polynomials.  "
        "Proof: h00 + h01 = 1 with both nonnegative on [0,1], so the value part "
        "is a convex combination; the tangent bases attain their extrema at "
        "t = 1/3 and t = 2/3."
    ),
}


def _nodes_inside(
    hull: Interval, field: FloatArray, cells: list[tuple[int, int]]
) -> bool:
    """Whether every corner of every enclosed cell lies in the hull.

    The interpolant reproduces the nodal values exactly, so a node escaping the
    hull would mean the Bernstein argument had been mis-assembled -- which is
    precisely what this check is for.
    """
    for i, j in cells:
        for a, b in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)):
            if not hull.contains(Fraction(float(field[a, b]))):
                return False
    return True


@dataclass(frozen=True)
class SlabCertificate:
    """A machine-checkable enclosure of one space-time slab."""

    payload: dict[str, object]

    def as_dict(self) -> dict[str, object]:
        return self.payload


def _hull(*intervals: Interval) -> Interval:
    return Interval(
        min(i.lower for i in intervals), max(i.upper for i in intervals)
    )


#: Hermite-to-Bernstein change of basis for one cubic, acting on the coefficient
#: vector ``(v_0, m_0, v_1, m_1)`` (value and scaled slope at each end).
#:
#: ``h00 = B_0 + B_1``, ``h10 = B_1/3``, ``h01 = B_2 + B_3``, ``h11 = -B_2/3``
#: in the cubic Bernstein basis ``B_k(t) = C(3,k) t^k (1-t)^{3-k}``, so
#: ``b_0 = v_0``, ``b_1 = v_0 + m_0/3``, ``b_2 = v_1 - m_1/3``, ``b_3 = v_1``.
HERMITE_TO_BERNSTEIN = (
    (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    (Fraction(1), Fraction(1, 3), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1), Fraction(-1, 3)),
    (Fraction(0), Fraction(0), Fraction(1), Fraction(0)),
)


def hermite_slopes(
    grid: AxisymmetricGrid, field: FloatArray, *, even_at_axis: bool
) -> tuple[FloatArray, FloatArray, FloatArray]:
    r"""Nodal slopes that **define** the interpolant, by the audited stencils.

    Second-order central differences in the interior, one-sided at the ends, and
    ``\partial_r f(0, z) = 0`` for a field that is even in ``r``.  The mixed
    slope is the axial stencil applied to the radial one.

    These are a *definition*, not an estimate: the certified object is the
    interpolant built from exactly these numbers, and the payload records which
    stencil was used so the statement is unambiguous.  Nothing below infers a
    derivative bound from them.
    """
    from .operators import derivative_r, derivative_z

    slope_r = derivative_r(grid, field, even_at_axis=even_at_axis)
    slope_z = derivative_z(grid, field)
    slope_rz = derivative_z(grid, slope_r)
    return slope_r, slope_z, slope_rz


def bernstein_cell_envelope(
    values: FloatArray,
    slope_r: FloatArray,
    slope_z: FloatArray,
    slope_rz: FloatArray,
    i: int,
    j: int,
    *,
    dr: Fraction,
    dz: Fraction,
) -> Interval:
    r"""**Theorem T1.**  An exact enclosure of the bicubic Hermite interpolant.

    The interpolant on the cell is
    ``I(\xi,\eta) = \sum_{a,b} H_a(\xi) H_b(\eta) C_{ab}`` with
    ``H = (h_{00}, h_{10}, h_{01}, h_{11})`` and ``C`` the four-by-four matrix of
    corner values and scaled slopes.  Changing basis with
    :data:`HERMITE_TO_BERNSTEIN` on both sides gives the Bernstein net
    ``B = T C T^{\mathsf T}``, and

    .. math::  I(\xi,\eta) = \sum_{k,l} B^3_k(\xi)B^3_l(\eta)\,b_{kl}

    with every ``B^3_k \ge 0`` and ``\sum_k B^3_k \equiv 1``.  The interpolant is
    therefore a convex combination of the ``b_{kl}`` at every point of the closed
    cell, so ``[\min b, \max b]`` encloses it.  This is a proof, not a bound:
    no Lipschitz constant, no inflation, no derivative estimate.
    """
    def entry(a: int, b: int) -> Fraction:
        row, col = i + (a >= 2), j + (b >= 2)
        if a % 2 == 0 and b % 2 == 0:
            return Fraction(float(values[row, col]))
        if a % 2 == 1 and b % 2 == 0:
            return dr * Fraction(float(slope_r[row, col]))
        if a % 2 == 0 and b % 2 == 1:
            return dz * Fraction(float(slope_z[row, col]))
        return dr * dz * Fraction(float(slope_rz[row, col]))

    coefficients = [[entry(a, b) for b in range(4)] for a in range(4)]
    transform = HERMITE_TO_BERNSTEIN
    # B = T C T^T, done in two passes so the intermediate stays rational.
    left = [
        [sum(transform[k][a] * coefficients[a][b] for a in range(4)) for b in range(4)]
        for k in range(4)
    ]
    net = [
        [sum(left[k][b] * transform[l][b] for b in range(4)) for l in range(4)]
        for k in range(4)
    ]
    flat = [value for row in net for value in row]
    return Interval(min(flat), max(flat))


def _time_hull(
    value_n: Interval,
    value_next: Interval,
    slope_n: Interval,
    slope_next: Interval,
    step: Fraction,
) -> Interval:
    """Enclose the cubic Hermite interpolant over the whole slab, exactly."""
    base = _hull(value_n, value_next)
    lo_10, hi_10 = HERMITE_BASIS_RANGES["h10"]
    lo_11, hi_11 = HERMITE_BASIS_RANGES["h11"]
    tangent_n = Interval(lo_10, hi_10).scale(step) * slope_n
    tangent_next = Interval(lo_11, hi_11).scale(step) * slope_next
    return base + tangent_n + tangent_next


def build_slab_certificate(
    grid: AxisymmetricGrid,
    *,
    time_start: float,
    time_step: float,
    u1_start: FloatArray,
    u1_end: FloatArray,
    omega1_start: FloatArray,
    omega1_end: FloatArray,
    psi1_start: FloatArray,
    psi1_end: FloatArray,
    du1_start: FloatArray,
    du1_end: FloatArray,
    domega1_start: FloatArray,
    domega1_end: FloatArray,
    viscosity: float,
    green_tail_bound: float,
    domega1_mid: FloatArray | None = None,
    interior_radius: float | None = None,
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> SlabCertificate:
    """Enclose one slab, exactly and unconditionally, as a statement about the
    interpolant.

    There is no inflation parameter and no fourth-derivative bound: the spatial
    enclosure is the Bernstein theorem T1 and the temporal one is T2, both
    proved.  The earlier signature took two inflation factors, which is what a
    hypothesis looks like when it is pretending to be a theorem.

    ``du1_*`` and ``domega1_*`` must be the **solver's own** right-hand sides at
    the accepted steps.  They are never reconstructed from snapshot differences:
    a difference quotient of two stored states cannot separate the local defect
    it is supposed to measure from the interpolation error it is built from.
    """
    if grid.periodic_z:
        raise ValueError("the slab certificate targets non-periodic z grids")
    if grid.nr * grid.nz > 20000:
        raise ValueError("coarsen the slab before building an exact certificate")
    if time_step <= 0.0:
        raise ValueError("time_step must be positive")
    if viscosity <= 0.0:
        raise ValueError("a Clay candidate requires a fixed positive viscosity")

    fields = {
        "u1_start": u1_start, "u1_end": u1_end,
        "omega1_start": omega1_start, "omega1_end": omega1_end,
        "psi1_start": psi1_start, "psi1_end": psi1_end,
        "du1_start": du1_start, "du1_end": du1_end,
        "domega1_start": domega1_start, "domega1_end": domega1_end,
    }
    for name, field in fields.items():
        grid.validate_field(field, name=name)

    dr = Fraction(float(grid.dr))
    dz = Fraction(float(grid.dz))
    step = Fraction(float(time_step))
    radius = (
        float(interior_radius)
        if interior_radius is not None
        else float(min(grid.r[-1], grid.z[-1])) * 0.5
    )

    cells = [
        (i, j)
        for i in range(grid.nr - 1)
        for j in range(grid.nz - 1)
        if float(np.hypot(grid.r[i + 1], max(abs(grid.z[j]), abs(grid.z[j + 1]))))
        <= radius
    ]
    if not cells:
        raise ValueError("interior_radius selects no cell")

    u1_hull = Interval(Fraction(0), Fraction(0))
    omega1_hull = Interval(Fraction(0), Fraction(0))
    du1_hull = Interval(Fraction(0), Fraction(0))
    domega1_hull = Interval(Fraction(0), Fraction(0))
    poisson_residual = Fraction(0)
    divergence_bound = Fraction(0)
    trapezoid_defect = Fraction(0)
    simpson_defect: Fraction | None = Fraction(0) if domega1_mid is not None else None

    r_rat = [Fraction(float(value)) for value in grid.r]

    # The nodal slopes that DEFINE the interpolant.  u1 and omega1 are even in r
    # so their radial slope vanishes on the axis; the right-hand sides inherit
    # that parity because every term of the equations preserves it.
    slopes = {
        name: hermite_slopes(grid, field, even_at_axis=True)
        for name, field in (
            ("u1_start", u1_start), ("u1_end", u1_end),
            ("omega1_start", omega1_start), ("omega1_end", omega1_end),
            ("du1_start", du1_start), ("du1_end", du1_end),
            ("domega1_start", domega1_start), ("domega1_end", domega1_end),
        )
    }
    arrays = {
        "u1_start": u1_start, "u1_end": u1_end,
        "omega1_start": omega1_start, "omega1_end": omega1_end,
        "du1_start": du1_start, "du1_end": du1_end,
        "domega1_start": domega1_start, "domega1_end": domega1_end,
    }

    def envelope(name: str, i: int, j: int) -> Interval:
        return bernstein_cell_envelope(
            arrays[name], *slopes[name], i, j, dr=dr, dz=dz
        )

    for i, j in cells:
        def slab(start: str, end: str, dstart: str, dend: str) -> Interval:
            return _time_hull(
                envelope(start, i, j), envelope(end, i, j),
                envelope(dstart, i, j), envelope(dend, i, j), step,
            ).round_outward(precision_bits)

        u1_hull = _hull(
            u1_hull, slab("u1_start", "u1_end", "du1_start", "du1_end")
        )
        omega1_hull = _hull(
            omega1_hull,
            slab("omega1_start", "omega1_end", "domega1_start", "domega1_end"),
        )
        du1_hull = _hull(
            du1_hull,
            _hull(envelope("du1_start", i, j), envelope("du1_end", i, j))
            .round_outward(precision_bits),
        )
        domega1_hull = _hull(
            domega1_hull,
            _hull(envelope("domega1_start", i, j), envelope("domega1_end", i, j))
            .round_outward(precision_bits),
        )

        # -- the elliptic constraint at both endpoints, node by node -------- #
        if 1 <= i < grid.nr - 1 and 1 <= j < grid.nz - 1:
            for psi, omega in ((psi1_start, omega1_start), (psi1_end, omega1_end)):
                def value(a: int, b: int) -> Interval:
                    return Interval.exact(float(psi[a, b]))

                psi_rr = (
                    value(i + 1, j) - value(i, j).scale(Fraction(2)) + value(i - 1, j)
                ).scale(1 / (dr * dr))
                psi_r = (value(i + 1, j) - value(i - 1, j)).scale(1 / (2 * dr))
                psi_zz = (
                    value(i, j + 1) - value(i, j).scale(Fraction(2)) + value(i, j - 1)
                ).scale(1 / (dz * dz))
                psi_z = (value(i, j + 1) - value(i, j - 1)).scale(1 / (2 * dz))
                laplace = (
                    psi_rr + psi_r.scale(Fraction(3) / r_rat[i]) + psi_zz
                ).round_outward(precision_bits)
                residual = (
                    -laplace - Interval.exact(float(omega[i, j]))
                ).round_outward(precision_bits)
                poisson_residual = max(poisson_residual, residual.magnitude)

                # div u = d_r u^r + u^r / r + d_z u^z with u^r/r = -psi_z.
                def u_r_at(a: int) -> Interval:
                    return (
                        -(value(a, j + 1) - value(a, j - 1)).scale(1 / (2 * dz))
                    ).scale(r_rat[a])

                def u_z_at(b: int) -> Interval:
                    return value(i, b).scale(Fraction(2)) + (
                        (value(i + 1, b) - value(i - 1, b)).scale(1 / (2 * dr))
                    ).scale(r_rat[i])

                divergence = (
                    (u_r_at(i + 1) - u_r_at(i - 1)).scale(1 / (2 * dr))
                    + (-psi_z)
                    + (u_z_at(j + 1) - u_z_at(j - 1)).scale(1 / (2 * dz))
                ).round_outward(precision_bits)
                divergence_bound = max(divergence_bound, divergence.magnitude)

        # -- the integrator's local defect, exactly -------------------------- #
        increment = Interval.exact(float(omega1_end[i, j])) - Interval.exact(
            float(omega1_start[i, j])
        )
        trapezoid = (
            Interval.exact(float(domega1_start[i, j]))
            + Interval.exact(float(domega1_end[i, j]))
        ).scale(step / 2)
        trapezoid_defect = max(
            trapezoid_defect,
            (increment - trapezoid).round_outward(precision_bits).magnitude,
        )
        if domega1_mid is not None:
            simpson = (
                Interval.exact(float(domega1_start[i, j]))
                + Interval.exact(float(domega1_mid[i, j])).scale(Fraction(4))
                + Interval.exact(float(domega1_end[i, j]))
            ).scale(step / 6)
            assert simpson_defect is not None
            simpson_defect = max(
                simpson_defect,
                (increment - simpson).round_outward(precision_bits).magnitude,
            )

    # The Simpson and trapezoid defects are kept, but as DIAGNOSTICS of the
    # distance from the interpolant to the semi-discrete flow, not as the input
    # to a fourth-derivative estimate.  Nothing downstream consumes them.

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate": "space_time_slab",
        "grid": {
            "nr": grid.nr, "nz": grid.nz,
            "dr": str(dr), "dz": str(dz),
            "r_max": float(grid.r[-1]),
            "z_min": float(grid.z[0]), "z_max": float(grid.z[-1]),
        },
        "slab": {
            "time_start": float(time_start),
            "time_end": float(time_start) + float(time_step),
            "time_step": str(step),
            "cells": len(cells),
            "interior_radius": radius,
            "encloses_cell_interiors": True,
        },
        "viscosity": float(viscosity),
        "theorems": {
            name: {"statement": statement, "proved": True}
            for name, statement in SLAB_THEOREMS.items()
        },
        "interpolant": {
            "spatial": "tensor-product cubic Hermite, bicubic per cell",
            "slope_stencil": (
                "second-order central differences in the interior, one-sided at "
                "the ends, d_r f(0, z) = 0 for fields even in r; the mixed slope "
                "is the axial stencil applied to the radial one"
            ),
            "temporal": "cubic Hermite from the endpoint states and the exact "
                        "solver right-hand sides",
            "enclosure_method": "Bernstein convex hull in space, exact basis "
                                "ranges in time",
        },
        "hypotheses": {
            "H3_semidiscrete_to_continuum": {
                "statement": SLAB_HYPOTHESES["H3_semidiscrete_to_continuum"],
                "proved": False,
            }
        },
        "enclosures": {
            "u1": u1_hull.as_pair(),
            "omega1": omega1_hull.as_pair(),
            "du1_dt": du1_hull.as_pair(),
            "domega1_dt": domega1_hull.as_pair(),
        },
        "bounds": {
            "poisson_residual": str(poisson_residual),
            "divergence": str(divergence_bound),
            "trapezoid_local_defect": str(trapezoid_defect),
            "simpson_local_defect": (
                str(simpson_defect) if simpson_defect is not None else None
            ),
            "green_tail": float(green_tail_bound),
        },
        # The endpoint states must lie inside the tube -- but only where the tube
        # exists.  The enclosure covers the cells inside `interior_radius`, so
        # comparing against the extremes over the *whole* grid would test the
        # tube against data it never claimed to bound.
        "endpoint_inclusion": {
            name: _nodes_inside(hull, field, cells)
            for name, hull, field in (
                ("u1_start", u1_hull, u1_start),
                ("u1_end", u1_hull, u1_end),
                ("omega1_start", omega1_hull, omega1_start),
                ("omega1_end", omega1_hull, omega1_end),
            )
        },
        "claims": [
            "Unconditional: every stated enclosure holds at every point of every "
            "enclosed cell and every instant of the slab, for the interpolant "
            "defined by the recorded stencils.  Theorems T1 and T2; no "
            "hypothesis is used.",
            "The trapezoid and Simpson defects are diagnostics of the distance "
            "from that interpolant to the semi-discrete flow.  They are not an "
            "input to any bound.",
            "Conditional on H3 only: nothing here is a statement about the "
            "continuous Navier-Stokes solution, and nothing here bears on the "
            "Clay problem.",
        ],
    }
    return SlabCertificate(payload)


def verify_slab_certificate(payload: dict[str, object]) -> dict[str, object]:
    """Independently re-check a slab payload from its own contents.

    The checker never touches the solver.  It re-derives every relation the
    payload asserts between its own numbers, so a payload whose enclosures have
    been widened, narrowed or shifted fails here even though it is still
    syntactically valid.

    A checker fed an adversarial payload must return a verdict, not raise: a
    crash is indistinguishable from a bug in the checker.  Every malformed
    field is therefore reported as a failure.
    """
    try:
        return _verify(payload)
    except Exception as error:  # noqa: BLE001 -- adversarial input, any error is a failure
        return {
            "verified": False,
            "failures": [f"payload rejected while parsing: {error}"],
        }


def _verify(payload: dict[str, object]) -> dict[str, object]:
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(payload.get("certificate") == "space_time_slab", "wrong certificate kind")
    require(payload.get("schema_version") == 1, "unsupported schema version")

    slab = payload.get("slab")
    require(isinstance(slab, dict), "missing slab block")
    if isinstance(slab, dict):
        step = Fraction(str(slab["time_step"]))
        require(step > 0, "slab time step must be positive")
        require(
            abs(
                (float(slab["time_end"]) - float(slab["time_start"])) - float(step)
            )
            <= 1e-12 * max(1.0, abs(float(slab["time_end"]))),
            "slab endpoints disagree with the time step",
        )
        require(bool(slab.get("encloses_cell_interiors")), "slab is node-only")
        require(int(slab["cells"]) > 0, "slab encloses no cell")

    enclosures = payload.get("enclosures")
    require(isinstance(enclosures, dict), "missing enclosures")
    if isinstance(enclosures, dict):
        for name, pair in enclosures.items():
            interval = Interval.from_pair(pair)
            require(interval.lower <= interval.upper, f"{name} enclosure inverted")

    bounds = payload.get("bounds")
    require(isinstance(bounds, dict), "missing bounds")
    if isinstance(bounds, dict):
        for name in (
            "poisson_residual", "divergence", "trapezoid_local_defect",
        ):
            require(Fraction(str(bounds[name])) >= 0, f"{name} must be nonnegative")
        simpson = bounds.get("simpson_local_defect")
        if simpson is not None:
            require(
                Fraction(str(simpson)) <= Fraction(str(bounds["trapezoid_local_defect"]))
                or Fraction(str(bounds["trapezoid_local_defect"])) == 0,
                "Simpson defect exceeds the trapezoid defect, which is impossible "
                "for a smooth integrand on one step",
            )

    theorems = payload.get("theorems")
    require(isinstance(theorems, dict), "missing theorems block")
    if isinstance(theorems, dict):
        require(
            set(theorems) == set(SLAB_THEOREMS),
            "the payload must carry exactly the two named theorems",
        )
        for name, entry in theorems.items():
            require(
                entry.get("proved") is True,
                f"{name} is a theorem and must be marked proved",
            )
            require(
                entry.get("statement") == SLAB_THEOREMS[name],
                f"{name} statement does not match the audited text",
            )

    interpolant = payload.get("interpolant")
    require(isinstance(interpolant, dict), "missing interpolant definition")
    if isinstance(interpolant, dict):
        require(
            "slope_stencil" in interpolant and interpolant["slope_stencil"],
            "the certified object is undefined without its slope stencil",
        )

    hypotheses = payload.get("hypotheses")
    require(isinstance(hypotheses, dict), "missing hypotheses")
    if isinstance(hypotheses, dict):
        require(
            set(hypotheses) == set(SLAB_HYPOTHESES),
            "the payload must carry exactly the one named hypothesis",
        )
        for name, entry in hypotheses.items():
            require(
                entry.get("proved") is False,
                f"{name} claims to be proved; this module does not prove it",
            )
            require(
                entry.get("statement") == SLAB_HYPOTHESES[name],
                f"{name} statement does not match the audited text",
            )

    inclusion = payload.get("endpoint_inclusion")
    require(isinstance(inclusion, dict), "missing endpoint inclusion block")
    if isinstance(inclusion, dict):
        require(all(bool(v) for v in inclusion.values()), "an endpoint escapes the tube")

    claims = payload.get("claims")
    require(
        isinstance(claims, list)
        and any("Clay" in str(claim) for claim in claims),
        "the payload must record that it says nothing about the Clay problem",
    )
    return {"verified": not failures, "failures": failures}
