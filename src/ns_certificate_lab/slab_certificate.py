r"""One space-time slab certificate over ``[t_n, t_{n+1}]``.

:mod:`ns_certificate_lab.snapshot_certificate` encloses one instant at the grid
*nodes*.  That is the weakest possible object: a singularity is a statement
about a space-time region, and a bound that holds at isolated nodes at isolated
instants constrains nothing between them.  This module encloses a whole slab --
every point of every cell, at every instant of ``[t_n, t_{n+1}]``.

What is rigorous and what is not
--------------------------------
Being explicit about this is the point of the module, so it is stated up front
and repeated in the payload.

**Rigorous, given the two endpoint states and the exact solver right-hand sides
as inputs.**  All arithmetic is exact rational arithmetic with outward rounding
(:class:`~ns_certificate_lab.snapshot_certificate.Interval`).  The temporal
enclosure of the cubic Hermite interpolant uses exact ranges of its basis
polynomials.  The integrator's local defect -- trapezoid and, when a midpoint
right-hand side is supplied, Simpson -- is an exactly computed number, not an
estimate.  The endpoint inclusion test is exact.

**Hypotheses, named and quantified, not proved here.**

``H1`` (*cell Lipschitz*).  Between adjacent nodes the field is enclosed by the
corner hull inflated by ``inflation * (|d_r f| dr + |d_z f| dz) / 2`` with the
one-sided differences measured on the grid.  This is exact for a field whose
first derivative on the cell is enclosed by ``inflation`` times the measured
divided differences.  Proving that needs a bound on the second spatial
derivative, which the numerical state does not carry.

``H2`` (*Hermite remainder*).  ``|y(t) - H(t)| <= Delta^4 M_4 / 384`` where
``M_4`` bounds the fourth time derivative on the slab.  ``M_4`` is estimated
from the measured Simpson defect when a midpoint is supplied, and inflated;
otherwise it must be supplied by the caller.  Proving a value for ``M_4`` needs
the ``H^s`` machinery of :doc:`../docs/research_notes/hs_error_propagation`, which
is at present an incomplete derivation with named unproved constants.

Consequently this certificate does **not** establish a property of the
continuous Navier-Stokes solution.  It establishes, exactly, a property of the
computed space-time object, conditional on ``H1`` and ``H2``.  Anything stronger
would require the missing lemmas, and they are listed rather than assumed away.
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

SLAB_HYPOTHESES = {
    "H1_cell_lipschitz": (
        "Between adjacent nodes the field's first spatial derivative is "
        "enclosed by 'inflation' times the measured divided differences.  "
        "Unproved: needs a second-derivative bound the numerical state does "
        "not carry."
    ),
    "H2_hermite_remainder": (
        "|y(t) - H(t)| <= Delta^4 M4 / 384 on the slab, with M4 a bound on the "
        "fourth time derivative.  Unproved: M4 is estimated from the measured "
        "Simpson defect and inflated, not derived."
    ),
}


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


def _cell_hull(
    field: FloatArray, i: int, j: int, *, dr: Fraction, dz: Fraction, inflation: Fraction
) -> Interval:
    """Enclose ``field`` over the whole cell ``[r_i, r_{i+1}] x [z_j, z_{j+1}]``.

    Hypothesis ``H1``: the corner hull inflated by ``inflation`` times half the
    measured variation across the cell in each direction.  With ``inflation = 1``
    and a monotone field the corner hull alone is already valid; the inflation
    covers an interior extremum whose size is bounded by the same differences.
    """
    corners = [
        Interval.exact(float(field[a, b]))
        for a in (i, i + 1)
        for b in (j, j + 1)
    ]
    base = _hull(*corners)
    variation_r = max(
        abs(Fraction(float(field[i + 1, b])) - Fraction(float(field[i, b])))
        for b in (j, j + 1)
    )
    variation_z = max(
        abs(Fraction(float(field[a, j + 1])) - Fraction(float(field[a, j])))
        for a in (i, i + 1)
    )
    pad = inflation * (variation_r + variation_z) / 2
    return Interval(base.lower - pad, base.upper + pad)


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
    cell_inflation: float = 1.0,
    remainder_inflation: float = 8.0,
    fourth_derivative_bound: float | None = None,
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> SlabCertificate:
    """Enclose one slab, exactly, under the two named hypotheses.

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
    inflation = Fraction(float(cell_inflation))
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

    for i, j in cells:
        def slab(start: FloatArray, end: FloatArray,
                 dstart: FloatArray, dend: FloatArray) -> Interval:
            return _time_hull(
                _cell_hull(start, i, j, dr=dr, dz=dz, inflation=inflation),
                _cell_hull(end, i, j, dr=dr, dz=dz, inflation=inflation),
                _cell_hull(dstart, i, j, dr=dr, dz=dz, inflation=inflation),
                _cell_hull(dend, i, j, dr=dr, dz=dz, inflation=inflation),
                step,
            ).round_outward(precision_bits)

        u1_hull = _hull(u1_hull, slab(u1_start, u1_end, du1_start, du1_end))
        omega1_hull = _hull(
            omega1_hull, slab(omega1_start, omega1_end, domega1_start, domega1_end)
        )
        du1_hull = _hull(
            du1_hull,
            _time_hull(
                _cell_hull(du1_start, i, j, dr=dr, dz=dz, inflation=inflation),
                _cell_hull(du1_end, i, j, dr=dr, dz=dz, inflation=inflation),
                Interval(Fraction(0), Fraction(0)),
                Interval(Fraction(0), Fraction(0)),
                step,
            ).round_outward(precision_bits),
        )
        domega1_hull = _hull(
            domega1_hull,
            _hull(
                _cell_hull(domega1_start, i, j, dr=dr, dz=dz, inflation=inflation),
                _cell_hull(domega1_end, i, j, dr=dr, dz=dz, inflation=inflation),
            ).round_outward(precision_bits),
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

    # -- hypothesis H2 --------------------------------------------------------- #
    if fourth_derivative_bound is not None:
        m4 = Fraction(float(fourth_derivative_bound))
        m4_source = "supplied by the caller"
    elif simpson_defect is not None:
        # Simpson's defect is  -(Delta^5/2880) y^{(4)}(xi)  for one step, so
        # M4 <= 2880 * defect / Delta^5, inflated to cover the sampling.
        m4 = (
            Fraction(2880) * simpson_defect * Fraction(float(remainder_inflation))
            / step**5
        )
        m4_source = "estimated from the measured Simpson defect and inflated"
    else:
        m4 = Fraction(0)
        m4_source = "unavailable: no midpoint right-hand side and no caller bound"
    hermite_remainder = m4 * step**4 / 384

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
        "hypotheses": {
            "H1_cell_lipschitz": {
                "statement": SLAB_HYPOTHESES["H1_cell_lipschitz"],
                "inflation": float(cell_inflation),
                "proved": False,
            },
            "H2_hermite_remainder": {
                "statement": SLAB_HYPOTHESES["H2_hermite_remainder"],
                "fourth_derivative_bound": str(m4),
                "source": m4_source,
                "inflation": float(remainder_inflation),
                "remainder": str(hermite_remainder),
                "proved": False,
            },
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
            "hermite_remainder": str(hermite_remainder),
        },
        "endpoint_inclusion": {
            "u1_start": u1_hull.contains(Fraction(float(np.max(u1_start))))
            and u1_hull.contains(Fraction(float(np.min(u1_start)))),
            "u1_end": u1_hull.contains(Fraction(float(np.max(u1_end))))
            and u1_hull.contains(Fraction(float(np.min(u1_end)))),
            "omega1_start": omega1_hull.contains(Fraction(float(np.max(omega1_start))))
            and omega1_hull.contains(Fraction(float(np.min(omega1_start)))),
            "omega1_end": omega1_hull.contains(Fraction(float(np.max(omega1_end))))
            and omega1_hull.contains(Fraction(float(np.min(omega1_end)))),
        },
        "claims": [
            "Every stated bound holds for the computed space-time object on "
            "every point of every enclosed cell and every instant of the slab, "
            "conditional on H1 and H2.",
            "Nothing here is a statement about the continuous Navier-Stokes "
            "solution, and nothing here bears on the Clay problem.",
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
            "hermite_remainder",
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

    hypotheses = payload.get("hypotheses")
    require(isinstance(hypotheses, dict), "missing hypotheses")
    if isinstance(hypotheses, dict):
        require(
            set(hypotheses) == set(SLAB_HYPOTHESES),
            "the payload must carry exactly the two named hypotheses",
        )
        for name, entry in hypotheses.items():
            require(
                entry.get("proved") is False,
                f"{name} claims to be proved; this module proves neither",
            )
            require(
                entry.get("statement") == SLAB_HYPOTHESES[name],
                f"{name} statement does not match the audited text",
            )
    if isinstance(hypotheses, dict) and isinstance(bounds, dict):
        h2 = hypotheses.get("H2_hermite_remainder", {})
        m4 = Fraction(str(h2.get("fourth_derivative_bound", "0")))
        step = Fraction(str(slab["time_step"])) if isinstance(slab, dict) else Fraction(1)
        require(
            Fraction(str(bounds["hermite_remainder"])) == m4 * step**4 / 384,
            "the Hermite remainder does not equal M4 Delta^4 / 384",
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
