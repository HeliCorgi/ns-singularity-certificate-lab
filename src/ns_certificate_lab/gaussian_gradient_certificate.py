r"""Gradient-form pressure pairing for the Gaussian--Hermite generation certificate.

The certified object
--------------------
The generation-rate certificate of :mod:`ns_certificate_lab.interval_gaussian`
bounds the pressure contribution to ``J(u_0)`` in the divergence form
``3\int p_h\,g\,dx`` with ``g = \nabla\cdot(|u|u)``: the discrete pressure
enters through its per-cell corner hull, the signed flux through its own
interval, and the certified lower bound is the sum of the worst corner
products.  Measured on the Gaussian--Hermite search candidate, the binding
width of that certificate is exactly this product — the pressure hull and the
flux enclosure are treated as independent, so the certificate pays
``|p|\cdot{\rm width}(g) + {\rm width}(p)\cdot|g|`` per cell even though the
integrand the float reference actually evaluates is the *paired* quantity.

This module re-certifies the same pressure contribution in the
integrated-by-parts (gradient) form the float reference already uses,

.. math::  P \;=\; -3\int_\Omega |u|\,\bigl(u^r\,\partial_r p_h
           + u^z\,\partial_z p_h\bigr)\,dx ,

with ``p_h`` the **multilinear interpolant** of the computed nodal pressure.
The gradient of that interpolant on a cell is an *exact rational linear
combination of the four nodal values*:

.. math::  \partial_r p_h(\xi,\eta) = \frac{(1-\eta)\,(p_{10}-p_{00})
           + \eta\,(p_{11}-p_{01})}{\Delta r},

linear in ``\eta`` and constant in ``\xi`` (and symmetrically for
``\partial_z p_h``).  Its range over any sub-box is therefore attained at the
two endpoints — a theorem of the interpolant, exactly like the corner-hull
range used for ``p_h`` itself — and the differentiation is applied to the
pressure *vector* symbolically, so the pressure and its gradient never receive
independent hulls.  The velocity factor ``|u|\,u`` is enclosed by the same
Gaussian--Hermite interval mechanics as every other analytic quantity
(:func:`~ns_certificate_lab.interval_gaussian.cell_velocity_enclosure`).

The three certificate paths built here
--------------------------------------
* ``path="plain"`` — the per-cell interval product
  ``[r\,|u|u^r]\cdot[\partial_r p_h] + [r\,|u|u^z]\cdot[\partial_z p_h]``;
* ``path="adaptive"`` — the same, with budgeted 2x2 subdivision of the cells
  on which the paired integrand is not sign-stable (its enclosure straddles
  zero), recursively up to a depth limit;
* ``path="taylor"`` — the pressure gradient carried as a Taylor model per
  cell: centre value (an exact rational) times the enclosed velocity moment,
  plus a remainder bounded by the within-cell gradient variation.  The
  Taylor-model contribution is intersected with the plain product, so it is
  never wider.

Every path also accumulates, cell by cell, the *natural* product-hull width
and the *correlated* (Taylor-model) width of the same contribution; their
difference is the measured dependency overestimation of the pairing.

What is proved and what is not
------------------------------
Unconditionally proved, in exact rational arithmetic: the stated enclosure of
the box integral ``-3\int_\Omega |u|\,u\cdot\nabla p_h\,dx`` for the *computed*
grid function ``p_h``, the upper bound on the viscous integral (interior
enclosure plus the Mills-ratio exterior tail), and the sensitivity integral.
The step from ``p_h`` to the true free-space pressure is the single named,
**unproved** hypothesis ``P1G`` below — the ``W^{1,\infty}`` analogue of the
``P1`` hypothesis of :mod:`ns_certificate_lab.l3_certificate`, belonging to
the same unproved circle.  The exterior pressure remainder is *not* enclosed,
exactly as in the divergence-form certificate, and the conditional claim is
restricted to the box contribution.

Nothing here is a singularity result and nothing here bears on the Clay
problem.  A positive initial generation rate of the critical norm would be a
screening criterion, not evidence of blow-up; a negative or failed margin is
evidence of nothing at all beyond the width of the enclosure.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .interval_gaussian import (
    cell_enclosure_gaussian,
    cell_velocity_enclosure,
    exterior_viscous_tail,
    gaussian_reference_field,
)
from .l3_certificate import GenerationCertificate, _cell_hull_from_nodes
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "GRADIENT_PRESSURE_HYPOTHESIS",
    "GRADIENT_CERTIFICATE_KIND",
    "CellPressureGradient",
    "build_gradient_generation_certificate",
    "compare_certificate_paths",
    "pressure_cell_gradient",
    "verify_gradient_generation_certificate",
]

GRADIENT_CERTIFICATE_KIND = "l3_generation_rate_gradient_form"

GRADIENT_PRESSURE_HYPOTHESIS = (
    "P1G (pressure-gradient discretisation): ||grad(p - p_h)||_{L^inf(Omega)} "
    "<= tau_star_grad, with p the free-space solution of -Laplacian p = "
    "d_i u_j d_j u_i and p_h the multilinear interpolant of the computed nodal "
    "pressure.  Unproved here: it needs an a posteriori W^{1,inf} bound on a "
    "Poisson solve, which belongs to the same unproved circle as P1."
)

_PATHS = ("plain", "adaptive", "taylor")

#: ``2\pi`` enclosed between two rationals (``2\pi = 6.283185307179586476...``).
_TWO_PI = Interval(
    Fraction(6283185307179586, 10**15), Fraction(6283185307179587, 10**15)
)
_TWO_PI_UPPER = _TWO_PI.upper
_ZERO = Interval(Fraction(0), Fraction(0))


# --------------------------------------------------------------------------- #
# the exact per-cell gradient of the multilinear pressure interpolant           #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CellPressureGradient:
    r"""The gradient of the multilinear ``p_h`` on one cell, exactly.

    ``d_r_low``/``d_r_high`` are the two radial edge difference quotients
    ``(p_{10}-p_{00})/\Delta r`` (at ``z = z_j``) and ``(p_{11}-p_{01})/\Delta
    r`` (at ``z = z_{j+1}``); ``d_z_low``/``d_z_high`` the axial ones at
    ``r = r_i`` and ``r = r_{i+1}``.  All four are exact rationals — the
    difference operator applied to the pressure vector symbolically, never a
    hull of anything.  ``\partial_r p_h`` is the linear interpolation of the
    first pair in ``\eta`` and ``\partial_z p_h`` of the second pair in
    ``\xi``, so the range over any sub-box is attained at the endpoints and
    the two range methods return **exact** ranges of the interpolant gradient,
    not enclosures with slack.
    """

    d_r_low: Fraction
    d_r_high: Fraction
    d_z_low: Fraction
    d_z_high: Fraction

    def radial_range(self, eta0: Fraction, eta1: Fraction) -> Interval:
        r"""Range of ``\partial_r p_h`` for ``\eta \in [\eta_0, \eta_1]``."""
        if not Fraction(0) <= eta0 <= eta1 <= Fraction(1):
            raise ValueError("the eta sub-range must sit inside [0, 1]")
        first = self.d_r_low + (self.d_r_high - self.d_r_low) * eta0
        second = self.d_r_low + (self.d_r_high - self.d_r_low) * eta1
        return Interval(min(first, second), max(first, second))

    def axial_range(self, xi0: Fraction, xi1: Fraction) -> Interval:
        r"""Range of ``\partial_z p_h`` for ``\xi \in [\xi_0, \xi_1]``."""
        if not Fraction(0) <= xi0 <= xi1 <= Fraction(1):
            raise ValueError("the xi sub-range must sit inside [0, 1]")
        first = self.d_z_low + (self.d_z_high - self.d_z_low) * xi0
        second = self.d_z_low + (self.d_z_high - self.d_z_low) * xi1
        return Interval(min(first, second), max(first, second))


def pressure_cell_gradient(
    pressure: FloatArray,
    i: int,
    j: int,
    dr: Fraction,
    dz: Fraction,
) -> CellPressureGradient:
    """The exact rational gradient data of cell ``(i, j)``.

    Every binary64 nodal value converts to a rational exactly, so nothing is
    approximated: the returned quotients are the difference operator applied
    to the stored pressure vector, in exact arithmetic.
    """
    if dr <= 0 or dz <= 0:
        raise ValueError("cell steps must be positive")
    p00 = Fraction(float(pressure[i, j]))
    p10 = Fraction(float(pressure[i + 1, j]))
    p01 = Fraction(float(pressure[i, j + 1]))
    p11 = Fraction(float(pressure[i + 1, j + 1]))
    return CellPressureGradient(
        d_r_low=(p10 - p00) / dr,
        d_r_high=(p11 - p01) / dr,
        d_z_low=(p01 - p00) / dz,
        d_z_high=(p11 - p10) / dz,
    )


def _intersect(first: Interval, second: Interval) -> Interval:
    """Intersection of two valid enclosures of the same quantity.

    Both arguments enclose the same true range, so the intersection is again a
    valid enclosure; an empty intersection can only come from a bug, and the
    defensive fallback returns the first argument rather than an invalid
    interval.
    """
    lower = max(first.lower, second.lower)
    upper = min(first.upper, second.upper)
    if lower > upper:  # both enclose the truth, so this cannot happen; be safe
        return first
    return Interval(lower, upper)


# --------------------------------------------------------------------------- #
# the builder                                                                   #
# --------------------------------------------------------------------------- #


def build_gradient_generation_certificate(
    family,
    grid: AxisymmetricGrid,
    *,
    viscosity: float,
    path: str = "plain",
    precision_bits: int = DEFAULT_PRECISION_BITS,
    exp_terms: int = 32,
    max_depth: int = 2,
    refine_budget: int = 4096,
) -> GenerationCertificate:
    r"""Build the gradient-form generation certificate on one grid.

    ``path`` selects the pairing strategy (see the module docstring):
    ``"plain"`` interval products, ``"adaptive"`` sign-stability subdivision
    with at most ``refine_budget`` extra leaf evaluations and ``max_depth``
    binary splits per direction, or ``"taylor"`` centre-plus-remainder
    pressure-gradient models intersected with the plain product.

    The payload records, besides the margin, the separated width accounting:
    the natural product-hull width, the correlated (Taylor-model) width, and
    their difference — the measured dependency overestimation of the pairing —
    together with the certified bound on the within-cell pressure-gradient
    variation and the exterior tails.  All statements about the true PDE
    pressure remain conditional on the named hypothesis ``P1G``; the enclosure
    of the discrete integral itself is unconditional.
    """
    from .l3_generation import AxisymmetricPressureSolver, l3_generation_rate
    from .l3_optimizer import require_clay_admissible

    require_clay_admissible(family)
    if viscosity <= 0.0:
        raise ValueError("a Clay candidate requires a fixed positive viscosity")
    if grid.periodic_z:
        raise ValueError("the generation certificate targets non-periodic z grids")
    if path not in _PATHS:
        raise ValueError(f"unknown path {path!r}; expected one of {_PATHS}")
    if max_depth < 0 or refine_budget < 0:
        raise ValueError("max_depth and refine_budget must be nonnegative")
    if not (float(grid.z[0]) < 0.0 < float(grid.z[-1])):
        raise ValueError(
            "the exterior tail bound needs a box that surrounds the origin: "
            "z_min < 0 < z_max"
        )

    field = gaussian_reference_field(family, grid)
    solver = AxisymmetricPressureSolver.build(grid)
    rate = l3_generation_rate(field, viscosity=viscosity, solver=solver)
    gradient = field.gradient()
    trace = solver.exterior_trace(field.u_r, field.u_theta, field.u_z)
    pressure = solver.solve(field.pressure_source(gradient), trace=trace)

    r_nodes = [Fraction(float(v)) for v in grid.r]
    z_nodes = [Fraction(float(v)) for v in grid.z]

    pressure_term = _ZERO                       # encloses -3 int |u| u . grad p_h
    viscous_upper = Fraction(0)
    sensitivity_upper = Fraction(0)
    natural_width = Fraction(0)
    correlated_width = Fraction(0)
    gradient_variation_bound = Fraction(0)
    gradient_spread_max = Fraction(0)
    cells = 0
    leaves = 0
    refined_cells = 0
    budget = refine_budget

    def leaf_terms(
        cell: CellPressureGradient,
        r_lo: Fraction, r_hi: Fraction, z_lo: Fraction, z_hi: Fraction,
        cell_r0: Fraction, cell_z0: Fraction, dr: Fraction, dz: Fraction,
    ) -> tuple:
        """Enclose one leaf box; returns everything the accumulators need."""
        r_box = Interval(r_lo, r_hi)
        z_box = Interval(z_lo, z_hi)
        v = cell_velocity_enclosure(
            family, r_box, z_box, terms=exp_terms, bits=precision_bits
        )
        g_r = cell.radial_range((z_lo - cell_z0) / dz, (z_hi - cell_z0) / dz)
        g_z = cell.axial_range((r_lo - cell_r0) / dr, (r_hi - cell_r0) / dr)
        rf_r = (r_box * (v.speed * v.u_r)).round_outward(precision_bits)
        rf_z = (r_box * (v.speed * v.u_z)).round_outward(precision_bits)
        area = (r_hi - r_lo) * (z_hi - z_lo)

        paired = rf_r * g_r + rf_z * g_z
        natural = (_TWO_PI * paired.scale(area)).round_outward(precision_bits)

        centre_r = (g_r.lower + g_r.upper) / 2
        centre_z = (g_z.lower + g_z.upper) / 2
        dev_r = (g_r.upper - g_r.lower) / 2
        dev_z = (g_z.upper - g_z.lower) / 2
        remainder = _TWO_PI_UPPER * area * (
            dev_r * rf_r.magnitude + dev_z * rf_z.magnitude
        )
        core = rf_r.scale(centre_r) + rf_z.scale(centre_z)
        correlated = (
            _TWO_PI * core.scale(area) + Interval(-remainder, remainder)
        ).round_outward(precision_bits)

        weight_up = _TWO_PI_UPPER * r_hi * area
        return v, natural, correlated, remainder, dev_r, dev_z, weight_up

    for i in range(grid.nr - 1):
        for j in range(grid.nz - 1):
            cells += 1
            cell_r0, cell_r1 = r_nodes[i], r_nodes[i + 1]
            cell_z0, cell_z1 = z_nodes[j], z_nodes[j + 1]
            dr = cell_r1 - cell_r0
            dz = cell_z1 - cell_z0
            cell = pressure_cell_gradient(pressure, i, j, dr, dz)
            gradient_spread_max = max(
                gradient_spread_max,
                abs(cell.d_r_high - cell.d_r_low),
                abs(cell.d_z_high - cell.d_z_low),
            )

            stack = [(cell_r0, cell_r1, cell_z0, cell_z1, 0)]
            cell_was_refined = False
            while stack:
                r_lo, r_hi, z_lo, z_hi, depth = stack.pop()
                v, natural, correlated, remainder, dev_r, dev_z, weight_up = (
                    leaf_terms(
                        cell, r_lo, r_hi, z_lo, z_hi, cell_r0, cell_z0, dr, dz
                    )
                )
                if (
                    path == "adaptive"
                    and depth < max_depth
                    and budget >= 4
                    and natural.contains_zero
                    and natural.upper > natural.lower
                ):
                    budget -= 4
                    cell_was_refined = True
                    r_mid = (r_lo + r_hi) / 2
                    z_mid = (z_lo + z_hi) / 2
                    stack.extend(
                        (
                            (r_lo, r_mid, z_lo, z_mid, depth + 1),
                            (r_mid, r_hi, z_lo, z_mid, depth + 1),
                            (r_lo, r_mid, z_mid, z_hi, depth + 1),
                            (r_mid, r_hi, z_mid, z_hi, depth + 1),
                        )
                    )
                    continue

                leaves += 1
                used = _intersect(natural, correlated) if path == "taylor" else natural
                pressure_term = pressure_term + used.scale(Fraction(-3))
                natural_width += 3 * (natural.upper - natural.lower)
                correlated_width += 3 * (correlated.upper - correlated.lower)
                gradient_variation_bound += 3 * remainder
                viscous_upper += weight_up * max(
                    v.viscous_integrand.upper, Fraction(0)
                )
                sensitivity_upper += weight_up * max(
                    v.speed_squared.upper, Fraction(0)
                )
            if cell_was_refined:
                refined_cells += 1
        pressure_term = pressure_term.round_outward(precision_bits)
        viscous_upper = Interval(viscous_upper, viscous_upper).round_outward(
            precision_bits
        ).upper
        sensitivity_upper = Interval(
            sensitivity_upper, sensitivity_upper
        ).round_outward(precision_bits).upper
        natural_width = Interval(natural_width, natural_width).round_outward(
            precision_bits
        ).upper
        correlated_width = Interval(
            correlated_width, correlated_width
        ).round_outward(precision_bits).upper
        gradient_variation_bound = Interval(
            gradient_variation_bound, gradient_variation_bound
        ).round_outward(precision_bits).upper

    tail = exterior_viscous_tail(
        family,
        r_nodes[-1],
        min(-z_nodes[0], z_nodes[-1]),
        terms=exp_terms,
        bits=precision_bits,
    )
    interior_viscous = viscous_upper
    viscous_total = viscous_upper + tail.viscous_tail

    nu = Fraction(float(viscosity))
    viscous_bound = -3 * nu * viscous_total     # a lower bound on V
    pressure_lower = pressure_term.lower        # a lower bound on the box P
    sensitivity_integral = 3 * sensitivity_upper
    margin = pressure_lower + viscous_bound
    threshold = (
        margin / sensitivity_integral if sensitivity_integral > 0 else Fraction(0)
    )

    distance = float(tail.distance)
    velocity_tail = float(tail.speed_coefficient) * math.exp(
        -float(tail.gamma) * distance * distance
    )

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate": GRADIENT_CERTIFICATE_KIND,
        "family": family.as_dict(),
        "grid": {
            "nr": grid.nr, "nz": grid.nz,
            "r_max": float(grid.r[-1]),
            "z_min": float(grid.z[0]), "z_max": float(grid.z[-1]),
        },
        "viscosity": float(viscosity),
        "path": path,
        "rigorous": {
            "enclosure_basis": "gaussian_hermite",
            "pressure_pairing": "gradient_form",
            "pressure_term_lower_bound": str(pressure_lower),
            "pressure_term_upper_bound": str(pressure_term.upper),
            "pressure_term_width": str(pressure_term.upper - pressure_term.lower),
            "viscous_lower_bound": str(viscous_bound),
            "viscous_integral_upper_bound": str(viscous_total),
            "interior_viscous_integral_upper_bound": str(interior_viscous),
            "exterior_viscous_tail": str(tail.viscous_tail),
            "exterior_flux_tail": str(tail.flux_tail),
            "sensitivity_integral_upper_bound": str(sensitivity_integral),
            "margin": str(margin),
            "pressure_error_threshold": str(threshold),
            "margin_is_positive": bool(margin > 0),
            "cells_enclosed": cells,
            "leaf_evaluations": leaves,
            "refined_cells": refined_cells,
            "refine_budget_left": budget,
            "natural_product_width": str(natural_width),
            "correlated_product_width": str(correlated_width),
            "dependency_overestimation": str(natural_width - correlated_width),
            "pressure_gradient_variation_bound": str(gradient_variation_bound),
            "pressure_gradient_spread_max": str(gradient_spread_max),
            "tail_gamma": str(tail.gamma),
            "tail_distance": str(tail.distance),
            "exterior_pressure_note": (
                "The gradient-form pressure enclosure covers the box only.  The "
                "exterior remainder -3 int_ext |u| u . grad p dx is NOT "
                "enclosed: bounding it needs a pressure decay estimate that "
                "belongs to the same unproved P1G circle, so the conditional "
                "claim is restricted to the box contribution of the pressure "
                "integral.  The exterior viscous contribution IS enclosed and "
                "folded in."
            ),
        },
        "hypotheses": {
            "P1G_pressure_gradient_discretisation": {
                "statement": GRADIENT_PRESSURE_HYPOTHESIS,
                "threshold": str(threshold),
                "proved": False,
            }
        },
        "floating_point_reference": rate.as_dict(),
        "diagnostics": {
            "pressure_solver_residual": float(rate.poisson_residual),
            "pressure_gradient_variation_bound": float(gradient_variation_bound),
            "pressure_gradient_spread_max": float(gradient_spread_max),
            "finite_box_viscous_tail_in_margin": float(3 * nu * tail.viscous_tail),
            "finite_box_flux_tail_scale": float(3 * tail.flux_tail),
            "gaussian_velocity_tail_at_box": velocity_tail,
            "quadrature_enclosure_width": float(correlated_width),
            "dependency_overestimation": float(natural_width - correlated_width),
            "pressure_term_width": float(pressure_term.upper - pressure_term.lower),
        },
        "claims": [
            "Rigorous and unconditional: the stated enclosure of the box "
            "integral -3 int |u| u . grad p_h dx for the computed discrete "
            "pressure p_h, the upper bound on the viscous integral (interior "
            "enclosure plus Mills-ratio exterior tail), and the sensitivity "
            "integral, all in exact rational arithmetic.",
            "The discrete pressure gradient is the exact per-cell derivative "
            "of the multilinear interpolant of the nodal pressure -- rational "
            "linear combinations of nodal values, never an independent hull.",
            "Conditional on P1G: if the free-space pressure gradient is within "
            "the stated threshold of the interpolant gradient in L^infinity on "
            "the box, then the box contribution to J(u_0) is positive whenever "
            "the margin is.",
            "The advective contribution is exactly zero by an algebraic "
            "identity, not by cancellation of computed numbers.",
            "Nothing here is a singularity result and nothing here bears on "
            "the Clay problem.  A positive initial generation rate of the "
            "critical norm is a screening criterion, not evidence of blow-up.",
        ],
    }
    return GenerationCertificate(payload)


# --------------------------------------------------------------------------- #
# the independent checker                                                       #
# --------------------------------------------------------------------------- #


def verify_gradient_generation_certificate(
    payload: dict[str, object],
) -> dict[str, object]:
    """Independently re-check a gradient-form payload from its own contents."""
    try:
        return _verify(payload)
    except Exception as error:  # noqa: BLE001 -- adversarial input is a failure
        return {"verified": False, "failures": [f"payload rejected: {error}"]}


def _verify(payload: dict[str, object]) -> dict[str, object]:
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(
        payload.get("certificate") == GRADIENT_CERTIFICATE_KIND,
        "wrong certificate kind",
    )
    require(payload.get("schema_version") == 1, "unsupported schema version")
    require(float(payload["viscosity"]) > 0.0, "viscosity must be positive")
    require(payload.get("path") in _PATHS, "unknown pairing path")

    block = payload.get("rigorous")
    require(isinstance(block, dict), "missing rigorous block")
    if isinstance(block, dict):
        pressure_lower = Fraction(str(block["pressure_term_lower_bound"]))
        pressure_upper = Fraction(str(block["pressure_term_upper_bound"]))
        viscous = Fraction(str(block["viscous_lower_bound"]))
        sensitivity = Fraction(str(block["sensitivity_integral_upper_bound"]))
        margin = Fraction(str(block["margin"]))
        threshold = Fraction(str(block["pressure_error_threshold"]))
        require(
            pressure_lower <= pressure_upper,
            "the pressure enclosure endpoints are out of order",
        )
        require(
            Fraction(str(block["pressure_term_width"]))
            == pressure_upper - pressure_lower,
            "the pressure width does not match its endpoints",
        )
        require(viscous <= 0, "the viscous contribution must be non-positive")
        require(sensitivity >= 0, "the sensitivity bound must be nonnegative")
        require(
            margin == pressure_lower + viscous,
            "the margin does not equal P_lower + V_lower",
        )
        require(
            bool(block["margin_is_positive"]) == (margin > 0),
            "the margin flag contradicts the margin",
        )
        if sensitivity > 0:
            require(
                threshold * sensitivity == margin,
                "the threshold does not equal margin / sensitivity",
            )
        require(
            (margin > 0) == (threshold > 0),
            "a positive threshold requires a positive margin",
        )
        total = Fraction(str(block["viscous_integral_upper_bound"]))
        interior = Fraction(str(block["interior_viscous_integral_upper_bound"]))
        tail = Fraction(str(block["exterior_viscous_tail"]))
        require(
            total == interior + tail,
            "interior viscous + exterior tail must equal the total",
        )
        require(tail > 0, "a Gaussian datum must carry a positive exterior tail")
        natural = Fraction(str(block["natural_product_width"]))
        correlated = Fraction(str(block["correlated_product_width"]))
        require(
            Fraction(str(block["dependency_overestimation"]))
            == natural - correlated,
            "the dependency accounting does not reconcile",
        )

    hypotheses = payload.get("hypotheses")
    require(isinstance(hypotheses, dict), "missing hypotheses")
    if isinstance(hypotheses, dict):
        require(
            set(hypotheses) == {"P1G_pressure_gradient_discretisation"},
            "the payload must carry exactly the one named hypothesis",
        )
        entry = hypotheses["P1G_pressure_gradient_discretisation"]
        require(entry.get("proved") is False, "P1G claims to be proved; it is not")
        require(
            entry.get("statement") == GRADIENT_PRESSURE_HYPOTHESIS,
            "the P1G statement does not match the audited text",
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


# --------------------------------------------------------------------------- #
# the four-path comparison harness                                              #
# --------------------------------------------------------------------------- #


def _baseline_width_diagnostics(
    family,
    grid: AxisymmetricGrid,
    pressure: FloatArray,
    *,
    precision_bits: int,
    exp_terms: int,
) -> dict[str, Fraction]:
    r"""Measure the width accounting of the divergence-form pairing.

    Reproduces the pressure-times-flux product of
    :func:`~ns_certificate_lab.interval_gaussian.
    build_gaussian_generation_certificate` cell by cell and records the
    natural product-hull width against the correlated (centred-pressure)
    width of the same contribution.  Diagnostic only: the certified baseline
    margin comes from the audited builder, never from this sweep.
    """
    r_nodes = [Fraction(float(v)) for v in grid.r]
    z_nodes = [Fraction(float(v)) for v in grid.z]
    p_lo, p_hi = _cell_hull_from_nodes(pressure)
    natural_width = Fraction(0)
    correlated_width = Fraction(0)
    for i in range(grid.nr - 1):
        for j in range(grid.nz - 1):
            r_box = Interval(r_nodes[i], r_nodes[i + 1])
            z_box = Interval(z_nodes[j], z_nodes[j + 1])
            enclosure = cell_enclosure_gaussian(
                family, r_box, z_box, terms=exp_terms, bits=precision_bits
            )
            weight = (
                _TWO_PI_UPPER
                * r_box.upper
                * (r_nodes[i + 1] - r_nodes[i])
                * (z_nodes[j + 1] - z_nodes[j])
            )
            pressure_box = Interval(
                Fraction(float(p_lo[i, j])), Fraction(float(p_hi[i, j]))
            )
            product = pressure_box * enclosure.flux
            natural_width += 3 * weight * (product.upper - product.lower)
            centre = (pressure_box.lower + pressure_box.upper) / 2
            deviation = (pressure_box.upper - pressure_box.lower) / 2
            centred = enclosure.flux.scale(centre) + Interval(
                -deviation * enclosure.flux.magnitude,
                deviation * enclosure.flux.magnitude,
            )
            correlated_width += 3 * weight * (centred.upper - centred.lower)
    return {
        "natural_product_width": natural_width,
        "correlated_product_width": correlated_width,
        "dependency_overestimation": natural_width - correlated_width,
    }


def compare_certificate_paths(
    family,
    grid: AxisymmetricGrid,
    *,
    viscosity: float,
    precision_bits: int = 48,
    exp_terms: int = 28,
    max_depth: int = 2,
    refine_budget: int = 4096,
    include_baseline: bool = True,
) -> dict[str, object]:
    """Run the divergence-form baseline and the three gradient-form paths.

    Measurement harness for one candidate on one grid: each entry of the
    returned dictionary carries the path's payload plus float summaries of the
    margin and the separated widths, so refinement studies can tabulate them
    directly.  Nothing here changes what any single certificate proves.
    """
    from .interval_gaussian import build_gaussian_generation_certificate
    from .l3_generation import AxisymmetricPressureSolver

    report: dict[str, object] = {"paths": {}}
    paths: dict[str, dict[str, object]] = report["paths"]  # type: ignore[assignment]

    if include_baseline:
        baseline = build_gaussian_generation_certificate(
            family,
            grid,
            viscosity=viscosity,
            precision_bits=precision_bits,
            exp_terms=exp_terms,
            subdivisions=1,
        ).as_dict()
        field = gaussian_reference_field(family, grid)
        solver = AxisymmetricPressureSolver.build(grid)
        trace = solver.exterior_trace(field.u_r, field.u_theta, field.u_z)
        pressure = solver.solve(
            field.pressure_source(field.gradient()), trace=trace
        )
        widths = _baseline_width_diagnostics(
            family, grid, pressure,
            precision_bits=precision_bits, exp_terms=exp_terms,
        )
        rigorous = baseline["rigorous"]
        paths["divergence_hull"] = {
            "payload": baseline,
            "margin": float(Fraction(str(rigorous["margin"]))),
            "pressure_lower_bound": float(
                Fraction(str(rigorous["pressure_lower_bound_discrete"]))
            ),
            "viscous_lower_bound": float(
                Fraction(str(rigorous["viscous_lower_bound"]))
            ),
            "natural_product_width": float(widths["natural_product_width"]),
            "correlated_product_width": float(widths["correlated_product_width"]),
            "dependency_overestimation": float(
                widths["dependency_overestimation"]
            ),
        }

    for path in _PATHS:
        payload = build_gradient_generation_certificate(
            family,
            grid,
            viscosity=viscosity,
            path=path,
            precision_bits=precision_bits,
            exp_terms=exp_terms,
            max_depth=max_depth,
            refine_budget=refine_budget,
        ).as_dict()
        rigorous = payload["rigorous"]
        paths[f"gradient_{path}"] = {
            "payload": payload,
            "margin": float(Fraction(str(rigorous["margin"]))),
            "pressure_lower_bound": float(
                Fraction(str(rigorous["pressure_term_lower_bound"]))
            ),
            "pressure_upper_bound": float(
                Fraction(str(rigorous["pressure_term_upper_bound"]))
            ),
            "pressure_width": float(Fraction(str(rigorous["pressure_term_width"]))),
            "viscous_lower_bound": float(
                Fraction(str(rigorous["viscous_lower_bound"]))
            ),
            "natural_product_width": float(
                Fraction(str(rigorous["natural_product_width"]))
            ),
            "correlated_product_width": float(
                Fraction(str(rigorous["correlated_product_width"]))
            ),
            "dependency_overestimation": float(
                Fraction(str(rigorous["dependency_overestimation"]))
            ),
            "diagnostics": payload["diagnostics"],
        }
    return report
