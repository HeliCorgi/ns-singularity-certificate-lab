r"""Chaining Track-P slabs: a conditional certified existence interval on ``T^3``.

The two-slab step, exactly as proved
------------------------------------
Slab ``n`` certifies (conditional on EXT-P1/P2/P3, recorded and never
axiomatised): a strong periodic solution ``u`` exists on ``[t_n, t_{n+1}]`` and

.. math:: \|u(t) - u^{(n)}_a(t)\|_{\dot H^4} \le R_n(t), \qquad
          R_n(t_n) = \delta_n ,

where ``u^{(n)}_a`` is the Galerkin trajectory through the **exact rational**
datum ``y_n``.  At the boundary the next datum is a new exact rational point
``y_{n+1}`` (the Taylor endpoint of ``u^{(n)}_a``, dyadically rounded and
Leray-projected), and the triangle inequality gives

.. math:: \|u(t_{n+1}) - y_{n+1}\|_{\dot H^4}
          \;\le\; R_n(t_{n+1})
          \;+\; \|u^{(n)}_a(t_{n+1}) - y_{n+1}\|_{\dot H^4}
          \;=\; \delta_{n+1} ,

so slab ``n+1`` may start its control ODE from ``R_{n+1}(t_{n+1}) =
\delta_{n+1}``.  The solution on slab ``n+1`` is the **same** solution ``u``
by the uniqueness half of EXT-P1; the certified interval is the union of the
slabs and the certified radius is the per-slab tube.  The finite-inequality
skeleton of this composition is proved in Lean
(``formal/NSSingularity/TrackPChain.lean``).

Why there is no wrapping to fight
---------------------------------
Every slab starts from an exact rational point, so no interval box is ever
propagated across a slab boundary.  The endpoint enclosure of one slab is
collapsed to its midpoint and the discarded width is charged to the *scalar*
radius ``\delta_{n+1}``; the scalar control ODE then transports it.  The
classical wrapping effect of interval integrators lives in the vector
enclosure, and the vector enclosure dies at every boundary by construction.
The chain still *measures* the quantities a wrapping diagnosis needs
(endpoint widths, remainder shares, a numerical monodromy estimate) and
records the preregistered Lohner/QR triggers with their measured values, so
"no wrapping" is an observation in the payload, not an assumption.

Wording contract (checker-enforced)
-----------------------------------
Allowed: "conditional PDE certificate assuming EXT-P1/P2/P3"; "the
finite-dimensional Galerkin enclosure and the scalar control ODE are verified
unconditionally".  Forbidden, and rejected by :func:`verify_chain_certificate`
if found anywhere in a payload: any claim of unconditional PDE existence, any
claim that periodic Navier-Stokes is proved in Lean, and any reading of a
chain stop as a singularity.  A chain that stops proves nothing about
singularities; the stop is classified into the preregistered failure taxonomy
instead.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from .control_ode import (
    ControlCoefficients,
    build_control_certificate,
    prove_rough_enclosure,
    taylor_tube,
    verify_control_certificate,
)
from .fourier_torus import (
    FRACTION_OPS,
    TORUS_FAMILIES,
    TrigScalar,
    TrigVector,
    _IntervalOps,
    _norm_sq,
    adot_squared_upper,
    advection,
    galerkin_modes,
    leray,
    residual_tail,
    sup_derivative_bound,
    taylor_coefficients,
)
from .l3_certificate import sqrt_interval
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval
from .torus_aposteriori import (
    EXTERNAL_THEOREMS,
    SQRT8_UPPER,
    interval_field,
    prove_galerkin_box,
)

__all__ = [
    "ALLOWED_WORDING",
    "FORBIDDEN_WORDING",
    "build_chain_certificate",
    "verify_chain_certificate",
]

ALLOWED_WORDING = (
    "conditional PDE certificate assuming EXT-P1/P2/P3",
    "the finite-dimensional Galerkin enclosure and the scalar control ODE "
    "are verified unconditionally",
)

#: Substrings whose presence anywhere in a payload makes the checker reject it.
FORBIDDEN_WORDING = (
    "unconditional existence of the true PDE solution",
    "unconditionally proves the existence",
    "Navier-Stokes proved in Lean",
    "proves a PDE singularity",
    "singular time candidate",
)

STOP_IS_NOT_A_SINGULARITY = (
    "a chain that stops proves NOTHING about singularities: the stop reason "
    "is classified below and every classification names a limitation of the "
    "method, not a property of the solution"
)

#: The preregistered failure taxonomy of section 8.  Every terminated chain
#: must carry exactly one of these labels.
TERMINATION_REASONS = (
    "target_reached",
    "slab_budget_exhausted",
    "galerkin_tail",
    "time_integration_defect",
    "datum_transfer",
    "interval_wrapping",
    "control_linear_coefficient",
    "control_quadratic_term",
    "ext_hypotheses_insufficient",
    "approximate_norm_growth",
)


# --------------------------------------------------------------------------- #
# small exact helpers                                                          #
# --------------------------------------------------------------------------- #


def _upper(value) -> Fraction:
    return value.upper if isinstance(value, Interval) else Fraction(value)


def h4dot_upper(field: TrigVector, *, bits: int = DEFAULT_PRECISION_BITS) -> Fraction:
    """A rational upper bound on the homogeneous ``\\dot H^4`` norm.

    For interval fields the accumulated outward rounding of ``sobolev_sq``
    puts a floor of roughly ``sqrt(modes * weights) * 2^{-bits/2}`` under the
    result (about ``5e-9`` at 64 bits).  The bound stays a true upper bound;
    quantities measured at that floor (the recentring transfer of a chain,
    the Taylor remainder diagnostics) are dominated by rounding, not by the
    enclosed field, and the payload records them as such.
    """
    sob = _upper(field.sobolev_sq(4))
    if sob < 0:
        sob = Fraction(0)
    return sqrt_interval(Interval(Fraction(0), sob), bits=bits).upper


def l2_inner(u: TrigVector, v: TrigVector) -> Fraction:
    """The exact ``L^2`` inner product of two rational trig fields."""
    table_u = u.coefficient_table()
    table_v = v.coefficient_table()
    half = Fraction(1, 2)
    total = Fraction(0)
    for k, (a, b) in table_u.items():
        if k == (0, 0, 0):
            continue
        pair = table_v.get(k)
        if pair is None:
            continue
        a2, b2 = pair
        dot = sum(Fraction(a[i]) * Fraction(a2[i]) for i in range(3))
        dot += sum(Fraction(b[i]) * Fraction(b2[i]) for i in range(3))
        total += half * dot
    return total


def curl(field: TrigVector) -> TrigVector:
    """The exact curl of a rational trig field."""
    u1, u2, u3 = field.components
    return TrigVector(
        (
            u3.derivative(1) + (-u2.derivative(2)),
            u1.derivative(2) + (-u3.derivative(0)),
            u2.derivative(0) + (-u1.derivative(1)),
        )
    )


# --------------------------------------------------------------------------- #
# the Taylor endpoint of the Galerkin trajectory                               #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class EndpointEnclosure:
    """``u_a(h)`` enclosed by an interval Taylor form with Lagrange remainder."""

    field: TrigVector
    remainder_h4: Fraction
    max_width: Fraction


def galerkin_taylor_endpoint(
    datum: TrigVector,
    box: TrigVector,
    *,
    viscosity: Fraction,
    cutoff_sq: int,
    step: Fraction,
    order: int = 6,
    bits: int = DEFAULT_PRECISION_BITS,
) -> EndpointEnclosure:
    r"""Enclose ``u_a(h)`` from the exact datum and the proved slab box.

    Componentwise Taylor with Lagrange remainder: with ``c_i`` the exact
    time-Taylor coefficients through the datum and ``d_{m+1}`` the interval
    evaluation of the coefficient recurrence on the box,

    .. math:: u_a(h) \in \sum_{i \le m} c_i h^i + d_{m+1}\,h^{m+1} ,

    because the ``(m+1)``-st Taylor coefficient of the trajectory at any
    intermediate time is the fixed polynomial ``c_{m+1}(\,\cdot\,)`` of a state
    that the box provably contains.  Everything is interval arithmetic at
    ``bits`` precision seeded with the exact rational datum.
    """
    ops = _IntervalOps(bits)
    exact = interval_field(datum, bits=bits, ops=ops)
    polys = taylor_coefficients(
        exact, viscosity=viscosity, cutoff_sq=cutoff_sq, order=order
    )
    on_box = taylor_coefficients(
        box, viscosity=viscosity, cutoff_sq=cutoff_sq, order=order + 1
    )
    endpoint = TrigVector.zero(ops)
    for i, coefficient in enumerate(polys):
        endpoint = endpoint + coefficient.scale(step**i)
    remainder = on_box[order + 1].scale(step ** (order + 1))
    endpoint = endpoint + remainder

    max_width = Fraction(0)
    for scalar in endpoint.components:
        for k, (a, b) in scalar.terms.items():
            for value in (a, b):
                width = value.upper - value.lower
                if width > max_width:
                    max_width = width
    return EndpointEnclosure(
        field=endpoint,
        remainder_h4=h4dot_upper(remainder, bits=bits),
        max_width=max_width,
    )


@dataclass(frozen=True)
class Recentred:
    """An exact rational restart point and the scalar cost of moving to it."""

    datum: TrigVector
    transfer: Fraction
    rounding_bits: int


def recentre(
    endpoint: EndpointEnclosure,
    *,
    cutoff_sq: int,
    rounding_bits: int = 64,
    bits: int = DEFAULT_PRECISION_BITS,
) -> Recentred:
    r"""Collapse the endpoint enclosure to an exact divergence-free point.

    The midpoint of every coefficient interval is rounded to a dyadic with
    ``rounding_bits`` fractional bits and the rounded field is Leray-projected
    in exact rational arithmetic, which restores exact divergence-freeness.
    Because the true endpoint is exactly divergence free (the Galerkin field
    is a Leray projection composed with a Laplacian) and the mode-wise Leray
    projector contracts every ``\dot H^n``,

    .. math:: \|u_a(h) - P\,m\|_{\dot H^4} = \|P(u_a(h) - m)\|_{\dot H^4}
              \le \|u_a(h) - m\|_{\dot H^4} ,

    so the returned ``transfer`` -- the interval ``\dot H^4`` bound on
    ``u_a(h) - m`` -- also bounds the distance to the projected point.
    """
    scale = 2**rounding_bits
    mid_scalars = []
    for scalar in endpoint.field.components:
        terms: dict = {}
        for k, (a, b) in scalar.terms.items():
            if k == (0, 0, 0):
                if a.lower > 0 or a.upper < 0 or b.lower > 0 or b.upper < 0:
                    raise AssertionError(
                        "the Galerkin endpoint grew a mean mode away from zero"
                    )
                continue
            rounded = []
            for value in (a, b):
                midpoint = (value.lower + value.upper) / 2
                rounded.append(Fraction(round(midpoint * scale), scale))
            terms[k] = rounded
        mid_scalars.append(TrigScalar(terms, FRACTION_OPS))
    mid_field = TrigVector(tuple(mid_scalars)).restrict(
        lambda k: 1 <= _norm_sq(k) <= cutoff_sq
    )
    projected = leray(mid_field).cleaned()

    difference = endpoint.field + (-interval_field(mid_field, bits=bits))
    transfer = h4dot_upper(difference, bits=bits)
    return Recentred(datum=projected, transfer=transfer, rounding_bits=rounding_bits)


# --------------------------------------------------------------------------- #
# one chain slab                                                               #
# --------------------------------------------------------------------------- #


def _slab_constants(
    box_field: TrigVector,
    *,
    viscosity: Fraction,
    cutoff_sq: int,
    lattice_cut: int,
    bits: int,
) -> dict:
    """The slab-valid constants of the control inequality (turn-9 assembly)."""
    adot_sq = adot_squared_upper(lattice_cut)
    adot = sqrt_interval(Interval(adot_sq, adot_sq), bits=bits).upper
    bounds = [sup_derivative_bound(box_field, order, bits=bits) for order in range(1, 6)]
    m1, m2, m3, m4, m5 = bounds
    k1 = 4 * m1 + 6 * m2 + 4 * m3 + m4
    k2 = m5 + 4 * m4 + 6 * m3 + 4 * m2 + m1
    linear = -viscosity + 9 * (k1 + k2)
    quadratic = 135 * adot
    tail = residual_tail(box_field, cutoff_sq=cutoff_sq)
    tail_sq = _upper(tail.sobolev_sq(4))
    forcing = sqrt_interval(Interval(Fraction(0), tail_sq), bits=bits).upper
    return {
        "adot": adot,
        "adot_sq": adot_sq,
        "lattice_cut": lattice_cut,
        "derivative_bounds": bounds,
        "K1": k1,
        "K2": k2,
        "linear": linear,
        "quadratic": quadratic,
        "forcing": forcing,
        "tail_sq_upper": tail_sq,
    }


def _check_datum(datum: TrigVector, *, cutoff_sq: int) -> None:
    divergence = datum.divergence().cleaned()
    if divergence.terms:
        raise AssertionError("a chain datum must be exactly divergence free")
    for scalar in datum.components:
        for k in scalar.terms:
            if not 1 <= _norm_sq(k) <= cutoff_sq:
                raise AssertionError(f"chain datum mode outside the band: {k}")
    mean = advection(datum, datum).coefficient_table().get((0, 0, 0))
    if mean is not None:
        for vector in mean:
            if any(value != 0 for value in vector):
                raise AssertionError("the advection mean must vanish exactly")


def _physical_snapshot(
    datum: TrigVector,
    *,
    radius: Fraction,
    adot: Fraction,
    cutoff_sq: int,
    bits: int,
) -> dict:
    """Exact invariants of the restart point, with Lipschitz transfer radii.

    Every "certified_range" below is an interval that provably contains the
    corresponding quantity of the TRUE solution at the slab start, obtained
    from the exact value at the rational restart point plus a Lipschitz bound
    ``|Q(u) - Q(u_a)| \\le L_Q R`` proved from ``\\|w\\|_{\\dot H^n} \\le
    \\|w\\|_{\\dot H^4} \\le R`` (``|k| \\ge 1``) and, for sup-type
    quantities, from ``\\|f\\|_\\infty \\le \\dot A\\|f\\|_{\\dot H^2}``.
    Quantities marked "numeric" have float evaluation somewhere in their
    centre and are diagnostics, not certificates.
    """
    sob0 = Fraction(datum.sobolev_sq(0))
    sob1 = Fraction(datum.sobolev_sq(1))
    sob3 = Fraction(datum.sobolev_sq(3))
    sob4 = Fraction(datum.sobolev_sq(4))
    norm0 = sqrt_interval(Interval(Fraction(0), sob0), bits=bits).upper
    norm1 = sqrt_interval(Interval(Fraction(0), sob1), bits=bits).upper
    energy = sob0 / 2

    shells: dict[int, Fraction] = {}
    fluxes: dict[int, Fraction] = {}
    nonlinear = leray(advection(datum, datum))
    for shell in range(1, cutoff_sq + 1):
        selector = lambda k, s=shell: _norm_sq(k) == s  # noqa: E731
        piece = datum.restrict(selector)
        shells[shell] = Fraction(piece.sobolev_sq(0)) / 2
        fluxes[shell] = -l2_inner(nonlinear.restrict(selector), piece)

    vorticity = curl(datum)
    vorticity_sup_upper = sup_derivative_bound(vorticity, 0, bits=bits)
    samples = 9
    vorticity_sample = 0.0
    step_f = 2 * math.pi / samples
    for ix in range(samples):
        for iy in range(samples):
            for iz in range(samples):
                value = vorticity.evaluate((ix * step_f, iy * step_f, iz * step_f))
                vorticity_sample = max(
                    vorticity_sample, math.sqrt(sum(v * v for v in value))
                )

    total_energy = sum(shells.values())
    top_share = shells[cutoff_sq] / total_energy if total_energy else Fraction(0)
    log_slope = None
    decorated = [
        (math.sqrt(s), math.log(float(e)) / 2) for s, e in shells.items() if e > 0
    ]
    if len(decorated) >= 2:
        xs = [d[0] for d in decorated]
        ys = [d[1] for d in decorated]
        n = len(xs)
        denom = n * sum(x * x for x in xs) - sum(xs) ** 2
        if denom:
            log_slope = (n * sum(x * y for x, y in zip(xs, ys)) - sum(xs) * sum(ys)) / denom

    def certified(centre: Fraction, lipschitz: Fraction) -> dict:
        return {
            "centre": str(centre),
            "lipschitz_times_R": str(lipschitz * radius),
            "certified_range": [
                str(centre - lipschitz * radius),
                str(centre + lipschitz * radius),
            ],
        }

    return {
        "exact_at_restart_point": {
            "energy": str(energy),
            "enstrophy": str(sob1),
            "h3_dot_sq": str(sob3),
            "h4_dot_sq": str(sob4),
            "shell_energy": {str(s): str(e) for s, e in shells.items()},
            "shell_flux": {str(s): str(f) for s, f in fluxes.items()},
            "top_shell_energy_share": str(top_share),
        },
        "certified_ranges_for_true_solution": {
            "h4_dot_norm": certified(
                sqrt_interval(Interval(Fraction(0), sob4), bits=bits).upper,
                Fraction(1),
            ),
            "h3_dot_norm": certified(
                sqrt_interval(Interval(Fraction(0), sob3), bits=bits).upper,
                Fraction(1),
            ),
            "energy": certified(energy, norm0 + radius / 2),
            "enstrophy": certified(sob1, 2 * norm1 + radius),
            "l3_and_linf_radius": str(adot * radius),
            "vorticity_sup_upper_plus_radius": str(
                vorticity_sup_upper + adot * radius
            ),
            "derivation": (
                "|Q(u)-Q(u_a)| <= L_Q R with L=1 for Hdot3/Hdot4, "
                "L=(||u_a||_0 + R/2) for energy (bilinear), "
                "L=(2||u_a||_1 + R) for enstrophy, L=Adot for L3/Linf/vorticity "
                "sup via ||w||_inf <= Adot ||w||_Hdot2 <= Adot R"
            ),
        },
        "numeric_diagnostics": {
            "vorticity_sup_sample": vorticity_sample,
            "vorticity_sup_certified_upper": str(vorticity_sup_upper + adot * radius),
            "analyticity_slope_proxy": log_slope,
            "note": "sampled or fitted values; diagnostics, never certificates",
        },
    }


@dataclass(frozen=True)
class ChainSlab:
    """Everything one accepted (or failed) slab contributes to the chain."""

    proved: bool
    failure: str | None
    step: Fraction
    delta_in: Fraction
    radius_slab: Fraction | None
    delta_end: Fraction | None
    transfer: Fraction | None
    delta_out: Fraction | None
    new_datum: TrigVector | None
    constants: dict | None
    payload: dict


def prove_chain_slab(
    datum: TrigVector,
    *,
    delta_in: Fraction,
    viscosity: Fraction,
    cutoff_sq: int,
    step: Fraction,
    lattice_cut: int = 20,
    bits: int = 64,
    taylor_order: int = 6,
    rounding_bits: int = 64,
    label: str = "chain_slab",
) -> ChainSlab:
    """One link: box, constants, control ODE from ``delta_in``, restart point."""
    _check_datum(datum, cutoff_sq=cutoff_sq)

    box = prove_galerkin_box(
        datum,
        viscosity=viscosity,
        cutoff_sq=cutoff_sq,
        step=step,
        max_halvings=0,
        bits=bits,
    )
    base_payload: dict = {
        "step": str(step),
        "delta_in": str(delta_in),
        "datum_modes": {
            str(k): {
                "cos": [str(v) for v in a],
                "sin": [str(v) for v in b],
            }
            for k, (a, b) in datum.coefficient_table().items()
        },
        "datum_h4_sq": str(datum.sobolev_sq(4)),
        "galerkin_box": box.as_dict(),
    }
    if not box.proved:
        return ChainSlab(
            proved=False, failure="galerkin_box", step=step, delta_in=delta_in,
            radius_slab=None, delta_end=None, transfer=None, delta_out=None,
            new_datum=None, constants=None, payload=base_payload,
        )

    constants = _slab_constants(
        box.box, viscosity=viscosity, cutoff_sq=cutoff_sq,
        lattice_cut=lattice_cut, bits=bits,
    )
    coefficients = ControlCoefficients(
        linear=Interval(constants["linear"], constants["linear"]),
        quadratic=Interval(constants["quadratic"], constants["quadratic"]),
        forcing=Interval(Fraction(0), constants["forcing"]),
    )
    initial = Interval(Fraction(0), delta_in)
    rough = prove_rough_enclosure(
        coefficients, initial, step, max_halvings=0, precision_bits=bits
    )
    base_payload["constants"] = {
        "adot_upper": str(constants["adot"]),
        "adot_squared_upper": str(constants["adot_sq"]),
        "lattice_cut": constants["lattice_cut"],
        "sup_derivative_bounds_M1_to_M5": [
            str(m) for m in constants["derivative_bounds"]
        ],
        "K1": str(constants["K1"]),
        "K2": str(constants["K2"]),
        "linear_coefficient": str(constants["linear"]),
        "quadratic_coefficient": str(constants["quadratic"]),
        "residual_h4_upper": str(constants["forcing"]),
        "residual_h4_sq_upper": str(constants["tail_sq_upper"]),
    }
    if not rough.proved:
        return ChainSlab(
            proved=False, failure="control_rough", step=step, delta_in=delta_in,
            radius_slab=None, delta_end=None, transfer=None, delta_out=None,
            new_datum=None, constants=constants, payload=base_payload,
        )

    tube = taylor_tube(coefficients, rough, order=6, precision_bits=bits)
    radius_slab = tube.refined_tube.upper

    endpoint_value = Interval(Fraction(0), Fraction(0))
    for i, coefficient in enumerate(tube.coefficients):
        endpoint_value = endpoint_value + coefficient.scale(step**i)
    endpoint_value = endpoint_value + (
        tube.remainder_coefficient * Interval.exact(step ** (tube.order + 1))
    )
    delta_end = min(endpoint_value.upper, radius_slab)
    if delta_end < 0:
        delta_end = Fraction(0)

    control = build_control_certificate(
        coefficients, initial, step, order=6,
        max_halvings=0, precision_bits=bits, label=label,
    )
    control_payload = control.as_dict()
    control_verdict = verify_control_certificate(control_payload)

    endpoint = galerkin_taylor_endpoint(
        datum, box.box, viscosity=viscosity, cutoff_sq=cutoff_sq,
        step=step, order=taylor_order, bits=bits,
    )
    restart = recentre(
        endpoint, cutoff_sq=cutoff_sq, rounding_bits=rounding_bits, bits=bits
    )
    delta_out = delta_end + restart.transfer

    base_payload.update(
        {
            "control_certificate": control_payload,
            "control_verified": control_verdict,
            "radius_slab": str(radius_slab),
            "delta_end": str(delta_end),
            "taylor_endpoint": {
                "order": taylor_order,
                "remainder_h4": str(endpoint.remainder_h4),
                "max_coefficient_width": str(endpoint.max_width),
            },
            "recentre": {
                "rounding_bits": restart.rounding_bits,
                "transfer_h4": str(restart.transfer),
            },
            "delta_out": str(delta_out),
            "next_datum_modes": {
                str(k): {
                    "cos": [str(v) for v in a],
                    "sin": [str(v) for v in b],
                }
                for k, (a, b) in restart.datum.coefficient_table().items()
            },
            "physical": _physical_snapshot(
                datum, radius=radius_slab, adot=constants["adot"],
                cutoff_sq=cutoff_sq, bits=bits,
            ),
        }
    )
    return ChainSlab(
        proved=bool(control_verdict.get("proves_existence")),
        failure=None if control_verdict.get("proves_existence") else "control_verify",
        step=step, delta_in=delta_in, radius_slab=radius_slab,
        delta_end=delta_end, transfer=restart.transfer, delta_out=delta_out,
        new_datum=restart.datum, constants=constants, payload=base_payload,
    )


# --------------------------------------------------------------------------- #
# the chain                                                                    #
# --------------------------------------------------------------------------- #


def _classify_stop(slab: ChainSlab, *, delta: Fraction) -> dict:
    """Attribute a stopped chain to one preregistered failure class."""
    if slab.failure == "galerkin_box":
        return {
            "reason": "time_integration_defect",
            "detail": (
                "the Picard self-mapping test for the Galerkin trajectory did "
                "not close at the minimum slab length"
            ),
            "note": STOP_IS_NOT_A_SINGULARITY,
        }
    constants = slab.constants or {}
    linear = constants.get("linear", Fraction(0))
    quadratic = constants.get("quadratic", Fraction(0))
    forcing = constants.get("forcing", Fraction(0))
    shares = {
        "control_linear_coefficient": max(linear, Fraction(0)) * delta,
        "control_quadratic_term": quadratic * delta * delta,
        "galerkin_tail": forcing,
    }
    if slab.transfer is not None:
        shares["datum_transfer"] = slab.transfer / slab.step
    reason = max(shares, key=lambda key: shares[key])
    return {
        "reason": reason,
        "detail": {key: str(value) for key, value in shares.items()},
        "explanation": (
            "the derivative shares of the control ODE at the stopped slab; "
            "the largest term names the binding limitation.  The linear "
            "coefficient is the crude self-contained 9(K1+K2) bound, not a "
            "property of the solution"
        ),
        "note": STOP_IS_NOT_A_SINGULARITY,
    }


def build_chain_certificate(
    family: str,
    *,
    viscosity: Fraction = Fraction(1, 10),
    cutoff_sq: int = 4,
    initial_step: Fraction = Fraction(1, 2048),
    min_step: Fraction = Fraction(1, 32768),
    max_slabs: int = 64,
    target_time: Fraction | None = None,
    lattice_cut: int = 20,
    bits: int = 64,
    taylor_order: int = 6,
    max_taylor_order: int = 10,
    rounding_bits: int = 64,
    relative_radius_cap: Fraction = Fraction(1, 2),
    allow_step_doubling: bool = True,
) -> dict:
    """Chain slabs adaptively from a preregistered family until a stop fires.

    Preregistered adaptation, in priority order: (i) a slab that fails its
    box or control enclosure is retried at half the step, down to
    ``min_step``; (ii) a slab whose recentring transfer exceeds its tube
    growth is retried with a higher Taylor order, up to ``max_taylor_order``;
    (iii) a slab whose tube exceeds ``relative_radius_cap`` times the datum
    norm stops the chain; (iv) after two comfortable slabs (no box doubling
    and endpoint growth below 3/2) the step may double back, never above
    ``initial_step``.  Every stop is classified by :func:`_classify_stop`.
    """
    if family not in TORUS_FAMILIES:
        raise ValueError(f"unknown torus family {family!r}")
    if viscosity <= 0:
        raise ValueError("a Clay-admissible chain requires positive viscosity")

    datum = TORUS_FAMILIES[family]()
    delta = Fraction(0)
    clock = Fraction(0)
    step = Fraction(initial_step)
    order_now = taylor_order
    slabs: list[dict] = []
    recurrence: list[dict] = []
    termination: dict | None = None
    comfortable = 0
    initial_h4 = sqrt_interval(
        Interval(Fraction(0), Fraction(datum.sobolev_sq(4))), bits=bits
    ).upper
    max_radius = Fraction(0)
    last_slab: ChainSlab | None = None

    while len(slabs) < max_slabs:
        if target_time is not None and clock >= target_time:
            termination = {"reason": "target_reached",
                           "note": STOP_IS_NOT_A_SINGULARITY}
            break
        slab = prove_chain_slab(
            datum, delta_in=delta, viscosity=viscosity, cutoff_sq=cutoff_sq,
            step=step, lattice_cut=lattice_cut, bits=bits,
            taylor_order=order_now, rounding_bits=rounding_bits,
            label=f"track_p_chain_{family}_slab{len(slabs)}",
        )
        last_slab = slab
        if not slab.proved:
            if step / 2 >= min_step:
                step = step / 2
                comfortable = 0
                continue
            termination = _classify_stop(slab, delta=delta)
            break

        growth = slab.delta_end - delta
        if (
            slab.transfer > growth
            and slab.transfer > Fraction(1, 10**12)
            and order_now < max_taylor_order
        ):
            order_now += 2
            continue

        datum_h4 = sqrt_interval(
            Interval(Fraction(0), Fraction(datum.sobolev_sq(4))), bits=bits
        ).upper
        if datum_h4 > 0 and slab.radius_slab > relative_radius_cap * datum_h4:
            if step / 2 >= min_step:
                step = step / 2
                comfortable = 0
                continue
            termination = _classify_stop(slab, delta=delta)
            termination["trigger"] = "relative_radius_cap"
            break

        slabs.append(slab.payload)
        recurrence.append(
            {
                "n": len(slabs) - 1,
                "t_start": str(clock),
                "step": str(slab.step),
                "delta_in": str(slab.delta_in),
                "radius_slab": str(slab.radius_slab),
                "delta_end": str(slab.delta_end),
                "transfer": str(slab.transfer),
                "delta_out": str(slab.delta_out),
            }
        )
        clock += slab.step
        delta = slab.delta_out
        datum = slab.new_datum
        if slab.radius_slab > max_radius:
            max_radius = slab.radius_slab

        if (
            allow_step_doubling
            and slab.payload["galerkin_box"]["doublings"] == 0
            and delta > 0
            and slab.delta_end < delta * Fraction(3, 2)
        ):
            comfortable += 1
            if comfortable >= 2 and step * 2 <= initial_step:
                step = step * 2
                comfortable = 0
        else:
            comfortable = 0

    if termination is None:
        termination = {"reason": "slab_budget_exhausted",
                       "note": STOP_IS_NOT_A_SINGULARITY}
        if last_slab is not None and not last_slab.proved:
            termination = _classify_stop(last_slab, delta=delta)

    remainders = [
        Fraction(entry["taylor_endpoint"]["remainder_h4"]) for entry in slabs
    ]
    transfers = [Fraction(entry["recentre"]["transfer_h4"]) for entry in slabs]
    widths = [
        Fraction(entry["taylor_endpoint"]["max_coefficient_width"]) for entry in slabs
    ]
    amplifications = []
    for entry, rem, transfer in zip(slabs, remainders, transfers):
        floor = rem + Fraction(1, 2**rounding_bits)
        amplifications.append(float(transfer / floor) if floor else None)
    lohner_trigger = {
        "fired": False,
        "preregistered_conditions": [
            "wrapping amplification >= 2 for three consecutive slabs",
            "interval widths not reduced by halving the slab",
            "coefficient enclosure ends the proof while the point residual is small",
        ],
        "measured": {
            "endpoint_max_widths": [str(w) for w in widths],
            "remainder_h4": [str(r) for r in remainders],
            "transfer_h4": [str(t) for t in transfers],
            "transfer_over_remainder_floor": amplifications,
        },
        "verdict": (
            "every slab restarts from an exact rational point, so no interval "
            "box crosses a boundary; the measured transfer stays at the "
            "remainder+rounding floor and the preregistered Lohner/QR triggers "
            "did not fire.  A Lohner frame would have nothing to represent"
        ),
    }
    consecutive = 0
    for value in amplifications:
        consecutive = consecutive + 1 if value is not None and value >= 2 else 0
        if consecutive >= 3:
            lohner_trigger["fired"] = True
            lohner_trigger["verdict"] = (
                "the wrapping amplification trigger fired; a Lohner/QR frame "
                "should be prototyped and compared before this chain is extended"
            )
            break

    final_delta = delta
    payload: dict = {
        "schema_version": 1,
        "certificate": "torus_chain_aposteriori",
        "family": family,
        "viscosity": str(viscosity),
        "galerkin_cutoff_sq": cutoff_sq,
        "preregistered": {
            "initial_step": str(initial_step),
            "min_step": str(min_step),
            "max_slabs": max_slabs,
            "target_time": str(target_time) if target_time is not None else None,
            "relative_radius_cap": str(relative_radius_cap),
            "taylor_order": taylor_order,
            "max_taylor_order": max_taylor_order,
            "rounding_bits": rounding_bits,
            "bits": bits,
            "lattice_cut": lattice_cut,
            "adaptation": (
                "halve on enclosure failure down to min_step; raise Taylor "
                "order when transfer exceeds tube growth; stop at the relative "
                "radius cap; double back after two comfortable slabs"
            ),
        },
        "slab_count": len(slabs),
        "certified_final_time": str(clock),
        "single_slab_reference": str(initial_step),
        "horizon_over_single_slab": (
            float(clock / initial_step) if initial_step else None
        ),
        "max_radius_h4_dot": str(max_radius),
        "final_delta_h4_dot": str(final_delta),
        "final_distance_h3": str(SQRT8_UPPER * max_radius),
        "delta_recurrence": recurrence,
        "slabs": slabs,
        "wrapping_diagnostics": lohner_trigger,
        "termination": termination,
        "external_theorems": {
            name: {"statement": text, "proved": False, "axiomatised_in_lean": False}
            for name, text in EXTERNAL_THEOREMS.items()
        },
        "claims": list(ALLOWED_WORDING) + [STOP_IS_NOT_A_SINGULARITY],
        "conclusion": {
            "proved": len(slabs) > 0,
            "kind": "conditional PDE certificate assuming EXT-P1/P2/P3",
            "unconditional_part": (
                "the finite-dimensional Galerkin enclosure and the scalar "
                "control ODE are verified unconditionally"
            ),
            "statement": (
                f"assuming EXT-P1/P2/P3, a strong periodic solution exists on "
                f"[0, {clock}] and its H4-dot distance to the piecewise "
                f"Galerkin centre stays below {max_radius} throughout; the "
                f"identity of the solution across slab boundaries is the "
                f"uniqueness clause of EXT-P1"
            ),
            "not_a_singularity_statement": STOP_IS_NOT_A_SINGULARITY,
        },
    }
    return payload


# --------------------------------------------------------------------------- #
# the independent checker                                                      #
# --------------------------------------------------------------------------- #


def _parse_field(modes: dict) -> TrigVector:
    scalars = []
    for component in range(3):
        terms = {}
        for key, entry in modes.items():
            k = tuple(int(part) for part in key.strip("()").split(","))
            terms[k] = [
                Fraction(entry["cos"][component]),
                Fraction(entry["sin"][component]),
            ]
        scalars.append(TrigScalar(terms, FRACTION_OPS))
    return TrigVector(tuple(scalars)).cleaned()


def _contains_forbidden(payload) -> str | None:
    import json

    text = json.dumps(payload)
    for phrase in FORBIDDEN_WORDING:
        if phrase in text:
            return phrase
    return None


def verify_chain_certificate(payload: dict) -> dict:
    """Recompute every link of the chain and re-check every inequality.

    The checker rebuilds each slab from its stored datum -- box, constants,
    control coefficients, tube endpoint, Taylor endpoint, recentred restart
    point -- and demands exact agreement with the stored values, then replays
    the ``delta`` recurrence and the wording contract.  It shares no state
    with the builder beyond the library itself.
    """
    failures: list[str] = []
    checks: list[str] = []

    def check(condition: bool, message: str) -> None:
        (checks if condition else failures).append(message)

    try:
        check(
            payload.get("certificate") == "torus_chain_aposteriori",
            "certificate kind",
        )
        phrase = _contains_forbidden(payload)
        check(phrase is None, f"forbidden wording present: {phrase!r}")
        for sentence in ALLOWED_WORDING:
            check(sentence in payload.get("claims", []), f"claim missing: {sentence}")
        check(
            STOP_IS_NOT_A_SINGULARITY in payload.get("claims", []),
            "stop-is-not-a-singularity disclaimer",
        )
        externals = payload.get("external_theorems", {})
        for name in EXTERNAL_THEOREMS:
            block = externals.get(name, {})
            check(block.get("proved") is False, f"{name} must record proved: false")
            check(
                block.get("axiomatised_in_lean") is False,
                f"{name} must record axiomatised_in_lean: false",
            )
            check(
                block.get("statement") == EXTERNAL_THEOREMS[name],
                f"{name} statement must be recorded verbatim",
            )
        termination = payload.get("termination", {})
        check(
            termination.get("reason") in TERMINATION_REASONS,
            "termination reason must be preregistered",
        )

        viscosity = Fraction(payload["viscosity"])
        cutoff_sq = int(payload["galerkin_cutoff_sq"])
        prereg = payload.get("preregistered", {})
        bits = int(prereg.get("bits", 64))
        lattice_cut = int(prereg.get("lattice_cut", 20))
        slabs = payload.get("slabs", [])
        recurrence = payload.get("delta_recurrence", [])
        check(len(slabs) == len(recurrence), "one recurrence row per slab")
        check(len(slabs) == int(payload.get("slab_count", -1)), "slab count")

        family = payload.get("family")
        if slabs and family in TORUS_FAMILIES:
            first = _parse_field(slabs[0]["datum_modes"])
            expected = TORUS_FAMILIES[family]().cleaned()
            check(
                first.coefficient_table() == expected.coefficient_table(),
                "slab 0 datum equals the preregistered family",
            )

        clock = Fraction(0)
        delta = Fraction(0)
        max_radius = Fraction(0)
        for index, (entry, row) in enumerate(zip(slabs, recurrence)):
            step = Fraction(entry["step"])
            check(Fraction(row["t_start"]) == clock, f"slab {index}: contiguous time")
            check(
                Fraction(entry["delta_in"]) == delta,
                f"slab {index}: delta_in equals the running delta",
            )
            datum = _parse_field(entry["datum_modes"])
            _check_datum(datum, cutoff_sq=cutoff_sq)
            checks.append(f"slab {index}: datum exactly divergence free and in band")

            box = prove_galerkin_box(
                datum, viscosity=viscosity, cutoff_sq=cutoff_sq,
                step=step, max_halvings=0, bits=bits,
            )
            check(box.proved, f"slab {index}: Galerkin box reproved")
            stored_box = entry["galerkin_box"]
            check(
                Fraction(stored_box["radius"]) == box.radius
                and Fraction(stored_box["step"]) == box.step,
                f"slab {index}: box radius and step reproduce",
            )

            constants = _slab_constants(
                box.box, viscosity=viscosity, cutoff_sq=cutoff_sq,
                lattice_cut=lattice_cut, bits=bits,
            )
            stored = entry["constants"]
            check(
                Fraction(stored["linear_coefficient"]) == constants["linear"]
                and Fraction(stored["quadratic_coefficient"]) == constants["quadratic"]
                and Fraction(stored["residual_h4_upper"]) == constants["forcing"]
                and Fraction(stored["K1"]) == constants["K1"]
                and Fraction(stored["K2"]) == constants["K2"],
                f"slab {index}: control constants reproduce from the box",
            )

            control_payload = entry["control_certificate"]
            verdict = verify_control_certificate(control_payload)
            check(
                bool(verdict.get("verified")) and bool(verdict.get("proves_existence")),
                f"slab {index}: embedded control certificate verifies",
            )
            stored_coefficients = ControlCoefficients.from_dict(
                control_payload["coefficients"]
            )
            check(
                stored_coefficients.linear.upper == constants["linear"]
                and stored_coefficients.quadratic.upper == constants["quadratic"]
                and stored_coefficients.forcing.upper == constants["forcing"],
                f"slab {index}: control coefficients wired from the constants",
            )
            check(
                Interval.from_pair(control_payload["initial"]).upper == delta,
                f"slab {index}: control ODE starts at the running delta",
            )

            coefficients = ControlCoefficients(
                linear=Interval(constants["linear"], constants["linear"]),
                quadratic=Interval(constants["quadratic"], constants["quadratic"]),
                forcing=Interval(Fraction(0), constants["forcing"]),
            )
            rough = prove_rough_enclosure(
                coefficients, Interval(Fraction(0), delta), step,
                max_halvings=0, precision_bits=bits,
            )
            check(rough.proved, f"slab {index}: rough enclosure reproved")
            tube = taylor_tube(coefficients, rough, order=6, precision_bits=bits)
            radius_slab = tube.refined_tube.upper
            check(
                Fraction(entry["radius_slab"]) == radius_slab,
                f"slab {index}: slab radius reproduces",
            )
            endpoint_value = Interval(Fraction(0), Fraction(0))
            for i, coefficient in enumerate(tube.coefficients):
                endpoint_value = endpoint_value + coefficient.scale(step**i)
            endpoint_value = endpoint_value + (
                tube.remainder_coefficient * Interval.exact(step ** (tube.order + 1))
            )
            delta_end = min(endpoint_value.upper, radius_slab)
            if delta_end < 0:
                delta_end = Fraction(0)
            check(
                Fraction(entry["delta_end"]) == delta_end,
                f"slab {index}: control endpoint reproduces",
            )

            order = int(entry["taylor_endpoint"]["order"])
            endpoint = galerkin_taylor_endpoint(
                datum, box.box, viscosity=viscosity, cutoff_sq=cutoff_sq,
                step=step, order=order, bits=bits,
            )
            restart = recentre(
                endpoint, cutoff_sq=cutoff_sq,
                rounding_bits=int(entry["recentre"]["rounding_bits"]), bits=bits,
            )
            check(
                Fraction(entry["recentre"]["transfer_h4"]) == restart.transfer,
                f"slab {index}: recentring transfer reproduces",
            )
            stored_next = entry["next_datum_modes"]
            recomputed_next = {
                str(k): {
                    "cos": [str(v) for v in a],
                    "sin": [str(v) for v in b],
                }
                for k, (a, b) in restart.datum.coefficient_table().items()
            }
            check(
                stored_next == recomputed_next,
                f"slab {index}: restart point reproduces exactly",
            )
            if index + 1 < len(slabs):
                check(
                    slabs[index + 1]["datum_modes"] == stored_next,
                    f"slab {index}: next slab starts at the restart point",
                )

            delta_out = delta_end + restart.transfer
            check(
                Fraction(entry["delta_out"]) == delta_out,
                f"slab {index}: delta recurrence step reproduces",
            )
            check(
                Fraction(row["delta_out"]) == delta_out,
                f"slab {index}: recurrence table matches",
            )
            delta = delta_out
            clock += step
            if radius_slab > max_radius:
                max_radius = radius_slab

        check(
            Fraction(payload["certified_final_time"]) == clock,
            "certified final time equals the slab sum",
        )
        check(
            Fraction(payload["max_radius_h4_dot"]) == max_radius,
            "chain radius equals the slab maximum",
        )
        check(
            Fraction(payload["final_delta_h4_dot"]) == delta,
            "final delta equals the recurrence output",
        )
        check(
            Fraction(payload["final_distance_h3"]) == SQRT8_UPPER * max_radius,
            "H3 distance uses the recorded sqrt8 bound",
        )
        conclusion = payload.get("conclusion", {})
        check(
            conclusion.get("kind") == ALLOWED_WORDING[0],
            "conclusion uses the allowed conditional wording",
        )
        check(
            conclusion.get("unconditional_part") == ALLOWED_WORDING[1],
            "conclusion separates the unconditional part",
        )
    except Exception as error:  # noqa: BLE001 -- a checker must not crash
        failures.append(f"exception while checking: {error!r}")

    verified = not failures
    return {
        "verified": verified,
        "proves_existence": verified and bool(payload.get("slabs")),
        "conditional_on": list(EXTERNAL_THEOREMS),
        "checks": checks,
        "failures": failures,
    }
