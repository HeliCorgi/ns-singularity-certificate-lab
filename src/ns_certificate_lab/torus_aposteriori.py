r"""The Track-P slab certificate: a rigorous a posteriori bound on ``T^3``.

What one certificate asserts
----------------------------
Take a finite Fourier datum ``u_0`` (rational coefficients, mean zero,
divergence free — hence ``C^\infty``), a rational viscosity ``\nu > 0`` and a
Galerkin band ``G = \{1 \le |k|^2 \le K^2\}``.  Let ``u_a`` be the Galerkin
trajectory through ``u_0``.  The certificate proves, in exact rational
arithmetic:

1. an enclosure ``B`` of ``u_a`` over the whole slab ``[0, h]`` by the Picard
   self-mapping test on the finite-dimensional quadratic Galerkin field — a
   theorem, not an estimate (see ``PICARD_SELF_MAPPING`` in
   :mod:`ns_certificate_lab.control_ode`);
2. rational bounds, valid on the whole slab, for every constant of the ``H^4``
   control inequality

   .. math::
      \tfrac{d^+}{dt}R \le \bigl(-\nu + 9(K_1+K_2)\bigr)R + 135\dot A\,R^2
        + \|e\|_{\dot H^4},\qquad R(0) = 0 ,

   where ``e = (I-P_G)P(u_a\cdot\nabla u_a)`` is the **exact** continuous
   Navier–Stokes residual of the trajectory (the Galerkin tail, a finite
   trigonometric polynomial — the periodic ``HS-5`` gap closes by
   construction), ``K_1, K_2`` are ``\ell^1`` bounds on derivatives of ``u_a``
   up to fifth order, and ``\dot A^2 \ge \sum_{k\ne0}|k|^{-4}`` is the one
   lattice constant behind the embedding
   ``\|f\|_\infty \le \dot A\|f\|_{\dot H^2}``;
3. a bound ``R(t)`` on the whole slab for any solution of that scalar Riccati
   inequality, by the control-ODE machinery of
   :mod:`ns_certificate_lab.control_ode`.

Combined with the three named classical external theorems below — recorded
faithfully in the payload and **never** inserted into Lean as axioms — the
conclusion is:

    **a strong periodic Navier–Stokes solution exists on all of ``[0, h]``
    and satisfies** ``\|u(t) - u_a(t)\|_{\dot H^4} \le R(t)``, hence
    ``\|u(t) - u_a(t)\|_{H^3} \le \sqrt8\,R(t)``.

This certifies *regularity near a computed trajectory*.  It is the opposite of
a singularity certificate, it is the correct first rigorous step, and it says
nothing about the Clay problem.

Why the constants close at ``n = 4`` and not ``n = 3``
------------------------------------------------------
Every commutator term of the ``H^4`` estimate is bounded using only the
Cauchy–Schwarz embedding ``\|f\|_\infty \le \dot A\|f\|_{\dot H^2}``: for the
cubic term ``\langle w\cdot\nabla w, w\rangle_4``, Leibniz factors with one or
two derivatives take the sup bound directly, and factors with three or four
derivatives put the sup on the *other* factor.  At ``n = 3`` the middle term
needs a Gagliardo–Nirenberg ``L^6\times L^3`` product, which is exactly what
the sharp Kato machinery of Morosi–Pizzocchero supplies (their ``G_3 = 0.438``;
see :doc:`../../docs/research_notes/a_posteriori_frameworks`).  The crude
self-contained constants used here cost slab length, not validity; the tradeoff
is recorded, not hidden.

The full derivation, constant by constant, is in
:doc:`../../docs/research_notes/track_p_periodic`.

External theorems (named, faithful, not axiomatised)
----------------------------------------------------
``EXT-P1`` — local existence and uniqueness of a strong solution in
``H^4(T^3)`` for mean-zero divergence-free ``H^4`` data (Kato's quasi-linear
theory; Temam, *Navier–Stokes Equations*, Ch. III, for periodic ``H^m``,
``m \ge 2``).
``EXT-P2`` — that solution has the regularity used by the energy estimate and
the Dini comparison: ``u \in C([0,T]; H^4) \cap C^1([0,T]; H^2)``, and
``t \mapsto \|w(t)\|_{\dot H^4}`` has a right Dini derivative bounded by the
estimate (standard for strong solutions).
``EXT-P3`` — the ``H^4`` continuation criterion: the strong solution extends
for as long as ``\|u\|_{H^4}`` stays finite, which the distance bound plus the
trajectory enclosure guarantees on the slab.

Each is classical; none is proved here; the payload carries all three with
``proved: false`` and the checker rejects a payload claiming otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

from .control_ode import (
    ControlCoefficients,
    build_control_certificate,
    verify_control_certificate,
)
from .fourier_torus import (
    TORUS_FAMILIES,
    TrigScalar,
    TrigVector,
    _IntervalOps,
    adot_squared_upper,
    advection,
    galerkin_modes,
    galerkin_rhs,
    leray,
    residual_tail,
    sup_derivative_bound,
)
from .l3_certificate import sqrt_interval
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

__all__ = [
    "EXTERNAL_THEOREMS",
    "SQRT8_UPPER",
    "TorusSlabCertificate",
    "build_torus_slab_certificate",
    "interval_field",
    "prove_galerkin_box",
    "verify_torus_slab_certificate",
]

#: ``\|w\|_{H^3} \le \sqrt8\,\|w\|_{\dot H^4}`` on mean-zero fields, because
#: ``(1+|k|^2)^3 \le (2|k|^2)^3 = 8|k|^6 \le 8|k|^8`` for ``|k| \ge 1``.
#: ``2828428/1000000 > 2\sqrt2``.
SQRT8_UPPER = Fraction(2828428, 1000000)

EXTERNAL_THEOREMS = {
    "EXT-P1": (
        "Local existence and uniqueness of a strong solution of the "
        "incompressible Navier-Stokes equations on the torus T^3 in H^4, for "
        "mean-zero divergence-free H^4 initial data and fixed positive "
        "viscosity.  Classical: Kato's quasi-linear theory; Temam, "
        "Navier-Stokes Equations, Ch. III (periodic H^m theory, m >= 2)."
    ),
    "EXT-P2": (
        "The strong solution of EXT-P1 satisfies u in C([0,T]; H^4) cap "
        "C^1([0,T]; H^2), the H^4 energy estimate for the difference with a "
        "smooth trajectory holds, and t -> ||w(t)||_{H^4-dot} admits a right "
        "Dini derivative bounded by that estimate.  Standard regularity of "
        "strong solutions."
    ),
    "EXT-P3": (
        "H^4 continuation: the strong solution extends in time for as long as "
        "its H^4 norm remains finite.  Classical continuation criterion for "
        "periodic strong solutions."
    ),
}


def interval_field(field: TrigVector, *, bits: int = DEFAULT_PRECISION_BITS,
                   ops: _IntervalOps | None = None) -> TrigVector:
    """The same field with degenerate interval coefficients."""
    interval_ops = _IntervalOps(bits) if ops is None else ops
    scalars = []
    for scalar in field.components:
        terms = {
            k: [Interval(Fraction(a), Fraction(a)), Interval(Fraction(b), Fraction(b))]
            for k, (a, b) in scalar.terms.items()
        }
        scalars.append(TrigScalar(terms, interval_ops))
    return TrigVector(tuple(scalars))


@dataclass(frozen=True)
class GalerkinBox:
    """A proved enclosure of the Galerkin trajectory over one slab."""

    proved: bool
    step: Fraction
    radius: Fraction
    doublings: int
    halvings: int
    box: TrigVector
    worst_margin: Fraction

    def as_dict(self) -> dict[str, object]:
        return {
            "proved": self.proved,
            "step": str(self.step),
            "radius": str(self.radius),
            "doublings": self.doublings,
            "halvings": self.halvings,
            "worst_margin": str(self.worst_margin),
        }


def prove_galerkin_box(
    datum: TrigVector,
    *,
    viscosity: Fraction,
    cutoff_sq: int,
    step: Fraction,
    max_doublings: int = 24,
    max_halvings: int = 8,
    bits: int = DEFAULT_PRECISION_BITS,
) -> GalerkinBox:
    r"""The Picard self-mapping test for the finite-dimensional Galerkin system.

    The box ``B`` is a coefficient-wise interval of uniform radius ``\rho``
    around the datum, over **every** canonical mode of the Galerkin band — also
    the modes on which the datum vanishes, because the image ``F(B)`` populates
    them and a box that omitted them could not contain its own image.  The test

    .. math::  Y_0 + [0, h]\,F(B) \subseteq B \quad\text{componentwise}

    is checked in exact rational interval arithmetic.  If it holds, the Picard
    operator maps ``C([0,h], B)`` into itself and is a contraction in the
    Bielecki norm (the field is a polynomial, hence Lipschitz on ``B``), so the
    Galerkin trajectory exists on the whole slab, is unique, and stays in ``B``.
    That is a theorem about a finite-dimensional ODE with no hypothesis left
    open — the same statement, verbatim, as ``PICARD_SELF_MAPPING`` in the
    scalar control module, applied coefficient-wise.

    The search doubles ``\rho`` on failure and then halves ``h``; the *test* is
    the theorem and the search is only a heuristic, so a failed search is
    reported as a failure rather than papered over.
    """
    ops = _IntervalOps(bits)
    modes = galerkin_modes(cutoff_sq)
    exact = interval_field(datum, bits=bits, ops=ops)

    rhs0 = galerkin_rhs(exact, viscosity=viscosity, cutoff_sq=cutoff_sq)
    scale = max(
        (
            max(pair[0].magnitude, pair[1].magnitude)
            for scalar in rhs0.components
            for pair in scalar.terms.values()
        ),
        default=Fraction(0),
    )
    step_now = Fraction(step)
    for halving in range(max_halvings + 1):
        radius = 2 * step_now * scale + Fraction(1, 1024)
        for doubling in range(max_doublings + 1):
            box_scalars = []
            for component in range(3):
                terms = {}
                base = exact.components[component].terms
                for k in modes:
                    centre = base.get(k, [ops.zero, ops.zero])
                    terms[k] = [
                        Interval(
                            centre[0].lower - radius, centre[0].upper + radius
                        ),
                        Interval(
                            centre[1].lower - radius, centre[1].upper + radius
                        ),
                    ]
                box_scalars.append(TrigScalar(terms, ops))
            box = TrigVector(tuple(box_scalars))
            image = galerkin_rhs(box, viscosity=viscosity, cutoff_sq=cutoff_sq)
            time_box = Interval(Fraction(0), step_now)
            worst = None
            contained = True
            image_table = image.coefficient_table()
            for component in range(3):
                base = exact.components[component].terms
                for k in modes:
                    centre = base.get(k, [ops.zero, ops.zero])
                    field_pair = image_table.get(k)
                    for slot in range(2):
                        increment = (
                            (time_box * field_pair[slot][component]).round_outward(bits)
                            if field_pair is not None
                            else ops.zero
                        )
                        reached = centre[slot] + increment
                        box_int = box.components[component].terms[k][slot]
                        margin = min(
                            reached.lower - box_int.lower,
                            box_int.upper - reached.upper,
                        )
                        worst = margin if worst is None else min(worst, margin)
                        if margin < 0:
                            contained = False
            if contained:
                return GalerkinBox(
                    proved=True, step=step_now, radius=radius,
                    doublings=doubling, halvings=halving, box=box,
                    worst_margin=worst if worst is not None else Fraction(0),
                )
            radius *= 2
        step_now /= 2
    return GalerkinBox(
        proved=False, step=step_now, radius=radius, doublings=max_doublings,
        halvings=max_halvings, box=exact,
        worst_margin=worst if worst is not None else Fraction(-1),
    )


@dataclass(frozen=True)
class TorusSlabCertificate:
    """A machine-checkable record of one Track-P slab."""

    payload: dict[str, object]

    def as_dict(self) -> dict[str, object]:
        return self.payload


def build_torus_slab_certificate(
    family: str,
    *,
    viscosity: Fraction = Fraction(1, 10),
    cutoff_sq: int = 4,
    step: Fraction = Fraction(1, 2048),
    lattice_cut: int = 20,
    bits: int = 64,
) -> TorusSlabCertificate:
    """Assemble one certificate: datum checks, slab box, constants, control ODE."""
    if family not in TORUS_FAMILIES:
        raise ValueError(f"unknown torus family {family!r}")
    if viscosity <= 0:
        raise ValueError("a Clay-admissible run requires fixed positive viscosity")
    datum = TORUS_FAMILIES[family]()

    # -- exact datum checks (rational arithmetic, no rounding) ---------------- #
    divergence = datum.divergence().cleaned()
    if divergence.terms:
        raise AssertionError("the datum is not exactly divergence free")
    outside = [
        k
        for scalar in datum.components
        for k in scalar.terms
        if not 1 <= (k[0] ** 2 + k[1] ** 2 + k[2] ** 2) <= cutoff_sq
    ]
    if outside:
        raise ValueError(f"datum modes outside the Galerkin band: {outside}")
    mean_advection = advection(datum, datum).coefficient_table().get((0, 0, 0))
    if mean_advection is not None:
        for vector in mean_advection:
            if any(value != 0 for value in vector):
                raise AssertionError("the advection mean must vanish exactly")

    # -- the slab enclosure of the trajectory --------------------------------- #
    box = prove_galerkin_box(
        datum, viscosity=viscosity, cutoff_sq=cutoff_sq, step=step, bits=bits
    )

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate": "torus_slab_aposteriori",
        "family": family,
        "viscosity": str(viscosity),
        "galerkin_cutoff_sq": cutoff_sq,
        "datum": {
            "modes": {
                str(k): {
                    "cos": [str(value) for value in a],
                    "sin": [str(value) for value in b],
                }
                for k, (a, b) in datum.coefficient_table().items()
            },
            "mean_zero": True,
            "divergence_free_exact": True,
            "smoothness": (
                "a finite trigonometric polynomial, hence C-infinity "
                "(Lean: TrackPFourier.lean)"
            ),
            "h4_norm_sq": str(datum.sobolev_sq(4)),
            "energy": str(datum.energy()),
            "enstrophy": str(datum.enstrophy()),
        },
        "distinction": {
            "fixed_band_trajectory": (
                "excluded by the finite-mode no-go (Track F): a trajectory that "
                "REMAINS in the band for all time cannot blow up"
            ),
            "finite_band_datum": (
                "NOT excluded: the true solution leaves the band immediately "
                "(the residual below is nonzero); Track P bounds its distance "
                "to the band-limited trajectory instead"
            ),
        },
        "galerkin_box": box.as_dict(),
    }

    if not box.proved:
        payload["conclusion"] = {
            "proved": False,
            "reason": "the Picard self-mapping test for the Galerkin trajectory "
                      "did not close at the requested slab length",
        }
        payload["claims"] = list(_CLAIMS)
        payload["external_theorems"] = _external_block()
        return TorusSlabCertificate(payload)

    # -- slab-valid constants of the control inequality ------------------------ #
    adot_sq = adot_squared_upper(lattice_cut)
    adot = sqrt_interval(Interval(adot_sq, adot_sq), bits=bits).upper
    derivative_bounds = [
        sup_derivative_bound(box.box, order, bits=bits) for order in range(1, 6)
    ]
    m1, m2, m3, m4, m5 = derivative_bounds
    k1 = 4 * m1 + 6 * m2 + 4 * m3 + m4
    k2 = m5 + 4 * m4 + 6 * m3 + 4 * m2 + m1
    linear = -viscosity + 9 * (k1 + k2)
    quadratic = 135 * adot

    tail = residual_tail(box.box, cutoff_sq=cutoff_sq)
    tail_sq = tail.sobolev_sq(4)
    forcing = sqrt_interval(
        Interval(Fraction(0), tail_sq.upper), bits=bits
    ).upper

    control = build_control_certificate(
        ControlCoefficients(
            linear=Interval(linear, linear),
            quadratic=Interval(quadratic, quadratic),
            forcing=Interval(Fraction(0), forcing),
        ),
        Interval(Fraction(0), Fraction(0)),
        box.step,
        order=6,
        precision_bits=bits,
        label=f"track_p_{family}",
    )
    control_payload = control.as_dict()
    control_verdict = verify_control_certificate(control_payload)
    tube_block = control_payload.get("tube")
    distance = (
        Interval.from_pair(tube_block["refined_tube"]).upper
        if isinstance(tube_block, dict)
        else None
    )

    payload["constants"] = {
        "adot_squared_upper": str(adot_sq),
        "adot_upper": str(adot),
        "lattice_cut": lattice_cut,
        "sup_derivative_bounds_M1_to_M5": [str(m) for m in derivative_bounds],
        "K1": str(k1),
        "K2": str(k2),
        "linear_coefficient": str(linear),
        "quadratic_coefficient": str(quadratic),
        "residual_h4_sq_upper": str(tail_sq.upper),
        "residual_h4_upper": str(forcing),
        "derivation": (
            "d+R/dt <= (-nu + 9(K1+K2)) R + 135 Adot R^2 + ||e||_{H4-dot}, "
            "R(0) = 0; K1 = 4M1+6M2+4M3+M4, K2 = M5+4M4+6M3+4M2+M1; the 9 is "
            "sqrt(81) from the multinomial weights of the Hdot4 form, the "
            "binomials are multi-Vandermonde, and 135 = 9 x 15 with "
            "15 = sum_{j=1..4} C(4,j).  Full derivation: "
            "docs/research_notes/track_p_periodic.md"
        ),
    }
    payload["residual_note"] = (
        "the residual is the EXACT continuous Navier-Stokes residual of the "
        "trajectory -- a finite trigonometric polynomial computed in rational "
        "arithmetic with no spatial discretisation; the periodic HS-5 gap "
        "closes by construction"
    )
    payload["control_certificate"] = control_payload
    payload["control_verified"] = control_verdict
    payload["conclusion"] = {
        "proved": bool(
            control_verdict.get("proves_existence") and distance is not None
        ),
        "conditional_on": ["EXT-P1", "EXT-P2", "EXT-P3"],
        "distance_h4_dot": str(distance) if distance is not None else None,
        "sqrt8_upper": str(SQRT8_UPPER),
        "distance_h3": str(SQRT8_UPPER * distance) if distance is not None else None,
        "statement": (
            "a strong periodic Navier-Stokes solution exists on the whole slab "
            "and its H4-dot distance to the certified Galerkin trajectory is at "
            "most distance_h4_dot at every instant; the H^3 distance is at most "
            "distance_h3"
        ),
    }
    payload["external_theorems"] = _external_block()
    payload["claims"] = list(_CLAIMS)
    return TorusSlabCertificate(payload)


def _external_block() -> dict[str, object]:
    return {
        name: {"statement": text, "proved": False, "axiomatised_in_lean": False}
        for name, text in EXTERNAL_THEOREMS.items()
    }


_CLAIMS = [
    "Rigorous and unconditional: the datum checks, the Picard enclosure of the "
    "finite-dimensional Galerkin trajectory over the slab, every constant of "
    "the control inequality as a bound valid on the whole slab, and the scalar "
    "control-ODE tube.",
    "Conditional on EXT-P1, EXT-P2, EXT-P3 (classical, named, never "
    "axiomatised in Lean): a strong periodic solution exists on the slab with "
    "the stated H4-dot and H3 distances to the trajectory.",
    "This certifies regularity near a computed trajectory.  It is not a "
    "singularity result, and nothing here bears on the Clay problem.",
]


def verify_torus_slab_certificate(payload: dict[str, object]) -> dict[str, object]:
    """Independently re-check a Track-P payload from its own numbers.

    The checker re-derives every relation among the recorded numbers — the
    coefficient assembly ``K_1, K_2 \\to`` linear coefficient, the ``135\\dot A``
    quadratic coefficient, the ``\\sqrt8`` conversion — and re-verifies the
    embedded control certificate with its own independent checker.  It does
    **not** re-run the interval convolutions that produced ``M_j`` and the
    residual bound; those are the generator's work, and the division of labour
    is recorded here so nobody mistakes consistency for recomputation.
    """
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
        payload.get("certificate") == "torus_slab_aposteriori",
        "wrong certificate kind",
    )
    require(payload.get("schema_version") == 1, "unsupported schema version")
    require(Fraction(str(payload["viscosity"])) > 0, "viscosity must be positive")

    datum = payload.get("datum")
    require(isinstance(datum, dict), "missing datum block")
    if isinstance(datum, dict):
        require(bool(datum.get("mean_zero")), "the datum must be mean zero")
        require(
            bool(datum.get("divergence_free_exact")),
            "the datum must be exactly divergence free",
        )

    externals = payload.get("external_theorems")
    require(isinstance(externals, dict), "missing external theorems")
    if isinstance(externals, dict):
        require(
            set(externals) == set(EXTERNAL_THEOREMS),
            "the payload must carry exactly the three named external theorems",
        )
        for name, entry in externals.items():
            require(
                entry.get("proved") is False,
                f"{name} claims to be proved; it is classical but not proved here",
            )
            require(
                entry.get("axiomatised_in_lean") is False,
                f"{name} must never be axiomatised in Lean",
            )
            require(
                entry.get("statement") == EXTERNAL_THEOREMS[name],
                f"{name} statement does not match the audited text",
            )

    conclusion = payload.get("conclusion")
    require(isinstance(conclusion, dict), "missing conclusion")

    constants = payload.get("constants")
    if isinstance(constants, dict) and isinstance(conclusion, dict):
        k1 = Fraction(str(constants["K1"]))
        k2 = Fraction(str(constants["K2"]))
        bounds = [Fraction(str(m)) for m in constants["sup_derivative_bounds_M1_to_M5"]]
        m1, m2, m3, m4, m5 = bounds
        require(
            k1 == 4 * m1 + 6 * m2 + 4 * m3 + m4,
            "K1 does not match its own derivative bounds",
        )
        require(
            k2 == m5 + 4 * m4 + 6 * m3 + 4 * m2 + m1,
            "K2 does not match its own derivative bounds",
        )
        nu = Fraction(str(payload["viscosity"]))
        require(
            Fraction(str(constants["linear_coefficient"]))
            == -nu + 9 * (k1 + k2),
            "the linear coefficient does not equal -nu + 9(K1+K2)",
        )
        adot = Fraction(str(constants["adot_upper"]))
        adot_sq = Fraction(str(constants["adot_squared_upper"]))
        require(adot * adot >= adot_sq, "adot_upper is not an upper square root")
        require(
            Fraction(str(constants["quadratic_coefficient"])) == 135 * adot,
            "the quadratic coefficient does not equal 135 adot",
        )

        control = payload.get("control_certificate")
        require(isinstance(control, dict), "missing control certificate")
        if isinstance(control, dict):
            inner = verify_control_certificate(control)
            require(bool(inner.get("verified")), "the control certificate fails")
            if conclusion.get("proved"):
                require(
                    bool(inner.get("proves_existence")),
                    "the conclusion claims a proof but the control slab did not close",
                )
            coefficients = control.get("coefficients", {})
            require(
                Interval.from_pair(coefficients["linear"]).upper
                == Fraction(str(constants["linear_coefficient"])),
                "the control ODE was fed a different linear coefficient",
            )
            require(
                Interval.from_pair(coefficients["quadratic"]).upper
                == Fraction(str(constants["quadratic_coefficient"])),
                "the control ODE was fed a different quadratic coefficient",
            )
            require(
                Interval.from_pair(coefficients["forcing"]).upper
                == Fraction(str(constants["residual_h4_upper"])),
                "the control ODE was fed a different residual bound",
            )

        if conclusion.get("proved"):
            sqrt8 = Fraction(str(conclusion["sqrt8_upper"]))
            require(sqrt8 * sqrt8 >= 8, "the sqrt8 constant is not an upper root")
            require(
                Fraction(str(conclusion["distance_h3"]))
                == sqrt8 * Fraction(str(conclusion["distance_h4_dot"])),
                "the H3 distance does not equal sqrt8 times the H4 distance",
            )

    claims = payload.get("claims")
    require(
        isinstance(claims, list)
        and any("Clay" in str(claim) for claim in claims)
        and any("not a singularity" in str(claim) for claim in claims),
        "the payload must record that it is not a singularity result and says "
        "nothing about the Clay problem",
    )
    return {"verified": not failures, "failures": failures}
