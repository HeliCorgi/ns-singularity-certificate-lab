r"""A rigorous scalar control ODE over one slab: a Picard box, then a Taylor tube.

A control ODE is the object that converts a *computed* trajectory into a
statement about an *interval of time*.  Every attempt in this repository to say
something about a time interval rather than about isolated instants has run
into the same wall: :mod:`ns_certificate_lab.snapshot_certificate` encloses one
instant, :mod:`ns_certificate_lab.slab_certificate` encloses one slab of the
computed object but only under two named and unproved hypotheses, and the
second of those hypotheses (the Hermite remainder ``H2``) is unfixable as it is
framed, because the two-point Hermite remainder ``Delta^4 M_4 / 384`` presumes
that the interpolant matches one and the same trajectory at both ends, which an
RK4 output does not.  This module is the replacement: instead of estimating how
far the interpolant strays from a trajectory, it *proves* that a trajectory
exists on the whole slab and stays inside an explicitly named box.

The honest form of the coefficients
-----------------------------------
The control ODE is written

.. math::   R'(t) = a(t)\,R(t) + b(t)\,R(t)^2 + e(t),\qquad R(t_n) \in [R_0] ,

with ``a``, ``b`` and ``e`` supplied as **interval constants over the slab**.
That is deliberate.  The coefficients of a real control ODE are not point
values: they are bounds on quantities like ``G_n D_n + K_n D_{n+1}`` assembled
from a numerical solution over the slab, and the only defensible thing to say
about such a quantity is that it lies between two numbers.  Handing the solver
point values and then claiming rigour would smuggle the entire question of
coefficient accuracy into a place where it cannot be seen.

Two stages, and why they are different in kind
----------------------------------------------
**Stage B, the rough enclosure, is a theorem, not an estimate.**  Given a
candidate box ``B`` containing ``[R_0]``, the single inclusion

.. math::   [R_0] + [0,h]\,F(B) \subseteq B ,\qquad
            F(B) = a B + b B^2 + e ,

proves that the initial value problem has a solution on the *whole* slab
``[t_n, t_n + h]`` and that the solution never leaves ``B``.  The proof is
Picard: the operator ``(T R)(t) = R(t_n) + \int_{t_n}^t f(s, R(s))\,ds`` maps
``C([t_n, t_n+h], B)`` into itself, because for any admissible selection the
integrand lies in ``F(B)`` and the integral over a subinterval of length
``t - t_n \in [0,h]`` therefore lies in ``[0,h]\,F(B)``.  The set ``B`` is
nonempty, convex, closed and bounded, so ``C([t_n,t_n+h], B)`` is a nonempty
closed subset of a Banach space; ``f`` is a polynomial in ``R`` hence Lipschitz
on the bounded ``B``; and in the Bielecki norm
``\|R\| = \sup e^{-2L(t-t_n)}|R(t)|`` the operator ``T`` is a contraction with
factor ``1/2`` **for every** ``h``, so no smallness condition is needed and no
hypothesis is left over.  (Bielecki, *Bull. Acad. Polon. Sci.* 4 (1956)
261-264.)  This is why the inclusion is tested rather than approximated: if it
holds the conclusion follows outright, and if it fails the honest report is a
failure.

**Stage C, the tight enclosure, is an interval Taylor method.**  The rough box
is what makes it legitimate: the Lagrange form of the Taylor remainder needs
the ``(p+1)``-st Taylor coefficient evaluated somewhere on the trajectory, and
"somewhere on the trajectory" is exactly what Stage B has pinned down to ``B``.
Running the coefficient recurrence a second time with ``a_0 = B`` therefore
encloses that remainder.  This is the standard Moore-Lohner construction
(Moore, *Interval Analysis*, 1966, ch. 11; Lohner 1987; Nedialkov, Jackson and
Corliss, *Appl. Math. Comput.* 105 (1999) 21-68), reduced to one scalar
equation so that every step can be read.

Why Bernstein and not a derivative bound
----------------------------------------
The Taylor polynomial still has to be *ranged* over ``[0,h]``, and evaluating a
polynomial in interval arithmetic by Horner or by term-by-term substitution
gives an overestimate that grows with the degree.  The Bernstein form avoids
the question entirely.  Converting the monomial coefficients ``c_k`` to

.. math::   b_j = \sum_{k=0}^{j} \frac{\binom{j}{k}}{\binom{p}{k}}\,c_k\,h^k

expresses the polynomial as ``\sum_j b_j B_j^p(t)`` with ``t = s/h`` and
``B_j^p`` the Bernstein basis on ``[0,1]``.  Those basis functions are
nonnegative and sum to one, so the polynomial's value at every point of the
interval is a **convex combination** of the ``b_j``, and

.. math::   \min_j b_j \le q(t) \le \max_j b_j \qquad \text{on } [0,1] .

That is the convex hull property of the Bernstein-Bezier form (Cargo and
Shisha, *J. Res. Nat. Bur. Standards* 70B (1966) 79-81; Garloff, *Freiburger
Intervall-Berichte* 85/9; Farin, *Curves and Surfaces for CAGD*, ch. 4).  It
needs no derivative bound, no subdivision and no hypothesis, which is precisely
the property the discarded Hermite remainder did not have.

What this module does **not** establish
---------------------------------------
It certifies a *scalar control ODE*.  The conclusion is conditional: **if** the
supplied coefficient intervals bound the true Sobolev-norm dynamics, **then**
the exact norm stays inside the stated tube for the stated time.  Establishing
the antecedent is HS-5, and HS-5 is not closed.  Nothing here is a statement
about the Navier-Stokes solution and nothing here bears on the Clay problem.

Separately, and worth stating because it is the most likely misquotation: the
``-\nu R_n`` damping term of the periodic-domain control ODE is **absent** on
``\mathbb{R}^3``.  On the torus, zero-mean fields have Laplacian spectrum
starting at ``|k| = 1``, which supplies that term; on ``\mathbb{R}^3`` the
spectrum is ``[0,\infty)`` and the term vanishes.  The existence interval
certified here is therefore strictly shorter than a periodic-domain analogue
would give, and a torus result quoted for the whole space would be wrong.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Sequence

from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

try:  # the parallel L3 module owns these; duplicate only if it has not landed
    from .l3_certificate import exp_interval
except ImportError:  # pragma: no cover - exercised only before the merge

    def exp_interval(value: Interval, *, terms: int = 40) -> Interval:
        """Fallback enclosure of ``exp``; merge with the L3 module's copy.

        ``exp`` is increasing, so only a lower bound at the left endpoint and an
        upper bound at the right endpoint are needed.  Arguments are halved
        until they are at most ``1/2``, where the truncated series has a
        geometric tail, and the result is squared back up.
        """
        if terms < 4:
            raise ValueError("at least four Taylor terms are required")

        def bounds_nonneg(y: Fraction) -> tuple[Fraction, Fraction]:
            squarings = 0
            while y > Fraction(1, 2):
                y /= 2
                squarings += 1
            partial = Fraction(0)
            term = Fraction(1)
            for k in range(terms):
                if k:
                    term = term * y / k
                partial += term
            ratio = y / (terms + 1)
            remainder = term * y / terms * ratio / (1 - ratio)
            low, high = partial, partial + remainder
            for _ in range(squarings):
                low, high = low * low, high * high
            return low, high

        low = (
            bounds_nonneg(value.lower)[0]
            if value.lower >= 0
            else 1 / bounds_nonneg(-value.lower)[1]
        )
        high = (
            bounds_nonneg(value.upper)[1]
            if value.upper >= 0
            else 1 / bounds_nonneg(-value.upper)[0]
        )
        return Interval(low, high)


__all__ = [
    "BERNSTEIN_HULL_PROPERTY",
    "CONTROL_CLAIMS",
    "CONTROL_LIMITATIONS",
    "ControlCertificate",
    "ControlCoefficients",
    "PICARD_SELF_MAPPING",
    "RoughEnclosure",
    "TaylorTube",
    "bernstein_coefficients",
    "build_control_certificate",
    "prove_rough_enclosure",
    "taylor_coefficients",
    "taylor_tube",
    "verify_control_certificate",
]


#: The statement the self-mapping test proves, recorded verbatim in every
#: payload so that a checker can refuse a payload whose theorem has been edited.
PICARD_SELF_MAPPING = (
    "If B is a compact interval containing the initial interval [R_0] and "
    "[R_0] + [0,h] F(B) is contained in B, where F(B) = a B + b B^2 + e is "
    "evaluated in interval arithmetic, then for every selection of coefficients "
    "inside a, b, e and every initial value in [R_0] the initial value problem "
    "has a solution on the whole of [t_n, t_n + h] and that solution stays in "
    "B.  The Picard operator maps C([t_n,t_n+h], B) into itself; B is nonempty, "
    "convex, closed and bounded; the right-hand side is a polynomial and hence "
    "Lipschitz on B; and in the Bielecki norm the operator is a contraction for "
    "every h, so the conclusion carries no smallness hypothesis."
)

#: The property that makes the polynomial range enclosure exact-by-construction.
BERNSTEIN_HULL_PROPERTY = (
    "The Bernstein basis of degree p on [0,1] is nonnegative and sums to one, "
    "so a polynomial written in that basis takes, at every point, a convex "
    "combination of its Bernstein coefficients.  Its range on [0,h] is "
    "therefore contained in the hull of the coefficients b_j, with no "
    "derivative bound, no subdivision and no hypothesis.  Cargo and Shisha, "
    "J. Res. Nat. Bur. Standards 70B (1966) 79-81."
)

#: Named limitations.  Each is recorded in the payload with ``closed: False``;
#: the checker compares the text verbatim and refuses a payload that drops one.
CONTROL_LIMITATIONS = {
    "HS5_coefficients_to_pde": (
        "The coefficient intervals a, b, e are an input, not a conclusion.  "
        "Deriving them from a computed trajectory so that they provably "
        "dominate the true Sobolev-norm dynamics is HS-5, and HS-5 is not "
        "closed.  What is certified here is the implication: if the control "
        "ODE dominates the true dynamics then the true norm stays inside the "
        "stated tube for the stated time."
    ),
    "C1_frozen_coefficients": (
        "The rough box holds for every measurable coefficient selection inside "
        "the stated intervals, because its proof uses nothing about the "
        "right-hand side beyond the bound F(B).  The tight Taylor tube is "
        "proved for selections that are constant on the slab: the Taylor "
        "recurrence differentiates the vector field, and a time-varying "
        "coefficient would contribute derivatives of that coefficient which "
        "the interval data does not carry."
    ),
    "R3_has_no_spectral_gap": (
        "The -nu R damping term of the periodic-domain control ODE is absent "
        "here.  On the torus the Laplacian spectrum on zero-mean fields starts "
        "at |k| = 1 and supplies that term; on R^3 the spectrum is [0,infinity) "
        "and it vanishes.  The existence interval certified here is therefore "
        "strictly shorter than a periodic-domain analogue would give, and "
        "quoting a torus result for the whole space would be wrong."
    ),
}

#: Claims recorded in every payload.  The Clay disclaimer is mandatory and the
#: checker refuses a payload without it.
CONTROL_CLAIMS = [
    "Rigorous and unconditional, given the coefficient intervals as inputs: "
    "the scalar initial value problem has a solution on the whole certified "
    "slab, that solution is unique for each admissible coefficient selection, "
    "and it stays inside the recorded rough box.",
    "Rigorous, given the same inputs and a coefficient selection constant on "
    "the slab: the solution stays inside the recorded Taylor tube, which is "
    "the Bernstein hull of the truncated Taylor polynomial plus an interval "
    "remainder evaluated on the rough box.",
    "Conditional on HS-5: if the coefficient intervals bound the true "
    "Sobolev-norm dynamics then the true norm obeys the same enclosures.  That "
    "step is not proved anywhere in this repository.",
    "This certifies a scalar control ODE.  It is not a regularity result, it "
    "is not a singularity result, and nothing here bears on the Clay problem.",
]


# --------------------------------------------------------------------------- #
# interval helpers the shared Interval class deliberately does not carry       #
# --------------------------------------------------------------------------- #


def _square(value: Interval) -> Interval:
    """The **exact** range of ``x -> x^2`` over an interval.

    ``value * value`` treats the two factors as independent and returns
    ``[-1, 1]`` for ``[-1, 1]``, which is not the range of any square.  The
    dependency matters here because the quadratic term is the entire reason a
    control ODE has a finite existence time, so the square is computed as a
    range rather than as a product.
    """
    lower, upper = value.lower, value.upper
    if lower >= 0:
        return Interval(lower * lower, upper * upper)
    if upper <= 0:
        return Interval(upper * upper, lower * lower)
    return Interval(Fraction(0), max(lower * lower, upper * upper))


def _subset(inner: Interval, outer: Interval) -> bool:
    return outer.lower <= inner.lower and inner.upper <= outer.upper


def _intersection(first: Interval, second: Interval) -> Interval:
    lower = max(first.lower, second.lower)
    upper = min(first.upper, second.upper)
    if lower > upper:
        raise ValueError("the two enclosures do not overlap")
    return Interval(lower, upper)


def _width(value: Interval) -> Fraction:
    return value.upper - value.lower


# --------------------------------------------------------------------------- #
# the coefficients                                                             #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ControlCoefficients:
    """The three interval constants of ``R' = a R + b R^2 + e`` on one slab.

    They are intervals rather than numbers because that is what a control ODE
    assembled from a numerical solution actually gives: bounds valid over the
    slab, not values at an instant.
    """

    linear: Interval
    quadratic: Interval
    forcing: Interval

    @staticmethod
    def from_floats(
        *, linear: float | Fraction, quadratic: float | Fraction,
        forcing: float | Fraction,
    ) -> "ControlCoefficients":
        """Degenerate intervals from point values, exactly (floats are dyadic)."""
        return ControlCoefficients(
            Interval.exact(linear), Interval.exact(quadratic), Interval.exact(forcing)
        )

    def field(self, value: Interval) -> Interval:
        """``F([x]) = a[x] + b[x]^2 + e``, the natural interval extension."""
        return self.linear * value + self.quadratic * _square(value) + self.forcing

    def lipschitz_bound(self, box: Interval) -> Fraction:
        """``|a| + 2|b| max|B|``, a Lipschitz constant for ``F`` on ``B``.

        This is what makes the Picard operator a contraction in the Bielecki
        norm and hence what makes the solution unique; it is recorded so that a
        reader can see the uniqueness half of the theorem has a number behind
        it, not only an assertion.
        """
        return self.linear.magnitude + 2 * self.quadratic.magnitude * box.magnitude

    def as_dict(self) -> dict[str, list[str]]:
        return {
            "linear": self.linear.as_pair(),
            "quadratic": self.quadratic.as_pair(),
            "forcing": self.forcing.as_pair(),
        }

    @staticmethod
    def from_dict(block: dict[str, Sequence[str]]) -> "ControlCoefficients":
        return ControlCoefficients(
            Interval.from_pair(block["linear"]),
            Interval.from_pair(block["quadratic"]),
            Interval.from_pair(block["forcing"]),
        )


# --------------------------------------------------------------------------- #
# stage B: the rough enclosure                                                 #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class RoughEnclosure:
    """The outcome of the self-mapping test, success or failure alike.

    A failure is a first-class result and carries the last box that was tried,
    the two exhausted budgets and a verdict in words.  Returning a widened box
    with a quiet flag would be the one outcome this module must never produce.
    """

    proved: bool
    initial: Interval
    box: Interval
    field_on_box: Interval
    image: Interval
    step: Fraction
    requested_step: Fraction
    doublings: int
    halvings: int
    inflation: Fraction
    attempts: int
    verdict: str


def prove_rough_enclosure(
    coefficients: ControlCoefficients,
    initial: Interval,
    step: float | Fraction,
    *,
    max_doublings: int = 40,
    max_halvings: int = 6,
    inflation_floor: Fraction = Fraction(1, 1024),
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> RoughEnclosure:
    r"""Search for a box ``B`` satisfying ``[R_0] + [0,h] F(B) \subseteq B``.

    The *search* is a heuristic; the *test* is the theorem.  That separation is
    what makes a failed search harmless: a box is only ever returned as proved
    when the inclusion has actually been checked in exact arithmetic, and the
    search failing means the function says so.

    The heuristic starts from ``rho = 2 h |F([R_0])|`` inflated by a floor
    proportional to ``1 + |[R_0]|``.  The leading term is the distance the
    solution can travel in the slab, doubled; the floor exists because
    ``F([R_0])`` can vanish at a momentarily stationary datum while the
    dynamics nearby does not, and because a search that starts from zero can
    never double its way anywhere.  On failure ``rho`` doubles, up to
    ``max_doublings`` times, and then ``h`` halves, up to ``max_halvings``
    times.  Doubling is not monotone for a quadratic field -- past a point a
    wider box has a strictly worse image, because ``b B^2`` grows quadratically
    while ``B`` grows linearly -- so exhausting the doublings is genuinely
    informative: it usually means the slab is too long, which is why the outer
    loop halves ``h`` rather than inflating further.

    Every box tested is rounded outward first and the field is then evaluated
    on the rounded box, so the inclusion that is checked is the inclusion for
    the box that is reported.  Rounding the image outward can only make the
    test harder, never easier.
    """
    h = Fraction(step)
    if h <= 0:
        raise ValueError("the slab length must be positive")
    if max_doublings < 0 or max_halvings < 0:
        raise ValueError("the retry budgets must be nonnegative")
    if inflation_floor <= 0:
        raise ValueError("the inflation floor must be positive")

    requested = h
    drift = coefficients.field(initial).magnitude
    floor = inflation_floor * (1 + initial.magnitude)
    attempts = 0
    box = field = image = initial
    rho = floor

    for halvings in range(max_halvings + 1):
        rho = 2 * h * drift + floor
        for doublings in range(max_doublings + 1):
            attempts += 1
            box = Interval(
                initial.lower - rho, initial.upper + rho
            ).round_outward(precision_bits)
            field = coefficients.field(box).round_outward(precision_bits)
            image = (
                initial + Interval(Fraction(0), h) * field
            ).round_outward(precision_bits)
            if _subset(image, box) and _subset(initial, box):
                return RoughEnclosure(
                    proved=True,
                    initial=initial,
                    box=box,
                    field_on_box=field,
                    image=image,
                    step=h,
                    requested_step=requested,
                    doublings=doublings,
                    halvings=halvings,
                    inflation=rho,
                    attempts=attempts,
                    verdict=(
                        f"self-mapping proved after {doublings} doubling(s) and "
                        f"{halvings} halving(s) of the slab length"
                    ),
                )
            rho *= 2
        if halvings < max_halvings:
            h = h / 2

    return RoughEnclosure(
        proved=False,
        initial=initial,
        box=box,
        field_on_box=field,
        image=image,
        step=h,
        requested_step=requested,
        doublings=max_doublings,
        halvings=max_halvings,
        inflation=rho,
        attempts=attempts,
        verdict=(
            f"no self-mapping box found: {max_doublings} doubling(s) exhausted at "
            f"each of {max_halvings + 1} slab lengths down to {h}.  Reported as a "
            "failure rather than inflated further, because a box that fails the "
            "inclusion proves nothing however wide it is"
        ),
    )


# --------------------------------------------------------------------------- #
# stage C: the Taylor coefficients, the Bernstein form and the tube            #
# --------------------------------------------------------------------------- #


def taylor_coefficients(
    coefficients: ControlCoefficients,
    seed: Interval,
    order: int,
    *,
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> list[Interval]:
    r"""The Taylor coefficients ``a_0 \dots a_{order}`` of the solution.

    Differentiating ``R' = a R + b R^2 + e`` and matching powers of ``s`` gives

    .. math::  (k+1)\,a_{k+1} = a\,a_k + b \sum_{i=0}^{k} a_i a_{k-i}
               + e\,[k = 0] ,

    which is the whole recurrence.  The convolution is assembled with the
    diagonal term ``a_{k/2}^2`` taken as an exact square and the off-diagonal
    pairs counted once and doubled; that is not an optimisation but a
    correctness point, since an interval straddling zero has a square strictly
    narrower than its self-product.

    Called with a thin ``seed`` this produces the polynomial part.  Called with
    the proved rough box as ``seed`` and ``order = p + 1`` it produces an
    enclosure of the Lagrange remainder coefficient, because for an autonomous
    right-hand side each ``a_k`` is a fixed polynomial function of the state and
    the state at the unknown intermediate point is known to lie in the box.
    """
    if order < 0:
        raise ValueError("the Taylor order must be nonnegative")
    zero = Interval(Fraction(0), Fraction(0))
    out = [seed]
    for k in range(order):
        convolution = zero
        for i in range(k + 1):
            if 2 * i == k:
                convolution = convolution + _square(out[i])
            elif 2 * i < k:
                convolution = convolution + (out[i] * out[k - i]).scale(Fraction(2))
        total = coefficients.linear * out[k] + coefficients.quadratic * convolution
        if k == 0:
            total = total + coefficients.forcing
        out.append(total.scale(Fraction(1, k + 1)).round_outward(precision_bits))
    return out


def bernstein_coefficients(
    monomial: Sequence[Interval],
    step: Fraction,
    *,
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> list[Interval]:
    r"""Convert monomial coefficients on ``[0,h]`` to Bernstein coefficients.

    With ``s = h t`` the polynomial ``\sum_k c_k s^k`` becomes
    ``\sum_k (c_k h^k) t^k`` on ``[0,1]``, and the monomial-to-Bernstein
    conversion ``t^k = \sum_{j \ge k} \binom{j}{k}/\binom{p}{k}\,B_j^p(t)``
    gives

    .. math::  b_j = \sum_{k=0}^{j} \frac{\binom{j}{k}}{\binom{p}{k}}\,c_k h^k .

    All the weights are nonnegative rationals, so the conversion is a single
    exact rational computation; the only division is by an integer binomial
    coefficient, done inside a :class:`~fractions.Fraction`.

    The point of doing this at all is stated in :data:`BERNSTEIN_HULL_PROPERTY`:
    the hull of the ``b_j`` encloses the polynomial's range on the whole
    interval, with no derivative bound anywhere.
    """
    degree = len(monomial) - 1
    if degree < 0:
        raise ValueError("a polynomial needs at least one coefficient")
    out: list[Interval] = []
    for j in range(degree + 1):
        accumulated = Interval(Fraction(0), Fraction(0))
        for k in range(j + 1):
            weight = Fraction(comb(j, k), comb(degree, k)) * step**k
            accumulated = accumulated + monomial[k].scale(weight)
        out.append(accumulated.round_outward(precision_bits))
    return out


@dataclass(frozen=True)
class TaylorTube:
    """The tight enclosure of the solution over the whole slab."""

    order: int
    step: Fraction
    coefficients: tuple[Interval, ...]
    remainder_coefficient: Interval
    bernstein: tuple[Interval, ...]
    polynomial_enclosure: Interval
    remainder: Interval
    tube: Interval
    refined_tube: Interval


def taylor_tube(
    coefficients: ControlCoefficients,
    rough: RoughEnclosure,
    *,
    order: int = 6,
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> TaylorTube:
    r"""Enclose the solution on ``[t_n, t_n + h]`` by an interval Taylor form.

    Requires a *proved* rough enclosure, and refuses to run without one.  That
    is not defensive programming: the remainder term is only an enclosure
    because the intermediate state is known to lie in the box, and without the
    box the whole construction would be an extrapolation dressed up as a bound.

    The result is

    .. math::  R(t_n + s) \in \sum_{k \le p} a_k s^k
               + [a_{p+1}]\,[0,h]^{p+1},\qquad s \in [0,h] ,

    with the polynomial part ranged by the Bernstein hull.  Because
    ``[0,h]^{p+1} = [0, h^{p+1}]`` contains zero, the remainder interval always
    contains zero, so the tube always contains the polynomial hull and in
    particular contains the initial interval at ``s = 0``.

    The tube is finally intersected with the rough box.  Both are valid
    enclosures of the same solution, so their intersection is one too, and on a
    long slab the box is sometimes the tighter of the pair.
    """
    if not rough.proved:
        raise ValueError(
            "a Taylor tube needs a proved rough enclosure; without one the "
            "remainder coefficient has nothing to be evaluated on"
        )
    if order < 1:
        raise ValueError("the Taylor order must be at least one")

    step = rough.step
    monomial = taylor_coefficients(
        coefficients, rough.initial, order, precision_bits=precision_bits
    )
    on_box = taylor_coefficients(
        coefficients, rough.box, order + 1, precision_bits=precision_bits
    )
    remainder_coefficient = on_box[order + 1]
    bernstein = bernstein_coefficients(
        monomial, step, precision_bits=precision_bits
    )
    polynomial = Interval.hull(bernstein)
    remainder = (
        remainder_coefficient * Interval(Fraction(0), step ** (order + 1))
    ).round_outward(precision_bits)
    tube = (polynomial + remainder).round_outward(precision_bits)
    return TaylorTube(
        order=order,
        step=step,
        coefficients=tuple(monomial),
        remainder_coefficient=remainder_coefficient,
        bernstein=tuple(bernstein),
        polynomial_enclosure=polynomial,
        remainder=remainder,
        tube=tube,
        refined_tube=_intersection(tube, rough.box),
    )


# --------------------------------------------------------------------------- #
# the certificate                                                              #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ControlCertificate:
    """A machine-checkable record of one certified control-ODE slab."""

    payload: dict[str, object]

    def as_dict(self) -> dict[str, object]:
        return self.payload


def build_control_certificate(
    coefficients: ControlCoefficients,
    initial: Interval,
    step: float | Fraction,
    *,
    time_start: float = 0.0,
    order: int = 6,
    max_doublings: int = 40,
    max_halvings: int = 6,
    inflation_floor: Fraction = Fraction(1, 1024),
    precision_bits: int = DEFAULT_PRECISION_BITS,
    label: str = "scalar_control_ode_slab",
) -> ControlCertificate:
    """Prove the self-mapping inclusion, then build the tube, then serialise.

    A failed self-mapping still produces a certificate.  The payload records the
    failure, carries a null tube, and verifies as a *consistent* record of a
    *failed* proof: :func:`verify_control_certificate` reports ``verified`` for
    internal consistency and ``proves_existence`` separately.  Refusing to emit
    anything on failure would make an unattempted slab indistinguishable from
    one that was attempted and did not close.
    """
    rough = prove_rough_enclosure(
        coefficients,
        initial,
        step,
        max_doublings=max_doublings,
        max_halvings=max_halvings,
        inflation_floor=inflation_floor,
        precision_bits=precision_bits,
    )

    tube_block: dict[str, object] | None = None
    picard_block: dict[str, object] | None = None
    if rough.proved:
        tube = taylor_tube(
            coefficients, rough, order=order, precision_bits=precision_bits
        )
        lipschitz = coefficients.lipschitz_bound(rough.box)
        growth = exp_interval(
            Interval.exact(lipschitz * rough.step)
        ).round_outward(precision_bits)
        picard_block = {
            "statement": PICARD_SELF_MAPPING,
            "lipschitz_bound_on_box": str(lipschitz),
            "lipschitz_times_step": str(lipschitz * rough.step),
            "gronwall_growth_factor": growth.as_pair(),
        }
        tube_block = {
            "order": tube.order,
            "bernstein_property": BERNSTEIN_HULL_PROPERTY,
            "taylor_coefficients": [item.as_pair() for item in tube.coefficients],
            "remainder_coefficient": tube.remainder_coefficient.as_pair(),
            "bernstein_coefficients": [item.as_pair() for item in tube.bernstein],
            "polynomial_enclosure": tube.polynomial_enclosure.as_pair(),
            "remainder": tube.remainder.as_pair(),
            "tube": tube.tube.as_pair(),
            "refined_tube": tube.refined_tube.as_pair(),
            "tube_width": str(_width(tube.tube)),
            "refined_tube_width": str(_width(tube.refined_tube)),
            "rough_box_width": str(_width(rough.box)),
            "tightness_ratio": float(
                Fraction(_width(tube.refined_tube), _width(rough.box))
            )
            if _width(rough.box) > 0
            else None,
        }

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate": "scalar_control_ode_slab",
        "label": label,
        "precision_bits": precision_bits,
        "coefficients": coefficients.as_dict(),
        "slab": {
            "time_start": float(time_start),
            "time_end": float(time_start) + float(rough.step),
            "requested_step": str(rough.requested_step),
            "step": str(rough.step),
            "halvings": rough.halvings,
        },
        "initial": initial.as_pair(),
        "self_mapping": {
            "statement": PICARD_SELF_MAPPING,
            "proved": rough.proved,
            "box": rough.box.as_pair(),
            "field_on_box": rough.field_on_box.as_pair(),
            "image": rough.image.as_pair(),
            "doublings": rough.doublings,
            "inflation": str(rough.inflation),
            "attempts": rough.attempts,
            "verdict": rough.verdict,
        },
        "picard": picard_block,
        "tube": tube_block,
        "limitations": {
            name: {"statement": text, "closed": False}
            for name, text in CONTROL_LIMITATIONS.items()
        },
        "claims": list(CONTROL_CLAIMS),
    }
    return ControlCertificate(payload)


def verify_control_certificate(payload: dict[str, object]) -> dict[str, object]:
    """Independently re-check a control payload from its own numbers alone.

    The checker never calls the solver.  It reads the coefficients, the box and
    the recorded coefficients of the Taylor and Bernstein forms, and re-derives
    every relation the payload asserts between them.  A payload whose box has
    been narrowed, whose ``F(B)`` no longer satisfies the inclusion, whose
    Bernstein enclosure is tighter than the hull of its own coefficients, or
    whose disclaimer has been removed fails here even though it is still
    syntactically valid JSON.

    Adversarial input gets a verdict, never an exception: a crash is
    indistinguishable from a bug in the checker, so every malformed field is
    reported as a failure instead.
    """
    try:
        return _verify(payload)
    except Exception as error:  # noqa: BLE001 -- adversarial input is a failure
        return {
            "verified": False,
            "proves_existence": False,
            "failures": [f"payload rejected while parsing: {error}"],
        }


def _verify(payload: dict[str, object]) -> dict[str, object]:
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(
        payload.get("certificate") == "scalar_control_ode_slab",
        "wrong certificate kind",
    )
    require(payload.get("schema_version") == 1, "unsupported schema version")

    coefficients = ControlCoefficients.from_dict(
        payload["coefficients"]  # type: ignore[arg-type]
    )
    for name, item in (
        ("a", coefficients.linear),
        ("b", coefficients.quadratic),
        ("e", coefficients.forcing),
    ):
        require(item.lower <= item.upper, f"coefficient {name} is inverted")

    slab = payload["slab"]
    step = Fraction(str(slab["step"]))  # type: ignore[index]
    requested = Fraction(str(slab["requested_step"]))  # type: ignore[index]
    halvings = int(slab["halvings"])  # type: ignore[index]
    require(step > 0, "the certified slab length must be positive")
    require(
        step * 2**halvings == requested,
        "the recorded slab length is not the requested one halved the recorded "
        "number of times",
    )
    require(
        abs(
            (float(slab["time_end"]) - float(slab["time_start"]))  # type: ignore[index]
            - float(step)
        )
        <= 1e-12 * max(1.0, abs(float(slab["time_end"]))),  # type: ignore[index]
        "the slab endpoints disagree with the certified step",
    )

    initial = Interval.from_pair(payload["initial"])  # type: ignore[arg-type]
    require(initial.lower <= initial.upper, "the initial enclosure is inverted")

    mapping = payload["self_mapping"]
    require(
        mapping.get("statement") == PICARD_SELF_MAPPING,  # type: ignore[union-attr]
        "the self-mapping statement does not match the audited text",
    )
    box = Interval.from_pair(mapping["box"])  # type: ignore[index]
    field = Interval.from_pair(mapping["field_on_box"])  # type: ignore[index]
    image = Interval.from_pair(mapping["image"])  # type: ignore[index]
    for name, item in (("box", box), ("field_on_box", field), ("image", image)):
        require(item.lower <= item.upper, f"the {name} enclosure is inverted")

    # -- the three relations the self-mapping claim actually rests on -------- #
    recomputed_field = coefficients.field(box)
    require(
        _subset(recomputed_field, field),
        "the recorded F(B) does not enclose a B + b B^2 + e on the recorded box",
    )
    recomputed_image = initial + Interval(Fraction(0), step) * field
    require(
        _subset(recomputed_image, image),
        "the recorded image does not enclose [R_0] + [0,h] F(B)",
    )
    claims_mapping = bool(mapping["proved"])  # type: ignore[index]
    holds = _subset(recomputed_image, box) and _subset(initial, box)
    require(
        claims_mapping == holds,
        "the self-mapping flag contradicts the inclusion its own numbers give",
    )
    if claims_mapping:
        require(holds, "the recorded F(B) does not satisfy the self-mapping inclusion")

    tube_block = payload.get("tube")
    picard = payload.get("picard")
    require(
        (tube_block is None) == (not claims_mapping),
        "a tube is present exactly when the self-mapping is claimed, and here "
        "the two disagree",
    )

    if claims_mapping and isinstance(picard, dict):
        lipschitz = Fraction(str(picard["lipschitz_bound_on_box"]))
        require(
            lipschitz >= coefficients.lipschitz_bound(box),
            "the recorded Lipschitz constant is below |a| + 2|b| max|B|",
        )
        require(
            Fraction(str(picard["lipschitz_times_step"])) == lipschitz * step,
            "the recorded L h does not equal the Lipschitz constant times the step",
        )
        growth = Interval.from_pair(picard["gronwall_growth_factor"])
        require(
            growth.upper >= 1 + lipschitz * step,
            "the recorded growth factor is below exp(L h) >= 1 + L h",
        )

    if claims_mapping and isinstance(tube_block, dict):
        order = int(tube_block["order"])
        require(order >= 1, "the Taylor order must be at least one")
        require(
            tube_block.get("bernstein_property") == BERNSTEIN_HULL_PROPERTY,
            "the Bernstein statement does not match the audited text",
        )
        monomial = [
            Interval.from_pair(pair)
            for pair in tube_block["taylor_coefficients"]  # type: ignore[union-attr]
        ]
        require(
            len(monomial) == order + 1,
            "the number of Taylor coefficients disagrees with the stated order",
        )
        require(
            _subset(monomial[0], initial) and _subset(initial, monomial[0]),
            "the zeroth Taylor coefficient is not the initial enclosure",
        )
        # Re-run the recurrence: each recorded coefficient must enclose the one
        # its own predecessors imply.
        implied = taylor_coefficients(coefficients, monomial[0], order)
        for index in range(1, order + 1):
            require(
                _subset(implied[index], monomial[index]),
                f"Taylor coefficient {index} does not enclose the value the "
                "recurrence gives from the recorded coefficients",
            )
        remainder_coefficient = Interval.from_pair(
            tube_block["remainder_coefficient"]  # type: ignore[arg-type]
        )
        implied_on_box = taylor_coefficients(coefficients, box, order + 1)
        require(
            _subset(implied_on_box[order + 1], remainder_coefficient),
            "the remainder coefficient does not enclose the recurrence evaluated "
            "on the rough box",
        )

        bernstein = [
            Interval.from_pair(pair)
            for pair in tube_block["bernstein_coefficients"]  # type: ignore[union-attr]
        ]
        require(
            len(bernstein) == order + 1,
            "the Bernstein form must have one coefficient per degree",
        )
        implied_bernstein = bernstein_coefficients(monomial, step)
        for index, item in enumerate(bernstein):
            require(item.lower <= item.upper, f"Bernstein coefficient {index} inverted")
            require(
                _subset(implied_bernstein[index], item),
                f"Bernstein coefficient {index} does not enclose the conversion "
                "of the recorded monomial coefficients",
            )

        polynomial = Interval.from_pair(
            tube_block["polynomial_enclosure"]  # type: ignore[arg-type]
        )
        hull = Interval.hull(bernstein)
        require(
            _subset(hull, polynomial),
            "the Bernstein enclosure is narrower than the hull of its own "
            "recorded coefficients, which the convex hull property forbids",
        )

        remainder = Interval.from_pair(tube_block["remainder"])  # type: ignore[arg-type]
        implied_remainder = remainder_coefficient * Interval(
            Fraction(0), step ** (order + 1)
        )
        require(
            _subset(implied_remainder, remainder),
            "the recorded remainder does not enclose [a_{p+1}] [0,h]^{p+1}",
        )
        require(
            remainder.contains_zero,
            "the remainder interval must contain zero, because [0,h]^{p+1} does",
        )

        tube = Interval.from_pair(tube_block["tube"])  # type: ignore[arg-type]
        refined = Interval.from_pair(tube_block["refined_tube"])  # type: ignore[arg-type]
        require(
            _subset(polynomial + remainder, tube),
            "the tube does not enclose the polynomial hull plus the remainder",
        )
        require(_subset(initial, tube), "the tube does not contain the initial data")
        require(
            _subset(refined, tube) and _subset(refined, box),
            "the refined tube is not contained in both the tube and the box",
        )
        require(
            _subset(initial, refined),
            "the refined tube does not contain the initial data, so the two "
            "enclosures it intersects are inconsistent",
        )
        require(
            Fraction(str(tube_block["tube_width"])) == tube.upper - tube.lower,
            "the recorded tube width is not the width of the recorded tube",
        )
        require(
            Fraction(str(tube_block["rough_box_width"])) == box.upper - box.lower,
            "the recorded box width is not the width of the recorded box",
        )

    limitations = payload.get("limitations")
    require(isinstance(limitations, dict), "missing limitations block")
    if isinstance(limitations, dict):
        require(
            set(limitations) == set(CONTROL_LIMITATIONS),
            "the payload must carry exactly the three named limitations",
        )
        for name, entry in limitations.items():
            require(
                entry.get("closed") is False,
                f"{name} claims to be closed; this module closes none of them",
            )
            require(
                entry.get("statement") == CONTROL_LIMITATIONS[name],
                f"the {name} statement does not match the audited text",
            )

    claims = payload.get("claims")
    require(
        isinstance(claims, list) and any("Clay" in str(item) for item in claims),
        "the payload must record that it says nothing about the Clay problem",
    )

    return {
        "verified": not failures,
        "proves_existence": bool(not failures and claims_mapping),
        "failures": failures,
    }
