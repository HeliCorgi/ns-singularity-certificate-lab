r"""Rigorous interval enclosures for the Gaussian--Hermite generator basis.

Why a second interval backend exists at all
-------------------------------------------
The flat-bump enclosures of :mod:`ns_certificate_lab.interval_generators` are
correct but catastrophically wide near the edge of a bump's support: the
second derivative of ``\chi(\sigma) = e^{-1/(1-\sigma)}`` carries a factor
``(1-\sigma)^{-4}``, and on any cell whose ``\sigma``-interval approaches ``1``
the only valid move is the crude sup bound, which inflates every downstream
quantity by tens of orders of magnitude.  Measured on the M-families the
divergence enclosure — an interval that must contain zero because the identity
is exact — came out at roughly ``\pm 1.45\times 10^2`` around zero, which makes
the certificate margins vacuous.

The Gaussian--Hermite basis removes the problem at its root.  A generator here
is ``P(s, z)\,e^{-\alpha s - \beta z^2}`` with ``s = r^2``, ``P`` a polynomial
with exact rational coefficients and ``\alpha, \beta`` positive rationals.
Differentiation **closes** on this family:

.. math::
   \partial_s\bigl(P e^{-\alpha s - \beta z^2}\bigr)
     = (\partial_sP - \alpha P)\,e^{-\alpha s - \beta z^2},\qquad
   \partial_z\bigl(P e^{-\alpha s - \beta z^2}\bigr)
     = (\partial_zP - 2\beta z P)\,e^{-\alpha s - \beta z^2},

so every partial the certificate needs is *another* polynomial — computed
exactly in :class:`fractions.Fraction` — times the *same* Gaussian.  The
Gaussian factor is enclosed sharply by monotonicity (its exponent is monotone
in ``s`` and in ``|z|``, so the range over a box is attained at box corners),
and each polynomial monomial ``s^mz^n`` is monotone in ``s`` on ``s\ge 0`` and
has known parity in ``z``, so it too is enclosed sharply.  There is no
``(1-\sigma)^{-4}`` anywhere, no interval division, and the only widening left
is the mild dependency loss of summing monomials and of multiplying the
polynomial box by the Gaussian box — both proportional to the cell width.

The price and how it is paid
----------------------------
A Gaussian is not compactly supported, so a certificate over a finite box must
bound the exterior contributions rigorously instead of getting them for free.
That is what :func:`gaussian_tail_moment` and :func:`exterior_viscous_tail`
provide: closed-form rational upper bounds on the exterior integrals, built
from the Mills-ratio bound

.. math::
   \int_R^\infty e^{-2\gamma\rho^2}\,d\rho \le \frac{e^{-2\gamma R^2}}{4\gamma R}
   \qquad(R > 0),

proved by inserting ``1 \le \rho/R`` into the integrand and integrating
``4\gamma\rho\,e^{-2\gamma\rho^2}`` exactly, together with repeated integration
by parts for higher moments and the calculus bound
``\sup_{t\ge0} t^me^{-at} = (m/(ae))^m`` for the polynomial prefactors.  The
tails decay like ``e^{-3\gamma R^2}`` and are folded into the certificate's
viscous and flux upper bounds, so being generous there costs a little margin
but never soundness.

Coordinate contract
-------------------
This module and the (parallel) ``gaussian_hermite`` symbolic module must agree
that ``s = r^2`` and that the stream generator produces

.. math::
   u^r = -r\,\partial_z\psi_1,\qquad u^z = 2\psi_1 + 2s\,\partial_s\psi_1,

which is the same audited recovery used by
:class:`~ns_certificate_lab.mixed_initial_data.MixedFamily` and is divergence
free as an algebraic identity.  :func:`as_gaussian_term` adapts anything that
carries ``(poly, alpha, beta)`` — a plain triple, or an object with those
attributes — so the wiring does not depend on the other module's class names.

Everything below is exact rational arithmetic with outward rounding
(:class:`~ns_certificate_lab.snapshot_certificate.Interval`); ``exp`` and
``sqrt`` are the audited series enclosures of
:mod:`ns_certificate_lab.l3_certificate`, never ``math.exp``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np
import numpy.typing as npt

from .interval_generators import CellEnclosure, GeneratorIntervals, divide, square
from .l3_certificate import exp_interval, sqrt_interval
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "ExteriorTailBound",
    "GaussianFamily",
    "GaussianTerm",
    "VelocityCellEnclosure",
    "as_gaussian_term",
    "build_gaussian_generation_certificate",
    "cell_enclosure_gaussian",
    "cell_velocity_enclosure",
    "clear_gaussian_cache",
    "exterior_viscous_tail",
    "gaussian_box",
    "gaussian_reference_field",
    "gaussian_tail_moment",
    "interval_poly_eval",
    "seed_gaussian_family",
]

_ZERO = Interval(Fraction(0), Fraction(0))
_ONE = Interval(Fraction(1), Fraction(1))

#: A strict upper bound on ``4\pi`` (``4\pi = 12.56637061435917295\ldots``).
_FOUR_PI_UPPER = Fraction(12566370614359173, 10**15)

#: The six partial-derivative slots, in the field order of
#: :class:`~ns_certificate_lab.interval_generators.GeneratorIntervals`.
_PARTIAL_KEYS = ("value", "ds", "dss", "dz", "dzz", "dsz")


def _round(value: Interval, bits: int) -> Interval:
    return value.round_outward(bits)


# --------------------------------------------------------------------------- #
# exact polynomial calculus in (s, z)                                          #
# --------------------------------------------------------------------------- #
#
# Polynomials are dictionaries {(m, n): Fraction} representing
# sum c_{mn} s^m z^n.  All manipulation below is exact; intervals enter only
# at evaluation time.  The derivative closure means the whole certificate
# needs nothing beyond addition, scaling, a monomial shift and the two formal
# derivatives.

_Poly = dict


def _iter_monomials(poly):
    """Yield ``(m, n, coefficient)`` from a dict or a canonical tuple."""
    if isinstance(poly, dict):
        for (m, n), coefficient in poly.items():
            yield int(m), int(n), Fraction(coefficient)
    else:
        for m, n, coefficient in poly:
            yield int(m), int(n), Fraction(coefficient)


def _poly_dict(poly) -> _Poly:
    out: _Poly = {}
    for m, n, coefficient in _iter_monomials(poly):
        if m < 0 or n < 0:
            raise ValueError("monomial exponents must be nonnegative")
        if coefficient != 0:
            out[(m, n)] = out.get((m, n), Fraction(0)) + coefficient
    return {key: value for key, value in out.items() if value != 0}


def _poly_add(first: _Poly, second: _Poly) -> _Poly:
    out = dict(first)
    for key, value in second.items():
        total = out.get(key, Fraction(0)) + value
        if total == 0:
            out.pop(key, None)
        else:
            out[key] = total
    return out


def _poly_scale(poly: _Poly, factor: Fraction) -> _Poly:
    if factor == 0:
        return {}
    return {key: value * factor for key, value in poly.items()}


def _poly_shift(poly: _Poly, dm: int, dn: int, factor: Fraction) -> _Poly:
    """Multiply by ``factor * s^dm * z^dn`` — exact, no widening."""
    if factor == 0:
        return {}
    return {(m + dm, n + dn): value * factor for (m, n), value in poly.items()}


def _poly_diff_s(poly: _Poly) -> _Poly:
    return {
        (m - 1, n): value * m for (m, n), value in poly.items() if m >= 1
    }


def _poly_diff_z(poly: _Poly) -> _Poly:
    return {
        (m, n - 1): value * n for (m, n), value in poly.items() if n >= 1
    }


# --------------------------------------------------------------------------- #
# the Gaussian--Hermite term and its adapters                                  #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class GaussianTerm:
    r"""One basis term ``P(s, z)\,e^{-\alpha s - \beta z^2}``, ``s = r^2``.

    The polynomial is stored as a **canonical sorted tuple** of
    ``(m, n, coefficient)`` rather than a dictionary so the term is hashable:
    the derivative closure and the one-dimensional Gaussian factors are cached
    per term across a grid sweep, and a hashable key is what makes the cache
    an ``O(n_r + n_z)`` affair instead of ``O(n_r n_z)``.

    ``\alpha`` and ``\beta`` must both be strictly positive: the tail bounds
    of :func:`exterior_viscous_tail` need genuine Gaussian decay in *both*
    coordinates, and a term without it would silently break every exterior
    integral.
    """

    monomials: tuple[tuple[int, int, Fraction], ...]
    alpha: Fraction
    beta: Fraction

    def __post_init__(self) -> None:
        if self.alpha <= 0 or self.beta <= 0:
            raise ValueError(
                "Gaussian decay rates alpha and beta must be strictly positive; "
                "the exterior tail bounds are meaningless without them"
            )
        for m, n, coefficient in self.monomials:
            if m < 0 or n < 0:
                raise ValueError("monomial exponents must be nonnegative")
            if not isinstance(coefficient, Fraction):
                raise ValueError("monomial coefficients must be Fractions")

    @property
    def poly(self) -> _Poly:
        return {(m, n): c for m, n, c in self.monomials}

    @staticmethod
    def from_triple(poly, alpha, beta) -> "GaussianTerm":
        """Build a term from a ``(poly_dict, alpha, beta)`` triple.

        ``float`` inputs convert exactly (every binary64 is a dyadic
        rational), so nothing is approximated on the way in.
        """
        canonical = tuple(
            sorted((m, n, c) for (m, n), c in _poly_dict(poly).items())
        )
        return GaussianTerm(
            monomials=canonical, alpha=Fraction(alpha), beta=Fraction(beta)
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "poly": [[m, n, str(c)] for m, n, c in self.monomials],
            "alpha": str(self.alpha),
            "beta": str(self.beta),
        }


def as_gaussian_term(entry) -> GaussianTerm:
    """Adapt anything carrying ``(poly, alpha, beta)`` to a :class:`GaussianTerm`.

    Accepted forms: a :class:`GaussianTerm` (returned unchanged), an object
    with ``poly``/``alpha``/``beta`` attributes (whatever class the symbolic
    ``gaussian_hermite`` module settles on), or a plain
    ``(poly_dict, alpha, beta)`` triple.  This is the whole wiring contract
    between the two modules, so it is deliberately dumb.
    """
    if isinstance(entry, GaussianTerm):
        return entry
    if hasattr(entry, "polynomial") and hasattr(entry, "alpha"):
        # gaussian_hermite.GaussianGenerator: an exact polynomial plus a float
        # amplitude.  The amplitude folds into the coefficients exactly, every
        # binary64 being a dyadic rational.
        amplitude = Fraction(float(getattr(entry, "amplitude", 1.0)))
        poly = {
            key: amplitude * Fraction(value)
            for key, value in entry.polynomial.coefficients.items()
        }
        return GaussianTerm.from_triple(poly, entry.alpha, entry.beta)
    if hasattr(entry, "poly") and hasattr(entry, "alpha") and hasattr(entry, "beta"):
        return GaussianTerm.from_triple(entry.poly, entry.alpha, entry.beta)
    poly, alpha, beta = entry
    return GaussianTerm.from_triple(poly, alpha, beta)


@dataclass(frozen=True)
class GaussianFamily:
    r"""A swirl generator and a stream generator, each a sum of Gaussian terms.

    The recovery convention is the audited one shared with
    :class:`~ns_certificate_lab.mixed_initial_data.MixedFamily`:
    ``u^\theta = r\,u_1``, ``u^r = -r\,\partial_z\psi_1``,
    ``u^z = 2\psi_1 + 2s\,\partial_s\psi_1`` with ``s = r^2``, which is
    divergence free for any ``\psi_1``.  ``clay_admissible`` is ``True``
    because a polynomial times a Gaussian is Schwartz, hence ``C^\infty`` —
    but, as everywhere in this repository, admissibility as a *datum* claims
    nothing about the Clay problem.
    """

    name: str
    swirl: tuple[GaussianTerm, ...]
    stream: tuple[GaussianTerm, ...]
    clay_admissible: bool = True

    def __post_init__(self) -> None:
        if not self.swirl and not self.stream:
            raise ValueError("a family needs at least one Gaussian term")

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "basis": "gaussian_hermite",
            "swirl": [term.as_dict() for term in self.swirl],
            "stream": [term.as_dict() for term in self.stream],
        }


def seed_gaussian_family(
    swirl_amplitude: float | Fraction = 10,
    stream_amplitude: float | Fraction = 1,
) -> GaussianFamily:
    r"""The Gaussian analogue of family M1: one swirl term, one stream term.

    ``u_1 = A_s\,z\,e^{-s-z^2}`` and ``\psi_1 = A_m\,z\,e^{-s-z^2}``.  The
    stream generator is **odd in z**, which is the parity the selection rule of
    :func:`~ns_certificate_lab.l3_generation.parity_selection_rule` requires
    for the pressure term to survive integration at all; an even stream term
    here would rebuild a family that is structurally unable to do the one
    thing the certificate measures.
    """
    swirl = GaussianTerm.from_triple(
        {(0, 1): Fraction(swirl_amplitude)}, Fraction(1), Fraction(1)
    )
    stream = GaussianTerm.from_triple(
        {(0, 1): Fraction(stream_amplitude)}, Fraction(1), Fraction(1)
    )
    return GaussianFamily(name="G-seed", swirl=(swirl,), stream=(stream,))


def _family_terms(family_or_triples):
    """Normalise a family-like object to two tuples of :class:`GaussianTerm`."""
    if hasattr(family_or_triples, "swirl") and hasattr(family_or_triples, "stream"):
        swirl, stream = family_or_triples.swirl, family_or_triples.stream
    elif isinstance(family_or_triples, dict):
        swirl = family_or_triples.get("swirl", ())
        stream = family_or_triples.get("stream", ())
    else:
        raise TypeError(
            "expected an object with 'swirl' and 'stream' term collections, or "
            "a dict with those keys; adapt raw triples with as_gaussian_term"
        )
    return (
        tuple(as_gaussian_term(entry) for entry in swirl),
        tuple(as_gaussian_term(entry) for entry in stream),
    )


# --------------------------------------------------------------------------- #
# the derivative closure, cached per term                                      #
# --------------------------------------------------------------------------- #
#
# Every partial of P e^{-alpha s - beta z^2} is Q e^{-alpha s - beta z^2} with
# Q computed EXACTLY in Fractions.  The six polynomials are cached against the
# (hashable) term because a grid sweep re-asks for them at every cell.

_DERIVATIVE_CACHE: dict[GaussianTerm, tuple[_Poly, ...]] = {}
_RADIAL_MONOMIAL_CACHE: dict[tuple, Interval] = {}
_AXIAL_MONOMIAL_CACHE: dict[tuple, Interval] = {}
_RADIAL_GAUSS_CACHE: dict[tuple, Interval] = {}
_AXIAL_GAUSS_CACHE: dict[tuple, Interval] = {}


def clear_gaussian_cache() -> None:
    """Empty every cache.  Only tests should need this."""
    _DERIVATIVE_CACHE.clear()
    _RADIAL_MONOMIAL_CACHE.clear()
    _AXIAL_MONOMIAL_CACHE.clear()
    _RADIAL_GAUSS_CACHE.clear()
    _AXIAL_GAUSS_CACHE.clear()


def _term_derivative_polys(term: GaussianTerm) -> tuple[_Poly, ...]:
    """``(value, ds, dss, dz, dzz, dsz)`` polynomial factors, exactly."""
    cached = _DERIVATIVE_CACHE.get(term)
    if cached is not None:
        return cached
    value = term.poly
    alpha, beta = term.alpha, term.beta

    def d_ds(poly: _Poly) -> _Poly:
        return _poly_add(_poly_diff_s(poly), _poly_scale(poly, -alpha))

    def d_dz(poly: _Poly) -> _Poly:
        return _poly_add(_poly_diff_z(poly), _poly_shift(poly, 0, 1, -2 * beta))

    ds = d_ds(value)
    dz = d_dz(value)
    polys = (value, ds, d_ds(ds), dz, d_dz(dz), d_dz(ds))
    _DERIVATIVE_CACHE[term] = polys
    return polys


# --------------------------------------------------------------------------- #
# sharp interval evaluation on a box                                           #
# --------------------------------------------------------------------------- #


def _s_power(s_box: Interval, m: int) -> Interval:
    """``s^m`` over ``s_box`` with ``s >= 0`` — sharp by monotonicity."""
    if m == 0:
        return _ONE
    key = (s_box.lower, s_box.upper, m)
    cached = _RADIAL_MONOMIAL_CACHE.get(key)
    if cached is None:
        cached = Interval(s_box.lower**m, s_box.upper**m)
        _RADIAL_MONOMIAL_CACHE[key] = cached
    return cached


def _z_power(z_box: Interval, n: int) -> Interval:
    """``z^n`` over ``z_box`` — sharp using the parity of ``n``.

    Odd powers are increasing on all of the line; even powers are functions of
    ``|z|``, whose minimum over the box is zero exactly when the box straddles
    the origin.  Either way the range is attained at explicit points, so no
    generic interval power (with its dependency widening) is ever used.
    """
    if n == 0:
        return _ONE
    key = (z_box.lower, z_box.upper, n)
    cached = _AXIAL_MONOMIAL_CACHE.get(key)
    if cached is not None:
        return cached
    if n % 2 == 1:
        result = Interval(z_box.lower**n, z_box.upper**n)
    else:
        high = max(abs(z_box.lower), abs(z_box.upper)) ** n
        low = (
            Fraction(0)
            if z_box.contains_zero
            else min(abs(z_box.lower), abs(z_box.upper)) ** n
        )
        result = Interval(low, high)
    _AXIAL_MONOMIAL_CACHE[key] = result
    return result


def interval_poly_eval(
    poly, s_box: Interval, z_box: Interval, bits: int = DEFAULT_PRECISION_BITS
) -> Interval:
    r"""Enclose ``P(s, z)`` over a box, with the monomial-sharp improvement.

    Two enclosures are computed and intersected — both are valid, so the
    intersection is valid and never worse than either:

    * **monomial sum**: each ``s^mz^n`` is enclosed *sharply* (``s^m`` by
      monotonicity on ``s \ge 0``, ``z^n`` by parity), and the products
      ``[s^m][z^n]`` are exact ranges because the two factors depend on
      different variables.  The only widening is the sum across monomials,
      which for the low-degree polynomials of the derivative closure is small;
    * **Horner in s** with interval ``z``-powers, which can be tighter when
      monomials of opposite sign share a power of ``s`` and the cancellation
      happens inside the scheme.

    Rounded outward at ``bits`` so denominators stay dyadic and bounded.
    """
    if s_box.lower < 0:
        raise ValueError("s = r^2 cannot be negative; the s-box is out of range")
    monomials = list(_iter_monomials(poly))
    if not monomials:
        return _ZERO

    direct = _ZERO
    by_degree: dict[int, Interval] = {}
    for m, n, coefficient in monomials:
        monomial = (_s_power(s_box, m) * _z_power(z_box, n)).scale(coefficient)
        direct = direct + monomial
        axial = _z_power(z_box, n).scale(coefficient)
        by_degree[m] = by_degree.get(m, _ZERO) + axial

    top = max(by_degree)
    horner = by_degree[top]
    for m in range(top - 1, -1, -1):
        horner = _round(horner * s_box, bits) + by_degree.get(m, _ZERO)

    lower = max(direct.lower, horner.lower)
    upper = min(direct.upper, horner.upper)
    if lower > upper:  # both enclose the range, so this cannot happen; be safe
        return _round(direct, bits)
    return _round(Interval(lower, upper), bits)


def _radial_gauss(
    alpha: Fraction, s_box: Interval, *, terms: int, bits: int
) -> Interval:
    r"""``e^{-\alpha s}`` over ``s_box``, sharp by monotonicity, cached."""
    key = (alpha, s_box.lower, s_box.upper, terms, bits)
    cached = _RADIAL_GAUSS_CACHE.get(key)
    if cached is None:
        exponent = Interval(-alpha * s_box.upper, -alpha * s_box.lower)
        cached = _round(exp_interval(exponent, terms=terms), bits)
        _RADIAL_GAUSS_CACHE[key] = cached
    return cached


def _axial_gauss(
    beta: Fraction, z_box: Interval, *, terms: int, bits: int
) -> Interval:
    r"""``e^{-\beta z^2}`` over ``z_box``, sharp through the parity of ``z^2``."""
    key = (beta, z_box.lower, z_box.upper, terms, bits)
    cached = _AXIAL_GAUSS_CACHE.get(key)
    if cached is None:
        squared = _z_power(z_box, 2)
        exponent = Interval(-beta * squared.upper, -beta * squared.lower)
        cached = _round(exp_interval(exponent, terms=terms), bits)
        _AXIAL_GAUSS_CACHE[key] = cached
    return cached


def gaussian_box(
    alpha,
    beta,
    s_box: Interval,
    z_box: Interval,
    *,
    terms: int = 32,
    bits: int = DEFAULT_PRECISION_BITS,
) -> Interval:
    r"""Enclose ``e^{-\alpha s - \beta z^2}`` over a box, sharply.

    The exponent is monotone decreasing in ``s`` and in ``|z|``, so the range
    of the Gaussian over the box is attained at explicit corner points and
    there is **no dependency problem at all** — this is the entire reason the
    Gaussian basis rescues the certificate from the flat-bump widths.

    Implemented as the product of the two one-dimensional factors
    ``e^{-\alpha s}`` and ``e^{-\beta z^2}``: because the factors depend on
    different variables and are positive, the product of their sharp
    enclosures is again the exact range (up to the two outward roundings), and
    the one-dimensional factors are cached so a grid sweep pays for
    ``O(n_r + n_z)`` exponential series rather than ``O(n_r n_z)`` — the same
    trade recorded in :mod:`ns_certificate_lab.interval_generators`, where the
    uncached version was measured at 209 ms per cell.
    """
    if s_box.lower < 0:
        raise ValueError("s = r^2 cannot be negative; the s-box is out of range")
    alpha, beta = Fraction(alpha), Fraction(beta)
    if alpha <= 0 or beta <= 0:
        raise ValueError("alpha and beta must be strictly positive")
    radial = _radial_gauss(alpha, s_box, terms=terms, bits=bits)
    axial = _axial_gauss(beta, z_box, terms=terms, bits=bits)
    return _round(radial * axial, bits)


def _term_intervals(
    term: GaussianTerm,
    s_box: Interval,
    z_box: Interval,
    *,
    terms: int,
    bits: int,
) -> GeneratorIntervals:
    """Enclose one term's value and five partials over the cell."""
    polys = _term_derivative_polys(term)
    gauss = gaussian_box(
        term.alpha, term.beta, s_box, z_box, terms=terms, bits=bits
    )
    return GeneratorIntervals(
        *(
            _round(interval_poly_eval(poly, s_box, z_box, bits=bits) * gauss, bits)
            for poly in polys
        )
    )


# --------------------------------------------------------------------------- #
# the cell enclosure                                                           #
# --------------------------------------------------------------------------- #


def cell_enclosure_gaussian(
    family_or_triples,
    r_box: Interval,
    z_box: Interval,
    *,
    terms: int = 32,
    bits: int = DEFAULT_PRECISION_BITS,
) -> CellEnclosure:
    r"""Enclose everything ``J`` needs over one cell, Gaussian--Hermite basis.

    The output contract is byte-identical to
    :func:`ns_certificate_lab.interval_generators.cell_enclosure` — the same
    :class:`~ns_certificate_lab.interval_generators.CellEnclosure` dataclass
    with the same fields — so the certificate builder consumes either backend
    without knowing which basis produced the numbers.  The final assembly
    below (gradient components, speed, the Kato-envelope flux intersection) is
    copied from the flat-bump version verbatim: the two backends differ only
    in how the six generator partials are enclosed, and sharing the assembly
    is what keeps the contract honest.
    """
    swirl_terms, stream_terms = _family_terms(family_or_triples)
    if r_box.lower < 0:
        raise ValueError("the radial box cannot extend below the axis")
    s_box = _round(square(r_box), bits)

    swirl = GeneratorIntervals.zero()
    for term in swirl_terms:
        swirl = swirl + _term_intervals(term, s_box, z_box, terms=terms, bits=bits)
    stream = GeneratorIntervals.zero()
    for term in stream_terms:
        stream = stream + _term_intervals(term, s_box, z_box, terms=terms, bits=bits)

    two_s = s_box.scale(Fraction(2))

    u_theta = _round(r_box * swirl.value, bits)
    u_r = _round(-(r_box * stream.dz), bits)
    u_z = _round(stream.value.scale(Fraction(2)) + _round(two_s * stream.ds, bits), bits)

    gradient = {
        "rr": -(stream.dz + _round(two_s * stream.dsz, bits)),
        "rt": swirl.value + _round(two_s * swirl.ds, bits),
        "rz": _round(
            r_box.scale(Fraction(2))
            * _round(stream.ds.scale(Fraction(4)) + _round(two_s * stream.dss, bits),
                     bits),
            bits,
        ),
        "tr": -swirl.value,
        "tt": -stream.dz,
        "tz": _ZERO,
        "zr": -_round(r_box * stream.dzz, bits),
        "zt": _round(r_box * swirl.dz, bits),
        "zz": stream.dz.scale(Fraction(2)) + _round(two_s * stream.dsz, bits),
    }
    gradient = {key: _round(value, bits) for key, value in gradient.items()}

    speed_squared = _round(
        _round(square(u_r), bits) + _round(square(u_theta), bits)
        + _round(square(u_z), bits),
        bits,
    )
    speed = Interval(
        max(speed_squared.lower, Fraction(0)), max(speed_squared.upper, Fraction(0))
    )
    speed = _round(sqrt_interval(speed, bits=bits), bits)
    gradient_squared = _round(
        sum((_round(square(v), bits) for v in gradient.values()), _ZERO), bits
    )

    # |V| integrand.  Kato gives |grad |u|| <= |grad u| pointwise, so
    # q(|grad u|^2 + |grad q|^2) <= 2 q |grad u|^2, which needs no division by q
    # and is therefore valid on the zero set as well.
    viscous = _round((speed * gradient_squared).scale(Fraction(2)), bits)

    # |g| = |u . grad |u|| <= |u| |grad u|, with no division, so it is valid on
    # the zero set too.  This is the Kato bound again.
    magnitude = _round(
        speed * _round(sqrt_interval(gradient_squared, bits=bits), bits), bits
    )
    envelope = Interval(-magnitude.upper, magnitude.upper)

    # The signed flux g = u . grad |u| = (u^r G_r + u^z G_z)/|u| with
    # G = (1/2) grad(|u|^2), which is smooth.  When the speed enclosure stays
    # away from zero the quotient is taken directly and is much tighter than the
    # Kato envelope; when it straddles zero the quotient does not exist and the
    # envelope is the only rigorous statement.  Intersecting the two is free and
    # never worse.
    g_r = _round(
        _round(u_r * gradient["rr"], bits) + _round(u_theta * gradient["rt"], bits)
        + _round(u_z * gradient["rz"], bits),
        bits,
    )
    g_z = _round(
        _round(u_r * gradient["zr"], bits) + _round(u_theta * gradient["zt"], bits)
        + _round(u_z * gradient["zz"], bits),
        bits,
    )
    numerator = _round(
        _round(u_r * g_r, bits) + _round(u_z * g_z, bits), bits
    )
    if speed.lower > 0:
        quotient = _round(divide(numerator, speed), bits)
        flux = Interval(
            max(quotient.lower, envelope.lower), min(quotient.upper, envelope.upper)
        )
    else:
        flux = envelope

    divergence = _round(gradient["rr"] + gradient["tt"] + gradient["zz"], bits)
    return CellEnclosure(
        speed_squared=speed_squared,
        speed=speed,
        gradient_squared=gradient_squared,
        viscous_integrand=viscous,
        flux_magnitude=Interval(Fraction(0), max(magnitude.upper, Fraction(0))),
        flux=flux,
        divergence=divergence,
    )


@dataclass(frozen=True)
class VelocityCellEnclosure:
    r"""Component-level enclosures over one cell, for pairings that need the
    velocity itself rather than only the assembled ``J`` integrands.

    The gradient-form pressure pairing of
    :mod:`ns_certificate_lab.gaussian_gradient_certificate` integrates
    ``-3\,|u|\,(u^r\partial_r p_h + u^z\partial_z p_h)`` and therefore needs
    rigorous cell enclosures of ``u^r``, ``u^z`` and ``|u|`` separately —
    :class:`~ns_certificate_lab.interval_generators.CellEnclosure` deliberately
    exposes only the assembled integrands.  The assembly here is byte-for-byte
    the one inside :func:`cell_enclosure_gaussian`; only the returned fields
    differ.
    """

    u_r: Interval
    u_theta: Interval
    u_z: Interval
    speed_squared: Interval
    speed: Interval
    gradient_squared: Interval
    viscous_integrand: Interval

    def as_dict(self) -> dict[str, list[str]]:
        return {
            name: getattr(self, name).as_pair() for name in self.__dataclass_fields__
        }


def cell_velocity_enclosure(
    family_or_triples,
    r_box: Interval,
    z_box: Interval,
    *,
    terms: int = 32,
    bits: int = DEFAULT_PRECISION_BITS,
) -> VelocityCellEnclosure:
    r"""Enclose the velocity components and the viscous integrand over a cell.

    Same Gaussian--Hermite mechanics as :func:`cell_enclosure_gaussian` — the
    six exact derivative polynomials per term, the sharp Gaussian box, the
    audited recovery ``u^\theta = r u_1``, ``u^r = -r\partial_z\psi_1``,
    ``u^z = 2\psi_1 + 2s\partial_s\psi_1`` — but returning the *components*.
    The viscous integrand ``2|u||\nabla u|^2`` (the Kato-bounded upper form)
    is computed in the same pass so a gradient-form sweep needs one cell
    evaluation, not two.
    """
    swirl_terms, stream_terms = _family_terms(family_or_triples)
    if r_box.lower < 0:
        raise ValueError("the radial box cannot extend below the axis")
    s_box = _round(square(r_box), bits)

    swirl = GeneratorIntervals.zero()
    for term in swirl_terms:
        swirl = swirl + _term_intervals(term, s_box, z_box, terms=terms, bits=bits)
    stream = GeneratorIntervals.zero()
    for term in stream_terms:
        stream = stream + _term_intervals(term, s_box, z_box, terms=terms, bits=bits)

    two_s = s_box.scale(Fraction(2))
    u_theta = _round(r_box * swirl.value, bits)
    u_r = _round(-(r_box * stream.dz), bits)
    u_z = _round(
        stream.value.scale(Fraction(2)) + _round(two_s * stream.ds, bits), bits
    )

    gradient = {
        "rr": -(stream.dz + _round(two_s * stream.dsz, bits)),
        "rt": swirl.value + _round(two_s * swirl.ds, bits),
        "rz": _round(
            r_box.scale(Fraction(2))
            * _round(stream.ds.scale(Fraction(4)) + _round(two_s * stream.dss, bits),
                     bits),
            bits,
        ),
        "tr": -swirl.value,
        "tt": -stream.dz,
        "tz": _ZERO,
        "zr": -_round(r_box * stream.dzz, bits),
        "zt": _round(r_box * swirl.dz, bits),
        "zz": stream.dz.scale(Fraction(2)) + _round(two_s * stream.dsz, bits),
    }
    gradient = {key: _round(value, bits) for key, value in gradient.items()}

    speed_squared = _round(
        _round(square(u_r), bits) + _round(square(u_theta), bits)
        + _round(square(u_z), bits),
        bits,
    )
    speed_squared = Interval(
        max(speed_squared.lower, Fraction(0)), max(speed_squared.upper, Fraction(0))
    )
    speed = _round(sqrt_interval(speed_squared, bits=bits), bits)
    gradient_squared = _round(
        sum((_round(square(v), bits) for v in gradient.values()), _ZERO), bits
    )
    viscous = _round((speed * gradient_squared).scale(Fraction(2)), bits)
    return VelocityCellEnclosure(
        u_r=u_r,
        u_theta=u_theta,
        u_z=u_z,
        speed_squared=speed_squared,
        speed=speed,
        gradient_squared=gradient_squared,
        viscous_integrand=viscous,
    )


def gaussian_reference_field(family, grid):
    """Public entry to the binary64 reference field (windowed to the box).

    Nothing rigorous consumes the result — it exists so the gradient-form
    certificate builder can produce the SAME discrete pressure ``p_h`` as
    :func:`build_gaussian_generation_certificate` without reaching into a
    private helper.  See :func:`_float_field` for the windowing rationale.
    """
    return _float_field(family, grid)


# --------------------------------------------------------------------------- #
# rigorous tail bounds outside the box                                         #
# --------------------------------------------------------------------------- #


def gaussian_tail_moment(
    gamma,
    R,
    m: int,
    *,
    terms: int = 32,
    bits: int = DEFAULT_PRECISION_BITS,
) -> Fraction:
    r"""A rational upper bound on ``\int_R^\infty \rho^m e^{-2\gamma\rho^2}d\rho``.

    Integration by parts, ``\rho^m = \rho^{m-1}\cdot\rho`` with
    ``\rho\,e^{-2\gamma\rho^2} = -\bigl(e^{-2\gamma\rho^2}/(4\gamma)\bigr)'``,
    gives the exact descent

    .. math::
       I_m(R) = \frac{R^{m-1}e^{-2\gamma R^2}}{4\gamma}
              + \frac{m-1}{4\gamma}\,I_{m-2}(R),

    which terminates at ``I_1(R) = e^{-2\gamma R^2}/(4\gamma)`` (exact) and at
    the **Mills-ratio bound**

    .. math::
       I_0(R) = \int_R^\infty e^{-2\gamma\rho^2}\,d\rho
              \le \frac{e^{-2\gamma R^2}}{4\gamma R}\qquad(R > 0),

    proved by inserting ``1 \le \rho/R`` into the integrand and integrating
    ``4\gamma\rho\,e^{-2\gamma\rho^2}`` exactly.  The positivity of ``R`` is a
    genuine hypothesis of the Mills bound and is enforced.  Replacing
    ``e^{-2\gamma R^2}`` by its series **upper** bound keeps every step an
    upper bound because all coefficients in the descent are positive.
    """
    gamma, R = Fraction(gamma), Fraction(R)
    if gamma <= 0:
        raise ValueError("gamma must be strictly positive")
    if R <= 0:
        raise ValueError(
            "the Mills-ratio bound needs R > 0; there is no finite bound at R = 0 "
            "for the m = 0 moment from this argument"
        )
    if m < 0 or int(m) != m:
        raise ValueError("the moment order m must be a nonnegative integer")
    m = int(m)

    point = -2 * gamma * R * R
    envelope = exp_interval(Interval(point, point), terms=terms).upper
    moments = [envelope / (4 * gamma * R), envelope / (4 * gamma)]
    for k in range(2, m + 1):
        moments.append(
            R ** (k - 1) * envelope / (4 * gamma)
            + Fraction(k - 1, 1) / (4 * gamma) * moments[k - 2]
        )
    return Interval(Fraction(0), moments[min(m, len(moments) - 1)]).round_outward(
        bits
    ).upper


def _sup_power_gaussian(
    j: int, a: Fraction, *, terms: int, bits: int
) -> Fraction:
    r"""A rational upper bound on ``\sup_{x\ge0} x^j e^{-a x^2}`` for ``a > 0``.

    Substituting ``t = x^2`` reduces it to the calculus bound

    .. math::  \sup_{t\ge0} t^\mu e^{-at} = \Bigl(\frac{\mu}{ae}\Bigr)^\mu
       \qquad(\mu = j/2,\ \text{attained at } t = \mu/a),

    so the supremum is ``(j/(2ae))^{j/2}``.  A rational **upper** bound of
    ``1/e`` comes from the audited series (``exp_interval`` at ``-1``); the
    base is then an upper bound of a nonnegative quantity and raising it to a
    positive power preserves the inequality.  Odd ``j`` needs a square root,
    taken with the outward :func:`sqrt_interval`.
    """
    if a <= 0:
        raise ValueError("the exponential rate must be strictly positive")
    if j == 0:
        return Fraction(1)
    inverse_e = exp_interval(
        Interval(Fraction(-1), Fraction(-1)), terms=terms
    ).upper
    base = Fraction(j) * inverse_e / (2 * a)
    if j % 2 == 0:
        value = base ** (j // 2)
    else:
        value = sqrt_interval(Interval(base**j, base**j), bits=bits).upper
    return Interval(Fraction(0), value).round_outward(bits).upper


@dataclass(frozen=True)
class ExteriorTailBound:
    r"""Rational upper bounds on the exterior integrals of the certificate.

    ``speed_coefficient`` and ``gradient_coefficient`` are the ``C_u`` and
    ``C_g`` of the global envelopes ``|u| \le C_u e^{-\gamma\rho^2}`` and
    ``|\nabla u| \le C_g e^{-\gamma\rho^2}`` (valid on **all** of space, hence
    in particular outside the box).  The three tails bound, over the exterior
    of the box:

    * ``viscous_tail``:   ``\int 2|u||\nabla u|^2\,dx``,
    * ``flux_tail``:      ``\int |u||\nabla u|\,dx`` — which dominates
      ``\int|\nabla\cdot(|u|u)|\,dx`` by the Kato bound, and is therefore the
      one folded into the certificate's flux integral,
    * ``cubed_flux_tail``: ``\int |u|^2|\nabla u|\,dx``, the cubic-order
      companion, recorded because the ``L^3`` bookkeeping asks for it.
    """

    gamma: Fraction
    distance: Fraction
    speed_coefficient: Fraction
    gradient_coefficient: Fraction
    viscous_tail: Fraction
    flux_tail: Fraction
    cubed_flux_tail: Fraction

    def as_dict(self) -> dict[str, str]:
        return {name: str(getattr(self, name)) for name in self.__dataclass_fields__}


def _pieces_sup(
    pieces, gamma: Fraction, *, terms: int, bits: int
) -> Fraction:
    r"""``\sup_{\mathbb R^3} |Q| e^{\gamma\rho^2}`` for a sum of basis pieces.

    Each piece is ``(k, P, \alpha, \beta)`` meaning ``r^kP(s,z)e^{-\alpha s
    -\beta z^2}`` with ``k \in \{0, 1\}``.  Splitting off ``e^{-\gamma\rho^2}``
    leaves ``r^{2m+k}e^{-(\alpha-\gamma)r^2}`` and ``|z|^ne^{-(\beta-\gamma)z^2}``
    per monomial, each bounded by :func:`_sup_power_gaussian`; the triangle
    inequality over monomials and pieces is generous but exact-rational and
    finite, which is all the tail needs.
    """
    total = Fraction(0)
    for k, poly, alpha, beta in pieces:
        residual_r = alpha - gamma
        residual_z = beta - gamma
        if residual_r <= 0 or residual_z <= 0:
            raise ValueError("gamma must leave residual decay in both coordinates")
        for m, n, coefficient in _iter_monomials(poly):
            total += (
                abs(coefficient)
                * _sup_power_gaussian(2 * m + k, residual_r, terms=terms, bits=bits)
                * _sup_power_gaussian(n, residual_z, terms=terms, bits=bits)
            )
    return Interval(Fraction(0), total).round_outward(bits).upper


def _sqrt_upper(value: Fraction, bits: int) -> Fraction:
    return sqrt_interval(Interval(value, value), bits=bits).upper


def exterior_viscous_tail(
    family_or_triples,
    r_max,
    z_max,
    *,
    terms: int = 32,
    bits: int = DEFAULT_PRECISION_BITS,
) -> ExteriorTailBound:
    r"""Bound the exterior contributions to the certificate integrals.

    Outside the box ``\{r \le R_{\max},\ |z| \le Z_{\max}\}`` every point has
    ``\rho^2 = r^2 + z^2 > d^2`` with ``d = \min(R_{\max}, Z_{\max})``, so the
    exterior is contained in the spherical shell ``\rho \ge d`` and the
    physical measure is bounded by the shell measure ``4\pi\rho^2\,d\rho``
    (an over-count of the meridional ``2\pi r\,dr\,dz`` measure, hence safe).
    With the global envelopes ``|u| \le C_ue^{-\gamma\rho^2}``,
    ``|\nabla u| \le C_ge^{-\gamma\rho^2}``:

    .. math::
       \int_{\rm ext} 2|u||\nabla u|^2\,dx
         \le 2C_uC_g^2\,4\pi\int_d^\infty \rho^2e^{-3\gamma\rho^2}d\rho ,
       \qquad
       \int_{\rm ext} |u||\nabla u|\,dx
         \le C_uC_g\,4\pi\int_d^\infty \rho^2e^{-2\gamma\rho^2}d\rho ,

    both closed by :func:`gaussian_tail_moment` (the cubic tail uses
    ``\gamma' = 3\gamma/2`` so that ``e^{-2\gamma'\rho^2} = e^{-3\gamma\rho^2}``).
    ``\gamma`` is half the weakest decay rate in the family, which leaves
    residual decay ``\alpha - \gamma \ge \alpha/2 > 0`` in each coordinate for
    the polynomial suprema.  The result decays like ``e^{-3\gamma d^2}``:
    generous, finite, and shrinking fast as the box grows.
    """
    swirl_terms, stream_terms = _family_terms(family_or_triples)
    every = (*swirl_terms, *stream_terms)
    if not every:
        raise ValueError("cannot bound the tail of an empty family")
    r_max, z_max = Fraction(r_max), Fraction(z_max)
    if r_max <= 0 or z_max <= 0:
        raise ValueError("the box must have positive extent in both directions")
    gamma = min(min(term.alpha, term.beta) for term in every) / 2
    distance = min(r_max, z_max)

    # Decompose each velocity and gradient component into (k, poly, alpha,
    # beta) pieces, k the explicit power of r in front.  These are the SAME
    # closed forms the cell enclosure assembles, so the envelope bounds the
    # identical quantities the interior enclosures measure.
    velocity: dict[str, list] = {"u_r": [], "u_t": [], "u_z": []}
    gradient: dict[str, list] = {
        key: [] for key in ("rr", "rt", "rz", "tr", "tt", "zr", "zt", "zz")
    }
    for term in swirl_terms:
        value, ds, _dss, dz, _dzz, _dsz = _term_derivative_polys(term)
        pair = (term.alpha, term.beta)
        velocity["u_t"].append((1, value, *pair))
        gradient["rt"].append(
            (0, _poly_add(value, _poly_shift(ds, 1, 0, Fraction(2))), *pair)
        )
        gradient["tr"].append((0, value, *pair))
        gradient["zt"].append((1, dz, *pair))
    for term in stream_terms:
        value, ds, dss, dz, dzz, dsz = _term_derivative_polys(term)
        pair = (term.alpha, term.beta)
        velocity["u_r"].append((1, dz, *pair))
        velocity["u_z"].append(
            (0, _poly_add(_poly_scale(value, Fraction(2)),
                          _poly_shift(ds, 1, 0, Fraction(2))), *pair)
        )
        rr = _poly_add(dz, _poly_shift(dsz, 1, 0, Fraction(2)))
        gradient["rr"].append((0, rr, *pair))
        gradient["zz"].append((0, rr, *pair))
        gradient["tt"].append((0, dz, *pair))
        gradient["zr"].append((1, dzz, *pair))
        gradient["rz"].append(
            (1, _poly_add(_poly_scale(ds, Fraction(8)),
                          _poly_shift(dss, 1, 0, Fraction(4))), *pair)
        )

    speed_sq = sum(
        _pieces_sup(pieces, gamma, terms=terms, bits=bits) ** 2
        for pieces in velocity.values()
    )
    gradient_sq = sum(
        _pieces_sup(pieces, gamma, terms=terms, bits=bits) ** 2
        for pieces in gradient.values()
    )
    speed_coefficient = _sqrt_upper(speed_sq, bits)
    gradient_coefficient = _sqrt_upper(gradient_sq, bits)

    cubic_moment = gaussian_tail_moment(
        3 * gamma / 2, distance, 2, terms=terms, bits=bits
    )
    quadratic_moment = gaussian_tail_moment(
        gamma, distance, 2, terms=terms, bits=bits
    )

    def _up(value: Fraction) -> Fraction:
        return Interval(Fraction(0), value).round_outward(bits).upper

    viscous_tail = _up(
        2 * speed_coefficient * gradient_coefficient**2 * _FOUR_PI_UPPER
        * cubic_moment
    )
    flux_tail = _up(
        speed_coefficient * gradient_coefficient * _FOUR_PI_UPPER
        * quadratic_moment
    )
    cubed_flux_tail = _up(
        speed_coefficient**2 * gradient_coefficient * _FOUR_PI_UPPER
        * cubic_moment
    )
    return ExteriorTailBound(
        gamma=gamma,
        distance=distance,
        speed_coefficient=speed_coefficient,
        gradient_coefficient=gradient_coefficient,
        viscous_tail=viscous_tail,
        flux_tail=flux_tail,
        cubed_flux_tail=cubed_flux_tail,
    )


# --------------------------------------------------------------------------- #
# floating-point evaluation, for the discrete pressure and the reference       #
# --------------------------------------------------------------------------- #


def _poly_float(poly: _Poly, s: FloatArray, z: FloatArray) -> FloatArray:
    total = np.zeros(np.broadcast(s, z).shape, dtype=np.float64)
    for m, n, coefficient in _iter_monomials(poly):
        total = total + float(coefficient) * s**m * z**n
    return total


def _float_partials(
    terms_list, r_mesh: FloatArray, z_mesh: FloatArray
) -> dict[str, FloatArray]:
    """Summed binary64 partials of a generator group — reference only.

    Nothing rigorous consumes these: they feed the discrete pressure solve
    (whose enclosure is the multilinear corner hull, a theorem) and the
    floating-point witness block of the payload.
    """
    s = r_mesh**2
    totals = {key: np.zeros(r_mesh.shape, dtype=np.float64) for key in _PARTIAL_KEYS}
    for term in terms_list:
        weight = np.exp(-(float(term.alpha) * s + float(term.beta) * z_mesh**2))
        for key, poly in zip(_PARTIAL_KEYS, _term_derivative_polys(term)):
            totals[key] = totals[key] + _poly_float(poly, s, z_mesh) * weight
    return totals


def _float_field(family, grid):
    """The :class:`~ns_certificate_lab.l3_generation.MixedField` in binary64.

    The generator partials are **windowed to the box** — zeroed on the outer
    radial ring and the two axial faces — before the field is assembled.  A
    Gaussian datum is Schwartz, not compactly supported, and the exterior
    representation behind
    :meth:`~ns_certificate_lab.l3_generation.AxisymmetricPressureSolver.exterior_trace`
    places its Dirichlet targets on those very faces: a source node sitting on
    a target produces a zero distance and a NaN.  Windowing is the same
    convention :meth:`~ns_certificate_lab.mixed_initial_data.MixedFamily.initial_state`
    applies, and it costs nothing rigorous: the float field feeds only the
    discrete pressure ``p_h`` (whose cell enclosure is the corner-hull theorem
    for *whatever* grid function ``p_h`` is) and the floating-point witness
    block, never an interval.
    """
    from .l3_generation import MixedField

    swirl_terms, stream_terms = _family_terms(family)
    r_mesh, z_mesh = grid.mesh()
    swirl = _float_partials(swirl_terms, r_mesh, z_mesh)
    stream = _float_partials(stream_terms, r_mesh, z_mesh)
    for group in (swirl, stream):
        for key in _PARTIAL_KEYS:
            group[key][-1, :] = 0.0
            group[key][:, 0] = 0.0
            group[key][:, -1] = 0.0
    s = r_mesh**2
    u_theta = r_mesh * swirl["value"]
    u_r = -r_mesh * stream["dz"]
    u_z = 2.0 * stream["value"] + 2.0 * s * stream["ds"]
    gradient = {
        "rr": -(stream["dz"] + 2.0 * s * stream["dsz"]),
        "rt": swirl["value"] + 2.0 * s * swirl["ds"],
        "rz": 2.0 * r_mesh * (4.0 * stream["ds"] + 2.0 * s * stream["dss"]),
        "tr": -swirl["value"],
        "tt": -stream["dz"],
        "tz": np.zeros_like(r_mesh),
        "zr": -r_mesh * stream["dzz"],
        "zt": r_mesh * swirl["dz"],
        "zz": 2.0 * stream["dz"] + 2.0 * s * stream["dsz"],
    }
    return MixedField(
        grid=grid, u_r=u_r, u_theta=u_theta, u_z=u_z, analytic_gradient=gradient
    )


# --------------------------------------------------------------------------- #
# the certificate                                                              #
# --------------------------------------------------------------------------- #


def build_gaussian_generation_certificate(
    family,
    grid,
    *,
    viscosity: float,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    exp_terms: int = 32,
    subdivisions: int = 1,
):
    r"""Build the generation certificate over the Gaussian--Hermite basis.

    Same payload schema, same single ``P1`` hypothesis and same verifier
    (:func:`~ns_certificate_lab.l3_certificate.verify_generation_certificate`)
    as the flat-bump builder — a checker must not care which basis produced
    the enclosures.  Two things differ, both recorded under ``"rigorous"``:

    * the interior cell enclosures come from :func:`cell_enclosure_gaussian`,
      whose widths shrink linearly with the cell instead of exploding at a
      support edge (there is no support edge);
    * the datum is **not** compactly supported, so the exterior contributions
      to the viscous integral and to ``\int|g|\,dx`` are bounded by
      :func:`exterior_viscous_tail` and **added** to the respective upper
      bounds.  Enlarging an upper bound is always sound: it can only shrink
      the margin and the pressure-error threshold.

    One honesty note, recorded in the payload: the discrete pressure bound
    remains a bound over the **box** only, exactly as in the flat-bump
    certificate.  For a compactly supported datum the exterior pressure term
    is zero identically; here the exterior remainder ``3\int_{\rm ext}p\,g``
    is *not* enclosed — controlling it needs a pressure decay estimate that is
    part of the same unproved ``P1`` circle — so the conditional claim is
    restricted to the box contribution of the pressure integral.

    ``subdivisions`` splits each grid cell before enclosure, as in the
    flat-bump builder: the residual dependency loss shrinks linearly with the
    subdivision width, at linear cost.
    """
    from .grid import AxisymmetricGrid  # noqa: F401  (documents the expected type)
    from .l3_certificate import (
        GenerationCertificate,
        PRESSURE_HYPOTHESIS,
        _cell_hull_from_nodes,
    )
    from .l3_generation import AxisymmetricPressureSolver, l3_generation_rate
    from .l3_optimizer import require_clay_admissible

    require_clay_admissible(family)
    if viscosity <= 0.0:
        raise ValueError("a Clay candidate requires a fixed positive viscosity")
    if grid.periodic_z:
        raise ValueError("the generation certificate targets non-periodic z grids")
    if subdivisions < 1:
        raise ValueError("subdivisions must be at least one")
    if not (float(grid.z[0]) < 0.0 < float(grid.z[-1])):
        raise ValueError(
            "the exterior tail bound needs a box that surrounds the origin: "
            "z_min < 0 < z_max"
        )

    field = _float_field(family, grid)
    solver = AxisymmetricPressureSolver.build(grid)
    rate = l3_generation_rate(field, viscosity=viscosity, solver=solver)

    gradient = field.gradient()
    trace = solver.exterior_trace(field.u_r, field.u_theta, field.u_z)
    pressure = solver.solve(field.pressure_source(gradient), trace=trace)

    # -- rigorous cell enclosures ------------------------------------------- #
    #
    # As in the flat-bump builder: the analytic quantities are enclosed by
    # interval evaluation over each cell, never by a hull of nodal values.
    # The computed pressure is the one exception -- it IS a grid function, and
    # the multilinear interpolant of a grid function has its cell range exactly
    # equal to the corner hull, so there the hull is a theorem.
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
                    enclosure = cell_enclosure_gaussian(
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

    # -- the exterior tails, the new obligation of a non-compact basis ------- #
    tail = exterior_viscous_tail(
        family,
        r_nodes[-1],
        min(-z_nodes[0], z_nodes[-1]),
        terms=exp_terms,
        bits=precision_bits,
    )
    interior_viscous = viscous_upper
    interior_flux = flux_upper
    viscous_upper = viscous_upper + tail.viscous_tail
    flux_upper = flux_upper + tail.flux_tail

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
            "enclosure_basis": "gaussian_hermite",
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
            "interior_viscous_integral_upper_bound": str(interior_viscous),
            "interior_flux_integral_upper_bound": str(interior_flux),
            "exterior_viscous_tail": str(tail.viscous_tail),
            "exterior_flux_tail": str(tail.flux_tail),
            "tail_gamma": str(tail.gamma),
            "tail_distance": str(tail.distance),
            "exterior_pressure_note": (
                "The discrete pressure lower bound covers the box only.  The "
                "exterior remainder 3 int_ext p g dx is NOT enclosed: bounding "
                "it needs a pressure decay estimate that belongs to the same "
                "unproved P1 circle, so the conditional claim is restricted to "
                "the box contribution of the pressure integral.  The exterior "
                "viscous and flux contributions ARE enclosed and folded in."
            ),
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
            "The enclosure basis is Gaussian-Hermite: polynomials in (r^2, z) "
            "times exp(-alpha r^2 - beta z^2), whose derivative closure keeps "
            "every enclosure sharp by monotonicity, and whose exterior "
            "contributions to the viscous and flux integrals are bounded by "
            "explicit Mills-ratio tail integrals and folded into the upper "
            "bounds.",
            "Conditional on P1: if the computed pressure is within the stated "
            "threshold of the true free-space pressure in L^infinity on the box, "
            "then the box contribution to J(u_0) is positive whenever the "
            "margin is.",
            "The advective contribution is exactly zero by an algebraic "
            "identity, not by cancellation of computed numbers.",
            "Nothing here is a singularity result and nothing here bears on the "
            "Clay problem.  A positive initial generation rate of the critical "
            "norm is a screening criterion, not evidence of blow-up.",
        ],
    }
    return GenerationCertificate(payload)
