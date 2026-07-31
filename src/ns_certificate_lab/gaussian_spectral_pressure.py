r"""Closed-form free-space pressure gradient for the Gaussian--Hermite basis.

What this module replaces
-------------------------
Both existing ``J`` certificates pair the velocity against a **discrete**
pressure: the divergence-form builder against the per-cell corner hull of the
computed nodal pressure, the gradient-form builder against the exact gradient
of its multilinear interpolant.  Either way the step from the computed
``p_h`` to the true free-space pressure is the named, *unproved* hypothesis
``P1``/``P1G``.  This module removes the discrete pressure from the certified
chain entirely: for the isotropic Gaussian--Hermite basis the free-space
pressure has a **finite closed form**, and the pressure term of ``J`` is
enclosed against that closed form directly.  No Poisson solve, no interpolant,
no discretisation hypothesis — the box pressure term becomes unconditional.

The closed form
---------------
For the audited recovery ``u^\theta = r\,u_1``, ``u^r = -r\,\partial_z\psi_1``,
``u^z = 2\psi_1 + 2s\,\partial_s\psi_1`` (``s = r^2``), writing ``f =
-\partial_z\psi_1``, ``g = u_1``, ``h = 2\psi_1 + 2s\,\partial_s\psi_1`` — so
that the Cartesian field is ``u = (xf - yg,\; yf + xg,\; h)`` — the pressure
source contracts to a purely meridional expression,

.. math::  \sigma \;=\; \partial_iu_j\,\partial_ju_i
   \;=\; 2f^2 + 4s\,ff_s + 4s^2f_s^2 - 2g^2 - 4s\,gg_s + 4s\,h_sf_z + h_z^2 ,

an algebraic identity checked against the cylindrical trace formula
``\sum_{ab}(\nabla u)_{ab}(\nabla u)_{ba}`` used by
:meth:`~ns_certificate_lab.l3_generation.MixedField.pressure_source` (the
angular dependence cancels exactly; a test pins the two against each other).
For generators with ``\alpha = \beta`` every product envelope is isotropic,
``\sigma = \sum P(s,z)\,e^{-A\rho^2}`` with ``\rho^2 = s + z^2``, and the
free-space pressure ``p = (-\Delta)^{-1}\sigma`` is computed through the heat
representation ``(-\Delta)^{-1} = \int_0^\infty e^{t\Delta}\,dt``:

* the heat flow acts **separably** and maps scaled Laguerre and Hermite
  functions to themselves.  With ``\lambda = 1 + 4At``, ``v = A/\lambda``,
  ``\ell_d(s;A) = A^{-d}L_d(As)`` and ``\hat h_j(z;A) = A^{-j/2}H_j(\sqrt Az)``
  (both have *rational* coefficients),

  .. math:: e^{t\Delta_2}\bigl[\ell_d(s;A)e^{-As}\bigr]
      = \lambda^{-(2d+1)}\ell_d(s;v)\,e^{-vs}, \qquad
      e^{t\partial_z^2}\bigl[\hat h_j(z;A)e^{-Az^2}\bigr]
      = \lambda^{-(2j+1)/2}\hat h_j(z;v)\,e^{-vz^2},

  proved from the Laguerre generating function ``\sum_d L_d(x)y^d =
  (1-y)^{-1}e^{-xy/(1-y)}`` and the Hermite generating function applied to the
  Gaussian translate ``e^{t\Delta}[e^{-A|x-a|^2}]`` (expand both sides in the
  shift ``a`` and match powers);
* the ``t``-integral of every resulting ``v``-monomial is **exact**: with
  ``\tau = 1/\lambda``,

  .. math:: \int_0^\infty \lambda^{-\nu}e^{-A\rho^2/\lambda}\,dt
     \;=\; \frac{1}{4A}\,g_{\nu-2}(A\rho^2), \qquad
     g_\mu(w) \;:=\; \int_0^1 \tau^\mu e^{-w\tau}\,d\tau ,

  so ``p`` is a finite sum of terms ``c\,s^kz^i\,g_\mu(A\rho^2)`` with exact
  rational ``c`` and half-integer ``\mu \ge -1/2``.  ``g_\mu`` is entire,
  positive, decreasing in ``w`` (sharp monotone box enclosures), and the family
  is closed under differentiation: ``g_\mu' = -g_{\mu+1}``, so ``\partial_sp``
  and ``\partial_zp`` are again finite sums of the same shape.  The base case
  reproduces the classical Newtonian potential of a Gaussian,
  ``(-\Delta)^{-1}e^{-A\rho^2} = \sqrt\pi\,\mathrm{erf}(\sqrt A\rho)/
  (4A^{3/2}\rho)``.

The representation is **self-verifying**: :meth:`SpectralPressure.
verify_laplace_identity` reduces ``\Delta_3 p + \sigma`` to zero by exact
rational arithmetic, using only the two-term recursion ``w\,g_\mu(w) =
\mu\,g_{\mu-1}(w) - e^{-w}`` (integration by parts of the defining integral)
and the linear independence of ``\{g_\mu\}\cup\{e^{-w}\}`` over polynomials.
The builder runs this check before certifying anything.

Enclosing ``g_\mu`` rigorously
------------------------------
For ``w`` below a switch point the alternating series ``g_\mu(w) =
\sum_k(-1)^kw^k/(k!\,(\mu+k+1))`` is summed in exact rationals; once the index
exceeds ``w`` the terms decrease, so the first omitted term bounds the
remainder.  The value at the **largest** needed ``\mu`` seeds the downward
recursion ``g_{\mu-1} = (w\,g_\mu + e^{-w})/\mu``, which is a sum of positive
quantities — no cancellation, so interval widths only shrink.  Above the
switch point the seed is ``g_{-1/2}(w) = (\sqrt\pi - 2I(\sqrt w))/\sqrt w``
with ``I(x) = \int_x^\infty e^{-\rho^2}d\rho`` enclosed by the elementary pair

.. math:: \frac{x\,e^{-x^2}}{2x^2+1} \;\le\; I(x) \;\le\; \frac{e^{-x^2}}{2x},

(the upper bound inserts ``1 \le \rho/x`` — the repository's Mills argument —
and the lower bound follows from ``I = e^{-x^2}/(2x) -
\tfrac12\int_x^\infty\rho^{-2}e^{-\rho^2}d\rho`` with ``\rho^{-2}\le x^{-2}``),
then the upward recursion ``g_\mu = (\mu\,g_{\mu-1} - e^{-w})/w``, whose
subtracted term ``e^{-w}`` is exponentially below ``\mu g_{\mu-1} \sim
\Gamma(\mu+1)w^{-\mu}`` — again no cancellation.

What is proved and what is not
------------------------------
Unconditionally proved, in exact rational arithmetic: the enclosure of the box
pressure integral ``-3\int_\Omega |u|\,u\cdot\nabla p\,dx`` for the **true**
free-space pressure, the bound on the exterior pressure remainder
``|{-3\int_{\rm ext}|u|u\cdot\nabla p}|`` (global decay of every ``g``-term
times the Gaussian exterior mass of ``|u|^2``), the viscous upper bound
(interior enclosure plus Mills-ratio tail), hence an **unconditional lower
bound on** ``J(u_0)``.  This certificate carries *no* pressure hypothesis: the
``hypotheses`` block is empty and the verifier enforces that.

Nothing here is a singularity result and nothing here bears on the Clay
problem.  A positive initial generation rate of the critical norm would be a
screening criterion, not evidence of blow-up; a negative margin is evidence of
nothing beyond the width of the enclosure.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .interval_gaussian import (
    _family_terms,
    _iter_monomials,
    _poly_add,
    _poly_diff_s,
    _poly_diff_z,
    _poly_scale,
    _poly_shift,
    _s_power,
    _z_power,
    cell_velocity_enclosure,
    exterior_viscous_tail,
    gaussian_tail_moment,
)
from .interval_generators import divide, square
from .l3_certificate import GenerationCertificate, exp_interval, sqrt_interval
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "SPECTRAL_CERTIFICATE_KIND",
    "SpectralPressure",
    "build_spectral_generation_certificate",
    "clear_spectral_cache",
    "isotropic_pressure_source",
    "verify_spectral_generation_certificate",
]

SPECTRAL_CERTIFICATE_KIND = "l3_generation_rate_spectral_pressure"

_PATHS = ("plain", "taylor", "adaptive")

#: ``\pi`` and ``2\pi`` enclosed between rationals.
_PI = Interval(
    Fraction(3141592653589793, 10**15), Fraction(3141592653589794, 10**15)
)
_TWO_PI = Interval(
    Fraction(6283185307179586, 10**15), Fraction(6283185307179587, 10**15)
)
_ZERO = Interval(Fraction(0), Fraction(0))

#: Switch point between the series and the erf-seeded evaluation of ``g_\mu``.
_G_SERIES_SWITCH = Fraction(30)


def _round(value: Interval, bits: int) -> Interval:
    return value.round_outward(bits)


def _floor_dyadic(value: Fraction, bits: int) -> Fraction:
    scale = Fraction(2) ** bits
    return Fraction(math.floor(value * scale), 1) / scale


def _ceil_dyadic(value: Fraction, bits: int) -> Fraction:
    scale = Fraction(2) ** bits
    return Fraction(math.ceil(value * scale), 1) / scale


# --------------------------------------------------------------------------- #
# the pressure source, purely meridional                                       #
# --------------------------------------------------------------------------- #


def _poly_mul(first: dict, second: dict) -> dict:
    out: dict = {}
    for (m1, n1), c1 in first.items():
        for (m2, n2), c2 in second.items():
            key = (m1 + m2, n1 + n2)
            total = out.get(key, Fraction(0)) + c1 * c2
            if total == 0:
                out.pop(key, None)
            else:
                out[key] = total
    return out


def _closure_d_s(poly: dict, rate: Fraction) -> dict:
    r"""``\partial_s(P e^{-A s - A z^2}) = (P_s - AP)e^{\ldots}`` — the factor."""
    return _poly_add(_poly_diff_s(poly), _poly_scale(poly, -rate))


def _closure_d_z(poly: dict, rate: Fraction) -> dict:
    return _poly_add(_poly_diff_z(poly), _poly_shift(poly, 0, 1, -2 * rate))


def isotropic_pressure_source(family) -> dict[Fraction, dict]:
    r"""``\sigma = \partial_iu_j\partial_ju_i`` as ``{A: P(s,z)}`` with
    ``\sigma = \sum_A P_A\,e^{-A\rho^2}``.

    Requires every generator to be isotropic (``\alpha = \beta``): the closed
    form of the Newtonian potential below needs spherical envelopes, and an
    anisotropic basis would silently produce a wrong pressure.  The honest
    retreat for anisotropic generators is the anisotropic heat calculus, which
    is not implemented; the guard raises instead.
    """
    swirl_terms, stream_terms = _family_terms(family)
    for term in (*swirl_terms, *stream_terms):
        if term.alpha != term.beta:
            raise ValueError(
                "the closed-form spectral pressure requires isotropic "
                "envelopes (alpha == beta for every generator); received "
                f"alpha={term.alpha}, beta={term.beta}"
            )

    f_env: dict[Fraction, dict] = {}
    g_env: dict[Fraction, dict] = {}
    h_env: dict[Fraction, dict] = {}
    for term in stream_terms:
        rate, poly = term.alpha, term.poly
        d_s, d_z = _closure_d_s(poly, rate), _closure_d_z(poly, rate)
        f_env[rate] = _poly_add(f_env.get(rate, {}), _poly_scale(d_z, Fraction(-1)))
        h_env[rate] = _poly_add(
            h_env.get(rate, {}),
            _poly_add(_poly_scale(poly, Fraction(2)),
                      _poly_shift(d_s, 1, 0, Fraction(2))),
        )
    for term in swirl_terms:
        g_env[term.alpha] = _poly_add(g_env.get(term.alpha, {}), term.poly)

    rates = sorted(set(f_env) | set(g_env) | set(h_env))
    f = {a: f_env.get(a, {}) for a in rates}
    g = {a: g_env.get(a, {}) for a in rates}
    h = {a: h_env.get(a, {}) for a in rates}
    fs = {a: _closure_d_s(f[a], a) for a in rates}
    fz = {a: _closure_d_z(f[a], a) for a in rates}
    gs = {a: _closure_d_s(g[a], a) for a in rates}
    hs = {a: _closure_d_s(h[a], a) for a in rates}
    hz = {a: _closure_d_z(h[a], a) for a in rates}

    source: dict[Fraction, dict] = {}
    for a1 in rates:
        for a2 in rates:
            acc = _poly_scale(_poly_mul(f[a1], f[a2]), Fraction(2))
            acc = _poly_add(
                acc, _poly_shift(_poly_mul(f[a1], fs[a2]), 1, 0, Fraction(4))
            )
            acc = _poly_add(
                acc, _poly_shift(_poly_mul(fs[a1], fs[a2]), 2, 0, Fraction(4))
            )
            acc = _poly_add(acc, _poly_scale(_poly_mul(g[a1], g[a2]), Fraction(-2)))
            acc = _poly_add(
                acc, _poly_shift(_poly_mul(g[a1], gs[a2]), 1, 0, Fraction(-4))
            )
            acc = _poly_add(
                acc, _poly_shift(_poly_mul(hs[a1], fz[a2]), 1, 0, Fraction(4))
            )
            acc = _poly_add(acc, _poly_mul(hz[a1], hz[a2]))
            if acc:
                combined = a1 + a2
                source[combined] = _poly_add(source.get(combined, {}), acc)
    return {a: poly for a, poly in source.items() if poly}


# --------------------------------------------------------------------------- #
# the Laguerre / Hermite decomposition and the exact t-integral                #
# --------------------------------------------------------------------------- #


def _laguerre_coefficients(d: int) -> dict[int, Fraction]:
    r"""``\ell_d(s;v) = \sum_k a_{dk}v^{k-d}s^k`` with ``a_{dk} =
    (-1)^k\binom{d}{k}/k!`` — the ``v``-free table."""
    return {
        k: Fraction((-1) ** k * math.comb(d, k), math.factorial(k))
        for k in range(d + 1)
    }


def _hermite_table(j_max: int) -> list[dict[int, int]]:
    r"""``\hat h_j(z;v) = \sum_i b_{ji}v^{(i-j)/2}z^i`` — integer tables from
    ``H_{j+1} = 2xH_j - 2jH_{j-1}``, so ``b_{j+1,i} = 2b_{j,i-1} -
    2j\,b_{j-1,i}``."""
    table: list[dict[int, int]] = [{0: 1}]
    if j_max >= 1:
        table.append({1: 2})
    for j in range(1, j_max):
        nxt: dict[int, int] = {}
        for i, c in table[j].items():
            nxt[i + 1] = nxt.get(i + 1, 0) + 2 * c
        for i, c in table[j - 1].items():
            nxt[i] = nxt.get(i, 0) - 2 * j * c
        table.append({i: c for i, c in nxt.items() if c})
    return table


def _decompose_radial(m: int, rate: Fraction) -> dict[int, Fraction]:
    """``s^m = sum_d out[d] * ell_d(s; rate)`` by triangular elimination."""
    rep = {m: Fraction(1)}
    out: dict[int, Fraction] = {}
    for d in range(m, -1, -1):
        table = _laguerre_coefficients(d)
        c = rep.get(d, Fraction(0)) / table[d]
        if c:
            out[d] = c
            for k, a_dk in table.items():
                rep[k] = rep.get(k, Fraction(0)) - c * a_dk * rate ** (k - d)
        rep = {k: value for k, value in rep.items() if value}
    if rep:  # pragma: no cover - triangularity makes this unreachable
        raise AssertionError("radial decomposition failed to terminate")
    return out


def _decompose_axial(
    n: int, rate: Fraction, table: list[dict[int, int]]
) -> dict[int, Fraction]:
    """``z^n = sum_j out[j] * hhat_j(z; rate)`` by triangular elimination."""
    rep = {n: Fraction(1)}
    out: dict[int, Fraction] = {}
    for j in range(n, -1, -1):
        c = rep.get(j, Fraction(0)) / Fraction(table[j][j])
        if c:
            out[j] = c
            for i, b_ji in table[j].items():
                rep[i] = rep.get(i, Fraction(0)) - c * b_ji * rate ** ((i - j) // 2)
        rep = {i: value for i, value in rep.items() if value}
    if rep:  # pragma: no cover - triangularity makes this unreachable
        raise AssertionError("axial decomposition failed to terminate")
    return out


_Rep = dict  # {(k, i, mu2): Fraction} -- one envelope's g-monomial terms


def _potential_terms(poly: dict, rate: Fraction) -> _Rep:
    r"""``(-\Delta)^{-1}[P(s,z)e^{-A\rho^2}]`` as g-monomials, exactly.

    Per source monomial ``s^mz^n``: decompose into ``\ell_d\hat h_j``
    components, evolve each with ``\lambda^{-(2d+j+3/2)}`` and coefficient
    polynomials in ``v``, and integrate each ``v``-monomial with
    ``\int_0^\infty\lambda^{-\nu}e^{-A\rho^2/\lambda}dt =
    g_{\nu-2}(A\rho^2)/(4A)``.  The result is
    ``\sum c\,s^kz^i\,g_\mu(A\rho^2)`` with ``\mu = d+k+(j+i)/2-1/2`` stored
    as ``mu2 = 2\mu``.
    """
    n_max = max((n for (_, n) in poly), default=0)
    table = _hermite_table(max(n_max, 1))
    out: _Rep = {}
    for (m, n), c in poly.items():
        radial = _decompose_radial(m, rate)
        axial = _decompose_axial(n, rate, table)
        for d, s_d in radial.items():
            laguerre = _laguerre_coefficients(d)
            for j, z_j in axial.items():
                base = c * s_d * z_j
                for k, a_dk in laguerre.items():
                    for i, b_ji in table[j].items():
                        q = (k - d) + (i - j) // 2
                        coefficient = base * a_dk * b_ji * rate**q / (4 * rate)
                        mu2 = 2 * d + 2 * k + j + i - 1
                        key = (k, i, mu2)
                        total = out.get(key, Fraction(0)) + coefficient
                        if total == 0:
                            out.pop(key, None)
                        else:
                            out[key] = total
    return out


def _rep_d_s(rep: _Rep, rate: Fraction) -> _Rep:
    r"""``\partial_s`` on g-monomials: ``\partial_s[s^kz^ig_\mu] =
    ks^{k-1}z^ig_\mu - As^kz^ig_{\mu+1}`` (``g_\mu' = -g_{\mu+1}``)."""
    out: _Rep = {}

    def add(key, value):
        total = out.get(key, Fraction(0)) + value
        if total == 0:
            out.pop(key, None)
        else:
            out[key] = total

    for (k, i, mu2), c in rep.items():
        if k >= 1:
            add((k - 1, i, mu2), k * c)
        add((k, i, mu2 + 2), -rate * c)
    return out


def _rep_d_z(rep: _Rep, rate: Fraction) -> _Rep:
    out: _Rep = {}

    def add(key, value):
        total = out.get(key, Fraction(0)) + value
        if total == 0:
            out.pop(key, None)
        else:
            out[key] = total

    for (k, i, mu2), c in rep.items():
        if i >= 1:
            add((k, i - 1, mu2), i * c)
        add((k, i + 1, mu2 + 2), -2 * rate * c)
    return out


# --------------------------------------------------------------------------- #
# rigorous enclosures of g_mu                                                  #
# --------------------------------------------------------------------------- #

_G_CACHE: dict[tuple, tuple[int, dict[int, Interval]]] = {}


def clear_spectral_cache() -> None:
    """Empty the g-function cache.  Only tests should need this."""
    _G_CACHE.clear()


def _g_series(mu2: int, w: Fraction, bits: int) -> Interval:
    r"""``g_\mu(w) = \sum_k(-1)^kw^k/(k!(\mu+k+1))`` for one ``\mu``, exact.

    The remainder after the last included term is bounded by the first omitted
    term once the index exceeds ``w`` (the term ratio is below ``w/(k+1)``),
    so the loop runs at least past ``w`` and stops when the term is below the
    working tolerance.
    """
    tolerance = Fraction(1, 2 ** (bits + 8))
    power = Fraction(1)  # w^k / k!
    total = Fraction(0)
    sign = 1
    k = 0
    while True:
        term = power * Fraction(2, mu2 + 2 * k + 2)
        if k > w and term < tolerance:
            return Interval(total - term, total + term)
        total += sign * term
        sign = -sign
        power = power * w / (k + 1)
        k += 1
        if k > 100000:  # pragma: no cover - defensive
            raise RuntimeError("g-series failed to converge")


def _g_base_erf(w: Fraction, bits: int, terms: int) -> Interval:
    r"""``g_{-1/2}(w) = (\sqrt\pi - 2I(\sqrt w))/\sqrt w`` for large ``w``,
    with the two-sided elementary bound on ``I`` from the module docstring."""
    x = sqrt_interval(Interval(w, w), bits=bits)
    envelope = exp_interval(Interval(-w, -w), terms=terms)
    upper = divide(envelope, x.scale(Fraction(2)))
    lower = divide(x * envelope, Interval(2 * w + 1, 2 * w + 1))
    integral = Interval(
        max(lower.lower, Fraction(0)), upper.upper
    )
    sqrt_pi = sqrt_interval(_PI, bits=bits)
    return _round(divide(sqrt_pi - integral.scale(Fraction(2)), x), bits)


def _g_family(
    w: Fraction, mu2_max: int, *, bits: int, terms: int
) -> dict[int, Interval]:
    r"""Enclosures of ``g_\mu(w)`` for every odd ``mu2 \le mu2\_max``, cached.

    Small ``w``: series at the top ``\mu``, then the downward recursion
    ``g_{\mu-1} = (w\,g_\mu + e^{-w})/\mu`` — positive quantities only.  Large
    ``w``: erf-seeded base, then the upward recursion ``g_\mu =
    (\mu g_{\mu-1} - e^{-w})/w`` — the subtracted term is exponentially
    smaller than the retained one.  Every value is intersected with the a
    priori range ``0 < g_\mu \le 1/(\mu+1)``.
    """
    key = (w, bits, terms)
    cached = _G_CACHE.get(key)
    if cached is not None and cached[0] >= mu2_max:
        return cached[1]
    envelope = exp_interval(Interval(-w, -w), terms=terms)
    values: dict[int, Interval] = {}

    def clamp(mu2: int, value: Interval) -> Interval:
        upper = Fraction(2, mu2 + 2)
        return _round(
            Interval(
                max(value.lower, Fraction(0)), min(value.upper, upper)
            ),
            bits,
        )

    if w < _G_SERIES_SWITCH:
        top = _g_series(mu2_max, w, bits)
        values[mu2_max] = clamp(mu2_max, top)
        for mu2 in range(mu2_max, -1, -2):
            if mu2 - 2 < -1:
                break
            mu = Fraction(mu2, 2)
            below = divide(
                (Interval(w, w) * values[mu2]) + envelope,
                Interval(mu, mu),
            )
            values[mu2 - 2] = clamp(mu2 - 2, below)
    else:
        values[-1] = clamp(-1, _g_base_erf(w, bits, terms))
        w_box = Interval(w, w)
        for mu2 in range(1, mu2_max + 1, 2):
            mu = Fraction(mu2, 2)
            above = divide(
                values[mu2 - 2].scale(mu) - envelope, w_box
            )
            values[mu2] = clamp(mu2, above)
    _G_CACHE[key] = (mu2_max, values)
    return values


# --------------------------------------------------------------------------- #
# the spectral pressure object                                                 #
# --------------------------------------------------------------------------- #


#: Highest derivative order kept in the symbolic lattice.  Order ``N`` Taylor
#: models need the lattice up to ``N + 1``.
_LATTICE_ORDER = 6


@dataclass(frozen=True)
class SpectralPressure:
    r"""The closed-form free-space pressure of one isotropic family.

    ``source`` is ``\sigma`` per envelope; ``derivatives[(a, b)]`` is the
    g-monomial representation of ``\partial_s^a\partial_z^b\,p`` per envelope
    (``\partial_rp = 2r\,\partial_sp``), kept for every ``a + b \le 6``.  All
    data are exact rationals; intervals appear only at evaluation time.
    """

    source: dict[Fraction, dict]
    derivatives: dict[tuple[int, int], dict[Fraction, _Rep]]

    @property
    def p(self) -> dict[Fraction, _Rep]:
        return self.derivatives[(0, 0)]

    @property
    def ps(self) -> dict[Fraction, _Rep]:
        return self.derivatives[(1, 0)]

    @property
    def pz(self) -> dict[Fraction, _Rep]:
        return self.derivatives[(0, 1)]

    @classmethod
    def from_family(cls, family) -> "SpectralPressure":
        source = isotropic_pressure_source(family)
        if not source:
            raise ValueError("the family has an identically zero pressure source")
        lattice: dict[tuple[int, int], dict[Fraction, _Rep]] = {
            (0, 0): {a: _potential_terms(poly, a) for a, poly in source.items()}
        }
        for order in range(1, _LATTICE_ORDER + 1):
            for a in range(order + 1):
                b = order - a
                if a > 0:
                    lattice[(a, b)] = {
                        rate: _rep_d_s(rep, rate)
                        for rate, rep in lattice[(a - 1, b)].items()
                    }
                else:
                    lattice[(a, b)] = {
                        rate: _rep_d_z(rep, rate)
                        for rate, rep in lattice[(a, b - 1)].items()
                    }
        return cls(source=source, derivatives=lattice)

    # -- exact self-verification -------------------------------------------- #

    def verify_laplace_identity(self) -> dict[str, object]:
        r"""Reduce ``\Delta_3p + \sigma`` to zero, exactly.

        ``\Delta_3 = 4s\partial_{ss} + 4\partial_s + \partial_{zz}`` on
        axisymmetric functions of ``(s, z)``.  The reduction repeatedly splits
        the top-``\mu`` coefficient polynomial as ``Q = (s+z^2)Q_1 + R`` (the
        remainder is ``Q`` at ``s = -z^2``) and rewrites ``(s+z^2)Q_1g_\mu``
        through ``w\,g_\mu = \mu g_{\mu-1} - e^{-w}``; the exponential parts
        must reassemble ``-\sigma`` and every division must be exact, because
        ``\{g_\mu\}\cup\{e^{-w}\}`` are linearly independent over polynomials
        (they have distinct decay laws in ``w``).  Any nonzero residual is
        reported and means the representation is wrong.
        """
        failures: list[str] = []
        for rate, rep in self.p.items():
            lap: _Rep = {}

            def add(target: _Rep, key, value):
                total = target.get(key, Fraction(0)) + value
                if total == 0:
                    target.pop(key, None)
                else:
                    target[key] = total

            for (k, i, mu2), c in self.derivatives[(2, 0)][rate].items():
                add(lap, (k + 1, i, mu2), 4 * c)
            for (k, i, mu2), c in self.derivatives[(1, 0)][rate].items():
                add(lap, (k, i, mu2), 4 * c)
            for (k, i, mu2), c in self.derivatives[(0, 2)][rate].items():
                add(lap, (k, i, mu2), c)

            exponential: dict = {}
            # Reduce level by level from the top.  The level set is re-derived
            # each pass because the reduction itself pushes coefficients onto
            # mu2 - 2, which may not have been populated initially; a static
            # level list would silently skip those and the check would not be
            # a check.
            while lap:
                mu2 = max(m2 for (_, _, m2) in lap)
                block = {
                    (k, i): c for (k, i, m2), c in lap.items() if m2 == mu2
                }
                for key in [key for key in lap if key[2] == mu2]:
                    del lap[key]
                if mu2 == -1:
                    if block:
                        failures.append(
                            f"envelope {rate}: residual g_(-1/2) coefficient "
                            f"{block}"
                        )
                    continue
                quotient, remainder = _divide_by_w(block)
                if remainder:
                    failures.append(
                        f"envelope {rate}: mu2={mu2} coefficient not divisible "
                        f"by (s+z^2): remainder {remainder}"
                    )
                    continue
                mu = Fraction(mu2, 2)
                for (k, i), c in quotient.items():
                    add(lap, (k, i, mu2 - 2), c * mu / rate)
                    add(exponential, (k, i), -c / rate)
            sigma = self.source[rate]
            residual = _poly_add(exponential, sigma)
            if residual:
                failures.append(
                    f"envelope {rate}: exponential part does not cancel the "
                    f"source: residual {residual}"
                )
        return {"verified": not failures, "failures": failures}

    # -- evaluation --------------------------------------------------------- #

    def _rep_box(
        self,
        rep_group: dict[Fraction, _Rep],
        s_box: Interval,
        z_box: Interval,
        *,
        bits: int,
        terms: int,
    ) -> Interval:
        """Interval enclosure of one representation over an (s, z) box.

        The inner loop works on raw rational endpoints: the monomial box is
        cached per ``(k, i)``, the g-box is nonnegative by construction, and
        the sign cases of the product are enumerated directly.  This is
        outward-exact — the endpoints computed are exactly the endpoints the
        generic interval product would produce — but several times faster,
        which is what makes the Taylor-model sweep affordable.
        """
        z_sq = _z_power(z_box, 2)
        zero = Fraction(0)
        total_low = zero
        total_high = zero
        monomials: dict[tuple[int, int], tuple[Fraction, Fraction]] = {}
        for rate, rep in rep_group.items():
            if not rep:
                continue
            mu2_max = max(mu2 for (_, _, mu2) in rep)
            w_low = _floor_dyadic(
                max(rate * (s_box.lower + z_sq.lower), zero), bits
            )
            w_high = _ceil_dyadic(rate * (s_box.upper + z_sq.upper), bits)
            g_at_low = _g_family(w_low, mu2_max, bits=bits, terms=terms)
            g_at_high = _g_family(w_high, mu2_max, bits=bits, terms=terms)
            for (k, i, mu2), c in rep.items():
                cached = monomials.get((k, i))
                if cached is None:
                    box = _s_power(s_box, k) * _z_power(z_box, i)
                    # Outward dyadic rounding keeps every later product and the
                    # running sums on denominators of bounded size, which is
                    # what keeps the gcd normalisation inside Fraction cheap.
                    cached = (
                        _floor_dyadic(box.lower, bits),
                        _ceil_dyadic(box.upper, bits),
                    )
                    monomials[(k, i)] = cached
                m_low, m_high = cached
                g_low = g_at_high[mu2].lower   # g decreasing in w
                g_high = g_at_low[mu2].upper
                if m_low >= zero:
                    low, high = m_low * g_low, m_high * g_high
                elif m_high <= zero:
                    low, high = m_low * g_high, m_high * g_low
                else:
                    low, high = m_low * g_high, m_high * g_high
                if c >= zero:
                    total_low += c * low
                    total_high += c * high
                else:
                    total_low += c * high
                    total_high += c * low
        return _round(Interval(total_low, total_high), bits)

    def gradient_boxes(
        self,
        r_box: Interval,
        z_box: Interval,
        *,
        path: str = "taylor",
        taylor_order: int = 4,
        bits: int = DEFAULT_PRECISION_BITS,
        terms: int = 28,
        crude_boxes: dict[tuple[int, int], Fraction] | None = None,
    ) -> tuple[Interval, Interval]:
        r"""Enclose ``\partial_sp`` and ``\partial_zp`` over one cell.

        ``plain``: direct monomial-times-``g`` boxes.  The triangle inequality
        across the g-monomials loses their cancellation, which costs two to
        three orders of magnitude of width (measured ~500--1000x on the
        preregistered candidate), so the default is a **Taylor model of order**
        ``taylor_order``: the derivatives of ``\partial p`` up to order
        ``N - 1`` are evaluated at the cell centre as *thin* intervals — exact
        cancellation, the whole point of the closed form — the polynomial part
        is evaluated in interval arithmetic over the centred cell, and only the
        order-``N`` Lagrange remainder pays the triangle-inequality toll, now
        multiplied by ``h^N``.  Each order contracts the width by roughly
        ``4\times`` at the standard grids; the result is intersected with the
        plain box, so it is never worse.  ``\partial_rp = 2r\,\partial_sp`` is
        assembled by the caller so the ``r``-factor pairs sharply with the
        measure.

        ``crude_boxes`` optionally supplies precomputed magnitudes of the
        order-``N+1`` derivatives valid on a superset of the cell (a coarse
        block): a bound over a superset is a bound over the cell, so this is
        sound, and it amortises the widest term loops across many cells.
        """
        if r_box.lower < 0:
            raise ValueError("the radial box cannot extend below the axis")
        order = taylor_order
        if order < 1 or order + 1 > _LATTICE_ORDER:
            raise ValueError(
                f"taylor_order must be between 1 and {_LATTICE_ORDER - 1}"
            )
        s_box = _round(square(r_box), bits)
        p_s = self._rep_box(self.ps, s_box, z_box, bits=bits, terms=terms)
        p_z = self._rep_box(self.pz, s_box, z_box, bits=bits, terms=terms)
        if path == "plain":
            return p_s, p_z

        centre_s = (s_box.lower + s_box.upper) / 2
        centre_z = (z_box.lower + z_box.upper) / 2
        s_thin = Interval(centre_s, centre_s)
        z_thin = Interval(centre_z, centre_z)
        half_s = (s_box.upper - s_box.lower) / 2
        half_z = (z_box.upper - z_box.lower) / 2

        # Thin derivative values at the centre, orders 1 .. N of p.
        thin: dict[tuple[int, int], Interval] = {}
        for total in range(1, order + 1):
            for a in range(total + 1):
                b = total - a
                thin[(a, b)] = self._rep_box(
                    self.derivatives[(a, b)], s_thin, z_thin,
                    bits=bits, terms=terms,
                )
        # Crude boxes of the order-(N+1) derivatives of p over the cell (or a
        # superset supplied by the caller).
        if crude_boxes is None:
            crude: dict[tuple[int, int], Fraction] = {}
            for a in range(order + 2):
                b = order + 1 - a
                box = self._rep_box(
                    self.derivatives[(a, b)], s_box, z_box,
                    bits=bits, terms=terms,
                )
                crude[(a, b)] = box.magnitude
        else:
            crude = crude_boxes

        def sym_power(half: Fraction, exponent: int) -> Interval:
            # Sharp range of Delta^exponent over [-half, half], by parity.
            if exponent == 0:
                return Interval(Fraction(1), Fraction(1))
            top = half**exponent
            if exponent % 2 == 0:
                return Interval(Fraction(0), top)
            return Interval(-top, top)

        def model(shift_s: int, shift_z: int) -> Interval:
            # Taylor model of d_s^{shift_s} d_z^{shift_z} p over the cell.
            total_poly = thin[(shift_s, shift_z)]
            for total in range(1, order):
                for a in range(total + 1):
                    b = total - a
                    coefficient = thin[(shift_s + a, shift_z + b)].scale(
                        Fraction(1, math.factorial(a) * math.factorial(b))
                    )
                    term = sym_power(half_s, a) * sym_power(half_z, b)
                    total_poly = total_poly + _round(coefficient * term, bits)
            remainder = Fraction(0)
            for a in range(order + 1):
                b = order - a
                remainder += (
                    crude[(shift_s + a, shift_z + b)]
                    * half_s**a * half_z**b
                    / (math.factorial(a) * math.factorial(b))
                )
            remainder = _ceil_dyadic(remainder, bits)
            return _round(total_poly + Interval(-remainder, remainder), bits)

        p_s = _intersect(p_s, model(1, 0))
        p_z = _intersect(p_z, model(0, 1))
        return p_s, p_z

    def gradient_value(
        self,
        r: Fraction | float,
        z: Fraction | float,
        *,
        bits: int = DEFAULT_PRECISION_BITS,
        terms: int = 28,
    ) -> tuple[Interval, Interval]:
        r"""Thin enclosures of ``(\partial_rp, \partial_zp)`` at one point."""
        r = Fraction(r)
        z = Fraction(z)
        p_s, p_z = self.gradient_boxes(
            Interval(r, r), Interval(z, z), path="plain", bits=bits, terms=terms
        )
        return p_s.scale(2 * r), p_z

    def pressure_value(
        self,
        r: Fraction | float,
        z: Fraction | float,
        *,
        bits: int = DEFAULT_PRECISION_BITS,
        terms: int = 28,
    ) -> Interval:
        r = Fraction(r)
        z = Fraction(z)
        return self._rep_box(
            self.p,
            Interval(r * r, r * r),
            Interval(z, z),
            bits=bits,
            terms=terms,
        )

    # -- the exterior bound -------------------------------------------------- #

    def exterior_gradient_bound(
        self, distance: Fraction, *, bits: int = DEFAULT_PRECISION_BITS
    ) -> Fraction:
        r"""A rational upper bound on ``\sup_{\rho\ge d}(|\partial_rp| +
        |\partial_zp|)``.

        Each g-monomial obeys ``g_\mu(w) \le \Gamma(\mu+1)/w^{\mu+1}``
        (``\gamma \le \Gamma``), so ``\rho^P g_\mu(A\rho^2) \le
        \Gamma(\mu+1)A^{-(\mu+1)}\rho^{P-2\mu-2}`` with ``P`` the total power
        of ``\rho`` dominating the monomial (``s^kz^i \le \rho^{2k+i}``, one
        more for the explicit ``r`` of ``\partial_rp``).  The construction
        guarantees ``P \le 2\mu+2`` for every term — asserted here — so the
        bound is decreasing in ``\rho`` and its value at ``\rho = d`` bounds
        the whole exterior.  ``\Gamma`` at half-integers is
        ``\sqrt\pi\,(2q)!/(4^qq!)``.
        """
        if distance <= 0:
            raise ValueError("the exterior bound needs a positive distance")
        sqrt_pi_upper = sqrt_interval(_PI, bits=bits).upper
        total = Fraction(0)
        for rep_group, extra in ((self.ps, 1), (self.pz, 0)):
            for rate, rep in rep_group.items():
                sqrt_rate_lower = sqrt_interval(
                    Interval(rate, rate), bits=bits
                ).lower
                for (k, i, mu2), c in rep.items():
                    power = 2 * k + i + extra
                    exponent = power - (mu2 + 2)
                    if exponent > 0:  # pragma: no cover - construction bound
                        raise AssertionError(
                            "a g-monomial grows relative to its decay; the "
                            "exterior bound derivation no longer applies"
                        )
                    q = (mu2 + 1) // 2
                    gamma_upper = (
                        sqrt_pi_upper
                        * Fraction(math.factorial(2 * q),
                                   4**q * math.factorial(q))
                    )
                    rate_power = rate**q * sqrt_rate_lower  # lower bound of A^(mu+1)
                    scale = Fraction(2) if extra else Fraction(1)
                    total += (
                        scale * abs(c) * gamma_upper
                        / rate_power
                        * distance**exponent
                    )
        return Interval(Fraction(0), total).round_outward(bits).upper


def _divide_by_w(block: dict) -> tuple[dict, dict]:
    r"""Divide a polynomial ``{(k, i): c}`` by ``(s + z^2)``, exactly.

    Descending synthetic division on the identity ``s^kz^i = (s+z^2)\,
    s^{k-1}z^i - s^{k-1}z^{i+2}``: each level-``k`` term moves one factor into
    the quotient and pushes an exact compensation to level ``k-1``, so after
    the sweep only ``s^0`` terms remain — the remainder, which is the
    polynomial at ``s = -z^2``.
    """
    work = {key: c for key, c in block.items() if c}
    quotient: dict = {}
    k_max = max((k for (k, _) in work), default=0)
    for k in range(k_max, 0, -1):
        for key in [key for key in work if key[0] == k]:
            c = work.pop(key)
            if not c:
                continue
            _, i = key
            q_key = (k - 1, i)
            quotient[q_key] = quotient.get(q_key, Fraction(0)) + c
            push = (k - 1, i + 2)
            work[push] = work.get(push, Fraction(0)) - c
    remainder = {key: c for key, c in work.items() if c}
    quotient = {key: c for key, c in quotient.items() if c}
    return quotient, remainder


def _intersect(first: Interval, second: Interval) -> Interval:
    lower = max(first.lower, second.lower)
    upper = min(first.upper, second.upper)
    if lower > upper:  # both enclose the truth; be safe
        return first
    return Interval(lower, upper)


# --------------------------------------------------------------------------- #
# the certificate                                                              #
# --------------------------------------------------------------------------- #


def build_spectral_generation_certificate(
    family,
    grid: AxisymmetricGrid,
    *,
    viscosity: float,
    path: str = "taylor",
    taylor_order: int = 4,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    exp_terms: int = 28,
    max_depth: int = 2,
    refine_budget: int = 8192,
    self_check: bool = True,
) -> GenerationCertificate:
    r"""Certify a lower bound on ``J(u_0)`` with the closed-form pressure.

    The pairing integrand is ``-3|u|(u^r\partial_rp + u^z\partial_zp)`` with
    ``\partial p`` the *exact* free-space pressure gradient of the datum —
    enclosed per cell through the g-monomial representation — so the payload
    carries **no pressure hypothesis**: the ``hypotheses`` block is empty and
    the verifier enforces it.  The exterior pressure remainder is enclosed
    too (``exterior_pressure_remainder``), so ``margin`` is an unconditional
    lower bound on ``J(u_0)``.
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
    if taylor_order < 1 or taylor_order + 1 > _LATTICE_ORDER:
        raise ValueError(
            f"taylor_order must be between 1 and {_LATTICE_ORDER - 1}"
        )
    if max_depth < 0 or refine_budget < 0:
        raise ValueError("max_depth and refine_budget must be nonnegative")
    if not (float(grid.z[0]) < 0.0 < float(grid.z[-1])):
        raise ValueError(
            "the exterior bounds need a box that surrounds the origin: "
            "z_min < 0 < z_max"
        )

    spectral = SpectralPressure.from_family(family)
    if self_check:
        report = spectral.verify_laplace_identity()
        if not report["verified"]:
            raise AssertionError(
                f"the spectral representation failed -Delta p = sigma: "
                f"{report['failures']}"
            )

    solver = AxisymmetricPressureSolver.build(grid)
    from .interval_gaussian import gaussian_reference_field

    field = gaussian_reference_field(family, grid)
    rate = l3_generation_rate(field, viscosity=viscosity, solver=solver)

    r_nodes = [Fraction(float(v)) for v in grid.r]
    z_nodes = [Fraction(float(v)) for v in grid.z]
    two_pi_upper = _TWO_PI.upper

    pressure_term = _ZERO
    viscous_upper = Fraction(0)
    natural_width = Fraction(0)
    width_from_pressure = Fraction(0)
    width_from_velocity = Fraction(0)
    gradient_width_max = Fraction(0)
    cells = 0
    leaves = 0
    refined_cells = 0
    budget = refine_budget
    box_path = "plain" if path == "plain" else "taylor"

    #: The crude Taylor remainders are magnitudes of order-(N+1) derivatives;
    #: bounding them over 4x4 blocks of cells is sound (a superset bound) and
    #: amortises the widest term loops.
    block = 4
    block_crude: dict[tuple[int, int], dict[tuple[int, int], Fraction]] = {}

    def crude_for(cell_i: int, cell_j: int):
        if box_path == "plain":
            return None
        key = (cell_i // block, cell_j // block)
        cached = block_crude.get(key)
        if cached is None:
            r_box = Interval(
                r_nodes[key[0] * block],
                r_nodes[min((key[0] + 1) * block, grid.nr - 1)],
            )
            z_box = Interval(
                z_nodes[key[1] * block],
                z_nodes[min((key[1] + 1) * block, grid.nz - 1)],
            )
            s_box = _round(square(r_box), precision_bits)
            cached = {}
            for a in range(taylor_order + 2):
                b = taylor_order + 1 - a
                box = spectral._rep_box(
                    spectral.derivatives[(a, b)], s_box, z_box,
                    bits=precision_bits, terms=exp_terms,
                )
                cached[(a, b)] = box.magnitude
            block_crude[key] = cached
        return cached

    def leaf_terms(r_lo, r_hi, z_lo, z_hi, crude):
        r_box = Interval(r_lo, r_hi)
        z_box = Interval(z_lo, z_hi)
        velocity = cell_velocity_enclosure(
            family, r_box, z_box, terms=exp_terms, bits=precision_bits
        )
        p_s, p_z = spectral.gradient_boxes(
            r_box, z_box, path=box_path, taylor_order=taylor_order,
            bits=precision_bits, terms=exp_terms, crude_boxes=crude,
        )
        s_box = _round(square(r_box), precision_bits)
        area = (r_hi - r_lo) * (z_hi - z_lo)
        flow_r = _round(
            s_box.scale(Fraction(2)) * (velocity.speed * velocity.u_r),
            precision_bits,
        )
        flow_z = _round(
            r_box * (velocity.speed * velocity.u_z), precision_bits
        )
        paired = flow_r * p_s + flow_z * p_z
        natural = _round((_TWO_PI * paired).scale(area), precision_bits)
        weight_up = two_pi_upper * r_hi * area
        return velocity, p_s, p_z, flow_r, flow_z, natural, weight_up

    for i in range(grid.nr - 1):
        for j in range(grid.nz - 1):
            cells += 1
            cell_crude = crude_for(i, j)
            stack = [(r_nodes[i], r_nodes[i + 1], z_nodes[j], z_nodes[j + 1], 0)]
            cell_was_refined = False
            while stack:
                r_lo, r_hi, z_lo, z_hi, depth = stack.pop()
                velocity, p_s, p_z, flow_r, flow_z, natural, weight_up = (
                    leaf_terms(r_lo, r_hi, z_lo, z_hi, cell_crude)
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
                pressure_term = pressure_term + natural.scale(Fraction(-3))
                natural_width += 3 * (natural.upper - natural.lower)
                area = (r_hi - r_lo) * (z_hi - z_lo)
                scale = 3 * two_pi_upper * area
                width_from_pressure += scale * (
                    flow_r.magnitude * (p_s.upper - p_s.lower)
                    + flow_z.magnitude * (p_z.upper - p_z.lower)
                )
                width_from_velocity += scale * (
                    p_s.magnitude * (flow_r.upper - flow_r.lower)
                    + p_z.magnitude * (flow_z.upper - flow_z.lower)
                )
                gradient_width_max = max(
                    gradient_width_max,
                    p_s.upper - p_s.lower,
                    p_z.upper - p_z.lower,
                )
                viscous_upper += weight_up * max(
                    velocity.viscous_integrand.upper, Fraction(0)
                )
            if cell_was_refined:
                refined_cells += 1
        pressure_term = pressure_term.round_outward(precision_bits)
        viscous_upper = Interval(viscous_upper, viscous_upper).round_outward(
            precision_bits
        ).upper
        natural_width = Interval(natural_width, natural_width).round_outward(
            precision_bits
        ).upper
        width_from_pressure = Interval(
            width_from_pressure, width_from_pressure
        ).round_outward(precision_bits).upper
        width_from_velocity = Interval(
            width_from_velocity, width_from_velocity
        ).round_outward(precision_bits).upper

    # -- the exterior contributions ----------------------------------------- #
    z_extent = min(-z_nodes[0], z_nodes[-1])
    tail = exterior_viscous_tail(
        family, r_nodes[-1], z_extent, terms=exp_terms, bits=precision_bits
    )
    distance = min(r_nodes[-1], z_extent)
    gradient_sup = spectral.exterior_gradient_bound(
        distance, bits=precision_bits
    )
    four_pi_upper = _TWO_PI.upper * 2
    speed_mass = (
        tail.speed_coefficient**2
        * four_pi_upper
        * gaussian_tail_moment(
            tail.gamma, distance, 2, terms=exp_terms, bits=precision_bits
        )
    )
    exterior_pressure = Interval(
        Fraction(0), 3 * gradient_sup * speed_mass
    ).round_outward(precision_bits).upper

    interior_viscous = viscous_upper
    viscous_total = viscous_upper + tail.viscous_tail
    nu = Fraction(float(viscosity))
    viscous_bound = -3 * nu * viscous_total
    pressure_lower = pressure_term.lower
    margin = pressure_lower - exterior_pressure + viscous_bound

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate": SPECTRAL_CERTIFICATE_KIND,
        "family": family.as_dict(),
        "grid": {
            "nr": grid.nr, "nz": grid.nz,
            "r_max": float(grid.r[-1]),
            "z_min": float(grid.z[0]), "z_max": float(grid.z[-1]),
        },
        "viscosity": float(viscosity),
        "path": path,
        "taylor_order": taylor_order,
        "rigorous": {
            "enclosure_basis": "gaussian_hermite",
            "pressure_pairing": "spectral_closed_form",
            "pressure_term_lower_bound": str(pressure_lower),
            "pressure_term_upper_bound": str(pressure_term.upper),
            "pressure_term_width": str(pressure_term.upper - pressure_term.lower),
            "exterior_pressure_remainder": str(exterior_pressure),
            "exterior_gradient_sup_bound": str(gradient_sup),
            "viscous_lower_bound": str(viscous_bound),
            "viscous_integral_upper_bound": str(viscous_total),
            "interior_viscous_integral_upper_bound": str(interior_viscous),
            "exterior_viscous_tail": str(tail.viscous_tail),
            "margin": str(margin),
            "margin_is_positive": bool(margin > 0),
            "cells_enclosed": cells,
            "leaf_evaluations": leaves,
            "refined_cells": refined_cells,
            "refine_budget_left": budget,
            "natural_product_width": str(natural_width),
            "width_from_pressure_gradient": str(width_from_pressure),
            "width_from_velocity_factor": str(width_from_velocity),
            "pressure_gradient_box_width_max": str(gradient_width_max),
            "tail_gamma": str(tail.gamma),
            "tail_distance": str(tail.distance),
            "laplace_identity_checked": bool(self_check),
        },
        "hypotheses": {},
        "floating_point_reference": rate.as_dict(),
        "diagnostics": {
            "pressure_term_width": float(
                pressure_term.upper - pressure_term.lower
            ),
            "width_from_pressure_gradient": float(width_from_pressure),
            "width_from_velocity_factor": float(width_from_velocity),
            "exterior_pressure_remainder": float(exterior_pressure),
            "pressure_gradient_box_width_max": float(gradient_width_max),
            "float_pressure_reference": float(rate.pressure),
            "float_pressure_in_enclosure": bool(
                pressure_term.lower
                <= Fraction(float(rate.pressure))
                <= pressure_term.upper
            ),
        },
        "claims": [
            "Rigorous and unconditional: the stated enclosure of the box "
            "integral -3 int |u| u . grad p dx for the TRUE free-space "
            "pressure of the datum (a finite closed form in Laguerre/Hermite "
            "heat components and the entire functions g_mu, verified against "
            "-Laplacian p = sigma by exact rational reduction), the bound on "
            "the exterior pressure remainder, and the viscous upper bound "
            "with its Mills-ratio exterior tail.",
            "No pressure hypothesis: this certificate replaces the discrete "
            "Poisson solve of the P1/P1G chain by the analytic free-space "
            "pressure, so the hypotheses block is empty and margin is an "
            "unconditional lower bound on J(u_0).",
            "The advective contribution is exactly zero by an algebraic "
            "identity, not by cancellation of computed numbers.",
            "Nothing here is a singularity result and nothing here bears on "
            "the Clay problem.  A positive initial generation rate of the "
            "critical norm is a screening criterion, not evidence of blow-up.",
        ],
    }
    return GenerationCertificate(payload)


# --------------------------------------------------------------------------- #
# the independent checker                                                      #
# --------------------------------------------------------------------------- #


def verify_spectral_generation_certificate(
    payload: dict[str, object],
) -> dict[str, object]:
    """Independently re-check a spectral-pressure payload from its contents."""
    try:
        return _verify(payload)
    except Exception as error:  # noqa: BLE001 - adversarial input is a failure
        return {"verified": False, "failures": [f"payload rejected: {error}"]}


def _verify(payload: dict[str, object]) -> dict[str, object]:
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(
        payload.get("certificate") == SPECTRAL_CERTIFICATE_KIND,
        "wrong certificate kind",
    )
    require(payload.get("schema_version") == 1, "unsupported schema version")
    require(float(payload["viscosity"]) > 0.0, "viscosity must be positive")
    require(payload.get("path") in _PATHS, "unknown pairing path")
    require(
        payload.get("hypotheses") == {},
        "the spectral certificate must carry no hypotheses",
    )

    block = payload.get("rigorous")
    require(isinstance(block, dict), "missing rigorous block")
    if isinstance(block, dict):
        pressure_lower = Fraction(str(block["pressure_term_lower_bound"]))
        pressure_upper = Fraction(str(block["pressure_term_upper_bound"]))
        exterior = Fraction(str(block["exterior_pressure_remainder"]))
        viscous = Fraction(str(block["viscous_lower_bound"]))
        margin = Fraction(str(block["margin"]))
        require(
            pressure_lower <= pressure_upper,
            "the pressure enclosure endpoints are out of order",
        )
        require(
            Fraction(str(block["pressure_term_width"]))
            == pressure_upper - pressure_lower,
            "the pressure width does not match its endpoints",
        )
        require(exterior >= 0, "the exterior remainder bound must be nonnegative")
        require(viscous <= 0, "the viscous contribution must be non-positive")
        require(
            margin == pressure_lower - exterior + viscous,
            "the margin does not equal P_lower - exterior + V_lower",
        )
        require(
            bool(block["margin_is_positive"]) == (margin > 0),
            "the margin flag contradicts the margin",
        )
        total = Fraction(str(block["viscous_integral_upper_bound"]))
        interior = Fraction(str(block["interior_viscous_integral_upper_bound"]))
        tail = Fraction(str(block["exterior_viscous_tail"]))
        require(
            total == interior + tail,
            "interior viscous + exterior tail must equal the total",
        )
        require(tail > 0, "a Gaussian datum must carry a positive exterior tail")
        require(
            block.get("laplace_identity_checked") is True,
            "the exact -Delta p = sigma self-check must have run",
        )

    claims = payload.get("claims")
    require(
        isinstance(claims, list) and any("Clay" in str(c) for c in claims),
        "the payload must record that it says nothing about the Clay problem",
    )
    require(
        isinstance(claims, list) and any("unconditional" in str(c) for c in claims),
        "the payload must state which bounds are unconditional",
    )
    return {"verified": not failures, "failures": failures}
