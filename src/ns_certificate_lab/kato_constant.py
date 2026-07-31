r"""The ``n = 3`` Kato-constant certificate (Lane K, Track P).

This module turns the derivation of
``docs/research_notes/kato_h3_constants.md`` into a machine-checkable
payload: a rigorous rational upper bound ``G3 <= 12 sqrt(A_4)`` for the
trilinear form

.. math::  |\langle P(v\cdot\nabla w),\, w\rangle_{\dot H^3}|
           \;\le\; G_3\, \|v\|_{\dot H^3}\, \|w\|_{\dot H^3}^2

on mean-zero divergence-free real fields on ``T^3``, in the repository's own
normalisation (section 1 of the note, reproduced verbatim in
:data:`CONVENTIONS` and embedded in every payload), together with

* the basic-inequality pair ``K_3``: coefficients ``4 sqrt(A_4)`` on
  ``||v||_3 ||w||_3`` and ``4 sqrt(A_6)`` on ``||v||_3 ||w||_4`` (section 7);
* the per-trajectory band constants :func:`c_kato_upper` (section 5a) and
  :func:`c_shift_upper` (section 5b) that the chain uses as the *linear*
  coefficient of the ``n = 3`` control inequality — exact finite band sums
  with no universal constant at all;
* a monotonicity witness (the bound at ``2 N`` is provably ``<=`` the bound
  at ``N``) and the four named crude steps recorded as sharpening levers;
* a literature comparison (Morosi–Pizzocchero) that is **verification only**
  — no external decimal is an input to any bound here, and the checker
  rejects a payload in which that role has been altered.

Everything is exact ``Fraction`` arithmetic; the only irrational operations
are outward-rounded square roots via
:func:`~ns_certificate_lab.l3_certificate.sqrt_interval`.  The checker
:func:`verify_kato_certificate` recomputes the lattice partial sums with its
own loop (an octant loop with symmetry multiplicities — deliberately a
*different* traversal from the builder's full-cube loop) and re-derives every
relation the payload asserts.

Scope: this is a constants certificate and nothing else.  No statement about
Navier–Stokes regularity, blow-up, or any Clay problem is made or implied;
the infinite-dimensional commutator inequality behind ``G3`` lives in the
audited note (section 9), is never axiomatised in Lean, and no payload field
claims otherwise.
"""

from __future__ import annotations

from fractions import Fraction
import sys

from .fourier_torus import TrigVector, _norm_sq, adot_squared_upper, euclidean_norm_upper
from .l3_certificate import sqrt_interval
from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

__all__ = [
    "CLAIMS",
    "CONVENTIONS",
    "DERIVATION_NOTE",
    "SCHEMA",
    "a4_upper",
    "a6_upper",
    "build_kato_certificate",
    "c_kato_upper",
    "c_shift_upper",
    "g3_upper",
    "verify_kato_certificate",
]

SCHEMA = "kato_h3_constant"

DERIVATION_NOTE = "docs/research_notes/kato_h3_constants.md"

#: Section 1 of the derivation note, verbatim (each table row as
#: ``object | convention``).  The payload embeds the full block and the
#: checker rejects any alteration: a constant is meaningless without the
#: normalisation it was derived in.
CONVENTIONS = [
    "torus | `T^3 = (R / 2 pi Z)^3`, coordinates `x`",
    "measure | normalised: `d mu = (2 pi)^{-3} dx`; `mu(T^3) = 1`",
    "Fourier basis | `e^{i k . x}`, `k` in `Z^3`; real fields: `u = sum_k uhat_k e^{ikx}`, `uhat_{-k} = conj(uhat_k)`",
    "repo real coefficients | canonical `k`: `u = sum a_k cos(kx) + b_k sin(kx)`, so `uhat_k = (a_k - i b_k)/2`",
    "Parseval | `<f, g>_{L^2(mu)} = sum_k fhat_k . conj(ghat_k)` (no `2 pi` factors)",
    "`H^n dot` norm | `||u||_n^2 = sum_{k != 0} |k|^{2n} |uhat_k|^2` `=` repo `sobolev_sq(n)` `=` `sum_{canonical} |k|^{2n} (|a_k|^2+|b_k|^2)/2`",
    "zero mode | absent (`uhat_0 = 0`); all sums skip `k = 0`",
    "reality | `uhat_{-k} = conj(uhat_k)`, equivalent to the `cos`/`sin` representation",
    "divergence-free | `k . uhat_k = 0` for every `k`",
    "Leray multiplier | `P_k = I - k k^T / |k|^2`, orthogonal projection per mode, self-adjoint, `P w = w` for divergence-free `w`",
    "convolution | `(fg)hat_k = sum_j fhat_j ghat_{k-j}` (exact; no measure factor in the normalised convention)",
    "advection | `(v . grad w)hat_k = sum_j i (vhat_{k-j} . j) what_j`",
    "lattice sums | `A_4 = sum_{p != 0} |p|^{-4}`, `A_6 = sum_{p != 0} |p|^{-6}`; rigorous rational upper bounds by exact partial sum + shell-count tail (`24 m^2 + 2 <= 26 m^2` points on the sup-shell `m`), as in `adot_squared_upper`",
]

#: What the payload does and does not assert.  The checker rejects a payload
#: whose claims block has been altered or dropped.
CLAIMS = (
    "This payload certifies rational upper bounds for trilinear constants in "
    "the repository's own normalisation (the embedded conventions block).",
    "It is a constants certificate only: no statement about Navier-Stokes "
    "regularity, blow-up, or the Clay problem is made or implied.",
    "The infinite-dimensional commutator inequality behind G3 is analysis "
    "over l^2(Z^3), lives in the audited derivation note, and is never "
    "axiomatised in Lean.",
    "The literature comparison is verification only; no external decimal is "
    "an input to any bound in this payload.",
)

G3_STATEMENT = (
    "| < P(v.grad w), w >_{H3dot} | <= G3 ||v||_{H3dot} ||w||_{H3dot}^2 for "
    "mean-zero divergence-free real fields on T^3, with G3 <= 12 sqrt(A_4) "
    "(sections 2-4 of the derivation note: symmetrised trilinear identity, "
    "discrete Kato-Ponce commutator, l^1 * l^2 Young + Cauchy-Schwarz)"
)

K3_STATEMENT = (
    "||P(v.grad w)||_3 <= ||v.grad w||_3 <= 4 sqrt(A_4) ||v||_3 ||w||_3 "
    "+ 4 sqrt(A_6) ||v||_3 ||w||_4, from |k|^3 <= 4 (|p|^3 + |j|^3) and "
    "mode-wise contractivity of P (section 7 of the derivation note)"
)

K3_COUPLING_NOTE = (
    "4 sqrt(A_4) multiplies ||v||_3 ||w||_3 and 4 sqrt(A_6) multiplies "
    "||v||_3 ||w||_4: K_n couples order n to order n+1 as expected; the "
    "chain does not use K_3 directly, it is recorded because the same "
    "certificate machinery covers it at no cost"
)

MONOTONICITY_STATEMENT = (
    "upper(2N) <= upper(N): the exact mass added between the cuts N and 2N "
    "is at most the shell-count tail dropped, sum_{m=N+1}^{2N} 26 m^{-2} <= "
    "26/N - 26/(2N) for A_4 (and the m^{-4} analogue for A_6), so doubling "
    "the cut provably improves the bound; the doubled-cut values are "
    "recorded outward-rounded to dyadics at the payload precision (still "
    "rigorous upper bounds), so the inequality recorded here is "
    "round_up(upper(2N)) <= upper(N), which implies upper(2N) <= upper(N)"
)

LITERATURE_ROLE = (
    "USED FOR NOTHING: order-of-magnitude cross-check of the derivation "
    "only; no external number is an input to any bound in this payload"
)

#: The four crude steps of the derivation, each a finite algebraic
#: inequality (Lean-able, section 9 of the note) and each a recorded lever:
#: sharpening any one of them tightens ``G3`` without touching the others.
SHARPENING_LEVERS = (
    {
        "name": "mvt_cube_bound",
        "statement": (
            "| |k|^3 - |j|^3 | = | |k| - |j| | (|k|^2 + |k||j| + |j|^2) "
            "<= 3 |p| max(|k|,|j|)^2 <= 3 |p| (|j| + |p|)^2"
        ),
    },
    {
        "name": "square_expansion",
        "statement": "(|j| + |p|)^2 = |j|^2 + 2 |j||p| + |p|^2, bounded termwise",
    },
    {
        "name": "am_gm",
        "statement": "|p|^2 |j|^2 <= (|p| |j|^3 + |p|^3 |j|)/2",
    },
    {
        "name": "l1_young",
        "statement": (
            "l^1 * l^2 Young plus Cauchy-Schwarz: sum_p |p| |vhat_p| "
            "<= sqrt(A_4) ||v||_3, and likewise in the j slot"
        ),
    },
)


# --------------------------------------------------------------------------- #
# exact serialisation                                                          #
# --------------------------------------------------------------------------- #


def _stringify(value: Fraction) -> str:
    """``str`` with the interpreter's int-digit guard lifted temporarily.

    The exact partial sums have lcm-sized denominators (thousands of digits
    at cut 40, tens of thousands at 80) and CPython refuses int-to-str
    conversions past ``sys.get_int_max_str_digits()`` by default.  The limit
    is raised only around the conversion and always restored — never
    lowered, never left mutated.
    """
    digits = (
        max(abs(value.numerator).bit_length(), value.denominator.bit_length()) // 3
        + 8
    )
    limit = sys.get_int_max_str_digits()
    if limit == 0 or digits <= limit:
        return str(value)
    sys.set_int_max_str_digits(digits)
    try:
        return str(value)
    finally:
        sys.set_int_max_str_digits(limit)


def _parse_fraction(text: object) -> Fraction:
    """``Fraction(str)`` with the same temporary guard as :func:`_stringify`."""
    if isinstance(text, Fraction):
        return text
    if not isinstance(text, str):
        raise TypeError(f"expected an exact fraction string, got {type(text).__name__}")
    limit = sys.get_int_max_str_digits()
    if limit == 0 or len(text) <= limit:
        return Fraction(text)
    sys.set_int_max_str_digits(len(text) + 8)
    try:
        return Fraction(text)
    finally:
        sys.set_int_max_str_digits(limit)


# --------------------------------------------------------------------------- #
# lattice sums                                                                 #
# --------------------------------------------------------------------------- #


def _sum_fractions(items) -> Fraction:
    """Pairwise (balanced) sum: keeps most additions between small fractions.

    A left-to-right fold drags the huge common denominator of the running
    total through every addition; merging pairwise confines the big-integer
    work to the top of the tree.  The value is exactly the same.
    """
    values = list(items)
    if not values:
        return Fraction(0)
    while len(values) > 1:
        merged = [values[i] + values[i + 1] for i in range(0, len(values) - 1, 2)]
        if len(values) % 2:
            merged.append(values[-1])
        values = merged
    return values[0]


def _partials(lattice_cut: int) -> tuple[Fraction, Fraction]:
    r"""Exact ``(sum |p|^{-4}, sum |p|^{-6})`` over ``0 < |p|_inf <= N``.

    The full-cube loop of :func:`adot_squared_upper`, grouped by the integer
    shell ``n = |p|^2`` so the expensive ``Fraction`` additions run over the
    few thousand distinct shell values instead of every lattice point.  The
    grouped value is *identical* to the pointwise sum (tested against
    ``adot_squared_upper`` directly).
    """
    if lattice_cut < 1:
        raise ValueError("the lattice cut must be at least one")
    counts: dict[int, int] = {}
    rng = range(-lattice_cut, lattice_cut + 1)
    for kx in rng:
        x_sq = kx * kx
        for ky in rng:
            xy_sq = x_sq + ky * ky
            for kz in rng:
                n = xy_sq + kz * kz
                if n:
                    counts[n] = counts.get(n, 0) + 1
    partial_4 = _sum_fractions(Fraction(c, n * n) for n, c in sorted(counts.items()))
    partial_6 = _sum_fractions(Fraction(c, n**3) for n, c in sorted(counts.items()))
    return partial_4, partial_6


def a4_upper(lattice_cut: int = 40) -> Fraction:
    r"""A rigorous rational upper bound on ``A_4 = \sum_{p\ne 0}|p|^{-4}``.

    The same object as :func:`~ns_certificate_lab.fourier_torus.adot_squared_upper`
    — exact partial sum over ``|p|_\infty \le N`` plus the shell-count tail
    ``26/N`` — re-exported here shell-grouped for speed; the test suite
    asserts exact equality of the two implementations.  Tail: the shell
    ``|p|_\infty = m`` has ``(2m+1)^3 - (2m-1)^3 = 24m^2 + 2 \le 26m^2``
    points, each with ``|p| \ge |p|_\infty = m``, so it contributes at most
    ``26 m^{-2}``, and ``\sum_{m>N} m^{-2} \le \int_N^\infty x^{-2}dx = 1/N``.
    """
    partial_4, _ = _partials(lattice_cut)
    return partial_4 + Fraction(26, lattice_cut)


def a6_upper(lattice_cut: int = 40) -> Fraction:
    r"""A rigorous rational upper bound on ``A_6 = \sum_{p\ne 0}|p|^{-6}``.

    Exact partial sum over ``|p|_\infty \le N`` plus the tail ``9/N^3``:

    * shell count: the sup-shell ``|p|_\infty = m`` holds
      ``(2m+1)^3 - (2m-1)^3 = 24m^2 + 2`` points, and ``24m^2 + 2 \le 26m^2``
      for ``m \ge 1`` because ``2 \le 2m^2``;
    * each such point has ``|p| \ge |p|_\infty = m``, hence
      ``|p|^{-6} \le m^{-6}``, so the shell contributes at most
      ``26 m^2 \cdot m^{-6} = 26\,m^{-4}``;
    * for ``m > N``, ``m^{-4} \le \int_{m-1}^{m} x^{-4}\,dx``, so
      ``\sum_{m>N} m^{-4} \le \int_N^\infty x^{-4}\,dx = 1/(3N^3)``;
    * therefore the tail is at most ``26/(3N^3) \le 9/N^3`` since
      ``26/3 \le 9``.
    """
    _, partial_6 = _partials(lattice_cut)
    return partial_6 + Fraction(9, lattice_cut**3)


def _sqrt_upper(value: Fraction, *, bits: int) -> Fraction:
    """Outward-rounded rational upper bound on ``sqrt(value)``."""
    return sqrt_interval(Interval(value, value), bits=bits).upper


def g3_upper(lattice_cut: int = 40, *, bits: int = DEFAULT_PRECISION_BITS) -> Fraction:
    """``G3 <= 12 sqrt(A_4)`` with both factors outward-rounded upper bounds."""
    return 12 * _sqrt_upper(a4_upper(lattice_cut), bits=bits)


# --------------------------------------------------------------------------- #
# per-trajectory band constants (sections 5a, 5b of the note)                  #
# --------------------------------------------------------------------------- #


def c_kato_upper(field: TrigVector, *, bits: int = DEFAULT_PRECISION_BITS) -> Fraction:
    r"""Section 5a: the band Kato constant ``C_kato(u_a)`` of the note.

    Implements exactly the inequality

    .. math::  |\langle P(u_a\cdot\nabla w),\, w\rangle_{\dot H^3}|
               \;\le\; C_{kato}(u_a)\, \|w\|_{\dot H^3}^2 ,\qquad
               C_{kato}(u_a) = 6 \sum_{\text{canonical } p}
               (|p| + |p|^3)\,\sqrt{|a_p|^2 + |b_p|^2}

    (note, section 5a): Steps 1–2 of the commutator bound rerun with
    ``\hat v`` supported on the band and the ``p``-sums kept exact, using
    ``|\hat u_{a,\pm p}| = \sqrt{|a_p|^2+|b_p|^2}/2`` and counting both
    signs.  Everything is a rational upper bound: ``|p| + |p|^3 =
    |p|\,(1 + |p|^2)`` is bounded by ``sqrt_upper(|p|^2)\,(1 + |p|^2)`` with
    ``|p|^2`` an exact integer, and the coefficient magnitude by
    :func:`~ns_certificate_lab.fourier_torus.euclidean_norm_upper` over the
    six components of ``(a_p, b_p)`` (interval coefficients — e.g. a Galerkin
    box field — contribute their interval magnitudes, so the result is valid
    over the whole slab box).
    """
    total = Fraction(0)
    for k, (a, b) in field.coefficient_table().items():
        if k == (0, 0, 0):
            continue
        n = _norm_sq(k)
        norm_upper = _sqrt_upper(Fraction(n), bits=bits)
        coefficient = euclidean_norm_upper(list(a) + list(b), bits=bits)
        total += norm_upper * (1 + n) * coefficient
    return 6 * total


def c_shift_upper(field: TrigVector, *, bits: int = DEFAULT_PRECISION_BITS) -> Fraction:
    r"""Section 5b: the band shift constant ``C_shift(u_a)`` of the note.

    Implements exactly the inequality

    .. math::  |\langle w\cdot\nabla u_a,\, w\rangle_{\dot H^3}|
               \;\le\; C_{shift}(u_a)\, \|w\|_{\dot H^3}^2 ,\qquad
               C_{shift}(u_a) = \sum_{\text{canonical } j}
               |j|\,(1 + |j|)^3\,\sqrt{|a_j|^2 + |b_j|^2}

    (note, section 5b): no cancellation is used — the advected factor is
    band-limited — via
    ``\sum_k |k|^6 |\hat w_{k-j}||\hat w_k| \le \|w\|_3\,
    \sqrt{\sum_m |m+j|^6 |\hat w_m|^2}`` and
    ``\sup_{|m|\ge 1} |m+j|^3/|m|^3 = (1+|j|)^3``.  With
    ``u = sqrt\_upper(|j|^2) \ge |j|`` the factor ``|j|(1+|j|)^3`` is bounded
    by ``u\,(1+u)^3`` (monotone in the root), so the result is a rational
    upper bound, valid verbatim for interval coefficients.
    """
    total = Fraction(0)
    for k, (a, b) in field.coefficient_table().items():
        if k == (0, 0, 0):
            continue
        n = _norm_sq(k)
        norm_upper = _sqrt_upper(Fraction(n), bits=bits)
        coefficient = euclidean_norm_upper(list(a) + list(b), bits=bits)
        total += norm_upper * (1 + norm_upper) ** 3 * coefficient
    return total


# --------------------------------------------------------------------------- #
# the certificate                                                              #
# --------------------------------------------------------------------------- #


def build_kato_certificate(lattice_cut: int = 40, bits: int = 64) -> dict:
    """Build the ``kato_h3_constant`` payload.

    Every number in the payload is an exact ``Fraction`` rendered as a
    string; the floats are display-only conveniences that the checker
    verifies against the exact values.  The monotonicity witness evaluates
    the same bounds at ``2 * lattice_cut`` and records that the larger cut
    gives the smaller (better) upper bound — a fact that is provable in
    general (see ``MONOTONICITY_STATEMENT``) and re-checked concretely.
    """
    if lattice_cut < 1:
        raise ValueError("the lattice cut must be at least one")
    if bits < 8:
        raise ValueError("refusing fewer than 8 bits of outward rounding")

    partial_4, partial_6 = _partials(lattice_cut)
    tail_4 = Fraction(26, lattice_cut)
    tail_6 = Fraction(9, lattice_cut**3)
    a4 = partial_4 + tail_4
    a6 = partial_6 + tail_6
    sqrt_a4 = _sqrt_upper(a4, bits=bits)
    sqrt_a6 = _sqrt_upper(a6, bits=bits)
    g3 = 12 * sqrt_a4

    double_cut = 2 * lattice_cut
    partial_4_d, partial_6_d = _partials(double_cut)
    # Outward-rounded dyadics: still rigorous upper bounds, but with bounded
    # string size (the exact doubled-cut fractions run to tens of thousands
    # of digits and would bloat the payload for no rigour gain).
    a4_double = (
        Interval.exact(partial_4_d + Fraction(26, double_cut))
        .round_outward(bits)
        .upper
    )
    a6_double = (
        Interval.exact(partial_6_d + Fraction(9, double_cut**3))
        .round_outward(bits)
        .upper
    )
    g3_double = 12 * _sqrt_upper(a4_double, bits=bits)
    if a4_double > a4 or a6_double > a6 or g3_double > g3:
        raise AssertionError(
            "monotonicity witness failed: doubling the lattice cut worsened "
            "a bound, which the shell-count tail makes impossible"
        )

    return {
        "schema": SCHEMA,
        "lattice_cut": lattice_cut,
        "precision_bits": bits,
        "claims": list(CLAIMS),
        "conventions": list(CONVENTIONS),
        "derivation": DERIVATION_NOTE,
        "a4": {
            "definition": "A_4 = sum_{p != 0} |p|^{-4}",
            "partial_region": f"exact rational sum over 0 < |p|_inf <= {lattice_cut}",
            "partial_sum": _stringify(partial_4),
            "tail": str(tail_4),
            "tail_formula": (
                "26/N: the sup-shell m holds 24 m^2 + 2 <= 26 m^2 points, each "
                "with |p| >= m, and sum_{m>N} m^{-2} <= 1/N"
            ),
            "upper": _stringify(a4),
            "upper_float": float(a4),
        },
        "a6": {
            "definition": "A_6 = sum_{p != 0} |p|^{-6}",
            "partial_region": f"exact rational sum over 0 < |p|_inf <= {lattice_cut}",
            "partial_sum": _stringify(partial_6),
            "tail": str(tail_6),
            "tail_formula": (
                "9/N^3: the sup-shell m contributes at most 26 m^2 . m^{-6} = "
                "26 m^{-4}, and sum_{m>N} m^{-4} <= 1/(3 N^3), with 26/3 <= 9"
            ),
            "upper": _stringify(a6),
            "upper_float": float(a6),
        },
        "g3": {
            "statement": G3_STATEMENT,
            "multiplier": 12,
            "sqrt_a4_upper": str(sqrt_a4),
            "upper": str(g3),
            "upper_float": float(g3),
        },
        "k3": {
            "statement": K3_STATEMENT,
            "multiplier": 4,
            "sqrt_a6_upper": str(sqrt_a6),
            "coeff_v3_w3": str(4 * sqrt_a4),
            "coeff_v3_w4": str(4 * sqrt_a6),
            "coupling_note": K3_COUPLING_NOTE,
        },
        "monotonicity": {
            "statement": MONOTONICITY_STATEMENT,
            "lattice_cut": lattice_cut,
            "double_cut": double_cut,
            "a4_upper_at_cut": _stringify(a4),
            "a4_upper_at_double_cut": str(a4_double),
            "a6_upper_at_cut": _stringify(a6),
            "a6_upper_at_double_cut": str(a6_double),
            "g3_upper_at_cut": str(g3),
            "g3_upper_at_double_cut": str(g3_double),
            "improves": True,
        },
        "sharpening_levers": [dict(lever) for lever in SHARPENING_LEVERS],
        "literature_comparison": {
            "verification_only": True,
            "role": LITERATURE_ROLE,
            "conversion": (
                "Morosi-Pizzocchero norms on (0, 2pi)^3 with unnormalised "
                "Lebesgue measure satisfy ||u||_{MP,n}^2 = (2 pi)^3 "
                "||u||_{ours,n}^2, so a trilinear constant converts as "
                "G_ours = (2 pi)^{3/2} G_MP"
            ),
            "conversion_factor_approx": 15.75,
            "mp_g3_approx": 0.438,
            "mp_g3_in_repo_normalisation_approx": 6.9,
            "comparison": (
                "our provable bound is about 7x cruder than the converted "
                "sharp value — the price of the four crude steps, each a "
                "recorded sharpening lever; agreement in order of magnitude "
                "checks the normalisation and is used for nothing else"
            ),
        },
    }


# --------------------------------------------------------------------------- #
# the independent checker                                                      #
# --------------------------------------------------------------------------- #


def _checker_partials(lattice_cut: int) -> tuple[Fraction, Fraction]:
    """The checker's own lattice loop — not the builder's.

    Traverses the closed octant ``0 <= kx, ky, kz <= N`` and weights each
    point by ``2^(#nonzero components)``, which enumerates the full cube
    ``|p|_inf <= N`` exactly once per point.  A deliberately different
    traversal from the builder's full-cube loop, so a bug in either shows up
    as a mismatch.
    """
    if lattice_cut < 1:
        raise ValueError("the lattice cut must be at least one")
    counts: dict[int, int] = {}
    for kx in range(lattice_cut + 1):
        weight_x = 2 if kx else 1
        x_sq = kx * kx
        for ky in range(lattice_cut + 1):
            weight_xy = weight_x * (2 if ky else 1)
            xy_sq = x_sq + ky * ky
            for kz in range(lattice_cut + 1):
                n = xy_sq + kz * kz
                if n == 0:
                    continue
                counts[n] = counts.get(n, 0) + weight_xy * (2 if kz else 1)
    partial_4 = _sum_fractions(Fraction(c, n * n) for n, c in sorted(counts.items()))
    partial_6 = _sum_fractions(Fraction(c, n**3) for n, c in sorted(counts.items()))
    return partial_4, partial_6


def verify_kato_certificate(payload: dict) -> dict:
    """Independently re-check a ``kato_h3_constant`` payload.

    The checker recomputes the ``A_4``/``A_6`` partial sums with its own
    octant loop (it never calls the builder's helper), re-derives the tails
    from their closed forms, re-runs the outward-rounded square roots at the
    recorded precision, checks the ``12``/``4`` multipliers, the
    monotonicity witness, the verbatim conventions block, the claims block,
    and the verification-only status of the literature comparison.
    Adversarial input gets a verdict, never an exception.
    """
    try:
        return _verify(payload)
    except Exception as error:  # noqa: BLE001 -- adversarial input is a failure
        return {
            "verified": False,
            "failures": [f"payload rejected while parsing: {error}"],
        }


def _verify(payload: dict) -> dict:
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(payload.get("schema") == SCHEMA, "wrong schema tag")
    require(
        payload.get("conventions") == CONVENTIONS,
        "the conventions block is missing or altered — the constant is "
        "meaningless without the normalisation it was derived in",
    )
    require(
        list(payload.get("claims", [])) == list(CLAIMS),
        "the claims block is missing or altered",
    )
    require(
        payload.get("derivation") == DERIVATION_NOTE,
        "the derivation reference is missing or altered",
    )

    lattice_cut = int(payload["lattice_cut"])
    bits = int(payload["precision_bits"])
    require(lattice_cut >= 1, "the lattice cut must be at least one")
    require(bits >= 8, "too few precision bits recorded")

    # -- A4 / A6: own loop, closed-form tails --------------------------------- #
    partial_4, partial_6 = _checker_partials(lattice_cut)
    a4_block = payload["a4"]
    a6_block = payload["a6"]
    require(
        _parse_fraction(a4_block["partial_sum"]) == partial_4,
        "the recorded A_4 partial sum does not equal the checker's own "
        "lattice sum",
    )
    require(
        _parse_fraction(a4_block["tail"]) == Fraction(26, lattice_cut),
        "the recorded A_4 tail is not 26/N",
    )
    a4 = _parse_fraction(a4_block["upper"])
    require(
        a4 == partial_4 + Fraction(26, lattice_cut),
        "the recorded A_4 upper bound is not partial sum + tail",
    )
    require(
        _parse_fraction(a6_block["partial_sum"]) == partial_6,
        "the recorded A_6 partial sum does not equal the checker's own "
        "lattice sum",
    )
    require(
        _parse_fraction(a6_block["tail"]) == Fraction(9, lattice_cut**3),
        "the recorded A_6 tail is not 9/N^3",
    )
    a6 = _parse_fraction(a6_block["upper"])
    require(
        a6 == partial_6 + Fraction(9, lattice_cut**3),
        "the recorded A_6 upper bound is not partial sum + tail",
    )
    for name, block, exact in (("A_4", a4_block, a4), ("A_6", a6_block, a6)):
        require(
            abs(Fraction(float(block["upper_float"])) - exact) <= Fraction(1, 10**6),
            f"the display float for {name} disagrees with the exact value",
        )

    # -- G3: the sqrt relation and the 12 ------------------------------------- #
    g3_block = payload["g3"]
    require(g3_block.get("multiplier") == 12, "the G3 multiplier is not 12")
    sqrt_a4 = _parse_fraction(g3_block["sqrt_a4_upper"])
    require(
        sqrt_a4 >= 0 and sqrt_a4 * sqrt_a4 >= a4,
        "the recorded sqrt(A_4) is not an upper root of the recorded A_4",
    )
    require(
        sqrt_a4 == sqrt_interval(Interval(a4, a4), bits=bits).upper,
        "the recorded sqrt(A_4) does not match outward rounding at the "
        "recorded precision",
    )
    g3 = _parse_fraction(g3_block["upper"])
    require(g3 == 12 * sqrt_a4, "the recorded G3 is not 12 * sqrt(A_4)")
    require(
        abs(Fraction(float(g3_block["upper_float"])) - g3) <= Fraction(1, 10**6),
        "the display float for G3 disagrees with the exact value",
    )
    require(g3_block.get("statement") == G3_STATEMENT, "the G3 statement was altered")

    # -- K3: the pair and the 4 ----------------------------------------------- #
    k3_block = payload["k3"]
    require(k3_block.get("multiplier") == 4, "the K3 multiplier is not 4")
    sqrt_a6 = _parse_fraction(k3_block["sqrt_a6_upper"])
    require(
        sqrt_a6 >= 0 and sqrt_a6 * sqrt_a6 >= a6,
        "the recorded sqrt(A_6) is not an upper root of the recorded A_6",
    )
    require(
        sqrt_a6 == sqrt_interval(Interval(a6, a6), bits=bits).upper,
        "the recorded sqrt(A_6) does not match outward rounding at the "
        "recorded precision",
    )
    require(
        _parse_fraction(k3_block["coeff_v3_w3"]) == 4 * sqrt_a4,
        "the K3 coefficient on ||v||_3 ||w||_3 is not 4 sqrt(A_4)",
    )
    require(
        _parse_fraction(k3_block["coeff_v3_w4"]) == 4 * sqrt_a6,
        "the K3 coefficient on ||v||_3 ||w||_4 is not 4 sqrt(A_6)",
    )
    coupling = str(k3_block.get("coupling_note", ""))
    require(
        "||v||_3 ||w||_3" in coupling and "||v||_3 ||w||_4" in coupling,
        "the K3 coupling note no longer names both norm pairings",
    )
    require(k3_block.get("statement") == K3_STATEMENT, "the K3 statement was altered")

    # -- monotonicity witness: recompute the doubled cut ---------------------- #
    mono = payload["monotonicity"]
    double_cut = int(mono["double_cut"])
    require(
        int(mono["lattice_cut"]) == lattice_cut and double_cut == 2 * lattice_cut,
        "the monotonicity witness does not use the cut and its double",
    )
    partial_4_d, partial_6_d = _checker_partials(double_cut)
    a4_double = (
        Interval.exact(partial_4_d + Fraction(26, double_cut))
        .round_outward(bits)
        .upper
    )
    a6_double = (
        Interval.exact(partial_6_d + Fraction(9, double_cut**3))
        .round_outward(bits)
        .upper
    )
    g3_double = 12 * sqrt_interval(Interval(a4_double, a4_double), bits=bits).upper
    require(
        _parse_fraction(mono["a4_upper_at_cut"]) == a4,
        "the monotonicity witness disagrees with the recorded A_4",
    )
    require(
        _parse_fraction(mono["a4_upper_at_double_cut"]) == a4_double,
        "the recorded doubled-cut A_4 does not equal the checker's own value",
    )
    require(
        _parse_fraction(mono["a6_upper_at_cut"]) == a6,
        "the monotonicity witness disagrees with the recorded A_6",
    )
    require(
        _parse_fraction(mono["a6_upper_at_double_cut"]) == a6_double,
        "the recorded doubled-cut A_6 does not equal the checker's own value",
    )
    require(
        _parse_fraction(mono["g3_upper_at_cut"]) == g3,
        "the monotonicity witness disagrees with the recorded G3",
    )
    require(
        _parse_fraction(mono["g3_upper_at_double_cut"]) == g3_double,
        "the recorded doubled-cut G3 does not equal the checker's own value",
    )
    require(
        a4_double <= a4 and a6_double <= a6 and g3_double <= g3,
        "monotonicity fails: the doubled cut did not improve every bound",
    )
    require(mono.get("improves") is True, "the monotonicity verdict was altered")

    # -- the crude steps stay recorded ---------------------------------------- #
    lever_names = [lever.get("name") for lever in payload.get("sharpening_levers", [])]
    require(
        lever_names == [lever["name"] for lever in SHARPENING_LEVERS],
        "the four sharpening levers (mvt_cube_bound, square_expansion, am_gm, "
        "l1_young) are no longer all recorded",
    )

    # -- the literature comparison must stay decorative ------------------------ #
    literature = payload["literature_comparison"]
    require(
        literature.get("verification_only") is True,
        "the literature comparison is no longer marked verification_only — "
        "an external number may have been promoted to an input",
    )
    require(
        "USED FOR NOTHING" in str(literature.get("role", "")),
        "the literature comparison role disclaimer was altered",
    )

    return {
        "verified": not failures,
        "failures": failures,
        "g3_upper": str(g3) if not failures else None,
        "g3_upper_float": float(g3) if not failures else None,
    }
