r"""The **smoothly** truncated coherent critical-spectrum family, and the one
continuum input the capacity lower bound needs.

Why this module exists
----------------------
The Proposition of ``docs/paper_lambda_dichotomy/theorem_statement.md`` is
conditional on Hypothesis (L\*), a capacity lower bound
``\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3`` for the **sharply**
truncated family ``\hat u_N(k)=P_kv_0/|k|^2`` on ``1\le|k|\le N``.  The
Proposition itself, however, only needs *some* family of real, mean-zero,
divergence-free trigonometric fields whose capacity ``K`` grows at least like
a positive power of ``N_0^2`` along a sequence with bounded logarithmic gaps.
That freedom is what this module exploits: it replaces the sharp indicator by
a smooth radial weight ``\chi``, for which the whole chain

    lattice family  ->  inner rescaling  ->  continuum profile  ->  duality

closes, with exactly one nontrivial continuum input.

The family
----------
For a band ``N`` and a seed ``v_0\in\mathbb Z^3\setminus\{0\}``,

.. math::
   \hat u_N(k)=\chi(|k|/N)\,\frac{P_kv_0}{|k|^2},\qquad k\ne0,
   \qquad P_k=I-\frac{k\otimes k}{|k|^2},

with ``\chi`` radial, ``\chi=1`` on ``[0,1/2]``, ``\mathrm{supp}\,\chi
\subset[0,1]``, ``0\le\chi\le1``.  The coefficients are real and even in
``k``, so ``u_N`` is a real cosine field; ``k\cdot P_kv_0=0``, so it is
divergence free.  Because ``\chi`` is *radial*, the band-symmetry identity
(Lemma 9 of the complete proof: ``\sum_{B}k_ik_jf(|k|^2)=\tfrac13\delta_{ij}
\sum_B|k|^2f(|k|^2)`` for any radial ``f\ge0`` on any sign-flip- and
permutation-invariant band) applies verbatim with
``f(|k|^2)=\chi(|k|/N)^2|k|^{-6}`` etc., so the exact moment laws survive:

.. math::
   H_0=\tfrac23\|v_0\|^2\sum_k\frac{\chi(|k|/N)^2}{|k|^4},\qquad
   H_1=\tfrac23\|v_0\|^2\sum_k\frac{\chi(|k|/N)^2}{|k|^2},\qquad
   u_N(0)=\tfrac23v_0\sum_k\frac{\chi(|k|/N)}{|k|^2}.

:func:`exact_moment_laws` verifies these with *zero* residual in exact
rational arithmetic.

The choice of ``\chi``
----------------------
``\chi`` must be evaluated at ``|k|/N``, and ``|k|`` is generically
irrational.  The family is therefore made exactly rational by taking ``\chi``
to be a **polynomial in ``r^2``**, since ``r^2=|k|^2/N^2\in\mathbb Q``:

.. math::
   \chi(r)=1-S\!\left(\frac{r^2-1/4}{3/4}\right),\qquad
   S(s)=126s^5-420s^6+540s^7-315s^8+70s^9,

with ``S`` clamped to ``[0,1]``.  ``S`` is the degree-9 smoothstep: it is the
normalised integral of ``s^4(1-s)^4``, hence increases monotonically from
``0`` to ``1`` on ``[0,1]`` with four vanishing derivatives at each end.  So
``\chi\in C^4(\mathbb R)``, ``\chi\equiv1`` on ``[0,1/2]``, ``\chi\equiv0``
on ``[1,\infty)``, ``0\le\chi\le1`` — and ``\chi(|k|/N)\in\mathbb Q`` for
every lattice point, which makes the *entire* family, its moments, and its
full nonlinear norm exactly rational.  (The sharp family is not rational in
this way; it is rational only because its weight is ``0`` or ``1``.)

The inner rescaling
-------------------
Write ``F(\xi)=\chi(|\xi|)P_\xi v_0/|\xi|^2`` on ``\mathbb R^3``; then
``\hat u_N(k)=N^{-2}F(k/N)`` exactly, because ``P_k=P_{k/N}``.  Let

.. math::
   V(y)=\int_{\mathbb R^3}F(\xi)e^{i\xi\cdot y}\,d\xi .

``F`` is compactly supported with an integrable ``|\xi|^{-2}`` singularity at
the origin, so ``F\in L^1`` and ``V`` is continuous, bounded, real and even.

**Radial reduction (exact, used by** :func:`radial_profiles` **).**  Because
``F`` is a radial weight times ``P_\xi v_0``,

.. math::
   V(y)=A(r)\,v_0+B(r)\,\hat y(\hat y\cdot v_0),\qquad r=|y|,
   \\
   A=\tfrac12(g+h),\quad B=\tfrac12(g-3h),
   \\
   g(r)=4\pi\!\int_0^1\!\chi(\rho)\,j_0(\rho r)\,d\rho,\qquad
   h(r)=4\pi\!\int_0^1\!\chi(\rho)\,\mathcal B(\rho r)\,d\rho,

where ``j_0(s)=\sin s/s`` and ``\mathcal B(s)=\int_0^1\mu^2\cos(s\mu)d\mu``.
Derivation: ``\int\chi|\xi|^{-2}e^{i\xi y}d\xi=g``; the tensor
``T_{ij}=\int\chi\,\xi_i\xi_j|\xi|^{-4}e^{i\xi y}d\xi=A_T\delta_{ij}
+B_T\hat y_i\hat y_j`` has trace ``3A_T+B_T=g`` and contraction
``A_T+B_T=h``; ``V=gv_0-T v_0``.  Every integrand above is **smooth on
``[0,1]``** — the ``|\xi|^{-2}`` singularity is cancelled by the ``\rho^2``
spherical Jacobian — so ``g`` and ``h`` are computed to machine precision by
Gauss–Legendre quadrature, and their ``r``-derivatives by differentiating
under the integral sign.

**Untruncated limit (closed form).**  For ``\chi\equiv1`` one has ``g=2\pi^2/r``
and ``h\equiv0`` (the second follows from
``\int_0^\infty[\sin s/s+2\cos s/s^2-2\sin s/s^3]ds=0``), so

.. math::
   V_\infty(y)=\frac{\pi^2}{|y|}\bigl(v_0+\hat y(\hat y\cdot v_0)\bigr).

Since ``\chi=1`` near ``0``, ``F-F_\infty=-(1-\chi)F_\infty`` is a smooth
symbol of order ``-2``, whose transform is rapidly decreasing away from the
origin; hence ``V(y)=V_\infty(y)+O(|y|^{-\infty})`` with all derivatives.
:func:`radial_profiles` confirms ``rA(r)\to\pi^2`` and ``rB(r)\to\pi^2``.

The decisive input (V-NONDEG)
-----------------------------
The duality step needs ``\mathbb P(V\cdot\nabla V)\not\equiv0``, i.e. ``V`` is
not a stationary Euler field.  ``\mathbb P W\equiv0`` iff ``W`` is a gradient
iff ``\nabla\times W\equiv0``, and **curl is local**, so it suffices to
exhibit one point where ``\nabla\times(V\cdot\nabla V)\ne0``.  ``V`` is
axisymmetric about ``v_0`` with no swirl; writing
``V=\alpha(r)\cos\theta\,\hat e_r-A(r)\sin\theta\,\hat e_\theta`` with
``\alpha=A+B=g-h`` and ``\theta`` the polar angle from ``v_0``,

.. math::
   \omega=\nabla\times V=W(r)\sin\theta\,\hat e_\varphi,\qquad
   W=\frac{\alpha-(rA)'}{r},
   \\
   \nabla\times(V\cdot\nabla V)=Z(r)\sin\theta\cos\theta\,\hat e_\varphi,
   \qquad
   Z=\frac{(rW\alpha)'-2WA}{r},

everything carrying an overall factor ``\|v_0\|^2``.  For ``V_\infty`` this
gives the closed form ``Z_\infty(r)=-12\pi^4\|v_0\|^2/r^4\ne0``, and the
asymptotics above force ``Z(r)=-12\pi^4\|v_0\|^2r^{-4}+O(r^{-\infty})``.
:func:`curl_coefficient` evaluates ``Z`` for the truncated ``\chi``;
:func:`asymptotic_curl_coefficient` is the closed form it must approach.

Independent routes to the same fact
-----------------------------------
* :func:`exact_nonlinear_mode` — the **exact rational** Fourier coefficient
  ``\widehat{u_N\cdot\nabla u_N}(k)=i\,b_k`` at a single mode.  Since
  ``\hat u_N(k)=N^{-2}F(k/N)``, the discrete convolution is a Riemann sum:
  ``b_k\to C(\zeta)`` with ``\zeta=k/N`` and
  ``C(\zeta)=\int(F(\zeta-b)\cdot b)F(b)\,db=-i\widehat{V\cdot\nabla V}(\zeta)``.
  No ``N``-power normalisation appears; the discrete number converges
  directly to the continuum one.
* :func:`convolution_coefficient` — the same ``C(\zeta)`` by direct 3-D
  quadrature of the convolution integral, split by a smooth partition of
  unity so that each half has a single singular point removed by the
  spherical Jacobian.

Labelling discipline
--------------------
Functions whose name starts with ``exact_`` return :class:`fractions.Fraction`
(or Python integers over a stated common denominator) and are **exact**.
Everything else is ``float`` (IEEE double).  No interval arithmetic is used
here; the exact lane is the certificate lane.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Iterable, Sequence

import numpy as np

__all__ = [
    "DEFAULT_V0",
    "SMOOTHSTEP_COEFFICIENTS",
    "asymptotic_curl_coefficient",
    "band_waves",
    "capacity_metrics",
    "chi_exact",
    "chi_float",
    "convolution_coefficient",
    "curl_coefficient",
    "dealias_grid",
    "exact_moment_laws",
    "exact_nonlinear_mode",
    "exact_nonlinear_norm_sq",
    "exact_scaled_coefficients",
    "meridional_profiles",
    "radial_profiles",
    "spectral_field",
    "vector_profile",
]

DEFAULT_V0: tuple[int, int, int] = (1, 2, 3)

#: Above this many decimal digits an exact rational is reported by its bit
#: sizes rather than digit-by-digit.  The value is *not* rounded: the
#: computation stays exact, only the transcript is abbreviated (CPython also
#: refuses ``str()`` on integers beyond ~4300 digits by default).
_RATIONAL_DIGIT_CAP = 1500


def _rational_record(value: Fraction) -> object:
    """Serialise an exact rational, abbreviating astronomically long ones."""

    if max(value.numerator.bit_length(), value.denominator.bit_length()) < int(
        _RATIONAL_DIGIT_CAP * 3.33
    ):
        return str(value)
    return {
        "exact_but_abbreviated": True,
        "numerator_bits": value.numerator.bit_length(),
        "denominator_bits": value.denominator.bit_length(),
        "sign": (value > 0) - (value < 0),
        "float": float(value),
    }

#: Coefficients of ``S(s)=126s^5-420s^6+540s^7-315s^8+70s^9`` (degree-9
#: smoothstep), indexed by the power ``5..9``.
SMOOTHSTEP_COEFFICIENTS: tuple[int, ...] = (126, -420, 540, -315, 70)

#: ``\chi=1`` for ``r^2\le`` this, ``\chi=0`` for ``r^2\ge1``.
_INNER_SQ = Fraction(1, 4)
_TRANSITION = Fraction(3, 4)


# --------------------------------------------------------------------------
# the radial weight chi
# --------------------------------------------------------------------------


def _smoothstep_numerator(numerator: int, denominator: int) -> int:
    r"""``S(p/q)\cdot q^9`` as an exact integer, for ``0\le p\le q``."""

    total = 0
    for index, coefficient in enumerate(SMOOTHSTEP_COEFFICIENTS):
        power = 5 + index
        total += coefficient * numerator**power * denominator ** (9 - power)
    return total


def chi_numerator(norm_sq: int, band: int) -> tuple[int, int]:
    r"""Return ``(C_k, D)`` with ``\chi(|k|/N)=C_k/D`` exactly.

    ``s=(|k|^2/N^2-1/4)/(3/4)=(4|k|^2-N^2)/(3N^2)``, so with
    ``q=3N^2`` and ``p=4|k|^2-N^2`` the weight is
    ``1-S(p/q)=(q^9-S(p/q)q^9)/q^9`` — a ratio of integers, and the same
    denominator ``D=q^9`` serves every lattice point of the band.
    """

    if band <= 0:
        raise ValueError("the band must be positive")
    q = 3 * band * band
    d = q**9
    p = 4 * norm_sq - band * band
    if p <= 0:
        return d, d
    if p >= q:
        return 0, d
    return d - _smoothstep_numerator(p, q), d


def chi_exact(norm_sq: int, band: int) -> Fraction:
    r"""``\chi(|k|/N)`` as an exact :class:`~fractions.Fraction`."""

    numerator, denominator = chi_numerator(norm_sq, band)
    return Fraction(numerator, denominator)


def chi_float(radius) -> np.ndarray:
    r"""``\chi(r)`` in IEEE double, vectorised over ``r``."""

    t = np.asarray(radius, dtype=np.float64) ** 2
    s = np.clip((t - 0.25) / 0.75, 0.0, 1.0)
    step = np.zeros_like(s)
    for index, coefficient in enumerate(SMOOTHSTEP_COEFFICIENTS):
        step = step + coefficient * s ** (5 + index)
    return 1.0 - step


# --------------------------------------------------------------------------
# the exact lattice family
# --------------------------------------------------------------------------


def band_waves(band: int) -> list[tuple[int, int, int]]:
    r"""Every lattice point with ``1\le|k|^2\le N^2`` (both signs).

    Points with ``|k|=N`` are retained even though ``\chi`` annihilates them;
    dropping them is left to the callers that care, so that the band is the
    same set as the sharp family's.
    """

    waves: list[tuple[int, int, int]] = []
    for kx in range(-band, band + 1):
        for ky in range(-band, band + 1):
            for kz in range(-band, band + 1):
                norm_sq = kx * kx + ky * ky + kz * kz
                if norm_sq == 0 or norm_sq > band * band:
                    continue
                waves.append((kx, ky, kz))
    return waves


def _lattice_lcm(waves: Iterable[tuple[int, int, int]]) -> int:
    value = 1
    for wave in waves:
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        value = value * norm_sq // math.gcd(value, norm_sq)
    return value


def exact_scaled_coefficients(
    band: int, v0: Sequence[int] = DEFAULT_V0
) -> tuple[int, dict[tuple[int, int, int], tuple[int, int, int]]]:
    r"""Return ``(M, A)`` with ``\hat u_N(k)=A_k/M`` exactly and ``A_k`` integral.

    ``M=D\,L^2`` where ``D`` is the common ``\chi`` denominator and
    ``L=\mathrm{lcm}\{|k|^2\}`` over the band, so
    ``A_k=C_k\,(v_0|k|^2-k(k\cdot v_0))\,(L/|k|^2)^2`` is an integer vector.
    Pure integer arithmetic downstream: no rounding anywhere, and no
    :class:`~fractions.Fraction` normalisation cost inside the hot loops.
    """

    waves = band_waves(band)
    lattice = _lattice_lcm(waves)
    _, chi_denominator = chi_numerator(1, band)
    scale = chi_denominator * lattice * lattice
    table: dict[tuple[int, int, int], tuple[int, int, int]] = {}
    for wave in waves:
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        weight, _ = chi_numerator(norm_sq, band)
        if weight == 0:
            continue
        dot = sum(c * v for c, v in zip(wave, v0))
        factor = weight * (lattice // norm_sq) ** 2
        vector = tuple((v0[i] * norm_sq - wave[i] * dot) * factor for i in range(3))
        if sum(c * v for c, v in zip(wave, vector)) != 0:
            raise AssertionError("exact coefficient is not divergence free")
        table[wave] = vector  # type: ignore[assignment]
    return scale, table


def exact_moment_laws(band: int, v0: Sequence[int] = DEFAULT_V0) -> dict[str, object]:
    r"""Exact verification of the family laws for the *smooth* weight.

    Checks, with residual required to be exactly zero:
    ``H_0=\tfrac23\|v_0\|^2\sum\chi^2|k|^{-4}``,
    ``H_1=\tfrac23\|v_0\|^2\sum\chi^2|k|^{-2}``,
    ``u_N(0)=\tfrac23v_0\sum\chi|k|^{-2}``,
    together with realness (the coefficient table is even in ``k`` and
    integral, so ``u`` is a pure cosine field) and ``k\cdot\hat u_N(k)=0``.
    """

    scale, table = exact_scaled_coefficients(band, v0)
    scale_sq = Fraction(scale) ** 2
    v0_norm_sq = sum(c * c for c in v0)

    h0 = Fraction(0)
    h1 = Fraction(0)
    chi_over_2 = Fraction(0)  # sum chi/|k|^2
    chi_sq_over_2 = Fraction(0)  # sum chi^2/|k|^2
    chi_sq_over_4 = Fraction(0)  # sum chi^2/|k|^4
    point_value = [Fraction(0), Fraction(0), Fraction(0)]
    hermitian = True
    divergence_free = True

    for wave, vector in table.items():
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        modulus_sq = Fraction(sum(c * c for c in vector)) / scale_sq
        h0 += modulus_sq
        h1 += norm_sq * modulus_sq
        weight = chi_exact(norm_sq, band)
        chi_over_2 += weight / norm_sq
        chi_sq_over_2 += weight * weight / norm_sq
        chi_sq_over_4 += weight * weight / (norm_sq * norm_sq)
        for i in range(3):
            point_value[i] += Fraction(vector[i], scale)
        mirror = table.get((-wave[0], -wave[1], -wave[2]))
        if mirror is None or mirror != vector:
            hermitian = False
        if sum(c * v for c, v in zip(wave, vector)) != 0:
            divergence_free = False

    prefactor = Fraction(2, 3) * v0_norm_sq
    h0_rhs = prefactor * chi_sq_over_4
    h1_rhs = prefactor * chi_sq_over_2
    point_rhs = [Fraction(2, 3) * chi_over_2 * v for v in v0]

    return {
        "label": "exact",
        "band": band,
        "v0": list(v0),
        "band_points": len(band_waves(band)),
        "active_points": len(table),
        "H0": _rational_record(h0),
        "H0_float": float(h0),
        "H1": _rational_record(h1),
        "H1_float": float(h1),
        "H0_law_rhs": _rational_record(h0_rhs),
        "H1_law_rhs": _rational_record(h1_rhs),
        "H0_residual": _rational_record(h0 - h0_rhs),
        "H1_residual": _rational_record(h1 - h1_rhs),
        "H0_law_exact": bool(h0 == h0_rhs),
        "H1_law_exact": bool(h1 == h1_rhs),
        "u_at_origin": [_rational_record(value) for value in point_value],
        "u_at_origin_float": [float(value) for value in point_value],
        "u_at_origin_law_rhs": [_rational_record(value) for value in point_rhs],
        "u_at_origin_residual": [
            _rational_record(value - rhs) for value, rhs in zip(point_value, point_rhs)
        ],
        "u_at_origin_law_exact": bool(
            all(value == rhs for value, rhs in zip(point_value, point_rhs))
        ),
        "real_cosine_field": bool(hermitian),
        "divergence_free": bool(divergence_free),
        "N0_sq": _rational_record(h1 / h0) if h0 else "0",
        "N0_sq_float": float(h1 / h0) if h0 else 0.0,
        "S_chi": _rational_record(chi_sq_over_2),
        "T_chi": _rational_record(chi_sq_over_4),
    }


def exact_nonlinear_mode(
    band: int, wave: Sequence[int], v0: Sequence[int] = DEFAULT_V0
) -> dict[str, object]:
    r"""One exact Fourier coefficient of the nonlinearity.

    Returns ``b_k=\sum_{l+m=k}(\hat u(l)\cdot m)\hat u(m)`` and
    ``P_kb_k`` as exact rationals, where
    ``\widehat{u_N\cdot\nabla u_N}(k)=i\,b_k`` and
    ``\widehat{\mathbb P(u_N\cdot\nabla u_N)}(k)=i\,P_kb_k``.

    Cost is ``O(|B_N|)``, one pass over the band, so this reaches bands far
    beyond the ``O(|B_N|^2)`` full contraction.  Since
    ``\hat u_N(k)=N^{-2}F(k/N)`` the sum is exactly a Riemann sum for
    ``C(k/N)=\int(F(k/N-b)\cdot b)F(b)\,db``, so ``b_k`` converges to
    ``-i\widehat{V\cdot\nabla V}(k/N)`` with **no** rescaling.
    """

    target = tuple(int(c) for c in wave)
    norm_sq = sum(c * c for c in target)
    if norm_sq == 0:
        raise ValueError("the mean mode carries no information here")
    scale, table = exact_scaled_coefficients(band, v0)

    accumulator = [0, 0, 0]
    for right, right_vector in table.items():
        left = (
            target[0] - right[0],
            target[1] - right[1],
            target[2] - right[2],
        )
        left_vector = table.get(left)
        if left_vector is None:
            continue
        weight = (
            left_vector[0] * right[0]
            + left_vector[1] * right[1]
            + left_vector[2] * right[2]
        )
        if weight == 0:
            continue
        accumulator[0] += weight * right_vector[0]
        accumulator[1] += weight * right_vector[1]
        accumulator[2] += weight * right_vector[2]

    scale_sq = Fraction(scale) ** 2
    b = [Fraction(component) / scale_sq for component in accumulator]
    dot = sum(Fraction(c) * value for c, value in zip(target, b))
    projected = [b[i] - Fraction(target[i]) * dot / norm_sq for i in range(3)]
    projected_norm_sq = sum(value * value for value in projected)

    return {
        "label": "exact",
        "band": band,
        "wave": list(target),
        "zeta": [c / band for c in target],
        "b_k": [_rational_record(value) for value in b],
        "b_k_float": [float(value) for value in b],
        "leray_b_k": [_rational_record(value) for value in projected],
        "leray_b_k_float": [float(value) for value in projected],
        "leray_b_k_norm_sq": _rational_record(projected_norm_sq),
        "leray_b_k_norm": math.sqrt(float(projected_norm_sq)),
        "is_exactly_zero": bool(all(value == 0 for value in projected)),
    }


def exact_nonlinear_norm_sq(
    band: int, v0: Sequence[int] = DEFAULT_V0
) -> dict[str, object]:
    r"""The full exact ``\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2`` and ``K``.

    ``O(|B_N|^2)`` ordered pairs; usable only at small bands.  Values are
    exact rationals with no truncation and no sampling.
    """

    scale, table = exact_scaled_coefficients(band, v0)
    items = list(table.items())

    accumulator: dict[tuple[int, int, int], list[int]] = {}
    for left, left_vector in items:
        c0, c1, c2 = left_vector
        for right, right_vector in items:
            weight = c0 * right[0] + c1 * right[1] + c2 * right[2]
            if weight == 0:
                continue
            key = (left[0] + right[0], left[1] + right[1], left[2] + right[2])
            current = accumulator.get(key)
            if current is None:
                accumulator[key] = [
                    weight * right_vector[0],
                    weight * right_vector[1],
                    weight * right_vector[2],
                ]
            else:
                current[0] += weight * right_vector[0]
                current[1] += weight * right_vector[1]
                current[2] += weight * right_vector[2]

    nonlinear_sq = Fraction(0)
    for key, vector in accumulator.items():
        norm_sq = key[0] ** 2 + key[1] ** 2 + key[2] ** 2
        if norm_sq == 0:
            if any(component != 0 for component in vector):
                raise AssertionError("the mean mode of u.grad u is nonzero")
            continue
        dot = key[0] * vector[0] + key[1] * vector[1] + key[2] * vector[2]
        modulus_sq = vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
        nonlinear_sq += Fraction(norm_sq * modulus_sq - dot * dot, norm_sq)
    nonlinear_sq /= Fraction(scale) ** 4

    scale_sq = Fraction(scale) ** 2
    h0 = Fraction(0)
    h1 = Fraction(0)
    for key, vector in items:
        norm_sq = key[0] ** 2 + key[1] ** 2 + key[2] ** 2
        modulus_sq = Fraction(sum(c * c for c in vector)) / scale_sq
        h0 += modulus_sq
        h1 += norm_sq * modulus_sq

    capacity = nonlinear_sq / h1**2
    n0_sq = h1 / h0
    return {
        "label": "exact",
        "band": band,
        "band_points": len(band_waves(band)),
        "active_points": len(items),
        "sumset_modes": len(accumulator),
        "H0": _rational_record(h0),
        "H1": _rational_record(h1),
        "nonlinear_sq": _rational_record(nonlinear_sq),
        "nonlinear_sq_float": float(nonlinear_sq),
        "nonlinear_sq_over_band_cubed": float(nonlinear_sq) / float(band**3),
        "K": _rational_record(capacity),
        "K_float": float(capacity),
        "N0_sq": _rational_record(n0_sq),
        "N0_sq_float": float(n0_sq),
        "K_over_N0_sq_float": float(capacity / n0_sq),
    }


# --------------------------------------------------------------------------
# the float / FFT lane
# --------------------------------------------------------------------------


def dealias_grid(band: int) -> int:
    r"""The smallest even grid on which ``u\cdot\nabla u`` is alias free.

    ``u`` is supported in ``|k|_\infty\le N`` and the quadratic product in
    ``|k|_\infty\le2N``; the repo's dealiasing guard
    (:func:`ns_certificate_lab.leray_response_relay.leray_advection`) demands
    ``2N<\lfloor G/2\rfloor``, so ``G\ge4N+2``.  Nothing is truncated: the
    *whole* product spectrum is represented, which is what a norm needs.
    """

    return 4 * band + 2


def spectral_field(
    grid: int, band: int, v0: Sequence[int] = DEFAULT_V0
) -> np.ndarray:
    r"""The smooth family as a Fourier array in the repo's torus convention.

    ``u(x)=\sum_k\hat u_k e^{ik\cdot x}``; the array index ``k`` is taken
    modulo ``grid``.  Coefficients are real and placed at both ``\pm k``, so
    the field is a real cosine field by construction.
    """

    if grid <= 2 * band:
        raise ValueError("the grid must resolve the band")
    seed = np.asarray(v0, dtype=np.float64)
    field = np.zeros((3, grid, grid, grid), dtype=np.complex128)
    for wave in band_waves(band):
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        weight = float(chi_exact(norm_sq, band))
        if weight == 0.0:
            continue
        k = np.asarray(wave, dtype=np.float64)
        coefficient = weight * (seed - k * (k @ seed) / norm_sq) / norm_sq
        index = tuple(int(c) % grid for c in wave)
        for component in range(3):
            field[(component, *index)] += coefficient[component]
    return field


def capacity_metrics(field: np.ndarray, band: int) -> dict[str, float]:
    r"""``H_0``, ``H_1``, ``N_0^2``, ``\|\mathbb P(u\cdot\nabla u)\|_2^2``, ``K``.

    Conventions are the paper's: ``H_0=\|u\|_2^2=\sum_k|\hat u_k|^2``,
    ``H_1=\|\nabla u\|_2^2``, ``K=\|\mathbb P(u\cdot\nabla u)\|_2^2/H_1^2``,
    all with the normalised measure ``(2\pi)^{-3}dx``.  IEEE double.
    """

    from .leray_response_relay import gradient_l2_squared, leray_advection

    h0 = float(np.sum(np.abs(field) ** 2))
    h1 = gradient_l2_squared(field)
    nonlinear = leray_advection(field, field)
    nonlinear_sq = float(np.vdot(nonlinear, nonlinear).real)
    n0_sq = h1 / h0
    capacity = nonlinear_sq / h1**2
    return {
        "label": "float",
        "band": band,
        "grid": int(field.shape[1]),
        "H0": h0,
        "H1": h1,
        "N0_sq": n0_sq,
        "z": math.log(n0_sq),
        "nonlinear_sq": nonlinear_sq,
        "nonlinear_sq_over_band_cubed": nonlinear_sq / float(band**3),
        "K": capacity,
        "K_over_N0_sq": capacity / n0_sq,
    }


# --------------------------------------------------------------------------
# the continuum profile V
# --------------------------------------------------------------------------

_SERIES_CUTOFF = 0.1


def _spherical_j0(s: np.ndarray) -> np.ndarray:
    s = np.abs(np.asarray(s, dtype=np.float64))
    out = np.empty_like(s)
    large = s > 1.0e-3
    big = s[large]
    out[large] = np.sin(big) / big
    small = s[~large]
    out[~large] = 1.0 - small**2 / 6.0 + small**4 / 120.0 - small**6 / 5040.0
    return out


def _spherical_j1(s: np.ndarray) -> np.ndarray:
    s = np.asarray(s, dtype=np.float64)
    out = np.empty_like(s)
    large = np.abs(s) > 1.0e-3
    big = s[large]
    out[large] = (np.sin(big) - big * np.cos(big)) / big**2
    small = s[~large]
    out[~large] = small / 3.0 - small**3 / 30.0 + small**5 / 840.0


    return out


def _kernel_b(s: np.ndarray) -> np.ndarray:
    r"""``\mathcal B(s)=\int_0^1\mu^2\cos(s\mu)d\mu``."""

    s = np.abs(np.asarray(s, dtype=np.float64))
    out = np.empty_like(s)
    large = s > _SERIES_CUTOFF
    big = s[large]
    out[large] = (
        np.sin(big) / big + 2.0 * np.cos(big) / big**2 - 2.0 * np.sin(big) / big**3
    )
    t = s[~large]
    out[~large] = (
        1.0 / 3.0
        - t**2 / 10.0
        + t**4 / 168.0
        - t**6 / 6480.0
        + t**8 / 443520.0
    )
    return out


def _kernel_b1(s: np.ndarray) -> np.ndarray:
    r"""``\mathcal B'(s)=-\int_0^1\mu^3\sin(s\mu)d\mu``."""

    s = np.asarray(s, dtype=np.float64)
    out = np.empty_like(s)
    large = np.abs(s) > _SERIES_CUTOFF
    big = s[large]
    out[large] = (
        np.cos(big) / big
        - 3.0 * np.sin(big) / big**2
        - 6.0 * np.cos(big) / big**3
        + 6.0 * np.sin(big) / big**4
    )
    t = s[~large]
    out[~large] = -t / 5.0 + t**3 / 42.0 - t**5 / 1080.0 + 8.0 * t**7 / 443520.0
    return out


def _kernel_b2(s: np.ndarray) -> np.ndarray:
    r"""``\mathcal B''(s)=-\int_0^1\mu^4\cos(s\mu)d\mu``."""

    s = np.abs(np.asarray(s, dtype=np.float64))
    out = np.empty_like(s)
    large = s > _SERIES_CUTOFF
    big = s[large]
    out[large] = -(
        np.sin(big) / big
        + 4.0 * np.cos(big) / big**2
        - 12.0 * np.sin(big) / big**3
        - 24.0 * np.cos(big) / big**4
        + 24.0 * np.sin(big) / big**5
    )
    t = s[~large]
    out[~large] = (
        -1.0 / 5.0 + 3.0 * t**2 / 42.0 - 5.0 * t**4 / 1080.0 + 56.0 * t**6 / 443520.0
    )
    return out


def _gauss_legendre_split(order: int) -> tuple[np.ndarray, np.ndarray]:
    r"""Nodes/weights on ``[0,1]``, split at ``\rho=1/2``.

    ``\chi`` is analytic on each half but only ``C^4`` across ``\rho=1/2``
    and ``\rho=1``; splitting at the interior break point restores spectral
    accuracy on both pieces.
    """

    nodes, weights = np.polynomial.legendre.leggauss(order)
    lower = 0.25 * nodes + 0.25, 0.25 * weights
    upper = 0.25 * nodes + 0.75, 0.25 * weights
    return (
        np.concatenate((lower[0], upper[0])),
        np.concatenate((lower[1], upper[1])),
    )


def radial_profiles(radius: float, order: int = 300) -> dict[str, float]:
    r"""``g,h`` and their first two ``r``-derivatives, by 1-D quadrature.

    ``g(r)=4\pi\int_0^1\chi j_0(\rho r)d\rho``,
    ``h(r)=4\pi\int_0^1\chi\mathcal B(\rho r)d\rho``, differentiated under
    the integral using ``j_0'=-j_1``, ``j_1'=\mathcal B``.  Integrands are
    smooth on each half of ``[0,1]``, so Gauss–Legendre converges
    spectrally; ``order=300`` per half is machine precision for
    ``r\lesssim100``.
    """

    if radius <= 0.0:
        raise ValueError("the radius must be positive")
    rho, weight = _gauss_legendre_split(order)
    chi = chi_float(rho)
    s = rho * radius
    four_pi = 4.0 * math.pi
    return {
        "r": float(radius),
        "g": four_pi * float(np.sum(weight * chi * _spherical_j0(s))),
        "g1": -four_pi * float(np.sum(weight * chi * rho * _spherical_j1(s))),
        "g2": -four_pi * float(np.sum(weight * chi * rho**2 * _kernel_b(s))),
        "h": four_pi * float(np.sum(weight * chi * _kernel_b(s))),
        "h1": four_pi * float(np.sum(weight * chi * rho * _kernel_b1(s))),
        "h2": four_pi * float(np.sum(weight * chi * rho**2 * _kernel_b2(s))),
    }


def vector_profile(radius: float, order: int = 300) -> dict[str, float]:
    r"""``V(y)=A(r)v_0+B(r)\hat y(\hat y\cdot v_0)`` for a unit seed.

    Also reports ``rA`` and ``rB``, which must both tend to ``\pi^2``
    (the closed-form untruncated profile
    ``V_\infty=\pi^2(v_0+\hat y(\hat y\cdot v_0))/|y|``).
    """

    p = radial_profiles(radius, order)
    a = 0.5 * (p["g"] + p["h"])
    b = 0.5 * (p["g"] - 3.0 * p["h"])
    return {
        "r": float(radius),
        "A": a,
        "B": b,
        "r_times_A": radius * a,
        "r_times_B": radius * b,
        "pi_squared": math.pi**2,
    }


def meridional_profiles(radius: float, order: int = 300) -> dict[str, float]:
    r"""The axisymmetric scalars ``A,\alpha,W,Z`` at one radius (unit seed).

    ``V=\alpha\cos\theta\,\hat e_r-A\sin\theta\,\hat e_\theta``,
    ``\omega=W\sin\theta\,\hat e_\varphi``,
    ``\nabla\times(V\cdot\nabla V)=Z\sin\theta\cos\theta\,\hat e_\varphi``.
    ``divergence_residual`` is ``r\alpha'+2\alpha-2A``, which the exact
    solenoidality of ``V`` forces to vanish; it is reported as a free
    consistency check on the whole radial construction.
    """

    p = radial_profiles(radius, order)
    a = 0.5 * (p["g"] + p["h"])
    a1 = 0.5 * (p["g1"] + p["h1"])
    a2 = 0.5 * (p["g2"] + p["h2"])
    alpha = p["g"] - p["h"]
    alpha1 = p["g1"] - p["h1"]
    r = float(radius)

    w = (alpha - a) / r - a1
    w1 = (alpha1 - a1) / r - (alpha - a) / r**2 - a2
    derivative = w * alpha + r * w1 * alpha + r * w * alpha1
    z = (derivative - 2.0 * w * a) / r
    return {
        "r": r,
        "A": a,
        "alpha": alpha,
        "W": w,
        "Z": z,
        "divergence_residual": r * alpha1 + 2.0 * alpha - 2.0 * a,
    }


def curl_coefficient(
    radius: float, v0: Sequence[float] = DEFAULT_V0, order: int = 300
) -> float:
    r"""``Z(r)``: the amplitude of ``\nabla\times(V\cdot\nabla V)``.

    ``\nabla\times(V\cdot\nabla V)=Z(r)\sin\theta\cos\theta\,\hat e_\varphi``
    with ``\theta`` measured from ``v_0``.  ``Z\ne0`` at a single radius
    already proves ``\mathbb P(V\cdot\nabla V)\not\equiv0``, because
    ``\mathbb PW\equiv0`` iff ``\nabla\times W\equiv0`` and the curl is a
    local operator.
    """

    seed_norm_sq = float(sum(float(c) * float(c) for c in v0))
    return seed_norm_sq * meridional_profiles(radius, order)["Z"]


def asymptotic_curl_coefficient(
    radius: float, v0: Sequence[float] = DEFAULT_V0
) -> float:
    r"""``-12\pi^4\|v_0\|^2/r^4``, the closed-form untruncated ``Z_\infty``."""

    seed_norm_sq = float(sum(float(c) * float(c) for c in v0))
    return -12.0 * math.pi**4 * seed_norm_sq / float(radius) ** 4


# --------------------------------------------------------------------------
# the continuum convolution coefficient C(zeta)
# --------------------------------------------------------------------------


def _symbol(xi: np.ndarray, v0: np.ndarray) -> np.ndarray:
    r"""``F(\xi)=\chi(|\xi|)P_\xi v_0/|\xi|^2``, zero outside the unit ball."""

    norm_sq = np.sum(xi * xi, axis=-1)
    out = np.zeros(xi.shape, dtype=np.float64)
    inside = (norm_sq > 0.0) & (norm_sq < 1.0)
    if not np.any(inside):
        return out
    x = xi[inside]
    n = norm_sq[inside]
    weight = chi_float(np.sqrt(n))
    out[inside] = weight[:, None] * (v0[None, :] - x * (x @ v0)[:, None] / n[:, None]) / n[:, None]
    return out


def _partition(b: np.ndarray, zeta: np.ndarray) -> np.ndarray:
    r"""``\varphi=d_\zeta^4/(d_0^4+d_\zeta^4)``: ``1`` at ``b=0``, ``0`` at ``b=\zeta``.

    The fourth power is what makes ``\varphi F(\zeta-b)`` vanish like
    ``|\zeta-b|^2`` at the far singularity (and ``(1-\varphi)F(b)`` like
    ``|b|^2`` at the near one), so each half of the split integrand is
    ``C^1`` and the product quadrature converges quickly.
    """

    near = np.sum(b * b, axis=-1) ** 2
    far = np.sum((zeta - b) * (zeta - b), axis=-1) ** 2
    return far / (near + far)


def _sphere_rule(polar: int, azimuth: int) -> tuple[np.ndarray, np.ndarray]:
    mu, mu_weight = np.polynomial.legendre.leggauss(polar)
    phi = 2.0 * math.pi * np.arange(azimuth) / azimuth
    phi_weight = np.full(azimuth, 2.0 * math.pi / azimuth)
    grid_mu, grid_phi = np.meshgrid(mu, phi, indexing="ij")
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - grid_mu**2))
    directions = np.stack(
        (sin_theta * np.cos(grid_phi), sin_theta * np.sin(grid_phi), grid_mu),
        axis=-1,
    ).reshape(-1, 3)
    return directions, np.outer(mu_weight, phi_weight).ravel()


def convolution_coefficient(
    zeta: Sequence[float],
    v0: Sequence[float] = DEFAULT_V0,
    radial: int = 60,
    polar: int = 60,
    azimuth: int = 60,
) -> np.ndarray:
    r"""``C(\zeta)=\int_{\mathbb R^3}(F(\zeta-b)\cdot b)F(b)\,db`` by quadrature.

    ``\widehat{V\cdot\nabla V}(\zeta)=i\,C(\zeta)``, so
    ``\mathbb P(V\cdot\nabla V)\not\equiv0`` follows from
    ``P_\zeta C(\zeta)\ne0`` at one ``\zeta``.  The integrand has two
    integrable ``|\cdot|^{-2}`` singularities, at ``b=0`` and ``b=\zeta``.
    They are separated by the smooth partition :func:`_partition` and each
    half is then integrated in spherical coordinates *centred on its own
    singularity*, where the ``\rho^2`` Jacobian cancels the ``|\cdot|^{-2}``
    exactly.  IEEE double; convergence is monitored by varying the three
    node counts.
    """

    zeta = np.asarray(zeta, dtype=np.float64)
    seed = np.asarray(v0, dtype=np.float64)
    rho, rho_weight = np.polynomial.legendre.leggauss(radial)
    rho = 0.5 * (rho + 1.0)
    rho_weight = 0.5 * rho_weight
    directions, direction_weight = _sphere_rule(polar, azimuth)

    points = (rho[:, None, None] * directions[None, :, :]).reshape(-1, 3)
    weight = (rho_weight[:, None] * direction_weight[None, :]).ravel()
    # F(rho omega) rho^2 = chi(rho) P_omega v0  -- the singular factor is gone.
    transverse = seed[None, :] - directions * (directions @ seed)[:, None]
    regular = (
        chi_float(rho)[:, None, None] * transverse[None, :, :]
    ).reshape(-1, 3)

    # near half: singularity at b = 0 only (the partition kills the far one)
    far_symbol = _symbol(zeta[None, :] - points, seed)
    near_weight = np.sum(far_symbol * points, axis=1) * _partition(points, zeta)
    near = np.sum((weight * near_weight)[:, None] * regular, axis=0)

    # far half: substitute b = zeta - e, singularity at e = 0 only
    shifted = zeta[None, :] - points
    far_weight = np.sum(regular * shifted, axis=1) * (
        1.0 - _partition(shifted, zeta)
    )
    far = np.sum((weight * far_weight)[:, None] * _symbol(shifted, seed), axis=0)

    return near + far
