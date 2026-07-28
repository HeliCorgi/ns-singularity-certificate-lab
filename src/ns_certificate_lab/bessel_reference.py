r"""Independent modified-Bessel reference oracle for the E-33 wall-response law.

Purpose
-------
Equation-audit entry **E-33** reduces the finite-cylinder problem
``-L5 psi = omega`` (E-25) to one modified Bessel equation of order one per
axial Fourier mode, and predicts that the wall-induced correction at a fixed
core radius is *exactly* proportional to

.. math::

   \rho(x) = \frac{K_1(x)}{I_1(x)}, \qquad x = kR .

Verifying that prediction numerically requires values of ``I_1`` and ``K_1``
that are computed **without** the finite-cylinder solvers, so that the
comparison is an oracle test rather than a tautology.  SciPy is not a
dependency of this project (see ``pyproject.toml``: NumPy only), so this module
implements the two functions from their defining series and integral, using
only :mod:`math` and :mod:`numpy`.

This module deliberately imports nothing from the rest of the package.  It
knows nothing about grids, solvers or experiments; it is a pure special-function
reference.

``I_nu`` by its ascending power series
--------------------------------------
For integer ``nu >= 0`` (DLMF 10.25.2)

.. math::

   I_\nu(x)=\sum_{m=0}^{\infty}
             \frac{(x/2)^{2m+\nu}}{m!\,(m+\nu)!},

which for ``nu = 1`` is the series named in the task specification,
``sum (x/2)^{2m+1} / (m! (m+1)!)``.  Every term is positive for ``x >= 0``, so
the sum has **no cancellation**: the relative error of the truncated sum is
bounded by the relative size of the discarded tail plus a few units in the last
place per accumulated term.  This is why the series is trustworthy here even
for ``x`` of order 60, where ``I_1(x) ~ 1e24``: the difficulty at large argument
is *overflow*, not accuracy.

Overflow is removed by folding ``exp(-x)`` into the first term and propagating
it through the recurrence

.. math::

   t_{m+1}=t_m\,\frac{(x/2)^2}{(m+1)(m+1+\nu)},

so that :func:`exp_scaled_i0` and :func:`exp_scaled_i1` accumulate
``exp(-x) I_nu(x)``, a quantity bounded by ``1`` for all ``x >= 0`` and
asymptotic to ``1/sqrt(2 pi x)``.  Nothing in the summation ever exceeds that
bound, so the scaled series is safe well past the range this project needs.

``K_nu`` by its integral representation
---------------------------------------
DLMF 10.32.9 gives, for ``x > 0``,

.. math::

   K_\nu(x)=\int_0^{\infty}e^{-x\cosh t}\cosh(\nu t)\,dt,

which for ``nu = 1`` is the representation named in the task specification.
The exponentially scaled form actually evaluated here is

.. math::

   e^{x}K_\nu(x)=\int_0^{\infty}
       e^{-x(\cosh t-1)}\cosh(\nu t)\,dt ,

whose integrand is bounded by ``cosh(nu t)`` near ``t = 0`` and decays
doubly exponentially, so the value stays ``O(1)`` (it tends to
``sqrt(pi/(2x))`` as ``x -> inf`` and grows only like ``1/x`` as ``x -> 0``).

**Quadrature choice, and why the error is not ``O(h^2)``.**  The integrand
``f(t) = exp(-x(cosh t - 1)) cosh(nu t)`` is *even* in ``t``, so

.. math::

   \int_0^\infty f = \tfrac12\int_{-\infty}^{\infty} f ,

and the uniform trapezoidal rule on the whole line,
``h * sum_{j in Z} f(jh)``, is the natural rule.  Halving it gives exactly the
one-sided rule used below, ``h * (f(0)/2 + sum_{j>=1} f(jh))``.  By Poisson
summation the whole-line trapezoidal error is
``sum_{n != 0} fhat(2 pi n / h)``, so it is governed by the decay of the
Fourier transform, i.e. by the strip of analyticity, and **not** by any Taylor
remainder.  Here ``cosh(u + iv)`` has real part ``cosh(u) cos(v)``, which is
positive and unbounded for ``|v| < pi/2``; hence ``f`` is analytic and
absolutely integrable on every strip ``|Im t| <= d < pi/2``, and the error is
bounded by ``C(x, d) exp(-2 pi d / h)``.  On the strip boundary the integrand is
bounded by ``exp(x)`` times an ``x``-independent factor, so the *relative* error
of the scaled integral is at worst of order ``exp(x) exp(-2 pi d / h)``.

With :data:`QUADRATURE_STEP` ``h = 1/32`` and ``d`` slightly below ``pi/2``
this bound is ``exp(x - 280)``, i.e. below ``1e-60`` for every argument this
project uses (``x <= 60``) and still below ``1e-16`` out to ``x ~ 240``.
:func:`k_quadrature_step_halving_defect` exposes the empirical check that this
bound is not vacuous; measured on this machine, halving ``h`` changes
``e^x K_1(x)`` by ``0`` to ``3e-15`` relative over ``x`` in ``[0.05, 60]``.

**Truncation.**  The tail is cut at ``t_max = arccosh(1 + L/x)`` with
:data:`QUADRATURE_EXPONENT_CUTOFF` ``L = 700``, so the first discarded sample
carries the factor ``exp(-700) ~ 1e-304``.  The discarded tail is bounded by
``(2/x) exp(-x cosh t_max) cosh(t_max) <= (2 L / x^2) exp(-L)``, which is below
``1e-290`` for every ``x >= 1e-3`` and is therefore never the limiting error.

The ratios
----------
:func:`k1_over_i1` never forms ``K_1`` or ``I_1`` separately.  It evaluates

.. math::

   \frac{K_1(x)}{I_1(x)}
   = e^{-2x}\,\frac{e^{x}K_1(x)}{e^{-x}I_1(x)} ,

in which both factors of the quotient are ``O(1)`` and only the explicit
``exp(-2x)`` is small.  For ``x = 60`` the result is ``~1e-53``, far above the
smallest normal double, while a naive ``bessel_k1(x)/bessel_i1(x)`` would divide
``1e-27`` by ``1e+24``; the scaled form is what keeps the ratio usable across
the whole range of ``kR`` that E-33 spans.

:func:`k0_over_k1` is the second ratio this project needs.  It is the only
special-function input to the exact modal transparent (Dirichlet-to-Neumann)
outer condition of ``docs/whole_space_transition.md`` (W-1),

.. math::

   \partial_r\hat\psi_k(R)
   +\Big[\frac2R+k\,\frac{K_0(kR)}{K_1(kR)}\Big]\hat\psi_k(R)=0 ,

and it is evaluated as the quotient of the two *exponentially scaled*
quadratures, ``(e^{x}K_0)/(e^{x}K_1)``.  The common factor ``e^{-x}`` cancels
identically, so no underflow can occur: at ``x = 700`` a naive
``bessel_k0(x)/bessel_k1(x)`` divides ``0.0`` by ``0.0``, while the scaled form
still returns a number within ``1e-3`` of its limit ``1``.  The quotient is
bounded in ``(0, 1)`` for every ``x > 0`` because ``K_0 < K_1`` pointwise
(DLMF 10.37.1), and it increases monotonically to ``1``; both facts are checked
by the tests, and the monotone bound is what makes the transparent bracket
``2/R + k K_0/K_1`` lie strictly between ``2/R`` and ``2/R + k``.

**Where the quadrature stops being the accurate branch.**  ``t_max =
arccosh(1 + L/x)`` shrinks like ``sqrt(2L/x)`` as ``x`` grows, so with a fixed
step the node count falls: ``x = 60`` gets 121 nodes, ``x = 1000`` gets 38 and
``x = 3000`` gets 22.  The Poisson-summation bound of the previous section is
``exp(x - 280)`` and is therefore vacuous well before that.  Measured on this
machine with :func:`k_quadrature_step_halving_defect`, ``e^x K_0`` is stable to
``6e-15`` at ``x = 240`` but only to ``3e-9`` at ``x = 1000``, ``3e-7`` at
``x = 1300`` and ``2e-3`` at ``x = 3000``.  The ratio is better conditioned
than either factor because the two quadratures share their nodes and their
errors partly cancel -- at ``x = 1300`` the ratio still matches the asymptotic
branch to ``4e-9`` -- but it degrades all the same.
:data:`LARGE_ARGUMENT_QUADRATURE_LIMIT` records the crossover, and consumers
that need ``K_0/K_1`` beyond it are expected to use the asymptotic branch and
to say so; :func:`ns_certificate_lab.transparent_boundary.outer_bracket` does
exactly that.  The functions themselves are left pure: they do not silently
switch formula behind the caller's back.

Its two asymptotic branches are exposed separately, because the transparent
condition's limits are exactly those branches:

* :func:`k0_over_k1_small_argument_asymptote` returns ``x(-log(x/2) - gamma)``,
  which follows from ``K_0(x) = -log(x/2) - gamma + O(x^2 log x)`` and
  ``K_1(x) = 1/x + O(x log x)`` (DLMF 10.30.2/10.30.3).  It vanishes as
  ``x -> 0``, so the transparent bracket tends to ``2/R``, i.e. to the exact
  condition satisfied by the ``k = 0`` decaying solution ``C/r^2``.
* :func:`k0_over_k1_large_argument_asymptote` returns
  ``1 - 1/(2x) + 3/(8x^2)``, obtained by dividing the standard expansions
  ``K_nu(x) ~ sqrt(pi/2x) e^{-x} (1 + (4 nu^2-1)/(8x)
  + (4 nu^2-1)(4 nu^2-9)/(128 x^2) + ...)`` (DLMF 10.40.2) for ``nu = 0`` and
  ``nu = 1`` and re-expanding.  Its limit ``1`` is the ``kR >> 1`` regime, where
  the bracket becomes ``2/R + k`` and the condition degenerates into pure
  outgoing exponential decay.

Validation
----------
``tests/test_wall_truncation_scaling.py`` checks this module against

* published values ``I_1(1) = 0.5651591039924850`` and
  ``K_1(1) = 0.6019072301972346`` (DLMF 10.25.2 / 10.31.1; Abramowitz & Stegun
  Table 9.8 lists the exponentially scaled forms ``e^{-1} I_1(1) ~ 0.20791042``
  and ``e^{1} K_1(1) ~ 1.6361535`` to fewer digits, so the test derives the
  scaled values from the two above rather than quoting a second table);
* two further points, ``K_1(1/2)`` and ``K_1(2)``, recomputed **inside the test
  module** from the independent ascending series of DLMF 10.31.1 rather than
  transcribed from a table, and ``I_1(2)``/``I_1(5)`` recomputed from the
  independent integral representation ``I_n(x) = (1/pi) int_0^pi e^{x cos t}
  cos(n t) dt`` (DLMF 10.32.3);
* the cross-product identity ``I_0(x) K_1(x) + I_1(x) K_0(x) = 1/x``
  (DLMF 10.28.2), which couples all four functions and is reproduced here to
  ``3e-16`` relative over ``x`` in ``[1e-2, 60]`` -- see
  :func:`wronskian_relative_defect`;
* the two asymptotic regimes of E-33(c).

``tests/test_transparent_boundary.py`` adds, for the ``K_0/K_1`` ratio,

* the published value ``K_0(1) = 0.4210244382407083`` (DLMF 10.31.1 /
  Abramowitz & Stegun Table 9.8), giving ``K_0(1)/K_1(1) = 0.699423...``, and a
  second point ``K_0(2)/K_1(2)`` recomputed **inside the test module** from the
  independent ascending series
  ``K_0(x) = -(log(x/2)+gamma) I_0(x) + sum_{m>=1} H_m (x^2/4)^m/(m!)^2``;
* the derivative identity ``K_1'(x) = -K_0(x) - K_1(x)/x`` (DLMF 10.29.2),
  evaluated by a fourth-order central difference of :func:`bessel_k1`, which is
  the identity the transparent condition is derived from and therefore the one
  that must hold if the derivation is to mean anything;
* monotonicity and the bound ``0 < K_0/K_1 < 1``;
* both asymptotic branches above.

``I_0`` was introduced only to make the cross-product identity available; the
``K_0`` branch is additionally the special-function content of the transparent
outer condition (W-1).
"""

from __future__ import annotations

import math

__all__ = [
    "EULER_MASCHERONI",
    "LARGE_ARGUMENT_QUADRATURE_LIMIT",
    "MINIMUM_K_ARGUMENT",
    "QUADRATURE_EXPONENT_CUTOFF",
    "QUADRATURE_STEP",
    "SERIES_RELATIVE_CUTOFF",
    "bessel_i0",
    "bessel_i1",
    "bessel_k0",
    "bessel_k1",
    "exp_scaled_i0",
    "exp_scaled_i1",
    "exp_scaled_k0",
    "exp_scaled_k1",
    "i1_over_argument",
    "k0_over_k1",
    "k0_over_k1_large_argument_asymptote",
    "k0_over_k1_small_argument_asymptote",
    "k1_over_i1",
    "k1_over_i1_large_argument_asymptote",
    "k1_over_i1_small_argument_asymptote",
    "k_quadrature_step_halving_defect",
    "wronskian_relative_defect",
]

QUADRATURE_STEP = 1.0 / 32.0
"""Trapezoidal step for the ``K_nu`` integral; see the module docstring."""

QUADRATURE_EXPONENT_CUTOFF = 700.0
"""Cut the ``K_nu`` integral where ``x (cosh t - 1)`` first exceeds this."""

SERIES_RELATIVE_CUTOFF = 1.0e-19
"""Stop the ``I_nu`` series once a term is this small relative to the sum."""

MINIMUM_K_ARGUMENT = 1.0e-10
"""Smallest argument accepted by the ``K_nu`` quadrature; see ``_exp_scaled_k``."""

EULER_MASCHERONI = 0.5772156649015328606065120900824
"""``gamma``, needed by the small-argument branch of ``K_0/K_1``."""

LARGE_ARGUMENT_QUADRATURE_LIMIT = 1.0e3
"""Argument beyond which the ``K_nu`` quadrature is no longer the better branch.

Below it the trapezoidal rule reproduces ``K_0/K_1`` to roundoff and the
three-term asymptote is the cruder of the two; above it the node count has
fallen far enough that the ordering reverses.  At exactly this argument the
asymptote's own remainder is about ``4e-10`` relative (it decays like
``0.375/x^3``) while the quadrature's step-halving defect is already ``3e-9``,
so switching here costs nothing and bounds the error for every larger argument.
See the module docstring for the measured table.
"""

_MAX_SERIES_TERMS = 4096
_MAX_QUADRATURE_NODES = 1_000_000


def _checked_argument(x: float, *, name: str, positive: bool) -> float:
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        raise TypeError(f"{name} must be a real number")
    value = float(x)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if positive:
        if value <= 0.0:
            raise ValueError(f"{name} must be strictly positive")
    elif value < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return value


def _exp_scaled_i(order: int, x: float) -> float:
    """Return ``exp(-x) * I_order(x)`` for ``order`` in ``{0, 1}``."""

    if x == 0.0:
        return 1.0 if order == 0 else 0.0
    half = 0.5 * x
    half_squared = half * half
    scale = math.exp(-x)
    term = scale if order == 0 else half * scale
    total = term
    for m in range(1, _MAX_SERIES_TERMS + 1):
        term *= half_squared / (m * (m + order))
        total += term
        # The series has strictly positive terms and the ratio
        # (x/2)^2 / (m (m+order)) falls below 1/2 once m exceeds x, so once a
        # term is negligible the whole remaining tail is too.
        if m > half and term <= SERIES_RELATIVE_CUTOFF * total:
            return total
    raise ArithmeticError(
        f"modified Bessel I_{order} series failed to converge at x={x!r}"
    )


def _exp_scaled_k(order: int, x: float, *, step: float) -> float:
    """Return ``exp(x) * K_order(x)`` for ``order`` in ``{0, 1}``."""

    if not step > 0.0:
        raise ValueError("quadrature step must be positive")
    if x < MINIMUM_K_ARGUMENT:
        # ``t_max`` grows like ``log(2 L / x)`` and ``cosh(t_max)`` like
        # ``L / x``; below this bound the truncation point would approach the
        # argument at which ``cosh`` itself overflows, so the quadrature would
        # silently stop meaning what its docstring says.  Callers that need
        # smaller arguments should use
        # :func:`k1_over_i1_small_argument_asymptote` and say so.
        raise ValueError(
            "x must be at least "
            f"{MINIMUM_K_ARGUMENT!r} for the K_nu quadrature; got {x!r}"
        )
    t_max = math.acosh(1.0 + QUADRATURE_EXPONENT_CUTOFF / x)
    node_count = int(math.ceil(t_max / step))
    if node_count > _MAX_QUADRATURE_NODES:
        raise ArithmeticError(f"K_{order} quadrature would need too many nodes")
    # One-sided form of the whole-line trapezoidal rule for the even integrand.
    total = 0.5  # the t=0 sample: exp(0) * cosh(0) = 1, halved
    for index in range(1, node_count + 1):
        t = index * step
        cosh_t = math.cosh(t)
        exponent = x * (cosh_t - 1.0)
        if exponent >= QUADRATURE_EXPONENT_CUTOFF:
            break
        total += math.exp(-exponent) * (cosh_t if order == 1 else 1.0)
    return step * total


def exp_scaled_i0(x: float) -> float:
    """Return ``exp(-x) * I_0(x)`` for ``x >= 0``."""

    return _exp_scaled_i(0, _checked_argument(x, name="x", positive=False))


def exp_scaled_i1(x: float) -> float:
    """Return ``exp(-x) * I_1(x)`` for ``x >= 0``."""

    return _exp_scaled_i(1, _checked_argument(x, name="x", positive=False))


def exp_scaled_k0(x: float, *, step: float = QUADRATURE_STEP) -> float:
    """Return ``exp(x) * K_0(x)`` for ``x > 0``."""

    return _exp_scaled_k(0, _checked_argument(x, name="x", positive=True), step=step)


def exp_scaled_k1(x: float, *, step: float = QUADRATURE_STEP) -> float:
    """Return ``exp(x) * K_1(x)`` for ``x > 0``."""

    return _exp_scaled_k(1, _checked_argument(x, name="x", positive=True), step=step)


def bessel_i0(x: float) -> float:
    """Return ``I_0(x)``.  Overflows to ``inf`` only beyond ``x ~ 713``."""

    return exp_scaled_i0(x) * math.exp(_checked_argument(x, name="x", positive=False))


def bessel_i1(x: float) -> float:
    """Return ``I_1(x)`` from the ascending power series."""

    return exp_scaled_i1(x) * math.exp(_checked_argument(x, name="x", positive=False))


def bessel_k0(x: float) -> float:
    """Return ``K_0(x)``.  Underflows to ``0`` only beyond ``x ~ 745``."""

    return exp_scaled_k0(x) * math.exp(-_checked_argument(x, name="x", positive=True))


def bessel_k1(x: float) -> float:
    """Return ``K_1(x)`` from the integral representation."""

    return exp_scaled_k1(x) * math.exp(-_checked_argument(x, name="x", positive=True))


def i1_over_argument(x: float) -> float:
    """Return ``I_1(x)/x``, which is finite at the origin with value ``1/2``.

    E-33(b) predicts that the wall-induced correction has the radial shape
    ``I_1(kr)/r = k * (I_1(kr)/(kr))``.  Writing it this way removes the
    removable singularity at the cylindrical axis, where the grid has a node,
    so the predicted shape can be evaluated at ``r = 0`` without a special
    case.  The series is the same one used by :func:`bessel_i1`, divided
    term by term by ``x``; it is *not* evaluated as a quotient.
    """

    value = _checked_argument(x, name="x", positive=False)
    half_squared = 0.25 * value * value
    term = 0.5
    total = term
    for m in range(1, _MAX_SERIES_TERMS + 1):
        term *= half_squared / (m * (m + 1))
        total += term
        if m > 0.5 * value and term <= SERIES_RELATIVE_CUTOFF * total:
            return total
    raise ArithmeticError(f"I_1(x)/x series failed to converge at x={value!r}")


def k1_over_i1(x: float) -> float:
    """Return ``K_1(x)/I_1(x)`` in the overflow-free scaled form.

    This is the exact ``R``-dependence of the E-33(b) wall correction.  See the
    module docstring for why the scaled quotient is used instead of dividing
    :func:`bessel_k1` by :func:`bessel_i1`.
    """

    value = _checked_argument(x, name="x", positive=True)
    return math.exp(-2.0 * value) * exp_scaled_k1(value) / exp_scaled_i1(value)


def k0_over_k1(x: float) -> float:
    """Return ``K_0(x)/K_1(x)``, the special-function part of the W-1 bracket.

    Both quadratures carry the same explicit factor ``e^{-x}``, so the ratio is
    formed from the exponentially scaled values and the factor cancels exactly.
    Nothing underflows: at ``x = 700`` the unscaled ``K_nu`` are both ``0.0`` in
    binary64, while this quotient is still accurate.  See the module docstring
    for the bounds ``0 < K_0/K_1 < 1`` and the monotonicity that the tests pin.
    """

    value = _checked_argument(x, name="x", positive=True)
    return exp_scaled_k0(value) / exp_scaled_k1(value)


def k0_over_k1_small_argument_asymptote(x: float) -> float:
    """Return ``x(-log(x/2) - gamma)``, the ``x << 1`` branch of ``K_0/K_1``.

    From ``K_0(x) = -log(x/2) - gamma + O(x^2 log x)`` and
    ``K_1(x) = 1/x + O(x log x)``.  It vanishes as ``x -> 0``, which is why the
    transparent bracket ``2/R + k K_0(kR)/K_1(kR)`` tends to ``2/R``.

    This branch is also the *evaluated* form below
    :data:`MINIMUM_K_ARGUMENT`, where the ``K_nu`` quadrature refuses to run;
    the transparent solver says so explicitly at its call site.
    """

    value = _checked_argument(x, name="x", positive=True)
    return value * (-math.log(0.5 * value) - EULER_MASCHERONI)


def k0_over_k1_large_argument_asymptote(x: float, *, terms: int = 3) -> float:
    """Return the ``x >> 1`` branch ``1 - 1/(2x) + 3/(8x^2)`` of ``K_0/K_1``.

    ``terms=1`` gives the bare limit ``1`` (the ``kR >> 1`` regime in which the
    transparent bracket degenerates to ``2/R + k``); ``terms=2`` adds
    ``-1/(2x)``; ``terms=3`` adds ``+3/(8x^2)``.  The coefficients come from
    dividing the DLMF 10.40.2 expansions of ``K_0`` and ``K_1`` and
    re-expanding, as recorded in the module docstring.
    """

    value = _checked_argument(x, name="x", positive=True)
    if terms not in (1, 2, 3):
        raise ValueError("terms must be 1, 2 or 3")
    series = 1.0
    if terms >= 2:
        series -= 0.5 / value
    if terms >= 3:
        series += 0.375 / (value * value)
    return series


def k1_over_i1_large_argument_asymptote(x: float, *, terms: int = 2) -> float:
    """Return the E-33(c) large-argument asymptote of ``K_1/I_1``.

    ``terms=1`` gives the leading form ``pi exp(-2x)``; ``terms=2`` adds the
    ``3/(4x)`` correction printed in E-33(c); ``terms=3`` adds the next
    coefficient ``9/(32 x^2)``, which follows from squaring the ``3/(8x)``
    term of the standard ``I_1``/``K_1`` expansions and is used by the tests to
    show that the E-33(c) remainder really is ``O(x^-2)``.
    """

    value = _checked_argument(x, name="x", positive=True)
    if terms not in (1, 2, 3):
        raise ValueError("terms must be 1, 2 or 3")
    series = 1.0
    if terms >= 2:
        series += 0.75 / value
    if terms >= 3:
        series += (9.0 / 32.0) / (value * value)
    return math.pi * math.exp(-2.0 * value) * series


def k1_over_i1_small_argument_asymptote(x: float) -> float:
    """Return the E-33(c) small-argument asymptote ``2/x^2`` of ``K_1/I_1``."""

    value = _checked_argument(x, name="x", positive=True)
    return 2.0 / (value * value)


def wronskian_relative_defect(x: float) -> float:
    """Return ``|x (I_0 K_1 + I_1 K_0) - 1|``, which must vanish.

    DLMF 10.28.2 states ``I_0(x) K_1(x) + I_1(x) K_0(x) = 1/x``.  Evaluated in
    the exponentially scaled variables the identity is
    ``x (S_{I0} S_{K1} + S_{I1} S_{K0}) = 1`` with ``S_{I} = e^{-x} I`` and
    ``S_{K} = e^{x} K``, so the check costs nothing and stays in range for every
    argument.  It ties the power series to the quadrature: a wrong factorial, a
    wrong quadrature weight or a wrong truncation in either one would break it.
    """

    value = _checked_argument(x, name="x", positive=True)
    product = (
        exp_scaled_i0(value) * exp_scaled_k1(value)
        + exp_scaled_i1(value) * exp_scaled_k0(value)
    )
    return abs(value * product - 1.0)


def k_quadrature_step_halving_defect(x: float, order: int = 1) -> float:
    """Return the relative change of ``e^x K_order(x)`` when ``h`` is halved.

    The module docstring argues that the trapezoidal error decays like
    ``exp(-2 pi d / h)``.  If that argument were wrong -- if the rule were only
    second order, say -- halving ``h`` would change the answer by a measurable
    amount.  This function makes the claim falsifiable at run time.
    """

    value = _checked_argument(x, name="x", positive=True)
    if order not in (0, 1):
        raise ValueError("order must be 0 or 1")
    coarse = _exp_scaled_k(order, value, step=QUADRATURE_STEP)
    fine = _exp_scaled_k(order, value, step=0.5 * QUADRATURE_STEP)
    return abs(coarse - fine) / abs(fine)
