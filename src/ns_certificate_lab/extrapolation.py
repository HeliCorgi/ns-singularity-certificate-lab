r"""Blind power-law extrapolation of a resolution series, with subset sensitivity.

Purpose
-------
A refinement series ``A(h_1) > A(h_2) > ...`` (or the reverse) is often
extrapolated to ``h -> 0`` through the one-term model

.. math::

   A(h) = A_\infty + C h^{p},

and the extrapolant ``A_inf`` is then quoted as "the" limit.  That step is only
legitimate when the series is already inside the asymptotic range of the model,
and the *only* honest way to test that from the series alone is to refit it on
different subsets and different windows and to look at how much the answer
moves.  This module does exactly that, and nothing else.

**Scope, stated honestly.**  This module

* fits ``A_inf + C h^p`` by three independent routes (an exact three-point
  solve, a nonlinear least-squares solve with a free exponent, and a
  fixed-exponent Richardson solve);
* reports the *spread* of the resulting extrapolants as its primary output;
* refuses to certify a limit when that spread exceeds preregistered
  thresholds.

It does **not**

* prove that ``A_inf + C h^p`` is the right model.  A series carrying a
  logarithm, two competing powers, or an unresolved feature can pass every
  test here and still extrapolate to the wrong number.  Passing the criteria
  below means "the extrapolants are mutually consistent", never "the limit is
  correct";
* know anything about any published value.  See the next section;
* replace a resolution gate.  Agreement of extrapolants across subsets is a
  necessary condition for quoting a limit, not a sufficient one.

Why there is no reference value anywhere in this API
----------------------------------------------------
Fitting a series with the intended answer held fixed -- solving only for ``C``
and ``p`` with ``A_inf`` pinned to an externally published number -- always
"works", because two free parameters can pass through the remaining points of
a short series regardless of whether the model applies.  Such a fit measures
nothing and yet reads as confirmation.  The project audit therefore forbids it.

A ban that lives only in prose is not a ban, so it is enforced structurally
here: **no function in this module accepts a reference, target, anchor or
expected value in any form.**  There is no keyword to pass one to, no attribute
to set, and :func:`blind_extrapolation_report` raises :class:`ValueError` if any
keyword argument at all is supplied, with a sharper message when the keyword
name is one of the anchoring words in :data:`ANCHORING_KEYWORD_STEMS`.  A fit
that used a published limit would have to be written somewhere else, where a
reviewer can see it.  For the same reason, no verdict produced here ever
compares ``A_inf`` against an external number: the verdict strings speak only
about internal consistency of the series with itself.

The three-point solve
---------------------
With three spacings ``h1 > h2 > h3`` and values ``A1, A2, A3``, subtracting the
model pairwise eliminates ``A_inf``,

.. math::

   A_1 - A_2 = C (h_1^p - h_2^p), \qquad A_2 - A_3 = C (h_2^p - h_3^p),

and dividing eliminates ``C`` as well, leaving one scalar equation for ``p``:

.. math::

   \frac{A_2 - A_3}{A_1 - A_2} = \frac{h_2^p - h_3^p}{h_1^p - h_2^p}
                               =: g(p).

Writing ``u = h_2/h_1`` and ``v = h_3/h_1`` (so ``0 < v < u < 1``) and dividing
through by ``h_1^p`` gives the scale-free form actually evaluated here,

.. math::

   g(p) = \frac{u^p - v^p}{1 - u^p},

which is what :func:`_exponent_equation_defect` computes, through ``expm1`` so
that neither difference is formed by cancelling two nearby numbers.  For an
equally-ratioed series ``h_k = h_1 r^{k-1}`` it collapses to ``g(p) = r^p``, so
the root is available in closed form and the numerical root here can be checked
against it.

Two facts about ``g`` shape the implementation.  First, ``g(p) > 0`` for every
``p > 0``, because ``0 < v < u < 1`` makes both the numerator and the
denominator positive.  A series whose two increments have opposite signs
therefore has **no** root at any exponent, and the correct answer is
``converged = False`` rather than a fitted number -- this is the normal outcome
for a non-monotone or noise-dominated series, not an error.  Second, ``g`` is
smooth and (for the spacing ratios that arise in refinement studies) monotone
decreasing, so a sign-change scan followed by bisection is a robust solver.
The scan is retained rather than assuming monotonicity: the number of sign
changes found is reported in the ``detail`` field, so a second root would be
visible instead of silently preferred.

Once ``p`` is known,

.. math::

   A_\infty = A_1 - \frac{A_1 - A_2}{1 - u^p}, \qquad
   C = \frac{A_1 - A_2}{h_1^p\,(1 - u^p)} ,

both of which are evaluated in the ``u``-scaled form above, so that ``A_inf``
never depends on the magnitude of ``h_1^p`` and stays accurate even when that
factor is tiny.

The least-squares solve
-----------------------
For ``n >= 3`` points the exponent is found by an outer scalar search with an
inner *linear* least squares for ``(A_inf, C)``, which is the standard
separable (variable-projection) treatment of this model: for fixed ``p`` the
model is linear in the two remaining parameters, so the inner problem has a
closed-form solution and the outer problem is one-dimensional.

The outer search is a dense scan on ``[0.05, 8]`` at spacing
:data:`COARSE_EXPONENT_STEP` followed by golden-section refinement inside the
bracketing interval.  A dense scan is used in preference to a bare
golden-section search over the whole interval because the residual as a
function of ``p`` need not be unimodal for noisy data, and a local search
started from an arbitrary point would then report whichever basin it happened
to land in.

The inner problem is posed in the scaled basis ``(h/h_max)^p``, whose entries
lie in ``(0, 1]`` for every admissible exponent.  In the raw basis ``h^p`` the
second column can be ``1e-20`` while the first is ``1``, and the normal
equations lose most of their significant digits at large ``p``; in the scaled
basis the design matrix is well conditioned for all ``p`` in range, and ``C``
is recovered afterwards by dividing out ``h_max^p``.  ``A_inf`` -- the quantity
this module exists to report -- is unaffected by that division.

If the best exponent on the coarse scan sits at an endpoint of ``[0.05, 8]``,
the data prefer an exponent outside the searched window and the fit is reported
as ``converged = False`` with no parameters, rather than as a fit pinned to the
window edge.

Preregistered acceptance thresholds
-----------------------------------
:data:`A_INF_SPREAD_TOLERANCE` ``= 0.05`` and :data:`P_SPREAD_TOLERANCE`
``= 0.2`` are fixed **here, in code, before application to any new series**,
together with the two minimum-count requirements
:data:`MINIMUM_CONVERGED_FITS` and :data:`MINIMUM_FREE_EXPONENT_FITS`.  They are
recorded in this module rather than chosen per experiment so that a series
cannot be declared asymptotic by widening a threshold after seeing the answer.
Changing either number is a change of the acceptance rule and must be reviewed
as such; the values are also echoed into every report under the ``thresholds``
key, so any saved evidence carries the rule that was applied to it.

:func:`subset_sensitivity` declares ``in_asymptotic_range = True`` only when
**both**

* ``A_inf_spread_relative = (max - min) / median(|A_inf|) <= 0.05`` over all
  converged fits, and
* ``p_max - p_min <= 0.2`` over the converged fits whose exponent is a *fitted*
  parameter,

and only when enough fits converged to make those two numbers mean anything.
The exponent range deliberately excludes the fixed-order fits: their exponent
was assumed, not measured, so including them would make ``p_max - p_min >= 1``
for every series and the criterion would be vacuous.  Their ``A_inf`` values
*are* included in the spread, because a fixed-order Richardson extrapolation is
a genuine alternative model of the same data and disagreement with it is
exactly the kind of model sensitivity this module is meant to expose.

**A ``False`` verdict means no limit estimate may be quoted from the series.**
Not the mean of the extrapolants, not the value from the finest window, not the
one that looks most plausible.  The series is a monotone approach that has not
entered the asymptotic range of the fitted model, and the reportable content is
the spread itself.  :func:`blind_extrapolation_report` returns exactly that
sentence in its ``verdict`` field.

Numerical conventions
---------------------
Inputs are validated and then sorted into order of *decreasing* ``h`` (coarsest
first), so "contiguous triple" and "contiguous window" are well defined however
the caller ordered the series.  Parameters of a fit that did not converge are
``None``, never ``NaN``: this project's evidence writers reject non-finite
numbers outright (``_integrity.require_finite_json``), so every dict returned
here is JSON-serializable as it stands, and a non-converged fit cannot be
mistaken for a number by arithmetic downstream.

Only :mod:`numpy` and :mod:`math` are used; SciPy is not a dependency of this
project.  Nothing in this module imports the rest of the package.
"""

from __future__ import annotations

import math
from typing import Any, Callable, Sequence

import numpy as np

__all__ = [
    "ANCHORING_KEYWORD_STEMS",
    "A_INF_SPREAD_TOLERANCE",
    "ASYMPTOTIC_RANGE_VERDICT",
    "COARSE_EXPONENT_STEP",
    "MAXIMUM_EXPONENT",
    "MINIMUM_CONVERGED_FITS",
    "MINIMUM_EXPONENT",
    "MINIMUM_FREE_EXPONENT_FITS",
    "NOT_ASYMPTOTIC_RANGE_VERDICT",
    "P_SPREAD_TOLERANCE",
    "blind_extrapolation_report",
    "fit_fixed_order",
    "fit_power_law_least_squares",
    "fit_power_law_three_point",
    "subset_sensitivity",
]

A_INF_SPREAD_TOLERANCE = 0.05
"""Preregistered bound on ``(max - min) / median(|A_inf|)`` across fits.

Fixed here, in code, before application to any new series.  Five percent is the
largest disagreement between independent extrapolants of the same data that is
still small enough for the common value to carry meaning at the precision these
studies quote.  Raising it is a change of the acceptance rule, not a tuning
knob; see the module docstring.
"""

P_SPREAD_TOLERANCE = 0.2
"""Preregistered bound on ``p_max - p_min`` across fits with a free exponent.

Fixed here, in code, before application to any new series.  The observed order
of a scheme is an integer or a half-integer in almost every case of interest,
so a window of ``0.2`` is tight enough to distinguish "all subsets see the same
order" from "the order drifts with resolution" while still tolerating the
scatter that finite-precision data produce.
"""

MINIMUM_EXPONENT = 0.05
"""Lower end of the searched exponent range.

Below this the model is nearly indistinguishable from a constant plus a
logarithm over any realistic range of ``h``, and ``C`` diverges.
"""

MAXIMUM_EXPONENT = 8.0
"""Upper end of the searched exponent range.

No discretization in this project is above eighth order; an apparent exponent
this large means the increments are dominated by something other than the
leading truncation term.
"""

COARSE_EXPONENT_STEP = 1.0e-3
"""Spacing of the dense scan that brackets the least-squares exponent."""

MINIMUM_CONVERGED_FITS = 3
"""Converged fits required before a spread is treated as informative.

Preregistered with the two tolerances.  A spread computed from one or two fits
is not evidence of agreement; with fewer than this many converged fits,
``in_asymptotic_range`` is ``False`` regardless of the numbers.
"""

MINIMUM_FREE_EXPONENT_FITS = 2
"""Converged free-exponent fits required before the exponent range is used.

Preregistered with the two tolerances.  One fitted exponent has a range of
zero, which would satisfy :data:`P_SPREAD_TOLERANCE` trivially.
"""

NOT_ASYMPTOTIC_RANGE_VERDICT = (
    "not_in_asymptotic_range: no limit estimate is quotable"
)
"""Verdict when the extrapolants disagree by more than the preregistered bounds."""

ASYMPTOTIC_RANGE_VERDICT = (
    "asymptotic_range_criteria_met: "
    "extrapolants agree within preregistered spread"
)
"""Verdict when every converged extrapolant agrees within the preregistered bounds.

It asserts internal consistency of the series with itself and nothing more.  It
is not a statement that the limit is correct, and it is never a comparison
against any external value.
"""

ANCHORING_KEYWORD_STEMS = frozenset(
    {
        "anchor",
        "expected",
        "gold",
        "known",
        "limit",
        "published",
        "reference",
        "target",
        "true",
        "truth",
    }
)
"""Keyword stems that would indicate an attempt to anchor a fit.

:func:`blind_extrapolation_report` rejects *every* extra keyword argument; these
stems only sharpen the error message, so the refusal cannot be sidestepped by
inventing a new spelling.
"""

_GOLDEN_SECTION_RATIO = 0.5 * (math.sqrt(5.0) - 1.0)
_GOLDEN_SECTION_ITERATIONS = 200
_GOLDEN_SECTION_TOLERANCE = 1.0e-12
_ROOT_SCAN_POINTS = 4001
_BISECTION_ITERATIONS = 200
_BISECTION_TOLERANCE = 1.0e-14


def _validated_series(
    h: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate a spacing/value series and return it sorted by decreasing ``h``.

    Sorting is part of the contract: "contiguous triple" and "contiguous
    window" in :func:`subset_sensitivity` refer to this order, so the set of
    subsets a series produces does not depend on how the caller happened to
    list it.
    """

    try:
        spacings = np.asarray(h, dtype=np.float64)
        magnitudes = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("h and values must be real numeric sequences") from exc

    if spacings.ndim != 1 or magnitudes.ndim != 1:
        raise ValueError("h and values must be one-dimensional")
    if spacings.size != magnitudes.size:
        raise ValueError(
            "h and values must have the same length; got "
            f"{spacings.size} and {magnitudes.size}"
        )
    if spacings.size < 3:
        raise ValueError(
            "at least three points are required to fit A_inf + C h^p; got "
            f"{spacings.size}"
        )
    if not np.all(np.isfinite(spacings)):
        raise ValueError("h must be finite")
    if not np.all(spacings > 0.0):
        raise ValueError("h must be strictly positive")
    if not np.all(np.isfinite(magnitudes)):
        raise ValueError("values must be finite")
    if np.unique(spacings).size != spacings.size:
        raise ValueError("h must contain no duplicated spacings")

    order = np.argsort(-spacings, kind="stable")
    return spacings[order], magnitudes[order]


def _failed_fit(model: str, n_points: int, detail: str, **extra: Any) -> dict[str, Any]:
    """Return the canonical record of a fit that produced no parameters."""

    record: dict[str, Any] = {
        "model": model,
        "n_points": int(n_points),
        "p_is_free": model != "fixed_order",
        "A_inf": None,
        "C": None,
        "p": None,
        "residual_norm": None,
        "converged": False,
        "detail": detail,
    }
    record.update(extra)
    return record


def _exponent_equation_defect(
    exponent: float, ratio_middle: float, ratio_fine: float, target: float
) -> float:
    """Return ``g(p) - target`` for the three-point exponent equation.

    ``g(p) = (u^p - v^p) / (1 - u^p)`` with ``u = h2/h1`` and ``v = h3/h1``.
    Both differences are formed with ``expm1`` rather than by subtracting two
    nearby powers, so the defect keeps its significant digits when ``p`` is
    small and the ratios are close to one.
    """

    log_middle = math.log(ratio_middle)
    log_fine = math.log(ratio_fine)
    # -expm1(x) > 0 for x < 0, and both exponents below are strictly negative.
    numerator = math.exp(exponent * log_middle) * (
        -math.expm1(exponent * (log_fine - log_middle))
    )
    denominator = -math.expm1(exponent * log_middle)
    return numerator / denominator - target


def _bracketed_roots(
    defect: Callable[[float], float],
) -> list[tuple[float, float]]:
    """Return every sign-change bracket of ``defect`` on the exponent range."""

    grid = np.linspace(MINIMUM_EXPONENT, MAXIMUM_EXPONENT, _ROOT_SCAN_POINTS)
    values = np.array([defect(float(point)) for point in grid])
    if not np.all(np.isfinite(values)):
        return []
    brackets: list[tuple[float, float]] = []
    for index in range(values.size - 1):
        left, right = float(values[index]), float(values[index + 1])
        if left == 0.0 or left * right < 0.0:
            brackets.append((float(grid[index]), float(grid[index + 1])))
    return brackets


def _bisect(defect: Callable[[float], float], lower: float, upper: float) -> float:
    """Bisect a sign-change bracket of ``defect`` to machine tolerance."""

    left, right = lower, upper
    left_value = defect(left)
    if left_value == 0.0:
        return left
    for _ in range(_BISECTION_ITERATIONS):
        middle = 0.5 * (left + right)
        if right - left <= _BISECTION_TOLERANCE * max(1.0, abs(middle)):
            break
        middle_value = defect(middle)
        if middle_value == 0.0:
            return middle
        if left_value * middle_value < 0.0:
            right = middle
        else:
            left, left_value = middle, middle_value
    return 0.5 * (left + right)


def fit_power_law_three_point(
    h: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
) -> dict[str, Any]:
    r"""Solve ``A(h) = A_inf + C h^p`` exactly on three points.

    Exactly three spacings and three values must be supplied.  ``A_inf`` and
    ``C`` are eliminated analytically and the remaining scalar equation for
    ``p`` is solved by a sign-change scan followed by bisection on
    ``[MINIMUM_EXPONENT, MAXIMUM_EXPONENT]``; see the module docstring for the
    algebra and for why a scan is used instead of assuming monotonicity.

    Returns a dict with ``"A_inf"``, ``"C"``, ``"p"`` and ``"converged"``,
    together with ``"model"``, ``"n_points"``, ``"p_is_free"``,
    ``"residual_norm"`` (the defect of the third equation, which the
    elimination does not impose and which is therefore an independent check on
    the solve) and ``"detail"``.  When ``"converged"`` is ``False`` the three
    parameters are ``None``.

    Non-convergence is a legitimate result, not an error.  ``g(p) > 0`` for
    every admissible exponent, so a series whose two increments have opposite
    signs -- a non-monotone or noise-dominated series -- has no root, and this
    function says so.
    """

    spacings, magnitudes = _validated_series(h, values)
    if spacings.size != 3:
        raise ValueError(
            f"fit_power_law_three_point requires exactly three points; got "
            f"{spacings.size}"
        )

    coarse, middle, fine = (float(value) for value in spacings)
    first, second, third = (float(value) for value in magnitudes)

    leading_increment = first - second
    if leading_increment == 0.0:
        return _failed_fit(
            "three_point",
            3,
            "degenerate: the two coarsest values coincide, so C is undetermined",
        )
    target = (second - third) / leading_increment

    ratio_middle = middle / coarse
    ratio_fine = fine / coarse

    def defect(exponent: float) -> float:
        return _exponent_equation_defect(exponent, ratio_middle, ratio_fine, target)

    brackets = _bracketed_roots(defect)
    if not brackets:
        return _failed_fit(
            "three_point",
            3,
            "no exponent in "
            f"[{MINIMUM_EXPONENT}, {MAXIMUM_EXPONENT}] reproduces the observed "
            f"increment ratio {target!r}",
        )

    exponent = _bisect(defect, *brackets[0])
    complement = -math.expm1(exponent * math.log(ratio_middle))  # 1 - (h2/h1)^p
    if not complement > 0.0:
        return _failed_fit(
            "three_point",
            3,
            "degenerate: 1 - (h2/h1)^p underflowed to zero",
        )

    limit = first - leading_increment / complement
    coarse_power = coarse**exponent
    coefficient = (
        leading_increment / (coarse_power * complement)
        if coarse_power > 0.0
        else math.inf
    )
    # The third equation was eliminated, so its defect is an independent check.
    residual = abs(
        limit
        + leading_increment * (ratio_fine**exponent) / complement
        - third
    )

    if not all(math.isfinite(value) for value in (limit, coefficient, exponent)):
        return _failed_fit(
            "three_point",
            3,
            "solve produced a non-finite parameter",
        )

    return {
        "model": "three_point",
        "n_points": 3,
        "p_is_free": True,
        "A_inf": float(limit),
        "C": float(coefficient),
        "p": float(exponent),
        "residual_norm": float(residual),
        "converged": True,
        "detail": f"bisection on {len(brackets)} sign-change bracket(s)",
    }


def _linear_fit_at_exponent(
    ratios: np.ndarray, magnitudes: np.ndarray, exponent: float
) -> tuple[float, float, float] | None:
    """Least-squares ``(A_inf, slope, residual_norm)`` in the scaled basis.

    ``ratios`` are ``h / h_max`` and ``slope`` multiplies ``ratios**exponent``,
    so the caller recovers ``C = slope / h_max**exponent``.  Returns ``None``
    when the basis column is constant, which makes the two columns collinear.
    """

    basis = ratios**exponent
    basis_mean = float(basis.mean())
    value_mean = float(magnitudes.mean())
    centered_basis = basis - basis_mean
    denominator = float(centered_basis @ centered_basis)
    if not denominator > 0.0 or not math.isfinite(denominator):
        return None
    slope = float(centered_basis @ (magnitudes - value_mean)) / denominator
    intercept = value_mean - slope * basis_mean
    residual = magnitudes - (intercept + slope * basis)
    return intercept, slope, float(np.linalg.norm(residual))


def _residual_norms_on_grid(
    ratios: np.ndarray, magnitudes: np.ndarray, exponents: np.ndarray
) -> np.ndarray:
    """Vectorized :func:`_linear_fit_at_exponent` residual norms over exponents.

    The residuals are formed explicitly rather than through the algebraic
    identity ``||r||^2 = S_yy - S_sy^2 / S_ss``: near the optimum that identity
    is a difference of two nearly equal numbers and loses every significant
    digit exactly where the scan has to be sharpest.
    """

    basis = ratios[None, :] ** exponents[:, None]
    basis_mean = basis.mean(axis=1, keepdims=True)
    centered_basis = basis - basis_mean
    value_mean = float(magnitudes.mean())
    denominator = np.einsum("ij,ij->i", centered_basis, centered_basis)
    usable = denominator > 0.0
    safe_denominator = np.where(usable, denominator, 1.0)
    slope = centered_basis @ (magnitudes - value_mean) / safe_denominator
    slope = np.where(usable, slope, 0.0)
    intercept = value_mean - slope * basis_mean[:, 0]
    residual = magnitudes[None, :] - (intercept[:, None] + slope[:, None] * basis)
    norms = np.linalg.norm(residual, axis=1)
    return np.where(usable, norms, np.inf)


def _golden_section_minimum(
    objective: Callable[[float], float], lower: float, upper: float
) -> float:
    """Minimize a scalar objective on ``[lower, upper]`` by golden section."""

    left, right = lower, upper
    inner_left = right - _GOLDEN_SECTION_RATIO * (right - left)
    inner_right = left + _GOLDEN_SECTION_RATIO * (right - left)
    value_left = objective(inner_left)
    value_right = objective(inner_right)
    for _ in range(_GOLDEN_SECTION_ITERATIONS):
        if right - left <= _GOLDEN_SECTION_TOLERANCE:
            break
        if value_left <= value_right:
            right, inner_right, value_right = inner_right, inner_left, value_left
            inner_left = right - _GOLDEN_SECTION_RATIO * (right - left)
            value_left = objective(inner_left)
        else:
            left, inner_left, value_left = inner_left, inner_right, value_right
            inner_right = left + _GOLDEN_SECTION_RATIO * (right - left)
            value_right = objective(inner_right)
    return 0.5 * (left + right)


def fit_power_law_least_squares(
    h: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
) -> dict[str, Any]:
    r"""Fit ``A(h) = A_inf + C h^p`` by separable least squares on ``n >= 3`` points.

    The exponent is found by a dense scan on
    ``[MINIMUM_EXPONENT, MAXIMUM_EXPONENT]`` at spacing
    :data:`COARSE_EXPONENT_STEP` followed by golden-section refinement inside
    the bracketing interval; for each candidate exponent the pair
    ``(A_inf, C)`` is obtained in closed form, since the model is linear in
    those two parameters once ``p`` is fixed.

    On exactly three points this reproduces :func:`fit_power_law_three_point`,
    because two linear parameters and one exponent fit three points exactly and
    the least-squares residual at the solution is zero.  That coincidence is a
    useful cross-check of the two independent implementations and is asserted
    in the tests.

    Returns ``"A_inf"``, ``"C"``, ``"p"``, ``"residual_norm"`` and
    ``"converged"``, plus ``"model"``, ``"n_points"``, ``"p_is_free"``,
    ``"exponent_at_search_boundary"`` and ``"detail"``.  A best exponent at
    either end of the searched range is reported as non-converged with ``None``
    parameters: the data prefer an exponent the preregistered window does not
    contain, and a value pinned to the edge would misrepresent that.
    """

    spacings, magnitudes = _validated_series(h, values)
    scale = float(spacings[0])
    ratios = spacings / scale

    point_count = int(
        round((MAXIMUM_EXPONENT - MINIMUM_EXPONENT) / COARSE_EXPONENT_STEP) + 1
    )
    exponents = np.linspace(MINIMUM_EXPONENT, MAXIMUM_EXPONENT, point_count)
    norms = _residual_norms_on_grid(ratios, magnitudes, exponents)
    index = int(np.argmin(norms))
    if not math.isfinite(float(norms[index])):
        return _failed_fit(
            "least_squares",
            spacings.size,
            "no exponent yields a solvable linear subproblem",
            exponent_at_search_boundary=False,
        )
    if index in (0, exponents.size - 1):
        return _failed_fit(
            "least_squares",
            spacings.size,
            "best exponent lies at the edge of "
            f"[{MINIMUM_EXPONENT}, {MAXIMUM_EXPONENT}]; the data prefer an "
            "exponent outside the preregistered search window",
            exponent_at_search_boundary=True,
        )

    lower = float(exponents[index - 1])
    upper = float(exponents[index + 1])

    def objective(exponent: float) -> float:
        fit = _linear_fit_at_exponent(ratios, magnitudes, exponent)
        return math.inf if fit is None else fit[2]

    exponent = _golden_section_minimum(objective, lower, upper)
    refined = _linear_fit_at_exponent(ratios, magnitudes, exponent)
    if refined is None:
        return _failed_fit(
            "least_squares",
            spacings.size,
            "refined exponent yields a collinear design matrix",
            exponent_at_search_boundary=False,
        )

    limit, slope, residual = refined
    scale_power = scale**exponent
    coefficient = slope / scale_power if scale_power > 0.0 else math.inf
    if not all(math.isfinite(value) for value in (limit, coefficient, exponent)):
        return _failed_fit(
            "least_squares",
            spacings.size,
            "fit produced a non-finite parameter",
            exponent_at_search_boundary=False,
        )

    return {
        "model": "least_squares",
        "n_points": int(spacings.size),
        "p_is_free": True,
        "A_inf": float(limit),
        "C": float(coefficient),
        "p": float(exponent),
        "residual_norm": float(residual),
        "converged": True,
        "exponent_at_search_boundary": False,
        "detail": (
            f"dense scan at step {COARSE_EXPONENT_STEP} then golden section"
        ),
    }


def fit_fixed_order(
    h: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
    *,
    order: float,
) -> dict[str, Any]:
    r"""Richardson extrapolation with the exponent *assumed* equal to ``order``.

    This is the classical extrapolation used when the order of the scheme is
    known a priori: ``A(h) = A_inf + C h^order`` is linear in its two remaining
    parameters, so ``(A_inf, C)`` follow from a single linear least squares.

    ``order`` is an assumption, not a measurement.  The returned record carries
    ``"p_is_free": False`` for that reason, and :func:`subset_sensitivity`
    excludes these exponents from its ``p_min``/``p_max`` range while keeping
    their ``A_inf`` in the spread -- see the module docstring.
    """

    spacings, magnitudes = _validated_series(h, values)
    if isinstance(order, bool) or not isinstance(order, (int, float)):
        raise ValueError("order must be a real number")
    exponent = float(order)
    if not math.isfinite(exponent):
        raise ValueError("order must be finite")
    if exponent <= 0.0:
        raise ValueError("order must be strictly positive")

    scale = float(spacings[0])
    ratios = spacings / scale
    fit = _linear_fit_at_exponent(ratios, magnitudes, exponent)
    if fit is None:
        return _failed_fit(
            "fixed_order",
            spacings.size,
            f"design matrix is collinear at the assumed order {exponent}",
            order=exponent,
        )

    limit, slope, residual = fit
    scale_power = scale**exponent
    coefficient = slope / scale_power if scale_power > 0.0 else math.inf
    if not all(math.isfinite(value) for value in (limit, coefficient)):
        return _failed_fit(
            "fixed_order",
            spacings.size,
            "fit produced a non-finite parameter",
            order=exponent,
        )

    return {
        "model": "fixed_order",
        "n_points": int(spacings.size),
        "p_is_free": False,
        "A_inf": float(limit),
        "C": float(coefficient),
        "p": exponent,
        "residual_norm": float(residual),
        "converged": True,
        "order": exponent,
        "detail": f"linear least squares at assumed order {exponent}",
    }


def _summarize(fits: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Reduce a fit table to the preregistered sensitivity summary."""

    converged = [fit for fit in fits if fit["converged"]]
    free_exponent = [fit for fit in converged if fit["p_is_free"]]
    limits = [float(fit["A_inf"]) for fit in converged]
    exponents = [float(fit["p"]) for fit in free_exponent]

    if limits:
        limit_min = min(limits)
        limit_max = max(limits)
        median_absolute = float(np.median(np.abs(np.asarray(limits))))
        if median_absolute > 0.0:
            spread: float | None = (limit_max - limit_min) / median_absolute
        elif limit_max == limit_min:
            spread = 0.0
        else:
            # Extrapolants straddle zero with a zero median: no meaningful
            # relative scale exists, so the spread is not assessable.
            spread = None
    else:
        limit_min = limit_max = median_absolute = None
        spread = None

    if exponents:
        exponent_min = min(exponents)
        exponent_max = max(exponents)
        exponent_spread: float | None = exponent_max - exponent_min
    else:
        exponent_min = exponent_max = exponent_spread = None

    enough_fits = (
        len(converged) >= MINIMUM_CONVERGED_FITS
        and len(free_exponent) >= MINIMUM_FREE_EXPONENT_FITS
    )
    in_asymptotic_range = bool(
        enough_fits
        and spread is not None
        and exponent_spread is not None
        and spread <= A_INF_SPREAD_TOLERANCE
        and exponent_spread <= P_SPREAD_TOLERANCE
    )

    return {
        "A_inf_min": limit_min,
        "A_inf_max": limit_max,
        "A_inf_median_absolute": median_absolute,
        "A_inf_spread_relative": spread,
        "p_min": exponent_min,
        "p_max": exponent_max,
        "p_spread": exponent_spread,
        "in_asymptotic_range": in_asymptotic_range,
        "converged_fit_count": len(converged),
        "free_exponent_fit_count": len(free_exponent),
        "total_fit_count": len(fits),
        "sufficient_fits_for_assessment": bool(enough_fits),
    }


def subset_sensitivity(
    h: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
) -> dict[str, Any]:
    r"""Refit the series on every subset and report how far the answer moves.

    The table contains, with the series sorted coarsest-first:

    * :func:`fit_power_law_three_point` on every contiguous triple;
    * :func:`fit_power_law_least_squares` on every contiguous window of three
      or more points.  The full set is the longest such window and therefore
      appears exactly once, flagged by ``"is_full_set"``;
    * :func:`fit_fixed_order` on the full set at assumed orders ``1`` and
      ``2``.

    Each record carries ``"indices"`` and ``"h_window"`` so that a disagreement
    can be traced back to the subset that produced it.

    The summary applies the thresholds preregistered in this module -- see
    :data:`A_INF_SPREAD_TOLERANCE` and :data:`P_SPREAD_TOLERANCE` -- and reports
    ``in_asymptotic_range``.  **A ``False`` verdict means no limit estimate may
    be quoted from this series**, including averages of the table and the value
    from the finest window; what is reportable is the spread itself.  Fits that
    did not converge are excluded from the summary but retained in the table.
    """

    spacings, magnitudes = _validated_series(h, values)
    count = int(spacings.size)
    fits: list[dict[str, Any]] = []

    def record(fit: dict[str, Any], start: int, stop: int, **extra: Any) -> None:
        fit = dict(fit)
        fit["indices"] = tuple(range(start, stop))
        fit["h_window"] = tuple(float(value) for value in spacings[start:stop])
        fit["is_full_set"] = (start == 0 and stop == count)
        fit.update(extra)
        fits.append(fit)

    for start in range(count - 2):
        stop = start + 3
        record(
            fit_power_law_three_point(
                spacings[start:stop], magnitudes[start:stop]
            ),
            start,
            stop,
        )

    for width in range(3, count + 1):
        for start in range(count - width + 1):
            stop = start + width
            record(
                fit_power_law_least_squares(
                    spacings[start:stop], magnitudes[start:stop]
                ),
                start,
                stop,
            )

    for order in (1, 2):
        record(fit_fixed_order(spacings, magnitudes, order=order), 0, count)

    return {
        "h": tuple(float(value) for value in spacings),
        "values": tuple(float(value) for value in magnitudes),
        "n_points": count,
        "fits": tuple(fits),
        "summary": _summarize(fits),
        "thresholds": {
            "A_inf_spread_relative_max": A_INF_SPREAD_TOLERANCE,
            "p_spread_max": P_SPREAD_TOLERANCE,
            "minimum_converged_fits": MINIMUM_CONVERGED_FITS,
            "minimum_free_exponent_fits": MINIMUM_FREE_EXPONENT_FITS,
            "exponent_search_range": (MINIMUM_EXPONENT, MAXIMUM_EXPONENT),
            "preregistered_in": "ns_certificate_lab.extrapolation",
        },
    }


def blind_extrapolation_report(
    h: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
    **forbidden_keywords: Any,
) -> dict[str, Any]:
    r"""Blind extrapolation entry point: subset sensitivity plus a verdict.

    This is the function experiments are expected to call.  It takes the
    spacings and the values and **nothing else**.  There is no way to supply a
    published limit, a target, an anchor or an expected value: every extra
    keyword argument raises :class:`ValueError`, so an anchored fit cannot be
    expressed through this API and would have to be written somewhere a
    reviewer can see it.  See the module docstring for why anchored fits are
    forbidden.

    The ``"verdict"`` field is one of

    * :data:`NOT_ASYMPTOTIC_RANGE_VERDICT` -- the extrapolants disagree by more
      than the preregistered spread, so the series has not entered the
      asymptotic range of the model and **no limit estimate is quotable from
      it**; or
    * :data:`ASYMPTOTIC_RANGE_VERDICT` -- every converged extrapolant agrees
      within the preregistered spread.

    Neither verdict is ever a comparison against an external number.  The
    second one asserts that the series is internally consistent with the fitted
    model, which is a necessary condition for quoting a limit and not a
    sufficient one.
    """

    if forbidden_keywords:
        offered = sorted(forbidden_keywords)
        anchoring = sorted(
            name
            for name in offered
            if any(stem in name.lower() for stem in ANCHORING_KEYWORD_STEMS)
        )
        message = (
            "blind_extrapolation_report takes only (h, values); refusing the "
            f"extra keyword argument(s) {offered}"
        )
        if anchoring:
            message += (
                f".  {anchoring} name(s) an external value: fitting this "
                "series with a published or expected limit held fixed is "
                "forbidden by the project audit, because two free parameters "
                "will absorb any short series and the fit then measures "
                "nothing"
            )
        raise ValueError(message)

    report = subset_sensitivity(h, values)
    report["verdict"] = (
        ASYMPTOTIC_RANGE_VERDICT
        if report["summary"]["in_asymptotic_range"]
        else NOT_ASYMPTOTIC_RANGE_VERDICT
    )
    return report
