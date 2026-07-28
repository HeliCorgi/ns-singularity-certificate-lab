"""Verification of the blind extrapolation and subset-sensitivity module.

Four things are checked here, in order of increasing importance.

1.  *Correctness on a case with a known answer.*  An exact power law
    ``A_inf + C h^p`` must be recovered by every model in the module, to
    roundoff, and must be declared to be in the asymptotic range.

2.  *Reproduction of the project audit's own arithmetic.*  The audit
    (``FABLE5_NEXT_TASK_AUDIT.md``, section P1-A) quotes three blind fits of the
    real vorticity-amplification series.  This module must reproduce them.  If
    it did not, the audit's diagnosis would rest on numbers nobody had checked.
    An independent brute-force grid search, written in this file from the raw
    normal equations rather than from the module's scaled and centered algebra,
    is used as the arbiter; the module is asserted against *both* the audit and
    the brute force, so a coincidental agreement with one of them cannot pass.

3.  *The audit's conclusion, encoded as a test.*  The real series must come out
    ``in_asymptotic_range = False``, and the verdict string must say that no
    limit estimate is quotable.

4.  *That the API cannot express an anchored fit.*  No function may accept a
    reference, target or anchor value, the published number ``20.5235`` must
    not appear anywhere in the module, and the blind entry point must reject
    every extra keyword argument.

No expected value in this file is produced by the module under test: the
synthetic series are written out in closed form, the audit's numbers are
transcribed from the audit, and the brute-force arbiter is implemented here.

Measured constants quoted in the assertions were obtained on this machine with
the pinned environment; the assertions keep visible margin so that they pin
behaviour rather than arithmetic noise.
"""

from __future__ import annotations

import inspect
import math
from pathlib import Path

import numpy as np
import pytest

import ns_certificate_lab.extrapolation as extrapolation_module
from ns_certificate_lab.extrapolation import (
    ANCHORING_KEYWORD_STEMS,
    A_INF_SPREAD_TOLERANCE,
    ASYMPTOTIC_RANGE_VERDICT,
    MAXIMUM_EXPONENT,
    MINIMUM_EXPONENT,
    NOT_ASYMPTOTIC_RANGE_VERDICT,
    P_SPREAD_TOLERANCE,
    blind_extrapolation_report,
    fit_fixed_order,
    fit_power_law_least_squares,
    fit_power_law_three_point,
    subset_sensitivity,
)

# The uniform-grid spacings of the Hou early-time resolution ladder.
SPACINGS = (1.0 / 64.0, 1.0 / 128.0, 1.0 / 192.0, 1.0 / 256.0)

# The measured grid-normalized vorticity amplifications at those spacings.
AMPLIFICATIONS = (6.11, 12.70, 15.63, 17.26)

# Transcribed from FABLE5_NEXT_TASK_AUDIT.md, section P1-A ("概算"), which
# quotes blind fits of A(h) = A_inf + C h^p to the series above:
#     all four points : A_inf ~ 27.38, p ~ 0.54
#     first three     : A_inf ~ 28.85, p ~ 0.49
#     last three      : A_inf ~ 24.60, p ~ 0.70
AUDIT_FULL_SET = (27.38, 0.54)
AUDIT_FIRST_TRIPLE = (28.85, 0.49)
AUDIT_LAST_TRIPLE = (24.60, 0.70)


def _brute_force_power_law(
    h: tuple[float, ...], values: tuple[float, ...], *, step: float = 1.0e-4
) -> tuple[float, float, float, float]:
    """Return ``(A_inf, C, p, residual_norm)`` from an exhaustive scan over ``p``.

    Deliberately independent of the module under test.  The exponent is swept
    on a uniform grid of resolution ``step`` over the same range the module
    searches, and for each candidate the linear pair ``(A_inf, C)`` is obtained
    from the raw ``2 x 2`` normal equations of the design ``[1, h^p]``:

        A_inf = (S_ss T_y  - S_s T_sy) / D,
        C     = (n   T_sy - S_s T_y ) / D,   D = n S_ss - S_s^2,

    with ``S_s = sum h^p``, ``S_ss = sum h^{2p}``, ``T_y = sum A`` and
    ``T_sy = sum h^p A``.  The module instead centers the data, rescales the
    basis by ``h_max^p`` and refines with golden section, so an error in either
    derivation cannot cancel against the other.  Candidates whose normal matrix
    is singular are skipped; residuals are always formed explicitly, so an
    ill-conditioned candidate cannot win by cancellation.
    """

    spacings = np.asarray(h, dtype=np.float64)
    magnitudes = np.asarray(values, dtype=np.float64)
    count = spacings.size

    steps = int(round((MAXIMUM_EXPONENT - MINIMUM_EXPONENT) / step))
    exponents = MINIMUM_EXPONENT + step * np.arange(steps + 1)

    basis = spacings[None, :] ** exponents[:, None]
    basis_sum = basis.sum(axis=1)
    square_sum = (basis * basis).sum(axis=1)
    value_sum = float(magnitudes.sum())
    cross_sum = basis @ magnitudes

    determinant = count * square_sum - basis_sum * basis_sum
    usable = determinant > 0.0
    safe = np.where(usable, determinant, 1.0)
    limit = (square_sum * value_sum - basis_sum * cross_sum) / safe
    coefficient = (count * cross_sum - basis_sum * value_sum) / safe
    residual = magnitudes[None, :] - (
        limit[:, None] + coefficient[:, None] * basis
    )
    norms = np.where(usable, np.linalg.norm(residual, axis=1), np.inf)

    index = int(np.argmin(norms))
    return (
        float(limit[index]),
        float(coefficient[index]),
        float(exponents[index]),
        float(norms[index]),
    )


def _fit_for(
    report: dict[str, object], model: str, indices: tuple[int, ...]
) -> dict[str, object]:
    matches = [
        fit
        for fit in report["fits"]
        if fit["model"] == model and fit["indices"] == indices
    ]
    assert len(matches) == 1, (model, indices)
    return matches[0]


def test_exact_power_law_is_recovered_by_every_model() -> None:
    """``A(h) = 20 + 100 h^2`` must be recovered to roundoff by all three models.

    The series is written out in closed form here, so the expected limit and
    exponent owe nothing to the module.  Measured on this machine: the
    three-point solves return ``A_inf`` exactly ``20.0`` with ``p`` in error by
    ``1.2e-12`` and ``1.2e-12``; the least-squares fits return ``A_inf`` in
    error by ``4e-15`` with ``p`` in error by ``8e-12``.
    """

    limit, coefficient, exponent = 20.0, 100.0, 2.0
    values = tuple(limit + coefficient * h**exponent for h in SPACINGS)

    for start in range(2):
        window = slice(start, start + 3)
        exact = fit_power_law_three_point(SPACINGS[window], values[window])
        assert exact["converged"]
        assert abs(exact["A_inf"] - limit) <= 1.0e-8
        assert abs(exact["p"] - exponent) <= 1.0e-6
        assert exact["C"] == pytest.approx(coefficient, rel=1.0e-6)
        # The third equation was eliminated by the solve, so its defect is an
        # independent check.  Measured: 3.6e-15.
        assert exact["residual_norm"] <= 1.0e-10

    for window in (slice(0, 3), slice(1, 4), slice(0, 4)):
        fitted = fit_power_law_least_squares(SPACINGS[window], values[window])
        assert fitted["converged"]
        assert abs(fitted["A_inf"] - limit) <= 1.0e-8
        assert abs(fitted["p"] - exponent) <= 1.0e-6
        assert fitted["C"] == pytest.approx(coefficient, rel=1.0e-6)
        assert fitted["residual_norm"] <= 1.0e-10

    # The assumed-order model is exact when the assumption is right, and wrong
    # when it is not.  Measured for order 1: A_inf = 19.9923.
    correct_order = fit_fixed_order(SPACINGS, values, order=2)
    assert correct_order["converged"]
    assert abs(correct_order["A_inf"] - limit) <= 1.0e-8
    assert correct_order["C"] == pytest.approx(coefficient, rel=1.0e-8)
    assert correct_order["p"] == 2.0
    assert correct_order["p_is_free"] is False

    wrong_order = fit_fixed_order(SPACINGS, values, order=1)
    assert wrong_order["converged"]
    assert wrong_order["A_inf"] < limit - 1.0e-4

    report = blind_extrapolation_report(SPACINGS, values)
    summary = report["summary"]
    assert summary["in_asymptotic_range"] is True
    assert report["verdict"] == ASYMPTOTIC_RANGE_VERDICT
    # Measured: spread 3.8e-4 (set by the deliberately wrong order-1 fit) and
    # exponent spread 9.2e-12.
    assert summary["A_inf_spread_relative"] <= A_INF_SPREAD_TOLERANCE
    assert summary["p_spread"] <= P_SPREAD_TOLERANCE
    assert summary["converged_fit_count"] == summary["total_fit_count"] == 7
    assert summary["free_exponent_fit_count"] == 5


def test_real_series_reproduces_the_audit_and_an_independent_brute_force() -> None:
    """The module must reproduce the audit's P1-A fits and a brute-force arbiter.

    Measured on this machine:

    ==================  =========  ========  =========  ========
    subset              A_inf      p         audit      audit p
    ==================  =========  ========  =========  ========
    all four points     27.37705   0.53743   27.38      0.54
    first three         28.85368   0.49360   28.85      0.49
    last three          24.60029   0.69709   24.60      0.70
    ==================  =========  ========  =========  ========

    Relative gaps against the audit: 1.0e-4, 1.3e-5 and 1.2e-5 on ``A_inf``;
    absolute gaps on ``p``: 0.0026, 0.0036 and 0.0029.  Both are far inside the
    2 percent / 0.05 tolerances the audit's two-significant-figure numbers
    deserve.

    The brute-force arbiter agrees with the module to 8.5e-5 relative on
    ``A_inf`` and 3.0e-5 absolute on ``p``, the latter being half its own grid
    resolution, as it must be.
    """

    windows = {
        (0, 1, 2, 3): AUDIT_FULL_SET,
        (0, 1, 2): AUDIT_FIRST_TRIPLE,
        (1, 2, 3): AUDIT_LAST_TRIPLE,
    }

    for indices, (audit_limit, audit_exponent) in windows.items():
        spacings = tuple(SPACINGS[index] for index in indices)
        values = tuple(AMPLIFICATIONS[index] for index in indices)

        fitted = fit_power_law_least_squares(spacings, values)
        assert fitted["converged"]

        brute_limit, brute_coefficient, brute_exponent, _ = _brute_force_power_law(
            spacings, values
        )

        # (a) The module agrees with the independent brute-force scan.  The
        # tolerance on the exponent is set by the arbiter's own 1e-4 grid.
        assert fitted["A_inf"] == pytest.approx(brute_limit, rel=5.0e-3)
        assert abs(fitted["p"] - brute_exponent) <= 5.0e-4
        assert fitted["C"] == pytest.approx(brute_coefficient, rel=5.0e-2)

        # (b) Both agree with the number the audit quotes.
        assert fitted["A_inf"] == pytest.approx(audit_limit, rel=0.02)
        assert abs(fitted["p"] - audit_exponent) <= 0.05
        assert brute_limit == pytest.approx(audit_limit, rel=0.02)
        assert abs(brute_exponent - audit_exponent) <= 0.05

        # (c) On three points the exact solve and the least-squares fit are the
        # same problem, so the two independent implementations must coincide
        # and the least-squares residual must vanish.  Measured: 8.2e-14.
        if len(indices) == 3:
            exact = fit_power_law_three_point(spacings, values)
            assert exact["converged"]
            assert exact["A_inf"] == pytest.approx(fitted["A_inf"], rel=1.0e-9)
            assert abs(exact["p"] - fitted["p"]) <= 1.0e-9
            assert fitted["residual_norm"] <= 1.0e-10
        else:
            # The four-point fit cannot be exact; if it were, the three subsets
            # would not disagree and there would be nothing to diagnose.
            # Measured: 0.0487.
            assert fitted["residual_norm"] > 1.0e-3

    # The same three fits, reached through the reporting entry point rather
    # than by calling the fitters directly.
    report = blind_extrapolation_report(SPACINGS, AMPLIFICATIONS)
    assert _fit_for(report, "least_squares", (0, 1, 2, 3))["A_inf"] == pytest.approx(
        AUDIT_FULL_SET[0], rel=0.02
    )
    assert _fit_for(report, "three_point", (0, 1, 2))["A_inf"] == pytest.approx(
        AUDIT_FIRST_TRIPLE[0], rel=0.02
    )
    assert _fit_for(report, "three_point", (1, 2, 3))["A_inf"] == pytest.approx(
        AUDIT_LAST_TRIPLE[0], rel=0.02
    )


def test_real_series_is_not_in_the_asymptotic_range() -> None:
    """The audit's conclusion, encoded: this series supports no limit estimate.

    Measured on this machine, over the seven fits of the full table:

    ===================  =========  ========
    model                A_inf      p
    ===================  =========  ========
    three_point (0,1,2)  28.85368   0.49360
    three_point (1,2,3)  24.60029   0.69709
    lsq (0,1,2)          28.85368   0.49360
    lsq (1,2,3)          24.60029   0.69709
    lsq (0,1,2,3)        27.37705   0.53743
    fixed order 1        20.52821   1 (assumed)
    fixed order 2        16.84418   2 (assumed)
    ===================  =========  ========

    ``A_inf`` spans ``[16.84, 28.85]`` about a median of ``24.60``, a relative
    spread of ``0.488`` -- nearly ten times the preregistered ``0.05``.  The
    free exponents span ``0.2035``, which also exceeds the preregistered
    ``0.2``, but only just; the verdict here is driven by the ``A_inf`` spread,
    and the assertions below say so separately rather than relying on the
    marginal criterion.

    The order-1 Richardson extrapolant, ``20.528``, happens to land within
    ``0.02`` percent of the externally published ``20.5235``.  That coincidence
    is recorded here precisely so that nobody mistakes it for evidence: it is
    one of seven mutually inconsistent extrapolants of a series that has not
    entered the asymptotic range, and it is not more privileged than the
    ``16.84`` produced by assuming second order instead of first.
    """

    report = blind_extrapolation_report(SPACINGS, AMPLIFICATIONS)
    summary = report["summary"]

    assert summary["in_asymptotic_range"] is False
    assert report["verdict"] == NOT_ASYMPTOTIC_RANGE_VERDICT
    assert report["verdict"].startswith("not_in_asymptotic_range")
    assert "no limit estimate is quotable" in report["verdict"]

    # The A_inf criterion fails on its own, by a wide margin.
    assert summary["A_inf_spread_relative"] > 5.0 * A_INF_SPREAD_TOLERANCE
    assert summary["A_inf_min"] < 20.0 < summary["A_inf_max"]
    assert summary["A_inf_max"] - summary["A_inf_min"] > 10.0

    # The exponent criterion also fails, but marginally; recorded as its own
    # assertion so that a future change to either criterion is visible.
    assert summary["p_spread"] > P_SPREAD_TOLERANCE
    assert summary["p_spread"] < 0.25

    assert summary["converged_fit_count"] == 7
    assert summary["free_exponent_fit_count"] == 5
    assert summary["sufficient_fits_for_assessment"] is True

    # Every fit in the table is monotone-consistent with the series but they do
    # not agree with each other; that disagreement is the reportable content.
    limits = sorted(fit["A_inf"] for fit in report["fits"] if fit["converged"])
    assert limits[0] < 17.0
    assert limits[-1] > 28.0

    # The thresholds travel with the report, so saved evidence carries the rule
    # that was applied to it.
    assert report["thresholds"]["A_inf_spread_relative_max"] == A_INF_SPREAD_TOLERANCE
    assert report["thresholds"]["p_spread_max"] == P_SPREAD_TOLERANCE
    assert report["thresholds"]["preregistered_in"] == "ns_certificate_lab.extrapolation"


def test_non_monotone_series_reports_non_convergence_without_crashing() -> None:
    """A triple with no admissible exponent must report it, not invent one.

    ``g(p) = (h2^p - h3^p)/(h1^p - h2^p)`` is strictly positive for every
    ``p > 0``, so a triple whose two increments have opposite signs has no root
    anywhere in the searched range.  The series below is the real one with its
    finest value perturbed downward, which makes the last triple non-monotone
    while leaving the first one untouched.

    Measured on this machine: the last triple's increment ratio is ``-0.1809``,
    the three-point solve reports no root, the least-squares fit on the same
    window drives its exponent to the edge of the search range and is also
    reported as non-converged, and the summary is computed from the remaining
    five fits.
    """

    values = (6.11, 12.70, 15.63, 15.10)
    report = blind_extrapolation_report(SPACINGS, values)

    failed_triple = _fit_for(report, "three_point", (1, 2, 3))
    assert failed_triple["converged"] is False
    assert failed_triple["A_inf"] is None
    assert failed_triple["C"] is None
    assert failed_triple["p"] is None
    assert failed_triple["residual_norm"] is None
    assert "no exponent" in failed_triple["detail"]

    failed_window = _fit_for(report, "least_squares", (1, 2, 3))
    assert failed_window["converged"] is False
    assert failed_window["A_inf"] is None
    assert failed_window["exponent_at_search_boundary"] is True

    # The untouched first triple is unaffected: non-convergence does not
    # contaminate its neighbours.
    intact = _fit_for(report, "three_point", (0, 1, 2))
    assert intact["converged"] is True
    assert intact["A_inf"] == pytest.approx(AUDIT_FIRST_TRIPLE[0], rel=0.02)

    summary = report["summary"]
    assert summary["total_fit_count"] == 7
    assert summary["converged_fit_count"] == 5
    assert summary["free_exponent_fit_count"] == 3
    # The summary is built only from converged fits: every extremum it reports
    # is attained by one of them.
    converged_limits = [
        fit["A_inf"] for fit in report["fits"] if fit["converged"]
    ]
    assert summary["A_inf_min"] == min(converged_limits)
    assert summary["A_inf_max"] == max(converged_limits)
    assert summary["in_asymptotic_range"] is False
    assert report["verdict"] == NOT_ASYMPTOTIC_RANGE_VERDICT

    # A directly called fitter reports the same failure, so the behaviour is
    # the fitter's and not the table builder's.
    direct = fit_power_law_three_point(SPACINGS[1:], values[1:])
    assert direct["converged"] is False
    assert direct["A_inf"] is None


def test_non_power_law_series_is_not_declared_asymptotic() -> None:
    """Monotone data from outside the model family must not pass the criteria.

    ``A(h) = 20 - 10 exp(-h^{-1/2})`` is monotone, bounded and converges to
    ``20``, but its approach is exponential rather than algebraic, so no single
    ``h^p`` describes it.  Every subset therefore reports a different exponent
    even though every subset happens to land near the right limit.

    This is the case the ``A_inf`` spread alone would wave through: measured on
    this machine the spread is ``8.0e-5``, comfortably inside the preregistered
    ``0.05``.  What catches it is the exponent range, measured at ``1.459``
    across free-exponent fits (``4.668`` on the first triple, ``6.127`` on the
    last, ``4.717`` on the full set).  Requiring *both* criteria is what makes
    the verdict honest here.
    """

    spacings = np.asarray(SPACINGS, dtype=np.float64)
    values = 20.0 - 10.0 * np.exp(-(spacings**-0.5))
    # The series really is monotone increasing towards 20 and never reaches it.
    assert np.all(np.diff(values) > 0.0)
    assert float(values[-1]) < 20.0

    report = blind_extrapolation_report(spacings, values)
    summary = report["summary"]

    assert summary["in_asymptotic_range"] is False
    assert report["verdict"] == NOT_ASYMPTOTIC_RANGE_VERDICT

    # The limit spread alone would have passed; the exponent spread is what
    # rejects the series.
    assert summary["A_inf_spread_relative"] <= A_INF_SPREAD_TOLERANCE
    assert summary["p_spread"] > 5.0 * P_SPREAD_TOLERANCE
    assert summary["p_min"] > 4.0
    assert summary["p_max"] < MAXIMUM_EXPONENT


def test_fitters_reject_malformed_series() -> None:
    """Every fitter validates its inputs identically and refuses bad ones.

    All five entry points share one validator, and the shape/finiteness checks
    run before any model-specific check, so the same malformed series produces
    the same message everywhere.  That is asserted rather than assumed: each
    case below is run through every entry point.
    """

    good_h = SPACINGS
    good_values = AMPLIFICATIONS

    def fitters(h, values):
        yield lambda: fit_power_law_three_point(h, values)
        yield lambda: fit_power_law_least_squares(h, values)
        yield lambda: fit_fixed_order(h, values, order=2)
        yield lambda: subset_sensitivity(h, values)
        yield lambda: blind_extrapolation_report(h, values)

    # Mismatched lengths.
    for call in fitters(good_h, good_values[:-1]):
        with pytest.raises(ValueError, match="same length"):
            call()

    # Fewer than three points.
    for call in fitters(good_h[:2], good_values[:2]):
        with pytest.raises(ValueError, match="at least three"):
            call()

    # Non-positive spacings.
    for bad in ((0.0, 0.5, 0.25, 0.125), (-0.25, 0.5, 0.3, 0.125)):
        for call in fitters(bad, good_values):
            with pytest.raises(ValueError, match="strictly positive"):
                call()

    # Non-finite spacings and values.
    for call in fitters((math.inf, 0.5, 0.25, 0.125), good_values):
        with pytest.raises(ValueError, match="h must be finite"):
            call()
    for bad_values in (
        (6.11, math.nan, 15.63, 17.26),
        (6.11, 12.70, math.inf, 17.26),
    ):
        for call in fitters(good_h, bad_values):
            with pytest.raises(ValueError, match="values must be finite"):
                call()

    # Duplicated spacings.
    for call in fitters((0.5, 0.25, 0.25, 0.125), good_values):
        with pytest.raises(ValueError, match="duplicated"):
            call()

    # Two-dimensional input.
    for call in fitters(np.reshape(good_h, (2, 2)), np.reshape(good_values, (2, 2))):
        with pytest.raises(ValueError, match="one-dimensional"):
            call()

    # Non-numeric input.
    for call in fitters(("a", "b", "c", "d"), good_values):
        with pytest.raises(ValueError, match="real numeric"):
            call()

    # The exact solve takes exactly three points, never more.
    with pytest.raises(ValueError, match="exactly three"):
        fit_power_law_three_point(good_h, good_values)

    # The assumed order must be a positive, finite, real number.
    for bad_order in (0.0, -1.0, math.nan, math.inf):
        with pytest.raises(ValueError):
            fit_fixed_order(good_h, good_values, order=bad_order)
    for bad_order in (True, "2", None):
        with pytest.raises(ValueError, match="real number"):
            fit_fixed_order(good_h, good_values, order=bad_order)


def test_blind_report_rejects_every_extra_keyword_argument() -> None:
    """The entry point takes ``(h, values)`` and refuses anything else.

    The point is not that the extra argument would be ignored -- it is that
    accepting one at all would make an anchored fit expressible through this
    API.  Names that suggest an external value get a sharper message, but the
    refusal does not depend on guessing the name: any keyword is rejected.
    """

    for name in ("reference", "target", "anchor", "published_value", "true_limit"):
        with pytest.raises(ValueError, match="forbidden by the project audit"):
            blind_extrapolation_report(SPACINGS, AMPLIFICATIONS, **{name: 20.5235})

    # A neutral name is refused too, just without the anchoring sentence.
    with pytest.raises(ValueError, match="takes only"):
        blind_extrapolation_report(SPACINGS, AMPLIFICATIONS, verbose=True)
    with pytest.raises(ValueError, match="takes only"):
        blind_extrapolation_report(SPACINGS, AMPLIFICATIONS, x=1, y=2)

    # Positionally the third argument does not exist either.
    with pytest.raises(TypeError):
        blind_extrapolation_report(SPACINGS, AMPLIFICATIONS, 20.5235)


def test_module_cannot_express_an_anchored_fit() -> None:
    """Structural guard: no signature and no literal admits an external limit.

    The prohibition on anchored fits is only real if it is enforced.  Two
    checks: every public entry point's signature is inspected for a parameter
    whose name suggests a reference value, and the module source is scanned for
    the published constant itself.  Neither check can be satisfied by a comment
    saying the right thing.
    """

    public = (
        fit_power_law_three_point,
        fit_power_law_least_squares,
        fit_fixed_order,
        subset_sensitivity,
        blind_extrapolation_report,
    )
    for function in public:
        parameters = inspect.signature(function).parameters
        for name, parameter in parameters.items():
            lowered = name.lower()
            if parameter.kind is inspect.Parameter.VAR_KEYWORD:
                # The only **kwargs in the module exists to reject arguments.
                assert function is blind_extrapolation_report
                continue
            for stem in ANCHORING_KEYWORD_STEMS:
                assert stem not in lowered, (function.__name__, name)
        # Nothing has a default that could smuggle a value in.
        for name, parameter in parameters.items():
            if parameter.default is not inspect.Parameter.empty:
                assert not isinstance(parameter.default, (int, float)), name

    assert list(inspect.signature(blind_extrapolation_report).parameters) == [
        "h",
        "values",
        "forbidden_keywords",
    ]

    source = Path(extrapolation_module.__file__).resolve().read_text(encoding="utf-8")
    # The externally published anchor must not appear, in any form, anywhere in
    # the module -- not as a constant, not as a default, not in a docstring.
    assert "20.5235" not in source
    assert "20.52" not in source

    # The preregistered thresholds are literals in the module, not values read
    # from configuration at run time.
    assert "A_INF_SPREAD_TOLERANCE = 0.05" in source
    assert "P_SPREAD_TOLERANCE = 0.2" in source
    assert A_INF_SPREAD_TOLERANCE == 0.05
    assert P_SPREAD_TOLERANCE == 0.2
    assert (MINIMUM_EXPONENT, MAXIMUM_EXPONENT) == (0.05, 8.0)


def test_report_is_serializable_and_windows_are_ordered() -> None:
    """Records carry their provenance and contain no non-finite numbers.

    This project's evidence writers reject NaN and infinity outright, so a fit
    table that used ``NaN`` for a failed fit could not be saved.  Every optional
    number here is ``None`` instead, which survives JSON.
    """

    report = blind_extrapolation_report(
        tuple(reversed(SPACINGS)), tuple(reversed(AMPLIFICATIONS))
    )

    # Input order does not matter: the series is sorted coarsest-first, so the
    # windows are the same ones the forward-ordered call produces.
    assert report["h"] == SPACINGS
    assert report["values"] == AMPLIFICATIONS

    windows = [(fit["model"], fit["indices"]) for fit in report["fits"]]
    assert windows == [
        ("three_point", (0, 1, 2)),
        ("three_point", (1, 2, 3)),
        ("least_squares", (0, 1, 2)),
        ("least_squares", (1, 2, 3)),
        ("least_squares", (0, 1, 2, 3)),
        ("fixed_order", (0, 1, 2, 3)),
        ("fixed_order", (0, 1, 2, 3)),
    ]
    assert sum(1 for fit in report["fits"] if fit["is_full_set"]) == 3

    for fit in report["fits"]:
        assert fit["h_window"] == tuple(SPACINGS[index] for index in fit["indices"])
        for key in ("A_inf", "C", "p", "residual_norm"):
            value = fit[key]
            assert value is None or math.isfinite(value)
        # A converged fit never carries a ``None`` parameter, and a failed one
        # never carries a number; the two states do not mix.
        carries_numbers = fit["A_inf"] is not None
        assert carries_numbers is bool(fit["converged"])
        assert (fit["p"] is not None) is carries_numbers
        assert (fit["C"] is not None) is carries_numbers

    for key, value in report["summary"].items():
        if isinstance(value, float):
            assert math.isfinite(value), key


def test_three_point_solve_matches_the_closed_form_for_a_geometric_ladder() -> None:
    """On an equally-ratioed ladder the exponent equation has a closed-form root.

    With ``h_k = h_1 r^{k-1}`` the ratio ``(A2-A3)/(A1-A2)`` equals ``r^p``
    exactly, so ``p = log(ratio)/log(r)`` without any root finding.  This pins
    the solver against algebra rather than against a second numerical method.
    """

    ratio = 0.5
    spacings = (0.08, 0.08 * ratio, 0.08 * ratio**2)
    for exponent in (0.3, 1.0, 1.5, 2.0, 3.7):
        limit, coefficient = -4.25, 13.5
        values = tuple(limit + coefficient * h**exponent for h in spacings)

        increment_ratio = (values[1] - values[2]) / (values[0] - values[1])
        closed_form = math.log(increment_ratio) / math.log(ratio)
        assert abs(closed_form - exponent) <= 1.0e-9

        fit = fit_power_law_three_point(spacings, values)
        assert fit["converged"]
        assert abs(fit["p"] - closed_form) <= 1.0e-9
        assert fit["A_inf"] == pytest.approx(limit, rel=1.0e-9, abs=1.0e-9)
        assert fit["C"] == pytest.approx(coefficient, rel=1.0e-8)
