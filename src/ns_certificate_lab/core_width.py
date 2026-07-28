r"""Core-width and points-per-scale diagnostics for snapshot fields (P0-D).

Purpose
-------
Amplification alone is not evidence that a structure is resolved.  A peak that
grows by a factor of a thousand while its core is only four cells wide is a
grid-scale object, and any convergence fit through such data is fitting the
discretization rather than the continuum problem.  The audit item P0-D requires
that this package never again report "no grid-scale saturation" from an
amplitude history alone.

This module answers one question about any snapshot field, in units the grid
itself supplies: *how many grid points resolve the narrowest structure that the
field actually contains?*  It measures

* the full width at half maximum of the peak along both coordinate directions
  (:func:`fwhm_1d`), with explicit truncation flags whenever a half-maximum
  crossing is never reached inside the sampled interval;
* the 10%-to-90% front thickness on a chosen side of the peak
  (:func:`front_thickness_10_90`) -- for a radial profile the interesting side
  is the one facing the axis, because that is where the shipped Hou snapshots
  put their steepest gradient;
* a sub-grid (three-point parabolic) peak location (:func:`subgrid_peak`), so
  that "the peak is six cells from the axis" is not confused with "the peak is
  at the sixth grid point";
* the pointwise gradient length scale ``|f| / |grad f|`` over the region where
  the field is not negligible (:func:`gradient_length_scale`);
* the fraction of discrete spectral energy sitting in the top third of the
  representable wavenumbers (:func:`high_frequency_tail`);
* an assembled report (:func:`peak_resolution_report`) and a hard precondition
  gate (:func:`fit_precondition`).

What this module is not
-----------------------
None of these numbers is a continuum statement.  A points-per-scale count is a
property of the sampled array: it says how the *data* is resolved, not how the
*solution* behaves.  In particular

* the high-frequency tail fractions are resolution indicators computed from a
  single snapshot on a finite grid.  They are not spectra of the continuum
  solution.  Aliasing folds unrepresented continuum content back into the
  represented modes, so a small tail fraction is consistent with a badly
  under-resolved field; only a *large* tail fraction is informative, and then
  only as a warning;
* the gradient length scale is built from three-point differences, which are
  themselves inaccurate exactly where the length scale is smallest.  It is a
  lower bound on trustworthiness, not a measurement of the true scale;
* nothing here is evidence about singularity formation, and nothing here
  licenses a convergence fit.  :func:`fit_precondition` can only *forbid* a
  fit; it never certifies one.

The preregistered threshold
---------------------------
:data:`PREREGISTERED_MIN_POINTS_PER_FRONT` is the minimum number of grid points
per 10%-90% front thickness below which convergence fits are forbidden.  It was
derived from :func:`manufactured_front_study` -- a closed-form ``tanh`` front
whose exact 10-90 thickness and exact derivative full width are known
analytically -- and *not* from any Hou snapshot.  The study is run over a
preregistered ladder of resolutions
(:data:`PREREGISTERED_DELTA_OVER_DR_VALUES`) and over
:data:`PREREGISTERED_FRONT_PHASES` sub-cell placements of the front, and the
reported error is the worst case over those placements: a threshold that only
held for a lucky alignment of the front with the mesh would be worthless.
The constant is the smallest integer number of points per front thickness at
which the worst-case relative front-thickness error stays at or below
:data:`PREREGISTERED_FRONT_ERROR_TOLERANCE`.

The constant is fixed in this file.  ``tests/test_core_width.py`` recomputes the
study and fails if the measured table does not support the stored value, so the
number cannot be quietly retuned to admit a data set that would otherwise be
rejected.  The order matters and is recorded here deliberately: the threshold
was fixed by manufactured fronts *before* being applied to any Hou data.

Conventions
-----------
Fields have shape ``(n_r, n_z)`` and are ``r``-major, following
:class:`~ns_certificate_lab.grid.AxisymmetricGrid`.  The axial direction is
periodic with the right endpoint omitted; axial profiles through a peak are
therefore unwrapped about that peak before their width is measured, so that a
peak sitting near index ``0`` is not reported as truncated when it is in fact
perfectly interior on the circle.  All arithmetic is binary64 and is not
outward rounded.
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid, FloatArray

__all__ = [
    "PREREGISTERED_DELTA_OVER_DR_VALUES",
    "PREREGISTERED_FRONT_ERROR_TOLERANCE",
    "PREREGISTERED_FRONT_PHASES",
    "PREREGISTERED_MIN_POINTS_PER_FRONT",
    "fit_precondition",
    "front_thickness_10_90",
    "fwhm_1d",
    "gradient_length_scale",
    "high_frequency_tail",
    "interpolate_to_common_grid",
    "manufactured_front_study",
    "peak_resolution_report",
    "subgrid_peak",
]

# --------------------------------------------------------------------------
# Preregistered constants.  See the module docstring and
# :func:`manufactured_front_study`.
# --------------------------------------------------------------------------

PREREGISTERED_FRONT_ERROR_TOLERANCE = 0.02
"""Largest relative front-thickness error accepted by the preregistration."""

PREREGISTERED_DELTA_OVER_DR_VALUES: tuple[float, ...] = (
    0.5,
    0.75,
    1.0,
    1.25,
    1.5,
    1.75,
    2.0,
    2.25,
    2.5,
    2.75,
    3.0,
    3.25,
    3.5,
    4.0,
    5.0,
    6.0,
    8.0,
    12.0,
    16.0,
    24.0,
    32.0,
)
"""Preregistered ladder of ``delta / dr`` values for the front study."""

PREREGISTERED_FRONT_PHASES = 8
"""Sub-cell placements of the manufactured front; the error is the worst case."""

PREREGISTERED_MIN_POINTS_PER_FRONT = 7
"""Minimum grid points per 10-90 front thickness required before a fit.

Derivation, reproduced and re-asserted by ``tests/test_core_width.py``.  Run
:func:`manufactured_front_study` over :data:`PREREGISTERED_DELTA_OVER_DR_VALUES`
and split the ladder at :data:`PREREGISTERED_FRONT_ERROR_TOLERANCE`:

* every case with ``delta / dr <= 2.75`` misses the tolerance, and the widest
  front among them spans 6.185 grid points (``delta / dr = 2.75``, worst-case
  relative error 2.362e-02);
* every case with ``delta / dr >= 3.0`` meets it, the narrowest of those
  spanning 6.708 points (worst-case relative error 1.767e-02).

Seven is therefore the smallest integer ``P`` for which *every* case in the
ladder whose measured 10-90 front spans at least ``P`` points meets the
tolerance; ``P = 6`` would admit the failing 6.185-point case.  The count is
expressed in points per measured front thickness -- the same quantity
:func:`peak_resolution_report` reports -- so no conversion is needed at the
gate.

This threshold gates convergence fits and nothing else, and it was fixed from
manufactured fronts before being applied to any Hou snapshot.
"""

_STUDY_DELTA = 0.05
"""Front half-width of the manufactured study, in units of the unit domain."""

_STUDY_CENTER = 0.5
"""Nominal front centre ``r0`` of the manufactured study."""

_SNAP_TOLERANCE = 1.0e-12
"""Cell-fraction below which an interpolation target is treated as a node."""

_LEVEL_HIGH = 0.9
_LEVEL_LOW = 0.1

_SIDES = ("decreasing_coordinate", "increasing_coordinate")


# --------------------------------------------------------------------------
# Validation helpers
# --------------------------------------------------------------------------


def _reject_complex(value: npt.ArrayLike, *, name: str) -> None:
    """Reject complex input instead of silently discarding the imaginary part."""

    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real-valued")
    if raw.dtype.kind == "O" and any(
        isinstance(item, complex) for item in raw.flat
    ):
        raise ValueError(f"{name} must be real-valued")


def _validate_profile(
    coordinates: npt.ArrayLike,
    values: npt.ArrayLike,
) -> tuple[FloatArray, FloatArray]:
    """Return a validated strictly increasing coordinate/value pair."""

    _reject_complex(coordinates, name="coordinates")
    _reject_complex(values, name="values")
    x = np.asarray(coordinates, dtype=np.float64)
    y = np.asarray(values, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("coordinates and values must be one-dimensional")
    if x.size != y.size:
        raise ValueError(
            f"coordinates and values must have equal length, got {x.size} "
            f"and {y.size}"
        )
    if x.size < 3:
        raise ValueError("at least three samples are required")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("coordinates and values must be finite")
    if np.any(np.diff(x) <= 0.0):
        raise ValueError("coordinates must be strictly increasing")
    return x, y


def _resolve_peak_index(
    peak_index: int | None,
    magnitude: FloatArray,
    *,
    required: bool = False,
) -> int:
    """Return a validated peak index, defaulting to ``argmax`` of ``|values|``."""

    if peak_index is None:
        if required:
            raise ValueError("peak_index is required")
        return int(np.argmax(magnitude))
    if isinstance(peak_index, bool) or not isinstance(
        peak_index, (int, np.integer)
    ):
        raise TypeError("peak_index must be an integer")
    index = int(peak_index)
    if not 0 <= index < magnitude.size:
        raise ValueError(
            f"peak_index {index} is outside the profile of length "
            f"{magnitude.size}"
        )
    return index


def _require_grid(grid: object, *, name: str = "grid") -> AxisymmetricGrid:
    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError(f"{name} must be an AxisymmetricGrid")
    return grid


def _validate_positive_fraction(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(
        value, (int, float, np.floating)
    ):
        raise TypeError(f"{name} must be a real number")
    number = float(value)
    if not np.isfinite(number) or not 0.0 < number <= 1.0:
        raise ValueError(f"{name} must lie in (0, 1]")
    return number


# --------------------------------------------------------------------------
# Crossing search
# --------------------------------------------------------------------------


def _descending_crossing(
    x: FloatArray,
    magnitude: FloatArray,
    start_index: int,
    level: float,
    step: int,
) -> tuple[float, bool, int]:
    """Locate the first fall to ``level`` walking away from ``start_index``.

    ``step`` is ``-1`` to walk toward smaller indices and ``+1`` to walk toward
    larger ones.  The returned coordinate is linearly interpolated between the
    two samples that bracket the crossing.  When no sample on that side is at
    or below ``level`` the interval boundary is returned and the truncation
    flag is set: the caller must be able to tell "the width is this" from "the
    width is at least this".

    The third return value is the index of the lower-index sample of the
    bracketing pair (or the boundary index when truncated), which the callers
    use to continue a search for a lower level without walking back past a
    crossing they have already found.
    """

    if step not in (-1, 1):
        raise ValueError("step must be -1 or +1")
    if magnitude[start_index] <= level:
        # The starting sample already sits at or below the level; the crossing
        # is the sample itself.  This happens for the 10% level of a front
        # that falls through both levels inside a single cell.
        return float(x[start_index]), False, int(start_index)

    if step == -1:
        below = np.flatnonzero(magnitude[: start_index + 1] <= level)
        if below.size == 0:
            return float(x[0]), True, 0
        lower = int(below[-1])
    else:
        below = np.flatnonzero(magnitude[start_index:] <= level)
        if below.size == 0:
            return float(x[-1]), True, int(magnitude.size - 1)
        lower = int(start_index + below[0]) - 1

    upper = lower + 1
    span = magnitude[upper] - magnitude[lower]
    if span == 0.0:  # pragma: no cover - excluded by the bracketing above
        raise ArithmeticError("degenerate bracketing interval in crossing search")
    fraction = (level - magnitude[lower]) / span
    fraction = min(1.0, max(0.0, float(fraction)))
    coordinate = float(x[lower] + fraction * (x[upper] - x[lower]))
    return coordinate, False, lower


# --------------------------------------------------------------------------
# One-dimensional width measurements
# --------------------------------------------------------------------------


def fwhm_1d(
    coordinates: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    peak_index: int | None = None,
) -> dict[str, Any]:
    """Return the full width at half maximum of ``|values|`` around its peak.

    Parameters
    ----------
    coordinates:
        Strictly increasing one-dimensional sample coordinates.
    values:
        Real samples of the same length.  The width is measured on ``|values|``.
    peak_index:
        Index of the peak to measure.  Defaults to ``argmax(|values|)``.  An
        explicit index is used when the caller wants the width of a specific
        peak (for example the global peak of a two-dimensional field, whose
        one-dimensional profile may contain a larger local maximum elsewhere).

    Returns
    -------
    dict
        ``fwhm``, ``left`` and ``right`` (the interpolated half-maximum
        crossings), ``peak_value`` (the magnitude at the peak),
        ``peak_coordinate``, and the flags ``truncated_left`` and
        ``truncated_right``.  A truncated side means the profile never fell to
        half maximum before the sampled interval ended; the corresponding
        boundary coordinate is used, so ``fwhm`` is then a *lower bound* and
        must not be reported as a width.
    """

    x, y = _validate_profile(coordinates, values)
    magnitude = np.abs(y)
    index = _resolve_peak_index(peak_index, magnitude)
    peak_value = float(magnitude[index])
    if not peak_value > 0.0:
        raise ValueError(
            "the peak magnitude must be positive to define a half maximum"
        )

    half = 0.5 * peak_value
    left, truncated_left, _ = _descending_crossing(x, magnitude, index, half, -1)
    right, truncated_right, _ = _descending_crossing(x, magnitude, index, half, 1)

    return {
        "fwhm": float(right - left),
        "left": float(left),
        "right": float(right),
        "half_maximum": half,
        "peak_value": peak_value,
        "signed_peak_value": float(y[index]),
        "peak_coordinate": float(x[index]),
        "peak_index": int(index),
        "truncated_left": bool(truncated_left),
        "truncated_right": bool(truncated_right),
    }


def front_thickness_10_90(
    coordinates: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    peak_index: int,
    side: str,
) -> dict[str, Any]:
    """Return the 10%-to-90%-of-peak front thickness on one side of a peak.

    Parameters
    ----------
    coordinates, values:
        As in :func:`fwhm_1d`.
    peak_index:
        Index of the peak whose flank is measured.  Required: a front is a
        property of a side of a specific peak, so there is no useful default.
    side:
        ``"decreasing_coordinate"`` walks toward smaller coordinates (for a
        radial profile this is the flank facing the axis) and
        ``"increasing_coordinate"`` walks toward larger ones.

    Returns
    -------
    dict
        ``thickness`` (the distance between the interpolated 90% and 10%
        crossings), the crossings themselves, the two levels, and the flags
        ``truncated_90``, ``truncated_10`` and ``truncated``.  As in
        :func:`fwhm_1d` a truncated crossing falls back to the interval
        boundary and turns ``thickness`` into a lower bound.
    """

    x, y = _validate_profile(coordinates, values)
    magnitude = np.abs(y)
    index = _resolve_peak_index(peak_index, magnitude, required=True)
    if not isinstance(side, str):
        raise TypeError("side must be a string")
    if side not in _SIDES:
        raise ValueError(f"side must be one of {_SIDES}, got {side!r}")
    peak_value = float(magnitude[index])
    if not peak_value > 0.0:
        raise ValueError(
            "the peak magnitude must be positive to define front levels"
        )

    step = -1 if side == "decreasing_coordinate" else 1
    level_high = _LEVEL_HIGH * peak_value
    level_low = _LEVEL_LOW * peak_value

    # Both levels are located by the same walk away from the peak.  The set of
    # samples at or below the 10% level is contained in the set at or below the
    # 90% level, so the 10% crossing is never found before the 90% one; when a
    # single under-resolved cell spans both levels the two crossings are
    # interpolated inside that one cell and the thickness is a fraction of dr,
    # which is exactly the signal this module exists to expose.
    coordinate_90, truncated_90, _ = _descending_crossing(
        x, magnitude, index, level_high, step
    )
    coordinate_10, truncated_10, _ = _descending_crossing(
        x, magnitude, index, level_low, step
    )

    thickness = abs(coordinate_10 - coordinate_90)
    return {
        "thickness": float(thickness),
        "coordinate_90": float(coordinate_90),
        "coordinate_10": float(coordinate_10),
        "level_90": level_high,
        "level_10": level_low,
        "peak_value": peak_value,
        "peak_coordinate": float(x[index]),
        "peak_index": int(index),
        "side": side,
        "truncated_90": bool(truncated_90),
        "truncated_10": bool(truncated_10),
        "truncated": bool(truncated_90 or truncated_10),
    }


def subgrid_peak(
    coordinates: npt.ArrayLike,
    values: npt.ArrayLike,
    peak_index: int,
) -> dict[str, Any]:
    """Refine a discrete peak with the parabola through its three samples.

    The parabola is built by divided differences through
    ``(x[i-1], y[i-1])``, ``(x[i], y[i])`` and ``(x[i+1], y[i+1])``, so it is
    correct for non-uniform coordinates as well; for uniform spacing it reduces
    to the familiar three-point vertex formula.  ``values`` is used with its
    sign, so the same routine refines a maximum and a minimum.

    Refinement is declined -- ``refined`` is ``False`` and the discrete sample
    is returned unchanged -- when ``peak_index`` is at either end of the
    profile, when the three samples are collinear (no vertex), or when the
    vertex falls outside the bracketing interval (in which case ``peak_index``
    is not the extremum and extrapolating would invent one).
    """

    x, y = _validate_profile(coordinates, values)
    index = _resolve_peak_index(peak_index, np.abs(y), required=True)
    discrete_coordinate = float(x[index])
    discrete_value = float(y[index])
    spacing = float(x[1] - x[0])

    result: dict[str, Any] = {
        "coordinate": discrete_coordinate,
        "value": discrete_value,
        "refined": False,
        "index": int(index),
        "discrete_coordinate": discrete_coordinate,
        "discrete_value": discrete_value,
        "offset_cells": 0.0,
        "curvature": 0.0,
        "reason": "",
    }

    if index == 0 or index == x.size - 1:
        result["reason"] = "peak_index is at a profile boundary"
        return result

    x0, x1, x2 = float(x[index - 1]), float(x[index]), float(x[index + 1])
    y0, y1, y2 = float(y[index - 1]), float(y[index]), float(y[index + 1])
    slope_01 = (y1 - y0) / (x1 - x0)
    slope_12 = (y2 - y1) / (x2 - x1)
    curvature = (slope_12 - slope_01) / (x2 - x0)
    result["curvature"] = float(curvature)
    if curvature == 0.0:
        result["reason"] = "the three samples are collinear"
        return result

    vertex = 0.5 * (x0 + x1) - 0.5 * slope_01 / curvature
    if not np.isfinite(vertex) or not x0 <= vertex <= x2:
        result["reason"] = "the parabolic vertex falls outside the stencil"
        return result

    value = (
        y0
        + slope_01 * (vertex - x0)
        + curvature * (vertex - x0) * (vertex - x1)
    )
    result["coordinate"] = float(vertex)
    result["value"] = float(value)
    result["refined"] = True
    result["offset_cells"] = float((vertex - x1) / spacing) if spacing else 0.0
    return result


# --------------------------------------------------------------------------
# Derivatives
# --------------------------------------------------------------------------


def _radial_derivative(grid: AxisymmetricGrid, field: FloatArray) -> FloatArray:
    """Centered radial difference, second-order one-sided at both ends."""

    dr = grid.dr
    out = np.empty_like(field)
    out[1:-1] = (field[2:] - field[:-2]) / (2.0 * dr)
    out[0] = (-3.0 * field[0] + 4.0 * field[1] - field[2]) / (2.0 * dr)
    out[-1] = (3.0 * field[-1] - 4.0 * field[-2] + field[-3]) / (2.0 * dr)
    return out


def _axial_derivative(grid: AxisymmetricGrid, field: FloatArray) -> FloatArray:
    """Centered axial difference; periodic via :func:`numpy.roll` when asked."""

    dz = grid.dz
    if grid.periodic_z:
        return (
            np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)
        ) / (2.0 * dz)
    out = np.empty_like(field)
    out[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2.0 * dz)
    out[:, 0] = (-3.0 * field[:, 0] + 4.0 * field[:, 1] - field[:, 2]) / (2.0 * dz)
    out[:, -1] = (
        3.0 * field[:, -1] - 4.0 * field[:, -2] + field[:, -3]
    ) / (2.0 * dz)
    return out


def _derivative_1d(x: FloatArray, y: FloatArray) -> FloatArray:
    """Centered difference of a uniformly sampled profile, one-sided at ends."""

    spacing = float(x[1] - x[0])
    out = np.empty_like(y)
    out[1:-1] = (y[2:] - y[:-2]) / (2.0 * spacing)
    out[0] = (-3.0 * y[0] + 4.0 * y[1] - y[2]) / (2.0 * spacing)
    out[-1] = (3.0 * y[-1] - 4.0 * y[-2] + y[-3]) / (2.0 * spacing)
    return out


def gradient_length_scale(
    grid: AxisymmetricGrid,
    field: npt.ArrayLike,
    *,
    threshold_fraction: float = 0.1,
) -> dict[str, Any]:
    """Return the smallest pointwise ``|f| / |grad f|`` over the active region.

    The active region is ``|f| >= threshold_fraction * max|f|``; outside it the
    ratio is dominated by the field's own zeros and carries no information
    about the structure that has to be resolved.  Points where the gradient
    vanishes are excluded for the same reason (their length scale is infinite,
    which is not a resolution constraint).

    The returned ``points_per_gradient_scale_r`` and
    ``points_per_gradient_scale_z`` are the minimum length scale divided by
    ``dr`` and ``dz``.  They are indicative only, for two reasons.  First, the
    differences used to build the gradient are themselves least accurate
    exactly where the length scale is smallest, so a small value is a warning
    and a large value is not a guarantee.  Second -- and this is easy to
    misread -- the minimum is frequently attained *at the edge of the active
    region*, where the field is passing through ``threshold_fraction`` of its
    peak, rather than at the narrowest feature: a field that decays smoothly
    to zero has a vanishing ``|f|/|grad f|`` wherever it vanishes.  The
    threshold is therefore part of the definition and is returned with the
    number; the reported location must be inspected before the minimum is
    interpreted as the width of anything.  :func:`peak_resolution_report` uses
    the FWHM and front measurements, not this one, for its points-per-scale
    counts.
    """

    grid = _require_grid(grid)
    _reject_complex(field, name="field")
    values = grid.validate_field(field, name="field")
    fraction = _validate_positive_fraction(
        threshold_fraction, name="threshold_fraction"
    )

    magnitude = np.abs(values)
    peak = float(np.max(magnitude))
    if not peak > 0.0:
        raise ValueError("field is identically zero; no length scale is defined")

    d_r = _radial_derivative(grid, values)
    d_z = _axial_derivative(grid, values)
    gradient_magnitude = np.sqrt(d_r * d_r + d_z * d_z)

    threshold = fraction * peak
    active = (magnitude >= threshold) & (gradient_magnitude > 0.0)
    if not np.any(active):
        raise ValueError(
            "no sample has both a significant magnitude and a nonzero gradient"
        )

    scale = np.full(grid.shape, np.inf, dtype=np.float64)
    np.divide(magnitude, gradient_magnitude, out=scale, where=active)
    index_r, index_z = np.unravel_index(int(np.argmin(scale)), grid.shape)
    minimum = float(scale[index_r, index_z])
    if not np.isfinite(minimum):  # pragma: no cover - excluded by ``active``
        raise ArithmeticError("minimum gradient length scale is not finite")

    return {
        "min_length_scale": minimum,
        "r": float(grid.r[index_r]),
        "z": float(grid.z[index_z]),
        "index_r": int(index_r),
        "index_z": int(index_z),
        "points_per_gradient_scale_r": minimum / grid.dr,
        "points_per_gradient_scale_z": minimum / grid.dz,
        "threshold_fraction": fraction,
        "threshold_value": float(threshold),
        "active_samples": int(np.count_nonzero(active)),
        "interpretation": (
            "|f|/|grad f| from three-point differences over the region "
            "|f| >= threshold_fraction * max|f|; a resolution indicator, not "
            "a continuum length scale"
        ),
    }


# --------------------------------------------------------------------------
# Discrete high-frequency content
# --------------------------------------------------------------------------


def _tail_cut(mode_count: int) -> int:
    """First index of the top third of ``mode_count`` non-negative modes."""

    if mode_count < 2:
        raise ValueError("at least two modes are required")
    highest = mode_count - 1
    return int(math.ceil(2.0 * highest / 3.0))


def _real_mode_weights(n_samples: int, mode_count: int) -> FloatArray:
    """Parseval weights for :func:`numpy.fft.rfft` of a real signal."""

    weights = np.ones(mode_count, dtype=np.float64)
    weights[1 : (n_samples + 1) // 2] = 2.0
    return weights


def _dct_ii(profile: FloatArray) -> FloatArray:
    """Orthonormal DCT-II of ``profile``, built from an even extension.

    ``X_k = sum_n x_n cos(pi (n + 1/2) k / N)`` is obtained from the length
    ``2N`` even extension ``v = [x_0..x_{N-1}, x_{N-1}..x_0]`` through

    ``X_k = 0.5 * Re( exp(-i pi k / (2N)) * FFT(v)_k )``,

    which follows from pairing the terms ``n`` and ``2N-1-n`` of the transform.
    The returned coefficients carry the orthonormal scaling
    (``sqrt(1/N)`` for ``k = 0`` and ``sqrt(2/N)`` otherwise), so Parseval's
    identity ``sum_k X_k^2 = sum_n x_n^2`` holds and energy fractions are
    well defined.  SciPy is deliberately not used; this package depends on
    NumPy alone.
    """

    n = int(profile.size)
    if n < 2:
        raise ValueError("at least two samples are required for a DCT")
    extended = np.concatenate((profile, profile[::-1]))
    spectrum = np.fft.rfft(extended)[:n]
    phase = np.exp(-1j * np.pi * np.arange(n) / (2.0 * n))
    coefficients = 0.5 * np.real(phase * spectrum)
    scale = np.full(n, math.sqrt(2.0 / n), dtype=np.float64)
    scale[0] = math.sqrt(1.0 / n)
    return coefficients * scale


def high_frequency_tail(
    grid: AxisymmetricGrid,
    field: npt.ArrayLike,
) -> dict[str, Any]:
    """Return the discrete energy fraction in the top third of the modes.

    Two independent indicators are computed.

    ``z``:
        every radial row is transformed with :func:`numpy.fft.rfft` (the axial
        direction is periodic, so this is the natural basis), the energy is
        summed with the Parseval weights that make the total equal the sum of
        squared samples, and the fraction carried by the top third of the
        non-negative wavenumbers is recorded.  The reported value is the
        maximum over rows together with the radius that attains it.
    ``r``:
        the radial profile through the axial index of the global ``|field|``
        peak is transformed with a DCT-II (:func:`_dct_ii`, implemented on an
        even extension with :mod:`numpy.fft`; the radial direction is not
        periodic, so a Fourier basis would manufacture a jump at ``r_max``),
        and the same top-third fraction is recorded.

    Both numbers are *resolution indicators for the sampled array*, not spectra
    of the continuum solution.  The discrete transform can only see the modes
    the grid represents: continuum content above the Nyquist wavenumber is
    aliased into the represented modes rather than appearing in the tail.  A
    large tail fraction is therefore evidence that the data is under-resolved,
    while a small one is not evidence that it is resolved.  The total energy
    includes the mean (mode zero), so a field dominated by a large mean will
    report a small fraction by construction.
    """

    grid = _require_grid(grid)
    _reject_complex(field, name="field")
    values = grid.validate_field(field, name="field")

    magnitude = np.abs(values)
    peak = float(np.max(magnitude))
    if not peak > 0.0:
        raise ValueError("field is identically zero; no tail fraction is defined")
    index_r, index_z = np.unravel_index(int(np.argmax(magnitude)), grid.shape)

    spectrum = np.abs(np.fft.rfft(values, axis=1)) ** 2
    mode_count_z = spectrum.shape[1]
    weights = _real_mode_weights(grid.nz, mode_count_z)
    energy = spectrum * weights[None, :]
    total = energy.sum(axis=1)
    cut_z = _tail_cut(mode_count_z)
    tail = energy[:, cut_z:].sum(axis=1)
    fraction_z = np.zeros(grid.nr, dtype=np.float64)
    np.divide(tail, total, out=fraction_z, where=total > 0.0)
    row = int(np.argmax(fraction_z))

    profile = values[:, index_z]
    coefficients = _dct_ii(profile)
    cut_r = _tail_cut(coefficients.size)
    energy_r = coefficients**2
    total_r = float(energy_r.sum())
    fraction_r = (
        float(energy_r[cut_r:].sum()) / total_r if total_r > 0.0 else 0.0
    )

    return {
        "z_tail_fraction": float(fraction_z[row]),
        "z_tail_fraction_r": float(grid.r[row]),
        "z_tail_fraction_index_r": row,
        "z_mode_count": int(mode_count_z),
        "z_cut_index": int(cut_z),
        "r_tail_fraction": float(fraction_r),
        "r_profile_index_z": int(index_z),
        "r_profile_z": float(grid.z[index_z]),
        "r_mode_count": int(coefficients.size),
        "r_cut_index": int(cut_r),
        "peak_index_r": int(index_r),
        "interpretation": (
            "discrete energy fractions in the top third of the representable "
            "modes; resolution indicators for the sampled array, not spectra "
            "of the continuum solution (aliasing folds unrepresented content "
            "back into the represented modes)"
        ),
    }


# --------------------------------------------------------------------------
# Assembled report and the fit gate
# --------------------------------------------------------------------------


def _axial_profile(
    grid: AxisymmetricGrid,
    magnitude: FloatArray,
    index_r: int,
    index_z: int,
) -> tuple[FloatArray, FloatArray, int]:
    """Return the axial profile through a peak, unwrapped when ``z`` is periodic.

    The axial direction omits its right endpoint, so a peak near index ``0`` is
    interior on the circle even though it sits at the edge of the stored array.
    Rolling the profile so that the peak lands in the middle and giving the
    samples the continued coordinates ``z_peak + (m - n_z//2) dz`` measures the
    physical width instead of reporting a spurious truncation.
    """

    row = magnitude[index_r, :]
    if not grid.periodic_z:
        return np.asarray(grid.z, dtype=np.float64), row, int(index_z)
    center = grid.nz // 2
    shift = center - int(index_z)
    rolled = np.roll(row, shift)
    offsets = (np.arange(grid.nz, dtype=np.float64) - center) * grid.dz
    return float(grid.z[index_z]) + offsets, rolled, center


def peak_resolution_report(
    grid: AxisymmetricGrid,
    field: npt.ArrayLike,
) -> dict[str, Any]:
    """Assemble every resolution diagnostic for the global ``|field|`` peak.

    The radial profile is taken at the peak's axial index and the axial profile
    at its radial index (unwrapped periodically, see :func:`_axial_profile`).
    The radial front is measured on the flank facing the axis, which is the
    steep side of the structures the shipped Hou snapshots contain.

    Returns
    -------
    dict
        Nested results from :func:`fwhm_1d`, :func:`front_thickness_10_90`,
        :func:`subgrid_peak`, :func:`gradient_length_scale` and
        :func:`high_frequency_tail`, together with the derived counts
        ``points_per_fwhm_r``, ``points_per_fwhm_z``, ``points_per_front`` and
        ``peak_to_axis_cells``, and the collected truncation flags.  Every
        truncation flag raised by any sub-measurement is propagated into
        ``truncation_flags`` and summarized by ``any_truncation``; a report
        with ``any_truncation`` true describes lower bounds, not widths.
    """

    grid = _require_grid(grid)
    _reject_complex(field, name="field")
    values = grid.validate_field(field, name="field")
    magnitude = np.abs(values)
    peak_value = float(np.max(magnitude))
    if not peak_value > 0.0:
        raise ValueError("field is identically zero; no peak to report")

    index_r, index_z = np.unravel_index(int(np.argmax(magnitude)), grid.shape)
    index_r = int(index_r)
    index_z = int(index_z)

    radial_profile = magnitude[:, index_z]
    axial_coordinates, axial_profile, axial_peak = _axial_profile(
        grid, magnitude, index_r, index_z
    )

    fwhm_r = fwhm_1d(grid.r, radial_profile, peak_index=index_r)
    fwhm_z = fwhm_1d(axial_coordinates, axial_profile, peak_index=axial_peak)
    front = front_thickness_10_90(
        grid.r,
        radial_profile,
        peak_index=index_r,
        side="decreasing_coordinate",
    )
    subgrid_r = subgrid_peak(grid.r, radial_profile, index_r)
    subgrid_z = subgrid_peak(axial_coordinates, axial_profile, axial_peak)
    gradient = gradient_length_scale(grid, values)
    tail = high_frequency_tail(grid, values)

    truncation_flags = {
        "fwhm_r_left": bool(fwhm_r["truncated_left"]),
        "fwhm_r_right": bool(fwhm_r["truncated_right"]),
        "fwhm_z_left": bool(fwhm_z["truncated_left"]),
        "fwhm_z_right": bool(fwhm_z["truncated_right"]),
        "front_90": bool(front["truncated_90"]),
        "front_10": bool(front["truncated_10"]),
    }

    return {
        "peak": {
            "index_r": index_r,
            "index_z": index_z,
            "r": float(grid.r[index_r]),
            "z": float(grid.z[index_z]),
            "value": float(values[index_r, index_z]),
            "magnitude": peak_value,
        },
        "subgrid_peak_r": subgrid_r,
        "subgrid_peak_z": subgrid_z,
        "fwhm_r": fwhm_r,
        "fwhm_z": fwhm_z,
        "points_per_fwhm_r": float(fwhm_r["fwhm"] / grid.dr),
        "points_per_fwhm_z": float(fwhm_z["fwhm"] / grid.dz),
        "peak_to_axis_cells": index_r,
        "front_to_axis": front,
        "points_per_front": float(front["thickness"] / grid.dr),
        "gradient_length_scale": gradient,
        "high_frequency_tail": tail,
        "truncation_flags": truncation_flags,
        "any_truncation": bool(any(truncation_flags.values())),
        "nr": grid.nr,
        "nz": grid.nz,
        "dr": grid.dr,
        "dz": grid.dz,
        "interpretation": (
            "points-per-scale counts for the sampled snapshot; they bound what "
            "the grid can represent and are not statements about the continuum "
            "solution or about singularity formation"
        ),
    }


def fit_precondition(
    report: dict[str, Any],
    *,
    min_points: int = PREREGISTERED_MIN_POINTS_PER_FRONT,
) -> dict[str, Any]:
    """Decide whether a report is resolved enough to permit a convergence fit.

    All four criteria must hold:

    1. at least ``min_points`` grid points across the 10-90 front facing the
       axis;
    2. at least ``min_points`` grid points across the radial full width at half
       maximum;
    3. no truncation flag anywhere in the report -- a truncated measurement is
       a lower bound, and a gate cannot be passed on a lower bound;
    4. the peak sits at least ``min_points`` cells from the axis, so that the
       flank being measured is not the axis boundary itself.

    The returned ``satisfied`` is ``False`` whenever any criterion fails, and
    ``reasons`` names each failure.  Passing this gate does not certify a fit;
    it only removes the preregistered prohibition on attempting one.
    """

    if not isinstance(report, dict):
        raise TypeError("report must be a dict produced by peak_resolution_report")
    if isinstance(min_points, bool) or not isinstance(
        min_points, (int, np.integer)
    ):
        raise TypeError("min_points must be an integer")
    min_points = int(min_points)
    if min_points < 1:
        raise ValueError("min_points must be at least one")

    required = (
        "points_per_front",
        "points_per_fwhm_r",
        "peak_to_axis_cells",
        "truncation_flags",
    )
    missing = [key for key in required if key not in report]
    if missing:
        raise ValueError(
            "report is missing required keys: " + ", ".join(sorted(missing))
        )
    flags = report["truncation_flags"]
    if not isinstance(flags, dict):
        raise ValueError("report['truncation_flags'] must be a dict of flags")

    points_per_front = float(report["points_per_front"])
    points_per_fwhm_r = float(report["points_per_fwhm_r"])
    peak_to_axis_cells = int(report["peak_to_axis_cells"])
    if not np.isfinite(points_per_front) or not np.isfinite(points_per_fwhm_r):
        raise ValueError("report contains non-finite points-per-scale counts")

    reasons: list[str] = []
    if points_per_front < min_points:
        reasons.append(
            f"points per 10-90 front thickness {points_per_front:.3f} is below "
            f"the preregistered minimum {min_points}"
        )
    if points_per_fwhm_r < min_points:
        reasons.append(
            f"points per radial FWHM {points_per_fwhm_r:.3f} is below the "
            f"preregistered minimum {min_points}"
        )
    for name in sorted(flags):
        if flags[name]:
            reasons.append(
                f"measurement {name} was truncated at the sampled boundary, so "
                f"it is a lower bound rather than a width"
            )
    if peak_to_axis_cells < min_points:
        reasons.append(
            f"the peak sits {peak_to_axis_cells} cells from the axis, below "
            f"the preregistered minimum {min_points}"
        )

    return {
        "points_per_front": points_per_front,
        "points_per_fwhm_r": points_per_fwhm_r,
        "peak_to_axis_cells": peak_to_axis_cells,
        "min_points": min_points,
        "any_truncation": bool(any(bool(flags[name]) for name in flags)),
        "satisfied": not reasons,
        "reasons": reasons,
        "interpretation": (
            "preregistered resolution gate for convergence fits; satisfied "
            "removes the prohibition on fitting, it does not validate a fit"
        ),
    }


# --------------------------------------------------------------------------
# Common-grid comparison
# --------------------------------------------------------------------------


def _catmull_rom_weights(t: FloatArray) -> FloatArray:
    """Cubic Catmull-Rom (Keys ``a = -1/2``) weights for parameters ``t``."""

    t2 = t * t
    t3 = t2 * t
    return 0.5 * np.stack(
        (
            -t3 + 2.0 * t2 - t,
            3.0 * t3 - 5.0 * t2 + 2.0,
            -3.0 * t3 + 4.0 * t2 + t,
            t3 - t2,
        ),
        axis=-1,
    )


def _stencil(
    targets: FloatArray,
    origin: float,
    spacing: float,
    count: int,
    *,
    periodic: bool,
) -> tuple[np.ndarray, FloatArray]:
    """Return four-point stencil indices and Catmull-Rom weights.

    Non-periodic axes clamp both the parameter and the stencil indices at the
    boundaries; periodic axes wrap.  Targets that coincide with a sample (to
    within :data:`_SNAP_TOLERANCE` of a cell) are snapped onto it, so restricting
    a field to its own grid is the identity rather than a roundoff-level
    perturbation of it.
    """

    s = (np.asarray(targets, dtype=np.float64) - origin) / spacing
    if not periodic:
        s = np.clip(s, 0.0, float(count - 1))
    base = np.floor(s)
    t = s - base
    high = t > 1.0 - _SNAP_TOLERANCE
    base = np.where(high, base + 1.0, base)
    t = np.where(high, 0.0, t)
    t = np.where(t < _SNAP_TOLERANCE, 0.0, t)

    offsets = np.arange(-1, 3, dtype=np.float64)
    raw = base[:, None] + offsets[None, :]
    if periodic:
        indices = np.mod(raw.astype(np.int64), count)
    else:
        indices = np.clip(raw.astype(np.int64), 0, count - 1)
    return indices, _catmull_rom_weights(t)


def _restrict(
    source: AxisymmetricGrid,
    field: FloatArray,
    target_r: FloatArray,
    target_z: FloatArray,
) -> FloatArray:
    """Separable Catmull-Rom interpolation of ``field`` onto a tensor target."""

    index_r, weight_r = _stencil(
        target_r, float(source.r[0]), source.dr, source.nr, periodic=False
    )
    partial = np.einsum("ik,ikj->ij", weight_r, field[index_r, :])
    index_z, weight_z = _stencil(
        target_z, float(source.z[0]), source.dz, source.nz, periodic=True
    )
    return np.einsum("jk,ijk->ij", weight_z, partial[:, index_z])


def interpolate_to_common_grid(
    coarse_grid: AxisymmetricGrid,
    coarse_field: npt.ArrayLike,
    fine_grid: AxisymmetricGrid,
    fine_field: npt.ArrayLike,
) -> dict[str, Any]:
    """Compare two fields on the coarse grid's physical coordinates.

    The common grid is the *coarse* one: the fine field is restricted onto it
    by separable Catmull-Rom cubic interpolation (third-order accurate, ``C^1``,
    interpolating), with periodic wraparound in ``z`` and index clamping at the
    radial boundaries.  The coarse field is passed through the same routine, so
    both sides of the comparison are produced by the same code path; because
    its targets coincide with its own samples the operation is the identity.

    Restricting the fine field to the coarse grid is a deliberate choice.
    Prolonging the coarse field to the fine grid would compare the fine data
    against an interpolant of the coarse data at points the coarse grid cannot
    represent, and the interpolation error would then be reported as a
    difference between the two solutions.

    Returns
    -------
    dict
        ``linf`` (maximum absolute difference), ``weighted_l2``
        (``sqrt(dr dz sum w_i r_i diff^2)`` with the trapezoid radial weight
        ``w`` halved at the first and last radial index), ``derivative_linf``
        (the larger of the radial and axial centered-derivative differences,
        both computed on the common grid), and ``argmax_r`` / ``argmax_z``, the
        coordinates where ``linf`` is attained.

    Raises
    ------
    ValueError
        If the two grids do not share ``r_max`` and the ``z`` period, or if
        either grid is not periodic in ``z``.
    """

    coarse_grid = _require_grid(coarse_grid, name="coarse_grid")
    fine_grid = _require_grid(fine_grid, name="fine_grid")
    if not coarse_grid.periodic_z or not fine_grid.periodic_z:
        raise ValueError(
            "interpolate_to_common_grid requires periodic_z=True on both grids"
        )
    coarse_r_max = float(coarse_grid.r[-1])
    fine_r_max = float(fine_grid.r[-1])
    if not np.isclose(coarse_r_max, fine_r_max, rtol=1.0e-12, atol=0.0):
        raise ValueError(
            f"grids must share r_max, got {coarse_r_max} and {fine_r_max}"
        )
    coarse_period = float(coarse_grid.z_period)
    fine_period = float(fine_grid.z_period)
    if not np.isclose(coarse_period, fine_period, rtol=1.0e-12, atol=0.0):
        raise ValueError(
            f"grids must share the z period, got {coarse_period} and "
            f"{fine_period}"
        )

    _reject_complex(coarse_field, name="coarse_field")
    _reject_complex(fine_field, name="fine_field")
    coarse = coarse_grid.validate_field(coarse_field, name="coarse_field")
    fine = fine_grid.validate_field(fine_field, name="fine_field")

    target_r = np.asarray(coarse_grid.r, dtype=np.float64)
    target_z = np.asarray(coarse_grid.z, dtype=np.float64)
    coarse_common = _restrict(coarse_grid, coarse, target_r, target_z)
    fine_common = _restrict(fine_grid, fine, target_r, target_z)

    difference = coarse_common - fine_common
    absolute = np.abs(difference)
    index_r, index_z = np.unravel_index(
        int(np.argmax(absolute)), coarse_grid.shape
    )

    radial_weight = np.ones(coarse_grid.nr, dtype=np.float64)
    radial_weight[0] = 0.5
    radial_weight[-1] = 0.5
    weighted = float(
        np.sum(
            radial_weight[:, None]
            * target_r[:, None]
            * difference
            * difference
        )
        * coarse_grid.dr
        * coarse_grid.dz
    )
    if weighted < 0.0:  # pragma: no cover - the summand is non-negative
        raise ArithmeticError("negative weighted squared difference")

    derivative_r = float(
        np.max(
            np.abs(
                _radial_derivative(coarse_grid, coarse_common)
                - _radial_derivative(coarse_grid, fine_common)
            )
        )
    )
    derivative_z = float(
        np.max(
            np.abs(
                _axial_derivative(coarse_grid, coarse_common)
                - _axial_derivative(coarse_grid, fine_common)
            )
        )
    )

    return {
        "linf": float(np.max(absolute)),
        "weighted_l2": float(np.sqrt(weighted)),
        "derivative_linf": max(derivative_r, derivative_z),
        "derivative_linf_r": derivative_r,
        "derivative_linf_z": derivative_z,
        "argmax_r": float(coarse_grid.r[index_r]),
        "argmax_z": float(coarse_grid.z[index_z]),
        "argmax_index_r": int(index_r),
        "argmax_index_z": int(index_z),
        "common_shape": coarse_grid.shape,
        "interpolation": (
            "separable Catmull-Rom cubic (Keys a=-1/2), periodic in z, "
            "clamped in r; the fine field is restricted onto the coarse grid"
        ),
    }


# --------------------------------------------------------------------------
# Preregistration study
# --------------------------------------------------------------------------


def _tanh_front(radius: FloatArray, center: float, delta: float) -> FloatArray:
    """Return ``0.5 (1 + tanh((center - r) / delta))``.

    The normalization matters: the profile falls from one to zero, so a
    fraction *of the peak value* is also a fraction of the full swing, and the
    10-90 thickness measured by :func:`front_thickness_10_90` is the classical
    10-90 rise distance whose exact value is ``2 atanh(0.8) delta``.
    """

    return 0.5 * (1.0 + np.tanh((center - radius) / delta))


def manufactured_front_study(
    *,
    delta_over_dr_values: Sequence[float] = PREREGISTERED_DELTA_OVER_DR_VALUES,
    n_z: int = 64,
) -> dict[str, Any]:
    r"""Measure width-diagnostic error against an exactly known front.

    The manufactured field is the ``tanh`` front
    ``f(r) = 0.5 (1 + tanh((r0 - r) / delta))`` on ``r`` in ``[0, 1]`` with
    ``delta = 0.05``, sampled on a uniform radial grid whose spacing is set by
    each requested ``delta / dr``.  Two continuum quantities are known in
    closed form and are recomputed here rather than quoted:

    * the 10%-to-90% front thickness is
      ``delta (atanh(0.8) - atanh(-0.8)) = 2 atanh(0.8) delta``;
    * the derivative ``f'(r) = -sech^2((r0 - r)/delta) / (2 delta)`` has full
      width at half maximum ``2 arccosh(sqrt(2)) delta``.

    For each resolution the front centre is placed at
    :data:`PREREGISTERED_FRONT_PHASES` equally spaced sub-cell offsets and the
    reported error is the **worst case** over those placements.  Without this,
    the measured error would oscillate with the accidental alignment of the
    front with the mesh and a threshold read off a single alignment would not
    transfer to data.

    The derivative profile is formed by the same three-point centered
    differences the rest of this module uses, so the reported ``fwhm_error``
    includes the differencing error a caller would actually incur; it is not
    an analytic derivative sampled on the grid.

    ``n_z`` only sets the axial size of the
    :class:`~ns_certificate_lab.grid.AxisymmetricGrid` that carries the radial
    coordinates.  The study is one-dimensional in ``r`` and no measured
    quantity depends on it; it exists so that the study samples the same grid
    object the rest of the package uses rather than a private coordinate
    array.

    Returns
    -------
    dict
        ``table`` maps each requested ``delta / dr`` to a dict with
        ``front_error`` and ``fwhm_error`` (worst-case relative errors), the
        measured and exact thicknesses, and ``points_per_front`` (the measured
        10-90 thickness in grid points).  The remaining top-level keys record
        the constants the study was run with, so the table can be audited
        without re-reading this source.

    Notes
    -----
    This study contains no Hou data and imports nothing from the experiment
    drivers.  It is the sole justification for
    :data:`PREREGISTERED_MIN_POINTS_PER_FRONT`.
    """

    if isinstance(n_z, bool) or not isinstance(n_z, (int, np.integer)):
        raise TypeError("n_z must be an integer")
    n_z = int(n_z)
    if n_z < 5:
        raise ValueError("n_z must be at least five")
    values = list(delta_over_dr_values)
    if not values:
        raise ValueError("delta_over_dr_values must not be empty")

    delta = _STUDY_DELTA
    exact_front = 2.0 * math.atanh(0.8) * delta
    exact_fwhm = 2.0 * math.acosh(math.sqrt(2.0)) * delta

    table: dict[float, dict[str, Any]] = {}
    for requested in values:
        if isinstance(requested, bool) or not isinstance(
            requested, (int, float, np.integer, np.floating)
        ):
            raise TypeError("delta_over_dr_values must contain real numbers")
        ratio = float(requested)
        if not np.isfinite(ratio) or ratio <= 0.0:
            raise ValueError("delta_over_dr_values must be positive and finite")

        n_r = int(round(ratio / delta)) + 1
        if n_r < 4:
            raise ValueError(
                f"delta/dr = {ratio} gives only {n_r} radial points; the study "
                "needs at least four"
            )
        grid = AxisymmetricGrid.uniform(
            nr=n_r, nz=n_z, r_max=1.0, z_min=0.0, z_max=1.0, periodic_z=True
        )
        radius = np.asarray(grid.r, dtype=np.float64)
        dr = grid.dr

        front_errors: list[float] = []
        fwhm_errors: list[float] = []
        thicknesses: list[float] = []
        widths: list[float] = []
        truncated = False
        for phase in range(PREREGISTERED_FRONT_PHASES):
            center = _STUDY_CENTER + phase * dr / PREREGISTERED_FRONT_PHASES
            profile = _tanh_front(radius, center, delta)
            # The closed-form references assume the front is fully developed
            # inside the domain; verify rather than assume.
            saturated = abs(float(profile[0]) - 1.0) > 1.0e-6
            decayed = abs(float(profile[-1])) > 1.0e-6
            if saturated or decayed:
                raise ValueError(
                    "the manufactured front is not fully contained in [0, 1]; "
                    "the closed-form references would not apply"
                )
            measured = front_thickness_10_90(
                radius, profile, peak_index=0, side="increasing_coordinate"
            )
            derivative = _derivative_1d(radius, profile)
            width = fwhm_1d(radius, derivative)
            truncated = truncated or bool(measured["truncated"])
            truncated = (
                truncated
                or bool(width["truncated_left"])
                or bool(width["truncated_right"])
            )
            thicknesses.append(float(measured["thickness"]))
            widths.append(float(width["fwhm"]))
            front_errors.append(
                abs(float(measured["thickness"]) - exact_front) / exact_front
            )
            fwhm_errors.append(abs(float(width["fwhm"]) - exact_fwhm) / exact_fwhm)

        worst = int(np.argmax(front_errors))
        table[ratio] = {
            "delta_over_dr_requested": ratio,
            "delta_over_dr_realized": float(delta / dr),
            "n_r": int(n_r),
            "dr": float(dr),
            "front_thickness": thicknesses[worst],
            "front_thickness_exact": exact_front,
            "front_error": float(front_errors[worst]),
            "points_per_front": float(thicknesses[worst] / dr),
            "points_per_front_exact": float(exact_front / dr),
            "derivative_fwhm": widths[int(np.argmax(fwhm_errors))],
            "derivative_fwhm_exact": exact_fwhm,
            "fwhm_error": float(max(fwhm_errors)),
            "truncated": bool(truncated),
        }

    return {
        "table": table,
        "delta": delta,
        "center": _STUDY_CENTER,
        "phases": PREREGISTERED_FRONT_PHASES,
        "n_z": n_z,
        "exact_front_thickness": exact_front,
        "exact_derivative_fwhm": exact_fwhm,
        "front_error_tolerance": PREREGISTERED_FRONT_ERROR_TOLERANCE,
        "preregistered_min_points_per_front": PREREGISTERED_MIN_POINTS_PER_FRONT,
        "interpretation": (
            "worst-case-over-phase relative error of the width diagnostics on "
            "a closed-form tanh front; the sole basis of "
            "PREREGISTERED_MIN_POINTS_PER_FRONT and independent of any "
            "simulation data"
        ),
    }
