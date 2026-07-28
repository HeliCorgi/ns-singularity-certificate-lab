"""Verification of the core-width and points-per-scale diagnostics (P0-D).

Every oracle used below is written out in closed form in this file:  the
Gaussian full width at half maximum ``2 sqrt(2 ln 2) sigma``, the ``tanh``
10-90 front thickness ``2 atanh(0.8) delta``, the ``sech^2`` derivative width
``2 arccosh(sqrt(2)) delta``, the vertex of an exact parabola, and the analytic
gradient of a closed-form field.  No production routine is used to derive the
value another production routine is checked against.

The preregistration test (:func:`test_manufactured_front_study_supports_the_
preregistered_threshold`) recomputes the manufactured-front study from scratch
and re-derives ``PREREGISTERED_MIN_POINTS_PER_FRONT`` from the measured table.
If the study ever stops supporting the stored constant this test fails, which
is the point: the constant may not be quietly retuned to admit a data set that
would otherwise be refused.

Measured constants quoted in the assertions were obtained on this machine with
the pinned environment.  The bounds keep margin so that they pin behaviour
rather than arithmetic noise.
"""

from __future__ import annotations

import ast
import inspect
import math
from pathlib import Path

import numpy as np
import pytest

import ns_certificate_lab.core_width as core_width_module
from ns_certificate_lab.core_width import (
    PREREGISTERED_DELTA_OVER_DR_VALUES,
    PREREGISTERED_FRONT_ERROR_TOLERANCE,
    PREREGISTERED_MIN_POINTS_PER_FRONT,
    _dct_ii,
    fit_precondition,
    front_thickness_10_90,
    fwhm_1d,
    gradient_length_scale,
    high_frequency_tail,
    interpolate_to_common_grid,
    manufactured_front_study,
    peak_resolution_report,
    subgrid_peak,
)
from ns_certificate_lab.grid import AxisymmetricGrid

GAUSSIAN_FWHM_FACTOR = 2.0 * math.sqrt(2.0 * math.log(2.0))
"""Full width at half maximum of ``exp(-x^2 / (2 sigma^2))`` in units of sigma."""

TANH_FRONT_FACTOR = 2.0 * math.atanh(0.8)
"""10-90 thickness of ``0.5 (1 + tanh(-x/delta))`` in units of delta."""

SECH2_FWHM_FACTOR = 2.0 * math.acosh(math.sqrt(2.0))
"""Full width at half maximum of ``sech^2(x/delta)`` in units of delta."""


def _grid(nr: int, nz: int, *, r_max: float = 1.0, z_max: float = 1.0):
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=r_max, z_min=0.0, z_max=z_max, periodic_z=True
    )


def _tanh_front(radius: np.ndarray, center: float, delta: float) -> np.ndarray:
    """``0.5 (1 + tanh((center - r)/delta))``: falls from one to zero at r=center.

    Written out here independently of the production helper of the same shape.
    With this normalization a fraction of the peak value is also a fraction of
    the full swing, so the 10-90 thickness measured by
    ``front_thickness_10_90`` is the classical rise distance
    ``delta (atanh(0.8) - atanh(-0.8)) = 2 atanh(0.8) delta``.
    """

    return 0.5 * (1.0 + np.tanh((center - radius) / delta))


# --------------------------------------------------------------------------
# 1-2.  Full width at half maximum
# --------------------------------------------------------------------------


def test_fwhm_matches_a_finely_sampled_gaussian() -> None:
    """FWHM of a Gaussian must reproduce ``2 sqrt(2 ln 2) sigma``.

    Measured on this machine: relative error 3.1e-07 with 4001 samples over
    twelve standard deviations, both truncation flags clear.  The requirement
    is 0.5%; the linear-interpolation error is four orders of magnitude below
    that, so this bound is not the thing under test -- the flags and the
    factor are.
    """

    sigma = 0.7
    x = np.linspace(-6.0 * sigma, 6.0 * sigma, 4001)
    y = np.exp(-0.5 * (x / sigma) ** 2)

    result = fwhm_1d(x, y)
    exact = GAUSSIAN_FWHM_FACTOR * sigma

    assert abs(result["fwhm"] - exact) <= 5.0e-3 * exact
    assert abs(result["fwhm"] - exact) <= 1.0e-5 * exact
    assert result["truncated_left"] is False
    assert result["truncated_right"] is False
    assert result["peak_value"] == pytest.approx(1.0, rel=1.0e-12)
    assert abs(result["peak_coordinate"]) <= 1.0e-12
    # The crossings straddle the peak symmetrically for a symmetric profile.
    assert result["left"] == pytest.approx(-0.5 * exact, abs=1.0e-6)
    assert result["right"] == pytest.approx(0.5 * exact, abs=1.0e-6)
    assert result["half_maximum"] == pytest.approx(0.5, rel=1.0e-12)

    # A sign flip must not change the width: the measurement is on |values|.
    flipped = fwhm_1d(x, -y)
    assert flipped["fwhm"] == result["fwhm"]
    assert flipped["signed_peak_value"] == -1.0 * result["signed_peak_value"]


def test_fwhm_flags_truncation_at_the_array_end() -> None:
    """A peak too close to the end must be reported as truncated, not narrow.

    The half-maximum crossing on the right lies outside the sampled interval,
    so the returned width is a lower bound.  Measured: fwhm 0.0789 against the
    true 0.1177, i.e. a silent answer would have understated the width by 33%.
    Silence is exactly what this module must never do.
    """

    sigma = 0.05
    x = np.linspace(0.0, 1.0, 201)
    y = np.exp(-0.5 * ((x - 0.98) / sigma) ** 2)

    result = fwhm_1d(x, y)
    exact = GAUSSIAN_FWHM_FACTOR * sigma

    assert result["truncated_right"] is True
    assert result["truncated_left"] is False
    assert result["right"] == float(x[-1])
    assert result["fwhm"] < exact
    # The left half is measured correctly even though the right one is not.
    assert result["peak_coordinate"] - result["left"] == pytest.approx(
        0.5 * exact, rel=1.0e-3
    )

    # Both sides truncated: a profile that never falls to half maximum at all.
    # (``2 + x`` peaks at 3 with a half maximum of 1.5, which the profile never
    # reaches; ``1 + x`` would touch its half maximum exactly at the left end
    # and is therefore not the degenerate case wanted here.)
    monotone = fwhm_1d(x, 2.0 + x)
    assert monotone["truncated_left"] is True
    assert monotone["truncated_right"] is True
    assert monotone["fwhm"] == pytest.approx(1.0, rel=1.0e-12)


# --------------------------------------------------------------------------
# 3.  10-90 front thickness
# --------------------------------------------------------------------------


def test_front_thickness_matches_the_analytic_tanh_front() -> None:
    """The 10-90 thickness of a resolved ``tanh`` front must be exact to <1%.

    ``delta / dr = 16`` (321 radial samples, ``delta = 0.05``).  Measured
    relative error 6.9e-04 for the mesh-aligned front; the requirement is 1%.
    The crossings land 3.8e-05 from their exact positions.
    """

    delta = 0.05
    n_r = int(round(16.0 / delta)) + 1
    radius = np.linspace(0.0, 1.0, n_r)
    assert radius[1] - radius[0] == pytest.approx(delta / 16.0, rel=1.0e-12)
    profile = _tanh_front(radius, 0.5, delta)

    result = front_thickness_10_90(
        radius, profile, peak_index=0, side="increasing_coordinate"
    )
    exact = TANH_FRONT_FACTOR * delta

    assert abs(result["thickness"] - exact) <= 1.0e-2 * exact
    assert abs(result["thickness"] - exact) <= 2.0e-3 * exact
    assert result["truncated"] is False
    assert result["truncated_90"] is False
    assert result["truncated_10"] is False
    # The crossings sit symmetrically about the front centre.
    assert result["coordinate_90"] == pytest.approx(
        0.5 - 0.5 * exact, abs=2.0e-4
    )
    assert result["coordinate_10"] == pytest.approx(
        0.5 + 0.5 * exact, abs=2.0e-4
    )
    assert result["level_90"] == pytest.approx(0.9 * float(profile[0]), rel=1e-12)
    assert result["level_10"] == pytest.approx(0.1 * float(profile[0]), rel=1e-12)

    # The other side of the same peak never falls, so it must be truncated
    # rather than reported as a zero-thickness front.
    inward = front_thickness_10_90(
        radius, profile, peak_index=0, side="decreasing_coordinate"
    )
    assert inward["truncated"] is True
    assert inward["thickness"] == 0.0


# --------------------------------------------------------------------------
# 4.  Sub-grid peak
# --------------------------------------------------------------------------


def test_subgrid_peak_recovers_an_exact_parabola_vertex() -> None:
    """Three samples of a parabola must return its vertex to roundoff.

    The samples are exact values of ``5 - 3 (x - 0.37)^2`` on a uniform mesh
    whose nodes miss the vertex.  Measured: coordinate error 6.1e-16, value
    error 0.0.
    """

    x = 0.1 * np.arange(7, dtype=np.float64)
    vertex, curvature, height = 0.37, -3.0, 5.0
    y = height + curvature * (x - vertex) ** 2

    discrete = int(np.argmax(np.abs(y)))
    assert discrete == 4  # the nearest node, x = 0.4, is not the vertex

    result = subgrid_peak(x, y, discrete)
    assert result["refined"] is True
    assert abs(result["coordinate"] - vertex) <= 1.0e-14
    assert abs(result["value"] - height) <= 1.0e-13
    assert result["curvature"] == pytest.approx(curvature, rel=1.0e-12)
    assert result["offset_cells"] == pytest.approx(-0.3, rel=1.0e-10)
    assert result["discrete_coordinate"] == float(x[discrete])
    assert result["discrete_value"] == float(y[discrete])

    # A boundary index cannot be refined: there is no bracketing stencil.
    for boundary in (0, x.size - 1):
        edge = subgrid_peak(x, y, boundary)
        assert edge["refined"] is False
        assert edge["coordinate"] == float(x[boundary])
        assert edge["value"] == float(y[boundary])
        assert edge["offset_cells"] == 0.0
        assert "boundary" in edge["reason"]

    # Collinear samples have no vertex; refinement is declined, not invented.
    # The coordinates are exact binary fractions so that the divided
    # differences cancel exactly and the collinear branch is the one taken.
    exact = np.arange(7.0)
    linear = subgrid_peak(exact, 2.0 - 0.5 * exact, 3)
    assert linear["refined"] is False
    assert linear["curvature"] == 0.0
    assert "collinear" in linear["reason"]

    # Nearly collinear samples do not reach that branch -- their curvature is
    # a rounding artefact -- and are caught by the out-of-stencil guard
    # instead.  Measured curvature: -1.2e-14 against a slope of -0.5.
    nearly = subgrid_peak(x, 2.0 - 0.5 * x, 3)
    assert nearly["refined"] is False
    assert abs(nearly["curvature"]) <= 1.0e-12
    assert "outside" in nearly["reason"]

    # A monotone stencil would put the vertex outside the three samples; the
    # discrete value is returned instead of an extrapolated fantasy.
    outside = subgrid_peak(
        np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 1.9]), 1
    )
    assert outside["refined"] is False
    assert "outside" in outside["reason"]


# --------------------------------------------------------------------------
# 5.  Gradient length scale
# --------------------------------------------------------------------------


def test_gradient_length_scale_matches_the_analytic_gradient() -> None:
    """``|f|/|grad f|`` must match the closed-form gradient within 10%.

    ``f = sin(2 pi z) (1 - r^2)`` on ``[0,1] x [0,1)``.  The analytic gradient
    is ``(-2 r sin(2 pi z), 2 pi cos(2 pi z)(1 - r^2))``.  The module's radial
    difference is exact on this quadratic profile and its axial difference
    carries the usual ``sin(2 pi dz)/(2 pi dz)`` factor, so the two agree to
    about 4e-4 relative.  Measured: 0.0236001 against the analytic 0.0235906,
    at ``(r, z) = (0.5625, 0.5234)``, giving 1.510 points per length scale
    radially and 3.021 axially.  Note where that minimum sits: not at the
    peak, but where the field falls through the 10% threshold.  This is the
    caveat documented on ``gradient_length_scale`` and the reason
    ``peak_resolution_report`` does not use this number as its width.
    """

    grid = _grid(65, 128)
    mesh_r, mesh_z = grid.mesh()
    field = np.sin(2.0 * np.pi * mesh_z) * (1.0 - mesh_r**2)

    d_r = -2.0 * mesh_r * np.sin(2.0 * np.pi * mesh_z)
    d_z = 2.0 * np.pi * np.cos(2.0 * np.pi * mesh_z) * (1.0 - mesh_r**2)
    magnitude = np.abs(field)
    gradient = np.hypot(d_r, d_z)
    active = (magnitude >= 0.1 * magnitude.max()) & (gradient > 0.0)
    safe = np.where(gradient > 0.0, gradient, 1.0)
    analytic_scale = np.where(active, magnitude / safe, np.inf)
    expected = float(analytic_scale.min())

    result = gradient_length_scale(grid, field)

    assert result["min_length_scale"] == pytest.approx(expected, rel=0.10)
    assert result["points_per_gradient_scale_r"] == pytest.approx(
        expected / grid.dr, rel=0.10
    )
    assert result["points_per_gradient_scale_z"] == pytest.approx(
        expected / grid.dz, rel=0.10
    )
    # The reported counts are the reported scale divided by the spacings; this
    # is an identity, checked so a future refactor cannot decouple them.
    assert result["points_per_gradient_scale_r"] == pytest.approx(
        result["min_length_scale"] / grid.dr, rel=1.0e-14
    )
    assert result["points_per_gradient_scale_z"] == pytest.approx(
        result["min_length_scale"] / grid.dz, rel=1.0e-14
    )

    # The reported location is sane: it lies on the grid, inside the active
    # region, and the analytic length scale there is the analytic minimum.
    assert 0.0 <= result["r"] <= 1.0
    assert 0.0 <= result["z"] < 1.0
    index_r, index_z = result["index_r"], result["index_z"]
    assert grid.r[index_r] == result["r"]
    assert grid.z[index_z] == result["z"]
    assert active[index_r, index_z]
    assert float(analytic_scale[index_r, index_z]) == pytest.approx(
        expected, rel=0.10
    )
    assert 0 < result["active_samples"] < grid.nr * grid.nz

    # The diagnostic must respond to a genuinely finer structure.  Quadrupling
    # the axial wavenumber shortens every axial scale by four; the measured
    # minimum falls from 0.02360 to 0.00796, a factor of 2.96.  (It is not
    # exactly four because the minimum is attained where the field crosses the
    # 10% threshold rather than at a fixed feature -- see the caveat in
    # ``gradient_length_scale``.)
    finer = np.sin(8.0 * np.pi * mesh_z) * (1.0 - mesh_r**2)
    tighter = gradient_length_scale(grid, finer)
    assert tighter["min_length_scale"] < 0.4 * result["min_length_scale"]

    # The threshold controls how much of the field is inspected.
    wide = gradient_length_scale(grid, field, threshold_fraction=0.9)
    assert wide["active_samples"] < result["active_samples"]
    assert wide["min_length_scale"] >= result["min_length_scale"]


# --------------------------------------------------------------------------
# 6.  High-frequency tail
# --------------------------------------------------------------------------


def test_high_frequency_tail_detects_a_nyquist_perturbation() -> None:
    """A 0.1% Nyquist ripple must raise the tail fraction by orders of magnitude.

    Measured in ``z``: 3.4e-32 for the smooth field, 4.1e-06 with a relative
    1e-3 Nyquist ripple.

    The radial indicator has a much higher floor -- 8.1e-08 for the same
    smooth field -- and the reason is worth recording: the radial profile is
    not periodic, so the DCT of a quadratic radial shape decays only
    algebraically and leaves genuine energy in the top third of the modes.  A
    1e-3 radial ripple lifts that floor by a factor of 14 only; the ripple
    used here is therefore 1e-2, which lifts it to 1.2e-04 (a factor of 1450).
    The asymmetry is a property of the two bases, not of the field.
    """

    grid = _grid(33, 64)
    mesh_r, mesh_z = grid.mesh()
    smooth = np.cos(2.0 * np.pi * mesh_z) * (1.0 - 0.3 * mesh_r**2)

    base = high_frequency_tail(grid, smooth)
    assert base["z_tail_fraction"] < 1.0e-6
    assert base["r_tail_fraction"] < 1.0e-6
    assert base["z_cut_index"] == 22  # top third of modes 0..32
    assert base["z_mode_count"] == grid.nz // 2 + 1
    assert base["r_mode_count"] == grid.nr
    assert base["r_cut_index"] == 22  # top third of modes 0..32

    axial_nyquist = (-1.0) ** np.arange(grid.nz)
    rippled = smooth + 1.0e-3 * axial_nyquist[None, :]
    ripple = high_frequency_tail(grid, rippled)
    assert ripple["z_tail_fraction"] > 1.0e-6
    assert ripple["z_tail_fraction"] > 100.0 * max(
        base["z_tail_fraction"], 1.0e-24
    )

    radial_nyquist = (-1.0) ** np.arange(grid.nr)
    radial_rippled = smooth + 1.0e-2 * radial_nyquist[:, None]
    radial_ripple = high_frequency_tail(grid, radial_rippled)
    assert radial_ripple["r_tail_fraction"] > 1.0e-5
    assert radial_ripple["r_tail_fraction"] > 100.0 * max(
        base["r_tail_fraction"], 1.0e-24
    )
    # The radial ripple leaves the axial indicator alone and vice versa: the
    # two indicators are not measuring the same thing twice.
    assert radial_ripple["z_tail_fraction"] < 1.0e-6
    assert ripple["r_tail_fraction"] < 1.0e-5


def test_dct_matches_a_direct_evaluation_and_conserves_energy() -> None:
    """The NumPy-only DCT-II is checked against its defining sum.

    The energy fractions reported by ``high_frequency_tail`` are meaningful
    only if the transform is orthonormal, so both the coefficients and
    Parseval's identity are verified against an independent ``O(N^2)``
    evaluation of ``X_k = sum_n x_n cos(pi (n + 1/2) k / N)``.  Measured over
    ``n = 4, 17, 64``: coefficient error at most 1.9e-14, Parseval defect at
    most 2.1e-14.
    """

    generator = np.random.default_rng(20260728)
    for n in (4, 17, 64):
        samples = generator.standard_normal(n)
        modes = np.arange(n)[:, None]
        positions = np.arange(n)[None, :]
        direct = (
            samples[None, :]
            * np.cos(np.pi * (positions + 0.5) * modes / n)
        ).sum(axis=1)
        scale = np.full(n, math.sqrt(2.0 / n))
        scale[0] = math.sqrt(1.0 / n)
        direct = direct * scale

        transformed = _dct_ii(samples)
        assert np.max(np.abs(transformed - direct)) <= 1.0e-12
        assert float(np.sum(transformed**2)) == pytest.approx(
            float(np.sum(samples**2)), rel=1.0e-12
        )


# --------------------------------------------------------------------------
# 7.  Common-grid comparison
# --------------------------------------------------------------------------


def _smooth_two_grid_field(grid) -> np.ndarray:
    """A closed-form smooth field, unit amplitude, written out here."""

    mesh_r, mesh_z = grid.mesh()
    angle = 2.0 * np.pi * mesh_z / float(grid.z_period)
    return (1.0 - 0.5 * mesh_r**2 + 0.2 * mesh_r**4) * np.cos(angle) + 0.4 * (
        np.exp(-mesh_r) * np.sin(2.0 * angle)
    )


def test_interpolate_to_common_grid_agrees_and_detects() -> None:
    """The same analytic field on two grids agrees; a perturbed one does not.

    Two grid pairs are used.  The nested pair (65x128 and 129x256) has every
    coarse node coincident with a fine node, so the restriction must be an
    exact node lookup: measured linf 0.0e+00, which also pins the snapping of
    interpolation parameters onto nodes.  The non-nested pair (65x128 and
    97x192) exercises the cubic itself: measured linf 1.9e-07, weighted L2
    8.5e-08 and derivative linf 3.7e-07, all far below the 1e-3 requirement.
    """

    coarse = _grid(65, 128, z_max=2.0 * np.pi)
    fine = _grid(129, 256, z_max=2.0 * np.pi)
    coarse_field = _smooth_two_grid_field(coarse)
    fine_field = _smooth_two_grid_field(fine)
    assert np.max(np.abs(coarse_field)) == pytest.approx(1.0, abs=0.5)

    nested = interpolate_to_common_grid(coarse, coarse_field, fine, fine_field)
    assert nested["linf"] < 1.0e-3
    assert nested["linf"] <= 1.0e-12
    assert nested["weighted_l2"] <= 1.0e-12
    assert nested["derivative_linf"] <= 1.0e-10
    assert nested["common_shape"] == coarse.shape

    odd = _grid(97, 192, z_max=2.0 * np.pi)
    odd_field = _smooth_two_grid_field(odd)
    interpolated = interpolate_to_common_grid(coarse, coarse_field, odd, odd_field)
    # Genuine cubic interpolation happens here: half the coarse radial nodes
    # and two thirds of the axial ones fall between fine samples.
    assert interpolated["linf"] > 0.0
    assert interpolated["linf"] < 1.0e-3
    assert interpolated["weighted_l2"] < 1.0e-3
    assert interpolated["derivative_linf"] < 1.0e-2

    # Detection: a deliberately different fine field must show up at its own
    # amplitude, not be smoothed away by the restriction.
    mesh_r, mesh_z = fine.mesh()
    perturbation = 0.1 * mesh_r * np.sin(
        4.0 * np.pi * mesh_z / float(fine.z_period)
    )
    detected = interpolate_to_common_grid(
        coarse, coarse_field, fine, fine_field + perturbation
    )
    assert 0.09 <= detected["linf"] <= 0.11
    assert detected["weighted_l2"] > 0.01
    assert detected["derivative_linf"] > 0.01
    # The perturbation peaks at r = r_max, and that is where it is reported.
    assert detected["argmax_r"] == pytest.approx(1.0, rel=1.0e-12)
    assert detected["argmax_index_r"] == coarse.nr - 1
    assert coarse.r[detected["argmax_index_r"]] == detected["argmax_r"]
    assert coarse.z[detected["argmax_index_z"]] == detected["argmax_z"]

    # A perturbation confined to the derivative is caught by derivative_linf
    # even though it is small in amplitude: 0.002 sin(16 * 2 pi z / L) has linf
    # 0.0020 but a coarse-grid derivative amplitude of 0.0288, a ratio of 14.4.
    wiggle = 0.002 * np.sin(32.0 * np.pi * mesh_z / float(fine.z_period))
    rough = interpolate_to_common_grid(
        coarse, coarse_field, fine, fine_field + wiggle
    )
    assert rough["linf"] <= 0.0021
    assert rough["derivative_linf"] > 10.0 * rough["linf"]


def test_interpolate_to_common_grid_rejects_incompatible_grids() -> None:
    """Mismatched domains must raise instead of comparing unrelated points."""

    coarse = _grid(65, 128)
    fine = _grid(129, 256)
    coarse_field = _smooth_two_grid_field(coarse)
    fine_field = _smooth_two_grid_field(fine)

    other_radius = _grid(129, 256, r_max=1.5)
    with pytest.raises(ValueError, match="r_max"):
        interpolate_to_common_grid(
            coarse, coarse_field, other_radius, _smooth_two_grid_field(other_radius)
        )

    other_period = _grid(129, 256, z_max=2.0)
    with pytest.raises(ValueError, match="z period"):
        interpolate_to_common_grid(
            coarse, coarse_field, other_period, _smooth_two_grid_field(other_period)
        )

    aperiodic = AxisymmetricGrid.uniform(
        nr=129, nz=256, r_max=1.0, z_min=0.0, z_max=1.0, periodic_z=False
    )
    with pytest.raises(ValueError, match="periodic_z"):
        interpolate_to_common_grid(
            coarse, coarse_field, aperiodic, np.zeros(aperiodic.shape)
        )

    with pytest.raises(ValueError, match="shape"):
        interpolate_to_common_grid(coarse, coarse_field[:-1], fine, fine_field)
    with pytest.raises(TypeError):
        interpolate_to_common_grid("not a grid", coarse_field, fine, fine_field)


# --------------------------------------------------------------------------
# 8.  The preregistration itself
# --------------------------------------------------------------------------


def test_manufactured_front_study_supports_the_preregistered_threshold() -> None:
    """Re-derive ``PREREGISTERED_MIN_POINTS_PER_FRONT`` from the study.

    The constant is fixed in ``core_width.py``.  This test recomputes the
    manufactured-front table and derives the constant again from the measured
    errors, so a future edit of the constant that the study does not support
    fails here.  Derivation: the widest front that misses the 2% tolerance
    spans 6.185 points (``delta/dr = 2.75``, error 2.362e-02) and the
    narrowest that meets it spans 6.708 points (``delta/dr = 3.0``, error
    1.767e-02), so seven is the smallest integer point count at which every
    ladder entry passes.

    Note on monotonicity.  The derivative-FWHM error decreases strictly along
    the whole dense ladder.  The front-thickness error does not, and the
    reason is understood: the two crossings sit ``2 atanh(0.8) delta/dr``
    cells apart, so their sub-cell phases differ by the fractional part of
    that number, and the leading errors of the two crossings therefore add or
    partly cancel depending on ``delta/dr`` itself.  Averaging or worst-casing
    over the *placement* of the front (which the study does, over eight
    sub-cell offsets) cannot remove this.  What is asserted instead is
    strict monotonicity along the doubling sub-ladder together with two-sided
    ``O((dr/delta)^2)`` envelopes on the dense ladder, which is the stronger
    statement anyway: measured coefficients lie in [0.128, 0.182].
    """

    study = manufactured_front_study()
    table = study["table"]
    tolerance = PREREGISTERED_FRONT_ERROR_TOLERANCE

    # The study's own closed-form references, recomputed here.
    delta = study["delta"]
    assert study["exact_front_thickness"] == pytest.approx(
        TANH_FRONT_FACTOR * delta, rel=1.0e-14
    )
    assert study["exact_derivative_fwhm"] == pytest.approx(
        SECH2_FWHM_FACTOR * delta, rel=1.0e-14
    )
    assert set(table) == {float(v) for v in PREREGISTERED_DELTA_OVER_DR_VALUES}

    for ratio, entry in table.items():
        assert entry["delta_over_dr_realized"] == pytest.approx(
            ratio, rel=1.0e-12
        )
        assert entry["truncated"] is False
        assert entry["front_thickness_exact"] == study["exact_front_thickness"]
        assert entry["points_per_front"] == pytest.approx(
            entry["front_thickness"] / entry["dr"], rel=1.0e-14
        )
        # Two-sided second-order envelope in dr/delta.
        assert entry["front_error"] <= 0.25 / ratio**2
        assert entry["front_error"] >= 0.09 / ratio**2

    # Strict monotonicity where it can be asserted honestly.
    dense = [table[float(v)] for v in PREREGISTERED_DELTA_OVER_DR_VALUES]
    for previous, following in zip(dense, dense[1:]):
        assert following["fwhm_error"] < previous["fwhm_error"]

    doubling = manufactured_front_study(
        delta_over_dr_values=(0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
    )["table"]
    ladder = [doubling[v] for v in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)]
    for previous, following in zip(ladder, ladder[1:]):
        assert following["front_error"] < previous["front_error"]
        assert following["fwhm_error"] < previous["fwhm_error"]
    # Second order: halving dr must shrink the error by roughly four.
    for previous, following in zip(ladder[1:], ladder[2:]):
        assert 2.5 <= previous["front_error"] / following["front_error"] <= 5.5
        assert 2.5 <= previous["fwhm_error"] / following["fwhm_error"] <= 5.5

    # --- the derivation of the constant, recomputed from the table ---
    failing = [
        entry["points_per_front"]
        for entry in table.values()
        if entry["front_error"] > tolerance
    ]
    passing = [
        entry["points_per_front"]
        for entry in table.values()
        if entry["front_error"] <= tolerance
    ]
    assert failing and passing
    derived = int(math.floor(max(failing))) + 1
    assert derived == PREREGISTERED_MIN_POINTS_PER_FRONT

    # Restated as the property the constant is used for: everything at or
    # above the threshold meets the tolerance, and the threshold is not
    # loose -- something just below it does not.
    for entry in table.values():
        if entry["points_per_front"] >= PREREGISTERED_MIN_POINTS_PER_FRONT:
            assert entry["front_error"] <= tolerance
    assert any(
        entry["points_per_front"] < PREREGISTERED_MIN_POINTS_PER_FRONT
        and entry["front_error"] > tolerance
        for entry in table.values()
    )
    assert max(failing) == pytest.approx(6.185, rel=1.0e-3)
    assert min(passing) == pytest.approx(6.708, rel=1.0e-3)
    assert study["preregistered_min_points_per_front"] == (
        PREREGISTERED_MIN_POINTS_PER_FRONT
    )


def test_manufactured_front_study_uses_no_simulation_data() -> None:
    """The preregistration must depend on nothing but closed-form fronts.

    The module may import only NumPy and the shared grid definition: no
    experiment driver, no solver, no SciPy (the DCT is built on an even
    extension with :mod:`numpy.fft` instead).  A preregistered threshold that
    could be reached through simulation output would not be a preregistration.
    """

    source = inspect.getsource(core_width_module)
    tree = ast.parse(source)
    imports: set[tuple[int, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update((0, alias.name) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add((node.level, node.module or ""))

    assert imports <= {
        (0, "__future__"),
        (0, "math"),
        (0, "typing"),
        (0, "numpy"),
        (0, "numpy.typing"),
        (1, "grid"),
    }
    imported = {name for _, name in imports}
    for forbidden in (
        "scipy",
        "poisson",
        "pde",
        "nonlinear_cylinder",
        "cartesian",
        "artifacts",
        "pathlib",
        "json",
    ):
        assert forbidden not in imported
        assert f"ns_certificate_lab.{forbidden}" not in imported

    text = Path(core_width_module.__file__).resolve().read_text(encoding="utf-8")
    assert "np.fft.rfft" in text  # the DCT is built from NumPy's own transform
    assert "outputs/" not in text
    assert ".npz" not in text


# --------------------------------------------------------------------------
# 9.  The fit gate
# --------------------------------------------------------------------------


def test_fit_precondition_accepts_a_well_resolved_peak() -> None:
    """A generously resolved Gaussian blob must pass every criterion.

    The blob has ``sigma_r = 0.1`` on ``dr = 1/128``, so its radial FWHM spans
    ``2 sqrt(2 ln 2) sigma_r / dr = 30.1`` points and its 10-90 front spans
    ``1.687 sigma_r / dr = 21.6`` points, and it sits 64 cells from the axis.
    Measured: 30.144 points per radial FWHM, 48.228 per axial FWHM, 21.617 per
    front, no truncation, and a sub-grid peak at exactly (0.5, 0.5).
    """

    grid = _grid(129, 256)
    mesh_r, mesh_z = grid.mesh()
    field = np.exp(
        -0.5 * ((mesh_r - 0.5) / 0.1) ** 2 - 0.5 * ((mesh_z - 0.5) / 0.08) ** 2
    )

    report = peak_resolution_report(grid, field)
    assert report["peak"]["index_r"] == 64
    assert report["any_truncation"] is False
    assert report["points_per_fwhm_r"] == pytest.approx(
        GAUSSIAN_FWHM_FACTOR * 0.1 / grid.dr, rel=1.0e-3
    )
    assert report["points_per_fwhm_z"] == pytest.approx(
        GAUSSIAN_FWHM_FACTOR * 0.08 / grid.dz, rel=1.0e-3
    )
    # 10-90 thickness of a Gaussian flank:
    #   sigma (sqrt(2 ln 10) - sqrt(2 ln(10/9))).
    # Measured 21.617 points against the exact 21.593, a relative error of
    # 1.1e-03 -- the linear-interpolation error of the two crossings at
    # dr/sigma = 0.078, consistent with the manufactured-front study.
    gaussian_front = 0.1 * (
        math.sqrt(2.0 * math.log(10.0)) - math.sqrt(2.0 * math.log(10.0 / 9.0))
    )
    assert report["points_per_front"] == pytest.approx(
        gaussian_front / grid.dr, rel=5.0e-3
    )
    assert report["subgrid_peak_r"]["coordinate"] == pytest.approx(0.5, abs=1e-4)
    assert report["subgrid_peak_z"]["coordinate"] == pytest.approx(0.5, abs=1e-4)

    gate = fit_precondition(report)
    assert gate["satisfied"] is True
    assert gate["reasons"] == []
    assert gate["min_points"] == PREREGISTERED_MIN_POINTS_PER_FRONT
    assert gate["any_truncation"] is False


def test_fit_precondition_refuses_the_shipped_hou_resolution() -> None:
    """The shipped 193x384 Hou snapshot must fail the gate on three counts.

    The numbers below were measured with this module from
    ``checkpoint_nr193_nz384_t004.npz`` (the last shipped early-time
    snapshot): the ``u1`` peak sits at radial index 6, its radial FWHM spans
    6.924 points, and the 10-90 front on its axis-facing flank spans 4.364
    points.  They are transcribed rather than re-read so that this test stays
    self-contained, and they are what the audit item P0-D records: the current
    data does *not* meet the precondition for a convergence fit, and no
    amplification history changes that.
    """

    hou_like = {
        "points_per_front": 4.364122291661473,
        "points_per_fwhm_r": 6.923776126352624,
        "peak_to_axis_cells": 6,
        "truncation_flags": {
            "fwhm_r_left": False,
            "fwhm_r_right": False,
            "fwhm_z_left": False,
            "fwhm_z_right": False,
            "front_90": False,
            "front_10": False,
        },
    }

    gate = fit_precondition(hou_like)
    assert gate["satisfied"] is False
    assert len(gate["reasons"]) == 3
    joined = " | ".join(gate["reasons"])
    assert "front" in joined and "FWHM" in joined and "axis" in joined
    assert gate["any_truncation"] is False

    # Every criterion is independently load-bearing.
    only_front = dict(hou_like, points_per_fwhm_r=30.0, peak_to_axis_cells=64)
    assert fit_precondition(only_front)["satisfied"] is False
    assert len(fit_precondition(only_front)["reasons"]) == 1

    only_axis = dict(hou_like, points_per_front=20.0, points_per_fwhm_r=30.0)
    assert fit_precondition(only_axis)["satisfied"] is False

    # A truncated measurement fails the gate even when every count is large:
    # a lower bound is not a width.
    truncated = dict(
        hou_like,
        points_per_front=20.0,
        points_per_fwhm_r=30.0,
        peak_to_axis_cells=64,
        truncation_flags=dict(hou_like["truncation_flags"], fwhm_r_left=True),
    )
    result = fit_precondition(truncated)
    assert result["satisfied"] is False
    assert result["any_truncation"] is True
    assert any("truncated" in reason for reason in result["reasons"])

    # The gate can be loosened, but only by naming a different threshold, and
    # the threshold used is returned alongside the verdict so that a loosened
    # run cannot be reported as the preregistered one.  At min_points = 5 the
    # front (4.364 points) still fails while the other two criteria pass; at
    # min_points = 4 everything passes -- which is precisely why the threshold
    # is preregistered at 7 and not chosen after looking at this data.
    loosened = fit_precondition(hou_like, min_points=5)
    assert loosened["satisfied"] is False
    assert loosened["min_points"] == 5
    assert len(loosened["reasons"]) == 1
    assert "front" in loosened["reasons"][0]

    post_hoc = fit_precondition(hou_like, min_points=4)
    assert post_hoc["satisfied"] is True
    assert post_hoc["min_points"] == 4
    assert post_hoc["min_points"] < PREREGISTERED_MIN_POINTS_PER_FRONT


def test_peak_resolution_report_unwraps_the_periodic_axial_profile() -> None:
    """A peak at axial index zero is interior on the circle, not truncated.

    The shipped Hou snapshots put their peak within a handful of cells of
    ``z = 0``.  Measuring the axial width on the stored array without
    unwrapping would report a truncated -- and therefore useless -- width for
    a structure that is perfectly resolved in ``z``.
    """

    grid = _grid(129, 256)
    mesh_r, mesh_z = grid.mesh()
    # Periodic distance to z = 0, so the peak straddles the array boundary.
    wrapped = np.minimum(mesh_z, 1.0 - mesh_z)
    field = np.exp(-0.5 * ((mesh_r - 0.5) / 0.1) ** 2 - 0.5 * (wrapped / 0.05) ** 2)

    report = peak_resolution_report(grid, field)
    assert report["peak"]["index_z"] == 0
    assert report["truncation_flags"]["fwhm_z_left"] is False
    assert report["truncation_flags"]["fwhm_z_right"] is False
    assert report["points_per_fwhm_z"] == pytest.approx(
        GAUSSIAN_FWHM_FACTOR * 0.05 / grid.dz, rel=1.0e-3
    )
    assert fit_precondition(report)["satisfied"] is True


# --------------------------------------------------------------------------
# 10.  Validation
# --------------------------------------------------------------------------


def test_invalid_inputs_raise() -> None:
    """Input guards, written independently of the module they check."""

    x = np.linspace(0.0, 1.0, 11)
    y = np.exp(-0.5 * ((x - 0.5) / 0.1) ** 2)

    with pytest.raises(ValueError, match="equal length"):
        fwhm_1d(x, y[:-1])
    with pytest.raises(ValueError, match="one-dimensional"):
        fwhm_1d(x, np.tile(y, (2, 1)))
    with pytest.raises(ValueError, match="strictly increasing"):
        fwhm_1d(x[::-1], y)
    with pytest.raises(ValueError, match="at least three"):
        fwhm_1d(x[:2], y[:2])
    with pytest.raises(ValueError, match="finite"):
        fwhm_1d(x, np.where(np.arange(x.size) == 3, np.nan, y))
    with pytest.raises(ValueError, match="positive"):
        fwhm_1d(x, np.zeros_like(y))
    with pytest.raises(ValueError, match="real-valued"):
        fwhm_1d(x, y.astype(np.complex128))
    with pytest.raises(ValueError, match="outside the profile"):
        fwhm_1d(x, y, peak_index=99)
    with pytest.raises(TypeError):
        fwhm_1d(x, y, peak_index=1.5)

    with pytest.raises(ValueError, match="side must be one of"):
        front_thickness_10_90(x, y, peak_index=5, side="left")
    with pytest.raises(TypeError):
        front_thickness_10_90(x, y, peak_index=5, side=None)
    with pytest.raises(ValueError, match="outside the profile"):
        front_thickness_10_90(x, y, peak_index=-1, side="increasing_coordinate")
    with pytest.raises(ValueError, match="outside the profile"):
        subgrid_peak(x, y, 11)

    grid = _grid(17, 32)
    mesh_r, mesh_z = grid.mesh()
    field = np.sin(2.0 * np.pi * mesh_z) * (1.0 - mesh_r**2)

    with pytest.raises(TypeError):
        gradient_length_scale("not a grid", field)
    with pytest.raises(ValueError, match="shape"):
        gradient_length_scale(grid, field[:-1])
    with pytest.raises(ValueError, match="finite"):
        gradient_length_scale(grid, np.where(field == field[0, 0], np.nan, field))
    with pytest.raises(ValueError, match="identically zero"):
        gradient_length_scale(grid, np.zeros(grid.shape))
    with pytest.raises(ValueError, match=r"\(0, 1\]"):
        gradient_length_scale(grid, field, threshold_fraction=0.0)
    with pytest.raises(ValueError, match=r"\(0, 1\]"):
        gradient_length_scale(grid, field, threshold_fraction=1.5)
    with pytest.raises(TypeError):
        gradient_length_scale(grid, field, threshold_fraction=True)
    with pytest.raises(ValueError, match="real-valued"):
        gradient_length_scale(grid, field.astype(np.complex128))

    with pytest.raises(ValueError, match="identically zero"):
        high_frequency_tail(grid, np.zeros(grid.shape))
    with pytest.raises(ValueError, match="identically zero"):
        peak_resolution_report(grid, np.zeros(grid.shape))
    with pytest.raises(TypeError):
        peak_resolution_report("not a grid", field)
    # The same call with valid arguments succeeds, so the guards above reject
    # the argument they name rather than the call shape.
    assert peak_resolution_report(grid, field)["nr"] == grid.nr

    with pytest.raises(TypeError):
        fit_precondition("not a report")
    with pytest.raises(ValueError, match="missing required keys"):
        fit_precondition({"points_per_front": 10.0})
    with pytest.raises(ValueError, match="must be a dict"):
        fit_precondition(
            {
                "points_per_front": 10.0,
                "points_per_fwhm_r": 10.0,
                "peak_to_axis_cells": 10,
                "truncation_flags": "none",
            }
        )
    with pytest.raises(ValueError, match="non-finite"):
        fit_precondition(
            {
                "points_per_front": float("nan"),
                "points_per_fwhm_r": 10.0,
                "peak_to_axis_cells": 10,
                "truncation_flags": {},
            }
        )
    with pytest.raises(ValueError, match="at least one"):
        fit_precondition(
            {
                "points_per_front": 10.0,
                "points_per_fwhm_r": 10.0,
                "peak_to_axis_cells": 10,
                "truncation_flags": {},
            },
            min_points=0,
        )
    with pytest.raises(TypeError):
        fit_precondition(
            {
                "points_per_front": 10.0,
                "points_per_fwhm_r": 10.0,
                "peak_to_axis_cells": 10,
                "truncation_flags": {},
            },
            min_points=7.5,
        )

    with pytest.raises(ValueError, match="must not be empty"):
        manufactured_front_study(delta_over_dr_values=())
    with pytest.raises(ValueError, match="positive and finite"):
        manufactured_front_study(delta_over_dr_values=(-1.0,))
    with pytest.raises(ValueError, match="at least five"):
        manufactured_front_study(delta_over_dr_values=(4.0,), n_z=4)
    with pytest.raises(TypeError):
        manufactured_front_study(delta_over_dr_values=(4.0,), n_z=8.0)
    with pytest.raises(ValueError, match="at least four"):
        manufactured_front_study(delta_over_dr_values=(0.1,))
