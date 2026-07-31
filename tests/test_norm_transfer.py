r"""Tests for the norm-transfer ladder and the extended axis conditions.

The reference profile pair is

    u1   = exp(-(r^2+z^2)),
    psi1 = (1 + r^2/2) exp(-(r^2+z^2)),

for which every quantity in Theorem NT has a closed form.  Writing
``G = exp(-(r^2+z^2))``, ``A_k = int_0^inf r^k e^{-2r^2} dr`` and
``B_m = int_R z^m e^{-2z^2} dz``, the values used below follow from

    A_3 = A_5 = 1/8,  A_7 = 3/16,  A_9 = 3/8,  A_11 = 15/16,
    B_0 = sqrt(pi/2),  B_2 = B_0/4,  B_4 = 3 B_0/16,

together with the exact derivatives

    d_r psi1  = -r (1 + r^2) G,
    d_z psi1  = -z (2 + r^2) G,
    L5 psi1   = (2 r^4 - 5 r^2 - 6 + 4 z^2 + 2 r^2 z^2) G.

Substituting and collecting powers of ``r^2`` gives the four rungs
``E_0(u1)^2 = B_0/8``, ``E_1(u1)^2 = 5 B_0/8``, ``E_1(psi1)^2 = 75 B_0/64``
and ``E_2(psi1)^2 = 409 B_0/64``.  These are exact numbers, not fits, so a
discretisation that drifts away from them under refinement is wrong rather
than merely coarse.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from ns_certificate_lab.axis import (
    check_axis_regularity,
    check_ladder_axis_regularity,
    required_axis_derivative_order,
)
from ns_certificate_lab.cartesian_validation import (
    RegularizedAxisymmetricProfile,
    UniformCartesianGrid,
    reconstruct_axisymmetric_regularized,
)
from ns_certificate_lab.critical_l3 import axisymmetric_volume_weights
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.norm_transfer import (
    apply_ladder_operator,
    cartesian_sobolev_norm,
    five_dimensional_algebra_warning,
    homogeneous_sobolev_norm_from_generators,
    isometry_defect,
    ladder_norm,
    mu_inner,
    mu_norm,
    mu_weights,
    rho_inner,
    rho_norm,
    rho_weights,
    sobolev_norm_from_generators,
)
from ns_certificate_lab.operators import derivative_r, derivative_z, laplacian_5d_formal


B0 = math.sqrt(math.pi / 2.0)
EXACT_E0_U1_SQUARED = B0 / 8.0
EXACT_E1_U1_SQUARED = 5.0 * B0 / 8.0
EXACT_E1_PSI1_SQUARED = 75.0 * B0 / 64.0
EXACT_E2_PSI1_SQUARED = 409.0 * B0 / 64.0
EXACT_HDOT0 = 2.0 * math.pi * (EXACT_E0_U1_SQUARED + EXACT_E1_PSI1_SQUARED)
EXACT_HDOT1 = 2.0 * math.pi * (EXACT_E1_U1_SQUARED + EXACT_E2_PSI1_SQUARED)


def _grid(nr: int = 161) -> AxisymmetricGrid:
    """Square-celled half-plane grid with the Gaussian tail fully resolved.

    ``r_max = 5`` puts the profile at ``exp(-50)`` on the wall, so truncation
    of the whole-space integral is far below every tolerance used here and any
    disagreement that appears is a disagreement about the theorem.
    """

    return AxisymmetricGrid.uniform(
        nr=nr,
        nz=2 * (nr - 1) + 1,
        r_max=5.0,
        z_min=-5.0,
        z_max=5.0,
        periodic_z=False,
    )


def _generators(grid: AxisymmetricGrid) -> tuple[np.ndarray, np.ndarray]:
    r, z = grid.mesh()
    gaussian = np.exp(-(r * r + z * z))
    return gaussian, (1.0 + 0.5 * r * r) * gaussian


def _profile() -> RegularizedAxisymmetricProfile:
    """The same pair as point evaluators, with hand-derived derivatives."""

    def gaussian(r: np.ndarray, z: np.ndarray) -> np.ndarray:
        return np.exp(-(r * r + z * z))

    return RegularizedAxisymmetricProfile(
        u1=gaussian,
        psi1=lambda r, z: (1.0 + 0.5 * r * r) * gaussian(r, z),
        dpsi1_dr=lambda r, z: -r * (1.0 + r * r) * gaussian(r, z),
        dpsi1_dz=lambda r, z: -z * (2.0 + r * r) * gaussian(r, z),
    )


def _trapezoid_weights(axis: np.ndarray) -> np.ndarray:
    column = np.full(axis.size, float(axis[1] - axis[0]), dtype=np.float64)
    column[0] *= 0.5
    column[-1] *= 0.5
    return column


# --------------------------------------------------------------------------- #
# the ladder against closed forms                                              #
# --------------------------------------------------------------------------- #


def test_ladder_rungs_match_the_closed_form_gaussian_values() -> None:
    """Each rung is checked against an exact integral, not against itself.

    A failure here is a failure of the discrete ladder to represent
    ``E_k``: a wrong weight, a wrong axis closure in ``L5``, or a wrong
    definition of the odd rungs as gradients.  Because the four values are hit
    simultaneously, an error common to all of them, such as a missing factor in
    the quadrature, cannot hide.
    """

    grid = _grid()
    u1, psi1 = _generators(grid)
    assert ladder_norm(grid, u1, 0) ** 2 == pytest.approx(
        EXACT_E0_U1_SQUARED, rel=1e-6
    )
    assert ladder_norm(grid, u1, 1) ** 2 == pytest.approx(
        EXACT_E1_U1_SQUARED, rel=2e-3
    )
    assert ladder_norm(grid, psi1, 1) ** 2 == pytest.approx(
        EXACT_E1_PSI1_SQUARED, rel=2e-3
    )
    assert ladder_norm(grid, psi1, 2) ** 2 == pytest.approx(
        EXACT_E2_PSI1_SQUARED, rel=2e-3
    )


def test_even_rung_zero_is_the_plain_mu_norm() -> None:
    """Rung zero must not secretly apply the operator once.

    An off-by-one in the ``k // 2`` bookkeeping would make every even rung one
    application of ``L5`` too deep, which for smooth data is a plausible-looking
    number rather than an obvious blow-up.
    """

    grid = _grid(41)
    u1, _ = _generators(grid)
    assert ladder_norm(grid, u1, 0) == pytest.approx(mu_norm(grid, u1), rel=0.0)


def test_apply_ladder_operator_is_repeated_application_of_the_whole_operator() -> None:
    """The contract is that ``L5^m`` is never expanded term by term.

    The individual terms of an expanded ``L5^2`` carry negative powers of ``r``
    that are separately unbounded at the axis; only the assembled operator is
    finite there.  Bit-for-bit agreement with a hand-written double application
    is the strongest available statement that no expansion has crept in.
    """

    grid = _grid(41)
    _, psi1 = _generators(grid)
    twice = laplacian_5d_formal(grid, laplacian_5d_formal(grid, psi1))
    assert np.array_equal(apply_ladder_operator(grid, psi1, 2), twice)
    assert np.array_equal(apply_ladder_operator(grid, psi1, 0), psi1)


# --------------------------------------------------------------------------- #
# the isometry itself                                                          #
# --------------------------------------------------------------------------- #


def test_l2_isometry_matches_a_fully_three_dimensional_quadrature() -> None:
    """The strongest independent check available: no cylindrical code at all.

    The field is rebuilt as three Cartesian components on a uniform 3D grid by
    :func:`reconstruct_axisymmetric_regularized`, which never divides by ``r``
    and never sees the half-plane grid, and ``int |u|^2 dV`` is taken with a
    plain tensor trapezoid rule.  For a Gaussian on a box this quadrature is
    accurate to roughly machine precision, so the comparison pins the constant
    ``2 pi`` and the ``r^3`` weight exactly.  A failure means the generator to
    velocity map and the ladder disagree about what field is being measured.
    """

    grid = _grid()
    u1, psi1 = _generators(grid)
    ladder = homogeneous_sobolev_norm_from_generators(grid, u1, psi1, 0)

    cartesian = UniformCartesianGrid.uniform(
        shape=(65, 65, 65),
        x_bounds=(-4.5, 4.5),
        y_bounds=(-4.5, 4.5),
        z_bounds=(-4.5, 4.5),
    )
    velocity = reconstruct_axisymmetric_regularized(cartesian, _profile())
    weights = (
        _trapezoid_weights(cartesian.x)[:, None, None]
        * _trapezoid_weights(cartesian.y)[None, :, None]
        * _trapezoid_weights(cartesian.z)[None, None, :]
    )
    volume_integral = float(np.sum(np.sum(velocity * velocity, axis=0) * weights))

    assert volume_integral == pytest.approx(EXACT_HDOT0, rel=1e-9)
    assert ladder.total_squared == pytest.approx(volume_integral, rel=1e-3)


def test_homogeneous_isometry_holds_at_the_two_cross_checked_rungs() -> None:
    """Ladder and physical routes agree, and both agree with the closed form.

    The two routes share no arithmetic: one applies ``L5`` and integrates
    against ``r^3``, the other differentiates ``u^r``, ``u^theta``, ``u^z`` and
    integrates against ``2 pi r``.  Agreement to the discretisation error, and
    with the exact value, is the numerical content of Theorem NT at ``n = 0``
    and ``n = 1``.
    """

    grid = _grid()
    u1, psi1 = _generators(grid)
    for n, exact in ((0, EXACT_HDOT0), (1, EXACT_HDOT1)):
        ladder = homogeneous_sobolev_norm_from_generators(grid, u1, psi1, n)
        physical = cartesian_sobolev_norm(grid, u1, psi1, n)
        assert ladder.total_squared == pytest.approx(exact, rel=2e-3)
        assert physical.total_squared == pytest.approx(exact, rel=2e-3)
        assert ladder.total_squared == pytest.approx(
            physical.total_squared, rel=2e-3
        )


def test_the_swirl_and_stream_halves_agree_separately() -> None:
    """The split must match term by term, not only in the sum.

    Two errors of opposite sign in the two halves would cancel in the total and
    leave the isometry looking healthy.  Comparing ``E_n(u1)`` against the
    swirl part of the physical integrand and ``E_{n+1}(psi1)`` against the
    meridional part removes that escape route.
    """

    grid = _grid()
    u1, psi1 = _generators(grid)
    ladder = homogeneous_sobolev_norm_from_generators(grid, u1, psi1, 1)
    physical = cartesian_sobolev_norm(grid, u1, psi1, 1)
    two_pi = 2.0 * math.pi
    assert two_pi * ladder.swirl_squared == pytest.approx(
        physical.swirl_part, rel=2e-3
    )
    assert two_pi * ladder.stream_squared == pytest.approx(
        physical.meridional_part, rel=2e-3
    )
    assert ladder.swirl_fraction == pytest.approx(
        ladder.swirl_squared / ladder.bracket, rel=1e-12
    )


def test_isometry_error_falls_at_second_order() -> None:
    """A constant residual would mean a wrong identity, not a coarse grid.

    Both routes are second-order accurate, so their difference must fall like
    ``h^2``.  If the identity itself were wrong the difference would tend to a
    nonzero limit and the observed order would collapse towards zero, which is
    exactly what this test refuses to accept.
    """

    spacings: list[float] = []
    defects: dict[int, list[float]] = {0: [], 1: []}
    for nr in (41, 81, 161):
        grid = _grid(nr)
        u1, psi1 = _generators(grid)
        spacings.append(grid.dr)
        for n in (0, 1):
            defects[n].append(isometry_defect(grid, u1, psi1, n))

    for n in (0, 1):
        values = defects[n]
        assert values[0] > values[1] > values[2] > 0.0
        for coarse, fine, h_coarse, h_fine in zip(
            values[:-1], values[1:], spacings[:-1], spacings[1:]
        ):
            order = math.log(coarse / fine) / math.log(h_coarse / h_fine)
            assert order > 1.8, (n, order)
        assert values[-1] < 1e-3


# --------------------------------------------------------------------------- #
# the two weights                                                              #
# --------------------------------------------------------------------------- #


def test_mu_and_rho_differ_by_exactly_one_half_on_a_known_profile() -> None:
    """The two weights are distinguished by a case with an exact ratio.

    For ``f = exp(-(r^2+z^2))`` the axial factor cancels and the ratio of the
    two inner products is ``A_3 / A_1 = (1/8)/(1/4) = 1/2`` exactly.  A test
    that only checked ``mu_inner`` against itself could not tell the two
    weights apart; this one fails loudly if either helper silently uses the
    other measure, and the specific value one half is what a transposed pair of
    calls would turn into two.
    """

    grid = _grid()
    r, z = grid.mesh()
    field = np.exp(-(r * r + z * z))
    ratio = mu_inner(grid, field, field) / rho_inner(grid, field, field)
    assert ratio == pytest.approx(0.5, rel=1e-3)
    assert mu_norm(grid, field) ** 2 == pytest.approx(
        mu_inner(grid, field, field), rel=1e-12
    )
    assert rho_norm(grid, field) ** 2 == pytest.approx(
        rho_inner(grid, field, field), rel=1e-12
    )


def test_physical_weights_agree_bit_for_bit_with_the_critical_l3_measure() -> None:
    """One physical measure in the repository, not two that nearly agree.

    ``critical_l3`` already owns a ``2 pi r dr dz`` rule.  If this module used a
    different endpoint convention the two would differ by half a cell at the
    wall, which is invisible on smooth data and fatal on data that has not
    decayed.  Exact equality is required because both are built from the same
    conventions on purpose.
    """

    grid = _grid(41)
    assert np.array_equal(
        2.0 * math.pi * rho_weights(grid), axisymmetric_volume_weights(grid)
    )
    assert np.array_equal(
        mu_weights(grid), rho_weights(grid) * (grid.r**2)[:, None]
    )


def test_integration_by_parts_holds_in_mu_and_fails_in_rho() -> None:
    """``<psi1,omega1>_mu = E_1(psi1)^2`` is the reason the weight is ``r^3``.

    ``L5`` is self-adjoint for ``r^3 dr dz`` and for no other power, so the
    identity that converts the pairing into a gradient norm holds in the ladder
    measure and is simply false in the physical one.  The second assertion is
    the important one: it records that the failure in the ``r`` weight is a
    gross ``O(1)`` failure, not a small perturbation that might be tolerated.
    """

    grid = _grid()
    _, psi1 = _generators(grid)
    omega1 = -laplacian_5d_formal(grid, psi1)

    mu_pairing = mu_inner(grid, psi1, omega1)
    assert mu_pairing == pytest.approx(ladder_norm(grid, psi1, 1) ** 2, rel=1e-3)
    assert mu_pairing == pytest.approx(EXACT_E1_PSI1_SQUARED, rel=2e-3)

    rho_pairing = rho_inner(grid, psi1, omega1)
    radial = derivative_r(grid, psi1, even_at_axis=True)
    axial = derivative_z(grid, psi1)
    rho_gradient = float(
        np.sum((radial * radial + axial * axial) * rho_weights(grid))
    )
    assert abs(rho_pairing - rho_gradient) / rho_gradient > 0.5


# --------------------------------------------------------------------------- #
# the inhomogeneous identity                                                   #
# --------------------------------------------------------------------------- #


def test_inhomogeneous_transfer_is_the_l2_norm_plus_the_homogeneous_part() -> None:
    """``H^n`` is assembled from two pieces that are separately checked.

    The reported bracket must reproduce ``||u||^2_{L^2} + ||u||^2_{Hdot^n}``
    with the ``L^2`` part written through ``<psi1,omega1>_mu``.  A failure means
    the rewrite that replaced ``E_1(psi1)^2`` by the pairing has gone wrong, and
    since the two agree only after an integration by parts, that is precisely
    the step most likely to be broken by a bad axis closure.
    """

    grid = _grid()
    u1, psi1 = _generators(grid)
    full = sobolev_norm_from_generators(grid, u1, psi1, 1)
    low = homogeneous_sobolev_norm_from_generators(grid, u1, psi1, 0)
    top = homogeneous_sobolev_norm_from_generators(grid, u1, psi1, 1)
    assert full.total_squared == pytest.approx(
        low.total_squared + top.total_squared, rel=1e-3
    )
    assert full.total_squared == pytest.approx(EXACT_HDOT0 + EXACT_HDOT1, rel=2e-3)
    assert full.norm == pytest.approx(math.sqrt(full.total_squared), rel=1e-12)
    assert full.dominant_term == "stream_top_squared"


def test_the_omega1_rewrite_of_the_top_term_is_exact_in_floating_point() -> None:
    """``E_{n-1}(omega1) = E_{n+1}(psi1)`` is algebra, not approximation.

    ``omega1 = -L5 psi1`` and ``L5`` is linear with sign-symmetric rounding, so
    the two evaluations must agree to the last bit rather than to the
    discretisation error.  A discrepancy would mean the two routes disagree
    about how many times ``L5`` has been applied.
    """

    grid = _grid(41)
    u1, psi1 = _generators(grid)
    for n in (1, 2, 3):
        full = sobolev_norm_from_generators(grid, u1, psi1, n)
        assert full.stream_top_squared == pytest.approx(
            ladder_norm(grid, psi1, n + 1) ** 2, rel=1e-13
        )
        assert full.swirl_top_squared == pytest.approx(
            ladder_norm(grid, u1, n) ** 2, rel=1e-13
        )


def test_reports_are_json_safe() -> None:
    """Every report has to survive canonical serialisation.

    The repository stores diagnostics as JSON with ``allow_nan=False``, so a
    dataclass that leaks a NumPy scalar or a non-finite value would fail at
    archive time rather than at computation time, which is the worst moment to
    discover it.
    """

    grid = _grid(41)
    u1, psi1 = _generators(grid)
    payloads = [
        sobolev_norm_from_generators(grid, u1, psi1, 2).as_dict(),
        homogeneous_sobolev_norm_from_generators(grid, u1, psi1, 2).as_dict(),
        cartesian_sobolev_norm(grid, u1, psi1, 1).as_dict(),
    ]
    for payload in payloads:
        text = json.dumps(payload, allow_nan=False, sort_keys=True)
        restored = json.loads(text)
        assert restored["n"] in (1, 2)
        for key, value in restored.items():
            assert isinstance(value, (int, float, str))
            if isinstance(value, float):
                assert math.isfinite(value)


# --------------------------------------------------------------------------- #
# rejection of bad input                                                       #
# --------------------------------------------------------------------------- #


def test_sobolev_transfer_rejects_order_zero_rather_than_guessing() -> None:
    """``H^0`` has no ``E_{n-1}(omega1)`` term and must not be improvised.

    Collapsing ``n = 0`` onto the ``L^2`` formula would double count the swirl
    term and produce a number too large by a factor that depends on the field,
    so the caller is made to say which norm is meant.
    """

    grid = _grid()
    u1, psi1 = _generators(grid)
    with pytest.raises(ValueError):
        sobolev_norm_from_generators(grid, u1, psi1, 0)
    # The honest n = 0 statement is available and is the L^2 norm.
    assert homogeneous_sobolev_norm_from_generators(
        grid, u1, psi1, 0
    ).total_squared == pytest.approx(EXACT_HDOT0, rel=2e-3)


def test_cross_check_refuses_rungs_it_cannot_independently_verify() -> None:
    """Silence about ``n >= 2`` is better than a second copy of the ladder.

    The physical route is genuinely independent only where the field itself is
    differentiated.  Returning a number for ``n = 2`` by applying ``L5`` again
    would look like a cross-check while testing nothing.
    """

    grid = _grid(41)
    u1, psi1 = _generators(grid)
    with pytest.raises(ValueError):
        cartesian_sobolev_norm(grid, u1, psi1, 2)


@pytest.mark.parametrize("bad", [-1, True, 1.5, "2"])
def test_ladder_rejects_a_rung_that_is_not_a_nonnegative_integer(bad: object) -> None:
    """``k`` indexes derivatives; a boolean or a float is a caller bug."""

    grid = _grid(41)
    u1, _ = _generators(grid)
    with pytest.raises(ValueError):
        ladder_norm(grid, u1, bad)  # type: ignore[arg-type]


def test_transfer_rejects_fields_of_the_wrong_shape_or_with_holes() -> None:
    """Shape and finiteness are checked before any quadrature runs.

    A transposed array would otherwise be broadcast into a plausible number,
    and a single NaN would propagate into a total that reads as ``nan`` only
    after the caller has stored it.
    """

    grid = _grid(41)
    u1, psi1 = _generators(grid)
    with pytest.raises(ValueError):
        mu_inner(grid, u1.T, psi1.T)
    with pytest.raises(ValueError):
        rho_inner(grid, u1[:-1], psi1[:-1])
    broken = psi1.copy()
    broken[3, 4] = np.nan
    with pytest.raises(ValueError):
        homogeneous_sobolev_norm_from_generators(grid, u1, broken, 1)


# --------------------------------------------------------------------------- #
# the axis conditions that the old check omitted                               #
# --------------------------------------------------------------------------- #


def _axis_test_grid() -> AxisymmetricGrid:
    """Coarse enough that a third derivative is not dominated by rounding.

    A one-sided third derivative amplifies floating-point noise by ``h^-3``, so
    the demonstration grid is kept coarse deliberately; the defect being
    detected is ``O(1)`` and does not need resolution.
    """

    return AxisymmetricGrid.uniform(
        nr=41, nz=33, r_max=5.0, z_min=-2.0, z_max=2.0, periodic_z=False
    )


def _polynomial_pair(grid: AxisymmetricGrid) -> tuple[np.ndarray, np.ndarray]:
    """An exactly axis-even profile: a polynomial in ``r^2`` times a Gaussian.

    Polynomials in ``r^2`` are annihilated exactly by the one-sided stencils
    used in the audit, so the clean case has no truncation error to confuse
    with a genuine parity defect.
    """

    r, z = grid.mesh()
    envelope = np.exp(-z * z)
    return envelope * (1.0 + 0.25 * r * r), envelope * (
        1.0 + 0.5 * r * r + 0.125 * r**4
    )


def test_extended_axis_check_catches_a_cubic_defect_the_old_check_missed() -> None:
    """This is the bug: ``d_r^3 psi1(0,z) != 0`` used to pass unnoticed.

    Adding ``c r^3`` to an axis-even field leaves ``d_r psi1(0,z)`` exactly
    zero, so the historical order-one check reports a clean field.  The field is
    nevertheless not the restriction of any smooth axisymmetric function and it
    breaks the ``H^4`` transfer, because ``L5`` turns the ``r^3`` term into an
    odd ``15 c r`` term which the next application of ``L5`` differentiates as
    if it were even.  A failure of the first assertion would mean the omission
    has been reintroduced; a failure of the second would mean the new check is
    firing on the wrong quantity.
    """

    grid = _axis_test_grid()
    _, clean = _polynomial_pair(grid)
    r, z = grid.mesh()
    defective = clean + 0.25 * r**3 * np.exp(-z * z)

    old_check = check_axis_regularity(
        grid, even_fields={"psi1": defective}, tolerance=1e-8
    )
    assert old_check.passed, old_check.to_dict()
    assert old_check.metrics["psi1.max_abs_dr_at_axis"] < 1e-8

    new_check = check_axis_regularity(
        grid, even_fields={"psi1": defective}, tolerance=1e-8, derivative_order=3
    )
    assert not new_check.passed
    assert new_check.violations == ("psi1.max_abs_dr3_at_axis",)
    # d_r^3 (c r^3) = 6c, so the residual is the defect itself, not a proxy.
    assert new_check.metrics["psi1.max_abs_dr3_at_axis"] == pytest.approx(
        6.0 * 0.25, rel=1e-6
    )


def test_extended_axis_check_still_passes_a_genuinely_even_field() -> None:
    """The new condition must not condemn correct data.

    A high-order one-sided stencil is noisy, and a check that flagged every
    smooth field would be worse than no check at all because it would be
    switched off.
    """

    grid = _axis_test_grid()
    u1, psi1 = _polynomial_pair(grid)
    report = check_axis_regularity(
        grid,
        even_fields={"u1": u1, "psi1": psi1},
        tolerance=1e-8,
        derivative_order=3,
    )
    assert report.passed, report.to_dict()
    assert set(report.metrics) == {
        "u1.max_abs_dr_at_axis",
        "u1.max_abs_dr3_at_axis",
        "psi1.max_abs_dr_at_axis",
        "psi1.max_abs_dr3_at_axis",
    }


def test_default_axis_report_is_unchanged_by_the_extension() -> None:
    """The historical API and its metric keys are frozen.

    Recorded artifacts contain ``max_abs_dr_at_axis`` values produced by the
    original stencil.  If the default call changed key names or numbers, every
    stored audit would become unreproducible, which is a worse failure than the
    omission being fixed.
    """

    grid = _axis_test_grid()
    u1, psi1 = _polynomial_pair(grid)
    r, _ = grid.mesh()
    report = check_axis_regularity(
        grid,
        even_fields={"psi1": psi1},
        odd_fields={"u_theta": r * u1},
        tolerance=1e-8,
    )
    assert list(report.metrics) == [
        "psi1.max_abs_dr_at_axis",
        "u_theta.max_abs_value_at_axis",
    ]
    assert report.passed, report.to_dict()


def test_odd_fields_gain_the_matching_even_order_conditions() -> None:
    """Parity forbids the even derivatives of an odd field, not only its value.

    ``u^theta = r u1`` vanishing at the axis is the order-zero condition; a
    field with a nonzero second radial derivative there is equally impossible
    and used to pass.
    """

    grid = _axis_test_grid()
    r, z = grid.mesh()
    honest = r * (1.0 + 0.25 * r * r) * np.exp(-z * z)
    defective = honest + 0.5 * r * r * np.exp(-z * z)
    good = check_axis_regularity(
        grid, odd_fields={"u_theta": honest}, tolerance=1e-8, derivative_order=2
    )
    assert good.passed, good.to_dict()
    bad = check_axis_regularity(
        grid, odd_fields={"u_theta": defective}, tolerance=1e-8, derivative_order=2
    )
    assert bad.violations == ("u_theta.max_abs_dr2_at_axis",)
    assert bad.metrics["u_theta.max_abs_dr2_at_axis"] == pytest.approx(1.0, rel=1e-6)


@pytest.mark.parametrize(
    ("sobolev_order", "expected"),
    [(0, 1), (1, 1), (2, 1), (3, 3), (4, 3), (5, 5), (6, 5)],
)
def test_required_axis_order_tracks_the_depth_of_the_ladder(
    sobolev_order: int, expected: int
) -> None:
    """The required order is set by how often ``L5`` is applied, not by taste.

    ``psi1`` is carried to rung ``n+1``, which applies ``L5`` exactly
    ``floor((n+1)/2)`` times, and each application consumes one more odd
    condition.  ``n = 4`` needing order three is the case that motivated the
    fix; a change to this table changes which fields the transfer accepts.
    """

    assert required_axis_derivative_order(sobolev_order) == expected


def test_ladder_axis_helper_applies_the_required_order_to_both_generators() -> None:
    """The convenience wrapper must not audit a weaker condition than needed."""

    grid = _axis_test_grid()
    u1, psi1 = _polynomial_pair(grid)
    r, z = grid.mesh()
    defective = psi1 + 0.25 * r**3 * np.exp(-z * z)
    assert check_ladder_axis_regularity(
        grid, u1=u1, psi1=psi1, sobolev_order=4, tolerance=1e-8
    ).passed
    report = check_ladder_axis_regularity(
        grid, u1=u1, psi1=defective, sobolev_order=4, tolerance=1e-8
    )
    assert not report.passed
    assert "psi1.max_abs_dr3_at_axis" in report.violations
    # At n = 2 the ladder never applies L5 twice, so the same field is accepted.
    assert check_ladder_axis_regularity(
        grid, u1=u1, psi1=defective, sobolev_order=2, tolerance=1e-8
    ).passed


def test_axis_check_rejects_bad_orders_and_undersized_grids() -> None:
    """A silently truncated stencil would report a meaningless residual."""

    grid = _axis_test_grid()
    _, psi1 = _polynomial_pair(grid)
    with pytest.raises(ValueError):
        check_axis_regularity(grid, even_fields={"psi1": psi1}, derivative_order=-1)
    with pytest.raises(ValueError):
        check_axis_regularity(grid, even_fields={"psi1": psi1}, derivative_order=1.0)
    with pytest.raises(ValueError):
        required_axis_derivative_order(-1)
    tiny = AxisymmetricGrid.uniform(
        nr=4, nz=8, r_max=1.0, z_min=-1.0, z_max=1.0, periodic_z=False
    )
    r, z = tiny.mesh()
    with pytest.raises(ValueError):
        check_axis_regularity(
            tiny,
            even_fields={"psi1": np.exp(-(r * r + z * z))},
            derivative_order=3,
        )


# --------------------------------------------------------------------------- #
# the standing warning                                                         #
# --------------------------------------------------------------------------- #


def test_five_dimensional_algebra_warning_keeps_all_five_obstructions() -> None:
    """Losing one reason is how a norm identity becomes an existence claim.

    Each of the five is fatal on its own, so the list is asserted item by item
    rather than by its length alone: a rewrite that quietly dropped the
    vortex-stretching sentence would leave a warning that still reads
    plausibly and is no longer true.
    """

    reasons = five_dimensional_algebra_warning()
    assert len(reasons) == 5
    assert all(isinstance(reason, str) and reason for reason in reasons)
    joined = " ".join(reasons).lower()
    for fragment in (
        "div_5 = 2 u^r / r",
        "vortex-stretching",
        "no five-dimensional pressure",
        "off by one derivative",
        "multiplier",
        "axis rather than infinity",
        "s > 5/2",
        "s > 3/2",
    ):
        assert fragment.lower() in joined, fragment
