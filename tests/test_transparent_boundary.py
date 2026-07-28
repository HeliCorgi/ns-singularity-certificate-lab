"""Acceptance tests for milestone W-A: the transparent outer boundary condition.

Five things are pinned here.

1. **The Dirichlet path is untouched.**  ``solve_streamfunction_poisson_outer``
   with its default ``boundary_condition="dirichlet"`` must return results that
   are *bitwise* equal to solver A's, field for field and metadata key for
   metadata key.  The transparent option is only allowed to exist if it costs
   the existing path nothing.

2. **Only the outermost radial row changes.**  The flux coefficients and the
   assembled tridiagonal rows of the transparent path are compared bitwise with
   solver A's on every row except the last, and the appended row is compared
   against an independent evaluation of the E-26 formula written out here.

3. **The condition and its two limits.**  ``K_0/K_1`` is checked against a
   published value, against an ascending series written out in this module, and
   against the derivative identity ``K_1' = -K_0 - K_1/x`` from which (W-1) is
   derived.  The bracket is then checked to equal ``2/R`` exactly at ``k = 0``
   -- with the ``k = 0`` decaying solution ``C/r^2`` satisfying the condition
   identically -- and to approach ``2/R + k`` as ``kR -> infinity``.

4. **The six W-B acceptance conditions**, at reduced resolution, through the
   experiment module.

5. **The refusals.**  Config validation, the compact-support premise, and the
   three fault injections.

This module imports solver A (to pin the bit-identity) and never imports solver
B, so it stays compatible with the guard in
``tests/test_poisson_cross_validation.py``.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from experiments import run_transparent_boundary as experiment
from experiments.run_transparent_boundary import (
    EXPERIMENT_ID,
    INJECTED_FAULTS,
    acceptance_report,
    build_grid,
    build_source,
    decaying_branch_derivatives,
    evaluate,
    manufactured_case_fields,
    manufactured_polynomial_coefficients,
    measure_support_leak,
    radial_point_count,
    run,
    source_moment,
    validate_config,
)
from ns_certificate_lab import bessel_reference
from ns_certificate_lab import transparent_boundary as transparent_module
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.poisson import (
    _radial_flux_coefficients,
    solve_streamfunction_poisson,
)
from ns_certificate_lab.transparent_boundary import (
    BRACKET_VARIANTS,
    assemble_radial_mode,
    assert_compact_support,
    modal_wavenumbers,
    outer_bracket,
    radial_flux_coefficients,
    solve_streamfunction_poisson_outer,
    transparent_condition_defect,
)

REPOSITORY = Path(__file__).resolve().parents[1]
SHIPPED_CONFIG = REPOSITORY / "configs" / "transparent_boundary.json"

EULER_MASCHERONI = 0.5772156649015328606065120900824

#: The E-33(e) zero-mode Dirichlet offset recorded in ``docs/equation_audit.md``
#: for the shipped source family at ``R = 1 -> 2``.
AUDITED_DIRICHLET_ZERO_MODE_OFFSET = -1.6969516537e-3


# --------------------------------------------------------------------------
# Independent reimplementations, written out here so that the oracle is
# checked against something that shares no code with it.
# --------------------------------------------------------------------------


def _i0_via_series(x: float, terms: int = 200) -> float:
    """``I_0(x) = sum (x/2)^{2m}/(m!)^2`` -- a self-contained series."""

    half = 0.5 * x
    term = 1.0
    total = 1.0
    for m in range(1, terms + 1):
        term *= half * half / (m * m)
        total += term
    return total


def _k0_via_ascending_series(x: float, terms: int = 120) -> float:
    """``K_0`` from DLMF 10.31.1, independent of :mod:`bessel_reference`.

    ``K_0(x) = -(ln(x/2)+gamma) I_0(x)
               + sum_{m>=1} H_m (x^2/4)^m / (m!)^2`` with
    ``H_m = sum_{j<=m} 1/j``.  Every quantity is recomputed here from scratch,
    including ``I_0``.
    """

    quarter = 0.25 * x * x
    term = 1.0
    harmonic = 0.0
    total = 0.0
    for m in range(1, terms + 1):
        term *= quarter / (m * m)
        harmonic += 1.0 / m
        total += harmonic * term
    return -(math.log(0.5 * x) + EULER_MASCHERONI) * _i0_via_series(x) + total


def _max_abs(values: np.ndarray) -> float:
    return float(np.max(np.abs(values)))


def _order(coarse: float, fine: float, ratio: float = 2.0) -> float:
    return math.log(coarse / fine) / math.log(ratio)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


# --------------------------------------------------------------------------
# 1. The K_0/K_1 oracle
# --------------------------------------------------------------------------


def test_k0_over_k1_reproduces_published_values() -> None:
    """DLMF 10.31.1 / A&S Table 9.8 values, and the ratio they imply."""

    assert bessel_reference.bessel_k0(1.0) == pytest.approx(
        0.4210244382407083, rel=1e-13
    )
    assert bessel_reference.bessel_k1(1.0) == pytest.approx(
        0.6019072301972346, rel=1e-13
    )
    expected = 0.4210244382407083 / 0.6019072301972346
    assert bessel_reference.k0_over_k1(1.0) == pytest.approx(expected, rel=1e-13)


def test_k0_over_k1_matches_an_independent_ascending_series() -> None:
    """The quadrature ``K_0`` must agree with a different formula entirely."""

    for argument in (0.25, 0.5, 1.0, 2.0, 4.0):
        series = _k0_via_ascending_series(argument)
        assert bessel_reference.bessel_k0(argument) == pytest.approx(
            series, rel=1e-12
        )
        ratio = series / bessel_reference.bessel_k1(argument)
        assert bessel_reference.k0_over_k1(argument) == pytest.approx(
            ratio, rel=1e-12
        )


def test_k0_over_k1_satisfies_the_derivative_identity_behind_w1() -> None:
    """``K_1'(x) = -K_0(x) - K_1(x)/x`` is the identity (W-1) is read off from.

    If it failed, the whole derivation of the transparent bracket would be
    wrong, so it is checked directly with a fourth-order central difference
    rather than assumed.
    """

    for value in (0.25, 0.5, 1.0, 2.0, 4.0):
        step = 1.0e-3 * value
        derivative = (
            bessel_reference.bessel_k1(value - 2.0 * step)
            - 8.0 * bessel_reference.bessel_k1(value - step)
            + 8.0 * bessel_reference.bessel_k1(value + step)
            - bessel_reference.bessel_k1(value + 2.0 * step)
        ) / (12.0 * step)
        predicted = (
            -bessel_reference.bessel_k0(value)
            - bessel_reference.bessel_k1(value) / value
        )
        assert derivative == pytest.approx(predicted, rel=1e-9)


def test_k0_over_k1_is_monotone_inside_the_unit_interval() -> None:
    """``0 < K_0/K_1 < 1`` and increasing: what bounds the transparent bracket."""

    ladder = [1e-6, 1e-3, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 60.0, 400.0]
    values = [bessel_reference.k0_over_k1(x) for x in ladder]
    assert all(0.0 < value < 1.0 for value in values)
    assert all(values[i] < values[i + 1] for i in range(len(values) - 1))
    # The naive unscaled quotient would be 0/0 at this argument; the scaled
    # form still returns the right number.
    assert bessel_reference.bessel_k0(760.0) == 0.0
    assert bessel_reference.bessel_k1(760.0) == 0.0
    assert bessel_reference.k0_over_k1(760.0) == pytest.approx(
        1.0 - 0.5 / 760.0, rel=1e-5
    )


def test_k0_over_k1_asymptotic_branches() -> None:
    for x in (1e-6, 1e-5, 1e-4):
        assert bessel_reference.k0_over_k1_small_argument_asymptote(
            x
        ) == pytest.approx(bessel_reference.k0_over_k1(x), rel=1e-3)
    for x in (20.0, 40.0, 60.0):
        assert bessel_reference.k0_over_k1_large_argument_asymptote(
            x
        ) == pytest.approx(bessel_reference.k0_over_k1(x), rel=1e-4)
    assert bessel_reference.k0_over_k1_large_argument_asymptote(10.0, terms=1) == 1.0
    with pytest.raises(ValueError):
        bessel_reference.k0_over_k1_large_argument_asymptote(10.0, terms=4)
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError):
            bessel_reference.k0_over_k1(bad)


def test_bessel_cross_product_identity_still_holds() -> None:
    """The new helpers must not have disturbed the existing oracle."""

    for x in (1e-2, 0.5, 1.0, 6.283185307179586, 60.0):
        assert bessel_reference.wronskian_relative_defect(x) <= 1e-14


# --------------------------------------------------------------------------
# 2. The bracket and its two limits
# --------------------------------------------------------------------------


def test_zero_mode_bracket_is_exactly_two_over_r() -> None:
    """W-1 at ``k = 0``: the bracket is ``2/R``, from an explicit branch."""

    for radius in (0.5, 1.0, 1.7, 3.0):
        assert outer_bracket(0.0, radius) == 2.0 / radius


def test_zero_mode_decaying_solution_satisfies_the_condition_identically() -> None:
    """``psi = C/r^2`` gives ``psi'(R) + (2/R) psi(R) = 0`` with no remainder."""

    for radius in (0.5, 1.0, 2.0):
        for amplitude in (1.0, -0.37, 12.5):
            derivative = -2.0 * amplitude / radius**3
            residual = derivative + outer_bracket(0.0, radius) * (
                amplitude / radius**2
            )
            assert residual == 0.0


def test_bracket_small_wavenumber_limit_is_two_over_r() -> None:
    """``k K_0(kR)/K_1(kR) -> 0`` faster than ``k``, so beta -> 2/R."""

    radius = 1.3
    previous = math.inf
    for k in (1.0, 1e-1, 1e-2, 1e-3, 1e-5, 1e-7):
        excess = outer_bracket(k, radius) - 2.0 / radius
        assert excess > 0.0
        assert excess / k < previous
        previous = excess / k
    # ``excess/k = K_0(kR)/K_1(kR) ~ kR(-log(kR/2)-gamma)``, so it vanishes
    # with ``k`` up to the logarithm; the bracket's limit really is 2/R.
    assert previous < 1e-5
    # Compared against the asymptote *without* the cancelling subtraction: at
    # k = 1e-7 the excess is 2e-13 next to 2/R = 1.5, so ``excess`` itself
    # carries only three significant digits and a tight comparison there would
    # be measuring binary64 cancellation rather than the limit.
    assert bessel_reference.k0_over_k1(1e-7 * radius) == pytest.approx(
        bessel_reference.k0_over_k1_small_argument_asymptote(1e-7 * radius),
        rel=1e-10,
    )
    assert previous == pytest.approx(
        bessel_reference.k0_over_k1(1e-7 * radius), rel=1e-2
    )


def test_bracket_large_wavenumber_limit_is_two_over_r_plus_k() -> None:
    """``K_0/K_1 -> 1``, so beta -> 2/R + k, approached like ``1/(2R)``."""

    radius = 1.3
    wavenumbers = (10.0, 100.0, 1000.0, 10000.0)
    gaps = []
    for k in wavenumbers:
        asymptote = 2.0 / radius + k
        bracket = outer_bracket(k, radius)
        assert bracket < asymptote
        gaps.append(abs(bracket - asymptote))
    # The shortfall tends to the constant 1/(2R), not to zero: it is the
    # -1/(2x) term of K_0/K_1 with x = kR, multiplied by k.  It approaches
    # that constant from below and monotonically.
    assert all(gaps[i] < gaps[i + 1] for i in range(len(gaps) - 1))
    assert all(gap < 0.5 / radius for gap in gaps)
    assert gaps[-1] == pytest.approx(0.5 / radius, rel=1e-3)
    # Relative to the bracket itself the gap vanishes, so beta/(2/R + k) -> 1.
    assert gaps[-1] / (2.0 / radius + wavenumbers[-1]) < 1e-4


def test_bracket_switches_to_the_documented_branches_outside_the_quadrature() -> None:
    """The two substitutions are made where documented, and are continuous."""

    limit = bessel_reference.LARGE_ARGUMENT_QUADRATURE_LIMIT
    radius = 1.0
    # Just below the limit the quadrature is used; just above, the asymptote.
    below = outer_bracket(limit * 0.999, radius)
    above = outer_bracket(limit * 1.001, radius)
    assert below == pytest.approx(
        2.0 / radius
        + limit * 0.999 * bessel_reference.k0_over_k1(limit * 0.999 * radius),
        rel=1e-15,
    )
    assert above == pytest.approx(
        2.0 / radius
        + limit
        * 1.001
        * bessel_reference.k0_over_k1_large_argument_asymptote(limit * 1.001 * radius),
        rel=1e-15,
    )
    # The switch is not a discontinuity: the two branches agree there to 1e-9.
    assert bessel_reference.k0_over_k1(limit) == pytest.approx(
        bessel_reference.k0_over_k1_large_argument_asymptote(limit), rel=1e-8
    )
    # And the quadrature really has degraded further out, which is the reason
    # the switch exists.
    assert bessel_reference.k_quadrature_step_halving_defect(3000.0, 0) > 1e-6
    assert bessel_reference.k_quadrature_step_halving_defect(240.0, 0) < 1e-13


def test_bracket_uses_the_magnitude_of_the_wavenumber() -> None:
    """NumPy stores negative frequencies; a signed bracket would be inward."""

    for k in (0.5, 3.0, 40.0):
        assert outer_bracket(-k, 1.1) == outer_bracket(k, 1.1)
        assert outer_bracket(k, 1.1) > 0.0


def test_bracket_variants_are_the_named_faults_and_nothing_else() -> None:
    assert BRACKET_VARIANTS == ("exact", *INJECTED_FAULTS)
    radius, k = 1.0, 2.0
    exact = outer_bracket(k, radius)
    assert outer_bracket(k, radius, variant="sign_flipped") == -exact
    assert outer_bracket(k, radius, variant="no_curvature_term") == pytest.approx(
        exact - 2.0 / radius, rel=1e-14
    )
    assert outer_bracket(k, radius, variant="frozen_ratio") == 2.0 / radius + k
    # The kR >> 1 misuse is an over-estimate for every finite kR.
    assert outer_bracket(k, radius, variant="frozen_ratio") > exact
    with pytest.raises(ValueError):
        outer_bracket(k, radius, variant="not_a_variant")
    with pytest.raises(ValueError):
        outer_bracket(k, 0.0)
    with pytest.raises(ValueError):
        outer_bracket(math.inf, 1.0)
    with pytest.raises(TypeError):
        outer_bracket(True, 1.0)


def test_frozen_ratio_is_bitwise_invisible_at_the_zero_mode() -> None:
    """A limitation the experiment reports rather than hides."""

    assert outer_bracket(0.0, 1.0, variant="frozen_ratio") == outer_bracket(0.0, 1.0)


# --------------------------------------------------------------------------
# 3. The Dirichlet path is bit-identical, and only the last row changes
# --------------------------------------------------------------------------


def _grid(nr: int, nz: int, r_max: float = 1.7, period: float = 2.0) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=r_max, z_min=0.0, z_max=period, periodic_z=True
    )


def _smooth_fields(grid: AxisymmetricGrid) -> tuple[np.ndarray, np.ndarray]:
    r, z = grid.mesh()
    q = 2.0 * np.pi * 2.0 / float(grid.z_period)
    omega = (1.0 - 0.4 * r**2) * np.cos(q * z) + 0.3 * r**2 - 0.11
    trace = 0.21 + 0.05 * np.cos(q * grid.z)
    return omega, trace


def test_dirichlet_path_is_bit_identical_to_solver_a() -> None:
    """The default option must cost the existing solver exactly nothing."""

    for nr, nz in ((17, 16), (33, 32)):
        grid = _grid(nr, nz)
        omega, trace = _smooth_fields(grid)
        for estimate in (False, True):
            expected = solve_streamfunction_poisson(
                grid, omega, trace, estimate_condition=estimate
            )
            actual = solve_streamfunction_poisson_outer(
                grid,
                omega,
                boundary_condition="dirichlet",
                outer_dirichlet=trace,
                estimate_condition=estimate,
            )
            assert np.array_equal(actual.psi1, expected.psi1)
            assert np.array_equal(
                actual.discrete_residual, expected.discrete_residual
            )
            assert np.array_equal(actual.pde_residual, expected.pde_residual)
            assert actual.metadata == expected.metadata
            # The default really is Dirichlet: omitting the keyword agrees.
            default = solve_streamfunction_poisson_outer(
                grid, omega, outer_dirichlet=trace, estimate_condition=estimate
            )
            assert np.array_equal(default.psi1, expected.psi1)


def test_dirichlet_path_rejects_transparent_only_arguments() -> None:
    grid = _grid(17, 16)
    omega, _ = _smooth_fields(grid)
    with pytest.raises(ValueError, match="support_radius is meaningful only"):
        solve_streamfunction_poisson_outer(
            grid, omega, boundary_condition="dirichlet", support_radius=1.0
        )
    with pytest.raises(ValueError, match="bracket_variant is meaningful only"):
        solve_streamfunction_poisson_outer(
            grid, omega, boundary_condition="dirichlet", bracket_variant="frozen_ratio"
        )


def test_flux_coefficients_extend_solver_a_without_touching_its_rows() -> None:
    """Interior rows and the axis row are solver A's, bitwise."""

    for nr in (17, 33, 65):
        grid = _grid(nr, 16)
        lower, upper = radial_flux_coefficients(grid)
        reference_lower, reference_upper = _radial_flux_coefficients(grid)
        assert lower.shape == (grid.nr,)
        assert np.array_equal(lower[:-1], reference_lower)
        assert np.array_equal(upper[:-1], reference_upper)
        # The axis row is the analytic limit with coefficient 8, unchanged.
        assert lower[0] == 0.0
        assert upper[0] == 8.0 / grid.dr**2

        # The appended row, from an independent evaluation of E-26 written out
        # here rather than taken from the module.
        radius = float(grid.r[-1])
        half_lower = radius - 0.5 * grid.dr
        half_upper = radius + 0.5 * grid.dr
        volume = (half_upper**4 - half_lower**4) / 4.0
        assert lower[-1] == pytest.approx(
            half_lower**3 / (grid.dr * volume), rel=1e-15
        )
        assert upper[-1] == pytest.approx(
            half_upper**3 / (grid.dr * volume), rel=1e-15
        )


def test_only_the_final_tridiagonal_row_differs_from_the_dirichlet_assembly() -> None:
    """The transparent option replaces one row and leaves the rest identical."""

    grid = _grid(33, 16)
    lower, upper = _radial_flux_coefficients(grid)
    for wave_number in (0.0, 1.7, -4.2):
        sub, diagonal, super_diagonal = assemble_radial_mode(grid, wave_number)
        dirichlet_diagonal = lower + upper + wave_number**2
        dirichlet_sub = -lower[1:]
        dirichlet_super = -upper[:-1]
        assert np.array_equal(diagonal[:-1], dirichlet_diagonal)
        assert np.array_equal(sub[:-1], dirichlet_sub)
        assert np.array_equal(super_diagonal[:-1], dirichlet_super)

        # The replaced row is (W-1h) with the ghost value eliminated.
        full_lower, full_upper = radial_flux_coefficients(grid)
        bracket = outer_bracket(wave_number, float(grid.r[-1]))
        assert sub[-1] == pytest.approx(-(full_lower[-1] + full_upper[-1]), rel=1e-15)
        assert diagonal[-1] == pytest.approx(
            full_lower[-1]
            + full_upper[-1]
            + wave_number**2
            + 2.0 * grid.dr * bracket * full_upper[-1],
            rel=1e-15,
        )
        # Strict diagonal dominance in the replaced row, which is what makes
        # the unpivoted Thomas elimination safe.
        assert diagonal[-1] > abs(sub[-1])


def test_transparent_solver_imports_nothing_unexpected() -> None:
    """The new module may reach solver A, the grid and the oracle -- no further."""

    tree = ast.parse(inspect.getsource(transparent_module))
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
        (1, ""),
        (1, "grid"),
        (1, "poisson"),
    }
    names = {name for _, name in imports}
    for forbidden in (
        "operators",
        "pde",
        "finite_cylinder_poisson",
        "realspace_poisson",
        "nonlinear_cylinder",
        "manufactured",
        "cartesian",
    ):
        assert forbidden not in names
        assert f"ns_certificate_lab.{forbidden}" not in names


# --------------------------------------------------------------------------
# 4. Solver guards
# --------------------------------------------------------------------------


def _compact_source(grid: AxisymmetricGrid, support_radius: float = 0.5) -> np.ndarray:
    profile = np.zeros(grid.nr, dtype=np.float64)
    inside = grid.r < support_radius
    scaled = grid.r[inside] / support_radius
    profile[inside] = (1.0 - scaled * scaled) ** 4
    axial = np.cos(2.0 * np.pi * grid.z / float(grid.z_period))
    return profile[:, None] * axial[None, :]


def test_transparent_path_rejects_bad_inputs() -> None:
    grid = _grid(17, 16, r_max=1.0)
    omega = _compact_source(grid)
    with pytest.raises(ValueError, match="requires an explicit"):
        solve_streamfunction_poisson_outer(grid, omega, boundary_condition="transparent")
    with pytest.raises(ValueError, match="homogeneous"):
        solve_streamfunction_poisson_outer(
            grid,
            omega,
            boundary_condition="transparent",
            outer_dirichlet=0.3,
            support_radius=0.5,
        )
    with pytest.raises(ValueError, match="must be real-valued"):
        solve_streamfunction_poisson_outer(
            grid,
            omega.astype(np.complex128),
            boundary_condition="transparent",
            support_radius=0.5,
        )
    broken = omega.copy()
    broken[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        solve_streamfunction_poisson_outer(
            grid, broken, boundary_condition="transparent", support_radius=0.5
        )
    with pytest.raises(ValueError, match="must have shape"):
        solve_streamfunction_poisson_outer(
            grid,
            omega[:-1],
            boundary_condition="transparent",
            support_radius=0.5,
        )
    with pytest.raises(ValueError, match="boundary_condition must be one of"):
        solve_streamfunction_poisson_outer(grid, omega, boundary_condition="robin")
    with pytest.raises(ValueError, match="bracket_variant must be one of"):
        solve_streamfunction_poisson_outer(
            grid,
            omega,
            boundary_condition="transparent",
            support_radius=0.5,
            bracket_variant="wrong",
        )
    with pytest.raises(TypeError):
        solve_streamfunction_poisson_outer(object(), omega)
    with pytest.raises(TypeError):
        solve_streamfunction_poisson_outer(
            grid, omega, estimate_condition="yes"  # type: ignore[arg-type]
        )
    nonperiodic = AxisymmetricGrid.uniform(
        nr=17, nz=16, r_max=1.0, z_min=0.0, z_max=1.0, periodic_z=False
    )
    with pytest.raises(ValueError, match="periodic"):
        solve_streamfunction_poisson_outer(
            nonperiodic,
            _compact_source(
                AxisymmetricGrid.uniform(
                    nr=17, nz=16, r_max=1.0, z_min=0.0, z_max=1.0, periodic_z=True
                )
            ),
            boundary_condition="transparent",
            support_radius=0.5,
        )


def test_transparent_metadata_records_the_choice_and_the_variant() -> None:
    grid = _grid(33, 16, r_max=1.0)
    omega = _compact_source(grid)
    solution = solve_streamfunction_poisson_outer(
        grid,
        omega,
        boundary_condition="transparent",
        support_radius=0.5,
        estimate_condition=False,
    )
    metadata = solution.metadata
    assert metadata["outer_boundary_condition"] == "transparent"
    assert metadata["bracket_variant"] == "exact"
    assert metadata["bracket_is_exact"] is True
    assert metadata["support_radius"] == 0.5
    assert metadata["wall_radius"] == 1.0
    assert metadata["zero_mode_bracket"] == 2.0
    assert metadata["zero_mode_bracket"] == metadata["zero_mode_bracket_closed_form"]
    assert "E-27" in metadata["outer_boundary_scope_warning"]
    assert metadata["axis_radial_coefficient_dimensionless"] == pytest.approx(8.0)
    # The interior discrete residual is the algebraic identity of the solve.
    assert metadata["discrete_residual_max_abs_interior"] < 1e-9
    assert metadata["discrete_residual_max_abs_boundary_row"] < 1e-9


def test_independent_condition_defect_is_second_order_and_not_an_identity() -> None:
    """The one-sided diagnostic is a different discretization, so it is not zero."""

    defects = []
    spacings = []
    for nr in (33, 65, 129):
        grid = _grid(nr, 16, r_max=1.0)
        omega = _compact_source(grid)
        solution = solve_streamfunction_poisson_outer(
            grid,
            omega,
            boundary_condition="transparent",
            support_radius=0.5,
            estimate_condition=False,
        )
        defect = transparent_condition_defect(grid, solution.psi1)
        assert np.array_equal(defect, solution.pde_residual[-1])
        defects.append(_max_abs(defect))
        spacings.append(grid.dr)
    assert min(defects) > 0.0
    for index in range(len(defects) - 1):
        assert 1.7 <= _order(defects[index], defects[index + 1]) <= 2.3


def test_compact_support_guard_refuses_and_reports() -> None:
    grid = _grid(33, 16, r_max=1.0)
    omega = _compact_source(grid, support_radius=0.5)
    assert assert_compact_support(grid, omega, support_radius=0.5) >= 1
    leaking = _compact_source(grid, support_radius=1.05)
    with pytest.raises(ValueError, match="not exactly zero outside"):
        assert_compact_support(grid, leaking, support_radius=0.5)
    with pytest.raises(ValueError, match="strictly inside the wall"):
        assert_compact_support(grid, omega, support_radius=1.0)
    with pytest.raises(ValueError, match="positive and finite"):
        assert_compact_support(grid, omega, support_radius=-1.0)
    with pytest.raises(TypeError):
        assert_compact_support(grid, omega, support_radius="0.5")


def test_modal_wavenumbers_match_solver_a_convention() -> None:
    grid = _grid(17, 16, period=3.0)
    expected = 2.0 * np.pi * np.fft.fftfreq(grid.nz, d=grid.dz)
    assert np.array_equal(modal_wavenumbers(grid), expected)


# --------------------------------------------------------------------------
# 5. The manufactured family
# --------------------------------------------------------------------------


def test_decaying_branch_derivatives_match_finite_differences() -> None:
    """The first line of this table is exactly where (W-1) comes from."""

    for k in (0.0, 0.5, 2.0 * math.pi):
        radius = 0.95
        step = 1e-4

        def value(x: float, wavenumber: float = k) -> float:
            return float(decaying_branch_derivatives(wavenumber, x)[0])

        first = (value(radius + step) - value(radius - step)) / (2.0 * step)
        second = (
            value(radius + step) - 2.0 * value(radius) + value(radius - step)
        ) / step**2
        exact = decaying_branch_derivatives(k, radius)
        assert first == pytest.approx(exact[1], rel=1e-6)
        assert second == pytest.approx(exact[2], rel=1e-6)
        # The branch is homogeneous: L_{5,k} u = 0.
        residual = exact[2] + 3.0 / radius * exact[1] - k * k * exact[0]
        assert abs(residual) <= 1e-10 * max(1.0, abs(exact[2]))


def test_manufactured_solution_satisfies_the_transparent_condition_exactly() -> None:
    """Outside its support the manufactured field *is* the decaying branch."""

    support = 0.95
    for k in (0.0, 2.0 * math.pi / 8.0, 2.0 * math.pi):
        coefficients = manufactured_polynomial_coefficients(k, support)
        assert coefficients.shape == (4,)
        for radius in (1.05, 1.25, 2.0):
            derivatives = decaying_branch_derivatives(k, radius)
            scale = decaying_branch_derivatives(k, support)[0]
            psi = derivatives[0] / scale
            slope = derivatives[1] / scale
            residual = slope + outer_bracket(k, radius) * psi
            assert abs(residual) <= 1e-12 * max(abs(slope), abs(psi))


def test_manufactured_source_is_exactly_zero_outside_its_support() -> None:
    grid = build_grid(
        points_per_unit_radius=64, wall_radius=1.25, axial_points=16, z_period=8.0
    )
    exact, omega = manufactured_case_fields(
        grid, support_radius=0.95, z_period=8.0, mode=1
    )
    assert _max_abs(omega[grid.r >= 0.95, :]) == 0.0
    assert _max_abs(omega) > 0.0
    # Continuity at the support edge: the last nonzero row is already small.
    inside = grid.r < 0.95
    assert abs(omega[inside][-1, 0]) < 0.05 * _max_abs(omega)
    assert _max_abs(exact) > 0.0


def test_manufactured_convergence_including_the_boundary_row() -> None:
    """W-B condition 3, run directly on the solver."""

    support = 0.95
    for period, mode in ((1.0, 0), (8.0, 1)):
        errors: list[float] = []
        boundary_errors: list[float] = []
        for ppur in (32, 64, 128):
            grid = build_grid(
                points_per_unit_radius=ppur,
                wall_radius=1.25,
                axial_points=16,
                z_period=period,
            )
            exact, omega = manufactured_case_fields(
                grid, support_radius=support, z_period=period, mode=mode
            )
            solution = solve_streamfunction_poisson_outer(
                grid,
                omega,
                boundary_condition="transparent",
                support_radius=support,
                estimate_condition=False,
            )
            errors.append(_max_abs(solution.psi1 - exact))
            boundary_errors.append(_max_abs(solution.psi1[-1] - exact[-1]))
        for series in (errors, boundary_errors):
            for index in range(len(series) - 1):
                assert 1.85 <= _order(series[index], series[index + 1]) <= 2.15
        # The boundary row is not accidentally exact; it carries real error.
        assert min(boundary_errors) > 0.0
        # A Dirichlet solve of the same problem does NOT converge to psi*,
        # which is what makes the comparison a test of the boundary row.
        grid = build_grid(
            points_per_unit_radius=128,
            wall_radius=1.25,
            axial_points=16,
            z_period=period,
        )
        exact, omega = manufactured_case_fields(
            grid, support_radius=support, z_period=period, mode=mode
        )
        dirichlet = solve_streamfunction_poisson(
            grid, omega, 0.0, estimate_condition=False
        ).psi1
        assert _max_abs(dirichlet - exact) > 100.0 * errors[-1]


# --------------------------------------------------------------------------
# 6. W-B conditions 1, 2 and 4 measured directly
# --------------------------------------------------------------------------


def _zero_mode_pair(ppur: int, radii: tuple[float, float]) -> dict[str, float]:
    support, exponent, core = 0.95, 8, 0.9
    transparent: dict[float, np.ndarray] = {}
    dirichlet: dict[float, np.ndarray] = {}
    for radius in radii:
        grid = build_grid(
            points_per_unit_radius=ppur,
            wall_radius=radius,
            axial_points=16,
            z_period=1.0,
        )
        omega = build_source(
            grid, support_radius=support, exponent=exponent, z_period=1.0, mode=0
        )
        mask = grid.r <= core + 1e-12
        transparent[radius] = solve_streamfunction_poisson_outer(
            grid,
            omega,
            boundary_condition="transparent",
            support_radius=support,
            estimate_condition=False,
        ).psi1[mask]
        dirichlet[radius] = solve_streamfunction_poisson(
            grid, omega, 0.0, estimate_condition=False
        ).psi1[mask]
    small, large = radii
    return {
        "transparent": _max_abs(transparent[small] - transparent[large]),
        "dirichlet": _max_abs(dirichlet[small] - dirichlet[large]),
        "dirichlet_signed": float(np.mean(dirichlet[small] - dirichlet[large])),
    }


def test_zero_mode_transparent_solution_is_r_independent_and_dirichlet_is_not() -> None:
    """W-B condition 1, including the audited Dirichlet number."""

    moment = source_moment(support_radius=0.95, exponent=8)
    closed_form = 0.5 * moment * (2.0**-2 - 1.0**-2)
    assert closed_form == pytest.approx(-1.6968880208e-3, rel=1e-9)

    measured: list[float] = []
    for ppur in (64, 128, 256):
        result = _zero_mode_pair(ppur, (1.0, 2.0))
        measured.append(result["transparent"])
        # The Dirichlet offset is the E-33(e) constant and does not shrink.
        assert result["dirichlet_signed"] == pytest.approx(closed_form, rel=1e-3)
        assert result["dirichlet"] > 1.6e-3
        assert result["transparent"] < 1.0e-5
    # ... while the transparent difference is O(dr^2) and vanishes.
    for index in range(len(measured) - 1):
        assert 1.85 <= _order(measured[index], measured[index + 1]) <= 2.15
    assert measured[-1] == pytest.approx(5.664e-8, rel=5e-2)
    # The finest Dirichlet value is the one recorded in docs/equation_audit.md.
    finest = _zero_mode_pair(256, (1.0, 2.0))
    assert finest["dirichlet_signed"] == pytest.approx(
        AUDITED_DIRICHLET_ZERO_MODE_OFFSET, rel=1e-6
    )


def test_long_wavelength_mode_is_where_the_transparent_condition_pays() -> None:
    """W-B condition 2: the ``kR << 1`` regime Dirichlet only reaches as R^-2."""

    support, exponent, core = 0.95, 8, 0.9
    period = 32.0  # kR = 0.196 at R = 1: deep in the algebraic regime
    ppur = 128
    transparent: dict[float, np.ndarray] = {}
    dirichlet: dict[float, np.ndarray] = {}
    for radius in (1.0, 2.0):
        grid = build_grid(
            points_per_unit_radius=ppur,
            wall_radius=radius,
            axial_points=16,
            z_period=period,
        )
        omega = build_source(
            grid, support_radius=support, exponent=exponent, z_period=period, mode=1
        )
        mask = grid.r <= core + 1e-12
        transparent[radius] = solve_streamfunction_poisson_outer(
            grid,
            omega,
            boundary_condition="transparent",
            support_radius=support,
            estimate_condition=False,
        ).psi1[mask]
        dirichlet[radius] = solve_streamfunction_poisson(
            grid, omega, 0.0, estimate_condition=False
        ).psi1[mask]
    transparent_gap = _max_abs(transparent[1.0] - transparent[2.0])
    dirichlet_gap = _max_abs(dirichlet[1.0] - dirichlet[2.0])
    assert dirichlet_gap > 1.0e-3
    assert transparent_gap < 1.0e-6
    assert dirichlet_gap / transparent_gap > 1.0e3


def test_transparent_small_radius_matches_a_large_radius_dirichlet_solve() -> None:
    """W-B condition 4, with the two terms of the split separated."""

    support, exponent, core = 0.95, 8, 0.9
    moment = source_moment(support_radius=support, exponent=exponent)
    reference_radius = 8.0
    wall_term = moment / (2.0 * reference_radius**2)
    residuals: list[float] = []
    for ppur in (64, 128):
        small_grid = build_grid(
            points_per_unit_radius=ppur,
            wall_radius=1.0,
            axial_points=16,
            z_period=1.0,
        )
        small_omega = build_source(
            small_grid,
            support_radius=support,
            exponent=exponent,
            z_period=1.0,
            mode=0,
        )
        transparent = solve_streamfunction_poisson_outer(
            small_grid,
            small_omega,
            boundary_condition="transparent",
            support_radius=support,
            estimate_condition=False,
        ).psi1[small_grid.r <= core + 1e-12]
        big_grid = build_grid(
            points_per_unit_radius=ppur,
            wall_radius=reference_radius,
            axial_points=16,
            z_period=1.0,
        )
        big_omega = build_source(
            big_grid,
            support_radius=support,
            exponent=exponent,
            z_period=1.0,
            mode=0,
        )
        dirichlet = solve_streamfunction_poisson(
            big_grid, big_omega, 0.0, estimate_condition=False
        ).psi1[big_grid.r <= core + 1e-12]
        total = _max_abs(transparent - dirichlet)
        residual = _max_abs(transparent - dirichlet - wall_term)
        # The total is dominated by the reference wall's own R^-2 error ...
        assert total == pytest.approx(wall_term, rel=0.15)
        # ... and what is left after removing that exact term is O(dr^2).
        assert residual < 0.3 * total
        residuals.append(residual)
    assert 1.85 <= _order(residuals[0], residuals[1]) <= 2.15


# --------------------------------------------------------------------------
# 7. Fault injection (W-B condition 5)
# --------------------------------------------------------------------------


def _fault_r_dependence(variant: str, period: float, mode: int) -> float:
    support, exponent, core = 0.95, 8, 0.9
    cores = []
    for radius in (1.0, 2.0):
        grid = build_grid(
            points_per_unit_radius=64,
            wall_radius=radius,
            axial_points=16,
            z_period=period,
        )
        omega = build_source(
            grid,
            support_radius=support,
            exponent=exponent,
            z_period=period,
            mode=mode,
        )
        cores.append(
            solve_streamfunction_poisson_outer(
                grid,
                omega,
                boundary_condition="transparent",
                support_radius=support,
                bracket_variant=variant,
                estimate_condition=False,
            ).psi1[grid.r <= core + 1e-12]
        )
    return _max_abs(cores[1] - cores[0])


def test_detects_sign_flipped_bracket() -> None:
    """Detected by the R-independence check, by three orders of magnitude."""

    for period, mode in ((1.0, 0), (8.0, 1)):
        baseline = _fault_r_dependence("exact", period, mode)
        faulted = _fault_r_dependence("sign_flipped", period, mode)
        assert faulted / baseline > 1.0e3


def test_detects_dropped_two_over_r_term() -> None:
    """Detected as a refusal: the k=0 modal matrix becomes singular.

    Dropping ``2/R`` leaves ``beta_0 = 0``, i.e. a pure discrete Neumann
    problem whose null space is the constant.  Since the source has a nonzero
    ``r^3``-weighted moment the system is inconsistent, and the unpivoted
    elimination stops on the final pivot.  Every axial case is affected,
    because the ``k = 0`` mode is present in every FFT.
    """

    for period, mode in ((1.0, 0), (8.0, 1)):
        with pytest.raises(np.linalg.LinAlgError, match="pivot"):
            _fault_r_dependence("no_curvature_term", period, mode)

    # Restricted to a k > 0 mode the matrix is still solvable, and there the
    # fault shows up as a number: the R-dependence grows by a large factor.
    support, exponent, core = 0.95, 8, 0.9
    k = 2.0 * math.pi / 8.0
    measurements: dict[str, float] = {}
    for variant in ("exact", "no_curvature_term"):
        cores = []
        for radius in (1.0, 2.0):
            grid = build_grid(
                points_per_unit_radius=64,
                wall_radius=radius,
                axial_points=16,
                z_period=8.0,
            )
            profile = experiment.radial_profile(
                grid.r, support_radius=support, exponent=exponent
            ).astype(np.complex128)
            solved = transparent_module.solve_radial_mode(
                grid, profile, k, bracket_variant=variant
            )
            cores.append(np.real(solved)[grid.r <= core + 1e-12])
        measurements[variant] = _max_abs(cores[1] - cores[0])
    assert measurements["no_curvature_term"] / measurements["exact"] > 100.0


def test_detects_k0_over_k1_frozen_to_one() -> None:
    """Detected on every k>0 case; invisible at k=0, and that is reported."""

    baseline = _fault_r_dependence("exact", 8.0, 1)
    faulted = _fault_r_dependence("frozen_ratio", 8.0, 1)
    assert faulted / baseline > 10.0

    baseline_long = _fault_r_dependence("exact", 32.0, 1)
    faulted_long = _fault_r_dependence("frozen_ratio", 32.0, 1)
    assert faulted_long / baseline_long > 10.0

    # At k = 0 the bracket is 2/R for both, so the solve is bitwise identical.
    # This is a real limitation of the fault, not a failure of the detector.
    assert _fault_r_dependence("frozen_ratio", 1.0, 0) == _fault_r_dependence(
        "exact", 1.0, 0
    )


def test_fault_variants_also_break_the_manufactured_convergence() -> None:
    """A second, independent place each fault is caught."""

    support = 0.95
    errors: dict[str, list[float]] = {}
    for variant in ("exact", "sign_flipped", "frozen_ratio"):
        errors[variant] = []
        for ppur in (64, 128):
            grid = build_grid(
                points_per_unit_radius=ppur,
                wall_radius=1.25,
                axial_points=16,
                z_period=8.0,
            )
            exact, omega = manufactured_case_fields(
                grid, support_radius=support, z_period=8.0, mode=1
            )
            solution = solve_streamfunction_poisson_outer(
                grid,
                omega,
                boundary_condition="transparent",
                support_radius=support,
                bracket_variant=variant,
                estimate_condition=False,
            )
            errors[variant].append(_max_abs(solution.psi1 - exact))
    assert 1.85 <= _order(*errors["exact"]) <= 2.15
    for variant in ("sign_flipped", "frozen_ratio"):
        assert _order(*errors[variant]) < 1.0
        assert errors[variant][-1] > 10.0 * errors["exact"][-1]


# --------------------------------------------------------------------------
# 8. Support-leak refusal (W-B condition 6)
# --------------------------------------------------------------------------


def test_support_leak_is_refused_and_a_compliant_source_is_not() -> None:
    config = _reduced_config()
    outcomes = measure_support_leak(config)
    assert outcomes["declared_support_reaches_wall_refused"] is True
    assert outcomes["undeclared_residue_refused"] is True
    assert outcomes["compliant_source_accepted"] is True
    assert outcomes["all_refusals_fired"] is True
    assert outcomes["max_abs_source_at_and_beyond_wall"] > 0.0
    assert "strictly inside the wall" in outcomes[
        "declared_support_reaches_wall_message"
    ]
    assert "not exactly zero outside" in outcomes["undeclared_residue_message"]


# --------------------------------------------------------------------------
# 9. Configuration
# --------------------------------------------------------------------------


def _reduced_config() -> dict[str, Any]:
    """A small but structurally complete config for the smoke run."""

    return {
        "schema_version": 1,
        "experiment": EXPERIMENT_ID,
        "description": "reduced transparent-boundary configuration for tests",
        "interpretation": (
            "smoke configuration; the shipped config carries the evidence"
        ),
        "physical_wall_warning": (
            "whole-space condition only; the Hou reproduction keeps the E-27 "
            "physical no-slip wall"
        ),
        "source_profile": {"support_radius": 0.95, "exponent": 8},
        "core_radius": 0.9,
        "axial_cases": [
            {"label": "k0", "z_period": 1.0, "mode": 0},
            {"label": "Lz32_m1", "z_period": 32.0, "mode": 1},
        ],
        "wall_radii": [1.0, 2.0],
        "axial_points": 16,
        "radial_resolutions": [32, 64, 128],
        "manufactured": {
            "wall_radius": 1.25,
            "radial_resolutions": [32, 64, 128],
            "axial_case_labels": ["k0", "Lz32_m1"],
        },
        "large_radius_dirichlet": {
            "transparent_wall_radius": 1.0,
            "dirichlet_wall_radii": [2.0, 8.0],
        },
        "fault_injection": {
            "points_per_unit_radius": 32,
            "wall_radii": [1.0, 2.0],
            "variants": list(INJECTED_FAULTS),
        },
        "support_leak": {
            "points_per_unit_radius": 32,
            "wall_radius": 1.0,
            "leaking_support_radius": 1.05,
        },
        "long_wavelength_max_kr": 0.5,
        "acceptance": {
            "max_transparent_r_independence_over_dr_squared": 0.05,
            "min_transparent_r_independence_observed_order": 1.85,
            "max_transparent_r_independence_observed_order": 2.15,
            "max_zero_mode_dirichlet_closed_form_relative_error": 1e-3,
            "min_dirichlet_over_transparent_improvement_factor": 100.0,
            "min_long_wavelength_improvement_factor": 1000.0,
            "min_manufactured_observed_order": 1.85,
            "max_manufactured_observed_order": 2.15,
            "min_manufactured_boundary_row_observed_order": 1.85,
            "max_manufactured_boundary_row_observed_order": 2.15,
            "max_manufactured_error_over_dr_squared": 100.0,
            "large_radius_discretization_constant": 0.05,
            "max_large_radius_dirichlet_split_factor": 1.5,
            "max_zero_mode_wall_corrected_residual_over_dr_squared": 0.05,
            "min_zero_mode_wall_corrected_observed_order": 1.85,
            "max_zero_mode_wall_corrected_observed_order": 2.15,
            "min_fault_detection_factor": 10.0,
            "max_outer_condition_defect_over_dr_squared": 0.5,
            "max_k0_over_k1_published_relative_error": 1e-12,
            "max_k0_over_k1_derivative_identity_relative_defect": 1e-9,
            "max_wronskian_relative_defect": 1e-14,
            "max_k0_over_k1_small_argument_relative_error": 1e-3,
            "max_k0_over_k1_large_argument_relative_error": 1e-4,
            "max_zero_mode_bracket_relative_defect": 1e-15,
        },
    }


def _shipped_config() -> dict[str, Any]:
    return json.loads(SHIPPED_CONFIG.read_text(encoding="utf-8"))


def test_shipped_config_is_valid_and_canonical() -> None:
    config = _shipped_config()
    assert validate_config(copy.deepcopy(config)) is not None
    assert config["experiment"] == EXPERIMENT_ID
    assert "E-27" in config["physical_wall_warning"]
    assert set(config["acceptance"]) == set(_reduced_config()["acceptance"])
    modes = {case["mode"] for case in config["axial_cases"]}
    assert 0 in modes and len(modes) > 1
    # The shipped sweep must reach the algebraic regime, where the Dirichlet
    # wall error decays only like R^-2 and the transparent condition matters.
    smallest_wall = min(config["wall_radii"])
    long_wavelength = [
        case
        for case in config["axial_cases"]
        if case["mode"]
        and 2.0 * math.pi * case["mode"] / case["z_period"] * smallest_wall
        <= config["long_wavelength_max_kr"]
    ]
    assert long_wavelength


def test_reduced_config_is_valid() -> None:
    assert validate_config(_reduced_config()) is not None


@pytest.mark.parametrize("key", sorted(experiment.TOP_LEVEL_KEYS))
def test_config_rejects_every_missing_key(key: str) -> None:
    config = _reduced_config()
    del config[key]
    with pytest.raises(ValueError):
        validate_config(config)


def test_config_rejects_unknown_key() -> None:
    config = _reduced_config()
    config["surprise"] = 1
    with pytest.raises(ValueError, match="unknown keys"):
        validate_config(config)


@pytest.mark.parametrize(
    "section",
    [
        "source_profile",
        "manufactured",
        "large_radius_dirichlet",
        "fault_injection",
        "support_leak",
        "acceptance",
    ],
)
def test_config_rejects_unknown_nested_keys(section: str) -> None:
    config = _reduced_config()
    config[section]["surprise"] = 1
    with pytest.raises(ValueError, match="unknown keys"):
        validate_config(config)


def test_config_rejects_structural_faults() -> None:
    cases: list[tuple[str, Any]] = [
        ("experiment", "something_else"),
        ("schema_version", 2),
        ("physical_wall_warning", "a warning without the audit reference"),
        ("wall_radii", [2.0, 1.0]),
        ("wall_radii", [1.0]),
        ("core_radius", 0.99),
        ("axial_points", 15),
        ("radial_resolutions", [64, 128]),
        ("radial_resolutions", [128, 64, 32]),
        ("long_wavelength_max_kr", 1.5),
    ]
    for key, value in cases:
        config = _reduced_config()
        config[key] = value
        with pytest.raises(ValueError):
            validate_config(config)


def test_config_rejects_vacuous_acceptance_tolerances() -> None:
    config = _reduced_config()
    config["acceptance"]["min_fault_detection_factor"] = 1.0
    with pytest.raises(ValueError, match="must exceed 1"):
        validate_config(config)

    config = _reduced_config()
    config["acceptance"]["large_radius_discretization_constant"] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        validate_config(config)


def test_config_rejects_a_source_that_reaches_the_smallest_wall() -> None:
    config = _reduced_config()
    config["source_profile"]["support_radius"] = 1.0
    with pytest.raises(ValueError, match="strictly smaller than the smallest wall"):
        validate_config(config)


def test_config_rejects_bad_axial_case_sets() -> None:
    config = _reduced_config()
    config["axial_cases"] = [{"label": "a", "z_period": 1.0, "mode": 1}]
    with pytest.raises(ValueError, match="exactly one axial case"):
        validate_config(config)

    config = _reduced_config()
    config["axial_cases"] = [
        {"label": "k0", "z_period": 1.0, "mode": 0},
        {"label": "k0", "z_period": 2.0, "mode": 0},
    ]
    with pytest.raises(ValueError):
        validate_config(config)


def test_config_rejects_bad_subsection_geometry() -> None:
    config = _reduced_config()
    config["manufactured"]["wall_radius"] = 0.5
    with pytest.raises(ValueError, match="strictly outside the source support"):
        validate_config(config)

    config = _reduced_config()
    config["large_radius_dirichlet"]["dirichlet_wall_radii"] = [0.5, 8.0]
    with pytest.raises(ValueError, match="must exceed the transparent"):
        validate_config(config)

    config = _reduced_config()
    config["large_radius_dirichlet"]["transparent_wall_radius"] = 1.7
    with pytest.raises(ValueError, match="must be one of"):
        validate_config(config)

    config = _reduced_config()
    config["fault_injection"]["wall_radii"] = [1.0, 3.0]
    with pytest.raises(ValueError, match="subset of wall_radii"):
        validate_config(config)

    config = _reduced_config()
    config["support_leak"]["leaking_support_radius"] = 0.5
    with pytest.raises(ValueError, match="must reach the wall"):
        validate_config(config)


def test_config_cannot_silently_drop_an_injected_fault() -> None:
    config = _reduced_config()
    config["fault_injection"]["variants"] = ["sign_flipped", "frozen_ratio"]
    with pytest.raises(ValueError, match="W-B condition 5 names all three"):
        validate_config(config)

    config = _reduced_config()
    config["fault_injection"]["variants"] = ["exact", "sign_flipped", "frozen_ratio"]
    with pytest.raises(ValueError):
        validate_config(config)


def test_radial_point_count_refuses_fractional_cells() -> None:
    assert radial_point_count(64, 1.5) == 97
    with pytest.raises(ValueError, match="whole cells"):
        radial_point_count(64, 1.01)


# --------------------------------------------------------------------------
# 10. End-to-end bundle
# --------------------------------------------------------------------------


def test_smoke_run_writes_a_verifiable_bundle(scratch_dir: Path) -> None:
    config = _reduced_config()
    summary = run(config, scratch_dir / "bundle")
    output = scratch_dir / "bundle"

    assert summary["accepted"] is True
    assert summary["experiment_id"] == EXPERIMENT_ID
    assert summary["outer_boundary_condition"] == "transparent"
    assert summary["outer_boundary_condition_baseline"] == "dirichlet"
    assert "E-27" in summary["physical_wall_warning"]
    assert set(summary["conditions_tested"]) == {
        "W-B.1",
        "W-B.2",
        "W-B.3",
        "W-B.4",
        "W-B.5",
        "W-B.6",
    }
    assert summary["limitations"] and summary["known_gaps"]

    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert set(manifest["files"]) == {
        "summary.json",
        "config.snapshot.json",
        "r_independence.csv",
        "manufactured_convergence.csv",
        "fault_injection.csv",
        "large_radius_agreement.csv",
    }
    for name, entry in manifest["files"].items():
        path = output / name
        assert _sha256(path) == entry["sha256"]
        assert path.stat().st_size == entry["bytes"]
        sidecar = output / (name + ".sha256")
        assert sidecar.read_text(encoding="ascii").strip() == entry["sha256"]
    assert (output / "manifest.json.sha256").read_text(
        encoding="ascii"
    ).strip() == _sha256(output / "manifest.json")

    snapshot = json.loads((output / "config.snapshot.json").read_text(encoding="utf-8"))
    assert snapshot == config
    provenance = summary["reproducibility"]["runtime_provenance"]
    assert len(provenance["source_fingerprint_sha256"]) == 64


def test_smoke_run_measures_every_wb_condition(scratch_dir: Path) -> None:
    config = _reduced_config()
    summary = run(config, scratch_dir / "bundle")
    checks = summary["acceptance"]["checks"]
    assert checks["all_passed"] is True
    for required in (
        "transparent_r_independence_at_discretization_level",
        "transparent_r_independence_orders_in_band",
        "dirichlet_zero_mode_offset_matches_e33e_closed_form",
        "transparent_beats_dirichlet_everywhere",
        "long_wavelength_improvement_is_large",
        "manufactured_orders_in_band",
        "manufactured_boundary_row_orders_in_band",
        "large_radius_agreement_within_predicted_split",
        "zero_mode_wall_corrected_residual_is_second_order",
        "every_injected_fault_is_detected",
        "every_injected_fault_is_detected_strongly",
        "support_leak_is_refused",
        "small_kr_bracket_tends_to_two_over_r",
        "large_kr_bracket_tends_to_two_over_r_plus_k",
    ):
        assert checks[required] is True, required

    analysis = summary["analysis"]
    verdicts = analysis["fault_verdicts"]
    for variant in INJECTED_FAULTS:
        assert verdicts[variant]["all_altered_cases_detected"] is True
    assert verdicts["frozen_ratio"]["cases_whose_bracket_is_bitwise_unchanged"] == [
        "k0"
    ]
    assert verdicts["no_curvature_term"]["cases_refused_by_the_solver"]
    limits = analysis["bracket_limits"]
    assert limits["zero_mode_bracket_is_bitwise_two_over_r"] is True
    assert limits["zero_mode_decaying_solution_residual"] == 0.0


def test_run_refuses_to_overwrite_a_nonempty_output_directory(
    scratch_dir: Path,
) -> None:
    output = scratch_dir / "bundle"
    output.mkdir()
    (output / "existing.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(FileExistsError, match="nonempty"):
        run(_reduced_config(), output)
    assert (output / "existing.txt").read_text(encoding="utf-8") == "keep me"


def test_run_validates_before_creating_any_output(scratch_dir: Path) -> None:
    config = _reduced_config()
    del config["core_radius"]
    output = scratch_dir / "bundle"
    with pytest.raises(ValueError):
        run(config, output)
    assert not output.exists()


def test_acceptance_fails_loudly_when_a_tolerance_is_violated() -> None:
    """The gate must be capable of failing; a green run is not automatic."""

    config = _reduced_config()
    analysis = evaluate(config)
    baseline = acceptance_report(analysis, config)
    assert baseline["checks"]["all_passed"] is True

    tightened = copy.deepcopy(config)
    tightened["acceptance"][
        "max_transparent_r_independence_over_dr_squared"
    ] = 1.0e-12
    report = acceptance_report(analysis, tightened)
    assert report["checks"]["transparent_r_independence_at_discretization_level"] is (
        False
    )
    assert report["checks"]["all_passed"] is False
