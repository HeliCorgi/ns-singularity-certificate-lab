"""Verification of solver C, the real-space finite-cylinder Poisson path.

Solver A is :mod:`ns_certificate_lab.poisson` (``r^3``-flux finite volume,
Fourier-diagonalized ``z``, Thomas elimination).  Solver C is
:mod:`ns_certificate_lab.realspace_poisson`: the same ``r^3``-flux radial
mathematics of E-26, transcribed independently, but with a three-point periodic
finite difference in ``z`` and a matrix-free preconditioned conjugate-gradient
inversion.  It uses no transform and no tridiagonal elimination, so it breaks
the axial common mode that solvers A and B share.

This module imports solver A for cross-checks and never imports the
direct-stencil sibling (solver B); the pairing guard in
``tests/test_poisson_cross_validation.py`` therefore stays green without being
weakened.

Every manufactured field below is written out in closed form in this file, and
every convergence order is recomputed here from the measured errors.  No
production code is used to derive an oracle.

Measured constants quoted in the assertions were obtained on this machine with
the pinned environment; the assertions keep at least an order of magnitude of
margin so that they pin behavior rather than arithmetic noise.
"""

from __future__ import annotations

import ast
import inspect
import math
from pathlib import Path

import numpy as np
import pytest

import ns_certificate_lab.realspace_poisson as realspace_module
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.poisson import solve_streamfunction_poisson
from ns_certificate_lab.realspace_poisson import (
    ConvergenceError,
    build_realspace_operator,
    solve_realspace_poisson,
)

R_MAX = 1.7
Z_PERIOD = 2.0 * np.pi
RESOLUTIONS = ((17, 32), (33, 64), (65, 128))


def _grid(nr: int, nz: int) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=R_MAX,
        z_min=0.0,
        z_max=Z_PERIOD,
        periodic_z=True,
    )


def _manufactured(
    grid: AxisymmetricGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``psi*``, ``omega* = -L5 psi*`` and the outer trace ``g(z)``.

    The analytic identity is ``L5(r**(2*k)) = 4*k*(k+1) * r**(2*k-2)``, since

        d_rr r**(2k) + (3/r) d_r r**(2k)
            = 2k(2k-1) r**(2k-2) + 6k r**(2k-2)
            = (4k**2 + 4k) r**(2k-2),

    so ``L5(r**2) = 8`` and ``L5(r**4) = 24 r**2``.  With

        psi*(r,z) = (1.1 + a r^2 + b r^4) + (c + d r^2 + e r^4) cos(q z),

    the radial part of ``L5 psi*`` is ``8a + 24 b r^2`` plus
    ``(8d + 24 e r^2) cos(q z)``, and ``d_zz`` contributes
    ``-q^2 (c + d r^2 + e r^4) cos(q z)``.  Hence

        omega* = -L5 psi*
               = -(8a + 24 b r^2)
                 + (q^2 (c + d r^2 + e r^4) - 8d - 24 e r^2) cos(q z).

    The outer trace is genuinely ``z``-dependent, so the Dirichlet elimination
    is exercised rather than trivially satisfied.
    """

    r, z = grid.mesh()
    q = 2.0 * np.pi * 3.0 / float(grid.z_period)
    a, b = -0.31, 0.07
    c, d, e = 0.22, -0.19, 0.035
    radial_constant = 1.1 + a * r**2 + b * r**4
    radial_mode = c + d * r**2 + e * r**4
    cosine = np.cos(q * z)

    psi = radial_constant + radial_mode * cosine
    omega = -(8.0 * a + 24.0 * b * r**2) + (
        q * q * radial_mode - 8.0 * d - 24.0 * e * r**2
    ) * cosine
    return psi, omega, psi[-1, :].copy()


def _radially_exact(
    grid: AxisymmetricGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a field whose radial profiles the ``r^3``-flux stencil is exact on.

    Only ``a + b r^2`` radial profiles appear, and the flux stencil reproduces
    ``L5(a + b r^2) = 8 b`` exactly (the telescoping face fluxes give
    ``2(r_{i+1/2}^4 - r_{i-1/2}^4) = 8 V_i``).  The axial factor is a single
    grid-resolved Fourier mode.  Solver A is therefore exact on this field to
    roundoff, so any solver-C error here is *entirely* the axial
    finite-difference truncation.
    """

    r, z = grid.mesh()
    q = 2.0 * np.pi * 3.0 / float(grid.z_period)
    radial_mode = 0.9 - 0.35 * r**2
    radial_constant = 0.4 + 0.21 * r**2
    axial = np.cos(q * z)

    psi = radial_mode * axial + radial_constant
    # L5 psi* = 8*(-0.35)*axial + 8*0.21 + radial_mode * d_zz(axial).
    l5_psi = 8.0 * (-0.35) * axial + 8.0 * 0.21 - q * q * radial_mode * axial
    return psi, -l5_psi, psi[-1, :].copy()


def _solve_a(
    grid: AxisymmetricGrid,
    omega: np.ndarray,
    boundary: np.ndarray,
) -> np.ndarray:
    return solve_streamfunction_poisson(
        grid,
        omega,
        boundary,
        estimate_condition=False,
    ).psi1


def _order(
    coarse_error: float,
    fine_error: float,
    coarse_h: float,
    fine_h: float,
) -> float:
    return math.log(coarse_error / fine_error) / math.log(coarse_h / fine_h)


def _max_abs(values: np.ndarray) -> float:
    return float(np.max(np.abs(values)))


def test_manufactured_solution_converges_at_second_order() -> None:
    """Solver C is second-order accurate against a closed-form solution.

    Measured on this machine: errors 9.319e-03, 2.330e-03, 5.826e-04 with
    observed orders 2.0000 and 1.9997, and iteration counts 34, 73, 146.
    """

    spacings: list[float] = []
    errors: list[float] = []
    iterations: list[int] = []

    for nr, nz in RESOLUTIONS:
        grid = _grid(nr, nz)
        exact, omega, boundary = _manufactured(grid)
        # The Dirichlet trace is genuinely z-dependent: peak-to-peak amplitude
        # 2*|0.22 - 0.19*R^2 + 0.035*R^4| = 0.0736 at R = 1.7.
        assert float(np.max(boundary) - np.min(boundary)) > 0.05

        result = solve_realspace_poisson(grid, omega, boundary)

        # The Dirichlet row is written back exactly, not merely approximated.
        assert _max_abs(result.psi1[-1] - boundary) == 0.0
        # The reported weighted residual honours the requested tolerance.
        assert result.weighted_relative_residual <= 1.0e-12
        # The unweighted maximum-norm residual is a much weaker statement than
        # the weighted one because the axis cell has volume dr^4/64; it is
        # reported separately and is checked here so that a weighted residual
        # cannot hide a bad axis row.  Measured: 1.3e-11, 6.2e-11, 6.6e-11.
        assert (
            result.metadata["algebraic_residual_max_abs_relative"] <= 1.0e-8
        )
        assert 0 < result.iterations <= 400

        spacings.append(grid.dr)
        errors.append(_max_abs(result.psi1 - exact))
        iterations.append(result.iterations)

    orders = [
        _order(errors[i], errors[i + 1], spacings[i], spacings[i + 1])
        for i in range(len(errors) - 1)
    ]
    assert min(orders) >= 1.80
    assert max(orders) <= 2.20
    # Measured error/dr^2 at the finest grid: 0.826.
    assert errors[-1] <= 2.0 * spacings[-1] ** 2
    # Performance pin: Jacobi-preconditioned CG must not silently degrade into
    # an O(N^2) iteration count.  Measured 34, 73, 146.
    assert iterations[-1] <= 400


def test_difference_from_solver_a_is_the_axial_discretization_gap() -> None:
    """Solver C must differ from solver A by ``O(dz^2)`` -- no more, no less.

    The two solvers share the ``r^3``-flux radial mathematics (independently
    transcribed), so their radial truncation errors cancel in the comparison.
    What does not cancel is the axial symbol: solver A applies the exact ``k^2``
    to each mode, solver C applies ``(4/dz^2) sin^2(k dz/2) = k^2 - k^4 dz^2/12
    + ...``.  The difference must therefore be a genuine ``O(dz^2)`` gap.

    Measured constants on this machine: ``D(h)/dz^2`` = 0.0912, 0.0903, 0.0901,
    i.e. a clean constant of about 0.09, with observed orders 2.014 and 2.004.
    The lower bound below is set at 1e-3 (about 90x below the measured
    constant) and the upper bound at 1.0 (about 11x above it).
    """

    axial_spacings: list[float] = []
    differences: list[float] = []

    for nr, nz in RESOLUTIONS:
        grid = _grid(nr, nz)
        _, omega, boundary = _manufactured(grid)

        psi_c = solve_realspace_poisson(grid, omega, boundary).psi1
        psi_a = _solve_a(grid, omega, boundary)

        axial_spacings.append(grid.dz)
        differences.append(_max_abs(psi_c - psi_a))

        # The zero axial mode sees no symbol gap at all: for k = 0 the exact
        # symbol and the finite-difference symbol are both zero.  The z-average
        # of the two solutions is therefore governed purely by the shared
        # radial stencil, and must agree to roundoff.  This is the sharpest
        # available statement of exactly what solver C does and does not share
        # with solver A.  Measured: 3.3e-16, 1.4e-14, 5.2e-15.
        assert _max_abs(psi_c.mean(axis=1) - psi_a.mean(axis=1)) <= 1.0e-12

    for spacing, difference in zip(axial_spacings, differences):
        # Upper bound: the gap is only the axial truncation error.
        assert difference <= 1.0 * spacing**2
        # CRITICAL LOWER BOUND.  Solver C exists to break the axial common mode
        # that solvers A and B share.  If a future refactor ever routes solver
        # C's axial term back through a transform (or through solver A), this
        # difference would collapse to roundoff and the third path would stop
        # being a third path.  Do not weaken this bound to make a refactor
        # green: an intentional merge must be reviewed as the loss of the
        # independent axial audit that it is.
        assert difference >= 1.0e-3 * spacing**2

    orders = [
        _order(
            differences[i],
            differences[i + 1],
            axial_spacings[i],
            axial_spacings[i + 1],
        )
        for i in range(len(differences) - 1)
    ]
    assert min(orders) >= 1.80
    assert max(orders) <= 2.20


def test_radially_exact_field_isolates_the_axial_path() -> None:
    """On a radially exact field, solver A is exact and solver C is not.

    This separates the two claims cleanly.  Solver A reproduces this field to
    roundoff because its radial stencil is exact on ``a + b r^2`` profiles and
    its axial treatment is exact on grid-resolved modes.  Solver C shares the
    first property but not the second, so its whole error here is the axial
    finite difference: measured 1.836e-02 and 4.559e-03, i.e. 0.476 dz^2 and
    0.473 dz^2, with observed order 2.01.
    """

    spacings: list[float] = []
    errors_c: list[float] = []

    for nr, nz in ((17, 32), (33, 64)):
        grid = _grid(nr, nz)
        exact, omega, boundary = _radially_exact(grid)

        psi_c = solve_realspace_poisson(grid, omega, boundary).psi1
        psi_a = _solve_a(grid, omega, boundary)

        # Solver A is exact here; this pins the manufactured field itself.
        assert _max_abs(psi_a - exact) <= 1.0e-12
        error_c = _max_abs(psi_c - exact)
        # Solver C's entire error is the axial truncation.
        assert 0.1 * grid.dz**2 <= error_c <= 2.0 * grid.dz**2
        assert abs(error_c - _max_abs(psi_c - psi_a)) <= 1.0e-12

        spacings.append(grid.dz)
        errors_c.append(error_c)

    order = _order(errors_c[0], errors_c[1], spacings[0], spacings[1])
    assert 1.80 <= order <= 2.20


def test_reduced_operator_is_symmetric_in_the_weighted_inner_product() -> None:
    """``<A u, v>_w == <u, A v>_w`` for random ``u`` and ``v``.

    The flux form guarantees ``a_i^+ V_i = F_i = a_{i+1}^- V_{i+1}``, so the
    weighted bilinear form telescopes into a sum over faces and is symmetric.
    The implementation stores one array of face coefficients ``F``, so the
    identity holds as an equality of the same floating-point number rather than
    only to rounding.  Measured relative asymmetry: 1.9e-16 and 0.0.
    """

    generator = np.random.default_rng(20260728)
    for nr, nz in ((17, 32), (33, 64)):
        grid = _grid(nr, nz)
        operator = build_realspace_operator(grid)
        assert operator.shape == (nr - 1, nz)

        for _ in range(3):
            u = generator.standard_normal(operator.shape)
            v = generator.standard_normal(operator.shape)

            left = operator.inner(operator.apply(u), v)
            right = operator.inner(u, operator.apply(v))
            scale = max(abs(left), abs(right))
            assert scale > 0.0
            assert abs(left - right) <= 1.0e-11 * scale

            # Positive definiteness on the same random draws.
            assert operator.inner(operator.apply(u), u) > 0.0

        # The face identity that produces the symmetry, checked directly.
        # a_i^+ V_i and a_{i+1}^- V_{i+1} are both F_i.
        upper_coefficient = operator.face_coefficient / operator.volume
        lower_coefficient = np.concatenate(
            ([0.0], operator.face_coefficient[:-1] / operator.volume[1:])
        )
        np.testing.assert_allclose(
            upper_coefficient[:-1] * operator.volume[:-1],
            lower_coefficient[1:] * operator.volume[1:],
            rtol=1.0e-14,
            atol=0.0,
        )

        # The zero vector is the only null vector: A is positive definite, so
        # a nonzero constant block has strictly positive energy (it is pinned
        # by the eliminated Dirichlet row).
        constant = np.ones(operator.shape)
        assert operator.inner(operator.apply(constant), constant) > 0.0


def test_axis_row_reproduces_the_e26b_coefficient() -> None:
    """The applied axis action must be exactly ``8 (psi_1 - psi_0) / dr^2``.

    The module never hard-codes ``8``; it evaluates the E-26a flux formula with
    ``r_{-1/2} = 0`` and ``V_0 = (dr/2)^4/4``.  This test probes the matrix-free
    apply with coordinate basis vectors, which is the operator that conjugate
    gradients actually inverts, not a separately assembled copy of it.
    """

    for nr, nz in ((17, 32), (33, 64)):
        grid = _grid(nr, nz)
        operator = build_realspace_operator(grid)
        axis_scale = 8.0 / grid.dr**2
        axial_scale = 2.0 / grid.dz**2
        column = 5

        # A = -(L5_r^h + d_zz^h), so row 0 of A is
        #   +8/dr^2 psi_0 - 8/dr^2 psi_1 + (2/dz^2) psi_0 - (1/dz^2)(neighbours)
        on_axis = np.zeros(operator.shape)
        on_axis[0, column] = 1.0
        applied = operator.apply(on_axis)
        assert applied[0, column] == pytest.approx(
            axis_scale + axial_scale, rel=1.0e-14
        )

        first_interior = np.zeros(operator.shape)
        first_interior[1, column] = 1.0
        applied = operator.apply(first_interior)
        assert applied[0, column] == pytest.approx(-axis_scale, rel=1.0e-14)

        # No other radial row couples into the axis row.
        second_interior = np.zeros(operator.shape)
        second_interior[2, column] = 1.0
        assert operator.apply(second_interior)[0, column] == 0.0

        # The coefficient 8 comes from the axis control volume
        # V_0 = (dr/2)^4 / 4 = dr^4 / 64 together with F_0 = (dr/2)^3 / dr,
        # not from a hard-coded constant.
        assert operator.volume[0] == pytest.approx(
            grid.dr**4 / 64.0, rel=1.0e-14
        )
        assert operator.face_coefficient[0] == pytest.approx(
            (0.5 * grid.dr) ** 3 / grid.dr, rel=1.0e-14
        )

        # The same coefficient is reported in the solve metadata.
        _, omega, boundary = _manufactured(grid)
        metadata = solve_realspace_poisson(grid, omega, boundary).metadata
        assert metadata["axis_radial_coefficient_dimensionless"] == (
            pytest.approx(8.0, rel=1.0e-12)
        )


def test_detects_sign_flipped_right_hand_side_against_solver_a() -> None:
    """A relative sign fault must not be absorbed by the C-vs-A comparison.

    Measured discrepancies: 1.0852 and 1.0734 against a solution scale of
    1.3200, i.e. the fault stays at the size of the solution itself instead of
    shrinking by a factor of four under refinement.
    """

    discrepancies: list[float] = []
    solution_scale = 0.0
    for nr, nz in ((17, 32), (33, 64)):
        grid = _grid(nr, nz)
        exact, omega, boundary = _manufactured(grid)
        solution_scale = max(solution_scale, _max_abs(exact))
        psi_a = _solve_a(grid, omega, boundary)
        psi_c = solve_realspace_poisson(grid, -omega, boundary).psi1
        discrepancies.append(_max_abs(psi_a - psi_c))

    for discrepancy in discrepancies:
        assert discrepancy > 0.5 * solution_scale

    coarse, fine = discrepancies
    # The fault does not converge away.  A genuine agreement would shrink by
    # about a factor of four here; measured change is 1.1%.
    assert fine >= 0.95 * coarse
    assert coarse / fine < 1.5


def test_detects_corrupted_dirichlet_trace() -> None:
    """Dropping the outer trace in solver C must be visible against solver A.

    Measured discrepancy: 0.8255 at both resolutions against a solution scale
    of 1.3200.
    """

    discrepancies: list[float] = []
    solution_scale = 0.0
    for nr, nz in ((17, 32), (33, 64)):
        grid = _grid(nr, nz)
        exact, omega, boundary = _manufactured(grid)
        solution_scale = max(solution_scale, _max_abs(exact))
        psi_a = _solve_a(grid, omega, boundary)
        # Injected fault: homogeneous outer data instead of the analytic trace.
        psi_c = solve_realspace_poisson(grid, omega, np.zeros(grid.nz)).psi1
        discrepancies.append(_max_abs(psi_a - psi_c))

    for discrepancy in discrepancies:
        assert discrepancy > 0.5 * solution_scale

    coarse, fine = discrepancies
    assert fine >= 0.95 * coarse
    assert coarse / fine < 1.5

    # A trace perturbed by a smooth O(1) function is also detected, and the
    # perturbation is carried into the interior rather than staying on the
    # boundary row.
    grid = _grid(33, 64)
    _, omega, boundary = _manufactured(grid)
    corrupted = boundary + 0.25 * np.cos(3.0 * grid.z)
    psi_clean = solve_realspace_poisson(grid, omega, boundary).psi1
    psi_dirty = solve_realspace_poisson(grid, omega, corrupted).psi1
    assert _max_abs(psi_dirty[-1] - corrupted) == 0.0
    assert _max_abs(psi_dirty[:-1] - psi_clean[:-1]) > 0.05


def test_detects_loose_tolerance_through_the_manufactured_gate() -> None:
    """A loose CG tolerance must fail the accuracy gate and say so honestly.

    Two independent things are checked, because either alone would be weak.

    1.  Honest reporting: the loose solve must *report* a weighted relative
        residual near its own loose tolerance, not a converged-looking one.
        Measured: 5.457e-03 for tol=1e-2, against 9.763e-13 for tol=1e-12.
    2.  The gate must bite: the loose solution's manufactured error must be
        visibly worse than the discretization error the tight solve achieves.
        Measured: 7.884e-01 loose versus 2.330e-03 tight, a factor of 338,
        while the second-order gate at this resolution is 2 dr^2 = 5.6e-03.

    The unweighted maximum-norm residual is also inspected, because the
    weighted norm gives the axis cell weight dr^4/64 and can therefore look
    small while the axis row is badly unconverged.  Measured for the loose
    solve: 7.09 relative to the right-hand side scale.
    """

    grid = _grid(33, 64)
    exact, omega, boundary = _manufactured(grid)

    tight = solve_realspace_poisson(grid, omega, boundary, tol=1.0e-12)
    loose = solve_realspace_poisson(grid, omega, boundary, tol=1.0e-2)

    tight_error = _max_abs(tight.psi1 - exact)
    loose_error = _max_abs(loose.psi1 - exact)

    # (1) honest reporting.
    assert tight.weighted_relative_residual <= 1.0e-12
    assert 1.0e-4 < loose.weighted_relative_residual <= 1.0e-2
    assert loose.iterations < tight.iterations
    assert tight.metadata["algebraic_residual_max_abs_relative"] <= 1.0e-8
    assert loose.metadata["algebraic_residual_max_abs_relative"] > 1.0e-3

    # (2) the manufactured gate rejects the loose solution and accepts the
    # tight one.
    second_order_gate = 2.0 * grid.dr**2
    assert tight_error <= second_order_gate
    assert loose_error > second_order_gate
    assert loose_error >= 50.0 * tight_error
    assert loose_error > 0.1


def test_rejects_invalid_inputs() -> None:
    """Input guards, written independently of the sibling solvers."""

    grid = _grid(17, 32)
    _, omega, boundary = _manufactured(grid)

    with pytest.raises(TypeError):
        solve_realspace_poisson("not a grid", omega, boundary)
    with pytest.raises(ValueError):
        solve_realspace_poisson(
            AxisymmetricGrid.uniform(nr=17, nz=32, r_max=R_MAX, periodic_z=False),
            np.zeros((17, 32)),
            0.0,
        )
    with pytest.raises(ValueError, match="real-valued"):
        solve_realspace_poisson(grid, omega.astype(np.complex128), boundary)
    with pytest.raises(ValueError, match="real-valued"):
        solve_realspace_poisson(
            grid, np.full(grid.shape, 1 + 2j, dtype=object), boundary
        )
    with pytest.raises(ValueError, match="real-valued"):
        solve_realspace_poisson(grid, omega, boundary.astype(np.complex128))

    broken = omega.copy()
    broken[3, 4] = np.nan
    with pytest.raises(ValueError, match="finite"):
        solve_realspace_poisson(grid, broken, boundary)
    with pytest.raises(ValueError, match="finite"):
        solve_realspace_poisson(grid, omega, boundary * np.inf)

    with pytest.raises(ValueError, match="shape"):
        solve_realspace_poisson(grid, omega[:-1], boundary)
    with pytest.raises(ValueError, match="shape"):
        solve_realspace_poisson(grid, omega, boundary[:-1])

    with pytest.raises(TypeError):
        solve_realspace_poisson(grid, omega, boundary, tol=True)
    with pytest.raises(ValueError):
        solve_realspace_poisson(grid, omega, boundary, tol=0.0)
    with pytest.raises(TypeError):
        solve_realspace_poisson(grid, omega, boundary, max_iterations=1.5)
    with pytest.raises(ValueError):
        solve_realspace_poisson(grid, omega, boundary, max_iterations=0)


def test_iteration_cap_raises_instead_of_returning_an_unconverged_field() -> None:
    """Exceeding the cap must raise, not silently return a partial solve."""

    grid = _grid(33, 64)
    _, omega, boundary = _manufactured(grid)
    with pytest.raises(ConvergenceError, match="did not reach"):
        solve_realspace_poisson(grid, omega, boundary, max_iterations=5)


def test_constant_trace_with_zero_source_returns_that_constant() -> None:
    """A discrete consistency check with an exactly known discrete solution.

    The flux operator annihilates constants in ``r`` and the periodic axial
    difference annihilates constants in ``z``, so ``psi = g`` solves the
    discrete system exactly when ``omega = 0`` and ``g`` is constant.
    """

    grid = _grid(33, 64)
    result = solve_realspace_poisson(grid, np.zeros(grid.shape), 2.5)
    assert _max_abs(result.psi1 - 2.5) <= 1.0e-12


def test_module_uses_no_transform_and_no_sibling_solver() -> None:
    """The independence claim must be enforced, not merely documented.

    Two separate checks: the import graph (AST) and the source text.  The text
    check is what pins the "no transform" claim, because a transform could be
    reached through the already-permitted ``numpy`` import without adding any
    new import statement.
    """

    source = inspect.getsource(realspace_module)
    tree = ast.parse(source)
    imports: set[tuple[int, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update((0, alias.name) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add((node.level, node.module or ""))

    assert imports <= {
        (0, "__future__"),
        (0, "dataclasses"),
        (0, "typing"),
        (0, "numpy"),
        (0, "numpy.typing"),
        (1, "grid"),
    }

    forbidden = {
        "operators",
        "pde",
        "poisson",
        "finite_cylinder_poisson",
        "manufactured",
        "cartesian",
        "cartesian_validation",
        "cartesian_candidate_adapter",
        "nonlinear_cylinder",
        "scipy",
    }
    imported_names = {name for _, name in imports}
    for name in forbidden:
        assert name not in imported_names
        assert f"ns_certificate_lab.{name}" not in imported_names

    # Source-text check: no transform is reachable from this module.  A
    # transform needs no new import statement -- it is reachable through the
    # already-permitted ``numpy`` import -- so the AST check above cannot see
    # it and this text check is the one that pins the claim.
    module_path = Path(realspace_module.__file__).resolve()
    text = module_path.read_text(encoding="utf-8")
    assert "fft" not in text.lower()
    # Positive check that the axial path really is the rolled finite
    # difference this module claims to use.
    assert "np.roll" in text

    # Structural check that the linear solver is the conjugate-gradient
    # iteration and not a tridiagonal elimination.  (The module docstring
    # mentions Thomas elimination when describing the sibling solvers it
    # deliberately does not use, so a raw text scan would be wrong here.)
    function_names = {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert "_preconditioned_conjugate_gradient" in function_names
    for name in function_names:
        lowered = name.lower()
        assert "tridiagonal" not in lowered
        assert "thomas" not in lowered
