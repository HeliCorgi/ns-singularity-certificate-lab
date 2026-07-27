"""Cross-validation of the two independent finite-cylinder Poisson solvers.

Solver A is :mod:`ns_certificate_lab.poisson`
(``solve_streamfunction_poisson``); solver B is
:mod:`ns_certificate_lab.finite_cylinder_poisson`
(``solve_finite_cylinder_poisson``).  Both discretize

    -L5 psi = -(d_rr + 3/r d_r + d_zz) psi = omega

on a finite cylinder with periodic ``z`` and explicit outer Dirichlet data, but
with different radial discretizations: A uses an ``r^3``-flux finite volume, B
uses the direct non-divergence stencil.  This is the ONLY test module that is
permitted to import both solver modules; every other test keeps one solver
isolated from the other.

Every manufactured field below is a closed-form analytic expression written out
in this file.  No production code is used to derive an oracle, and every
convergence order is recomputed here from the measured errors.
"""

from __future__ import annotations

import ast
import inspect
import math
from pathlib import Path

import numpy as np

import ns_certificate_lab.finite_cylinder_poisson as finite_cylinder_module
from ns_certificate_lab.finite_cylinder_poisson import solve_finite_cylinder_poisson
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.poisson import solve_streamfunction_poisson

R_MAX = 1.7
Z_PERIOD = 2.0 * np.pi


def _grid(nr: int, nz: int) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=R_MAX,
        z_min=0.0,
        z_max=Z_PERIOD,
        periodic_z=True,
    )


def _cv1_fields(
    grid: AxisymmetricGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``psi*``, ``omega* = -L5 psi*`` and the outer trace ``g(z)``.

    The analytic identity used here is ``L5(r**(2*k)) = 4*k*(k+1)*r**(2*k-2)``,
    i.e. ``L5(r**2) = 8`` and ``L5(r**4) = 24*r**2``; it is verified directly:
    ``d_rr r**4 + (3/r) d_r r**4 = 12 r**2 + 12 r**2 = 24 r**2``.  With

        psi*(r,z) = (1.1 + a r^2 + b r^4) + (c + d r^2 + e r^4) cos(q z)

    the radial part of ``L5 psi*`` is ``(8a + 24 b r^2)`` plus
    ``(8d + 24 e r^2) cos(q z)``, and ``d_zz`` contributes
    ``-q^2 (c + d r^2 + e r^4) cos(q z)``.  Hence

        omega* = -L5 psi*
               = -(8a + 24 b r^2)
                 + (q^2 (c + d r^2 + e r^4) - 8d - 24 e r^2) cos(q z).
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


def _cv2_fields(
    grid: AxisymmetricGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a field whose radial profiles are exact for both stencils.

    Only ``a + b r^2`` radial profiles appear, and both radial discretizations
    reproduce ``L5(a + b r^2) = 8 b`` exactly, while the axial dependence is a
    pair of grid-exact Fourier modes.  Both solvers must therefore return
    ``psi*`` to roundoff.
    """

    r, z = grid.mesh()
    period = float(grid.z_period)
    q1 = 2.0 * np.pi * 3.0 / period
    q2 = 2.0 * np.pi * 5.0 / period

    radial_mode = 0.9 - 0.35 * r**2
    radial_constant = 0.4 + 0.21 * r**2
    axial = np.cos(q1 * z) + 0.4 * np.sin(q2 * z)
    # d_zz of the axial factor, differentiated analytically.
    axial_zz = -q1 * q1 * np.cos(q1 * z) - 0.4 * q2 * q2 * np.sin(q2 * z)

    psi = radial_mode * axial + radial_constant
    # L5 psi* = [radial L5 of the r-profiles] + [d_zz of the axial factor]:
    #   radial L5 (0.9 - 0.35 r^2) = 8 * (-0.35),
    #   radial L5 (0.4 + 0.21 r^2) = 8 * 0.21.
    l5_psi = 8.0 * (-0.35) * axial + 8.0 * 0.21 + radial_mode * axial_zz
    omega = -l5_psi
    return psi, omega, psi[-1, :].copy()


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


def _solve_b(
    grid: AxisymmetricGrid,
    omega: np.ndarray,
    boundary: np.ndarray,
) -> np.ndarray:
    return solve_finite_cylinder_poisson(
        omega,
        grid,
        outer_boundary=boundary,
        condition_mode_indices=(),
    ).psi


def _order(coarse_error: float, fine_error: float, coarse_h: float, fine_h: float) -> float:
    return math.log(coarse_error / fine_error) / math.log(coarse_h / fine_h)


def _max_abs(values: np.ndarray) -> float:
    return float(np.max(np.abs(values)))


def test_independent_solvers_agree_at_second_order_without_collapsing() -> None:
    """CV-1: the two solvers must agree at ``O(dr^2)`` and only at ``O(dr^2)``."""

    spacings: list[float] = []
    differences: list[float] = []
    errors_a: list[float] = []
    errors_b: list[float] = []

    for nr in (17, 33, 65):
        grid = _grid(nr, 2 * (nr - 1))
        exact, omega, boundary = _cv1_fields(grid)
        assert _max_abs(boundary) > 0.0
        # The trace is genuinely z-dependent: its peak-to-peak amplitude is
        # 2*|0.22 + (-0.19)*R^2 + 0.035*R^4| = 0.0736 at R = 1.7.
        assert float(np.max(boundary) - np.min(boundary)) > 0.05

        # Identical right-hand side, identical Dirichlet trace, identical grid.
        psi_a = _solve_a(grid, omega, boundary)
        psi_b = _solve_b(grid, omega, boundary)

        spacings.append(grid.dr)
        differences.append(_max_abs(psi_a - psi_b))
        errors_a.append(_max_abs(psi_a - exact))
        errors_b.append(_max_abs(psi_b - exact))

    # The bounds are checked before any order is computed so that a collapsed
    # difference reports the explanatory assertion below instead of a division
    # by zero inside the order formula.
    for spacing, difference in zip(spacings, differences):
        # Upper bound: the two discretizations differ only by their truncation
        # error, measured at D(h) ~= 0.115 * dr^2 on this machine.
        assert difference <= 1.0 * spacing**2
        # CRITICAL LOWER BOUND.  The two solvers are meant to be genuinely
        # independent discretizations, so their difference must remain a real
        # O(dr^2) truncation gap.  If a future refactor ever makes one solver
        # call the other's stencil (or makes them share an assembly routine),
        # D(h) would collapse to roundoff and this clause would fail.  A
        # passing agreement test would then be vacuous.  Never weaken or
        # remove this bound to make a refactor green: an intentional merge of
        # the two implementations must be reviewed as the loss of independent
        # cross-validation that it is.
        assert difference >= 1.0e-3 * spacing**2

    difference_orders = [
        _order(differences[i], differences[i + 1], spacings[i], spacings[i + 1])
        for i in range(len(differences) - 1)
    ]
    assert min(difference_orders) >= 1.80
    assert max(difference_orders) <= 2.20

    # Each solver is independently second-order accurate against the analytic
    # field, so the agreement above is not two implementations sharing a bug.
    for errors in (errors_a, errors_b):
        orders = [
            _order(errors[i], errors[i + 1], spacings[i], spacings[i + 1])
            for i in range(len(errors) - 1)
        ]
        assert min(orders) >= 1.85
        assert errors[-1] <= 2.0 * spacings[-1] ** 2


def test_radially_exact_field_pins_shared_solver_conventions() -> None:
    """CV-2: on a field both stencils integrate exactly, they must agree to roundoff.

    Any disagreement in FFT normalization, wavenumber convention, Nyquist
    handling, Dirichlet semantics, ``(n_r, n_z)`` index order, operator sign or
    axis coefficient would show up here far above ``1e-12``.
    """

    for nr, nz in ((17, 32), (33, 64)):
        grid = _grid(nr, nz)
        exact, omega, boundary = _cv2_fields(grid)
        psi_a = _solve_a(grid, omega, boundary)
        psi_b = _solve_b(grid, omega, boundary)

        assert _max_abs(psi_a - psi_b) <= 1.0e-12
        assert _max_abs(psi_a - exact) <= 1.0e-12
        assert _max_abs(psi_b - exact) <= 1.0e-12


def test_detects_sign_flipped_right_hand_side_between_solvers() -> None:
    """CV-3(i): a relative sign fault must not be absorbed by the comparison."""

    discrepancies: list[float] = []
    solution_scale = 0.0
    for nr in (17, 33):
        grid = _grid(nr, 2 * (nr - 1))
        exact, omega, boundary = _cv1_fields(grid)
        solution_scale = max(solution_scale, _max_abs(exact))
        psi_a = _solve_a(grid, omega, boundary)
        psi_b = _solve_b(grid, -omega, boundary)  # injected relative sign fault
        discrepancies.append(_max_abs(psi_a - psi_b))

    for discrepancy in discrepancies:
        assert discrepancy > 0.5 * solution_scale

    coarse, fine = discrepancies
    # The fault does not converge away.  A genuine agreement would shrink by
    # about a factor of four under this refinement; here the discrepancy stays
    # at the size of the solution itself.  (Measured on this machine: 1.0800
    # then 1.0722, a 0.7% change, so the honest pin is "does not decrease at a
    # convergent rate" rather than a literal non-decrease.)
    assert fine >= 0.95 * coarse
    assert coarse / fine < 1.5


def test_detects_dropped_dirichlet_trace_between_solvers() -> None:
    """CV-3(ii): dropping the outer trace in one solver must be detected."""

    discrepancies: list[float] = []
    solution_scale = 0.0
    for nr in (17, 33):
        grid = _grid(nr, 2 * (nr - 1))
        exact, omega, boundary = _cv1_fields(grid)
        solution_scale = max(solution_scale, _max_abs(exact))
        psi_a = _solve_a(grid, omega, boundary)
        # Injected fault: solver B is given homogeneous outer data instead of
        # the analytic trace g(z).
        psi_b = _solve_b(grid, omega, np.zeros(grid.nz))
        discrepancies.append(_max_abs(psi_a - psi_b))

    for discrepancy in discrepancies:
        assert discrepancy > 0.5 * solution_scale

    coarse, fine = discrepancies
    assert fine >= 0.95 * coarse
    assert coarse / fine < 1.5


def test_detects_axis_row_perturbation_in_the_comparison_metric() -> None:
    """CV-3(iii): the comparison metric must see a fault confined to ``r=0``."""

    grid = _grid(17, 32)
    _, omega, boundary = _cv1_fields(grid)
    psi_a = _solve_a(grid, omega, boundary)
    psi_b = _solve_b(grid, omega, boundary)

    baseline = _max_abs(psi_a - psi_b)
    assert baseline < 1.0e-2  # the honest agreement is far below the fault

    corrupted = psi_a.copy()
    corrupted[0, :] += 0.01  # injected axis-row fault, post-solve
    assert _max_abs(corrupted - psi_b) > 1.0e-2


def test_finite_cylinder_solver_does_not_import_operator_pde_or_sibling_modules() -> None:
    """Mirror of the solver-A import guard, applied to solver B.

    Independence is only meaningful if solver B never reaches into the shared
    cylindrical operators, the PDE driver, the sibling Poisson solver, the
    manufactured-solution helpers or the Cartesian modules.
    """

    source = inspect.getsource(finite_cylinder_module)
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
        (0, "functools"),
        (0, "math"),
        (0, "typing"),
        (0, "numpy"),
        (0, "numpy.typing"),
        (1, "grid"),
    }

    forbidden = {
        "operators",
        "pde",
        "poisson",
        "manufactured",
        "cartesian",
        "cartesian_validation",
        "cartesian_candidate_adapter",
    }
    imported_names = {name for _, name in imports}
    for name in forbidden:
        assert name not in imported_names
        assert f"ns_certificate_lab.{name}" not in imported_names


def test_cross_validation_is_the_only_module_importing_both_solvers() -> None:
    """The isolation claim above is only credible if it is enforced somewhere."""

    this_file = Path(__file__).resolve()
    for path in sorted(this_file.parent.glob("test_*.py")):
        if path == this_file:
            continue
        text = path.read_text(encoding="utf-8")
        imports_a = "ns_certificate_lab.poisson" in text
        imports_b = "ns_certificate_lab.finite_cylinder_poisson" in text
        assert not (imports_a and imports_b), path.name
