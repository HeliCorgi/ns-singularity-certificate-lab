from __future__ import annotations

import ast
import inspect
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.run_poisson_manufactured import (
    evaluate as evaluate_poisson_experiment,
    run as run_poisson_experiment,
)
from ns_certificate_lab._integrity import (
    sha256_file,
    strict_json_loads,
    verify_digest,
)
from ns_certificate_lab.grid import AxisymmetricGrid
import ns_certificate_lab.poisson as poisson_module
from ns_certificate_lab.poisson import (
    independent_physical_poisson_residual,
    solve_streamfunction_poisson,
)


def _manufactured(
    grid: AxisymmetricGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Independent calculus oracle with modes 0, 1 and 2.

    For ``r**(2*q)``, the radial part of L5 is
    ``4*q*(q+1)*r**(2*q-2)``.  This identity is used directly below rather
    than applying any production discrete operator.
    """

    r, z = grid.mesh()
    p1 = 1.0 + 0.2 * r**2 + 0.1 * r**4
    p2 = 0.3 * (1.0 - r**2 + 0.25 * r**4)
    psi = 0.15 * r**2 + p1 * np.cos(z) + p2 * np.sin(2.0 * z)
    radial_l5_p1 = 1.6 + 2.4 * r**2
    radial_l5_p2 = 0.3 * (-8.0 + 6.0 * r**2)
    omega = (
        -1.2
        + (p1 - radial_l5_p1) * np.cos(z)
        + (4.0 * p2 - radial_l5_p2) * np.sin(2.0 * z)
    )
    return psi, omega, psi[-1].copy()


def _test_oracle_residual(
    grid: AxisymmetricGrid,
    psi: np.ndarray,
    omega: np.ndarray,
    expected_boundary: np.ndarray,
) -> np.ndarray:
    """Test-local physical-space residual, independent of the solver module."""

    axial = (
        np.roll(psi, -1, axis=1)
        - 2.0 * psi
        + np.roll(psi, 1, axis=1)
    ) / grid.dz**2
    radial = np.zeros_like(psi)
    radial[0] = 8.0 * (psi[1] - psi[0]) / grid.dr**2
    radial[1:-1] = (
        (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) / grid.dr**2
        + 3.0
        * (psi[2:] - psi[:-2])
        / (2.0 * grid.dr * grid.r[1:-1, None])
    )
    residual = np.empty_like(psi)
    residual[:-1] = -radial[:-1] - axial[:-1] - omega[:-1]
    residual[-1] = psi[-1] - expected_boundary
    return residual


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=np.float64) ** 2)))


def test_multimode_nonzero_boundary_manufactured_solution_converges() -> None:
    spacings: list[float] = []
    solution_errors: list[float] = []
    independent_residuals: list[float] = []
    for nr in (17, 33, 65):
        grid = AxisymmetricGrid.uniform(
            nr=nr,
            nz=2 * (nr - 1),
            r_max=1.0,
            periodic_z=True,
        )
        exact, omega, boundary = _manufactured(grid)
        assert np.max(np.abs(boundary)) > 1.0
        solution = solve_streamfunction_poisson(grid, omega, boundary)
        spacings.append(grid.dr)
        solution_errors.append(_rms(solution.psi1 - exact))
        oracle = _test_oracle_residual(
            grid,
            solution.psi1,
            omega,
            boundary,
        )
        independent_residuals.append(_rms(oracle[:-1]))

        assert np.array_equal(solution.psi1[-1], boundary)
        assert solution.metadata["outer_boundary_max_abs_defect"] == 0.0
        assert (
            solution.metadata["axis_radial_coefficient_dimensionless"]
            == 8.0
        )
        assert (
            solution.metadata["discrete_residual_max_abs_interior"]
            < 5.0e-10
        )
        assert np.allclose(solution.pde_residual, oracle, rtol=0.0, atol=2e-12)
        assert math.isfinite(
            solution.metadata["zero_mode_condition_number_inf"]
        )

    solution_orders = [
        math.log(coarse / fine) / math.log(h_coarse / h_fine)
        for coarse, fine, h_coarse, h_fine in zip(
            solution_errors[:-1],
            solution_errors[1:],
            spacings[:-1],
            spacings[1:],
        )
    ]
    residual_orders = [
        math.log(coarse / fine) / math.log(h_coarse / h_fine)
        for coarse, fine, h_coarse, h_fine in zip(
            independent_residuals[:-1],
            independent_residuals[1:],
            spacings[:-1],
            spacings[1:],
        )
    ]
    assert max(solution_errors) < 5.0e-4
    assert min(solution_orders) > 1.95
    assert min(residual_orders) > 1.95


def test_periodic_fourier_solve_is_translation_covariant() -> None:
    grid = AxisymmetricGrid.uniform(nr=25, nz=48, r_max=1.0)
    r, z = grid.mesh()
    exact = (1.0 + 0.25 * r**2) * np.cos(3.0 * z)
    omega = (
        9.0 * (1.0 + 0.25 * r**2) - 2.0
    ) * np.cos(3.0 * z)
    reference = solve_streamfunction_poisson(grid, omega, exact[-1])
    shift = 11
    shifted = solve_streamfunction_poisson(
        grid,
        np.roll(omega, shift, axis=1),
        np.roll(exact[-1], shift),
    )
    assert np.allclose(
        shifted.psi1,
        np.roll(reference.psi1, shift, axis=1),
        rtol=0.0,
        atol=2.0e-13,
    )
    assert _rms(reference.psi1 - exact) < 2.0e-13


def test_scalar_nonzero_boundary_and_unused_finite_outer_rhs_row() -> None:
    grid = AxisymmetricGrid.uniform(nr=21, nz=32, r_max=1.0)
    r, _ = grid.mesh()
    exact = np.broadcast_to(0.5 + r**2, grid.shape).copy()
    omega = np.full(grid.shape, -8.0)
    reference = solve_streamfunction_poisson(grid, omega, 1.5)
    assert np.max(np.abs(reference.psi1 - exact)) < 3.0e-13
    assert np.array_equal(reference.psi1[-1], np.full(grid.nz, 1.5))

    changed_outer_row = omega.copy()
    changed_outer_row[-1] = np.linspace(-1.0e100, 1.0e100, grid.nz)
    repeated = solve_streamfunction_poisson(
        grid,
        changed_outer_row,
        1.5,
        estimate_condition=False,
    )
    assert np.array_equal(repeated.psi1, reference.psi1)
    assert repeated.metadata["zero_mode_condition_number_inf"] is None
    assert "not evaluated" in repeated.metadata["condition_estimate_interpretation"]


def test_zero_mode_condition_formula_matches_independent_dense_matrix() -> None:
    for nr in (4, 5, 17):
        grid = AxisymmetricGrid.uniform(nr=nr, nz=8, r_max=1.7)
        unknown_r = grid.r[:-1]
        lower_face = np.maximum(unknown_r - 0.5 * grid.dr, 0.0)
        upper_face = unknown_r + 0.5 * grid.dr
        volume = 0.25 * (upper_face**4 - lower_face**4)
        lower = lower_face**3 / (grid.dr * volume)
        upper = upper_face**3 / (grid.dr * volume)
        lower[0] = 0.0
        upper[0] = 8.0 / grid.dr**2
        matrix = np.diag(lower + upper)
        matrix += np.diag(-lower[1:], k=-1)
        matrix += np.diag(-upper[:-1], k=1)
        solution = solve_streamfunction_poisson(
            grid,
            np.zeros(grid.shape),
            np.zeros(grid.nz),
        )
        assert solution.metadata["zero_mode_condition_number_inf"] == (
            pytest.approx(float(np.linalg.cond(matrix, p=np.inf)))
        )


def test_axis_limit_and_overall_sign_are_fixed_by_independent_oracle() -> None:
    grid = AxisymmetricGrid.uniform(nr=33, nz=32, r_max=1.0)
    r, _ = grid.mesh()
    exact = np.broadcast_to(r**2, grid.shape).copy()
    omega = np.full(grid.shape, -8.0)
    good = independent_physical_poisson_residual(
        grid,
        psi1=exact,
        omega1=omega,
        outer_dirichlet=np.ones(grid.nz),
    )
    assert np.max(np.abs(good)) < 3.0e-12

    sign_fault = _test_oracle_residual(
        grid,
        exact,
        -omega,
        np.ones(grid.nz),
    )
    assert _rms(sign_fault[:-1]) > 15.0

    axis_fault = exact.copy()
    axis_fault[0] += 0.01
    axis_defect = _test_oracle_residual(
        grid,
        axis_fault,
        omega,
        np.ones(grid.nz),
    )
    assert np.max(np.abs(axis_defect[0])) > 80.0


def test_boundary_fault_is_not_hidden_in_a_pde_norm() -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=24, r_max=1.0)
    exact, omega, boundary = _manufactured(grid)
    solution = solve_streamfunction_poisson(grid, omega, boundary)
    corrupted = solution.psi1.copy()
    corrupted[-1, 3] += 0.125
    defect = _test_oracle_residual(grid, corrupted, omega, boundary)
    assert np.max(np.abs(defect[-1])) == pytest.approx(0.125)


def test_sign_flipped_rhs_solves_the_wrong_problem_and_is_rejected() -> None:
    grid = AxisymmetricGrid.uniform(nr=25, nz=48, r_max=1.0)
    exact, omega, boundary = _manufactured(grid)
    wrong = solve_streamfunction_poisson(grid, -omega, boundary)
    # A small algebraic residual only confirms the supplied (wrong-sign)
    # discrete equation.  The independent analytic target rejects it.
    assert wrong.metadata["discrete_residual_max_abs_interior"] < 2.0e-10
    assert _rms(wrong.psi1 - exact) > 0.1
    correct_equation_defect = _test_oracle_residual(
        grid,
        wrong.psi1,
        omega,
        boundary,
    )
    assert _rms(correct_equation_defect[:-1]) > 1.0


@pytest.mark.parametrize(
    ("rhs_factory", "boundary_factory", "message"),
    [
        (
            lambda grid: np.zeros((grid.nr - 1, grid.nz)),
            lambda grid: np.zeros(grid.nz),
            "omega1 must have shape",
        ),
        (
            lambda grid: np.full(grid.shape, np.nan),
            lambda grid: np.zeros(grid.nz),
            "omega1 must contain only finite",
        ),
        (
            lambda grid: np.zeros(grid.shape),
            lambda grid: np.zeros(grid.nz - 1),
            "outer_dirichlet must be scalar or have shape",
        ),
        (
            lambda grid: np.zeros(grid.shape),
            lambda grid: np.full(grid.nz, np.inf),
            "outer_dirichlet must contain only finite",
        ),
        (
            lambda grid: np.full(grid.shape, 1.0j),
            lambda grid: np.zeros(grid.nz),
            "omega1 must be real-valued",
        ),
        (
            lambda grid: np.zeros(grid.shape),
            lambda grid: np.full(grid.nz, 1.0j),
            "outer_dirichlet must be real-valued",
        ),
    ],
)
def test_invalid_inputs_are_rejected(
    rhs_factory: object,
    boundary_factory: object,
    message: str,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=24, r_max=1.0)
    rhs = rhs_factory(grid)  # type: ignore[operator]
    boundary = boundary_factory(grid)  # type: ignore[operator]
    with pytest.raises(ValueError, match=message):
        solve_streamfunction_poisson(grid, rhs, boundary)


def test_nonperiodic_z_is_rejected() -> None:
    grid = AxisymmetricGrid.uniform(
        nr=17,
        nz=25,
        r_max=1.0,
        periodic_z=False,
    )
    with pytest.raises(ValueError, match="periodic_z=True"):
        solve_streamfunction_poisson(
            grid,
            np.zeros(grid.shape),
            np.zeros(grid.nz),
        )


def test_manufactured_control_supports_declared_nonstandard_period() -> None:
    repository = Path(__file__).resolve().parents[1]
    config = strict_json_loads(
        (repository / "configs" / "poisson_manufactured.json").read_text(
            encoding="utf-8"
        ),
        label="test Poisson config",
    )
    assert isinstance(config, dict)
    config["z_min"] = -0.75
    config["z_max"] = 3.25
    records, _, orders = evaluate_poisson_experiment(config)
    assert len(records) == 3
    assert min(orders) > 1.9


def test_solver_does_not_import_existing_cylindrical_operator_or_pde_modules() -> None:
    source = inspect.getsource(poisson_module)
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


def test_manufactured_experiment_writes_checksummed_evidence(
    scratch_dir: Path,
) -> None:
    repository = Path(__file__).resolve().parents[1]
    config = strict_json_loads(
        (repository / "configs" / "poisson_manufactured.json").read_text(
            encoding="utf-8"
        ),
        label="test Poisson config",
    )
    assert isinstance(config, dict)
    output = scratch_dir / "poisson-evidence"
    summary = run_poisson_experiment(config, output)
    assert summary["accepted_as_smooth_poisson_control"]
    assert min(summary["observed_orders"]) > 1.95
    assert min(summary["independent_pde_residual_orders"]) > 1.95
    assert summary["checks"]["independent_pde_residual_order_passes"]
    assert summary["checks"]["finest_independent_pde_residual_passes"]

    for name in (
        "config.snapshot.json",
        "diagnostics.csv",
        "diagnostics.json",
        "summary.json",
        "manifest.json",
    ):
        verify_digest(output / name)
    manifest = strict_json_loads(
        (output / "manifest.json").read_text(encoding="utf-8"),
        label="test Poisson manifest",
    )
    assert isinstance(manifest, dict)
    assert set(manifest["files"]) == {
        "config.snapshot.json",
        "diagnostics.csv",
        "diagnostics.json",
        "finest_fields.npz",
        "summary.json",
    }
    for name, entry in manifest["files"].items():
        path = output / name
        assert path.stat().st_size == entry["bytes"]
        assert sha256_file(path) == entry["sha256"]

    with np.load(output / "finest_fields.npz", allow_pickle=False) as archive:
        assert set(archive.files) == {
            "r",
            "z",
            "psi1_numerical",
            "psi1_exact",
            "omega1",
            "outer_dirichlet",
            "discrete_residual",
            "independent_pde_residual",
        }
        assert all(np.all(np.isfinite(archive[name])) for name in archive.files)

    marker = output / "prior-evidence.txt"
    marker.write_text("preserve\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run_poisson_experiment(config, output)
    assert marker.read_text(encoding="utf-8") == "preserve\n"


@pytest.mark.parametrize(
    ("field", "wrong_value"),
    [
        ("experiment_id", "wrong_experiment"),
        ("coordinate_system", "formal five-dimensional space"),
        ("equation", "+(d_rr + 3/r d_r + d_zz) psi1 = omega1"),
    ],
)
def test_manufactured_experiment_rejects_mislabeled_equations(
    scratch_dir: Path,
    field: str,
    wrong_value: str,
) -> None:
    repository = Path(__file__).resolve().parents[1]
    config = strict_json_loads(
        (repository / "configs" / "poisson_manufactured.json").read_text(
            encoding="utf-8"
        ),
        label="test Poisson config",
    )
    assert isinstance(config, dict)
    config[field] = wrong_value
    with pytest.raises(ValueError, match=f"{field} must equal"):
        run_poisson_experiment(config, scratch_dir / f"bad-{field}")


def test_independent_residual_acceptance_gate_cannot_be_omitted(
    scratch_dir: Path,
) -> None:
    repository = Path(__file__).resolve().parents[1]
    config = strict_json_loads(
        (repository / "configs" / "poisson_manufactured.json").read_text(
            encoding="utf-8"
        ),
        label="test Poisson config",
    )
    assert isinstance(config, dict)
    acceptance = config["acceptance"]
    assert isinstance(acceptance, dict)
    acceptance["minimum_independent_residual_order"] = 9.0
    summary = run_poisson_experiment(
        config,
        scratch_dir / "rejected-independent-residual",
    )
    assert not summary["checks"]["independent_pde_residual_order_passes"]
    assert not summary["accepted_as_smooth_poisson_control"]
