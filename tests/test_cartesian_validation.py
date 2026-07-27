"""Tests for the independent uniform-Cartesian primitive audit path."""

from __future__ import annotations

import ast
import inspect

import numpy as np
import pytest

import ns_certificate_lab.cartesian_candidate_adapter as candidate_adapter_module
import ns_certificate_lab.cartesian_validation as cartesian_module
from ns_certificate_lab.artifacts import (
    CandidateDescription,
    load_candidate,
    save_candidate,
)
from ns_certificate_lab.cartesian_candidate_adapter import (
    reconstruct_loaded_candidate_on_cartesian,
)
from ns_certificate_lab.cartesian_validation import (
    RegularizedAxisymmetricProfile,
    UniformCartesianGrid,
    audit_cartesian_reconstruction,
    cartesian_curl,
    cartesian_divergence,
    cartesian_gradient,
    cartesian_vector_laplacian,
    primitive_ns_residual,
    reconstruct_axisymmetric_regularized,
)
from ns_certificate_lab.grid import AxisymmetricGrid

_ANALYTIC_AUDIT_TOLERANCES = {
    "divergence_rms_tolerance": 0.005,
    "divergence_max_tolerance": 0.01,
    "curl_rms_tolerance": 0.006,
    "curl_max_tolerance": 0.015,
}

_SAVED_CANDIDATE_TOLERANCES = {
    "divergence_rms_tolerance": 0.005,
    "divergence_max_tolerance": 0.01,
    "curl_rms_tolerance": 0.015,
    "curl_max_tolerance": 0.04,
}


def _candidate_description() -> CandidateDescription:
    return CandidateDescription(
        representation="nodal little-endian float64 arrays",
        coordinate_system="axisymmetric cylindrical half-plane (r,z)",
        units={
            "r": "dimensionless length",
            "z": "dimensionless length",
            "u1": "dimensionless transformed swirl",
            "omega1": "dimensionless transformed azimuthal vorticity",
            "psi1": "dimensionless transformed streamfunction",
        },
        normalization="dimensionless manufactured normalization",
        physical_time=0.2,
        viscosity=0.05,
        basis_convention=(
            "array axes (r,z); E-18a velocity; E-18b full vorticity"
        ),
    )


def _manufactured_nodal_arrays(
    grid: AxisymmetricGrid,
    *,
    omega1_sign: float = 1.0,
) -> dict[str, np.ndarray]:
    """Direct transcription of the local manufactured closed forms."""

    r, z = grid.mesh()
    r_squared = r * r
    q = (1.0 - r_squared) ** 2
    a = np.exp(-0.2)
    b = np.exp(-0.4)
    return {
        "psi1": a * q * np.cos(z),
        "omega1": (
            omega1_sign
            * a
            * (17.0 - 26.0 * r_squared + r_squared * r_squared)
            * np.cos(z)
        ),
        "u1": b * q * np.sin(2.0 * z),
    }


def _save_manufactured_candidate(
    archive,
    *,
    nr: int,
    nz: int,
    omega1_sign: float = 1.0,
):
    grid = AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=1.0,
        periodic_z=True,
    )
    save_candidate(
        archive,
        grid=grid,
        fields=_manufactured_nodal_arrays(grid, omega1_sign=omega1_sign),
        config={"test": "saved Cartesian end-to-end audit"},
        seed=2718,
        description=_candidate_description(),
    )
    return load_candidate(archive)


def _manufactured_cartesian_oracle(
    grid: UniformCartesianGrid,
) -> tuple[np.ndarray, np.ndarray]:
    """Direct E-18a/E-18b closed forms, independent of both adapter paths."""

    x, y, z = grid.mesh()
    radius_squared = x * x + y * y
    radius = np.sqrt(radius_squared)
    q = (1.0 - radius_squared) ** 2
    q_r = -4.0 * radius + 4.0 * radius * radius_squared
    a = np.exp(-0.2)
    b = np.exp(-0.4)
    u1 = b * q * np.sin(2.0 * z)
    u1_r = b * q_r * np.sin(2.0 * z)
    u1_z = 2.0 * b * q * np.cos(2.0 * z)
    psi1 = a * q * np.cos(z)
    psi1_r = a * q_r * np.cos(z)
    psi1_z = -a * q * np.sin(z)
    omega1 = (
        a
        * (17.0 - 26.0 * radius_squared + radius_squared * radius_squared)
        * np.cos(z)
    )
    velocity = np.stack(
        (
            -x * psi1_z - y * u1,
            -y * psi1_z + x * u1,
            2.0 * psi1 + radius * psi1_r,
        ),
        axis=0,
    )
    vorticity = np.stack(
        (
            -x * u1_z - y * omega1,
            -y * u1_z + x * omega1,
            2.0 * u1 + radius * u1_r,
        ),
        axis=0,
    )
    return velocity, vorticity


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=np.float64) ** 2)))


def _orders(errors: list[float], spacings: list[float]) -> list[float]:
    return [
        float(np.log(errors[index] / errors[index + 1])
              / np.log(spacings[index] / spacings[index + 1]))
        for index in range(len(errors) - 1)
    ]


def _periodic_manufactured(
    grid: UniformCartesianGrid,
    *,
    time: float = 0.2,
    viscosity: float = 0.07,
) -> dict[str, np.ndarray]:
    """Analytic divergence-free field; no numerical derivative is called."""

    x, y, z = grid.mesh()
    amplitude = np.exp(-time)
    sin_x, cos_x = np.sin(x), np.cos(x)
    sin_2x, cos_2x = np.sin(2.0 * x), np.cos(2.0 * x)
    sin_y, cos_y = np.sin(y), np.cos(y)
    sin_z, cos_z = np.sin(z), np.cos(z)

    velocity = np.stack(
        (
            amplitude * sin_2x * cos_y * cos_z,
            -2.0 * amplitude * cos_2x * sin_y * cos_z,
            np.zeros_like(x),
        ),
        axis=0,
    )
    velocity_t = -velocity
    pressure = sin_x * sin_y * sin_z
    pressure_gradient = np.stack(
        (
            cos_x * sin_y * sin_z,
            sin_x * cos_y * sin_z,
            sin_x * sin_y * cos_z,
        ),
        axis=0,
    )
    vorticity = np.stack(
        (
            -2.0 * amplitude * cos_2x * sin_y * sin_z,
            -amplitude * sin_2x * cos_y * sin_z,
            5.0 * amplitude * sin_2x * sin_y * cos_z,
        ),
        axis=0,
    )
    laplacian = -6.0 * velocity
    advection = np.stack(
        (
            2.0 * amplitude**2 * sin_2x * cos_2x * cos_z**2,
            4.0 * amplitude**2 * sin_y * cos_y * cos_z**2,
            np.zeros_like(x),
        ),
        axis=0,
    )
    body_force = velocity_t + advection + pressure_gradient - viscosity * laplacian
    return {
        "velocity": velocity,
        "velocity_t": velocity_t,
        "pressure": pressure,
        "vorticity": vorticity,
        "laplacian": laplacian,
        "advection": advection,
        "pressure_gradient": pressure_gradient,
        "body_force": body_force,
    }


def _periodic_grid(count: int) -> UniformCartesianGrid:
    return UniformCartesianGrid.uniform(
        shape=(count, count, count),
        x_bounds=(-np.pi, np.pi),
        y_bounds=(-np.pi, np.pi),
        z_bounds=(-np.pi, np.pi),
        periodic=(True, True, True),
    )


def test_cartesian_operators_and_primitive_residual_converge() -> None:
    metrics = {
        "divergence": [],
        "curl": [],
        "laplacian": [],
        "advection": [],
        "pressure_gradient": [],
        "viscous_term": [],
        "primitive_residual": [],
    }
    spacings: list[float] = []
    for count in (12, 24, 48):
        grid = _periodic_grid(count)
        exact = _periodic_manufactured(grid)
        residual = primitive_ns_residual(
            grid,
            velocity=exact["velocity"],
            velocity_t=exact["velocity_t"],
            pressure=exact["pressure"],
            viscosity=0.07,
            body_force=exact["body_force"],
        )
        metrics["divergence"].append(_rms(residual.divergence))
        metrics["curl"].append(
            _rms(cartesian_curl(grid, exact["velocity"]) - exact["vorticity"])
        )
        metrics["laplacian"].append(
            _rms(
                cartesian_vector_laplacian(grid, exact["velocity"])
                - exact["laplacian"]
            )
        )
        metrics["advection"].append(
            _rms(residual.advection - exact["advection"])
        )
        metrics["pressure_gradient"].append(
            _rms(residual.pressure_gradient - exact["pressure_gradient"])
        )
        metrics["viscous_term"].append(
            _rms(residual.viscous + 0.07 * exact["laplacian"])
        )
        metrics["primitive_residual"].append(
            _rms(residual.defect_against_body_force)
        )
        spacings.append(max(grid.spacings))

        assert np.allclose(
            residual.total,
            residual.time_derivative
            + residual.advection
            + residual.pressure_gradient
            + residual.viscous,
        )
        assert np.array_equal(residual.time_derivative, exact["velocity_t"])
        assert np.array_equal(residual.body_force, exact["body_force"])
        assert np.allclose(
            residual.defect_against_body_force,
            residual.total - residual.body_force,
        )

    for name, errors in metrics.items():
        assert errors[-1] < errors[0], (name, errors)
        assert min(_orders(errors, spacings)) > 1.8, (name, errors)


def test_nonperiodic_boundary_closures_are_exact_for_quadratics() -> None:
    """Exercise every one-sided first/second derivative boundary branch."""

    grid = UniformCartesianGrid.uniform(
        shape=(7, 8, 9),
        x_bounds=(-0.7, 1.1),
        y_bounds=(-1.3, 0.9),
        z_bounds=(-0.4, 1.5),
        periodic=(False, False, False),
    )
    x, y, z = grid.mesh()
    scalar = x * x + 2.0 * y * y - 0.5 * z * z + x * y
    expected_gradient = np.stack(
        (2.0 * x + y, 4.0 * y + x, -z),
        axis=0,
    )
    velocity = np.stack(
        (
            x * x + 2.0 * x * y + 3.0 * z,
            -y * y + x * z,
            z * z + y * z - 4.0 * x,
        ),
        axis=0,
    )
    expected_divergence = 2.0 * x + y + 2.0 * z
    expected_curl = np.stack(
        (z - x, np.full_like(x, 7.0), z - 2.0 * x),
        axis=0,
    )
    expected_laplacian = np.stack(
        (
            np.full_like(x, 2.0),
            np.full_like(x, -2.0),
            np.full_like(x, 2.0),
        ),
        axis=0,
    )

    assert np.max(
        np.abs(cartesian_gradient(grid, scalar) - expected_gradient)
    ) < 2.0e-13
    assert np.max(
        np.abs(cartesian_divergence(grid, velocity) - expected_divergence)
    ) < 2.0e-13
    assert np.max(
        np.abs(cartesian_curl(grid, velocity) - expected_curl)
    ) < 2.0e-13
    assert np.max(
        np.abs(
            cartesian_vector_laplacian(grid, velocity)
            - expected_laplacian
        )
    ) < 2.0e-12


def _axisymmetric_candidate(
    grid: UniformCartesianGrid,
) -> tuple[np.ndarray, np.ndarray]:
    """Smooth regularized swirl profile and its independent exact curl."""

    profile = RegularizedAxisymmetricProfile(
        u1=lambda r, z: np.exp(-(r * r)) * np.sin(z),
        psi1=lambda r, z: np.exp(-(r * r)) * np.cos(z),
        dpsi1_dr=lambda r, z: -2.0 * r * np.exp(-(r * r)) * np.cos(z),
        dpsi1_dz=lambda r, z: -np.exp(-(r * r)) * np.sin(z),
    )
    velocity = reconstruct_axisymmetric_regularized(grid, profile)

    # Analytic full curl, derived directly in Cartesian coordinates:
    # omega_x=-x u1_z-y omega1, omega_y=-y u1_z+x omega1,
    # omega_z=2u1+r u1_r, with omega1=-L5 psi1.
    x, y, z = grid.mesh()
    radius_squared = x * x + y * y
    exponential = np.exp(-radius_squared)
    u1_z = exponential * np.cos(z)
    omega1 = (9.0 - 4.0 * radius_squared) * exponential * np.cos(z)
    omega_z = 2.0 * (1.0 - radius_squared) * exponential * np.sin(z)
    vorticity = np.stack(
        (
            -x * u1_z - y * omega1,
            -y * u1_z + x * omega1,
            omega_z,
        ),
        axis=0,
    )
    return velocity, vorticity


def _axisymmetric_grid() -> UniformCartesianGrid:
    return UniformCartesianGrid.uniform(
        shape=(33, 33, 32),
        x_bounds=(-0.8, 0.8),
        y_bounds=(-0.8, 0.8),
        z_bounds=(0.0, 2.0 * np.pi),
        periodic=(False, False, True),
    )


def test_axisymmetric_candidate_reconstruction_passes_end_to_end_audit() -> None:
    grid = _axisymmetric_grid()
    velocity, vorticity = _axisymmetric_candidate(grid)
    report = audit_cartesian_reconstruction(
        grid,
        velocity=velocity,
        expected_vorticity=vorticity,
        **_ANALYTIC_AUDIT_TOLERANCES,
    )
    assert report.passed, report
    assert report.divergence_max < 0.01
    assert report.curl_defect_max < 0.015


def test_saved_candidate_arrays_pass_uniform_cartesian_end_to_end_audit(
    scratch_dir,
) -> None:
    archive = scratch_dir / "cartesian-audit-candidate.npz"
    loaded = _save_manufactured_candidate(
        archive,
        nr=65,
        nz=128,
    )
    cartesian_grid = UniformCartesianGrid.uniform(
        shape=(25, 25, 64),
        x_bounds=(-0.6, 0.6),
        y_bounds=(-0.6, 0.6),
        z_bounds=(0.0, 2.0 * np.pi),
        periodic=(False, False, True),
    )
    reconstructed = reconstruct_loaded_candidate_on_cartesian(
        loaded,
        cartesian_grid,
    )
    report = audit_cartesian_reconstruction(
        cartesian_grid,
        velocity=reconstructed.velocity,
        expected_vorticity=reconstructed.vorticity_e18b,
        **_SAVED_CANDIDATE_TOLERANCES,
        interior_margin=2,
    )
    assert report.all_points_finite
    assert report.passed, report


def test_saved_candidate_adapter_converges_to_direct_cartesian_closed_form(
    scratch_dir,
) -> None:
    cartesian_grid = UniformCartesianGrid.uniform(
        shape=(17, 19, 48),
        x_bounds=(-0.55, 0.55),
        y_bounds=(-0.55, 0.55),
        z_bounds=(0.0, 2.0 * np.pi),
        periodic=(False, False, True),
    )
    exact_velocity, exact_vorticity = _manufactured_cartesian_oracle(
        cartesian_grid
    )
    velocity_errors: list[float] = []
    vorticity_errors: list[float] = []
    spacings: list[float] = []
    for nr, nz in ((33, 64), (65, 128)):
        loaded = _save_manufactured_candidate(
            scratch_dir / f"candidate-{nr}-{nz}.npz",
            nr=nr,
            nz=nz,
        )
        reconstructed = reconstruct_loaded_candidate_on_cartesian(
            loaded,
            cartesian_grid,
        )
        velocity_errors.append(
            _rms(reconstructed.velocity - exact_velocity)
        )
        vorticity_errors.append(
            _rms(reconstructed.vorticity_e18b - exact_vorticity)
        )
        spacings.append(max(loaded.grid.dr, loaded.grid.dz))

    assert velocity_errors[1] < velocity_errors[0]
    assert vorticity_errors[1] < vorticity_errors[0]
    assert velocity_errors[1] < 2.0e-4
    assert vorticity_errors[1] < 7.0e-4
    assert _orders(velocity_errors, spacings)[0] > 1.7
    assert _orders(vorticity_errors, spacings)[0] > 1.7


def test_saved_candidate_omega1_sign_fault_is_rejected_after_reload(
    scratch_dir,
) -> None:
    loaded = _save_manufactured_candidate(
        scratch_dir / "wrong-omega-sign.npz",
        nr=65,
        nz=128,
        omega1_sign=-1.0,
    )
    cartesian_grid = UniformCartesianGrid.uniform(
        shape=(25, 25, 64),
        x_bounds=(-0.6, 0.6),
        y_bounds=(-0.6, 0.6),
        z_bounds=(0.0, 2.0 * np.pi),
        periodic=(False, False, True),
    )
    reconstructed = reconstruct_loaded_candidate_on_cartesian(
        loaded,
        cartesian_grid,
    )
    report = audit_cartesian_reconstruction(
        cartesian_grid,
        velocity=reconstructed.velocity,
        expected_vorticity=reconstructed.vorticity_e18b,
        **_SAVED_CANDIDATE_TOLERANCES,
        interior_margin=2,
    )
    assert report.divergence_passed
    assert not report.curl_passed
    assert not report.passed


def test_independent_audit_rejects_cylindrical_radial_sign_fault() -> None:
    grid = _axisymmetric_grid()
    velocity, vorticity = _axisymmetric_candidate(grid)
    x, y, z = grid.mesh()
    psi_z = -np.exp(-(x * x + y * y)) * np.sin(z)
    u1 = np.exp(-(x * x + y * y)) * np.sin(z)
    wrong = velocity.copy()
    # Wrong u^r sign while preserving the swirl component.
    wrong[0] = x * psi_z - y * u1
    wrong[1] = y * psi_z + x * u1
    report = audit_cartesian_reconstruction(
        grid,
        velocity=wrong,
        expected_vorticity=vorticity,
        **_ANALYTIC_AUDIT_TOLERANCES,
    )
    assert not report.passed
    assert not report.divergence_passed


def test_independent_audit_rejects_component_transform_fault() -> None:
    grid = _axisymmetric_grid()
    velocity, vorticity = _axisymmetric_candidate(grid)
    wrong = velocity.copy()
    wrong[[0, 1]] = wrong[[1, 0]]
    report = audit_cartesian_reconstruction(
        grid,
        velocity=wrong,
        expected_vorticity=vorticity,
        **_ANALYTIC_AUDIT_TOLERANCES,
    )
    assert not report.passed
    assert not report.curl_passed


def test_independent_audit_rejects_vorticity_sign_fault() -> None:
    grid = _axisymmetric_grid()
    velocity, vorticity = _axisymmetric_candidate(grid)
    report = audit_cartesian_reconstruction(
        grid,
        velocity=velocity,
        expected_vorticity=-vorticity,
        **_ANALYTIC_AUDIT_TOLERANCES,
    )
    assert not report.passed
    assert report.divergence_passed
    assert not report.curl_passed


def test_independent_audit_rejects_divergence_pollution() -> None:
    grid = _axisymmetric_grid()
    velocity, vorticity = _axisymmetric_candidate(grid)
    x, _, _ = grid.mesh()
    wrong = velocity.copy()
    wrong[0] += 0.2 * x
    report = audit_cartesian_reconstruction(
        grid,
        velocity=wrong,
        expected_vorticity=vorticity,
        **_ANALYTIC_AUDIT_TOLERANCES,
    )
    assert not report.passed
    assert not report.divergence_passed


def test_localized_curl_defect_is_rejected_by_maximum_gate() -> None:
    grid = _axisymmetric_grid()
    velocity, vorticity = _axisymmetric_candidate(grid)
    localized_fault = vorticity.copy()
    localized_fault[0, 16, 16, 10] += 0.2
    report = audit_cartesian_reconstruction(
        grid,
        velocity=velocity,
        expected_vorticity=localized_fault,
        **_ANALYTIC_AUDIT_TOLERANCES,
    )
    assert report.curl_rms_passed
    assert not report.curl_max_passed
    assert not report.passed


def test_periodic_z_seam_fault_is_not_hidden_by_interior_margin() -> None:
    grid = _axisymmetric_grid()
    velocity, vorticity = _axisymmetric_candidate(grid)
    seam_fault = velocity.copy()
    seam_fault[0, :, :, 0] += 0.3
    report = audit_cartesian_reconstruction(
        grid,
        velocity=seam_fault,
        expected_vorticity=vorticity,
        **_ANALYTIC_AUDIT_TOLERANCES,
        interior_margin=2,
    )
    assert not report.curl_rms_passed
    assert not report.passed


def _imports_in(module) -> set[tuple[int, str]]:
    tree = ast.parse(inspect.getsource(module))
    imports: set[tuple[int, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update((0, alias.name) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add((node.level, node.module or ""))
    return imports


def test_cartesian_validation_has_no_cylindrical_operator_dependency() -> None:
    checker_imports = _imports_in(cartesian_module)
    adapter_imports = _imports_in(candidate_adapter_module)
    assert checker_imports <= {
        (0, "__future__"),
        (0, "dataclasses"),
        (0, "typing"),
        (0, "numpy"),
        (0, "numpy.typing"),
    }
    assert adapter_imports <= {
        (0, "__future__"),
        (0, "dataclasses"),
        (0, "numpy"),
        (0, "numpy.typing"),
        (1, "artifacts"),
        (1, "cartesian_validation"),
    }


@pytest.mark.parametrize(
    "tolerance_name",
    (
        "divergence_rms_tolerance",
        "divergence_max_tolerance",
        "curl_rms_tolerance",
        "curl_max_tolerance",
    ),
)
def test_cartesian_audit_requires_each_positive_tolerance(
    tolerance_name: str,
) -> None:
    grid = _periodic_grid(4)
    tolerances = {
        "divergence_rms_tolerance": 1.0,
        "divergence_max_tolerance": 1.0,
        "curl_rms_tolerance": 1.0,
        "curl_max_tolerance": 1.0,
    }
    tolerances[tolerance_name] = 0.0
    with pytest.raises(ValueError, match=tolerance_name):
        audit_cartesian_reconstruction(
            grid,
            velocity=np.zeros((3, *grid.shape)),
            expected_vorticity=np.zeros((3, *grid.shape)),
            **tolerances,
        )


@pytest.mark.parametrize(
    ("field_name", "values"),
    (
        ("velocity", np.zeros((3, 4, 4, 3))),
        ("velocity_t", np.zeros((3, 4, 4, 3))),
        ("pressure", np.zeros((4, 4, 3))),
        ("body_force", np.zeros((3, 4, 4, 3))),
    ),
)
def test_primitive_residual_rejects_wrong_shapes(
    field_name: str,
    values: np.ndarray,
) -> None:
    grid = UniformCartesianGrid.uniform(
        shape=(4, 4, 4),
        x_bounds=(-1.0, 1.0),
        y_bounds=(-1.0, 1.0),
        z_bounds=(-1.0, 1.0),
    )
    arguments: dict[str, object] = {
        "velocity": np.zeros((3, 4, 4, 4)),
        "velocity_t": np.zeros((3, 4, 4, 4)),
        "pressure": np.zeros((4, 4, 4)),
        "viscosity": 0.1,
        "body_force": np.zeros((3, 4, 4, 4)),
    }
    arguments[field_name] = values
    with pytest.raises(ValueError, match="shape"):
        primitive_ns_residual(grid, **arguments)


def test_cartesian_validation_rejects_nonfinite_inputs() -> None:
    grid = _periodic_grid(4)
    velocity = np.zeros((3, 4, 4, 4))
    velocity[1, 2, 1, 3] = np.nan
    with pytest.raises(ValueError, match="finite"):
        cartesian_divergence(grid, velocity)
    with pytest.raises(ValueError, match="nonnegative"):
        primitive_ns_residual(
            grid,
            velocity=np.zeros((3, 4, 4, 4)),
            velocity_t=np.zeros((3, 4, 4, 4)),
            pressure=np.zeros((4, 4, 4)),
            viscosity=-0.1,
        )
