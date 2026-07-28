"""Tests for the independent Cartesian audit of saved Hou snapshots.

The real ``outputs/hou_early_time_v1`` checkpoints are read (never written) so
that the fault injections act on genuine saved data rather than on a
manufactured stand-in.  Audit grids are deliberately small here; the shipped
config uses 64x64x128 and is run by the orchestrator, not by this module.
"""

from __future__ import annotations

import ast
import copy
import dataclasses
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import experiments.run_hou_snapshot_cartesian_audit as audit_module
from experiments.run_hou_snapshot_cartesian_audit import (
    DENOMINATORS,
    EXPECTED_COORDINATE_SYSTEM,
    EXPECTED_EXPERIMENT_ID,
    PRIMITIVE_RESIDUAL_GAP,
    audit_cartesian_arrays,
    build_cartesian_grid,
    check_box_is_inscribed,
    primitive_residual_from_snapshot_pair,
    reconstruct_checkpoint,
    resolution_comparison,
    run,
    validate_config,
)
from ns_certificate_lab._integrity import (
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    verify_digest,
)
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.cartesian_validation import UniformCartesianGrid

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "hou_snapshot_cartesian_audit.json"
CHECKPOINT_DIR = REPOSITORY_ROOT / "outputs" / "hou_early_time_v1" / "checkpoints"
NR65_T000 = "outputs/hou_early_time_v1/checkpoints/checkpoint_nr65_nz128_t000.npz"

# Small audit grids: the fault injections only need enough resolution to
# separate a broken component map from the audit's own truncation error.
FAULT_SHAPE = (32, 32, 64)
GRID_SPEC: dict[str, Any] = {
    "half_width": 0.7,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_periodic": True,
}
AUDIT_KWARGS: dict[str, float | int] = {
    "interior_margin": 2,
    "relative_divergence_rms_gate": 0.02,
    "relative_curl_defect_rms_gate": 0.05,
    "relative_divergence_max_reference": 0.25,
    "relative_curl_defect_max_reference": 1.5,
}

# Measured on the real nr65 t000 checkpoint at the small smoke shape
# (24, 24, 48) with interior_margin 2: relative divergence RMS 1.27e-3 and
# relative curl-defect RMS 4.08e-3 against the gradient row-sum denominator.
# The recorded bounds below are deliberately about an order of magnitude
# looser so that they pin the order of magnitude, not the exact build.
SMOKE_RELATIVE_DIVERGENCE_BOUND = 1.0e-2
SMOKE_RELATIVE_CURL_BOUND = 3.0e-2


def _shipped_config() -> dict[str, Any]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="hou snapshot cartesian audit config",
    )
    assert isinstance(value, dict)
    return value


def _smoke_config() -> dict[str, Any]:
    config = _shipped_config()
    config["checkpoint_paths"] = [NR65_T000]
    config["cartesian_grid"] = {
        **config["cartesian_grid"],
        "shape": [24, 24, 48],
        "coarse_shape": [12, 12, 24],
    }
    return config


def _relative(metrics: dict[str, Any], name: str, region: str = "interior") -> float:
    value = metrics["regions"][region]["relative"][name]
    assert value is not None
    return float(value)


# --------------------------------------------------------------------------
# configuration contract
# --------------------------------------------------------------------------


def test_shipped_config_is_valid_and_names_the_real_checkpoints() -> None:
    config = _shipped_config()
    validate_config(config)
    assert config["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert config["coordinate_system"] == EXPECTED_COORDINATE_SYSTEM
    for entry in config["checkpoint_paths"]:
        assert (REPOSITORY_ROOT / entry).is_file(), entry
    # The default set is t000 and t004 (= T_1) at all three resolutions.
    assert len(config["checkpoint_paths"]) == 6
    assert sum("t000" in entry for entry in config["checkpoint_paths"]) == 3
    assert sum("t004" in entry for entry in config["checkpoint_paths"]) == 3


@pytest.mark.parametrize(
    "key",
    (
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "coordinate_system",
        "source_experiment",
        "checkpoint_paths",
        "cartesian_grid",
        "hard_gates",
        "record_only_thresholds",
    ),
)
def test_config_rejects_missing_key(key: str) -> None:
    config = _shipped_config()
    del config[key]
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(config)


def test_config_rejects_unknown_key() -> None:
    config = _shipped_config()
    config["extra_setting"] = 1
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(config)


def test_config_rejects_mislabeled_experiment_id() -> None:
    config = _shipped_config()
    config["experiment_id"] = "hou_early_time_v1"
    with pytest.raises(ValueError, match="experiment_id must be"):
        validate_config(config)


def test_config_rejects_mislabeled_coordinate_system() -> None:
    config = _shipped_config()
    config["coordinate_system"] = "axisymmetric cylindrical half-plane (r,z)"
    with pytest.raises(ValueError, match="coordinate_system"):
        validate_config(config)


def test_config_rejects_unsupported_schema_version() -> None:
    config = _shipped_config()
    config["schema_version"] = 2
    with pytest.raises(ValueError, match="schema_version"):
        validate_config(config)


@pytest.mark.parametrize(
    "entry",
    (
        "C:/absolute/checkpoint.npz",
        "/absolute/checkpoint.npz",
        "../outside/checkpoint.npz",
    ),
)
def test_config_rejects_escaping_checkpoint_path(entry: str) -> None:
    config = _shipped_config()
    config["checkpoint_paths"] = [entry]
    with pytest.raises(ValueError, match="repository-relative"):
        validate_config(config)


def test_config_rejects_non_npz_checkpoint_path() -> None:
    config = _shipped_config()
    config["checkpoint_paths"] = ["outputs/hou_early_time_v1/summary.json"]
    with pytest.raises(ValueError, match=r"\.npz"):
        validate_config(config)


def test_config_rejects_duplicate_checkpoint_paths() -> None:
    config = _shipped_config()
    config["checkpoint_paths"] = [NR65_T000, NR65_T000]
    with pytest.raises(ValueError, match="unique"):
        validate_config(config)


def test_config_rejects_coarse_shape_that_is_not_coarser() -> None:
    config = _shipped_config()
    config["cartesian_grid"]["coarse_shape"] = config["cartesian_grid"]["shape"]
    with pytest.raises(ValueError, match="strictly coarser"):
        validate_config(config)


def test_config_rejects_interior_margin_that_empties_the_norm_region() -> None:
    config = _shipped_config()
    config["cartesian_grid"]["coarse_shape"] = [6, 6, 12]
    config["cartesian_grid"]["interior_margin"] = 3
    with pytest.raises(ValueError, match="empty norm region"):
        validate_config(config)


def test_config_rejects_half_width_larger_than_the_unit_cylinder() -> None:
    config = _shipped_config()
    config["cartesian_grid"]["half_width"] = 1.5
    with pytest.raises(ValueError, match="half_width"):
        validate_config(config)


def test_audit_box_corners_must_stay_inside_the_source_cylinder() -> None:
    candidate = load_candidate(CHECKPOINT_DIR / "checkpoint_nr65_nz128_t000.npz")
    check_box_is_inscribed(candidate, GRID_SPEC)
    with pytest.raises(ValueError, match="leave the source cylinder"):
        check_box_is_inscribed(candidate, {**GRID_SPEC, "half_width": 0.95})


def test_cli_help_text_is_printable_on_a_legacy_console() -> None:
    """argparse prints the module docstring; a non-ASCII glyph crashes cp932."""

    assert audit_module.__doc__ is not None
    audit_module.__doc__.encode("ascii")
    with pytest.raises(SystemExit) as exit_info:
        audit_module.main(["--help"])
    assert exit_info.value.code == 0


def test_audit_never_imports_a_cylindrical_operator_module() -> None:
    """The whole point of this experiment is a second, independent path."""

    tree = ast.parse(inspect.getsource(audit_module))
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            modules.add(node.module or "")
    forbidden = {
        "ns_certificate_lab.operators",
        "ns_certificate_lab.pde",
        "ns_certificate_lab.poisson",
        "ns_certificate_lab.finite_cylinder_poisson",
        "ns_certificate_lab.nonlinear_cylinder",
        "ns_certificate_lab.cartesian",
    }
    assert not (modules & forbidden), sorted(modules & forbidden)


# --------------------------------------------------------------------------
# smoke run on a real checkpoint
# --------------------------------------------------------------------------


def test_smoke_run_on_real_nr65_t000_checkpoint(scratch_dir: Path) -> None:
    config = _smoke_config()
    output_dir = scratch_dir / "audit"
    summary = run(config, output_dir)

    assert summary["schema_version"] == 1
    assert summary["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert summary["source_experiment"] == "hou_early_time_v1"
    assert set(summary["denominators"]) == set(DENOMINATORS)
    assert summary["primitive_navier_stokes_residual_gap"] == PRIMITIVE_RESIDUAL_GAP
    assert summary["limitations"]
    assert summary["acceptance_checks"]["all_metrics_finite"]
    assert summary["acceptance_checks"]["all_points_finite"]
    assert summary["acceptance_checks"]["relative_divergence_rms_within_gate"]
    assert summary["acceptance_checks"]["relative_curl_defect_rms_within_gate"]
    assert summary["accepted_as_independent_cartesian_snapshot_audit"] is True

    (record,) = summary["checkpoints"]
    assert record["checkpoint"] == NR65_T000
    assert record["archive_sha256"] == sha256_file(REPOSITORY_ROOT / NR65_T000)
    assert record["source"] == {
        **record["source"],
        "nr": 65,
        "nz": 128,
        "physical_time": 0.0,
        "field_names": ["omega1", "psi1", "u1"],
    }
    assert set(record["evaluations"]) == {"fine", "coarse"}

    fine = record["evaluations"]["fine"]
    assert fine["shape"] == [24, 24, 48]
    assert set(fine["scales"]) == set(DENOMINATORS)
    for region in ("full", "interior"):
        block = fine["regions"][region]
        assert all(math.isfinite(value) for value in block["absolute"].values())
        finite_relative = [
            value for value in block["relative"].values() if value is not None
        ]
        # Only the signal-scale azimuthal ratio may be null (omega1 == 0 here).
        assert len(finite_relative) == len(block["relative"]) - 1
        assert all(math.isfinite(value) for value in finite_relative)
        for name in ("divergence_max", "curl_defect_max", "azimuthal_defect_max"):
            location = block["locations"][name]
            assert {"index", "x", "y", "z", "radius", "value"} <= set(location)
            assert len(location["index"]) == 3

    # omega1 vanishes identically in the E-29 initial datum, so r^2*omega1 is
    # exactly zero and no ratio against it is reported.
    assert fine["azimuthal_signal_degenerate"] is True
    assert (
        fine["regions"]["interior"]["relative"][
            "azimuthal_defect_max_over_signal_scale"
        ]
        is None
    )

    relative_divergence = _relative(fine, "divergence_rms_over_gradient_scale")
    relative_curl = _relative(fine, "curl_defect_rms_over_gradient_scale")
    assert 0.0 < relative_divergence < SMOKE_RELATIVE_DIVERGENCE_BOUND
    assert 0.0 < relative_curl < SMOKE_RELATIVE_CURL_BOUND

    refinement = record["audit_grid_refinement"]
    assert refinement["spacing_ratio"] > 1.0
    assert refinement["divergence_decreases_under_refinement"] is True
    assert refinement["curl_defect_decreases_under_refinement"] is True

    table = summary["resolution_comparison"]
    assert [entry["physical_time"] for entry in table] == [0.0]
    assert table[0]["rows"][0]["nr"] == 65

    names = sorted(path.name for path in output_dir.iterdir())
    assert names == [
        "config.snapshot.json",
        "config.snapshot.json.sha256",
        "manifest.json",
        "manifest.json.sha256",
        "metrics.csv",
        "summary.json",
        "summary.json.sha256",
    ]
    verify_digest(output_dir / "manifest.json")
    verify_digest(output_dir / "summary.json")
    verify_digest(output_dir / "config.snapshot.json")
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    for name, entry in manifest["files"].items():
        target = output_dir / name
        assert entry["sha256"] == sha256_file(target)
        assert entry["bytes"] == target.stat().st_size
    assert manifest["audited_inputs"][NR65_T000]["archive_sha256"] == (
        sha256_file(REPOSITORY_ROOT / NR65_T000)
    )
    snapshot = (output_dir / "config.snapshot.json").read_bytes()
    assert summary["reproducibility"]["config_sha256"] == sha256_bytes(snapshot)


def test_run_refuses_a_nonempty_output_directory(scratch_dir: Path) -> None:
    output_dir = scratch_dir / "occupied"
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(_smoke_config(), output_dir)
    assert (output_dir / "existing.txt").read_text(encoding="utf-8") == "keep me"


def test_run_refuses_to_write_inside_the_audited_source_directory() -> None:
    target = REPOSITORY_ROOT / "outputs" / "hou_early_time_v1" / "cartesian_audit"
    with pytest.raises(ValueError, match="audited source evidence directory"):
        run(_smoke_config(), target)
    assert not target.exists()


def test_resolution_comparison_groups_by_time_and_orders_by_resolution() -> None:
    def stub(nr: int, moment: float) -> dict[str, Any]:
        relative = {
            "divergence_rms_over_gradient_scale": 1.0,
            "divergence_max_over_gradient_scale": 2.0,
            "curl_defect_rms_over_gradient_scale": 3.0,
            "curl_defect_max_over_gradient_scale": 4.0,
            "azimuthal_defect_max_over_cancellation_scale": 5.0,
        }
        return {
            "checkpoint": f"nr{nr}-t{moment}",
            "source": {"nr": nr, "nz": 2 * nr, "physical_time": moment},
            "evaluations": {"fine": {"regions": {"interior": {"relative": relative}}}},
            "audit_grid_refinement": {
                "relative_divergence_rms_observed_order": 2.0,
                "relative_curl_defect_rms_observed_order": 2.0,
            },
        }

    table = resolution_comparison(
        [stub(193, 1.0), stub(65, 0.0), stub(129, 1.0), stub(65, 1.0), stub(129, 0.0)]
    )
    assert [entry["physical_time"] for entry in table] == [0.0, 1.0]
    assert [row["nr"] for row in table[0]["rows"]] == [65, 129]
    assert [row["nr"] for row in table[1]["rows"]] == [65, 129, 193]


# --------------------------------------------------------------------------
# fault injection on the real snapshot
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def real_reconstruction():
    """Real nr65 t000 checkpoint reconstructed once on the fault-audit box."""

    candidate = load_candidate(CHECKPOINT_DIR / "checkpoint_nr65_nz128_t000.npz")
    grid = build_cartesian_grid(GRID_SPEC, FAULT_SHAPE)
    arrays = reconstruct_checkpoint(candidate, grid)
    baseline = audit_cartesian_arrays(
        grid,
        velocity=arrays.velocity,
        vorticity_e18b=arrays.vorticity_e18b,
        **AUDIT_KWARGS,
    )
    return candidate, grid, arrays, baseline


def test_uncorrupted_real_snapshot_passes_the_same_audit(real_reconstruction) -> None:
    _, _, _, baseline = real_reconstruction
    assert _relative(baseline, "divergence_rms_over_gradient_scale") < 0.02
    assert _relative(baseline, "curl_defect_rms_over_gradient_scale") < 0.05
    assert baseline["cartesian_validation_gate"]["divergence_rms_passed"]
    assert baseline["cartesian_validation_gate"]["curl_rms_passed"]


def test_detects_sign_flipped_u_y_component(real_reconstruction) -> None:
    """A wrong swirl sign in u_y breaks E-18b by far more than 10x."""

    _, grid, arrays, baseline = real_reconstruction
    corrupted = arrays.velocity.copy()
    corrupted[1] = -corrupted[1]
    metrics = audit_cartesian_arrays(
        grid,
        velocity=corrupted,
        vorticity_e18b=arrays.vorticity_e18b,
        **AUDIT_KWARGS,
    )
    curl_ratio = _relative(metrics, "curl_defect_rms_over_gradient_scale") / _relative(
        baseline, "curl_defect_rms_over_gradient_scale"
    )
    divergence_ratio = _relative(
        metrics, "divergence_rms_over_gradient_scale"
    ) / _relative(baseline, "divergence_rms_over_gradient_scale")
    assert curl_ratio > 10.0, curl_ratio
    assert divergence_ratio > 10.0, divergence_ratio
    assert not metrics["cartesian_validation_gate"]["curl_rms_passed"]
    assert not metrics["cartesian_validation_gate"]["divergence_rms_passed"]
    # The azimuthal agreement collapses from machine precision to O(1).
    assert (
        _relative(metrics, "azimuthal_defect_max_over_cancellation_scale") > 0.1
    )


def test_detects_swapped_component_mapping(real_reconstruction) -> None:
    """u_x <-> u_y is a wrong Cartesian component map, not a sign error."""

    _, grid, arrays, baseline = real_reconstruction
    corrupted = arrays.velocity.copy()
    corrupted[[0, 1]] = corrupted[[1, 0]]
    metrics = audit_cartesian_arrays(
        grid,
        velocity=corrupted,
        vorticity_e18b=arrays.vorticity_e18b,
        **AUDIT_KWARGS,
    )
    assert _relative(metrics, "curl_defect_rms_over_gradient_scale") / _relative(
        baseline, "curl_defect_rms_over_gradient_scale"
    ) > 10.0
    assert _relative(metrics, "divergence_rms_over_gradient_scale") / _relative(
        baseline, "divergence_rms_over_gradient_scale"
    ) > 10.0
    assert not metrics["cartesian_validation_gate"]["curl_rms_passed"]


def _axis_perturbations(
    candidate,
    *,
    amplitude_fraction: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an odd-in-r psi1 kink and a matched even-in-r control pair.

    ``phi(z) = cos(2 pi z)`` is a smooth periodic axial profile.  The odd
    perturbation ``c*r*phi`` violates E-16c (``psi1_r(0,z) = c*phi != 0``);
    the even perturbation ``c*(1-r^2)^2*phi`` satisfies it.  For the even case
    the matching ``omega1`` correction is available in closed form because
    ``L5[(1-r^2)^2] = -16 + 24 r^2`` and ``omega1 = -L5 psi1``.
    """

    radius = candidate.grid.r[:, None]
    axial = candidate.grid.z[None, :]
    phi = np.cos(2.0 * np.pi * axial)
    phi_zz = -((2.0 * np.pi) ** 2) * phi
    amplitude = amplitude_fraction * float(np.max(np.abs(candidate.fields["u1"])))
    odd_psi1 = amplitude * radius * phi
    even_psi1 = amplitude * (1.0 - radius**2) ** 2 * phi
    even_omega1 = -amplitude * (
        (-16.0 + 24.0 * radius**2) * phi + (1.0 - radius**2) ** 2 * phi_zz
    )
    return odd_psi1, even_psi1, even_omega1


def _near_axis_curl_defect_rms(
    grid: UniformCartesianGrid,
    velocity: np.ndarray,
    vorticity_e18b: np.ndarray,
    *,
    cells: float = 3.0,
) -> float:
    """RMS curl defect on the column of cells closest to the axis."""

    from ns_certificate_lab.cartesian_validation import cartesian_curl

    x, y, _ = grid.mesh()
    mask = np.hypot(x, y) <= cells * grid.spacings[0]
    defect = cartesian_curl(grid, velocity) - vorticity_e18b
    return float(np.sqrt(np.mean(defect[:, mask] ** 2)))


def test_detects_axis_condition_break_odd_radial_psi1(real_reconstruction) -> None:
    r"""An odd-in-``r`` ``psi1`` kink is invisible to E-02 but breaks E-18b.

    Why this is a genuine axis-condition violation and not an arbitrary edit:

    * The E-18a map makes ``div u`` identically zero for *every* ``(u1,psi1)``
      pair, so the divergence check is structurally blind to any change of
      ``psi1``.  It must stay flat, and it does.
    * ``psi1 -> psi1 + c*r*phi(z)`` gives ``psi1_r(0,z) = c*phi(z) != 0``,
      violating the E-16c pole condition.  The consequence is that
      ``L5(c*r*phi) = 3*c*phi/r + c*r*phi''`` is unbounded at the axis, so
      **no finite stored omega1 array can satisfy ``-L5 psi1 = omega1``** for
      the perturbed streamfunction.  The reconstructed Cartesian curl acquires
      a direction-dependent, non-vanishing value on the axis, which is exactly
      the loss of Cartesian smoothness that E-16 exists to forbid.
    * The full-curl check therefore must reject it, and the companion test
      below shows an *axis-regular* perturbation of the same amplitude, with
      its analytically matched omega1, is accepted.
    """

    candidate, grid, clean, baseline = real_reconstruction
    odd_psi1, _, _ = _axis_perturbations(candidate)
    corrupted = dataclasses.replace(
        candidate,
        fields={
            **candidate.fields,
            "psi1": candidate.fields["psi1"] + odd_psi1,
        },
    )
    arrays = reconstruct_checkpoint(corrupted, grid)
    metrics = audit_cartesian_arrays(
        grid,
        velocity=arrays.velocity,
        vorticity_e18b=arrays.vorticity_e18b,
        **AUDIT_KWARGS,
    )
    curl_ratio = _relative(metrics, "curl_defect_rms_over_gradient_scale") / _relative(
        baseline, "curl_defect_rms_over_gradient_scale"
    )
    divergence_ratio = _relative(
        metrics, "divergence_rms_over_gradient_scale"
    ) / _relative(baseline, "divergence_rms_over_gradient_scale")
    assert curl_ratio > 10.0, curl_ratio
    assert not metrics["cartesian_validation_gate"]["curl_rms_passed"]
    # E-02 is an identity of the component map, so it cannot see this fault.
    assert divergence_ratio < 2.0, divergence_ratio
    assert metrics["cartesian_validation_gate"]["divergence_rms_passed"]
    # The near-axis column degrades too: the perturbed curl no longer has a
    # single limit on the axis.  The interpolated baseline defect is already
    # largest there, so this ratio is smaller than the global one, but it is
    # the companion control below that makes it meaningful.
    axis_ratio = _near_axis_curl_defect_rms(
        grid, arrays.velocity, arrays.vorticity_e18b
    ) / _near_axis_curl_defect_rms(grid, clean.velocity, clean.vorticity_e18b)
    assert axis_ratio > 2.0, axis_ratio


def test_accepts_axis_regular_perturbation_with_matched_omega1(
    real_reconstruction,
) -> None:
    """Control for the test above: the audit is not merely psi1-sensitive."""

    candidate, grid, clean, baseline = real_reconstruction
    _, even_psi1, even_omega1 = _axis_perturbations(candidate)
    perturbed = dataclasses.replace(
        candidate,
        fields={
            **candidate.fields,
            "psi1": candidate.fields["psi1"] + even_psi1,
            "omega1": candidate.fields["omega1"] + even_omega1,
        },
    )
    arrays = reconstruct_checkpoint(perturbed, grid)
    metrics = audit_cartesian_arrays(
        grid,
        velocity=arrays.velocity,
        vorticity_e18b=arrays.vorticity_e18b,
        **AUDIT_KWARGS,
    )
    curl_ratio = _relative(metrics, "curl_defect_rms_over_gradient_scale") / _relative(
        baseline, "curl_defect_rms_over_gradient_scale"
    )
    assert curl_ratio < 1.25, curl_ratio
    assert metrics["cartesian_validation_gate"]["curl_rms_passed"]
    assert metrics["cartesian_validation_gate"]["divergence_rms_passed"]
    axis_ratio = _near_axis_curl_defect_rms(
        grid, arrays.velocity, arrays.vorticity_e18b
    ) / _near_axis_curl_defect_rms(grid, clean.velocity, clean.vorticity_e18b)
    assert axis_ratio < 1.1, axis_ratio


def test_detects_vorticity_oracle_sign_fault(real_reconstruction) -> None:
    """A flipped E-18b oracle must be rejected by the curl check alone."""

    _, grid, arrays, baseline = real_reconstruction
    metrics = audit_cartesian_arrays(
        grid,
        velocity=arrays.velocity,
        vorticity_e18b=-arrays.vorticity_e18b,
        **AUDIT_KWARGS,
    )
    assert _relative(metrics, "curl_defect_rms_over_gradient_scale") / _relative(
        baseline, "curl_defect_rms_over_gradient_scale"
    ) > 10.0
    assert not metrics["cartesian_validation_gate"]["curl_rms_passed"]
    assert metrics["cartesian_validation_gate"]["divergence_rms_passed"]


def test_audit_rejects_a_degenerate_zero_velocity_field() -> None:
    grid = build_cartesian_grid(GRID_SPEC, (8, 8, 8))
    with pytest.raises(ValueError, match="no denominator"):
        audit_cartesian_arrays(
            grid,
            velocity=np.zeros((3, 8, 8, 8)),
            vorticity_e18b=np.zeros((3, 8, 8, 8)),
            **AUDIT_KWARGS,
        )


# --------------------------------------------------------------------------
# the primitive residual: manufactured pair only
# --------------------------------------------------------------------------

_VISCOSITY = 0.07


def _manufactured_pair_fields(
    grid: UniformCartesianGrid,
    time: float,
) -> dict[str, np.ndarray]:
    """Closed-form decaying Taylor-Green-like field with an exact pressure."""

    x, y, z = grid.mesh()
    amplitude = math.exp(-time)
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
    laplacian = -6.0 * velocity
    advection = np.stack(
        (
            2.0 * amplitude**2 * sin_2x * cos_2x * cos_z**2,
            4.0 * amplitude**2 * sin_y * cos_y * cos_z**2,
            np.zeros_like(x),
        ),
        axis=0,
    )
    return {
        "velocity": velocity,
        "velocity_t": velocity_t,
        "pressure": pressure,
        "body_force": (
            velocity_t + advection + pressure_gradient - _VISCOSITY * laplacian
        ),
    }


def _periodic_grid(count: int) -> UniformCartesianGrid:
    return UniformCartesianGrid.uniform(
        shape=(count, count, count),
        x_bounds=(-np.pi, np.pi),
        y_bounds=(-np.pi, np.pi),
        z_bounds=(-np.pi, np.pi),
        periodic=(True, True, True),
    )


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=np.float64) ** 2)))


def test_primitive_residual_pair_reproduces_a_manufactured_body_force() -> None:
    """The pair assembly recovers f on a field where p and u_t are exact."""

    grid = _periodic_grid(32)
    midpoint = 0.25
    step = 0.01
    early = _manufactured_pair_fields(grid, midpoint - 0.5 * step)
    late = _manufactured_pair_fields(grid, midpoint + 0.5 * step)
    exact = _manufactured_pair_fields(grid, midpoint)

    assembled = primitive_residual_from_snapshot_pair(
        grid,
        velocity_early=early["velocity"],
        velocity_late=late["velocity"],
        time_early=midpoint - 0.5 * step,
        time_late=midpoint + 0.5 * step,
        midpoint_pressure=exact["pressure"],
        viscosity=_VISCOSITY,
        body_force=exact["body_force"],
    )
    assert assembled.midpoint_time == pytest.approx(midpoint)
    assert assembled.time_step == pytest.approx(step)
    residual = assembled.residual
    assert np.allclose(
        residual.total,
        residual.time_derivative
        + residual.advection
        + residual.pressure_gradient
        + residual.viscous,
    )
    assert np.allclose(
        assembled.velocity_time_derivative,
        (late["velocity"] - early["velocity"]) / step,
    )
    signal = _rms(exact["body_force"])
    assert _rms(residual.defect_against_body_force) / signal < 0.02
    # The manufactured field is analytically divergence free; what remains is
    # the O(h^2) level of a wavenumber-2 field on 32 points across 2*pi,
    # which is about (2h)^2/6 = 2.6e-2 relative.
    assert _rms(residual.divergence) / _rms(exact["velocity"]) < 0.05


def test_primitive_residual_pair_is_second_order_in_checkpoint_spacing() -> None:
    """Isolate the temporal error by differencing against the exact midpoint."""

    from ns_certificate_lab.cartesian_validation import primitive_ns_residual

    grid = _periodic_grid(16)
    midpoint = 0.25
    exact = _manufactured_pair_fields(grid, midpoint)
    # Same grid and same spatial stencils, but the exact midpoint velocity and
    # the exact time derivative.  Differencing against this cancels every
    # spatial truncation term and leaves only the temporal error of the pair.
    exact_residual = primitive_ns_residual(
        grid,
        velocity=exact["velocity"],
        velocity_t=exact["velocity_t"],
        pressure=exact["pressure"],
        viscosity=_VISCOSITY,
    )

    errors: list[float] = []
    steps = [0.2, 0.1, 0.05]
    for step in steps:
        early = _manufactured_pair_fields(grid, midpoint - 0.5 * step)
        late = _manufactured_pair_fields(grid, midpoint + 0.5 * step)
        assembled = primitive_residual_from_snapshot_pair(
            grid,
            velocity_early=early["velocity"],
            velocity_late=late["velocity"],
            time_early=midpoint - 0.5 * step,
            time_late=midpoint + 0.5 * step,
            midpoint_pressure=exact["pressure"],
            viscosity=_VISCOSITY,
        )
        errors.append(_rms(assembled.residual.total - exact_residual.total))

    assert errors[0] > errors[1] > errors[2] > 0.0
    orders = [
        math.log(errors[index] / errors[index + 1])
        / math.log(steps[index] / steps[index + 1])
        for index in range(len(errors) - 1)
    ]
    assert min(orders) > 1.8, (orders, errors)


def test_primitive_residual_pair_requires_an_increasing_time_pair() -> None:
    grid = _periodic_grid(8)
    fields = _manufactured_pair_fields(grid, 0.0)
    for time_early, time_late in ((0.5, 0.5), (0.5, 0.4)):
        with pytest.raises(ValueError, match="strictly greater"):
            primitive_residual_from_snapshot_pair(
                grid,
                velocity_early=fields["velocity"],
                velocity_late=fields["velocity"],
                time_early=time_early,
                time_late=time_late,
                midpoint_pressure=fields["pressure"],
                viscosity=_VISCOSITY,
            )


def test_primitive_residual_pair_rejects_mismatched_shapes() -> None:
    grid = _periodic_grid(8)
    fields = _manufactured_pair_fields(grid, 0.0)
    with pytest.raises(ValueError, match="shape"):
        primitive_residual_from_snapshot_pair(
            grid,
            velocity_early=fields["velocity"],
            velocity_late=np.zeros((3, 8, 8, 7)),
            time_early=0.0,
            time_late=0.1,
            midpoint_pressure=fields["pressure"],
            viscosity=_VISCOSITY,
        )


def test_primitive_residual_gap_is_declared_not_silently_dropped() -> None:
    """The one thing a single checkpoint cannot support stays visible."""

    assert PRIMITIVE_RESIDUAL_GAP["status"].startswith("NOT evaluated")
    reasons = " ".join(PRIMITIVE_RESIDUAL_GAP["missing_inputs"]).lower()
    assert "pressure" in reasons
    assert "adjacent" in reasons
    # No shipped code path feeds a real checkpoint into the pair assembly.
    source = inspect.getsource(audit_module)
    body = source.split("def primitive_residual_from_snapshot_pair", 1)[1]
    assert "load_candidate" not in body


def test_shipped_config_and_docs_agree_on_the_gap() -> None:
    config = _shipped_config()
    assert "primitive" in config["interpretation"].lower()
    assert "gap" in config["interpretation"].lower()
    assert copy.deepcopy(PRIMITIVE_RESIDUAL_GAP) == PRIMITIVE_RESIDUAL_GAP
