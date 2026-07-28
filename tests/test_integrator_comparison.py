r"""Tests for the Gate 1 cross-integrator comparison experiment."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "experiments"))

from ns_certificate_lab._integrity import strict_json_loads, verify_digest
from ns_certificate_lab.nonlinear_cylinder import STEP_STREAM_FIELDS

import run_integrator_comparison as experiment


def _smoke_config() -> dict:
    return {
        "experiment_id": "integrator_comparison_smoke",
        "description": "smoke",
        "grid": {"nr": 17, "nz": 32},
        "initial_amplitude": 3000.0,
        "amplitude_scale": 1.0,
        "viscosity_schedule": [[0.0, 5.0e-3]],
        "t_final": 4.0e-5,
        "requested_time_steps": [2.0e-6, 1.0e-6],
        "methods": ["heun", "ssprk3", "rk4"],
        "diagnostic_stride": 10,
        "seed": 7,
        "acceptance": {
            "pairwise_amplification_relative_tolerance": 1.0e-3,
            "pairwise_max_u1_relative_tolerance": 1.0e-3,
            "pairwise_field_linf_relative_tolerance": 5.0e-3,
            "pairwise_differences_shrink_with_dt": True,
            "argmax_within_one_cell": True,
            "max_poisson_algebraic_residual_relative": 1.0e-12,
            "max_wall_u1_abs": 1.0e-12,
        },
        "preregistration_note": "smoke",
    }


def test_e29b_common_reference_matches_the_closed_form() -> None:
    assert experiment.E29B_MAX_CARTESIAN_VORTICITY_EXACT == pytest.approx(
        7569.6226982, rel=1.0e-9
    )


def test_config_validation_rejects_bad_configs() -> None:
    config = _smoke_config()
    experiment.validate_config(config)

    broken = _smoke_config()
    del broken["acceptance"]
    with pytest.raises(ValueError, match="missing keys"):
        experiment.validate_config(broken)

    broken = _smoke_config()
    broken["methods"] = ["heun", "leapfrog"]
    with pytest.raises(ValueError, match="unknown integrator"):
        experiment.validate_config(broken)

    broken = _smoke_config()
    broken["methods"] = ["ssprk3", "rk4"]
    with pytest.raises(ValueError, match="requires heun"):
        experiment.validate_config(broken)

    broken = _smoke_config()
    broken["requested_time_steps"] = [1.0e-6, 2.0e-6]
    with pytest.raises(ValueError, match="strictly decreasing"):
        experiment.validate_config(broken)

    broken = _smoke_config()
    broken["acceptance"]["extra"] = 1.0
    with pytest.raises(ValueError, match="exactly the keys"):
        experiment.validate_config(broken)


def test_shipping_config_is_valid() -> None:
    config = strict_json_loads(
        (REPOSITORY_ROOT / "configs" / "integrator_comparison.json").read_text(
            encoding="utf-8"
        ),
        label="shipping config",
    )
    experiment.validate_config(config)
    assert config["requested_time_steps"] == [6.0e-7, 3.0e-7]
    assert set(config["methods"]) == {"heun", "ssprk3", "rk4"}


@pytest.fixture(scope="module")
def smoke_bundle(tmp_path_factory):
    output_dir = tmp_path_factory.mktemp("integrator_comparison") / "bundle"
    summary = experiment.run(_smoke_config(), output_dir)
    return output_dir, summary


def test_smoke_run_passes_and_writes_a_verified_bundle(smoke_bundle) -> None:
    output_dir, summary = smoke_bundle
    assert summary["accepted_as_gate1_comparison"] is True
    assert all(summary["acceptance_checks"].values())
    # Manifest digests verify and cover every payload file.
    verify_digest(output_dir / "manifest.json")
    manifest = json.loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8")
    )
    for name, entry in manifest["files"].items():
        path = output_dir / name
        assert path.is_file()
        assert path.stat().st_size == entry["bytes"]
    listed = set(manifest["files"])
    on_disk = {
        p.relative_to(output_dir).as_posix()
        for p in output_dir.rglob("*")
        if p.is_file()
        and p.name not in ("manifest.json", "manifest.json.sha256")
    }
    assert on_disk == listed


def test_smoke_run_reports_both_normalizations(smoke_bundle) -> None:
    _, summary = smoke_bundle
    for level in summary["levels"]:
        for metrics in level["runs"].values():
            grid_norm = metrics["amplification_grid_normalized"]
            common = metrics["amplification_common_reference"]
            assert grid_norm > 0.0 and common > 0.0
            # Same final absolute value, different denominators.
            reconstructed = (
                metrics["final_max_cartesian_vorticity"]
                / metrics["initial_max_cartesian_vorticity_common_reference"]
            )
            assert common == pytest.approx(reconstructed, rel=1.0e-12)
            assert (
                metrics["initial_max_cartesian_vorticity_discrete"]
                != metrics[
                    "initial_max_cartesian_vorticity_common_reference"
                ]
            )


def test_smoke_run_stores_full_step_streams(smoke_bundle) -> None:
    output_dir, summary = smoke_bundle
    with np.load(output_dir / "step_streams.npz") as data:
        names = set(data.files)
        for level in summary["levels"]:
            dt_label = f"{level['requested_dt']:.3e}"
            for method, metrics in level["runs"].items():
                for field in STEP_STREAM_FIELDS:
                    key = f"{field}__{method}_dt{dt_label}"
                    assert key in names
                    assert data[key].size == metrics["steps"]


def test_zero_tolerance_fails_the_acceptance(smoke_bundle, tmp_path) -> None:
    config = _smoke_config()
    config["acceptance"]["pairwise_amplification_relative_tolerance"] = 0.0
    _, checks = experiment.evaluate(config)
    assert checks["pairwise_amplification_within_tolerance"] is False


def test_prepare_output_refuses_nonempty_directories(tmp_path) -> None:
    target = tmp_path / "occupied"
    target.mkdir()
    (target / "existing.txt").write_text("evidence", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing"):
        experiment._prepare_output(target)
