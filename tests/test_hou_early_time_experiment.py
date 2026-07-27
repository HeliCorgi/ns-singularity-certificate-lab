"""Acceptance tests for the early-time Hou finite-cylinder experiment.

The real E-30 target time is deliberately never integrated here: these tests
exercise the configuration contract, the artifact bundle and the E-29b initial
norms on grids small enough for continuous integration.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from experiments.run_hou_early_time import (
    E29B_MAX_ABS_U1,
    E29B_MAX_CARTESIAN_VORTICITY,
    E29_AMPLITUDE,
    build_grid,
    evolve_resolution,
    initial_norms,
    run,
    validate_config,
)
from ns_certificate_lab._integrity import strict_json_loads
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.nonlinear_cylinder import load_checkpoint


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "hou_early_time.json"


def _config() -> dict[str, Any]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="hou early time test config",
    )
    assert isinstance(value, dict)
    return value


def _smoke_config() -> dict[str, Any]:
    """A few steps on tiny grids: the full E-30 target is never run in tests."""

    config = _config()
    config["resolutions"] = [[17, 32], [25, 48]]
    config["t_final"] = 2.0e-6
    config["snapshot_times"] = [0.0, 1.0e-6, 2.0e-6]
    config["max_time_step"] = 2.5e-7
    config["max_steps"] = 200
    config["diagnostic_stride"] = 1
    config["acceptance"] = {
        # Loosened only because the smoke grids are far too coarse to resolve
        # the E-29b maxima; the shipped config keeps the 1 percent bound.
        "maximum_circulation_growth_ratio": 1.0e-3,
        "maximum_energy_growth_ratio": 1.0e-8,
        "maximum_initial_norm_relative_error": 0.5,
        "maximum_odd_symmetry_defect_ratio": 1.0e-6,
    }
    return config


@pytest.fixture(scope="module")
def generated_run(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    config = _smoke_config()
    output = tmp_path_factory.mktemp("hou-early-time") / "evidence"
    summary = run(config, output)
    return config, output, summary


# ------------------------------------------------------ config validation


def test_shipped_config_is_valid_and_targets_the_audited_protocol() -> None:
    config = _config()
    validate_config(config)
    assert config["experiment_id"] == "hou_early_time_v1"
    assert config["r_max"] == 1.0
    assert config["z_period"] == 1.0
    assert config["amplitude_scale"] == 1.0
    assert config["viscosity_schedule"] == [[0.0, 5e-4], [0.00227375, 5e-3]]
    assert config["cfl_coefficient"] == 0.1
    assert len(config["resolutions"]) >= 2
    assert config["t_final"] == pytest.approx(0.002191729)
    assert config["snapshot_times"][-1] == pytest.approx(config["t_final"])


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"experiment_id": "something_else"}, "experiment_id must be"),
        ({"schema_version": 2}, "schema_version must be 1"),
        ({"resolutions": [[129, 256]]}, "at least two resolutions"),
        ({"resolutions": [[129, 256], [65, 128]]}, "strictly increasing"),
        ({"resolutions": [[129, 256], [193]]}, "\\[nr, nz\\] integer pair"),
        ({"amplitude_scale": 0.0}, "amplitude_scale must be positive"),
        ({"cfl_coefficient": -1.0}, "cfl_coefficient must be positive"),
        ({"max_steps": 0}, "max_steps must be a positive integer"),
        ({"diagnostic_stride": 1.5}, "diagnostic_stride must be a positive"),
        ({"t_final": 0.0}, "t_final must be positive"),
        ({"snapshot_times": []}, "snapshot_times must be a nonempty list"),
        ({"snapshot_times": [0.001, 0.0005]}, "strictly increasing"),
        ({"snapshot_times": [0.0, 0.001]}, "last snapshot time must be t_final"),
        ({"snapshot_times": [0.0, 9.0]}, "snapshot times must lie"),
        ({"viscosity_schedule": [[0.001, 5e-4]]}, "first viscosity switch time"),
        ({"viscosity_schedule": [[0.0, -1.0]]}, "viscosity values must be"),
        ({"units": {}}, "units must be a nonempty object"),
        ({"description": "  "}, "description must be a nonempty string"),
        ({"seed": "abc"}, "seed must be an integer"),
    ],
)
def test_invalid_config_values_are_rejected(
    mutation: dict[str, Any],
    message: str,
) -> None:
    config = _config()
    config.update(mutation)
    with pytest.raises(ValueError, match=message):
        validate_config(config)


def test_missing_and_unknown_config_keys_are_rejected() -> None:
    missing = _config()
    del missing["viscosity_schedule"]
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(missing)

    unknown = _config()
    unknown["filter_strength"] = 0.1
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(unknown)

    bad_acceptance = _config()
    bad_acceptance["acceptance"] = dict(bad_acceptance["acceptance"])
    bad_acceptance["acceptance"]["extra"] = 1.0
    with pytest.raises(ValueError, match="acceptance object has missing or unknown"):
        validate_config(bad_acceptance)

    missing_acceptance = _config()
    missing_acceptance["acceptance"] = dict(missing_acceptance["acceptance"])
    del missing_acceptance["acceptance"]["maximum_energy_growth_ratio"]
    with pytest.raises(ValueError, match="acceptance object has missing or unknown"):
        validate_config(missing_acceptance)


# ------------------------------------------------------------ E-29b norms


def test_e29b_initial_norms_are_recovered_on_a_moderate_grid() -> None:
    config = _config()
    grid = build_grid(config, nr=129, nz=256)
    norms = initial_norms(grid, amplitude_scale=1.0)
    assert norms["expected_max_abs_u1"] == E29B_MAX_ABS_U1
    assert norms["expected_max_cartesian_vorticity"] == E29B_MAX_CARTESIAN_VORTICITY
    # Measured on (129,256): 3265.5332 and 7554.6095, i.e. 1.39e-4 and 1.98e-3
    # relative.  Both are grid maxima of smooth fields and converge from below.
    assert norms["max_abs_u1_relative_error"] < 0.01
    assert norms["max_cartesian_vorticity_relative_error"] < 0.01
    assert norms["max_abs_u1"] < E29B_MAX_ABS_U1
    assert norms["max_cartesian_vorticity"] < E29B_MAX_CARTESIAN_VORTICITY

    # The norms scale linearly with amplitude_scale.
    half = initial_norms(grid, amplitude_scale=0.5)
    assert half["max_abs_u1"] == pytest.approx(0.5 * norms["max_abs_u1"])
    assert half["max_abs_u1_relative_error"] == pytest.approx(
        norms["max_abs_u1_relative_error"]
    )
    assert E29_AMPLITUDE == 12000.0


def test_initial_vorticity_norm_converges_under_refinement() -> None:
    config = _config()
    errors = [
        initial_norms(
            build_grid(config, nr=nr, nz=nz),
            amplitude_scale=1.0,
        )["max_cartesian_vorticity_relative_error"]
        for nr, nz in ((65, 128), (129, 256))
    ]
    assert errors[0] > errors[1]
    assert errors[1] < 0.5 * errors[0]


# ------------------------------------------------------- end-to-end smoke


def test_smoke_run_produces_a_complete_accepted_bundle(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, summary = generated_run
    assert summary["experiment_id"] == "hou_early_time_v1"
    assert summary["accepted_as_early_time_observation"] is True
    expected_checks = {
        "all_resolutions_completed",
        "all_diagnostics_finite",
        "energy_non_increasing",
        "circulation_max_principle",
        "early_swirl_maximum_decreases",
        "odd_symmetry_preserved",
        "initial_norms_match_e29b",
        "cross_solver_elliptic_agreement_recorded",
    }
    assert set(summary["acceptance_checks"]) == expected_checks
    assert all(summary["acceptance_checks"].values())
    assert "not a reproduction" in summary["interpretation"]
    assert summary["derived_reference_norms"]["max_abs_u1"] == E29B_MAX_ABS_U1
    assert summary["reproducibility"]["seed"] == config["seed"]
    assert summary["reproducibility"]["runtime_provenance"][
        "source_fingerprint_sha256"
    ]

    limitations = " ".join(summary["limitations"])
    assert "uniform fixed grid" in limitations
    assert "numerical observation" in limitations
    assert "not a reproduction claim" in limitations
    assert "monitored" in limitations

    resolutions = summary["resolutions"]
    assert [(item["nr"], item["nz"]) for item in resolutions] == [
        tuple(pair) for pair in config["resolutions"]
    ]
    for item in resolutions:
        assert item["completed"] is True
        assert item["failure"] is None
        assert item["all_diagnostics_finite"] is True
        assert item["final_time"] == pytest.approx(config["t_final"])
        assert item["early_swirl_maximum_decreases"] is True
        assert item["maximum_energy_growth_ratio"] <= 0.0
        assert math.isfinite(item["final_amplification"])
        assert math.isfinite(
            item["snapshot_cross_solver_psi_max_abs_difference"]
        )
        assert item["maximum_wall_u1_abs"] == 0.0


def test_smoke_run_artifacts_are_finite_and_checksummed(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, _ = generated_run
    with (output / "diagnostics.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    required = {
        "nr",
        "nz",
        "time",
        "dt",
        "viscosity",
        "energy",
        "enstrophy",
        "max_abs_u1",
        "max_abs_omega1",
        "max_cartesian_vorticity",
        "amplification",
        "argmax_u1_r",
        "argmax_u1_z",
        "circulation_max",
        "divergence_residual_max",
        "odd_symmetry_defect",
        "axis_parity_defect",
    }
    assert required <= rows[0].keys()
    assert all(
        math.isfinite(float(value)) for row in rows for value in row.values()
    )

    with (output / "snapshots.csv").open(newline="", encoding="utf-8") as stream:
        snapshots = list(csv.DictReader(stream))
    assert len(snapshots) == len(config["resolutions"]) * len(config["snapshot_times"])
    assert {"psi_cross_solver_max_abs_difference", "checkpoint"} <= snapshots[
        0
    ].keys()
    assert all(
        math.isfinite(float(row["psi_cross_solver_max_abs_difference"]))
        for row in snapshots
    )

    with np.load(output / "trajectories.npz", allow_pickle=False) as archive:
        assert any(name.startswith("psi1_final_") for name in archive.files)
        assert all(np.all(np.isfinite(archive[name])) for name in archive.files)

    manifest = strict_json_loads(
        (output / "manifest.json").read_text(encoding="utf-8"),
        label="hou early time manifest",
    )
    assert {
        "config.snapshot.json",
        "diagnostics.csv",
        "snapshots.csv",
        "summary.json",
        "trajectories.npz",
    } <= set(manifest["files"])
    assert any(name.startswith("checkpoints/") for name in manifest["files"])
    for name, record in manifest["files"].items():
        payload = (output / name).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == record["sha256"]
        assert len(payload) == record["bytes"]
    assert (output / "manifest.json.sha256").read_text(
        encoding="ascii"
    ).strip() == hashlib.sha256((output / "manifest.json").read_bytes()).hexdigest()


def test_snapshot_checkpoints_reload_as_v2_candidates(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, _ = generated_run
    archives = sorted((output / "checkpoints").glob("*.npz"))
    assert len(archives) == len(config["resolutions"]) * len(config["snapshot_times"])
    artifact = load_candidate(archives[-1])
    assert set(artifact.fields) == {"u1", "omega1", "psi1"}
    assert artifact.schema.endswith("/v2")
    assert artifact.description is not None
    assert artifact.provenance is not None
    checkpoint = load_checkpoint(archives[-1])
    # The stored psi1 is the elliptic solution of the stored omega1, so
    # re-constraining the checkpoint must reproduce it bit for bit.
    assert np.array_equal(checkpoint.state.psi1, artifact.fields["psi1"])
    assert checkpoint.time == pytest.approx(config["snapshot_times"][-1])
    assert checkpoint.viscosity == pytest.approx(config["viscosity_schedule"][0][1])
    assert checkpoint.metadata["experiment_id"] == "hou_early_time_v1"
    assert np.all(np.isfinite(checkpoint.state.u1))
    assert float(np.max(np.abs(checkpoint.state.u1[-1]))) == 0.0


def test_nonempty_output_directory_is_refused(scratch_dir: Path) -> None:
    config = _smoke_config()
    occupied = scratch_dir / "occupied"
    occupied.mkdir()
    marker = occupied / "preserve.txt"
    marker.write_text("keep\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(config, occupied)
    assert marker.read_text(encoding="utf-8") == "keep\n"
    assert list(occupied.iterdir()) == [marker]


def test_max_steps_records_an_explicit_incomplete_resolution() -> None:
    config = _smoke_config()
    config["max_steps"] = 3
    result = evolve_resolution(config, nr=17, nz=32)
    assert result.completed is False
    assert result.failure is not None
    assert "before t_final" in result.failure
    assert result.step_count == 3
