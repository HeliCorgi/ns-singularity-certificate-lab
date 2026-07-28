"""Tests for the P0-D core-width audit of the saved Hou snapshots.

The heavy assertions run against a *synthetic* evidence bundle built in
``tmp_path`` by the experiment module's own writer, so the audit's format
contract is exercised end to end without depending on shipped simulation data.
A final test reads the real ``outputs/core_width_audit_v1`` bundle when it is
present and skips when it is not, so a checkout without the evidence still
passes.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pytest

import experiments.run_core_width_audit as audit_module
from experiments.run_core_width_audit import (
    AMPLIFICATION_STATEMENT,
    AUDITED_FIELDS,
    E29B_MAX_CARTESIAN_VORTICITY_EXACT,
    EXPECTED_EXPERIMENT_ID,
    MANIFEST_SCHEMA,
    _prepare_output,
    continuum_initial_vorticity_maximum,
    evidence_directory,
    run,
    validate_config,
    verify_evidence_bundle,
    write_minimal_evidence_bundle,
)
from ns_certificate_lab._integrity import (
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    verify_digest,
)
from ns_certificate_lab.core_width import PREREGISTERED_MIN_POINTS_PER_FRONT

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "core_width_audit.json"
REAL_OUTPUT_DIR = REPOSITORY_ROOT / "outputs" / "core_width_audit_v1"

SYNTHETIC_RESOLUTIONS = ((17, 32), (33, 64))
SYNTHETIC_TIMES = (0.0, 0.001)
EXPECTED_ARTIFACTS = [
    "amplification_renormalization.csv",
    "common_grid_differences.csv",
    "config.snapshot.json",
    "config.snapshot.json.sha256",
    "fit_precondition.csv",
    "manifest.json",
    "manifest.json.sha256",
    "peak_resolution_reports.csv",
    "summary.json",
    "summary.json.sha256",
]


def _shipped_config() -> dict[str, Any]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="core width audit config",
    )
    assert isinstance(value, dict)
    return value


def _synthetic_bundle(scratch_dir: Path) -> Path:
    bundle = scratch_dir / "synthetic_hou_v0"
    write_minimal_evidence_bundle(
        bundle,
        experiment_id="synthetic_core_width_source_v0",
        resolutions=SYNTHETIC_RESOLUTIONS,
        snapshot_times=SYNTHETIC_TIMES,
    )
    return bundle


def _synthetic_config(bundle: Path) -> dict[str, Any]:
    config = _shipped_config()
    config["input_evidence_dirs"] = [
        bundle.resolve().relative_to(REPOSITORY_ROOT).as_posix()
    ]
    config["snapshot_labels"] = ["t000", "t001"]
    config["initial_snapshot_label"] = "t000"
    config["final_snapshot_label"] = "t001"
    config["continuum_reference"] = {
        "search_points_per_axis": 201,
        "refinement_stages": 2,
    }
    return config


# --------------------------------------------------------------------------
# configuration contract
# --------------------------------------------------------------------------


def test_shipped_config_is_valid_and_names_the_real_bundles() -> None:
    config = _shipped_config()
    validate_config(config)
    assert config["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert config["audited_fields"] == list(AUDITED_FIELDS)
    assert config["snapshot_labels"] == ["t000", "t001", "t002", "t003", "t004"]
    assert config["initial_snapshot_label"] == "t000"
    assert config["final_snapshot_label"] == "t004"
    for entry in config["input_evidence_dirs"]:
        directory = evidence_directory(entry)
        assert directory.is_dir(), entry
        assert (directory / "manifest.json").is_file(), entry


@pytest.mark.parametrize(
    "key",
    (
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "input_evidence_dirs",
        "snapshot_labels",
        "initial_snapshot_label",
        "final_snapshot_label",
        "audited_fields",
        "continuum_reference",
    ),
)
def test_config_rejects_missing_key(key: str) -> None:
    config = _shipped_config()
    del config[key]
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(config)


def test_config_rejects_unknown_key() -> None:
    config = _shipped_config()
    config["min_points_per_front"] = 7
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(config)


def test_config_rejects_mislabeled_experiment_id() -> None:
    config = _shipped_config()
    config["experiment_id"] = "hou_early_time_v1"
    with pytest.raises(ValueError, match="experiment_id must be"):
        validate_config(config)


@pytest.mark.parametrize(
    "entry",
    (
        "C:/absolute/outputs",
        "/absolute/outputs",
        "../outside/outputs",
    ),
)
def test_config_rejects_escaping_input_evidence_dir(entry: str) -> None:
    config = _shipped_config()
    config["input_evidence_dirs"] = [entry]
    with pytest.raises(ValueError, match="repository-relative"):
        validate_config(config)


def test_config_rejects_duplicate_input_evidence_dirs() -> None:
    config = _shipped_config()
    config["input_evidence_dirs"] = ["outputs/hou_early_time_v1"] * 2
    with pytest.raises(ValueError, match="unique"):
        validate_config(config)


def test_config_rejects_partial_audited_field_set() -> None:
    config = _shipped_config()
    config["audited_fields"] = ["u1"]
    with pytest.raises(ValueError, match="audited_fields must be exactly"):
        validate_config(config)


def test_config_rejects_unsorted_snapshot_labels() -> None:
    config = _shipped_config()
    config["snapshot_labels"] = ["t004", "t000"]
    with pytest.raises(ValueError, match="unique and sorted"):
        validate_config(config)


def test_config_rejects_out_of_range_continuum_search() -> None:
    config = _shipped_config()
    config["continuum_reference"]["search_points_per_axis"] = 3
    with pytest.raises(ValueError, match="search_points_per_axis"):
        validate_config(config)


def test_config_does_not_restate_the_preregistered_threshold() -> None:
    """The module constant is the single source of truth (P0-D)."""

    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert "PREREGISTERED_MIN_POINTS_PER_FRONT" in text
    config = _shipped_config()
    flat = json.dumps(config, sort_keys=True)
    for forbidden in (
        "min_points_per_front",
        "minimum_points_per_front",
        "min_points",
    ):
        assert f'"{forbidden}"' not in flat


def test_run_rejects_an_unknown_input_evidence_directory(scratch_dir: Path) -> None:
    config = _shipped_config()
    config["input_evidence_dirs"] = ["outputs/no_such_bundle_v0"]
    with pytest.raises(ValueError, match="input evidence directory is missing"):
        run(config, scratch_dir / "audit")


def test_cli_help_text_is_printable_on_a_legacy_console() -> None:
    assert audit_module.__doc__ is not None
    audit_module.__doc__.encode("ascii")
    with pytest.raises(SystemExit) as exit_info:
        audit_module.main(["--help"])
    assert exit_info.value.code == 0


# --------------------------------------------------------------------------
# input integrity
# --------------------------------------------------------------------------


def test_verify_evidence_bundle_accepts_the_writer_it_ships_with(
    scratch_dir: Path,
) -> None:
    bundle = _synthetic_bundle(scratch_dir)
    verify_digest(bundle / "manifest.json")
    record = verify_evidence_bundle(bundle)
    assert record["verified"] is True
    assert record["experiment_id"] == "synthetic_core_width_source_v0"
    assert record["checkpoint_count"] == len(SYNTHETIC_RESOLUTIONS) * len(
        SYNTHETIC_TIMES
    )
    assert set(record["checkpoints"]) == {
        (nr, nz, f"t{index:03d}")
        for nr, nz in SYNTHETIC_RESOLUTIONS
        for index in range(len(SYNTHETIC_TIMES))
    }


def test_verify_evidence_bundle_aborts_on_a_corrupted_payload(
    scratch_dir: Path,
) -> None:
    bundle = _synthetic_bundle(scratch_dir)
    target = bundle / "checkpoints" / "checkpoint_nr17_nz32_t000.npz"
    payload = bytearray(target.read_bytes())
    payload[-1] ^= 0xFF
    target.write_bytes(bytes(payload))
    with pytest.raises(ValueError, match="does not match"):
        verify_evidence_bundle(bundle)


def test_verify_evidence_bundle_aborts_on_a_tampered_manifest(
    scratch_dir: Path,
) -> None:
    bundle = _synthetic_bundle(scratch_dir)
    manifest = bundle / "manifest.json"
    manifest.write_bytes(manifest.read_bytes() + b" ")
    with pytest.raises(ValueError, match="manifest integrity check failed"):
        verify_evidence_bundle(bundle)


# --------------------------------------------------------------------------
# smoke run on a synthetic bundle
# --------------------------------------------------------------------------


def test_smoke_run_on_a_synthetic_bundle(scratch_dir: Path) -> None:
    bundle = _synthetic_bundle(scratch_dir)
    config = _synthetic_config(bundle)
    output_dir = scratch_dir / "audit"
    summary = run(config, output_dir)

    # -- the bundle itself ------------------------------------------------
    assert sorted(path.name for path in output_dir.iterdir()) == EXPECTED_ARTIFACTS
    verify_digest(output_dir / "manifest.json")
    verify_digest(output_dir / "summary.json")
    verify_digest(output_dir / "config.snapshot.json")
    manifest = strict_json_loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8"),
        label="core width audit manifest",
    )
    assert manifest["schema"] == MANIFEST_SCHEMA
    assert manifest["experiment_id"] == EXPECTED_EXPERIMENT_ID
    for name, entry in manifest["files"].items():
        target = output_dir / name
        assert entry["sha256"] == sha256_file(target)
        assert entry["bytes"] == target.stat().st_size
    assert summary["reproducibility"]["config_sha256"] == sha256_bytes(
        (output_dir / "config.snapshot.json").read_bytes()
    )

    # -- preregistration is taken from the module, not the config ---------
    assert summary["preregistration"]["min_points_per_front"] == (
        PREREGISTERED_MIN_POINTS_PER_FRONT
    )
    assert summary["amplification_normalization_statement"] == (
        AMPLIFICATION_STATEMENT
    )
    assert "Hou published values are never fit targets" in AMPLIFICATION_STATEMENT
    assert summary["limitations"]
    assert summary["acceptance_checks"] == {
        "all_input_manifests_verified": True,
        "all_requested_snapshots_reported_finite": True,
        "final_snapshot_fit_precondition_recorded_for_every_resolution": True,
        "convergence_fit_precondition_recorded": True,
        "adjacent_resolution_differences_complete": True,
    }
    assert summary["accepted_as_core_width_audit"] is True

    # -- per snapshot reports and the gate --------------------------------
    assert summary["resolutions"] == [
        {"nr": nr, "nz": nz} for nr, nz in SYNTHETIC_RESOLUTIONS
    ]
    assert len(summary["snapshots"]) == len(SYNTHETIC_RESOLUTIONS) * len(
        SYNTHETIC_TIMES
    )
    for snapshot in summary["snapshots"]:
        assert set(snapshot["fields"]) == set(AUDITED_FIELDS)
        for name in AUDITED_FIELDS:
            block = snapshot["fields"][name]
            report = block["report"]
            assert block["report_is_finite"] is True
            for key in (
                "points_per_fwhm_r",
                "points_per_fwhm_z",
                "points_per_front",
            ):
                assert math.isfinite(float(report[key]))
            assert isinstance(report["peak_to_axis_cells"], int)
            assert set(report["truncation_flags"]) == {
                "fwhm_r_left",
                "fwhm_r_right",
                "fwhm_z_left",
                "fwhm_z_right",
                "front_90",
                "front_10",
            }
            assert "high_frequency_tail" in report
            assert "gradient_length_scale" in report
            assert "subgrid_peak_r" in report and "subgrid_peak_z" in report
            gate = block["fit_precondition"]
            assert isinstance(gate["satisfied"], bool)
            assert gate["min_points"] == PREREGISTERED_MIN_POINTS_PER_FRONT
            if not gate["satisfied"]:
                assert gate["reasons"]

    # -- the gate verdict is recorded for every resolution at T_final -----
    final_rows = [
        snapshot
        for snapshot in summary["snapshots"]
        if snapshot["snapshot_label"] == config["final_snapshot_label"]
    ]
    assert sorted((row["nr"], row["nz"]) for row in final_rows) == sorted(
        SYNTHETIC_RESOLUTIONS
    )

    # -- and the aggregate boolean is computed, not asserted --------------
    assert isinstance(
        summary["convergence_fit_precondition_satisfied_anywhere"], bool
    )
    assert summary["convergence_fit_precondition_satisfied_anywhere"] == any(
        ladder["ladder_satisfied"]
        for ladder in summary["fit_precondition_ladders"]
    )
    for ladder in summary["fit_precondition_ladders"]:
        assert ladder["ladder_satisfied"] == (
            len(ladder["verdicts"]) >= 2
            and all(item["satisfied"] for item in ladder["verdicts"])
        )

    # -- both normalizations, side by side and consistent -----------------
    amplification = summary["amplification_renormalization"]
    assert len(amplification) == len(SYNTHETIC_RESOLUTIONS)
    denominators = set()
    for entry in amplification:
        discrete = entry["initial_max_cartesian_vorticity_discrete"]
        reference = entry["common_continuum_reference"]
        final = entry["final_max_cartesian_vorticity"]
        assert discrete > 0.0
        assert reference == pytest.approx(E29B_MAX_CARTESIAN_VORTICITY_EXACT)
        assert entry["amplification_grid_normalized"] == pytest.approx(
            final / discrete
        )
        assert entry["amplification_common_reference"] == pytest.approx(
            final / reference
        )
        assert entry["discrete_initial_over_continuum_reference"] == pytest.approx(
            discrete / reference
        )
        offset = entry["initial_max_location_offset"]
        assert offset["distance"] == pytest.approx(
            math.hypot(offset["delta_r"], offset["delta_z_wrapped"])
        )
        assert abs(offset["delta_z_wrapped"]) <= 0.5 + 1.0e-12
        denominators.add(discrete)
    # The grid-normalized denominator is resolution dependent; the common
    # reference is a single number shared by every row.  That is the whole
    # point of reporting both.
    assert len(denominators) == len(SYNTHETIC_RESOLUTIONS)
    assert len({entry["common_continuum_reference"] for entry in amplification}) == 1

    # -- adjacent-resolution common-grid differences ----------------------
    differences = summary["common_grid_differences"]
    assert len(differences) == (
        (len(SYNTHETIC_RESOLUTIONS) - 1) * len(SYNTHETIC_TIMES) * len(AUDITED_FIELDS)
    )
    for entry in differences:
        assert (entry["coarse_nr"], entry["fine_nr"]) == (
            SYNTHETIC_RESOLUTIONS[0][0],
            SYNTHETIC_RESOLUTIONS[1][0],
        )
        assert entry["common_shape"] == list(SYNTHETIC_RESOLUTIONS[0])
        for key in ("linf", "weighted_l2", "derivative_linf"):
            assert math.isfinite(float(entry[key]))
            assert float(entry[key]) >= 0.0
    assert {entry["field"] for entry in differences} == set(AUDITED_FIELDS)

    # -- the CSV tables mirror the JSON -----------------------------------
    import csv

    def _rows(name: str) -> list[dict[str, str]]:
        with (output_dir / name).open(encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream))

    assert len(_rows("peak_resolution_reports.csv")) == len(
        summary["snapshots"]
    ) * len(AUDITED_FIELDS)
    assert len(_rows("fit_precondition.csv")) == len(summary["snapshots"]) * len(
        AUDITED_FIELDS
    )
    assert len(_rows("common_grid_differences.csv")) == len(differences)
    assert len(_rows("amplification_renormalization.csv")) == len(amplification)


def test_run_refuses_a_nonempty_output_directory(scratch_dir: Path) -> None:
    bundle = _synthetic_bundle(scratch_dir)
    output_dir = scratch_dir / "occupied"
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(_synthetic_config(bundle), output_dir)
    assert (output_dir / "existing.txt").read_text(encoding="utf-8") == "keep me"


def test_prepare_output_refuses_nonempty_and_accepts_empty(
    scratch_dir: Path,
) -> None:
    empty = scratch_dir / "empty"
    _prepare_output(empty)
    assert empty.is_dir()
    _prepare_output(empty)
    (empty / "stray.bin").write_bytes(b"x")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        _prepare_output(empty)
    occupied_file = scratch_dir / "not_a_directory"
    occupied_file.write_text("x", encoding="utf-8")
    with pytest.raises(FileExistsError, match="not a directory"):
        _prepare_output(occupied_file)


def test_run_refuses_to_write_inside_an_audited_source_directory() -> None:
    config = _shipped_config()
    target = REPOSITORY_ROOT / "outputs" / "hou_early_time_v1" / "core_width_audit"
    with pytest.raises(ValueError, match="audited source evidence directory"):
        run(config, target)
    assert not target.exists()


# --------------------------------------------------------------------------
# the continuum reference is computed, not quoted
# --------------------------------------------------------------------------


def test_continuum_initial_maximum_reproduces_the_closed_form() -> None:
    result = continuum_initial_vorticity_maximum(
        search_points_per_axis=801, refinement_stages=3
    )
    assert result["relative_difference_to_closed_form"] < 1.0e-9
    assert result["r"] == pytest.approx(1.0 / math.sqrt(37.0), abs=1.0e-6)
    assert result["z"] == pytest.approx(0.0, abs=1.0e-6)
    assert result["closed_form_value"] == pytest.approx(
        E29B_MAX_CARTESIAN_VORTICITY_EXACT
    )
    assert len(result["stages"]) == 3
    windows = [
        stage["r_window"][1] - stage["r_window"][0] for stage in result["stages"]
    ]
    assert windows[0] > windows[1] > windows[2] > 0.0


# --------------------------------------------------------------------------
# the shipped evidence bundle, when it is present
# --------------------------------------------------------------------------


def test_real_bundle_records_the_257x512_final_fit_precondition_verdict() -> None:
    if not (REAL_OUTPUT_DIR / "summary.json").is_file():
        pytest.skip("outputs/core_width_audit_v1 is not in this checkout")
    verify_digest(REAL_OUTPUT_DIR / "manifest.json")
    verify_digest(REAL_OUTPUT_DIR / "summary.json")
    manifest = strict_json_loads(
        (REAL_OUTPUT_DIR / "manifest.json").read_text(encoding="utf-8"),
        label="core width audit manifest",
    )
    assert manifest["schema"] == MANIFEST_SCHEMA
    for name, entry in manifest["files"].items():
        target = REAL_OUTPUT_DIR / name
        assert entry["sha256"] == sha256_file(target)
        assert entry["bytes"] == target.stat().st_size

    summary = strict_json_loads(
        (REAL_OUTPUT_DIR / "summary.json").read_text(encoding="utf-8"),
        label="core width audit summary",
    )
    assert summary["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert summary["amplification_normalization_statement"] == (
        AMPLIFICATION_STATEMENT
    )
    assert summary["preregistration"]["min_points_per_front"] == (
        PREREGISTERED_MIN_POINTS_PER_FRONT
    )
    assert isinstance(
        summary["convergence_fit_precondition_satisfied_anywhere"], bool
    )
    assert summary["convergence_fit_precondition_satisfied_anywhere"] == any(
        ladder["ladder_satisfied"]
        for ladder in summary["fit_precondition_ladders"]
    )

    final_label = "t004"
    rows = [
        snapshot
        for snapshot in summary["snapshots"]
        if snapshot["nr"] == 257
        and snapshot["nz"] == 512
        and snapshot["snapshot_label"] == final_label
    ]
    assert len(rows) == 1, "the 257x512 T_1 snapshot must be audited exactly once"
    for name in AUDITED_FIELDS:
        gate = rows[0]["fields"][name]["fit_precondition"]
        assert isinstance(gate["satisfied"], bool)
        assert gate["min_points"] == PREREGISTERED_MIN_POINTS_PER_FRONT
        assert math.isfinite(float(gate["points_per_front"]))
        if not gate["satisfied"]:
            assert gate["reasons"]

    # Every audited resolution must carry a verdict at T_1: no silent gaps.
    resolutions = {(entry["nr"], entry["nz"]) for entry in summary["resolutions"]}
    recorded = {
        (snapshot["nr"], snapshot["nz"])
        for snapshot in summary["snapshots"]
        if snapshot["snapshot_label"] == final_label
    }
    assert recorded == resolutions
    assert len(summary["amplification_renormalization"]) == len(resolutions)
