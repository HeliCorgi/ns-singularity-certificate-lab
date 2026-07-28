"""Tests for the frozen-coefficient von Neumann audit experiment.

The expensive path (every recorded row of every shipped Heun bundle) is run by
the orchestrator, not by this module.  What is exercised here is the config
contract, the exclusion and verdict arithmetic on a SYNTHETIC bundle whose two
auditable rows are constructed to fall on opposite sides of the Heun
frozen-coefficient threshold, the refusal to overwrite evidence, and -- if it
has been produced -- the integrity of the real bundle.

Synthetic operating points, and why their verdicts are not guessed
-----------------------------------------------------------------
With ``viscosity = 0`` the symbol is purely imaginary, ``a = dt*(|c_r|/dr +
|c_z|/dz)`` at the worst wavenumber, and Heun gives ``|G|^2 = 1 + a^4/4``:

* ``dt = 1e-12`` on the 9x16 synthetic grid gives ``a = 2.4e-11`` and
  ``a^4/4 = 8e-44``, which is below the binary64 resolution of 1, so the scan
  returns exactly 1.0 and the row passes at tolerance 1e-12;
* ``dt = 0.05`` gives ``a = 1.2``, ``|G| = 1.232`` and a predictor stage
  ``|1 + i a| = 1.562``, so the row fails by a wide margin.

Both numbers are recomputed inside the tests from the same closed forms rather
than pasted, so a change to the stability polynomial cannot silently pass.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import pytest

import experiments.run_von_neumann_audit as audit_module
from experiments.run_von_neumann_audit import (
    AUDIT_N_THETA,
    EXPECTED_EXPERIMENT_ID,
    INTERPRETATION,
    SHIPPED_REFERENCE_POINTS,
    UNVERIFIED_PHRASE,
    VERDICT_UNVERIFIED,
    VERDICT_VERIFIED,
    _prepare_output,
    evaluate,
    read_source_rows,
    run,
    validate_config,
    verify_bundle,
)
from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_file,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "von_neumann_audit.json"
REAL_OUTPUT_DIR = REPOSITORY_ROOT / "outputs" / "von_neumann_audit_v1"

SYNTHETIC_NR = 9
SYNTHETIC_NZ = 16
SYNTHETIC_DR = 1.0 / (SYNTHETIC_NR - 1)
SYNTHETIC_DZ = 1.0 / SYNTHETIC_NZ
STABLE_DT = 1.0e-12
FAILING_DT = 0.05
SYNTHETIC_SPEED = 1.0

DIAGNOSTIC_COLUMNS = (
    "nr",
    "nz",
    "step",
    "time",
    "dt",
    "viscosity",
    "max_abs_u_r",
    "max_abs_u_z",
)


def _inviscid_heun_amplification(dt: float) -> float:
    """``|G| = (1 + a^4/4)^(1/2)`` at the worst wavenumber, ``nu = 0``."""

    a = dt * (SYNTHETIC_SPEED / SYNTHETIC_DR + SYNTHETIC_SPEED / SYNTHETIC_DZ)
    return math.sqrt(1.0 + a**4 / 4.0)


def _inviscid_euler_amplification(dt: float) -> float:
    """``|1 + i a|`` at the worst wavenumber, ``nu = 0``."""

    a = dt * (SYNTHETIC_SPEED / SYNTHETIC_DR + SYNTHETIC_SPEED / SYNTHETIC_DZ)
    return math.hypot(1.0, a)


def _synthetic_rows(*, include_failing_row: bool = True) -> list[dict[str, Any]]:
    """A dt == 0 record plus one clearly passing and one clearly failing row."""

    rows = [
        {
            "nr": SYNTHETIC_NR,
            "nz": SYNTHETIC_NZ,
            "step": 0,
            "time": 0.0,
            "dt": 0.0,
            "viscosity": 0.0,
            "max_abs_u_r": 0.0,
            "max_abs_u_z": 0.0,
        },
        {
            "nr": SYNTHETIC_NR,
            "nz": SYNTHETIC_NZ,
            "step": 25,
            "time": 2.5e-11,
            "dt": STABLE_DT,
            "viscosity": 0.0,
            "max_abs_u_r": SYNTHETIC_SPEED,
            "max_abs_u_z": SYNTHETIC_SPEED,
        },
    ]
    if include_failing_row:
        rows.append(
            {
                "nr": SYNTHETIC_NR,
                "nz": SYNTHETIC_NZ,
                "step": 50,
                "time": 1.25,
                "dt": FAILING_DT,
                "viscosity": 0.0,
                "max_abs_u_r": SYNTHETIC_SPEED,
                "max_abs_u_z": SYNTHETIC_SPEED,
            }
        )
    return rows


def _write_synthetic_bundle(
    directory: Path,
    rows: list[dict[str, Any]],
    *,
    experiment_id: str = "synthetic_bundle_v1",
) -> Path:
    """Write a diagnostics.csv plus the manifest the experiments would write."""

    directory.mkdir(parents=True, exist_ok=True)
    diagnostics = directory / "diagnostics.csv"
    with diagnostics.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(DIAGNOSTIC_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": 1,
        "experiment_id": experiment_id,
        "files": {
            "diagnostics.csv": {
                "sha256": sha256_file(diagnostics),
                "bytes": diagnostics.stat().st_size,
            }
        },
    }
    write_with_digest(directory / "manifest.json", canonical_json_bytes(manifest))
    return directory


def _relative(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY_ROOT).as_posix()


def _synthetic_config(
    bundles: list[tuple[str, Path]],
    *,
    stride: int = 25,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "experiment_id": EXPECTED_EXPERIMENT_ID,
        "description": "synthetic smoke bundle for the von Neumann audit",
        "interpretation": (
            "A failing row is " + UNVERIFIED_PHRASE + "; this synthetic bundle "
            "exists only to exercise the audit machinery."
        ),
        "sources": [
            {
                "label": label,
                "bundle": _relative(path),
                "stride": stride,
                "r_max": 1.0,
                "z_period": 1.0,
                "grid": None,
                "extra_run_key_columns": [],
            }
            for label, path in bundles
        ],
        "methods": ["heun"],
        "tolerance": 1.0e-12,
        "n_theta": AUDIT_N_THETA,
    }


def _shipped_config() -> dict[str, Any]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="von neumann audit config",
    )
    assert isinstance(value, dict)
    return value


# --------------------------------------------------------------------------
# configuration contract
# --------------------------------------------------------------------------


def test_shipped_config_is_valid_and_names_the_real_bundles() -> None:
    config = _shipped_config()
    validate_config(config)
    assert config["experiment_id"] == EXPECTED_EXPERIMENT_ID
    labels = [source["label"] for source in config["sources"]]
    assert labels == [
        "hou_early_time_v1",
        "hou_early_time_v2_hires",
        "hou_time_refinement_v1",
    ]
    for source in config["sources"]:
        bundle = REPOSITORY_ROOT / source["bundle"]
        assert (bundle / "diagnostics.csv").is_file(), source["bundle"]
        assert (bundle / "manifest.json").is_file(), source["bundle"]
        # The stride limitation is a fact about the shipped bundles, not a
        # tuning knob: every audited run recorded one row every 25 steps.
        snapshot = strict_json_loads(
            (bundle / "config.snapshot.json").read_text(encoding="utf-8"),
            label="source config snapshot",
        )
        assert int(snapshot["diagnostic_stride"]) == int(source["stride"])


@pytest.mark.parametrize(
    "key",
    (
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "sources",
        "methods",
        "tolerance",
        "n_theta",
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


def test_config_rejects_unsupported_schema_version() -> None:
    config = _shipped_config()
    config["schema_version"] = 2
    with pytest.raises(ValueError, match="schema_version"):
        validate_config(config)


def test_config_rejects_interpretation_without_the_required_phrase() -> None:
    config = _shipped_config()
    config["interpretation"] = "A failing row means the run was unstable."
    with pytest.raises(ValueError, match="exact phrase"):
        validate_config(config)


def test_config_rejects_empty_sources() -> None:
    config = _shipped_config()
    config["sources"] = []
    with pytest.raises(ValueError, match="nonempty list"):
        validate_config(config)


@pytest.mark.parametrize(
    "bundle",
    ("C:/absolute/bundle", "/absolute/bundle", "../outside/bundle"),
)
def test_config_rejects_escaping_bundle_path(bundle: str) -> None:
    config = _shipped_config()
    config["sources"][0]["bundle"] = bundle
    with pytest.raises(ValueError, match="repository-relative"):
        validate_config(config)


def test_config_rejects_duplicate_source_labels() -> None:
    config = _shipped_config()
    config["sources"][1]["label"] = config["sources"][0]["label"]
    with pytest.raises(ValueError, match="labels must be unique"):
        validate_config(config)


def test_config_rejects_unknown_source_key() -> None:
    config = _shipped_config()
    config["sources"][0]["unexpected"] = 1
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(config)


@pytest.mark.parametrize("stride", (0, -1, 2.5, True))
def test_config_rejects_bad_stride(stride: Any) -> None:
    config = _shipped_config()
    config["sources"][0]["stride"] = stride
    with pytest.raises(ValueError, match="stride"):
        validate_config(config)


@pytest.mark.parametrize("grid", ([65], [65, 128, 256], [2, 128], "65x128"))
def test_config_rejects_bad_grid(grid: Any) -> None:
    config = _shipped_config()
    config["sources"][2]["grid"] = grid
    with pytest.raises(ValueError, match="grid must be null"):
        validate_config(config)


@pytest.mark.parametrize("methods", (["ssprk3"], ["heun", "rk4"], [], "heun"))
def test_config_rejects_methods_other_than_heun_alone(methods: Any) -> None:
    config = _shipped_config()
    config["methods"] = methods
    with pytest.raises(ValueError, match="methods must be exactly"):
        validate_config(config)


@pytest.mark.parametrize("tolerance", (-1.0e-12, "1e-12", True))
def test_config_rejects_bad_tolerance(tolerance: Any) -> None:
    config = _shipped_config()
    config["tolerance"] = tolerance
    with pytest.raises(ValueError, match="tolerance"):
        validate_config(config)


@pytest.mark.parametrize("n_theta", (181, 720, 4097))
def test_config_rejects_an_n_theta_that_audit_snapshot_would_not_scan(
    n_theta: int,
) -> None:
    config = _shipped_config()
    config["n_theta"] = n_theta
    with pytest.raises(ValueError, match="scan_amplification default"):
        validate_config(config)


def test_config_rejects_duplicate_run_key_columns() -> None:
    config = _shipped_config()
    config["sources"][2]["extra_run_key_columns"] = ["level_index", "level_index"]
    with pytest.raises(ValueError, match="repeat a column"):
        validate_config(config)


def test_module_interpretation_carries_the_preregistered_phrase() -> None:
    assert UNVERIFIED_PHRASE == "stability-unverified, not unstable"
    assert UNVERIFIED_PHRASE in INTERPRETATION
    # The verdict vocabulary must never collapse to "stable"/"unstable".
    assert VERDICT_UNVERIFIED == "stability-unverified"
    assert VERDICT_VERIFIED == "stability-verified-at-recorded-rows"
    assert VERDICT_VERIFIED.endswith("at-recorded-rows")
    for verdict in (VERDICT_VERIFIED, VERDICT_UNVERIFIED):
        assert verdict not in {"stable", "unstable"}
    assert SHIPPED_REFERENCE_POINTS["radial_binding"]["passes"] is False
    assert SHIPPED_REFERENCE_POINTS["axial_binding"]["passes"] is True
    assert (
        SHIPPED_REFERENCE_POINTS["radial_binding"][
            "heun_full_step_max_amplification"
        ]
        == 1.000003460085721
    )
    assert (
        SHIPPED_REFERENCE_POINTS["axial_binding"]["heun_full_step_max_amplification"]
        == 1.0
    )
    assert "numerical_stability_audit.md" in SHIPPED_REFERENCE_POINTS["source"]


def test_cli_help_text_is_printable_on_a_legacy_console() -> None:
    """argparse prints the module docstring; a non-ASCII glyph crashes cp932."""

    assert audit_module.__doc__ is not None
    audit_module.__doc__.encode("ascii")
    with pytest.raises(SystemExit) as exit_info:
        audit_module.main(["--help"])
    assert exit_info.value.code == 0


# --------------------------------------------------------------------------
# integrity of the audited inputs
# --------------------------------------------------------------------------


def test_verify_bundle_rejects_a_tampered_diagnostics_file(scratch_dir: Path) -> None:
    bundle = _write_synthetic_bundle(scratch_dir / "bundle", _synthetic_rows())
    verify_bundle(bundle)
    diagnostics = bundle / "diagnostics.csv"
    diagnostics.write_text(
        diagnostics.read_text(encoding="utf-8") + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_bundle(bundle)


def test_verify_bundle_rejects_a_missing_manifest_sidecar(scratch_dir: Path) -> None:
    bundle = _write_synthetic_bundle(scratch_dir / "bundle", _synthetic_rows())
    (bundle / "manifest.json.sha256").unlink()
    with pytest.raises(ValueError, match="missing data or checksum sidecar"):
        verify_bundle(bundle)


def test_read_source_rows_counts_the_dt_zero_exclusion(scratch_dir: Path) -> None:
    bundle = _write_synthetic_bundle(scratch_dir / "bundle", _synthetic_rows())
    config = _synthetic_config([("synthetic", bundle)])
    rows, excluded = read_source_rows(config["sources"][0])
    assert excluded == 1
    assert len(rows) == 3
    assert [row["excluded"] for row in rows] == [True, False, False]
    auditable = [row for row in rows if not row["excluded"]]
    assert auditable[0]["dr"] == pytest.approx(SYNTHETIC_DR)
    assert auditable[0]["dz"] == pytest.approx(SYNTHETIC_DZ)


# --------------------------------------------------------------------------
# smoke run on the synthetic bundle
# --------------------------------------------------------------------------


def test_smoke_run_on_synthetic_bundle(scratch_dir: Path) -> None:
    bundle = _write_synthetic_bundle(scratch_dir / "bundle", _synthetic_rows())
    config = _synthetic_config([("synthetic", bundle)])
    output_dir = scratch_dir / "audit"
    summary = run(config, output_dir)

    assert summary["schema_version"] == 1
    assert summary["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert UNVERIFIED_PHRASE in summary["interpretation"]
    assert summary["accepted_as_von_neumann_audit"] is True
    assert summary["acceptance_checks"] == {
        "all_input_manifests_verified": True,
        "all_audited_rows_finite": True,
        "per_run_verdict_recorded": True,
        "summary_reports_verified_and_unverified_booleans": True,
        "interpretation_states_unverified_not_unstable": True,
    }

    # Acceptance check 4: the booleans exist and are computed from the data.
    assert summary["any_run_stability_verified_at_recorded_rows"] is False
    assert summary["all_runs_stability_unverified"] is True

    totals = summary["totals"]
    assert totals["runs"] == 1
    assert totals["rows_audited"] == 2
    assert totals["rows_excluded_dt_zero"] == 1
    assert totals["failing_rows"] == 1

    (source,) = summary["sources"]
    assert source["integrity"]["manifest_digest_verified"] is True
    assert source["integrity"]["diagnostics_matches_manifest"] is True
    assert source["integrity"]["diagnostics_sha256"] == sha256_file(
        bundle / "diagnostics.csv"
    )
    assert source["rows_excluded_dt_zero"] == 1

    (record,) = summary["runs"]
    assert record["nr"] == SYNTHETIC_NR
    assert record["nz"] == SYNTHETIC_NZ
    assert record["rows_audited"] == 2
    assert record["rows_excluded_dt_zero"] == 1
    assert record["failing_rows"] == 1
    assert record["passing_rows"] == 1
    assert record["verdict"] == VERDICT_UNVERIFIED
    assert record["min_dt"] == pytest.approx(STABLE_DT)
    assert record["max_dt"] == pytest.approx(FAILING_DT)

    # The worst row is the constructed failing one, and its two amplification
    # numbers agree with the closed forms recomputed here.
    worst = record["worst_row"]
    assert worst["dt"] == pytest.approx(FAILING_DT)
    assert worst["heun_passes"] is False
    assert worst["heun_max_amplification"] == pytest.approx(
        _inviscid_heun_amplification(FAILING_DT), rel=1.0e-9
    )
    assert record["worst_predictor_row"]["predictor_stage"] == "euler"
    assert record["worst_predictor_row"][
        "predictor_max_amplification"
    ] == pytest.approx(_inviscid_euler_amplification(FAILING_DT), rel=1.0e-9)

    # The strided extrapolation is labelled, reported and reproducible.
    assert "not a bound" in record["strided_growth_extrapolation_note"]
    expected_log_sum = 25.0 * (
        math.log(_inviscid_heun_amplification(STABLE_DT))
        + math.log(_inviscid_heun_amplification(FAILING_DT))
    )
    assert record["strided_log_growth_sum_not_a_bound"] == pytest.approx(
        expected_log_sum, rel=1.0e-9
    )
    assert record["strided_growth_extrapolation_not_a_bound"] == pytest.approx(
        math.exp(expected_log_sum), rel=1.0e-9
    )

    # Per-row verdicts split exactly as constructed.
    with (output_dir / "row_audits.csv").open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2
    stable_row, failing_row = rows
    assert float(stable_row["dt"]) == pytest.approx(STABLE_DT)
    assert stable_row["heun_passes"] == "True"
    assert stable_row["row_verdict"] == "stability-verified-in-model"
    assert float(stable_row["heun_max_amplification"]) == pytest.approx(1.0)
    assert failing_row["heun_passes"] == "False"
    assert failing_row["row_verdict"] == VERDICT_UNVERIFIED
    assert float(failing_row["heun_max_amplification"]) > 1.2
    assert float(failing_row["predictor_max_amplification"]) > 1.5

    # The bundle is written with the shipped conventions and verifies.
    names = sorted(path.name for path in output_dir.iterdir())
    assert names == [
        "config.snapshot.json",
        "config.snapshot.json.sha256",
        "manifest.json",
        "manifest.json.sha256",
        "row_audits.csv",
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
    assert manifest["audited_inputs"][_relative(bundle)]["diagnostics_sha256"] == (
        sha256_file(bundle / "diagnostics.csv")
    )
    assert summary["limitations"]
    assert any("STRIDE" in item for item in summary["limitations"])
    assert summary["runtime"]["scans_performed"] == 4


def test_run_booleans_are_computed_from_the_data(scratch_dir: Path) -> None:
    """A passing bundle and a failing bundle must move both booleans."""

    passing = _write_synthetic_bundle(
        scratch_dir / "passing",
        _synthetic_rows(include_failing_row=False),
        experiment_id="synthetic_passing_v1",
    )
    failing = _write_synthetic_bundle(
        scratch_dir / "failing",
        _synthetic_rows(),
        experiment_id="synthetic_failing_v1",
    )
    config = _synthetic_config([("pass", passing), ("fail", failing)])
    _, runs, _, checks = evaluate(config)
    verdicts = {record["run_id"]: record["verdict"] for record in runs}
    assert sorted(verdicts.values()) == [VERDICT_UNVERIFIED, VERDICT_VERIFIED]
    assert checks["per_run_verdict_recorded"] is True

    summary = run(config, scratch_dir / "audit")
    assert summary["any_run_stability_verified_at_recorded_rows"] is True
    assert summary["all_runs_stability_unverified"] is False
    assert summary["totals"]["runs_stability_verified_at_recorded_rows"] == 1
    assert summary["totals"]["runs_stability_unverified"] == 1


def test_prepare_output_refuses_a_nonempty_directory(scratch_dir: Path) -> None:
    output_dir = scratch_dir / "occupied"
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        _prepare_output(output_dir)
    assert (output_dir / "existing.txt").read_text(encoding="utf-8") == "keep me"

    fresh = scratch_dir / "fresh"
    _prepare_output(fresh)
    assert fresh.is_dir()
    _prepare_output(fresh)  # empty is fine


def test_run_refuses_a_nonempty_output_directory(scratch_dir: Path) -> None:
    bundle = _write_synthetic_bundle(scratch_dir / "bundle", _synthetic_rows())
    output_dir = scratch_dir / "occupied"
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(_synthetic_config([("synthetic", bundle)]), output_dir)
    assert (output_dir / "existing.txt").read_text(encoding="utf-8") == "keep me"


def test_run_refuses_to_write_inside_an_audited_bundle(scratch_dir: Path) -> None:
    bundle = _write_synthetic_bundle(scratch_dir / "bundle", _synthetic_rows())
    target = bundle / "von_neumann"
    with pytest.raises(ValueError, match="audited source evidence directory"):
        run(_synthetic_config([("synthetic", bundle)]), target)
    assert not target.exists()


# --------------------------------------------------------------------------
# the real bundle, if it has been produced
# --------------------------------------------------------------------------


def test_real_bundle_verifies_and_records_every_per_run_verdict() -> None:
    if not (REAL_OUTPUT_DIR / "summary.json").is_file():
        pytest.skip("outputs/von_neumann_audit_v1 has not been produced")

    verify_digest(REAL_OUTPUT_DIR / "manifest.json")
    verify_digest(REAL_OUTPUT_DIR / "summary.json")
    verify_digest(REAL_OUTPUT_DIR / "config.snapshot.json")
    manifest = strict_json_loads(
        (REAL_OUTPUT_DIR / "manifest.json").read_text(encoding="utf-8"),
        label="von neumann audit manifest",
    )
    for name, entry in manifest["files"].items():
        target = REAL_OUTPUT_DIR / name
        assert entry["sha256"] == sha256_file(target), name
        assert entry["bytes"] == target.stat().st_size, name

    summary = strict_json_loads(
        (REAL_OUTPUT_DIR / "summary.json").read_text(encoding="utf-8"),
        label="von neumann audit summary",
    )
    assert summary["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert UNVERIFIED_PHRASE in summary["interpretation"]
    assert isinstance(summary["any_run_stability_verified_at_recorded_rows"], bool)
    assert isinstance(summary["all_runs_stability_unverified"], bool)
    assert summary["accepted_as_von_neumann_audit"] is True

    required_run_fields = {
        "run_id",
        "bundle",
        "nr",
        "nz",
        "rows_audited",
        "rows_excluded_dt_zero",
        "failing_rows",
        "passing_rows",
        "min_dt",
        "max_dt",
        "worst_row",
        "worst_predictor_row",
        "verdict",
        "strided_growth_extrapolation_not_a_bound",
        "strided_growth_extrapolation_note",
    }
    assert summary["runs"]
    for record in summary["runs"]:
        assert required_run_fields <= set(record), record.get("run_id")
        assert record["verdict"] in (VERDICT_VERIFIED, VERDICT_UNVERIFIED)
        assert math.isfinite(record["worst_row"]["heun_max_amplification"])
        assert math.isfinite(
            record["worst_predictor_row"]["predictor_max_amplification"]
        )

    # The booleans in the shipped summary are the ones the verdicts imply.
    verdicts = [record["verdict"] for record in summary["runs"]]
    assert summary["any_run_stability_verified_at_recorded_rows"] == (
        VERDICT_VERIFIED in verdicts
    )
    assert summary["all_runs_stability_unverified"] == all(
        verdict == VERDICT_UNVERIFIED for verdict in verdicts
    )
    assert summary["totals"]["rows_audited"] == sum(
        record["rows_audited"] for record in summary["runs"]
    )

    with (REAL_OUTPUT_DIR / "row_audits.csv").open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == summary["totals"]["rows_audited"]
    for row in rows:
        assert math.isfinite(float(row["heun_max_amplification"]))
        assert math.isfinite(float(row["predictor_max_amplification"]))
        assert row["row_verdict"] in (
            "stability-verified-in-model",
            VERDICT_UNVERIFIED,
        )
