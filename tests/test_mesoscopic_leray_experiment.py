from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import shutil
import sys

import pytest

from experiments.run_mesoscopic_leray_cloud import (
    CONFIG_SCHEMA,
    PAYLOAD_NAMES,
    _build_scaling_fits,
    _load_config,
    main,
    verify_mesoscopic_leray_bundle,
)
from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_file,
    write_with_digest,
)


def _small_config(path: Path) -> None:
    payload = {
        "schema": CONFIG_SCHEMA,
        "seed": 20260801,
        "viscosity": "1/40",
        "energy_constant": 1.0,
        "tau": 0.25,
        "scales": [8, 12, 16, 24],
        "gammas": [0.4],
        "relative_width_fractions": [0.1],
        "exact_target_limit": 0,
        "strata_per_axis": 1,
        "exact_energy_pair_limit": 0,
        "local_fft_maximum_working_bytes": 200000000,
        "required_duhamel_ratio": 0.5,
        "small_galerkin": {
            "scale": 2,
            "width": 1,
            "grid_size": 32,
            "tau": 0.1,
            "steps": 2,
        },
        "carrier_two_stage": {
            "scale": 2,
            "grid_size": 32,
            "tau": 0.1,
            "time_multiples": 2.0,
            "steps": 2,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


@pytest.fixture(scope="module")
def generated_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("mesoscopic_bundle")
    config = root / "config.json"
    output = root / "bundle"
    _small_config(config)
    old_argv = sys.argv
    try:
        sys.argv = [
            "run_mesoscopic_leray_cloud",
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ]
        main()
    finally:
        sys.argv = old_argv
    return output


def _copy_bundle(source: Path, destination: Path) -> Path:
    shutil.copytree(source, destination)
    return destination


def _resign_payload(bundle: Path, name: str, payload: bytes) -> None:
    path = bundle / name
    write_with_digest(path, payload)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"][name] = {
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    write_with_digest(manifest_path, canonical_json_bytes(manifest))


def test_bundle_verifies_and_records_both_galerkin_diagnostics(
    generated_bundle: Path,
) -> None:
    verdict = verify_mesoscopic_leray_bundle(generated_bundle)
    assert verdict["verified"]
    assert verdict["payload_count"] == len(PAYLOAD_NAMES) == 8
    assert verdict["power_row_count"] == 4
    assert verdict["relative_row_count"] == 4

    summary = json.loads(
        (generated_bundle / "summary.json").read_text(encoding="utf-8")
    )
    assert summary["sublinear_verdict"]["classification"].startswith("REJECTED")
    assert summary["two_stage_comparison"]["joint_strict_finite_search"] == {
        "status": "INFEASIBLE_IN_STATED_FINITE_ALPHABET",
        "exhaustive": True,
        "orientations_tested": 16,
        "found": False,
        "joint_objective_j_n": None,
        "reason": "no admissible graph, so no optimization score is manufactured",
    }
    assert summary["two_stage_galerkin"]["pathway_contaminated"] is True
    assert summary["two_stage_galerkin"]["grandchild_binary64_resolved"] is True

    comparison = json.loads(
        (generated_bundle / "galerkin_comparison.json").read_text(
            encoding="utf-8"
        )
    )
    assert comparison["one_stage_mesoscopic"]["initial_child_energy"] == 0.0
    assert "pathway-contaminated" in comparison["two_stage_carrier"][
        "interpretation"
    ]
    svg = (generated_bundle / "shell_energy.svg").read_text(encoding="utf-8")
    assert 'id="cross-talk-energy"' in svg
    assert 'id="intended-grandchild-energy"' in svg

    fits = json.loads(
        (generated_bundle / "scaling_fits.json").read_text(encoding="utf-8")
    )
    assert fits["groups"][0]["fits"]["gain_g"]["input_scales"] == [
        8,
        12,
        16,
        24,
    ]
    assert fits["groups"][0]["fits"]["gain_g"]["eligibility_class"] == (
        "diagnostic_core_geometry_not_relay_acceptance"
    )
    assert fits["groups"][0]["fits"]["gain_g"]["status"] == (
        "fit_last_four_eligible_points"
    )


def test_checksum_tampering_is_rejected(
    generated_bundle: Path, tmp_path: Path
) -> None:
    bundle = _copy_bundle(generated_bundle, tmp_path / "checksum_tamper")
    summary = bundle / "summary.json"
    summary.write_bytes(summary.read_bytes() + b" ")
    with pytest.raises(ValueError, match="checksum"):
        verify_mesoscopic_leray_bundle(bundle)


def test_resigned_sublinear_acceptance_is_semantically_rejected(
    generated_bundle: Path, tmp_path: Path
) -> None:
    bundle = _copy_bundle(generated_bundle, tmp_path / "semantic_tamper")
    path = bundle / "mesoscopic_scaling.csv"
    reader = csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")))
    rows = list(reader)
    rows[0]["asymptotic_screen_pass"] = "True"
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=reader.fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    _resign_payload(bundle, path.name, stream.getvalue().encode("utf-8"))
    with pytest.raises(ValueError, match="sublinear gamma"):
        verify_mesoscopic_leray_bundle(bundle)


@pytest.mark.parametrize(
    ("field", "replacement", "message"),
    (
        ("critical_energy", "999", "critical energy"),
        ("parent_norm_squared", "999", "parent norm"),
        ("effective_mode_count", "999999999", "effective mode count"),
        ("heat_factor_h", "999", "heat factor"),
        ("child_forcing_fraction", "0.123", "child forcing fraction"),
        ("energy_cancellation_pairing", "1", "cancellation residual"),
        ("parent_divergence_relative", "0.001", "legacy structural residual"),
        ("nonlinear_reality_relative", "0.001", "structural residual"),
    ),
)
def test_resigned_independent_row_invariant_tampering_is_rejected(
    generated_bundle: Path,
    tmp_path: Path,
    field: str,
    replacement: str,
    message: str,
) -> None:
    bundle = _copy_bundle(generated_bundle, tmp_path / f"row_tamper_{field}")
    path = bundle / "mesoscopic_scaling.csv"
    reader = csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")))
    rows = list(reader)
    rows[0][field] = replacement
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=reader.fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    _resign_payload(bundle, path.name, stream.getvalue().encode("utf-8"))
    with pytest.raises(ValueError, match=message):
        verify_mesoscopic_leray_bundle(bundle)


def test_resigned_wrong_current_source_hash_is_rejected(
    generated_bundle: Path, tmp_path: Path
) -> None:
    bundle = _copy_bundle(generated_bundle, tmp_path / "source_tamper")
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    first_source = sorted(manifest["source_files"])[0]
    manifest["source_files"][first_source] = "0" * 64
    write_with_digest(manifest_path, canonical_json_bytes(manifest))
    with pytest.raises(ValueError, match="source checksum mismatch"):
        verify_mesoscopic_leray_bundle(bundle)


def test_driver_refuses_to_overwrite_existing_evidence(
    tmp_path: Path,
) -> None:
    config = tmp_path / "config.json"
    output = tmp_path / "occupied"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("preserve", encoding="utf-8")
    _small_config(config)
    old_argv = sys.argv
    try:
        sys.argv = [
            "run_mesoscopic_leray_cloud",
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ]
        with pytest.raises(ValueError, match="new or empty"):
            main()
    finally:
        sys.argv = old_argv
    assert marker.read_text(encoding="utf-8") == "preserve"


def test_config_requires_registered_energy_limit_and_two_stage_block(
    tmp_path: Path,
) -> None:
    config = tmp_path / "config.json"
    _small_config(config)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload.pop("exact_energy_pair_limit")
    config.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="missing or unknown"):
        _load_config(config)

    _small_config(config)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["carrier_two_stage"]["time_multiples"] = 0.0
    config.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="time_multiples"):
        _load_config(config)


def test_core_gain_fit_is_diagnostic_when_full_relay_geometry_fails() -> None:
    rows = []
    for scale in (16, 32, 48, 64):
        rows.append(
            {
                "family": "power_width",
                "family_parameter": 0.7,
                "gamma": 0.7,
                "base_scale": scale,
                "core_fit_eligible": True,
                "relay_fit_eligible": False,
                "full_output_bands_disjoint": False,
                "forcing_metric_method": "overlapping-band-projections-only",
                "support_mode_count": float(scale**2),
                "a_unit": float(scale**2),
                "a_critical": float(scale),
                "gain_g": float(scale) ** 0.2,
                "duhamel_energy_ratio": float(scale) ** -0.9,
                "off_chain_main_ratio": 1.0,
                "outside_child_full_main_ratio": 1.0,
                "self_interaction_cross_interaction_ratio": 1.0,
            }
        )
    fits = _build_scaling_fits(rows)["groups"][0]["fits"]
    assert fits["gain_g"]["status"] == "fit_last_four_eligible_points"
    assert fits["gain_g"]["exponent"] == pytest.approx(0.2)
    assert fits["gain_g"]["eligibility_class"] == (
        "diagnostic_core_geometry_not_relay_acceptance"
    )
    assert fits["outside_child_full_main_ratio"]["status"] == (
        "insufficient_eligible_points"
    )
    assert fits["self_interaction_cross_interaction_ratio"]["status"] == (
        "insufficient_eligible_points"
    )

    for row in rows:
        row["forcing_metric_method"] = "exact-sparse-full"
    exact_full_fits = _build_scaling_fits(rows)["groups"][0]["fits"]
    assert exact_full_fits["outside_child_full_main_ratio"]["status"] == (
        "fit_last_four_eligible_points"
    )
    assert exact_full_fits["self_interaction_cross_interaction_ratio"][
        "status"
    ] == "fit_last_four_eligible_points"
