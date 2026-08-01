from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from experiments.run_leray_relay_discovery import (
    CONFIG_SCHEMA,
    _load_config,
    _simple_response_orbit_passes,
    main,
    verify_leray_relay_bundle,
)
from ns_certificate_lab._integrity import canonical_json_bytes, write_with_digest


def _small_config(path: Path) -> None:
    value = {
        "schema": CONFIG_SCHEMA,
        "seed": 20260801,
        "viscosity": "1/40",
        "exact_scale": 1,
        "exact_parent_sine": "1",
        "exact_parent_cosine": "1",
        "exact_child_cosine": "1/8",
        "search_grid_size": 48,
        "recursive_grid_size": 48,
        "carrier": 4,
        "envelope": 2,
        "random_polarization_samples": 1,
        "chirp_values": [0.0],
        "recursive_harmonics": [2],
        "baseline_polarization_real": [0.63377081, -0.12052672, -0.23055982],
        "baseline_polarization_imag": [-0.71460668, -0.13277247, -0.04856754],
    }
    path.write_text(json.dumps(value), encoding="utf-8")


def test_discovery_bundle_verifies_and_tampering_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "config.json"
    output = tmp_path / "bundle"
    _small_config(config)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_leray_relay_discovery",
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ],
    )
    main()
    verdict = verify_leray_relay_bundle(output)
    assert verdict["verified"]
    assert verdict["payload_count"] == 6

    summary = output / "summary.json"
    summary.write_bytes(summary.read_bytes() + b" ")
    with pytest.raises(ValueError, match="checksum"):
        verify_leray_relay_bundle(output)


def test_driver_rejects_nonempty_output_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "config.json"
    output = tmp_path / "occupied"
    output.mkdir()
    (output / "keep.txt").write_text("preserve", encoding="utf-8")
    _small_config(config)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_leray_relay_discovery",
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ],
    )
    with pytest.raises(ValueError, match="new or empty"):
        main()
    assert (output / "keep.txt").read_text(encoding="utf-8") == "preserve"


def test_bundle_rejects_a_self_consistent_but_wrong_source_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "config.json"
    output = tmp_path / "bundle"
    _small_config(config)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_leray_relay_discovery",
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ],
    )
    main()
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    first_source = sorted(manifest["source_files"])[0]
    manifest["source_files"][first_source] = "0" * 64
    write_with_digest(manifest_path, canonical_json_bytes(manifest))
    with pytest.raises(ValueError, match="source checksum mismatch"):
        verify_leray_relay_bundle(output)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("exact_child_cosine", "-1/8", "positive B\\*C\\*D"),
        ("exact_child_cosine", "0", "positive B\\*C\\*D"),
        ("viscosity", "10", "post-viscous"),
    ],
)
def test_config_rejects_a_nonpositive_exact_relay(
    tmp_path: Path, field: str, value: str, message: str
) -> None:
    config = tmp_path / "config.json"
    _small_config(config)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload[field] = value
    config.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        _load_config(config)


def test_config_rejects_zero_polarization_and_unfitted_recursive_grid(
    tmp_path: Path,
) -> None:
    config = tmp_path / "config.json"
    _small_config(config)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["baseline_polarization_real"] = [0.0, 0.0, 0.0]
    payload["baseline_polarization_imag"] = [0.0, 0.0, 0.0]
    config.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="nonzero"):
        _load_config(config)

    _small_config(config)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["recursive_harmonics"] = [2, 4, 8]
    config.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="Nyquist"):
        _load_config(config)


def test_response_orbit_screen_requires_every_full_flux() -> None:
    passing = [
        {
            "normalized_injection": 1.0,
            "off_chain_forcing_ratio": 0.5,
            "full_populated_flux_tested": True,
            "full_populated_child_flux": 0.25,
        },
        {
            "normalized_injection": 0.6,
            "off_chain_forcing_ratio": 0.75,
            "full_populated_flux_tested": True,
            "full_populated_child_flux": 0.1,
        },
    ]
    assert _simple_response_orbit_passes(passing)

    untested = [dict(row) for row in passing]
    untested[-1]["full_populated_flux_tested"] = False
    untested[-1]["full_populated_child_flux"] = None
    assert not _simple_response_orbit_passes(untested)

    nonpositive = [dict(row) for row in passing]
    nonpositive[-1]["full_populated_child_flux"] = 0.0
    assert not _simple_response_orbit_passes(nonpositive)
