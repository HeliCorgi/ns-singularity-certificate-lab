"""End-to-end integrity checks for the Zeno algebraic pilot bundle."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest

from experiments.run_zeno_packet_relay_pilot import (
    main,
    verify_zeno_packet_relay_bundle,
)


ROOT = Path(__file__).resolve().parents[1]


def _small_config(path: Path) -> None:
    config = json.loads(
        (ROOT / "configs" / "zeno_packet_relay_pilot_v2.json").read_text(
            encoding="utf-8"
        )
    )
    config["levels"] = [2, 4]
    config["random_orientation_samples"] = 8
    path.write_text(json.dumps(config), encoding="utf-8")


def test_v2_bundle_verifies_and_tampering_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "config.json"
    output = tmp_path / "bundle"
    _small_config(config)
    monkeypatch.chdir(ROOT)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_zeno_packet_relay_pilot",
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ],
    )
    main()

    verdict = verify_zeno_packet_relay_bundle(output)
    assert verdict == {
        "verified": True,
        "schema": "ns-certificate-lab/zeno-packet-relay-pilot/v2",
        "payload_count": 5,
        "source_file_count": 6,
    }
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["pressure_checks"]["moment_eigenvalues"] == pytest.approx(
        [0.0, 0.5, 0.5]
    )
    assert not summary["pressure_checks"][
        "exact_pressure_and_multipole_remainder_tested"
    ]
    assert "source_fingerprint_sha256" in summary["provenance"]

    (output / "summary.json").write_bytes(
        (output / "summary.json").read_bytes() + b"\n"
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_zeno_packet_relay_bundle(output)


def test_driver_rejects_nonempty_output_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "config.json"
    output = tmp_path / "bundle"
    output.mkdir()
    (output / "sentinel").write_text("preserve", encoding="utf-8")
    _small_config(config)
    monkeypatch.chdir(ROOT)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_zeno_packet_relay_pilot",
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ],
    )
    with pytest.raises(ValueError, match="new or empty"):
        main()
    assert (output / "sentinel").read_text(encoding="utf-8") == "preserve"
