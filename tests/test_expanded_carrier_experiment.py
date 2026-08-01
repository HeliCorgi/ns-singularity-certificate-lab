from __future__ import annotations

import json
from pathlib import Path
import shutil

import pytest

from experiments.run_expanded_carrier_search import (
    PAYLOAD_NAMES,
    run,
    verify_expanded_carrier_bundle,
)
from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_file,
    write_with_digest,
)


@pytest.fixture(scope="module")
def expanded_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    output = tmp_path_factory.mktemp("expanded_carrier") / "bundle"
    run(Path("configs/expanded_carrier_search_v1.json"), output)
    return output


def test_expanded_bundle_recomputes_search_and_independent_verifier(
    expanded_bundle: Path,
) -> None:
    verdict = verify_expanded_carrier_bundle(expanded_bundle)
    assert verdict == {
        "verified": True,
        "payload_count": 4,
        "structural_records_verified": 256,
        "strict_exact_witnesses_verified": 2,
    }
    summary = json.loads(
        (expanded_bundle / "summary.json").read_text(encoding="utf-8")
    )
    assert summary["strict_hit"]["relay_two_signed_flux"] == "25/486"
    assert summary["strict_hit"]["next_signed_flux"] == "125/10206"
    assert summary["strict_hit"]["full_energy_defect"] == "0"
    assert summary["strict_hit"]["simple_binary_recursion_rejected"] is True
    assert summary["quantitative_obstruction"][
        "low_to_intended_power_ratio"
    ] == "459/106"
    assert summary["relay_quality_screen"]["passing_orientation_count"] == 0
    assert summary["independent_verification"]["verified"] is True


def test_expanded_bundle_checksum_tampering_is_rejected(
    expanded_bundle: Path, tmp_path: Path
) -> None:
    copied = tmp_path / "checksum_tamper"
    shutil.copytree(expanded_bundle, copied)
    path = copied / "summary.json"
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="checksum"):
        verify_expanded_carrier_bundle(copied)


def test_expanded_bundle_resigned_semantic_tampering_is_rejected(
    expanded_bundle: Path, tmp_path: Path
) -> None:
    copied = tmp_path / "semantic_tamper"
    shutil.copytree(expanded_bundle, copied)
    path = copied / "independent_verification.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["strict_exact_witnesses_verified"] = 999
    write_with_digest(path, canonical_json_bytes(payload))
    manifest_path = copied / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"][path.name] = {
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    write_with_digest(manifest_path, canonical_json_bytes(manifest))
    with pytest.raises(ValueError, match="independent verification"):
        verify_expanded_carrier_bundle(copied)


def test_expanded_driver_refuses_to_overwrite_evidence(tmp_path: Path) -> None:
    output = tmp_path / "occupied"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("preserve", encoding="utf-8")
    with pytest.raises(ValueError, match="new or empty"):
        run(Path("configs/expanded_carrier_search_v1.json"), output)
    assert marker.read_text(encoding="utf-8") == "preserve"


def test_expanded_driver_requires_reviewed_scope(tmp_path: Path) -> None:
    config = json.loads(
        Path("configs/expanded_carrier_search_v1.json").read_text(
            encoding="utf-8"
        )
    )
    config["reviewed_experiment_specification"]["reviewed"] = False
    path = tmp_path / "unreviewed.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    with pytest.raises(ValueError, match="reviewed finite specification"):
        run(path, tmp_path / "output")
