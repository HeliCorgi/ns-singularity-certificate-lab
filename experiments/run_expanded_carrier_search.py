"""Write and independently verify the expanded exact carrier certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    digest_sidecar,
    sha256_file,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)
from ns_certificate_lab.exact_carrier_record_verifier import (
    verify_serialized_expanded_carrier_certificate,
)
from ns_certificate_lab.expanded_carrier_search import (
    ExpandedCarrierSearchConfig,
    search_expanded_known_relay,
)
from ns_certificate_lab.provenance import collect_runtime_provenance


CONFIG_SCHEMA = "ns-certificate-lab/expanded-carrier-search-config/v1"
MANIFEST_SCHEMA = "ns-certificate-lab/expanded-carrier-search-manifest/v1"
SUMMARY_SCHEMA = "ns-certificate-lab/expanded-carrier-search-summary/v1"
STATUS = "EXACT FINITE CARRIER SEARCH / NOT A PDE OR SINGULARITY PROOF"
ROOT = Path(__file__).resolve().parents[1]
PAYLOAD_NAMES = (
    "config.snapshot.json",
    "expanded_carrier_certificate.json",
    "independent_verification.json",
    "summary.json",
)
SOURCE_PATHS = (
    "src/ns_certificate_lab/_integrity.py",
    "src/ns_certificate_lab/provenance.py",
    "src/ns_certificate_lab/fourier_torus.py",
    "src/ns_certificate_lab/torus_chain.py",
    "src/ns_certificate_lab/exact_carrier_search.py",
    "src/ns_certificate_lab/expanded_carrier_search.py",
    "src/ns_certificate_lab/exact_carrier_record_verifier.py",
    "experiments/run_expanded_carrier_search.py",
    "tests/test_exact_carrier_search.py",
    "tests/test_expanded_carrier_search.py",
    "tests/test_exact_carrier_record_verifier.py",
    "tests/test_expanded_carrier_experiment.py",
)


def _load_config(path: Path) -> dict[str, Any]:
    value = strict_json_loads(path.read_text(encoding="utf-8"), label="expanded config")
    required = {
        "schema",
        "reviewed_experiment_specification",
        "wave_bound",
        "polarization_component_bound",
        "maximum_ordered_wave_pair_tests",
        "maximum_polarization_pair_tests",
        "maximum_saved_structural_records",
        "maximum_saved_unique_hits",
        "top_partial_limit",
    }
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("expanded config has missing or unknown fields")
    if value["schema"] != CONFIG_SCHEMA:
        raise ValueError("expanded config schema is invalid")
    review = value["reviewed_experiment_specification"]
    if not isinstance(review, dict) or set(review) != {
        "reviewed",
        "requested_by_user",
        "scope_fixed_before_execution",
        "recorded_date",
        "estimated_runtime_seconds",
        "claim_boundary",
    }:
        raise ValueError("expanded reviewed experiment specification is invalid")
    if any(
        review.get(name) is not True
        for name in ("reviewed", "requested_by_user", "scope_fixed_before_execution")
    ) or review.get("claim_boundary") != "finite exact carrier algebra only":
        raise ValueError("expanded search requires a reviewed finite specification")
    ExpandedCarrierSearchConfig(
        wave_bound=value["wave_bound"],
        polarization_component_bound=value["polarization_component_bound"],
        maximum_ordered_wave_pair_tests=value[
            "maximum_ordered_wave_pair_tests"
        ],
        maximum_polarization_pair_tests=value["maximum_polarization_pair_tests"],
        maximum_saved_structural_records=value[
            "maximum_saved_structural_records"
        ],
        maximum_saved_unique_hits=value["maximum_saved_unique_hits"],
        top_partial_limit=value["top_partial_limit"],
    )
    return value


def _search_config(config: dict[str, Any]) -> ExpandedCarrierSearchConfig:
    return ExpandedCarrierSearchConfig(
        wave_bound=int(config["wave_bound"]),
        polarization_component_bound=int(
            config["polarization_component_bound"]
        ),
        maximum_ordered_wave_pair_tests=int(
            config["maximum_ordered_wave_pair_tests"]
        ),
        maximum_polarization_pair_tests=int(
            config["maximum_polarization_pair_tests"]
        ),
        maximum_saved_structural_records=int(
            config["maximum_saved_structural_records"]
        ),
        maximum_saved_unique_hits=int(config["maximum_saved_unique_hits"]),
        top_partial_limit=int(config["top_partial_limit"]),
    )


def _summary(
    config: dict[str, Any],
    certificate: dict[str, object],
    verification: dict[str, object],
) -> dict[str, object]:
    preferred = certificate["strict_hits"][0]
    return {
        "schema": SUMMARY_SCHEMA,
        "status": STATUS,
        "claim_boundary": certificate["claim_boundary"],
        "reviewed_experiment_specification": config[
            "reviewed_experiment_specification"
        ],
        "search_counts": certificate["counts"],
        "rejection_histogram": certificate["rejection_histogram"],
        "strict_hit": {
            "second_left": preferred["second_left"],
            "second_right": preferred["second_right"],
            "relay_two_signed_flux": preferred["second_relay"]["interaction"][
                "signed_parent_to_child_flux"
            ],
            "next_signed_flux": preferred["next_interaction"][
                "signed_parent_to_child_flux"
            ],
            "all_unintended_parent_pairs_target_shell_zero": preferred[
                "exact_predicates"
            ]["all_unintended_parent_pairs_target_shell_zero"],
            "full_energy_defect": preferred["full_gadget"][
                "full_field_energy_defect"
            ],
            "simple_binary_recursion_rejected": preferred["full_gadget"][
                "simple_binary_recursion_rejected"
            ],
        },
        "quantitative_obstruction": preferred["four_parent_forcing_budget"],
        "relay_quality_screen": certificate["relay_quality_screen"],
        "parametric_rational_family": certificate["parametric_rational_family"],
        "continuous_polarization_classification": certificate[
            "continuous_polarization_classification"
        ],
        "independent_verification": verification,
        "provenance": collect_runtime_provenance(ROOT),
    }


def _manifest(config_path: Path, output_dir: Path) -> dict[str, object]:
    return {
        "schema": MANIFEST_SCHEMA,
        "status": STATUS,
        "source_config": {
            "requested_path": str(config_path.resolve()),
            "snapshot_sha256": sha256_file(output_dir / "config.snapshot.json"),
        },
        "source_files": {
            name: sha256_file(ROOT / name) for name in SOURCE_PATHS
        },
        "files": {
            name: {
                "bytes": (output_dir / name).stat().st_size,
                "sha256": sha256_file(output_dir / name),
            }
            for name in PAYLOAD_NAMES
        },
    }


def verify_expanded_carrier_bundle(output_dir: Path) -> dict[str, object]:
    output_dir = Path(output_dir)
    data_names = (*PAYLOAD_NAMES, "manifest.json")
    expected = set(data_names)
    expected.update(digest_sidecar(Path(name)).name for name in data_names)
    if not output_dir.is_dir() or {path.name for path in output_dir.iterdir()} != expected:
        raise ValueError("expanded bundle file inventory is invalid")
    for name in data_names:
        verify_digest(output_dir / name)
    manifest = strict_json_loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8"),
        label="expanded manifest",
    )
    if not isinstance(manifest, dict) or manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("expanded manifest schema is invalid")
    if set(manifest) != {"schema", "status", "source_config", "source_files", "files"}:
        raise ValueError("expanded manifest has missing or unknown fields")
    if manifest["status"] != STATUS or set(manifest["files"]) != set(PAYLOAD_NAMES):
        raise ValueError("expanded manifest inventory is inconsistent")
    for name, record in manifest["files"].items():
        path = output_dir / name
        if record != {"bytes": path.stat().st_size, "sha256": sha256_file(path)}:
            raise ValueError("expanded manifest payload checksum is inconsistent")
    if set(manifest["source_files"]) != set(SOURCE_PATHS):
        raise ValueError("expanded manifest source inventory is inconsistent")
    for name, digest in manifest["source_files"].items():
        if sha256_file(ROOT / name) != digest:
            raise ValueError(f"expanded source checksum mismatch for {name}")
    config = _load_config(output_dir / "config.snapshot.json")
    if manifest["source_config"]["snapshot_sha256"] != sha256_file(
        output_dir / "config.snapshot.json"
    ):
        raise ValueError("expanded config snapshot hash is inconsistent")
    certificate = strict_json_loads(
        (output_dir / "expanded_carrier_certificate.json").read_text(
            encoding="utf-8"
        ),
        label="expanded carrier certificate",
    )
    recomputed = search_expanded_known_relay(_search_config(config))
    if certificate != recomputed:
        raise ValueError("expanded search certificate does not recompute")
    verification = verify_serialized_expanded_carrier_certificate(certificate)
    saved_verification = strict_json_loads(
        (output_dir / "independent_verification.json").read_text(encoding="utf-8"),
        label="expanded independent verification",
    )
    if saved_verification != verification:
        raise ValueError("expanded independent verification is inconsistent")
    summary = strict_json_loads(
        (output_dir / "summary.json").read_text(encoding="utf-8"),
        label="expanded summary",
    )
    if not isinstance(summary, dict) or summary.get("schema") != SUMMARY_SCHEMA:
        raise ValueError("expanded summary schema is invalid")
    if summary.get("search_counts") != certificate["counts"] or summary.get(
        "independent_verification"
    ) != verification:
        raise ValueError("expanded summary is inconsistent")
    return {
        "verified": True,
        "payload_count": len(PAYLOAD_NAMES),
        "structural_records_verified": verification[
            "structural_records_verified"
        ],
        "strict_exact_witnesses_verified": verification[
            "strict_exact_witnesses_verified"
        ],
    }


def run(config_path: Path, output_dir: Path) -> dict[str, object]:
    config = _load_config(config_path)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("output directory must be new or empty")
    output_dir.mkdir(parents=True, exist_ok=True)
    certificate = search_expanded_known_relay(_search_config(config))
    verification = verify_serialized_expanded_carrier_certificate(certificate)
    payloads = {
        "config.snapshot.json": config,
        "expanded_carrier_certificate.json": certificate,
        "independent_verification.json": verification,
        "summary.json": _summary(config, certificate, verification),
    }
    for name, payload in payloads.items():
        write_with_digest(output_dir / name, canonical_json_bytes(payload))
    write_with_digest(
        output_dir / "manifest.json",
        canonical_json_bytes(_manifest(config_path, output_dir)),
    )
    return verify_expanded_carrier_bundle(output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()
    print(json.dumps(run(arguments.config, arguments.output_dir), sort_keys=True))


if __name__ == "__main__":
    main()
