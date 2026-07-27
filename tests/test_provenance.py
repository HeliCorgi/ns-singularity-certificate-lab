from __future__ import annotations

import hashlib
import json
from pathlib import Path

from ns_certificate_lab.provenance import (
    PROVENANCE_SCHEMA_V1,
    PROVENANCE_SCHEMA_V2,
    THREAD_ENVIRONMENT_VARIABLES,
    collect_runtime_provenance,
    source_fingerprint,
    validate_runtime_provenance,
)


def test_runtime_provenance_is_strict_and_source_hash_is_stable(
    scratch_dir: Path,
) -> None:
    first = validate_runtime_provenance(collect_runtime_provenance())
    scratch_file = scratch_dir / "excluded-generated-data.txt"
    scratch_file.write_text("this test artifact is not an execution input", encoding="utf-8")
    second = validate_runtime_provenance(collect_runtime_provenance())

    assert first["source_fingerprint_sha256"] == second["source_fingerprint_sha256"]
    assert (
        first["numpy_build_config_sha256"]
        == second["numpy_build_config_sha256"]
    )
    assert len(first["source_fingerprint_sha256"]) == 64
    assert first["python_version"]
    assert first["numpy_version"]
    assert first["platform"]
    assert first["schema"] == PROVENANCE_SCHEMA_V2
    assert first["machine"]
    assert first["architecture"]
    assert first["byteorder"] in {"little", "big"}
    assert first["git_head"]
    assert isinstance(first["git_dirty"], (bool, type(None)))
    assert "outputs" not in first["source_fingerprint_scope"]
    assert set(first["thread_environment"]) == set(THREAD_ENVIRONMENT_VARIABLES)
    build_config = first["numpy_build_config"]
    assert set(build_config) == {"source", "configuration"}
    assert build_config["configuration"]
    encoded = json.dumps(
        build_config,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    assert (
        first["numpy_build_config_sha256"]
        == hashlib.sha256(encoded).hexdigest()
    )
    if build_config["source"] == "numpy.__config__.CONFIG":
        dependencies = build_config["configuration"].get("Build Dependencies")
        assert isinstance(dependencies, dict)
        assert {"blas", "lapack"}.issubset(dependencies)


def test_source_fingerprint_reports_nonempty_stable_input_set() -> None:
    digest, file_count, scope = source_fingerprint()
    assert len(digest) == 64
    assert file_count > 0
    assert {
        "pyproject.toml",
        "Dockerfile",
        ".dockerignore",
        "requirements",
        "src",
        "experiments",
        "configs",
        "scripts",
    }.issubset(scope)


def test_source_fingerprint_ignores_editable_install_metadata(
    scratch_dir: Path,
) -> None:
    project = scratch_dir / "project"
    package = project / "src" / "example"
    package.mkdir(parents=True)
    (project / "pyproject.toml").write_text(
        "[project]\nname='example'\n",
        encoding="utf-8",
    )
    (package / "__init__.py").write_text("VALUE = 1\n", encoding="utf-8")
    before = source_fingerprint(project)

    egg_info = project / "src" / "example.egg-info"
    egg_info.mkdir()
    (egg_info / "PKG-INFO").write_text(
        "Metadata-Version: 2.1\nName: example\n",
        encoding="utf-8",
    )
    after = source_fingerprint(project)
    assert before == after


def test_runtime_provenance_records_thread_environment(
    monkeypatch,
) -> None:
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    monkeypatch.delenv("OPENBLAS_NUM_THREADS", raising=False)
    provenance = validate_runtime_provenance(collect_runtime_provenance())
    assert provenance["thread_environment"]["OMP_NUM_THREADS"] == "3"
    assert provenance["thread_environment"]["OPENBLAS_NUM_THREADS"] is None


def test_legacy_v1_runtime_provenance_remains_readable() -> None:
    current = collect_runtime_provenance()
    legacy_only = {
        key: value
        for key, value in current.items()
        if key
        not in {
            "machine",
            "architecture",
            "byteorder",
            "thread_environment",
            "numpy_build_config",
            "numpy_build_config_sha256",
        }
    }
    legacy_only["schema"] = PROVENANCE_SCHEMA_V1
    validated = validate_runtime_provenance(legacy_only)
    assert validated == legacy_only
    assert validated["schema"] == PROVENANCE_SCHEMA_V1


def test_recorded_v1_artifact_runtime_provenance_remains_readable() -> None:
    manifest_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "manufactured_v5"
        / "manufactured_candidate.manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validated = validate_runtime_provenance(manifest["provenance"])
    assert validated["schema"] == PROVENANCE_SCHEMA_V1
    assert validated["source_fingerprint_scope"] == [
        "pyproject.toml",
        "src",
        "experiments",
        "configs",
    ]


def test_runtime_provenance_v2_rejects_build_config_tampering() -> None:
    provenance = collect_runtime_provenance()
    provenance["numpy_build_config"]["configuration"]["tampered"] = True
    try:
        validate_runtime_provenance(provenance)
    except ValueError as exc:
        assert "checksum" in str(exc)
    else:
        raise AssertionError("tampered NumPy build configuration was accepted")


def test_runtime_provenance_validation_returns_a_deep_copy() -> None:
    provenance = collect_runtime_provenance()
    validated = validate_runtime_provenance(provenance)
    validated["thread_environment"]["OMP_NUM_THREADS"] = "changed"
    validated["numpy_build_config"]["configuration"]["changed"] = True
    assert (
        provenance["thread_environment"]["OMP_NUM_THREADS"]
        != "changed"
    )
    assert "changed" not in provenance["numpy_build_config"]["configuration"]
