from __future__ import annotations

from pathlib import Path

from ns_certificate_lab.provenance import (
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
    assert len(first["source_fingerprint_sha256"]) == 64
    assert first["python_version"]
    assert first["numpy_version"]
    assert first["platform"]
    assert first["git_head"]
    assert isinstance(first["git_dirty"], (bool, type(None)))
    assert "outputs" not in first["source_fingerprint_scope"]


def test_source_fingerprint_reports_nonempty_stable_input_set() -> None:
    digest, file_count, scope = source_fingerprint()
    assert len(digest) == 64
    assert file_count > 0
    assert {"pyproject.toml", "src", "experiments", "configs"}.issubset(scope)


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
