from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import zipfile

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "build_release_bundle.py"
SPEC = importlib.util.spec_from_file_location("build_release_bundle", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
release_bundle = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release_bundle)


def _archive_content(path: Path) -> dict[str, bytes]:
    with zipfile.ZipFile(path) as archive:
        assert archive.testzip() is None
        return {name: archive.read(name) for name in archive.namelist()}


def _git_show(relative: str) -> bytes:
    result = subprocess.run(
        ["git", "-C", str(PROJECT_ROOT), "show", f"HEAD:{relative}"],
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        raise AssertionError(f"missing committed release evidence: {relative}")
    return result.stdout


def test_release_bundle_is_deterministic_small_and_self_checking(
    scratch_dir: Path,
) -> None:
    first = release_bundle.build_release_bundle(
        scratch_dir / "release-a.zip",
        project_root=PROJECT_ROOT,
    )
    second = release_bundle.build_release_bundle(
        scratch_dir / "release-b.zip",
        project_root=PROJECT_ROOT,
    )
    assert first.read_bytes() == second.read_bytes()

    content = _archive_content(first)
    prefix = f"{release_bundle.ARCHIVE_ROOT}/"
    assert all(name.startswith(prefix) for name in content)
    relative_content = {
        name.removeprefix(prefix): data
        for name, data in content.items()
    }

    assert "src/ns_certificate_lab/provenance.py" in relative_content
    assert "scripts/build_release_bundle.py" in relative_content
    assert "requirements/constraints-windows-cpython311.txt" in relative_content
    assert "requirements/constraints-container.txt" in relative_content
    assert "Dockerfile" in relative_content
    assert "RELEASE_INDEX.json" in relative_content
    assert "SHA256SUMS" in relative_content

    forbidden_parts = {
        ".git",
        ".pytest_cache",
        "__pycache__",
        ".venv",
        "venv",
        "_scratch",
    }
    for relative in relative_content:
        assert not forbidden_parts.intersection(Path(relative).parts)
    index = json.loads(relative_content["RELEASE_INDEX.json"])
    assert index["schema"] == release_bundle.INDEX_SCHEMA
    assert index["archive_root"] == release_bundle.ARCHIVE_ROOT
    expected_outputs = {
        family: path.relative_to(PROJECT_ROOT).as_posix()
        for family, path in release_bundle._select_latest_outputs(  # noqa: SLF001
            PROJECT_ROOT
        ).items()
    }
    assert index["selected_outputs"] == expected_outputs
    for selected in index["selected_outputs"].values():
        assert any(
            relative.startswith(selected + "/")
            for relative in relative_content
        )
    bundled_output_directories = {
        "/".join(Path(relative).parts[:2])
        for relative in relative_content
        if relative.startswith("outputs/") and relative != "outputs/README.md"
    }
    assert bundled_output_directories == set(expected_outputs.values())

    indexed_paths = set()
    for entry in index["files"]:
        data = relative_content[entry["path"]]
        indexed_paths.add(entry["path"])
        assert entry["bytes"] == len(data)
        assert entry["sha256"] == hashlib.sha256(data).hexdigest()
    assert indexed_paths == (
        set(relative_content) - {"RELEASE_INDEX.json", "SHA256SUMS"}
    )

    checksum_entries = {}
    for line in relative_content["SHA256SUMS"].decode("ascii").splitlines():
        digest, relative = line.split("  ", maxsplit=1)
        checksum_entries[relative] = digest
    assert set(checksum_entries) == set(relative_content) - {"SHA256SUMS"}
    for relative, digest in checksum_entries.items():
        assert digest == hashlib.sha256(relative_content[relative]).hexdigest()

    for selected in index["selected_outputs"].values():
        manifest_name = f"{selected}/manifest.json"
        if manifest_name not in relative_content:
            continue
        manifest = json.loads(relative_content[manifest_name])
        for name, metadata in manifest["files"].items():
            data = relative_content[f"{selected}/{name}"]
            assert metadata["bytes"] == len(data)
            assert metadata["sha256"] == hashlib.sha256(data).hexdigest()


def test_release_bundle_refuses_overwrite_and_external_destination(
    scratch_dir: Path,
) -> None:
    destination = release_bundle.build_release_bundle(
        scratch_dir / "release.zip",
        project_root=PROJECT_ROOT,
    )
    with pytest.raises(FileExistsError, match="overwrite"):
        release_bundle.build_release_bundle(
            destination,
            project_root=PROJECT_ROOT,
        )
    with pytest.raises(ValueError, match="inside the repository"):
        release_bundle.build_release_bundle(
            PROJECT_ROOT.parent / "outside-release.zip",
            project_root=PROJECT_ROOT,
        )


def test_release_bundle_rejects_manifest_payload_hash_mismatch(
    scratch_dir: Path,
) -> None:
    output = scratch_dir / "outputs" / "synthetic_v1"
    output.mkdir(parents=True)
    payload_path = output / "summary.json"
    payload_path.write_bytes(b'{"value":1}\n')
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "files": {
                    "summary.json": {
                        "bytes": payload_path.stat().st_size,
                        "sha256": hashlib.sha256(
                            payload_path.read_bytes()
                        ).hexdigest(),
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    payload_path.write_bytes(b'{"value":2}\n')
    with pytest.raises(ValueError, match="manifest hash mismatch"):
        release_bundle._verify_included_manifests(  # noqa: SLF001
            [
                ("outputs/synthetic_v1/manifest.json", manifest_path),
                ("outputs/synthetic_v1/summary.json", payload_path),
            ]
        )


def test_evidence_git_attributes_preserve_bytes_and_trim_exports() -> None:
    paths = [
        "outputs/baseline_v5/summary.json",
        "outputs/manufactured_v5/diagnostics.json",
        "outputs/manufactured_v4/diagnostics.json",
        "certificates/README.md",
    ]
    result = subprocess.run(
        [
            "git",
            "-C",
            str(PROJECT_ROOT),
            "check-attr",
            "text",
            "export-ignore",
            "--",
            *paths,
        ],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    attributes: dict[tuple[str, str], str] = {}
    for line in result.stdout.splitlines():
        path, attribute, value = line.rsplit(": ", maxsplit=2)
        attributes[(path, attribute)] = value
    for path in paths:
        assert attributes[(path, "text")] == "unset"
    assert (
        attributes[("outputs/manufactured_v4/diagnostics.json", "export-ignore")]
        == "set"
    )
    assert (
        attributes[("outputs/manufactured_v5/diagnostics.json", "export-ignore")]
        == "unset"
    )


def test_committed_evidence_hashes_match_raw_git_blobs() -> None:
    committed_attributes = _git_show(".gitattributes")
    if b"/outputs/** -text" not in committed_attributes:
        pytest.skip("raw evidence attribute has not been committed yet")

    committed_directories = subprocess.run(
        [
            "git",
            "-C",
            str(PROJECT_ROOT),
            "ls-tree",
            "-d",
            "--name-only",
            "HEAD:outputs",
        ],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.splitlines()
    selected_directories = []
    for family in (
        *release_bundle.LATEST_OUTPUT_FAMILIES,
        *release_bundle.OPTIONAL_OUTPUT_FAMILIES,
    ):
        versions = []
        for name in committed_directories:
            match = release_bundle._OUTPUT_VERSION.fullmatch(name)  # noqa: SLF001
            if match is None or match.group("family") != family:
                continue
            version_text = match.group("version")
            versions.append(
                (
                    int(version_text) if version_text is not None else 0,
                    name,
                )
            )
        if versions:
            selected_directories.append(
                "outputs/" + max(versions, key=lambda item: item)[1]
            )

    for relative_directory in selected_directories:
        listing = subprocess.run(
            [
                "git",
                "-C",
                str(PROJECT_ROOT),
                "ls-tree",
                "-r",
                "--name-only",
                "HEAD",
                relative_directory,
            ],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        ).stdout.splitlines()
        assert listing, f"latest evidence is not committed: {relative_directory}"

        for relative in listing:
            if relative.endswith(".sha256"):
                expected = _git_show(relative).decode("ascii").strip()
                assert expected == hashlib.sha256(
                    _git_show(relative[: -len(".sha256")])
                ).hexdigest()

        experiment_manifest = f"{relative_directory}/manifest.json"
        if experiment_manifest in listing:
            manifest = json.loads(_git_show(experiment_manifest))
            for name, metadata in manifest["files"].items():
                data = _git_show(f"{relative_directory}/{name}")
                assert metadata["bytes"] == len(data)
                assert metadata["sha256"] == hashlib.sha256(data).hexdigest()

        for relative in listing:
            if relative.endswith(".manifest.json"):
                manifest = json.loads(_git_show(relative))
                archive_path = (
                    relative.rsplit("/", maxsplit=1)[0]
                    + "/"
                    + manifest["archive"]
                )
                assert manifest["archive_sha256"] == hashlib.sha256(
                    _git_show(archive_path)
                ).hexdigest()
