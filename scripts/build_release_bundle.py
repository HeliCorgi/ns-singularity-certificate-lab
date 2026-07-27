#!/usr/bin/env python3
"""Build a deterministic, review-sized release ZIP inside this repository.

The bundle contains source/review material and only the newest small evidence
directory in each preregistered output family.  Preserved older evidence stays
in Git but is intentionally not copied into the distribution archive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
from typing import Iterable
import zipfile

ARCHIVE_ROOT = "ns-singularity-certificate-lab"
INDEX_SCHEMA = "ns-certificate-lab/release-index/v1"
LATEST_OUTPUT_FAMILIES = ("manufactured", "baseline", "time_convergence")
OPTIONAL_OUTPUT_FAMILIES = ("poisson_manufactured", "nonlinear_control")
ROOT_FILES = frozenset(
    {
        ".dockerignore",
        ".gitattributes",
        ".gitignore",
        "AGENTS.md",
        "Dockerfile",
        "PLAN.md",
        "README.md",
        "REFERENCES.md",
        "SPEC.md",
        "STATUS.md",
        "pyproject.toml",
    }
)
SOURCE_DIRECTORIES = frozenset(
    {
        ".github",
        "certificates",
        "configs",
        "docs",
        "experiments",
        "requirements",
        "scripts",
        "src",
        "tests",
    }
)
EXCLUDED_PARTS = frozenset(
    {
        ".git",
        ".pytest_cache",
        ".mypy_cache",
        ".ruff_cache",
        ".hypothesis",
        ".tox",
        ".nox",
        ".venv",
        "venv",
        "__pycache__",
        "_scratch",
        "build",
        "dist",
    }
)
EXCLUDED_SUFFIXES = frozenset({".pyc", ".pyo"})
_OUTPUT_VERSION = re.compile(r"^(?P<family>[a-z_]+?)(?:_v(?P<version>[0-9]+))?$")
_FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_excluded(relative: Path) -> bool:
    return (
        any(
            part in EXCLUDED_PARTS or part.endswith(".egg-info")
            for part in relative.parts
        )
        or relative.suffix.lower() in EXCLUDED_SUFFIXES
    )


def _walk_directory(root: Path, directory: str) -> list[tuple[str, Path]]:
    base = root / directory
    if not base.exists():
        raise ValueError(f"release source directory is missing: {directory}")
    if base.is_symlink() or not base.is_dir():
        raise ValueError(f"release source directory is invalid: {directory}")
    files: list[tuple[str, Path]] = []
    for path in sorted(base.rglob("*"), key=lambda item: item.as_posix()):
        relative = path.relative_to(root)
        if _is_excluded(relative):
            continue
        if path.is_symlink():
            raise ValueError(
                f"release bundle refuses symbolic link: {relative.as_posix()}"
            )
        if path.is_file():
            files.append((relative.as_posix(), path))
    return files


def _select_latest_outputs(root: Path) -> dict[str, Path]:
    output_root = root / "outputs"
    if not output_root.is_dir():
        raise ValueError("outputs directory is missing")
    selected: dict[str, Path] = {}
    for family in (*LATEST_OUTPUT_FAMILIES, *OPTIONAL_OUTPUT_FAMILIES):
        versions: list[tuple[int, str, Path]] = []
        for path in output_root.iterdir():
            if not path.is_dir() or path.is_symlink():
                continue
            match = _OUTPUT_VERSION.fullmatch(path.name)
            if match is None or match.group("family") != family:
                continue
            version_text = match.group("version")
            version = int(version_text) if version_text is not None else 0
            versions.append((version, path.name, path))
        if not versions and family in OPTIONAL_OUTPUT_FAMILIES:
            continue
        if not versions:
            raise ValueError(f"no output directory found for family {family!r}")
        selected[family] = max(versions, key=lambda item: (item[0], item[1]))[2]
    return selected


def _collect_payload(root: Path) -> tuple[list[tuple[str, Path]], dict[str, str]]:
    payload: list[tuple[str, Path]] = []
    for name in sorted(ROOT_FILES):
        path = root / name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"required release source file is missing: {name}")
        payload.append((name, path))
    for directory in sorted(SOURCE_DIRECTORIES):
        payload.extend(_walk_directory(root, directory))

    output_readme = root / "outputs" / "README.md"
    if not output_readme.is_file() or output_readme.is_symlink():
        raise ValueError("outputs/README.md is missing or invalid")
    payload.append(("outputs/README.md", output_readme))

    selected = _select_latest_outputs(root)
    for path in selected.values():
        payload.extend(_walk_directory(root, f"outputs/{path.name}"))

    by_name: dict[str, Path] = {}
    for relative, path in payload:
        if relative in by_name:
            raise ValueError(f"duplicate release path: {relative}")
        by_name[relative] = path
    required_source = "src/ns_certificate_lab/provenance.py"
    if required_source not in by_name:
        raise ValueError(f"required release source file is missing: {required_source}")
    selected_names = {
        family: path.relative_to(root).as_posix()
        for family, path in sorted(selected.items())
    }
    return sorted(by_name.items()), selected_names


def _verify_included_sidecars(payload: Iterable[tuple[str, Path]]) -> None:
    by_name = dict(payload)
    for relative, sidecar in by_name.items():
        if not relative.startswith("outputs/") or not relative.endswith(".sha256"):
            continue
        target_relative = relative[: -len(".sha256")]
        target = by_name.get(target_relative)
        if target is None:
            raise ValueError(f"checksum sidecar has no bundled target: {relative}")
        try:
            expected = sidecar.read_text(encoding="ascii").strip()
        except UnicodeDecodeError as exc:
            raise ValueError(f"checksum sidecar is not ASCII: {relative}") from exc
        if (
            len(expected) != 64
            or any(character not in "0123456789abcdef" for character in expected)
            or expected != _sha256_file(target)
        ):
            raise ValueError(f"checksum sidecar mismatch: {relative}")


def _strict_json(path: Path, *, label: str) -> object:
    def reject_duplicate_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    def reject_constant(token: str) -> object:
        raise ValueError(f"{label} contains non-finite constant {token}")

    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=reject_duplicate_pairs,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is not strict UTF-8 JSON") from exc


def _valid_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _verify_included_manifests(payload: Iterable[tuple[str, Path]]) -> None:
    """Verify internal file hashes in experiment and candidate manifests."""

    by_name = dict(payload)
    for relative, manifest_path in by_name.items():
        if not relative.startswith("outputs/") or not relative.endswith(
            "manifest.json"
        ):
            continue
        manifest = _strict_json(manifest_path, label=relative)
        if not isinstance(manifest, dict):
            raise ValueError(f"output manifest must be an object: {relative}")
        parent = relative.rsplit("/", maxsplit=1)[0]

        if relative.endswith("/manifest.json"):
            files = manifest.get("files")
            if not isinstance(files, dict) or not files:
                raise ValueError(f"output manifest has no files table: {relative}")
            for name, metadata in files.items():
                if (
                    not isinstance(name, str)
                    or not name
                    or "/" in name
                    or "\\" in name
                    or name in {".", ".."}
                    or not isinstance(metadata, dict)
                    or set(metadata) != {"bytes", "sha256"}
                    or isinstance(metadata.get("bytes"), bool)
                    or not isinstance(metadata.get("bytes"), int)
                    or metadata["bytes"] < 0
                    or not _valid_sha256(metadata.get("sha256"))
                ):
                    raise ValueError(
                        f"output manifest file entry is invalid: {relative}"
                    )
                target_relative = f"{parent}/{name}"
                target = by_name.get(target_relative)
                if target is None:
                    raise ValueError(
                        f"output manifest target is not bundled: {target_relative}"
                    )
                if (
                    target.stat().st_size != metadata["bytes"]
                    or _sha256_file(target) != metadata["sha256"]
                ):
                    raise ValueError(
                        f"output manifest hash mismatch: {target_relative}"
                    )
            continue

        # Candidate manifests use an archive/archive_sha256 pair rather than a
        # general files table.  Verify that outer archive before distribution.
        archive_name = manifest.get("archive")
        archive_sha256 = manifest.get("archive_sha256")
        if (
            not isinstance(archive_name, str)
            or not archive_name
            or "/" in archive_name
            or "\\" in archive_name
            or archive_name in {".", ".."}
            or not _valid_sha256(archive_sha256)
        ):
            raise ValueError(
                f"candidate manifest archive entry is invalid: {relative}"
            )
        archive_relative = f"{parent}/{archive_name}"
        archive_path = by_name.get(archive_relative)
        if archive_path is None or _sha256_file(archive_path) != archive_sha256:
            raise ValueError(
                f"candidate archive hash mismatch: {archive_relative}"
            )


def _git_metadata(root: Path) -> tuple[str | None, bool | None]:
    try:
        head = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--verify", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
        )
        status = subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "status",
                "--porcelain",
                "--untracked-files=all",
            ],
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
        )
    except (FileNotFoundError, OSError, subprocess.SubprocessError):
        return None, None
    if head.returncode != 0 or status.returncode != 0:
        return None, None

    generated_prefixes = ("outputs/", "dist/", "tests/_scratch/")
    relevant_status = []
    for line in status.stdout.splitlines():
        path_text = line[3:].strip().replace("\\", "/")
        if not path_text.startswith(generated_prefixes):
            relevant_status.append(line)
    return head.stdout.strip().lower(), bool(relevant_status)


def _canonical_json(value: object) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def _zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name, date_time=_FIXED_ZIP_TIME)
    info.create_system = 3
    info.compress_type = zipfile.ZIP_STORED
    info.external_attr = (0o100644 & 0xFFFF) << 16
    return info


def _resolve_destination(root: Path, output: str | Path) -> Path:
    candidate = Path(output)
    if not candidate.is_absolute():
        candidate = root / candidate
    destination = candidate.resolve()
    if root != destination.parent and root not in destination.parents:
        raise ValueError("release bundle output must stay inside the repository")
    if destination.suffix.lower() != ".zip":
        raise ValueError("release bundle output must use the .zip extension")
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite release bundle: {destination}")
    return destination


def build_release_bundle(
    output: str | Path,
    *,
    project_root: str | Path | None = None,
    require_clean_source: bool = False,
) -> Path:
    """Build a deterministic archive and refuse unsafe paths or overwrites."""

    root = (
        Path(project_root).resolve()
        if project_root is not None
        else Path(__file__).resolve().parents[1]
    )
    if not (root / "pyproject.toml").is_file() or not (root / "src").is_dir():
        raise ValueError("project_root is not an ns-certificate-lab checkout")
    destination = _resolve_destination(root, output)
    payload, selected_outputs = _collect_payload(root)
    _verify_included_sidecars(payload)
    _verify_included_manifests(payload)
    git_head, source_dirty = _git_metadata(root)
    if require_clean_source and source_dirty is not False:
        raise ValueError(
            "release source is dirty or Git cleanliness could not be verified"
        )

    entries = []
    content: dict[str, bytes] = {}
    for relative, path in payload:
        data = path.read_bytes()
        content[relative] = data
        entries.append(
            {
                "path": relative,
                "bytes": len(data),
                "sha256": _sha256_bytes(data),
            }
        )
    index = {
        "schema": INDEX_SCHEMA,
        "archive_root": ARCHIVE_ROOT,
        "source_git_head": git_head,
        "source_worktree_dirty_excluding_generated": source_dirty,
        "selected_outputs": selected_outputs,
        "files": entries,
    }
    index_bytes = _canonical_json(index)
    content["RELEASE_INDEX.json"] = index_bytes
    checksum_lines = [
        f"{_sha256_bytes(data)}  {relative}"
        for relative, data in sorted(content.items())
    ]
    content["SHA256SUMS"] = ("\n".join(checksum_lines) + "\n").encode("ascii")

    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(destination, mode="x") as archive:
            for relative, data in sorted(content.items()):
                archive.writestr(
                    _zip_info(f"{ARCHIVE_ROOT}/{relative}"),
                    data,
                )
    except Exception:
        if destination.exists():
            destination.unlink()
        raise
    return destination


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="dist/ns-singularity-certificate-lab-release.zip",
        help="ZIP destination inside the repository (must not already exist)",
    )
    parser.add_argument(
        "--require-clean-source",
        action="store_true",
        help="reject non-generated worktree changes or unavailable Git state",
    )
    return parser.parse_args()


def main() -> None:
    arguments = _parse_arguments()
    path = build_release_bundle(
        arguments.output,
        require_clean_source=arguments.require_clean_source,
    )
    print(path)


if __name__ == "__main__":
    main()
