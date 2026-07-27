"""Stable runtime and source provenance for generated research artifacts."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import platform as platform_module
import subprocess
from typing import Any, Iterator

import numpy as np

PROVENANCE_SCHEMA = "ns-certificate-lab/runtime-provenance/v1"
SOURCE_FINGERPRINT_ALGORITHM = "sha256-path-and-content-v1"
SOURCE_SCOPE = ("pyproject.toml", "src", "experiments", "configs")
GIT_ABSENT = "unborn/not-a-git-checkout"
_EXCLUDED_DIRECTORY_NAMES = frozenset(
    {
        "__pycache__",
        ".pytest_cache",
        ".mypy_cache",
        ".ruff_cache",
        ".git",
        "outputs",
    }
)
_EXCLUDED_SUFFIXES = frozenset({".pyc", ".pyo"})


def discover_project_root() -> Path:
    """Find the source checkout containing this installed module."""

    module_path = Path(__file__).resolve()
    for parent in module_path.parents:
        if (parent / "pyproject.toml").is_file() and (parent / "src").is_dir():
            return parent
    # A wheel may not retain pyproject.toml.  Hash the installed package tree
    # rather than pretending that unavailable checkout inputs were inspected.
    return module_path.parent


def _iter_source_files(project_root: Path) -> Iterator[tuple[str, Path]]:
    root = project_root.resolve()
    candidates: list[Path] = []
    if root.name == "ns_certificate_lab" and (root / "__init__.py").is_file():
        candidates.extend(
            path
            for path in root.rglob("*")
            if path.is_file() or path.is_symlink()
        )
    else:
        for entry in SOURCE_SCOPE:
            path = root / entry
            if path.is_file() or path.is_symlink():
                candidates.append(path)
            elif path.is_dir():
                candidates.extend(
                    nested
                    for nested in path.rglob("*")
                    if nested.is_file() or nested.is_symlink()
                )

    for path in sorted(candidates, key=lambda item: item.as_posix()):
        relative_parts = path.relative_to(root).parts
        if any(
            part in _EXCLUDED_DIRECTORY_NAMES or part.endswith(".egg-info")
            for part in relative_parts
        ):
            continue
        if path.suffix.lower() in _EXCLUDED_SUFFIXES:
            continue
        yield path.relative_to(root).as_posix(), path


def source_fingerprint(
    project_root: Path | None = None,
) -> tuple[str, int, list[str]]:
    """Hash stable relative paths and bytes of executable project inputs.

    Generated ``outputs/`` data is never traversed.  Symbolic links are hashed
    by link text and are not followed outside the project.
    """

    root = (project_root or discover_project_root()).resolve()
    digest = hashlib.sha256()
    digest.update((SOURCE_FINGERPRINT_ALGORITHM + "\0").encode("ascii"))
    files = list(_iter_source_files(root))
    for relative, path in files:
        relative_bytes = relative.encode("utf-8")
        digest.update(len(relative_bytes).to_bytes(8, "big"))
        digest.update(relative_bytes)
        if path.is_symlink():
            payload = ("symlink:" + os.readlink(path)).encode("utf-8")
        else:
            payload = path.read_bytes()
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    scope = (
        list(SOURCE_SCOPE)
        if root.name != "ns_certificate_lab"
        else ["installed-package:ns_certificate_lab"]
    )
    return digest.hexdigest(), len(files), scope


def _run_git(root: Path, *arguments: str) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(
            ["git", "-C", str(root), *arguments],
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
        )
    except (FileNotFoundError, OSError, subprocess.SubprocessError):
        return None


def _git_state(root: Path) -> tuple[str, bool | None]:
    inside = _run_git(root, "rev-parse", "--is-inside-work-tree")
    if inside is None or inside.returncode != 0 or inside.stdout.strip() != "true":
        return GIT_ABSENT, None
    head_result = _run_git(root, "rev-parse", "--verify", "HEAD")
    head = (
        head_result.stdout.strip().lower()
        if head_result is not None and head_result.returncode == 0
        else GIT_ABSENT
    )
    status = _run_git(root, "status", "--porcelain", "--untracked-files=all")
    dirty = (
        bool(status.stdout.strip())
        if status is not None and status.returncode == 0
        else None
    )
    return head, dirty


def collect_runtime_provenance(
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Collect runtime, VCS and stable source-input provenance."""

    root = (project_root or discover_project_root()).resolve()
    fingerprint, file_count, scope = source_fingerprint(root)
    git_head, git_dirty = _git_state(root)
    return {
        "schema": PROVENANCE_SCHEMA,
        "python_version": platform_module.python_version(),
        "python_implementation": platform_module.python_implementation(),
        "numpy_version": np.__version__,
        "platform": platform_module.platform(),
        "git_head": git_head,
        "git_dirty": git_dirty,
        "source_fingerprint_algorithm": SOURCE_FINGERPRINT_ALGORITHM,
        "source_fingerprint_sha256": fingerprint,
        "source_fingerprint_file_count": file_count,
        "source_fingerprint_scope": scope,
    }


def validate_runtime_provenance(value: Any) -> dict[str, Any]:
    """Strictly validate and defensively copy a provenance object."""

    if not isinstance(value, dict):
        raise ValueError("runtime provenance must be an object")
    required = {
        "schema",
        "python_version",
        "python_implementation",
        "numpy_version",
        "platform",
        "git_head",
        "git_dirty",
        "source_fingerprint_algorithm",
        "source_fingerprint_sha256",
        "source_fingerprint_file_count",
        "source_fingerprint_scope",
    }
    if set(value) != required or value.get("schema") != PROVENANCE_SCHEMA:
        raise ValueError("runtime provenance has an invalid schema")
    for key in (
        "python_version",
        "python_implementation",
        "numpy_version",
        "platform",
        "git_head",
    ):
        if not isinstance(value.get(key), str) or not value[key].strip():
            raise ValueError(f"runtime provenance {key} must be a nonempty string")
    dirty = value.get("git_dirty")
    if dirty is not None and not isinstance(dirty, bool):
        raise ValueError("runtime provenance git_dirty must be boolean or null")
    if (
        value.get("source_fingerprint_algorithm") != SOURCE_FINGERPRINT_ALGORITHM
        or not isinstance(value.get("source_fingerprint_sha256"), str)
        or len(value["source_fingerprint_sha256"]) != 64
        or any(
            character not in "0123456789abcdef"
            for character in value["source_fingerprint_sha256"]
        )
        or isinstance(value.get("source_fingerprint_file_count"), bool)
        or not isinstance(value.get("source_fingerprint_file_count"), int)
        or value["source_fingerprint_file_count"] < 1
    ):
        raise ValueError("runtime provenance source fingerprint is invalid")
    scope = value.get("source_fingerprint_scope")
    if (
        not isinstance(scope, list)
        or not scope
        or any(not isinstance(item, str) or not item for item in scope)
        or any(item == "outputs" or item.startswith("outputs/") for item in scope)
    ):
        raise ValueError("runtime provenance source scope is invalid")
    return {
        key: list(item) if key == "source_fingerprint_scope" else item
        for key, item in value.items()
    }
