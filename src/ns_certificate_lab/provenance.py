"""Stable runtime and source provenance for generated research artifacts."""

from __future__ import annotations

import contextlib
import copy
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform as platform_module
import subprocess
import sys
from typing import Any, Iterator

import numpy as np

PROVENANCE_SCHEMA_V1 = "ns-certificate-lab/runtime-provenance/v1"
PROVENANCE_SCHEMA_V2 = "ns-certificate-lab/runtime-provenance/v2"
PROVENANCE_SCHEMA = PROVENANCE_SCHEMA_V2
SOURCE_FINGERPRINT_ALGORITHM = "sha256-path-and-content-v1"
SOURCE_SCOPE = (
    "pyproject.toml",
    "Dockerfile",
    ".dockerignore",
    "requirements",
    "src",
    "experiments",
    "configs",
    "scripts",
)
GIT_ABSENT = "unborn/not-a-git-checkout"
THREAD_ENVIRONMENT_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
    "OMP_DYNAMIC",
    "MKL_DYNAMIC",
)
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
_PROVENANCE_V1_KEYS = frozenset(
    {
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
)
_PROVENANCE_V2_KEYS = _PROVENANCE_V1_KEYS | frozenset(
    {
        "machine",
        "architecture",
        "byteorder",
        "thread_environment",
        "numpy_build_config",
        "numpy_build_config_sha256",
    }
)


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


def _normalize_json_value(value: Any) -> Any:
    """Return a deterministic JSON-compatible copy of NumPy build metadata."""

    if isinstance(value, dict):
        return {
            str(key): _normalize_json_value(value[key])
            for key in sorted(value, key=lambda item: str(item))
        }
    if isinstance(value, (list, tuple)):
        return [_normalize_json_value(item) for item in value]
    if isinstance(value, np.generic):
        return _normalize_json_value(value.item())
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    return str(value)


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _collect_numpy_build_config() -> dict[str, Any]:
    """Capture NumPy's build/BLAS/LAPACK report in a stable JSON structure."""

    configuration = getattr(np.__config__, "CONFIG", None)
    if isinstance(configuration, dict) and configuration:
        return {
            "source": "numpy.__config__.CONFIG",
            "configuration": _normalize_json_value(configuration),
        }

    # NumPy versions without CONFIG still expose show_config().  Normalize line
    # endings and trailing whitespace so repeated collection is byte-stable.
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        np.show_config()
    lines = [line.rstrip() for line in stream.getvalue().splitlines()]
    text = "\n".join(lines).strip()
    return {
        "source": "numpy.show_config",
        "configuration": {"text": text},
    }


def collect_runtime_provenance(
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Collect runtime, VCS and stable source-input provenance."""

    root = (project_root or discover_project_root()).resolve()
    fingerprint, file_count, scope = source_fingerprint(root)
    git_head, git_dirty = _git_state(root)
    numpy_build_config = _collect_numpy_build_config()
    return {
        "schema": PROVENANCE_SCHEMA,
        "python_version": platform_module.python_version(),
        "python_implementation": platform_module.python_implementation(),
        "numpy_version": np.__version__,
        "platform": platform_module.platform(),
        "machine": platform_module.machine() or "unknown",
        "architecture": platform_module.architecture()[0] or "unknown",
        "byteorder": sys.byteorder,
        "thread_environment": {
            name: os.environ.get(name)
            for name in THREAD_ENVIRONMENT_VARIABLES
        },
        "numpy_build_config": numpy_build_config,
        "numpy_build_config_sha256": hashlib.sha256(
            _canonical_json_bytes(numpy_build_config)
        ).hexdigest(),
        "git_head": git_head,
        "git_dirty": git_dirty,
        "source_fingerprint_algorithm": SOURCE_FINGERPRINT_ALGORITHM,
        "source_fingerprint_sha256": fingerprint,
        "source_fingerprint_file_count": file_count,
        "source_fingerprint_scope": scope,
    }


def validate_runtime_provenance(value: Any) -> dict[str, Any]:
    """Strictly validate v2 or a legacy v1 object and defensively copy it."""

    if not isinstance(value, dict):
        raise ValueError("runtime provenance must be an object")
    schema = value.get("schema")
    if schema == PROVENANCE_SCHEMA_V1:
        required = _PROVENANCE_V1_KEYS
    elif schema == PROVENANCE_SCHEMA_V2:
        required = _PROVENANCE_V2_KEYS
    else:
        raise ValueError("runtime provenance has an invalid schema")
    if set(value) != required:
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

    if schema == PROVENANCE_SCHEMA_V2:
        for key in ("machine", "architecture"):
            if not isinstance(value.get(key), str) or not value[key].strip():
                raise ValueError(
                    f"runtime provenance {key} must be a nonempty string"
                )
        if value.get("byteorder") not in {"little", "big"}:
            raise ValueError("runtime provenance byteorder is invalid")
        thread_environment = value.get("thread_environment")
        if (
            not isinstance(thread_environment, dict)
            or set(thread_environment) != set(THREAD_ENVIRONMENT_VARIABLES)
            or any(
                item is not None and not isinstance(item, str)
                for item in thread_environment.values()
            )
        ):
            raise ValueError("runtime provenance thread_environment is invalid")

        numpy_build_config = value.get("numpy_build_config")
        if (
            not isinstance(numpy_build_config, dict)
            or set(numpy_build_config) != {"source", "configuration"}
            or numpy_build_config.get("source")
            not in {"numpy.__config__.CONFIG", "numpy.show_config"}
            or not isinstance(numpy_build_config.get("configuration"), dict)
            or not numpy_build_config["configuration"]
        ):
            raise ValueError("runtime provenance NumPy build config is invalid")
        try:
            encoded_build_config = _canonical_json_bytes(numpy_build_config)
        except (TypeError, ValueError):
            raise ValueError(
                "runtime provenance NumPy build config is invalid"
            ) from None
        build_config_sha256 = value.get("numpy_build_config_sha256")
        if (
            not isinstance(build_config_sha256, str)
            or len(build_config_sha256) != 64
            or any(
                character not in "0123456789abcdef"
                for character in build_config_sha256
            )
            or build_config_sha256
            != hashlib.sha256(encoded_build_config).hexdigest()
        ):
            raise ValueError(
                "runtime provenance NumPy build config checksum is invalid"
            )
    return copy.deepcopy(value)
