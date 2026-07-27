"""Small integrity and canonical-serialization helpers."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import uuid
from typing import Any


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    text = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return (text + "\n").encode("utf-8")


def require_finite_json(value: Any, *, label: str = "JSON value") -> None:
    """Recursively reject non-JSON types and non-finite floating values."""

    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{label} contains a non-finite number")
        return
    if isinstance(value, list):
        for item in value:
            require_finite_json(item, label=label)
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError(f"{label} contains a non-string object key")
            require_finite_json(item, label=label)
        return
    raise ValueError(f"{label} contains a value that is not standard JSON data")


def strict_json_loads(data: str, *, label: str = "JSON data") -> Any:
    """Decode JSON while rejecting NaN, Infinity and exponent overflow."""

    def reject_constant(token: str) -> None:
        raise ValueError(f"{label} contains forbidden constant {token}")

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, item in pairs:
            if key in result:
                raise ValueError(f"{label} contains duplicate object key {key!r}")
            result[key] = item
        return result

    try:
        value = json.loads(
            data,
            parse_constant=reject_constant,
            object_pairs_hook=unique_object,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} is not valid JSON") from exc
    require_finite_json(value, label=label)
    return value


def atomic_write_bytes(path: Path, data: bytes) -> None:
    """Replace ``path`` atomically using a temporary sibling file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def digest_sidecar(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".sha256")


def write_with_digest(path: Path, data: bytes) -> None:
    atomic_write_bytes(path, data)
    atomic_write_bytes(
        digest_sidecar(path),
        (sha256_bytes(data) + "\n").encode("ascii"),
    )


def verify_digest(path: Path) -> None:
    sidecar = digest_sidecar(path)
    if not path.is_file() or not sidecar.is_file():
        raise ValueError(f"missing data or checksum sidecar for {path.name}")
    expected = sidecar.read_text(encoding="ascii").strip()
    if len(expected) != 64 or any(character not in "0123456789abcdef" for character in expected):
        raise ValueError(f"invalid checksum sidecar for {path.name}")
    observed = sha256_file(path)
    if observed != expected:
        raise ValueError(f"checksum mismatch for {path.name}")
