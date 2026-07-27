"""Integrity-protected JSON and CSV diagnostic output."""

from __future__ import annotations

import csv
import io
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from ._integrity import (
    canonical_json_bytes,
    require_finite_json,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)

SCHEMA = "ns-certificate-lab/diagnostics/v1"
Scalar = str | int | float | bool | None


def _is_nonfinite_csv_numeric(value: Scalar) -> bool:
    if value is None or isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return not math.isfinite(float(value))
    stripped = value.strip()
    if not stripped:
        return False
    try:
        return not math.isfinite(float(stripped))
    except ValueError:
        return False


def write_diagnostics_json(
    path: str | Path,
    *,
    records: Sequence[Mapping[str, Any]],
    metadata: Mapping[str, Any] | None = None,
) -> Path:
    """Write canonical diagnostic records and a checksum sidecar."""

    destination = Path(path)
    if destination.suffix.lower() != ".json":
        raise ValueError("JSON diagnostics must use the .json extension")
    payload = {
        "schema": SCHEMA,
        "metadata": dict(metadata or {}),
        "records": [dict(record) for record in records],
    }
    try:
        require_finite_json(payload, label="diagnostics")
        encoded = canonical_json_bytes(payload)
    except (TypeError, ValueError) as exc:
        raise ValueError("diagnostics must contain finite JSON values") from exc
    write_with_digest(destination, encoded)
    return destination


def load_diagnostics_json(path: str | Path) -> dict[str, Any]:
    destination = Path(path)
    verify_digest(destination)
    try:
        payload = strict_json_loads(
            destination.read_text(encoding="utf-8"),
            label="diagnostics",
        )
    except UnicodeDecodeError as exc:
        raise ValueError("diagnostics are not valid UTF-8 JSON") from exc
    if (
        not isinstance(payload, dict)
        or payload.get("schema") != SCHEMA
        or not isinstance(payload.get("metadata"), dict)
        or not isinstance(payload.get("records"), list)
        or not all(isinstance(record, dict) for record in payload["records"])
    ):
        raise ValueError("diagnostics have an invalid schema")
    return payload


def write_diagnostics_csv(
    path: str | Path,
    *,
    records: Sequence[Mapping[str, Scalar]],
) -> Path:
    """Write tabular scalar diagnostics and a checksum sidecar."""

    destination = Path(path)
    if destination.suffix.lower() != ".csv":
        raise ValueError("CSV diagnostics must use the .csv extension")
    rows = [dict(record) for record in records]
    if not rows:
        raise ValueError("at least one CSV diagnostic record is required")
    keys = [key for row in rows for key in row]
    if not keys or any(not isinstance(key, str) or not key for key in keys):
        raise ValueError("CSV diagnostic keys must be nonempty strings")
    fieldnames = sorted(set(keys))
    allowed = (str, int, float, bool, type(None))
    if any(not isinstance(value, allowed) for row in rows for value in row.values()):
        raise ValueError("CSV diagnostics may contain only scalar values")
    if any(_is_nonfinite_csv_numeric(value) for row in rows for value in row.values()):
        raise ValueError("CSV diagnostics may not contain NaN or Infinity")

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    write_with_digest(destination, buffer.getvalue().encode("utf-8"))
    return destination


def load_diagnostics_csv(path: str | Path) -> list[dict[str, str]]:
    destination = Path(path)
    verify_digest(destination)
    try:
        with destination.open("r", encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
    except UnicodeDecodeError as exc:
        raise ValueError("CSV diagnostics are not valid UTF-8") from exc
    for row in rows:
        for value in row.values():
            if value is None or not value.strip():
                continue
            if _is_nonfinite_csv_numeric(value):
                raise ValueError("CSV diagnostics contain a non-finite numeric value")
    return rows
