"""Experiment configuration and deterministic seed recording."""

from __future__ import annotations

from pathlib import Path
import random
from typing import Any, Mapping

import numpy as np

from ._integrity import (
    canonical_json_bytes,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)
from .provenance import (
    collect_runtime_provenance,
    validate_runtime_provenance,
)

SCHEMA_V1 = "ns-certificate-lab/run-config/v1"
SCHEMA_V2 = "ns-certificate-lab/run-config/v2"
SCHEMA = SCHEMA_V2


def seed_everything(seed: int) -> np.random.Generator:
    """Seed Python and NumPy, and return the preferred explicit generator."""

    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    value = int(seed)
    if value < 0 or value > 2**32 - 1:
        raise ValueError("seed must lie in [0, 2^32 - 1]")
    random.seed(value)
    np.random.seed(value)
    return np.random.default_rng(value)


def save_run_config(
    path: str | Path,
    *,
    parameters: Mapping[str, Any],
    seed: int,
    provenance: Mapping[str, Any] | None = None,
) -> Path:
    """Write canonical JSON plus a SHA-256 sidecar.

    Provenance is collected automatically unless a validated shared run
    snapshot is supplied.
    """

    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    seed_value = int(seed)
    if seed_value < 0 or seed_value > 2**32 - 1:
        raise ValueError("seed must lie in [0, 2^32 - 1]")
    payload = {
        "schema": SCHEMA_V2,
        "seed": seed_value,
        "parameters": dict(parameters),
        "provenance": validate_runtime_provenance(
            collect_runtime_provenance()
            if provenance is None
            else dict(provenance)
        ),
    }
    try:
        encoded = canonical_json_bytes(payload)
    except (TypeError, ValueError) as exc:
        raise ValueError("configuration must contain finite JSON values") from exc
    destination = Path(path)
    if destination.suffix.lower() != ".json":
        raise ValueError("run configuration must use the .json extension")
    write_with_digest(destination, encoded)
    return destination


def load_run_config(path: str | Path) -> dict[str, Any]:
    destination = Path(path)
    verify_digest(destination)
    try:
        payload = strict_json_loads(
            destination.read_text(encoding="utf-8"),
            label="run configuration",
        )
    except UnicodeDecodeError as exc:
        raise ValueError("run configuration is not valid UTF-8 JSON") from exc
    schema = payload.get("schema") if isinstance(payload, dict) else None
    if (
        not isinstance(payload, dict)
        or schema not in {SCHEMA_V1, SCHEMA_V2}
        or isinstance(payload.get("seed"), bool)
        or not isinstance(payload.get("seed"), int)
        or not 0 <= payload["seed"] <= 2**32 - 1
        or not isinstance(payload.get("parameters"), dict)
    ):
        raise ValueError("run configuration has an invalid schema")
    if schema == SCHEMA_V2:
        if set(payload) != {"schema", "seed", "parameters", "provenance"}:
            raise ValueError("run configuration v2 has missing or unknown keys")
        payload["provenance"] = validate_runtime_provenance(payload["provenance"])
        payload["legacy_provenance_absent"] = False
    else:
        payload["provenance"] = None
        payload["legacy_provenance_absent"] = True
    return payload
