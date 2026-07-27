"""Explicit candidate arrays with strict manifests and runtime provenance."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
from typing import Any, Mapping
import uuid

import numpy as np
import numpy.typing as npt

from ._integrity import (
    atomic_write_bytes,
    canonical_json_bytes,
    digest_sidecar,
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    verify_digest,
)
from .grid import AxisymmetricGrid, FloatArray
from .provenance import (
    collect_runtime_provenance,
    validate_runtime_provenance,
)

SCHEMA_V1 = "ns-certificate-lab/candidate/v1"
SCHEMA_V2 = "ns-certificate-lab/candidate/v2"
SCHEMA = SCHEMA_V2
_FIELD_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_RESERVED_NAMES = frozenset({"r", "z"})
_CANONICAL_DTYPE = np.dtype("<f8")


@dataclass(frozen=True)
class CandidateDescription:
    """Required interpretation of every v2 candidate array bundle."""

    representation: str
    coordinate_system: str
    units: dict[str, str]
    normalization: str
    physical_time: float
    viscosity: float
    basis_convention: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "representation": self.representation,
            "coordinate_system": self.coordinate_system,
            "units": dict(self.units),
            "normalization": self.normalization,
            "physical_time": self.physical_time,
            "viscosity": self.viscosity,
            "basis_convention": self.basis_convention,
        }

    @classmethod
    def from_value(
        cls,
        value: "CandidateDescription | Mapping[str, Any]",
        *,
        field_names: set[str],
    ) -> "CandidateDescription":
        if isinstance(value, cls):
            raw = value.to_dict()
        elif isinstance(value, Mapping):
            raw = dict(value)
        else:
            raise ValueError("candidate description must be an object")
        required_keys = {
            "representation",
            "coordinate_system",
            "units",
            "normalization",
            "physical_time",
            "viscosity",
            "basis_convention",
        }
        if set(raw) != required_keys:
            raise ValueError("candidate description has missing or unknown keys")
        for key in (
            "representation",
            "coordinate_system",
            "normalization",
            "basis_convention",
        ):
            if not isinstance(raw[key], str) or not raw[key].strip():
                raise ValueError(f"candidate description {key} must be nonempty text")
        units = raw["units"]
        required_units = {"r", "z", *field_names}
        if (
            not isinstance(units, dict)
            or any(
                not isinstance(key, str)
                or not isinstance(unit, str)
                or not unit.strip()
                for key, unit in units.items()
            )
            or not required_units.issubset(units)
        ):
            raise ValueError(
                "candidate description units must name r, z and every stored field"
            )
        physical_time = raw["physical_time"]
        viscosity = raw["viscosity"]
        if (
            isinstance(physical_time, bool)
            or not isinstance(physical_time, (int, float))
            or not np.isfinite(physical_time)
            or physical_time < 0.0
        ):
            raise ValueError("candidate physical_time must be finite and nonnegative")
        if (
            isinstance(viscosity, bool)
            or not isinstance(viscosity, (int, float))
            or not np.isfinite(viscosity)
            or viscosity <= 0.0
        ):
            raise ValueError("candidate viscosity must be finite and positive")
        return cls(
            representation=raw["representation"],
            coordinate_system=raw["coordinate_system"],
            units={key: units[key] for key in sorted(units)},
            normalization=raw["normalization"],
            physical_time=float(physical_time),
            viscosity=float(viscosity),
            basis_convention=raw["basis_convention"],
        )


@dataclass(frozen=True)
class CandidateArtifact:
    grid: AxisymmetricGrid
    fields: dict[str, FloatArray]
    config: dict[str, Any]
    seed: int
    metadata: dict[str, Any]
    schema: str
    description: CandidateDescription | None
    provenance: dict[str, Any] | None
    legacy_provenance_absent: bool


@dataclass(frozen=True)
class CandidatePaths:
    archive: Path
    manifest: Path
    manifest_checksum: Path


def candidate_paths(path: str | Path) -> CandidatePaths:
    archive = Path(path)
    if archive.suffix.lower() != ".npz":
        raise ValueError("candidate archive must use the .npz extension")
    manifest = archive.with_suffix(".manifest.json")
    return CandidatePaths(
        archive=archive,
        manifest=manifest,
        manifest_checksum=digest_sidecar(manifest),
    )


def _array_digest(array: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(array)
    return sha256_bytes(contiguous.tobytes(order="C"))


def _json_copy(value: Mapping[str, Any], *, label: str) -> dict[str, Any]:
    try:
        encoded = canonical_json_bytes(dict(value))
        decoded = strict_json_loads(encoded.decode("utf-8"), label=label)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be finite JSON data") from exc
    if not isinstance(decoded, dict):
        raise ValueError(f"{label} must be a JSON object")
    return decoded


def save_candidate(
    path: str | Path,
    *,
    grid: AxisymmetricGrid,
    fields: Mapping[str, npt.ArrayLike],
    config: Mapping[str, Any],
    seed: int,
    description: CandidateDescription | Mapping[str, Any],
    metadata: Mapping[str, Any] | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> CandidatePaths:
    """Save a v2 candidate with canonical float64 arrays and provenance.

    Candidate fields must already be explicit arrays; optimizer or neural
    network state alone is not accepted.  Every array is serialized as
    canonical little-endian IEEE-754 binary64.  The required description makes
    coordinates, units, normalization and physical interpretation explicit.
    ``provenance=None`` collects automatically; callers producing several
    files in one run may pass one pre-run snapshot to every writer.
    """

    paths = candidate_paths(path)
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    seed_value = int(seed)
    if seed_value < 0 or seed_value > 2**63 - 1:
        raise ValueError("seed must lie in [0, 2^63 - 1]")
    if not fields:
        raise ValueError("at least one candidate field is required")

    arrays: dict[str, np.ndarray] = {
        "r": np.ascontiguousarray(grid.r, dtype=_CANONICAL_DTYPE),
        "z": np.ascontiguousarray(grid.z, dtype=_CANONICAL_DTYPE),
    }
    for name, values in sorted(fields.items()):
        if not isinstance(name, str) or not _FIELD_NAME.fullmatch(name):
            raise ValueError(f"invalid field name: {name!r}")
        if name in _RESERVED_NAMES:
            raise ValueError(f"reserved field name: {name}")
        validated = grid.validate_field(values, name=name)
        arrays[name] = np.ascontiguousarray(validated, dtype=_CANONICAL_DTYPE)

    field_names = set(fields)
    description_value = CandidateDescription.from_value(
        description,
        field_names=field_names,
    )
    config_copy = _json_copy(config, label="config")
    metadata_copy = _json_copy(metadata or {}, label="metadata")
    provenance_value = validate_runtime_provenance(
        collect_runtime_provenance()
        if provenance is None
        else dict(provenance)
    )

    paths.archive.parent.mkdir(parents=True, exist_ok=True)
    temporary = paths.archive.parent / f".{paths.archive.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("wb") as stream:
            np.savez_compressed(stream, **arrays)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, paths.archive)
    finally:
        if temporary.exists():
            temporary.unlink()

    array_entries = {
        name: {
            "dtype": _CANONICAL_DTYPE.str,
            "shape": list(array.shape),
            "sha256": _array_digest(array),
        }
        for name, array in arrays.items()
    }
    manifest = {
        "schema": SCHEMA_V2,
        "archive": paths.archive.name,
        "archive_sha256": sha256_file(paths.archive),
        "grid": {
            "periodic_z": grid.periodic_z,
            "z_period": grid.z_period,
        },
        "arrays": array_entries,
        "field_names": sorted(field_names),
        "config": config_copy,
        "seed": seed_value,
        "metadata": metadata_copy,
        "description": description_value.to_dict(),
        "provenance": provenance_value,
    }
    manifest_bytes = canonical_json_bytes(manifest)
    atomic_write_bytes(paths.manifest, manifest_bytes)
    atomic_write_bytes(
        paths.manifest_checksum,
        (sha256_bytes(manifest_bytes) + "\n").encode("ascii"),
    )
    return paths


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _load_manifest(paths: CandidatePaths) -> dict[str, Any]:
    verify_digest(paths.manifest)
    try:
        value = strict_json_loads(
            paths.manifest.read_text(encoding="utf-8"),
            label="candidate manifest",
        )
    except UnicodeDecodeError as exc:
        raise ValueError("candidate manifest is not valid UTF-8 JSON") from exc
    if not isinstance(value, dict) or value.get("schema") not in {
        SCHEMA_V1,
        SCHEMA_V2,
    }:
        raise ValueError("unsupported or missing candidate schema")
    if value["schema"] == SCHEMA_V2:
        required_keys = {
            "schema",
            "archive",
            "archive_sha256",
            "grid",
            "arrays",
            "field_names",
            "config",
            "seed",
            "metadata",
            "description",
            "provenance",
        }
        if set(value) != required_keys:
            raise ValueError("candidate v2 manifest has missing or unknown keys")
    if value.get("archive") != paths.archive.name:
        raise ValueError("manifest archive name does not match candidate path")
    if (
        not _valid_sha256(value.get("archive_sha256"))
        or value["archive_sha256"] != sha256_file(paths.archive)
    ):
        raise ValueError("candidate archive checksum mismatch")
    return value


def load_candidate(path: str | Path) -> CandidateArtifact:
    """Load a strict v2 candidate or a provenance-free legacy v1 candidate."""

    paths = candidate_paths(path)
    if not paths.archive.is_file():
        raise ValueError("candidate archive is missing")
    manifest = _load_manifest(paths)
    schema = manifest["schema"]
    entries = manifest.get("arrays")
    field_names = manifest.get("field_names")
    if not isinstance(entries, dict) or not isinstance(field_names, list):
        raise ValueError("candidate manifest has invalid array metadata")
    if (
        any(not isinstance(name, str) for name in field_names)
        or any(not _FIELD_NAME.fullmatch(name) for name in field_names)
        or any(name in _RESERVED_NAMES for name in field_names)
        or len(field_names) != len(set(field_names))
        or field_names != sorted(field_names)
    ):
        raise ValueError("candidate manifest has invalid field names")
    expected_names = {"r", "z", *field_names}
    if set(entries) != expected_names:
        raise ValueError("candidate manifest array set is inconsistent")

    try:
        with np.load(paths.archive, allow_pickle=False) as archive:
            if set(archive.files) != expected_names:
                raise ValueError("candidate archive array set is inconsistent")
            arrays = {name: np.asarray(archive[name]).copy() for name in expected_names}
    except (OSError, ValueError, KeyError) as exc:
        if isinstance(exc, ValueError) and str(exc).startswith("candidate archive array"):
            raise
        raise ValueError("candidate archive is not a valid non-pickle NPZ") from exc

    for name, array in arrays.items():
        entry = entries.get(name)
        if not isinstance(entry, dict):
            raise ValueError(f"missing metadata for array {name}")
        if schema == SCHEMA_V2 and set(entry) != {"dtype", "shape", "sha256"}:
            raise ValueError(f"candidate v2 array metadata has unknown keys for {name}")
        if (
            entry.get("dtype") != array.dtype.str
            or entry.get("shape") != list(array.shape)
            or entry.get("sha256") != _array_digest(array)
        ):
            raise ValueError(f"candidate array metadata mismatch for {name}")
        if array.dtype != _CANONICAL_DTYPE:
            raise ValueError(
                f"candidate array {name} must use little-endian float64"
            )
        if not np.issubdtype(array.dtype, np.floating) or not np.all(np.isfinite(array)):
            raise ValueError(f"candidate array {name} must be finite floating point")

    grid_metadata = manifest.get("grid")
    if (
        not isinstance(grid_metadata, dict)
        or not isinstance(grid_metadata.get("periodic_z"), bool)
        or (schema == SCHEMA_V2 and set(grid_metadata) != {"periodic_z", "z_period"})
    ):
        raise ValueError("candidate grid metadata is invalid")
    periodic_z = grid_metadata["periodic_z"]
    z_period = grid_metadata.get("z_period")
    if (
        periodic_z
        and (
            isinstance(z_period, bool)
            or not isinstance(z_period, (int, float))
            or not np.isfinite(z_period)
            or z_period <= 0.0
        )
    ) or (not periodic_z and z_period is not None):
        raise ValueError("candidate grid period metadata is invalid")

    description: CandidateDescription | None
    provenance: dict[str, Any] | None
    legacy_provenance_absent = schema == SCHEMA_V1
    if schema == SCHEMA_V2:
        description = CandidateDescription.from_value(
            manifest.get("description"),
            field_names=set(field_names),
        )
        provenance = validate_runtime_provenance(manifest.get("provenance"))
    else:
        description = None
        provenance = None

    grid = AxisymmetricGrid(
        r=np.asarray(arrays.pop("r"), dtype=np.float64),
        z=np.asarray(arrays.pop("z"), dtype=np.float64),
        periodic_z=periodic_z,
        z_period=z_period,
    )
    for name, array in arrays.items():
        if array.shape != grid.shape:
            raise ValueError(f"candidate field {name} does not match the grid")

    seed = manifest.get("seed")
    config = manifest.get("config")
    metadata = manifest.get("metadata")
    if (
        isinstance(seed, bool)
        or not isinstance(seed, int)
        or seed < 0
        or seed > 2**63 - 1
        or not isinstance(config, dict)
        or not isinstance(metadata, dict)
    ):
        raise ValueError("candidate run metadata is invalid")
    config = _json_copy(config, label="config")
    metadata = _json_copy(metadata, label="metadata")
    fields_out = {name: value.copy() for name, value in arrays.items()}
    return CandidateArtifact(
        grid=grid,
        fields=fields_out,
        config=config,
        seed=seed,
        metadata=metadata,
        schema=schema,
        description=description,
        provenance=provenance,
        legacy_provenance_absent=legacy_provenance_absent,
    )
