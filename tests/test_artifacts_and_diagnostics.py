from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from ns_certificate_lab.artifacts import (
    SCHEMA_V1 as CANDIDATE_SCHEMA_V1,
    SCHEMA_V2 as CANDIDATE_SCHEMA_V2,
    CandidateDescription,
    load_candidate,
    save_candidate,
)
from ns_certificate_lab.config import (
    SCHEMA_V1 as CONFIG_SCHEMA_V1,
    load_run_config,
    save_run_config,
    seed_everything,
)
from ns_certificate_lab.diagnostics import (
    load_diagnostics_csv,
    load_diagnostics_json,
    write_diagnostics_csv,
    write_diagnostics_json,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.manufactured import manufactured_fields


def _description(*field_names: str) -> CandidateDescription:
    return CandidateDescription(
        representation="nodal little-endian float64 arrays",
        coordinate_system="axisymmetric cylindrical half-plane (r,z)",
        units={
            "r": "dimensionless length",
            "z": "dimensionless length",
            **{name: "dimensionless transformed field" for name in field_names},
        },
        normalization="dimensionless manufactured normalization",
        physical_time=0.2,
        viscosity=0.01,
        basis_convention=(
            "array axes (r,z); u^r=-r psi1_z; "
            "u^z=2 psi1+r psi1_r; -L5 psi1=omega1"
        ),
    )


def _write_json_with_valid_digest(path: Path, value: object) -> None:
    encoded = (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")
    path.write_bytes(encoded)
    path.with_suffix(path.suffix + ".sha256").write_text(
        hashlib.sha256(encoded).hexdigest() + "\n",
        encoding="ascii",
    )


def test_candidate_roundtrip_is_exact_and_records_configuration(
    scratch_dir: Path,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=32, r_max=1.0)
    fields = manufactured_fields(grid)
    path = scratch_dir / "candidate.npz"
    save_candidate(
        path,
        grid=grid,
        fields={
            "u1": fields.u1,
            "omega1": fields.omega1,
            "psi1": fields.psi1,
        },
        config={"viscosity": 0.01, "method": "manufactured"},
        seed=1729,
        description=_description("u1", "omega1", "psi1"),
        metadata={"claim": "smooth test field, not a singularity"},
    )
    loaded = load_candidate(path)

    assert loaded.seed == 1729
    assert loaded.schema == CANDIDATE_SCHEMA_V2
    assert not loaded.legacy_provenance_absent
    assert loaded.provenance is not None
    assert loaded.provenance["python_version"]
    assert loaded.provenance["numpy_version"] == np.__version__
    assert loaded.provenance["source_fingerprint_file_count"] > 0
    assert all(
        not item.startswith("outputs")
        for item in loaded.provenance["source_fingerprint_scope"]
    )
    assert loaded.description == _description("u1", "omega1", "psi1")
    assert loaded.config == {"method": "manufactured", "viscosity": 0.01}
    assert loaded.metadata["claim"].startswith("smooth")
    assert np.array_equal(loaded.grid.r, grid.r)
    assert np.array_equal(loaded.grid.z, grid.z)
    assert loaded.grid.periodic_z
    assert set(loaded.fields) == {"u1", "omega1", "psi1"}
    for name in loaded.fields:
        assert loaded.fields[name].dtype == np.dtype("<f8")
        assert np.array_equal(
            loaded.fields[name],
            getattr(fields, name),
        )


def test_invalid_candidate_archive_is_rejected(scratch_dir: Path) -> None:
    path = scratch_dir / "invalid.npz"
    path.write_bytes(b"not an npz archive")
    with pytest.raises(ValueError, match="missing|checksum|archive"):
        load_candidate(path)


def test_candidate_archive_tampering_is_detected(scratch_dir: Path) -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=16, r_max=1.0)
    fields = manufactured_fields(grid)
    path = scratch_dir / "candidate.npz"
    save_candidate(
        path,
        grid=grid,
        fields={"psi1": fields.psi1},
        config={},
        seed=1,
        description=_description("psi1"),
    )
    with path.open("ab") as stream:
        stream.write(b"tamper")
    with pytest.raises(ValueError, match="checksum"):
        load_candidate(path)


def test_malformed_manifest_field_names_are_rejected_as_validation_error(
    scratch_dir: Path,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=16, r_max=1.0)
    fields = manufactured_fields(grid)
    path = scratch_dir / "candidate.npz"
    paths = save_candidate(
        path,
        grid=grid,
        fields={"psi1": fields.psi1},
        config={},
        seed=1,
        description=_description("psi1"),
    )
    manifest = json.loads(paths.manifest.read_text(encoding="utf-8"))
    manifest["field_names"] = [["unhashable", "hostile"]]
    encoded = (
        json.dumps(
            manifest,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")
    paths.manifest.write_bytes(encoded)
    paths.manifest_checksum.write_text(
        hashlib.sha256(encoded).hexdigest() + "\n",
        encoding="ascii",
    )
    with pytest.raises(ValueError, match="field names"):
        load_candidate(path)


def test_v2_rejects_noncanonical_float32_even_with_matching_hashes(
    scratch_dir: Path,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=16, r_max=1.0)
    fields = manufactured_fields(grid)
    path = scratch_dir / "candidate.npz"
    paths = save_candidate(
        path,
        grid=grid,
        fields={"psi1": fields.psi1},
        config={},
        seed=1,
        description=_description("psi1"),
    )
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name].copy() for name in archive.files}
    arrays["psi1"] = arrays["psi1"].astype("<f4")
    with path.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    manifest = json.loads(paths.manifest.read_text(encoding="utf-8"))
    changed = arrays["psi1"]
    manifest["archive_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest["arrays"]["psi1"] = {
        "dtype": changed.dtype.str,
        "shape": list(changed.shape),
        "sha256": hashlib.sha256(
            np.ascontiguousarray(changed).tobytes(order="C")
        ).hexdigest(),
    }
    _write_json_with_valid_digest(paths.manifest, manifest)
    with pytest.raises(ValueError, match="little-endian float64"):
        load_candidate(path)


def test_v2_rejects_resigned_invalid_runtime_provenance(
    scratch_dir: Path,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=16, r_max=1.0)
    fields = manufactured_fields(grid)
    path = scratch_dir / "candidate.npz"
    paths = save_candidate(
        path,
        grid=grid,
        fields={"psi1": fields.psi1},
        config={},
        seed=1,
        description=_description("psi1"),
    )
    manifest = json.loads(paths.manifest.read_text(encoding="utf-8"))
    manifest["provenance"]["git_dirty"] = "unknown"
    _write_json_with_valid_digest(paths.manifest, manifest)
    with pytest.raises(ValueError, match="git_dirty"):
        load_candidate(path)


def test_legacy_v1_candidate_loads_with_explicit_provenance_absence(
    scratch_dir: Path,
) -> None:
    grid = AxisymmetricGrid.uniform(nr=9, nz=16, r_max=1.0)
    fields = manufactured_fields(grid)
    path = scratch_dir / "candidate.npz"
    paths = save_candidate(
        path,
        grid=grid,
        fields={"psi1": fields.psi1},
        config={},
        seed=1,
        description=_description("psi1"),
    )
    manifest = json.loads(paths.manifest.read_text(encoding="utf-8"))
    manifest["schema"] = CANDIDATE_SCHEMA_V1
    manifest.pop("description")
    manifest.pop("provenance")
    _write_json_with_valid_digest(paths.manifest, manifest)

    loaded = load_candidate(path)
    assert loaded.schema == CANDIDATE_SCHEMA_V1
    assert loaded.provenance is None
    assert loaded.description is None
    assert loaded.legacy_provenance_absent


def test_run_config_records_seed_and_reproduces_rng(scratch_dir: Path) -> None:
    path = scratch_dir / "run.json"
    save_run_config(path, parameters={"nr": 17, "viscosity": 0.01}, seed=42)
    loaded = load_run_config(path)
    first = seed_everything(loaded["seed"]).standard_normal(8)
    second = seed_everything(42).standard_normal(8)
    assert loaded["parameters"] == {"nr": 17, "viscosity": 0.01}
    assert loaded["provenance"]["source_fingerprint_sha256"]
    assert not loaded["legacy_provenance_absent"]
    assert np.array_equal(first, second)


def test_legacy_v1_run_config_loads_with_explicit_provenance_absence(
    scratch_dir: Path,
) -> None:
    path = scratch_dir / "run.json"
    save_run_config(path, parameters={"nr": 17}, seed=42)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["schema"] = CONFIG_SCHEMA_V1
    payload.pop("provenance")
    _write_json_with_valid_digest(path, payload)
    loaded = load_run_config(path)
    assert loaded["provenance"] is None
    assert loaded["legacy_provenance_absent"]


def test_json_diagnostics_roundtrip_and_tamper_detection(
    scratch_dir: Path,
) -> None:
    path = scratch_dir / "diagnostics.json"
    records = [
        {"resolution": 17, "divergence_rms": 1.25e-3},
        {"resolution": 33, "divergence_rms": 3.1e-4},
    ]
    write_diagnostics_json(
        path,
        records=records,
        metadata={"experiment": "manufactured"},
    )
    loaded = load_diagnostics_json(path)
    assert loaded["records"] == records

    original = path.read_text(encoding="utf-8")
    path.write_text(original.replace("0.00125", "0.00001"), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        load_diagnostics_json(path)


def test_csv_diagnostics_roundtrip(scratch_dir: Path) -> None:
    path = scratch_dir / "diagnostics.csv"
    write_diagnostics_csv(
        path,
        records=[
            {"resolution": 17, "passed": False},
            {"resolution": 33, "passed": True},
        ],
    )
    records = load_diagnostics_csv(path)
    assert records == [
        {"passed": "False", "resolution": "17"},
        {"passed": "True", "resolution": "33"},
    ]


@pytest.mark.parametrize("bad_value", [float("nan"), float("inf"), float("-inf")])
def test_diagnostic_and_config_writers_reject_nonfinite_values(
    scratch_dir: Path,
    bad_value: float,
) -> None:
    with pytest.raises(ValueError, match="finite|NaN|Infinity"):
        write_diagnostics_json(
            scratch_dir / "bad.json",
            records=[{"value": bad_value}],
        )
    with pytest.raises(ValueError, match="NaN|Infinity"):
        write_diagnostics_csv(
            scratch_dir / "bad.csv",
            records=[{"value": bad_value}],
        )
    with pytest.raises(ValueError, match="finite"):
        save_run_config(
            scratch_dir / "bad-config.json",
            parameters={"value": bad_value},
            seed=1,
        )


def test_checksum_valid_loaders_reject_nonfinite_json_and_csv(
    scratch_dir: Path,
) -> None:
    diagnostic_path = scratch_dir / "diagnostics.json"
    write_diagnostics_json(diagnostic_path, records=[{"value": 1.0}])
    diagnostic_payload = json.loads(diagnostic_path.read_text(encoding="utf-8"))
    diagnostic_payload["records"][0]["value"] = float("nan")
    _write_json_with_valid_digest(diagnostic_path, diagnostic_payload)
    with pytest.raises(ValueError, match="forbidden|non-finite"):
        load_diagnostics_json(diagnostic_path)

    config_path = scratch_dir / "config.json"
    save_run_config(config_path, parameters={"value": 1.0}, seed=1)
    config_payload = json.loads(config_path.read_text(encoding="utf-8"))
    config_payload["parameters"]["value"] = float("inf")
    _write_json_with_valid_digest(config_path, config_payload)
    with pytest.raises(ValueError, match="forbidden|non-finite"):
        load_run_config(config_path)

    csv_path = scratch_dir / "diagnostics.csv"
    write_diagnostics_csv(csv_path, records=[{"value": 1.0}])
    csv_data = b"value\nNaN\n"
    csv_path.write_bytes(csv_data)
    csv_path.with_suffix(".csv.sha256").write_text(
        hashlib.sha256(csv_data).hexdigest() + "\n",
        encoding="ascii",
    )
    with pytest.raises(ValueError, match="non-finite"):
        load_diagnostics_csv(csv_path)


def test_checksum_valid_run_config_rejects_duplicate_json_keys(
    scratch_dir: Path,
) -> None:
    path = scratch_dir / "config.json"
    save_run_config(path, parameters={"nr": 17}, seed=1)
    encoded = path.read_bytes().replace(b'"seed":1', b'"seed":1,"seed":2')
    assert encoded.count(b'"seed"') == 2
    path.write_bytes(encoded)
    path.with_suffix(".json.sha256").write_text(
        hashlib.sha256(encoded).hexdigest() + "\n",
        encoding="ascii",
    )
    with pytest.raises(ValueError, match="duplicate object key"):
        load_run_config(path)
