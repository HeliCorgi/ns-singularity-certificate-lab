#!/usr/bin/env python
r"""Frozen-coefficient von Neumann audit of the RECORDED operating points of the
shipped Heun evidence bundles (P0-A items 2-4).

What this experiment does
-------------------------
``ns_certificate_lab.von_neumann`` answers, for a single frozen operating point
``(max|u^r|, max|u^z|, dr, dz, dt, nu)``, whether Heun with centered differences
amplifies any resolved Fourier mode of the model advection--diffusion symbol.
Until now that machinery had never been pointed at the numbers the shipped runs
actually ran with.  This experiment does exactly that: it re-reads the recorded
``diagnostics.csv`` of each shipped Heun bundle, replays every recorded row
through :func:`ns_certificate_lab.von_neumann.audit_snapshot`, and attaches a
per-row and a per-run verdict to it.

Vocabulary, fixed here and used everywhere below
------------------------------------------------
A row or a run that fails this audit is **stability-unverified**, NOT
**unstable**.  The frozen-coefficient model discards the variable-coefficient
commutators, the axis row (E-17 coefficient-8 limit), the wall rows
(E-27/E-31), the per-stage elliptic solve and the zeroth-order production
terms, so a failure is grounds for suspicion and for re-running at a smaller
step, and never a proof that a recorded trajectory grew.  A pass is equally
limited: it verifies the model at the sampled wavenumbers, not the solver.

The stride limitation, stated up front
--------------------------------------
The audited bundles record a diagnostic row only every ``stride`` accepted
steps (25 for every bundle audited here).  Every statement produced by this
experiment is therefore a statement about the RECORDED rows only.  The steps
between two recorded rows are NOT covered by any number in this bundle; they
carry their own ``dt``, their own ``max|u^r|`` and their own ``max|u^z|``, none
of which were written down.  The per-step ``step_stream`` instrumentation added
later closes this gap for future runs; it cannot close it retroactively, so the
per-run verdict is named ``stability-verified-at-recorded-rows`` rather than
``stability-verified``.

The cumulative growth question
------------------------------
The honest cumulative worst case would be the product of the per-step
amplification factors over every accepted step, and those factors were not
recorded.  What this experiment reports instead is
``exp(sum over recorded rows of ln(maxG) * stride)``, stored under the field
name ``strided_growth_extrapolation_not_a_bound``.  The name is the warning:
that expression assumes each recorded row's worst case persisted unchanged for
the whole stride, which is an assumption and not a bound in either direction.

This experiment reads ``outputs/`` and never writes into an audited bundle.
It is a numerical observation about a linear model problem.  It is not evidence
for or against singularity formation.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import json
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
import sys
import time
from typing import Any, Mapping, Sequence

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.von_neumann import METHODS, audit_snapshot, scan_amplification

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPECTED_EXPERIMENT_ID = "von_neumann_audit_v1"
SUMMARY_SCHEMA_VERSION = 1
MANIFEST_SCHEMA = "ns-certificate-lab/von-neumann-audit-manifest/v1"

#: The only diagnostics file this experiment reads out of an audited bundle.
DIAGNOSTICS_NAME = "diagnostics.csv"

#: Columns every audited diagnostics file must carry.
REQUIRED_ROW_COLUMNS: tuple[str, ...] = (
    "step",
    "time",
    "dt",
    "viscosity",
    "max_abs_u_r",
    "max_abs_u_z",
)

#: Per-run verdicts.  Deliberately NOT called "stable"/"unstable".
VERDICT_VERIFIED = "stability-verified-at-recorded-rows"
VERDICT_UNVERIFIED = "stability-unverified"

#: ``audit_snapshot`` does not expose ``n_theta``; it always uses the module
#: default of :func:`scan_amplification`.  The config must agree with that
#: default so the recorded ``n_theta`` is the one that was actually scanned.
AUDIT_N_THETA = int(
    inspect.signature(scan_amplification).parameters["n_theta"].default
)

#: Acceptance check 5 requires this exact phrase in the summary interpretation.
UNVERIFIED_PHRASE = "stability-unverified, not unstable"

INTERPRETATION = (
    "Frozen-coefficient von Neumann audit applied to the recorded operating "
    "points of shipped Heun runs. A failing row or run is "
    f"{UNVERIFIED_PHRASE}: the model freezes the measured velocity maxima over "
    "the whole grid and omits the variable-coefficient commutators, the axis "
    "row (E-17 coefficient-8 limit), the wall rows (E-27/E-31), the per-stage "
    "elliptic solve and the zeroth-order production terms, so a failure is "
    "grounds for suspicion and for a smaller step and never a proof that a "
    "recorded trajectory grew. A passing run is verified only AT THE RECORDED "
    "ROWS: the audited bundles record one diagnostic row every stride accepted "
    "steps and the steps in between were never written down, so no number here "
    "covers them. Nothing in this bundle is evidence for or against "
    "singularity formation."
)

#: Headline reference points transcribed from the prose audit so that a reader
#: can calibrate the per-row numbers below against the already-published pair.
#: These are NOT recomputed here; they are the shipped self-consistent reading.
SHIPPED_REFERENCE_POINTS: dict[str, Any] = {
    "source": "docs/numerical_stability_audit.md §1.2",
    "operating_point": (
        "the self-consistent reading of outputs/hou_early_time_v1 "
        "(min dt = 2.7587e-7, max advective CFL = 0.10023) at tolerance 1e-12"
    ),
    "radial_binding": {
        "reading": "the recorded advective CFL is carried by the radial direction",
        "heun_full_step_max_amplification": 1.000003460085721,
        "passes": False,
    },
    "axial_binding": {
        "reading": "the recorded advective CFL is carried by the axial direction",
        "heun_full_step_max_amplification": 1.0,
        "passes": True,
    },
    "note": (
        "the verdict of the shipped operating point turns in the fifth to "
        "sixth significant digit and depends on which direction is assumed to "
        "carry the recorded CFL; the per-run numbers in this bundle use the "
        "recorded max|u^r| and max|u^z| of each row simultaneously and are "
        "therefore at least as pessimistic as the radial-binding reading"
    ),
}

STRIDED_EXTRAPOLATION_NOTE = (
    "exp(sum over recorded rows of ln(max|G|) * stride). The per-step "
    "amplification factors of the steps between two recorded rows were never "
    "recorded, so this expression assumes each recorded row's frozen worst "
    "case persisted unchanged for its whole stride. That is a strided "
    "extrapolation, not a bound, in either direction."
)

LIMITATIONS: tuple[str, ...] = (
    "STRIDE. Only the recorded diagnostic rows are audited; the audited "
    "bundles record one row every 25 accepted steps, and the roughly 24 steps "
    "between two records carry their own dt, max|u^r| and max|u^z|, none of "
    "which were written down. No verdict here covers them. The per-step "
    "step_stream instrumentation closes this gap for future runs only.",
    "FROZEN COEFFICIENTS. The measured velocity maxima are frozen over the "
    "whole grid. Frozen-coefficient analysis discards the commutator terms, so "
    "it is neither an upper nor a lower bound for the variable-coefficient "
    "operator that the solver actually applied.",
    "BOUNDARY ROWS. The axis row (E-17 coefficient-8 limit) and the wall rows "
    "(E-27/E-31) are boundary modifications of the stencil. No Fourier mode "
    "diagonalizes them and nothing in this bundle says anything about them; "
    "the frozen 3*nu/r first-derivative term is a proxy, not the axis row.",
    "OMITTED OPERATORS. The per-stage elliptic solve -L5 psi1 = omega1 and the "
    "zeroth-order production terms 2 u1 psi1_z and d_z(u1^2) are not modelled, "
    "which makes the model optimistic in those directions.",
    "SAMPLED SUPREMUM. Each scan maximizes over a finite "
    f"{AUDIT_N_THETA}x{AUDIT_N_THETA} wavenumber sample containing theta = 0 "
    "and theta = +/-pi exactly, so every reported max|G| is a lower bound for "
    "the true frozen-coefficient supremum.",
    "ACCEPTANCE. 'accepted' in this bundle means the audit executed correctly "
    "on every recorded row. It does NOT mean the audited runs are stable.",
)


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------


def _nonempty_text(value: Any, *, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _finite_positive(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite positive number")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a finite positive number")
    return result


def _bundle_directory(entry: Any) -> Path:
    """Return the resolved directory of a repository-relative bundle path."""

    text = _nonempty_text(entry, name="source bundle")
    candidate = Path(text)
    rooted = (
        candidate.is_absolute()
        or bool(candidate.anchor)
        or PurePosixPath(text).is_absolute()
        or bool(PureWindowsPath(text).anchor)
    )
    if rooted or ".." in candidate.parts:
        raise ValueError(
            "source bundle entries must be repository-relative paths without "
            "parent traversal"
        )
    resolved = (REPOSITORY_ROOT / candidate).resolve()
    try:
        resolved.relative_to(REPOSITORY_ROOT)
    except ValueError as error:
        raise ValueError(
            "source bundle entries must resolve inside this repository"
        ) from error
    return resolved


def _validate_source(source: Any) -> None:
    keys = {
        "label",
        "bundle",
        "stride",
        "r_max",
        "z_period",
        "grid",
        "extra_run_key_columns",
    }
    if not isinstance(source, dict) or set(source) != keys:
        raise ValueError("each source has missing or unknown keys")
    _nonempty_text(source["label"], name="source label")
    _bundle_directory(source["bundle"])
    stride = source["stride"]
    if isinstance(stride, bool) or not isinstance(stride, int) or stride < 1:
        raise ValueError("source stride must be an integer of at least 1")
    _finite_positive(source["r_max"], name="source r_max")
    _finite_positive(source["z_period"], name="source z_period")
    grid = source["grid"]
    if grid is not None:
        if (
            not isinstance(grid, list)
            or len(grid) != 2
            or any(
                isinstance(item, bool) or not isinstance(item, int) or item < 3
                for item in grid
            )
        ):
            raise ValueError(
                "source grid must be null or a pair of integers of at least 3"
            )
    columns = source["extra_run_key_columns"]
    if not isinstance(columns, list) or any(
        not isinstance(item, str) or not item for item in columns
    ):
        raise ValueError("extra_run_key_columns must be a list of column names")
    if len(set(columns)) != len(columns):
        raise ValueError("extra_run_key_columns must not repeat a column")


def validate_config(config: Mapping[str, Any]) -> None:
    """Reject any config with missing, unknown or out-of-range entries."""

    required = {
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "sources",
        "methods",
        "tolerance",
        "n_theta",
    }
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError("von neumann audit config has missing or unknown keys")
    if config["schema_version"] != SUMMARY_SCHEMA_VERSION:
        raise ValueError("schema_version must be 1")
    if config["experiment_id"] != EXPECTED_EXPERIMENT_ID:
        raise ValueError(f"experiment_id must be {EXPECTED_EXPERIMENT_ID}")
    for key in ("description", "interpretation"):
        _nonempty_text(config[key], name=key)
    if UNVERIFIED_PHRASE not in config["interpretation"]:
        raise ValueError(
            "config interpretation must contain the exact phrase "
            f"{UNVERIFIED_PHRASE!r}"
        )

    sources = config["sources"]
    if not isinstance(sources, list) or not sources:
        raise ValueError("sources must be a nonempty list")
    for source in sources:
        _validate_source(source)
    labels = [source["label"] for source in sources]
    if len(set(labels)) != len(labels):
        raise ValueError("source labels must be unique")
    bundles = [source["bundle"] for source in sources]
    if len(set(bundles)) != len(bundles):
        raise ValueError("source bundles must be unique")

    methods = config["methods"]
    if not isinstance(methods, list) or methods != ["heun"]:
        raise ValueError(
            "methods must be exactly [\"heun\"]: every audited bundle was "
            "produced with Heun and the per-row schema records the Heun full "
            "step together with its explicit-Euler predictor stage"
        )
    if any(method not in METHODS for method in methods):
        raise ValueError("methods contains a name von_neumann does not support")

    tolerance = config["tolerance"]
    if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)):
        raise ValueError("tolerance must be a nonnegative finite number")
    if not math.isfinite(float(tolerance)) or float(tolerance) < 0.0:
        raise ValueError("tolerance must be a nonnegative finite number")

    n_theta = config["n_theta"]
    if isinstance(n_theta, bool) or not isinstance(n_theta, int):
        raise ValueError("n_theta must be an integer")
    if n_theta != AUDIT_N_THETA:
        raise ValueError(
            "n_theta must equal the von_neumann.scan_amplification default of "
            f"{AUDIT_N_THETA}: audit_snapshot does not expose n_theta, so any "
            "other value would be recorded but not scanned"
        )


# --------------------------------------------------------------------------
# reading the audited bundles
# --------------------------------------------------------------------------


def verify_bundle(directory: Path) -> dict[str, Any]:
    """Verify one audited bundle's manifest and its diagnostics digest.

    The manifest sidecar is checked first, then the diagnostics file is hashed
    and compared against the hash the verified manifest records for it.  Any
    mismatch raises; nothing downstream runs on unverified bytes.
    """

    manifest_path = directory / "manifest.json"
    diagnostics_path = directory / DIAGNOSTICS_NAME
    if not directory.is_dir():
        raise ValueError(f"audited bundle is missing: {directory}")
    if not diagnostics_path.is_file():
        raise ValueError(f"audited bundle has no {DIAGNOSTICS_NAME}: {directory}")
    verify_digest(manifest_path)
    manifest = strict_json_loads(
        manifest_path.read_text(encoding="utf-8"),
        label=f"manifest of {directory.name}",
    )
    if not isinstance(manifest, dict) or not isinstance(manifest.get("files"), dict):
        raise ValueError(f"manifest of {directory.name} has no files mapping")
    entry = manifest["files"].get(DIAGNOSTICS_NAME)
    if not isinstance(entry, dict) or not isinstance(entry.get("sha256"), str):
        raise ValueError(
            f"manifest of {directory.name} does not record {DIAGNOSTICS_NAME}"
        )
    observed = sha256_file(diagnostics_path)
    if observed != entry["sha256"]:
        raise ValueError(
            f"checksum mismatch for {DIAGNOSTICS_NAME} of {directory.name}: "
            f"manifest records {entry['sha256']}, file hashes to {observed}"
        )
    return {
        "manifest_sha256": sha256_file(manifest_path),
        "manifest_digest_verified": True,
        "diagnostics_sha256": observed,
        "diagnostics_matches_manifest": True,
        "experiment_id": manifest.get("experiment_id"),
    }


def _row_float(row: Mapping[str, str], name: str, *, where: str) -> float:
    text = row.get(name)
    if text is None or text == "":
        raise ValueError(f"{where}: column {name!r} is missing or empty")
    try:
        value = float(text)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{where}: column {name!r} is not a number") from error
    if not math.isfinite(value):
        raise ValueError(f"{where}: column {name!r} is not finite")
    return value


def _row_grid(
    row: Mapping[str, str],
    source: Mapping[str, Any],
    *,
    where: str,
) -> tuple[int, int]:
    grid = source["grid"]
    if grid is not None:
        return int(grid[0]), int(grid[1])
    for name in ("nr", "nz"):
        if row.get(name) in (None, ""):
            raise ValueError(
                f"{where}: source declares no grid and the row has no "
                f"{name!r} column"
            )
    nr = int(float(row["nr"]))
    nz = int(float(row["nz"]))
    if nr < 3 or nz < 3:
        raise ValueError(f"{where}: recorded grid {nr}x{nz} is too small")
    return nr, nz


def read_source_rows(source: Mapping[str, Any]) -> tuple[list[dict[str, Any]], int]:
    """Return the auditable rows of one bundle and the dt == 0 exclusion count.

    Rows whose recorded ``dt`` is exactly zero are step-0 records: no step was
    taken, so there is no operating point to audit.  They are excluded and
    counted rather than silently dropped.  A negative or non-finite ``dt`` is
    an error, not an exclusion.
    """

    directory = _bundle_directory(source["bundle"])
    path = directory / DIAGNOSTICS_NAME
    stride = int(source["stride"])
    r_max = float(source["r_max"])
    z_period = float(source["z_period"])
    extra_columns: Sequence[str] = tuple(source["extra_run_key_columns"])

    rows: list[dict[str, Any]] = []
    excluded = 0
    with path.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        header = reader.fieldnames or []
        missing = [name for name in REQUIRED_ROW_COLUMNS if name not in header]
        if missing:
            raise ValueError(
                f"{path.name} of {directory.name} lacks columns {missing}"
            )
        for name in extra_columns:
            if name not in header:
                raise ValueError(
                    f"{path.name} of {directory.name} lacks the declared run "
                    f"key column {name!r}"
                )
        for index, raw in enumerate(reader):
            where = f"{source['label']} row {index}"
            dt = _row_float(raw, "dt", where=where)
            if dt < 0.0:
                raise ValueError(f"{where}: recorded dt is negative")
            nr, nz = _row_grid(raw, source, where=where)
            extra = {name: str(raw[name]) for name in extra_columns}
            key_parts = [f"nr{nr}", f"nz{nz}"]
            key_parts.extend(f"{name}={value}" for name, value in extra.items())
            run_id = f"{source['label']}/" + "/".join(key_parts)
            if dt == 0.0:
                excluded += 1
                rows.append(
                    {
                        "run_id": run_id,
                        "nr": nr,
                        "nz": nz,
                        "extra": extra,
                        "excluded": True,
                    }
                )
                continue
            viscosity = _row_float(raw, "viscosity", where=where)
            speed_r = _row_float(raw, "max_abs_u_r", where=where)
            speed_z = _row_float(raw, "max_abs_u_z", where=where)
            if viscosity < 0.0 or speed_r < 0.0 or speed_z < 0.0:
                raise ValueError(f"{where}: viscosity and speeds must be >= 0")
            rows.append(
                {
                    "run_id": run_id,
                    "label": str(source["label"]),
                    "bundle": str(source["bundle"]),
                    "nr": nr,
                    "nz": nz,
                    "extra": extra,
                    "stride": stride,
                    "dr": r_max / float(nr - 1),
                    "dz": z_period / float(nz),
                    "step": _row_float(raw, "step", where=where),
                    "time": _row_float(raw, "time", where=where),
                    "dt": dt,
                    "viscosity": viscosity,
                    "max_abs_u_r": speed_r,
                    "max_abs_u_z": speed_z,
                    "excluded": False,
                }
            )
    if not any(not row["excluded"] for row in rows):
        raise ValueError(
            f"{directory.name} has no auditable row with a positive dt"
        )
    return rows, excluded


# --------------------------------------------------------------------------
# the audit itself
# --------------------------------------------------------------------------


def audit_row(
    row: Mapping[str, Any],
    *,
    tolerance: float,
    methods: Sequence[str],
) -> dict[str, Any]:
    """Audit one recorded operating point with :func:`audit_snapshot`."""

    result = audit_snapshot(
        max_abs_u_r=float(row["max_abs_u_r"]),
        max_abs_u_z=float(row["max_abs_u_z"]),
        dr=float(row["dr"]),
        dz=float(row["dz"]),
        dt=float(row["dt"]),
        viscosity=float(row["viscosity"]),
        tolerance=float(tolerance),
        methods=tuple(methods),
    )
    if int(result["n_theta"]) != AUDIT_N_THETA:
        raise ValueError(
            "audit_snapshot scanned "
            f"{result['n_theta']} wavenumbers, expected {AUDIT_N_THETA}"
        )
    heun = result["methods"]["heun"]
    record = {
        "run_id": row["run_id"],
        "label": row["label"],
        "bundle": row["bundle"],
        "nr": int(row["nr"]),
        "nz": int(row["nz"]),
        "run_key_extra": ";".join(
            f"{name}={value}" for name, value in row["extra"].items()
        ),
        "step": float(row["step"]),
        "time": float(row["time"]),
        "dt": float(row["dt"]),
        "viscosity": float(row["viscosity"]),
        "max_abs_u_r": float(row["max_abs_u_r"]),
        "max_abs_u_z": float(row["max_abs_u_z"]),
        "dr": float(result["dr"]),
        "dz": float(result["dz"]),
        "advective_cfl_r": float(result["advective_cfl_r"]),
        "advective_cfl_z": float(result["advective_cfl_z"]),
        "viscous_number": float(result["viscous_number"]),
        "first_derivative_coefficient_r": float(
            result["first_derivative_coefficient_r"]
        ),
        "heun_max_amplification": float(heun["max_amplification"]),
        "heun_argmax_theta_r": float(heun["argmax_theta_r"]),
        "heun_argmax_theta_z": float(heun["argmax_theta_z"]),
        "heun_passes": bool(heun["passes"]),
        "predictor_stage": str(heun["predictor_stage"]),
        "predictor_max_amplification": float(heun["predictor_max_amplification"]),
        "predictor_argmax_theta_r": float(heun["predictor_argmax_theta_r"]),
        "predictor_argmax_theta_z": float(heun["predictor_argmax_theta_z"]),
        "row_verdict": str(result["verdict"]),
    }
    return record


_FINITE_ROW_FIELDS: tuple[str, ...] = (
    "heun_max_amplification",
    "heun_argmax_theta_r",
    "heun_argmax_theta_z",
    "predictor_max_amplification",
    "predictor_argmax_theta_r",
    "predictor_argmax_theta_z",
    "advective_cfl_r",
    "advective_cfl_z",
    "viscous_number",
)


def _row_is_finite(record: Mapping[str, Any]) -> bool:
    return all(math.isfinite(float(record[name])) for name in _FINITE_ROW_FIELDS)


def _strided_extrapolation(
    records: Sequence[Mapping[str, Any]],
    stride: int,
) -> tuple[float, float | None]:
    """Return the strided log-growth sum and its exponential, or ``None``.

    ``None`` is returned for the exponential when it would overflow binary64,
    which keeps the summary strictly JSON-representable instead of emitting an
    ``Infinity`` token.
    """

    total = 0.0
    for record in records:
        magnitude = float(record["heun_max_amplification"])
        if magnitude <= 0.0:
            raise ValueError("a recorded amplification factor was not positive")
        total += math.log(magnitude) * float(stride)
    if not math.isfinite(total):
        raise ValueError("the strided log-growth sum is not finite")
    try:
        value = math.exp(total)
    except OverflowError:
        return total, None
    return total, (value if math.isfinite(value) else None)


def summarize_run(
    run_id: str,
    records: Sequence[Mapping[str, Any]],
    *,
    excluded: int,
    stride: int,
    source: Mapping[str, Any],
) -> dict[str, Any]:
    """Aggregate one run's audited rows into a per-run verdict."""

    if not records:
        raise ValueError(f"run {run_id} has no audited rows")
    worst = max(records, key=lambda item: float(item["heun_max_amplification"]))
    worst_predictor = max(
        records, key=lambda item: float(item["predictor_max_amplification"])
    )
    failing = [item for item in records if not bool(item["heun_passes"])]
    log_sum, extrapolation = _strided_extrapolation(records, stride)
    verdict = VERDICT_UNVERIFIED if failing else VERDICT_VERIFIED
    return {
        "run_id": run_id,
        "label": str(source["label"]),
        "bundle": str(source["bundle"]),
        "nr": int(records[0]["nr"]),
        "nz": int(records[0]["nz"]),
        "run_key_extra": str(records[0]["run_key_extra"]),
        "dr": float(records[0]["dr"]),
        "dz": float(records[0]["dz"]),
        "stride": int(stride),
        "rows_audited": len(records),
        "rows_excluded_dt_zero": int(excluded),
        "failing_rows": len(failing),
        "passing_rows": len(records) - len(failing),
        "min_dt": min(float(item["dt"]) for item in records),
        "max_dt": max(float(item["dt"]) for item in records),
        "min_viscosity": min(float(item["viscosity"]) for item in records),
        "max_viscosity": max(float(item["viscosity"]) for item in records),
        "max_advective_cfl_r": max(
            float(item["advective_cfl_r"]) for item in records
        ),
        "max_advective_cfl_z": max(
            float(item["advective_cfl_z"]) for item in records
        ),
        "worst_row": {
            "step": float(worst["step"]),
            "time": float(worst["time"]),
            "dt": float(worst["dt"]),
            "viscosity": float(worst["viscosity"]),
            "max_abs_u_r": float(worst["max_abs_u_r"]),
            "max_abs_u_z": float(worst["max_abs_u_z"]),
            "heun_max_amplification": float(worst["heun_max_amplification"]),
            "heun_argmax_theta_r": float(worst["heun_argmax_theta_r"]),
            "heun_argmax_theta_z": float(worst["heun_argmax_theta_z"]),
            "heun_passes": bool(worst["heun_passes"]),
        },
        "worst_predictor_row": {
            "step": float(worst_predictor["step"]),
            "time": float(worst_predictor["time"]),
            "dt": float(worst_predictor["dt"]),
            "predictor_stage": str(worst_predictor["predictor_stage"]),
            "predictor_max_amplification": float(
                worst_predictor["predictor_max_amplification"]
            ),
            "predictor_argmax_theta_r": float(
                worst_predictor["predictor_argmax_theta_r"]
            ),
            "predictor_argmax_theta_z": float(
                worst_predictor["predictor_argmax_theta_z"]
            ),
        },
        "verdict": verdict,
        "verdict_meaning": (
            f"{VERDICT_UNVERIFIED} means the frozen-coefficient worst case was "
            "not verified at one or more recorded rows; it is grounds for "
            "suspicion and for a smaller step, not a proof that the run was "
            f"unstable. {VERDICT_VERIFIED} covers the recorded rows only."
        ),
        "strided_log_growth_sum_not_a_bound": log_sum,
        "strided_growth_extrapolation_not_a_bound": extrapolation,
        "strided_growth_extrapolation_note": STRIDED_EXTRAPOLATION_NOTE,
    }


def evaluate(
    config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, bool]]:
    """Audit every configured bundle without writing any evidence.

    Returns ``(source_records, run_records, row_records, checks)``.
    """

    validate_config(config)
    tolerance = float(config["tolerance"])
    methods = [str(name) for name in config["methods"]]

    source_records: list[dict[str, Any]] = []
    run_records: list[dict[str, Any]] = []
    row_records: list[dict[str, Any]] = []
    manifests_verified = True

    for source in config["sources"]:
        directory = _bundle_directory(source["bundle"])
        integrity = verify_bundle(directory)
        manifests_verified = manifests_verified and bool(
            integrity["manifest_digest_verified"]
            and integrity["diagnostics_matches_manifest"]
        )
        rows, excluded_total = read_source_rows(source)

        grouped: dict[str, list[dict[str, Any]]] = {}
        exclusions: dict[str, int] = {}
        order: list[str] = []
        for row in rows:
            run_id = str(row["run_id"])
            if run_id not in grouped:
                grouped[run_id] = []
                exclusions[run_id] = 0
                order.append(run_id)
            if row["excluded"]:
                exclusions[run_id] += 1
                continue
            record = audit_row(row, tolerance=tolerance, methods=methods)
            grouped[run_id].append(record)
            row_records.append(record)

        run_ids: list[str] = []
        for run_id in order:
            run_records.append(
                summarize_run(
                    run_id,
                    grouped[run_id],
                    excluded=exclusions[run_id],
                    stride=int(source["stride"]),
                    source=source,
                )
            )
            run_ids.append(run_id)

        source_records.append(
            {
                "label": str(source["label"]),
                "bundle": str(source["bundle"]),
                "stride": int(source["stride"]),
                "r_max": float(source["r_max"]),
                "z_period": float(source["z_period"]),
                "grid": source["grid"],
                "extra_run_key_columns": list(source["extra_run_key_columns"]),
                "integrity": integrity,
                "rows_audited": sum(len(grouped[key]) for key in order),
                "rows_excluded_dt_zero": int(excluded_total),
                "run_ids": run_ids,
            }
        )

    checks = {
        "all_input_manifests_verified": bool(manifests_verified),
        "all_audited_rows_finite": all(
            _row_is_finite(record) for record in row_records
        ),
        "per_run_verdict_recorded": all(
            record.get("verdict") in (VERDICT_VERIFIED, VERDICT_UNVERIFIED)
            for record in run_records
        )
        and bool(run_records),
        "summary_reports_verified_and_unverified_booleans": True,
        "interpretation_states_unverified_not_unstable": (
            UNVERIFIED_PHRASE in INTERPRETATION
        ),
    }
    return source_records, run_records, row_records, checks


# --------------------------------------------------------------------------
# evidence bundle
# --------------------------------------------------------------------------

ROW_FIELDNAMES: tuple[str, ...] = (
    "bundle",
    "label",
    "run_id",
    "nr",
    "nz",
    "run_key_extra",
    "step",
    "time",
    "dt",
    "viscosity",
    "max_abs_u_r",
    "max_abs_u_z",
    "dr",
    "dz",
    "advective_cfl_r",
    "advective_cfl_z",
    "viscous_number",
    "first_derivative_coefficient_r",
    "heun_max_amplification",
    "heun_argmax_theta_r",
    "heun_argmax_theta_z",
    "heun_passes",
    "predictor_stage",
    "predictor_max_amplification",
    "predictor_argmax_theta_r",
    "predictor_argmax_theta_z",
    "row_verdict",
)


def _prepare_output(output_dir: Path) -> None:
    if output_dir.exists():
        if not output_dir.is_dir():
            raise FileExistsError(f"output path is not a directory: {output_dir}")
        if any(output_dir.iterdir()):
            raise FileExistsError(
                f"output directory {output_dir} is not empty; refusing to "
                "overwrite existing evidence"
            )
    else:
        output_dir.mkdir(parents=True)


def _inside_repository(path: Path) -> bool:
    try:
        path.resolve().relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def _refuse_audited_destination(config: Mapping[str, Any], output_dir: Path) -> None:
    resolved = output_dir.resolve()
    for source in config["sources"]:
        audited = _bundle_directory(source["bundle"])
        if resolved == audited or audited in resolved.parents:
            raise ValueError(
                "refusing to write inside an audited source evidence directory"
            )


def run(config: Mapping[str, Any], output_dir: Path) -> dict[str, Any]:
    """Audit every configured bundle and write a checksummed evidence bundle."""

    validate_config(config)
    if not _inside_repository(output_dir):
        raise ValueError("output directory must remain inside this repository")
    _refuse_audited_destination(config, output_dir)
    provenance = collect_runtime_provenance()
    config_bytes = canonical_json_bytes(dict(config))
    _prepare_output(output_dir)

    started = time.perf_counter()
    sources, runs, rows, checks = evaluate(config)
    elapsed = time.perf_counter() - started

    any_verified = any(record["verdict"] == VERDICT_VERIFIED for record in runs)
    all_unverified = all(record["verdict"] == VERDICT_UNVERIFIED for record in runs)

    summary: dict[str, Any] = {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "experiment_id": config["experiment_id"],
        "description": config["description"],
        "interpretation": INTERPRETATION,
        "config_interpretation": config["interpretation"],
        "method": {
            "auditor": (
                "ns_certificate_lab.von_neumann.audit_snapshot with "
                f"methods={list(config['methods'])}, tolerance="
                f"{float(config['tolerance'])!r}, n_theta={AUDIT_N_THETA}"
            ),
            "operating_point": (
                "the recorded max_abs_u_r and max_abs_u_z of each diagnostics "
                "row are frozen over the whole grid and fed in simultaneously, "
                "with dr = r_max/(nr-1) and dz = z_period/nz"
            ),
            "predictor_stage": (
                "the explicit-Euler predictor of Heun is scanned separately, "
                "because a stage can overflow even when the completed step "
                "does not"
            ),
            "exclusions": (
                "rows whose recorded dt is exactly zero are step-0 records "
                "with no step taken; they are excluded and counted"
            ),
            "run_verdict_rule": (
                f"{VERDICT_VERIFIED} if every audited row of the run passes at "
                f"the configured tolerance, otherwise {VERDICT_UNVERIFIED}"
            ),
        },
        "shipped_reference_points": dict(SHIPPED_REFERENCE_POINTS),
        "sources": sources,
        "runs": runs,
        "totals": {
            "runs": len(runs),
            "rows_audited": len(rows),
            "rows_excluded_dt_zero": sum(
                int(record["rows_excluded_dt_zero"]) for record in sources
            ),
            "failing_rows": sum(int(record["failing_rows"]) for record in runs),
            "runs_stability_unverified": sum(
                1 for record in runs if record["verdict"] == VERDICT_UNVERIFIED
            ),
            "runs_stability_verified_at_recorded_rows": sum(
                1 for record in runs if record["verdict"] == VERDICT_VERIFIED
            ),
        },
        "any_run_stability_verified_at_recorded_rows": bool(any_verified),
        "all_runs_stability_unverified": bool(all_unverified),
        "acceptance_checks": checks,
        "acceptance_meaning": (
            "'accepted' means this audit executed correctly on every recorded "
            "row of every verified input bundle. It does NOT mean the audited "
            "runs are stable; see the per-run verdicts."
        ),
        "accepted_as_von_neumann_audit": bool(all(checks.values())),
        "limitations": list(LIMITATIONS),
        "runtime": {
            "audit_wall_time_seconds": float(elapsed),
            "scans_performed": 2 * len(rows),
            "scans_note": (
                "one full-step and one predictor-stage scan per audited row, "
                f"each over a {AUDIT_N_THETA}x{AUDIT_N_THETA} wavenumber sample"
            ),
        },
        "reproducibility": {
            "config_sha256": sha256_bytes(config_bytes),
            "runtime_provenance": provenance,
            "command": (
                "python experiments/run_von_neumann_audit.py --config "
                "configs/von_neumann_audit.json --output-dir "
                "outputs/von_neumann_audit_v1"
            ),
        },
    }

    write_with_digest(output_dir / "config.snapshot.json", config_bytes)
    with (output_dir / "row_audits.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ROW_FIELDNAMES))
        writer.writeheader()
        writer.writerows(
            {name: record[name] for name in ROW_FIELDNAMES} for record in rows
        )
    write_with_digest(output_dir / "summary.json", canonical_json_bytes(summary))

    artifact_names = [
        "config.snapshot.json",
        "config.snapshot.json.sha256",
        "row_audits.csv",
        "summary.json",
        "summary.json.sha256",
    ]
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "experiment_id": config["experiment_id"],
        "files": {
            name: {
                "sha256": sha256_file(output_dir / name),
                "bytes": (output_dir / name).stat().st_size,
            }
            for name in sorted(artifact_names)
        },
        "audited_inputs": {
            record["bundle"]: {
                "manifest_sha256": record["integrity"]["manifest_sha256"],
                "diagnostics_sha256": record["integrity"]["diagnostics_sha256"],
            }
            for record in sources
        },
    }
    write_with_digest(output_dir / "manifest.json", canonical_json_bytes(manifest))
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "von_neumann_audit.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "von_neumann_audit_v1",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    try:
        config = strict_json_loads(
            config_path.read_text(encoding="utf-8"),
            label="von neumann audit config",
        )
        if not isinstance(config, dict):
            raise ValueError("config must be a JSON object")
        summary = run(config, output_dir)
    except (OSError, ValueError, FileExistsError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary["acceptance_checks"], indent=2, sort_keys=True))
    print(
        json.dumps(
            {
                "any_run_stability_verified_at_recorded_rows": summary[
                    "any_run_stability_verified_at_recorded_rows"
                ],
                "all_runs_stability_unverified": summary[
                    "all_runs_stability_unverified"
                ],
                "runs": {
                    record["run_id"]: record["verdict"] for record in summary["runs"]
                },
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if summary["accepted_as_von_neumann_audit"] else 2


if __name__ == "__main__":
    sys.exit(main())
