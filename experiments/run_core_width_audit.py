#!/usr/bin/env python
r"""Core-width / points-per-scale audit of the saved Hou early-time snapshots.

This experiment is the P0-D application step.  It never integrates anything.
It re-reads the *already saved* checkpoints of ``outputs/hou_early_time_v1``
and ``outputs/hou_early_time_v2_hires`` under checksum verification and applies
:mod:`ns_certificate_lab.core_width` to them, so that the question "is the
structure whose amplitude we keep quoting actually resolved by the grid it
lives on?" is answered with machine-readable numbers instead of an assertion.

What is measured, per snapshot and per resolution, on *two* fields -- the
stored swirl ``u1`` and the pointwise magnitude of the E-18b Cartesian
vorticity ``|omega|`` rebuilt from ``(u1, omega1)``:

* radial and axial full width at half maximum and the matching
  points-per-FWHM counts;
* the 10%-to-90% front thickness on the flank facing the axis and its
  points-per-``dr`` count;
* how many cells separate the peak from the axis;
* the sub-grid (three-point parabolic) peak location and value;
* the smallest pointwise gradient length scale ``|f|/|grad f|`` on the active
  region and its points-per-cell counts;
* the discrete high-frequency tail fractions in ``r`` and ``z``;
* every truncation flag raised by any of the above.

Each report is then passed through :func:`core_width.fit_precondition`, whose
threshold :data:`core_width.PREREGISTERED_MIN_POINTS_PER_FRONT` was fixed from
manufactured ``tanh`` fronts before any Hou snapshot was looked at.  The gate
can only *forbid* a convergence fit; passing it certifies nothing.

Two further tables are produced.

``common-grid differences``
    For every pair of adjacent resolutions that share a snapshot time, both
    fields are brought onto the coarser grid with
    :func:`core_width.interpolate_to_common_grid` (Catmull-Rom restriction of
    the fine field, never prolongation of the coarse one) and the ``linf``,
    radially weighted ``L2`` and derivative-``linf`` differences are recorded.
    These are differences between two unconverged discrete solutions.  They are
    not error estimates, and no order of accuracy is fitted from them here.

``P1-B amplification renormalization``
    For every resolution the summary reports side by side (a) the *discrete*
    initial maximum Cartesian vorticity on that grid, (b) the common continuum
    reference ``24000 pi 37^{-1/2} (36/37)^18`` evaluated in code from that
    closed form, (c) the absolute final maximum Cartesian vorticity at ``T_1``,
    (d) the grid-normalized amplification ``(c)/(a)`` and (e) the
    common-reference amplification ``(c)/(b)``, together with (f) the ``(r,z)``
    location of the discrete initial maximum and its distance from the
    continuum location.  The continuum location is *computed* here by nested
    grid refinement of the closed-form E-29 datum
    (:func:`continuum_initial_vorticity_maximum`); it is not quoted.

The audit reads ``outputs/`` and never writes into any audited directory.
Nothing here is evidence for or against singularity formation.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any, Mapping, Sequence

import numpy as np

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.core_width import (
    PREREGISTERED_DELTA_OVER_DR_VALUES,
    PREREGISTERED_FRONT_ERROR_TOLERANCE,
    PREREGISTERED_FRONT_PHASES,
    PREREGISTERED_MIN_POINTS_PER_FRONT,
    fit_precondition,
    interpolate_to_common_grid,
    peak_resolution_report,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    cartesian_vorticity,
    constrain_state,
    save_checkpoint,
)
from ns_certificate_lab.provenance import collect_runtime_provenance

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPECTED_EXPERIMENT_ID = "core_width_audit_v1"
SUMMARY_SCHEMA_VERSION = 1
MANIFEST_SCHEMA = "ns-certificate-lab/core-width-audit-manifest/v1"
BUNDLE_MANIFEST_KEYS = frozenset({"experiment_id", "files", "schema_version"})

#: The two fields the audit is required to report on.  ``u1`` is what the
#: solver stores; the Cartesian vorticity magnitude is what the amplification
#: numbers are built from, and it is one derivative rougher than ``u1``.
AUDITED_FIELDS: tuple[str, ...] = ("u1", "cartesian_vorticity_magnitude")

#: E-29 initial datum, closed form (docs/equation_audit.md, E-29b).
E29_INITIAL_AMPLITUDE = 12000.0
E29_AXIAL_DENOMINATOR_COEFFICIENT = 12.5
E29_RADIAL_EXPONENT = 18

#: The P1-B common reference: the continuum maximum Cartesian vorticity of the
#: E-29 datum, ``24000 pi 37^{-1/2} (36/37)^18``.  Evaluated from the closed
#: form here so that no decimal is copied between files.
E29B_MAX_CARTESIAN_VORTICITY_EXACT = (
    24000.0 * math.pi / math.sqrt(37.0) * (36.0 / 37.0) ** 18
)

#: Required verbatim in every summary this experiment writes.
AMPLIFICATION_STATEMENT = (
    "amplification numbers under both normalizations are unconverged numerical "
    "observations; the fit precondition gate "
    "(PREREGISTERED_MIN_POINTS_PER_FRONT) decides whether convergence fits may "
    "be attempted, and Hou published values are never fit targets"
)

_CHECKPOINT_NAME = re.compile(
    r"^checkpoints/checkpoint_nr(?P<nr>[0-9]+)_nz(?P<nz>[0-9]+)_"
    r"(?P<label>t[0-9]{3})\.npz$"
)
_SNAPSHOT_LABEL = re.compile(r"^t[0-9]{3}$")

LIMITATIONS: tuple[str, ...] = (
    "Every points-per-scale count is a property of the SAMPLED ARRAY.  It "
    "bounds what the grid can represent; it is not a statement about the "
    "continuum solution, and a large count is not a certificate of accuracy.",
    "The high-frequency tail fractions are computed from single snapshots on "
    "finite grids.  Aliasing folds unrepresented continuum content back into "
    "the represented modes, so a small tail fraction is consistent with badly "
    "under-resolved data; only a large one is informative, and then only as a "
    "warning.",
    "The gradient length scale is built from three-point differences, which "
    "are least accurate exactly where the length scale is smallest, and its "
    "minimum is frequently attained at the edge of the active region rather "
    "than at the narrowest feature.  It is recorded, not gated on.",
    "The common-grid differences compare two unconverged discrete solutions "
    "restricted onto the coarser mesh.  They include the restriction error of "
    "the Catmull-Rom operator and are not error estimates; no convergence "
    "order is fitted from them in this experiment.",
    "The 129x256 resolution appears in both audited bundles.  The audit checks "
    "that the two copies are byte identical and then audits one of them; it "
    "does not treat them as independent samples.",
    "fit_precondition can only forbid a convergence fit.  Passing it removes a "
    "preregistered prohibition; it never validates a fit, and it says nothing "
    "about time discretization, boundary treatment or the elliptic solve.",
    "Binary64 floating-point arithmetic with tolerances only.  There is no "
    "interval enclosure and nothing here is evidence for or against "
    "singularity formation.",
)


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------


def _nonempty_text(value: Any, *, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _bounded_int(value: Any, *, name: str, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer in [{low}, {high}]")
    if not low <= value <= high:
        raise ValueError(f"{name} must be an integer in [{low}, {high}]")
    return int(value)


def evidence_directory(entry: Any) -> Path:
    """Resolve one configured input directory, refusing anything escaping."""

    text = _nonempty_text(entry, name="input_evidence_dirs entry")
    candidate = Path(text)
    rooted = (
        candidate.is_absolute()
        or bool(candidate.anchor)
        or PurePosixPath(text).is_absolute()
        or PureWindowsPath(text).is_absolute()
        or bool(PureWindowsPath(text).anchor)
    )
    if rooted or ".." in candidate.parts:
        raise ValueError(
            "input_evidence_dirs entries must be repository-relative paths "
            "without parent traversal"
        )
    resolved = (REPOSITORY_ROOT / candidate).resolve()
    try:
        resolved.relative_to(REPOSITORY_ROOT)
    except ValueError as error:
        raise ValueError(
            "input_evidence_dirs entries must resolve inside this repository"
        ) from error
    return resolved


def _display_path(path: Path) -> str:
    """Repository-relative POSIX name when possible, absolute POSIX otherwise."""

    try:
        return path.resolve().relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def validate_config(config: Mapping[str, Any]) -> None:
    """Reject any config with missing, unknown or out-of-range entries."""

    required = {
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "input_evidence_dirs",
        "snapshot_labels",
        "initial_snapshot_label",
        "final_snapshot_label",
        "audited_fields",
        "continuum_reference",
    }
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError("core width audit config has missing or unknown keys")
    if config["schema_version"] != SUMMARY_SCHEMA_VERSION:
        raise ValueError("schema_version must be 1")
    if config["experiment_id"] != EXPECTED_EXPERIMENT_ID:
        raise ValueError(f"experiment_id must be {EXPECTED_EXPERIMENT_ID}")
    for key in ("description", "interpretation"):
        _nonempty_text(config[key], name=key)

    directories = config["input_evidence_dirs"]
    if not isinstance(directories, list) or not directories:
        raise ValueError("input_evidence_dirs must be a nonempty list")
    if len(set(map(repr, directories))) != len(directories):
        raise ValueError("input_evidence_dirs entries must be unique")
    for entry in directories:
        evidence_directory(entry)

    labels = config["snapshot_labels"]
    if not isinstance(labels, list) or not labels:
        raise ValueError("snapshot_labels must be a nonempty list")
    for label in labels:
        if not isinstance(label, str) or not _SNAPSHOT_LABEL.fullmatch(label):
            raise ValueError("snapshot_labels entries must look like 't000'")
    if len(set(labels)) != len(labels) or labels != sorted(labels):
        raise ValueError("snapshot_labels must be unique and sorted ascending")
    for key in ("initial_snapshot_label", "final_snapshot_label"):
        label = config[key]
        if label not in labels:
            raise ValueError(f"{key} must appear in snapshot_labels")
    if config["initial_snapshot_label"] >= config["final_snapshot_label"]:
        raise ValueError(
            "initial_snapshot_label must precede final_snapshot_label"
        )

    fields = config["audited_fields"]
    if not isinstance(fields, list) or list(fields) != list(AUDITED_FIELDS):
        raise ValueError(
            "audited_fields must be exactly " + repr(list(AUDITED_FIELDS))
        )

    reference = config["continuum_reference"]
    reference_keys = {"search_points_per_axis", "refinement_stages"}
    if not isinstance(reference, dict) or set(reference) != reference_keys:
        raise ValueError("continuum_reference has missing or unknown keys")
    _bounded_int(
        reference["search_points_per_axis"],
        name="continuum_reference.search_points_per_axis",
        low=101,
        high=4001,
    )
    _bounded_int(
        reference["refinement_stages"],
        name="continuum_reference.refinement_stages",
        low=1,
        high=8,
    )


# --------------------------------------------------------------------------
# input verification
# --------------------------------------------------------------------------


def verify_evidence_bundle(directory: Path) -> dict[str, Any]:
    """Verify one saved evidence bundle before any payload is read.

    ``manifest.json`` is checked against its own checksum sidecar and then
    every file it records is re-hashed and re-sized.  Any mismatch aborts the
    audit loudly: reading a payload whose digest does not match its manifest
    would make every number downstream unattributable.
    """

    if not directory.is_dir():
        raise ValueError(f"input evidence directory is missing: {directory}")
    manifest_path = directory / "manifest.json"
    try:
        verify_digest(manifest_path)
    except ValueError as error:
        raise ValueError(
            f"manifest integrity check failed for {directory.name}: {error}"
        ) from error
    manifest = strict_json_loads(
        manifest_path.read_text(encoding="utf-8"),
        label=f"{directory.name} manifest",
    )
    if not isinstance(manifest, dict) or set(manifest) != BUNDLE_MANIFEST_KEYS:
        raise ValueError(
            f"{directory.name}/manifest.json has missing or unknown keys"
        )
    if manifest["schema_version"] != 1:
        raise ValueError(f"{directory.name}/manifest.json schema_version must be 1")
    experiment_id = _nonempty_text(
        manifest["experiment_id"], name=f"{directory.name} manifest experiment_id"
    )
    files = manifest["files"]
    if not isinstance(files, dict) or not files:
        raise ValueError(f"{directory.name}/manifest.json lists no files")

    digests: dict[str, str] = {}
    for name in sorted(files):
        entry = files[name]
        if not isinstance(entry, dict) or set(entry) != {"sha256", "bytes"}:
            raise ValueError(
                f"{directory.name}/manifest.json entry {name!r} is malformed"
            )
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise ValueError(
                f"{directory.name}/manifest.json entry {name!r} escapes the bundle"
            )
        target = directory.joinpath(*relative.parts)
        if not target.is_file():
            raise ValueError(
                f"{directory.name}/manifest.json records a missing file: {name}"
            )
        size = target.stat().st_size
        digest = sha256_file(target)
        if digest != entry["sha256"] or size != entry["bytes"]:
            raise ValueError(
                f"integrity failure in {directory.name}: {name} does not match "
                "its recorded sha256/size"
            )
        digests[name] = digest

    checkpoints: dict[tuple[int, int, str], dict[str, Any]] = {}
    for name, digest in digests.items():
        match = _CHECKPOINT_NAME.fullmatch(name)
        if match is None:
            continue
        key = (int(match.group("nr")), int(match.group("nz")), match.group("label"))
        checkpoints[key] = {
            "name": name,
            "path": directory.joinpath(*PurePosixPath(name).parts),
            "sha256": digest,
        }
    if not checkpoints:
        raise ValueError(
            f"{directory.name} contains no checkpoint archives to audit"
        )
    return {
        "directory": _display_path(directory),
        "experiment_id": experiment_id,
        "manifest_sha256": sha256_file(manifest_path),
        "file_count": len(digests),
        "checkpoint_count": len(checkpoints),
        "verified": True,
        "checkpoints": checkpoints,
    }


# --------------------------------------------------------------------------
# continuum reference maximization
# --------------------------------------------------------------------------


def _e29_axial_profile(z: np.ndarray) -> np.ndarray:
    """``f(z) = sin(2 pi z) / (1 + 12.5 sin^2(pi z))``."""

    sine = np.sin(np.pi * z)
    return np.sin(2.0 * np.pi * z) / (
        1.0 + E29_AXIAL_DENOMINATOR_COEFFICIENT * sine * sine
    )


def _e29_axial_profile_derivative(z: np.ndarray) -> np.ndarray:
    r"""Closed-form ``f'(z)``.

    With ``s = sin(pi z)``, ``c = cos(pi z)`` and ``D = 1 + 12.5 s^2``,

    ``f' = 2 pi [cos(2 pi z) D - 25 s^2 c^2] / D^2``,

    obtained from the quotient rule using ``sin(2 pi z) = 2 s c``.  The
    derivative is analytic, not a finite difference, so the reported continuum
    maximum carries no differencing error.
    """

    sine = np.sin(np.pi * z)
    cosine = np.cos(np.pi * z)
    denominator = 1.0 + E29_AXIAL_DENOMINATOR_COEFFICIENT * sine * sine
    numerator = np.cos(2.0 * np.pi * z) * denominator - (
        2.0 * E29_AXIAL_DENOMINATOR_COEFFICIENT * sine * sine * cosine * cosine
    )
    return 2.0 * np.pi * numerator / (denominator * denominator)


def continuum_initial_vorticity_magnitude(
    r: np.ndarray,
    z: np.ndarray,
    *,
    amplitude: float = E29_INITIAL_AMPLITUDE,
) -> np.ndarray:
    r"""Cartesian vorticity magnitude of the closed-form E-29 initial datum.

    ``u1(0,r,z) = A (1-r^2)^18 f(z)`` with ``omega1(0) = 0``, so the E-18b
    components reduce to ``omega_r = -r u1_z``, ``omega_theta = 0`` and
    ``omega_z = 2 u1 + r u1_r``.  Both survive in closed form:

    * ``r u1_z = A r (1-r^2)^18 f'(z)``;
    * ``2 u1 + r u1_r = A f(z) (1-r^2)^17 (2 - 38 r^2)``, using
      ``r u1_r = -36 A r^2 (1-r^2)^17 f(z)``.

    The returned array has shape ``(r.size, z.size)``.
    """

    radius = np.asarray(r, dtype=np.float64)[:, None]
    axial = np.asarray(z, dtype=np.float64)[None, :]
    one_minus = 1.0 - radius * radius
    profile = _e29_axial_profile(axial)
    derivative = _e29_axial_profile_derivative(axial)
    radial_component = (
        amplitude * radius * one_minus**E29_RADIAL_EXPONENT * derivative
    )
    axial_component = (
        amplitude
        * profile
        * one_minus ** (E29_RADIAL_EXPONENT - 1)
        * (2.0 - 2.0 * (E29_RADIAL_EXPONENT + 1) * radius * radius)
    )
    return np.sqrt(radial_component**2 + axial_component**2)


def continuum_initial_vorticity_maximum(
    *,
    search_points_per_axis: int,
    refinement_stages: int,
    amplitude: float = E29_INITIAL_AMPLITUDE,
) -> dict[str, Any]:
    r"""Locate the continuum initial vorticity maximum by nested refinement.

    The search is a plain tensor scan of
    :func:`continuum_initial_vorticity_magnitude` over
    ``r in [0, 1]`` and ``z in [-1/2, 1/2]`` -- a full axial period placed so
    that the expected maximum at ``z = 0`` is interior rather than sitting on a
    window edge.  Each stage evaluates ``search_points_per_axis`` samples per
    axis, takes the discrete argmax and shrinks the window to the two cells on
    either side of it, clamped to ``r in [0, 1]``.  With ``N`` samples and ``S``
    stages the final window is ``(4/(N-1))^(S-1) / (N-1)`` wide, so the located
    coordinates are resolved far below the grid spacing of any audited
    snapshot.

    The result is compared against the closed form
    ``24000 pi 37^{-1/2} (36/37)^18``, which is derivable by hand: the radial
    term ``A r (1-r^2)^18 f'(z)`` is maximized at ``f'(0) = 2 pi`` and at
    ``r^2 = 1/37`` (from ``d/dr[r (1-r^2)^18] = (1-r^2)^17 (1 - 37 r^2)``),
    where the axial term vanishes because ``f(0) = 0``.  The agreement is
    *recorded*, not assumed; a search that disagreed with the closed form would
    show up in ``relative_difference_to_closed_form``.
    """

    points = _bounded_int(
        search_points_per_axis, name="search_points_per_axis", low=101, high=4001
    )
    stages = _bounded_int(
        refinement_stages, name="refinement_stages", low=1, high=8
    )
    r_low, r_high = 0.0, 1.0
    z_low, z_high = -0.5, 0.5
    history: list[dict[str, Any]] = []
    best_r = best_z = best_value = 0.0
    for stage in range(stages):
        radius = np.linspace(r_low, r_high, points, dtype=np.float64)
        axial = np.linspace(z_low, z_high, points, dtype=np.float64)
        magnitude = continuum_initial_vorticity_magnitude(
            radius, axial, amplitude=amplitude
        )
        index_r, index_z = np.unravel_index(
            int(np.argmax(magnitude)), magnitude.shape
        )
        best_r = float(radius[index_r])
        best_z = float(axial[index_z])
        best_value = float(magnitude[index_r, index_z])
        dr = float(radius[1] - radius[0])
        dz = float(axial[1] - axial[0])
        history.append(
            {
                "stage": int(stage),
                "r_window": [r_low, r_high],
                "z_window": [z_low, z_high],
                "dr": dr,
                "dz": dz,
                "r": best_r,
                "z": best_z,
                "value": best_value,
            }
        )
        r_low = max(0.0, best_r - 2.0 * dr)
        r_high = min(1.0, best_r + 2.0 * dr)
        z_low = best_z - 2.0 * dz
        z_high = best_z + 2.0 * dz

    closed_form = amplitude / E29_INITIAL_AMPLITUDE * (
        E29B_MAX_CARTESIAN_VORTICITY_EXACT
    )
    return {
        "amplitude": float(amplitude),
        "search_points_per_axis": points,
        "refinement_stages": stages,
        "search_domain": {"r": [0.0, 1.0], "z": [-0.5, 0.5]},
        "r": best_r,
        "z": best_z,
        "value": best_value,
        "closed_form_value": float(closed_form),
        "closed_form_expression": "24000*pi*37^(-1/2)*(36/37)^18",
        "closed_form_r": float(1.0 / math.sqrt(37.0)),
        "closed_form_z": 0.0,
        "relative_difference_to_closed_form": float(
            abs(best_value - closed_form) / closed_form
        ),
        "r_difference_to_closed_form": float(best_r - 1.0 / math.sqrt(37.0)),
        "z_difference_to_closed_form": float(best_z),
        "stages": history,
        "method": (
            "nested tensor-grid refinement of the closed-form E-29 Cartesian "
            "vorticity magnitude sqrt((r u1_z)^2 + (2 u1 + r u1_r)^2) with an "
            "analytic axial derivative; the window shrinks to two cells on "
            "either side of the discrete argmax at every stage"
        ),
    }


# --------------------------------------------------------------------------
# snapshot measurement
# --------------------------------------------------------------------------


def vorticity_magnitude(
    grid: AxisymmetricGrid,
    u1: np.ndarray,
    omega1: np.ndarray,
) -> np.ndarray:
    """Return the pointwise norm of the E-18b Cartesian vorticity."""

    components = cartesian_vorticity(grid, u1, omega1)
    total = np.zeros(grid.shape, dtype=np.float64)
    for component in components:
        total += component * component
    return np.sqrt(total)


def audited_fields(
    grid: AxisymmetricGrid,
    fields: Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Return the two audited arrays keyed by :data:`AUDITED_FIELDS`."""

    return {
        "u1": np.asarray(fields["u1"], dtype=np.float64),
        "cartesian_vorticity_magnitude": vorticity_magnitude(
            grid, fields["u1"], fields["omega1"]
        ),
    }


def _report_is_finite(report: Mapping[str, Any]) -> bool:
    stack: list[Any] = [report]
    while stack:
        value = stack.pop()
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float)):
            if not math.isfinite(float(value)):
                return False
        elif isinstance(value, dict):
            stack.extend(value.values())
        elif isinstance(value, (list, tuple)):
            stack.extend(value)
    return True


def measure_snapshot(
    grid: AxisymmetricGrid,
    fields: Mapping[str, np.ndarray],
) -> dict[str, dict[str, Any]]:
    """Run the full width report and the fit gate on both audited fields."""

    measurements: dict[str, dict[str, Any]] = {}
    arrays = audited_fields(grid, fields)
    for name in AUDITED_FIELDS:
        report = peak_resolution_report(grid, arrays[name])
        gate = fit_precondition(report)
        measurements[name] = {
            "report": report,
            "fit_precondition": gate,
            "report_is_finite": bool(_report_is_finite(report)),
        }
    return measurements


def _wrapped(delta: float, period: float) -> float:
    """Signed representative of ``delta`` in ``[-period/2, period/2)``."""

    return float(delta - period * round(delta / period))


# --------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------


def evaluate(config: Mapping[str, Any]) -> dict[str, Any]:
    """Verify every input bundle, audit every snapshot and assemble the tables."""

    validate_config(config)
    labels = [str(label) for label in config["snapshot_labels"]]
    initial_label = str(config["initial_snapshot_label"])
    final_label = str(config["final_snapshot_label"])

    bundles = [
        verify_evidence_bundle(evidence_directory(entry))
        for entry in config["input_evidence_dirs"]
    ]

    # ------------------------------------------------------------------
    # Resolve the resolution ladder.  A resolution appearing in more than
    # one bundle is audited once, after checking that every shared archive
    # is byte identical; two copies of the same bytes are not two samples.
    # ------------------------------------------------------------------
    owner: dict[tuple[int, int], dict[str, Any]] = {}
    duplicates: list[dict[str, Any]] = []
    for bundle in bundles:
        for (nr, nz, label), entry in bundle["checkpoints"].items():
            key = (nr, nz)
            record = owner.setdefault(
                key,
                {
                    "nr": nr,
                    "nz": nz,
                    "source_evidence_dir": bundle["directory"],
                    "entries": {},
                },
            )
            if record["source_evidence_dir"] == bundle["directory"]:
                record["entries"][label] = entry
                continue
            previous = record["entries"].get(label)
            identical = previous is not None and previous["sha256"] == entry["sha256"]
            duplicates.append(
                {
                    "nr": nr,
                    "nz": nz,
                    "snapshot_label": label,
                    "audited_source": record["source_evidence_dir"],
                    "duplicate_source": bundle["directory"],
                    "byte_identical": bool(identical),
                }
            )
            if not identical:
                raise ValueError(
                    f"resolution {nr}x{nz} snapshot {label} differs between "
                    f"{record['source_evidence_dir']} and {bundle['directory']}; "
                    "the audit refuses to pick one silently"
                )

    resolutions = sorted(owner)
    if not resolutions:
        raise ValueError("no resolutions were discovered in the audited bundles")

    # ------------------------------------------------------------------
    # Load and measure every requested snapshot.
    # ------------------------------------------------------------------
    snapshots: list[dict[str, Any]] = []
    arrays: dict[tuple[int, int, str], dict[str, Any]] = {}
    for key in resolutions:
        record = owner[key]
        for label in labels:
            entry = record["entries"].get(label)
            if entry is None:
                raise ValueError(
                    f"resolution {key[0]}x{key[1]} has no checkpoint for "
                    f"snapshot label {label}"
                )
            candidate = load_candidate(entry["path"])
            if candidate.description is None:
                raise ValueError(f"checkpoint lacks a v2 description: {entry['name']}")
            grid = candidate.grid
            measurements = measure_snapshot(grid, candidate.fields)
            arrays[(key[0], key[1], label)] = {
                "grid": grid,
                "fields": audited_fields(grid, candidate.fields),
            }
            snapshots.append(
                {
                    "source_evidence_dir": record["source_evidence_dir"],
                    "checkpoint": entry["name"],
                    "archive_sha256": entry["sha256"],
                    "nr": int(grid.nr),
                    "nz": int(grid.nz),
                    "dr": float(grid.dr),
                    "dz": float(grid.dz),
                    "snapshot_label": label,
                    "physical_time": float(candidate.description.physical_time),
                    "viscosity": float(candidate.description.viscosity),
                    "fields": measurements,
                }
            )

    # ------------------------------------------------------------------
    # Fit-precondition ladders.  A convergence fit is a statement about a
    # resolution SERIES at one time for one quantity, so the ladder-level
    # verdict is the conjunction over the whole series; a single resolution
    # passing the gate cannot license a fit on its own.
    # ------------------------------------------------------------------
    ladders: list[dict[str, Any]] = []
    for label in labels:
        for name in AUDITED_FIELDS:
            rows = [
                snapshot
                for snapshot in snapshots
                if snapshot["snapshot_label"] == label
            ]
            verdicts = [
                {
                    "nr": snapshot["nr"],
                    "nz": snapshot["nz"],
                    "satisfied": bool(
                        snapshot["fields"][name]["fit_precondition"]["satisfied"]
                    ),
                    "points_per_front": float(
                        snapshot["fields"][name]["report"]["points_per_front"]
                    ),
                    "points_per_fwhm_r": float(
                        snapshot["fields"][name]["report"]["points_per_fwhm_r"]
                    ),
                    "peak_to_axis_cells": int(
                        snapshot["fields"][name]["report"]["peak_to_axis_cells"]
                    ),
                }
                for snapshot in sorted(rows, key=lambda item: item["nr"])
            ]
            ladders.append(
                {
                    "snapshot_label": label,
                    "physical_time": float(rows[0]["physical_time"]),
                    "field": name,
                    "resolutions_in_ladder": len(verdicts),
                    "verdicts": verdicts,
                    "ladder_satisfied": bool(
                        len(verdicts) >= 2
                        and all(item["satisfied"] for item in verdicts)
                    ),
                }
            )

    satisfied_anywhere = any(ladder["ladder_satisfied"] for ladder in ladders)
    final_ladders = [
        ladder for ladder in ladders if ladder["snapshot_label"] == final_label
    ]
    satisfied_at_final = any(
        ladder["ladder_satisfied"] for ladder in final_ladders
    )

    # ------------------------------------------------------------------
    # Adjacent-resolution common-grid differences.
    # ------------------------------------------------------------------
    differences: list[dict[str, Any]] = []
    expected_pairs = 0
    for label in labels:
        present = [key for key in resolutions if (key[0], key[1], label) in arrays]
        for coarse_key, fine_key in zip(present, present[1:]):
            expected_pairs += len(AUDITED_FIELDS)
            coarse = arrays[(coarse_key[0], coarse_key[1], label)]
            fine = arrays[(fine_key[0], fine_key[1], label)]
            physical_time = next(
                snapshot["physical_time"]
                for snapshot in snapshots
                if snapshot["snapshot_label"] == label
            )
            for name in AUDITED_FIELDS:
                result = interpolate_to_common_grid(
                    coarse["grid"],
                    coarse["fields"][name],
                    fine["grid"],
                    fine["fields"][name],
                )
                scale = max(
                    float(np.max(np.abs(coarse["fields"][name]))),
                    float(np.max(np.abs(fine["fields"][name]))),
                )
                differences.append(
                    {
                        "snapshot_label": label,
                        "physical_time": float(physical_time),
                        "field": name,
                        "coarse_nr": coarse_key[0],
                        "coarse_nz": coarse_key[1],
                        "fine_nr": fine_key[0],
                        "fine_nz": fine_key[1],
                        "linf": float(result["linf"]),
                        "weighted_l2": float(result["weighted_l2"]),
                        "derivative_linf": float(result["derivative_linf"]),
                        "derivative_linf_r": float(result["derivative_linf_r"]),
                        "derivative_linf_z": float(result["derivative_linf_z"]),
                        "argmax_r": float(result["argmax_r"]),
                        "argmax_z": float(result["argmax_z"]),
                        "field_max_abs": scale,
                        "linf_over_field_max_abs": (
                            float(result["linf"] / scale) if scale > 0.0 else None
                        ),
                        "common_shape": list(result["common_shape"]),
                        "interpolation": str(result["interpolation"]),
                    }
                )

    # ------------------------------------------------------------------
    # P1-B amplification renormalization.
    # ------------------------------------------------------------------
    continuum = continuum_initial_vorticity_maximum(
        search_points_per_axis=int(
            config["continuum_reference"]["search_points_per_axis"]
        ),
        refinement_stages=int(config["continuum_reference"]["refinement_stages"]),
    )
    amplification: list[dict[str, Any]] = []
    for key in resolutions:
        initial = arrays[(key[0], key[1], initial_label)]
        final = arrays[(key[0], key[1], final_label)]
        grid = initial["grid"]
        initial_field = initial["fields"]["cartesian_vorticity_magnitude"]
        final_field = final["fields"]["cartesian_vorticity_magnitude"]
        index_r, index_z = np.unravel_index(
            int(np.argmax(initial_field)), grid.shape
        )
        discrete_initial = float(initial_field[index_r, index_z])
        if not discrete_initial > 0.0:
            raise ValueError(
                f"resolution {key[0]}x{key[1]} has a vanishing initial maximum "
                "Cartesian vorticity; no grid normalization is defined"
            )
        final_maximum = float(np.max(final_field))
        location_r = float(grid.r[index_r])
        location_z = float(grid.z[index_z])
        delta_r = location_r - float(continuum["r"])
        delta_z = _wrapped(location_z - float(continuum["z"]), float(grid.z_period))
        initial_time = next(
            snapshot["physical_time"]
            for snapshot in snapshots
            if snapshot["nr"] == key[0] and snapshot["snapshot_label"] == initial_label
        )
        final_time = next(
            snapshot["physical_time"]
            for snapshot in snapshots
            if snapshot["nr"] == key[0] and snapshot["snapshot_label"] == final_label
        )
        amplification.append(
            {
                "nr": key[0],
                "nz": key[1],
                "dr": float(grid.dr),
                "dz": float(grid.dz),
                "source_evidence_dir": owner[key]["source_evidence_dir"],
                "initial_snapshot_label": initial_label,
                "initial_physical_time": float(initial_time),
                "final_snapshot_label": final_label,
                "final_physical_time": float(final_time),
                "initial_max_cartesian_vorticity_discrete": discrete_initial,
                "common_continuum_reference": float(
                    E29B_MAX_CARTESIAN_VORTICITY_EXACT
                ),
                "final_max_cartesian_vorticity": final_maximum,
                "amplification_grid_normalized": float(
                    final_maximum / discrete_initial
                ),
                "amplification_common_reference": float(
                    final_maximum / E29B_MAX_CARTESIAN_VORTICITY_EXACT
                ),
                "discrete_initial_over_continuum_reference": float(
                    discrete_initial / E29B_MAX_CARTESIAN_VORTICITY_EXACT
                ),
                "initial_max_location": {
                    "r": location_r,
                    "z": location_z,
                    "index_r": int(index_r),
                    "index_z": int(index_z),
                },
                "continuum_initial_max_location": {
                    "r": float(continuum["r"]),
                    "z": float(continuum["z"]),
                },
                "initial_max_location_offset": {
                    "delta_r": float(delta_r),
                    "delta_z_wrapped": float(delta_z),
                    "distance": float(math.hypot(delta_r, delta_z)),
                    "delta_r_in_cells": float(delta_r / grid.dr),
                    "delta_z_in_cells": float(delta_z / grid.dz),
                },
            }
        )

    checks = {
        "all_input_manifests_verified": bool(
            len(bundles) == len(config["input_evidence_dirs"])
            and all(bundle["verified"] for bundle in bundles)
        ),
        "all_requested_snapshots_reported_finite": bool(
            len(snapshots) == len(resolutions) * len(labels)
            and all(
                snapshot["fields"][name]["report_is_finite"]
                for snapshot in snapshots
                for name in AUDITED_FIELDS
            )
        ),
        "final_snapshot_fit_precondition_recorded_for_every_resolution": bool(
            sorted(
                (snapshot["nr"], snapshot["nz"])
                for snapshot in snapshots
                if snapshot["snapshot_label"] == final_label
                and all(
                    isinstance(
                        snapshot["fields"][name]["fit_precondition"]["satisfied"],
                        bool,
                    )
                    for name in AUDITED_FIELDS
                )
            )
            == resolutions
        ),
        "convergence_fit_precondition_recorded": bool(
            isinstance(satisfied_anywhere, bool)
            and satisfied_anywhere
            == any(ladder["ladder_satisfied"] for ladder in ladders)
        ),
        "adjacent_resolution_differences_complete": bool(
            len(differences) == expected_pairs and expected_pairs > 0
        ),
    }

    return {
        "bundles": [
            {key: value for key, value in bundle.items() if key != "checkpoints"}
            for bundle in bundles
        ],
        "duplicate_resolution_checks": duplicates,
        "resolutions": [{"nr": nr, "nz": nz} for nr, nz in resolutions],
        "snapshots": snapshots,
        "fit_precondition_ladders": ladders,
        "convergence_fit_precondition_satisfied_anywhere": bool(satisfied_anywhere),
        "convergence_fit_precondition_satisfied_at_final_snapshot": bool(
            satisfied_at_final
        ),
        "common_grid_differences": differences,
        "continuum_initial_maximum": continuum,
        "amplification_renormalization": amplification,
        "acceptance_checks": checks,
    }


# --------------------------------------------------------------------------
# CSV tables
# --------------------------------------------------------------------------


def _report_rows(snapshots: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for snapshot in snapshots:
        for name in AUDITED_FIELDS:
            report = snapshot["fields"][name]["report"]
            flags = report["truncation_flags"]
            gradient = report["gradient_length_scale"]
            tail = report["high_frequency_tail"]
            rows.append(
                {
                    "source_evidence_dir": snapshot["source_evidence_dir"],
                    "nr": snapshot["nr"],
                    "nz": snapshot["nz"],
                    "dr": snapshot["dr"],
                    "dz": snapshot["dz"],
                    "snapshot_label": snapshot["snapshot_label"],
                    "physical_time": snapshot["physical_time"],
                    "field": name,
                    "peak_r": report["peak"]["r"],
                    "peak_z": report["peak"]["z"],
                    "peak_index_r": report["peak"]["index_r"],
                    "peak_index_z": report["peak"]["index_z"],
                    "peak_magnitude": report["peak"]["magnitude"],
                    "subgrid_peak_r": report["subgrid_peak_r"]["coordinate"],
                    "subgrid_peak_r_value": report["subgrid_peak_r"]["value"],
                    "subgrid_peak_r_offset_cells": report["subgrid_peak_r"][
                        "offset_cells"
                    ],
                    "subgrid_peak_r_refined": report["subgrid_peak_r"]["refined"],
                    "subgrid_peak_z": report["subgrid_peak_z"]["coordinate"],
                    "subgrid_peak_z_value": report["subgrid_peak_z"]["value"],
                    "subgrid_peak_z_offset_cells": report["subgrid_peak_z"][
                        "offset_cells"
                    ],
                    "subgrid_peak_z_refined": report["subgrid_peak_z"]["refined"],
                    "fwhm_r": report["fwhm_r"]["fwhm"],
                    "points_per_fwhm_r": report["points_per_fwhm_r"],
                    "fwhm_z": report["fwhm_z"]["fwhm"],
                    "points_per_fwhm_z": report["points_per_fwhm_z"],
                    "front_thickness_10_90": report["front_to_axis"]["thickness"],
                    "points_per_front": report["points_per_front"],
                    "peak_to_axis_cells": report["peak_to_axis_cells"],
                    "min_gradient_length_scale": gradient["min_length_scale"],
                    "points_per_gradient_scale_r": gradient[
                        "points_per_gradient_scale_r"
                    ],
                    "points_per_gradient_scale_z": gradient[
                        "points_per_gradient_scale_z"
                    ],
                    "gradient_scale_r": gradient["r"],
                    "gradient_scale_z": gradient["z"],
                    "z_tail_fraction": tail["z_tail_fraction"],
                    "z_tail_fraction_r": tail["z_tail_fraction_r"],
                    "r_tail_fraction": tail["r_tail_fraction"],
                    "any_truncation": report["any_truncation"],
                    "truncated_fwhm_r_left": flags["fwhm_r_left"],
                    "truncated_fwhm_r_right": flags["fwhm_r_right"],
                    "truncated_fwhm_z_left": flags["fwhm_z_left"],
                    "truncated_fwhm_z_right": flags["fwhm_z_right"],
                    "truncated_front_90": flags["front_90"],
                    "truncated_front_10": flags["front_10"],
                }
            )
    return rows


def _fit_rows(snapshots: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for snapshot in snapshots:
        for name in AUDITED_FIELDS:
            gate = snapshot["fields"][name]["fit_precondition"]
            rows.append(
                {
                    "source_evidence_dir": snapshot["source_evidence_dir"],
                    "nr": snapshot["nr"],
                    "nz": snapshot["nz"],
                    "snapshot_label": snapshot["snapshot_label"],
                    "physical_time": snapshot["physical_time"],
                    "field": name,
                    "min_points": gate["min_points"],
                    "points_per_front": gate["points_per_front"],
                    "points_per_fwhm_r": gate["points_per_fwhm_r"],
                    "peak_to_axis_cells": gate["peak_to_axis_cells"],
                    "any_truncation": gate["any_truncation"],
                    "satisfied": gate["satisfied"],
                    "reason_count": len(gate["reasons"]),
                    "reasons": " | ".join(gate["reasons"]),
                }
            )
    return rows


def _difference_rows(
    differences: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            "snapshot_label": entry["snapshot_label"],
            "physical_time": entry["physical_time"],
            "field": entry["field"],
            "coarse_nr": entry["coarse_nr"],
            "coarse_nz": entry["coarse_nz"],
            "fine_nr": entry["fine_nr"],
            "fine_nz": entry["fine_nz"],
            "linf": entry["linf"],
            "weighted_l2": entry["weighted_l2"],
            "derivative_linf": entry["derivative_linf"],
            "derivative_linf_r": entry["derivative_linf_r"],
            "derivative_linf_z": entry["derivative_linf_z"],
            "argmax_r": entry["argmax_r"],
            "argmax_z": entry["argmax_z"],
            "field_max_abs": entry["field_max_abs"],
            "linf_over_field_max_abs": (
                "" if entry["linf_over_field_max_abs"] is None
                else entry["linf_over_field_max_abs"]
            ),
        }
        for entry in differences
    ]


def _amplification_rows(
    amplification: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            "nr": entry["nr"],
            "nz": entry["nz"],
            "source_evidence_dir": entry["source_evidence_dir"],
            "initial_snapshot_label": entry["initial_snapshot_label"],
            "initial_physical_time": entry["initial_physical_time"],
            "final_snapshot_label": entry["final_snapshot_label"],
            "final_physical_time": entry["final_physical_time"],
            "initial_max_cartesian_vorticity_discrete": entry[
                "initial_max_cartesian_vorticity_discrete"
            ],
            "common_continuum_reference": entry["common_continuum_reference"],
            "final_max_cartesian_vorticity": entry["final_max_cartesian_vorticity"],
            "amplification_grid_normalized": entry["amplification_grid_normalized"],
            "amplification_common_reference": entry[
                "amplification_common_reference"
            ],
            "discrete_initial_over_continuum_reference": entry[
                "discrete_initial_over_continuum_reference"
            ],
            "initial_max_r": entry["initial_max_location"]["r"],
            "initial_max_z": entry["initial_max_location"]["z"],
            "continuum_initial_max_r": entry["continuum_initial_max_location"]["r"],
            "continuum_initial_max_z": entry["continuum_initial_max_location"]["z"],
            "delta_r": entry["initial_max_location_offset"]["delta_r"],
            "delta_z_wrapped": entry["initial_max_location_offset"][
                "delta_z_wrapped"
            ],
            "distance": entry["initial_max_location_offset"]["distance"],
            "delta_r_in_cells": entry["initial_max_location_offset"][
                "delta_r_in_cells"
            ],
            "delta_z_in_cells": entry["initial_max_location_offset"][
                "delta_z_in_cells"
            ],
        }
        for entry in amplification
    ]


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path.name}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------
# evidence bundle
# --------------------------------------------------------------------------


def _prepare_output(output_dir: Path) -> None:
    if output_dir.exists():
        if not output_dir.is_dir():
            raise FileExistsError(f"output path is not a directory: {output_dir}")
        if any(output_dir.iterdir()):
            raise FileExistsError(
                f"refusing to overwrite nonempty output directory: {output_dir}"
            )
    else:
        output_dir.mkdir(parents=True)


def _inside_repository(path: Path) -> bool:
    try:
        path.resolve().relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def _refuse_audited_sources(config: Mapping[str, Any], output_dir: Path) -> None:
    resolved = output_dir.resolve()
    for entry in config["input_evidence_dirs"]:
        source = evidence_directory(entry)
        if resolved == source or source in resolved.parents:
            raise ValueError(
                "refusing to write inside an audited source evidence directory"
            )


def run(config: Mapping[str, Any], output_dir: Path) -> dict[str, Any]:
    """Audit every snapshot and write a checksummed evidence bundle."""

    validate_config(config)
    if not _inside_repository(output_dir):
        raise ValueError("output directory must remain inside this repository")
    _refuse_audited_sources(config, output_dir)
    provenance = collect_runtime_provenance()
    config_bytes = canonical_json_bytes(dict(config))
    _prepare_output(output_dir)

    payload = evaluate(config)
    checks = payload["acceptance_checks"]

    summary: dict[str, Any] = {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "experiment_id": config["experiment_id"],
        "description": config["description"],
        "interpretation": config["interpretation"],
        "amplification_normalization_statement": AMPLIFICATION_STATEMENT,
        "preregistration": {
            "min_points_per_front": int(PREREGISTERED_MIN_POINTS_PER_FRONT),
            "front_error_tolerance": float(PREREGISTERED_FRONT_ERROR_TOLERANCE),
            "front_phases": int(PREREGISTERED_FRONT_PHASES),
            "delta_over_dr_values": [
                float(value) for value in PREREGISTERED_DELTA_OVER_DR_VALUES
            ],
            "source": (
                "ns_certificate_lab.core_width; the threshold was fixed by "
                "manufactured tanh fronts before any Hou snapshot was measured "
                "and is not restated in this experiment's config"
            ),
        },
        "method": {
            "reports": (
                "ns_certificate_lab.core_width.peak_resolution_report applied "
                "to the stored u1 array and to the pointwise magnitude of the "
                "E-18b Cartesian vorticity rebuilt from (u1, omega1)"
            ),
            "gate": (
                "ns_certificate_lab.core_width.fit_precondition with the "
                "module's own PREREGISTERED_MIN_POINTS_PER_FRONT"
            ),
            "ladder_rule": (
                "a convergence fit is a statement about a resolution SERIES at "
                "one time for one quantity, so ladder_satisfied is the "
                "conjunction of the gate over every resolution in that series; "
                "one resolution passing on its own cannot license a fit"
            ),
            "common_grid": (
                "ns_certificate_lab.core_width.interpolate_to_common_grid "
                "restricts the finer field onto the coarser grid with a "
                "separable Catmull-Rom cubic; the coarse field passes through "
                "the same code path, where it is the identity"
            ),
            "amplification": (
                "P1-B: both normalizations are reported side by side for every "
                "resolution, the grid-normalized one with a resolution "
                "dependent denominator and the common-reference one with the "
                "single continuum value below"
            ),
            "continuum_reference": (
                "computed in code from the closed form "
                "24000*pi*37^(-1/2)*(36/37)^18; its (r,z) location is located "
                "numerically by nested grid refinement of the closed-form "
                "E-29 datum, not quoted"
            ),
        },
        "audited_fields": list(AUDITED_FIELDS),
        "input_bundles": payload["bundles"],
        "duplicate_resolution_checks": payload["duplicate_resolution_checks"],
        "resolutions": payload["resolutions"],
        "snapshots": payload["snapshots"],
        "fit_precondition_ladders": payload["fit_precondition_ladders"],
        "convergence_fit_precondition_satisfied_anywhere": payload[
            "convergence_fit_precondition_satisfied_anywhere"
        ],
        "convergence_fit_precondition_satisfied_at_final_snapshot": payload[
            "convergence_fit_precondition_satisfied_at_final_snapshot"
        ],
        "convergence_fit_precondition_definition": (
            "there exists at least one (snapshot label, audited field) whose "
            "complete resolution ladder satisfies fit_precondition at every "
            "resolution audited; see fit_precondition_ladders"
        ),
        "common_grid_differences": payload["common_grid_differences"],
        "continuum_initial_maximum": payload["continuum_initial_maximum"],
        "common_reference_normalization": {
            "value": float(E29B_MAX_CARTESIAN_VORTICITY_EXACT),
            "closed_form": "24000*pi*37^(-1/2)*(36/37)^18",
            "source": "E-29b, docs/equation_audit.md",
        },
        "amplification_renormalization": payload["amplification_renormalization"],
        "acceptance_checks": checks,
        "accepted_as_core_width_audit": bool(all(checks.values())),
        "limitations": list(LIMITATIONS),
        "reproducibility": {
            "config_sha256": sha256_bytes(config_bytes),
            "runtime_provenance": provenance,
            "command": (
                "python experiments/run_core_width_audit.py --config "
                "configs/core_width_audit.json --output-dir "
                "outputs/core_width_audit_v1"
            ),
        },
    }

    write_with_digest(output_dir / "config.snapshot.json", config_bytes)
    write_with_digest(output_dir / "summary.json", canonical_json_bytes(summary))
    _write_csv(
        output_dir / "peak_resolution_reports.csv", _report_rows(payload["snapshots"])
    )
    _write_csv(output_dir / "fit_precondition.csv", _fit_rows(payload["snapshots"]))
    _write_csv(
        output_dir / "common_grid_differences.csv",
        _difference_rows(payload["common_grid_differences"]),
    )
    _write_csv(
        output_dir / "amplification_renormalization.csv",
        _amplification_rows(payload["amplification_renormalization"]),
    )

    artifact_names = [
        "amplification_renormalization.csv",
        "common_grid_differences.csv",
        "config.snapshot.json",
        "config.snapshot.json.sha256",
        "fit_precondition.csv",
        "peak_resolution_reports.csv",
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
            bundle["directory"]: {
                "experiment_id": bundle["experiment_id"],
                "manifest_sha256": bundle["manifest_sha256"],
                "file_count": bundle["file_count"],
            }
            for bundle in payload["bundles"]
        },
        "audited_checkpoints": {
            snapshot["checkpoint"]: snapshot["archive_sha256"]
            for snapshot in payload["snapshots"]
        },
    }
    write_with_digest(output_dir / "manifest.json", canonical_json_bytes(manifest))
    return summary


# --------------------------------------------------------------------------
# synthetic bundle writer (test support)
# --------------------------------------------------------------------------


def synthetic_state(
    grid: AxisymmetricGrid,
    *,
    time: float,
    amplitude: float = 1.0,
):
    r"""Return a constrained state built from a smooth closed-form datum.

    ``u1 = A e^{-40 t} r^2 (1-r^2)^2 sin(2 pi z)`` peaks off the axis at
    ``r = 1/sqrt(2)`` and vanishes at the wall, and
    ``omega1 = A e^{-40 t} (1-r^2)^2 cos(2 pi z)`` is even in ``r``.  Neither
    solves anything; the pair exists so that a bundle with the real schema, the
    real checkpoint writer and a real elliptic solve can be produced for tests
    without shipping simulation data.
    """

    radius, axial = grid.mesh()
    decay = math.exp(-40.0 * float(time))
    scale = float(amplitude) * decay
    u1 = scale * radius**2 * (1.0 - radius**2) ** 2 * np.sin(2.0 * np.pi * axial)
    omega1 = scale * (1.0 - radius**2) ** 2 * np.cos(2.0 * np.pi * axial)
    return constrain_state(grid, u1, omega1)


def write_minimal_evidence_bundle(
    directory: Path,
    *,
    experiment_id: str,
    resolutions: Sequence[tuple[int, int]],
    snapshot_times: Sequence[float],
    amplitude: float = 1.0,
    viscosity: float = 5.0e-4,
    seed: int = 0,
) -> Path:
    """Write a minimal evidence bundle in exactly the audited layout.

    The checkpoints go through :func:`nonlinear_cylinder.save_checkpoint`, so
    they are genuine schema-v2 candidates, and ``manifest.json`` is written with
    the same per-file ``sha256``/``bytes`` shape and the same
    :func:`write_with_digest` call that the shipped experiments use.  This is
    the fixture :mod:`tests.test_core_width_audit` audits; keeping it here
    guarantees the fixture cannot drift away from the format the verifier
    demands.
    """

    directory = Path(directory)
    _prepare_output(directory)
    checkpoint_dir = directory / "checkpoints"
    checkpoint_dir.mkdir()
    provenance = collect_runtime_provenance()
    for nr, nz in resolutions:
        grid = AxisymmetricGrid.uniform(
            nr=int(nr), nz=int(nz), r_max=1.0, z_min=0.0, z_max=1.0, periodic_z=True
        )
        for index, time in enumerate(snapshot_times):
            state = synthetic_state(grid, time=float(time), amplitude=amplitude)
            save_checkpoint(
                checkpoint_dir / f"checkpoint_nr{nr}_nz{nz}_t{index:03d}.npz",
                grid=grid,
                state=state,
                time=float(time),
                viscosity=float(viscosity),
                seed=int(seed),
                config={"experiment_id": experiment_id},
                provenance=provenance,
                metadata={"experiment_id": experiment_id, "synthetic": True},
            )
    names = sorted(
        path.relative_to(directory).as_posix()
        for path in directory.rglob("*")
        if path.is_file()
    )
    manifest = {
        "schema_version": 1,
        "experiment_id": experiment_id,
        "files": {
            name: {
                "sha256": sha256_file(directory.joinpath(*PurePosixPath(name).parts)),
                "bytes": directory.joinpath(
                    *PurePosixPath(name).parts
                ).stat().st_size,
            }
            for name in names
        },
    }
    write_with_digest(directory / "manifest.json", canonical_json_bytes(manifest))
    return directory


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "core_width_audit.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "core_width_audit_v1",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    try:
        config = strict_json_loads(
            config_path.read_text(encoding="utf-8"),
            label="core width audit config",
        )
        if not isinstance(config, dict):
            raise ValueError("config must be a JSON object")
        summary = run(config, output_dir)
    except (OSError, ValueError, FileExistsError) as exc:
        parser.error(str(exc))
    print(
        canonical_json_bytes(
            {
                "acceptance_checks": summary["acceptance_checks"],
                "convergence_fit_precondition_satisfied_anywhere": summary[
                    "convergence_fit_precondition_satisfied_anywhere"
                ],
                "convergence_fit_precondition_satisfied_at_final_snapshot": summary[
                    "convergence_fit_precondition_satisfied_at_final_snapshot"
                ],
            }
        ).decode("utf-8"),
        end="",
    )
    return 0 if summary["accepted_as_core_width_audit"] else 2


if __name__ == "__main__":
    sys.exit(main())
