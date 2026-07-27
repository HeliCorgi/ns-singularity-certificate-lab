#!/usr/bin/env python
r"""Independent uniform-Cartesian audit of saved Hou finite-cylinder snapshots.

Each audited artifact is a schema-v2 candidate written by
``experiments/run_hou_early_time.py``.  It stores only the transformed
axisymmetric nodal arrays ``(u1, omega1, psi1)`` on the cylindrical
half-plane.  This experiment reloads such a checkpoint with checksum
verification, maps it onto a uniform Cartesian box with
:mod:`ns_certificate_lab.cartesian_candidate_adapter`, and then evaluates it
with the stencils owned by :mod:`ns_certificate_lab.cartesian_validation`
alone.  Neither module imports ``operators.py``, ``pde.py`` or any cylindrical
finite-difference result, so the three quantities reported here are an
independent re-derivation rather than a restatement of the producing solver.

What is well defined from a *single* snapshot, and is therefore evaluated:

* **(a) E-02 divergence.**  The E-18a component map makes
  :math:`\nabla\cdot u\equiv 0` an identity for *every* pair
  ``(u1, psi1)``, so the measured Cartesian divergence is purely the
  interpolation plus truncation error of this audit.
* **(b) E-18b full curl.**  ``curl u`` recomputed on the Cartesian box is
  compared against the complete E-18b vorticity oracle
  ``omega=(-x u1_z - y omega1, -y u1_z + x omega1, 2 u1 + r u1_r)``
  reconstructed from the same stored arrays.  Its ``x`` and ``y`` components
  encode the elliptic relation ``-L5 psi1 = omega1`` and its ``z`` component
  encodes the swirl relation, so this is a real cross-check, not an identity.
* **(c) Cylindrical agreement.**  The azimuthal projection
  ``r*dot(curl u, e_theta) = -y (curl u)_x + x (curl u)_y`` is compared against
  ``r^2 omega1``, which is recovered from the adapter's own E-18b output as
  ``-y omega_x + x omega_y``.  This is the ``omega^theta = r*omega1`` half of
  E-18b in an axis-safe, division-free form.
* **(d) Relative versions** of all of the above, with the denominators
  documented in :data:`DENOMINATORS` and reported alongside every ratio.

What is *not* well defined from a single snapshot, and is therefore **not**
faked here: the primitive momentum residual
``u_t + (u.grad)u + grad p - nu*Laplacian(u)`` needs a pressure field and a
time derivative, and a checkpoint stores neither.
:func:`primitive_residual_from_snapshot_pair` assembles exactly that residual
from a closely spaced pair of Cartesian velocity fields and an independently
supplied pressure, using the existing ``cartesian_validation`` term machinery.
It is unit-tested on a manufactured analytic field whose pressure and time
derivative are known in closed form.  Applying it to the real Hou checkpoints
remains an explicit gap, recorded as :data:`PRIMITIVE_RESIDUAL_GAP` in every
summary rather than silently dropped.

This experiment reads ``outputs/`` and never writes into it.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
import sys
from typing import Any, Iterable, Mapping

import numpy as np
import numpy.typing as npt

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    write_with_digest,
)
from ns_certificate_lab.artifacts import CandidateArtifact, load_candidate
from ns_certificate_lab.cartesian_candidate_adapter import (
    CandidateCartesianArrays,
    reconstruct_loaded_candidate_on_cartesian,
)
from ns_certificate_lab.cartesian_validation import (
    PrimitiveNSResidual,
    UniformCartesianGrid,
    audit_cartesian_reconstruction,
    cartesian_curl,
    cartesian_divergence,
    cartesian_gradient,
    primitive_ns_residual,
)
from ns_certificate_lab.provenance import collect_runtime_provenance

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPECTED_EXPERIMENT_ID = "hou_snapshot_cartesian_audit_v1"
EXPECTED_COORDINATE_SYSTEM = (
    "uniform Cartesian (x,y,z) box [-a,a]^2 x [z_min,z_max) inscribed in the "
    "finite cylinder of the source checkpoints"
)
SUMMARY_SCHEMA_VERSION = 1
MANIFEST_SCHEMA = "ns-certificate-lab/hou-snapshot-cartesian-audit-manifest/v1"

#: Every relative metric divides an absolute defect by one of these documented
#: scales.  All of them are global maxima over the *whole* audit box (never
#: over the trimmed interior), so a single denominator normalizes both the
#: full-box and the interior-margin norms of one checkpoint.
DENOMINATORS: dict[str, str] = {
    "velocity_gradient_row_sum_max": (
        "max over grid points of max_i sum_j |d_j u_i|, the largest absolute "
        "row sum of the reconstructed Cartesian velocity gradient.  This is "
        "the audit's default cancellation denominator: every quantity here is "
        "a signed combination of first derivatives of u, so no such "
        "combination can be resolved below this scale times the relative "
        "rounding and truncation level."
    ),
    "divergence_term_sum_max": (
        "max over grid points of |d_x u_x| + |d_y u_y| + |d_z u_z|, i.e. the "
        "sum of the magnitudes of exactly the three terms that cancel in "
        "div u.  It is the tightest honest denominator for the divergence "
        "identity and is never larger than the row-sum scale."
    ),
    "curl_term_sum_max": (
        "max over grid points and over the three curl components of "
        "|d_j u_k| + |d_k u_j| for the cyclic pair (j,k) entering that "
        "component, i.e. the magnitudes of exactly the two terms differenced "
        "in each curl component."
    ),
    "vorticity_e18b_max_abs": (
        "max over grid points and components of |omega_E18b|.  Signal scale, "
        "not a cancellation scale: it says how large the compared quantity is "
        "rather than how much cancellation the comparison required."
    ),
    "azimuthal_cancellation_max": (
        "max over grid points of |y|*(|(curl u)_x| + |omega_x|) + "
        "|x|*(|(curl u)_y| + |omega_y|), the cancellation scale of the "
        "r-weighted azimuthal projection -y*(.)_x + x*(.)_y applied to both "
        "sides of the comparison."
    ),
    "azimuthal_signal_max": (
        "max over grid points of |r^2 * omega1| recovered as "
        "|-y omega_x + x omega_y| from the adapter's E-18b output.  This is "
        "identically zero whenever omega1 vanishes (as it does in the E-29 "
        "initial datum), so ratios against it are reported as null and "
        "flagged degenerate instead of dividing by zero."
    ),
    "velocity_max_abs": (
        "max over grid points and components of |u|.  Reported for scale "
        "context only; no relative metric divides by it."
    ),
}

LIMITATIONS: tuple[str, ...] = (
    "The adapter samples the cylindrical nodal arrays with BILINEAR "
    "interpolation.  Bilinear sampling of a field whose radial structure is "
    "(1-r^2)^18 with amplitude 12000 is only first-order accurate in the "
    "sampled values near the axis, and the audit's own second-order Cartesian "
    "stencils then differentiate that interpolant.  Every defect reported "
    "here therefore includes an interpolation error that no part of this "
    "experiment separates from a genuine inconsistency of the snapshot.",
    "The uniform Cartesian audit box is generally COARSER than the source "
    "cylindrical grid (for example 64 points across 1.4 in x versus 193 "
    "points across 1.0 in r), so the audit is deliberately re-examined on two "
    "audit resolutions and the change is recorded.  Defects that shrink under "
    "audit-grid refinement are audit truncation, not snapshot defects; that "
    "distinction is recorded, not assumed.",
    "The Cartesian divergence is an identity of the E-18a component map for "
    "every (u1, psi1) whatsoever.  Passing the divergence gate therefore "
    "certifies the reconstruction arithmetic and the audit resolution; it "
    "certifies nothing about whether the snapshot solves any equation.",
    "Pointwise maxima of the defects are dominated by a few near-axis points "
    "and do not converge monotonically under audit-grid refinement; only the "
    "RMS norms do.  The hard gates are therefore RMS gates and the maxima are "
    "recorded without gating.",
    "Floating-point binary64 arithmetic with tolerances only.  There is no "
    "interval enclosure, no discretization error bound, and nothing here is "
    "evidence for or against singularity formation.",
)

#: The one thing a single checkpoint cannot support.  Recorded verbatim in
#: every summary so that the omission is visible rather than inferred.
PRIMITIVE_RESIDUAL_GAP: dict[str, Any] = {
    "quantity": "u_t + (u.grad)u + grad p - nu*Laplacian(u)",
    "status": "NOT evaluated on any real Hou checkpoint by this experiment",
    "implemented": (
        "primitive_residual_from_snapshot_pair() assembles the residual at "
        "the temporal midpoint of a checkpoint pair using the existing "
        "cartesian_validation.primitive_ns_residual term machinery"
    ),
    "verified_on": (
        "a manufactured analytic Cartesian field whose pressure, time "
        "derivative and body force are known in closed form; the pair-based "
        "residual reproduces the exact-midpoint residual to second order in "
        "the checkpoint spacing (tests/test_hou_snapshot_cartesian_audit.py)"
    ),
    "missing_inputs": [
        "an independent pressure solve: a checkpoint stores (u1, omega1, "
        "psi1) only, and the kinematic pressure is never written, so grad p "
        "is unavailable at any snapshot time",
        "an adjacent-in-time checkpoint pair whose spacing is small compared "
        "with the integrator step: the shipped snapshot spacing is 5e-4 while "
        "the integrator cap is 1e-6, about 500 steps apart, so a two-point "
        "difference over that interval is not a controlled approximation of "
        "u_t for this amplitude",
    ],
    "consequence": (
        "the primitive momentum residual of the real snapshots is an open "
        "gap; no number in this summary should be read as bounding it"
    ),
}

_SQRT_TWO = math.sqrt(2.0)
_TINY = float(np.finfo(np.float64).tiny)


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


def _finite_number(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _shape(value: Any, *, name: str) -> tuple[int, int, int]:
    if (
        not isinstance(value, list)
        or len(value) != 3
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
    ):
        raise ValueError(f"{name} must be a list of three integers")
    if any(item < 4 for item in value):
        raise ValueError(f"{name} entries must each be at least 4")
    return (int(value[0]), int(value[1]), int(value[2]))


def _checkpoint_path(entry: Any) -> Path:
    text = _nonempty_text(entry, name="checkpoint_paths entry")
    candidate = Path(text)
    # Reject anything absolute under either flavour of path syntax so that a
    # config written on one platform cannot escape on another.
    rooted = (
        candidate.is_absolute()
        or bool(candidate.anchor)
        or PurePosixPath(text).is_absolute()
        or PureWindowsPath(text).is_absolute()
        or bool(PureWindowsPath(text).anchor)
    )
    if rooted or ".." in candidate.parts:
        raise ValueError(
            "checkpoint_paths entries must be repository-relative paths "
            "without parent traversal"
        )
    if candidate.suffix.lower() != ".npz":
        raise ValueError("checkpoint_paths entries must name .npz candidates")
    resolved = (REPOSITORY_ROOT / candidate).resolve()
    try:
        resolved.relative_to(REPOSITORY_ROOT)
    except ValueError as error:
        raise ValueError(
            "checkpoint_paths entries must resolve inside this repository"
        ) from error
    return resolved


def validate_config(config: Mapping[str, Any]) -> None:
    """Reject any config with missing, unknown or out-of-range entries."""

    required = {
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "coordinate_system",
        "source_experiment",
        "checkpoint_paths",
        "cartesian_grid",
        "hard_gates",
        "record_only_thresholds",
    }
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError(
            "hou snapshot cartesian audit config has missing or unknown keys"
        )
    if config["schema_version"] != SUMMARY_SCHEMA_VERSION:
        raise ValueError("schema_version must be 1")
    if config["experiment_id"] != EXPECTED_EXPERIMENT_ID:
        raise ValueError(f"experiment_id must be {EXPECTED_EXPERIMENT_ID}")
    for key in ("description", "interpretation", "source_experiment"):
        _nonempty_text(config[key], name=key)
    if config["coordinate_system"] != EXPECTED_COORDINATE_SYSTEM:
        raise ValueError("coordinate_system must equal the audited canonical value")

    paths = config["checkpoint_paths"]
    if not isinstance(paths, list) or not paths:
        raise ValueError("checkpoint_paths must be a nonempty list")
    if len(set(map(repr, paths))) != len(paths):
        raise ValueError("checkpoint_paths entries must be unique")
    for entry in paths:
        _checkpoint_path(entry)

    grid_spec = config["cartesian_grid"]
    grid_keys = {
        "half_width",
        "z_min",
        "z_max",
        "z_periodic",
        "shape",
        "coarse_shape",
        "interior_margin",
    }
    if not isinstance(grid_spec, dict) or set(grid_spec) != grid_keys:
        raise ValueError("cartesian_grid has missing or unknown keys")
    half_width = _finite_positive(grid_spec["half_width"], name="half_width")
    if half_width > 1.0:
        raise ValueError(
            "half_width above 1.0 cannot be inscribed in any unit-radius "
            "source cylinder"
        )
    z_min = _finite_number(grid_spec["z_min"], name="z_min")
    z_max = _finite_number(grid_spec["z_max"], name="z_max")
    if z_max <= z_min:
        raise ValueError("z_max must exceed z_min")
    if not isinstance(grid_spec["z_periodic"], bool):
        raise ValueError("z_periodic must be a boolean")
    fine = _shape(grid_spec["shape"], name="shape")
    coarse = _shape(grid_spec["coarse_shape"], name="coarse_shape")
    if any(c >= f for c, f in zip(coarse, fine)):
        raise ValueError(
            "coarse_shape must be strictly coarser than shape on every axis"
        )
    margin = grid_spec["interior_margin"]
    if isinstance(margin, bool) or not isinstance(margin, int) or margin < 0:
        raise ValueError("interior_margin must be a nonnegative integer")
    if any(size <= 2 * margin for size in (fine[0], fine[1], coarse[0], coarse[1])):
        raise ValueError("interior_margin leaves an empty norm region")

    gates = config["hard_gates"]
    expected_gates = {
        "maximum_relative_divergence_rms",
        "maximum_relative_curl_defect_rms",
    }
    if not isinstance(gates, dict) or set(gates) != expected_gates:
        raise ValueError("hard_gates has missing or unknown keys")
    for key in expected_gates:
        _finite_positive(gates[key], name=f"hard_gates.{key}")

    thresholds = config["record_only_thresholds"]
    expected_thresholds = {
        "reference_relative_divergence_max",
        "reference_relative_curl_defect_max",
        "reference_relative_azimuthal_defect_max",
    }
    if not isinstance(thresholds, dict) or set(thresholds) != expected_thresholds:
        raise ValueError("record_only_thresholds has missing or unknown keys")
    for key in expected_thresholds:
        _finite_positive(thresholds[key], name=f"record_only_thresholds.{key}")


# --------------------------------------------------------------------------
# Cartesian grid construction
# --------------------------------------------------------------------------


def build_cartesian_grid(
    grid_spec: Mapping[str, Any],
    shape: tuple[int, int, int],
) -> UniformCartesianGrid:
    """Build the audit box ``[-a,a]^2 x [z_min, z_max)`` at ``shape``."""

    half_width = float(grid_spec["half_width"])
    return UniformCartesianGrid.uniform(
        shape=shape,
        x_bounds=(-half_width, half_width),
        y_bounds=(-half_width, half_width),
        z_bounds=(float(grid_spec["z_min"]), float(grid_spec["z_max"])),
        periodic=(False, False, bool(grid_spec["z_periodic"])),
    )


def check_box_is_inscribed(
    candidate: CandidateArtifact,
    grid_spec: Mapping[str, Any],
) -> None:
    """Refuse an audit box whose corners leave the source cylinder."""

    half_width = float(grid_spec["half_width"])
    r_max = float(candidate.grid.r[-1])
    corner_radius = half_width * _SQRT_TWO
    if corner_radius > r_max * (1.0 + 1.0e-12):
        raise ValueError(
            "Cartesian audit box corners leave the source cylinder: "
            f"half_width*sqrt(2)={corner_radius!r} exceeds r_max={r_max!r}"
        )


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------


def _rms(values: npt.NDArray[np.float64]) -> float:
    array = np.asarray(values, dtype=np.float64)
    return float(np.sqrt(np.mean(array * array)))


def _ratio(numerator: float, denominator: float) -> float | None:
    if not math.isfinite(denominator) or denominator <= _TINY:
        return None
    value = numerator / denominator
    return float(value) if math.isfinite(value) else None


def _location(
    field: npt.NDArray[np.float64],
    mesh: tuple[
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ],
    *,
    offset: tuple[int, int, int],
    vector: bool,
) -> dict[str, Any]:
    """Locate the maximum magnitude of ``field`` in whole-box coordinates."""

    index = np.unravel_index(int(np.argmax(np.abs(field))), field.shape)
    component = int(index[0]) if vector else None
    spatial = tuple(int(value) for value in (index[1:] if vector else index))
    absolute = tuple(int(value) + int(shift) for value, shift in zip(spatial, offset))
    x, y, z = mesh
    return {
        "component": component,
        "index": list(absolute),
        "x": float(x[absolute]),
        "y": float(y[absolute]),
        "z": float(z[absolute]),
        "radius": float(math.hypot(float(x[absolute]), float(y[absolute]))),
        "value": float(field[index]),
    }


def audit_cartesian_arrays(
    grid: UniformCartesianGrid,
    *,
    velocity: npt.ArrayLike,
    vorticity_e18b: npt.ArrayLike,
    interior_margin: int,
    relative_divergence_rms_gate: float,
    relative_curl_defect_rms_gate: float,
    relative_divergence_max_reference: float,
    relative_curl_defect_max_reference: float,
) -> dict[str, Any]:
    """Evaluate (a)-(d) on already-Cartesian arrays and return JSON metrics.

    Only :mod:`ns_certificate_lab.cartesian_validation` stencils are used.
    The caller owns the reconstruction, so fault-injection tests can corrupt
    the Cartesian arrays and re-enter exactly this code path.
    """

    state = grid.validate_vector(velocity, name="velocity")
    oracle = grid.validate_vector(vorticity_e18b, name="vorticity_e18b")
    if isinstance(interior_margin, bool) or not isinstance(interior_margin, int):
        raise ValueError("interior_margin must be a nonnegative integer")
    margin = int(interior_margin)
    if margin < 0:
        raise ValueError("interior_margin must be a nonnegative integer")

    x, y, z = grid.mesh()
    divergence = cartesian_divergence(grid, state)
    curl = cartesian_curl(grid, state)
    curl_defect = curl - oracle

    azimuthal_curl = -y * curl[0] + x * curl[1]
    azimuthal_oracle = -y * oracle[0] + x * oracle[1]
    azimuthal_defect = azimuthal_curl - azimuthal_oracle

    # grad[i][j] = d_j u_i, built only from the public Cartesian gradient.
    grad = np.stack(
        [cartesian_gradient(grid, state[component]) for component in range(3)],
        axis=0,
    )
    gradient_row_sum_max = float(np.max(np.abs(grad).sum(axis=1)))
    divergence_term_sum_max = float(
        np.max(np.abs(grad[0, 0]) + np.abs(grad[1, 1]) + np.abs(grad[2, 2]))
    )
    curl_term_sum_max = float(
        np.max(
            np.maximum(
                np.maximum(
                    np.abs(grad[2, 1]) + np.abs(grad[1, 2]),
                    np.abs(grad[0, 2]) + np.abs(grad[2, 0]),
                ),
                np.abs(grad[1, 0]) + np.abs(grad[0, 1]),
            )
        )
    )
    azimuthal_cancellation_max = float(
        np.max(
            np.abs(y) * (np.abs(curl[0]) + np.abs(oracle[0]))
            + np.abs(x) * (np.abs(curl[1]) + np.abs(oracle[1]))
        )
    )
    azimuthal_signal_max = float(np.max(np.abs(azimuthal_oracle)))
    vorticity_scale = float(np.max(np.abs(oracle)))
    velocity_scale = float(np.max(np.abs(state)))
    if gradient_row_sum_max <= _TINY:
        raise ValueError(
            "reconstructed velocity gradient vanishes; the relative audit "
            "has no denominator"
        )
    scales = {
        "velocity_max_abs": velocity_scale,
        "vorticity_e18b_max_abs": vorticity_scale,
        "velocity_gradient_row_sum_max": gradient_row_sum_max,
        "divergence_term_sum_max": divergence_term_sum_max,
        "curl_term_sum_max": curl_term_sum_max,
        "azimuthal_cancellation_max": azimuthal_cancellation_max,
        "azimuthal_signal_max": azimuthal_signal_max,
    }
    signal_degenerate = azimuthal_signal_max <= 1.0e-9 * max(
        azimuthal_cancellation_max,
        _TINY,
    )

    regions: dict[str, Any] = {}
    for label, active_margin in (("full", 0), ("interior", margin)):
        if active_margin and any(
            size <= 2 * active_margin
            for size, periodic in zip(grid.shape, grid.periodic)
            if not periodic
        ):
            raise ValueError("interior_margin leaves an empty norm region")
        window = tuple(
            slice(None)
            if periodic or active_margin == 0
            else slice(active_margin, -active_margin)
            for periodic in grid.periodic
        )
        vector_window = (slice(None), *window)
        offset = tuple(
            0 if periodic or active_margin == 0 else active_margin
            for periodic in grid.periodic
        )
        region_divergence = divergence[window]
        region_curl_defect = curl_defect[vector_window]
        region_azimuthal = azimuthal_defect[window]
        absolute = {
            "divergence_rms": _rms(region_divergence),
            "divergence_max": float(np.max(np.abs(region_divergence))),
            "curl_defect_rms": _rms(region_curl_defect),
            "curl_defect_max": float(np.max(np.abs(region_curl_defect))),
            "azimuthal_defect_rms": _rms(region_azimuthal),
            "azimuthal_defect_max": float(np.max(np.abs(region_azimuthal))),
        }
        relative = {
            "divergence_rms_over_gradient_scale": _ratio(
                absolute["divergence_rms"], gradient_row_sum_max
            ),
            "divergence_max_over_gradient_scale": _ratio(
                absolute["divergence_max"], gradient_row_sum_max
            ),
            "divergence_rms_over_divergence_term_scale": _ratio(
                absolute["divergence_rms"], divergence_term_sum_max
            ),
            "divergence_max_over_divergence_term_scale": _ratio(
                absolute["divergence_max"], divergence_term_sum_max
            ),
            "curl_defect_rms_over_gradient_scale": _ratio(
                absolute["curl_defect_rms"], gradient_row_sum_max
            ),
            "curl_defect_max_over_gradient_scale": _ratio(
                absolute["curl_defect_max"], gradient_row_sum_max
            ),
            "curl_defect_rms_over_curl_term_scale": _ratio(
                absolute["curl_defect_rms"], curl_term_sum_max
            ),
            "curl_defect_max_over_curl_term_scale": _ratio(
                absolute["curl_defect_max"], curl_term_sum_max
            ),
            "curl_defect_max_over_vorticity_scale": _ratio(
                absolute["curl_defect_max"], vorticity_scale
            ),
            "azimuthal_defect_rms_over_cancellation_scale": _ratio(
                absolute["azimuthal_defect_rms"], azimuthal_cancellation_max
            ),
            "azimuthal_defect_max_over_cancellation_scale": _ratio(
                absolute["azimuthal_defect_max"], azimuthal_cancellation_max
            ),
            "azimuthal_defect_max_over_signal_scale": (
                None
                if signal_degenerate
                else _ratio(absolute["azimuthal_defect_max"], azimuthal_signal_max)
            ),
        }
        regions[label] = {
            "margin": int(active_margin),
            "absolute": absolute,
            "relative": relative,
            "locations": {
                "divergence_max": _location(
                    region_divergence, (x, y, z), offset=offset, vector=False
                ),
                "curl_defect_max": _location(
                    region_curl_defect, (x, y, z), offset=offset, vector=True
                ),
                "azimuthal_defect_max": _location(
                    region_azimuthal, (x, y, z), offset=offset, vector=False
                ),
            },
        }

    # Re-enter the cartesian_validation gate itself.  Its tolerances are
    # absolute, so the configured RELATIVE thresholds are converted with the
    # measured gradient denominator; the converted values are recorded.  The
    # RMS tolerances come from the hard gates and the maximum tolerances from
    # the record-only references, matching how the verdict is assembled.
    divergence_rms_tolerance = (
        float(relative_divergence_rms_gate) * gradient_row_sum_max
    )
    divergence_max_tolerance = (
        float(relative_divergence_max_reference) * gradient_row_sum_max
    )
    curl_rms_tolerance = float(relative_curl_defect_rms_gate) * gradient_row_sum_max
    curl_max_tolerance = (
        float(relative_curl_defect_max_reference) * gradient_row_sum_max
    )
    report = audit_cartesian_reconstruction(
        grid,
        velocity=state,
        expected_vorticity=oracle,
        divergence_rms_tolerance=divergence_rms_tolerance,
        divergence_max_tolerance=divergence_max_tolerance,
        curl_rms_tolerance=curl_rms_tolerance,
        curl_max_tolerance=curl_max_tolerance,
        interior_margin=margin,
    )

    return {
        "shape": list(grid.shape),
        "spacings": list(grid.spacings),
        "periodic": [bool(value) for value in grid.periodic],
        "interior_margin": margin,
        "scales": scales,
        "azimuthal_signal_degenerate": bool(signal_degenerate),
        "regions": regions,
        "cartesian_validation_gate": {
            "divergence_rms_tolerance": divergence_rms_tolerance,
            "divergence_max_tolerance": divergence_max_tolerance,
            "curl_rms_tolerance": curl_rms_tolerance,
            "curl_max_tolerance": curl_max_tolerance,
            "divergence_rms": report.divergence_rms,
            "divergence_max": report.divergence_max,
            "curl_defect_rms": report.curl_defect_rms,
            "curl_defect_max": report.curl_defect_max,
            "divergence_rms_passed": bool(report.divergence_rms_passed),
            "divergence_max_passed": bool(report.divergence_max_passed),
            "curl_rms_passed": bool(report.curl_rms_passed),
            "curl_max_passed": bool(report.curl_max_passed),
            "all_points_finite": bool(report.all_points_finite),
        },
    }


def reconstruct_checkpoint(
    candidate: CandidateArtifact,
    grid: UniformCartesianGrid,
) -> CandidateCartesianArrays:
    """Reconstruct one loaded candidate with the existing independent adapter."""

    return reconstruct_loaded_candidate_on_cartesian(candidate, grid)


def _observed_order(
    coarse_value: float | None,
    fine_value: float | None,
    spacing_ratio: float,
) -> float | None:
    if (
        coarse_value is None
        or fine_value is None
        or coarse_value <= 0.0
        or fine_value <= 0.0
        or spacing_ratio <= 1.0
    ):
        return None
    value = math.log(coarse_value / fine_value) / math.log(spacing_ratio)
    return float(value) if math.isfinite(value) else None


def audit_checkpoint(
    path: Path,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Load one checkpoint, reconstruct it and evaluate every audit metric."""

    grid_spec = config["cartesian_grid"]
    gates = config["hard_gates"]
    thresholds = config["record_only_thresholds"]
    candidate = load_candidate(path)
    check_box_is_inscribed(candidate, grid_spec)
    if candidate.description is None:
        raise ValueError(f"checkpoint lacks a v2 description: {path}")

    evaluations: dict[str, dict[str, Any]] = {}
    for label, shape_key in (("fine", "shape"), ("coarse", "coarse_shape")):
        shape = _shape(grid_spec[shape_key], name=shape_key)
        grid = build_cartesian_grid(grid_spec, shape)
        arrays = reconstruct_checkpoint(candidate, grid)
        evaluations[label] = audit_cartesian_arrays(
            grid,
            velocity=arrays.velocity,
            vorticity_e18b=arrays.vorticity_e18b,
            interior_margin=int(grid_spec["interior_margin"]),
            relative_divergence_rms_gate=float(
                gates["maximum_relative_divergence_rms"]
            ),
            relative_curl_defect_rms_gate=float(
                gates["maximum_relative_curl_defect_rms"]
            ),
            relative_divergence_max_reference=float(
                thresholds["reference_relative_divergence_max"]
            ),
            relative_curl_defect_max_reference=float(
                thresholds["reference_relative_curl_defect_max"]
            ),
        )

    fine = evaluations["fine"]["regions"]["interior"]["relative"]
    coarse = evaluations["coarse"]["regions"]["interior"]["relative"]
    spacing_ratio = max(evaluations["coarse"]["spacings"]) / max(
        evaluations["fine"]["spacings"]
    )
    refinement = {
        "spacing_ratio": float(spacing_ratio),
        "relative_divergence_rms_coarse": coarse[
            "divergence_rms_over_gradient_scale"
        ],
        "relative_divergence_rms_fine": fine["divergence_rms_over_gradient_scale"],
        "relative_divergence_rms_observed_order": _observed_order(
            coarse["divergence_rms_over_gradient_scale"],
            fine["divergence_rms_over_gradient_scale"],
            spacing_ratio,
        ),
        "relative_curl_defect_rms_coarse": coarse[
            "curl_defect_rms_over_gradient_scale"
        ],
        "relative_curl_defect_rms_fine": fine["curl_defect_rms_over_gradient_scale"],
        "relative_curl_defect_rms_observed_order": _observed_order(
            coarse["curl_defect_rms_over_gradient_scale"],
            fine["curl_defect_rms_over_gradient_scale"],
            spacing_ratio,
        ),
    }
    refinement["divergence_decreases_under_refinement"] = bool(
        refinement["relative_divergence_rms_fine"] is not None
        and refinement["relative_divergence_rms_coarse"] is not None
        and refinement["relative_divergence_rms_fine"]
        < refinement["relative_divergence_rms_coarse"]
    )
    refinement["curl_defect_decreases_under_refinement"] = bool(
        refinement["relative_curl_defect_rms_fine"] is not None
        and refinement["relative_curl_defect_rms_coarse"] is not None
        and refinement["relative_curl_defect_rms_fine"]
        < refinement["relative_curl_defect_rms_coarse"]
    )

    manifest_path = path.with_suffix(".manifest.json")
    return {
        "checkpoint": path.relative_to(REPOSITORY_ROOT).as_posix(),
        "archive_sha256": sha256_file(path),
        "manifest_sha256": (
            sha256_file(manifest_path) if manifest_path.is_file() else None
        ),
        "source": {
            "nr": int(candidate.grid.nr),
            "nz": int(candidate.grid.nz),
            "dr": float(candidate.grid.dr),
            "dz": float(candidate.grid.dz),
            "r_max": float(candidate.grid.r[-1]),
            "periodic_z": bool(candidate.grid.periodic_z),
            "physical_time": float(candidate.description.physical_time),
            "viscosity": float(candidate.description.viscosity),
            "field_names": sorted(candidate.fields),
            "experiment_id": candidate.metadata.get("experiment_id"),
        },
        "evaluations": evaluations,
        "audit_grid_refinement": refinement,
    }


# --------------------------------------------------------------------------
# the deliberately unfinished primitive residual
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class SnapshotPairResidual:
    """Midpoint primitive residual assembled from two Cartesian snapshots."""

    residual: PrimitiveNSResidual
    midpoint_time: float
    time_step: float
    velocity_midpoint: npt.NDArray[np.float64]
    velocity_time_derivative: npt.NDArray[np.float64]


def primitive_residual_from_snapshot_pair(
    grid: UniformCartesianGrid,
    *,
    velocity_early: npt.ArrayLike,
    velocity_late: npt.ArrayLike,
    time_early: float,
    time_late: float,
    midpoint_pressure: npt.ArrayLike,
    viscosity: float,
    body_force: npt.ArrayLike | None = None,
) -> SnapshotPairResidual:
    r"""Assemble the primitive residual at the midpoint of a checkpoint pair.

    Given Cartesian velocities at ``time_early`` and ``time_late`` and a
    kinematic pressure at the temporal midpoint, this forms

    ``u_t ~ (u_late - u_early) / dt`` and ``u ~ (u_late + u_early) / 2``,

    both of which are second-order accurate at ``t = (t_early+t_late)/2``, and
    hands them to :func:`cartesian_validation.primitive_ns_residual`.  No new
    momentum-term stencil is introduced here.

    **This function is not applicable to the shipped Hou checkpoints.**  They
    store no pressure at all, and their 5e-4 spacing is roughly 500 integrator
    steps, so neither input above is available at the accuracy the formula
    assumes.  See :data:`PRIMITIVE_RESIDUAL_GAP`.
    """

    early = grid.validate_vector(velocity_early, name="velocity_early")
    late = grid.validate_vector(velocity_late, name="velocity_late")
    start = _finite_number(time_early, name="time_early")
    stop = _finite_number(time_late, name="time_late")
    step = stop - start
    if step <= 0.0:
        raise ValueError("time_late must be strictly greater than time_early")
    velocity_t = (late - early) / step
    velocity_mid = 0.5 * (late + early)
    residual = primitive_ns_residual(
        grid,
        velocity=velocity_mid,
        velocity_t=velocity_t,
        pressure=midpoint_pressure,
        viscosity=viscosity,
        body_force=body_force,
    )
    return SnapshotPairResidual(
        residual=residual,
        midpoint_time=0.5 * (start + stop),
        time_step=float(step),
        velocity_midpoint=velocity_mid,
        velocity_time_derivative=velocity_t,
    )


# --------------------------------------------------------------------------
# evaluation and evidence
# --------------------------------------------------------------------------


def _interior_relative(record: Mapping[str, Any], name: str) -> float | None:
    return record["evaluations"]["fine"]["regions"]["interior"]["relative"][name]


def _finite_values(value: Any) -> Iterable[float]:
    if isinstance(value, bool) or value is None:
        return ()
    if isinstance(value, (int, float)):
        return (float(value),)
    if isinstance(value, dict):
        return [item for entry in value.values() for item in _finite_values(entry)]
    if isinstance(value, list):
        return [item for entry in value for item in _finite_values(entry)]
    return ()


def resolution_comparison(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group the fine-grid relative metrics by snapshot time, ordered by nr."""

    times = sorted({record["source"]["physical_time"] for record in records})
    table: list[dict[str, Any]] = []
    for moment in times:
        rows = sorted(
            (
                record
                for record in records
                if record["source"]["physical_time"] == moment
            ),
            key=lambda record: record["source"]["nr"],
        )
        table.append(
            {
                "physical_time": float(moment),
                "rows": [
                    {
                        "nr": row["source"]["nr"],
                        "nz": row["source"]["nz"],
                        "checkpoint": row["checkpoint"],
                        "relative_divergence_rms": _interior_relative(
                            row, "divergence_rms_over_gradient_scale"
                        ),
                        "relative_divergence_max": _interior_relative(
                            row, "divergence_max_over_gradient_scale"
                        ),
                        "relative_curl_defect_rms": _interior_relative(
                            row, "curl_defect_rms_over_gradient_scale"
                        ),
                        "relative_curl_defect_max": _interior_relative(
                            row, "curl_defect_max_over_gradient_scale"
                        ),
                        "relative_azimuthal_defect_max": _interior_relative(
                            row, "azimuthal_defect_max_over_cancellation_scale"
                        ),
                        "audit_grid_divergence_order": row["audit_grid_refinement"][
                            "relative_divergence_rms_observed_order"
                        ],
                        "audit_grid_curl_defect_order": row["audit_grid_refinement"][
                            "relative_curl_defect_rms_observed_order"
                        ],
                    }
                    for row in rows
                ],
            }
        )
    return table


def evaluate(
    config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, bool], dict[str, bool]]:
    """Audit every configured checkpoint without writing evidence."""

    validate_config(config)
    records = [
        audit_checkpoint(_checkpoint_path(entry), config)
        for entry in config["checkpoint_paths"]
    ]
    gates = config["hard_gates"]
    thresholds = config["record_only_thresholds"]

    divergence_rms = [
        _interior_relative(record, "divergence_rms_over_gradient_scale")
        for record in records
    ]
    curl_rms = [
        _interior_relative(record, "curl_defect_rms_over_gradient_scale")
        for record in records
    ]
    divergence_max = [
        _interior_relative(record, "divergence_max_over_gradient_scale")
        for record in records
    ]
    curl_max = [
        _interior_relative(record, "curl_defect_max_over_gradient_scale")
        for record in records
    ]
    azimuthal_max = [
        _interior_relative(record, "azimuthal_defect_max_over_cancellation_scale")
        for record in records
    ]
    checks = {
        "all_checkpoints_audited": len(records) == len(config["checkpoint_paths"]),
        "all_metrics_finite": all(
            math.isfinite(value) for value in _finite_values(records)
        ),
        "all_points_finite": all(
            evaluation["cartesian_validation_gate"]["all_points_finite"]
            for record in records
            for evaluation in record["evaluations"].values()
        ),
        "relative_divergence_rms_within_gate": all(
            value is not None
            and value <= float(gates["maximum_relative_divergence_rms"])
            for value in divergence_rms
        ),
        "relative_curl_defect_rms_within_gate": all(
            value is not None
            and value <= float(gates["maximum_relative_curl_defect_rms"])
            for value in curl_rms
        ),
        "cartesian_validation_rms_gates_passed": all(
            evaluation["cartesian_validation_gate"]["divergence_rms_passed"]
            and evaluation["cartesian_validation_gate"]["curl_rms_passed"]
            for record in records
            for evaluation in record["evaluations"].values()
        ),
        "audit_defects_decrease_under_refinement": all(
            record["audit_grid_refinement"]["divergence_decreases_under_refinement"]
            and record["audit_grid_refinement"][
                "curl_defect_decreases_under_refinement"
            ]
            for record in records
        ),
    }
    record_only = {
        "relative_divergence_max_below_reference": all(
            value is not None
            and value <= float(thresholds["reference_relative_divergence_max"])
            for value in divergence_max
        ),
        "relative_curl_defect_max_below_reference": all(
            value is not None
            and value <= float(thresholds["reference_relative_curl_defect_max"])
            for value in curl_max
        ),
        "relative_azimuthal_defect_max_below_reference": all(
            value is not None
            and value
            <= float(thresholds["reference_relative_azimuthal_defect_max"])
            for value in azimuthal_max
        ),
    }
    return records, checks, record_only


def _metric_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        for label, evaluation in sorted(record["evaluations"].items()):
            for region, block in sorted(evaluation["regions"].items()):
                rows.append(
                    {
                        "checkpoint": record["checkpoint"],
                        "nr": record["source"]["nr"],
                        "nz": record["source"]["nz"],
                        "physical_time": record["source"]["physical_time"],
                        "audit_grid": label,
                        "audit_shape": "x".join(
                            str(value) for value in evaluation["shape"]
                        ),
                        "region": region,
                        **{
                            key: block["absolute"][key]
                            for key in sorted(block["absolute"])
                        },
                        **{
                            key: (
                                "" if block["relative"][key] is None
                                else block["relative"][key]
                            )
                            for key in sorted(block["relative"])
                        },
                        **{
                            f"scale_{key}": evaluation["scales"][key]
                            for key in sorted(evaluation["scales"])
                        },
                    }
                )
    return rows


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


def run(config: Mapping[str, Any], output_dir: Path) -> dict[str, Any]:
    """Audit every checkpoint and write a checksummed evidence bundle."""

    validate_config(config)
    if not _inside_repository(output_dir):
        raise ValueError("output directory must remain inside this repository")
    resolved_output = output_dir.resolve()
    outputs_root = (REPOSITORY_ROOT / "outputs" / "hou_early_time_v1").resolve()
    if resolved_output == outputs_root or outputs_root in resolved_output.parents:
        raise ValueError(
            "refusing to write inside the audited source evidence directory"
        )
    provenance = collect_runtime_provenance()
    config_bytes = canonical_json_bytes(dict(config))
    _prepare_output(output_dir)

    records, checks, record_only = evaluate(config)
    summary: dict[str, Any] = {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "experiment_id": config["experiment_id"],
        "source_experiment": config["source_experiment"],
        "interpretation": config["interpretation"],
        "method": {
            "reconstruction": (
                "ns_certificate_lab.cartesian_candidate_adapter: private "
                "second-order (r,z) stencils and bilinear sampling of the "
                "reloaded nodal arrays, producing E-18a velocity and the "
                "E-18b full-vorticity oracle"
            ),
            "checker": (
                "ns_certificate_lab.cartesian_validation: second-order "
                "Cartesian stencils owned by that module alone; no import of "
                "operators.py, pde.py or any cylindrical result"
            ),
            "divergence": "E-02 physical 3D divergence of the reconstructed velocity",
            "curl": "E-18b full curl defect curl(u) - omega_E18b",
            "cylindrical_agreement": (
                "-y (curl u)_x + x (curl u)_y compared with r^2 omega1 "
                "recovered as -y omega_x + x omega_y from the same E-18b "
                "oracle; this is omega^theta = r omega1 without dividing by r"
            ),
            "audit_grid_refinement": (
                "every checkpoint is evaluated on both the fine and the "
                "coarse audit shape so that audit truncation can be "
                "distinguished from a snapshot defect"
            ),
            "norm_region": (
                "interior_margin trims the non-periodic x and y axes only; "
                "the periodic z seam is always retained"
            ),
        },
        "denominators": dict(DENOMINATORS),
        "cartesian_grid": dict(config["cartesian_grid"]),
        "hard_gates": dict(config["hard_gates"]),
        "record_only_thresholds": dict(config["record_only_thresholds"]),
        "checkpoints": records,
        "resolution_comparison": resolution_comparison(records),
        "acceptance_checks": checks,
        "record_only_checks": record_only,
        "accepted_as_independent_cartesian_snapshot_audit": bool(all(checks.values())),
        "primitive_navier_stokes_residual_gap": dict(PRIMITIVE_RESIDUAL_GAP),
        "limitations": list(LIMITATIONS),
        "reproducibility": {
            "config_sha256": sha256_bytes(config_bytes),
            "runtime_provenance": provenance,
            "command": (
                "python experiments/run_hou_snapshot_cartesian_audit.py "
                "--config configs/hou_snapshot_cartesian_audit.json "
                "--output-dir outputs/hou_snapshot_cartesian_audit_v1"
            ),
        },
    }

    write_with_digest(output_dir / "config.snapshot.json", config_bytes)
    write_with_digest(output_dir / "summary.json", canonical_json_bytes(summary))
    rows = _metric_rows(records)
    with (output_dir / "metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    artifact_names = [
        "config.snapshot.json",
        "config.snapshot.json.sha256",
        "metrics.csv",
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
            record["checkpoint"]: {
                "archive_sha256": record["archive_sha256"],
                "manifest_sha256": record["manifest_sha256"],
            }
            for record in records
        },
    }
    write_with_digest(output_dir / "manifest.json", canonical_json_bytes(manifest))
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "hou_snapshot_cartesian_audit.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "hou_snapshot_cartesian_audit_v1",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    try:
        config = strict_json_loads(
            config_path.read_text(encoding="utf-8"),
            label="hou snapshot cartesian audit config",
        )
        if not isinstance(config, dict):
            raise ValueError("config must be a JSON object")
        summary = run(config, output_dir)
    except (OSError, ValueError, FileExistsError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary["acceptance_checks"], indent=2, sort_keys=True))
    return 0 if summary["accepted_as_independent_cartesian_snapshot_audit"] else 2


if __name__ == "__main__":
    sys.exit(main())
