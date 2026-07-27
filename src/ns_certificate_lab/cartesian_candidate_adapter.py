"""Adapter from saved axisymmetric nodal candidates to Cartesian audit arrays.

This is deliberately an adapter, not part of the Cartesian checker.  It may
interpret the repository's candidate artifact format, but it does not import
or call ``operators.py`` and does not consume a precomputed cylindrical
velocity or derivative.  Its private second-order cylindrical stencils and
bilinear interpolator are owned by this module.

The returned vorticity is the complete E-18b oracle reconstructed from
``(u1, omega1)``:

``omega_x=-x*u1_z-y*omega1``,
``omega_y=-y*u1_z+x*omega1``,
``omega_z=2*u1+r*u1_r``.

It is intended for cross-checking candidate arrays after a verified save/load
round trip.  It does not certify interpolation or truncation errors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from .artifacts import CandidateArtifact
from .cartesian_validation import UniformCartesianGrid, VectorArray

ScalarArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class CandidateCartesianArrays:
    """Cartesian velocity and full-vorticity oracle from one loaded candidate."""

    velocity: VectorArray
    vorticity_e18b: VectorArray


def _radial_first(values: ScalarArray, spacing: float) -> ScalarArray:
    """Independent second-order radial derivative, including both boundaries."""

    if values.ndim != 2 or values.shape[0] < 3:
        raise ValueError("candidate radial derivative requires a two-dimensional array")
    result = np.empty_like(values)
    result[1:-1] = (values[2:] - values[:-2]) / (2.0 * spacing)
    result[0] = (-3.0 * values[0] + 4.0 * values[1] - values[2]) / (
        2.0 * spacing
    )
    result[-1] = (
        3.0 * values[-1] - 4.0 * values[-2] + values[-3]
    ) / (2.0 * spacing)
    return result


def _axial_first(candidate: CandidateArtifact, values: ScalarArray) -> ScalarArray:
    """Independent second-order axial derivative."""

    spacing = candidate.grid.dz
    if candidate.grid.periodic_z:
        return (
            np.roll(values, -1, axis=1) - np.roll(values, 1, axis=1)
        ) / (2.0 * spacing)
    if values.shape[1] < 3:
        raise ValueError("candidate axial derivative requires at least three points")
    result = np.empty_like(values)
    result[:, 1:-1] = (values[:, 2:] - values[:, :-2]) / (2.0 * spacing)
    result[:, 0] = (
        -3.0 * values[:, 0] + 4.0 * values[:, 1] - values[:, 2]
    ) / (2.0 * spacing)
    result[:, -1] = (
        3.0 * values[:, -1] - 4.0 * values[:, -2] + values[:, -3]
    ) / (2.0 * spacing)
    return result


def _bilinear_sample(
    candidate: CandidateArtifact,
    values: ScalarArray,
    radius: ScalarArray,
    axial: ScalarArray,
) -> ScalarArray:
    """Sample a cylindrical nodal array at Cartesian ``(r,z)`` points."""

    source_grid = candidate.grid
    if values.shape != source_grid.shape or not np.all(np.isfinite(values)):
        raise ValueError("candidate interpolation source is invalid")
    radial_tolerance = 32.0 * np.finfo(np.float64).eps * max(
        1.0,
        float(source_grid.r[-1]),
    )
    if (
        np.any(radius < -radial_tolerance)
        or np.any(radius > source_grid.r[-1] + radial_tolerance)
    ):
        raise ValueError("Cartesian audit grid extends beyond candidate r_max")
    radial_coordinate = np.clip(radius, 0.0, source_grid.r[-1]) / source_grid.dr
    radial_lower = np.floor(radial_coordinate).astype(np.intp)
    radial_lower = np.clip(radial_lower, 0, source_grid.nr - 2)
    radial_fraction = radial_coordinate - radial_lower

    if source_grid.periodic_z:
        assert source_grid.z_period is not None
        wrapped = (
            np.mod(axial - source_grid.z[0], source_grid.z_period)
            + source_grid.z[0]
        )
        axial_coordinate = (wrapped - source_grid.z[0]) / source_grid.dz
        axial_lower = np.floor(axial_coordinate).astype(np.intp) % source_grid.nz
        axial_upper = (axial_lower + 1) % source_grid.nz
        axial_fraction = axial_coordinate - np.floor(axial_coordinate)
    else:
        axial_tolerance = 32.0 * np.finfo(np.float64).eps * max(
            1.0,
            abs(float(source_grid.z[0])),
            abs(float(source_grid.z[-1])),
        )
        if (
            np.any(axial < source_grid.z[0] - axial_tolerance)
            or np.any(axial > source_grid.z[-1] + axial_tolerance)
        ):
            raise ValueError("Cartesian audit grid extends beyond candidate z range")
        axial_coordinate = (
            np.clip(axial, source_grid.z[0], source_grid.z[-1])
            - source_grid.z[0]
        ) / source_grid.dz
        axial_lower = np.floor(axial_coordinate).astype(np.intp)
        axial_lower = np.clip(axial_lower, 0, source_grid.nz - 2)
        axial_upper = axial_lower + 1
        axial_fraction = axial_coordinate - axial_lower

    radial_upper = radial_lower + 1
    value_00 = values[radial_lower, axial_lower]
    value_10 = values[radial_upper, axial_lower]
    value_01 = values[radial_lower, axial_upper]
    value_11 = values[radial_upper, axial_upper]
    lower_axial = value_00 + radial_fraction * (value_10 - value_00)
    upper_axial = value_01 + radial_fraction * (value_11 - value_01)
    return lower_axial + axial_fraction * (upper_axial - lower_axial)


def reconstruct_loaded_candidate_on_cartesian(
    candidate: CandidateArtifact,
    grid: UniformCartesianGrid,
) -> CandidateCartesianArrays:
    """Reconstruct velocity and the E-18b vorticity oracle from loaded arrays.

    Required fields are ``u1``, ``psi1``, and ``omega1``.  No derivative or
    velocity stored elsewhere in the artifact is trusted.
    """

    missing = {"u1", "psi1", "omega1"} - set(candidate.fields)
    if missing:
        raise ValueError(
            "candidate is missing Cartesian audit fields: "
            + ", ".join(sorted(missing))
        )
    u1 = candidate.grid.validate_field(candidate.fields["u1"], name="u1")
    psi1 = candidate.grid.validate_field(candidate.fields["psi1"], name="psi1")
    omega1 = candidate.grid.validate_field(candidate.fields["omega1"], name="omega1")

    # These stencils are local to this adapter and do not call the production
    # cylindrical differentiation path.
    psi1_r = _radial_first(psi1, candidate.grid.dr)
    psi1_z = _axial_first(candidate, psi1)
    u1_r = _radial_first(u1, candidate.grid.dr)
    u1_z = _axial_first(candidate, u1)

    x, y, z = grid.mesh()
    radius = np.hypot(x, y)
    sampled_u1 = _bilinear_sample(candidate, u1, radius, z)
    sampled_psi1 = _bilinear_sample(candidate, psi1, radius, z)
    sampled_psi1_r = _bilinear_sample(candidate, psi1_r, radius, z)
    sampled_psi1_z = _bilinear_sample(candidate, psi1_z, radius, z)
    sampled_u1_r = _bilinear_sample(candidate, u1_r, radius, z)
    sampled_u1_z = _bilinear_sample(candidate, u1_z, radius, z)
    sampled_omega1 = _bilinear_sample(candidate, omega1, radius, z)

    velocity = np.stack(
        (
            -x * sampled_psi1_z - y * sampled_u1,
            -y * sampled_psi1_z + x * sampled_u1,
            2.0 * sampled_psi1 + radius * sampled_psi1_r,
        ),
        axis=0,
    )
    vorticity_e18b = np.stack(
        (
            -x * sampled_u1_z - y * sampled_omega1,
            -y * sampled_u1_z + x * sampled_omega1,
            2.0 * sampled_u1 + radius * sampled_u1_r,
        ),
        axis=0,
    )
    if not np.all(np.isfinite(velocity)) or not np.all(np.isfinite(vorticity_e18b)):
        raise ValueError("candidate Cartesian reconstruction produced non-finite values")
    return CandidateCartesianArrays(
        velocity=velocity,
        vorticity_e18b=vorticity_e18b,
    )
