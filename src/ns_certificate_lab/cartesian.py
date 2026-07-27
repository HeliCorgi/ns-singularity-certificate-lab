"""Independent physical-space reconstruction checks.

This module deliberately keeps its radial and axial first-difference
implementation separate from :mod:`ns_certificate_lab.operators`.  It is an
audit path for reconstructing the physical three-dimensional velocity and for
checking the azimuthal curl identity

    omega^theta = partial_z u^r - partial_r u^z.

It is not a second time integrator and it does not turn a discrete candidate
into a proof.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from .artifacts import CandidateArtifact
from .grid import AxisymmetricGrid, FloatArray


@dataclass(frozen=True)
class CartesianVelocity:
    """Velocity and coordinates on a cylindrical tensor-product sampling."""

    x: npt.NDArray[np.float64]
    y: npt.NDArray[np.float64]
    z: npt.NDArray[np.float64]
    u_x: npt.NDArray[np.float64]
    u_y: npt.NDArray[np.float64]
    u_z: npt.NDArray[np.float64]
    theta: npt.NDArray[np.float64]


def _first_r_independent(grid: AxisymmetricGrid, field: FloatArray) -> FloatArray:
    """Fourth-order radial derivative, independent of the production operator."""

    if grid.nr < 5:
        raise ValueError("independent fourth-order audit requires nr >= 5")
    out = np.empty_like(field)
    h = grid.dr
    out[2:-2] = (
        field[:-4] - 8.0 * field[1:-3] + 8.0 * field[3:-1] - field[4:]
    ) / (12.0 * h)
    out[0] = (
        -25.0 * field[0]
        + 48.0 * field[1]
        - 36.0 * field[2]
        + 16.0 * field[3]
        - 3.0 * field[4]
    ) / (12.0 * h)
    out[1] = (
        -3.0 * field[0]
        - 10.0 * field[1]
        + 18.0 * field[2]
        - 6.0 * field[3]
        + field[4]
    ) / (12.0 * h)
    out[-2] = (
        -field[-5]
        + 6.0 * field[-4]
        - 18.0 * field[-3]
        + 10.0 * field[-2]
        + 3.0 * field[-1]
    ) / (12.0 * h)
    out[-1] = (
        25.0 * field[-1]
        - 48.0 * field[-2]
        + 36.0 * field[-3]
        - 16.0 * field[-4]
        + 3.0 * field[-5]
    ) / (12.0 * h)
    return out


def _first_z_independent(grid: AxisymmetricGrid, field: FloatArray) -> FloatArray:
    """Fourth-order axial derivative, independent of the production operator."""

    h = grid.dz
    if grid.periodic_z:
        return (
            np.roll(field, 2, axis=1)
            - 8.0 * np.roll(field, 1, axis=1)
            + 8.0 * np.roll(field, -1, axis=1)
            - np.roll(field, -2, axis=1)
        ) / (12.0 * h)
    out = np.empty_like(field)
    out[:, 2:-2] = (
        field[:, :-4]
        - 8.0 * field[:, 1:-3]
        + 8.0 * field[:, 3:-1]
        - field[:, 4:]
    ) / (12.0 * h)
    out[:, 0] = (
        -25.0 * field[:, 0]
        + 48.0 * field[:, 1]
        - 36.0 * field[:, 2]
        + 16.0 * field[:, 3]
        - 3.0 * field[:, 4]
    ) / (12.0 * h)
    out[:, 1] = (
        -3.0 * field[:, 0]
        - 10.0 * field[:, 1]
        + 18.0 * field[:, 2]
        - 6.0 * field[:, 3]
        + field[:, 4]
    ) / (12.0 * h)
    out[:, -2] = (
        -field[:, -5]
        + 6.0 * field[:, -4]
        - 18.0 * field[:, -3]
        + 10.0 * field[:, -2]
        + 3.0 * field[:, -1]
    ) / (12.0 * h)
    out[:, -1] = (
        25.0 * field[:, -1]
        - 48.0 * field[:, -2]
        + 36.0 * field[:, -3]
        - 16.0 * field[:, -4]
        + 3.0 * field[:, -5]
    ) / (12.0 * h)
    return out


def recover_meridional_velocity_independent(
    grid: AxisymmetricGrid,
    psi1: npt.ArrayLike,
) -> tuple[FloatArray, FloatArray]:
    """Recover ``(u^r,u^z)`` without calling the production differentiators."""

    stream = grid.validate_field(psi1, name="psi1")
    radius = grid.r[:, None]
    stream_r = _first_r_independent(grid, stream)
    stream_z = _first_z_independent(grid, stream)
    u_r = -radius * stream_z
    u_z = 2.0 * stream + radius * stream_r
    u_r[0] = 0.0
    return u_r, u_z


def omega_theta_from_velocity_independent(
    grid: AxisymmetricGrid,
    u_r: npt.ArrayLike,
    u_z: npt.ArrayLike,
) -> FloatArray:
    """Compute physical ``omega^theta`` from velocity via an audit stencil."""

    radial = grid.validate_field(u_r, name="u_r")
    axial = grid.validate_field(u_z, name="u_z")
    omega_theta = (
        _first_z_independent(grid, radial)
        - _first_r_independent(grid, axial)
    )
    # A smooth axisymmetric vector field has odd omega^theta.  Record its
    # represented-axis value exactly; nearby points still test the stencil.
    omega_theta[0] = 0.0
    return omega_theta


def curl_consistency_defect(
    grid: AxisymmetricGrid,
    *,
    u_r: npt.ArrayLike,
    u_z: npt.ArrayLike,
    omega1: npt.ArrayLike,
) -> FloatArray:
    """Return ``curl(u)^theta - r omega1`` using the independent curl path."""

    transformed = grid.validate_field(omega1, name="omega1")
    physical = omega_theta_from_velocity_independent(grid, u_r, u_z)
    return physical - grid.r[:, None] * transformed


def omega_theta_from_cartesian_independent(
    grid: AxisymmetricGrid,
    velocity: CartesianVelocity,
) -> npt.NDArray[np.float64]:
    """Recover ``omega^theta`` after Cartesian reconstruction.

    The Cartesian components are first projected back onto the physical radial
    direction at each sampled angle.  Separate fourth-order stencils then
    evaluate ``partial_z u^r - partial_r u^z``.  This path is intentionally
    redundant: it catches component-map and sign errors that an elliptic
    ``-L5 psi1`` check alone cannot see.
    """

    expected = (grid.nr, grid.nz, velocity.theta.size)
    if velocity.theta.ndim != 1 or velocity.theta.size < 1:
        raise ValueError("Cartesian velocity has invalid theta coordinates")
    for name, values in (
        ("u_x", velocity.u_x),
        ("u_y", velocity.u_y),
        ("u_z", velocity.u_z),
    ):
        if values.shape != expected or not np.all(np.isfinite(values)):
            raise ValueError(f"Cartesian velocity {name} must have shape {expected}")
    output = np.empty(expected, dtype=np.float64)
    for index, angle in enumerate(velocity.theta):
        radial = (
            velocity.u_x[:, :, index] * np.cos(angle)
            + velocity.u_y[:, :, index] * np.sin(angle)
        )
        output[:, :, index] = (
            _first_z_independent(grid, radial)
            - _first_r_independent(grid, velocity.u_z[:, :, index])
        )
        output[0, :, index] = 0.0
    return output


def cartesian_curl_consistency_defect(
    grid: AxisymmetricGrid,
    *,
    velocity: CartesianVelocity,
    omega1: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """Return the 3D reconstruction/curl defect against ``r omega1``."""

    transformed = grid.validate_field(omega1, name="omega1")
    omega_theta = omega_theta_from_cartesian_independent(grid, velocity)
    target = np.broadcast_to(
        (grid.r[:, None] * transformed)[:, :, None],
        omega_theta.shape,
    )
    return omega_theta - target


def reconstruct_cartesian_velocity(
    grid: AxisymmetricGrid,
    *,
    u1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    theta: npt.ArrayLike,
) -> CartesianVelocity:
    """Reconstruct the physical velocity on selected azimuthal angles.

    Returned arrays have shape ``(n_r, n_z, n_theta)``.  The map is

    ``x=r cos(theta)``, ``y=r sin(theta)``,
    ``u_x=u^r cos(theta)-u^theta sin(theta)``,
    ``u_y=u^r sin(theta)+u^theta cos(theta)``.
    """

    swirl = grid.validate_field(u1, name="u1")
    angles = np.asarray(theta, dtype=np.float64)
    if angles.ndim != 1 or angles.size < 1 or not np.all(np.isfinite(angles)):
        raise ValueError("theta must be a nonempty finite one-dimensional array")

    u_r, u_z_cylindrical = recover_meridional_velocity_independent(grid, psi1)
    radius = grid.r[:, None]
    u_theta = radius * swirl
    cosine = np.cos(angles)[None, None, :]
    sine = np.sin(angles)[None, None, :]
    radius_3d = radius[:, :, None]
    z_3d = np.broadcast_to(
        grid.z[None, :, None],
        (grid.nr, grid.nz, angles.size),
    ).copy()
    return CartesianVelocity(
        x=np.broadcast_to(radius_3d * cosine, z_3d.shape).copy(),
        y=np.broadcast_to(radius_3d * sine, z_3d.shape).copy(),
        z=z_3d,
        u_x=u_r[:, :, None] * cosine - u_theta[:, :, None] * sine,
        u_y=u_r[:, :, None] * sine + u_theta[:, :, None] * cosine,
        u_z=np.broadcast_to(u_z_cylindrical[:, :, None], z_3d.shape).copy(),
        theta=angles.copy(),
    )


def reconstruct_candidate_cartesian_velocity(
    candidate: CandidateArtifact,
    *,
    theta: npt.ArrayLike,
) -> CartesianVelocity:
    """Reconstruct a loaded candidate, rejecting missing required fields."""

    missing = {"u1", "psi1"} - set(candidate.fields)
    if missing:
        names = ", ".join(sorted(missing))
        raise ValueError(f"candidate is missing Cartesian reconstruction fields: {names}")
    return reconstruct_cartesian_velocity(
        candidate.grid,
        u1=candidate.fields["u1"],
        psi1=candidate.fields["psi1"],
        theta=theta,
    )
