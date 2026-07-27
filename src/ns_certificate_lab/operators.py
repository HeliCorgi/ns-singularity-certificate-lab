"""Second-order finite differences for axisymmetric transformed variables.

The symbol ``laplacian_5d_formal`` names the scalar differential expression

    d_rr f + (3/r) d_r f + d_zz f.

It is a useful algebraic shorthand for the transformed axisymmetric equations;
it is *not* used as a five-dimensional incompressibility condition.  Physical
incompressibility is always checked with the three-dimensional cylindrical
divergence implemented by :func:`divergence_physical`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid, FloatArray


def derivative_r(
    grid: AxisymmetricGrid,
    field: npt.ArrayLike,
    *,
    even_at_axis: bool = False,
) -> FloatArray:
    """Second-order radial first derivative.

    When ``even_at_axis`` is true the exact regularity condition ``f_r(0)=0``
    is imposed at the first grid line.  No outer-boundary condition is imposed;
    a second-order backward stencil is used there.
    """

    f = grid.validate_field(field)
    out = np.empty_like(f)
    h = grid.dr
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * h)
    out[-1] = (3.0 * f[-1] - 4.0 * f[-2] + f[-3]) / (2.0 * h)
    if even_at_axis:
        out[0] = 0.0
    else:
        out[0] = (-3.0 * f[0] + 4.0 * f[1] - f[2]) / (2.0 * h)
    return out


def derivative_rr(
    grid: AxisymmetricGrid,
    field: npt.ArrayLike,
    *,
    even_at_axis: bool = False,
) -> FloatArray:
    """Second-order radial second derivative."""

    f = grid.validate_field(field)
    out = np.empty_like(f)
    h2 = grid.dr * grid.dr
    out[1:-1] = (f[2:] - 2.0 * f[1:-1] + f[:-2]) / h2
    out[-1] = (
        2.0 * f[-1] - 5.0 * f[-2] + 4.0 * f[-3] - f[-4]
    ) / h2
    if even_at_axis:
        # Even reflection supplies the ghost value f(-h)=f(h).
        out[0] = 2.0 * (f[1] - f[0]) / h2
    else:
        out[0] = (2.0 * f[0] - 5.0 * f[1] + 4.0 * f[2] - f[3]) / h2
    return out


def derivative_z(grid: AxisymmetricGrid, field: npt.ArrayLike) -> FloatArray:
    """Second-order axial first derivative."""

    f = grid.validate_field(field)
    h = grid.dz
    if grid.periodic_z:
        return (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2.0 * h)
    out = np.empty_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * h)
    out[:, 0] = (-3.0 * f[:, 0] + 4.0 * f[:, 1] - f[:, 2]) / (2.0 * h)
    out[:, -1] = (
        3.0 * f[:, -1] - 4.0 * f[:, -2] + f[:, -3]
    ) / (2.0 * h)
    return out


def derivative_zz(grid: AxisymmetricGrid, field: npt.ArrayLike) -> FloatArray:
    """Second-order axial second derivative."""

    f = grid.validate_field(field)
    h2 = grid.dz * grid.dz
    if grid.periodic_z:
        return (
            np.roll(f, -1, axis=1) - 2.0 * f + np.roll(f, 1, axis=1)
        ) / h2
    out = np.empty_like(f)
    out[:, 1:-1] = (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / h2
    out[:, 0] = (
        2.0 * f[:, 0] - 5.0 * f[:, 1] + 4.0 * f[:, 2] - f[:, 3]
    ) / h2
    out[:, -1] = (
        2.0 * f[:, -1]
        - 5.0 * f[:, -2]
        + 4.0 * f[:, -3]
        - f[:, -4]
    ) / h2
    return out


def laplacian_5d_formal(
    grid: AxisymmetricGrid,
    field: npt.ArrayLike,
) -> FloatArray:
    """Evaluate ``f_rr + 3 f_r/r + f_zz`` for an axis-even scalar.

    Smooth axis-even fields obey ``lim_(r->0) f_r/r = f_rr(0)``, so the
    radial part at the axis is ``4 f_rr(0)``.
    """

    f = grid.validate_field(field)
    f_r = derivative_r(grid, f, even_at_axis=True)
    f_rr = derivative_rr(grid, f, even_at_axis=True)
    out = f_rr + derivative_zz(grid, f)
    out[1:] += 3.0 * f_r[1:] / grid.r[1:, None]
    out[0] += 3.0 * f_rr[0]
    return out


def recover_velocity(
    grid: AxisymmetricGrid,
    psi1: npt.ArrayLike,
) -> tuple[FloatArray, FloatArray]:
    """Recover physical meridional velocity from ``psi1 = psi^theta/r``.

    The conventions are

    ``u^r = -r psi1_z`` and ``u^z = 2 psi1 + r psi1_r``.
    """

    psi = grid.validate_field(psi1, name="psi1")
    r_column = grid.r[:, None]
    u_r = -r_column * derivative_z(grid, psi)
    u_z = 2.0 * psi + r_column * derivative_r(
        grid, psi, even_at_axis=True
    )
    # Enforce the exact consequence of odd radial parity at the represented
    # axis, avoiding signed-zero noise in serialized candidates.
    u_r[0] = 0.0
    return u_r, u_z


def divergence_physical(
    grid: AxisymmetricGrid,
    u_r: npt.ArrayLike,
    u_z: npt.ArrayLike,
) -> FloatArray:
    """Evaluate the physical 3D cylindrical divergence.

    Away from the axis this is ``u^r_r + u^r/r + u^z_z``.  For a smooth odd
    radial velocity, the axis limit is ``2 u^r_r(0) + u^z_z(0)``.
    """

    radial = grid.validate_field(u_r, name="u_r")
    axial = grid.validate_field(u_z, name="u_z")
    radial_r = derivative_r(grid, radial, even_at_axis=False)
    out = radial_r + derivative_z(grid, axial)
    out[1:] += radial[1:] / grid.r[1:, None]
    out[0] += radial_r[0]
    return out


def elliptic_defect(
    grid: AxisymmetricGrid,
    psi1: npt.ArrayLike,
    omega1: npt.ArrayLike,
) -> FloatArray:
    """Return the defect in ``-L psi1 = omega1`` as ``omega1 + L psi1``."""

    omega = grid.validate_field(omega1, name="omega1")
    return omega + laplacian_5d_formal(grid, psi1)
