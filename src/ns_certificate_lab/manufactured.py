"""Independent analytic fields for method-of-manufactured-solutions tests.

No finite-difference routine is called here.  The explicit derivatives provide
an independent oracle for velocity recovery, the elliptic relation, and the
forced PDE residual.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .grid import AxisymmetricGrid, FloatArray


@dataclass(frozen=True)
class ManufacturedFields:
    psi1: FloatArray
    omega1: FloatArray
    u1: FloatArray
    u_r: FloatArray
    u_z: FloatArray
    dpsi1_dr: FloatArray
    dpsi1_dz: FloatArray
    du1_dr: FloatArray
    du1_dz: FloatArray
    du1_dt: FloatArray
    domega1_dr: FloatArray
    domega1_dz: FloatArray
    domega1_dt: FloatArray
    laplacian_u1: FloatArray
    laplacian_omega1: FloatArray

    def forcing_u1(self, viscosity: float) -> FloatArray:
        """Exact forcing making the displayed fields solve the ``u1`` PDE."""

        return (
            self.du1_dt
            + self.u_r * self.du1_dr
            + self.u_z * self.du1_dz
            - 2.0 * self.u1 * self.dpsi1_dz
            - viscosity * self.laplacian_u1
        )

    def forcing_omega1(self, viscosity: float) -> FloatArray:
        """Exact forcing making the displayed fields solve the ``omega1`` PDE."""

        # For u1 = b q sin(2z), d_z(u1^2) = 2 b^2 q^2 sin(4z).
        swirl_source = self._swirl_square_z()
        return (
            self.domega1_dt
            + self.u_r * self.domega1_dr
            + self.u_z * self.domega1_dz
            - swirl_source
            - viscosity * self.laplacian_omega1
        )

    def _swirl_square_z(self) -> FloatArray:
        # Recover it from two independently stored exact factors.  Division is
        # deliberately avoided at zeros of sin(2z).
        return 2.0 * self._b_squared_q_squared() * self._sin_four_z()

    # Private arrays are attached by the factory despite the frozen public
    # dataclass.  Keeping these helpers out of the public artifact avoids
    # mistaking manufactured-only factors for candidate fields.
    def _b_squared_q_squared(self) -> FloatArray:
        return object.__getattribute__(self, "_mms_b2q2")

    def _sin_four_z(self) -> FloatArray:
        return object.__getattribute__(self, "_mms_sin4z")


def manufactured_fields(
    grid: AxisymmetricGrid,
    *,
    time: float = 0.2,
) -> ManufacturedFields:
    r"""Return smooth periodic fields satisfying the elliptic relation exactly.

    The analytic definitions are

    .. math::

       q(r)&=(1-r^2)^2,\quad a=e^{-t},\quad b=e^{-2t},\\
       \psi_1&=a q(r)\cos z,\\
       \omega_1&=a(17-26r^2+r^4)\cos z=-L\psi_1,\\
       u_1&=b q(r)\sin(2z).

    They are not asserted to solve the unforced Navier--Stokes equations.
    Exact forcing is supplied for convergence tests.
    """

    if not np.isfinite(time):
        raise ValueError("time must be finite")
    r, z = grid.mesh()
    r2 = r * r
    r4 = r2 * r2
    q = 1.0 - 2.0 * r2 + r4
    q_r = -4.0 * r + 4.0 * r * r2
    a = float(np.exp(-time))
    b = float(np.exp(-2.0 * time))
    sin_z = np.sin(z)
    cos_z = np.cos(z)
    sin_2z = np.sin(2.0 * z)
    cos_2z = np.cos(2.0 * z)

    psi1 = a * q * cos_z
    dpsi1_dr = a * q_r * cos_z
    dpsi1_dz = -a * q * sin_z
    radial_omega = 17.0 - 26.0 * r2 + r4
    omega1 = a * radial_omega * cos_z
    u1 = b * q * sin_2z

    u_r = a * r * q * sin_z
    u_z = a * (2.0 - 8.0 * r2 + 6.0 * r4) * cos_z
    du1_dr = b * q_r * sin_2z
    du1_dz = 2.0 * b * q * cos_2z
    du1_dt = -2.0 * u1
    domega1_dr = a * (-52.0 * r + 4.0 * r * r2) * cos_z
    domega1_dz = -a * radial_omega * sin_z
    domega1_dt = -omega1
    laplacian_u1 = b * (-20.0 + 32.0 * r2 - 4.0 * r4) * sin_2z
    laplacian_omega1 = a * (-225.0 + 50.0 * r2 - r4) * cos_z

    fields = ManufacturedFields(
        psi1=psi1,
        omega1=omega1,
        u1=u1,
        u_r=u_r,
        u_z=u_z,
        dpsi1_dr=dpsi1_dr,
        dpsi1_dz=dpsi1_dz,
        du1_dr=du1_dr,
        du1_dz=du1_dz,
        du1_dt=du1_dt,
        domega1_dr=domega1_dr,
        domega1_dz=domega1_dz,
        domega1_dt=domega1_dt,
        laplacian_u1=laplacian_u1,
        laplacian_omega1=laplacian_omega1,
    )
    object.__setattr__(fields, "_mms_b2q2", b * b * q * q)
    object.__setattr__(fields, "_mms_sin4z", np.sin(4.0 * z))
    return fields
