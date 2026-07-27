"""Term-by-term residuals for the axisymmetric swirl formulation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid, FloatArray
from .operators import (
    derivative_r,
    derivative_z,
    laplacian_5d_formal,
    recover_velocity,
)


@dataclass(frozen=True)
class EquationTerms:
    """Signed decomposition of one equation.

    ``residual = time + advection_r + advection_z - source - diffusion``.
    The source and diffusion members are stored with their right-hand-side
    signs, making each physical contribution inspectable before assembly.
    """

    time: FloatArray
    advection_r: FloatArray
    advection_z: FloatArray
    source: FloatArray
    diffusion: FloatArray
    residual: FloatArray

    def as_dict(self) -> dict[str, FloatArray]:
        return {
            "time": self.time,
            "advection_r": self.advection_r,
            "advection_z": self.advection_z,
            "source": self.source,
            "diffusion": self.diffusion,
            "residual": self.residual,
        }


@dataclass(frozen=True)
class PDEResiduals:
    """Residual terms for the ``u1`` and ``omega1`` equations."""

    u1: EquationTerms
    omega1: EquationTerms
    u_r: FloatArray
    u_z: FloatArray


def _assemble(
    *,
    time: FloatArray,
    advection_r: FloatArray,
    advection_z: FloatArray,
    source: FloatArray,
    diffusion: FloatArray,
) -> EquationTerms:
    residual = time + advection_r + advection_z - source - diffusion
    return EquationTerms(
        time=time,
        advection_r=advection_r,
        advection_z=advection_z,
        source=source,
        diffusion=diffusion,
        residual=residual,
    )


def compute_pde_residuals(
    grid: AxisymmetricGrid,
    *,
    u1: npt.ArrayLike,
    omega1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    du1_dt: npt.ArrayLike,
    domega1_dt: npt.ArrayLike,
    viscosity: float,
) -> PDEResiduals:
    r"""Compute all terms of the transformed axisymmetric equations.

    With ``L = d_rr + 3 d_r/r + d_zz``, the audited convention is

    .. math::

       u_{1,t}+u^r u_{1,r}+u^z u_{1,z}
          &=2u_1\psi_{1,z}+\nu L u_1,\\
       \omega_{1,t}+u^r\omega_{1,r}+u^z\omega_{1,z}
          &=\partial_z(u_1^2)+\nu L\omega_1.

    Time derivatives are inputs so that stored trajectories and independent
    differentiators can be audited without hiding a time scheme in this
    routine.
    """

    if not np.isfinite(viscosity) or viscosity < 0.0:
        raise ValueError("viscosity must be nonnegative and finite")
    swirl = grid.validate_field(u1, name="u1")
    vorticity = grid.validate_field(omega1, name="omega1")
    stream = grid.validate_field(psi1, name="psi1")
    swirl_t = grid.validate_field(du1_dt, name="du1_dt")
    vorticity_t = grid.validate_field(domega1_dt, name="domega1_dt")

    u_r, u_z = recover_velocity(grid, stream)
    u1_r = derivative_r(grid, swirl, even_at_axis=True)
    u1_z = derivative_z(grid, swirl)
    omega1_r = derivative_r(grid, vorticity, even_at_axis=True)
    omega1_z = derivative_z(grid, vorticity)

    u1_terms = _assemble(
        time=swirl_t.copy(),
        advection_r=u_r * u1_r,
        advection_z=u_z * u1_z,
        source=2.0 * swirl * derivative_z(grid, stream),
        diffusion=float(viscosity) * laplacian_5d_formal(grid, swirl),
    )
    omega1_terms = _assemble(
        time=vorticity_t.copy(),
        advection_r=u_r * omega1_r,
        advection_z=u_z * omega1_z,
        source=derivative_z(grid, swirl * swirl),
        diffusion=float(viscosity) * laplacian_5d_formal(grid, vorticity),
    )
    return PDEResiduals(
        u1=u1_terms,
        omega1=omega1_terms,
        u_r=u_r,
        u_z=u_z,
    )
