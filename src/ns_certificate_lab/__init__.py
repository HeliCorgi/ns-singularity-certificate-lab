"""Auditable numerical primitives for axisymmetric Navier--Stokes studies.

The package deliberately provides diagnostics and candidate serialization, not
a claim of singularity formation.  Arrays use the convention ``(n_r, n_z)``.
"""

from .axis import AxisRegularityReport, check_axis_regularity
from .grid import AxisymmetricGrid
from .operators import (
    divergence_physical,
    elliptic_defect,
    laplacian_5d_formal,
    recover_velocity,
)
from .pde import EquationTerms, PDEResiduals, compute_pde_residuals

__all__ = [
    "AxisRegularityReport",
    "AxisymmetricGrid",
    "EquationTerms",
    "PDEResiduals",
    "check_axis_regularity",
    "compute_pde_residuals",
    "divergence_physical",
    "elliptic_defect",
    "laplacian_5d_formal",
    "recover_velocity",
]

__version__ = "0.1.0"
