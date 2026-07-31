"""Auditable numerical primitives for axisymmetric Navier--Stokes studies.

The package deliberately provides diagnostics and candidate serialization, not
a claim of singularity formation.  Arrays use the convention ``(n_r, n_z)``.
"""

from .axis import (
    AxisRegularityReport,
    check_axis_regularity,
    check_ladder_axis_regularity,
    required_axis_derivative_order,
)
from .critical_l3 import (
    CriticalGrowthReport,
    CriticalL3Report,
    analyze_critical_l3,
    anisotropic_l3_scaling_factor,
    axisymmetric_l3_cubed,
    axisymmetric_l3_norm,
    critical_growth_series,
    state_critical_l3,
)
from .grid import AxisymmetricGrid
from .norm_transfer import (
    CartesianCrossCheck,
    HomogeneousSobolevTransfer,
    SobolevTransfer,
    cartesian_sobolev_norm,
    five_dimensional_algebra_warning,
    homogeneous_sobolev_norm_from_generators,
    isometry_defect,
    ladder_norm,
    mu_inner,
    mu_norm,
    rho_inner,
    rho_norm,
    sobolev_norm_from_generators,
)
from .operators import (
    divergence_physical,
    elliptic_defect,
    laplacian_5d_formal,
    recover_velocity,
)
from .pde import EquationTerms, PDEResiduals, compute_pde_residuals

__all__ = [
    "CriticalGrowthReport",
    "CriticalL3Report",
    "analyze_critical_l3",
    "anisotropic_l3_scaling_factor",
    "axisymmetric_l3_cubed",
    "axisymmetric_l3_norm",
    "critical_growth_series",
    "state_critical_l3",
    "AxisRegularityReport",
    "AxisymmetricGrid",
    "CartesianCrossCheck",
    "EquationTerms",
    "HomogeneousSobolevTransfer",
    "PDEResiduals",
    "SobolevTransfer",
    "cartesian_sobolev_norm",
    "check_axis_regularity",
    "check_ladder_axis_regularity",
    "compute_pde_residuals",
    "divergence_physical",
    "elliptic_defect",
    "five_dimensional_algebra_warning",
    "homogeneous_sobolev_norm_from_generators",
    "isometry_defect",
    "ladder_norm",
    "laplacian_5d_formal",
    "mu_inner",
    "mu_norm",
    "recover_velocity",
    "required_axis_derivative_order",
    "rho_inner",
    "rho_norm",
    "sobolev_norm_from_generators",
]

__version__ = "0.1.0"
