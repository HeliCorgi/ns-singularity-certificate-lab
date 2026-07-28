"""Critical-``L^3`` diagnostics for axisymmetric velocity fields.

The endpoint regularity obstruction makes the global ``L^3`` norm a required
quantity for any whole-space singularity study.  This module provides finite-
volume diagnostics only.  It does not turn a finite-cylinder computation into a
whole-space result, and it does not certify truncation or floating-point error.

Arrays use the repository convention ``(n_r, n_z)``.  The physical volume
measure is ``2*pi*r dr dz``.  Periodic ``z`` grids omit the endpoint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class CriticalL3Report:
    """Critical-norm and concentration diagnostics for one snapshot.

    ``radial_scale`` and ``axial_scale`` are RMS widths of the critical density
    ``|u|^3``.  The radial center is fixed to the symmetry axis.  This is the
    appropriate diagnostic for an axis-focused Hou-like core; it is not a
    point-centered diagnostic for an off-axis vortex ring.
    """

    l3_cubed: float
    l3_norm: float
    max_velocity: float
    center_z: float
    center_resultant_ratio: float
    radial_scale: float
    axial_scale: float
    isotropic_scale: float
    critical_scaling_product: float
    outer_radial_fraction: float
    shell_edges: tuple[float, ...]
    shell_masses: tuple[float, ...]
    shell_fractions: tuple[float, ...]
    shell_entropy: float
    effective_shell_count: float
    maximum_shell_fraction: float


@dataclass(frozen=True)
class CriticalGrowthReport:
    """Finite-difference growth indicators for an ordered time series.

    These slopes are diagnostics, not blow-up exponents.  They are meaningful
    only after space, time, domain and solver refinement.
    """

    times: tuple[float, ...]
    l3_norms: tuple[float, ...]
    scaling_products: tuple[float, ...]
    effective_shell_counts: tuple[float, ...]
    log_l3_slopes: tuple[float, ...]
    log_scaling_product_slopes: tuple[float, ...]
    shell_count_slopes: tuple[float, ...]


def _validated_components(
    grid: AxisymmetricGrid,
    u_r: npt.ArrayLike,
    u_theta: npt.ArrayLike,
    u_z: npt.ArrayLike,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    radial = grid.validate_field(u_r, name="u_r")
    swirl = grid.validate_field(u_theta, name="u_theta")
    axial = grid.validate_field(u_z, name="u_z")
    return radial, swirl, axial


def axisymmetric_volume_weights(grid: AxisymmetricGrid) -> FloatArray:
    """Return tensor-product quadrature weights for ``2*pi*r dr dz``.

    The radial rule is trapezoidal, including the outer wall.  The periodic
    axial rule is the rectangle rule because the endpoint is omitted.  A
    nonperiodic axial grid uses trapezoidal endpoint weights.
    """

    radial_weights = np.full(grid.nr, grid.dr, dtype=np.float64)
    radial_weights[0] *= 0.5
    radial_weights[-1] *= 0.5
    radial_weights *= 2.0 * math.pi * grid.r

    axial_weights = np.full(grid.nz, grid.dz, dtype=np.float64)
    if not grid.periodic_z:
        axial_weights[0] *= 0.5
        axial_weights[-1] *= 0.5
    weights = radial_weights[:, None] * axial_weights[None, :]
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise FloatingPointError("invalid axisymmetric volume weights")
    return weights


def velocity_magnitude(
    grid: AxisymmetricGrid,
    u_r: npt.ArrayLike,
    u_theta: npt.ArrayLike,
    u_z: npt.ArrayLike,
) -> FloatArray:
    """Return ``sqrt((u^r)^2 + (u^theta)^2 + (u^z)^2)``."""

    radial, swirl, axial = _validated_components(grid, u_r, u_theta, u_z)
    magnitude = np.sqrt(radial * radial + swirl * swirl + axial * axial)
    if not np.all(np.isfinite(magnitude)):
        raise FloatingPointError("velocity magnitude is not finite")
    return magnitude


def axisymmetric_l3_cubed(
    grid: AxisymmetricGrid,
    u_r: npt.ArrayLike,
    u_theta: npt.ArrayLike,
    u_z: npt.ArrayLike,
) -> float:
    """Approximate ``int |u|^3 dx`` with the physical cylindrical measure."""

    magnitude = velocity_magnitude(grid, u_r, u_theta, u_z)
    value = float(np.sum(magnitude**3 * axisymmetric_volume_weights(grid)))
    if not math.isfinite(value) or value < 0.0:
        raise FloatingPointError("computed L3 cubed is invalid")
    return value


def axisymmetric_l3_norm(
    grid: AxisymmetricGrid,
    u_r: npt.ArrayLike,
    u_theta: npt.ArrayLike,
    u_z: npt.ArrayLike,
) -> float:
    """Approximate the physical ``L^3`` norm on the represented domain."""

    return axisymmetric_l3_cubed(grid, u_r, u_theta, u_z) ** (1.0 / 3.0)


def anisotropic_l3_scaling_factor(
    amplitude: float,
    radial_scale: float,
    axial_scale: float,
) -> float:
    """Return the exact anisotropic critical factor ``A^3 L_r^2 L_z``."""

    values = (float(amplitude), float(radial_scale), float(axial_scale))
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("amplitude and scales must be positive and finite")
    return values[0] ** 3 * values[1] ** 2 * values[2]


def rescaled_l3_cubed(
    base_l3_cubed: float,
    *,
    amplitude: float,
    radial_scale: float,
    axial_scale: float,
) -> float:
    """Apply the exact change-of-variable factor to a base ``L^3`` integral."""

    base = float(base_l3_cubed)
    if not math.isfinite(base) or base < 0.0:
        raise ValueError("base_l3_cubed must be finite and nonnegative")
    return base * anisotropic_l3_scaling_factor(
        amplitude,
        radial_scale,
        axial_scale,
    )


def _periodic_displacement(grid: AxisymmetricGrid, center: float) -> FloatArray:
    coordinates = grid.z - float(center)
    if not grid.periodic_z:
        return coordinates
    period = float(grid.z_period)
    return (coordinates + 0.5 * period) % period - 0.5 * period


def _critical_center_z(
    grid: AxisymmetricGrid,
    density: FloatArray,
    weights: FloatArray,
) -> tuple[float, float]:
    axial_mass = np.sum(density * weights, axis=0)
    total = float(np.sum(axial_mass))
    if total <= 0.0:
        return float(grid.z[0]), 0.0
    if not grid.periodic_z:
        return float(np.sum(grid.z * axial_mass) / total), 1.0

    period = float(grid.z_period)
    angles = 2.0 * math.pi * (grid.z - grid.z[0]) / period
    cosine = float(np.sum(axial_mass * np.cos(angles)))
    sine = float(np.sum(axial_mass * np.sin(angles)))
    resultant = math.hypot(cosine, sine)
    ratio = resultant / total
    if ratio <= 1.0e-12:
        # Circular mean is undefined for nearly uniform or symmetric mass.
        # Use the maximum-mass column as a deterministic diagnostic center.
        return float(grid.z[int(np.argmax(axial_mass))]), float(ratio)
    angle = math.atan2(sine, cosine) % (2.0 * math.pi)
    center = float(grid.z[0] + period * angle / (2.0 * math.pi))
    return center, float(ratio)


def _entropy_and_effective_count(fractions: FloatArray) -> tuple[float, float]:
    positive = fractions[fractions > 0.0]
    if positive.size == 0:
        return 0.0, 0.0
    entropy = float(-np.sum(positive * np.log(positive)))
    return entropy, float(math.exp(entropy))


def analyze_critical_l3(
    grid: AxisymmetricGrid,
    u_r: npt.ArrayLike,
    u_theta: npt.ArrayLike,
    u_z: npt.ArrayLike,
    *,
    shell_count: int = 10,
    outer_radial_cells: int = 4,
) -> CriticalL3Report:
    """Analyze critical mass, anisotropic widths and dyadic shell spreading.

    Shells are centered on the symmetry axis at the critical-density axial
    center.  The first shell is ``rho < L``; subsequent finite shells are
    ``2^(j-1)L <= rho < 2^j L``; the final shell is the remaining represented
    domain.  Here ``L = (L_r^2 L_z)^(1/3)``.
    """

    if isinstance(shell_count, bool) or shell_count < 2:
        raise ValueError("shell_count must be an integer at least two")
    if isinstance(outer_radial_cells, bool) or not 1 <= outer_radial_cells < grid.nr:
        raise ValueError("outer_radial_cells must lie in [1, nr-1]")

    magnitude = velocity_magnitude(grid, u_r, u_theta, u_z)
    density = magnitude**3
    weights = axisymmetric_volume_weights(grid)
    weighted_density = density * weights
    total = float(np.sum(weighted_density))
    if total <= 0.0:
        raise ValueError("critical diagnostics require a nonzero velocity field")

    center_z, resultant_ratio = _critical_center_z(grid, density, weights)
    axial_displacement = _periodic_displacement(grid, center_z)
    radial_coordinate = grid.r[:, None]
    radial_variance = float(np.sum(radial_coordinate**2 * weighted_density) / total)
    axial_variance = float(
        np.sum(axial_displacement[None, :] ** 2 * weighted_density) / total
    )
    radial_scale = math.sqrt(max(0.0, radial_variance))
    axial_scale = math.sqrt(max(0.0, axial_variance))
    if radial_scale <= 0.0 or axial_scale <= 0.0:
        raise ValueError("critical density has a degenerate RMS scale")

    amplitude = float(np.max(magnitude))
    isotropic_scale = (radial_scale**2 * axial_scale) ** (1.0 / 3.0)
    scaling_product = anisotropic_l3_scaling_factor(
        amplitude,
        radial_scale,
        axial_scale,
    )

    rho = np.sqrt(radial_coordinate**2 + axial_displacement[None, :] ** 2)
    finite_edges = isotropic_scale * 2.0 ** np.arange(shell_count - 1, dtype=np.float64)
    shell_indices = np.searchsorted(finite_edges, rho, side="right")
    masses = np.bincount(
        shell_indices.ravel(),
        weights=weighted_density.ravel(),
        minlength=shell_count,
    ).astype(np.float64)
    fractions = masses / total
    entropy, effective_count = _entropy_and_effective_count(fractions)

    outer_start = grid.nr - int(outer_radial_cells)
    outer_fraction = float(np.sum(weighted_density[outer_start:]) / total)
    edges = tuple(float(value) for value in finite_edges)
    return CriticalL3Report(
        l3_cubed=total,
        l3_norm=total ** (1.0 / 3.0),
        max_velocity=amplitude,
        center_z=center_z,
        center_resultant_ratio=resultant_ratio,
        radial_scale=radial_scale,
        axial_scale=axial_scale,
        isotropic_scale=isotropic_scale,
        critical_scaling_product=scaling_product,
        outer_radial_fraction=outer_fraction,
        shell_edges=edges,
        shell_masses=tuple(float(value) for value in masses),
        shell_fractions=tuple(float(value) for value in fractions),
        shell_entropy=entropy,
        effective_shell_count=effective_count,
        maximum_shell_fraction=float(np.max(fractions)),
    )


def state_critical_l3(
    grid: AxisymmetricGrid,
    *,
    u1: npt.ArrayLike,
    u_r: npt.ArrayLike,
    u_z: npt.ArrayLike,
    shell_count: int = 10,
    outer_radial_cells: int = 4,
) -> CriticalL3Report:
    """Analyze a repository state using ``u^theta = r*u1``."""

    swirl_rate = grid.validate_field(u1, name="u1")
    u_theta = grid.r[:, None] * swirl_rate
    return analyze_critical_l3(
        grid,
        u_r,
        u_theta,
        u_z,
        shell_count=shell_count,
        outer_radial_cells=outer_radial_cells,
    )


def critical_growth_series(
    times: Sequence[float],
    reports: Sequence[CriticalL3Report],
) -> CriticalGrowthReport:
    """Return local logarithmic growth slopes for an ordered report series."""

    if len(times) != len(reports) or len(times) < 2:
        raise ValueError("times and reports must have the same length >= 2")
    moments = np.asarray(times, dtype=np.float64)
    if not np.all(np.isfinite(moments)) or np.any(np.diff(moments) <= 0.0):
        raise ValueError("times must be finite and strictly increasing")

    l3 = np.asarray([item.l3_norm for item in reports], dtype=np.float64)
    products = np.asarray(
        [item.critical_scaling_product for item in reports], dtype=np.float64
    )
    shell_counts = np.asarray(
        [item.effective_shell_count for item in reports], dtype=np.float64
    )
    if np.any(l3 <= 0.0) or np.any(products <= 0.0):
        raise ValueError("reports contain nonpositive critical quantities")
    dt = np.diff(moments)
    return CriticalGrowthReport(
        times=tuple(float(value) for value in moments),
        l3_norms=tuple(float(value) for value in l3),
        scaling_products=tuple(float(value) for value in products),
        effective_shell_counts=tuple(float(value) for value in shell_counts),
        log_l3_slopes=tuple(float(value) for value in np.diff(np.log(l3)) / dt),
        log_scaling_product_slopes=tuple(
            float(value) for value in np.diff(np.log(products)) / dt
        ),
        shell_count_slopes=tuple(float(value) for value in np.diff(shell_counts) / dt),
    )
