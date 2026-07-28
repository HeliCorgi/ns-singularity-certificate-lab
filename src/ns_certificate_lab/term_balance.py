"""PDE term-balance diagnostics for axisymmetric swirl checkpoints.

This module classifies *resolved numerical balances* in the transformed
``u1`` and ``omega1`` equations.  It does not prove an asymptotic regime.  In
particular, sparse checkpoint differentiation can dominate the reported
residual, so every report includes a time-derivative sensitivity indicator.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
import numpy.typing as npt

from .critical_l3 import axisymmetric_volume_weights
from .grid import AxisymmetricGrid, FloatArray
from .pde import EquationTerms, compute_pde_residuals


@dataclass(frozen=True)
class BalanceMetrics:
    """Norms and cancellation ratios for one transformed PDE equation."""

    time_rms: float
    advection_rms: float
    source_rms: float
    diffusion_rms: float
    inviscid_block_rms: float
    residual_rms: float
    normalized_residual: float
    advection_source_cancellation: float
    time_derivative_sensitivity: float
    dominant_balance: str


@dataclass(frozen=True)
class SnapshotBalanceReport:
    """Global and critical-core balances at one checkpoint."""

    time: float
    viscosity: float
    core_fraction: float
    u1_global: BalanceMetrics
    omega1_global: BalanceMetrics
    u1_core: BalanceMetrics
    omega1_core: BalanceMetrics


@dataclass(frozen=True)
class BalanceGateDecision:
    """Pre-registered decision for promoting a numerical balance claim."""

    passed: bool
    reasons: tuple[str, ...]


def balance_gate(
    metrics: BalanceMetrics,
    *,
    maximum_normalized_residual: float = 0.10,
    maximum_time_derivative_sensitivity: float = 0.20,
) -> BalanceGateDecision:
    """Reject balance labels unsupported by trajectory resolution.

    The defaults are deliberately strict and may be changed only before a run,
    not after inspecting its result.
    """

    thresholds = (
        float(maximum_normalized_residual),
        float(maximum_time_derivative_sensitivity),
    )
    if any(not math.isfinite(value) or value < 0.0 for value in thresholds):
        raise ValueError("gate thresholds must be finite and nonnegative")
    reasons: list[str] = []
    if metrics.normalized_residual > thresholds[0]:
        reasons.append("normalized_residual")
    if metrics.time_derivative_sensitivity > thresholds[1]:
        reasons.append("time_derivative_sensitivity")
    if metrics.dominant_balance in {"unresolved_or_other", "zero"}:
        reasons.append("unresolved_classification")
    return BalanceGateDecision(passed=not reasons, reasons=tuple(reasons))


def _validate_times(times: Sequence[float], count: int) -> FloatArray:
    values = np.asarray(times, dtype=np.float64)
    if values.shape != (count,):
        raise ValueError(f"times must have shape ({count},)")
    if count < 2:
        raise ValueError("at least two snapshots are required")
    if not np.all(np.isfinite(values)) or np.any(np.diff(values) <= 0.0):
        raise ValueError("times must be finite and strictly increasing")
    return values


def _lagrange_derivative_weights(nodes: FloatArray, x0: float) -> FloatArray:
    """Derivative weights of the Lagrange interpolant at ``x0``."""

    count = int(nodes.size)
    weights = np.zeros(count, dtype=np.float64)
    for j in range(count):
        total = 0.0
        for m in range(count):
            if m == j:
                continue
            product = 1.0 / (nodes[j] - nodes[m])
            for ell in range(count):
                if ell == j or ell == m:
                    continue
                product *= (x0 - nodes[ell]) / (nodes[j] - nodes[ell])
            total += product
        weights[j] = total
    return weights


def time_derivatives(
    times: Sequence[float],
    fields: Sequence[npt.ArrayLike],
) -> tuple[FloatArray, ...]:
    """Differentiate a nonuniform time series with local polynomial stencils.

    Two snapshots use a first-order secant.  Three or more use local three-point
    quadratic interpolation, including one-sided endpoint stencils.
    """

    arrays = [np.asarray(field, dtype=np.float64) for field in fields]
    if not arrays:
        raise ValueError("fields must be nonempty")
    shape = arrays[0].shape
    if any(array.shape != shape for array in arrays):
        raise ValueError("all fields must have the same shape")
    if any(not np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("fields must be finite")
    values = _validate_times(times, len(arrays))
    if len(arrays) == 2:
        derivative = (arrays[1] - arrays[0]) / (values[1] - values[0])
        return derivative.copy(), derivative.copy()

    output: list[FloatArray] = []
    for index in range(len(arrays)):
        if index == 0:
            chosen = np.array([0, 1, 2], dtype=int)
        elif index == len(arrays) - 1:
            chosen = np.array([len(arrays) - 3, len(arrays) - 2, len(arrays) - 1])
        else:
            chosen = np.array([index - 1, index, index + 1], dtype=int)
        weights = _lagrange_derivative_weights(values[chosen], float(values[index]))
        derivative = sum(
            float(weight) * arrays[int(position)]
            for weight, position in zip(weights, chosen, strict=True)
        )
        output.append(np.asarray(derivative, dtype=np.float64))
    return tuple(output)


def derivative_sensitivity(
    times: Sequence[float],
    fields: Sequence[npt.ArrayLike],
    derivatives: Sequence[npt.ArrayLike] | None = None,
) -> tuple[FloatArray, ...]:
    """Return pointwise disagreement with the nearest two-point secant.

    This is not an error bound.  It flags when checkpoint spacing is too coarse
    for a stable term-balance interpretation.
    """

    arrays = [np.asarray(field, dtype=np.float64) for field in fields]
    values = _validate_times(times, len(arrays))
    derived = time_derivatives(values, arrays) if derivatives is None else tuple(
        np.asarray(item, dtype=np.float64) for item in derivatives
    )
    if len(derived) != len(arrays):
        raise ValueError("derivatives length must match fields")
    output: list[FloatArray] = []
    for index, derivative in enumerate(derived):
        if index == 0:
            other = 1
        elif index == len(arrays) - 1:
            other = index - 1
        else:
            left_dt = values[index] - values[index - 1]
            right_dt = values[index + 1] - values[index]
            other = index - 1 if left_dt <= right_dt else index + 1
        secant = (arrays[other] - arrays[index]) / (values[other] - values[index])
        output.append(np.asarray(derivative - secant, dtype=np.float64))
    return tuple(output)


def _weighted_rms(
    grid: AxisymmetricGrid,
    field: npt.ArrayLike,
    mask: npt.ArrayLike | None,
) -> float:
    values = grid.validate_field(field, name="term")
    weights = axisymmetric_volume_weights(grid)
    if mask is not None:
        selected = np.asarray(mask, dtype=bool)
        if selected.shape != grid.shape:
            raise ValueError("mask must match the grid")
        weights = weights * selected
    denominator = float(np.sum(weights))
    if denominator <= 0.0:
        raise ValueError("mask has zero physical volume")
    result = math.sqrt(float(np.sum(values * values * weights)) / denominator)
    if not math.isfinite(result):
        raise FloatingPointError("non-finite weighted RMS")
    return result


def _classification(time: float, inviscid: float, diffusion: float, residual: float) -> str:
    scale = time + inviscid + diffusion
    if scale <= 0.0:
        return "zero"
    if residual / scale > 0.25:
        return "unresolved_or_other"
    values = np.array([time, inviscid, diffusion], dtype=np.float64)
    maximum = float(np.max(values))
    active = values >= maximum / 4.0
    if bool(np.all(active)):
        return "triple_balance"
    if active[0] and active[1] and not active[2]:
        return "time_inviscid"
    if not active[0] and active[1] and active[2]:
        return "quasi_steady_viscous"
    if active[0] and not active[1] and active[2]:
        return "time_diffusion"
    return "single_block_or_mixed"


def summarize_terms(
    grid: AxisymmetricGrid,
    terms: EquationTerms,
    *,
    time_sensitivity: npt.ArrayLike,
    mask: npt.ArrayLike | None = None,
) -> BalanceMetrics:
    """Summarize one equation using physical-volume-weighted RMS norms."""

    advection = terms.advection_r + terms.advection_z
    inviscid = advection - terms.source
    t_norm = _weighted_rms(grid, terms.time, mask)
    a_norm = _weighted_rms(grid, advection, mask)
    s_norm = _weighted_rms(grid, terms.source, mask)
    d_norm = _weighted_rms(grid, terms.diffusion, mask)
    i_norm = _weighted_rms(grid, inviscid, mask)
    r_norm = _weighted_rms(grid, terms.residual, mask)
    sensitivity_norm = _weighted_rms(grid, time_sensitivity, mask)
    total = t_norm + a_norm + s_norm + d_norm
    epsilon = np.finfo(np.float64).tiny
    return BalanceMetrics(
        time_rms=t_norm,
        advection_rms=a_norm,
        source_rms=s_norm,
        diffusion_rms=d_norm,
        inviscid_block_rms=i_norm,
        residual_rms=r_norm,
        normalized_residual=r_norm / max(total, epsilon),
        advection_source_cancellation=i_norm / max(a_norm + s_norm, epsilon),
        time_derivative_sensitivity=sensitivity_norm / max(t_norm, epsilon),
        dominant_balance=_classification(t_norm, i_norm, d_norm, r_norm),
    )


def critical_core_mask(
    grid: AxisymmetricGrid,
    *,
    center_z: float,
    radial_scale: float,
    axial_scale: float,
    radius_multiplier: float = 2.0,
) -> FloatArray:
    """Elliptic core mask based on critical-density RMS scales."""

    if not math.isfinite(radius_multiplier) or radius_multiplier <= 0.0:
        raise ValueError("radius_multiplier must be positive and finite")
    if radial_scale <= 0.0 or axial_scale <= 0.0:
        raise ValueError("scales must be positive")
    dz = grid.z - float(center_z)
    if grid.periodic_z:
        period = float(grid.z_period)
        dz = (dz + 0.5 * period) % period - 0.5 * period
    radius_squared = (
        (grid.r[:, None] / float(radial_scale)) ** 2
        + (dz[None, :] / float(axial_scale)) ** 2
    )
    return np.asarray(radius_squared <= radius_multiplier**2, dtype=bool)


def snapshot_balance_report(
    grid: AxisymmetricGrid,
    *,
    time: float,
    viscosity: float,
    u1: npt.ArrayLike,
    omega1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    du1_dt: npt.ArrayLike,
    domega1_dt: npt.ArrayLike,
    du1_dt_sensitivity: npt.ArrayLike,
    domega1_dt_sensitivity: npt.ArrayLike,
    core_mask: npt.ArrayLike,
) -> SnapshotBalanceReport:
    """Compute global and core balances for one checkpoint."""

    residuals = compute_pde_residuals(
        grid,
        u1=u1,
        omega1=omega1,
        psi1=psi1,
        du1_dt=du1_dt,
        domega1_dt=domega1_dt,
        viscosity=viscosity,
    )
    mask = np.asarray(core_mask, dtype=bool)
    if mask.shape != grid.shape:
        raise ValueError("core_mask must match grid")
    weights = axisymmetric_volume_weights(grid)
    core_fraction = float(np.sum(weights * mask) / np.sum(weights))
    return SnapshotBalanceReport(
        time=float(time),
        viscosity=float(viscosity),
        core_fraction=core_fraction,
        u1_global=summarize_terms(
            grid, residuals.u1, time_sensitivity=du1_dt_sensitivity
        ),
        omega1_global=summarize_terms(
            grid, residuals.omega1, time_sensitivity=domega1_dt_sensitivity
        ),
        u1_core=summarize_terms(
            grid, residuals.u1, time_sensitivity=du1_dt_sensitivity, mask=mask
        ),
        omega1_core=summarize_terms(
            grid, residuals.omega1, time_sensitivity=domega1_dt_sensitivity, mask=mask
        ),
    )
