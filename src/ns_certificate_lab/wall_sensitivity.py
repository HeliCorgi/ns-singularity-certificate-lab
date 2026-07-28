"""Exact modal wall-correction formulas for the regularized axisymmetric Poisson problem.

For one axial Fourier mode, let ``psi_inf`` solve the free-radial problem

    -[d_rr + (3/r)d_r - k^2] psi_inf = f,

with regularity at ``r=0`` and decay as ``r -> infinity``.  Let ``psi_R`` solve
that equation on ``0 <= r <= R`` with the artificial condition ``psi_R(R)=0``.
Their difference is homogeneous and regular at the axis.

For ``k > 0``:

    psi_R(r) - psi_inf(r)
      = -psi_inf(R) [I1(k r)/r] / [I1(k R)/R].

For ``k = 0``, the regular homogeneous solution is constant:

    psi_R(r) - psi_inf(r) = -psi_inf(R).

If the source vanishes outside ``R``, the free boundary value is also explicit:

    psi_inf(R) = (1/(2 R^2)) int rho^3 f(rho) d rho,               k = 0,

    psi_inf(R) = [K1(kR)/R] int [I1(k rho)/rho] f(rho) rho^3 d rho, k > 0.

These formulas explain the transition from algebraic zero-mode wall influence to
exponentially small fixed-nonzero-mode influence as ``kR`` grows.
"""

from __future__ import annotations

import math

import numpy as np

from .free_space_poisson import _i1e_cephes, _k1e_cephes


def _uniform_trapezoid_weights(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if r.ndim != 1 or r.size < 2:
        raise ValueError("r must be one-dimensional with at least two points")
    if r[0] != 0.0 or np.any(np.diff(r) <= 0.0):
        raise ValueError("r must start at zero and be strictly increasing")
    dr = np.diff(r)
    if not np.allclose(dr, dr[0], rtol=1e-12, atol=1e-14):
        raise ValueError("the prototype requires a uniform radial grid")
    weights = np.full(r.size, float(dr[0]))
    weights[[0, -1]] *= 0.5
    return weights


def regular_mode_ratio(
    r: np.ndarray,
    k: float,
    *,
    outer_radius: float | None = None,
) -> np.ndarray:
    """Return the regular homogeneous mode normalized to one at ``R``.

    For ``k=0`` the ratio is identically one.  For ``k>0`` it is

        [I1(k r)/r] / [I1(k R)/R].

    The axis value is evaluated by the removable limit ``I1(kr)/r -> k/2``.
    """

    r = np.asarray(r, dtype=float)
    if r.ndim != 1 or r.size < 2:
        raise ValueError("r must be a one-dimensional grid")
    if r[0] != 0.0 or np.any(np.diff(r) <= 0.0):
        raise ValueError("r must start at zero and be strictly increasing")
    if not (math.isfinite(k) and k >= 0.0):
        raise ValueError("k must be finite and nonnegative")
    radius = float(r[-1] if outer_radius is None else outer_radius)
    if not (math.isfinite(radius) and radius > 0.0):
        raise ValueError("outer_radius must be finite and positive")
    if np.any(r > radius + 1e-14):
        raise ValueError("all radial nodes must lie within outer_radius")

    if k == 0.0:
        return np.ones_like(r)

    x = k * r
    x_radius = k * radius
    denominator_scaled = float(_i1e_cephes(np.array([x_radius]))[0])
    if denominator_scaled <= 0.0:
        raise ArithmeticError("invalid scaled I1 value at outer boundary")

    ratio = np.empty_like(r)
    positive = r > 0.0
    ratio[positive] = (
        (radius / r[positive])
        * (_i1e_cephes(x[positive]) / denominator_scaled)
        * np.exp(x[positive] - x_radius)
    )
    ratio[0] = k * radius * math.exp(-x_radius) / (2.0 * denominator_scaled)
    return ratio


def predicted_zero_dirichlet_correction(
    r: np.ndarray,
    k: float,
    free_boundary_amplitude: float,
    *,
    outer_radius: float | None = None,
) -> np.ndarray:
    """Predict ``psi_R - psi_inf`` for one Fourier-mode amplitude."""

    if not math.isfinite(free_boundary_amplitude):
        raise ValueError("free_boundary_amplitude must be finite")
    return -float(free_boundary_amplitude) * regular_mode_ratio(
        r,
        k,
        outer_radius=outer_radius,
    )


def zero_mode_moment(source_radial: np.ndarray, r: np.ndarray) -> float:
    """Return ``int rho^3 f(rho) d rho`` by composite trapezoid."""

    source = np.asarray(source_radial, dtype=float)
    r = np.asarray(r, dtype=float)
    if source.shape != r.shape:
        raise ValueError("source_radial and r must have the same shape")
    if not np.all(np.isfinite(source)):
        raise ValueError("source_radial must be finite")
    weights = _uniform_trapezoid_weights(r)
    return float(np.sum(weights * r**3 * source))


def predicted_free_boundary_from_source(
    source_radial: np.ndarray,
    r: np.ndarray,
    k: float,
) -> float:
    """Predict the free-radial solution amplitude at the last radial node.

    The supplied source is interpreted as zero beyond ``r[-1]``.
    """

    source = np.asarray(source_radial, dtype=float)
    r = np.asarray(r, dtype=float)
    if source.shape != r.shape:
        raise ValueError("source_radial and r must have the same shape")
    if not np.all(np.isfinite(source)):
        raise ValueError("source_radial must be finite")
    if not (math.isfinite(k) and k >= 0.0):
        raise ValueError("k must be finite and nonnegative")
    radius = float(r[-1])
    weights = _uniform_trapezoid_weights(r)

    if k == 0.0:
        return zero_mode_moment(source, r) / (2.0 * radius**2)

    x = k * r
    i1_over_r = np.empty_like(r)
    positive = r > 0.0
    i1_over_r[positive] = (
        _i1e_cephes(x[positive]) * np.exp(x[positive]) / r[positive]
    )
    i1_over_r[0] = 0.5 * k
    coefficient = float(np.sum(weights * r**3 * source * i1_over_r))

    x_radius = k * radius
    k1_over_radius = (
        float(_k1e_cephes(np.array([x_radius]))[0])
        * math.exp(-x_radius)
        / radius
    )
    return coefficient * k1_over_radius
