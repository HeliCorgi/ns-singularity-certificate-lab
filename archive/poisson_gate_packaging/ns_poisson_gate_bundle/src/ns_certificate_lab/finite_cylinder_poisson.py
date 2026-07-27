r"""Independent finite-cylinder solver for ``-L5 psi = omega``.

This module solves

.. math::

    -(\partial_{rr} + 3 r^{-1}\partial_r + \partial_{zz})\psi = \omega

on ``0 <= r <= R`` with periodic ``z`` and an explicitly supplied Dirichlet
value at ``r=R``.  The cylindrical axis is included.  Smooth even fields obey

.. math::

    (\partial_{rr} + 3 r^{-1}\partial_r)\psi\big|_{r=0}
    = 4\partial_{rr}\psi(0),

which gives the second-order discrete axis row
``8 * (psi[1] - psi[0]) / dr**2``.

The solve is deliberately independent of :mod:`ns_certificate_lab.operators`:
periodicity is diagonalized with a discrete Fourier transform and every radial
mode is solved by a local tridiagonal implementation.  A separate physical-
space stencil is provided for residual auditing.

This is a floating-point numerical primitive, not an interval proof and not a
full-space Green-function solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid

FloatArray = npt.NDArray[np.float64]
ComplexArray = npt.NDArray[np.complex128]


@dataclass(frozen=True)
class PoissonSolveDiagnostics:
    """Diagnostics attached to one finite-cylinder Poisson solve."""

    fourier_matrix_residual_max: float
    physical_cross_stencil_defect_rms: float
    physical_cross_stencil_defect_max: float
    boundary_error_max: float
    condition_numbers: Mapping[int, float]


@dataclass(frozen=True)
class PoissonSolveResult:
    """Solution and numerical diagnostics."""

    psi: FloatArray
    diagnostics: PoissonSolveDiagnostics


def _validate_periodic_grid(grid: AxisymmetricGrid) -> None:
    if not grid.periodic_z:
        raise ValueError("finite-cylinder Poisson solver requires periodic_z=True")
    if grid.z_period is None or not np.isfinite(grid.z_period) or grid.z_period <= 0.0:
        raise ValueError("periodic grid must have a positive finite z_period")


def _validate_boundary(boundary: npt.ArrayLike | float, grid: AxisymmetricGrid) -> FloatArray:
    array = np.asarray(boundary, dtype=np.float64)
    if array.ndim == 0:
        array = np.full(grid.nz, float(array), dtype=np.float64)
    if array.shape != (grid.nz,):
        raise ValueError(f"outer_boundary must have shape {(grid.nz,)}, got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError("outer_boundary must contain only finite values")
    return array


def axial_wavenumbers(grid: AxisymmetricGrid) -> FloatArray:
    """Return angular Fourier wavenumbers matching NumPy's FFT convention."""

    _validate_periodic_grid(grid)
    return 2.0 * np.pi * np.fft.fftfreq(grid.nz, d=grid.dz)


def radial_mode_tridiagonal(
    grid: AxisymmetricGrid,
    *,
    q_squared: float,
) -> tuple[ComplexArray, ComplexArray, ComplexArray]:
    """Assemble the radial tridiagonal system for one axial Fourier mode.

    Arrays use the standard tridiagonal convention: ``lower[j]`` is row
    ``j+1``, column ``j`` and ``upper[j]`` is row ``j``, column ``j+1``.
    The final row is the outer Dirichlet identity equation.
    """

    _validate_periodic_grid(grid)
    if not np.isfinite(q_squared) or q_squared < 0.0:
        raise ValueError("q_squared must be finite and non-negative")

    nr = grid.nr
    dr = grid.dr
    inv_dr2 = 1.0 / (dr * dr)

    lower = np.zeros(nr - 1, dtype=np.complex128)
    diagonal = np.zeros(nr, dtype=np.complex128)
    upper = np.zeros(nr - 1, dtype=np.complex128)

    # Axis limit: -4 psi_rr + q^2 psi, with
    # psi_rr(0) ~= 2(psi_1 - psi_0)/dr^2.
    diagonal[0] = 8.0 * inv_dr2 + q_squared
    upper[0] = -8.0 * inv_dr2

    for i in range(1, nr - 1):
        r_i = grid.r[i]
        radial_skew = 3.0 / (2.0 * r_i * dr)
        lower[i - 1] = -inv_dr2 + radial_skew
        diagonal[i] = 2.0 * inv_dr2 + q_squared
        upper[i] = -inv_dr2 - radial_skew

    # Explicit outer Dirichlet row.  The coupling from row nr-2 to the
    # boundary remains in upper[nr-2]; the boundary row itself has no lower
    # coupling.
    lower[-1] = 0.0
    diagonal[-1] = 1.0
    return lower, diagonal, upper


def assemble_radial_mode_matrix(
    grid: AxisymmetricGrid,
    *,
    q_squared: float,
) -> ComplexArray:
    """Return a dense radial matrix for auditing and condition estimates."""

    lower, diagonal, upper = radial_mode_tridiagonal(grid, q_squared=q_squared)
    matrix = np.diag(diagonal)
    matrix += np.diag(lower, k=-1)
    matrix += np.diag(upper, k=1)
    return np.asarray(matrix, dtype=np.complex128)


def _solve_tridiagonal(
    lower: ComplexArray,
    diagonal: ComplexArray,
    upper: ComplexArray,
    rhs: ComplexArray,
) -> ComplexArray:
    """Solve a complex tridiagonal system by guarded Thomas elimination."""

    n = diagonal.size
    if lower.shape != (n - 1,) or upper.shape != (n - 1,) or rhs.shape != (n,):
        raise ValueError("incompatible tridiagonal array shapes")

    b = np.array(diagonal, dtype=np.complex128, copy=True)
    c = np.array(upper, dtype=np.complex128, copy=True)
    d = np.array(rhs, dtype=np.complex128, copy=True)
    a = np.asarray(lower, dtype=np.complex128)

    scale = max(1.0, float(np.max(np.abs(b))), float(np.max(np.abs(a))), float(np.max(np.abs(c))))
    pivot_floor = 64.0 * np.finfo(np.float64).eps * scale

    for i in range(1, n):
        pivot = b[i - 1]
        if abs(pivot) <= pivot_floor:
            raise np.linalg.LinAlgError(f"near-zero tridiagonal pivot at row {i - 1}")
        multiplier = a[i - 1] / pivot
        b[i] -= multiplier * c[i - 1]
        d[i] -= multiplier * d[i - 1]

    if abs(b[-1]) <= pivot_floor:
        raise np.linalg.LinAlgError(f"near-zero tridiagonal pivot at row {n - 1}")

    solution = np.empty(n, dtype=np.complex128)
    solution[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        if abs(b[i]) <= pivot_floor:
            raise np.linalg.LinAlgError(f"near-zero tridiagonal pivot at row {i}")
        solution[i] = (d[i] - c[i] * solution[i + 1]) / b[i]
    return solution


def apply_finite_cylinder_system(
    psi: npt.ArrayLike,
    grid: AxisymmetricGrid,
) -> FloatArray:
    """Apply the discrete PDE rows and outer Dirichlet identity row.

    Rows ``0 .. nr-2`` contain the second-order approximation to ``-L5 psi``.
    The final radial row returns ``psi(R,z)`` itself, matching the identity row
    used by the solver.  This physical-space implementation does not invoke the
    Fourier assembly used in :func:`solve_finite_cylinder_poisson`.
    """

    _validate_periodic_grid(grid)
    field = grid.validate_field(psi, name="psi")
    dr = grid.dr
    dz = grid.dz
    inv_dr2 = 1.0 / (dr * dr)
    inv_dz2 = 1.0 / (dz * dz)

    dzz = (np.roll(field, -1, axis=1) - 2.0 * field + np.roll(field, 1, axis=1)) * inv_dz2
    result = np.empty_like(field)

    radial_axis = 8.0 * (field[1, :] - field[0, :]) * inv_dr2
    result[0, :] = -(radial_axis + dzz[0, :])

    r = grid.r[1:-1, None]
    drr = (field[2:, :] - 2.0 * field[1:-1, :] + field[:-2, :]) * inv_dr2
    dr_centered = (field[2:, :] - field[:-2, :]) / (2.0 * dr)
    radial = drr + 3.0 * dr_centered / r
    result[1:-1, :] = -(radial + dzz[1:-1, :])

    result[-1, :] = field[-1, :]
    return result


def _condition_number_sample(
    grid: AxisymmetricGrid,
    *,
    mode_indices: tuple[int, ...],
) -> dict[int, float]:
    q = axial_wavenumbers(grid)
    result: dict[int, float] = {}
    for raw_index in mode_indices:
        index = int(raw_index) % grid.nz
        matrix = assemble_radial_mode_matrix(grid, q_squared=float(q[index] * q[index]))
        result[index] = float(np.linalg.cond(matrix))
    return result


def solve_finite_cylinder_poisson(
    omega: npt.ArrayLike,
    grid: AxisymmetricGrid,
    *,
    outer_boundary: npt.ArrayLike | float,
    condition_mode_indices: tuple[int, ...] = (0, 1),
) -> PoissonSolveResult:
    """Solve ``-L5 psi = omega`` with periodic ``z`` and outer Dirichlet data.

    Parameters
    ----------
    omega:
        Right-hand side sampled on the full ``(nr, nz)`` grid.  Values on the
        outer row are ignored because that row is replaced by the Dirichlet
        equation.
    grid:
        Uniform axisymmetric grid including ``r=0`` and ``r=R``.  ``z`` must be
        periodic and omit its right endpoint.
    outer_boundary:
        Scalar or length-``nz`` array specifying ``psi(R,z)``.
    condition_mode_indices:
        Fourier mode indices for which dense 2-norm condition estimates are
        recorded.  This diagnostic can be disabled with ``()``.
    """

    _validate_periodic_grid(grid)
    rhs_physical = grid.validate_field(omega, name="omega")
    boundary = _validate_boundary(outer_boundary, grid)

    rhs_hat = np.fft.fft(rhs_physical, axis=1)
    boundary_hat = np.fft.fft(boundary)
    psi_hat = np.empty_like(rhs_hat, dtype=np.complex128)
    q = axial_wavenumbers(grid)

    fourier_residual_max = 0.0
    for mode in range(grid.nz):
        lower, diagonal, upper = radial_mode_tridiagonal(
            grid,
            q_squared=float(q[mode] * q[mode]),
        )
        mode_rhs = np.array(rhs_hat[:, mode], dtype=np.complex128, copy=True)
        mode_rhs[-1] = boundary_hat[mode]
        mode_solution = _solve_tridiagonal(lower, diagonal, upper, mode_rhs)
        psi_hat[:, mode] = mode_solution

        check = diagonal * mode_solution
        check[1:] += lower * mode_solution[:-1]
        check[:-1] += upper * mode_solution[1:]
        fourier_residual_max = max(
            fourier_residual_max,
            float(np.max(np.abs(check - mode_rhs))),
        )

    psi_complex = np.fft.ifft(psi_hat, axis=1)
    imaginary_max = float(np.max(np.abs(psi_complex.imag)))
    real_scale = max(1.0, float(np.max(np.abs(psi_complex.real))))
    if imaginary_max > 2048.0 * np.finfo(np.float64).eps * real_scale:
        raise ArithmeticError(
            "inverse FFT produced an unexpectedly large imaginary component: "
            f"{imaginary_max:.3e}"
        )
    psi = np.asarray(psi_complex.real, dtype=np.float64)

    system_rhs = np.array(rhs_physical, dtype=np.float64, copy=True)
    system_rhs[-1, :] = boundary
    physical_residual = apply_finite_cylinder_system(psi, grid) - system_rhs
    interior_residual = physical_residual[:-1, :]
    boundary_error = psi[-1, :] - boundary

    diagnostics = PoissonSolveDiagnostics(
        fourier_matrix_residual_max=fourier_residual_max,
        physical_cross_stencil_defect_rms=float(np.sqrt(np.mean(interior_residual * interior_residual))),
        physical_cross_stencil_defect_max=float(np.max(np.abs(interior_residual))),
        boundary_error_max=float(np.max(np.abs(boundary_error))),
        condition_numbers=_condition_number_sample(
            grid,
            mode_indices=condition_mode_indices,
        ),
    )
    return PoissonSolveResult(psi=psi, diagnostics=diagnostics)
