r"""Independent finite-cylinder solver for ``-L5 psi1 = omega1``.

The regularized axisymmetric streamfunction equation is

.. math::

   -\left(\partial_{rr} + \frac{3}{r}\partial_r
           + \partial_{zz}\right)\psi_1 = \omega_1.

This module deliberately does not import the cylindrical finite-difference
operators used elsewhere in the project.  The axial direction is represented
by NumPy Fourier coefficients.  The radial operator is a conservative
finite-volume discretization of

.. math::

   r^{-3}\partial_r(r^3\partial_r\psi_1).

The control volume touching the axis gives the exact discrete limit
``8 * (psi1[1] - psi1[0]) / dr**2``.  Every Fourier mode is solved by a
locally implemented Thomas algorithm; no external sparse solver is required.

The returned ``discrete_residual`` measures the algebraic equation solved by
the Fourier/finite-volume method.  ``pde_residual`` is intentionally computed
through a second, physical-space centered-difference path.  Neither is a
continuum error bound, and analytic manufactured-solution error must be
measured separately.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid, FloatArray

ComplexArray = npt.NDArray[np.complex128]


@dataclass(frozen=True)
class PoissonSolution:
    """Result and diagnostics for one finite-cylinder elliptic solve."""

    psi1: FloatArray
    discrete_residual: FloatArray
    pde_residual: FloatArray
    metadata: dict[str, Any]


def _outer_trace(
    grid: AxisymmetricGrid,
    outer_dirichlet: npt.ArrayLike,
) -> FloatArray:
    raw_trace = np.asarray(outer_dirichlet)
    if np.iscomplexobj(raw_trace) or (
        raw_trace.dtype.kind == "O"
        and any(isinstance(value, complex) for value in raw_trace.flat)
    ):
        raise ValueError("outer_dirichlet must be real-valued")
    trace = np.asarray(outer_dirichlet, dtype=np.float64)
    if trace.ndim == 0:
        trace = np.full(grid.nz, float(trace), dtype=np.float64)
    elif trace.shape == (grid.nz,):
        trace = trace.copy()
    else:
        raise ValueError(
            f"outer_dirichlet must be scalar or have shape ({grid.nz},)"
        )
    if not np.all(np.isfinite(trace)):
        raise ValueError("outer_dirichlet must contain only finite values")
    return trace


def _radial_flux_coefficients(
    grid: AxisymmetricGrid,
) -> tuple[FloatArray, FloatArray]:
    """Return lower/upper positive coefficients of the radial ``L5`` part.

    Only the ``nr - 1`` unknown rows are returned; the outermost radial point
    is prescribed by Dirichlet data.
    """

    radial_unknowns = grid.r[:-1]
    half_lower = np.maximum(radial_unknowns - 0.5 * grid.dr, 0.0)
    half_upper = radial_unknowns + 0.5 * grid.dr
    volumes = 0.25 * (half_upper**4 - half_lower**4)
    lower = half_lower**3 / (grid.dr * volumes)
    upper = half_upper**3 / (grid.dr * volumes)
    # Keep the analytically derived axis row literal, rather than trusting
    # cancellation of small powers in the general finite-volume expression.
    lower[0] = 0.0
    upper[0] = 8.0 / grid.dr**2
    return lower, upper


def _solve_tridiagonal(
    lower: FloatArray,
    diagonal: FloatArray,
    upper: FloatArray,
    right_hand_side: ComplexArray,
) -> ComplexArray:
    """Solve a real tridiagonal matrix against one complex right-hand side."""

    size = diagonal.size
    if (
        lower.shape != (size - 1,)
        or upper.shape != (size - 1,)
        or right_hand_side.shape != (size,)
    ):
        raise ValueError("invalid tridiagonal system shape")
    modified_diagonal = diagonal.astype(np.complex128, copy=True)
    modified_rhs = right_hand_side.astype(np.complex128, copy=True)
    upper_complex = upper.astype(np.complex128, copy=False)
    scale = max(
        1.0,
        float(np.max(np.abs(diagonal))),
        float(np.max(np.abs(lower), initial=0.0)),
        float(np.max(np.abs(upper), initial=0.0)),
    )
    pivot_floor = 32.0 * np.finfo(np.float64).eps * scale

    for row in range(1, size):
        if abs(modified_diagonal[row - 1]) <= pivot_floor:
            raise np.linalg.LinAlgError("unsafe pivot in radial Thomas solve")
        multiplier = lower[row - 1] / modified_diagonal[row - 1]
        modified_diagonal[row] -= multiplier * upper_complex[row - 1]
        modified_rhs[row] -= multiplier * modified_rhs[row - 1]
    if abs(modified_diagonal[-1]) <= pivot_floor:
        raise np.linalg.LinAlgError("unsafe final pivot in radial Thomas solve")

    solution = np.empty(size, dtype=np.complex128)
    solution[-1] = modified_rhs[-1] / modified_diagonal[-1]
    for row in range(size - 2, -1, -1):
        solution[row] = (
            modified_rhs[row] - upper_complex[row] * solution[row + 1]
        ) / modified_diagonal[row]
    return solution


def _zero_mode_condition_number_inf(nr: int) -> float:
    """Return one zero-mode matrix's induced infinity-norm condition number.

    The zero-mode radial matrix is an irreducible M-matrix.  The exact discrete
    solution of ``A y = 1`` is ``y_i = (R**2-r_i**2)/8``, so positivity gives
    ``||A^-1||_inf = R**2/8``.  The axis row has the largest absolute row sum,
    ``||A||_inf = 16/dr**2``.  Therefore the unweighted coordinate-basis
    condition number is exactly ``2*(nr-1)**2``.  It is not a weighted
    coercivity estimate.
    """

    return float(2 * (nr - 1) ** 2)


def _spectral_discrete_residual(
    grid: AxisymmetricGrid,
    *,
    psi1: FloatArray,
    omega1: FloatArray,
    outer_trace: FloatArray,
    lower_flux: FloatArray,
    upper_flux: FloatArray,
) -> tuple[FloatArray, float]:
    """Apply the exact discrete operator used by the solver."""

    psi_hat = np.fft.fft(psi1, axis=1)
    omega_hat = np.fft.fft(omega1, axis=1)
    wave_numbers = 2.0 * np.pi * np.fft.fftfreq(grid.nz, d=grid.dz)
    residual_hat = np.zeros_like(psi_hat)
    last_unknown = grid.nr - 2
    for mode, wave_number in enumerate(wave_numbers):
        radial_values = psi_hat[:, mode]
        applied = (
            (lower_flux + upper_flux + wave_number**2)
            * radial_values[:-1]
        )
        applied[1:] -= lower_flux[1:] * radial_values[:-2]
        applied[:-1] -= upper_flux[:-1] * radial_values[1:-1]
        applied[last_unknown] -= (
            upper_flux[last_unknown] * radial_values[-1]
        )
        residual_hat[:-1, mode] = applied - omega_hat[:-1, mode]
    boundary_hat = psi_hat[-1] - np.fft.fft(outer_trace)
    residual_hat[-1] = boundary_hat
    residual_complex = np.fft.ifft(residual_hat, axis=1)
    imaginary_leakage = float(np.max(np.abs(residual_complex.imag)))
    return np.asarray(residual_complex.real, dtype=np.float64), imaginary_leakage


def independent_physical_poisson_residual(
    grid: AxisymmetricGrid,
    *,
    psi1: npt.ArrayLike,
    omega1: npt.ArrayLike,
    outer_dirichlet: npt.ArrayLike,
) -> FloatArray:
    """Compute ``-L5 psi1 - omega1`` with an independent physical-space path.

    Centered, periodically wrapped second differences are used in ``z``.
    Interior radial rows use the direct non-divergence form.  The axis row
    uses the even-field limit with coefficient eight.  The final row stores
    the Dirichlet trace defect rather than a PDE residual.
    """

    if not grid.periodic_z:
        raise ValueError("the Poisson audit requires a periodic z grid")
    psi = grid.validate_field(psi1, name="psi1")
    omega = grid.validate_field(omega1, name="omega1")
    trace = _outer_trace(grid, outer_dirichlet)
    dz2 = (
        np.roll(psi, -1, axis=1)
        - 2.0 * psi
        + np.roll(psi, 1, axis=1)
    ) / grid.dz**2
    radial = np.empty_like(psi)
    radial[0] = 8.0 * (psi[1] - psi[0]) / grid.dr**2
    r_interior = grid.r[1:-1, None]
    radial[1:-1] = (
        (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) / grid.dr**2
        + 3.0
        * (psi[2:] - psi[:-2])
        / (2.0 * grid.dr * r_interior)
    )
    residual = np.empty_like(psi)
    residual[:-1] = -radial[:-1] - dz2[:-1] - omega[:-1]
    residual[-1] = psi[-1] - trace
    return residual


def solve_streamfunction_poisson(
    grid: AxisymmetricGrid,
    omega1: npt.ArrayLike,
    outer_dirichlet: npt.ArrayLike,
    *,
    estimate_condition: bool = True,
) -> PoissonSolution:
    """Solve ``-L5 psi1 = omega1`` with periodic ``z`` and outer Dirichlet data.

    Parameters
    ----------
    grid:
        Uniform :class:`AxisymmetricGrid`.  ``periodic_z`` must be true.
    omega1:
        Finite nodal right-hand side with shape ``grid.shape``.  The outer
        radial row is not used by the PDE solve because it is a boundary row,
        but it must still be finite.
    outer_dirichlet:
        Finite scalar or one-dimensional array of length ``grid.nz``.
    estimate_condition:
        Whether to include the cached dense zero-mode Euclidean infinity-norm
        condition number in the unweighted coordinate basis.  Time
        integrators may disable this optional diagnostic; all residual arrays
        are still computed.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    if not isinstance(estimate_condition, bool):
        raise TypeError("estimate_condition must be boolean")
    if not grid.periodic_z:
        raise ValueError("solve_streamfunction_poisson requires periodic_z=True")
    raw_omega = np.asarray(omega1)
    if np.iscomplexobj(raw_omega) or (
        raw_omega.dtype.kind == "O"
        and any(isinstance(value, complex) for value in raw_omega.flat)
    ):
        raise ValueError("omega1 must be real-valued")
    omega = grid.validate_field(omega1, name="omega1")
    trace = _outer_trace(grid, outer_dirichlet)
    lower_flux, upper_flux = _radial_flux_coefficients(grid)
    omega_hat = np.fft.fft(omega, axis=1)
    trace_hat = np.fft.fft(trace)
    wave_numbers = 2.0 * np.pi * np.fft.fftfreq(grid.nz, d=grid.dz)
    unknown_count = grid.nr - 1
    psi_hat = np.empty(grid.shape, dtype=np.complex128)
    psi_hat[-1] = trace_hat

    for mode, wave_number in enumerate(wave_numbers):
        diagonal = lower_flux + upper_flux + wave_number**2
        lower = -lower_flux[1:]
        upper = -upper_flux[:-1]
        right_hand_side = omega_hat[:unknown_count, mode].copy()
        right_hand_side[-1] += upper_flux[-1] * trace_hat[mode]
        psi_hat[:unknown_count, mode] = _solve_tridiagonal(
            lower,
            diagonal,
            upper,
            right_hand_side,
        )

    psi_complex = np.fft.ifft(psi_hat, axis=1)
    solution_imaginary_leakage = float(np.max(np.abs(psi_complex.imag)))
    solution_scale = max(1.0, float(np.max(np.abs(psi_complex.real))))
    solution_relative_imaginary_leakage = (
        solution_imaginary_leakage / solution_scale
    )
    if solution_relative_imaginary_leakage > 512.0 * np.finfo(np.float64).eps:
        raise FloatingPointError("unexpected imaginary leakage in real Poisson solve")
    psi = np.asarray(psi_complex.real, dtype=np.float64)
    # Assign the requested boundary trace exactly in physical space.
    psi[-1] = trace
    if not np.all(np.isfinite(psi)):
        raise FloatingPointError("Poisson solve produced a non-finite field")

    discrete_residual, residual_imaginary_leakage = (
        _spectral_discrete_residual(
            grid,
            psi1=psi,
            omega1=omega,
            outer_trace=trace,
            lower_flux=lower_flux,
            upper_flux=upper_flux,
        )
    )
    pde_residual = independent_physical_poisson_residual(
        grid,
        psi1=psi,
        omega1=omega,
        outer_dirichlet=trace,
    )
    interior_discrete = discrete_residual[:-1]
    interior_pde = pde_residual[:-1]
    residual_scale = max(
        1.0,
        float(np.max(np.abs(omega))),
        float(np.max(np.abs(psi)))
        * (grid.dr**-2 + grid.dz**-2),
    )
    residual_relative_imaginary_leakage = (
        residual_imaginary_leakage / residual_scale
    )
    zero_mode_condition = (
        _zero_mode_condition_number_inf(grid.nr)
        if estimate_condition
        else None
    )
    metadata: dict[str, Any] = {
        "equation": "-(d_rr + 3/r d_r + d_zz) psi1 = omega1",
        "radial_discretization": "second-order r^3-flux finite volume",
        "rhs_interpretation": (
            "nodal omega1 values approximate radial control-volume averages; "
            "for smooth even fields this replacement is second-order "
            "consistent, not exact cell integration"
        ),
        "axial_discretization": "NumPy FFT Fourier coefficients",
        "linear_solver": "per-mode Thomas tridiagonal elimination",
        "nr": grid.nr,
        "nz": grid.nz,
        "dr": grid.dr,
        "dz": grid.dz,
        "z_period": grid.z_period,
        "axis_radial_coefficient": float(upper_flux[0]),
        "axis_radial_coefficient_dimensionless": float(
            upper_flux[0] * grid.dr**2
        ),
        "outer_boundary_max_abs_defect": float(
            np.max(np.abs(psi[-1] - trace))
        ),
        "discrete_residual_rms_interior": float(
            math.sqrt(np.mean(interior_discrete**2))
        ),
        "discrete_residual_max_abs_interior": float(
            np.max(np.abs(interior_discrete))
        ),
        "independent_pde_residual_rms_interior": float(
            math.sqrt(np.mean(interior_pde**2))
        ),
        "independent_pde_residual_max_abs_interior": float(
            np.max(np.abs(interior_pde))
        ),
        "zero_mode_condition_number_inf": zero_mode_condition,
        "condition_estimate_interpretation": (
            "Unweighted coordinate-basis induced infinity-norm condition "
            "number of the zero-Fourier-mode radial matrix only; not a "
            "weighted coercivity bound"
            if estimate_condition
            else "not evaluated for this call"
        ),
        "solution_imaginary_leakage_max_abs": solution_imaginary_leakage,
        "residual_imaginary_leakage_max_abs": residual_imaginary_leakage,
        "solution_relative_imaginary_leakage": (
            solution_relative_imaginary_leakage
        ),
        "residual_relative_imaginary_leakage": (
            residual_relative_imaginary_leakage
        ),
        "interpretation": (
            "finite-dimensional elliptic solve diagnostics; not a continuum "
            "error bound and not singularity evidence"
        ),
    }
    return PoissonSolution(
        psi1=psi,
        discrete_residual=discrete_residual,
        pde_residual=pde_residual,
        metadata=metadata,
    )
