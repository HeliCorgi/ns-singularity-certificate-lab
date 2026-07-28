r"""Exact modal transparent (Dirichlet-to-Neumann) outer condition for solver A.

Scope and warning
-----------------
This module implements the **whole-space** outer condition designed in
``docs/whole_space_transition.md`` (equation W-1).  It is *not* the Hou physical
wall.  Equation-audit entry **E-27** fixes ``psi_1(t, 1, z) = 0`` as a
*physical* no-flow boundary for the Hou reproduction; that is a statement about
a real cylinder wall, and replacing it by the condition below would change the
problem being reproduced.  The two are therefore never selected implicitly:
:func:`solve_streamfunction_poisson_outer` defaults to
``boundary_condition="dirichlet"`` and the transparent path must be requested by
name, and every experiment summary records which one was used.

What the condition is
---------------------
Per axial Fourier mode ``k`` the radial operator is
``L_{5,k} = d_rr + (3/r) d_r - k^2`` and E-33(a)/(b) give the decaying
homogeneous solution ``psi = K_1(kr)/r``.  Outside the support of ``omega_1``
the whole-space solution *is* a multiple of that branch, so

.. math::

   \frac{d}{dr}\Big[\frac{K_1(kr)}{r}\Big]
   =-\frac{k K_0(kr)}{r}-\frac{2K_1(kr)}{r^2}

(using ``K_1'(x) = -K_0(x) - K_1(x)/x``, DLMF 10.29.2) and hence, at any radius
``R`` beyond the support,

.. math::

   \partial_r\hat\psi_k(R)+\beta_k\,\hat\psi_k(R)=0,
   \qquad
   \beta_k=\frac2R+k\,\frac{K_0(kR)}{K_1(kR)} .
   \tag{W-1}

Two limits fix the interpretation and are checked by the tests.

* ``k -> 0``: ``K_0(x)/K_1(x) -> 0`` like ``x(-log(x/2)-gamma)``, so
  ``k K_0(kR)/K_1(kR) -> 0`` and ``beta_0 = 2/R`` **exactly**.  The ``k = 0``
  decaying solution is ``psi = C/r^2``, for which
  ``psi'(R) + (2/R) psi(R) = -2C/R^3 + 2C/R^3 = 0``: the condition is satisfied
  identically, not approximately.  :func:`outer_bracket` returns ``2/R`` for
  ``k = 0`` from an explicit branch, never from a limit of the quadrature.
* ``kR -> infinity``: ``K_0/K_1 -> 1``, so ``beta_k -> 2/R + k`` and the
  condition degenerates to outgoing exponential decay ``psi' ~ -k psi``.

Under the premise that ``omega_1`` vanishes for ``r >= R`` the truncation error
of (W-1) is **exactly zero**; only the discretization error of the boundary row
survives.  That premise is not optional and is enforced by
:func:`assert_compact_support`, which refuses -- it does not warn -- when the
source reaches the wall.

Discretization of the boundary row, and why it is second order
--------------------------------------------------------------
The interior rows and the axis row are E-26's, unchanged, and are taken from
solver A itself (:func:`radial_flux_coefficients` extends
``poisson._radial_flux_coefficients`` by one row and is pinned bitwise against
it by ``tests/test_transparent_boundary.py``).  With ``r_i = i dr``,
``r_N = R``, ``V_i = (r_{i+1/2}^4 - r_{i-1/2}^4)/4`` and
``a_i^{\pm} = r_{i\pm1/2}^3/(dr V_i)``, those rows are

.. math::

   -a_i^-\hat\psi_{i-1}+(a_i^-+a_i^++k^2)\hat\psi_i-a_i^+\hat\psi_{i+1}
   =\hat\omega_i ,
   \qquad
   (a_0^-,a_0^+)=(0,\,8/dr^2).

**Ghost cell.**  Only the outermost row changes.  Because ``omega_1`` vanishes
for ``r >= R``, the continuum solution continues to exist and to be smooth
beyond ``R`` -- it is ``A K_1(kr)/r`` there -- so the node ``r_{N+1} = R + dr``
is not fictitious, and the *same* finite-volume row may be written at ``i = N``
over the full control volume ``[R - dr/2, R + dr/2]``:

.. math::

   -a_N^-\hat\psi_{N-1}+(a_N^-+a_N^++k^2)\hat\psi_N-a_N^+\hat\psi_{N+1}
   =\hat\omega_N .

The ghost value is eliminated with the *centered* second-order form of (W-1) at
``r = R``,

.. math::

   \frac{\hat\psi_{N+1}-\hat\psi_{N-1}}{2\,dr}+\beta_k\hat\psi_N=0
   \quad\Longrightarrow\quad
   \hat\psi_{N+1}=\hat\psi_{N-1}-2\,dr\,\beta_k\hat\psi_N ,

which gives the single replaced row actually assembled here:

.. math::

   -(a_N^-+a_N^+)\hat\psi_{N-1}
   +\big(a_N^-+a_N^++k^2+2\,dr\,\beta_k a_N^+\big)\hat\psi_N
   =\hat\omega_N .
   \tag{W-1h}

**Order.**  Two error sources enter (W-1h).  The finite-volume row itself is
E-26's and is second-order consistent for smooth even fields.  The ghost
elimination is exact up to

.. math::

   \hat\psi(R+dr)-\big[\hat\psi(R-dr)-2\,dr\,\beta_k\hat\psi(R)\big]
   =\tfrac{dr^3}{3}\hat\psi'''(R)+O(dr^5),

because ``beta_k psi(R) = -psi'(R)`` holds *exactly* for the continuum
solution.  Multiplying by ``a_N^+ = O(dr^{-2})`` leaves a local truncation error
``O(dr)`` in that one row -- the boundary row is genuinely one order worse
locally, and saying otherwise would be false.  It is nevertheless **second
order globally**: the assembled matrix is an irreducibly diagonally dominant
M-matrix (row ``N`` is strictly dominant because ``beta_k > 0``), its rows carry
the scaling ``1/V_i = O(dr^{-1} r^{-3})`` of the flux form, and consequently the
column of ``A^{-1}`` belonging to the boundary row is ``O(dr)`` in the max norm
-- the discrete Green's function of a flux-form operator picks up one factor of
the cell width at the boundary.  An ``O(dr)`` defect confined to that row
therefore contributes ``O(dr^2)`` to the solution, matching the interior order.
This is the standard ghost-cell Neumann/Robin argument, and it is *measured*
rather than assumed: ``experiments/run_transparent_boundary.py`` reports the
observed order of the max-norm error over the whole cylinder **including the
boundary node**, against a manufactured solution that satisfies (W-1) exactly.

Half-cell alternative (rejected).  Integrating ``r^3`` times the equation over
``[R - dr/2, R]`` and substituting ``psi'(R) = -beta_k psi(R)`` also gives a
tridiagonal row, but its cell-average approximations are centred on the
half-cell while the unknown sits at its edge, which costs an extra ``O(dr)``
*and* makes the row's volume differ from the interior convention.  The ghost
form above keeps ``V_N`` in the interior family, so ``a_N^{\pm}`` come from the
very same expression as every other row.

Fault-injection variants
------------------------
:func:`outer_bracket` accepts a ``variant`` from :data:`BRACKET_VARIANTS`.
Anything other than ``"exact"`` is a deliberately wrong condition used by the
acceptance experiment to show that the checks detect it.  The variant is a
required part of the solver metadata, so a faulted run can never be mistaken
for a good one in a summary.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import numpy.typing as npt

from . import bessel_reference
from .grid import AxisymmetricGrid, FloatArray
from .poisson import (
    PoissonSolution,
    _radial_flux_coefficients,
    _solve_tridiagonal,
    independent_physical_poisson_residual,
    solve_streamfunction_poisson,
)

__all__ = [
    "BOUNDARY_CONDITIONS",
    "BRACKET_VARIANTS",
    "DIRICHLET",
    "TRANSPARENT",
    "assemble_radial_mode",
    "assert_compact_support",
    "modal_wavenumbers",
    "outer_bracket",
    "radial_flux_coefficients",
    "solve_radial_mode",
    "solve_streamfunction_poisson_outer",
    "transparent_condition_defect",
]

DIRICHLET = "dirichlet"
TRANSPARENT = "transparent"

BOUNDARY_CONDITIONS: tuple[str, ...] = (DIRICHLET, TRANSPARENT)
"""The outer conditions this module can impose.  Dirichlet is the default."""

BRACKET_VARIANTS: tuple[str, ...] = (
    "exact",
    "sign_flipped",
    "no_curvature_term",
    "frozen_ratio",
)
"""``"exact"`` is (W-1); the rest are the named faults of W-B condition 5."""

ComplexArray = npt.NDArray[np.complex128]

_EPS = float(np.finfo(np.float64).eps)


def modal_wavenumbers(grid: AxisymmetricGrid) -> FloatArray:
    """Return the signed axial wavenumbers in NumPy FFT storage order.

    The same expression solver A uses, kept here so that the boundary row and
    the interior rows cannot drift apart in their wavenumber convention.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    if not grid.periodic_z:
        raise ValueError("the transparent outer condition requires periodic z")
    return np.asarray(
        2.0 * np.pi * np.fft.fftfreq(grid.nz, d=grid.dz), dtype=np.float64
    )


def outer_bracket(
    wavenumber: float, radius: float, *, variant: str = "exact"
) -> float:
    r"""Return the W-1 bracket ``2/R + |k| K_0(|k|R)/K_1(|k|R)``.

    ``wavenumber`` may be negative: the radial problem depends on ``k`` only
    through ``k^2`` and through the decay rate ``|k|``, and NumPy's FFT stores
    negative frequencies in the upper half of the spectrum, so the absolute
    value is taken.  Using the signed value there would impose an *inward*
    growing condition on half of the modes.

    ``k = 0`` is handled by an explicit branch that returns ``2/R`` exactly;
    ``K_0/K_1`` is never evaluated at zero argument.  Outside the range where
    the oracle's quadrature is the accurate branch, the corresponding
    *documented* asymptote is used instead, and each substitution is counted in
    the metadata of :func:`solve_streamfunction_poisson_outer` rather than made
    silently:

    * below :data:`~ns_certificate_lab.bessel_reference.MINIMUM_K_ARGUMENT` the
      quadrature refuses to run at all, and
      ``x(-log(x/2) - gamma)`` is used; there the whole term is below
      ``1e-9 |k|`` and cannot matter next to ``2/R``;
    * above
      :data:`~ns_certificate_lab.bessel_reference.LARGE_ARGUMENT_QUADRATURE_LIMIT`
      the quadrature has too few nodes left and ``1 - 1/(2x) + 3/(8x^2)`` is
      the more accurate of the two (see that constant's documentation).  Modes
      that far into ``kR >> 1`` are damped by ``e^{-2kR}`` and contribute
      nothing to the core, but leaving a knowingly degraded formula in place
      would be a trap for a later, finer grid.
    """

    if variant not in BRACKET_VARIANTS:
        raise ValueError(
            f"variant must be one of {BRACKET_VARIANTS!r}, got {variant!r}"
        )
    if isinstance(wavenumber, bool) or not isinstance(wavenumber, (int, float)):
        raise TypeError("wavenumber must be a real number")
    if isinstance(radius, bool) or not isinstance(radius, (int, float)):
        raise TypeError("radius must be a real number")
    magnitude = abs(float(wavenumber))
    outer_radius = float(radius)
    if not math.isfinite(magnitude):
        raise ValueError("wavenumber must be finite")
    if not math.isfinite(outer_radius) or outer_radius <= 0.0:
        raise ValueError("radius must be positive and finite")

    curvature = 2.0 / outer_radius
    if magnitude == 0.0:
        decay = 0.0
    elif variant == "frozen_ratio":
        # The kR >> 1 asymptote K_0/K_1 -> 1 applied everywhere: the misuse
        # W-B condition 5 asks to detect.
        decay = magnitude
    else:
        argument = magnitude * outer_radius
        if argument < bessel_reference.MINIMUM_K_ARGUMENT:
            ratio = bessel_reference.k0_over_k1_small_argument_asymptote(argument)
        elif argument > bessel_reference.LARGE_ARGUMENT_QUADRATURE_LIMIT:
            ratio = bessel_reference.k0_over_k1_large_argument_asymptote(argument)
        else:
            ratio = bessel_reference.k0_over_k1(argument)
        decay = magnitude * ratio

    if variant == "no_curvature_term":
        return decay
    bracket = curvature + decay
    if variant == "sign_flipped":
        return -bracket
    return bracket


def radial_flux_coefficients(
    grid: AxisymmetricGrid,
) -> tuple[FloatArray, FloatArray]:
    """Return E-26's ``a^-``/``a^+`` for **all** ``nr`` rows, wall included.

    Solver A's own helper stops one row short because its outermost radial
    point is prescribed Dirichlet data.  The transparent path makes that point
    an unknown, so one further row of the identical expression is needed.  The
    first ``nr - 1`` entries are therefore *taken from solver A itself* rather
    than recomputed, which is what makes "only the final row changes" a
    structural fact instead of a claim; a test pins the appended row against an
    independent evaluation of the same E-26 formula.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    lower_interior, upper_interior = _radial_flux_coefficients(grid)
    outer_radius = float(grid.r[-1])
    half_lower = outer_radius - 0.5 * grid.dr
    half_upper = outer_radius + 0.5 * grid.dr
    volume = 0.25 * (half_upper**4 - half_lower**4)
    lower = np.concatenate(
        (lower_interior, np.array([half_lower**3 / (grid.dr * volume)]))
    )
    upper = np.concatenate(
        (upper_interior, np.array([half_upper**3 / (grid.dr * volume)]))
    )
    return lower, upper


def assemble_radial_mode(
    grid: AxisymmetricGrid,
    wavenumber: float,
    *,
    bracket_variant: str = "exact",
    flux: tuple[FloatArray, FloatArray] | None = None,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return ``(sub, diagonal, super)`` for one mode's ``nr``-row system.

    Rows ``0 .. nr-2`` are E-26's, byte for byte; row ``nr-1`` is (W-1h).  The
    matrix is an irreducibly diagonally dominant M-matrix whenever the bracket
    is positive, with strict dominance in the final row -- that is what makes
    the unpivoted Thomas elimination of solver A safe here.  A fault variant
    that makes the bracket zero (``"no_curvature_term"`` at ``k = 0``) destroys
    exactly that property and the elimination refuses; see the experiment.
    """

    lower_flux, upper_flux = (
        radial_flux_coefficients(grid) if flux is None else flux
    )
    radius = float(grid.r[-1])
    bracket = outer_bracket(wavenumber, radius, variant=bracket_variant)
    diagonal = lower_flux + upper_flux + float(wavenumber) ** 2
    diagonal[-1] += 2.0 * grid.dr * bracket * upper_flux[-1]
    sub_diagonal = -lower_flux[1:].copy()
    sub_diagonal[-1] = -(lower_flux[-1] + upper_flux[-1])
    super_diagonal = -upper_flux[:-1].copy()
    return sub_diagonal, diagonal, super_diagonal


def solve_radial_mode(
    grid: AxisymmetricGrid,
    modal_source: npt.ArrayLike,
    wavenumber: float,
    *,
    bracket_variant: str = "exact",
    flux: tuple[FloatArray, FloatArray] | None = None,
) -> ComplexArray:
    """Solve one axial mode's radial system under the transparent condition."""

    source = np.asarray(modal_source, dtype=np.complex128)
    if source.shape != (grid.nr,):
        raise ValueError(f"modal_source must have shape ({grid.nr},)")
    sub_diagonal, diagonal, super_diagonal = assemble_radial_mode(
        grid, wavenumber, bracket_variant=bracket_variant, flux=flux
    )
    return _solve_tridiagonal(
        sub_diagonal, diagonal, super_diagonal, source.copy()
    )


def assert_compact_support(
    grid: AxisymmetricGrid,
    omega1: FloatArray,
    *,
    support_radius: float,
) -> int:
    """Refuse a source that is not compactly supported strictly inside ``R``.

    (W-1) is an *exact* condition only when the solution beyond ``R`` is the
    pure decaying branch, i.e. only when ``omega_1`` vanishes there.  A source
    that reaches the wall would still produce numbers, and they would silently
    mean something else, so this is an error rather than a warning -- the same
    stance ``experiments/run_wall_truncation_scaling.py`` takes for E-33.

    Returns the number of interior grid rows on which the source is exactly
    zero, which must be at least one: otherwise compact support is not
    observable on this grid at all.
    """

    outer_radius = float(grid.r[-1])
    if isinstance(support_radius, bool) or not isinstance(
        support_radius, (int, float)
    ):
        raise TypeError("support_radius must be a real number")
    radius = float(support_radius)
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError("support_radius must be positive and finite")
    if radius >= outer_radius:
        raise ValueError(
            "the transparent outer condition requires omega1 to be compactly "
            f"supported strictly inside the wall: support radius {radius!r} is "
            f"not smaller than wall radius {outer_radius!r}"
        )
    outside = grid.r >= radius
    interior_zero_rows = int(np.count_nonzero(outside[:-1]))
    if interior_zero_rows < 1:
        raise ValueError(
            "no interior grid node lies between the source support and the "
            "wall, so compact support is not observable on this grid"
        )
    residue = float(np.max(np.abs(omega1[outside, :])))
    if residue != 0.0:
        raise ValueError(
            "omega1 is not exactly zero outside its declared support radius "
            f"{radius!r}: max |omega1| there is {residue!r}; the transparent "
            "outer condition is not valid for this source"
        )
    return interior_zero_rows


def transparent_condition_defect(
    grid: AxisymmetricGrid,
    psi1: FloatArray,
    *,
    variant: str = "exact",
) -> FloatArray:
    r"""Return an *independent* second-order defect of (W-1) at ``r = R``.

    The solve eliminates a ghost node with the centered difference
    ``(psi_{N+1} - psi_{N-1})/(2 dr)``.  This diagnostic instead uses the
    one-sided second-order derivative

    .. math::

       \hat\psi'(R)\approx
       \frac{3\hat\psi_N-4\hat\psi_{N-1}+\hat\psi_{N-2}}{2\,dr},

    and reports ``psi'(R) + beta_k psi(R)`` per mode, transformed back to
    physical space.  It is a different discretization of the same condition, so
    it is a consistency check and not an algebraic identity: it must be
    ``O(dr^2)``, not zero.  If it were zero the diagnostic would merely be
    restating the solve.
    """

    wave_numbers = modal_wavenumbers(grid)
    radius = float(grid.r[-1])
    psi_hat = np.fft.fft(np.asarray(psi1, dtype=np.float64), axis=1)
    one_sided = (
        3.0 * psi_hat[-1] - 4.0 * psi_hat[-2] + psi_hat[-3]
    ) / (2.0 * grid.dr)
    brackets = np.array(
        [
            outer_bracket(float(value), radius, variant=variant)
            for value in wave_numbers
        ],
        dtype=np.float64,
    )
    defect_hat = one_sided + brackets * psi_hat[-1]
    return np.asarray(np.fft.ifft(defect_hat).real, dtype=np.float64)


def _modal_discrete_residual(
    grid: AxisymmetricGrid,
    *,
    psi1: FloatArray,
    omega1: FloatArray,
    lower_flux: FloatArray,
    upper_flux: FloatArray,
    brackets: FloatArray,
) -> tuple[FloatArray, float]:
    """Apply the exact discrete operator the transparent path assembles."""

    psi_hat = np.fft.fft(psi1, axis=1)
    omega_hat = np.fft.fft(omega1, axis=1)
    wave_numbers = modal_wavenumbers(grid)
    residual_hat = np.zeros_like(psi_hat)
    for mode, wave_number in enumerate(wave_numbers):
        values = psi_hat[:, mode]
        applied = (lower_flux + upper_flux + wave_number**2) * values
        applied[1:] -= lower_flux[1:] * values[:-1]
        applied[:-1] -= upper_flux[:-1] * values[1:]
        # (W-1h): the ghost value psi_{N+1} = psi_{N-1} - 2 dr beta psi_N.
        applied[-1] -= upper_flux[-1] * (
            values[-2] - 2.0 * grid.dr * brackets[mode] * values[-1]
        )
        residual_hat[:, mode] = applied - omega_hat[:, mode]
    residual_complex = np.fft.ifft(residual_hat, axis=1)
    leakage = float(np.max(np.abs(residual_complex.imag)))
    return np.asarray(residual_complex.real, dtype=np.float64), leakage


def _reject_real_array(values: npt.ArrayLike, *, name: str) -> None:
    raw = np.asarray(values)
    if np.iscomplexobj(raw) or (
        raw.dtype.kind == "O"
        and any(isinstance(value, complex) for value in raw.flat)
    ):
        raise ValueError(f"{name} must be real-valued")


def solve_streamfunction_poisson_outer(
    grid: AxisymmetricGrid,
    omega1: npt.ArrayLike,
    *,
    boundary_condition: str = DIRICHLET,
    outer_dirichlet: npt.ArrayLike = 0.0,
    support_radius: float | None = None,
    bracket_variant: str = "exact",
    estimate_condition: bool = True,
) -> PoissonSolution:
    """Solve ``-L5 psi1 = omega1`` with a selectable outer condition.

    ``boundary_condition="dirichlet"`` (the default, and the only behaviour any
    existing caller can get) delegates verbatim to
    :func:`ns_certificate_lab.poisson.solve_streamfunction_poisson`, so that
    path is bit-identical by construction and pinned by
    ``tests/test_transparent_boundary.py``.

    ``boundary_condition="transparent"`` imposes (W-1).  It requires
    ``support_radius``: the condition is exact only when ``omega1`` vanishes
    outside that radius, and the premise is checked, not assumed.  There is no
    outer trace in that case -- the condition is homogeneous -- so passing a
    nonzero ``outer_dirichlet`` is rejected instead of ignored.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    if not isinstance(estimate_condition, bool):
        raise TypeError("estimate_condition must be boolean")
    if boundary_condition not in BOUNDARY_CONDITIONS:
        raise ValueError(
            f"boundary_condition must be one of {BOUNDARY_CONDITIONS!r}, got "
            f"{boundary_condition!r}"
        )
    if bracket_variant not in BRACKET_VARIANTS:
        raise ValueError(
            f"bracket_variant must be one of {BRACKET_VARIANTS!r}, got "
            f"{bracket_variant!r}"
        )

    if boundary_condition == DIRICHLET:
        if support_radius is not None:
            raise ValueError(
                "support_radius is meaningful only for the transparent outer "
                "condition; the Dirichlet path does not assume compact support"
            )
        if bracket_variant != "exact":
            raise ValueError(
                "bracket_variant is meaningful only for the transparent outer "
                "condition"
            )
        return solve_streamfunction_poisson(
            grid,
            omega1,
            outer_dirichlet,
            estimate_condition=estimate_condition,
        )

    if not grid.periodic_z:
        raise ValueError("the transparent outer condition requires periodic_z=True")
    trace = np.asarray(outer_dirichlet, dtype=np.float64)
    if trace.ndim != 0 or float(trace) != 0.0:
        raise ValueError(
            "the transparent outer condition is homogeneous; outer_dirichlet "
            "must be left at its default 0.0"
        )
    if support_radius is None:
        raise ValueError(
            "the transparent outer condition requires an explicit "
            "support_radius: it is exact only when omega1 vanishes for r >= R, "
            "and that premise is verified rather than assumed"
        )

    _reject_real_array(omega1, name="omega1")
    omega = grid.validate_field(omega1, name="omega1")
    interior_zero_rows = assert_compact_support(
        grid, omega, support_radius=support_radius
    )

    radius = float(grid.r[-1])
    lower_flux, upper_flux = radial_flux_coefficients(grid)
    wave_numbers = modal_wavenumbers(grid)
    brackets = np.array(
        [
            outer_bracket(float(value), radius, variant=bracket_variant)
            for value in wave_numbers
        ],
        dtype=np.float64,
    )
    arguments = np.abs(wave_numbers) * radius
    small_argument_modes = int(
        np.count_nonzero(
            (np.abs(wave_numbers) > 0.0)
            & (arguments < bessel_reference.MINIMUM_K_ARGUMENT)
        )
    )
    large_argument_modes = int(
        np.count_nonzero(
            arguments > bessel_reference.LARGE_ARGUMENT_QUADRATURE_LIMIT
        )
    )

    omega_hat = np.fft.fft(omega, axis=1)
    psi_hat = np.empty(grid.shape, dtype=np.complex128)
    for mode, wave_number in enumerate(wave_numbers):
        psi_hat[:, mode] = solve_radial_mode(
            grid,
            omega_hat[:, mode],
            float(wave_number),
            bracket_variant=bracket_variant,
            flux=(lower_flux, upper_flux),
        )

    psi_complex = np.fft.ifft(psi_hat, axis=1)
    solution_imaginary_leakage = float(np.max(np.abs(psi_complex.imag)))
    solution_scale = max(1.0, float(np.max(np.abs(psi_complex.real))))
    solution_relative_imaginary_leakage = solution_imaginary_leakage / solution_scale
    if solution_relative_imaginary_leakage > 512.0 * _EPS:
        raise FloatingPointError(
            "unexpected imaginary leakage in real transparent Poisson solve"
        )
    psi = np.asarray(psi_complex.real, dtype=np.float64)
    if not np.all(np.isfinite(psi)):
        raise FloatingPointError(
            "transparent Poisson solve produced a non-finite field"
        )

    discrete_residual, residual_imaginary_leakage = _modal_discrete_residual(
        grid,
        psi1=psi,
        omega1=omega,
        lower_flux=lower_flux,
        upper_flux=upper_flux,
        brackets=brackets,
    )
    # The interior rows are checked by solver A's independent physical-space
    # path.  Its final row is a Dirichlet-trace defect, which is meaningless
    # here, so it is handed the solved trace (making that row identically zero)
    # and then overwritten with the modal condition defect.  The transparent
    # condition is diagonal in z, not local in z, so there is no physical-space
    # stencil for it -- that is exactly the incompatibility recorded for
    # solver C in docs/whole_space_transition.md section 3, note 2.
    pde_residual = independent_physical_poisson_residual(
        grid,
        psi1=psi,
        omega1=omega,
        outer_dirichlet=psi[-1],
    )
    boundary_defect = transparent_condition_defect(
        grid, psi, variant=bracket_variant
    )
    pde_residual = pde_residual.copy()
    pde_residual[-1] = boundary_defect

    interior_discrete = discrete_residual[:-1]
    interior_pde = pde_residual[:-1]
    residual_scale = max(
        1.0,
        float(np.max(np.abs(omega))),
        float(np.max(np.abs(psi))) * (grid.dr**-2 + grid.dz**-2),
    )
    solution_max = float(np.max(np.abs(psi)))
    defect_scale = max(
        float(np.finfo(np.float64).tiny), solution_max / max(grid.dr, _EPS)
    )

    metadata: dict[str, Any] = {
        "equation": "-(d_rr + 3/r d_r + d_zz) psi1 = omega1",
        "outer_boundary_condition": TRANSPARENT,
        "outer_boundary_reference": (
            "docs/whole_space_transition.md (W-1); modal Dirichlet-to-Neumann "
            "condition d_r psi_hat_k(R) + [2/R + |k| K_0(|k|R)/K_1(|k|R)] "
            "psi_hat_k(R) = 0"
        ),
        "outer_boundary_scope_warning": (
            "this is the WHOLE-SPACE outer condition; it is NOT the E-27 "
            "physical no-slip wall used by the Hou reproduction, and the two "
            "must not be interchanged"
        ),
        "outer_boundary_row_discretization": (
            "ghost cell at r=R+dr eliminated by the centered second-order form "
            "of W-1; interior r^3-flux rows and the axis row (coefficient 8) "
            "are unchanged"
        ),
        "bracket_variant": bracket_variant,
        "bracket_is_exact": bracket_variant == "exact",
        "wall_radius": radius,
        "support_radius": float(support_radius),
        "interior_zero_source_rows": interior_zero_rows,
        "zero_mode_bracket": float(brackets[0]),
        "zero_mode_bracket_closed_form": 2.0 / radius,
        "max_abs_bracket": float(np.max(np.abs(brackets))),
        "min_abs_bracket": float(np.min(np.abs(brackets))),
        "small_argument_asymptote_modes": small_argument_modes,
        "large_argument_asymptote_modes": large_argument_modes,
        "max_modal_bessel_argument": float(np.max(arguments)),
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
        "axis_radial_coefficient_dimensionless": float(upper_flux[0] * grid.dr**2),
        "outer_condition_defect_max_abs": float(np.max(np.abs(boundary_defect))),
        "outer_condition_defect_relative": float(
            np.max(np.abs(boundary_defect)) / defect_scale
        ),
        "discrete_residual_rms_interior": float(
            math.sqrt(np.mean(interior_discrete**2))
        ),
        "discrete_residual_max_abs_interior": float(np.max(np.abs(interior_discrete))),
        "discrete_residual_max_abs_boundary_row": float(
            np.max(np.abs(discrete_residual[-1]))
        ),
        "independent_pde_residual_rms_interior": float(
            math.sqrt(np.mean(interior_pde**2))
        ),
        "independent_pde_residual_max_abs_interior": float(
            np.max(np.abs(interior_pde))
        ),
        "zero_mode_condition_number_inf": None,
        "condition_estimate_interpretation": (
            "solver A's closed-form zero-mode condition number assumes the "
            "Dirichlet outer row and does not describe the transparent matrix; "
            "it is deliberately not reported here"
        ),
        "solution_imaginary_leakage_max_abs": solution_imaginary_leakage,
        "residual_imaginary_leakage_max_abs": residual_imaginary_leakage,
        "solution_relative_imaginary_leakage": solution_relative_imaginary_leakage,
        "residual_relative_imaginary_leakage": (
            residual_imaginary_leakage / residual_scale
        ),
        "interpretation": (
            "finite-dimensional elliptic solve diagnostics; the transparent "
            "condition removes the wall-truncation error under its compact "
            "support premise but leaves the O(dr^2) boundary discretization "
            "error, and none of this is a continuum error bound or singularity "
            "evidence"
        ),
    }
    return PoissonSolution(
        psi1=psi,
        discrete_residual=discrete_residual,
        pde_residual=pde_residual,
        metadata=metadata,
    )
