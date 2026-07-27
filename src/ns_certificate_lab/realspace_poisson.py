r"""Third, real-space finite-cylinder solver for ``-L5 psi1 = omega1``.

Purpose
-------
Two solvers for this boundary-value problem already exist in the package:
``ns_certificate_lab.poisson`` (solver A, ``r^3``-flux finite volume) and
``ns_certificate_lab.finite_cylinder_poisson`` (solver B, direct
non-divergence stencil).  Their radial stencils differ, but their *entire*
axial treatment is identical: both diagonalize the periodic direction with
NumPy's discrete Fourier transform, both take their axial wavenumbers from the
same NumPy sample-frequency helper, and both invert the resulting radial
systems with a locally written Thomas elimination.
A defect in either of those two shared components would move both solvers the
same way and would therefore be invisible to their cross-validation.

This module is solver C.  It shares neither of those components:

* the periodic axial direction is discretized with the ordinary three-point
  centered difference, applied with :func:`numpy.roll` -- no transform of any
  kind is used anywhere in this file;
* the linear system is inverted by a matrix-free Jacobi-preconditioned
  conjugate-gradient iteration -- no tridiagonal elimination is performed;
* the radial coefficients are transcribed here directly from the equation
  audit (E-26a/E-26b), not imported or copied from either sibling solver.

What remains shared is stated plainly:
:class:`~ns_certificate_lab.grid.AxisymmetricGrid` is imported as the common
data structure that fixes the ``(n_r, n_z)`` layout, the uniform spacings and
the field-shape validation.  That dependency is unavoidable -- it is what makes
the three solvers comparable on *the same* grid -- and it is not an
independent check of itself.  All three solvers also share binary64 arithmetic
and the E-26 discretization mathematics.

Equation and boundary conditions (E-13, E-16, E-25)
---------------------------------------------------
.. math::

   -\mathcal L_5\psi_1
   = -\left(\partial_{rr}+\frac{3}{r}\partial_r+\partial_{zz}\right)\psi_1
   = \omega_1,
   \qquad 0\le r\le R,\quad z\in[0,L_z)\ \text{periodic},

with the outer Dirichlet trace ``psi1(R, z) = g(z)`` and the axis condition
``d_r psi1(0, z) = 0`` implied by evenness (E-16c).  For an even smooth field
the axis limit is ``(L5 f)(0,z) = 4 f_rr(0,z) + f_zz(0,z)`` (E-17); the radial
part of that limit is what the axis control volume reproduces below.

Radial discretization (E-26a, E-26b), transcribed
-------------------------------------------------
With ``r_i = i*dr``, ``r_{i+1/2} = (i+1/2)*dr`` and the control volumes

.. math::

   V_i=\int_{r_{i-1/2}}^{r_{i+1/2}}r^3\,dr
       =\frac{r_{i+1/2}^4-r_{i-1/2}^4}{4},

the interior row of the radial operator ``r^{-3} d_r (r^3 d_r psi)`` is

.. math::

   (\mathcal L_{5,r}^{h}\psi)_i=
   \frac{r_{i+1/2}^3(\psi_{i+1}-\psi_i)/\Delta r
         -r_{i-1/2}^3(\psi_i-\psi_{i-1})/\Delta r}{V_i}.

The axis cell has ``r_{-1/2} = 0`` and ``r_{1/2} = dr/2``, hence
``V_0 = (dr/2)^4/4 = dr^4/64`` and

.. math::

   (\mathcal L_{5,r}^{h}\psi)_0=\frac{8(\psi_1-\psi_0)}{\Delta r^2},

which is E-26b.  This module does *not* hard-code the value ``8``: it evaluates
the general flux formula and then checks that the resulting axis coefficient
agrees with ``8/dr^2`` to a relative tolerance of ``1e-12``, raising if it does
not.  The coefficient is therefore pinned by the audit formula rather than
asserted by fiat.

Axial discretization
--------------------
``(d_zz psi)_{ij} = (psi_{i,j+1} - 2 psi_{ij} + psi_{i,j-1}) / dz^2`` with the
index wrapped periodically.  This is the component that breaks the common mode
shared by solvers A and B.  Consequently solver C is expected to differ from
solver A on the *same* grid by an ``O(dz^2)`` amount even though their radial
stencils agree: solver A applies the exact symbol ``k^2`` to each axial mode
while this module applies ``(4/dz^2) sin^2(k dz/2) = k^2 - k^4 dz^2/12 + ...``.
That gap is a feature, not a defect; it is what makes the comparison
non-vacuous, and ``tests/test_realspace_poisson.py`` bounds it from both sides.

Symmetry of the reduced system
------------------------------
Let the unknowns be rows ``i = 0 .. n_r-2``; the outer row is prescribed by the
Dirichlet trace and is eliminated into the right-hand side.  Define the weighted
inner product

.. math::

   \langle u,v\rangle_w=\sum_{i,j} w_i\,u_{ij}v_{ij},
   \qquad w_i = V_i\,\Delta z,

which is the Riemann sum of ``\int\int u v\, r^3\,dr\,dz`` for the control
volumes above.  Write the face coefficient

.. math::

   F_i=\frac{r_{i+1/2}^3}{\Delta r},

so the audit's coefficients are ``a_i^+ = F_i / V_i`` and
``a_{i+1}^- = F_i / V_{i+1}``.  The key identity is that *both* products equal
the same face quantity,

.. math::

   a_i^+ V_i = F_i = a_{i+1}^- V_{i+1},

which is exactly the statement that the flux leaving cell ``i`` through the
face at ``r_{i+1/2}`` is the flux entering cell ``i+1`` through that same face.
Abel summation in ``i`` (using ``F_{-1} = 0`` at the axis and ``u = v = 0`` on
the eliminated outer row) and periodic summation by parts in ``j`` give, for
``A = -(\mathcal L_{5,r}^{h} + d_{zz}^{h})``,

.. math::

   \langle A u, v\rangle_w
   = \Delta z\sum_j\sum_{i=0}^{n_r-2}F_i(u_{i+1,j}-u_{ij})(v_{i+1,j}-v_{ij})
     + \frac{1}{\Delta z}\sum_i V_i\sum_j
       (u_{i,j+1}-u_{ij})(v_{i,j+1}-v_{ij}).

The right-hand side is manifestly symmetric in ``u`` and ``v`` and
non-negative for ``u = v``.  It vanishes only when ``u`` is constant in both
indices, and the radial term's outer boundary value is zero, so ``u = 0``.  The
reduced operator is therefore symmetric positive definite in
``\langle\cdot,\cdot\rangle_w``, which is what licenses conjugate gradients.
Note that the operator is *not* symmetric in the unweighted coordinate basis;
the weight is essential and is the same ``r^3 dr`` measure that E-26 warns must
not be confused with an unweighted condition number.

Because the implementation applies the operator through the single array ``F``
of face coefficients, the identity ``a_i^+ V_i = a_{i+1}^- V_{i+1}`` holds in
floating point as an equality of the *same* stored number, not merely to
rounding.  The symmetry test in the test module verifies this numerically
rather than trusting the argument.

Limitations
-----------
All arithmetic is binary64 and is not outward rounded.  The returned residual
is the algebraic residual of the discrete system defined above; it is not a
continuum error bound, and it says nothing about the artificial outer boundary
at ``r = R`` (E-25).  Nothing here is evidence about singularity formation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid, FloatArray

__all__ = [
    "ConvergenceError",
    "RealspaceOperator",
    "RealspacePoissonSolution",
    "build_realspace_operator",
    "solve_realspace_poisson",
]

_AXIS_COEFFICIENT_TOLERANCE = 1.0e-12


class ConvergenceError(RuntimeError):
    """Raised when the conjugate-gradient iteration exceeds its cap.

    The exception message carries the iteration count and the achieved
    weighted relative residual so that a failure is reported with numbers
    rather than as a bare "did not converge".
    """


def _reject_complex(value: npt.ArrayLike, *, name: str) -> None:
    """Reject complex input instead of silently discarding the imaginary part.

    ``numpy.asarray(..., dtype=float64)`` would drop the imaginary part behind
    a :class:`numpy.ComplexWarning`.  Object arrays holding Python ``complex``
    scalars are caught separately because their dtype kind is ``"O"``.
    """

    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real-valued")
    if raw.dtype.kind == "O" and any(
        isinstance(item, complex) for item in raw.flat
    ):
        raise ValueError(f"{name} must be real-valued")


def _outer_trace(
    grid: AxisymmetricGrid,
    outer_dirichlet: npt.ArrayLike,
) -> FloatArray:
    """Normalize the outer Dirichlet data to a finite length-``n_z`` array."""

    _reject_complex(outer_dirichlet, name="outer_dirichlet")
    trace = np.asarray(outer_dirichlet, dtype=np.float64)
    if trace.ndim == 0:
        trace = np.full(grid.nz, float(trace), dtype=np.float64)
    elif trace.shape == (grid.nz,):
        trace = np.array(trace, dtype=np.float64, copy=True)
    else:
        raise ValueError(
            f"outer_dirichlet must be scalar or have shape ({grid.nz},), "
            f"got shape {trace.shape}"
        )
    if not np.all(np.isfinite(trace)):
        raise ValueError("outer_dirichlet must contain only finite values")
    return trace


@dataclass(frozen=True)
class RealspaceOperator:
    """Matrix-free reduced operator ``A = -(L5_r^h + d_zz^h)`` on the unknowns.

    The unknown block has shape ``(n_r - 1, n_z)``: the outer radial row is
    prescribed by Dirichlet data and is not an unknown.  :meth:`apply` treats
    that eliminated row as zero, which is the homogeneous operator whose
    inverse conjugate gradients builds; the inhomogeneous Dirichlet data is
    moved to the right-hand side by :func:`solve_realspace_poisson`.
    """

    grid: AxisymmetricGrid
    face_coefficient: FloatArray
    """``F_i = r_{i+1/2}^3 / dr`` for the faces ``i = 0 .. n_r-2``."""

    volume: FloatArray
    """``V_i = (r_{i+1/2}^4 - r_{i-1/2}^4) / 4`` for ``i = 0 .. n_r-2``."""

    weight: FloatArray
    """``w_i = V_i * dz``, the discrete ``r^3 dr dz`` measure."""

    diagonal: FloatArray
    """Main diagonal of ``A``; constant along ``z`` and strictly positive."""

    @property
    def shape(self) -> tuple[int, int]:
        return (self.grid.nr - 1, self.grid.nz)

    def _validate(self, field: npt.ArrayLike, *, name: str) -> FloatArray:
        array = np.asarray(field, dtype=np.float64)
        if array.shape != self.shape:
            raise ValueError(
                f"{name} must have shape {self.shape}, got {array.shape}"
            )
        return array

    def apply(self, interior: npt.ArrayLike) -> FloatArray:
        """Return ``A u`` for the unknown block ``u``.

        Radial rows use the flux form of E-26a; the axis row degenerates to
        ``8 (u_1 - u_0) / dr^2`` automatically because ``r_{-1/2} = 0``.  The
        axial term is the periodic three-point centered difference.
        """

        field = self._validate(interior, name="interior")
        extended = np.zeros(self.grid.shape, dtype=np.float64)
        extended[:-1] = field
        # Face differences: entry i is (u_{i+1} - u_i) across the face at
        # r_{i+1/2}.  The last entry uses the eliminated (zero) outer row.
        face_difference = extended[1:] - extended[:-1]
        flux = self.face_coefficient[:, None] * face_difference
        divergence = np.empty_like(field)
        divergence[0] = flux[0]  # F_{-1} = 0 at the axis
        divergence[1:] = flux[1:] - flux[:-1]
        radial = divergence / self.volume[:, None]
        axial = (
            np.roll(field, -1, axis=1)
            - 2.0 * field
            + np.roll(field, 1, axis=1)
        ) / self.grid.dz**2
        return -(radial + axial)

    def inner(self, left: npt.ArrayLike, right: npt.ArrayLike) -> float:
        """Return the weighted inner product ``<left, right>_w``."""

        a = self._validate(left, name="left")
        b = self._validate(right, name="right")
        return float(np.sum(self.weight[:, None] * a * b))

    def norm(self, field: npt.ArrayLike) -> float:
        """Return the weighted Euclidean norm ``sqrt(<field, field>_w)``."""

        value = self.inner(field, field)
        if value < 0.0:
            raise ArithmeticError("negative weighted squared norm")
        return float(np.sqrt(value))


def build_realspace_operator(grid: AxisymmetricGrid) -> RealspaceOperator:
    """Assemble the E-26 flux coefficients for the Dirichlet-reduced system.

    The radial coefficients are computed from the grid coordinates only.  The
    axis coefficient produced by the general flux formula is cross-checked
    against the closed form ``8/dr^2`` of E-26b.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    if not grid.periodic_z:
        raise ValueError("build_realspace_operator requires periodic_z=True")

    dr = grid.dr
    dz = grid.dz
    # r_{i+1/2} for i = 0 .. n_r-2, taken from the stored coordinates so that
    # the faces are consistent with the grid rather than with an assumed
    # spacing formula.
    upper_face = 0.5 * (grid.r[:-1] + grid.r[1:])
    lower_face = np.concatenate(([0.0], upper_face[:-1]))
    volume = 0.25 * (upper_face**4 - lower_face**4)
    if not np.all(volume > 0.0):
        raise ValueError("radial control volumes must be positive")
    face_coefficient = upper_face**3 / dr

    axis_coefficient = float(face_coefficient[0] / volume[0])
    closed_form_axis = 8.0 / dr**2
    axis_defect = abs(axis_coefficient - closed_form_axis) / closed_form_axis
    if axis_defect > _AXIS_COEFFICIENT_TOLERANCE:
        raise ArithmeticError(
            "axis control volume does not reproduce the E-26b coefficient "
            f"8/dr^2: relative defect {axis_defect:.3e}"
        )

    diagonal = np.empty_like(volume)
    diagonal[0] = face_coefficient[0] / volume[0]
    diagonal[1:] = (face_coefficient[1:] + face_coefficient[:-1]) / volume[1:]
    diagonal += 2.0 / dz**2
    if not np.all(diagonal > 0.0):
        raise ArithmeticError("Jacobi diagonal must be strictly positive")

    return RealspaceOperator(
        grid=grid,
        face_coefficient=np.asarray(face_coefficient, dtype=np.float64),
        volume=np.asarray(volume, dtype=np.float64),
        weight=np.asarray(volume * dz, dtype=np.float64),
        diagonal=np.asarray(diagonal, dtype=np.float64),
    )


@dataclass(frozen=True)
class RealspacePoissonSolution:
    """Solution and iteration diagnostics for one real-space solve."""

    psi1: FloatArray
    iterations: int
    weighted_relative_residual: float
    metadata: dict[str, Any]


def _preconditioned_conjugate_gradient(
    operator: RealspaceOperator,
    right_hand_side: FloatArray,
    *,
    tol: float,
    max_iterations: int,
) -> tuple[FloatArray, int, float]:
    """Solve ``A x = b`` by Jacobi-preconditioned conjugate gradients.

    The stopping test uses the recursively updated residual because that is
    free, but every candidate acceptance is confirmed against the freshly
    recomputed residual ``b - A x``.  If the fresh residual disagrees, it
    replaces the recursion and the iteration continues; the value returned is
    always the freshly recomputed one, so the reported residual cannot be an
    artifact of accumulated drift.
    """

    inverse_diagonal = 1.0 / operator.diagonal[:, None]
    solution = np.zeros_like(right_hand_side)
    right_hand_side_norm = operator.norm(right_hand_side)
    if right_hand_side_norm == 0.0:
        return solution, 0, 0.0

    residual = right_hand_side.copy()
    preconditioned = residual * inverse_diagonal
    direction = preconditioned.copy()
    residual_dot = operator.inner(residual, preconditioned)
    achieved = 1.0

    for iteration in range(1, max_iterations + 1):
        curvature_field = operator.apply(direction)
        curvature = operator.inner(direction, curvature_field)
        if not curvature > 0.0:
            raise ArithmeticError(
                "non-positive curvature in the weighted inner product; the "
                "reduced operator is not positive definite as assumed"
            )
        step = residual_dot / curvature
        solution += step * direction
        residual -= step * curvature_field

        if operator.norm(residual) <= tol * right_hand_side_norm:
            true_residual = right_hand_side - operator.apply(solution)
            achieved = operator.norm(true_residual) / right_hand_side_norm
            if achieved <= tol:
                return solution, iteration, achieved
            # Residual replacement: restart the recursion from the exact
            # residual rather than reporting a drifted one as converged.
            residual = true_residual
            preconditioned = residual * inverse_diagonal
            direction = preconditioned.copy()
            residual_dot = operator.inner(residual, preconditioned)
            continue

        preconditioned = residual * inverse_diagonal
        updated_dot = operator.inner(residual, preconditioned)
        if not updated_dot > 0.0:
            raise ArithmeticError(
                "non-positive preconditioned residual product; the Jacobi "
                "preconditioner is not positive definite as assumed"
            )
        direction = preconditioned + (updated_dot / residual_dot) * direction
        residual_dot = updated_dot

    achieved = (
        operator.norm(right_hand_side - operator.apply(solution))
        / right_hand_side_norm
    )
    raise ConvergenceError(
        f"conjugate gradients did not reach the weighted relative residual "
        f"tolerance {tol:.3e} within {max_iterations} iterations; "
        f"achieved {achieved:.3e}"
    )


def solve_realspace_poisson(
    grid: AxisymmetricGrid,
    omega1: npt.ArrayLike,
    outer_dirichlet: npt.ArrayLike,
    *,
    tol: float = 1.0e-12,
    max_iterations: int = 5000,
) -> RealspacePoissonSolution:
    """Solve ``-L5 psi1 = omega1`` in real space by preconditioned CG.

    Parameters
    ----------
    grid:
        Uniform :class:`~ns_certificate_lab.grid.AxisymmetricGrid` including
        ``r=0`` and ``r=r_max``.  ``periodic_z`` must be true and the axial
        grid must omit its right endpoint (enforced by the grid class).
    omega1:
        Real, finite nodal right-hand side of shape ``grid.shape``.  The outer
        radial row is not used by the solve because it is replaced by the
        Dirichlet equation, but it must still be finite.
    outer_dirichlet:
        Real, finite scalar or array of shape ``(grid.nz,)`` giving
        ``psi1(r_max, z)``.
    tol:
        Target weighted relative residual, measured as
        ``||b - A x||_w / ||b||_w`` with ``w_i = V_i dz``.
    max_iterations:
        Hard cap on conjugate-gradient iterations.  Exceeding it raises
        :class:`ConvergenceError` rather than returning an unconverged field.

    Notes
    -----
    The nodal right-hand side is used as an approximation of the control-volume
    average, exactly as in E-26; this is second-order consistent for smooth
    even fields but is not exact cell integration.  The weighted residual norm
    gives the axis cell a very small weight because ``V_0 = dr^4/64``; the
    unweighted maximum-norm residual is therefore also reported in the
    metadata and should be inspected alongside it.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    if not grid.periodic_z:
        raise ValueError("solve_realspace_poisson requires periodic_z=True")
    if isinstance(tol, bool) or not isinstance(tol, (int, float)):
        raise TypeError("tol must be a real number")
    tol = float(tol)
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be positive and finite")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int):
        raise TypeError("max_iterations must be an integer")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least one")

    _reject_complex(omega1, name="omega1")
    omega = grid.validate_field(omega1, name="omega1")
    trace = _outer_trace(grid, outer_dirichlet)

    operator = build_realspace_operator(grid)
    last_unknown = grid.nr - 2

    # Dirichlet elimination: the face at r_{n_r-3/2} couples the last unknown
    # row to the prescribed outer value, so that coupling moves to the
    # right-hand side as +a^+ g (E-26c).
    right_hand_side = np.array(omega[:-1], dtype=np.float64, copy=True)
    outer_coupling = (
        operator.face_coefficient[last_unknown] / operator.volume[last_unknown]
    )
    right_hand_side[last_unknown] += outer_coupling * trace

    interior, iterations, achieved = _preconditioned_conjugate_gradient(
        operator,
        right_hand_side,
        tol=tol,
        max_iterations=max_iterations,
    )

    psi = np.empty(grid.shape, dtype=np.float64)
    psi[:-1] = interior
    psi[-1] = trace
    if not np.all(np.isfinite(psi)):
        raise FloatingPointError(
            "real-space Poisson solve produced a non-finite field"
        )

    algebraic_residual = right_hand_side - operator.apply(interior)
    right_hand_side_scale = max(
        1.0, float(np.max(np.abs(right_hand_side)))
    )
    metadata: dict[str, Any] = {
        "equation": "-(d_rr + 3/r d_r + d_zz) psi1 = omega1",
        "radial_discretization": (
            "second-order r^3-flux finite volume transcribed from E-26a/E-26b"
        ),
        "axial_discretization": (
            "three-point periodic centered difference applied with numpy.roll; "
            "no transform is used"
        ),
        "linear_solver": (
            "matrix-free Jacobi-preconditioned conjugate gradients in the "
            "V_i*dz weighted inner product"
        ),
        "rhs_interpretation": (
            "nodal omega1 values approximate radial control-volume averages; "
            "second-order consistent for smooth even fields, not exact cell "
            "integration"
        ),
        "shared_dependencies": (
            "AxisymmetricGrid data layout, binary64 arithmetic and the E-26 "
            "discretization mathematics remain common to all three solvers"
        ),
        "nr": grid.nr,
        "nz": grid.nz,
        "dr": grid.dr,
        "dz": grid.dz,
        "z_period": grid.z_period,
        "tolerance": tol,
        "max_iterations": max_iterations,
        "iterations": iterations,
        "weighted_relative_residual": achieved,
        "axis_radial_coefficient": float(
            operator.face_coefficient[0] / operator.volume[0]
        ),
        "axis_radial_coefficient_dimensionless": float(
            operator.face_coefficient[0] / operator.volume[0] * grid.dr**2
        ),
        "axis_cell_volume": float(operator.volume[0]),
        "outer_boundary_max_abs_defect": float(
            np.max(np.abs(psi[-1] - trace))
        ),
        "algebraic_residual_max_abs": float(np.max(np.abs(algebraic_residual))),
        "algebraic_residual_max_abs_relative": float(
            np.max(np.abs(algebraic_residual)) / right_hand_side_scale
        ),
        "algebraic_residual_rms": float(
            np.sqrt(np.mean(algebraic_residual**2))
        ),
        "interpretation": (
            "finite-dimensional elliptic solve diagnostics; not a continuum "
            "error bound, not an outer-boundary-truncation bound, and not "
            "singularity evidence"
        ),
    }

    return RealspacePoissonSolution(
        psi1=psi,
        iterations=iterations,
        weighted_relative_residual=achieved,
        metadata=metadata,
    )
