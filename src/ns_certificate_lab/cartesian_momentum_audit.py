r"""Pressure-free Cartesian momentum audit and a scoped pressure recovery.

This module closes the half of the primitive Navier--Stokes check that does
not need a pressure field.  For the exact incompressible equations the
momentum residual

.. math::

   R(u) := \partial_t u + (u\cdot\nabla)u - \nu\Delta u

must equal :math:`-\nabla p`, so its curl vanishes identically.  Two *distinct*
discrete forms of that statement are evaluated here:

``curl R``
    The curl of the assembled momentum residual.  The nonlinear term enters as
    ``curl((u.grad)u)``: one curl applied to the advection field.

``S := omega_t + (u.grad)omega - (omega.grad)u - nu*Laplacian(omega)``
    The Cartesian vorticity-transport residual with ``omega = curl u``.  The
    nonlinear term enters as the difference of two products of first
    derivatives.

The two agree in the continuum whenever ``div u = 0``, but they are *different*
discretizations, so reporting both is a genuine cross-check rather than a
restatement (TM-14).  Neither needs a pressure, which is exactly why they can
be applied to real checkpoints that store none.

Every stencil used below is taken from
:mod:`ns_certificate_lab.cartesian_validation` through its public surface.  In
particular :func:`cartesian_scalar_laplacian` is expressed through the public
``cartesian_vector_laplacian`` (which is componentwise) rather than by adding a
private second-derivative of its own, and the least-squares pressure recovery
uses the exact algebraic transpose of the module's own ``cartesian_gradient``
(verified against it in the test suite) rather than a new discretization.  This
module imports no cylindrical operator, PDE or Poisson module.

Relativization (TM-09).  Every defect below is a signed sum of terms that
cancel.  A raw magnitude is therefore meaningless on its own, so each residual
is reported together with the pointwise sum of the magnitudes of exactly the
terms that were summed to produce it, and the ratio is taken against the
maximum of that sum *over the same region* in which the defect norm is taken.
Denominators taken over a different region than the numerator would flatter the
interior numbers, so they are not used here.

Nothing in this module is an interval enclosure, a discretization error bound,
or evidence about singularity formation.  It is floating-point binary64
arithmetic with tolerances.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np
import numpy.typing as npt

from .cartesian_validation import (
    UniformCartesianGrid,
    cartesian_curl,
    cartesian_divergence,
    cartesian_gradient,
    cartesian_vector_laplacian,
)

ScalarArray = npt.NDArray[np.float64]
VectorArray = npt.NDArray[np.float64]

_TINY = float(np.finfo(np.float64).tiny)

#: Documented meaning of every denominator this module divides by.  These are
#: cancellation scales in the sense of TM-09: the pointwise sum of the
#: magnitudes of exactly the terms that were added together, maximized over the
#: region in which the numerator norm was taken.
DENOMINATORS: dict[str, str] = {
    "curl_momentum_term_sum_max": (
        "max over the region and over the three components of "
        "|curl(u_t)| + |curl((u.grad)u)| + |nu*curl(Laplacian u)|, i.e. the "
        "magnitudes of exactly the three curled terms whose signed sum is the "
        "reported curl-of-momentum residual"
    ),
    "vorticity_transport_term_sum_max": (
        "max over the region and over the three components of "
        "|omega_t| + |(u.grad)omega| + |(omega.grad)u| + |nu*Laplacian omega|, "
        "i.e. the magnitudes of exactly the four terms whose signed sum is the "
        "reported vorticity-transport residual"
    ),
    "momentum_term_sum_max": (
        "max over the region and over the three components of "
        "|u_t| + |(u.grad)u| + |nu*Laplacian u|.  This is the scale of the "
        "momentum residual R itself; R is NOT a defect (it equals -grad p), so "
        "this is the denominator of the pressure-consistency diagnostic only."
    ),
    "pressure_consistency_term_sum_max": (
        "max over the region and over the three components of "
        "|u_t| + |(u.grad)u| + |nu*Laplacian u| + |grad p_recovered|, the "
        "cancellation scale of R + grad p"
    ),
    "velocity_gradient_row_sum_max": (
        "max over the region of max_i sum_j |d_j u_i|, the largest absolute "
        "row sum of the Cartesian velocity gradient.  Recorded as the common "
        "first-derivative scale of the reconstructed field; the divergence "
        "ratio divides by it."
    ),
    "divergence_term_sum_max": (
        "max over the region of |d_x u_x| + |d_y u_y| + |d_z u_z|, the "
        "cancellation scale of div u"
    ),
    "vorticity_max_abs": (
        "max over the region and components of |curl u|.  A signal scale, not "
        "a cancellation scale."
    ),
    "velocity_max_abs": (
        "max over the region and components of |u|.  Scale context only."
    ),
}


# --------------------------------------------------------------------------
# regions
# --------------------------------------------------------------------------


def region_window(
    grid: UniformCartesianGrid,
    margin: int,
) -> tuple[slice, slice, slice]:
    """Return the slice that trims ``margin`` points off every *non-periodic* axis.

    Periodic axes keep every point including their wrap seam, exactly as the
    snapshot audit does.
    """

    if isinstance(margin, bool) or not isinstance(margin, (int, np.integer)):
        raise ValueError("interior_margin must be a nonnegative integer")
    active = int(margin)
    if active < 0:
        raise ValueError("interior_margin must be a nonnegative integer")
    if active and any(
        size <= 2 * active
        for size, periodic in zip(grid.shape, grid.periodic)
        if not periodic
    ):
        raise ValueError("interior_margin leaves an empty norm region")
    return tuple(  # type: ignore[return-value]
        slice(None) if periodic or active == 0 else slice(active, -active)
        for periodic in grid.periodic
    )


def region_offset(grid: UniformCartesianGrid, margin: int) -> tuple[int, int, int]:
    """Index shift from a trimmed region back to whole-box indices."""

    active = int(margin)
    return tuple(  # type: ignore[return-value]
        0 if periodic or active == 0 else active for periodic in grid.periodic
    )


# --------------------------------------------------------------------------
# derived Cartesian operators, all built from the public audit stencils
# --------------------------------------------------------------------------


def cartesian_scalar_laplacian(
    grid: UniformCartesianGrid,
    scalar: npt.ArrayLike,
) -> ScalarArray:
    """Return ``Laplacian(scalar)`` using only the public vector Laplacian.

    ``cartesian_vector_laplacian`` is componentwise, so applying it to the
    vector whose three components are all ``scalar`` returns three identical
    copies of the scalar Laplacian.  No new second-derivative stencil is
    introduced.
    """

    values = grid.validate_scalar(scalar, name="scalar field")
    stacked = np.stack((values, values, values), axis=0)
    return cartesian_vector_laplacian(grid, stacked)[0]


def cartesian_directional_derivative(
    grid: UniformCartesianGrid,
    *,
    advecting: npt.ArrayLike,
    field: npt.ArrayLike,
) -> VectorArray:
    """Return ``(a.grad) f`` for vector ``a`` and vector ``f``."""

    carrier = grid.validate_vector(advecting, name="advecting field")
    target = grid.validate_vector(field, name="advected field")
    result = np.empty_like(target)
    for component in range(3):
        gradient = cartesian_gradient(grid, target[component])
        result[component] = sum(
            (carrier[direction] * gradient[direction] for direction in range(3)),
            start=np.zeros(grid.shape, dtype=np.float64),
        )
    return result


def _check_viscosity(viscosity: float) -> float:
    value = float(viscosity)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("viscosity must be finite and nonnegative")
    return value


# --------------------------------------------------------------------------
# the two pressure-free residuals
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class MomentumTerms:
    """The three pressure-free terms of ``R`` and their signed sum.

    ``viscous`` stores ``-nu*Laplacian(u)`` so that
    ``residual = time_derivative + advection + viscous`` exactly, matching the
    sign convention of :class:`cartesian_validation.PrimitiveNSResidual`.
    ``R`` is *not* a defect: for an exact solution it equals ``-grad p``.
    """

    time_derivative: VectorArray
    advection: VectorArray
    viscous: VectorArray
    residual: VectorArray
    term_sum: VectorArray


def momentum_terms(
    grid: UniformCartesianGrid,
    *,
    velocity: npt.ArrayLike,
    velocity_t: npt.ArrayLike,
    viscosity: float,
    viscous_sign: float = 1.0,
) -> MomentumTerms:
    """Assemble ``R = u_t + (u.grad)u - nu*Laplacian(u)`` term by term.

    ``viscous_sign`` exists solely for the sign-flip fault injection: the
    physical value is ``+1.0`` and any other value deliberately breaks the
    equation being evaluated.
    """

    state = grid.validate_vector(velocity, name="velocity")
    rate = grid.validate_vector(velocity_t, name="velocity_t")
    nu = _check_viscosity(viscosity)
    sign = float(viscous_sign)
    if not math.isfinite(sign):
        raise ValueError("viscous_sign must be finite")
    advection = cartesian_directional_derivative(grid, advecting=state, field=state)
    viscous = -sign * nu * cartesian_vector_laplacian(grid, state)
    residual = rate + advection + viscous
    term_sum = np.abs(rate) + np.abs(advection) + np.abs(viscous)
    return MomentumTerms(
        time_derivative=rate,
        advection=advection,
        viscous=viscous,
        residual=residual,
        term_sum=term_sum,
    )


@dataclass(frozen=True)
class CurlMomentumResidual:
    """``curl R``, which must vanish because ``R = -grad p``."""

    curl_time_derivative: VectorArray
    curl_advection: VectorArray
    curl_viscous: VectorArray
    residual: VectorArray
    term_sum: VectorArray
    momentum: MomentumTerms


def curl_momentum_residual(
    grid: UniformCartesianGrid,
    *,
    velocity: npt.ArrayLike,
    velocity_t: npt.ArrayLike,
    viscosity: float,
    viscous_sign: float = 1.0,
) -> CurlMomentumResidual:
    """Return ``curl(u_t) + curl((u.grad)u) - nu*curl(Laplacian u)``.

    The curl is applied to each momentum term separately so that the pointwise
    cancellation scale of the sum is available (TM-09).  ``viscous_sign`` is
    the fault-injection hook documented on :func:`momentum_terms`.
    """

    terms = momentum_terms(
        grid,
        velocity=velocity,
        velocity_t=velocity_t,
        viscosity=viscosity,
        viscous_sign=viscous_sign,
    )
    curl_time = cartesian_curl(grid, terms.time_derivative)
    curl_advection = cartesian_curl(grid, terms.advection)
    curl_viscous = cartesian_curl(grid, terms.viscous)
    residual = curl_time + curl_advection + curl_viscous
    term_sum = np.abs(curl_time) + np.abs(curl_advection) + np.abs(curl_viscous)
    return CurlMomentumResidual(
        curl_time_derivative=curl_time,
        curl_advection=curl_advection,
        curl_viscous=curl_viscous,
        residual=residual,
        term_sum=term_sum,
        momentum=terms,
    )


@dataclass(frozen=True)
class VorticityTransportResidual:
    """``S = omega_t + (u.grad)omega - (omega.grad)u - nu*Laplacian(omega)``."""

    vorticity: VectorArray
    time_derivative: VectorArray
    advection: VectorArray
    stretching: VectorArray
    viscous: VectorArray
    residual: VectorArray
    term_sum: VectorArray


def vorticity_transport_residual(
    grid: UniformCartesianGrid,
    *,
    velocity: npt.ArrayLike,
    vorticity: npt.ArrayLike,
    vorticity_t: npt.ArrayLike,
    viscosity: float,
    viscous_sign: float = 1.0,
) -> VorticityTransportResidual:
    """Assemble the Cartesian vorticity-transport residual term by term.

    ``vorticity`` is passed in rather than recomputed so that the caller owns
    the choice of how ``omega`` and ``omega_t`` were formed from the saved
    states, and so that a fault injection can corrupt exactly one term.

    ``viscous_sign`` exists solely for the sign-flip fault injection: the
    physical value is ``+1.0`` and any other value deliberately breaks the
    equation.  ``stretching`` stores ``-(omega.grad)u`` and ``viscous`` stores
    ``-viscous_sign*nu*Laplacian(omega)`` so that ``residual`` is the plain sum
    of the four stored fields.
    """

    state = grid.validate_vector(velocity, name="velocity")
    omega = grid.validate_vector(vorticity, name="vorticity")
    omega_t = grid.validate_vector(vorticity_t, name="vorticity_t")
    nu = _check_viscosity(viscosity)
    sign = float(viscous_sign)
    if not math.isfinite(sign):
        raise ValueError("viscous_sign must be finite")
    advection = cartesian_directional_derivative(grid, advecting=state, field=omega)
    stretching = -cartesian_directional_derivative(
        grid, advecting=omega, field=state
    )
    viscous = -sign * nu * cartesian_vector_laplacian(grid, omega)
    residual = omega_t + advection + stretching + viscous
    term_sum = (
        np.abs(omega_t) + np.abs(advection) + np.abs(stretching) + np.abs(viscous)
    )
    return VorticityTransportResidual(
        vorticity=omega,
        time_derivative=omega_t,
        advection=advection,
        stretching=stretching,
        viscous=viscous,
        residual=residual,
        term_sum=term_sum,
    )


# --------------------------------------------------------------------------
# least-squares pressure recovery (scoped diagnostic, see the module docstring)
# --------------------------------------------------------------------------


def adjoint_first_derivative(
    grid: UniformCartesianGrid,
    field: npt.ArrayLike,
    axis: int,
) -> ScalarArray:
    r"""Return ``D_axis^T q`` for the audit module's own first-derivative ``D``.

    This is the exact algebraic transpose of the matrix that
    :func:`cartesian_validation.cartesian_gradient` applies along ``axis``: the
    same coefficients, transposed, and nothing else.  It is not an independent
    discretization and it is verified against ``cartesian_gradient`` by the
    inner-product identity ``<D u, q> == <u, D^T q>`` in the test suite.

    Periodic axis: ``D`` is antisymmetric, so ``D^T = -D``.
    Non-periodic axis: the interior rows contribute the transposed centered
    pair and the two one-sided closure rows contribute their three
    coefficients each.
    """

    if isinstance(axis, bool) or axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    direction = int(axis)
    values = grid.validate_scalar(field, name="adjoint operand")
    spacing = grid.spacings[direction]
    if grid.periodic[direction]:
        return -(
            np.roll(values, -1, axis=direction)
            - np.roll(values, 1, axis=direction)
        ) / (2.0 * spacing)

    count = grid.shape[direction]
    if count < 4:
        raise ValueError("the adjoint closure needs at least four points per axis")
    scale = 1.0 / (2.0 * spacing)
    result = np.zeros_like(values)

    def window(start: Any, stop: Any = None) -> tuple[Any, ...]:
        index: list[Any] = [slice(None)] * 3
        index[direction] = start if stop is None else slice(start, stop)
        return tuple(index)

    interior = values[window(1, count - 1)]
    result[window(2, count)] += scale * interior
    result[window(0, count - 2)] -= scale * interior

    low = values[window(0)]
    result[window(0)] += -3.0 * scale * low
    result[window(1)] += 4.0 * scale * low
    result[window(2)] += -1.0 * scale * low

    high = values[window(count - 1)]
    result[window(count - 1)] += 3.0 * scale * high
    result[window(count - 2)] += -4.0 * scale * high
    result[window(count - 3)] += 1.0 * scale * high
    return result


def gradient_null_space(grid: UniformCartesianGrid) -> tuple[ScalarArray, ...]:
    r"""Return an orthonormal basis of ``ker(grad_h)`` on this grid.

    Along a non-periodic axis the one-sided closure rows force a null vector to
    be constant.  Along a periodic axis with an even number of points the
    centered wrapped difference additionally annihilates the alternating
    (checkerboard) mode.  The kernel of the full gradient is the tensor product
    of the per-axis kernels, so it has dimension ``2**k`` where ``k`` is the
    number of periodic axes with an even point count.
    """

    axis_bases: list[list[ScalarArray]] = []
    for direction in range(3):
        count = grid.shape[direction]
        vectors = [np.ones(count, dtype=np.float64)]
        if grid.periodic[direction] and count % 2 == 0:
            vectors.append((-1.0) ** np.arange(count, dtype=np.float64))
        axis_bases.append(vectors)

    basis: list[ScalarArray] = []
    for first in axis_bases[0]:
        for second in axis_bases[1]:
            for third in axis_bases[2]:
                vector = (
                    first[:, None, None] * second[None, :, None] * third[None, None, :]
                )
                norm = float(np.sqrt(np.sum(vector * vector)))
                basis.append(np.asarray(vector / norm, dtype=np.float64))
    return tuple(basis)


def _project_out(field: ScalarArray, basis: tuple[ScalarArray, ...]) -> ScalarArray:
    result = field
    for vector in basis:
        result = result - float(np.sum(result * vector)) * vector
    return result


@dataclass(frozen=True)
class PressureRecovery:
    """Least-squares kinematic pressure recovered from ``R = -grad p``.

    ``pressure`` is determined only up to the kernel of the discrete gradient
    (see :func:`gradient_null_space`); the returned field is the representative
    orthogonal to that kernel.  ``consistency`` is ``R + grad p_recovered``,
    which the exact continuum relation would make zero.
    """

    pressure: ScalarArray
    pressure_gradient: VectorArray
    consistency: VectorArray
    term_sum: VectorArray
    iterations: int
    converged: bool
    normal_equation_relative_residual: float
    null_space_dimension: int


def recover_pressure(
    grid: UniformCartesianGrid,
    momentum_residual: npt.ArrayLike,
    *,
    momentum_term_sum: npt.ArrayLike | None = None,
    tolerance: float = 1.0e-10,
    max_iterations: int = 20000,
) -> PressureRecovery:
    r"""Recover ``p`` from ``R = -grad p`` in the discrete least-squares sense.

    The exact relation is ``grad p = -R``.  That linear system is
    overdetermined (three equations per unknown), so ``p`` is chosen to
    minimize ``||grad_h p + R||_2^2`` over the whole box, whose normal
    equations are the discrete Poisson problem

    ``sum_a D_a^T D_a p = - sum_a D_a^T R_a``,

    i.e. ``div grad p = -div R`` written with the audit module's own gradient
    and its exact transpose.  No boundary condition is invented: the
    least-squares formulation supplies the natural one implied by the audit's
    one-sided closures.

    **Scope.**  The audit box is periodic in ``z`` but is an artificial cut in
    ``x`` and ``y``: the true pressure of the source problem is not determined
    there by data inside the box.  The recovered field is therefore a *scoped
    diagnostic*, never a validated pressure, and only its interior consistency
    is meaningful.  Callers must record it as record-only.
    """

    residual_field = grid.validate_vector(
        momentum_residual, name="momentum_residual"
    )
    limit = float(tolerance)
    if not math.isfinite(limit) or limit <= 0.0:
        raise ValueError("tolerance must be positive and finite")
    if (
        isinstance(max_iterations, bool)
        or not isinstance(max_iterations, int)
        or max_iterations < 1
    ):
        raise ValueError("max_iterations must be a positive integer")

    basis = gradient_null_space(grid)

    def normal_operator(field: ScalarArray) -> ScalarArray:
        gradient = cartesian_gradient(grid, field)
        return sum(
            (
                adjoint_first_derivative(grid, gradient[axis], axis)
                for axis in range(3)
            ),
            start=np.zeros(grid.shape, dtype=np.float64),
        )

    rhs = -sum(
        (
            adjoint_first_derivative(grid, residual_field[axis], axis)
            for axis in range(3)
        ),
        start=np.zeros(grid.shape, dtype=np.float64),
    )
    rhs = _project_out(rhs, basis)
    rhs_norm = float(np.sqrt(np.sum(rhs * rhs)))

    pressure = np.zeros(grid.shape, dtype=np.float64)
    iterations = 0
    converged = True
    relative = 0.0
    if rhs_norm > _TINY:
        current = rhs.copy()
        direction = current.copy()
        squared = float(np.sum(current * current))
        converged = False
        for step in range(1, int(max_iterations) + 1):
            applied = _project_out(normal_operator(direction), basis)
            curvature = float(np.sum(direction * applied))
            if not math.isfinite(curvature) or curvature <= 0.0:
                iterations = step
                break
            alpha = squared / curvature
            pressure = pressure + alpha * direction
            current = current - alpha * applied
            updated = float(np.sum(current * current))
            iterations = step
            if math.sqrt(updated) <= limit * rhs_norm:
                converged = True
                break
            direction = current + (updated / squared) * direction
            squared = updated
        relative = (
            float(np.sqrt(np.sum(current * current))) / rhs_norm
            if rhs_norm > _TINY
            else 0.0
        )
    pressure = _project_out(pressure, basis)
    gradient = cartesian_gradient(grid, pressure)
    consistency = residual_field + gradient
    if momentum_term_sum is None:
        term_sum = np.abs(residual_field) + np.abs(gradient)
    else:
        term_sum = (
            grid.validate_vector(momentum_term_sum, name="momentum_term_sum")
            + np.abs(gradient)
        )
    return PressureRecovery(
        pressure=pressure,
        pressure_gradient=gradient,
        consistency=consistency,
        term_sum=term_sum,
        iterations=int(iterations),
        converged=bool(converged),
        normal_equation_relative_residual=float(relative),
        null_space_dimension=len(basis),
    )


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def rms(values: npt.NDArray[np.float64]) -> float:
    array = np.asarray(values, dtype=np.float64)
    return float(np.sqrt(np.mean(array * array)))


def ratio(numerator: float, denominator: float) -> float | None:
    """Return ``numerator/denominator`` or ``None`` for a degenerate scale."""

    if not math.isfinite(denominator) or denominator <= _TINY:
        return None
    value = numerator / denominator
    return float(value) if math.isfinite(value) else None


def locate_maximum(
    field: npt.NDArray[np.float64],
    mesh: tuple[ScalarArray, ScalarArray, ScalarArray],
    *,
    offset: tuple[int, int, int],
    vector: bool,
) -> dict[str, Any]:
    """Locate the largest magnitude of ``field`` in whole-box coordinates."""

    index = np.unravel_index(int(np.argmax(np.abs(field))), field.shape)
    component = int(index[0]) if vector else None
    spatial = tuple(int(value) for value in (index[1:] if vector else index))
    absolute = tuple(int(value) + int(shift) for value, shift in zip(spatial, offset))
    x, y, z = mesh
    return {
        "component": component,
        "index": list(absolute),
        "x": float(x[absolute]),
        "y": float(y[absolute]),
        "z": float(z[absolute]),
        "radius": float(math.hypot(float(x[absolute]), float(y[absolute]))),
        "value": float(field[index]),
    }


def defect_report(
    grid: UniformCartesianGrid,
    *,
    defect: npt.NDArray[np.float64],
    term_sum: npt.NDArray[np.float64],
    window: tuple[slice, slice, slice],
    offset: tuple[int, int, int],
    mesh: tuple[ScalarArray, ScalarArray, ScalarArray],
    vector: bool = True,
) -> dict[str, Any]:
    """Absolute norms, the region denominator and the relative ratios."""

    selector = (slice(None), *window) if vector else window
    region_defect = defect[selector]
    region_terms = term_sum[selector]
    denominator = float(np.max(np.abs(region_terms)))
    absolute_rms = rms(region_defect)
    absolute_max = float(np.max(np.abs(region_defect)))
    return {
        "rms": absolute_rms,
        "max": absolute_max,
        "term_sum_max": denominator,
        "relative_rms": ratio(absolute_rms, denominator),
        "relative_max": ratio(absolute_max, denominator),
        "argmax": locate_maximum(region_defect, mesh, offset=offset, vector=vector),
        "term_sum_argmax": locate_maximum(
            region_terms, mesh, offset=offset, vector=vector
        ),
    }


def velocity_scales(
    grid: UniformCartesianGrid,
    velocity: npt.ArrayLike,
    vorticity: npt.ArrayLike,
    *,
    window: tuple[slice, slice, slice],
) -> dict[str, float]:
    """Report the first-derivative and signal scales of one reconstruction."""

    state = grid.validate_vector(velocity, name="velocity")
    omega = grid.validate_vector(vorticity, name="vorticity")
    gradient = np.stack(
        [cartesian_gradient(grid, state[component]) for component in range(3)],
        axis=0,
    )
    selector = (slice(None), slice(None), *window)
    region_gradient = gradient[selector]
    vector_window = (slice(None), *window)
    return {
        "velocity_max_abs": float(np.max(np.abs(state[vector_window]))),
        "vorticity_max_abs": float(np.max(np.abs(omega[vector_window]))),
        "velocity_gradient_row_sum_max": float(
            np.max(np.abs(region_gradient).sum(axis=1))
        ),
        "divergence_term_sum_max": float(
            np.max(
                np.abs(gradient[0, 0][window])
                + np.abs(gradient[1, 1][window])
                + np.abs(gradient[2, 2][window])
            )
        ),
    }


def divergence_report(
    grid: UniformCartesianGrid,
    velocity: npt.ArrayLike,
    *,
    window: tuple[slice, slice, slice],
    offset: tuple[int, int, int],
    mesh: tuple[ScalarArray, ScalarArray, ScalarArray],
) -> dict[str, Any]:
    """E-02 divergence of the reconstruction, relative to its own term sum."""

    state = grid.validate_vector(velocity, name="velocity")
    gradient = np.stack(
        [cartesian_gradient(grid, state[component]) for component in range(3)],
        axis=0,
    )
    term_sum = (
        np.abs(gradient[0, 0]) + np.abs(gradient[1, 1]) + np.abs(gradient[2, 2])
    )
    return defect_report(
        grid,
        defect=cartesian_divergence(grid, state),
        term_sum=term_sum,
        window=window,
        offset=offset,
        mesh=mesh,
        vector=False,
    )


def observed_order(
    coarse_value: float | None,
    fine_value: float | None,
    spacing_ratio: float,
) -> float | None:
    """Return ``log(coarse/fine)/log(spacing_ratio)`` when it is well defined."""

    if (
        coarse_value is None
        or fine_value is None
        or not math.isfinite(float(coarse_value))
        or not math.isfinite(float(fine_value))
        or float(coarse_value) <= 0.0
        or float(fine_value) <= 0.0
        or spacing_ratio <= 1.0
    ):
        return None
    value = math.log(float(coarse_value) / float(fine_value)) / math.log(spacing_ratio)
    return float(value) if math.isfinite(value) else None


def finite_values(value: Any) -> list[float]:
    """Flatten every numeric leaf of a JSON-shaped object."""

    if isinstance(value, bool) or value is None:
        return []
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, Mapping):
        return [item for entry in value.values() for item in finite_values(entry)]
    if isinstance(value, (list, tuple)):
        return [item for entry in value for item in finite_values(entry)]
    return []
