r"""Independent uniform-Cartesian validation of primitive Navier--Stokes fields.

This module is deliberately self-contained apart from NumPy.  In particular,
it does not import the cylindrical grid, ``operators.py``, ``pde.py``, or any
cylindrical finite-difference result.  It is therefore suitable as a
cross-implementation audit path for an axisymmetric candidate reconstructed in
physical Cartesian coordinates.

Arrays use the component-first convention ``(3, n_x, n_y, n_z)``.  The
primitive momentum residual is

.. math::

   R = \partial_t u + (u\cdot\nabla)u + \nabla p
       - \nu\Delta u.

Here ``p`` is kinematic pressure (physical pressure divided by constant
density), ``velocity_t`` is supplied at the same physical time as
``velocity``, and ``viscosity`` is the kinematic viscosity.  Pressure has no
time derivative in this residual.  Incompressibility is reported separately
as ``divergence``.  An optional force per unit mass ``f`` is retained only as a
manufactured-solution oracle: ``defect_against_body_force = R - f``; it never
changes the meaning of the unforced ``total`` residual.

The second-order stencils are intentionally different from the existing
cylindrical implementation.  Periodic axes use centered wrapped differences;
non-periodic axes use centered interior differences and second-order one-sided
boundary closures.  These floating-point checks are not interval proofs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import numpy.typing as npt

ScalarArray = npt.NDArray[np.float64]
VectorArray = npt.NDArray[np.float64]
ProfileEvaluator = Callable[[ScalarArray, ScalarArray], npt.ArrayLike]


def _finite_positive(value: float, *, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be positive and finite") from error
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


@dataclass(frozen=True)
class UniformCartesianGrid:
    """Uniform three-dimensional Cartesian grid.

    ``periods[i]`` must be supplied for a periodic axis and must equal
    ``n_i * h_i`` because the endpoint is omitted.  It must be ``None`` for a
    non-periodic axis.  Prefer :meth:`uniform`, which constructs these
    conventions without ambiguity.
    """

    x: ScalarArray
    y: ScalarArray
    z: ScalarArray
    periodic: tuple[bool, bool, bool] = (False, False, False)
    periods: tuple[float | None, float | None, float | None] = (
        None,
        None,
        None,
    )

    def __post_init__(self) -> None:
        if len(self.periodic) != 3 or len(self.periods) != 3:
            raise ValueError("periodic and periods must each contain three entries")
        coordinates: list[ScalarArray] = []
        normalized_periods: list[float | None] = []
        for axis, (raw, is_periodic, raw_period) in enumerate(
            zip((self.x, self.y, self.z), self.periodic, self.periods)
        ):
            values = np.array(raw, dtype=np.float64, copy=True)
            if values.ndim != 1 or values.size < 4:
                raise ValueError(
                    f"Cartesian axis {axis} must be one-dimensional with at least "
                    "four points"
                )
            if not np.all(np.isfinite(values)):
                raise ValueError(f"Cartesian axis {axis} must contain finite values")
            increments = np.diff(values)
            if np.any(increments <= 0.0):
                raise ValueError(
                    f"Cartesian axis {axis} coordinates must be strictly increasing"
                )
            if not np.allclose(
                increments,
                increments[0],
                rtol=1.0e-12,
                atol=1.0e-14,
            ):
                raise ValueError(f"Cartesian axis {axis} must be uniform")
            if is_periodic:
                period = _finite_positive(raw_period, name=f"periods[{axis}]")
                inferred = float(values.size * increments[0])
                if not np.isclose(period, inferred, rtol=1.0e-11, atol=1.0e-13):
                    raise ValueError(
                        f"periodic Cartesian axis {axis} must omit its endpoint: "
                        "period must equal n * spacing"
                    )
                normalized_periods.append(period)
            else:
                if raw_period is not None:
                    raise ValueError(
                        f"periods[{axis}] must be None for a non-periodic axis"
                    )
                normalized_periods.append(None)
            values.setflags(write=False)
            coordinates.append(values)

        object.__setattr__(self, "x", coordinates[0])
        object.__setattr__(self, "y", coordinates[1])
        object.__setattr__(self, "z", coordinates[2])
        object.__setattr__(self, "periodic", tuple(bool(v) for v in self.periodic))
        object.__setattr__(self, "periods", tuple(normalized_periods))

    @classmethod
    def uniform(
        cls,
        *,
        shape: tuple[int, int, int],
        x_bounds: tuple[float, float],
        y_bounds: tuple[float, float],
        z_bounds: tuple[float, float],
        periodic: tuple[bool, bool, bool] = (False, False, False),
    ) -> "UniformCartesianGrid":
        """Construct a grid, omitting the upper endpoint on periodic axes."""

        if len(shape) != 3 or len(periodic) != 3:
            raise ValueError("shape and periodic must each contain three entries")
        axes: list[ScalarArray] = []
        periods: list[float | None] = []
        for axis, (count, bounds, is_periodic) in enumerate(
            zip(shape, (x_bounds, y_bounds, z_bounds), periodic)
        ):
            if isinstance(count, bool) or int(count) != count or int(count) < 4:
                raise ValueError(f"shape[{axis}] must be an integer >= 4")
            lower, upper = (float(bounds[0]), float(bounds[1]))
            if not (
                np.isfinite(lower)
                and np.isfinite(upper)
                and upper > lower
            ):
                raise ValueError(
                    f"bounds for Cartesian axis {axis} must be finite and increasing"
                )
            axes.append(
                np.linspace(
                    lower,
                    upper,
                    int(count),
                    endpoint=not is_periodic,
                    dtype=np.float64,
                )
            )
            periods.append(upper - lower if is_periodic else None)
        return cls(
            x=axes[0],
            y=axes[1],
            z=axes[2],
            periodic=periodic,
            periods=tuple(periods),
        )

    @property
    def shape(self) -> tuple[int, int, int]:
        return (int(self.x.size), int(self.y.size), int(self.z.size))

    @property
    def spacings(self) -> tuple[float, float, float]:
        return (
            float(self.x[1] - self.x[0]),
            float(self.y[1] - self.y[0]),
            float(self.z[1] - self.z[0]),
        )

    def mesh(self) -> tuple[ScalarArray, ScalarArray, ScalarArray]:
        """Return ``X, Y, Z`` arrays with shape :attr:`shape`."""

        return np.meshgrid(self.x, self.y, self.z, indexing="ij")

    def validate_scalar(
        self,
        field: npt.ArrayLike,
        *,
        name: str = "scalar field",
    ) -> ScalarArray:
        values = np.asarray(field, dtype=np.float64)
        if values.shape != self.shape:
            raise ValueError(f"{name} must have shape {self.shape}, got {values.shape}")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must contain only finite values")
        return values

    def validate_vector(
        self,
        field: npt.ArrayLike,
        *,
        name: str = "vector field",
    ) -> VectorArray:
        values = np.asarray(field, dtype=np.float64)
        expected = (3, *self.shape)
        if values.shape != expected:
            raise ValueError(f"{name} must have shape {expected}, got {values.shape}")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must contain only finite values")
        return values


def _validate_axis(axis: int) -> int:
    if isinstance(axis, bool) or axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    return int(axis)


def _first_derivative(
    grid: UniformCartesianGrid,
    field: npt.ArrayLike,
    axis: int,
) -> ScalarArray:
    """Second-order first derivative owned solely by this audit module."""

    direction = _validate_axis(axis)
    values = grid.validate_scalar(field)
    spacing = grid.spacings[direction]
    if grid.periodic[direction]:
        return (
            np.roll(values, -1, axis=direction)
            - np.roll(values, 1, axis=direction)
        ) / (2.0 * spacing)

    result = np.empty_like(values)
    interior = [slice(None)] * 3
    forward = [slice(None)] * 3
    backward = [slice(None)] * 3
    interior[direction] = slice(1, -1)
    forward[direction] = slice(2, None)
    backward[direction] = slice(None, -2)
    result[tuple(interior)] = (
        values[tuple(forward)] - values[tuple(backward)]
    ) / (2.0 * spacing)

    edge = [slice(None)] * 3
    one = [slice(None)] * 3
    two = [slice(None)] * 3
    edge[direction] = 0
    one[direction] = 1
    two[direction] = 2
    result[tuple(edge)] = (
        -3.0 * values[tuple(edge)]
        + 4.0 * values[tuple(one)]
        - values[tuple(two)]
    ) / (2.0 * spacing)

    edge[direction] = -1
    one[direction] = -2
    two[direction] = -3
    result[tuple(edge)] = (
        3.0 * values[tuple(edge)]
        - 4.0 * values[tuple(one)]
        + values[tuple(two)]
    ) / (2.0 * spacing)
    return result


def _second_derivative(
    grid: UniformCartesianGrid,
    field: npt.ArrayLike,
    axis: int,
) -> ScalarArray:
    """Second-order second derivative owned solely by this audit module."""

    direction = _validate_axis(axis)
    values = grid.validate_scalar(field)
    spacing_squared = grid.spacings[direction] ** 2
    if grid.periodic[direction]:
        return (
            np.roll(values, -1, axis=direction)
            - 2.0 * values
            + np.roll(values, 1, axis=direction)
        ) / spacing_squared

    result = np.empty_like(values)
    interior = [slice(None)] * 3
    forward = [slice(None)] * 3
    backward = [slice(None)] * 3
    interior[direction] = slice(1, -1)
    forward[direction] = slice(2, None)
    backward[direction] = slice(None, -2)
    result[tuple(interior)] = (
        values[tuple(forward)]
        - 2.0 * values[tuple(interior)]
        + values[tuple(backward)]
    ) / spacing_squared

    edge = [slice(None)] * 3
    one = [slice(None)] * 3
    two = [slice(None)] * 3
    three = [slice(None)] * 3
    edge[direction] = 0
    one[direction] = 1
    two[direction] = 2
    three[direction] = 3
    result[tuple(edge)] = (
        2.0 * values[tuple(edge)]
        - 5.0 * values[tuple(one)]
        + 4.0 * values[tuple(two)]
        - values[tuple(three)]
    ) / spacing_squared

    edge[direction] = -1
    one[direction] = -2
    two[direction] = -3
    three[direction] = -4
    result[tuple(edge)] = (
        2.0 * values[tuple(edge)]
        - 5.0 * values[tuple(one)]
        + 4.0 * values[tuple(two)]
        - values[tuple(three)]
    ) / spacing_squared
    return result


def cartesian_gradient(
    grid: UniformCartesianGrid,
    scalar: npt.ArrayLike,
) -> VectorArray:
    """Return ``(partial_x, partial_y, partial_z) scalar``."""

    values = grid.validate_scalar(scalar)
    return np.stack(
        [_first_derivative(grid, values, axis) for axis in range(3)],
        axis=0,
    )


def cartesian_divergence(
    grid: UniformCartesianGrid,
    velocity: npt.ArrayLike,
) -> ScalarArray:
    """Return the physical three-dimensional divergence ``∇·u``."""

    values = grid.validate_vector(velocity, name="velocity")
    return sum(
        (_first_derivative(grid, values[axis], axis) for axis in range(3)),
        start=np.zeros(grid.shape, dtype=np.float64),
    )


def cartesian_curl(
    grid: UniformCartesianGrid,
    velocity: npt.ArrayLike,
) -> VectorArray:
    """Return all three Cartesian components of ``∇×u``."""

    values = grid.validate_vector(velocity, name="velocity")
    return np.stack(
        (
            _first_derivative(grid, values[2], 1)
            - _first_derivative(grid, values[1], 2),
            _first_derivative(grid, values[0], 2)
            - _first_derivative(grid, values[2], 0),
            _first_derivative(grid, values[1], 0)
            - _first_derivative(grid, values[0], 1),
        ),
        axis=0,
    )


def cartesian_vector_laplacian(
    grid: UniformCartesianGrid,
    velocity: npt.ArrayLike,
) -> VectorArray:
    """Return the componentwise Cartesian vector Laplacian ``Δu``."""

    values = grid.validate_vector(velocity, name="velocity")
    result = np.empty_like(values)
    for component in range(3):
        result[component] = sum(
            (
                _second_derivative(grid, values[component], axis)
                for axis in range(3)
            ),
            start=np.zeros(grid.shape, dtype=np.float64),
        )
    return result


@dataclass(frozen=True)
class PrimitiveNSResidual:
    """Term-by-term primitive momentum residual on a Cartesian grid.

    Every vector-valued member has shape ``(3, n_x, n_y, n_z)``.  ``viscous``
    stores ``-nu * Laplacian(u)``.  ``total`` is the direct sum of the four
    primitive Navier--Stokes terms and therefore always denotes the unforced
    residual ``R``.  ``defect_against_body_force`` is ``R-f`` and is useful
    only for a forced manufactured solution.
    """

    time_derivative: VectorArray
    advection: VectorArray
    pressure_gradient: VectorArray
    viscous: VectorArray
    total: VectorArray
    body_force: VectorArray
    defect_against_body_force: VectorArray
    divergence: ScalarArray


def primitive_ns_residual(
    grid: UniformCartesianGrid,
    *,
    velocity: npt.ArrayLike,
    velocity_t: npt.ArrayLike,
    pressure: npt.ArrayLike,
    viscosity: float,
    body_force: npt.ArrayLike | None = None,
) -> PrimitiveNSResidual:
    """Compute every Cartesian component of the primitive NS residual.

    No temporal differencing is hidden here: callers must supply
    ``velocity_t``.  This keeps a temporal discretization independently
    auditable.  If ``body_force`` is omitted, the unforced equation is tested.
    """

    state = grid.validate_vector(velocity, name="velocity")
    time_derivative = grid.validate_vector(
        velocity_t,
        name="velocity_t",
    ).copy()
    kinematic_pressure = grid.validate_scalar(pressure, name="pressure")
    nu = float(viscosity)
    if not np.isfinite(nu) or nu < 0.0:
        raise ValueError("viscosity must be finite and nonnegative")
    if body_force is None:
        physical_force = np.zeros_like(state)
    else:
        physical_force = grid.validate_vector(body_force, name="body_force")

    advection = np.empty_like(state)
    for component in range(3):
        advection[component] = sum(
            (
                state[direction]
                * _first_derivative(grid, state[component], direction)
                for direction in range(3)
            ),
            start=np.zeros(grid.shape, dtype=np.float64),
        )
    pressure_gradient = cartesian_gradient(grid, kinematic_pressure)
    viscous = -nu * cartesian_vector_laplacian(grid, state)
    total = (
        time_derivative
        + advection
        + pressure_gradient
        + viscous
    )
    return PrimitiveNSResidual(
        time_derivative=time_derivative,
        advection=advection,
        pressure_gradient=pressure_gradient,
        viscous=viscous,
        total=total,
        body_force=physical_force.copy(),
        defect_against_body_force=total - physical_force,
        divergence=cartesian_divergence(grid, state),
    )


@dataclass(frozen=True)
class RegularizedAxisymmetricProfile:
    """Point evaluators for ``u1``, ``psi1`` and two exact profile derivatives.

    Evaluators receive two identically shaped arrays ``(r, z)``.  Supplying
    derivatives explicitly permits coefficient-space or analytic candidates
    to be reconstructed without using a cylindrical finite-difference result.
    For a nodal candidate, interpolation and derivative certification should be
    a separately audited adapter.
    """

    u1: ProfileEvaluator
    psi1: ProfileEvaluator
    dpsi1_dr: ProfileEvaluator
    dpsi1_dz: ProfileEvaluator


def _evaluate_profile(
    evaluator: ProfileEvaluator,
    radius: ScalarArray,
    axial: ScalarArray,
    *,
    name: str,
) -> ScalarArray:
    raw = np.asarray(evaluator(radius, axial), dtype=np.float64)
    try:
        values = np.broadcast_to(raw, radius.shape)
    except ValueError as error:
        raise ValueError(
            f"axisymmetric evaluator {name} is not broadcastable to {radius.shape}"
        ) from error
    if not np.all(np.isfinite(values)):
        raise ValueError(f"axisymmetric evaluator {name} returned non-finite values")
    return np.asarray(values, dtype=np.float64)


def reconstruct_axisymmetric_regularized(
    grid: UniformCartesianGrid,
    profile: RegularizedAxisymmetricProfile,
) -> VectorArray:
    r"""Reconstruct a regularized axisymmetric profile directly on ``(x,y,z)``.

    The nonsingular component map is

    .. math::

       u_x=-x\psi_{1,z}-y u_1,\quad
       u_y=-y\psi_{1,z}+x u_1,\quad
       u_z=2\psi_1+r\psi_{1,r}.

    It remains well-defined on the axis and does not divide by ``r``.
    """

    x, y, z = grid.mesh()
    radius = np.hypot(x, y)
    u1 = _evaluate_profile(profile.u1, radius, z, name="u1")
    psi1 = _evaluate_profile(profile.psi1, radius, z, name="psi1")
    psi1_r = _evaluate_profile(
        profile.dpsi1_dr,
        radius,
        z,
        name="dpsi1_dr",
    )
    psi1_z = _evaluate_profile(
        profile.dpsi1_dz,
        radius,
        z,
        name="dpsi1_dz",
    )
    return np.stack(
        (
            -x * psi1_z - y * u1,
            -y * psi1_z + x * u1,
            2.0 * psi1 + radius * psi1_r,
        ),
        axis=0,
    )


@dataclass(frozen=True)
class CartesianAuditReport:
    """Scalar acceptance summary for a reconstructed physical field."""

    divergence_rms: float
    divergence_max: float
    curl_defect_rms: float
    curl_defect_max: float
    divergence_rms_tolerance: float
    divergence_max_tolerance: float
    curl_rms_tolerance: float
    curl_max_tolerance: float
    interior_margin: int
    all_points_finite: bool

    @property
    def divergence_rms_passed(self) -> bool:
        return self.divergence_rms <= self.divergence_rms_tolerance

    @property
    def divergence_max_passed(self) -> bool:
        return self.divergence_max <= self.divergence_max_tolerance

    @property
    def divergence_passed(self) -> bool:
        return self.divergence_rms_passed and self.divergence_max_passed

    @property
    def curl_rms_passed(self) -> bool:
        return self.curl_defect_rms <= self.curl_rms_tolerance

    @property
    def curl_max_passed(self) -> bool:
        return self.curl_defect_max <= self.curl_max_tolerance

    @property
    def curl_passed(self) -> bool:
        return self.curl_rms_passed and self.curl_max_passed

    @property
    def passed(self) -> bool:
        return (
            self.all_points_finite
            and self.divergence_rms_passed
            and self.divergence_max_passed
            and self.curl_rms_passed
            and self.curl_max_passed
        )


def audit_cartesian_reconstruction(
    grid: UniformCartesianGrid,
    *,
    velocity: npt.ArrayLike,
    expected_vorticity: npt.ArrayLike,
    divergence_rms_tolerance: float,
    divergence_max_tolerance: float,
    curl_rms_tolerance: float,
    curl_max_tolerance: float,
    interior_margin: int = 0,
) -> CartesianAuditReport:
    """Audit divergence and full curl against an independent vorticity oracle.

    Finiteness is checked over *all* points before norms are computed.
    ``interior_margin`` may then omit a fixed number of points on each
    *non-periodic* side from the reported RMS and maximum norms, separating
    boundary-closure error from interior consistency.  Periodic axes retain
    every point, including their wrap seam.  The margin never hides a NaN or
    infinity because whole-array finiteness is checked first.

    RMS and maximum tolerances are distinct required gates.  This prevents a
    localized defect from being diluted by the number of grid points.
    """

    state = grid.validate_vector(velocity, name="velocity")
    target = grid.validate_vector(expected_vorticity, name="expected_vorticity")
    divergence_rms_limit = _finite_positive(
        divergence_rms_tolerance,
        name="divergence_rms_tolerance",
    )
    divergence_max_limit = _finite_positive(
        divergence_max_tolerance,
        name="divergence_max_tolerance",
    )
    curl_rms_limit = _finite_positive(
        curl_rms_tolerance,
        name="curl_rms_tolerance",
    )
    curl_max_limit = _finite_positive(
        curl_max_tolerance,
        name="curl_max_tolerance",
    )
    if (
        isinstance(interior_margin, bool)
        or not isinstance(interior_margin, (int, np.integer))
        or int(interior_margin) < 0
    ):
        raise ValueError("interior_margin must be a nonnegative integer")
    margin = int(interior_margin)
    if margin and any(
        size <= 2 * margin
        for size, is_periodic in zip(grid.shape, grid.periodic)
        if not is_periodic
    ):
        raise ValueError("interior_margin leaves an empty norm region")
    divergence = cartesian_divergence(grid, state)
    curl_defect = cartesian_curl(grid, state) - target
    region = tuple(
        slice(None)
        if is_periodic or margin == 0
        else slice(margin, -margin)
        for is_periodic in grid.periodic
    )
    vector_region = (slice(None), *region)
    norm_divergence = divergence[region]
    norm_curl_defect = curl_defect[vector_region]
    return CartesianAuditReport(
        divergence_rms=float(
            np.sqrt(np.mean(norm_divergence * norm_divergence))
        ),
        divergence_max=float(np.max(np.abs(norm_divergence))),
        curl_defect_rms=float(
            np.sqrt(np.mean(norm_curl_defect * norm_curl_defect))
        ),
        curl_defect_max=float(np.max(np.abs(norm_curl_defect))),
        divergence_rms_tolerance=divergence_rms_limit,
        divergence_max_tolerance=divergence_max_limit,
        curl_rms_tolerance=curl_rms_limit,
        curl_max_tolerance=curl_max_limit,
        interior_margin=margin,
        all_points_finite=True,
    )
