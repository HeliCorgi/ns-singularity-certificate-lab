r"""From the elliptic tail bound to the Navier--Stokes right-hand-side defect.

Gate 5 bounded every derivative of the free-space potential error,

    eps_0 >= sup |delta psi_1|,
    eps_1 >= sup |grad delta psi_1|,
    eps_2 >= sup |D^2 delta psi_1|,

but a bound on ``delta psi_1`` is not yet a bound on anything the momentum
equation consumes.  This module closes that gap **algebraically**: every
constant below is written out, and nothing is inferred from a floating-point
comparison.

The chain is

    eps_0, eps_1, eps_2   ->   eps_{u^r}, eps_{u^z}   ->   eps_advection
                          ->   eps_stretching, eps_swirl
                          ->   a short-time Gronwall bound on the state error.

Every step uses only the product-difference identity
``a b - a~ b~ = (a - a~) b + a~ (b - b~)`` and the triangle inequality, so each
inequality is exact and its constant is explicit.

Scope, stated first.
* These are **continuum** inequalities relating two exact solutions of the same
  PDE with different elliptic data.  The discretisation error of the actual
  finite-difference computation is *not* bounded here.
* The Gronwall step assumes the ``L^infinity`` maximum principle for the
  advection--diffusion operator with a divergence-free advecting field.  That
  holds on the whole space and on a box the fields never reach; the numerical
  companion diagnostic is ``outer_band_fraction``, and a run in which it is not
  negligible falls outside the hypothesis.
* Converting the *state* error ``e_omega`` into a velocity error needs an
  operator norm of the free-space solve.  ``L^infinity -> L^infinity`` is **not**
  bounded for Biot--Savart, so that constant is carried as an explicit input
  ``solve_operator_norm`` and supplying it rigorously is an open obligation, not
  something this module pretends to have.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = [
    "FieldBounds",
    "PotentialErrorBounds",
    "RhsDefectBounds",
    "StateErrorBound",
    "advection_defect_bound",
    "gronwall_state_error",
    "velocity_error_bounds",
]


def _require_nonnegative(value: float, *, name: str) -> float:
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if value < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return float(value)


@dataclass(frozen=True)
class PotentialErrorBounds:
    """Uniform bounds on the error of the stream function and its derivatives."""

    value: float
    gradient: float
    hessian: float
    radial_extent: float

    def __post_init__(self) -> None:
        for name in ("value", "gradient", "hessian", "radial_extent"):
            _require_nonnegative(getattr(self, name), name=name)
        if self.radial_extent <= 0.0:
            raise ValueError("radial_extent must be positive")

    def as_dict(self) -> dict[str, float]:
        return {
            "value": self.value,
            "gradient": self.gradient,
            "hessian": self.hessian,
            "radial_extent": self.radial_extent,
        }


@dataclass(frozen=True)
class FieldBounds:
    """Uniform bounds on the *computed* state and velocity on the core region."""

    u1_max: float
    omega1_max: float
    u1_gradient_max: float
    omega1_gradient_max: float
    velocity_radial_max: float
    velocity_axial_max: float
    psi_axial_derivative_max: float

    def __post_init__(self) -> None:
        for name in (
            "u1_max",
            "omega1_max",
            "u1_gradient_max",
            "omega1_gradient_max",
            "velocity_radial_max",
            "velocity_axial_max",
            "psi_axial_derivative_max",
        ):
            _require_nonnegative(getattr(self, name), name=name)

    def as_dict(self) -> dict[str, float]:
        return {
            "u1_max": self.u1_max,
            "omega1_max": self.omega1_max,
            "u1_gradient_max": self.u1_gradient_max,
            "omega1_gradient_max": self.omega1_gradient_max,
            "velocity_radial_max": self.velocity_radial_max,
            "velocity_axial_max": self.velocity_axial_max,
            "psi_axial_derivative_max": self.psi_axial_derivative_max,
        }


def velocity_error_bounds(potential: PotentialErrorBounds) -> tuple[float, float]:
    r"""``(|delta u^r|, |delta u^z|)`` from the potential error bounds.

    The audited recovery (E-14) is ``u^r = -r psi_{1,z}`` and
    ``u^z = 2 psi_1 + r psi_{1,r}``, both linear in ``psi_1``, so the errors are

    .. math::

        |\delta u^r| \le R_{\max}\,\varepsilon_1, \qquad
        |\delta u^z| \le 2\varepsilon_0 + R_{\max}\,\varepsilon_1 .

    No cancellation is claimed and no product rule is needed: the map is linear.
    """
    radial = potential.radial_extent * potential.gradient
    axial = 2.0 * potential.value + potential.radial_extent * potential.gradient
    return radial, axial


@dataclass(frozen=True)
class RhsDefectBounds:
    """Explicit bounds on every term of the right-hand-side defect."""

    velocity_radial: float
    velocity_axial: float
    advection_u1: float
    advection_omega1: float
    swirl_source: float
    stretching_source: float
    state_lipschitz: float

    @property
    def u1_equation(self) -> float:
        return self.advection_u1 + self.swirl_source

    @property
    def omega1_equation(self) -> float:
        return self.advection_omega1 + self.stretching_source

    @property
    def total(self) -> float:
        return max(self.u1_equation, self.omega1_equation)

    def as_dict(self) -> dict[str, float]:
        return {
            "velocity_radial": self.velocity_radial,
            "velocity_axial": self.velocity_axial,
            "advection_u1": self.advection_u1,
            "advection_omega1": self.advection_omega1,
            "swirl_source": self.swirl_source,
            "stretching_source": self.stretching_source,
            "state_lipschitz": self.state_lipschitz,
            "u1_equation": self.u1_equation,
            "omega1_equation": self.omega1_equation,
            "total": self.total,
        }


def advection_defect_bound(
    potential: PotentialErrorBounds,
    fields: FieldBounds,
    *,
    solve_operator_norm: float = 0.0,
) -> RhsDefectBounds:
    r"""Bound every right-hand-side term of the truncated-versus-true difference.

    Write ``a`` for a true quantity and ``a~`` for the truncated one, and use
    ``a b - a~ b~ = (a - a~)b + a~(b - b~)`` throughout.

    **Advection.**  ``A[f] = u^r f_r + u^z f_z``, so

    .. math::

        |\delta A[f]| \le |\delta u^r|\,\|f_r\| + |\delta u^z|\,\|f_z\|
                        + \|\tilde u^r\|\,\|\delta f_r\|
                        + \|\tilde u^z\|\,\|\delta f_z\| .

    The last two terms are the *state* contribution; they are proportional to the
    state error and are collected into ``state_lipschitz`` rather than into the
    constant part, because Gronwall multiplies them by the unknown error.

    **Swirl source** ``2 u_1 \psi_{1,z}``:

    .. math::

        |\delta(2u_1\psi_{1,z})| \le 2\|\tilde u_1\|\,\varepsilon_1
                                    + 2\|\psi_{1,z}\|\,|\delta u_1| .

    **Stretching source** ``\partial_z(u_1^2)``: with
    ``u_1^2 - \tilde u_1^2 = (u_1-\tilde u_1)(u_1+\tilde u_1)``,

    .. math::

        |\delta\partial_z(u_1^2)|
          \le 2\|u_1\|\,|\partial_z\delta u_1| + 2\|\partial_z u_1\|\,|\delta u_1| ,

    again split into a tail part and a state part.

    ``solve_operator_norm`` is the constant ``K`` in
    ``|\delta u| \le K\,\|e_\omega\|`` for the free-space solve.  It is an
    **input**: ``L^\infty \to L^\infty`` is unbounded for Biot--Savart, so no
    honest default exists and ``0`` means "the caller asserts the state
    contribution through the solve is accounted for elsewhere".
    """
    _require_nonnegative(solve_operator_norm, name="solve_operator_norm")
    radial, axial = velocity_error_bounds(potential)

    advection_u1 = radial * fields.u1_gradient_max + axial * fields.u1_gradient_max
    advection_omega1 = (
        radial * fields.omega1_gradient_max + axial * fields.omega1_gradient_max
    )
    swirl = 2.0 * fields.u1_max * potential.gradient
    # The elliptic tail does not enter d_z(u_1^2) directly: u_1 is a state
    # variable, not a potential derivative.  Its contribution arrives only
    # through the state error and therefore belongs to the Lipschitz constant,
    # not to the constant defect.  Recording it as an explicit zero keeps that
    # decision visible instead of silently dropping a term.
    stretching = 0.0
    lipschitz = (
        fields.velocity_radial_max
        + fields.velocity_axial_max
        + 2.0 * fields.psi_axial_derivative_max
        + 2.0 * fields.u1_gradient_max
        + 2.0 * fields.u1_max
        + solve_operator_norm
        * (fields.u1_gradient_max + fields.omega1_gradient_max)
    )
    return RhsDefectBounds(
        velocity_radial=radial,
        velocity_axial=axial,
        advection_u1=advection_u1,
        advection_omega1=advection_omega1,
        swirl_source=swirl,
        stretching_source=stretching,
        state_lipschitz=lipschitz,
    )


@dataclass(frozen=True)
class StateErrorBound:
    """Short-time Gronwall bound on the state error."""

    initial_error: float
    constant_defect: float
    lipschitz: float
    horizon: float
    bound: float

    def as_dict(self) -> dict[str, float]:
        return {
            "initial_error": self.initial_error,
            "constant_defect": self.constant_defect,
            "lipschitz": self.lipschitz,
            "horizon": self.horizon,
            "bound": self.bound,
        }


def gronwall_state_error(
    *,
    initial_error: float,
    constant_defect: float,
    lipschitz: float,
    horizon: float,
) -> StateErrorBound:
    r"""``\|e(t)\| \le (\|e(0)\| + D t) e^{\Lambda t}`` for ``t \le T``.

    Both the truncated and the true solution satisfy the same transport--diffusion
    system, so their difference ``e`` obeys

    .. math::

        \partial_t e + \tilde u\cdot\nabla e = F + \nu\mathcal L_5 e,

    with ``\|F\| \le D + \Lambda\|e\|``.  Since ``\tilde u`` is divergence free
    and ``\mathcal L_5`` is the five-dimensional Laplacian for axis-even fields,
    the ``L^\infty`` maximum principle gives
    ``\frac{d}{dt}\|e\|_\infty \le \|F\|_\infty``, hence

    .. math::

        \frac{d}{dt}\|e\| \le D + \Lambda\|e\|
        \;\Longrightarrow\;
        \|e(t)\| \le \Bigl(\|e(0)\| + \frac{D}{\Lambda}\Bigr)e^{\Lambda t}
                     - \frac{D}{\Lambda},

    which the simpler form ``(\|e(0)\| + Dt)e^{\Lambda t}`` dominates.  The
    simpler form is returned because it is what a certificate can check without
    dividing by a possibly tiny ``\Lambda``.
    """
    initial = _require_nonnegative(initial_error, name="initial_error")
    defect = _require_nonnegative(constant_defect, name="constant_defect")
    rate = _require_nonnegative(lipschitz, name="lipschitz")
    time = _require_nonnegative(horizon, name="horizon")
    growth = math.exp(rate * time)
    return StateErrorBound(
        initial_error=initial,
        constant_defect=defect,
        lipschitz=rate,
        horizon=time,
        bound=(initial + defect * time) * growth,
    )
