r"""Frozen-coefficient von Neumann audit of the explicit time steppers.

Why this module exists (P0-A)
-----------------------------
The production integrator :mod:`ns_certificate_lab.nonlinear_cylinder` advances
``(u_1, \omega_1)`` with explicit Heun/RK2 and discretizes every spatial
derivative with a centered difference.  For the constant-coefficient scalar
model problem

.. math::

   q_t + c\,q_z = 0,
   \qquad
   (q_z)_j \approx \frac{q_{j+1}-q_{j-1}}{2\Delta z},

the semi-discrete Fourier symbol is ``lambda = -i c sin(theta)/dz``, which is
*purely imaginary*.  Heun's linear stability function is
``G(z) = 1 + z + z^2/2``, and on the imaginary axis

.. math::

   |G(ia)|^2 = \left(1-\frac{a^2}{2}\right)^2 + a^2 = 1 + \frac{a^4}{4} > 1
   \qquad\text{for every } a \neq 0 .

So Heun with centered differences and *no viscosity* amplifies every resolved
mode at every nonzero CFL number: the explicit method has no imaginary-axis
stability interval at all.  Viscosity moves the symbol off the imaginary axis
and can restore ``|G| \le 1``, but only if the parabolic damping at the
worst-case wavenumber beats the quartic advective growth.  Whether it actually
does so at a given operating point is an arithmetic question, and this module
is what answers it, so that recorded snapshots can be audited after the fact
instead of being trusted.

Contrast with the other methods tabulated here: SSPRK3 and RK4 *do* have
imaginary-axis intervals (``|a| \le \sqrt{3}`` and ``|a| \le 2\sqrt{2}``
respectively), which is why they are listed -- they are what a fix would look
like.  Their stability functions are stated below and every boundary claim in
the tests is recomputed from :func:`scan_amplification` rather than asserted
from these numbers.

What this analysis does *not* cover
-----------------------------------
Von Neumann analysis is exact only for a linear, constant-coefficient,
*periodic* problem.  Every one of those three words fails somewhere in the
production solver, and the failures are listed here rather than buried:

* **Variable coefficients.**  ``u^r`` and ``u^z`` vary over the grid.  This
  module freezes them at their measured maxima, which is a model: it is
  neither an upper bound nor a lower bound for the variable-coefficient
  operator, because frozen-coefficient analysis discards the commutator terms
  entirely.
* **The axis row.**  ``\mathcal L_5 = \partial_{rr} + (3/r)\partial_r +
  \partial_{zz}`` degenerates at ``r = 0``; the audited discretization uses the
  coefficient-8 limit ``8(f_1-f_0)/\Delta r^2`` (E-17 and E-26b in
  ``docs/equation_audit.md``).  That row is a *boundary modification* of the
  stencil, not a translation-invariant one, so no Fourier mode diagonalizes it
  and nothing below says anything about it.  The frozen ``3\nu/r`` first
  derivative modelled by ``first_derivative_coefficient_r`` is the nearest
  honest proxy, evaluated at the worst interior radius ``r \approx \Delta r``;
  it is a proxy, not the axis row.
* **The wall rows.**  The no-slip trace ``u_1(r_{max},\cdot) = 0`` and the
  Thom-type wall vorticity are algebraic constraints on the outermost radial
  row (E-27, E-31).  They are boundary modifications for the same reason, and
  GKS-type boundary stability is out of scope here.
* **The elliptic coupling.**  Every Heun stage re-solves ``-\mathcal L_5\psi_1
  = \omega_1`` and recovers ``(u^r,u^z)``.  That is a nonlocal operator applied
  inside the stage; it is not represented at all by the advection--diffusion
  symbol below.
* **The zeroth-order terms.**  ``2u_1\psi_{1,z}`` and ``\partial_z(u_1^2)``
  are reaction/production terms.  They are omitted, which makes the model
  optimistic in that direction.

Consequently a *failing* audit means the run is **stability-unverified**: the
frozen-coefficient worst case is sufficient grounds for suspicion and for
re-running at a smaller step, and it is not a proof that the scheme blew up.
A *passing* audit is equally limited -- it verifies the model, not the solver.
Nothing in this module is evidence for or against singularity formation.

Sign conventions
----------------
The symbol is written so that the semi-discrete system is ``q_t = lambda q``:

.. math::

   \lambda = -i\left(\frac{c_r\sin\theta_r}{\Delta r}
                     +\frac{c_z\sin\theta_z}{\Delta z}\right)
             -\nu\left(\frac{2-2\cos\theta_r}{\Delta r^2}
                       +\frac{2-2\cos\theta_z}{\Delta z^2}\right)
             -i\,\kappa_r\frac{\sin\theta_r}{\Delta r},

with ``\kappa_r`` the frozen radial first-derivative coefficient.  The advective
and ``\kappa_r`` contributions are written with the *same* sign, so that they
add.  That is deliberate: ``u^r`` takes both signs over the meridional plane
while ``3\nu/r`` has a fixed sign, so both relative orientations occur
somewhere in the domain and the additive one is the conservative choice.
:func:`scan_amplification` therefore feeds the *magnitudes* of the advection
speeds into the symbol; :func:`advection_diffusion_symbol` is the lower-level
routine and uses whatever signs it is given.

All arithmetic is binary64 and is not outward rounded.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Mapping

import numpy as np
import numpy.typing as npt

__all__ = [
    "FrozenCoefficientSymbol",
    "METHODS",
    "advection_diffusion_symbol",
    "audit_snapshot",
    "reference_mode_amplification",
    "reference_propagate",
    "scan_amplification",
    "stability_polynomial",
    "theta_grid",
]

FloatArray = npt.NDArray[np.float64]
ComplexArray = npt.NDArray[np.complex128]

# For a *linear autonomous* problem q' = lambda q every classical explicit
# Runge-Kutta method of order p with s = p stages reduces to the degree-p
# Taylor polynomial of the exponential.  That is why SSPRK3 and RK4 appear here
# as degrees 3 and 4: the statement is specific to the linear model problem and
# says nothing about their nonlinear behaviour (SSP coefficients, for example,
# are invisible to it).  "euler" is listed under its own name because it is the
# predictor stage of Heun, and auditing the predictor separately matters: a
# stage can overflow even when the completed step does not.
_TAYLOR_DEGREE: Mapping[str, int] = {
    "euler": 1,
    "heun": 2,
    "ssprk3": 3,
    "rk4": 4,
}

METHODS: tuple[str, ...] = ("euler", "heun", "ssprk3", "rk4")

_PREDICTOR_STAGE: Mapping[str, str] = {"heun": "euler"}


def _require_method(method: object) -> str:
    """Return a validated method name or raise :class:`ValueError`."""

    if not isinstance(method, str):
        raise TypeError(f"method must be a string, got {type(method).__name__}")
    if method not in _TAYLOR_DEGREE:
        raise ValueError(
            f"unknown method {method!r}; supported methods are "
            f"{', '.join(sorted(_TAYLOR_DEGREE))}"
        )
    return method


_REAL_TYPES = (int, float, np.floating, np.integer)


def _as_real(value: object, *, name: str) -> float:
    """Return ``value`` as a float, rejecting bools and non-numbers.

    ``bool`` is excluded explicitly because it is a subclass of ``int`` and
    would otherwise be accepted as the step size ``1.0``.
    """

    if isinstance(value, bool) or not isinstance(value, _REAL_TYPES):
        raise TypeError(f"{name} must be a real number")
    return float(value)


def _require_positive(value: object, *, name: str) -> float:
    """Return ``value`` as a positive finite float or raise."""

    number = _as_real(value, name=name)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite, got {number!r}")
    return number


def _require_nonnegative(value: object, *, name: str) -> float:
    """Return ``value`` as a nonnegative finite float or raise."""

    number = _as_real(value, name=name)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(
            f"{name} must be nonnegative and finite, got {number!r}"
        )
    return number


def _require_finite(value: object, *, name: str) -> float:
    """Return ``value`` as a finite float of either sign or raise."""

    number = _as_real(value, name=name)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {number!r}")
    return number


def _require_count(value: object, *, name: str, minimum: int) -> int:
    """Return ``value`` as an integer at least ``minimum`` or raise."""

    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    count = int(value)
    if count < minimum:
        raise ValueError(f"{name} must be at least {minimum}, got {count}")
    return count


def stability_polynomial(method: str, z: complex | npt.ArrayLike) -> ComplexArray:
    r"""Return the linear stability function ``G(z)`` of an explicit method.

    ``G`` is defined by ``q_{n+1} = G(\Delta t\,\lambda) q_n`` for the scalar
    model problem ``q' = \lambda q``:

    ======== =========================================================
    method   ``G(z)``
    ======== =========================================================
    euler    ``1 + z``
    heun     ``1 + z + z^2/2``
    ssprk3   ``1 + z + z^2/2 + z^3/6``
    rk4      ``1 + z + z^2/2 + z^3/6 + z^4/24``
    ======== =========================================================

    The last two coincide with the degree-3 and degree-4 Taylor polynomials of
    ``exp`` only because the model problem is linear and autonomous; on a
    nonlinear problem the three-stage SSP method and classical RK4 differ from
    each other and from these polynomials.

    The value is evaluated by Horner's rule
    ``1 + z(1 + \tfrac{z}{2}(1 + \tfrac{z}{3}(\cdots)))``, which is the same
    polynomial exactly and is better conditioned than summing the monomials.

    Returns a ``complex128`` array (0-d for scalar input).  Non-finite input is
    rejected rather than propagated, because a NaN amplification factor would
    silently pass a ``<=`` acceptance test.
    """

    name = _require_method(method)
    values = np.asarray(z, dtype=np.complex128)
    if not np.all(np.isfinite(values)):
        raise ValueError("z must contain only finite values")

    degree = _TAYLOR_DEGREE[name]
    result = np.ones_like(values)
    for term in range(degree, 0, -1):
        result = 1.0 + values * result / float(term)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError(
            "stability polynomial evaluation produced a non-finite value"
        )
    return result


@dataclass(frozen=True)
class FrozenCoefficientSymbol:
    r"""Validated frozen-coefficient data for the centered advection--diffusion symbol.

    The fields are the coefficients of the *model* operator, not of the
    production operator; see the module docstring for the list of production
    features this model omits.
    """

    advection_r: float
    advection_z: float
    viscosity: float
    dr: float
    dz: float
    first_derivative_coefficient_r: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "advection_r", _require_finite(self.advection_r, name="advection_r")
        )
        object.__setattr__(
            self, "advection_z", _require_finite(self.advection_z, name="advection_z")
        )
        object.__setattr__(
            self, "viscosity", _require_nonnegative(self.viscosity, name="viscosity")
        )
        object.__setattr__(self, "dr", _require_positive(self.dr, name="dr"))
        object.__setattr__(self, "dz", _require_positive(self.dz, name="dz"))
        object.__setattr__(
            self,
            "first_derivative_coefficient_r",
            _require_finite(
                self.first_derivative_coefficient_r,
                name="first_derivative_coefficient_r",
            ),
        )

    def evaluate(
        self,
        theta_r: npt.ArrayLike,
        theta_z: npt.ArrayLike,
    ) -> ComplexArray:
        """Return ``lambda(theta_r, theta_z)`` with NumPy broadcasting."""

        angle_r = np.asarray(theta_r, dtype=np.float64)
        angle_z = np.asarray(theta_z, dtype=np.float64)
        if not (np.all(np.isfinite(angle_r)) and np.all(np.isfinite(angle_z))):
            raise ValueError("theta_r and theta_z must contain only finite values")

        sin_r = np.sin(angle_r)
        sin_z = np.sin(angle_z)
        imaginary = (
            self.advection_r * sin_r / self.dr
            + self.advection_z * sin_z / self.dz
            + self.first_derivative_coefficient_r * sin_r / self.dr
        )
        real = -self.viscosity * (
            (2.0 - 2.0 * np.cos(angle_r)) / self.dr**2
            + (2.0 - 2.0 * np.cos(angle_z)) / self.dz**2
        )
        return (real - 1j * imaginary).astype(np.complex128, copy=False)


def advection_diffusion_symbol(
    theta_r: npt.ArrayLike,
    theta_z: npt.ArrayLike,
    *,
    advection_r: float,
    advection_z: float,
    viscosity: float,
    dr: float,
    dz: float,
    first_derivative_coefficient_r: float = 0.0,
) -> ComplexArray:
    r"""Frozen-coefficient symbol of the centered advection--diffusion stencil.

    .. math::

       \lambda = -i\left(\frac{c_r\sin\theta_r}{\Delta r}
                         +\frac{c_z\sin\theta_z}{\Delta z}\right)
                 -\nu\left(\frac{2-2\cos\theta_r}{\Delta r^2}
                           +\frac{2-2\cos\theta_z}{\Delta z^2}\right)
                 -i\,\kappa_r\frac{\sin\theta_r}{\Delta r}.

    The first two groups are the exact symbols of ``-c\,\partial`` with the
    three-point centered first difference and of ``\nu\,\partial^2`` with the
    three-point centered second difference on a *periodic* uniform grid.

    ``first_derivative_coefficient_r`` (``\kappa_r``) models the frozen
    ``3\nu/r`` first-derivative part of ``\mathcal L_5 = \partial_{rr} +
    (3/r)\partial_r + \partial_{zz}``.  At the worst *interior* radius
    ``r \approx \Delta r`` it takes the value ``3\nu/\Delta r``, so the term
    contributes ``3\nu\sin\theta_r/\Delta r^2``.

    This is a periodic frozen-coefficient MODEL.  The actual axis row is the
    coefficient-8 limit ``8(f_1-f_0)/\Delta r^2`` of E-17/E-26b, and the
    outermost radial row carries the E-27/E-31 wall constraints.  Both are
    boundary modifications of the stencil; von Neumann analysis does not cover
    either of them, and no value returned here is a statement about them.
    """

    symbol = FrozenCoefficientSymbol(
        advection_r=advection_r,
        advection_z=advection_z,
        viscosity=viscosity,
        dr=dr,
        dz=dz,
        first_derivative_coefficient_r=first_derivative_coefficient_r,
    )
    return symbol.evaluate(theta_r, theta_z)


def theta_grid(n_theta: int = 721) -> FloatArray:
    r"""Return ``n_theta`` samples of ``[-\pi, \pi]``.

    The grid contains ``0`` and ``\pm\pi`` as *exact* samples.
    ``n_theta`` must be odd, which is what makes ``theta = 0`` (the constant
    mode) a sample point; ``theta = \pm\pi`` (the sawtooth mode, the shortest
    wave the grid carries) are the endpoints.  Both matter: the sawtooth is
    where the diffusive damping is largest and where a centered advection
    symbol vanishes, and the constant mode is where both vanish, so a scan that
    missed either could report a spuriously small maximum.

    The grid is built by mirroring ``numpy.linspace(0, pi, (n_theta+1)//2)``.
    That construction has a property the tests rely on: for ``m`` and ``4m-3``
    style refinements the coarse samples are *bitwise* a subset of the fine
    ones, because the fine step is the coarse step divided by a power of two
    and binary64 division by a power of two is exact.  The scan maximum is
    therefore monotone nondecreasing under such refinement, rather than merely
    approximately so.
    """

    count = _require_count(n_theta, name="n_theta", minimum=3)
    if count % 2 == 0:
        raise ValueError(
            "n_theta must be odd so that theta = 0 is a sample point, "
            f"got {count}"
        )
    half = np.linspace(0.0, np.pi, (count + 1) // 2, dtype=np.float64)
    grid = np.concatenate((-half[:0:-1], half))
    if grid.shape != (count,):
        raise ArithmeticError("theta grid construction produced the wrong length")
    # The three claims above are pinned here rather than trusted.
    if not (
        grid[0] == -np.pi
        and grid[-1] == np.pi
        and grid[(count - 1) // 2] == 0.0
    ):
        raise ArithmeticError(
            "theta grid does not contain -pi, 0 and pi as exact samples"
        )
    return grid


def _normalize_methods(methods: Iterable[str] | str) -> tuple[str, ...]:
    """Validate a method collection, preserving order and dropping duplicates."""

    if isinstance(methods, str):
        raise TypeError(
            "methods must be a collection of method names, not a single string"
        )
    try:
        candidates = list(methods)
    except TypeError as error:  # pragma: no cover - defensive
        raise TypeError("methods must be iterable") from error
    if not candidates:
        raise ValueError("methods must contain at least one method name")
    ordered: list[str] = []
    for candidate in candidates:
        name = _require_method(candidate)
        if name not in ordered:
            ordered.append(name)
    return tuple(ordered)


def scan_amplification(
    *,
    methods: Iterable[str],
    dt: float,
    dr: float,
    dz: float,
    advection_r: float,
    advection_z: float,
    viscosity: float,
    n_theta: int = 721,
    include_radial_first_derivative: bool = True,
) -> dict[str, Any]:
    r"""Maximize ``|G(\Delta t\,\lambda)|`` over the discrete wavenumber square.

    ``theta_r`` and ``theta_z`` each run over ``theta_grid(n_theta)``, which
    always contains ``0`` and ``\pm\pi`` exactly.  The scan is a maximum over a
    *finite sample* of a continuous function; it is therefore a lower bound for
    the true supremum over ``[-\pi,\pi]^2``.  Refining ``n_theta`` can only
    raise it (see :func:`theta_grid`), and the tests bound the gap on a case
    whose maximizer is interior.

    The magnitudes ``|advection_r|`` and ``|advection_z|`` are used, together
    with a nonnegative ``kappa_r``, so that the advective and ``3\nu/r``
    imaginary contributions add; see the module docstring for why that
    orientation is the conservative one.

    When ``include_radial_first_derivative`` is true, ``kappa_r`` is set to
    ``3*viscosity/dr``, the frozen worst case at the smallest interior radius.

    Returns
    -------
    dict
        ``"methods"`` maps each requested method name to
        ``{"max_amplification", "argmax_theta_r", "argmax_theta_z"}``, plus
        ``"predictor_max_amplification"`` and the predictor argmax for
        ``"heun"`` (its explicit-Euler predictor stage).  The top level also
        carries ``"advective_cfl_r" = dt|c_r|/dr``,
        ``"advective_cfl_z" = dt|c_z|/dz`` and
        ``"viscous_number" = 4*nu*dt/min(dr,dz)^2``, the grid data, and a
        scope string.  Everything is a plain JSON-representable value.
    """

    names = _normalize_methods(methods)
    step = _require_positive(dt, name="dt")
    spacing_r = _require_positive(dr, name="dr")
    spacing_z = _require_positive(dz, name="dz")
    nu = _require_nonnegative(viscosity, name="viscosity")
    speed_r = abs(_require_finite(advection_r, name="advection_r"))
    speed_z = abs(_require_finite(advection_z, name="advection_z"))
    if not isinstance(include_radial_first_derivative, bool):
        raise TypeError("include_radial_first_derivative must be a bool")

    kappa_r = 3.0 * nu / spacing_r if include_radial_first_derivative else 0.0
    symbol = FrozenCoefficientSymbol(
        advection_r=speed_r,
        advection_z=speed_z,
        viscosity=nu,
        dr=spacing_r,
        dz=spacing_z,
        first_derivative_coefficient_r=kappa_r,
    )

    angles = theta_grid(n_theta)
    scaled = step * symbol.evaluate(angles[:, None], angles[None, :])

    def _summarize(name: str) -> tuple[float, float, float]:
        magnitude = np.abs(stability_polynomial(name, scaled))
        flat = int(np.argmax(magnitude))
        row, column = np.unravel_index(flat, magnitude.shape)
        return (
            float(magnitude[row, column]),
            float(angles[row]),
            float(angles[column]),
        )

    per_method: dict[str, dict[str, float]] = {}
    for name in names:
        peak, at_r, at_z = _summarize(name)
        entry: dict[str, float] = {
            "max_amplification": peak,
            "argmax_theta_r": at_r,
            "argmax_theta_z": at_z,
        }
        predictor = _PREDICTOR_STAGE.get(name)
        if predictor is not None:
            stage_peak, stage_r, stage_z = _summarize(predictor)
            entry["predictor_stage"] = predictor
            entry["predictor_max_amplification"] = stage_peak
            entry["predictor_argmax_theta_r"] = stage_r
            entry["predictor_argmax_theta_z"] = stage_z
        per_method[name] = entry

    return {
        "methods": per_method,
        "dt": step,
        "dr": spacing_r,
        "dz": spacing_z,
        "advection_r": speed_r,
        "advection_z": speed_z,
        "viscosity": nu,
        "n_theta": int(n_theta),
        "include_radial_first_derivative": include_radial_first_derivative,
        "first_derivative_coefficient_r": float(kappa_r),
        "advective_cfl_r": step * speed_r / spacing_r,
        "advective_cfl_z": step * speed_z / spacing_z,
        "viscous_number": 4.0 * nu * step / min(spacing_r, spacing_z) ** 2,
        "model": (
            "frozen-coefficient periodic advection-diffusion with centered "
            "first and second differences; the axis row (E-17 coefficient-8 "
            "limit), the wall rows (E-27/E-31), the elliptic coupling and the "
            "zeroth-order production terms are NOT modelled"
        ),
        "scope": (
            "maximum over a finite wavenumber sample, hence a lower bound for "
            "the true frozen-coefficient supremum; not a bound on the "
            "variable-coefficient solver"
        ),
    }


def audit_snapshot(
    *,
    max_abs_u_r: float,
    max_abs_u_z: float,
    dr: float,
    dz: float,
    dt: float,
    viscosity: float,
    tolerance: float = 1.0e-12,
    methods: Iterable[str] = ("heun",),
) -> dict[str, Any]:
    r"""Audit one recorded snapshot against the frozen-coefficient worst case.

    ``max_abs_u_r`` and ``max_abs_u_z`` are the measured meridional velocity
    maxima of the snapshot; they are frozen over the whole grid, which is the
    worst case for the *model* and is not an upper bound for the actual
    variable-coefficient operator.

    A method ``"passes"`` when its worst-case one-step amplification is at most
    ``1 + tolerance``.

    Interpretation, stated exactly
    ------------------------------
    A run whose worst accepted step fails this audit is
    **stability-unverified**, not **unstable**.  The frozen-coefficient worst
    case is sufficient grounds for suspicion -- and for re-running the segment
    at a smaller step -- but it is *not* a proof that the scheme's solution
    grows, because the model discards the variable-coefficient commutators, the
    axis and wall boundary rows, the elliptic coupling and the zeroth-order
    terms.  Conversely a pass verifies the model and not the solver.

    Returns
    -------
    dict
        The :func:`scan_amplification` result, with ``"passes"`` added to each
        method entry, plus top-level ``"tolerance"``, ``"passes"`` (the
        conjunction over the requested methods) and ``"verdict"``, which is
        either ``"stability-verified-in-model"`` or ``"stability-unverified"``.
    """

    speed_r = _require_nonnegative(max_abs_u_r, name="max_abs_u_r")
    speed_z = _require_nonnegative(max_abs_u_z, name="max_abs_u_z")
    accept = _require_nonnegative(tolerance, name="tolerance")

    result = scan_amplification(
        methods=methods,
        dt=dt,
        dr=dr,
        dz=dz,
        advection_r=speed_r,
        advection_z=speed_z,
        viscosity=viscosity,
    )

    threshold = 1.0 + accept
    all_pass = True
    for entry in result["methods"].values():
        passes = bool(entry["max_amplification"] <= threshold)
        entry["passes"] = passes
        all_pass = all_pass and passes

    result["max_abs_u_r"] = speed_r
    result["max_abs_u_z"] = speed_z
    result["tolerance"] = accept
    result["passes"] = all_pass
    result["verdict"] = (
        "stability-verified-in-model" if all_pass else "stability-unverified"
    )
    result["verdict_meaning"] = (
        "a failing verdict means the frozen-coefficient worst case was not "
        "verified; it is grounds for suspicion and for a smaller step, not a "
        "proof that the scheme is unstable"
    )
    return result


# ---------------------------------------------------------------------------
# Discrete reference propagators.  VALIDATION ONLY.
#
# These exist so that the symbol path above can be checked against an actual
# time integration.  Two independent paths are provided and neither calls the
# other:
#
#   * :func:`reference_mode_amplification` runs the Butcher tableau on a single
#     complex mode amplitude driven by the discrete symbol;
#   * :func:`reference_propagate` runs the same tableau on a real (or complex)
#     array with ``numpy.roll`` centered differences and never forms a symbol.
#
# Neither calls :func:`stability_polynomial`, so a wrong stability polynomial
# cannot make its own verification pass.  They are not production integrators
# and know nothing about the cylinder geometry.
# ---------------------------------------------------------------------------


def _advance(
    method: str,
    values: Any,
    operator: Any,
    dt: float,
) -> Any:
    """Apply one explicit step of ``method`` to ``values`` under ``operator``.

    ``operator`` is any callable implementing the semi-discrete right-hand
    side; it is used identically by the scalar-mode and the array propagators,
    which is what makes those two paths share their time stepping and differ
    only in their spatial operator.

    The tableaux are written in their standard forms rather than as
    polynomials, so that reproducing :func:`stability_polynomial` is a genuine
    check.  SSPRK3 uses the Shu--Osher convex-combination form.
    """

    if method == "euler":
        return values + dt * operator(values)
    if method == "heun":
        k1 = operator(values)
        k2 = operator(values + dt * k1)
        return values + 0.5 * dt * (k1 + k2)
    if method == "ssprk3":
        stage1 = values + dt * operator(values)
        stage2 = 0.75 * values + 0.25 * (stage1 + dt * operator(stage1))
        return (1.0 / 3.0) * values + (2.0 / 3.0) * (
            stage2 + dt * operator(stage2)
        )
    if method == "rk4":
        k1 = operator(values)
        k2 = operator(values + 0.5 * dt * k1)
        k3 = operator(values + 0.5 * dt * k2)
        k4 = operator(values + dt * k3)
        return values + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    raise ValueError(f"unknown method {method!r}")  # pragma: no cover - guarded


def reference_mode_amplification(
    method: str,
    *,
    wavenumber_index: int,
    n_points: int,
    dt: float,
    advection: float,
    viscosity: float,
    n_steps: int,
    period: float = 1.0,
) -> complex:
    r"""Return the accumulated amplitude of one Fourier mode after ``n_steps``.

    The mode is ``e^{2\pi i k x}`` sampled on a periodic grid of ``n_points``
    covering ``[0, period)``, so ``\Delta x = period/n_points`` and
    ``\theta = 2\pi k/n_{points}``.  Its amplitude obeys the scalar ODE
    ``a' = \lambda a`` with ``\lambda`` the discrete symbol of the centered
    stencil, obtained from :func:`advection_diffusion_symbol` at
    ``theta_r = 0`` (which annihilates both radial terms exactly).

    The starting amplitude is ``1``, so the returned complex number *is* the
    accumulated amplification factor; the per-step factor is its ``n_steps``-th
    root.  The Butcher tableau is applied stage by stage, so this routine is
    independent of :func:`stability_polynomial`.

    VALIDATION ONLY.  This is not a solver.
    """

    name = _require_method(method)
    count = _require_count(n_points, name="n_points", minimum=3)
    index = _require_count(wavenumber_index, name="wavenumber_index", minimum=0)
    if index > count // 2:
        raise ValueError(
            f"wavenumber_index must not exceed n_points//2 = {count // 2}, "
            f"got {index}"
        )
    step = _require_positive(dt, name="dt")
    steps = _require_count(n_steps, name="n_steps", minimum=1)
    length = _require_positive(period, name="period")
    speed = _require_finite(advection, name="advection")
    nu = _require_nonnegative(viscosity, name="viscosity")

    spacing = length / count
    theta = 2.0 * np.pi * index / count
    symbol = complex(
        advection_diffusion_symbol(
            0.0,
            theta,
            advection_r=0.0,
            advection_z=speed,
            viscosity=nu,
            dr=spacing,
            dz=spacing,
        )
    )

    amplitude = 1.0 + 0.0j
    for _ in range(steps):
        amplitude = _advance(name, amplitude, lambda value: symbol * value, step)
        if not math.isfinite(amplitude.real) or not math.isfinite(amplitude.imag):
            raise FloatingPointError(
                "reference mode amplification became non-finite; the requested "
                "step is violently unstable for this method"
            )
    return complex(amplitude)


def reference_propagate(
    method: str,
    initial_values: npt.ArrayLike,
    *,
    dt: float,
    dx: float,
    advection: float,
    viscosity: float,
    n_steps: int,
) -> npt.NDArray[Any]:
    r"""Advance a periodic 1D array under ``q_t + c q_x = \nu q_{xx}``.

    Space is discretized with the same centered differences the production
    solver uses, applied with :func:`numpy.roll`:

    .. math::

       (q_x)_j = \frac{q_{j+1}-q_{j-1}}{2\Delta x},
       \qquad
       (q_{xx})_j = \frac{q_{j+1}-2q_j+q_{j-1}}{\Delta x^2}.

    No Fourier transform and no symbol appears anywhere in this function: it is
    the array path against which the symbol path is checked.  Real input is
    kept real and complex input is kept complex, so a single mode can be
    propagated exactly as a complex array.

    VALIDATION ONLY.  This is not a solver.
    """

    name = _require_method(method)
    step = _require_positive(dt, name="dt")
    spacing = _require_positive(dx, name="dx")
    steps = _require_count(n_steps, name="n_steps", minimum=1)
    speed = _require_finite(advection, name="advection")
    nu = _require_nonnegative(viscosity, name="viscosity")

    raw = np.asarray(initial_values)
    dtype = np.complex128 if np.iscomplexobj(raw) else np.float64
    values = np.array(raw, dtype=dtype, copy=True)
    if values.ndim != 1:
        raise ValueError(
            f"initial_values must be one-dimensional, got {values.ndim} dimensions"
        )
    if values.size < 3:
        raise ValueError("initial_values must have at least three points")
    if not np.all(np.isfinite(values)):
        raise ValueError("initial_values must contain only finite values")

    def operator(field: npt.NDArray[Any]) -> npt.NDArray[Any]:
        forward = np.roll(field, -1)
        backward = np.roll(field, 1)
        return (
            -speed * (forward - backward) / (2.0 * spacing)
            + nu * (forward - 2.0 * field + backward) / spacing**2
        )

    # A violently unstable request overflows to infinity, which NumPy reports
    # as a warning by default.  It is promoted to an exception here so that no
    # caller can mistake a saturated array for a computed one.
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        try:
            for _ in range(steps):
                values = _advance(name, values, operator, step)
                if not np.all(np.isfinite(values)):
                    raise FloatingPointError(
                        "reference propagation became non-finite; the "
                        "requested step is violently unstable for this method"
                    )
        except FloatingPointError as error:
            raise FloatingPointError(
                "reference propagation became non-finite; the requested step "
                f"is violently unstable for this method ({error})"
            ) from error
    return values
