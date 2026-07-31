r"""Gaussian--Hermite mixed initial data: polynomial-times-Gaussian generators.

The chi-bump generators of :mod:`ns_certificate_lab.mixed_initial_data` buy
compact support at a real price: ``\chi`` has enormous high derivatives near
the edge of its support, so the finite-difference cross-check of the analytic
gradient only enters its asymptotic range around ``nr = 97`` and converges at
an order of about ``1.55`` even there (see
``test_analytic_and_difference_gradients_agree_under_refinement``).  This
module is the opposite trade.  Nothing here is compactly supported, but the
generators are *Hermite-type* functions

.. math::

   u_1 = P(s, z)\,e^{-\alpha s - \beta z^2}, \qquad
   \psi_1 = Q(s, z)\,e^{-\alpha s - \beta z^2}, \qquad s = r^2 ,

with rational ``\alpha, \beta > 0`` and ``P, Q`` polynomials with rational
coefficients.  The velocity is recovered through the same audited map as the
chi families,

.. math::

   u^\theta = r\,u_1,\qquad u^r = -r\,\partial_z\psi_1,\qquad
   u^z = 2\psi_1 + 2s\,\partial_s\psi_1 ,

which is divergence free as an algebraic identity for **any** stream generator
(the cancellation is spelled out in :class:`~ns_certificate_lab.l3_generation.
MixedField`), and everything is a function of ``s = r^2``, never of ``r``, so
the Cartesian field is ``C^\infty`` across the axis by exactly the argument in
the :mod:`~ns_certificate_lab.mixed_initial_data` docstring.

Why this family earns its keep
------------------------------
The class is **closed under the derivatives the gradient needs**:

.. math::

   \partial_s\bigl(s^mz^n E\bigr) = (m\,s^{m-1} - \alpha s^m)\,z^n E, \qquad
   \partial_z\bigl(s^mz^n E\bigr) = (n\,z^{n-1} - 2\beta z^{n+1})\,s^m E ,

with ``E = e^{-\alpha s - \beta z^2}``.  Every partial derivative of a
generator is therefore *another* polynomial of the same kind, computed by an
exact recursion on rational coefficients and only then evaluated in floating
point.  The exact velocity gradient costs a few dictionary operations, no
automatic differentiation, and the derivatives stay mild: the finite-difference
cross-check converges at clean second order from ``nr = 33``, in sharp contrast
to the chi bump.

Admissibility and the tail
--------------------------
Gaussian--Hermite data are Schwartz class: ``C^\infty`` with all derivatives
decaying faster than any polynomial.  They are legitimate smooth, decaying
initial data in the sense the repository's guard requires, so
``clay_admissible`` is ``True`` — these are not the finite-``C^k`` spline
surrogates that :func:`~ns_certificate_lab.l3_optimizer.require_clay_admissible`
exists to refuse.  What is lost against the chi families is *compact* support,
and everything compactness used to certify for free — the validity of the
exterior pressure representation, the outer gap of the box — is instead
controlled by an explicit, rigorous tail bound:
:meth:`GaussianMixedFamily.tail_bound` returns a certified upper bound on
``|u|`` outside any sphere, and :meth:`GaussianMixedFamily.support_radius`
inverts it.  The bound is a theorem about the datum, not a sampling of it; its
derivation is in the :meth:`~GaussianMixedFamily.tail_bound` docstring.

The axial-parity selection rule applies verbatim: the pressure contribution to
the ``L^3`` generation rate survives integration over a ``z``-symmetric domain
only when the stream generator is **odd** in ``z``, which for a Hermite
generator means ``Q`` contains only odd powers of ``z`` — a property read off
the exponent dictionary rather than sampled.  All stream generators of the
preregistered basis below are odd.

Nothing in this module bears on the Clay problem; the family is search
infrastructure for the ``\mathrm{Re}_{\rm crit}`` objective of
:mod:`ns_certificate_lab.l3_optimizer`, and no result obtained with it should
be described as progress on, or evidence about, global regularity.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
import math
from typing import Mapping

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "GaussianPolynomial",
    "GaussianGenerator",
    "GaussianMixedFamily",
    "gaussian_search_basis",
]


def _power(base: float, exponent: float) -> float:
    """``base**exponent`` with the convention ``x**0 = 1`` even at ``x = 0``.

    The tail bound needs ``0^0 = 1`` in three places (``a``, ``b/2`` and their
    sum), and Python's ``0.0**0.0`` already returns ``1.0``; the helper exists
    so the convention is *stated* rather than relied on silently.
    """
    return 1.0 if exponent == 0.0 else float(base) ** exponent


class GaussianPolynomial:
    r"""A polynomial in ``(s, z)`` with exact rational coefficients.

    Stored as ``{(m, n): coefficient}`` for the monomial ``s^m z^n``.  The
    point of keeping coefficients as :class:`~fractions.Fraction` is that the
    derivative recursions

    .. math::

       \partial_s\colon (m, n) \mapsto m\,(m{-}1, n) - \alpha\,(m, n), \qquad
       \partial_z\colon (m, n) \mapsto n\,(m, n{-}1) - 2\beta\,(m, n{+}1)

    are then *exact*: the certified gradient of a family is produced by
    integer and rational arithmetic and floating point enters only at the final
    evaluation, so there is no accumulation of rounding across repeated
    differentiation.  Instances are treated as immutable; every operation
    returns a new polynomial and zero coefficients are dropped on construction.
    """

    __slots__ = ("_coefficients",)

    def __init__(
        self, coefficients: Mapping[tuple[int, int], Fraction | int | float | str]
    ) -> None:
        cleaned: dict[tuple[int, int], Fraction] = {}
        for key, value in coefficients.items():
            m, n = key
            if not isinstance(m, int) or not isinstance(n, int) or m < 0 or n < 0:
                raise ValueError("monomial exponents must be nonnegative integers")
            coefficient = Fraction(value)
            if coefficient != 0:
                cleaned[(m, n)] = coefficient
        self._coefficients = cleaned

    @property
    def coefficients(self) -> dict[tuple[int, int], Fraction]:
        return dict(self._coefficients)

    @property
    def is_zero(self) -> bool:
        return not self._coefficients

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GaussianPolynomial):
            return NotImplemented
        return self._coefficients == other._coefficients

    def __repr__(self) -> str:  # pragma: no cover - debugging aid only
        return f"GaussianPolynomial({self._coefficients!r})"

    # -- exact dictionary arithmetic ----------------------------------------- #

    def add(self, other: "GaussianPolynomial") -> "GaussianPolynomial":
        total = dict(self._coefficients)
        for key, value in other._coefficients.items():
            total[key] = total.get(key, Fraction(0)) + value
        return GaussianPolynomial(total)

    def scale(self, factor: Fraction | int) -> "GaussianPolynomial":
        factor = Fraction(factor)
        return GaussianPolynomial(
            {key: factor * value for key, value in self._coefficients.items()}
        )

    def mul_s(self) -> "GaussianPolynomial":
        return GaussianPolynomial(
            {(m + 1, n): value for (m, n), value in self._coefficients.items()}
        )

    def mul_z(self) -> "GaussianPolynomial":
        return GaussianPolynomial(
            {(m, n + 1): value for (m, n), value in self._coefficients.items()}
        )

    def d_s(self, alpha: Fraction) -> "GaussianPolynomial":
        r"""The polynomial ``P_s - \alpha P``, so that
        ``\partial_s(P\,E) = (P_s - \alpha P)\,E``."""
        alpha = Fraction(alpha)
        total: dict[tuple[int, int], Fraction] = {}
        for (m, n), value in self._coefficients.items():
            if m > 0:
                key = (m - 1, n)
                total[key] = total.get(key, Fraction(0)) + m * value
            key = (m, n)
            total[key] = total.get(key, Fraction(0)) - alpha * value
        return GaussianPolynomial(total)

    def d_z(self, beta: Fraction) -> "GaussianPolynomial":
        r"""The polynomial ``P_z - 2\beta z P``, so that
        ``\partial_z(P\,E) = (P_z - 2\beta zP)\,E``."""
        beta = Fraction(beta)
        total: dict[tuple[int, int], Fraction] = {}
        for (m, n), value in self._coefficients.items():
            if n > 0:
                key = (m, n - 1)
                total[key] = total.get(key, Fraction(0)) + n * value
            key = (m, n + 1)
            total[key] = total.get(key, Fraction(0)) - 2 * beta * value
        return GaussianPolynomial(total)

    # -- evaluation ---------------------------------------------------------- #

    def evaluate(self, s: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        """Vectorised float evaluation of the bare polynomial (no envelope)."""
        s_array = np.asarray(s, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(s_array, z_array).shape
        total = np.zeros(shape, dtype=np.float64)
        for (m, n), value in self._coefficients.items():
            total = total + float(value) * s_array**m * z_array**n
        return total

    @property
    def axial_parity(self) -> str:
        """``'odd'``, ``'even'`` or ``'none'`` — the parity in ``z``.

        Read off the exponent dictionary, not sampled: a polynomial is odd in
        ``z`` exactly when every monomial carries an odd power of ``z``.  The
        zero polynomial is reported ``'even'``; it contributes nothing, so the
        label is inert.
        """
        parities = {n % 2 for (_, n) in self._coefficients}
        if parities <= {0}:
            return "even"
        if parities == {1}:
            return "odd"
        return "none"


@dataclass(frozen=True)
class GaussianGenerator:
    r"""One Hermite-type generator ``a\,P(s,z)\,e^{-\alpha s-\beta z^2}``.

    The amplitude is carried as a separate float field, exactly as
    :class:`~ns_certificate_lab.mixed_initial_data.GeneratorComponent` does,
    so :meth:`GaussianMixedFamily.with_amplitudes` is a
    :func:`dataclasses.replace` and the polynomial's rational coefficients are
    never touched by the optimiser: the *shape* stays exact while the search
    moves only the amplitudes.
    """

    polynomial: GaussianPolynomial
    alpha: Fraction
    beta: Fraction
    amplitude: float = 1.0

    def __post_init__(self) -> None:
        polynomial = self.polynomial
        if not isinstance(polynomial, GaussianPolynomial):
            polynomial = GaussianPolynomial(polynomial)
            object.__setattr__(self, "polynomial", polynomial)
        object.__setattr__(self, "alpha", Fraction(self.alpha))
        object.__setattr__(self, "beta", Fraction(self.beta))
        if self.alpha <= 0 or self.beta <= 0:
            raise ValueError("the Gaussian widths alpha and beta must be positive")
        if not math.isfinite(self.amplitude):
            raise ValueError("amplitude must be finite")

    def _envelope(self, s: FloatArray, z: FloatArray) -> FloatArray:
        return np.exp(-float(self.alpha) * s - float(self.beta) * z * z)

    def evaluate(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        s = r_array * r_array
        return (
            self.amplitude
            * self.polynomial.evaluate(s, z_array)
            * self._envelope(s, z_array)
        )

    def partials(self, r: npt.ArrayLike, z: npt.ArrayLike) -> dict[str, FloatArray]:
        r"""``value``, ``ds``, ``dss``, ``dz``, ``dzz``, ``dsz`` — all exact.

        Each partial is obtained by the *polynomial* recursion first —
        ``\partial_s`` maps ``P`` to ``P_s - \alpha P`` and ``\partial_z`` maps
        it to ``P_z - 2\beta zP`` — and evaluated only afterwards, so the six
        arrays are values of closed-form Hermite functions, not difference
        quotients of anything.  The keys and their meaning match
        :meth:`GeneratorComponent.partials
        <ns_certificate_lab.mixed_initial_data.GeneratorComponent.partials>`
        exactly, which is what lets :class:`GaussianMixedFamily` reuse the
        audited gradient formulas verbatim.
        """
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        s = r_array * r_array
        envelope = self._envelope(s, z_array)
        value = self.polynomial
        ds = value.d_s(self.alpha)
        dz = value.d_z(self.beta)
        a = self.amplitude
        return {
            "value": a * value.evaluate(s, z_array) * envelope,
            "ds": a * ds.evaluate(s, z_array) * envelope,
            "dss": a * ds.d_s(self.alpha).evaluate(s, z_array) * envelope,
            "dz": a * dz.evaluate(s, z_array) * envelope,
            "dzz": a * dz.d_z(self.beta).evaluate(s, z_array) * envelope,
            "dsz": a * ds.d_z(self.beta).evaluate(s, z_array) * envelope,
        }

    @property
    def axial_parity(self) -> str:
        """The parity of the generator in ``z``, read off the exponents."""
        return self.polynomial.axial_parity

    def as_dict(self) -> dict[str, object]:
        """A JSON-serialisable record; fractions are rendered as strings so the
        exact rational data survive a round trip through a results file."""
        return {
            "amplitude": float(self.amplitude),
            "alpha": str(self.alpha),
            "beta": str(self.beta),
            "polynomial": {
                f"{m},{n}": str(value)
                for (m, n), value in sorted(self.polynomial.coefficients.items())
            },
            "axial_parity": self.axial_parity,
        }


@dataclass(frozen=True)
class GaussianMixedFamily:
    r"""A mixed family over Gaussian--Hermite generators.

    Duck-type compatible with
    :class:`~ns_certificate_lab.mixed_initial_data.MixedFamily` at every point
    :func:`~ns_certificate_lab.l3_optimizer.evaluate_shape` touches:
    ``with_amplitudes``, ``field``, ``reference_length``, ``name``,
    ``amplitudes``, ``as_dict`` and the ``clay_admissible`` flag.
    """

    name: str
    swirl: tuple[GaussianGenerator, ...]
    stream: tuple[GaussianGenerator, ...]
    reference_length: float

    #: ``True``: Gaussian--Hermite data are Schwartz class, hence legitimate
    #: smooth decaying initial data, not finite-``C^k`` surrogates.  The guard
    #: in :func:`~ns_certificate_lab.l3_optimizer.require_clay_admissible`
    #: accepts them.  This says nothing beyond admissibility of the *datum*;
    #: no claim about the Clay problem is made or implied anywhere here.
    clay_admissible: bool = True

    def __post_init__(self) -> None:
        if not self.swirl and not self.stream:
            raise ValueError("a family needs at least one generator")
        if self.reference_length <= 0.0:
            raise ValueError("reference_length must be positive")

    # -- generators ---------------------------------------------------------- #

    def partials(
        self, components: tuple[GaussianGenerator, ...], r: npt.ArrayLike,
        z: npt.ArrayLike,
    ) -> dict[str, FloatArray]:
        """Summed exact partials of a generator group, keys as in
        :meth:`GaussianGenerator.partials`."""
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        rb = np.broadcast_to(r_array, shape)
        zb = np.broadcast_to(z_array, shape)
        keys = ("value", "ds", "dss", "dz", "dzz", "dsz")
        total = {key: np.zeros(shape, dtype=np.float64) for key in keys}
        for component in components:
            piece = component.partials(rb, zb)
            for key in keys:
                total[key] = total[key] + piece[key]
        return total

    def u1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        return self.partials(self.swirl, r, z)["value"]

    def psi1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        return self.partials(self.stream, r, z)["value"]

    # -- the field and its exact gradient ------------------------------------ #

    def exact_gradient(self, grid: AxisymmetricGrid) -> dict[str, FloatArray]:
        r"""The nine physical components of ``\nabla u``, analytically.

        These are byte-for-byte the formulas of
        :meth:`MixedFamily.exact_gradient
        <ns_certificate_lab.mixed_initial_data.MixedFamily.exact_gradient>`,
        because the generators expose the same six partials with the same
        meaning; only the mechanism *behind* the partials differs (polynomial
        recursion here, chi-bump chain rule there).  In particular the
        divergence ``rr + tt + zz`` cancels algebraically, so the field is
        divergence free to rounding, not to truncation.
        """
        r_mesh, z_mesh = grid.mesh()
        swirl = self.partials(self.swirl, r_mesh, z_mesh)
        stream = self.partials(self.stream, r_mesh, z_mesh)
        s = r_mesh**2
        return {
            "rr": -(stream["dz"] + 2.0 * s * stream["dsz"]),
            "rt": swirl["value"] + 2.0 * s * swirl["ds"],
            "rz": 2.0 * r_mesh * (4.0 * stream["ds"] + 2.0 * s * stream["dss"]),
            "tr": -swirl["value"],
            "tt": -stream["dz"],
            "tz": np.zeros_like(r_mesh),
            "zr": -r_mesh * stream["dzz"],
            "zt": r_mesh * swirl["dz"],
            "zz": 2.0 * stream["dz"] + 2.0 * s * stream["dsz"],
        }

    def field(self, grid: AxisymmetricGrid):
        r"""The :class:`~ns_certificate_lab.l3_generation.MixedField`, with the
        analytic gradient attached so nothing downstream ever differentiates a
        Gaussian by finite differences on the certified path."""
        from .l3_generation import MixedField

        r_mesh, z_mesh = grid.mesh()
        swirl = self.partials(self.swirl, r_mesh, z_mesh)
        stream = self.partials(self.stream, r_mesh, z_mesh)
        u_theta = r_mesh * swirl["value"]
        u_r = -r_mesh * stream["dz"]
        u_z = 2.0 * stream["value"] + 2.0 * (r_mesh**2) * stream["ds"]
        return MixedField(
            grid=grid, u_r=u_r, u_theta=u_theta, u_z=u_z,
            analytic_gradient=self.exact_gradient(grid),
        )

    def velocity(
        self, r: npt.ArrayLike, z: npt.ArrayLike
    ) -> tuple[FloatArray, FloatArray, FloatArray]:
        """``(u^r, u^theta, u^z)`` at arbitrary meridional points.

        Off-grid evaluation exists so the tail bound can be *tested* against
        the actual field on spheres that need not be grid-aligned; the bound
        itself never samples anything.
        """
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        rb = np.broadcast_to(r_array, shape)
        swirl = self.partials(self.swirl, r_array, z_array)
        stream = self.partials(self.stream, r_array, z_array)
        u_theta = rb * swirl["value"]
        u_r = -rb * stream["dz"]
        u_z = 2.0 * stream["value"] + 2.0 * (rb**2) * stream["ds"]
        return u_r, u_theta, u_z

    # -- the rigorous tail bound --------------------------------------------- #

    def _tail_terms(self) -> list[tuple[float, float, float, float]]:
        r"""Monomial data ``(|c|, a, b, gamma)`` for every velocity term.

        Each velocity component is a finite sum of terms ``c\,s^a|z|^b
        e^{-\alpha s-\beta z^2}`` with half-integer ``a``: the swirl gives
        ``u^\theta = r\,P E`` (one factor ``r = s^{1/2}``), the stream gives
        ``u^r = -r\,(P_z - 2\beta zP)E`` and ``u^z = (2P + 2s(P_s - \alpha
        P))E``, all produced by the same exact recursions as the gradient, so
        the tail bound bounds precisely the field the solver sees.
        """
        terms: list[tuple[float, float, float, float]] = []

        def collect(
            poly: GaussianPolynomial, amplitude: float, shift: float, gamma: float
        ) -> None:
            for (m, n), value in poly.coefficients.items():
                coefficient = abs(amplitude * float(value))
                if coefficient > 0.0:
                    terms.append((coefficient, m + shift, float(n), gamma))

        for generator in self.swirl:
            gamma = float(min(generator.alpha, generator.beta))
            collect(generator.polynomial, generator.amplitude, 0.5, gamma)
        for generator in self.stream:
            gamma = float(min(generator.alpha, generator.beta))
            collect(
                generator.polynomial.d_z(generator.beta),
                generator.amplitude, 0.5, gamma,
            )
            axial = generator.polynomial.scale(2).add(
                generator.polynomial.d_s(generator.alpha).mul_s().scale(2)
            )
            collect(axial, generator.amplitude, 0.0, gamma)
        return terms

    def tail_bound(self, radius: float) -> float:
        r"""A certified upper bound on ``\sup_{|x|\ge\rho}|u(x)|``.

        **The exact bound implemented.**  Every velocity component is a finite
        sum of terms ``c\,s^a|z|^b e^{-\alpha s-\beta z^2}`` with ``s = r^2``
        and half-integer ``a`` (see :meth:`_tail_terms`).  Two elementary facts
        bound one term on the sphere ``s + z^2 = \rho'^2``.  First, with ``t =
        z^2`` the constrained maximum of ``s^a t^{b/2}`` subject to ``s + t =
        \rho'^2`` is, by weighted AM--GM,

        .. math::

           \max s^a t^{b/2} = K\,\rho'^{\,2a+b}, \qquad
           K = \frac{a^a\,(b/2)^{b/2}}{(a+b/2)^{a+b/2}}, \quad 0^0 := 1 .

        Second, ``\alpha s + \beta t \ge \gamma(s+t) = \gamma\rho'^2`` with
        ``\gamma = \min(\alpha,\beta)``, an identity (no loss) whenever
        ``\alpha = \beta`` as in the preregistered basis.  Hence on that sphere
        ``|{\rm term}| \le |c|\,K\rho'^{\,p}e^{-\gamma\rho'^2}`` with ``p = 2a
        + b``, and since ``\rho^pe^{-\gamma\rho^2}`` is decreasing once
        ``\rho^2 \ge p/(2\gamma)`` (its calculus maximum), the supremum over
        the *entire exterior* ``|x| \ge \rho`` is the same expression evaluated
        at ``\max(\rho, \sqrt{p/(2\gamma)})``.  Summing over terms and using
        ``|u| \le |u^r| + |u^\theta| + |u^z|`` gives the bound.

        This is the sharp form of the cruder split ``s^mz^n e^{-\alpha s/2 -
        \beta z^2/2} \le (2m/\alpha e)^m (n/\beta e)^{n/2}`` times ``e^{-\gamma
        \rho^2/2}``: keeping the constrained maximum *on the sphere* retains
        the full decay rate ``e^{-\gamma\rho^2}`` instead of half of it, which
        at the standard box is the difference between a useful truncation bound
        and a vacuous one.  It is a true bound at every radius, not an
        asymptotic, and it is monotone nonincreasing in ``radius``.
        """
        if not math.isfinite(radius) or radius < 0.0:
            raise ValueError("radius must be finite and nonnegative")
        total = 0.0
        radius_sq = radius * radius
        for coefficient, a, b, gamma in self._tail_terms():
            half = 0.5 * b
            weight = a + half
            k_factor = _power(a, a) * _power(half, half) / _power(weight, weight)
            p = 2.0 * weight
            rho_sq = max(radius_sq, p / (2.0 * gamma)) if p > 0.0 else radius_sq
            total += (
                coefficient * k_factor * _power(rho_sq, 0.5 * p)
                * math.exp(-gamma * rho_sq)
            )
        return total

    def support_radius(self, threshold: float) -> float:
        r"""The radius beyond which the **analytic** tail bound certifies
        ``|u| < threshold``.

        Compact support does not exist here, so this is the honest replacement
        for :attr:`MixedFamily.support_radius
        <ns_certificate_lab.mixed_initial_data.MixedFamily.support_radius>`:
        outside the returned sphere the field is *provably*, not empirically,
        below the threshold.  Computed by bisection on the monotone
        :meth:`tail_bound`, always returning the side on which the bound has
        been evaluated and found ``\le`` threshold, so conservativeness is
        preserved by construction.  An identically zero family returns ``0``.
        """
        if not math.isfinite(threshold) or threshold <= 0.0:
            raise ValueError("threshold must be positive and finite")
        if self.tail_bound(0.0) <= threshold:
            return 0.0
        low, high = 0.0, 1.0
        while self.tail_bound(high) > threshold:
            low = high
            high *= 2.0
        for _ in range(100):
            midpoint = 0.5 * (low + high)
            if self.tail_bound(midpoint) <= threshold:
                high = midpoint
            else:
                low = midpoint
        return high

    # -- bookkeeping --------------------------------------------------------- #

    @property
    def amplitudes(self) -> FloatArray:
        return np.array(
            [g.amplitude for g in (*self.swirl, *self.stream)], dtype=np.float64
        )

    def with_amplitudes(self, values: npt.ArrayLike) -> "GaussianMixedFamily":
        """A copy with every generator amplitude replaced, in order.

        Same contract as :meth:`MixedFamily.with_amplitudes
        <ns_certificate_lab.mixed_initial_data.MixedFamily.with_amplitudes>`:
        the optimiser's design variables are exactly these amplitudes, swirl
        first, and every point of the amplitude space is automatically
        ``C^\\infty``, divergence free and Schwartz.
        """
        array = np.asarray(values, dtype=np.float64).ravel()
        total = len(self.swirl) + len(self.stream)
        if array.size != total:
            raise ValueError(f"expected {total} amplitudes, received {array.size}")
        swirl = tuple(
            replace(g, amplitude=float(a))
            for g, a in zip(self.swirl, array[: len(self.swirl)])
        )
        stream = tuple(
            replace(g, amplitude=float(a))
            for g, a in zip(self.stream, array[len(self.swirl):])
        )
        return replace(self, swirl=swirl, stream=stream)

    def scaled(self, *, swirl: float = 1.0, stream: float = 1.0) -> "GaussianMixedFamily":
        """A copy with the two generator groups rescaled independently."""
        return replace(
            self,
            swirl=tuple(
                replace(g, amplitude=g.amplitude * swirl) for g in self.swirl
            ),
            stream=tuple(
                replace(g, amplitude=g.amplitude * stream) for g in self.stream
            ),
        )

    @property
    def swirl_amplitude(self) -> float:
        return max((abs(g.amplitude) for g in self.swirl), default=0.0)

    @property
    def meridional_amplitude(self) -> float:
        return max((abs(g.amplitude) for g in self.stream), default=0.0)

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "reference_length": self.reference_length,
            "clay_admissible": self.clay_admissible,
            "swirl_amplitude": self.swirl_amplitude,
            "meridional_amplitude": self.meridional_amplitude,
            "swirl": [g.as_dict() for g in self.swirl],
            "stream": [g.as_dict() for g in self.stream],
        }


# --------------------------------------------------------------------------- #
# the preregistered search basis                                               #
# --------------------------------------------------------------------------- #


def gaussian_search_basis(
    *,
    alpha_narrow: Fraction = Fraction(3),
    beta_narrow: Fraction = Fraction(3),
    alpha_wide: Fraction = Fraction(5, 2),
    beta_wide: Fraction = Fraction(5, 2),
    reference_length: float = 1.0,
) -> GaussianMixedFamily:
    r"""The preregistered eight-amplitude Gaussian--Hermite basis.

    Four swirl and four stream generators, amplitudes all zero (the optimiser
    owns them), over two Gaussian envelopes ``E_1 = e^{-3s-3z^2}`` and ``E_2 =
    e^{-\frac52 s-\frac52 z^2}``:

    ========  =====================  =====================
    slot      polynomial             envelope
    ========  =====================  =====================
    swirl 1   ``z``                  ``E_1``
    swirl 2   ``sz``                 ``E_1``
    swirl 3   ``z^3``                ``E_1``
    swirl 4   ``sz``                 ``E_2``
    stream 1  ``z``                  ``E_1``
    stream 2  ``sz``                 ``E_1``
    stream 3  ``z^3``                ``E_1``
    stream 4  ``z``                  ``E_2``
    ========  =====================  =====================

    **Every stream polynomial contains only odd powers of** ``z``: the parity
    selection rule of :mod:`ns_certificate_lab.mixed_initial_data` kills the
    pressure term of the generation rate for an even stream generator on a
    ``z``-symmetric domain, so an even stream slot would be a wasted design
    variable, structurally incapable of the one thing the basis is for.  The
    swirl parity is irrelevant to the rule and the odd choice simply matches
    the strain-flow geometry of the chi families.

    **Why the widths are what they are.**  The standard box (``r_{\max} = 3``,
    ``z_{\max} = 3.6``) contains the ball of radius ``3``, so every exterior
    point lies at distance at least ``3`` from the origin and the truncation
    committed by working on the box is controlled by ``tail\_bound(3)``.  The
    naive width bookkeeping fails here: the bare exponential ``e^{-\alpha s}``
    with ``\alpha = 1`` is only ``e^{-9} \approx 1.2\times10^{-4}`` at ``r =
    3``, nowhere near the ``10^{-7}`` target, so the widest admissible Gaussian
    must be much narrower than ``e^{-s-z^2}``; and even ``\alpha = \beta = 2``
    leaves the *bound* — polynomial prefactors included — at
    ``1.5\times10^{-5}`` of the peak speed, still two orders of magnitude too
    slack.  With ``\alpha = \beta = 3`` for ``E_1`` and ``\alpha = \beta =
    \frac52`` for ``E_2``, measured with unit amplitudes on the standard
    ``41\times81`` grid, the certified tail bound at radius ``3`` is
    ``2.7\times10^{-8}`` absolute against a peak speed of ``1.13``, a
    **relative truncation level of ``2.4\times10^{-8}``**, below the
    ``10^{-7}`` requirement for both envelopes (the ``E_2`` terms alone
    contribute ``2.4\times10^{-8}`` absolute, which is why ``\frac52`` and not
    ``\frac94``: at ``\frac94`` the ``E_2`` share alone is ``1.8\times10^{-7}``
    relative and breaches the target).  Both envelopes are isotropic in
    ``(s, z^2)`` (``\alpha = \beta``) so the ``\gamma = \min(\alpha,\beta)``
    step of the tail bound is lossless, and all widths are rational so the
    generator data stay exact.

    ``reference_length`` is ``1.0``: the ``E_1`` scale ``1/\sqrt3`` and the
    ``E_2`` scale ``\sqrt{2/5}`` bracket it, and the critical-Reynolds
    objective carries a single explicit power of this length, so a documented
    order-one convention is what matters, not a fitted value.
    """
    e1 = (alpha_narrow, beta_narrow)
    e2 = (alpha_wide, beta_wide)
    z1 = GaussianPolynomial({(0, 1): 1})
    sz = GaussianPolynomial({(1, 1): 1})
    z3 = GaussianPolynomial({(0, 3): 1})
    swirl = (
        GaussianGenerator(polynomial=z1, alpha=e1[0], beta=e1[1], amplitude=0.0),
        GaussianGenerator(polynomial=sz, alpha=e1[0], beta=e1[1], amplitude=0.0),
        GaussianGenerator(polynomial=z3, alpha=e1[0], beta=e1[1], amplitude=0.0),
        GaussianGenerator(polynomial=sz, alpha=e2[0], beta=e2[1], amplitude=0.0),
    )
    stream = (
        GaussianGenerator(polynomial=z1, alpha=e1[0], beta=e1[1], amplitude=0.0),
        GaussianGenerator(polynomial=sz, alpha=e1[0], beta=e1[1], amplitude=0.0),
        GaussianGenerator(polynomial=z3, alpha=e1[0], beta=e1[1], amplitude=0.0),
        GaussianGenerator(polynomial=z1, alpha=e2[0], beta=e2[1], amplitude=0.0),
    )
    return GaussianMixedFamily(
        name="gaussian-hermite-basis", swirl=swirl, stream=stream,
        reference_length=reference_length,
    )
