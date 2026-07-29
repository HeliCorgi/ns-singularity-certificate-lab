r"""Interval evaluation of the smooth generators over a **cell**, not at a node.

The rule the certificates must obey is that a derivative bound may never be
inferred from stored nodal values.  For the initial datum that rule is easy to
keep, because the datum is not stored: it is a finite combination of explicit
generators, so every quantity is a closed-form composition of ``exp``, powers
and rational functions of ``(s, z)`` with ``s = r^2``.  Evaluating that
composition in interval arithmetic on the box
``[r_i, r_{i+1}]\times[z_j, z_{j+1}]`` returns an enclosure valid at **every
point of the cell interior**, derivatives included, with no Lipschitz
hypothesis and no inflation factor.

Everything factors
------------------
A generator is ``a\,X(s)\,Z(z)``, so every partial factors as
``a\,X^{(m)}(s)\,Z^{(n)}(z)`` and only six one-dimensional interval functions
are needed: ``X, X', X''`` and ``Z, Z', Z''``.  With

.. math::
   \chi(\sigma) = e^{-1/(1-\sigma)},\quad
   \chi' = -\frac{\chi}{(1-\sigma)^2},\quad
   \chi'' = \chi\Bigl(\frac1{(1-\sigma)^4} - \frac2{(1-\sigma)^3}\Bigr),

``X(s) = \chi(s/R^2)`` and ``Z(z) = \chi(\zeta^2)\,\Pi(\zeta)`` with
``\zeta = (z-z_0)/W``, the rest is interval arithmetic.

The dependency problem
----------------------
Interval arithmetic over-estimates when a variable occurs more than once, and
these expressions repeat ``\sigma`` and ``\zeta`` several times.  Two things
control it.  ``\chi`` is evaluated by **monotonicity** — it is strictly
decreasing — so its own enclosure is sharp; and the caller subdivides cells,
which shrinks the over-estimate linearly in the subdivision width.  No attempt
is made to hide the widening: :func:`cell_enclosure` returns intervals whose
width the caller can inspect and narrow.

Only unit generator powers
--------------------------
``radial_power`` and ``axial_power`` are required to be ``1`` here.  General
powers would need ``\chi^{p}`` and two derivatives of it, whose interval form
has three more repeated occurrences of ``\chi``; the certificate would still be
valid but much wider, and the search basis does not use them.  A non-unit power
raises rather than silently widening.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

__all__ = [
    "CellEnclosure",
    "GeneratorIntervals",
    "cell_enclosure",
    "clear_generator_cache",
    "chi_interval",
    "chi_first_interval",
    "chi_second_interval",
    "divide",
    "square",
]

_ZERO = Interval(Fraction(0), Fraction(0))
_ONE = Interval(Fraction(1), Fraction(1))


def _round(value: Interval, bits: int) -> Interval:
    return value.round_outward(bits)


def square(value: Interval) -> Interval:
    """``x^2``, using the sign structure so the result is sharp."""
    if value.lower >= 0:
        return Interval(value.lower**2, value.upper**2)
    if value.upper <= 0:
        return Interval(value.upper**2, value.lower**2)
    return Interval(Fraction(0), max(value.lower**2, value.upper**2))


def divide(numerator: Interval, denominator: Interval) -> Interval:
    """``x/y`` for a denominator that does not contain zero."""
    if denominator.contains_zero:
        raise ZeroDivisionError("interval division by an interval containing zero")
    reciprocal = Interval(
        min(1 / denominator.lower, 1 / denominator.upper),
        max(1 / denominator.lower, 1 / denominator.upper),
    )
    return numerator * reciprocal


def chi_interval(argument: Interval, *, terms: int, bits: int) -> Interval:
    r"""``\chi(\sigma)`` on an interval, sharp by monotonicity.

    ``\chi`` is strictly decreasing on ``[0,1)`` and zero on ``[1,\infty)``, so
    the range on ``[a,b]`` is ``[\chi(b),\chi(a)]``.  Using that rather than
    generic interval arithmetic removes the dependency problem from the one
    place it would hurt most.
    """
    from .l3_certificate import exp_interval

    low = max(argument.lower, Fraction(0))
    high = argument.upper
    if low >= 1:
        return _ZERO
    upper_point = -1 / (1 - low)
    upper = exp_interval(Interval(upper_point, upper_point), terms=terms).upper
    if high >= 1:
        return _round(Interval(Fraction(0), upper), bits)
    lower_point = -1 / (1 - high)
    lower = exp_interval(Interval(lower_point, lower_point), terms=terms).lower
    return _round(Interval(lower, upper), bits)


def _gap(argument: Interval) -> Interval:
    """``1 - \\sigma`` clipped to the support, guaranteed positive."""
    low = max(argument.lower, Fraction(0))
    high = min(argument.upper, Fraction(1))
    if high >= 1:
        high = Fraction(1)
    return Interval(1 - high, 1 - low)


def chi_first_interval(argument: Interval, *, terms: int, bits: int) -> Interval:
    r"""``\chi'(\sigma) = -\chi/(1-\sigma)^2``."""
    if argument.lower >= 1:
        return _ZERO
    chi = chi_interval(argument, terms=terms, bits=bits)
    gap = _gap(argument)
    if gap.contains_zero:
        # chi decays faster than any power, so chi/(1-s)^2 -> 0 as s -> 1.
        # sup_{0<g<=G} e^{-1/g}/g^2 is attained at g = 1/2 with value 4/e^2 < 1.
        magnitude = max(chi.magnitude * 4, Fraction(1))
        return _round(Interval(-magnitude, Fraction(0)), bits)
    return _round(-divide(chi, square(gap)), bits)


def chi_second_interval(argument: Interval, *, terms: int, bits: int) -> Interval:
    r"""``\chi''(\sigma) = \chi\bigl((1-\sigma)^{-4} - 2(1-\sigma)^{-3}\bigr)``."""
    if argument.lower >= 1:
        return _ZERO
    chi = chi_interval(argument, terms=terms, bits=bits)
    gap = _gap(argument)
    if gap.contains_zero:
        # sup_{0<g<=1} e^{-1/g}(g^{-4} + 2g^{-3}) is below 12; use it as a crude
        # but valid bound rather than dividing by an interval containing zero.
        magnitude = max(chi.magnitude * 12, Fraction(12))
        return _round(Interval(-magnitude, magnitude), bits)
    inverse = divide(_ONE, gap)
    quad = square(square(inverse))
    cube = inverse * square(inverse)
    return _round(chi * (quad - cube.scale(Fraction(2))), bits)


@dataclass(frozen=True)
class GeneratorIntervals:
    """Enclosures of one generator group and its partials over a cell."""

    value: Interval
    ds: Interval
    dss: Interval
    dz: Interval
    dzz: Interval
    dsz: Interval

    def __add__(self, other: "GeneratorIntervals") -> "GeneratorIntervals":
        return GeneratorIntervals(
            *(getattr(self, k) + getattr(other, k)
              for k in ("value", "ds", "dss", "dz", "dzz", "dsz"))
        )

    @staticmethod
    def zero() -> "GeneratorIntervals":
        return GeneratorIntervals(_ZERO, _ZERO, _ZERO, _ZERO, _ZERO, _ZERO)


#: Cache for the one-dimensional factors.
#:
#: A grid sweep evaluates the same radial box against every axial box and the
#: same axial box against every radial box, so the transcendental work is
#: ``O(n_r + n_z)`` per component rather than ``O(n_r n_z)`` once the factors are
#: cached.  Without this the exponential series dominates everything: measured at
#: 209 ms per cell uncached, which makes even a coarse certificate impractical.
_RADIAL_CACHE: dict[tuple, tuple[Interval, Interval, Interval]] = {}
_AXIAL_CACHE: dict[tuple, tuple[Interval, Interval, Interval]] = {}


def clear_generator_cache() -> None:
    """Empty the factor caches.  Only tests should need this."""
    _RADIAL_CACHE.clear()
    _AXIAL_CACHE.clear()


def _radial_factors(
    support: Fraction, r_box: Interval, *, terms: int, bits: int
) -> tuple[Interval, Interval, Interval]:
    r"""``(X, X_s, X_{ss})`` for ``X(s) = \chi(s/R^2)``, ``s = r^2``."""
    key = (support, r_box.lower, r_box.upper, terms, bits)
    cached = _RADIAL_CACHE.get(key)
    if cached is not None:
        return cached
    radius_sq = support * support
    sigma = square(r_box).scale(1 / radius_sq)
    value = chi_interval(sigma, terms=terms, bits=bits)
    first = chi_first_interval(sigma, terms=terms, bits=bits).scale(1 / radius_sq)
    second = chi_second_interval(sigma, terms=terms, bits=bits).scale(
        1 / (radius_sq * radius_sq)
    )
    _RADIAL_CACHE[key] = (value, first, second)
    return value, first, second


def _axial_factors(
    component, z_box: Interval, *, terms: int, bits: int
) -> tuple[Interval, Interval, Interval]:
    r"""``(Z, Z_z, Z_{zz})`` for ``Z(z) = \chi(\zeta^2)\,\Pi(\zeta)``."""
    width = Fraction(float(component.axial_support))
    centre = Fraction(float(component.axial_center))
    concentration = Fraction(float(component.axial_concentration))
    key = (
        width, centre, concentration, bool(component.odd_axial),
        z_box.lower, z_box.upper, terms, bits,
    )
    cached = _AXIAL_CACHE.get(key)
    if cached is not None:
        return cached

    zeta = Interval(z_box.lower - centre, z_box.upper - centre).scale(1 / width)
    argument = square(zeta)
    bump = chi_interval(argument, terms=terms, bits=bits)
    bump_1 = chi_first_interval(argument, terms=terms, bits=bits)
    bump_2 = chi_second_interval(argument, terms=terms, bits=bits)
    # d/dzeta chi(zeta^2) = chi'(zeta^2) 2 zeta ; second derivative
    # chi''(zeta^2) 4 zeta^2 + chi'(zeta^2) 2 .
    d_bump = bump_1 * zeta.scale(Fraction(2))
    d2_bump = bump_2 * argument.scale(Fraction(4)) + bump_1.scale(Fraction(2))

    if component.odd_axial:
        c = concentration
        denominator = _ONE + argument.scale(c)
        profile = divide(zeta, denominator)
        profile_1 = divide(_ONE - argument.scale(c), square(denominator))
        profile_2 = divide(
            (-zeta.scale(Fraction(2) * c))
            * (Interval(Fraction(3), Fraction(3)) - argument.scale(c)),
            denominator * square(denominator),
        )
    else:
        profile, profile_1, profile_2 = _ONE, _ZERO, _ZERO

    value = _round(bump * profile, bits)
    first = _round((d_bump * profile + bump * profile_1).scale(1 / width), bits)
    second = _round(
        (d2_bump * profile + (d_bump * profile_1).scale(Fraction(2))
         + bump * profile_2).scale(1 / (width * width)),
        bits,
    )
    _AXIAL_CACHE[key] = (value, first, second)
    return value, first, second


def _component_intervals(
    component, r_box: Interval, z_box: Interval, *, terms: int, bits: int
) -> GeneratorIntervals:
    if component.radial_power != 1.0 or component.axial_power != 1.0:
        raise ValueError(
            "the interval generator evaluation supports unit bump powers only; "
            "a non-unit power would widen every enclosure and the search basis "
            "does not use one"
        )
    chi, chi_1, chi_2 = _radial_factors(
        Fraction(float(component.radial_support)), r_box, terms=terms, bits=bits
    )
    z_value, z_first, z_second = _axial_factors(
        component, z_box, terms=terms, bits=bits
    )
    amplitude = Fraction(float(component.amplitude))
    return GeneratorIntervals(
        value=_round((chi * z_value).scale(amplitude), bits),
        ds=_round((chi_1 * z_value).scale(amplitude), bits),
        dss=_round((chi_2 * z_value).scale(amplitude), bits),
        dz=_round((chi * z_first).scale(amplitude), bits),
        dzz=_round((chi * z_second).scale(amplitude), bits),
        dsz=_round((chi_1 * z_first).scale(amplitude), bits),
    )


@dataclass(frozen=True)
class CellEnclosure:
    """Enclosures of the velocity, its gradient and the ``J`` integrands."""

    speed_squared: Interval
    speed: Interval
    gradient_squared: Interval
    viscous_integrand: Interval
    flux_magnitude: Interval
    flux: Interval
    divergence: Interval

    def as_dict(self) -> dict[str, list[str]]:
        return {
            name: getattr(self, name).as_pair() for name in self.__dataclass_fields__
        }


def cell_enclosure(
    family,
    r_box: Interval,
    z_box: Interval,
    *,
    terms: int = 32,
    bits: int = DEFAULT_PRECISION_BITS,
) -> CellEnclosure:
    r"""Enclose everything ``J`` needs over one cell.

    The gradient components are the analytic ones of
    :meth:`~ns_certificate_lab.mixed_initial_data.MixedFamily.exact_gradient`;
    the divergence enclosure is returned so a caller can check that the
    *interval* divergence contains zero, which it must, the identity being
    exact.
    """
    swirl = GeneratorIntervals.zero()
    for component in family.swirl:
        swirl = swirl + _component_intervals(
            component, r_box, z_box, terms=terms, bits=bits
        )
    stream = GeneratorIntervals.zero()
    for component in family.stream:
        stream = stream + _component_intervals(
            component, r_box, z_box, terms=terms, bits=bits
        )

    s_box = _round(square(r_box), bits)
    two_s = s_box.scale(Fraction(2))

    u_theta = _round(r_box * swirl.value, bits)
    u_r = _round(-(r_box * stream.dz), bits)
    u_z = _round(stream.value.scale(Fraction(2)) + _round(two_s * stream.ds, bits), bits)

    gradient = {
        "rr": -(stream.dz + _round(two_s * stream.dsz, bits)),
        "rt": swirl.value + _round(two_s * swirl.ds, bits),
        "rz": _round(
            r_box.scale(Fraction(2))
            * _round(stream.ds.scale(Fraction(4)) + _round(two_s * stream.dss, bits),
                     bits),
            bits,
        ),
        "tr": -swirl.value,
        "tt": -stream.dz,
        "tz": _ZERO,
        "zr": -_round(r_box * stream.dzz, bits),
        "zt": _round(r_box * swirl.dz, bits),
        "zz": stream.dz.scale(Fraction(2)) + _round(two_s * stream.dsz, bits),
    }
    gradient = {key: _round(value, bits) for key, value in gradient.items()}

    speed_squared = _round(
        _round(square(u_r), bits) + _round(square(u_theta), bits)
        + _round(square(u_z), bits),
        bits,
    )
    speed = Interval(
        max(speed_squared.lower, Fraction(0)), max(speed_squared.upper, Fraction(0))
    )
    from .l3_certificate import sqrt_interval

    speed = _round(sqrt_interval(speed, bits=bits), bits)
    gradient_squared = _round(
        sum((_round(square(v), bits) for v in gradient.values()), _ZERO), bits
    )

    # |V| integrand.  Kato gives |grad |u|| <= |grad u| pointwise, so
    # q(|grad u|^2 + |grad q|^2) <= 2 q |grad u|^2, which needs no division by q
    # and is therefore valid on the zero set as well.
    viscous = _round((speed * gradient_squared).scale(Fraction(2)), bits)

    # |g| = |u . grad |u|| <= |u| |grad u|, with no division, so it is valid on
    # the zero set too.  This is the Kato bound again.
    magnitude = _round(
        speed * _round(sqrt_interval(gradient_squared, bits=bits), bits), bits
    )
    envelope = Interval(-magnitude.upper, magnitude.upper)

    # The signed flux g = u . grad |u| = (u^r G_r + u^z G_z)/|u| with
    # G = (1/2) grad(|u|^2), which is smooth.  When the speed enclosure stays
    # away from zero the quotient is taken directly and is much tighter than the
    # Kato envelope; when it straddles zero the quotient does not exist and the
    # envelope is the only rigorous statement.  Intersecting the two is free and
    # never worse.
    g_r = _round(
        _round(u_r * gradient["rr"], bits) + _round(u_theta * gradient["rt"], bits)
        + _round(u_z * gradient["rz"], bits),
        bits,
    )
    g_z = _round(
        _round(u_r * gradient["zr"], bits) + _round(u_theta * gradient["zt"], bits)
        + _round(u_z * gradient["zz"], bits),
        bits,
    )
    numerator = _round(
        _round(u_r * g_r, bits) + _round(u_z * g_z, bits), bits
    )
    if speed.lower > 0:
        quotient = _round(divide(numerator, speed), bits)
        flux = Interval(
            max(quotient.lower, envelope.lower), min(quotient.upper, envelope.upper)
        )
    else:
        flux = envelope

    divergence = _round(gradient["rr"] + gradient["tt"] + gradient["zz"], bits)
    return CellEnclosure(
        speed_squared=speed_squared,
        speed=speed,
        gradient_squared=gradient_squared,
        viscous_integrand=viscous,
        flux_magnitude=Interval(Fraction(0), max(magnitude.upper, Fraction(0))),
        flux=flux,
        divergence=divergence,
    )
