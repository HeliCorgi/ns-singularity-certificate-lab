r"""An exact rational three-mode Leray relay witness on ``T^3``.

The field

.. math::

   u = B e_3\sin(s p\cdot x)+C e_2\cos(s q\cdot x)
       +D n\cos(s c\cdot x),

with ``p=(1,1,0)``, ``q=(1,0,1)``, ``c=p+q=(2,1,1)`` and
``n=p cross q=(1,-1,-1)``, is exactly divergence free.  Direct rational
Fourier algebra gives

.. math::

   P_c\mathbb P((u\cdot\nabla)u)=-sBC n\cos(sc\cdot x)/3,
   \qquad \Pi_c=sBCD/2.

The difference mode ``p-q`` cancels exactly.  The witness proves only that a
true Navier--Stokes triad can have the desired signed transfer.  Scaling a
fixed number of modes does not provide the localized multiplicity needed for
a Zeno cascade; :func:`fixed_cardinality_scaling` records that obstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .fourier_torus import TrigVector, advection, leray
from .torus_chain import l2_inner

__all__ = [
    "ExactRelayMetrics",
    "build_exact_relay_triad",
    "exact_relay_metrics",
    "fixed_cardinality_scaling",
]


Wavevector = tuple[int, int, int]


def _fraction(value: Fraction | int) -> Fraction:
    result = Fraction(value)
    if result.denominator == 0:
        raise ValueError("invalid rational coefficient")
    return result


def _norm_squared(wave: Wavevector) -> int:
    return sum(component * component for component in wave)


def build_exact_relay_triad(
    *,
    scale: int = 1,
    parent_sine: Fraction | int = 1,
    parent_cosine: Fraction | int = 1,
    child_cosine: Fraction | int = Fraction(1, 8),
) -> TrigVector:
    """Return the exact rational parent-parent-child trigonometric field."""

    if isinstance(scale, bool) or not isinstance(scale, int) or scale < 1:
        raise ValueError("scale must be a positive integer")
    b = _fraction(parent_sine)
    c_coefficient = _fraction(parent_cosine)
    d = _fraction(child_cosine)
    p = (scale, scale, 0)
    q = (scale, 0, scale)
    child = (2 * scale, scale, scale)
    direction = (1, -1, -1)
    return TrigVector.from_modes(
        [
            (p, (0, 0, 0), (0, 0, b)),
            (q, (0, c_coefficient, 0), (0, 0, 0)),
            (
                child,
                tuple(d * component for component in direction),
                (0, 0, 0),
            ),
        ]
    )


@dataclass(frozen=True)
class ExactRelayMetrics:
    """Exact fractions for the signed relay and every immediate loss."""

    scale: int
    parent_energy: Fraction
    child_energy: Fraction
    child_flux: Fraction
    parent_flux: Fraction
    child_viscous_loss: Fraction
    parent_viscous_loss: Fraction
    child_net: Fraction
    off_chain_nonlinear_l2_squared: Fraction
    formula_off_chain_l2_squared: Fraction
    total_nonlinear_energy_defect: Fraction
    difference_mode_present: bool
    child_coefficient_matches_formula: bool

    def as_dict(self) -> dict[str, int | str | bool]:
        return {
            "scale": self.scale,
            "parent_energy": str(self.parent_energy),
            "child_energy": str(self.child_energy),
            "child_flux": str(self.child_flux),
            "parent_flux": str(self.parent_flux),
            "child_viscous_loss": str(self.child_viscous_loss),
            "parent_viscous_loss": str(self.parent_viscous_loss),
            "child_net": str(self.child_net),
            "off_chain_nonlinear_l2_squared": str(
                self.off_chain_nonlinear_l2_squared
            ),
            "formula_off_chain_l2_squared": str(
                self.formula_off_chain_l2_squared
            ),
            "total_nonlinear_energy_defect": str(
                self.total_nonlinear_energy_defect
            ),
            "difference_mode_present": self.difference_mode_present,
            "child_coefficient_matches_formula": self.child_coefficient_matches_formula,
        }


def exact_relay_metrics(
    *,
    viscosity: Fraction | int,
    scale: int = 1,
    parent_sine: Fraction | int = 1,
    parent_cosine: Fraction | int = 1,
    child_cosine: Fraction | int = Fraction(1, 8),
) -> ExactRelayMetrics:
    """Evaluate the exact field and independently compare the closed formulas."""

    nu = _fraction(viscosity)
    if nu <= 0:
        raise ValueError("viscosity must be positive")
    b = _fraction(parent_sine)
    c_coefficient = _fraction(parent_cosine)
    d = _fraction(child_cosine)
    field = build_exact_relay_triad(
        scale=scale,
        parent_sine=b,
        parent_cosine=c_coefficient,
        child_cosine=d,
    )
    divergence = field.divergence().cleaned()
    if divergence.terms:
        raise AssertionError("the relay witness lost exact divergence-freeness")
    nonlinear = leray(advection(field, field)).cleaned()

    parent_norm_squared = 2 * scale * scale
    child_norm_squared = 6 * scale * scale
    parent = field.restrict(lambda wave: _norm_squared(wave) == parent_norm_squared)
    child = field.restrict(lambda wave: _norm_squared(wave) == child_norm_squared)
    parent_nonlinear = nonlinear.restrict(
        lambda wave: _norm_squared(wave) == parent_norm_squared
    )
    child_nonlinear = nonlinear.restrict(
        lambda wave: _norm_squared(wave) == child_norm_squared
    )
    parent_flux = -l2_inner(parent_nonlinear, parent)
    child_flux = -l2_inner(child_nonlinear, child)
    parent_energy = Fraction(parent.sobolev_sq(0)) / 2
    child_energy = Fraction(child.sobolev_sq(0)) / 2
    parent_viscous = nu * Fraction(parent.sobolev_sq(1))
    child_viscous = nu * Fraction(child.sobolev_sq(1))

    populated = set(field.coefficient_table())
    off_chain = nonlinear.restrict(lambda wave: wave not in populated)
    difference = (0, scale, -scale)
    # coefficient_table canonicalises signs internally, so direct membership of
    # either orientation is the robust exact check.
    difference_present = difference in nonlinear.coefficient_table() or tuple(
        -value for value in difference
    ) in nonlinear.coefficient_table()

    expected_child = -scale * b * c_coefficient / 3
    child_wave = (2 * scale, scale, scale)
    child_pair = nonlinear.coefficient_table().get(child_wave)
    expected_vector = (expected_child, -expected_child, -expected_child)
    child_matches = child_pair is not None and tuple(child_pair[0]) == expected_vector and all(
        value == 0 for value in child_pair[1]
    )
    formula_flux = scale * b * c_coefficient * d / 2
    if child_flux != formula_flux:
        raise AssertionError("exact child flux disagrees with the closed formula")
    formula_off = (
        Fraction(3, 8)
        * scale
        * scale
        * d
        * d
        * (b * b + c_coefficient * c_coefficient)
    )
    off_norm = Fraction(off_chain.sobolev_sq(0))
    if off_norm != formula_off:
        raise AssertionError("exact off-chain norm disagrees with the closed formula")

    return ExactRelayMetrics(
        scale=scale,
        parent_energy=parent_energy,
        child_energy=child_energy,
        child_flux=child_flux,
        parent_flux=parent_flux,
        child_viscous_loss=child_viscous,
        parent_viscous_loss=parent_viscous,
        child_net=child_flux - child_viscous,
        off_chain_nonlinear_l2_squared=off_norm,
        formula_off_chain_l2_squared=formula_off,
        total_nonlinear_energy_defect=l2_inner(nonlinear, field),
        difference_mode_present=difference_present,
        child_coefficient_matches_formula=child_matches,
    )


def fixed_cardinality_scaling() -> dict[str, str]:
    r"""Return the obstruction exponents for a critically weighted fixed triad.

    At wavenumber ``N``, a fixed number of modes with total energy ``N^-1``
    has amplitude ``N^-1/2``.  The cubic transfer is then
    ``N * (N^-1/2)^3 = N^-1/2``, whereas child viscosity is
    ``N^2 * N^-1 = N``.  Their ratio is ``N^-3/2`` and tends to zero.
    """

    return {
        "mode_count": "N^0",
        "energy": "N^-1",
        "coefficient_amplitude": "N^-1/2",
        "signed_flux": "N^-1/2",
        "viscous_loss": "N^1",
        "flux_to_viscosity_ratio": "N^-3/2",
        "necessary_mode_count": "(nu^2/c_E) * N^3 up to fixed constants",
        "verdict": "REJECTED as a fixed-cardinality scale iteration",
    }
