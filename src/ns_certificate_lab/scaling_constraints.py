"""Conditional power-law constraints for Navier--Stokes singularity candidates.

This module studies a deliberately narrow asymptotic model.  Let
``tau = T - t`` and suppose a localized core has velocity amplitude
``A(t) ~ tau**(-alpha)`` and axisymmetric length scales
``L_r(t) ~ tau**beta_r`` and ``L_z(t) ~ tau**beta_z``.  Write
``B = 2*beta_r + beta_z`` for the core-volume exponent and
``gamma = max(beta_r, beta_z)`` for the thinnest represented length scale.

The checks below are necessary or conditional checks, not a proof of blow-up.
They assume a nondegenerate profile with nonzero L2/L3 mass and, for the
viscous-dissipation estimate, a nonzero derivative in a direction with scale
``tau**gamma``.  The dominant-balance classifier also assumes no unmeasured
leading cancellation between convection and pressure.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import math
from typing import Iterable


class BalanceFamily(str, Enum):
    """Leading-order term balance in the one-effective-gradient model."""

    TIME_INERTIAL = "time_inertial"
    TIME_VISCOUS = "time_viscous"
    INERTIAL_VISCOUS = "inertial_viscous"
    TRIPLE_PARABOLIC = "triple_parabolic"
    REQUIRES_LEADING_CANCELLATION = "requires_leading_cancellation"


class CandidateFamily(str, Enum):
    """Generic finite-energy, critical-growth families surviving the checks."""

    EULER_LIKE_TYPE_II = "euler_like_type_ii"
    ANISOTROPIC_PARABOLIC = "anisotropic_parabolic"
    QUASI_STEADY_VISCOUS_INERTIAL = "quasi_steady_viscous_inertial"
    NONE = "none"


@dataclass(frozen=True)
class PowerLawScaling:
    """Power-law exponents for one localized axisymmetric core."""

    alpha: float
    beta_r: float
    beta_z: float

    def __post_init__(self) -> None:
        values = (self.alpha, self.beta_r, self.beta_z)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("all scaling exponents must be finite")
        if self.alpha < 0.0 or self.beta_r <= 0.0 or self.beta_z <= 0.0:
            raise ValueError("alpha must be nonnegative and length exponents positive")

    @property
    def volume_exponent(self) -> float:
        """Return ``B = 2 beta_r + beta_z``."""

        return 2.0 * self.beta_r + self.beta_z

    @property
    def thinnest_scale_exponent(self) -> float:
        """Return ``gamma = max(beta_r, beta_z)``.

        For ``tau -> 0``, the larger positive exponent is the smaller physical
        length scale.
        """

        return max(self.beta_r, self.beta_z)

    def lp_growth_exponent(self, p: float) -> float:
        """Return ``g_p`` in ``||u||_p ~ tau**(-g_p)`` for the core model."""

        p_value = float(p)
        if not math.isfinite(p_value) or p_value <= 0.0:
            raise ValueError("p must be positive and finite")
        return self.alpha - self.volume_exponent / p_value

    def serrin_time_exponent(self, p: float) -> float:
        """Return the time singularity exponent in the critical Serrin norm.

        For ``p > 3``, the scale-critical time exponent is
        ``q = 2p/(p-3)``.  The integral of ``||u||_p**q`` is locally finite
        when this returned exponent is strictly less than one.
        """

        p_value = float(p)
        if not math.isfinite(p_value) or p_value <= 3.0:
            raise ValueError("Serrin p must be finite and greater than three")
        q_value = 2.0 * p_value / (p_value - 3.0)
        return q_value * self.lp_growth_exponent(p_value)


@dataclass(frozen=True)
class ScalingConstraintReport:
    """Margins and classifications for one power-law candidate."""

    scaling: PowerLawScaling
    volume_exponent: float
    thinnest_scale_exponent: float
    energy_margin: float
    critical_l3_margin: float
    serrin_infinity_margin: float
    dissipation_margin: float
    time_term_exponent: float
    inertial_term_exponent: float
    viscous_term_exponent: float
    balance_family: BalanceFamily
    candidate_family: CandidateFamily
    finite_energy_compatible: bool
    critical_l3_divergent: bool
    escapes_classical_serrin_scaling: bool
    finite_dissipation_compatible: bool
    generic_balance_compatible: bool
    passes_conditional_core_gates: bool

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable dictionary."""

        result = asdict(self)
        result["balance_family"] = self.balance_family.value
        result["candidate_family"] = self.candidate_family.value
        return result


def _leading_balance(
    time_exponent: float,
    inertial_exponent: float,
    viscous_exponent: float,
    *,
    tolerance: float,
) -> BalanceFamily:
    exponents = {
        "time": time_exponent,
        "inertial": inertial_exponent,
        "viscous": viscous_exponent,
    }
    maximum = max(exponents.values())
    leaders = {
        name for name, value in exponents.items() if maximum - value <= tolerance
    }
    if leaders == {"time", "inertial", "viscous"}:
        return BalanceFamily.TRIPLE_PARABOLIC
    if leaders == {"time", "inertial"}:
        return BalanceFamily.TIME_INERTIAL
    if leaders == {"time", "viscous"}:
        return BalanceFamily.TIME_VISCOUS
    if leaders == {"inertial", "viscous"}:
        return BalanceFamily.INERTIAL_VISCOUS
    return BalanceFamily.REQUIRES_LEADING_CANCELLATION


def analyze_power_law_scaling(
    scaling: PowerLawScaling,
    *,
    tolerance: float = 1.0e-10,
) -> ScalingConstraintReport:
    """Evaluate conditional necessary constraints for a localized core.

    Margin conventions:

    - ``energy_margin = B - 2 alpha``; nonnegative is compatible with bounded
      core energy.
    - ``critical_l3_margin = 3 alpha - B``; positive means the modeled core
      L3 norm diverges.
    - ``serrin_infinity_margin = 2 alpha - 1``; nonnegative is necessary to
      avoid the ``L^2_t L^infty_x`` member of the classical Serrin class.
    - ``dissipation_margin = 1 - (2 alpha + 2 gamma - B)``; positive is
      compatible with finite energy dissipation under the nondegenerate
      thinnest-direction gradient assumption.
    """

    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("tolerance must be finite and nonnegative")

    alpha = scaling.alpha
    volume = scaling.volume_exponent
    gamma = scaling.thinnest_scale_exponent

    energy_margin = volume - 2.0 * alpha
    critical_l3_margin = 3.0 * alpha - volume
    serrin_infinity_margin = 2.0 * alpha - 1.0
    dissipation_power = 2.0 * alpha + 2.0 * gamma - volume
    dissipation_margin = 1.0 - dissipation_power

    time_exponent = alpha + 1.0
    inertial_exponent = 2.0 * alpha + gamma
    viscous_exponent = alpha + 2.0 * gamma
    balance = _leading_balance(
        time_exponent,
        inertial_exponent,
        viscous_exponent,
        tolerance=tolerance,
    )

    finite_energy = energy_margin >= -tolerance
    critical_l3 = critical_l3_margin > tolerance
    serrin_escape = serrin_infinity_margin >= -tolerance and critical_l3_margin >= -tolerance
    finite_dissipation = dissipation_margin > tolerance
    generic_balance = balance is not BalanceFamily.REQUIRES_LEADING_CANCELLATION

    family = CandidateFamily.NONE
    if finite_energy and critical_l3 and serrin_escape and finite_dissipation:
        if balance is BalanceFamily.TIME_INERTIAL and gamma < 0.5 - tolerance:
            family = CandidateFamily.EULER_LIKE_TYPE_II
        elif balance is BalanceFamily.TRIPLE_PARABOLIC:
            family = CandidateFamily.ANISOTROPIC_PARABOLIC
        elif balance is BalanceFamily.INERTIAL_VISCOUS and gamma > 0.5 + tolerance:
            family = CandidateFamily.QUASI_STEADY_VISCOUS_INERTIAL

    passes = (
        finite_energy
        and critical_l3
        and serrin_escape
        and finite_dissipation
        and family is not CandidateFamily.NONE
    )

    return ScalingConstraintReport(
        scaling=scaling,
        volume_exponent=volume,
        thinnest_scale_exponent=gamma,
        energy_margin=energy_margin,
        critical_l3_margin=critical_l3_margin,
        serrin_infinity_margin=serrin_infinity_margin,
        dissipation_margin=dissipation_margin,
        time_term_exponent=time_exponent,
        inertial_term_exponent=inertial_exponent,
        viscous_term_exponent=viscous_exponent,
        balance_family=balance,
        candidate_family=family,
        finite_energy_compatible=finite_energy,
        critical_l3_divergent=critical_l3,
        escapes_classical_serrin_scaling=serrin_escape,
        finite_dissipation_compatible=finite_dissipation,
        generic_balance_compatible=generic_balance,
        passes_conditional_core_gates=passes,
    )


def isotropic_type_ii_wedge_contains(
    alpha: float,
    beta: float,
    *,
    tolerance: float = 1.0e-10,
) -> bool:
    """Return whether exponents lie in the generic isotropic Type-II wedge.

    Under the module assumptions the wedge is

    ``alpha + beta = 1`` and ``2/5 <= beta < 1/2``.
    """

    values = (float(alpha), float(beta), float(tolerance))
    if not all(math.isfinite(value) for value in values) or tolerance < 0.0:
        raise ValueError("alpha, beta and tolerance must be finite")
    return (
        abs(alpha + beta - 1.0) <= tolerance
        and beta >= 0.4 - tolerance
        and beta < 0.5 - tolerance
    )


def generic_family_bounds(gamma: float) -> dict[str, float | str | bool]:
    """Return analytic volume-exponent bounds for a generic balance family.

    The bounds include finite core energy, divergent core L3, finite total
    dissipation, and geometry ``B <= 3 gamma``.  Open endpoints are identified
    explicitly in the returned dictionary.
    """

    value = float(gamma)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("gamma must be positive and finite")

    if 0.4 <= value < 0.5:
        return {
            "family": CandidateFamily.EULER_LIKE_TYPE_II.value,
            "alpha": 1.0 - value,
            "volume_lower": 2.0 - 2.0 * value,
            "volume_upper": 3.0 * value,
            "lower_open": False,
            "upper_open": False,
            "feasible": True,
        }
    if math.isclose(value, 0.5, rel_tol=0.0, abs_tol=1.0e-12):
        return {
            "family": CandidateFamily.ANISOTROPIC_PARABOLIC.value,
            "alpha": 0.5,
            "volume_lower": 1.0,
            "volume_upper": 1.5,
            "lower_open": True,
            "upper_open": True,
            "feasible": True,
        }
    if 0.5 < value < 1.0:
        return {
            "family": CandidateFamily.QUASI_STEADY_VISCOUS_INERTIAL.value,
            "alpha": value,
            "volume_lower": 4.0 * value - 1.0,
            "volume_upper": 3.0 * value,
            "lower_open": True,
            "upper_open": True,
            "feasible": True,
        }
    return {
        "family": CandidateFamily.NONE.value,
        "alpha": math.nan,
        "volume_lower": math.nan,
        "volume_upper": math.nan,
        "lower_open": True,
        "upper_open": True,
        "feasible": False,
    }


def scan_power_law_scalings(
    alphas: Iterable[float],
    beta_rs: Iterable[float],
    beta_zs: Iterable[float],
    *,
    tolerance: float = 1.0e-10,
) -> list[ScalingConstraintReport]:
    """Return all grid points that pass the conditional generic core gates."""

    reports: list[ScalingConstraintReport] = []
    for alpha in alphas:
        for beta_r in beta_rs:
            for beta_z in beta_zs:
                report = analyze_power_law_scaling(
                    PowerLawScaling(float(alpha), float(beta_r), float(beta_z)),
                    tolerance=tolerance,
                )
                if report.passes_conditional_core_gates:
                    reports.append(report)
    return reports
