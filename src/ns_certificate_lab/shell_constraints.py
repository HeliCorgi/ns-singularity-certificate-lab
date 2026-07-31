r"""Historical positive-slope conditions for a diverging-bandwidth candidate.

The fixed-finite-mode no-go (``docs/research_notes/track_f_finite_mode_nogo.md``)
closes exactly one search class: velocity fields whose Fourier support stays
inside a *fixed finite* set.  Its contrapositive says that any surviving Track-F
candidate must have **diverging bandwidth**, ``N(t) → ∞`` as ``t → T``.  This
module implements the exponent bookkeeping for that surviving class.

Derivation and citations: ``docs/research_notes/track_f_shell_constraints.md``.
The short version.  Model the dyadic shell energies near ``T`` by

    E_j(t) = A(t) (lambda_j / N(t))^beta   for lambda_j <= N(t),  0 above,
    N(t) = (T-t)^(-gamma),   A(t) = (T-t)^sigma,   beta > 0,

with ``lambda_j = 2^j``.  Four necessary conditions then become inequalities on
``(gamma, sigma, beta)``:

* **bounded energy** — ``sum_j E_j`` converges and stays bounded, so
  ``beta > 0`` and ``sigma >= 0``;
* **finite total dissipation** — ``int_0^T sum_j lambda_j^2 E_j dt < infinity``
  requires ``sigma - 2 gamma > -1``;
* **critical-norm divergence** — a bounded ``B^{1/2}_{2,1}`` norm bounds
  ``L^3`` and the endpoint regularity theorem then forbids a singularity, so
  ``sigma < gamma``;
* **diverging bandwidth** — ``gamma > 0``, the contrapositive of the
  fixed-finite-mode no-go.

Together these force ``0 < gamma < 1`` and ``max(0, 2 gamma - 1) < sigma <
gamma``.  Everything outside that open triangle is excluded.

Correction (2026-08-01).  On the torus the lowest nonzero shell is fixed, so
the low-frequency sum has only ``J+1`` terms when ``N=2**J``.  Consequently
``beta > 0`` is not necessary in the full finite-floor model.  This module is
retained as the historical **positive-slope branch** and for artifact
compatibility.  The corrected piecewise classifier, including the surviving
``beta=-1, sigma=gamma`` logarithmic boundary, is implemented in
``zeno_packet_relay.classify_finite_floor_shell_exponents`` and derived in
``track_f_shell_constraints_finite_floor_erratum.md``.

Scope.  Within ``beta > 0`` these are necessary conditions on an assumed power law, not an
existence statement: no candidate is constructed, and a point inside the
feasible triangle is not a solution of anything.  The critical-norm condition
quotes Escauriaza--Seregin--Sverak; it is neither proved nor formalized here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

__all__ = [
    "ShellExponents",
    "ShellVerdict",
    "feasible_sigma_interval",
    "feasible_sigma_lower_is_attained",
    "scan_shell_exponents",
    "classify_shell_exponents",
    "MAXIMUM_BANDWIDTH_EXPONENT",
]

#: ``gamma`` must be strictly below this value; see the note, Proposition 3.
MAXIMUM_BANDWIDTH_EXPONENT = 1.0


@dataclass(frozen=True)
class ShellExponents:
    """One power-law point ``(gamma, sigma, beta)`` of the shell ansatz."""

    bandwidth: float
    peak_energy: float
    spectral_slope: float

    def __post_init__(self) -> None:
        for name, value in (
            ("bandwidth", self.bandwidth),
            ("peak_energy", self.peak_energy),
            ("spectral_slope", self.spectral_slope),
        ):
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")

    # -- the four necessary conditions, each as a named predicate ---------- #

    @property
    def bandwidth_diverges(self) -> bool:
        """``gamma > 0``: the contrapositive of the fixed-finite-mode no-go."""
        return self.bandwidth > 0.0

    @property
    def spectrum_summable(self) -> bool:
        """Whether this point belongs to the historical positive-slope branch.

        This is not a necessary condition in the corrected finite-floor torus
        model.  It is preserved for compatibility with preregistered scans.
        """
        return self.spectral_slope > 0.0

    @property
    def energy_bounded(self) -> bool:
        """``sigma >= 0``: the total energy does not diverge as ``t -> T``."""
        return self.peak_energy >= 0.0

    @property
    def dissipation_integrable(self) -> bool:
        """``sigma - 2 gamma > -1``: the total dissipation is finite."""
        return self.peak_energy - 2.0 * self.bandwidth > -1.0

    @property
    def critical_norm_diverges(self) -> bool:
        """``sigma < gamma``: the critical Besov norm actually blows up."""
        return self.peak_energy < self.bandwidth

    @property
    def sobolev_exponent_threshold(self) -> float:
        """Smallest ``s`` with ``||u||_{H^s} -> infinity``, namely ``sigma/(2 gamma)``.

        Useful as a design number: it says which Sobolev norm a candidate would
        first have to be tracked in.  Meaningless when ``gamma <= 0``.
        """
        if self.bandwidth <= 0.0:
            return math.inf
        return self.peak_energy / (2.0 * self.bandwidth)

    @property
    def dissipation_time_exponent(self) -> float:
        """``sigma - 2 gamma``; the enstrophy behaves like ``(T-t)`` to this power."""
        return self.peak_energy - 2.0 * self.bandwidth

    @property
    def critical_norm_time_exponent(self) -> float:
        """``(sigma - gamma)/2``; the ``B^{1/2}_{2,1}`` norm grows at this rate."""
        return 0.5 * (self.peak_energy - self.bandwidth)


@dataclass(frozen=True)
class ShellVerdict:
    """Classification of one exponent point."""

    exponents: ShellExponents
    failed_conditions: tuple[str, ...]

    @property
    def feasible(self) -> bool:
        return not self.failed_conditions

    @property
    def verdict(self) -> str:
        return "conditionally_admissible" if self.feasible else "excluded"

    def as_dict(self) -> dict[str, object]:
        return {
            "bandwidth": self.exponents.bandwidth,
            "peak_energy": self.exponents.peak_energy,
            "spectral_slope": self.exponents.spectral_slope,
            "failed_conditions": list(self.failed_conditions),
            "feasible": self.feasible,
            "verdict": self.verdict,
            "dissipation_time_exponent": self.exponents.dissipation_time_exponent,
            "critical_norm_time_exponent": (
                self.exponents.critical_norm_time_exponent
            ),
            "sobolev_exponent_threshold": (
                self.exponents.sobolev_exponent_threshold
            ),
        }


def classify_shell_exponents(exponents: ShellExponents) -> ShellVerdict:
    """Apply the four necessary conditions and name every one that fails."""
    failed: list[str] = []
    if not exponents.bandwidth_diverges:
        failed.append("bandwidth_diverges")
    if not exponents.spectrum_summable:
        failed.append("spectrum_summable")
    if not exponents.energy_bounded:
        failed.append("energy_bounded")
    if not exponents.dissipation_integrable:
        failed.append("dissipation_integrable")
    if not exponents.critical_norm_diverges:
        failed.append("critical_norm_diverges")
    return ShellVerdict(exponents=exponents, failed_conditions=tuple(failed))


def feasible_sigma_interval(bandwidth: float) -> tuple[float, float] | None:
    """The open interval of admissible ``sigma`` for a given ``gamma``.

    Returns ``None`` when no ``sigma`` survives, which happens exactly for
    ``gamma <= 0`` and for ``gamma >= 1``.  The upper endpoint is ``gamma``
    (critical-norm divergence) and the lower endpoint is
    ``max(0, 2 gamma - 1)`` (bounded energy and finite dissipation).

    The upper endpoint is never attained.  The lower endpoint **is** attained
    exactly when ``gamma < 1/2``, because there it comes from the non-strict
    bounded-energy condition ``sigma >= 0`` rather than from the strict
    dissipation condition; see :func:`feasible_sigma_lower_is_attained`.
    """
    if bandwidth <= 0.0 or bandwidth >= MAXIMUM_BANDWIDTH_EXPONENT:
        return None
    lower = max(0.0, 2.0 * bandwidth - 1.0)
    if lower >= bandwidth:
        return None
    return (lower, bandwidth)


def feasible_sigma_lower_is_attained(bandwidth: float) -> bool:
    """Whether the lower endpoint of :func:`feasible_sigma_interval` is admissible.

    It is, exactly when ``2 gamma - 1 < 0``: the endpoint is then ``sigma = 0``,
    which satisfies the non-strict energy bound and the strict dissipation
    bound.  For ``gamma >= 1/2`` the endpoint is ``2 gamma - 1`` itself and the
    dissipation condition excludes it.
    """
    interval = feasible_sigma_interval(bandwidth)
    if interval is None:
        return False
    return 2.0 * bandwidth - 1.0 < 0.0


def scan_shell_exponents(
    bandwidths: Sequence[float],
    peak_energies: Sequence[float],
    spectral_slopes: Sequence[float],
) -> list[ShellVerdict]:
    """Classify a product grid of exponent points."""
    verdicts: list[ShellVerdict] = []
    for bandwidth in bandwidths:
        for peak_energy in peak_energies:
            for slope in spectral_slopes:
                verdicts.append(
                    classify_shell_exponents(
                        ShellExponents(
                            bandwidth=float(bandwidth),
                            peak_energy=float(peak_energy),
                            spectral_slope=float(slope),
                        )
                    )
                )
    return verdicts


def smooth_force_shell_bound(
    wavenumber: float,
    order: int,
    constant: float,
) -> float:
    r"""``||f_j||_{L^2} <= constant * lambda_j^{-order}`` for a smooth force.

    A Clay-admissible force is ``C^infinity`` on ``T^3 x [0,infinity)``, so on
    the compact time interval ``[0,T]`` its Fourier coefficients decay faster
    than any polynomial.  The consequence recorded in the note is that the force
    contributes ``O(lambda_j^{-m})`` to the shell-``j`` energy budget for every
    ``m``: **the force cannot pump the cascade**, only reshape the ``O(1)``
    shells.  A Track-F candidate with diverging bandwidth therefore has to build
    its high-wavenumber energy by the same nonlinear flux an unforced candidate
    would.
    """
    if wavenumber <= 0.0:
        raise ValueError("wavenumber must be positive")
    if order < 0:
        raise ValueError("order must be nonnegative")
    if constant < 0.0:
        raise ValueError("constant must be nonnegative")
    return constant * wavenumber ** (-float(order))
