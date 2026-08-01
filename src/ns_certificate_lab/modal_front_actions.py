r"""Exact finite-Fourier checks for modal front identities.

For the unforced periodic Navier--Stokes equation write

.. math:: \partial_t\hat u_k=\mathcal N_k-\nu |k|^2\hat u_k,

where ``N=-P(u.grad u)``.  With ``e_k`` the paired real Fourier ``L2``
contribution and ``a_k=<u_k,N_k>``, define

.. math:: H_r=\sum x_k^r e_k,\quad T_r=\sum x_k^r a_k,\quad x_k=|k|^2.

Then ``(1/2) Hdot_r=T_r-nu H_{r+1}``.  The radial bandwidth
``N_r^2=H_{r+1}/H_r`` obeys the exact covariance identity described by
:func:`modal_growth_identity`.  The identity is a regularity research tool,
not a proof that its action is integrable for arbitrary PDE solutions.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .fourier_torus import TrigVector, advection, leray

__all__ = [
    "ModalGrowthIdentity",
    "h3_bandwidth_factorization",
    "modal_growth_identity",
]


def _norm_squared(wave: tuple[int, int, int]) -> int:
    return sum(component * component for component in wave)


def _mode_energy(pair) -> Fraction:
    cosine, sine = pair
    return Fraction(1, 2) * sum(
        Fraction(value) * Fraction(value) for vector in (cosine, sine) for value in vector
    )


def _mode_inner(left, right) -> Fraction:
    left_cosine, left_sine = left
    right_cosine, right_sine = right
    return Fraction(1, 2) * sum(
        Fraction(a) * Fraction(b)
        for first, second in (
            (left_cosine, right_cosine),
            (left_sine, right_sine),
        )
        for a, b in zip(first, second)
    )


def _moments(field: TrigVector, maximum_order: int):
    if maximum_order < 0:
        raise ValueError("maximum_order must be nonnegative")
    if field.divergence().cleaned().terms:
        raise ValueError("field must be exactly divergence free")
    table = field.coefficient_table()
    mean_pair = table.get((0, 0, 0))
    if mean_pair is not None and any(
        Fraction(value) != 0 for vector in mean_pair for value in vector
    ):
        raise ValueError("field must have zero spatial mean")
    nonlinear_growth = -leray(advection(field, field)).cleaned()
    growth_table = nonlinear_growth.coefficient_table()
    energies: dict[tuple[int, int, int], Fraction] = {}
    growth: dict[tuple[int, int, int], Fraction] = {}
    for wave, pair in table.items():
        if wave == (0, 0, 0):
            continue
        energy = _mode_energy(pair)
        if energy <= 0:
            continue
        energies[wave] = energy
        growth_pair = growth_table.get(
            wave,
            ([Fraction(0)] * 3, [Fraction(0)] * 3),
        )
        growth[wave] = _mode_inner(pair, growth_pair)
    h = [Fraction(0) for _ in range(maximum_order + 1)]
    t = [Fraction(0) for _ in range(maximum_order + 1)]
    for wave, energy in energies.items():
        x = Fraction(_norm_squared(wave))
        for order in range(maximum_order + 1):
            h[order] += x**order * energy
            t[order] += x**order * growth[wave]
    return energies, growth, h, t


@dataclass(frozen=True)
class ModalGrowthIdentity:
    """Exact rational components of the modal covariance diagnostic."""

    order: int
    h_r: Fraction
    h_r1: Fraction
    h_r2: Fraction
    t_r: Fraction
    t_r1: Fraction
    bandwidth_squared: Fraction
    spectral_variance_relative: Fraction
    modal_growth_variance: Fraction
    covariance: Fraction
    covariance_front_term: Fraction
    direct_front_term: Fraction
    viscous_front_term: Fraction
    log_bandwidth_derivative: Fraction
    cauchy_schwarz_gap: Fraction
    square_completion_upper: Fraction
    square_completion_gap: Fraction

    def as_dict(self) -> dict[str, int | str]:
        return {
            "order": self.order,
            "h_r": str(self.h_r),
            "h_r1": str(self.h_r1),
            "h_r2": str(self.h_r2),
            "t_r": str(self.t_r),
            "t_r1": str(self.t_r1),
            "bandwidth_squared": str(self.bandwidth_squared),
            "spectral_variance_relative": str(self.spectral_variance_relative),
            "modal_growth_variance": str(self.modal_growth_variance),
            "covariance": str(self.covariance),
            "covariance_front_term": str(self.covariance_front_term),
            "direct_front_term": str(self.direct_front_term),
            "viscous_front_term": str(self.viscous_front_term),
            "log_bandwidth_derivative": str(self.log_bandwidth_derivative),
            "cauchy_schwarz_gap": str(self.cauchy_schwarz_gap),
            "square_completion_upper": str(self.square_completion_upper),
            "square_completion_gap": str(self.square_completion_gap),
        }


def modal_growth_identity(
    field: TrigVector,
    *,
    order: int,
    viscosity: Fraction | int,
) -> ModalGrowthIdentity:
    """Evaluate the covariance identity in exact rational arithmetic."""

    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    nu = Fraction(viscosity)
    if nu <= 0:
        raise ValueError("viscosity must be positive")
    energies, growth, h, t = _moments(field, order + 2)
    h_r, h_r1, h_r2 = h[order], h[order + 1], h[order + 2]
    t_r, t_r1 = t[order], t[order + 1]
    if h_r <= 0 or h_r1 <= 0:
        raise ValueError("the requested modal moments must be positive")
    mu = h_r1 / h_r
    mean_growth = t_r / h_r
    covariance = t_r1 / h_r - mu * mean_growth
    front_covariance = covariance / mu
    direct_front = t_r1 / h_r1 - t_r / h_r
    if front_covariance != direct_front:
        raise AssertionError("the exact covariance front identity failed")
    variance_x = h_r2 / h_r - mu * mu
    delta = variance_x / (mu * mu)
    variance_growth = Fraction(0)
    for wave, energy in energies.items():
        x = Fraction(_norm_squared(wave))
        probability = x**order * energy / h_r
        relative_growth = growth[wave] / energy
        variance_growth += probability * (relative_growth - mean_growth) ** 2
    cauchy_gap = variance_x * variance_growth - covariance * covariance
    if cauchy_gap < 0:
        raise AssertionError("the exact covariance Cauchy--Schwarz gap is negative")
    viscous = nu * mu * delta
    log_derivative = direct_front - viscous
    upper = variance_growth / (4 * nu * mu)
    square_gap = upper - log_derivative
    if square_gap < 0:
        raise AssertionError("the square-completion upper bound failed")
    return ModalGrowthIdentity(
        order=order,
        h_r=h_r,
        h_r1=h_r1,
        h_r2=h_r2,
        t_r=t_r,
        t_r1=t_r1,
        bandwidth_squared=mu,
        spectral_variance_relative=delta,
        modal_growth_variance=variance_growth,
        covariance=covariance,
        covariance_front_term=front_covariance,
        direct_front_term=direct_front,
        viscous_front_term=viscous,
        log_bandwidth_derivative=log_derivative,
        cauchy_schwarz_gap=cauchy_gap,
        square_completion_upper=upper,
        square_completion_gap=square_gap,
    )


def h3_bandwidth_factorization(field: TrigVector) -> dict[str, str]:
    r"""Check ``H3=H0 N0^2 N1^2 N2^2`` exactly for a finite field."""

    _, _, h, _ = _moments(field, 3)
    if any(value <= 0 for value in h):
        raise ValueError("H0 through H3 must be positive")
    bandwidths = [h[index + 1] / h[index] for index in range(3)]
    reconstructed = h[0]
    for value in bandwidths:
        reconstructed *= value
    if reconstructed != h[3]:
        raise AssertionError("the H3 bandwidth factorization failed")
    return {
        "H0": str(h[0]),
        "H1": str(h[1]),
        "H2": str(h[2]),
        "H3": str(h[3]),
        "N0_squared": str(bandwidths[0]),
        "N1_squared": str(bandwidths[1]),
        "N2_squared": str(bandwidths[2]),
        "reconstructed_H3": str(reconstructed),
    }
