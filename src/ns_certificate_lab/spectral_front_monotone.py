r"""Exact spectral front-gap ledger for the closable bandwidth monotone.

For the unforced periodic Navier--Stokes equation in repository conventions,
with ``e_k`` the paired real modal energy, ``a_k`` the modal growth inner
product, and (new here) ``n_k`` the paired real modal energy of the Leray
nonlinearity ``N = -P((u.grad)u)``, define the moments

.. math:: H_r=\sum x_k^re_k,\quad T_r=\sum x_k^ra_k,\quad G_r=\sum x_k^rn_k,
          \qquad x_k=|k|^2 .

The exact front identity of ``modal_front_actions`` bounds the covariance by
``sqrt(V_r Var(g))``, which is sharper but *not closable*: ``Var(g)`` is
dominated by spectral near-nodes and is bounded by no Sobolev norm.  The
closable route applies the modal Cauchy--Schwarz ``a_k^2 <= e_k n_k``
*before* dividing by ``e_k``:

.. math:: |\mathrm{Cov}_{p_r}(x,g)|\le\sqrt{V_rG_r/H_r},

and, after a square completion in ``sqrt(V_r)``,

.. math:: \frac{d}{dt}\log N_r^2\ \le\ \frac{G_r}{2\nu H_{r+1}} .

The gap between the two sides decomposes into two nonnegative pieces (a
Cauchy--Schwarz gap and a square-completion gap); their *sum* is exactly
rational even though each piece separately contains a square root, so this
module certifies the decomposition through the rational total together with
the squared Cauchy--Schwarz certificate ``Cov^2 <= V_r G_r/H_r``.

Companion unconditional bound (Lemma K): the nonlinear front wavenumber
``K = \|P(u.grad u)\|_2^2/\|\nabla u\|_2^4`` of a band-limited field obeys
``K <= S_N`` with ``S_N=\sum_{0<|k|^2\le N^2}|k|^{-2}`` (an exactly rational
lattice sum ``~ 4 pi N``), with no phase or coherence assumption.

Everything in the exact lane is `fractions.Fraction` arithmetic on finite
trigonometric fields; the float lane provides the amplitude-delocalization
diagnostic ``mu_N = M_eff/N^3`` for sparse mesoscopic parents.  This module
is a regularity research tool; nothing in it is a PDE regularity proof.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Mapping

import numpy as np
import numpy.typing as npt

from .fourier_torus import TrigVector, advection, leray

__all__ = [
    "FrontGapIdentity",
    "front_gap_identity",
    "front_wavenumber",
    "inverse_square_lattice_sum",
    "lemma_k_certificate",
    "sparse_parent_delocalization",
]


def _norm_squared(wave: tuple[int, int, int]) -> int:
    return sum(component * component for component in wave)


def _pair_energy(pair) -> Fraction:
    cosine, sine = pair
    return Fraction(1, 2) * sum(
        Fraction(value) * Fraction(value)
        for vector in (cosine, sine)
        for value in vector
    )


def _pair_inner(left, right) -> Fraction:
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


def _ledger(field: TrigVector, maximum_order: int):
    """Return modal ``(e_k, a_k, n_k)`` and moments ``(H, T, G)`` exactly."""

    if field.divergence().cleaned().terms:
        raise ValueError("field must be exactly divergence free")
    table = field.coefficient_table()
    mean_pair = table.get((0, 0, 0))
    if mean_pair is not None and any(
        Fraction(value) != 0 for vector in mean_pair for value in vector
    ):
        raise ValueError("field must have zero spatial mean")
    nonlinear = -leray(advection(field, field)).cleaned()
    growth_table = nonlinear.coefficient_table()
    zero_pair = ([Fraction(0)] * 3, [Fraction(0)] * 3)
    energies: dict[tuple[int, int, int], Fraction] = {}
    growth: dict[tuple[int, int, int], Fraction] = {}
    power: dict[tuple[int, int, int], Fraction] = {}
    for wave, pair in table.items():
        if wave == (0, 0, 0):
            continue
        energy = _pair_energy(pair)
        if energy <= 0:
            continue
        growth_pair = growth_table.get(wave, zero_pair)
        energies[wave] = energy
        growth[wave] = _pair_inner(pair, growth_pair)
        power[wave] = _pair_energy(growth_pair)
        # Modal Cauchy--Schwarz certificate |a_k|^2 <= e_k n_k, exact.
        if growth[wave] ** 2 > energies[wave] * power[wave]:
            raise AssertionError("the modal Cauchy--Schwarz certificate failed")
    h = [Fraction(0) for _ in range(maximum_order + 1)]
    t = [Fraction(0) for _ in range(maximum_order + 1)]
    g = [Fraction(0) for _ in range(maximum_order + 1)]
    for wave, energy in energies.items():
        x = Fraction(_norm_squared(wave))
        for order in range(maximum_order + 1):
            h[order] += x**order * energy
            t[order] += x**order * growth[wave]
            g[order] += x**order * power[wave]
    return energies, growth, power, h, t, g


@dataclass(frozen=True)
class FrontGapIdentity:
    """Exact rational front-gap ledger at one moment order."""

    order: int
    viscosity: Fraction
    h_r: Fraction
    h_r1: Fraction
    h_r2: Fraction
    t_r: Fraction
    t_r1: Fraction
    g_r: Fraction
    bandwidth_squared: Fraction
    covariance: Fraction
    variance_direct: Fraction
    log_bandwidth_derivative: Fraction
    closable_upper: Fraction
    gap_total: Fraction
    saturation_deficit: Fraction
    cauchy_schwarz_certificate: Fraction

    def as_dict(self) -> dict[str, int | str]:
        return {
            "order": self.order,
            "viscosity": str(self.viscosity),
            "h_r": str(self.h_r),
            "h_r1": str(self.h_r1),
            "h_r2": str(self.h_r2),
            "t_r": str(self.t_r),
            "t_r1": str(self.t_r1),
            "g_r": str(self.g_r),
            "bandwidth_squared": str(self.bandwidth_squared),
            "covariance": str(self.covariance),
            "variance_direct": str(self.variance_direct),
            "log_bandwidth_derivative": str(self.log_bandwidth_derivative),
            "closable_upper": str(self.closable_upper),
            "gap_total": str(self.gap_total),
            "saturation_deficit": str(self.saturation_deficit),
            "cauchy_schwarz_certificate": str(self.cauchy_schwarz_certificate),
        }


def front_gap_identity(
    field: TrigVector,
    *,
    order: int,
    viscosity: Fraction | int,
) -> FrontGapIdentity:
    """Evaluate the closable front-gap ledger in exact rational arithmetic."""

    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    nu = Fraction(viscosity)
    if nu <= 0:
        raise ValueError("viscosity must be positive")
    energies, growth, _power, h, t, g = _ledger(field, order + 2)
    h_r, h_r1, h_r2 = h[order], h[order + 1], h[order + 2]
    t_r, t_r1 = t[order], t[order + 1]
    g_r = g[order]
    if h_r <= 0 or h_r1 <= 0:
        raise ValueError("the requested modal moments must be positive")
    mu = h_r1 / h_r
    covariance = t_r1 / h_r - mu * (t_r / h_r)
    # Direct variance sum (never the moment difference): TM-09 discipline,
    # although in exact arithmetic both agree and the agreement is asserted.
    variance_direct = Fraction(0)
    for wave, energy in energies.items():
        x = Fraction(_norm_squared(wave))
        variance_direct += x**order * energy * (x - mu) ** 2
    variance_direct /= h_r
    if variance_direct != h_r2 / h_r - mu * mu:
        raise AssertionError("exact variance bookkeeping disagreement")
    log_derivative = 2 * (covariance / mu - nu * variance_direct / mu)
    closable_upper = g_r / (2 * nu * h_r1)
    gap_total = closable_upper - log_derivative
    # (I.4) as a rational identity: the two square roots cancel in the sum.
    identity_rhs = (
        -2 * covariance / mu
        + 2 * nu * variance_direct / mu
        + g_r / (2 * nu * mu * h_r)
    )
    if gap_total != identity_rhs:
        raise AssertionError("the exact gap decomposition identity failed")
    # Cauchy--Schwarz certificate for the closable covariance bound.
    certificate = variance_direct * g_r / h_r - covariance * covariance
    if covariance > 0 and certificate < 0:
        raise AssertionError("the closable covariance certificate failed")
    if gap_total < 0:
        raise AssertionError("the closable upper bound (I.3) failed")
    if closable_upper > 0:
        deficit = gap_total / closable_upper
    else:
        deficit = Fraction(0)
    return FrontGapIdentity(
        order=order,
        viscosity=nu,
        h_r=h_r,
        h_r1=h_r1,
        h_r2=h_r2,
        t_r=t_r,
        t_r1=t_r1,
        g_r=g_r,
        bandwidth_squared=mu,
        covariance=covariance,
        variance_direct=variance_direct,
        log_bandwidth_derivative=log_derivative,
        closable_upper=closable_upper,
        gap_total=gap_total,
        saturation_deficit=deficit,
        cauchy_schwarz_certificate=certificate,
    )


def front_wavenumber(field: TrigVector) -> Fraction:
    """Return ``K = G_0/(2 H_1^2)``, the nonlinear front wavenumber."""

    _, _, _, h, _, g = _ledger(field, 1)
    if h[1] <= 0:
        raise ValueError("the field must have positive enstrophy")
    return g[0] / (2 * h[1] * h[1])


def inverse_square_lattice_sum(maximum_norm_squared: int) -> Fraction:
    """Return ``sum_{0<|k|^2<=X} |k|^{-2}`` over the integer lattice, exactly."""

    if (
        isinstance(maximum_norm_squared, bool)
        or not isinstance(maximum_norm_squared, int)
        or maximum_norm_squared < 1
    ):
        raise ValueError("maximum_norm_squared must be a positive integer")
    bound = math.isqrt(maximum_norm_squared)
    axis = np.arange(-bound, bound + 1, dtype=np.int64)
    kx, ky, kz = np.meshgrid(axis, axis, axis, indexing="ij")
    norms = (kx * kx + ky * ky + kz * kz).ravel()
    counts = np.bincount(norms, minlength=maximum_norm_squared + 1)
    total = Fraction(0)
    for value in range(1, maximum_norm_squared + 1):
        if counts[value]:
            total += Fraction(int(counts[value]), value)
    return total


def lemma_k_certificate(field: TrigVector) -> dict[str, str | int]:
    """Certify ``K <= S_N`` exactly for a band-limited field."""

    table = field.coefficient_table()
    occupied = [
        _norm_squared(wave)
        for wave, pair in table.items()
        if wave != (0, 0, 0) and _pair_energy(pair) > 0
    ]
    if not occupied:
        raise ValueError("the field has no occupied modes")
    x_max = max(occupied)
    lattice_sum = inverse_square_lattice_sum(x_max)
    k_value = front_wavenumber(field)
    margin = lattice_sum - k_value
    if margin < 0:
        raise AssertionError("the Lemma K lattice certificate failed")
    return {
        "band_norm_squared": x_max,
        "front_wavenumber_K": str(k_value),
        "lattice_sum_S": str(lattice_sum),
        "margin": str(margin),
        "four_pi_n_float": 4.0 * math.pi * math.sqrt(float(x_max)),
    }


def sparse_parent_delocalization(
    parent: Mapping[tuple[int, int, int], npt.ArrayLike],
    *,
    scale: int,
) -> dict[str, float]:
    """Return the delocalization diagnostics of a sparse Hermitian parent.

    ``parent`` maps positive-representative wavevectors to complex coefficient
    triples; the Hermitian negatives are implied.  ``M_eff`` is the amplitude
    (Rényi-1/2) effective mode count over the *full* (both-sign) support and
    ``mu_N = M_eff/scale^3`` is the delocalization ratio whose positive floor
    the mesoscopic no-go requires.  Binary64 diagnostic.
    """

    if isinstance(scale, bool) or not isinstance(scale, int) or scale < 1:
        raise ValueError("scale must be a positive integer")
    amplitudes = []
    for wave, coefficient in parent.items():
        if wave == (0, 0, 0):
            raise ValueError("sparse parents must not carry a mean mode")
        vector = np.asarray(coefficient, dtype=np.complex128)
        if vector.shape != (3,):
            raise ValueError("sparse coefficients must be three-vectors")
        amplitudes.append(float(np.linalg.norm(vector)))
    if not amplitudes or not any(value > 0.0 for value in amplitudes):
        raise ValueError("the sparse parent is empty")
    amplitude_sum = 2.0 * float(np.sum(amplitudes))
    energy_sum = 2.0 * float(np.sum(np.square(amplitudes)))
    support = 2 * sum(1 for value in amplitudes if value > 0.0)
    effective = amplitude_sum**2 / energy_sum
    return {
        "support_mode_count": float(support),
        "effective_mode_count": effective,
        "effective_over_support": effective / support,
        "mu_delocalization": effective / float(scale**3),
    }
