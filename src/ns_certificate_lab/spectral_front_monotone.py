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

**Convention warning (verification sprint V1, workstream C).** :func:`_ledger`
iterates over the *field's* coefficient table, so the ``G_r`` it returns — and
therefore the ``G_r``, ``K`` and deficit of :func:`front_gap_identity` and
:func:`front_wavenumber` — is the **in-support** moment
``G_r^{in}=\sum_{k:e_k>0}x_k^rn_k``.  That is legitimate and in fact sharper in
(I.2) (``a_k=0`` off the support of ``u``), but it is *not* the published
``G_0=\tfrac12\|\mathbb P(u\cdot\nabla u)\|_2^2`` that Lemma K bounds and that
``KD\le\|u\|_\infty^2`` refers to.  Use :func:`full_nonlinear_power` or
``front_defect_decomposition(..., convention="full")`` for the published
object; the difference is the *leakage* ``G_r^{out}``.  See
``docs/research_notes/verification_sprint_v1/lambda_O9_defect_decomposition.md``
§4, which also records that the deficit reported by
:func:`front_gap_identity` is minimised over ``nu`` at ``nu_*=Cov/(2V_r)`` with
the viscosity-free value ``1-Cov^2H_r/(V_rG_r)`` — the pilot's 0.865/0.740 are
artefacts of ``nu=1/40``, not saturation measurements.

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
from .snapshot_certificate import Interval

__all__ = [
    "FrontDefectDecomposition",
    "FrontGapIdentity",
    "front_defect_decomposition",
    "front_gap_identity",
    "front_wavenumber",
    "full_nonlinear_power",
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


def full_nonlinear_power(field: TrigVector, maximum_order: int) -> list[Fraction]:
    r"""Return ``G_r=\sum_k x_k^r n_k`` over **all** modes of ``N``, exactly.

    :func:`_ledger` deliberately restricts every moment to the support of the
    field, so its ``G_r`` is the *sharper* in-support moment
    ``G_r^{in}=\sum_{k:e_k>0}x_k^rn_k`` (legitimate in (I.2), because
    ``a_k=0`` off the support of ``u``).  The published monotone instead uses
    the full ``G_0=\tfrac12\|\mathbb P(u\cdot\nabla u)\|_2^2`` — that is the
    quantity Lemma K bounds and the one the Serrin domination
    ``KD\le\|u\|_\infty^2`` refers to.  The two differ by the *leakage*
    ``G_r^{out}``, the nonlinear power that the convolution deposits on modes
    the field does not occupy.  Both are returned by callers so the
    convention is never silently swapped.
    """

    if maximum_order < 0:
        raise ValueError("maximum_order must be nonnegative")
    nonlinear = -leray(advection(field, field)).cleaned()
    totals = [Fraction(0) for _ in range(maximum_order + 1)]
    for wave, pair in nonlinear.coefficient_table().items():
        if wave == (0, 0, 0):
            continue
        power = _pair_energy(pair)
        if power == 0:
            continue
        x = Fraction(_norm_squared(wave))
        for order in range(maximum_order + 1):
            totals[order] += x**order * power
    return totals


@dataclass(frozen=True)
class FrontDefectDecomposition:
    """Exact defect decomposition of the chain (I.1) -> (I.2) -> (I.3).

    Every field is a string of an exact ``Fraction`` unless the name ends in
    ``_lower`` / ``_upper``, which are the endpoints of a rigorous rational
    enclosure of an irrational (square-root) quantity.
    """

    order: int
    viscosity: Fraction
    convention: str
    h_r: Fraction
    h_r1: Fraction
    h_r2: Fraction
    t_r: Fraction
    t_r1: Fraction
    g_in: Fraction
    g_full: Fraction
    g_r: Fraction
    bandwidth_squared: Fraction
    variance: Fraction
    covariance: Fraction
    absolute_covariance: Fraction
    delta_sign: Fraction
    rational_remainder: Fraction
    gap_total: Fraction
    closable_upper: Fraction
    saturation_deficit: Fraction
    optimal_viscosity: Fraction | None
    optimal_deficit: Fraction
    delta_cs_modal: tuple[Fraction, Fraction]
    delta_cs_vector: tuple[Fraction, Fraction]
    delta_square_completion: tuple[Fraction, Fraction]
    split_residual: tuple[Fraction, Fraction]
    per_mode_cs_defect: dict[tuple[int, int, int], Fraction]

    def as_dict(self) -> dict[str, object]:
        def pair(value: tuple[Fraction, Fraction]) -> list[str]:
            return [str(value[0]), str(value[1])]

        return {
            "order": self.order,
            "viscosity": str(self.viscosity),
            "convention": self.convention,
            "h_r": str(self.h_r),
            "h_r1": str(self.h_r1),
            "h_r2": str(self.h_r2),
            "t_r": str(self.t_r),
            "t_r1": str(self.t_r1),
            "g_in_support": str(self.g_in),
            "g_full": str(self.g_full),
            "g_r_used": str(self.g_r),
            "bandwidth_squared": str(self.bandwidth_squared),
            "variance": str(self.variance),
            "covariance": str(self.covariance),
            "absolute_covariance": str(self.absolute_covariance),
            "delta_sign": str(self.delta_sign),
            "rational_remainder": str(self.rational_remainder),
            "gap_total": str(self.gap_total),
            "closable_upper": str(self.closable_upper),
            "saturation_deficit": str(self.saturation_deficit),
            "optimal_viscosity": (
                None if self.optimal_viscosity is None
                else str(self.optimal_viscosity)
            ),
            "optimal_deficit": str(self.optimal_deficit),
            "delta_cs_modal_enclosure": pair(self.delta_cs_modal),
            "delta_cs_vector_enclosure": pair(self.delta_cs_vector),
            "delta_square_completion_enclosure": pair(self.delta_square_completion),
            "split_residual_enclosure": pair(self.split_residual),
            "per_mode_cs_defect": {
                str(list(wave)): str(value)
                for wave, value in sorted(self.per_mode_cs_defect.items())
            },
        }


def _sqrt_enclosure(value: Fraction, *, bits: int) -> tuple[Fraction, Fraction]:
    from .l3_certificate import sqrt_interval

    if value < 0:
        raise ValueError("cannot enclose the square root of a negative rational")
    result = sqrt_interval(Interval(value, value), bits=bits)
    return (result.lower, result.upper)


def front_defect_decomposition(
    field: TrigVector,
    *,
    order: int,
    viscosity: Fraction | int,
    convention: str = "full",
    bits: int = 96,
) -> FrontDefectDecomposition:
    r"""Split the (I.3) gap into four explicitly nonnegative defects.

    With ``w_k=x_k^r|x_k-\mu|/H_r``, ``C=\mathrm{Cov}_{p_r}(x,g)``,
    ``A=\sum_kw_k|a_k|``, ``B=\sum_kw_k\sqrt{e_kn_k}`` and
    ``S=\sqrt{V_rG_r/H_r}`` the chain ``C\le A\le B\le S`` gives

    ``\hat\Delta_{sign}=(2/\mu)(A-C)``,
    ``\hat\Delta_{CS,modal}=(2/\mu)(B-A)``,
    ``\hat\Delta_{CS,vector}=(2/\mu)(S-B)``,
    ``\hat\Delta_{SC}=(2\nu/\mu)(\sqrt{V_r}-\sqrt{G_r/H_r}/(2\nu))^2``,

    and the telescoping identity

    ``G_r/(2\nu H_{r+1})-\tfrac{d}{dt}\log N_r^2
      =\hat\Delta_{sign}+\hat\Delta_{CS,modal}+\hat\Delta_{CS,vector}
       +\hat\Delta_{SC}``.

    ``\hat\Delta_{sign}`` and the *sum* of the last three,
    ``R=(2/\mu)[\nu V_r+G_r/(4\nu H_r)-A]``, are exactly rational and are
    checked as exact ``Fraction`` equalities; the individual last three carry
    square roots and are returned as rigorous rational enclosures whose sum is
    verified to contain ``R``.

    ``convention`` selects ``G_r``: ``"full"`` (the published monotone's
    ``G_0=\tfrac12\|\mathbb P(u\cdot\nabla u)\|_2^2``) or ``"in_support"``
    (the sharper moment used by :func:`front_gap_identity`).
    """

    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    if convention not in ("full", "in_support"):
        raise ValueError("convention must be 'full' or 'in_support'")
    nu = Fraction(viscosity)
    if nu <= 0:
        raise ValueError("viscosity must be positive")
    energies, growth, power, h, t, g = _ledger(field, order + 2)
    h_r, h_r1, h_r2 = h[order], h[order + 1], h[order + 2]
    if h_r <= 0 or h_r1 <= 0:
        raise ValueError("the requested modal moments must be positive")
    g_in = g[order]
    g_full = full_nonlinear_power(field, order)[order]
    if g_full < g_in:
        raise AssertionError("the full nonlinear power is below its restriction")
    g_r = g_full if convention == "full" else g_in

    mu = h_r1 / h_r
    covariance = t[order + 1] / h_r - mu * (t[order] / h_r)
    variance = Fraction(0)
    absolute = Fraction(0)
    defects: dict[tuple[int, int, int], Fraction] = {}
    modal_terms: list[tuple[Fraction, Fraction]] = []
    for wave, energy in energies.items():
        x = Fraction(_norm_squared(wave))
        weight = x**order * abs(x - mu) / h_r
        variance += x**order * energy * (x - mu) ** 2
        absolute += weight * abs(growth[wave])
        defects[wave] = energy * power[wave] - growth[wave] ** 2
        if defects[wave] < 0:
            raise AssertionError("the modal Cauchy--Schwarz certificate failed")
        modal_terms.append((weight, energy * power[wave]))
    variance /= h_r
    if variance != h_r2 / h_r - mu * mu:
        raise AssertionError("exact variance bookkeeping disagreement")

    log_derivative = 2 * (covariance / mu - nu * variance / mu)
    closable_upper = g_r / (2 * nu * h_r1)
    gap_total = closable_upper - log_derivative
    delta_sign = 2 * (absolute - covariance) / mu
    remainder = 2 * (nu * variance + g_r / (4 * nu * h_r) - absolute) / mu
    if delta_sign < 0:
        raise AssertionError("the sign defect is negative")
    if remainder < 0:
        raise AssertionError("the rational remainder defect is negative")
    if gap_total != delta_sign + remainder:
        raise AssertionError("the exact defect telescoping failed")

    # Enclosures for the three radical defects.
    b_lower = Fraction(0)
    b_upper = Fraction(0)
    for weight, product in modal_terms:
        low, high = _sqrt_enclosure(product, bits=bits)
        b_lower += weight * low
        b_upper += weight * high
    s_lower, s_upper = _sqrt_enclosure(variance * g_r / h_r, bits=bits)
    v_lower, v_upper = _sqrt_enclosure(variance, bits=bits)
    q_lower, q_upper = _sqrt_enclosure(g_r / h_r, bits=bits)
    modal = (2 * (b_lower - absolute) / mu, 2 * (b_upper - absolute) / mu)
    vector = (2 * (s_lower - b_upper) / mu, 2 * (s_upper - b_lower) / mu)
    low_gap = v_lower - q_upper / (2 * nu)
    high_gap = v_upper - q_lower / (2 * nu)
    squares = sorted(value * value for value in (low_gap, high_gap))
    inner = Fraction(0) if low_gap <= 0 <= high_gap else squares[0]
    completion = (2 * nu * inner / mu, 2 * nu * squares[1] / mu)
    split = (
        modal[0] + vector[0] + completion[0] - remainder,
        modal[1] + vector[1] + completion[1] - remainder,
    )
    if not (split[0] <= 0 <= split[1]):
        raise AssertionError("the radical defect split does not enclose R")
    if modal[1] < 0 or vector[1] < 0 or completion[1] < 0:
        raise AssertionError("a radical defect enclosure is strictly negative")

    deficit = gap_total / closable_upper if closable_upper > 0 else Fraction(0)
    if covariance > 0 and variance > 0:
        optimal_viscosity = covariance / (2 * variance)
        optimal_deficit = 1 - covariance * covariance * h_r / (variance * g_r)
    else:
        optimal_viscosity = None
        optimal_deficit = Fraction(1)
    if optimal_deficit < 0:
        raise AssertionError("the viscosity-optimised deficit is negative")
    return FrontDefectDecomposition(
        order=order,
        viscosity=nu,
        convention=convention,
        h_r=h_r,
        h_r1=h_r1,
        h_r2=h_r2,
        t_r=t[order],
        t_r1=t[order + 1],
        g_in=g_in,
        g_full=g_full,
        g_r=g_r,
        bandwidth_squared=mu,
        variance=variance,
        covariance=covariance,
        absolute_covariance=absolute,
        delta_sign=delta_sign,
        rational_remainder=remainder,
        gap_total=gap_total,
        closable_upper=closable_upper,
        saturation_deficit=deficit,
        optimal_viscosity=optimal_viscosity,
        optimal_deficit=optimal_deficit,
        delta_cs_modal=modal,
        delta_cs_vector=vector,
        delta_square_completion=completion,
        split_residual=split,
        per_mode_cs_defect=defects,
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
