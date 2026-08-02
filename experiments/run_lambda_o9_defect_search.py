r"""Lambda O-9: defect decomposition and adversarial minimisation of the
saturation deficit of the closable front bound (I.3).

Verification sprint V1, workstream C.  Three lanes, never mixed:

* **exact lane** (`fractions.Fraction`) -- the four-defect telescoping
  ``G_r/(2 nu H_{r+1}) - d/dt log N_r^2
   = Delta_sign + Delta_CS_modal + Delta_CS_vector + Delta_SC``
  on the repository's exact rational fields, the closed-form relay family,
  and the rationalised near-minimisers of the float search.  Every number
  reported from this lane is an exact rational (or a rigorous rational
  enclosure of a square root).
* **float lane** (binary64, numpy only) -- adversarial minimisation of the
  viscosity-optimised deficit over divergence-free zero-mean fields
  band-limited to ``|k|_inf <= B``, by analytic-gradient ascent from random
  and structured multi-starts.  Float minimisers are *never* reported as
  bounds; they are rationalised and re-evaluated in the exact lane, and only
  those exact values are certified upper bounds on the infimum.
* **closed form** -- the three-mode Leray relay family, whose deficit is
  derived in closed rational form and verified exactly against the ledger.

Key exact reduction used throughout (derived in the companion note): for a
fixed field the deficit is minimised over ``nu`` at
``nu_*=Cov/(2V_r)`` with value

    d_*(u) = 1 - Cov^2 H_r/(V_r G_r),

a viscosity-free, amplitude-free, scale-invariant rational functional.  The
rigidity lemma of obligation O-8/O-9 is exactly the statement
``sup_u Cov^2H_r/(V_rG_r) < 1``.

Research diagnostic.  Nothing here is a PDE regularity proof.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from ns_certificate_lab._integrity import canonical_json_bytes, write_with_digest
from ns_certificate_lab.exact_leray_relay import build_exact_relay_triad
from ns_certificate_lab.fourier_torus import (
    TrigVector,
    family_P1,
    family_P2,
    family_P3,
)
from ns_certificate_lab.spectral_front_monotone import (
    front_defect_decomposition,
    front_gap_identity,
    full_nonlinear_power,
)

OUTPUT_SCHEMA = "ns-certificate-lab/lambda-o9-defect-search/v1"
STATUS = (
    "EXACT RATIONAL DEFECT TELESCOPING + BINARY64 ADVERSARIAL SEARCH / "
    "NOT A PROOF"
)


# --------------------------------------------------------------------------- #
# float lane: band-limited divergence-free parameterisation                    #
# --------------------------------------------------------------------------- #


class Band:
    r"""Divergence-free zero-mean fields with ``1 <= |k|_inf <= band``.

    Each canonical wavevector ``k`` carries four real parameters, the
    components of the cosine and sine coefficient vectors in the *integer*
    orthogonal basis ``t1 = k x e_j`` (first ``j`` with ``k x e_j != 0``) and
    ``t2 = k x t1`` of ``k^\perp``.  Divergence-freeness is therefore
    structural, and a rational parameter vector maps to an exactly
    divergence-free rational field.
    """

    def __init__(self, band: int) -> None:
        if band < 1:
            raise ValueError("band must be at least one")
        self.band = band
        grid = 4 * band + 1
        self.grid = grid
        axis = np.fft.fftfreq(grid, d=1.0 / grid).astype(np.int64)
        kx, ky, kz = np.meshgrid(axis, axis, axis, indexing="ij")
        self.k = np.stack([kx, ky, kz]).astype(np.float64)
        self.ksq = (kx * kx + ky * ky + kz * kz).astype(np.float64)
        self.inv_ksq = 1.0 / np.where(self.ksq == 0.0, 1.0, self.ksq)
        modes = []
        for a in range(-band, band + 1):
            for b in range(-band, band + 1):
                for c in range(-band, band + 1):
                    wave = (a, b, c)
                    if wave == (0, 0, 0):
                        continue
                    for component in wave:
                        if component > 0:
                            modes.append(wave)
                            break
                        if component < 0:
                            break
        self.modes = sorted(modes)
        self.basis: list[tuple[np.ndarray, np.ndarray, float, float]] = []
        for wave in self.modes:
            vector = np.array(wave, dtype=np.int64)
            first = None
            for j in range(3):
                unit = np.zeros(3, dtype=np.int64)
                unit[j] = 1
                candidate = np.cross(vector, unit)
                if np.any(candidate):
                    first = candidate
                    break
            second = np.cross(vector, first)
            self.basis.append(
                (
                    first,
                    second,
                    float(np.linalg.norm(first)),
                    float(np.linalg.norm(second)),
                )
            )
        self.parameter_count = 4 * len(self.modes)
        self.index = [tuple(int(c) % grid for c in w) for w in self.modes]
        self.conjugate = [tuple((-int(c)) % grid for c in w) for w in self.modes]
        self.in_band = np.zeros((grid, grid, grid), dtype=bool)
        for positive, negative in zip(self.index, self.conjugate):
            self.in_band[positive] = True
            self.in_band[negative] = True

    # -- transforms ---------------------------------------------------------- #

    def coefficients(self, theta: np.ndarray) -> np.ndarray:
        grid = self.grid
        hat = np.zeros((3, grid, grid, grid), dtype=np.complex128)
        table = theta.reshape(len(self.modes), 4)
        for m, (first, second, n1, n2) in enumerate(self.basis):
            cosine = table[m, 0] * first / n1 + table[m, 1] * second / n2
            sine = table[m, 2] * first / n1 + table[m, 3] * second / n2
            positive, negative = self.index[m], self.conjugate[m]
            for component in range(3):
                hat[component][positive] += (
                    cosine[component] - 1j * sine[component]
                ) / 2
                hat[component][negative] += (
                    cosine[component] + 1j * sine[component]
                ) / 2
        return hat

    def to_parameters(self, hat: np.ndarray) -> np.ndarray:
        """``d/dtheta`` of ``<F, u(theta)>`` for the real field with FT ``hat``."""

        out = np.zeros((len(self.modes), 4))
        for m, (first, second, n1, n2) in enumerate(self.basis):
            positive = self.index[m]
            value = np.array([hat[c][positive] for c in range(3)])
            real, imaginary = np.real(value), np.imag(value)
            out[m, 0] = real @ first / n1
            out[m, 1] = real @ second / n2
            out[m, 2] = -(imaginary @ first) / n1
            out[m, 3] = -(imaginary @ second) / n2
        return out.ravel()

    def physical(self, hat: np.ndarray) -> np.ndarray:
        return np.fft.ifftn(hat, axes=(-3, -2, -1)) * (self.grid**3)

    def spectral(self, field: np.ndarray) -> np.ndarray:
        return np.fft.fftn(field, axes=(-3, -2, -1)) / (self.grid**3)

    def leray(self, hat: np.ndarray) -> np.ndarray:
        radial = sum(self.k[i] * hat[i] for i in range(3))
        return np.array(
            [hat[i] - self.k[i] * radial * self.inv_ksq for i in range(3)]
        )

    def nonlinearity(self, hat: np.ndarray) -> np.ndarray:
        field = self.physical(hat).real
        gradient = self.physical(1j * self.k[:, None] * hat[None, :]).real
        product = np.einsum("ixyz,ijxyz->jxyz", field, gradient)
        return -self.leray(self.spectral(product))


@dataclass
class FloatRecord:
    """Binary64 diagnostics of one field.  Never a certified bound."""

    h0: float
    h1: float
    h2: float
    t0: float
    t1: float
    g_full: float
    g_band: float

    @property
    def mu(self) -> float:
        return self.h1 / self.h0

    @property
    def spread(self) -> float:
        return self.h0 * self.h2 / (self.h1 * self.h1) - 1.0

    def deficit(self, *, full: bool) -> float:
        weight = self.h2 - self.h1 * self.h1 / self.h0
        power = self.g_full if full else self.g_band
        if weight <= 0.0 or power <= 0.0 or self.t1 <= 0.0:
            return 1.0
        return max(0.0, 1.0 - self.t1 * self.t1 / (weight * power))


def measure(band: Band, hat: np.ndarray, nonlinear: np.ndarray) -> FloatRecord:
    energy = np.sum(np.abs(hat) ** 2, axis=0)
    power = np.sum(np.abs(nonlinear) ** 2, axis=0)
    growth = np.sum(np.real(hat * np.conj(nonlinear)), axis=0)
    return FloatRecord(
        h0=float(np.sum(energy)),
        h1=float(np.sum(band.ksq * energy)),
        h2=float(np.sum(band.ksq**2 * energy)),
        t0=float(np.sum(growth)),
        t1=float(np.sum(band.ksq * growth)),
        g_full=float(np.sum(power)),
        g_band=float(np.sum(power[band.in_band])),
    )


def ratio_and_gradient(
    band: Band, theta: np.ndarray, *, full: bool
) -> tuple[float, np.ndarray, FloatRecord]:
    r"""Return ``J = T_1^2/((H_2-H_1^2/H_0)G_0)``, its gradient, and the ledger.

    ``J = 1 - d_*`` is the viscosity-optimised saturation ratio; the gradient
    is analytic (validated against central differences in the test suite).
    """

    hat = band.coefficients(theta)
    nonlinear = band.nonlinearity(hat)
    record = measure(band, hat, nonlinear)
    effective = nonlinear if full else np.where(band.in_band, nonlinear, 0.0)
    power = record.g_full if full else record.g_band
    weight = record.h2 - record.h1 * record.h1 / record.h0
    if weight <= 0.0 or power <= 0.0:
        return 0.0, np.zeros_like(theta), record
    ratio = record.t1 * record.t1 / (weight * power)
    field = band.physical(hat).real
    gradient_u = band.physical(1j * band.k[:, None] * hat[None, :]).real

    def adjoint(test_hat: np.ndarray) -> np.ndarray:
        """``grad_u <N(u), psi>`` for divergence-free ``psi`` (physical form)."""

        test = band.physical(test_hat).real
        gradient_psi = band.physical(
            1j * band.k[:, None] * test_hat[None, :]
        ).real
        first = -np.einsum("mjxyz,jxyz->mxyz", gradient_u, test)
        second = np.einsum("ixyz,imxyz->mxyz", field, gradient_psi)
        return first + second

    d_t1 = band.spectral(
        band.physical(band.ksq * nonlinear).real
    ) + band.spectral(adjoint(band.ksq * hat))
    d_power = 2.0 * band.spectral(adjoint(effective))
    d_weight = 2.0 * (band.ksq - record.mu) ** 2 * hat
    d_ratio = (2 * record.t1 / (weight * power)) * d_t1 - (
        record.t1 * record.t1 / (weight * weight * power * power)
    ) * (power * d_weight + weight * d_power)
    return ratio, band.to_parameters(d_ratio), record


def spread_gradient(band: Band, theta: np.ndarray) -> np.ndarray:
    """Gradient of ``V_0/mu^2 = H_0H_2/H_1^2 - 1`` in the parameters."""

    hat = band.coefficients(theta)
    record = measure(band, hat, band.nonlinearity(hat))
    d_hat = (
        (record.h2 * 2.0 * hat + record.h0 * 2.0 * band.ksq**2 * hat)
        / (record.h1**2)
        - 2.0 * record.h0 * record.h2 * (2.0 * band.ksq * hat) / (record.h1**3)
    )
    return band.to_parameters(d_hat)


def ascend(
    band: Band,
    theta: np.ndarray,
    *,
    full: bool,
    steps: int,
    trace: list[tuple[float, float]] | None = None,
) -> tuple[np.ndarray, float]:
    """Projected gradient ascent of ``J`` on the unit sphere of parameters."""

    theta = theta / np.linalg.norm(theta)
    ratio, gradient, record = ratio_and_gradient(band, theta, full=full)
    if trace is not None:
        trace.append((record.spread, 1.0 - ratio))
    step = 0.5
    for _ in range(steps):
        projected = gradient - (gradient @ theta) * theta
        norm = float(np.linalg.norm(projected))
        if norm < 1e-15:
            break
        direction = projected / norm
        improved = False
        for _ in range(45):
            candidate = theta + step * direction
            candidate /= np.linalg.norm(candidate)
            value, new_gradient, new_record = ratio_and_gradient(
                band, candidate, full=full
            )
            if value > ratio:
                theta, ratio, gradient, record = (
                    candidate,
                    value,
                    new_gradient,
                    new_record,
                )
                step *= 1.7
                improved = True
                break
            step *= 0.4
        if trace is not None:
            trace.append((record.spread, 1.0 - ratio))
        if not improved:
            break
    return theta, ratio


def descend_constrained(
    band: Band,
    theta: np.ndarray,
    *,
    target: float,
    steps: int,
    penalty: float = 5.0,
) -> np.ndarray:
    r"""Minimise ``d_*`` under a soft constraint ``V_0/mu^2 = target``.

    The penalised objective ``-J + penalty*(log(V_0/mu^2)-log target)^2`` is
    well scaled because ``J in [0,1]``; the constraint residual achieved is
    reported alongside every row, so a poorly enforced constraint is visible
    rather than hidden.
    """

    theta = theta / np.linalg.norm(theta)

    def objective(vector: np.ndarray) -> tuple[float, np.ndarray]:
        ratio, gradient, record = ratio_and_gradient(band, vector, full=True)
        spread = max(record.spread, 1e-300)
        mismatch = math.log(spread) - math.log(target)
        value = -ratio + penalty * mismatch * mismatch
        grad = -gradient + (
            2.0 * penalty * mismatch / spread
        ) * spread_gradient(band, vector)
        return value, grad

    value, gradient = objective(theta)
    step = 0.3
    for _ in range(steps):
        projected = gradient - (gradient @ theta) * theta
        norm = float(np.linalg.norm(projected))
        if norm < 1e-15:
            break
        direction = -projected / norm
        improved = False
        for _ in range(40):
            candidate = theta + step * direction
            candidate /= np.linalg.norm(candidate)
            new_value, new_gradient = objective(candidate)
            if new_value < value:
                theta, value, gradient = candidate, new_value, new_gradient
                step *= 1.7
                improved = True
                break
            step *= 0.4
        if not improved:
            break
    return theta


def relay_seed(band: Band, child: float) -> np.ndarray | None:
    r"""The exact relay triad embedded in the band parameterisation.

    ``u = e_3 sin(p.x) + e_2 cos(q.x) + child*(1,-1,-1)cos((p+q).x)`` with
    ``p=(1,1,0)``, ``q=(1,0,1)``: a known near-saturating configuration whose
    deficit is available in closed form, used as a structured multi-start.
    """

    wanted = {
        (1, 1, 0): (np.zeros(3), np.array([0.0, 0.0, 1.0])),
        (1, 0, 1): (np.array([0.0, 1.0, 0.0]), np.zeros(3)),
        (2, 1, 1): (child * np.array([1.0, -1.0, -1.0]), np.zeros(3)),
    }
    theta = np.zeros((len(band.modes), 4))
    for wave, (cosine, sine) in wanted.items():
        if wave not in band.modes:
            return None
        m = band.modes.index(wave)
        first, second, n1, n2 = band.basis[m]
        u1, u2 = first / n1, second / n2
        theta[m, 0] = cosine @ u1
        theta[m, 1] = cosine @ u2
        theta[m, 2] = sine @ u1
        theta[m, 3] = sine @ u2
    return theta.ravel()


# --------------------------------------------------------------------------- #
# exact lane                                                                   #
# --------------------------------------------------------------------------- #


def rationalise(
    band: Band, theta: np.ndarray, *, keep: int, bits: int
) -> TrigVector | None:
    """Sparsify and dyadically rationalise a float parameter vector.

    Dyadic rounding (a common denominator ``2^bits``) keeps the exact
    arithmetic of the downstream ledger tractable; ``limit_denominator``
    produces coprime denominators whose products blow up.
    """

    table = theta.reshape(len(band.modes), 4)
    strength = np.sum(table * table, axis=1)
    order = np.argsort(strength)[::-1][:keep]
    scale = float(np.max(np.abs(table[order]))) or 1.0
    unit = Fraction(1, 2**bits)
    modes = []
    for m in order:
        first, second, n1, n2 = band.basis[m]
        lambdas = [
            unit
            * round(
                float(table[m, j]) / scale / (n1 if j % 2 == 0 else n2) * 2**bits
            )
            for j in range(4)
        ]
        cosine = tuple(
            lambdas[0] * int(first[i]) + lambdas[1] * int(second[i])
            for i in range(3)
        )
        sine = tuple(
            lambdas[2] * int(first[i]) + lambdas[3] * int(second[i])
            for i in range(3)
        )
        if all(value == 0 for value in cosine + sine):
            continue
        modes.append((band.modes[m], cosine, sine))
    if not modes:
        return None
    return TrigVector.from_modes(modes)


def exact_certificate(field: TrigVector, *, label: str) -> dict[str, Any] | None:
    """Exact ``d_*`` (viscosity-optimised deficit) for both G conventions.

    ``d_*`` is defined only when ``Cov > 0``; since ``N(-u) = N(u)`` and
    ``a_k(-u) = -a_k(u)``, the orientation with positive covariance is the
    admissible one and is selected here (the float objective ``T_1^2/(WG)`` is
    sign-blind, so a minimiser can come out with either orientation).
    """

    probe = front_defect_decomposition(
        field, order=0, viscosity=Fraction(1), convention="full"
    )
    orientation = "u"
    if probe.covariance < 0:
        field = field.scale(Fraction(-1))
        orientation = "-u"
    out: dict[str, Any] = {"label": label, "orientation": orientation}
    for convention in ("full", "in_support"):
        try:
            record = front_defect_decomposition(
                field, order=0, viscosity=Fraction(1), convention=convention
            )
        except (ValueError, AssertionError) as error:
            out[convention] = {"status": f"unavailable: {error}"}
            continue
        out[convention] = {
            "status": "exact",
            "optimal_viscosity": (
                None
                if record.optimal_viscosity is None
                else str(record.optimal_viscosity)
            ),
            "optimal_deficit": str(record.optimal_deficit),
            "optimal_deficit_float": float(record.optimal_deficit),
            "spectral_spread": str(record.variance / record.bandwidth_squared**2),
            "spectral_spread_float": float(
                record.variance / record.bandwidth_squared**2
            ),
            "g_in_support": str(record.g_in),
            "g_full": str(record.g_full),
            "leakage_fraction": (
                str((record.g_full - record.g_in) / record.g_full)
                if record.g_full > 0
                else None
            ),
        }
    return out


def relay_family_closed_form(
    b: Fraction, c: Fraction, d: Fraction, scale: int
) -> dict[str, str]:
    r"""Closed-form deficit of the exact three-mode Leray relay.

    With ``p=s(1,1,0)``, ``q=s(1,0,1)``, ``child=p+q``, ``n=p x q/s^2`` and
    ``u = B e_3 sin(p.x) + C e_2 cos(q.x) + D n cos(child.x)``:

    ``H_0=(B^2+C^2+3D^2)/2``, ``H_1=s^2(B^2+C^2+9D^2)``,
    ``H_2=2s^4(B^2+C^2+27D^2)``, ``T_0=0``, ``T_1=2s^3BCD``,
    ``n_child=s^2B^2C^2/6``, ``n_p=3s^2C^2D^2/8``, ``n_q=3s^2B^2D^2/8``,
    ``n_{p-q}=0`` (exact cancellation), ``n_{p+child}=3s^2B^2D^2/8``,
    ``n_{q+child}=3s^2C^2D^2/8``, hence with ``P=B^2+C^2``

    ``d_*^{in}   = 1 - B^2C^2(P+3D^2)/(P[B^2C^2+(9/4)D^2P])``,
    ``d_*^{full} = 1 - B^2C^2(P+3D^2)/(P[B^2C^2+(9/2)D^2P])``,

    both independent of ``s`` (exact scale invariance of the r=0 bound), and
    ``V_0/mu^2 = H_0H_2/H_1^2-1 = 12PD^2/(P+9D^2)^2``.
    """

    p = b * b + c * c
    t = d * d
    numerator = b * b * c * c * (p + 3 * t)
    in_support = 1 - numerator / (p * (b * b * c * c + Fraction(9, 4) * t * p))
    full = 1 - numerator / (p * (b * b * c * c + Fraction(9, 2) * t * p))
    spread = 12 * p * t / (p + 9 * t) ** 2
    return {
        "B": str(b),
        "C": str(c),
        "D": str(d),
        "scale": str(scale),
        "deficit_in_support": str(in_support),
        "deficit_in_support_float": repr(float(in_support)),
        "deficit_full": str(full),
        "deficit_full_float": repr(float(full)),
        "spectral_spread": str(spread),
        "spectral_spread_float": repr(float(spread)),
        "ratio_full_over_spread": repr(float(full / spread)),
    }


def equality_conditions(field: TrigVector) -> dict[str, Any]:
    """Extract which saturation conditions hold at a (near) minimiser."""

    from ns_certificate_lab.spectral_front_monotone import _ledger, _norm_squared

    energies, growth, power, h, _t, _g = _ledger(field, 2)
    mu = h[1] / h[0]
    total_weight = Fraction(0)
    aligned_weight = Fraction(0)
    modal_defect = Fraction(0)
    modal_scale = Fraction(0)
    rows = []
    for wave, energy in sorted(energies.items()):
        x = Fraction(_norm_squared(wave))
        weight = abs(x - mu) * abs(growth[wave])
        total_weight += weight
        sign_ok = (x - mu) * growth[wave] >= 0
        if sign_ok:
            aligned_weight += weight
        defect = energy * power[wave] - growth[wave] ** 2
        modal_defect += abs(x - mu) * defect
        modal_scale += abs(x - mu) * energy * power[wave]
        rows.append(
            {
                "wave": list(wave),
                "x_minus_mu": str(x - mu),
                "energy": str(energy),
                "growth": str(growth[wave]),
                "power": str(power[wave]),
                "relative_modal_cs_defect": (
                    repr(float(defect / (energy * power[wave])))
                    if energy * power[wave] > 0
                    else None
                ),
                "sign_aligned": bool(sign_ok),
                "vector_cs_ratio": (
                    repr(float(power[wave] / ((x - mu) ** 2 * energy)))
                    if (x - mu) != 0 and energy > 0
                    else None
                ),
            }
        )
    return {
        "bandwidth_squared": str(mu),
        "sign_aligned_weight_fraction": (
            repr(float(aligned_weight / total_weight))
            if total_weight > 0
            else None
        ),
        "aggregate_relative_modal_cs_defect": (
            repr(float(modal_defect / modal_scale)) if modal_scale > 0 else None
        ),
        "modes": rows,
    }


# --------------------------------------------------------------------------- #
# driver                                                                       #
# --------------------------------------------------------------------------- #


def stage_telescoping() -> list[dict[str, Any]]:
    fields: dict[str, TrigVector] = {
        "relay_triad_s1": build_exact_relay_triad(scale=1),
        "relay_triad_s2": build_exact_relay_triad(scale=2),
        "relay_triad_s3": build_exact_relay_triad(scale=3),
        "relay_triad_s1_d_1_32": build_exact_relay_triad(
            scale=1, child_cosine=Fraction(1, 32)
        ),
        "family_P1": family_P1(),
        "family_P2": family_P2(),
        "family_P3": family_P3(),
    }
    records: list[dict[str, Any]] = []
    for name, field in fields.items():
        entry: dict[str, Any] = {"field": name, "rows": []}
        for convention in ("full", "in_support"):
            for viscosity in (Fraction(1, 40), Fraction(1, 10), Fraction(1)):
                for order in (0, 1, 2):
                    decomposition = front_defect_decomposition(
                        field,
                        order=order,
                        viscosity=viscosity,
                        convention=convention,
                    )
                    row = decomposition.as_dict()
                    if convention == "in_support" and order <= 2:
                        legacy = front_gap_identity(
                            field, order=order, viscosity=viscosity
                        )
                        row["matches_front_gap_identity"] = (
                            legacy.gap_total == decomposition.gap_total
                        )
                    entry["rows"].append(row)
        records.append(entry)
        print(f"telescoping ok: {name}", flush=True)
    return records


def stage_closed_form() -> dict[str, Any]:
    rows = []
    for scale in (1, 2):
        for d in (
            Fraction(1),
            Fraction(1, 2),
            Fraction(1, 8),
            Fraction(1, 32),
            Fraction(1, 128),
            Fraction(1, 512),
        ):
            closed = relay_family_closed_form(Fraction(1), Fraction(1), d, scale)
            field = build_exact_relay_triad(scale=scale, child_cosine=d)
            checks = {}
            for convention, key in (
                ("full", "deficit_full"),
                ("in_support", "deficit_in_support"),
            ):
                record = front_defect_decomposition(
                    field, order=0, viscosity=Fraction(1), convention=convention
                )
                checks[f"ledger_{key}"] = str(record.optimal_deficit)
                checks[f"{key}_matches_ledger"] = str(
                    record.optimal_deficit
                ) == closed[key]
                checks[f"{key}_optimal_viscosity"] = str(record.optimal_viscosity)
            rows.append({**closed, **checks})
    asymmetric = [
        relay_family_closed_form(Fraction(1), Fraction(4), Fraction(1, 64), 1),
        relay_family_closed_form(Fraction(3), Fraction(2), Fraction(1, 64), 1),
    ]
    return {
        "family": (
            "u = B e3 sin(p.x) + C e2 cos(q.x) + D n cos((p+q).x), "
            "p=s(1,1,0), q=s(1,0,1), n=(1,-1,-1)"
        ),
        "rows": rows,
        "asymmetric_amplitudes": asymmetric,
        "limit": (
            "d_full = D^2[(9/2)P^2-3B^2C^2]/(P[B^2C^2+(9/2)D^2P]) -> 0 as D->0; "
            "at B=C the leading terms are d_full=15D^2/P, d_in=6D^2/P and "
            "V_0/mu^2=12D^2/P, so d_full/(V_0/mu^2) -> 5/4 and "
            "d_in/(V_0/mu^2) -> 1/2"
        ),
        "two_mode_classification": stage_two_mode(),
    }


def stage_two_mode() -> dict[str, Any]:
    r"""Exact closed form for every two-mode divergence-free field.

    Let ``supp(u)={\pm k_1,\pm k_2}``.  If ``k_2=m k_1`` the field is
    ``u=f(k_1.x)`` with ``k_1.f=0``, so ``(u.\nabla)u\equiv0`` and every
    moment vanishes.  Otherwise the sums ``l+m`` with ``l,m in supp(u)`` land
    only on ``{0, +-2k_1, +-2k_2, +-(k_1+k_2), +-(k_1-k_2)}``, and
    non-collinearity excludes every one of these from being ``+-k_1`` or
    ``+-k_2``; the self-interaction ``N_{2k}=-i(k.\hat u_k)P_{2k}\hat u_k``
    vanishes by divergence-freeness in any case.  Hence ``a_k\equiv0``,
    ``Cov=T_1=0`` and ``d_*=1`` exactly: **the deficit of every two-mode
    field is maximal.**  A closing triad ``k_1+k_2=k_3`` is the smallest
    structure with a nonzero deficit gradient.
    """

    from ns_certificate_lab.spectral_front_monotone import _ledger

    checks = []
    for label, modes in (
        (
            "non_collinear",
            [((1, 0, 0), (0, 1, 0), (0, 0, 1)), ((0, 1, 0), (1, 0, 0), (0, 0, 2))],
        ),
        (
            "non_collinear_generic",
            [
                ((1, 1, 0), (0, 0, 1), (1, -1, 0)),
                ((2, 0, 1), (0, 1, 0), (1, 0, -2)),
            ],
        ),
        (
            "collinear",
            [((1, 0, 0), (0, 1, 0), (0, 0, 1)), ((2, 0, 0), (0, 1, 1), (0, 2, -1))],
        ),
    ):
        field = TrigVector.from_modes(modes)
        _e, _a, _n, _h, t, g = _ledger(field, 2)
        checks.append(
            {
                "case": label,
                "waves": [list(m[0]) for m in modes],
                "T_0": str(t[0]),
                "T_1": str(t[1]),
                "G_0_in_support": str(g[0]),
                "G_0_full": str(full_nonlinear_power(field, 0)[0]),
                "deficit": "1",
            }
        )
    return {
        "statement": (
            "every two-mode divergence-free zero-mean field has T_1 = 0 and "
            "therefore d_* = 1 (maximal deficit); collinear pairs have "
            "N == 0 identically"
        ),
        "checks": checks,
    }


def stage_search(
    bands: tuple[int, ...], starts: int, steps: int, seed: int
) -> dict[str, Any]:
    results: dict[str, Any] = {"bands": [], "cloud": {}}
    rng = np.random.default_rng(seed)
    best_thetas: dict[int, np.ndarray] = {}
    for size in bands:
        band = Band(size)
        entry: dict[str, Any] = {
            "band_inf_norm": size,
            "canonical_modes": len(band.modes),
            "parameters": band.parameter_count,
        }
        seeds: list[np.ndarray] = []
        for child in (0.5, 0.125, 0.03, 0.008):
            seed_vector = relay_seed(band, child)
            if seed_vector is not None:
                seeds.append(seed_vector)
        previous = max(
            (key for key in best_thetas if key < size), default=None
        )
        if previous is not None:
            older = Band(previous)
            embedded = np.zeros((len(band.modes), 4))
            table = best_thetas[previous].reshape(len(older.modes), 4)
            for m, wave in enumerate(older.modes):
                embedded[band.modes.index(wave)] = table[m]
            seeds.append(embedded.ravel())
        for full in (True, False):
            label = "full" if full else "in_band"
            best = (1.0, None)
            trace: list[tuple[float, float]] = []
            started = time.time()
            for index in range(starts + len(seeds)):
                if index < len(seeds):
                    theta = seeds[index] + 1e-3 * rng.normal(
                        size=band.parameter_count
                    )
                else:
                    theta = rng.normal(size=band.parameter_count)
                theta, ratio = ascend(
                    band, theta, full=full, steps=steps, trace=trace
                )
                if 1.0 - ratio < best[0]:
                    best = (1.0 - ratio, theta)
            entry[f"deficit_min_{label}"] = best[0]
            entry[f"seconds_{label}"] = time.time() - started
            if full and best[1] is not None:
                best_thetas[size] = best[1]
                hat = band.coefficients(best[1])
                record = measure(band, hat, band.nonlinearity(hat))
                entry["minimiser_spread"] = record.spread
                entry["minimiser_leakage_fraction"] = (
                    record.g_full - record.g_band
                ) / record.g_full
                entry["minimiser_energy_neutrality_residual"] = abs(
                    record.t0
                ) / max(abs(record.t1), 1e-300)
            key = f"band{size}_{label}"
            results["cloud"][key] = [
                [float(s), float(d)] for s, d in trace if d > 0.0 and s > 0.0
            ]
            print(
                f"band |k|_inf<={size} G={label}: d_min={best[0]:.6e} "
                f"({entry[f'seconds_{label}']:.1f}s)",
                flush=True,
            )
        results["bands"].append(entry)
    results["_best_thetas"] = best_thetas
    return results


def stage_constrained(
    bands: tuple[int, ...], targets: tuple[float, ...], starts: int, steps: int,
    seed: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(seed)
    for size in bands:
        band = Band(size)
        for target in targets:
            # Relay child amplitude (B=C=1) that hits V_0/mu^2 = target:
            # 24 D^2/(2+9D^2)^2 = target, solved by Newton from the small-D
            # root D^2 ~ target/6.
            square = target / 6.0
            for _ in range(80):
                value = 24.0 * square - target * (2.0 + 9.0 * square) ** 2
                slope = 24.0 - 18.0 * target * (2.0 + 9.0 * square)
                square = max(1e-18, square - value / slope)
            seeds = [
                vector
                for vector in (
                    relay_seed(band, child)
                    for child in (math.sqrt(square), 0.4, 0.1, 0.02)
                )
                if vector is not None
            ]
            best = None

            def consider(vector: np.ndarray, tag: str) -> None:
                nonlocal best
                hat = band.coefficients(vector)
                record = measure(band, hat, band.nonlinearity(hat))
                deficit = record.deficit(full=True)
                if record.spread <= 0.0 or deficit <= 0.0:
                    return
                ratio = deficit / record.spread
                if best is None or ratio < best[0]:
                    best = (ratio, deficit, record.spread, tag)

            for index in range(starts + len(seeds)):
                if index < len(seeds):
                    start = seeds[index]
                    consider(start / np.linalg.norm(start), "seed")
                    theta = start + 1e-3 * rng.normal(size=band.parameter_count)
                else:
                    theta = rng.normal(size=band.parameter_count)
                theta = descend_constrained(
                    band, theta, target=target, steps=steps
                )
                consider(
                    theta, "descent_from_seed" if index < len(seeds) else "descent"
                )
            if best is None:
                continue
            rows.append(
                {
                    "band_inf_norm": size,
                    "target_spread": target,
                    "achieved_spread": best[2],
                    "constraint_log_residual": math.log(best[2] / target),
                    "deficit_full": best[1],
                    "deficit_over_spread": best[0],
                    "origin": best[3],
                }
            )
            print(
                f"  constrained band{size} v={target:g}: "
                f"V/mu^2={best[2]:.4e} d={best[1]:.4e} ratio={best[0]:.4f} "
                f"[{best[3]}]",
                flush=True,
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--starts", type=int, default=10)
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260802)
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    started = time.time()
    telescoping = stage_telescoping()
    closed_form = stage_closed_form()
    search = stage_search(
        (1, 2, 3), arguments.starts, arguments.steps, arguments.seed
    )
    best_thetas = search.pop("_best_thetas")
    constrained = stage_constrained(
        (2, 3),
        (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0),
        max(2, arguments.starts // 3),
        arguments.steps,
        arguments.seed + 1,
    )

    certificates: list[dict[str, Any]] = []
    conditions: list[dict[str, Any]] = []
    for size, theta in sorted(best_thetas.items()):
        band = Band(size)
        widest = min(len(band.modes), 64)
        for keep, bits in ((8, 14), (24, 16), (widest, 16), (widest, 20)):
            field = rationalise(band, theta, keep=keep, bits=bits)
            if field is None:
                continue
            label = f"band{size}_keep{keep}_dyadic2^{bits}"
            try:
                certificate = exact_certificate(field, label=label)
            except (ValueError, AssertionError) as error:
                certificates.append({"label": label, "status": str(error)})
                continue
            certificates.append(certificate)
            print(
                f"  exact certificate {label}: "
                f"d*_full={certificate['full'].get('optimal_deficit_float')}",
                flush=True,
            )
    for name, field in (
        ("relay_triad_D=1/8", build_exact_relay_triad(scale=1)),
        (
            "relay_triad_D=1/128",
            build_exact_relay_triad(scale=1, child_cosine=Fraction(1, 128)),
        ),
    ):
        certificates.append(exact_certificate(field, label=name))
        conditions.append({"field": name, **equality_conditions(field)})

    envelope: dict[str, list[list[float]]] = {}
    for key, points in search["cloud"].items():
        bins: dict[int, float] = {}
        for spread, deficit in points:
            index = int(math.floor(4.0 * math.log10(spread)))
            bins[index] = min(bins.get(index, 1.0), deficit)
        envelope[key] = [
            [10.0 ** (index / 4.0), value] for index, value in sorted(bins.items())
        ]

    summary = {
        "schema": OUTPUT_SCHEMA,
        "status": STATUS,
        "elapsed_seconds": time.time() - started,
        "reduction": {
            "statement": (
                "min over nu>0 of the (I.3) saturation deficit equals "
                "d_*(u) = 1 - Cov^2 H_r/(V_r G_r), attained at "
                "nu_* = Cov/(2 V_r); d_* is viscosity-free, invariant under "
                "u -> alpha u and under u -> lambda u(lambda x), hence no "
                "critical normalisation can change its range."
            )
        },
        "exact_telescoping": telescoping,
        "closed_form_relay_family": closed_form,
        "adversarial_search": search,
        "constrained_search": constrained,
        "lower_envelope": envelope,
        "exact_certificates": certificates,
        "equality_conditions": conditions,
    }
    write_with_digest(
        arguments.output_dir / "summary.json", canonical_json_bytes(summary)
    )
    print(json.dumps({"elapsed_seconds": summary["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
