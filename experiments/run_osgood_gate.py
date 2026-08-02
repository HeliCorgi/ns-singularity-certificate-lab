"""Osgood-gate computations for the spectral front monotone.

Decides the single question of Verification Sprint V1's Osgood gate: can the
front wavenumber ``K = ||P(u.grad u)||_2^2 / ||grad u||_2^4`` of divergence
free trigonometric fields be bounded by ``C(1 + log N_0^2)`` (plus an
integrable remainder), which via the Osgood closure lemma would upgrade the
Lambda dichotomy to unconditional periodic regularity?

Static lane: builds the coherent critical-spectrum family
``u_hat(k) = |k|^{-2} P_k v_0`` (all phases aligned) on bands ``1 <= |k| <= N``
and measures ``K``, ``N_0^2``, and the ratio ``K/(1+log N_0^2)`` as ``N``
grows, against control families (single shell, Beltrami perturbation,
mesoscopic two-box parent).  An exact rational anchor at ``N = 4`` is computed
with the ``fourier_torus`` machinery and compared against the float path.

Dynamic lane: evolves the coherent field by the dealiased Fourier-Galerkin
Navier-Stokes ODE over half a parabolic time ``0.5 N^{-2}`` and tracks
``K(t), N_0^2(t), D(t), z(t)`` plus the (I.4) gap, with a step-halving
control, to test whether the large-``K`` configuration is dynamically
transient or persistent.

Binary64 except where labelled exact.  Diagnostic only; no PDE theorem.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    write_with_digest,
)
from ns_certificate_lab.fourier_torus import TrigVector, advection, leray
from ns_certificate_lab.leray_response_relay import (
    gradient_l2_squared,
    leray_advection,
    mean_energy,
)
from ns_certificate_lab.mesoscopic_cloud_scaling import (
    MesoscopicCloudConfig,
    build_sparse_parent,
)

SEED_VECTOR = np.array([1.0, 2.0, 3.0])


def _integer_frequencies(grid: int) -> np.ndarray:
    return np.rint(np.fft.fftfreq(grid, d=1.0 / grid)).astype(np.int64)


def _frequency_mesh(grid: int):
    f = _integer_frequencies(grid).astype(np.float64)
    return np.meshgrid(f, f, f, indexing="ij")


def _band_waves(band: int) -> list[tuple[int, int, int]]:
    """Positive-representative waves with 1 <= |k|^2 <= band^2."""

    waves = []
    for kx in range(-band, band + 1):
        for ky in range(-band, band + 1):
            for kz in range(-band, band + 1):
                n = kx * kx + ky * ky + kz * kz
                if n == 0 or n > band * band:
                    continue
                if (kx, ky, kz) < (0, 0, 0):
                    continue
                waves.append((kx, ky, kz))
    return waves


def coherent_field(grid: int, band: int) -> np.ndarray:
    """Coherent critical-spectrum field: u_hat(k) = |k|^{-2} P_k v0, phases 1."""

    field = np.zeros((3, grid, grid, grid), dtype=np.complex128)
    for wave in _band_waves(band):
        k = np.asarray(wave, dtype=np.float64)
        nsq = float(k @ k)
        e = SEED_VECTOR - k * (k @ SEED_VECTOR) / nsq
        if np.linalg.norm(e) < 1.0e-9:
            continue
        coefficient = e / nsq
        idx = tuple(int(c) % grid for c in wave)
        nidx = tuple(int(-c) % grid for c in wave)
        for comp in range(3):
            field[(comp, *idx)] += coefficient[comp]
            field[(comp, *nidx)] += coefficient[comp]
    return field


def single_shell_field(grid: int, band: int) -> np.ndarray:
    field = np.zeros((3, grid, grid, grid), dtype=np.complex128)
    target = band * band
    for wave in _band_waves(band):
        if sum(c * c for c in wave) != target:
            continue
        k = np.asarray(wave, dtype=np.float64)
        e = SEED_VECTOR - k * (k @ SEED_VECTOR) / target
        if np.linalg.norm(e) < 1.0e-9:
            continue
        idx = tuple(int(c) % grid for c in wave)
        nidx = tuple(int(-c) % grid for c in wave)
        for comp in range(3):
            field[(comp, *idx)] += e[comp]
            field[(comp, *nidx)] += e[comp]
    return field


def beltrami_perturbed_field(grid: int, band: int, delta: float) -> np.ndarray:
    """ABC-type Beltrami core plus a small coherent perturbation."""

    field = np.zeros((3, grid, grid, grid), dtype=np.complex128)

    def add(wave, cvec, svec):
        idx = tuple(int(c) % grid for c in wave)
        nidx = tuple(int(-c) % grid for c in wave)
        for comp in range(3):
            coeff = 0.5 * (cvec[comp] - 1.0j * svec[comp])
            field[(comp, *idx)] += coeff
            field[(comp, *nidx)] += np.conjugate(coeff)

    add((1, 0, 0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))
    add((0, 1, 0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    add((0, 0, 1), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0))
    return field + delta * coherent_field(grid, band)


def measure(field: np.ndarray) -> dict[str, float]:
    grid = field.shape[1]
    kx, ky, kz = _frequency_mesh(grid)
    wsq = kx * kx + ky * ky + kz * kz
    density = np.sum(np.abs(field) ** 2, axis=0)
    h0 = 0.5 * float(np.sum(density))
    h1 = 0.5 * float(np.sum(wsq * density))
    h2 = 0.5 * float(np.sum(wsq * wsq * density))
    physical = np.fft.ifftn(field, axes=(1, 2, 3)) * float(grid**3)
    sup = float(np.max(np.sqrt(np.sum(np.abs(physical) ** 2, axis=0))))
    nonlinear = leray_advection(field, field)
    nsq = float(np.vdot(nonlinear, nonlinear).real)
    gradsq = gradient_l2_squared(field)
    k_value = nsq / gradsq**2
    n0sq = h1 / h0
    # Leray retention: ||P(u.grad u)|| / ||u.grad u|| (projection loss only).
    raw = _raw_advection(field)
    rawsq = float(np.vdot(raw, raw).real)
    return {
        "H0": h0,
        "H1": h1,
        "H2": h2,
        "sup_u": sup,
        "grad_sq": gradsq,
        "nonlinear_sq": nsq,
        "K": k_value,
        "N0_sq": n0sq,
        "z": math.log(n0sq),
        "K_over_1plusz": k_value / (1.0 + math.log(n0sq)),
        "leray_retention": math.sqrt(nsq / rawsq) if rawsq > 0 else 0.0,
        "lemma_k_bound_sup": sup * sup / gradsq,
    }


def _raw_advection(field: np.ndarray) -> np.ndarray:
    grid = field.shape[1]
    freqs = _frequency_mesh(grid)
    scale = float(grid**3)
    physical = np.fft.ifftn(field, axes=(1, 2, 3)) * scale
    out = np.empty_like(field)
    for comp in range(3):
        acc = np.zeros((grid, grid, grid), dtype=np.complex128)
        for direction in range(3):
            dhat = 1.0j * freqs[direction] * field[comp]
            dphys = np.fft.ifftn(dhat, axes=(0, 1, 2)) * scale
            acc += physical[direction] * dphys
        out[comp] = np.fft.fftn(acc, axes=(0, 1, 2)) / scale
    return out


def exact_anchor(band: int) -> dict[str, str | float]:
    """Exact rational K for the coherent family at a small band."""

    modes = []
    for wave in _band_waves(band):
        k = wave
        nsq = sum(c * c for c in k)
        v0 = (Fraction(1), Fraction(2), Fraction(3))
        kdotv = sum(Fraction(c) * v for c, v in zip(k, v0))
        e = tuple(
            v - Fraction(c) * kdotv / nsq for c, v in zip(k, v0)
        )
        if all(component == 0 for component in e):
            continue
        cvec = tuple(2 * component / nsq for component in e)
        modes.append((k, cvec, (0, 0, 0)))
    field = TrigVector.from_modes(modes)
    if field.divergence().cleaned().terms:
        raise AssertionError("exact coherent field lost divergence-freeness")
    nonlinear = leray(advection(field, field)).cleaned()
    nl_sq = Fraction(0)
    for _, pair in nonlinear.coefficient_table().items():
        cosine, sine = pair
        nl_sq += Fraction(1, 2) * sum(
            Fraction(v) ** 2 for vec in (cosine, sine) for v in vec
        )
    h0 = Fraction(0)
    h1 = Fraction(0)
    for wave, pair in field.coefficient_table().items():
        if wave == (0, 0, 0):
            continue
        cosine, sine = pair
        e_k = Fraction(1, 2) * sum(
            Fraction(v) ** 2 for vec in (cosine, sine) for v in vec
        )
        h0 += e_k
        h1 += sum(c * c for c in wave) * e_k
    # _pair_energy(C, S) = (|C|^2+|S|^2)/2 is exactly the L^2 mean of
    # C cos(k.x) + S sin(k.x), so nl_sq = ||N||_2^2 and h1 = ||grad u||_2^2
    # with no further convention factor; hence K = nl_sq/h1^2.
    k_exact = nl_sq / h1**2
    n0_sq = h1 / h0
    return {
        "band": band,
        "K_exact": str(k_exact),
        "K_exact_float": float(k_exact),
        "N0_sq_exact": str(n0_sq),
        "N0_sq_exact_float": float(n0_sq),
    }


def _rhs(state, *, viscosity, wsq, mask):
    retained = np.asarray(state * mask[None, ...], dtype=np.complex128)
    nonlinear = leray_advection(retained, retained)
    return np.asarray(
        (-nonlinear - viscosity * wsq[None, ...] * retained) * mask[None, ...],
        dtype=np.complex128,
    )


def persistence_run(
    *,
    band: int,
    cutoff: int,
    grid: int,
    viscosity: float,
    amplitude: float,
    steps: int,
    samples: int,
) -> list[dict[str, float]]:
    if 2 * cutoff >= grid // 2:
        raise ValueError("grid too small for dealiased products")
    field = amplitude * coherent_field(grid, band)
    kx, ky, kz = _frequency_mesh(grid)
    wsq = kx * kx + ky * ky + kz * kz
    mask = (
        (np.abs(kx) <= cutoff) & (np.abs(ky) <= cutoff) & (np.abs(kz) <= cutoff)
    )
    total_time = 0.5 / float(band * band)
    dt = total_time / steps
    records = []
    state = field
    sample_every = max(1, steps // samples)
    for step in range(steps + 1):
        if step % sample_every == 0 or step == steps:
            m = measure(state)
            m["time"] = step * dt
            m["time_over_parabolic"] = step * dt * band * band
            records.append(m)
        if step == steps:
            break
        k1 = _rhs(state, viscosity=viscosity, wsq=wsq, mask=mask)
        k2 = _rhs(state + 0.5 * dt * k1, viscosity=viscosity, wsq=wsq, mask=mask)
        k3 = _rhs(state + 0.5 * dt * k2, viscosity=viscosity, wsq=wsq, mask=mask)
        k4 = _rhs(state + dt * k3, viscosity=viscosity, wsq=wsq, mask=mask)
        state = np.asarray(
            (state + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)) * mask[None, ...],
            dtype=np.complex128,
        )
        state[:, 0, 0, 0] = 0.0
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    static_rows = []
    for band, grid in ((4, 48), (8, 80), (16, 144), (24, 208), (32, 272)):
        row = {"family": "coherent_critical", "band": band, "grid": grid}
        row.update(measure(coherent_field(grid, band)))
        static_rows.append(row)
        print(
            f"coherent band={band}: K={row['K']:.4f} N0^2={row['N0_sq']:.3f} "
            f"K/(1+z)={row['K_over_1plusz']:.4f} "
            f"retention={row['leray_retention']:.3f}",
            flush=True,
        )
    for band, grid in ((8, 80), (16, 144)):
        row = {"family": "single_shell", "band": band, "grid": grid}
        row.update(measure(single_shell_field(grid, band)))
        static_rows.append(row)
        print(f"single_shell band={band}: K={row['K']:.6f}", flush=True)
    row = {"family": "beltrami_perturbed", "band": 8, "grid": 80, "delta": 0.05}
    row.update(measure(beltrami_perturbed_field(80, 8, 0.05)))
    static_rows.append(row)
    print(f"beltrami+0.05*coherent: K={row['K']:.6f}", flush=True)
    config = MesoscopicCloudConfig(base_scale=16, gamma=1.0, width_override=3)
    sparse = build_sparse_parent(config)
    grid = 144
    cloud = np.zeros((3, grid, grid, grid), dtype=np.complex128)
    for wave, coefficient in sparse.items():
        idx = tuple(int(c) % grid for c in wave)
        nidx = tuple(int(-c) % grid for c in wave)
        vec = np.asarray(coefficient, dtype=np.complex128)
        for comp in range(3):
            cloud[(comp, *idx)] += vec[comp]
            cloud[(comp, *nidx)] += np.conjugate(vec[comp])
    row = {"family": "mesoscopic_two_box", "band": 16, "grid": grid}
    row.update(measure(cloud))
    static_rows.append(row)
    print(f"mesoscopic two-box N=16: K={row['K']:.6f}", flush=True)

    anchor = exact_anchor(4)
    float_anchor = next(
        r for r in static_rows
        if r["family"] == "coherent_critical" and r["band"] == 4
    )
    anchor["float_K"] = float_anchor["K"]
    anchor["relative_difference"] = abs(
        anchor["K_exact_float"] - float_anchor["K"]
    ) / anchor["K_exact_float"]
    print(
        f"exact anchor band=4: K={anchor['K_exact_float']:.6f} "
        f"(float {float_anchor['K']:.6f}, "
        f"rel diff {anchor['relative_difference']:.2e})",
        flush=True,
    )

    persistence = {}
    for label, amplitude, steps in (
        ("moderate", 4.0, 192),
        ("strong", 16.0, 768),
    ):
        records = persistence_run(
            band=8,
            cutoff=20,
            grid=96,
            viscosity=1.0 / 40.0,
            amplitude=amplitude,
            steps=steps,
            samples=12,
        )
        control = persistence_run(
            band=8,
            cutoff=20,
            grid=96,
            viscosity=1.0 / 40.0,
            amplitude=amplitude,
            steps=2 * steps,
            samples=12,
        )
        drift = abs(records[-1]["K"] - control[-1]["K"]) / max(
            records[-1]["K"], 1.0e-30
        )
        persistence[label] = {
            "amplitude": amplitude,
            "records": records,
            "step_halving_relative_drift_K": drift,
        }
        print(
            f"persistence {label}: K(0)={records[0]['K']:.4f} -> "
            f"K(0.5/N^2)={records[-1]['K']:.4f} "
            f"min K={min(r['K'] for r in records):.4f} "
            f"dt-drift={drift:.2e}",
            flush=True,
        )

    summary = {
        "schema": "ns-certificate-lab/osgood-gate/v1",
        "status": "BINARY64 + ONE EXACT ANCHOR / NOT A PROOF",
        "static_rows": static_rows,
        "exact_anchor": anchor,
        "persistence": persistence,
    }
    write_with_digest(
        args.output_dir / "summary.json", canonical_json_bytes(summary)
    )
    print("written", args.output_dir / "summary.json")


if __name__ == "__main__":
    main()
