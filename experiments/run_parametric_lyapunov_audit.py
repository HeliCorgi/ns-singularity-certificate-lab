"""VERIFICATION SPRINT v1 / WORKSTREAM D.

Audit of Proposition 1 ("tree no-go") and the cyclic holonomy invariant C of
``docs/research_notes/ideas_2026_08_01/idea_parametric_resonance.md``.

The audited object is the CHILD LINEAR SYSTEM

    xdot = A(t) x,   A(t) = -nu * D + K(t),   D = diag(|k_i|^2),

with K(t) zero-diagonal, supported on a coupling graph G, and each edge
e = {i,j} carrying constants (beta_fwd, beta_bwd) modulated by scalar
functions of time.  Proposition 1 claims: G a forest  =>  the top Floquet
exponent is <= -nu * min_i |k_i|^2 for ALL admissible modulations.

This script does not generate ideas.  It only (i) validates the corrected
proposition numerically over an exhaustive/randomized grid, (ii) exhibits
counterexamples for the hypotheses the note left implicit, and (iii) re-derives
and stress-tests the 3-cycle holonomy condition, including the note's own
C = -1239/128 gadget in exact rational arithmetic.

ARITHMETIC LABELS (mandatory repo convention):
  * EXACT (fractions.Fraction): all beta-coefficients of the B.5 gadget, the
    edge ratios r_e, the cycle invariant C, det K, c1 = sum omega_e^2, the
    balanced-cycle construction, and the diagonal metric Theta.
  * FLOAT (numpy float64): every matrix exponential, monodromy, spectral
    radius, Floquet exponent, and the delta-scans.  No float quantity is
    load-bearing for a "PROVEN" label anywhere in the companion note.

scipy is NOT available in this repo; expm is a hand-rolled batched
scaling-and-squaring Taylor evaluation (see ``expm_batch``).

Run:
    .venv/Scripts/python.exe experiments/run_parametric_lyapunov_audit.py
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import numpy as np

SCHEMA = "ns-certificate-lab/parametric-common-lyapunov-audit/v1"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPOSITORY_ROOT / "outputs" / "verification_sprint_v1" / "parametric"

# Absolute tolerance below which a Floquet margin is called "zero".  The
# float64 round-off floor observed on the exactly-balanced null controls is
# ~1e-9 in absolute margin (see outputs JSON), so 1e-7 is a safe, honest
# threshold; every table also reports the raw max margin.
MARGIN_TOL = 1.0e-7


# --------------------------------------------------------------------------
# batched matrix exponential (numpy only, no scipy)
# --------------------------------------------------------------------------


def _inf_norm(mats: np.ndarray) -> np.ndarray:
    """max_i sum_j |A_ij| for a stack of matrices, floored away from zero."""
    return np.maximum(np.abs(mats).sum(axis=-1).max(axis=-1), 1e-300)


def expm_batch(
    mats: np.ndarray, taylor_terms: int = 22
) -> tuple[np.ndarray, np.ndarray]:
    """exp(M) = exp(logscale) * E with ||E||_inf = 1.  FLOAT.

    Scaling-and-squaring + Taylor, with the partial result renormalized to
    unit infinity-norm after the Taylor step and after every squaring, the log
    of the discarded factor being carried additively.  This makes the routine
    immune to both overflow (growing generators, where ||exp(A)|| can exceed
    1e308) and underflow (strongly damped generators).  A single global
    squaring exponent is used for the whole batch so the routine stays
    branch-free; with ||A/2^s||_inf <= 1/4 the truncated Taylor series at 22
    terms is below float64 round-off.
    """
    mats = np.asarray(mats, dtype=np.float64)
    if mats.size == 0:
        return mats.copy(), np.zeros(mats.shape[:-2])
    dim = mats.shape[-1]
    eye = np.eye(dim)
    nrm = float(np.max(np.abs(mats).sum(axis=-1)))
    squarings = 0
    if nrm > 0.25:
        squarings = int(math.ceil(math.log2(nrm / 0.25)))
    scaled = mats / (2.0**squarings)
    out = np.broadcast_to(eye, mats.shape).copy()
    term = out.copy()
    for order in range(1, taylor_terms + 1):
        term = (term @ scaled) / order
        out = out + term
    scale = _inf_norm(out)
    out = out / scale[..., None, None]
    logscale = np.log(scale)
    for _ in range(squarings):
        out = out @ out
        step = _inf_norm(out)
        out = out / step[..., None, None]
        logscale = 2.0 * logscale + np.log(step)
    return out, logscale


def spectral_radius_batch(mats: np.ndarray) -> np.ndarray:
    """rho(M) for a stack of matrices.  FLOAT."""
    eigs = np.linalg.eigvals(mats)
    return np.max(np.abs(eigs), axis=-1)


# --------------------------------------------------------------------------
# monodromy of a piecewise-constant modulated child system
# --------------------------------------------------------------------------


def build_generators(
    edges: Sequence[tuple[int, int]],
    beta_fwd: np.ndarray,
    beta_bwd: np.ndarray,
    mod_fwd: np.ndarray,
    mod_bwd: np.ndarray,
    ksq: np.ndarray,
    nu: float,
) -> np.ndarray:
    """A_j = -nu D + K_j for every batch element and every switching interval.

    Shapes: beta_* (B, E); mod_* (B, S, E); ksq (n,).  Returns (B, S, n, n).
    Edge (i, j) means x_j is driven by x_i with beta_fwd, x_i by x_j with
    beta_bwd, i.e. K[j, i] = beta_fwd * m_fwd and K[i, j] = beta_bwd * m_bwd.
    """
    batch, steps, n_edges = mod_fwd.shape
    dim = len(ksq)
    gen = np.zeros((batch, steps, dim, dim), dtype=np.float64)
    for idx, (i, j) in enumerate(edges):
        gen[:, :, j, i] = beta_fwd[:, None, idx] * mod_fwd[:, :, idx]
        gen[:, :, i, j] = beta_bwd[:, None, idx] * mod_bwd[:, :, idx]
    diag = -nu * np.asarray(ksq, dtype=np.float64)
    gen[:, :, np.arange(dim), np.arange(dim)] += diag[None, None, :]
    return gen


def monodromy(
    generators: np.ndarray, durations: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """M = exp(A_S t_S) ... exp(A_1 t_1) as (M_normalized, log_scale).

    The true monodromy is ``exp(log_scale)[:, None, None] * M_normalized`` with
    ``||M_normalized||_inf = 1``; the split is what makes strongly growing or
    strongly decaying products representable in float64.  FLOAT.
    """
    batch, steps, dim, _ = generators.shape
    scaled = generators * durations[:, :, None, None]
    flat, shift = expm_batch(scaled.reshape(batch * steps, dim, dim))
    blocks = flat.reshape(batch, steps, dim, dim)
    shifts = shift.reshape(batch, steps)
    out = blocks[:, 0]
    logs = shifts[:, 0].copy()
    for step in range(1, steps):
        out = blocks[:, step] @ out
        logs += shifts[:, step]
        nrm = _inf_norm(out)
        out = out / nrm[:, None, None]
        logs += np.log(nrm)
    return out, logs


def log_spectral_radius(generators: np.ndarray, durations: np.ndarray) -> np.ndarray:
    """log rho(M) for the piecewise-constant monodromy.  FLOAT."""
    mono, logs = monodromy(generators, durations)
    rho = spectral_radius_batch(mono)
    return np.log(np.maximum(rho, 1e-300)) + logs


def floquet_margin(
    edges: Sequence[tuple[int, int]],
    beta_fwd: np.ndarray,
    beta_bwd: np.ndarray,
    mod_fwd: np.ndarray,
    mod_bwd: np.ndarray,
    durations: np.ndarray,
    ksq: np.ndarray,
    nu: float,
) -> tuple[np.ndarray, np.ndarray]:
    """(margin, chi) where chi = log rho(M)/T and margin = chi + nu*min|k|^2.

    margin > 0 is a counterexample to the Proposition-1 claim.  FLOAT.
    """
    gen = build_generators(edges, beta_fwd, beta_bwd, mod_fwd, mod_bwd, ksq, nu)
    period = durations.sum(axis=1)
    chi = log_spectral_radius(gen, durations) / period
    return chi + nu * float(np.min(ksq)), chi


# --------------------------------------------------------------------------
# random family samplers
# --------------------------------------------------------------------------

GRAPHS: dict[str, tuple[int, tuple[tuple[int, int], ...]]] = {
    "edge_2node": (2, ((0, 1),)),
    "path_3node": (3, ((0, 1), (1, 2))),
    "star_4node": (4, ((0, 1), (0, 2), (0, 3))),
    "path_4node": (4, ((0, 1), (1, 2), (2, 3))),
    "cycle_3node": (3, ((0, 1), (1, 2), (2, 0))),
}


def _loguniform(rng: np.random.Generator, shape, lo: float, hi: float) -> np.ndarray:
    return np.exp(rng.uniform(math.log(lo), math.log(hi), size=shape))


def sample_betas(
    rng: np.random.Generator,
    batch: int,
    n_edges: int,
    *,
    elliptic: bool,
    mag_lo: float,
    mag_hi: float,
) -> tuple[np.ndarray, np.ndarray]:
    """beta_fwd, beta_bwd.  elliptic=True enforces beta_fwd*beta_bwd < 0 (H2)."""
    fwd = _loguniform(rng, (batch, n_edges), mag_lo, mag_hi)
    fwd *= rng.choice(np.array([-1.0, 1.0]), size=(batch, n_edges))
    ratio = _loguniform(rng, (batch, n_edges), 0.05, 20.0)
    if elliptic:
        bwd = -np.sign(fwd) * np.abs(fwd) * ratio
    else:
        bwd = np.sign(fwd) * np.abs(fwd) * ratio
    return fwd, bwd


def sample_modulation(
    rng: np.random.Generator,
    batch: int,
    steps: int,
    n_edges: int,
    modulation_class: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (mod_fwd, mod_bwd), shape (B, S, E).

    class_a  : K(t) = m(t) K0, m piecewise constant, one scalar for all edges.
    class_b  : common smooth (sinusoidal) scalar sampled onto the S intervals.
    class_c  : fully independent per-edge m_e(t).
    class_d  : (H1)-BREAKING -- independent modulation of the two DIRECTIONS of
               each edge.  Not realizable by a single parent amplitude; used
               only to show (H1) is load-bearing.
    """
    if modulation_class == "class_a":
        scal = rng.uniform(-3.0, 3.0, size=(batch, steps, 1))
        mod = np.broadcast_to(scal, (batch, steps, n_edges)).copy()
        return mod, mod
    if modulation_class == "class_b":
        amp = rng.uniform(0.2, 3.0, size=(batch, 1, 1))
        off = rng.uniform(-2.0, 2.0, size=(batch, 1, 1))
        phase = rng.uniform(0.0, 2.0 * math.pi, size=(batch, 1, 1))
        grid = (np.arange(steps, dtype=np.float64) + 0.5) / steps
        wave = off + amp * np.sin(2.0 * math.pi * grid[None, :, None] + phase)
        return (
            np.broadcast_to(wave, (batch, steps, n_edges)).copy(),
            np.broadcast_to(wave, (batch, steps, n_edges)).copy(),
        )
    if modulation_class == "class_c":
        mod = rng.uniform(-3.0, 3.0, size=(batch, steps, n_edges))
        # inject exact zeros (staggered / duty-cycled pumps) about 25% of the time
        mask = rng.random((batch, steps, n_edges)) < 0.25
        mod = np.where(mask, 0.0, mod)
        return mod, mod
    if modulation_class == "class_d":
        fwd = rng.uniform(-3.0, 3.0, size=(batch, steps, n_edges))
        bwd = rng.uniform(-3.0, 3.0, size=(batch, steps, n_edges))
        return fwd, bwd
    raise ValueError(f"unknown modulation class {modulation_class!r}")


def sample_ksq(rng: np.random.Generator, batch: int, dim: int, mode: str) -> np.ndarray:
    if mode == "single_shell":
        return np.full(dim, 9.0)
    return np.array([2.0, 6.0, 9.0, 14.0][:dim], dtype=np.float64)


def run_family(
    rng: np.random.Generator,
    graph: str,
    modulation_class: str,
    *,
    elliptic: bool,
    balanced_cycle: bool,
    steps: int,
    batch: int,
    nu: float,
    shell_mode: str,
    mag_lo: float,
    mag_hi: float,
    dur_lo: float,
    dur_hi: float,
) -> dict[str, Any]:
    dim, edges = GRAPHS[graph]
    n_edges = len(edges)
    fwd, bwd = sample_betas(
        rng, batch, n_edges, elliptic=elliptic, mag_lo=mag_lo, mag_hi=mag_hi
    )
    if balanced_cycle:
        # force C = r_a r_b r_c = -1 exactly (up to float64): choose r_a, r_b
        # freely negative and set r_c = -1/(r_a r_b).
        ratio = bwd / fwd  # all negative when elliptic
        ratio[:, 2] = -1.0 / (ratio[:, 0] * ratio[:, 1])
        bwd = ratio * fwd
    mod_f, mod_b = sample_modulation(rng, batch, steps, n_edges, modulation_class)
    durations = _loguniform(rng, (batch, steps), dur_lo, dur_hi)
    ksq = sample_ksq(rng, batch, dim, shell_mode)
    margin, chi = floquet_margin(
        edges, fwd, bwd, mod_f, mod_b, durations, ksq, nu
    )
    finite = np.isfinite(margin)
    margin = margin[finite]
    chi = chi[finite]
    worst = int(np.argmax(margin))
    record: dict[str, Any] = {
        "graph": graph,
        "topology": "forest" if graph != "cycle_3node" else "cycle",
        "modulation_class": modulation_class,
        "hypotheses": {
            "H1_reciprocal_shared_parent": modulation_class != "class_d",
            "H2_elliptic_edges": bool(elliptic),
            "H3_balanced_cycles": graph != "cycle_3node" or bool(balanced_cycle),
        },
        "switching_intervals": steps,
        "samples": int(margin.size),
        "nu": nu,
        "shell_mode": shell_mode,
        "ksq": ksq.tolist(),
        "coupling_magnitude_box": [mag_lo, mag_hi],
        "edge_ratio_box": [0.05, 20.0],
        "interval_length_box": [dur_lo, dur_hi],
        "max_margin": float(np.max(margin)),
        "max_relative_margin": float(np.max(margin) / (nu * float(np.min(ksq)))),
        "median_margin": float(np.median(margin)),
        "min_margin": float(np.min(margin)),
        "frac_positive_margin": float(np.mean(margin > MARGIN_TOL)),
        "max_chi": float(np.max(chi)),
        "worst_case": {
            "beta_fwd": fwd[finite][worst].tolist(),
            "beta_bwd": bwd[finite][worst].tolist(),
            "durations": durations[finite][worst].tolist(),
            "mod_fwd": mod_f[finite][worst].tolist(),
            "mod_bwd": mod_b[finite][worst].tolist(),
        },
        "verdict": (
            "NO COUNTEREXAMPLE (margin <= tol)"
            if float(np.max(margin)) <= MARGIN_TOL
            else "COUNTEREXAMPLE FOUND (margin > tol)"
        ),
    }
    return record


# --------------------------------------------------------------------------
# exhaustive small-case grid (2-node and 3-path, piecewise constant)
# --------------------------------------------------------------------------


def run_exhaustive_grid(
    graph: str,
    modulation_class: str,
    *,
    steps: int,
    nu: float,
    elliptic: bool,
    beta_mag: Sequence[float],
    mod_vals: Sequence[float],
    dur_vals: Sequence[float],
    chunk: int = 120_000,
) -> dict[str, Any]:
    """Deterministic full-tensor grid over (beta, modulation, durations).

    Every point of the cartesian product
        beta_mag^(2E) x mod_vals^(S E or S) x dur_vals^S
    is evaluated; no RNG anywhere.  FLOAT.
    """
    import itertools

    dim, edges = GRAPHS[graph]
    n_edges = len(edges)
    sign = 1.0 if elliptic else -1.0

    beta_combos = list(
        itertools.product(itertools.product(beta_mag, beta_mag), repeat=n_edges)
    )
    fwd_tab = np.array([[bb[0] for bb in combo] for combo in beta_combos])
    bwd_tab = np.array([[-sign * bb[1] for bb in combo] for combo in beta_combos])

    if modulation_class == "class_a":
        mod_tab = np.array(
            [
                [[v] * n_edges for v in seq]
                for seq in itertools.product(mod_vals, repeat=steps)
            ],
            dtype=np.float64,
        )
    else:
        mod_tab = np.array(
            list(
                itertools.product(
                    list(itertools.product(mod_vals, repeat=n_edges)), repeat=steps
                )
            ),
            dtype=np.float64,
        )
    dur_tab = np.array(list(itertools.product(dur_vals, repeat=steps)), dtype=np.float64)

    n_beta, n_mod, n_dur = len(fwd_tab), len(mod_tab), len(dur_tab)
    total = n_beta * n_mod * n_dur
    ksq = np.full(dim, 9.0)

    best = -np.inf
    best_index = -1
    for start in range(0, total, chunk):
        idx = np.arange(start, min(start + chunk, total))
        i_beta, rem = np.divmod(idx, n_mod * n_dur)
        i_mod, i_dur = np.divmod(rem, n_dur)
        margin, _ = floquet_margin(
            edges,
            fwd_tab[i_beta],
            bwd_tab[i_beta],
            mod_tab[i_mod],
            mod_tab[i_mod],
            dur_tab[i_dur],
            ksq,
            nu,
        )
        local = int(np.argmax(margin))
        if float(margin[local]) > best:
            best = float(margin[local])
            best_index = int(idx[local])
    i_beta, rem = divmod(best_index, n_mod * n_dur)
    i_mod, i_dur = divmod(rem, n_dur)
    return {
        "kind": "exhaustive_grid",
        "graph": graph,
        "topology": "forest" if graph != "cycle_3node" else "cycle",
        "modulation_class": modulation_class,
        "switching_intervals": steps,
        "nu": nu,
        "elliptic_edges": elliptic,
        "grid": {
            "beta_magnitudes": list(beta_mag),
            "modulation_values": list(mod_vals),
            "interval_lengths": list(dur_vals),
            "ksq": ksq.tolist(),
            "beta_points": n_beta,
            "modulation_points": n_mod,
            "duration_points": n_dur,
        },
        "evaluations": int(total),
        "max_margin": best,
        "max_relative_margin": best / (nu * float(ksq.min())),
        "argmax": {
            "beta_fwd": fwd_tab[i_beta].tolist(),
            "beta_bwd": bwd_tab[i_beta].tolist(),
            "mod": mod_tab[i_mod].tolist(),
            "durations": dur_tab[i_dur].tolist(),
        },
        "verdict": (
            "NO COUNTEREXAMPLE" if best <= MARGIN_TOL else "COUNTEREXAMPLE FOUND"
        ),
    }


# --------------------------------------------------------------------------
# EXACT rational analysis of the note's B.5 / B.6 cyclic gadget
# --------------------------------------------------------------------------

# Raw exact Leray coefficients quoted in idea_parametric_resonance.md B.5:
#   forward:  -8/5, -1/24, 32/45      backward: 2/5, 7/5, -59/72
GADGET_FWD = (Fraction(-8, 5), Fraction(-1, 24), Fraction(32, 45))
GADGET_BWD = (Fraction(2, 5), Fraction(7, 5), Fraction(-59, 72))
NOTE_R = (Fraction(-9, 4), Fraction(-7, 3), Fraction(-59, 32))
NOTE_C = Fraction(-1239, 128)


def cycle_invariants_exact(
    fwd: Sequence[Fraction], bwd: Sequence[Fraction]
) -> dict[str, Any]:
    """EXACT rational invariants of a 3-cycle generator K.

    Edges a={0,1}, b={1,2}, c={2,0} with forward K[1,0], K[2,1], K[0,2].
    C = (product of backward)/(product of forward); det K = F + B = F(1+C);
    c1 = sum of squared elliptic frequencies = -sum_e beta_fwd*beta_bwd.
    """
    ratios = tuple(b / f for f, b in zip(fwd, bwd))
    invariant = ratios[0] * ratios[1] * ratios[2]
    prod_f = fwd[0] * fwd[1] * fwd[2]
    prod_b = bwd[0] * bwd[1] * bwd[2]
    det_k = prod_f + prod_b
    omega_sq = tuple(-f * b for f, b in zip(fwd, bwd))
    c1 = omega_sq[0] + omega_sq[1] + omega_sq[2]
    return {
        "r_edges_exact": [str(r) for r in ratios],
        "r_edges_float": [float(r) for r in ratios],
        "C_exact": str(invariant),
        "C_float": float(invariant),
        "one_plus_C_exact": str(1 + invariant),
        "forward_product_exact": str(prod_f),
        "backward_product_exact": str(prod_b),
        "det_K_exact": str(det_k),
        "det_K_float": float(det_k),
        "omega_sq_exact": [str(w) for w in omega_sq],
        "omega_float": [math.sqrt(float(w)) if w > 0 else float("nan") for w in omega_sq],
        "c1_exact": str(c1),
        "c1_float": float(c1),
        "all_edges_elliptic": all(w > 0 for w in omega_sq),
    }


def steady_cycle_growth(c1: float, det_k: float) -> dict[str, float]:
    """max Re lambda for a traceless elliptic 3-cycle.

    char poly: lam^3 + c1 lam - det K, c1 > 0 => one real root mu_r and one
    complex pair (-mu_r/2 +- i sigma).  FLOAT (root of a cubic).
    """
    roots = np.roots([1.0, 0.0, c1, -det_k])
    real_root = float(roots[np.argmin(np.abs(roots.imag))].real)
    return {
        "mu_real": real_root,
        "sigma_cycle": float(np.max(roots.real)),
        "spectrum_real": [float(z.real) for z in roots],
        "spectrum_imag": [float(z.imag) for z in roots],
        "det_K_sign": float(np.sign(det_k)),
    }


def staggered_chi(
    fwd: Sequence[float], bwd: Sequence[float], deltas: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """chi_par(delta) = log rho(exp(Gc d) exp(Gb d) exp(Ga d)) / (3 delta).

    FLOAT.  Edge order a, b, c; nu = 0 (the viscous shift is exactly
    -nu|k|^2 on a single shell, so it is added analytically later).
    """
    edges = GRAPHS["cycle_3node"][1]
    batch = deltas.size
    fwd_arr = np.tile(np.asarray(fwd, dtype=np.float64), (batch, 1))
    bwd_arr = np.tile(np.asarray(bwd, dtype=np.float64), (batch, 1))
    # one interval per edge; edge e active alone on interval e
    mod = np.zeros((batch, 3, 3))
    for e in range(3):
        mod[:, e, e] = 1.0
    durations = np.repeat(deltas[:, None], 3, axis=1)
    gen = build_generators(
        edges, fwd_arr, bwd_arr, mod, mod, np.zeros(3), 0.0
    )
    log_rho = log_spectral_radius(gen, durations)
    return log_rho / (3.0 * deltas), np.exp(np.minimum(log_rho, 700.0))


def holonomy_stress(rng: np.random.Generator, n_samples: int) -> dict[str, Any]:
    """Scan sigma_cycle and max_delta chi_par as a function of C at fixed omega.

    Shows (i) growth vanishes only at C = -1, (ii) there is no magnitude
    threshold for EXISTENCE of growth, (iii) the sign of det K = F(1+C)
    changes the rate by a factor ~2 at fixed |det K|.  FLOAT.
    """
    omegas = (0.8, 0.24152, 0.76336)  # corrected B.5 elliptic frequencies
    c_values = np.concatenate(
        [
            -np.geomspace(50.0, 1.0000001, 60),
            -np.geomspace(0.9999999, 0.002, 60),
        ]
    )
    deltas = np.geomspace(1e-3, 20.0, 260)
    rows = []
    for c_val in c_values:
        # canonical gauge: r_a = r_b = -1, r_c = C
        f = [omegas[0], omegas[1], omegas[2] / math.sqrt(-c_val)]
        b = [-omegas[0], -omegas[1], -omegas[2] * math.sqrt(-c_val)]
        inv = cycle_invariants_exact(
            [Fraction(x).limit_denominator(10**12) for x in f],
            [Fraction(x).limit_denominator(10**12) for x in b],
        )
        c1 = float(inv["c1_float"])
        det_k = float(inv["det_K_float"])
        steady = steady_cycle_growth(c1, det_k)
        chi, _ = staggered_chi(f, b, deltas)
        best = int(np.argmax(chi))
        rows.append(
            {
                "C": float(c_val),
                "one_plus_C": float(1.0 + c_val),
                "det_K": det_k,
                "c1": c1,
                "sigma_cycle": steady["sigma_cycle"],
                "mu_real": steady["mu_real"],
                "sigma_over_absdetK_c1": steady["sigma_cycle"] * c1 / abs(det_k)
                if det_k != 0
                else float("nan"),
                "chi_par_max": float(chi[best]),
                "delta_star": float(deltas[best]),
                "chi_over_sigma": float(chi[best] / steady["sigma_cycle"])
                if steady["sigma_cycle"] > 0
                else float("nan"),
                "chi_small_delta": float(chi[0]),
                "chi_min_over_scan": float(np.min(chi)),
            }
        )
    # random 3-cycles: is chi_par ever exactly zero at C != -1?
    zero_hits = 0
    n_rand = n_samples
    r_ab = -_loguniform(rng, (n_rand, 2), 0.05, 20.0)
    om = _loguniform(rng, (n_rand, 3), 0.05, 5.0)
    c_rand = r_ab[:, 0] * r_ab[:, 1] * (-_loguniform(rng, n_rand, 0.05, 20.0))
    min_chis = []
    deltas_rand = np.geomspace(1e-3, 20.0, 120)
    for idx in range(min(n_rand, 400)):
        c_val = float(c_rand[idx])
        f = [om[idx, 0], om[idx, 1], om[idx, 2] / math.sqrt(-c_val)]
        b = [-om[idx, 0], -om[idx, 1], -om[idx, 2] * math.sqrt(-c_val)]
        chi, rho = staggered_chi(f, b, deltas_rand)
        min_chis.append(float(np.min(chi)))
        if float(np.min(rho)) < 1.0 - 1e-9:
            zero_hits += 1
    return {
        "omega_used": list(omegas),
        "C_scan": rows,
        "random_cycles_tested": int(min(n_rand, 400)),
        "random_min_chi_over_delta_scan": {
            "min": float(np.min(min_chis)),
            "median": float(np.median(min_chis)),
            "count_below_1e-9": int(np.sum(np.asarray(min_chis) < -1e-9)),
        },
        "count_rho_below_one": zero_hits,
        "note": (
            "det M = exp(delta * tr(Ga+Gb+Gc)) = 1 for traceless generators, "
            "so rho(M) >= 1 always; chi_par >= 0 is structural, not evidence."
        ),
    }


def gadget_report() -> dict[str, Any]:
    """Reproduce and correct the note's B.5/B.6 numbers.  EXACT + FLOAT."""
    exact = cycle_invariants_exact(GADGET_FWD, GADGET_BWD)
    note_c = NOTE_R[0] * NOTE_R[1] * NOTE_R[2]
    # gauge check: r_e from the note vs raw ratios must differ by (theta_i/theta_j)^2
    raw_r = [Fraction(b, 1) / f for f, b in zip(GADGET_FWD, GADGET_BWD)]
    gauge_factors = [nr / rr for nr, rr in zip(NOTE_R, raw_r)]
    gauge_product = gauge_factors[0] * gauge_factors[1] * gauge_factors[2]

    c1 = float(exact["c1_float"])
    det_k = float(exact["det_K_float"])
    steady = steady_cycle_growth(c1, det_k)

    # the note's own (inconsistent) numbers, for side-by-side
    note_omega_c = 1.0796
    note_c1 = 0.8**2 + 0.24152**2 + note_omega_c**2
    note_det = -0.58192
    steady_note = steady_cycle_growth(note_c1, note_det)

    deltas = np.geomspace(1e-3, 30.0, 900)
    chi, rho = staggered_chi(
        [float(x) for x in GADGET_FWD], [float(x) for x in GADGET_BWD], deltas
    )
    best = int(np.argmax(chi))
    # commutator-free floor: chi_par(delta -> 0) = sigma_cycle / 3
    return {
        "source": "idea_parametric_resonance.md sections B.5 / B.6",
        "exact_invariants": exact,
        "note_quoted_C_exact": str(NOTE_C),
        "note_r_product_exact": str(note_c),
        "C_matches_note": exact["C_exact"] == str(NOTE_C) == str(note_c),
        "gauge_check": {
            "note_r_over_raw_r": [str(g) for g in gauge_factors],
            "product_must_be_one": str(gauge_product),
            "diagonal_gauge_consistent": gauge_product == 1,
            "implied_theta_ratios": [float(g) ** 0.5 for g in gauge_factors],
        },
        "corrected_from_raw_coefficients": {
            "omega_a": math.sqrt(float(Fraction(-1) * GADGET_FWD[0] * GADGET_BWD[0])),
            "omega_b": math.sqrt(float(Fraction(-1) * GADGET_FWD[1] * GADGET_BWD[1])),
            "omega_c": math.sqrt(float(Fraction(-1) * GADGET_FWD[2] * GADGET_BWD[2])),
            "c1": c1,
            "det_K": det_k,
            "sigma_cycle": steady["sigma_cycle"],
            "spectrum_real": steady["spectrum_real"],
            "spectrum_imag": steady["spectrum_imag"],
        },
        "note_quoted_values": {
            "omega": [0.8, 0.24152, note_omega_c],
            "det_K": note_det,
            "sigma_cycle": 0.149013,
            "chi_par": 0.052773,
            "delta_star": 0.6593,
            "rho_M": 1.1100,
        },
        "note_reproduced_from_its_own_numbers": {
            "c1": note_c1,
            "sigma_cycle": steady_note["sigma_cycle"],
            "spectrum_real": steady_note["spectrum_real"],
            "spectrum_imag": steady_note["spectrum_imag"],
        },
        "discrepancy": {
            "omega_c_note_over_corrected": note_omega_c
            / math.sqrt(float(Fraction(-1) * GADGET_FWD[2] * GADGET_BWD[2])),
            "det_K_note_over_corrected": note_det / det_k,
            "comment": (
                "Both ratios equal sqrt(2) to 4 digits: a single spurious "
                "sqrt(2) on edge c.  C, r_e and the C=-1 null control are "
                "unaffected (C is gauge- and edge-scaling-invariant); "
                "omega_c, det K and sigma_cycle are not."
            ),
        },
        "staggered_scan": {
            "delta_grid": [float(deltas[0]), float(deltas[-1]), int(deltas.size)],
            "chi_par_max": float(chi[best]),
            "delta_star": float(deltas[best]),
            "rho_at_delta_star": float(rho[best]),
            "chi_par_small_delta_limit": float(chi[0]),
            "sigma_cycle_over_3": steady["sigma_cycle"] / 3.0,
            "chi_over_sigma_cycle": float(chi[best]) / steady["sigma_cycle"],
            "holonomy_excess_over_trotter_floor": float(chi[best])
            / (steady["sigma_cycle"] / 3.0)
            - 1.0,
            "chi_min_over_scan": float(np.min(chi)),
            "rho_min_over_scan": float(np.min(rho)),
        },
    }


def balanced_null_control() -> dict[str, Any]:
    """EXACT construction of a C = -1 3-cycle + FLOAT null control.

    Take r_a = -2, r_b = -3, r_c = -1/6 so C = r_a r_b r_c = -1 exactly.
    Then Theta = diag(theta) with theta_i/theta_j = -beta_ji/beta_ij closes
    around the cycle and Theta K(t) + K(t)^T Theta = 0 for every modulation.
    """
    fwd = (Fraction(3, 2), Fraction(5, 7), Fraction(11, 4))
    ratios = (Fraction(-2), Fraction(-3), Fraction(-1, 6))
    bwd = tuple(r * f for r, f in zip(ratios, fwd))
    exact = cycle_invariants_exact(fwd, bwd)

    # Theta by propagation.  Edges a={0,1}, b={1,2}, c={2,0}; forward entries
    # are K[1,0], K[2,1], K[0,2], backward entries K[0,1], K[1,2], K[2,0].
    # The metric condition Theta K + K^T Theta = 0 reads, entrywise,
    #     theta_i K_ij + theta_j K_ji = 0,  i.e. theta_i/theta_j = -K_ji/K_ij.
    # (The note's proof sketch writes theta_i/theta_j = -beta_ij/beta_ji, the
    #  reciprocal; both are positive under (H2) so existence is unaffected,
    #  but the metric itself is the transpose of the note's.)
    theta = [Fraction(1), Fraction(0), Fraction(0)]
    theta[1] = -theta[0] * bwd[0] / fwd[0]  # theta_0 K[0,1] + theta_1 K[1,0] = 0
    theta[2] = -theta[1] * bwd[1] / fwd[1]  # theta_1 K[1,2] + theta_2 K[2,1] = 0
    closure = theta[2] * bwd[2] + theta[0] * fwd[2]  # must be 0 exactly

    rng = np.random.default_rng(20260802)
    edges = GRAPHS["cycle_3node"][1]
    batch, steps = 40000, 4
    fwd_arr = np.tile(np.array([float(x) for x in fwd]), (batch, 1))
    bwd_arr = np.tile(np.array([float(x) for x in bwd]), (batch, 1))
    mod, _ = sample_modulation(rng, batch, steps, 3, "class_c")
    durations = _loguniform(rng, (batch, steps), 0.05, 4.0)
    ksq = np.full(3, 9.0)
    nu = 0.07
    margin, chi = floquet_margin(
        edges, fwd_arr, bwd_arr, mod, mod, durations, ksq, nu
    )
    return {
        "construction": "exact rational C = -1 cycle (r = -2, -3, -1/6)",
        "beta_fwd_exact": [str(x) for x in fwd],
        "beta_bwd_exact": [str(x) for x in bwd],
        "exact_invariants": exact,
        "theta_exact": [str(t) for t in theta],
        "cycle_closure_residual_exact": str(closure),
        "theta_closes_exactly": closure == 0,
        "condition_number_theta": float(
            max(float(t) for t in theta) / min(float(t) for t in theta)
        ),
        "null_control_float": {
            "samples": batch,
            "switching_intervals": steps,
            "modulation_class": "class_c (fully independent per-edge, incl. zeros)",
            "nu": nu,
            "ksq": ksq.tolist(),
            "max_abs_margin": float(np.max(np.abs(margin))),
            "max_margin": float(np.max(margin)),
            "max_abs_chi_plus_nu_ksq": float(np.max(np.abs(chi + nu * 9.0))),
            "passes_1e-12": bool(np.max(np.abs(margin)) < 1e-12),
            "passes_1e-10": bool(np.max(np.abs(margin)) < 1e-10),
        },
    }


def theta_certificate_check(rng: np.random.Generator, batch: int = 20000) -> dict[str, Any]:
    """Directly verify the common-Lyapunov certificate on random forests.

    For a 4-node path with (H1)+(H2) build Theta by edge-wise propagation and
    check (i) Theta K(t) + K(t)^T Theta = 0 to round-off at random modulation
    values, (ii) the resulting decay bound d/dt ||x||^2_Theta <= -2 nu min|k|^2
    ||x||^2_Theta, and (iii) the Euclidean overshoot factor sqrt(kappa(Theta)).
    FLOAT (the exact-rational version is in ``balanced_null_control``).
    """
    dim, edges = GRAPHS["path_4node"]
    fwd, bwd = sample_betas(rng, batch, len(edges), elliptic=True, mag_lo=0.05, mag_hi=20.0)
    theta = np.ones((batch, dim))
    for idx, (i, j) in enumerate(edges):
        # theta_i K[i,j] + theta_j K[j,i] = 0 with K[j,i]=fwd, K[i,j]=bwd
        theta[:, j] = -theta[:, i] * bwd[:, idx] / fwd[:, idx]
    mod = rng.uniform(-4.0, 4.0, size=(batch, 1, len(edges)))
    ksq = np.array([2.0, 6.0, 9.0, 14.0])
    nu = 0.3
    gen = build_generators(edges, fwd, bwd, mod, mod, ksq, nu)[:, 0]
    kmat = gen + nu * np.diag(ksq)[None]
    theta_mat = theta[:, :, None] * np.eye(dim)[None]
    resid = theta_mat @ kmat + np.transpose(kmat, (0, 2, 1)) @ theta_mat
    scale = np.maximum(np.abs(theta_mat @ kmat).max(axis=(1, 2)), 1e-300)
    sym = theta_mat @ gen + np.transpose(gen, (0, 2, 1)) @ theta_mat
    # generalized eigenvalues of sym vs -2 Theta: worst decay rate
    inv_sqrt = 1.0 / np.sqrt(theta)
    whitened = inv_sqrt[:, :, None] * sym * inv_sqrt[:, None, :]
    rates = np.max(np.linalg.eigvalsh((whitened + np.transpose(whitened, (0, 2, 1))) / 2.0), axis=-1)
    return {
        "graph": "path_4node",
        "samples": batch,
        "max_relative_residual_ThetaK_plus_KTTheta": float(
            np.max(np.abs(resid).max(axis=(1, 2)) / scale)
        ),
        "certificate_holds": bool(np.max(np.abs(resid).max(axis=(1, 2)) / scale) < 1e-12),
        "worst_instantaneous_rate": float(np.max(rates)),
        "required_bound_minus_2nu_min_ksq": float(-2.0 * nu * ksq.min()),
        "rate_bound_respected": bool(
            np.all(rates <= -2.0 * nu * ksq.min() + 1e-9)
        ),
        "max_kappa_theta": float(np.max(theta.max(axis=1) / theta.min(axis=1))),
        "max_euclidean_overshoot": float(
            np.sqrt(np.max(theta.max(axis=1) / theta.min(axis=1)))
        ),
    }


def trace_identity_check() -> dict[str, Any]:
    """sum_i chi_i = -nu tr D exactly (det M = exp(-nu tr D T)).

    Consequence: chi_max >= -nu * mean_i |k_i|^2 for ANY zero-diagonal K.
    Combined with Prop 1 (chi_max <= -nu min|k|^2) this PINS chi_max =
    -nu|k|^2 on a single shell.  FLOAT verification.
    """
    rng = np.random.default_rng(7)
    edges = GRAPHS["path_3node"][1]
    batch, steps = 20000, 3
    fwd, bwd = sample_betas(rng, batch, 2, elliptic=True, mag_lo=0.05, mag_hi=30.0)
    mod, _ = sample_modulation(rng, batch, steps, 2, "class_c")
    durations = _loguniform(rng, (batch, steps), 0.05, 3.0)
    ksq = np.array([2.0, 6.0, 9.0])
    nu = 0.5
    gen = build_generators(edges, fwd, bwd, mod, mod, ksq, nu)
    mono, logs = monodromy(gen, durations)
    period = durations.sum(axis=1)
    dim = ksq.size
    logdet = np.log(np.maximum(np.abs(np.linalg.det(mono)), 1e-300)) + dim * logs
    predicted = -nu * ksq.sum() * period
    eigs = np.linalg.eigvals(mono)
    chi_max = (np.log(np.max(np.abs(eigs), axis=-1)) + logs) / period
    return {
        "identity": "log|det M| = -nu * tr(D) * T",
        "max_abs_relative_error": float(
            np.max(np.abs((logdet - predicted) / predicted))
        ),
        "chi_max_min_observed": float(np.min(chi_max)),
        "lower_bound_mean_shell": float(-nu * ksq.mean()),
        "lower_bound_respected": bool(np.all(chi_max >= -nu * ksq.mean() - 1e-9)),
        "upper_bound_min_shell": float(-nu * ksq.min()),
        "upper_bound_respected": bool(np.all(chi_max <= -nu * ksq.min() + MARGIN_TOL)),
    }


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seed", type=int, default=20260802)
    parser.add_argument("--batch", type=int, default=60000)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args(argv)

    started = time.time()
    rng = np.random.default_rng(args.seed)
    batch = 4000 if args.quick else args.batch

    randomized: list[dict[str, Any]] = []
    forest_graphs = ["edge_2node", "path_3node", "star_4node", "path_4node"]
    for graph in forest_graphs:
        for modulation_class in ("class_a", "class_b", "class_c"):
            for steps in (2, 3, 4):
                for shell_mode in ("single_shell", "mixed_shell"):
                    randomized.append(
                        run_family(
                            rng,
                            graph,
                            modulation_class,
                            elliptic=True,
                            balanced_cycle=False,
                            steps=steps,
                            batch=batch,
                            nu=0.5 if shell_mode == "mixed_shell" else 0.05,
                            shell_mode=shell_mode,
                            mag_lo=0.02,
                            mag_hi=50.0,
                            dur_lo=0.02,
                            dur_hi=6.0,
                        )
                    )

    # control families: hypotheses deliberately broken
    controls: list[dict[str, Any]] = []
    for graph in ("edge_2node", "path_3node"):
        for steps in (2, 3, 4):
            controls.append(
                run_family(
                    rng,
                    graph,
                    "class_d",
                    elliptic=True,
                    balanced_cycle=False,
                    steps=steps,
                    batch=batch,
                    nu=0.05,
                    shell_mode="single_shell",
                    mag_lo=0.02,
                    mag_hi=50.0,
                    dur_lo=0.02,
                    dur_hi=6.0,
                )
            )
        for steps in (2, 3):
            controls.append(
                run_family(
                    rng,
                    graph,
                    "class_c",
                    elliptic=False,
                    balanced_cycle=False,
                    steps=steps,
                    batch=batch,
                    nu=0.05,
                    shell_mode="single_shell",
                    mag_lo=0.02,
                    mag_hi=50.0,
                    dur_lo=0.02,
                    dur_hi=6.0,
                )
            )
    for modulation_class in ("class_a", "class_c"):
        for steps in (3, 4):
            controls.append(
                run_family(
                    rng,
                    "cycle_3node",
                    modulation_class,
                    elliptic=True,
                    balanced_cycle=False,
                    steps=steps,
                    batch=batch,
                    nu=0.05,
                    shell_mode="single_shell",
                    mag_lo=0.02,
                    mag_hi=50.0,
                    dur_lo=0.02,
                    dur_hi=6.0,
                )
            )
            controls.append(
                run_family(
                    rng,
                    "cycle_3node",
                    modulation_class,
                    elliptic=True,
                    balanced_cycle=True,
                    steps=steps,
                    batch=batch,
                    nu=0.05,
                    shell_mode="single_shell",
                    mag_lo=0.02,
                    mag_hi=50.0,
                    dur_lo=0.02,
                    dur_hi=6.0,
                )
            )

    exhaustive: list[dict[str, Any]] = []
    if not args.quick:
        wide_beta = (0.25, 1.0, 4.0, 16.0)
        wide_mod = (-2.0, -0.5, 0.0, 1.0, 3.0)
        wide_dur = (0.1, 0.7, 3.0)
        for steps in (2, 3, 4):
            exhaustive.append(
                run_exhaustive_grid(
                    "edge_2node",
                    "class_a",
                    steps=steps,
                    nu=0.05,
                    elliptic=True,
                    beta_mag=wide_beta,
                    mod_vals=wide_mod,
                    dur_vals=wide_dur,
                )
            )
        # 3-node path: coarser per-axis grid, fully independent per-edge
        # modulation (class c), still an exhaustive cartesian product.
        path_beta = (0.25, 1.0, 4.0)
        path_mod = (-2.0, 0.0, 1.0, 3.0)
        path_dur = (0.3, 2.0)
        for steps in (2, 3):
            exhaustive.append(
                run_exhaustive_grid(
                    "path_3node",
                    "class_c",
                    steps=steps,
                    nu=0.05,
                    elliptic=True,
                    beta_mag=path_beta,
                    mod_vals=path_mod,
                    dur_vals=path_dur,
                )
            )
        # 3-cycle at the same grid: the unbalanced cycle MUST break the bound.
        exhaustive.append(
            run_exhaustive_grid(
                "cycle_3node",
                "class_c",
                steps=3,
                nu=0.05,
                elliptic=True,
                beta_mag=(0.5, 2.0),
                mod_vals=(0.0, 1.0, 3.0),
                dur_vals=(0.3, 2.0),
            )
        )

    report = {
        "schema": SCHEMA,
        "audited_claim": (
            "Proposition 1 of idea_parametric_resonance.md: for the child "
            "system xdot = (K(t) - nu D) x with forest coupling graph, the top "
            "Floquet exponent is <= -nu min_i |k_i|^2 for ALL admissible "
            "modulations."
        ),
        "arithmetic_labels": {
            "exact_fraction": [
                "cycle invariant C",
                "edge ratios r_e",
                "det K",
                "c1 = sum omega_e^2",
                "diagonal metric Theta",
                "gauge-consistency check on the B.5 gadget",
            ],
            "float64": [
                "all matrix exponentials / monodromies",
                "spectral radii and Floquet exponents",
                "delta scans and sigma_cycle (cubic root)",
            ],
        },
        "margin_definition": "margin = log rho(M)/T + nu*min_i|k_i|^2 ; >0 kills Prop 1",
        "margin_tolerance": MARGIN_TOL,
        "randomized_forest_families": randomized,
        "control_families": controls,
        "exhaustive_grids": exhaustive,
        "balanced_null_control": balanced_null_control(),
        "theta_certificate": theta_certificate_check(rng),
        "trace_identity": trace_identity_check(),
        "b5_b6_gadget": gadget_report(),
        "holonomy_stress": holonomy_stress(rng, 400),
        "runtime": {
            "seconds": None,
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "platform": platform.platform(),
            "scipy": "absent (numpy-only expm)",
            "seed": args.seed,
            "batch": batch,
        },
    }

    max_forest = max(r["max_margin"] for r in randomized)
    forest_grids = [g for g in exhaustive if g["topology"] == "forest"]
    cycle_grids = [g for g in exhaustive if g["topology"] == "cycle"]
    max_exhaustive = max((g["max_margin"] for g in forest_grids), default=-math.inf)
    report["summary"] = {
        "forest_families_evaluated": len(randomized),
        "forest_samples_total": sum(r["samples"] for r in randomized),
        "exhaustive_evaluations_total": sum(g["evaluations"] for g in exhaustive),
        "exhaustive_forest_evaluations": sum(g["evaluations"] for g in forest_grids),
        "max_margin_over_all_forest_families": max_forest,
        "max_margin_over_exhaustive_forest_grids": max_exhaustive,
        "max_margin_over_exhaustive_cycle_grid": max(
            (g["max_margin"] for g in cycle_grids), default=None
        ),
        "prop1_counterexample_found_under_H1_H2": bool(
            max(max_forest, max_exhaustive) > MARGIN_TOL
        ),
        "control_H1_broken_max_margin": max(
            (c["max_margin"] for c in controls if c["modulation_class"] == "class_d"),
            default=None,
        ),
        "control_H2_broken_max_margin": max(
            (
                c["max_margin"]
                for c in controls
                if not c["hypotheses"]["H2_elliptic_edges"]
            ),
            default=None,
        ),
        "control_unbalanced_cycle_max_margin": max(
            (
                c["max_margin"]
                for c in controls
                if c["graph"] == "cycle_3node"
                and not c["hypotheses"]["H3_balanced_cycles"]
            ),
            default=None,
        ),
        "control_balanced_cycle_max_margin": max(
            (
                c["max_margin"]
                for c in controls
                if c["graph"] == "cycle_3node" and c["hypotheses"]["H3_balanced_cycles"]
            ),
            default=None,
        ),
    }
    report["runtime"]["seconds"] = round(time.time() - started, 2)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / "parametric_common_lyapunov_audit.json"
    target.write_text(json.dumps(report, indent=2, sort_keys=False), encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(f"\nwrote {target}")
    print(f"runtime {report['runtime']['seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
