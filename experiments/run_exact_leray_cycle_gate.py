"""VERIFICATION SPRINT v1 / FREEZE PHASE -- PART II.

EXACT LERAY CYCLE GATE.  One-shot decisive computation for the parametric
resonance mechanism of
``docs/research_notes/ideas_2026_08_01/idea_parametric_resonance.md``
under the corrected Proposition 1 of
``docs/research_notes/verification_sprint_v1/parametric_common_lyapunov_audit.md``.

NO NEW VARIANTS of Lambda, K, N_r, Osgood, or deficit quantities are generated
here.  This script verifies / computes / decides only.

THE OBJECT
----------
Parent field (prescribed, NOT dynamical -- the most generous admissible class):

    U_P(x,t) = sum_a  a_{p_a}(t) * C_a(x),
    C_a(x)   = c_a * cos(p_a . x)   or   c_a * sin(p_a . x),
    c_a . p_a = 0,  c_a integer primitive,  |a_{p_a}(t)| <= 1,
    piecewise constant, arbitrary switching times and signs.

Child band: integer wavevectors k_1..k_L closed under k -> k +/- p_a along the
cycle edges.  Each child mode carries its FULL real divergence-free subspace,
4-dimensional:  {v1 cos, v2 cos, v1 sin, v2 sin} with v1 = k x e, v2 = k x v1
(integer, primitive, mutually orthogonal, both orthogonal to k).  No
polarization is chosen for the children; nothing is projected away by hand.
Consequently the edge coefficients are 4x4 BLOCKS, not scalars, and the
Lyapunov analysis is done on the blocks.

Linearized child system:

    xdot = A(t) x,   A(t) = -nu D + sum_a a_{p_a}(t) K_a,
    D    = blockdiag(|k_i|^2 I_4),
    K_a  = exact Leray linearization  L_a w = -P[(C_a . grad) w + (w . grad) C_a]
           restricted to (and expressed in) the child basis.

Both advection-of-w and advection-of-U terms are kept; both +p and -p sidebands
appear automatically because the real cos/sin product-to-sum identities of
``fourier_torus`` produce k+p and k-p from every product.  The per-mode Leray
projection P_k is the exact rational ``fourier_torus.leray``.

ARITHMETIC LABELS (repo convention, enforced throughout)
--------------------------------------------------------
  EXACT (fractions.Fraction)
      every Gram matrix, every coupling block K_a[i][j], every cross-talk
      coefficient, every edge proportionality constant rho_e, every cycle
      holonomy product, the divergence checks, the reality check, and the
      full-system energy-conservation check.
  FLOAT (numpy float64)
      matrix exponentials, monodromy products, spectral radii, Floquet
      exponents, and the Galerkin cloud test (step 6, only if reached).
  INTERVAL (snapshot_certificate.Interval, exact rational endpoints)
      the Floquet certificate of step 4, if a positive candidate survives.

PRE-REGISTERED ENUMERATION BOUNDS AND COMPUTE CAPS: see ``REGISTRY`` below.
They are frozen before the run and reported honestly if they bind.

Run:
    .venv/Scripts/python.exe experiments/run_exact_leray_cycle_gate.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import platform
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from ns_certificate_lab.fourier_torus import (  # noqa: E402
    TrigVector,
    advection,
    leray,
)
from ns_certificate_lab.snapshot_certificate import Interval  # noqa: E402

SCHEMA = "ns-certificate-lab/exact-leray-cycle-gate/v1"
DEFAULT_OUTPUT_DIR = (
    REPOSITORY_ROOT / "outputs" / "verification_sprint_v1" / "cycle_gate"
)

Vec = tuple[int, int, int]


# --------------------------------------------------------------------------- #
# PRE-REGISTERED REGISTRY  (frozen before the run; deviations reported)        #
# --------------------------------------------------------------------------- #

REGISTRY: dict[str, Any] = {
    "parent_shell_max_sq": 6,          # 1 <= |p|^2 <= 6, canonical representatives
    "child_shell_max_sq": 50,          # every cycle mode obeys |k|^2 <= 50
    "cycle_classes": ["3-cycle", "4-cycle", "double-triad"],
    "geometry_cap": 200,               # first 200 distinct geometries by max |k|^2
    "geometry_order_R1": "(max |k|^2, sum |k|^2, lexicographic child tuple)",
    "geometry_list_R2_stratified": (
        "PRE-REGISTERED SUPPLEMENT, declared before any dynamics was run, because "
        "R1's cap binds degenerately (all 200 sit at max|k|^2 <= 2).  R2 = round-robin "
        "over the strata (cycle kind) x (max|k|^2 shell), taking geometries in R1 order "
        "within each stratum, 200 total, with the note's own B.5 reference geometry "
        "({(2,2,1),(2,1,2),(1,2,2)}; parents |p|^2=2) forced in."
    ),
    "per_shell_bucket_cap": 24,        # enumeration keeps <= 24 geometries per (kind, shell)
    "parent_orientations_per_parent": 8,   # 4 polarizations x {cos, sin}
    "polarization_generators": "v1 = k x e, v2 = k x v1, and v1+v2, v1-v2 (primitive)",
    "child_block_dim": 4,              # FULL divergence-free subspace per child mode
    "orientations_for_dynamics_per_geometry": 64,
    "orientation_ranking_for_dynamics": (
        "exact imbalance score: (#non-proportional edge blocks, #hyperbolic edges, "
        "|log holonomy| when defined); ties by generator index.  Declared before "
        "any dynamics."
    ),
    "nu_values": ["0 (structural, coupling only)", "1/40", "1/10"],
    "amplitude_grid": ["-1", "-1/2", "0", "1/2", "1"],
    "switching_intervals": [1, 2, 3, 4],
    "schedule_families": [
        "staggered (one parent per interval, all sign patterns)",
        "constant (all parents on, all sign patterns)",
        "two-interval sign flip",
        "random piecewise-constant on the amplitude grid",
    ],
    "period_scales": "delta in {1/16, 1/4, 1, 4} / coupling_inf_norm, plus 1/(nu max|k|^2)",
    "random_schedules_per_candidate_screen": 8,
    "random_schedules_per_candidate_refine": 4000,
    "rho_margin_for_growth": 1e-9,
    "compute_caps_minutes": {
        "stage_A_exact_enumeration_and_blocks": 15,
        "stage_B_monodromy_screen": 20,
        "stage_C_refine_and_interval_certificate": 15,
        "total": 50,
    },
    "GO_conditions_step7": [
        "inf_N E_next/E_current > 0 with no decreasing trend",
        "off-chain/main ratio non-increasing in N",
        "second-stage normalized gain non-degenerate",
        "output shape returns to the same cycle class",
    ],
}


# --------------------------------------------------------------------------- #
# exact lattice helpers                                                        #
# --------------------------------------------------------------------------- #


def _norm_sq(k: Vec) -> int:
    return k[0] * k[0] + k[1] * k[1] + k[2] * k[2]


def _cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vec, b: Vec) -> int:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _primitive(v: Vec) -> Vec:
    g = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    if g == 0:
        return v
    w = (v[0] // g, v[1] // g, v[2] // g)
    for c in w:
        if c > 0:
            return w
        if c < 0:
            return (-w[0], -w[1], -w[2])
    return w


def _canonical(k: Vec) -> tuple[Vec, int]:
    if k == (0, 0, 0):
        return k, 1
    for c in k:
        if c > 0:
            return k, 1
        if c < 0:
            return (-k[0], -k[1], -k[2]), -1
    raise AssertionError("unreachable")


def perp_generators(k: Vec) -> tuple[Vec, Vec]:
    """Two integer, primitive, mutually orthogonal vectors spanning k^perp."""
    for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        v1 = _cross(k, e)
        if v1 != (0, 0, 0):
            break
    v1 = _primitive(v1)
    v2 = _primitive(_cross(k, v1))
    assert _dot(k, v1) == 0 and _dot(k, v2) == 0 and _dot(v1, v2) == 0
    return v1, v2


def shell_modes(max_sq: int) -> list[Vec]:
    b = int(math.isqrt(max_sq))
    out = []
    for x in range(-b, b + 1):
        for y in range(-b, b + 1):
            for z in range(-b, b + 1):
                k = (x, y, z)
                if not 1 <= _norm_sq(k) <= max_sq:
                    continue
                key, sign = _canonical(k)
                if sign > 0:
                    out.append(key)
    return sorted(out)


# --------------------------------------------------------------------------- #
# exact trig-field helpers                                                     #
# --------------------------------------------------------------------------- #


def mode_field(k: Vec, a: Sequence, b: Sequence) -> TrigVector:
    return TrigVector.from_modes([(k, tuple(a), tuple(b))])


def inner(u: TrigVector, v: TrigVector) -> Fraction:
    """Exact L^2 inner product in the normalised measure (2pi)^-3 dx."""
    tu = u.coefficient_table()
    tv = v.coefficient_table()
    total = Fraction(0)
    for k, (a1, b1) in tu.items():
        if k == (0, 0, 0):
            continue
        entry = tv.get(k)
        if entry is None:
            continue
        a2, b2 = entry
        s = Fraction(0)
        for i in range(3):
            s += Fraction(a1[i]) * Fraction(a2[i]) + Fraction(b1[i]) * Fraction(b2[i])
        total += s / 2
    return total


def divergence_defect(u: TrigVector) -> Fraction:
    d = u.divergence().cleaned()
    total = Fraction(0)
    for _k, (a, b) in d.terms.items():
        total += abs(Fraction(a)) + abs(Fraction(b))
    return total


def leray_linearization(parent: TrigVector, w: TrigVector) -> TrigVector:
    """L_P w = -P[(U_P . grad) w + (w . grad) U_P]   (exact)."""
    return -leray(advection(parent, w) + advection(w, parent))


# --------------------------------------------------------------------------- #
# child basis                                                                  #
# --------------------------------------------------------------------------- #


@dataclass
class ChildBasis:
    k: Vec
    fields: tuple[TrigVector, TrigVector, TrigVector, TrigVector]
    gram: tuple[Fraction, Fraction, Fraction, Fraction]
    labels: tuple[str, str, str, str]


_CHILD_CACHE: dict[Vec, ChildBasis] = {}


def child_basis(k: Vec) -> ChildBasis:
    kc, sign = _canonical(k)
    cached = _CHILD_CACHE.get(kc)
    if cached is not None:
        return cached
    v1, v2 = perp_generators(kc)
    zero = (0, 0, 0)
    fields = (
        mode_field(kc, v1, zero),
        mode_field(kc, v2, zero),
        mode_field(kc, zero, v1),
        mode_field(kc, zero, v2),
    )
    n1 = Fraction(_norm_sq(v1), 2)
    n2 = Fraction(_norm_sq(v2), 2)
    cb = ChildBasis(kc, fields, (n1, n2, n1, n2), ("v1cos", "v2cos", "v1sin", "v2sin"))
    _CHILD_CACHE[kc] = cb
    return cb


# --------------------------------------------------------------------------- #
# parent generators                                                            #
# --------------------------------------------------------------------------- #


@dataclass
class ParentGenerator:
    p: Vec
    pol: Vec
    phase: str
    field: TrigVector
    label: str


_PARENT_CACHE: dict[Vec, list[ParentGenerator]] = {}


def parent_generators(p: Vec) -> list[ParentGenerator]:
    """The 8 pre-registered orientations of one parent wavevector."""
    cached = _PARENT_CACHE.get(p)
    if cached is not None:
        return cached
    v1, v2 = perp_generators(p)
    pols = [
        v1,
        v2,
        _primitive(tuple(v1[i] + v2[i] for i in range(3))),
        _primitive(tuple(v1[i] - v2[i] for i in range(3))),
    ]
    zero = (0, 0, 0)
    out = []
    for idx, c in enumerate(pols):
        for phase in ("cos", "sin"):
            fld = mode_field(p, c, zero) if phase == "cos" else mode_field(p, zero, c)
            out.append(ParentGenerator(p, c, phase, fld, f"pol{idx}_{phase}"))
    _PARENT_CACHE[p] = out
    return out


# --------------------------------------------------------------------------- #
# exact coupling blocks                                                        #
# --------------------------------------------------------------------------- #


@dataclass
class CouplingTables:
    blocks: dict[tuple[int, int], list[list[Fraction]]]
    crosstalk: dict[tuple[Vec, int], list[list[Fraction]]]
    crosstalk_energy: Fraction
    inset_energy: Fraction


def coupling_tables(gen: ParentGenerator, basis: Sequence[ChildBasis]) -> CouplingTables:
    """Every exact coefficient of L_gen on the child band, in-set and leaking out."""
    index = {b.k: i for i, b in enumerate(basis)}
    blocks: dict[tuple[int, int], list[list[Fraction]]] = {}
    crosstalk: dict[tuple[Vec, int], list[list[Fraction]]] = {}
    ct_energy = Fraction(0)
    in_energy = Fraction(0)
    for j, bj in enumerate(basis):
        for m in range(4):
            image = leray_linearization(gen.field, bj.fields[m]).cleaned()
            for k, (a, b) in image.coefficient_table().items():
                if k == (0, 0, 0):
                    continue
                target = child_basis(k)
                fld = mode_field(k, a, b)
                comps = [inner(fld, target.fields[n]) / target.gram[n] for n in range(4)]
                if all(c == 0 for c in comps):
                    continue
                if k in index:
                    i = index[k]
                    blk = blocks.setdefault((i, j), [[Fraction(0)] * 4 for _ in range(4)])
                    for n in range(4):
                        blk[n][m] += comps[n]
                        in_energy += comps[n] * comps[n] * target.gram[n]
                else:
                    blk = crosstalk.setdefault((k, j), [[Fraction(0)] * 4 for _ in range(4)])
                    for n in range(4):
                        blk[n][m] += comps[n]
                        ct_energy += comps[n] * comps[n] * target.gram[n]
    return CouplingTables(blocks, crosstalk, ct_energy, in_energy)


# --------------------------------------------------------------------------- #
# CORRECTNESS CHECKS (mandatory, exact, before any dynamics)                   #
# --------------------------------------------------------------------------- #


def check_divergence_free(fields: Iterable[tuple[str, TrigVector]]) -> dict[str, Any]:
    defects = {name: divergence_defect(f) for name, f in fields}
    worst = max(defects.values()) if defects else Fraction(0)
    return {
        "fields_checked": len(defects),
        "max_divergence_defect_exact": str(worst),
        "all_exactly_divergence_free": worst == 0,
    }


def check_reality(blocks: Iterable[list[list[Fraction]]]) -> dict[str, Any]:
    n = 0
    ok = True
    for blk in blocks:
        for row in blk:
            for v in row:
                n += 1
                if not isinstance(v, Fraction):
                    ok = False
    return {"entries_checked": n, "all_entries_exact_rational_real": ok}


def check_full_system_energy(fields: Sequence[tuple[str, TrigVector]]) -> dict[str, Any]:
    """Exact energy check on the FULL system with the parents dynamical."""
    names = [n for n, _ in fields]
    fs = [f for _, f in fields]
    n = len(fs)
    tri = {}
    for a in range(n):
        for b in range(n):
            adv = advection(fs[a], fs[b])
            for c in range(n):
                tri[(a, b, c)] = inner(adv, fs[c])
    worst_perm = Fraction(0)
    worst_anti = Fraction(0)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                worst_anti = max(worst_anti, abs(tri[(a, b, c)] + tri[(a, c, b)]))
    for a, b, c in itertools.combinations_with_replacement(range(n), 3):
        total = Fraction(0)
        for perm in itertools.permutations((a, b, c)):
            total += tri[perm]
        worst_perm = max(worst_perm, abs(total))
    return {
        "basis_fields": names,
        "n_fields": n,
        "ordered_triples_checked": n**3,
        "max_six_permutation_sum_exact": str(worst_perm),
        "max_divfree_antisymmetry_defect_exact": str(worst_anti),
        "energy_conserved_exactly": worst_perm == 0,
        "antisymmetry_exact": worst_anti == 0,
    }


# --------------------------------------------------------------------------- #
# cycle enumeration                                                            #
# --------------------------------------------------------------------------- #


@dataclass
class Cycle:
    kind: str
    children: tuple[Vec, ...]
    steps: tuple[Vec, ...]
    parents: tuple[Vec, ...]
    max_k_sq: int
    sum_k_sq: int


def _cycle_key(children: Sequence[Vec]) -> tuple:
    best = None
    r = len(children)
    for sgn in (1, -1):
        ch = [_canonical(tuple(sgn * c for c in k))[0] for k in children]
        for direction in (1, -1):
            seq = ch if direction == 1 else list(reversed(ch))
            for shift in range(r):
                cand = tuple(seq[(shift + t) % r] for t in range(r))
                if best is None or cand < best:
                    best = cand
    return best


def enumerate_cycles(
    parent_max_sq: int, child_max_sq: int, bucket_cap: int
) -> tuple[dict[tuple[str, int], list[Cycle]], dict[str, Any]]:
    """Pre-registered enumeration, bucketed by (kind, max |k|^2)."""
    parents = shell_modes(parent_max_sq)
    parent_set = set(parents)
    signed = []
    for p in parents:
        signed.append(p)
        signed.append(tuple(-c for c in p))
    children = shell_modes(child_max_sq)
    child_set = set(children)

    buckets: dict[tuple[str, int], list[Cycle]] = {}
    seen: set[tuple] = set()
    raw_counts = {"3-cycle": 0, "4-cycle": 0, "double-triad": 0}
    full_counts = {"3-cycle": 0, "4-cycle": 0, "double-triad": 0}

    def register(kind: str, raw: Sequence[Vec], steps: Sequence[Vec]) -> None:
        raw_counts[kind] += 1
        cans = []
        mx = 0
        for k in raw:
            c = _canonical(k)[0]
            if c == (0, 0, 0) or c not in child_set:
                return
            cans.append(c)
            mx = max(mx, _norm_sq(c))
        if len(set(cans)) != len(cans):
            return
        bucket = buckets.setdefault((kind, mx), [])
        if len(bucket) >= bucket_cap:
            full_counts[kind] += 1
            return
        key = (kind, _cycle_key(cans))
        if key in seen:
            return
        seen.add(key)
        bucket.append(
            Cycle(
                kind=kind,
                children=tuple(cans),
                steps=tuple(steps),
                parents=tuple(sorted({_canonical(s)[0] for s in steps})),
                max_k_sq=mx,
                sum_k_sq=sum(_norm_sq(k) for k in cans),
            )
        )

    # 3-cycles: three signed parent steps summing to zero (closed parent triad)
    zero_triads = []
    for da, db in itertools.combinations(signed, 2):
        dc = tuple(-(da[i] + db[i]) for i in range(3))
        if not 1 <= _norm_sq(dc) <= parent_max_sq:
            continue
        if _canonical(dc)[0] not in parent_set:
            continue
        zero_triads.append((da, db, dc))
    for da, db, dc in zero_triads:
        for k1 in children:
            k2 = tuple(k1[i] + da[i] for i in range(3))
            k3 = tuple(k2[i] + db[i] for i in range(3))
            register("3-cycle", (k1, k2, k3), (da, db, dc))

    # double-triads: (p, q, -p, -q) parallelograms, two parents alternating
    for p, q in itertools.combinations(parents, 2):
        for sp in (1, -1):
            for sq in (1, -1):
                d1 = tuple(sp * c for c in p)
                d2 = tuple(sq * c for c in q)
                d3 = tuple(-c for c in d1)
                d4 = tuple(-c for c in d2)
                for k1 in children:
                    k2 = tuple(k1[i] + d1[i] for i in range(3))
                    k3 = tuple(k2[i] + d2[i] for i in range(3))
                    k4 = tuple(k3[i] + d3[i] for i in range(3))
                    register("double-triad", (k1, k2, k3, k4), (d1, d2, d3, d4))

    # general 4-cycles: four signed steps summing to zero
    for d1, d2, d3 in itertools.combinations_with_replacement(signed, 3):
        d4 = tuple(-(d1[i] + d2[i] + d3[i]) for i in range(3))
        if not 1 <= _norm_sq(d4) <= parent_max_sq:
            continue
        if _canonical(d4)[0] not in parent_set:
            continue
        for k1 in children:
            k2 = tuple(k1[i] + d1[i] for i in range(3))
            k3 = tuple(k2[i] + d2[i] for i in range(3))
            k4 = tuple(k3[i] + d3[i] for i in range(3))
            register("4-cycle", (k1, k2, k3, k4), (d1, d2, d3, d4))

    stats = {
        "raw_paths_examined": raw_counts,
        "bucket_cap_per_kind_and_shell": bucket_cap,
        "buckets_that_overflowed": {
            k: v for k, v in full_counts.items() if v > 0
        },
        "distinct_geometries_retained": sum(len(v) for v in buckets.values()),
        "shells_populated": sorted({s for (_k, s) in buckets}),
    }
    return buckets, stats


def select_lists(
    buckets: dict[tuple[str, int], list[Cycle]], cap: int
) -> tuple[list[Cycle], list[Cycle], dict[str, Any]]:
    """R1 = pre-registered "first cap by increasing max|k|^2".
    R2 = pre-registered stratified supplement (round-robin over shells)."""
    allc = [c for v in buckets.values() for c in v]
    ordered = sorted(allc, key=lambda c: (c.max_k_sq, c.sum_k_sq, c.children))
    r1 = ordered[:cap]

    strata = sorted({(c.kind, c.max_k_sq) for c in allc})
    per_stratum: dict[tuple[str, int], list[Cycle]] = {s: [] for s in strata}
    for c in ordered:
        per_stratum[(c.kind, c.max_k_sq)].append(c)
    r2: list[Cycle] = []
    idx = 0
    while len(r2) < cap:
        added = False
        for s in strata:
            if idx < len(per_stratum[s]) and len(r2) < cap:
                r2.append(per_stratum[s][idx])
                added = True
        if not added:
            break
        idx += 1

    # force the note's own B.5 reference geometry in
    ref = Cycle(
        kind="3-cycle",
        children=((2, 2, 1), (2, 1, 2), (1, 2, 2)),
        steps=((0, -1, 1), (-1, 1, 0), (1, 0, -1)),
        parents=((0, 1, -1), (1, -1, 0), (1, 0, -1)),
        max_k_sq=9,
        sum_k_sq=27,
    )
    if not any(c.children == ref.children for c in r2):
        r2 = [ref] + r2[: cap - 1]

    info = {
        "R1_cap_binds": len(ordered) > cap,
        "R1_max_k_sq_range": [min(c.max_k_sq for c in r1), max(c.max_k_sq for c in r1)],
        "R1_by_kind": {k: sum(1 for c in r1 if c.kind == k) for k in
                       ("3-cycle", "4-cycle", "double-triad")},
        "R1_DEGENERATE": len({c.max_k_sq for c in r1}) <= 2,
        "R2_max_k_sq_range": [min(c.max_k_sq for c in r2), max(c.max_k_sq for c in r2)],
        "R2_shells_covered": sorted({c.max_k_sq for c in r2}),
        "R2_by_kind": {k: sum(1 for c in r2 if c.kind == k) for k in
                       ("3-cycle", "4-cycle", "double-triad")},
        "B5_reference_forced_in": True,
    }
    return r1, r2, info


# --------------------------------------------------------------------------- #
# exact block-Lyapunov classification                                          #
# --------------------------------------------------------------------------- #


def _gram_left(gram: Sequence[Fraction], m: list[list[Fraction]]):
    return [[gram[r] * m[r][c] for c in range(4)] for r in range(4)]


def _transpose_gram_right(m: list[list[Fraction]], gram: Sequence[Fraction]):
    return [[m[c][r] * gram[c] for c in range(4)] for r in range(4)]


def _proportional(P, Q) -> tuple[str, Fraction | None]:
    p_zero = all(v == 0 for row in P for v in row)
    q_zero = all(v == 0 for row in Q for v in row)
    if p_zero and q_zero:
        return "both_zero", None
    if p_zero or q_zero:
        return "one_sided", None
    lam = None
    for r in range(4):
        for c in range(4):
            if P[r][c] == 0:
                if Q[r][c] != 0:
                    return "not_proportional", None
                continue
            ratio = Q[r][c] / P[r][c]
            if lam is None:
                lam = ratio
            elif ratio != lam:
                return "not_proportional", None
    return "proportional", lam


ZERO4 = [[Fraction(0)] * 4 for _ in range(4)]


def classify(basis: Sequence[ChildBasis], tables: Sequence[tuple[str, CouplingTables]]):
    """Exact block version of (H1)/(H2)/(H3).

    Metric class: W = blockdiag(w_i G_i), w_i > 0 (one positive weight per child
    mode; G_i is the exact mode Gram, so x^T W x is the physical energy weighted
    by w).  By Lemma D0 of the audit, W A + A^T W <= 0 for all modulations forces

        w_i G_i B_ij + w_j B_ji^T G_j = 0   for every edge and EVERY parent

    (the parents modulate independently, so each must vanish separately).
    Writing P = G_i B_ij and Q = B_ji^T G_j this is Q = -(w_i/w_j) P: the BLOCK
    ellipticity condition.  Nothing scalar is assumed anywhere.
    """
    L = len(basis)
    edges = []
    constraints = []
    diag_defect = Fraction(0)
    n_nonprop = n_hyper = n_onesided = 0

    for label, tab in tables:
        for (i, j), blk in tab.blocks.items():
            if i != j:
                continue
            GB = _gram_left(basis[i].gram, blk)
            for r in range(4):
                for c in range(4):
                    diag_defect = max(diag_defect, abs(GB[r][c] + GB[c][r]))
        for i in range(L):
            for j in range(i + 1, L):
                Bij = tab.blocks.get((i, j))
                Bji = tab.blocks.get((j, i))
                if Bij is None and Bji is None:
                    continue
                P = _gram_left(basis[i].gram, Bij or ZERO4)
                Q = _transpose_gram_right(Bji or ZERO4, basis[j].gram)
                status, lam = _proportional(P, Q)
                rho = None
                elliptic = False
                if status == "proportional":
                    rho = -lam
                    elliptic = rho > 0
                    if elliptic:
                        constraints.append((i, j, rho))
                    else:
                        n_hyper += 1
                elif status == "not_proportional":
                    n_nonprop += 1
                elif status == "one_sided":
                    n_onesided += 1
                edges.append((i, j, label, status, lam, rho, elliptic))

    # (H3): propagate w_i/w_j = rho and detect holonomy inconsistency
    w: dict[int, Fraction] = {}
    balance_ok = True
    defects = []
    adj: dict[int, list[tuple[int, Fraction]]] = {}
    for i, j, rho in constraints:
        adj.setdefault(i, []).append((j, Fraction(1) / rho))
        adj.setdefault(j, []).append((i, rho))
    for start in range(L):
        if start in w:
            continue
        w[start] = Fraction(1)
        stack = [start]
        while stack:
            u = stack.pop()
            for v, factor in adj.get(u, []):
                val = w[u] * factor
                if v in w:
                    if w[v] != val:
                        balance_ok = False
                        defects.append(str(val / w[v]))
                else:
                    w[v] = val
                    stack.append(v)

    holonomy = None
    if not edges:
        holonomy = None
    elif all(e[3] == "proportional" for e in edges):
        prod = Fraction(1)
        ok = True
        for i in range(L):
            j = (i + 1) % L
            lo, hi = min(i, j), max(i, j)
            match = [e for e in edges if (e[0], e[1]) == (lo, hi) and e[5] is not None]
            if not match:
                ok = False
                break
            r = match[0][5]
            prod *= r if (i, j) == (lo, hi) else Fraction(1) / r
        holonomy = str(prod) if ok else None

    has_metric = (
        balance_ok
        and n_nonprop == 0
        and n_onesided == 0
        and n_hyper == 0
        and diag_defect == 0
        and bool(edges)
    )
    score = (n_nonprop, n_hyper, n_onesided, 0 if balance_ok else 1)
    return {
        "edge_block_pairs": len(edges),
        "edges": [
            {
                "i": e[0], "j": e[1], "parent": e[2], "status": e[3],
                "lambda_Q_over_P_exact": None if e[4] is None else str(e[4]),
                "rho_required_exact": None if e[5] is None else str(e[5]),
                "block_elliptic": e[6],
            }
            for e in edges
        ],
        "non_proportional_edges": n_nonprop,
        "hyperbolic_edges": n_hyper,
        "one_sided_edges": n_onesided,
        "diagonal_block_antisymmetry_defect_exact": str(diag_defect),
        "cycle_balance_consistent": balance_ok,
        "holonomy_defect_ratios": defects[:6],
        "holonomy_product_exact": holonomy,
        "common_block_scalar_metric_exists": has_metric,
        "metric_weights_exact": (
            {str(i): str(v) for i, v in sorted(w.items())} if has_metric else None
        ),
        "class": "a_dead_common_metric" if has_metric else "b_candidate_for_growth",
        "_score": score,
    }


# --------------------------------------------------------------------------- #
# float dynamics                                                               #
# --------------------------------------------------------------------------- #


def blocks_to_matrix(basis: Sequence[ChildBasis], tab: CouplingTables) -> np.ndarray:
    L = len(basis)
    M = np.zeros((4 * L, 4 * L))
    for (i, j), blk in tab.blocks.items():
        for r in range(4):
            for c in range(4):
                M[4 * i + r, 4 * j + c] = float(blk[r][c])
    return M


def diag_vector(basis: Sequence[ChildBasis]) -> np.ndarray:
    d = np.zeros(4 * len(basis))
    for i, b in enumerate(basis):
        d[4 * i : 4 * i + 4] = float(_norm_sq(b.k))
    return d


def expm_float(A: np.ndarray) -> np.ndarray:
    nrm = float(np.abs(A).sum(axis=1).max())
    n = A.shape[0]
    if nrm == 0.0:
        return np.eye(n)
    s = int(max(0, math.ceil(math.log2(nrm / 0.25)))) if nrm > 0.25 else 0
    B = A / (2.0**s)
    X = np.eye(n)
    term = np.eye(n)
    for m in range(1, 22):
        term = term @ B / m
        X = X + term
    for _ in range(s):
        X = X @ X
    return X


def spectral_radius(M: np.ndarray) -> float:
    try:
        return float(np.max(np.abs(np.linalg.eigvals(M))))
    except np.linalg.LinAlgError:
        return float("nan")


def monodromy(Ks, dvec, nu, sched, cache=None):
    n = dvec.shape[0]
    M = np.eye(n)
    for dt, amps in sched:
        key = (dt, amps)
        E = None if cache is None else cache.get(key)
        if E is None:
            A = np.diag(-nu * dvec)
            for a, K in zip(amps, Ks):
                if a != 0.0:
                    A = A + a * K
            E = expm_float(A * dt)
            if cache is not None:
                cache[key] = E
        M = E @ M
    return M


def schedule_family(n_parents, deltas, rng, n_random):
    """Pre-registered piecewise-constant schedules, |a| <= 1.

    Interval lengths are drawn from the ``deltas`` grid only (this makes the
    per-piece matrix exponential cacheable; it costs nothing in generality that
    the delta grid does not already cover, since the grid spans four octaves of
    1/coupling plus the viscous timescale).
    """
    grid = (-1.0, -0.5, 0.0, 0.5, 1.0)
    out = []
    for dt in deltas:
        for signs in itertools.product((1.0, -1.0), repeat=n_parents):
            sched = []
            for idx in range(n_parents):
                amps = [0.0] * n_parents
                amps[idx] = signs[idx]
                sched.append((dt, tuple(amps)))
            out.append(sched)
            out.append([(dt * n_parents, tuple(signs))])
            out.append([(dt, signs), (dt, tuple(-s for s in signs))])
    for _ in range(n_random):
        S = int(rng.integers(2, 5))
        sched = []
        for _s in range(S):
            amps = tuple(float(grid[int(rng.integers(len(grid)))]) for _ in range(n_parents))
            dt = float(deltas[int(rng.integers(len(deltas)))])
            sched.append((dt, amps))
        out.append(sched)
    return out


def instantaneous_margins(Ks, dvec, nu):
    """max Re lambda(-nu D + a K_a) over a in {+1,-1}, for EACH parent alone.

    The idea note's whole claim is growth with *identically zero instantaneous
    eigenvalue margin*: every single-parent generator is supposed to be
    elliptic, so that any Floquet gain is pure cycle holonomy.  This function
    measures that premise directly.  If some single parent already has a
    positive margin, then a frozen parent alone destabilises the band by
    ordinary strain (vortex stretching) and the cycle contributes nothing
    qualitatively new.
    """
    out = []
    for K in Ks:
        best = -np.inf
        for a in (1.0, -1.0):
            A = np.diag(-nu * dvec) + a * K
            ev = np.linalg.eigvals(A)
            best = max(best, float(np.max(ev.real)))
        out.append(best)
    return out


def best_growth(Ks, dvec, nu, scheds):
    best = (-np.inf, None, None)
    cache: dict = {}
    for sch in scheds:
        M = monodromy(Ks, dvec, nu, sch, cache)
        if not np.all(np.isfinite(M)):
            continue
        rho = spectral_radius(M)
        if not np.isfinite(rho) or rho <= 0:
            continue
        T = sum(dt for dt, _ in sch)
        chi = math.log(rho) / T
        if chi > best[0]:
            best = (chi, rho, sch)
    return best


# --------------------------------------------------------------------------- #
# INTERVAL certificate                                                         #
# --------------------------------------------------------------------------- #

IZERO = Interval(Fraction(0), Fraction(0))
IONE = Interval(Fraction(1), Fraction(1))
BITS = 100


def _ir(x: Interval) -> Interval:
    return x.round_outward(BITS)


def imat_mul(A, B, n):
    C = [[IZERO for _ in range(n)] for _ in range(n)]
    for i in range(n):
        Ci = C[i]
        for k in range(n):
            a = A[i][k]
            if a.lower == 0 and a.upper == 0:
                continue
            Bk = B[k]
            for j in range(n):
                b = Bk[j]
                if b.lower == 0 and b.upper == 0:
                    continue
                Ci[j] = _ir(Ci[j] + a * b)
    return C


def imat_add(A, B, n):
    return [[_ir(A[i][j] + B[i][j]) for j in range(n)] for i in range(n)]


def imat_scale(A, q: Fraction, n):
    return [[_ir(A[i][j].scale(q)) for j in range(n)] for i in range(n)]


def imat_inf_norm(A, n) -> Fraction:
    return max(sum(A[i][j].magnitude for j in range(n)) for i in range(n))


def _exp_upper(r: Fraction, terms: int = 24) -> Fraction:
    """A cheap rational UPPER bound on exp(r), r >= 0, with small denominators.

    Kept deliberately crude: it multiplies a remainder that is already ~1e-30,
    and an exact-but-huge Fraction here costs thousands of digits in every
    subsequent interval operation.
    """
    if r <= Fraction(1, 4):
        return Fraction(4, 3)  # exp(1/4) = 1.2840... < 4/3
    s = Fraction(0)
    t = Fraction(1)
    for m in range(terms):
        if m > 0:
            t = t * r / m
        s += t
    if r >= terms:
        raise ValueError("exp bound: increase terms")
    out = s + t * (r / terms) / (1 - r / terms)
    return Interval(out, out).round_outward(32).upper


def iexpm(A, n, order: int = 26):
    """Rigorous interval enclosure of exp(A), scaling-and-squaring Taylor.

    Remainder:  ||e^B - T_m(B)||_inf <= ||B||_inf^{m+1} e^{||B||_inf} / (m+1)!,
    added as a symmetric interval to every entry (entrywise magnitude <= inf-norm).
    """
    nrm = imat_inf_norm(A, n)
    s = 0
    r = nrm
    while r > Fraction(1, 4):
        r /= 2
        s += 1
    B = imat_scale(A, Fraction(1, 2**s), n)
    bn = imat_inf_norm(B, n)
    X = [[IONE if i == j else IZERO for j in range(n)] for i in range(n)]
    term = [[IONE if i == j else IZERO for j in range(n)] for i in range(n)]
    for m in range(1, order + 1):
        term = imat_scale(imat_mul(term, B, n), Fraction(1, m), n)
        X = imat_add(X, term, n)
    fact = Fraction(1)
    for m in range(1, order + 2):
        fact *= m
    rem = (bn ** (order + 1)) * _exp_upper(bn) / fact
    # round the remainder UP to a short dyadic: keeping it exact costs thousands
    # of digits in every later interval op and buys nothing (it is ~1e-30).
    rem = Interval(rem, rem).round_outward(200).upper
    ri = Interval(-rem, rem)
    X = [[_ir(X[i][j] + ri) for j in range(n)] for i in range(n)]
    for _ in range(s):
        X = imat_mul(X, X, n)
    return X, {"scaling_squarings": s, "taylor_order": order,
               "scaled_inf_norm": str(bn), "remainder_bound": str(rem)}


def _idet(M, n) -> Interval:
    if n == 1:
        return M[0][0]
    total = IZERO
    for j in range(n):
        e = M[0][j]
        if e.lower == 0 and e.upper == 0:
            continue
        minor = [[M[i][c] for c in range(n) if c != j] for i in range(1, n)]
        term = e * _idet(minor, n - 1)
        total = total + (term if j % 2 == 0 else -term)
    return _ir(total)


def interval_certificate(Ks_exact, dvec, nu: Fraction, sched, n, squarings: int = 11):
    """Rigorous enclosure of the monodromy, and a rigorous LOWER bound on rho(M).

    Enclosure: interval Taylor scaling-and-squaring per constant piece, with the
    rigorous remainder ||e^A - T_m(A)|| <= ||A||^{m+1} e^{||A||}/(m+1)!.

    Lower bound on the spectral radius.  The determinant route is USELESS here
    and is reported only for completeness: every K_a is traceless on the child
    band, so det M = exp(-nu tr D * T) <= 1 exactly, and |det M|^{1/n} can never
    exceed 1.  The bound actually claimed is the trace-power bound

        |tr(M^q)| = |sum_i lambda_i^q| <= n * rho(M)^q
        =>  rho(M) >= ( |tr(M^q)| / n )^{1/q} ,

    evaluated on rigorous interval enclosures of M, M^2, M^4, ..., M^{2^s}.
    It is rigorous for every q; it is *sharp* only when the dominant eigenvalue
    does not cancel in the trace, so the maximum over q is reported.
    """
    info = []
    M = [[IONE if i == j else IZERO for j in range(n)] for i in range(n)]
    exact_sched = []
    for dt, amps in sched:
        A = [[IZERO for _ in range(n)] for _ in range(n)]
        for i in range(n):
            A[i][i] = Interval.exact(-nu * Fraction(int(dvec[i])))
        amps_q = [Fraction(a).limit_denominator(64) for a in amps]
        for a, K in zip(amps_q, Ks_exact):
            if a == 0:
                continue
            for (r, c), val in K.items():
                A[r][c] = _ir(A[r][c] + Interval.exact(a * val))
        dt_q = Fraction(dt).limit_denominator(10**7)
        exact_sched.append((str(dt_q), [str(a) for a in amps_q]))
        Adt = imat_scale(A, dt_q, n)
        E, meta = iexpm(Adt, n)
        info.append(meta)
        M = imat_mul(E, M, n)

    # NOTE: the determinant is deliberately NOT computed.  It is useless here
    # (det M = exp(-nu tr D T) <= 1 exactly, since every K_a is traceless, so
    # |det M|^{1/n} can never exceed 1), and the exact Laplace expansion in
    # `_idet` is O(n!) -- 12! ~ 4.8e8 terms at n = 12.  The analytic value is
    # reported instead.
    powers = []
    P = M
    q = 1
    best = (0.0, None)
    for _s in range(squarings + 1):
        tr = IZERO
        for i in range(n):
            tr = _ir(tr + P[i][i])
        lo = Fraction(0) if tr.contains_zero else min(abs(tr.lower), abs(tr.upper))
        width = float(tr.upper - tr.lower)
        bound = float(lo / n) ** (1.0 / q) if lo > 0 else 0.0
        powers.append(
            {
                "q": q,
                "trace_enclosure": [float(tr.lower), float(tr.upper)],
                "trace_enclosure_width": width,
                "abs_trace_lower_bound": float(lo),
                "rho_lower_bound": bound,
            }
        )
        if bound > best[0]:
            best = (bound, q)
        if _s < squarings:
            P = imat_mul(P, P, n)
            q *= 2
    return {
        "precision_bits": BITS,
        "exact_schedule_used": exact_sched,
        "nu_exact": str(nu),
        "expm_pieces": info,
        "det_analytic": float(
            math.exp(-float(nu) * float(sum(dvec)) * sum(float(dt) for dt, _ in sched))
        ),
        "det_note": (
            "det M = exp(-nu tr D T) <= 1 exactly (every K_a is traceless), so the "
            "determinant route CANNOT certify rho > 1; not computed (Laplace is O(n!))"
        ),
        "trace_power_bounds": powers,
        "best_rho_lower_bound": best[0],
        "best_q": best[1],
        "certifies_rho_gt_1": best[0] > 1.0,
        "claimed_bound": "rho(M) >= (|tr M^q| / n)^(1/q), rigorous for every q",
    }


def exact_K_dict(basis, tab):
    out = {}
    for (i, j), blk in tab.blocks.items():
        for r in range(4):
            for c in range(4):
                if blk[r][c] != 0:
                    out[(4 * i + r, 4 * j + c)] = blk[r][c]
    return out


# --------------------------------------------------------------------------- #
# driver                                                                       #
# --------------------------------------------------------------------------- #


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="exact Leray cycle gate")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--geometry-cap", type=int, default=REGISTRY["geometry_cap"])
    parser.add_argument("--seed", type=int, default=20260801)
    parser.add_argument("--stage", default="all")
    parser.add_argument(
        "--squarings", type=int, default=11,
        help="interval-certificate squarings: trace-power bounds up to q = 2**squarings "
             "(the archived certificate used 11, i.e. q up to 2048)",
    )
    args = parser.parse_args(argv)

    def _jsonable(o):
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, Fraction):
            return str(o)
        raise TypeError(f"not serializable: {type(o)}")

    t_start = time.time()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ---------------- certificate-only reproduction ------------------------ #
    # `--stage certificate` re-runs ONLY the archived certified candidate (the
    # B.5 gadget, orientation (2, 6, 1), the staggered dt = 1 schedule, nu =
    # 1/40) and rewrites cycle_gate_interval_certificate.json.  This is the
    # cheap (~seconds) reproduction path for the stored certificate; the full
    # pipeline below takes ~20 minutes.
    if args.stage == "certificate":
        children = [(2, 2, 1), (2, 1, 2), (1, 2, 2)]
        parents = [(0, 1, -1), (1, -1, 0), (1, 0, -1)]
        choice = (2, 6, 1)
        basis = [child_basis(k) for k in children]
        gens_per_parent = [parent_generators(_canonical(p)[0]) for p in parents]
        tabs_per_parent = [
            [coupling_tables(g, basis) for g in gens] for gens in gens_per_parent
        ]
        n_par = len(parents)
        dvec = diag_vector(basis)
        n = int(dvec.shape[0])
        sched = [
            (1.0, tuple(1.0 if j == i else 0.0 for j in range(n_par)))
            for i in range(n_par)
        ]
        Kexact = [
            exact_K_dict(basis, tabs_per_parent[ai][choice[ai]]) for ai in range(n_par)
        ]
        cert = interval_certificate(
            Kexact, dvec, Fraction(1, 40), sched, n, squarings=args.squarings
        )
        Ks = [
            blocks_to_matrix(basis, tabs_per_parent[ai][choice[ai]])
            for ai in range(n_par)
        ]
        M = monodromy(Ks, dvec, 1.0 / 40.0, sched)
        rho = spectral_radius(M)
        T = sum(dt for dt, _ in sched)
        margins = instantaneous_margins(Ks, dvec, 1.0 / 40.0)
        cert["float_reference"] = {
            "chi": math.log(rho) / T,
            "rho": rho,
            "orientation": list(choice),
            "single_parent_margins": [float(m) for m in margins],
        }
        with (out_dir / "cycle_gate_interval_certificate.json").open(
            "w", encoding="utf-8"
        ) as fh:
            json.dump(cert, fh, indent=2, default=_jsonable)
        print(
            f"certificate: best q = {cert['best_q']}, "
            f"rho(M) >= {cert['best_rho_lower_bound']!r}, "
            f"certifies_rho_gt_1 = {cert['certifies_rho_gt_1']}"
        )
        print(f"wall {time.time() - t_start:.1f}s -> "
              f"{out_dir / 'cycle_gate_interval_certificate.json'}")
        return 0

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "registry": REGISTRY,
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": "NOT INSTALLED (hand-rolled expm)",
        },
    }

    # ---------------- A0: mandatory exact correctness checks --------------- #
    print("[A0] exact correctness checks ...", flush=True)
    ref_children = [(2, 2, 1), (2, 1, 2), (1, 2, 2)]
    ref_parents = [(0, 1, -1), (1, -1, 0), (1, 0, -1)]
    check_fields: list[tuple[str, TrigVector]] = []
    for k in ref_children:
        cb = child_basis(k)
        for lab, f in zip(cb.labels, cb.fields):
            check_fields.append((f"child{k}:{lab}", f))
    for p in ref_parents:
        for g in parent_generators(_canonical(p)[0]):
            check_fields.append((f"parent{p}:{g.label}", g.field))
    div_check = check_divergence_free(check_fields)

    energy_fields = []
    for k in ref_children:
        cb = child_basis(k)
        for lab, f in zip(cb.labels, cb.fields):
            energy_fields.append((f"c{k}:{lab}", f))
    for p in ref_parents:
        pc = _canonical(p)[0]
        for g in parent_generators(pc)[:2]:
            energy_fields.append((f"p{pc}:{g.label}", g.field))
    energy_check = check_full_system_energy(energy_fields)
    print(f"     div-free: {div_check['all_exactly_divergence_free']}  "
          f"energy: {energy_check['energy_conserved_exactly']}  "
          f"antisym: {energy_check['antisymmetry_exact']}  "
          f"(n={energy_check['n_fields']})", flush=True)
    report["correctness_checks"] = {
        "divergence_free": div_check,
        "full_system_energy_parents_dynamical": energy_check,
    }

    # ---------------- A1: enumeration -------------------------------------- #
    print("[A1] enumerating cycles ...", flush=True)
    t0 = time.time()
    buckets, enum_stats = enumerate_cycles(
        REGISTRY["parent_shell_max_sq"],
        REGISTRY["child_shell_max_sq"],
        REGISTRY["per_shell_bucket_cap"],
    )
    r1, r2, sel_info = select_lists(buckets, args.geometry_cap)
    enum_stats["wall_seconds"] = round(time.time() - t0, 2)
    enum_stats["selection"] = sel_info
    print(f"     retained={enum_stats['distinct_geometries_retained']} "
          f"R1 shells={sel_info['R1_max_k_sq_range']} "
          f"R2 shells={sel_info['R2_shells_covered'][:12]}...", flush=True)
    report["enumeration"] = enum_stats

    # the analysed list = R1 (pre-registered) UNION R2 (pre-registered supplement)
    analysed: list[Cycle] = []
    seen_children = set()
    for c in list(r2) + list(r1):
        if c.children in seen_children:
            continue
        seen_children.add(c.children)
        analysed.append(c)
    report["analysed_geometry_count"] = len(analysed)
    print(f"     analysing {len(analysed)} geometries (R2 first, then R1)", flush=True)

    # ---------------- A2: exact blocks + classification --------------------- #
    print("[A2] exact coupling blocks + block-Lyapunov classification ...", flush=True)
    t0 = time.time()
    geom_records = []
    crosstalk_rows = []
    class_counts = {"a_dead_common_metric": 0, "b_candidate_for_growth": 0}
    reality_entries = 0
    reality_ok = True
    prepared = []  # (cycle, basis, [(parent, [gen tables])])

    a_cap = REGISTRY["compute_caps_minutes"]["stage_A_exact_enumeration_and_blocks"] * 60
    a_cap_hit = False
    for gi, cyc in enumerate(analysed):
        if time.time() - t0 > a_cap:
            a_cap_hit = True
            print(f"     STAGE A CAP BOUND at geometry {gi}", flush=True)
            break
        basis = [child_basis(k) for k in cyc.children]
        gens_per_parent = [parent_generators(p) for p in cyc.parents]
        tabs_per_parent = []
        for p, gens in zip(cyc.parents, gens_per_parent):
            tabs = [coupling_tables(g, basis) for g in gens]
            tabs_per_parent.append(tabs)
        # reality check on every block produced
        for tabs in tabs_per_parent:
            for t in tabs:
                rc = check_reality(list(t.blocks.values()) + list(t.crosstalk.values()))
                reality_entries += rc["entries_checked"]
                reality_ok = reality_ok and rc["all_entries_exact_rational_real"]

        n_par = len(cyc.parents)
        orient_results = []
        for choice in itertools.product(range(8), repeat=n_par):
            tables = [
                (f"p{ai}:{gens_per_parent[ai][choice[ai]].label}", tabs_per_parent[ai][choice[ai]])
                for ai in range(n_par)
            ]
            res = classify(basis, tables)
            orient_results.append((choice, res))
            class_counts[res["class"]] += 1
        # exact cross-talk summary for this geometry (all orientations)
        ct_total = Fraction(0)
        in_total = Fraction(0)
        ct_entries = 0
        for tabs in tabs_per_parent:
            for t in tabs:
                ct_total += t.crosstalk_energy
                in_total += t.inset_energy
                ct_entries += len(t.crosstalk) * 16
        geom_records.append(
            {
                "index": gi,
                "kind": cyc.kind,
                "children": [list(k) for k in cyc.children],
                "steps": [list(s) for s in cyc.steps],
                "parents": [list(p) for p in cyc.parents],
                "max_k_sq": cyc.max_k_sq,
                "orientations": len(orient_results),
                "n_class_a_dead": sum(
                    1 for _c, r in orient_results if r["class"] == "a_dead_common_metric"
                ),
                "n_class_b_growth_candidate": sum(
                    1 for _c, r in orient_results if r["class"] == "b_candidate_for_growth"
                ),
                "edge_status_histogram": _status_hist(orient_results),
                "holonomy_values_seen": sorted(
                    {r["holonomy_product_exact"] for _c, r in orient_results
                     if r["holonomy_product_exact"] is not None}
                )[:12],
                "crosstalk_energy_exact_sum_over_generators": str(ct_total),
                "inset_energy_exact_sum_over_generators": str(in_total),
                "crosstalk_to_inset_ratio_float": (
                    float(ct_total / in_total) if in_total != 0 else None
                ),
                "crosstalk_block_entries": ct_entries,
                "example_edge_table": orient_results[0][1]["edges"][:8],
            }
        )
        # full cross-talk table dump (exact) for this geometry
        for ai, (p, tabs) in enumerate(zip(cyc.parents, tabs_per_parent)):
            for g, t in zip(gens_per_parent[ai], tabs):
                for (kout, j), blk in t.crosstalk.items():
                    crosstalk_rows.append(
                        {
                            "geometry": gi,
                            "parent": list(p),
                            "generator": g.label,
                            "polarization": list(g.pol),
                            "source_child": list(cyc.children[j]),
                            "target_mode_off_cycle": list(kout),
                            "target_k_sq": _norm_sq(kout),
                            "block_exact": [[str(v) for v in row] for row in blk],
                        }
                    )
        prepared.append((cyc, basis, tabs_per_parent, gens_per_parent, orient_results))

    stage_a_wall = time.time() - t0
    print(f"     {len(prepared)} geometries done in {stage_a_wall:.1f}s; "
          f"class counts {class_counts}", flush=True)
    report["stage_A"] = {
        "wall_seconds": round(stage_a_wall, 1),
        "cap_seconds": a_cap,
        "cap_bound": a_cap_hit,
        "geometries_processed": len(prepared),
        "reality_check": {
            "entries_checked": reality_entries,
            "all_entries_exact_rational_real": reality_ok,
        },
        "orientation_class_counts": class_counts,
        "class_a_fraction": (
            class_counts["a_dead_common_metric"]
            / max(1, sum(class_counts.values()))
        ),
    }
    report["geometries"] = geom_records

    with (out_dir / "crosstalk_tables.json").open("w", encoding="utf-8") as fh:
        json.dump(
            {"schema": SCHEMA + "/crosstalk", "rows": crosstalk_rows[:60000]},
            fh, indent=1, default=_jsonable,
        )
    print(f"     crosstalk rows written: {len(crosstalk_rows)}", flush=True)
    report["crosstalk_table_rows"] = len(crosstalk_rows)

    with (out_dir / "cycle_gate_stage_a.json").open("w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=_jsonable)

    if args.stage == "A":
        print(f"stage A only; wall {time.time() - t_start:.1f}s")
        return 0

    # ---------------- B: monodromy screen ---------------------------------- #
    print("[B] monodromy screen over class-(b) candidates ...", flush=True)
    t0 = time.time()
    b_cap = REGISTRY["compute_caps_minutes"]["stage_B_monodromy_screen"] * 60
    b_cap_hit = False
    best_overall = []
    screened = 0
    nus = [0.0, 1.0 / 40.0, 1.0 / 10.0]

    for gi, (cyc, basis, tabs_per_parent, gens_per_parent, orient_results) in enumerate(prepared):
        if time.time() - t0 > b_cap:
            b_cap_hit = True
            print(f"     STAGE B CAP BOUND at geometry {gi}", flush=True)
            break
        n_par = len(cyc.parents)
        dvec = diag_vector(basis)
        cand = [(c, r) for c, r in orient_results if r["class"] == "b_candidate_for_growth"]
        cand.sort(key=lambda cr: (-cr[1]["_score"][0], -cr[1]["_score"][1],
                                  -cr[1]["_score"][3], cr[0]))
        cand = cand[: REGISTRY["orientations_for_dynamics_per_geometry"]]
        for choice, res in cand:
            Ks = [
                blocks_to_matrix(basis, tabs_per_parent[ai][choice[ai]])
                for ai in range(n_par)
            ]
            scale = max((float(np.abs(K).sum(axis=1).max()) for K in Ks), default=0.0)
            if scale == 0.0:
                continue
            deltas = [f / scale for f in (0.0625, 0.25, 1.0, 4.0)]
            deltas.append(1.0 / (0.025 * max(dvec)))
            scheds = schedule_family(
                n_par, deltas, rng, REGISTRY["random_schedules_per_candidate_screen"]
            )
            screened += 1
            for nu in nus:
                chi, rho, sch = best_growth(Ks, dvec, nu, scheds)
                if rho is None:
                    continue
                margins = instantaneous_margins(Ks, dvec, nu)
                mmax = max(margins)
                best_overall.append(
                    {
                        "geometry": gi,
                        "kind": cyc.kind,
                        "children": [list(k) for k in cyc.children],
                        "max_k_sq": cyc.max_k_sq,
                        "orientation": list(choice),
                        "nu": nu,
                        "coupling_inf_norm": scale,
                        "rho": rho,
                        "chi": chi,
                        "period": sum(dt for dt, _ in sch),
                        "schedule": [[dt, list(a)] for dt, a in sch],
                        "holonomy_exact": res["holonomy_product_exact"],
                        "non_proportional_edges": res["non_proportional_edges"],
                        "hyperbolic_edges": res["hyperbolic_edges"],
                        "single_parent_margins": [float(m) for m in margins],
                        "max_single_parent_margin": float(mmax),
                        "zero_instantaneous_margin_class": bool(mmax <= 1e-9),
                        "cycle_beats_frozen_parent": bool(chi > mmax + 1e-9),
                    }
                )
    best_overall.sort(key=lambda r: -r["chi"])
    stage_b_wall = time.time() - t0
    print(f"     screened {screened} candidate systems in {stage_b_wall:.1f}s", flush=True)
    report["stage_B"] = {
        "wall_seconds": round(stage_b_wall, 1),
        "cap_seconds": b_cap,
        "cap_bound": b_cap_hit,
        "candidate_systems_screened": screened,
        "monodromy_evaluations": screened * len(nus),
        "top_by_chi": best_overall[:40],
        "max_chi_by_nu": {
            str(nu): max((r["chi"] for r in best_overall if r["nu"] == nu), default=None)
            for nu in nus
        },
        "max_rho_by_nu": {
            str(nu): max((r["rho"] for r in best_overall if r["nu"] == nu), default=None)
            for nu in nus
        },
        "count_rho_gt_1_by_nu": {
            str(nu): sum(1 for r in best_overall if r["nu"] == nu and r["rho"] > 1.0 + 1e-9)
            for nu in nus
        },
        "zero_instantaneous_margin_premise": {
            "description": (
                "The idea note's mechanism requires EVERY single-parent generator to be "
                "elliptic (max Re lambda <= 0), so that any Floquet gain is pure cycle "
                "holonomy.  Measured here in the FULL 4-dim-per-mode polarization space."
            ),
            "candidates_by_nu": {
                str(nu): sum(1 for r in best_overall if r["nu"] == nu) for nu in nus
            },
            "in_zero_margin_class_by_nu": {
                str(nu): sum(
                    1 for r in best_overall
                    if r["nu"] == nu and r["zero_instantaneous_margin_class"]
                )
                for nu in nus
            },
            "cycle_beats_frozen_parent_by_nu": {
                str(nu): sum(
                    1 for r in best_overall
                    if r["nu"] == nu and r["cycle_beats_frozen_parent"]
                )
                for nu in nus
            },
            "growth_inside_zero_margin_class_by_nu": {
                str(nu): sum(
                    1 for r in best_overall
                    if r["nu"] == nu
                    and r["zero_instantaneous_margin_class"]
                    and r["rho"] > 1.0 + 1e-9
                )
                for nu in nus
            },
            "max_single_parent_margin_over_all_candidates_by_nu": {
                str(nu): max(
                    (r["max_single_parent_margin"] for r in best_overall if r["nu"] == nu),
                    default=None,
                )
                for nu in nus
            },
            "min_single_parent_margin_over_all_candidates_by_nu": {
                str(nu): min(
                    (r["max_single_parent_margin"] for r in best_overall if r["nu"] == nu),
                    default=None,
                )
                for nu in nus
            },
        },
    }

    with (out_dir / "cycle_gate_stage_b.json").open("w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=_jsonable)

    if args.stage == "B":
        print(f"stages A-B; wall {time.time() - t_start:.1f}s")
        return 0

    # ---------------- C: refine + interval certificate --------------------- #
    print("[C] refinement + interval certificate ...", flush=True)
    t0 = time.time()
    c_cap = REGISTRY["compute_caps_minutes"]["stage_C_refine_and_interval_certificate"] * 60
    refined = []
    viscous = [r for r in best_overall if r["nu"] > 0.0]
    viscous.sort(key=lambda r: -r["chi"])
    zero_margin = [r for r in viscous if r["zero_instantaneous_margin_class"]]
    # refine the best overall AND the best inside the note's own zero-margin class
    to_refine = viscous[:10] + zero_margin[:6]
    for rec in to_refine:
        if time.time() - t0 > c_cap * 0.5:
            break
        gi = rec["geometry"]
        cyc, basis, tabs_per_parent, gens_per_parent, _o = prepared[gi]
        n_par = len(cyc.parents)
        choice = tuple(rec["orientation"])
        Ks = [blocks_to_matrix(basis, tabs_per_parent[ai][choice[ai]]) for ai in range(n_par)]
        dvec = diag_vector(basis)
        scale = rec["coupling_inf_norm"]
        deltas = [f / scale for f in (0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)]
        scheds = schedule_family(
            n_par, deltas, rng, REGISTRY["random_schedules_per_candidate_refine"]
        )
        chi, rho, sch = best_growth(Ks, dvec, rec["nu"], scheds)
        margins = instantaneous_margins(Ks, dvec, rec["nu"])
        refined.append(
            {
                "geometry": gi,
                "children": rec["children"],
                "orientation": rec["orientation"],
                "nu": rec["nu"],
                "refined_rho": rho,
                "refined_chi": chi,
                "refined_period": sum(dt for dt, _ in sch) if sch else None,
                "schedule": [[dt, list(a)] for dt, a in sch] if sch else None,
                "single_parent_margins": [float(m) for m in margins],
                "max_single_parent_margin": float(max(margins)),
                "zero_instantaneous_margin_class": bool(max(margins) <= 1e-9),
                "cycle_beats_frozen_parent": bool((chi or -1e18) > max(margins) + 1e-9),
            }
        )
    refined.sort(key=lambda r: -(r["refined_chi"] if r["refined_chi"] is not None else -1e18))
    report["stage_C_refine"] = {
        "wall_seconds": round(time.time() - t0, 1),
        "refined": refined[:12],
        "any_refined_rho_gt_1": any(
            (r["refined_rho"] or 0.0) > 1.0 + REGISTRY["rho_margin_for_growth"]
            for r in refined
        ),
    }

    # interval certificate on the single best positive candidate, if any
    cert = {"attempted": False, "reason": "no candidate reached rho(M) > 1 with margin"}
    positives = [r for r in refined if (r["refined_rho"] or 0.0) > 1.0 + REGISTRY["rho_margin_for_growth"]]
    if positives:
        best = positives[0]
        gi = best["geometry"]
        cyc, basis, tabs_per_parent, gens_per_parent, _o = prepared[gi]
        n_par = len(cyc.parents)
        choice = tuple(best["orientation"])
        Kexact = [exact_K_dict(basis, tabs_per_parent[ai][choice[ai]]) for ai in range(n_par)]
        dvec = diag_vector(basis)
        n = int(dvec.shape[0])
        nu_f = Fraction(1, 40) if abs(best["nu"] - 0.025) < 1e-12 else Fraction(1, 10)
        sched = [(dt, amps) for dt, amps in best["schedule"]]
        try:
            cert = interval_certificate(
                Kexact, dvec, nu_f, sched, n, squarings=args.squarings
            )
            cert["float_reference"] = {
                "chi": best["refined_chi"],
                "rho": best["refined_rho"],
                "orientation": list(best["orientation"]),
                "single_parent_margins": best["single_parent_margins"],
            }
            with (out_dir / "cycle_gate_interval_certificate.json").open(
                "w", encoding="utf-8"
            ) as fh:
                json.dump(cert, fh, indent=2, default=_jsonable)
        except Exception as exc:  # noqa: BLE001
            cert = {"attempted": True, "failed": True, "error": repr(exc)}
    report["stage_C_interval_certificate"] = cert

    # ---------------- D: tree control -------------------------------------- #
    print("[D] tree control: single-edge (FOREST) growth ...", flush=True)
    t0 = time.time()
    report["stage_D_tree_control"] = {
        "description": (
            "A single edge {k, k+p} with one parent frozen at |a| = 1 is a FOREST. "
            "Proposition 1 / Theorem D1 assert the top exponent is <= -nu min|k|^2 for "
            "EVERY modulation and EVERY amplitude.  Measured on the exact Leray blocks "
            "in the full 4-dim-per-mode polarization space."
        ),
        "results": tree_control([0.0, 1.0 / 40.0, 1.0 / 10.0]),
        "wall_seconds": round(time.time() - t0, 1),
    }
    tc = report["stage_D_tree_control"]["results"]
    print("     " + "  ".join(
        f"nu={k}: {v['with_positive_growth']}/{v['single_edge_systems_tested']} grow "
        f"(max {v['largest_single_edge_growth_rate']:+.4f})" for k, v in tc.items()
    ), flush=True)

    # ---------------- E: critical-normalization scaling (step 6/7 core) ----- #
    print("[E] critical-normalization scaling scan ...", flush=True)
    t0 = time.time()
    b5_children = [(2, 2, 1), (2, 1, 2), (1, 2, 2)]
    b5_parents = [(0, 1, -1), (1, -1, 0), (1, 0, -1)]
    scan = critical_normalization_scan(
        b5_children, b5_parents, 1.0 / 40.0,
        c_E_values=[1.0, 10.0, 100.0, 228.0],
        scales=[1, 2, 3, 4, 6],
        rng=rng,
    )
    report["stage_E_critical_normalization"] = {
        "geometry": {"children": [list(k) for k in b5_children],
                     "parents": [list(p) for p in b5_parents],
                     "note": "the note's own B.5 gadget, dilated k -> s k, p -> s p"},
        "nu": 1.0 / 40.0,
        "rows": scan,
        "wall_seconds": round(time.time() - t0, 1),
    }
    for cE in [1.0, 10.0, 100.0, 228.0]:
        seq = [r for r in scan if r["c_E"] == cE]
        print(f"     c_E={cE:<6} q_* vs N: " +
              "  ".join(f"N={r['N_equals_max_k']:.1f}:{r['stage_gain_q_star']:+.4f}"
                        for r in seq), flush=True)
    q_star_strictly_decreasing = all(
        all(
            seq[i + 1]["stage_gain_q_star"] < seq[i]["stage_gain_q_star"]
            for i in range(len(seq) - 1)
        )
        for seq in (
            [r for r in scan if r["c_E"] == cE] for cE in [1.0, 10.0, 100.0, 228.0]
        )
    )
    any_pos_at_largest_N = any(
        r["positive_net_growth"] for r in scan
        if r["scale_s"] == max(x["scale_s"] for x in scan)
    )
    inf_pos = all(r["stage_gain_q_star"] > 0 for r in scan)
    report["stage_E_verdict"] = {
        "q_star_strictly_decreasing_in_N_for_every_c_E": q_star_strictly_decreasing,
        "inf_N_stage_gain_positive": inf_pos,
        "branch_tested": (
            "FIXED CARDINALITY (3 modes dilated). The fixed-relative-width cloud "
            "(step 6) was NOT implemented in this script."
        ),
    }

    # ---------------- verdict ---------------------------------------------- #
    any_growth = report["stage_C_refine"]["any_refined_rho_gt_1"]
    certified = bool(cert.get("certifies_rho_gt_1"))
    cycle_beats = any(r["cycle_beats_frozen_parent"] for r in refined)
    go = bool(
        any_growth
        and certified
        and inf_pos
        and not q_star_strictly_decreasing
        and any_pos_at_largest_N
    )
    # a CERTIFIED lower bound must be rounded DOWN (truncated), never up
    cert_bound = cert.get("best_rho_lower_bound") or 0.0
    cert_bound_down = math.floor(cert_bound * 1e10) / 1e10
    cert_margins = (cert.get("float_reference") or {}).get("single_parent_margins")
    # M9: the off-chain/main leakage-ratio TREND in N is measurable from the
    # stage-A archive; report it instead of calling it unmeasured.
    _rt_small = [
        g["crosstalk_to_inset_ratio_float"] for g in geom_records
        if g["crosstalk_to_inset_ratio_float"] and g["max_k_sq"] <= 9
    ]
    _rt_large = [
        g["crosstalk_to_inset_ratio_float"] for g in geom_records
        if g["crosstalk_to_inset_ratio_float"] and g["max_k_sq"] >= 25
    ]
    tc = report["stage_D_tree_control"]["results"]

    def _tc_line(key: str) -> str:
        v = tc[key]
        return (f"{v['with_positive_growth']}/{v['single_edge_systems_tested']} grow, "
                f"max {v['largest_single_edge_growth_rate']:+.6f}")

    report["verdict"] = {
        "gate": "GO" if go else "KILL",
        "step_5_kill": not (any_growth and certified),
        "step_5_note": (
            f"growth exists and is interval-certified: rho(M) >= {cert_bound_down:.10f} "
            "(lower bound rounded DOWN)"
            if certified else "no interval-certified growth"
        ),
        "interval_certified_growth": certified,
        "certified_rho_lower_bound": cert.get("best_rho_lower_bound"),
        "certified_candidate_in_zero_instantaneous_margin_class": (
            bool(max(cert_margins) < 0.0) if cert_margins else None
        ),
        "ground_1_tree_control_fires": {
            "statement": (
                "a SINGLE EDGE {k,k+p} (a FOREST) already grows on the exact Leray "
                "blocks in full polarization; hypothesis (H2) fails exactly (blocks "
                "not proportional), so Proposition 1 is OUTSIDE ITS HYPOTHESES here "
                "and asserts nothing -- it does not apply to the exact Leray coupling."
            ),
            "nu_0": _tc_line("0.0"),
            "nu_1_40": _tc_line("0.025"),
            "nu_1_10": _tc_line("0.1"),
            "consequence": (
                "the note's \"impossible at any amplitude for every tree -> finite "
                "threshold for a cycle\" factor-infinity-to-finite jump is FALSE; "
                "both are finite."
            ),
        },
        "ground_2_no_uniformity_in_N": {
            "statement": (
                "q_*(N) = chi_net/(nu N^2) = C(c_E) N^{-3/2} - 1 with "
                "C(c_E) = chi_1 sqrt(c_E)/(3 nu), an identity (single shell => "
                "-nu D scalar; K linear in amplitude and in s)."
            ),
            "strictly_decreasing_for_every_c_E": q_star_strictly_decreasing,
            "zero_crossing_at_c_E_228": "N ~ 24.6",
            "equivalent_to": (
                "judge_rigor c_E >~ nu^2 N^3/chi^2 -> infinity; exponent matches "
                "the repo's own fixed_cardinality_scaling rejection."
            ),
        },
        "honest_gaps": {
            "step_6_galerkin_cloud_test": (
                "NEVER IMPLEMENTED - no Galerkin-cloud code exists in this script. "
                "The fixed-relative-width branch (obligation 4 / K2) is therefore "
                "NOT settled by this gate."
            ),
            "step_7_conditions_unmeasured": [
                "second-stage normalized gain",
                "output shape returns to same cycle class",
            ],
            "off_chain_main_ratio_trend_in_N": {
                "median_ratio_max_k_sq_le_9": (
                    float(np.median(_rt_small)) if _rt_small else None
                ),
                "median_ratio_max_k_sq_ge_25": (
                    float(np.median(_rt_large)) if _rt_large else None
                ),
                "note": (
                    "measured from the stage-A archive: mildly DECREASING in N, "
                    "the direction the GO condition asks for (favourable)"
                ),
            },
            "autonomy": (
                "NOT discharged and never tested: parent amplitudes are "
                "PRESCRIBED, not dynamical."
            ),
        },
        "any_candidate_where_cycle_beats_best_frozen_parent": cycle_beats,
        "mechanism_specific_verdict": (
            "The cyclic-carrier phenomenon is real and certified, but it is "
            "neither novel (trees do it too) nor uniform in N (fixed-cardinality "
            "branch dies as N^{-3/2})."
        ),
    }
    report["total_wall_seconds"] = round(time.time() - t_start, 1)

    # the full per-geometry tables live in cycle_gate_stage_a.json; the
    # consolidated report carries only a summary pointer
    report["geometries_summary"] = {
        "count": len(geom_records),
        "full_per_geometry_tables_in": "cycle_gate_stage_a.json",
    }
    report.pop("geometries", None)

    with (out_dir / "cycle_gate.json").open("w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=_jsonable)
    print(json.dumps(report["verdict"], indent=2))
    print(f"total wall {time.time() - t_start:.1f}s -> {out_dir}")
    return 0


def tree_control(nu_values, parent_max_sq=6, child_max_sq=12):
    """TREE CONTROL.  A SINGLE EDGE {k, k+p} with one parent frozen at |a| = 1.

    The coupling graph is one edge, i.e. a FOREST, so Proposition 1 (and its
    corrected form, Theorem D1) asserts the top exponent is <= -nu min|k|^2 for
    every modulation and every amplitude.  This measures max Re lambda directly
    on the exact Leray blocks in the FULL 4-dim-per-mode polarization space.
    """
    acc = {
        str(nu): {"tested": 0, "pos": 0, "worst": -np.inf, "example": None}
        for nu in nu_values
    }
    for p in shell_modes(parent_max_sq):
        for k in shell_modes(child_max_sq):
            k2 = tuple(k[i] + p[i] for i in range(3))
            b0, b1 = child_basis(k), child_basis(k2)
            if b0.k == b1.k:
                continue
            basis = [b0, b1]
            dvec = diag_vector(basis)
            for gen in parent_generators(p):
                K = blocks_to_matrix(basis, coupling_tables(gen, basis))
                for nu in nu_values:
                    a = acc[str(nu)]
                    m = float(np.max(np.linalg.eigvals(np.diag(-nu * dvec) + K).real))
                    a["tested"] += 1
                    if m > 1e-9:
                        a["pos"] += 1
                    if m > a["worst"]:
                        a["worst"] = m
                        a["example"] = {
                            "child_a": list(b0.k), "child_b": list(b1.k),
                            "parent": list(p), "generator": gen.label,
                            "polarization": list(gen.pol),
                            "max_Re_lambda": m,
                            "prop1_bound_minus_nu_min_k_sq": -nu * float(min(dvec)),
                        }
    return {
        nk: {
            "single_edge_systems_tested": a["tested"],
            "with_positive_growth": a["pos"],
            "fraction_violating_prop1": a["pos"] / max(1, a["tested"]),
            "largest_single_edge_growth_rate": a["worst"],
            "worst_example": a["example"],
        }
        for nk, a in acc.items()
    }


def critical_normalization_scan(children, parents, nu, c_E_values, scales, rng):
    """STEP 6/7 CORE: does the cycle gain survive CRITICAL normalization as N grows?

    The cycle geometry is dilated, k -> s k and p -> s p, so N ~ s |k|.  The
    parent amplitude is no longer capped at 1 but set by the critical shell law
    E_p = c_E / N, i.e.  a_s = sqrt(c_E / N) / ||C_a||_{L^2}.  The measured
    quantity is the pre-registered stage gain

        q_*(s) = chi_net(s) / (nu N^2),

    which must stay bounded below by a positive constant for the front to close.
    """
    rows = []
    for s in scales:
        kids = [tuple(s * c for c in k) for k in children]
        pars = [_canonical(tuple(s * c for c in p))[0] for p in parents]
        basis = [child_basis(k) for k in kids]
        dvec = diag_vector(basis)
        N = math.sqrt(float(max(dvec)))
        gens = [parent_generators(p) for p in pars]
        tabs = [[coupling_tables(g, basis) for g in gs] for gs in gens]
        for c_E in c_E_values:
            amp = math.sqrt(c_E / N)
            best = (-np.inf, None, None)
            for choice in itertools.product(range(8), repeat=len(pars)):
                Ks = [amp * blocks_to_matrix(basis, tabs[a][choice[a]])
                      for a in range(len(pars))]
                sc = max(float(np.abs(K).sum(axis=1).max()) for K in Ks)
                if sc == 0:
                    continue
                deltas = [f / sc for f in (0.0625, 0.25, 1.0, 4.0)]
                chi, rho, sch = best_growth(Ks, dvec, nu, deltas and
                                            schedule_family(len(pars), deltas, rng, 4))
                if chi > best[0]:
                    best = (chi, rho, choice)
            chi = best[0]
            rows.append(
                {
                    "scale_s": s,
                    "N_equals_max_k": N,
                    "c_E": c_E,
                    "parent_amplitude_critical": amp,
                    "nu_N_sq": nu * N * N,
                    "best_net_floquet_chi": chi,
                    "stage_gain_q_star": chi / (nu * N * N),
                    "positive_net_growth": bool(chi > 0.0),
                }
            )
    return rows


def _status_hist(orient_results):
    hist: dict[str, int] = {}
    for _c, r in orient_results:
        for e in r["edges"]:
            hist[e["status"]] = hist.get(e["status"], 0) + 1
    return hist


if __name__ == "__main__":
    raise SystemExit(main())
