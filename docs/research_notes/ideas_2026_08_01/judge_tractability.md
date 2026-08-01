# JUDGE REPORT — 14 candidates, NS certificate lab (2026-08-01)

## 1. Score matrix (0–5)

Axes: **Cl**=distance-to-Clay · **Nv**=novelty (lit + repo) · **Cs**=closability · **NG**=no-go evasion soundness · **Cp**=computability · **Iv**=interval feasibility · **Ln**=Lean-formalizability · **Fl**=falsifiability · **Cn**=genuine singularity/regularity connection

| # | candidate | Cl | Nv | Cs | NG | Cp | Iv | Ln | Fl | Cn | Σ |
|---|---|--|--|--|--|--|--|--|--|--|--|
| 1 | **monotone-dichotomy** (Λ) | 4 | 4 | 5 | 5 | 5 | 5 | 5 | 4 | 5 | **42** |
| 2 | **critical-quantity** (μ_N) | 4 | 4 | 4 | 5 | 4 | 4 | 5 | 5 | 5 | **40** |
| 3 | **parametric-resonance** | 3 | 5 | 3 | 5 | 5 | 4 | 5 | 5 | 3 | **38** |
| 4 | **helicity-ledger** | 3 | 3 | 4 | 4 | 4 | 4 | 5 | 5 | 3 | **35** |
| 5 | **moving-front-ode** | 3 | 4 | 3 | 4 | 4 | 3 | 4 | 5 | 4 | **34** |
| 6 | forcing-controller | 2 | 3 | 4 | 4 | 5 | 4 | 5 | 3 | 2 | 32 |
| 7 | pressure-feedback | 2 | 5 | 1 | 4 | 3 | 4 | 3 | 4 | 3 | 29 |
| 8 | lagrangian-degeneracy | 3 | 4 | 3 | 3 | 2 | 2 | 4 | 3 | 5 | 29 |
| 9 | two-scale-core | 2 | 3 | 3 | 4 | 3 | 3 | 4 | 4 | 3 | 29 |
| 10 | aniso-multiscale | 3 | 4 | 2 | 4 | 2 | 3 | 3 | 3 | 4 | 28 |
| 11 | log-periodic-multitype | 2 | 3 | 2 | 3 | 5 | 3 | 4 | 4 | 2 | 28 |
| 12 | space-fourier-mismatch | 3 | 4 | 2 | 4 | 2 | 2 | 2 | 4 | 4 | 27 |
| 13 | front-flow (seed) | 2 | 2 | 2 | 2 | 4 | 2 | 3 | 4 | 4 | 25 |
| 14 | multi-center (NSDL) | 2 | 3 | 2 | 3 | 2 | 3 | 2 | 3 | 3 | 23 |

One-line rationales (extremes only): **Λ** — every object already exists as exact rationals in `modal_front_actions`/`fourier_torus`, and (I.3) is scale-invariant at r=0 with a dominated-by-Serrin RHS. **μ_N** — derives the repo's own M_eff≳N³ floor a second time from the *flux* side, and (B.7) shows the frozen-parent no-go is the true evolution to O(τ). **parametric** — Prop 1 (tree no-go for *all* modulations) is the only genuinely new structural theorem in the portfolio and retro-explains every registered relay failure. **helicity** — D1 predicts an *exact zero* on the registered 9:1 cross-talk gate. **moving-front** — g*→2/5 contradicts the repo's γ=1/2 target with a measurable χ(N) exponent. Lowest scores are driven by needing new solvers (aniso 𝔓^δ, lagrangian tracer/SVD, space-Fourier Gevrey blob, NSDL 192³ dyad) rather than by wrong math.

## 2. FATAL flaws (exact equation / collision)

- **front-flow (13)** — **FATAL, two.** (a) S2.5: `e≡e_c` and `F=χξe^{3/2}` are inconsistent; capacity `F(1)ξ⁻²` lies *below* requirement `F(1)(1−β(ξ−1))` for all `1<ξ<φ`, so ξ_max=2 is unreachable and the "predicted profile" is not a solution. (b) T.2: a steady front flow at constant `a` **is literally a Leray profile** (`û=N⁻²Û(k/N)`); with the critical wake `U≍|x|⁻¹`, local energy and `∫∫|∇u|²≍∫λ⁻¹dt` are both finite ⟹ **Tsai1998 ⟹ u≡0 on ℝ³**. Survives only as torus DSS, gated by T.4.
- **pressure-feedback (7)** — **FATAL to the realization.** V3: along the *entire* locked family `min_θ p₂^mat = 1.036 > 1/2`; the material blob falls below √(ντ) and is viscously arrested. The exhibited family solves a closed 5-eq/6-unknown *model*, not NS.
- **multi-center (14)** — **FATAL to the mechanism as derived.** B.1 gives `ŝ·e₁ ≡ 0 ⟹ ḋ = 0` at line-vortex order: the dyad never collapses. The repair (generation j+1 = internal sub-dyad of the *same* tubes at the *same* Γ) is asserted; Kelvin conserves Γ per material circuit, not per sub-dyad.
- **log-periodic (11)** — **near-fatal.** `r = p−q` has `|r|²=2`, identical to the parent modulus, so 𝒜₀ and 𝒜₁ are *not* radially disjoint; λ₀=λ₁=√2 is an rms artifact of a bimodal set, and stage-0↔stage-1 same-radius cross-talk is outside LP-9's scope. Plus η≤1/14 ⟹ 50× loss in M_eff.
- **lagrangian (8) vs moving-front (5) vs helicity (4)** — **mutual contradiction:** V4 kills γ∈(1/2,1) via `𝒱≍νN²/λ_c→∞`; Lens 13 argues g∈(2/5,1/2); Lens 7 T3 pins γ=1/2 exactly. At most one is right; all three rest on order-of-magnitude estimates of the same viscous term.
- **aniso (10)** — non-fatal but load-bearing: `M_eff/N³ ≍ τ^{2a₁−a₂−a₃}=τ^{3/16}→0` *would* kill it; the three-point non-applicability argument (no critical wake / no empty child / no frozen parent) is unproven. Limit system is Prandtl-type, Sobolev-ill-posed (self-flagged).
- **forcing-controller (6)** — self-defeating: its own decoupling theorem proves `f` cannot reach the front, so (D) reduces to (B) and obligation 6 (flux positivity) is untouched; §E.10 concedes cross-talk unfixed.
- **helicity (4)** — T1 is Waleffe's 1992 helical-triad identity (repo-novel, not literature-novel); T2 needs (H3) sign-definite `H_k` for all k,t, unverifiable. Residual: `A₁+B₁`, `A₂+B₂` stay bichromatic and are untested.
- **Λ (1), μ_N (2), parametric (3)** — no fatal equation found. Λ's exposure is literature novelty (∫KD criteria are adjacent to known Serrin refinements) + O-3/O-4 bridges; parametric's is K2 (𝒞→−1 in the width-saturating limit) and that staggering costs 65% of the cycle gain (χ_par=0.0528 vs σ_cycle=0.149).

## 3. Top 5

**1. monotone-dichotomy · 2. critical-quantity · 3. parametric-resonance · 4. helicity-ledger · 5. moving-front-ode.**
1–2 are the only PROOF CANDIDATEs whose central objects are *already* exact-rational repo computations; 3 is the best new mechanism-side theorem; 4 has the sharpest one-hour test; 5 has the sharpest scaling prediction.

## 4. Answers to the tractability questions

- **Advanceable THIS SESSION:** Λ (1) — `modal_front_actions._moments` + `fourier_torus` (which already ships `_FractionOps` *and* `_IntervalOps`) supply H_r, T_r exactly; only G_r, V_r, Cov, Γ^CS, Γ^SC are new (~120 lines). All five test-field builders exist. Runner-up: parametric P1 (pure `Fraction`, `fourier_torus.{advection,leray}` + `exact_carrier_search` enumeration, ~20 min).
- **Shortest path to an exact-rational / interval-verified statement:** Λ's **Lemma K** (`K ≤ 4πN+C`, a finite lattice sum) and identity **(I.4)** (exact equality with two nonnegativity assertions). Both are finite sums + Cauchy–Schwarz + square completion — Lean-able in the `MesoscopicDuhamelNoGo.lean` style *today*, with interval enclosure essentially free via the existing ops backend.
- **Sharpest falsifiable prediction:** helicity **D1** — the 2-colouring `s(A₁)=s(B₂)=+`, `s(A₂)=s(B₁)=−` at `|k|²=2` makes both diagonal cross channels **exactly 0** (coefficient `s|ℓ|−s'|m|`) while intended channels take the maximal `2√2`. Binary, exact, and it flips a registered 16/16 negative result. Runner-up: parametric P2's null control (`𝒞=−1 ⟹ ρ(M)=1±10⁻¹²` vs `𝒞=−1239/128 ⟹ ρ=1.1100`).

## 5. Best single ≤1h pilot

**Run Lens 11 §F: `src/ns_certificate_lab/spectral_front_monotone.py` — exact-rational front-gap ledger.**

`front_gap_identity(field, r, ν) → Fraction` record of `{H_r,H_{r+1},H_{r+2}, T_r,T_{r+1}, G_r, V_r=Σp_r(x−μ)², Cov, Γ^CS_r, Γ^SC_r, G_r/(4νH_{r+1})}`, asserting **(I.4) as an exact equality** and `Γ^CS,Γ^SC ≥ 0`; plus `front_wavenumber → K=G₀/2H₁²`, `lemma_k_margin → 4πN−K`, and `renyi_half_perplexity → M_eff`.

Inputs (all already exact-rational, zero floats, zero grids): `exact_leray_relay.build_exact_relay_triad()`; both strict orientations from `exact_carrier_search.search_exact_carrier_gadget()`; `carrier_two_stage_galerkin.build_partial_carrier_parent()`; `mesoscopic_cloud_scaling.build_sparse_parent()` at η∈{0.10,0.15,0.20}, N∈{16,32,48,64}. Compute V_r as `Σp_r(x−μ)²`, never as `H_{r+2}/H_r−μ²` (TM-09).

Why this one: it is the only pilot with ~0.9 completion probability that needs *no* new numerics, *no* threat-model exposure, and produces an exact-arithmetic verdict either way — and because `renyi_half_perplexity = M_eff = e^{H_{1/2}}` is precisely Lens 10's `μ_N N³`, **the same hour instruments the #1 and #2 candidates on the same field set**.

Pre-registered kills: (K1) any `Γ^CS<0`, `Γ^SC<0`, or (I.4) failing exact equality ⟹ Lens 11 derivation is wrong, abandon. (K2) `lemma_k_margin<0` on any band-limited field ⟹ re-derive the lattice constant. (K3) `𝔡 = (Γ^CS₀+Γ^SC₀)/(G₀/2νH₁) → 0` along the N-sweep ⟹ Λ demoted to a consistency diagnostic. (K4, Lens 10) `μ_N → 0` across η∈(0,1/3) as N grows ⟹ **the only structurally uncondemned mesoscopic family fails the delocalization floor and the whole cloud lane dies by a new, independent argument** — the single most valuable negative outcome available in one hour.