# VERIFICATION SPRINT v1 — WORKSTREAM D
## Parametric resonance, Proposition 1 (common diagonal Lyapunov) — audit

**Audited document:** `docs/research_notes/ideas_2026_08_01/idea_parametric_resonance.md`
(Lens 6), Propositions 1 and 2, sections B.2–B.6.
**Judge criticisms consulted:** `judge_rigor.md` (parametric = "decisive gap, not
fatal"; `c_E ≳ ν²N³/χ_par²→∞`), `judge_novelty.md` ("Prop 1 … is a new exact
no-go"; "staggered pump is imposed, not autonomous"), `judge_tractability.md`
("no fatal equation found … exposure is K2 and that staggering costs 65 % of the
cycle gain").
**Script:** `experiments/run_parametric_lyapunov_audit.py`
**Data:** `outputs/verification_sprint_v1/parametric/parametric_common_lyapunov_audit.json`
**No new ideas were generated.** This note verifies, corrects, or kills existing claims only.

---

## 0. Verdict summary

| item | status |
|---|---|
| Prop 1 as a statement about **forest topology alone** | **CORRECTED** — topology is neither necessary nor sufficient; one *unstated* hypothesis (H1, reciprocal modulation) is load-bearing, and forest is only a sufficient proxy for the real one (H3, cycle balance) |
| Prop 1's proof sketch (`θ_i/θ_j = −β_ij/β_ji`) | **CORRIGENDUM** — indices transposed; the correct metric is `θ_i/θ_j = −β_ji/β_ij` (existence unaffected) |
| Prop 1 under the corrected hypotheses (H1)+(H2)+(H3) | **PROVEN** (finite linear algebra, all three modulation classes (a)/(b)/(c) at once) |
| Prop 2 (`C ≠ −1 ⇒ Re λ > 0`) | **CONFIRMED**, and sharpened (factor-2 rate asymmetry in `sign det K`) |
| `C = −1` null control (`ρ = 1 ± 10⁻¹²`) | **CONFIRMED** — and upgraded from numerics to a *theorem* |
| B.5 exact invariants `r_e`, `C = −1239/128` | **CONFIRMED exactly** (gauge-consistency verified in `Fraction`) |
| B.6 `ω_c = 1.0796`, `det K = −0.58192` | **ERRATUM** — both carry a spurious factor `√2` |
| B.6 "the growth is pure holonomy" | **DEFLATED** — ≈ 92 % of `χ_par` is the arithmetic mean of the generators |
| B.6 `ρ(M) > 1` as evidence | **WEAKENED** — `ρ(M) ≥ 1` is structural (traceless generators), so only the null control carries information |

---

## 1. The audited claim, formalized

The child band is the linear non-autonomous system

    (1)   ẋ = A(t) x ,   A(t) = −ν D + K(t) ,   D = diag(|k₁|², …, |k_n|²),

with `K(t)` real, zero-diagonal, supported on a coupling graph `G` on the child
modes `{k₁,…,k_n}`, and (the note's model, §B.2/§G.1)

    (2)   K_ij(t) = β_ij X_{P(ij)}(t) ,   β_ij β_ji < 0  (cascade edges).

`x` is written in **energy coordinates** so that `E_i = x_i²/2`; the note's
triad reduction (§B.1) fixes `(β_B,β_C,β_D) = √(2/3)·s·(−½,−½,1)`.

> **Claim under audit (Prop 1 as written).** If `G` is a **forest**, then a
> constant positive diagonal `Θ` makes `K(t)` `Θ`-antisymmetric for all `t`,
> whence `‖x(t)‖_Θ ≤ e^{−ν min_i|k_i|² t}‖x(0)‖_Θ`, i.e. the top Floquet
> exponent obeys `χ_max ≤ −ν min_i |k_i|²` for **every** modulation — any
> amplitude, waveform, period.

Equivalent Lyapunov formulation (task item 2): find a constant `W = diag(w_i) ≻ 0`
and `δ > 0` with

    (3)   W A(t) + A(t)ᵀ W ⪯ −2 δ W    for all t.

Since `A = −νD + K`, (3) splits as

    (4)   W A + Aᵀ W = −2ν W D + (W K + Kᵀ W).

`W D` is diagonal positive, so `−2νWD ⪯ −2ν(min_i|k_i|²) W` always. Therefore

    (5)   **(3) holds with δ = ν·min_i|k_i|²  ⇐  W K(t) + K(t)ᵀ W ⪯ 0 for all t.**

The whole question is the *sign-semidefiniteness of the coupling part in a
constant diagonal metric*. Topology enters only through the solvability of the
edge conditions.

---

## 2. When does a common positive diagonal `W` exist?

### 2.1 Edge-wise reduction

`W K + Kᵀ W` is symmetric with zero diagonal (`K` has zero diagonal) and
`(i,j)` entry

    (6)   (W K + Kᵀ W)_ij = w_i K_ij(t) + w_j K_ji(t).

`W K + Kᵀ W` has **zero diagonal** (`K_ii = 0`), and a symmetric matrix with zero
diagonal is negative semidefinite **iff it is zero**: `M ⪯ 0` with `M_ii = 0`
forces every `2×2` principal minor `[[0, M_ij],[M_ij, 0]]` (eigenvalues `±|M_ij|`)
to be `⪯ 0`, i.e. `M_ij = 0`. Hence:

> **Lemma D0.** Let `W ≻ 0` be diagonal and `δ = ν min_i|k_i|²`.
> (i) If all children lie on one shell (`D = |k|² I`), condition (3) is
> *equivalent* to `W K(t) + K(t)ᵀ W = 0` — exact `W`-antisymmetry.
> (ii) In general (3) leaves slack `2ν W(D − min_i|k_i|² I) ⪰ 0` on the
> non-minimal shells, but that slack is **fixed** while `K` is not: since
> Prop 1 quantifies over *all* modulation amplitudes and
> `W K(cm) + K(cm)ᵀ W = c (W K(m) + K(m)ᵀ W)`, dividing by `c` and letting
> `c → ∞` forces `W K(m) + K(m)ᵀ W ⪯ 0`, hence `= 0`.
>
> **There is no weaker "merely dissipative" diagonal certificate to exploit.**
> Any constant-diagonal proof of Prop 1 must go through exact antisymmetrization.

So the problem reduces, edge by edge, to the `2×2` block condition

    (7)   w_i K_ij(t) + w_j K_ji(t) = 0     for every edge {i,j} and every t,

i.e.

    (8)   w_i / w_j = − K_ji(t) / K_ij(t) .

The right-hand side of (8) must be **(a) positive** and **(b) independent of `t`**.
These are exactly the two hypotheses the note leaves implicit.

### 2.2 The three hypotheses that are actually used

> **(H1) Reciprocal (shared-parent) modulation.** For every edge `e = {i,j}`
> there is a *single* scalar `m_e(t)` with
> `K_ij(t) = β_ij m_e(t)` **and** `K_ji(t) = β_ji m_e(t)`.
> Then `K_ji/K_ij ≡ β_ji/β_ij` is time-independent.

> **(H2) Elliptic edges.** `β_ij β_ji < 0` for every edge, so
> `ρ_ij := −β_ji/β_ij > 0`.  *(Stated in the note as "cascade edges".)*

> **(H3) Cycle balance.** Around every cycle `γ = (i₁ i₂ … i_r i₁)` of `G`,
> `∏_{s} ρ_{i_s i_{s+1}} = 1`. A forest has no cycles, so (H3) is vacuous —
> **this, and only this, is what "forest" buys.**  *(Implicit in the note, via
> the forest hypothesis.)*

Of the three, only **(H1) is absent from the note in any form**, and it is the
one whose failure produces genuine unbounded growth on a tree (§4.2).

Note that (H1) is *not* the statement that the modulation is slow, periodic, or
of fixed shape. It is a **structural reciprocity**: both directions of one edge
are driven by the same parent amplitude. In the Fourier NS model this is
automatic — the edge `k_i → k_j` needs the parent `p = k_j − k_i` and the edge
`k_j → k_i` needs `−p`, and for a real field `û_{−p} = conj(û_p)`, one complex
degree of freedom. (In a real-amplitude reduction with frozen parent
polarization and phase, `m_e(t) = X_p(t)`; see §7 for the residual complex-phase
gap.) **But it is a hypothesis, and it is the one that can fail** (§4).

### 2.3 Theorem D1 (corrected Proposition 1) and its proof

> **Theorem D1.** Assume (H1), (H2), (H3) for the system (1)–(2). Then there is
> a constant `Θ = diag(θ_i) ≻ 0`, depending only on the constants `β_ij` and not
> on the modulations `m_e(·)`, such that
>
>     Θ K(t) + K(t)ᵀ Θ = 0    for all t,
>
> and consequently, for every locally integrable modulation `m_e : ℝ → ℝ` — any
> amplitude, any waveform, any period, any sign, independent across edges —
>
>     (9)   d/dt ‖x‖²_Θ = −2ν xᵀ Θ D x ,
>     (10)  e^{−2ν max_i|k_i|² (t−t₀)} ≤ ‖x(t)‖²_Θ / ‖x(t₀)‖²_Θ ≤ e^{−2ν min_i|k_i|² (t−t₀)} .
>
> Hence every Lyapunov/Floquet exponent of (1) lies in
> `[ −ν max_i|k_i|² , −ν min_i|k_i|² ]`, and in the Euclidean norm
>
>     (11)  ‖x(t)‖₂ ≤ √(κ(Θ)) · e^{−ν min_i|k_i|²(t−t₀)} ‖x(t₀)‖₂ ,
>       κ(Θ) = max_i θ_i / min_i θ_i .

*Proof.* Build `Θ` by propagation. Pick a spanning forest of `G` and a root in
each component; set `θ_root = 1`; for a tree edge `{i,j}` with `j` the child, set
`θ_j = θ_i / ρ_ij`, so `θ_i/θ_j = ρ_ij = −β_ji/β_ij > 0` by (H2). Every non-tree
edge closes a cycle, and (H3) says the product of the required ratios around
that cycle is 1, i.e. the already-assigned `θ_i, θ_j` satisfy (8) automatically.
So (8) holds on all of `G`, and by (H1) the required ratio is time-independent,
so a single constant `Θ` works for all `t`. Then (7) holds edgewise for every `t`
and every modulation, so `ΘK(t)+K(t)ᵀΘ = 0`. Differentiating `‖x‖²_Θ = xᵀΘx`
along (1),

    d/dt xᵀΘx = xᵀ(ΘA + AᵀΘ)x = xᵀ(ΘK+KᵀΘ)x − 2ν xᵀΘDx = −2ν Σ_i θ_i |k_i|² x_i²,

which is (9). Bounding `min_i|k_i|² ≤ |k_i|² ≤ max_i|k_i|²` gives (10) by
Grönwall; (11) follows from `min θ · ‖x‖₂² ≤ ‖x‖²_Θ ≤ max θ · ‖x‖₂²`. Floquet
exponents of a periodic modulation are Lyapunov exponents of (1), so they are
trapped by (10). ∎

**Both the upper and the lower bound in (10) are new relative to the note**, and
the lower bound matters: it says the child band cannot be over-damped either.

### 2.4 The three modulation classes are all closed by Theorem D1

| class | definition | status |
|---|---|---|
| (a) fixed edge ratios | `K(t) = m(t) K₀`, one scalar `m` | **PROVEN** (special case of (H1) with `m_e ≡ m`) |
| (b) common scalar modulation | one waveform, e.g. `m(t) = m₀ + m₁cos Ωt`, shared by all edges | **PROVEN** (same) |
| (c) fully independent edge modulations | `m_e(t)` arbitrary and independent per edge, including exact zeros (staggered/duty-cycled pumps) | **PROVEN** |

The reason is structural and worth stating plainly: **`Θ` is built from the
constants `β_ij` only.** The modulations never enter it. This is why (c) — the
class the note's own §B.6 staggered pump lives in — costs nothing extra, and it
is why the "no modulation, any amplitude/waveform/period" phrasing of Prop 1 is
*correct as stated for these classes*. The note's escape in §B.4 is therefore
**not** "modulation beats the metric"; it is "an unbalanced cycle has no metric".

### 2.5 Necessity: (H2) and (H3) cannot be dropped

- **(H2) is necessary, and its failure is not a large-amplitude effect.** If
  `β_ij β_ji > 0` on some edge, (8) demands `w_i/w_j < 0`, impossible for
  `W ≻ 0`. The *claim* fails too, and immediately: freezing `m_e ≡ 1`, the
  isolated `2×2` block has eigenvalues
  `−ν k̄² ± √( ν²Δ² + β_ijβ_ji )`, `k̄² = (|k_i|²+|k_j|²)/2`,
  `Δ = (|k_j|²−|k_i|²)/2`. Since `√(ν²Δ² + β_ijβ_ji) > ν|Δ|` for **any** nonzero
  coupling, `χ_max > −νk̄² + ν|Δ| = −ν min(|k_i|²,|k_j|²)`: the Prop-1 bound is
  violated at *every* nonzero amplitude, on a **single edge**, i.e. a forest.
  Confirmed numerically (§4.2, control H2-broken: `+589.1`, and 94–99.9 % of
  samples positive).
- **(H3) is necessary and forest is only sufficient.** Around a 3-cycle,
  `ρ_a ρ_b ρ_c = (−β_ba/β_ab)(−β_cb/β_bc)(−β_ac/β_ca) = −C` in the note's own
  notation, so (H3) ⟺ `C = −1` ⟺ `det K = 0`. **A balanced cycle is therefore
  just as dead as a tree.** The note's §B.4 sentence "Proposition 1 fails exactly
  when the coupling graph carries a cycle" is an overstatement: it fails exactly
  when the graph carries an **unbalanced** cycle. §4.2 verifies the balanced
  3-cycle at `+4.6 × 10⁻¹³` over 40 000 fully independent modulations, against
  `+561.8` for the unbalanced one on the same grid.

### 2.6 Corrigendum to the note's proof sketch

The note writes `θ_i/θ_j = −β_ij/β_ji`. From (8) the correct relation is the
**reciprocal**, `θ_i/θ_j = −β_ji/β_ij`. Both are positive under (H2), so the
*existence* argument is unaffected, but the metric produced by the note's formula
is the transposed one. Sanity check on the repo triad (§B.2): the `(X_C, X_D)`
block has `K_{CD} = β_C X_B = −(s/2)√(2/3) X_B` and `K_{DC} = β_D X_B = s√(2/3) X_B`,
so the correct `θ_C/θ_D = −K_{DC}/K_{CD} = 2`, `κ(Θ) = 2`, and
`E(t) ≤ 2 e^{−2λ_C t} E(0)` — reproducing the note's own "×2 overshoot" line.
(`κ(Θ)` is invariant under transposing the ratio in this 2-node case, which is
why the error is invisible in the triad.)

### 2.7 Two corollaries the note does not state

**(C1) Trace identity / floor.** `tr A(t) = −ν tr D` for every `t` because `K`
is zero-diagonal. Hence `det M = exp(−ν tr D · T)` for the monodromy over a
period `T`, so `Σ_i χ_i = −ν Σ_i |k_i|²` and

    (12)   χ_max ≥ −ν · mean_i |k_i|²    for ANY zero-diagonal K, forest or not.

Combined with Theorem D1 this **pins** the answer on a single shell: if all
`|k_i|² = |k|²` and (H1)–(H3) hold, then `χ_max = −ν|k|²` *exactly*, not merely
`≤`. Verified to `1.6 × 10⁻¹¹` relative (§6). This also means the
`min` in Prop 1 is sharp only when the children are on one shell; off a single
shell the true bound is somewhere in `[−ν mean|k|², −ν min|k|²]` and Prop 1 is
not tight.

**(C2) Transient overshoot is not `×2` in general — it is exponential in tree
depth.** From (11) the Euclidean overshoot is `√κ(Θ)`, and along a path
`i₁ … i_r`, `κ(Θ) = ∏_e max(ρ_e, ρ_e^{−1})` (worst case). On a 4-node path with
per-edge ratios drawn from `[1/20, 20]` the audit measures
`κ(Θ) = 5.950 × 10³`, i.e. an Euclidean overshoot of **≈ 77.1×**, with the Floquet
exponent still exactly `−ν min|k|²`. This is a direct, quantified instance of the
note's own threat **TM-P1** ("transient non-normality mistaken for Floquet
growth") and it is *much* worse than the "≤ ×2" the note quotes from the triad.
Any numerical relay experiment on a deep tree can show two orders of magnitude of
energy growth with zero Floquet gain.

---

## 3. Search protocol

All searches evaluate the **Floquet margin**

    (13)   margin := log ρ(M) / T + ν·min_i |k_i|² ,    M = ∏_{j=S..1} exp(A_j t_j),

with `T = Σ_j t_j`. `margin > 0` is a counterexample to the Prop-1 bound;
`margin = 0` means the bound is attained.

*Numerics (FLOAT, float64).* `scipy` is absent. `expm` is a hand-rolled batched
scaling-and-squaring Taylor evaluation (22 terms at `‖A/2^s‖_∞ ≤ 1/4`) that
renormalises the partial product to unit `∞`-norm after the Taylor step and
after every squaring, carrying `log` of the discarded factor additively; the
monodromy product is renormalised the same way. This makes the pipeline immune to
both overflow (the growth controls reach `log ρ ≈ 10³`) and underflow. Validated
against an eigendecomposition reference to `≤ 3 × 10⁻¹³` relative on `3×3`
matrices with `‖A‖` up to `10³`.

*Arithmetic split (repo convention).*
**EXACT `Fraction`:** the cycle invariant `C`, edge ratios `r_e`, `det K`,
`c₁ = Σ ω_e²`, the balanced-cycle construction and its `Θ`, and the
gauge-consistency check on the B.5 gadget.
**FLOAT:** every matrix exponential, monodromy, spectral radius, Floquet
exponent, cubic root, and `δ`-scan. No float quantity supports a PROVEN label.

*Tolerance.* `MARGIN_TOL = 10⁻⁷` absolute. The measured float64 round-off floor
on the exactly-balanced null controls is `4.6 × 10⁻¹³`; on the randomized forest
families with couplings up to `10³` and `T` up to `24` it is `≈ 3 × 10⁻¹⁰`
(relative to `ν min|k|² = 0.45`, i.e. `7 × 10⁻¹⁰`). Nothing in the audit sits
between `10⁻⁹` and `10⁻⁷`, so the tolerance is not load-bearing.

### 3.1 Domains

**Randomized families** (72 forest families × 4000–60000 samples each):

| axis | domain |
|---|---|
| graphs | `edge_2node`, `path_3node`, `star_4node`, `path_4node` |
| `\|β_fwd\|` | log-uniform on `[0.02, 50]`, random sign |
| edge ratio `\|β_bwd/β_fwd\|` | log-uniform on `[0.05, 20]` (so `\|β_bwd\|` up to `10³`) |
| modulation values | class (a)/(c): uniform `[−3, 3]`, with 25 % exact zeros in (c); class (b): `m₀ + m₁ sin(2πt/T + φ)`, `m₁ ∈ [0.2,3]`, `m₀ ∈ [−2,2]`, `φ ∈ [0,2π)` |
| switching intervals per period `S` | 2, 3, 4 |
| interval lengths | log-uniform on `[0.02, 6]` (so `T ∈ [0.04, 24]`) |
| shells | single shell `\|k_i\|² = 9` (`ν = 0.05`) and mixed `\|k_i\|² ∈ {2,6,9,14}` (`ν = 0.5`) |

**Exhaustive cartesian grids** (no RNG; every point of the full product
`β-magnitudes^{2E} × modulation-values^{S·E} × interval-lengths^{S}` evaluated):

| grid | `β` magnitudes | modulation values | interval lengths | `S` | points |
|---|---|---|---|---|---|
| `edge_2node`, class (a) | `{0.25,1,4,16}` | `{−2,−0.5,0,1,3}` | `{0.1,0.7,3}` | 2 / 3 / 4 | 3 600 / 54 000 / 810 000 |
| `path_3node`, class (c) | `{0.25,1,4}` | `{−2,0,1,3}` | `{0.3,2}` | 2 / 3 | 82 944 / 2 654 208 |
| `cycle_3node`, class (c) *(control)* | `{0.5,2}` | `{0,1,3}` | `{0.3,2}` | 3 | 10 077 696 |

Sign of `β_fwd` ranges over both values through the `β`-magnitude product, and
`β_bwd = −sign(β_fwd)·(magnitude)` under (H2). Forest grid total: **3 604 752**
exhaustive evaluations; grand total including the cycle control: 13 682 448.
Wall clock for the whole audit (all randomized families + all grids + the
`δ`-scans): **348 s** on one CPU core.

**Controls with a hypothesis deliberately broken** (these *must* produce
positive margins or the harness is not sensitive):

- **class (d) = (H1)-broken:** the two directions of each edge get *independent*
  modulations `m⁺_e(t) ≠ m⁻_e(t)`. Graphs `edge_2node`, `path_3node`.
- **(H2)-broken:** `β_ij β_ji > 0` (hyperbolic edge), forest graphs.
- **(H3)-broken:** unbalanced 3-cycle, `C ≠ −1`.
- **(H3)-satisfied cycle (null control):** exact `C = −1`, fully independent
  per-edge modulation.

---

## 4. Results — Proposition 1

*(numbers from `parametric_common_lyapunov_audit.json`; `margin` per (13))*

### 4.1 Forest families under (H1)+(H2) — no counterexample

| family group | evaluations | max margin | verdict |
|---|---|---|---|
| 72 randomized forest families, classes (a)/(b)/(c), `S ∈ {2,3,4}`, single + mixed shell | 4 320 000 | `+6.82 × 10⁻¹⁰` | no counterexample |
| exhaustive grids, `edge_2node` + `path_3node` | 3 604 752 | `+1.64 × 10⁻¹²` | no counterexample |

Per graph and modulation class (randomized, 360 000 samples each):

| graph | class (a) | class (b) | class (c) |
|---|---|---|---|
| `edge_2node` | `4.10e−10` | `2.91e−10` | `2.78e−10` |
| `path_3node` | `3.67e−10` | `4.38e−10` | `3.36e−10` |
| `star_4node` | `4.34e−10` | `4.64e−10` | `3.06e−10` |
| `path_4node` | `3.44e−10` | `6.82e−10` | `3.32e−10` |

Exhaustive grids:

| grid | `S` | evaluations | max margin |
|---|---|---|---|
| `edge_2node` (a) | 2 / 3 / 4 | 3 600 / 54 000 / 810 000 | `1.639e−12` / `1.640e−12` / `1.641e−12` |
| `path_3node` (c) | 2 / 3 | 82 944 / 2 654 208 | `1.895e−13` / `1.895e−13` |

Every value sits at the float64 round-off floor, and on single-shell families the
margin is `0` to round-off *from both sides* — exactly as corollary (C1)
predicts (`χ_max = −ν|k|²` is pinned, not merely bounded). Note the exhaustive
grids are *cleaner* than the randomized families (`10⁻¹²`–`10⁻¹³` vs `10⁻¹⁰`)
because they do not reach the extreme edge-ratio spread (`κ(Θ)` up to `6 × 10³`)
that costs conditioning in the randomized draws.

Direct certificate check (`theta_certificate`, 4-node path, 20 000 random
`β`-draws and modulation values): the constructed `Θ` satisfies
`‖ΘK + KᵀΘ‖ / ‖ΘK‖ ≤ 4.438 × 10⁻¹⁶`, and the worst instantaneous
`Θ`-metric rate equals `−2ν min|k|²` exactly. The certificate is not merely
consistent with the data — it is the mechanism.

Trace identity (`trace_identity`, 3-node path, mixed shell, `ν = 0.5`):
`log|det M| = −ν tr D · T` to `1.6 × 10⁻¹¹` relative; the observed
`min χ_max = −2.8317` against the theoretical floor `−ν mean|k|² = −2.8333` and
the Prop-1 ceiling `−ν min|k|² = −1.0`. Both bounds respected.

### 4.2 Controls — the hypotheses are load-bearing

| control | hypothesis broken | max margin | frac. of samples with margin > 0 | interpretation |
|---|---|---|---|---|
| `class_d`, `edge_2node` & `path_3node`, `S ∈ {2,3,4}` | (H1) | **`+534.0`** | 0.60 – 0.998 | **(H1) is load-bearing.** Per-direction modulation of a *single* edge already produces unbounded Floquet growth on a forest. |
| `β_ij β_ji > 0`, `edge_2node` & `path_3node` | (H2) | **`+589.1`** | 0.94 – 0.999 | (H2) is load-bearing (and elementary: a hyperbolic edge). |
| unbalanced 3-cycle, classes (a)/(c), `S ∈ {3,4}` | (H3) | **`+561.8`** | 0.95 – 0.9997 | (H3) is load-bearing; this is the note's own escape channel. |
| exhaustive `cycle_3node` grid, 10 077 696 pts | (H3) | **`+4.50`** | — | same, on a deterministic grid |
| **balanced 3-cycle, `C = −1` exactly**, classes (a)/(c), `S ∈ {3,4}` | none | `+4.41 × 10⁻⁹` | **0.0000** | **A cycle obeying (H3) is as dead as a tree.** Forest is sufficient, not necessary. |
| **balanced 3-cycle, exact-rational `Θ` null control**, 40 000 independent-modulation samples | none | `+4.61 × 10⁻¹³` | **0.0000** | passes the note's own `ρ = 1 ± 10⁻¹²` bar; `Θ = (1,2,6)` exact, cycle-closure residual `= 0` exactly in `Fraction` |

The `class_d` control deserves a sentence, because it is the only way Prop 1 can
be attacked and it is *not* what the note's §B.4 does. Take one edge with

    K(t) = [[0, −a m⁺(t)], [b m⁻(t), 0]],  a, b > 0,

and alternate `(m⁺,m⁻) = (1, ε)` then `(ε, 1)`. Each phase is elliptic (purely
imaginary instantaneous eigenvalues, `±i√(abε)`), yet the two rotations have
orthogonal high-eccentricity axes and their product is hyperbolic — the classic
switched-elliptic instability. This is a genuine counterexample to Prop 1 *as
literally stated* (a forest, `β_ijβ_ji < 0`, "any waveform"), and the only thing
that saves the proposition is (H1), which the note never writes down. In the
Fourier NS model (H1) does hold, because both directions of an edge are mediated
by the single parent amplitude `û_p` (with `û_{−p} = conj û_p`) — but that is a
fact about NS, not about forests, and it belongs in the statement.

---

## 5. Results — the 3-cycle, `C`, and the holonomy claim

### 5.1 `C` is a gauge invariant; the individual `r_e` are not

Diagonal rescaling `x → Θx` sends `K_ij → (θ_i/θ_j)K_ij`, hence

- `r_e = β_bwd/β_fwd  →  (θ_i/θ_j)² r_e`  (**not** invariant),
- `ω_e := √(−β_fwd β_bwd)` (the elliptic frequency of the isolated edge) — **invariant**,
- `F = K₂₁K₃₂K₁₃`, `B = K₁₂K₂₃K₃₁`, hence `C = B/F` and `det K = F + B = F(1+C)` — **invariant**.

This resolves an apparent inconsistency in §B.5/§B.6 of the note. The quoted
"exact forward/backward Leray coefficients" `(−8/5, −1/24, 32/45)` /
`(2/5, 7/5, −59/72)` give raw ratios `(−1/4, −168/5, −295/256)`, which do **not**
equal the quoted `r = (−9/4, −7/3, −59/32)`. But the audit verifies **exactly**
(`Fraction`) that the two sets differ by a legitimate diagonal gauge:

    r_note / r_raw = (9, 5/72, 8/5),   product = 1 exactly
    ⇒ implied θ-ratios (3, 0.263523…, 1.264911…),  product = 1 exactly,

and both sets give the same `C = −1239/128` exactly. **The note's `C` is
confirmed.**

**Complete invariant parametrization.** A 3-cycle generator has 6 real
parameters and a 2-dimensional diagonal gauge group, leaving 4 continuous
invariants plus a discrete one:

    (14)   (ω_a, ω_b, ω_c ; C ; sign F),   |F| = ω_aω_bω_c/√(−C),  det K = F(1+C).

Under (H2) all `r_e < 0`, hence `C = r_ar_br_c < 0` automatically. Canonical
gauge: `r_a = r_b = −1`, `r_c = C` — all the imbalance can be pushed onto one edge.
The audit uses (14) as the scan coordinate.

### 5.2 ERRATUM in §B.6: a spurious `√2` on edge `c`

From the note's own §B.5 coefficients the audit computes, in exact rational
arithmetic,

    ω_a² = 16/25,  ω_b² = 7/120,  ω_c² = 236/405
    ω    = (0.800000, 0.241523, 0.763358)
    c₁   = Σω_e² = 20753/16200 = 1.281049…
    det K = −1111/2700 = −0.411481…
    spectrum(K) = { 0.150054 ± 1.161291 i , −0.300107 } ,  σ_cycle = 0.150054

The note reports `ω_c = 1.0796` and `det K = −0.58192`. The audit finds

    ω_c(note)/ω_c(correct)   = 1.41428
    det K(note)/det K(correct) = 1.41421

i.e. **both are `√2` too large, and by exactly the same factor** — the signature
of a single edge-`c` amplitude carrying a spurious `√2` (rescaling `f_c, b_c` by
`λ` scales `ω_c` and `det K` by `λ` and leaves `r_c` and `C` untouched). The
note's own internal arithmetic is self-consistent given its `ω_c` (re-deriving
its `σ_cycle = 0.149006 ± 1.389416 i` from `c₁ = 1.86387`, `det K = −0.58192`),
so this is a data error upstream of the spectral computation, not an algebra
error. **What is unaffected:** `C`, all `r_e`, `1+C`, Prop 2, and the `C = −1`
null control (all invariant under an edge rescaling). **What is affected:**
`ω_c`, `det K`, `σ_cycle`, `χ_par`, `δ*`, `ρ(M)`, and hence the physical
thresholds `X_P > 60.4 ν s` / `X_P > 170.5 ν s` quoted in §B.6.

Numerically the damage to `σ_cycle` is small by coincidence (`0.150054` vs the
note's `0.149013`, `+0.7 %`) because the `√2` inflates `c₁` and `|det K|` in
partly compensating directions; the damage to `χ_par` and `δ*` is larger
(§5.4).

### 5.3 Proposition 2 — CONFIRMED and sharpened

Under (H2), `c₁ := −Σ_e β_fwd β_bwd = Σ_e ω_e² > 0`, and `tr K = 0`, so the
characteristic polynomial is

    (15)   λ³ + c₁ λ − det K = 0 ,  c₁ > 0 .

With `c₁ > 0` the cubic has **exactly one real root** `μ_r` and one complex pair
`(−μ_r/2 ± iσ)`. `μ_r` is the unique real solution of `μ³ + c₁μ = det K`, so
`sign μ_r = sign det K`, and

    (16)   σ_cycle := max Re λ(K) = { μ_r          if det K > 0
                                    { |μ_r|/2      if det K < 0
                                    { 0            iff det K = 0 iff C = −1.

This confirms the note's Prop 2 ("`C ≠ −1 ⇒` an eigenvalue with `Re > 0`, for
either sign of `det K`") and adds the piece the note omits — **the rate is not
symmetric in the sign**. Linearising (15) for small `|det K|`,

    (17)   σ_cycle ≃ |det K| / c₁ · { 1    if det K > 0
                                    { 1/2  if det K < 0 .

The audit's 120-point `C`-scan (at fixed `ω = (0.8, 0.24152, 0.76336)`) measures
`σ_cycle·c₁/|det K| → 0.99999999999999989` as `C → −1⁺` (`det K > 0`) and
`→ 0.49999999999999956` as `C → −1⁻` (`det K < 0`).
**The factor-2 asymmetry is confirmed to 15 digits.** The B.5
gadget sits on the `det K < 0` side, i.e. the *slow* branch: it gives away a
factor 2 relative to a same-`|det K|` cycle with the opposite orientation. That
is a concrete, checkable design lever the note does not use.

### 5.4 The corrected precise cycle condition

Putting (14), (16), (17) together, with `K = X_P K̂` (parent amplitude factored
out) and children on a **single shell** `|k_i|² = |k|²` — so that `−νD = −ν|k|²I`
is scalar and the Floquet exponents of `A` are those of `K` rigidly shifted:

> **Corrected cycle condition.** For an elliptic 3-cycle (all `r_e < 0`, so
> `C < 0`):
>
> 1. **Existence of growth is a codimension-one condition, with no magnitude
>    threshold in `C`:** the steady cycle has `σ_cycle > 0` **iff `C ≠ −1`**,
>    and `σ_cycle` is a continuous, strictly monotone function of `1 + C` on each
>    side of `−1`, vanishing linearly there. The *sign* of `1 + C` does not
>    decide existence; it decides the rate, via the factor 2 in (17).
> 2. **The magnitude threshold is physical, not structural.** Net growth of the
>    child band requires
>
>        (18)   σ_cycle(K̂) · X_P > ν|k|²,
>
>    i.e., in the small-imbalance regime and using `|F̂| = ω̂_aω̂_bω̂_c/√(−C)`,
>
>        (19)   |1 + C| / √(−C) · ω̂_aω̂_bω̂_c / (ω̂_a²+ω̂_b²+ω̂_c²) · X_P
>                 >  ν|k|² · { 1  if det K > 0 ; 2  if det K < 0 } .
>
>    Only through (19) does the *magnitude* of `C` enter, and it enters weakly.
>    At fixed `ω̂_e`, `|det K| = ω̂_aω̂_bω̂_c·|1+C|/√(−C) ≈ ω̂_aω̂_bω̂_c·√(−C)` for
>    large `|C|`, so `σ_cycle` grows only like `√|C|` while `|det K| ≪ c₁^{3/2}`,
>    and like `|C|^{1/6}` beyond (there `μ_r ≈ (det K)^{1/3}`). Measured on the
>    `C`-scan at fixed `ω`: `C = −2.07 → σ_cycle = 0.0427`;
>    `C = −10.2 → 0.1542`; `C = −50 → 0.3078` — a 24-fold increase in `|C|` buys
>    a 7.2-fold increase in growth rate. **There is no `|C|` threshold and no
>    `|C|` windfall.**
> 3. **Staggered (zero-instantaneous-margin) version.** With one edge active per
>    third of the period, `χ_par(δ) = log ρ(M)/(3δ)`,
>    `M = e^{G_cδ}e^{G_bδ}e^{G_aδ}`:
>    - `det M = 1` because every `G_e` is traceless, so **`ρ(M) ≥ 1` always** and
>      `χ_par ≥ 0` is *structural*, not evidence of anything;
>    - `χ_par(δ) → σ_cycle/3` as `δ → 0` (Trotter/averaging limit, `Ḡ = K/3`);
>    - `χ_par ≡ 0` for all `δ` **iff `C = −1`** (then every `G_e` is
>      `Θ`-antisymmetric for one common `Θ`, each `e^{G_eδ}` is `Θ`-orthogonal,
>      and so is the product — a *theorem*, replacing the note's numerical null
>      control);
>    - for `C ≠ −1`, `χ_par(δ) > 0` on the whole scanned range `δ ∈ [10⁻³, 30]`
>      (400 random cycles: `min_δ χ_par = 3.43 × 10⁻⁸`, median `8.5 × 10⁻⁵`,
>      zero instances of `χ_par < −10⁻⁹`, zero instances of `ρ(M) < 1`). No
>      re-entrant elliptic windows were found, but this is a **numerical**
>      statement over a finite `δ`-grid, not a proof.

### 5.5 "Pure holonomy" — DEFLATED

The note frames the staggered gain as holonomy: "three elliptic rotations …
compose to a hyperbolic map". Quantitatively, the audit finds for the B.5 gadget
(corrected coefficients, `δ`-grid of 900 points on `[10⁻³, 30]`):

| quantity | value | note's value |
|---|---|---|
| `σ_cycle` (steady cycle, all edges on) | `0.1500537` | `0.149013` |
| `χ_par(δ → 0)` | `0.0500261` | — |
| `σ_cycle / 3` (commutator-free Trotter floor) | `0.0500179` | — |
| `max_δ χ_par` | `0.0544523` | `0.052773` |
| `δ*` | `0.92928` | `0.6593` |
| `ρ(M(δ*))` | `1.163933` | `1.1100` |
| `χ_par / σ_cycle` | `0.362885` | `0.354` |
| **holonomy excess over the Trotter floor `σ_cycle/3`** | **`+8.87 %`** | (`+6.2 %` from the note's own numbers) |
| `min_δ χ_par` over `δ ∈ [10⁻³, 30]`, 900 pts | `+2.73 × 10⁻⁵` | — |
| `min_δ ρ(M)` over the same grid | `1.000150` | — |

`χ_par(δ→0)` reproduces `σ_cycle/3` to `2 × 10⁻⁵`. The optimal stagger beats the
commutator-free floor by only ≈ 8.9 %. **So ≈ 92 % of the staggered
Floquet gain is simply the arithmetic mean of the three generators — a
first-order averaging effect, not holonomy.** Across the whole `C`-scan
`χ_par/σ_cycle` stays in `[0.347, 0.393]` against the floor `1/3 = 0.3333`, i.e.
the genuine non-commutative contribution is between `4 %` and `18 %` everywhere.

This does not falsify §B.6 — the phenomenon "positive Floquet exponent with
identically zero instantaneous eigenvalue margin" is **real and confirmed** — but
it reframes it. The correct statement is the classical one: *the spectral
abscissa is not an average*; `max Re λ(Σ G_e/3) > (1/3)Σ max Re λ(G_e) = 0`
already at first order. Holonomy is a small correction on top. The note's own
comparison `χ_par/σ_cycle = 0.354` versus "staggering costs 65 % of the cycle
gain" (judge_tractability) should be read as: staggering costs exactly the
factor `1/3` you would expect from duty-cycling, plus a few percent back from
non-commutativity.

Finally, **`ρ(M) > 1` is not by itself evidence**: `det M = 1` forces `ρ(M) ≥ 1`
for any traceless staggered generator set. The informative quantities are (i) the
null control `C = −1 ⇒ ρ = 1` (which the audit upgrades to a theorem) and
(ii) the *size* of `χ_par` relative to `ν|k|²`.

---

## 6. Corrected minimal Proposition 1

> **Proposition 1′ (minimal hypotheses).** Let `ẋ = (K(t) − νD)x` on `ℝⁿ`,
> `D = diag(|k_i|²)`, `K(t)` real and zero-diagonal with support graph `G`.
> Assume
>
> **(H1)** *reciprocal modulation*: for every edge `e = {i,j}` of `G` there is one
> locally integrable scalar `m_e : ℝ → ℝ` and constants `β_ij, β_ji` with
> `K_ij(t) = β_ij m_e(t)`, `K_ji(t) = β_ji m_e(t)`;
> **(H2)** *elliptic edges*: `β_ij β_ji < 0`, so `ρ_ij := −β_ji/β_ij > 0`;
> **(H3)** *balance*: `∏_{e ∈ γ} ρ_e^{±1} = 1` around every cycle `γ` of `G`
> (vacuous if `G` is a forest).
>
> Then `Θ = diag(θ_i) ≻ 0` exists, depending only on `{β_ij}`, with
> `ΘK(t) + K(t)ᵀΘ = 0` for all `t`, and for every choice of the `m_e`
>
>     −ν max_i|k_i|²  ≤  χ_min ≤ χ_max  ≤  −ν min_i|k_i|² ,
>     ‖x(t)‖₂ ≤ √(κ(Θ)) e^{−ν min_i|k_i|²(t−t₀)} ‖x(t₀)‖₂ .
>
> If in addition all `|k_i|² = |k|²` (single shell), then `χ_max = −ν|k|²`
> **exactly** (trace identity (12)).
>
> Each of (H1), (H2), (H3) is necessary: dropping any one of them admits explicit
> counterexamples with `χ_max > 0` (§4.2 — measured margins `+534`, `+589`,
> `+562` respectively, against `+7 × 10⁻¹⁰` when all three hold).

*Proof:* §2.3 (existence of `Θ`), §2.5 (necessity), §2.7 (single-shell pinning). ∎

**Status of the modulation classes:**
(a) fixed edge ratios `K(t)=m(t)K₀` — **PROVEN**;
(b) common scalar modulation (any waveform, any period) — **PROVEN**;
(c) fully independent per-edge modulations `m_e(t)` — **PROVEN**.
There is no residual "unproved general-modulation case" *within the class where
(H1) holds*: the `W`-based proof of §2.3 closes all of it at once, because `Θ`
never sees the modulation. The unproved/false region is precisely (H1)-violating
per-direction modulation, and there the proposition is **false** (§4.2), so the
right move is to keep (H1) as a hypothesis and discharge it separately from the
NS convolution structure (§7).

**Lean-ability.** Proposition 1′ is finite linear algebra plus Grönwall, in the
style of `formal/NSSingularity/MesoscopicDuhamelNoGo.lean`: construct `Θ` by
forest propagation, prove `ΘK + KᵀΘ = 0` entrywise, then a single differential
inequality. The `C = −1` upgrade (a balanced cycle is dead) is the same proof
with one extra `∏ρ_e = 1` hypothesis. This supports the note's own obligation
G.1 and the judges' "new exact no-go" reading — **with the corrected hypothesis
list**.

---

## 7. Residual gaps (honestly recorded, not audited here)

1. **(H1) from NS.** The audit *assumes* the reciprocal structure and shows it is
   load-bearing. It does **not** re-derive it from `N_k = −iP_k Σ_{ℓ+m=k}(m·û_ℓ)û_m`.
   The real-amplitude/frozen-polarization reduction used in §B.1–B.2 satisfies
   (H1) by inspection; the **complex-parent-phase** case (where
   `K_ij ∝ conj(û_p)`, `K_ji ∝ û_p`) requires the anti-Hermitian version
   `θ_i β_ij = −θ_j conj(β_ji)`, which is again time-independent — so (H1)
   plausibly survives — but this was not verified against the exact `fourier_torus`
   kernel. **Registered as the top follow-up.** If it fails, Prop 1 does not
   apply to the repo's relays and the whole "retro-explains every relay failure"
   claim (§E.14) collapses.
2. **The 512-orientation table** (§B.5: 128 with `C = −1`, 192 "live and every
   edge stable") was **not** re-derived; it needs the exact `fourier_torus`
   enumeration, which is a separate workstream.
3. **Judge criticism not addressed here.** `judge_rigor`'s
   `c_E ≳ ν²N³/χ_par² → ∞` and the fixed-cardinality objection (kill condition
   K2) are *scaling* objections downstream of Prop 1/Prop 2 and are untouched by
   this audit. Nothing found here rescues them, and the `√2` erratum of §5.2
   moves the quoted amplitude thresholds `X_P > 60.4νs` / `170.5νs`, which feed
   directly into that estimate — they should be recomputed before being quoted
   again.
4. **`δ`-scan is finite.** "`χ_par > 0` for all `δ`" is asserted only on
   `[10⁻³, 30]` over 900 (gadget) / 120 (random) points.
5. **Non-single-shell cycles.** Everything in §5 assumes `D = |k|²I`. Off a
   single shell the shift is not rigid and `σ_cycle` versus `ν min|k|²` is no
   longer the right comparison.

---

## 8. Bottom line

Proposition 1 is **not a theorem about forests**. It is a theorem about a
constant diagonal metric, and the metric exists iff (H1) reciprocity, (H2) edge
ellipticity, and (H3) cycle balance hold. Forest topology is one sufficient
condition for (H3) and is neither necessary (balanced cycles are equally dead)
nor sufficient on its own (a forest with per-direction modulation grows without
bound). Under the corrected hypotheses the proposition is **proven for all three
modulation classes at once** — the strongest form the note claims, and by an
argument that is Lean-able today. The cycle invariant `C = −1239/128` and
Proposition 2 survive exactly; the escape channel is real. Three quantitative
claims in §B.6 need correction: a spurious `√2` on edge `c`, the "pure holonomy"
framing (≈ 92 % of the staggered gain is first-order averaging), and the
evidential weight of `ρ(M) > 1` (structural, since `det M = 1`).
