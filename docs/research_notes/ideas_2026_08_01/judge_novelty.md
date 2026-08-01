# JUDGE REPORT — 14 candidates, NS certificate lab (2026-08-01)

**Axis key** (0–5): **D**=distance-to-Clay · **N**=novelty (lit + repo) · **C**=closability · **E**=no-go-evasion soundness · **P**=computability · **I**=interval feasibility · **L**=Lean-formalizability · **F**=falsifiability · **G**=genuine connection to singularity/regularity. (Per-axis prose for 14×9 cells exceeds the word cap; the matrix carries the scores, the note carries the binding axis.)

| # | candidate | D | N | C | E | P | I | L | F | G | binding-axis note |
|---|---|---|---|---|---|---|---|---|---|---|---|
|1|aniso-multiscale|3|3|2|3|3|3|3|3|4|**E/C**: limit system is Prandtl, Sobolev-ill-posed; the growth mechanism *is* the ill-posedness. Anisotropic dynamic rescaling already specified in repo `future_search` §4 → partial reinvention; exponent polytope + ceiling a₁−a₃<1/6 is new.|
|2|multi-center (NSDL)|3|4|2|4|3|3|2|4|4|**N/E**: closed-form skew-dyad strain σ(θ,φ) and *capacity-floor ⟹ n_j=O(1)* are new couplings; lacunary gate ρ(1−η)>1+ρ⁻¹ is exactly checkable. Overlaps Moffatt–Kimura (rejected as M-1) and Pelz-type vortex cascades.|
|3|space-fourier-mismatch|3|3|2|3|3|2|2|4|3|**C**: its own δ=0 theorem kills all exact power laws; the repair makes the front a viscous *equilibrium*, so nothing drives it. c_E(τ)→∞ abandons critical normalization used by every no-go it cites.|
|4|pressure-feedback|3|3|2|4|4|3|2|4|4|**N**: "restricted Euler = spherical depolarization" overlaps Chevillard–Meneveau/tetrad closures; the c₁²c₃²=−c₂²(c₁²+c₃²) impossibility and eigenframe-locking family are new.|
|5|lagrangian-degeneracy|3|4|3|3|3|2|4|4|5|**N/G**: γ=1/Θ bridges Lagrangian SVD geometry to Fourier bandwidth — genuinely new coupling. L1/L3 are cheap, Lean-able, and L3 permanently closes the Euler-perturbative lane.|
|6|**parametric-resonance**|3|**5**|3|3|5|4|4|5|3|**N**: Prop 1 (every tree/forest relay has Floquet exponent ≤ −min ν\|k\|², *for all* modulations) is a new exact no-go that retro-explains every repo relay failure; cycle invariant 𝒞 with ρ(M)>1 at identically-zero instantaneous margin is a new escape channel.|
|7|helicity-ledger|2|**2**|4|3|5|4|5|5|3|**N penalty**: T1's coefficient (s\|ℓ\|−s′\|m\|) *is* Waleffe 1992; \|H_N\|≤2NE_N *is* Kraichnan/Brissaud realizability; Beltrami depletion classical. Only D1 (2-colouring the carrier alphabet) is new.|
|8|log-periodic-multitype|3|4|3|3|5|3|4|4|3|**N**: LP-9 (odd-prime divisibility ⟹ exact sumset miss at *every* dyadic level) is new and elementary; log-4 period with log-2 second harmonic is a real discriminator. Executes repo eq (6.11), so novelty-vs-repo is only moderate.|
|9|forcing-controller|2|3|4|4|4|4|4|4|2|**G**: mathematically clean (open-loop necessity via ‖∂_t f‖≳2νN²‖P_{≤K₀}B‖→∞ is new), but targets (D) and explicitly does not improve the hard step. Hygiene, not mechanism.|
|10|critical-quantity (μ_N)|3|4|5|4|5|4|5|5|4|**C**: elementary, Lean-able; two *independent* derivations (flux side + response side) of the same μ≳1 floor. Ceiling: Theorem D is near-vacuous (μ_*≤5.42).|
|11|**monotone-dichotomy (Λ)**|4|4|5|5|5|4|5|4|5|**C/E**: only candidate with a genuinely provable-now core (I.1–I.4, Lemma K) and *no* equation-level flaw found; ∫KD criterion strictly stronger than Serrin (∞,2) and vorticity (3,2); pressure-blind, so LG-4/LG-9 pathologies cannot recur.|
|12|two-scale-core|2|3|4|4|4|3|3|4|3|**G**: results are negative (LD drain law, Lemma-7 impotence). Valuable hygiene; V4/V6 selection of point-packets is close to standard intermittency heuristics.|
|13|moving-front-ode|3|3|3|3|4|3|4|4|4|**N/E**: front-ODE closure overlaps Katz–Pavlović/Cheskidov dyadic-model front analysis. Front–BKM identity ∫‖ω‖_∞dt=(θ/χ√m)log N and "BKM is never numerically observable" are genuinely new and methodologically important.|
|14|front-flow (seed)|3|2|2|**1**|4|3|2|4|4|**E**: see FATAL. The flow itself is just dynamic rescaling written in Fourier — low novelty; it is connective tissue, not a candidate.|

## FATAL flaws (exact equation / no-go collision)

- **front-flow, S2 closure — FATAL.** `e≡e_c` and `F=χξe^{3/2}` are over-determined: capacity `F_cap/F(1)=ξ⁻²` (slope −2) vs requirement `F_req/F(1)=1−β(ξ−1)` (slope −β, β=√2ν/(χ√c_E)≤1 at threshold). Cap < req on all `1<ξ<φ=1.618`, **for every β**. `ξ_max≥2` handoff is unreachable; the "predicted profile" is not a solution.
- **front-flow, T.2 — FATAL on ℝ³.** `∂_sΨ=0` at constant `a` **is literally** a Leray profile (`û=N⁻²Û(k/N)`). Critical wake ⟹ `U≍|x|⁻¹`: local energy bounded and `∫∫|∇u|²≍∫λ⁻¹dt<∞` ⟹ **Tsai1998 ⟹ u≡0**. Escape only via lattice-DSS, an unfiltered class (Seregin/Chae–Wolf).
- **helicity-ledger — novelty-fatal.** T1 = Waleffe (1992) verbatim; capacity bound = Kraichnan realizability. V1 already refuted (2^{J/2} ≫ J).
- **aniso-multiscale — near-fatal.** `M^eff/N³ ≲ τ^{2a₁−a₂−a₃}=τ^{3/16}→0`: L-11 *would* kill it; the escape rests on three unproven hypothesis failures. Plus Prandtl ill-posedness.
- **pressure-feedback — near-fatal.** Its own V3: `min_θ p₂^mat = 1.036 > 1/2` along the *entire* solved family ⟹ the locked states are inconsistent with a materially advected profile; needs an unbuilt Burgers-renewed profile.
- **parametric-resonance — decisive gap (not fatal).** The staggered pump is imposed, not autonomous; the exhibited 6-mode gadget is fixed-cardinality, already rejected by `fixed_cardinality_scaling` (N^{−3/2}). Survival = its own pre-registered K2.
- **log-periodic — structural limit.** LP-7 closure needs irrational `δ/γ=(1−3^{3/4})/(2√2)`: no exact rational carrier can ever witness the closure; P-LP1 can only reach G1–G3.
- **moving-front-ode — accounting error.** §E.8 applies the L-11 constant at critical normalization (`M_eff=mN³`, `E_N=c_E/N`) while the mechanism runs `E_N=c_E N^{a−1}`, a≈0.78.

## Top 5 (novelty-weighted)

1. **parametric-resonance** — highest novelty; supplies a *new no-go* (Prop 1) valuable even if the mechanism dies, plus a new multiplicative (Floquet/holonomy) escape from L-11's additive Duhamel frame.
2. **monotone-dichotomy (Λ)** — the only candidate whose core is provable *today*; no flaw found; strongest closability/Lean/evasion profile; regularity-side, so it is orthogonal insurance.
3. **critical-quantity (μ_N)** — new interpolation `‖u‖_{X^{-1}}≤3μ_*(‖u‖₂‖∇u‖₂)^{1/2}`; its kill condition can retire the entire mesoscopic cloud lane.
4. **lagrangian-degeneracy** — genuinely new Lagrangian↔Fourier bridge; L1/L3 are cheap, Lean-able, permanent.
5. **multi-center (NSDL)** — best *mechanism* candidate: exact strain law, capacity floor forcing bounded multiplicity, support-based (not amplitude-based) cross-talk gate, and a proof that S6(a) is unfixable.

## Best 1-hour pilot

**Parametric-resonance P1+P2: exact cyclic-carrier enumeration + monodromy with the 𝒞=−1 null control.**
Enumerate every closed child triangle `{k₁,k₂,k₃}` on one shell with parents `p_a+p_b+p_c=0`; compute `r_a,r_b,r_c`, `𝒞=r_ar_br_c`, `det K=β_aβ_bβ_c(1+𝒞)` in exact `Fraction` arithmetic via `fourier_torus.{TrigVector,advection,leray}` (verified present), enumeration style from `exact_carrier_search.canonical_waves_in_box`/`primitive_polarizations`, independent re-check through `exact_carrier_record_verifier`. Then RK4 the staggered monodromy `M=e^{G_cδ}e^{G_bδ}e^{G_aδ}` and report `ρ(M^n)^{1/n}`, n=1,2,3.
**Why it wins:** ~20 min, exact rational (no float in the load-bearing claim), a *mandatory null control* (`𝒞=−1 ⟹ ρ=1±10⁻¹²`), a second independent code path (TM-14), and both outcomes are publishable — either a growth channel with identically zero instantaneous eigenvalue margin (something no tree can do), or `𝒞≡−1` for all physical Leray triads, which extends Prop 1 into a sweeping new no-go over the repo's entire relay history.
**Runner-up:** `critical-quantity` §F — `μ_N` across `η∈{0.10,0.15,0.20}`, `N∈{16,32,48,64}` using `mesoscopic_cloud_scaling`/`mesoscopic_local_fft` unchanged; `μ_N→0` would kill the only structurally uncondemned mesoscopic family by an independent argument.