# ADVERSARIAL AUDIT — `front_flow_seed.md` (S1–S7)

Auditor: adversarial math auditor. Date 2026-08-01.
Every claim below was re-derived from scratch; nothing is taken on the seed's word.
Labels: **[C]** confirmed, **[E]** error (with correction), **[G]** gap (unjustified but not
wrong), **[R]** refuted (claim is false as stated).

Conventions used throughout (repo's, `constraint_map.md` §4):
`u(x)=Σ_k û_k e^{ik·x}`, `∂_t û_k = 𝒩_k − ν|k|²û_k`, `𝒩_k = −iP_kΣ_{ℓ+m=k}(m·û_ℓ)û_m`,
`P_k = I − k⊗k/|k|²`, `E_N = ½‖u_N‖² = c_E/N`.

---

## S1. Carrier-frame front flow

### S1.1 Term-by-term substitution — **[C] CONFIRMED, every exponent**

Ansatz `û(k,t) = N(t)^{−2}Ψ(k/N(t),s)`, `ds = N²dt`, `ξ = k/N`.

```
∂_t û = −2N^{−3}Ṅ Ψ                       (from ∂_t N^{−2})
        + N^{−2}∇_ξΨ · d(k/N)/dt          d(k/N)/dt = −(Ṅ/N)ξ  ⟹ −N^{−3}Ṅ ξ·∇_ξΨ
        + N^{−2}Ψ_s (ds/dt = N²)          ⟹ +Ψ_s
      = Ψ_s − (Ṅ/N³)(2Ψ + ξ·∇_ξΨ).            ✔ matches seed exactly
−ν|k|²û = −ν N²|ξ|² · N^{−2}Ψ = −ν|ξ|²Ψ.       ✔ O(N⁰)
```

Nonlinear term, exponent audit (this is the one that had to be checked digit by digit):
`(m·û_ℓ)û_m = (Nμ · N^{−2}Ψ(λ))·N^{−2}Ψ(μ) = N^{1−4}(μ·Ψ(λ))Ψ(μ) = N^{−3}(μ·Ψ(λ))Ψ(μ)`
with `λ=ℓ/N, μ=(k−ℓ)/N=ξ−λ`. The Riemann sum `Σ_ℓ → N³∫d³λ` gives exactly
`N³·N^{−3} = N⁰`:

```
𝒩_k → −i P_ξ ∫ ((ξ−λ)·Ψ(λ)) Ψ(ξ−λ) d³λ = −𝒬(Ψ,Ψ)(ξ).
```

So the **Riemann-sum factor N³ is correct and it is exactly what makes 𝒩 O(N⁰)**; the
`N^{−3}` from the three amplitude/derivative factors (`N^{+1}` from `m`, `N^{−2}` twice)
is what it must cancel. `P_k = P_ξ` is scale-invariant ✔. The **sign of 𝒬 is correct**:
with `𝒬 := +iP_ξ∫(...)`, `𝒩 = −𝒬`, and the seed's flow carries `−𝒬` ✔.

Assembling `Ψ_s − a(2Ψ+ξ∇Ψ) = −ν|ξ|²Ψ − 𝒬`:
```
∂_sΨ = a(s)(2Ψ + ξ·∇_ξΨ) − ν|ξ|²Ψ − 𝒬(Ψ,Ψ),   a = Ṅ/N³.   ✔
```

### S1.2 `a = Ṅ/N³ = −½ d(N^{−2})/dt` and the front ODE — **[C]**

`d(N^{−2})/dt = −2N^{−3}Ṅ` ⟹ `−½ d(N^{−2})/dt = Ṅ/N³ = a` ✔.
`N^{−2}(t) = N₀^{−2} − 2∫₀ᵗ a dt` ✔ (direct integration).
`a ≥ a₋ > 0 ⟹ T ≤ N₀^{−2}/(2a₋)` ✔. `a→a₊ ⟹ N = (2a₊(T−t))^{−1/2}`, γ = 1/2 ✔.
`dN/ds = (dN/dt)/N² = aN` ✔ exponential. `T−t(s) = ∫_s^∞ N^{−2}ds′ < ∞` ✔ (`dt = N^{−2}ds`,
`N^{−2} ~ e^{−2a₊s}`). Doubling period `S = log2/a₊` from `e^{a₊S}=2` ✔.

### S1.3 **[E] The RG step over-determines τ: `τ = log 2 / a` is forced, τ is not free**

The doubling pullback is the time-`S` flow of the dilation generator `(2+ξ∂_ξ)`, i.e.
`D_μ` with `μ = e^{aΔs}`; requiring `μ = 2` gives `aΔs = log 2`. The frozen-parent Duhamel
step runs for physical `Δt = τN^{−2}`, i.e. `Δs = τ`. **Hence a·τ = log 2.** The pilot (S7)
scans `(η,c_E,τ,ν)` as if τ were independent; it is not — once you demand a doubling
pullback, `a₊ = log2/τ` is an *output identity*, and combined with S3's `a ∝ √c_E`
(below) it forces `τ ∝ c_E^{−1/2}`. Corollary prediction: `T = τN₀^{−2}/(2log2)`.

### S1.4 **[E]/[G] S1 and S6 are mutually exclusive regimes**

The Riemann sum `Σ_ℓ → N³∫dλ` requires the profile to be *dense* on the lattice —
`M_N ≍ N³`, i.e. box half-width `ηN` (γ = 1). S6's sumset-miss / cross-talk geometry is
an argument about *sparse* integer carriers. Both cannot hold at the same N. The seed
never reconciles them. (This is not fatal — the repo's own (6.5)–(6.6) resolve it by using
finitely many carriers `a∈𝒜` each carrying a *continuum profile* `U_a(ξ)` on a cell `Ω_η`
— but then the "Riemann sum" is per-cell, and the S1 flow must be written cell-wise, which
changes `𝒬` into the repo's (6.6) with the `Σ_{a+b=r}` bookkeeping. The seed's identification
"𝒬 restricted to carrier cells IS eq (6.6)" is **[C] correct** modulo this rewriting.)

### S1.5 **[G] `a` is an eigenvalue, not a parameter**

Nothing in S1 says how `a(s)` is determined. In the front flow `a` is fixed by the
normalization constraint (fixing `∫e dξ`-type moment or the position of the front), i.e.
it is a nonlinear eigenvalue. This is the usual modulation/gauge condition and is missing.

### S1.6 **[G] "attracting periodic orbit found by forward integration" is unsupported**

Self-similar blowup profiles are generically saddles. Forward RG integration *with energy
renormalization at each stage* is a power method: it removes the leading (amplitude/scaling)
unstable direction — which is exactly what (6.10)'s projective normalization does — but no
argument is offered that the remaining shape spectrum is contracting. Repo PO-08 lists
exactly this as strategy-less. Claim is a hope, not a derivation.

---

## S2. 1D shell reduction

### S2.1 Sphere-average bookkeeping and the coefficient 2 — **[C] EXACT**

With `g(ξ) := ∮_{S²}|Ψ(ξω)|²dω` and `e = ξ²g`:
```
dilation contribution to ∂_s e = ξ²∮ 2Re(Ψ̄·a(2Ψ+ξ∂_ξΨ))dω = a(4ξ²g + ξ³g′)
a(2e + ξe′)   with e′ = 2ξg + ξ²g′  ⟹  a(2ξ²g + 2ξ²g + ξ³g′) = a(4ξ²g + ξ³g′).  ✔ IDENTICAL
```
Structural reason (worth recording): conjugating the generator `(4+ξ∂_ξ)` acting on `|Ψ|²`
by the Jacobian `ξ²` gives `ξ²(4+ξ∂_ξ)ξ^{−2} = 2 + ξ∂_ξ`. **The "2" is correct.**
Viscous: `ξ²∮2Re(Ψ̄(−ν|ξ|²Ψ)) = −2νξ²e` ✔.
Nonlinear written as `−∂_ξF` with `F(0)=F(∞)=0`: legitimate **as a definition**
(`F(ξ):=−∫₀^ξ(nonlinear)`); `F(∞)=F(0)` follows from continuum energy neutrality
`∫Re(Ψ̄·𝒬)d³ξ=0`, and `F(0)=0` from regularity at ξ=0. **[C]** as a definition; F is *not*
determined by e — the closure is a genuine ansatz, correctly labelled.

### S2.2 Dilation term vanishes on `e_c = 2c_Eξ^{−2}` — **[C] EXACTLY**

`2e_c + ξe_c′ = 4c_Eξ^{−2} + ξ(−4c_Eξ^{−3}) = 0` ✔ identically. Good, non-trivial: the
critical spectrum is the kernel of the dilation operator on the shell density.

### S2.3 Factor-2 in the front-energy formula — **[E] but self-cancelling**

Seed writes "front energy = `N^{−1}∫e dξ`". With `E = ½Σ_k|û_k|²` the Riemann sum gives
`E = ½N^{−4}·N³∫|Ψ|²d³ξ = (2N)^{−1}∫e dξ`. **Seed is off by a factor 2.**
Consequence chain if corrected: dyadic-shell density `ε(k)=2c_E/k²` (from `∫_N^{2N}ε=c_E/N`),
so `e_c = 4c_Eξ^{−2}`, `F′=−8νc_E`, `F(1)=χ(4c_E)^{3/2}=8χc_E^{3/2}`, `ξ_max = 1+χ√c_E/ν`,
threshold `χ√c_E ≥ ν`. But `χ` is defined *relative to the same e-normalization*
(`χ_new = χ_old/√2`), so **the threshold `χ√c_E ≥ √2ν` is invariant** ✔. Net: bookkeeping
slip with no structural consequence, but any numerical comparison of χ against a measured
flux must state which normalization is in force.

### S2.4 `F′=−4νc_E`, `ξ_max`, threshold — **[C] algebra**

On the wake (dilation ≡ 0, steady): `0 = −2νξ²e_c − F′ ⟹ F′ = −2νξ²·2c_Eξ^{−2} = −4νc_E` ✔.
`F(ξ)=F(1)−4νc_E(ξ−1)` ✔. `F(1)=χ·1·(2c_E)^{3/2}=2√2χc_E^{3/2}` ✔.
`ξ_max = 1 + F(1)/(4νc_E) = 1 + (2√2/4)χ√c_E/ν = 1 + (χ/√2)√c_E/ν` ✔.
`ξ_max ≥ 2 ⟺ χ√c_E ≥ √2ν` ✔.

### S2.5 **[E] MAJOR — the two closures are mutually inconsistent; the "predicted measurable profile" is not a solution**

The derivation uses `e ≡ e_c = 2c_Eξ^{−2}` (to kill the dilation term and get
`F′=−4νc_E`) *and simultaneously* `F = χξe^{3/2}` with the same `e`. But
`χξ e_c^{3/2} = χ(2c_E)^{3/2}ξ^{−2} = F(1)ξ^{−2}`, whose derivative is `−2F(1)ξ^{−3}`, which
equals `−4νc_E` only at ξ=1 and only if `F(1)=2νc_E` ⟺ `χ√c_E = ν/√2` — **a factor 2 BELOW
the required threshold `χ√c_E ≥ √2ν`.** The system is over-determined; the line
`e(ξ) = [(F(1)−4νc_E(ξ−1))/(χξ)]^{2/3}` is inconsistent with the equation used to derive
`F(1)−4νc_E(ξ−1)`, and is **not** a "predicted measurable profile".

Quantitative form of the damage. Write `β := 4νc_E/F(1) = √2 ν/(χ√c_E)` (threshold ⟺ β≤1).
- *required* flux to sustain the critical wake: `F_req(ξ)/F(1) = 1 − β(ξ−1)`
- *capacity* of the closure on the critical wake: `F_cap(ξ)/F(1) = ξ^{−2}`

At β = 1 (exact threshold): `ξ^{−2} = 2−ξ ⟺ ξ³−2ξ²+1=0 ⟺ (ξ−1)(ξ²−ξ−1)=0`, roots
ξ = 1 and **ξ = φ = 1.618…** (golden ratio). For **1 < ξ < 1.618 the capacity is strictly
below the flux required to hold the critical wake** (e.g. ξ=1.2: cap 0.694 < req 0.800).
So the critical wake cannot be maintained on the very first sub-octave, and `ξ_max = 2` is
*not* reached by any solution of the closed system.

**Correction / what is actually true.** Only one of two readings survives:
(a) *Necessary-condition reading* (charitable, and the only defensible one): `F_req` is the
flux the critical wake demands, `F_cap = χξe^{3/2}` is what a shape can deliver; then
`χ√c_E ≥ √2ν` is **necessary** but **not sufficient**, and the extra sufficiency condition
`F_cap(ξ) ≥ F_req(ξ) for all 1<ξ<ξ_max` **fails for every β** (both curves start at F(1),
capacity has slope `−2F(1)`, requirement has slope `−βF(1)`, and β≤1<2).
(b) *Honest steady-state reading*: keep the dilation term and solve
`a(2e+ξe′) − 2νξ²e − d/dξ(χξe^{3/2}) = 0`, a nonlinear first-order ODE with `a` as
eigenvalue. `e_c` is **not** a solution of this ODE. The wake must overshoot
(`e > e_c` for ξ slightly >1) to carry the required flux; that overshoot is precisely the
`s`-dependent "front bump" and it re-activates the dilation term. **This is the real
mathematical content and it is not derived anywhere in the seed.**

### S2.6 Dimensional consistency with the repo capacity law — **[C] and it pins γ=1**

Repo: `|Π_N| ≤ C N √(M_N) E_N^{3/2}` with `E_N=c_E/N` ⟹ `Π ≤ C√(M_N) c_E^{3/2} N^{−1/2}`.
Seed: `Π_phys = N·F`, `F(1)=2√2χ c_E^{3/2}` ⟹ `Π_phys = 2√2 χ c_E^{3/2} N`.
Matching: `χ ≲ C √(M_N) N^{−3/2}`. With `M_N ≍ (ηN)³`: **`χ ≲ Cη^{3/2}`, scale-independent.**
✔ Consistent — and it shows the closure constant χ can be an N-independent shape number
**only when `M_N ≍ N³`** (γ=1). For γ<1, `χ ≲ Cη N^{3(γ−1)/2} → 0`, matching the L-11
rejection of every γ<1. Good internal coherence.
*Caveat*: the conversion `Π_phys = N·F` carries the same factor-2 slip as S2.3
(`Π_phys = NF/2` in the repo's `E=½‖u‖²` convention); c_E- and N-scalings are unaffected.

---

## S3. c_E-collapse

### S3.1 The collapse is real, and stronger than the seed's dimensional counting — **[C], upgraded**

The front flow has an **exact one-parameter scaling group**: if `Ψ` solves with `(a,ν)`,
then `Φ = μΨ`, `σ = s/μ` solves the same equation with `(a,ν) → (μa, μν)`; equivalently
```
∂_σΦ = (a/μ)(2Φ+ξ∇Φ) − (ν/μ)|ξ|²Φ − 𝒬(Φ,Φ).
```
Since `e ~ |Ψ|² ~ c_E`, `μ = √c_E`, so **the whole c_E-family is the single family with
effective viscosity `ν_eff = ν/√c_E`, effective rate `a_eff = a/√c_E`, and s-clock
rescaled by √c_E.** Every "channel ratio is a c_E-independent shape functional" statement
follows immediately. ✔

### S3.2 **The dilation channel: [C], but only because `a ∝ √c_E` — a consequence the seed omits**

The dilation term in the e-equation is `a·(2e+ξe′) ~ a c_E`, which is `c_E^{3/2}` **iff
`a ∝ √c_E`**. This is forced by the scaling group above (`a_eff` must be an O(1) shape
number). Therefore:
- **`a₊ ∝ √c_E`**, so `T = N₀^{−2}/(2a₊) ∝ N₀^{−2} c_E^{−1/2}`,
- doubling period `S = log2/a₊ ∝ c_E^{−1/2}`,
- with S1.3, **`τ = log2/a₊ ∝ c_E^{−1/2}`.** The pilot must scan τ *along* this curve.
Answer to the audit question: yes, the dilation channel is covered, but only after you
accept that `a` is not c_E-independent — which contradicts any reading of S3 in which
c_E is tuned "at fixed dynamics".

### S3.3 **[E]/[R] "large c_E is legal / the problem reduces to shape closure only" is a re-labelling, not a reduction**

`c_E ≥ 2ν²/χ_shape²` ⟺ `ν_eff = ν/√c_E ≤ χ_shape/√2`. So "choose c_E large" **is exactly**
"take the effective-Reynolds-number limit". The residual question — does the front flow at
small `ν_eff` admit a recurrent profile with positive net forward flux — is the
inviscid/high-Re limit, i.e. the discretely-self-similar *Euler* blowup problem. Calling
this "shape closure only" understates it by exactly the amount that makes the problem hard.
The reduction is correct as algebra, **[C]**, but the claim that it makes the lane cheap is
**[R]**.

### S3.4 Does large c_E break smoothness or the energy budget? — **[C] no, with three sharp caveats**

- *Smoothness*: a band-limited field on 𝕋³ is a trig polynomial, hence real-analytic. No
  obstruction. On ℝ³ this fails (Z-01/Z-02 kinematic incompatibility) — the seed is a torus
  object only.
- *Energy*: `Σ_{j≥0} c_E/(2^jN₀) = 2c_E/N₀` ✔ arithmetic correct; small by taking N₀ large.
  **[C]**, provided the wake's IR edge is the *initial* shell N₀ (not k=1); if the wake
  extended to k=1 the total would be `≈2c_E`, O(1). The seed's statement is correct as
  written but only because the cascade is UV-directed.
- *Critical norm*: `‖u₀‖²_{Ḣ^{1/2}} ≍ N₀E_{N₀} = c_E`, i.e. `‖u₀‖_{Ḣ^{1/2}} ≍ √c_E`. The
  Fujita–Kato small-data threshold is `‖u₀‖_{Ḣ^{1/2}} ≲ ν`. **The threshold
  `χ√c_E ≥ √2ν` is, up to the O(1) factor χ, exactly the statement that the datum is
  supercritical in Ḣ^{1/2}.** So S3(i) is legitimate but *content-free*: it says "pick data
  outside the small-data theory". **[C] + [G]** — this deserves to be stated in the seed
  because it shows the threshold is not an artifact of the closure.
- *Dissipation budget*: `ν∫₀ᵀ‖∇u‖²dt ≈ 2νc_E∫₀ᵀN(t)dt = 2νc_E/(a₊N₀)`; requiring this
  ≤ available energy `2c_E/N₀` gives **`ν ≲ a₊`** ✔ — the *same* threshold again
  (`a₊ ∝ √c_E`). Three independent routes to one condition: strong internal consistency. **[C]**

### S3.5 **[G] `𝒜 ∝ c_E` is a first-order-Duhamel fact and cannot be extrapolated to 𝒜 ≈ 1/2**

`D_N = 2c_E N^{−1}H_N G_N²` is linear in c_E only in the *frozen-parent* (first Duhamel
iterate) approximation. The repo's "required c_E = 228/546/1902" numbers are obtained by
linear extrapolation from `D_N ~ 2×10^{−3}` to `D_N = 1/2` — a factor 250–2000 extrapolation
of a perturbative formula into the fully nonlinear regime. **Those numbers are not
predictions.** (The one full-Galerkin cross-check, ratio 0.968, was at `D ≈ 2×10^{−3}`,
deep in the linear regime, and says nothing about `D ≈ 1/2`.)

---

## S4. Flux → L³ per-octave lemma

### S4.1 Hölder/Bernstein chain — **[C], constants check out**

```
|Π_N| = |⟨P_N ℙ((u·∇)u), u_N⟩| ≤ ‖(u·∇)u‖_{3/2} ‖P_N ℙ u_N‖_3        (Hölder 2/3+1/3=1)
‖(u·∇)u‖_{3/2} ≤ ‖u_band‖_3 ‖∇u_band‖_3                             (Hölder 1/3+1/3=2/3)
‖∇u_band‖_3 ≤ κN‖u_band‖_3                                          (Bernstein, band-limited)
‖P_N ℙ u_N‖_3 ≤ C‖u_N‖_3                                            (Riesz/LP bdd on L³, 1<3<∞)
⟹ |Π_N| ≤ CκN‖u_band‖_3³.   ✔
```

### S4.2 Is band-locality circular? — **[G] no, but it hides the sweeping commutator**

Not circular: it is an extra hypothesis, not a consequence of the conclusion. But it is
doing real work. The dangerous piece is the low-high (sweeping) contribution
`⟨(u_{low}·∇)u_N, u_N⟩`, which is *exactly zero* without projections
(`= ∫u_low·∇(|u_N|²/2) = 0`) but only *commutator-small* once `P_N` is inserted. For a field
that is genuinely band-limited to an annulus `[N/κ₀,κN]` (trig polynomial on 𝕋³) the chain
above is valid as written. For a Littlewood–Paley decomposition of a general field the
constant C hides a Bony paraproduct/commutator estimate. **On the torus lane this is fine;
it must not be exported to ℝ³ packets (Z-01/Z-02).**

### S4.3 The c_E and N exponents — **[C] after fixing a notational inconsistency**

Bernstein: `‖u_N‖₃ ≤ CN^{1/2}‖u_N‖₂ = CN^{1/2}√(2c_E/N) = C√(2c_E)` — **N-independent**
(L³ is critical) and `‖u_N‖₃³ ≍ c_E^{3/2}` ✔, saturating the inequality's c_E-scaling.
The seed writes both `Π_N ≥ q c_E N` and `Π_N ~ χ(2c_E)^{3/2}N`; these force
`q = 2^{3/2}χ√c_E`, i.e. **q is not a pure shape constant** — it carries `√c_E`, consistent
with `a ∝ √c_E` (S3.2). With that q, `‖u_band‖₃³ ≥ (q/Cκ)c_E = 2^{3/2}χc_E^{3/2}/(Cκ)` ✔.
The seed's display `≥ (q/C)c_E^{3/2}` is only right under this reading. **[E] notational**,
resolvable.

### S4.4 **[E] The per-octave summation step is missing — but it is repairable, and the repair is clean**

`‖u‖₃³` is **not** additive over Littlewood–Paley pieces; the seed asserts
`‖u‖₃³ ≳ const·#octaves` with no argument. Repair:
```
‖u‖₃ ≳ c‖(Σ_j|u_j|²)^{1/2}‖₃          (LP square-function lower bound, 1<3<∞)
‖(Σ_j|u_j|²)^{1/2}‖₃³ = ‖Σ_j|u_j|²‖_{3/2}^{3/2} ≥ ∫Σ_j|u_j|³ = Σ_j‖u_j‖₃³
```
the middle step by the pointwise elementary inequality `(Σa_j)^p ≥ Σa_j^p` for `p=3/2 ≥ 1`,
`a_j = |u_j|² ≥ 0`. **So `‖u‖₃³ ≥ c³ Σ_j‖u_j‖₃³` with a universal constant.** ✔ This is a
genuine (small) new lemma and is Lean-able as stated. With it, S4's conclusion
`‖u(t)‖₃³ ≳ c_E^{3/2} log₂(N(t)/N₀) → ∞` **[C]**, *conditional on* the flux hypothesis.

### S4.5 Wake persistence (heat multiplier positivity) — **[C], and the constant is exactly scale-free**

Claim: each completed octave retains a constant energy fraction. Verify: octave j completes
at `t_j` with `N(t_j)=N_j`; the residual lifetime is `T−t_j = N_j^{−2}/(2a₊)` (exact under
`N=(2a₊(T−t))^{−1/2}`), so the total heat decay factor is
```
exp(−ν N_j²(T−t_j)) = exp(−ν/(2a₊)),
```
**independent of j** ✔ — the retained fraction is a fixed number, exactly because γ=1/2 makes
the remaining lifetime one parabolic time of the shell's own scale. This is correct and is
the sharpest part of the seed. Positivity/realness of `e^{−ν|k|²Δt}` ⟹ no phase rotation ✔.
Combined with the S3.4 budget `ν ≲ a₊`, the retained fraction is `≥ e^{−1/2} ≈ 0.61`.
**[G]**: this bounds only the *linear* depletion; nonlinear export out of the wake is not
bounded, and the whole statement presumes the octave's shape does not decohere.

### S4.6 **[G] "removes proof-chain steps 8–9" is conditional on the hardest open step**

S4 derives ESS-divergence *from* a sustained flux lower bound `Π_N ≥ qc_EN`. The constraint
map names precisely this as "the single largest unproven step across both the
shell-constraint and Zeno-candidate lines" (Bernstein gives only an upper bound; the *sign*
is open). So S4 does not remove obligations; it *re-routes* 8–9 onto step 5. That is real
value, and the seed says so — but the headline "follows automatically" should read
"follows conditionally".

---

## S5. Helicity ledger

### S5.1 Linear polarization has exactly zero helicity — **[C], pointwise**

`u = A e cos(k·x+φ)`, `e·k=0`, `e,k,A` real. Then `ω = ∇×u = −A(k×e)sin(k·x+φ)` and
```
u·ω = −A²(e·(k×e)) cos·sin ≡ 0     since (k×e) ⊥ e.
```
So the helicity **density** vanishes identically, not merely its integral ✔. A fortiori all
current linearly-polarized relay alphabets have H ≡ 0. **[C]**

### S5.2 Beltrami depletion — **[C]**

`∇×u = λu ⟹ (u·∇)u = ∇(|u|²/2) − u×(∇×u) = ∇(|u|²/2) − λ(u×u) = ∇(|u|²/2)`, and
`ℙ∇(anything) = 0`. So `B(u,u) = 0` exactly ✔.

### S5.3 Helicity capacity `|H_N| ≤ 2NE_N = 2c_E` — **[C]**

`|ω̂_k| = |k×û_k| ≤ N|û_k|` on the shell ⟹ `|H_N| ≤ N‖u_N‖₂² = 2NE_N = 2c_E` ✔.

### S5.4 **[E] "coherent transfer requires helicity-neutrality" is a non sequitur**

The implications proved are: *maximal* helicity (Beltrami, `|h| := |H_N|/(2NE_N) = 1`)
⟹ zero transfer; linear polarization ⟹ `h = 0`. Neither gives "transfer ⟹ h = 0".
**Correct statement**: the flux efficiency obeys `χ = χ(h)` with `χ(±1) = 0`, so the
threshold `χ√c_E ≥ √2ν` requires `|h|` bounded away from 1 — nothing more. The design rule
"use zero-helicity (linearly polarized) carriers" is **not** implied; "use non-Beltrami
carriers" is.

### S5.5 **[R] The suggested no-go for single-sign helical critical cascades FAILS via the route proposed, but a marginal balance exists via a different route**

- *Proposed route (Cauchy–Schwarz)*: `|H| ≤ ‖u‖₂‖ω‖₂`. Here `‖u‖₂ ≍ 2√(c_E/N₀)`,
  `‖ω‖₂ ≍ √(2c_E N)` ⟹ bound `≍ 2√2 c_E √(N/N₀) = 2√2 c_E 2^{J/2}`. The claimed accumulation
  is `2c_E J`. **Linear vs exponential in J — the bound is slack by an exponential factor
  and constrains nothing.** The suggested no-go is refuted. **[R]**
- *The route that does bind*: helicity is not free — `dH/dt = −2ν∫ω·(∇×ω)`, so
  `|dH/dt| ≤ 2ν‖ω‖₂‖∇ω‖₂ ≍ 2√2 νc_E N²`. Since `N² = (2a₊(T−t))^{−1}`,
  `∫₀^{T}|dH/dt|dt ≍ (√2 νc_E/a₊)·2log N` — **log-divergent, i.e. the same order as the
  demanded accumulation `2c_E log₂N`.** Comparing coefficients:
  single-sign accumulation is affordable iff `ν/a₊ ≳ 1/(√2 log2) ≈ 1.02`.
  But S3.4's energy budget requires `ν/a₊ ≲ 1`. **The two budgets meet with essentially
  zero margin.** Verdict: no theorem, but a sharp heuristic that a *single-sign* maximally
  helical critical cascade is exactly marginal, while an alternating/neutral helicity code
  removes the tension entirely. This *supports* the seed's design rule by a completely
  different argument than the one given.

### S5.6 **[C]** "exact rational helicity columns at trivial cost" — true; `H = Σ_k i k·(û_k × conj(û_k))`-type sums are rational for the rational alphabets already in `exact_carrier_search.py`.

---

## S6. Cross-talk removal by scale-stagger

### S6.1 **[E] The translation-split suppression exponent is wrong / conflated**

The seed says power suppression `~ (W|Δx|)^{−3}`. Two different (both standard) mechanisms
exist and neither gives that exponent:
- *Coherent kernel decay.* The cross child coefficient is `Σ_m g(m)e^{−im·Δx}` over a box of
  half-width W. For a sharp box (Dirichlet kernel) the ratio to the aligned value is
  `~ (W|Δx|)^{−1}` per dimension, i.e. **`(W|Δx|)^{−3}` in AMPLITUDE, `(W|Δx|)^{−6}` in
  POWER.** For a Fejér/triangular envelope it is `(W|Δx|)^{−2}` per dimension, hence
  `(W|Δx|)^{−12}` in power. For a C^∞ envelope it is faster than any power.
- *Incoherent floor.* Random phases give `E|Σg|² = Σ|g|²` vs the coherent `M·Σ|g|²`, so the
  power ratio floors at `M^{−1} = W^{−3}` — **`W^{−3}` in power, with no `|Δx|` dependence**
  once `|Δx| ≳ 1/W`.
**Correct statement**: power suppression `≍ max{(W|Δx|)^{−6}, W^{−3}}` (box envelope),
saturating at the incoherent floor `W^{−3}` for `|Δx| ≳ W^{−1/2}`. The seed's
`(W|Δx|)^{−3}` is the incoherent floor with a spurious `|Δx|` attached.

### S6.2 The stage-2 failure claim — **[C], the symmetry argument is airtight**

Child₁, born from two parents both centered at x₀, inherits the phase `e^{−ik·x₀}`, i.e. it
is localized at x₀; child₂ at x₁. The stage-2 pairing child₁⊗child₂ therefore sits at
separation `|Δx| = |x₁−x₀|` with the *same* box width scale (2W after one convolution), so it
suffers the *same* suppression law as the diagonal cross-talk it was meant to kill. The
symmetry is exact and (a) is correctly **REJECTED** ✔. (Escape the seed names — "redesign
the stage-2 pairing" — would require the stage-2 partners to be co-located, i.e. undoing the
split.)

### S6.3 **[E] "L-periodic multi-type cycle with Πλ_j = 2 per super-period" — misquoted and, for equal ratios, arithmetically impossible**

- The repo's (6.11) requires `Π_j λ_j > 1`, **not `= 2`**. Misquote.
- If all L ratios are equal, `λ^L = 2 ⟹ λ = 2^{1/L}`, **irrational for L ≥ 2**, hence
  impossible on ℤ³ (`λN ∈ ℕ` for all N in the chain is required by (6.8)/lattice
  compatibility). So the "L types each advancing by 2^{1/L}" reading is dead.
- The only realizable version is *unequal* rational ratios with rational product, e.g.
  `λ₀ = 4/3, λ₁ = 3/2` — which is exactly **two interleaved doubling ladders**
  `{2^j}` and `{ρ2^j}`. Under that reading each type still doubles per stage, so the
  "super-period" language is vacuous but the construction is coherent.

### S6.4 **Is the scale-stagger sumset-miss satisfiable? — [C] YES, but only for small η, and ρ=1/2 is degenerate. Hand computation follows.**

Carriers (repo alphabet, `constraint_map.md` §"Second-stage / cross-talk"):
relay-1 `p=(1,1,0)`, `q=(1,0,1)` at scale N, child `(2,1,1)N`;
relay-2 `p′=(0,1,−1)`, `q′=(1,0,−1)` at scale ρN, children `(1,1,−2)ρN`, `(1,1,1)ρN`.
Diagonal cross channels are `N(p±ρp′)`, `N(p±ρq′)`, `N(q±ρp′)`, `N(q±ρq′)` (± because real
fields carry ±k). Boxes have ℓ^∞ half-width `ηN`; a cross output box has half-width `2ηN`
(one convolution broadens W→2W). **Miss condition**: ℓ^∞ distance from a cross output centre
to every tagged centre `> 3ηN`.

Computed ℓ^∞ margins (units of N; `m` = min over all tagged centres at all scales `2^j`,
`ρ2^j`):

| ρ | cross outputs (units of N) | binding neighbour | margin m | ⟹ requires |
|---|---|---|---|---|
| **1/2** | (1,1.5,−.5), (1.5,1,−.5), (1,.5,.5), (1.5,0,.5), (.5,1,.5), (1,−.5,1.5), (.5,0,1.5) | (1,1,0)/(1,0,1)/(2,1,1) | **1/2** | **η < 1/6** |
| **2/3** | (1,5/3,−2/3), (5/3,1,−2/3), (1,2/3,1/3), (5/3,0,1/3), (1,1/3,2/3), (1/3,1,2/3), … | (1,1,0) at (0,1/3,1/3); ρ(1,1,1)=(2/3,2/3,2/3) | **1/3** | **η < 1/9** |
| **3/4** | (1,1.75,−.75), (1.75,1,−.75), **(1,.75,.25)**, … | (1,1,0) at (0,−1/4,1/4) | **1/4** | **η < 1/12** |

Worked example (the binding row): `q + ρp′ = (1,0,1) + ¾(0,1,−1) = (1, ¾, ¼)`; distance to
the parent box `(1,1,0)` is `‖(0,−¼,¼)‖_∞ = ¼`. Similarly for ρ=2/3,
`q+ρp′ = (1,⅔,⅓)`, distance to `(1,1,0)` is `⅓`, and to `ρ(1,1,1)=(⅔,⅔,⅔)` also `⅓`.

**Conclusions.**
1. The sumset-miss **is satisfiable** — the seed's proposal is not empty. ✔
2. But it forces `η ≤ 1/9` (ρ=2/3) or `η ≤ 1/12` (ρ=3/4), well inside the surviving window
   `η<1/3`. Since the transfer floor scales like `η³`, this costs a factor `27–64` in `D_N`
   relative to `η=1/3`, i.e. **required c_E rises by the same factor** — from the repo's
   measured `c_E ≈ 228` at η=0.2 to `≈1.3×10³` (ρ=2/3) or `≈3.2×10³` (ρ=3/4), i.e.
   `ν_eff = ν/√c_E` down by another factor ~4. Scale-stagger buys geometric cleanliness by
   *worsening* the very budget it was invoked to protect.
3. **ρ=1/2 is degenerate and must not be counted as a success.** It has the best margin
   (1/2) but `ρ = 2^{−1}` means "relay-2 is relay-1 one octave earlier" — it is not a
   stagger at all, it is the ordinary parent/grandparent adjacency of a single cascade, and
   the cross-talk it "avoids" is exactly the cross-talk every cascade has. Genuine stagger
   requires `ρ ∉ 2^ℤ`.
4. Number of cross channels grows as (#tagged boxes)², so the integer-programming check the
   seed proposes is necessary; the three-row table above is the ρ-selection, not the full test.

---

## S7. Doubling pullback bookkeeping

### S7.1 **[E] "child band 2k → k subsampling" is the wrong discretization and loses 7/8 of the child**

The map `û_new(k) := α û_old(2k)` samples the child only on the **even sublattice** `2ℤ³`.
The child band is the sumset of two parent boxes and generically populates **all** parities,
so subsampling discards ~7/8 of the child modes/energy. It is exact only if the child is
supported on `2ℤ³`, i.e. if the field is exactly half-period periodic — and periodicity
cannot spontaneously increase under NS, so this cannot be arranged dynamically without
making the field periodic at *every* scale from the start (a dense singular set, which
conflicts with CKN for suitable weak solutions).
**The repo's own map is right and the seed's paraphrase is wrong**: (6.8) is
`(𝔗_λ U)_a(ξ) = λ² V_tot(λ(a+ξ))` — a **continuum dilation of the profile**, evaluated by
*resampling/interpolation*, not lattice subsampling. Any pilot must resample `Ψ` (or run the
next stage at scale 2N with the same `(a,Ω_η)` geometry and compare shapes), never subsample.

### S7.2 The continuum energy bookkeeping is correct — **[C]**

With `û(k)=N^{−2}Ψ_j(k/N)` and the next stage at `2N`, the pullback is
`Ψ_{j+1}(ξ) = 4 Ψ̃(2ξ)` (from `(2N)^{−2}Ψ_{j+1}(ξ) = û(2Nξ) = N^{−2}Ψ̃(2ξ)`). Then
`e_{j+1}(ξ) = 16ξ²∮|Ψ̃(2ξω)|² = 4 e_{Ψ̃}(2ξ)`, so `∫e_{j+1}dξ = 2∫e_{Ψ̃}dη`, and
`E_{j+1} = (1/(2·2N))∫e_{j+1} = (1/(2N))∫e_{Ψ̃} = E[Ψ̃]`. **The pullback is exactly
energy-preserving in profile variables** ✔ — as it must be, being a relabelling. The factor
`λ` in the repo's (6.9) is the same statement.

### S7.3 **[E] "𝒜_j → 𝒜₊ > 0 (ANY positive constant)" is the wrong success criterion as stated**

The repo is explicit right below (6.9): *"a critical doubling relay needs selected child
energy ratio at least 1/2"*. `𝒜₊ > 0` is sufficient only via the S3 c_E-lever
(`D_N ∝ c_E` ⟹ rescale c_E until `𝒜 ≥ 1/2`), and that lever is exactly the perturbative
extrapolation flagged in S3.5. So the honest criterion is: **`𝒜₊ > 0` at fixed c_E is a
screen; `𝒜 ≥ 1/2` at the self-consistent c_E is the gate**, and the two are connected only
by an unvalidated linear extrapolation. Also `ρ_j → 1` must be measured in a **critical**
(not L²) norm — see S7.5.

### S7.4 **[G] τ is not a free scan parameter** — see S1.3: `a₊τ = log2`, and `a₊ ∝ √c_E`
forces `τ ∝ c_E^{−1/2}`. The stated 4-parameter scan `(η,c_E,τ,ν)` has (at most) 3
independent directions, and one of the three published constraint combinations
(`χ√c_E ≥ √2ν`, `ν ≲ a₊`, `a₊τ = log2`) is redundant with the others.

### S7.5 **[G] The front profile is unbounded in L², so no L²-norm orbit exists**

The wake's IR edge sits at `ξ = N₀/N → 0`, so `∫e dξ ≍ 2c_E N/N₀ → ∞`: `Ψ` does **not**
converge in the natural `∫e dξ` (L²-profile) topology. Periodicity/attraction can only be
sought (i) modulo the wake tail, in a topology weighted so that `e_c = 2c_Eξ^{−2}` has finite
norm (e.g. `sup_ξ ξ²e(ξ)`, i.e. the L^{3,∞}/critical topology), or (ii) locally in ξ near
the active front. A naive Galerkin/L²-overlap pilot will not detect the orbit. This is a
concrete, actionable correction to the S7 recipe.

---

## The NRS/Tsai no-go: does the torus escape, and does the Riemann-sum limit re-import it?

### T.1 **The continuum steady state IS a Leray profile — exactly, not approximately**

Leray ansatz `u(x,t) = λ(t)^{−1}U(x/λ)`, `λ = √(2a(T−t))`. Its Fourier transform is
`û(k,t) = λ²Û(λk)`, and with `N := 1/λ` this is **literally** `û(k,t) = N^{−2}Û(k/N)` —
the S1 ansatz with `Ψ ≡ Û` s-independent, and `Ṅ/N³ = a` constant. Therefore:

> **A steady state (`∂_sΨ = 0`) of the S1 front flow at constant `a` is precisely a backward
> self-similar (Leray) profile of NS. NRS1996 and Tsai1998 apply to it directly.**

### T.2 The critical wake sits exactly at the NRS borderline — and Tsai still kills it

`e_c ∝ ξ^{−2}` ⟹ `∮|Ψ|² ∝ ξ^{−4}` ⟹ `|Û(ξ)| ≍ √c_E ξ^{−2}` ⟹ `U(x) ≍ √c_E |x|^{−1}`.
- `|x|^{−1} ∈ L^{3,∞}(ℝ³) \ L³(ℝ³)`: `∫|x|^{−3}dx` diverges logarithmically at **both** ends.
  So NRS's hypothesis `U ∈ L³` **fails, but only logarithmically** — which is the same
  logarithm as S4's `‖u‖₃³ ≍ log N`. NRS is evaded by a hair. **[C]**
- **Tsai1998 is not evaded.** Its hypothesis is finite *local* energy on `Q₁(0,T)`:
  `ess sup_t ∫_{B₁}|u|² + ∫∫_{Q₁}|∇u|² < ∞`. Compute for `u ≍ 1/|x|` outside the core
  `|x| ≲ λ(t)`: `∫_{B₁}|u|² ≍ ∫_λ^1 r^{−2}r²dr + λ^{−2}λ³ ≍ 1` (bounded), and
  `∫_{B₁}|∇u|² ≍ ∫_λ^1 r^{−4}r²dr ≍ λ^{−1}`, so
  `∫_{T−1}^{T}λ(t)^{−1}dt = ∫(2a(T−t))^{−1/2}dt < ∞`. **Both hypotheses hold ⟹ u ≡ 0.**
  So the exactly-steady front flow with the critical wake is **dead on ℝ³**, notwithstanding
  the L³ log-escape. This is the single most important thing the seed does not say.

### T.3 Why the torus/discrete version escapes — the precise statement

`ℤ³` admits only the discrete dilation subgroup `2^ℤ` (`λN ∈ ℕ`, repo (6.8)). Hence the
S1 flow *on the lattice* cannot have a nontrivial continuously-self-similar solution; the
best available object is an **`s`-periodic orbit with period `S = log2/a₊`**, i.e. a
**discretely self-similar (DSS) solution with ratio 2**, `u(x,t) = 2u(2x, T−4(T−t))`-type.
NRS1996 and Tsai1998 both assume *exact* self-similarity and **do not apply to DSS
solutions**. That, and only that, is the escape. Additionally the torus supplies a hard IR
cutoff `|k| ≥ 1` (and in practice `|k| ≥ N₀`), so the wake is truncated and `‖u(t)‖₃` is
finite at each `t < T`, diverging only logarithmically as `t↑T` — consistent with ESS being
*saturated* rather than violated.

Two riders the seed omits:
- The repo's own `future_search.md` §3.2 already warns: *"Must filter against [Seregin2024]
  Liouville-type restrictions on discrete self-similar solutions — periodicity does NOT
  automatically evade backward-self-similar non-existence theorems."* DSS is a live, not a
  safe, class (Chae–Wolf, Bradshaw–Tsai results on DSS exist and must be checked).
- The escape is only via the *torus/lattice*, so this lane is CLAY-B/D (periodic), **not**
  CLAY-A/TARGET-U. Exporting it to ℝ³ re-imports Z-01/Z-02 (compact-support + exact band
  limitation are incompatible) *and* T.2.

### T.4 **Does the Riemann-sum limit re-import the no-go? — Partly YES, and this is the sharpest falsifiable criterion in the whole audit**

The Riemann sum `Σ_ℓ → N³∫d³λ` is *precisely* the erasure of the lattice. Two facts:

1. **Good news**: the front flow at fixed `(a,ν)` does **not** have continuous dilation
   symmetry — `D_μ` maps `(a,ν)`-solutions to `(a/μ²,ν)`-solutions after rescaling `s`
   (the `ν|ξ|²` term breaks it, fixing the scale `ξ ≍ √(a/ν)`). So an `s`-periodic orbit is
   *not* automatically a symmetry-generated relative equilibrium reducible to a fixed point.
   The escape is not empty on formal grounds.

2. **Bad news**: in the continuum flow the period `S` is not quantized — `log2/a₊` is
   meaningful only because the *lattice* forces the pullback ratio to be an integer power of
   2. But the Riemann-sum approximation is valid exactly when `Ψ` varies on scale `ξ ≍ 1`,
   i.e. when the lattice is invisible to the leading dynamics; the DSS and exactly-SS objects
   then differ only at `O(1/N)` — **unless the orbit's `s`-dependence over one period is
   `O(1)` in the critical topology.** Therefore:

> **Falsifiable criterion (the gate this lane must pass).** The RG pilot must exhibit an
> `s`-periodic orbit whose profile variation over one period,
> `osc := sup_{s∈[0,S]} ‖Ψ(·,s) − Ψ̄‖_crit / ‖Ψ̄‖_crit` (critical/`sup_ξ ξ²e` norm),
> is bounded **below** by an `N`-independent constant. If `osc → 0` as `N` grows, the orbit
> is converging to a steady state of the continuum front flow = a Leray profile, and
> **Tsai1998 kills it** (T.2). "Nearly time-independent renormalized profile" is a *failure*
> signature here, not a success signature — the opposite of the usual RG intuition, and the
> opposite of what a naive reading of S1's "attracting periodic orbit" would celebrate.

Recommend this be pre-registered as a kill condition alongside `ρ_j` decreasing and
`𝒜_j → 0`.

---

## Summary tables

### CONFIRMED
S1 full substitution incl. `N³` Riemann factor, `N^{−3}` exponent chain, sign of `𝒬`, `P_ξ`
scale-invariance; `a = Ṅ/N³ = −½d(N^{−2})/dt`; `N^{−2}=N₀^{−2}−2∫a`; `γ=1/2`; `dN/ds=aN`;
`T−t=∫N^{−2}ds<∞`; `S=log2/a₊`. S2 sphere-average bookkeeping and the coefficient 2 (exact);
`−2νξ²e`; `−∂_ξF` as a definition; **exact** vanishing of the dilation term on `e_c`;
`F′=−4νc_E`; `ξ_max = 1+(χ/√2)√c_E/ν`; `χ√c_E ≥ √2ν`; capacity match `χ ≲ Cη^{3/2}`
(pins γ=1). S3 the c_E-collapse as an *exact* scaling group `(Ψ,s,a,ν)→(μΨ,s/μ,μa,μν)`;
`Σc_E/N_j = 2c_E/N₀`; smoothness unharmed on 𝕋³; three independent routes to the same
threshold. S4 Hölder/Bernstein/Riesz chain; N-independence of `‖u_N‖₃`; wake retention factor
`e^{−ν/(2a₊)}` scale-free. S5 linear polarization ⟹ pointwise zero helicity; Beltrami ⟹
`B(u,u)=0`; `|H_N| ≤ 2c_E`. S6 stage-2 failure of translation-split (symmetry exact);
sumset-miss satisfiable. S7 continuum pullback energy bookkeeping exactly preserving.

### ERRORS (with corrections)
1. **S2.5 (major)** `e ≡ e_c` and `F = χξe^{3/2}` are mutually inconsistent; the "predicted
   measurable profile" is not a solution. Capacity `F(1)ξ^{−2}` lies **below** requirement
   `F(1)(1−β(ξ−1))` for `1<ξ<φ=1.618` at threshold β=1. Fix: `χ√c_E ≥ √2ν` is necessary, not
   sufficient; the true steady problem is `a(2e+ξe′)−2νξ²e−(χξe^{3/2})′=0` with `a` an
   eigenvalue, and the wake must overshoot `e_c`.
2. **S1.3 / S7.4** τ is over-determined: `a₊τ = log2`; with `a₊ ∝ √c_E`, `τ ∝ c_E^{−1/2}`.
3. **S3.2** `a ∝ √c_E` is required for the dilation channel to be `c_E^{3/2}`; the seed
   treats `a` as c_E-independent. Consequences: `T ∝ N₀^{−2}c_E^{−1/2}`, `S ∝ c_E^{−1/2}`.
4. **S3.3** "reduces to shape closure only" — the reduction is `ν_eff = ν/√c_E`, i.e. the
   high-Reynolds limit; it renames the difficulty rather than removing it.
5. **S3.5** `D_N ∝ c_E` is first-Duhamel-only; the repo's `c_E ≈ 228/546/1902` are 250–2000×
   extrapolations of a perturbative formula and are not predictions.
6. **S2.3** front energy is `(2N)^{−1}∫e dξ`, not `N^{−1}∫e dξ` (and `Π_phys = NF/2`);
   harmless because χ's normalization absorbs it, but must be stated.
7. **S4.3/S4.4** `q` is not a shape constant (`q = 2^{3/2}χ√c_E`); the per-octave summation
   `‖u‖₃³ ≳ Σ_j‖u_j‖₃³` is asserted without proof — **repaired here** via the LP square
   function plus `(Σa_j)^{3/2} ≥ Σa_j^{3/2}`.
8. **S5.4** "coherent transfer requires helicity-neutrality" is a non sequitur; correct rule
   is `|h| = |H_N|/(2NE_N)` bounded away from 1 (non-Beltrami), not `h = 0`.
9. **S6.1** suppression exponent: `(W|Δx|)^{−3}` in **amplitude** / `(W|Δx|)^{−6}` in power
   for a box envelope, with an incoherent floor at `W^{−3}` in power; the seed's
   "`(W|Δx|)^{−3}` in power" mixes the two.
10. **S6.3** `Πλ_j = 2` misquotes repo (6.11) (`>1`); equal-ratio L-cycles need `λ=2^{1/L}`,
    irrational hence lattice-incompatible for L≥2; only unequal rational ratios (interleaved
    ladders) work.
11. **S7.1** `2k→k` subsampling loses ~7/8 of the child; the correct map is the repo's (6.8)
    continuum dilation with resampling.
12. **S7.3** `𝒜₊ > 0` is a screen, not the gate; the repo's own gate is child energy ratio
    `≥ 1/2`.

### REFUTED
- **S5.5** the proposed helicity no-go via `|H| ≤ ‖u‖₂‖ω‖₂` fails by an exponential margin
  (`2c_E J` vs `2√2 c_E 2^{J/2}`). A *different* balance (viscous helicity flux
  `≍ (νc_E/a₊)log N` vs demanded `2c_E log₂N`) is marginal at `ν ≈ a₊` — supportive of the
  design rule, but by another route, and not a theorem.
- **S3's** implicit "large c_E is a free lever": it is exactly the small-`ν_eff` limit and is
  simultaneously the Fujita–Kato supercriticality requirement — legitimate but content-free
  as a *reduction*.

### GAPS
G1 `a(s)` has no modulation/gauge condition (nonlinear eigenvalue, unspecified).
G2 No argument that the RG shape direction is contracting (PO-08 unresolved); forward
integration finding an attractor is a hope.
G3 S1 (dense, `M≍N³`, Riemann sum) and S6 (sparse integer carriers) are mutually exclusive
regimes; the reconciliation must go through the repo's per-cell (6.5)–(6.6) form.
G4 S4's band-locality hides the low–high sweeping commutator; valid for exactly band-limited
torus fields, must not be exported to ℝ³ packets (Z-01/Z-02).
G5 S4 is conditional on the sustained **signed** flux lower bound — the repo's own
"single largest unproven step". S4 re-routes obligations 8–9 onto 5; it does not remove them.
G6 The profile is unbounded in the L²-profile topology (`∫e dξ ≍ 2c_E N/N₀`); orbits must be
sought in a critical topology (`sup_ξ ξ²e`) or locally in ξ.
G7 Scale-stagger forces `η ≤ 1/9`–`1/12`, raising required `c_E` by 27–64×; the geometric
cleanliness is bought from the same budget it protects. No accounting given.
G8 **DSS class is not safe**: Seregin-type / Chae–Wolf / Bradshaw–Tsai restrictions on
discretely self-similar solutions are untested against this ansatz.
G9 **The T.4 gate** (`osc` bounded below, N-independently) is not in the seed and must be
pre-registered: an orbit whose s-dependence vanishes in the continuum limit is a Leray
profile and is killed by Tsai1998.
G10 The lane is CLAY-B/D (torus) only; nothing here transfers to ℝ³/TARGET-U.
