# EXT-P1/P2/P3, audited: complete statements, the minimal integral form, repairs, and closure verdicts

**Status (2026-07-31, turn 11; assembler's record).**  Three independent
adversarial audits examined the EXT-P1 paper proof of
[`ext_dependencies.md`](ext_dependencies.md) §§1–3 together with its §5
imports from [`track_p_periodic.md`](track_p_periodic.md), on three lenses:
(A1) compactness, topologies and limit exchange; (A2) inequalities, constants
and normalisation; (A3) statement completeness and what the chain actually
consumes.  Verdicts: A1 "gaps", A2 "sound", A3 "gaps" — eleven raw findings,
deduplicating to **nine distinct items, every one of severity minor**; no
fatal or major defect was found on any audited axis, and every displayed
constant of §3 was independently re-derived by A2 (`270A = 1110.78`,
`135A = 555.39`, `L_k`, `9Y₀/(2ν)`, `16Y₀/K²`, `8A`, `(3.e.1)`, the exact
lattice partial sum `16.27443269…`).  This document does five things: it
fixes the three external statements as complete mathematics, with EXT-P2
replaced by the minimal *integral* form the chain composition needs (§1); it
maps every hypothesis and conclusion onto the repository conventions (§2,
against [`kato_h3_constants.md`](kato_h3_constants.md) §1); it resolves every
audit finding — repairs written out in full, with erratum edits applied to
`ext_dependencies.md` §3 marked `[turn-11 audit erratum]` (§3); it issues
per-theorem closure verdicts (§4); and it records exactly what a
`proved: true` flag would and would not mean (§5).

**The re-audit stage has veto over every verdict below**, and the material
that is *new this turn* (the EXT-P2-INT delivery proof, the EXT-P3 proof, the
`32A` smoothing estimate, the integral comparison lemma) is explicitly
labelled as unaudited new writing with highest re-audit priority.  Wording
contract unchanged and repeated: the only permitted summary sentences remain
"**conditional PDE certificate assuming EXT-P1/P2/P3**" and "**the
finite-dimensional Galerkin enclosure and the scalar control ODE are verified
unconditionally**".  Nothing here is a singularity statement, nothing here
bears on the Clay problem, and **no payload flag is changed by this
document** — flag decisions are the coordinator's alone (§5).

## 0. Notation

Fixed by [`kato_h3_constants.md`](kato_h3_constants.md) §1 and
[`track_p_periodic.md`](track_p_periodic.md) §5.1, restated once:
`T³ = (ℝ/2πℤ)³` with the normalised measure `(2π)⁻³dx`; real, mean-zero
fields `u = Σ_{k≠0} û_k e^{ik·x}`, `û_{−k} = conj(û_k)`; divergence-free
means `k·û_k = 0`; Parseval `⟨f,g⟩ = Σ_k f̂_k·conj(ĝ_k)`; homogeneous norms
`‖f‖_n² = Σ_k |k|^{2n}|f̂_k|²`, monotone in `n` (spectral gap `|k| ≥ 1`);
`Λ^n` is the multiplier `|k|^n`; `P_k = I − kkᵀ/|k|²` the Leray multiplier;
`B_k(v,w) = Σ_{m+n=k}(v̂_m·n)ŵ_n` so `(v·∇w)^̂_k = iB_k(v,w)`;
`A² = Σ_{k≠0}|k|⁻⁴ = A₄`, `4.034 ≤ A ≤ 4.114`, `‖f‖_∞ ≤ A‖f‖₂`;
`H⁴_σ` = mean-zero divergence-free `H⁴(T³;ℝ³)`;
`‖f‖₄ ≤ ‖f‖_{H⁴} ≤ 4‖f‖₄` and `‖f‖_{H³} ≤ √8‖f‖₃ ≤ 2.8285‖f‖₃` on this
class.  Numbered displays `(3.x.y)` refer to `ext_dependencies.md` §3;
`§5.x` refers to `track_p_periodic.md`.

## 1. Complete statements

### 1.1 Definition (strong solution), with both previously missing converses

Let `0 ≤ t₀ < t₁`.  A **strong solution** on `[t₀,t₁]` is a mean-zero
divergence-free `u ∈ C([t₀,t₁]; H⁴_σ) ∩ C¹([t₀,t₁]; H²)` with

```
∂_t u = ν Δu − P(u·∇u)      in H², pointwise in t ∈ [t₀,t₁] .      (Leray form)
```

Three equivalent presentations, each direction now proved:

**(a) Pressure form.**  *Forward* (proved in (3.d), "Pressure via Leray"):
with `p̂₀ = 0`, `p̂_k = i k·(u·∇u)^̂_k/|k|²`, the pair `(u,p)` satisfies
`∂_t u + u·∇u + ∇p = νΔu` classically in space, `p ∈ C([t₀,t₁];H⁴)`.
*Converse* (repairs G7): if `u` is in the regularity class above,
mean-zero and divergence-free, and `∂_t u + u·∇u + ∇p = νΔu` holds for some
distributional `p`, apply `P`: `P ∂_t u = ∂_t u` (the divergence-free
mean-zero fields form a closed subspace containing the `C¹`-curve `u`, so
`∂_t u(t)` lies in it), `P∇p = 0` mode by mode (`P_k(ik p̂_k) = i p̂_k P_k k
= 0`), and `P` is bounded on `H²`; the Leray form follows.  Applying `I − P`
to the same equation then identifies `∇p = −(I−P)(u·∇u)`, i.e. `p` is the
recovered pressure up to the mean normalisation.  Consequence: the
uniqueness clause of EXT-P1★ below covers pressure-form competitors in the
same `u`-regularity class — the class in which the cited classical sources
(Temam Ch. III; Majda–Bertozzi Ch. 3) state their theorems.

**(b) Coefficient form.**  *Forward* (repairs G2; this is the converse the
ledger asserted but never proved): every strong solution satisfies, for every
`k ≠ 0`,

```
û_k(t) = û_k(t₀) − ∫_{t₀}^t [ ν|k|² û_k(s) + i P_k B_k(u,u)(s) ] ds .   (3.d.1)
```

Proof in three lines: `f ↦ f̂_k` is a bounded linear functional on `H²`
(componentwise `|f̂_k| ≤ ‖f‖₀`), so it commutes with the `H²`-valued Bochner
FTC identity `u(t) = u(t₀) + ∫ ∂_t u(s) ds` supplied by `C¹([t₀,t₁];H²)`;
substituting the Leray form pointwise in `s` and using
`(u·∇u)^̂_k = iB_k(u,u)` — an absolutely convergent convolution for
`H⁴ ⊂ C²` fields, by the Cauchy–Schwarz bound `|B_k(v,w)| ≤ ‖v‖₀‖w‖₁`
displayed before (3.d.1) — gives the claim.  *Converse*: (3.d.1) for all `k`
plus `u ∈ C([t₀,t₁];H⁴_σ)` implies `u` is a strong solution — via the steps
of the ledger's (3.d) that are valid at bare `C(H⁴)` regularity: integrand
continuity giving (3.d.2), the bilinear estimate (3.d.3), the compact-range
tail argument, and the `C¹H²` upgrade `[turn-11 re-audit erratum: the
absolute-continuity-of-`‖u‖₄²` step is NOT on this route — it uses the
`L²Ḣ⁵`/`L²Ḣ³` budgets, available only a posteriori via the patching remark
once strongness is established]`; these apply verbatim to any such field,
not only the constructed one.  The ledger §1 sentence "equivalent characterisation" is therefore now
true as stated; the missing direction is the display above.

### 1.2 EXT-P1★ (local existence and uniqueness — final form)

> **Theorem EXT-P1★.**  Fix `ν > 0` and `u₀ ∈ H⁴_σ` (equivalently:
> `û₀,₀ = 0`, `k·û₀ₖ = 0` for all `k`, `Σ|k|⁸|û₀ₖ|² < ∞`).  Set
>
> ```
> T* := 1/(270 A ‖u₀‖₄) ≥ 1/(1110.78 · ‖u₀‖_{H⁴}) ,     A ≤ 4.114 ,
> ```
>
> **depending only on `‖u₀‖₄` and uniform in `ν > 0`** (the viscous term is
> sign-definite and only improves every estimate; this repairs the
> misleading "depending on `ν` and `‖u₀‖_{H⁴}`" of the ledger §1.1 — G4).
> If `u₀ = 0`, read `T* = +∞`: `u ≡ 0` is the (unique) solution and every
> clause below holds trivially `[turn-11 re-audit erratum: the formula
> divides by ‖u₀‖₄; the degenerate case is covered by §3's explicit branch
> and uniqueness, and the patching remark reads the same convention at any
> restart time with ‖u(s)‖₄ = 0]`.
> Then:
>
> *(E) Existence and regularity.*  There is a strong solution `u` on
> `[0,T*]` with `u(0) = u₀`, and it satisfies
> * **time regularity**: `u ∈ C¹([0,T*]; H²)`, equation pointwise in `H²`;
>   additionally `∂_t u ∈ L²(0,T*; Ḣ³)` with the explicit budget (3.d.4);
> * **space regularity**: `u(t) ∈ H⁴_σ ⊂ C²(T³)` for every `t`;
>   `sup_{[0,T*]}‖u(t)‖₄ ≤ 2‖u₀‖₄`; `u ∈ L²(0,T*; Ḣ⁵)` with
>   `2ν∫₀^{T*}‖u‖₅² dt ≤ 9‖u₀‖₄²`;
> * **continuity at the datum**: `u(t) → u₀` *strongly in `H⁴`* as
>   `t ↓ 0` (proved, not assumed: the weak-plus-norm-convergence argument of
>   (3.d) includes `t = 0`);
> * **pressure**: the recovered `p ∈ C([0,T*]; H⁴)`, mean-zero,
>   `‖p(t)‖₄ ≤ 8A‖u(t)‖₄²`, and `(u,p)` solves the pressure form
>   classically in space; conversely any pressure-form solution in the class
>   is of this shape (§1.1(a)).
>
> *(U) Uniqueness clause.*  If `u, v` are strong solutions on any `[t₀,t₁]`
> (same `ν`) with `u(t₀) = v(t₀)`, then `u ≡ v` on `[t₀,t₁]`, with the
> quantitative stability
> `‖u(t)−v(t)‖₀² ≤ ‖u(t₀)−v(t₀)‖₀² · exp(2A(t−t₀)·sup_{[t₀,t]}‖v‖₄)`
> — inequality (3.e.1).  **The uniqueness class is exactly
> `C([t₀,t₁];H⁴_σ) ∩ C¹([t₀,t₁];H²)`** — uniqueness among *strong* solutions
> in this sense.  Not asserted (and not consumed anywhere in the chain):
> uniqueness among Leray weak solutions or mild solutions.
>
> *(S) Time shift.*  The equation is autonomous: `s ↦ u(s−t₀)` is a
> bijection between strong solutions on `[0, t₁−t₀]` and strong solutions on
> `[t₀,t₁]`, preserving data, and **both the existence clause and the
> uniqueness clause hold with any initial time** (repairs G6: the ledger
> recorded the shift for uniqueness only, while the chain and EXT-P3 apply
> the existence clause at `t₀ > 0`).
>
> *(M) Maximal solution.*  All strong solutions with datum `u₀` patch (on
> overlapping intervals, licensed by (U)) to a unique maximal one on
> `[0, T_max)`, `T_max ≥ T* > 0`, which coincides with every strong solution
> with datum `u₀` on that solution's interval.

*Proof.*  `ext_dependencies.md` §3 (Theorem 1 and (3.a)–(3.e)), as repaired
by the `[turn-11 audit erratum]` edits (Lemma 3.a.0 and the (3.d.2)
rewording) and by §1.1 above; audited by A1/A2/A3 with the nine findings of
§3 below as the complete defect list.  ∎

> **Corollary P1-3 (the `n = 3` instantiation).**  The repository's data are
> exact rational finite trigonometric polynomials — `C^∞`, exactly mean-zero
> and exactly divergence-free — so the `H⁴_σ` hypothesis is *free*: EXT-P1★
> applies verbatim to every chain datum `y_n`.  The same solution satisfies
> `u ∈ C([0,T*]; H³_σ)` with `‖u‖_{H³} ≤ √8‖u‖₃ ≤ √8‖u‖₄`, and every `n = 3`
> tube statement of the turn-11 chain concerns *this* solution; the
> uniqueness class remains the `H⁴` class of EXT-P1★.  No separate `H³`
> well-posedness theory is stated, proved, or consumed.

### 1.3 EXT-P2-INT (the minimal integral form that replaces the Dini clause)

The payload EXT-P2 carries a pointwise right-Dini clause (iii) whose
mollification upgrade the ledger deliberately did not reproduce (§1.2 honesty
note there).  The audit (A3, finding G9) confirmed the chain needs strictly
less.  Here is the minimal form, stated exactly, then proved for the EXT-P1★
solution.  **No Dini derivative appears anywhere in this subsection.**

> **EXT-P2-INT (`n ∈ {3,4}`).**
> *Hypotheses.*  `u` a strong solution on a slab `I = [s₀,s₁]` (no relation
> between `u(s₀)` and `u_a(s₀)` is imposed — later slabs restart with
> `δ_n > 0`); `u_a` a Galerkin trajectory on the band `G`: a finite
> trigonometric polynomial in `x`, mean-zero, exactly divergence-free, with
> real-analytic-in-`t` coefficients solving
> `∂_t u_a = νΔu_a − P_G P(u_a·∇u_a)`; `w := u − u_a`;
> `W_n(t) := ‖w(t)‖_n`; `e := (I−P_G)P(u_a·∇u_a)` the exact residual.
> *Conclusions.*
> (i) `w ∈ C(I; H⁴_σ)`; in particular `t ↦ W_n(t)` is **continuous** on `I`
> for `n ≤ 4`.
> (ii) `t ↦ W_n(t)²` is **absolutely continuous** on `I`, and for all
> `s₀ ≤ s ≤ t ≤ s₁` the **integral energy inequality** holds:
>
> ```
> W_n(t)² ≤ W_n(s)² + ∫_s^t 2 [ a_n W_n(τ)² + b_n W_n(τ)³ + ε_n W_n(τ) ] dτ ,
> ```
>
> with any constants dominating the slab-exact quantities, in particular the
> certificate's interval enclosures:
> `n = 4`: `a₄ ≥ sup_I(−ν + 9(K₁+K₂))`, `b₄ ≥ 135A`,
> `ε₄ ≥ sup_I ‖e(τ)‖₄` (§5.4 constants);
> `n = 3`: `a₃ ≥ sup_I(−ν + C_kato(u_a) + C_shift(u_a))`, `b₃ ≥ G3`,
> `ε₃ ≥ sup_I ‖e(τ)‖₃` (`kato_h3_constants.md` §§4–6 constants).

The scalar comparison then runs through the following lemma — **the
`ChainAnalysis` Lean formalisation target** (a statement about two scalar
functions; its formalisation status is tracked by the coordinator and not
asserted here) — and through nothing else.

> **Lemma C (integral comparison; repairs the second half of G8).**  Let
> `W : [s,T] → [0,∞)` be continuous with `W²` absolutely continuous and
>
> ```
> W(t)² ≤ W(σ)² + ∫_σ^t 2 [ aW² + bW³ + εW ](τ) dτ     for all s ≤ σ ≤ t ≤ T ,
> ```
>
> where `a ∈ ℝ`, `b ≥ 0`, `ε ≥ 0`.  Let `R : [s,T] → ℝ` be `C¹` with
> `R' = aR + bR² + ε'` for some `ε' ≥ ε`, and `R(s) ≥ W(s)`.  Then
> `R ≥ 0` and `W ≤ R` on `[s,T]`.

*Proof.*  *`R ≥ 0`*: `R(s) ≥ W(s) ≥ 0`; if `R(t₂) < 0`, let
`τ₂ := sup{t ≤ t₂ : R(t) ≥ 0}`, so `R(τ₂) = 0` and `R < 0` on `(τ₂,t₂]`;
there `R' = aR + bR² + ε' ≥ aR`, so `(R e^{−at})' ≥ 0` and
`R(t₂)e^{−at₂} ≥ R(τ₂)e^{−aτ₂} = 0` — contradiction.

*A.e. differential inequality for `W²`*: `W²` is AC, so `(W²)'` exists a.e.
and integrates back; the two-sided integral hypothesis at pairs `(σ,t)`
shrinking to a point `τ` of differentiability gives
`(W²)'(τ) ≤ 2[aW² + bW³ + εW](τ)` at every such `τ` (the integrand is
continuous, so every point is a Lebesgue point).

*Comparison with the `η`-retreat* (the retreat is what handles the crossing
at `W = R = 0`, where `z ↦ √z` is not Lipschitz — the first-slab regime
`R(s) = W(s) = 0`): suppose `W(t₁) > R(t₁)` for some `t₁`.  Let
`τ := sup{t ∈ [s,t₁] : W(t) ≤ R(t)}` — the set contains `s` and is closed
(`W, R` continuous), so `W(τ) = R(τ)` (equality: `≤` by membership-closure,
`≥` by right-limits) and `W > R ≥ 0` on `(τ, t₁]`.  Fix
`η ∈ (0, t₁−τ)`.  On `[τ+η, t₁]`, `m := min W > 0` (continuous positive
function on a compact set), so `W = √(W²)` is AC there (composition with
`√·`, Lipschitz on `[m², max W²]`), with
`W' = (W²)'/(2W) ≤ aW + bW² + ε ≤ aW + bW² + ε'` a.e.  Then `D := W − R` is
AC on `[τ+η, t₁]` with `D' ≤ aD + b(W+R)D` a.e.; with
`Φ(t) := exp(−∫_{τ+η}^t (a + b(W+R)))` (AC, bounded away from `0` and `∞`
uniformly in `η` by `E := |a|(t₁−τ) + b(t₁−τ)·max_{[τ,t₁]}(W+R) < ∞`),
`(ΦD)' ≤ 0` a.e. with `ΦD` AC, so `ΦD` is nonincreasing and
`D(t₁) ≤ D(τ+η) · e^{E}`.  Let `η ↓ 0`: `D(τ+η) → D(τ) = 0` by continuity,
so `D(t₁) ≤ 0` — contradicting `W(t₁) > R(t₁)`.  ∎

**Patching remark (used twice below).**  Let `u` be a strong solution on a
compact `I = [s₀,s₁]`.  Then `u` enjoys on `I` every property the §3
construction delivers on `[0,T*]`-type intervals, namely
`u ∈ L²(I;Ḣ⁵)`, `∂_t u ∈ L²(I;Ḣ³)`, and absolute continuity of
`t ↦ ‖u(t)‖₄²`.  *Proof*: `μ := max_I ‖u‖₄ < ∞` (continuity, compactness);
for each `s ∈ I`, EXT-P1★ time-shifted at `s` with datum `u(s)` constructs a
solution on `J_s := [s, min(s₁, s + 1/(270Aμ))]` carrying all three
properties ((3.c.4), (3.d.4), and the AC-of-`F` paragraph of (3.d)); the
uniqueness clause (U) identifies it with `u|_{J_s}`.  Finitely many `J_s`
of the uniform length `1/(270Aμ)` cover `I`; each of the three properties is
additive over a finite cover of overlapping closed intervals.  ∎

> **Theorem (delivery: the EXT-P1★ solution satisfies EXT-P2-INT).**  Under
> the EXT-P2-INT hypotheses, conclusions (i) and (ii) hold for `w = u − u_a`,
> for both `n = 4` and `n = 3`.

*Proof.*

*Regularity inventory.*  By the patching remark, on `I`:
`u ∈ C(I;H⁴_σ) ∩ C¹(I;H²) ∩ L²(I;Ḣ⁵)`, `∂_t u ∈ L²(I;Ḣ³)`.  The Galerkin
trajectory `u_a` is a band-limited trigonometric polynomial with
real-analytic coefficients on the compact slab, hence bounded with all time
derivatives in every `Ḣ^m`.  So `w` inherits the full inventory:
`w ∈ C(I;H⁴_σ) ∩ C¹(I;H²) ∩ L²(I;Ḣ⁵)`, `∂_t w ∈ L²(I;Ḣ³)`, and `w` is
mean-zero and divergence-free.  Conclusion (i) is immediate (`‖·‖_n`
continuous on `Ḣ⁴` for `n ≤ 4` by norm monotonicity and the triangle
inequality).

*Coefficientwise error equation.*  Subtracting the Galerkin system from the
coefficient form (3.d.1)/(3.d.2) of `u` (§1.1(b)) and using bilinearity
`B_k(u,u) − B_k(u_a,u_a) = B_k(u_a,w) + B_k(w,u_a) + B_k(w,w)`:

```
ŵ_k' = −ν|k|² ŵ_k − i P_k [ B_k(u_a,w) + B_k(w,u_a) + B_k(w,w) ] − ê_k ,
ê_k  = i (1 − 1_{k∈G}) P_k B_k(u_a,u_a)      (the exact tail of §5.3) ,
```

valid for every `k ≠ 0` and every `t ∈ I` — the coefficientwise §5.3 error
equation.

*(ii), step 1: absolute continuity of `W_n²` (repairs the first half of G8,
clause (ii) of the ledger's honesty note).*  This is the ledger's
AC-of-`F` argument verbatim with `c_k` replaced by `ŵ_k`.  Write
`W_n² = Σ_k F_k`, `F_k := |k|^{2n}|ŵ_k|²`; each `F_k ∈ C¹` with
`F_k' = 2|k|^{2n}Re(ŵ_k'·conj(ŵ_k))` and the majorants

```
n = 4 :  |F_k'| ≤ 2 (|k|⁵|ŵ_k|)(|k|³|ŵ_k'|) =: g_k ,
         Σ_k ∫_I g_k ≤ 2 (∫_I ‖w‖₅²)^{1/2} (∫_I ‖∂_t w‖₃²)^{1/2} < ∞ ,
n = 3 :  |F_k'| ≤ 2 (|k|³|ŵ_k|)(|k|³|ŵ_k'|) =: g_k ,
         Σ_k ∫_I g_k ≤ 2 (∫_I ‖w‖₃²)^{1/2} (∫_I ‖∂_t w‖₃²)^{1/2} < ∞ ,
```

by Cauchy–Schwarz in the product of counting and Lebesgue measure and the
regularity inventory.  Fubini–Tonelli with the summable majorant gives
`W_n²(t) = W_n²(s) + ∫_s^t Σ_k F_k'`, `Σ_k F_k' ∈ L¹(I)`: `W_n²` is AC with

```
(W_n²)'(t) = 2 Σ_k |k|^{2n} Re( ŵ_k'(t) · conj(ŵ_k(t)) )     for a.e. t .
```

*(ii), step 2: the a.e. differential inequality, `n = 4`.*  Let `t` be a
point where the identity above holds, `Σ_k g_k(t) < ∞` (a.e., since
integrable), **and `w(t) ∈ Ḣ⁵`** (a.e., since `w ∈ L²(I;Ḣ⁵)` — this is the
"a.e.-`H⁵` insertion" A3's finding G8 demanded be stated).  Substituting the
error equation, the sum splits into five families, *each absolutely
convergent at such `t`*:

* viscous: `−2νΣ|k|^{10}|ŵ_k|² = −2ν‖w‖₅² ≤ −2ν W₄²` (finite by
  `w(t) ∈ Ḣ⁵`; spectral gap);
* transport `⟨u_a·∇w, w⟩₄`: the §5.4 bound `9K₁‖w‖₄²`.  The `β = 0`
  Leibniz term pairs to zero by the classical cancellation with
  `div u_a = 0` — applied to `f = ∂^α w`, which lies in `H¹` precisely
  because `w(t) ∈ Ḣ⁵`; all other terms carry the exact band majorants
  `M_j`;
* stretching `⟨w·∇u_a, w⟩₄`: the §5.4 bound `9K₂‖w‖₄²` (band majorants
  `N_m`, Frobenius `ℓ¹` factor included);
* cubic `⟨w·∇w, w⟩₄`: the §5.4 bound `135A‖w‖₄³`; its `β = 0` cancellation
  uses `div w = 0` and again `∂^α w ∈ H¹ ⟸ w(t) ∈ Ḣ⁵`;
* residual: `|⟨e, w⟩₄| ≤ ‖e(t)‖₄ W₄` (Cauchy–Schwarz; `e` band-limited and
  exact).

Hence for a.e. `t ∈ I`

```
(W₄²)' ≤ 2 [ (−ν + 9(K₁+K₂)) W₄² + 135A·W₄³ + ‖e(t)‖₄ W₄ ]
       ≤ 2 [ a₄ W₄² + b₄ W₄³ + ε₄ W₄ ] ,
```

and integrating the AC function `W₄²` gives conclusion (ii) at `n = 4`.

*(ii), step 3: the a.e. differential inequality, `n = 3` — the tube runs
under `H⁴` regularity, with no `H⁵` insertion at all.*  At `n = 3` every
pairing closes at the regularity `w(t) ∈ H⁴_σ`, which holds at **every**
`t ∈ I`, not merely a.e.:

* viscous: `−2νΣ|k|⁸|ŵ_k|² = −2ν‖w‖₄² ≤ −2ν W₃²` — finite at every `t`;
* transport: `|⟨P(u_a·∇w), w⟩₃| ≤ C_kato(u_a) W₃²`
  (`kato_h3_constants.md` §5(a)).  The cancellation-first commutator
  argument subtracts `⟨u_a·∇Λ³w, Λ³w⟩ = 0`, which needs
  `Λ³w ∈ H¹ ⟺ w(t) ∈ Ḣ⁴` — available at every `t`; the commutator
  majorants are finite band sums against `‖w‖₃`;
* stretching: `|⟨P(w·∇u_a), w⟩₃| ≤ C_shift(u_a) W₃²`
  (`kato_h3_constants.md` §5(b)); a direct absolutely convergent double sum
  (band `ℓ¹` factor times the exact `(1+|j|)³` shift bound), needing only
  `w(t) ∈ Ḣ³`;
* cubic: `|⟨P(w·∇w), w⟩₃| ≤ G3‖w‖₃³` (`kato_h3_constants.md` §4, with
  `v := w`).  Two convergence checks, both supplied by `H⁴`: the subtracted
  cancellation term `⟨w·∇Λ³w, Λ³w⟩ = 0` needs `Λ³w ∈ H¹ ⟸ w(t) ∈ Ḣ⁴`, and
  its double sum converges absolutely via the `ℓ¹` factor
  `Σ|ŵ_p| ≤ A‖w‖₂` and the `ℓ²` pairing of `|j|⁴|ŵ_j|` (needs `w ∈ Ḣ⁴`)
  with `|k|³|ŵ_k|`; the commutator majorants themselves
  (`ℓ¹ ∗ ℓ²` with `Σ|p||ŵ_p| ≤ A‖w‖₃`) close at `Ḣ³`;
* residual: `|⟨e, w⟩₃| ≤ ‖e(t)‖₃ W₃`.

Hence for a.e. `t ∈ I` (a.e. only because of step 1's derivative, not for
any regularity reason)

```
(W₃²)' ≤ 2 [ (−ν + C_kato(u_a) + C_shift(u_a)) W₃² + G3·W₃³ + ε₃ W₃ ] ,
```

and integration gives (ii) at `n = 3`.  ∎

**Relation to the payload EXT-P2 (finding G9, recorded for the
coordinator).**  The payload wording carries the pointwise right-Dini clause
(iii); the chain's comparison needs only EXT-P2-INT plus Lemma C — the
mismatch is in the safe direction (the chain assumes more than it uses).
The Dini upgrade (classical mollification, Temam III / Majda–Bertozzi §3.2)
remains **unproved in this repository and deliberately so**: it is
superseded, not repaired.  Whether the payloads re-point EXT-P2 at the
minimal form is a coordinator decision; this document only fixes the minimal
form and proves it.

### 1.4 EXT-P3★ (the continuation principle)

> **Theorem EXT-P3★ (`H⁴` continuation criterion).**  Let `u` be the maximal
> strong solution of EXT-P1★ on `[0, T_max)`.  If `T ∈ (0,∞)`,
> `T ≤ T_max`, and `M := sup_{t<T} ‖u(t)‖_{H⁴} < ∞`, then `T_max > T`.
> Equivalently: if `T_max < ∞` then
> `limsup_{t↑T_max} ‖u(t)‖_{H⁴} = ∞`.

*Proof* (repairs G5; the route is: uniform local time → restart →
uniqueness on a **common interval** → Cauchy property → gluing on a
**positive-length overlap**, which eliminates the missing point-concatenation
lemma instead of proving it).

If `T < T_max` there is nothing to prove, so assume `T = T_max < ∞` and
derive a contradiction.  Set `M₄ := sup_{t<T}‖u(t)‖₄ ≤ M < ∞`; if `M₄ = 0`
then `u ≡ 0` extends globally, so assume `M₄ > 0`.

1. *Uniform local time.*  By EXT-P1★ clauses (E)+(S), for every
   `t₀ ∈ [0,T)` there is a strong solution `v` on
   `[t₀, t₀ + T*(t₀)]` with datum `u(t₀)`, where
   `T*(t₀) = 1/(270A‖u(t₀)‖₄) ≥ 1/(270A·M₄) > 0` — a positive lower bound
   **uniform over `t₀ < T`**.
2. *Restart.*  Pick `t₀ ∈ (T − 1/(540A·M₄), T)`.  Then
   `t₁ := t₀ + 1/(270A·M₄) > t₀ + 2(T − t₀) = 2T − t₀ > T`, and `v` is a
   strong solution on `[t₀, t₁]` with `v(t₀) = u(t₀)`.
3. *Identification and the Cauchy property.*  For every `s ∈ (t₀, T)`, both
   `u|_{[t₀,s]}` and `v|_{[t₀,s]}` are strong solutions with the same datum
   at `t₀`; the uniqueness clause (U) — applied on the **common interval**
   `[t₀,s]`, which is all it can do and all it is asked to do — gives
   `u = v` on `[t₀, T)`.  Since `v ∈ C([t₀,t₁]; H⁴)`, it follows that
   `u(t) = v(t) → v(T)` strongly in `H⁴` as `t ↑ T`; in particular
   `(u(t_j))_j` is Cauchy in `H⁴` for every sequence `t_j ↑ T`.
4. *Gluing on an overlap.*  Define `ũ := u` on `[0, t₀]` and `ũ := v` on
   `[t₀, t₁]`.  By step 3, `ũ = u` on all of `[0, s]` for every
   `s ∈ (t₀, T)`, and `ũ = v` on `[t₀, t₁]`; the two closed intervals
   `[0,s]` and `[t₀,t₁]` overlap in the positive-length `[t₀, s]` and cover
   `[0,t₁]`.  Membership of `C(·;H⁴_σ)`, membership of `C¹(·;H²)` and the
   pointwise-in-`t` equation are all local properties, and every point of
   `[0,t₁]` has a relative neighbourhood inside `[0,s]` or inside
   `[t₀,t₁]`, on which `ũ` coincides with the strong solution `u` or `v`
   respectively.  Hence `ũ` is a strong solution on `[0,t₁]` with
   `ũ(0) = u₀` — **no one-sided-derivative matching at an isolated junction
   point is needed anywhere** (the ledger's §1.3 remark glued at the single
   point `t₀`; the audited gap G5 was exactly that unproved step, and this
   proof removes it rather than filling it).
5. *Contradiction.*  By clause (M), `T_max ≥ t₁ > T = T_max`.  ∎

> **Corollary P3-3 (`H³` continuation criterion) [new this turn — not part
> of the audited §3 text; highest re-audit priority].**  Under the same
> hypotheses, if `T ∈ (0,∞)`, `T ≤ T_max` and
> `M₃ := sup_{t<T} ‖u(t)‖₃ < ∞`, then `T_max > T`.

*Proof.*

*Step 1 (a smoothing estimate, `Λ⁴`-commutator form; every ingredient is a
repo-standard step).*  For mean-zero divergence-free `v ∈ H⁵_σ`:

```
|⟨v·∇v, v⟩₄| ≤ 32 A ‖v‖₃ ‖v‖₄² ,        32A ≤ 131.65 .
```

Proof: `⟨v·∇v, v⟩₄ = ⟨Λ⁴(v·∇v), Λ⁴v⟩₀ = ⟨[Λ⁴, v·∇]v, Λ⁴v⟩₀`, since
`⟨v·∇Λ⁴v, Λ⁴v⟩₀ = 0` (classical cancellation with `div v = 0`, applied to
`f = Λ⁴v ∈ H¹ ⟸ v ∈ Ḣ⁵`; the double sum converges absolutely via the `ℓ¹`
factor `Σ|v̂_p| ≤ A‖v‖₂` and the `ℓ²` pairing of `|j|⁵|v̂_j|` with
`|k|⁴|v̂_k|`, i.e. `‖v‖₅‖v‖₄`).  The commutator coefficient is
`([Λ⁴,v·∇]v)^̂_k = Σ_j i(|k|⁴−|j|⁴)(v̂_{k−j}·j) v̂_j`, and with `p = k−j`,

```
| |k|⁴ − |j|⁴ | ≤ 4|p| max(|k|,|j|)³ ≤ 4|p|(|j|+|p|)³ ≤ 16 |p| (|j|³ + |p|³)
```

(MVT on `x⁴`, then the convexity split `(a+b)³ ≤ 4(a³+b³)` already used at
(3.d.3)).  So the coefficient is dominated by
`16 Σ_j [(|p||v̂_p|)(|j|⁴|v̂_j|) + (|p|⁴|v̂_p|)(|j||v̂_j|)]`; Young
`ℓ¹ ∗ ℓ² ⊂ ℓ²` with the §5.1 embeddings `Σ|p||v̂_p| ≤ A‖v‖₃`,
`Σ|j||v̂_j| ≤ A‖v‖₃` gives `‖[Λ⁴,v·∇]v‖₀ ≤ 32A‖v‖₃‖v‖₄`, and Cauchy–Schwarz
against `‖Λ⁴v‖₀ = ‖v‖₄` finishes.

*Step 2 (Grönwall for `Y₄ := ‖u‖₄²` along the solution).*  Fix `τ < T`.  By
the patching remark (§1.3), `F := Y₄` is AC on `[0,τ]` with, for a.e. `t`
(where additionally `u(t) ∈ Ḣ⁵`, a.e. by `L²Ḣ⁵` on each patch),

```
F' = −2ν‖u‖₅² − 2⟨u·∇u, u⟩₄ ≤ 2·32A‖u‖₃·F ≤ 64A·M₃·F ,   64A ≤ 263.30 ,
```

the identity because the a.e. derivative `Σ_k F_k'` rearranges, at such `t`,
into the viscous series (finite by `u(t) ∈ Ḣ⁵`) plus the trilinear pairing
(absolutely convergent: `Σ(|k|³|B_k(u,u)|)(|k|⁵|c_k|) ≤ ‖u·∇u‖₃‖u‖₅ ≤
8A‖u‖₄²‖u‖₅ < ∞` by (3.d.3)), and the inequality by Step 1.  Since `F` is AC
with an a.e. linear differential inequality, `(F e^{−64AM₃t})' ≤ 0` a.e.
with the product AC, so

```
sup_{t ≤ τ} ‖u(t)‖₄² ≤ ‖u₀‖₄² · e^{64A·M₃·τ} ≤ ‖u₀‖₄² · e^{64A·M₃·T} ,
```

a bound **uniform in `τ < T`**.  Hence
`sup_{t<T}‖u(t)‖_{H⁴} ≤ 4·sup_{t<T}‖u(t)‖₄ < ∞`.

*Step 3.*  Apply Theorem EXT-P3★.  ∎

**Consumption note (the two chains).**  The `n = 4` chain consumes EXT-P3★
exactly as in the ledger §1.3 usage note: on a certified slab
`‖u‖₄ ≤ ‖u_a‖₄ + R < ∞`, so `T_max` cannot land inside the slab.  The
turn-11 `n = 3` chain's tube bounds only `‖u − u_a‖₃`; the corresponding
existence clause therefore consumes **Corollary P3-3** (via
`sup_slab ‖u‖₃ ≤ sup‖u_a‖₃ + sup R₃ < ∞`), not the bare `H⁴` criterion —
without P3-3 the `n = 3` lane would have an unnamed gap in its existence
clause.  This dependency is hereby named and discharged at the paper level,
subject to re-audit.

## 2. Correspondence table

Every hypothesis and conclusion of §1 against the repository conventions —
the authoritative conventions table is
[`kato_h3_constants.md`](kato_h3_constants.md) **§1**, row names cited in
the middle column.

| object in §1 statements | repository convention (`kato_h3_constants.md` §1 row) | where computed / enforced |
|---|---|---|
| torus, measure | `T³ = (ℝ/2πℤ)³`, `dμ = (2π)⁻³dx` ("torus", "measure") — no `2π` factors in any norm or pairing | all of `torus_aposteriori` / `torus_chain` rational formulas |
| field data `û_k` | canonical real coefficients `û_k = (a_k − ib_k)/2` ("repo real coefficients"); reality `û_{−k} = conj(û_k)` automatic in the `cos`/`sin` representation ("reality") | `TrigVector` |
| mean-zero hypothesis | zero mode absent, `û₀ = 0` ("zero mode") | construction of the data; `galerkin_modes` excludes `k = 0` |
| divergence-free hypothesis `k·û_k = 0` | "divergence-free" row; exact rational check | `TrigVector.divergence().cleaned()` must vanish (`_check_datum`) |
| `‖·‖_n`, Parseval | "`H^n dot` norm" = `sobolev_sq(n)`; "Parseval" | `sobolev_sq`, `hdot_upper` (turn-11 generalisation of `h4dot_upper`) |
| `H⁴`/`Ḣ⁴`, `H³`/`Ḣ³` conversion | `‖f‖₄ ≤ ‖f‖_{H⁴} ≤ 4‖f‖₄`; `‖f‖_{H³} ≤ √8‖f‖₃` (spectral gap `\|k\| ≥ 1`) | statement conversions in §1.2, §1.4; `track_p_periodic.md` §6 |
| Leray projector `P`, `P_k` | "Leray multiplier" `P_k = I − kkᵀ/\|k\|²`, self-adjoint per mode, `Pw = w` for divergence-free `w` | `leray()`; per-mode Lean brick `TrackPFourier.leray` |
| convolution / advection coefficients `B_k` | "convolution" and "advection" rows (no measure factor); `(v·∇w)^̂_k = iB_k(v,w)` | exact product-to-sum in `torus_aposteriori` |
| Galerkin band `G_N`, band `G` of `u_a` | `G_N = {1 ≤ \|k\|² ≤ N²}` = `galerkin_modes(N²)`; the chain band is the certificate's `G` | `galerkin_modes`; slab configs |
| residual `e` | `e = (I−P_G)P(u_a·∇u_a)`, exact tail on `(G+G)∖G` (§5.3) | `residual_tail`, measured in `Ḣ⁴` or `Ḣ³` (`_slab_constants`, `_slab_constants_h3`) |
| lattice constant `A` | `A² = A₄ = Σ\|p\|⁻⁴` ("lattice sums"): exact partial sum + `26/N` shell tail; `A ≤ 4.114` | `adot_squared_upper` |
| `T* = 1/(270A‖u₀‖₄)` | scalar; uniform in `ν`; **not** the slab length (the slabs are chosen by the Picard box and control ODE — the §3 caution stands) | never conflated with `h` in any config |
| `n = 4` tube constants `a₄, b₄, ε₄` | `a₄ = −ν + 9(K₁+K₂)` (band `ℓ¹` majorants `M_j, N_m`), `b₄ = 135A ≤ 555.4`, `ε₄ = sup‖e‖₄` | `_slab_constants`; `ControlCoefficients` |
| `n = 3` tube constants `a₃, b₃, ε₃` | `a₃ = −ν + C_kato(u_a) + C_shift(u_a)` (exact band sums, §5 of `kato_h3_constants.md`), `b₃ = G3 ≤ 12√A₄`, `ε₃ = sup‖e‖₃` | `_slab_constants_h3`; `kato_constant.py` (`c_kato_upper`, `c_shift_upper`) |
| uniqueness class in (U) | strong solutions `C(·;H⁴_σ) ∩ C¹(·;H²)` — exactly the class the chain's boundary identification uses (`track_p_chain.md` §1); no weak/mild uniqueness consumed anywhere | chain composition step |
| integral comparison (Lemma C) | scalar; no function spaces | the `ChainAnalysis` Lean target (coordinator-tracked) |
| datum class of the chain | exact rational finite trig polynomials, exactly divergence-free after `leray` — `C^∞`, so the `H⁴_σ` hypothesis of EXT-P1★ is free (Corollary P1-3) | `recentre` (Leray-projected exact restart points) |

Numerical anchors re-verified by audit A2 and reused above:
`A² ≤ 16.9245` (partial sum `16.27443269…` at `|k|_∞ ≤ 40` plus `26/40`),
`A ∈ [4.034, 4.114]`, `135A ≤ 555.39`, `270A ≤ 1110.78`, `√8 ≤ 2.8285`;
new this turn: `32A ≤ 131.65`, `64A ≤ 263.30`.

## 3. Line-audit results and resolutions

Deduplication: A1 and A2 found the same two leading items independently.
`Gn` labels are used throughout this note and in the erratum markers.

| id | raw findings | severity | location | issue (compressed) | resolution |
|---|---|---|---|---|---|
| G1 | A1.1 ≡ A2.1 | minor | `ext_dependencies.md` §3(a) (trilinear identity), §3(b) (`β = 0`, `Pu_N = u_N`), §3(c) (constraint inheritance) | invariance of `div u_N = 0` and reality along the Galerkin flow used three times, never proved | **REPAIRED** — Lemma 3.a.0 (linear ODE `s_k' = −ν\|k\|²s_k` with zero datum for the divergence; conjugation symmetry `Θ` + Picard uniqueness for reality) inserted in §3(a) with full proof, plus marked pointers at all three use sites `[turn-11 audit erratum]` |
| G2 | A1.2 ≡ A2.2 | minor | `ext_dependencies.md` §1, Definition, "equivalent characterisation" | only one direction of the strong-solution ⟺ (3.d.1) equivalence proved | **REPAIRED** — the missing direction (strong ⟹ (3.d.1)) proved in §1.1(b) of this note (bounded functional through the Bochner FTC + absolutely convergent convolution identity); the ledger sentence is true as stated with this note as its proof record; no §3 edit needed |
| G3 | A1.3 | minor | `ext_dependencies.md` §3(d), sentence after (3.d.2) | forward reference to `∂_t u` before its `H²`-valued construction (apparent circularity) | **REPAIRED** — sentence reworded to an explicit forward pointer, `[turn-11 audit erratum]`; no mathematics changed |
| G4 | A2.3 | minor | `ext_dependencies.md` §1.1 | "`T` depending only on `ν` and `‖u₀‖_{H⁴}`" — `T*` is in fact `ν`-independent | **REPAIRED** — EXT-P1★ (§1.2 here) states "depending only on `‖u₀‖₄`, uniform in `ν > 0`"; ledger §1 wording superseded by this note (no §3 erratum: the defect is in §1 prose, outside the erratum mandate, and Theorem 1 in §3 is correctly stated) |
| G5 | A3.1 | minor | `ext_dependencies.md` §1.3 remark, §2 row 14 | EXT-P3 restart remark glues two solutions at a single point; the needed concatenation lemma nowhere stated | **REPAIRED** — full proof of EXT-P3★ in §1.4 of this note; the point-gluing is *eliminated* (uniqueness on the common interval `[t₀,s]` identifies `u` with the restart `v`, so the gluing happens on a positive-length overlap where strongness is local); the Cauchy property of `u(t)` as `t ↑ T` falls out of the same identification |
| G6 | A3.2 | minor | `ext_dependencies.md` §1.1 (autonomy sentence) | time-shift invariance recorded for the uniqueness clause only; the chain and EXT-P3 apply the *existence* clause at `t₀ > 0` | **REPAIRED** — EXT-P1★ clause (S) states the shift for both clauses, with the shift map displayed |
| G7 | A3.3 | minor | `ext_dependencies.md` §1, Definition | pressure-form ⟹ Leray-form converse unstated (uniqueness-class identification with the cited sources) | **REPAIRED** — one-line converse proved in §1.1(a) of this note (`P∂_t u = ∂_t u`, `P∇p = 0`), plus the identification of `p` with the recovered pressure |
| G8 | A3.4 | minor | `ext_dependencies.md` §1.2 honesty note (parenthesis) | the claim "(ii) alone suffices for the comparison" left two sub-steps unwritten: (a) a.e.-`H⁵` validity of the §5.4 algebra, (b) the AC-comparison lemma with the `η`-retreat at the `W = R = 0` crossing | **REPAIRED** — both written in full in §1.3 of this note: the delivery theorem states the a.e.-`Ḣ⁵` insertion explicitly at `n = 4` (and shows `n = 3` needs none), and Lemma C proves the comparison with the `η`-retreat |
| G9 | A3.5 | minor | payload EXT-P2 (iii) vs actual consumption | chain conditions on the pointwise-Dini clause while needing only the integral form — a safe-direction over-consumption | **RECORDED, not repaired** — the minimal consumption is now the fixed statement EXT-P2-INT + Lemma C (§1.3); re-pointing the payload wording is the coordinator's call; the Dini upgrade itself stays deliberately unproved (named gap, §4) |

Erratum edits applied to `ext_dependencies.md` (all marked
`[turn-11 audit erratum]`, none touching any payload, checker, or
coordinator-owned file): Lemma 3.a.0 with proof in §3(a); three one-line
pointers at the G1 use sites in §3(a)/(b)/(c); the G3 rewording in §3(d).
Nothing else in the ledger was modified.

## 4. Closure verdicts (per theorem)

Vocabulary: **"closed (audited paper proof)"** is issued only where zero
unproved steps remain after the §3 repairs *and* the underlying text has
been through the three-auditor adversarial audit.  Material written for the
first time in this note is at best **"closed at the paper level, new this
turn"** — complete to this assembler's knowledge, but by definition not yet
independently audited.  **The re-audit stage has veto over every line of
this table**, and closure here moves no flag (§5).

| theorem | verdict | remaining named gaps |
|---|---|---|
| **EXT-P1★** (§1.2; = ledger §3 Theorem 1 + repairs G1–G4, G6, G7) | **closed (audited paper proof)** — the §3 Galerkin proof survived all three adversarial lenses with only the mechanically repairable findings above, and every repair is now written and applied | none |
| **Corollary P1-3** | closed (immediate specialisation; two-line proof in §1.2) | none |
| **EXT-P2** (payload wording, with Dini clause (iii)) | **open** | named gap **G-DINI**: the pointwise right-Dini upgrade (classical mollification, Temam III / Majda–Bertozzi §3.2) — deliberately unproved and superseded by EXT-P2-INT; it is retained only because the payload text carries it |
| **EXT-P2-INT** + Lemma C (§1.3) | **closed at the paper level, new this turn** — zero unproved steps identified in the argument as written; the AC step, the a.e.-`Ḣ⁵` insertion (`n = 4`), the every-`t` `H⁴` closure (`n = 3`) and the `η`-retreat are all displayed; **unaudited: highest re-audit priority** | none identified; re-audit pending |
| **EXT-P3★** (`H⁴`, §1.4) | **closed at the paper level, new this turn** — the overlap-gluing proof removes G5's missing lemma rather than assuming it; **unaudited: re-audit priority** | none identified; re-audit pending |
| **Corollary P3-3** (`H³`, §1.4) | **closed at the paper level, new this turn** — the `32A` smoothing estimate and the patched Grönwall are new mathematics not present in the audited §3 text; **unaudited: highest re-audit priority**, and the `n = 3` chain's existence clause depends on it (consumption note, §1.4) | none identified; re-audit pending |

Anti-bias note, recorded on purpose: the desire to write "closed"
everywhere was a live failure mode for this assembly, so the table
deliberately splits "audited paper proof" (only EXT-P1★ qualifies — the
audited object was the §3 proof, and the repairs are exactly the auditors'
own suggested repairs, executed) from "new this turn" (everything §1.3–§1.4
adds).  If the re-audit stage rejects any new-this-turn argument, the
corresponding verdict reverts to **open** with the rejected step as its
named gap, and nothing downstream may cite it meanwhile.

## 5. What `proved: true` would mean — and who decides

If, after the re-audit stage, the coordinator were to flip any of the
EXT-P1/P2/P3 payload flags to `proved: true` on the strength of this
document, that flag would mean exactly this and no more: **an audited paper
proof exists in this repository — classical analysis, written and
independently audited, and NOT Lean-formalised.**  It would *not* mean
machine-checked; the infinite-dimensional steps (the ledger's
L-EXT-P1-c1 … e1 obligations, and everything in §1.3–§1.4 here) live outside
the repository's verified layer, and the finite-dimensional Lean obligations
of `ext_dependencies.md` §3.6 remain open targets regardless.  **This
document does not flip any flag, and this assembler has no authority to**:
the decision — including whether an audited-paper-proof tier is representable
in the payloads at all, and any closure-annotated wording variant — belongs
to the coordinator alone, with the re-audit stage's veto ahead of it.  Until
then every payload keeps `proved: false` and `axiomatised_in_lean: false`,
`verify_torus_slab_certificate` keeps enforcing both, and the only true
summary sentences remain: the Track-P certificates are a **conditional PDE
certificate assuming EXT-P1/P2/P3**, and **the finite-dimensional Galerkin
enclosure and the scalar control ODE are verified unconditionally**.
Finally, and independent of every verdict above: **Lean axiomatisation of
EXT-P1/P2/P3 remains forbidden** — an audited paper proof is a citation with
provenance, not an axiom, and the three checked reasons of
`ext_dependencies.md` §4 ("Why the axiom-free rule forbids the shortcut")
apply verbatim to the repaired statements.  Nothing in this document is a
singularity statement, and nothing in it bears on the Clay problem.

## 6. Re-audit round and the coordinator's disposition (turn 11, final)

**Re-audit verdicts** (two fresh adversarial re-auditors, each actively
attempting refutation; full transcripts in the session workflow record):

| theorem | re-auditor A | re-auditor B (statement-vs-use) |
|---|---|---|
| EXT-P1★ | **closed** — survives independent re-derivation, zero unproved steps | **closed** — chain consumes (U)/(S)/(M) exactly as delivered |
| Corollary P1-3 | closed | closed |
| EXT-P2 (payload Dini clause) | **not closed** (G-DINI, correctly recorded) | not closed — *and consumed by nothing* |
| EXT-P2-INT + Lemma C | **closed** (paper level; delivery proof and η-retreat verified line by line) | **closed** — the slab constants instantiate the stated sups verbatim |
| EXT-P3★ (`H⁴`) | **closed** (overlap gluing verified) | closed for the `n = 4` lane |
| Corollary P3-3 (`H³`) | **closed** (`32A ≤ 131.65` re-derived from scratch) | closed; consumption by the `n = 3` lane correctly named |

Both re-auditors independently re-verified the `32A` smoothing estimate and
re-confirmed `270A = 1110.78`, `135A = 555.39`.  Re-auditor A additionally
verified the three `kato_h3_constants.md` derivations (`G₃ ≤ 12√A₄`,
`C_kato`, `C_shift`) that the `n = 3` delivery branch imports — closing the
provenance gap it had itself flagged; the §4 EXT-P2-INT row should be read
with this addendum.

**Minor findings of the re-audit round, all applied in place:** (i) the
`u₀ = 0` convention in EXT-P1★ (`[turn-11 re-audit erratum]` in §1.2);
(ii) the converse's proof route in §1.1(b) no longer names the
AC-of-`‖u‖₄²` step (erratum in §1.1); (iii) every `d⁺/dt` display in
`kato_h3_constants.md` §6 and in the chain payload derivation strings is
annotated as shorthand for the delivered integral form (EXT-P2-INT +
Lemma C); (iv) reissued `n = 3` payloads name Corollary P3-3 as the
continuation form they consume.

**Coordinator's flag decision.**  Under §5's protocol and the user's turn-11
instruction (`proved: true` only after completed audit; closure = complete
audited paper proof or faithful reduction, zero unproved steps, gaps named):

* the audited records `EXTERNAL_THEOREMS_AUDITED`
  (`torus_aposteriori.py`) now carry **`proved: true`** for
  EXT-P1/EXT-P2(-INT)/EXT-P3 **with mandatory closure metadata** — method
  `audited_paper_proof`, this document as provenance,
  `lean_formalised: false`, `axiomatised_in_lean: false` — and the open
  G-DINI clause recorded open and unconsumed;
* chain payloads are **reissued** by
  `torus_chain.reissue_chain_certificate`: the conclusion kind becomes the
  checker-enforced sentence "unconditional PDE certificate modulo audited
  classical theorems …", the quantitative content is untouched, and the
  independent checker re-verifies every link of every reissued payload and
  **rejects mixtures** (`proved: true` without closure metadata, closure
  with a conditional kind, a closed Dini clause, a missing P3-3 note on an
  `n = 3` payload);
* v1 payloads remain valid historical artifacts under the v1 wording, and
  the legacy `EXTERNAL_THEOREMS` block is unchanged;
* **Lean axiomatisation remains forbidden**; every Lean theorem in the
  repository remains unconditionally true and none encodes any of this
  analysis.  Nothing in this section is a singularity statement, and
  nothing in it bears on the Clay problem.
