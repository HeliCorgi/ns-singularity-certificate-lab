# Referee report B — adversarial independent re-derivation

**Documents under review.**
Primary: `docs/paper_lambda_dichotomy/complete_proof.md`, `theorem_statement.md`.
Cross-read: `dependency_and_gap_audit.md`, and the source materials
`docs/research_notes/verification_sprint_v1/spectral_front_final_theorem_and_nogo.md`,
`spectral_front_osgood_gate.md`, `lambda_O9_defect_decomposition.md`,
`coherent_family_certificates.md`.

**Method.** Every identity in (i)–(v) of the referee brief was rederived from
scratch in independent notation, with every constant tracked; the band-symmetry
identity was verified by explicit integer enumeration in exact `Fraction`
arithmetic at `N = 2` (and the family laws at `N = 2, 3, 4`); the lattice-sum
and certificate numbers were recomputed independently up to `N = 140`; the
defect decomposition was re-checked on an independently generated random
divergence-free band-limited field. Scripts are throwaway (scratchpad); all
numbers below are reproducible from the formulas quoted.

**Summary of the re-derivation.** *(i)–(iv) reproduce with no discrepancy in any
constant.* I found **no** sign error, no factor of 2, and no arithmetic slip in
(I.1), the closable bound, the square completion, the two-part defect, the
action identity, or the family laws. The defects reported below are of a
different kind: one **circular** positivity argument that everything else rests
on, one **false corollary** for which I exhibit an explicit counterexample, a
**non-exhaustive** dichotomy, and several overclaims and unproved steps.

---

## PART A — line-by-line re-derivation (what checked out)

### A.1 (i) `d/dt log N₀² = 2 Cov/μ − 2ν V/μ` — CONFIRMED, constants exact

From `½Ḣ₀ = T₀ − νH₁`, `½Ḣ₁ = T₁ − νH₂`:

```
d/dt log N₀² = Ḣ₁/H₁ − Ḣ₀/H₀ = 2[T₁/H₁ − T₀/H₀] − 2ν[H₂/H₁ − H₁/H₀].
```

First bracket `= (T₁ − μT₀)/H₁ = (H₀/H₁)·Cov = Cov/μ` with `μ = H₁/H₀` and
`Cov = (T₁ − μT₀)/H₀`. Second bracket `= (H₀H₂ − H₁²)/(H₀H₁)`; expanding
`V = Σp(x−μ)² = H₂/H₀ − μ²` gives `H₀H₂ − H₁² = H₀²V`, so the bracket is
`H₀V/H₁ = V/μ`. Hence exactly the displayed identity. **Matches Lemma 2 term for
term, including the factor 2 on both terms and the placement of `μ` in the
denominator.** Verified numerically (residual `0.0` in binary64) on a random
divergence-free field, together with `T₀ = 0` (`|T₀|/|T₁| = 2.6e−16`).

### A.2 (ii) closable bound and square completion — CONFIRMED, constants exact

Modal Cauchy–Schwarz on `ℂ³`: `|α_k| = |Re⟨û_k, N̂_k⟩| ≤ √(ε_k η_k)`. Vector
Cauchy–Schwarz on `(|x_k−μ|√ε_k)` and `(√η_k)`, using
`Σ(x_k−μ)²ε_k = H₀V` (verified exactly, residual 0):

```
|Cov| ≤ H₀⁻¹ √(H₀V) √(Σηₖ) = √(V‖N‖₂²/H₀).            (= Lemma 3)
```

Square completion: `sup_{w≥0}[aw − νw²] = a²/(4ν)` at `a = √(‖N‖²/H₀)`, so

```
d/dt log N₀² ≤ (2/μ)·‖N‖²/(4νH₀) = ‖N‖²/(2νμH₀) = ‖N‖²/(2νH₁) = KD/(2ν),
```

the last step because `KD = (‖N‖²/D²)·D = ‖N‖²/H₁`. **Every constant agrees with
Lemma 4 (`2/μ`, `1/(4ν)`, `1/(2ν)`).**

Defect: with `S := √(VG/H₀)`, `G := ‖N‖²`,

```
(2/μ)(S − Cov) + (2ν/μ)[√V − √(G/H₀)/(2ν)]²
   = (2/μ)S − 2Cov/μ + 2νV/μ − (2/μ)S + G/(2νμH₀)
   = KD/(2ν) − d/dt log N₀² ,
```

so `Γ^CS + Γ^SC` is exactly the gap and the two radicals cancel in the sum (the
individual terms are irrational; only the sum is rational). **Matches Lemma 4
exactly**; independently re-verified numerically (residual `0.0`, both parts
`≥ 0`).

### A.3 (iii) action representation — CONFIRMED, constants exact

`N = ∂ₜu + νAu` ⇒ `‖N‖² = ‖∂ₜu‖² + 2ν⟨∂ₜu,Au⟩ + ν²‖Au‖²`, and
`2⟨∂ₜu,Au⟩ = Ḣ₁ = Ḋ`. Dividing by `D > 0` and integrating,

```
∫₀^{T'} KD = ∫₀^{T'} ‖∂ₜu‖²/D + ν² ∫₀^{T'} N₁² + ν log(D(T')/D(0)),
```

using `‖Au‖²/D = H₂/H₁ = N₁²` and `‖N‖²/D = KD`. **Identical to Lemma 6,
including the single factor `ν` (not `2ν`) on the logarithm.** The identity
itself is correct; the *corollary* drawn from it is not — see CRITICAL-2.

### A.4 (iv) band symmetry and the family laws — CONFIRMED exactly at N = 2, 3, 4

Explicit integer enumeration of `B₂ = {k ∈ ℤ³ : 1 ≤ |k| ≤ 2}` (32 points),
exact `Fraction` arithmetic, shell by shell:

| shell `|k|²` | # pts | `Σ k_i²` (i = 1,2,3) | `Σ k_ik_j`, `i≠j` | `m·#/3` |
|---|---|---|---|---|
| 1 | 6 | 2, 2, 2 | 0 | 2 |
| 2 | 12 | 8, 8, 8 | 0 | 8 |
| 3 | 8 | 8, 8, 8 | 0 | 8 |
| 4 | 6 | 8, 8, 8 | 0 | 8 |

`Σ_{B₂} k_ik_j f(|k|²) = δ_ij·(1/3)Σ_{B₂}|k|²f(|k|²)` verified exactly for
`f = 1, |k|⁻², |k|⁻⁴, |k|⁻⁶` **and for the sign-changing weight `f = |k|²−3`**
(so the hypothesis `f ≥ 0` in Lemma 9 is superfluous — MINOR-1). Sample values:
`Σ|k|²f = 32` (`f=|k|⁻²`), diagonal `32/3` each; `Σ|k|²f = 97/6` (`f=|k|⁻⁴`),
diagonal `97/18` each; `Σ|k|²f = 739/72` (`f=|k|⁻⁶`), diagonal `739/216` each.

Family laws, exact, `v₀ = (1,2,3)`, `‖v₀‖² = 14`, `(2/3)‖v₀‖² = 28/3`:

| N | `S_N` (exact) | `T_N` (exact) | `H₀ − (28/3)T_N` | `H₁ − (28/3)S_N` | `u_N(0) − (2/3)S_N v₀` | `H₁/H₀ = S_N/T_N` |
|---|---|---|---|---|---|---|
| 2 | `97/6` | `739/72` | **0** | **0** | **(0,0,0)** | exact equality |
| 3 | `149/5` | `134443/10800` | **0** | **0** | **(0,0,0)** | exact equality |
| 4 | `4888669/120120` | `231195938461/17314577280` | **0** | **0** | **(0,0,0)** | exact equality |

My `S₄`, `T₄` reproduce `coherent_family_certificates.md` Table 1 exactly.
I also confirm `|u_N(0)|²/H₁ = (2/3)S_N` **as an equality** at `N = 2,3,4`, so
Lemma 10's `‖u_N‖∞²/H₁ ≥ (2/3)S_N` is correct and sharp at `x = 0`.
Derivations reproduce the paper's: `|P_kv₀|² = ‖v₀‖² − (k·v₀)²/|k|²`,
`Σ(k·v₀)²|k|⁻⁶ = ⅓‖v₀‖²T_N`, `Σ(k·v₀)²|k|⁻⁴ = ⅓‖v₀‖²S_N`,
`Σ k(k·v₀)|k|⁻⁴ = ⅓ S_N v₀`. **No constant discrepancy.**

### A.5 (v) exponential-class conclusion — arithmetic CONFIRMED, covering step NOT (MAJOR-4)

`K(u_N) = ‖N‖²/H₁² ≥ c₀N³/((28/3)S_N)² ≥ c₀′N` (uses `S_N ≤ 4πN + C₀`), and
`N₀² = S_N/T_N ≤ (4πN + C₀)/T₁ = O(N)` (uses `T_N ≥ T₁ = 6`), so
`K(u_N) ≥ c₀″N₀² = c₀″e^{s_N}`; then `Φ(s_N) ≥ c₀″e^{s_N}` and
`Φ(s) ≥ ½c₀″e^{s}` once `s_{N+1} − s_N ≤ log 2`, giving `∫^∞ ds/Φ < ∞`.
**The arithmetic is right**; the *covering* of `[s_{N₀},∞)` by `[s_N, s_{N+1}]`
is asserted, not proved (MAJOR-4).

### A.6 Numerical claims re-checked against repository outputs

Independent recomputation (`|k| ≤ 140`, exact-count shell sums):

| claim (location) | independent value | verdict |
|---|---|---|
| `4π(N−2) ≤ S_N`, `N ≥ 2` (L11) | `min_{2≤N≤140}(S_N − 4π(N−2)) = +15.565` | **holds** |
| `S_N ≤ 4πN + C₀` (L11) | `max_{2≤N≤140}(S_N − 4πN) = −7.899` → `C₀ = 0` suffices | holds, but `C₀` never exhibited (MAJOR-5) |
| `T_∞ ≤ 16.6` (thm stmt, Constant dependence) | `T_140 = 16.44255`, tail `≤ 4π/139.5`, so `T_∞ ≤ 16.5326` | **holds** |
| `N₀² = S_N/T_N ≍ N` | `N₀²/N = 0.76199, 0.76125, 0.76147, 0.76099` at `N = 4,8,16,32` | holds |
| `N₀²` at `N = 4,6,8,16,32` | `3.0479, 4.6387, 6.0900, 12.1835, 24.3516` | matches `osgood_gate` §2 table and certificates Table 1b |
| `s_{N+1} − s_N → 0` | strictly decreasing `0.419 → 0.051` over `N = 2…20`, `0.0166` at `N = 60` | holds |
| `s_N` monotone (needed, not claimed) | strictly increasing for every consecutive `N ∈ [2,60]` | true, unstated (MAJOR-4) |
| `K(u₄) = 0.7884107043…`, `K(u₆) = 1.4344718…`, `K(u₈) = 2.0372430…` | certificates Table 2: `0.7884107043392768`, `1.4344718079662586`, `2.0372430105569106` | **match** |
| float-pipeline agreement `≤ 1.3e−15` | certificates Table 3: `0.0`, `4.64e−16`, **`1.31e−15`** | `1.31 > 1.30`; write `≤ 1.4e−15` (MINOR-2) |
| sweeping bound captures `91.2 %–92.1 %` for `4 ≤ N ≤ 32` | `0.92104, 0.91693, 0.91509, 0.91265, 0.91176` at `N = 4,6,8,16,32` | range is `91.18–92.10 %`, monotonically **decreasing**, sampled at 5 points only (MINOR-3) |
| `K/N₀²` increases `0.259 → 0.396` | `0.2586696704 → 0.3961869815` | **match** |

Free-standing check of the audit's reduction (§3): with
`‖(v₀·∇)u‖₂² = Σ(v₀·k)²|P_kv₀|²|k|⁻⁴ = ⅓‖v₀‖⁴S_N − Σ(v₀·k)⁴|k|⁻⁶ ≤ ⅓‖v₀‖⁴S_N`,
so the stated bound `≤ ⅔‖v₀‖⁴S_N` is correct but loose by a factor 2 (at `N = 4`
the true value `1089.27` is `0.41×` even the `⅓` bound). Harmless.

---

## CRITICAL

### CRITICAL-1 — the proof that `H₀ > 0` is circular; **forward** uniqueness is the wrong tool

*Location.* `complete_proof.md` Lemma 0, lines 43–49 ("Positivity: … hence
everything"). Inherited from `spectral_front_final_theorem_and_nogo.md` §0 and
`spectral_front_osgood_gate.md` §1 ("backward uniqueness … is not needed"), and
asserted as settled in `dependency_and_gap_audit.md` §4.2.

*Defect.* Two separate problems.

1. *Wrong direction.* "if `H₀(t₁) = 0` then `u(t₁) = 0` and by forward uniqueness
   (F1) `u ≡ 0` on `[t₁,∞)`" is true but vacuous: it constrains the solution
   **after** `t₁` and says nothing about `[0,t₁)`. To rule out `H₀(t₁) = 0` for
   `u₀ ≠ 0` one needs **backward** uniqueness. The audit's boast (§2) that the
   Main Theorem "does not use … backward uniqueness" is therefore either false or
   the proof is incomplete.
2. *Circularity.* The offered substitute — "`H₀(t) ≥ H₀(0)exp(−2ν∫₀ᵗN₀²)` … the
   set `{H₀>0}` is open and closed" — is exactly circular. From
   `Ḣ₀ = −2νH₁ = −2νN₀²H₀` one gets the **identity**
   `log H₀(T) = log H₀(0) − 2ν∫₀^T N₀²`, i.e. `H₀(T) = 0 ⟺ ∫₀^T N₀² = ∞`. The
   quoted Grönwall inequality is that identity; it cannot exclude its own
   failure case. Continuity of `N₀²` **on the open set** `{H₀>0}` gives no
   integrability up to its boundary, so `{H₀>0}` is *not* shown closed: `H₀` is
   non-increasing, so `H₀(t_n)>0` for `t_n ↑ t*` with `H₀(t*) = 0` is entirely
   consistent with every fact used.

*Why CRITICAL.* `H₀ > 0` (hence `H₁ > 0`, `μ > 0`, `D > 0`) is presupposed by
Lemmas 1–8 without exception: it is the denominator in `p(k)`, `Cov`, `V`,
`log N₀²` and in the division by `D` in Lemma 6.

*Minimal fix (avoids backward uniqueness, but forces a re-ordering).* Prove
Lemma 7 (`KD ≤ ‖u‖_∞²`) first — it needs no positivity, only
`‖ℙ(u·∇u)‖₂ ≤ ‖u‖_∞‖∇u‖₂` — then bootstrap: fix `T' < T_max`, put
`M := sup_{[0,T']}‖u‖_∞ < ∞` (finite since `H^m ↪ L^∞` for `m > 5/2`) and
`t_* := sup{t ≤ T' : H₀ > 0 on [0,t]}`. On `[0,t_*)` Lemmas 0–4 apply, so
`log N₀²(t) ≤ log N₀²(0) + M²t/(2ν)`; hence `∫₀^{t_*}N₀² < ∞` and
`H₀(t_*) = H₀(0)e^{−2ν∫₀^{t_*}N₀²} > 0`, so `t_* = T'`. (Alternatively cite
backward uniqueness and amend the audit.) **Either repair breaks
`dependency_and_gap_audit.md` §1's claim "no lemma uses a later one".**

### CRITICAL-2 — the "In particular" clause of Main Theorem (c) is FALSE; explicit counterexample

*Location.* `theorem_statement.md` (c), lines 73–76: "In particular the
regularity criterion of (b1) is equivalent (given the energy equality) to the
a-priori finiteness of the `Ḣ¹`-bandwidth action
`∫(‖∂ₜu‖₂²/D + ν²N₁²)dt`." Same claim in
`spectral_front_final_theorem_and_nogo.md` §4 ("Reading (fixed)") — there hedged
as "given `log D` bounded, which follows from either side plus energy", which is
also false — and as the "free corollary" of `spectral_front_osgood_gate.md` §3.

*Counterexample (exact, elementary).* Take `u₀ = (cos x₂, 0, 0)` on
`𝕋³ = (ℝ/2πℤ)³`: zero mean, divergence-free (`∂₁u₁ = 0`), `C^∞`, `u₀ ≠ 0`, and
`(u₀·∇)u₀ = u₁∂₁u₀ = 0`. Hence `𝒩 ≡ 0` and `u(t) = e^{−νt}u₀` is the (unique,
global) strong solution, with `Δu₀ = −u₀`. Then, normalising `H₀(0) = 1`:

```
K ≡ 0,        ∫₀^∞ KD dt = 0 < ∞                        → case (b1) holds
D(t) = e^{−2νt},  ‖∂ₜu‖₂² = ν²e^{−2νt},  N₁² = H₂/H₁ ≡ 1
∫₀^{T'}‖∂ₜu‖²/D = ν²T',   ν²∫₀^{T'}N₁² = ν²T',   ν log(D(T')/D(0)) = −2ν²T'
```

Lemma 6 checks exactly (`0 = ν²T' + ν²T' − 2ν²T'` ✓), but the
`Ḣ¹`-bandwidth action is `2ν²T' → ∞`. **So (b1) holds while
`∫(‖∂ₜu‖²/D + ν²N₁²)dt = +∞`: the claimed equivalence fails, and it fails for
*every* decaying global solution**, because `log D(T') → −∞` at rate `−2νT'`
while the paper's parenthetical "(given the energy equality)" only bounds `D`
from **above**.

The converse direction is also unjustified: if the action is finite but
`D(T') → ∞`, then `∫KD = A + ν log(D/D₀) → ∞`, i.e. case (b2). So neither
implication is proved.

*Minimal fix.* Delete "equivalent" and state the honest content: **(c) is an
identity**, and `∫KD < ∞ ⟺ ∫(‖∂ₜu‖²/D + ν²N₁²)dt < ∞` **holds only on the set
where `log D` is bounded above and below**, which is not implied by the energy
equality (it fails as `t → ∞` for every solution decaying to rest). Restrict the
corollary to `T_max < ∞` *and* add a proof that `inf_{[0,T_max)}D > 0` there, or
drop the corollary. Correct the same sentence in the two source notes.

---

## MAJOR

### MAJOR-1 — the dichotomy (b) is **not exhaustive**; "Exactly one" is unproved

`theorem_statement.md` (b), lines 60–64 (and `…final_theorem_and_nogo.md` §2 T2,
`complete_proof.md` Lemma 5). The two alternatives are mutually exclusive (the
proof establishes only this), but the case

```
T_max = ∞  and  ∫₀^∞ KD dt = ∞
```

satisfies neither (1) nor (2) and is nowhere excluded. Lemma 5's last sentence
("Exclusivity: if `T_max < ∞`, the same argument rules out `∫KD < ∞`") proves
exclusivity, not exhaustiveness. *Minimal fix:* state (b) as
"`∫₀^{T_max}KD dt < ∞ ⟹ T_max = ∞` and `u ∈ L^∞H¹`; equivalently
`T_max < ∞ ⟹ ∫₀^{T_max}KD dt = ∞`", and replace "Exactly one of the following
holds" by "At most one"; or prove that every global strong solution on `𝕋³` has
finite action (plausible via eventual exponential decay, but that is a further
argument, not in the document).

### MAJOR-2 — the Proposition is quoted as refuting clause (e), but (e) contains an `R` that the Proposition never touches

*Locations.* `theorem_statement.md` (e), lines 94–96: "*(This clause is a proven
implication; the following proposition shows its hypothesis cannot hold
uniformly.)*"; and the Proposition's conclusion, lines 120–123: "clause (e) of
the Main Theorem cannot be activated by **any field-inequality route**"; echoed
in `complete_proof.md` and `dependency_and_gap_audit.md` §5.

*Defect.* Clause (e)'s hypothesis is `KD ≤ Φ(z)D + R` with `R ≥ 0`,
`∫R dt < ∞`. The Proposition assumes and refutes only `K(u) ≤ Φ(log N₀²(u))`,
i.e. the `R ≡ 0` case, and only as a **static, all-fields** inequality. It says
nothing about (a) a solution-dependent `R` with a uniform budget, (b) an
inequality restricted to fields that actually occur along a Navier–Stokes
trajectory, or (c) inequalities of the form `K ≤ Φ(z) + Ψ(V/μ²)` (which the
`lambda_O9` note explicitly keeps alive as conjecture O-8′). The source
`spectral_front_osgood_gate.md` is *scrupulous* here (§1: "The only meaningful
gate is the **uniform** bridge"; §5: "**Scope limitation, stated plainly** … §5
establishes the *persistence mechanism* … not a dynamic realisation of the
excess itself"). **The compression drops this caveat and the parenthetical in
(e) overstates the source.**

*Minimal fix.* In (e) write "…shows its hypothesis **with `R ≡ 0`** cannot hold
as a uniform field inequality"; in the Proposition's conclusion replace "any
field-inequality route" by "any uniform pointwise field inequality
`K ≤ Φ(log N₀²)` valid on all zero-mean divergence-free trigonometric fields",
and add one sentence recording that a trajectory-restricted or `R`-assisted
bridge is untouched.

### MAJOR-3 — "(equality cases characterised)" is asserted in the Theorem and never proved

`theorem_statement.md` (a), line 58. `complete_proof.md` Lemma 4 proves only
`Γ^CS ≥ 0` and `Γ^SC ≥ 0`; there is no equality analysis anywhere in the
proof document. The source *does* have it
(`…final_theorem_and_nogo.md` (I.4): `Γ^CS = 0` iff modal alignment
`a_k² = e_kn_k` plus `a_k ∝ (x_k−μ)x_k^r e_k`; `Γ^SC = 0` iff
`V_r = G_r/(4ν²H_r)`; and `lambda_O9…md` §7 shows the combination is the
nonlinear eigen-equation `ℙ((u·∇)u) = 2ν(Δ+μ)u`). Also note that with the
**full** `G = ‖𝒩‖₂²` used in the paper (rather than the source's `G^in`),
`Γ^CS = 0` additionally forces **zero leakage** `η_k = 0` off `supp(u)`
(`lambda_O9…md` §7(2)) — a condition the source's statement does not carry.
*Minimal fix:* import the characterisation (with the leakage clause) into
Lemma 4, or delete the parenthetical from (a).

### MAJOR-4 — the Proposition's covering step presumes `s_N` is monotone

`complete_proof.md` lines 233–237: "for `s ∈ [s_N, s_{N+1}]`, monotonicity and
`s_{N+1} − s_N → 0` give …". Lemma 11 proves only `s_{N+1} − s_N → 0`; nothing
establishes `s_N ≤ s_{N+1}`, and without it `⋃_N[s_N,s_{N+1}]` need not cover a
half-line, so `Φ(s) ≥ ½c₀″e^s` "for all large `s`" does not follow. Same gap in
`…final_theorem_and_nogo.md` §7. (I verified numerically that `s_N` **is**
strictly increasing for all consecutive `2 ≤ N ≤ 60`, so the statement is true.)
*Minimal fix — two lines, and it avoids needing monotonicity at all:* since
`s_N → ∞` and `s_{N+1} − s_N → 0`, for `s` large let `N` be the **first** index
with `s_N > s`; then `s_{N−1} ≤ s` and `s − s_{N−1} ≤ s_N − s_{N−1} → 0`, so
`Φ(s) ≥ Φ(s_{N−1}) ≥ c₀″e^{s_{N−1}} ≥ ½c₀″e^{s}`. (Alternatively prove
monotonicity from `ΔT_N ≤ N⁻²ΔS_N`, `T_N ≥ 6`, `S_N ≤ 4πN`.)

### MAJOR-5 — Lemma 11 is a sketch containing literal ellipses; "all constants are explicit" is false

`complete_proof.md` lines 217–227 contains, verbatim, unfinished formulas:
`|k|⁻²(1 + (√3/2)·(2/|k|)·…)⁻¹` and `∫…dx = 4π(N ∓ …)`. This is not a proof of
`4π(N−2) ≤ S_N ≤ 4πN + C₀`, and `C₀` — declared "explicit" in Lemma 11 and again
in `theorem_statement.md` §"Constant dependence" — is never exhibited. (Numerically
`S_N − 4πN ≤ −7.9` for `2 ≤ N ≤ 140`, so `C₀ = 0` works and the parenthetical
"(first-shell corrections)" is unnecessary.) Also `N₀² ∈ [c₁N, c₂N]` fails at
`N = 2` (the lower bound `4π(N−2)/T_∞` is `0`); the claim needs "for `N ≥ N₁`"
with `N₁` named. *Minimal fix:* write out the cube comparison
(`|x| ≤ |k| + √3/2` on the unit cube at `k`), state `C₀ = 0` with the verified
range, and restrict the ratio bound to `N ≥ 4`.

### MAJOR-6 — audit claims that are refuted by CRITICAL-1/2

`dependency_and_gap_audit.md` §4.2 ("`H₀ = 0` … proven impossible … via the
Grönwall two-sidedness of `H₀` and forward uniqueness"), §5 ("There are no other
gaps: the Main Theorem's proof chain … is complete"), §1 ("no lemma uses a later
one"), and §2 ("the Main Theorem does not use … backward uniqueness"). At least
one of these must be withdrawn once CRITICAL-1 is repaired. §5's gap table must
also record the (c)-corollary defect. *Minimal fix:* rewrite §4.2 with the
bootstrap of CRITICAL-1, redraw the graph with the `L7 → L0` edge, and add a row
to §5.

### MAJOR-7 — two mutually inconsistent ledger conventions in the cited sources (paper is correct; provenance is not)

`spectral_front_osgood_gate.md` preamble declares `D = ‖∇u‖₂² = 2H₁` and
`K = G₀^full/(2H₁²)`, and its §2 T2 proof asserts `G₀^full = ‖𝒩‖₂²`; but
`lambda_O9_defect_decomposition.md` §1 declares `H_r = Σ|k|^{2r}|û_k|²`,
`G_r = Σ|k|^{2r}|N̂_k|²` (full lattice) while its §4 declares
`G₀ = ½‖ℙ(u·∇u)‖₂²`. These cannot all hold: the pair differs from the complex
full-lattice pair by exactly a factor 2 in **both** `H_r` and `G_r`.
`coherent_family_certificates.md` uses a third statement (`H₁ = ‖∇u‖₂² = D`),
which is the one `complete_proof.md` adopts and which my exact recomputation of
`H₁ = (28/3)S_N` confirms. **The paper is self-consistent and correct**; the
hazard is that anyone tracing a number back through the source notes will be off
by 2. *Minimal fix:* add one convention-map line to §"Conventions" of
`complete_proof.md` (`repo H_r^paired = ½H_r^paper`, `G₀^paired = ½‖𝒩‖₂²`) and
correct the two source preambles.

---

## MINOR

1. **Lemma 9**: hypothesis "`f ≥ 0` radial" — nonnegativity is unused and
   unnecessary (verified with a sign-changing `f`). Say "`f` radial" (finite
   sums, so no convergence hypothesis needed either).
2. **`≤ 1.3 × 10⁻¹⁵`** (complete_proof line 244; audit §3) — the actual maximum
   in `coherent_family_certificates.md` Table 3 is `1.31e−15`. Write `≤ 1.4e−15`
   or `≈ 1.3e−15`.
3. **"near-constant 91.2 %–92.1 % … for `4 ≤ N ≤ 32`"** (complete_proof line
   249) — the measured values are `0.92104, 0.91693, 0.91509, 0.91265, 0.91176`,
   i.e. `91.18–92.10 %`, **monotonically decreasing**, at 5 sampled `N`. Write
   "`≥ 91.1 %` at the five tested bands `N ∈ {4,6,8,16,32}`, decreasing slowly".
4. **`g_k = α_k/ε_k` is `0/0` on empty modes** (theorem statement (a); source
   T1). `p(k) = 0` there so `Cov_p(x,g)` is well defined under the convention
   `g_k := 0` when `ε_k = 0`, but the convention should be stated. Also
   `Cov_p(x,g) = Σp(x−μ)g` only because `E_p[g] = T₀/H₀ = 0` (Lemma 1); the
   theorem statement uses this before Lemma 1 is available.
5. **`Ω ∈ C¹` is false for merely nondecreasing `Φ`** (Lemma 8, line 167;
   `osgood_gate` §1). A monotone `Φ` may be discontinuous. Fix: either assume
   `Φ` continuous, or note `Ω` is locally Lipschitz and increasing, `Ω∘z` is
   absolutely continuous, and the chain-rule inequality holds a.e. — which is
   all that is used.
6. **Lemma 8 does not literally deliver (b1)**: it gives `z ≤ z_*`, hence
   `T_max = ∞` and `u ∈ L^∞H¹`, but "(b1)" as stated in the theorem *is*
   `∫KD < ∞`. Add the one line
   `∫₀^∞KD ≤ Φ(z_*)∫₀^∞D + ∫₀^∞R ≤ Φ(z_*)H₀(0)/(2ν) + ∫R < ∞`.
7. **`(2/q + 3/p = 1/2 < 1)` at `(p,q) = (6,∞)`** — arithmetic correct, but the
   text should say `q = ∞`, `p = 6` explicitly (as written, `(p,q) = (6,∞)` with
   `2/q + 3/p` is fine; `theorem_statement`'s "Serrin `(∞,2)` action" in (d) uses
   the *opposite* `(q,p)` ordering in the same document — unify the convention).
8. **`C_S` is the homogeneous constant.** `theorem_statement.md` "Constant
   dependence" and audit F3 say "`H¹(𝕋³) ↪ L⁶`"; Lemma 7 needs
   `‖u‖_{L⁶} ≤ C_S‖∇u‖_{L²}` — the **`Ḣ¹`** (mean-zero) inequality. Write `Ḣ¹`.
9. **Mean-zero preservation is never stated.** The setting *assumes* the
   `k ≠ 0` expansion; that `∫u dx` is conserved (hence `= 0`) under
   `∂ₜu = νΔu − ℙ(u·∇u)` should be one line, since every use of `C_S` depends on
   it. (See "external inputs" below — it *is* available.)
10. **Lemma 5's parenthetical "(mean-zero Sobolev embedding)" is misapplied.**
    At L5 the *inhomogeneous* `H¹ ↪ L⁶` suffices (`H₀ ≤ H₀(0)` bounds the full
    `H¹` norm), and it needs no mean-zero. Mean-zero is needed only at L7.
11. **"the two external classical inputs F1, F2"** (complete_proof line 5) versus
    the audit's **three** (F1, F2, F3). Also the audit's dependency graph routes
    F3 only into L7, while L5 also invokes Sobolev.
12. **"(the solution is global and smooth for `t > 0`)"** (theorem statement
    (b1)) consumes parabolic instantaneous smoothing `H^m → C^∞`, which is a
    fourth external input not listed in the audit's "No other external result is
    consumed".
13. **`‖∇u‖_{L³}` is called "the critical vorticity-class action"** (d). The
    scaling is right (`2/2 + 3/3 = 2`), but converting `‖∇u‖_{L³}` to `‖ω‖_{L³}`
    needs Calderón–Zygmund (an unlisted input). Either say "velocity-gradient
    class" or list the CZ constant.
14. **`ℙ` self-adjointness** is used implicitly in Lemma 1
    (`⟨u, ℙ(u·∇u)⟩ = ⟨ℙu, u·∇u⟩`); the text cites only `ℙu = u`.
15. **Lemma 0's proof contains typographic debris**: "`(r ≤ 1 ≤ m−… indeed …)`"
    (line 39). Clean up.
16. **`v₀ ∈ ℤ³∖{0}`** (Proposition) — integrality is irrelevant to every stated
    conclusion (only to exact-`Fraction` certification); say so, or drop to
    `v₀ ∈ ℝ³∖{0}`. Note also that (L\*) is certified only for `v₀ = (1,2,3)`
    while the Proposition quantifies over all `v₀`.
17. **The sweeping-pairing bound needs `(v₀·∇)u ≢ 0`** (complete_proof lines
    245–249) — true for the family (some occupied `k` has `v₀·k ≠ 0`) but the
    division is unremarked.
18. **`Q(v₀)/N²` is not `≍`.** The audit §3 calls `Q(v₀)` "`≍ N²`-sized
    numerically"; the measured `Q/N²` is `640.8, 804.3, 850.7, 969.2, 1029.5`
    over `N = 4…32` — still increasing by `60 %`, with no plateau. This only
    *helps* (L\*), so nothing is wrong, but "`≍`" claims a two-sided bound the
    data do not show. Similarly `‖ℙ(u·∇u)‖₂/N^{3/2}` rises `42.2 → 62.9` without
    stabilising.
19. **Scale non-invariance of the Proposition's hypothesis.** `K` has dimensions
    of length², `N₀²` of length⁻², so `K ≤ Φ(log N₀²)` is not invariant under
    `u ↦ λu(λ·)`; on `𝕋³` only integer `λ` is available, so nothing breaks, but
    the statement is intrinsically a torus statement and should say so.
20. **`u₀ ≠ 0` and `T_max` bookkeeping.** In Lemma 5, the hypothesis
    `∫₀^{T_max}KD < ∞` refers to the *pre-extension* `T_max`; the conclusion
    `u ∈ L^∞(0,T_max;H¹)` after re-defining `T_max = ∞` needs the standard
    continuation bootstrap, which is elided.

---

## Stress test of the two external inputs (as requested)

**(1) `L^∞L⁶ ⟹` regularity — attribution and use.** *Attribution is correct.*
`(p,q) = (6,∞)` gives `3/p + 2/q = 1/2 < 1`, strictly **sub**critical, so this
is the classical Prodi–Serrin/Ladyzhenskaya regime with strict inequality; it
does **not** need the ESS endpoint `L^∞L³`, and the audit correctly says so.
*Use is valid but heavier than necessary and slightly mis-motivated.* What
Lemma 5 actually produces is a uniform bound on the **full** `H¹` norm
(`H₁ ≤ e^{z_*}H₀(0)` together with `H₀ ≤ H₀(0)`), and for a strong `H^m`
solution `sup_{t<T_max}‖u(t)‖_{H¹} < ∞` already forces `T_max = ∞` by the
standard `H¹` local-well-posedness continuation (existence time bounded below by
`c ν³/‖u(t)‖_{H¹}⁴`) plus `H^m`-persistence — no Serrin theorem required. I
recommend keeping F2 but adding a remark that the elementary route suffices; it
reduces the external-input surface by one item.

**(2) Is mean-zero needed, and is it available?** *Needed:* only for the
**homogeneous** inequality `‖u‖_{L⁶} ≤ C_S‖∇u‖_{L²}` used in **Lemma 7** (the
second bound of Main (d)); constants are counterexamples without it. *Not
needed* at Lemma 5, where the inhomogeneous `H¹(𝕋³) ↪ L⁶` applies to the already
uniformly bounded full `H¹` norm — the paper's parenthetical there is
misleading (MINOR-10). *Available:* yes. Zero mean is part of the setting
(`k ≠ 0`), and it is propagated: `d/dt ∫u dx = ν∫Δu − ∫ℙ(u·∇u) = 0` since the
`k = 0` Fourier coefficient of both `Δu` and `ℙ(u·∇u)` vanishes (the latter
because `∫(u·∇)u = ∫∇·(u⊗u) = 0` and `ℙ` preserves that). **This one-line
propagation argument is nowhere in the documents and should be added.**

**(3) Is forward uniqueness the right tool where used?** *No.* It is used at
exactly one place — Lemma 0's positivity of `H₀` — and it is the wrong tool
there: it constrains the solution forward from a hypothetical zero, whereas the
statement to be proved is that no such zero occurs at a *later* time than a
nonzero datum, i.e. a **backward** statement. See CRITICAL-1 for the two
admissible repairs (bootstrap through Lemma 7, which keeps backward uniqueness
out; or cite backward uniqueness and amend the audit). Forward uniqueness *is*
used correctly and harmlessly elsewhere as part of framework input F1 (maximal
solution well-defined).

---

## Summary table

| # | Severity | Location | One-line statement |
|---|---|---|---|
| C1 | CRITICAL | `complete_proof.md` L0, lines 43–49 | positivity of `H₀` is circular; forward uniqueness is the wrong direction |
| C2 | CRITICAL | `theorem_statement.md` (c), lines 73–76 | "(b1) ⟺ finite `Ḣ¹`-bandwidth action" is false; counterexample `e^{−νt}(cos x₂,0,0)` |
| M1 | MAJOR | `theorem_statement.md` (b) | "Exactly one" — the case `T_max = ∞`, `∫KD = ∞` is not excluded |
| M2 | MAJOR | `theorem_statement.md` (e) + Proposition | Proposition refutes only the `R ≡ 0`, all-fields bridge; source caveats dropped |
| M3 | MAJOR | `theorem_statement.md` (a) | "equality cases characterised" is nowhere proved in `complete_proof.md` |
| M4 | MAJOR | `complete_proof.md` Prop., lines 233–237 | covering of `[s_{N₀},∞)` presumes unproved monotonicity of `s_N` |
| M5 | MAJOR | `complete_proof.md` L11 | proof contains literal ellipses; `C₀` never exhibited; `c₁N` fails at small `N` |
| M6 | MAJOR | `dependency_and_gap_audit.md` §1/§2/§4.2/§5 | four audit claims refuted or invalidated by C1/C2 |
| M7 | MAJOR | source notes | two mutually inconsistent factor-2 ledger conventions across the cited sources |
| 1–20 | MINOR | see above | missing hypotheses, misattributed embeddings, numeric rounding, wording |

**Nothing in this report bears on the truth of the underlying mathematics of
(i)–(iv): all identities and every constant reproduced exactly under independent
re-derivation, and every quoted certificate value reproduced from the repository
outputs.** The defects are concentrated in (a) one foundational positivity
argument, (b) one corollary drawn from a correct identity, and (c) the
compression of the source notes' carefully hedged scope statements.

**REFEREE: MAJOR-REVISIONS**
