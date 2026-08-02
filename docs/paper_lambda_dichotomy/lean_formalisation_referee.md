# Adversarial referee report — `SpectralFrontIdentities.lean` and `BandSymmetry.lean`

**Date:** 2026-08-02
**Mode:** read-only. No file was edited; `lake build` was **not** run (another
agent was building). Compilation and axiom-freedom are taken from the author's
report and are *not* independently re-verified here. This report checks
**statements**, not proofs.

**Files under review** (both untracked-new on branch `deps/security-bumps`):

* `formal/NSSingularity/SpectralFrontIdentities.lean` (24.6 kB, 22 declarations)
* `formal/NSSingularity/BandSymmetry.lean` (13.1 kB, 20 declarations)

**Also changed:** `formal/NSSingularity.lean` (+2 imports),
`formal/AxiomAudit.lean` (+2 imports, +31 `#print axioms` lines in the existing
comment style). Both diffs are additive only; nothing pre-existing was touched.
`grep -E "sorry|admit|^axiom|native_decide|implemented_by"` over the two new
files returns only the two docstring sentences asserting their absence.

**Sources of truth:** `docs/paper_lambda_dichotomy/complete_proof.md`
(Lemmas 2, 3, 4, 9) and `docs/paper_lambda_dichotomy/theorem_statement.md`.

---

## Verdict

**No misstatement was found in the mathematical content of either file.**
Every constant, sign, factor and inequality direction I checked matches the
paper. Both of the author's declared deviations are genuine strengthenings and
are correctly documented. No theorem is vacuous.

Six defects were found, **all in the module docstrings** — i.e. in the
scope-honesty layer, which is exactly the layer a reader of a formalisation
trusts most. Two of them (D1, D3) could mislead a reader about what the paper
depends on and about how much of Lemma 4 is formalized. None of them affects a
theorem statement.

---

## 1. Vacuity — every theorem instantiated

### 1.1 A single witness discharges (S1)–(S8)

Take `ι = Fin 2`, `s = univ`, and

```
x     = (1, 4)      eps  = (1, 1)      alpha = (1, -1)     eta = (1, 1)
nu    = 1
```

Then all the defining hypotheses are satisfied with

```
H0 = 2      H1 = 5      H2 = 17     T0 = 0      T1 = -3     G  = 2
mu = 5/2    Vvar = 9/4  Cov = -3/2  Hd0 = -10   Hd1 = -40
```

Checks, all exact:

| fact | check |
|---|---|
| `hH0pos`, `hH1pos`, `hnu` | `2 > 0`, `5 > 0`, `1 > 0` |
| `heps`, `heta` | all entries `= 1 ≥ 0` |
| `hCS` | `alpha i ^2 = 1 ≤ eps i * eta i = 1` (saturated) |
| Gram (S1) | `H0*H2 - H1^2 = 34 - 25 = 9 = 4 * 9/4 = H0^2 * Vvar` |
| centred cov (S2) | `Σ (x-mu)*alpha = -3/2 - 3/2·... = -3 = T1 - mu*T0` |
| **(I.1) / (S3)** | LHS `-40/5 - (-10/2) = -8 + 5 = -3`; RHS `2(-3/2)/(5/2) - 2·1·(9/4)/(5/2) = -6/5 - 9/5 = -3` ✔ |
| **(I.2) / (S5)** | `Cov² = 9/4 ≤ Vvar·G/H0 = (9/4)(2)/2 = 9/4` ✔ (saturated) |
| **(I.3) / (S7)** | `-3 ≤ G/(2·nu·H1) = 2/10 = 1/5` ✔ |
| **(I.4) / (S8)** | gap `= 1/5 - (-3) = 16/5`; `gapCS = (2/(5/2))(√(9/4·1) + 3/2) = (4/5)(3) = 12/5`; `gapSC = (2/(5/2))(3/2 - 1/2)² = 4/5`; sum `= 16/5` ✔ |

The witness has `T0 = 0`, so it *also* discharges
`log_bandwidth_derivative_identity_of_energyNeutral` and
`log_bandwidth_derivative_le_of_energyNeutral` (with `Cov = T1/H0 = -3/2` ✔).
It is non-degenerate in every relevant slot: `Vvar > 0`, `Cov ≠ 0`, `G > 0`,
and **both** defect parts are strictly positive (`12/5` and `4/5`), so (S8) is
not being verified on a collapsed configuration.

`sum_centred_sq_eq`, `centred_variance_eq_moment_ratio`,
`moment_gram_eq_sq_mul_variance`, `variance_nonneg`, `centred_covariance_eq`,
`sq_sum_centred_pairing_le`, `covariance_sq_le_variance_mul_action`,
`abs_covariance_le_sqrt`, `spectral_front_gapCS_nonneg`,
`spectral_front_gapSC_nonneg`, `spectral_front_gap_decomposition`,
`spectral_front_defect_decomposition` are all instantiated by the same witness.

`mul_sqrt_sub_mul_le` (`hW : 0 ≤ W`, `hnu : 0 < nu`) is trivially non-vacuous.

### 1.2 Lemma K

`s = Fin 2`, `x = (1,4)`, `c = (1,1)`, `supNorm = 2`, `gradNorm = √5`,
`nlNorm = 2√5`. Then `hsup : 2 ≤ Σ|c| = 2` ✔; `hgrad : 5 = Σ x·c² = 5` ✔;
`hnl : 2√5 ≤ 2·√5` ✔; positivity ✔. Conclusion `20/25 = 4/5 ≤ Σ x⁻¹ = 5/4` ✔.
Non-vacuous, and non-degenerate (the conclusion is a strict inequality).

### 1.3 Band symmetry

Take the genuine lattice ball `s = {k ∈ ℤ³ : 1 ≤ |k|² ≤ 2}` (18 points: the 6
of norm 1 and the 12 of norm √2) with `w k = (|k|² : ℝ)⁻¹`, i.e.
`f m = ((m : ℤ) : ℝ)⁻¹`. This is *exactly* the intended object (`B_N` with
`N² = 2`), and it satisfies all four closure hypotheses, since a ball is
invariant under the full signed-permutation group.

* `Σ k₁²w = 2·1 + 8·(1/2) = 6`; `Σ|k|²w = 6·1 + 12·1 = 18`; `3·6 = 18` ✔
  (`three_mul_sum_sq_fst_eq_sum_normSq`).
* `Σ k₁k₂w`: the eight `(±1,±1,0)` cancel in sign-pairs; everything else has a
  zero factor. `= 0` ✔.
* `Σ k₁²w = Σ k₂²w = Σ k₃²w = 6` by symmetry ✔.

The author's smaller witness (the six norm-1 points) also works and gives
`3·2 = 6 = Σ|k|²w`.

**No theorem in either file could not be instantiated. Nothing is vacuous.**

---

## 2. Faithfulness — constant-by-constant

### 2.1 (I.1) / Lemma 2 — the factor 2 is right

Paper Lemma 2: `d/dt log N₀² = 2·Cov/μ − 2νV/μ`, with
`Cov = Σ(x_k−μ)α_k/H₀`, `V = Σ p(k)(x_k−μ)²`, `μ = N₀² = H₁/H₀`.

Lean `log_bandwidth_derivative_identity`:
`Hd1/H1 − Hd0/H0 = 2*Cov/mu − 2*nu*Vvar/mu`.

Independently re-derived:

```
Hd1/H1 − Hd0/H0 = 2(T1/H1 − T0/H0) − 2ν(H2/H1 − H1/H0)
T1/H1 − T0/H0   = (T1 − μT0)/H1 = Cov·H0/H1 = Cov/μ
H2/H1 − H1/H0   = (H0H2 − H1²)/(H0H1) = H0²Vvar/(H0H1) = Vvar/μ
```

⇒ `= 2Cov/μ − 2νVvar/μ`. **Exact match, factor 2 confirmed on both terms.**
The left-hand side `Hd1/H1 − Hd0/H0` is the correct algebraic shadow of
`d/dt log(H1/H0)`.

`Cov` and `Vvar` carry **the same `1/H₀` normalisation as the paper**
(`hCov : Cov = (Σ (x−mu)·alpha)/H0`, `hVvar : Vvar = (Σ (x−mu)²·eps)/H0`;
paper: `Cov = Σ(x−μ)α/H₀`, `V = Σ p(x−μ)² = Σ(x−μ)²ε/H₀`). ✔

### 2.2 (I.2) / Lemma 3 — direction and normalisation

Paper: `|Cov| ≤ √(V‖𝒩‖₂²/H₀)`. Lean `abs_covariance_le_sqrt`:
`|Cov| ≤ Real.sqrt (Vvar * G / H0)` with `hG : G = Σ eta = ‖𝒩‖₂²`. Lean parses
`Vvar * G / H0` as `(Vvar*G)/H0` ✔. Squared form
`Cov ^ 2 ≤ Vvar * G / H0` ✔. Direction ✔.

The Cauchy–Schwarz input is the **per-mode** hypothesis
`hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i`, which is the paper's
`|α_k| ≤ |û_k||𝒩̂_k| = √(ε_kη_k)` before summation. It is *not* the aggregated
conclusion. See §3.

### 2.3 (S6) — the square completion

`a·√W − ν·W ≤ a²/(4ν)`. Paper: "maximising the concave function
`√V ↦ α√V − νV`", giving `α²/(4ν)`. ✔ Constant `4ν` correct.

`hW : 0 ≤ W` is **genuinely necessary** and correctly retained: with `a = 0`
and `W < 0`, `Real.sqrt W = 0` and the LHS is `−νW > 0 = ` RHS. The author's
omission of `0 ≤ a` is separately assessed in §3.

### 2.4 (I.3) / Lemma 4 — the constant is `G/(2·nu·H1)`

Paper Lemma 4: `d/dt log N₀² ≤ ‖𝒩‖₂²/(2νH₁) = KD/(2ν)`.
Lean: `Hd1/H1 − Hd0/H0 ≤ G / (2 * nu * H1)`. ✔

Re-derived: `(2/μ)(√(G/H₀)·√Vvar − νVvar) ≤ (2/μ)·(G/H₀)/(4ν) = G/(2νμH₀)`,
and `H₁ = μH₀` ⇒ `= G/(2νH₁)`. **Constant confirmed exactly.** The chain
against the paper's `KD/(2ν)` also closes: `K = ‖𝒩‖₂²/‖∇u‖₂⁴ = G/H₁²`,
`D = H₁`, so `KD/(2ν) = G/(2νH₁)`. ✔

### 2.5 (I.4) / Lemma 4 defect — both parts, exactly

Paper:
`Γ^CS = (2/μ)[√(V‖𝒩‖²/H₀) − Cov]`, `Γ^SC = (2ν/μ)[√V − (1/(2ν))√(‖𝒩‖²/H₀)]²`.

Lean, with `g := G/H0`:
`gapCS mu Cov Vvar g = 2/mu * (√(Vvar*g) − Cov)` ✔,
`gapSC mu nu Vvar g = 2*nu/mu * (√Vvar − √g/(2*nu))^2` ✔.
Both match term for term, including the `1/(2ν)` inside the square.

Re-derived the decomposition independently. With `A = √Vvar`, `B = √(G/H₀)`,
`G = B²H₀`, `H₁ = μH₀`:

```
LHS = B²/(2νμ) − 2Cov/μ + 2νA²/μ
RHS = (2/μ)(AB − Cov) + (2ν/μ)(A² − AB/ν + B²/(4ν²))
    = 2AB/μ − 2Cov/μ + 2νA²/μ − 2AB/μ + B²/(2νμ)
```

The `2AB/μ` terms cancel — the radicals really do drop out — and LHS = RHS
identically. **The paper's claim that the decomposition is an exact rational
identity is faithfully rendered.** ✔

`spectral_front_gapSC_nonneg` needs only `0 < mu`, `0 < nu` (it is a square) —
correct, and correctly not over-hypothesised. `spectral_front_gapCS_nonneg` is
precisely (I.2) and correctly consumes `hCS`. ✔

### 2.6 Lemma K

`lemmaK_bound : nlNorm ^ 2 / gradNorm ^ 4 ≤ ∑ i ∈ s, (x i)⁻¹`, i.e.
`K ≤ S_N` in the paper's unpaired normalisation. Re-derived:
`‖𝒩‖ ≤ ‖u‖_∞‖∇u‖₂` and `‖u‖_∞ ≤ Σ|û_k| = Σ|k|⁻¹·|k||û_k| ≤ √(S_N)·‖∇u‖₂`
give `K = ‖𝒩‖²/‖∇u‖₂⁴ ≤ S_N`. ✔ The docstring's claim that "for `x i = |k|²`
over a band the right-hand side is the lattice sum `S_N` of the paper" is
correct (`Σ(|k|²)⁻¹ = Σ|k|⁻²`). But see **D1** — Lemma K is not in the paper.

### 2.7 Lemma 9 / band symmetry

Paper: `Σ_{k≠0} k_i k_j f(k) = δ_ij · (1/3) · Σ_{k≠0}|k|² f(k)` for radial `f`.

Lean `band_symmetry` states the three off-diagonals `= 0` (the `δ_ij = 0`
case), the three diagonals equal, and `3·Σk₁²w = Σ|k|²w` (the `δ_ij = 1` case,
rearranged). Modulo `i ↔ j` symmetry of `k_ik_j`, this **is** Lemma 9. ✔

Radiality is faithful: `hradial : ∀ k ∈ s, w k = f (k.1^2 + k.2.1^2 + k.2.2^2)`
is `f(k) = φ(|k|²)`, imposed only on `s`, and the four invariances are
*derived* from it (`radial_flipFst` etc.) rather than assumed — as the
docstring claims. ✔

Sums over `s` may include `k = 0`; harmless, since `0` contributes `0` to every
sum here, and the intended `B_N` excludes it anyway.

### 2.8 Nothing else silently changed

`d/dt log N₀² = Ḣ₁/H₁ − Ḣ₀/H₀`, the ledger `Ḣ₀ = 2T₀ − 2νH₁`,
`Ḣ₁ = 2T₁ − 2νH₂` (paper Lemma 1), `Σ η_k = ‖𝒩‖₂²`, `D = H₁`,
`K = ‖𝒩‖₂²/‖∇u‖₂⁴` — every constant in these bridging relations is reproduced
correctly wherever it appears as a hypothesis.

---

## 3. Hidden strength — none found

**Checked and clean:**

* **Energy neutrality is not smuggled in.** `log_bandwidth_derivative_identity`
  and `log_bandwidth_derivative_le` leave `T0` free (`hT0 : T0 = Σ alpha`, no
  constraint). `T0 = 0` appears only as the explicit named hypothesis
  `hneutral` in the two `_of_energyNeutral` corollaries, and those are proved
  *from* the general forms, not the other way round. The paper *derives*
  `T₀ = 0` (Lemma 1); Lean *assumes* it. That is a weakening, correctly
  disclosed, not hidden strength.
* **The Cauchy–Schwarz input is per-mode, not the conclusion.** `hCS` is
  `∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i`. The aggregated bound (I.2) is
  *derived* (`sq_sum_centred_pairing_le` →
  `covariance_sq_le_variance_mul_action`). No circularity.
* **The Gram identity is not double-assumed.** `hvar : H0*H2 − H1^2 = H0^2*Vvar`
  is a hypothesis in the two "bare algebra" lemmas
  (`log_bandwidth_derivative_identity`, `spectral_front_gap_decomposition`),
  but in both packaged theorems (`log_bandwidth_derivative_le`,
  `spectral_front_defect_decomposition`) it is discharged from the sum
  definitions via `moment_gram_eq_sq_mul_variance`. No hypothesis is assumed
  that the file elsewhere proves.
* **`Hd0`/`Hd1` are unconstrained reals** pinned only by the two ledger
  equations. They cannot secretly encode a derivative.
* **`hmupos` in the two `gap*_nonneg` lemmas** is a hypothesis, not a derived
  fact — a weakening, harmless, and discharged from `hH0pos`/`hH1pos` in the
  packaged version.
* **`hflip2` in `band_symmetry` is redundant** (`flipSnd = swapFstSnd ∘ flipFst
  ∘ swapFstSnd`), i.e. the hypothesis set is *stronger than necessary* — but
  strictly weakens the theorem, and the intended object satisfies it. Cosmetic
  only (D6).

**Author's declared deviations — both confirmed as genuine strengthenings:**

1. `mul_sqrt_sub_mul_le` without `0 ≤ a`. **Correct.** The completion
   `ν(√W − a/(2ν))² ≥ 0` expands to `νW − a√W + a²/(4ν) ≥ 0` for *every* real
   `a`; the sign of `a` never enters. Dropping the hypothesis strictly
   strengthens the lemma. `hW : 0 ≤ W` is retained and is necessary (§2.3).
2. `log_bandwidth_derivative_le` without `T0 = 0`. **Correct.** Because `Cov`
   is defined as the *centred* covariance `(Σ(x−μ)α)/H₀ = (T₁ − μT₀)/H₀`, both
   (I.1) and (I.3) hold for arbitrary `T₀`; neutrality is needed only to
   rewrite `Cov = T₁/H₀`. The paper's literal form is preserved as the
   corollary, which genuinely consumes `hneutral`. Same for (S8). This is the
   right way to have done it.

---

## 4. Scope honesty

The docstrings are, on the whole, unusually candid — better than the repository
norm. `SpectralFrontIdentities.lean` states plainly that there is "no PDE, no
Fourier transform, no Leray projection, no notion of derivative, and no
Navier–Stokes solution anywhere in this file", that `Hd0`/`Hd1` are plain
reals, that `hCS` is an input, that `T0 = 0` is never assumed silently, and
that Lemma K's `hsup`/`hnl` are unformalised analytic inputs.
`BandSymmetry.lean` states that all sums are finite, that no absolute-summability
or rearrangement argument for the infinite lattice is made, and that the
connection to `H₀`/`H₁` "is made in the paper, not in Lean".

**A reader would not mistake these for PDE theorems** — with the exceptions D1
and D3 below.

### Defects

**D1 (moderate) — the file claims a paper dependency that does not exist.**
`SpectralFrontIdentities.lean` line 11–14: "Its Lemmas 2, 3 and 4 consume four
algebraic facts about the modal moments … **plus the 'Lemma K' bandwidth bound
on the front functional `K`**."

`complete_proof.md` contains **no Lemma K**. `grep -n "Lemma K"` over
`docs/paper_lambda_dichotomy/` returns nothing; the bound `K ≤ S_N` appears
only in `docs/candidates/CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md:101`
and in `STATUS.md`. Lemmas 2–4 of the paper do not consume it, and no lemma in
`complete_proof.md` uses it: `S_N` occurs there only inside Lemma 11's lattice
bounds on `S_N^χ`, an unrelated role. Lemma K is therefore a *fifth* result
imported from a different document, presented as if the paper depended on it.
The mathematics of `lemmaK_bound` is fine (§2.6); the provenance sentence is
not. **Fix:** cite the candidate document for Lemma K, and drop it from the
list of things the paper's Lemmas 2–4 consume.

**D2 (moderate) — normalisation label collision.** The docstring attributes
"(I.1)–(I.4)" to "the repository". The repository's own (I.1)
(`CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md`) is
`½ d/dt log N_r² = Cov/μ − νV_r/μ` in a **paired** modal convention with
`K = G₀/(2H₁²)`, `D = 2H₁`, `G₀ = ½‖ℙ(u·∇u)‖₂²`. The Lean statement is the
**unpaired** convention of `complete_proof.md` Lemma 2 (factor 2, no ½,
`K = ‖𝒩‖₂²/‖∇u‖₂⁴`, `D = H₁`). The two do reconcile — paired `H₁` is half the
unpaired one, and `G₀` is half `‖𝒩‖²` — and the Lean form is the correct one
for the designated source. But a reader who follows the docstring's citation
to the repository's "(I.1)" will find a spurious factor of 2 and have to
re-derive the reconciliation. **Fix:** cite `complete_proof.md` Lemmas 2–4 by
number and state that the unpaired full-lattice convention is used.

**D3 (minor) — (S8) is billed as more of Lemma 4 than it is.** "Together with
(S7) this is Lemma 4 of the paper, minus every analytic ingredient." Lemma 4
also concludes "**hence `Λ` is non-increasing**", obtained by integrating the
pointwise inequality in `t`. Nothing in the file integrates anything; the
monotone `Λ` — the headline object of the whole paper — is not formalized in
any form. "Minus every analytic ingredient" is doing a lot of silent work
here. **Fix:** say explicitly that the time-integration step and the monotone
`Λ` are not formalized.

**D4 (minor) — one unformalised analytic input is unlisted.** The Lemma K
docstring names `hsup` (Fourier inversion) and `hnl` (bilinear + Leray) as the
analytic inputs. `hgrad : gradNorm ^ 2 = ∑ x·c²` is a third one — it is
Parseval for `‖∇u‖₂²`. It is an equality hypothesis, so it cannot make the
theorem wrong, but the list reads as exhaustive and is not.

**D5 (minor) — `BandSymmetry` overstates what (c) drives.**
`three_mul_sum_sq_fst_eq_sum_normSq` is called "the identity behind the paper's
exact family laws `H₀ = (2/3)‖v₀‖²T_N` and `H₁ = (2/3)‖v₀‖²S_N` (Lemma 10)".
Lemma 10 needs `|P_kv₀|² = ‖v₀‖² − (k·v₀)²/|k|²` with
`(k·v₀)² = Σ_{i,j} k_ik_j v_iv_j`, so it consumes the **off-diagonal**
vanishing (a) as well as (b)/(c); and the contraction against `v₀` — the step
that actually produces the `2/3` — is not formalized at all. (c) alone is not
"the engine". **Fix:** say (a)+(b)+(c) are the ingredients and that the
contraction with `v₀` is done in the paper, not in Lean.

**D6 (cosmetic).** `hflip2` is redundant (see §3). The "What is formalized"
list omits a bullet for (S4) although the section header `(S4)–(S5)` exists,
and omits `sum_centred_sq_eq` and `variance_nonneg` from the (S1) bullet.

### The analytic bridge, stated in full

Neither file proves, and neither claims to prove: the sup-norm/Fourier-inversion
bound; `L²`-boundedness of the Leray projection; existence or continuity of
`Ḣ₀`, `Ḣ₁`, or `d/dt log N₀²`; Parseval identifying `Σ x^r ε` with
`‖A^{r/2}u‖₂²`; the ledger `Ḣ_r = 2T_r − 2νH_{r+1}` (Lemma 1); energy
neutrality `T₀ = 0` (Lemma 1); positivity `H₀ > 0` (Lemma 0b); the passage from
finite Finset sums to infinite-lattice/PDE quantities; the time integration
yielding the monotone `Λ`; the dichotomy (Lemma 5); the action representation
(Lemma 6); the dominations (Lemma 7); Osgood closure (Lemma 8); and everything
in Theorem O beyond Lemma 9's finite algebra (Lemmas 10, 11, the Capacity
Theorem). With D1/D3 corrected, the docstrings say so adequately.

---

## 5. Band symmetry — hypotheses, intended object, and (c)

**Is the symmetry hypothesis strong enough?** Yes, and it is used correctly:

| conclusion | hypothesis used | mechanism | correct? |
|---|---|---|---|
| `Σ k₁k₂w = 0` | `flipFst` closure | `(−k₁)k₂ = −k₁k₂`, `w` invariant | ✔ |
| `Σ k₁k₃w = 0` | `flipFst` closure | `(−k₁)k₃ = −k₁k₃` | ✔ |
| `Σ k₂k₃w = 0` | `flipSnd` closure | `k₁(−k₂)k₃ ↦ −k₂k₃` | ✔ |
| `Σ k₁²w = Σ k₂²w` | `swapFstSnd` closure | `Finset.sum_nbij'` along an involution | ✔ |
| `Σ k₁²w = Σ k₃²w` | `swapFstThd` closure | same | ✔ |

Note the hypotheses are *closure* conditions `∀ k ∈ s, σ k ∈ s`, not
`σ '' s = s`; combined with involutivity this is equivalent, and
`sum_comp_involution` uses exactly `hmem`/`hinv` on both sides. Correct.

**Radial-weight hypothesis.** Used correctly and minimally: `hradial` is
required only on `s`, and each `radial_*` lemma discharges
`w (σ k) = w k` by rewriting through `f` and closing the integer identity with
`ring` (e.g. `(−k₁)² + k₂² + k₃² = k₁² + k₂² + k₃²`). Radiality is genuinely
*needed* — without it the off-diagonal sums do not vanish — and it is not
over-used.

**Does (c) follow from (a)+(b)?** (c) follows from **(b) alone**:
`3Σk₁²w = Σk₁²w + Σk₂²w + Σk₃²w = Σ(k₁²+k₂²+k₃²)w`. The Lean proof indeed uses
only `h12` and `h13` plus `Finset.sum_add_distrib`. Correct; (a) is not needed
for (c). (But (a) *is* needed for Lemma 10 — see D5.)

**Is it satisfied by the intended object?** Yes. `B_N = {k : 1 ≤ |k|² ≤ N²}`
and the smooth-truncation support `{k ≠ 0 : χ(|k|/N) ≠ 0}` are both invariant
under the full hyperoctahedral group, hence under all four maps; and the
weights the paper feeds in — `χ(|k|/N)²|k|^{-6}`, `χ(|k|/N)²|k|^{-4}`,
`χ(|k|/N)|k|^{-4}` — are radial in exactly the `f(|k|²)` sense, with
`f : ℤ → ℝ` well-typed since `|k|²` is an integer. Concretely instantiated in
§1.3.

**One genuine limitation, correctly disclosed.** The paper's Lemma 9 is stated
on the infinite lattice under absolute summability; Lean proves the finite
case. The docstring says so. Worth recording in the paper's favour: for the
*smooth truncated* family actually used by Theorem O, `χ` is supported in
`|k| ≤ N`, so all three weights are **finitely supported** (`complete_proof.md`
Lemma 10: "both finitely supported, hence absolutely summable"). The finite
version is therefore fully adequate for the only application in the paper. The
generality gap is real but not load-bearing.

---

## 6. Summary of required changes

None to any theorem statement. Docstring edits only:

1. **D1** — remove the claim that the paper's Lemmas 2–4 consume Lemma K;
   attribute Lemma K to `docs/candidates/CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md`.
2. **D2** — cite `complete_proof.md` Lemmas 2–4 rather than the repository's
   "(I.1)–(I.4)", and state that the unpaired full-lattice convention is used.
3. **D3** — state that the time integration and the monotone `Λ` are not
   formalized.
4. **D4** — add `hgrad` to Lemma K's list of unformalised analytic inputs.
5. **D5** — correct the attribution of Lemma 10's `2/3` to (a)+(b)+(c) plus an
   unformalised contraction with `v₀`.
6. **D6** — optional: note `hflip2`'s redundancy; complete the (S1)/(S4)
   bullet lists.

---

## 7. Caveats on this report

* `lake build` was not run and the axiom audit was not re-executed. If the
  author's build report is wrong, this report says nothing about it. Everything
  above concerns **statements**, which I read directly.
* I did not check proof *tactics* for admissibility beyond confirming the
  absence of `sorry`, `admit`, `axiom`, `native_decide` and `implemented_by`.
* The author reports checking non-vacuity in a scratch file that was then
  deleted; I could not inspect it, so §1 is an independent pen-and-paper
  reconstruction, not a confirmation of theirs.
