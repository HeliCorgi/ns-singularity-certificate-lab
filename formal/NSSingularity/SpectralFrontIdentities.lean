/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import Mathlib

/-!
# Spectral-front identities (I.1)–(I.4) and Lemma K: the finite-sum algebra

`docs/paper_lambda_dichotomy/complete_proof.md` proves a bandwidth–dissipation
dichotomy for the periodic Navier–Stokes equations.  Its Lemmas 2, 3 and 4
consume four algebraic facts about the modal moments — the identities the
repository calls (I.1)–(I.4).  **This file is exactly those facts, and nothing
else.**

The file additionally formalizes the finite core of the bandwidth bound that
`docs/candidates/CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md` calls
"Lemma K".  Provenance, stated precisely because an earlier draft of this
docstring got it wrong: `complete_proof.md` states no Lemma K and no lemma of
that paper consumes one; the bound belongs to the candidate document, and it is
formalized here only because it is the same finite Cauchy–Schwarz inequality
with different weights.

Everything here is finite-sum real algebra.  There is no PDE, no Fourier
transform, no Leray projection, no notion of derivative, and no Navier–Stokes
solution anywhere in this file.

## What is formalized

Let `s : Finset ι` carry per-mode reals `eps i ≥ 0` (modal energy `|û_k|²`),
`alpha i` (modal growth pairing `Re⟨û_k, N̂_k⟩`), `eta i ≥ 0` (modal nonlinear
power `|N̂_k|²`) and `x i` (the weight `|k|²`), with moments
`H0 = Σ eps`, `H1 = Σ x·eps`, `H2 = Σ x²·eps`, `T0 = Σ alpha`, `T1 = Σ x·alpha`,
`G = Σ eta`, and `mu = H1/H0`, `Cov = (Σ (x-mu)·alpha)/H0`,
`Vvar = (Σ (x-mu)²·eps)/H0`.

* **(S1)** `centred_variance_eq_moment_ratio` and `moment_gram_eq_sq_mul_variance`
  — the centred variance identity in both forms.
* **(S2)** `centred_covariance_eq` — `Σ (x-mu)·alpha = T1 - mu·T0`.
* **(S3)** `log_bandwidth_derivative_identity` — identity (I.1), the logarithmic
  bandwidth derivative identity, in its general form, and
  `log_bandwidth_derivative_identity_of_energyNeutral`, its `T0 = 0` corollary.
* **(S5)** `covariance_sq_le_variance_mul_action` and `abs_covariance_le_sqrt`
  — identity (I.2) in squared and un-squared form.
* **(S6)** `mul_sqrt_sub_mul_le` — the square-completion inequality.
* **(S7)** `log_bandwidth_derivative_le` — identity (I.3), and
  `log_bandwidth_derivative_le_of_energyNeutral`, the paper's `T0 = 0` form.
* **(S8)** `spectral_front_gap_decomposition`, `spectral_front_gapCS_nonneg`,
  `spectral_front_gapSC_nonneg`, `spectral_front_defect_decomposition` —
  identity (I.4), the exact gap decomposition with both parts nonnegative.
* **(S9)** `sq_sum_abs_le_sum_inv_mul_sum_mul_sq` and `lemmaK_bound` — Lemma K.
  (Lemma K lives here rather than in its own file because it is the same finite
  Cauchy–Schwarz inequality as (I.2), applied to different weights.)

## What is deliberately **not** proved

* `Hd0` and `Hd1` are **plain real numbers** constrained by the two named
  ledger hypotheses `Hd0 = 2·T0 - 2·nu·H1` and `Hd1 = 2·T1 - 2·nu·H2`.  No
  derivative of any moment is defined, and nothing here asserts that a
  Navier–Stokes solution satisfies those relations.  That is Lemma 1 of the
  paper, an analytic statement, and it is not formalized.
* The modal Cauchy–Schwarz input `alpha i ^ 2 ≤ eps i * eta i` is a **named
  hypothesis** (`hCS`).  In the paper it comes from
  `|Re⟨û_k, N̂_k⟩| ≤ |û_k| |N̂_k|` on `ℂ³`; here it is an input.
* Energy neutrality `T0 = 0` is a **named hypothesis**, never assumed silently.
  The general (`T0` arbitrary) form of (I.1) is proved first, and `T0 = 0` is a
  corollary of it.
* In Lemma K, `hsup : supNorm ≤ Σ |c i|` (Fourier inversion) and
  `hnl : nlNorm ≤ supNorm * gradNorm` (the bilinear bound together with the
  boundedness of the Leray projection on `L²`) are the **unformalised analytic
  inputs**.  They appear as hypotheses; no Fourier inversion theorem and no
  projection bound is proved in this file.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/

namespace NSSingularity

noncomputable section

variable {ι : Type*}

/-! ## (S1) The centred variance identity -/

/-- Expanding the centred second moment into the three raw moments.  This is
the step the repository computes as a direct centred sum rather than as a
moment difference. -/
theorem sum_centred_sq_eq (s : Finset ι) (x eps : ι → ℝ) {H0 H1 H2 mu : ℝ}
    (hH0 : H0 = ∑ i ∈ s, eps i)
    (hH1 : H1 = ∑ i ∈ s, x i * eps i)
    (hH2 : H2 = ∑ i ∈ s, x i ^ 2 * eps i) :
    ∑ i ∈ s, (x i - mu) ^ 2 * eps i = H2 - 2 * mu * H1 + mu ^ 2 * H0 := by
  have hpt : ∀ i ∈ s, (x i - mu) ^ 2 * eps i
      = x i ^ 2 * eps i - 2 * mu * (x i * eps i) + mu ^ 2 * eps i := by
    intro i _
    ring
  rw [Finset.sum_congr rfl hpt, Finset.sum_add_distrib, Finset.sum_sub_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, hH0, hH1, hH2]

/-- **(S1), ratio form.**  The centred variance equals the moment difference:
`(Σ (x - mu)² eps)/H0 = H2/H0 - mu²` when `mu = H1/H0`. -/
theorem centred_variance_eq_moment_ratio (s : Finset ι) (x eps : ι → ℝ)
    {H0 H1 H2 mu : ℝ}
    (hH0 : H0 = ∑ i ∈ s, eps i)
    (hH1 : H1 = ∑ i ∈ s, x i * eps i)
    (hH2 : H2 = ∑ i ∈ s, x i ^ 2 * eps i)
    (hH0pos : 0 < H0) (hmu : mu = H1 / H0) :
    (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0 = H2 / H0 - mu ^ 2 := by
  rw [sum_centred_sq_eq s x eps hH0 hH1 hH2, hmu]
  field_simp
  ring

/-- **(S1), Gram form.**  `H0·H2 - H1² = H0²·Vvar`.  This is the shape used in
the second bracket of (I.1). -/
theorem moment_gram_eq_sq_mul_variance (s : Finset ι) (x eps : ι → ℝ)
    {H0 H1 H2 mu Vvar : ℝ}
    (hH0 : H0 = ∑ i ∈ s, eps i)
    (hH1 : H1 = ∑ i ∈ s, x i * eps i)
    (hH2 : H2 = ∑ i ∈ s, x i ^ 2 * eps i)
    (hH0pos : 0 < H0) (hmu : mu = H1 / H0)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0) :
    H0 * H2 - H1 ^ 2 = H0 ^ 2 * Vvar := by
  rw [hVvar, sum_centred_sq_eq s x eps hH0 hH1 hH2, hmu]
  field_simp
  ring

/-- The centred variance is nonnegative when every modal energy is. -/
theorem variance_nonneg (s : Finset ι) (x eps : ι → ℝ) {H0 mu Vvar : ℝ}
    (hH0pos : 0 < H0)
    (heps : ∀ i ∈ s, 0 ≤ eps i)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0) :
    0 ≤ Vvar := by
  rw [hVvar]
  refine div_nonneg (Finset.sum_nonneg fun i hi => ?_) hH0pos.le
  exact mul_nonneg (sq_nonneg _) (heps i hi)

/-! ## (S2) The covariance form -/

/-- **(S2).**  The centred pairing is the raw pairing minus the mean weight
times the total pairing: `Σ (x - mu)·alpha = T1 - mu·T0`. -/
theorem centred_covariance_eq (s : Finset ι) (x alpha : ι → ℝ) {T0 T1 mu : ℝ}
    (hT0 : T0 = ∑ i ∈ s, alpha i)
    (hT1 : T1 = ∑ i ∈ s, x i * alpha i) :
    ∑ i ∈ s, (x i - mu) * alpha i = T1 - mu * T0 := by
  have hpt : ∀ i ∈ s, (x i - mu) * alpha i = x i * alpha i - mu * alpha i := by
    intro i _
    ring
  rw [Finset.sum_congr rfl hpt, Finset.sum_sub_distrib, ← Finset.mul_sum, hT0, hT1]

/-! ## (S3) Identity (I.1): the logarithmic bandwidth derivative identity -/

/-- **(S3), general form — identity (I.1).**  With `Hd0`, `Hd1` plain reals
constrained by the two ledger relations, `mu = H1/H0`,
`Cov = (T1 - mu·T0)/H0` and `H0·H2 - H1² = H0²·Vvar`:

`Hd1/H1 - Hd0/H0 = 2·Cov/mu - 2·nu·Vvar/mu`.

No energy-neutrality assumption is made here; `T0` is arbitrary.  Nothing in
this statement is a derivative: `Hd0` and `Hd1` are real numbers. -/
theorem log_bandwidth_derivative_identity
    {H0 H1 H2 T0 T1 Hd0 Hd1 nu mu Cov Vvar : ℝ}
    (hH0pos : 0 < H0) (hH1pos : 0 < H1)
    (hmu : mu = H1 / H0)
    (hCov : Cov = (T1 - mu * T0) / H0)
    (hvar : H0 * H2 - H1 ^ 2 = H0 ^ 2 * Vvar)
    (hHd0 : Hd0 = 2 * T0 - 2 * nu * H1)
    (hHd1 : Hd1 = 2 * T1 - 2 * nu * H2) :
    Hd1 / H1 - Hd0 / H0 = 2 * Cov / mu - 2 * nu * Vvar / mu := by
  have hH0ne : H0 ≠ 0 := ne_of_gt hH0pos
  have hH1ne : H1 ≠ 0 := ne_of_gt hH1pos
  have hH2 : H2 = (H0 ^ 2 * Vvar + H1 ^ 2) / H0 := by
    field_simp
    linarith [hvar]
  subst hHd0 hHd1 hCov hmu
  rw [hH2]
  field_simp
  ring

/-- **(S3), energy-neutral corollary.**  With the additional named hypothesis
`T0 = 0` — energy neutrality, Lemma 1 of the paper, an analytic input that is
*not* proved here — the covariance is simply `T1/H0`. -/
theorem log_bandwidth_derivative_identity_of_energyNeutral
    {H0 H1 H2 T0 T1 Hd0 Hd1 nu mu Cov Vvar : ℝ}
    (hH0pos : 0 < H0) (hH1pos : 0 < H1)
    (hmu : mu = H1 / H0)
    (hneutral : T0 = 0)
    (hCov : Cov = T1 / H0)
    (hvar : H0 * H2 - H1 ^ 2 = H0 ^ 2 * Vvar)
    (hHd0 : Hd0 = 2 * T0 - 2 * nu * H1)
    (hHd1 : Hd1 = 2 * T1 - 2 * nu * H2) :
    Hd1 / H1 - Hd0 / H0 = 2 * Cov / mu - 2 * nu * Vvar / mu := by
  refine log_bandwidth_derivative_identity hH0pos hH1pos hmu ?_ hvar hHd0 hHd1
  rw [hCov, hneutral]
  ring

/-! ## (S4)–(S5) Identity (I.2): the modal Cauchy–Schwarz bound -/

/-- The finite Cauchy–Schwarz step behind (I.2), before dividing by `H0`.

`hCS` is the modal Cauchy–Schwarz **input**: in the paper it is
`|Re⟨û_k, N̂_k⟩| ≤ |û_k| |N̂_k|`, which is not proved here. -/
theorem sq_sum_centred_pairing_le (s : Finset ι) (x eps alpha eta : ι → ℝ)
    {mu : ℝ}
    (heps : ∀ i ∈ s, 0 ≤ eps i) (heta : ∀ i ∈ s, 0 ≤ eta i)
    (hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i) :
    (∑ i ∈ s, (x i - mu) * alpha i) ^ 2
      ≤ (∑ i ∈ s, (x i - mu) ^ 2 * eps i) * ∑ i ∈ s, eta i := by
  refine Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul s
    (f := fun i => (x i - mu) ^ 2 * eps i) (g := eta)
    (fun i hi => mul_nonneg (sq_nonneg _) (heps i hi)) heta ?_
  intro i hi
  have h : (x i - mu) ^ 2 * alpha i ^ 2 ≤ (x i - mu) ^ 2 * (eps i * eta i) :=
    mul_le_mul_of_nonneg_left (hCS i hi) (sq_nonneg _)
  calc ((x i - mu) * alpha i) ^ 2 = (x i - mu) ^ 2 * alpha i ^ 2 := by ring
    _ ≤ (x i - mu) ^ 2 * (eps i * eta i) := h
    _ = (x i - mu) ^ 2 * eps i * eta i := by ring

/-- **(S5), squared form — identity (I.2).**  `Cov² ≤ Vvar·G/H0`. -/
theorem covariance_sq_le_variance_mul_action (s : Finset ι)
    (x eps alpha eta : ι → ℝ) {H0 G mu Cov Vvar : ℝ}
    (hH0pos : 0 < H0)
    (hG : G = ∑ i ∈ s, eta i)
    (hCov : Cov = (∑ i ∈ s, (x i - mu) * alpha i) / H0)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0)
    (heps : ∀ i ∈ s, 0 ≤ eps i) (heta : ∀ i ∈ s, 0 ≤ eta i)
    (hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i) :
    Cov ^ 2 ≤ Vvar * G / H0 := by
  have hkey := sq_sum_centred_pairing_le s x eps alpha eta (mu := mu) heps heta hCS
  have hH0ne : H0 ≠ 0 := ne_of_gt hH0pos
  have hrhs : Vvar * G / H0
      = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) * (∑ i ∈ s, eta i) / H0 ^ 2 := by
    rw [hVvar, hG]
    field_simp
  rw [hCov, div_pow, hrhs]
  exact div_le_div_of_nonneg_right hkey (sq_nonneg H0)

/-- **(S5), un-squared form.**  `|Cov| ≤ √(Vvar·G/H0)`. -/
theorem abs_covariance_le_sqrt (s : Finset ι) (x eps alpha eta : ι → ℝ)
    {H0 G mu Cov Vvar : ℝ}
    (hH0pos : 0 < H0)
    (hG : G = ∑ i ∈ s, eta i)
    (hCov : Cov = (∑ i ∈ s, (x i - mu) * alpha i) / H0)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0)
    (heps : ∀ i ∈ s, 0 ≤ eps i) (heta : ∀ i ∈ s, 0 ≤ eta i)
    (hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i) :
    |Cov| ≤ Real.sqrt (Vvar * G / H0) := by
  have hsq := covariance_sq_le_variance_mul_action s x eps alpha eta hH0pos hG
    hCov hVvar heps heta hCS
  have h : |Cov| ^ 2 ≤ Vvar * G / H0 := by rwa [sq_abs]
  exact Real.le_sqrt_of_sq_le h

/-! ## (S6) The square-completion inequality -/

/-- **(S6).**  `a·√W - nu·W ≤ a²/(4·nu)` for `0 ≤ W` and `0 < nu`.

The paper states this for `a ≥ 0`; nonnegativity of `a` is in fact not needed,
so it is omitted rather than carried as an unused hypothesis. -/
theorem mul_sqrt_sub_mul_le {a W nu : ℝ} (hW : 0 ≤ W) (hnu : 0 < nu) :
    a * Real.sqrt W - nu * W ≤ a ^ 2 / (4 * nu) := by
  have hnune : nu ≠ 0 := ne_of_gt hnu
  have hsq : Real.sqrt W ^ 2 = W := Real.sq_sqrt hW
  have hkey : nu * (Real.sqrt W - a / (2 * nu)) ^ 2
      = nu * Real.sqrt W ^ 2 - a * Real.sqrt W + a ^ 2 / (4 * nu) := by
    field_simp
    ring
  rw [hsq] at hkey
  have hnn : 0 ≤ nu * (Real.sqrt W - a / (2 * nu)) ^ 2 :=
    mul_nonneg hnu.le (sq_nonneg _)
  linarith [hkey ▸ hnn]

/-! ## (S7) Identity (I.3) -/

/-- **(S7) — identity (I.3).**  Combining (I.1), (I.2) and the square
completion:

`Hd1/H1 - Hd0/H0 ≤ G/(2·nu·H1)`.

The constant is exactly `G/(2·nu·H1)`: the square completion produces
`G/(2·nu·mu·H0)`, and `H1 = mu·H0` converts it.

Energy neutrality `T0 = 0` is *not* needed in this form, because `Cov` is the
**centred** covariance `(Σ (x - mu)·alpha)/H0`; it is needed only to rewrite
`Cov` as `T1/H0`, which is
`log_bandwidth_derivative_le_of_energyNeutral` below. -/
theorem log_bandwidth_derivative_le (s : Finset ι) (x eps alpha eta : ι → ℝ)
    {H0 H1 H2 T0 T1 G Hd0 Hd1 nu mu Cov Vvar : ℝ}
    (hH0 : H0 = ∑ i ∈ s, eps i)
    (hH1 : H1 = ∑ i ∈ s, x i * eps i)
    (hH2 : H2 = ∑ i ∈ s, x i ^ 2 * eps i)
    (hT0 : T0 = ∑ i ∈ s, alpha i)
    (hT1 : T1 = ∑ i ∈ s, x i * alpha i)
    (hG : G = ∑ i ∈ s, eta i)
    (hH0pos : 0 < H0) (hH1pos : 0 < H1) (hnu : 0 < nu)
    (heps : ∀ i ∈ s, 0 ≤ eps i) (heta : ∀ i ∈ s, 0 ≤ eta i)
    (hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i)
    (hmu : mu = H1 / H0)
    (hCov : Cov = (∑ i ∈ s, (x i - mu) * alpha i) / H0)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0)
    (hHd0 : Hd0 = 2 * T0 - 2 * nu * H1)
    (hHd1 : Hd1 = 2 * T1 - 2 * nu * H2) :
    Hd1 / H1 - Hd0 / H0 ≤ G / (2 * nu * H1) := by
  have hH0ne : H0 ≠ 0 := ne_of_gt hH0pos
  have hmupos : 0 < mu := by rw [hmu]; exact div_pos hH1pos hH0pos
  have hvar : H0 * H2 - H1 ^ 2 = H0 ^ 2 * Vvar :=
    moment_gram_eq_sq_mul_variance s x eps hH0 hH1 hH2 hH0pos hmu hVvar
  have hCov' : Cov = (T1 - mu * T0) / H0 := by
    rw [hCov, centred_covariance_eq s x alpha hT0 hT1]
  have hident : Hd1 / H1 - Hd0 / H0 = 2 * Cov / mu - 2 * nu * Vvar / mu :=
    log_bandwidth_derivative_identity hH0pos hH1pos hmu hCov' hvar hHd0 hHd1
  -- the two analytic ingredients, as finite-sum facts
  have hVnn : 0 ≤ Vvar := variance_nonneg s x eps hH0pos heps hVvar
  have hGnn : 0 ≤ G := by
    rw [hG]; exact Finset.sum_nonneg fun i hi => heta i hi
  have hgnn : (0 : ℝ) ≤ G / H0 := div_nonneg hGnn hH0pos.le
  have hCovle : Cov ≤ Real.sqrt (G / H0) * Real.sqrt Vvar := by
    have habs := abs_covariance_le_sqrt s x eps alpha eta hH0pos hG hCov hVvar
      heps heta hCS
    have hrw : Vvar * G / H0 = Vvar * (G / H0) := by ring
    calc Cov ≤ |Cov| := le_abs_self Cov
      _ ≤ Real.sqrt (Vvar * G / H0) := habs
      _ = Real.sqrt (G / H0) * Real.sqrt Vvar := by
          rw [hrw, Real.sqrt_mul hVnn, mul_comm]
  -- the square completion
  have hcompl : Real.sqrt (G / H0) * Real.sqrt Vvar - nu * Vvar
      ≤ Real.sqrt (G / H0) ^ 2 / (4 * nu) := mul_sqrt_sub_mul_le hVnn hnu
  rw [Real.sq_sqrt hgnn] at hcompl
  have hstep : Cov - nu * Vvar ≤ G / H0 / (4 * nu) := by linarith
  have hfinal : 2 * Cov / mu - 2 * nu * Vvar / mu ≤ G / (2 * nu * H1) := by
    have hH1eq : H1 = mu * H0 := by
      rw [hmu]; field_simp
    have h2 : 2 / mu * (Cov - nu * Vvar) ≤ 2 / mu * (G / H0 / (4 * nu)) :=
      mul_le_mul_of_nonneg_left hstep (by positivity)
    have hleft : 2 / mu * (Cov - nu * Vvar) = 2 * Cov / mu - 2 * nu * Vvar / mu := by
      field_simp
    have hright : 2 / mu * (G / H0 / (4 * nu)) = G / (2 * nu * H1) := by
      rw [hH1eq]
      field_simp
      ring
    linarith [hleft ▸ hright ▸ h2]
  linarith [hident ▸ hfinal]

/-- **(S7), the paper's form.**  Under the named energy-neutrality hypothesis
`T0 = 0` — Lemma 1 of the paper, an analytic input that is *not* proved here —
the covariance is `T1/H0` and (I.3) reads

`Hd1/H1 - Hd0/H0 ≤ G/(2·nu·H1)`. -/
theorem log_bandwidth_derivative_le_of_energyNeutral (s : Finset ι)
    (x eps alpha eta : ι → ℝ)
    {H0 H1 H2 T0 T1 G Hd0 Hd1 nu mu Cov Vvar : ℝ}
    (hH0 : H0 = ∑ i ∈ s, eps i)
    (hH1 : H1 = ∑ i ∈ s, x i * eps i)
    (hH2 : H2 = ∑ i ∈ s, x i ^ 2 * eps i)
    (hT0 : T0 = ∑ i ∈ s, alpha i)
    (hT1 : T1 = ∑ i ∈ s, x i * alpha i)
    (hG : G = ∑ i ∈ s, eta i)
    (hH0pos : 0 < H0) (hH1pos : 0 < H1) (hnu : 0 < nu)
    (heps : ∀ i ∈ s, 0 ≤ eps i) (heta : ∀ i ∈ s, 0 ≤ eta i)
    (hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i)
    (hmu : mu = H1 / H0)
    (hneutral : T0 = 0)
    (hCov : Cov = T1 / H0)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0)
    (hHd0 : Hd0 = 2 * T0 - 2 * nu * H1)
    (hHd1 : Hd1 = 2 * T1 - 2 * nu * H2) :
    Hd1 / H1 - Hd0 / H0 ≤ G / (2 * nu * H1) := by
  refine log_bandwidth_derivative_le (Cov := Cov) s x eps alpha eta hH0 hH1 hH2
    hT0 hT1 hG hH0pos hH1pos hnu heps heta hCS hmu ?_ hVvar hHd0 hHd1
  rw [centred_covariance_eq s x alpha hT0 hT1, hneutral, hCov]
  ring

/-! ## (S8) Identity (I.4): the exact gap decomposition -/

/-- The Cauchy–Schwarz part of the defect. -/
def gapCS (mu Cov Vvar g : ℝ) : ℝ := 2 / mu * (Real.sqrt (Vvar * g) - Cov)

/-- The square-completion part of the defect. -/
def gapSC (mu nu Vvar g : ℝ) : ℝ :=
  2 * nu / mu * (Real.sqrt Vvar - Real.sqrt g / (2 * nu)) ^ 2

/-- **(S8)(ii) — identity (I.4).**  The exact decomposition of the gap between
the (I.3) bound and the (I.1) value.  The radicals cancel: this is an
*identity*, not an inequality, and it holds without any Cauchy–Schwarz input. -/
theorem spectral_front_gap_decomposition
    {H0 H1 H2 T0 T1 G Hd0 Hd1 nu mu Cov Vvar : ℝ}
    (hH0pos : 0 < H0) (hH1pos : 0 < H1) (hnu : 0 < nu)
    (hVnn : 0 ≤ Vvar) (hGnn : 0 ≤ G)
    (hmu : mu = H1 / H0)
    (hCov : Cov = (T1 - mu * T0) / H0)
    (hvar : H0 * H2 - H1 ^ 2 = H0 ^ 2 * Vvar)
    (hHd0 : Hd0 = 2 * T0 - 2 * nu * H1)
    (hHd1 : Hd1 = 2 * T1 - 2 * nu * H2) :
    G / (2 * nu * H1) - (Hd1 / H1 - Hd0 / H0)
      = gapCS mu Cov Vvar (G / H0) + gapSC mu nu Vvar (G / H0) := by
  have hH0ne : H0 ≠ 0 := ne_of_gt hH0pos
  have hnune : nu ≠ 0 := ne_of_gt hnu
  have hmupos : 0 < mu := by rw [hmu]; exact div_pos hH1pos hH0pos
  have hmune : mu ≠ 0 := ne_of_gt hmupos
  have hgnn : (0 : ℝ) ≤ G / H0 := div_nonneg hGnn hH0pos.le
  have hident : Hd1 / H1 - Hd0 / H0 = 2 * Cov / mu - 2 * nu * Vvar / mu :=
    log_bandwidth_derivative_identity hH0pos hH1pos hmu hCov hvar hHd0 hHd1
  obtain ⟨A, hA0, hAV⟩ : ∃ A : ℝ, 0 ≤ A ∧ A ^ 2 = Vvar :=
    ⟨Real.sqrt Vvar, Real.sqrt_nonneg _, Real.sq_sqrt hVnn⟩
  obtain ⟨B, hB0, hBG⟩ : ∃ B : ℝ, 0 ≤ B ∧ B ^ 2 = G / H0 :=
    ⟨Real.sqrt (G / H0), Real.sqrt_nonneg _, Real.sq_sqrt hgnn⟩
  have hsA : Real.sqrt Vvar = A := by rw [← hAV, Real.sqrt_sq hA0]
  have hsB : Real.sqrt (G / H0) = B := by rw [← hBG, Real.sqrt_sq hB0]
  have hsAB : Real.sqrt (Vvar * (G / H0)) = A * B := by
    rw [Real.sqrt_mul hVnn, hsA, hsB]
  have hGB : G = B ^ 2 * H0 := by
    rw [hBG]; field_simp
  have hH1eq : H1 = mu * H0 := by rw [hmu]; field_simp
  rw [gapCS, gapSC, hident, hsAB, hsA, hsB, ← hAV, hGB, hH1eq]
  field_simp
  ring

/-- **(S8)(i), square-completion part.**  `0 ≤ GammaSC`. -/
theorem spectral_front_gapSC_nonneg {nu mu Vvar g : ℝ}
    (hmupos : 0 < mu) (hnu : 0 < nu) :
    0 ≤ gapSC mu nu Vvar g := by
  rw [gapSC]
  positivity

/-- **(S8)(i), Cauchy–Schwarz part.**  `0 ≤ GammaCS`, which is exactly the
content of (I.2). -/
theorem spectral_front_gapCS_nonneg (s : Finset ι) (x eps alpha eta : ι → ℝ)
    {H0 G mu Cov Vvar : ℝ}
    (hH0pos : 0 < H0) (hmupos : 0 < mu)
    (hG : G = ∑ i ∈ s, eta i)
    (hCov : Cov = (∑ i ∈ s, (x i - mu) * alpha i) / H0)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0)
    (heps : ∀ i ∈ s, 0 ≤ eps i) (heta : ∀ i ∈ s, 0 ≤ eta i)
    (hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i) :
    0 ≤ gapCS mu Cov Vvar (G / H0) := by
  have habs := abs_covariance_le_sqrt s x eps alpha eta hH0pos hG hCov hVvar
    heps heta hCS
  have hrw : Vvar * G / H0 = Vvar * (G / H0) := by ring
  rw [hrw] at habs
  have hCovle : Cov ≤ Real.sqrt (Vvar * (G / H0)) :=
    le_trans (le_abs_self Cov) habs
  rw [gapCS]
  have h2 : (0 : ℝ) ≤ 2 / mu := by positivity
  exact mul_nonneg h2 (by linarith)

/-- **(S8), packaged.**  Both defect parts are nonnegative and their sum is
exactly the gap.

Scope, stated precisely.  Together with (S7) this is the *pointwise* content of
Lemma 4 of the paper: the differential inequality at one instant, with `Hd0`
and `Hd1` supplied as plain reals rather than as derivatives.  Lemma 4's actual
conclusion — that the monotone `Λ` is non-increasing, which is the paper's
headline object — additionally requires integrating that inequality in time,
and **no time integration is formalized anywhere in this file**. -/
theorem spectral_front_defect_decomposition (s : Finset ι)
    (x eps alpha eta : ι → ℝ)
    {H0 H1 H2 T0 T1 G Hd0 Hd1 nu mu Cov Vvar : ℝ}
    (hH0 : H0 = ∑ i ∈ s, eps i)
    (hH1 : H1 = ∑ i ∈ s, x i * eps i)
    (hH2 : H2 = ∑ i ∈ s, x i ^ 2 * eps i)
    (hT0 : T0 = ∑ i ∈ s, alpha i)
    (hT1 : T1 = ∑ i ∈ s, x i * alpha i)
    (hG : G = ∑ i ∈ s, eta i)
    (hH0pos : 0 < H0) (hH1pos : 0 < H1) (hnu : 0 < nu)
    (heps : ∀ i ∈ s, 0 ≤ eps i) (heta : ∀ i ∈ s, 0 ≤ eta i)
    (hCS : ∀ i ∈ s, alpha i ^ 2 ≤ eps i * eta i)
    (hmu : mu = H1 / H0)
    (hCov : Cov = (∑ i ∈ s, (x i - mu) * alpha i) / H0)
    (hVvar : Vvar = (∑ i ∈ s, (x i - mu) ^ 2 * eps i) / H0)
    (hHd0 : Hd0 = 2 * T0 - 2 * nu * H1)
    (hHd1 : Hd1 = 2 * T1 - 2 * nu * H2) :
    0 ≤ gapCS mu Cov Vvar (G / H0) ∧ 0 ≤ gapSC mu nu Vvar (G / H0) ∧
      G / (2 * nu * H1) - (Hd1 / H1 - Hd0 / H0)
        = gapCS mu Cov Vvar (G / H0) + gapSC mu nu Vvar (G / H0) := by
  have hmupos : 0 < mu := by rw [hmu]; exact div_pos hH1pos hH0pos
  have hVnn : 0 ≤ Vvar := variance_nonneg s x eps hH0pos heps hVvar
  have hGnn : 0 ≤ G := by
    rw [hG]; exact Finset.sum_nonneg fun i hi => heta i hi
  have hvar : H0 * H2 - H1 ^ 2 = H0 ^ 2 * Vvar :=
    moment_gram_eq_sq_mul_variance s x eps hH0 hH1 hH2 hH0pos hmu hVvar
  have hCov' : Cov = (T1 - mu * T0) / H0 := by
    rw [hCov, centred_covariance_eq s x alpha hT0 hT1]
  refine ⟨spectral_front_gapCS_nonneg s x eps alpha eta hH0pos hmupos hG hCov
      hVvar heps heta hCS,
    spectral_front_gapSC_nonneg hmupos hnu,
    spectral_front_gap_decomposition hH0pos hH1pos hnu hVnn hGnn hmu hCov' hvar
      hHd0 hHd1⟩

/-! ## (S9) Lemma K

The finite core is the same Cauchy–Schwarz inequality as (I.2), with the
splitting `|c i| = (x i)^{-1/2} · ((x i)^{1/2} |c i|)`.  The packaged version
carries the paper's two analytic inputs as named hypotheses.
-/

/-- **(S9), finite core.**  `(Σ |c|)² ≤ (Σ x⁻¹)·(Σ x·c²)` for positive
weights `x`. -/
theorem sq_sum_abs_le_sum_inv_mul_sum_mul_sq (s : Finset ι) (c x : ι → ℝ)
    (hx : ∀ i ∈ s, 0 < x i) :
    (∑ i ∈ s, |c i|) ^ 2 ≤ (∑ i ∈ s, (x i)⁻¹) * ∑ i ∈ s, x i * c i ^ 2 := by
  refine Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul s
    (f := fun i => (x i)⁻¹) (g := fun i => x i * c i ^ 2)
    (fun i hi => inv_nonneg.2 (hx i hi).le)
    (fun i hi => mul_nonneg (hx i hi).le (sq_nonneg _)) ?_
  intro i hi
  rw [sq_abs, ← mul_assoc, inv_mul_cancel₀ (ne_of_gt (hx i hi)), one_mul]

/-- **(S9), packaged Lemma K.**  With `supNorm`, `gradNorm`, `nlNorm` standing
for `‖u‖_∞`, `‖∇u‖₂` and `‖ℙ(u·∇u)‖₂`, and with the two **unformalised
analytic inputs**

* `hsup : supNorm ≤ Σ |c i|` — Fourier inversion, and
* `hnl : nlNorm ≤ supNorm * gradNorm` — the bilinear bound together with the
  `L²`-boundedness of the Leray projection,

as named hypotheses, the front functional obeys
`nlNorm²/gradNorm⁴ ≤ Σ (x i)⁻¹`.  For `x i = |k|²` over a band the right-hand
side is the lattice sum `S_N` of the paper.  Neither `hsup` nor `hnl` is proved
here; no Fourier transform occurs in this file. -/
theorem lemmaK_bound (s : Finset ι) (c x : ι → ℝ)
    {supNorm gradNorm nlNorm : ℝ}
    (hx : ∀ i ∈ s, 0 < x i)
    (hsup : supNorm ≤ ∑ i ∈ s, |c i|)
    (hgrad : gradNorm ^ 2 = ∑ i ∈ s, x i * c i ^ 2)
    (hnl : nlNorm ≤ supNorm * gradNorm)
    (hgradpos : 0 < gradNorm) (hsupnn : 0 ≤ supNorm) (hnlnn : 0 ≤ nlNorm) :
    nlNorm ^ 2 / gradNorm ^ 4 ≤ ∑ i ∈ s, (x i)⁻¹ := by
  have hcore := sq_sum_abs_le_sum_inv_mul_sum_mul_sq s c x hx
  rw [← hgrad] at hcore
  have hSnn : 0 ≤ ∑ i ∈ s, (x i)⁻¹ :=
    Finset.sum_nonneg fun i hi => inv_nonneg.2 (hx i hi).le
  have h1 : nlNorm ^ 2 ≤ supNorm ^ 2 * gradNorm ^ 2 := by
    nlinarith [hnl, hnlnn, mul_nonneg hsupnn hgradpos.le]
  have h2 : supNorm ^ 2 ≤ (∑ i ∈ s, |c i|) ^ 2 := by
    nlinarith [hsup, hsupnn]
  have h3 : supNorm ^ 2 ≤ (∑ i ∈ s, (x i)⁻¹) * gradNorm ^ 2 := le_trans h2 hcore
  have h4 : nlNorm ^ 2 ≤ (∑ i ∈ s, (x i)⁻¹) * gradNorm ^ 4 := by
    nlinarith [h1, h3, sq_nonneg gradNorm, hgradpos]
  rw [div_le_iff₀ (by positivity)]
  exact h4

end

end NSSingularity
