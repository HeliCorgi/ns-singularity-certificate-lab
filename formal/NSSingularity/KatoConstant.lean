/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under the Apache 2.0 licence.
-/
import Mathlib
import NSSingularity.GaussianTransfer

/-!
# The `n = 3` Kato constant — the finite algebra

`docs/research_notes/kato_h3_constants.md` derives, in the repository's own
normalisation, a computable constant `G3` with

  `|⟨P(v·∇w), w⟩_{Ḣ³}| ≤ G3 ‖v‖_{Ḣ³} ‖w‖_{Ḣ³}²`

for mean-zero divergence-free fields on `T³`, and its section 9 lists the
finitely many elementary inequalities the derivation consumes.  This file is
that list, verbatim, as unconditional theorems about real numbers and vectors:

* `cube_diff_bound` — `|x³ - y³| ≤ 3 |x - y| (max x y)²` for `0 ≤ x, y`, the
  scalar brick behind the cubic-difference continuity step; delegated to the
  existing `abs_cube_sub_cube_le` of `GaussianTransfer.lean`, restated in the
  factor order the note uses.
* `am_gm_split` — the AM–GM split `p² j² ≤ (p j³ + p³ j)/2`, which is how the
  derivation trades the mixed weight `|p|²|j|²` for the two extreme weights
  that the lattice sums `A₄`, `A₆` can absorb.
* `shifted_ratio_bound` — `‖m + j‖³ ≤ (1 + ‖j‖)³ ‖m‖³` for `‖m‖ ≥ 1`: on the
  mean-zero torus every retained wavenumber has `|m| ≥ 1`, and this is the
  `sup_{|m| ≥ 1}` ratio bound that controls the weight transported from the
  shifted mode `m + j` back to `m`.
* `inv_pow_four_succ_le_telescope`, `inv_pow_tail_bound` — the shell-tail
  comparison: termwise `1/(m+1)⁴ ≤ (1/m³ - 1/(m+1)³)/3`, hence
  `Σ_{m=N+1}^{M} m⁻⁴ ≤ 1/(3N³)` for every finite `M`.  This is the tail bound
  the certificate charges against the exactly-summed head of the lattice sum
  `A₄` (and, with the shell count `26 m²`, of `A₆`).
* `g3_assembly_mono`, `g3_of_a4` — the checker's final assembly step
  `G3_upper = 12 · sqrt(A4_upper)` is monotone in the enclosure of `A₄`, so an
  upper bound on `A₄` yields an upper bound on `G3`.  Trivial, and stated so
  that the monotonicity appears in the axiom audit as a named theorem.

## What is not formalised, and why

The commutator estimate that consumes these bricks — the passage from the
trilinear form on `ℓ²(ℤ³)` Fourier coefficients to the constant `G3` — is
infinite-dimensional analysis.  Per the repository's rules it is **not**
axiomatised in Lean; it lives in the research note, and the audit trail for it
is the note plus the interval-arithmetic checker
(`src/ns_certificate_lab/kato_constant.py`).  Nothing here mentions
Navier–Stokes, no PDE statement is formalised, and nothing in this file, alone
or in combination with the rest of the library, is a statement about the Clay
problem.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/

namespace NSSingularity

noncomputable section

open Set

/-! ## K1 — the cubic-difference brick -/

/-- **K1.**  `|x³ - y³| ≤ 3 |x - y| (max x y)²` for nonnegative `x`, `y` —
the note's cubic-difference continuity brick, delegated to the existing
`abs_cube_sub_cube_le` and restated in the factor order of
`kato_h3_constants.md` §9. -/
theorem cube_diff_bound {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) :
    |x ^ 3 - y ^ 3| ≤ 3 * |x - y| * max x y ^ 2 :=
  calc |x ^ 3 - y ^ 3| ≤ 3 * max x y ^ 2 * |x - y| := abs_cube_sub_cube_le hx hy
    _ = 3 * |x - y| * max x y ^ 2 := by ring

/-! ## K2 — the AM–GM split -/

/-- **K2.**  The AM–GM split `p² j² ≤ (p j³ + p³ j)/2` for `0 ≤ p, j`: the
mixed weight is dominated by the mean of the two extreme weights, which is the
step that lets the derivation replace `|p|²|j|²` by lattice sums it can
actually evaluate.  The whole proof is `p j (p - j)² ≥ 0`. -/
theorem am_gm_split {p j : ℝ} (hp : 0 ≤ p) (hj : 0 ≤ j) :
    p ^ 2 * j ^ 2 ≤ (p * j ^ 3 + p ^ 3 * j) / 2 := by
  nlinarith [mul_nonneg (mul_nonneg hp hj) (sq_nonneg (p - j))]

/-! ## K3 — the shifted-ratio bound -/

/-- **K3.**  For `‖m‖ ≥ 1` (every retained wavenumber of the mean-zero torus),
`‖m + j‖³ ≤ (1 + ‖j‖)³ ‖m‖³`.  This is the `sup_{|m| ≥ 1}` ratio bound of the
note: the triangle inequality gives `‖m + j‖ ≤ ‖m‖ + ‖j‖ ≤ (1 + ‖j‖) ‖m‖`
because `‖j‖ · 1 ≤ ‖j‖ · ‖m‖`, and cubing preserves the inequality.  Stated in
an arbitrary normed abelian group since nothing about `ℤ³` enters. -/
theorem shifted_ratio_bound {E : Type*} [NormedAddCommGroup E] {m j : E}
    (hm : 1 ≤ ‖m‖) :
    ‖m + j‖ ^ 3 ≤ (1 + ‖j‖) ^ 3 * ‖m‖ ^ 3 := by
  have hj : (0 : ℝ) ≤ ‖j‖ := norm_nonneg j
  have h1 : ‖m + j‖ ≤ (1 + ‖j‖) * ‖m‖ := by
    have htri : ‖m + j‖ ≤ ‖m‖ + ‖j‖ := norm_add_le m j
    nlinarith [mul_nonneg hj (by linarith : (0 : ℝ) ≤ ‖m‖ - 1)]
  calc ‖m + j‖ ^ 3 ≤ ((1 + ‖j‖) * ‖m‖) ^ 3 :=
        pow_le_pow_left₀ (norm_nonneg _) h1 3
    _ = (1 + ‖j‖) ^ 3 * ‖m‖ ^ 3 := mul_pow _ _ _

/-! ## K4 — the inverse-fourth-power tail -/

/-- **K4, one step.**  The telescoping brick: for `1 ≤ m`,
`1/(m+1)⁴ ≤ 1/(3m³) - 1/(3(m+1)³)`.  Clearing denominators, the difference of
the two sides is `3 (m+1)³ (6m² + 4m + 1)` up to a positive factor, which is
manifestly positive. -/
theorem inv_pow_four_succ_le_telescope {m : ℕ} (hm : 1 ≤ m) :
    (1 : ℝ) / ((m : ℝ) + 1) ^ 4
      ≤ 1 / (3 * (m : ℝ) ^ 3) - 1 / (3 * ((m : ℝ) + 1) ^ 3) := by
  have hx : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
  have hx0 : (0 : ℝ) < (m : ℝ) := by linarith
  have hx1 : (0 : ℝ) < (m : ℝ) + 1 := by linarith
  rw [div_sub_div _ _ (by positivity) (by positivity),
    div_le_div_iff₀ (by positivity) (by positivity)]
  have hprod : (0 : ℝ) < ((m : ℝ) + 1) ^ 3 * (6 * (m : ℝ) ^ 2 + 4 * (m : ℝ) + 1) := by
    positivity
  nlinarith [hprod]

/-- **K4, the tail.**  For `1 ≤ N` and every finite `M`,
`Σ_{m = N+1}^{M} 1/m⁴ ≤ 1/(3N³)`: summing the telescoping brick collapses the
partial sum to `1/(3N³) - 1/(3M³)` and the subtracted term is dropped.  This
is the finite-sum form the interval checker uses — the head of the lattice sum
`A₄` is summed exactly and every finite continuation of the tail is bounded by
`1/(3N³)`, uniformly in `M`. -/
theorem inv_pow_tail_bound {N : ℕ} (hN : 1 ≤ N) (M : ℕ) :
    ∑ m ∈ Finset.Icc (N + 1) M, (1 : ℝ) / (m : ℝ) ^ 4
      ≤ 1 / (3 * (N : ℝ) ^ 3) := by
  rcases le_or_gt N M with hNM | hMN
  · have key : ∀ K : ℕ, N ≤ K →
        ∑ m ∈ Finset.Icc (N + 1) K, (1 : ℝ) / (m : ℝ) ^ 4
          ≤ 1 / (3 * (N : ℝ) ^ 3) - 1 / (3 * (K : ℝ) ^ 3) := by
      intro K hK
      induction K, hK using Nat.le_induction with
      | base =>
        rw [Finset.Icc_eq_empty (by omega), Finset.sum_empty]
        simp
      | succ K hK ih =>
        rw [Finset.sum_Icc_succ_top (by omega : N + 1 ≤ K + 1)]
        have hK1 : 1 ≤ K := le_trans hN hK
        have htel := inv_pow_four_succ_le_telescope hK1
        have hcast : ((K + 1 : ℕ) : ℝ) = (K : ℝ) + 1 := by push_cast; ring
        rw [hcast]
        linarith
    have hdrop : (0 : ℝ) ≤ 1 / (3 * (M : ℝ) ^ 3) := by positivity
    linarith [key M hNM]
  · rw [Finset.Icc_eq_empty (by omega), Finset.sum_empty]
    positivity

/-! ## K5 — the assembly step -/

/-- **K5 (i).**  The checker's final assembly `G3_upper = 12 · sqrt(A4_upper)`
is monotone in the enclosure: `a ≤ a'` gives `12 √a ≤ 12 √a'`.  The
nonnegativity hypothesis mirrors the checker's precondition on the lattice sum
(mathlib's `Real.sqrt` does not need it, since `√` of a negative number is
`0`), and is deliberately kept in the statement. -/
theorem g3_assembly_mono {a a' : ℝ} (_ha : 0 ≤ a) (h : a ≤ a') :
    12 * Real.sqrt a ≤ 12 * Real.sqrt a' := by
  have := Real.sqrt_le_sqrt h
  linarith

/-- **K5 (ii).**  The assembly inequality in the exact shape the checker
emits: an upper enclosure `A₄ ≤ bound` for the lattice sum yields the upper
enclosure `12 √A₄ ≤ 12 √bound` for the assembled constant.  This is the whole
machine-checked content of the step `G3_upper = 12 · sqrt_upper(A4_upper)`;
the value of `bound` is produced by interval arithmetic outside Lean. -/
theorem g3_of_a4 {A4 bound : ℝ} (h : A4 ≤ bound) :
    12 * Real.sqrt A4 ≤ 12 * Real.sqrt bound := by
  have := Real.sqrt_le_sqrt h
  linarith

end

end NSSingularity
