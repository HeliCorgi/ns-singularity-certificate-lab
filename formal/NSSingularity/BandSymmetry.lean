/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import Mathlib

/-!
# The band symmetry lemma: finite symmetric lattice sums

Lemma 9 of `docs/paper_lambda_dichotomy/complete_proof.md` states that for a
radial weight `f(k) = φ(|k|²)` on `ℤ³ \ {0}`,

`Σ_k k_i k_j f(k) = δ_{ij} · (1/3) · Σ_k |k|² f(k)`.

This file formalizes exactly that, **over a finite index set** `s : Finset (ℤ × ℤ × ℤ)`
that is closed under the coordinate sign flips and the coordinate transpositions.

## What is formalized

* `sum_offDiagonal_12_eq_zero`, `sum_offDiagonal_13_eq_zero`,
  `sum_offDiagonal_23_eq_zero` — off-diagonal vanishing.  The proof is that the
  sign flip in one of the two coordinates is an involution of `s` which
  preserves the weight and negates the summand, so the sum equals its own
  negation.
* `sum_sq_fst_eq_sum_sq_snd`, `sum_sq_fst_eq_sum_sq_thd` — diagonal equality,
  by the coordinate-swap bijections.
* `three_mul_sum_sq_fst_eq_sum_normSq` — the identity
  `3 · Σ k₁² w(k) = Σ (k₁² + k₂² + k₃²) w(k)`, which is the engine behind the
  paper's exact family laws `H₀ = (2/3)‖v₀‖² T_N` and `H₁ = (2/3)‖v₀‖² S_N`
  (Lemma 10).
* `band_symmetry` — the three statements packaged together, which is Lemma 9
  in the form the paper uses it.

Radiality is imposed exactly as the paper does: the hypothesis is that
`w k = f (k.1² + k.2.1² + k.2.2²)` for some `f : ℤ → ℝ`.  The invariance of `w`
under the sign flips and the transpositions is **derived** from that, not
assumed separately.

## What is deliberately **not** proved

* All sums here are over a **finite** `Finset`.  No lattice-point asymptotics,
  no absolute-summability argument for the infinite lattice `ℤ³ \ {0}`, and no
  rearrangement theorem is used or claimed.  The paper's Lemma 9 is stated for
  an absolutely summable weight on the infinite lattice; the passage from the
  finite symmetric set to that statement is not formalized.
* Nothing here mentions Fourier series, the Leray projection, or any
  Navier–Stokes quantity.  The connection to `H₀` and `H₁` is made in the
  paper, not in Lean.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/

namespace NSSingularity

noncomputable section

/-! ## Two reindexing helpers -/

/-- Reindexing a finite sum along an involution of the index set. -/
theorem sum_comp_involution {ι : Type*} {s : Finset ι} (F G : ι → ℝ)
    (σ : ι → ι) (hmem : ∀ k ∈ s, σ k ∈ s) (hinv : ∀ k ∈ s, σ (σ k) = k)
    (hFG : ∀ k ∈ s, F k = G (σ k)) :
    ∑ k ∈ s, F k = ∑ k ∈ s, G k :=
  Finset.sum_nbij' (i := σ) (j := σ) hmem hmem hinv hinv hFG

/-- A finite sum whose summand is negated by an involution of the index set
vanishes. -/
theorem sum_eq_zero_of_neg_involution {ι : Type*} {s : Finset ι} (F : ι → ℝ)
    (σ : ι → ι) (hmem : ∀ k ∈ s, σ k ∈ s) (hinv : ∀ k ∈ s, σ (σ k) = k)
    (hneg : ∀ k ∈ s, F (σ k) = -F k) :
    ∑ k ∈ s, F k = 0 := by
  have h1 : ∑ k ∈ s, F k = ∑ k ∈ s, (fun k => -F k) k := by
    refine sum_comp_involution F (fun k => -F k) σ hmem hinv ?_
    intro k hk
    rw [hneg k hk, neg_neg]
  simp only [Finset.sum_neg_distrib] at h1
  linarith

/-! ## The three coordinate involutions -/

/-- Sign flip in the first coordinate. -/
def flipFst (k : ℤ × ℤ × ℤ) : ℤ × ℤ × ℤ := (-k.1, k.2.1, k.2.2)

/-- Sign flip in the second coordinate. -/
def flipSnd (k : ℤ × ℤ × ℤ) : ℤ × ℤ × ℤ := (k.1, -k.2.1, k.2.2)

/-- Transposition of the first two coordinates. -/
def swapFstSnd (k : ℤ × ℤ × ℤ) : ℤ × ℤ × ℤ := (k.2.1, k.1, k.2.2)

/-- Transposition of the first and third coordinates. -/
def swapFstThd (k : ℤ × ℤ × ℤ) : ℤ × ℤ × ℤ := (k.2.2, k.2.1, k.1)

theorem flipFst_flipFst (k : ℤ × ℤ × ℤ) : flipFst (flipFst k) = k := by
  simp [flipFst]

theorem flipSnd_flipSnd (k : ℤ × ℤ × ℤ) : flipSnd (flipSnd k) = k := by
  simp [flipSnd]

theorem swapFstSnd_swapFstSnd (k : ℤ × ℤ × ℤ) :
    swapFstSnd (swapFstSnd k) = k := by
  simp [swapFstSnd]

theorem swapFstThd_swapFstThd (k : ℤ × ℤ × ℤ) :
    swapFstThd (swapFstThd k) = k := by
  simp [swapFstThd]

/-! ## Radiality: the weight is invariant under all four maps -/

/-- A radial weight is invariant under the sign flip in the first
coordinate. -/
theorem radial_flipFst {s : Finset (ℤ × ℤ × ℤ)} {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2))
    (hflip : ∀ k ∈ s, flipFst k ∈ s) {k : ℤ × ℤ × ℤ} (hk : k ∈ s) :
    w (flipFst k) = w k := by
  rw [hradial _ (hflip k hk), hradial k hk]
  simp only [flipFst]
  congr 1
  ring

/-- A radial weight is invariant under the sign flip in the second
coordinate. -/
theorem radial_flipSnd {s : Finset (ℤ × ℤ × ℤ)} {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2))
    (hflip : ∀ k ∈ s, flipSnd k ∈ s) {k : ℤ × ℤ × ℤ} (hk : k ∈ s) :
    w (flipSnd k) = w k := by
  rw [hradial _ (hflip k hk), hradial k hk]
  simp only [flipSnd]
  congr 1
  ring

/-- A radial weight is invariant under the transposition of the first two
coordinates. -/
theorem radial_swapFstSnd {s : Finset (ℤ × ℤ × ℤ)} {w : ℤ × ℤ × ℤ → ℝ}
    {f : ℤ → ℝ}
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2))
    (hswap : ∀ k ∈ s, swapFstSnd k ∈ s) {k : ℤ × ℤ × ℤ} (hk : k ∈ s) :
    w (swapFstSnd k) = w k := by
  rw [hradial _ (hswap k hk), hradial k hk]
  simp only [swapFstSnd]
  congr 1
  ring

/-- A radial weight is invariant under the transposition of the first and third
coordinates. -/
theorem radial_swapFstThd {s : Finset (ℤ × ℤ × ℤ)} {w : ℤ × ℤ × ℤ → ℝ}
    {f : ℤ → ℝ}
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2))
    (hswap : ∀ k ∈ s, swapFstThd k ∈ s) {k : ℤ × ℤ × ℤ} (hk : k ∈ s) :
    w (swapFstThd k) = w k := by
  rw [hradial _ (hswap k hk), hradial k hk]
  simp only [swapFstThd]
  congr 1
  ring

/-! ## (a) Off-diagonal vanishing -/

/-- **Lemma 9(a), the `(1,2)` entry.**  `Σ_{k ∈ s} k₁ k₂ w(k) = 0`.

The sign flip in the first coordinate is an involution of `s` that preserves a
radial `w` and negates the summand. -/
theorem sum_offDiagonal_12_eq_zero {s : Finset (ℤ × ℤ × ℤ)}
    {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hflip : ∀ k ∈ s, flipFst k ∈ s)
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2)) :
    ∑ k ∈ s, (k.1 : ℝ) * (k.2.1 : ℝ) * w k = 0 := by
  refine sum_eq_zero_of_neg_involution _ flipFst hflip
    (fun k _ => flipFst_flipFst k) ?_
  intro k hk
  have hw : w (flipFst k) = w k := radial_flipFst hradial hflip hk
  simp only [flipFst] at hw ⊢
  rw [hw]
  push_cast
  ring

/-- **Lemma 9(a), the `(1,3)` entry.**  `Σ_{k ∈ s} k₁ k₃ w(k) = 0`. -/
theorem sum_offDiagonal_13_eq_zero {s : Finset (ℤ × ℤ × ℤ)}
    {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hflip : ∀ k ∈ s, flipFst k ∈ s)
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2)) :
    ∑ k ∈ s, (k.1 : ℝ) * (k.2.2 : ℝ) * w k = 0 := by
  refine sum_eq_zero_of_neg_involution _ flipFst hflip
    (fun k _ => flipFst_flipFst k) ?_
  intro k hk
  have hw : w (flipFst k) = w k := radial_flipFst hradial hflip hk
  simp only [flipFst] at hw ⊢
  rw [hw]
  push_cast
  ring

/-- **Lemma 9(a), the `(2,3)` entry.**  `Σ_{k ∈ s} k₂ k₃ w(k) = 0`. -/
theorem sum_offDiagonal_23_eq_zero {s : Finset (ℤ × ℤ × ℤ)}
    {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hflip : ∀ k ∈ s, flipSnd k ∈ s)
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2)) :
    ∑ k ∈ s, (k.2.1 : ℝ) * (k.2.2 : ℝ) * w k = 0 := by
  refine sum_eq_zero_of_neg_involution _ flipSnd hflip
    (fun k _ => flipSnd_flipSnd k) ?_
  intro k hk
  have hw : w (flipSnd k) = w k := radial_flipSnd hradial hflip hk
  simp only [flipSnd] at hw ⊢
  rw [hw]
  push_cast
  ring

/-! ## (b) Diagonal equality -/

/-- **Lemma 9(b).**  `Σ k₁² w(k) = Σ k₂² w(k)`, by the transposition of the
first two coordinates. -/
theorem sum_sq_fst_eq_sum_sq_snd {s : Finset (ℤ × ℤ × ℤ)}
    {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hswap : ∀ k ∈ s, swapFstSnd k ∈ s)
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2)) :
    ∑ k ∈ s, (k.1 : ℝ) ^ 2 * w k = ∑ k ∈ s, (k.2.1 : ℝ) ^ 2 * w k := by
  refine sum_comp_involution _ _ swapFstSnd hswap
    (fun k _ => swapFstSnd_swapFstSnd k) ?_
  intro k hk
  have hw : w (swapFstSnd k) = w k := radial_swapFstSnd hradial hswap hk
  simp only [swapFstSnd] at hw ⊢
  rw [hw]

/-- **Lemma 9(b).**  `Σ k₁² w(k) = Σ k₃² w(k)`, by the transposition of the
first and third coordinates. -/
theorem sum_sq_fst_eq_sum_sq_thd {s : Finset (ℤ × ℤ × ℤ)}
    {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hswap : ∀ k ∈ s, swapFstThd k ∈ s)
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2)) :
    ∑ k ∈ s, (k.1 : ℝ) ^ 2 * w k = ∑ k ∈ s, (k.2.2 : ℝ) ^ 2 * w k := by
  refine sum_comp_involution _ _ swapFstThd hswap
    (fun k _ => swapFstThd_swapFstThd k) ?_
  intro k hk
  have hw : w (swapFstThd k) = w k := radial_swapFstThd hradial hswap hk
  simp only [swapFstThd] at hw ⊢
  rw [hw]

/-! ## (c) The third identity -/

/-- **Lemma 9(c).**  `3 · Σ k₁² w(k) = Σ (k₁² + k₂² + k₃²) w(k)`.

This is the identity behind the paper's exact family laws
`H₀ = (2/3)‖v₀‖² T_N` and `H₁ = (2/3)‖v₀‖² S_N` (Lemma 10).  The sum runs over
a **finite** symmetric set; no lattice-point asymptotics are claimed. -/
theorem three_mul_sum_sq_fst_eq_sum_normSq {s : Finset (ℤ × ℤ × ℤ)}
    {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hswap12 : ∀ k ∈ s, swapFstSnd k ∈ s)
    (hswap13 : ∀ k ∈ s, swapFstThd k ∈ s)
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2)) :
    3 * ∑ k ∈ s, (k.1 : ℝ) ^ 2 * w k
      = ∑ k ∈ s, ((k.1 : ℝ) ^ 2 + (k.2.1 : ℝ) ^ 2 + (k.2.2 : ℝ) ^ 2) * w k := by
  have h12 := sum_sq_fst_eq_sum_sq_snd hswap12 hradial
  have h13 := sum_sq_fst_eq_sum_sq_thd hswap13 hradial
  have hsplit : ∑ k ∈ s, ((k.1 : ℝ) ^ 2 + (k.2.1 : ℝ) ^ 2 + (k.2.2 : ℝ) ^ 2) * w k
      = (∑ k ∈ s, (k.1 : ℝ) ^ 2 * w k) + (∑ k ∈ s, (k.2.1 : ℝ) ^ 2 * w k)
        + ∑ k ∈ s, (k.2.2 : ℝ) ^ 2 * w k := by
    rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun k _ => by ring
  rw [hsplit, ← h12, ← h13]
  ring

/-- **Lemma 9, packaged.**  For a finite index set closed under the sign flips
in the first two coordinates and under the two transpositions, and a radial
weight, the off-diagonal moments vanish, the diagonal moments agree, and three
times a diagonal moment is the full squared-norm moment. -/
theorem band_symmetry {s : Finset (ℤ × ℤ × ℤ)} {w : ℤ × ℤ × ℤ → ℝ} {f : ℤ → ℝ}
    (hflip1 : ∀ k ∈ s, flipFst k ∈ s)
    (hflip2 : ∀ k ∈ s, flipSnd k ∈ s)
    (hswap12 : ∀ k ∈ s, swapFstSnd k ∈ s)
    (hswap13 : ∀ k ∈ s, swapFstThd k ∈ s)
    (hradial : ∀ k ∈ s, w k = f (k.1 ^ 2 + k.2.1 ^ 2 + k.2.2 ^ 2)) :
    (∑ k ∈ s, (k.1 : ℝ) * (k.2.1 : ℝ) * w k = 0 ∧
      ∑ k ∈ s, (k.1 : ℝ) * (k.2.2 : ℝ) * w k = 0 ∧
      ∑ k ∈ s, (k.2.1 : ℝ) * (k.2.2 : ℝ) * w k = 0) ∧
    (∑ k ∈ s, (k.1 : ℝ) ^ 2 * w k = ∑ k ∈ s, (k.2.1 : ℝ) ^ 2 * w k ∧
      ∑ k ∈ s, (k.1 : ℝ) ^ 2 * w k = ∑ k ∈ s, (k.2.2 : ℝ) ^ 2 * w k) ∧
    3 * ∑ k ∈ s, (k.1 : ℝ) ^ 2 * w k
      = ∑ k ∈ s, ((k.1 : ℝ) ^ 2 + (k.2.1 : ℝ) ^ 2 + (k.2.2 : ℝ) ^ 2) * w k :=
  ⟨⟨sum_offDiagonal_12_eq_zero hflip1 hradial,
    sum_offDiagonal_13_eq_zero hflip1 hradial,
    sum_offDiagonal_23_eq_zero hflip2 hradial⟩,
   ⟨sum_sq_fst_eq_sum_sq_snd hswap12 hradial,
    sum_sq_fst_eq_sum_sq_thd hswap13 hradial⟩,
   three_mul_sum_sq_fst_eq_sum_normSq hswap12 hswap13 hradial⟩

end

end NSSingularity
