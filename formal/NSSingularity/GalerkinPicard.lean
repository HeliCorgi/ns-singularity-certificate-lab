/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import Mathlib
import NSSingularity.TimeDependentGalerkin

/-!
# The Galerkin Picard–Lindelöf lemma — the finite-dimensional core of EXT-P1

The Track-P dependency table names three external classical theorems.  EXT-P1
is local existence and uniqueness in `H^4` for periodic Navier–Stokes
(Kato/Temam); it is a PDE statement, it is recorded faithfully in the research
notes and payloads, and it is **not** proved here — encoding it as a Lean
axiom is forbidden in this repository, and rightly so.  What *is* proved here
is the finite-dimensional half of the same dependency row: the Galerkin
system that the certificate integrates is an ODE `u' = A u + B u u` in a
finite-dimensional (here: arbitrary Banach) space, with `A` the projected
dissipative linear part and `B` the projected advective bilinear part, and
for *that* system existence and uniqueness are unconditional theorems, proved
below through mathlib's Picard–Lindelöf machinery.  The Python Picard box
test (`prove_rough_enclosure` and the `torus_chain` slab loop) instantiates
these lemmas concretely with interval arithmetic.  In the project's fixed
wording: the finite-dimensional Galerkin enclosure and the scalar control ODE
are verified unconditionally, and everything a payload claims about the PDE
itself is a conditional PDE certificate assuming EXT-P1/P2/P3.  Nothing in
this file is a statement about the Clay problem.

## What is formalised

* `quadraticField` — the autonomous field `F u = A u + B u u`.
* `quadratic_field_lipschitzOnWith` — `F` is Lipschitz on the closed ball of
  radius `r` about `x₀` with the explicit constant
  `‖A‖ + 2 ‖B‖ (‖x₀‖ + r)`: the linear part contributes its operator norm
  globally, and the quadratic part is only *locally* Lipschitz, through the
  telescoping `B x x - B y y = B x (x - y) + B (x - y) y`, which is why the
  ball radius must appear in the constant and a global Picard–Lindelöf
  argument is impossible for a genuinely quadratic field.
* `quadratic_ode_local_solution` — through every datum `(t₀, x₀)` there is a
  local solution on `[t₀ - ε, t₀ + ε]` with `ε > 0` exhibited explicitly from
  the operator norms (`ε = 1/(L+1)` for the sup bound `L` of the field on the
  unit ball about the datum).  Existence, so `CompleteSpace` is required.
* `quadratic_ode_local_solution_hasDerivAt` — the same solution restated with
  a genuine two-sided derivative at every interior time, so "solution on a
  nontrivial open interval" holds in the strongest pointwise sense.
* `quadratic_ode_unique` — two solutions with the same datum that stay in a
  common ball agree on the whole interval, by the Lipschitz bound and
  mathlib's Grönwall-based uniqueness.

## What is *not* formalised, and why

* **EXT-P1 itself.**  The passage from the Galerkin family to the PDE — the
  compactness/limit argument of Kato/Temam in `H^4`, uniqueness of the strong
  PDE solution, and the continuation criterion (EXT-P3) — needs Sobolev
  spaces on `T³` and the Fourier isometry, which mathlib does not yet have in
  usable form.  These remain named external theorems; no Lean statement below
  mentions Navier–Stokes, and none is weakened to pretend otherwise.
* **The specific `A`, `B` of a mode set.**  That the concrete projected
  operators of a Galerkin family satisfy the operator-norm bounds fed into
  these lemmas is finite arithmetic performed by the Python layer.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/

namespace NSSingularity

noncomputable section

open Metric Set
open scoped NNReal

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The autonomous quadratic Galerkin field `F u = A u + B u u`: `A` is the
projected linear (dissipative) part, `B` the projected bilinear (advective)
part.  The time-dependent generalisation is `galerkinField` in
`TimeDependentGalerkin.lean`; this is its autonomous specialisation, which is
what the torus chain integrates on each slab. -/
def quadraticField (A : E →L[ℝ] E) (B : E →L[ℝ] E →L[ℝ] E) : E → E :=
  fun u => A u + B u u

@[simp] theorem quadraticField_apply (A : E →L[ℝ] E) (B : E →L[ℝ] E →L[ℝ] E) (u : E) :
    quadraticField A B u = A u + B u u := rfl

omit [NormedSpace ℝ E] in
/-- Points of a closed ball are bounded by centre plus radius. -/
private theorem norm_le_of_mem_closedBall' {x x₀ : E} {r : ℝ}
    (hx : x ∈ closedBall x₀ r) : ‖x‖ ≤ ‖x₀‖ + r := by
  have h : ‖x - x₀‖ ≤ r := by simpa [dist_eq_norm] using hx
  calc ‖x‖ = ‖x₀ + (x - x₀)‖ := by rw [add_sub_cancel]
    _ ≤ ‖x₀‖ + ‖x - x₀‖ := norm_add_le _ _
    _ ≤ ‖x₀‖ + r := by linarith

/-- The telescoping bound for the quadratic term on a ball: the constant is
`2‖B‖` times the outer radius of the ball and cannot be made global. -/
private theorem quadratic_term_bound (B : E →L[ℝ] E →L[ℝ] E) {x y x₀ : E} {r : ℝ}
    (hx : x ∈ closedBall x₀ r) (hy : y ∈ closedBall x₀ r) :
    ‖B x x - B y y‖ ≤ 2 * ‖B‖ * (‖x₀‖ + r) * ‖x - y‖ := by
  have hxn : ‖x‖ ≤ ‖x₀‖ + r := norm_le_of_mem_closedBall' hx
  have hyn : ‖y‖ ≤ ‖x₀‖ + r := norm_le_of_mem_closedBall' hy
  have hsplit : B x x - B y y = B x (x - y) + B (x - y) y := by
    have e1 : B x (x - y) = B x x - B x y := map_sub _ _ _
    have e2 : B (x - y) y = B x y - B y y := by rw [map_sub]; simp
    rw [e1, e2]; abel
  have h1 : ‖B x (x - y)‖ ≤ ‖B‖ * ‖x‖ * ‖x - y‖ := B.le_opNorm₂ x (x - y)
  have h2 : ‖B (x - y) y‖ ≤ ‖B‖ * ‖x - y‖ * ‖y‖ := B.le_opNorm₂ (x - y) y
  have hB : (0 : ℝ) ≤ ‖B‖ := ContinuousLinearMap.opNorm_nonneg B
  have hd : (0 : ℝ) ≤ ‖x - y‖ := norm_nonneg _
  have h1' : ‖B x (x - y)‖ ≤ ‖B‖ * (‖x₀‖ + r) * ‖x - y‖ := by
    refine h1.trans ?_
    exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hxn hB) hd
  have h2' : ‖B (x - y) y‖ ≤ ‖B‖ * (‖x₀‖ + r) * ‖x - y‖ := by
    refine h2.trans ?_
    calc ‖B‖ * ‖x - y‖ * ‖y‖
        ≤ ‖B‖ * ‖x - y‖ * (‖x₀‖ + r) :=
          mul_le_mul_of_nonneg_left hyn (mul_nonneg hB hd)
      _ = ‖B‖ * (‖x₀‖ + r) * ‖x - y‖ := by ring
  calc ‖B x x - B y y‖ = ‖B x (x - y) + B (x - y) y‖ := by rw [hsplit]
    _ ≤ ‖B x (x - y)‖ + ‖B (x - y) y‖ := norm_add_le _ _
    _ ≤ 2 * ‖B‖ * (‖x₀‖ + r) * ‖x - y‖ := by linarith

/-- **GP1, the local Lipschitz bound with its explicit constant.**  The
quadratic field is Lipschitz on the closed ball of radius `r` about `x₀` with
constant `‖A‖ + 2 ‖B‖ (‖x₀‖ + r)`: the linear part contributes `‖A‖`, and the
telescoped quadratic part contributes `2‖B‖` times the outer radius `‖x₀‖ + r`
of the ball.  The hypothesis `0 ≤ r` is what makes the constant nonnegative,
so that its `Real.toNNReal` coercion is lossless. -/
theorem quadratic_field_lipschitzOnWith (A : E →L[ℝ] E) (B : E →L[ℝ] E →L[ℝ] E)
    (x₀ : E) {r : ℝ} (hr : 0 ≤ r) :
    LipschitzOnWith (Real.toNNReal (‖A‖ + 2 * ‖B‖ * (‖x₀‖ + r)))
      (quadraticField A B) (closedBall x₀ r) := by
  have hxr : (0 : ℝ) ≤ ‖x₀‖ + r := add_nonneg (norm_nonneg _) hr
  have hC0 : (0 : ℝ) ≤ ‖A‖ + 2 * ‖B‖ * (‖x₀‖ + r) :=
    add_nonneg (ContinuousLinearMap.opNorm_nonneg A)
      (mul_nonneg (mul_nonneg (by norm_num) (ContinuousLinearMap.opNorm_nonneg B)) hxr)
  refine LipschitzOnWith.of_dist_le_mul fun x hx y hy => ?_
  have hquad := quadratic_term_bound B hx hy
  have hlin : ‖A x - A y‖ ≤ ‖A‖ * ‖x - y‖ := by
    rw [← map_sub]; exact A.le_opNorm _
  have hsplit : quadraticField A B x - quadraticField A B y
      = (A x - A y) + (B x x - B y y) := by
    simp only [quadraticField]; abel
  calc dist (quadraticField A B x) (quadraticField A B y)
      = ‖(A x - A y) + (B x x - B y y)‖ := by rw [dist_eq_norm, hsplit]
    _ ≤ ‖A x - A y‖ + ‖B x x - B y y‖ := norm_add_le _ _
    _ ≤ ‖A‖ * ‖x - y‖ + 2 * ‖B‖ * (‖x₀‖ + r) * ‖x - y‖ := by linarith
    _ = (‖A‖ + 2 * ‖B‖ * (‖x₀‖ + r)) * ‖x - y‖ := by ring
    _ = (Real.toNNReal (‖A‖ + 2 * ‖B‖ * (‖x₀‖ + r)) : ℝ) * dist x y := by
        rw [Real.coe_toNNReal _ hC0, dist_eq_norm]

/-- **GP2, local existence.**  Through every datum `(t₀, x₀)` the quadratic
system `u' = A u + B u u` has a local solution on `[t₀ - ε, t₀ + ε]`, with the
interval half-width exhibited explicitly: `ε = 1 / (L + 1)` where
`L = ‖B‖ (‖x₀‖ + 1)² + ‖A‖ (‖x₀‖ + 1)` is the sup bound of the field on the
unit ball about the datum.  Obtained from `galerkin_local_solution`
(`TimeDependentGalerkin.lean`) at zero forcing, i.e. from mathlib's
time-dependent `IsPicardLindelof`.  This, instantiated by the Python interval
layer for the concrete projected operators of a mode family, is the
finite-dimensional Galerkin half of the EXT-P1 dependency row; the PDE half
remains an external classical theorem and is not claimed. -/
theorem quadratic_ode_local_solution [CompleteSpace E]
    (A : E →L[ℝ] E) (B : E →L[ℝ] E →L[ℝ] E) (t₀ : ℝ) (x₀ : E) :
    ∃ ε : ℝ, 0 < ε ∧ ∃ f : ℝ → E, f t₀ = x₀ ∧
      ∀ t ∈ Icc (t₀ - ε) (t₀ + ε),
        HasDerivWithinAt f (quadraticField A B (f t)) (Icc (t₀ - ε) (t₀ + ε)) t := by
  classical
  have hx1 : (0 : ℝ) ≤ ‖x₀‖ + 1 := by have := norm_nonneg x₀; linarith
  set Lr : ℝ := ‖B‖ * (‖x₀‖ + 1) ^ 2 + ‖A‖ * (‖x₀‖ + 1) with hLr
  have hLr0 : 0 ≤ Lr := by
    rw [hLr]
    exact add_nonneg
      (mul_nonneg (ContinuousLinearMap.opNorm_nonneg B) (sq_nonneg _))
      (mul_nonneg (ContinuousLinearMap.opNorm_nonneg A) hx1)
  set Kr : ℝ := 2 * ‖B‖ * (‖x₀‖ + 1) + ‖A‖ with hKr
  have hKr0 : 0 ≤ Kr := by
    rw [hKr]
    exact add_nonneg
      (mul_nonneg (mul_nonneg (by norm_num) (ContinuousLinearMap.opNorm_nonneg B)) hx1)
      (ContinuousLinearMap.opNorm_nonneg A)
  have hLr1 : (0 : ℝ) < Lr + 1 := by linarith
  have hε : 0 < 1 / (Lr + 1) := one_div_pos.mpr hLr1
  have hg : ContinuousOn (fun _ : ℝ => (0 : E))
      (Icc (t₀ - 1 / (Lr + 1)) (t₀ + 1 / (Lr + 1))) := continuousOn_const
  have hgb : ∀ t ∈ Icc (t₀ - 1 / (Lr + 1)) (t₀ + 1 / (Lr + 1)),
      ‖(fun _ : ℝ => (0 : E)) t‖ ≤ ((0 : ℝ≥0) : ℝ) := by
    intro t _; simp
  have hK : 2 * ‖B‖ * (‖x₀‖ + ((1 : ℝ≥0) : ℝ)) + ‖A‖ ≤ ((Real.toNNReal Kr : ℝ≥0) : ℝ) := by
    rw [Real.coe_toNNReal _ hKr0, NNReal.coe_one, hKr]
  have hL : ((0 : ℝ≥0) : ℝ) + ‖B‖ * (‖x₀‖ + ((1 : ℝ≥0) : ℝ)) ^ 2
      + ‖A‖ * (‖x₀‖ + ((1 : ℝ≥0) : ℝ)) ≤ ((Real.toNNReal Lr : ℝ≥0) : ℝ) := by
    rw [Real.coe_toNNReal _ hLr0, NNReal.coe_zero, NNReal.coe_one, hLr]
    linarith
  have hstep : ((Real.toNNReal Lr : ℝ≥0) : ℝ) * (1 / (Lr + 1)) ≤ ((1 : ℝ≥0) : ℝ) := by
    rw [Real.coe_toNNReal _ hLr0, NNReal.coe_one, mul_one_div]
    rw [div_le_one hLr1]
    linarith
  obtain ⟨α, hα₀, hα⟩ :=
    galerkin_local_solution B A (fun _ : ℝ => (0 : E)) x₀ hε hg hgb hK hL hstep
  refine ⟨1 / (Lr + 1), hε, α, hα₀, fun t ht => ?_⟩
  have h' := hα t ht
  have heq : (fun _ : ℝ => (0 : E)) t + B (α t) (α t) + A (α t)
      = quadraticField A B (α t) := by
    show (0 : E) + B (α t) (α t) + A (α t) = quadraticField A B (α t)
    simp only [quadraticField_apply]
    abel
  rw [heq] at h'
  exact h'

/-- **GP2', interior form.**  The same local solution with a genuine two-sided
`HasDerivAt` at every time of the open interval `(t₀ - ε, t₀ + ε)`, which
contains the datum time.  This is the "local solution on a nontrivial
interval" statement in its strongest pointwise sense. -/
theorem quadratic_ode_local_solution_hasDerivAt [CompleteSpace E]
    (A : E →L[ℝ] E) (B : E →L[ℝ] E →L[ℝ] E) (t₀ : ℝ) (x₀ : E) :
    ∃ ε : ℝ, 0 < ε ∧ ∃ f : ℝ → E, f t₀ = x₀ ∧
      ∀ t ∈ Ioo (t₀ - ε) (t₀ + ε), HasDerivAt f (quadraticField A B (f t)) t := by
  obtain ⟨ε, hε, f, hf₀, hf⟩ := quadratic_ode_local_solution A B t₀ x₀
  refine ⟨ε, hε, f, hf₀, fun t ht => ?_⟩
  exact (hf t (Ioo_subset_Icc_self ht)).hasDerivAt (Icc_mem_nhds ht.1 ht.2)

/-- **GP3, uniqueness.**  Two solutions of `u' = A u + B u u` with the same
datum that both remain in a common closed ball agree on the whole interval.
The Lipschitz constant of GP1 on that ball feeds mathlib's Grönwall-based
`ODE_solution_unique_of_mem_Icc_right`.  Confinement to *some* ball is not a
restriction for the certificate: the Picard box of the Python layer is
exactly such a ball, and the enclosure inside it is the thing being
certified.  This is the uniqueness half of the finite-dimensional Galerkin
statement; the uniqueness half of EXT-P1 — for the PDE — is not claimed. -/
theorem quadratic_ode_unique
    (A : E →L[ℝ] E) (B : E →L[ℝ] E →L[ℝ] E) {x₀ : E} {r : ℝ} (hr : 0 ≤ r)
    {t₀ t₁ : ℝ} {f g : ℝ → E}
    (hf : ContinuousOn f (Icc t₀ t₁)) (hg : ContinuousOn g (Icc t₀ t₁))
    (hf' : ∀ t ∈ Ico t₀ t₁, HasDerivWithinAt f (quadraticField A B (f t)) (Ici t) t)
    (hg' : ∀ t ∈ Ico t₀ t₁, HasDerivWithinAt g (quadraticField A B (g t)) (Ici t) t)
    (hfm : ∀ t ∈ Ico t₀ t₁, f t ∈ closedBall x₀ r)
    (hgm : ∀ t ∈ Ico t₀ t₁, g t ∈ closedBall x₀ r)
    (h₀ : f t₀ = g t₀) :
    EqOn f g (Icc t₀ t₁) :=
  ODE_solution_unique_of_mem_Icc_right
    (v := fun _ x => quadraticField A B x) (s := fun _ => closedBall x₀ r)
    (fun _ _ => quadratic_field_lipschitzOnWith A B x₀ hr)
    hf hf' hfm hg hg' hgm h₀

end

end NSSingularity
