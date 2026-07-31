/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under the Apache 2.0 licence.
-/
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.ODE.Gronwall
import Mathlib.Analysis.ODE.ExistUnique
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import NSSingularity.CertificateLayer

/-!
# The control-ODE layer of an a posteriori error bound

An a posteriori error bound for a computed trajectory has three layers.  The
*residual* layer measures, by interval arithmetic, how badly the computed field
fails to solve the equation.  The *estimate* layer turns that residual into a
differential inequality for a scalar quantity `R t` — the size of the error in
whatever norm the estimate closes in.  The *ODE* layer integrates that
differential inequality.  This file is the third layer, and nothing else.

Every theorem here is a statement about real-valued functions satisfying a
differential inequality.  None of them says anything about Navier–Stokes: the
inequality is a hypothesis, and supplying it is the business of
`docs/research_notes/hs_error_propagation.md`, whose obligations `HS-1`, `HS-2a`,
`HS-2b`, `HS-4` and `HS-5` remain open.  What *is* closed here is `HS-6`.

## Why there is no linear damping

The natural comparison object for a spectral computation is Morosi and
Pizzocchero's control ODE for the Sobolev norm `R_n` of the error,

```
d⁺R_n/dt  ≥  -ν R_n + (G_n D_n + K_n D_{n+1}) R_n + G_n R_n² + ε_n .
```

The `-ν R_n` term is the whole reason their smallness and global-existence
corollaries work: it is a genuine linear damping which, for small enough data,
beats the quadratic term forever.  It comes from the Poincaré inequality on the
torus, where the admissible wavenumbers are `k ∈ ℤ^d \ {0}` and hence `|k| ≥ 1`,
so the Laplacian has a spectral gap and `‖∇v‖ ≥ ‖v‖` on mean-zero fields.

On `ℝ³` the spectrum of the Laplacian is `[0, ∞)` with no gap.  There is no
Poincaré inequality on the whole space, `‖∇v‖ ≥ c‖v‖` is false for every `c > 0`,
and the `-ν R_n` term simply is not there.  So the control ODE available on `ℝ³`
is the *undamped* one,

```
R' ≤ a R + b R² + e ,
```

and its solution blows up in finite time whatever the data.  Every bound in this
file is therefore local in time by construction, and the existence interval it
certifies is strictly shorter than the periodic analogue would give.  That is not
a defect of the proof; it is the geometry of the problem.  It is also why nothing
here can be chained up to a putative blow-up time: the majorant blows up first.

## What is proved

* `nonpos_of_deriv_le_mul_of_pos` — the Chaplygin–Dini one-sided comparison
  lemma.  If `Z t₀ ≤ 0` and `Z' t ≤ M Z t` *at those times where `Z t > 0`*, then
  `Z ≤ 0` throughout.  This is the reusable half of Morosi–Pizzocchero's
  machinery, and everything else in the file is a corollary of it or of the same
  fencing theorem.  The hypothesis is deliberately conditional on `Z t > 0`,
  because that is all that a comparison argument ever gets to use, and demanding
  the inequality everywhere would exclude the applications.
* `gronwall_variable_coefficient` — **HS-6**, the variable-coefficient Grönwall
  step, in the form that needs no measure theory: the caller supplies
  antiderivatives.  No sign condition on `a` or `e`.
* `gronwall_variable_coefficient_integral` — **HS-6** in the literal form (3.2)
  of `hs_error_propagation.md`, with the antiderivatives written as interval
  integrals, obtained from the previous theorem by the fundamental theorem of
  calculus.
* `riccati_comparison` — the full quadratic comparison: a subsolution of
  `R' = a R + b R² + e` never overtakes a supersolution that starts above it.  No
  a priori bound on `R` or `S` is assumed; the Lipschitz constant the argument
  needs is manufactured from continuity on a compact interval.
* `le_quadratic_bound` — the explicit blow-up majorant `R₀/(1 - b R₀ t)` for
  `R' ≤ b R²`, valid on the whole interval where the majorant is finite.  This is
  the quantitative statement of "the undamped control ODE gives a short interval
  and no more".
* `clampTo`, `roughEnclosure_isPicardLindelof`, `exists_roughEnclosure_solution`,
  `roughEnclosure_solution_unique` — the rough-enclosure (Picard self-mapping)
  lemma in scalar form.

## The rough enclosure, and what it replaces

`slab_certificate.py` currently bounds a solution over a time slab through a
hypothesis called `H2`: a two-point Hermite remainder `|y - H| ≤ Δ⁴ M₄/384`.  That
bound is a theorem only when the interpolant `H` matches *one* trajectory at both
ends.  In the certificate the right-hand endpoint is an RK4 output and the
right-hand slope belongs to a different trajectory, so the hypothesis of the
Hermite remainder is not met and `H2` is not repairable as framed.

`exists_roughEnclosure_solution` is the replacement, and it is a theorem with no
hypothesis beyond finitely many checkable inequalities.  Given a box `B = [lo,hi]`
containing the initial value `y₀`, a Lipschitz bound for the field on `B`, and an
interval enclosure `[Flo, Fhi] ⊇ F(B)` satisfying the *self-mapping* conditions

```
y₀ + h Fhi ≤ hi ,      lo ≤ y₀ + h Flo ,
```

there exists a solution on the whole slab `[0,h]` and it stays in `B`.  Those two
inequalities are exactly the statement `y₀ + [0,h]·F(B) ⊆ B`, and they are exactly
what an interval-arithmetic certificate already computes.

The proof does not use a fixed-point argument by hand.  It clamps the field to
`B`, which makes it globally Lipschitz and globally bounded without changing it on
`B`, hands the clamped field to mathlib's `IsPicardLindelof`, and then shows that
the resulting solution never leaves `B` — so that the clamp was inert and the
solution solves the original equation.  The enclosure is proved by fencing the
solution between the two affine functions `y₀ + Flo t` and `y₀ + Fhi t`, which is
where the self-mapping inequalities are consumed.

## Scope

`CertificateLayer.lean` proves the constant-coefficient Grönwall step
(`norm_le_simple_gronwall`, `gronwallBound_le_simple`); this file extends it and
does not restate it.  There is no `sorry`, no `admit` and no project-specific
axiom.  Nothing here licenses any claim about the Clay problem.
-/

namespace NSSingularity

noncomputable section

open Set Metric MeasureTheory
open scoped NNReal Topology

/-! ## The Chaplygin–Dini one-sided comparison lemma -/

section Dini

/-- **The one-sided comparison lemma.**  A function that starts nonpositive and
whose right derivative is dominated by `M Z` *whenever it is positive* never
becomes positive.

The conditional form of the growth hypothesis is the point.  In every application
the inequality `Z' ≤ M Z` is only available where `Z > 0`: it is obtained by
bounding a coefficient `c t` that multiplies `Z`, and `c t Z ≤ M Z` requires
`Z ≥ 0`.  Demanding the inequality at all times would make the lemma unusable.

The proof is the classical strict-boundary trick.  For each `δ > 0` the curve
`B t = δ exp((M+1)(t - t₀))` is a *strict* supersolution: where `Z` touches `B`
one has `Z = B > 0`, hence `Z' ≤ M Z < (M+1) Z = B'`, and mathlib's fencing
theorem `image_le_of_deriv_right_lt_deriv_boundary` applies.  Letting `δ → 0`
gives `Z ≤ 0`.  The extra `+1` in the exponent is what buys the strictness that
the fencing theorem demands and that the ODE itself does not supply. -/
theorem nonpos_of_deriv_le_mul_of_pos {Z Z' : ℝ → ℝ} {M t₀ T : ℝ}
    (hZ : ContinuousOn Z (Icc t₀ T))
    (hZ' : ∀ t ∈ Ico t₀ T, HasDerivWithinAt Z (Z' t) (Ici t) t)
    (h₀ : Z t₀ ≤ 0)
    (hgrow : ∀ t ∈ Ico t₀ T, 0 < Z t → Z' t ≤ M * Z t) :
    ∀ t ∈ Icc t₀ T, Z t ≤ 0 := by
  intro t ht
  have key : ∀ δ : ℝ, 0 < δ → Z t ≤ δ * Real.exp ((M + 1) * (t - t₀)) := by
    intro δ hδ
    have hBderiv : ∀ x : ℝ,
        HasDerivAt (fun y : ℝ => δ * Real.exp ((M + 1) * (y - t₀)))
          ((M + 1) * (δ * Real.exp ((M + 1) * (x - t₀)))) x := by
      intro x
      have h1 : HasDerivAt (fun y : ℝ => (M + 1) * (y - t₀)) (M + 1) x := by
        simpa using ((hasDerivAt_id x).sub_const t₀).const_mul (M + 1)
      have h2 := h1.exp.const_mul δ
      have heq : δ * (Real.exp ((M + 1) * (x - t₀)) * (M + 1))
          = (M + 1) * (δ * Real.exp ((M + 1) * (x - t₀))) := by ring
      rw [← heq]
      exact h2
    refine image_le_of_deriv_right_lt_deriv_boundary
      (B := fun y : ℝ => δ * Real.exp ((M + 1) * (y - t₀)))
      (B' := fun y : ℝ => (M + 1) * (δ * Real.exp ((M + 1) * (y - t₀))))
      hZ hZ' ?_ hBderiv ?_ ht
    · simpa using h₀.trans hδ.le
    · intro x hx hcontact
      have hBpos : 0 < δ * Real.exp ((M + 1) * (x - t₀)) := by positivity
      have hZpos : 0 < Z x := hcontact ▸ hBpos
      have hexp : (M + 1) * Z x = M * Z x + Z x := by ring
      rw [← hcontact]
      linarith [hgrow x hx hZpos]
  rcases le_or_gt (Z t) 0 with hle | hpos
  · exact hle
  exfalso
  have hCpos : 0 < Real.exp ((M + 1) * (t - t₀)) := Real.exp_pos _
  have hCne : Real.exp ((M + 1) * (t - t₀)) ≠ 0 := ne_of_gt hCpos
  have h1 := key (Z t / (2 * Real.exp ((M + 1) * (t - t₀)))) (by positivity)
  have h2 : Z t / (2 * Real.exp ((M + 1) * (t - t₀))) * Real.exp ((M + 1) * (t - t₀))
      = Z t / 2 := by
    field_simp
  rw [h2] at h1
  linarith

end Dini

/-! ## HS-6 — Grönwall with a time-dependent coefficient -/

section VariableGronwall

/-- **HS-6, the antiderivative form.**  If `R' t ≤ a t · R t + e t` on `[t₀, T]`
and `R t₀ ≤ R₀`, then

`R t ≤ exp(A t) · (R₀ + Φ t)`,

where `A` is any antiderivative of `a` vanishing at `t₀` and `Φ` any
antiderivative of `s ↦ exp(-A s) e s` vanishing at `t₀`.

This is the statement `hs_error_propagation.md` §3.5 needs and
`CertificateLayer.norm_le_simple_gronwall` does not supply: mathlib's Grönwall
inequality, and hence F-19, has a *constant* rate `K`, whereas the `H^s` estimate
(3.1) produces the rate `C_s M(t)` with `M(t) = ‖u‖_{H^{s+1}} + ‖ũ‖_{H^{s+1}}`,
which is a genuine function of time.

Antiderivatives are taken as hypotheses rather than built as integrals so that no
integrability or measurability side condition appears, and so that a caller with
a closed-form `A` — the usual case when `a` is a polynomial in `t` produced by a
certificate — pays nothing for the fundamental theorem of calculus.  The integral
form is `gronwall_variable_coefficient_integral`.

No sign hypothesis is imposed on `a` or on `e`.  The proof is the integrating
factor: `Y = exp(-A) R - Φ` has `Y' = exp(-A)(R' - a R - e) ≤ 0`, so `Y` is
nonincreasing, and `Y t₀ = R t₀`. -/
theorem gronwall_variable_coefficient {R R' a e A Φ : ℝ → ℝ} {R₀ t₀ T : ℝ}
    (hA : ∀ t, HasDerivAt A (a t) t) (hA₀ : A t₀ = 0)
    (hΦ : ∀ t, HasDerivAt Φ (Real.exp (-A t) * e t) t) (hΦ₀ : Φ t₀ = 0)
    (hR : ContinuousOn R (Icc t₀ T))
    (hR' : ∀ t ∈ Ico t₀ T, HasDerivWithinAt R (R' t) (Ici t) t)
    (h₀ : R t₀ ≤ R₀)
    (hbound : ∀ t ∈ Ico t₀ T, R' t ≤ a t * R t + e t) :
    ∀ t ∈ Icc t₀ T, R t ≤ Real.exp (A t) * (R₀ + Φ t) := by
  have hAdiff : Differentiable ℝ A := fun t => (hA t).differentiableAt
  have hΦdiff : Differentiable ℝ Φ := fun t => (hΦ t).differentiableAt
  have hAcont : Continuous A := hAdiff.continuous
  have hΦcont : Continuous Φ := hΦdiff.continuous
  have hYcont : ContinuousOn (fun s => Real.exp (-A s) * R s - Φ s) (Icc t₀ T) :=
    (((Real.continuous_exp.comp hAcont.neg).continuousOn).mul hR).sub hΦcont.continuousOn
  have hYderiv : ∀ t ∈ Ico t₀ T,
      HasDerivWithinAt (fun s => Real.exp (-A s) * R s - Φ s)
        (Real.exp (-A t) * (R' t - a t * R t - e t)) (Ici t) t := by
    intro t ht
    have h1 : HasDerivWithinAt (fun s => Real.exp (-A s)) (Real.exp (-A t) * -(a t)) (Ici t) t :=
      ((hA t).neg.exp).hasDerivWithinAt
    have heq : Real.exp (-A t) * (R' t - a t * R t - e t)
        = Real.exp (-A t) * -(a t) * R t + Real.exp (-A t) * R' t - Real.exp (-A t) * e t := by
      ring
    rw [heq]
    exact (h1.mul (hR' t ht)).sub ((hΦ t).hasDerivWithinAt)
  have hYle : ∀ t ∈ Icc t₀ T,
      Real.exp (-A t) * R t - Φ t ≤ Real.exp (-A t₀) * R t₀ - Φ t₀ := by
    intro t ht
    refine image_le_of_deriv_right_le_deriv_boundary
      (B := fun _ : ℝ => Real.exp (-A t₀) * R t₀ - Φ t₀) (B' := fun _ : ℝ => (0 : ℝ))
      hYcont hYderiv le_rfl continuousOn_const (fun x _ => hasDerivWithinAt_const x _ _) ?_ ht
    intro x hx
    have hexp : (0 : ℝ) < Real.exp (-A x) := Real.exp_pos _
    have hneg : R' x - a x * R x - e x ≤ 0 := by linarith [hbound x hx]
    have := mul_le_mul_of_nonneg_left hneg hexp.le
    simpa using this
  intro t ht
  have hstart : Real.exp (-A t₀) * R t₀ - Φ t₀ ≤ R₀ := by
    rw [hA₀, hΦ₀, neg_zero, Real.exp_zero, one_mul, sub_zero]
    exact h₀
  have hmain : Real.exp (-A t) * R t - Φ t ≤ R₀ := (hYle t ht).trans hstart
  have hexp : (0 : ℝ) < Real.exp (A t) := Real.exp_pos _
  have hid : Real.exp (A t) * ((Real.exp (-A t) * R t - Φ t) + Φ t) = R t := by
    have hne : Real.exp (A t) ≠ 0 := ne_of_gt hexp
    rw [Real.exp_neg]
    field_simp
    ring
  calc R t = Real.exp (A t) * ((Real.exp (-A t) * R t - Φ t) + Φ t) := hid.symm
    _ ≤ Real.exp (A t) * (R₀ + Φ t) :=
        mul_le_mul_of_nonneg_left (by linarith) hexp.le

/-- **HS-6, the integral form.**  Literally inequality (3.2) of
`hs_error_propagation.md`:

`R t ≤ exp(∫_{t₀}^{t} a) · (R₀ + ∫_{t₀}^{t} exp(-∫_{t₀}^{s} a) e s ds)`.

Continuity of `a` and `e` is what the fundamental theorem of calculus wants; it
is not needed by the underlying comparison argument, which is why
`gronwall_variable_coefficient` is stated separately and is the theorem to reach
for when the antiderivative is known in closed form. -/
theorem gronwall_variable_coefficient_integral {R R' a e : ℝ → ℝ} {R₀ t₀ T : ℝ}
    (ha : Continuous a) (he : Continuous e)
    (hR : ContinuousOn R (Icc t₀ T))
    (hR' : ∀ t ∈ Ico t₀ T, HasDerivWithinAt R (R' t) (Ici t) t)
    (h₀ : R t₀ ≤ R₀)
    (hbound : ∀ t ∈ Ico t₀ T, R' t ≤ a t * R t + e t) :
    ∀ t ∈ Icc t₀ T,
      R t ≤ Real.exp (∫ s in t₀..t, a s) *
        (R₀ + ∫ s in t₀..t, Real.exp (-∫ τ in t₀..s, a τ) * e s) := by
  have hAderiv : ∀ t : ℝ, HasDerivAt (fun u : ℝ => ∫ s in t₀..u, a s) (a t) t := fun t =>
    intervalIntegral.integral_hasDerivAt_right (ha.intervalIntegrable _ _)
      (ha.stronglyMeasurableAtFilter _ _) ha.continuousAt
  have hAdiff : Differentiable ℝ fun u : ℝ => ∫ s in t₀..u, a s :=
    fun t => (hAderiv t).differentiableAt
  have hAcont : Continuous (fun u : ℝ => ∫ s in t₀..u, a s) := hAdiff.continuous
  have hgcont : Continuous (fun s : ℝ => Real.exp (-∫ τ in t₀..s, a τ) * e s) :=
    (Real.continuous_exp.comp hAcont.neg).mul he
  have hΦderiv : ∀ t : ℝ,
      HasDerivAt (fun u : ℝ => ∫ s in t₀..u, Real.exp (-∫ τ in t₀..s, a τ) * e s)
        (Real.exp (-∫ τ in t₀..t, a τ) * e t) t := fun t =>
    intervalIntegral.integral_hasDerivAt_right (hgcont.intervalIntegrable _ _)
      (hgcont.stronglyMeasurableAtFilter _ _) hgcont.continuousAt
  exact gronwall_variable_coefficient (A := fun u : ℝ => ∫ s in t₀..u, a s)
    (Φ := fun u : ℝ => ∫ s in t₀..u, Real.exp (-∫ τ in t₀..s, a τ) * e s)
    hAderiv (by simp) hΦderiv (by simp) hR hR' h₀ hbound

end VariableGronwall

/-! ## The quadratic control ODE -/

section Riccati

/-- **Riccati comparison.**  A subsolution of `R' = a R + b R² + e` that starts
below a supersolution stays below it, on the whole interval where both are
defined.

This is the comparison the undamped control ODE of the introduction needs: `a`
carries the linear part of the estimate, `b` the quadratic advection term, and
`e` the residual.  Only `b ≥ 0` is assumed; `a` and `e` are unrestricted in sign,
and neither function is assumed nonnegative.

No a priori bound on `R` or on `S` appears in the hypotheses.  One is needed —
the difference `Z = R - S` satisfies `Z' ≤ (a + b(R+S)) Z` and the coefficient
has to be bounded before `nonpos_of_deriv_le_mul_of_pos` applies — but it is
manufactured inside the proof from continuity on the compact interval `[t₀, T]`.
That is the honest form of "the Riccati comparison is local": the constant exists
but is not effective, so this theorem certifies the *ordering* of the two
trajectories and never a quantitative rate. -/
theorem riccati_comparison {R R' S S' : ℝ → ℝ} {a b e t₀ T : ℝ}
    (hb : 0 ≤ b)
    (hR : ContinuousOn R (Icc t₀ T)) (hS : ContinuousOn S (Icc t₀ T))
    (hR' : ∀ t ∈ Ico t₀ T, HasDerivWithinAt R (R' t) (Ici t) t)
    (hS' : ∀ t ∈ Ico t₀ T, HasDerivWithinAt S (S' t) (Ici t) t)
    (hsub : ∀ t ∈ Ico t₀ T, R' t ≤ a * R t + b * R t ^ 2 + e)
    (hsuper : ∀ t ∈ Ico t₀ T, a * S t + b * S t ^ 2 + e ≤ S' t)
    (h₀ : R t₀ ≤ S t₀) :
    ∀ t ∈ Icc t₀ T, R t ≤ S t := by
  obtain ⟨CR, hCR⟩ := isCompact_Icc.exists_bound_of_continuousOn hR
  obtain ⟨CS, hCS⟩ := isCompact_Icc.exists_bound_of_continuousOn hS
  have hZcont : ContinuousOn (fun t => R t - S t) (Icc t₀ T) := hR.sub hS
  have hZderiv : ∀ t ∈ Ico t₀ T,
      HasDerivWithinAt (fun t => R t - S t) (R' t - S' t) (Ici t) t :=
    fun t ht => (hR' t ht).sub (hS' t ht)
  have hgrow : ∀ t ∈ Ico t₀ T, 0 < R t - S t →
      R' t - S' t ≤ (a + b * (CR + CS)) * (R t - S t) := by
    intro t ht hpos
    have htI : t ∈ Icc t₀ T := Ico_subset_Icc_self ht
    have hRb : R t ≤ CR := (le_abs_self _).trans (by simpa [Real.norm_eq_abs] using hCR t htI)
    have hSb : S t ≤ CS := (le_abs_self _).trans (by simpa [Real.norm_eq_abs] using hCS t htI)
    have hid : (a + b * (R t + S t)) * (R t - S t)
        = a * R t + b * R t ^ 2 - (a * S t + b * S t ^ 2) := by ring
    have h1 : R' t - S' t ≤ (a + b * (R t + S t)) * (R t - S t) := by
      linarith [hsub t ht, hsuper t ht]
    have hcoef : a + b * (R t + S t) ≤ a + b * (CR + CS) := by
      nlinarith [mul_nonneg hb (by linarith : (0 : ℝ) ≤ CR + CS - (R t + S t))]
    exact h1.trans (mul_le_mul_of_nonneg_right hcoef hpos.le)
  have hnonpos :=
    nonpos_of_deriv_le_mul_of_pos hZcont hZderiv (sub_nonpos.2 h₀) hgrow
  intro t ht
  linarith [hnonpos t ht]

/-- **The quadratic blow-up majorant.**  For `R' ≤ b R²` with `R t₀ ≤ R₀` the
explicit solution `R₀/(1 - b R₀ t)` is an upper bound on the whole interval where
it is finite.

This is the undamped control ODE in its purest form, and the statement makes
explicit what the missing `-ν R` term of Morosi–Pizzocchero would have bought:
without it the majorant is a hyperbola with a pole at `t = 1/(b R₀)`, and no
choice of data pushes the pole to infinity.  The bound is proved for all
`t ∈ [0, T]` with `b R₀ T < 1`, which is the largest interval on which the
right-hand side means anything.

No sign condition on `R` is needed — the comparison does not use one — and the
proof is `riccati_comparison` with `a = e = 0` against the exact solution. -/
theorem le_quadratic_bound {R R' : ℝ → ℝ} {b R₀ T : ℝ}
    (hb : 0 ≤ b) (hR₀ : 0 ≤ R₀) (hTb : b * R₀ * T < 1)
    (hR : ContinuousOn R (Icc 0 T))
    (hR' : ∀ t ∈ Ico (0 : ℝ) T, HasDerivWithinAt R (R' t) (Ici t) t)
    (h₀ : R 0 ≤ R₀)
    (hbound : ∀ t ∈ Ico (0 : ℝ) T, R' t ≤ b * R t ^ 2) :
    ∀ t ∈ Icc (0 : ℝ) T, R t ≤ R₀ / (1 - b * R₀ * t) := by
  have hden : ∀ t ∈ Icc (0 : ℝ) T, 0 < 1 - b * R₀ * t := by
    intro t ht
    have hbR : (0 : ℝ) ≤ b * R₀ := mul_nonneg hb hR₀
    have : b * R₀ * t ≤ b * R₀ * T := mul_le_mul_of_nonneg_left ht.2 hbR
    linarith
  have hScont : ContinuousOn (fun t => R₀ / (1 - b * R₀ * t)) (Icc 0 T) := by
    refine ContinuousOn.div continuousOn_const (by fun_prop) ?_
    intro t ht
    exact ne_of_gt (hden t ht)
  have hSderiv : ∀ t ∈ Ico (0 : ℝ) T,
      HasDerivWithinAt (fun t => R₀ / (1 - b * R₀ * t))
        (b * (R₀ / (1 - b * R₀ * t)) ^ 2) (Ici t) t := by
    intro t ht
    have hne : (1 - b * R₀ * t) ≠ 0 := ne_of_gt (hden t (Ico_subset_Icc_self ht))
    have hd : HasDerivAt (fun s : ℝ => 1 - b * R₀ * s) (-(b * R₀)) t := by
      simpa using ((hasDerivAt_id t).const_mul (b * R₀)).const_sub 1
    have hq := (hasDerivAt_const t R₀).div hd hne
    have heq : (0 * (1 - b * R₀ * t) - R₀ * -(b * R₀)) / (1 - b * R₀ * t) ^ 2
        = b * (R₀ / (1 - b * R₀ * t)) ^ 2 := by
      field_simp
      ring
    rw [← heq]
    exact hq.hasDerivWithinAt
  have hstart : R 0 ≤ R₀ / (1 - b * R₀ * 0) := by simpa using h₀
  exact riccati_comparison (a := 0) (e := 0) hb hR hScont hR' hSderiv
    (fun t ht => by simpa using hbound t ht) (fun t _ => by simp) hstart

end Riccati

/-! ## The rough enclosure: a hypothesis-free replacement for `H2` -/

section RoughEnclosure

/-- Clamp a real number into `[lo, hi]`.  Composing a vector field with this map
leaves it unchanged on the box and makes it globally Lipschitz and globally
bounded off it, which is the whole trick behind `exists_roughEnclosure_solution`.
-/
def clampTo (lo hi x : ℝ) : ℝ := max lo (min hi x)

/-- The clamp lands in the box. -/
theorem clampTo_mem {lo hi : ℝ} (h : lo ≤ hi) (x : ℝ) : clampTo lo hi x ∈ Icc lo hi :=
  ⟨le_max_left _ _, max_le h (min_le_left _ _)⟩

/-- The clamp is the identity on the box, so composing with it changes nothing
where it matters. -/
theorem clampTo_eq_self {lo hi x : ℝ} (hx : x ∈ Icc lo hi) : clampTo lo hi x = x := by
  simp only [clampTo, min_eq_right hx.2, max_eq_right hx.1]

/-- The clamp is a contraction, so it does not spoil a Lipschitz constant. -/
theorem abs_clampTo_sub_clampTo_le (lo hi x y : ℝ) :
    |clampTo lo hi x - clampTo lo hi y| ≤ |x - y| := by
  have h1 : |min hi x - min hi y| ≤ |x - y| :=
    calc |min hi x - min hi y| ≤ max |hi - hi| |x - y| := abs_min_sub_min_le_max hi x hi y
      _ = |x - y| := by rw [sub_self, abs_zero]; exact max_eq_right (abs_nonneg _)
  have h2 : |max lo (min hi x) - max lo (min hi y)| ≤ |min hi x - min hi y| :=
    calc |max lo (min hi x) - max lo (min hi y)| ≤ max |lo - lo| |min hi x - min hi y| :=
          abs_max_sub_max_le_max lo (min hi x) lo (min hi y)
      _ = |min hi x - min hi y| := by rw [sub_self, abs_zero]; exact max_eq_right (abs_nonneg _)
  exact h2.trans h1

/-- **The Picard–Lindelöf conditions for the clamped field.**  The field
`x ↦ F (clamp x)` satisfies all four of mathlib's conditions on *every* ball, with
the Lipschitz constant and the sup bound that `F` has on the box, because the
clamp never lets its argument leave the box.

The fourth condition, `mul_max_le`, is precisely the self-mapping condition
`L h ≤ a` of the introduction; here it is met by choosing the ball radius, and the
sharper two-sided version is what `exists_roughEnclosure_solution` then uses to
prove the enclosure. -/
theorem roughEnclosure_isPicardLindelof (F : ℝ → ℝ) {lo hi y₀ h : ℝ} {a L K : ℝ≥0}
    (hlohi : lo ≤ hi) (hh : 0 < h)
    (hlip : LipschitzOnWith K F (Icc lo hi))
    (hLb : ∀ x ∈ Icc lo hi, |F x| ≤ (L : ℝ))
    (hstep : (L : ℝ) * h ≤ (a : ℝ)) :
    IsPicardLindelof (fun (_ : ℝ) (x : ℝ) => F (clampTo lo hi x))
      (⟨0, ⟨le_rfl, hh.le⟩⟩ : Icc (0 : ℝ) h) y₀ a 0 L K where
  lipschitzOnWith := by
    intro t _
    refine LipschitzOnWith.of_dist_le_mul fun x _ y _ => ?_
    have hcx := clampTo_mem hlohi x
    have hcy := clampTo_mem hlohi y
    calc dist (F (clampTo lo hi x)) (F (clampTo lo hi y))
        ≤ (K : ℝ) * dist (clampTo lo hi x) (clampTo lo hi y) := hlip.dist_le_mul _ hcx _ hcy
      _ ≤ (K : ℝ) * dist x y := by
          rw [Real.dist_eq, Real.dist_eq]
          exact mul_le_mul_of_nonneg_left (abs_clampTo_sub_clampTo_le lo hi x y) K.coe_nonneg
  continuousOn := fun _ _ => continuousOn_const
  norm_le := by
    intro t _ x _
    rw [Real.norm_eq_abs]
    exact hLb _ (clampTo_mem hlohi _)
  mul_max_le := by
    simpa [max_eq_left hh.le] using hstep

/-- **The rough enclosure, and the replacement for the slab certificate's `H2`.**

Hypotheses, all of them finitely checkable by interval arithmetic:

* `y₀ ∈ [lo, hi]` — the initial value lies in the box;
* `F` is Lipschitz with constant `K` on the box, and `|F| ≤ L` there;
* `[Flo, Fhi]` encloses `F` on the box;
* the *self-mapping* inequalities `y₀ + h Fhi ≤ hi` and `lo ≤ y₀ + h Flo`.

Conclusion: a solution exists on the entire slab `[0, h]`, starts at `y₀`, and
never leaves the box.

The last two hypotheses say exactly `y₀ + [0,h]·F(B) ⊆ B`.  Nothing is assumed
about a numerical integrator, an interpolant, or a fourth derivative; in
particular the two-point Hermite remainder that the current `slab_certificate.py`
calls `H2` — which is a theorem only when one trajectory is interpolated at both
ends, and therefore not applicable to an RK4 endpoint paired with a foreign slope
— is not used and not needed.

The proof clamps `F` to the box.  The clamped field agrees with `F` on the box,
is globally Lipschitz and globally bounded, so mathlib's Picard–Lindelöf theorem
applies on the ball of radius `L h` about `y₀` with no smallness condition at all.
The resulting solution is then fenced between the affine functions `y₀ + Flo t`
and `y₀ + Fhi t` — this is where the two self-mapping inequalities are spent, and
convexity of the box in `t` is what lets a one-sided enclosure of `F` do the work
of a symmetric one.  Being inside the box, the solution never sees the clamp, so
it solves the original equation. -/
theorem exists_roughEnclosure_solution (F : ℝ → ℝ) {lo hi y₀ h Flo Fhi : ℝ} {L K : ℝ≥0}
    (hh : 0 < h) (hy₀ : y₀ ∈ Icc lo hi)
    (hlip : LipschitzOnWith K F (Icc lo hi))
    (hLb : ∀ x ∈ Icc lo hi, |F x| ≤ (L : ℝ))
    (hFlo : ∀ x ∈ Icc lo hi, Flo ≤ F x)
    (hFhi : ∀ x ∈ Icc lo hi, F x ≤ Fhi)
    (hself_hi : y₀ + h * Fhi ≤ hi)
    (hself_lo : lo ≤ y₀ + h * Flo) :
    ∃ y : ℝ → ℝ, y 0 = y₀ ∧ (∀ t ∈ Icc (0 : ℝ) h, y t ∈ Icc lo hi) ∧
      ∀ t ∈ Icc (0 : ℝ) h, HasDerivWithinAt y (F (y t)) (Icc (0 : ℝ) h) t := by
  have hlohi : lo ≤ hi := hy₀.1.trans hy₀.2
  have hstep : (L : ℝ) * h ≤ ((L * Real.toNNReal h : ℝ≥0) : ℝ) := by
    rw [NNReal.coe_mul, Real.coe_toNNReal _ hh.le]
  obtain ⟨α, hα₀, hα⟩ :=
    (roughEnclosure_isPicardLindelof F hlohi hh hlip hLb
      hstep).exists_eq_forall_mem_Icc_hasDerivWithinAt₀
  have hα₀' : α 0 = y₀ := hα₀
  have hαcont : ContinuousOn α (Icc (0 : ℝ) h) := fun t ht => (hα t ht).continuousWithinAt
  have hIci : ∀ t ∈ Ico (0 : ℝ) h,
      HasDerivWithinAt α (F (clampTo lo hi (α t))) (Ici t) t := by
    intro t ht
    refine (hα t (Ico_subset_Icc_self ht)).mono_of_mem_nhdsWithin ?_
    exact (Set.ordConnected_Icc).mem_nhdsGE (Ico_subset_Icc_self ht)
      (right_mem_Icc.2 hh.le) ht.2
  -- The upper fence.
  have hup : ∀ t ∈ Icc (0 : ℝ) h, α t ≤ y₀ + Fhi * t := by
    intro t ht
    refine image_le_of_deriv_right_le_deriv_boundary
      (B := fun s : ℝ => y₀ + Fhi * s) (B' := fun _ : ℝ => Fhi)
      hαcont hIci ?_ (by fun_prop) ?_ ?_ ht
    · rw [hα₀']; simp
    · intro x _
      have : HasDerivAt (fun s : ℝ => y₀ + Fhi * s) Fhi x := by
        simpa using ((hasDerivAt_id x).const_mul Fhi).const_add y₀
      exact this.hasDerivWithinAt
    · intro x _
      exact hFhi _ (clampTo_mem hlohi _)
  -- The lower fence, obtained by fencing `-α` from above.
  have hdown : ∀ t ∈ Icc (0 : ℝ) h, y₀ + Flo * t ≤ α t := by
    intro t ht
    have hneg : -α t ≤ -y₀ - Flo * t := by
      refine image_le_of_deriv_right_le_deriv_boundary
        (f := fun s => -α s) (f' := fun s => -F (clampTo lo hi (α s)))
        (B := fun s : ℝ => -y₀ - Flo * s) (B' := fun _ : ℝ => -Flo)
        hαcont.neg (fun x hx => (hIci x hx).neg) ?_ (by fun_prop) ?_ ?_ ht
      · rw [hα₀']; simp
      · intro x _
        have : HasDerivAt (fun s : ℝ => -y₀ - Flo * s) (-Flo) x := by
          simpa using ((hasDerivAt_id x).const_mul Flo).const_sub (-y₀)
        exact this.hasDerivWithinAt
      · intro x _
        exact neg_le_neg (hFlo _ (clampTo_mem hlohi _))
    linarith
  -- The enclosure itself: convexity of the box in `t` turns the one-sided
  -- enclosure of `F` into a two-sided enclosure of the solution.
  have hmem : ∀ t ∈ Icc (0 : ℝ) h, α t ∈ Icc lo hi := by
    intro t ht
    obtain ⟨ht0, hth⟩ := ht
    have hlow := hdown t ⟨ht0, hth⟩
    have hhigh := hup t ⟨ht0, hth⟩
    constructor
    · rcases le_or_gt 0 Flo with hs | hs
      · nlinarith [mul_nonneg hs ht0, hy₀.1]
      · nlinarith [mul_nonneg (le_of_lt (neg_pos.2 hs)) (sub_nonneg.2 hth)]
    · rcases le_or_gt 0 Fhi with hs | hs
      · nlinarith [mul_nonneg hs (sub_nonneg.2 hth)]
      · nlinarith [mul_nonneg (le_of_lt (neg_pos.2 hs)) ht0, hy₀.2]
  refine ⟨α, hα₀', hmem, ?_⟩
  intro t ht
  have hd := hα t ht
  rwa [clampTo_eq_self (hmem t ht)] at hd

/-- **Uniqueness on the slab.**  Two solutions that start together and both stay
in the box agree throughout it.  Together with `exists_roughEnclosure_solution`
this is the full statement the slab certificate wants: on `[0, h]` the initial
value problem has one and only one solution, and it is inside the box.

Note the hypotheses: uniqueness holds *among solutions that remain in the box*.
That is not a restriction in the intended use, since the existence theorem
produces such a solution and the enclosure is the thing being certified. -/
theorem roughEnclosure_solution_unique (F : ℝ → ℝ) {lo hi h : ℝ} {K : ℝ≥0}
    (hlip : LipschitzOnWith K F (Icc lo hi)) {y z : ℝ → ℝ}
    (hy : ContinuousOn y (Icc 0 h)) (hz : ContinuousOn z (Icc 0 h))
    (hy' : ∀ t ∈ Ico (0 : ℝ) h, HasDerivWithinAt y (F (y t)) (Ici t) t)
    (hz' : ∀ t ∈ Ico (0 : ℝ) h, HasDerivWithinAt z (F (z t)) (Ici t) t)
    (hym : ∀ t ∈ Ico (0 : ℝ) h, y t ∈ Icc lo hi)
    (hzm : ∀ t ∈ Ico (0 : ℝ) h, z t ∈ Icc lo hi)
    (h₀ : y 0 = z 0) :
    EqOn y z (Icc 0 h) :=
  ODE_solution_unique_of_mem_Icc_right (v := fun _ x => F x) (s := fun _ => Icc lo hi)
    (fun _ _ => hlip) hy hy' hym hz hz' hzm h₀

end RoughEnclosure

end

end NSSingularity
