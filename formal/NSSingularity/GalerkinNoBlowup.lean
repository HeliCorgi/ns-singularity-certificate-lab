/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
/-
Stage 1, item **F-6** of `docs/formalization_map.md`: the **Galerkin energy
bound**, i.e. the finite-dimensional core of the Track-F no-go theorem
(`docs/research_notes/track_f_finite_mode_nogo.md`, Theorem 1(i)).

Mathematical content.  A Track-F ansatz for Clay statement (C)/(D) prescribes a
velocity field `u` and defines the force as the residual
`f = ∂_t u + (u·∇)u - νΔu + ∇p`.  If `u(t)` is constrained to a *fixed*
finite-dimensional space `V` of smooth divergence-free fields — every
trigonometric-polynomial ansatz, every finite symbolic search space — then
projecting the momentum equation onto `V` gives a closed ordinary differential
equation

  `c' = g(t) + B(c,c) + A c`

in which
* `B` is the projected advection term and satisfies `⟪x, B x x⟫ = 0`, because
  `∫ u·(u·∇)u = 0` for divergence-free periodic `u` (Lemma 1 of the note; the
  cancellation is machine-verified in exact integer arithmetic by
  `src/ns_certificate_lab/galerkin_obstruction.py`),
* `A` is the viscous part and satisfies `⟪x, A x⟫ ≤ 0`,
* `g` carries the projected force, with `‖g(t)‖ ≤ F t`.

This file proves that such a trajectory obeys

  `‖c(t)‖ ≤ ‖c(0)‖ + ∫₀ᵗ F`,                                        (F-6)

so its norm cannot blow up in finite time.  Since all norms on a
finite-dimensional space are equivalent, that removes the whole ansatz class
from Track F.

Representation choices.
* The state space is an arbitrary real inner product space `E`, not a concrete
  `EuclideanSpace`.  Nothing in the argument uses finite-dimensionality: the
  bound is an energy estimate, and it is the *use* made of it (equivalence of
  norms) that needs `dim E < ∞`.  Keeping `E` abstract makes the Lean statement
  strictly stronger and keeps it free of any coordinate bookkeeping.
* The hypothesis is stated as the differential inequality
  `⟪u t, u' t⟫ ≤ ‖u t‖ * F t` (`norm_le_of_energy_inequality`), with the
  Galerkin structure supplied separately (`inner_galerkin_le`).  This
  separation matters: the analytic content is Grönwall, the PDE content is the
  two structural identities, and only the latter is specific to Navier–Stokes.
* The differential inequality is integrated through
  `t ↦ √(⟪u t, u t⟫ + ε) - ∫₀ᵗ F` for `ε > 0` and then `ε ↓ 0`.  The
  regularization is not cosmetic: `‖·‖` is not differentiable at `0`, and the
  trajectory is allowed to vanish.

Scope caveat.  This file is a statement about a *given* curve `u` in a *given*
inner product space.  It does not construct the Navier–Stokes Galerkin system,
does not prove `⟪x, B x x⟫ = 0` for the actual advection term (that is Lemma 1,
proved on paper and machine-checked in exact arithmetic outside Lean), and does
not perform the ODE continuation of Theorem 1(iii) (obligation F-7).  It is
therefore not connected to `ClayStatement.lean` and proves nothing about the
Millennium problem.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib

namespace NSSingularity

noncomputable section

open Set MeasureTheory Filter Topology RealInnerProductSpace

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-! ### The two structural hypotheses -/

/-- **Energy neutrality** of a quadratic nonlinearity: the state is orthogonal
to the value the nonlinearity takes at it.  For the Galerkin truncation of the
incompressible Navier–Stokes advection term this is exactly
`∫ u·(u·∇)u = 0`, which holds for every divergence-free periodic field. -/
def EnergyNeutral (B : E → E → E) : Prop := ∀ x : E, ⟪x, B x x⟫ = 0

/-- **Dissipativity** of a linear part: `⟪x, A x⟫ ≤ 0`.  For `A = νΔ` restricted
to a space of trigonometric polynomials this is `-ν‖∇u‖² ≤ 0`. -/
def Dissipative (A : E → E) : Prop := ∀ x : E, ⟪x, A x⟫ ≤ 0

/-! ### The differential inequality -/

/-- **F-6, analytic core.**  A curve whose radial energy production is bounded
by `‖u‖ · F` grows at most like `∫F`:

  `‖u b‖ ≤ ‖u 0‖ + ∫₀ᵇ F`.

No inner-product structure beyond `⟪u, u'⟫ ≤ ‖u‖ F` is used, and no
finite-dimensionality; nonnegativity of `F` is not assumed either, since it
follows on `Ioo 0 b` from the hypothesis at points where `u` vanishes and is
otherwise not needed. -/
theorem norm_le_of_energy_inequality
    {u u' : ℝ → E} {F : ℝ → ℝ} {b : ℝ} (hb : 0 ≤ b)
    (hu : ContinuousOn u (Icc 0 b))
    (hu' : ∀ t ∈ Ioo 0 b, HasDerivAt u (u' t) t)
    (hF : ContinuousOn F (Icc 0 b))
    (hF₀ : ∀ t ∈ Icc 0 b, 0 ≤ F t)
    (hkey : ∀ t ∈ Ioo 0 b, ⟪u t, u' t⟫ ≤ ‖u t‖ * F t) :
    ‖u b‖ ≤ ‖u 0‖ + ∫ s in (0 : ℝ)..b, F s := by
  have huIcc : uIcc (0 : ℝ) b = Icc 0 b := uIcc_of_le hb
  have hFint : IntegrableOn F (Icc 0 b) volume := hF.integrableOn_Icc
  set G : ℝ → ℝ := fun t => ∫ s in (0 : ℝ)..t, F s with hGdef
  have hGcont : ContinuousOn G (Icc 0 b) := by
    have h := intervalIntegral.continuousOn_primitive_interval
      (a := (0 : ℝ)) (b := b) (μ := volume) (f := F) (by rwa [huIcc])
    rwa [huIcc] at h
  have hGderiv : ∀ t ∈ Ioo (0 : ℝ) b, HasDerivAt G (F t) t := by
    intro t ht
    have hmem : Icc (0 : ℝ) b ∈ 𝓝 t := Icc_mem_nhds ht.1 ht.2
    have hcont : ContinuousAt F t := hF.continuousAt hmem
    have hint : IntervalIntegrable F volume 0 t := by
      refine (hFint.mono_set ?_).intervalIntegrable
      rw [uIcc_of_le ht.1.le]
      exact Icc_subset_Icc le_rfl ht.2.le
    have hmeas : StronglyMeasurableAtFilter F (𝓝 t) volume :=
      (hF.mono Ioo_subset_Icc_self).stronglyMeasurableAtFilter isOpen_Ioo t ht
    exact intervalIntegral.integral_hasDerivAt_right hint hmeas hcont
  -- The regularized comparison function, one `ε` at a time.
  have step : ∀ ε : ℝ, 0 < ε → ‖u b‖ ≤ Real.sqrt (‖u 0‖ ^ 2 + ε) + G b := by
    intro ε hε
    set q : ℝ → ℝ := fun t => ⟪u t, u t⟫ + ε with hqdef
    have hqpos : ∀ t, 0 < q t := fun t => by
      have := real_inner_self_nonneg (x := u t)
      simp only [hqdef]
      linarith
    set y : ℝ → ℝ := fun t => Real.sqrt (q t) - G t with hydef
    have hqcont : ContinuousOn q (Icc 0 b) := (hu.inner hu).add continuousOn_const
    have hycont : ContinuousOn y (Icc 0 b) :=
      (Real.continuous_sqrt.comp_continuousOn hqcont).sub hGcont
    have hyderiv : ∀ t ∈ Ioo (0 : ℝ) b,
        HasDerivAt y (2 * ⟪u t, u' t⟫ / (2 * Real.sqrt (q t)) - F t) t := by
      intro t ht
      have hq : HasDerivAt q (2 * ⟪u t, u' t⟫) t := by
        have h := ((hu' t ht).inner ℝ (hu' t ht)).add_const ε
        have hcomm : ⟪u t, u' t⟫ + ⟪u' t, u t⟫ = 2 * ⟪u t, u' t⟫ := by
          rw [real_inner_comm (u' t) (u t)]; ring
        rwa [hcomm] at h
      exact (hq.sqrt (hqpos t).ne').sub (hGderiv t ht)
    have hynonpos : ∀ t ∈ Ioo (0 : ℝ) b,
        2 * ⟪u t, u' t⟫ / (2 * Real.sqrt (q t)) - F t ≤ 0 := by
      intro t ht
      have hs : 0 < Real.sqrt (q t) := Real.sqrt_pos.mpr (hqpos t)
      have hnorm : ‖u t‖ ≤ Real.sqrt (q t) := by
        have hsq : ‖u t‖ ^ 2 ≤ q t := by
          simp only [hqdef, real_inner_self_eq_norm_sq]
          linarith
        calc ‖u t‖ = Real.sqrt (‖u t‖ ^ 2) := (Real.sqrt_sq (norm_nonneg _)).symm
          _ ≤ Real.sqrt (q t) := Real.sqrt_le_sqrt hsq
      have hFt : 0 ≤ F t := hF₀ t (Ioo_subset_Icc_self ht)
      have h1 : 2 * ⟪u t, u' t⟫ ≤ 2 * (Real.sqrt (q t) * F t) := by
        have := hkey t ht
        nlinarith [mul_le_mul_of_nonneg_right hnorm hFt]
      have h2 : 2 * ⟪u t, u' t⟫ / (2 * Real.sqrt (q t)) ≤ F t := by
        rw [div_le_iff₀ (by linarith)]
        nlinarith
      linarith
    have hanti : AntitoneOn y (Icc 0 b) := by
      refine antitoneOn_of_deriv_nonpos (convex_Icc 0 b) hycont ?_ ?_
      · intro t ht
        rw [interior_Icc] at ht
        exact (hyderiv t ht).differentiableAt.differentiableWithinAt
      · intro t ht
        rw [interior_Icc] at ht
        rw [(hyderiv t ht).deriv]
        exact hynonpos t ht
    have hy : y b ≤ y 0 := hanti (left_mem_Icc.mpr hb) (right_mem_Icc.mpr hb) hb
    have hG0 : G 0 = 0 := by simp [hGdef]
    have hqb : q b = ‖u b‖ ^ 2 + ε := by simp [hqdef]
    have hq0 : q 0 = ‖u 0‖ ^ 2 + ε := by simp [hqdef]
    have hub : ‖u b‖ ≤ Real.sqrt (q b) := by
      calc ‖u b‖ = Real.sqrt (‖u b‖ ^ 2) := (Real.sqrt_sq (norm_nonneg _)).symm
        _ ≤ Real.sqrt (q b) := Real.sqrt_le_sqrt (by rw [hqb]; linarith)
    simp only [hydef, hG0, sub_zero, hq0] at hy
    linarith
  -- Now let `ε ↓ 0` along `ε = δ²`.
  refine le_of_forall_pos_le_add fun δ hδ => ?_
  have h := step (δ ^ 2) (by positivity)
  have hle : Real.sqrt (‖u 0‖ ^ 2 + δ ^ 2) ≤ ‖u 0‖ + δ := by
    calc Real.sqrt (‖u 0‖ ^ 2 + δ ^ 2)
        ≤ Real.sqrt ((‖u 0‖ + δ) ^ 2) := by
          refine Real.sqrt_le_sqrt ?_
          nlinarith [norm_nonneg (u 0), hδ.le]
      _ = ‖u 0‖ + δ := Real.sqrt_sq (by positivity)
  linarith

/-! ### The Galerkin structure supplies the hypothesis -/

/-- The energy production of an energy-neutral, dissipative quadratic ODE is
controlled by the forcing alone:
`⟪x, g + B x x + A x⟫ ≤ ‖x‖ ‖g‖`.  This is the only place where the
Navier–Stokes structure enters. -/
theorem inner_galerkin_le
    {B : E → E → E} {A : E → E} (hB : EnergyNeutral B) (hA : Dissipative A)
    (x g : E) : ⟪x, g + B x x + A x⟫ ≤ ‖x‖ * ‖g‖ := by
  have hsplit : ⟪x, g + B x x + A x⟫ = ⟪x, g⟫ + ⟪x, B x x⟫ + ⟪x, A x⟫ := by
    rw [inner_add_right, inner_add_right]
  rw [hsplit, hB x]
  have h1 : ⟪x, g⟫ ≤ ‖x‖ * ‖g‖ := real_inner_le_norm x g
  have h2 := hA x
  linarith

/-- **F-6, packaged.**  A trajectory of the energy-neutral dissipative
quadratic system `u' = g + B(u,u) + A u` satisfies
`‖u b‖ ≤ ‖u 0‖ + ∫₀ᵇ F` whenever `‖g t‖ ≤ F t`.

Nonnegativity of `F` is not a hypothesis: it follows from `‖g t‖ ≤ F t`. -/
theorem galerkin_norm_le
    {B : E → E → E} {A : E → E} {g u : ℝ → E} {F : ℝ → ℝ} {b : ℝ} (hb : 0 ≤ b)
    (hB : EnergyNeutral B) (hA : Dissipative A)
    (hg : ∀ t, ‖g t‖ ≤ F t)
    (hu : ContinuousOn u (Icc 0 b))
    (hode : ∀ t ∈ Ioo (0 : ℝ) b,
      HasDerivAt u (g t + B (u t) (u t) + A (u t)) t)
    (hF : ContinuousOn F (Icc 0 b)) :
    ‖u b‖ ≤ ‖u 0‖ + ∫ s in (0 : ℝ)..b, F s := by
  refine norm_le_of_energy_inequality hb hu hode hF
    (fun t _ => le_trans (norm_nonneg (g t)) (hg t)) ?_
  intro t ht
  refine le_trans (inner_galerkin_le hB hA (u t) (g t)) ?_
  exact mul_le_mul_of_nonneg_left (hg t) (norm_nonneg _)

/-- **F-6, uniform bound.**  The same trajectory is bounded on the whole of
`[0,T]` by the single constant `‖u 0‖ + ∫₀ᵀ F`. -/
theorem galerkin_norm_le_of_mem
    {B : E → E → E} {A : E → E} {g u : ℝ → E} {F : ℝ → ℝ} {T : ℝ} (hT : 0 ≤ T)
    (hB : EnergyNeutral B) (hA : Dissipative A)
    (hg : ∀ t, ‖g t‖ ≤ F t)
    (hu : ContinuousOn u (Icc 0 T))
    (hode : ∀ t ∈ Ioo (0 : ℝ) T,
      HasDerivAt u (g t + B (u t) (u t) + A (u t)) t)
    (hF : ContinuousOn F (Icc 0 T)) :
    ∀ t ∈ Icc (0 : ℝ) T, ‖u t‖ ≤ ‖u 0‖ + ∫ s in (0 : ℝ)..T, F s := by
  intro t ht
  have hsub : Icc (0 : ℝ) t ⊆ Icc 0 T := Icc_subset_Icc le_rfl ht.2
  have hlocal : ‖u t‖ ≤ ‖u 0‖ + ∫ s in (0 : ℝ)..t, F s :=
    galerkin_norm_le ht.1 hB hA hg (hu.mono hsub)
      (fun s hs => hode s ⟨hs.1, lt_of_lt_of_le hs.2 ht.2⟩) (hF.mono hsub)
  have hfi : IntervalIntegrable F volume 0 T := by
    have hInt : IntegrableOn F (uIcc (0 : ℝ) T) volume := by
      rw [uIcc_of_le hT]
      exact hF.integrableOn_Icc
    exact hInt.intervalIntegrable
  have hmono : (∫ s in (0 : ℝ)..t, F s) ≤ ∫ s in (0 : ℝ)..T, F s :=
    intervalIntegral.integral_mono_interval le_rfl ht.1 ht.2
      (Eventually.of_forall fun s => le_trans (norm_nonneg (g s)) (hg s)) hfi
  linarith

/-- **F-6, no finite-time blow-up.**  Under the same hypotheses the norm of the
trajectory does not tend to `+∞` as `t → T` from the left.  This is the
statement a Track-F ansatz would have to contradict. -/
theorem galerkin_not_tendsto_atTop
    {B : E → E → E} {A : E → E} {g u : ℝ → E} {F : ℝ → ℝ} {T : ℝ} (hT : 0 < T)
    (hB : EnergyNeutral B) (hA : Dissipative A)
    (hg : ∀ t, ‖g t‖ ≤ F t)
    (hu : ContinuousOn u (Icc 0 T))
    (hode : ∀ t ∈ Ioo (0 : ℝ) T,
      HasDerivAt u (g t + B (u t) (u t) + A (u t)) t)
    (hF : ContinuousOn F (Icc 0 T)) :
    ¬ Tendsto (fun t => ‖u t‖) (𝓝[<] T) atTop := by
  intro hdiv
  set C : ℝ := ‖u 0‖ + ∫ s in (0 : ℝ)..T, F s with hC
  have hbound := galerkin_norm_le_of_mem hT.le hB hA hg hu hode hF
  have hev : ∀ᶠ t in 𝓝[<] T, C < ‖u t‖ := hdiv.eventually_gt_atTop C
  have hnear : ∀ᶠ t in 𝓝[<] T, t ∈ Ioo (0 : ℝ) T := by
    have h1 : ∀ᶠ t in 𝓝[<] T, (0 : ℝ) < t :=
      eventually_nhdsWithin_of_eventually_nhds (eventually_gt_nhds hT)
    have h2 : ∀ᶠ t in 𝓝[<] T, t < T := eventually_mem_nhdsWithin.mono fun _ h => h
    exact h1.and h2
  obtain ⟨t, htC, htmem⟩ := (hev.and hnear).exists
  exact absurd (hbound t (Ioo_subset_Icc_self htmem)) (not_le.mpr htC)

end

end NSSingularity
