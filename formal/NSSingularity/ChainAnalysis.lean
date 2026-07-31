/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import Mathlib
import NSSingularity.ControlODE

/-!
# Chain analysis — the finite-statement halves of the EXT audit

The EXT audit lane replaces the blanket external hypotheses EXT-P1/P2/P3 of
the Track-P chain by narrower, named obligations.  Two of those obligations
have a *scalar, finite-statement* half that lives entirely inside real
analysis on an interval, and this file proves those halves unconditionally:

* **EXT-P2's scalar half** — the comparison of an *integral-inequality*
  subsolution with a genuine ODE supersolution.  The `H³` energy argument
  produces, for the scalar error size `φ`, only the integral inequality
  `φ t ≤ φ s + ∫ₛᵗ f(φ τ) dτ` (no differentiability of `φ` is available a
  priori), while the certificate's majorant `R` solves `R' = f(R)` exactly.
  `integral_comparison` shows `φ ≤ R` on the slab, and
  `integral_riccati_comparison` instantiates it for the Riccati field
  `f r = a r + b r² + e` of the control inequality, manufacturing the
  one-sided Lipschitz constant from compactness.
* **EXT-P3's gluing and endpoint logic** — `glued_continuous` (two continuous
  pieces agreeing at the joint give a continuous piecewise function on the
  union interval) and the endpoint extension: a function on `[t₀, T)` into a
  complete metric space whose values are Cauchy along `T` extends to a
  continuous function on `[t₀, T]` (`exists_continuousOn_Icc_extension`),
  with a uniform-modulus sufficient condition
  (`cauchy_map_of_uniform_modulus`, combined form
  `exists_continuousOn_Icc_extension_of_modulus`).
* **The discharge shape** — `cond_to_uncond`, the propositional record that
  closing every hypothesis of a conditional payload converts the implication
  into its conclusion.  Trivial, and stated so that the certificate-dependency
  step is a named, audited theorem rather than prose.

## The mathematics `integral_comparison` does and does not do

The proof upgrades the integral inequality to a right-Dini slope bound — for
`r > f(φ x)`, continuity of `f ∘ φ` bounds the integrand by `r` on a small
right interval, so `(φ z - φ x)/(z - x) < r` there — and then fences `φ`
under the strict supersolutions `R + δ e^{(L+1)(t-t₀)}` using mathlib's
`image_le_of_liminf_slope_right_lt_deriv_boundary'`, exactly the
strict-boundary trick of `ControlODE.nonpos_of_deriv_le_mul_of_pos`.  The
one-sided Lipschitz hypothesis is consumed only at contact points, and only
for the pair `(φ x, R x)` — which is why the Riccati instance can discharge
it from compactness bounds on the two trajectories.

**Retreat, stated precisely.**  The lemma hypothesises the integral
inequality *between every pair* `s ≤ t` in the slab, not merely from the left
endpoint.  This is what the energy argument actually provides (it integrates
the estimate over an arbitrary subinterval) and it is what the Dini upgrade
needs; a version assuming only `φ t ≤ φ t₀ + ∫_{t₀}^t f(φ)` is *not* proved
here, and no claim is made about it.  No absolute continuity, a.e.
differentiability, or measure-theoretic apparatus is assumed of `φ` beyond
continuity — the integral inequality itself carries the regularity.

## What is not formalised, and why

Supplying the integral inequality for the true Navier–Stokes error, the
uniform modulus for the `H³` trajectory near a continuation time, and the
identification of the extended function with a solution are PDE statements:
they are exactly the analytic halves of the EXT obligations, they are
recorded in the research notes, and they are **never** Lean axioms.  Nothing
here mentions Navier–Stokes; every theorem is an unconditional statement
about real functions, metric spaces, or propositions.  Nothing in this file,
alone or in combination with the rest of the library, is a statement about
the Clay problem.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/

namespace NSSingularity

noncomputable section

open Set Filter MeasureTheory intervalIntegral
open scoped Topology

/-! ## A1 — the integral-inequality comparison (EXT-P2, scalar half) -/

section IntegralComparison

/-- **A1, the general comparison.**  Let `φ` be continuous on `[t₀, t₁]` and
satisfy the *integral* inequality `φ t ≤ φ s + ∫ₛᵗ F(φ τ) dτ` for all
`t₀ ≤ s ≤ t ≤ t₁`, let `R` solve `R' = F(R)` (right derivative on
`[t₀, t₁)`), and let `F` be one-sided Lipschitz along the pair, i.e.
`F(φ t) - F(R t) ≤ L (φ t - R t)` whenever `R t ≤ φ t`.  If `φ t₀ ≤ R t₀`
then `φ ≤ R` on the whole slab.

The subsolution is *not* assumed differentiable anywhere: the integral
inequality plus continuity of `F ∘ φ` yields the right-Dini slope bound that
mathlib's fencing theorem consumes, and the strict supersolutions
`R + δ e^{(L+1)(t-t₀)}` (for `δ → 0⁺`) buy the strictness the fencing
demands.  The Lipschitz hypothesis is only ever used at contact points, where
`R t ≤ φ t` holds automatically. -/
theorem integral_comparison {φ R F : ℝ → ℝ} {L t0 t1 : ℝ}
    (hFc : Continuous F)
    (hφ : ContinuousOn φ (Icc t0 t1))
    (hR : ContinuousOn R (Icc t0 t1))
    (hR' : ∀ t ∈ Ico t0 t1, HasDerivWithinAt R (F (R t)) (Ici t) t)
    (hlip : ∀ t ∈ Ico t0 t1, R t ≤ φ t → F (φ t) - F (R t) ≤ L * (φ t - R t))
    (h0 : φ t0 ≤ R t0)
    (hint : ∀ s t : ℝ, t0 ≤ s → s ≤ t → t ≤ t1 →
      φ t ≤ φ s + ∫ τ in s..t, F (φ τ)) :
    ∀ t ∈ Icc t0 t1, φ t ≤ R t := by
  have hci : ContinuousOn (fun τ => F (φ τ)) (Icc t0 t1) :=
    hFc.comp_continuousOn hφ
  -- The Dini upgrade: the integral inequality bounds every right slope.
  have hslope : ∀ x ∈ Ico t0 t1, ∀ r, F (φ x) < r →
      ∃ᶠ z in 𝓝[>] x, slope φ x z < r := by
    intro x hx r hr
    have hr'lt : F (φ x) < (F (φ x) + r) / 2 := by linarith
    have hr'r : (F (φ x) + r) / 2 < r := by linarith
    have hcw : ContinuousWithinAt (fun τ => F (φ τ)) (Icc t0 t1) x :=
      hci x (Ico_subset_Icc_self hx)
    have h1 : ∀ᶠ τ in 𝓝[Icc t0 t1] x, F (φ τ) < (F (φ x) + r) / 2 :=
      hcw.eventually_lt_const hr'lt
    have h2 : 𝓝[Icc x t1] x ≤ 𝓝[Icc t0 t1] x :=
      nhdsWithin_mono x (Icc_subset_Icc hx.1 le_rfl)
    have h3 : ∀ᶠ τ in 𝓝[Icc x t1] x, F (φ τ) < (F (φ x) + r) / 2 :=
      h1.filter_mono h2
    rw [nhdsWithin_Icc_eq_nhdsGE hx.2] at h3
    obtain ⟨u, hu, hsub⟩ := mem_nhdsGE_iff_exists_Ico_subset.mp h3
    have hxv : x < min u t1 := lt_min hu hx.2
    have hall : ∀ z ∈ Ioo x (min u t1), slope φ x z < r := by
      intro z hz
      have hxz : x < z := hz.1
      have hzt1 : z ≤ t1 := le_of_lt (lt_of_lt_of_le hz.2 (min_le_right u t1))
      have hzu : z < u := lt_of_lt_of_le hz.2 (min_le_left u t1)
      have hmem : ∀ τ ∈ Icc x z, F (φ τ) ≤ (F (φ x) + r) / 2 := fun τ hτ =>
        le_of_lt (hsub ⟨hτ.1, lt_of_le_of_lt hτ.2 hzu⟩)
      have hInt : IntervalIntegrable (fun τ => F (φ τ)) volume x z := by
        apply ContinuousOn.intervalIntegrable
        rw [uIcc_of_le hxz.le]
        exact hci.mono (Icc_subset_Icc hx.1 hzt1)
      have hmono : (∫ τ in x..z, F (φ τ)) ≤ ∫ τ in x..z, (F (φ x) + r) / 2 :=
        intervalIntegral.integral_mono_on hxz.le hInt intervalIntegrable_const hmem
      have hconst : (∫ τ in x..z, (F (φ x) + r) / 2) = (z - x) * ((F (φ x) + r) / 2) := by
        rw [intervalIntegral.integral_const, smul_eq_mul]
      have hφz : φ z - φ x ≤ (z - x) * ((F (φ x) + r) / 2) := by
        have hI := hint x z hx.1 hxz.le hzt1
        rw [hconst] at hmono
        linarith
      rw [slope_def_field, div_lt_iff₀ (by linarith : (0 : ℝ) < z - x)]
      nlinarith [mul_lt_mul_of_pos_right hr'r (by linarith : (0 : ℝ) < z - x)]
    have hIoo : Ioo x (min u t1) ∈ 𝓝[>] x := Ioo_mem_nhdsGT hxv
    exact Filter.Eventually.frequently
      (mem_of_superset hIoo fun z hz => hall z hz)
  -- The strict fencing under `R + δ e^{(L+1)(t-t₀)}`, for every `δ > 0`.
  have key : ∀ δ : ℝ, 0 < δ → ∀ t ∈ Icc t0 t1,
      φ t ≤ R t + δ * Real.exp ((L + 1) * (t - t0)) := by
    intro δ hδ t ht
    have hBderiv : ∀ y : ℝ,
        HasDerivAt (fun s : ℝ => δ * Real.exp ((L + 1) * (s - t0)))
          ((L + 1) * (δ * Real.exp ((L + 1) * (y - t0)))) y := by
      intro y
      have h1 : HasDerivAt (fun s : ℝ => (L + 1) * (s - t0)) (L + 1) y := by
        simpa using ((hasDerivAt_id y).sub_const t0).const_mul (L + 1)
      have h2 := h1.exp.const_mul δ
      have heq : δ * (Real.exp ((L + 1) * (y - t0)) * (L + 1))
          = (L + 1) * (δ * Real.exp ((L + 1) * (y - t0))) := by ring
      rw [← heq]
      exact h2
    refine image_le_of_liminf_slope_right_lt_deriv_boundary'
      (f := φ) (f' := fun y => F (φ y))
      (B := fun y => R y + δ * Real.exp ((L + 1) * (y - t0)))
      (B' := fun y => F (R y) + (L + 1) * (δ * Real.exp ((L + 1) * (y - t0))))
      hφ hslope ?_ ?_ ?_ ?_ ht
    · have h1 : δ * Real.exp ((L + 1) * (t0 - t0)) = δ := by simp
      rw [h1]
      linarith
    · exact hR.add (Continuous.continuousOn (by fun_prop))
    · intro y hy
      exact (hR' y hy).add ((hBderiv y).hasDerivWithinAt)
    · intro y hy hcontact
      have hE : 0 < Real.exp ((L + 1) * (y - t0)) := Real.exp_pos _
      have hδE : 0 < δ * Real.exp ((L + 1) * (y - t0)) := by positivity
      have hRy : R y ≤ φ y := by rw [hcontact]; linarith
      have hlipy := hlip y hy hRy
      have hgap : φ y - R y = δ * Real.exp ((L + 1) * (y - t0)) := by
        rw [hcontact]; ring
      rw [hgap] at hlipy
      have hexpand : (L + 1) * (δ * Real.exp ((L + 1) * (y - t0)))
          = L * (δ * Real.exp ((L + 1) * (y - t0)))
            + δ * Real.exp ((L + 1) * (y - t0)) := by ring
      rw [hexpand]
      linarith
  -- `δ → 0⁺`.
  intro t ht
  rcases le_or_gt (φ t) (R t) with hle | hpos
  · exact hle
  exfalso
  have hEpos : 0 < Real.exp ((L + 1) * (t - t0)) := Real.exp_pos _
  have h1 := key ((φ t - R t) / (2 * Real.exp ((L + 1) * (t - t0))))
    (div_pos (by linarith) (by positivity)) t ht
  have h2 : (φ t - R t) / (2 * Real.exp ((L + 1) * (t - t0)))
      * Real.exp ((L + 1) * (t - t0)) = (φ t - R t) / 2 := by
    field_simp
  rw [h2] at h1
  linarith

/-- **A1, the Riccati instance.**  The target field of the control
inequality: `f r = a r + b r² + e` with `0 ≤ b` (no sign condition on `a` or
`e` is needed for the *comparison*; the certificate's `e ≥ 0` matters only
for the size of the majorant, not for its validity).  An integral-inequality
subsolution `φ` never overtakes an exact solution `R` of the Riccati ODE that
starts above it.  The one-sided Lipschitz constant
`L = |a| + b (C_φ + C_R)` is manufactured from compactness bounds on the two
trajectories, exactly as in `ControlODE.riccati_comparison`; as there, the
constant exists but is not effective, so this theorem certifies the ordering
and never a quantitative rate. -/
theorem integral_riccati_comparison {φ R : ℝ → ℝ} {a b e t0 t1 : ℝ}
    (hb : 0 ≤ b)
    (hφ : ContinuousOn φ (Icc t0 t1))
    (hR : ContinuousOn R (Icc t0 t1))
    (hR' : ∀ t ∈ Ico t0 t1,
      HasDerivWithinAt R (a * R t + b * R t ^ 2 + e) (Ici t) t)
    (h0 : φ t0 ≤ R t0)
    (hint : ∀ s t : ℝ, t0 ≤ s → s ≤ t → t ≤ t1 →
      φ t ≤ φ s + ∫ τ in s..t, (a * φ τ + b * φ τ ^ 2 + e)) :
    ∀ t ∈ Icc t0 t1, φ t ≤ R t := by
  obtain ⟨Cp, hCp⟩ := isCompact_Icc.exists_bound_of_continuousOn hφ
  obtain ⟨CR, hCR⟩ := isCompact_Icc.exists_bound_of_continuousOn hR
  refine integral_comparison (F := fun r => a * r + b * r ^ 2 + e)
    (L := |a| + b * (Cp + CR)) (by fun_prop) hφ hR hR' ?_ h0 hint
  intro y hy hle
  have hyI : y ∈ Icc t0 t1 := Ico_subset_Icc_self hy
  have hpb : φ y ≤ Cp :=
    (le_abs_self _).trans (by simpa [Real.norm_eq_abs] using hCp y hyI)
  have hRb : R y ≤ CR :=
    (le_abs_self _).trans (by simpa [Real.norm_eq_abs] using hCR y hyI)
  have hDG : 0 ≤ (φ y - R y) * ((|a| - a) + b * ((Cp + CR) - (φ y + R y))) :=
    mul_nonneg (by linarith)
      (add_nonneg (by linarith [le_abs_self a]) (mul_nonneg hb (by linarith)))
  nlinarith [hDG]

end IntegralComparison

/-! ## A2 — gluing and endpoint extension (EXT-P3, finite halves) -/

section GluingExtension

/-- **A2 (i), gluing.**  Two continuous pieces that agree at the joint glue to
a continuous piecewise function on the union interval: `u₁` continuous on
`[t₀, t₁]`, `u₂` continuous on `[t₁, t₂]`, `u₁ t₁ = u₂ t₁`, then
`t ↦ if t ≤ t₁ then u₁ t else u₂ t` is continuous on `[t₀, t₂]`.  This is the
piecewise-centre shape of `TrackPChain.two_slab_composition`, upgraded from a
pointwise bound to continuity of the glued trajectory — the topological half
of the continuation step, with mathlib's `ContinuousOn.if` doing the frontier
bookkeeping at `t₁`. -/
theorem glued_continuous {E : Type*} [TopologicalSpace E]
    {u1 u2 : ℝ → E} {t0 t1 t2 : ℝ}
    (h01 : t0 ≤ t1) (h12 : t1 ≤ t2)
    (hu1 : ContinuousOn u1 (Icc t0 t1)) (hu2 : ContinuousOn u2 (Icc t1 t2))
    (hglue : u1 t1 = u2 t1) :
    ContinuousOn (fun t => if t ≤ t1 then u1 t else u2 t) (Icc t0 t2) := by
  refine ContinuousOn.if ?_ ?_ ?_
  · rintro a ⟨-, hfr⟩
    have hfr' : a ∈ frontier (Iic t1) := hfr
    rw [frontier_Iic] at hfr'
    obtain rfl : a = t1 := hfr'
    exact hglue
  · have hcl : closure {a : ℝ | a ≤ t1} = Iic t1 := closure_Iic t1
    have hset : Icc t0 t2 ∩ closure {a : ℝ | a ≤ t1} = Icc t0 t1 := by
      rw [hcl]
      ext x
      simp only [mem_inter_iff, mem_Icc, mem_Iic]
      constructor
      · rintro ⟨⟨hx0, -⟩, hx1⟩
        exact ⟨hx0, hx1⟩
      · rintro ⟨hx0, hx1⟩
        exact ⟨⟨hx0, hx1.trans h12⟩, hx1⟩
    rw [hset]
    exact hu1
  · have hcl : closure {a : ℝ | ¬a ≤ t1} = Ici t1 := by
      have h : {a : ℝ | ¬a ≤ t1} = Ioi t1 := by
        ext x
        simp [not_le]
      rw [h, closure_Ioi]
    have hset : Icc t0 t2 ∩ closure {a : ℝ | ¬a ≤ t1} = Icc t1 t2 := by
      rw [hcl]
      ext x
      simp only [mem_inter_iff, mem_Icc, mem_Ici]
      constructor
      · rintro ⟨⟨-, hx2⟩, hx1⟩
        exact ⟨hx1, hx2⟩
      · rintro ⟨hx1, hx2⟩
        exact ⟨⟨h01.trans hx1, hx2⟩, hx1⟩
    rw [hset]
    exact hu2

/-- **A2 (ii), the uniform modulus gives the Cauchy property.**  If `u` has a
uniform continuity modulus on `[t₀, T)` — for every `ε > 0` a `δ > 0` with
`dist (u s) (u t) < ε` whenever `|s - t| < δ` — then the image filter of `u`
along `T⁻` is Cauchy.  This is the finite bridge from "uniformly continuous
on the half-open slab" to the completeness argument of
`exists_continuousOn_Icc_extension`: the interval `(T - δ, T)` is mapped into
an `ε`-set. -/
theorem cauchy_map_of_uniform_modulus {E : Type*} [MetricSpace E]
    {u : ℝ → E} {t0 T : ℝ} (h : t0 < T)
    (hmod : ∀ ε > 0, ∃ δ > 0, ∀ s ∈ Ico t0 T, ∀ t ∈ Ico t0 T,
      |s - t| < δ → dist (u s) (u t) < ε) :
    Cauchy (Filter.map u (𝓝[<] T)) := by
  rw [Metric.cauchy_iff]
  refine ⟨Filter.map_neBot, ?_⟩
  intro ε hε
  obtain ⟨δ, hδ, hmod'⟩ := hmod ε hε
  have hδ' : 0 < min δ (T - t0) := lt_min hδ (by linarith)
  refine ⟨u '' Ioo (T - min δ (T - t0)) T,
    Filter.image_mem_map (Ioo_mem_nhdsLT (by linarith)), ?_⟩
  rintro x ⟨s, hs, rfl⟩ y ⟨t, ht, rfl⟩
  have hmin1 : min δ (T - t0) ≤ δ := min_le_left _ _
  have hmin2 : min δ (T - t0) ≤ T - t0 := min_le_right _ _
  have hsI : s ∈ Ico t0 T := ⟨by linarith [hs.1], hs.2⟩
  have htI : t ∈ Ico t0 T := ⟨by linarith [ht.1], ht.2⟩
  refine hmod' s hsI t htI ?_
  rw [abs_sub_lt_iff]
  constructor <;> linarith [hs.1, hs.2, ht.1, ht.2]

/-- **A2 (iii), the endpoint extension.**  A continuous function on the
half-open slab `[t₀, T)` into a complete metric space whose values are Cauchy
along `T⁻` extends to a continuous function on the closed slab `[t₀, T]`,
agreeing with the original on `[t₀, T)`.  Completeness supplies the endpoint
value, and mathlib's `extendFrom`/`continuousOn_extendFrom` supply the
extension; nothing here says what equation, if any, the extended function
satisfies — identifying the extension with a solution is the analytic half of
EXT-P3 and stays outside Lean. -/
theorem exists_continuousOn_Icc_extension {E : Type*} [MetricSpace E]
    [CompleteSpace E] {u : ℝ → E} {t0 T : ℝ} (h : t0 < T)
    (hu : ContinuousOn u (Ico t0 T))
    (hcauchy : Cauchy (Filter.map u (𝓝[<] T))) :
    ∃ v : ℝ → E, ContinuousOn v (Icc t0 T) ∧ EqOn u v (Ico t0 T) := by
  obtain ⟨y, hy⟩ := CompleteSpace.complete hcauchy
  have hyT : Filter.Tendsto u (𝓝[<] T) (𝓝 y) := hy
  refine ⟨extendFrom (Ico t0 T) u, ?_,
    fun x hx => (extendFrom_extends hu x hx).symm⟩
  refine continuousOn_extendFrom ?_ ?_
  · rw [closure_Ico h.ne]
  · intro x hx
    rcases lt_or_eq_of_le hx.2 with hlt | heq
    · exact ⟨u x, hu x ⟨hx.1, hlt⟩⟩
    · subst heq
      exact ⟨y, hyT.mono_left (nhdsWithin_mono x Ico_subset_Iio_self)⟩

/-- **A2 (iv), the combined form.**  Uniform modulus on `[t₀, T)` plus
completeness of the target give the continuous extension to `[t₀, T]`
directly — (ii) composed with (iii). -/
theorem exists_continuousOn_Icc_extension_of_modulus {E : Type*}
    [MetricSpace E] [CompleteSpace E] {u : ℝ → E} {t0 T : ℝ} (h : t0 < T)
    (hu : ContinuousOn u (Ico t0 T))
    (hmod : ∀ ε > 0, ∃ δ > 0, ∀ s ∈ Ico t0 T, ∀ t ∈ Ico t0 T,
      |s - t| < δ → dist (u s) (u t) < ε) :
    ∃ v : ℝ → E, ContinuousOn v (Icc t0 T) ∧ EqOn u v (Ico t0 T) :=
  exists_continuousOn_Icc_extension h hu (cauchy_map_of_uniform_modulus h hmod)

end GluingExtension

/-! ## A3 — the discharge shape -/

section Discharge

/-- **A3.**  The certificate-dependency step, as a named theorem: a payload
conditional on three hypotheses, together with the three hypotheses, yields
the payload's conclusion.  Logically this is nothing but `modus ponens`
through a conjunction; it is stated (and audited) so that the *shape* of the
EXT discharge — closing every item of the external list converts the
conditional certificate into its unconditional conclusion — is machine
checked rather than prose.  No instance of `H1`, `H2`, `H3` or `C` is
supplied anywhere in this repository's Lean code: in particular the EXT
obligations themselves are recorded in the research notes and are never Lean
axioms, so this theorem discharges nothing by itself. -/
theorem cond_to_uncond {H1 H2 H3 C : Prop}
    (hpayload : H1 ∧ H2 ∧ H3 → C) (h1 : H1) (h2 : H2) (h3 : H3) : C :=
  hpayload ⟨h1, h2, h3⟩

end Discharge

end

end NSSingularity
