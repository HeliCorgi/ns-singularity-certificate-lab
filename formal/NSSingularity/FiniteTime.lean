/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
/-
Stage 1, item **F-2** of `docs/formalization_map.md`:
the **finite physical time** condition for a dynamically rescaled trajectory.

Mathematical content.  A dynamic-rescaling ansatz replaces physical time `t`
by a rescaled time `s`, the two being linked by the scale function `L` through

  `dt/ds = L(s)^2`,   i.e.   `t(s) = t₀ + ∫_{s₀}^{s} L(σ)^2 dσ`.

Rescaled time runs to `+∞` at a putative singularity, so the assertion that
the singularity happens at a *finite* physical time is exactly the assertion
that the improper integral `∫_{s₀}^{∞} L^2` converges.  The obligation
formalized here is therefore

  IF `L^2` is integrable on `[s₀,∞)`
  THEN `t(·)` has a finite limit `T = t₀ + ∫_{(s₀,∞)} L^2` as `s → ∞`,
       `t` is monotone, and `t(s) < T` for every `s`

together with the converse-flavoured statement that if `L^2` is *not*
integrable then `t(s) → ∞`, so integrability is not merely sufficient but
is the exact dividing line.

Representation choices.
* `physicalTime` is an **interval integral** `∫ σ in s₀..s, L σ ^ 2` rather
  than a set integral.  This is the form in which the fundamental theorem of
  calculus is stated in mathlib (`intervalIntegral.integral_hasDerivAt_right`),
  so `hasDerivAt_physicalTime` below can certify that the definition really
  does solve `dt/ds = L²` with `t(s₀) = t₀`; it also fixes the sign for
  `s < s₀` without a case split.
* The limit `blowupTime` is a **set integral over `Set.Ioi s₀`**, because
  that is the form of mathlib's improper-integral convergence theorem
  `MeasureTheory.intervalIntegral_tendsto_integral_Ioi`.  The two forms are
  reconciled by exactly that lemma, which is the single place where the real
  content (dominated convergence for an `AECover`) enters.
* The integrand is written `L σ ^ 2`, so nonnegativity of the rate is free
  (`sq_nonneg`) and **no positivity hypothesis on `L` is needed** for
  monotonicity, for the finite limit, or for the bound `t(s) ≤ T`.
  Positivity of `L` is used only for the *strict* statements
  (`physicalTime_strictMonoOn`, `physicalTime_lt_blowupTime`).
* Lean's Bochner integral is junk-valued (`= 0`) at non-integrable
  integrands.  Every statement below that could be affected by this carries
  the relevant integrability hypothesis explicitly; the two hypotheses used
  are `IntegrableOn (scaleRate L) (Set.Ici s₀)` (the F-2 hypothesis proper)
  and the weaker local form `∀ b, s₀ ≤ b → IntervalIntegrable (scaleRate L)
  volume s₀ b`, which is all that monotonicity needs and which also survives
  in the non-integrable case.

Scope caveat.  This file is a statement about a *given* function `L`.  It
does not construct `L`, does not connect `L` to any Navier–Stokes solution,
and does not verify the integrability hypothesis for any concrete candidate:
supplying `IntegrableOn (fun σ => L σ ^ 2) (Set.Ici s₀)` as a rigorous
(interval-arithmetic) bound is precisely the obligation a future numerical
certificate would have to discharge.  See `formal/README.md`.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib

namespace NSSingularity

noncomputable section

open MeasureTheory Set Filter Topology

variable {t₀ s₀ : ℝ} {L : ℝ → ℝ}

/-! ### The rescaling rate and the physical-time map -/

/-- The rescaling rate `dt/ds = L(σ)²` attached to a scale function `L`. -/
def scaleRate (L : ℝ → ℝ) (σ : ℝ) : ℝ := L σ ^ 2

/-- The rate is nonnegative for *every* `L`: it is a square. -/
theorem scaleRate_nonneg (L : ℝ → ℝ) (σ : ℝ) : 0 ≤ scaleRate L σ := sq_nonneg _

/-- For a strictly positive scale function the rate is strictly positive. -/
theorem scaleRate_pos (hL : ∀ σ, 0 < L σ) (σ : ℝ) : 0 < scaleRate L σ :=
  pow_pos (hL σ) 2

/-- **The physical-time map.**  `t(s) = t₀ + ∫_{s₀}^{s} L(σ)² dσ`, the
solution of `dt/ds = L²` normalized by `t(s₀) = t₀` (see
`hasDerivAt_physicalTime` and `physicalTime_base`). -/
def physicalTime (t₀ s₀ : ℝ) (L : ℝ → ℝ) (s : ℝ) : ℝ :=
  t₀ + ∫ σ in s₀..s, scaleRate L σ

/-- **The candidate blow-up time** `T = t₀ + ∫_{(s₀,∞)} L(σ)² dσ`.  Under the
F-2 hypothesis this is the finite limit of `physicalTime` (see
`tendsto_physicalTime`); without it the set integral is junk-valued and the
statement to use instead is `tendsto_physicalTime_atTop`. -/
def blowupTime (t₀ s₀ : ℝ) (L : ℝ → ℝ) : ℝ :=
  t₀ + ∫ σ in Ioi s₀, scaleRate L σ

@[simp] theorem physicalTime_base (t₀ s₀ : ℝ) (L : ℝ → ℝ) :
    physicalTime t₀ s₀ L s₀ = t₀ := by
  simp [physicalTime]

/-- The definition is the right one: `physicalTime` really does solve the
rescaling ODE `dt/ds = L(s)²`.  (Continuity of `L` is assumed only here; it
is not used anywhere else in this file.) -/
theorem hasDerivAt_physicalTime (hL : Continuous L) (t₀ s₀ s : ℝ) :
    HasDerivAt (physicalTime t₀ s₀ L) (scaleRate L s) s := by
  have hc : Continuous (scaleRate L) := hL.pow 2
  exact (intervalIntegral.integral_hasDerivAt_right (a := s₀) (b := s)
    (hc.intervalIntegrable _ _) (hc.stronglyMeasurableAtFilter _ _)
    hc.continuousAt).const_add t₀

/-! ### Integrability bookkeeping -/

/-- Integrability of the rate on the half-line `[s₀,∞)` gives interval
integrability on every `[s₀,b]` with `b ≥ s₀`.  This is the only form of the
F-2 hypothesis that the monotonicity arguments consume. -/
theorem intervalIntegrable_scaleRate
    (hint : IntegrableOn (scaleRate L) (Ici s₀)) {b : ℝ} (hb : s₀ ≤ b) :
    IntervalIntegrable (scaleRate L) volume s₀ b := by
  refine (hint.mono_set ?_).intervalIntegrable
  rw [uIcc_of_le hb]
  exact fun x hx => hx.1

/-! ### Monotonicity

Both monotonicity statements are proved from the *local* hypothesis
`hint : ∀ b, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b`, which is
implied by the F-2 hypothesis (`intervalIntegrable_scaleRate`) but also holds
in the non-integrable case treated at the end of the file. -/

/-- Forward differences of the physical-time map are nonnegative: the rate is
a square, so `∫_{s₀}^{s}` is nonnegative for `s ≥ s₀` — no integrability is
needed for this direction. -/
theorem le_physicalTime_of_le (t₀ : ℝ) {s : ℝ} (hs : s₀ ≤ s) :
    t₀ ≤ physicalTime t₀ s₀ L s := by
  have := intervalIntegral.integral_nonneg (μ := volume) (f := scaleRate L) hs
    (fun u _ => scaleRate_nonneg L u)
  simpa [physicalTime] using this

/-- Below the base point the physical time never exceeds `t₀`.  (For `s < s₀`
the interval integral is `-∫_{s}^{s₀}`, or the junk value `0` if the rate is
not integrable there; either way it is `≤ 0`.) -/
theorem physicalTime_le_of_le_base (t₀ : ℝ) {s : ℝ} (hs : s ≤ s₀) :
    physicalTime t₀ s₀ L s ≤ t₀ := by
  have h : 0 ≤ ∫ σ in s..s₀, scaleRate L σ :=
    intervalIntegral.integral_nonneg (μ := volume) (f := scaleRate L) hs
      (fun u _ => scaleRate_nonneg L u)
  have hsymm : (∫ σ in s₀..s, scaleRate L σ) = -∫ σ in s..s₀, scaleRate L σ :=
    intervalIntegral.integral_symm s s₀
  simp only [physicalTime, hsymm]
  linarith

/-- **Monotonicity on the rescaled half-line.** -/
theorem physicalTime_monotoneOn
    (hint : ∀ b : ℝ, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b)
    (t₀ : ℝ) : MonotoneOn (physicalTime t₀ s₀ L) (Ici s₀) := by
  intro a ha b hb hab
  simp only [physicalTime, add_le_add_iff_left]
  exact intervalIntegral.integral_mono_interval le_rfl ha hab
    (Eventually.of_forall fun u => scaleRate_nonneg L u) (hint b hb)

/-- The comparison form of monotonicity that is actually used below: it also
covers base points `a < s₀`, where `physicalTime` need not be monotone
because the rate may fail to be integrable to the left of `s₀`. -/
theorem physicalTime_le_physicalTime
    (hint : ∀ b : ℝ, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b)
    (t₀ : ℝ) {a b : ℝ} (hab : a ≤ b) (hb : s₀ ≤ b) :
    physicalTime t₀ s₀ L a ≤ physicalTime t₀ s₀ L b := by
  rcases le_or_gt s₀ a with ha | ha
  · exact physicalTime_monotoneOn hint t₀ ha hb hab
  · exact (physicalTime_le_of_le_base t₀ ha.le).trans (le_physicalTime_of_le t₀ hb)

/-- **Strict monotonicity where `L > 0`.**  Here the positivity of `L` is
genuinely used: it forces the rate to be nonzero on the whole of `(a,b]`, a
set of positive Lebesgue measure, hence `∫_a^b L² > 0`. -/
theorem physicalTime_strictMonoOn
    (hint : ∀ b : ℝ, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b)
    (hpos : ∀ σ, 0 < L σ) (t₀ : ℝ) :
    StrictMonoOn (physicalTime t₀ s₀ L) (Ici s₀) := by
  intro a ha b hb hab
  have hii : IntervalIntegrable (scaleRate L) volume a b :=
    (hint a ha).symm.trans (hint b hb)
  have hIoc : IntegrableOn (scaleRate L) (Ioc a b) volume :=
    (intervalIntegrable_iff_integrableOn_Ioc_of_le hab.le).mp hii
  have hsupp : Function.support (scaleRate L) ∩ Ioc a b = Ioc a b := by
    rw [inter_eq_right]
    exact fun x _ => Function.mem_support.mpr (scaleRate_pos hpos x).ne'
  have hlt : 0 < ∫ σ in a..b, scaleRate L σ := by
    rw [intervalIntegral.integral_of_le hab.le,
      setIntegral_pos_iff_support_of_nonneg_ae
        (ae_of_all _ fun u => scaleRate_nonneg L u) hIoc, hsupp, Real.volume_Ioc]
    exact ENNReal.ofReal_pos.mpr (by linarith)
  have hadd : (∫ σ in s₀..a, scaleRate L σ) + ∫ σ in a..b, scaleRate L σ
      = ∫ σ in s₀..b, scaleRate L σ :=
    intervalIntegral.integral_add_adjacent_intervals (hint a ha) hii
  simp only [physicalTime]
  linarith

/-! ### The main theorem: infinite rescaled time, finite physical time -/

/-- **F-2, convergence.**  If `L²` is integrable on `[s₀,∞)` then the physical
time converges, as rescaled time `s → ∞`, to the finite value
`blowupTime t₀ s₀ L = t₀ + ∫_{(s₀,∞)} L²`.

This is where the analytic content enters, through mathlib's improper
integral theorem `MeasureTheory.intervalIntegral_tendsto_integral_Ioi`. -/
theorem tendsto_physicalTime (hint : IntegrableOn (scaleRate L) (Ici s₀))
    (t₀ : ℝ) :
    Tendsto (physicalTime t₀ s₀ L) atTop (𝓝 (blowupTime t₀ s₀ L)) := by
  have hIoi : IntegrableOn (scaleRate L) (Ioi s₀) := hint.mono_set Ioi_subset_Ici_self
  have h := MeasureTheory.intervalIntegral_tendsto_integral_Ioi
    (b := fun s : ℝ => s) (l := atTop) s₀ hIoi tendsto_id
  exact h.const_add t₀

/-- **F-2, the bound.**  Physical time never reaches the limit value from
above: `t(s) ≤ T` for every rescaled time `s`. -/
theorem physicalTime_le_blowupTime (hint : IntegrableOn (scaleRate L) (Ici s₀))
    (t₀ : ℝ) (s : ℝ) : physicalTime t₀ s₀ L s ≤ blowupTime t₀ s₀ L := by
  refine ge_of_tendsto (tendsto_physicalTime hint t₀) ?_
  filter_upwards [eventually_ge_atTop s, eventually_ge_atTop s₀] with x hx hx₀
  exact physicalTime_le_physicalTime
    (fun b hb => intervalIntegrable_scaleRate hint hb) t₀ hx hx₀

/-- **F-2, the strict bound.**  If moreover `L > 0` everywhere, the physical
time stays *strictly* below the limit: `t(s) < T` for every `s`.  So the
rescaled trajectory occupies exactly the physical time interval
`[t₀, T)` with `T < ∞`. -/
theorem physicalTime_lt_blowupTime (hint : IntegrableOn (scaleRate L) (Ici s₀))
    (hpos : ∀ σ, 0 < L σ) (t₀ : ℝ) (s : ℝ) :
    physicalTime t₀ s₀ L s < blowupTime t₀ s₀ L := by
  have hloc : ∀ b : ℝ, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b :=
    fun b hb => intervalIntegrable_scaleRate hint hb
  set a := max s s₀ with ha
  have ha₀ : s₀ ≤ a := le_max_right _ _
  have h1 : physicalTime t₀ s₀ L s ≤ physicalTime t₀ s₀ L a :=
    physicalTime_le_physicalTime hloc t₀ (le_max_left _ _) ha₀
  have h2 : physicalTime t₀ s₀ L a < physicalTime t₀ s₀ L (a + 1) :=
    physicalTime_strictMonoOn hloc hpos t₀ ha₀ (by simp; linarith) (by linarith)
  have h3 : physicalTime t₀ s₀ L (a + 1) ≤ blowupTime t₀ s₀ L :=
    physicalTime_le_blowupTime hint t₀ (a + 1)
  linarith

/-- **F-2, packaged.**  Integrability of `L²` on `[s₀,∞)` means: infinite
rescaled time corresponds to a finite physical time interval.  There is a
finite `T` which the physical time approaches but never attains, and the
physical time is monotone on `[s₀,∞)`.

The witness is `T = blowupTime t₀ s₀ L = t₀ + ∫_{(s₀,∞)} L²`. -/
theorem exists_finite_blowupTime (hint : IntegrableOn (scaleRate L) (Ici s₀))
    (t₀ : ℝ) :
    ∃ T : ℝ, Tendsto (physicalTime t₀ s₀ L) atTop (𝓝 T) ∧
      (∀ s : ℝ, physicalTime t₀ s₀ L s ≤ T) ∧
      MonotoneOn (physicalTime t₀ s₀ L) (Ici s₀) :=
  ⟨blowupTime t₀ s₀ L, tendsto_physicalTime hint t₀,
    physicalTime_le_blowupTime hint t₀,
    physicalTime_monotoneOn (fun _ hb => intervalIntegrable_scaleRate hint hb) t₀⟩

/-! ### The converse: non-integrable rate ⟹ infinite physical time -/

/-- **F-2, converse.**  If the rate is locally integrable above `s₀` but *not*
integrable on `(s₀,∞)`, then physical time runs to `+∞` along the rescaled
trajectory.  Together with `tendsto_physicalTime` this shows the F-2
hypothesis is exactly the dividing line between a finite and an infinite
physical blow-up time — it is not merely a convenient sufficient condition.

The local hypothesis cannot be dropped: without it the interval integrals
would be junk-valued and nothing could be said. -/
theorem tendsto_physicalTime_atTop
    (hloc : ∀ b : ℝ, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b)
    (hnot : ¬ IntegrableOn (scaleRate L) (Ioi s₀)) (t₀ : ℝ) :
    Tendsto (physicalTime t₀ s₀ L) atTop atTop := by
  rw [tendsto_atTop]
  intro C
  obtain ⟨s₁, hs₁, hCs₁⟩ : ∃ s₁, s₀ ≤ s₁ ∧ C ≤ physicalTime t₀ s₀ L s₁ := by
    by_contra hcon
    push Not at hcon
    refine hnot ?_
    refine integrableOn_Ioi_of_intervalIntegral_norm_bounded (C - t₀) s₀
      (b := fun s : ℝ => s) (l := atTop) (fun i => ?_) tendsto_id ?_
    · rcases le_or_gt s₀ i with hi | hi
      · exact (intervalIntegrable_iff_integrableOn_Ioc_of_le hi).mp (hloc i hi)
      · rw [Ioc_eq_empty (by simp [hi.le])]
        exact integrableOn_empty
    · filter_upwards [eventually_ge_atTop s₀] with i hi
      have hnorm : ∀ x : ℝ, ‖scaleRate L x‖ = scaleRate L x :=
        fun x => Real.norm_of_nonneg (scaleRate_nonneg L x)
      simp only [hnorm]
      have := hcon i hi
      simp only [physicalTime] at this
      linarith
  filter_upwards [eventually_ge_atTop s₁, eventually_ge_atTop s₀] with x hx hx₀
  exact hCs₁.trans (physicalTime_le_physicalTime hloc t₀ hx hx₀)

end

end NSSingularity
