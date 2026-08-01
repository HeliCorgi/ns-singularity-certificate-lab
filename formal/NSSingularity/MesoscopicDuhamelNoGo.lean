/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import Mathlib
import NSSingularity.FiniteModeNoGo

/-!
# Phase-independent mesoscopic Duhamel no-go: finite algebraic core

This file formalizes only the nonnegative-real algebra behind the estimate

`(V / U)^2 <= 2 * kappa^2 * tau^2 * cE * Meff / N^3`.

The hypotheses deliberately expose every analytic input: the energy
normalization, the bilinear bound, and the short-time Duhamel bound.  No
Fourier transform, Leray projection, Navier--Stokes solution, empty-child
geometry, or PDE-to-finite-algebra bridge is defined or proved here.

The final finite-support lemma reuses F-13's Cauchy--Schwarz inequality to
show that the usual coefficient effective count is at most the support
cardinality.  It is still a finite-sum statement, not a PDE estimate.
-/

namespace NSSingularity

noncomputable section

/-! ## MD1 -- abstract phase-independent ratio estimate -/

/-- **MD1.**  The finite algebraic core of the phase-independent
empty-child estimate.

`U` is the parent amplitude scale, `B` a bilinear source bound, and `V` the
resulting short-time Duhamel scale.  The theorem assumes, rather than proves,
the three application-specific inputs

* `U^2 = 2*cE/N`,
* `B <= kappa*N*sqrt(Meff)*U^2`, and
* `V <= (tau/N^2)*B`.

Thus this theorem has no conclusion about a Navier--Stokes solution unless a
separate Fourier/PDE bridge proves those hypotheses. -/
theorem emptyChild_duhamel_ratio_sq_le
    {U B V N kappa tau cE Meff : ℝ}
    (hU : 0 < U) (hN : 0 < N)
    (_hB : 0 ≤ B) (hV : 0 ≤ V)
    (hkappa : 0 ≤ kappa) (htau : 0 ≤ tau)
    (_hcE : 0 ≤ cE) (hMeff : 0 ≤ Meff)
    (henergy : U ^ 2 = 2 * cE / N)
    (hbilinear : B ≤ kappa * N * Real.sqrt Meff * U ^ 2)
    (hduhamel : V ≤ tau / N ^ 2 * B) :
    (V / U) ^ 2 ≤ 2 * kappa ^ 2 * tau ^ 2 * cE * Meff / N ^ 3 := by
  have htime : 0 ≤ tau / N ^ 2 := div_nonneg htau (sq_nonneg N)
  have hsqrt : 0 ≤ Real.sqrt Meff := Real.sqrt_nonneg Meff
  have hsource : 0 ≤ kappa * N * Real.sqrt Meff * U ^ 2 := by positivity
  have hVB :
      V ≤ (tau / N ^ 2) * (kappa * N * Real.sqrt Meff * U ^ 2) :=
    hduhamel.trans (mul_le_mul_of_nonneg_left hbilinear htime)
  have hR :
      0 ≤ (tau / N ^ 2) * (kappa * N * Real.sqrt Meff * U ^ 2) :=
    mul_nonneg htime hsource
  have hsq :
      V ^ 2 ≤ ((tau / N ^ 2) *
        (kappa * N * Real.sqrt Meff * U ^ 2)) ^ 2 := by
    exact (sq_le_sq₀ hV hR).2 hVB
  calc
    (V / U) ^ 2 = V ^ 2 / U ^ 2 := by ring
    _ ≤ ((tau / N ^ 2) *
        (kappa * N * Real.sqrt Meff * U ^ 2)) ^ 2 / U ^ 2 :=
      div_le_div_of_nonneg_right hsq (sq_nonneg U)
    _ = 2 * kappa ^ 2 * tau ^ 2 * cE * Meff / N ^ 3 := by
      field_simp [ne_of_gt hU, ne_of_gt hN]
      rw [Real.sq_sqrt hMeff]
      have henergyN : U ^ 2 * N = 2 * cE :=
        (eq_div_iff (ne_of_gt hN)).mp henergy
      calc
        tau ^ 2 * N * kappa ^ 2 * Meff * U ^ 2 =
            tau ^ 2 * kappa ^ 2 * Meff * (U ^ 2 * N) := by ring
        _ = tau ^ 2 * kappa ^ 2 * Meff * (2 * cE) := by rw [henergyN]
        _ = tau ^ 2 * kappa ^ 2 * 2 * cE * Meff := by ring

/-! ## MD2 -- replacing the effective count by any certified upper bound -/

/-- If `Meff <= M`, MD1 immediately yields the weaker support-scale bound
with `M`.  This isolates the only monotonicity step needed after a separate
effective-count estimate. -/
theorem emptyChild_duhamel_ratio_sq_le_of_effectiveCount_le
    {U B V N kappa tau cE Meff M : ℝ}
    (hU : 0 < U) (hN : 0 < N)
    (hB : 0 ≤ B) (hV : 0 ≤ V)
    (hkappa : 0 ≤ kappa) (htau : 0 ≤ tau)
    (hcE : 0 ≤ cE) (hMeff : 0 ≤ Meff)
    (hMeffM : Meff ≤ M)
    (henergy : U ^ 2 = 2 * cE / N)
    (hbilinear : B ≤ kappa * N * Real.sqrt Meff * U ^ 2)
    (hduhamel : V ≤ tau / N ^ 2 * B) :
    (V / U) ^ 2 ≤ 2 * kappa ^ 2 * tau ^ 2 * cE * M / N ^ 3 := by
  calc
    (V / U) ^ 2 ≤ 2 * kappa ^ 2 * tau ^ 2 * cE * Meff / N ^ 3 :=
      emptyChild_duhamel_ratio_sq_le hU hN hB hV hkappa htau hcE hMeff
        henergy hbilinear hduhamel
    _ ≤ 2 * kappa ^ 2 * tau ^ 2 * cE * M / N ^ 3 := by
      apply div_le_div_of_nonneg_right _ (by positivity)
      exact mul_le_mul_of_nonneg_left hMeffM (by positivity)

/-! ## MD3 -- finite coefficient effective count -/

/-- The coefficient effective count `(sum |c_i|)^2 / sum c_i^2` on a finite
support.  Division is total in Lean; the useful comparison theorem below
assumes positive squared mass, so its denominator is genuinely nonzero. -/
def finiteEffectiveModeCount {ι : Type*} (s : Finset ι) (c : ι → ℝ) : ℝ :=
  (∑ i ∈ s, |c i|) ^ 2 / ∑ i ∈ s, c i ^ 2

/-- A positive-energy finite coefficient family has nonnegative effective
mode count. -/
theorem finiteEffectiveModeCount_nonneg {ι : Type*} (s : Finset ι)
    (c : ι → ℝ) (henergy : 0 < ∑ i ∈ s, c i ^ 2) :
    0 ≤ finiteEffectiveModeCount s c := by
  exact div_nonneg (sq_nonneg _) henergy.le

/-- **MD3.**  F-13's finite Cauchy--Schwarz inequality gives
`Meff <= |support|`.  No infinite Fourier series or PDE norm is involved. -/
theorem finiteEffectiveModeCount_le_card {ι : Type*} (s : Finset ι)
    (c : ι → ℝ) (henergy : 0 < ∑ i ∈ s, c i ^ 2) :
    finiteEffectiveModeCount s c ≤ (s.card : ℝ) := by
  rw [finiteEffectiveModeCount, div_le_iff₀ henergy]
  exact sq_sum_abs_le_card_mul_sum_sq s c

/-- Combining MD1 with finite Cauchy--Schwarz replaces the effective count by
the coefficient support cardinality.  The bilinear hypothesis is still an
input; this theorem does not derive it from Navier--Stokes. -/
theorem emptyChild_duhamel_ratio_sq_le_card
    {ι : Type*} (s : Finset ι) (c : ι → ℝ)
    {U B V N kappa tau cE : ℝ}
    (hU : 0 < U) (hN : 0 < N)
    (hB : 0 ≤ B) (hV : 0 ≤ V)
    (hkappa : 0 ≤ kappa) (htau : 0 ≤ tau)
    (hcE : 0 ≤ cE)
    (hcoeffEnergy : 0 < ∑ i ∈ s, c i ^ 2)
    (henergy : U ^ 2 = 2 * cE / N)
    (hbilinear : B ≤ kappa * N *
      Real.sqrt (finiteEffectiveModeCount s c) * U ^ 2)
    (hduhamel : V ≤ tau / N ^ 2 * B) :
    (V / U) ^ 2 ≤
      2 * kappa ^ 2 * tau ^ 2 * cE * (s.card : ℝ) / N ^ 3 := by
  exact emptyChild_duhamel_ratio_sq_le_of_effectiveCount_le
    hU hN hB hV hkappa htau hcE
    (finiteEffectiveModeCount_nonneg s c hcoeffEnergy)
    (finiteEffectiveModeCount_le_card s c hcoeffEnergy)
    henergy hbilinear hduhamel

end

end NSSingularity
