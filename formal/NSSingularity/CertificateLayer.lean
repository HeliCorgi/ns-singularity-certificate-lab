/-
Stage 1, items **F-17**, **F-18**, **F-19** and the packaged Clay restriction of
`docs/final_target.md` §4.

This file is the Lean side of the certificate chain that
`src/ns_certificate_lab/tail_propagation.py` implements:

  eps_psi, eps_dpsi, eps_d2psi  ->  eps_u  ->  eps_advection  ->  state error.

Each arrow is a finite inequality between nonnegative reals.  That is a
deliberate choice of scope: formalizing the analytic Green integral in one go is
not attempted, and would not be the useful thing to attempt.  What a numerical
certificate can supply is *finitely many nonnegative upper bounds*, and what a
proof needs from it is that those bounds compose correctly.  This file is
exactly that composition layer.

* **F-17** — velocity error from potential error.  The recovery
  `u^r = -r psi_z`, `u^z = 2 psi + r psi_r` is linear, so the bound needs no
  product rule and loses nothing.
* **F-18** — advection error from velocity error and state gradient bounds,
  through the product-difference identity `ab - a'b' = (a-a')b + a'(b-b')`.
* **F-19** — the short-time Gronwall step, obtained from mathlib's
  `norm_le_gronwallBound_of_norm_deriv_right_le` and simplified to the form a
  certificate can check without dividing by a possibly tiny rate.

Scope caveats.
* Every statement is about *given* real numbers.  Nothing here computes a bound,
  and nothing here connects to a PDE solution: supplying the inputs is the job
  of the interval certificate, and that remains a separate obligation.
* The `L^infinity` maximum principle that justifies the differential inequality
  fed to F-19 is **not** formalized; F-19 takes the differential inequality as a
  hypothesis.
* The Clay restriction at the end is still a statement about a coefficient
  trajectory, not about `ClayStatement.ClayPeriodicBreakdown`.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib
import NSSingularity.GalerkinNoBlowup
import NSSingularity.FiniteModeNoGo
import NSSingularity.GreenAndCascade

namespace NSSingularity

noncomputable section

open Set Filter Topology

/-! ## F-17 — potential error to velocity error -/

section VelocityError

variable {R eps0 eps1 : ℝ}

/-- **F-17, radial.**  `|delta u^r| = |r| |d_z delta psi| <= R eps1`. -/
theorem velocity_radial_error_le
    (hR : 0 ≤ R) (hr : |r| ≤ R) (hz : |dz| ≤ eps1) :
    |(-r) * dz| ≤ R * eps1 := by
  rw [abs_mul, abs_neg]
  exact mul_le_mul hr hz (abs_nonneg _) hR

/-- **F-17, axial.**  `|delta u^z| = |2 delta psi + r d_r delta psi|
<= 2 eps0 + R eps1`. -/
theorem velocity_axial_error_le
    (hR : 0 ≤ R) (hr : |r| ≤ R) (hpsi : |p| ≤ eps0) (hdr : |dr| ≤ eps1) :
    |2 * p + r * dr| ≤ 2 * eps0 + R * eps1 := by
  refine (abs_add_le _ _).trans (add_le_add ?_ ?_)
  · rw [abs_mul]
    simpa using mul_le_mul_of_nonneg_left hpsi (by norm_num : (0:ℝ) ≤ |(2:ℝ)|)
  · rw [abs_mul]
    exact mul_le_mul hr hdr (abs_nonneg _) hR

end VelocityError

/-! ## F-18 — the product-difference identity and the advection error -/

section AdvectionError

/-- The identity every bound in this file rests on. -/
theorem product_difference (a b a' b' : ℝ) :
    a * b - a' * b' = (a - a') * b + a' * (b - b') := by ring

/-- **F-18, one product.**  `|ab - a'b'| <= eps_a |b| + |a'| eps_b`. -/
theorem product_error_le {a b a' b' epsa epsb : ℝ}
    (ha : |a - a'| ≤ epsa) (hb : |b - b'| ≤ epsb) :
    |a * b - a' * b'| ≤ epsa * |b| + |a'| * epsb := by
  rw [product_difference]
  refine (abs_add_le _ _).trans (add_le_add ?_ ?_)
  · rw [abs_mul]
    exact mul_le_mul_of_nonneg_right ha (abs_nonneg _)
  · rw [abs_mul]
    exact mul_le_mul_of_nonneg_left hb (abs_nonneg _)

/-- **F-18, the advection term.**  With `A = u^r f_r + u^z f_z`,

`|delta A| <= eps_ur G_r + U_r eps_fr + eps_uz G_z + U_z eps_fz`. -/
theorem advection_error_le
    {ur uz fr fz ur' uz' fr' fz' epsur epsuz epsfr epsfz Gr Gz Ur Uz : ℝ}
    (hur : |ur - ur'| ≤ epsur) (huz : |uz - uz'| ≤ epsuz)
    (hfr : |fr - fr'| ≤ epsfr) (hfz : |fz - fz'| ≤ epsfz)
    (hGr : |fr| ≤ Gr) (hGz : |fz| ≤ Gz)
    (hUr : |ur'| ≤ Ur) (hUz : |uz'| ≤ Uz)
    (hepsur : 0 ≤ epsur) (hepsuz : 0 ≤ epsuz)
    (hepsfr : 0 ≤ epsfr) (hepsfz : 0 ≤ epsfz) :
    |(ur * fr + uz * fz) - (ur' * fr' + uz' * fz')|
      ≤ (epsur * Gr + Ur * epsfr) + (epsuz * Gz + Uz * epsfz) := by
  have hsplit :
      (ur * fr + uz * fz) - (ur' * fr' + uz' * fz')
        = (ur * fr - ur' * fr') + (uz * fz - uz' * fz') := by ring
  rw [hsplit]
  refine (abs_add_le _ _).trans (add_le_add ?_ ?_)
  · refine (product_error_le hur hfr).trans (add_le_add ?_ ?_)
    · exact mul_le_mul_of_nonneg_left hGr hepsur
    · exact mul_le_mul_of_nonneg_right hUr hepsfr
  · refine (product_error_le huz hfz).trans (add_le_add ?_ ?_)
    · exact mul_le_mul_of_nonneg_left hGz hepsuz
    · exact mul_le_mul_of_nonneg_right hUz hepsfz

end AdvectionError

/-! ## F-19 — the short-time Gronwall step -/

section Gronwall

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- `gronwallBound` is dominated by the certificate-friendly form
`(delta + eps t) exp(K t)` for nonnegative data.  The proof is
`exp(x) - 1 <= x exp(x)` for `x >= 0`.  Nonnegativity of `delta` turns out to be
unnecessary and is kept only because every caller has it. -/
theorem gronwallBound_le_simple {δ K ε t : ℝ}
    (_hδ : 0 ≤ δ) (hK : 0 ≤ K) (hε : 0 ≤ ε) (ht : 0 ≤ t) :
    gronwallBound δ K ε t ≤ (δ + ε * t) * Real.exp (K * t) := by
  rcases eq_or_lt_of_le hK with hK0 | hKpos
  · simp [gronwallBound, ← hK0, Real.exp_zero]
  · have hne : K ≠ 0 := ne_of_gt hKpos
    rw [gronwallBound_of_K_ne_0 hne]
    have hx : 0 ≤ K * t := mul_nonneg hK ht
    have hpos : 0 < Real.exp (K * t) := Real.exp_pos _
    have hkey : Real.exp (K * t) - 1 ≤ (K * t) * Real.exp (K * t) := by
      have h1 : (1 : ℝ) - K * t ≤ Real.exp (-(K * t)) := by
        have := Real.add_one_le_exp (-(K * t))
        linarith
      have h2 := mul_le_mul_of_nonneg_right h1 hpos.le
      rw [Real.exp_neg, inv_mul_cancel₀ (ne_of_gt hpos)] at h2
      nlinarith [h2]
    have hεK : 0 ≤ ε / K := div_nonneg hε hKpos.le
    have hstep : (ε / K) * (Real.exp (K * t) - 1) ≤ ε * t * Real.exp (K * t) := by
      have hmul := mul_le_mul_of_nonneg_left hkey hεK
      have hid : (ε / K) * ((K * t) * Real.exp (K * t)) = ε * t * Real.exp (K * t) := by
        field_simp
      linarith [hmul, hid.le, hid.ge]
    nlinarith [hpos, hstep]

/-- **F-19.**  A curve whose derivative satisfies `‖f'‖ ≤ K‖f‖ + eps` and whose
initial norm is at most `delta` obeys `‖f t‖ ≤ (delta + eps t) exp(K t)`.

This is the form the tail-propagation certificate checks: it never divides by
`K`, so a certificate with a tiny Lipschitz constant stays checkable. -/
theorem norm_le_simple_gronwall {f f' : ℝ → E} {δ K ε a b : ℝ}
    (_hab : a ≤ b) (hδ : 0 ≤ δ) (hK : 0 ≤ K) (hε : 0 ≤ ε)
    (hf : ContinuousOn f (Icc a b))
    (hf' : ∀ x ∈ Ico a b, HasDerivWithinAt f (f' x) (Ici x) x)
    (ha : ‖f a‖ ≤ δ)
    (bound : ∀ x ∈ Ico a b, ‖f' x‖ ≤ K * ‖f x‖ + ε) :
    ∀ x ∈ Icc a b, ‖f x‖ ≤ (δ + ε * (x - a)) * Real.exp (K * (x - a)) := by
  intro x hx
  have hmain := norm_le_gronwallBound_of_norm_deriv_right_le hf hf' ha bound x hx
  exact hmain.trans (gronwallBound_le_simple hδ hK hε (by linarith [hx.1]))

end Gronwall

/-! ## The packaged Clay restriction

`breakdown_time_set_empty` refutes a breakdown for a trajectory supplied with
loose hypotheses.  Packaging those hypotheses into a structure makes the class
the theorem excludes explicit, and makes it impossible to drop one silently.
-/

section ClayRestriction

variable (E : Type*) [NormedAddCommGroup E] [InnerProductSpace ℝ E] [CompleteSpace E]

open RealInnerProductSpace

/-- A **fixed finite-bandwidth Track-F candidate**, with every hypothesis named.

This is the class the no-go excludes.  It is *not* all of Track F: candidates
whose Fourier bandwidth diverges as `t → T` are outside it, and their necessary
conditions live in `ShellAdmissible` (F-16). -/
structure FixedBandwidthCandidate where
  /-- The coefficient trajectory of the ansatz. -/
  trajectory : ℝ → E
  /-- The projected advection term. -/
  advection : E →L[ℝ] E →L[ℝ] E
  /-- The projected viscous term. -/
  viscous : E →L[ℝ] E
  /-- The projected force. -/
  force : ℝ → E
  /-- A uniform bound on the projected force. -/
  forceBound : ℝ
  energy_neutral : EnergyNeutral fun x y => advection x y
  dissipative : Dissipative fun x => viscous x
  force_bounded : ∀ t, ‖force t‖ ≤ forceBound
  force_continuous : ContinuousOn force (Ici (0 : ℝ))
  trajectory_continuous : ContinuousOn trajectory (Ici (0 : ℝ))
  evolution : ∀ t ∈ Ici (0 : ℝ),
    HasDerivAt trajectory
      (force t + advection (trajectory t) (trajectory t) + viscous (trajectory t)) t

variable {E}

/-- **The packaged no-go.**  A fixed finite-bandwidth candidate never witnesses a
breakdown, at any time.

What this does not say, and what `formal/README.md` records: it is a statement
about the coefficient trajectory, and the bridge to
`ClayStatement.ClayPeriodicBreakdown` still needs the Fourier isomorphism, the
analytic half of the trilinear identity, and Navier–Stokes local uniqueness. -/
theorem FixedBandwidthCandidate.breakdown_times_empty
    (candidate : FixedBandwidthCandidate E) :
    {T : ℝ | IsBreakdownCandidate candidate.trajectory T} = ∅ :=
  breakdown_time_set_empty candidate.energy_neutral candidate.dissipative
    candidate.force_bounded candidate.force_continuous
    candidate.trajectory_continuous candidate.evolution

/-- The same conclusion in the form a reader of Clay (D) wants: at no positive
time does the candidate fail to reach its state. -/
theorem FixedBandwidthCandidate.reaches_every_time
    (candidate : FixedBandwidthCandidate E) {T : ℝ} (hT : 0 < T) :
    ∃ L : E, Tendsto candidate.trajectory (𝓝[<] T) (𝓝 L) := by
  by_contra hno
  have hmem : T ∈ {T : ℝ | IsBreakdownCandidate candidate.trajectory T} := ⟨hT, hno⟩
  rw [candidate.breakdown_times_empty] at hmem
  exact hmem

end ClayRestriction

end

end NSSingularity
