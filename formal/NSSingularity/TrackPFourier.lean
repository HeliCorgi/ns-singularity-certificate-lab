/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
/-
Track P — the periodic lane — algebraic and finite-dimensional core.

## The mathematics this file sits under

Track P controls the distance between the true solution of Navier–Stokes on the
torus `T^3` (fixed viscosity `nu > 0`, mean-zero divergence-free data,
normalised measure) and a finite Fourier Galerkin trajectory `u_a`.  The error
`w = u - u_a` obeys

  `dw/dt = nu Δw - P(w·∇u_a) - P(u_a·∇w) - P(w·∇w) - e`,

where `e = (I - P_G) P(u_a·∇u_a)` is *exactly* the Galerkin tail.  Because
`u_a` is a finite trigonometric polynomial, its continuous PDE residual is
again a finite trigonometric polynomial, supported on the modes of `(G+G) \ G`
and computable exactly in rational arithmetic on Fourier coefficients.  That is
the entire reason the periodic lane closes where the `R^3` lane of
`ControlODE.lean` does not: no spatial interpolation and no
discrete-to-continuous recovery ever enters, and the Poincaré gap `|k| ≥ 1` of
the mean-zero torus supplies the linear damping that the whole space lacks.
The `H^4` energy estimate — with the crude `ℓ^1` embedding
`‖f‖_∞ ≤ Adot ‖f‖_2`, which is the reason `n = 4` and not `n = 3`: at `n = 3`
one commutator term needs a sharp Gagliardo–Nirenberg constant, at `n = 4`
every term closes with `Adot` alone — yields the Riccati control inequality

  `d⁺R/dt ≤ (-nu + 9(K1 + K2)) R + 135 Adot R² + ‖e‖_4`,  `R(0) = 0`,

whose constants `K1`, `K2` (sup norms of derivatives of `u_a`) and `Adot` (a
lattice sum) are finite Fourier data.  Conditional on the named external
classical theorems EXT-P1 (Kato/Temam local `H^4` existence and uniqueness),
EXT-P2 (regularity sufficient to run the estimate and the Dini comparison) and
EXT-P3 (the `H^4` continuation criterion) — recorded faithfully in the research
notes and **never as Lean axioms** — the comparison of `ControlODE.lean` bounds
`‖u - u_a‖_{Ḣ^4}` by the majorant on the whole slab.  This is an a posteriori
*distance* bound for a surrogate trajectory.  It is not a singularity result,
and nothing in this file, alone or in combination with the rest of the
library, is a statement about the Clay problem.

## The fixed-band / finite-band distinction

The finite-mode no-go (`FiniteModeNoGo.lean`, packaged in
`CertificateLayer.lean`) excludes trajectories that *remain* in a fixed finite
Fourier band for all time.  A finite-band *initial datum* does not satisfy
that hypothesis: the Galerkin tail `e` is generically nonzero, so the true
solution leaves the band immediately, and nothing about it is excluded.
Track P works precisely in the unexcluded region — finite-band datum,
infinite-band true solution, distance controlled.  `FixedBandTrajectory` and
`FiniteBandDatum` state the two hypotheses so that they cannot be conflated:
the implication one way is a one-liner, and the converse is refuted by an
explicit counterexample (`u t = (1, t)` against the first-coordinate axis).
`FixedBandwidthCandidate.fixedBand_scope` then restates the scope of the
existing no-go through the new predicate: the excluded class carries the
fixed-band hypothesis, and in particular has a finite-band datum — but the
counterexample shows the converse inclusion fails, which is why the no-go says
nothing about Track P's objects.

## What is formalised

* `leray` and its five properties — the single-mode Leray multiplier
  `a ↦ a - (⟪k,a⟫/⟪k,k⟫) k`, which is what the Leray projector `P` *is* on the
  Fourier coefficient of the mode `k`: annihilation by `k`
  (`inner_leray_eq_zero`), fixing of transverse amplitudes
  (`leray_eq_self_of_inner_eq_zero`), idempotence (`leray_leray`),
  self-adjointness (`inner_leray_left_eq_inner_leray_right`) and the
  contraction property (`norm_leray_le`, from Pythagoras).  All five hold with
  no hypothesis on `k`: at `k = 0` the multiplier degenerates to the identity
  because Lean's `0/0 = 0`, and every statement remains true.
* `slotDivergence_cosMode`, `slotDivergence_sinMode` — a single Fourier mode
  `x ↦ cos(k·x) a` (or `sin`) with transverse amplitude `k·a = 0` has
  pointwise divergence zero, stated as a genuine divergence theorem on
  `Fin 3 → ℝ` through slot derivatives, in the same style as the cylindrical
  divergence of `L3Generation.lean`.  The collapse is structural: each slot
  derivative produces `-sin(k·x) kᵢaᵢ`, and the sum dies by transversality.
* `contDiff_trigPolynomial` — a finite sum of cosine and sine modes is
  `C^∞`.  This is the certificate that finite Fourier data are smooth
  classical fields, so the Galerkin trajectory is a *bona fide* strong
  solution of its own (projected) equation and its PDE residual is a classical
  object, not a distribution.
* `FixedBandTrajectory`, `FiniteBandDatum`,
  `FixedBandTrajectory.finiteBandDatum`,
  `exists_finiteBandDatum_not_fixedBandTrajectory`,
  `FixedBandwidthCandidate.fixedBand_scope` — the distinction of the
  paragraph above, stated so it cannot be misread.
* `weighted_sum_succ_mono`, `weighted_sum_pow_mono` — the homogeneous-norm
  ladder: with weights `w ≥ 1` (on the mean-zero torus, `w = |k|² ≥ 1`) the
  weighted coefficient sums are monotone in the exponent, which is
  `‖f‖_m ≤ ‖f‖_n` for `m ≤ n` in the form the control inequality uses.
* `trackP_slab_error_le` — the composition shell: the Track-P control
  inequality is exactly a Riccati subsolution hypothesis with datum error
  zero (the Galerkin datum *is* the datum), so `ControlODE.riccati_comparison`
  bounds the error by any supersolution on the slab.  This is the point where
  the two layers — the estimate layer above, the ODE layer of
  `ControlODE.lean` — are shown to compose.

## What is not formalised, and why

* **The `H^4` energy estimate itself.**  Passing from the error equation to
  the Riccati inequality needs Sobolev spaces on `T^3`, the Fourier isometry,
  and integration by parts on the torus; mathlib has none of these in usable
  form, and taking the estimate as a hypothesis is exactly what
  `trackP_slab_error_le` does.  The paper derivation (Leibniz commutators,
  the multinomial weight `Σ c_α = 81`, the `9 = √81` Cauchy–Schwarz, the
  `135 = 9·15` advection constant) lives in the research notes; per the task
  instruction, the identity `Σ_{|α|=4} 4!/α! = 81` is deliberately not
  restated in Lean because no Lean statement here consumes it.
* **The Galerkin tail identity** `e = (I - P_G) P(u_a·∇u_a)` as an operator
  statement — same missing torus function spaces.  Its finite *coefficient*
  ingredients are here: the Leray multiplier (F1) is the mode-wise `P`, and
  the trilinear resonance algebra is already `advectionForm_eq_zero` in
  `FiniteModeNoGo.lean`.
* **The `Adot` lattice-sum enclosure** (exact sum over `|k|_∞ ≤ N` plus the
  `26/N` tail).  It is finite rational arithmetic plus a comparison of a tail
  with `Σ m⁻²`, checked by the Python layer; formalising it is a computation,
  not an idea, and it resisted honest effort within this session's budget
  without adding anything the composition shell needs.
* **EXT-P1/2/3** (Kato/Temam local existence, regularity, `H^4`
  continuation).  External classical theorems, recorded as named conditions
  in the research notes.  Encoding them as Lean axioms is forbidden in this
  repository, and rightly so; they therefore appear nowhere in this file, and
  every theorem below is unconditionally true as stated.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib
import NSSingularity.CertificateLayer
import NSSingularity.ControlODE

namespace NSSingularity

noncomputable section

open Set Filter Topology
open scoped RealInnerProductSpace

/-! ## F1 — the Leray multiplier, finite algebra

On the mode `k` the Leray projector acts on the Fourier amplitude `a` as
`a ↦ a - (⟪k,a⟫/⟪k,k⟫) k`, the orthogonal projection onto the plane
transverse to `k`.  Everything the error estimate needs from `P` — that it
fixes divergence-free data, kills gradients, is self-adjoint, and never
increases a norm — is this finite-dimensional algebra, mode by mode.  The
statements are given in an arbitrary real inner product space because the
dimension never enters. -/

section LerayMultiplier

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The Leray multiplier on the mode `k`: orthogonal projection of the
amplitude `a` onto the complement of `k`.  At `k = 0` the coefficient is
`0/0 = 0` and the multiplier is the identity, which is the right degeneration:
the mean mode has no divergence constraint. -/
def leray (k a : E) : E := a - (⟪k, a⟫ / ⟪k, k⟫) • k

/-- **F1 (i).**  The projected amplitude is transverse to its wavevector: this
is `div (P u) = 0` mode by mode, and it needs no hypothesis on `k`. -/
theorem inner_leray_eq_zero (k a : E) : ⟪k, leray k a⟫ = 0 := by
  rcases eq_or_ne k 0 with rfl | hk
  · simp [leray]
  · have hkk : ⟪k, k⟫ ≠ (0 : ℝ) := inner_self_ne_zero.mpr hk
    simp only [leray, inner_sub_right, real_inner_smul_right]
    rw [div_mul_cancel₀ _ hkk, sub_self]

/-- **F1 (ii).**  A transverse amplitude is fixed: `P` restricted to
divergence-free data is the identity, which is why `⟨P f, w⟩ = ⟨f, w⟩` in the
energy estimate. -/
theorem leray_eq_self_of_inner_eq_zero {k a : E} (h : ⟪k, a⟫ = 0) :
    leray k a = a := by
  simp [leray, h]

/-- **F1 (iii).**  Idempotence, the algebraic signature of a projection.
Immediate from (i) and (ii): the output of `leray k` is transverse, and
transverse amplitudes are fixed. -/
theorem leray_leray (k a : E) : leray k (leray k a) = leray k a :=
  leray_eq_self_of_inner_eq_zero (inner_leray_eq_zero k a)

/-- **F1 (iv).**  Self-adjointness.  This is what lets the estimate move `P`
off the nonlinearity and onto the divergence-free error `w`, where it
disappears by (ii).  Real symmetry of the inner product is the only input. -/
theorem inner_leray_left_eq_inner_leray_right (k a b : E) :
    ⟪leray k a, b⟫ = ⟪a, leray k b⟫ := by
  simp only [leray, inner_sub_left, inner_sub_right, real_inner_smul_left,
    real_inner_smul_right]
  rw [real_inner_comm a k]
  ring

/-- **F1 (v).**  The contraction property, from Pythagoras: the amplitude
splits orthogonally as `leray k a + (⟪k,a⟫/⟪k,k⟫) k`, so dropping the second
summand can only shrink the norm.  This is the reason `P` costs nothing in any
of the estimate's constants. -/
theorem norm_leray_le (k a : E) : ‖leray k a‖ ≤ ‖a‖ := by
  have hperp : ⟪leray k a, (⟪k, a⟫ / ⟪k, k⟫) • k⟫ = 0 := by
    rw [real_inner_smul_right,
      show ⟪leray k a, k⟫ = ⟪k, leray k a⟫ from real_inner_comm _ _,
      inner_leray_eq_zero, mul_zero]
  have hsum : leray k a + (⟪k, a⟫ / ⟪k, k⟫) • k = a := by
    simp only [leray]
    abel
  have hsq : ‖a‖ ^ 2 = ‖leray k a‖ ^ 2 + ‖(⟪k, a⟫ / ⟪k, k⟫) • k‖ ^ 2 := by
    calc ‖a‖ ^ 2 = ‖leray k a + (⟪k, a⟫ / ⟪k, k⟫) • k‖ ^ 2 := by rw [hsum]
      _ = ‖leray k a‖ ^ 2 + 2 * ⟪leray k a, (⟪k, a⟫ / ⟪k, k⟫) • k⟫
          + ‖(⟪k, a⟫ / ⟪k, k⟫) • k‖ ^ 2 := norm_add_sq_real _ _
      _ = ‖leray k a‖ ^ 2 + ‖(⟪k, a⟫ / ⟪k, k⟫) • k‖ ^ 2 := by
          rw [hperp]; ring
  have hle : ‖leray k a‖ ^ 2 ≤ ‖a‖ ^ 2 := by
    have := sq_nonneg ‖(⟪k, a⟫ / ⟪k, k⟫) • k‖
    linarith
  calc ‖leray k a‖ = Real.sqrt (‖leray k a‖ ^ 2) :=
        (Real.sqrt_sq (norm_nonneg _)).symm
    _ ≤ Real.sqrt (‖a‖ ^ 2) := Real.sqrt_le_sqrt hle
    _ = ‖a‖ := Real.sqrt_sq (norm_nonneg _)

end LerayMultiplier

/-! ## F2 — single-mode divergence

A single Fourier mode with transverse amplitude is divergence free as a
classical field on `Fin 3 → ℝ`, not merely as coefficient algebra.  The
divergence is assembled from slot derivatives exactly as the cylindrical
divergence of `L3Generation.lean`; each slot contributes
`∓ sin/cos (k·x) · kᵢ aᵢ` by the one-variable chain rule along an affine
phase, and the sum collapses by `k·a = 0`.  The cancellation is structural —
no regularity of anything beyond the trigonometric functions is used. -/

section SingleModeDivergence

/-- Partial derivative of a scalar field on `Fin 3 → ℝ` in coordinate slot
`i`, through `Function.update`: freeze the other two coordinates and
differentiate the remaining one-variable function. -/
def slotDeriv (f : (Fin 3 → ℝ) → ℝ) (i : Fin 3) (x : Fin 3 → ℝ) : ℝ :=
  deriv (fun t => f (Function.update x i t)) (x i)

/-- The divergence of a vector field on `Fin 3 → ℝ`, as the sum of the three
diagonal slot derivatives. -/
def slotDivergence (V : (Fin 3 → ℝ) → Fin 3 → ℝ) (x : Fin 3 → ℝ) : ℝ :=
  ∑ i, slotDeriv (fun y => V y i) i x

/-- The cosine mode `x ↦ cos(k·x) a` with wavevector `k` and amplitude `a`. -/
def cosMode (k a : Fin 3 → ℝ) : (Fin 3 → ℝ) → Fin 3 → ℝ :=
  fun x i => Real.cos (∑ j, k j * x j) * a i

/-- The sine mode `x ↦ sin(k·x) a`. -/
def sinMode (k a : Fin 3 → ℝ) : (Fin 3 → ℝ) → Fin 3 → ℝ :=
  fun x i => Real.sin (∑ j, k j * x j) * a i

/-- Freezing all slots but `i` makes the phase affine with slope `k i`. -/
private theorem sum_mul_update (k x : Fin 3 → ℝ) (i : Fin 3) (t : ℝ) :
    ∑ j, k j * Function.update x i t j
      = k i * t + ∑ j ∈ Finset.univ.erase i, k j * x j := by
  rw [← Finset.add_sum_erase _ (fun j => k j * Function.update x i t j)
    (Finset.mem_univ i)]
  congr 1
  · rw [Function.update_self]
  · exact Finset.sum_congr rfl fun j hj => by
      rw [Function.update_of_ne (Finset.ne_of_mem_erase hj)]

/-- The one-variable core: the phase seen through slot `i` has derivative
`k i`. -/
private theorem hasDerivAt_slotPhase (k x : Fin 3 → ℝ) (i : Fin 3) :
    HasDerivAt (fun t => ∑ j, k j * Function.update x i t j) (k i) (x i) := by
  have hfun : (fun t => ∑ j, k j * Function.update x i t j)
      = fun t => k i * t + ∑ j ∈ Finset.univ.erase i, k j * x j := by
    funext t
    exact sum_mul_update k x i t
  have haff : HasDerivAt
      (fun t : ℝ => k i * t + ∑ j ∈ Finset.univ.erase i, k j * x j)
      (k i * 1) (x i) :=
    ((hasDerivAt_id (x i)).const_mul (k i)).add_const _
  rw [hfun]
  simpa using haff

/-- At the base point the frozen phase is the phase. -/
private theorem slotPhase_at (k x : Fin 3 → ℝ) (i : Fin 3) :
    ∑ j, k j * Function.update x i (x i) j = ∑ j, k j * x j := by
  rw [Function.update_eq_self]

/-- **F2, cosine.**  A cosine mode with transverse amplitude is divergence
free: each slot derivative is `-sin(k·x) kᵢ aᵢ`, and the sum is
`-sin(k·x) (k·a) = 0`.  This is the pointwise fact that makes finite Fourier
data with Leray-projected amplitudes admissible velocity fields. -/
theorem slotDivergence_cosMode (k a : Fin 3 → ℝ)
    (hka : ∑ i, k i * a i = 0) (x : Fin 3 → ℝ) :
    slotDivergence (cosMode k a) x = 0 := by
  have hterm : ∀ i : Fin 3, slotDeriv (fun y => cosMode k a y i) i x
      = -Real.sin (∑ j, k j * x j) * (k i * a i) := by
    intro i
    have hcos : HasDerivAt
        (fun t => Real.cos (∑ j, k j * Function.update x i t j))
        (-Real.sin (∑ j, k j * x j) * k i) (x i) := by
      have h := (hasDerivAt_slotPhase k x i).cos
      rwa [slotPhase_at] at h
    have hmul := hcos.mul_const (a i)
    have hd : slotDeriv (fun y => cosMode k a y i) i x
        = -Real.sin (∑ j, k j * x j) * k i * a i := by
      simp only [slotDeriv, cosMode]
      exact hmul.deriv
    rw [hd]; ring
  simp only [slotDivergence]
  calc ∑ i, slotDeriv (fun y => cosMode k a y i) i x
      = ∑ i, -Real.sin (∑ j, k j * x j) * (k i * a i) :=
        Finset.sum_congr rfl fun i _ => hterm i
    _ = -Real.sin (∑ j, k j * x j) * ∑ i, k i * a i := by
        rw [Finset.mul_sum]
    _ = 0 := by rw [hka, mul_zero]

/-- **F2, sine.**  Same collapse with `cos(k·x)` in place of `-sin(k·x)`. -/
theorem slotDivergence_sinMode (k a : Fin 3 → ℝ)
    (hka : ∑ i, k i * a i = 0) (x : Fin 3 → ℝ) :
    slotDivergence (sinMode k a) x = 0 := by
  have hterm : ∀ i : Fin 3, slotDeriv (fun y => sinMode k a y i) i x
      = Real.cos (∑ j, k j * x j) * (k i * a i) := by
    intro i
    have hsin : HasDerivAt
        (fun t => Real.sin (∑ j, k j * Function.update x i t j))
        (Real.cos (∑ j, k j * x j) * k i) (x i) := by
      have h := (hasDerivAt_slotPhase k x i).sin
      rwa [slotPhase_at] at h
    have hmul := hsin.mul_const (a i)
    have hd : slotDeriv (fun y => sinMode k a y i) i x
        = Real.cos (∑ j, k j * x j) * k i * a i := by
      simp only [slotDeriv, sinMode]
      exact hmul.deriv
    rw [hd]; ring
  simp only [slotDivergence]
  calc ∑ i, slotDeriv (fun y => sinMode k a y i) i x
      = ∑ i, Real.cos (∑ j, k j * x j) * (k i * a i) :=
        Finset.sum_congr rfl fun i _ => hterm i
    _ = Real.cos (∑ j, k j * x j) * ∑ i, k i * a i := by
        rw [Finset.mul_sum]
    _ = 0 := by rw [hka, mul_zero]

end SingleModeDivergence

/-! ## F3 — finite Fourier data are `C^∞`

The Galerkin trajectory at each instant is a finite sum of cosine and sine
modes, so it is a smooth classical field; consequently every derivative the
`H^4` estimate takes of `u_a` is a classical derivative, and the PDE residual
`e` is a classical field too.  The proof is compositional: the phase is a
continuous linear functional, `cos`/`sin` are `C^∞`, and finite sums preserve
smoothness. -/

section TrigSmoothness

/-- **F3.**  A finite trigonometric polynomial with vector amplitudes is
`C^∞`.  This is the "finite Fourier data are smooth" certificate: the
statement quantifies over an arbitrary finite mode set and arbitrary cosine
and sine amplitudes, which is exactly the class of Track-P surrogate
trajectories. -/
theorem contDiff_trigPolynomial {ι : Type*} (s : Finset ι)
    (k : ι → Fin 3 → ℝ) (a b : ι → Fin 3 → ℝ) :
    ContDiff ℝ (⊤ : ℕ∞) fun x : Fin 3 → ℝ =>
      ∑ m ∈ s, (Real.cos (∑ j, k m j * x j) • a m
        + Real.sin (∑ j, k m j * x j) • b m) := by
  refine ContDiff.sum fun m _ => ?_
  have hproj : ∀ j : Fin 3, ContDiff ℝ (⊤ : ℕ∞) fun x : Fin 3 → ℝ => x j :=
    fun j => (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : Fin 3 => ℝ)
      j).contDiff
  have hphase : ContDiff ℝ (⊤ : ℕ∞) fun x : Fin 3 → ℝ => ∑ j, k m j * x j :=
    ContDiff.sum fun j _ => contDiff_const.mul (hproj j)
  exact (hphase.cos.smul contDiff_const).add (hphase.sin.smul contDiff_const)

end TrigSmoothness

/-! ## F4 — the fixed-band / finite-band distinction

Two predicates that must never be conflated.  The finite-mode no-go consumes
the first; Track P's objects satisfy only the second, and the counterexample
shows the second does not imply the first — so the no-go excludes nothing
about a finite-band datum whose true evolution leaves the band, which is the
generic case because the Galerkin tail is generically nonzero. -/

section BandDistinction

variable {E : Type*} [AddCommGroup E] [Module ℝ E]

/-- The hypothesis of the finite-mode no-go: the trajectory *remains* in the
subspace `V` for all time. -/
def FixedBandTrajectory (u : ℝ → E) (V : Submodule ℝ E) : Prop :=
  ∀ t, u t ∈ V

/-- The hypothesis Track P actually has: the *datum* lies in `V`.  Nothing is
said about later times. -/
def FiniteBandDatum (u : ℝ → E) (V : Submodule ℝ E) : Prop :=
  u 0 ∈ V

/-- **F4 (i).**  A fixed-band trajectory has, in particular, a finite-band
datum: evaluate at `t = 0`. -/
theorem FixedBandTrajectory.finiteBandDatum {u : ℝ → E} {V : Submodule ℝ E}
    (h : FixedBandTrajectory u V) : FiniteBandDatum u V :=
  h 0

end BandDistinction

/-- **F4 (ii).  The converse fails.**  The curve `t ↦ (1, t)` starts on the
first-coordinate axis and leaves it at every later time: a finite-band datum
with no fixed-band trajectory.  This two-dimensional toy is the exact logical
shape of the Track-P situation — the Galerkin datum lies in the band, the
Galerkin tail pushes the true solution out — reduced to its minimum. -/
theorem exists_finiteBandDatum_not_fixedBandTrajectory :
    ∃ (u : ℝ → ℝ × ℝ) (V : Submodule ℝ (ℝ × ℝ)),
      FiniteBandDatum u V ∧ ¬ FixedBandTrajectory u V := by
  refine ⟨fun t : ℝ => ((1 : ℝ), t), LinearMap.ker (LinearMap.snd ℝ ℝ ℝ),
    ?_, ?_⟩
  · simp [FiniteBandDatum, LinearMap.mem_ker]
  · intro hfix
    have h1 := hfix 1
    simp only [LinearMap.mem_ker, LinearMap.snd_apply] at h1
    exact one_ne_zero h1

/-- **F4, the scope of the existing no-go restated.**  A packaged
fixed-bandwidth candidate whose trajectory satisfies `FixedBandTrajectory`
has a finite-band datum *and* an empty breakdown set: the first conjunct is
F4 (i), the second is the packaged no-go of `CertificateLayer.lean`,
delegated verbatim.  Read with `exists_finiteBandDatum_not_fixedBandTrajectory`
this states the scope exactly: the excluded class sits strictly inside the
finite-band-datum class, and the no-go is silent about the rest of it — in
particular about a Galerkin datum whose true Navier–Stokes evolution leaves
the band, which is Track P's whole subject. -/
theorem FixedBandwidthCandidate.fixedBand_scope
    {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    [CompleteSpace E]
    (candidate : FixedBandwidthCandidate E) {V : Submodule ℝ E}
    (hband : FixedBandTrajectory candidate.trajectory V) :
    FiniteBandDatum candidate.trajectory V ∧
      {T : ℝ | IsBreakdownCandidate candidate.trajectory T} = ∅ :=
  ⟨hband.finiteBandDatum, candidate.breakdown_times_empty⟩

/-! ## F5 — the homogeneous-norm ladder

On the mean-zero torus every retained wavenumber has `|k| ≥ 1`, so the
spectral weights `w = |k|²` satisfy `w ≥ 1` and the weighted coefficient sums
`Σ wⁿ c` are monotone in `n`.  This is `‖f‖_m ≤ ‖f‖_n` for `m ≤ n` in exactly
the finite-sum form the control inequality uses — it is how `-ν‖w‖_5²` is
absorbed into `-ν‖w‖_4²` in the viscous term, and how the final `H^3`
conclusion follows from the `Ḣ^4` bound. -/

section CoefficientLadder

/-- **F5, one step.**  With weights `≥ 1` and nonnegative coefficients, one
more power of the weight can only increase the sum: termwise
`wⁿ c ≤ wⁿ⁺¹ c`. -/
theorem weighted_sum_succ_mono {ι : Type*} (s : Finset ι) (w c : ι → ℝ)
    (n : ℕ) (hw : ∀ i ∈ s, 1 ≤ w i) (hc : ∀ i ∈ s, 0 ≤ c i) :
    ∑ i ∈ s, w i ^ n * c i ≤ ∑ i ∈ s, w i ^ (n + 1) * c i := by
  refine Finset.sum_le_sum fun i hi => ?_
  have h1 := hw i hi
  have h0 : (0 : ℝ) ≤ w i := zero_le_one.trans h1
  have hpow : w i ^ n ≤ w i ^ (n + 1) := by
    rw [pow_succ]
    nlinarith [pow_nonneg h0 n]
  exact mul_le_mul_of_nonneg_right hpow (hc i hi)

/-- **F5, the ladder.**  The full monotone ladder `m ≤ n` implies
`Σ wᵐ c ≤ Σ wⁿ c`, by iterating the single step.  With `w k = |k|²` and
`c k = |û_k|²` this is `‖f‖_{Ḣ^m} ≤ ‖f‖_{Ḣ^n}` for mean-zero finite Fourier
data. -/
theorem weighted_sum_pow_mono {ι : Type*} (s : Finset ι) (w c : ι → ℝ)
    {m n : ℕ} (hmn : m ≤ n) (hw : ∀ i ∈ s, 1 ≤ w i)
    (hc : ∀ i ∈ s, 0 ≤ c i) :
    ∑ i ∈ s, w i ^ m * c i ≤ ∑ i ∈ s, w i ^ n * c i := by
  induction n, hmn using Nat.le_induction with
  | base => exact le_rfl
  | succ n _ ih => exact ih.trans (weighted_sum_succ_mono s w c n hw hc)

end CoefficientLadder

/-! ## F6 — the slab composition shell

The Track-P control inequality is a Riccati subsolution statement for
`R t ≥ ‖w t‖_{Ḣ^4}` with datum error zero — the Galerkin datum *is* the
datum.  The instantiation below is deliberately tiny: it shows the estimate
layer's output plugs directly into `ControlODE.riccati_comparison`, with
`a = -ν + 9(K1 + K2)`, `b = 135·Adot` and `e = sup ‖e(t)‖_4` read off the
derivation in the module docstring.  Supplying the subsolution hypothesis for
the true Navier–Stokes error is the business of EXT-P1/2/3 plus the interval
certificate, and remains outside Lean. -/

section SlabShell

/-- **F6.**  A Riccati subsolution starting from zero never overtakes a
supersolution starting nonnegative: the Track-P error bound on the slab
`[0, h]`, conditional on the control inequality, is an instance of the
existing quadratic comparison.  Nothing here is specific to Navier–Stokes —
that is the point: the PDE content is entirely in the hypothesis `hsub`. -/
theorem trackP_slab_error_le {R R' S S' : ℝ → ℝ} {a b e h : ℝ}
    (hb : 0 ≤ b)
    (hR : ContinuousOn R (Icc 0 h)) (hS : ContinuousOn S (Icc 0 h))
    (hR' : ∀ t ∈ Ico (0 : ℝ) h, HasDerivWithinAt R (R' t) (Ici t) t)
    (hS' : ∀ t ∈ Ico (0 : ℝ) h, HasDerivWithinAt S (S' t) (Ici t) t)
    (hR0 : R 0 = 0) (hS0 : 0 ≤ S 0)
    (hsub : ∀ t ∈ Ico (0 : ℝ) h, R' t ≤ a * R t + b * R t ^ 2 + e)
    (hsuper : ∀ t ∈ Ico (0 : ℝ) h, a * S t + b * S t ^ 2 + e ≤ S' t) :
    ∀ t ∈ Icc (0 : ℝ) h, R t ≤ S t :=
  riccati_comparison hb hR hS hR' hS' hsub hsuper (hR0.le.trans hS0)

end SlabShell

end

end NSSingularity
