/-
Stage 1, item **F-3** of `docs/formalization_map.md`:
the velocity-recovery formulas **E-14** and the divergence cancellation
**E-15** of `docs/equation_audit.md`.

Mathematical content.  Let `ψ₁ = ψ₁(r,z)` be the normalized angular stream
function of an axisymmetric flow.  E-14 recovers the meridional velocity
components by

  `uʳ(r,z) = -r ∂_z ψ₁(r,z)`,
  `u^z(r,z) = 2 ψ₁(r,z) + r ∂_r ψ₁(r,z)`,

and E-15 asserts that the *physical* (three-dimensional, axisymmetric)
divergence of the recovered field vanishes:

  `∂_r uʳ + uʳ / r + ∂_z u^z = 0`   for `r ≠ 0`.

Term by term,

  `∂_r uʳ = -∂_z ψ₁ - r ∂_r∂_z ψ₁`,
  `uʳ / r = -∂_z ψ₁`,
  `∂_z u^z = 2 ∂_z ψ₁ + r ∂_z∂_r ψ₁`,

so the sum vanishes **only** because the two mixed second partials agree.
That Schwarz/Clairaut symmetry is the entire mathematical content of E-15;
it is isolated below in `mixed_partial_comm`, which is the single place
where the `C²` hypothesis is consumed.

Representation choices (kept deliberately close to `ClayStatement.lean`).
* The meridional half-plane is modelled by the *plain product* `ℝ × ℝ`,
  with `p.1` the radial coordinate `r` and `p.2` the axial coordinate `z`.
  `ClayStatement.lean` uses `EuclideanSpace ℝ (Fin 3)` with
  `coordVector i = EuclideanSpace.single i 1`; here the ambient space is
  two-dimensional and *not* a Euclidean coordinate space, because `(r,z)`
  are curvilinear coordinates and no inner product on them is meaningful
  for E-15.  `ℝ × ℝ` carries the same (equivalent) norm topology, so
  `fderiv` and `ContDiff` mean exactly the same thing, while `(1,0)` and
  `(0,1)` are literally the coordinate directions.
* Exactly as in `ClayStatement.lean`, partial derivatives are Fréchet
  derivatives evaluated on a coordinate direction:
  `partialR f p = fderiv ℝ f p (1,0)`, `partialZ f p = fderiv ℝ f p (0,1)`.
  `fderiv` is junk-valued at points of non-differentiability, which is
  harmless because every statement below carries a `ContDiff ℝ 2`
  hypothesis in the same proposition.
* Division `uR ψ p / p.1` is Lean's total division (`x / 0 = 0`).  The main
  theorem therefore carries the hypothesis `p.1 ≠ 0`, matching the `r ≠ 0`
  side condition of E-15.  The axis is treated separately by
  `uROverR`, the continuous extension of `uʳ/r`, for which the identity
  holds at *every* point including `r = 0` (see
  `divergence_of_recovered_velocity_eq_zero'`).

Scope caveat.  Everything here is an identity about a *chosen coordinate
representation of a scalar function*.  Nothing in this file mentions the
Navier–Stokes equations, axisymmetry of an actual 3D field, or the
cylindrical-to-Cartesian correspondence (E-18/E-24); see `formal/README.md`.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib

namespace NSSingularity

noncomputable section

open ContinuousLinearMap

/-- The radial coordinate direction `∂/∂r` in the meridional plane `ℝ × ℝ`. -/
def dirR : ℝ × ℝ := (1, 0)

/-- The axial coordinate direction `∂/∂z` in the meridional plane `ℝ × ℝ`. -/
def dirZ : ℝ × ℝ := (0, 1)

/-- Radial partial derivative of a scalar function on the meridional plane.
Same convention as `NSSingularity.partialDeriv` in `ClayStatement.lean`:
a Fréchet derivative evaluated on a coordinate direction. -/
def partialR (f : ℝ × ℝ → ℝ) (p : ℝ × ℝ) : ℝ := fderiv ℝ f p dirR

/-- Axial partial derivative of a scalar function on the meridional plane. -/
def partialZ (f : ℝ × ℝ → ℝ) (p : ℝ × ℝ) : ℝ := fderiv ℝ f p dirZ

/-- **E-14, radial component.**  `uʳ = -r ∂_z ψ₁`. -/
def uR (ψ : ℝ × ℝ → ℝ) (p : ℝ × ℝ) : ℝ := -(p.1 * partialZ ψ p)

/-- **E-14, axial component.**  `u^z = 2 ψ₁ + r ∂_r ψ₁`. -/
def uZ (ψ : ℝ × ℝ → ℝ) (p : ℝ × ℝ) : ℝ := 2 * ψ p + p.1 * partialR ψ p

/-- The continuous extension of `uʳ / r` to the whole plane.  For `r ≠ 0`
this equals `uR ψ p / p.1` (`uROverR_eq_div`); at `r = 0` it is the limit
value `-∂_z ψ₁(0,z)`, which is what the physical axisymmetric divergence
uses on the axis (cf. E-16). -/
def uROverR (ψ : ℝ × ℝ → ℝ) (p : ℝ × ℝ) : ℝ := -(partialZ ψ p)

section Regularity

variable {ψ : ℝ × ℝ → ℝ}

/-- A `C²` function has a differentiable Fréchet-derivative map. -/
theorem differentiable_fderiv_of_contDiff_two (hψ : ContDiff ℝ 2 ψ) :
    Differentiable ℝ (fderiv ℝ ψ) :=
  (hψ.fderiv_right (m := 1) (by norm_num)).differentiable (by norm_num)

/-- Differentiating `q ↦ Dψ(q) w` in the direction `v` produces the second
Fréchet derivative `D²ψ(p) v w`.  This is pure bookkeeping (application of
a fixed vector is a continuous linear map), stated once so that the
`fderiv`-of-`fderiv` shuffling does not clutter the real proofs. -/
theorem hasFDerivAt_apply_const (hψ : ContDiff ℝ 2 ψ) (w p : ℝ × ℝ) :
    HasFDerivAt (fun q : ℝ × ℝ => fderiv ℝ ψ q w)
      ((fderiv ℝ (fderiv ℝ ψ) p).flip w) p := by
  have h1 : HasFDerivAt (fderiv ℝ ψ) (fderiv ℝ (fderiv ℝ ψ) p) p :=
    (differentiable_fderiv_of_contDiff_two hψ p).hasFDerivAt
  simpa using h1.clm_apply (hasFDerivAt_const w p)

theorem fderiv_apply_const (hψ : ContDiff ℝ 2 ψ) (v w p : ℝ × ℝ) :
    fderiv ℝ (fun q : ℝ × ℝ => fderiv ℝ ψ q w) p v
      = fderiv ℝ (fderiv ℝ ψ) p v w := by
  rw [(hasFDerivAt_apply_const hψ w p).fderiv]
  simp

/-- **The Schwarz / Clairaut step — the only place smoothness is used.**

For `C²` data the two mixed second partials in E-15 agree,
`∂_r ∂_z ψ₁ = ∂_z ∂_r ψ₁`.  Everything else in this file is the product
rule.  The mathlib input is
`ContDiffAt.isSymmSndFDerivAt : ContDiffAt 𝕜 n f x → minSmoothness 𝕜 2 ≤ n →
IsSymmSndFDerivAt 𝕜 f x`, whose conclusion unfolds to
`∀ v w, fderiv 𝕜 (fderiv 𝕜 f) x v w = fderiv 𝕜 (fderiv 𝕜 f) x w v`. -/
theorem mixed_partial_comm (hψ : ContDiff ℝ 2 ψ) (p : ℝ × ℝ) :
    partialR (partialZ ψ) p = partialZ (partialR ψ) p := by
  have hsymm : IsSymmSndFDerivAt ℝ ψ p :=
    hψ.contDiffAt.isSymmSndFDerivAt (by simp)
  have h1 : partialR (partialZ ψ) p = fderiv ℝ (fderiv ℝ ψ) p dirR dirZ :=
    fderiv_apply_const hψ dirR dirZ p
  have h2 : partialZ (partialR ψ) p = fderiv ℝ (fderiv ℝ ψ) p dirZ dirR :=
    fderiv_apply_const hψ dirZ dirR p
  rw [h1, h2]
  exact hsymm dirR dirZ

end Regularity

section Expansion

variable {ψ : ℝ × ℝ → ℝ}

/-- Product rule for the radial component: `∂_r uʳ = -∂_z ψ₁ - r ∂_r∂_z ψ₁`
(first bracket of E-15). -/
theorem partialR_uR (hψ : ContDiff ℝ 2 ψ) (p : ℝ × ℝ) :
    partialR (uR ψ) p = -(partialZ ψ p) - p.1 * partialR (partialZ ψ) p := by
  have hg : HasFDerivAt (fun q : ℝ × ℝ => fderiv ℝ ψ q dirZ)
      ((fderiv ℝ (fderiv ℝ ψ) p).flip dirZ) p := hasFDerivAt_apply_const hψ dirZ p
  have hf : HasFDerivAt (fun q : ℝ × ℝ => q.1) (ContinuousLinearMap.fst ℝ ℝ ℝ) p :=
    hasFDerivAt_fst
  have h : HasFDerivAt (uR ψ)
      (-(p.1 • (fderiv ℝ (fderiv ℝ ψ) p).flip dirZ
        + fderiv ℝ ψ p dirZ • ContinuousLinearMap.fst ℝ ℝ ℝ)) p := (hf.mul hg).neg
  have hR : partialR (uR ψ) p
      = -(p.1 * fderiv ℝ (fderiv ℝ ψ) p dirR dirZ + fderiv ℝ ψ p dirZ) := by
    show fderiv ℝ (uR ψ) p dirR = _
    rw [h.fderiv]
    simp [dirR]
  have hcomm : partialR (partialZ ψ) p = fderiv ℝ (fderiv ℝ ψ) p dirR dirZ :=
    fderiv_apply_const hψ dirR dirZ p
  have hz : partialZ ψ p = fderiv ℝ ψ p dirZ := rfl
  linear_combination hR + hz + p.1 * hcomm

/-- Product rule for the axial component:
`∂_z u^z = 2 ∂_z ψ₁ + r ∂_z∂_r ψ₁` (third bracket of E-15). -/
theorem partialZ_uZ (hψ : ContDiff ℝ 2 ψ) (p : ℝ × ℝ) :
    partialZ (uZ ψ) p = 2 * partialZ ψ p + p.1 * partialZ (partialR ψ) p := by
  have hψ1 : HasFDerivAt ψ (fderiv ℝ ψ p) p :=
    (hψ.differentiable (by norm_num) p).hasFDerivAt
  have hg : HasFDerivAt (fun q : ℝ × ℝ => fderiv ℝ ψ q dirR)
      ((fderiv ℝ (fderiv ℝ ψ) p).flip dirR) p := hasFDerivAt_apply_const hψ dirR p
  have hf : HasFDerivAt (fun q : ℝ × ℝ => q.1) (ContinuousLinearMap.fst ℝ ℝ ℝ) p :=
    hasFDerivAt_fst
  have h : HasFDerivAt (uZ ψ)
      ((2 : ℝ) • fderiv ℝ ψ p
        + (p.1 • (fderiv ℝ (fderiv ℝ ψ) p).flip dirR
          + fderiv ℝ ψ p dirR • ContinuousLinearMap.fst ℝ ℝ ℝ)) p :=
    (hψ1.const_mul (2 : ℝ)).add (hf.mul hg)
  have hZ : partialZ (uZ ψ) p
      = 2 * fderiv ℝ ψ p dirZ + p.1 * fderiv ℝ (fderiv ℝ ψ) p dirZ dirR := by
    show fderiv ℝ (uZ ψ) p dirZ = _
    rw [h.fderiv]
    simp [dirZ]
  have hcomm : partialZ (partialR ψ) p = fderiv ℝ (fderiv ℝ ψ) p dirZ dirR :=
    fderiv_apply_const hψ dirZ dirR p
  have hz : partialZ ψ p = fderiv ℝ ψ p dirZ := rfl
  linear_combination hZ - 2 * hz - p.1 * hcomm

end Expansion

/-- **F-3 / E-15, axis-inclusive form.**

With `uʳ/r` replaced by its continuous extension `uROverR ψ = -∂_z ψ₁`, the
axisymmetric divergence of the recovered velocity field vanishes at *every*
point of the meridional plane, including the axis `r = 0`.  This is the
sharpest form: no `r ≠ 0` hypothesis is needed. -/
theorem divergence_of_recovered_velocity_eq_zero' (ψ : ℝ × ℝ → ℝ)
    (hψ : ContDiff ℝ 2 ψ) (p : ℝ × ℝ) :
    partialR (uR ψ) p + uROverR ψ p + partialZ (uZ ψ) p = 0 := by
  rw [partialR_uR hψ p, partialZ_uZ hψ p, ← mixed_partial_comm hψ p]
  show -(partialZ ψ p) - p.1 * partialR (partialZ ψ) p + -(partialZ ψ p)
      + (2 * partialZ ψ p + p.1 * partialR (partialZ ψ) p) = 0
  ring

/-- Off the axis the continuous extension really is the quotient `uʳ/r`. -/
theorem uROverR_eq_div (ψ : ℝ × ℝ → ℝ) (p : ℝ × ℝ) (hr : p.1 ≠ 0) :
    uROverR ψ p = uR ψ p / p.1 := by
  rw [uROverR, uR]
  field_simp

/-- **F-3 / E-15, the literal statement.**

For `C²` normalized stream function `ψ₁` and any point with `r ≠ 0`, the
velocity field recovered by E-14 satisfies the physical axisymmetric
divergence-free condition

  `∂_r uʳ + uʳ/r + ∂_z u^z = 0`.

The `r ≠ 0` hypothesis is needed only because Lean's `/` is total
(`x / 0 = 0`); see `divergence_of_recovered_velocity_eq_zero'` for the form
that also covers the axis. -/
theorem divergence_of_recovered_velocity_eq_zero (ψ : ℝ × ℝ → ℝ)
    (hψ : ContDiff ℝ 2 ψ) (p : ℝ × ℝ) (hr : p.1 ≠ 0) :
    partialR (uR ψ) p + uR ψ p / p.1 + partialZ (uZ ψ) p = 0 := by
  rw [← uROverR_eq_div ψ p hr]
  exact divergence_of_recovered_velocity_eq_zero' ψ hψ p

end

end NSSingularity
