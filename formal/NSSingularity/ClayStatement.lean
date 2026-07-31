/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
/-
Stage 0 of `LEAN4_VERIFICATION_POLICY.md`: fix the Clay Navier–Stokes
statements as Lean `Prop`s.

This file contains DEFINITIONS ONLY — no proofs, no `sorry`, no axioms.
Its purpose is to pin down, in machine-checkable form, exactly which
propositions the research program ultimately targets, so that an
independent mathematician can review the correspondence with the official
problem statement (C. Fefferman, "Existence and smoothness of the
Navier–Stokes equation", Clay Mathematics Institute).  The natural-language
correspondence table lives in `formal/README.md`.

Conventions:
* Physical space is `ℝ³` as `EuclideanSpace ℝ (Fin 3)`.
* Time-dependent fields are curried `ℝ → Space → _` with time first.
* Derivatives are Fréchet derivatives evaluated on coordinate directions;
  a field that is not differentiable makes `fderiv` junk-valued, which is
  harmless here because every use below is guarded by a smoothness
  hypothesis in the same statement.
* The official statement requires solutions smooth on `ℝ³ × [0,∞)`;
  we encode this as joint `ContDiffOn` on `Set.Ici 0 ×ˢ Set.univ`.

Faithfulness caveats (reviewed in `formal/README.md`):
* Fefferman's conditions (4)/(5) bound derivatives of all orders by
  `C_{α,K} (1+|x|)^{-K}`; we quantify over the total derivative order
  `n : ℕ` via `iteratedFDeriv`, which dominates the multi-index form.
* The pressure is required smooth on the same half-space; no decay is
  imposed on `p`, matching the official text.
-/
import Mathlib

namespace NSSingularity

noncomputable section

open MeasureTheory Set

/-- Physical space `ℝ³`. -/
abbrev Space : Type := EuclideanSpace ℝ (Fin 3)

/-- Time-dependent velocity field (time first). -/
abbrev VelocityField : Type := ℝ → Space → Space

/-- Time-dependent scalar pressure field. -/
abbrev PressureField : Type := ℝ → Space → ℝ

/-- The `i`-th coordinate unit vector of `ℝ³`. -/
def coordVector (i : Fin 3) : Space := EuclideanSpace.single i 1

/-- Partial derivative of a scalar field along the `i`-th coordinate. -/
def partialDeriv (i : Fin 3) (f : Space → ℝ) (x : Space) : ℝ :=
  fderiv ℝ f x (coordVector i)

/-- Divergence of a (static) vector field on `ℝ³`. -/
def divergence (u : Space → Space) (x : Space) : ℝ :=
  ∑ i, partialDeriv i (fun y => u y i) x

/-- Scalar Laplacian on `ℝ³`. -/
def laplacianScalar (f : Space → ℝ) (x : Space) : ℝ :=
  ∑ i, partialDeriv i (fun y => partialDeriv i f y) x

/-- A static vector field is divergence-free. -/
def DivergenceFree (u : Space → Space) : Prop :=
  ∀ x : Space, divergence u x = 0

/-- Smooth (C^∞) on all of `ℝ³`. -/
def SmoothStatic (u : Space → Space) : Prop :=
  ContDiff ℝ (⊤ : ℕ∞) u

/-- Fefferman condition (4): every spatial derivative decays faster than
any inverse power of `|x|`.  We quantify over the total derivative order
`n`, which dominates the multi-index formulation. -/
def RapidlyDecaying (u : Space → Space) : Prop :=
  ∀ n K : ℕ, ∃ C : ℝ, ∀ x : Space,
    ‖iteratedFDeriv ℝ n u x‖ ≤ C * (1 + ‖x‖) ^ (-(K : ℝ))

/-- Admissible Clay initial datum on `ℝ³`: smooth, divergence-free,
rapidly decaying (official conditions (1) and (4)). -/
def AdmissibleData (u₀ : Space → Space) : Prop :=
  SmoothStatic u₀ ∧ DivergenceFree u₀ ∧ RapidlyDecaying u₀

/-- Joint C^∞ smoothness of a time-dependent vector field on
`[0,∞) × ℝ³` (official condition (6)/(7) side condition). -/
def SmoothOnHalfTime (u : VelocityField) : Prop :=
  ContDiffOn ℝ (⊤ : ℕ∞) (fun q : ℝ × Space => u q.1 q.2)
    (Set.Ici (0 : ℝ) ×ˢ Set.univ)

/-- Joint C^∞ smoothness of a scalar field on `[0,∞) × ℝ³`. -/
def SmoothScalarOnHalfTime (p : PressureField) : Prop :=
  ContDiffOn ℝ (⊤ : ℕ∞) (fun q : ℝ × Space => p q.1 q.2)
    (Set.Ici (0 : ℝ) ×ˢ Set.univ)

/-- Fefferman condition (5) for a force: every space-time derivative is
bounded by `C (1 + |x| + t)^{-K}` on the physical half-time domain.
We quantify over joint derivative orders through the curried map. -/
def AdmissibleForce (f : VelocityField) : Prop :=
  SmoothOnHalfTime f ∧
  ∀ n K : ℕ, ∃ C : ℝ, ∀ t : ℝ, 0 ≤ t → ∀ x : Space,
    ‖iteratedFDeriv ℝ n (fun q : ℝ × Space => f q.1 q.2) (t, x)‖
      ≤ C * (1 + ‖x‖ + t) ^ (-(K : ℝ))

/-- The incompressible Navier–Stokes momentum equation with viscosity `ν`
and force `f`, imposed pointwise for `t > 0` (official equation (1)). -/
def SatisfiesMomentum (ν : ℝ) (u : VelocityField) (p : PressureField)
    (f : VelocityField) : Prop :=
  ∀ t : ℝ, 0 < t → ∀ x : Space, ∀ i : Fin 3,
    deriv (fun s => u s x i) t
      + ∑ j, u t x j * partialDeriv j (fun y => u t y i) x
      = ν * laplacianScalar (fun y => u t y i) x
        - partialDeriv i (fun y => p t y) x
        + f t x i

/-- Incompressibility at every nonnegative time (official equation (2)). -/
def SatisfiesIncompressibility (u : VelocityField) : Prop :=
  ∀ t : ℝ, 0 ≤ t → DivergenceFree (u t)

/-- Bounded energy (official condition (7)). -/
def BoundedEnergy (u : VelocityField) : Prop :=
  ∃ C : ℝ, ∀ t : ℝ, 0 ≤ t → (∫ x : Space, ‖u t x‖ ^ 2) ≤ C

/-- A Clay-admissible global smooth solution for datum `u₀`, force `f`,
viscosity `ν`. -/
def IsGlobalSmoothSolution (ν : ℝ) (u₀ : Space → Space) (f : VelocityField)
    (u : VelocityField) (p : PressureField) : Prop :=
  SmoothOnHalfTime u ∧ SmoothScalarOnHalfTime p ∧
  (∀ x : Space, u 0 x = u₀ x) ∧
  SatisfiesMomentum ν u p f ∧
  SatisfiesIncompressibility u ∧
  BoundedEnergy u

/-- **Clay statement (A)** — whole-space existence and smoothness:
for every admissible datum and every positive viscosity there is a global
smooth bounded-energy solution with zero force. -/
def ClayWholeSpaceRegularity : Prop :=
  ∀ ν : ℝ, 0 < ν → ∀ u₀ : Space → Space, AdmissibleData u₀ →
    ∃ u : VelocityField, ∃ p : PressureField,
      IsGlobalSmoothSolution ν u₀ (fun _ _ => 0) u p

/-- **Clay statement (C)** — whole-space breakdown: there exist a positive
viscosity, an admissible datum, and an admissible force for which no global
smooth bounded-energy solution exists. -/
def ClayWholeSpaceBreakdown : Prop :=
  ∃ ν : ℝ, 0 < ν ∧ ∃ u₀ : Space → Space, AdmissibleData u₀ ∧
    ∃ f : VelocityField, AdmissibleForce f ∧
      ∀ u : VelocityField, ∀ p : PressureField,
        ¬ IsGlobalSmoothSolution ν u₀ f u p

/-- **This repository's preferred target** — unforced whole-space
breakdown, STRICTLY STRONGER than (C): the force is identically zero. -/
def UnforcedWholeSpaceBreakdown : Prop :=
  ∃ ν : ℝ, 0 < ν ∧ ∃ u₀ : Space → Space, AdmissibleData u₀ ∧
    ∀ u : VelocityField, ∀ p : PressureField,
      ¬ IsGlobalSmoothSolution ν u₀ (fun _ _ => 0) u p

section Periodic

/-- Spatial periodicity with respect to the integer lattice `ℤ³`
(official condition (8), period 1 in each coordinate). -/
def SpacePeriodic (u : Space → Space) : Prop :=
  ∀ x : Space, ∀ i : Fin 3, u (x + coordVector i) = u x

/-- Periodic admissible datum: smooth, divergence-free, `ℤ³`-periodic
(official conditions (10) and (8); no decay is required in the periodic
setting). -/
def AdmissiblePeriodicData (u₀ : Space → Space) : Prop :=
  SmoothStatic u₀ ∧ DivergenceFree u₀ ∧ SpacePeriodic u₀

/-- A global smooth spatially periodic solution (bounded energy is
automatic on the torus, but we keep the official phrasing: `u` smooth on
`ℝ³ × [0,∞)` and periodic at every time). -/
def IsGlobalSmoothPeriodicSolution (ν : ℝ) (u₀ : Space → Space)
    (f : VelocityField) (u : VelocityField) (p : PressureField) : Prop :=
  SmoothOnHalfTime u ∧ SmoothScalarOnHalfTime p ∧
  (∀ x : Space, u 0 x = u₀ x) ∧
  (∀ t : ℝ, 0 ≤ t → SpacePeriodic (u t)) ∧
  SatisfiesMomentum ν u p f ∧
  SatisfiesIncompressibility u

/-- **Clay statement (B)** — periodic existence and smoothness. -/
def ClayPeriodicRegularity : Prop :=
  ∀ ν : ℝ, 0 < ν → ∀ u₀ : Space → Space, AdmissiblePeriodicData u₀ →
    ∃ u : VelocityField, ∃ p : PressureField,
      IsGlobalSmoothPeriodicSolution ν u₀ (fun _ _ => 0) u p

/-- A `ℤ³`-periodic smooth force with all space-time derivatives bounded
uniformly on the half-time domain (periodic analogue of condition (5),
following the official statement (9)). -/
def AdmissiblePeriodicForce (f : VelocityField) : Prop :=
  SmoothOnHalfTime f ∧ (∀ t : ℝ, 0 ≤ t → SpacePeriodic (f t)) ∧
  ∀ n : ℕ, ∃ C : ℝ, ∀ t : ℝ, 0 ≤ t → ∀ x : Space,
    ‖iteratedFDeriv ℝ n (fun q : ℝ × Space => f q.1 q.2) (t, x)‖ ≤ C

/-- **Clay statement (D)** — periodic breakdown. -/
def ClayPeriodicBreakdown : Prop :=
  ∃ ν : ℝ, 0 < ν ∧ ∃ u₀ : Space → Space, AdmissiblePeriodicData u₀ ∧
    ∃ f : VelocityField, AdmissiblePeriodicForce f ∧
      ∀ u : VelocityField, ∀ p : PressureField,
        ¬ IsGlobalSmoothPeriodicSolution ν u₀ f u p

end Periodic

end

end NSSingularity
