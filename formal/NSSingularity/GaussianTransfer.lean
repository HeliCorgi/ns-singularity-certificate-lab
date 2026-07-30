/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under the Apache 2.0 licence.
-/
import Mathlib
import NSSingularity.ControlODE

/-!
# The Gaussian transfer layer

The optimisation half of this repository works on finite-`C^k` spline
surrogates, because a spline has finitely many coefficients and interval
arithmetic can certify statements about finitely many numbers.  A spline
surrogate is *never* a Clay candidate — that is an absolute rule of the
project, not a modelling convenience: the Clay statement quantifies over
smooth fields, and a surrogate optimum acquires smooth meaning only through a
*transfer*.  Mollify the spline to a smooth field, and prove that the
generation functional `J` moves by less than the margin with which the
surrogate optimum was certified; only the mollified smooth field is ever a
candidate, and only after that transfer bound is proved.  This file supplies
the layers of the transfer that are pure algebra or finite inequalities,
together with the torus specialisation of the control ODE that the Python
payload cites by a single theorem name.

## `G1` — closure of the Gaussian–Hermite class under differentiation

The smooth side of the transfer lives in the span of the functions
`x ↦ p(x) · exp(-α x²)` with `p` a polynomial.  That this class is closed
under differentiation is usually left as a slogan; here it is a theorem with
the witness polynomial exhibited:

* `hasDerivAt_poly_mul_gaussian` — the derivative of `p(x) exp(-α x²)` is
  `(p'(x) - 2 α x p(x)) exp(-α x²)`, another member of the class;
* `gaussianDerivPoly`, `gaussianDerivPoly_eval` — the witness
  `q = p' - C(2α)·X·p` and the identity `q(x) = p'(x) - 2 α x p(x)`;
* `hasDerivAt_poly_gaussian` — the closure restated with the witness in the
  derivative slot, which is the citable form.

The downstream consequence is that every space derivative of a mollified
candidate is obtained by *coefficient algebra* — exact in the certified
arithmetic — and never by numerical differentiation.  No sign condition is
imposed on `α`: positivity matters for integrability, which is not at stake in
a pointwise derivative computation, and demanding it would only weaken the
algebra.

## `G2` — the `J`-continuity bricks

"`J` is continuous along a mollification" is the licence to transfer a
surrogate optimum, and its analytic content reduces to three pointwise
inequalities applied under an integral sign:

* `abs_cube_sub_cube_le` — the scalar brick,
  `|a³ - b³| ≤ 3 (max a b)² |a - b|` for nonnegative `a`, `b`;
* `abs_norm_cube_sub_norm_cube_le` — the brick at the integrand of `J`,
  `| ‖u‖³ - ‖v‖³ | ≤ 3 (‖u‖ + ‖v‖)² ‖u - v‖`, obtained from the scalar brick
  through the reverse triangle inequality;
* `norm_smul_sub_smul_le` — the vector brick for the nonlinearity `|u| u`
  appearing in the mollification bound,
  `‖ ‖u‖•u - ‖v‖•v ‖ ≤ (‖u‖ + ‖v‖) ‖u - v‖`.

Everything here is pointwise.  The passage from these bricks to continuity of
the integral functional needs a dominated-convergence layer that is out of
scope for the same reason as the integration by parts in `L3Generation.lean`;
the honest division of labour is recorded there and not repeated by weakening
anything below.

## `G3` — the torus control bound

`ControlODE.lean` explains why the whole-space control ODE has no linear
damping: the Laplacian on `ℝ³` has no spectral gap.  On the torus the
admissible wavenumbers are bounded away from zero, the Poincaré inequality
holds, and the viscous term contributes a genuine `-ν R` to the linear
coefficient, against which the estimate layer contributes `+a₀ R`; the net
linear coefficient `-ν + a₀` carries no fixed sign.  `torus_control_bound` is
`riccati_comparison` with that substitution performed, `R 0 = 0` as the
Track-P payload has it, and the Track-P sign conventions (`ν > 0`, `b ≥ 0`,
`e ≥ 0`) recorded in the hypotheses, so that the Python payload can cite one
theorem name instead of a substitution recipe.

## Scope

This file contains no `sorry`, no `admit`, and introduces no axioms.  The
viscosity is fixed and positive throughout the project, and nothing here,
alone or in combination with the rest of the library, is a statement about
the Clay problem.
-/

namespace NSSingularity

noncomputable section

open Set Polynomial

/-! ## G1 — the Gaussian–Hermite class is closed under differentiation -/

section GaussianClosure

/-- **G1, the derivative computation.**  For every polynomial `p` and every
`α`,

`d/dx ( p(x) · exp(-α x²) ) = ( p'(x) - 2 α x p(x) ) · exp(-α x²)`.

The derivative of a polynomial times a Gaussian is a polynomial times the
*same* Gaussian, with the polynomial in closed form.  This is the product
rule, the polynomial derivative, and the chain rule through `exp` — nothing
else — but it is the identity that makes symbolic differentiation of a
Gaussian–Hermite candidate exact rather than approximate. -/
theorem hasDerivAt_poly_mul_gaussian (p : Polynomial ℝ) (α x : ℝ) :
    HasDerivAt (fun t => p.eval t * Real.exp (-(α * t ^ 2)))
      ((p.derivative.eval x - 2 * α * x * p.eval x) * Real.exp (-(α * x ^ 2))) x := by
  have hp : HasDerivAt (fun t : ℝ => p.eval t) (p.derivative.eval x) x := p.hasDerivAt x
  have hsq : HasDerivAt (fun t : ℝ => t ^ 2) (2 * x) x := by
    simpa using hasDerivAt_pow 2 x
  have hexp : HasDerivAt (fun t : ℝ => Real.exp (-(α * t ^ 2)))
      (Real.exp (-(α * x ^ 2)) * -(α * (2 * x))) x := ((hsq.const_mul α).neg).exp
  have heq : (p.derivative.eval x - 2 * α * x * p.eval x) * Real.exp (-(α * x ^ 2))
      = p.derivative.eval x * Real.exp (-(α * x ^ 2))
        + p.eval x * (Real.exp (-(α * x ^ 2)) * -(α * (2 * x))) := by ring
  rw [heq]
  exact hp.mul hexp

/-- The witness of the closure: differentiation acts on the polynomial part of
`p(x) exp(-α x²)` as `p ↦ p' - 2αX·p`.  Keeping this as a named polynomial —
rather than as an evaluated expression — is the point of `G1`: the certificate
manipulates coefficient vectors, and this is the operator it iterates to reach
higher derivatives without ever leaving the class. -/
def gaussianDerivPoly (α : ℝ) (p : Polynomial ℝ) : Polynomial ℝ :=
  p.derivative - C (2 * α) * X * p

/-- The witness evaluates to the coefficient the product rule produces:
`q(x) = p'(x) - 2 α x p(x)`.  Pure `Polynomial.eval` algebra. -/
theorem gaussianDerivPoly_eval (α : ℝ) (p : Polynomial ℝ) (x : ℝ) :
    (gaussianDerivPoly α p).eval x = p.derivative.eval x - 2 * α * x * p.eval x := by
  simp [gaussianDerivPoly, mul_assoc]

/-- **G1, the closure.**  The derivative of `p(x) exp(-α x²)` is
`q(x) exp(-α x²)` with `q = gaussianDerivPoly α p` — the same statement as
`hasDerivAt_poly_mul_gaussian`, but with the new class member named, which is
what "closed under differentiation" means as a theorem rather than a slogan. -/
theorem hasDerivAt_poly_gaussian (α : ℝ) (p : Polynomial ℝ) (x : ℝ) :
    HasDerivAt (fun t => p.eval t * Real.exp (-(α * t ^ 2)))
      ((gaussianDerivPoly α p).eval x * Real.exp (-(α * x ^ 2))) x := by
  rw [gaussianDerivPoly_eval]
  exact hasDerivAt_poly_mul_gaussian p α x

end GaussianClosure

/-! ## G2 — the `J`-continuity bricks

The transfer bound integrates these three inequalities against the
mollification error.  They are stated for a general normed group — the inner
product structure of the target plays no role at this layer, and taking the
general form keeps the bricks usable for both the velocity and any derived
field.  The surrogate these bricks are applied to is never a Clay candidate;
only the mollified smooth field is, and only once the transfer bound assembled
from them is proved. -/

section JContinuity

variable {E : Type*} [NormedAddCommGroup E]

/-- **G2 (i), the scalar brick.**  For nonnegative `a` and `b`,

`|a³ - b³| ≤ 3 (max a b)² |a - b|`.

The factorisation `a³ - b³ = (a² + ab + b²)(a - b)` does all the work: each of
the three terms of the symmetric factor is at most `(max a b)²` precisely
because `a` and `b` are nonnegative, which is why the hypotheses are not
decorative.  The exponent-three Lipschitz constant on a bounded set is the
entire analytic mechanism of `J`-continuity. -/
theorem abs_cube_sub_cube_le {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    |a ^ 3 - b ^ 3| ≤ 3 * max a b ^ 2 * |a - b| := by
  have hM : 0 ≤ max a b := le_trans ha (le_max_left a b)
  have haM : a ≤ max a b := le_max_left a b
  have hbM : b ≤ max a b := le_max_right a b
  have hfac : a ^ 3 - b ^ 3 = (a ^ 2 + a * b + b ^ 2) * (a - b) := by ring
  have hcoef : |a ^ 2 + a * b + b ^ 2| = a ^ 2 + a * b + b ^ 2 :=
    abs_of_nonneg (by nlinarith [mul_nonneg ha hb, sq_nonneg a, sq_nonneg b])
  have hbound : a ^ 2 + a * b + b ^ 2 ≤ 3 * max a b ^ 2 := by
    nlinarith [mul_le_mul haM haM ha hM, mul_le_mul hbM hbM hb hM,
      mul_le_mul haM hbM hb hM]
  calc |a ^ 3 - b ^ 3| = |a ^ 2 + a * b + b ^ 2| * |a - b| := by rw [hfac, abs_mul]
    _ ≤ 3 * max a b ^ 2 * |a - b| := by
        rw [hcoef]
        exact mul_le_mul_of_nonneg_right hbound (abs_nonneg _)

/-- **G2 (ii), the brick at the integrand of `J`.**

`| ‖u‖³ - ‖v‖³ | ≤ 3 (‖u‖ + ‖v‖)² ‖u - v‖`.

This is the scalar brick at `a = ‖u‖`, `b = ‖v‖`, with the reverse triangle
inequality `| ‖u‖ - ‖v‖ | ≤ ‖u - v‖` converting the scalar increment into the
vector one, and `max ‖u‖ ‖v‖ ≤ ‖u‖ + ‖v‖` — true because norms are
nonnegative — trading the maximum for the sum that the mollification estimates
actually carry.  Integrated against a mollification error, this is the
statement that `J` of the mollified field converges to `J` of the surrogate,
which is the transfer licence. -/
theorem abs_norm_cube_sub_norm_cube_le (u v : E) :
    |‖u‖ ^ 3 - ‖v‖ ^ 3| ≤ 3 * (‖u‖ + ‖v‖) ^ 2 * ‖u - v‖ := by
  have h1 := abs_cube_sub_cube_le (norm_nonneg u) (norm_nonneg v)
  have h2 : |‖u‖ - ‖v‖| ≤ ‖u - v‖ := abs_norm_sub_norm_le u v
  have hM : (0 : ℝ) ≤ max ‖u‖ ‖v‖ := le_trans (norm_nonneg u) (le_max_left _ _)
  have hmax : max ‖u‖ ‖v‖ ≤ ‖u‖ + ‖v‖ :=
    max_le (le_add_of_nonneg_right (norm_nonneg v)) (le_add_of_nonneg_left (norm_nonneg u))
  have hsq : max ‖u‖ ‖v‖ ^ 2 ≤ (‖u‖ + ‖v‖) ^ 2 := by
    nlinarith [mul_self_le_mul_self hM hmax]
  calc |‖u‖ ^ 3 - ‖v‖ ^ 3| ≤ 3 * max ‖u‖ ‖v‖ ^ 2 * |‖u‖ - ‖v‖| := h1
    _ ≤ 3 * (‖u‖ + ‖v‖) ^ 2 * ‖u - v‖ :=
        mul_le_mul (by linarith) h2 (abs_nonneg _) (by positivity)

variable [NormedSpace ℝ E]

/-- **G2 (iii), the vector brick.**  For the nonlinearity `|u| u` that the
mollification bound differences,

`‖ ‖u‖•u - ‖v‖•v ‖ ≤ (‖u‖ + ‖v‖) ‖u - v‖`.

The proof is the telescoping `‖u‖•u - ‖v‖•v = ‖u‖•(u - v) + (‖u‖ - ‖v‖)•v`:
the first term is bounded by `‖u‖ ‖u - v‖` because scalar multiplication is
norm-multiplicative, the second by `‖u - v‖ ‖v‖` through the reverse triangle
inequality, and the sum is the claimed constant.  The sharper-looking constant
`2 max ‖u‖ ‖v‖` bounds this one from above, so nothing downstream is lost by
proving the sum form, and the sum is what the telescoping actually yields. -/
theorem norm_smul_sub_smul_le (u v : E) :
    ‖‖u‖ • u - ‖v‖ • v‖ ≤ (‖u‖ + ‖v‖) * ‖u - v‖ := by
  have hsplit : ‖u‖ • u - ‖v‖ • v = ‖u‖ • (u - v) + (‖u‖ - ‖v‖) • v := by
    rw [smul_sub, sub_smul]
    abel
  rw [hsplit]
  calc ‖‖u‖ • (u - v) + (‖u‖ - ‖v‖) • v‖
      ≤ ‖‖u‖ • (u - v)‖ + ‖(‖u‖ - ‖v‖) • v‖ := norm_add_le _ _
    _ = ‖u‖ * ‖u - v‖ + |‖u‖ - ‖v‖| * ‖v‖ := by
        rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs,
          abs_of_nonneg (norm_nonneg u)]
    _ ≤ ‖u‖ * ‖u - v‖ + ‖u - v‖ * ‖v‖ :=
        add_le_add le_rfl
          (mul_le_mul_of_nonneg_right (abs_norm_sub_norm_le u v) (norm_nonneg v))
    _ = (‖u‖ + ‖v‖) * ‖u - v‖ := by ring

end JContinuity

/-! ## G3 — the torus control bound -/

section TorusControl

/-- **G3.  The Track-P torus control bound.**

`riccati_comparison` specialised to the coefficient shape the torus estimate
produces.  On the torus the Poincaré inequality gives the viscous term a
genuine linear damping `-ν R` — the term whose absence on `ℝ³` is documented
at length in `ControlODE.lean` — and the estimate layer contributes a linear
growth `a₀ R`, so the net linear coefficient is `-ν + a₀`, of no fixed sign.
The quadratic coefficient `b ≥ 0` is the advection term, the constant `e ≥ 0`
is the residual, and `R 0 = 0` because the error of the computed trajectory
vanishes at the initial time.  A subsolution of the resulting Riccati
inequality never overtakes a supersolution that starts nonnegative.

The hypotheses `_hν` and `_he` are not consumed by the comparison — the
underlying theorem is indifferent to both signs — but they pin the Track-P
shape (fixed positive viscosity, nonnegative residual) so that citing this
theorem asserts the intended regime and not an accidental generalisation.
The proof is one line by design: the point of this corollary is the citable
name, so that the Python payload references a theorem whose statement already
contains the substitution instead of a recipe for performing it. -/
theorem torus_control_bound {R R' S S' : ℝ → ℝ} {ν a₀ b e T : ℝ}
    (_hν : 0 < ν) (hb : 0 ≤ b) (_he : 0 ≤ e)
    (hR : ContinuousOn R (Icc 0 T)) (hS : ContinuousOn S (Icc 0 T))
    (hR' : ∀ t ∈ Ico (0 : ℝ) T, HasDerivWithinAt R (R' t) (Ici t) t)
    (hS' : ∀ t ∈ Ico (0 : ℝ) T, HasDerivWithinAt S (S' t) (Ici t) t)
    (hsub : ∀ t ∈ Ico (0 : ℝ) T, R' t ≤ (-ν + a₀) * R t + b * R t ^ 2 + e)
    (hsuper : ∀ t ∈ Ico (0 : ℝ) T, (-ν + a₀) * S t + b * S t ^ 2 + e ≤ S' t)
    (h₀ : R 0 = 0) (hS₀ : 0 ≤ S 0) :
    ∀ t ∈ Icc (0 : ℝ) T, R t ≤ S t :=
  riccati_comparison hb hR hS hR' hS' hsub hsuper (h₀.le.trans hS₀)

end TorusControl

end

end NSSingularity
