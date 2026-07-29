/-
The algebraic core of the `L^3` generation identity, and of the pure-swirl no-go
that rules out the most obvious family of candidate initial data.

## The mathematics this file sits under

For `u` solving incompressible Navier-Stokes on `R^3` with the decaying Riesz
pressure `p = R_i R_j (u_i u_j)`, and with `|u| + |grad u| <= C <x>^{-4}` locally
uniformly in time, the quantity `F(t) = int |u|^3` obeys the classical `L^p`
energy identity at `p = 3`,

  `F'(t) = 3 int p div(|u| u) - 3 nu int ( |u| |grad u|^2 + |u| |grad |u||^2 )`.

Three features of that identity decide everything downstream, and all three are
*pointwise algebra* dressed up as analysis:

* the transport contribution vanishes identically for every divergence-free
  field, because `|u| u . (u.grad) u = (1/3) u . grad(|u|^3)` is a pure
  divergence -- the cubic chain rule and nothing else;
* the viscous contribution is `-3 nu` times an integral of two manifestly
  nonnegative densities, hence never positive;
* therefore a positive `F'(0)` can only be produced by the pressure term.

The pure-swirl no-go then says the pressure term is *also* zero for the natural
first guess at initial data.  If `u_0 = u^theta(r,z) e_theta` with
`u^theta = r g(r^2,z)` (the smoothness-at-the-axis form), then the Cartesian
field is `g(x^2+y^2, z) . (-y, x, 0)`; the pressure inherited from a rotationally
equivariant map is axisymmetric; and a pure swirl is pointwise orthogonal to the
gradient of *any* axisymmetric function.  So `F'(0) <= 0`, with equality only for
`u_0 = 0`.

## What is formalised here

Everything below is a genuine theorem about real numbers, about curves in a real
inner product space, or about scalar fields on `R^3` written as `R x R x R`.
Nothing below is an axiom, and nothing below is weakened to make it cheap.

* `hasDerivAt_norm_cube` -- the chain rule at exponent three,
  `d/dt ||f t||^3 = 3 ||f t|| <f t, f'>`, valid *including at the zeros of* `f`,
  which is exactly the point: `||.||` itself is not differentiable at the origin
  but its cube is, and it is the cube that appears in `F`.
* `hasDerivAt_cube`, `transport_eq_one_third_deriv` -- the transport
  cancellation.  The second is the statement worth citing: the transport
  integrand `|u| <u, (u.grad)u>` is *exactly one third* of a total derivative
  along the flow, so it integrates to zero against a divergence-free field.
* `regSpeed`, `eps_le_regSpeed`, `norm_le_regSpeed`, `regSpeed_le_norm_add`,
  `contDiff_regSpeed`, `hasDerivAt_regSpeed` -- the regularised speed
  `s_eps = sqrt(|v|^2 + eps^2)`.  This is the device that lets the interval
  certificate differentiate the speed without ever dividing by `|u|`; the
  sandwich `|v| <= s_eps <= |v| + eps` is what makes the regularisation error
  linear in `eps`.
* `viscous_contribution_nonpos`, `kato_split_le`,
  `positive_generation_forces_pressure` -- the sign structure.  `kato_split_le`
  is the sharp form the certificate actually uses: given the Kato bound
  `|G| <= q g` with `q = |u|`, `g = |grad u|` and `G = |u| grad |u|`, the split
  density `q g^2 + |G|^2/q` never exceeds `2 q g^2`, with the quotient defined to
  be `0` at `q = 0` so that no division by the speed survives into the checkable
  arithmetic.
* `generation_pos_iff_reynolds_gt` -- the scaling criterion.  Under
  `u -> A u(x/L)` the pressure term scales as `A^4 L^2` and the viscous term as
  `nu A^3 L`, so generation is positive precisely when the Reynolds number
  `A L / nu` exceeds the shape-dependent threshold `V_1 / P_1`.
* `swirl_cartesianDiv_eq_zero` -- a pure-swirl field is divergence free.
* `swirl_dot_grad_axisymmetric_eq_zero` -- **the heart of the no-go**: a pure
  swirl is pointwise orthogonal to the gradient of any axisymmetric scalar.  This
  is the reason the pressure term dies; it is an identity, not an estimate.
* `eventually_eq_of_fderiv_eq_zero_on`,
  `eq_of_locallyConstant_of_ne_zero`,
  `eq_zero_of_locallyConstant_of_tendsto_cocompact`,
  `pure_swirl_equality_case` -- the equality case.  A field whose derivative
  vanishes wherever the field itself does not is locally constant there; the set
  where it equals that constant is then open and closed; connectedness of `R^3`
  makes it everything; and decay at infinity forces the constant to be zero.

## What is *not* formalised, and why

* **The integration by parts.**  Passing from the pointwise identities above to
  `F'(t) = ...` needs a divergence theorem on `R^3` for a field with only
  polynomial decay, together with differentiation under the integral sign.
  Mathlib has neither in a form that applies to this hypothesis class, and
  assuming them as hypotheses would move the entire content of the identity into
  an unproved premise.  So the file stops at the integrands.
* **The Riesz pressure.**  `p = R_i R_j (u_i u_j)` and, above all, its
  equivariance under rotations about `e_z` -- which is what makes `p`
  axisymmetric for axisymmetric data -- would need the Riesz transform on
  `R^3` and its commutation with the orthogonal group.  Neither is available.
  `swirl_dot_grad_axisymmetric_eq_zero` therefore takes the axisymmetry of the
  scalar as a hypothesis rather than deriving it: that is the honest division,
  and it is the analytic half that is missing, not the algebraic half.
* **The decay class and its propagation.**  `|u| + |grad u| <= C <x>^{-4}`
  locally uniformly in time is a hypothesis on the solution, not a conclusion;
  the Schwartz class is *not* propagated (the Brandolese localisation
  obstruction), so no shortcut through Schwartz data is available.  In this file
  decay enters only through `Tendsto u (cocompact _) (nhds 0)`, which is all the
  equality case needs.
* **The parity selection rule.**  That the pressure integrand has `z`-parity
  `-sigma_psi`, so that it vanishes on a symmetric domain unless `psi_1` is odd
  in `z`, is a statement about an integral over a symmetric domain and belongs
  to the missing measure-theoretic layer.  The pointwise parity bookkeeping on
  its own proves nothing.
* **`||u||_{Hdot^n}` norm transfer.**  The exact isometry between the
  three-dimensional Sobolev norms and the `r^3`-weighted cylindrical ones is a
  measure-theoretic statement about `int_{r>0} f g r^3 dr dz`; it is out of
  scope for the same reason as the integration by parts.

Nothing in this file, alone or in combination with the rest of the library, is
a statement about the Clay problem.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib

namespace NSSingularity

noncomputable section

open Set Filter Topology
open scoped RealInnerProductSpace

/-! ## The chain rule at exponent three, and the transport cancellation

The `L^3` functional differentiates through the norm cube.  The only delicate
point is that `‖·‖` is not differentiable at the origin while `‖·‖ ^ 3` is, so
the identity must not be obtained by writing `‖v‖ ^ 3 = ‖v‖ ^ 2 * ‖v‖` and
differentiating the second factor.  Mathlib's `hasFDerivAt_norm_rpow` handles the
origin correctly for every real exponent `p > 1`, and `p = 3` is the case wanted
here. -/

section NormCube

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- **T1.  The chain rule at exponent three.**

`d/dt ‖f t‖ ^ 3 = 3 ‖f t‖ ⟪f t, f' t⟫`, with no hypothesis excluding the zeros
of `f`.  This is the derivative that produces the whole `L^3` generation
identity: every term of `F'(t)` is this coefficient contracted against one of
the three terms of the momentum equation. -/
theorem hasDerivAt_norm_cube {f : ℝ → E} {f' : E} {t : ℝ} (hf : HasDerivAt f f' t) :
    HasDerivAt (fun x => ‖f x‖ ^ 3) (3 * ‖f t‖ * ⟪f t, f'⟫) t := by
  have hrpow : HasFDerivAt (fun v : E => ‖v‖ ^ (3 : ℝ))
      ((3 * ‖f t‖ ^ ((3 : ℝ) - 2)) • innerSL ℝ (f t)) (f t) :=
    hasFDerivAt_norm_rpow (f t) (by norm_num)
  have hcomp : HasDerivAt (fun x => ‖f x‖ ^ (3 : ℝ))
      (((3 * ‖f t‖ ^ ((3 : ℝ) - 2)) • innerSL ℝ (f t)) f') t :=
    hrpow.comp_hasDerivAt t hf
  have hpt : ∀ x : ℝ, ‖f x‖ ^ (3 : ℝ) = ‖f x‖ ^ (3 : ℕ) := fun x => by
    rw [show ((3 : ℝ)) = ((3 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
  have hval : ((3 * ‖f t‖ ^ ((3 : ℝ) - 2)) • innerSL ℝ (f t)) f'
      = 3 * ‖f t‖ * ⟪f t, f'⟫ := by
    rw [show ((3 : ℝ) - 2) = 1 by norm_num, Real.rpow_one]
    simp [mul_assoc]
  simpa only [hpt, hval] using hcomp

/-- **T3, the scalar chain rule.**  `d/dt (f t) ^ 3 = 3 (f t) ^ 2 f'`.

Stated separately because the transport cancellation is used along a *scalar*
transported quantity as often as along the velocity itself. -/
theorem hasDerivAt_cube {f : ℝ → ℝ} {f' t : ℝ} (hf : HasDerivAt f f' t) :
    HasDerivAt (fun x => f x ^ 3) (3 * f t ^ 2 * f') t := by
  have h := hf.pow 3
  norm_num at h
  exact h

/-- **T3.  The transport cancellation.**

The transport integrand `‖u‖ ⟪u, (u . grad) u⟫` is *exactly one third* of the
derivative of `‖u‖ ^ 3` along the flow.  That is the entire reason the transport
term of the `L^3` identity vanishes: `(1/3) u . grad (|u|^3)` is a pure
divergence for a divergence-free `u`, so it integrates to zero, and the
cancellation is exact rather than an estimate.

Here `f` is the value of the field along the flow line and `f'` its derivative
along it, so `⟪f t, f'⟫` is the pointwise `u . (u . grad) u`. -/
theorem transport_eq_one_third_deriv {f : ℝ → E} {f' : E} {t : ℝ}
    (hf : HasDerivAt f f' t) :
    ‖f t‖ * ⟪f t, f'⟫ = (1 / 3) * deriv (fun x => ‖f x‖ ^ 3) t := by
  rw [(hasDerivAt_norm_cube hf).deriv]
  ring

end NormCube

/-! ## The regularised speed

The certificate never differentiates `‖u‖`, because that fails at the zeros of
`u` and, worse, is numerically hopeless near them.  It differentiates
`s_eps = sqrt(‖u‖ ^ 2 + eps ^ 2)` instead, which is smooth everywhere, bounded
below by `eps`, and within `eps` of the true speed.  The last two facts are what
turn the regularisation into a linear-in-`eps` perturbation of the identity
rather than an uncontrolled change of problem. -/

section RegularisedSpeed

variable {E : Type*} [NormedAddCommGroup E]

/-- The regularised speed `s_eps(v) = sqrt(‖v‖ ^ 2 + eps ^ 2)`. -/
def regSpeed (ε : ℝ) (v : E) : ℝ := Real.sqrt (‖v‖ ^ 2 + ε ^ 2)

/-- The regularised speed is bounded below by the regularisation parameter, so
dividing by it is always legitimate. -/
theorem eps_le_regSpeed {ε : ℝ} (hε : 0 ≤ ε) (v : E) : ε ≤ regSpeed ε v := by
  have h : ε ^ 2 ≤ ‖v‖ ^ 2 + ε ^ 2 := by nlinarith [sq_nonneg ‖v‖]
  calc ε = Real.sqrt (ε ^ 2) := (Real.sqrt_sq hε).symm
    _ ≤ Real.sqrt (‖v‖ ^ 2 + ε ^ 2) := Real.sqrt_le_sqrt h

/-- The regularised speed never underestimates the true speed. -/
theorem norm_le_regSpeed (ε : ℝ) (v : E) : ‖v‖ ≤ regSpeed ε v := by
  have h : ‖v‖ ^ 2 ≤ ‖v‖ ^ 2 + ε ^ 2 := by nlinarith [sq_nonneg ε]
  calc ‖v‖ = Real.sqrt (‖v‖ ^ 2) := (Real.sqrt_sq (norm_nonneg v)).symm
    _ ≤ Real.sqrt (‖v‖ ^ 2 + ε ^ 2) := Real.sqrt_le_sqrt h

/-- The regularised speed overestimates the true speed by at most `eps`.  With
`norm_le_regSpeed` this is the sandwich `‖v‖ ≤ s_eps ≤ ‖v‖ + eps`, which makes
every regularisation error linear in `eps`. -/
theorem regSpeed_le_norm_add {ε : ℝ} (hε : 0 ≤ ε) (v : E) :
    regSpeed ε v ≤ ‖v‖ + ε := by
  have h : ‖v‖ ^ 2 + ε ^ 2 ≤ (‖v‖ + ε) ^ 2 := by nlinarith [norm_nonneg v]
  calc regSpeed ε v ≤ Real.sqrt ((‖v‖ + ε) ^ 2) := Real.sqrt_le_sqrt h
    _ = ‖v‖ + ε := Real.sqrt_sq (by positivity)

variable [InnerProductSpace ℝ E]

/-- **T2.**  The regularised speed is `C^∞` on all of `E`.  Unregularised, `‖·‖`
is not even once differentiable at the origin, which is precisely the reason the
regularisation exists. -/
theorem contDiff_regSpeed {ε : ℝ} (hε : 0 < ε) :
    ContDiff ℝ (⊤ : ℕ∞) (fun v : E => regSpeed ε v) := by
  have hne : ∀ v : E, ‖v‖ ^ 2 + ε ^ 2 ≠ 0 := fun v => by
    have := pow_pos hε 2
    nlinarith [sq_nonneg ‖v‖]
  exact ((contDiff_norm_sq ℝ).add contDiff_const).sqrt hne

/-- **T2, the derivative.**  Along a curve, `d/dt s_eps(f t) = ⟪f t, f'⟫ /
s_eps(f t)`, with no side condition: the denominator is at least `eps > 0`.

This is the identity the interval certificate propagates in place of the
unregularised speed equation. -/
theorem hasDerivAt_regSpeed {ε : ℝ} (hε : 0 < ε) {f : ℝ → E} {f' : E} {t : ℝ}
    (hf : HasDerivAt f f' t) :
    HasDerivAt (fun x => regSpeed ε (f x)) (⟪f t, f'⟫ / regSpeed ε (f t)) t := by
  have hpos : (0 : ℝ) < ‖f t‖ ^ 2 + ε ^ 2 := by
    have := pow_pos hε 2
    nlinarith [sq_nonneg ‖f t‖]
  have hsq : HasDerivAt (fun x => ‖f x‖ ^ 2 + ε ^ 2) (2 * ⟪f t, f'⟫) t := by
    simpa using hf.norm_sq.add_const (ε ^ 2)
  have hroot := hsq.sqrt (ne_of_gt hpos)
  have hsqrt_pos : (0 : ℝ) < Real.sqrt (‖f t‖ ^ 2 + ε ^ 2) := Real.sqrt_pos.mpr hpos
  simp only [regSpeed]
  convert hroot using 1
  field_simp

end RegularisedSpeed

/-! ## Signs, the Kato split, and the scaling criterion

The viscous term of the `L^3` identity is `-3 nu` times an integral of
nonnegative densities, so it can only push `F` down.  Consequently a positive
generation rate is *by elimination* a statement about the pressure term, and
that is the logical step `positive_generation_forces_pressure` records.

`kato_split_le` is the form the interval certificate consumes.  The two viscous
densities are `|u| |grad u|^2` and `|u| |grad |u||^2`, and the second is the
awkward one: `grad |u| = (u . grad u)/|u|` is not evaluable at the zeros of `u`.
The certificate therefore never forms `grad |u|`.  It carries
`G = |u| grad |u| = (1/2) grad(|u|^2)`, which is a polynomial in the state, and
recovers the density as `|G|^2 / q` with `q = |u|`, declared to be `0` when
`q = 0`.  With `g = |grad u|` the Kato inequality `|grad |u|| <= |grad u|` is
exactly `|G| <= q g`, and the lemma then says the assembled viscous density
`q g^2 + |G|^2/q` never exceeds `2 q g^2`. -/

section Signs

/-- **T4.**  The viscous contribution to `F'` is never positive.  This pins the
sign convention: `a` stands for the speed, `b` and `c` for the two nonnegative
gradient densities. -/
theorem viscous_contribution_nonpos {ν a b c : ℝ} (hν : 0 < ν)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) :
    -3 * ν * (a * b + c) ≤ 0 := by
  nlinarith [mul_nonneg ha hb]

/-- **T4, sharp.**  The Kato split without division by the speed.

If `q ≥ 0` and `‖G‖ ≤ q g`, then `q g ^ 2 + ‖G‖ ^ 2 / q ≤ 2 q g ^ 2`, where the
quotient is *defined* to be `0` at `q = 0`.  Read `q = |u|`, `g = |grad u|` and
`G = |u| grad |u|`; the hypothesis is then the Kato inequality and the left-hand
side is the exact viscous density of the `L^3` identity.  The point is not the
factor two but that the left-hand side is a total function of `(q, g, G)`: the
certificate can evaluate it on an interval straddling `q = 0` without a case
split and without a removable singularity. -/
theorem kato_split_le {F : Type*} [NormedAddCommGroup F] {q g : ℝ} {G : F}
    (hq : 0 ≤ q) (hG : ‖G‖ ≤ q * g) :
    q * g ^ 2 + (if 0 < q then ‖G‖ ^ 2 / q else 0) ≤ 2 * (q * g ^ 2) := by
  rcases eq_or_lt_of_le hq with h | h
  · rw [← h]
    norm_num
  · rw [if_pos h]
    have hGn : (0 : ℝ) ≤ ‖G‖ := norm_nonneg G
    have hsq : ‖G‖ ^ 2 ≤ (q * g) ^ 2 := by nlinarith
    have hdiv : ‖G‖ ^ 2 / q ≤ q * g ^ 2 := by
      rw [div_le_iff₀ h]
      nlinarith
    linarith

/-- **The elimination step.**  If the generation rate splits as pressure plus
viscosity and the viscous part is nonpositive, then a positive generation rate
forces a positive pressure part.

Trivial as arithmetic, and the whole strategic content of the `L^3` identity:
no divergence-free field generates `L^3` norm through transport, and none
generates it through viscosity, so a candidate has to earn its growth from the
pressure. -/
theorem positive_generation_forces_pressure {Fdot P V : ℝ}
    (hsplit : Fdot = P + V) (hV : V ≤ 0) (hpos : 0 < Fdot) : 0 < P := by
  linarith [hsplit ▸ hpos]

/-- **The scaling criterion.**

Under `u -> A u(x/L)` the pressure contribution scales as `A ^ 4 L ^ 2` and the
viscous contribution as `nu A ^ 3 L`, with shape constants `P₁ > 0` and `V₁`.
Generation is positive exactly when the Reynolds number `A L / nu` exceeds the
shape-dependent threshold `V₁ / P₁`, which is the definition of `Re_crit`.

The physical amplitude `A` here and the swirl-variable amplitude `A₁` of the
repository's parametrisation give the same number, since `A = A₁ L`. -/
theorem generation_pos_iff_reynolds_gt {A L ν P₁ V₁ : ℝ}
    (hA : 0 < A) (hL : 0 < L) (hν : 0 < ν) (hP : 0 < P₁) :
    0 < A ^ 4 * L ^ 2 * P₁ - ν * A ^ 3 * L * V₁ ↔ V₁ / P₁ < A * L / ν := by
  have hfac : A ^ 4 * L ^ 2 * P₁ - ν * A ^ 3 * L * V₁
      = A ^ 3 * L * (A * L * P₁ - ν * V₁) := by ring
  have hpos : (0 : ℝ) < A ^ 3 * L := by positivity
  rw [hfac, mul_pos_iff_of_pos_left hpos, div_lt_div_iff₀ hP hν]
  constructor <;> intro h <;> linarith

end Signs

/-! ## Pure swirl

A swirl field smooth across the axis has `u^theta = r g(r^2, z)`, so in Cartesian
coordinates, using `e_theta = (-y, x, 0) / r`,

  `u(x,y,z) = g(x^2 + y^2, z) . (-y, x, 0)`.

Two identities follow, and between them they are the no-go.  The first is that
such a field is divergence free, so it is admissible initial data.  The second is
that it is pointwise orthogonal to the gradient of any function of `(r^2, z)` --
in particular to the gradient of the pressure, which is axisymmetric because the
map `u -> R_i R_j (u_i u_j)` commutes with rotations about `e_z`.  The pressure
term of the `L^3` identity is therefore identically zero, not merely small, and
by `positive_generation_forces_pressure` no pure swirl generates `L^3` norm.

Only the algebra is here.  That the Riesz pressure of an axisymmetric field is
axisymmetric is the missing analytic half; `swirl_dot_grad_axisymmetric_eq_zero`
takes it as a hypothesis on the scalar. -/

section PureSwirl

/-- Partial derivative in the first Cartesian slot. -/
def slotX (V : ℝ × ℝ × ℝ → ℝ) (p : ℝ × ℝ × ℝ) : ℝ :=
  deriv (fun t => V (t, p.2.1, p.2.2)) p.1

/-- Partial derivative in the second Cartesian slot. -/
def slotY (V : ℝ × ℝ × ℝ → ℝ) (p : ℝ × ℝ × ℝ) : ℝ :=
  deriv (fun t => V (p.1, t, p.2.2)) p.2.1

/-- Partial derivative in the third Cartesian slot. -/
def slotZ (V : ℝ × ℝ × ℝ → ℝ) (p : ℝ × ℝ × ℝ) : ℝ :=
  deriv (fun t => V (p.1, p.2.1, t)) p.2.2

/-- First component of the pure-swirl field `g(r^2, z) . (-y, x, 0)`. -/
def swirlX (g : ℝ → ℝ → ℝ) (p : ℝ × ℝ × ℝ) : ℝ :=
  -p.2.1 * g (p.1 ^ 2 + p.2.1 ^ 2) p.2.2

/-- Second component of the pure-swirl field. -/
def swirlY (g : ℝ → ℝ → ℝ) (p : ℝ × ℝ × ℝ) : ℝ :=
  p.1 * g (p.1 ^ 2 + p.2.1 ^ 2) p.2.2

/-- Third component of the pure-swirl field.  A swirl has no axial part. -/
def swirlZ (_g : ℝ → ℝ → ℝ) (_p : ℝ × ℝ × ℝ) : ℝ := 0

/-- An axisymmetric scalar: a function of `r ^ 2 = x ^ 2 + y ^ 2` and `z`.  Using
`r ^ 2` rather than `r` is not a convenience -- it is the condition that makes the
field smooth across the axis. -/
def axisym (φ : ℝ → ℝ → ℝ) (p : ℝ × ℝ × ℝ) : ℝ :=
  φ (p.1 ^ 2 + p.2.1 ^ 2) p.2.2

/-- The Cartesian divergence of a field given by its three components. -/
def cartesianDiv (Vx Vy Vz : ℝ × ℝ × ℝ → ℝ) (p : ℝ × ℝ × ℝ) : ℝ :=
  slotX Vx p + slotY Vy p + slotZ Vz p

/-- The radial chain rule used twice below: `d/dx (x ^ 2 + y ^ 2) = 2x`. -/
private theorem hasDerivAt_radiusSq_fst (x y : ℝ) :
    HasDerivAt (fun t : ℝ => t ^ 2 + y ^ 2) (2 * x) x := by
  simpa using (hasDerivAt_pow 2 x).add_const (y ^ 2)

/-- The radial chain rule in the second slot: `d/dy (x ^ 2 + y ^ 2) = 2y`. -/
private theorem hasDerivAt_radiusSq_snd (x y : ℝ) :
    HasDerivAt (fun t : ℝ => x ^ 2 + t ^ 2) (2 * y) y := by
  simpa using (hasDerivAt_pow 2 y).const_add (x ^ 2)

/-- **T5.  A pure swirl is divergence free.**

`d_x(-y g) = -y . 2x . g_s` and `d_y(x g) = x . 2y . g_s` cancel, and the axial
component is identically zero.  Note that no regularity in `z` is used and no
condition on `g` beyond differentiability in its first argument: the
cancellation is structural. -/
theorem swirl_cartesianDiv_eq_zero (g dg : ℝ → ℝ → ℝ)
    (hg : ∀ s z : ℝ, HasDerivAt (fun σ => g σ z) (dg s z) s) (p : ℝ × ℝ × ℝ) :
    cartesianDiv (swirlX g) (swirlY g) (swirlZ g) p = 0 := by
  obtain ⟨x, y, z⟩ := p
  have hgx := (hg (x ^ 2 + y ^ 2) z).comp x (hasDerivAt_radiusSq_fst x y)
  have hgy := (hg (x ^ 2 + y ^ 2) z).comp y (hasDerivAt_radiusSq_snd x y)
  have hX : HasDerivAt (fun t : ℝ => -y * g (t ^ 2 + y ^ 2) z)
      (-y * (dg (x ^ 2 + y ^ 2) z * (2 * x))) x := hgx.const_mul (-y)
  have hY : HasDerivAt (fun t : ℝ => x * g (x ^ 2 + t ^ 2) z)
      (x * (dg (x ^ 2 + y ^ 2) z * (2 * y))) y := hgy.const_mul x
  have hZ : HasDerivAt (fun _ : ℝ => (0 : ℝ)) 0 z := hasDerivAt_const z 0
  simp only [cartesianDiv, slotX, slotY, slotZ, swirlX, swirlY, swirlZ]
  rw [hX.deriv, hY.deriv, hZ.deriv]
  ring

/-- **T5, the heart of the no-go.**  A pure swirl is pointwise orthogonal to the
gradient of every axisymmetric scalar.

`grad phi(r^2, z) = (2x phi_s, 2y phi_s, phi_z)` and the swirl is
`g . (-y, x, 0)`, so the inner product is `g (-2xy + 2xy) phi_s = 0`.  The axial
derivative of `phi` never enters, because the swirl has no axial component --
which is why no regularity of `phi` in `z` is assumed.

Applied with `phi` the pressure of an axisymmetric field, this says
`u_0 . grad p = 0` pointwise, so the pressure term of the `L^3` identity is zero
and, by `positive_generation_forces_pressure`, no pure swirl generates. -/
theorem swirl_dot_grad_axisymmetric_eq_zero (g φ dφ : ℝ → ℝ → ℝ)
    (hφ : ∀ s z : ℝ, HasDerivAt (fun σ => φ σ z) (dφ s z) s) (p : ℝ × ℝ × ℝ) :
    swirlX g p * slotX (axisym φ) p
      + swirlY g p * slotY (axisym φ) p
      + swirlZ g p * slotZ (axisym φ) p = 0 := by
  obtain ⟨x, y, z⟩ := p
  have hXc := (hφ (x ^ 2 + y ^ 2) z).comp x (hasDerivAt_radiusSq_fst x y)
  have hYc := (hφ (x ^ 2 + y ^ 2) z).comp y (hasDerivAt_radiusSq_snd x y)
  have hX : HasDerivAt (fun t : ℝ => φ (t ^ 2 + y ^ 2) z)
      (dφ (x ^ 2 + y ^ 2) z * (2 * x)) x := hXc
  have hY : HasDerivAt (fun t : ℝ => φ (x ^ 2 + t ^ 2) z)
      (dφ (x ^ 2 + y ^ 2) z * (2 * y)) y := hYc
  simp only [slotX, slotY, slotZ, swirlX, swirlY, swirlZ, axisym]
  rw [hX.deriv, hY.deriv]
  ring

end PureSwirl

/-! ## The equality case

The no-go gives `F'(0) ≤ 0`, and the equality case has to rule out a nonzero
field achieving `F'(0) = 0` for trivial reasons.  The argument is:
`grad u_0 = 0` on the open set where `u_0 ≠ 0`, so `u_0` is locally constant
there; the set on which `u_0` equals a fixed nonzero value is then both open --
by local constancy -- and closed -- by continuity -- hence all of `R^3` by
connectedness; and a nonzero constant does not tend to zero at infinity.

Each of those four steps is a separate theorem below.  The first is a bridge
from `fderiv = 0` to local constancy, so that a caller who has the derivative
condition never has to produce the topological hypothesis by hand. -/

section EqualityCase

/-- A function whose Fréchet derivative vanishes on an open set is locally
constant at each point of that set.  The mathlib input is
`IsOpen.isOpen_inter_preimage_of_fderiv_eq_zero`, which says that the trace of a
level set on the open set is open. -/
theorem eventually_eq_of_fderiv_eq_zero_on
    {E G : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [NormedAddCommGroup G] [NormedSpace ℝ G]
    {U : Set E} (hU : IsOpen U) {u : E → G} (hdiff : DifferentiableOn ℝ u U)
    (hzero : U.EqOn (fderiv ℝ u) 0) {x : E} (hx : x ∈ U) :
    ∀ᶠ y in 𝓝 x, u y = u x := by
  have hopen := hU.isOpen_inter_preimage_of_fderiv_eq_zero hdiff hzero {u x}
  have hmem : x ∈ U ∩ u ⁻¹' {u x} := ⟨hx, rfl⟩
  filter_upwards [hopen.mem_nhds hmem] with y hy
  exact hy.2

/-- **T6, the clopen step.**  On a preconnected space, a continuous function that
is locally constant wherever it is nonzero, and is nonzero somewhere, is constant.

The set where the function takes the value `u x₀ ≠ 0` is closed by continuity and
open by local constancy -- local constancy is available at each of its points
precisely because the value there is nonzero -- so preconnectedness makes it
everything. -/
theorem eq_of_locallyConstant_of_ne_zero
    {X G : Type*} [TopologicalSpace X] [PreconnectedSpace X]
    [TopologicalSpace G] [Zero G] [T1Space G]
    {u : X → G} (hu : Continuous u)
    (hloc : ∀ x, u x ≠ 0 → ∀ᶠ y in 𝓝 x, u y = u x)
    {x₀ : X} (hx₀ : u x₀ ≠ 0) (x : X) : u x = u x₀ := by
  have hclosed : IsClosed {x : X | u x = u x₀} :=
    IsClosed.preimage hu isClosed_singleton
  have hopen : IsOpen {x : X | u x = u x₀} := by
    rw [isOpen_iff_mem_nhds]
    intro y hy
    have hne : u y ≠ 0 := by
      rw [show u y = u x₀ from hy]
      exact hx₀
    filter_upwards [hloc y hne] with w hw
    show u w = u x₀
    rw [hw]
    exact hy
  have hclopen : IsClopen {x : X | u x = u x₀} := ⟨hclosed, hopen⟩
  have huniv : {x : X | u x = u x₀} = univ := hclopen.eq_univ ⟨x₀, rfl⟩
  have hx : x ∈ {x : X | u x = u x₀} := huniv ▸ mem_univ x
  exact hx

/-- **T6, with decay.**  Adding decay at infinity to the previous theorem forces
the function to vanish identically.

Decay is expressed as convergence to zero along the cocompact filter, which is
the weakest form of `|u| -> 0 as |x| -> infinity` that makes sense on a general
space, and is implied by any of the pointwise decay rates used elsewhere in the
project.  `NoncompactSpace` is what makes the cocompact filter nontrivial; on a
compact space the statement is false and should be. -/
theorem eq_zero_of_locallyConstant_of_tendsto_cocompact
    {X G : Type*} [TopologicalSpace X] [PreconnectedSpace X] [NoncompactSpace X]
    [TopologicalSpace G] [Zero G] [T1Space G]
    {u : X → G} (hu : Continuous u)
    (hloc : ∀ x, u x ≠ 0 → ∀ᶠ y in 𝓝 x, u y = u x)
    (hdecay : Tendsto u (cocompact X) (𝓝 0)) (x : X) : u x = 0 := by
  by_contra hx
  have hconst : u = fun _ => u x := by
    funext y
    exact eq_of_locallyConstant_of_ne_zero hu hloc hx y
  rw [hconst] at hdecay
  exact hx (tendsto_const_nhds_iff.mp hdecay)

/-- **T6, the equality case of the pure-swirl no-go.**

A differentiable field on `R^3` whose derivative vanishes wherever the field does
not, and which decays at infinity, is identically zero.  This is the step that
upgrades `F'(0) ≤ 0` to `F'(0) < 0 unless u_0 = 0`: vanishing of the viscous
density forces `grad u_0 = 0` on `{u_0 ≠ 0}`, and the only decaying field with
that property is the zero field. -/
theorem pure_swirl_equality_case
    {u : EuclideanSpace ℝ (Fin 3) → EuclideanSpace ℝ (Fin 3)}
    (hu : Differentiable ℝ u)
    (hzero : ∀ x, u x ≠ 0 → fderiv ℝ u x = 0)
    (hdecay : Tendsto u (cocompact (EuclideanSpace ℝ (Fin 3))) (𝓝 0)) :
    ∀ x, u x = 0 := by
  set U : Set (EuclideanSpace ℝ (Fin 3)) := {x | u x ≠ 0} with hUdef
  have hUopen : IsOpen U := by
    rw [hUdef]
    exact isOpen_ne_fun hu.continuous continuous_const
  have hEqOn : U.EqOn (fderiv ℝ u) 0 := fun y hy => hzero y hy
  have hloc : ∀ x, u x ≠ 0 → ∀ᶠ y in 𝓝 x, u y = u x := fun x hx =>
    eventually_eq_of_fderiv_eq_zero_on hUopen hu.differentiableOn hEqOn hx
  exact fun x =>
    eq_zero_of_locallyConstant_of_tendsto_cocompact hu.continuous hloc hdecay x

end EqualityCase

end

end NSSingularity
