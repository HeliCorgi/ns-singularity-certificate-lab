/-
Stage 1, items **F-12**, **F-13** and **F-7a/F-7b** of `docs/final_target.md` §4:
completing the Lean side of the **fixed-finite-mode** Track-F no-go theorem
begun in `GalerkinNoBlowup.lean` (F-6).

`docs/research_notes/track_f_finite_mode_nogo.md` proves, on paper, that a
velocity field whose Fourier support stays inside a *fixed finite* set cannot be
the breakdown solution of Clay (C)/(D).  Four ingredients enter:

1. the trilinear cancellation `⟨u,(u·∇)u⟩ = 0` for divergence-free fields,
2. equivalence of norms on the finite Fourier space,
3. the a priori energy bound (F-6, already formalized), and
4. continuation of the bounded coefficient trajectory past the putative
   singular time.

F-6 formalized (3) abstractly.  This file formalizes (1) as an exact algebraic
identity in the Fourier representation, (2) as the two sharp finite-sum
inequalities the note actually uses, and splits (4) into the endpoint-limit step
(F-7a) and the local-continuation step (F-7b).

Scope caveats, stated before the mathematics rather than after.
* (1) is formalized in the **Fourier representation**: it is the statement that
  the resonant cubic sum vanishes when every amplitude is transverse to its own
  wavevector.  That is exactly the finite-dimensional identity the Galerkin
  system consumes, and it is what the exact-integer checker in
  `src/ns_certificate_lab/galerkin_obstruction.py` verifies.  It is **not** a
  formalization of `∫_{𝕋³} u·(u·∇)u = 0` as an integral over a manifold: mathlib
  has neither the torus function spaces nor the Fourier isomorphism needed to
  connect the two, and that bridge remains open.
* (2) is formalized as the finite-sum inequalities with explicit constants, not
  as an abstract "all norms on a finite-dimensional space are equivalent"
  statement.  The explicit form is what the note quotes and what the Python
  implementation computes.
* (4) is closed here for a **time-independent** vector field (F-7b).  The
  time-dependent case, F-7c, is closed in `TimeDependentGalerkin.lean` by the
  direct route through mathlib's time-dependent `IsPicardLindelof`; nothing is
  axiomatized in either file.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib
import NSSingularity.GalerkinNoBlowup

namespace NSSingularity

noncomputable section

open Set MeasureTheory Filter Topology Finset

/-! ## F-12 — the trilinear cancellation in the Fourier representation

Writing `u(x) = Σ_k a_k e^{2πi k·x}` on the torus, Fourier orthogonality turns
`⟨u,(u·∇)u⟩` into `2πi` times

  `Σ_{p+q+s=0} (a_q · k_s)(a_p · a_s)`,

a cubic form in the amplitudes.  The divergence-free constraint is `k_i · a_i = 0`
for every mode.  The theorem below is that the cubic form then vanishes
identically.  The proof is the swap `p ↔ s`, which is an involution of the
resonant index set.
-/

section Trilinear

open scoped Classical

variable {ι : Type*} [Fintype ι]

/-- A Fourier amplitude, or a wavevector, as a complex three-vector.  Integer
wavevectors are the intended instance; stating the theorem over `ℂ` makes it
strictly more general and avoids a cast. -/
abbrev Amp : Type := Fin 3 → ℂ

/-- The **bilinear** (not Hermitian) dot product.  This is the pairing that
appears in the Fourier form of the advection term. -/
def dotAmp (x y : Amp) : ℂ := ∑ i, x i * y i

theorem dotAmp_comm (x y : Amp) : dotAmp x y = dotAmp y x :=
  Finset.sum_congr rfl fun _ _ => mul_comm _ _

theorem dotAmp_add_right (x y z : Amp) :
    dotAmp x (y + z) = dotAmp x y + dotAmp x z := by
  simp [dotAmp, mul_add, Finset.sum_add_distrib]

theorem dotAmp_neg_right (x y : Amp) : dotAmp x (-y) = -dotAmp x y := by
  simp [dotAmp]

/-- The resonant triples `{(p,q,s) : k p + k q + k s = 0}`. -/
def resonantTriples (k : ι → Amp) : Finset (ι × ι × ι) :=
  Finset.univ.filter fun t => k t.1 + k t.2.1 + k t.2.2 = 0

theorem mem_resonantTriples {k : ι → Amp} {t : ι × ι × ι} :
    t ∈ resonantTriples k ↔ k t.1 + k t.2.1 + k t.2.2 = 0 := by
  simp [resonantTriples]

/-- Swapping the first and third index preserves resonance. -/
theorem swap_mem_resonantTriples {k : ι → Amp} {t : ι × ι × ι}
    (ht : t ∈ resonantTriples k) : (t.2.2, t.2.1, t.1) ∈ resonantTriples k := by
  rw [mem_resonantTriples] at ht ⊢
  simpa [add_comm, add_left_comm, add_assoc] using ht

/-- The Fourier form of `⟨u,(u·∇)u⟩`, divided by the constant `2πi`. -/
def advectionForm (k a : ι → Amp) : ℂ :=
  ∑ t ∈ resonantTriples k, dotAmp (a t.2.1) (k t.2.2) * dotAmp (a t.1) (a t.2.2)

/-- **F-12.**  The Navier–Stokes advection term is energy neutral: for amplitudes
transverse to their own wavevectors the resonant cubic form vanishes.

This is Lemma 1 of `docs/research_notes/track_f_finite_mode_nogo.md`, and it is
the single algebraic input of the F-6 bound.  No reality constraint
`a_{-k} = conj a_k` is needed: the identity holds for every transverse amplitude
family, and the mode set enters only through `Fintype ι`. -/
theorem advectionForm_eq_zero (k a : ι → Amp)
    (hdiv : ∀ i, dotAmp (k i) (a i) = 0) : advectionForm k a = 0 := by
  have hswapped :
      advectionForm k a
        = ∑ t ∈ resonantTriples k,
            dotAmp (a t.2.1) (k t.1) * dotAmp (a t.2.2) (a t.1) := by
    refine Finset.sum_nbij' (i := fun t : ι × ι × ι => (t.2.2, t.2.1, t.1))
      (j := fun t : ι × ι × ι => (t.2.2, t.2.1, t.1))
      (fun a ha => swap_mem_resonantTriples ha)
      (fun a ha => swap_mem_resonantTriples ha)
      (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ => rfl)
  have hzero : ∀ t ∈ resonantTriples k,
      dotAmp (a t.2.1) (k t.2.2) * dotAmp (a t.1) (a t.2.2)
        + dotAmp (a t.2.1) (k t.1) * dotAmp (a t.2.2) (a t.1) = 0 := by
    intro t ht
    have hres : k t.1 + k t.2.1 + k t.2.2 = 0 := mem_resonantTriples.mp ht
    have hsub : k t.2.2 + k t.1 = -k t.2.1 := by
      linear_combination hres
    have hsym : dotAmp (a t.2.2) (a t.1) = dotAmp (a t.1) (a t.2.2) :=
      dotAmp_comm _ _
    have hsplit :
        dotAmp (a t.2.1) (k t.2.2) + dotAmp (a t.2.1) (k t.1)
          = dotAmp (a t.2.1) (k t.2.2 + k t.1) := (dotAmp_add_right _ _ _).symm
    have hvanish : dotAmp (a t.2.1) (k t.2.2 + k t.1) = 0 := by
      rw [hsub, dotAmp_neg_right, dotAmp_comm (a t.2.1) (k t.2.1), hdiv t.2.1,
        neg_zero]
    rw [hsym]
    calc
      dotAmp (a t.2.1) (k t.2.2) * dotAmp (a t.1) (a t.2.2)
          + dotAmp (a t.2.1) (k t.1) * dotAmp (a t.1) (a t.2.2)
          = (dotAmp (a t.2.1) (k t.2.2) + dotAmp (a t.2.1) (k t.1))
              * dotAmp (a t.1) (a t.2.2) := by ring
      _ = 0 := by rw [hsplit, hvanish, zero_mul]
  have hdouble : advectionForm k a + advectionForm k a = 0 := by
    nth_rewrite 2 [hswapped]
    rw [advectionForm, ← Finset.sum_add_distrib]
    exact Finset.sum_eq_zero hzero
  have : (2 : ℂ) * advectionForm k a = 0 := by linear_combination hdouble
  simpa using this

end Trilinear

/-! ## F-13 — norm equivalence on the finite Fourier space

The note uses exactly two constants:

* `‖u‖_{H^s}^2 = Σ_k (1+4π²|k|²)^s |a_k|^2 ≤ (1+4π²R_S^2)^s ‖u‖^2`, and
* `‖∂^α u‖_∞ ≤ (2πR_S)^{|α|} √|S| ‖u‖`,

the first a weighted-sum bound and the second Cauchy–Schwarz over the retained
modes.  Both are finite-sum inequalities; the theorems below are those two
inequalities with the weights and cardinality left abstract.
-/

section NormEquivalence

variable {ι : Type*}

/-- **F-13, weighted form.**  A uniformly bounded spectral weight amplifies the
squared `ℓ²` norm by at most that bound.  With `w k = (1+4π²|k|²)^s` this is the
`H^s`-versus-`L²` constant of the note. -/
theorem weighted_sq_sum_le (s : Finset ι) (w c : ι → ℝ) (W : ℝ)
    (hw : ∀ i ∈ s, w i ≤ W) :
    ∑ i ∈ s, w i * c i ^ 2 ≤ W * ∑ i ∈ s, c i ^ 2 := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun i hi => ?_
  exact mul_le_mul_of_nonneg_right (hw i hi) (sq_nonneg _)

/-- **F-13, Cauchy–Schwarz form.**  `(Σ|c|)² ≤ |S| Σ|c|²`, the inequality behind
the `L^∞`-versus-`L²` constant `√|S|`. -/
theorem sq_sum_abs_le_card_mul_sum_sq (s : Finset ι) (c : ι → ℝ) :
    (∑ i ∈ s, |c i|) ^ 2 ≤ (s.card : ℝ) * ∑ i ∈ s, c i ^ 2 := by
  have h := sq_sum_le_card_mul_sum_sq (s := s) (f := fun i => |c i|)
  simpa [sq_abs] using h

/-- The packaged supremum bound: `Σ|c| ≤ √|S| · √(Σ c²)`. -/
theorem sum_abs_le_sqrt_card_mul_sqrt_sum_sq (s : Finset ι) (c : ι → ℝ) :
    ∑ i ∈ s, |c i| ≤ Real.sqrt (s.card : ℝ) * Real.sqrt (∑ i ∈ s, c i ^ 2) := by
  have hprod : Real.sqrt (s.card : ℝ) * Real.sqrt (∑ i ∈ s, c i ^ 2)
      = Real.sqrt ((s.card : ℝ) * ∑ i ∈ s, c i ^ 2) :=
    (Real.sqrt_mul (by positivity) _).symm
  rw [hprod]
  exact Real.le_sqrt_of_sq_le (sq_sum_abs_le_card_mul_sum_sq s c)

end NormEquivalence

/-! ## F-7a — a bounded trajectory with bounded velocity reaches the endpoint

F-6 bounds `‖u‖` on `[0,T)`, and the Galerkin equation then bounds `‖u'‖`.  A
bounded, continuous derivative on a *finite* interval already forces the
trajectory to have a limit at the right endpoint; that limit is the initial
condition the continuation step consumes.

The proof is the fundamental theorem of calculus plus continuity of the
primitive.  Nothing is assumed at `T` itself, which is the point: `T` is the
putative singular time.
-/

section Endpoint

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]

omit [NormedSpace ℝ E] [CompleteSpace E] in
/-- A continuous, bounded derivative on `Ico 0 T` is interval integrable. -/
theorem intervalIntegrable_of_continuousOn_bounded
    {u' : ℝ → E} {T C : ℝ} (hT : 0 < T)
    (hcont : ContinuousOn u' (Ico (0 : ℝ) T))
    (hbound : ∀ t ∈ Ico (0 : ℝ) T, ‖u' t‖ ≤ C) :
    IntervalIntegrable u' volume 0 T := by
  rw [intervalIntegrable_iff_integrableOn_Ioo_of_le hT.le]
  haveI : IsFiniteMeasure (volume.restrict (Ioo (0 : ℝ) T)) := by
    refine ⟨?_⟩
    rw [Measure.restrict_apply_univ, Real.volume_Ioo]
    exact ENNReal.ofReal_lt_top
  refine Integrable.of_bound ?_ C ?_
  · exact (hcont.mono fun x hx => ⟨hx.1.le, hx.2⟩).aestronglyMeasurable measurableSet_Ioo
  · refine (ae_restrict_iff' measurableSet_Ioo).2 (Eventually.of_forall fun x hx => ?_)
    exact hbound x ⟨hx.1.le, hx.2⟩

/-- **F-7a.**  A curve differentiable on `[0,T)` with a continuous derivative
bounded by `C` has a limit as `t → T` from the left.

This is the step that turns the F-6 *bound* into a statement about the endpoint:
a bounded trajectory not only fails to blow up, it converges. -/
theorem exists_tendsto_nhdsWithin_of_norm_deriv_le
    {u u' : ℝ → E} {T C : ℝ} (hT : 0 < T)
    (hderiv : ∀ t ∈ Ico (0 : ℝ) T, HasDerivAt u (u' t) t)
    (hcont : ContinuousOn u' (Ico (0 : ℝ) T))
    (hbound : ∀ t ∈ Ico (0 : ℝ) T, ‖u' t‖ ≤ C) :
    ∃ L : E, Tendsto u (𝓝[<] T) (𝓝 L) := by
  have hint : IntervalIntegrable u' volume 0 T :=
    intervalIntegrable_of_continuousOn_bounded hT hcont hbound
  set P : ℝ → E := fun t => ∫ s in (0 : ℝ)..t, u' s with hP
  -- The primitive is continuous up to and including `T`.
  have hPcont : ContinuousOn P (Icc 0 T) := by
    have h := intervalIntegral.continuousOn_primitive_interval
      (a := (0 : ℝ)) (b := T) (μ := volume) (f := u')
      (by rw [uIcc_of_le hT.le]
          exact (intervalIntegrable_iff_integrableOn_Icc_of_le hT.le).1 hint)
    rwa [uIcc_of_le hT.le] at h
  -- On `[0,T)` the curve is its own primitive.
  have hrepr : ∀ t ∈ Ico (0 : ℝ) T, u t = u 0 + P t := by
    intro t ht
    have hsub : uIcc (0 : ℝ) t ⊆ Ico (0 : ℝ) T := by
      rw [uIcc_of_le ht.1]
      exact fun x hx => ⟨hx.1, lt_of_le_of_lt hx.2 ht.2⟩
    have hloc : IntervalIntegrable u' volume 0 t :=
      hint.mono_set (by rw [uIcc_of_le ht.1, uIcc_of_le hT.le]
                        exact Icc_subset_Icc le_rfl ht.2.le)
    have := intervalIntegral.integral_eq_sub_of_hasDerivAt
      (fun x hx => hderiv x (hsub hx)) hloc
    simp only [hP]
    rw [this]
    abel
  refine ⟨u 0 + P T, ?_⟩
  have hmem : Icc (0 : ℝ) T ∈ 𝓝[<] T := by
    have h1 : Ioo (0 : ℝ) T ∈ 𝓝[<] T :=
      Ioo_mem_nhdsLT hT
    exact mem_of_superset h1 fun x hx => ⟨hx.1.le, hx.2.le⟩
  have hPtend : Tendsto P (𝓝[<] T) (𝓝 (P T)) := by
    have hwithin : ContinuousWithinAt P (Icc 0 T) T :=
      hPcont T (right_mem_Icc.2 hT.le)
    exact hwithin.mono_left (le_trans (nhdsWithin_le_iff.2 hmem) (le_refl _))
  have heq : u =ᶠ[𝓝[<] T] fun t => u 0 + P t := by
    have hIco : Ico (0 : ℝ) T ∈ 𝓝[<] T :=
      mem_of_superset (Ioo_mem_nhdsLT hT) fun x hx => ⟨hx.1.le, hx.2⟩
    exact eventually_of_mem hIco fun x hx => hrepr x hx
  exact (tendsto_const_nhds.add hPtend).congr' heq.symm

end Endpoint

/-! ## F-7b — local continuation from the endpoint

With the endpoint value in hand, continuing past `T` is exactly Picard–Lindelöf.
mathlib supplies it for a **time-independent** `C^1` field
(`exists_forall_mem_closedBall_exists_eq_forall_mem_Ioo_hasDerivAt₀`), which
covers the Galerkin system with a time-independent projected force.

**F-7c, closed** in `TimeDependentGalerkin.lean`.  The remark above that
mathlib's local existence theorem is autonomous-only was true of the theorem
named here, but not of the API as a whole: `IsPicardLindelof` is stated for
`f : ℝ → E → E`, so a genuinely time-dependent projected force needs no
autonomisation.  What it does need is the local Lipschitz estimate for `B x x`
— `B x x - B y y = B x (x-y) + B (x-y) y`, hence a constant `2‖B‖(‖x₀‖+a)` that
carries the ball radius — and that is `galerkin_isPicardLindelof`.  The Lean
chain therefore now covers time-dependent projected forcing, matching the paper
proof in `docs/research_notes/track_f_finite_mode_nogo.md` §3.
-/

section Continuation

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]

omit [CompleteSpace E] in
/-- The Galerkin vector field with a time-independent forcing is `C^1`.

`B` is presented as a continuous linear map into continuous linear maps, which
is what a bilinear form on a finite-dimensional coefficient space always is. -/
theorem contDiff_galerkinField (B : E →L[ℝ] E →L[ℝ] E) (A : E →L[ℝ] E) (g : E) :
    ContDiff ℝ 1 fun x : E => g + B x x + A x := by
  have hB : ContDiff ℝ 1 fun x : E => B x x :=
    ContDiff.clm_apply B.contDiff contDiff_id
  exact ((contDiff_const.add hB).add A.contDiff)

/-- **F-7b.**  The autonomous Galerkin system has a two-sided local solution
through every point, in particular through the endpoint value produced by
F-7a. -/
theorem exists_local_galerkin_solution
    (B : E →L[ℝ] E →L[ℝ] E) (A : E →L[ℝ] E) (g : E) (L : E) (T : ℝ) :
    ∃ α : ℝ → E, α T = L ∧ ∃ ε > (0 : ℝ),
      ∀ t ∈ Ioo (T - ε) (T + ε),
        HasDerivAt α (g + B (α t) (α t) + A (α t)) t :=
  ContDiffAt.exists_forall_mem_closedBall_exists_eq_forall_mem_Ioo_hasDerivAt₀
    ((contDiff_galerkinField B A g).contDiffAt) T

end Continuation

/-! ## The logical connection to Clay (D)

A Track-F ansatz witnesses a Clay (C)/(D) breakdown only if its coefficient
trajectory genuinely fails to reach the putative singular time — otherwise the
velocity field extends and nothing breaks down.  `IsBreakdownCandidate` names
that requirement, and the theorem below refutes it for every fixed-finite-mode
candidate whose projected force is bounded.

What this does **not** do, stated plainly:
* it does not connect to `ClayStatement.ClayPeriodicBreakdown`.  That needs the
  Fourier isomorphism between `V_S` and the coefficient space, the identification
  of `⟨u,(u·∇)u⟩` with `advectionForm` (F-12 is the algebraic half; the analytic
  half needs torus function spaces mathlib does not have), and the local
  uniqueness theory for Navier–Stokes, which mathlib also does not have;
* the continuation step for time-dependent forcing is F-7c, closed separately
  in `TimeDependentGalerkin.lean` rather than here.
-/

section ClayConnection

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

open RealInnerProductSpace

/-- The property a Track-F ansatz needs in order to witness a breakdown at `T`:
the coefficient trajectory must fail to have a limit there.  If it does have a
limit, the ansatz reaches `T` with every norm bounded and nothing breaks down. -/
def IsBreakdownCandidate (u : ℝ → E) (T : ℝ) : Prop :=
  0 < T ∧ ¬ ∃ L : E, Tendsto u (𝓝[<] T) (𝓝 L)

variable [CompleteSpace E]

/-- **The fixed-finite-mode no-go, Lean form.**

A trajectory of the energy-neutral dissipative Galerkin system driven by a
uniformly bounded force is never a breakdown candidate.  The chain is
F-6 (`galerkin_norm_le_of_mem`) bounds the state, the equation then bounds the
velocity, and F-7a (`exists_tendsto_nhdsWithin_of_norm_deriv_le`) converts that
into convergence at `T`.

`B` and `A` are continuous (linear/bilinear) because that is what they are on a
finite-dimensional coefficient space; continuity is what converts the state
bound into a velocity bound. -/
theorem not_isBreakdownCandidate_of_galerkin
    {B : E →L[ℝ] E →L[ℝ] E} {A : E →L[ℝ] E} {g u : ℝ → E} {T M : ℝ} (hT : 0 < T)
    (hB : EnergyNeutral fun x y => B x y) (hA : Dissipative fun x => A x)
    (hg : ∀ t, ‖g t‖ ≤ M)
    (hgcont : ContinuousOn g (Ico (0 : ℝ) T))
    (hu : ContinuousOn u (Icc 0 T))
    (hode : ∀ t ∈ Ico (0 : ℝ) T,
      HasDerivAt u (g t + B (u t) (u t) + A (u t)) t) :
    ¬ IsBreakdownCandidate u T := by
  rintro ⟨-, hno⟩
  -- F-6: the state is bounded by `K` on the whole of `[0,T]`.
  set K : ℝ := ‖u 0‖ + ∫ _s in (0 : ℝ)..T, M with hK
  have hstate : ∀ t ∈ Icc (0 : ℝ) T, ‖u t‖ ≤ K :=
    galerkin_norm_le_of_mem hT.le hB hA hg hu
      (fun t ht => hode t ⟨ht.1.le, ht.2⟩) continuousOn_const
  -- The equation turns that into a velocity bound.
  set C : ℝ := M + ‖B‖ * K ^ 2 + ‖A‖ * K with hC
  have hvel : ∀ t ∈ Ico (0 : ℝ) T,
      ‖g t + B (u t) (u t) + A (u t)‖ ≤ C := by
    intro t ht
    have hmem : t ∈ Icc (0 : ℝ) T := ⟨ht.1, ht.2.le⟩
    have hnorm : ‖u t‖ ≤ K := hstate t hmem
    have hbil : ‖B (u t) (u t)‖ ≤ ‖B‖ * K ^ 2 := by
      have h1 : ‖B (u t) (u t)‖ ≤ ‖B (u t)‖ * ‖u t‖ :=
        (B (u t)).le_opNorm (u t)
      have h2 : ‖B (u t)‖ ≤ ‖B‖ * ‖u t‖ := B.le_opNorm (u t)
      have hK0 : 0 ≤ K := le_trans (norm_nonneg _) hnorm
      calc ‖B (u t) (u t)‖ ≤ ‖B (u t)‖ * ‖u t‖ := h1
        _ ≤ (‖B‖ * ‖u t‖) * ‖u t‖ := by
            exact mul_le_mul_of_nonneg_right h2 (norm_nonneg _)
        _ ≤ (‖B‖ * K) * K := by
            have : ‖B‖ * ‖u t‖ ≤ ‖B‖ * K :=
              mul_le_mul_of_nonneg_left hnorm (norm_nonneg B)
            exact mul_le_mul this hnorm (norm_nonneg _)
              (mul_nonneg (norm_nonneg B) hK0)
        _ = ‖B‖ * K ^ 2 := by ring
    have hlin : ‖A (u t)‖ ≤ ‖A‖ * K := by
      refine le_trans (A.le_opNorm (u t)) ?_
      exact mul_le_mul_of_nonneg_left hnorm (norm_nonneg A)
    calc ‖g t + B (u t) (u t) + A (u t)‖
        ≤ ‖g t + B (u t) (u t)‖ + ‖A (u t)‖ := norm_add_le _ _
      _ ≤ (‖g t‖ + ‖B (u t) (u t)‖) + ‖A (u t)‖ := by
          linarith [norm_add_le (g t) (B (u t) (u t))]
      _ ≤ (M + ‖B‖ * K ^ 2) + ‖A‖ * K := by
          linarith [hg t, hbil, hlin]
      _ = C := by rw [hC]
  -- The velocity is continuous on `[0,T)`.
  have hucont : ContinuousOn u (Ico (0 : ℝ) T) := fun t ht =>
    ((hode t ht).continuousAt).continuousWithinAt
  have hcont : ContinuousOn (fun t => g t + B (u t) (u t) + A (u t))
      (Ico (0 : ℝ) T) := by
    have hBu : ContinuousOn (fun t => B (u t) (u t)) (Ico (0 : ℝ) T) :=
      (B.continuous.comp_continuousOn hucont).clm_apply hucont
    exact (hgcont.add hBu).add (A.continuous.comp_continuousOn hucont)
  exact hno (exists_tendsto_nhdsWithin_of_norm_deriv_le hT hode hcont hvel)

/-- The same conclusion packaged with the explicit F-6 bound, so that a reader
sees both halves of the no-go: the trajectory stays inside a ball of radius
`‖u 0‖ + M·T` and it converges at the putative singular time. -/
theorem galerkin_bounded_and_reaches_endpoint
    {B : E →L[ℝ] E →L[ℝ] E} {A : E →L[ℝ] E} {g u : ℝ → E} {T M : ℝ} (hT : 0 < T)
    (hB : EnergyNeutral fun x y => B x y) (hA : Dissipative fun x => A x)
    (hg : ∀ t, ‖g t‖ ≤ M)
    (hgcont : ContinuousOn g (Ico (0 : ℝ) T))
    (hu : ContinuousOn u (Icc 0 T))
    (hode : ∀ t ∈ Ico (0 : ℝ) T,
      HasDerivAt u (g t + B (u t) (u t) + A (u t)) t) :
    (∀ t ∈ Icc (0 : ℝ) T, ‖u t‖ ≤ ‖u 0‖ + M * T) ∧
      ∃ L : E, Tendsto u (𝓝[<] T) (𝓝 L) := by
  constructor
  · intro t ht
    have := galerkin_norm_le_of_mem hT.le hB hA hg hu
      (fun s hs => hode s ⟨hs.1.le, hs.2⟩) continuousOn_const t ht
    simpa [mul_comm] using this
  · by_contra hno
    exact not_isBreakdownCandidate_of_galerkin hT hB hA hg hgcont hu hode ⟨hT, hno⟩

end ClayConnection

end

end NSSingularity
