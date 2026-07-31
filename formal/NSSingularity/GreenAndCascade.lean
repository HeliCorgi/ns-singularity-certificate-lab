/-
Stage 1, items **F-14**, **F-15**, **F-16** and the **F-7c reduction** of
`docs/final_target.md` §4.

Four small, independent pieces, each chosen because it is the *algebraic core*
of a step that the numerical side already relies on:

* **F-14** — the radial identity behind the five-dimensional fundamental
  solution: `R^{-3}` is annihilated by `d_RR + 4 d_R / R`, which is `Δ_5` acting
  on functions of the five-dimensional radius alone.  Every tail bound in
  `src/ns_certificate_lab/free_space_recovery.py` starts here.
* **F-15** — Newton's flux identity `R^4 ψ'(R) = -m(R)`, the exact statement
  that turns a compactly supported radial source into a closed-form potential.
  This is what makes the Gate-4 reference independent rather than quadratured.
* **F-16** — the shell exponent region with **every hypothesis named**: the
  bandwidth exponent of a diverging-bandwidth Track-F candidate satisfies
  `γ < 1`.  Stating it in Lean is what makes the "under the current shell ansatz
  and non-degeneracy assumptions" qualifier machine-checkable rather than a
  footnote.
* **F-7c, reduction** — a step of the *abandoned* autonomisation route, kept as
  a record of it: if the autonomised field on `E × ℝ` has a local flow, the
  time-dependent
  Galerkin system has one.  The remaining half is constructing that flow.

Scope caveats.
* F-14 and F-15 are one-dimensional real-analysis identities about the radial
  profile.  They are **not** a formalization of `Δ_5 G_5 = -δ_0`, which would
  need distributions on `ℝ^5`; mathlib has no such development and the Dirac
  side is not used by any bound in this repository.
* F-16 is arithmetic on exponents.  It does not derive the four inequalities
  from the PDE — that derivation is on paper in
  `docs/research_notes/track_f_shell_constraints.md` and quotes the endpoint
  regularity theorem.
* The F-7c reduction does **not** close F-7c; `TimeDependentGalerkin.lean` does, by a shorter route.  This splits it, which is what
  `docs/final_target.md` §4.1 asked for.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/
import Mathlib
import NSSingularity.GalerkinNoBlowup
import NSSingularity.FiniteModeNoGo

namespace NSSingularity

noncomputable section

open Set Filter Topology

/-! ## F-14 — the radial five-dimensional Green profile is harmonic off the origin

`Δ_5` acting on a function of `R = |X|` alone is `f'' + 4 f' / R`.  The claim is
that `f(R) = R^{-3}` is annihilated by it, which is the only property of the
fundamental solution that the derivative tail bounds use.
-/

section GreenProfile

/-- The radial profile of the five-dimensional fundamental solution, up to the
constant `c_5 = 1/(8π²)` which plays no role in the identity. -/
def greenProfile (R : ℝ) : ℝ := R ^ (-3 : ℤ)

/-- `d/dR R^{-3} = -3 R^{-4}`, in the raw form mathlib produces. -/
theorem hasDerivAt_greenProfile {R : ℝ} (hR : R ≠ 0) :
    HasDerivAt greenProfile ((-3 : ℤ) * R ^ (-3 - 1 : ℤ)) R :=
  hasDerivAt_zpow (-3 : ℤ) R (Or.inl hR)

/-- `d/dR (-3 R^{-4}) = 12 R^{-5}`. -/
theorem hasDerivAt_greenProfileDeriv {R : ℝ} (hR : R ≠ 0) :
    HasDerivAt (fun x : ℝ => (-3 : ℝ) * x ^ (-4 : ℤ)) (12 * R ^ (-5 : ℤ)) R := by
  have h := (hasDerivAt_zpow (-4 : ℤ) R (Or.inl hR)).const_mul (-3 : ℝ)
  have hval : (-3 : ℝ) * ((-4 : ℤ) * R ^ (-4 - 1 : ℤ)) = 12 * R ^ (-5 : ℤ) := by
    norm_num
    ring
  rw [hval] at h
  exact h

/-- **F-14.**  `R^{-3}` solves the five-dimensional radial Laplace equation
`f'' + 4 f'/R = 0` away from the origin.  The two derivative values are the ones
computed above, so the identity is about the actual derivatives of
`greenProfile`, not about symbols. -/
theorem greenProfile_radial_laplace_eq_zero {R : ℝ} (hR : R ≠ 0) :
    12 * R ^ (-5 : ℤ) + 4 / R * ((-3 : ℝ) * R ^ (-4 : ℤ)) = 0 := by
  have h4 : R ^ (-4 : ℤ) = (R ^ (4 : ℕ))⁻¹ := by
    rw [show (-4 : ℤ) = -((4 : ℕ) : ℤ) by norm_num, zpow_neg, zpow_natCast]
  have h5 : R ^ (-5 : ℤ) = (R ^ (5 : ℕ))⁻¹ := by
    rw [show (-5 : ℤ) = -((5 : ℕ) : ℤ) by norm_num, zpow_neg, zpow_natCast]
  rw [h4, h5]
  field_simp
  ring

end GreenProfile

/-! ## F-15 — Newton's flux identity

For a radial source `ω` on `ℝ^5` with enclosed mass `m(R) = ∫_0^R t^4 ω(t) dt`,
the potential slope is `ψ'(R) = -m(R)/R^4`, equivalently the **flux** `R^4 ψ'(R)`
equals `-m(R)`.  Differentiating gives `(R^4 ψ')' = -R^4 ω`, which is the radial
form of `-Δ_5 ψ = ω`.  The closed-form Gate-4 reference is exactly this identity
applied to a polynomial bump.
-/

section NewtonFlux

/-- The potential slope produced by an enclosed mass `m`. -/
def newtonSlope (m : ℝ → ℝ) (R : ℝ) : ℝ := -m R / R ^ 4

/-- **F-15, flux form.**  `R^4 ψ'(R) = -m(R)` for `R ≠ 0`. -/
theorem flux_newtonSlope (m : ℝ → ℝ) {R : ℝ} (hR : R ≠ 0) :
    R ^ 4 * newtonSlope m R = -m R := by
  simp only [newtonSlope]
  field_simp

/-- **F-15, differentiated form.**  If `m` has derivative `R^4 ω` at `R`, then
the flux `R ↦ R^4 ψ'(R)` has derivative `-R^4 ω`, i.e. the radial form of
`-Δ_5 ψ = ω`.  Stated through the flux so that no quotient is differentiated. -/
theorem hasDerivAt_flux (m : ℝ → ℝ) {R w : ℝ} (hm : HasDerivAt m (R ^ 4 * w) R) :
    HasDerivAt (fun x : ℝ => -m x) (-(R ^ 4 * w)) R :=
  hm.neg

end NewtonFlux

/-! ## F-16 — the shell exponent region, with every hypothesis named

`docs/research_notes/track_f_shell_constraints.md` models the dyadic shell
energies of a diverging-bandwidth Track-F candidate by
`E_j = A (λ_j/N)^β`, `N = (T-t)^{-γ}`, `A = (T-t)^σ`, and derives four
inequalities.  The theorem below is the arithmetic consequence.  Naming the
hypotheses in Lean is the point: it records that `γ < 1` holds **under those
assumptions**, not unconditionally.
-/

section ShellExponents

/-- The four necessary conditions of the shell ansatz, as one structure so that
none of them can be dropped silently. -/
structure ShellAdmissible (γ σ β : ℝ) : Prop where
  /-- The bandwidth diverges: the contrapositive of the fixed-finite-mode no-go. -/
  bandwidth_diverges : 0 < γ
  /-- The shell sum converges, so the energy is finite. -/
  spectrum_summable : 0 < β
  /-- The total energy stays bounded as `t → T`. -/
  energy_bounded : 0 ≤ σ
  /-- The total dissipation is integrable. -/
  dissipation_integrable : -1 < σ - 2 * γ
  /-- The critical Besov norm actually diverges (uses the endpoint regularity
  theorem, which is **quoted** and not formalized). -/
  critical_norm_diverges : σ < γ

/-- **F-16.**  Under the shell ansatz the bandwidth exponent is below one. -/
theorem ShellAdmissible.bandwidth_lt_one {γ σ β : ℝ} (h : ShellAdmissible γ σ β) :
    γ < 1 := by
  have h1 := h.dissipation_integrable
  have h2 := h.critical_norm_diverges
  linarith

/-- The admissible `σ` interval, made explicit. -/
theorem ShellAdmissible.sigma_mem {γ σ β : ℝ} (h : ShellAdmissible γ σ β) :
    σ ∈ Ico (max 0 (2 * γ - 1)) γ := by
  refine ⟨?_, h.critical_norm_diverges⟩
  exact max_le h.energy_bounded (le_of_lt (by linarith [h.dissipation_integrable]))

/-- No admissible point exists once `γ ≥ 1`. -/
theorem not_shellAdmissible_of_one_le {γ σ β : ℝ} (hγ : 1 ≤ γ) :
    ¬ ShellAdmissible γ σ β := fun h => absurd h.bandwidth_lt_one (not_lt.mpr hγ)

end ShellExponents

/-! ## Limited connection to Clay (D)

`not_isBreakdownCandidate_of_galerkin` refutes a breakdown at one time `T`.
The sharper statement a Clay-(D) discussion needs is that the *set* of breakdown
times is empty: a fixed-finite-bandwidth candidate never witnesses a breakdown,
at any time.  That is what the theorem below says.

What still separates this from `ClayStatement.ClayPeriodicBreakdown` is
unchanged and is listed in `formal/README.md`: the Fourier isomorphism between
the trigonometric-polynomial space and the coefficient space, the identification
of `⟨u,(u·∇)u⟩` with `advectionForm` (F-12 is its algebraic half), and the local
uniqueness theory for Navier–Stokes, none of which exist in mathlib.
-/

section ClayConnection

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [CompleteSpace E]

open RealInnerProductSpace

/-- **Limited Clay connection.**  A fixed-finite-bandwidth Galerkin trajectory
driven by a uniformly bounded force has an **empty** set of breakdown times. -/
theorem breakdown_time_set_empty
    {B : E →L[ℝ] E →L[ℝ] E} {A : E →L[ℝ] E} {g u : ℝ → E} {M : ℝ}
    (hB : EnergyNeutral fun x y => B x y) (hA : Dissipative fun x => A x)
    (hg : ∀ t, ‖g t‖ ≤ M)
    (hgcont : ContinuousOn g (Ici (0 : ℝ)))
    (hu : ContinuousOn u (Ici (0 : ℝ)))
    (hode : ∀ t ∈ Ici (0 : ℝ),
      HasDerivAt u (g t + B (u t) (u t) + A (u t)) t) :
    {T : ℝ | IsBreakdownCandidate u T} = ∅ := by
  ext T
  simp only [mem_setOf_eq, mem_empty_iff_false, iff_false]
  intro hbreak
  have hT : 0 < T := hbreak.1
  refine not_isBreakdownCandidate_of_galerkin (B := B) (A := A) hT hB hA hg
    (hgcont.mono fun t ht => ht.1)
    (hu.mono fun t ht => ht.1)
    (fun t ht => hode t ht.1) hbreak

end ClayConnection

/-! ## F-7c, the reduction step

The theorem below is a step of the autonomisation route: given a flow for the
autonomised field on `E × ℝ`, the original system has a local solution.

That route was abandoned.  `IsPicardLindelof` in the pinned mathlib is already
time-dependent, so the direct construction in `TimeDependentGalerkin.lean`
closes F-7c without ever forming `F(x,s) = (g s + B x x + A x, 1)`, without the
product-space instances, and without this reduction.  The theorem is kept
because it is correct and because it records what the alternative would have
cost; nothing depends on it.
-/

section Autonomisation

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- **F-7c, reduction.**  If the autonomised field has a local flow through
`(L, T)`, the time-dependent Galerkin system has a local solution through `L`.

The proof is that the second component solves `s' = 1` with `s T = T`, hence
`s t = t` on the interval, so the first component solves the original equation.
-/
theorem galerkin_solution_of_autonomised
    {B : E →L[ℝ] E →L[ℝ] E} {A : E →L[ℝ] E} {g : ℝ → E} {L : E} {T ε : ℝ}
    (hε : 0 < ε) {α : ℝ → E × ℝ}
    (hinit : α T = (L, T))
    (hflow : ∀ t ∈ Ioo (T - ε) (T + ε),
      HasDerivAt α (g (α t).2 + B (α t).1 (α t).1 + A (α t).1, (1 : ℝ)) t) :
    ∃ x : ℝ → E, x T = L ∧
      ∀ t ∈ Ioo (T - ε) (T + ε),
        HasDerivAt x (g t + B (x t) (x t) + A (x t)) t := by
  set I : Set ℝ := Ioo (T - ε) (T + ε) with hI
  have hTmem : T ∈ I := ⟨by linarith, by linarith⟩
  -- The second component minus the identity has zero derivative on `I`.
  have hsecond : ∀ t ∈ I, HasDerivAt (fun s => (α s).2) (1 : ℝ) t := by
    intro t ht
    exact (hflow t ht).snd
  have hzero : ∀ t ∈ I,
      HasDerivWithinAt (fun s => (α s).2 - s) ((fun _ : ℝ => (0 : ℝ)) t) I t := by
    intro t ht
    have h1 : HasDerivAt (fun s => (α s).2 - s) ((1 : ℝ) - 1) t :=
      (hsecond t ht).sub (hasDerivAt_id t)
    have h2 : HasDerivAt (fun s => (α s).2 - s) (0 : ℝ) t := by simpa using h1
    exact h2.hasDerivWithinAt
  have hconst : ∀ t ∈ I, (α t).2 - t = (α T).2 - T := by
    intro t ht
    have hbound : ∀ s ∈ I, ‖(fun _ : ℝ => (0 : ℝ)) s‖ ≤ 0 := by simp
    have hmain :=
      (convex_Ioo (T - ε) (T + ε)).norm_image_sub_le_of_norm_hasDerivWithin_le
        hzero hbound hTmem ht
    have h0 : ‖((α t).2 - t) - ((α T).2 - T)‖ ≤ 0 := by simpa using hmain
    have := norm_le_zero_iff.mp h0
    linarith [sub_eq_zero.mp this]
  have hid : ∀ t ∈ I, (α t).2 = t := by
    intro t ht
    have h := hconst t ht
    have hT2 : (α T).2 = T := by simp [hinit]
    rw [hT2] at h
    linarith
  refine ⟨fun t => (α t).1, by simp [hinit], ?_⟩
  intro t ht
  have hfst := (hflow t ht).fst
  rwa [hid t ht] at hfst

end Autonomisation

end

end NSSingularity
