/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import Mathlib

/-!
# Track P — the chain-composition skeleton

`src/ns_certificate_lab/torus_chain.py` chains Track-P slab certificates: slab
`n` certifies `‖u(t) - c_n(t)‖ ≤ R_n(t)` on `[t_n, t_{n+1}]` *given* the datum
bound `‖u(t_n) - c_n(t_n)‖ ≤ δ_n`, and the next datum radius is produced by the
exact recurrence `δ_{n+1} = R_n(t_{n+1}) + ‖c_n(t_{n+1}) - y_{n+1}‖`, where
`y_{n+1} = c_{n+1}(t_{n+1})` is the exact rational recentring of the endpoint.
This file is the finite-inequality skeleton of that composition.  Everything
here is an unconditional statement about an arbitrary function `u : ℝ → E`
into a normed group and arbitrary centre/radius data: the *supply* of the
per-slab conditional tubes for the true Navier–Stokes solution is the business
of the analytic estimate layer together with the named external classical
theorems EXT-P1/P2/P3, which are recorded in the research notes and **never as
Lean axioms**.  In the project's fixed wording: the payload downstream of this
file is a conditional PDE certificate assuming EXT-P1/P2/P3, while the
finite-dimensional Galerkin enclosure and the scalar control ODE are verified
unconditionally.  Nothing in this file, alone or in combination with the rest
of the library, is a statement about the Clay problem, and a chain that stops
proves nothing about singularities.

## What is formalised

* `two_slab_composition` — the two-slab step, exactly as the Python docstring
  states it: a certified tube on `[t₀, t₁]`, a transfer bound `η` at `t₁`, the
  budget inequality `R₁ t₁ + η ≤ δ₂`, and a *conditional* tube on `[t₁, t₂]`
  whose condition is the datum bound at `t₁`, compose to a piecewise tube on
  `[t₀, t₂]`, with piecewise centre and radius given literally by
  `if t ≤ t₁ then … else …`.
* `transfer_triangle` — the three-term triangle inequality
  `‖u - y‖ ≤ ‖u - c‖ + ‖c - m‖ + ‖m - y‖` behind the recentring budget: tube
  endpoint (`u` against the Galerkin centre `c`), Taylor remainder (`c` against
  the evaluated Taylor endpoint `m`), and dyadic rounding plus Leray projection
  (`m` against the new exact rational datum `y`).
* `ChainLink`, `LinkCertified`, `LinkComposable`, `chain_composition` — the
  `n`-slab induction.  A link carries its slab `[tStart, tEnd]`, its centre and
  radius functions and its datum radius `δ`; `LinkCertified` is the per-slab
  *conditional* tube (the certificate's output), and `LinkComposable` is the
  contiguity `L'.tStart = L.tEnd` together with the delta recurrence
  `L.radius L.tEnd + ‖L.centre L.tEnd - L'.centre L'.tStart‖ ≤ L'.δ`, which is
  the Lean transcription of the Python `delta_out = delta_end + transfer`
  (with `≤` so that the exact rational equality is a special case).  The
  conclusion is the per-slab tube on every link of the chain — precisely the
  Python docstring's "the certified interval is the union of the slabs and the
  certified radius is the per-slab tube".  This per-link form is the honest
  piecewise statement: at a shared boundary time *both* adjacent bounds hold,
  which is strictly stronger than any single `if`/`else` selection, and
  `two_slab_composition` exhibits the `if`/`else` shape in the base case.
* `chain_composition_union` — the same conclusion re-packaged over the union:
  every time covered by some slab admits a link whose tube bounds `u` there.
* `discrete_gronwall` — the discrete Grönwall inequality
  `x (n+1) ≤ A x n + B → x n ≤ Aⁿ x 0 + B Σ_{i<n} Aⁱ`, the a priori shape of
  the chained datum radii when every slab is run with the same contraction
  factor `A` and per-slab injection `B`.
* `piecewise_radius_le_max`, `le_foldr_max`, `chain_radius_le_foldr_max` — the
  chained bound never exceeds the maximum of the per-slab radius bounds:
  pointwise for the two-slab `if`/`else` radius, and through
  `List.foldr max 0` for a whole chain.

## What is not formalised, and why

* **The per-slab tube for the true solution.**  `LinkCertified` is a
  hypothesis here.  Discharging it for periodic Navier–Stokes needs the `H^4`
  energy estimate and EXT-P1/P2/P3 (Kato/Temam local existence and uniqueness,
  regularity, `H^4` continuation) — external classical theorems that this
  repository records faithfully in docs and payloads and never axiomatises.
  In particular the fact that the solution continued across a slab boundary is
  the *same* solution is the uniqueness half of EXT-P1 and enters nowhere
  below: the skeleton needs only one fixed function `u`.
* **The numerical values.**  Radii, transfers and deltas are produced by the
  Python interval layer; this file proves only that the finitely many
  inequalities it emits compose as claimed.

This file contains no `sorry`, no `admit`, and introduces no axioms.
-/

namespace NSSingularity

open Set
open scoped Nat

section TwoSlab

variable {E : Type*} [NormedAddCommGroup E]

/-- **C1, the two-slab step.**  A certified tube on `[t₀, t₁]`, a transfer
bound at the boundary, the budget inequality `R₁ t₁ + η ≤ δ₂`, and a
*conditional* tube on `[t₁, t₂]` — conditional exactly on the datum bound that
the first three items establish — compose to a piecewise tube on `[t₀, t₂]`.
The piecewise centre and radius are literally `if t ≤ t₁ then … else …`; at
`t = t₁` the first branch is chosen, and the datum bound at `t₁` (the
`hstart` of the proof) is exactly the quantity the Python chain calls
`δ_{n+1}` when `c₂ t₁` is the recentred datum.  Nothing about differential
equations enters: this is the finite-inequality shell into which the
conditional PDE certificate assuming EXT-P1/P2/P3 is plugged. -/
theorem two_slab_composition {u c₁ c₂ : ℝ → E} {R₁ R₂ : ℝ → ℝ}
    {t₀ t₁ t₂ : ℝ} {η δ₂ : ℝ}
    (h01 : t₀ ≤ t₁) (_h12 : t₁ ≤ t₂)
    (h1 : ∀ t ∈ Icc t₀ t₁, ‖u t - c₁ t‖ ≤ R₁ t)
    (htr : ‖c₁ t₁ - c₂ t₁‖ ≤ η)
    (hd : R₁ t₁ + η ≤ δ₂)
    (h2 : ‖u t₁ - c₂ t₁‖ ≤ δ₂ → ∀ t ∈ Icc t₁ t₂, ‖u t - c₂ t‖ ≤ R₂ t) :
    ∀ t ∈ Icc t₀ t₂,
      ‖u t - (if t ≤ t₁ then c₁ t else c₂ t)‖ ≤ if t ≤ t₁ then R₁ t else R₂ t := by
  have hsplit : u t₁ - c₂ t₁ = (u t₁ - c₁ t₁) + (c₁ t₁ - c₂ t₁) := by abel
  have htri : ‖u t₁ - c₂ t₁‖ ≤ ‖u t₁ - c₁ t₁‖ + ‖c₁ t₁ - c₂ t₁‖ := by
    rw [hsplit]; exact norm_add_le _ _
  have h1e := h1 t₁ ⟨h01, le_rfl⟩
  have hstart : ‖u t₁ - c₂ t₁‖ ≤ δ₂ := by linarith
  have htube₂ := h2 hstart
  intro t ht
  by_cases hcase : t ≤ t₁
  · simp only [if_pos hcase]
    exact h1 t ⟨ht.1, hcase⟩
  · simp only [if_neg hcase]
    exact htube₂ t ⟨(not_le.1 hcase).le, ht.2⟩

/-- **C2, the transfer triangle.**  The recentring budget at a slab boundary
splits into three finite quantities: the tube endpoint `‖u - c‖` (bounded by
`R_n(t_{n+1})`), the Taylor-endpoint remainder `‖c - m‖`, and the rounding
plus projection displacement `‖m - y‖` to the new exact rational datum.  The
Python layer bounds each term by interval arithmetic; this is the inequality
that lets it add them. -/
theorem transfer_triangle (u c m y : E) :
    ‖u - y‖ ≤ ‖u - c‖ + ‖c - m‖ + ‖m - y‖ := by
  have hsplit : u - y = (u - c) + (c - m) + (m - y) := by abel
  calc ‖u - y‖ = ‖(u - c) + (c - m) + (m - y)‖ := by rw [hsplit]
    _ ≤ ‖(u - c) + (c - m)‖ + ‖m - y‖ := norm_add_le _ _
    _ ≤ ‖u - c‖ + ‖c - m‖ + ‖m - y‖ := by
        have := norm_add_le (u - c) (c - m); linarith

end TwoSlab

/-! ## C3 — the `n`-slab chain -/

section Chain

/-- One slab link of a Track-P chain: the slab `[tStart, tEnd]`, the certified
centre trajectory, the certified radius, and the datum radius `δ` at `tStart`.
Pure data — the hypotheses live in `LinkCertified` and `LinkComposable` so
that they cannot be smuggled. -/
structure ChainLink (E : Type*) where
  /-- Left endpoint of the slab. -/
  tStart : ℝ
  /-- Right endpoint of the slab. -/
  tEnd : ℝ
  /-- The certified centre trajectory (the Galerkin surrogate of the slab). -/
  centre : ℝ → E
  /-- The certified radius (the control-ODE supersolution of the slab). -/
  radius : ℝ → ℝ
  /-- The datum radius: the certificate's tube is conditional on
  `‖u tStart - centre tStart‖ ≤ δ`. -/
  δ : ℝ

variable {E : Type*} [NormedAddCommGroup E]

/-- The per-slab certificate, as a *conditional* statement: the slab is
nondegenerate, and IF the datum bound holds at `tStart` THEN the tube holds on
the whole slab.  For the true periodic solution this conditional is exactly
what the slab certificate provides under EXT-P1/P2/P3; here it is a
hypothesis. -/
def LinkCertified (u : ℝ → E) (L : ChainLink E) : Prop :=
  L.tStart ≤ L.tEnd ∧
    (‖u L.tStart - L.centre L.tStart‖ ≤ L.δ →
      ∀ t ∈ Icc L.tStart L.tEnd, ‖u t - L.centre t‖ ≤ L.radius t)

/-- Contiguity and the delta recurrence between consecutive links: the next
slab starts where this one ends, and the next datum radius absorbs the tube
endpoint plus the recentring transfer — the Lean form of the Python
`delta_out = delta_end + transfer`, with `≤` so the exact equality is a
special case. -/
def LinkComposable (L L' : ChainLink E) : Prop :=
  L'.tStart = L.tEnd ∧
    L.radius L.tEnd + ‖L.centre L.tEnd - L'.centre L'.tStart‖ ≤ L'.δ

/-- **C3, the chain composition.**  If every link of the chain carries the
conditional per-slab certificate, consecutive links are composable, and the
datum bound holds at the head of the chain, then the tube holds on *every*
slab of the chain.  This is the `n`-fold iteration of `two_slab_composition`:
the induction transports the datum bound across each boundary through the
triangle inequality and the delta recurrence, so only the *first* datum bound
is ever assumed.  The conclusion is the piecewise bound on the whole union of
slabs in its strongest form: at a boundary time shared by two links both
tubes hold. -/
theorem chain_composition (u : ℝ → E) :
    ∀ links : List (ChainLink E),
      (∀ L ∈ links, LinkCertified u L) →
      links.IsChain LinkComposable →
      (∀ L₀ ∈ links.head?, ‖u L₀.tStart - L₀.centre L₀.tStart‖ ≤ L₀.δ) →
      ∀ L ∈ links, ∀ t ∈ Icc L.tStart L.tEnd, ‖u t - L.centre t‖ ≤ L.radius t := by
  intro links
  induction links with
  | nil => intro _ _ _ L hL; simp at hL
  | cons L₀ rest ih =>
    intro hcert hchain hd M hM
    have hd₀ : ‖u L₀.tStart - L₀.centre L₀.tStart‖ ≤ L₀.δ := hd L₀ rfl
    have hcL : LinkCertified u L₀ := hcert L₀ (by simp)
    have htube₀ := hcL.2 hd₀
    rcases List.mem_cons.1 hM with rfl | hM'
    · exact htube₀
    · cases rest with
      | nil => simp at hM'
      | cons L₁ rest' =>
        have hrel : LinkComposable L₀ L₁ := (List.isChain_cons_cons.1 hchain).1
        have hchain' : (L₁ :: rest').IsChain LinkComposable :=
          (List.isChain_cons_cons.1 hchain).2
        have hd₁ : ‖u L₁.tStart - L₁.centre L₁.tStart‖ ≤ L₁.δ := by
          have hend := htube₀ L₀.tEnd ⟨hcL.1, le_rfl⟩
          have hrec := hrel.2
          have hsplit : u L₀.tEnd - L₁.centre L₀.tEnd
              = (u L₀.tEnd - L₀.centre L₀.tEnd) + (L₀.centre L₀.tEnd - L₁.centre L₀.tEnd) := by
            abel
          have htri : ‖u L₀.tEnd - L₁.centre L₀.tEnd‖
              ≤ ‖u L₀.tEnd - L₀.centre L₀.tEnd‖ + ‖L₀.centre L₀.tEnd - L₁.centre L₀.tEnd‖ := by
            rw [hsplit]; exact norm_add_le _ _
          rw [hrel.1] at hrec ⊢
          linarith
        have hcert' : ∀ L ∈ L₁ :: rest', LinkCertified u L := fun L hL =>
          hcert L (List.mem_cons_of_mem _ hL)
        have hd' : ∀ La ∈ (L₁ :: rest').head?,
            ‖u La.tStart - La.centre La.tStart‖ ≤ La.δ := by
          intro La hLa
          have h' : some L₁ = some La := hLa
          obtain rfl := Option.some.inj h'
          exact hd₁
        exact ih hcert' hchain' hd' M hM'

/-- **C3, union form.**  Every time instant covered by some slab of a
certified chain admits a link whose tube bounds `u` there — the "certified
interval is the union of the slabs" statement, with the certified radius being
the per-slab tube. -/
theorem chain_composition_union (u : ℝ → E) (links : List (ChainLink E))
    (hcert : ∀ L ∈ links, LinkCertified u L)
    (hchain : links.IsChain LinkComposable)
    (hd : ∀ L₀ ∈ links.head?, ‖u L₀.tStart - L₀.centre L₀.tStart‖ ≤ L₀.δ) :
    ∀ t : ℝ, (∃ L ∈ links, t ∈ Icc L.tStart L.tEnd) →
      ∃ L ∈ links, t ∈ Icc L.tStart L.tEnd ∧ ‖u t - L.centre t‖ ≤ L.radius t := by
  intro t ⟨L, hL, ht⟩
  exact ⟨L, hL, ht, chain_composition u links hcert hchain hd L hL t ht⟩

end Chain

/-! ## C4 — the discrete Grönwall inequality -/

section DiscreteGronwall

/-- **C4, discrete Grönwall.**  If each step contracts by `A ≥ 0` and injects
at most `B`, then after `n` steps
`x n ≤ Aⁿ x 0 + B · Σ_{i<n} Aⁱ`.  This is the a priori shape of the chained
datum radii `δ_n` when every slab certifies the same contraction factor and
per-slab transfer budget; the Python chain measures the actual per-slab values
instead, and this inequality is the reason a uniformly contracting chain has
uniformly bounded deltas. -/
theorem discrete_gronwall {x : ℕ → ℝ} {A B : ℝ} (hA : 0 ≤ A)
    (h : ∀ n, x (n + 1) ≤ A * x n + B) (n : ℕ) :
    x n ≤ A ^ n * x 0 + B * ∑ i ∈ Finset.range n, A ^ i := by
  induction n with
  | zero => simp
  | succ n ih =>
    have hmul : A * x n ≤ A * (A ^ n * x 0 + B * ∑ i ∈ Finset.range n, A ^ i) :=
      mul_le_mul_of_nonneg_left ih hA
    calc x (n + 1) ≤ A * x n + B := h n
      _ ≤ A * (A ^ n * x 0 + B * ∑ i ∈ Finset.range n, A ^ i) + B := by linarith
      _ = A ^ (n + 1) * x 0 + B * ∑ i ∈ Finset.range (n + 1), A ^ i := by
          rw [geom_sum_succ, pow_succ]; ring

end DiscreteGronwall

/-! ## C5 — the chained bound never exceeds the per-slab maximum -/

section RadiusMax

variable {E : Type*} [NormedAddCommGroup E]

/-- **C5, pointwise.**  The two-slab piecewise radius is at most the pointwise
maximum of the branch radii. -/
theorem piecewise_radius_le_max {R₁ R₂ : ℝ → ℝ} (t₁ t : ℝ) :
    (if t ≤ t₁ then R₁ t else R₂ t) ≤ max (R₁ t) (R₂ t) := by
  split
  · exact le_max_left _ _
  · exact le_max_right _ _

/-- Every member of a list is at most the `foldr max` of the list over any
base value. -/
theorem le_foldr_max {l : List ℝ} {x : ℝ} (b : ℝ) (hx : x ∈ l) :
    x ≤ l.foldr max b := by
  induction l with
  | nil => simp at hx
  | cons a l ih =>
    rcases List.mem_cons.1 hx with rfl | hx'
    · exact le_max_left _ _
    · exact le_trans (ih hx') (le_max_right _ _)

/-- **C5, chain form.**  If each link's radius is bounded on its slab by a
per-link constant `M L`, then the chained tube is bounded everywhere by the
maximum of those constants (as a `foldr max 0`; the base `0` is free because
the left side is a norm).  This is the statement that chaining never
manufactures a radius larger than the worst slab. -/
theorem chain_radius_le_foldr_max {u : ℝ → E} {links : List (ChainLink E)}
    {M : ChainLink E → ℝ}
    (htube : ∀ L ∈ links, ∀ t ∈ Icc L.tStart L.tEnd, ‖u t - L.centre t‖ ≤ L.radius t)
    (hM : ∀ L ∈ links, ∀ t ∈ Icc L.tStart L.tEnd, L.radius t ≤ M L) :
    ∀ L ∈ links, ∀ t ∈ Icc L.tStart L.tEnd,
      ‖u t - L.centre t‖ ≤ (links.map M).foldr max 0 := by
  intro L hL t ht
  exact ((htube L hL t ht).trans (hM L hL t ht)).trans
    (le_foldr_max 0 (List.mem_map.2 ⟨L, hL, rfl⟩))

end RadiusMax

/-! ## C6 — the Lagrange endpoint bound -/

section TaylorEndpoint

/-- **C6, the Lagrange endpoint bound.**  For a real function that is
`C^{m+1}` on the slab `[t₀, t₀ + h]`, the distance between the true endpoint
value and the degree-`m` Taylor evaluation from `t₀` is at most
`M · h^{m+1} / (m+1)!` whenever `M` dominates the `(m+1)`-st derivative on the
open slab.  This is the remainder inequality behind the recentring transfer:
the Python chain evaluates the Taylor endpoint of the Galerkin trajectory and
charges exactly this bound (computed by interval arithmetic on the analytic
coefficient recursion) into `δ_{n+1}`.  Specialisation of mathlib's
`taylor_mean_remainder_lagrange_iteratedDeriv`. -/
theorem taylor_endpoint_remainder_bound {f : ℝ → ℝ} {t₀ h M : ℝ} {m : ℕ}
    (hh : 0 < h)
    (hf : ContDiffOn ℝ (m + 1) f (Icc t₀ (t₀ + h)))
    (hM : ∀ s ∈ Ioo t₀ (t₀ + h), |iteratedDeriv (m + 1) f s| ≤ M) :
    |f (t₀ + h) - taylorWithinEval f m (Icc t₀ (t₀ + h)) t₀ (t₀ + h)|
      ≤ M * h ^ (m + 1) / (m + 1)! := by
  have hle : t₀ ≤ t₀ + h := by linarith
  have hne : t₀ ≠ t₀ + h := ne_of_lt (by linarith)
  have huIcc : uIcc t₀ (t₀ + h) = Icc t₀ (t₀ + h) := uIcc_of_le hle
  have huIoo : uIoo t₀ (t₀ + h) = Ioo t₀ (t₀ + h) := uIoo_of_le hle
  obtain ⟨x', hx', heq⟩ := taylor_mean_remainder_lagrange_iteratedDeriv hne
    (by rw [huIcc]; exact hf)
  rw [huIcc] at heq
  rw [huIoo] at hx'
  rw [heq]
  have hsub : t₀ + h - t₀ = h := by ring
  rw [hsub, abs_div, abs_mul, abs_pow, abs_of_pos hh, Nat.abs_cast]
  have hpow : (0 : ℝ) ≤ h ^ (m + 1) := pow_nonneg hh.le _
  have hnum : |iteratedDeriv (m + 1) f x'| * h ^ (m + 1) ≤ M * h ^ (m + 1) :=
    mul_le_mul_of_nonneg_right (hM x' hx') hpow
  have hinv : (0 : ℝ) ≤ (((m + 1)! : ℝ))⁻¹ := by positivity
  calc |iteratedDeriv (m + 1) f x'| * h ^ (m + 1) / ((m + 1)! : ℝ)
      = |iteratedDeriv (m + 1) f x'| * h ^ (m + 1) * (((m + 1)! : ℝ))⁻¹ := by
        rw [div_eq_mul_inv]
    _ ≤ M * h ^ (m + 1) * (((m + 1)! : ℝ))⁻¹ :=
        mul_le_mul_of_nonneg_right hnum hinv
    _ = M * h ^ (m + 1) / ((m + 1)! : ℝ) := by rw [div_eq_mul_inv]

end TaylorEndpoint

end NSSingularity
