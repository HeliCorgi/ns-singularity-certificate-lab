# The EXT-P1/P2/P3 ledger: exact statements, dependencies, and an unaudited Galerkin proof of EXT-P1

**Status: a dependency document, plus one complete — but UNAUDITED — paper
proof.**  Track P and the Track-P chain issue a *conditional PDE certificate
assuming EXT-P1/P2/P3*; the finite-dimensional Galerkin enclosure and the
scalar control ODE are verified unconditionally, and nothing else is.  This
note does three things and refuses a fourth.  It expands the three payload
one-liners of `torus_aposteriori.EXTERNAL_THEOREMS` into fully quantified
statements (§1); it maps every ingredient of those statements onto what exists
today in mathlib, on paper, and in this repository's Lean tree (§2); and it
writes out a complete paper proof of EXT-P1 by the Galerkin method, every
inequality displayed (§3).  What it does **not** do: change a single payload
flag.  EXT-P1/P2/P3 keep `proved: false` and `axiomatised_in_lean: false` in
every payload, the checkers keep rejecting anything else, and §4 records why an
audit of §3 would still not move those flags.  Nothing here is a singularity
result, and nothing here bears on the Clay problem.

Notation is that of [`track_p_periodic.md`](track_p_periodic.md), fixed once:
`T³ = (ℝ/2πℤ)³` with the normalised measure `(2π)⁻³dx`; fields are real,
mean-zero, divergence-free; Fourier expansion `u = Σ_{k≠0} û_k e^{ik·x}` with
`û_{−k} = conj(û_k)`, `k·û_k = 0`; Parseval `⟨f,g⟩ = Σ_k f̂_k·conj(ĝ_k)`;
homogeneous norms `‖f‖_n² = Σ_k |k|^{2n}|f̂_k|²`, monotone in `n` on mean-zero
fields because `|k| ≥ 1`.  The inhomogeneous `H⁴` norm is equivalent on this
class: `‖f‖₄ ≤ ‖f‖_{H⁴} ≤ 4‖f‖₄`, since `(1+|k|²)⁴ ≤ (2|k|²)⁴ = 16|k|⁸` for
`|k| ≥ 1`.  `A` is the lattice constant of
[`track_p_periodic.md`](track_p_periodic.md) §5.1: `A² = Σ_{k≠0}|k|⁻⁴ ≤
16.9245`, `A ≤ 4.114`, and `‖f‖_∞ ≤ A‖f‖₂` for mean-zero `f`.
`H⁴_σ` abbreviates "mean-zero divergence-free `H⁴(T³;ℝ³)`".

## 1. Exact statements

The one-line forms in `EXTERNAL_THEOREMS` are the payload versions and they do
not change; the checker (`verify_torus_slab_certificate`) compares them
verbatim.  The expansions below are conservative refinements: they quantify
the same claims, in the same function spaces, with nothing added that the
named classical sources do not supply.  Sources, with the repository's
recall-honesty convention: Kato's quasi-linear theory; Temam, *Navier–Stokes
Equations*, Ch. III (periodic `H^m` theory, `m ≥ 2`; the theorem *number* is
[recalled], not re-checked); Majda–Bertozzi, *Vorticity and Incompressible
Flow*, Ch. 3 (`H^m`, `m > 5/2`, covering `m = 4`).

**Definition (strong solution).**  Let `0 ≤ t₀ < t₁`.  A *strong solution* on
`[t₀, t₁]` is a mean-zero divergence-free
`u ∈ C([t₀,t₁]; H⁴_σ) ∩ C¹([t₀,t₁]; H²)` such that

```
∂_t u = ν Δu − P(u·∇u)      in H², pointwise in t ∈ [t₀,t₁] ,
```

where `P` is the Leray projector (`P_k = I − kkᵀ/|k|²` mode by mode).  For
such `u` the pressure is recovered by `∇p = −(I−P)(u·∇u)`, i.e.
`p̂_k = i k·(u·∇u)^̂_k / |k|²`, `p̂₀ = 0`, and the pair `(u, p)` satisfies
`∂_t u + u·∇u + ∇p = νΔu` classically in space (`H⁴ ⊂ C²` by the `A`
embedding).  Equivalent characterisation used in §3: `u ∈ C([t₀,t₁]; H⁴_σ)`
and every Fourier coefficient satisfies the integral equation (3.d.1) below.

### 1.1 EXT-P1, expanded

> **EXT-P1 (local existence and uniqueness).**
> *Hypotheses.*  `ν > 0` fixed; `u₀ ∈ H⁴(T³;ℝ³)` with `∫_{T³} u₀ = 0` and
> `div u₀ = 0` in the distributional sense (equivalently `k·û₀ₖ = 0` for all
> `k` and `û₀,₀ = 0`).
> *Existence.*  There exist `T > 0`, depending only on `ν` and
> `‖u₀‖_{H⁴}` — the §3 proof gives the explicit admissible value
> `T = T* := 1/(270 A ‖u₀‖₄) ≥ 1/(1110.78·‖u₀‖_{H⁴})` — and a strong
> solution `u` on `[0, T]` with `u(0) = u₀`, together with the recovered
> pressure `p ∈ C([0,T]; H⁴)` (normalised mean-zero).
> *Uniqueness clause.*  If `u` and `v` are strong solutions on `[t₀, t₁]`
> (same `ν`) with `u(t₀) = v(t₀)`, then `u(t) = v(t)` for every
> `t ∈ [t₀, t₁]`.  The equation is autonomous, so the clause holds on every
> time interval, not only ones containing `0`.
> *Maximal solution.*  By the uniqueness clause, all strong solutions with
> datum `u₀` patch to a unique maximal one on `[0, T_max)`,
> `T_max ≥ T* > 0`.

**How the chain certificate uses the uniqueness clause.**  The slab chain
(`torus_chain.py`; the finite-inequality skeleton of the composition is
`formal/NSSingularity/TrackPChain.lean`) certifies, conditionally, a strong
solution on
`[t_n, t_{n+1}]` for each slab.  At the boundary `t_{n+1}` the next slab's
certificate speaks about the strong solution with datum `u(t_{n+1})` — an
object supplied by EXT-P1 applied at time `t_{n+1}`.  The uniqueness clause,
in its time-shifted form, is what identifies that solution with the
continuation of the slab-`n` solution: both are strong solutions on
`[t_{n+1}, t_{n+2}]` agreeing at `t_{n+1}`, so they coincide, and "the
certified solution on the union of the slabs" is well defined.  Without the
uniqueness half of EXT-P1 the chain would certify a *family* of per-slab
solutions with no license to call them one solution.

### 1.2 EXT-P2, expanded

> **EXT-P2 (regularity to run the estimate).**
> *Hypotheses.*  `u` the strong solution of EXT-P1 on `[0,T]`; `u_a` a
> Galerkin trajectory (a finite trigonometric polynomial in `x` with
> real-analytic-in-`t` coefficients, mean-zero, divergence-free);
> `w = u − u_a`.
> *Conclusions.*
> (i) `u ∈ C([0,T]; H⁴) ∩ C¹([0,T]; H²)`; moreover (standard, and delivered
> by the §3 construction) `u ∈ L²(0,T; H⁵)` and `∂_t u ∈ L²(0,T; H³)`.
> (ii) The `H⁴` energy estimate for `w` is justified: `t ↦ ‖w(t)‖₄²` is
> absolutely continuous with `(d/dt)‖w‖₄² = 2⟨∂_t w, w⟩₄` for a.e. `t`, the
> pairing being the absolutely convergent series
> `2Σ_k |k|⁸ Re(ŵ_k' · conj(ŵ_k))` (finite for a.e. `t` by (i)).
> (iii) `t ↦ ‖w(t)‖₄` admits at every `t` a right Dini derivative bounded by
> the §5.4 estimate of [`track_p_periodic.md`](track_p_periodic.md),
> `d⁺W/dt ≤ (−ν + 9(K₁+K₂))W + 135A·W² + ‖e(t)‖₄`, so the Chaplygin/Dini
> comparison with the control ODE applies and yields `W ≤ R` on the slab.

**Honesty note (a deliberate retreat, recorded).**  The §3 machinery of this
note delivers (i) and (ii) for the solution it constructs, and (ii) already
suffices to run the comparison in integral form (an absolutely continuous
`W²` with an a.e. differential inequality admits the same comparison
conclusion).  The *pointwise right-Dini* form (iii), which is the wording the
payload carries, is the classical statement and needs the standard
mollification argument (Temam Ch. III; Majda–Bertozzi §3.2) on top of (ii);
that upgrade is **not** reproduced in this note.  EXT-P2 therefore remains an
external theorem in full, and no part of §3 is offered as a proof of it.

### 1.3 EXT-P3, expanded

> **EXT-P3 (H⁴ continuation criterion).**
> *Hypotheses.*  `u` the maximal strong solution of EXT-P1 on `[0, T_max)`.
> *Conclusion.*  If `T ≤ T_max` and `sup_{t < T} ‖u(t)‖_{H⁴} =: M < ∞`,
> then `T_max > T`.  Equivalently: if `T_max < ∞` then
> `limsup_{t ↑ T_max} ‖u(t)‖_{H⁴} = ∞`.

**Remark (route to a proof, not offered as one).**  Given the quantitative
EXT-P1 of §3 — local time bounded below by `1/(270A·M')` whenever the datum
norm is at most `M'` — EXT-P3 follows by a three-line restart: pick
`t₀ < T` with `T − t₀ < 1/(540A·M)`, apply EXT-P1 at `t₀` with datum
`u(t₀)` (norm `≤ M`), and glue with the uniqueness clause; the glued solution
is strong past `T`.  This note does **not** promote that remark to a proved
status: §3's mandate is EXT-P1 alone, and EXT-P3 stays external with
`proved: false` like the others.

**How the slab certificate uses EXT-P3.**  On a certified slab,
`‖u(t)‖₄ ≤ ‖u_a(t)‖₄ + R(t)`, with `‖u_a‖₄` bounded by the Picard box and
`R` by the control tube; both are finite on the whole slab, so EXT-P3 forbids
`T_max` from landing inside the slab and the existence interval reaches `h`.

## 2. Dependency table

One row per ingredient.  "Mathlib today" is checked against the vendored tree
under `formal/.lake/packages/mathlib` (this session), not against memory.
"Finite Fourier sums" asks whether every object in the ingredient is a finite
trigonometric sum / finite-dimensional vector (so that the repository's exact
rational style applies instance by instance); "genuinely ∞-dim" marks the
steps that intrinsically live in the infinite-dimensional function space.

**Honest baseline on mathlib.**  Mathlib has Picard–Lindelöf ODE theory
(`Mathlib/Analysis/ODE/PicardLindelof.lean`: structure `IsPicardLindelof`,
theorem `IsPicardLindelof.exists_eq_forall_mem_Icc_eq_picard`; uniqueness in
`Mathlib/Analysis/ODE/ExistUnique.lean`: `ODE_solution_unique_of_mem_Icc` and
variants), linear Grönwall comparison
(`Mathlib/Analysis/ODE/Gronwall.lean`: `gronwallBound`,
`norm_le_gronwallBound_of_norm_deriv_right_le`,
`le_gronwallBound_of_liminf_deriv_right_le`), the Ascoli theorem
(`Mathlib/Topology/UniformSpace/Ascoli.lean`), Banach–Alaoglu
(`Mathlib/Analysis/Normed/Module/WeakDual.lean`), `p`-series summability
(`Mathlib/Analysis/PSeries.lean`) and the `tsum`/`lp` infrastructure.  It has
**no** `H^s(T³)` Sobolev spaces, no Leray projector on function spaces, no
vector-valued `L²(0,T;H^s)` parabolic framework, and no Navier–Stokes
anything; Fourier analysis on the torus is essentially one-dimensional
(`AddCircle`).  Lane L's `formal/NSSingularity/GalerkinPicard.lean` is the
instantiation of that ODE theory for a quadratic coefficient field: it landed
in the working tree this turn with `quadraticField`,
`quadratic_field_lipschitzOnWith`, `quadratic_ode_local_solution` (the lemma
this note cross-references), `quadratic_ode_local_solution_hasDerivAt` and
`quadratic_ode_unique`, declaring no `sorry`, `admit` or axioms; this session
did not run a Lean build (the repository verifies builds in CI), so rows
citing it say "in tree" rather than "verified built".

| # | ingredient (step in §3) | available in mathlib today | provable on paper | Lean-addable near-term | needs only finite Fourier sums | genuinely ∞-dim |
|---|---|---|---|---|---|---|
| 1 | (a) Picard–Lindelöf for the Galerkin coefficient ODE | yes: `IsPicardLindelof`, `exists_eq_forall_mem_Icc_eq_picard`, `ODE_solution_unique_of_mem_Icc` | §3(a); classical | yes — `GalerkinPicard.lean` (`quadratic_ode_local_solution`, `quadratic_ode_unique`; lane L, in tree, CI build pending) instantiates it for the quadratic field | yes | no |
| 2 | (a) `L²` energy identity and global extension at fixed `N` | abstract core already in this repo: `GalerkinNoBlowup.lean` (`EnergyNeutral`, `galerkin_norm_le`, `galerkin_not_tendsto_atTop`); mathlib has no packaged maximal-interval/escape API | §3(a) | yes — glue existing repo lemmas to the ODE layer | yes | no |
| 3 | (b) `H⁴` self-estimate algebra (Leibniz, Vandermonde, `9 = √81`, `135A`) | nothing (no `H^s(T³)`) | [`track_p_periodic.md`](track_p_periodic.md) §5.2/§5.4, instantiated in §3(b) | plausible but heavy: per-`N` it is a finite rational inequality; bricks exist in `TrackPFourier.lean` (`leray`, `slotDivergence_cosMode`, `contDiff_trigPolynomial`, `weighted_sum_pow_mono`) | yes (each instance; `∀N` is a schema over finite sums) | no |
| 4 | (b) lattice constant `A² = Σ_{k≠0}|k|⁻⁴ ≤ 16.9245` | summability tools: `Mathlib/Analysis/PSeries.lean`, `tsum` | [`track_p_periodic.md`](track_p_periodic.md) §5.1 (exact partial sum + `26/N` tail) | yes | no — one infinite lattice series with an elementary tail majorant | no (a single scalar series) |
| 5 | (b) scalar Riccati bound `Y ≤ 4Y₀` on `[0,T*]` via `φ = Y^{−1/2}` | FTC/MVT yes; nonlinear scalar comparison not packaged; related scalar machinery in repo `ControlODE.lean` | §3(b), complete | yes — scalar real analysis only | scalar (no Fourier at all) | no |
| 6 | (c) per-coefficient Lipschitz bounds + Arzelà–Ascoli + diagonal extraction | Ascoli: `Mathlib/Topology/UniformSpace/Ascoli.lean`; diagonal argument standard | §3(c), complete | plausible | no — countably many modes at once | yes (first ∞-dim step) |
| 7 | (c) uniform `H⁴` tail bound → `C([0,T*];H³)` convergence; Fatou → `L^∞H⁴`, `L²Ḣ⁵` | `tsum`/`liminf` tools only | §3(c), complete | with effort | no | yes |
| 8 | (d) lattice bilinear estimates (`ℓ¹∗ℓ²` Young; `‖v·∇w‖₃ ≤ 8A‖v‖₄‖w‖₄`) | `lp` machinery; no `Ḣ^s(T³)` product calculus | §3(d), complete | possible as pure `ℓ²(ℤ³)` statements | no — infinite convolutions | yes |
| 9 | (d) limit passage into the coefficient ODEs; norm-continuity (absolute continuity of `‖u‖₄²`) → `u ∈ C H⁴`; pressure recovery | FTC and dominated convergence available; no Lions–Magenes-type lemma; per-mode Leray is `TrackPFourier.leray` | §3(d), complete (coefficientwise, self-contained) | per-mode parts yes; the summation/interchange steps with effort | per-mode yes; summation no | yes |
| 10 | (e) `L²` Grönwall uniqueness (trilinear cancellation + linear Grönwall) | Grönwall half: `norm_le_gronwallBound_of_norm_deriv_right_le` is the right shape; torus integration by parts for `⟨u·∇w,w⟩ = 0`: absent | §3(e), complete | Grönwall glue yes; the cancellation needs the function-space layer | no | yes |
| 11 | (P2-i) `u ∈ C¹([0,T];H²)`, `L²H⁵`, `∂_t u ∈ L²H³` | — | delivered for the §3 solution (rows 7, 9) | with row-9 effort | no | yes |
| 12 | (P2-ii/iii) AC energy identity for `w = u − u_a`; pointwise Dini upgrade; Chaplygin comparison | comparison-adjacent lemmas only (`Gronwall.lean`) | AC form: same class as row 9.  Pointwise Dini: **open here** (classical mollification, Temam III / Majda–Bertozzi §3.2).  Chaplygin: scalar | scalar comparison yes; Dini upgrade no | no | yes |
| 13 | (P3) quantitative local time `T*(M) = 1/(270A·M)` | — | §3(b) corollary | tied to rows 3–5 | mixed | no |
| 14 | (P3) restart + gluing via uniqueness → continuation | — | three-line remark in §1.3, **not audited** | needs rows 9–10 first | no | yes |

Reading the table honestly: every ingredient that is finite-dimensional
(rows 1–5, 13) is either already covered by repository Lean files, planned in
lane L, or a bounded amount of rational algebra; every genuinely
infinite-dimensional row (6–12, 14) is classical on paper but has **no**
mathlib substrate to build on today, because the function spaces themselves
are missing.  That asymmetry is the entire reason EXT-P1/P2/P3 are carried as
external theorems rather than formalisation targets of the current milestone.

## 3. A complete paper proof of EXT-P1 by the Galerkin method — UNAUDITED

Everything in this section is classical mathematics written out with explicit
constants; nothing in it has been machine-checked, and §3.6 records its
status.  Fix `ν > 0` and `u₀ ∈ H⁴_σ`, write `Y₀ := ‖u₀‖₄²`.  If `u₀ = 0`
the solution `u ≡ 0` satisfies everything trivially, so assume `Y₀ > 0`.

**Theorem 1 (quantitative EXT-P1).**  Let

```
T* := 1 / (270 A √Y₀) = 1 / (270 A ‖u₀‖₄) ,       A ≤ 4.114 .
```

Then there is a strong solution `u` on `[0, T*]` (in the sense of §1) with
`u(0) = u₀`,

```
sup_{[0,T*]} ‖u(t)‖₄ ≤ 2‖u₀‖₄ ,      ∫₀^{T*} ‖u(t)‖₅² dt ≤ 9Y₀/(2ν) ,
```

and pressure `p ∈ C([0,T*]; H⁴)`.  Moreover any two strong solutions on any
`[t₀,t₁]` that agree at `t₀` coincide on `[t₀,t₁]` (§3(e), with the explicit
Grönwall constant).  Since `‖u₀‖₄ ≤ ‖u₀‖_{H⁴}`, the time `T*` admits the
lower bound `1/(270A‖u₀‖_{H⁴})` stated in §1.1.

*Caution.*  `T*` is the local-existence time of EXT-P1 and nothing else.  The
slab lengths of the certificates are chosen by the Picard box and the control
ODE, not by `T*`; the two must not be conflated.

### 3(a) The Galerkin system: local Picard, then global at fixed N

For `N ≥ 1` let `G_N = {k ∈ ℤ³ : 1 ≤ |k|² ≤ N²}` (the band of
`galerkin_modes(N²)`) and let `P_N` be the sharp Fourier cutoff onto
`span{e^{ik·x} : k ∈ G_N}` — self-adjoint in every `⟨·,·⟩_n` (a real
diagonal multiplier) and commuting with derivatives and with `P`.  The
Galerkin system is

```
∂_t u_N = ν Δ u_N − P_N P (u_N·∇u_N) ,      u_N(0) = P_N u₀ .          (3.a.1)
```

In coefficients: for `k ∈ G_N`,

```
d/dt û_k = −ν|k|² û_k − i P_k B_k(u_N,u_N) ,
B_k(v,w) := Σ_{m+n=k} (v̂_m·n) ŵ_n            (finite sum over G_N here). (3.a.2)
```

The right-hand side is a polynomial (linear plus quadratic) in the finitely
many real coefficients, hence `C^∞` and locally Lipschitz on every ball, so
Picard–Lindelöf gives a unique local solution through any datum.  This is
row 1 of the table: mathlib's `IsPicardLindelof` +
`ODE_solution_unique_of_mem_Icc`, instantiated for abstract quadratic fields
by `GalerkinPicard.lean` (`quadratic_ode_local_solution`, in tree, CI build
pending); the repository's
rational-arithmetic engine proves the same statement *with explicit boxes* for
specific data (`prove_galerkin_box`, the `PICARD_SELF_MAPPING` theorem of
`control_ode.py` applied coefficient-wise).

**Lemma 3.a.0 (constraint invariance along the Galerkin flow) [turn-11 audit
erratum — repairs the audited gap that `div u_N(t) = 0` and reality of
`u_N(t)` were used for `t > 0` but proved only at `t = 0`; see
`ext_p1_p2_p3_audit.md` §3, item G1].**  For every `N` and every `t` in the
Galerkin solution's interval of existence, `k·û_{N,k}(t) = 0` and
`û_{N,−k}(t) = conj(û_{N,k}(t))` for all `k ∈ G_N`; in particular `u_N(t)`
stays real, mean-zero and divergence-free, and `P u_N(t) = u_N(t)`, for as
long as it exists.

*Proof.*  *Divergence.*  `k·P_k v = 0` for every `v ∈ ℂ³` (`P_k` projects
onto `k^⊥`), so `s_k := k·û_{N,k}` satisfies, by (3.a.2), the **linear**
scalar ODE `s_k' = −ν|k|² s_k` with datum
`s_k(0) = k·(P_N u₀)^̂_k = k·û₀ₖ = 0`; by uniqueness for locally Lipschitz
ODEs (the row-1 shape: `ODE_solution_unique_of_mem_Icc`, instantiated in
`GalerkinPicard.lean` as `quadratic_ode_unique`), `s_k ≡ 0`.  *Reality.*  The
conjugation symmetry `Θ : (v_k)_{k∈G_N} ↦ (conj(v_{−k}))_{k∈G_N}` maps
solutions of (3.a.2) to solutions: conjugating the `−k` equation and
substituting `m ↦ −m`, `n ↦ −n` in the convolution gives

```
conj(B_{−k}(v,v)) = Σ_{m+n=k} (conj(v_{−m})·(−n)) conj(v_{−n}) = −B_k(Θv, Θv) ,
```

and `P_{−k} = P_k` is a real matrix, so
`d/dt (Θv)_k = conj(−ν|k|² v_{−k} − i P_{−k} B_{−k}(v,v))
= −ν|k|² (Θv)_k − i P_k B_k(Θv, Θv)` — the `k` equation again.  The datum
`P_N u₀` is the coefficient vector of a real field, hence `Θ`-fixed, so
Picard uniqueness gives `Θ ∘ û_N = û_N` on the whole interval, i.e.
`û_{N,−k} = conj(û_{N,k})`.  ∎

*Global extension at fixed `N`.*  Pairing (3.a.1) with `u_N` in `⟨·,·⟩₀`:
`P_N` and `P` drop (self-adjoint, `P_N u_N = u_N`, `P u_N = u_N`), and the
exact trilinear identity `⟨u_N·∇u_N, u_N⟩₀ = ½∫u_N·∇|u_N|² = −½∫(div u_N)|u_N|² = 0`
holds — an identity of finite trigonometric polynomials, whose hypotheses
`div u_N(t) = 0` and reality of `u_N(t)` are supplied for every `t` by
Lemma 3.a.0 [turn-11 audit erratum].  Hence

```
d/dt ‖u_N‖₀² = −2ν‖u_N‖₁² ≤ 0 ,      ‖u_N(t)‖₀ ≤ ‖P_N u₀‖₀ ≤ ‖u₀‖₀ .
```

On the band, `‖u_N‖₄ ≤ N⁴‖u_N‖₀`, so the coefficient vector stays in a fixed
compact ball; a maximal ODE solution with locally Lipschitz field either is
global or leaves every compact set, so the Galerkin solution is global in
time.  (The abstract core of exactly this argument is already formalised in
`formal/NSSingularity/GalerkinNoBlowup.lean`: `EnergyNeutral`,
`galerkin_norm_le`, `galerkin_not_tendsto_atTop` — row 2.)  All norms below
are therefore defined for all `t ≥ 0`, and `t ↦ û_k(t)` is real-analytic; in
particular `Y_N(t) := ‖u_N(t)‖₄²` is `C¹`, and **no Dini derivative is needed
anywhere at the Galerkin level**.

### 3(b) The uniform-in-N a priori bound, with the §5.4 algebra

Pair (3.a.1) with `u_N` in `⟨·,·⟩₄`.  As in (a), `P_N` and `P` drop at no
cost, and the viscous term gives exactly `−ν‖u_N‖₅²`:

```
½ d/dt Y_N = −ν‖u_N‖₅² − ⟨u_N·∇u_N, u_N⟩₄ .                          (3.b.1)
```

The nonlinear pairing is the **cubic term** of
[`track_p_periodic.md`](track_p_periodic.md) §5.4, instantiated with
`w := u_N` — the same commutator algebra, verbatim; nothing new is derived
here.  For the record, the instantiated chain of inequalities: expand
`⟨u_N·∇u_N, u_N⟩₄ = Σ_{|α|=4} c_α ⟨∂^α(u_N·∇u_N), ∂^α u_N⟩` with
`c_α = 4!/α!` (§5.2); Leibniz over `β ≤ α`; the `β = 0` term pairs to zero
because `div u_N = 0` (Lemma 3.a.0 [turn-11 audit erratum]); each surviving
`|β| = j` term is bounded, using only
the `A` embedding of §5.1, by

```
j = 1, 2 :  ‖∂^β u_N‖_∞ · ‖∇∂^{α−β}u_N‖ · ‖∂^α u_N‖
            ≤ (A‖u_N‖_{j+2}) · ‖u_N‖_{5−j} · ‖∂^α u_N‖ ≤ A‖u_N‖₄² ‖∂^α u_N‖ ,
j = 3, 4 :  ‖|∇∂^{α−β}u_N|_F‖_∞ · ‖∂^β u_N‖ · ‖∂^α u_N‖
            ≤ (A‖u_N‖_{7−j}) · ‖u_N‖_j · ‖∂^α u_N‖ ≤ A‖u_N‖₄² ‖∂^α u_N‖ ;
```

the Vandermonde identity collapses the `β`-sum to the weights `binom(4,j)`
with `Σ_{j=1}^4 binom(4,j) = 15`, and Cauchy–Schwarz over `α` with the
weights `c_α` supplies `9 = √81 = (Σ_α c_α)^{1/2}`:

```
|⟨u_N·∇u_N, u_N⟩₄| ≤ 9·15·A ‖u_N‖₄³ = 135 A ‖u_N‖₄³ .               (3.b.2)
```

Where is `9(K₁+K₂)`?  Nowhere in this self-estimate — and that is a statement
about the algebra, not a discrepancy.  The §5.4 machinery has two
instantiations: when one Leibniz factor is a *known band-limited field* the
sup is taken from the exact rational `ℓ¹` majorants `M_j, N_m` and the
transport/stretching constants `9K₁, 9K₂` come out; when every factor is the
*unknown* field the sup is taken from the `A` embedding and the cubic
constant `135A` comes out.  The certificate's difference estimate
(`w = u − u_a`) uses both; the a priori self-estimate uses only the second.
Same Leibniz expansion, same Vandermonde collapse, same `√81` — one algebra,
two instantiations.

From (3.b.1)–(3.b.2), with `Y_N = ‖u_N‖₄²`:

```
d/dt Y_N ≤ −2ν‖u_N‖₅² + 270 A Y_N^{3/2}
         ≤ −2ν Y_N + 270 A Y_N^{3/2} .                                (3.b.3)
```

So the form that actually comes out is `Y' ≤ C₁ Y^{3/2} + C₂ Y` with
`C₁ = 270A ≤ 1110.78` and `C₂ = −2ν ≤ 0`; the linear term *helps* and is
discarded for the uniform bound: `Y_N' ≤ C₁ Y_N^{3/2}`, with
`Y_N(0) = ‖P_N u₀‖₄² ≤ Y₀` (sharp cutoffs contract every `‖·‖_n`).

*The uniform bound, by hand.*  Claim: `Y_N(t) ≤ 4Y₀` for all
`t ∈ [0, T*]`, `T* = 1/(C₁√Y₀)`, for every `N`.  Suppose not:
`Y_N(t₁) > 4Y₀` for some `t₁ ≤ T*`.  Let
`t₀ := sup{t ≤ t₁ : Y_N(t) ≤ Y₀}`; the set contains `0`, and by continuity
`Y_N(t₀) = Y₀` and `Y_N > Y₀ > 0` on `(t₀, t₁]`.  On `[t₀, t₁]` the function
`φ := Y_N^{−1/2}` is `C¹` (positivity) with

```
φ' = −½ Y_N^{−3/2} Y_N' ≥ −½ C₁ ,
```

so `φ(t₁) ≥ φ(t₀) − ½C₁(t₁−t₀) ≥ Y₀^{−1/2} − ½C₁T* = Y₀^{−1/2}(1 − ½) = ½Y₀^{−1/2}`,
i.e. `Y_N(t₁) ≤ 4Y₀` — contradiction.  Hence

```
sup_{[0,T*]} ‖u_N(t)‖₄ ≤ 2√Y₀      for every N ,                      (3.b.4)
```

with `T* = 1/(270A√Y₀)` **independent of `N`**, depending on the datum only
through `‖u₀‖₄`.  (This scalar argument is row 5: nothing but the fundamental
theorem of calculus on a positive `C¹` function.)

*Corollary (uniform `Ḣ⁵` budget).*  Integrating the first line of (3.b.3)
over `[0,T*]` and using (3.b.4):

```
2ν ∫₀^{T*} ‖u_N‖₅² dt ≤ Y_N(0) + 270A·T*·(4Y₀)^{3/2} = Y₀ + 8Y₀ = 9Y₀ ,
∫₀^{T*} ‖u_N‖₅² dt ≤ 9Y₀/(2ν) .                                       (3.b.5)
```

### 3(c) Compactness: coefficientwise Arzelà–Ascoli, diagonal, tail upgrade

*The per-coefficient Lipschitz bound, written out.*  For `k ∈ G_N`,
Cauchy–Schwarz over the free convolution index in (3.a.2) gives
`|B_k(u_N,u_N)| ≤ ‖u_N‖₀‖u_N‖₁ ≤ Y_N ≤ 4Y₀` on `[0,T*]`, and
`|û_k| ≤ ‖u_N‖₄/|k|⁴ ≤ 2√Y₀/|k|⁴`; `|P_k v| ≤ |v|`.  Hence, for every
`N ≥ |k|` and `t ∈ [0,T*]`,

```
|d/dt û_{N,k}(t)| ≤ ν|k|²·(2√Y₀/|k|⁴) + 4Y₀ = 2ν√Y₀/|k|² + 4Y₀ =: L_k
                  ≤ 2ν√Y₀ + 4Y₀ .                                     (3.c.1)
```

*Extraction.*  For fixed `k`, the family `{t ↦ û_{N,k}(t)}_{N ≥ |k|}` is
uniformly bounded (by `2√Y₀`) and `L_k`-Lipschitz on `[0,T*]`, hence
equicontinuous; Arzelà–Ascoli yields a uniformly convergent subsequence.
Enumerate `ℤ³∖{0}` and diagonalise: there is one subsequence `N_j → ∞` and
functions `c_k : [0,T*] → ℂ³` with

```
û_{N_j,k} → c_k   uniformly on [0,T*],   for every k ≠ 0 .            (3.c.2)
```

Each `c_k` is `L_k`-Lipschitz and inherits the exact constraints pointwise:
`c_{−k} = conj(c_k)` (reality), `k·c_k = 0` (divergence-free),
`c_k(0) = û₀ₖ` (since `û_{N,k}(0) = û₀ₖ` for `N ≥ |k|`).  (The first two
hold for every `û_{N_j,k}(t)` by Lemma 3.a.0 and are closed linear
conditions on the coefficients, hence pass to the pointwise limit
[turn-11 audit erratum].)

*Uniform `H⁴` control of the limit.*  For every `K` and `t`,
`Σ_{|k|≤K} |k|⁸|c_k(t)|² = lim_j Σ_{|k|≤K} |k|⁸|û_{N_j,k}(t)|² ≤ 4Y₀` by
(3.b.4); taking `K → ∞` (monotone limit of partial sums — the elementary
Fatou step),

```
u(t) := Σ_{k≠0} c_k(t) e^{ik·x}  ∈ H⁴_σ ,    ‖u(t)‖₄² ≤ 4Y₀ ,  t ∈ [0,T*]. (3.c.3)
```

The same partial-sum argument applied to (3.b.5) (integrate the uniformly
convergent finite partial sums, then let `K → ∞` monotonically) gives

```
∫₀^{T*} ‖u(t)‖₅² dt ≤ 9Y₀/(2ν) .                                      (3.c.4)
```

*Upgrade to `C([0,T*]; H³)` convergence.*  The uniform `H⁴` bound controls
all `H³` tails at once: `Σ_{|k|>K}|k|⁶|v̂_k|² ≤ K^{−2}Σ|k|⁸|v̂_k|²`, so for
every `K`, using `|a−b|² ≤ 2|a|²+2|b|²` on the tail,

```
‖u_{N_j}(t) − u(t)‖₃² ≤ K⁶ Σ_{|k|≤K} |û_{N_j,k}(t) − c_k(t)|²
                        + (2·4Y₀ + 2·4Y₀)/K²  .                        (3.c.5)
```

Given `ε > 0` choose `K` with `16Y₀/K² < ε/2`, then `j` large enough that the
finite sum (uniform in `t` by (3.c.2)) is `< ε/2`.  Hence
`u_{N_j} → u` in `C([0,T*]; Ḣ³)`, and `u ∈ C([0,T*]; H³)` (uniform limit of
continuous functions; `‖·‖_{H³} ≤ √8‖·‖₃` on mean-zero fields).  Together
with (3.c.3) this gives weak-`H⁴` control along continuous-in-time
coefficients: `u ∈ C_w([0,T*]; H⁴)` (coefficientwise continuity plus a
`t`-uniform `H⁴` bound), which is what (d) upgrades to strong `H⁴`
continuity.

### 3(d) Limit passage: the limit is a strong solution; pressure by Leray

*The coefficient integral equations.*  Fix `k ≠ 0` and `N_j ≥ |k|`.
Integrating (3.a.2),

```
û_{N_j,k}(t) = û₀ₖ − ∫₀ᵗ [ ν|k|² û_{N_j,k}(s) + i P_k B_k(u_{N_j},u_{N_j})(s) ] ds .
```

The bilinear form obeys the Cauchy–Schwarz bound `|B_k(v,w)| ≤ ‖v‖₀‖w‖₁`
(sum over the free index `m`), hence the difference bound

```
|B_k(u_{N_j},u_{N_j}) − B_k(u,u)|
  ≤ ‖u_{N_j}−u‖₀‖u_{N_j}‖₁ + ‖u‖₀‖u_{N_j}−u‖₁
  ≤ 4√Y₀ · ‖u_{N_j}−u‖₃  → 0    uniformly on [0,T*] ,
```

using `‖·‖₀, ‖·‖₁ ≤ ‖·‖₃` and (3.c.5).  (For the limit field the series
`B_k(u,u) = Σ_{m+n=k}(û_m·n)û_n` converges absolutely, by the same
Cauchy–Schwarz.)  Passing to the limit in the integral equation (uniform
convergence of the integrand):

```
c_k(t) = û₀ₖ − ∫₀ᵗ [ ν|k|² c_k(s) + i P_k B_k(u,u)(s) ] ds .           (3.d.1)
```

The integrand is continuous in `s` (bilinear bound plus `C H³` continuity of
`u`), so each `c_k ∈ C¹([0,T*])` with

```
c_k'(t) = −ν|k|² c_k(t) − i P_k B_k(u,u)(t) ,    k·c_k(t) = 0 .        (3.d.2)
```

Equation (3.d.2) is, coefficientwise, the projected Navier–Stokes system;
the paragraph "`u ∈ C¹([0,T*]; H²)`" below constructs `∂_t u` as an
`H²`-valued object and upgrades (3.d.2) to the identity
`∂_t u = νΔu − P(u·∇u)` in `H²` — no property of that object is used before
it is constructed [turn-11 audit erratum: forward pointer replacing an
apparent forward reference; no mathematics changed].  Equivalently, testing
the time-integrated equation against an arbitrary divergence-free
trigonometric polynomial recovers the weak form, and the modes span.  The product `u·∇u` is classical:
`H⁴ ⊂ C²` by the `A` embedding (`Σ|k|²|û_k| ≤ A‖u‖₄`), so `u(t)·∇u(t)` is a
continuous function whose Fourier coefficients are `iB_k(u,u)`.

*A bilinear lattice estimate (used twice below).*  For mean-zero `v, w` with
finite right-hand sides,

```
‖v·∇w‖₃ ≤ 4A(‖v‖₃‖w‖₃ + ‖v‖₂‖w‖₄) ≤ 8A ‖v‖₄ ‖w‖₄ .                   (3.d.3)
```

Proof: `|k|³ ≤ (|m|+|n|)³ ≤ 4(|m|³+|n|³)` for `m+n = k` (convexity of
`x³`), so `|k|³|(v·∇w)^̂_k|` is dominated by four times the convolution of
`(|m|³|v̂_m|)` with `(|n||ŵ_n|)` plus the convolution of `(|v̂_m|)` with
`(|n|⁴|ŵ_n|)`; Young `ℓ¹ ∗ ℓ² ⊂ ℓ²` (triangle inequality in `ℓ²`) and the
`ℓ¹` embeddings `Σ|n||ŵ_n| ≤ A‖w‖₃`, `Σ|v̂_m| ≤ A‖v‖₂` give (3.d.3).  In
particular `‖∂_t u‖₃ ≤ ν‖u‖₅ + 8A‖u‖₄²` coefficientwise, so by (3.c.4)

```
∂_t u ∈ L²(0,T*; Ḣ³)  with  ∫₀^{T*}Σ_k|k|⁶|c_k'|² dt
   ≤ 2ν²·(9Y₀/2ν) + 2T*·(32AY₀)²  < ∞ .                               (3.d.4)
```

*Norm continuity: `u ∈ C([0,T*]; H⁴)`.*  This is the one step where the
classical proofs invoke a Lions–Magenes-type lemma; on the torus it is
self-contained.  Let `F(t) := ‖u(t)‖₄² = Σ_k F_k(t)`, `F_k := |k|⁸|c_k|²`.
Each `F_k ∈ C¹` with `F_k' = 2|k|⁸Re(conj(c_k)·c_k')`, and

```
|F_k'| ≤ 2·(|k|⁵|c_k|)(|k|³|c_k'|) =: g_k ,
Σ_k ∫₀^{T*} g_k dt ≤ 2 (∫Σ|k|^{10}|c_k|²)^{1/2} (∫Σ|k|⁶|c_k'|²)^{1/2} < ∞
```

by Cauchy–Schwarz in the product of counting and Lebesgue measure, (3.c.4)
and (3.d.4).  Therefore `Σ_k |F_k(t) − F_k(0)| ≤ Σ_k∫|F_k'| < ∞` and the
sum–integral interchange (Fubini–Tonelli with the summable majorant `Σg_k`)
gives

```
F(t) = F(0) + ∫₀ᵗ Σ_k F_k'(s) ds ,      Σ_k F_k' ∈ L¹(0,T*) ,
```

so `F` is absolutely continuous — in particular continuous.  Now take
`t_n → t`: the coefficients converge (`c_k` continuous) and the `H⁴` bounds
are uniform, so `u(t_n) ⇀ u(t)` weakly in `Ḣ⁴`; together with
`‖u(t_n)‖₄ → ‖u(t)‖₄` (continuity of `F`) this forces strong convergence in
the Hilbert space:
`‖u(t_n)−u(t)‖₄² = F(t_n) − 2Re⟨u(t_n),u(t)⟩₄ + F(t) → 0`.  Hence

```
u ∈ C([0,T*]; H⁴_σ) ,    u(0) = u₀   attained strongly in H⁴ .
```

*`u ∈ C¹([0,T*]; H²)`.*  The range `{u(t)}` is compact in `Ḣ⁴` (continuous
image of `[0,T*]`), so its `H⁴` tails are uniformly small: cover by finitely
many `ε`-balls and truncate the centres.  Then, from (3.d.2),
`Σ_{|k|>K}|k|⁴|c_k'|² ≤ 2ν²Σ_{|k|>K}|k|⁸|c_k|² + 2K^{−2}(8A·4Y₀)²`
is uniformly small in `t` for large `K` (the second term by
`Σ|k|⁶|B̂_k|² ≤ ‖u·∇u‖₃²` and (3.d.3)), while every finite partial sum of
`Σ|k|⁴|c_k'(t)|²` is continuous in `t`.  Hence `t ↦ u'(t)` is continuous into
`Ḣ²`, and the coefficientwise FTC upgrades to
`u(t) = u(0) + ∫₀ᵗ u'(s)ds` in `H²`: `u ∈ C¹([0,T*]; H²)` and the equation
holds in `H²` pointwise in `t` — `u` is a strong solution in the §1 sense.
(This paragraph also delivers clause (i) of EXT-P2 *for this constructed
solution*; it does not discharge EXT-P2, per §1.2.)

*Pressure via Leray.*  Define `p̂₀ = 0` and
`p̂_k := i k·(u·∇u)^̂_k / |k|²`.  Then `∇p = −(I−P)(u·∇u)` mode by mode,
`|p̂_k| ≤ |(u·∇u)^̂_k| / |k|`, so `‖p‖_{Ḣ⁴} ≤ ‖u·∇u‖₃ ≤ 8A‖u‖₄² < ∞`, and
`t ↦ p(t)` is continuous into `H⁴` by the bilinear continuity
`‖u·∇u − v·∇v‖₃ ≤ 8A(‖u‖₄+‖v‖₄)‖u−v‖₄` from (3.d.3).  Adding
`(I−P)(u·∇u) + ∇p = 0` to the projected equation yields the unprojected
Navier–Stokes system `∂_t u + u·∇u + ∇p = νΔu` classically in `x`,
completing the existence half of Theorem 1.

### 3(e) Uniqueness: L² Grönwall for the difference

Let `u, v` be strong solutions on `[t₀,t₁]` with `u(t₀) = v(t₀)` and set
`w := u − v ∈ C([t₀,t₁]; H⁴_σ) ∩ C¹([t₀,t₁]; H²)`.  Subtracting the two
equations and splitting the bilinear difference exactly,

```
∂_t w = νΔw − P(u·∇w) − P(w·∇v) .
```

The function `h(t) := ‖w(t)‖₀²` is `C¹` with `h' = 2⟨∂_t w, w⟩₀` (calculus of
a `C¹` curve in `L²`).  Term by term, with all fields in `H⁴ ⊂ C²` so that
every integration by parts is classical:

* `⟨νΔw, w⟩₀ = −ν‖w‖₁² ≤ 0` (Parseval);
* `⟨P f, w⟩₀ = ⟨f, w⟩₀` (`P` self-adjoint, `Pw = w`);
* `⟨u·∇w, w⟩₀ = ½∫ u·∇|w|² = −½∫ (div u)|w|² = 0` (periodic boundary,
  `div u = 0`);
* `|⟨w·∇v, w⟩₀| ≤ ‖|∇v|_F‖_∞ ‖w‖₀² ≤ A‖v‖₃ ‖w‖₀² ≤ A‖v‖₄ ‖w‖₀²`, since
  the mode-`k` contribution to `∇v` has Frobenius norm `|k||v̂_k|` and
  `Σ|k||v̂_k| ≤ A‖v‖₃` by the §5.1 Cauchy–Schwarz.

Hence, with `M := sup_{[t₀,t₁]} ‖v(t)‖₄` (finite by continuity),

```
h'(t) ≤ 2A‖v(t)‖₄ · h(t) ≤ 2AM · h(t) ,
```

and `(h e^{−2AM(t−t₀)})' ≤ 0` gives the explicit Grönwall inequality

```
‖u(t) − v(t)‖₀² ≤ ‖u(t₀) − v(t₀)‖₀² · exp( 2A (t−t₀) · sup_{[t₀,t]}‖v‖₄ ) ,
A ≤ 4.114 .                                                            (3.e.1)
```

With `u(t₀) = v(t₀)` the right side vanishes, so `w(t) = 0` in `L²` for all
`t`; all coefficients agree, hence `u = v` in `H⁴`.  ∎  (Inequality (3.e.1)
is also the quantitative stability statement behind the chain's boundary
step; the constant is stated in terms of the `H⁴` norm through
`‖v‖₃ ≤ ‖v‖₄ ≤ ‖v‖_{H⁴}`.)

### 3.6 Audit box

```
STATUS: UNAUDITED.  This §3 proof is classical mathematics written down by
this repository; it has not been independently audited, and none of its
infinite-dimensional steps is formalised.  EXT-P1 therefore keeps
proved: false and axiomatised_in_lean: false in EVERY payload, and the
checkers keep enforcing that, until BOTH of the following happen:
  (1) an independent audit of this section, recorded with auditor and date;
  (2) formalisation of its finite-dimensional parts in Lean.
Even then the flag semantics do not change by fiat: proved means proved
inside the repository's verified layer, and the infinite-dimensional steps
below remain outside it.

Per-step Lean formalisation obligations:
  finite-dimensional (near-term):
    L-EXT-P1-a1  Picard local existence/uniqueness for the Galerkin
                 coefficient ODE — GalerkinPicard.lean, lemmas
                 quadratic_ode_local_solution / quadratic_ode_unique (lane L,
                 in tree this turn, CI build pending; mathlib
                 IsPicardLindelof + ODE_solution_unique_of_mem_Icc).
                 Remaining glue: apply the abstract quadratic field to the
                 concrete band-limited Navier-Stokes coefficient vector.
    L-EXT-P1-a2  L2 energy neutrality + global extension at fixed N —
                 abstract core EXISTS (GalerkinNoBlowup.lean); glue to the
                 ODE layer open.
    L-EXT-P1-b1  the per-N H4 self-estimate (3.b.2) as a finite rational
                 inequality schema — open; bricks in TrackPFourier.lean.
    L-EXT-P1-b2  the scalar bound Y <= 4*Y0 on [0,T*] (the phi = Y^{-1/2}
                 argument) — open; scalar real analysis only.
    L-EXT-P1-b3  the lattice constant A (one infinite series with an
                 elementary tail) — open; mathlib PSeries/tsum suffice.
  infinite-dimensional (OPEN, no mathlib substrate — the function spaces
  themselves do not exist there):
    L-EXT-P1-c1  Arzela–Ascoli + diagonal extraction over all modes.
    L-EXT-P1-c2  uniform H4 tail bound, C([0,T*];H3) upgrade, Fatou to
                 L-infinity H4 and L2 H5-dot.
    L-EXT-P1-d1  lattice bilinear estimates (l1*l2 Young; (3.d.3)).
    L-EXT-P1-d2  limit passage into the coefficient ODEs (3.d.1)-(3.d.2).
    L-EXT-P1-d3  absolute continuity of the H4-dot norm and strong H4
                 continuity; C1-H2 regularity; pressure summation.
    L-EXT-P1-e1  trilinear cancellation in L2 and the Gronwall closure
                 (mathlib norm_le_gronwallBound_of_norm_deriv_right_le is
                 the right scalar shape).
Wording contract unchanged: the only true sentences remain "conditional PDE
certificate assuming EXT-P1/P2/P3" and "the finite-dimensional Galerkin
enclosure and the scalar control ODE are verified unconditionally".
```

## 4. What an audit would — and would not — change

Suppose §3 were independently audited tomorrow and the audit passed.

**What stays, all of it.**
* `verify_torus_slab_certificate` keeps requiring, for each of EXT-P1/P2/P3:
  `proved: false`, `axiomatised_in_lean: false`, and the statement text
  byte-identical to `EXTERNAL_THEOREMS`.  An audited paper proof is still a
  paper proof; the `proved` flag is reserved for the repository's verified
  layer (exact rational checkers and Lean), and §3 lives in neither.
* Every conditional wording stays: certificate conclusions keep
  `conditional_on: ["EXT-P1", "EXT-P2", "EXT-P3"]`, the chain checker keeps
  its `ALLOWED_WORDING`/`FORBIDDEN_WORDING` lists, and the only permitted
  summary sentences remain the two of the wording contract.  This holds
  *until formalisation*, not until audit: an audit upgrades confidence in
  the mathematics, not the status of the machine-checked chain.
* The division of labour stays: unconditional means the finite-dimensional
  Galerkin enclosure and the scalar control ODE, nothing more.

**What an audit would actually buy.**  Provenance: the ledger (this note,
STATUS.md at the coordinator's discretion) could record "EXT-P1 paper proof
audited by X on date D", and the §3.6 finite-dimensional obligations become
well-posed Lean targets with an audited reference text — the difference
between formalising folklore and formalising a fixed document.  Nothing in
any payload changes until those Lean targets are discharged, and even full
discharge of the finite-dimensional list leaves EXT-P1 external, because the
infinite-dimensional steps (L-EXT-P1-c1 … e1) are the theorem.

**Why the axiom-free rule forbids the shortcut.**  The tempting move — `axiom
extP1 : ...` in Lean, then "prove" the conditional theorems unconditionally —
is forbidden by `LEAN4_VERIFICATION_POLICY.md` and by the standing rule
recorded in `torus_aposteriori.py` ("recorded faithfully in the payload and
never inserted into Lean as axioms"), for three reasons that are checked, not
aspirational.  First, an axiom does not appear in the *statement* of the
theorems that use it: downstream results would silently carry an unverifiable
hypothesis, converting "conditional PDE certificate assuming EXT-P1/P2/P3"
into exactly the forbidden unconditional claim, invisibly.  Second, the
trusted base would grow by a sentence of informal mathematics whose
faithfulness nothing checks — the same failure mode as a [recalled] constant,
but embedded in the kernel's trust boundary where `#print axioms` is the only
witness.  Third, the repository's honesty architecture depends on the
checker being *able* to reject a payload that claims too much
(`proved is False` is a hard test in `verify_torus_slab_certificate`); an
axiomatised EXT-P1 would make the Lean artifact assert more than the payload
is permitted to say, and the two layers would disagree about what has been
established.  A citation is a citation.  The path for EXT-P1 into Lean is
the §3.6 obligation list, in order, or nothing.
