# Morosi–Pizzocchero and the whole-space obstruction

**Status: a verdict document with its evidence attached.**  It records, as
faithfully as we can state it, the one existing framework that turns a computed
approximate solution into a certified lower bound on the existence time of the
exact solution; it maps that framework onto this repository's objects; and it
says why it must not be ported here.  The verdict is negative and the reason is
structural, so the document is written to be re-checkable rather than believed.

Every factual claim carries one of three labels, and they are used strictly.

* **[verified]** — read from the source in this session (2026-07-29), by
  fetching the arXiv abstract page or the API listing.  The exact route is named.
* **[corroborated]** — a secondary source consulted in this session states it,
  but the primary source was not read.
* **[recalled]** — stated from memory, not re-checked against the source in this
  session.  Treat these as leads, not as citations.  Anything load-bearing that
  is only **[recalled]** is flagged where it is used.

Nothing in this document is a claim about the Navier–Stokes problem itself.

## 1. The framework, faithfully

### 1.1 Setting

Morosi and Pizzocchero work on the torus `T^d`, in the Sobolev spaces
`H^n_{Σ0}(T^d)` of divergence-free, **zero-mean** vector fields, for real
`n > d/2 + 1`.  In applications `d = 3` and `n = 3`.  **[verified]**: the
abstract of arXiv:1104.3832 states "on a torus `T^d` in the functional setting
of the Sobolev spaces `H^n(T^d)` of divergence free, zero mean vector fields on
`T^d`, for `n > d/2+1`".

The approximate solution is assumed regular enough to be differentiated twice
above the working level: `u_a ∈ C(H^5) ∩ C¹(H^3)` for the `d = 3, n = 3`
applications.  **[recalled]** — the two-above-working-level shape is certain, the
exact function-space pair as written is not re-checked.

### 1.2 The objects

Write `ℙ` for the Leray projection and `P(v, w) = ℙ((v·∇)w)` for the projected
bilinear form.  Given an exact solution `u` of

```
∂_t u = νΔu + P(u, u) + f ,      u(0) = u₀ ,
```

and any sufficiently regular field `u_a` — the **approximate solution**, which
need not solve anything — define:

* the **differential residual**
  ```
  e(u_a) := ∂_t u_a - νΔu_a - P(u_a, u_a) - f ,
  ```
  everything by which `u_a` fails to solve the equation, with
  `‖e(u_a)(t)‖_m ≤ ε_m(t)`;
* the **datum error** `δ_m ≥ ‖u_a(0) - u₀‖_m`;
* the **growth estimator** `D_m(t) ≥ ‖u_a(t)‖_m`, a known function bounding the
  approximate solution's own norm — note that "Sobolev norm estimator", the
  phrase used elsewhere in this repository, is **not** their terminology;
* the **existence interval** `T_c`, the certified lower bound on the time of
  existence of the exact solution.

Two constants enter, and they are the two inequalities of §1.2 of
[`hs5_function_space.md`](hs5_function_space.md):

```
‖P(v, w)‖_n ≤ K_n ‖v‖_n ‖w‖_{n+1}          the "basic" inequality
|⟨P(v, w), w⟩_n| ≤ G_n ‖v‖_n ‖w‖_n²        the Kato inequality
```

`K_n` loses a derivative and `G_n` does not; it is `G_n` that closes the argument
in one unknown.  **[verified]** that these two inequalities are the subject of
two dedicated papers, by title: arXiv:1007.4412 "On the constants in a basic
inequality for the Euler and Navier-Stokes equations" and arXiv:1009.2051 "On the
constants in a Kato inequality for the Euler and Navier-Stokes equations", both
Morosi–Pizzocchero; and arXiv:1511.00533 "New results on the constants in some
inequalities for the Navier-Stokes quadratic nonlinearity".  **[verified]** that
the abstract of 1009.2051 says "numerical values of our upper and lower bounds
are reported for `d=3` and some values of `n`" without displaying them.

### 1.3 The control ODE

The whole framework reduces to one scalar differential inequality.  If `R_n` is
locally Lipschitz, satisfies

```
d⁺R_n/dt  ≥  -ν R_n + ( G_n D_n(t) + K_n D_{n+1}(t) ) R_n + G_n R_n² + ε_n(t) ,
R_n(0) ≥ δ_n ,                                                          (1.3)
```

on `[0, T_c)`, with `d⁺` the upper right Dini derivative, then the exact solution
exists on `[0, T_c)` and

```
‖u(t) - u_a(t)‖_n ≤ R_n(t)     for all t ∈ [0, T_c) .
```

**[verified]** from the abstract of 1104.3832: "this allows to infer a lower
bound `T_c` on the time of existence of the exact solution `u` analyzing a
posteriori any approximate solution `u_a`, and also to construct a function `R_n`
such that `‖u(t) - u_a(t)‖_n ≤ R_n(t)` for all `t ∈ [0,T_c)`.  Both `T_c` and
`R_n` are determined solving suitable 'control inequalities', depending on the
error of `u_a`."  The exact right-hand side of (1.3) as written is
**[recalled]**.

The mechanism has three parts and it is worth separating them, because only two
of them survive the move to `ℝ³`.

1. **Derivation.**  Set `W_n(t) = ‖u(t) - u_a(t)‖_n`.  Differentiating and using
   the two inequalities of §1.2 term by term produces exactly the right-hand
   side of (1.3) as an **upper** bound for `d⁺W_n/dt`.  The `-ν R_n` comes from
   the viscous term, the `G_nD_n + K_nD_{n+1}` from the cross terms between the
   error and the approximate solution, the `G_nR_n²` from the error against
   itself, and `ε_n` from the residual.
2. **Comparison.**  Since `W_n` satisfies the inequality with `≤` and `R_n` with
   `≥`, and `R_n(0) ≥ W_n(0)`, a Chaplygin/Dini comparison lemma gives
   `W_n ≤ R_n` for as long as `R_n` exists.  This is where the Dini derivative
   is needed: `W_n` is a norm and is not differentiable everywhere.
3. **Continuation.**  As long as `R_n` stays finite, `‖u‖_n` stays finite, and
   the local existence theory in `H^n` with `n > d/2 + 1` continues the solution.
   Blow-up of `R_n` at time `T_c` is the only thing that stops the argument, and
   it does not imply blow-up of `u`.

### 1.4 Citations

| paper | content | label |
|---|---|---|
| arXiv:1104.3832, Morosi & Pizzocchero, "On approximate solutions of the incompressible Euler and Navier-Stokes equations" | the general theory; Definitions 3.4, 4.4, 4.6 | title, authors, abstract, torus setting, control inequalities, `T_c`, `R_n` **[verified]**; the definition numbers **[recalled]** |
| arXiv:1310.5642, Morosi, Pernici & Pizzocchero, "A posteriori estimates for Euler and Navier-Stokes equations" | the applications; Def 2.2 | title, authors **[verified]**; that it treats the Behr–Nečas–Wu datum with a Taylor polynomial of order 52 and has a section "the Euler equations on `T³`, with the Behr–Nečas–Wu initial datum" **[verified]** from the abstract and section list; the definition number **[recalled]** |
| arXiv:1405.3421, Morosi & Pizzocchero, "Smooth solutions of the Euler and Navier-Stokes equations from the a posteriori analysis…" | the `C^∞` variant; Prop 4.4 | title **[verified]**; abstract **[verified]** as presenting "a variant, based on a `C^∞` formulation of the Cauchy problem" giving "a bound on the Sobolev distance of any order between the exact and the approximate solution"; the proposition number **[recalled]** |
| arXiv:1007.4412 / arXiv:1009.2051 / arXiv:1511.00533 | the constants `K_n` and `G_n` | titles and authors **[verified]**; all numerical values **[recalled]** |

### 1.5 Reference constants for `d = 3`

**All numbers in this subsection are [recalled].**  None of them was read from a
paper in this session — the arXiv abstract pages do not display them and the PDF
text was not extractable in this environment.  They are recorded because their
*magnitudes* are what matters for the verdict, and are flagged so that nobody
quotes them as sourced.

```
G_3  = 0.438        upper bound on the Kato constant
G_3⁻ = 0.114        lower bound on the sharp Kato constant
K_3  = 0.323        upper bound on the basic-inequality constant
K_33⁺ = 0.320       a later, slightly improved upper bound
K_3⁻ = 0.179        lower bound on the sharp constant
```

The gap between upper and lower bounds is a factor of roughly 3 to 4 — about 3.8
for `G_3` and about 1.8 for `K_3` — and that uncertainty propagates directly into
`T_c`, since (1.3) is a Riccati inequality whose blow-up time depends on the
coefficient of the quadratic term.  Any port of this framework inherits the
uncertainty; no port can be sharper than the constants it uses.

### 1.6 Benchmark outcomes

**All numbers here are [recalled] except where noted.**

For the **Navier–Stokes** equations on `T³` with the Behr–Nečas–Wu datum and a
Galerkin approximate solution using 150 Fourier modes:

```
ν = 0    T_c = 0.0667              (the Euler case)
ν = 7    T_c = 0.2386
ν = 8    T_c = +∞ , i.e. global existence
```

against a naive criterion — the elementary smallness condition that does not use
the approximate solution at all — needing `ν ≥ 67.58` for the same conclusion.
The a posteriori analysis therefore improves the global-existence threshold by a
factor of about **8.4**.  **[corroborated]**: a web search on 2026-07-29 returned
the statement that Morosi and Pizzocchero, analysing the Behr–Nečas–Wu datum by
Galerkin approximation, obtained "a lower bound `T > 0.066` for the time of
existence in `H^3` for the Euler equation", and that "for the Navier–Stokes
equations they showed global existence when the viscosity coefficient is above an
explicit threshold".  That corroborates `0.0667` and the shape of the `ν = 8`
result; it does not corroborate `0.2386`, `8`, `67.58` or the mode count.

For the **Euler** equations on `T³` with the same datum and a degree-52 Taylor
approximant in time:

```
T_c = 0.242        certified
T ≈ 0.32           conjectured blow-up time
```

**[verified]** that the degree-52 Taylor polynomial for the Behr–Nečas–Wu datum
appears in 1310.5642; the two times are **[recalled]**.

The last line is the one to internalise.  Even at degree 52, with a datum chosen
for its numerical tractability, on a compact domain, with no truncation error and
no discrete-to-continuum gap, **the certified interval does not reach the
conjectured singularity**.  It stops at about three quarters of the way.  That is
the state of the art for this style of argument, and it is a strong prior on what
any descendant of it can achieve.

## 2. Mapping onto this repository

| Morosi–Pizzocchero object | counterpart here | where | notes |
|---|---|---|---|
| approximate solution `u_a` | the integrated full-nonlinear trajectory of `PicardLadder` | `picard_continuation.py`, class `PicardLadder` | the ladder also carries Picard levels 0, 1, 2, but those are *iterates*, a different kind of approximation: they are compared against the full solution, not residual-controlled against the PDE |
| differential residual `e(u_a)`, with `‖e(u_a)‖_m ≤ ε_m` | **partial, and in the wrong norms** | `slab_certificate.py` (integrator local defect, trapezoid `2.35e-6`, Simpson `8.25e-10`); `ContinuationCheckpoint.poisson_residual`, `.relative_divergence`, `.curl_defect` | every measured quantity is a discrete `L^∞` or grid ratio at nodes.  Converting any of them into `‖·‖_{H^n}` is exactly `HS-5`, the largest gap.  There is **no** `ε_n` in this repository |
| datum error `δ_m` | zero for the current runs | the datum is evaluated exactly on the grid; see `hs_error_propagation.md` §4, row `R_embed` | it becomes nonzero the moment the truncated box is accounted for, at which point it is `R_trunc` and needs `HS-4`, which has no MP analogue |
| growth estimator `D_m(t)` | **no counterpart at the right order** | `ContinuationCheckpoint.max_u1`, `.max_omega1`, `.max_speed`, `.kinetic_energy`, `.enstrophy`, `.l3_norm` | these are `L^∞`, `L²` and `L³` quantities.  `D_3` and `D_4` would be `H^3` and `H^4` norms of the computed field, which `norm_transfer.py` will make computable via Theorem NT — that module does not exist yet |
| control inequality / control ODE (1.3) | **no counterpart** | — | there is no `control_ode.py` in `src/ns_certificate_lab`; verified by search on 2026-07-29 |
| existence interval `T_c` | **no counterpart; do not conflate** | `ContinuationResult.reached_tau`, `.stopped_reason` | these record where the *computation* stopped — resolution, front points, departure thresholds.  `T_c` records where *existence of the exact solution* is certified.  Reporting one as the other would be the single most misleading thing this repository could do |
| exact–approximate distance `‖u - u_a‖_n ≤ R_n` | **no counterpart** | `ContinuationCheckpoint.distance_to_first_picard`, `.distance_to_second_picard` | these are distances between two *computed* trajectories.  They diagnose whether the run has left the first Picard iterate.  They say nothing about any exact solution |
| the constants `G_n`, `K_n` | **no counterpart, and none available** | — | MP's values are lattice sums on `T^d`; see §3 |
| the blow-up criterion used for continuation | present, but a different one | `critical_l3.py`, `l3_certificate.py`, `l3_generation.py`, `docs/research_notes/critical_l3_obstruction.md` | MP continue in `H^n` and invoke Beale–Kato–Majda; this repository works with the Escauriaza–Seregin–Šverák endpoint `L^∞_t L³_x` criterion.  The two are compatible but not interchangeable |
| the Grönwall/comparison step | present, constant coefficient only | `NSSingularity.norm_le_simple_gronwall` in `formal/NSSingularity/CertificateLayer.lean` | the variable-coefficient form is `HS-6`; mathlib has `norm_le_gronwallBound_of_norm_deriv_right_le` for constant coefficients |
| the domain | `T^d`, compact, finite volume | `ℝ³`, whole space, truncated to a box | this single row is the source of every entry in the right-hand column of §3 |

The shape of the mapping is worth stating plainly: this repository has an
approximate solution and a great deal of *diagnostic* information about it, and
essentially none of the *certified* information that (1.3) consumes.  What is
missing is not the framework — it is `ε_n`, `D_n` and `δ_n` as rigorous `H^n`
quantities.  Two of those three are downstream of `HS-5`.

## 3. What transfers and what does not

| directly reusable on `ℝ³` | new lemma needed for the whole space |
|---|---|
| **The Riccati/Chaplygin comparison structure and the Dini-derivative formulation.**  Nothing in the comparison lemma knows about the domain.  A locally Lipschitz `R` satisfying the reversed inequality with `R(0) ≥ W(0)` dominates `W`, on any domain, and the Dini derivative is needed for the same reason — `W` is a norm | **The Poincaré / spectral gap.**  The `-ν R_n` term in (1.3) exists **only** because on `T^d` with zero mean the Fourier index runs over `Z^d \ {0}`, so `\|k\| ≥ 1` and `‖∇v‖_n ≥ ‖v‖_n`.  On `ℝ³` the spectrum of `-Δ` is `[0, ∞)`, the sharp constant is `0`, and the term **vanishes entirely**.  This is the severest obstruction and it is discussed below the table |
| **The term-by-term derivation of the ODE from `W_n = ‖u - u_a‖_n`.**  Differentiating the norm, splitting the bilinear difference by `u·∇u - u_a·∇u_a = e·∇u + u_a·∇e`, and applying the two inequalities of §1.2 is a calculation in a Hilbert space, not on a domain.  It is the same calculation as §2 and §3 of `hs_error_propagation.md` | **The Kato constant as a discrete lattice sum.**  MP's numerical values for `G_n` and `K_n` are obtained by bounding and then evaluating sums over `Z^d`.  On `ℝ³` the corresponding object is a continuous supremum over `ℝ³ × ℝ³`, and the entire computational machinery — truncation of the sum, interval evaluation of the finitely many remaining terms, the tail bound — has no analogue.  The constants would have to be re-derived from scratch |
| **The Beale–Kato–Majda criterion.**  Valid on `ℝ³`; it is one of the standard continuation criteria there.  Its role — turning "the `H^n` norm did not blow up" into "the solution continues" — is domain-independent | **The zero-mean scaffolding `Σ0`.**  `H^n_{Σ0}` is a Hilbert space on which the homogeneous and inhomogeneous norms are equivalent with computable constants, and on which `Λ^n` is boundedly invertible.  On `ℝ³` "zero mean" is not a constraint one can impose — it is not preserved, not meaningful for `L²` fields without extra decay, and carries no information without finite volume.  Every place MP use it, a different argument is needed |
| **The threshold `n > d/2 + 1` and its role.**  The reason for it — `n > d/2` gives only the derivative-losing product estimate, `n > d/2 + 1` buys the Kato inequality that closes the ODE in `R_n` alone — is a Sobolev-embedding argument, identical on `T^d` and `ℝ³`.  See `hs5_function_space.md` §1.2 | **The exact finite-support Galerkin residual.**  On `T^d` a Galerkin approximant is a trigonometric polynomial; its residual is computable exactly, in closed form, with finitely many terms, because the nonlinearity of a band-limited field is band-limited.  On `ℝ³` there is no such basis, the approximant is a grid function, and the residual is only known at nodes.  This is precisely `HS-5` and it is the reason MP never needed anything like it |
| **The absence of boundary terms.**  Every integration by parts in the derivation is boundary-free on `T^d` because it is compact without boundary, and boundary-free on `ℝ³` because the fields decay.  Nothing is lost here | **The Leray projection in the Sobolev inner product — a partial port, not free.**  `ℙ` is a Fourier multiplier on `ℝ³` too, so it remains an orthogonal projection in every `H^n` inner product, and that much transfers.  What does not transfer is its interaction with the low frequencies: on `T^d` the zero mode is simply removed, whereas on `ℝ³` the low-frequency behaviour of `ℙ` is exactly where the homogeneous and inhomogeneous norms part company.  And this repository does not use `ℙ` at all — the pressure is eliminated by the stream-function/swirl formulation — so porting means either introducing it or redoing every estimate in the `(u₁, ψ₁, ω₁)` variables, which is a different calculation |
| | **The domain-truncation budget.**  MP have no truncation error: the domain is compact and the computation covers it.  Here the computation covers a box inside `ℝ³`, so there is a fifth error source with no counterpart in (1.3) — `R_trunc`, bounded pointwise by the multipole tails of `free_space_recovery.py` and `green_derivative_tail_bounds.md`, and requiring `HS-4` to become an `H^n` quantity |

### The `-ν R_n` term, in detail

This is the load-bearing point of the whole document, so it is spelled out.

On `T^d` with zero mean, the viscous term contributes `-ν‖∇e‖²_n` to the
derivative of `‖e‖²_n`, and the spectral gap converts it into `-ν‖e‖²_n`, which
is the `-ν R_n` of (1.3).  Consider what that term does to the Riccati structure.
Dropping the `D` terms for clarity, the control ODE is

```
R′ = -ν R + G R² + ε .
```

The right-hand side has a **negative** root — a stable equilibrium — whenever
`ν² > 4Gε`, that is whenever the residual is small relative to the viscosity.
Starting below that equilibrium, `R` never reaches it, `R` stays bounded for all
time, and `T_c = +∞`.  **That is the entire mechanism behind every
global-existence-from-smallness corollary in the programme**, including the
`ν = 8` benchmark of §1.6.  It is also the mechanism behind the decay
corollaries, which are just the statement that `R` relaxes toward the equilibrium.

On `ℝ³` the term is gone.  The ODE is

```
R′ = (G D_n + K D_{n+1}) R + G R² + ε ,
```

with every coefficient nonnegative.  `R` is strictly increasing whenever `ε > 0`,
and a Riccati equation with a positive quadratic coefficient and no negative
linear term blows up in finite time from **any** positive initial value and any
positive `ε`.  So:

* `T_c` is always finite;
* no smallness of `δ_n` or `ε_n`, relative to `ν` or to anything else, produces
  global existence;
* there is no decay corollary, because there is no equilibrium to decay to;
* the only remaining content is a finite-time bound whose length is governed by
  `1/(G R(0))`-type quantities, that is, by how small the residual and datum
  error are — which is exactly the regime in which the Euler benchmark of §1.6
  reached only three quarters of the conjectured blow-up time on a *compact*
  domain with an *exact* residual.

Replacing the missing gap is not a technical exercise.  There is no inequality on
`ℝ³` that plays the role of `‖∇v‖ ≥ ‖v‖`; the failure is a fact about the
spectrum of the Laplacian on an infinite-volume space, not a gap in the
literature.  Any whole-space version of this programme must be built on a
different mechanism — decay of the heat semigroup in a scaling-critical norm, for
instance — and that is a different framework, not a port.

## 4. The verdict

> **Do not port Morosi–Pizzocchero to this repository.**

The `-ν R_n` term is load-bearing — it is what every quotable corollary of the
programme rests on — and it does not exist on `ℝ³`.  What would remain after the
port is the comparison lemma and the term-by-term derivation, both of which this
repository already has in a more directly usable form in
`hs_error_propagation.md` §2–§3 and in `NSSingularity.norm_le_simple_gronwall`.
The cost of the port would be re-deriving `G_n` and `K_n` on `ℝ³` from scratch,
which is a research project in its own right, in exchange for a finite-time bound
that the existing derivation already produces.  The trade is bad.

**If a residual-to-existence framework is wanted on `ℝ³`, the better template is
Brunk–Giesselmann–Tscherpel, arXiv:2509.25105, "A posteriori existence of strong
solutions to the Navier-Stokes equations in 3D".**  **[verified]**: title,
authors and abstract read on 2026-07-29.  Three reasons it is a better template:

1. **It uses only negative Sobolev norms of the residual.**  The abstract states
   "The a posteriori criterion that, if satisfied, verifies existence of strong
   solutions, involves only negative Sobolev norms of the residual".  The
   specific norms `L²(W^{-1,2})` and `L³(W^{-1,3})` are **[recalled]**.  This
   matters enormously here: a negative-order norm of a residual is *easier* to
   bound from nodal data than a positive-order one, because it is dual to a
   space of test functions and does not require differentiating the discrete
   residual.  It would make `HS-5` a materially easier lemma than the `H^n`
   version.
2. **Its blow-up criterion is `L^∞(L³)`.**  The abstract names "a version of the
   celebrated blow-up criterion in the critical space `L^∞(L³)` by Iskauriaza,
   Serëgin and Shverak (2003)", built on "a conditional stability estimate in
   `L²` and `L³`".  That is **the same critical norm** this repository already
   works in — see `critical_l3_obstruction.md`, `critical_l3.py`,
   `l3_generation.py`, and the `L³` energy identity
   `F′ = 3∫p div(|u|u) - 3ν∫(|u||∇u|² + |u||∇|u||²)`.  The alignment is not
   coincidental; it is the natural criterion for whole-space work because it is
   scaling-critical.
3. **The criterion is fully computable without extra assumptions.**  The abstract
   states "A posteriori error estimates allow us to derive a fully computable
   criterion without imposing any extra assumptions on the solution".  That is
   precisely the property `HS-5` currently lacks.

Two limitations, recorded so that the recommendation is not oversold.  It is
**periodic** — **[verified]**, the abstract says "subject to periodic boundary
conditions" — so the domain problem is not solved by adopting it, only relocated
to a place where the fix is a truncation argument rather than a missing spectral
gap.  And its own constants are described as crude — **[recalled]**, from the
paper body, not from the abstract, which says only "While limited to short time
intervals".  It is a template for the *structure* of the argument, not a source
of numbers.

**No whole-space version of the Morosi–Pizzocchero programme exists, by them or
by anyone we can find.**  **[verified]** on 2026-07-29 by querying the arXiv API
for `au:"Pizzocchero" AND cat:math.AP`, which returned 17 entries spanning 2003
to 2026.  A correction to the figure this note was drafted with: the count is 17
on that query, not 14, and it is not literally true that every entry is on `T^d`.
The accurate statement is finer and, for our purposes, stronger:

* every entry belonging to the **a posteriori / approximate-solution programme**
  for fluids is on a torus — 1104.3832, 1310.5642, 1405.3421, 1203.6865
  (Behr–Nečas–Wu power series, `T³`), 0909.3707 (`H¹` setting, `T³`), 1304.2972
  and 1402.0487 (Reynolds expansions), 1511.00533, 1905.13722 (MHD), 2602.06933
  (MHD global stability), 1007.4412 and 1009.2051 (the constants);
* the only whole-space entries are the abstract inequality-constant papers —
  math/0501229 and 0902.0708 on multiplication in Sobolev spaces, and 1611.00734
  on fractional Gagliardo–Nirenberg and Sobolev constants — which are **inputs**
  to the constants, not applications of the framework;
* math-ph/0309016 and 0709.1670 are the abstract semilinear-evolution
  antecedents, in Banach spaces, and are the closest thing to a
  domain-independent statement of the method.

The domain labels in that list come from a machine summary of the API response
rather than from reading each abstract, so treat the per-entry classification as
**[corroborated]** and the overall pattern as solid.  The relevant conclusion —
that twenty-three years of work by the originators produced no whole-space
version — stands.

## 5. The other frameworks, and why each is or is not the template

All five are periodic.  That is not a coincidence and it is the point of this
section: every existing computer-assisted framework for three-dimensional
Navier–Stokes lives on a torus, because every one of them needs either a Fourier
lattice, a spectral gap, or a compact domain, and usually all three.

**Chernyshenko–Constantin–Robinson–Titi**, "A posteriori regularity of the
three-dimensional Navier–Stokes equations from numerical computations", J. Math.
Phys. **48** (2007) 065204.  **[corroborated]** by a web search on 2026-07-29 for
the title, authors, journal and volume.  The structure: a sufficiently smooth
strong solution is robust under small perturbations of the datum and forcing;
Galerkin approximations converge to it if it exists; therefore existence of a
sufficiently regular strong solution can be verified by sufficiently refined
numerical computation.  **Not the right template**, because the verification is
conditional on a rigorous bound for the numerical error, which is assumed
available rather than produced.  That is `HS-5` again, moved one step earlier in
the argument.

**Dashti–Robinson**, "An a posteriori condition on the numerical approximations
of the Navier–Stokes equations for the existence of a strong solution", SIAM J.
Numer. Anal. **46** (2008) 3136–3150, arXiv:math/0701341.  **[verified]** title
and arXiv identifier from the arXiv listing surfaced in the same search;
**[recalled]** that the setting is periodic.  Closest in spirit to
Brunk–Giesselmann–Tscherpel: a computable condition on the numerical solution
implying existence of a strong solution.  Worth reading before writing anything
in this direction, as the earlier and simpler statement of the same idea.

**Robinson–Sadowski**, "Numerical verification of regularity in the
three-dimensional Navier–Stokes equations", arXiv:math/0701268; and the
bounded-sets-of-initial-data version in Asymptotic Analysis (2008).
**[verified]** title, authors and the abstract's content — a numerical procedure
verifying regularity for **bounded sets** of initial conditions, proved to
terminate in finite time under the assumption that the solutions really are
regular.  **[recalled]** that the domain is periodic; the abstract does not say.
Interesting for a different reason from the others: it certifies a whole *set* of
data rather than one trajectory, which is closer to what a search programme
wants.  Not a template for the residual problem.

**van den Berg–Breden–Lessard–van Veen**, "Spontaneous periodic orbits in the
Navier–Stokes flow", arXiv:1902.00384.  **[corroborated]** on 2026-07-29:
constructive existence proofs for periodic orbits of the forced autonomous
Navier–Stokes equations **on the three-torus**, by posing a zero-finding problem
on a Banach space of **geometrically decaying Fourier coefficients** and applying
a Newton–Kantorovich theorem, with applications to Taylor–Green forcing.  **Not
the right template**, and the reason is visible in that description: the method is
a fixed-point argument in a weighted `ℓ¹` space of Fourier coefficients.  It needs
a discrete lattice and geometric decay, neither of which exists on `ℝ³`.  The
radii-polynomial technique itself is excellent and is the right tool for a
different question — validated continuation of a *finite-dimensional* reduction —
but not for a whole-space PDE residual.

**Cyranka–Zgliczyński**, rigorous integration of dissipative PDEs by
self-consistent bounds on Galerkin projections: "Existence of globally attracting
solutions for one-dimensional viscous Burgers equation with nonautonomous forcing
— a computer assisted proof", SIAM J. Appl. Dyn. Syst. **14**-2 (2015);
"Stabilizing effect of large average initial velocity in forced dissipative PDEs
invariant with respect to Galilean transformations", J. Differential Equations
**261**-8 (2016) 4648–4708; and with Mucha and Titi, "Stabilizing the long-time
behavior of the forced Navier–Stokes and damped Euler systems by large mean
flow", Physica D **369** (2018) 18–29.  **[corroborated]** titles and venues on
2026-07-29.  The method rigorously encloses a Galerkin projection and controls the
tail by self-consistent bounds.  **Not the right template** for the same reason as
the previous entry, and more sharply: "the tail" is defined as the part of the
Fourier series above the projection, so the whole apparatus is a statement about a
lattice.  On `ℝ³` the analogous object would be a high-frequency remainder in a
continuous transform, for which the self-consistent bound machinery does not
exist.

The common structure of all five failures is worth naming once.  Each framework
converts an infinite-dimensional problem into a finite computation by a *discrete
spectral decomposition with a controlled tail*.  On a torus that decomposition is
canonical, orthogonal, and exact on band-limited fields.  On `ℝ³` it does not
exist.  Every whole-space attempt must instead convert by *spatial truncation with
a controlled tail*, which is what `free_space_recovery.py`,
`green_derivative_tail_bounds.md` and `HS-4` are for, and which is a genuinely
different technique with genuinely different failure modes.  That, rather than any
particular missing constant, is why none of these five is the template.
