# Track P: a rigorous a posteriori lane on the periodic torus

**Status: a derivation document with every constant exposed.**  It records the
H⁴ control inequality that Track P runs on, with each constant reduced to a
lattice sum, a binomial identity or a Cauchy–Schwarz step that a referee can
re-check on paper.  Nothing here uses an external numerical constant: the
Morosi–Pizzocchero machinery analysed in
[`a_posteriori_frameworks.md`](a_posteriori_frameworks.md) is the comparison
point, not an ingredient.  Nothing here is a singularity result, and nothing
here bears on the Clay problem; §6 states exactly what is concluded and under
which named external theorems.  Viscosity is fixed and positive throughout, as
everywhere in this repository.

Notation, fixed once.  The domain is `T³ = (ℝ/2πℤ)³` with the **normalised**
measure `(2π)⁻³ dx`.  Fields are real, mean-zero and divergence-free.  Fourier
expansion `f = Σ_k f̂_k e^{ik·x}` over `k ∈ ℤ³`, Parseval
`⟨f, g⟩ = Σ_k f̂_k · conj(ĝ_k)`, and the **homogeneous** norms

```
‖f‖_n² = Σ_k |k|^{2n} |f̂_k|² .
```

Mean zero means the sum runs over `k ≠ 0`, so `|k| ≥ 1` and the norms are
monotone: `‖f‖_m ≤ ‖f‖_n` for `m ≤ n`.  That single inequality — the spectral
gap — is the reason this lane exists.

## 1. Why a periodic lane, honestly

The Clay problem statements (B) and (D) are the periodic ones, so a periodic
lane is not a detour from the official problem; it is one of its two official
settings.  What makes the torus worth a dedicated track is that **both**
whole-space obstructions identified in
[`a_posteriori_frameworks.md`](a_posteriori_frameworks.md) disappear here, and
they disappear for structural reasons rather than by extra work.

First, the spectral gap.  On `ℝ³` the control ODE loses its `-ν R` damping term
because the Laplacian's spectrum reaches down to zero (§3 of that note, "The
`-ν R_n` term, in detail").  On the torus with mean-zero fields, `|k| ≥ 1`
gives `‖∇w‖_n ≥ ‖w‖_n`, the viscous term donates a genuine `-ν R`, and the
Riccati structure recovers the stable equilibrium that every
global-existence-from-smallness corollary rests on.

Second, and more important for this repository: **the residual is exactly
computable.**  The approximate solution on this lane is a Galerkin trajectory —
a finite trigonometric polynomial at every instant — and the nonlinearity of a
band-limited field is band-limited.  Its *continuous* PDE residual is therefore
a finite Fourier polynomial whose coefficients are computed exactly in rational
arithmetic.  There is no spatial interpolation, no grid, no
discrete-to-continuous recovery step.  The lemma that
[`hs_error_propagation.md`](hs_error_propagation.md) §4 names HS-5 — turning a
nodal residual into a Sobolev norm of the true residual, "the largest single
gap" — **closes on this lane by construction**, because the measured object and
the analytical object are the same finite list of rational numbers.

Stated plainly, because it is the load-bearing claim of the whole track: on
Track P the approximate solution is an exact trig polynomial, and its
continuous PDE residual is computed exactly in Fourier coefficients.  The `ℝ³`
lane keeps its HS-5 gap, its truncation budget (HS-4) and its missing spectral
gap.  The two lanes must never be conflated, and no result from this note may
be quoted for the whole-space programme.

One clarification about the verdict of
[`a_posteriori_frameworks.md`](a_posteriori_frameworks.md) §4, "do not port
Morosi–Pizzocchero".  That verdict was about porting the framework to `ℝ³`, and
it stands.  Track P is not a port and not an application of their constants: it
is a self-contained rederivation of the same *shape* of argument on the domain
where the shape is native, using deliberately crude constants that this
repository can verify from scratch (§5.6 records the cost).

## 2. Fixed band versus finite band: what the no-go does and does not exclude

The repository already holds a finite-mode exclusion theorem:
[`track_f_finite_mode_nogo.md`](track_f_finite_mode_nogo.md) proves (Theorem 1
and Corollary 2, with the abstract core formalised in
`formal/NSSingularity/GalerkinNoBlowup.lean` and
`formal/NSSingularity/FiniteModeNoGo.lean`) that a trajectory which **remains
in a fixed finite Fourier band for all time** cannot blow up: the energy bound
plus finite-dimensional norm equivalence extends it past any putative singular
time.

It is essential not to misread that theorem as covering Track P's data.  The
hypothesis is that the *trajectory* stays in the band.  A **finite-band initial
datum** does not satisfy it: the moment the nonlinearity acts, the true
solution populates modes outside the band — the Galerkin tail `e` of §5.3 is
generically nonzero, and it is precisely the part of `P(u_a·∇u_a)` that the
band cannot hold.  So the true solution launched from a trig-polynomial datum
leaves the band immediately, the no-go says nothing about it, and no
circularity arises when Track P certifies facts about that solution.

The distinction is being formalised in parallel:
`formal/NSSingularity/TrackPFourier.lean` defines `FixedBandTrajectory` (the
excluded object: band-limited for all `t`) and `FiniteBandDatum` (Track P's
object: band-limited at `t = 0` only), and exhibits the counterexample showing
the second does not imply the first — a datum whose exact solution has, for
arbitrarily small `t > 0`, energy outside the initial band.  Track P works in
exactly the unexcluded region: finite-band datum, infinite-band true solution,
and the distance between the true solution and the band-limited surrogate
controlled by the inequality of §5.

## 3. The data

The initial data are real, mean-zero, divergence-free finite trigonometric
polynomials with rational coefficients.  A finite trig sum is `C^∞`, so these
data are legitimate periodic Clay-class data as they stand — no mollification,
no closure, no surrogate status.  (Contrast the `C^k` spline rule recorded in
[`gaussian_hermite_basis.md`](gaussian_hermite_basis.md): splines are tools,
never candidates.  Trig polynomials need no such rule because they are already
smooth.)

Three families, whose concrete wavevectors and rational amplitudes live in
`src/ns_certificate_lab/fourier_torus.py` (written in parallel with this note;
this section records intent, the module records numbers):

* **P1, the helical triad.**  Three wavevectors summing to zero, each carrying
  a helical amplitude orthogonal to its wavevector.  The minimal genuinely
  nonlinear periodic datum: the convolution of the triad with itself is small,
  finite and analytically transparent, so every quantity in §5 can be audited
  by hand on this family before being trusted on the others.
* **P2, two connected triads sharing a mode.**  The minimal cascade coupling:
  energy has a path from one triad to the other through the shared leg, the
  interaction set `(G+G) \ G` is materially richer, and the Galerkin tail
  stops being a single-shell object.  This family is what exercises the
  residual bookkeeping.
* **P3, symmetry-broken and pressure-active.**  Amplitudes chosen so that the
  gradient part `(I−P)(u_a·∇u_a)` is nonzero, and *measured* to be nonzero.
  Gate 8 taught the repository that a symmetry can silently kill the pressure
  channel and make a whole search vacuous
  ([`l3_positive_generation_search.md`](l3_positive_generation_search.md) §2);
  P3 exists so that the periodic lane provably does not live in the
  pressure-free subspace.

## 4. The exact operations

Everything below is rational arithmetic on finitely many coefficients.  That is
the entire point of the lane, so the operations are listed with the reason each
one is exact.

**Representation.**  A field is stored as cos/sin pairs over a half-lattice of
represented modes: `u = Σ_{k∈G⁺} (a_k cos(k·x) + b_k sin(k·x))` with rational
vectors `a_k, b_k ∈ ℚ³`.  Reality is built into the representation;
divergence-freeness and mean-zeroness are per-mode rational linear constraints
(`k·a_k = k·b_k = 0`, `0 ∉ G⁺`).

**Convolution by product-to-sum.**  Products of modes reduce by the exact
identities `cos A cos B = ½[cos(A−B) + cos(A+B)]`,
`sin A sin B = ½[cos(A−B) − cos(A+B)]`,
`sin A cos B = ½[sin(A−B) + sin(A+B)]`.  The bilinear term of a band-limited
field is computed **exactly**, term by term, in rationals.  There is no grid,
hence no aliasing and no dealiasing error: the quantity that a pseudospectral
computation would call dealiasing error appears here instead as the
exactly-computed Galerkin tail of §5.3 — named, bounded and carried through the
estimate rather than filtered away.

**Leray projection.**  A per-mode rational matrix: `P_k = I − k kᵀ/(k·k)`.
Orthogonal projection in every `‖·‖_n` inner product, commuting with all
derivatives, self-adjoint mode by mode.

**Norms and invariants.**  With the normalised measure,
`‖u‖_n² = ½ Σ_{k∈G⁺} (k·k)ⁿ (|a_k|² + |b_k|²)` — a rational quadratic form for
every integer `n`.  Energy `½‖u‖₀²` and enstrophy `‖u‖₁²` are the `n = 0, 1`
cases.  The `ℓ¹`-type majorants `M_j`, `N_m` of §5.4 are rational sums over the
same coefficients.

**Time-Taylor coefficients.**  The Galerkin right-hand side is quadratic, so
the Taylor coefficients of the trajectory obey the exact recursion

```
(m+1) u^{(m+1)} = ν Δ u^{(m)} − Σ_{i+j=m} P_G P(u^{(i)}·∇u^{(j)}) ,
```

evaluated in rationals.  This recursion integrates the trajectory (an interval
Taylor step whose Lagrange remainder is evaluated on the Picard rough box, the
scalar template of `control_ode.py` lifted to the coefficient vector) and,
run once more on the box, encloses every slab quantity of §6.

## 5. The H⁴ control inequality

### 5.1 The embedding constant `A`

For mean-zero `f`, Cauchy–Schwarz against the weight `|k|⁻²` gives the sup
bound

```
‖f‖_∞ ≤ Σ_{k≠0} |f̂_k| = Σ_{k≠0} |k|⁻² · |k|²|f̂_k|
      ≤ ( Σ_{k≠0} |k|⁻⁴ )^{1/2} ‖f‖₂  =:  A ‖f‖₂ .
```

`A² = Σ_{k≠0} |k|⁻⁴` is bounded rigorously by an exact lattice sum plus a tail.
The shell `|k|_∞ = m` contains `(2m+1)³ − (2m−1)³ = 24m² + 2 ≤ 26 m²` points,
and `|k| ≥ |k|_∞`, so

```
Σ_{|k|_∞ = m} |k|⁻⁴ ≤ (24m² + 2) m⁻⁴ ≤ 26 m⁻² ,
Σ_{|k|_∞ > N} |k|⁻⁴ ≤ Σ_{m>N} 26 m⁻² ≤ 26/N .
```

With the exact rational sum over `|k|_∞ ≤ 40` (evaluated in
`fractions.Fraction`; its float image is `16.27443269…`) this gives

```
A² ≤ 16.2745 + 26/40 = 16.9245 ,      A ≤ 4.114 ,
```

and the partial sum alone shows `A ≥ 4.034`, so the bound is tight to about
2%.  A referee can reproduce the whole constant with a thirty-line script and
the two displayed inequalities; nothing else enters.

### 5.2 The inner product and the multinomial weights

The working inner product is the homogeneous H⁴ form written through fourth
derivatives:

```
⟨f, g⟩₄ = Σ_{|α|=4} c_α ⟨∂^α f, ∂^α g⟩ ,      c_α = 4!/α! ,
```

which reproduces `‖f‖₄²` exactly because
`Σ_{|α|=4} c_α k^{2α} = (k₁² + k₂² + k₃²)⁴ = |k|⁸` (the multinomial theorem).
Two numbers from this choice recur below and are worth naming now:

```
Σ_{|α|=4} c_α = 3⁴ = 81           (set every k_i² = 1) ,
```

and its square root `9`, which is the price of Cauchy–Schwarz over the index
`α`: for any nonnegative `y_α`,

```
Σ_α c_α y_α ≤ (Σ_α c_α)^{1/2} (Σ_α c_α y_α²)^{1/2} = 9 (Σ_α c_α y_α²)^{1/2} .
```

The other combinatorial identity is the multi-dimensional Vandermonde
convolution: for `|α| = 4`,

```
Σ_{|β|=j, β≤α} binom(α,β) = binom(4,j) ,
```

which collapses every Leibniz sum below into four binomial weights, and
`Σ_{j=1}^{4} binom(4,j) = 2⁴ − 1 = 15`.

### 5.3 The error equation and the exact residual

Let `u` be the true strong solution (§6 supplies it), `u_a` the Galerkin
trajectory on the symmetric band `G`, and `w = u − u_a`.  Both are mean-zero
and divergence-free.  Subtracting the projected equations and splitting the
bilinear difference exactly (`(u_a+w)·∇(u_a+w) − u_a·∇u_a
= u_a·∇w + w·∇u_a + w·∇w`),

```
∂_t w = ν Δw − P(w·∇u_a) − P(u_a·∇w) − P(w·∇w) − e ,

e := ∂_t u_a − ν Δu_a + P(u_a·∇u_a) = (I − P_G) P(u_a·∇u_a) .
```

The second line is the crucial identity: because `u_a` solves the Galerkin
system `∂_t u_a = ν Δu_a − P_G P(u_a·∇u_a)` exactly, its continuous residual is
**exactly** the Galerkin tail — a finite Fourier polynomial supported on the
modes of `(G+G) \ G`, computed exactly in rational arithmetic by the
product-to-sum convolution of §4.  This is the sentence in which HS-5 closes on
the periodic lane: the residual is not estimated, sampled or recovered; it is
written down.

### 5.4 The energy estimate, term by term

Pair the error equation with `w` in `⟨·,·⟩₄`.  Throughout, `P` drops out at no
cost: `P` is self-adjoint and commutes with derivatives, so
`⟨P f, w⟩₄ = ⟨f, P w⟩₄ = ⟨f, w⟩₄` since `w` is divergence-free.

**Viscous term — helps, exactly.**
`⟨ν Δw, w⟩₄ = −ν ‖w‖₅² ≤ −ν ‖w‖₄²` by the spectral gap.  This is the term the
whole-space lane does not have.

**Transport term `⟨u_a·∇w, w⟩₄` — commutator only.**  Leibniz:
`∂^α(u_a·∇w) = Σ_{β≤α} binom(α,β) (∂^β u_a)·∇ ∂^{α−β} w`.  The `β = 0` term
pairs to zero: `⟨u_a·∇∂^α w, ∂^α w⟩ = −½∫(div u_a)|∂^α w|² = 0`.  For
`|β| = j ≥ 1`, Cauchy–Schwarz over the contracted index and Hölder
`L^∞ × L² × L²` give

```
|⟨(∂^β u_a)·∇∂^{α−β}w, ∂^α w⟩| ≤ ‖∂^β u_a‖_∞ · ‖∇∂^{α−β}w‖ · ‖∂^α w‖ ,
```

with `‖∇∂^{α−β}w‖ ≤ ‖w‖_{5−j} ≤ ‖w‖₄` (Parseval; `5−j ≤ 4`).  For a finite
Fourier field the sup is majorised by the exact rational sum

```
M_j := Σ_{k∈G⁺} |k|_∞^j (|a_k| + |b_k|)  ≥  max_{|β|=j} ‖∂^β u_a‖_∞ ,
```

using `|k^β| ≤ |k|_∞^{|β|}` mode by mode.  Vandermonde collapses the `β`-sum,
Cauchy–Schwarz over `α` supplies the `9`:

```
|⟨u_a·∇w, w⟩₄| ≤ 9 K₁ ‖w‖₄² ,      K₁ := Σ_{j=1}^{4} binom(4,j) M_j .
```

**Stretching term `⟨w·∇u_a, w⟩₄` — no cancellation, fifth derivatives.**  The
same Leibniz expansion, but now every term including `β = 0` survives, and the
sup falls on the *matrix* `∇∂^{α−β}u_a`.  For `|β| = j`, `0 ≤ j ≤ 4`,

```
|⟨(∂^β w)·∇∂^{α−β}u_a, ∂^α w⟩| ≤ ‖ |∇∂^{α−β}u_a|_F ‖_∞ · ‖∂^β w‖ · ‖∂^α w‖ ,
```

with `‖∂^β w‖ ≤ ‖w‖_j ≤ ‖w‖₄`.  The Frobenius contraction is the one place a
Euclidean factor of `k` is unavoidable: the mode-`k` contribution to
`∇∂^γ u_a` is the rank-one matrix `k ⊗ (k^γ a_k)` (times a trig function), with
Frobenius norm `|k| |k^γ| |a_k|`, and `|k| ≤ |k|_∞` fails off the axes by up to
`√3`.  The correct — still exactly computable, still rational — majorant
carries one `ℓ¹` factor:

```
N_m := Σ_{k∈G⁺} |k|₁ |k|_∞^{m−1} (|a_k| + |b_k|)
     ≥ max_{|γ|=m−1} ‖ |∇∂^γ u_a|_F ‖_∞ ,
```

using `|k| ≤ |k|₁` and `|k^γ| ≤ |k|_∞^{|γ|}`.  (Any implementation that
majorises this constant with `|k|_∞` powers alone is wrong by up to `√3`;
`N_m ≤ 3 M_m`, so the correction changes no structure, only a constant.)  With
`m = 5−j` — note `j = 0` needs `N₅`, i.e. **fifth** derivatives of `u_a`,
exactly available for finite Fourier data — the same Vandermonde and `9` give

```
|⟨w·∇u_a, w⟩₄| ≤ 9 K₂ ‖w‖₄² ,      K₂ := Σ_{j=0}^{4} binom(4,j) N_{5−j} .
```

**Cubic term `⟨w·∇w, w⟩₄` — closes with the `A` embedding alone.**  The
`β = 0` term again pairs to zero (`div w = 0`).  Each surviving commutator term
with `|β| = j` is a product of two `w`-factors and `∂^α w`, and at `n = 4`
every one of them closes with **only** the embedding of §5.1:

* `j = 1, 2`: put the sup on the low factor.
  `‖∂^β w‖_∞ ≤ Σ_k |k^β||ŵ_k| ≤ Σ_k |k|^j |ŵ_k| ≤ A ‖w‖_{j+2} ≤ A ‖w‖₄`
  (the §5.1 Cauchy–Schwarz with weight `|k|⁻²`, and `j+2 ≤ 4`),
  and `‖∇∂^{α−β}w‖ ≤ ‖w‖_{5−j} ≤ ‖w‖₄`.
* `j = 3, 4`: put the sup on the gradient factor.
  `‖ |∇∂^{α−β}w|_F ‖_∞ ≤ Σ_k |k|^{5−j}|ŵ_k| ≤ A ‖w‖_{7−j} ≤ A ‖w‖₄`
  (since `7−j ≤ 4`), and `‖∂^β w‖ ≤ ‖w‖_j ≤ ‖w‖₄`.

Every term is `≤ A ‖w‖₄² · ‖∂^α w‖`; Vandermonde sums the weights to
`Σ_{j=1}^4 binom(4,j) = 15`, and Cauchy–Schwarz over `α` supplies the `9`:

```
|⟨w·∇w, w⟩₄| ≤ 9 · 15 · A ‖w‖₄³ = 135 A ‖w‖₄³ ≤ 555.4 ‖w‖₄³ .
```

**Residual term.**  `|⟨e, w⟩₄| ≤ ‖e‖₄ ‖w‖₄`, Cauchy–Schwarz in the H⁴ inner
product, with `‖e‖₄` the exact rational norm of the exactly computed tail.

**Assembled.**  With `W(t) = ‖w(t)‖₄` and the `ε`-regularisation
`y_ε = (W² + ε)^{1/2}` handling the non-differentiability of the norm at zero
(the same device as [`track_f_finite_mode_nogo.md`](track_f_finite_mode_nogo.md)
§3, proof of (i)),

```
d⁺W/dt ≤ ( −ν + 9(K₁ + K₂) ) W + 135 A W² + ‖e(t)‖₄ .
```

### 5.5 Why `n = 4` and not `n = 3`

Run the cubic term's bookkeeping at `n = 3`.  The commutator terms are
`(∂^β w)·∇∂^{α−β}w` with `|α| = 3`, `j = 1, 2, 3`.  The `j = 1` term closes
(`‖∂w‖_∞ ≤ A‖w‖₃`, gradient factor `‖w‖₃`); so does `j = 3` (`‖∂³w‖ = ‖w‖₃`,
`‖∇w‖_∞ ≤ A‖w‖₃`).  The middle term `j = 2` does not: putting the sup on
`∂²w` needs `A‖w‖₄`, putting it on `∇∂w` also needs `A‖w‖₄`, and both exceed
the working level.  The classical escape is the Hölder split `L³ × L⁶` plus
Gagliardo–Nirenberg, which is exactly the sharp machinery Morosi and
Pizzocchero built: their Kato constant `G₃ = 0.438` (recorded, flagged as
**[recalled]**, in [`a_posteriori_frameworks.md`](a_posteriori_frameworks.md)
§1.5) is the optimised form of precisely this step.  At `n = 4` the middle
terms acquire one derivative of slack on each side and the crude `ℓ¹` embedding
closes everything.  That is the entire reason Track P works at `n = 4`:
self-containedness is bought with one extra derivative.

### 5.6 What the crude constants cost, on the record

The cost is not hidden and it is not an oversight.  Three line items:

* The cubic coefficient is `135 A ≤ 555.4`, against Morosi–Pizzocchero's
  `G₃ = 0.438` at `n = 3` — three orders of magnitude.  The blow-up time of a
  Riccati inequality scales inversely with its quadratic coefficient, so at
  comparable residuals Track P's certified slabs are of order a thousand times
  shorter than what the sharp constants would certify.
* `K₁, K₂` majorise sups by `ℓ¹` sums of coefficients; for concentrated
  spectra this loses a further modest factor, and `N_m ≤ 3 M_m` records the
  `ℓ¹`-on-one-index correction of §5.4.
* Working at `n = 4` weights the residual by `|k|⁸`.  The tail lives on
  `(G+G) \ G`, whose modes reach `2·max|k|`, so choosing `n = 4` over `n = 3`
  inflates `‖e‖_n` by up to `(2 max|k|)²` relative.

All three are the recorded price of a derivation in which **every constant is a
lattice sum, a binomial coefficient or a Cauchy–Schwarz step** — checkable by a
referee in an afternoon, with no dependence on constants this repository has
only ever seen as [recalled] values.  Sharpening any of the three is a quality
upgrade that changes no structure; the trade was made knowingly, interval
length for auditability.

## 6. The control ODE and the conclusion

**The slab enclosure.**  Fix a slab `[0, h]`.  The Galerkin trajectory itself
needs no external theorem: its coefficient vector obeys a polynomial ODE, and
the Stage-B self-mapping inclusion of
`src/ns_certificate_lab/control_ode.py` (`prove_rough_enclosure`, the
Picard–Bielecki theorem recorded verbatim in `PICARD_SELF_MAPPING`) lifted to
the coefficient vector proves existence on the whole slab inside an explicit
rough box.  Interval evaluation of the rational formulas of §4 and §5.4 on that
box yields slab-uniform interval bounds for `M_j`, `N_m` (hence `K₁, K₂`) and
`‖e(t)‖₄` — the interval coefficients that `ControlCoefficients` expects.

**The Riccati inequality.**  With `a ⊇ −ν + 9(K₁+K₂)`, `b ⊇ 135 A` and
`ε ⊇ sup_{[0,h]} ‖e(t)‖₄` as intervals, the control ODE

```
d⁺R/dt = a R + b R² + ε ,      R(0) = 0
```

is certified by `build_control_certificate` exactly as in the scalar module:
rough box, then Bernstein-hulled Taylor tube.  **`R(0) = 0` is not an
assumption but a construction: the Galerkin datum IS the datum.**  The initial
condition of the true solution is defined to be `u_a(0)`, so the datum error is
identically zero — the `R_embed = 0` situation of
[`hs_error_propagation.md`](hs_error_propagation.md) §4, here exact rather than
grid-exact.

**The comparison.**  `W(t) = ‖w(t)‖₄` satisfies the differential inequality of
§5.4 with `≤`; `R` satisfies the equation with these (interval-dominating)
coefficients; `R(0) = 0 ≥ W(0) = 0`.  The Chaplygin/Dini comparison lemma —
the domain-independent part of the Morosi–Pizzocchero mechanism, §1.3 of
[`a_posteriori_frameworks.md`](a_posteriori_frameworks.md) — gives `W ≤ R` for
as long as both exist.

**The named external theorems.**  The argument consumes three classical facts
about periodic Navier–Stokes, recorded here faithfully and by name.  They are
**never inserted as Lean axioms** — the rule of `LEAN4_VERIFICATION_POLICY.md`
and of [`track_f_finite_mode_nogo.md`](track_f_finite_mode_nogo.md) §6 applies
verbatim: a citation is a citation, and anything on a Lean final path needs
either a faithful formalisation or an explicitly audited theorem interface.

> **EXT-P1 (local existence and uniqueness).**  For mean-zero, divergence-free
> `u₀ ∈ H⁴(T³)` and fixed `ν > 0` there exist `T* > 0` and a unique strong
> solution `u` of the Navier–Stokes equations on `[0, T*)` with
> `u(0) = u₀`, `u ∈ C([0,T]; H⁴)` for every `T < T*`.  Source: the Kato-style
> local theory; Temam, *Navier–Stokes Equations: Theory and Numerical
> Analysis*, Ch. III (the "Thm 3.7-ish" locus — the theorem *number* is
> [recalled], not re-checked this session), or Majda–Bertozzi, *Vorticity and
> Incompressible Flow*, Ch. 3, whose `H^m`, `m > 5/2` local theory covers
> `m = 4` on the torus.
>
> **EXT-P2 (regularity to run the estimate).**  That strong solution has
> `u ∈ C([0,T]; H⁴) ∩ C¹([0,T]; H²)`, enough regularity to justify the pairing
> `d/dt‖w‖₄² = 2⟨∂_t w, w⟩₄` (by the standard mollification/Galerkin
> justification of energy identities at the top level) and the Dini-derivative
> manipulation of §5.4.  Standard within the same sources.
>
> **EXT-P3 (continuation criterion).**  The solution of EXT-P1 extends beyond
> `T` whenever `sup_{t<T} ‖u(t)‖_{H⁴} < ∞`: the local existence time of EXT-P1
> depends only on the `H⁴` norm of the datum, so a bounded norm allows the
> solution to be restarted past `T`.  Standard; same sources.

**The conclusion.**  Suppose the slab certificate closes: the Picard box for
the Galerkin trajectory holds on `[0, h]` and the control-ODE tube `R(t)` is
finite on `[0, h]`.  Then, conditional on EXT-P1/P2/P3 and on nothing else:
the true strong solution with datum `u_a(0)` exists on the whole of `[0, h]`
(EXT-P3 plus `‖u‖₄ ≤ ‖u_a‖₄ + R < ∞` forbids earlier breakdown), and

```
‖u(t) − u_a(t)‖_{Ḣ⁴} ≤ R(t)          for all t ∈ [0, h] ,
‖u(t) − u_a(t)‖_{H³} ≤ √8 · R(t) ≤ 2.8285 · R(t) ,
```

the second line because `(1 + |k|²)³ ≤ 8|k|⁶ ≤ 8|k|⁸` for `|k| ≥ 1`.  This is
an existence-and-proximity statement over one short slab.  It is **not** a
singularity result, it is the opposite of one (§8), and it says nothing about
the Clay problem.

## 7. The mapping table

The same table as [`a_posteriori_frameworks.md`](a_posteriori_frameworks.md)
§2, now with the periodic lane's entries filled in.  `fourier_torus.py` is
landing in parallel with this note; its rows name the module's objects by role.

| Morosi–Pizzocchero object | Track P counterpart | where |
|---|---|---|
| approximate solution `u_a` | the Galerkin trajectory: band-limited coefficient ODE integrated by the exact Taylor recursion | `src/ns_certificate_lab/fourier_torus.py` |
| datum error `δ_m` | **zero by construction** — the Galerkin datum is the datum | §6; no code needed |
| differential residual `e(u_a)`, `‖e‖_m ≤ ε_m` | the **exact** Galerkin tail `(I−P_G)P(u_a·∇u_a)`, finite support `(G+G)\G`, rational coefficients; `‖e‖₄` an exact rational form | `fourier_torus.py` (tail assembly via the product-to-sum convolution) |
| growth estimator `D_m(t)` | slab-uniform interval bounds `M_j`, `N_m` evaluated on the Picard rough box of the coefficient ODE | `fourier_torus.py` box + §5.4 formulas |
| Sobolev norm estimator | the exact rational `H^n` quadratic forms `½Σ(k·k)ⁿ(\|a_k\|²+\|b_k\|²)` | `fourier_torus.py` |
| the constants `G_n`, `K_n` | **not used.**  Replaced by the self-contained `A ≤ 4.114`, `9 = √81`, `binom(4,j)`, `K₁`, `K₂`, `135A`; MP's `G₃ = 0.438` recorded for comparison only | §5; [`a_posteriori_frameworks.md`](a_posteriori_frameworks.md) §1.5 |
| control inequality / control ODE | `R' = aR + bR² + ε`, `R(0)=0`, certified rough-box-then-Taylor-tube | `src/ns_certificate_lab/control_ode.py`: `ControlCoefficients`, `prove_rough_enclosure`, `taylor_tube`, `build_control_certificate`, `verify_control_certificate` |
| existence interval `T_c` | the slab `[0, h]`, iterable in principle (one slab in v1, §8) | §6 |
| exact–approximate distance `‖u−u_a‖_n ≤ R_n` | `R(t)`, the refined Taylor tube; `H³` via the factor `√8` | `control_ode.py` payload, §6 |

## 8. Limitations

**Short slabs.**  §5.6 quantifies it: the quadratic coefficient is roughly a
thousand times Morosi–Pizzocchero's, and Riccati blow-up time scales inversely
with it.  The slabs this lane certifies will be short, and their brevity is the
recorded, deliberate price of constants a referee can check from nothing.

**One slab in v1.**  Chaining slabs is not a loop over §6: restarting from the
endpoint interval of one tube as the datum interval of the next loses the shape
of the reachable set, and the wrapping effect compounds exponentially across
slabs.  Doing it honestly needs Lohner-type moving frames — the machinery
already cited in `control_ode.py`'s Moore–Lohner discussion, and the same
dependency/wrapping-effect suppression obligation recorded against interval
verification in `docs/proof_obligations.md` PO-13.  Until that machinery
exists, v1 certifies one slab and says so.  Note also that `R(0) = 0` holds
only on the *first* slab; every subsequent slab starts from a nonzero interval,
so chaining changes the datum-error row of §7 as well as the frame problem.

**Small Galerkin cutoff.**  The families of §3 live on small bands, both to
keep the exact convolution auditable and because `‖e‖₄`'s `|k|⁸` weighting
punishes wide bands.  A small band means a large tail unless the datum is
gentle or the viscosity strong; the lane's first honest numbers will reflect
that.

**The direction of the certificate, stated so it cannot be misquoted.**  A
certificate of existence-near-a-trajectory is the **opposite** of a singularity
certificate.  Track P certifies *regularity, locally in time*: the true
solution exists on the slab and stays near a smooth surrogate.  That is the
correct first rigorous step for any a posteriori programme — it is what
Morosi–Pizzocchero's framework does, and what every framework surveyed in
[`a_posteriori_frameworks.md`](a_posteriori_frameworks.md) §5 does — and it
must not be oversold.  No accumulation of Track P slabs can certify a blow-up:
near a putative singular time the constants `K₁, K₂, ‖e‖₄` grow and the slabs
shrink to nothing, exactly as §5 of
[`hs_error_propagation.md`](hs_error_propagation.md) explains for the
whole-space estimate.  Nothing in this note is a claim about the Clay problem,
in either direction.
