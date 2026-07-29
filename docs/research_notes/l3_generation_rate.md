# The generation rate of the critical `L³` norm, and the pure-swirl no-go

**Status: a complete derivation of a classical identity, plus one new-to-this-repository
corollary and a search objective.**  Nothing here is new mathematics.  The identity of §1
is the `p = 3` member of the classical `L^p` energy identity for Navier–Stokes — the same
computation that stands behind the endpoint theory of Escauriaza, Seregin and Šverák and
behind every textbook proof that `‖u(t)‖_{L³}` cannot grow under a smallness hypothesis.
It is written out here in full only because the repository needs it in a form a referee
can check line by line, with its hypothesis class stated honestly rather than waved at,
and because the corollary in §4 kills a search lane that this repository had been running
for eighteen computations.

Every step below is marked *proved* or *required lemma*.  Every required lemma carries a
name, a precise statement, and the constant it would supply.  §8 says bluntly what none of
this establishes, and §9 collects the lemmas.  Nothing here bears on the Clay problem, and
no claim about it is made anywhere in this note.

---

## 1. Statement and hypothesis class

### 1.1 The identity

Let `ν > 0` be fixed and let `u` solve the incompressible Navier–Stokes equations on `ℝ³`,

```
∂_t u + (u·∇)u + ∇p = νΔu ,        ∇·u = 0 ,
```

on a compact time interval `J`.  Put

```
F(t) = ∫_{ℝ³} |u(x,t)|³ dx = ‖u(t)‖_{L³}³ .
```

> **Theorem (the `L³` generation identity).**  Under the hypothesis class `H` of §1.2,
> `F ∈ C¹(J)` and
> ```
> F'(t) = 3∫ p ∇·(|u|u) dx  -  3ν∫ ( |u| |∇u|² + |u| |∇|u||² ) dx .      (1.1)
> ```
> Equivalently, in the form the numerics use,
> ```
> F'(t) = -3∫ |u| (u·∇p) dx  -  3ν∫ ( |u| |∇u|² + |u| |∇|u||² ) dx .     (1.2)
> ```

Write `P` for the first integral (with its factor `3`) and `V` for the second (with its
factor `-3ν`), so `F' = P + V`.

This is the `p = 3` case of the family

```
d/dt ∫|u|^p = p∫ p_press ∇·(|u|^{p-2}u) dx
              - pν∫ ( |u|^{p-2}|∇u|² + (p-2)|u|^{p-2}|∇|u||² ) dx ,
```

valid for `p ≥ 2` under the same kind of hypotheses, where `p_press` denotes the pressure —
written out in full here only because this is the one display in which the exponent `p` and
the pressure would otherwise collide.  Two features of the exponent are
worth naming because they are exactly why `p = 3` is the interesting case.  At `p = 2` the
pressure term is `2∫ p_press ∇·u = 0` and the identity degenerates to the energy identity;
`p = 3` is the first exponent at which the pressure genuinely enters.  And `p = 3` is the
exponent at which `‖u‖_{L^p}` is invariant under the Navier–Stokes scaling
`u ↦ λu(λx, λ²t)`, which is what makes `F` a critical quantity rather than a
units-dependent one.  The coefficient `(p-2)` equals `1` at `p = 3`, which is why the two
viscous integrands appear with equal weight.

### 1.2 The hypothesis class `H`

```
(H1)  u ∈ C(J; H³(ℝ³)) ∩ C¹(J; L²(ℝ³)) ;
(H2)  ∇·u(·,t) = 0 pointwise for every t ∈ J ;
(H3)  |u(x,t)| + |∇u(x,t)| ≤ C_J ⟨x⟩^{-4}  for all x ∈ ℝ³ and all t ∈ J ,
      with one constant C_J for the whole compact interval  ( ⟨x⟩ = (1+|x|²)^{1/2} ) ;
(H4)  p = R_i R_j (u_i u_j) , the decaying Riesz pressure.
```

`(H1)` gives enough regularity for the equation to hold in `L²` at every time and for the
integrations by parts of §2 to be legitimate; `H³(ℝ³) ↪ C^{1,1/2}` so `u` and `∇u` are
genuinely continuous functions, and `Δu ∈ H¹ ⊂ L²`.

`(H3)` is the working hypothesis and it is the one that has to be argued for.  It is used
in exactly four places: to put `u(t)` in `L¹` (which `(H1)` alone does *not* give — `H³` is
contained in `L² ∩ L^∞` and hence in every `L^q` with `q ≥ 2`, but not in `L¹`); to make
`|u||∇u|²` and `(1+|u|)|∇u|²` integrable; to kill the three boundary terms of §2; and to
run the equality case of §4.3.  The `L¹` membership is the one that would be easy to
overlook, and `⟨x⟩^{-4}` is exactly the borderline that supplies it: `∫_{ℝ³}⟨x⟩^{-4}dx`
converges because `4 > 3`, while `⟨x⟩^{-3}` would not.

### 1.3 Why "Schwartz" is the wrong class

It is tempting to state (1.1) for Schwartz solutions and be done with it.  That would be a
statement about the empty set at positive times, and the corollary of §4 would be a
statement about the empty set at positive times *for exactly the data it is about*.

**The Schwartz class is not propagated by Navier–Stokes.**  Take `u_0` smooth,
divergence free and compactly supported.  The time derivative at `t = 0` is

```
∂_t u |_{t=0} = νΔu_0 - (u_0·∇)u_0 - ∇p_0 .
```

The first two terms are compactly supported.  The third is not, and its tail is exactly
computable.  Since `-Δp_0 = ∂_i u_j ∂_j u_i = ∂_i∂_j(u_{0,i}u_{0,j})`, the source of the
pressure Poisson equation is a second derivative, so its low moments cancel:

```
∫ ∂_i∂_j(u_iu_j) dx = 0 ,              (monopole: a total derivative)
∫ x_k ∂_i∂_j(u_iu_j) dx = 0 ,          (dipole:  ∂_i∂_j x_k = 0)
∫ x_k x_l ∂_i∂_j(u_iu_j) dx
   = ∫ (δ_{ik}δ_{jl} + δ_{il}δ_{jk}) u_iu_j dx = 2 M_{kl} ,   M = ∫ u⊗u dx .
```

So the leading far-field term of the pressure is the quadrupole,

```
p_0(x) = (1/4π|x|³) ( 3 M_{ij} n_i n_j - tr M ) + O(|x|^{-4}) ,   n = x/|x| ,
```

and this vanishes for every direction `n` if and only if `M` is a multiple of the
identity.  Consequently `∇p_0 = O(|x|^{-4})` and **not** `o(|x|^{-4})` unless
`M ∝ I`.  For small `t > 0` the solution therefore carries a tail
`u(t) ≈ -t∇p_0 ~ |x|^{-4}`: compact support is destroyed instantly, and so is Schwartz
decay.  This is the Brandolese localisation obstruction, and `|x|^{-4}` is its sharp
exponent in three dimensions — which is precisely why `(H3)` is stated with that exponent
and not with a better one.

**The sharp instance is the very datum §4 is about.**  Let `u_0 = u^θ(r,z)e_θ` be a
nonzero pure swirl.  In Cartesian components `u_0 = u^θ(-sinθ, cosθ, 0)`, so

```
M_{11} = ∫ (u^θ)² sin²θ dx = ½∫(u^θ)² dx ,
M_{22} = ∫ (u^θ)² cos²θ dx = ½∫(u^θ)² dx ,
M_{12} = -∫ (u^θ)² sinθ cosθ dx = 0 ,      M_{i3} = 0 ,
```

the off-diagonal and third-row entries vanishing on the `θ` integration.  Hence

```
M = diag(a, a, 0) ,     a = ½∫(u^θ)² dx > 0 ,
```

which is **not** a multiple of the identity, and

```
3 M_{ij}n_in_j - tr M = a(1 - 3n_3²) ,
```

nonzero except on the cone `n_3² = 1/3`.  So a pure swirl generates an exact `|x|^{-3}`
pressure tail and an exact `|x|^{-4}` velocity tail at `t > 0`.  A theorem stated for
Schwartz solutions would say nothing at all about the evolution of a pure swirl, which is
the one thing §4 is for.  This is not pedantry: it is the difference between a hypothesis
that holds and one that is empty.

> **Required lemma LG-1 (propagation of the decay class).**  There is a class of initial
> data — smooth, divergence free, `⟨x⟩^{-4}`-decaying together with its first derivatives —
> and a `T > 0` such that the solution satisfies `(H3)` on `J = [0,T]` with a single
> constant `C_J`.
>
> *Unproved constant:* `C_J`.  *What supplies it:* Brandolese's localisation theory for
> the Navier–Stokes equations in `ℝ^n`, which establishes both that `|x|^{-(n+1)}` decay
> is propagated and that it is generically optimal.  *Status here:* assumed, not proved,
> and not proved anywhere in this repository.  **Everything below is conditional on it.**
> It is the only genuinely unproved input to §§2–7.

### 1.4 The pressure gauge must be fixed, and `(H4)` fixes it

Taking the divergence of the momentum equation gives `-Δp = ∂_i u_j ∂_j u_i`, which
determines `p` only up to an additive function harmonic in `x` (with arbitrary time
dependence).  For a velocity known only to lie in `L³` this ambiguity is real and not
removable by fiat — it is the elementary half of Serrin's observation that the pressure
of a distributional solution is not a function of the velocity alone — and it is fatal to
(1.1): replacing `p` by `p + h(x,t)` with `h` harmonic changes the right-hand side by
`3∫ h ∇·(|u|u) dx`, which is nonzero for perfectly ordinary `h` and `u`.  **The identity
is a statement about a gauge, not about `u` alone.**

`(H4)` removes the ambiguity.  With `u(t) ∈ L⁶` (which `(H1)` gives), Calderón–Zygmund
gives `p(t) = R_iR_j(u_iu_j) ∈ L³` with `‖p‖_{L³} ≤ C_CZ ‖u‖_{L⁶}²`, and an `L³` harmonic
function on `ℝ³` is zero by the mean-value property.  So `p` is the *unique* decaying
solution, and it is what the note means by "the pressure" everywhere below.

> **LG-2 (gauge uniqueness) — *proved*.**  Under `(H1)`, the equation `-Δp = ∂_iu_j∂_ju_i`
> has exactly one solution in `L³(ℝ³)` for each `t`, namely `R_iR_j(u_iu_j)`.  Two such
> solutions differ by a harmonic function in `L³`, hence by zero.

---

## 2. The derivation, step by step

Throughout §2 fix `t ∈ J` and write `u = u(·,t)`.  For `ε > 0` set

```
s_ε = ( |u|² + ε² )^{1/2} ,        F_ε(t) = ∫ ( s_ε³ - ε³ ) dx .
```

### 2.1 Why the `- ε³` is necessary, and why the regularisation is needed at all

**The subtraction is not cosmetic.**  `s_ε → ε` as `|x| → ∞`, so `s_ε³` does not decay and
`∫ s_ε³ dx = +∞` for every `ε > 0`.  With the constant removed,

```
0 ≤ s_ε³ - ε³ = (s_ε - ε)(s_ε² + s_ε ε + ε²) ≤ |u| ( |u|² + 3|u|ε + 3ε² )
              = |u|³ + 3ε|u|² + 3ε²|u| ,
```

using `s_ε ≤ |u| + ε` twice.  So `F_ε` is finite as soon as `u ∈ L¹ ∩ L² ∩ L³`, and the
`L¹` requirement is genuinely there — it is carried by the `3ε²|u|` term and by nothing
else.  `F` itself needs only `u ∈ L³`.  This is the first of the two places where `(H3)`
does real work rather than decorative work.

Letting `ε ↓ 0` in the same bound, with `ε ≤ 1`, gives the dominating function
`|u|³ + 3|u|² + 3|u| ∈ L¹`, and `s_ε³ - ε³ → |u|³` pointwise, so

```
F_ε(t) → F(t)   for every t ∈ J,   by dominated convergence.               (2.1)
```

**The regularisation is not needed for the time derivative.**  The map `w ↦ |w|³` on `ℝ³`
is `C²`, with gradient `3|w|w` and Hessian `3(|w|I + w⊗w/|w|)` whose norm is at most
`6|w| → 0` at the origin.  There is no non-differentiability to repair in `t`.  What is
not smooth is the *spatial* object `∂_k(|u|u_j)`, which involves `∂_k|u| = (u·∂_ku)/|u|`
and is merely bounded and measurable — `|u|` is Lipschitz but not `C¹` across the zero set
of `u`.  The integration by parts of §2.5 is the only step that needs `ε`, and everything
else is carried along for the ride.

### 2.2 Differentiating under the integral — *proved* (LG-3)

`∂_t(s_ε³) = 3s_ε² ∂_t s_ε = 3s_ε² (u·∂_t u)/s_ε = 3 s_ε (u·∂_t u)`, and this holds for
`ε = 0` as well with the value `3|u|(u·∂_tu)`.

> **LG-3 — *proved*.**  For every `ε ≥ 0`, `F_ε ∈ C¹(J)` and
> `F_ε'(t) = 3∫ s_ε (u·∂_t u) dx`.

*Proof.*  Write `Δ_h = (u(t+h) - u(t))/h`.  By the fundamental theorem of calculus applied
to `σ ↦ (|u(t)+σh Δ_h|² + ε²)^{3/2}`,

```
( F_ε(t+h) - F_ε(t) ) / h = 3∫∫_0^1 s_ε(u_σ) ( u_σ · Δ_h ) dσ dx ,
     u_σ = u(t) + σ h Δ_h .
```

`(H1)` gives `Δ_h → ∂_t u` in `L²`, and `|s_ε(w)w| ≤ |w|² + ε|w|`, so `s_ε(u_σ)u_σ` is
bounded in `L²` uniformly in `σ` and `h` small (using `u ∈ C(J; H³) ⊂ C(J; L² ∩ L^∞)`) and
converges to `s_ε(u)u` in `L²` as `h → 0`.  The pairing therefore converges to
`3∫ s_ε(u)(u·∂_tu)`.  Continuity of `t ↦ F_ε'(t)` follows the same way from
`u ∈ C¹(J;L²)`. ∎

Substituting the equation `∂_t u = νΔu - (u·∇)u - ∇p`, which holds in `L²` at each `t`,
splits `F_ε'` into three pieces:

```
F_ε' = T_ε + P_ε + V_ε ,
T_ε = -3∫ s_ε u·(u·∇)u dx ,   P_ε = -3∫ s_ε (u·∇p) dx ,   V_ε = 3ν∫ s_ε (u·Δu) dx .
```

### 2.3 The transport term is exactly zero — *proved* (LG-4)

No limit is taken here.  `T_ε = 0` for every `ε > 0`, and for `ε = 0`.

Since `u·(u·∇)u = ½ u·∇|u|² = ½ u·∇(s_ε²)` (the constant `ε²` differentiates away) and
`∇(s_ε³) = (3/2) s_ε ∇(s_ε²)`,

```
s_ε u·(u·∇)u = ½ s_ε u·∇(s_ε²) = ⅓ u·∇(s_ε³) = ⅓ u·∇(s_ε³ - ε³) = ⅓ ∇·( u (s_ε³ - ε³) ) ,
```

the last step by `∇·u = 0`.  At `ε = 0` this is the familiar
`|u| u·(u·∇)u = ⅓ u·∇(|u|³)`.  Integrating over `B_R` and applying the divergence theorem,

```
|T_ε| = | ∮_{|x|=R} (s_ε³ - ε³)(u·n) dS |
      ≤ 4πR² · ( C_J³R^{-12} + 3εC_J²R^{-8} + 3ε²C_JR^{-4} ) · C_JR^{-4}
      = O(ε² R^{-6}) → 0 .
```

> **LG-4 — *proved*.**  For every divergence-free field in the class `(H3)` and every
> `ε ≥ 0`, `T_ε = 0`.

This is the structural fact that §3 turns into a statement about the search: **advection
cannot change the critical norm at all**, at any time, for any divergence-free field.  It
is not an estimate with a small constant; it is an exact cancellation, and it holds
before any regularisation or limit.

### 2.4 The pressure term — *proved* (LG-5)

Integrate by parts on `B_R`, using `∇·u = 0` so that `∇·(s_ε u) = u·∇s_ε`:

```
P_ε = -3∫_{B_R} s_ε (u·∇p) dx
    = 3∫_{B_R} p ∇·(s_ε u) dx - 3∮_{|x|=R} p s_ε (u·n) dS
    = 3∫_{B_R} p (u·∇s_ε) dx - 3∮_{|x|=R} p s_ε (u·n) dS .
```

**The trap and how it is defused.**  `s_ε` does not decay: `s_ε → ε ≠ 0` at infinity, so
the surface integrand is *not* small because of `s_ε`.  It is small because of the
explicit factor `u`, which decays like `R^{-4}`.  With the sharp pressure tail
`|p| ≤ C R^{-3}` of §1.3,

```
| ∮_{|x|=R} p s_ε (u·n) dS |  ≤  4πR² · C R^{-3} · (ε + C_J R^{-4}) · C_J R^{-4}
                              =  O(ε R^{-5}) + O(R^{-9})  →  0 .           (2.2)
```

The estimate is stated in this form because it is sharp and because it exhibits exactly
which factor does the work.  It is worth recording that the conclusion does not depend on
having the sharp tail at all.  Since `p(t) ∈ L³` by LG-2, the radial function
`h(R) = ∮_{|x|=R}|p|³dS` is in `L¹(0,∞)`, so `liminf_{R→∞} R h(R) = 0`; along a sequence
`R_n → ∞` realising it, Hölder gives
`∮_{|x|=R_n}|p| dS ≤ (4πR_n²)^{2/3} h(R_n)^{1/3} = o(R_n)`, and the surface term is
bounded by `sup_{|x|=R_n}(s_ε|u|) · o(R_n) = O(ε R_n^{-3}) → 0`.  **So the boundary term
vanishes using only `p ∈ L³` and `(H3)`, along a subsequence of radii, which is all the
divergence theorem needs.**  No pointwise pressure-decay lemma is required, and none is
listed in §9.

> **LG-5 — *proved*.**  For every `ε > 0`, `P_ε = 3∫ p (u·∇s_ε) dx`, the integral being
> absolutely convergent.

### 2.5 The viscous term — *proved* (LG-6)

This is the step that needs `ε`.  Integrating by parts on `B_R`,

```
V_ε = 3ν∫_{B_R} s_ε u_j Δu_j dx
    = -3ν∫_{B_R} ∂_k(s_ε u_j) ∂_k u_j dx + 3ν∮_{|x|=R} s_ε u_j ∂_ku_j n_k dS .
```

`s_ε` is smooth in `u` for `ε > 0`, with `∂_k s_ε = (u·∂_k u)/s_ε`, so

```
∂_k(s_ε u_j) ∂_k u_j = (u·∂_ku)(u_j∂_ku_j)/s_ε + s_ε ∂_ku_j∂_ku_j
                     = (u·∂_ku)²/s_ε + s_ε |∇u|² .
```

At `ε = 0` the first factor `∂_k|u| = (u·∂_ku)/|u|` is undefined on `{u = 0}` and `|u|u`
is only Lipschitz, so this manipulation is not available; that, and nothing else, is what
the regularisation buys.  The boundary term is

```
| 3ν ∮_{|x|=R} s_ε u_j ∂_ku_j n_k dS | ≤ 3ν · 4πR² · (ε + C_JR^{-4}) · C_JR^{-4} · C_JR^{-4}
                                       = O(ε R^{-6}) → 0 .
```

> **LG-6 — *proved*.**  For every `ε > 0`,
> ```
> V_ε = -3ν ∫ ( s_ε |∇u|² + Σ_k (u·∂_ku)²/s_ε ) dx ,
> ```
> and **both integrands are pointwise nonnegative**, so `V_ε ≤ 0`.

### 2.6 The `ε → 0` limit, term by term — *proved* (LG-7)

Each integrand converges pointwise **everywhere**, not merely almost everywhere.  On
`{u ≠ 0}` this is continuity of `ε ↦ s_ε` at `0`.  On `{u = 0}` the numerator of every
quotient vanishes identically, so for every `ε > 0` the value is already `0`, and the
limit is `0`, which is the conventional value assigned to `(u·∂_ku)²/|u|` there.  **The
zero set of `u` need not be Lebesgue-null and the argument does not need it to be**: it is
handled by the numerator, not by measure theory.  Likewise `u·∇s_ε = u_k(u·∂_ku)/s_ε`
vanishes identically on `{u = 0}` for every `ε`.

The three dominating functions, valid for `0 < ε ≤ 1`, using `|u| ≤ s_ε ≤ |u| + ε`:

| term | integrand | dominating function | why integrable |
|---|---|---|---|
| pressure | `3 p (u·∇s_ε)` | `3\|p\|\|u\|\|∇u\|` | `‖p‖_{L³}‖ \|u\|\|∇u\| ‖_{L^{3/2}} < ∞`, and `‖p‖_{L³} ≤ C_CZ‖u‖_{L⁶}²` |
| viscous-A | `3ν Σ_k (u·∂_ku)²/s_ε` | `3ν \|u\| \|∇u\|²` | `≤ 3νC_J³⟨x⟩^{-12}` |
| viscous-B | `3ν s_ε \|∇u\|²` | `3ν (\|u\|+1) \|∇u\|²` | `≤ 3ν(C_J+1)C_J²⟨x⟩^{-8}`; also finite from `∇u ∈ L²` and `u ∈ L^∞` |

The pressure bound uses `|u_k(u·∂_ku)/s_ε| ≤ |u|²|∇u|/s_ε ≤ |u||∇u|` since `s_ε ≥ |u|`;
viscous-A uses `Σ_k(u·∂_ku)² ≤ |u|²|∇u|²` and again `s_ε ≥ |u|`.  Both `⟨x⟩^{-12}` and
`⟨x⟩^{-8}` are integrable on `ℝ³`.  The Calderón–Zygmund bound is quoted because it is
what one would use outside the decay class; inside the class the integrability is
elementary.

**Why monotone convergence does not apply.**  `s_ε` is strictly increasing in `ε`.  So as
`ε ↓ 0` the viscous-B integrand `s_ε|∇u|²` **decreases** to `|u||∇u|²` while the viscous-A
integrand `(u·∂_ku)²/s_ε` **increases** to `(u·∂_ku)²/|u|`.  One piece goes up, the other
goes down; their sum is not monotone in `ε`, and neither is the pressure integrand in
either direction.  Each viscous piece separately could be treated by monotone convergence
— the decreasing one only after producing an integrable majorant at `ε = 1`, which is
dominated convergence in disguise — but the identity is a statement about the sum and the
sum has no monotonicity to exploit.  Dominated convergence handles all three uniformly
and is in any case what §2.7 needs.

Finally `(u·∂_ku)²/|u| = (|u|∂_k|u|)²/|u| = |u| (∂_k|u|)²`, so viscous-A assembles to
`|u||∇|u||²`, giving the second integrand of (1.1).

> **LG-7 — *proved*.**  `P_ε → 3∫p(u·∇|u|)dx = 3∫ p ∇·(|u|u) dx` and
> `V_ε → -3ν∫(|u||∇u|² + |u||∇|u||²)dx` as `ε ↓ 0`, for each fixed `t ∈ J`.

### 2.7 The step that needs care: `lim F_ε' ≠ (lim F_ε)'` in general — *proved* (LG-8)

At this point one has `F_ε → F` pointwise and `F_ε' → G` pointwise, where `G` is the
right-hand side of (1.1).  **That is not enough**, and the gap is not a technicality.
Take

```
F_n(t) = t / (1 + n t²) .
```

Then `sup_t |F_n(t)| = 1/(2√n) → 0`, so `F_n → 0` *uniformly*; yet `F_n'(0) = 1` for every
`n`, so `F_n'(0) → 1 ≠ 0 = (lim F_n)'(0)`.  Uniform convergence of the functions is
irrelevant.  What is needed is uniform convergence of the **derivatives**, and that is what
the following explicit rate supplies.

> **LG-8 (the `O(ε)` rate) — *proved*.**  For every `t ∈ J` and every `0 < ε ≤ 1`,
> ```
> | F_ε'(t) - G(t) | ≤ 3ε ( ‖p(t)‖_{L²} ‖∇u(t)‖_{L²} + 2ν ‖∇u(t)‖_{L²}² ) .     (2.3)
> ```

*Proof.*  Write `w_k = u·∂_k u`, so `|w_k| ≤ |u||∇u|`, and note
`| |u| - s_ε | = ε²/(s_ε + |u|) ≤ ε²/s_ε ≤ ε` since `s_ε ≥ ε`.

*Pressure.*  On `{u = 0}` both integrands vanish identically and contribute nothing.  On
`{u ≠ 0}`,

```
| u·∇s_ε - u·∇|u| | = | u_k w_k ( |u| - s_ε ) / ( s_ε |u| ) |
                    ≤ |u| · |u||∇u| · ε / ( s_ε |u| ) = ε |u||∇u|/s_ε ≤ ε |∇u| ,
```

using `|u| ≤ s_ε`.  Hence the pressure difference is at most
`3ε∫|p||∇u| ≤ 3ε‖p‖_{L²}‖∇u‖_{L²}`.

*Viscous-B.*  `| s_ε - |u| | ≤ ε` gives `3ν ε ‖∇u‖_{L²}²`.

*Viscous-A.*  On `{u ≠ 0}`,
`Σ_k w_k²( |u| - s_ε )/( s_ε|u| ) ≤ |u|²|∇u|² ε/(s_ε|u|) = ε|u||∇u|²/s_ε ≤ ε|∇u|²`,
giving another `3νε‖∇u‖_{L²}²`; on `{u = 0}` both sides are `0`. ∎

The bracket in (2.3) is bounded **uniformly on the compact interval `J`**: `u ∈ C(J;H³)`
makes `sup_J‖∇u‖_{L²}` finite, and `‖p‖_{L²} ≤ C_CZ ‖u‖_{L⁴}²` with `H³ ↪ L⁴` makes
`sup_J‖p‖_{L²}` finite.  Therefore `F_ε' → G` **uniformly on `J`**.  Combined with
`F_ε → F` pointwise (2.1), the classical theorem on differentiating a limit gives: `F` is
differentiable on `J` with `F' = G`.  Each `F_ε'` is continuous by LG-3, and a uniform
limit of continuous functions is continuous, so `G` is continuous and `F ∈ C¹(J)`.  This
completes the proof of (1.1). ∎

---

## 3. The three structural consequences

Read off (1.1) directly.

**Advection contributes exactly nothing.**  By LG-4 the transport term is identically
zero for *every* divergence-free field, not approximately zero and not zero on average.
So the mechanism that everyone pictures when they picture blow-up — nonlinear advective
steepening — is by itself incapable of increasing the critical norm at any instant.
Whatever advection does to the critical norm, it does it indirectly, by reshaping the
field so that the pressure term becomes favourable later.

**The viscous term is non-positive, always.**  Both integrands in `V` are products of
`|u| ≥ 0` with a square.  Moreover Kato's inequality `|∇|u|| ≤ |∇u|` pointwise a.e. pins
`V` between two multiples of the same quantity,

```
-6ν∫|u||∇u|² dx  ≤  V  ≤  -3ν∫|u||∇u|² dx ,
```

so `V` is comparable to `-3ν∫|u||∇u|²` up to a factor of two and can never be made small
by arranging cancellation between its two pieces.

**Only the pressure channel can be positive.**  `F'(0) > 0` requires `P > 0`, and `P` is
the *only* candidate.  Its integrand, in the form (1.2), is `-3|u|(u·∇p)`, and in the
form (1.1) it is `3p ∇·(|u|u) = 3p (u·∇|u|)`.

The geometric reading of the second form is worth stating because it is what makes §4
obvious in hindsight.  `∇·(|u|u) = u·∇|u|` is the rate at which a fluid particle carries
its own speed along.  **The critical norm can only grow if the flow advects its own speed
in correlation with the pressure**: `u` must have a component along `∇|u|`, and that
component must be aligned with regions of the right pressure sign.  A field whose velocity
is everywhere perpendicular to the gradient of its own speed has `P = 0` exactly, whatever
its amplitude, whatever its Reynolds number.

---

## 4. The pure-swirl no-go

### 4.1 Statement

> **Corollary (pure-swirl no-go).**  Let
> ```
> u_0 = u^θ(r,z) e_θ ,        u^θ(r,z) = r g(r², z) ,
> ```
> with `g` smooth and `u_0` in the class `(H3)`.  Then
> ```
> J(u_0) := F'(0) = -3ν∫ ( |u_0||∇u_0|² + |u_0||∇|u_0||² ) dx ≤ 0 ,
> ```
> with equality if and only if `u_0 ≡ 0`.

**The axis condition is a hypothesis, not decoration.**  The Cartesian field is
`u_0 = g(x²+y², z)·(-y, x, 0)`.  Written as `u^θ = h(r,z)` with a general `h`, the field is
discontinuous on the axis unless `h(0,z) = 0`, and it fails to be `C¹` unless `h` is `r`
times a smooth function of `(r², z)`.  Without `u^θ = r g(r²,z)` the datum is not in
`(H1)`, the identity does not apply, and there is nothing to prove or disprove.

### 4.2 The proof, in the strong form

The weak form of the argument computes `∇·(|u_0|u_0)` and finds it zero:
`|u_0|u_0 = |u^θ|u^θ e_θ =: G(r,z) e_θ`, whose cylindrical divergence is
`r^{-1}∂_θ G = 0`.  That is correct, but it is the wrong proof, because it computes the
integrand of (1.1) — the form that came out of an integration by parts and an `ε → 0`
limit.  The following argument works on (1.2) directly and touches neither.

> **LG-9 (rotational equivariance of the Riesz pressure) — *proved*.**  Let `R` be a
> rotation and let `(R·u)(x) = R u(R^{-1}x)` be the induced action on vector fields.  Then
> the decaying pressure satisfies `p[R·u] = p[u] ∘ R^{-1}`.

*Proof.*  The source `-Δp = ∂_iu_j∂_ju_i = tr((∇u)²)` is a scalar invariant of `∇u`, and
`∇(R·u)(x) = R (∇u)(R^{-1}x) R^{-1}`, so `tr((∇(R·u))²)(x) = tr((∇u)²)(R^{-1}x)`: the
source transforms as a scalar.  `Δ` commutes with rotations, so `p[u]∘R^{-1}` solves the
transformed Poisson equation and lies in `L³`; by LG-2 it is *the* solution. ∎

Now let `R_φ` be the rotation by `φ` about `e_z`.  A pure swirl is axisymmetric in exactly
the sense `R_φ·u_0 = u_0` for every `φ`.  By LG-9, `p_0 ∘ R_φ^{-1} = p_0`, i.e.

```
∂_θ p_0 ≡ 0 .
```

Therefore, pointwise and everywhere,

```
u_0 · ∇p_0 = u^r ∂_r p_0 + u^θ r^{-1} ∂_θ p_0 + u^z ∂_z p_0
           =   0        +        0           +    0        = 0 ,
```

the first and third terms because `u^r = u^z = 0` and the middle because `p_0` is
axisymmetric.  The pressure integrand of (1.2) is `-3|u_0|(u_0·∇p_0) ≡ 0` **as a
function**, so `P = 0`.  This is a pointwise statement about the datum; it does not use
the integration by parts of §2.4, it does not use the dominated convergence of §2.6, and
it does not use the `ε → 0` interchange of §2.7.  The corollary is therefore as robust as
the trivially checkable facts that `u^r = u^z = 0` and that `∂_θ p_0 = 0`.

The same conclusion in the language of §3: a pure swirl has `|u_0| = |u^θ(r,z)|`, a
function of `(r,z)` only, so `∇|u_0|` lies in the meridional plane while `u_0` is purely
azimuthal.  **A pure swirl is exactly perpendicular to the gradient of its own speed, at
every point.**  It is the exact zero of the only channel that can ever be positive.

Hence `J(u_0) = V ≤ 0` by LG-6.

### 4.3 The equality case, done properly

The one-line version — "`∫|u||∇u|² = 0` gives `∇u = 0` gives `u` constant gives `u = 0`" —
is a **non-sequitur**, and it is worth saying why so that nobody re-derives it.  The
integrand `|u||∇u|²` vanishes at every point of `{u = 0}` *regardless of the value of
`∇u` there*.  Vanishing of the integral gives no information about `∇u` on the zero set,
and the zero set of a nonzero decaying field is large — it contains the whole complement
of any compact support.  A correct argument has to be topological.

> **LG-10 (equality case) — *proved*.**  If `u_0` is continuous, `u_0(x) → 0` as
> `|x| → ∞`, and `∫|u_0||∇u_0|² dx = 0`, then `u_0 ≡ 0`.

*Proof.*  The integrand is continuous and nonnegative, so `|u_0||∇u_0|² ≡ 0`.  Let
`Ω = {x : u_0(x) ≠ 0}`, an open set on which `∇u_0 = 0`; so `u_0` is locally constant on
`Ω`.  Let `C` be a connected component of `Ω`, say `u_0 ≡ c` on `C` with `c ≠ 0`.  `C` is
open.  `C` is also **closed** in `ℝ³`: if `x ∈ ∂C` then `u_0(x) = c ≠ 0` by continuity, so
`x ∈ Ω`, and continuity gives a ball `B ∋ x` with `B ⊂ Ω`; `B` is connected and meets `C`,
so `B ∪ C` is a connected subset of `Ω` containing `C`, whence `B ⊂ C` by maximality of
the component, and `x ∈ C`.  A nonempty subset of `ℝ³` that is both open and closed is all
of `ℝ³`, since `ℝ³` is connected.  Then `u_0 ≡ c ≠ 0` everywhere, contradicting
`u_0(x) → 0`.  So `Ω` has no components, i.e. `Ω = ∅`. ∎

Under `(H3)` the decay hypothesis holds, so `J(u_0) = 0` forces `u_0 ≡ 0`, and for every
nonzero pure swirl `J(u_0) < 0` **strictly**.

---

## 5. What the no-go does not exclude, and the correct framing

**It is a statement at one instant.**  `J(u_0) = F'(0)`.  It says the critical norm is
strictly decreasing at `t = 0`.  It says nothing whatever about `F(t)` for `t > 0`.

**Pure swirl is not preserved by the flow.**  The axisymmetric vorticity `ω^θ` obeys

```
∂_t ω^θ + (u·∇)ω^θ - (u^r/r) ω^θ = ν( Δ - r^{-2} ) ω^θ + r^{-1} ∂_z( (u^θ)² ) ,
```

and the source term `r^{-1}∂_z((u^θ)²)` does not vanish for a generic swirl.  Starting from
`ω^θ(0) = 0` — which is what "pure swirl" means — the equation gives
`ω^θ = t r^{-1}∂_z((u^θ)²) + O(t²)`, so meridional flow appears **immediately** at first
order in `t`, and with it a nonzero `∇·(|u|u)`.  The only way to avoid it is
`∂_z u^θ ≡ 0`, the columnar case; and within the hypothesis class `(H3)` that case is
vacuous, because a `z`-independent `u^θ` that tends to zero as `|z| → ∞` is identically
zero.  (The genuinely columnar swirl `u^θ = u^θ(r)` is a global smooth solution of
Navier–Stokes for trivial reasons — it obeys a linear heat equation and the pressure
merely balances the centrifugal force — but it has infinite energy and is not in the
class.)

**The pressure channel switches on at `t = 0+`, and its sign is not determined.**  For
`t > 0` the field is no longer a pure swirl, `P(t) ≠ 0` in general, and nothing in this
note fixes whether `P(t) > |V(t)|`.

**So what is this result?**  It is a **kill**, not a support.  It removes a search lane; it
does not provide evidence for any conjecture about singularities.  Concretely, the family
of purely azimuthal initial data has `P ≡ 0` identically, which by §7 means its critical
Reynolds number is `+∞`: **no amount of increasing `Re` can make `J` positive for a pure
swirl**, because the quantity that `Re` multiplies is zero.

This explains the outcome of Gate 7 completely, in the following exact sense.  All three
preregistered whole-space families S, A and H are pure swirls (`u^r = u^z = 0` in the
datum, meridional flow generated only through the elliptic solve at order `t`).  The
eighteen runs swept `Re ∈ {10, 25, 50, 100, 200, 400}` looking for growth of the critical
norm.  By the corollary the initial rate was negative in every one of them **before any
computation was performed**, and by the scaling analysis of §7 the sweep in `Re` was
searching along a direction that could not possibly change that sign.  The uniform decay
recorded in `outputs/tau_continuation_gate7_v1/` was not an experimental finding about
Navier–Stokes; it was a consequence of the choice of data class.  Leaving that class
requires meridional flow **in the datum**, which is what
`ns_certificate_lab.mixed_initial_data` exists to supply.

---

## 6. The axial parity selection rule

Work in the generator variables used throughout the repository:

```
u^θ = r u₁(r²,z) ,     u^r = -r ∂_z ψ₁(r²,z) ,     u^z = 2ψ₁ + r ∂_r ψ₁ .
```

The meridional part is divergence free for any `ψ₁`, and the swirl part is `θ`-independent,
so the sum is exactly divergence free with no projection.  Let `σ_u, σ_ψ ∈ {+1,-1}` be the
parities of `u₁` and `ψ₁` under `z ↦ -z` (writing `f(r,-z) = σ f(r,z)`).

> **LG-11 (parity selection rule) — *proved*.**  On a `z`-symmetric domain, the pressure
> contribution `P` vanishes identically unless `σ_ψ = -1`, that is, unless the **stream
> generator is odd in `z`**.  The swirl parity `σ_u` does not enter.

*Proof.*  The component parities are `σ(u^θ) = σ_u`, `σ(u^r) = -σ_ψ` (one `∂_z` applied to
`ψ₁`), `σ(u^z) = σ_ψ`.

`|u|² = (u^r)² + (u^θ)² + (u^z)²` is a sum of squares, hence **always even**, hence `|u|`
is even, whatever the generator parities.

The pressure source `∂_iu_j∂_ju_i = tr((∇u)²)` in the orthonormal cylindrical frame is

```
tr((∇u)²) = G_rr² + G_θθ² + G_zz² + 2 G_rθ G_θr + 2 G_rz G_zr + 2 G_θz G_zθ ,
```

with

| entry | value | parity |
|---|---|---|
| `G_rr` | `∂_r u^r` | `-σ_ψ` |
| `G_θθ` | `u^r/r` | `-σ_ψ` |
| `G_zz` | `∂_z u^z` | `-σ_ψ` |
| `G_rθ` | `∂_r u^θ` | `σ_u` |
| `G_θr` | `-u^θ/r` | `σ_u` |
| `G_rz` | `∂_r u^z` | `σ_ψ` |
| `G_zr` | `∂_z u^r` | `σ_ψ` |
| `G_θz` | `0` | — |
| `G_zθ` | `∂_z u^θ` | `-σ_u` |

Every term of the trace is either a square or one of the two cross products
`G_rθG_θr` and `G_rzG_zr`, each a product of two factors of **equal** parity.  So the
pressure source is **always even**, and by uniqueness of the decaying solution (LG-2) on a
`z`-symmetric configuration, `p` is even.

The pressure integrand of (1.2) is `-3|u|(u^r∂_rp + u^z∂_zp)`.  `|u|` is even; `∂_rp` is
even and `∂_zp` is odd; so `u^r∂_rp` has parity `-σ_ψ` and `u^z∂_zp` has parity
`σ_ψ·(-1) = -σ_ψ`.  The whole integrand therefore has parity `-σ_ψ`.  If `σ_ψ = +1` the
integrand is odd in `z` and integrates to zero on a `z`-symmetric domain. ∎

Two remarks.  The rule is one-directional: `σ_ψ = -1` does not make `P` positive, it only
removes the forced cancellation.  And the viscous integrand `|u|(|∇u|² + |∇|u||²)` is a
product of even factors and is therefore **always even**, so it never cancels — parity can
only kill the helpful term, never the harmful one.

An odd `ψ₁` is also the physically sensible choice, independently of the parity count: it
gives `u^r` even and `u^z` odd, the axisymmetric strain flow that converges radially at
every height and ejects along both `±z`.  An even `ψ₁` gives inflow on one side and
outflow on the other, and the two halves cancel — which is exactly what the parity
computation says.

---

## 7. Scaling, and the search objective

### 7.1 The two-parameter family

Fix a shape `û` and set `u_{A,L}(x) = A û(x/L)`.  Then `∇u ~ (A/L)`, `p ~ A²`,
`dx ~ L³`, and

```
P  ~  A⁴ L² ,        V  ~  -ν A³ L ,        ‖u‖_{L³}³ ~ A³L³ ,
```

so

```
P / |V|  ~  A L / ν  =:  Re .                                          (7.1)
```

> **LG-12 (the sign criterion) — *proved*.**  Along the two-parameter family, `J > 0` if
> and only if `Re > Re_crit(shape)`, where
> ```
> Re_crit  =  ( A L / ν ) · |V| / P                                     (7.2)
> ```
> evaluated at any one representative of the family; the value is independent of which
> representative is used, so it is a function of the shape alone.

There is one further exact symmetry, which removes a search direction for free.  Under
`u ↦ -u` the pressure `p = R_iR_j(u_iu_j)` is unchanged (it is quadratic) while
`∇·(|u|u)` changes sign; the viscous integrands are even.  So `P` is **odd** and `V` is
**even** under the sign flip, and for any shape with `P ≠ 0` exactly one of `±u` has
`P > 0`.  The overall sign is read off, never searched.

### 7.2 Consequence (a): the critical rescaling cannot manufacture growth

Take the Navier–Stokes scaling that preserves the critical norm:

```
A ↦ λA ,     L ↦ L/λ .
```

Then `‖u‖_{L³}³ ~ A³L³` is invariant, as it must be; and

```
P ~ A⁴L² ↦ λ⁴ · λ^{-2} · (A⁴L²) = λ² P ,
V ~ νA³L ↦ λ³ · λ^{-1} · (νA³L) = λ² V .
```

**Both scale by the same factor `λ²`, so the ratio `P/|V|` is exactly invariant**, and so
is `Re = AL/ν`.  A shape that fails `J > 0` cannot be rescued by zooming in and amplifying,
however violently.  Only the shape moves `Re_crit`.  This is not a heuristic: it is the
statement that `Re_crit` is a scale-invariant functional, and it is what makes the search
in §7.3 a search over an eight-dimensional sphere of amplitudes rather than over a
noncompact family.

### 7.3 The objective, and the amplitude-variable trap

The search objective is

```
minimise   Re_crit(shape)  =  A L · ( |V| / ν ) / P                     (7.3)
```

over a fixed compactly supported generator basis, with `A` the **physical** velocity
amplitude `max|u|` and `L` the reference length of the basis.  `Re_crit` does not mention
`ν` (the `ν` in `|V|/ν` cancels the one in the definition), is invariant under
`u ↦ λu`, and is `+∞` exactly when `P ≤ 0`.

**One power of `L`, not two, because `A` is the physical amplitude.**  This is the same
Reynolds number the repository already uses.  `nondimensional.py` and the Gate 7 note write
`Re = A_{u₁} L² / ν`, where `A_{u₁}` is the amplitude of the *swirl variable* `u₁ = u^θ/r`.
Since `u^θ = r u₁`, the physical amplitude is `A_u = A_{u₁} L`, and

```
A_u L / ν  =  ( A_{u₁} L ) L / ν  =  A_{u₁} L² / ν .
```

Same number, different variable.  The trap is to pair one convention's amplitude with the
other convention's power of `L`; the two must be changed together, and a formula containing
`max|u|` together with `L²` is off by a factor of `L`.  Any implementation of (7.3) should
be checked against this identity before its output is compared with a Gate 7 Reynolds
number.

---

## 8. Limitations, stated bluntly

**`J > 0` is a statement at one instant, and the endpoint theorem is a statement about a
limit.**  What Escauriaza–Seregin–Šverák requires, to conclude that `T` is a singular
time, is `limsup_{t→T} ‖u(t)‖_{L³} = ∞`.  What `J(u_0) > 0` gives is
`d/dt‖u(t)‖_{L³}³ > 0` at the single instant `t = 0`.  **Nothing connects the two.**  There
is no monotonicity argument, no differential inequality of the form `F' ≥ cF^α` with
`α > 1`, and no reason from (1.1) to expect one: the pressure channel is a correlation
between two objects that both evolve, and correlations do not obey autonomous ODEs.
Finding a datum with `J > 0` would be a necessary-looking first step in a programme whose
second step does not exist.

**The pressure channel is sign-indefinite for `t > 0`.**  `P(t) = 3∫p(u·∇|u|)dx` has no
definite sign and no known lower bound in terms of `F` alone.  A datum can start with
`P > 0` and have `P < 0` an instant later; the viscous term meanwhile is negative
throughout and, by §3, is bounded below by `-6ν∫|u||∇u|²`, which grows as the field
steepens.  Nothing rules out — and nothing establishes — that `P` loses this race.

**The regime where `P` wins is the regime the numerics are worst at.**  By (7.1), `J > 0`
requires `Re > Re_crit`.  Large `Re` means thin structures, large gradients, and a small
viscous length; that is exactly where the fixed grid, the truncated box and the finite
integrator are least trustworthy, and where the `H1`/`H2` slab hypotheses of
`slab_certificate.py` and the discrete-to-continuous gap `HS-5` of
[`hs_error_propagation.md`](hs_error_propagation.md) bite hardest.  A positive `J` computed
at large `Re` on a coarse grid is a statement about the grid.

**Even `J > 0` on a whole interval gives bounded growth.**  If `F' > 0` on `[0,T]` then
`F(T) > F(0)`.  That is a finite increase.  Divergence needs a mechanism that survives
`T → T_*`, and a positive derivative at any finite collection of times does not supply one.

**Pre-registered kill condition.**  Recorded here, before the optimiser of
`ns_certificate_lab.l3_optimizer` is run to convergence on the fixed basis:

> If the best `Re_crit` attained over the preregistered eight-dimensional search basis does
> not fall **below `10²`**, this lane is dead and will be closed with a negative result.

Two things make this a real bar rather than a movable one.  First, §7.2 forbids the obvious
rescue: `Re_crit` is scale invariant, so no amount of rescaling, zooming or amplifying can
lower it, and enlarging the amplitude search range is provably futile.  Only enlarging the
*basis* — a different, and declarable, experiment — can change the number.  Second, there
is a reference point to calibrate against: a generic band-limited field gives
`|P/V| ≈ 1.07e-3`, that is, `Re_crit ≈ 9e2`.  Beating that by an order of magnitude is the
minimum that would distinguish a designed shape from an arbitrary one.  A best value that
merely lands near `9e2` means the optimiser has found nothing that a random field does not
already do.

**Nothing here bears on the Clay problem.**  The identity is classical.  The corollary is
an exclusion of a class of initial data at a single instant.  The objective in §7 is a
numerical search target whose success condition is `J(u_0) > 0` for one datum at one
time — which, as the first paragraph of this section says, is not a step toward a
singularity proof, only a step out of a region where a singularity is provably not starting
to form.  No claim is made, and none should be read in, that any of this constitutes
progress on the existence and smoothness problem.

---

## 9. Summary of the steps

| id | statement | status |
|---|---|---|
| LG-1 | the decay class `(H3)` is propagated on a compact interval | **required lemma**; supplies `C_J`; Brandolese, not proved here |
| LG-2 | the decaying pressure is unique, so `(H4)` is a well-posed gauge | *proved* |
| LG-3 | `F_ε ∈ C¹` and `F_ε' = 3∫s_ε(u·∂_tu)`, for every `ε ≥ 0` | *proved* |
| LG-4 | the transport term is **exactly** zero, for every `ε ≥ 0` | *proved* |
| LG-5 | pressure integration by parts; boundary term vanishes | *proved* |
| LG-6 | viscous integration by parts; two nonnegative integrands | *proved* |
| LG-7 | `ε → 0` by dominated convergence, pointwise **everywhere** | *proved* |
| LG-8 | the `O(ε)` rate (2.3), uniform on `J`; hence `F ∈ C¹`, `F' = G` | *proved* |
| LG-9 | rotational equivariance of the Riesz pressure | *proved* |
| LG-10 | the equality case, by connectedness | *proved* |
| LG-11 | the axial parity selection rule | *proved* |
| LG-12 | the scaling criterion `J > 0 ⟺ Re > Re_crit`, and its invariance | *proved* |

Constants that are named but not computed anywhere in this repository:

| symbol | what it bounds | supplied by | status |
|---|---|---|---|
| `C_J` | `\|u\| + \|∇u\| ≤ C_J⟨x⟩^{-4}` on `J` | LG-1 | assumed |
| `C_CZ(q)` | `‖R_iR_jf‖_{L^q} ≤ C_CZ‖f‖_{L^q}`, `1<q<∞` | Calderón–Zygmund | classical, not proved here |
| `sup_J‖∇u‖_{L²}`, `sup_J‖p‖_{L²}` | the rate in (2.3) | `(H1)` plus `C_CZ` | finite, not quantified |
| `Re_crit(shape)` | the sign of `J` | numerics | measured, never certified |

## 10. References

* L. Escauriaza, G. A. Seregin, V. Šverák, *`L_{3,∞}`-solutions of the Navier–Stokes
  equations and backward uniqueness*, Russian Math. Surveys **58** (2003), 211–250.
* L. Brandolese, *On the localization of symmetric and asymmetric solutions of the
  Navier–Stokes equations in `ℝⁿ`*, C. R. Acad. Sci. Paris **332** (2001), 125–130; and
  *Space-time decay of Navier–Stokes flows invariant under rotations*, Math. Ann. **329**
  (2004), 685–706.
* J. Serrin, *On the interior regularity of weak solutions of the Navier–Stokes equations*,
  Arch. Rational Mech. Anal. **9** (1962), 187–195, for the non-uniqueness of the pressure
  of a distributional solution.
* Companion notes in this repository:
  [`hs_error_propagation.md`](hs_error_propagation.md) (the `H^s` error chain and the
  discrete-to-continuous gap `HS-5`),
  [`critical_l3_obstruction.md`](critical_l3_obstruction.md) (the endpoint obstruction and
  the rescaled-candidate no-go),
  [`tau_continuation_gate7.md`](tau_continuation_gate7.md) (the eighteen runs whose uniform
  decay §5 explains).
* Implementation: `src/ns_certificate_lab/l3_generation.py` (the identity and the no-go),
  `src/ns_certificate_lab/mixed_initial_data.py` (data with meridional flow, and the parity
  rule), `src/ns_certificate_lab/l3_optimizer.py` (the objective of §7),
  `src/ns_certificate_lab/l3_certificate.py` (the interval certificate for `J > 0`,
  conditional on the named pressure hypothesis `P1`).
