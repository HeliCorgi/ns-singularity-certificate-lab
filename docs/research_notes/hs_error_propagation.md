# Moving the error estimate from `L^∞` to `H^s`, `s > 5/2`

**Status: an incomplete derivation.**  Every step that is complete is marked
*proved*.  Every step that is not is marked *required lemma* and carries a name,
a precise statement and the constant it would supply.  Nothing below licenses
calling the existing certificate an unconditional PDE proof, and §7 says exactly
what the certificate does and does not currently establish.

## 1. Why `L^∞` is the wrong space

The certificate layer (`formal/NSSingularity/CertificateLayer.lean`, F-17–F-19)
propagates the error `e = u - ũ` in `L^∞` through the Grönwall step

```
d/dt ‖e‖_∞ ≤ C ‖e‖_∞ + ‖R‖_∞ .
```

That inequality is not available for Navier–Stokes, and the reason is structural
rather than technical.

1. **The pressure is nonlocal.**  Eliminating `p` gives `u_t + P(u·∇u) = νΔu`
   with `P` the Leray projector.  `P` is bounded on `L^q` for `1 < q < ∞` and on
   every `H^s`, but **not** on `L^∞`.  So the very first step of writing the
   equation for `e` leaves `L^∞`.
2. **The advection term loses a derivative.**  `(e·∇)ũ` needs `∇ũ`, and
   `‖∇ũ‖_∞` is not controlled by `‖ũ‖_∞`.  Any `L^∞` closure has to assume a
   gradient bound that the norm itself does not carry.
3. **`L^∞` is not an algebra under the product with a derivative.**  What the
   estimate needs is a space where `‖fg‖ ≲ ‖f‖‖g‖` *and* the norm dominates the
   gradient.  `H^s` with `s > d/2` does both.

In the axisymmetric-with-swirl formulation the ambient dimension for the
operator `L₅` is `d = 5`, so the algebra threshold is `s > 5/2`.  That is where
the `5/2` in the title comes from; it is not a stylistic choice.

## 2. The equation for the error

Let `u` solve Navier–Stokes on `ℝ³` and let `ũ` be the computed field, extended
by zero outside the truncated box.  Write

```
e = u - ũ ,        R = ũ_t + P(ũ·∇ũ) - νΔũ .
```

`R` is the **residual**: everything by which the computed field fails to solve
the equation.  Subtracting,

```
e_t + P(u·∇u) - P(ũ·∇ũ) = νΔe - R .
```

Applying the product-difference identity already proved in Lean
(`NSSingularity.product_difference`, `ab - ãb̃ = (a-ã)b + ã(b-b̃)`) to the
bilinear term:

```
u·∇u - ũ·∇ũ = e·∇u + ũ·∇e .                                       (2.1)
```

*Proved.*  This is the same identity F-18 uses, applied to the vector bilinear
form instead of the scalar product.  It is exact, not an estimate.

So

```
e_t = -P(e·∇u) - P(ũ·∇e) + νΔe - R .                              (2.2)
```

## 3. The `H^s` energy step, term by term

Pair (2.2) with `e` in `H^s` — that is, apply `Λ^s = (I-Δ)^{s/2}`, take the
`L²` inner product with `Λ^s e`, and use `d/dt ‖e‖_{H^s}² = 2⟨Λ^s e_t, Λ^s e⟩`.

### 3.1 The pressure term — *proved*

`P` is an orthogonal projection on `L²` and commutes with `Λ^s`, being a Fourier
multiplier.  Hence

```
⟨Λ^s P F, Λ^s e⟩ = ⟨Λ^s F, Λ^s P e⟩ = ⟨Λ^s F, Λ^s e⟩
```

whenever `div e = 0`.  **The pressure contributes nothing** — provided `e` is
divergence free, which is §3.2.

### 3.2 The divergence-free condition — *required lemma* `HS-1`

`u` is divergence free exactly; `ũ` is divergence free only up to the discrete
divergence measured by the certificate.  So `div e = -div ũ =: δ`, and §3.1
acquires a correction

```
|⟨Λ^s P F, Λ^s e⟩ - ⟨Λ^s F, Λ^s e⟩| ≤ C_P ‖F‖_{H^s} ‖δ‖_{H^{s-1}} .
```

> **Required lemma HS-1.**  There is `C_P = C_P(s)` such that for `F ∈ H^s` and
> `e ∈ H^s` with `div e = δ`,
> `|⟨Λ^s (I-P) F, Λ^s e⟩| ≤ C_P ‖F‖_{H^s} ‖δ‖_{H^{s-1}}`.
>
> *Unproved constant:* `C_P(s)`.  *What supplies it:* `(I-P)F = ∇Δ^{-1} div F`
> and the `H^s` boundedness of the Riesz transforms; the `Δ^{-1} div` gains one
> derivative, which is where `H^{s-1}` comes from.
>
> *What the numerics currently give:* `relative_divergence` in
> `outputs/tau_continuation_gate7_v1/checkpoints.csv`, a discrete `L^∞` ratio,
> **not** an `H^{s-1}` norm of `δ`.  Converting the one into the other is itself
> an unproved step (`HS-5`).

### 3.3 The advection terms — *required lemma* `HS-2`

The two terms behave differently and must not be lumped together.

**`ũ·∇e`, the transport term.**  Integrating by parts,

```
⟨Λ^s(ũ·∇e), Λ^s e⟩ = ⟨[Λ^s, ũ·∇]e, Λ^s e⟩ - (1/2)∫ (div ũ) |Λ^s e|² ,
```

so the leading term cancels and only the **commutator** survives.  This is the
Kato–Ponce structure, and it is the only reason the estimate closes at all.

> **Required lemma HS-2a (Kato–Ponce commutator).**  For `s > 0`,
> `‖[Λ^s, f·∇]g‖_{L²} ≤ C_{KP}(s) (‖∇f‖_{L^∞} ‖g‖_{H^s} + ‖f‖_{H^s} ‖∇g‖_{L^∞})`.
>
> *Unproved constant:* `C_{KP}(s)`.  *What supplies it:* the standard
> Kato–Ponce/Coifman–Meyer commutator estimate.  It is classical; it is listed
> because this repository has not proved it and mathlib does not have it.

**`e·∇u`, the stretching term.**  No cancellation is available: this term is
genuinely `‖e‖ · ‖∇u‖`.

> **Required lemma HS-2b (algebra property).**  For `s > 5/2`,
> `‖fg‖_{H^s} ≤ C_A(s) ‖f‖_{H^s} ‖g‖_{H^s}`, and `‖∇u‖_{H^s} ≤ ‖u‖_{H^{s+1}}`.
>
> *Unproved constant:* `C_A(s)`.  *What supplies it:* the Sobolev embedding
> `H^s ↪ L^∞` for `s > d/2` plus Leibniz.  This is why the estimate is stated
> with `‖u‖_{H^{s+1}}` on the right and not `‖u‖_{H^s}`: **the estimate loses one
> derivative**, and that loss is not removable by any rearrangement of these
> steps.

### 3.4 The viscous term — *proved, and it has a sign*

```
ν⟨Λ^s Δe, Λ^s e⟩ = -ν ‖∇Λ^s e‖_{L²}² ≤ 0 .
```

This is the one term that helps, and it is exact.  It is worth keeping rather
than discarding, because it absorbs part of the commutator: by Young,

```
C ‖∇Λ^s e‖ ‖Λ^s e‖ ≤ ν ‖∇Λ^s e‖² + (C²/4ν) ‖Λ^s e‖² ,
```

which trades a gradient for a factor `1/ν`.  With `ν` fixed and positive — as it
is in every run of this repository — that trade is legitimate.  **It would not
be legitimate under a vanishing or time-dependent viscosity protocol**, which is
one more reason those are excluded from the candidate computations.

### 3.5 Assembling — *conditional*

Under `HS-1`, `HS-2a` and `HS-2b`, and writing
`M(t) = ‖u‖_{H^{s+1}} + ‖ũ‖_{H^{s+1}}`,

```
d/dt ‖e‖_{H^s} ≤ C_s M(t) ‖e‖_{H^s} + ‖R‖_{H^s} + C_P ‖δ‖_{H^{s-1}} M(t)   (3.1)
```

with `C_s = C_s(C_{KP}, C_A, ν)`.  This is the target inequality of the session
instruction, with the two extra pieces the honest derivation produces: the
divergence defect and its coefficient.

Grönwall then gives, on `[0,T]`,

```
‖e(T)‖_{H^s} ≤ e^{C_s ∫₀^T M} ( ‖e(0)‖_{H^s}
                + ∫₀^T (‖R‖_{H^s} + C_P M ‖δ‖_{H^{s-1}}) dt ) .        (3.2)
```

*The Grönwall step itself is proved* — `NSSingularity.norm_le_simple_gronwall`
in `CertificateLayer.lean` is exactly this, for a constant coefficient.  The
time-dependent-coefficient version is `HS-6` below.

## 4. The residual, decomposed

`R` is not one thing.  Bounding it means bounding four separate objects, and
conflating them is how an `L^∞` certificate comes to look stronger than it is.

| piece | what it is | current status |
|---|---|---|
| `R_disc` | the discrete operator minus the continuous one | measured at nodes only; `HS-5` |
| `R_time` | the integrator's local defect | measured exactly (`slab_certificate.json`, trapezoid `2.35e-6`, Simpson `8.25e-10`) |
| `R_trunc` | truncating `ℝ³` to the box | bounded by the Green tail (`green_derivative_tail_bounds.md`) in `L^∞`, not `H^s`; `HS-4` |
| `R_embed` | the initial embedding error `e(0)` | zero for these runs — the datum is evaluated exactly on the grid |

> **Required lemma HS-4 (truncation in `H^s`).**  The multipole tail bounds of
> `free_space_recovery.py` are pointwise.  Converting them to `‖R_trunc‖_{H^s}`
> needs a bound on `s` derivatives of the tail, which the Green-kernel
> differentiation constants `A_m` do supply for integer `m` up to 6 — but `s`
> here is fractional and the conversion is not written.
>
> *Unproved constant:* the `H^s` tail constant.

> **Required lemma HS-5 (discrete to continuous).**  Every measured quantity in
> this repository is a grid quantity: `relative_divergence`, `poisson_residual`,
> the front counts.  Turning "the discrete residual is ≤ ρ at the nodes" into
> "‖R_disc‖_{H^s} ≤ ϱ(ρ, h)" needs an inverse inequality plus a bound on a
> higher derivative of the true solution.
>
> *Unproved constant:* `ϱ`.  **This is the largest single gap.**  It is the step
> that would make any of this a statement about the PDE rather than about the
> computation, and it is not done.

> **Required lemma HS-6 (Grönwall with a time-dependent coefficient).**  The
> Lean Grönwall step is proved for constant `K`.  (3.2) needs
> `∫₀^T M(t) dt`.  mathlib has `norm_le_gronwallBound_of_norm_deriv_right_le`
> for constant coefficients; the variable-coefficient form is a short extension
> but is not written.

## 5. What `s > 5/2` costs

The estimate needs `u ∈ H^{s+1}` on the right, i.e. `H^{7/2+}`.  For the
computed field that is a statement about a grid function and is meaningless
without `HS-5`.  For the true solution it is exactly the regularity whose
persistence is in question.  **The estimate is therefore conditional on the
quantity it would be used to control** — which is normal for a local-in-time
stability estimate and is precisely why (3.2) can only ever certify a *short*
interval, never a blow-up.

This is worth stating plainly: even fully proved, (3.2) could never establish a
Clay singularity.  It could establish that a computed trajectory stays close to
*some* exact solution over a short interval.  Chaining such intervals up to a
putative blow-up time fails because `∫M` diverges there by construction.

## 6. What would have to change in the Lean layer

`CertificateLayer.lean` currently proves the `L^∞` chain F-17–F-19.  Those
theorems are correct as stated — they are statements about real-valued
functions satisfying a differential inequality, not claims about Navier–Stokes.
What would change is their *interpretation*:

* F-17 (potential error → velocity error) survives unchanged; it is algebra.
* F-18 (product difference) survives unchanged and is reused as (2.1).
* F-19 (Grönwall) needs the variable-coefficient extension `HS-6`.
* A new layer would be needed for §3.1–3.4, and it needs `HS-1`, `HS-2a`,
  `HS-2b`, none of which mathlib currently has.

Estimated missing mathlib API: Kato–Ponce, the Leray projector's `H^s`
boundedness, and fractional Sobolev multiplication. That is a substantial body
of analysis, not a formalisation detail.

## 7. What the existing certificate does and does not establish

Stated once, flatly, so that no reader has to reconstruct it:

* `snapshot_certificate.py` and `slab_certificate.py` produce **exact rational
  enclosures of computed quantities**.  Under the two hypotheses `H1` and `H2`
  named in `slab_certificate.py`, they bound those quantities over whole cells
  and whole time slabs.
* They are **not** a proof about the continuous Navier–Stokes equations.  The
  step that would connect them is `HS-5`, which is not done.
* The `L^∞` Grönwall chain in Lean is **not** an unconditional PDE stability
  proof, and the repository does not claim it is.  It is a proved implication
  whose hypotheses — a pointwise differential inequality with a given constant —
  are supplied by hand, not derived from Navier–Stokes.

## 8. Summary of the required lemmas

| id | statement | supplies | status |
|---|---|---|---|
| HS-1 | `(I-P)` correction bounded by `‖δ‖_{H^{s-1}}` | `C_P(s)` | not proved here |
| HS-2a | Kato–Ponce commutator | `C_{KP}(s)` | classical, not proved here |
| HS-2b | `H^s` algebra for `s > 5/2` | `C_A(s)` | classical, not proved here |
| HS-4 | pointwise tail → `H^s` tail | `H^s` truncation constant | not written |
| HS-5 | discrete residual → `‖R‖_{H^s}` | `ϱ(ρ,h)` | **largest gap** |
| HS-6 | Grönwall with variable coefficient | — | short, not written |

Proved here and reusable: (2.1), §3.1, §3.4, and the constant-coefficient
Grönwall step already in Lean.
