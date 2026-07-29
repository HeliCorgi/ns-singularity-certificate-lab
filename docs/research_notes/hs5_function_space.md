# HS-5: fixing the function space before implementing anything

**Status: a decision document.**  It settles which norm the error estimate is
measured in, and why, before any of `norm_transfer.py`, the extended axis
checker, or the control-ODE layer is written.  Nothing here is a new theorem
about Navier–Stokes.  One thing here *is* a theorem — Theorem NT of §3 — and it
is an exact isometry rather than an estimate, which is the reason the whole
decision is available at zero cost in constants.

## 0. What this supersedes, and why

This note **supersedes the choice of function space made in
[`hs_error_propagation.md`](hs_error_propagation.md)**.  That note chose `H^s(ℝ⁵)`
with real `s > 5/2`, on the stated ground that "in the axisymmetric-with-swirl
formulation the ambient dimension for the operator `L₅` is `d = 5`, so the
algebra threshold is `s > 5/2`" (§1 of that note).  The reasoning was that if
`L₅` behaves like a five-dimensional Laplacian then the Sobolev calculus should
be done in five dimensions.

That inference is wrong, and it is wrong in a way that would have cost a great
deal of implementation effort before it surfaced.  The five-dimensional lift is
an exact tool for the **linear** operator `L₅` and for nothing else.  The moment
the nonlinearity is written down, five dimensions stops being a convenience and
becomes an obstruction: the lifted field is not divergence free, so the transport
term is compressible and the antisymmetry that makes every energy estimate work
is lost; there is no five-dimensional pressure; and the two generators `u₁` and
`ψ₁` do not sit at the same differential order, so a "common `H^s`" is silently
off by one derivative in the direction that hurts.  §5 gives this in full, with
the one open lemma (`NT-N1`) that anyone who still wants the five-dimensional
lane must supply first.

**The replacement decision.**  The error is measured in `H^n(ℝ³)`, `n` an
**integer**, `n = 3` or `n = 4`, on the **physical three-dimensional velocity**.
The five-dimensional lift survives only as a proof device inside Theorem NT, and
appears nowhere in the statement of any estimate.

Everything else in `hs_error_propagation.md` stands, with `s` replaced by the
integer `n` and the ambient dimension replaced by `3`.  In particular §2 (the
error equation), §3.1 (the pressure term contributes nothing to a divergence-free
error), §3.4 (the viscous term has a sign, and buys a factor `1/ν` by Young), §5
(the estimate is conditional on the quantity it would control, so it can certify
a short interval and never a blow-up) and §7 (what the existing certificate does
and does not establish) are unaffected by the change of dimension.  The required
lemmas `HS-1`, `HS-2a`, `HS-2b`, `HS-4`, `HS-5`, `HS-6` are unaffected in
substance and are restated with their new indices in §6.  §8 of that note — its
summary table of required lemmas — is superseded by §6 here.

## 1. The decision

> The true PDE error is measured in `H^n(ℝ³)` for **integer** `n = 3` or `n = 4`,
> on the **physical** three-dimensional velocity field `u`.  Not in `H^s` of the
> five-dimensional lift, and not in any weighted meridional norm.

Three separate choices are being made, and each has its own reason.

### 1.1 The physical field, not a generator

The Clay statement, the Escauriaza–Seregin–Šverák endpoint criterion used by
`critical_l3.py`, the Beale–Kato–Majda criterion, the Sobolev embeddings, and
every continuation theorem in the literature are statements about `u` on `ℝ³`.
A bound on `u₁` or `ψ₁` in some weighted half-plane norm is a bound on a
*coordinate representation*, and converting it costs either a constant or a
derivative unless the conversion is exact.  Theorem NT makes it exact, which is
why the choice of the physical field costs nothing.  Without Theorem NT it would
cost a constant at every level of the hierarchy, and those constants compound.

### 1.2 Why `n > d/2 + 1 = 5/2`, so `n ≥ 3`

With `d = 3` the two thresholds that matter are `d/2 = 3/2` and `d/2 + 1 = 5/2`,
and the distance between them is exactly one derivative.  They buy different
things and confusing them is how an unclosable hierarchy gets mistaken for a
closed one.

**`n > d/2` buys the product estimate.**  For `n > d/2`, `H^n` is an algebra and
the Leray-projected bilinear form `P(v, w) = ℙ((v·∇)w)` satisfies

```
‖P(v, w)‖_n ≤ K_n ‖v‖_n ‖w‖_{n+1} .                                     (1.1)
```

The `n+1` is not an artefact of a lazy proof.  `P(v, w)` differentiates `w` once,
and no rearrangement of Leibniz recovers it, because there is nothing to
integrate by parts against: (1.1) is a bound on a *norm*, not on a pairing.  Used
in a Grönwall argument, (1.1) gives `d/dt‖e‖_n ≲ ‖e‖_{n+1} · (…)`, so controlling
level `n` requires level `n+1`, which requires level `n+2`.  The hierarchy never
closes.  This is exactly the loss recorded at `HS-2b` in
`hs_error_propagation.md` — "the estimate loses one derivative, and that loss is
not removable by any rearrangement of these steps" — and the reason it is not
removable *at that threshold* is that `n > d/2` is simply not enough regularity.

**`n > d/2 + 1` buys the Kato inequality.**  The extra derivative converts the
norm bound into a bound on the *pairing against the second argument*:

```
|⟨P(v, w), w⟩_n| ≤ G_n ‖v‖_n ‖w‖_n²  .                                  (1.2)
```

Here the top-order term does cancel — it is the same cancellation as
`⟨(v·∇)w, w⟩ = -½∫(div v)|w|²` at the level of `Λ^n`, leaving only a commutator —
and the commutator can be absorbed because `∇v ∈ L^∞` is now available from
`H^{n-1} ↪ L^∞`, which needs `n - 1 > d/2`, that is `n > d/2 + 1`.  (1.2) is
what closes the control ODE **in `R_n` alone**, with no reference to level
`n+1` on the right-hand side of the quadratic term.  That is the whole
difference, and it is the reason the threshold in Morosi–Pizzocchero and in every
descendant of that argument is `n > d/2 + 1` and not `n > d/2`; see
[`a_posteriori_frameworks.md`](a_posteriori_frameworks.md) §1.

For `d = 3` the least integer above `5/2` is `n = 3`.

### 1.3 Why an integer

A fractional `H^s` norm is defined by the Fourier multiplier `Λ^s = (I - Δ)^{s/2}`,
or equivalently by a Gagliardo double integral.  Both are **nonlocal**.  A grid
of nodes carries no information about either without a global argument, and
an interval-arithmetic enclosure of a nonlocal quantity from finitely many nodes
is not something this repository can produce: `snapshot_certificate.py` and
`slab_certificate.py` enclose *local* expressions over cells, exactly, by
rational arithmetic with outward rounding.  For integer `n`,

```
‖u‖²_{H^n(ℝ³)} = Σ_{|α| ≤ n} ‖∂^α u‖²_{L²(ℝ³)}
```

is a finite sum of `L²` norms of derivatives.  Each summand is an integral of a
square of a local differential expression, which is precisely the shape an
interval quadrature over cells can bound.  That is the entire argument for
integrality, and it is a practical argument, not an aesthetic one: fractional `s`
would make `HS-5` — already the largest gap — unapproachable rather than merely
hard.

### 1.4 Why `n = 3` *or* `n = 4`

`n = 3` is the minimum the Kato inequality permits, and it is the level at which
Morosi–Pizzocchero state all their `d = 3` applications, so it is the level at
which any borrowed constant is meaningful.  It is also the level at which the
existing axis checker is already adequate (§7).

`n = 4` buys one spare derivative.  That spare derivative is wanted in three
places: the Kato–Ponce commutator (`HS-2a`) is cleaner when `∇u ∈ H^{n-1}` with
room to spare; the truncation tail (`HS-4`) is bounded through the Green-kernel
differentiation constants `A_m`, which `green_derivative_tail_bounds.md` supplies
for integer `m` up to 6, so `n = 4` is comfortably inside the range while
fractional `s` was not; and the discrete-to-continuum step (`HS-5`, `H3`) needs a
bound on a derivative one order above the norm, so working at `n = 3` means
needing four derivatives of the true solution anyway.

The price of `n = 4` is an extra axis condition — `∂_r³ψ₁(0, z) = 0` — which
nothing currently checks.  That is §7, and it is the reason this document exists
before the implementation rather than after it.

## 2. The five norms, kept strictly apart

Five different norms are in play, they are all called "the Sobolev norm" in
casual speech, and three of them differ by powers of `r`.  Each is written here
with its measure spelled out.  Nothing below should ever be written without its
measure.

Throughout, cylindrical coordinates are `x = (r cos θ, r sin θ, z)`, the
meridional half-plane is `H = {(r, z) : r > 0, z ∈ ℝ}`, and

```
dρ = r dr dz            the physical meridional measure
dμ = r³ dr dz           the ladder measure
⟨f, g⟩_μ = ∫_H f g r³ dr dz ,      ‖f‖_μ² = ⟨f, f⟩_μ
⟨f, g⟩_ρ = ∫_H f g r  dr dz ,      ‖f‖_ρ² = ⟨f, f⟩_ρ
L₅ = ∂_rr + (3/r)∂_r + ∂_zz = r^{-3} ∂_r(r³ ∂_r ·) + ∂_zz
ω₁ = -L₅ψ₁ ,   u^θ = r u₁ ,   u^r = -r ∂_zψ₁ ,   u^z = 2ψ₁ + r ∂_rψ₁
```

`L₅` is `laplacian_5d_formal` in `operators.py`; the velocity recovery is
`PicardLadder.velocity` in `picard_continuation.py`.

### (i) The physical three-dimensional Sobolev norm

```
‖u‖²_{H^n(ℝ³)} = Σ_{|α| ≤ n} ∫_{ℝ³} |∂^α u|² dx ,     dx = r dr dθ dz .
```

This is the norm the decision of §1 selects.  It is a norm of a **vector field**,
so the Cartesian components are differentiated, not the cylindrical ones — the
distinction matters because `e_θ` is not a constant field and differentiating it
produces the `1/r` terms that make the axis analysis nontrivial.

Two conventions for `H^n` circulate and they are *not* equal, only equivalent.
Theorem NT is stated for

```
‖u‖²_{H^n} := ‖u‖²_{L²} + ‖u‖²_{Ḣ^n} ,
```

which is the Fourier definition with weight `(1 + |ξ|^{2n})`.  This is the
convention in which the constant in Theorem NT is exactly 1.  Under the
full-sum convention above there are intermediate orders and the identity acquires
equivalence constants.  **Use the two-term convention**, and say so in every
docstring, because the point of Theorem NT is that no constant is lost.

### (ii) The two weighted meridional norms

```
‖f‖_ρ² = ∫_H f² r dr dz          (physical: dx = 2π dρ after the θ integration)
‖f‖_μ² = ∫_H f² r³ dr dz         (ladder: dV₅ = |S³| dμ, and ‖V_f‖²_{L²(ℝ³)} = 2π‖f‖²_μ)
```

**Confusing these two is the single easiest error in this entire project.**  They
differ by `r²`, they are both natural, both appear in the literature, and both
appear in this repository's own history.  A norm written `∫ f² r dr dz` in a
docstring is `ρ`; a norm written `∫ f² r³ dr dz` is `μ`.  The relation is exactly

```
‖f‖_μ = ‖r f‖_ρ ,                                                        (2.1)
```

so multiplication by `r` is an **isometric isomorphism** `L²(μ) → L²(ρ)` — not a
bounded map on either space alone.  §4 is entirely about why `μ` is the right one
and `ρ` is a trap.

### (iii) The five-dimensional lift norm

For `X = (X′, z) ∈ ℝ⁴ × ℝ = ℝ⁵` write `r = |X′|` and define the lift
`T₅f(X) = f(|X′|, z)`.  Then

```
‖T₅f‖²_{L²(ℝ⁵)} = |S³| ∫_H f² r³ dr dz = 2π² ‖f‖²_μ ,
Δ_{ℝ⁵} T₅f = T₅(L₅f) .                                                  (2.2)
```

The second identity is why `L₅` is called the five-dimensional Laplacian: the
radial part of the `ℝ⁴` Laplacian in `|X′|` is `∂_rr + (3/r)∂_r`.  This lift is
used **only inside the proof of Theorem NT**, to identify `L₅` as a genuine
Laplacian and thereby to get the spectral calculus of `L₅^{k/2}` for free.  It is
not used to state any estimate, for the reasons in §5.

### (iv) The generator norms

For integer `k ≥ 0`,

```
E_k(f) = ‖L₅^{k/2} f‖_μ                    k even
E_k(f) = ‖∇ L₅^{(k-1)/2} f‖_μ              k odd,   |∇g|² = g_r² + g_z²
```

where `∇` here is the **meridional** gradient `(∂_r, ∂_z)`, taken in the `μ`
inner product.  The two cases are consistent because integration by parts in
`L²(μ)` gives

```
⟨-L₅f, f⟩_μ = ‖∇f‖²_μ ,                                                  (2.3)
```

which is exactly the divergence-form structure `L₅ = r^{-3}∂_r(r³∂_r·) + ∂_zz`
and is the first of the four reasons in §4.

**Never expand `L₅^m` term by term for `n ≥ 3`.**  The individual terms of the
expansion carry negative powers of `r` and diverge separately at the axis; only
the assembled operator is finite.  This is not a numerical-conditioning remark,
it is a statement about the integrals: the expansion of `L₅²` contains
`9 r^{-2} f_r²`-type terms whose `μ`-integrals need not converge for admissible
`f`, while `‖L₅f‖_μ` does.  Any implementation must apply `laplacian_5d_formal`
repeatedly, never a symbolic expansion of its square.

### (v) The recovered Cartesian velocity

`reconstruct_cartesian_velocity` in `cartesian.py` produces `(u_x, u_y, u_z)` on
a Cartesian grid from `(u₁, ψ₁)`.  This is the object the Clay statement is about,
the object on which `‖u‖_{H^n(ℝ³)}` of (i) is literally defined, and the object
on which the truncation tail bounds of `free_space_recovery.py` and
`green_derivative_tail_bounds.md` are stated.  It is also the only one of the five
that is a vector field on a flat space with no weight and no axis, which is why it
is the right place to state the final estimate and the wrong place to compute.

### The dictionary

| object | lives in | measure | where in the repository |
|---|---|---|---|
| `u`, physical velocity | `H^n(ℝ³)`, `n = 3` or `4` | `dx = r dr dθ dz` | `cartesian.reconstruct_cartesian_velocity` |
| `u₁`, swirl generator | `L²(μ)`, and `E_n(u₁)` | `dμ = r³ dr dz` | state array on `AxisymmetricGrid` |
| `ψ₁`, stream generator | `E_{n+1}(ψ₁)` | `dμ` | `WholeSpaceEllipticSolver.solve` output |
| `ω₁ = -L₅ψ₁` | `E_{n-1}(ω₁)` | `dμ` | state array; `PicardLadder` integrates it |
| `u^θ = r u₁`, `u^r`, `u^z` | components of `u`; `L²(ρ)` per component | `dρ = r dr dz` | `PicardLadder.velocity` |
| `T₅u₁`, `T₅ψ₁` | `L²(ℝ⁵)` | `dV₅ = |S³| r³ dr dz` | proof device only; no code |
| elliptic source mass | `∫ω₁ dV₅` | `dμ` | `whole_space_gate.py` lines 402, 410, 873 |
| radial quadrature weights | — | `dμ` | `free_space_poisson.py` line 329 |

## 3. Theorem NT

> **Theorem NT (norm transfer).**  Let `u` be an axisymmetric-with-swirl vector
> field on `ℝ³` with generators `u₁`, `ψ₁` in the sense
> `u = V_{u₁} + curl V_{ψ₁}`, where `V_f = r f e_θ`, and let `ω₁ = -L₅ψ₁`.  Then
> for every integer `n ≥ 0`,
>
> ```
> ‖u‖²_{Ḣ^n(ℝ³)} = 2π [ E_n(u₁)² + E_{n+1}(ψ₁)² ]                        (NT-Ḣ)
> ‖u‖²_{H^n(ℝ³)}  = 2π [ ‖u₁‖²_μ + E_n(u₁)² + ⟨ψ₁, ω₁⟩_μ + E_{n-1}(ω₁)² ] (NT-H)
> ```
>
> with `‖u‖²_{H^n} := ‖u‖²_{L²} + ‖u‖²_{Ḣ^n}`.
>
> **These are equalities.**  The transfer is an exact isometry with constant 1;
> no constant is lost, in either direction, at any order.

That last sentence is the reason the decision of §1 is free.  A merely equivalent
pair of norms would force a constant `C_n` into the control ODE at every use, and
those constants multiply through the Grönwall exponential, so an equivalence with
`C_n = 2` at `n = 4` would cost a factor `2^{…}` in the certified interval.  An
isometry costs nothing.

The two forms are consistent: `ω₁ = -L₅ψ₁` gives
`L₅^{(n+1)/2}ψ₁ = -L₅^{(n-1)/2}ω₁`, hence `E_{n+1}(ψ₁) = E_{n-1}(ω₁)`, and
`‖u‖²_{L²} = 2π[‖u₁‖²_μ + ⟨ψ₁, ω₁⟩_μ]` is (NT-H) at the bottom order.  The
`⟨ψ₁, ω₁⟩_μ` term is the meridional kinetic energy written the way the solver
already produces it — the elliptic solve returns `ψ₁` from `ω₁`, so the pairing
is available without differentiating anything.

### The two cases written out

For `n = 3`:

```
‖u‖²_{Ḣ³} = 2π [ ‖∇L₅u₁‖²_μ + ‖L₅²ψ₁‖²_μ ]
‖u‖²_{H³} = 2π [ ‖u₁‖²_μ + ‖∇L₅u₁‖²_μ + ⟨ψ₁, ω₁⟩_μ + ‖L₅ω₁‖²_μ ]
```

For `n = 4`:

```
‖u‖²_{Ḣ⁴} = 2π [ ‖L₅²u₁‖²_μ + ‖∇L₅²ψ₁‖²_μ ]
‖u‖²_{H⁴} = 2π [ ‖u₁‖²_μ + ‖L₅²u₁‖²_μ + ⟨ψ₁, ω₁⟩_μ + ‖∇L₅ω₁‖²_μ ]
```

Read off the implementation requirement: `n = 3` needs `L₅` applied twice to `ω₁`
(equivalently three times to `ψ₁`) and once to `u₁` followed by a meridional
gradient; `n = 4` needs `L₅` twice on `u₁` and `L₅` once on `ω₁` followed by a
gradient.  Nothing beyond repeated `laplacian_5d_formal`, `derivative_r`,
`derivative_z`, and a `μ`-weighted quadrature.  That is the whole of
`norm_transfer.py`.

### The proof chain, by name

The proof is a chain of five exact steps, each of which is either classical or a
one-line computation.  They are named here so that the implementation and any
future Lean layer can be checked step by step against them.

1. **Polar factorisation.**  For an axisymmetric integrand the `θ` integral
   factors out as `2π`, reducing every `ℝ³` integral to a `ρ`-weighted integral
   over the half-plane.  This is where the `2π` in Theorem NT comes from, and it
   is the only place it comes from.
2. **The two Laplacian intertwinings.**  `Δ_{ℝ³} V_f = V_{L₅f}` and
   `Δ_{ℝ⁵} T₅f = T₅(L₅f)`.  The first is the direct computation
   `Δ(A e_θ) = (A_rr + A_r/r - A/r² + A_zz)e_θ` with `A = r f`, which collapses
   to `r(f_rr + 3f_r/r + f_zz)`; the `-A/r²` term is exactly what turns the
   three-dimensional radial Laplacian into the five-dimensional one, and this is
   the honest origin of the "five dimensions" in the whole formulation.  The
   second is (2.2).
3. **The weight identity `‖f‖_μ = ‖r f‖_ρ`**, equation (2.1).  This is what
   converts a `ρ`-weighted statement about the physical component `u^θ = r u₁`
   into a `μ`-weighted statement about the generator `u₁`, with no loss.
4. **Swirl/meridional orthogonality.**  `V_{u₁}` is azimuthal and `curl V_{ψ₁}`
   is meridional, so they are pointwise orthogonal; and `Δ` preserves the
   splitting by step 2, so the orthogonality survives every application of
   `Λ^n`.  Hence the squared norms add with no cross term at any order.  This is
   the step that makes Theorem NT a sum of two independent contributions rather
   than a quadratic form with off-diagonal entries.
5. **The decomposition `u = V_{u₁} + curl V_{ψ₁}`**, which is the standard
   axisymmetric-with-swirl representation, valid for smooth solenoidal
   axisymmetric fields (Liu–Wang, `REFERENCES.md` `[LiuWang2009]`, text
   verified).  Its validity at finite regularity, and the axis conditions it
   requires, are `WH-1` and `WH-1f` of §6.

Step 2's second identity is used only to get the functional calculus of
`L₅^{k/2}` — self-adjointness, non-positivity, and the spectral definition of
half-integer powers — cheaply.  It is not used to formulate anything.

## 4. Why `r³` and not `r`

Four reasons, of which the first two are structural and the last two are
practical.  Any one of them would settle it; together they close the question.

**(a) `L₅` is self-adjoint and non-positive on `L²(μ)`, and is not even symmetric
on `L²(ρ)`.**  In divergence form `L₅ = r^{-3}∂_r(r³∂_r·) + ∂_zz`, so for
decaying `f, g`

```
⟨L₅f, g⟩_μ = -∫_H (f_r g_r + f_z g_z) r³ dr dz = ⟨f, L₅g⟩_μ ,
```

which is symmetry and, at `g = f`, non-positivity — equation (2.3).  In `ρ` the
same integration by parts leaves a residue:

```
⟨L₅f, g⟩_ρ - ⟨f, L₅g⟩_ρ = 2 ∫_H (f_r g - f g_r) dr dz ≠ 0 .
```

The operator that is symmetric for `ρ` is the three-dimensional radial Laplacian
`∂_rr + (1/r)∂_r + ∂_zz`, not `L₅`.  Since every one of `E_k`, the spectral
powers `L₅^{k/2}`, and the identity `E_{n+1}(ψ₁) = E_{n-1}(ω₁)` rests on
self-adjointness, choosing `ρ` does not merely change constants — it destroys the
functional calculus that the norm definition uses.

**(b) Multiplication by `r` is an isometry `L²(μ) → L²(ρ)` and is unbounded on
either space alone.**  By (2.1) the map `f ↦ rf` is an isometric isomorphism
between the two weighted spaces.  It is *not* a bounded operator `L²(ρ) → L²(ρ)`
(take mass out to large `r`), nor is its inverse bounded `L²(ρ) → L²(ρ)` (take
mass to the axis).  So the factor `r` in `u^θ = r u₁` is free in `μ` and
uncontrollable in `ρ`.  Since that factor is present in the definition of the
swirl generator and cannot be removed, the measure must be the one that absorbs
it.

**(c) The repository already computes in `μ`.**  The `r³` weight is not a
proposal, it is the existing convention:

```
free_space_poisson.py:329      radial_weights = base_weights * (r**3)
whole_space_gate.py:402,410    weighted = values * grid.r[:, None] ** 3     (mass, L¹ source)
whole_space_gate.py:873        weights = grid.r[:, None] ** 3               (relative error)
wall_sensitivity.py:129,163    np.sum(weights * r**3 * source)
```

`whole_space_gate.py` documents the reason in place: `dV₅ = |S³| r³ dr dz`, and
the monopole boundary condition of the free-space elliptic solver consumes
exactly `∫ω₁ dV₅`.  Switching to `ρ` would silently invalidate the boundary
condition of the solver on which every whole-space run depends, and the failure
would appear as a slow drift in the far field rather than as an error.

**(d) The `ρ`-weighted norm of `u₁` is off by exactly one derivative.**  This is
the sharpest reason and it deserves its numbers.  For an azimuthal field
`A e_θ`,

```
‖∇(A e_θ)‖²_{L²(ℝ³)} = ∫_{ℝ³} ( |∇A|² + A²/r² ) dx ,
```

the `A²/r²` being the same term that produced `L₅` in step 2 of §3.  Taking
`A = u^θ = r u₁` gives `A²/r² = u₁²`, hence exactly

```
‖u₁‖²_ρ = (2π)^{-1} ∫_{ℝ³} |u^θ|²/r² dx ≤ (2π)^{-1} ‖∇u‖²_{L²(ℝ³)} .      (4.1)
```

So `‖u₁‖_ρ` is controlled by the **gradient** of `u`, never by `u` itself.  It is
a first-order quantity wearing a zeroth-order name.  And it is genuinely not
bounded by `‖u‖_{L²}`: take

```
u^θ = r g(z) χ(r/ε)
```

with `χ` a fixed cut-off.  Then `‖u‖²_{L²} ~ ∫ r² g² χ² r dr dz = O(ε⁴)` while
`‖u₁‖²_ρ = ∫ g²χ² r dr dz = O(ε²)`, so the ratio `‖u₁‖²_ρ / ‖u‖²_{L²}` is of
order `ε^{-2}` and diverges as the field concentrates on the axis.  In `μ` the
same field gives `‖u₁‖²_μ = ∫ g²χ² r³ dr dz = O(ε⁴)` and the ratio is exactly
`(2π)^{-1}`, as Theorem NT says it must be.

Using `ρ` would therefore not merely change constants; it would shift the whole
ladder by one derivative in the direction that hurts, and the shift would be
invisible until a run with a concentrating core — which is precisely the run this
repository is built to perform.

## 5. The obstruction: why the five-dimensional `H^s` lane is blocked

The five-dimensional lift is exact for `L₅` and for nothing else.  Five separate
obstructions block its use as the space in which the *nonlinear* estimate is
stated.  They are independent; fixing any one leaves the others.

### (a) The natural five-dimensional lift is not divergence free

Take the natural lift of the meridional velocity to `ℝ⁵ = ℝ⁴_{X′} × ℝ_z`,

```
U(X′, z) = ( u^r(r, z) X′/r ,  u^z(r, z) ) ,     r = |X′| .
```

The divergence of a radial field `f(r) X′/r` on `ℝ⁴` is `f′(r) + 3f(r)/r`, so

```
div₅ U = ∂_r u^r + 3u^r/r + ∂_z u^z .
```

Physical incompressibility on `ℝ³` is `∂_r u^r + u^r/r + ∂_z u^z = 0`.  Hence

```
div₅ U = 2 u^r / r  ≠ 0 .                                                (5.1)
```

**The five-dimensional transport is therefore compressible**, and the single
identity on which every energy estimate for Navier–Stokes rests,

```
⟨(U·∇)U, U⟩ = -½ ∫ (div₅ U) |U|² ,
```

no longer vanishes.  What it leaves behind is not a small perturbation.  By
(5.1) it is `-∫ (u^r/r)|U|²`, and `u^r/r` is precisely the coefficient of the
axisymmetric vortex-stretching term — the quantity whose sign and size is the
entire difficulty of the axisymmetric problem, the thing that makes `ω₁` grow.
So the five-dimensional formulation does not remove the difficulty; it removes
the *machinery* and leaves the difficulty sitting in a term that the machinery
can no longer absorb.  This is the decisive obstruction.

### (b) There is no five-dimensional pressure

The pressure is defined by the Leray projection on `ℝ³`, `p = R_iR_j(u_iu_j)`.
There is no field on `ℝ⁵` whose gradient is the lift of `∇p`: the azimuthal
structure and the `r`-weights are wrong, and (5.1) says the lifted flow is not in
the kernel of any five-dimensional Leray projection anyway.  Consequently the
"five-dimensional Navier–Stokes system" does not exist.  What exists is a system
of scalar transport–diffusion equations for `u₁` and `ω₁` with a nonlocal source,
which is a different object with different structure — in particular with none of
the orthogonality that makes `§3.1` of `hs_error_propagation.md` (the pressure
contributes nothing) true.

### (c) `u₁` and `ψ₁` do not sit at the same order

Theorem NT says it plainly: `u₁` enters `‖u‖_{Ḣ^n}` at order `n` and `ψ₁` at
order `n+1`.  A five-dimensional formulation that puts the pair `(u₁, ψ₁)` in a
common `H^s(ℝ⁵)` therefore either wastes a derivative on `u₁` or under-resolves
`ψ₁` by one.  Both are silent.  The second is the one that occurs in practice,
because the natural instinct is to bound both generators at the level the norm
names, and the second is the unsafe direction: `ψ₁` is what generates the
meridional velocity that does the advecting, so under-resolving it under-resolves
the nonlinearity itself.

### (d) Multiplication by `r` is not an `H^s(ℝ⁵)` multiplier, and the obstruction
is at the **axis**, not at infinity

The transfer between physical components and generators is multiplication by `r`
in one direction and by `1/r` in the other, and both directions are needed.  On
`ℝ⁵`, `r = |X′|` is Lipschitz but not `C²` at the axis `{X′ = 0}`, a
codimension-four subspace: `∇r = X′/r` is discontinuous there and `|D²r| ~ 1/r`.

What is actually verified here is the following, and it is worth being precise
because the folklore statement ("not a multiplier for `s ≥ 2`") is not what the
computation gives.  The first bad term, `(D²r) f ~ f/r`, *is* saved: the Hardy
inequality in codimension four gives `∫_{ℝ⁵} |f|²/r² dX ≤ ∫ |∇f|² dX`.  The next
one, `f/r²`, needs the Rellich inequality `∫ |f|²/r⁴ ≤ C ∫ |Δf|²`, whose constant
degenerates **exactly** in codimension four.  So the argument stops after one
derivative, and in particular the inverse map `f ↦ f/r` is already unbounded at
the first order.  The sharp threshold in `s` is not settled here, and settling it
is part of `NT-N1`.

What *is* settled is the qualitative point: the failure is **local at the axis**.
No truncation at infinity removes it — compact support does not help — and
truncating near the axis destroys the axis regularity on which the entire
representation `u = V_{u₁} + curl V_{ψ₁}` depends.  There is no cheap way out.

### (e) The embedding threshold is worse in five dimensions

`H^s(ℝ⁵) ↪ L^∞` needs `s > 5/2`; `H^n(ℝ³) ↪ L^∞` needs `n > 3/2`.  For exactly
the same pointwise information the five-dimensional lane costs a full extra
derivative, and in this repository every derivative costs an axis condition (§7),
a quadrature order, and a term in the discrete-to-continuum estimate.  Even if
(a)–(d) were all resolved, (e) alone would make the five-dimensional lane the
more expensive of the two.

### The open lemma

> **`NT-N1` (OPEN).**  A weighted commutator estimate for the dilation field
> `X′·∇_{X′}` on `H^s(ℝ⁵)`.  Concretely: identify the range of real `s` and a
> constant `C_{N1}(s)` such that, for `f` in the axis-regular subspace (the image
> of `T₅` on smooth even generators),
>
> ```
> ‖ [Λ^s_{ℝ⁵}, X′·∇_{X′}] f ‖_{L²(ℝ⁵)} ≤ C_{N1}(s) ‖f‖_{H^s(ℝ⁵)}
> ```
>
> together with the two-sided bound
> `C^{-1}‖f‖_{H^s(ℝ⁵)} ≤ ‖r f‖_{H^s(ℝ⁵)} ≤ C‖f‖_{H^s(ℝ⁵)}` on that subspace.
>
> `X′·∇_{X′} = r ∂_r` is the Euler operator, and it is the commutator that
> controls the interaction between multiplication by `r` and the fractional
> derivative.  Without it, none of (c) or (d) can be repaired.

**Anyone proposing to close the argument "in five dimensions" must supply `NT-N1`
in writing, with its constant and its range of `s`, before any code is written.**
This is recorded as a hard gate, not a preference.  The reason is that the
five-dimensional lane is attractive precisely because `L₅` looks like a
Laplacian, and that attraction survives the first four obstructions by sheer
plausibility; `NT-N1` is the smallest concrete object whose absence proves the
lane is not open.

## 6. Dependency table for HS-5

This supersedes §8 of `hs_error_propagation.md`.  Indices `s` there become
integer `n ∈ {3, 4}` here, and the ambient dimension is 3.

Status vocabulary, used strictly:

* **proved** — proved in this project, in prose, with the argument recorded.
* **classical, unformalised** — a standard theorem of analysis, true and
  citable, which neither this repository nor mathlib has a formal proof of.
* **open** — no proof is known to us, or the statement is not yet pinned down.
  These are the ones that block.

| id | statement | supplies | depended on by | status |
|---|---|---|---|---|
| **NT** | `‖u‖²_{H^n(ℝ³)} = 2π[‖u₁‖²_μ + E_n(u₁)² + ⟨ψ₁,ω₁⟩_μ + E_{n-1}(ω₁)²]`, and the `Ḣ^n` form | the exact isometry, constant 1, between physical and generator norms | everything: HS-1, HS-2a/b, HS-4, HS-5, H3, the control ODE | **proved** (§3), exact isometry |
| **WH-1** | Whitney even-function theorem with parameters: `f(r,z)` smooth and even in `r` implies `f = g(r², z)` with `g` smooth | legitimacy of `u^θ = r u₁` with smooth `u₁`, hence the axis and parity structure | NT, WH-1f, the parity selection rule, `axis.py` | **classical, unformalised** (Whitney 1943, with `z` as parameter); not in mathlib |
| **WH-1f** | finite-order version: `f ∈ C^{2m}` even in `r` implies `g ∈ C^m` with `‖g‖_{C^m} ≤ C_m‖f‖_{C^{2m}}`, and the `H^n` analogue | the quantitative axis conditions `∂_r^{2k+1}f(0,z) = 0` and the constants for NT at finite regularity | NT at finite regularity, §7's checker, HS-5 | **classical, unformalised**; the constants `C_m` are not written down anywhere we have read |
| **NT-N1** | weighted commutator estimate for `X′·∇_{X′}` on `H^s(ℝ⁵)`, plus the two-sided `r`-multiplier bound | the only route to a five-dimensional formulation | the 5D lane, and nothing else | **OPEN — the blocker for the 5D lane** (§5); the 3D lane does not need it |
| **HS-1** | `\|⟨Λ^n(I-ℙ)F, Λ^n e⟩\| ≤ C_P ‖F‖_{H^n}‖δ‖_{H^{n-1}}`, `δ = div e` | `C_P(n)`; the cost of the computed field not being exactly divergence free | the control ODE's pressure term | **classical, unformalised**; follows from `(I-ℙ)F = ∇Δ^{-1}div F` and Riesz boundedness |
| **HS-2a** | Kato–Ponce: `‖[Λ^n, f·∇]g‖_{L²} ≤ C_{KP}(‖∇f‖_∞‖g‖_{H^n} + ‖f‖_{H^n}‖∇g‖_∞)` | `C_{KP}(n)`; the transport term's only surviving piece | the control ODE's transport term; the Kato inequality (1.2) | **classical, unformalised**; not in mathlib |
| **HS-2b** | `H^n(ℝ³)` algebra for `n > 3/2`: `‖fg‖_{H^n} ≤ C_A‖f‖_{H^n}‖g‖_{H^n}` | `C_A(n)`; the stretching term `e·∇u` | the control ODE's quadratic term | **classical, unformalised**; Morosi–Pizzocchero's constants papers are on `T^d`, so their numbers do not transfer |
| **HS-4** | pointwise multipole tail bounds → `‖R_trunc‖_{H^n}` | the `H^n` truncation constant | the residual decomposition; the domain-truncation budget | **not written**; the Green-kernel constants `A_m` in `green_derivative_tail_bounds.md` cover integer `m ≤ 6`, so `n ≤ 4` is inside range — this is now tractable, which it was not for fractional `s` |
| **HS-5** | discrete nodal residual → `‖R_disc‖_{H^n}`: `ϱ(ρ, h)` | the step that makes any of this a statement about the PDE rather than the computation | every claim about `u` as opposed to `ũ` | **OPEN — the largest single gap**; needs an inverse inequality plus a bound on a derivative one order above `n` |
| **HS-6** | Grönwall with a time-dependent coefficient | the `∫₀^T M(t) dt` form of the exponential | the Lean certificate layer | **classical, short, not written**; mathlib has the constant-coefficient form, `norm_le_gronwallBound_of_norm_deriv_right_le` |
| **H3** | semi-discrete → continuum: the finite-difference operators' consistency error measured in `H^n`, not at nodes | the other half of what `HS-5` needs; distinct from it | HS-5 | **open**; overlaps HS-5 but is separable — `H3` is about the operator, `HS-5` is about the norm of the residual |
| **H2′** | Picard/Banach rough enclosure replacing the slab certificate's `H2` | a hypothesis-free slab enclosure: if `[Y_n] + [0,h]F(B) ⊂ int B` then the Picard operator maps `C([t_n, t_n+h], B)` into itself and the IVP has a unique solution on the slab staying in `B` | `slab_certificate.py`, whose docstring currently points at `hs_error_propagation.md` for `H2` | **proved**; the Hermite-remainder form of `H2` is unfixable as framed, because `u1_end` is the RK4 output rather than the same trajectory's value; not yet implemented |

Reading of the table: the **only** genuinely open items on the chosen
three-dimensional lane are `HS-5` and `H3`, and they are two halves of the same
discrete-to-continuum step.  `NT-N1` is open but is not on the chosen lane — it
blocks only the alternative that §5 rejects.  Everything else is classical
analysis that has to be written down, formalised, or both, but which nobody has
to invent.  That is a considerably better position than §8 of
`hs_error_propagation.md` described, and the improvement is entirely due to
Theorem NT removing the transfer constants and to integrality removing the
fractional calculus.

## 7. Bug notice: the axis checker is one condition short at `n = 4`

`axis.check_axis_regularity` in `src/ns_certificate_lab/axis.py` tests, for each
even field, `max|∂_r f|` at the axis, and for each odd field `max|f|` at the
axis.  Its helper `_axis_derivative_independent` computes a one-sided **first**
derivative only, at fourth order with five points or third order with four.

That is the `k = 0` condition.  The axis conditions that Theorem NT's finite-
regularity form requires are

```
∂_r^{2k+1} f(0, z) = 0     for every odd order 2k+1 the norm sees .
```

At `n = 3` the only such order is `2k+1 = 1`, and the existing check is
therefore adequate.  At `n = 4` it is not: the norm additionally requires

```
∂_r³ ψ₁(0, z) = 0 ,
```

**and nothing in the repository checks this.**  A datum that passes
`check_axis_regularity` today may have a nonzero third radial derivative of `ψ₁`
at the axis, in which case the `n = 4` transfer is being applied outside its
hypotheses and the resulting norm is not the `H⁴(ℝ³)` norm of anything.

Two honest caveats on the index range, recorded rather than papered over.  The
exact cut-off — whether it is `2k+1 ≤ n-2` or `2k+1 ≤ n-1` — differs precisely at
`n = 4`, which is the case that matters, and the two candidate statements were
both in circulation when this note was written.  The condition stated above is
the one consistent with `n = 3` being covered by the existing check and `n = 4`
requiring one more.  **It should be pinned down by `WH-1f` and not by assertion**,
and until it is, the `n = 4` lane should be treated as carrying an unverified
hypothesis.  Second, `∂_r³` at the axis needs a one-sided stencil the module does
not have: five points give the first derivative to fourth order, but a useful
third derivative needs at least six, so extending the checker means extending
`_axis_derivative_independent` as well, not merely calling it again.

A parallel task is extending the checker.  This note records the gap so that the
extension has a written specification to satisfy rather than an inferred one, and
so that any `n = 4` result produced before the extension lands is known to be
provisional.
