# Gate 8 — searching for a positive initial generation rate

Companion to [`l3_generation_rate.md`](l3_generation_rate.md), which derives the
identity and the pure-swirl no-go.  This note records what the search actually
found, including the parts that did not work.

**Nothing here bears on the Clay problem.**  A positive initial rate of the
critical norm is a screening criterion for a search lane, not evidence of
blow-up.

## 1. What the no-go explains

Gate 7 followed three whole-space families to `τ = 1` across a factor of forty
in Reynolds number and found the critical norm `‖u‖₃` decaying in **all eighteen
runs**.  That was not bad luck.  All three families — S, A and H — are pure
swirls, and for a pure swirl

* the advective contribution to `d/dt‖u‖₃³` is zero, as it is for *every*
  divergence-free field;
* the pressure contribution is zero **identically**, because `p` is axisymmetric
  and `u₀·∇p = u^θ r⁻¹∂_θ p = 0` pointwise;
* the viscous contribution is a sum of two nonnegative integrands.

So `J(u₀) := F'(0) ≤ 0` was forced, with equality only for `u₀ ≡ 0`.

**What the theorem covers, precisely.**  It guarantees only that the initial
generation rate is non-positive.  Pure swirl is not preserved: vortex stretching
generates meridional flow immediately, the pressure channel switches on at
`t = 0⁺` with indefinite sign, and the theorem then says nothing.  The sustained
decrease Gate 7 observed across `0 < τ ≤ 1` in all eighteen runs is a
**numerical observation** — consistent with the theorem, not implied by it.

This is a **kill**, not support: a pure swirl has `u ⊥ ∇|u|` identically, so it
is the exact zero of the only channel that can produce growth.

## 2. Leaving the pure-swirl lane

`mixed_initial_data.py` builds data with both swirl and meridional flow from two
generators `u₁(r²,z)` and `ψ₁(r²,z)`, giving

```
u^θ = r u₁ ,   u^r = −r ∂_z ψ₁ ,   u^z = 2ψ₁ + r ∂_r ψ₁ = 2ψ₁ + 2s ∂_s ψ₁ .
```

Divergence-freedom is an algebraic identity for **any** `ψ₁`, and every partial
derivative is available in closed form because each generator factors as
`a X(s) Z(z)`.  Both meridional components are therefore built from analytic
derivatives, and the measured discrete divergence is `3e-16` relative — the
difference operator's truncation error and nothing else.

### The parity selection rule

Before running anything, a parity argument decides whether the pressure term can
be nonzero at all.  With `σ_u, σ_ψ` the `z`-parities of the generators:

| object | parity |
|---|---|
| `u^θ` | `σ_u` |
| `u^r` | `−σ_ψ` |
| `u^z` | `σ_ψ` |
| `\|u\|` | **always even** (a sum of squares) |
| `∂_iu_j∂_ju_i` | **always even** (squares, plus two cross terms of equal-parity factors) |
| `p` | **always even** |
| pressure integrand `−3\|u\|(u^r∂_rp + u^z∂_zp)` | `−σ_ψ` |

On a `z`-symmetric domain the pressure term therefore **vanishes identically
unless `ψ₁` is odd in `z`**.  The swirl parity never enters.

This was found the hard way: the first M1 and M2 shipped with *even* stream
generators and returned `P = 1.5e-19` and `1.3e-17`, i.e. zero.  The rule is now
derived in the module docstring, implemented as `parity_selection_rule`, and
checked numerically against the actual integral, which is the only way to catch
a sign slip in the parity argument itself.

An odd `ψ₁` is also the physically right choice: it gives `u^r` even and `u^z`
odd, the axisymmetric strain flow that converges radially at every height and
ejects along both `±z`.  An even `ψ₁` produces inflow on one side and outflow on
the other, and the halves cancel.

## 3. The objective

Two exact symmetries collapse the search space.

**`P` is odd and `V` is even under `u ↦ −u`.**  `p` is quadratic in `u` so it is
unchanged, while `∇·(|u|u)` flips.  Measured: flipping the sign of a family
gives `P = ∓4.227834e-4` and `V = −4.952322e-2` **to every printed digit**.  So
for any datum with `P ≠ 0`, one of the two signs works, and the search reads it
off rather than exploring it.

**Amplitude is not a search direction.**  `P ~ A⁴L²` and `|V| ~ νA³L`, so
`P/|V| ~ AL/ν`.  Under the Navier–Stokes scaling `u ↦ λu(λx)` that preserves
`‖u‖₃`, `A ↦ λA` and `L ↦ L/λ`, so the ratio is **invariant**: `J>0` cannot be
manufactured by rescaling inside a critical family.

What remains is the shape, and the scale-free objective is the **critical
Reynolds number**

```
Re_crit(shape) = A L (|V|/ν) / P ,        J > 0  ⟺  Re > Re_crit .
```

One power of `L`, because `A` is the physical velocity amplitude.  This is the
same number as the `Re = A_{u₁}L²/ν` of `nondimensional.py`, whose amplitude is
the swirl *variable*, since `A = A_{u₁}L`.

The objective is deliberately **not** the maximum vorticity: `max|ω|` is not
invariant under the Navier–Stokes scaling, so maximising it is a statement about
units.

## 4. The search, and its preregistered kill condition

Eight amplitudes over a fixed basis of four swirl and four stream generators,
two of each centred away from `z = 0` so the symmetric subspace is not
invariant.  Multi-start projected-gradient descent on `log Re_crit`, with the
amplitude vector projected onto the unit sphere after each step — exactly,
because scale is not a direction, so rescaled duplicates are literally the same
point.

**Kill condition, recorded in `l3_optimizer.py` before the search ran:** a
generic band-limited divergence-free field gives `|P/V| ≈ 1e-3`, i.e.
`Re_crit ≈ 1e3`.  If shape optimisation cannot bring the best `Re_crit` below
`1e2`, the shape-factor lane is dead — and §3 forbids rescuing it by rescaling.

### Result

| stage | grid | `Re_crit` |
|---|---|---|
| hand-designed seed (M2, flipped) | 65×129 | 1.27e4 |
| optimiser, coarse | 49×97 | 1.13e3 |
| the same shape re-verified | 97×193 | 2.04e3 |
| optimiser, medium | 65×129 | 1.30e3 |
| the same shape re-verified | 97×193 | 1.37e3 |
| the same shape re-verified | 129×257 | **1.41e3** |

Two findings.

* **`J > 0` is achieved.**  The optimised shape has `P > 0` and gives a strictly
  positive initial generation rate for every `Re > ≈1.4e3`, confirmed directly:
  at half the critical Reynolds number the rate is negative and at twice it is
  positive, with the crossing where the formula puts it.
* **The preregistered kill condition fires.**  `1.41e3` is more than ten times
  the `1e2` bar.  Shape optimisation bought a factor of nine over the
  hand-designed seed and roughly nothing over a generic field.  By the rule
  recorded in advance, **the shape-factor lane is declared dead**.

A third finding, methodological: the coarse-grid optimum read `1.13e3` and
verified at `2.04e3` — the optimiser overfits the grid.  Any objective computed
through a discretisation must be re-verified at higher resolution before it is
believed, and the table above does that at every stage.

## 5. Domain and convergence

| quantity | behaviour |
|---|---|
| `P` | 2.31e-6, 2.33e-6, 2.27e-6, 2.25e-6 at `nr = 65, 97, 129, 161` |
| `V/ν` | 4.88e-2, 5.06e-2, 5.17e-2, 5.22e-2 — still creeping, `O(h^{1.5})` |
| transport residual | 4e-9 to 2e-8, i.e. 0.2–0.9 % of `P` |
| `Re_crit` across boxes `(2.5,3.0)` to `(4.2,5.0)` | 2093, 2042, 2059, 2045 |

The transport residual is the load-bearing check: the identity says it is
exactly zero, so a residual comparable to `P` would mean the quadrature cannot
resolve `P` either.  At 0.2–0.9 % it does not.  Domain independence is good to
2 %.

The finite-difference gradient converges to the analytic one at order ≈1.55 and
only from `nr ≈ 97` upward; below that the edge of the `χ` bump is unresolved and
the comparison is not in its asymptotic range.  That is why the certified
gradient is the analytic one.

## 6. The rigorous certificate — and why it does not close

`l3_certificate.py` and `interval_generators.py` build the certificate in exact
rational arithmetic, enclosing every analytic quantity by **interval evaluation
over each cell**, never by a hull of nodal values.  `exp` is bounded by a
truncated series with a rigorous geometric remainder; `χ` is enclosed by
monotonicity, which is sharp.

The structure is right and the independent checker passes.  **The bounds are
useless.**  Measured on the optimised shape at 25×49:

| bound | value | true value |
|---|---|---|
| viscous integral upper bound | 1.44e53 | 4.4e-2 |
| flux integral upper bound | 1.85e26 | — |
| divergence excursion | 1.45e2 | 0 |

Sixty orders of magnitude too loose.  The cause is the dependency problem
compounded near the edge of the bump support, where `χ''` carries `(1−σ)^{-4}`
and the interval evaluation cannot see the cancellation against `χ`'s
faster-than-polynomial decay.  The divergence enclosure straddling `±1.5e2`
around an exactly zero quantity is the clean diagnostic: the widening is in the
gradient components themselves, before any integration.

**Consequently no candidate is promoted.**  The rule was that `J > 0` must be
proved in intervals before promotion; it is not, so it is not.

What would be needed: adaptive subdivision concentrated at the support edge
(the widening is local, so uniform refinement is the wrong instrument), or a
generator basis whose derivatives do not carry `(1−σ)^{-4}` — a bump built from
a rational or a compactly supported spline rather than `exp(−1/(1−s))` would
lose `C^∞` but keep `C^k` for any fixed `k`, which is all the identity needs.
Both are separate experiments and neither is done.

## 7. Limitations

* `J = F'(0) > 0` is a statement about a single instant.  Escauriaza–Seregin–
  Šverák needs `limsup_{t→T}‖u(t)‖₃ = ∞`, and nothing here connects the two.
* The pressure channel is sign-indefinite for `t > 0`, with no monotonicity
  argument available.
* The large-`Re` regime where `P` wins is exactly where the numerics are least
  trustworthy.
* Even `J > 0` throughout a short interval gives bounded growth, not divergence.
* Viscosity is fixed and positive everywhere; no time-dependent or two-stage
  protocol appears.


## 8. Gate 9 update (2026-07-30)

Three changes recorded after this note was written.

* **The certificate basis changed.**  Section 6 diagnosed the flat-bump
  enclosures as structurally unfixable by subdivision; the basis is now
  Gaussian–Hermite (`gaussian_hermite.py`,
  [`gaussian_hermite_basis.md`](gaussian_hermite_basis.md)).  Measured effect:
  the divergence enclosure went from ±1.45e2 around an exact zero (gradient
  scale 5e-2) to ±0.38 on a 1/8-cell at gradient scale ~1, and the rigorous
  viscous lower bound from 47 orders of magnitude off to within a factor 7 of
  the true value.
* **The retried `J > 0` interval certificate still does not close.**  With the
  enclosure catastrophe gone, the remaining obstacle is the *product* of the
  discrete-pressure corner hull with the signed-flux enclosure: margin
  −2464 → −266 → −47 → −17 under refinement (9×17 → 33×65, sub=2), converging
  roughly linearly toward the barely-positive float value (+3.5e-3 at the
  favourable seed), i.e. about 2⁷ more linear refinement than is practical.
  The exterior Gaussian tails are *not* the obstacle (≤3.6e-6 viscous at box
  extent 4).  Recorded next step: evaluate the pressure term as
  `−3∫|u|u·∇p_h` with `∇p_h` an exact linear combination of nodal values,
  instead of `3∫p_h g` with the wide `g` enclosure.  No candidate is
  promoted; that remains the rule.
* **The spline rule became a guard.**  `require_clay_admissible` in
  `l3_optimizer.py` refuses to score or certify any family with
  `clay_admissible = False`; finite-C^k surrogates are tools (checker
  development, HS-5 prototyping, optimisation, mollification sources), never
  candidates.
