# The Gaussian–Hermite generator basis

**Status: a design decision with its reason attached.**  It records why the
certificate basis for the axisymmetric generator families changed from
flat-bump cutoffs to Gaussian–Hermite generators, what properties the new basis
must have and how each is met, and the standing rules for the two bases it
displaces or demotes.  Nothing here is a result about Navier–Stokes; it is
about making interval enclosures of initial-datum functionals tight enough to
be worth computing.

## 1. Why the basis changed

The flat-bump construction evaluates every generator through
`χ(σ) = e^{−1/(1−σ)}`, and its second derivative is

```
χ'' = χ ( (1−σ)^{−4} − 2(1−σ)^{−3} ) ,
```

so every gradient enclosure near the support edge multiplies an interval
containing a fourth-order pole by an interval containing an essential zero.
Interval arithmetic cannot see that cancellation — the two factors are treated
as independent — and the result, measured in
[`l3_positive_generation_search.md`](l3_positive_generation_search.md) §6, was
a certificate whose structure was right and whose bounds were useless: a
viscous integral truly worth `4.4e-2` enclosed below `1.44e53`, an exactly-zero
divergence enclosed in `±1.45e2`.  Roughly sixty orders of magnitude of
dependency widening, concentrated where `χ''` carries `(1−σ)^{−4}`.  By the
promotion rule recorded there, no candidate was promoted, and by the diagnosis
recorded there, the fix is not subdivision alone — the widening is structural
to the basis.

## 2. The new generators, and why they cannot repeat the failure

The Gaussian–Hermite generators are

```
u₁, ψ₁  of the form  P(r², z) · e^{−α r² − β z²} ,
```

with `P` a polynomial and `α, β > 0`.  The family is **closed under
differentiation with polynomial coefficient recursion**: writing `s = r²`,

```
∂_s [ Q e^{−αs−βz²} ] = ( ∂_s Q − α Q ) e^{−αs−βz²} ,
∂_z [ Q e^{−αs−βz²} ] = ( ∂_z Q − 2βz Q ) e^{−αs−βz²} ,
```

so every derivative of every order is (a polynomial with explicitly computed
coefficients) times the *same* exponential — the pure-Gaussian case generates
the Hermite polynomials, whence the name.  Interval evaluation of any
derivative on a box is therefore **polynomial interval arithmetic plus one
monotone `exp` call per box**: the exponent `−αs − βz²` is linear in `s` and an
exact square in `z` (the exact-range square of `control_ode.py`'s `_square`,
not a self-product), its range on a box is exact, and `exp` is monotone, so the
single transcendental evaluation is sharp end to end via the existing
`exp_interval` truncated-series enclosure.  There is no pole, no
`(1−σ)^{−4}`-against-essential-zero cancellation for the arithmetic to miss,
and hence no dependency catastrophe: the only widening left is the ordinary
polynomial dependency, which subdivision does shrink linearly.

## 3. The required properties, one by one

| property required of a certificate basis | how the Gaussian–Hermite basis meets it |
|---|---|
| Cartesian `C^∞` | generators are functions of `(r², z)`, so the assembled field is smooth across the axis by the same argument as the flat-bump construction (`mixed_initial_data.py`); here the generators are moreover entire |
| divergence-free | the same algebraic identity as the flat-bump construction: `u^θ = r u₁`, `u^r = −r ∂_z ψ₁`, `u^z = 2ψ₁ + 2s ∂_s ψ₁`, an identity for **any** `ψ₁` |
| rapid decay, finite energy, finite `L³` | the Gaussian factor dominates every polynomial, so every derivative decays faster than any power and every `L^p` and Sobolev integral converges |
| closed-form derivatives | the polynomial recursion of §2, exact at every order |
| analytic tail bounds | in `z`, the Mills ratio `∫_x^∞ e^{−βt²} dt ≤ e^{−βx²}/(2βx)` plus the by-parts moment recursion `∫_x^∞ t^m e^{−βt²} dt = x^{m−1}e^{−βx²}/(2β) + ((m−1)/(2β)) ∫_x^∞ t^{m−2} e^{−βt²} dt`; in `s = r²`, the radial measure `r dr = ½ ds` makes tails **exactly** `∫_S^∞ s^m e^{−cs} ds = (m!/c^{m+1}) e^{−cS} Σ_{i≤m} (cS)^i/i!` — closed form, one exponential |
| rational parameterisation | `α`, `β` and the coefficients of `P` are rational, and the recursion of §2 preserves rationality at every derivative order; the single transcendental per box is `exp`, enclosed by the audited series routine |
| independent pressure-source computation | unchanged: the existing 3D solver (`AxisymmetricPressureSolver` and the free-space Poisson machinery) consumes field values and is indifferent to the basis that produced them |

Truncation is taken at the standard evaluation box already used by the Gate-8
grids, with everything outside it bounded analytically by the tail rows above —
the tail is a bound that gets *added*, not an error that gets ignored.

## 4. The standing of the old bases

**Flat-bump families remain, for regression only.**  They stay in the
repository so that the Gaussian–Hermite results can be compared against the
Gate-8 numbers on identical footing, and because their float-level evaluation
is perfectly sound — only their *interval* evaluation is poisoned.  They are no
longer the certificate basis.  If a flat-bump family is ever to be certified
again, the natural-interval route through `χ`, `χ'`, `χ''` is **forbidden** —
§1 is the standing evidence — and the recorded plan is a Taylor-model route:
substitute `y = 1/(1−s)` and evaluate `y^m e^{−y}` **as a single unit** (a
bounded function with its maximum at `y = m` and monotone on either side), so
that the cancellation the natural extension cannot see is performed inside one
enclosure rather than between two.  That route is a plan, not an
implementation, and nothing may cite it as done.

**The `C^k` spline rule, restated.**  Finite-`C^k` splines are surrogates:
permitted for checker development, HS-5 prototyping, optimisation, and as
mollification sources — never Clay candidates, because the Clay statements
require `C^∞` data and a finite-`C^k` object is not in the class.  The rule is
enforced mechanically: `require_clay_admissible` in
`src/ns_certificate_lab/l3_optimizer.py` refuses to score or certify any family
flagged `clay_admissible = False`, and nothing downstream of that guard may
describe a surrogate as a candidate.  The Gaussian–Hermite basis needs no such
flag — its members are genuinely `C^∞` with genuine decay — which is precisely
why it can be the certificate basis and a spline cannot.
