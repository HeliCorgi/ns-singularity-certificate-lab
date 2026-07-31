# The `n = 3` Kato constant in this repository's conventions

**Status (2026-07-31, turn 11).**  A self-contained derivation, from the
repository's own normalisation with no literature decimals anywhere in the
chain of inequalities, of a computable constant `G3` with

```
|  < P(v.grad w), w >_{H3dot}  |   <=   G3 ||v||_{H3dot} ||w||_{H3dot}^2
```

for mean-zero divergence-free real fields on `T^3`, together with the exact
per-trajectory band constants that replace the crude `9(K1+K2)` linear
coefficient of the turn-9 control inequality.  Certificate implementation:
`src/ns_certificate_lab/kato_constant.py`.  Literature comparison at the end
is *verification only* — no external number is an input to anything here.

## 1. Conventions (all of them)

| object | convention |
|---|---|
| torus | `T^3 = (R / 2 pi Z)^3`, coordinates `x` |
| measure | normalised: `d mu = (2 pi)^{-3} dx`; `mu(T^3) = 1` |
| Fourier basis | `e^{i k . x}`, `k` in `Z^3`; real fields: `u = sum_k uhat_k e^{ikx}`, `uhat_{-k} = conj(uhat_k)` |
| repo real coefficients | canonical `k`: `u = sum a_k cos(kx) + b_k sin(kx)`, so `uhat_k = (a_k - i b_k)/2` |
| Parseval | `<f, g>_{L^2(mu)} = sum_k fhat_k . conj(ghat_k)` (no `2 pi` factors) |
| `H^n dot` norm | `||u||_n^2 = sum_{k != 0} |k|^{2n} |uhat_k|^2` `=` repo `sobolev_sq(n)` `=` `sum_{canonical} |k|^{2n} (|a_k|^2+|b_k|^2)/2` |
| zero mode | absent (`uhat_0 = 0`); all sums skip `k = 0` |
| reality | `uhat_{-k} = conj(uhat_k)`, equivalent to the `cos`/`sin` representation |
| divergence-free | `k . uhat_k = 0` for every `k` |
| Leray multiplier | `P_k = I - k k^T / |k|^2`, orthogonal projection per mode, self-adjoint, `P w = w` for divergence-free `w` |
| convolution | `(fg)hat_k = sum_j fhat_j ghat_{k-j}` (exact; no measure factor in the normalised convention) |
| advection | `(v . grad w)hat_k = sum_j i (vhat_{k-j} . j) what_j` |
| lattice sums | `A_4 = sum_{p != 0} |p|^{-4}`, `A_6 = sum_{p != 0} |p|^{-6}`; rigorous rational upper bounds by exact partial sum + shell-count tail (`24 m^2 + 2 <= 26 m^2` points on the sup-shell `m`), as in `adot_squared_upper` |

Throughout, `p = k - j` is the `v`-frequency, `V_p = |p|^3 |vhat_p|`,
`W_k = |k|^3 |what_k|`, so `||V||_{l^2} = ||v||_3`, `||W||_{l^2} = ||w||_3`.

## 2. The symmetrised trilinear identity

For smooth mean-zero divergence-free real `v, w`, using `P w = w` and
self-adjointness of `P` (mode-wise, commuting with `|k|^{2n}`):

```
T := < P(v.grad w), w >_3 = < v.grad w, w >_3
   = sum_{k,j} |k|^6 i (vhat_{k-j} . j) (what_j . conj(what_k)) .
```

Taking the real part (T is real), using reality of `vhat`, `what` and the
divergence-freeness of `v` in the form `vhat_{k-j} . k = vhat_{k-j} . j`:

```
T = 1/2 sum_{k,j} i (vhat_{k-j} . j) ( |k|^6 - |j|^6 ) (what_j . conj(what_k)) .
```

At `n = 0` the factor `|k|^0 - |j|^0` vanishes: this is the classical
cancellation `<v.grad w, w>_{L^2} = 0`, and the identity above is its `H^3dot`
generalisation.  Everything below bounds this symmetrised sum.

## 3. Two dead ends, recorded

The natural target form `G3 <= sqrt( sup_k G(k) )` with `G(k)` a lattice sum
in the other frequency **does not converge in this kernel**, and the failure
is worth recording so nobody re-derives it:

- fix `k` and sum the squared kernel over `j`: for `|j| -> infinity` the
  summand behaves like `|j|/(2|k|^3)` squared — divergent;
- fix the `v`-frequency `p` and take `sup_k` of the kernel, then sum over
  `p`: the sup does *not* decay in `|p|` — in the regime `|k| = 1`,
  `|j| ~ |p|` the kernel grows like `|p|/2`, so `sum_p sup_k (...)^2`
  divergences as well.

The largeness always sits in whichever slot the arrangement refuses to
weight.  The convergent route is the classical one, made discrete and
explicit: a commutator estimate that uses the cancellation *first*.

## 4. The discrete Kato–Ponce commutator bound

Write `Lambda^3` for the multiplier `|k|^3`.  Then
`T = < [Lambda^3, v.grad] w, Lambda^3 w >_{L^2}` (the term
`< v.grad Lambda^3 w, Lambda^3 w >` vanishes by the `n = 0` cancellation),
so `|T| <= || [Lambda^3, v.grad] w ||_{L^2} ||w||_3`, and the commutator has
coefficients

```
( [Lambda^3, v.grad] w )hat_k = sum_j i ( |k|^3 - |j|^3 ) (vhat_{k-j} . j) what_j .
```

**Step 1 (algebra).**  `| |k|^3 - |j|^3 | = | |k| - |j| | (|k|^2 + |k||j| + |j|^2)
<= |p| . 3 max(|k|,|j|)^2 <= 3 |p| (|j| + |p|)^2`, and expanding
`(|j|+|p|)^2 = |j|^2 + 2|j||p| + |p|^2` gives three kernels; the middle one is
split by AM–GM, `|p|^2 |j|^2 <= (|p| |j|^3 + |p|^3 |j|)/2`.  Altogether the
coefficient is bounded by

```
3 (|p| |vhat_p|) (|j|^3 |what_j|)          (t1)
+ 3 (|p| |vhat_p|) (|j|^3 |what_j|)  + 3 (|p|^3 |vhat_p|) (|j| |what_j|)   (t2, after AM-GM)
+ 3 (|p|^3 |vhat_p|) (|j| |what_j|)        (t3)
```

summed over `j` (with `p = k - j`), i.e. by the convolution of an `l^1`
sequence with an `l^2` sequence in each term.

**Step 2 (Young `l^1 * l^2` + Cauchy–Schwarz).**  With
`sum_p |p| |vhat_p| <= sqrt(A_4) ||v||_3` and
`sum_j |j| |what_j| <= sqrt(A_4) ||w||_3` (Cauchy–Schwarz against the
`|.|^3`-weighted `l^2` norms; these are exactly the `adot` mechanics):

```
|| [Lambda^3, v.grad] w ||_{L^2}
   <= 6 sqrt(A_4) ||v||_3 ||w||_3     (t1 + first half of t2, l^1 in p)
    + 6 sqrt(A_4) ||v||_3 ||w||_3     (t3 + second half of t2, l^1 in j)
   = 12 sqrt(A_4) ||v||_3 ||w||_3 .
```

**Theorem (G3).**  In the conventions of section 1,

```
G3  <=  12 sqrt(A_4) ,
```

with `A_4 <= exact lattice sum over |p|_inf <= N  +  26/N` (rigorous rational
upper bound, monotone improving in `N`; the same object as
`adot_squared_upper`).  Computed by the certificate:
`N = 20`: `A_4 <= 17.32295`, `G3 <= 49.9451`;
`N = 40`: `A_4 <= 16.92444`, `G3 <= 49.3672`;
`N = 60`: `A_4 <= 16.79301`, `G3 <= 49.1752` — monotone as required.

The `12` is the price of the three crude steps (MVT bound on
`|k|^3 - |j|^3`, the `(|j|+|p|)^2` expansion, AM–GM); each is a finite
algebraic inequality (Lean-able) and each is a recorded sharpening lever.
The known-sharper machinery of Morosi–Pizzocchero (their `G_3` in their
normalisation) is a *comparison point*, section 8 — not an input.

## 5. Per-trajectory band constants (the actual chain linear coefficient)

For the chain, `v = u_a` is a finite trigonometric polynomial with interval
coefficients from the slab box, so both mixed terms of the difference
equation are bounded by **exact finite band sums**, with no universal
constant at all:

**(a) `< P(u_a . grad w), w >_3`** — rerun Step 1–2 with `vhat` supported on
the band and keep the `p`-sums exact:

```
| < P(u_a.grad w), w >_3 |  <=  C_kato(u_a) ||w||_3^2 ,
C_kato(u_a) = 6 sum_{p in spec(u_a)} ( |p| + |p|^3 ) |uhat_{a,p}|
            = 6 sum_{canonical p} ( |p| + |p|^3 ) sqrt(|a_p|^2 + |b_p|^2) .
```

(The last equality uses `|uhat_{-p}| = |uhat_p| = sqrt(|a|^2+|b|^2)/2` and
counts both signs.)

**(b) `< P(w . grad u_a), w >_3`** — no cancellation is available and none is
needed, because the advected factor is band-limited:

```
| < w.grad u_a, w >_3 |
  <= sum_{j in spec(u_a)} |j| |uhat_{a,j}| sum_k |k|^6 |what_{k-j}| |what_k|
  <= [ sum_{canonical j} |j| (1+|j|)^3 sqrt(|a_j|^2+|b_j|^2) ] ||w||_3^2
  =: C_shift(u_a) ||w||_3^2 ,
```

using `sum_k |k|^6 |what_{k-j}||what_k| <= ||w||_3 . sqrt( sum_m |m+j|^6 |what_m|^2 )`
and `sup_{m != 0} |m+j|^3/|m|^3 = (1+|j|)^3` (attained at `|m| = 1` parallel
to `j`; `|m| >= 1` on the mean-zero lattice).

## 6. The `n = 3` control inequality

For `w = u - u_a`, subtracting the equations and pairing with `w` in
`H^3dot` (the analytic regularity to run this is the business of the audited
EXT layer; every constant below is exact rational or a lattice certificate):

```
d+/dt R  <=  ( -nu + C_kato(u_a) + C_shift(u_a) ) R  +  G3 R^2  +  eps_3 ,
R = ||w||_{H3dot} ,   eps_3 = || (I - P_G) P(u_a.grad u_a) ||_{H3dot} ,
```

with `C_kato`, `C_shift` evaluated on the slab box (interval coefficients),
`G3 = 12 sqrt(A_4)` from the certificate, and the viscous term retained as
`-nu R` via the spectral gap `|k| >= 1`, exactly as at `n = 4`.  The `d+/dt`
display is shorthand: the analytic form actually delivered by the audited
EXT layer is the INTEGRAL energy inequality (EXT-P2-INT) together with the
comparison Lemma C of `ext_p1_p2_p3_audit.md` section 1.3 — no Dini
derivative is delivered or consumed (the pointwise-Dini upgrade is the named
open gap G-DINI, consumed by nothing).  The shape
`a R + b R^2 + e` is unchanged, so the entire `control_ode` machinery, the
Picard box, the Taylor endpoint and the recentring transfer (now measured in
`H^3dot`) are reused verbatim; `||w||_{H^3} <= sqrt8 ||w||_{H3dot}` converts
to the inhomogeneous norm as before.

What changed against turn 9's `n = 4` inequality: the linear coefficient
`9(K1+K2)` (an `l^1`-type bound, ~1.24e3 for P1) becomes
`C_kato + C_shift` (exact band sums, tens); the quadratic `135 Adot ~ 562`
becomes `G3 <= 48.8`.  The measured effect on the certified horizon is the
subject of the turn-11 chain rerun and is *not* asserted here.

## 7. The basic inequality `K_3` (same method)

Without the cancellation (bounding the full advection norm):
`|k|^3 <= 4(|p|^3 + |j|^3)` gives

```
|| v.grad w ||_3  <=  4 sqrt(A_4) ||v||_3 ||w||_3  +  4 sqrt(A_6) ||v||_3 ||w||_4 ,
```

and `||P(v.grad w)||_3 <= ||v.grad w||_3` since `P` contracts mode-wise.
`A_6 = sum |p|^{-6} <= exact sum + 9/N^3` (shell tail
`26 m^2 . m^{-6} = 26 m^{-4}`, `sum_{m>N} <= 26/(3N^3) <= 9/N^3`).  The
`||w||_4` factor is the expected one (`K_n` couples `n` to `n+1`); the chain
does not use `K_3` directly, it is recorded because the same certificate
machinery covers it at no cost.

## 8. Literature comparison (verification only)

Morosi–Pizzocchero work on `(0, 2pi)^3` with the *unnormalised* Lebesgue
measure; their norms satisfy `||u||_{MP,n}^2 = (2 pi)^3 ||u||_{ours,n}^2`
for the same Fourier data, so a trilinear constant converts as
`G_ours = (2 pi)^{3/2} G_MP ~ 15.75 G_MP`.  Their sharp `G_3 ~ 0.438`
therefore corresponds to `~ 6.9` in our normalisation, against our provable
`48.8`: our constant is about `7x` cruder, the price of the four elementary
steps in section 4, each individually sharpenable.  This comparison checks
the *order of magnitude* of the derivation and is used for nothing else; if
the converted numbers had disagreed wildly, that would have signalled a
normalisation error — they do not.

## 9. Lean-able finite parts

`|x^3 - y^3| <= 3 |x - y| max(x,y)^2` for `0 <= x, y`; the AM–GM split
`p^2 j^2 <= (p j^3 + p^3 j)/2`; `(j+p)^2 = j^2 + 2jp + p^2`;
`sup_{|m| >= 1} (1+|j|)^3`-ratio bound; the shell-count tail bounds for
`A_4`/`A_6`; and the checker's final assembly
`G3_upper = 12 * sqrt_upper(A4_upper)`.  The infinite-dimensional commutator
inequality itself is analysis over `l^2(Z^3)` and is **not** inserted into
Lean as an axiom; the Lean file `KatoConstant.lean` carries the finite
algebra and the audit note carries the rest.
