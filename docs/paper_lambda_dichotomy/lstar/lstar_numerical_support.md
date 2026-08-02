# Numerical support for (L\*) on the **smoothly** truncated coherent family

Producer: [`experiments/run_smooth_family_capacity.py`](../../../experiments/run_smooth_family_capacity.py)
· module: [`src/ns_certificate_lab/smooth_coherent_family.py`](../../../src/ns_certificate_lab/smooth_coherent_family.py)
· tests: [`tests/test_smooth_coherent_family.py`](../../../tests/test_smooth_coherent_family.py)
· data: [`outputs/lstar/smooth_family_capacity.json`](../../../outputs/lstar/smooth_family_capacity.json)

Run: 752.7 s wall (pre-registered cap 35 min, not binding).
Content digest `sha256:527b3d0b084cdfa2ac516fd0a8ecdd1c1fefb8f19466bc554e8f536a72e7d4b3`.

Companion analytic document: [`lstar_proof_main.md`](lstar_proof_main.md).
Compatibility: that document's Definition 0.1 asks for `chi in C^inf`; the
`chi` used here is `C^4` (and a rational `chi in C^p` is available for every
`p`, see section 0). Everything measured below is therefore evidence for the
`C^p` version of the statement; if the analytic route genuinely needs
`C^inf`, a non-polynomial `chi` must be used and the exact-rational lane
(tables A, B, D) is lost, while the float lanes (C, E, F) carry over
unchanged.

**No PDE theorem is asserted here.** Every row is labelled `exact`
(`fractions.Fraction`, no rounding) or `float` (IEEE double). Tables only;
the accompanying claims are labelled PROVEN / NUMERICAL / OPEN.

> **Status (2026-08-02).** The capacity bound these tables corroborate is
> now **PROVEN** analytically, with the sharp exponent \(N^3\), for the
> smoothly truncated family, uniformly in the admissible cutoff and in
> \(v_0\) ([`lstar_proof_main.md`](lstar_proof_main.md) Thm 7.1(3);
> adversarially refereed, [`lstar_referee.md`](lstar_referee.md), verdict
> (i)). Nothing below is load-bearing for that proof — in particular the
> (V-NONDEG) rows corroborate but do not carry Theorem 6.7, which is
> proven by a \(\chi\)-independent closed-form identity. The constants
> \(c_0,N_*\) of the proven bound are **non-effective**. The paper's
> literal Hypothesis (L\*), for the **sharply** truncated family, remains
> **OPEN** and is used by nothing; this note changes the family and does
> not prove it.

---

## 0. The family and the weight

| item | value |
|---|---|
| field | `u_hat_N(k) = chi(\|k\|/N) P_k v0 / \|k\|^2`, `1 <= \|k\|^2 <= N^2`, `P_k = I - k k^T/\|k\|^2` |
| seed | `v0 = (1,2,3)`, `\|v0\|^2 = 14` |
| weight | `chi(r) = 1 - S(s)`, `s = clamp((r^2 - 1/4)/(3/4), 0, 1)`, `S(s) = 126 s^5 - 420 s^6 + 540 s^7 - 315 s^8 + 70 s^9` |
| properties of `chi` | `chi = 1` on `[0,1/2]`; `chi = 0` on `[1,inf)`; `0 <= chi <= 1`; nonincreasing; `chi in C^4` |
| why this `chi` | `S` is the degree-9 smoothstep (normalised integral of `s^4(1-s)^4`); it is a **polynomial in `r^2`**, so `chi(\|k\|/N) in Q` for every lattice point and the whole family is exactly rational |
| convention | `u(x) = sum_k u_hat_k e^{ik.x}`, measure `(2pi)^{-3}dx`; `H0 = \|u\|_2^2`, `H1 = \|grad u\|_2^2`, `N0^2 = H1/H0`, `K = \|P(u.grad u)\|_2^2 / H1^2` — identical to `theorem_statement.md` and to `run_osgood_gate.py` |
| generalisation | replacing `S` by the degree-`(2p+1)` smoothstep gives a rational `chi in C^p` for any `p`; only a fixed finite `p` is used below |

---

## A. Exact moment laws — residual is the rational `0`, not "small"

Laws checked (Lemma 9 of `complete_proof.md` holds for **any** radial weight
on the sign-flip- and permutation-invariant band, so it applies verbatim to
`f(\|k\|^2) = chi(\|k\|/N)^2 \|k\|^{-4}` etc.):

```
H0     = (2/3)|v0|^2 sum_k chi(|k|/N)^2 / |k|^4
H1     = (2/3)|v0|^2 sum_k chi(|k|/N)^2 / |k|^2
u_N(0) = (2/3) v0    sum_k chi(|k|/N)   / |k|^2
```

| N | lattice points | active (`chi != 0`) | `H0` residual | `H1` residual | `u_N(0)` residual | real cosine field | `k . u_hat = 0` | `N0^2` (float of exact) |
|---:|---:|---:|:--|:--|:--|:--|:--|---:|
| 4 | 256 | 250 | `0` | `0` | `0,0,0` | yes | yes | 2.316030 |
| 6 | 924 | 894 | `0` | `0` | `0,0,0` | yes | yes | 3.455252 |
| 8 | 2 108 | 2 102 | `0` | `0` | `0,0,0` | yes | yes | 4.585959 |
| 10 | 4 168 | 4 138 | `0` | `0` | `0,0,0` | yes | yes | 5.718747 |
| 12 | 7 152 | 7 122 | `0` | `0` | `0,0,0` | yes | yes | 6.852077 |

`u_N(0)` is exactly proportional to `v0` (components in ratio `1 : 2 : 3`),
as the law demands. Lane label: **exact**.

---

## B. Exact capacity, full `O(|B_N|^2)` contraction

No truncation, no sampling: all ordered pairs of active band points.
`K` and `||P(u.grad u)||_2^2` are stored as full rationals in the JSON
(519 decimal digits for `K` at `N = 8`; the transcript abbreviates larger
ones by bit length without rounding).

| N | `\|\|P(u.grad u)\|\|_2^2 / N^3` | `K` | `N0^2` | `K / N0^2` | label |
|---:|---:|---:|---:|---:|:--|
| 4 | 678.284307 | 0.61657681 | 2.316030 | 0.266221 | exact |
| 6 | 1027.307364 | 1.13589428 | 3.455252 | 0.328744 | exact |
| 8 | 1241.509295 | 1.67053905 | 4.585959 | 0.364273 | exact |
| 10 | 1384.356255 | 2.20774755 | 5.718747 | 0.386054 | exact |
| 12 | 1485.779047 | 2.74636641 | 6.852077 | 0.400808 | exact |
| 16 | 1619.652353 | 3.82541413 | 9.119441 | 0.419479 | exact |

The sharp family is exactly certified only to `N = 8` in the paper's
appendix; the smooth family reaches `N = 16` because its weight is rational.

---

## C. Dealiased FFT continuation, and the exact cross-check

Grid `= 4N + 2`, the smallest even grid clearing the repository's dealiasing
guard (`u` lives in `|k|_inf <= N`, the product in `|k|_inf <= 2N`, and
`leray_advection` demands `2N < floor(G/2)`). The **whole** product spectrum
is represented — nothing is truncated, which is what a norm requires.

| N | grid | `\|\|P(u.grad u)\|\|_2^2 / N^3` | `K` | `N0^2` | `K / N0^2` | `H1 / N` | label |
|---:|---:|---:|---:|---:|---:|---:|:--|
| 8 | 34 | 1241.509295 | 1.67053905 | 4.585959 | 0.364273 | 77.106620 | float |
| 16 | 66 | 1619.652353 | 3.82541413 | 9.119441 | 0.419479 | 82.306032 | float |
| 24 | 98 | 1761.491502 | 5.98587403 | 13.655282 | 0.438356 | 84.039239 | float |
| 32 | 130 | 1835.406669 | 8.14717671 | 18.191651 | 0.447853 | 84.905842 | float |
| 40 | 162 | 1880.711599 | 10.30870339 | 22.728223 | 0.453564 | 85.425804 | float |
| 48 | 194 | 1911.311477 | 12.47030632 | 27.264894 | 0.457376 | 85.772445 | float |

**Exact vs float** (lane B against lane C at the shared bands):

| N | relative difference in `\|\|P(u.grad u)\|\|_2^2` | relative difference in `K` |
|---:|---:|---:|
| 8 | 5.5e-16 | 5.3e-16 |
| 16 | 7.0e-16 | 8.1e-16 |

Agreement at the level of double-precision roundoff, at a band where the
exact contraction enumerates `2102^2` and `17176^2` ordered pairs
respectively.

**Two-point Richardson limits** (fit `L - c/N` on `N = 40, 48`; a
diagnostic, not a proof):

| quantity | limit |
|---|---:|
| `\|\|P(u_N.grad u_N)\|\|_2^2 / N^3` | **2064.31** |
| `K / N0^2` | **0.47644** |
| `H1 / N` | 87.5057 |
| `N0^2 / N` | 0.56708 |
| sharp-cutoff reference `K/N0^2` (from `complete_proof.md`) | 0.259 (`N=4`) -> 0.396 (`N=32`) |

Reading: both target sequences increase monotonically toward **positive**
constants. `K/N0^2 -> 0.476` means `K >~ 0.476 e^{z}` with `z = log N0^2`,
i.e. the **sharp** exponent, not merely some `a > 0`. The smooth family's
`K/N0^2` runs slightly above the sharp family's at every comparable band.

---

## D. Exact single Fourier coefficients of the nonlinearity

`u_hat_N(k) = N^{-2} F(k/N)` **exactly** (because `P_k = P_{k/N}`), where
`F(xi) = chi(|xi|) P_xi v0 / |xi|^2`. The discrete convolution
`b_k = sum_{l+m=k} (u_hat(l).m) u_hat(m)` is therefore *literally* a Riemann
sum, with spacing `1/N`, for

```
C(zeta) = int_{R^3} (F(zeta-b).b) F(b) db  =  -i * FT[V.grad V](zeta),
```

so `b_k -> C(k/N)` with **no `N`-power normalisation**. `|P_k b_k|` below is
exact rational arithmetic throughout; `zeta_outer` is deliberately outside the
band (`|zeta| = 1.118`), where the sharp-cutoff picture says nothing.

| N | `zeta = (1/2,1/2,1/2)` | `zeta = (1/4,1/2,1/4)` | `zeta = (1,1/2,0)` |
|---:|---:|---:|---:|
| 8 | 6.6569313637 | 27.0252497212 | 1.1554654517 |
| 16 | 6.9250757721 | 30.5573481635 | 1.1561889652 |
| 24 | 7.0221271844 | 31.7466591193 | 1.1561903251 |
| 32 | 7.0717405781 | 32.3426883926 | 1.1561904559 |
| 40 | 7.1017732957 | 32.7006261834 | 1.1561904798 |
| 48 | 7.1218836284 | 32.9393567138 | 1.1561904817 |
| Richardson `L - c/N` | **7.2224352916** | **34.1330093663** | **1.1561904909** |

All entries are exactly nonzero as rationals (`is_exactly_zero = false`), and
each column increases monotonically. Convergence is `O(1/N)` for the two
inner frequencies and far faster for `zeta_outer` — there `|zeta| > 1` makes
`F(b)` vanish at the far singular point, leaving only one integrable
singularity in the convolution. Lane label: **exact**.

---

## E. The same coefficient by direct continuum quadrature

`C(zeta)` integrated in 3-D. The two integrable `|.|^{-2}` singularities (at
`b = 0` and `b = zeta`) are separated by the smooth partition
`phi = d_zeta^4/(d_0^4 + d_zeta^4)`, and each half is then taken in spherical
coordinates centred on **its own** singularity, where the `rho^2` Jacobian
cancels the singular factor exactly. Node counts `(radial, polar, azimuth)`:

| orders | increment from previous (`zeta = (1/2,1/2,1/2)`) |
|---|---:|
| 20,20,20 | — |
| 30,30,30 | 1.04e-04 |
| 40,40,40 | 9.70e-07 |
| 60,60,60 | 5.39e-08 |
| 80,80,80 | 6.18e-09 |

| `zeta` | `\|P_zeta C(zeta)\|` (float, quadrature) | discrete Richardson limit (from exact lane D) | relative difference |
|---|---:|---:|---:|
| `(1/2,1/2,1/2)` | 7.2230770878 | 7.2224352916 | 8.9e-05 |
| `(1/4,1/2,1/4)` | 34.1337627827 | 34.1330093663 | 2.2e-05 |
| `(1,1/2,0)` | 1.1561903241 | 1.1561904909 | 1.4e-07 |

Two completely disjoint code paths — exact rational lattice sums, and
floating-point 3-D quadrature of a continuum integral — agree to 1e-4 or
better on a quantity that would be `0` if `V` were a stationary Euler field.

---

## F. The continuum profile `V` and the decisive fact (V-NONDEG)

Radial reduction (all integrands smooth on `[0,1]`; the `|xi|^{-2}`
singularity is removed by the spherical Jacobian):

```
V(y) = A(r) v0 + B(r) yhat (yhat.v0),     A = (g+h)/2,  B = (g-3h)/2
g(r) = 4pi int_0^1 chi(rho) j0(rho r) drho
h(r) = 4pi int_0^1 chi(rho) B(rho r) drho,     B(s) = int_0^1 mu^2 cos(s mu) dmu
```

In the frame with polar angle `theta` measured from `v0`, with
`alpha = A + B = g - h`:

```
V     = alpha cos(theta) e_r  -  A sin(theta) e_theta
omega = W(r) sin(theta) e_phi,                 W = (alpha - (rA)')/r
curl(V.grad V) = Z(r) sin(theta) cos(theta) e_phi,   Z = ((rW alpha)' - 2WA)/r
```

with an overall `||v0||^2`. Untruncated closed form (`chi == 1`):
`g = 2pi^2/r`, `h == 0`, hence `V_inf = pi^2 (v0 + yhat(yhat.v0))/|y|` and
`Z_inf = -12 pi^4 ||v0||^2 / r^4`.

| r | `Z(r)·\|\|v0\|\|^2` | `-12 pi^4 \|\|v0\|\|^2 / r^4` | ratio | `r A(r)` | `r B(r)` | solenoidality residual `r alpha' + 2alpha - 2A` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | -3.158e-01 | -4.189e+06 | 0.000000 | | | -3.6e-15 |
| 0.50 | -1.250e+00 | -2.618e+05 | 0.000005 | | | -1.1e-14 |
| 1.00 | -4.796e+00 | -1.636e+04 | 0.000293 | 6.3217 | 0.1358 | -7.1e-15 |
| 2.00 | -1.624e+01 | -1.023e+03 | 0.015877 | 11.1351 | 0.9984 | 0.0e+00 |
| 4.00 | **-3.332e+01** | -6.392e+01 | 0.521166 | 13.5094 | 5.6797 | -1.8e-15 |
| 8.00 | -8.172e+00 | -3.995e+00 | 2.045374 | 8.2118 | 11.6421 | -1.3e-15 |
| 16.00 | -5.662e-01 | -2.497e-01 | 2.267395 | 9.3837 | 10.3719 | -2.0e-15 |
| 32.00 | -1.622e-02 | -1.561e-02 | 1.039593 | 9.8587 | 9.8797 | 4.4e-16 |
| 48.00 | -2.949e-03 | -3.083e-03 | 0.956539 | | | 9.4e-16 |
| 64.00 | -9.915e-04 | -9.754e-04 | 1.016516 | 9.8693 | 9.8699 | 1.1e-15 |
| 96.00 | -1.932e-04 | -1.927e-04 | 1.002815 | 9.8693 | 9.8699 | 1.7e-16 |

`pi^2 = 9.869604`. `r A` and `r B` both converge to it, confirming the
claimed `V = V_inf + o(1/|y|)` tail; the `ratio` column oscillates around `1`
with decaying amplitude, confirming `Z ~ -12 pi^4 ||v0||^2 / r^4`.

**Honest note on sign.** `Z` is *not* one-signed: the compactly supported
`chi` produces oscillatory corrections and `Z` changes sign near
`r ~ 10.2, 12.2, 18.2, 20.1`. Nothing in the argument needs a sign — only
`Z != 0` at **one** radius, and `|Z(4)| = 33.3` with a quadrature error near
machine epsilon.

**Independent implementation.** `V` is rebuilt in Cartesian coordinates from
`A` and `B` with the seed along `e_1` (a frame sharing nothing with the
derivation above), `V.grad V` and its curl taken by central differences:

| point | relative difference from the axisymmetric reduction |
|---|---:|
| `(2, 1, 0.5)` | 9.8e-07 |
| `(4, -2, 1)` | 4.0e-07 |

(at the `1e-4` finite-difference truncation floor).

---

## G. Status ledger

| step | status | note |
|---|---|---|
| `K = \|\|N\|\|^2/H1^2`, `H1 ≍ N`, `N0^2 ≍ N`; so `\|\|N\|\|^2 >= cN^{2+a}` gives `K >= c' e^{as}` and `int ds/Phi < inf` | **PROVEN** | re-derived independently. Table C gives `H1/N -> 87.51` and `N0^2/N -> 0.567`, so along the dyadic sequence `N_j = 2^j N_0` the gaps `s_{j+1} - s_j` tend to `log 2` -- bounded above and below, which is exactly what the Proposition's monotonicity step needs; no crude Lemma-11-style lattice bound is required |
| moment laws survive the smooth radial weight | **PROVEN** + exactly verified | Lemma 9 needs only that the band is sign-flip/permutation invariant and the weight radial. Table A: every residual is the rational `0` |
| `u_hat_N(k) = N^{-2} F(k/N)` | **PROVEN** (elementary; `P_k = P_{k/N}`) | this, not Poisson summation, is what the argument needs |
| Poisson-summation periodization `u_N(x) = N sum_m V(N(x+2pi m))` with `m != 0` terms `O(1/N)` uniformly on `\|x\| <= pi` | **FALSE**, in the naive form and in every regularised form | `V(y) ~ pi^2(v0 + yhat(yhat.v0))/\|y\|` decays only like `\|y\|^{-1}` with a sign-definite spherical mean in the `v0` direction, and `sum_{m in Z^3} \|m\|^{-1} = inf`: the periodization diverges. Because the divergent terms are *eventually positive*, no summation method (Abel included) rescues it — see `lstar_proof_main.md` Remark 3.4. The route avoids the identity entirely (next row); nothing here or there depends on it |
| duality lower bound, Fourier-side | **PROVEN** | with `psi_N(x) = N^{3/2} Psi(Nx)`, `Psi` divergence-free `C_c^inf`: `psi_hat_N(k) = N^{-3/2} tilde-Psi(k/N)`, `\|psi_N\|_{L^2(T^3)}^2 = (2pi)^{-3}\|Psi\|_{L^2(R^3)}^2`, and the pairing is `N^{3/2}(2pi)^{-3}<V.grad V, Psi> (1+o(1))` by a Riemann sum. Hence `\|P(u_N.grad u_N)\|_2^2 >= c N^3`: **exponent `a = 1`, the sharp one** |
| Riemann-sum convergence of that pairing | **PROVEN modulo standard** | the integrand `C(zeta) conj(tilde-Psi(zeta))` is bounded, compactly supported in `\|zeta\| <= 2`, and continuous off `{0}` (`C` is degree-0 homogeneous near `0`), hence Riemann integrable; the lattice sum omits only the singular cells, contributing `O(N^{-2})` |
| `V_inf = pi^2(v0 + yhat(yhat.v0))/\|y\|` | **PROVEN modulo standard** | needs only `FT[\|xi\|^{-2}] = 2pi^2/\|y\|`; the tensor `T_ij` is fixed by its trace (`= g`) and its `yhat yhat` contraction (`= 0`, from `int_0^inf [sin s/s + 2cos s/s^2 - 2 sin s/s^3] ds = 0`). Table F confirms `rA, rB -> pi^2` |
| `V = V_inf + O(\|y\|^{-p})`, `p` large | **PROVEN modulo standard** | `(1-chi)F_inf` is a `C^p` symbol of order `-2` vanishing near `0`; its transform is `O(\|y\|^{-p})` away from the origin. `chi in C^p` for any `p` is available inside the rational family |
| **(V-NONDEG): `P(V.grad V) != 0`** | **PROVEN** (given the two rows above) + verified 3 independent ways | `ker P` = curl-free, and curl is **local**, so one point suffices. `curl(V.grad V) = Z(r) sin cos e_phi` with `Z -> -12 pi^4\|\|v0\|\|^2/r^4 != 0`. Numerics: lane F (`\|Z(4)\| = 33.3`), lane E (`\|P_zeta C(zeta)\| = 7.223`), lane D (exact rational, same number) |
| an `N`-uniform lower bound proved *from the exact lattice lane alone* | **OPEN** | lane D certifies exact nonzero values at finitely many `N`; upgrading it to a proof needs a quantitative Riemann-sum error bound, which is not attempted here. The analytic route above does not need it |
| (L\*) for the **sharp** family | **OPEN, and untouched** | this note changes the family; it does not prove the paper's stated (L\*) |
