# Continuum-to-lattice shadowing scaffolding for the doubling Poincare map

**Verification sprint v1, workstream E — 2026-08-02.**
**Label: BINARY64 VERIFICATION DIAGNOSTIC / NOT A PROOF.**
No PDE theorem, no singularity, no Clay statement. Every number below is a
binary64 measurement produced by
[`experiments/run_poincare_shadowing.py`](../../../experiments/run_poincare_shadowing.py);
raw records with digests live in
`outputs/verification_sprint_v1/shadowing/summary.json`; the tabulated rows are
[`poincare_shadowing_scaling.csv`](poincare_shadowing_scaling.csv).
Nothing in this note is exact-rational; there is no interval enclosure.

## 0. What is being audited, and what is already dead

Pilot A of the front-flow portfolio
([`discovery_2026_08_01_front_flow_portfolio.md` §5](../discovery_2026_08_01_front_flow_portfolio.md))
forward-integrated the lattice doubling map over
`c_E in {1,30,100,228,500} x {front-only, keep-sea}` at `N_0 = 4`, `W = 2` and
found **no positive-gain plateau and no recurrent shape anywhere**: nine of ten
runs collapsed into the binary64 noise floor within 6–11 stages, and the
strongest run finished *decaying* at gain `9.4e-5`.

Consequently **there is no periodic orbit in the scanned box, hence no orbit,
no monodromy, and no Floquet radius**. An orbit-shadowing theorem needs all
three. This workstream therefore does not attempt to shadow anything; it
builds and measures the *scaffolding* such a theorem would require, and issues
an honest verdict on the "orbit shadowing proof candidate".

The scaffolding has three measurable pieces:

| piece | question | section |
|---|---|---|
| `R_N`, `I_N` | do lattice and continuum data map into each other cleanly, and at what order? | §1 |
| consistency of `Q` | how far is the lattice Leray convolution from the continuum `Q` at resolution `N`? | §2 |
| consistency and contraction of `P` | how far is the one-stage lattice Poincare map from its `2N` refinement, and is `P` a contraction? | §3, §4 |

The pre-registered acceptance rule was fixed before any number was produced
(§5).

## 1. The restriction and interpolation operators

### 1.1 Data model

Repository conventions throughout: `u(x) = sum_k u_hat(k) e^{i k.x}` on `T^3`
with `||u||_2^2 = sum_k |u_hat(k)|^2`;
`N_k = -i P_k sum_{l+m=k} (m . u_l) u_m` with `P_k = I - k (x) k / |k|^2`;
critical energy law `E_N = c_E / N`.

A **cell profile** is a finite carrier set `A subset Z^3` together with
maps `U_a : Omega_eta -> C^3`, `Omega_eta = [-eta, eta]^3`, with
`U_{-a}(-xi) = conj(U_a(xi))` (reality) and the **continuum transversality**
condition

```
(a + xi) . U_a(xi) = 0            for all xi in Omega_eta.          (T)
```

Cells are assumed separated: `2 eta N < min_{a != b} |N a - N b|_inf`, so every
`k` lies in at most one cell and the decomposition `k = N a + q`,
`q/N in Omega_eta`, is unique.

### 1.2 Restriction `R_N` (continuum -> lattice)

```
(R_N U)(k) := N^{-2} U_a(q/N)     if k = N a + q with a in A, q/N in Omega_eta,
(R_N U)(k) := 0                   otherwise.
```

This is exactly the front-flow ansatz `u_hat(k) = N^{-2} Psi(k/N)` written
cell-wise, and it is the amplitude convention under which the repository's own
`E_N = c_E/N` normalization is `N`-independent in profile variables:

```
sum_k |u_hat(k)|^2 = N^{-4} sum_q |U(q/N)|^2  ~  N^{-4} . N^3 int |U|^2
                   = N^{-1} int |U|^2  =  2 c_E / N   <=>   int |U|^2 = 2 c_E.
```

Two exact structural facts, both used below:

* **`R_N` preserves reality exactly**, by `U_{-a}(-xi) = conj(U_a(xi))`.
* **`R_N` preserves divergence-freeness exactly, for every `N`**: with
  `k = N(a + xi)`, `k . (R_N U)(k) = N^{-1} (a + xi) . U_a(xi) = 0` by (T).
  Transversality is a continuum condition and survives sampling untouched.

### 1.3 Interpolation `I_N` (lattice -> continuum)

Let `w(y) = prod_{i=1}^{3} max(0, 1 - |y_i|)` be the tensor-product hat kernel.

```
(I_N u)_a(xi) := sum_{q in Z^3} w(N xi - q) . N^2 u_hat(N a + q),   xi in Omega_eta,
```

i.e. on each lattice cube `q/N + [0, 1/N]^3` the profile is the trilinear
interpolant of the eight corner values `N^2 u_hat(N a + q')`.

### 1.4 `R_N I_N = id` on lattice data

`w` is nodally exact: `w(N xi - q)|_{xi = q'/N} = delta_{q q'}`. Hence
`(I_N u)_a(q/N) = N^2 u_hat(N a + q)` and

```
(R_N I_N u)(N a + q) = N^{-2} . N^2 u_hat(N a + q) = u_hat(N a + q),
```

with both sides zero off the cells. So **`R_N I_N = id`** on lattice data
supported in the cells: `R_N` is onto that subspace, `I_N` is a right inverse,
and `I_N R_N` is the projection onto continuous piecewise-trilinear profiles.

### 1.5 `I_N R_N - id` on `C^1` profiles

On a cube `C_q = q/N + [0,1/N]^3` the trilinear interpolant is a *convex
combination* of the eight corner values of `U`, so for every `xi in C_q`

```
| (I_N R_N U)(xi) - U(xi) |  <=  max_{corners c} | U(c) - U(xi) |
                             <=  diam(C_q) . ||grad U||_inf
                             =   sqrt(3) N^{-1} ||grad U||_inf .
```

Hence, uniformly over the cell,

```
|| (I_N R_N - id) U ||_{L^inf(Omega_eta)}  <=  sqrt(3) ||grad U||_inf / N,
|| (I_N R_N - id) U ||_{L^2(Omega_eta)}    <=  sqrt(3) |Omega_eta|^{1/2} ||grad U||_inf / N.
```

**First order in `1/N` on `C^1`.** The bound is attained in order for merely
Lipschitz profiles (trilinear interpolation reproduces affine data exactly, so
the leading term is the modulus of continuity of `grad U`). For `U in C^2` the
classical tensor-product bound upgrades this to
`(1/8) N^{-2} sum_i ||d^2_i U||_inf` — second order — and the fixed-relative
two-box family of §2 is smooth except at the Fejer kink, which is why the
measured exponent there lands between 1 and 2.

### 1.6 The one obstruction `I_N` carries: it breaks the Leray constraint

`R_N` preserves (T) exactly; `I_N` does **not**. Since `(a + .) . U ≡ 0`, the
interpolated divergence is a pure commutator:

```
(a + xi) . (I_N R_N U)(xi)
  = sum_q w(N xi - q) [ (a + xi) - (a + q/N) ] . U(q/N)
  = sum_q w(N xi - q) (xi - q/N) . U(q/N),
```

so

```
| (a + xi) . (I_N R_N U)(xi) |  <=  sqrt(3) ||U||_inf / N .            (D)
```

The reconstructed continuum profile therefore carries a divergence defect of
size `O(N^{-1})` relative to `||U||_inf` — first order, *not* exact. A
shadowing argument that reconstructs a continuum profile from lattice data
must either carry (D) as an explicit error term or re-project after `I_N`
(re-projection is available exactly, since `P_{a+xi}` is a bounded pointwise
projector, but it changes the profile at the same `O(N^{-1})` order). This is
recorded as a scaffolding obligation, not a defect of the measurements below,
which use `R_N` on both sides and never invoke `I_N`.

## 2. Consistency of the lattice Leray convolution against continuum `Q`

### 2.1 What the comparison measures

With `u_hat(k) = N^{-2} Psi(k/N)`, the repo's lattice nonlinearity at `k = N xi`
is

```
N_k = -i P_k sum_l (m . u_l) u_m
    = -i P_xi N^{-3} sum_l (mu . Psi(lambda)) Psi(mu)
    ~ -i P_xi int ((xi - lambda) . Psi(lambda)) Psi(xi - lambda) d^3 lambda
    = -Q(Psi, Psi)(xi),
```

the Riemann factor `N^3` cancelling the amplitude factor `N^{-3}` exactly
(audit item S1, confirmed). The lattice coefficient is therefore an `N^0`
approximation of `-Q(xi)` with **no rescaling**, so blocks computed at `N` and
at `2N` are directly comparable — and this is confirmed empirically: the
child-band centre coefficient norm is
`0.5773503, 0.5764007, 0.5763499, 0.5763388` at `N = 8, 16, 32, 64`.

The doubling pairing is exactly `2k <-> k`: a coarse child target
`k = N(p+q) + o` and the fine target `k' = 2N(p+q) + 2o = 2k` correspond to the
*same* `xi = (p+q) + o/N`. Comparing `Q_N` at `k` with `Q_{2N}` at `2k` is
therefore comparison of the two lattice quadratures of one continuum integral
at a common set of continuum points, i.e. `R_N` vs `R_{2N}` composed with the
doubling pairing — precisely the `R/I` comparison asked for.

### 2.2 Family and protocol

Fixed-relative two-box cloud, `p = (1,1,0)`, `q = (1,0,1)`, `eta = 3/16`,
`width_override = floor(3N/16)`, Fejer weight, modewise Leray polarization,
`c_E = 1`, `nu = 1/40`. Backend: the exact zero-padded local-FFT convolution
`ns_certificate_lab.mesoscopic_local_fft` (every coefficient of
`P((u.grad)u)`, linear convolution, no aliasing, binary64 coefficient
arithmetic).

Two corrections were necessary to make this a genuine continuum-limit test,
and both are load-bearing:

1. **The phase slope must act on `xi`, not on the integer offset.** The
   module's default `phase_slope` multiplies integer offsets, so the profile's
   phase gradient in `xi` grows like `N` and the underlying continuum profile
   is *not* fixed. The runs below use
   `phase_slope = CONTINUUM_PHASE_SLOPE / N` with
   `CONTINUUM_PHASE_SLOPE = 16 x (0.173, -0.119, 0.071)`, so that the `N = 16`
   member reproduces the repository default exactly. The unrescaled family is
   run as an explicit **null control**.
2. **`floor(3N/16)` realises `eta = 3/16` exactly only when `16 | N`.** For
   `N in {8, 24}` the coarse and fine members have *different* relative widths
   (`0.1250` vs `0.1875`, and `0.1667` vs `0.1875`), so those pairs compare two
   different continuum profiles and cannot measure a discretization rate. They
   are **measured, recorded, and excluded from the fit** — never deleted.

For `16 | N` the family is an exact restriction: the Fejer weight
`prod_i (1 - |o_i|/W)` equals the fixed continuum tent `prod_i (1 - |xi_i|/eta)`
sampled at `xi = o/N`, and it vanishes precisely at `|o_i| = W`, which the box
`|o|_inf <= W-1` excludes. So `R_N U` and `R_{2N} U` are *exact* restrictions of
one and the same continuum tent profile and the only residual error is the
lattice quadrature of `Q`.

### 2.3 Results

`err(N) = || Q_N - Q_{2N} ||_2 / || Q_{2N} ||_2` on the common child band
(`Q_N` at `k = N(p+q)+o`, `Q_{2N}` at `2k`), plain Euclidean norm because both
sides are pointwise values of a profile at the *same* `xi` points.

| `N -> 2N` | rel. width coarse/fine | paired modes | `err(N)` | in fit |
|---|---|---|---|---|
| 8 -> 16 | 0.1250 / 0.1875 | 1 | 1.647e-03 | no (geometry mismatch) |
| 16 -> 32 | 0.1875 / 0.1875 | 125 | **5.143e-02** | yes |
| 24 -> 48 | 0.1667 / 0.1875 | 343 | 1.777e-01 | no (geometry mismatch) |
| 32 -> 64 | 0.1875 / 0.1875 | 1331 | **1.534e-02** | yes |
| 48 -> 96 | 0.1875 / 0.1875 | 4913 | **7.117e-03** | yes |
| 64 -> 128 | 0.1875 / 0.1875 | 12167 | **4.079e-03** | yes |

Log-log fit `err ~ C N^{-sigma}` on the four geometry-matched points:

```
sigma = 1.8249,   log C = 2.1117  (C = 8.26),   R^2 = 0.99937,   4 points.
```

Successive two-point slopes are `1.744` (16->32), `1.894` (32->48),
`1.933` (48->64) — drifting up toward 2, consistent with a second-order
quadrature degraded by the Fejer kink, and bracketed by the §1.5 operator
bounds (`1` on `C^1`, `2` on `C^2`).

**TM-22 label: `sigma = 1.82` is a diagnostic exponent, not a proven
convergence rate.** It is a four-point binary64 fit over a single octave-and-a-
half of `N` on one profile family, with no interval enclosure and no
resolution-adequacy precondition of the kind TM-22 requires before a
convergence fit may be trusted. Nothing below uses it as a rate; it is used
only once, explicitly flagged, as an order-of-magnitude cross-check in §3.3.

**Null control (mandatory, and it fires correctly).** Re-running the identical
protocol with the module's *unrescaled* integer phase slope — so the continuum
profile is no longer fixed as `N` grows — gives

```
0.2451 (16), 0.4994 (32), 0.7317 (48), 0.9307 (64),
sigma = -0.9686 (i.e. the error GROWS like N^{+0.97}), R^2 = 0.99840.
```

The control separates cleanly from the signal, so the decay in the table is a
property of the fixed continuum profile and not an artefact of the pairing,
the backend, or the norm.

**Two recorded failures, kept per repository discipline.** The rows
`8 -> 16` and `24 -> 48` compare two *different* continuum profiles, because
`floor(3N/16)/N` equals `3/16` only when `16 | N`. `24 -> 48` accordingly shows
a **larger** error (0.178) than its neighbours; `8 -> 16` shows a spuriously
tiny error (1.6e-03) because at `N=8` the box degenerates to a single mode per
carrier and the comparison has exactly one paired coefficient. Both are
recorded in the CSV, excluded from the fit, and neither is deleted.

## 3. One-stage Poincare map consistency

### 3.1 The operator

`P` is one stage of `run_renormalized_cascade`, applied as a *map* with a
**frozen** integrator so that `P` is well defined: evolve the full Fourier
Galerkin NS ODE (RK4, exact padded products) for the parabolic stage time
`tau / N_0^2` with the step count fixed by the base state and then reused
verbatim for every argument; apply the lattice doubling pullback
`w(k) = 4 v(2k)`; optionally drop the sub-front sea; renormalize the energy
back to `c_E / N_0`. Freezing the step count matters: with the cascade's
adaptive CFL rule a perturbed argument would silently change the integrator,
and the resulting jump would be reported as map sensitivity.

Both stages run for the *same renormalized time*: `ds = N^2 dt`, so
`Delta s = N_0^2 . tau / N_0^2 = tau` at either scale. `drop_below = None`
here, because a fixed *lattice* sea floor is a different `xi` floor at the two
scales (`3/4` at `N_0=4` versus `3/8` at `N_0=8`) and would inject an `O(1)`
mismatch that has nothing to do with discretization.

### 3.2 Pairing and results

Profiles are compared as `Psi_N(k/N) = N^2 u_hat(k)` on the common band
`|k|_inf <= min(K_4, K_8/2)`, `K_N = cutoff_N / 2` being the pullback support,
with `Psi_8` read at `2k`. Two pairs were run:

* **prescribed** — `(N_0=4, W=2, grid 64)` vs `(N_0=8, W=2, grid 80)`, the
  configuration named in the task. Note that `W=2` at both scales *halves* the
  relative cell width (`1/4` -> `1/8`), so this pair does **not** restrict one
  continuum profile.
* **profile-matched** — `(N_0=4, W=2, grid 64)` vs `(N_0=8, W=3, grid 96)`,
  which keeps `(W-1)/N_0 = 1/4` fixed and is the honest `R_N` comparison.
  `cutoff = 20` needs `grid > 80`, hence 96.

| pair | `c_E` | steps | `||P_4 - P_8|| / ||P_8||` | shape overlap |
|---|---|---|---|---|
| prescribed | 1 | 16 / 16 | 0.6493 | 0.99492 |
| profile-matched | 1 | 16 / 16 | **0.8355** | 0.99832 |
| prescribed | 100 | 63 / 21 | 0.7690 | 0.75243 |
| profile-matched | 100 | 63 / 48 | **0.7043** | 0.94128 |

**The map-level consistency error is `O(1)` — between 0.65 and 0.84 — at every
setting.** The near-unit shape overlaps show the discrepancy is almost entirely
in the band-restricted *amplitude*, not the direction.

### 3.3 Where the `O(1)` comes from (two candidate causes, not separated)

The lattice normalization gives the exact identity
`sum_k |Psi_N(k/N)|^2 = 2 c_E N^3`, which the runs reproduce to all displayed
digits (`||Psi_4|| = sqrt(128) = 11.3137` at `c_E=1`; `sqrt(12800) = 113.1371`
at `c_E=100`). Using it, the fraction of the *fine* profile's squared norm
that lives on the paired (coarse-resolvable) `xi` points is

| pair | `c_E` | paired / total fine `L^2` | value for a resolved profile |
|---|---|---|---|
| profile-matched | 1 | **3.73 %** | 12.5 % (= 1/8) |
| profile-matched | 100 | **5.14 %** | 12.5 % |
| prescribed | 1 | 99.91 % | 12.5 % |
| prescribed | 100 | 92.15 % | 12.5 % |

Neither matches. In the profile-matched pair the `N_0=8` pullback image puts
96 % of its energy *off* the coarse lattice: it is a grid-scale object that the
coarse `xi` lattice cannot see — precisely the TM-22 unresolved-peak signature.
In the prescribed pair the fine image is instead almost entirely *on* the even
sublattice, which is the opposite pathology and reflects that `W=2` at `N_0=8`
is a narrower cell, not the same profile.

Two causes are consistent with the data and **this budget cannot separate
them**:

1. **Plain under-resolution at `N_0 = 4`.** Extrapolating the §2.3 law
   `err ~ 8.26 N^{-1.82}` down to `N = 4` gives `0.66` — numerically almost
   exactly the measured map error. (Stated only as an order-of-magnitude
   coincidence check; per TM-22 an extrapolation outside the fitted window is
   *not* a prediction, and this one runs the wrong way, toward coarser `N`.)
2. **The decimation pullback is not a discretization of the continuum
   dilation.** Audit item S7.1 of
   [`audit_front_flow_seed.md`](../ideas_2026_08_01/audit_front_flow_seed.md)
   states that `w(k) = 4 v(2k)` samples only the even sublattice and discards
   ~7/8 of a generic child, whereas the repository's own (6.8),
   `(T_lambda U)_a(xi) = lambda^2 V_tot(lambda(a+xi))`, is a continuum dilation
   requiring *resampling/interpolation*. The 3.73 % figure is a direct,
   quantitative instance of that criticism.

Separating them needs `N_0 >= 16`, i.e. `cutoff = 34` and grid `>= 144`, with
its `2N_0 = 32` partner at `cutoff = 66` and grid `>= 288`. A `288^3` RK4 stage
is roughly `27x` the cost of the `96^3` stage already run and is far outside
both this sprint's budget and the task's explicit 64/80 resolution ceiling.
**Recorded as an open obligation, not attempted.**

## 4. Contraction margin along the measured trajectory

At `(N_0 = 4, W = 2, grid 64)`, front-only (`drop_below = 3`), stage-0 state
`u = ` the critical-energy angle-box parent, with the step count frozen as in
§3.1. Three random perturbations `delta` per `c_E`, each verified real
(reality defect exactly `0.0`) and divergence-free (defect `< 9e-17`), with
`||delta|| = 1e-6 ||u||`, and

```
L = || P(u + delta) - P(u) || / || delta || .
```

| `c_E` | steps | `L` (3 trials) | margin `1 - L` |
|---|---|---|---|
| 1 | 16 | 14.7081, 14.5956, 15.0138 | **-13.71 … -14.01** |
| 100 | 63 | 1.7335, 1.8060, 1.7869 | **-0.73 … -0.81** |

**`P` is expanding, not contracting, at the stage-0 state at both `c_E`.**
Note `||P(u)|| = ||u||` here (both are normalized to the same critical energy
`c_E/N_0`), so `L` is exactly the relative shape sensitivity and needs no
further scaling.

The conclusion is one-sided and therefore robust: three random directions can
only *understate* the top singular value, so the true local operator norm
satisfies `||DP(u)|| >= max_trial L > 1`. No amount of extra sampling can turn
this into a contraction. Note also that the map is expanding *while* the
pullback gain is decaying (pilot A: final gains `1e-10` to `1e-4`) — the
projective renormalization re-amplifies whatever survives the pullback by
`1 / sqrt(gain)`, which is exactly the mechanism that makes `L > 1` compatible
with a collapsing orbit.

Scope: `L` is measured at the stage-0 state only. "Along the measured
trajectory" cannot mean more than this, because pilot A's trajectories collapse
into the binary64 noise floor within 6–11 stages and there is no orbit to
follow.

## 5. Pre-registered acceptance rule and verdict

**Rule (fixed before any number in this note was produced).** The orbit
shadowing proof candidate survives only if

```
contraction margin (1 - L)  >  map-level consistency error.
```

This is the standard finite-resolution shadowing requirement: if the lattice
map differs from its refinement by `eps` and the map contracts with factor
`L < 1`, the shadowing distance is `eps / (1 - L)`; the rule demands that this
stay below the perturbation scale.

**Measured.**

```
map-level consistency error (worst, profile-matched)  =  0.8355
local Lipschitz factor (worst)                        = 15.0138
contraction margin 1 - L                              = -14.0138
rule satisfied                                        =  NO
```

The rule fails at both `c_E`, by a wide margin, and it fails in the worst
possible way: the margin is not merely too small, it is **negative** — there is
no contraction to trade against consistency error at all. Independently, the
premise of the rule is empty: pilot A found no periodic orbit in the scanned
box, so there is no orbit, no monodromy and no Floquet radius, and hence
nothing to shadow.

### 5.1 Verdict

> **KILLED — "orbit shadowing proof candidate" (front-flow lane, lattice
> doubling map).** No orbit exists in the scanned box (pilot A); the one-stage
> renormalized map is expanding (`L = 1.73`–`15.01`, margin `1 - L` from
> `-0.73` to `-14.01`); and the map-level consistency error is `O(1)`
> (`0.65`–`0.84`) at the only doubling this budget permits. The pre-registered
> rule `(1 - L) > eps` fails at both `c_E` by more than an order of magnitude.

### 5.2 What survives

The kill is of the *proof candidate*, not of the scaffolding. Three pieces
survive as reusable, and one as a corrected fact:

1. **`R_N` / `I_N` with sharp orders** (§1): `R_N I_N = id` exactly; `R_N`
   preserves reality and divergence-freeness exactly at every `N`;
   `|| I_N R_N - id || <= sqrt(3) ||grad U||_inf / N` on `C^1`; and the
   quantified obstruction (D) — `I_N` breaks the Leray constraint at first
   order, `|(a+xi) . I_N R_N U| <= sqrt(3) ||U||_inf / N`.
2. **The lattice Leray convolution IS consistent with continuum `Q`** at a
   measured diagnostic exponent `sigma = 1.82` (`R^2 = 0.9994`), with a null
   control that correctly fails. The *algebra* of the lane is not the problem.
3. **The `2k <-> k` doubling pairing** is the correct `R_N`-level comparison,
   and the exact profile identity `sum_k |Psi_N|^2 = 2 c_E N^3` gives a free
   resolution diagnostic (paired-sublattice `L^2` fraction versus `1/8`).
4. **Corrected fact for the portfolio.** §5 of
   [`discovery_2026_08_01_front_flow_portfolio.md`](../discovery_2026_08_01_front_flow_portfolio.md)
   lists "interpolating (non-decimating) pullback" as merely one of several
   open realizations. This measurement upgrades it: at `N_0 = 4 -> 8` the
   decimated pullback image is 96 % invisible to the coarse `xi` lattice, so
   **no orbit search built on the decimating pullback can be given a continuum
   reading**, whatever it finds. Obligation 6 of that note's §8 should be read
   as a prerequisite, not an option.

### 5.3 Open obligations handed on (not attempted here)

* Implement the interpolating pullback `(T_2 U)_a(xi) = 4 V_tot(2(a+xi))` via
  `I_N` and re-measure §3 with it. Carry (D) explicitly, or re-project.
* Re-run §3 at `N_0 >= 16` (grid `>= 144` / `>= 288`) to separate
  under-resolution from the decimation defect. Costed at ~27x the `96^3` stage.
* Replace the `L^2` profile metric by a critical one (audit S7.5: the front
  profile is unbounded in `L^2`, so `sup_xi xi^2 e(xi)` or a local-in-`xi`
  metric is required before any orbit or contraction statement is meaningful).
  Every number in §3 and §4 is in `L^2` and inherits this caveat.

## 6. Reproduction

```
.venv/Scripts/python.exe experiments/run_poincare_shadowing.py \
    --output-dir outputs/verification_sprint_v1/shadowing \
    --csv docs/research_notes/verification_sprint_v1/poincare_shadowing_scaling.csv
.venv/Scripts/python.exe -m pytest tests/test_poincare_shadowing.py
```

Runtime ~13 min on the reference machine (numpy only; scipy is not installed
and is not used). RNG seed `20260802`, fixed in the script. The digest sidecar
`summary.json.sha256` fixes the raw record.

*Binding non-claims: no singularity, no regularity theorem, no Clay statement,
no continuum enclosure. `sigma = 1.82` is a diagnostic exponent, not a proven
rate. Every quantity is binary64 and labelled `float` in the CSV.*
