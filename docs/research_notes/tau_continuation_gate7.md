# Gate 7 — leaving the Picard region, and what was found there

Run: `outputs/tau_continuation_gate7_v1/`.
Driver: `experiments/run_tau_continuation.py`, config `configs/tau_continuation_gate7.json`.

**Headline: the gate passed and every candidate was rejected.**  All three
preregistered families reach `τ = A t = 1` and leave the first Picard iterate by
17–36 %, so the computations are genuinely nonlinear.  None of them concentrates.
The critical `L³` norm decays monotonically in every one of the eighteen runs,
and the dyadic shell count *increases*, meaning the field spreads.  Nothing here
bears on the Clay problem, and nothing here is a claim about it.

## 1. Why the previous thirty-two runs never left the linear regime

Gate 6 ended with thirty-two runs all inside the first Picard iterate.  The
reclassification in §1 of the driver explains it in one number.

The datum is `A_nom · χ((r/R)²) · χ((z/Z)²) · ζ/(1+cζ²)`.  Its peak is not
`A_nom`: the radial bump contributes `χ(0) = e^{-1} ≈ 0.368` and the odd axial
profile peaks at about `0.13`.  So

```
A_eff = max |u₁(0)| ≈ 0.048 · A_nom .
```

The sweep's "amplitude" overstated the field the equations see by a factor of
about **twenty-three**.  Running it to `t = 0.024` therefore reached

```
τ_max = A_eff · t = 0.0233  (32 points, min 0.00163)
```

not the `O(1)` it looked like.  Rerunning `deduplicate_settings` on the sweep in
`(Re, aspect, c, τ)` coordinates:

| quantity | value |
|---|---|
| points | 32 |
| distinct computations | 32 (no exact duplicates) |
| distinct dimensionless shapes | **6** |
| distinct Reynolds numbers | 16 |
| Reynolds range | 8.69 – 279.83 |
| dimensionless time reached | 0.00163 – **0.0233** |

So the sweep was not wasteful in the duplicate sense — no two points were the
same computation — but it varied only six shapes, and it stopped **43× short in
dimensionless time**.  The `A²` scaling of `max|ω₁|` that Gate 6 diagnosed after
the fact was not a subtle finding; it was the inevitable consequence of never
getting past `τ = 0.023`.

The per-point table is `legacy_reclassified.csv`.

## 2. The single parameter

`nondimensional.py` carries the substitution `r = Lρ, z = Lζ, τ = At, u₁ = AU,
ω₁ = (A/L)W, ψ₁ = ALΨ` through the audited system.  Every term of the swirl
equation carries `A²` except the viscous one, which carries `νA/L²`; every term
of the vorticity equation carries `A²/L` except the viscous one, `νA/L³`.
Dividing leaves exactly one parameter:

```
Re = A L² / ν .
```

Amplitude, length and viscosity are therefore **not** three search directions.
`Setting.dimensionless_key` makes that checkable and `deduplicate_settings`
enforces it.

The scaling exponents (`E ~ A²L⁵`, `‖u‖_{L³} ~ AL²`, `|ω| ~ A`, width `~ L`,
`t ~ A⁻¹`) are in `ScalingExponents` and are cross-checked against each other in
`tests/test_gate7_modules.py` rather than merely asserted.

## 3. How far `τ` has to go, measured rather than assumed

Pure swirl advects nothing.  The meridional velocity that does the advecting is
generated through the elliptic solve, at order `τ`, and the elliptic kernel
carries `1/(8π²)`.  So the advection-to-source ratio grows quadratically.
Measured on the flagship:

| `τ` | advection/source | `C = ratio/τ²` | distance to 1st Picard | `τ` for 10 % |
|---|---|---|---|---|
| 0.025 | 4.40e-4 | 0.705 | 3.43e-4 | 0.377 |
| 0.05 | 1.73e-3 | 0.691 | 1.35e-3 | 0.381 |
| 0.1 | 6.64e-3 | 0.664 | 5.22e-3 | 0.388 |
| 0.2 | 2.47e-2 | 0.617 | 1.95e-2 | 0.403 |
| 0.4 | 8.58e-2 | 0.536 | 6.84e-2 | 0.432 |
| 0.7 | 2.16e-1 | 0.441 | 1.73e-1 | 0.476 |
| 1.0 | 4.01e-1 | 0.401 | 2.98e-1 | 0.499 |

The quadratic law holds cleanly at small `τ` and softens as the response
saturates.  The extrapolation `τ_* = τ√(0.1/ratio)` predicts 0.38–0.50 and the
measured crossing of the 10 % departure bar is at `τ ≈ 0.5`.  **`τ = O(1)` is
the right target after all** — the earlier estimate of `τ ≈ 9` was an artefact
of measuring `τ` with the nominal amplitude.

## 4. The Picard ladder

The distance to the Picard iterates is *measured*, not inferred.
`PicardLadder` integrates

* level 0 — pure diffusion (vorticity stays identically zero),
* level 1 — nonlinear terms frozen at level 0,
* level 2 — nonlinear terms frozen at level 1,
* the full nonlinear solution,

**together**, with the same integrator and the same accepted steps.  There is no
interpolation between trajectories and no differencing of stored snapshots; the
exact solver right-hand side is available at every accepted step because it is
computed there.

Two structural facts the ladder makes visible.  Level 0's vorticity is exactly
zero, so its velocity is exactly zero — which is why level 1 reduces to
"vorticity grows linearly from the stretching source `∂_z(u₁²)`", the familiar
first-Picard picture.  And level 2 is the first iterate that feels advection at
all.

## 5. The preregistered departure gate

Thresholds were recorded in the config before any run
(`preregistration.json` in the output directory).  Evaluated on the flagship at
the finest resolution:

| check | value | bar | verdict |
|---|---|---|---|
| distance to 1st Picard | 3.61e-1 | ≥ 0.10 | pass |
| advection measurable | 4.89e-1 | ≥ 0.02 | pass |
| resolution stability | 5.87e-2 | ≤ 0.20 | pass |
| integrator agreement (RK4 vs SSPRK3) | 4.41e-9 | ≤ 0.05 | pass |
| time-step agreement (`dt` halved) | 8.54e-12 | ≤ 0.05 | pass |
| front points, radial | 10.0 | ≥ 10 | pass |
| front points, axial | 11.0 | ≥ 10 | pass |
| domain tail below signal | 1.08e-90 | ≤ 1e-6 | pass |
| divergence free | 1.65e-3 | ≤ 0.05 | pass |

### The one amendment, and why it is not a moved threshold

The first pass failed exactly two checks, both front-resolution: 6.5 radial and
7.0 axial points on the 97×193 grid against a bar of 10.  **The bar was not
moved.**  A grid was added (145×289) and the gate re-evaluated there, with the
resolution-stability comparison shifted up accordingly so it remains a genuine
two-grid comparison.  The amendment, its reason and the fact that no threshold
changed are recorded in `configs/tau_continuation_gate7.json` under `amendments`
and copied into `summary.json`.

Spending more computation to meet an unchanged bar is not the thing the
preregistration rule forbids; changing the bar after seeing the number is, and
that did not happen.

## 6. The Reynolds continuation — every candidate rejected

Preregistered series `Re ∈ {10, 25, 50, 100, 200, 400}` × families S, A, H, each
to `τ = 1`.  Final values at `τ = 1` (medium grid):

| run | dist. 1st Picard | adv/src | visc/src | `L³` growth | Δ shell count | min width ratio |
|---|---|---|---|---|---|---|
| S_Re10 | 0.0094 | 0.0101 | 2.9794 | 0.1945 | +0.3150 | 1.2026 |
| S_Re25 | 0.0554 | 0.0594 | 1.8076 | 0.4378 | +0.2333 | 1.0504 |
| S_Re50 | 0.1262 | 0.1662 | 1.0677 | 0.6278 | +0.1534 | 1.0116 |
| S_Re100 | 0.2001 | 0.3491 | 0.6734 | 0.7722 | +0.0987 | 0.9998 |
| S_Re200 | 0.2578 | 0.5006 | 0.3719 | 0.8652 | +0.0594 | 0.9977 |
| S_Re400 | 0.2947 | 0.6078 | 0.2009 | 0.9191 | +0.0486 | 0.9981 |
| A_Re10 | 0.0027 | 0.0030 | 3.4524 | 0.1159 | +0.5613 | 0.8215 |
| A_Re25 | 0.0186 | 0.0178 | 2.4297 | 0.3135 | +0.3817 | 0.9218 |
| A_Re50 | 0.0420 | 0.0502 | 1.8379 | 0.5117 | +0.1963 | 0.9593 |
| A_Re100 | 0.0800 | 0.1073 | 1.0654 | 0.6895 | +0.0558 | 0.9792 |
| A_Re200 | 0.1249 | 0.2215 | 0.6155 | 0.8166 | **−0.0298** | 0.9881 |
| A_Re400 | 0.1677 | 0.3284 | 0.3586 | 0.8953 | **−0.0612** | 0.9925 |
| H_Re10 | 6.20e-4 | 0.0013 | 12.1340 | 0.0366 | +0.0985 | 1.5673 |
| H_Re25 | 0.0063 | 0.0080 | 3.6579 | 0.1199 | +0.1174 | 1.2111 |
| H_Re50 | 0.0290 | 0.0312 | 2.0395 | 0.2476 | +0.0732 | 1.0891 |
| H_Re100 | 0.0902 | 0.1117 | 1.7755 | 0.4233 | +0.0889 | 1.0440 |
| H_Re200 | 0.1848 | 0.2406 | 1.4885 | 0.6051 | +0.0750 | 1.0360 |
| H_Re400 | 0.2975 | 0.4011 | 1.1299 | 0.7514 | +0.0646 | 1.0202 |

Read across the table:

* **The departure grows monotonically with `Re`,** exactly as it should: the
  nonlinear term is `Re` times the viscous one in dimensionless variables.  At
  `Re = 400` all three families are outside the first Picard iterate.
* **`L³` decays in every run.**  Growth ratios rise toward 1 as `Re` grows
  (S: 0.772, 0.865, 0.919 at `Re` = 100, 200, 400) but the increments shrink
  (0.093, 0.054), so the sequence appears to saturate **below** 1.  Extrapolating
  from three points is not a proof, but the direction is unambiguous and it is
  the opposite of what a critical-concentration candidate needs.
* **The shell count increases in sixteen of eighteen runs.**  The exceptions are
  family A at `Re ≥ 200`, whose asymmetric quadrupole does draw the field
  inward slightly (−0.03, −0.06).  That is the only sign of concentration
  anywhere in the sweep, and it is far too small to matter.
* **Widths are flat.**  The narrowest ratio in the whole sweep is A_Re10 at
  0.82, and that is a diffusive artefact at the most viscous setting, not a
  nonlinear contraction.

### Promotion decision: **no candidate promoted**

The preregistered criteria are

* *critical concentration*: `L³` growth ≥ 1.05 **and** Δ shell count ≤ −0.05
  **and** narrowest width ratio ≤ 0.90;
* *multiscale*: Δ shell count ≥ 0.30 **and** outer-shell growth ≥ 2.0.

The best run in the sweep, A_Re400, satisfies exactly one of the six conditions
(Δ shell count = −0.0612 ≤ −0.05) and fails the other two of its criterion:
`L³` growth 0.8953 (needs ≥ 1.05), narrowest width ratio 0.9925 (needs ≤ 0.90).
Every rejection carries its numbers in `summary.json` under `promotions`.

**The rejection is the result.**  Three explicit whole-space families, followed
into the genuinely nonlinear regime across a factor of forty in Reynolds number,
all decay.

## 7. What the families are

Preregistered in `initial_data.py` before any run.

* **S** — the symmetric baseline, byte-identical to the legacy `SwirlFamily`, so
  every earlier result remains comparable.  Odd in `z`.
* **A** — axially asymmetric: two bumps at different centres with different
  amplitudes, widths and parities.  Built specifically so `u₁²` is not even in
  `z` and the axial quadrupole does **not** degenerate.
* **H** — Hou-shaped: concentrated near the axis (bump raised to the power 2.5,
  which preserves `C^∞` since `χ(s)^k = exp(-k/(1-s))`) and thin in `z`.  It
  keeps the two features of the Hou setup that matter and imports **neither**
  the finite-cylinder wall condition **nor** the two-stage viscosity protocol.

### Which moments degenerate, and why — analytically

At leading order `ω₁ = t ∂_z(u₁²)`, so the multipole moments are integrals of a
`z` derivative against `1`, `z`, `z²`, `r²` with `dV₅ = |S³| r³ dr dz`:

| moment | value | degenerate? |
|---|---|---|
| monopole | `∫ ∂_z(u₁²) dV₅ = 0` | **always** — the `z` integral of a `z` derivative |
| axial dipole | `−∫ u₁² dV₅` | **never** |
| axial quadrupole | `−2∫ z u₁² dV₅` | iff `u₁²` even in `z`: **S and H yes, A no** |
| transverse quadrupole | `¼∫ r² ∂_z(u₁²) dV₅ = 0` | **always** |

This is the structural reason Gate 6's multipole boundary comparisons kept
degenerating, and it is why family A exists.  Measured: A's normalised axial
quadrupole is −0.308; S and H give −5.6e-18 and 3.2e-19, i.e. zero.

The monopole is not permanently zero — it grows at `O(t²)` through advection,
whose five-dimensional divergence `2u^r/r` is nonzero — but at first order it
vanishes for every compactly supported datum, so it can never serve as a
first-order boundary discriminator.

## 8. The space-time slab certificate

`slab_certificate.py` encloses one accepted step over **whole cells at every
instant**, not nodes at two instants.

* Temporal enclosure: cubic Hermite from the two endpoint states and the two
  exact solver right-hand sides.  The basis ranges are exact rationals —
  `h₀₀ + h₀₁ ≡ 1` with both in `[0,1]`, so the value part is a convex
  combination and lies in the endpoint hull; the tangent bases have exact
  extrema `±4/27`.  The tests check these against the polynomials rather than
  trusting them.
* Spatial enclosure: corner hull inflated by the measured cell variation.
* All arithmetic in `fractions.Fraction` with outward rounding.

Measured on the flagship at `τ = 0.4`, 706 cells:

| bound | value |
|---|---|
| Poisson residual | 1.57e-15 |
| divergence | 2.13e-4 |
| trapezoid local defect | 2.35e-6 |
| Simpson local defect | 8.25e-10 |
| Hermite remainder | 1.54e-7 |

The independent checker re-derives every relation the payload asserts between
its own numbers and rejects: an inverted enclosure, a hypothesis marked
`proved`, a fourth-derivative bound inconsistent with the stated remainder, and
a payload with the Clay disclaimer removed.  Fed garbage it returns a verdict
rather than raising, because a crash is indistinguishable from a checker bug.

### The two hypotheses, stated rather than hidden

`H1` (cell Lipschitz) and `H2` (Hermite remainder) are **not proved**.  They are
named in the payload, marked `proved: false`, and the checker rejects any
payload claiming otherwise.  `H2`'s constant `M₄` is estimated from the measured
Simpson defect and inflated — an estimate, not a derivation.

Consequently the certificate establishes, exactly, a property of the *computed*
space-time object, conditional on `H1` and `H2`.  It is not a statement about
the continuous Navier–Stokes solution.  The missing step is `HS-5` of
[`hs_error_propagation.md`](hs_error_propagation.md).

## 9. Related work in this session

* [`hs_error_propagation.md`](hs_error_propagation.md) — moving the error
  estimate from `L^∞` to `H^s`, `s > 5/2`.  An **incomplete** derivation: §3.1,
  §3.4 and the product-difference identity are proved; `HS-1`, `HS-2a`, `HS-2b`,
  `HS-4`, `HS-5`, `HS-6` are named with their statements and their unproved
  constants.  §7 states flatly what the existing `L^∞` certificate does and does
  not establish.
* `formal/NSSingularity/TimeDependentGalerkin.lean` — F-7c closed by the direct
  route; mathlib's `IsPicardLindelof` was already time-dependent, so the
  autonomisation was never needed.

## 10. Limitations

* This is a resolved numerical continuation, not a proof.  Nothing here
  establishes or refutes a Clay singularity, and no such claim is made.
* The departure gate certifies distance to the Picard iterates of the *discrete*
  system on a given grid.  That the computation has left linear response is
  measured; that the continuous solution has is not.
* The `L³` saturation argument in §6 extrapolates from three Reynolds numbers.
  It is a direction, not a theorem.
* Viscosity is fixed and positive in every run.  No time-dependent or two-stage
  viscosity protocol appears anywhere in the candidate computations.
