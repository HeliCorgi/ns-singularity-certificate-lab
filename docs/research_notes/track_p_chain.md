# Track P chain: composing slab certificates into a certified interval

**Status (2026-07-31, turn 10).**  Conditional PDE certificates assuming
EXT-P1/P2/P3; the finite-dimensional Galerkin enclosure and the scalar control
ODE are verified unconditionally.  Nothing here bears on the Clay problem, and
a chain that stops proves nothing about singularities.

Code: `src/ns_certificate_lab/torus_chain.py`.
Lean skeleton: `formal/NSSingularity/TrackPChain.lean`.
Preregistration: `configs/track_p_chain.json` (thresholds fixed before any
run).  Outputs: `outputs/track_p_chain_v1/`.

## 1. The two-slab theorem

Fix `t0 <= t1 <= t2`.  Suppose slab `n` gives (conditional on EXT-P1/P2/P3)

```
||u(t) - c1(t)||_{H4dot} <= R1(t)   on [t0, t1],  R1(t0) = delta_0 ,
```

where `c1` is the Galerkin trajectory through the exact rational datum `y1`.
Let `y2` be any exact rational point with `||c1(t1) - y2||_{H4dot} <= eta`
(the recentring), and let `delta_1 >= R1(t1) + eta`.  Then by the triangle
inequality

```
||u(t1) - y2||_{H4dot} <= R1(t1) + eta <= delta_1 ,
```

so the slab-2 certificate — whose control ODE starts at `R2(t1) = delta_1` —
applies to the *same* solution `u` (uniqueness clause of EXT-P1), giving the
piecewise bound on `[t0, t2]`.  Induction over finitely many slabs is
immediate; the certified interval is the union of the slabs and the certified
radius is the per-slab tube.  The finite-inequality skeleton (two-slab
composition, list induction, transfer triangle, discrete Gronwall, union
bound) is in `TrackPChain.lean`; the comparison between the Dini derivative
of `||u - c_n||` and the scalar tube is the same EXT-P2-conditioned step as in
the single-slab certificate, now started from `delta_n` instead of `0`.
The control inequality itself needs no smallness of the initial radius: its
constants depend only on the slab enclosure of the centre trajectory.

## 2. Exact restart points, and why there is nothing to wrap

The classical wrapping effect afflicts interval integrators that propagate a
*box* through the flow: each step re-encloses an increasingly non-axis-aligned
set in axis-aligned coordinates, and the width inflates exponentially.  The
chain avoids the phenomenon *structurally* rather than fighting it:

1. every slab starts from an **exact rational, exactly divergence-free**
   datum `y_n` (width zero — there is no box to wrap);
2. the endpoint of the centre trajectory is enclosed once per slab by an
   interval Taylor form with Lagrange remainder,

   ```
   u_a(h) in  sum_{i<=m} c_i h^i  +  d_{m+1} h^{m+1} ,
   ```

   with `c_i` the exact time-Taylor coefficients through `y_n` (interval
   arithmetic seeded with a degenerate interval) and `d_{m+1}` the
   coefficient recurrence evaluated on the proved Picard box — valid because
   the `(m+1)`-st Taylor coefficient of the trajectory at the unknown
   intermediate time is a fixed polynomial of a state the box contains;
3. the enclosure is collapsed to its midpoint, dyadically rounded
   (`rounding_bits`), Leray-projected back to exact divergence-freeness
   (a mode-wise orthogonal projection, so it contracts every `H^ndot` and the
   projection step never increases the error), and the entire discarded
   width — remainder, rounding, projection — is charged to the **scalar**
   radius as the transfer

   ```
   delta_{n+1} = R_n(t_{n+1}) + ||u_a(t_{n+1}) - y_{n+1}||_{H4dot} .
   ```

The vector enclosure therefore dies at every boundary; only a scalar crosses.
What the scalar loses is directional information — it treats the error as an
isotropic ball, which is exactly the crudeness measured by the linear
coefficient of the control ODE, and exactly what the failure classification
reports (section 5).

**Measured floors.**  The recorded per-slab transfer sits at `~1.4e-8`
(64-bit runs): this is the accumulated outward-rounding slop of the interval
`sobolev_sq`, not the true remainder (the order-6 Taylor remainder at
`h = 1/2048` is around `1e-23`; the coefficient-width diagnostics show it).
Against per-slab tube radii of `1e-3` and up, the transfer contributes less
than `1e-4` of the recurrence and is irrelevant to the horizon.

**Preregistered Lohner/QR triggers** (fixed in the config before any run):
wrapping amplification `>= 2` on three consecutive slabs; interval widths not
reduced by halving the slab; the coefficient enclosure ending the proof while
the point residual is small.  The measured amplification is `~1.00` on every
slab of every chain — the transfer is *at* the remainder+rounding floor — so
the triggers did not fire and no Lohner/QR frame was introduced.  This is the
honest outcome of the turn-10 instruction: measure first, introduce the frame
only if the measurements demand it.  A Lohner frame would represent the
orientation of a propagated box; here no box is propagated.

## 3. The growth law and the horizon ceiling

On each slab the tube obeys `R' <= a R + b R^2 + eps` with
`a = -nu + 9(K1 + K2)`, `b = 135 Adot`, `eps = ||tail||_{H4dot}`.  While
`b R^2` is subdominant the endpoint recurrence is essentially

```
delta_{n+1} ~ e^{a h} delta_n + (eps/a)(e^{a h} - 1) ,
```

so the chain grows by the factor `e^{a h}` per slab and the reachable horizon
is capped near the blow-up time of the Riccati majorant,

```
T* ~ (1/a) * log( a^2 / (b * eps) )        (delta_0 = 0) .
```

For P1 at `cutoff_sq = 4` the slab constants give `a ~ 1.24e3`,
`b ~ 5.6e2`, `eps ~ 1.9`, hence `T* ~ 5.7e-3` — and the chains do approach
this ceiling before the relative-radius cap stops them.  Three facts keep the
interpretation honest:

- the ceiling is a property of the **crude self-contained constants**
  (`9(K1+K2)` is an l1-type bound that overestimates the true logarithmic
  norm of the linearised flow by orders of magnitude), not of the solution;
  the trajectory itself is *decaying* in every norm we track;
- halving slabs approaches `T*` asymptotically but cannot cross it; the
  productive lever is `a` (sharper Kato-type constants, the recorded
  Morosi-Pizzocchero tradeoff) or `eps` (a larger Galerkin band — with the
  measured caveat that at fixed `h` the uniform Picard box inflates the tail
  bound when the band grows);
- accordingly every horizon-limited stop classifies as
  `control_linear_coefficient`, and the classification is part of the
  payload, checker-enforced, with the note that a stop proves nothing about
  singularities.

## 4. What one chain certifies

Conditional on EXT-P1/P2/P3 (recorded faithfully, never axiomatised, with
the uniqueness clause of EXT-P1 identifying the solution across boundaries):
a strong periodic solution exists on `[0, T_N]` and stays within the
per-slab tube of the piecewise Galerkin centre in `H4dot` (and within
`sqrt8` times that in `H3`).  Unconditionally: every Picard box, every
constant assembly, every scalar tube, every restart identity and the whole
`delta` recurrence are verified in exact rational arithmetic by
`verify_chain_certificate`, which rebuilds each link from its stored datum
and demands exact agreement.

Alongside the tube the chain tracks physical quantities at every restart
point — energy, enstrophy, `H3dot`/`H4dot`, shell energies and exact shell
fluxes, top-shell share, a numeric analyticity-slope proxy, a numeric
vorticity sup with a certified upper bound — each with a Lipschitz transfer
radius `|Q(u) - Q(u_a)| <= L_Q R` proved from `||w||_{Hndot} <= R` (`|k| >= 1`)
and `||w||_inf <= Adot ||w||_{H2dot}`:

| quantity | Lipschitz bound `L_Q` | status |
|---|---|---|
| `H4dot`, `H3dot` norm | `1` | exact centre, certified radius |
| energy | `||u_a||_0 + R/2` | exact centre, certified radius |
| enstrophy | `2||u_a||_1 + R` | exact centre, certified radius |
| global `L3` (normalised measure) | `Adot` (via `L3 <= Linf`) | numeric centre, certified radius |
| vorticity sup | `Adot` (via `||curl w||_{H2dot} <= ||w||_{H3dot}`) | numeric/upper centre, certified radius |
| shell energy | `sqrt(2 E_s(u_a)) + R/2` | exact centre, certified radius |
| shell flux, analyticity slope | — | numeric diagnostics only |

No tracked quantity grows along the certified chains (the data are decaying
Galerkin trajectories at these viscosities); nothing is called "certified
growth", and nothing here may be called that unless a *lower* bound from the
tube demonstrates it.

## 5. Failure taxonomy (preregistered)

Every stopped chain carries exactly one of: `galerkin_tail`,
`time_integration_defect`, `datum_transfer`, `interval_wrapping`,
`control_linear_coefficient`, `control_quadratic_term`,
`ext_hypotheses_insufficient`, `approximate_norm_growth` — or the two
successful stops `target_reached` / `slab_budget_exhausted`.  The classifier
compares the derivative shares `a*delta`, `b*delta^2`, `eps`, `transfer/h`
at the stopped slab and names the largest; the checker rejects unregistered
reasons.  "証明区間終了" は特異点の主張ではない — the payload carries this
as a machine-checked sentence.

## 6. Results (turn-10 runs, `outputs/track_p_chain_v1/`)

Thirteen chains, all built from the preregistered config, all verified by the
full-recomputation checker, all stops classified
`control_linear_coefficient` (as predicted in section 3), no Lohner trigger
fired anywhere:

| chain | slabs | certified `T` | vs single slab `1/2048` | max `R` (`H4dot`) |
|---|---|---|---|---|
| P1, `nu = 1/4`   | 13 | `79/16384  ~ 4.82e-3` | 9.9x  | 0.92 |
| P1, `nu = 1/10`  | 13 | `5/1024    ~ 4.88e-3` | 10.0x | 0.93 |
| P1, `nu = 1/40`  | 12 | `5/1024    ~ 4.88e-3` | 10.0x | 0.92 |
| P1, `nu = 1/100` | 12 | `5/1024    ~ 4.88e-3` | 10.0x | 0.92 |
| P2, `nu = 1/4`   | 14 | `45/16384  ~ 2.75e-3` | 5.6x  | 1.62 |
| P2, `nu = 1/10..1/100` | 15 | `91/32768 ~ 2.78e-3` | 5.7x | 1.74-1.75 |
| P3, `nu = 1/4`   | 12 | `33/16384  ~ 2.01e-3` | 4.1x  | 2.57 |
| P3, `nu = 1/10..1/100` | 13 | `67/32768 ~ 2.05e-3` | 4.2x | 2.92-2.93 |
| **long**: P1, `nu = 1/10`, `h = 1/8192` | **41** | `321/65536 ~ 4.90e-3` | **10.0x** | 0.92 |

All five preregistered acceptance criteria hold: every chain has at least 4
slabs; at least one chain exceeds 8x the reference single-slab length (the
P1 chains and the long chain reach 10.0x); the long run exceeds the 32-slab
target with 41 chained slabs; the checker verifies every emitted chain; and
every termination is classified.  The P1 horizons at four different
viscosities agree to three digits — direct evidence that the stop is the
`9(K1+K2)` Riccati ceiling of section 3 (viscosity enters `a` only through
the negligible `-nu`), not any property of the solutions, which decay in
every tracked quantity.  P2 and P3 stop earlier for the same reason at
larger `a` (their data are larger in the `M_j` sums).  The adaptive
controller's endgame is visible in the non-dyadic final times: the last one
or two slabs run at half or quarter step before the relative-radius cap
(`R > 0.5 ||datum||`) ends the chain.
