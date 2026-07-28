# Early Hou snapshots: conditional singular-scaling gate result

## Status

This is a negative diagnostic result for the available early-time snapshots. It
does not exclude later Hou dynamics and does not prove global regularity.

The input was the previously exported critical-`L^3` snapshot table with five
times for each of three uniform grids. The fit scanner varied the proposed
singular time over

```text
T_last + 1e-6 <= T <= T_last + 2e-3
```

and used every suffix window containing at least three points. For every
candidate time and window it fitted

\[
A(t)\sim (T-t)^{-\alpha},\qquad
L_r(t)\sim(T-t)^{\beta_r},\qquad
L_z(t)\sim(T-t)^{\beta_z}.
\]

The amplitude was the maximum physical velocity and the widths were RMS widths
of the physical critical density `|u|^3`.

## Result

No fitted window on any grid passed the conditional finite-energy,
critical-growth, finite-dissipation and generic-balance gates.

| grid | gate-passing fits | best scanned `alpha` | best `beta_r` | best `beta_z` |
|---|---:|---:|---:|---:|
| 65 x 128 | 0 | 0.0743 | -0.0125 | -0.0101 |
| 129 x 256 | 0 | 0.1250 | -0.0116 | -0.00865 |
| 193 x 384 | 0 | 0.1390 | -0.0114 | -0.00824 |

The best post-selected fits put the proposed singular time almost immediately
after the final sample. Even then, both fitted width exponents are negative.
Thus the measured critical-density widths are expanding, not shrinking.

The fitted amplitudes also have `alpha < 1/2`, which lies inside the
`L^2_t L^infinity_x` Serrin-regular range for this simple global power-law
model.

## Interpretation

The early snapshots show velocity-amplitude growth but do not show a resolved
shrinking critical core. They therefore do not enter any of the three generic
single-core families:

1. Euler-like Type-II time/inertial balance;
2. anisotropic parabolic triple balance;
3. quasi-steady viscous/inertial balance.

This is consistent with the direct critical-`L^3` observation that the global
represented-domain `L^3` norm decreases over the early interval and the RMS
critical widths increase.

The conclusion is limited to the available early interval. A later transition
could in principle produce shrinking widths, critical-norm growth, a shell
cascade or an outer-field contribution. Such a transition must be demonstrated
with resolved, converged data rather than inferred from maximum vorticity.

## Required next experiment

Do not fit a blow-up law to the current early snapshots. The next useful run
must provide, at every accepted step or sufficiently dense checkpoint:

- physical maximum velocity;
- at least three independent radial and axial scale definitions;
- global and shell-decomposed `L^3` mass;
- energy and dissipation;
- time, inertial-pressure and viscous PDE-block norms;
- cancellation ratio for the inertial-pressure block;
- full space/time/domain/integrator refinement metadata.

A fit may be promoted only when all measured scale definitions are positive,
shrink under refinement, and place the candidate in a stable exponent family
across preregistered windows and candidate-time ranges.
