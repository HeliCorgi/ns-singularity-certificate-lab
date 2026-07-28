# Transformed-PDE Term-Balance Progress

## Status

This work does **not** solve the Navier–Stokes Millennium Problem and does not
produce a singularity candidate.  It adds a falsifiable diagnostic for deciding
whether a resolved checkpoint series exhibits one of the term balances required
by the Type-II scaling families already recorded in the repository.

## Implemented diagnostic

For each checkpoint series, the code constructs nonuniform three-point time
derivatives and evaluates the audited transformed equations

\[
 u_{1,t}+u^r u_{1,r}+u^z u_{1,z}
 =2u_1\psi_{1,z}+\nu\mathcal L_5u_1,
\]

\[
 \omega_{1,t}+u^r\omega_{1,r}+u^z\omega_{1,z}
 =\partial_z(u_1^2)+\nu\mathcal L_5\omega_1.
\]

It reports physical-volume-weighted RMS norms of

- the time derivative;
- advection;
- source/stretching;
- diffusion;
- the combined inviscid block `advection - source`;
- the complete residual.

The same quantities are computed globally and in an elliptic core defined by
two RMS widths of the physical critical density `|u|^3`.

A separate time-derivative sensitivity compares the quadratic checkpoint
derivative with the nearest two-point secant.  This is not an error bound, but
it prevents sparse checkpoint output from being mistaken for a resolved PDE
balance.

The pre-registered promotion gate requires

- normalized residual at most `0.10`;
- time-derivative sensitivity at most `0.20`;
- a non-unresolved balance classification.

## Test result

The new module and the previously imported critical-L3, scaling and free-space
Poisson modules pass the focused suite:

```text
44 passed
```

The full historical repository suite was started.  It produced no failures
before the execution environment timed out near the middle, so this report does
not claim full-suite completion.

## Application to the stored Hou early-time checkpoints

Fifteen checkpoints were analyzed: five times at each of `65x128`, `129x256`
and `193x384`.

### Resolution-stable observation in the middle of the early interval

At `t=0.0005` and `t=0.001`, all three resolutions classify both transformed
equations as a time--inviscid balance **before** applying the strict promotion
gate.

At `t=0.001`, the finest-grid critical-core values are:

| quantity | `u1` equation | `omega1` equation |
|---|---:|---:|
| normalized checkpoint residual | 0.0622 | 0.0668 |
| time-derivative sensitivity | 0.430 | 0.621 |
| `||advection-source||/(||advection||+||source||)` | 0.352 | 0.243 |
| `||diffusion||/(||advection||+||source||)` | 2.85e-4 | 3.85e-4 |

The corresponding ratios differ only slightly across the three grids.  Thus,
within the represented discrete finite-cylinder model, the early middle window
is overwhelmingly inviscid in term magnitude and contains substantial
advection/source cancellation.

This is consistent with the previously observed fact that viscosity is too
small in this window to identify its sign from the Hou run itself.  It does not
show an Euler singularity or a Navier–Stokes singularity.

### Negative result: no balance claim passes the promotion gate

Across the four positive-time checkpoints and two equations at each of the
three resolutions, the strict gate passes

```text
0 / 8 equation-time cases per resolution.
```

The main reason is sparse-time differentiation.  The derivative-sensitivity
ratio is already about `0.31--0.58` at `t=0.0005`; it reaches about
`0.69--0.96` around `t=0.0015`.  At the final stored time, normalized residuals
rise to about

- `0.318` for the finest-grid core `u1` equation;
- `0.444` for the finest-grid core `omega1` equation.

Therefore the later part of the stored early series cannot distinguish a true
change of PDE balance from checkpoint-interpolation error.

### Relation to the previous critical-L3 and scaling results

The previous diagnostics found, over this same early interval,

- decreasing represented-domain physical `L3` norm;
- increasing rather than shrinking critical-density RMS widths;
- no accepted Type-II or anisotropic power-law fit.

The new result adds:

- the term magnitudes are initially inviscid dominated;
- viscosity is not the missing mechanism in the early window;
- the existing checkpoint cadence is too sparse to classify the later balance.

Thus the early Hou window is not promoted to any of the three surviving
single-core scaling families.

## Consequence for the next run

The next nonlinear experiment should not merely extend the final time.  It must
write sufficiently dense checkpoints, or compute independent stage-level time
derivatives, so that the term-balance gate can pass or fail for mathematical
rather than sampling reasons.

Minimum required additions:

1. Record dense checkpoints with a cadence selected by a derivative-convergence
   study, not by output convenience.
2. Save the solver-stage RHS and also estimate the time derivative independently
   from adjacent accepted states.
3. Require agreement of those two paths under time-step refinement.
4. Measure the three PDE blocks in the same physical core used by the critical
   `L3` and shell diagnostics.
5. Apply the analysis after the viscosity switch and in the nonperiodic/free-
   space setup.
6. Promote an Euler-type, triple, or quasi-steady-viscous balance only when the
   residual, derivative sensitivity, spatial resolution and cross-integrator
   gates pass simultaneously.

## Files

- `src/ns_certificate_lab/term_balance.py`
- `tests/test_term_balance.py`
- `experiments/analyze_term_balance.py`
- `outputs/term_balance_old_snapshot/term_balance.csv`
- `outputs/term_balance_old_snapshot/term_balance.json`
- `outputs/term_balance_old_snapshot/manifest.json`

## Current conclusion

The stored early Hou data show a resolution-stable, strongly inviscid term
magnitude hierarchy in the middle of the interval, but the time sampling is too
sparse to validate the balance in the later interval.  Combined with decreasing
critical `L3` and expanding widths, this early window is not a whole-space
singularity candidate.
