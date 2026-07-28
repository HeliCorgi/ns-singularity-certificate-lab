# Type-II / anisotropic Navier--Stokes search progress

## What was added

- conditional exponent-constraint engine;
- classification of generic localized power-law cores;
- candidate singular-time and fit-window scanner;
- tests using exact synthetic power laws;
- application to the available early Hou critical-`L^3` snapshots.

## Conditional generic families

Let

\[
A\sim(T-t)^{-\alpha},\quad
L_r\sim(T-t)^{\beta_r},\quad
L_z\sim(T-t)^{\beta_z},
\]

\[
B=2\beta_r+\beta_z,\qquad
\gamma=\max(\beta_r,\beta_z).
\]

Under finite energy, finite energy dissipation, critical-norm growth, one
effective thinnest derivative scale and no leading-order cancellation, the
surviving families are:

1. Euler-like Type-II:
   \[
   2/5\le\gamma<1/2,\quad
   \alpha=1-\gamma,\quad
   2-2\gamma\le B\le3\gamma.
   \]

2. Anisotropic parabolic:
   \[
   \alpha=\gamma=1/2,\quad 1<B<3/2.
   \]

3. Quasi-steady viscous-inertial:
   \[
   1/2<\gamma<1,\quad
   \alpha=\gamma,\quad
   4\gamma-1<B<3\gamma.
   \]

The isotropic Type-II subcase is

\[
2/5\le\beta<1/2,\qquad \alpha=1-\beta.
\]

These are conditional search regions, not solutions.

## Early Hou result

Scanning candidate singular times and all suffix windows with at least three
points produced zero gate-passing fits at all three resolutions.

| grid | passing fits | best alpha | best beta_r | best beta_z |
|---|---:|---:|---:|---:|
| 65 x 128 | 0 | 0.0743 | -0.0125 | -0.0101 |
| 129 x 256 | 0 | 0.1250 | -0.0116 | -0.00865 |
| 193 x 384 | 0 | 0.1390 | -0.0114 | -0.00824 |

The critical-density widths expand during the available early interval. The
early observed amplitude growth is therefore not a Type-II, anisotropic or
multiscale singular-core candidate under this gate.

## Test status

The new focused suite passed:

```text
23 passed
```

The full historical suite was started but exceeded the execution timeout, so a
full-suite success is not claimed here.

## Next scientific step

Collect later, resolved whole-space-compatible snapshots and measure the three
PDE blocks and inertial-pressure cancellation. Family I has an Euler-like
leading profile and is subject to known nonexistence results under additional
decay/integrability assumptions; those assumptions must be audited before it
can be used as the main target.
