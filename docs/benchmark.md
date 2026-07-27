# Non-singular baseline experiment

## Purpose

This experiment is a negative control.  It verifies stability, refinement,
energy diagnostics, outer-boundary sensitivity, and a minimal false-positive
guard on a field known analytically to remain smooth.  It is not a candidate
search.

## Exact field

On a cylinder periodic in \(z\), take

\[
u^r=u^z=0,\qquad u^\theta(r,t)=r u_1(r,t),\qquad
\omega_1=\omega^\theta/r=0 .
\]

The radial pressure gradient balances the centripetal acceleration.  The
azimuthal equation reduces to

\[
\partial_tu_1=\nu\left(\partial_{rr}+\frac3r\partial_r\right)u_1.
\]

For \(a>0,\sigma>0\), the smooth Gaussian solution is

\[
u_1(r,t)=a\left(\frac{\sigma^2}{\sigma^2+4\nu t}\right)^2
\exp\left[-\frac{r^2}{\sigma^2+4\nu t}\right].
\]

The factor with exponent two is the heat-kernel scaling for four radial
dimensions.  This use of the formal four-radial-dimensional heat operator is
only an algebraic identity for \(u_1\); physical divergence and energy remain
three-dimensional.

The exact field is regular at the axis because \(u_1\) is even in \(r\), hence
\(u^\theta=r u_1\) is odd.  The energy per unit axial length is

\[
\frac{E}{L_z}=\pi\int_0^{R} r^3u_1(r,t)^2\,dr.
\]

The periodic cylinder has finite total energy for its finite axial period.
However, the same \(z\)-independent field extended to all of
\(\mathbb R^3\) has infinite total energy.  Consequently this experiment is a
local numerical control for the transformed operator and diagnostics, not a
finite-energy whole-space benchmark for the main Clay problem.

The recorded control uses \(L_0=U_0=1\), time scale \(L_0/U_0=1\),
\(\nu=0.05\), hence \(Re=U_0L_0/\nu=20\), and a documented periodic axial
length \(L_z=2\pi\).  Energy is reported per unit dimensionless \(z\)-length.

## Independent discretization

`experiments/run_baseline.py` uses a tridiagonal second-order radial stencil
and Crank--Nicolson time stepping.  It intentionally does not call the
differentiation routines in `src/ns_certificate_lab/operators.py`.  At the
axis it uses the analytic even-field limit

\[
\left(\partial_{rr}+3r^{-1}\partial_r\right)f(0)=4f_{rr}(0).
\]

The convergence runs use the analytic outer Dirichlet trace.  A separate pair
of runs imposes homogeneous outer data at two radii and compares the common
interior.  This separates discretization convergence from finite-domain
sensitivity.

## Recorded rejection logic

The diagnostic refuses even to attempt a reciprocal blow-up-time fit unless
the resolved physical-vorticity peak grows monotonically and materially over a
trailing window.  The Gaussian control should fail that growth gate and be
classified as bounded or decaying.  This guard can catch a simple false
positive; it cannot prove regularity and is not sufficient to validate a
future singularity candidate.

## Fixed-grid temporal convergence

`experiments/run_time_convergence.py` isolates the Crank--Nicolson time error
from a changing spatial grid.  The radial grid is fixed at
\(n_r=513\), \(R=5\), and the final time is \(1\); the requested and actual
time steps are exactly \(0.5,0.25,0.125\), corresponding to \(2,4,8\) steps.
The analytic outer trace is used for the convergence comparison.

The weighted relative \(L^2\) norm uses the physical swirl-energy weight
\(r^3\,dr\).  The recorded values are:

| \(\Delta t\) | relative \(L^2\) error | maximum absolute error | final energy per unit \(z\) | final max vorticity | auxiliary inner boundary sensitivity |
|---:|---:|---:|---:|---:|---:|
| 0.5 | \(8.586076\times10^{-4}\) | \(1.626945\times10^{-3}\) | 0.2725504862 | 1.3856349981 | \(2.155502\times10^{-8}\) |
| 0.25 | \(2.048943\times10^{-4}\) | \(3.882867\times10^{-4}\) | 0.2726696354 | 1.3881123155 | \(5.382548\times10^{-9}\) |
| 0.125 | \(4.382275\times10^{-5}\) | \(8.493974\times10^{-5}\) | 0.2726992741 | 1.3887190094 | \(2.031933\times10^{-9}\) |

The adjacent analytic-error orders are \(2.067119\) and \(2.225128\).
On the same grid,

\[
\frac{\lVert u_{\Delta t}-u_{\Delta t/2}\rVert}
     {\lVert u_{\Delta t/2}-u_{\Delta t/4}\rVert}
\]

gives order \(2.020029\); the corresponding normalized differences are
\(6.537258\times10^{-4}\) and \(1.611782\times10^{-4}\).
Raw analytic errors still contain the fixed spatial discretization error and
can eventually plateau, so the raw order is not interpreted as a standalone
continuum result.

The initial energy per unit \(z\) is \(0.3926990819\) for every run and
decreases in every recorded history.  The maximum physical-vorticity value
over each complete history is its initial value \(2.0\); there is no growth
signal.  Auxiliary boundary sensitivity is measured independently with
homogeneous outer traces on radii \(3\) and \(4\), compared only inside
\(r\leq1.5\).  The all-recorded-time maximum happens at the final time here,
so both saved metrics have the tabulated value.  This auxiliary comparison
does not directly bound the boundary error of the main \(R=5\) run.
These checks establish the expected finite-grid behavior of this smooth
control only.  They neither prove regularity of general solutions nor bear
evidence of a singularity.

## Run

```console
python experiments/run_baseline.py --config configs/baseline.json --output-dir outputs/baseline_replay
python -m experiments.run_time_convergence --config configs/baseline_time_convergence.json --output-dir outputs/time_convergence_replay
```

The script refuses to overwrite a non-empty evidence directory.  Use a new
`--output-dir` for a replay.
