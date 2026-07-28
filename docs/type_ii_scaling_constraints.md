# Conditional Type-II and anisotropic scaling constraints

## Scope

This document does not prove a Navier--Stokes singularity.  It classifies a
restricted power-law model and identifies exponent ranges that are not
immediately excluded by finite energy, finite energy dissipation, critical
Serrin scaling, and a generic leading-order PDE balance.

Let

\[
\tau=T-t,\qquad
A(t)\asymp \tau^{-\alpha},\qquad
L_r(t)\asymp\tau^{\beta_r},\qquad
L_z(t)\asymp\tau^{\beta_z}.
\]

Define

\[
B=2\beta_r+\beta_z,
\qquad
\gamma=\max(\beta_r,\beta_z).
\]

Here `B` is the modeled core-volume exponent and `gamma` is the thinnest
length-scale exponent.  The derivation assumes a nondegenerate localized
profile.  In particular, it assumes nonzero profile mass and a nonzero
derivative in a direction with scale `tau^gamma`.

## Necessary margins

The modeled core energy scales as

\[
E_{\rm core}\asymp A^2L_r^2L_z\asymp\tau^{B-2\alpha}.
\]

Finite global energy therefore requires

\[
B\ge 2\alpha.
\]

The modeled critical norm scales as

\[
\|u\|_3^3\asymp A^3L_r^2L_z\asymp\tau^{B-3\alpha}.
\]

A core-driven escape from the endpoint `L^3` obstruction requires

\[
B<3\alpha.
\]

The classical critical Serrin pair with `p=infinity`, `q=2` rules out
`alpha<1/2` in this model.  More generally,

\[
\|u\|_p\asymp\tau^{-\alpha+B/p},\qquad
q=\frac{2p}{p-3}.
\]

Avoiding every classical critical Serrin class requires

\[
\alpha\ge\frac12,\qquad B\le3\alpha.
\]

Under the nondegenerate thinnest-gradient assumption,

\[
\|\nabla u\|_2^2
\asymp
A^2L_r^2L_zL_{\min}^{-2}
\asymp
\tau^{-2\alpha-2\gamma+B}.
\]

Finite total energy dissipation requires

\[
2\alpha+2\gamma-B<1.
\]

Equality produces a logarithmically divergent dissipation integral and is not
accepted.

## Generic term balance

Using the thinnest effective derivative scale, the modeled singular exponents
of the PDE blocks are

\[
E_t=\alpha+1,
\qquad
E_I=2\alpha+\gamma,
\qquad
E_\nu=\alpha+2\gamma.
\]

`E_I` represents the inertial-pressure block.  If exactly one block has the
largest exponent, that leading term must vanish or cancel for structural
reasons.  Without such a measured cancellation, at least two leading exponents
must agree.

This gives three generic families.

### Family I: Euler-like Type-II core

For `gamma < 1/2`, time and inertia can balance:

\[
\alpha+\gamma=1.
\]

Combining energy, geometry and critical growth yields

\[
\boxed{
\frac25\le\gamma<\frac12,
\qquad
\alpha=1-\gamma,
\qquad
2-2\gamma\le B\le3\gamma.
}
\]

For isotropic scales, `B=3 gamma`, so

\[
\boxed{
\frac25\le\beta<\frac12,
\qquad
\alpha=1-\beta.
}
\]

The viscous term is lower order.  A candidate in this family must therefore
exhibit an Euler-like inner balance plus a valid viscous/outer matching.

### Family II: anisotropic parabolic core

At

\[
\alpha=\gamma=\frac12,
\]

all three modeled blocks balance.  Critical growth is possible only through
anisotropic volume collapse:

\[
\boxed{1<B<\frac32.}
\]

The isotropic endpoint `B=3/2` has bounded critical `L^3` and is excluded.

### Family III: quasi-steady viscous-inertial core

For `gamma > 1/2`, inertia and viscosity can balance with time lower order:

\[
\alpha=\gamma.
\]

The remaining constraints give

\[
\boxed{
\frac12<\gamma<1,
\qquad
\alpha=\gamma,
\qquad
4\gamma-1<B<3\gamma.
}
\]

This family must exhibit a resolved near-cancellation of the inertial-pressure
and viscous blocks while the time derivative is asymptotically smaller.

## Escape routes not covered

The classification is conditional and does not exclude:

- leading convection-pressure cancellation;
- different derivative scales for different velocity components;
- multiple interacting cores;
- a shell cascade rather than one localized profile;
- global `L^3` growth carried by an outer field;
- logarithmic corrections;
- non-power-law or intermittent dynamics.

Such behavior must be measured explicitly rather than inferred from a maximum
vorticity curve.

## Required numerical observables

Future candidate runs should fit and cross-check:

- amplitude exponent `alpha` from several independent amplitudes;
- `beta_r` and `beta_z` from FWHM, quantile width and derivative scales;
- `B=2 beta_r+beta_z`;
- global and shell-decomposed `L^3` mass;
- energy and dissipation scaling;
- each PDE block norm;
- inertial-pressure cancellation ratio;
- fit sensitivity to the proposed singular time and fitting window.

No exponent fit is admissible until the relevant scale is resolved, the
space/time/domain convergence gates pass, and multiple definitions of the
scale agree.

## Machine implementation

- `src/ns_certificate_lab/scaling_constraints.py`
- `tests/test_scaling_constraints.py`
- `experiments/scan_scaling_constraints.py`

The scanner enumerates conditional exponent triples.  Its output is a search
prior, not a list of Navier--Stokes solutions.
