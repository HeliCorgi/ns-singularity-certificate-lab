# Numerical foundation

## What is implemented

The initial code is an audit layer, not a production Navier--Stokes solver.
Arrays have shape `(n_r, n_z)` on a uniform grid containing the axis \(r=0\).
The default axial grid is periodic and omits the duplicated right endpoint.

`src/ns_certificate_lab/operators.py` supplies second-order finite differences:

- \(\partial_r,\partial_{rr},\partial_z,\partial_{zz}\);
- the scalar expression
  \(L_5=\partial_{rr}+3r^{-1}\partial_r+\partial_{zz}\);
- recovery of \(u^r,u^z\) from \(\psi_1\);
- the physical three-dimensional cylindrical divergence; and
- the defect in \(-L_5\psi_1=\omega_1\).

`src/ns_certificate_lab/cartesian.py` supplies a deliberately separate audit
route.  It uses independent fourth-order first-derivative stencils to recover
the meridional velocity, maps

\[
(u^r,u^\theta,u^z)\mapsto
(u_x,u_y,u_z)
\]

on user-selected azimuthal angles, projects the stored Cartesian components
back to the radial direction, and checks

\[
\omega^\theta=\partial_z u^r-\partial_r u^z=r\omega_1.
\]

This catches component-map, curl-sign, and elliptic-sign errors without
calling the production \(L_5\) routine.  It is not yet a finite-difference
calculation on a uniform \(x,y,z\) grid and does not evaluate all primitive
Cartesian Navier--Stokes residual components.

The name \(L_5\) is shorthand for a scalar differential expression.  It is
never used as an incompressibility condition.  The physical check is

\[
\partial_r u^r+\frac{u^r}{r}+\partial_z u^z=0.
\]

At the axis, even scalar fields use

\[
\lim_{r\to0}\left(f_{rr}+\frac3r f_r\right)=4f_{rr}(0),
\]

and smooth odd \(u^r\) uses

\[
\lim_{r\to0}\left(\partial_r u^r+\frac{u^r}{r}\right)
=2\partial_r u^r(0).
\]

The outer radial line uses one-sided diagnostic stencils.  No physical
outer-boundary condition is silently imposed by these operator functions.

## PDE residual decomposition

`compute_pde_residuals` accepts the fields and independently supplied time
derivatives.  It returns, without combining away the evidence:

- time derivative;
- radial advection;
- axial advection;
- right-hand-side source;
- viscous contribution; and
- signed total residual.

For each equation the convention is

\[
R=(\text{time})+(\text{radial advection})+(\text{axial advection})
-(\text{source})-(\text{diffusion}).
\]

This makes a sign or coefficient failure attributable to an individual term.
The routine does not estimate time derivatives internally, so a future time
integrator can be audited independently.

## Independent manufactured oracle

`src/ns_certificate_lab/manufactured.py` contains closed-form polynomials and
trigonometric fields.  Every derivative used as an oracle is written
analytically and does not call the finite-difference implementation.  The
fields satisfy the elliptic relation exactly but solve a **forced** transformed
PDE; the exact forcing is evaluated from the closed forms.

The refinement ladder \(17\times32\), \(33\times64\), \(65\times128\) checks:

- first and second derivatives;
- velocity recovery;
- physical divergence;
- elliptic compatibility;
- stored-field Cartesian reconstruction and independent azimuthal curl;
- both forced PDE residuals; and
- axis parity.

## Non-singular evolution control

The control in `experiments/run_baseline.py` uses a separate tridiagonal
Crank--Nicolson implementation; it does not reuse the package operator.  Its
derivation and acceptance checks are in [benchmark.md](benchmark.md).

## Deliberate current limitations

- There is no nonlinear production time integrator.
- There is no elliptic solver; the initial milestone checks an explicit
  \((\psi_1,\omega_1)\) pair.
- There is no independent full \(x,y,z\) divergence, full curl, or primitive
  Cartesian momentum-residual implementation; the present independent curl
  check operates on the cylindrical sampling after component reconstruction.
- There is no adaptive mesh, spectral solver, neural network, dynamic
  rescaling, interval arithmetic, or candidate search.
- Second-order consistency and checksum integrity are not truncation-error
  bounds, interval enclosures, or proofs.
- A one-sided finite-domain stencil does not certify a whole-space boundary
  condition.

These omissions are intentional gates: candidate discovery begins only after
the audited foundation has an independent reproduction.
