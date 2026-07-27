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
calling the production \(L_5\) routine.

`src/ns_certificate_lab/cartesian_validation.py` is a second, more independent
route on a uniform three-dimensional \((x,y,z)\) grid.  Its local second-order
stencils compute

\[
\nabla\cdot u,\qquad \nabla\times u,\qquad \Delta u,
\]

the three Cartesian pressure-gradient components, the advective term, the
unforced primitive residual, and its forced defect

\[
R_0=u_t+(u\cdot\nabla)u+\nabla p-\nu\Delta u,\qquad
R_f=R_0-f.
\]

It does not import the cylindrical grid, `operators.py`, `pde.py`, or the older
Cartesian audit.  Periodic manufactured fields test every component against
closed-form derivatives over three refinements.

For the artifact-level path,
`src/ns_certificate_lab/cartesian_candidate_adapter.py` accepts a candidate
after `save_candidate`/`load_candidate`.  It differentiates the stored
\(u_1,\psi_1\) arrays with stencils local to that adapter, bilinearly samples
them at \(r=\sqrt{x^2+y^2}\), and reconstructs

\[
\begin{aligned}
u_x&=-x\psi_{1,z}-y u_1,&
u_y&=-y\psi_{1,z}+x u_1,&
u_z&=2\psi_1+r\psi_{1,r},\\
\omega_x&=-x u_{1,z}-y\omega_1,&
\omega_y&=-y u_{1,z}+x\omega_1,&
\omega_z&=2u_1+r u_{1,r}.
\end{aligned}
\]

The uniform-grid checker then recomputes divergence and full curl without
consulting any cylindrical derivative result.  Thus the audit is independent
of the production differentiation path, although it remains a finite
second-order calculation with interpolation error, not an exact continuum
identity.

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

A separate periodic Cartesian manufactured field on \(12^3,24^3,48^3\)
checks three-component divergence, full curl, vector Laplacian, advection,
pressure gradient, viscous contribution, and the forced primitive residual.
Every reported error decreases at approximately second order.  The saved
axisymmetric-candidate route is also exercised after an actual artifact
round-trip; injected radial-sign, component-map, vorticity-sign, and
divergence faults must be rejected.

## Non-singular evolution control

The control in `experiments/run_baseline.py` uses a separate tridiagonal
Crank--Nicolson implementation; it does not reuse the package operator.  Its
derivation and acceptance checks are in [benchmark.md](benchmark.md).

`experiments/run_time_convergence.py` holds that spatial grid fixed at
\(n_r=513\) and changes only the requested time step through
\(0.5,0.25,0.125\).  It reports both raw analytic-solution errors and the
step-doubling difference
\(\lVert u_{\Delta t}-u_{\Delta t/2}\rVert\).  Raw errors retain the fixed
spatial error; the difference ratio is included because the common leading
spatial error largely cancels.

## Legacy Poisson prototype decision

Only the prescribed files in the old prototype archive were inspected
read-only.  Its periodic-\(z\), outer-Dirichlet finite-cylinder matrix contains
useful axis-row and manufactured-solution design ideas, but depends on the old
grid/operator conventions and leaves conditioning and several boundary tests
unverified.  No code was copied.  The detailed audit and future acceptance
criteria are in [legacy_reuse_review.md](legacy_reuse_review.md).

## Deliberate current limitations

- There is no nonlinear production time integrator.
- There is no elliptic solver; the initial milestone checks an explicit
  \((\psi_1,\omega_1)\) pair.  The next minimum gate is a newly written
  independent finite-cylinder Poisson solver, not a copy of the legacy code.
- The uniform Cartesian audit uses finite differences and bilinear candidate
  interpolation on a bounded box; it is not a truncation-error enclosure.
- Primitive pressure and forcing are supplied to the residual checker.  There
  is not yet an independent pressure recovery or projection solver.
- There is no adaptive mesh, spectral solver, neural network, dynamic
  rescaling, interval arithmetic, or candidate search.
- Second-order consistency and checksum integrity are not truncation-error
  bounds, interval enclosures, or proofs.
- A one-sided finite-domain stencil does not certify a whole-space boundary
  condition.

These omissions are intentional gates: candidate discovery begins only after
the audited foundation has an independent reproduction.
