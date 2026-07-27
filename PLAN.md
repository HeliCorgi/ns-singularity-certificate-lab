# Research plan

## Scope of this repository

The present milestone builds an audited, reproducible foundation and runs one
deliberately non-singular benchmark.  It does **not** launch a large candidate
search and does not claim a finite-time singularity.

## Phase 0 — audited equations (completed for the implemented scope)

- Fix the physical domain, conventions, dimensions, and regularity class.
- Derive the axisymmetric-with-swirl system from the three-dimensional PDE.
- Audit every sign, radial coefficient, recovery formula, and axis condition.
- Map known exclusion and continuation theorems onto the proposed search space.

Exit gate: no implemented formula has status `未確認`, `不整合`, or `誤り` in
`docs/equation_audit.md`.

Status: reached for the formulas currently implemented.  Adding a new solver
still requires a new equation-audit entry before that solver can become an
implementation premise.

## Phase 1 — independently testable numerical primitives (completed)

- Implement finite-difference operators on a cell-inclusive cylindrical grid.
- Recover meridional velocity from the stream function.
- Compute physical divergence, the elliptic defect, and each PDE term
  separately.
- Enforce explicit candidate schemas, hashes, configurations, and seeds.
- Verify with analytic manufactured fields and intentional fault injection.
- On a uniform Cartesian \((x,y,z)\) grid, independently compute the
  three-component divergence, full curl, vector Laplacian, and every primitive
  Navier--Stokes momentum-residual term without calling the cylindrical
  production operators.
- Save and reload an axisymmetric candidate, reconstruct its Cartesian
  velocity and E-18b vorticity through a separate stencil/interpolation
  adapter, and reject injected sign, component-map, vorticity, and divergence
  faults.

Exit gate: the full test suite passes, observed refinement errors decrease at
the documented rate, and corrupt artifacts are rejected.

Status: reached on the recorded local replay.  This is a consistency gate, not
a proof of any continuum claim.

## Phase 2 — non-singular controls (completed)

- Evolve a smooth, exactly controlled axisymmetric field on modest CPU grids.
- Check stability, refinement, energy behavior, boundary sensitivity, and the
  distinction between a steep transient and sustained norm growth.
- Store data and plots with a machine-readable run summary.
- Hold the \(513\)-point radial grid fixed and repeat the same smooth control at
  \(\Delta t,\Delta t/2,\Delta t/4\), recording analytic errors,
  step-doubling order, energy, physical-vorticity maxima, and boundary
  sensitivity.

Exit gate: replay from a fresh environment gives the same diagnostics within
declared tolerances.

Status: the recorded run passed its numerical acceptance checks; the CI
workflow now includes this replay.  It remains a finite-grid negative control.

## Phase 2.5 — independent elliptic solve gate (next minimum step)

- Design a finite-cylinder solver for
  \(-(\partial_{rr}+3r^{-1}\partial_r+\partial_{zz})\psi_1=\omega_1\)
  with periodic \(z\) and an explicitly declared outer radial boundary.
- Implement it anew in the current conventions; do not copy the audited legacy
  prototype.
- Directly test the axis row, global sign, outer boundary rows, and a
  manufactured solution over at least three refinements.
- Record conditioning, residual, and boundary/truncation limitations before
  using the solver in any trajectory code.

Exit gate: an implementation independent of the existing diagnostic
operators passes non-circular manufactured and fault-injection tests.  Do not
start candidate discovery at this phase.

## Phase 3 — preregistered candidate discovery (future)

- Implement dynamically rescaled Type-II-compatible coordinates with separate
  radial and axial scales.
- Explore equilibria, periodic orbits, and slowly drifting trajectories in
  rescaled time.
- If machine learning is used, export every retained field to explicit spectral
  coefficients and recompute all diagnostics without the network.
- Preregister rejection thresholds, refinement ladders, and all reported runs.

Exit gate: a candidate survives independent discretizations, precision levels,
domain sizes, and withheld diagnostic tests.  This remains numerical evidence.

## Phase 4 — validated numerics (future)

- Replace floating-point residual estimates with interval enclosures.
- Bound truncation, interpolation, nonlinear convolution, spectral tails, and
  domain truncation.
- Prove an orbit/profile and its local invariant structure with a posteriori
  estimates.

Exit gate: a reproducible computer-assisted theorem with all hypotheses mapped
to the original three-dimensional equations.

## Phase 5 — analytic bridge (future)

- Prove that smooth finite-energy initial data enter the validated neighborhood.
- Prove that rescaled time corresponds to finite physical time.
- Prove divergence of a regularity-controlling physical norm, rather than a
  coordinate artifact.
- Formalize the theorem and machine-checkable certificate where practical.

Only after all applicable obligations in `docs/proof_obligations.md` are
discharged could a singularity claim be considered.
