# Research plan

## Scope of this repository

The present milestone builds an audited, reproducible foundation and runs one
deliberately non-singular benchmark.  It does **not** launch a large candidate
search and does not claim a finite-time singularity.

## Phase 0 — audited equations (current)

- Fix the physical domain, conventions, dimensions, and regularity class.
- Derive the axisymmetric-with-swirl system from the three-dimensional PDE.
- Audit every sign, radial coefficient, recovery formula, and axis condition.
- Map known exclusion and continuation theorems onto the proposed search space.

Exit gate: no implemented formula has status `未確認`, `不整合`, or `誤り` in
`docs/equation_audit.md`.

## Phase 1 — independently testable numerical primitives (current)

- Implement finite-difference operators on a cell-inclusive cylindrical grid.
- Recover meridional velocity from the stream function.
- Compute physical divergence, the elliptic defect, and each PDE term
  separately.
- Enforce explicit candidate schemas, hashes, configurations, and seeds.
- Verify with analytic manufactured fields and intentional fault injection.

Exit gate: the full test suite passes, observed refinement errors decrease at
the documented rate, and corrupt artifacts are rejected.

## Phase 2 — non-singular controls (current)

- Evolve a smooth, exactly controlled axisymmetric field on modest CPU grids.
- Check stability, refinement, energy behavior, boundary sensitivity, and the
  distinction between a steep transient and sustained norm growth.
- Store data and plots with a machine-readable run summary.

Exit gate: replay from a fresh environment gives the same diagnostics within
declared tolerances.

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
