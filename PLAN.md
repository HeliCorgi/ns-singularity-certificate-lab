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

## Phase 2.5 — independent elliptic solve gate (completed)

- Two independent implementations now exist: the committed `poisson.py`
  (\(r^3\)-flux finite volume, Fourier in \(z\), Thomas solve) and the
  integrated bundle solver `finite_cylinder_poisson.py` (non-divergence-form
  direct differences, Fourier in \(z\), independent Thomas solve).
- Both pass non-circular manufactured tests at three refinements with
  observed order ~2, and both reject sign/axis/boundary faults.
- A dedicated cross-validation suite (`tests/test_poisson_cross_validation.py`)
  pins the shared conventions to roundoff on radially-exact fields (CV-2),
  measures the genuine \(O(\Delta r^2)\) inter-solver disagreement with an
  explicit lower bound that fails if independence ever collapses (CV-1), and
  runs paired fault injections (CV-3).
- Documented limitation: the axial Fourier treatment and the grid class are
  shared between the two solvers, so independence holds for the radial
  operator, boundary bookkeeping, and linear algebra only.  A non-Fourier
  axial path remains an open verification gate.

Status: reached.  Conditioning and boundary/truncation limitations are
recorded; no whole-space claim is made.

## Phase 2.6 — nonlinear wall-bounded production solver (completed)

- Hou finite-cylinder setup audited from the arXiv LaTeX sources
  (v1, v2, and the deferred method paper arXiv:2102.06663) into
  `docs/hou_setup_audit.md` and equation-audit entries E-27–E-31.
- Production Heun/RK2 integrator for the full nonlinear
  \((u_1,\omega_1,\psi_1)\) system with no-slip/no-flow wall conditions
  (single Dirichlet \(\psi_1(1,z)=0\) elliptic solve plus second-order
  Thom-type wall vorticity, E-31), adaptive CFL time step, and the
  two-stage viscosity protocol (E-30).
- Verified by forced manufactured convergence (space and time, order ~2),
  zero-field/small-amplitude/symmetry/circulation invariants, five fault
  injections, and restart fidelity.

Status: reached.  This validates the discretization, not any physical claim.

## Phase 2.7 — Hou early-time run on uniform grids (evidence stage)

- Integrate the audited E-29 datum at three uniform resolutions to
  \(t=T_1=0.002191729\) with \(\nu=5\times10^{-4}\), full-period \(z\),
  monitored (not imposed) odd symmetry.
- Record amplification trajectories, \((R(t),Z(t))\), energy/enstrophy,
  circulation, divergence residuals, and the independent solver-B elliptic
  cross-check at every snapshot.
- Compare against the published 1536² amplification 20.5235 at \(T_1\)
  while stating plainly that a uniform fixed grid cannot resolve Hou's
  adaptive scales; discrepancies are evidence, not failure to be hidden.

Exit gate: multi-resolution trends recorded with raw checkpoints; either a
resolution-consistent growth trend or a documented negative result.  Do not
start candidate discovery from this phase alone.

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
