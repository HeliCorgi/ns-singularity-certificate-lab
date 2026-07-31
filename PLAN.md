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

Status amendment (2026-07-28): the shipped Heun runs of this phase are
reclassified **stability-unverified** by the frozen-coefficient von Neumann
audit (`docs/numerical_stability_audit.md`); their amplification numbers must
not feed candidate decisions unless corroborated by the Phase 2.8
cross-integrator comparison, and no convergence fit may be attempted while
the `core_width.fit_precondition` gate fails on their snapshots.

## Phase 2.8 — numerical-stability and resolution gates (P0, 2026-07-28)

Mandated by FABLE5_NEXT_TASK_AUDIT.md before any mid/late-time Hou run,
blow-up-time fit, or AI candidate search:

- **Gate 1 (stability)**: frozen-coefficient von Neumann audit of Heun with
  centered differences (`von_neumann.py`); pre/predictor/post CFL recorded
  per accepted step with optional predictor-stage rejection; SSPRK3/RK4
  cross-checking integrators sharing the identical spatial discretization
  (`run_integrator_comparison.py`).  Heun-only amplification is barred from
  candidate decisions.
- **Gate 2 (time convergence)**: dt, dt/2, dt/4 at fixed grid
  (`outputs/hou_time_refinement_v1`, to be re-examined under the new
  full-step streaming gates).
- **Gate 3 (space convergence)**: resolution ladder with common-grid profile
  differences and the preregistered points-per-scale precondition
  (`core_width.py`, minimum 7 points per 10--90 front) before any fit;
  blind extrapolation only (`extrapolation.py`, no external anchors) — the
  current ladder is *not in its asymptotic range*.
- **Gate 4 (whole-space transition)**: non-periodic \(z\), compact support
  in \(r\) and \(z\), a free-space elliptic path, independent
  \(R_{\max}\)/\(Z_{\max}\) enlargement, low-wavenumber stress tests
  (`docs/whole_space_transition.md` §7).  Not yet implemented; the W-A
  transparent radial condition covers the radial component only, and all
  existing wall results are periodic-z radial-wall sensitivity observations.

Exit gate: Gates 1–3 measured and recorded (pass or documented failure);
Gate 4 implemented and passed.  Until then the Hou mechanism is not called
an \(\mathbb R^3\) candidate.

## Phase 2.85 — Gate 4, the linear whole-space elliptic gate (completed, 2026-07-29)

Gate 4 of Phase 2.8 is implemented in `src/ns_certificate_lab/whole_space_gate.py`
and passed all twenty preregistered acceptance checks in
`outputs/whole_space_gate4_v1`.

- Non-periodic `z` on a finite box, with the axial second difference
  diagonalised by a dense **discrete sine transform built from `numpy.sin`** —
  the non-Fourier independent axial path this repository has listed as open
  since 2026-07-28.
- An **exact closed-form free-space reference**: `L5` is the axisymmetric
  five-dimensional Laplacian, so Newton's theorem solves a compactly supported
  radial bump in closed form, and `z`-translation invariance makes superpositions
  exact too.  No quadrature, no discretisation on the reference side.
- Independent `R_max` and `Z_max` enlargement, with the truncation component
  isolated by common-interior differencing because the direct comparison
  saturates at the discretisation floor.
- The periodic zero axial mode is shown to overstate the far field by exactly
  `2R/L`, and the non-periodic Dirichlet operator has **no zero mode at all**.
- An a posteriori bound on the continuum truncation error, dominating every
  measured monopole row.
- An independent Cartesian audit of the recovered velocity: `div u = 0` and
  `curl u = omega1 (-y, x, 0)`, both second order.

Exit gate: reached, with the scope stated precisely: **the gate certifies the
free-space potential *value* only.**  The derivative tail needed for velocity
recovery and the nonlinear coupling were not verified by this phase; they are
Phase 2.86.

## Phase 2.86 — Gate 5: derivative tails and the nonlinear entry point (completed, 2026-07-29)

Gate 4 certified the free-space potential *value*.  Phase 2.86 promotes that to a
differentiated entry point.

- **Explicit derivative tail bounds** from analytic differentiation of the
  five-dimensional Green kernel: `|D^k(psi - psi_multipole)| <= A_{k+1} I_1 /
  d^{4+k}` for the monopole truncation and `<= (1/2) A_{k+2} I_2 / d^{5+k}` for
  the dipole one, with `A_m` the exact homogeneity constants of `D^m G_5`.  No
  maximum principle is used; the interior propagation uses the classical
  interior elliptic estimate `(nk/rho)^k` instead, which the maximum principle
  cannot supply for `k >= 1`.
- **Two independent evaluation paths**: finite differences of the solver output
  (the production path) and an analytic-Green quadrature that differentiates the
  kernel and never touches the solver.  The second converges at order 4 and
  agrees with the closed-form reference to `1.6e-8`.
- **A free-space velocity recovery API** returning `psi1`, its first and second
  derivatives, `u^r`, `u^z`, and both the boundary and interior derivative tail
  bounds, all second order and with exact axis regularity.
- **A small-amplitude nonlinear whole-space run** from smooth compactly
  supported divergence-free pure-swirl data, refined independently in `dt`,
  `dr`, `dz`, `R_max`, `Z_max`, integrator and outer boundary order.
- **A finite dyadic cascade model** showing that low-mode-only smooth forcing
  drives high-shell amplitude by 26 orders of magnitude with exactly zero direct
  injection, so the earlier reading "forcing gives Track F no advantage" was
  wrong and has been corrected everywhere.

Exit gate: reached, with three sub-gates recorded as **uninformative rather than
passed** — at this amplitude the time refinement, the domain enlargement and the
boundary-order comparison all produce identical answers, because the temporal
error is far below the spatial one and the field never reaches the boundary.
They must be repeated at an amplitude that exercises them.

## Phase 2.87 — Gate 6: mid-amplitude calibration and continuation (2026-07-29)

- **Nonlinear tail propagation** (`tail_propagation.py`): explicit constants from
  the potential error bounds through the velocity error to every right-hand-side
  term, using only the product-difference identity, plus a short-time Gronwall
  step.  The Lean side is the certificate layer F-17/F-18/F-19.
- **An explicit initial-data family** (`initial_data.py`): smooth, compactly
  supported, axis-regular, exactly divergence-free pure swirl with the radial
  factor a function of `r^2`.  Fixed positive viscosity; no two-stage protocol.
- **Multipole hierarchy** to quadrupole, and a `tail_bound_available` flag for
  the boxes where no multipole bound exists at all.
- **Domain expansion** with preregistered triggers and before/after invariants.
- **An interval-arithmetic snapshot certificate** in exact rational arithmetic,
  with an independent checker: PO-05/06/07/13 are no longer design-only.

Two preregistered criteria **failed and were not retuned**:
the core boundary-condition difference is `7.9e-3` times the Richardson
discretisation error rather than the required `8`, and the amplitude
continuation never left the quadratic-response regime (`max|omega1|` scales as
`A^2` to better than `5e-5` across a factor of ten).  Both failures are
quantified in `outputs/whole_space_gate6_v1` and neither is hidden.

Status: the entry point is instrumented; **no candidate was promoted**, and the
three best were rejected with named reasons.

## Phase 2.9 — Track F fixed-finite-mode class closed (completed, 2026-07-29)

`START_NEW_SESSION_NAVIER_STOKES.md` section 6 "candidate A" proposed a symbolic
search over low-order divergence-free Fourier ansaetze for a residual force that
stays smooth across a finite-time singularity (Clay statements (C)/(D)).  **That
search space is provably empty and the search was therefore not run.**

- `docs/research_notes/track_f_finite_mode_nogo.md` proves: if the velocity of a
  Track-F ansatz stays inside a fixed finite-dimensional divergence-free space
  and the residual force is in \(L^1_tL^2_x\), then the energy identity bounds
  every norm of the ansatz and the projected ODE continues it smoothly past the
  putative singular time.  Contrapositive: **a Clay (C)/(D) counterexample must
  have unbounded Fourier bandwidth as \(t	o T^-\).**
- The one algebraic input, \(\int u\cdot(u\cdot
abla)u = 0\), is verified in
  **exact integer arithmetic** for each scanned mode set by
  `src/ns_certificate_lab/galerkin_obstruction.py` (every monomial coefficient of
  the cubic form over \(\mathbb Z[i]\)); no floating point is involved.
- The a priori bound itself is proved in Lean 4 as **F-6**
  (`formal/NSSingularity/GalerkinNoBlowup.lean`, no `sorry`/`admit`/axiom).
- Two necessary conditions for *any* Track-F counterexample follow without
  finite-dimensionality: bounded energy and finite total dissipation up to `T`.
  Clay condition (7) is therefore automatic, and a design that expects the
  energy to diverge is wrong from the start.
- The evidence run is `outputs/track_f_finite_mode_scan_v1` (10 families, all
  rejected, all 11 preregistered acceptance checks passed).

Open remainder: the ODE continuation step is proved on paper only (**F-7**), and
nothing here constrains ansaetze whose bandwidth grows as \(t	o T^-\).

Status: reached.  This is an exclusion theorem for a search class, not a step
toward a singularity.

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

## Phase 2.88 — Gate 7: leaving the Picard region (2026-07-29)

The single parameter is `Re = A L^2 / nu`; amplitude, length and viscosity are
not three search directions and `deduplicate_settings` now enforces that.

Reclassifying the Gate 6 sweep in `(Re, aspect, c, tau)` showed why all
thirty-two runs stayed inside the first Picard iterate: the nominal amplitude
overstated the field's peak by about twenty-three, so the sweep reached
`tau <= 0.0233`, forty-three times short of the `tau = 1` it needed.

`PicardLadder` integrates levels 0, 1, 2 and the full solution together, with
the same integrator and the same accepted steps, so the distance to the iterates
is measured rather than inferred, and the solver's exact right-hand side is
stored at every accepted step rather than reconstructed from snapshots.

Eighteen runs (Re = 10..400 x families S, A, H) reached `tau = 1`.  The
preregistered departure gate passed on all nine checks.  **No candidate was
promoted:** the critical `L^3` norm decays in every run and the dyadic shell
count grows in sixteen of eighteen, i.e. the fields spread rather than
concentrate.  The one amendment — adding a 145x289 grid — changed no threshold
and is recorded with its reason in the config.

F-7c is closed in Lean by the direct route: mathlib's `IsPicardLindelof` was
already time-dependent, so the autonomisation on `E x R` was never needed.

The `H^s` derivation (`docs/research_notes/hs_error_propagation.md`) is
deliberately incomplete, with HS-1..HS-6 named.  HS-5 — discrete residual to
`||R||_{H^s}` — is the step that would make any certificate a statement about
the PDE, and it is not done.  Until it is, the `L^infty` certificate is not
presented as an unconditional PDE proof.

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
