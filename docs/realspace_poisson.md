# Real-space `-L5` Poisson path (solver C)

## Purpose

Two solvers for the finite-cylinder problem E-25 already exist:

| | module | radial | axial | linear solve |
|---|---|---|---|---|
| A | `ns_certificate_lab.poisson` | `r^3`-flux finite volume (E-26a) | discrete Fourier transform | Thomas elimination |
| B | `ns_certificate_lab.finite_cylinder_poisson` | direct non-divergence stencil | discrete Fourier transform | Thomas elimination |
| C | `ns_certificate_lab.realspace_poisson` | `r^3`-flux finite volume (E-26a) | three-point periodic difference | preconditioned conjugate gradients |

A and B were built to cross-validate each other, and their **radial** stencils
really are different. Their **axial** treatment is not: both diagonalize `z`
with NumPy's discrete Fourier transform, both take the axial wavenumbers from
the same NumPy sample-frequency helper, and both invert the resulting radial
systems with a locally written Thomas elimination. A defect in the wavenumber
convention, in the Nyquist mode, in the real/imaginary bookkeeping, or in the
tridiagonal elimination itself would move A and B in exactly the same way and
would be invisible to `tests/test_poisson_cross_validation.py`. This limitation
is recorded in `docs/finite_cylinder_poisson.md`.

Solver C exists to break that common mode. It is not a better solver; on
grid-resolved axial modes it is strictly *less* accurate than A and B. It is a
solver whose failure modes are different, which is the only property that makes
it useful as an audit path for snapshots.

## Method

The problem is E-25 exactly as A and B solve it:

```
-(d_rr + (3/r) d_r + d_zz) psi1 = omega1,   0 <= r <= R,
z periodic on [0, L_z) with the right endpoint omitted,
psi1(R, z) = g(z),   d_r psi1(0, z) = 0  (evenness, E-16c).
```

**Radial.** The E-26a control-volume form, transcribed into
`realspace_poisson.py` directly from `docs/equation_audit.md` rather than
copied from `poisson.py`:

```
V_i = (r_{i+1/2}^4 - r_{i-1/2}^4) / 4,     F_i = r_{i+1/2}^3 / dr,
(L5_r^h psi)_i = [ F_i (psi_{i+1} - psi_i) - F_{i-1} (psi_i - psi_{i-1}) ] / V_i,
```

with `F_{-1} = 0` at the axis. The axis cell has `r_{-1/2} = 0`,
`r_{1/2} = dr/2`, hence `V_0 = dr^4/64` and the row collapses to
`8 (psi_1 - psi_0) / dr^2`, which is E-26b. **The module does not hard-code
`8`.** It evaluates the general formula and then checks that
`F_0 / V_0 == 8/dr^2` to a relative tolerance of `1e-12`, raising
`ArithmeticError` otherwise, so the audited coefficient is derived rather than
asserted.

**Axial.** `(psi_{i,j+1} - 2 psi_{ij} + psi_{i,j-1}) / dz^2`, wrapped with
`numpy.roll`. No transform of any kind appears in the module; a test reads the
module source and asserts that the substring `fft` does not occur in it, because
a transform is reachable through the already-permitted `numpy` import and would
therefore be invisible to an import-graph check.

**Linear solve.** Matrix-free Jacobi-preconditioned conjugate gradients in the
weighted inner product below, with a relative-residual tolerance (default
`1e-12`), an iteration cap, and an explicit `ConvergenceError` when the cap is
reached. No tridiagonal elimination, no dense assembly, no SciPy. The stopping
test uses the cheap recursively updated residual, but every acceptance is
confirmed against a freshly recomputed `b - A x`; if they disagree the fresh
residual replaces the recursion and iteration continues. The value reported in
`weighted_relative_residual` is always the freshly recomputed one, so it cannot
be an artifact of accumulated drift.

## Symmetry

The outer row is prescribed, so the unknowns are rows `i = 0 .. n_r-2` and the
Dirichlet data moves to the right-hand side as `+ a_i^+ g` (E-26c). Define

```
<u, v>_w = sum_{i,j} w_i u_ij v_ij,     w_i = V_i * dz,
```

the Riemann sum of the `r^3 dr dz` measure over the control volumes. E-26c's
coefficients are `a_i^+ = F_i / V_i` and `a_{i+1}^- = F_i / V_{i+1}`, so

```
a_i^+ V_i = F_i = a_{i+1}^- V_{i+1}.
```

This is the discrete statement that the flux leaving cell `i` through the face
at `r_{i+1/2}` is the flux entering cell `i+1` through the same face. With it,
Abel summation in `i` (using `F_{-1} = 0` and `u = v = 0` on the eliminated
outer row) and periodic summation by parts in `j` give, for
`A = -(L5_r^h + d_zz^h)`,

```
<A u, v>_w
  = dz * sum_j sum_{i=0}^{n_r-2} F_i (u_{i+1,j} - u_ij)(v_{i+1,j} - v_ij)
  + (1/dz) * sum_i V_i sum_j (u_{i,j+1} - u_ij)(v_{i,j+1} - v_ij).
```

Both terms are manifestly symmetric and non-negative on the diagonal. The
energy vanishes only for a field constant in both indices whose outer boundary
value is zero, i.e. only for zero. `A` is therefore symmetric positive definite
in `<.,.>_w`, which is what licenses conjugate gradients. The operator is
**not** symmetric in the unweighted coordinate basis; the `r^3` weight is
essential, and it is the same measure that E-26 warns must not be confused with
an unweighted condition number.

Because the implementation applies the operator through a single stored array
of face coefficients `F`, the identity above holds in floating point as an
equality of the *same* number, not merely to rounding. The test verifies this
numerically anyway (`<A u, v>_w` vs `<u, A v>_w` on random vectors), and a
mutation that scales the lower coefficients by 1.05 produces a relative
asymmetry of `1.5e-02` against a gate of `1e-11`.

## Measured behaviour

At `R = 1.7`, `L_z = 2*pi`, resolutions `(17,32)`, `(33,64)`, `(65,128)`:

| quantity | coarse | medium | fine | observed order |
|---|---|---|---|---|
| manufactured max error, C | 9.319e-03 | 2.330e-03 | 5.826e-04 | 2.000, 2.000 |
| manufactured max error, A | 5.804e-03 | 1.460e-03 | 3.656e-04 | 1.991, 1.997 |
| max abs difference C vs A | 3.515e-03 | 8.704e-04 | 2.171e-04 | 2.014, 2.004 |
| that difference divided by `dz^2` | 0.0912 | 0.0903 | 0.0901 | — |
| CG iterations at `tol = 1e-12` | 34 | 73 | 146 | — |

The C-vs-A gap is a clean `0.09 dz^2`, exactly what the theory predicts: A
applies the exact axial symbol `k^2` while C applies
`(4/dz^2) sin^2(k dz/2) = k^2 - k^4 dz^2/12 + ...`. Two consequences are pinned
by the tests:

- The **zero axial mode** sees no symbol gap at all, so the `z`-average of the
  two solutions agrees to roundoff (measured `3e-16` to `1e-14`). This is the
  sharpest available statement of exactly what C shares with A: the radial
  mathematics, and nothing else.
- On a field that is a single grid-resolved axial mode times an `a + b r^2`
  radial profile, solver A is exact to roundoff while solver C is not. C's
  entire error there is the axial truncation, measured at `0.476 dz^2` and
  `0.473 dz^2` with observed order 2.01.

## What remains shared

Stated plainly, because independence claims are only worth what their exclusions
are worth:

- **`AxisymmetricGrid`.** All three solvers import it. It fixes the `(n_r, n_z)`
  layout, the uniform spacings, the `r[0] = 0` requirement, the
  endpoint-omitting periodic `z` convention, and `validate_field`. This is
  unavoidable — it is what makes the three solvers comparable on *the same*
  grid — and it is not an independent check of itself. A defect in the grid
  class is invisible to all three.
- **binary64 arithmetic**, not outward rounded, in all three.
- **The E-26 discretization mathematics.** C's radial stencil is an independent
  *transcription* of the same formulas solver A implements, not an independent
  *discretization*. A vs C therefore says nothing about whether E-26a is the
  right radial scheme; that is what solver B's different radial stencil is for.
- **The manufactured-solution family** used by the tests, which is written out
  in closed form in each test module.

## What this does and does not establish

Establishes, at the tested resolutions only:

- The axial Fourier path and the Thomas elimination shared by A and B are not
  the sole reason those two agree. A completely different axial discretization
  and a completely different linear solver reproduce the same solution to the
  predicted `O(dz^2)`.
- The reduced operator is symmetric positive definite in the `r^3 dr dz`
  measure, verified matrix-free rather than assumed.
- The E-26b axis coefficient `8` follows from the control-volume formula and is
  reproduced by the operator that conjugate gradients actually inverts.

Does **not** establish:

- Anything about the continuum. The reported residual is the algebraic residual
  of the discrete system; it is not an error bound.
- Anything about the artificial outer boundary at `r = R`. E-25's Dirichlet
  trace is a declared finite-domain condition, not a decay condition inherited
  from the `R^3` problem. No domain-truncation bound is provided.
- Anything about E-26a being the correct radial scheme, since C shares it.
- Anything about singularity formation.

## Notes and known gaps

- **The weighted residual is not a maximum-norm guarantee.** The axis cell has
  weight `V_0 = dr^4/64`, so `||b - A x||_w` sees it almost not at all. At
  `tol = 1e-12` the unweighted relative max-norm residual is measured at
  `1.3e-11` to `6.6e-11`, which is fine; but for a deliberately loosened
  `tol = 1e-2` it is **7.09**, i.e. larger than the right-hand side itself.
  `solve_realspace_poisson` therefore reports
  `algebraic_residual_max_abs_relative` in its metadata alongside the weighted
  figure, and the tests check both. Do not quote the weighted residual alone.
- **Global error gates are blind to the axis coefficient.** A mutation that
  scales the axis face coefficient by 1.5 (turning the row into `12/dr^2`)
  leaves the manufactured error at `9.278e-03` instead of `9.319e-03` and keeps
  the C-vs-A difference inside its band. Only the direct basis-vector probe of
  the axis row catches it. That is why
  `tests/test_realspace_poisson.py::test_axis_row_reproduces_the_e26b_coefficient`
  probes the matrix-free apply directly and must not be replaced by a global
  accuracy check.
- **Pairwise import isolation.** `tests/test_poisson_cross_validation.py`
  enforces that no other test module imports both A and B.
  `tests/test_realspace_poisson.py` imports A and C only, so that guard stays
  green untouched. There is currently **no** guard against a future module
  importing both B and C; if a B-vs-C comparison is ever added it should be
  given the same treatment.
- Solver C is not wired into any production trajectory path. It is an audit
  path for snapshots.
