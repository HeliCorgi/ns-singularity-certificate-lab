# Finite-cylinder `-L5` Poisson gate

## Scope

This milestone implements the floating-point boundary-value problem

\[
-\left(\partial_{rr}+\frac{3}{r}\partial_r+\partial_{zz}\right)\psi_1
=\omega_1,
\qquad 0\le r\le R,
\]

with periodic `z`, regular even behavior at the axis, and explicitly supplied
Dirichlet data at `r=R`.

It is an independent elliptic numerical gate. It is not a whole-space solver,
does not bound the effect of the artificial outer boundary, and does not use
interval arithmetic.

## Discretization

The periodic axial direction is diagonalized with NumPy's discrete Fourier
transform. For each axial angular wavenumber `q`, the radial problem is
tridiagonal.

At an interior radial point,

\[
-\mathcal L_5\psi_i =
\left(-\frac1{\Delta r^2}+\frac{3}{2r_i\Delta r}\right)\psi_{i-1}
+\left(\frac2{\Delta r^2}+q^2\right)\psi_i
+\left(-\frac1{\Delta r^2}-\frac{3}{2r_i\Delta r}\right)\psi_{i+1}.
\]

For a smooth even function at the axis,

\[
\left(\partial_{rr}+\frac3r\partial_r\right)\psi\big|_{r=0}
=4\partial_{rr}\psi(0).
\]

With an even ghost point, `psi(-dr)=psi(dr)`, this gives

\[
(-\mathcal L_5\psi)_0
\approx \frac{8}{\Delta r^2}(\psi_0-\psi_1)+q^2\psi_0.
\]

The final radial row is the identity equation `psi(R,z)=g(z)`.

### Known structural limitation

The first interior row `i=1` sits at `r_1 = dr`, so its lower coefficient is

\[
-\frac1{\Delta r^2}+\frac{3}{2 r_1\Delta r}
= -\frac1{\Delta r^2}+\frac{3}{2\Delta r^2}
= +\frac{1}{2\Delta r^2} > 0 .
\]

Two consequences follow and are recorded here rather than hidden:

- The matrix has a **positive off-diagonal entry**, so it is **not an
  M-matrix**, and there is **no discrete maximum principle** for this
  discretization. Row `i=1` is also **not diagonally dominant**: for `q=0` the
  diagonal is `2/dr^2` while the off-diagonal absolute row sum is
  `0.5/dr^2 + 2.5/dr^2 = 3/dr^2`.
- The per-mode linear solve is a plain Thomas elimination with only a pivot
  floor guard. There is **no pivoting**, so no unconditional backward-stability
  argument is available from the code.

Empirically the solver is second-order accurate and the sampled dense
condition numbers stay moderate at the tested resolutions, and the residual and
manufactured-solution gates pass. That is measured behavior at those
resolutions only. The scheme is **not** "verified stable", and this limitation
must not be described as resolved without a separate proof or a change of
discretization. `tests/test_finite_cylinder_poisson.py` pins the row-`1`
coefficient so this known limitation cannot change meaning silently.

## Independence

The solver does not call the existing cylindrical diagnostic operators. It
uses Fourier diagonalization and a local complex tridiagonal solver.

A separate physical-space implementation applies the discrete PDE rows and the
outer identity row. This supplies a second residual path. The continuum error
is measured independently against an analytic polynomial-times-cosine
manufactured solution with nonzero outer boundary data.

### Scope of the A/B independence claim

This solver (B) and `ns_certificate_lab.poisson` (A) differ **only** in their
radial stencils. Their axial treatment is identical: both diagonalize `z` with
NumPy's discrete Fourier transform, both take the axial wavenumbers from the
same NumPy sample-frequency helper, and both invert the resulting radial
systems with a locally written Thomas elimination. A defect in the wavenumber
convention, the Nyquist mode, the real/imaginary bookkeeping, or the
tridiagonal elimination would move A and B identically and would not be
detected by `tests/test_poisson_cross_validation.py`. The A/B cross-validation
must therefore be described as a *radial-stencil* cross-validation, not as a
full independence result.

`ns_certificate_lab.realspace_poisson` (solver C) is a third path added to
break that common mode: three-point periodic finite differences in `z` and a
matrix-free preconditioned conjugate-gradient solve, with no transform and no
tridiagonal elimination anywhere. It shares the `r^3`-flux radial mathematics
with solver A, so it does **not** replace the A/B radial comparison. See
`docs/realspace_poisson.md`.

## Acceptance gate

The initial gate requires:

- direct inspection tests for the axis coefficient `8`;
- direct inspection of the outer identity row;
- a global-sign fault test;
- nonzero analytic outer boundary data;
- at least three refinements;
- decreasing analytic error with observed order at least `1.85`;
- small Fourier-system and physical-space matrix residuals;
- rejected nonperiodic, nonfinite, and shape-invalid inputs;
- selected dense 2-norm condition estimates recorded as diagnostics.

## Limitations

- All arithmetic is binary64 and is not outward rounded.
- Fourier differentiation is used only in the solve; the residual audit uses a
  second-order physical-space periodic stencil.
- Dense condition numbers are sampled for selected modes and are not a
  stability theorem.
- The radial matrix is not an M-matrix and row `i=1` is not diagonally
  dominant; see "Known structural limitation" above.
- The finite outer radius is an artificial boundary. No full-space Green tail
  or domain-truncation bound is provided.
- Passing this gate permits later trajectory-code design review only. It is not
  evidence for a singularity.
