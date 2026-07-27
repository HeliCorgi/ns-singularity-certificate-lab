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

## Independence

The solver does not call the existing cylindrical diagnostic operators. It
uses Fourier diagonalization and a local complex tridiagonal solver.

A separate physical-space implementation applies the discrete PDE rows and the
outer identity row. This supplies a second residual path. The continuum error
is measured independently against an analytic polynomial-times-cosine
manufactured solution with nonzero outer boundary data.

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
- The finite outer radius is an artificial boundary. No full-space Green tail
  or domain-truncation bound is provided.
- Passing this gate permits later trajectory-code design review only. It is not
  evidence for a singularity.
