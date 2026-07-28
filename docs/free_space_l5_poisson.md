# Free-radial Green solver for `-L5 psi1 = omega1`

## Purpose

The finite-cylinder Hou solver cannot by itself test the Clay whole-space problem.
This module removes the artificial radial wall from the elliptic recovery step.
It is an intermediate gate toward an `R^3` axisymmetric solver.

It does **not** yet remove the `z` periodicity rigorously.  The implementation
zero-pads the source in `z`; convergence under increasing padding measures the
period-image sensitivity.

## Five-dimensional interpretation

For an axisymmetric function in the first four coordinates of `R^5`,

```text
L5 = d_rr + (3/r) d_r + d_zz
```

is the ordinary five-dimensional Laplacian.  After Fourier transforming in `z`,

```text
-[d_rr + (3/r)d_r - k^2] psi_hat = omega_hat.
```

The regular and decaying homogeneous solutions for `k > 0` are

```text
I1(k r)/r,    K1(k r)/r.
```

Their weighted Wronskian is `r^3 W = -1`, giving

```text
G_k(r,rho) = I1(k r_<) K1(k r_>) / (r_< r_>).
```

For `k = 0`,

```text
G_0(r,rho) = 1 / [2 max(r,rho)^2].
```

Thus

```text
psi_hat(r,k) = integral G_k(r,rho) omega_hat(rho,k) rho^3 d rho.
```

The radial decay condition at infinity is built into `K1`; no outer radial
Dirichlet wall is imposed.

## Numerical status

The current implementation uses:

- FFT in the padded `z` box;
- explicit Green matrices in `r`;
- composite trapezoid quadrature with the `rho^3` measure;
- scaled `I1/K1` approximations to avoid overflow;
- an optional SciPy special-function backend for cross-checking.

The dependency-free Bessel approximation is not interval-certified.  It is an
independent floating-point prototype, not a proof artifact.

## Manufactured tests

1. `z`-independent radial Gaussian;
2. a single periodic Fourier mode with radial Gaussian profile;
3. a five-dimensional Gaussian on a truncated `r,z` box, with padding factors
   `1,2,4` to verify decreasing period-image error;
4. Bessel-product comparison against SciPy when available;
5. symmetry, positivity and input-fault tests.

## Remaining work before a whole-space claim

- replace periodized `z` with a rigorously controlled free-space transform or
  a certified image-tail bound;
- bound radial quadrature truncation and Bessel approximation errors;
- cross-check against direct five-dimensional Green quadrature on small grids;
- propagate the solver through the nonlinear time step;
- use compactly supported smooth `R^3` initial data;
- demonstrate convergence in time, space, radial extent and axial extent;
- add interval arithmetic and a machine-checkable certificate.
