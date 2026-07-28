# Whole-space Poisson and low-frequency wall progress

## Scope

This milestone does not solve the Navier–Stokes Millennium Problem.  It closes
one numerical and logical gap between the periodic finite-cylinder Hou setup
and a possible `R^3` computation.

## Exact result added

For one axial Fourier mode, define

```text
A_k psi = -[d_rr + (3/r)d_r - k^2] psi.
```

Let `psi_inf` solve `A_k psi=f`, regular at `r=0` and decaying at infinity.
Let `psi_R` solve the same equation on `0<=r<=R` with `psi_R(R)=0`.
Then

```text
psi_R(r)-psi_inf(r)
 = -psi_inf(R) [I1(kr)/r]/[I1(kR)/R],   k>0,

psi_R(r)-psi_inf(r)
 = -psi_inf(R),                          k=0.
```

For a source supported inside the wall,

```text
psi_inf(R) = [1/(2R^2)] integral rho^3 f(rho) d rho,  k=0.
```

Thus the zero mode has an algebraic `R^-2` tail.  Fixed nonzero modes have a
wall coefficient containing `K1(kR)/I1(kR)`, which is exponentially small for
large `kR`.

## New solver

`src/ns_certificate_lab/free_space_poisson.py` implements

```text
-L5 psi = source,
L5 = d_rr + (3/r)d_r + d_zz,
```

using the exact free-radial Fourier Green kernel

```text
G_0(r,rho) = 1/[2 max(r,rho)^2],
G_k(r,rho) = I1(k r_<) K1(k r_>)/(r_< r_>), k>0.
```

The radial wall is removed analytically.  The `z` direction remains periodic
after zero padding, so padding convergence is a period-image diagnostic rather
than a rigorous free-space-z certificate.

## Manufactured gate

Five-dimensional Gaussian test, interior relative `L2` error:

| z padding | padded period | relative error |
|---:|---:|---:|
| 1 | 4 | 1.18936e-2 |
| 2 | 8 | 4.64317e-3 |
| 4 | 16 | 4.19483e-3 |

The first doubling reduces image error by a factor of 2.56.  The later plateau
is dominated by radial quadrature and source truncation; it is not a rigorous
error bound.

## Low-frequency wall gate

A localized radial source was tested in the first periodic Fourier mode while
increasing the axial period.

| axial period | `kR` | relative finite-wall error in `r<=2` |
|---:|---:|---:|
| 1 | 25.1327 | 9.0581e-4 |
| 2 | 12.5664 | 2.7636e-4 |
| 4 | 6.2832 | 2.2521e-4 |
| 8 | 3.1416 | 6.3546e-3 |
| 16 | 1.5708 | 3.2936e-2 |
| zero mode | 0 | 9.1256e-2 |

At large `kR`, the true wall correction falls below finite-difference error, so
the modal-shape prediction is not resolved.  At periods 8 and 16, the exact
homogeneous correction predicts the finite/free difference with relative
defects `7.95e-3` and `1.12e-3`.  For the zero mode the predicted constant shift
has defect `3.37e-4`.

## Scientific conclusion

The earlier small wall sensitivity at short axial period is explained by the
Fourier spectral gap.  It is not evidence that the mechanism is independent of
the wall in `R^3`.

A whole-space transition must explicitly resolve:

1. the zero axial mode;
2. the continuous low-frequency limit `k -> 0`;
3. nonperiodic `z` image/tail error;
4. radial and axial domain enlargement independently;
5. the nonlinear time evolution using the free-space elliptic recovery.

## Verification performed

```text
13 focused tests passed
```

These include Gaussian manufactured solutions, Bessel cross-checks, exact zero-
mode moment formulas, finite/free wall-correction shapes, low-frequency stress
tests, and fault/input validation.

The complete legacy test suite was started and reached 112 passed tests before
the execution time limit.  It did not report a failure before timeout.  This is
not a claim that the full current `fable5-mainline` suite passes, because the
latest GitHub branch could not be fetched in this execution environment.

## Next proof-relevant step

Integrate the free-radial solver into a nonperiodic-z time-step prototype with
compactly supported smooth initial data.  Before candidate search, establish a
four-way convergence table in

```text
dt, dr, dz, Rmax, Zmax, z-padding/free-space-tail method.
```

Every run must record global and shell-decomposed `L3`, the three PDE term
blocks, divergence, energy balance, and elliptic tail error.  Only candidates
that escape the critical `L3` obstruction under all refinements should be
promoted.
