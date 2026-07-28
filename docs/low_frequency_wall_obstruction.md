# Low-frequency obstruction to inferring whole-space behavior from a periodic cylinder

## Exact modal correction

For one axial Fourier mode, write

```text
-[d_rr + (3/r)d_r - k^2] psi = f.
```

Let `psi_inf` be regular at the axis and decaying at radial infinity.  Let
`psi_R` satisfy the same equation on `0 <= r <= R` with `psi_R(R)=0`.
The difference is homogeneous.

For `k > 0`, the regular homogeneous solution is `I1(kr)/r`, hence

```text
psi_R(r) - psi_inf(r)
  = -psi_inf(R) [I1(kr)/r] / [I1(kR)/R].
```

For `k = 0`, the regular homogeneous solution is constant:

```text
psi_R(r) - psi_inf(r) = -psi_inf(R).
```

If the source is zero beyond `R`,

```text
psi_inf(R) = M/(2R^2),
M = integral rho^3 f(rho) d rho,
```

for the zero mode.  Its wall influence is therefore algebraic.

For fixed `k > 0`, the exterior free solution is proportional to `K1(kr)/r`.
The wall-correction coefficient contains `K1(kR)/I1(kR)`, which is exponentially
small as `kR` grows.

## Consequence for the Hou cylinder

A short periodic `z` interval enforces a large minimum nonzero wavenumber.  For
period one, `k_min = 2 pi`; at radial wall `R=1`, `k_min R = 2 pi`.  Nonzero-mode
radial wall effects can therefore be small even when the corresponding
nonperiodic or low-frequency whole-space field would have a much longer tail.

Thus:

> Small wall sensitivity in a fixed short-period cylinder is not evidence of
> radial-wall independence in `R^3`.

The zero mode and the limit `k -> 0` must be tested separately.  Increasing the
axial period is a necessary stress test before moving to a nonperiodic solver.

## Numerical gate

`experiments/run_low_frequency_wall_gate.py` compares:

- the free-radial Green solver;
- the finite-cylinder zero-Dirichlet solver;
- the exact homogeneous wall-correction shape.

It uses the first periodic mode for periods `1,2,4,8,16`, followed by the zero
mode.  The observed wall error must rise as `kR` enters the low-frequency
regime, and the exact modal correction must explain the difference once the
wall effect exceeds finite-difference error.

## Proof relevance

This result does not construct or exclude a singularity.  It closes a logical
hole in the transition plan: periodic-cylinder wall experiments cannot be
promoted to whole-space evidence without controlling the continuous low-
frequency spectrum and the algebraic zero-mode tail.
