# Finite-cylinder Poisson gate implementation report

## Implemented milestone

This bundle implements the next recorded repository gate:

\[
-\left(\partial_{rr}+\frac{3}{r}\partial_r+\partial_{zz}\right)\psi_1=\omega_1
\]

on a finite cylinder with periodic `z`, regular axis behavior, and explicit
nonzero outer Dirichlet data.

## Files

- `src/ns_certificate_lab/finite_cylinder_poisson.py`
- `tests/test_finite_cylinder_poisson.py`
- `experiments/run_poisson_gate.py`
- `configs/poisson_gate.json`
- `docs/finite_cylinder_poisson.md`
- `outputs/poisson_gate_v1/`

The implementation intentionally uses a new module name and does not overwrite
an existing `poisson.py`.

## Numerical method

- periodic `z`: discrete Fourier transform;
- radial direction: newly assembled second-order tridiagonal system;
- axis limit: coefficient `8` in the negative operator row;
- outer radius: explicit Dirichlet identity row;
- solve: guarded complex Thomas algorithm;
- independent audit: separate physical-space second-order periodic stencil.

## Recorded results

Manufactured solution resolutions:

| `(nr,nz)` | weighted relative L2 error | observed order | Fourier matrix residual max | cross-stencil defect max |
|---|---:|---:|---:|---:|
| `(17,32)` | `1.9421867e-3` | — | `9.95e-13` | `5.67355e-2` |
| `(33,64)` | `5.0027631e-4` | `1.9569` | `1.03e-11` | `1.42805e-2` |
| `(65,128)` | `1.2697393e-4` | `1.9782` | `1.17e-10` | `3.57619e-3` |

The cross-stencil defect is not expected to be at roundoff because the solver
uses Fourier axial differentiation while the audit path uses an independent
second-order physical-space axial stencil. Its observed orders were about
`1.9902` and `1.9976`.

Selected raw dense condition estimates increased with refinement. These values
are representation-dependent because the PDE rows and unscaled Dirichlet
identity row have different magnitudes; they are diagnostics, not stability
bounds.

## Validation performed

Against the previously uploaded repository snapshot:

```text
53 passed
```

This consists of the snapshot's existing tests plus 11 new Poisson-gate tests.
The experiment's seven acceptance checks passed. Python compilation and all
four manifest payload hashes were also checked.

## Important integration note

The public GitHub branch is newer than the uploaded snapshot. The new files use
unique paths and only depend on the stable `AxisymmetricGrid` interface, so they
are designed to be copied into the current repository without replacing the
newer Cartesian and time-convergence work. The full current GitHub test suite
must still be run locally after copying. Do not update `STATUS.md` to claim the
gate is complete until that current-branch replay succeeds.

## Reproduction

From the repository root:

```powershell
$env:PYTHONPATH='src'
$env:PYTHONDONTWRITEBYTECODE='1'
python -m pytest -q
python -m experiments.run_poisson_gate --config configs/poisson_gate.json --output-dir outputs/poisson_gate_local
```

Choose a fresh output directory; the experiment refuses to overwrite an
existing evidence directory.

## What this does not establish

- no full-space boundary-tail estimate;
- no interval arithmetic;
- no nonlinear time evolution;
- no singularity candidate;
- no evidence that the Millennium problem is solved.
