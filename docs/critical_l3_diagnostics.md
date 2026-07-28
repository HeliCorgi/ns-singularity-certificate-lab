# Critical-`L^3` diagnostics

## Purpose

A whole-space finite-time singularity must evade the endpoint critical-regularity obstruction.  In particular, a standard one-scale representation

\[
u(x,t)=L(t)^{-1}U((x-x_*(t))/L(t),s)
\]

satisfies

\[
\|u(t)\|_{L^3(\mathbb R^3)}=\|U(s)\|_{L^3(\mathbb R^3)}.
\]

A uniformly `L^3`-bounded global rescaled orbit therefore cannot be the final blow-up mechanism.  The numerical search must distinguish Type-II amplitude growth, anisotropic concentration, growth of the rescaled critical norm, a multiscale cascade, and an outer-tail contribution.

The module `ns_certificate_lab.critical_l3` adds the required finite-domain diagnostics.  It does not prove the endpoint regularity theorem and does not convert a cylinder computation into an `R^3` result.

## Measured quantities

For an axisymmetric physical velocity

\[
|u|=\sqrt{(u^r)^2+(u^\theta)^2+(u^z)^2},
\]

the represented-domain critical mass is evaluated with the Cartesian volume measure

\[
M_3=2\pi\int\!\int |u(r,z)|^3 r\,dr\,dz.
\]

The state adapter uses `u^theta = r u1`; it never substitutes `u1` for physical swirl.

The diagnostic amplitude and widths are preregistered as

- `A = max |u|`,
- `L_r`: RMS radius of the density `|u|^3` about the symmetry axis,
- `L_z`: periodic RMS axial width about the circular critical-density center,
- `L = (L_r^2 L_z)^(1/3)`,
- `Q = A^3 L_r^2 L_z`.

`Q` is an exact scaling factor only for an explicitly supplied anisotropic rescaling with stable dimensionless profile.  When estimated from RMS widths, it is a numerical shape diagnostic.

## Multiscale shells

The code partitions represented critical mass into axis-centered dyadic shells around the axial center:

\[
\rho=\sqrt{r^2+d_{\rm per}(z,z_*)^2}.
\]

It records shell masses, fractions, Shannon entropy, effective shell count, maximum shell fraction and the fraction in the outer radial cells.  Increasing shell entropy alone is not proof of a cascade; it must survive space, time, domain and solver refinement.

An off-axis axisymmetric maximum is a ring, not a point.  The current shell geometry is intended for an axis-focused Hou-like core.  A separate toroidal-shell diagnostic is required before interpreting off-axis concentration.

## Post-processing checkpoints

Example:

```powershell
python -m experiments.analyze_critical_l3 `
  outputs/hou_early_time/checkpoints/checkpoint_nr257_nz512_t000.npz `
  outputs/hou_early_time/checkpoints/checkpoint_nr257_nz512_t001.npz `
  --output-dir outputs/critical_l3_nr257
```

The output contains:

- `critical_l3_snapshots.csv`,
- `critical_l3_summary.json`,
- `manifest.json`.

## Interpretation gate

A run remains only a finite-cylinder concentration observation unless all of the following are addressed:

1. nonperiodic `z` and a free-space elliptic path;
2. radial and axial domain enlargement;
3. resolved outer-tail contribution;
4. space/time/integrator convergence of `M3`, `Q` and shell masses;
5. a critical quantity that grows rather than merely `||omega||_infinity`;
6. interval or equivalent validated bounds before any proof claim.

Uniform boundedness of a trustworthy global `L^3` estimate rejects the candidate as a whole-space finite-time singularity mechanism.
