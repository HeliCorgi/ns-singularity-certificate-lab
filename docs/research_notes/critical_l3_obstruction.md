# Critical \(L^3\) Obstruction and Revised Singularity Search Target

## Status

This note does **not** solve the three-dimensional Navier–Stokes Millennium Problem.

It records a rigorous obstruction that materially changes the current search strategy:

> A globally defined, one-scale rescaled candidate with uniformly bounded rescaled \(L^3\) norm cannot produce a finite-time singularity in \(\mathbb R^3\).

The result follows from the scale invariance of the \(L^3\) norm and the Escauriaza–Seregin–Šverák endpoint regularity theorem.

---

## 1. Physical equation

Consider the unforced incompressible Navier–Stokes equations on \(\mathbb R^3\):

\[
\partial_t u + (u\cdot\nabla)u + \nabla p = \nu\Delta u,
\qquad
\nabla\cdot u=0.
\]

Assume that \(u\) is smooth for \(0\le t<T\).

---

## 2. Isotropic dynamic rescaling

Let

\[
u(x,t)
=
\frac{1}{L(t)}
U\!\left(
\frac{x-x_*(t)}{L(t)},s(t)
\right),
\qquad
\frac{ds}{dt}=\frac{1}{L(t)^2},
\]

where \(L(t)>0\).

Set

\[
y=\frac{x-x_*(t)}{L(t)}.
\]

For every fixed \(t<T\),

\[
\begin{aligned}
\|u(t)\|_{L^3(\mathbb R^3)}^3
&=
\int_{\mathbb R^3}
L(t)^{-3}
\left|
U\!\left(
\frac{x-x_*(t)}{L(t)},s(t)
\right)
\right|^3
\,dx\\
&=
\int_{\mathbb R^3}|U(y,s(t))|^3\,dy\\
&=
\|U(s(t))\|_{L^3(\mathbb R^3)}^3.
\end{aligned}
\]

Therefore,

\[
\boxed{
\|u(t)\|_{L^3(\mathbb R^3)}
=
\|U(s(t))\|_{L^3(\mathbb R^3)}.
}
\]

This identity is exact. It does not depend on numerical resolution or an asymptotic approximation.

---

## 3. No-go theorem for bounded rescaled orbits

### Theorem

Suppose the physical solution has the above global rescaled representation near \(T\), and

\[
\sup_{s\ge s_0}
\|U(s)\|_{L^3(\mathbb R^3)}
<\infty.
\]

Then

\[
u\in L^\infty((t_0,T);L^3(\mathbb R^3)).
\]

By the endpoint \(L^\infty_tL^3_x\) regularity theorem of Escauriaza, Seregin and Šverák, \(u\) is regular at \(T\). Hence \(T\) is not a singular time.

### Consequence

None of the following can by itself establish a Navier–Stokes blow-up in \(\mathbb R^3\):

- a stationary rescaled profile \(U\) with finite \(L^3\) norm;
- a periodic rescaled orbit with uniformly bounded \(L^3\) norm;
- a quasiperiodic rescaled orbit with uniformly bounded \(L^3\) norm;
- a bounded, localized, one-scale Type-I core whose complete physical velocity is represented by the rescaling above.

A singularity candidate must escape this obstruction.

---

## 4. Anisotropic scaling

Consider a localized anisotropic representation

\[
u(x,t)
=
A(t)
U\!\left(
\frac{x_1-x_{*,1}(t)}{L_r(t)},
\frac{x_2-x_{*,2}(t)}{L_r(t)},
\frac{x_3-x_{*,3}(t)}{L_z(t)},
s(t)
\right).
\]

Then

\[
\boxed{
\|u(t)\|_{L^3}^3
=
A(t)^3 L_r(t)^2L_z(t)
\|U(s(t))\|_{L^3}^3.
}
\]

Thus, if \(U\) remains uniformly bounded in \(L^3\), a necessary condition for a singularity is

\[
\boxed{
A(t)^3L_r(t)^2L_z(t)\longrightarrow\infty
}
\]

along some sequence approaching the alleged singular time.

For standard isotropic parabolic scaling,

\[
A=L^{-1},
\qquad
L_r=L_z=L,
\]

and therefore

\[
A^3L_r^2L_z=1.
\]

So a uniformly \(L^3\)-bounded standard one-scale profile is excluded.

---

## 5. What can escape the obstruction?

At least one of the following must occur:

1. **Type-II amplitude growth**

   \[
   A^3L_r^2L_z\to\infty.
   \]

2. **Unbounded rescaled critical norm**

   \[
   \|U(s)\|_{L^3}\to\infty.
   \]

3. **A multiscale cascade**

   No single bounded profile captures all scales; critical \(L^3\) mass accumulates over an increasing number of shells.

4. **A nonlocal outer contribution**

   The local core is bounded in critical mass, but the outer field causes the global \(L^3\) norm to diverge.

5. **Failure of the proposed global representation**

   A locally fitted rescaling does not represent the full physical solution.

These are not optional refinements. One is necessary for a whole-space finite-time singularity.

---

## 6. Immediate numerical diagnostics

The repository should add the following quantities to every candidate run.

### Global critical norm

\[
M_3(t)=\int_{\mathbb R^3}|u(x,t)|^3\,dx.
\]

### Rescaled critical norm

\[
\widetilde M_3(s)
=
\int |U(y,s)|^3\,dy.
\]

### Critical scaling product

\[
Q(t)
=
A(t)^3L_r(t)^2L_z(t).
\]

### Shell decomposition

For a core scale \(L(t)\), compute

\[
M_{3,j}(t)
=
\int_{2^jL(t)\le |x-x_*(t)|<2^{j+1}L(t)}
|u(x,t)|^3\,dx.
\]

This distinguishes:

- one bounded core;
- a growing outer tail;
- a multiscale cascade;
- critical mass spread over increasingly many shells.

### Expanding rescaled domains

Compute

\[
\int_{|y|<R}|U(y,s)|^3\,dy
\]

for several increasing \(R\). A fixed small rescaled box is insufficient, because the missing critical mass may move outward in \(y\).

---

## 7. Required acceptance rule

A candidate may be promoted from “Hou-like concentration” to a “whole-space singularity candidate” only if it exhibits reproducible growth in a critical quantity.

At minimum, one of these must hold under space, time, domain and solver refinement:

\[
\|u(t)\|_{L^3}\to\infty,
\]

or

\[
A^3L_r^2L_z\to\infty,
\]

or a quantitatively resolved multiscale shell accumulation implying divergence of the global \(L^3\) norm.

Growth only in

- \(\|\omega\|_\infty\),
- \(\|u_1\|_\infty\),
- a fitted inverse time scale,
- or a small local residual

is not sufficient.

---

## 8. Implication for the current Hou program

The current early-time Hou reproduction remains scientifically useful. It validates a concentration mechanism and the numerical infrastructure.

However, a standard one-scale interpretation with velocity amplitude proportional to \(L^{-1}\) and a uniformly bounded rescaled profile is not an admissible final blow-up mechanism in \(\mathbb R^3\).

The next whole-space experiments must determine whether the Hou-like mechanism develops:

- Type-II critical growth;
- anisotropic critical growth;
- a logarithmically growing rescaled \(L^3\) norm;
- a multiscale cascade;
- or an outer-tail contribution.

If none occurs and the global \(L^3\) norm stays uniformly bounded, the candidate should be rejected as a possible whole-space singularity mechanism.

---

## 9. Lean 4 targets

Suggested formalization identifiers:

- `F-4`: isotropic \(L^3\)-scaling identity;
- `F-5`: anisotropic \(L^3\)-scaling identity;
- `F-6`: bounded scaling product plus bounded profile implies bounded physical \(L^3\);
- `F-7`: conditional bridge from the endpoint regularity theorem to exclusion of the rescaled candidate class.

`F-4` through `F-6` are algebraic/change-of-variable results.

`F-7` requires a faithful formalization of the endpoint Navier–Stokes regularity theorem or an explicitly audited theorem interface. It must not be inserted as an unproved project-specific axiom in the final proof path.

---

## 10. References

- L. Escauriaza, G. A. Seregin, V. Šverák, *\(L_{3,\infty}\)-solutions of the Navier–Stokes equations and backward uniqueness*, Russian Mathematical Surveys 58 (2003), 211–250.
- Charles L. Fefferman, *Existence and Smoothness of the Navier–Stokes Equation*, official Clay Millennium Problem description.
