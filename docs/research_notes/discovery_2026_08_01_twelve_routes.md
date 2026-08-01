# Discovery portfolio 2026-08-01 — twelve routes beyond the audited dead ends

**Status:** discovery record; every surviving route is `AUDIT REQUIRED`
**Claim boundary:** no route below is a proof, a singularity, or a numerical
singularity candidate.  `REJECTED` routes are retained because their failure
equations constrain the surviving designs.

## 0. Starting point and notation

The common source is the three-dimensional incompressible system, never the
formal five-dimensional scalar operator used elsewhere in the repository:

\[
 \partial_tu+\mathbb P(u\cdot\nabla u)=\nu\Delta u+\mathbb Pf,
 \qquad \nabla\cdot u=0,
 \qquad -\Delta p=\partial_i\partial_j(u_i u_j).
 \tag{0.1}
\]

Here \(\tau=T-t\), \(s=-\log\tau\), and all norm scalings are physical
three-dimensional scalings.  Clay labels are used as follows: (A)/(B) are
unforced regularity on \(\mathbb R^3/\mathbb T^3\), while (C)/(D) are forced
breakdown on \(\mathbb R^3/\mathbb T^3\).  A construction with \(f=0\) would be
stronger than the corresponding forced statement.

For a packet of amplitude \(A=\tau^{-\alpha}\), three spatial widths
\(L_i=\tau^{\beta_i}\), \(B=\beta_1+\beta_2+\beta_3\), and derivative scale
\(D=\tau^{-m}\), the dimensional audit used repeatedly below is

\[
\begin{array}{c|c}
\text{quantity}&\tau\text{-power}\\ \hline
\|u\|_2^2&B-2\alpha\\
\|\nabla u\|_2^2&B-2\alpha-2m\\
\|u\|_3&B/3-\alpha\\
\|\omega\|_\infty&-(\alpha+m)\\
\|(u\cdot\nabla)u\|_2,\ \|\nabla p\|_2&B/2-2\alpha-m\\
\|p\|_2&B/2-2\alpha
\end{array}
\tag{0.2}
\]

when one localized component dominates.  Packet counts multiply squared
norms and cubed \(L^3\) norms.  Every table below separately states when that
rule is being used.

---

## 1. Zeno critical-packet relay

**Label:** `FORMAL ANSATZ / SYMBOLIC CANDIDATE`
**Nonstandard feature:** a relative traveling front leaves a scale-critical
wake; there is no bounded rescaled fixed profile.

### A. Clay target

Clay (D), \(\mathbb T^3\), fixed \(\nu>0\), mean-zero smooth divergence-free
initial data (finite Fourier support is allowed), and a smooth force supported
only on a fixed finite low-mode set, with every time derivative bounded and a
smooth extension past \(T\).  The desired high modes must be produced by (0.1),
not placed in the force.

### B. Central equations

Let \(\lambda_j=2^j\), choose a divergence-free
\(W\in C_c^\infty(B(0,\rho))\), and put

\[
 x_j=x_*+C\lambda_j^{-1}e_1,\quad C>3\rho,\qquad
 w_j(x)=\lambda_jR_jW(\lambda_jR_j^T(x-x_j)).
 \tag{1.1}
\]

The supports are disjoint and

\[
 \|w_j\|_2^2=\lambda_j^{-1}\|W\|_2^2,\quad
 \|w_j\|_3^3=\|W\|_3^3,\quad
 \|\nabla w_j\|_2^2=\lambda_j\|\nabla W\|_2^2.
 \tag{1.2}
\]

The compact-support scaffold gives exact physical-space additivity, but it is
not exactly annular in Fourier space: compact support and exact band limitation
cannot hold simultaneously.  Thus \(\lambda_j^{-1}\) in (1.2) is total packet
energy.  Identifying it with an exact Littlewood--Paley shell energy requires
either a band-limited Schwartz template plus quantitative core-tail bounds, or
explicit LP leakage/cross-term estimates for the compact template.  That bridge
is an open obligation, not part of the finite-sum identity.

Let \(a_c\) be the one-cell rescaled time and set the resolved activation
times

\[
 t_j=T-{4a_c\over3}\lambda_j^{-2},\qquad
 \Delta t_j=a_c\lambda_j^{-2},
 \tag{1.3}
\]

and require a genuine NS cascade cell to leave a persistent parent wake while
creating the child:

\[
 \Phi_{a_c}(V+B_{\rm low})
 =V_{\rm wake}+2R_*V(2R_*^T(\,\cdot-y_*))+r_*,
 \qquad \sum_j\|r_{*,j}\|_{L^3}<\infty.
\tag{1.4}
\]

These are threshold times, not Heaviside creation events.  A classical PDE
orbit must contain a convergent family of smooth, super-algebraically small
future seeds, or one global relative-front field; adding a nonzero packet at
\(t_j\) would be time-discontinuous.

Equivalently, if \(J=\log_2N\), the rescaled shell amplitudes obey the
relative-periodic condition

\[
 b_{j+1}(s+2\log2)=b_j(s),\qquad b_j(s)=q(j-J(s)).
 \tag{1.5}
\]

The finite torus floor changes the shell audit.  With
\(E_j=A(\lambda_j/N)^\beta\), the previously excluded endpoint

\[
 \boxed{\beta=-1,\qquad \sigma=\gamma,\qquad
        E_j=\lambda_j^{-1},\qquad N=\tau^{-\gamma}}
 \tag{1.6}
\]

has bounded energy but one constant critical contribution per occupied scale.
Bernstein supplies only the scale-local absolute upper bound

\[
 |\Pi_N|\le C_{\rm NL}N^{5/2}E_N^{3/2}
 =C_{\rm NL}c_E^{3/2}N,
 \qquad E_N=c_EN^{-1}.
 \tag{1.7}
\]

Positive sign is a separate hypothesis.  Assume the true projected cell,
after child viscosity and all wake/off-chain loss, has a signed lower margin
\(q_*N>0\).  Dividing it by the child energy \(E_{2N}=c_E/(2N)\) gives the
continuous front envelope

\[
 \dot J=\kappa N^2,\qquad \kappa=2q_*/c_E,
 \boxed{\dot N=kN^3},\qquad k=(\log2)\kappa,
 \qquad N(t)=[N_0^{-2}-2k(t-t_0)]^{-1/2}.
 \tag{1.8}
\]

The discrete staircase \(N=2^J\) is not differentiable; (1.8) interpolates its
activation times.  It matches (1.3) when \(a_c=3/(8k)\).  Its unproved input
is the positive signed lower margin, not the elementary integration.

### C. Scaling

The ODE fixes \(\gamma=1/2\).  The table concerns the whole packet wake.

| quantity | scaling |
|---|---:|
| energy | \(O(1)\) |
| enstrophy / dissipation rate | \(\tau^{-1/2}\) |
| global \(L^3\) | \((\log(1/\tau))^{1/3}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining total dissipation | \(O(\tau^{1/2})\) |
| nonlinear and pressure-gradient \(L^2\) | \(O(\tau^{-3/4})\) |
| pressure \(L^2\) | \(O(\tau^{-1/4})\) |
| physical stage time | \(N^{-2}\asymp\tau\) |
| Fourier bandwidth | \(N\asymp\tau^{-1/2}\) |

### D. Feedback loop

\[
 \Pi_j-D_{j+1}-L_{\rm wake}-L_{\rm off}>0
 \Longrightarrow \dot N>0
 \Longrightarrow \text{child anisotropic moment }M_j
 \Longrightarrow -\nabla^2p_{j\to j+1}>0
 \Longrightarrow \text{child strain/phase lock}
 \Longrightarrow \Pi_{j+1}>0.
 \tag{1.9}
\]

For an explicit divergence-free Schwartz transmitter
\(W_a=(a\times x)e^{-|x|^2/2}\), \(a=e_1\),
\(M/\operatorname{tr}M=\operatorname{diag}(0,1/2,1/2)\).  At separation
\(d e_1\), the leading tensor is

\[
 {4\pi d^5\nabla^2p\over\operatorname{tr}M}
 =\operatorname{diag}(-12,6,6).
\tag{1.10}
\]

Thus \(-\nabla^2p\) has normalized favourable strength 12 along \(e_1\).
On the torus the exact kernel is
\(G_{\mathbb T^3}=(4\pi|x|)^{-1}+H_{\mathbb T^3}\) locally.  Taking
\(d=C\lambda_j^{-1}\) makes \(M_jd^{-5}\asymp\lambda_j^4\), but the
support/separation ratio is fixed: increasing \(C\) reduces the multipole
remainder while weakening the relay as \(C^{-5}\).  The exact-kernel tradeoff
is untested.  Pressure is a diagnostic part of the same Leray interaction as
the triad flux, not a second energy source.

### E. Obstruction audit

Energy is summable, total dissipation is finite, and the constructed lower
bound—not merely a Besov upper bound—would make \(L^3\) diverge.  The bandwidth
is not fixed; the packets are general 3D, not pure swirl.  The wake makes (1.5)
unlike a bounded Leray fixed profile, but a relative-front ancient-limit audit
is still required before claiming that every DSS Liouville hypothesis is
avoided.  At \(\gamma=1/2\), \(\sqrt\tau\|u\|_\infty\asymp1\); weak-\(L^3\)
and CKN local-energy endpoint criteria also remain to be checked.  Periodicity is
the Clay domain rather than a wall approximation.  Every finite Galerkin
truncation stops the front and is therefore only a screening model.  Smooth
forcing requires the high-frequency residual of (1.4) to vanish exactly outside
the fixed forced set and all time derivatives to extend through \(T\); it may
not be hidden in \(f\).

### F. One-hour falsification test

Use four dyadic annuli with 30--300 true helical modes per shell.  Optimize
parent phases, rotations, centers, and fixed low modes.  Promote only if the
normalized post-viscous flux is positive on at least three successive shell
transfers and off-chain leakage decreases with shell count.  Kill if two
successive optimized transfers have nonpositive margin or the phase cone is
linearly unstable.  Binary64 is enough for discovery; promotion requires
interval complex arithmetic for every triad and leakage sum.

### G. Proof chain

Explicit triad graph -> one-cell interval theorem -> invariant phase cone ->
summable leakage -> scale iteration -> smooth low force and datum -> classical
solution for every \(t<T\) -> local packet \(L^3\) lower bounds despite heat
tails -> continuation failure -> Clay (D).

---

## 2. Logarithmic endpoint coherent front

**Label:** `FORMAL ANSATZ`
**Nonstandard feature:** one moving front, no persistent packet wake; a
logarithm is the smallest correction that beats viscosity at the shell
endpoint.

### A. Clay target

Clay (D) on \(\mathbb T^3\), fixed \(\nu>0\), smooth finite-Fourier datum; force
is either zero or the fixed low-mode controller of Route 11.

### B. Central equations

For a front packet with energy \(e(N)\), define its exact dimensionless
coherence coefficient \(\chi\) by

\[
 \dot{\log N}=\chi N^{5/2}e^{1/2}-\nu\delta N^2,
 \qquad \delta\ge0.
 \tag{2.1}
\]

The endpoint choice

\[
 \boxed{e(N)=N^{-1}(\log N)^{2a},\qquad a>0}
 \tag{2.2}
\]

This logarithm belongs only to the active/front-dominated packet.  Multiplying
every member of the Route-1 \(\beta=-1\) wake by the same
\((\log N)^{2a}\) would make total energy diverge and is explicitly excluded.
An energy-safe log-wake hybrid would instead require a scale-dependent
adiabatic factor \((1+j)^a\), losing exact relative periodicity.

corresponds to a Bernstein-saturating packet

\[
 u_N(x)=N(\log N)^aW(N(x-x_N)),\qquad
 \|u_N\|_3\asymp(\log N)^a.
 \tag{2.3}
\]

Equation (2.1) becomes

\[
 \dot N=N^3[\chi(\log N)^a-\nu\delta].
 \tag{2.4}
\]

If \(\chi\ge\chi_*>0\) and \(\delta\le C\), then

\[
 \tau\sim[2\chi_*N^2(\log N)^a]^{-1},\qquad
 N\sim\tau^{-1/2}(\log(1/\tau))^{-a/2}.
 \tag{2.5}
\]

### C. Scaling

Write \(\mathcal L=\log(1/\tau)\).

| quantity | scaling |
|---|---:|
| front energy | \(\tau^{1/2}\mathcal L^{5a/2}\) |
| enstrophy / dissipation rate | \(\tau^{-1/2}\mathcal L^{3a/2}\) |
| global \(L^3\) | \(\mathcal L^a\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining dissipation | \(O(\tau^{1/2}\mathcal L^{3a/2})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-3/2}\mathcal L^{a/2}\) |
| pressure pointwise | \(\tau^{-1}\mathcal L^{2a}\) |
| viscous term pointwise | \(\tau^{-3/2}\mathcal L^{-a/2}\) |
| physical time / bandwidth | \(\tau\), \(\tau^{-1/2}\mathcal L^{-a/2}\) |

### D. Feedback loop

\[
 \chi>0\to \log^aN\text{ amplitude gain}
 \to {\text{nonlinear}\over\text{viscous}}\sim\log^aN
 \to \dot N\sim N^3\log^aN
 \to N\uparrow\to\log^aN\uparrow.
 \tag{2.6}
\]

### E. Obstruction audit

Energy vanishes and total dissipation is finite, while \(L^3\) and BKM both
diverge.  The front is infinite-band and not self-similar or globally DSS.
The construction uses true 3D Fourier interactions.  Its exact collision with
known theory is concentrated in one missing assertion: persistent
\(\chi_*>0\); energy cancellation alone cannot imply it.

### F. One-hour falsification test

Search 10,000 low-shell phase seeds, evolve the best 24^3 states for one
turnover, and repeat at 32^3.  Measure \(H_0,\ldots,H_3,\chi,\delta\) every
step.  Success means \(\chi>0.05\), creation of the next shell, front-identity
residual below 5%, and energy-cancellation defect below \(10^{-12}\) at both
resolutions.  Kill if \(\chi\) loses sign before one transfer.  Floating point
screens phases; exact/interval sums are needed only after promotion.

### G. Proof chain

Coherent packet family -> positive phase-cone flux lower bound -> comparison
for (2.4) -> finite-time front -> critical-norm divergence -> compactness of
infinite-band approximants -> smooth forcing/datum -> noncontinuation -> Clay
(D).

---

## 3. Branching Floquet tree

**Label:** `FORMAL ANSATZ`
**Nonstandard feature:** rescaled evolution returns to several smaller copies,
not to itself.

### A. Clay target

Initially Clay (D) on \(\mathbb T^3\) with a smooth low-mode force; a
free-space Schwartz/unforced matching would strengthen it toward (C).

### B. Central equations

Use \(u=\tau^{-\alpha}U((x-X)/\tau^\beta,s)\),
\(\alpha+\beta=1\), and let the resolved lobe count be
\(K(\tau)=\tau^{-\kappa}\).  Bounded total energy uniquely requires

\[
 K\tau^{3\beta-2\alpha}=O(1)
 \quad\Longrightarrow\quad
 \boxed{\kappa=5\beta-2},\qquad {2\over5}<\beta<{1\over2}.
 \tag{3.1}
\]

After one rescaled period \(S\), a parent must split according to

\[
 U_{\rm parent}(y,S)=e^{(1-\beta)S}
 \sum_{r=1}^{b}R_rU_0\!\left(R_r^T{y-d_r\over\rho}\right)+r_S,
 \quad \rho=e^{-\beta S},\quad b=e^{\kappa S},
 \tag{3.2}
\]

with \(b\) an integer and a contraction on \(r_S\).  The child centers obey

\[
 a_j'=\beta a_j+V_j(U),\qquad X_j=\tau^\beta a_j(s).
 \tag{3.3}
\]

A concrete dyadic choice is

\[
 \beta={9\over20},\quad\alpha={11\over20},\quad
 \kappa={1\over4},\quad b=2,\quad S=4\log2.
 \tag{3.4}
\]

### C. Scaling

The table includes all \(K\) children.

| quantity | scaling |
|---|---:|
| total energy | \(O(1)\) |
| total enstrophy / dissipation rate | \(\tau^{-9/10}\) |
| global \(L^3\) | \(\tau^{-11/60}\) |
| per-core vorticity maximum | \(\tau^{-1}\) |
| remaining total dissipation | \(O(\tau^{1/10})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-31/20}\) |
| pressure pointwise | \(\tau^{-11/10}\) |
| physical time | \(\tau\) |
| Fourier bandwidth / lobe count | \(\tau^{-9/20}\), \(\tau^{-1/4}\) |

### D. Feedback loop

\[
 \text{parent strain}\to\text{two-lobe instability}
 \to\text{cross-pressure separation}
 \to\text{equal-energy children}
 \to K\uparrow
 \to\text{more pressure interaction}\to\text{next split}.
 \tag{3.5}
\]

### E. Obstruction audit

Energy and total dissipation pass; \(L^3\) diverges; bandwidth and lobe count
both diverge.  Equation (3.2) is outside fixed-profile and global periodic-DSS
no-go hypotheses.  Occupied volume is
\(K\tau^{3\beta}=\tau^{11/10}\to0\), so packing is not the immediate
obstruction.  The integer \(K\) must be defined through smooth lobe resolution,
not inserted as a discontinuous PDE state.

### F. One-hour falsification test

Solve one parent-to-two-child space-time BVP in a divergence-free
Gaussian--Hermite basis with 100--300 unknowns.  Success requires two separated
lobes, energy split \(1/2+1/2\), decreasing PDE and phase residuals under basis
refinement.  Kill if pressure recombines the lobes or the residual converges
only by moving energy into the truncation edge.  Float Newton--Krylov suffices;
the final BVP needs radii polynomials.

### G. Proof chain

One branching-cell theorem -> contraction of the branching renormalization
operator -> locally finite infinite tree -> pressure/tail matching -> smooth
force and datum -> energy/dissipation -> \(L^3\) divergence -> finite physical
time -> Clay connection.

---

## 4. Local \(3/5\)--\(2/5\) pressure-Hessian Floquet lens

**Label:** `FORMAL ANSATZ / PROOF CANDIDATE`; the globally decaying periodic
profile version is `REJECTED`.
**Nonstandard feature:** periodic only on a growing rescaled core; a moving
outer energy flux prevents reduction to a global DSS profile.

### A. Clay target

Clay (C), \(\mathbb R^3\), preferably \(f=0\), fixed \(\nu>0\), with a
divergence-free Schwartz datum.  It is fully Cartesian and non-axisymmetric.

### B. Central equations

Set

\[
 A=\nabla u,\quad S={A+A^T\over2},\quad
 \xi={\omega\over|\omega|},\quad
 \alpha=\xi^TS\xi,\quad \eta=(I-\xi\xi^T)S\xi,
 \quad P=\nabla^2p.
 \tag{4.1}
\]

Differentiating (0.1) and differentiating \(\alpha\) along a vortex line gives
the exact identity

\[
\boxed{
 D_t\alpha=|\eta|^2-\alpha^2-\xi^TP\xi
\nu\xi^T\Delta S\xi
{2\nu\over|\omega|}
 ((I-\xi\xi^T)\Delta\omega)\cdot S\xi .}
\tag{4.2}
\]

For \(Q=\operatorname{tr}(A^2)=|S|^2-|\omega|^2/2\), the directional pressure
Hessian is

\[
 \xi^TP(x)\xi=-{Q(x)\over3}
 +\operatorname{p.v.}\int_{\mathbb R^3}
 {3(\widehat z\cdot\xi)^2-1\over4\pi|z|^3}Q(x-z)\,dz.
 \tag{4.3}
\]

An equatorial \(Q>0\) belt can therefore contribute negatively to
\(P_{\xi\xi}\).  The quantitative lens condition is

\[
 P_{\xi\xi}\le-(2+\delta)\alpha^2,\qquad
 R_\nu\ge-\frac\delta2\alpha^2,
 \tag{4.4}
\]

where \(R_\nu\) denotes the last two terms of (4.2).  It implies
\(D_t\alpha\ge(1+\delta/2)\alpha^2\).

For

\[
 u=\tau^{-\alpha_*}U((x-X)/\tau^{\beta_*},s),\qquad
 \alpha_*+\beta_*=1,
 \tag{4.5}
\]

the rescaled PDE is

\[
 U_s+\alpha_*U+\beta_*y\cdot\nabla U-c\cdot\nabla U
 +\mathbb P(U\cdot\nabla U)
 =\nu e^{-(1-2\beta_*)s}\Delta U.
 \tag{4.6}
\]

A nonzero asymptotically periodic finite-energy core has zero mean dilation
only when \(\alpha_*=3\beta_*/2\).  Thus

\[
 \boxed{\alpha_*={3\over5},\qquad\beta_*={2\over5}}.
 \tag{4.7}
\]

It cannot be a globally localized periodic \(U\).  The surviving boundary
condition is periodicity only in \(|y|<R(s)\), \(R(s)\to\infty\), with the
moving-boundary flux

\[
 \mathcal F_R=\int_{\partial B_R}
 \left[\left(U+{2\over5}y\right){|U|^2\over2}+PU
 -\nu e^{-s/5}\nabla{|U|^2\over2}-R'n{|U|^2\over2}\right]\cdot n.
 \tag{4.8}
\]

The mean of (4.8) must balance core dissipation over each Floquet cycle.

### C. Scaling

| quantity | scaling |
|---|---:|
| core energy | \(O(1)\) |
| enstrophy / dissipation rate | \(\tau^{-4/5}\) |
| global \(L^3\) | \(\tau^{-1/5}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining total dissipation | \(O(\tau^{1/5})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-8/5}\) |
| pressure pointwise / Hessian | \(\tau^{-6/5}\), \(\tau^{-2}\) |
| physical/rescaled time | \(\tau\), \(s=-\log\tau\) |
| Fourier bandwidth | \(\tau^{-2/5}\) |

### D. Feedback loop

\[
 Q_{\rm equator}>0\overset{(4.3)}\longrightarrow
 P_{\xi\xi}<-2\alpha^2\overset{(4.2)}\longrightarrow
 \dot\alpha>\alpha^2\longrightarrow |\omega|\uparrow
 \longrightarrow Q\text{ satellites reform}
 \longrightarrow Q_{\rm equator}>0.
 \tag{4.9}
\]

The outer flux (4.8) must close the last arrow without making the profile
globally periodic.

### E. Obstruction audit

The scale is Type II, has bounded energy, finite dissipation, divergent \(L^3\),
and unbounded bandwidth.  It is not pure swirl.  A globally decaying periodic
solution of the inviscid limiting rescaled equation falls under known
discrete-self-similar/profile nonexistence hypotheses, so that version is
explicitly rejected; local Floquet plus nonperiodic outer matching is the only
surviving version.  A periodic-box pilot would be invalid unless image effects
are removed with free-space pressure.

### F. One-hour falsification test

Optimize a divergence-free Gaussian--Hermite core plus 4--6 equatorial
satellites using analytic/free-space pressure at two basis orders.  Require
the margin in (4.4) to exceed 0.1 on an open core for one turnover and require
the measured outer flux not to vanish as the matching radius grows.  Kill if
\(-P_{\xi\xi}/\alpha^2\le1\), the vorticity direction leaves the favourable
eigenspace immediately, or the optimized tail converges to a global periodic
profile.  Float optimization is enough; the final open-set sign requires
interval quadrature.

### G. Proof chain

Explicit lens family -> free-space pressure-Hessian interval bound -> local
Floquet orbit -> moving outer matching -> viscous Fredholm correction ->
trapping tube from Schwartz data -> Riccati comparison -> energy/tail control
-> noncontinuation -> Clay (C).

---

## 5. Carrier/envelope helical cloud

**Label:** `FORMAL ANSATZ`
**Nonstandard feature:** the occupied spatial diameter and the Fourier
wavelength shrink with different exponents.

### A. Clay target

Clay (D) on \(\mathbb T^3\), fixed \(\nu>0\), finite-Fourier smooth datum and
fixed low-mode smooth forcing.  A free-space localization could instead aim at
(C).

### B. Central equations

Let the cloud envelope and carrier scale be

\[
 L_e=\tau^{2/5},\qquad N=\tau^{-1/2},\qquad A=\tau^{-1/2},
 \qquad K=(L_eN)^3=\tau^{-3/10}.
 \tag{5.1}
\]

Place \(K\) carrier cells of width \(N^{-1}\) in the envelope and write them
with exact curls so that divergence is not merely asymptotic:

\[
 u_c=\sum_{q=1}^{K}\nabla\times
 \left[N^{-1}A\,\chi(N(x-x_q))
 \Re\{a_qe^{iNk_q\cdot(x-x_q)}\}\right],
 \qquad k_q\cdot a_q=0.
 \tag{5.2}
\]

Each cell uses both helicities.  The leading carrier balance and child
creation condition are

\[
 \partial_\theta W_q+\mathbb P_\zeta(W_q\cdot\nabla_\zeta W_q)
 =\nu\Delta_\zeta W_q+\Pi_{q-1\to q}-\Pi_{q\to q+1},
 \tag{5.3}
\]

where \(\theta=N^2(t-t_q)\).  A single isolated periodic carrier cell is
impossible: taking the inner product in (5.3) with no incoming flux gives
\(\nu\|\nabla_\zeta W_q\|_2^2=0\).  Hence the flux-through term is a necessary,
falsifiable part of the ansatz.

Neighbouring cell moments satisfy

\[
 M_q\asymp A^2N^{-3}=\tau^{1/2},\qquad
 d_q\asymp N^{-1},\qquad
 \nabla^2p_{q\to q+1}\asymp M_qd_q^{-5}=\tau^{-2},
 \tag{5.4}
\]

which matches the time derivative of cell strain \(AN\asymp\tau^{-1}\).

### C. Scaling

The table includes \(K\) nearly orthogonal or disjoint cells.

| quantity | scaling |
|---|---:|
| total energy | \(\tau^{1/5}\) |
| enstrophy / dissipation rate | \(\tau^{-4/5}\) |
| global \(L^3\) | \(\tau^{-1/10}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining total dissipation | \(O(\tau^{1/5})\) |
| nonlinear and pressure-gradient \(L^2\) | \(\tau^{-9/10}\) |
| pressure \(L^2\) | \(\tau^{-2/5}\) |
| physical cell time | \(N^{-2}=\tau\) |
| bandwidth / envelope diameter | \(\tau^{-1/2}\), \(\tau^{2/5}\) |

### D. Feedback loop

\[
 \text{heterochiral flux}\to\text{new carrier cells}
 \to K,N\uparrow\to\text{neighbour pressure Hessian }\tau^{-2}
 \to\text{phase/strain lock}\to\text{heterochiral flux}.
 \tag{5.5}
\]

### E. Obstruction audit

The spatial core is not tied to \(N^{-1}\), so one-scale profile theorems do
not apply directly.  Energy vanishes, dissipation is finite, \(L^3\) diverges,
and the active band is unbounded.  It is neither pure swirl nor a fixed
Galerkin system.  The exact no-incoming-flux identity above rejects the most
tempting closed-cell simplification.  Front-resolution tests must grow both
the number of carrier cells and modes per cell.

### F. One-hour falsification test

Use 2--4 shells and \(K=4,8,16\) curl-built cells.  Measure incoming/outgoing
shell flux, neighbour pressure-Hessian sign, and off-cloud leakage at 64^3 and
96^3.  Success requires a positive scale-normalized flux and a pressure sign
stable under doubling \(K\).  Kill if removing the imposed incoming flux does
not collapse the cell (implementation defect), or if the true outgoing flux is
always below viscous loss.  Float is adequate except for a promoted triad sign.

### G. Proof chain

Flux-through cell BVP -> helical phase cone -> pressure-neighbour lower bound
-> envelope packing and modulation -> summable leakage -> smooth forcing/datum
-> energy/dissipation and critical lower bound -> finite-time front -> Clay.

---

## 6. Opposite-helicity capacitor

**Label:** `FORMAL ANSATZ`

### A. Clay target

Clay (C), \(\mathbb R^3\), \(f=0\), fixed \(\nu>0\), divergence-free Schwartz
data.  The field is general Cartesian 3D.

### B. Central equations

With helical projectors

\[
 \mathcal P^\pm={1\over2}(I\pm|D|^{-1}\operatorname{curl}),
 \qquad u=v^++v^-,
 \tag{6.1}
\]

assume both fields lie in \([N,(1+\varepsilon)N]\).  Then

\[
 \operatorname{curl}v^\pm=\pm Nv^\pm+r^\pm,\qquad
 \|r^\pm\|_2\le\varepsilon N\|v^\pm\|_2,
 \tag{6.2}
\]

and the Lamb vector has the exact leading cancellation identity

\[
 \boxed{u\times\omega=-2N v^+\times v^-+u\times(r^++r^-).}
 \tag{6.3}
\]

Thus each chirality is almost Beltrami while their cross-interaction is large.
The local helicity obeys

\[
 \partial_t(u\cdot\omega)+\nabla\cdot
 \left[(u\cdot\omega)u+(p-|u|^2/2)\omega\right]
 =\nu\Delta(u\cdot\omega)
 -2\nu\partial_ju_i\partial_j\omega_i.
 \tag{6.4}
\]

Choose \(\|v^+\|_2=\|v^-\|_2\), so signed helicity nearly cancels while
\(|u\times\omega|\) does not.  Use the \(3/5\)--\(2/5\) local scale of Route 4.

### C. Scaling

| quantity | scaling |
|---|---:|
| energy | \(O(1)\) |
| enstrophy / dissipation rate | \(\tau^{-4/5}\) |
| global \(L^3\) | \(\tau^{-1/5}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining dissipation | \(O(\tau^{1/5})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-8/5}\) |
| pressure pointwise | \(\tau^{-6/5}\) |
| physical time / bandwidth | \(\tau\), \(\tau^{-2/5}\) |
| signed packet helicities | \(H^\pm\sim\pm\tau^{-2/5}\), net \(o(H^+)\) |

### D. Feedback loop

\[
 H^+\simeq-H^-\to |u\times\omega|\simeq2N|v^+\times v^-|
 \to\text{strain}\to N,|H^\pm|\uparrow
 \to\text{larger cross-helicity interaction}.
 \tag{6.5}
\]

### E. Obstruction audit

The full field is deliberately not Beltrami, so Beltrami depletion does not
apply.  Energy/dissipation/ESS scalings pass and bandwidth diverges.  Pure
swirl and fixed-band exclusions are irrelevant.  The direct obstruction is
off-annulus convolution: if it destroys the opposite-helicity overlap in less
than one turnover, the capacitor is dead.

### F. One-hour falsification test

Construct localized helical packets at 64^3, 96^3 and 128^3, optimize relative
phase/translation, and monitor net helicity, Lamb-vector norm, forward flux and
off-band energy for one turnover.  Promote if chirality purity exceeds 0.9,
net helicity is below 5% of either signed part, normalized forward flux exceeds
0.2, and leakage stays below 50%.  Kill on rapid homochiralization or inverse
transfer.  Float screens; interval complex triads certify a promoted sign.

### G. Proof chain

Localized helical basis -> bound (6.2) -> invariant opposite-helicity phase
cone -> rescaled orbit -> leakage/tail control -> Schwartz connection ->
critical and BKM divergence -> noncontinuation -> Clay (C).

---

## 7. Coherent vortex-line fold: velocity cancellation, strain addition

**Label:** `FORMAL ANSATZ`

### A. Clay target

Clay (C), \(\mathbb R^3\), unforced fixed-viscosity flow from smooth,
divergence-free, rapidly decaying closed vortex tubes.

### B. Central equations

The exact Biot--Savart strain is

\[
 S(x)={3\over8\pi}\operatorname{p.v.}\int
 {[(z\times\omega)\otimes z+z\otimes(z\times\omega)]\over|z|^5}\,dy,
 \qquad z=x-y.
 \tag{7.1}
\]

Seek a quadrupolar fold arrangement for which

\[
 \sum_j u_j(x_*)=O(\Gamma/r),\qquad
 \xi_*^T\sum_jS_j(x_*)\xi_*\ge c_S{\Gamma\over r^2},
 \tag{7.2}
\]

so the odd velocity pieces cancel but differentiated even strain pieces add.
The vortex-line identity

\[
 \nabla\cdot(|\omega|\xi)=0,
 \qquad \partial_\ell\log|\omega|=-\nabla\cdot\xi
 \tag{7.3}
\]

makes the required growing fold curvature explicit.  Choose

\[
 r=\tau^{1/2},\quad L_{\rm seg}=\tau^{7/20},\quad
 M=\tau^{-1/10},\quad \Omega=\tau^{-1},\quad
 \Gamma=\Omega r^2=\Gamma_0.
 \tag{7.4}
\]

Then the effective amplitude equation is

\[
 \dot\Omega=\left(c_S-{c_D\nu\over\Gamma_0}\right)\Omega^2
 +R_{\rm reconnect}+R_{\rm geometry}.
 \tag{7.5}
\]

### C. Scaling

| quantity | scaling |
|---|---:|
| velocity amplitude | \(\tau^{-1/2}\) |
| energy | \(\tau^{1/4}|\log\tau|\) |
| enstrophy / dissipation rate | \(\tau^{-3/4}\) |
| global \(L^3\) | \(\tau^{-1/12}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining dissipation | \(O(\tau^{1/4})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-3/2}\) |
| pressure pointwise / Hessian | \(\tau^{-1}\), \(\tau^{-2}\) |
| physical time / bandwidth / fold count | \(\tau\), \(\tau^{-1/2}\), \(\tau^{-1/10}\) |

### D. Feedback loop

\[
 M\uparrow\to\text{velocity cancellation + strain addition}
 \to\dot\Omega\ge c\Omega^2\to r\downarrow
 \to\text{more folds fit}\to M\uparrow.
 \tag{7.6}
\]

### E. Obstruction audit

Energy and total dissipation pass; \(L^3\), BKM, and bandwidth diverge.  The
price paid to evade vortex-direction coherence criteria is explicit:
\(\int\kappa\,d\ell\asymp\tau^{-1/4}\to\infty\).  The design is not pure
swirl or a one-core profile.  Viscous reconnection is the immediate physical
kill mechanism, not a detail to postpone.

### F. One-hour falsification test

Generate volume-preserving tubes with \(M=4,8,16\), evaluate free-space
Biot--Savart on 128^3, and measure normalized velocity cancellation and strain
addition in (7.2).  Kill if the strain cancels with velocity or reconnection
time is shorter than stretching time.  Float geometry search is enough;
promoted kernel signs need interval cubature.

### G. Proof chain

Explicit folded tubes -> interval Biot--Savart sign -> packing/curvature
persistence -> comparison in (7.5) -> reconnection and diffusion bounds ->
energy/local-energy control -> critical divergence -> Schwartz embedding ->
Clay (C).

---

## 8. Lagrangian viscous Cauchy defect

**Label:** singular-tube version `REJECTED`; material-cover regularity lemma
`PROOF CANDIDATE`.

### A. Clay target

The rejected branch aimed at (C) on \(\mathbb R^3\), unforced, fixed
\(\nu>0\), Schwartz data.  Its reusable branch aims at unforced regularity (A)
by showing that sufficiently coherent material tubes cannot realize BKM growth.

### B. Central equations

Let

\[
 \dot X(a,t)=u(X(a,t),t),\quad F=\nabla_aX,\quad\det F=1,
 \quad G=F^{-1}F^{-T},\quad \Omega=\omega\circ X.
 \tag{8.1}
\]

Piola's identity yields

\[
 (\Delta_x\omega)\circ X
 =\partial_{a_\alpha}(G_{\alpha\beta}\partial_{a_\beta}\Omega).
 \tag{8.2}
\]

The Cauchy-defect variable \(Z=F^{-1}\Omega\) removes stretching exactly:

\[
 \boxed{\partial_tZ=\nu F^{-1}\nabla_a\cdot
       [G\nabla_a(FZ)].}
 \tag{8.3}
\]

Try to hide diffusion along the most stretched material direction with

\[
 s(F)=(\tau^{-q},\tau^{q/2},\tau^{q/2}),\qquad
 G\sim\operatorname{diag}(\tau^{2q},\tau^{-q},\tau^{-q}).
 \tag{8.4}
\]

Volume preservation forces transverse radius \(r\sim\tau^{q/2}\), so the
transverse diffusion time is \(r^2\sim\tau^q\), while
\(|\Omega|\sim\tau^{-q}\).

### C. Scaling and rejection

For one stretched tube, ignoring logarithmic slender-body corrections:

| quantity | scaling |
|---|---:|
| material volume | \(O(1)\) |
| velocity / energy | \(\tau^{-q/2}\), \(\tau^{-q}\) |
| enstrophy / dissipation rate | \(\tau^{-2q}\) |
| global \(L^3\) | \(\tau^{-q/2}\) |
| vorticity maximum | \(\tau^{-q}\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-3q/2}\) |
| pressure pointwise | \(\tau^{-q}\) |
| physical time / transverse bandwidth | \(\tau\), \(\tau^{-q/2}\) |

The required inequalities are mutually incompatible:

\[
 \underbrace{q\ge1}_{\rm BKM},\qquad
 \underbrace{q<1}_{\rm diffusion\ slower\ than\ remaining\ time},\qquad
 \underbrace{q<1/2}_{\rm finite\ dissipation},\qquad
 \underbrace{q\le0}_{\rm bounded\ energy}.
 \tag{8.5}
\]

### D. Mechanism and failure

\[
 F\text{ stretches}\to F\omega_0\text{ grows}
 \to\text{cross-section shrinks}
 \to G_\perp\text{ grows}
 \to\text{transverse diffusion accelerates}.
 \tag{8.6}
\]

The last arrow reverses the desired singular feedback.  For regularity, it is
exactly the useful arrow.

### E. Obstruction audit

Equation (8.5), not a numerical failure, rejects the fixed material-tube
singularity.  Selecting a time-dependent nonmaterial high-vorticity subset
reintroduces boundary flux and loses the closure (8.3).  The possible reuse is
a covering theorem: every BKM-scale high-vorticity set would have to contain
enough coherent material tubes for transverse Poincare dissipation to violate
the energy inequality.  Such a cover is not currently proved.

### F. One-hour falsification test

No PDE run is needed.  Verify the rational feasibility problem (8.5), then
test the proposed regularity reuse on saved smooth trajectories by measuring
material singular values and transverse gradients.  Kill the simple cover
lemma if high vorticity consistently jumps between material tubes faster than
one stretching time.  Exponent feasibility is exact rational arithmetic;
trajectory screening may be float.

### G. Proof chain for the reused route

Material-tube cover -> transverse Poincare inequality -> lower bound from
\(G_\perp\) -> BKM growth implies nonintegrable dissipation -> energy
contradiction -> control outside the cover by a critical norm -> Sobolev
continuation -> Clay (A).

---

## 9. Adaptive analyticity impedance

**Label:** `PROOF CANDIDATE / SYMBOLIC CANDIDATE`

### A. Clay target

Clay (B), unforced \(\mathbb T^3\), arbitrary smooth mean-zero
divergence-free datum and fixed \(\nu>0\).  The opposite branch supplies a
quantitative necessary mechanism for a (D) front.

### B. Central equations

Define the analytic energy and its first two weighted moments by

\[
 E_\rho={1\over2}\sum_ke^{2\rho|k|}|\widehat u_k|^2,\quad
 M_\rho=\sum_k|k|e^{2\rho|k|}|\widehat u_k|^2,\quad
 D_\rho=\sum_k|k|^2e^{2\rho|k|}|\widehat u_k|^2.
 \tag{9.1}
\]

If \(B_k\) is the exact Leray nonlinearity, let

\[
 \mathcal T_\rho=\sum_ke^{2\rho|k|}
 \Re(\overline{\widehat u_k}\cdot B_k(u,u)).
 \tag{9.2}
\]

Then, with a time-dependent radius,

\[
 \boxed{\dot E_\rho=\dot\rho M_\rho-\nu D_\rho+\mathcal T_\rho.}
 \tag{9.3}
\]

Choose the radius dynamically:

\[
 \boxed{\dot\rho=-{(\mathcal T_\rho)_+\over M_\rho}
                  +{\nu\over2}{D_\rho\over M_\rho}.}
 \tag{9.4}
\]

It follows exactly that

\[
 \dot E_\rho\le-\frac\nu2D_\rho.
 \tag{9.5}
\]

The new phase-sensitive impedance is

\[
 \boxed{\mathfrak I_\rho={ (\mathcal T_\rho)_+\over\rho M_\rho}}.
 \tag{9.6}
\]

Since \(\dot\rho\ge-\mathfrak I_\rho\rho\), finite
\(\int_0^T\mathfrak I_\rho dt\) keeps \(\rho(T^-)>0\), hence gives analytic
continuation.  Every singular branch must therefore have
\(\int_0^T\mathfrak I_\rho dt=\infty\).

### C. Scaling

A candidate saturating obstruction uses front exponents
\((\gamma,\sigma)=(9/20,1/4)\), so \(A=\tau^{-11/20}\).

| quantity | obstruction-front scaling |
|---|---:|
| front energy | \(\tau^{1/4}\) |
| enstrophy / dissipation rate | \(\tau^{-13/20}\) |
| global \(L^3\) | \(\tau^{-1/10}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining dissipation | \(O(\tau^{7/20})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-31/20}\) |
| pressure pointwise | \(\tau^{-11/10}\) |
| physical time / bandwidth | \(\tau\), \(\tau^{-9/20}\) |
| impedance | \(\mathfrak I_\rho\asymp\tau^{-1}\) |

### D. Regularity/singularity dichotomy

\[
 \int\mathfrak I<\infty\to\rho(T)>0\to\text{regular continuation},
 \tag{9.7a}
\]

or

\[
 \text{coherent triad phase}\to\mathcal T_\rho>0
 \to\mathfrak I\sim\tau^{-1}\to\rho\downarrow0
 \to\text{front amplitude grows}\to\mathcal T_\rho\uparrow.
 \tag{9.7b}
\]

### E. Obstruction audit

Equations (9.3)--(9.6) use exact NS coefficients and retain phase information
lost by ordinary norm estimates.  The construction does not claim that
analyticity-radius methods themselves are new; the candidate object is the
adaptive positive-flux impedance.  A proof of (B) still requires a universal
integrable bound on it.  A finite Fourier run can only test the identity and
phase sensitivity, not establish regularity or blow-up.

### F. One-hour falsification test

On a 64^3--96^3 pseudospectral trajectory compute (9.1)--(9.6) for
\(\rho N\in[0.25,2]\), and compare the original phases with shellwise scrambled
phases.  Kill the proposed discriminant if \(\mathfrak I\) is insensitive to
phase or if the radius ODE reaches zero despite an integrable measured
impedance.  Finite sums verify (9.3) to roundoff; a promoted algebraic identity
can be checked exactly.

### G. Proof chain

Finite-sum identity -> limit and radius-ODE existence -> positivity of radius
under \(L^1_t\) impedance -> analytic-to-Sobolev continuation -> universal
phase-defect estimate making \(\mathfrak I\in L^1_t\) -> Clay (B).

---

## 10. Exact Fourier-moment front and viscosity-barrier hierarchy

**Label:** `SYMBOLIC CANDIDATE / PROOF CANDIDATE`

### A. Clay target

The barrier branch targets unforced regularity (B); failure of the barrier in
a controlled phase cone supplies a possible forced breakdown route (D).
Initial data are smooth, mean-zero and divergence-free on \(\mathbb T^3\),
with fixed \(\nu>0\).

### B. Central equations

For \(r\ge0\), define

\[
 H_r=\sum_k|k|^{2r}|\widehat u_k|^2,\qquad
 {1\over2}\dot H_r=T_r-\nu H_{r+1}+F_r,\qquad
 N_r^2={H_{r+1}\over H_r}.
 \tag{10.1}
\]

Direct differentiation gives the exact moving-front identity

\[
\boxed{
 {d\over dt}\log N_r=
 {T_{r+1}+F_{r+1}\over H_{r+1}}
 -{T_r+F_r\over H_r}
 -\nu\left({H_{r+2}\over H_{r+1}}-{H_{r+1}\over H_r}\right).}
\tag{10.2}
\]

Log-convexity makes the last bracket

\[
 D_r=N_r^2\delta_r\ge0,\qquad
 \delta_r={H_{r+2}H_r\over H_{r+1}^2}-1.
 \tag{10.3}
\]

If the nonlinear difference \(G_r\) satisfies, for \(r=0,1,2\),

\[
 G_r\le(1-\eta)\nu D_r+g_r,\qquad g_r\in L^1_t,\quad\eta>0,
 \tag{10.4}
\]

then

\[
 \mathcal M_r=\log N_r+\eta\nu\int_0^tD_r-\int_0^tg_r
 \tag{10.5}
\]

is nonincreasing.  Bounded \(N_0,N_1,N_2\) gives
\(H_3=H_0N_0^2N_1^2N_2^2\), hence continuation.

Conversely, for a localized coherent front with energy \(e=N^{-q}\),
\(0\le q<1\), (10.2) motivates

\[
 \dot N\sim\chi N^{(7-q)/2},
 \qquad \boxed{N\sim\tau^{-2/(5-q)}}.
 \tag{10.6}
\]

### C. Scaling of the violating branch

Let \(\gamma=2/(5-q)\).

| quantity | scaling |
|---|---:|
| front energy | \(\tau^{q\gamma}\) |
| enstrophy / dissipation rate | \(\tau^{-\gamma(2-q)}\) |
| global \(L^3\) | \(\tau^{-\gamma(1-q)/2}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining dissipation | \(O(\tau^{(1+q)/(5-q)})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-2(4-q)/(5-q)}\) |
| pressure pointwise | \(\tau^{-2(3-q)/(5-q)}\) |
| physical time / bandwidth | \(\tau\), \(\tau^{-2/(5-q)}\) |

### D. Dichotomy

\[
 \text{regularity via (10.4)}
 \quad\text{or}\quad
 \bigl[G_r\text{ beats viscosity}\bigr]
 \ \vee\ \bigl[\delta_r\to0\text{ and the spectrum becomes thin}\bigr].
 \tag{10.7}
\]

The first violating branch is Routes 1--2; the thin-spectrum branch becomes a
separate rigidity target.

### E. Obstruction audit

The identity is not a self-consistency model and is independent of
axisymmetry.  The power-law violating branch has bounded energy, finite total
dissipation, divergent \(L^3\), and moving bandwidth for \(q<1\).  But (10.4)
is not known for arbitrary NS solutions, while \(\chi>0\) in (10.6) is not
implied by the exact identity.  Galerkin data only diagnose which branch is
active before cutoff.

### F. One-hour falsification test

Measure \(G_r/(\nu D_r)\), \(\delta_r\), and the exact residual of (10.2) on
existing low-Re periodic runs and a short phase-optimized run.  Kill the simple
barrier (10.4) if the ratio is unbounded even in the benchmark; then test
whether adding a measured phase-defect term yields an integrable \(g_r\).
Float is sufficient for rejection; log-convexity and finite Fourier identities
are exact-algebra targets.

### G. Proof chain

Moment identity -> log-convexity -> universal phase-defect inequality (10.4)
-> three bounded moment ratios -> \(H^3\) bound -> strong continuation -> Clay
(B); or invariant positive-flux cone -> comparison (10.6) -> Clay (D).

---

## 11. Smooth low-frequency phase control

**Label:** `SYMBOLIC CANDIDATE / AUDIT REQUIRED`

### A. Clay target

Clay (D), \(\mathbb T^3\), fixed \(\nu>0\), smooth finite-Fourier initial datum,
and a \(C^\infty_{x,t}\) divergence-free force supported on one fixed finite
low-mode set \(L\).

### B. Central equations

Split \(u=b^*+v\), prescribe low coefficients
\(b_k^*=\rho_k(t)e^{i\phi_k(t)}\), and put \(P_Hf=0\).  The only allowed force is

\[
 \boxed{f_k=\dot b_k^*+\nu|k|^2b_k^*
 +iP_k\sum_{\ell+m=k}(m\cdot\widehat u_\ell)\widehat u_m,
 \qquad k\in L.}
 \tag{11.1}
\]

For this to extend smoothly through \(T\), the high-high backreaction must
satisfy

\[
 \boxed{G_k(v,v)=P_k\sum_{\substack{\ell+m=k\\\ell,m\notin L}}
 (m\cdot\widehat v_\ell)\widehat v_m
 =g_k(t)\in C^\infty([0,T+\epsilon]),\qquad k\in L.}
 \tag{11.2}
\]

The strongest and easiest-to-certify condition is \(G_k\equiv0\).  In a
helical basis, the controller holds triad phases

\[
 \Theta_{k\ell m}=\phi_\ell+\phi_m-\phi_k+arg C_{k\ell m}
 \tag{11.3}
\]

inside a forward-flux cone near \(\pi/2\), while all actual front transfer is
high-high.

### C. Scaling

Use the high field of Route 2; the finite low reservoir only adds \(O(1)\)
energy.

| quantity | high-field scaling |
|---|---:|
| total energy | \(O(1)\) |
| high enstrophy / dissipation | \(\tau^{-1/2}\mathcal L^{3a/2}\) |
| global \(L^3\) | \(\mathcal L^a\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining dissipation | finite |
| nonlinear / pressure gradient | \(\tau^{-3/2}\mathcal L^{a/2}\) |
| pressure | \(\tau^{-1}\mathcal L^{2a}\) |
| physical time / bandwidth | \(\tau\), \(\tau^{-1/2}\mathcal L^{-a/2}\) |
| directly forced high modes | exactly zero |

### D. Feedback loop

\[
 f_L\to\Theta\approx\pi/2\to\chi>0\to N\uparrow
 \to G_L(v,v)=g_L\text{ remains smooth}
 \to f_L\text{ remains smooth}\to\Theta\approx\pi/2.
 \tag{11.4}
\]

### E. Obstruction audit

Smooth-force high-frequency decay is met by exact zero direct forcing, not by
an asymptotic excuse.  A bounded low-mode strain can move \(N\) only
exponentially, so low-high transfer alone cannot blow up; high-high local
transfer is mandatory.  The front evades fixed-band and Galerkin no-go, and
the high \(L^3\) norm evades ESS.  Condition (11.2) is the sharp unresolved
collision: arbitrary inverse-designed residual forcing would not be Clay
smooth.

### F. One-hour falsification test

On a finite exact helical triad graph solve
\(G_L(v,v)=0\) while maximizing the forward coefficient \(\chi\), using random
phase search followed by polynomial least squares or a Groebner calculation.
Kill if \(G_L=0\) algebraically forces \(\chi\le0\).  A floating witness is
discovery evidence only; a promoted witness needs exact algebraic numbers or
interval complex coefficients.

### G. Proof chain

Finite phase-algebra witness -> invariant phase cone -> exact/smooth low force
-> high-front PDE theorem -> force extension past \(T\) -> energy/dissipation
-> critical divergence -> strong noncontinuation -> Clay (D).

---

## 12. Nested outer/inner strain relay

**Label:** disjoint-core two-way version `REJECTED`; overlapping Reynolds-stress
variant is absorbed into Route 5.

### A. Clay target

The rejected ansatz aimed at (C), \(\mathbb R^3\), unforced, fixed \(\nu>0\),
from Schwartz data.

### B. Central equations

Take an outer \((\alpha_o,\beta_o)=(3/5,2/5)\) core and an inner

\[
 \alpha_i={11\over20},\qquad\beta_i={9\over20},\qquad
 L_i\ll L_o.
 \tag{12.1}
\]

In inner coordinates the outer gradient is leading order:

\[
 V_s+\alpha_iV+\beta_i y\cdot\nabla V
 +(\Sigma y)\cdot\nabla V+\Sigma V+\mathbb P(V\cdot\nabla V)
 =\nu e^{-s/10}\Delta V+O(e^{-s/20}),
 \tag{12.2}
\]

where \(\nabla u_o(X,t)=\tau^{-1}\Sigma(s)\).  Thus outer-to-inner coupling
works.  But the reverse pressure Hessian at outer distance is

\[
 {\nabla^2p_i(L_o)\over\nabla^2p_o}
 \asymp {E_i/L_o^5\over\tau^{-2}}
 =\tau^{5(\beta_i-\beta_o)}=\tau^{1/4}\to0.
 \tag{12.3}
\]

Leading feedback would require \(K\asymp\tau^{-1/4}\) inner cores, whereas
disjoint packing supplies at most

\[
 K_{\max}\asymp(L_o/L_i)^3=\tau^{-3/20}.
 \tag{12.4}
\]

The required effective packing dimension is five, larger than physical
dimension three.

### C. Scaling of one inner core

| quantity | scaling |
|---|---:|
| energy | \(\tau^{1/4}\) |
| enstrophy / dissipation rate | \(\tau^{-13/20}\) |
| global \(L^3\) | \(\tau^{-1/10}\) |
| vorticity maximum | \(\tau^{-1}\) |
| remaining dissipation | \(O(\tau^{7/20})\) |
| nonlinear and pressure-gradient pointwise | \(\tau^{-31/20}\) |
| pressure pointwise | \(\tau^{-11/10}\) |
| physical time / inner bandwidth | \(\tau\), \(\tau^{-9/20}\) |

### D. Failed feedback loop

\[
 \text{outer strain}\to\text{inner growth}
 \to \underbrace{\text{inner pressure return}}_{\tau^{1/4}\text{ too small}}
 \not\to\text{outer strain}.
 \tag{12.5}
\]

### E. Obstruction audit

The one-way multiscale coupling passes energy, dissipation, ESS and bandwidth
checks, but it is not a closed singular mechanism.  Multiplying disjoint inner
cores enough to repair it violates packing before any known theorem is needed.
Overlapping oscillatory carriers can have a Reynolds stress without disjoint
packing and are therefore retained only as Route 5.

### F. One-hour falsification test

Compute the free-space pressure Hessian of two-scale Gaussian packets for
\(L_i/L_o=2^{-m}\).  Confirmation of the \((L_i/L_o)^5\) law kills the disjoint
two-way route; failure signals a multipole or overlap term that must be isolated.
Float is sufficient for exponent rejection, while the asymptotic kernel can be
proved symbolically.

### G. Reuse chain

Discard disjoint cores -> construct overlapping carrier family -> derive
nonzero averaged Reynolds stress -> close envelope pressure equation -> Route
5.  No direct Clay chain remains for the rejected version.

---

## 13. Coverage map and genuinely new equations

The requested discovery categories are covered as follows.

| requested category | routes |
|---|---|
| whole-space anisotropic/multiscale self-similarity | 4, 7, 12 |
| multi-center concentration | 1, 3, 5, 7 |
| different spatial and Fourier scales | 5, 12 |
| pressure as positive nonlocal feedback | 1, 4, 5, 7 |
| Lagrangian deformation-gradient degeneration | 8 |
| stretching/diffusion phase synchronization | 4, 5, 7, 8 |
| helicity/local helicity subsystem | 6, 11 |
| finite Fourier datum to infinite bandwidth | 1, 2, 10, 11 |
| smooth low-frequency phase control | 11 |
| fixed point replaced by periodic/relative/branching orbit | 1, 3, 4 |
| new critical quantity | 9, 10 |
| monotone quantity or exclusive regularity dichotomy | 8, 9, 10 |

Routes 4--10 are independent of the Hou/axisymmetric reduction and start from
general Cartesian, Lagrangian, or periodic Fourier NS.  The principal new
derived objects are:

1. the finite-floor \(\beta=-1\) critical staircase (1.6);
2. the cubic Zeno bandwidth law (1.8);
3. the branching renormalization boundary condition (3.2);
4. the directional pressure-lens Riccati identity (4.2)--(4.4);
5. the exact viscous Cauchy-defect PDE (8.3);
6. the adaptive analyticity energy/radius pair (9.3)--(9.6);
7. the exact Fourier-moment front identity (10.2);
8. the low-mode forcing cancellation condition (11.2).

Thus the portfolio contains more than five worked derivations and more than
three departures from a single self-similar core.

## 14. Comparative score

Scores are 1 (weak) to 5 (strong).  They rank research priority, not truth.
For “Clay distance” and “connection,” a higher score means fewer currently
visible logical bridges.  `Int` means interval-arithmetic suitability and
`Lean` means suitability after the analytic content has been proved.

| # | route | Clay | novel | closes | no-go | compute | Int | Lean | falsify | connects | total |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Zeno packet relay | 4 | 5 | 3 | 5 | 4 | 4 | 3 | 5 | 4 | **37** |
| 2 | logarithmic endpoint front | 4 | 4 | 3 | 5 | 5 | 4 | 3 | 5 | 3 | **36** |
| 11 | smooth low-frequency phase control | 4 | 4 | 3 | 5 | 5 | 4 | 4 | 5 | 2 | **36** |
| 9 | adaptive analyticity impedance | 3 | 4 | 3 | 4 | 5 | 4 | 4 | 5 | 3 | **35** |
| 10 | Fourier-moment barrier/front | 3 | 4 | 3 | 4 | 5 | 4 | 5 | 5 | 2 | **35** |
| 5 | carrier/envelope cloud | 3 | 5 | 2 | 5 | 4 | 3 | 2 | 5 | 3 | **32** |
| 3 | branching Floquet tree | 3 | 5 | 2 | 5 | 3 | 2 | 2 | 4 | 3 | **29** |
| 4 | local pressure Floquet lens | 3 | 4 | 2 | 4 | 3 | 3 | 2 | 5 | 3 | **29** |
| 6 | opposite-helicity capacitor | 3 | 4 | 2 | 4 | 4 | 3 | 2 | 5 | 2 | **29** |
| 7 | vortex-line fold | 3 | 5 | 2 | 4 | 3 | 3 | 2 | 5 | 2 | **29** |
| 8 | Lagrangian Cauchy defect | 2 | 4 | 3 | 2 | 4 | 4 | 4 | 5 | 2 | **30** |
| 12 | nested disjoint relay | 1 | 4 | 1 | 2 | 5 | 4 | 3 | 5 | 1 | **26** |

The top three research experiments are therefore:

1. **Route 1:** it is the only route born from a concrete correction to the
   repository's obstruction geometry and it supplies an actual critical
   \(L^3\) lower-bound architecture.  Decisive hole: a true NS cascade cell.
2. **Route 2:** it is the smallest active-front modification that makes the
   nonlinear/viscous ratio diverge while retaining finite total dissipation.
   Decisive hole: persistent positive \(\chi\).
3. **Route 11:** it uses exactly the freedom granted by Clay (D) without
   directly forcing high modes.  Decisive hole: simultaneous forward flux and
   smooth/zero high-high backreaction at the controlled low modes.

Routes 9 and 10 are the leading regularity-side fallbacks: failure to construct
Routes 1/2 should be translated into an impedance or viscosity-barrier bound,
not discarded as an unstructured negative search.

## 15. Pilot executed for the top route

The reproducible CPU-small pilot is

```text
python -m experiments.run_zeno_packet_relay_pilot \
  --config configs/zeno_packet_relay_pilot_v2.json \
  --output-dir outputs/zeno_packet_relay_pilot_v2
```

with seed `20260801`.  It records config, runtime provenance, SHA-256 sidecars,
and a manifest.  Its checked scope is deliberately narrow:

- levels \(J=8,16,24,32\): energy remained below 2, enstrophy divided by
  bandwidth approached 2, and both critical packet mass and the finite-floor
  modeled Besov sum grew as \(J+1\), as two separate kinematic checks rather
  than simultaneous exact properties of one template;
- the explicit rank-two Gaussian vortex gave normalized leading pressure
  stretch 12, and 20,000 seeded physical orientations found
  \(11.9976793898\);
- the leading tensor's \(d^{-5}\) invariant had zero relative spread, a
  tautological regression of the homogeneous formula rather than an exact
  packet-pressure/remainder test;
- the integrated front ODE had formal blow-up time 1 and
  \(N\sqrt{T-t}\) relative spread `2.22e-16`;
- the dedicated unit/bundle-integrity suite passed 22 tests.

The preserved v1 run used an unrealizable rank-one energy moment and is
rejected as a physical packet witness; v2 does not overwrite it.  V2 supports
the finite-sum algebra, a realizable leading-tensor sign, and ODE implementation
only.
It does **not** test the positive NS flux lower bound, invariant phase cone, or
existence of an orbit.  Those are the shortest kill conditions, not conclusions
smuggled into the pilot.

## 16. Handoff to the verification phase

The next phase should proceed in this order:

1. implement one true helical shell-to-shell transfer using the exact Leray
   coefficients and measure post-viscous margin plus all off-chain leakage;
2. search simultaneously for the Route-11 constraint \(G_L(v,v)=0\);
3. repeat at two larger mode sets and stop before the front reaches the top
   third of the retained spectrum;
4. if the normalized margin decays or changes sign, mark Routes 1/2/11
   `REJECTED` and feed the measurements into Routes 9/10;
5. if it persists, enclose the finite triad graph with interval complex
   arithmetic and prove an invariant phase cone;
6. only then solve the one-cell relative/branching BVP and bound leakage;
7. prove local \(L^3\) lower bounds with parabolic tails, not exact support
   disjointness;
8. prove that the low-mode force extends smoothly beyond \(T\);
9. construct the infinite-band PDE limit and check local energy/pressure;
10. audit the exact bridge to the selected Clay statement.

The complete unverified top-route specification is saved separately in
`docs/candidates/CANDIDATE_SOLUTION_ZENO_PACKET_RELAY.md`.  The finite-floor
correction and the retained historical error are recorded in
`docs/research_notes/track_f_shell_constraints_finite_floor_erratum.md` and
`docs/research_notes/track_f_shell_constraints.md` respectively.
