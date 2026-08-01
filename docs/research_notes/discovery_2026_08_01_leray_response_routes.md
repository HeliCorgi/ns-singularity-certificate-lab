# Discovery portfolio 2026-08-01 — exact Leray response and front actions

**Status: discovery record; all surviving mechanisms require audit**

This is the second portfolio. It does not repeat the finite-floor discovery.
It starts from the open obligation left by the first portfolio: construct or
exclude a positive signed transfer for the true three-dimensional Leray
nonlinearity. All Fourier formulas use

\[
\partial_t\widehat u_k=\mathcal N_k-\nu|k|^2\widehat u_k,
\qquad
\mathcal N_k=-iP_k\sum_{\ell+m=k}
(m\cdot\widehat u_\ell)\widehat u_m.
\tag{0.1}
\]

For physical-space scaling tables, \(\tau=T-t\). “Pressure scale” is the
scale of the pressure gradient, not a second energy source independent of the
Leray term.

## 1. Equiangular difference-cancelling Leray relay

**Label: SYMBOLIC CANDIDATE**

### A. Clay target

Clay (D) on \(\mathbb T^3\), with a fixed low-mode smooth controller; Clay (B)
if it can be initiated without forcing. The datum is an exact rational,
finite Fourier polynomial and hence smooth and divergence free. Viscosity is
fixed and positive.

### B. Central equations

Set

\[
p=(1,1,0),\quad q=(1,0,1),\quad c=p+q,\quad
n=p\times q=(1,-1,-1),
\tag{1.1}
\]

\[
u=B e_3\sin(p\cdot x)+C e_2\cos(q\cdot x)+Dn\cos(c\cdot x).
\tag{1.2}
\]

Exact Leray convolution yields

\[
P_{p-q}\mathbb P(u\cdot\nabla u)=0,\qquad
\Pi_c={BCD\over2},
\tag{1.3}
\]

\[
D_{\nu,c}=9\nu D^2,\qquad
\|N_{\rm off}\|_2^2={3D^2\over8}(B^2+C^2).
\tag{1.4}
\]

At \(B=C=1,D=1/8,\nu=1/40\), the exact post-viscous child margin is
\(151/2560>0\).

### C. Scaling

Under integer dilation \(s\) with a fixed number of modes and critical energy:

| quantity | scale |
|---|---:|
| energy | \(s^{-1}\) |
| enstrophy/dissipation | \(s\) |
| global \(L^3\) upper scale | \(s^{-1/2}\) |
| vorticity coefficient scale | \(s^{1/2}\) |
| signed nonlinear flux | \(s^{-1/2}\) |
| pressure/nonlinear \(L^2\) | \(s^0\) |
| physical viscous time | \(s^{-2}\) |
| bandwidth | \(s\) |

Thus flux/viscosity is \(s^{-3/2}\).

### D. Feedback loop

\[
\text{equiangular parents}
\xrightarrow{P_{p-q}N=0}
\text{one-sided sum mode}
\xrightarrow{\Pi_c>0}
\text{child energy}.
\]

This is one signed step, not a self-replicating loop.

### E. Obstruction audit

The exact sign survives pressure/Leray and energy cancellation. Fixed
cardinality does not survive viscosity and is **REJECTED** at high scale.
The three wavevectors are coplanar; a five-mode connected wrapper is needed
for a genuinely rank-three graph. No ESS or finite-time conclusion follows.

### F. One-hour falsification

Use exact rational convolution for the three- and five-mode graphs, then a
Taylor expansion through \(h=1/4\). Kill a proposed wrapper if off-graph
energy reaches half the child gain or the phase cone closes before the child
reaches its prescribed energy.

### G. Proof chain

Exact kernel → rank-three wrapper → mode-cloud amplification → phase cone →
interval energy budget → infinite-band limit → actual \(L^3\) lower bound →
Clay connection.

## 2. Adjoint Leray-response orbit

**Label: NUMERICAL CANDIDATE / SYMBOLIC CANDIDATE**

### A. Clay target

Clay (D), periodic, fixed viscosity, finite-Fourier datum plus fixed low-mode
forcing. The high-frequency response is entirely nonlinear.

### B. Central equations

For parent band \(P_j\) and child band \(C_{j+1}\), define

\[
g_j(v)=-P_{C_{j+1}}\mathbb P((v\cdot\nabla)v),
\tag{2.1}
\]

\[
\mathcal R_j(v)=
\sqrt{2E_{j+1}}\,{g_j(v)\over\|g_j(v)\|_2}.
\tag{2.2}
\]

Then

\[
\langle\mathcal R_j(v),g_j(v)\rangle
=\sqrt{2E_{j+1}}\|g_j(v)\|_2>0.
\tag{2.3}
\]

Search for a relative-periodic projective orbit

\[
\mathcal S\mathcal R_{j+L-1}\cdots\mathcal R_j(v_j)
=e^{i\phi}Rv_j+r_j.
\tag{2.4}
\]

### C. Scaling

For a localized \(N^3\)-mode block with \(E_N=N^{-1}\):

| quantity | \(N\)-scale | \(\tau\)-scale if \(N=\tau^{-1/2}\) |
|---|---:|---:|
| energy | \(N^{-1}\) | \(\tau^{1/2}\) |
| enstrophy/dissipation | \(N\) | \(\tau^{-1/2}\) |
| \(L^3{}^3\) per core | \(1\) | \(1\) |
| vorticity | \(N^2\) | \(\tau^{-1}\) |
| nonlinear/pressure \(L^2\) | \(N^{3/2}\) | \(\tau^{-3/4}\) |
| signed flux capacity | \(N\) | \(\tau^{-1/2}\) |
| physical step | \(N^{-2}\) | \(\tau\) |
| bandwidth | \(N\) | \(\tau^{-1/2}\) |

### D. Feedback loop

\[
v_j\to g_j(v_j)\to
\text{optimally aligned }v_{j+1}\to
g_{j+1}(v_{j+1}).
\]

The first arrow has an exact sign. Recurrence of shape and magnitude is the
unknown feedback.

### E. Obstruction audit

It is not a fixed-point self-similar ansatz and can be log-periodic. The
implemented simple Fejer orbit loses normalized injection and has low forcing
larger than child forcing; that simple family is **REJECTED**. A more general
phase code remains open.

### F. One-hour falsification

Iterate the dealiased response map for three doublings. Success requires
normalized injection and child forcing fraction bounded below and off-chain
ratio below one. Two consecutive decay factors below \(1/2\) kill the chosen
template.

### G. Proof chain

Finite response orbit → rationalization → invariant phase cone → interval
stage budget → summable errors → smooth pre-seeds → PDE limit → Clay.

## 3. Phase-coded coherent \(N^3\)-mode cloud

**Label: SYMBOLIC CANDIDATE / AUDIT REQUIRED**

### A. Clay target

Clay (D) on \(\mathbb T^3\), with the same force restrictions as Route 2.
The complete specification is saved in the phase-coded Leray cloud candidate
document under docs/candidates.

### B. Central equations

With \(M_N\asymp N^3\), \(E_N=c_E/N\), each coefficient is \(N^{-2}\).
The support must have fixed **relative** thickness \(\rho N\).  Absolute
thickness \(O(1)\) contains only \(O(N^2)\) lattice modes and fails the
multiplicity requirement.
Define

\[
\chi_N={\Pi_N\over N\sqrt{M_N}E_N^{3/2}}.
\tag{3.1}
\]

For phases \(\theta(q)\), low and high outputs are respectively weighted
autocorrelation and convolution:

\[
\mathcal C(\ell)=\sum_q
e^{i(\theta(q+\ell)-\theta(q))}K_-(q,\ell),
\tag{3.2}
\]

\[
\mathcal H(\ell)=\sum_q
e^{i(\theta(q)+\theta(\ell-q))}K_+(q,\ell).
\tag{3.3}
\]

Required:

\[
\inf_N\chi_N>0,\qquad
\|\mathcal C\|\le\varepsilon\|\mathcal H\|,\quad\varepsilon<1.
\tag{3.4}
\]

### C. Scaling

| quantity | scale for \(N=\tau^{-1/2}\) |
|---|---:|
| total energy | \(O(1)\) |
| enstrophy/dissipation | \(\tau^{-1/2}\) |
| global \(L^3{}^3\) | \(\log(1/\tau)\) |
| vorticity | \(\tau^{-1}\) |
| nonlinear/pressure \(L^2\) | \(\tau^{-3/4}\) |
| signed flux | \(\tau^{-1/2}\) |
| physical time | \(\tau\) |
| bandwidth | \(\tau^{-1/2}\) |
| active modes | \(\tau^{-3/2}\) |

### D. Feedback loop

\[
\text{phase autocorrelation cancellation}
\to\text{coherent forward convolution}
\to\text{child cloud}
\to\text{renormalized phase code}.
\]

### E. Obstruction audit

The \(N^3\) multiplicity is the minimum compatible with viscosity. Independent
random phases have only a heuristic \(N^{-1/2}\) baseline and require a
separate finite-size test. For the requested mesoscopic width
\(W=N^\gamma\), critical normalization and an empty child give the exact
frozen-response identity

\[
{E_C(\tau N^{-2})\over E_P(0)}
={2c_E\over N}H_NG_N^2.
\tag{3.5}
\]

More strongly, set
\(M_N^{\rm eff}=(\sum_k|\widehat u_k|)^2/\|u\|_2^2\le M_N\).  For
\(|k|\le\kappa N\) on the support,

\[
{E_C(\tau N^{-2})\over E_P(0)}
\le2\kappa^2\tau^2c_E\,{M_N^{\rm eff}\over N^3}
\le2\kappa^2\tau^2c_E\,{M_N\over N^3}.
\tag{3.6}
\]

This is a phase-independent Bernstein/heat-semigroup bound. A relay therefore
requires amplitude-effective delocalisation \(M_N^{\rm eff}\gtrsim N^3\).
Hence every \(M_N=o(N^3)\), not only the affine ansatz, fails this one-stage
critical empty-child test. Capacity saturation gives the matching exponent
\(N^{3\gamma-3}\to0\) for every \(\gamma<1\). Therefore \(G_N\) growth at
\(\gamma>2/3\) is not a success condition; only fixed-relative width
\(W=\eta N\), \(0<\eta<1/3\), survives this screen. Exact compact support is
not claimed. A filled box has an additional shape obstruction: convolution
broadens width \(W\) to \(2W\), so its child-spill/core ratio need not decay.
Only \(\gamma=1\) matches \(W_{2N}\sim2W_N\) in a doubling relay. The
surviving limit is the global-overlap-aware continuum operator
\(\mathfrak T_2\) in candidate equations (6.4)--(6.10), whose projective
periodic points replace independent mode-phase optimization. Type-I,
weak-\(L^3\), CKN, and parabolic-tail audits remain.

### F. One-hour falsification

Discretize the global continuum response at two resolutions and optimize the
two-step soft minimum together with profile covariance. Kill if the selected
child energy stays below \(1/2\), the projective two-step covariance decreases
under refinement, global-combination leakage exceeds selected energy, or
independent sparse and zero-padded convolutions disagree.

### G. Proof chain

Exact cloud family → uniform coherence → phase cone → interval budget →
localization → smooth force/data → infinite solution → \(L^3\) divergence →
Clay (D).

## 4. Perron multi-child relay

**Label: FORMAL ANSATZ**

### A. Clay target

Clay (D), periodic and fixed viscosity. Unlike Route 3, this route recruits
secondary outputs into finitely many child types rather than cancelling them.

### B. Central equations

Partition outputs into types \(C^\alpha\). For normalized shapes \(v^\alpha\),

\[
T_{\alpha\beta\gamma}
=-\langle v^\alpha,
P_{C^\alpha}\mathbb P((v^\beta\cdot\nabla)v^\gamma)\rangle.
\tag{4.1}
\]

For amplitudes \(z_\alpha\),

\[
\dot E_\alpha=
\sum_{\beta,\gamma}T_{\alpha\beta\gamma}z_\alpha z_\beta z_\gamma
-2\nu d_\alpha N^2E_\alpha.
\tag{4.2}
\]

Seek a donor/child cone \(K\) and positive vector \(e\) with

\[
\mathscr P(T[z,z])\ge\lambda e,\qquad z\in K,\quad\lambda>0.
\tag{4.3}
\]

### C. Scaling

With finitely many types but \(N^3\) modes in each type, energy is \(N^{-1}\),
dissipation and flux are \(N\), nonlinear/pressure \(L^2\) is \(N^{3/2}\),
the step is \(N^{-2}\), and \(N=\tau^{-1/2}\).

### D. Feedback loop

\[
\text{parent type vector}
\to\text{full output vector}
\to\text{positive Perron cone}
\to\text{next type vector}.
\]

### E. Obstruction audit

Energy cancellation forbids every component from growing simultaneously;
donors must be explicit. Branch count cannot outrun physical packing.
Unpopulated off modes begin with quadratic energy in time, not instantaneous
energy loss.

### F. One-hour falsification

Build the complete transfer tensor for 2–6 types and solve a cone feasibility
problem. Kill if every cone violates energy cancellation or loses more than
half its child energy within one step.

### G. Proof chain

Exact tensor → donor/child cone → interval invariance → uniform scale family →
front law → PDE limit → Clay.

## 5. Modal radial-growth variance monotone

**Label: PROOF CANDIDATE / AUDIT REQUIRED**

### A. Clay target

Clay (B), unforced \(\mathbb T^3\), arbitrary smooth mean-zero divergence-free
datum and fixed \(\nu>0\). The intended output is global regularity if a new
modal action can be bounded universally.

### B. Central equations and derivation

Let

\[
x_k=|k|^2,\quad e_k=|\widehat u_k|^2,\quad
a_k=\Re(\overline{\widehat u_k}\cdot\mathcal N_k),
\tag{5.1}
\]

\[
H_r=\sum_kx_k^re_k,\qquad T_r=\sum_kx_k^ra_k,\qquad
N_r^2={H_{r+1}\over H_r}.
\tag{5.2}
\]

Then

\[
{1\over2}\dot H_r=T_r-\nu H_{r+1}.
\tag{5.3}
\]

Define

\[
p_{r,k}={x_k^re_k\over H_r},\quad
\mu_r=\mathbb E_{p_r}x=N_r^2,\quad
g_k={a_k\over e_k}
\tag{5.4}
\]

on occupied modes. Direct division gives the exact identity

\[
{T_{r+1}\over H_{r+1}}-{T_r\over H_r}
={\operatorname{Cov}_{p_r}(x,g)\over\mu_r}.
\tag{5.5}
\]

Set

\[
\delta_r={H_{r+2}H_r\over H_{r+1}^2}-1
={\operatorname{Var}_{p_r}(x)\over\mu_r^2},
\qquad
\sigma_r^2=\operatorname{Var}_{p_r}(g).
\tag{5.6}
\]

Then

\[
{d\over dt}\log N_r
={\operatorname{Cov}(x,g)\over\mu_r}-\nu\mu_r\delta_r
\le \sigma_r\sqrt{\delta_r}-\nu N_r^2\delta_r,
\tag{5.7}
\]

and square completion gives

\[
\boxed{
{d\over dt}\log N_r\le{\sigma_r^2\over4\nu N_r^2}.}
\tag{5.8}
\]

Thus

\[
\boxed{
\mathcal M_r(t)=\log N_r(t)
-\int_0^t{\sigma_r(s)^2\over4\nu N_r(s)^2}\,ds}
\tag{5.9}
\]

is nonincreasing.

### C. Scaling

For a dyadic critical wake \(E_j\asymp\lambda_j^{-1}\),
\(1\le\lambda_j\le N=\tau^{-1/2}\), and relative front-growth dispersion
\(\varepsilon\):

| quantity | scale |
|---|---:|
| energy | \(O(1)\) |
| enstrophy/dissipation | \(\tau^{-1/2}\) |
| \(L^3{}^3\) | \(\log(1/\tau)\) in a wake |
| vorticity | \(\tau^{-1}\) |
| nonlinear/pressure \(L^2\) | \(\tau^{-3/4}\) |
| bandwidth factors | \(N_0\asymp N^{1/2}\), \(N_1,N_2\asymp N\) |
| \(\sigma_0\) | \(\varepsilon N^{3/2}=\varepsilon\tau^{-3/4}\) |
| \(\sigma_1,\sigma_2\) | \(\varepsilon N^2=\varepsilon\tau^{-1}\) |
| each action density | \(O(\varepsilon^2N^2/\nu)=O(\varepsilon^2/(\nu\tau))\) |
| physical time | \(\tau\) |
| bandwidth | \(\tau^{-1/2}\) |

### D. Regularity loop

If the actions are integrable for \(r=0,1,2\), all \(N_r\) are bounded and

\[
H_3=H_0N_0^2N_1^2N_2^2
\tag{5.10}
\]

is bounded. Standard \(H^3\) continuation then gives regularity. Conversely,
a singular solution must make at least one action nonintegrable.

### E. Obstruction audit

The monotone formula is exact and does not assume the old Route-10 barrier.
The missing theorem is a universal action estimate. Tiny newly born modes can
make \(g_k=a_k/e_k\) large.  Put \(\bar g_r=T_r/H_r\), sum only over
\(e_k>0\), and assign zero contribution when \(e_k=a_k=0\).  The exact
variance, which avoids naming \(g_k\) but still divides by \(e_k\), is

\[
\sigma_r^2={1\over H_r}
\sum_{e_k>0}x_k^r{(a_k-\bar g_re_k)^2\over e_k}.
\tag{5.11}
\]

A genuinely division-free sufficient upper bound is

\[
\boxed{
\sigma_r^2\le {1\over H_r}\sum_kx_k^r
|\mathcal N_k-\bar g_r\widehat u_k|^2.}
\tag{5.12}
\]

The formulas are first finite-sum identities.  Passing to the full PDE uses
spectral cutoffs on each classical interval \([0,T']\), then
\(T'\uparrow T\); zero crossings and infinite second moments use the
extended-value convention.  In particular,

\[
N_r(t)\le N_r(0)\exp\!\left(
\int_0^t{\sigma_r^2\over4\nu N_r^2}\,ds\right).
\tag{5.13}
\]

### F. One-hour falsification

On dealiased 32, 48, and 64 grids, compare (5.5) with the direct front
difference, verify the square-completion slack, and test cutoff convergence
under phase scrambling. Kill the practical diagnostic if tiny modes dominate
without resolution convergence.

### G. Proof chain

Exact identity → robust zero-mode convention → universal action bound for
\(r=0,1,2\) → bounded bandwidth factors → bounded \(H^3\) → Clay (B).

## 6. \(H^3\)-quantile front action

**Label: PROOF CANDIDATE / AUDIT REQUIRED**

### A. Clay target

Clay (B), unforced periodic data. This is a regularity route designed not to
miss a front whose ordinary energy tends to zero.

### B. Central equations

Use the strict cutoff

\[
\phi_m(z)={z^m\over1+z^m},\qquad m\ge8\text{ a fixed even integer},
\tag{6.1}
\]

\[
H_{r,K}=\sum_k\phi_m(|k|/K)x_k^re_k,\qquad
T_{r,K}=\sum_k\phi_m(|k|/K)x_k^ra_k.
\tag{6.2}
\]

For a nonzero mean-zero smooth solution and fixed
\(\theta\in[\theta_0,1-\theta_0]\), \(0<\theta_0<1/2\), define the unique
\(K_{r,\theta}>0\) by

\[
H_{r,K_{r,\theta}}=\theta H_r.
\tag{6.3}
\]

The front density is strictly positive:

\[
B_{r,K}={m\over2}\sum_k
\phi_m(1-\phi_m)x_k^re_k>0.
\tag{6.4}
\]

As \(K\) increases, the left side of (6.3) is continuous and strictly
decreases from \(H_r\) to zero, with
\(\partial_{\log K}H_{r,K}=-2B_{r,K}<0\).  Thus the implicit-function theorem
applies on every classical interval.  Differentiating (6.3), with all cutoff
moments evaluated at the same current \(K=K_{r,\theta}\), gives the exact
front ODE

\[
\boxed{
\dot{\log K}_{r,\theta}
={T_{r,K}-\theta T_r
-\nu(H_{r+1,K}-\theta H_{r+1})
\over B_{r,K}}.}
\tag{6.5}
\]

Let \(\mathfrak Q_{r,\theta}\) be the positive part of the right-hand side.

### C. Scaling

| front | \(K\) | \(\dot{\log K}\) | energy/enstrophy note |
|---|---:|---:|---|
| \(K=\tau^{-\gamma}\) | \(\tau^{-\gamma}\) | \(\gamma/\tau\) | arbitrary power front |
| Zeno | \(\tau^{-1/2}\) | \(1/(2\tau)\) | finite dissipation |
| log endpoint | \(\tau^{-1/2}\mathcal L^{-a/2}\) | \((2\tau)^{-1}(1-a/\mathcal L)\) | vanishing front energy |

For a core \(U=\tau^{-a},L=\tau^b\), energy, enstrophy, \(L^3{}^3\),
vorticity, nonlinear \(L^2\), pressure \(L^2\), and bandwidth scale as
\(\tau^{3b-2a},\tau^{b-2a},\tau^{3b-3a},
\tau^{-a-b},\tau^{b/2-2a},\tau^{b/2-2a},\tau^{-b}\).

### D. Regularity loop

At \(r=3\), put

\[
C_m:=\sup_{z\ge0}{z^6\over1+z^m}<\infty.
\tag{6.6a}
\]

Then

\[
(1-\theta)H_3
\le C_mK_{3,\theta}^6H_0.
\tag{6.6}
\]

Moreover

\[
\log K(t)\le\log K(0)+\int_0^t\mathfrak Q_{3,\theta}(s)\,ds.
\tag{6.6b}
\]

Therefore

\[
\int_0^T\mathfrak Q_{3,\theta}dt<\infty
\Longrightarrow K_{3,\theta}\text{ bounded}
\Longrightarrow H_3\text{ bounded}
\Longrightarrow\text{regularity}.
\tag{6.7}
\]

### E. Obstruction audit

An energy quantile \(r=0\) would miss the critical wake; \(r=3\) is essential.
Finite Galerkin saturation is not a PDE conclusion. A compact cutoff would
create a separate spectral-jump branch; the strict sigmoid removes exact
\(B=0\) for nonzero solutions but gives no quantitative lower bound on
\(B/H_r\). A lacunary spectrum can make that ratio arbitrarily small, so
front-density collapse remains inside the action.  The finite-sum identity
passes to the PDE by spectral cutoff on \([0,T']\) before \(T'\uparrow T\).

### F. One-hour falsification

Root-solve \(K_{3,\theta}\) at \(\theta=1/4,1/2,3/4\), compare (6.5) with a
time-refined finite difference, and inject a lacunary manufactured spectrum.
Kill if the identity fails to converge or benign smooth states drive
\(B/H_3\) to the floating-point floor.

### G. Proof chain

Exact implicit ODE → universal positive-front action bound → bounded quantile
→ bounded \(H^3\) → Clay (B).

## 7. Material stretching versus Fourier slippage

**Label: FORMAL ANSATZ**

### A. Clay target

Clay (A) or (B), unforced, with either whole-space Littlewood–Paley blocks or
periodic Fourier blocks.

### B. Central equations

Apply \(P_j\) to

\[
D_t\omega=\omega\cdot\nabla u+\nu\Delta\omega.
\tag{7.1}
\]

Because transport and projection do not commute,

\[
{1\over2}{d\over dt}\|\omega_j\|_2^2
+\nu\|\nabla\omega_j\|_2^2=S_j+C_j,
\tag{7.2}
\]

\[
S_j=\langle\omega_j,P_j(\omega\cdot\nabla u)\rangle,
\quad
C_j=\langle\omega_j,[u\cdot\nabla,P_j]\omega\rangle.
\tag{7.3}
\]

### C. Scaling

For \(U=\tau^{-a},L=\tau^b,N=L^{-1}\):

| quantity | exponent |
|---|---:|
| energy | \(3b-2a\) |
| enstrophy/dissipation | \(b-2a\) |
| \(L^3{}^3\) | \(3b-3a\) |
| vorticity | \(-a-b\) |
| nonlinear/pressure \(L^2\) | \(b/2-2a\) |
| \(S_j/E_j^\omega,C_j/E_j^\omega\) | \(-a-b\) |
| physical time | \(1\) |
| bandwidth | \(-b\) |

### D. Dichotomy

If for all high \(j\),

\[
(S_j)_++(C_j)_+\le\nu\|\nabla\omega_j\|_2^2,
\tag{7.4}
\]

the high-shell enstrophy cannot grow. Hence any singularity must contain a
sequence with either stretching or slippage exceeding half the viscous term.

### E. Obstruction audit

A crude commutator bound returns to \(\|\nabla u\|_\infty\) and adds nothing
to BKM. The route requires a phase-sensitive bound. Sharp blocks may be
unstable under small spectral shifts, so smooth-block robustness is mandatory.

### F. One-hour falsification

Compute \(S_j\) and \(C_j\) separately in a dealiased run, verify that their
sum equals total shell production, then test Galilean shifts and phase
scrambles. Kill the diagnostic if the branch assignment is not stable.

### G. Proof chain

Exact commutator identity → phase-sensitive estimate on one branch → high
enstrophy bound → finite low-shell remainder → continuation.

## 8. Lagrangian deformation-inhomogeneity action

**Label: FORMAL ANSATZ**

### A. Clay target

Clay (A) or (B), unforced. The domain can be \(\mathbb R^3\) with decaying data
or \(\mathbb T^3\).

### B. Central equations

Let \(X(a,t)\) be the flow map,

\[
F=\nabla_aX,\quad G=F^{-1}F^{-T},\quad
\Omega=\omega\circ X,\quad Z=F^{-1}\Omega.
\tag{8.1}
\]

The exact viscous Cauchy-defect equation is

\[
\partial_tZ=
\nu F^{-1}\partial_{a_\alpha}
\left(G_{\alpha\beta}\partial_{a_\beta}(FZ)\right).
\tag{8.2}
\]

For \(Y=\frac12\int|Z|^2da\), integration by parts gives

\[
\dot Y+\nu D_G=\nu R_F,
\quad
D_G=\int G_{\alpha\beta}
\partial_\alpha Z\cdot\partial_\beta Z\,da,
\tag{8.3}
\]

where \(R_F\) is the sum of the three terms containing
\(\partial F\) or \(\partial F^{-T}\). Define

\[
\mathfrak J_F=\nu{(R_F-D_G)_+\over Y}.
\tag{8.4}
\]

### C. Scaling

For label wavelength \(\ell_Z\),

\[
{|R_F|\over D_G}
\sim\kappa(F)\ell_Z|F^{-1}\nabla_aF|
+\kappa(F)^2\ell_Z^2|F^{-1}\nabla_aF|^2.
\tag{8.5}
\]

The physical core scaling table is the generic table of Route 7. The new
dimensionless variable is
\(\ell_Z|F^{-1}\nabla_aF|\); time is \(\tau\), physical bandwidth is
\(\tau^{-b}\), and label bandwidth may differ by the singular values of \(F\).

### D. Regularity loop

\[
\sup_{t<T}\|F\|_\infty<\infty,\qquad
\int_0^T\mathfrak J_Fdt<\infty
\tag{8.6}
\]

implies bounded \(Y\) and

\[
\|\omega\|_2\le\|F\|_\infty(2Y)^{1/2}.
\tag{8.7}
\]

Bounded \(H^1\) velocity gives a Serrin-class continuation. A singularity must
therefore have unbounded deformation or nonintegrable material-wrinkling
action.

### E. Obstruction audit

No Jacobian collapse is possible because \(\det F=1\). If \(F\) is spatially
constant in label space, \(R_F=0\) and the action decays; coherent affine
strain alone cannot use this route.

### F. One-hour falsification

Track particles and \(F_t=(\nabla u\circ X)F\) for one turnover. Verify (8.3)
and a manufactured constant-\(F\) case. Kill if benign flows make
\(\mathfrak J_F\) diverge under particle refinement.

### G. Proof chain

Exact material identity → quantitative commutator estimate → action
integrability or deformation alternative → bounded enstrophy → continuation.

## 9. Curved slow-axis 2D3C cigar

**Label: FORMAL ANSATZ**

### A. Clay target

Clay (C), \(\mathbb R^3\), unforced, fixed viscosity, with a divergence-free
Schwartz datum. The axisymmetric version is rejected; only a curved and
twisted Cartesian cigar remains.

### B. Central equations

Let

\[
L_\perp=N^{-1},\qquad L_\parallel=N^{-4/5},\qquad U=N^{5/4},
\tag{9.1}
\]

and \(y=(x_1/L_\perp,x_2/L_\perp,x_3/L_\parallel)\). From a Schwartz vector
potential

\[
A_N=UL_\perp(a_1(y),a_2(y),\psi(y)),
\tag{9.2}
\]

define \(u=\nabla\times A_N\). With
\(\varepsilon=L_\perp/L_\parallel=N^{-1/5}\),

\[
{u\over U}=
(\partial_2\psi-\varepsilon\partial_3a_2,\,
\varepsilon\partial_3a_1-\partial_1\psi,\,
\partial_1a_2-\partial_2a_1).
\tag{9.3}
\]

The physical 3D divergence is identically zero. Require the leading 2D3C
profile \(h=(\partial_2\psi,-\partial_1\psi)\), \(q=\partial_1a_2-\partial_2a_1\)
to satisfy

\[
h\cdot\nabla_\perp h+\nabla_\perp\pi=0,\qquad
h\cdot\nabla_\perp q=0.
\tag{9.4}
\]

This cancels the fast transverse Euler term. The surviving slow coupling
\(U^2/L_\parallel\) balances \(U/\tau\) when

\[
N=\tau^{-20/41}.
\tag{9.5}
\]

### C. Scaling

| quantity | \(\tau\)-scale |
|---|---:|
| energy | \(\tau^{6/41}\) |
| enstrophy/dissipation | \(\tau^{-34/41}\) |
| global \(L^3{}^3\) | \(\tau^{-19/41}\) |
| vorticity | \(\tau^{-45/41}\) |
| uncancelled transverse nonlinear/pressure \(L^2\) | \(\tau^{-42/41}\) |
| slow residual/time derivative \(L^2\) | \(\tau^{-38/41}\) |
| physical time | \(\tau\) |
| Fourier bandwidth | \(\tau^{-20/41}\) |

Total dissipation is finite; the BKM integral diverges.

### D. Feedback loop

\[
\text{curvature/twist}
\to\text{slow-axis transport}
\to\text{transverse concentration}
\to\text{stronger curvature coupling}.
\]

### E. Obstruction audit

Anisotropic scaling is not an NS symmetry. Any failure of (9.4) dominates the
desired residual by \(\tau^{-4/41}\). The straight axisymmetric version would
have \(\Gamma=ru^\theta\sim N^{1/4}\), contradicting the swirl maximum
principle, and is **REJECTED**.

### F. One-hour falsification

At \(N=8,16,32\), compute the free-space Leray RHS and project it onto the
ansatz tangent. Kill if the normalized tangent coefficient is nonpositive,
the normal residual does not decrease, or the exact initial \(L^3\) production
is nonpositive.

### G. Proof chain

Localized profile solving (9.4) → curved-coordinate expansion → positive
modulation ODE → residual contraction → Schwartz PDE orbit → \(L^3\)
divergence → Clay (C).

## 10. Realizable near-field pressure cage

**Label: FORMAL ANSATZ**

### A. Clay target

Clay (C), \(\mathbb R^3\), unforced and fixed viscosity, with compactly
generated divergence-free transmitter and receiver fields.

### B. Central equations

Let

\[
\psi=b(x_1/L_1)b(x_2/L_2)b(x_3/L_3),\qquad
W=\nabla\times(\psi e_3).
\tag{10.1}
\]

Its energy moment is physically realizable:

\[
M(W)=\operatorname{diag}(m_1,m_2,0),\qquad
{m_1\over m_2}=\left({L_1\over L_2}\right)^2.
\tag{10.2}
\]

For transmitter direction \(e_1\) and receiver direction \(\xi=e_2\), the
far-field sign diagnostic is

\[
-4\pi d^5\,\xi^T\nabla^2p\,\xi=12m_1-9m_2>0
\tag{10.3}
\]

when \(m_1>3m_2/4\). The pilot must nevertheless solve the exact free-space
pressure, not substitute (10.3) as a near-field premise.

Place two transmitters at \(\pm c\ell e_1\) around a compact receiver whose
core strain is \(S=\operatorname{diag}(-s/2,s,-s/2)\). Along its vorticity
direction,

\[
D_t\alpha=|\eta|^2-\alpha^2-\xi^T(\nabla^2p)\xi+R_\nu.
\tag{10.4}
\]

### C. Scaling

With \(\ell=\tau^{5/11}\), \(U=\tau^{-6/11}\):

| quantity | \(\tau\)-scale |
|---|---:|
| energy | \(\tau^{3/11}\) |
| enstrophy/dissipation | \(\tau^{-7/11}\) |
| global \(L^3{}^3\) | \(\tau^{-3/11}\) |
| vorticity/strain | \(\tau^{-1}\) |
| nonlinear/pressure-gradient \(L^2\) | \(\tau^{-19/22}\) |
| pressure Hessian | \(\tau^{-2}\) |
| physical time | \(\tau\) |
| bandwidth | \(\tau^{-5/11}\) |

### D. Feedback loop

\[
\text{transmitter Reynolds stress}
\to-\xi^T\nabla^2p\,\xi>0
\to\text{receiver strain}
\to\text{vorticity amplification}
\to\text{stronger pressure cage}.
\]

### E. Obstruction audit

No arbitrary rank-one moment is used. Increasing separation improves a
multipole remainder but weakens the effect as \(c^{-5}\). At \(c=O(1)\), all
self, cross, and rotation terms must come from exact pressure. A favorable
Hessian alone does not imply positive global \(L^3\) production.

### F. One-hour falsification

Use a padded 128-cube free-space solve over aspect ratio, separation, and
rotation. Require a positive Riccati margin on a whole receiver ball and
positive exact \(L^3\) production. Kill if the best margin is nonpositive or
the vorticity direction exits before one turnover.

### G. Proof chain

Explicit curl fields → exact pressure sign on a ball → invariant orientation
cone → Riccati lower bound → viscosity/tail control → PDE orbit → Clay (C).

## 11. Viscous geometric-phase strain ratchet

**Label: FORMAL ANSATZ**

### A. Clay target

Clay (C), unforced whole space. A forced high-frequency pulse train is not
allowed because it would violate Clay smooth-force decay.

### B. Central equations

Along a trajectory,

\[
\dot F=(\nabla u\circ X)F,\qquad \det F=1,
\tag{11.1}
\]

\[
{d\over dt}\left(F^{-1}\omega(X,t)\right)
=\nu F^{-1}\Delta\omega(X,t).
\tag{11.2}
\]

Use two noncommuting trace-free strains

\[
S_1=\begin{pmatrix}1&0&0\\0&-1&0\\0&0&0\end{pmatrix},
\qquad
S_2=\begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix}.
\tag{11.3}
\]

The zero-mean four-pulse commutator has monodromy

\[
\mathcal M(a)=
e^{-aS_2}e^{-aS_1}e^{aS_2}e^{aS_1}.
\tag{11.4}
\]

At \(a=1\), its spectral radius is approximately \(5.446\), despite zero
mean strain.

### C. Scaling

Use the same \(\ell=\tau^{5/11},U=\tau^{-6/11}\) scale as Route 10:

| quantity | \(\tau\)-scale |
|---|---:|
| energy | \(\tau^{3/11}\) |
| enstrophy/dissipation | \(\tau^{-7/11}\) |
| \(L^3{}^3\) | \(\tau^{-3/11}\) |
| vorticity | \(\tau^{-1}\) |
| nonlinear/pressure \(L^2\) | \(\tau^{-19/22}\) |
| pulse time | \(\ell/U=\tau\) |
| viscous defect fraction | \(\nu\tau/\ell^2=\nu\tau^{1/11}\) |
| bandwidth | \(\tau^{-5/11}\) |

### D. Feedback loop

\[
\text{noncommuting strain cycle}
\to\text{geometric monodromy gain}
\to\text{vorticity alignment}
\to\text{stronger next strain cycle}.
\]

### E. Obstruction audit

\(\det F=1\) forbids Jacobian collapse or particle coalescence. Shrinking
cores require material renewal and boundary flux. Kinematic monodromy does not
show that NS pressure self-generates the pulses.

### F. One-hour falsification

Run one compact strain/receiver cycle. Measure monodromy gain, the exact
viscous defect in (11.2), material residence time, pressure-Hessian mismatch,
and integrated \(L^3\) production. Kill if gain is at most one or the receiver
leaves before four pulses.

### G. Proof chain

Self-generated pulse cell → invariant material residence → viscous-corrected
gain → periodic scale orbit → finite physical time → critical norm divergence
→ Clay (C).

## 12. Collapsing vortex necklace

**Label: FORMAL ANSATZ**

### A. Clay target

Clay (C), unforced \(\mathbb R^3\), fixed viscosity, from a smooth
divergence-free finite-energy datum. No ancestral multiscale wake is used.

### B. Central equations

Place \(K\) rotated curl-built beads on a ring:

\[
x_j=R(\cos(2\pi j/K),\sin(2\pi j/K),0),
\tag{12.1}
\]

\[
u=\sum_{j=1}^KU R_jW(R_j^T(x-x_j)/\ell).
\tag{12.2}
\]

The reduced target system is

\[
\dot K={1\over5}{U\over\ell},\qquad
\dot\ell=-{3\over5}{U\over K},\qquad
\dot U={3\over5}{U^2\over K\ell},
\tag{12.3}
\]

\[
\dot R=-{2\over5}U,\qquad R=K\ell.
\tag{12.4}
\]

It has

\[
R=\tau^{2/5},\quad
\ell=\tau^{3/5},\quad
U=\tau^{-3/5},\quad
K=\tau^{-1/5}.
\tag{12.5}
\]

### C. Scaling

| quantity | \(\tau\)-scale |
|---|---:|
| total energy | \(\tau^{2/5}\) |
| enstrophy/dissipation | \(\tau^{-4/5}\) |
| global \(L^3{}^3\) | \(\tau^{-1/5}\) |
| vorticity | \(\tau^{-6/5}\) |
| nonlinear/pressure \(L^2\) | \(\tau^{-1}\) |
| physical time | \(\tau\) |
| bead bandwidth | \(\tau^{-3/5}\) |
| center count | \(\tau^{-1/5}\) |

### D. Feedback loop

\[
\text{ring contraction}
\to\text{more bead interactions}
\to\text{angular mode creation}
\to\text{stronger contraction}.
\]

### E. Obstruction audit

The earlier guess \(R=\tau^{1/4},\ell=\tau^{1/2},U=\tau^{-1/2}\) is
**REJECTED** because \(|\dot R|>U\). The repaired scale satisfies material
speed but requires pressure and viscosity to cancel to relative \(1/K\).
\(K(t)\) is not literally a continuously varying integer; it must be an
angular spectral centroid. Viscous tails destroy exact bead disjointness.

### F. One-hour falsification

For \(K=4,8,16\), project the exact NS tangent onto
\((R,\ell,U,K_{\rm eff})\). Kill if coefficients have wrong sign, vanish with
\(K\), fail the \(1/K\) cancellation, or reconnection occurs before the
envelope time.

### G. Proof chain

Explicit bead family → tangent modulation law → pressure sign → angular
spectral-centroid dynamics → tail/local-energy control → finite-time collapse
→ \(L^3\) divergence → Clay (C).

## 13. Coverage and new derived equations

| requested direction | routes |
|---|---|
| anisotropic whole-space scaling | 9 |
| multi-center concentration | 10, 12 |
| different physical/Fourier scales | 3, 8, 9 |
| positive pressure feedback | 10 |
| Lagrangian deformation | 8, 11 |
| stretching/diffusion synchronization | 7, 11 |
| helicity/phase subsystem | 1–4 |
| finite datum to infinite bandwidth | 1–4 |
| smooth low-mode control | 2–4 |
| periodic/log-periodic renormalisation | 2, 11 |
| new critical quantity | 5, 6 |
| new monotone/dichotomy | 5–8 |

Routes 1–8 start directly from general three-dimensional Cartesian Fourier,
vorticity, or Lagrangian NS. New derivations include:

1. the exact difference-cancelling triad (1.3);
2. the mode-multiplicity necessity \(M_N\gtrsim N^3\);
3. the adjoint response identity (2.3);
4. the modal covariance identity and monotone (5.5)--(5.9);
5. the \(H^3\)-quantile front ODE (6.5);
6. the stretching/slippage commutator split (7.2);
7. the deformation-inhomogeneity action (8.3)--(8.5);
8. the anisotropic cigar exponents;
9. the realizable pressure-cage sign;
10. the strain-commutator monodromy and necklace ODE.

## 14. Comparative score

Scores are 1–5 and rank research priority, not truth. Higher Clay and
connection scores mean fewer visible bridges.

| # | route | Clay | novel | closes | no-go | compute | interval | Lean | falsify | connects | total |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | modal-growth variance | 4 | 5 | 4 | 4 | 5 | 5 | 5 | 5 | 4 | **41** |
| 6 | \(H^3\) quantile front | 4 | 5 | 4 | 4 | 5 | 4 | 5 | 5 | 4 | **40** |
| 3 | phase-coded Leray cloud | 4 | 5 | 3 | 5 | 5 | 4 | 3 | 5 | 5 | **39** |
| 1 | exact equiangular relay | 3 | 5 | 4 | 3 | 5 | 5 | 5 | 5 | 3 | **38** |
| 2 | adjoint response orbit | 4 | 4 | 3 | 4 | 5 | 4 | 3 | 5 | 3 | **35** |
| 7 | stretching/slippage | 3 | 4 | 3 | 4 | 5 | 4 | 4 | 5 | 3 | **35** |
| 4 | Perron multi-child | 3 | 4 | 3 | 4 | 4 | 4 | 3 | 5 | 3 | **33** |
| 10 | pressure cage | 3 | 5 | 2 | 4 | 4 | 3 | 2 | 5 | 3 | **31** |
| 9 | curved 2D3C cigar | 3 | 5 | 2 | 5 | 3 | 3 | 2 | 5 | 3 | **31** |
| 8 | deformation action | 3 | 5 | 3 | 4 | 3 | 3 | 3 | 4 | 3 | **31** |
| 11 | strain ratchet | 3 | 5 | 2 | 4 | 3 | 3 | 2 | 5 | 2 | **29** |
| 12 | vortex necklace | 2 | 5 | 2 | 4 | 4 | 3 | 2 | 5 | 2 | **29** |

The top three are Route 5, Route 6, and Route 3. Route 5 supplies the strongest
new exact mathematical object; Route 6 converts a front into a continuation
quantity without losing vanishing energy; Route 3 is the strongest remaining
singularity construction.

## 15. Pilot executed

The reproducible pilot is

    python -m experiments.run_leray_relay_discovery
      --config configs/leray_relay_discovery_v1.json
      --output-dir outputs/leray_relay_discovery_v1

It performs:

1. exact rational evaluation of the three-mode field and all off modes;
2. exact covariance identities for \(r=0,1,2\) and the \(H^3\) factorization;
3. a seeded complex-polarization/chirp search for two response stages;
4. a three-stage, dealiased recursive response ladder;
5. source hashes, runtime provenance, payload sidecars, and independent bundle
   verification.

The exact child sign survives. The fixed-cardinality continuation is
analytically rejected. The simple Fejer response iteration also fails the
hard-coded scale-uniformity/off-chain screen. The retained next experiment is
the fixed-relative \(N^3\)-mode boundary cloud, not another run of the rejected
simple orbit.

The review-driven mesoscopic audit is replayed separately:

    python -m experiments.run_mesoscopic_leray_cloud
      --config configs/mesoscopic_leray_cloud_v1.json
      --output-dir outputs/mesoscopic_leray_cloud_v1

It writes the \((N,\gamma)\) table, fixed-relative controls, exact finite
carrier certificate, frozen/full one-stage comparison, two-stage
carrier-shell history and plot, scaling fits only where finite geometry is
admissible, hashes, and a semantic verifier.  Its main negative conclusions
are: every \(\gamma<1\) fails the critical empty-child exponent even under
ideal coherence; the first finite two-relay gadget fails exact cross-talk; and
the resolved Galerkin grandchild is contaminated by cross channels.

## 16. Verification handoff

1. Implement weighted autocorrelation/convolution objectives (5.2)--(5.3) for
   the surviving fixed-relative family.
2. Search a larger exact carrier alphabet for a strict two-relay graph before
   assigning a joint \(J_N\) objective.
3. Compare optimized, quadratic-chirp, complementary-code, and random phases.
4. Stop if normalized coherence decays twice or low/high forcing exceeds one.
5. Rationalize only a surviving graph.
6. Enclose the full Leray tensor and all output bands with intervals.
7. Test a positive invariant phase cone over a child-filling interval.
8. In parallel, test modal action and quantile identities on benign DNS.
9. If the cloud fails, try to prove a universal bound using Routes 5–6.
10. Advance neither side to a Clay claim until the infinite-dimensional bridge
    and the relevant continuation theorem are closed.
