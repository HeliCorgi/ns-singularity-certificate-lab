# CANDIDATE SOLUTION — phase-coded coherent Leray cloud

**Status: unverified solution candidate; not yet a PDE solution candidate**

**Label: SYMBOLIC CANDIDATE / AUDIT REQUIRED**

This records a candidate mechanism, not a singularity claim. Its new hard
input is an exact rational, positive signed transfer for the true
three-dimensional Leray nonlinearity. Its unresolved step is to amplify that
finite algebra into a uniformly coherent cloud of order \(N^3\) modes.

## 1. Clay target

The primary target is Clay (D) on \(\mathbb T^3\): fixed \(\nu>0\), smooth
mean-zero divergence-free initial datum, and a \(C^\infty\) force whose Fourier
coefficients vanish above one fixed low cutoff. The force may lock the initial
low-frequency phase but may not inject energy directly into the moving
high-frequency front. If the autonomous cloud can be started without that
controller, it would be the zero-force special case of Clay (D), equivalently
a counterexample to the regularity assertion (B), not evidence for (B).

## 2. Exact signed-transfer kernel

Let

\[
p=(1,1,0),\quad q=(1,0,1),\quad c=p+q=(2,1,1),
\quad n=p\times q=(1,-1,-1).
\tag{2.1}
\]

For integer \(s\ge1\), define

\[
u_s(x)=B e_3\sin(sp\cdot x)+C e_2\cos(sq\cdot x)
       +D n\cos(sc\cdot x).
\tag{2.2}
\]

All coefficients are perpendicular to their wavevectors. With
\(\mathcal B(u,u)=\mathbb P((u\cdot\nabla)u)\), exact rational convolution
gives

\[
\begin{array}{c|c|c}
k&a_k^{\mathcal B}&b_k^{\mathcal B}\\ \hline
sp&0&-sCD\,n/2\\
sq&-sBD\,n/2&0\\
sc&-sBC\,n/3&0\\
s(c+p)&sBD\,n/2&0\\
s(c+q)&0&-sCD\,n/2.
\end{array}
\tag{2.3}
\]

The difference sideband vanishes exactly:

\[
P_{s(p-q)}\mathcal B(u_s,u_s)=0.
\tag{2.4}
\]

For \(\Pi_k=-\langle P_k\mathcal B(u_s,u_s),P_ku_s\rangle\),

\[
\Pi_p=-{sBCD\over4},\qquad
\Pi_q=-{sBCD\over4},\qquad
\boxed{\Pi_c={sBCD\over2}},
\tag{2.5}
\]

so \(\Pi_p+\Pi_q+\Pi_c=0\) exactly. The child energy, viscosity, and immediate
off-graph forcing are

\[
E_c={3D^2\over4},\qquad D_{\nu,c}=9\nu s^2D^2,
\tag{2.6}
\]

\[
\|P_{\rm off}\mathcal B(u_s,u_s)\|_2^2
={3s^2D^2\over8}(B^2+C^2).
\tag{2.7}
\]

Thus the instantaneous child margin is positive when

\[
BC>18\nu sD.
\tag{2.8}
\]

For \(s=B=C=1\), \(D=1/8\), \(\nu=1/40\),

\[
\Pi_c={1\over16},\quad
D_{\nu,c}={9\over2560},\quad
\Pi_c-D_{\nu,c}={151\over2560}>0.
\tag{2.9}
\]

These are exact fractions recomputed by the repository's rational cos/sin
convolution and Leray projector.

## 3. Why the finite graph itself is rejected

If a shell of wavenumber \(N\) contains \(M_N\) modes and energy \(E_N\),

\[
\|u_N\|_\infty\lesssim\sqrt{M_N}\|u_N\|_2,\qquad
|\Pi_N|\lesssim N\sqrt{M_N}E_N^{3/2}.
\tag{3.1}
\]

For \(E_N=c_E/N\),

\[
|\Pi_N|\lesssim C\sqrt{M_N}c_E^{3/2}N^{-1/2},
\qquad D_{\nu,N}\asymp\nu c_EN,
\tag{3.2}
\]

and hence

\[
{\Pi_N\over D_{\nu,N}}
\lesssim {C\sqrt{c_E}\sqrt{M_N}\over\nu N^{3/2}}.
\tag{3.3}
\]

A fixed-cardinality graph has ratio \(O(N^{-3/2})\) and is **REJECTED** as a
scale iteration. A necessary multiplicity law is

\[
\boxed{M_N\gtrsim {\nu^2\over c_E}N^3.}
\tag{3.4}
\]

The implicit constant includes the Bernstein/Leray constant.  Suppressing the
fixed \(c_E\) recovers the shorter \(M_N\gtrsim\nu^2N^3\) notation.

## 4. Coherent cloud ansatz

Let \(N_j=2^jN_0\). Choose finite lattice clouds
\(\Lambda_j^P,\Lambda_j^Q,\Lambda_j^C\) in disjoint
**fixed-relative-thickness** boxes or annular sectors around
\(N_j,N_j,2N_j\). Their absolute thickness is \(\eta N_j\), for fixed
sufficiently small \(\eta>0\), and

\[
|\Lambda_j^P|+|\Lambda_j^Q|+|\Lambda_j^C|
\asymp M_j\asymp\varrho_jN_j^3,
\qquad \varrho_j\asymp\eta^3.
\tag{4.1}
\]

An annulus of absolute thickness \(O(1)\) has only \(O(N_j^2)\) lattice
points and cannot meet (3.4); that thinner version is rejected.

Set

\[
u_j(x)=\sum_{k\in\Lambda_j}
\left(\widehat u_j(k)e^{ik\cdot x}
+\overline{\widehat u_j(k)}e^{-ik\cdot x}\right),
\tag{4.2}
\]

\[
k\cdot\widehat u_j(k)=0,\qquad
E_j={1\over2}\sum_{k\in\pm\Lambda_j}|\widehat u_j(k)|^2
={c_E\over N_j}.
\tag{4.3}
\]

A typical coefficient has size

\[
|\widehat u_j(k)|\asymp\sqrt{E_j/M_j}\asymp N_j^{-2}.
\tag{4.4}
\]

One quadratic pair contributes \(N_j^{-3}\) to a target coefficient. If
\(\Theta(N_j^3)\) pairs arrive with the same phase,

\[
|\widehat{\mathcal B}_j(k)|\asymp1,\qquad
\Pi_{j\to j+1}\asymp N_j^3N_j^{-2}\asymp N_j.
\tag{4.5}
\]

## 4a. Mesoscopic thickening audit

The intermediate falsification family uses

\[
W_N=\lfloor N^\gamma\rfloor,\qquad
M_N\asymp W_N^3,qquad 0<\gamma<1,
\tag{4a.1}
\]

with two angle-\(O(1)\) parent boxes around \(Np\) and \(Nq\). Their affine
phases \(x_0\cdot k+\alpha_A\) and \(x_0\cdot k+\alpha_B\) align every pair
entering a fixed sum mode. Write \(\phi_N\) for the unit-\(L^2\) shape and
normalize the physical parent by

\[
u_N=s_N\phi_N,qquad s_N^2=\|u_N\|_2^2={2c_E\over N}.
\tag{4a.2}
\]

At full Bernstein-capacity saturation,

\[
A_N^{\rm unit}
:=\|P_C\mathbb P((\phi_N\cdot\nabla)\phi_N)\|_2
\asymp N\sqrt{M_N}\asymp N^{1+3\gamma/2}.
\tag{4a.3}
\]

This is a hypothesis, not a consequence of mode count. The forcing of the
critically normalized field is instead

\[
A_N^{\rm crit}=s_N^2A_N^{\rm unit}
\asymp N^{3\gamma/2},\qquad
G_N={A_N^{\rm crit}\over N^2\|u_N\|_2^2}
={A_N^{\rm unit}\over N^2}
\asymp N^{3\gamma/2-1}.
\tag{4a.4}
\]

For an initially empty child and \(t=\tau N^{-2}\), the frozen-parent heat
response is exact mode by mode. If \(H_N\) denotes the child-forcing-weighted
mean of

\[
h_{\nu,\tau}(r)^2,\qquad
h_{\nu,\tau}(r)={1-e^{-\nu\tau r^2}\over\nu r^2},
\tag{4a.5}
\]

then

\[
\boxed{
{E_{\rm child}(\tau N^{-2})\over E_{\rm parent}(0)}
={\|v\|_2^2\over\|u_N\|_2^2}
={2c_E\over N}H_NG_N^2.}
\tag{4a.6}
\]

This implication is phase-independent, not only a capacity heuristic. Set
\(M_N^{\rm eff}=(\sum_k|\widehat u_k|)^2/\|u_N\|_2^2\le M_N\). If
\(|k|\le\kappa N\) on the parent support, then

\[
\|P_CB(u_N,u_N)\|_2
\le\kappa N\sqrt{M_N^{\rm eff}}\|u_N\|_2^2,\qquad
{E_{\rm child}\over E_{\rm parent}}
\le2\kappa^2\tau^2c_E{M_N^{\rm eff}\over N^3}
\le2\kappa^2\tau^2c_E{M_N\over N^3}.
\tag{4a.7}
\]

Consequently amplitude-effective multiplicity
\(M_N^{\rm eff}\gtrsim N^3\) is necessary; \(M_N=o(N^3)\), and hence every
\(\gamma<1\), is universally insufficient for a one-parabolic-time
empty-child relay. Ideal coherence
merely attains the matching exponent \(N^{3\gamma-3}\to0\). Thus increasing
\(G_N\) for \(\gamma>2/3\) is not a relay success. The requested mesoscopic
family is a diagnostic approach to the boundary, while the only surviving
critical-width possibility is

\[
W_N=\eta N+O(1),\qquad 0<\eta<1/3,
\tag{4a.8}
\]

for which the empty-child ratio may remain scale-independent. Every finite
fit must also enforce the conservative separation condition
\(N>3(W_N-1)\); overlapping rows are geometry failures, not asymptotic data.

The independent-random-phase heuristic predicts only

\[
\Pi_j^{\rm random}\asymp N_j^{-1/2}.
\tag{4.6}
\]

The candidate therefore needs a coherence gain \(N_j^{3/2}\) over that
heuristic baseline.  This exponent has not yet been promoted to a theorem or
an independently converged experiment.

## 5. Phase code and low-sideband cancellation

Write positive-carrier coefficients as

\[
\widehat u_j(N_je_1+q)
=N_j^{-2}w(q/N_j)h_j(q)e^{i\theta_j(q)},
\qquad h_j(q)\perp N_je_1+q.
\tag{5.1}
\]

Reality supplies the negative-carrier conjugates. The low difference branch
is a weighted autocorrelation,

\[
\mathcal C_j(\ell)=
\sum_q e^{i(\theta_j(q+\ell)-\theta_j(q))}
K_-(q,\ell)w_{q+\ell}w_q,
\tag{5.2}
\]

whereas the forward branch is a weighted convolution,

\[
\mathcal H_j(\ell)=
\sum_q e^{i(\theta_j(q)+\theta_j(\ell-q))}
K_+(q,\ell)w_qw_{\ell-q}.
\tag{5.3}
\]

The exact kernel (2.1)--(2.4) shows that \(K_-\) can be made parallel to the
difference wavevector and removed by the Leray projection without killing
\(K_+\). The cloud must satisfy

\[
\|P_{\rm low}\mathcal B(u_j,u_j)\|_2
\le\varepsilon_*\|P_{\rm child}\mathcal B(u_j,u_j)\|_2,
\qquad\varepsilon_*<1,
\tag{5.4}
\]

and

\[
\chi_j:=
{\Pi_{j\to j+1}\over N_j\sqrt{M_j}E_j^{3/2}}
\ge\chi_*>0.
\tag{5.5}
\]

Quadratic chirps, complementary phase pairs, and finite cyclic
constant-amplitude/zero-autocorrelation codes are seeds for \(\theta_j\).
The Leray tensor weights in (5.2)--(5.3) must be included; scalar sequence
autocorrelation alone is insufficient.

There is a separate support obstruction. A filled profile of scaled width
one has a quadratic convolution of scaled width two. The part outside a child
core of the original width is therefore generally an order-one shape
fraction, not \(O(W_N/N)\). Sublinear widths cannot repair this by going to
large \(N\): \(W_{2N}/W_N\to2^\gamma<2\). For a doubling cascade, only the
fixed-relative boundary \(\gamma=1\) has

\[
W_{2N}=2W_N+O(1),
\tag{5.6}
\]

so that the full convolution support can become the next parent support.
The recorded `off_core`, `child_spill`, and `outside_child_full` ratios must
therefore be kept distinct. Calling their sum \(O(N^{\gamma-1})\) is false
for this filled-box profile.

## 6. Adjoint response and log-periodic renormalisation

For a parent block \(v\), define the true Leray best response

\[
\mathcal R_j(v)=
-\sqrt{2E_{j+1}}\,
{P_{\Lambda_{j+1}}\mathcal B(v,v)
\over\|P_{\Lambda_{j+1}}\mathcal B(v,v)\|_2}.
\tag{6.1}
\]

Then

\[
\left\langle\mathcal R_j(v),
-P_{\Lambda_{j+1}}\mathcal B(v,v)\right\rangle
=\sqrt{2E_{j+1}}\,
\|P_{\Lambda_{j+1}}\mathcal B(v,v)\|_2>0.
\tag{6.2}
\]

The sign is algebraic; shape recurrence is not. A fixed point is not required.
The proposed object is an \(L\)-periodic orbit in the reality-preserving
projective phase space:

\[
\mathcal S_{j+L,j}
\mathcal R_{j+L-1}\circ\cdots\circ\mathcal R_j(v_j)
=\mathcal Q_\varphi Rv_j+r_j,\qquad
\|r_j\|\le\epsilon\|v_j\|.
\tag{6.3}
\]

Here \(\mathcal S\) renormalizes frequency and \(R\) is a cubic-lattice
symmetry.  On one chosen representative of each \(\{k,-k\}\) pair,
\((\mathcal Q_\varphi\widehat v)(k)=e^{i\varphi}\widehat v(k)\) and the
negative representative receives \(e^{-i\varphi}\); hence
\(\widehat v(-k)=\overline{\widehat v(k)}\) is preserved.  Ordinary scalar
multiplication of the full field by \(e^{i\varphi}\) is not permitted unless
\(\varphi=0\) or \(\pi\).  Invariance of the proposed orbit under this
reality map is an unproved condition, not an NS symmetry being assumed.

### 6.1 Fixed-relative continuum response operator

The surviving \(\gamma=1\) problem has an \(N\to\infty\) equation that no
longer contains the lattice scale. Let
\(\mathcal A=-\mathcal A\subset\mathbb Z^3\) be finite, assume the cells
\(a+\Omega_\eta\) are pairwise disjoint for
\(\Omega_\eta=[-\eta,\eta]^3\), and set

\[
\widehat u_N(Na+q)=N^{-2}U_a(q/N),\qquad
(a+\xi)\cdot U_a(\xi)=0,\qquad
U_{-a}(-\xi)=\overline{U_a(\xi)}.
\tag{6.4}
\]

Here the profiles are bounded and Riemann integrable, with discontinuity and
support boundaries of measure zero; smoother profiles may be used in the
eventual candidate. These hypotheses are sufficient for the displayed
Riemann-sum limit and prevent carrier cells from being counted twice.
Then \(N\|u_N\|_2^2\to\mathcal E(U)\), where

\[
\mathcal E(U)=\sum_{a\in\mathcal A}
\int_{\Omega_\eta}|U_a(\xi)|^2\,d\xi.
\tag{6.5}
\]

For an output carrier \(r\in\mathcal A+\mathcal A\), the exact discrete
Leray convolution has the Riemann-sum limit

\[
\mathcal Q_r(U,U)(\zeta)=i\mathbb P_{r+\zeta}
\sum_{a+b=r}\int
 \bigl[U_a(\xi)\cdot(b+\zeta-\xi)\bigr]
 U_b(\zeta-\xi)\,d\xi,
\tag{6.6}
\]

with the integrand zero outside its two profile supports. The
pair-sum cells \(r+2\Omega_\eta\) need not be disjoint. Consequently the
physical field must be combined at the same scaled wave number before any
norm or selection:

\[
\mathcal Q_{\rm tot}(x)=
\sum_{r\in\mathcal A+\mathcal A}\mathcal Q_r(U,U)(x-r),
\qquad
V_{\rm tot}(x)=-h_{\nu,\tau}(|x|)\mathcal Q_{\rm tot}(x).
\tag{6.7}
\]

Use \(\mathbb P_0=0\) and \(h_{\nu,\tau}(0)=\tau\). If the next scale is
\(N'=\lambda N\), require the lattice compatibilities
\(\lambda N\in\mathbb N\); a nonzero clean carrier relay will normally also
need \(\lambda a\in\mathcal A+\mathcal A\) for every selected \(a\).
The globally combined child is pulled
back to the next parent coordinates by the quadratic renormalisation map

\[
(\mathfrak T_\lambda U)_a(\xi)
=\lambda^2V_{\rm tot}\bigl(\lambda(a+\xi)\bigr).
\tag{6.8}
\]

Indeed \(N^{-2}V_{\rm tot}(\lambda(a+\xi))
=(\lambda N)^{-2}(\mathfrak T_\lambda U)_a(\xi)\), including the required
coefficient normalization. Its selected energy satisfies

\[
\mathcal E(\mathfrak T_\lambda U)
=\lambda\sum_{a\in\mathcal A}
\int_{\lambda(a+\Omega_\eta)}|V_{\rm tot}(x)|^2\,dx.
\tag{6.9}
\]

Thus a critical doubling relay needs selected child energy ratio at least
\(1/2\). On the domain
\(\mathcal E(U)>0\), \(\mathcal E(\mathfrak T_2U)>0\), shape closure is the
projective fixed/periodic-orbit problem

\[
\widetilde{\mathfrak T}_2^{\,L}U
=\mathcal Q_\varphi R U,\qquad
\widetilde{\mathfrak T}_2U
:=\sqrt{{\mathcal E(U)\over\mathcal E(\mathfrak T_2U)}}\,
\mathfrak T_2U.
\tag{6.10}
\]

On the unit-energy slice this reduces to division by
\(\sqrt{\mathcal E(\mathfrak T_2U)}\). For the doubling map, leakage consists
of the globally combined energy of \(V_{\rm tot}\) outside
\(\bigcup_{a\in\mathcal A}2(a+\Omega_\eta)\). Individual \(r\)-components
are not called orthogonal leakage channels when their supports overlap.
Equations (6.6)--(6.10), rather
than independent per-mode phase maximisation, are the next simultaneous
two-stage search problem.

A single carrier set and a doubling are not mandatory. For an
\(L\)-periodic carrier cycle
\(\mathcal A_0,\ldots,\mathcal A_{L-1}\), define the same globally combined
map \(\mathfrak T_j\) with scale ratio \(\lambda_j>1\) from stage \(j\) to
\(j+1\). The log-periodic closure equation is

\[
\widetilde{\mathfrak T}_{L-1}\circ\cdots\circ
\widetilde{\mathfrak T}_{0}U^{(0)}
=\mathcal Q_\varphi R U^{(0)},\qquad
\prod_{j=0}^{L-1}\lambda_j>1.
\tag{6.11}
\]

This multi-type version is needed when exact integer children have different
carrier geometry from their parents. Every intermediate output is still
combined globally before the next map, and lattice compatibility must hold at
each physical scale.

## 7. Interval budget and moving front

Let \(d_*>0\) be the child viscous shape coefficient, assume
\(M_j=\varrho_jN_j^3\) with \(\varrho_j\ge\varrho_->0\), and let all low, wake, and
off-chain losses be bounded by \(\ell_*N_j\). The required margin is

\[
q_*=\chi_*\sqrt{\varrho_-}\,c_E^{3/2}
-\nu d_*c_E-\ell_*>0.
\tag{7.1}
\]

Until the first time at which the child reaches \(c_E/(2N_j)\), require the
pointwise differential inequality

\[
\dot E_{j+1}(t)\ge q_*N_j.
\tag{7.2}
\]

Equivalently, a time-averaged version must be stated on every prefix ending
before that first hitting time.  Integration then gives

\[
\Delta t_j\le {c_E\over2q_*}N_j^{-2},\qquad
\sum_j\Delta t_j<\infty,
\tag{7.3}
\]

and the continuous activation envelope obeys

\[
\dot N=kN^3,\qquad N(t)\asymp(T-t)^{-1/2}.
\tag{7.4}
\]

## 8. Physical localisation and critical norm

The \(N^3\) coherently phased modes should form a physical core of diameter
\(N^{-1}\), amplitude \(N\), and volume \(N^{-3}\). The required core-tail
estimates are

\[
\|u_j\|_2^2\asymp N_j^{-1},\qquad
\int_{B_j}|u_j|^3dx\ge c_3,
\tag{8.1}
\]

\[
\|u-u_j\|_{L^3(B_j)}
\le\frac12\|u_j\|_{L^3(B_j)}.
\tag{8.2}
\]

They imply the displayed global lower bound only if the completed wake balls
\(B_j\) are pairwise disjoint at each time and each old core lower bound
persists uniformly until \(T\). Under those additional obligations,

\[
\|u(t)\|_3^3\gtrsim J(t)\asymp\log{1\over T-t}.
\tag{8.3}
\]

Exact compact support is not assumed and cannot coexist with exact finite-band
support.

## 9. Scaling table

With \(\tau=T-t\), \(N=\tau^{-1/2}\):

| quantity | cloud scale | \(\tau\)-scale |
|---|---:|---:|
| total energy | \(O(1)\) | \(O(1)\) |
| front energy | \(N^{-1}\) | \(\tau^{1/2}\) |
| enstrophy | \(N\) | \(\tau^{-1/2}\) |
| global \(L^3{}^3\) | \(\log N\) | \(\log(1/\tau)\) |
| \(\|\omega\|_\infty\) | \(N^2\) | \(\tau^{-1}\) |
| dissipation rate | \(N\) | \(\tau^{-1/2}\) |
| nonlinear \(L^2\) capacity | \(N^{3/2}\) | \(\tau^{-3/4}\) |
| pressure-gradient \(L^2\) upper scale | \(N^{3/2}\) | \(\tau^{-3/4}\) |
| physical time remaining | \(N^{-2}\) | \(\tau\) |
| Fourier bandwidth | \(N\) | \(\tau^{-1/2}\) |
| active mode count | \(N^3\) | \(\tau^{-3/2}\) |

Thus \(\int_0^T N(t)\,dt<\infty\), whereas
\(\int_0^T N(t)^2dt=\infty\).

## 10. Smooth forcing

This candidate deliberately imposes the stronger finite-band controller
condition

\[
P_{>K_0}f(t)\equiv0
\tag{10.1}
\]

for one fixed \(K_0\).  Clay (D) itself only requires a smooth periodic force
with the stated uniform derivative bounds; it does **not** require a Fourier
cutoff.  For every \(q,s\ge0\), we additionally require

\[
\sup_{t\ge0}\|\partial_t^qf(t)\|_{H^s}<\infty.
\tag{10.2}
\]

The force may select the initial phase basin but cannot repair a high-mode
residual. Making it fixed or flat near \(T\) is one sufficient way, not an
official Clay requirement, to obtain a smooth extension past \(T\).

## 11. Obstruction audit

| obstruction | candidate response |
|---|---|
| bounded energy | \(\sum_jc_E/N_j<\infty\) |
| finite dissipation | \(\int N(t)dt<\infty\) |
| ESS endpoint | actual local-core lower bound (8.3) |
| fixed finite bandwidth | \(M_N,N\to\infty\) |
| finite Galerkin ODE | every cutoff stops the cloud and is only a pilot |
| pure swirl decrease | general Cartesian 3D cloud |
| one-scale self-similar no-go | wake plus projective periodic orbit |
| smooth-force high tail | high front is autonomous |
| pressure double counting | pressure appears once through \(\mathbb P\) |
| compact/band-limited conflict | only core-tail bounds are requested |
| Type-I, weak-\(L^3\), CKN | still open |

## 12. Pilot support and retained failures

The exact three-mode fractions in (2.9) are verified. A dealiased Fejer packet
pilot implements (6.1) and confirms positive parent-parent injection,
divergence near roundoff, and nonlinear energy cancellation near roundoff.
Best-response alignment one is an algebraic consequence of the definition,
not independent evidence.

The mesoscopic angled-box pilot separately enforces affine phase, physical
Leray polarization, Hermitian symmetry, critical energy normalization, an
initially empty child, and heat-Duhamel response.  Its decisive analytic
result is the phase-independent bound (4a.7): even before asking whether the
optimistic capacity law is saturated, every
\(W_N=N^\gamma\), \(\gamma<1\), has child-energy ratio
\(O(N^{3\gamma-3})\) and is therefore **REJECTED** as a critical relay.  A
full-coefficient zero-padded local-FFT measurement at
\(N=64,\gamma=4/5,W_N=27\) gives
\(G_N=0.963601\), \(\|v\|_2/\|u\|_2=0.157979\), and
\(E_{\rm child}/E_{\rm parent}=0.0249574\).  Thus the nearly order-one
instantaneous gain still produces a small empty child.  A
small full Galerkin check at \(N=4,W=2\) gives
\(D_{\rm frozen}=2.04074\times10^{-3}\) and
\(D_{\rm full}=1.97616\times10^{-3}\); parent evolution retains about 96.8%
of the already small frozen response.

An exact two-relay carrier search also found a partial gadget but rejected it:
the intended next-stage outputs carry only \(222/2483\) of its complete
relay-cross power, diagonal cross-talk lies on the target shell, and the two
aligned grandchildren interact exactly to zero.  With the known relay fixed,
all 16 compatible orientations in the stated finite alphabet fail the strict
cross-talk screen.  This negative result is exhaustive only in that finite
alphabet.  A full two-stage Galerkin run resolves the same failure dynamically:
the intended children and diagonal cross-talk receive respectively
\(5.74749\times10^{-4}\) and \(5.68886\times10^{-4}\) of the initial parent
energy, while the contaminated grandchild receives only
\(6.35542\times10^{-8}\).

The following variants are retained as negative results:

1. fixed three- or five-mode dilation: **REJECTED** by (3.3);
2. existing P1 phase: **REJECTED** because all initial shell fluxes vanish;
3. existing P2 phase: **REJECTED** because its top occupied shell donates;
4. simple single-polarization response iteration: **REJECTED** by decreasing
   normalized injection and excessive low forcing;
5. independent random phases: retained as a **SPECULATION**-level negative
   baseline from (4.6), pending a finite-size scaling test;
6. every sublinear mesoscopic width: **REJECTED** by the empty-child exponent
   \(3\gamma-3<0\), even if \(G_N\) itself grows;
7. the first finite two-relay carrier: **REJECTED** by exact diagonal
   cross-talk and zero simple recursion.

The fixed-relative \(N^3\) boundary cloud remains unproved; exact finite
convolutions test constants and geometry but do not establish persistence.
At \(N=64\), the \(W_N=\lfloor0.2N\rfloor=12\) control has
\(E_{\rm child}/E_{\rm parent}=0.00219411\); reaching \(1/2\) by amplitude
rescaling alone would require \(c_E\simeq227.9\) in the frozen model and says
nothing about nonlinear persistence.

## 13. Unproved lemmas

1. Lattice clouds with \(M_N\asymp N^3\) satisfying (5.4)--(5.5) exist.
2. Their phase code lies in a positive invariant cone.
3. The continuum projective response map (6.10) has a finite-period orbit and
   its lattice Riemann sums have uniform error.
4. The interval budget (7.2) closes with every output band included.
5. Cloud localisation gives (8.1)--(8.2) with viscous tails.
6. Future seeds form one smooth convergent initial datum.
7. The low-band controller satisfies (10.1)--(10.2), or is unnecessary.
8. Galerkin limits give one classical solution on every \(t<T\).
9. Local energy and pressure bounds pass to the infinite-band limit.
10. Actual \(L^3\) divergence connects faithfully to Clay (D).

## 14. Shortest falsification experiment

Fix a rational relative width such as \(\eta=3/16\), choose the smallest
carrier graph surviving the exact cross-talk search, and discretize
\(\mathfrak T_2\) at three profile resolutions. Optimize both stages together.
Measure

\[
\mathcal A_1={\mathcal E(\mathfrak T_2U)\over\mathcal E(U)},\qquad
\mathcal A_2={\mathcal E(\mathfrak T_2^2U)\over
 \mathcal E(\mathfrak T_2U)},\qquad
\rho_{\rm shape}=
{|\langle U,\widetilde{\mathfrak T}_2^{\,L}U\rangle|
\over\|U\|_2\|\widetilde{\mathfrak T}_2^{\,L}U\|_2},
\tag{14.1}
\]

and

\[
\mathcal L=
{\int_{\mathbb R^3\setminus\cup_a2(a+\Omega_\eta)}
 |V_{\rm tot}(x)|^2\,dx
\over
 \int_{\cup_a2(a+\Omega_\eta)}|V_{\rm tot}(x)|^2\,dx}.
\tag{14.2}
\]

Kill this family if:

1. no exact finite carrier graph survives diagonal cross-talk;
2. \(\min(\mathcal A_1,\mathcal A_2)<1\) at every tested critical-energy
   constant that remains numerically resolvable;
3. \(\rho_{\rm shape}\) decreases under two successive refinements;
4. \(\mathcal L\ge1\);
5. the full Galerkin response loses more than half of the frozen response; or
6. independent global sparse and zero-padded convolutions disagree.

Float arithmetic is sufficient for discovery. An accepted graph must then be
rationalized and enclosed by interval complex arithmetic.

## 15. Final proof chain

1. Specify one exact phase-coded \(N^3\)-mode cloud family.
2. Prove divergence, reality, energy law, and lattice closure.
3. Prove low-sideband cancellation and uniform positive coherence.
4. Prove an invariant phase cone or finite-period response orbit.
5. Enclose the full interval flux budget.
6. Construct smooth pre-seeds and, if used, a fixed low-band force.
7. Pass finite clouds to a classical solution on every \(t<T\).
8. Prove local-core \(L^3\) lower bounds and finite total dissipation.
9. Prove \(\|u(t_n)\|_3\to\infty\) as \(t_n\uparrow T\).
10. Use uniqueness and continuity of a hypothetical global smooth solution to
    contradict that divergence and conclude Clay (D).

## 16. Present verdict

The exact sign problem for one true Leray triad is solved. The finite-mode and
simple response-map escalations are not. The decisive question is now:

> Does the global-overlap-aware continuum Leray map \(\mathfrak T_2\) admit a
> physically real projective periodic orbit with a strict two-stage energy
> margin and leakage below the selected child energy?

Until that is answered, this is a symbolic candidate rather than a PDE
candidate.
