# LENS 13 — The Moving Fourier Front $N(t)$ as a State Variable

**STATUS: FORMAL ANSATZ.** (Sub-results: the *Front–BKM identity* §D.3 and the *capacity floor* §B.5 are **PROOF CANDIDATE**, finite algebra, Lean-able. The flux closures are ansätze; no PDE theorem is claimed.)

Conventions: $\mathbb T^3$, $u=\sum_k\hat u_ke^{ik\cdot x}$, $k\in\mathbb Z^3$ (so $N$ is dimensionless),
$\partial_t\hat u_k=\mathcal N_k-\nu|k|^2\hat u_k+\hat f_k$, $\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$, $P_k=I-k\otimes k/|k|^2$, $h_{\nu,\tau}(r)=(1-e^{-\nu\tau r^2})/(\nu r^2)$. Clock $\tau=T-t$.

---

## A. Clay target

**Primary: (D)** — breakdown on $\mathbb T^3$, smooth forcing *allowed*. **Forcing: none used** ($f\equiv0$), so a success also refutes **(B)**. Forcing is declared inadmissible-as-a-crutch by F-N4 (smooth $f$ has super-polynomially decaying $\hat f_k$: it cannot feed the front) and F-N1/N2 (energy/dissipation stay bounded with forcing anyway). Initial data: $u_0\in C^\infty(\mathbb T^3)$, $\nabla\cdot u_0=0$, mean zero, band-limited to $|k|\le N_0$. Viscosity: $\nu>0$ fixed, treated **exactly** (no inviscid limit, no vanishing-$\nu$ argument). Singularity signature: $\limsup_{t\uparrow T}\|u\|_{H^m}=\infty$ certified through $\int_0^T\|\omega\|_\infty=\infty$ and $\|u\|_{L^3}\to\infty$, **never** through energy or enstrophy (dead diagnostics).

**Derived scope restriction (§E.5):** this lens *cannot* target (C)/TARGET-U through the axisymmetric $\mathbb R^3$ route at $g=1/2$, because the critical wake forces a pointwise Type-I bound and collides with CSTY2009. This is why (D)/(B) is the target and why §B.7 abandons the exactly-critical wake.

---

## B. Central mathematics: the front ODE, derived

### B.1 State variables
Wake law with one free slowly-varying amplitude $\Lambda$:
$$E_j:=\tfrac12\|u_{\lambda_j}\|_2^2=\frac{c_E\Lambda(\lambda_j)}{\lambda_j},\qquad \lambda_j=2^j\le N(t),\qquad \Lambda(\lambda)=\lambda^a,\ a\ge0\ (\;a=0\Rightarrow\Lambda\equiv1,\text{ repo-critical}).$$
Continuum density $e(k)=2c_E\Lambda(k)k^{-2}$, so $\int_{\lambda}^{2\lambda}e=c_E\Lambda(\lambda)/\lambda$ up to a fixed constant absorbed in $c_E$. Front $N(t):=\sup\{k:e(k,t)\ge\tfrac12\cdot2c_E\Lambda k^{-2}\}$: the half-amplitude leading edge. Occupation $M_N=mN^3$, $m\in(0,4\pi/3]$ ($\gamma=1$, fixed-relative width — the only regime not pre-rejected by §6 of the constraint map).

Dimensions (torus radius $1$, $k$ integer, everything a power of $T^{-1}$): $[c_E]=T^{-2}$, $[\nu]=T^{-1}$, $[\Pi]=T^{-3}$, $[c_E^{3/2}]=T^{-3}$ ✓.

### B.2 Front kinematics (exact, given the wake law)
Top-octave energy $\mathcal E_{\rm top}=2c_E\Lambda(N)/N$. Advancing the front by $dN$ at the prescribed wake level costs
$$\dot{\mathcal E}_{\rm fill}=2c_E\Lambda(N)N^{-2}\dot N \quad[\,T^{-3}\,]\ ✓ .$$
Viscous drain of the top octave: $2\nu\int_{N/2}^{N}k^2e\,dk\cdot 2 = 4\nu c_E\Lambda(N)N/2\cdot 2 \asymp 4\nu c_E\Lambda(N)N$. Band budget (advection is *exactly* energy-neutral, $\int\partial_k\Pi\,dk=0$, so all gain is flux $\Pi_N$ through the front):
$$\boxed{\,2c_E\Lambda N^{-2}\dot N=\Pi_N-4\nu c_E\Lambda N\,}\tag{FL}$$

### B.3 Capacity bound (repo L-02, re-derived)
$|\Pi_N|=|\langle P_NB(u,u),P_Nu\rangle|\le\|u_N\|_\infty\|\nabla u_N\|_2\|u_N\|_2\le\sqrt{M_N}\,\|u_N\|_2\cdot N\|u_N\|_2\cdot\|u_N\|_2$
$=N\sqrt{M_N}(2c_E\Lambda/N)^{3/2}=2\sqrt2\,c_E^{3/2}\Lambda^{3/2}\sqrt m\,N.$ Write
$$\Pi_N=\chi\cdot2\sqrt2\,\sqrt m\,c_E^{3/2}\Lambda(N)^{3/2}N,\qquad \chi\in[-1,1]\ \text{(net coherence efficiency)}.\tag{CAP}$$

### B.4 The front ODE
Insert (CAP) into (FL) and divide by $2c_E\Lambda N^{-2}$:
$$\boxed{\ \dot N=\Big(\sqrt2\,\chi\sqrt m\,\sqrt{c_E\Lambda(N)}-2\nu\Big)N^3\ }\tag{ODE}$$
Sign check: $\chi>0$ (forward transfer) $\Rightarrow$ front advances iff coherence beats viscosity. Dimension: $[\sqrt{c_E}]=T^{-1}=[\nu]$ ✓.

At $\Lambda\equiv1$ this is exactly the repo's $\dot N=kN^3$, $k=\tfrac12(X\sqrt{c_E}-2\nu)$, and reproduces seed S2's threshold $\chi\sqrt{c_E}\gtrsim\sqrt2\nu$ ✓ and seed S1's $\gamma=1/2$ ✓ — a two-way consistency check of the derivation.

### B.5 Capacity floor on the front exponent (PROOF CANDIDATE)
If $\dot N\le\kappa N^p$ then $N$ reaches $\infty$ only if $p>1$ and $N\sim[(p-1)\kappa\tau]^{-1/(p-1)}$, i.e.
$$g=\frac1{p-1}.$$
(CAP) gives $p\le3+a/2$. Hence **$g\ge\dfrac{2}{4+a}$, with equality iff the capacity is saturated.** For the critical wake $a=0$: $g\ge1/2$ — the repo's floor, here derived from *capacity*, not from BKM. Independent second derivation, and strictly stronger (it holds even where BKM is uninformative).

### B.6 The three closures — solved
Eddy time at the front: $t_N=1/(N\|u_N\|_2)=N^{-(1+a)/2}/\sqrt{2c_E}$.

**(i) Random phase.** $\Pi^{\rm RP}_N=\mathcal E_{\rm top}/t_N=\sqrt2c_E^{3/2}\Lambda^{3/2}N^{-1/2}$ (equivalently $\chi_{\rm RP}=\tfrac12 N^{-3/2}/\sqrt m$: the $1/\sqrt{M_N}$ cancellation). (FL) $\Rightarrow$
$$\dot N=\tfrac1{\sqrt2}\sqrt{c_E}\,N^{3/2+a/2}-2\nu N^3 .$$
Since $3/2+a/2<3$ for all $a<3$, **$\dot N<0$ for $N>N_*=\big(\sqrt{c_E}/(2\sqrt2\nu)\big)^{2/(3-a)}$**: the front is pinned at a finite dissipation scale. At $a=0$, $N_*^3=c_E/(8\nu^2)$, matching the Kolmogorov cutoff $(\varepsilon/\nu^3)^{1/4}$ with $\varepsilon=4\nu c_EN_*$. **REJECTED — V1.** *No random-phase front blows up, for any admissible wake.*

**(ii) Coherent capacity, critical wake ($a=0$, $\chi$ const).** $\dot N=kN^3$, $k=\sqrt2\chi\sqrt{mc_E}-2\nu$.
Blow-up iff $\chi\sqrt{mc_E}>\sqrt2\nu$; then $N=(2k\tau)^{-1/2}$, $g=1/2$ **exactly**, for *every* $\chi$ above threshold. Integrals: $\int N\,dt<\infty$ ✓; $\int N^2dt=\int(2k\tau)^{-1}d\tau=\infty$ **logarithmically** ✓ (marginal). Type I: $\sqrt\tau\|u\|_\infty\to2\sqrt{2c_E}/\sqrt{2k}$, a finite constant — **pointwise Type I**. Wake retention over remaining time: $e^{-2\nu N_j^2\tau_j}=e^{-2\nu C^2}$, a constant fraction per octave ✓ (this *is* seed S4's retention claim, now derived).
**Verdict: survives on $\mathbb T^3$, DEMOTED.** Two costs: a hard, permanent, $N$-independent coherence threshold $\chi\ge\sqrt2\nu/\sqrt{mc_E}$; and Type-I status (§E.5).

**(iii) Intermittent packet / log-critical and mildly supercritical wake ($a>0$).** $\sqrt\Lambda=N^{a/2}\to\infty$ *swallows the viscous term*: for $N>N_{\rm th}=(2\nu/(\sqrt2\chi\sqrt{mc_E}))^{2/a}$, (ODE) becomes $\dot N\simeq\sqrt2\chi\sqrt{mc_E}\,N^{3+a/2}$, hence
$$N(\tau)=\Big[\tfrac{4+a}{2}\sqrt2\chi\sqrt{mc_E}\,\tau\Big]^{-2/(4+a)},\qquad \boxed{g=\frac{2}{4+a}\in(2/5,\,1/2]\ \text{ for }a\in[0,1).}$$
The marginal case $\Lambda=\log N$ ("$a=0^+$", the repo's log-critical $\beta=-1,\sigma=\gamma$ branch) gives $N\simeq(\chi\sqrt{c_E}\tau)^{-1/2}(\tfrac12\log\tfrac1\tau)^{-1/4}$: $g=1/2$ with a $(\log)^{-1/4}$ correction and required coherence $\chi_{\rm req}=\sqrt2\nu/\sqrt{mc_E\log N}\to0$.

### B.7 Required coherence and the optimal front exponent $g^*$ — the deliverable
Capacity-saturated propagation ties $a$ and $g$: $a=(2-4g)/g$. Solving (ODE) for $\chi$:
$$\boxed{\ \chi_{\rm req}(g;N)=\frac{2\nu+\;g\,C^{-1/g}N^{1/g-2}}{\sqrt2\sqrt{m}\sqrt{c_E}\,N^{a/2}}\ \xrightarrow[\ \text{saturated}\ ]{}\ \frac{\sqrt2\,\nu}{\sqrt{mc_E}}\,N^{-a/2}=\frac{\sqrt2\,\nu}{\sqrt{mc_E}}\,N^{-(2-4g)/(2g)} }$$
* $g>1/2$ ($a<0$, sub-critical wake): $\chi_{\rm req}\sim N^{|a|/2}\to\infty>1$. **REJECTED — V2.**
* $g>1/2$ via *unsaturated tuned* capacity at $\Lambda\equiv1$: needs $\chi\sqrt{2mc_E}-2\nu=2gC^{-1/g}N^{1/g-2}\downarrow0$ at an exact rate. Any $\liminf\chi\sqrt{2mc_E}>2\nu$ forces $\frac{d}{dt}N^{-2}=$ const $<0$, i.e. $g=1/2$; any deficit reverses the front. Codimension-$\infty$; and the wake is destroyed, $e^{-2\nu C^2\tau^{1-2g}}\to0$. **REJECTED — V3 (non-generic + wake erosion).**
* $g=1/2$ ($a=0$): $\chi_{\rm req}=\sqrt2\nu/\sqrt{mc_E}$, a **positive constant floor**, the only $g$ where viscosity never becomes negligible.
* $g<1/2$ ($0<a<1$): $\chi_{\rm req}\to0$ like $N^{-a/2}$. **Viscosity ceases to be an obstruction at the front.**
* $a\ge1$: energy $\sum_jc_E2^{j(a-1)}=\infty$. **REJECTED — V4.**

**Trade-off.** Lowering $g$ below $1/2$ buys vanishing coherence demand and pays in *energy*: $E_{\rm tot}=2c_E/(1-2^{a-1})$ blows up as $a\to1^-$. Fix the energy budget $E_{\rm tot}\le E_0$ (so $c_E=E_0(1-2^{a-1})/2$) and minimise the **binding** requirement, the demand at the entry scale $N_0=2^{J_0}$:
$$\chi_{\rm req}(a)=\frac{2\nu}{\sqrt{mE_0(1-2^{a-1})}\;N_0^{a/2}} .$$
$\partial_a\log[(1-2^{a-1})N_0^{a}]=0\Rightarrow \dfrac{2^{a-1}}{1-2^{a-1}}=J_0$, hence
$$\boxed{\,a^*=1-\log_2\!\Big(1+\tfrac1{J_0}\Big),\qquad g^*=\frac{2}{5-\log_2(1+1/J_0)},\qquad \chi_{\min}=\frac{2\nu\sqrt{1+J_0}}{\sqrt{mE_0}\;N_0^{a^*/2}}\, }$$
| $N_0$ | $J_0$ | $a^*$ | $g^*$ | $\chi_{\min}\sqrt{mE_0}/2\nu$ |
|---|---|---|---|---|
| 64 | 6 | 0.7776 | **0.4186** | 0.161 |
| 1024 | 10 | 0.8625 | **0.4113** | 0.0546 |
| $\to\infty$ | | $\to1$ | $\to\mathbf{2/5}$ | $\sim\sqrt{J_0}\,N_0^{-1/2}$ |

**The program should aim at $g\approx0.41$–$0.42$ (infimum $2/5$), not at $g=1/2$.** The repo's window $1/2\le g<1$ is an artefact of imposing the *exactly*-critical wake; a mildly supercritical wake $E_N=c_EN^{a-1}$, $a\approx0.78$, opens $2/5\le g<1/2$, kills the viscous coherence threshold, and makes Type-II automatic. In repo shell notation $(\gamma,\sigma,\beta)=\big(\tfrac2{4+a},\tfrac{2(1-a)}{4+a},a-1\big)=(0.419,0.093,-0.222)$ — inside $0<\gamma<1$, $\max(0,2\gamma-1)=0\le\sigma<\gamma$, with $\beta\in(-1,0)$ (allowed after the 2026-08-01 erratum).

**Invariant (independent of $a,g$ on the saturated branch):** $\Delta t_{\rm octave}/t_N=N^{(1+a)/2-1/g}=N^{-3/2}$. The front always outruns eddy turnover by $N^{3/2}$ — transfer must be *ballistic/phase-coded*, never turbulent. This is the structural reason V1 dies.

---

## C. Scaling table ($\tau=T-t$, saturated branch, $a\in(0,1)$, $g=\frac2{4+a}$; numbers for $a=0.7776$)

| quantity | law in $N$ | exponent in $\tau$ | $a=0.778$ |
|---|---|---|---|
| Fourier front $N$ | — | $\tau^{-g}$ | $\tau^{-0.419}$ |
| energy $E$ | $2c_E/(1-2^{a-1})$ | $\tau^{0}$ (bounded) | const |
| enstrophy $Z$ | $c_EN^{1+a}$ | $\tau^{-g(1+a)}$ | $\tau^{-0.744}$ |
| dissipation rate $\varepsilon=4\nu c_EN^{1+a}$ | $N^{1+a}$ | $\tau^{-0.744}$ ($\int<\infty$ ✓) | $\tau^{-0.744}$ |
| $\|u\|_{L^3}^3$ | $\asymp N^{3a/2}$ | $\tau^{-3ag/2}$ | $\tau^{-0.489}$ |
| $\|u\|_\infty$ | $\sqrt{2c_E}N^{1+a/2}$ | $\tau^{-g(1+a/2)}$ | $\tau^{-0.581}$ |
| $\sqrt\tau\|u\|_\infty$ (Type-I test) | — | $\tau^{-0.081}\to\infty$ | **Type II** ✓ |
| $\|\omega\|_\infty$ | $\theta\sqrt{2c_E}N^{2+a/2}$ | $\tau^{-g(2+a/2)}=\tau^{-1}$ | $\tau^{-1}$ (BKM log-div.) |
| nonlinear $\|P(u\!\cdot\!\nabla u)\|_2$ | $\asymp N^{3/2+a}$ | $\tau^{-g(3/2+a)}$ | $\tau^{-0.953}$ |
| pressure $\|\nabla p\|_2$ | same order | same | $\tau^{-0.953}$ |
| time remaining $\tau$ | $\asymp N^{-(4+a)/2}$ | — | $N^{-2.39}$ |
| active mode count $M_N$ | $mN^3$ | $\tau^{-3g}$ | $\tau^{-1.256}$ |
| eddy time $t_N$ | $N^{-(1+a)/2}$ | — | — |

Note $\|\nabla p\|_2\asymp\|P(u\!\cdot\!\nabla u)\|_2$: pressure is **never** subleading. All net gain is a cancellation *pattern*, never a magnitude — consistent with $\langle x,B(x,x)\rangle=0$.

---

## D. The closed feedback loop

$$\underbrace{N}_{\text{front}}\ \xrightarrow{\ \|u_N\|_2=(2c_E\Lambda(N)/N)^{1/2}\ }\ \underbrace{\text{amplitude}}_{}\ \xrightarrow{\ \Pi_N=2\sqrt2\chi\sqrt m\,c_E^{3/2}\Lambda^{3/2}N\ }\ \underbrace{\text{flux}}_{}$$
$$\xrightarrow{\ \dot N=(\sqrt2\chi\sqrt{mc_E\Lambda}-2\nu)N^3\ }\ \underbrace{\dot N>0}_{}\ \xrightarrow{\ \text{octave }j\text{ filled to }c_E\lambda_j^{a-1}\ }\ \underbrace{\text{wake}}_{}\ \xrightarrow{\ e^{-2\nu N_j^2\tau_j}=e^{-2\nu C^2\tau_j^{a/(4+a)}}\to1\ }\ \underbrace{\Lambda(N)=N^a}_{\text{re-supplied}}\ \circlearrowleft$$

**D.1 Retention (derived, and the decisive advantage of $g<1/2$).** Octave $j$ born at $\tau_j$ decays by $\exp(-2\nu N_j^2\tau_j)=\exp(-2\nu C^2\tau_j^{1-2g})$. $1-2g=a/(4+a)>0$ $\Rightarrow$ **retention $\to1$: the wake is asymptotically loss-free.** ($g=1/2$: constant fraction lost per octave; $g>1/2$: wake annihilated — the exact failing equation of V3.)

**D.2 Closure of $\Lambda$.** $\Lambda=N^a$ is the ansatz's one free structural input. It is *pinned* — not derived — by three inequalities: (a) $\sum_j\Lambda(\lambda_j)/\lambda_j<\infty\Rightarrow a<1$; (b) Type-II escape $\Rightarrow a>0$; (c) capacity/BKM coincidence $\Rightarrow g=2/(4+a)$. Whether the front *self-generates* $\Lambda=N^a$ is the pilot's RG-recurrence test (§F).

**D.3 Front–BKM identity (PROOF CANDIDATE, exact algebra).** On the saturated branch with $\|\omega\|_\infty=\theta\sqrt{2c_E\Lambda}N^2$,
$$\frac{\dot N}{N}=\sqrt2\chi\sqrt{mc_E\Lambda}N^2\ \Rightarrow\ \|\omega\|_\infty=\frac{\theta}{\chi\sqrt m}\frac{\dot N}{N}\ \Rightarrow\ \boxed{\int_0^t\|\omega\|_\infty ds=\frac{\theta}{\chi\sqrt m}\,\log\frac{N(t)}{N(0)}}$$
independent of $\nu,c_E,\Lambda,a,g$. **Consequences.** (1) BKM divergence $\iff N\to\infty$ — nothing more is needed. (2) The BKM integral is *always* only logarithmic in $N$: $I_{\rm BKM}=10$ requires $N/N_0=e^{10\chi\sqrt m/\theta}$. **No numerical run can ever observe BKM divergence**; it must be certified structurally. (3) Design rule: the sup-norm coherence $\theta$ must not be parametrically smaller than the flux coherence $\chi\sqrt m$.

**D.4 Capacity $=$ BKM.** $\int\|\omega\|_\infty dt\propto\int\tau^{-g(2+a/2)}d\tau=\infty\iff g\ge\frac2{4+a}$ — *identically* the capacity floor of §B.5. Saturating capacity is exactly being BKM-marginal.

---

## E. Obstruction audit (exact collision points)

1. **Energy bound / finite dissipation (F-N1/N2, Leray).** Never used as a signature. $E\to2c_E/(1-2^{a-1})<\infty$; $\int_0^T\varepsilon\,dt\propto\int\tau^{-g(1+a)}d\tau<\infty\iff g(1+a)<1\iff a<2$ ✓ (we have $a<1$, margin $\ge0.26$).
2. **ESS $L^\infty_tL^3_x$.** Collision point: needs $\|u\|_3\to\infty$. Here $\|u\|_3^3\asymp\sum_j(\lambda_j^{1/2}\|u_j\|_2)^3=(2c_E)^{3/2}\sum_j\lambda_j^{3a/2}\asymp N^{3a/2}\to\infty$ **polynomially** (vs. only $\log N$ at $a=0$). Requires the reverse-Bernstein/occupation lemma (seed S4) — obligation 5.
3. **Serrin pairs, all of them.** $\|u\|_p\asymp N^{1-3/p+a/2}$, $q=2p/(p-3)$, so $q(1-3/p+a/2)=2+\frac{ap}{p-3}$, and $\int\|u\|_p^qdt=\infty\iff g\big(2+\frac{ap}{p-3}\big)\ge1\iff\frac{2p}{p-3}\ge1$, **true for every $p>3$ with strict margin** (endpoint $p=\infty$: $2g(1+a/2)=1.163>1$ ✓). Every admissible pair diverges — not just the tracked ones.
4. **Fixed-finite-bandwidth no-go (F-6/F-7, F-$\alpha$1, VR-L-011).** $N(t)\to\infty$ is the *state variable*; the ansatz is nowhere confined to a fixed mode set. Collision avoided by construction; and $\gamma=g<1$ satisfies VR-L-019.
5. **CSTY2009 Type-I exclusion (axisym $\mathbb R^3$).** THE collision point of the $a=0$ branch: there $\sqrt\tau\|u\|_\infty\to2\sqrt{c_E/k}$, exactly Type I, so an axisymmetric $\mathbb R^3$ realisation would be excluded. **Evasion, stated precisely:** (i) domain is $\mathbb T^3$, where CSTY (an $\mathbb R^3$ axisymmetric theorem) has no statement; (ii) more importantly, $a>0$ gives $\sqrt\tau\|u\|_\infty\asymp\tau^{-ga/2}=\tau^{-0.081}\to\infty$ — the Type-I bound (2.1) is violated with a definite exponent, and the weighted bound (2.2) fails likewise since $\|u\|_\infty$ exceeds $C_*|t|^{-1/2}$ globally in the front region, not merely locally.
6. **KNSS Liouville / bounded ancient solutions.** Rescaling at the blow-up point yields an ancient limit. KNSS (b) kills swirl-retaining ancient solutions obeying $|u|\le C/r$. Our field is *not* pointwise $C/r$: the rescaled profile has $\|u\|_\infty$ growing like $\tau^{-ga/2}$ relative to the parabolic scale, i.e. the natural rescaling $\lambda=\sqrt\tau$ does **not** produce a bounded limit — the correct rescaling is $\lambda=N^{-1}$, under which the limit object is the front-flow profile $\Psi$ of seed S1, not an ancient NS solution in the KNSS class. Collision point flagged as **open**: nontriviality/boundedness of that limit is unproven (obligation 8).
7. **NRS/Tsai one-scale self-similar no-go.** Our object is **not** one-scale: $E_N=c_EN^{a-1}$ with $a\ne0$ means the profile amplitude is *not* scale-invariant; the front flow is $s$-periodic (discrete doubling), not stationary. Exact collision point: NRS requires $u=(2a(T-t))^{-1/2}U(x/\sqrt{2a(T-t)})$ with $U$ time-independent and $U\in L^3$; here the similarity exponent is $g=0.419\ne1/2$ and the profile is $s$-dependent, so the profile equation of NRS is never satisfied.
8. **Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}_N/N^3$.** Collision point: we require $M^{\rm eff}_N=mN^3$, $m=\Theta(1)$ — i.e. $\gamma=1$ exactly, fixed-relative width $\eta\in(0,1/3)$. Then $D_N\le2\kappa^2\tau^2c_Em$ is $N$-independent, so the no-go is *not* violated, it is *saturated at the only value of $\gamma$ it permits*. Note the amplitude-effective requirement: $M^{\rm eff}\asymp M$ demands Fourier-amplitude delocalisation, which is compatible with (indeed implied by) *physical-space* packet intermittency. No conflict.
9. **Diagonal cross-talk gate.** Unaddressed by this lens and inherited as an obligation: with $m=\Theta(1)$ the front shell is densely occupied, so the $9{:}1$ leakage of the sparse 4-mode alphabet is not directly applicable, but no positive result replaces it. Registered as obligation 4, **open**.
10. **Galerkin global existence.** Any truncation of this ansatz is globally regular; the pilot therefore measures *rates*, never blow-up. Explicitly acknowledged (TM-01/22).
11. **Smooth-forcing high-frequency decay (F-N4).** $f\equiv0$; nothing to evade.
12. **Pure-swirl $L^3$ no-go (VR-L-016).** The ansatz is a generic $\mathbb T^3$ Fourier front, not a pure swirl; the pressure channel $P=3\int p\,\nabla\!\cdot\!(|u|u)$ is not identically zero by symmetry. Not applicable — but also not yet shown positive.
13. **Front-resolution threat model (TM-22, $\ge7$ points per front; TM-10 power-law fitting).** The pilot below fits *no* blow-up time and reports no $T$; per D.3 BKM divergence is unobservable numerically, so the pilot measures only $\chi$ and its $N$-exponent, at $\ge4$ scales with exact-rational signs.

---

## F. Minimal falsification experiment ($\le1$ h)

**Question:** does the achieved coherence $\chi(N)$ decay *slower* than the required $N^{-a/2}$, for some $a<1$, with a **sign-stable positive** flux?

**Design.** Reuse `mesoscopic_cloud_scaling.build_sparse_parent` / `exact_sparse_leray_convolution` / `measure_mesoscopic_cloud`, and `mesoscopic_local_fft.measure_local_fft_cloud` (zero-padded, $K=4W-3$, no wrap aliasing). Fixed-relative width $\eta=0.2$ ($\gamma=1$), $N\in\{16,32,48,64\}$, $a\in\{0,0.25,0.5,0.75\}$, $\nu=1/40$, $\tau=1/4$. For each $(N,a)$: set shell amplitude $\|u_N\|_2^2=2c_EN^{a-1}$, compute the **signed** transfer $\Pi_N=\langle P_CB(u_N,u_N),\,\hat c\rangle$ into the doubled band, and
$$\chi_{\rm meas}=\frac{\Pi_N}{2\sqrt2\sqrt m\,c_E^{3/2}N^{3a/2}N},\qquad \mathrm{VM}:=\frac{\sqrt2\,\chi_{\rm meas}\sqrt{mc_E}\,N^{a/2}}{2\nu}\ \ (\text{viscous margin}).$$

**Variables:** $N$, $a$, $\eta$, phase code. **Success:** $\chi_{\rm meas}>0$ at all four $N$ for at least one $a<1$, with fitted exponent $\beta_\chi>-a/2-1/4$ **and** $\mathrm{VM}$ increasing in $N$ at $\ge3$ of 4 scales. **Kill:** sign of $\Pi_N$ indefinite at any $N\ge32$ (this is the repo's standing "positive flux margin" gap and would kill the lens outright), or $\beta_\chi<-a/2-1/4$ across all $a$.

**Arithmetic:** the **sign** of $\Pi_N$ and of the doubling-pullback overlap must be **exact rational** — route through `exact_carrier_search` / independent `exact_carrier_record_verifier`. Magnitudes, exponent regressions, and $\mathrm{VM}$ may be binary64. Resolution: local-FFT backend at $N\le64$ suffices; $\ge7$ points across the front width $W=\eta N=12.8$ ✓ (TM-22).

**Second, free check (5 min, symbolic):** verify $\int\|\omega\|_\infty dt=\frac{\theta}{\chi\sqrt m}\log(N/N_0)$ (D.3) as finite algebra; it is Lean-able in the style of `MesoscopicDuhamelNoGo.lean`.

**Reported ceiling:** "front-coherence exponent measured at stated $\eta,\nu,\tau,N$" — never a blow-up claim.

---

## G. Proof chain (10 obligations)

1. **(FL)** — front energy budget as a rigorous inequality for a band-limited $\mathbb T^3$ field (finite Fourier algebra; Lean-able).
2. **(CAP)** — $|\Pi_N|\le2\sqrt2\sqrt{m}c_E^{3/2}\Lambda^{3/2}N$ with explicit $m$; already essentially `MesoscopicDuhamelNoGo`.
3. **Positive flux margin:** $\exists\chi_0>0$, a phase code with $\Pi_N\ge\chi_0\cdot(\text{CAP})$ for all $N$ in the chain. *The single decisive open step of the whole repo.*
4. **Off-chain leakage summability:** cross-talk + spill energy per octave $\le\epsilon_j$, $\sum\epsilon_j<\tfrac12$ (diagonal cross-talk gate at $\gamma=1$).
5. **Occupation / reverse Bernstein:** $\|u_j\|_3\ge c\lambda_j^{1/2}\|u_j\|_2$ on the chain $\Rightarrow\|u\|_3^3\gtrsim N^{3a/2}$ (seed S4).
6. **Wake retention:** $\prod_j e^{-2\nu N_j^2\tau_j}\ge c>0$ for $g<1/2$ (D.1).
7. **Front ODE $\Rightarrow$ $N\to\infty$ at finite $T$** with $T-t_0=\int N^{-(4+a)/2}\!\cdot$const$<\infty$, $g=2/(4+a)$.
8. **RG recurrence:** the front flow of seed S1 has an attracting $s$-periodic profile reproducing $\Lambda=N^a$ (forward RG integration, not optimisation).
9. **Entry:** a smooth band-limited $u_0$ enters the basin of that orbit (PO-09, strategy-less; the acknowledged hardest analytic step).
10. **BKM $\Rightarrow$ non-extendability:** $\int_0^T\|\omega\|_\infty=\frac{\theta}{\chi\sqrt m}\log N\to\infty$ (D.3) $\Rightarrow$ $T$ maximal $\Rightarrow$ **Clay (D)** (and with $f\equiv0$, refutation of (B)).

---

## Rejected sub-variants (kept, with the failing equation)

| id | variant | exact failing equation |
|---|---|---|
| **V1** | random-phase front, any wake | $\dot N=\tfrac1{\sqrt2}\sqrt{c_E}N^{3/2+a/2}-2\nu N^3<0$ for $N>(\sqrt{c_E}/2\sqrt2\nu)^{2/(3-a)}$ |
| **V2** | $g>1/2$, capacity-saturated (needs $a<0$) | $\chi_{\rm req}=\frac{\sqrt2\nu}{\sqrt{mc_E}}N^{|a|/2}\to\infty>1$ |
| **V3** | $g>1/2$, tuned sub-capacity, $\Lambda\equiv1$ | requires $\chi\sqrt{2mc_E}-2\nu=2gC^{-1/g}N^{1/g-2}\downarrow0$ exactly; and wake retention $e^{-2\nu C^2\tau^{1-2g}}\to0$ |
| **V4** | supercritical wake $a\ge1$ | $E=\sum_j c_E2^{j(a-1)}=\infty$ |
| **V5** | $g<1/2$ with exactly-critical wake $a=0$ | $\chi_{\rm req}=\frac{2\nu+gC^{-1/g}N^{1/g-2}}{\sqrt{2mc_E}}\to\infty$ since $1/g-2>0$ |
| **V6 (demoted, not rejected)** | $g=1/2$, $a=0$ (current repo target) | survives on $\mathbb T^3$ but has a permanent threshold $\chi\ge\sqrt2\nu/\sqrt{mc_E}$ and is pointwise Type I ($\sqrt\tau\|u\|_\infty=2\sqrt{c_E/k}$), colliding with CSTY2009 on axisym $\mathbb R^3$ |
