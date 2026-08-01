# LENS 5 — Lagrangian deformation-gradient degeneracy: the $\Theta$-selection front

**Working name:** CLFA — *Chirally-biased Lagrangian Fold-Annihilation front.*
**Status: FORMAL ANSATZ**, containing two **PROOF CANDIDATE** lemmas (L1, L3) and three **REJECTED** sub-variants (V1, V2, V4) kept below with their exact failing equations.

---

## A. Clay target

- **Target: TARGET-U** (unforced $\mathbb R^3$ breakdown; implies **(C)**). Not (A)/(B).
- **Domain:** $\mathbb R^3$, whole space. (Repo rule: no import of Track-P periodic closures.)
- **Forcing:** none. $f\equiv0$, so F-N1/F-N2/F-N4 are structurally inapplicable.
- **Initial data:** $u_0\in C_c^\infty(\mathbb R^3;\mathbb R^3)$, $\nabla\!\cdot u_0=0$, **deliberately non-axisymmetric** and with non-zero, single-signed local helicity density in the seed region. Non-axisymmetry is a *load-bearing design constraint* (see E: CSTY, KNSS).
- **Viscosity:** $\nu>0$ fixed and treated **non-perturbatively**. The whole point of §B/E is that at the selected front $\nu$ enters at order one; the "Euler first, viscosity as a correction" strategy is proved inadmissible here (V3).
- **State variable:** the deformation gradient $F(a,t)=\partial X/\partial a$ of the *deterministic* NS flow map $\dot X=u(X,t)$, $\det F\equiv1$.

---

## B. Central mathematics (derived)

### B.1 Exact singular-value ODE

$F$ satisfies $\dot F=AF$, $A_{ij}=\partial_j u_i$, **for NS exactly** (only the Cauchy formula $\omega=F\omega_0$ needs $\nu=0$; the flow map does not). Write the SVD $F=R\Sigma Q^{T}$, $\Sigma=\operatorname{diag}(\sigma_1\ge\sigma_2\ge\sigma_3>0)$, $\sigma_1\sigma_2\sigma_3=1$. With $\Omega_R=R^T\dot R$, $\Omega_Q=Q^T\dot Q$ (both antisymmetric) and $\hat A=R^TAR$, the relation $\dot F=AF$ becomes $\Omega_R\Sigma+\dot\Sigma-\Sigma\Omega_Q=\hat A\Sigma$. Its diagonal ($\Omega$'s have zero diagonal) gives

$$\boxed{\ \dot\sigma_i=\sigma_i\,(r_i\cdot S\,r_i)\ },\qquad S=\tfrac12(A+A^T),\ r_i=Re_i. \tag{B1}$$

(The antisymmetric part drops since $r\cdot Wr=0$.) Check: $\sum_i r_i\!\cdot\!Sr_i=\operatorname{tr}S=0\Leftrightarrow\frac{d}{dt}\det F=0$. ✓ Dimensions: $[S]=T^{-1}$. ✓

Define the **compressive rate** $\lambda_c:=-r_3\cdot S r_3$ and $n:=r_3$. Then

$$\sigma_3(a,t)=\exp\Big(-\!\int_0^t\!\lambda_c\,ds\Big),\qquad \delta(t):=\epsilon\,\sigma_3(t),\qquad N(t):=1/\delta(t). \tag{B2}$$

**Geometric meaning of $n$.** The off-diagonal of the same SVD relation gives $(\Omega_R)_{ij}=(\hat A_{ij}\sigma_j^2+\hat A_{ji}\sigma_i^2)/(\sigma_j^2-\sigma_i^2)$, hence for $\sigma_3\ll\sigma_2$

$$\dot n=-(I-n\otimes n)A^{T}n+O(\sigma_3^2/\sigma_2^2). \tag{B3}$$

This is *exactly* the direction equation of a materially transported gradient ($\frac{d}{dt}\nabla\theta=-A^T\nabla\theta$). Indeed $\nabla_x(q_3\!\cdot\!a)=\sigma_3^{-1}n$: **$n$ is the normal of the folded material sheet and $1/\sigma_3$ is the label-gradient amplification.** This is the analyticity-strip link: the Fourier envelope of a field whose level sets are separated by $\delta$ is $e^{-\delta_a|k|}$ with $\delta_a\asymp\delta=\epsilon\sigma_3$, so $\delta_a(t)\asymp1/N(t)$ **is** $\sigma_3\to0$.

### B.2 The $\lambda_c$ Riccati and the shape number $\Theta$

Along a trajectory, $\dot A=-A^2-H+\nu\Delta A$ with $H_{ij}=\partial_i\partial_jp$, $\operatorname{tr}H=-\operatorname{tr}A^2$. Symmetrising with $A=S+W$, $W_{ij}=-\tfrac12\epsilon_{ijk}\omega_k$, $(W^2)_{ij}=\tfrac14(\omega_i\omega_j-|\omega|^2\delta_{ij})$:

$$\dot S=-S^2+\tfrac14\big(|\omega|^2I-\omega\otimes\omega\big)-H+\nu\Delta S. \tag{B4}$$

(Trace check: $-|S|^2+\tfrac12|\omega|^2+\operatorname{tr}A^2=0$ since $\operatorname{tr}A^2=|S|^2-\tfrac12|\omega|^2$. ✓)

Now $\dot\lambda_c=-\big[2\dot n\cdot Sn+n\cdot\dot Sn\big]$. Using (B3), $P^\perp Sn=Sn+\lambda_c n$ and $n\cdot Wn=0$:
$2\dot n\cdot Sn=-2|Sn|^2+2\lambda_c^2+2(Wn)\!\cdot\!(Sn)$. Substituting (B4):

$$\boxed{\ \dot\lambda_c=3|Sn|^2-2\lambda_c^2-2(Wn)\!\cdot\!(Sn)-\tfrac14\big(|\omega|^2-(\omega\!\cdot\!n)^2\big)+n\!\cdot\!Hn-\nu\,n\!\cdot\!(\Delta S)n\ } \tag{B5}$$

**Exact validation.** Pure linear strain $u=Sx$, $S=\operatorname{diag}(\alpha,0,-\alpha)$, $\omega=0$, $H=-S^2$: RHS $=3\alpha^2-2\alpha^2+0-0-\alpha^2-0=0=\dot\lambda_c$. ✓ Exact, non-trivial cancellation.

Normalise by $\lambda_c^2$ and define the **scale-free shape number**

$$\Theta:=\underbrace{3\rho-2}_{\text{misalignment}}+\underbrace{\chi}_{\text{tilt}}-\underbrace{\Xi}_{\text{rotation brake}}+\underbrace{\Pi}_{\text{pressure}}-\underbrace{\mathcal V}_{\text{viscous}},$$
$$\rho=\frac{|Sn|^2}{\lambda_c^2}\ge1,\quad \chi=\frac{-2(Wn)\!\cdot\!(Sn)}{\lambda_c^2},\quad \Xi=\frac{|\omega|^2-(\omega\!\cdot\!n)^2}{4\lambda_c^2}\ge0,\quad \Pi=\frac{n\!\cdot\!Hn}{\lambda_c^2},\quad \mathcal V=\frac{\nu\,n\!\cdot\!(\Delta S)n}{\lambda_c^2}. \tag{B6}$$

Then $\dot\lambda_c=\Theta\lambda_c^2$. Under NS scaling $u\mapsto\lambda u(\lambda x,\lambda^2t)$ every one of $\rho,\chi,\Xi,\Pi,\mathcal V$ is invariant ($S\mapsto\lambda^2S$, $H\mapsto\lambda^4H$, $\Delta S\mapsto\lambda^4\Delta S$). **$\Theta$ is a pure shape functional** — the exact Lagrangian analogue of the seed's $\chi_{\rm shape}$/$c_E$-collapse (S3): the whole lane reduces to one dimensionless number.

If $\Theta\to\Theta_\star$ then $\lambda_c\sim1/(\Theta_\star\tau)$, $\tau=T-t$, and by (B2)

$$\sigma_3\asymp\tau^{\gamma},\qquad N\asymp\tau^{-\gamma},\qquad \boxed{\gamma=1/\Theta_\star}. \tag{B7}$$

Sign/consistency: $|Sn|^2\ge(n\!\cdot\!Sn)^2=\lambda_c^2$, so $3\rho-2\ge1$; the aligned inviscid baseline is $\Theta=1$, i.e. $\gamma=1$. **Misalignment of the sheet normal with the strain eigenframe is exactly the quantity that pushes $\gamma$ below 1.** Note also $\Pi_{\rm linear\ strain}=-1$: a globally linear strain field is Riccati-neutral; a *positive* normal pressure curvature $\partial_n^2p>0$ (a pressure minimum in the sheet plane — the vortex-core signature) is required to drive $\Theta$ up.

### B.3 The monotone Lagrangian quantity

From $\dot\sigma_3=-\lambda_c\sigma_3$ and $\dot\lambda_c=\Theta\lambda_c^2$:

$$\mathrm{Re}_L:=\frac{\sigma_3^2\lambda_c}{\nu}\quad\Longrightarrow\quad \boxed{\ \frac{d}{dt}\log \mathrm{Re}_L=\lambda_c\,(\Theta-2)\ }. \tag{B8}$$

$\mathrm{Re}_L$ is the **local cell Reynolds number of the compressive front** ($[\sigma_3^2\lambda_c]=L^2T^{-1}$). Consequences, all exact:
- $\Theta\le2$ everywhere $\Rightarrow$ $\mathrm{Re}_L$ non-increasing (a genuine monotone Lagrangian functional).
- $\mathrm{Re}_L\to0$ $\Rightarrow$ the front is viscously slaved; no collapse can be sustained.
- A collapse with $\int^T\!\lambda_c=\infty$ and $\mathrm{Re}_L$ bounded above and below forces $\langle\Theta\rangle=2$ in the $\lambda_c\,dt$-average. **$\Theta_\star=2$, i.e. $\gamma=1/2$, is not chosen — it is *selected* by (B8).**

### B.4 Geometry pinned by the ordering + BKM

$\sigma_1\sigma_2=1/\sigma_3=\tau^{-\gamma}$; write $\sigma_2=\tau^{m}$, $\sigma_1=\tau^{-\gamma-m}$. Two constraints:
1. **$\Xi=O(1)$**: $\Xi\asymp|\omega|^2/(4\lambda_c^2)$ and $|\omega|\asymp\sigma_1=\tau^{-\gamma-m}$ (stretching), $\lambda_c\asymp\gamma/\tau$. $\Xi$ bounded $\Leftrightarrow\gamma+m=1$. If $\gamma+m>1$ then $\Xi\to\infty$ and $\Theta\to-\infty$: collapse self-quenches. So $m=1-\gamma$: $\ \sigma=(\tau^{-1},\ \tau^{1-\gamma},\ \tau^{\gamma})$.
2. **Ordering $\sigma_2\ge\sigma_3$**: $\tau^{1-\gamma}\ge\tau^{\gamma}\Leftrightarrow\gamma\ge1/2$.

Combined with $\gamma=1/\Theta_\star$: $\Theta_\star\in(1,2]$, i.e. $\gamma\in[1/2,1)$ — **exactly the repo's surviving Zeno window**, derived here independently from Lagrangian algebra.

Circulation check: $|\omega|\,\sigma_2\sigma_3=\tau^{-1}\cdot\tau^{1-\gamma}\cdot\tau^{\gamma}=O(1)$ — Kelvin-consistent. ✓

### B.5 Where viscosity enters at order one (the decisive computation)

$n\cdot\Delta S\,n\asymp N^2\lambda_c$ in magnitude, with sign $+$ (the compressive extremum is a minimum of $n\!\cdot\!Sn$), so

$$\mathcal V\asymp\frac{\nu N^2}{\lambda_c}=\frac{\nu}{\gamma}\,\tau^{\,1-2\gamma}. \tag{B9}$$

- $\gamma>1/2$: $\mathcal V\to\infty$, $\Theta\to-\infty$. **The entire interval $\gamma\in(1/2,1)$ is viscously destroyed** (sub-variant V4, REJECTED below).
- $\gamma=1/2$: $\mathcal V\to2\nu\,c_{\rm shape}$, a *finite order-one constant*. This is the honest answer to the lens question: viscosity is neither negligible nor dominant at $\gamma=1/2$; it contributes a fixed subtraction to $\Theta$.
- $\gamma<1/2$: excluded by B.4.

**Selection theorem (heuristic level):** $\gamma=1/2$ uniquely, and the closure condition is

$$\boxed{\ \big\langle 3\rho-2+\chi-\Xi+\Pi\big\rangle=2+\mathcal V_\star,\qquad \mathcal V_\star=\nu\lim\frac{n\cdot\Delta Sn}{\lambda_c^2}>0.\ } \tag{B10}$$

### B.6 The mechanism proper: fold–annihilation

At $\gamma=1/2$: $\sigma_1/\sigma_3=\tau^{-3/2}$ while the front's *outer extent* is $\ell\asymp\tau^{1/2}$ (see C). A material line therefore folds $\asymp\tau^{-3/2}/\tau^{-1/2}\!\cdot\!$-many times inside the front, producing structure at scale $\tau$, i.e. bandwidth $\tau^{-1}\gg N$. If those folds persisted, $\gamma=1$ (excluded). They do not: adjacent folds are **antiparallel** and are annihilated by viscosity at the scale $\sqrt{\nu\tau}\asymp\sqrt\nu\,\delta$ — the *same* scale as the front. The surviving vorticity is the **chirality-biased residue**: exactly antiparallel folds cancel completely, so the net $\omega$ is proportional to the helicity/chirality imbalance of the fold packing. This is the precise sense in which the Cauchy formula fails at order one (B.5, V3), and it is why the seed's helicity ledger (S5) becomes a *design rule*, not a diagnostic: zero-helicity fold packing $\Rightarrow$ complete annihilation $\Rightarrow$ $\Theta$ collapses.

---

## C. Scaling table ($\tau=T-t$, $\gamma=1/2$, critical wake $E_N=c_E/N$)

| quantity | law in $\tau$ | source |
|---|---|---|
| $\sigma_1,\sigma_2,\sigma_3$ | $\tau^{-1},\ \tau^{1/2},\ \tau^{1/2}$ | B.4 |
| $\lambda_c$ | $\tfrac12\tau^{-1}$ | (B7) |
| Fourier bandwidth $N=1/\delta=1/(\epsilon\sigma_3)$ | $\tau^{-1/2}$ | (B2) |
| analyticity strip $\delta_a$ | $\tau^{1/2}$ (floor $\sqrt{\nu\tau}$) | B.1 |
| $\|\omega\|_\infty=\|\nabla u\|_\infty$ | $\tau^{-1}$ | B.4 |
| $\|u\|_\infty$ | $\tau^{-1/2}$ (**Type-I marginal**) | $\|\nabla u\|_\infty/N$ |
| front extent $\ell_1\!=\!\ell_2\!=\!\ell_3$ | $\tau^{1/2}$; active volume $2c_E\tau^{3/2}$ | $U^2V_N=2c_E/N$ |
| front energy $E_N$ | $c_E\tau^{1/2}\to0$ | critical wake |
| front enstrophy $N^2E_N$ | $c_E\tau^{-1/2}$ | — |
| dissipation rate $2\nu N^2E_N$ | $2\nu c_E\tau^{-1/2}$; $\int_0^T\!<\infty$ ✓ | integrable iff $\gamma<1$ |
| $\|\mathbb P(u\!\cdot\!\nabla u)\|_{L^2}$ | $\tau^{-3/4}$ | $U\!\cdot\!NU\!\cdot\!V_N^{1/2}$ |
| $\|\nabla p\|_{L^2}$ | $\tau^{-3/4}$ (same order; $\mathbb P$ bounded) | — |
| $\nu\|\Delta u\|_{L^2}$ | $\nu\sqrt{2c_E}\,\tau^{-3/4}$ | **all three balance** ✓ |
| global $\|u\|_{L^3}^3$ | $\asymp$ #octaves $\asymp\tfrac12\log(1/\tau)\to\infty$ | seed S4; per-octave $U^3V_N=O(1)$ ✓ |
| $\int^t\|\omega\|_\infty$ | $\log(1/\tau)\to\infty$ (BKM ✓, log-slow) | — |
| physical time remaining | $\tau=e^{-s}$, $s=\log(1/\tau)$; $\int e^{-s}ds<\infty$ ✓ | — |
| active mode count $M_N^{\rm eff}$ | $\asymp\tfrac43\pi N^3$ (front size $=1/N$ $\Rightarrow$ maximal Fourier delocalisation) | E-6 below |
| $\mathrm{Re}_L$ | constant (marginal) | (B8) |

The triple balance at $\tau^{-3/4}$ of nonlinear / pressure / viscous $L^2$ norms is an independent confirmation that $\gamma=1/2$ is the distinguished exponent.

---

## D. Closed feedback loop (every arrow a formula)

$$\lambda_c \xrightarrow{\ \sigma_3=e^{-\int\lambda_c}\ }\ \delta=\epsilon\sigma_3\ \xrightarrow{\ N=1/\delta\ }\ \text{fold multiplicity } \mathcal F=\sigma_1/(N^{-1}\sigma_3^{-1})$$
$$\xrightarrow[\text{annihilation}]{\ \text{viscous, scale }\sqrt{\nu\tau}\ }\ \omega_{\rm net}=h\,\Gamma_0/(\sigma_2\sigma_3)\ \xrightarrow[\text{Biot--Savart}]{\ \Pi=n\cdot Hn/\lambda_c^2,\ \rho=|Sn|^2/\lambda_c^2\ }\ \Theta=3\rho-2+\chi-\Xi+\Pi-\mathcal V$$
$$\xrightarrow{\ \dot\lambda_c=\Theta\lambda_c^2\ }\ \lambda_c .$$

Gain condition of the loop: $\Theta=2$ (B10). Regulator: (B8), $\frac{d}{dt}\log\mathrm{Re}_L=\lambda_c(\Theta-2)$ — the loop is **self-stabilising at $\Theta=2$**: $\Theta>2$ raises $\mathrm{Re}_L$, thickening the front relative to $\sqrt{\nu\tau}$, which *reduces* the annihilation efficiency and lowers $h$ hence $\Pi$; $\Theta<2$ does the reverse. $h\in(0,1]$ is the chirality bias of the fold packing (S5); $h=0$ (achiral, exactly antiparallel folds) $\Rightarrow\omega_{\rm net}=0\Rightarrow\Pi\to-1\Rightarrow\Theta<0$: **the loop is dead without chirality.** This is the single new design rule this lens contributes.

---

## E. Obstruction audit (exact collision points)

**E-1. Energy bound / finite dissipation (Leray; F-N1/N2; VR).** No collision: $E_N=c_E\tau^{1/2}\to0$, total dissipation $\propto\int_0^T\!\tau^{-1/2}d\tau<\infty$. The blow-up signature is *never* energy or enstrophy; it is $\|u\|_{L^3}\to\infty$ (log) and $\int\|\omega\|_\infty=\infty$ (log). Collision point avoided *because* $\gamma=1/2<1$; $\gamma\ge1$ would make $\int\nu N^2E_N\,dt$ divergent.

**E-2. ESS $L^\infty_tL^3_x$ (S3).** Must be violated, and is, only logarithmically: $\|u\|_3^3\asymp\tfrac12\log(1/\tau)$, via seed S4's per-octave lemma with the exact per-octave value $U^3V_N=O(1)$ computed in C. This is the weakest quantitative margin in the whole scheme and is the first thing an interval computation must certify.

**E-3. Fixed-finite-bandwidth no-go (F-6/F-7, F-α1, VR-L-011).** No collision: $N(t)=\tau^{-1/2}\to\infty$; the ansatz is *defined* by unbounded bandwidth. The Galerkin theorem's hypothesis (fixed mode set $V_S$) fails at every $t$.

**E-4. One-scale self-similar no-go (NRS1996 $L^3$; Tsai1998 local energy).** **Live collision.** $\gamma=1/2$ *is* the self-similar exponent, and $\sigma_2=\sigma_3=\tau^{1/2}$ makes the core isotropic. Evasion must be exact, not vague: the object here is **not** a time-independent profile $U(y)$. The Lagrangian data $(\sigma_1,\sigma_2,\sigma_3)=(\tau^{-1},\tau^{1/2},\tau^{1/2})$ has $\sigma_1/\sigma_3=\tau^{-3/2}$, which is *not* a similarity invariant — a genuine backward-self-similar solution has $\nabla u=\tau^{-1}G(y)$ with $G$ frozen, forcing $\Theta\equiv2$ pointwise and $\mathcal F\equiv$ const. Here $\mathcal F\to\infty$, so the profile in similarity variables $y=x/\sqrt{2a\tau}$, $s=\log(1/\tau)$ **must be $s$-dependent**. The concrete claim is that it is $s$-*periodic* (discrete self-similar / log-periodic), which is exactly the recurrent object of the orchestrator's front flow S1 ($S=\log2/a_+$). NRS/Tsai hypotheses require $\partial_sU=0$; they are not satisfied. **Residual threat: Seregin2024 Liouville-type restrictions on DSS solutions** — not yet converted to a filter in-repo; this is the single largest unresolved obstruction for this idea and must be checked before any promotion beyond FORMAL ANSATZ.

**E-5. CSTY2009 axisymmetric Type-I exclusion.** **Live collision.** Table C gives $\sqrt\tau\|u\|_\infty\to$ const, i.e. exactly Type-I marginal, and B.4 *proves* Type-II is unreachable in this lens: Type-II would need $\|\nabla u\|_\infty\asymp\tau^{-b}$, $b>\gamma+1/2$, while $\lambda_c\asymp\gamma/\tau$; then $\Xi\asymp\tau^{2-2b}\to\infty$ and $\Theta\to-\infty$ unless $\omega\parallel n$, but $\omega\parallel r_3$ means vorticity along the compressed axis, which decays. So $b=1$ is forced. **Therefore the candidate cannot evade CSTY by being Type-II; it evades only by failing CSTY's hypothesis of axisymmetry.** Design constraint §A (non-axisymmetric, chiral) is load-bearing. Any axisymmetric realisation of this mechanism is *dead on arrival*.

**E-6. Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM_N^{\rm eff}/N^3$ (L-11/L-11a).** **Direct collision, and this is the idea's strongest point.** The required relay ratio is $D_N=E_{2N}/E_N=1/2$ at critical wake. The bound needs $M_N^{\rm eff}\gtrsim N^3$. Here the front's *physical* extent is $\ell\asymp\tau^{1/2}$ and its bandwidth is $N\asymp\tau^{-1/2}$, so $\ell=1/N$ exactly: the front is a bump whose width equals its own inverse bandwidth, hence its Fourier coefficients fill the whole ball $|k|\le\kappa N$ with comparable magnitudes and $M_N^{\rm eff}\asymp\frac43\pi N^3$ — **the L-11 floor is saturated, not violated.** In the cloud lane's own language this is the fixed-relative-width regime $W_N=\eta N$ (the *only* structurally uncondemned regime per the constraint map), reached here not by construction but as a forced consequence of $\gamma=1/2$ + critical energy. The residual obligation is the same one the cloud lane left open: a *scale-independent positive floor* for $D_N$, which in this lens is literally (B10), $\langle\Theta\rangle=2$.

**E-7. Diagonal cross-talk gate.** Not directly applicable — there is no discrete carrier alphabet. Its analogue is off-front leakage: energy deposited outside the $\sqrt\tau$-ball. The regulator is (B8); leakage manifests as $\Theta<2$ and monotone $\mathrm{Re}_L$ decay. The pilot (F) measures this directly instead of a sumset-miss condition.

**E-8. KNSS2009 bounded ancient Liouville.** At Type-I marginality the rescaled limit is a nonzero bounded ancient mild solution. KNSS(a) needs axisymmetry+no swirl; KNSS(b) needs axisymmetry+swirl+$|u|\le C/r$. Neither hypothesis holds: the candidate is non-axisymmetric by construction. The general (non-symmetric) bounded-ancient Liouville problem is open, so this obstruction does not bite — but equally it provides no support.

**E-9. Pure-swirl $L^3$ no-go (VR-L-016).** Not applicable: the mechanism's driver is $\Pi=n\!\cdot\!Hn/\lambda_c^2>0$, an explicitly *non-zero* pressure channel, whereas the pure-swirl kill is precisely $u\!\cdot\!\nabla p_0\equiv0$. Non-axisymmetry guarantees $\partial_\theta p\not\equiv0$.

**E-10. Smooth-forcing high-frequency decay (F-N4).** N/A (unforced).

**E-11. BKM (Kato–Ponce).** **This is the lens's main positive lemma, not an obstruction.** See L1 below.

**E-12. Front-resolution threat model (TM-01/04/20/22).** $\delta=\tau^{1/2}$ with viscous floor $\sqrt{\nu\tau}$; PREREGISTERED_MIN_POINTS_PER_FRONT $=7$ caps the admissible $\tau$ for any given grid. Tracer SVD is exposed to TM-09 (catastrophic cancellation) once $\sigma_1/\sigma_3\gtrsim10^{8}$ in binary64 — the pilot must log $\operatorname{cond}(F)$ and stop before that.

### Named lemmas

**L1 (PROOF CANDIDATE, viscosity-free).** Let $u$ be a strong solution on $[0,T)$, $u\in C([0,T);H^m)$, $m>5/2$, so the flow map is $C^1$. If for some label $a$ and some $T<\infty$, $\liminf_{t\uparrow T}\sigma_3(a,t)=0$ (equivalently $\int_0^T\lambda_c\,dt=+\infty$), then $\int_0^T\|\nabla u\|_{L^\infty}dt=\infty$, and hence by BKM–Kato–Ponce together with the Kozono–Taniuchi logarithmic inequality $\|\nabla u\|_\infty\le C(1+\|\omega\|_{BMO}\log(e+\|u\|_{H^m}))$, also $\int_0^T\|\omega\|_{L^\infty}dt=\infty$ — i.e. $T$ is a genuine breakdown time. *Proof sketch:* $\lambda_c\le\|S(t)\|_{L^\infty}\le\|\nabla u\|_{L^\infty}$, so $\int\lambda_c=\infty\Rightarrow\int\|\nabla u\|_\infty=\infty$; if $\int\|\omega\|_\infty<\infty$ then BKM extends the solution, bounding $\|u\|_{H^m}$ on $[0,T]$ and making $\|\nabla u\|_\infty$ bounded — contradiction. **Crucially this uses only $\dot F=AF$, which is exact for NS.** Lean-able: it is finite-dimensional ODE algebra plus two cited inequalities.

**L3 (PROOF CANDIDATE, negative).** *The Cauchy formula cannot be treated perturbatively at any collapsing front.* Exactly, $\frac{d}{dt}\big(F^{-1}\omega\big)=\nu F^{-1}\Delta\omega$ (from $\frac{d}{dt}F^{-1}=-F^{-1}A$ and $\dot\omega=A\omega+\nu\Delta\omega$). Hence $|\dot\zeta|\le\nu\sigma_3^{-1}\|\Delta\omega\|_\infty$, and with $\sigma_3=\tau^\gamma$, $\|\Delta\omega\|_\infty\asymp N^2\|\omega\|_\infty=\tau^{-2\gamma-1}$,
$$\int^T\!\!\nu\,\sigma_3^{-1}\|\Delta\omega\|_\infty\,dt\ \asymp\ \nu\!\int^T\!\!\tau^{-3\gamma-1}d\tau\ =\ \infty\quad\text{for every }\gamma>0 .$$
So the Euler-first strategy is inadmissible; the correct representation is Constantin–Iyer stochastic ($u=\mathbb E\,\mathbb P[(\nabla A)^T(u_0\!\circ\!A)]$), in which $F^W$ still obeys $\dot F^W=A(X^W)F^W$ with $\det F^W=1$, so **(B1), (B5), (B8) and L1 survive verbatim path-by-path.**

---

## REJECTED sub-variants (kept, with exact failing equation)

- **V1 — Restricted-Euler / Vieillefosse closure** ($H=-\tfrac13\operatorname{tr}(A^2)I$). Fails: this is the aligned, isotropic-pressure, $\Xi$-free case, i.e. $\rho=1,\chi=\Xi=0,\Pi=1-3\rho+2=\ldots$ collapsing to $\Theta_{\rm RE}=1$, hence **$\gamma=1/\Theta=1$**, exactly the boundary excluded by VR-L-019/F-α (and by $\int\nu N^2E_N\,dt=\infty$). *Failing equation:* $\Theta_{\rm RE}=1\Rightarrow\gamma=1$.
- **V2 — Pancake branch** $\sigma_1=\sigma_2=\sigma_3^{-1/2}$. Then $\|\omega\|_\infty\le\sigma_1\|\omega_0\|_\infty=\tau^{-\gamma/2}$ and BKM needs $\gamma/2\ge1$, i.e. $\gamma\ge2$, contradicting $\gamma<1$. *Failing equation:* $\gamma/2\ge1$ vs $\gamma<1$.
- **V3 — Perturbative viscous Cauchy formula.** REJECTED by L3: $\nu\int\tau^{-3\gamma-1}d\tau=\infty$ for all $\gamma>0$.
- **V4 — Type-II branch $\gamma\in(1/2,1)$.** REJECTED by (B9): $\mathcal V\asymp(\nu/\gamma)\tau^{1-2\gamma}\to\infty\Rightarrow\Theta\to-\infty\Rightarrow\dot\lambda_c<0$. *Failing equation:* $\mathcal V=\nu N^2/\lambda_c\to\infty$ for $\gamma>1/2$.

---

## F. Minimal falsification experiment ($\le$1 h)

**Claim under test:** there exist Lagrangian trajectories in a non-axisymmetric, chirality-biased Leray flow along which $\Theta$ (B6) is sustained in $(1,2]$ with $\mathrm{Re}_L$ (B8) not decaying.

**Design.** Grid $128^3$ (fallback $64^3$), dealiased, $\nu$ from the repo's mesoscopic configs, $t\in[0,\tau_p]$ with $\tau_p=\tau N^{-2}$-style parabolic units; RK4$\times$16 as in `mesoscopic_galerkin`. Initial data: (i) `leray_response_relay.helical_fejer_packet` (chiral, $h\ne0$), (ii) `mesoscopic_galerkin.build_angle_box_parent` (achiral control, expected $h\approx0$), (iii) the repo's Gaussian vortex $W_a=(a\times x)e^{-|x|^2/2}$ superposed antiparallel pair.

**Per-step diagnostics.** Compute $A=\nabla u$ spectrally; $S,W$; pressure Hessian exactly in Fourier, $\hat H_{ij}(k)=\dfrac{k_ik_jk_lk_m}{|k|^2}\widehat{u_lu_m}(k)$ (from $\hat p=-\frac{k_lk_m}{|k|^2}\widehat{u_lu_m}$); $\Delta S$ spectrally. Advect $10^4$ tracers, integrate $\dot F=AF$ with the same RK4, SVD each $F$, extract $n=r_3$, $\lambda_c$, then $\rho,\chi,\Xi,\Pi,\mathcal V,\Theta$ and $\mathrm{Re}_L$.

**Variables logged:** $\sigma_{1,2,3}$, $\operatorname{cond}F$, $\lambda_c$, the five $\Theta$-components, $\Theta$, $\mathrm{Re}_L$, local helicity $u\!\cdot\!\omega$, $\|\omega\|_\infty$, points-per-front.

**Success criterion (pre-registered).** A set of tracers of non-vanishing measure (>1% of tracers, contiguous in label space) with $\lambda_c$-weighted mean $\langle\Theta\rangle\in(1,2]$ over the last half of the run **and** $\mathrm{Re}_L$ non-decaying, **and** $\langle\Theta\rangle$ measurably larger for the chiral datum than for the achiral control (the loop's own design rule, D).

**Kill conditions (pre-registered).** (a) $\langle\Theta\rangle<1$ for every tracer in every datum ($\Rightarrow\gamma>1$, excluded — mechanism dead). (b) $\mathrm{Re}_L$ decays exponentially for all tracers. (c) $\Xi$ grows without bound on the top-$\lambda_c$ tracers (rotation brake wins). (d) chiral and achiral data are statistically indistinguishable in $\langle\Theta\rangle$ (kills the fold-annihilation loop D specifically, even if $\Theta$ is large).

**Arithmetic.** *Float is sufficient* for $\Theta$ statistics. **Exact rational arithmetic is required** for three identity checks, run on `exact_leray_relay.build_exact_relay_triad` (exact rational trig field): (i) $\operatorname{tr}H+\operatorname{tr}(A^2)=0$; (ii) $\sum_i r_i\!\cdot\!Sr_i=0$; (iii) the pure-linear-strain validation of (B5) giving exactly $0$. Interval arithmetic is *not* needed at this stage (it belongs to PO-13).

**Repo modules to reuse.** `leray_response_relay.leray_project`, `.leray_advection`, `.helical_fejer_packet`, `.spectral_inner`; `mesoscopic_galerkin.build_angle_box_parent`, `.run_small_mesoscopic_galerkin` (time-stepper harness); `mesoscopic_local_fft.local_fft_leray_coefficients` (exact convolution cross-check); `exact_leray_relay` + `exact_carrier_record_verifier` for the rational identity checks. New code needed: a tracer/$F$-integrator and an SVD diagnostic module (~200 lines).

---

## G. Proof chain to TARGET-U (10 obligations)

1. **(L1)** Formalise: strong solution + $\liminf\sigma_3(a,t)=0$ at finite $T$ $\Rightarrow$ $\int_0^T\|\omega\|_\infty=\infty$ $\Rightarrow$ non-extendability. *(Lean-able now; depends only on ODE algebra + BKM + Kozono–Taniuchi.)*
2. **(B1),(B5),(B8)** Formalise the exact identities $\dot\sigma_i=\sigma_i(r_i\!\cdot\!Sr_i)$, the $\lambda_c$ evolution, and $\frac{d}{dt}\log\mathrm{Re}_L=\lambda_c(\Theta-2)$ for $C^2$ velocity fields.
3. **(L3)** Formalise $\frac{d}{dt}(F^{-1}\omega)=\nu F^{-1}\Delta\omega$ and the divergence of its bound — closing off the Euler-perturbative lane permanently.
4. Construct the **front flow** (seed S1) in Lagrangian variables and prove that an $s$-periodic orbit corresponds to $\langle\Theta\rangle=2$, $\gamma=1/2$.
5. **Existence** of that periodic orbit: forward RG integration to an attractor, then a computer-assisted fixed point of the period map (PO-04, radii-polynomial $Y+Z(r)<r$).
6. **Spectral tail / truncation** bounds for the front profile in a Gevrey/weighted-$\ell^1$ space (PO-05/06/07), with the $\sqrt{\nu\tau}$ viscous floor as the tail scale.
7. **Chirality lemma:** $h>0$ is preserved by the period map (the loop-D design rule as a theorem, not an assumption).
8. **Seregin2024 DSS clearance:** prove the constructed log-periodic object is outside the hypotheses of every DSS Liouville theorem (E-4), or exhibit the extra integrability failure explicitly.
9. **Nonlinear stability + entry** (PO-08/PO-09): unstable manifold of the periodic orbit intersects the image of an explicit $C_c^\infty$ non-axisymmetric datum. *Currently strategy-less repo-wide; the codimension count from step 5's Floquet spectrum is the first deliverable.*
10. **Assembly:** $T-t(s)=\int_s^\infty e^{-s'}ds'<\infty$ (PO-10) + step 1 $\Rightarrow$ $\limsup_{t\uparrow T}\|u\|_{H^m}=\infty$ (PO-11) + PO-12 (no coordinate artifact; automatic here, Cartesian throughout) $\Rightarrow$ TARGET-U.

**Honest bottom line.** Steps 1–3 are genuinely new, cheap, and Lean-able and are worth doing regardless of whether the mechanism survives. Step 5 is the same computer-assisted-proof wall every lane hits. Steps 8 and 9 are where this idea most plausibly dies.
