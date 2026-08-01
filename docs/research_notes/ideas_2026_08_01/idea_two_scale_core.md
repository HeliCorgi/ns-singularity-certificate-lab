# LENS 12 — Two-scale core: the O(1) sea, the inertial bridge, and the autonomy of the inner front

**Status: FORMAL ANSATZ** (mechanism). The decoupling estimates §B4–B6 are separately **PROOF CANDIDATE** (finite-Fourier algebra + one integration by parts; Lean-able). Six sub-variants died during derivation and are kept below marked **REJECTED** with their exact failing equation.

Conventions are the repo's: $\mathbb T^3$, $u=\sum_k\hat u_ke^{ik\cdot x}$, $\partial_t\hat u_k=\mathcal N_k-\nu|k|^2\hat u_k+\hat f_k$, $\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$, $P_k=I-k\otimes k/|k|^2$, critical shell normalization $E_N=\tfrac12\|u_N\|_2^2=c_E/N$, $h_{\nu,\tau}(r)=(1-e^{-\nu\tau r^2})/(\nu r^2)$.

---

## A. Clay target

**Primary: negation of (B)** — finite-time breakdown on $\mathbb T^3=[0,2\pi)^3$, **no forcing**. **Fallback: (D)** — same with a smooth force $f$ of *finite* Fourier support $|k|\le K_0=O(1)$, time-independent, mean-zero (admissible: $C^\infty$, and trivially satisfies F‑N4 since $\hat f_k\equiv0$ for $|k|>K_0$; only the F‑N5 *indirect* channel is invoked). The decoupling theorem of §D is exactly the statement that $f$ is irrelevant, so (D) is a strictly weaker corollary and the honest target is (B).

Initial data: $u_0\in C^\infty(\mathbb T^3;\mathbb R^3)$, $\nabla\!\cdot u_0=0$, mean zero, **not axisymmetric** (see E-13). Regularity class $C([0,T);H^m)\cap C^1([0,T);H^{m-2})$, $m>5/2$. Viscosity $\nu>0$ **fixed and never rescaled**; the free parameters are the critical amplitude $c_E$ and the entry scale $N_0$. Total energy $E(0)\approx 2c_E/N_0$ is made *small* by starting deep even though $c_E$ is large — this is the S3 c_E-collapse, used here in the only place it is legal.

Breakdown signature: $\limsup_{t\uparrow T}\|u(t)\|_{L^3}=\infty$ and $\int_0^T\|\omega\|_\infty dt=\infty$. Energy and enstrophy are *not* used (F‑N1/N2).

---

## B. Central mathematics

### B0. The two blocks

$u=U+v$: outer **sea** $U=P_{\le K_0}u$, $\tfrac12\|U\|_2^2=O(1)$, strain $S_{\rm out}=\mathrm{sym}\nabla U$, $S_0:=\|S_{\rm out}\|_\infty=O(1)$; inner **front** $v$ at scale $1/N(t)$. Between them the **bridge** $E_k=c_E/k$ (dyadic), i.e. $\mathcal E(k):=dE/dk=2c_E k^{-2}$.

**REJECTED V0 (direct outer drive).** Linear stretching of the front by the sea: $\frac{d}{dt}\|\omega_{\rm fr}\|\le (S_0-\nu N^2)\|\omega_{\rm fr}\|$. Fails for every $N>(S_0/\nu)^{1/2}$. Exact failing equation: $S_0-\nu N^2<0$.

**REJECTED V1 (single helper at scale $m$, critical energy).** Helper amplitude $\sqrt{2c_E/m}$ spread over the box gives strain $S_m=\sqrt{2c_Em}$; need $S_m\ge\nu N^2$, i.e. $m\ge\nu^2N^4/(2c_E)$, but $m\le N$. Fails for $N^3>2c_E/\nu^2$. Exact failing equation: $\nu^2N^4/(2c_E)>N$.

### B1. The bridge is a linear drain, not a conveyor — **LD no-go**

Shell budget $\partial_t\mathcal E(k)+\partial_k\Pi(k)=-2\nu k^2\mathcal E(k)$. On the critical bridge $\mathcal E=2c_Ek^{-2}$ the dissipation density is **$k$-independent**: $2\nu k^2\cdot2c_Ek^{-2}=4\nu c_E$. Hence, quasi-statically,
$$\boxed{\ \Pi(k)=\Pi(1)-4\nu c_E\,(k-1).\ }\tag{B1}$$
Dimension check: $[\nu c_E]=[\text{energy}]/[\text{time}]/[k]$ ✓; sign: monotone decreasing ✓. With $\Pi(1)\le\varepsilon_0=O(1)$ (the sea's turnover power):

**REJECTED V2 (strain conveyor from the sea).** $\Pi(N)=\varepsilon_0-4\nu c_E(N-1)<0$ for all $N>1+\varepsilon_0/(4\nu c_E)=O(1)$. The sea can feed *at most $O(1)$ octaves*. The "power-law inertial bridge as strain conveyor" in the lens brief is **false in the sea$\to$front direction**.

### B2. What survives: octave-local conveyance

The conveyor works, but only *one octave at a time*, and only if each octave is **point-localized**. Let octave $j$ occupy volume fraction $\phi_j$ with amplitude $A_j$; critical normalization gives $\phi_jA_j^2=2c_E/N_j$. Its strain is $S_j=A_jN_j$. Then
$$\frac{S_j}{\nu N_{j+1}^2}=\frac{A_jN_j}{4\nu N_j^2}=\frac{A_j}{4\nu N_j},\qquad A_j=\sqrt{\tfrac{2c_E}{N_j\phi_j}}.$$
This is $N$-independent **iff $\phi_j=c_\ell^3N_j^{-3}$**, i.e. the octave is a *single three-dimensionally localized packet of physical extent $\ell_j=c_\ell/N_j$*. Then $A_j=\alpha N_j$, $\alpha:=\sqrt{2c_E/c_\ell^3}$, and
$$\boxed{\ \frac{S_j}{\nu N_{j+1}^2}=\frac{\alpha}{4\nu}=\frac{1}{4\nu}\sqrt{\frac{2c_E}{c_\ell^3}}\ \ \text{(scale-free)},\qquad \frac{S_j}{\nu N_{j+n}^2}=\frac{\alpha}{4^n\nu}.}\tag{B2}$$
**Conveyor range:** $n\le R_{\rm conv}=\log_4(\alpha/\nu)$ octaves. The sea sits at $n=\log_2 N$, so it is out of range by an exponential factor — this is (B1) again, now octave-resolved.

**REJECTED V4 (delocalized / box-filling front).** $\phi=O(1)$ gives $A=\sqrt{2c_E/N}\to0$, $\|u\|_\infty=O(1)$, and $\|u_{N_j}\|_3^3=\phi A^3=(2c_E/N_j)^{3/2}$, summable $\Rightarrow \|u\|_{L^3}$ bounded $\Rightarrow$ ESS $\Rightarrow$ regular. Note this member *passes* the $M^{\rm eff}_N\gtrsim N^3$ screen of §6 — the screen alone does not select it out; ESS does.

**REJECTED V6 (filament / sheet fronts).** Filament ($\phi=N^{-2}$): $\|u_N\|_3^3=\phi A^3=(2c_E)^{3/2}N^{-1/2}\to0$. Sheet ($\phi=N^{-1}$): $=(2c_E)^{3/2}N^{-1}\to0$. Both give convergent octave sums $\Rightarrow$ bounded $L^3\Rightarrow$ regular. **Only the point-like packet gives $\|u_{N_j}\|_3^3=c_\ell^{-3/2}(2c_E)^{3/2}=$ const per octave.** Its Fourier ball has radius $N/c_\ell$, i.e. relative width $\eta=1/c_\ell$; taking $c_\ell=3.5$ puts it inside the *only* structurally uncondemned window $\eta\in(0,1/3)$ of the mesoscopic lane.

### B3. Inner problem and the front flow (S1, with the IR condition derived)

Front ansatz $\hat v(k,t)=N(t)^{-2}\Psi(\xi,s)$, $\xi=k/N$, $ds=N^2dt$, $a:=\dot N/N^3$. Substituting ($\dot N=aN^3$):
$$\partial_t\hat v=\Psi_s-a(2\Psi+\xi\!\cdot\!\nabla_\xi\Psi),\quad \nu|k|^2\hat v=\nu\xi^2\Psi,\quad \mathcal N\to-\mathcal Q(\Psi,\Psi),$$
$$\partial_s\Psi=a(2\Psi+\xi\!\cdot\!\nabla_\xi\Psi)-\nu\xi^2\Psi-\mathcal Q(\Psi,\Psi).\tag{B3}$$
Every term is $O(N^0)$ **precisely because** $E_N=c_E/N$: $\sum_k|\hat v_k|^2=N^{-1}\!\int|\Psi|^2d^3\xi$, so $\int|\Psi|^2=2c_E$. The dimensionless nonlinear/viscous ratio in (B3) is $\mathcal Q/(\nu\xi^2\Psi)\sim\sqrt{c_E}/\nu$ — **$N$-free**, the same threshold as S2's $\chi\sqrt{c_E}\ge\sqrt2\,\nu$ and as (B2).

Shell density $e(\xi,s)=\xi^2\oint|\Psi(\xi\omega)|^2d\omega$; from (B3), $\partial_se=a(2e+\xi e')-2\nu\xi^2e-\partial_\xi F$. Writing $a(2e+\xi e')=\partial_\xi(a\xi e)+ae$ puts it in conservation form with **$\xi$-flux $J=F-a\xi e$**; characteristics $d\xi/ds=-a\xi$ move modes *downward* in $\xi$ — the wake is produced by **dilation drift**, not by inertial flux.

**Wake law (derived, not assumed).** Below the front base $\xi_0$ set $F\equiv0$ (inert wake, see M2). Stationarity gives $a\xi e'=-2ae+2\nu\xi^2e$, hence
$$\boxed{\ e(\xi)=e(\xi_0)\Big(\frac{\xi_0}{\xi}\Big)^{2}\exp\!\Big[-\frac{\nu}{a}(\xi_0^2-\xi^2)\Big].}\tag{B4}$$
As $\xi\to0$: $e\to 2c_E\xi^{-2}e^{-\nu\xi_0^2/a}$ — **the critical wake $e_c=2c_E\xi^{-2}$ of S2 is an output of (B3), not an input**, with a scale-independent retention $\theta=e^{-\nu\xi_0^2/a}>0$. Equivalently $|\Psi|\sim\xi^{-2}$, i.e. $\hat u_k\sim k^{-2}$ *independent of $N$* ✓ (the wake is frozen in physical $k$, as it must be). Cross-check: $\Phi=\int_{K_0/N}e\,d\xi\approx 2c_EN/K_0$, so wake+front energy $=\Phi/(2N)=c_E/K_0=O(1)$, finite ✓.

### B4. Matching at the front base $\xi_0$ — the outer boundary condition

**(M1) Infrared spectral condition.** $\Psi(\xi,s)\to\mathcal C(\hat\xi,s)\,\xi^{-2}$ as $\xi\downarrow0$, with the power law truncated at $\xi=K_0/N$ where the sea's $O(1)$ spectrum takes over. **As $N\to\infty$ the matching point recedes to $\xi=0$**: the sea is expelled to the boundary of the inner domain and survives only as the (integrable) IR cutoff of $\int\xi^2e\,d\xi$.

**(M2) Flux node.** $F(\xi_0^+,s)=0$: the wake is inertially inert (phases frozen, heat-only), so no radial flux crosses the base. This *replaces* S2's $F(1)=\chi(2c_E)^{3/2}>0$.

**REJECTED V5 (S2 steady-wake closure).** $F(1)=\chi(2c_E)^{3/2}>0$ requires an inertial source at $\xi<1$; by (B1) that source would have to be the sea, which delivers $\Pi(N)<0$. Hence S2's front-extent formula $\xi_{\max}=1+(\chi/\sqrt2)\sqrt{c_E}/\nu$ is superseded; the correct BC is (M2) and the front extent is fixed by shape recurrence.

**(M3) RDT matching term — the actual outer BC.** Expand the sea about the front centre $x_*(t)$: $U(x)=U(x_*)+S(t)(x-x_*)+O(|x-x_*|^2)$, and choose the frame $\dot x_*=U(x_*,t)$ (this removes sweeping exactly; sweeping is energy-neutral, see B5). The residual linear strain enters (B3) as a rapid-distortion term. In physical time it is $O(S_0)$; in the $s$-clock $\partial_s=N^{-2}\partial_t$, so
$$\boxed{\ \partial_s\Psi=a(2\Psi+\xi\!\cdot\!\nabla_\xi\Psi)-\nu\xi^2\Psi-\mathcal Q(\Psi,\Psi)\;+\;\frac{1}{N^{2}}\Big[(S^{\!\top}\xi)\!\cdot\!\nabla_\xi\Psi-P_\xi S\Psi\Big]+O(N^{-3}).}\tag{B5}$$
The entire sea enters the inner problem through **one traceless symmetric $3\times3$ matrix $S(t)$ (5 numbers) at order $N^{-2}$**; the curvature of $U$ enters at $O(N^{-3})$ ($|x-x_*|^2\sim N^{-2}$, one $\nabla$ giving $N$, one $N^{-2}$ from the clock).

### B5. Sea $\to$ front energy flux (exact)

Split $\mathcal N_k$ by which leg is the sea. The advecting-leg term is $-iP_k\sum_{|\ell|\le K_0}((k-\ell)\!\cdot\!\hat U_\ell)\hat v_{k-\ell}$ ($=-P\,\widehat{(U\!\cdot\!\nabla)v}$); the advected-leg term is $-iP_k\sum_{|m|\le K_0}(m\!\cdot\!\hat v_{k-m})\hat U_m$ ($=-P\,\widehat{(v\!\cdot\!\nabla)U}$). Energy transfer into the front band:
$$\Theta_{sf}=\mathrm{Re}\!\!\sum_{|k|\sim N}\!\!\hat v_k^{*}\!\cdot\!\mathcal N_k^{sf} =-\!\int\! v\!\cdot\!(U\!\cdot\!\nabla)v-\!\int\! v\!\cdot\!(v\!\cdot\!\nabla)U =-\tfrac12\!\int\! U\!\cdot\!\nabla|v|^2-\!\int\! v\!\otimes\! v:\nabla U,$$
and the first integral is **exactly zero** ($\nabla\!\cdot U=0$, periodic). Hence
$$\boxed{\ \Theta_{sf}(t)=-\int_{\mathbb T^3}(v\otimes v):S_{\rm out}\,d^3x,\qquad |\Theta_{sf}|\le S_0\|v\|_2^2=\frac{2S_0c_E}{N}.\ }\tag{B6}$$
*Sweeping transfers exactly zero energy; only the traceless strain couples.* Compare with the front's own budget rate $|dE_{\rm fr}/dt|=c_E|\dot N|/N^2=c_EaN$:
$$\frac{|\Theta_{sf}|}{c_EaN}\le\frac{2S_0}{aN^2}=O(N^{-2}).\tag{B7}$$
Lifetime-integrated sea-induced strain on the front: $\int_{t_J}^{T}S_0\,dt=S_0(T-t_J)=S_0/(2aN_J^2)=O(N_J^{-2})$, whereas viscosity accumulates $\nu N^2\Delta t=\nu\log2/a=O(1)$ per doubling. **The sea is weaker than viscosity by $N^{-2}$ — and viscosity is only marginally beaten, by the front's own strain (B2).**

### B6. Front $\to$ sea back-reaction (explicit integral)

Pressure: $\hat p_k=-\frac{k_ik_j}{|k|^2}\widehat{(u_iu_j)}_k$. The front's contribution to the low band, with $v$ localized at $x_*$ on scale $\ell=c_\ell/N$:
$$P_{\le K_0}p^{(v)}(x)=-\!\!\sum_{0<|k|\le K_0}\!\!\frac{k_ik_j}{|k|^2}e^{ik\cdot x}\!\!\int_{\mathbb T^3}\!\!e^{-ik\cdot y}v_i(y)v_j(y)\,d^3y =-\!\!\sum_{0<|k|\le K_0}\!\!\frac{k_ik_j}{|k|^2}e^{ik\cdot(x-x_*)}\big[\mathfrak M_{ij}+O(|k|\ell\,\mathfrak M)\big],$$
$$\mathfrak M_{ij}:=\int v_iv_j\,d^3y,\qquad \mathrm{tr}\,\mathfrak M=\|v\|_2^2=\frac{2c_E}{N}.$$
**The front acts on the sea as a point quadrupole of strength $\mathfrak M=O(c_E/N)$, with the first correction $O(c_E/N^2)$.** Hence, with $C_{K_0}:=(\sum_{0<|k|\le K_0}|k|^2)^{1/2}\asymp(4\pi/5)^{1/2}K_0^{5/2}$,
$$\|\nabla P_{\le K_0}p^{(v)}\|_{L^2}\le C_{K_0}\frac{2c_E}{N}=O(N^{-1}).\tag{B8}$$
Integrating over the *entire remaining life* of the cascade, using $N(t)=(2a(T-t))^{-1/2}$:
$$\|\delta U\|_{L^2}\le\int_{t_J}^{T}\!\!2c_EC_{K_0}\sqrt{2a}\,(T-t)^{1/2}dt=\frac{4}{3}c_EC_{K_0}\sqrt{2a}\,(T-t_J)^{3/2} =\boxed{\ \frac{2\,c_EC_{K_0}}{3\,a\,N_J^{3}}=O(N_J^{-3}).}\tag{B9}$$
The direct nonlinear back-reaction is the exact antisymmetric partner of (B6), $\Theta_{fs}=-\Theta_{sf}$, and integrates to $\int_{t_J}^T|\Theta_{fs}|dt\le\frac{2S_0c_E}{3aN_J^3}$ — **same order $O(N_J^{-3})$**.

---

## C. Scaling table ($\tau=T-t$, $N=(2a\tau)^{-1/2}$, $\gamma=1/2$; $A=\alpha N$, $\ell=c_\ell/N$, $\alpha=\sqrt{2c_E/c_\ell^3}$)

| quantity | law | $\tau$-exponent |
|---|---|---|
| total energy $E(t)$ | $\to E(T)>0$, bounded (Leray) | $0$ |
| front energy $c_E/N$ | $c_E\sqrt{2a}\,\tau^{1/2}$ | $+1/2$ |
| enstrophy $\|\omega\|_2^2$ | $\asymp2c_EN$ | $-1/2$ |
| dissipation rate $\nu\|\nabla u\|_2^2$ | $\asymp2\nu c_EN$ | $-1/2$ (so $\int_0^T<\infty$ ✓) |
| global $\|u\|_{L^3}^3$ | $c_3\log_2N+O(1)$, $c_3=c_\ell^{-3/2}(2c_E)^{3/2}\theta$ | $\log(1/\tau)$ |
| $\|u\|_{\dot H^{1/2}}^2$ | $\asymp2c_E\theta\log_2N$ | $\log(1/\tau)$ |
| $\|u\|_{L^\infty}$ | $\alpha N$ | $-1/2$ |
| $\|\omega\|_{L^\infty}$ | $(\alpha/c_\ell)N^2$ | $-1$ |
| $\int_0^t\|\omega\|_\infty$ | $\propto\log(1/\tau)$ | divergent ✓ (BKM) |
| $\sqrt{\tau}\|u\|_\infty$ | $\to\alpha/\sqrt{2a}$ | $0$ (**Type-I marginal**) |
| nonlinear $\|P(u\!\cdot\!\nabla)u\|_2$ | $\alpha^2c_\ell^{1/2}N^{3/2}$ | $-3/4$ |
| pressure $\|\nabla p\|_2$ | $\le$ same | $-3/4$ |
| viscous $\nu\|\Delta u\|_2$ | $\nu\alpha c_\ell^{-1/2}N^{3/2}$ | $-3/4$ (**same order** — ratio $\alpha/\nu$, $N$-free) |
| physical time remaining | $\tau=1/(2aN^2)$ | — |
| Fourier bandwidth | $\kappa N$ | $-1/2$ |
| active mode count $M^{\rm eff}$ | $(N/c_\ell)^3$ | $-3/2$ |

---

## D. Closed feedback loop (every arrow a formula)

$$
\underbrace{\Psi_*(\xi,s)}_{\text{shape}}\ \xrightarrow{\ \mathcal Q(\Psi,\Psi)\ }\ \underbrace{F(\xi)>0\text{ on }[\xi_0,\xi_{\max}]}_{\text{forward flux}}\ \xrightarrow{\ \text{(B2): }S_j/\nu N_{j+1}^2=\alpha/4\nu\ }\ \underbrace{\text{octave }j{+}1\text{ built}}_{A_{j+1}=\alpha N_{j+1}}
$$
$$
\xrightarrow{\ d\xi/ds=-a\xi\ }\ \underbrace{e(\xi)=e(\xi_0)(\xi_0/\xi)^2e^{-\frac{\nu}{a}(\xi_0^2-\xi^2)}}_{\text{(B4) wake, retention }\theta=e^{-\nu\xi_0^2/a}}\ \xrightarrow{\ \text{S4 per-octave}\ }\ \|u\|_3^3\ge c_3\log_2N
$$
$$
\xrightarrow{\ a=\log2/S>0\ }\ N^{-2}(t)=N_0^{-2}-2\!\!\int_0^t\!\!a\ \Rightarrow\ T\le\frac{N_0^{-2}}{2a_-}<\infty\ \xrightarrow{\ \text{doubling pullback}\ }\ \Psi_*(\xi,s+S)=\Psi_*(\xi,s)\ \circlearrowleft
$$
Two **damping** arrows close the loop from outside and are shown small:
$$\text{sea}\to\text{front}:\ \frac{1}{N^2}\big[(S^{\!\top}\xi)\!\cdot\!\nabla_\xi\Psi-P_\xi S\Psi\big]\ \ (\text{B5}),\qquad \Theta_{sf}=-\!\int\!(v\otimes v):S_{\rm out},\ \ \frac{|\Theta_{sf}|}{c_EaN}\le\frac{2S_0}{aN^2};$$
$$\text{front}\to\text{sea}:\ \hat p^{(v)}_k=-\frac{k_ik_j\mathfrak M_{ij}}{|k|^2}e^{-ik\cdot x_*}+O(c_EN^{-2}),\qquad \|\delta U\|_2\le\frac{2c_EC_{K_0}}{3aN_J^3}.$$
**Decoupling statement.** The front is autonomous to relative accuracy $O(N_J^{-2})$ over its entire remaining lifetime; the sea is perturbed by $O(N_J^{-3})$ over the same interval. Both are summable in $J$. Consequently the inner problem (B3) with IR condition (M1) and flux node (M2) is a *closed* renormalization-group flow, and **the only surviving obstruction is shape recurrence: does (B3) admit an $s$-periodic $\Psi_*$ with $a>0$?** (S3's $\chi_{\rm shape}>0$, re-derived here by an independent route.)

**Verdict on cloud-candidate unproved Lemma 7 (low-band controller).** **CONFIRMED — and strengthened from "unnecessary" to "impotent."** By (B6)–(B7) a low-band controller can inject at most $2S_0c_E/N$ of energy per unit time against a front budget $c_EaN$, ratio $O(N^{-2})$; and by (B1) any attempt to route power up a critical bridge dies at $k\approx\varepsilon_0/(4\nu c_E)$. The corollary is a *negative* one for the cloud lane: **no low-band engineering can repair a front with $\chi_{\rm shape}\le0$.** Lemma 7 is removed from the obligation list and replaced by nothing.

---

## E. Obstruction audit (exact collision points)

- **E-1 Energy bound / F‑N1.** No collision: $E(t)$ bounded by construction (table row 1). Not used as signature.
- **E-2 Finite dissipation / F‑N2.** No collision: $\nu\!\int_0^T\!\|\nabla u\|_2^2dt=2\nu c_E\!\int_0^{\tau_0}\!(2a\tau)^{-1/2}d\tau<\infty$. Not used as signature.
- **E-3 ESS $L^\infty_tL^3_x$.** **Must be violated and is** — but only logarithmically, $\|u\|_3^3\ge c_3\log_2N$. Collision point: $c_3=c_\ell^{-3/2}(2c_E)^{3/2}\theta$ with $\theta=e^{-\nu\xi_0^2/a}$. If wake retention $\theta\to0$ (front consumes its wake) the candidate dies — see **REJECTED V3**: a wake-consuming travelling front has $O(1)$ active octaves, $\|u\|_3^3=O(1)$, ESS $\Rightarrow$ regular.
- **E-4 Serrin $(p,q)=(\infty,2)$.** $\|u\|_\infty\asymp\alpha N\asymp\tau^{-1/2}\Rightarrow\int_0^T\|u\|_\infty^2dt\asymp\int\tau^{-1}d\tau=\infty$ ✓ (marginally). This is the closest non-endpoint escape and it is *exactly* logarithmic; V4/V6 die here.
- **E-5 Fixed-finite-bandwidth no-go (F‑α1, VR‑L‑011, `galerkin_not_tendsto_atTop`).** Evaded: the *inner* support $[\xi_0,\xi_{\max}]$ is fixed but the physical band $N\cdot[\xi_0,\xi_{\max}]$ and mode count $(N/c_\ell)^3$ diverge like $\tau^{-1/2},\tau^{-3/2}$. Collision point: the no-go is about a *fixed* set $S\subset\mathbb Z^3$; ours is $s$-dependent with $|S(s)|\to\infty$.
- **E-6 Galerkin global existence.** Same evasion as E-5; the truncated system is globally regular for every fixed cutoff and the limit is not uniform — this is why the pilot (§F) may *never* be read as evidence of blowup.
- **E-7 Pure-swirl $L^3$ no-go (VR‑L‑016, LG‑9).** Collision point: for a rotationally equivariant $\Psi$ about a fixed axis, $u\cdot\nabla p_0\equiv0$, $\chi_{\rm shape}\le0$ and $F(\xi)\le0$. **Design rule: $\Psi_*$ must break axisymmetry; add the exact rational $J(u_0)=3\!\int\!p\nabla\!\cdot(|u|u)$ column to the pilot.**
- **E-8 One-scale backward self-similar (NRS1996 / Tsai1998).** **This is the sharpest collision.** With $a=$ const and $\Psi$ *$s$-independent*, $\hat u_k=N^{-2}\Psi(k/N)$ is an exact Leray profile with $\gamma=1/2$ and NRS gives $\Psi\equiv0$ under $U\in L^3$. Evasion is **two-fold and must be stated together**: (i) $\mathbb Z^3$ breaks continuous dilation, so the recurrent object is $S$-periodic with $S=\log2/a$, i.e. **discretely self-similar with ratio $2$**, not self-similar; (ii) decisively, the rescaled field $V(y,s)=\sqrt{2a\tau}\,u(x_*+\sqrt{2a\tau}\,y,t)$ has $\|V\|_{L^3}^3=\|u\|_{L^3}^3\ge c_3\log_2N\to\infty$, so **the NRS hypothesis $U\in L^3$ fails by exactly the same logarithm that supplies the ESS violation.** Tsai's local-energy version: $\int_{B_1}|V|^2$ *is* finite ($=O(c_E)$, only the top octave sits in $B_1$), so Tsai is evaded **only** by (i). Honest flag: Seregin2024/2026 DSS-Liouville preprints are the live residual threat here; the repo forbids excluding by abstract, and equally forbids *assuming* they do not apply.
- **E-9 CSTY2009 Type-I exclusion.** $\sqrt{\tau}\|u\|_\infty\to\alpha/\sqrt{2a}$: the candidate **is Type I**. CSTY does not apply — it is an *axisymmetric $\mathbb R^3$* theorem and our field is neither. Collision point is therefore fully loaded onto E-13: axisymmetry must be broken *structurally*, not accidentally.
- **E-10 KNSS ancient-solution Liouville.** Type-I $\Rightarrow$ the rescaled limit is a nonzero bounded ancient mild solution; ours is $V$, $S$-periodic, with wake octave $n=J-j$ sitting at $|y|\sim2^{n}$ with amplitude $\alpha2^{-n}$, i.e. **$|V(y,s)|\lesssim\alpha/|y|$ — exactly the KNSS(b) decay $C/r$.** KNSS(b) then forces $V=0$ **if $V$ is axisymmetric**. Collision point named precisely: *the candidate is killed by KNSS the instant the front shape is axisymmetric with swirl.* Non-axisymmetry is not cosmetic; it is the sole barrier at two independent no-gos (E-7, E-10).
- **E-11 Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}_N/N^3$.** We sit at $\gamma=1$: $M^{\rm eff}_N=(N/c_\ell)^3$, so $D_N\le2\kappa^2\tau^2c_E/c_\ell^3$, $N$-independent ✓ (screen passed, not violated). Required half-transfer $D_N=1/2$ gives the **necessary** condition $\ c_E\ge c_\ell^3/(4\kappa^2\tau^2)$; with $\kappa=2$, $c_\ell=3.5$, $\tau=1/4$ this is $c_E\ge 107$ — the same order as the registered "$c_E\approx228$ required" row. Legal by S3 (total energy $2c_E/N_0$ small for large $N_0$), but it is an *upper* bound: sufficiency still requires the shape.
- **E-12 Diagonal cross-talk gate.** No collision in its original form: there is **one** structure per octave, hence no two same-scale relays and no diagonal parent pairs $A_1{+}B_2$. The surviving analogue is front($N$)$\times$wake($N/2$) $\to$ bands $N/2$ and $3N/2$, which is *not* leakage but is also *not small*; it is fully contained inside $\mathcal Q$ in (B3) provided the inner domain extends to $\xi\ll1$. **Own unproved lemma:** (M2) $F(\xi_0)=0$ is an *output to be measured*, not an input; predicted $F(\xi_0)=O(\nu)$.
- **E-13 Front-resolution threat model (TM-22, PREREGISTERED_MIN_POINTS_PER_FRONT=7).** The front is $\ell=3.5/N$ wide, so a grid with $\ge7$ points across $\ell$ needs $\Delta x\le0.5/N$, i.e. $\ge 2N$ points per direction *inside the packet* — the pilot below fixes this by working in Fourier space with the packet's own band, never by fitting a physical-space peak.
- **E-14 Smooth-forcing high-frequency decay (F‑N4/N5).** No collision: $\hat f$ supported in $|k|\le K_0$, never injects into high shells; and by the decoupling theorem it does not drive them indirectly either — which is why the honest target is unforced (B), and why F‑N5's toy indirect cascade is *specifically contradicted* by (B1) at critical normalization.

---

## F. Minimal falsification experiment (SFD‑1, $\le$1 hour)

**Question tested:** is the inner front autonomous (Lemma-7 verdict) and is the bridge a drain?

**Configuration.** $\mathbb T^3$, grid $128^3$, 2/3 dealiasing. Sea: $|k|\le K_0=2$, random div-free, energy $1/2$, giving a measured $S_{\rm out}$. Front: single localized real div-free packet, $\ell=3.5/N$, energy $c_E/N$, $N\in\{8,12,16,20,24\}$. $\nu\in\{1/40,1/80\}$, $c_E\in\{1,10\}$.

**Variables / success criteria** (fit $\log(\cdot)$ vs $\log N$, $\ge4$ points, require $R^2\ge0.95$):
1. **V1 (sea$\to$front).** Measure $\Theta_{sf}$ two ways — spectral $\mathrm{Re}\sum_{|k|\sim N}\hat v_k^*\mathcal N_k^{sf}$ and the closed form $-\!\int(v\otimes v):S_{\rm out}$ — agreement to $\le10^{-12}$ relative (**convention check**), and fitted exponent of $|\Theta_{sf}|/(c_EaN)$ in $[-2.3,-1.7]$.
2. **V2 (front$\to$sea).** $\|\nabla P_{\le K_0}p^{(v)}\|_2$ exponent in $[-1.3,-0.7]$.
3. **V3 (quadrupole).** Relative error of $\hat p^{(v)}_k$ vs $-k_ik_j\mathfrak M_{ij}|k|^{-2}e^{-ik\cdot x_*}$ decays with fitted exponent in $[-1.3,-0.7]$.
4. **V4 (LD law).** Build the critical bridge $E_k=c_E/k$, $1\le k\le N$; measure exact Leray flux $\Pi(k)$; fit $\Pi(k)-\Pi(1)$ vs $k$: slope must be **negative and linear**, matching $-4\nu c_E$ within $20\%$.

**Kill conditions (pre-registered).** V1 exponent $>-1.5$ $\Rightarrow$ the front is *not* autonomous, Lemma 7 is resurrected, **this lens dies**. V4 slope not negative-linear $\Rightarrow$ (B1) is wrong and V2's rejection is void. V1's two-way identity mismatch $>10^{-10}$ $\Rightarrow$ sign/convention bug, halt.

**Arithmetic.** V1–V4 exponents: binary64 is sufficient (they are ratios, not cancellations). **Exact rational required** for: the identity $\Theta_{sf}=-\int(v\otimes v):S_{\rm out}$ on a tiny configuration ($N=4$, $K_0=1$, $\le9$ wavevectors) — this is an algebraic convention check on $\mathcal N_k=-iP_k\sum(m\cdot\hat u_\ell)\hat u_m$ and must not be validated in floats; and the $\mathfrak M_{ij}$ trace identity.

**Reuse (repo API).** `leray_response_relay.leray_project`, `.leray_advection` (dealiased $P((u\!\cdot\!\nabla)u)$), `.spectral_inner`, `.gradient_l2_squared`; `mesoscopic_local_fft.local_fft_leray_coefficients` (zero-padded exact linear convolution — required for $\Pi(k)$ in V4, no wrap aliasing); `mesoscopic_cloud_scaling.exact_sparse_leray_convolution` (sparse sea block); `mesoscopic_galerkin.build_angle_box_parent` (adapt: replace the two angled boxes by one localized packet, keep the `scale > 3*(width-1)` disjointness guard as $c_\ell>3$); `exact_leray_relay` + `exact_carrier_record_verifier` for the exact-rational leg. **Do not** reuse `leray_response_relay.relay_stage`'s Fejér orbit family (rejected route).

**Permitted claim ceiling.** At most: "sea/front coupling exponents consistent with the derived $N^{-2}/N^{-1}$ decoupling at the stated resolution." Never a blowup claim.

---

## G. Proof chain (10 obligations)

1. **PO‑03** Local existence/uniqueness, $C([0,T);H^m)\cap C^1([0,T);H^{m-2})$, $m>5/2$ on $\mathbb T^3$ (classical; Lean gap F‑δ noted).
2. **Front-flow well-posedness.** (B3) locally well-posed in a weighted $\ell^1_\nu$/Gevrey space with IR weight enforcing $\Psi\sim\xi^{-2}$ (M1); $\mathcal Q$ bounded via the Banach-algebra estimate (PO‑07 style).
3. **Recurrent shape.** Existence of $\Psi_*$ with $\Psi_*(\cdot,s+S)=\Psi_*(\cdot,s)$, $S=\log2/a$, $a\ge a_->0$, $F>0$ on $[\xi_0,\xi_{\max}]$, non-axisymmetric. Route: forward RG integration (S1/S7) then radii-polynomial $Y+Z(r)<r$ interval fixed point (PO‑04, PO‑13).
4. **Decoupling theorem** (this document, §B4–B6): (B5)+(B6)+(B8)+(B9). Finite-Fourier algebra + one integration by parts; formalizable.
5. **Structural stability.** Floquet spectrum of $D\mathcal F(\Psi_*)$: finitely many unstable directions, no marginal modes beyond the gauge group (translation, dilation, rotation, phase) $\Rightarrow$ the $O(N^{-2})$-perturbed system (B5) shadows $\Psi_*$ for all $s\ge s_0$ (PO‑08).
6. **Wake persistence & per-octave $L^3$.** (B4) gives retention $\theta=e^{-\nu\xi_0^2/a}>0$; S4's Hölder+Bernstein per-octave lemma gives $\|u_{N_j}\|_3^3\ge c_3>0$; hence $\|u(t)\|_3^3\ge c_3\log_2N(t)-C$.
7. **Finite blow-up time.** $N^{-2}(t)=N_0^{-2}-2\int_0^t a\Rightarrow T\le N_0^{-2}/(2a_-)<\infty$ (PO‑10).
8. **Entry (PO‑09).** Explicit $u_0\in C^\infty(\mathbb T^3)$ (trigonometric polynomial + controlled tail) whose orbit meets the stable manifold of $\Psi_*$. Reduced by (5) to a **finite-codimension shooting problem**; this is the chain's weakest link and the decoupling theorem is what makes it finite-codimension rather than infinite.
9. **Norm divergence & non-extendability (PO‑11, PO‑12).** $\limsup_{t\uparrow T}\|u\|_{L^3}=\infty$ (from 6), $\int_0^T\|\omega\|_\infty dt=\infty$ (table), $\Rightarrow$ by Kato–Ponce/BKM no extension in $H^m$ past $T$.
10. **Clay.** Negation of (B) on $\mathbb T^3$, unforced; (D) as corollary. Whole-space (A)/(C) transfer is *plausible but not claimed*: the decoupling theorem says the front never sees the domain, but HS‑5 forbids importing the periodic closure — a separate obligation.

---

### Retained REJECTED sub-variants (with failing equations)

| id | variant | exact failing equation |
|---|---|---|
| V0 | outer $O(1)$ strain drives the front | $S_0-\nu N^2<0$ for $N>(S_0/\nu)^{1/2}$ |
| V1 | single helper at scale $m$, critical energy | $m\ge\nu^2N^4/(2c_E)>N$ |
| V2 | sea-fed strain conveyor along $E_k=c_E/k$ | $\Pi(N)=\varepsilon_0-4\nu c_E(N-1)<0$ |
| V3 | wake-consuming travelling front | $\|u\|_3^3=O(1)\Rightarrow$ ESS $\Rightarrow$ regular |
| V4 | delocalized/box-filling front | $\|u_{N_j}\|_3^3=(2c_E/N_j)^{3/2}$, summable |
| V5 | S2 steady-wake closure $F(1)=\chi(2c_E)^{3/2}>0$ | needs source at $\xi<1$, excluded by V2 |
| V6 | filament ($\phi=N^{-2}$) or sheet ($\phi=N^{-1}$) front | $\|u_N\|_3^3=(2c_E)^{3/2}N^{-1/2}$ resp. $N^{-1}\to0$ |
