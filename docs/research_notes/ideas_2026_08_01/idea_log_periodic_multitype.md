# LENS 8 — Log-periodic renormalisation beyond period 1: the Butterfly carrier cycle with odd-denominator scale stagger

**Status: FORMAL ANSATZ**, containing two **SYMBOLIC CANDIDATE** sub-results proved here in closed form:
(LP-4) an explicit rational polarization witness realising a *both-channels-nonzero* carrier pair, and
(LP-9) an exact divisibility theorem making diagonal cross-talk miss every tagged carrier box at every scale.
Direct continuation of eq. (6.11) of `CANDIDATE_SOLUTION_PHASE_CODED_LERAY_CLOUD.md`.

---

## A. Clay target

**Target: (B)** — global regularity on $\mathbb T^3$, **no forcing** (a *negation* candidate: we construct a mechanism whose success would refute (B); its failure is a regularity datum). Because no forcing is used, (D) follows a fortiori and F-N1/F-N2/F-N4 are vacuous here.

- Domain: $\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3$, $u(x)=\sum_k\hat u_ke^{ik\cdot x}$, $k\in\mathbb Z^3$, zero mean, $\mathbb P_0=0$.
- Forcing: none. Admissibility question does not arise.
- Initial regularity/decay: $u_0$ a **trigonometric polynomial** (hence $C^\omega$, finite energy, all Clay (7) conditions automatic), supported on the $N_0$-scale carrier boxes with $N_0\gg1$. Total energy $\mathcal E_{\rm tot}=4c_E/N_0$, which is small even for large $c_E$ (S3 $c_E$-collapse).
- Viscosity: $\nu>0$ fixed, treated **exactly** through the repo heat factor $h_{\nu,\tau}(r)=(1-e^{-\nu\tau r^2})/(\nu r^2)$ at each stage; never taken to $0$.

---

## B. Central mathematics

### B1. Why period must exceed 1: the harmonic channel is $O(\eta)$-null

Let a carrier occupy $N(a+\Omega_\eta)$, $\Omega_\eta=[-\eta,\eta]^3$, $\gamma_W=1$ (fixed-relative width — the *only* regime not condemned by the mesoscopic no-go). For a **single** linearly polarized mode $u=e\cos(Na\cdot x)$, $e\perp a$,
$$(u\cdot\nabla)u=-N(e\cdot a)\,e\,\sin(Na\cdot x)\cos(Na\cdot x)\equiv 0 ,\tag{LP-1}$$
because $e\cdot a=0$ exactly. For a box, $U_a(\xi)\perp(a+\xi)$ gives $U_a\cdot a=-U_a\cdot\xi=O(\eta)$, so the harmonic channel $a+a\to 2a$ has amplitude $O(\eta)$ relative to a mixed channel $a+b$, $b\neq\pm a$.

**Consequence.** A period-1 (single carrier set, $\lambda=2$) relay $\mathcal A\to2\mathcal A$ can only use harmonic channels unless $\mathcal A+\mathcal A\supset 2\mathcal A$ *by mixed pairs*. For the known relay-1 alphabet $\mathcal A_0=\{\pm p,\pm q\}$, $p=(1,1,0)$, $q=(1,0,1)$: $\mathcal A_0+\mathcal A_0=\{0,\pm2p,\pm2q,\pm(p{+}q),\pm(p{-}q)\}$ and $2p,2q$ arise **only** as harmonics. Hence period-1 doubling on this alphabet is $O(\eta)$-suppressed. **This forces $L\ge2$.**

### B2. The Butterfly cycle (composition condition $\prod\lambda_j=2$, derived)

Take the stage map on carrier *pairs*
$$\mathfrak B:\{\pm a,\pm b\}\longmapsto\{\pm(a+b),\pm(a-b)\},\qquad
\mathfrak B^2=2\cdot\mathrm{Id}.\tag{LP-2}$$
Both output carriers come from **mixed** channels, so both are $O(1)$. $\mathfrak B^2=2\,\mathrm{Id}$ is the composition condition: $L=2$, super-period ratio $\Lambda=2$, **automatically**, with no free parameter. Instantiated on relay-1's alphabet:
$$\mathcal A_0=\{\pm(1,1,0),\pm(1,0,1)\}\ \xrightarrow{\ \mathfrak B\ }\
\mathcal A_1=\{\pm C_1,\pm r\}=\{\pm(2,1,1),\pm(0,1,-1)\}\ \xrightarrow{\ \mathfrak B\ }\ 2\mathcal A_0 .\tag{LP-3}$$
Here $C_1=p+q$ is exactly relay-1's child, and $r=p-q$ is exactly **relay-2's first parent** — the two registered exact relays are the two stages of one cycle, as the lens requires. Applying the lattice rotation $R:(x,y,z)\mapsto(z,y,-x)$ (order 4, $R\in O_h$) gives $Rp=r$, $Rq=s=(1,0,-1)$ and $R\mathcal A_1=\{\pm(1,1,-2),\pm(1,-1,0)\}$, whose first element is exactly relay-2's child $C_2$. So **relay-2 $=R\cdot$relay-1**, and (6.11)'s $\mathcal Q_\varphi R$ closure is realised with a *concrete* $R$.

**Per-stage ratios.** Using the energy-weighted (rms) carrier radius $\bar N_j=(\,|\mathcal A_j|^{-1}\sum_{a\in\mathcal A_j}|a|^2)^{1/2}N$: $\bar N_0^2=2N^2$, $\bar N_1^2=(6+2)/2\,N^2=4N^2$, $\bar N_2=8N^2$. Hence
$$\lambda_0=\lambda_1=\sqrt2,\qquad \prod_{j=0}^{1}\lambda_j=2 .$$
**Crucially the lattice unit does not change between stage 0 and stage 1** — only after the super-period does $N\mapsto2N$. This is the mechanism by which an *irrational* per-stage scale ratio $\sqrt2$ is realised on $\mathbb Z^3$; it is impossible at $L=1$. The cycle is genuinely period-2 because $\mathcal A_0\cap\mathcal A_1=\varnothing$ (isotropic-radius pair vs. split-radius pair), so $U(\cdot,s)\neq U(\cdot,s+S/2)$ with strictly positive support distance.

**Per-stage critical energy ratio.** $E_{\bar N}=c_E/\bar N$ gives $E_1/E_0=1/\lambda_0=1/\sqrt2\approx0.7071$ per stage, $1/2$ per super-period (consistent with (6.9)). Within stage 1 the split between the two carriers is fixed by criticality:
$$\frac{E(r)}{E(C_1)}=\frac{|C_1|}{|r|}=\frac{\sqrt6}{\sqrt2}=\sqrt3 .\tag{LP-4a}$$

### B3. Exact polarization algebra; the both-channels witness

With $u=\alpha e_a\cos(Na\cdot x)+\beta e_b\cos(Nb\cdot x)$, $e_a\perp a$, $e_b\perp b$, the Fourier convention $\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$ gives, with $X:=(e_a\!\cdot\! b)\,e_b$, $Y:=(e_b\!\cdot\! a)\,e_a$,
$$\mathcal N_{\pm(a+b)}\propto P_{a+b}(X+Y),\qquad
\mathcal N_{\pm(a-b)}\propto P_{a-b}(Y-X).\tag{LP-4}$$
*Sign/consistency check against the repo.* For relay-1's registered polarizations $e_p=e_3$, $e_q=e_2$: $X=(e_p\!\cdot\!q)e_q=(0,1,0)$, $Y=(e_q\!\cdot\!p)e_p=(0,0,1)$, $Y-X=(0,-1,1)=-(p-q)$, so $P_{p-q}(Y-X)=0$ — **exactly reproducing the repo identity $P_{p-q}\mathcal N=0$**. And $P_c(X+Y)=(0,1,1)-\tfrac26(2,1,1)=\tfrac23(-1,1,1)\neq0$, transverse to $c$. Dimensions: $[\mathcal N]=[N][\hat u]^2$, matching $[\partial_t\hat u]$ under $t=\tau N^{-2}$, $[\hat u]\sim N^{-1}$. ✓

**The repo's exact searches gated on *difference-zero*** (216/256 orientations were discarded for "difference nonzero", L-14a). The Butterfly **inverts that gate**. An explicit rational witness with both channels alive:
$$\boxed{\ e_p=(0,0,1),\qquad e_q=(1,0,-1)\ }\quad\Longrightarrow\quad
P_c(X+Y)=\tfrac13(1,-1,-1),\qquad P_{p-q}(Y-X)=(-1,1,1).\tag{LP-5}$$
(Check: $X=(1,0,-1)$, $Y=(0,0,1)$, $X+Y=(1,0,0)$, $Y-X=(-1,0,2)$; $(-1,1,1)\cdot(0,1,-1)=0$ ✓, $(1,-1,-1)\cdot(2,1,1)=0$ ✓.)

**Tuning the critical split.** Parametrise $e_p=(0,0,1)$, $e_q=\gamma(1,0,-1)/\sqrt2+\delta(0,1,0)$, $\gamma^2+\delta^2=1$. A short computation using $\kappa:=\sqrt2\gamma+2\delta$ gives the closed forms
$$P_{p-q}(Y-X)=\tfrac{\gamma}{\sqrt2}(-1,1,1),\qquad
P_c(X+Y)=u\,(1,-1,-1),\quad u=\frac{\gamma-2\sqrt2\,\delta}{3\sqrt2}.\tag{LP-6}$$
Including the Duhamel weights $h_{\nu,\tau}(\sqrt2)$, $h_{\nu,\tau}(\sqrt6)$ (children at $|k|=\sqrt2N$, $\sqrt6N$), the **stage-closure coefficient relation** is
$$\frac{E(r)}{E(C_1)}=\frac{9\gamma^2}{(\gamma-2\sqrt2\delta)^2}\cdot
\left(\frac{h_{\nu,\tau}(\sqrt2)}{h_{\nu,\tau}(\sqrt6)}\right)^{\!2}\stackrel{!}{=}\sqrt3 .\tag{LP-7}$$
In the inviscid-stage limit $h$-ratio $\to1$ this is solved exactly by
$$\delta/\gamma=\frac{1-3^{3/4}}{2\sqrt2}=-0.45241\ldots\tag{LP-8}$$
So the critical energy allocation is **attainable**, at an irrational polarization angle — the point where rational arithmetic must hand over to interval arithmetic.

### B4. Scale stagger and the exact sumset-miss theorem

Run a second Butterfly tower (the $R$-rotated copy, i.e. relay-2's alphabet $\{\pm r,\pm s\}$) at scale $\rho N$, $\rho=P/Q\in(1,2)$ rational, $\gcd(P,Q)=1$. Put $N=QM$, so tower-I centres are $M\,Qa$ ($a\in\mathcal A^{\rm I}:=\{\pm p,\pm q,\pm C_1,\pm r\}$) and tower-II centres $M\,Pb$ ($b\in\mathcal A^{\rm II}:=R\mathcal A^{\rm I}$). Box half-widths $\eta QM$, $\eta PM$. Diagonal cross-talk outputs sit at $M(Qa\pm Pb)$ with half-width $\eta(P+Q)M$; tagged carrier boxes at all tower levels sit at $2^mMQa'$ and $2^mMPb'$, $m\in\mathbb Z$.

> **Theorem (LP-9, exact cross-talk exclusion).** If $P$ and $Q$ each possess an odd prime factor (and $\gcd(P,Q)=1$), then for all $m\in\mathbb Z$, all $a\in\mathcal A^{\rm I}$, $b\in\mathcal A^{\rm II}$, $a'\in\mathcal A^{\rm I}$, $b'\in\mathcal A^{\rm II}$,
> $$Qa\pm Pb\ \neq\ 2^mQa',\qquad Qa\pm Pb\ \neq\ 2^mPb'.$$
>
> *Proof.* Multiply by $2^{|m|}$ so all entries are integers. From $2^{|m|}(Qa\pm Pb)=Qa'$ (case $m<0$; $m\ge0$ is identical with the factor on the other side) we get $Q\mid 2^{|m|}Pb_i$ for every $i$, hence $Q\mid2^{|m|}b_i$. Writing $Q=2^\alpha Q_{\rm odd}$ forces $Q_{\rm odd}\mid b_i$ for all $i$. Every $b\in\mathcal A^{\rm II}$ has a component equal to $\pm1$ (indeed $r,s,C_2,(1,-1,0)$ all do), so $Q_{\rm odd}=1$, contradicting the hypothesis. The second family is symmetric, using that every $a\in\mathcal A^{\rm I}$ has a component $\pm1$. $\square$

Because the non-coincidences are between integer (or dyadic-rational) vectors, they upgrade to **disjoint supports** once
$$\eta<\eta_*(P,Q):=\min_{m,X,\text{centre}}\frac{\|X-2^m\!\cdot\text{centre}\|_\infty}{P+Q+2^m\max(P,Q)},\tag{LP-10}$$
a minimum over a *finite* window: for $m\ge m_+$ geometric separation $2^mQ-2(Q{+}P)>\eta(\cdot)$ takes over, and for $m\le m_-$ the shrinking wake boxes lie inside $\{\|x\|_\infty\le2^mQ(1+\eta)\}$ while the blob lies outside $\{\|x\|_\infty\ge X_{\min}-\eta(P+Q)\}$. A crude but valid bound is $\eta_*\ge1/(P+3Q)$.

**Enumeration of small $\rho$.** Require $P,Q$ coprime, each with an odd prime factor, $1<P/Q<2$ (WLOG: $\rho$ and $2/\rho$ label the same interleaved system). Minimising $P+3Q$ (which maximises the admissible $\eta_*$, hence $M^{\rm eff}\propto\eta^3$):

| $\rho$ | $(P,Q)$ | $P+3Q$ | $\eta_*\gtrsim$ |
|---|---|---|---|
| **5/3** | (5,3) | **14** | **1/14** |
| 6/5 | (6,5) | 21 | 1/21 |
| 7/5 | (7,5) | 22 | 1/22 |
| 9/5 | (9,5) | 24 | 1/24 |
| 7/6 | (7,6) | 25 | 1/25 |

$\rho=3/2$ and $\rho=3/4$ (the seed's example) have $Q=2$ or $4$, a pure power of $2$: LP-9's hypothesis fails and exclusion is no longer automatic at every $m$ (it must then be verified by finite enumeration; for $\rho=3/4$ enumeration does pass at $|m|\le2$, but there is no all-$m$ guarantee). **Selected: $\rho=5/3$, $\eta\le1/14$.**

**Honest cost accounting.** One tower with $\eta\le1/3$ has $M^{\rm eff}/N^3\le 4(2\eta)^3=32/27$; two staggered towers at $\eta=1/14$ give $2\cdot32/14^3=0.0233$ — a $50\times$ loss in effective mode count. The trade is deliberate: by S3's $c_E$-collapse **$c_E$ is free** (all channel ratios are $c_E$-independent; total energy $4c_E/N_0$ is made small by starting deep), whereas uncontrolled cross-talk is **not** compensable — the registered same-scale gadget is leakage-dominated at $\|B_{\rm low}\|^2/\|B_{\rm child}\|^2=459/106$ and $\|B_{\rm intended}\|^2/\|B_{\rm cross}\|^2=222/2483$. Therefore: **the stagger is required exactly when a single Butterfly pair cannot meet the per-stage floor $1/\sqrt2$**, and then it supplies extra carriers at exactly-zero cross-talk cost.

### B5. Sub-variants REJECTED during this derivation (kept, with failing equations)

- **R1 — "child$_1$ + $\rho\,$child$_2$ closes the cycle" (the seed's literal $L=2$ stagger closure).** Requires $C_1\pm\rho C_2\in\{2p,2q,2\rho r,2\rho s\}$. Component-wise: $C_1+\rho C_2=(2{+}\rho,1{+}\rho,1{-}2\rho)$; $=2p\Rightarrow\rho=0$; $=2q\Rightarrow\rho=0$; $=2\rho r\Rightarrow2{+}\rho=0$; $=2\rho s\Rightarrow\rho=2$ but then $1{+}\rho=3\neq0$. $C_1-\rho C_2=(2{-}\rho,1{-}\rho,1{+}2\rho)$; $=2p,2q\Rightarrow\rho=0$; $=2\rho r\Rightarrow\rho=2$ but $1{-}\rho=-1\neq2\rho=4$; $=2\rho s\Rightarrow\rho=1$ but $2{-}\rho=1\neq2\rho=2$. **No $\rho>0$ works — REJECTED.** This is what forced the Butterfly (LP-2) instead.
- **R2 — period-1 doubling $\mathcal A_0\to2\mathcal A_0$.** Fails by (LP-1): amplitude $O(\eta)$, i.e. $\ge14\times$ down in amplitude, $\ge200\times$ in energy at $\eta=1/14$. **REJECTED as a primary channel.**
- **R3 — same-scale ($\rho=1$) four-parent two-relay gadget.** $\|B_{\rm low}\|_2^2/\|B_{\rm child}\|_2^2=459/106$, child fraction $106/565$. **REJECTED (registered).**
- **R4 — translation split (S6a).** Suppression $(W|\Delta x|)^{-3}$ applies equally to the intended stage-2 pair. **REJECTED.**
- **R5 — relay-1's registered polarization $(e_3,e_2)$ inside the Butterfly.** $P_{p-q}(Y-X)=P_{p-q}(0,-1,1)=0$ exactly: the difference carrier is never populated, $\mathcal A_1$ degenerates to one carrier, and by R2 the cycle stalls. **REJECTED** — replaced by (LP-5).

---

## C. Scaling table

Clock: $\tau=T-t$; front $N=(2a\tau)^{-1/2}$ from $ds=N^2dt$, $N=N_0e^{as}$; every entry carries a multiplicative $S$-periodic factor $\Pi(\log\tau)$ of period $\log4$ (see D). $\gamma_N=1/2$; box exponent $\gamma_W=1$.

| Quantity | Exponent in $\tau$ | Status |
|---|---|---|
| Front-shell energy $E_N=c_E/N$ | $\tau^{+1/2}$ | $\to0$ |
| Total energy $\sum_j2c_E/N_j$ | $\tau^{0}$, $\le4c_E/N_{\min}$ | bounded ✓ (required) |
| Enstrophy $\sum_jN_j^2\|u_j\|_2^2\simeq2c_EN$ | $\tau^{-1/2}$ | diverges |
| Dissipation rate $\nu\!\cdot\!$enstrophy | $\nu\,\tau^{-1/2}$ | $\int_0^T\!<\infty$ ✓ |
| Global $\|u\|_{L^3}^3\simeq c\,\#\text{octaves}$ | $\tfrac{c}{2}\log(1/\tau)$ | $\to\infty$ (ESS) ✓ |
| $\|\omega\|_\infty\simeq N^2$ | $\tau^{-1}$ | $\int\|\omega\|_\infty dt=\infty$ (log) ✓ BKM |
| $\|u\|_\infty\simeq N$ | $\tau^{-1/2}$ | Type-I boundary |
| Nonlinear term $\|\mathbb PB(u,u)\|_2\simeq N\|u\|_\infty\|u\|_2$ | $\tau^{-3/4}$ | |
| Pressure term $\|\nabla p\|_2$ (Riesz, same order) | $\tau^{-3/4}$ | |
| Physical time remaining $\tau=\Theta_L\tau_*N^{-2}\!\cdot\!\tfrac43$ | $N^{-2}$ | $\Theta_2=1+\tfrac12=\tfrac32$ |
| Fourier bandwidth $\kappa N$, $\kappa=\sqrt6+\eta$ | $\tau^{-1/2}$ | |
| Active mode count $M_N=2\!\cdot\!4(2\eta N)^3$ | $\tau^{-3/2}$ | $M^{\rm eff}/N^3=0.0233$ |

Per-super-period time budget: $\Delta t=\tau_*(\bar N_0^{-2}+\bar N_1^{-2})=\Theta_2\tau_*N^{-2}$, $\Theta_L=\sum_{j<L}(\prod_{i<j}\lambda_i)^{-2}$; total $T-t_0=\tfrac43\Theta_2\tau_*N_0^{-2}=2\tau_*N_0^{-2}<\infty$.

---

## D. Closed feedback loop (every arrow a formula)

$$\underbrace{u_j\ \text{on}\ \mathcal A_j\ \text{at}\ N_j,\ E=c_E/\bar N_j}_{\text{state}}
\xrightarrow{\ \text{(i)}\ }\underbrace{\mathcal Q_r(U,U)}_{\text{(6.6)}}
\xrightarrow{\ \text{(ii)}\ }\underbrace{V=-N^{-2}h_{\nu,\tau_*}(|\cdot|/N)\mathcal Q_{\rm tot}}_{\text{Duhamel}}
\xrightarrow{\ \text{(iii)}\ }\underbrace{\mathfrak T_jU}_{\text{(6.8)}}\xrightarrow{\ \text{(iv)}\ }u_{j+1}$$

(i) $\mathcal N_{a\pm b}=-iP_{a\pm b}\big[(e_a\!\cdot\!b)e_b\pm(e_b\!\cdot\!a)e_a\big]\cdot(\text{amp})$ — **both** signs $O(1)$ by (LP-5); harmonic channels $O(\eta)$ by (LP-1); diagonal cross-channels **identically zero on every tagged box** by LP-9 + (LP-10).
(ii) Weight $h_{\nu,\tau_*}(\sqrt2)$ on the $r$-child, $h_{\nu,\tau_*}(\sqrt6)$ on the $C_1$-child; $0<h\le\tau_*$, phases preserved ($h>0$ real multiplier).
(iii) Renormalise: $\mathcal E(\mathfrak T U)/\mathcal E(U)\ge1/\lambda_j=1/\sqrt2$ is the stage closure; split enforced by (LP-7).
(iv) $\mathcal A_{j+2}=2\mathcal A_j$, $N\mapsto2N$, $\tau\mapsto\tau/4$; the loop closes after $L=2$ maps, up to $\mathcal Q_\varphi R$ with $R^4=\mathrm{Id}$.

**Log-periodic observable (the testable signature).** $U(y,s+S)=U(y,s)$, $S=\log2/a$, and $\tau=N^{-2}/(2a)$ give $s=-\tfrac1{2a}\log(2a\tau)+{\rm const}$, so with $W$ of period $2aS=2\log2$:
$$\tau\,\|\omega(t)\|_\infty=W\big(\log(1/\tau)\big),\qquad \sqrt\tau\,\|u(t)\|_\infty=A\big(\log(1/\tau)\big),\qquad \text{period }\log4 .\tag{LP-11}$$
Because $\mathcal A_0\neq\mathcal A_1$ with disjoint support, $W$ is **not** $\log2$-periodic: its Fourier spectrum in $\log(1/\tau)$ has a fundamental at $2\pi/\log4$ *and* a comparable second harmonic at $2\pi/\log2$. A period-1 cascade would show only $2\pi/\log2$. **This is the discriminating measurement.**

---

## E. Obstruction audit (collision points named)

1. **Energy bound (F-N1/N2, Leray).** Total energy $\le4c_E/N_{\min}$; we never use energy as a signature. No collision.
2. **Finite dissipation.** $\nu\int_0^T\!\text{enstrophy}\,dt\asymp\nu\int_0^T\tau^{-1/2}d\tau<\infty$. Consistent, not violated.
3. **ESS $L^\infty_tL^3_x$.** We *must* diverge. Collision point: if the wake decayed geometrically per octave, $\|u\|_3$ would stay bounded. Evasion (S4): the octave at $N_j$ is exposed to viscosity for elapsed time $O(\tau_*N_j^{-2})$, so its retention factor is $e^{-C\nu\tau_*}$, **a constant independent of $j$**; hence $\|u\|_3^3\gtrsim c\,\#\text{octaves}\to\infty$.
4. **Fixed-finite-bandwidth no-go (F-α1 / VR-L-011).** Bandwidth $\kappa N(t)\to\infty$; the trajectory is in no fixed $V_S$. No collision.
5. **Pure-swirl $L^3$ no-go (VR-L-016).** Its proof uses $\partial_\theta p_0\equiv0$ from axisymmetry. Our field is non-axisymmetric (carriers $p,q$ share no common rotation axis with matched polarizations) and lives on $\mathbb T^3$; and our signature is not $F'(0)$. Not applicable.
6. **One-scale self-similar no-go (NRS1996 / Tsai1998).** Hypothesis: $U$ **$s$-independent**. Collision point exactly here. Evasion certificate is *exact and support-based*: $\mathcal A_0\cap\mathcal A_1=\varnothing$ ($\{\pm p,\pm q\}$ vs $\{\pm C_1,\pm r\}$), so $\|U(\cdot,s)-U(\cdot,s+S/2)\|_2^2=\|U(\cdot,s)\|^2+\|U(\cdot,s+S/2)\|^2>0$ with no cancellation. Not a Leray profile.
7. **Seregin DSS-Liouville (2024/2026).** Our object is DSS with factor 2 *twisted by $R$ and $\mathcal Q_\varphi$*, on $\mathbb T^3$ where the scaling group is broken below $|k|=1$ (the tower is finite-depth-from-below, DSS only asymptotically). Seregin's hypotheses are $\mathbb R^3$ + local mixed-norm extra integrability, unverified here and preprint-status in the repo. **Flagged as a live theoretical risk, not evaded by construction.**
8. **Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}_N/N^3$.** We sit at $\gamma_W=1$, the only uncondemned regime. Collision arithmetic, $\kappa^2=(\sqrt6+\eta)^2\approx6.6$, $M^{\rm eff}/N^3=0.0233$, per-stage requirement $D\ge1/\sqrt2$:
$$c_E\ \ge\ \frac{1/\sqrt2}{2(6.6)\tau_*^2(0.0233)}=\frac{2.30}{\tau_*^2}\qquad(\tau_*=\tfrac14\Rightarrow c_E\ge36.8).$$
Legal by S3 ($c_E$ free; total energy $4c_E/N_0$). **This is a feasibility window, not a proof.**
9. **Diagonal cross-talk gate.** THE new content: exactly zero by LP-9. Named collision point: $Qa\pm Pb=2^mQa'$ requires $Q_{\rm odd}\mid b_i$, impossible since every $b\in\mathcal A^{\rm II}$ has a unit component and $Q_{\rm odd}\ge3$.
10. **CSTY2009 axisym Type-I exclusion.** We *are* Type-I ($\sqrt\tau\|u\|_\infty=A$ bounded). CSTY's hypothesis is an **axisymmetric** strong solution on $\mathbb R^3$; its proof runs De Giorgi–Nash–Moser on the swirl equation and has no known non-axisymmetric analogue. **This is the candidate's single largest exposure**: any future general-3D Type-I exclusion kills it outright. Pre-registered.
11. **KNSS ancient-solution Liouville.** The rescaling limit is the two-sided tower, non-axisymmetric, so neither KNSS(a) (no-swirl) nor KNSS(b) ($|u|\le C/r$ + axisymmetry) applies. Noted honestly: the ancient limit does satisfy $|u(x)|\asymp C/|x|$ (amplitude $N$ at scale $1/|x|$) — borderline, saved only by the axisymmetry hypothesis.
12. **Galerkin global existence.** Growing band; not applicable.
13. **Smooth-forcing high-frequency decay (F-N4).** No forcing (target (B)). Vacuous.
14. **Front-resolution threat model (TM-22, $\ge7$ points per front scale).** Front width $2\eta N\ge N/7$; $\ge7$ modes per box requires $N\ge49$ — the pilot's resolution floor.
15. **Z-01/Z-02 kinematic incompatibility.** Carriers are exact Fourier boxes, not compact-support packets: exact band-limitation holds by construction, no Littlewood–Paley leakage.
16. **VR-L-019 (bandwidth exponent $<1$).** That exponent is $\gamma_N$ (time), here $1/2\in(0,1)$ ✓ — distinct from the box exponent $\gamma_W=1$. Do not conflate.

---

## F. Minimal falsification experiment ($\le1$ h)

**P-LP1 (exact rational, minutes) — invert the difference gate.** Reuse `expanded_carrier_search.canonical_waves_in_box(2)`, `primitive_polarizations(wave, 2)`, `projected_mixed_channels`. Variables: wavevector pair $(a,b)$ with $|a|^2=|b|^2=2$, $\|k\|_\infty\le2$; primitive polarizations, component bound 2 then 3. **Gates (all in `fractions.Fraction`):** (G1) $P_{a+b}(X+Y)\neq0$ **and** $P_{a-b}(Y-X)\neq0$; (G2) both signed fluxes $\Pi_{a+b}>0$, $\Pi_{a-b}>0$; (G3) Butterfly closure — the stage-1 pair $\{a{+}b,a{-}b\}$ has both second-stage channels $\Pi_{2a},\Pi_{2b}>0$; (G4) energy split within $[\sqrt3/2,2\sqrt3]$ of (LP-4a) before tuning. Independently re-verify hits with `exact_carrier_record_verifier.verify_serialized_expanded_carrier_certificate`. **Success:** $\ge1$ orientation passing G1–G3. **Kill:** zero hits at component bound 3 ⟹ the Butterfly has no exact carrier realisation in that scope; lane closed at that scope. (LP-5 already guarantees G1; G2/G3 are genuinely open.)

**P-LP2 (exact integer, seconds) — sumset-miss enumerator.** For all $P,Q\le15$ coprime, $1<P/Q<2$: enumerate the finite window $m\in[m_-,m_+]$ and compute $\eta_*(P,Q)$ from (LP-10). **Success:** $\eta_*(5/3)>0$ and $\ge1/14$; LP-9 confirmed by brute force. **Kill:** $\eta_*=0$ for all admissible $\rho$ ⟹ theorem or bookkeeping wrong.

**P-LP3 (float, $\le1$ h) — two-tower Galerkin cross-talk null test.** $Q=3,P=5,M=1$ ($N=3$, $\rho N=5$), cutoff $\ge2\sqrt6\cdot5\approx25$, grid $64^3$, RK4$\times16$, $t=2\tau_*N^{-2}$. Reuse `leray_response_relay.leray_advection` / `leray_project`, `mesoscopic_galerkin.run_small_mesoscopic_galerkin` harness, `mesoscopic_local_fft.measure_local_fft_cloud` for tagged-channel splits. Measure: intended child energies per tower vs. energy landing on tagged boxes from diagonal pairs. **Success:** tagged-box cross-talk $\le10^{-14}$ relative (contrast: unstaggered gadget gave $5.689\times10^{-4}$ cross vs $5.747\times10^{-4}$ intended — near parity). **Kill:** $>10^{-10}$ relative.

**Arithmetic discipline.** P-LP1, P-LP2 exact rational/integer (mandatory — these are the symbolic claims). P-LP3 float is permitted (it is a null test with a $10^{4}$ margin). The polarization angle (LP-8) is irrational and must be handled by interval arithmetic at certificate stage, never by float equality.

---

## G. Proof chain to Clay (B) (10 obligations)

1. **O1.** Exact carrier witness: a polarization pair passing P-LP1 G1–G3, certified in rational arithmetic with an independent verifier. *(status: LP-5 gives G1; G2/G3 open)*
2. **O2.** Continuum profile version: existence of box profiles $U_a$ on $\Omega_\eta$ realising (LP-7) with $\lambda_j=\sqrt2$; Riemann-sum limit (6.4)–(6.6) justified.
3. **O3.** Exact cross-talk exclusion at $\rho=5/3$: LP-9 + explicit $\eta_*$ (P-LP2), promoted to a support-disjointness lemma for the full sumset.
4. **O4.** Per-stage positive flux margin net of viscosity: $q_*>0$ uniformly in $j$ — *the single largest unproven step across the whole repo lane*.
5. **O5.** Existence of an attracting $L=2$ periodic orbit of the S1 front flow $\partial_s\Psi=a(2\Psi+\xi\!\cdot\!\nabla_\xi\Psi)-\nu|\xi|^2\Psi-\mathcal Q(\Psi,\Psi)$, obtained by forward RG integration (no nonconvex optimisation), satisfying (6.11) with $\Lambda=2$, $R$ as in B2.
6. **O6.** Interval-arithmetic enclosure of O5 (radii-polynomial / Y+Z(r)<r), i.e. PO-04/PO-13 for this object.
7. **O7.** Off-chain leakage summability: $\sum_j\ell_j<\infty$ with $\ell_j$ from the *non-tagged* sumset (LP-9 removes the diagonal family exactly; harmonic and wake-sweeping remain, the latter absorbed into a Galilean/translation gauge).
8. **O8.** Finite physical time: $T-t_0=\tfrac43\Theta_2\tau_*N_0^{-2}<\infty$ with rigorous stage-rate bounds (PO-10).
9. **O9.** Norm divergence: $\limsup_{t\uparrow T}\|u(t)\|_{L^3}=\infty$ via the S4 per-octave lemma (Hölder+Bernstein+constant per-octave retention) — a limsup over a tail interval, not a fitted exponent (PO-11).
10. **O10.** Entry: a trigonometric-polynomial $u_0$ whose true solution enters the stable manifold of the O5 orbit (PO-09; currently strategy-less repo-wide). Then non-extendability past $T$ contradicts Clay (B).

**Nothing here is a theorem about NS.** O1–O3 are finite algebra and are the only parts this lens claims to have advanced; O4, O5, O10 remain open, and obstruction #10 (a hypothetical non-axisymmetric Type-I exclusion) would kill the whole construction.
