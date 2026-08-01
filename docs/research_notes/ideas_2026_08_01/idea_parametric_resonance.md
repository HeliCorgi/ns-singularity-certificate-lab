# Parametric Resonance on a Cyclic Carrier: escaping the instantaneous-margin requirement

**LENS 6** — vortex stretching / viscous diffusion phase synchronization.
**STATUS: SYMBOLIC CANDIDATE.** (Core algebra — the tree no-go, the cycle invariant
$\mathcal C$, the exact value $\mathcal C=-1239/128$, and Floquet positivity with zero
instantaneous eigenvalue margin — is computed in exact rational arithmetic with the
repo's own `fourier_torus` kernel. Nothing here is a PDE theorem.)

---

## A. Clay target

- **Statement:** negation of **CLAY-B** (global regularity, $\mathbb T^3$, no forcing);
  this implies **CLAY-D** as a corollary. Whole-space (A)/(C)/TARGET-U is *not* targeted
  (see E.8 — the $\mathbb R^3$ port collides with CSTY unless made non-axisymmetric).
- **Domain:** $\mathbb T^3$, convention $u(x)=\sum_k\hat u_ke^{ik\cdot x}$,
  $\partial_t\hat u_k=\mathcal N_k-\nu|k|^2\hat u_k$,
  $\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$, $P_k=I-k\otimes k/|k|^2$.
- **Forcing:** none ($f\equiv0$). F-N4 is therefore vacuous.
- **Initial regularity/decay:** $u_0$ a finite trigonometric polynomial (hence $C^\omega$,
  mean zero, divergence free). Finite-band *datum*, infinite-band *trajectory* — the
  VR-L-012 gap, not the VR-L-011 excluded class.
- **Viscosity:** $\nu>0$ fixed, never sent to $0$. $\nu$ enters only through the
  ratio $\nu|k|^2/\chi$, $\chi$ = Floquet exponent of the carrier network.

---

## B. Central mathematics

### B.1 The exact triad the repo already owns

For $p=(1,1,0)$, $q=(1,0,1)$, $c=p+q=(2,1,1)$, $n=p\times q=(1,-1,-1)$ and
$u=Be_3\sin(sp\cdot x)+Ce_2\cos(sq\cdot x)+Dn\cos(sc\cdot x)$, `exact_leray_relay`
gives $\Pi_c=sBCD/2$. Running the same rational kernel on all three legs
(`leray(advection(u,u))`, projected on the populated polarizations) yields the **full**
amplitude system, which the repo had not recorded:

$$\dot B=\alpha_B CD-2\nu s^2B,\quad \dot C=\alpha_C BD-2\nu s^2C,\quad
\dot D=\alpha_D BC-6\nu s^2D,$$
$$\boxed{\alpha_B=\alpha_C=-\tfrac s2,\qquad \alpha_D=+\tfrac s3},\qquad
\alpha_B+\alpha_C+3\alpha_D=0 \text{ exactly.}$$

(Verified for $s=1,2,3,4$; the last identity is detailed energy conservation, since
$E_B=B^2/4,\ E_C=C^2/4,\ E_D=3D^2/4$. Sign check: $\alpha_D>0$, so
$\Pi_c=\tfrac32 D\dot D|_{\rm nl}=\tfrac s2BCD>0$ — energy climbs, as advertised.
Dimension check: $[\alpha]=[k]$, so $\alpha\propto s$ ✓.)

**Energy coordinates.** $X_B=B/\sqrt2$, $X_C=C/\sqrt2$, $X_D=D\sqrt{3/2}$ give
$E_i=X_i^2/2$ and $\dot X_i=\beta_iX_jX_k$ with
$(\beta_B,\beta_C,\beta_D)=\sqrt{2/3}\,s\,(-\tfrac12,-\tfrac12,1)$, $\sum\beta_i=0$.

### B.2 Freezing the parent: the doublet is a damped rotation (Adler form)

Freeze $B$, linearize in $(X_C,X_D)$. Writing $\lambda_C=\nu|sq|^2=2\nu s^2$,
$\lambda_D=\nu|sc|^2=6\nu s^2$, $\bar\lambda=4\nu s^2$, $\Delta=(\lambda_D-\lambda_C)/2=2\nu s^2$:

$$\frac{d}{dt}\binom{X_C}{X_D}=\begin{pmatrix}-\lambda_C & -\tfrac{s}{2\sqrt3}B\\
\tfrac{s}{\sqrt3}B & -\lambda_D\end{pmatrix}\binom{X_C}{X_D}.$$

The off-diagonal **product is negative**, $-s^2B^2/6<0$, so the block is *antisymmetrisable*:
with $X'=X_C 2^{1/4},\ Y'=X_D2^{-1/4}$ the coupling becomes exactly antisymmetric with
rate $\omega=sB/\sqrt6$, and in polar coordinates $(\|v\|_*,\varphi)$

$$\boxed{\frac{d}{dt}\log\|v\|_*=-\bar\lambda+\Delta\cos2\varphi,\qquad
\dot\varphi=\omega(t)-\Delta\sin2\varphi.}$$

This is an Adler phase-synchronisation equation: the *only* growth resource is parking the
phase at $\varphi=0$ (all energy in the least-damped leg), and even there
$-\bar\lambda+\Delta=-\lambda_C<0$.

> **Proposition 1 (tree no-go — new, exact).** Let the child band obey
> $\dot x=(K(t)-\Lambda)x$, $\Lambda=\mathrm{diag}(\nu|k_i|^2)$, $K$ zero-diagonal with
> $K_{ij}(t)=\beta_{ij}X_{P(ij)}(t)$, $\beta_{ij}\beta_{ji}<0$ (cascade edges: parent is the
> lowest leg). If the coupling graph is a **forest**, a constant positive diagonal metric
> $\Theta$ makes $K(t)$ antisymmetric for **all** $t$, whence
> $\|x(t)\|_\Theta\le e^{-\min_i\nu|k_i|^2\,t}\|x(0)\|_\Theta$.
> **No modulation — any amplitude, waveform, period — produces net growth.**
> For the repo triad: $E(t)\le 2e^{-2\lambda_Ct}E(0)$, transient overshoot $\le\!\times2$.

*Proof sketch.* On a forest one may solve $\theta_i/\theta_j=-\beta_{ij}/\beta_{ji}>0$
edge-by-edge without consistency conditions; $\Theta^{1/2}K\Theta^{-1/2}$ is then
antisymmetric and $\tfrac{d}{dt}\|x\|_\Theta^2=-2x^{\!\top}\Theta\Lambda x$. $\square$

### B.3 REJECTED sub-variants (kept, with the failing equation)

- **REJECTED — scalar parametric ansatz.** $\dot D=\sigma(t)D-\nu|k|^2D$ with
  $\sigma=\sigma_0+\sigma_1\cos\Omega t$. Floquet exponent of a *scalar* linear ODE is the
  time average: $\mu=\sigma_0-\nu|k|^2$. **Gain exactly zero.** Parametric resonance needs
  $\ge2$ dimensions.
- **REJECTED — Mathieu on the doublet.** Eliminating $X_D$ and removing damping
  ($\xi=B^{1/2}w$) gives the Hill equation
  $\ddot w+[(\omega_0^2-\Delta^2)+\varepsilon((2\omega_0^2-\tfrac{\Omega^2}{2})\cos\Omega t-\Delta\Omega\sin\Omega t)]w=0$,
  $\omega_0=sB_0/\sqrt6$, $B=B_0(1+\varepsilon\cos\Omega t)$. At principal resonance
  $\Omega=2\omega_{\rm eff}=2\sqrt{\omega_0^2-\Delta^2}$ the pump amplitude collapses to
  $\tilde\varepsilon=2\varepsilon\Delta\omega_0$, so
  $\mu_w=\tilde\varepsilon/4\omega_{\rm eff}=\varepsilon\Delta\omega_0/2\omega_{\rm eff}\le\varepsilon\Delta/2$
  in the perturbative regime $\omega_0\gg\Delta$. Net growth needs
  $\mu_w>\bar\lambda=2\Delta$, i.e. $\boxed{\varepsilon>4}$ — outside Mathieu validity, and
  killed outright by Proposition 1 at **any** $\varepsilon$. *Failing inequality:*
  $\varepsilon\Delta/2<2\Delta$ for all $\varepsilon\le4$.
- **REJECTED — harvesting non-normal transient growth.** The symmetric part of the doublet
  has $\lambda_{\max}=-\bar\lambda+\sqrt{\Delta^2+s^2B^2/48}>0$ iff $|B|>24\nu s$
  (numerically confirmed: $\lambda_{\max}^{\rm sym}=0.0000$ at $B=24$, $s=\nu=1$). But the
  *eigenvalues* stay at $\mathrm{Re}=-\bar\lambda$ for every $B$, and Proposition 1 caps the
  Floquet exponent at $-\lambda_C$. Transients do not compound.

### B.4 The escape: cycles, not chains

Proposition 1 fails exactly when the coupling graph carries a **cycle**. For a 3-cycle
$1\to2\to3\to1$ mediated by parents $p_a,p_b,p_c$ with $p_a+p_b+p_c=0$, define the edge
ratios $r_e=\beta_{\rm bwd}/\beta_{\rm fwd}$ and the

$$\textbf{cycle invariant}\qquad \mathcal C:=r_ar_br_c,\qquad
\det K=\beta_a\beta_b\beta_c\,(1+\mathcal C).$$

Antisymmetrisability $\iff\mathcal C=-1$ $\iff\det K=0$. Since $\mathrm{tr}\,K=0$, the
eigenvalues $\{\mu_r,\rho\pm i\sigma\}$ satisfy $\mu_r+2\rho=0$ and
$\mu_r(\rho^2+\sigma^2)=\det K$, so

> **Proposition 2.** $\mathcal C\ne-1\Rightarrow K$ has an eigenvalue with
> $\mathrm{Re}>0$ (namely $\max(\mu_r,-\mu_r/2)>0$), for **either** sign of $\det K$.

### B.5 A concrete integer cyclic carrier (exact)

$$k_1=(2,2,1),\ k_2=(2,1,2),\ k_3=(1,2,2)\quad (|k|^2=9),$$
$$p_a=(0,-1,1),\ p_b=(-1,1,0),\ p_c=(1,0,-1)\quad (|p|^2=2),\ p_a+p_b+p_c=0,$$
$k_1+p_a=k_2$, $k_2+p_b=k_3$, $k_3+p_c=k_1$. Children on **one** shell; parents form a
**closed** triad. Exhaustive exact search over the $2^3\!\times\!2^3\!\times\!2^3\!\times\!2^3=512$
polarization/phase orientations (rational kernel, validated against the known triad —
recovers $r=-1/2$ exactly):

| class | count |
|---|---|
| all six couplings nonzero | 512 |
| $\mathcal C=-1$ exactly (**provably dead**, Prop. 1) | 128 |
| $\mathcal C\ne-1$ (parametrically live) | 384 |
| live **and** every edge individually stable ($r_e<0$: **zero instantaneous margin**) | **192** |

Best of the marginal class (child pols $(-5,4,2),(0,2,-1),(-8,2,2)$; parent pols
$(0,1,1),(-1,-1,0),(0,-1,0)$; children $\cos$, parents $\sin$):

$$r=\Big(-\tfrac94,\,-\tfrac73,\,-\tfrac{59}{32}\Big),\qquad
\boxed{\mathcal C=-\tfrac{1239}{128}=-9.6797},\qquad 1+\mathcal C=-\tfrac{1111}{128}.$$

Exact forward/backward Leray coefficients $-\tfrac85,-\tfrac1{24},\tfrac{32}{45}$ /
$\tfrac25,\tfrac75,-\tfrac{59}{72}$. In energy coordinates at $X_P=1$, $s=1$:
$\mathrm{tr}\,K=0$, $\det K=-0.58192$, spectrum
$\{0.149013\pm1.389381i,\,-0.298027\}$, so $\sigma_{\rm cycle}=0.149013\,sX_P$.

### B.6 The genuinely oscillatory version (phase-staggered pump)

Run the three parents **sequentially**: $p_a$ on for $[0,T_p/3)$, then $p_b$, then $p_c$.
Then at every instant $K(t)$ is a **single** 2-node block with $r_e<0$, i.e. an *elliptic*
generator: its eigenvalues are $\pm0.8i$, $\pm0.2415i$, $\pm1.0796i$ — **$\max\mathrm{Re}=0$
exactly, at all times**. Yet the monodromy $M=e^{G_c\delta}e^{G_b\delta}e^{G_a\delta}$ has

$$\rho(M)=1.1100>1,\qquad \chi_{\rm par}=\frac{\log\rho(M)}{3\delta}=0.052773\,sX_P
\quad(\delta^*=0.6593/sX_P),$$

i.e. **net growth per period with the instantaneous relay margin identically zero**
($\chi_{\rm par}/\sigma_{\rm cycle}=0.354$). The growth is pure holonomy: three elliptic
rotations with incommensurate eccentricities ($\prod\rho_e\ne1\iff\mathcal C\ne-1$)
compose to a hyperbolic map. Null control: at $\mathcal C=-1$ the same computation gives
$\rho=1.0000000$ and $\chi_{\rm par}=0$ to machine precision.

**Net-growth condition (physical):**
$$\boxed{\ \chi_{\rm par}\,sX_P>\nu|k|^2\ }\quad\Longleftrightarrow\quad X_P>170.5\,\nu s
\ \ \text{(this gadget)},$$
versus $X_P>60.4\nu s$ for the steady cycle, versus **impossible at any amplitude** for
every tree (Prop. 1). *That factor-$\infty$-to-finite jump is the whole idea.*

### B.7 Front-flow consequence (seed S1 becomes a periodic orbit)

At capacity, $\chi_{\rm par}sX_P\asymp\hat\chi\,\kappa\sqrt{2c_EM_{\rm eff}N}$ and the pump
period is set by $\omega_{\rm edge}\asymp\nu N^2$:
$$T_{\rm pump}=\frac{c_p}{\nu N^2}\ \ (c_p=O(1))\ \Longrightarrow\ S_{\rm pump}=N^2T_{\rm pump}=\frac{c_p}{\nu}=O(1)\ \text{in the parabolic clock,}$$
exactly the regime the lens specifies. The RG doubling period is $S_{\rm dbl}=\log2/a_+$
(S1). Closure of the orbit demands **commensurability**:
$$\boxed{\ a_+=\frac{\nu\log 2}{c_p\,m},\quad m\in\mathbb N\ }$$
— the front speed is **quantised**. S1's fixed point $\Psi_*$ is replaced by
$\Psi(\xi,s+S_{\rm pump})=\Psi(\xi,s)$ with $S_{\rm dbl}=mS_{\rm pump}$. The instantaneous
margin $q(t)=\Pi_c(t)/(2\nu|c|^2E_{\rm child}(t))-1$ oscillates in sign $6$ times per pump
period (once per elliptic quarter-turn per active edge); only
$\langle\Pi_c\rangle/\langle2\nu|c|^2E\rangle=1+\chi_{\rm par}/(\nu|k|^2)>1$ holds.

---

## C. Scaling table

Clock $\tau=T-t$; $N=(2a_+\tau)^{-1/2}$ ($\gamma=1/2$, Type-I boundary, matching S1/STATUS);
shell law $E_j=c_E/\lambda_j$ ($\beta=-1,\sigma=\gamma$, the post-erratum branch);
$M_{\rm eff}=\theta N^3$, $\theta=O(1)$.

| quantity | law in $N$ | law in $\tau$ |
|---|---|---|
| energy $E$ | $\le 2c_E/N_0$ | $\tau^{0}$ (bounded ✓ F-N1) |
| enstrophy $\|\omega\|_2^2$ | $\asymp c_EN$ | $\tau^{-1/2}$ |
| global $\|u\|_{L^3}^3$ | $\asymp e^{-\pi m}\log N$ | $\asymp\tfrac12 e^{-\pi m}|\log\tau|$ (**$\to\infty$**) |
| $\|\omega\|_\infty$ | $\sqrt{2\theta c_E}\,N^2$ | $\tau^{-1}$ |
| $\int_0^t\|\omega\|_\infty$ | $\asymp\log N$ | $|\log\tau|\to\infty$ ✓ BKM |
| dissipation rate $\nu\|\nabla u\|_2^2$ | $\asymp\nu c_EN$ | $\tau^{-1/2}$; $\int_0^T<\infty$ ✓ |
| nonlinear $\|\mathbb PB(u,u)\|_2$ | $2c_E\sqrt\theta N^{3/2}$ | $\tau^{-3/4}$ |
| pressure $\|\nabla p\|_2$ | same order $N^{3/2}$ | $\tau^{-3/4}$ |
| physical time remaining | $N^{-2}/2a_+$ | $\tau$ ✓ finite |
| Fourier bandwidth $\kappa N$ | $\kappa N$ | $\tau^{-1/2}$ |
| active mode count $M_{\rm eff}$ | $\theta N^3$ | $\tau^{-3/2}$ |
| pump period $T_{\rm pump}$ | $c_p/\nu N^2$ | $\propto\tau$ (∞ periods before $T$) |
| $\sqrt\tau\|u\|_\infty$ | $\sqrt{2\theta c_E/2a_+}$ | $\tau^0$ (**Type-I marginal**) |

Consistency: the BKM requirement $\|\omega\|_\infty\gtrsim N^2$ forces
$M_{\rm eff}\gtrsim N^3/2c_E$ — the *same* $N^3$ floor as L-02/L-11. Nothing is free.

---

## D. Closed feedback loop

$$
\underbrace{\text{parent triad }p_a,p_b,p_c \text{ at }N,\ \sum p_e=0}_{X_P=\sqrt{2c_E/3N}}
\xrightarrow{\ \text{stagger: }p_e\text{ on for }T_{\rm pump}/3\ }
\underbrace{K_e(t)\ \text{elliptic},\ \max\mathrm{Re}\,\lambda(K_e)=0}_{\text{no instantaneous margin}}
$$
$$
\xrightarrow{\ \mathcal C=r_ar_br_c\ne-1\ }
\underbrace{\rho(M)=\rho(e^{G_c\delta}e^{G_b\delta}e^{G_a\delta})>1}_{\chi_{\rm par}=\log\rho/3\delta}
\xrightarrow{\ \chi_{\rm par}sX_P>\nu N^2\ }
\underbrace{E_{\rm child}(nT_{\rm pump})=e^{2(\chi_{\rm par}sX_P-\nu N^2)nT_{\rm pump}}E_{\rm child}(0)}_{\text{multiplicative, not additive}}
$$
$$
\xrightarrow{\ E_{\rm child}\uparrow c_E/2N\ }
\underbrace{\text{child triangle at }|k|=3s\text{ re-tagged as parent triad at }N'=2N}_{\text{doubling pullback, }S_{\rm dbl}=mS_{\rm pump}}
\xrightarrow{\ a_+=\dot N/N^3\ }
\underbrace{N^{-2}(t)=N_0^{-2}-2a_+t}_{T=N_0^{-2}/2a_+<\infty}
$$
$$
\xrightarrow{\ \text{S4 wake persistence, retained fraction }e^{-\nu/2a_+}=e^{-\pi m}\ }
\underbrace{\|u\|_{L^3}^3\gtrsim e^{-\pi m}\log N\to\infty}_{\text{ESS-required divergence}}\ \longrightarrow\ \text{breakdown at }T.
$$
Every arrow is a displayed formula above. The loop closes on itself: the child triangle at
step $j$ is the parent triad at step $j+1$, because $\{k_1,k_2,k_3\}$ has
$k_i-k_j\in\{\pm p_e\}$ — the triangle is simultaneously a closed sum-zero difference set.

---

## E. Obstruction audit

1. **Energy bound / finite dissipation (F-N1/N2, Leray).** Table C: $E=\tau^0$,
   $\int_0^T\nu\|\nabla u\|^2\asymp\int\tau^{-1/2}d\tau<\infty$. We never use these as the
   signature. **Not violated; not evaded — respected.**
2. **ESS $L^3$ endpoint (U-X1, VR-N-002).** Requires $\limsup\|u\|_3=\infty$. Delivered by
   S4 per-octave persistence with retained fraction $e^{-\nu/2a_+}=e^{-\pi m}$; $m=1$ is
   forced (larger $m$ costs $e^{-\pi m}$ in the prefactor). **Collision point:** if retired
   octaves decay faster than $e^{-\nu N_j^2(T-t_j)}$ the sum converges and ESS kills us —
   this is pilot kill-condition K3.
3. **Fixed-finite-bandwidth no-go (F-α1, F-6/F-7, VR-L-011, `galerkin_not_tendsto_atTop`).**
   Hypothesis: trajectory confined to a **fixed** mode set $S$. **Collision point:** our
   trajectory exits $\{|k|\le\Lambda\}$ at $t_\Lambda=T-1/(2a_+\Lambda^2)<T$ for every
   $\Lambda$. Datum is finite-band (VR-L-012 gap), trajectory is not. Evaded by hypothesis.
4. **Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM_{\rm eff}/N^3$.**
   This is an exact phase-independent inequality; we do **not** contradict it. Two exact
   collision points: (i) *width* — we use only fixed-relative width $W=\eta N$,
   $\eta\in(0,1/3)$, so $M_{\rm eff}=\theta N^3$ and the bound reads $D_N\le2\kappa^2\tau^2c_E\theta$,
   $N$-independent; with $\kappa=\tfrac13,\tau=\tfrac14,\theta=\tfrac{4\pi}{3}\kappa^3f$ this needs
   $c_E\ge232/f$ — matching the repo's own measured "$c_E\approx228$ at $\eta=0.20$" row, and
   legal by S3 ($\sum_jc_E/N_j=2c_E/N_0$ small for $N_0$ large). (ii) *frame* — the L-11
   derivation sets $v(0)=0$ and freezes $f_N$, bounding the **additive** one-step Duhamel
   response. The cyclic mechanism is **multiplicative**: $v(0)\ne0$ and $f$ depends linearly
   on $v$, so the governing quantity is the Floquet exponent $\chi_{\rm par}sX_P$ versus
   $\nu N^2$, not $D_N$ versus $1/2$. $D_N$ is simply not the right functional for a
   Floquet relay. *This is the single load-bearing evasion of the document.*
5. **Diagonal cross-talk gate (leakage $2483/222\approx9{:}1$).** The gate fires because two
   parent pairs are simultaneously active, producing $A_1{+}B_2$ on the target shell.
   **Collision point:** the staggered pump has $X_{p_a}(t)X_{p_b}(t)\equiv0$ pointwise in
   $t$ — the diagonal channel is **exactly zero**, not merely suppressed. Unlike S6(a)
   (translation split) this does not kill the intended pair, because the intended cycle is
   itself sequential. Unlike S6(b) it needs no scale-stagger or sumset-miss integer program.
6. **Pure-swirl $L^3$ no-go (VR-L-016).** Hypothesis $u_0=u^\theta(r,z)e_\theta$. Our
   polarizations are $(-5,4,2),(0,2,-1),(-8,2,2)$ on $\mathbb T^3$ — not azimuthal, pressure
   channel $P\not\equiv0$. Not applicable.
7. **One-scale self-similar (NRS1996 / Tsai1998).** Both require *exact continuous* backward
   self-similarity, $u=(2a\tau)^{-1/2}U(x/\sqrt{2a\tau})$ with $U$ $s$-independent.
   **Collision point:** our $\Psi(\xi,s)$ is $s$-**periodic** with period $S_{\rm pump}=c_p/\nu>0$
   and $\partial_s\Psi\not\equiv0$; invariance holds only under the discrete subgroup
   $\lambda=2^{1/m}$. NRS/Tsai hypotheses fail. **Open risk (not evaded):** Seregin2024
   Liouville results for *discretely* self-similar solutions — flagged in `future_search` §3.2
   — must be audited before promotion. Registered as PO-audit item, not claimed evaded.
8. **CSTY Type-I exclusion.** Hypotheses: $v$ **axisymmetric**, on $\mathbb R^3$.
   **Collision point:** our object is on $\mathbb T^3$ and non-axisymmetric. But table C gives
   $\sqrt\tau\|u\|_\infty\to\sqrt{c_E/a_+}$ — **Type-I marginal**. Therefore any future
   $\mathbb R^3$ port must be non-axisymmetric or it dies immediately. Registered.
9. **KNSS ancient Liouville.** Both branches assume axisymmetry (no-swirl $\Rightarrow$
   axis-constant; with-swirl $+|u|\le C/r\Rightarrow0$). **Collision point:** our rescaled
   limit is a non-axisymmetric, genuinely $s$-periodic (hence non-stationary) bounded
   ancient object. Neither branch applies. Non-axisymmetric extensions are open.
10. **Galerkin global existence.** Same as (3).
11. **Smooth-forcing high-frequency decay (F-N4).** $f\equiv0$. Vacuous.
12. **Front-resolution threat model (TM-22, ≥7 points per front; TM-20 integrator
    amplification; TM-04 spectral pile-up).** Pilot mandates ≥7 grid points per front scale,
    ≥7 diagnostic samples per pump period (no stride thinning, TM-21), RK4 with ≥16 substeps
    per period and a doubling check (TM-20), $3N$ dealiasing (TM-03).
13. **NEW threat TM-P1 (mine): transient non-normality mistaken for Floquet growth.** The
    doublet already overshoots by $\times2$ in energy with zero Floquet gain. Kill rule:
    report $\rho(M)$ of the **monodromy over $\ge3$ full periods**, not $\|M\|$, and require
    the per-period ratio to be constant to $5\%$.
14. **NEW no-go (mine): Proposition 1.** Every chain/tree relay in the repo's entire history
    is now excluded for *all* modulations, not merely for steady ones. This retro-explains
    the $\gamma<1$ failures, the "no strict hit" of L-12, and the absence of any positive
    $q_*$: the searches were structurally acyclic.

---

## F. Minimal falsification experiment (≤1 hour)

**Variables.** Child shell $|k|^2$, parent shell $|p|^2$, 512 orientations
(polarization index $\times$ phase), $\nu$, $c_E$, stagger duty cycle $\delta$, grid $G$.

**P1 — exact rational (≤10 min).** Reuse `fourier_torus.{TrigVector,advection,leray}` and the
enumeration style of `exact_carrier_search.canonical_waves_in_box` /
`primitive_polarizations`. For every closed child triangle on a single shell with
$\sum_ep_e=0$: compute exact $r_a,r_b,r_c$, $\mathcal C$, $\det K$ as `Fraction`s.
*Must be exact rational.* Independent re-verification via
`exact_carrier_record_verifier.verify_serialized_strict_orientation_records` (different code
path: complex conjugate-symmetric convolution).
**Success S1:** $\ge1$ orientation with all $r_e<0$ (zero instantaneous margin) and
$|1+\mathcal C|\ge1/2$. *Already met in miniature:* 192/512, best $\mathcal C=-1239/128$.

**P2 — monodromy (≤10 min, float OK, exact-seeded).** Build $G_a,G_b,G_c$ in energy
coordinates; RK4 the staggered period; report $\rho(M^n)^{1/n}$ for $n=1,2,3$.
**Null control (mandatory):** an orientation with exact $\mathcal C=-1$ must give
$\rho=1\pm10^{-12}$.
**Success S2:** $\chi_{\rm par}>0$, agreeing with the exact-algebra prediction $0.0528\,sX_P$
to 3 digits, and $\rho(M^3)^{1/3}=\rho(M)$ to $5\%$ (kills TM-P1).

**P3 — full Galerkin cross-check (≤30 min).** Reuse the pattern of
`mesoscopic_galerkin.run_small_mesoscopic_galerkin` and
`leray_response_relay.{leray_project,leray_advection}` (dealiased $\mathbb P((u\cdot\nabla)u)$).
Seed the 6-mode field ($s=1$, $N=3$), $\nu$ chosen so $\nu|k|^2=\tfrac12\chi_{\rm par}sX_P$,
integrate $3T_{\rm pump}$ on $G=48^3,64^3,96^3$ with $16$ and $32$ RK4 substeps.
**Success S3:** child-triangle energy grows geometrically over $\ge3$ periods, ratio stable
to $5\%$ across all six (grid, substep) combinations, and the $\mathcal C=-1$ control decays.

**Kill conditions (pre-registered).**
- **K1:** $\rho(M)\le1+10^{-6}$ for every all-$r_e<0$ orientation $\Rightarrow$ the marginal
  class is secretly antisymmetrisable; mechanism dead.
- **K2:** the exact search over the fixed-relative-width family $W=\eta N$ shows the
  shape-averaged $\mathcal C\to-1$ as $M_{\rm eff}/N^3\to\theta>0$ (cycle invariant degenerates
  in the capacity-saturating limit) $\Rightarrow$ **decisive kill** — the mechanism would
  exist only at fixed cardinality, which `fixed_cardinality_scaling` already rejects
  ($N^{-3/2}$).
- **K3:** the retained wake fraction per octave is not bounded below by a constant
  $\Rightarrow$ $L^3$ stays bounded $\Rightarrow$ ESS kills it.
- **K4:** P3 growth disappears or shifts $>5\%$ under refinement $\Rightarrow$ TM-01/TM-20.

**Arithmetic split.** Exact `Fraction`: $r_e$, $\mathcal C$, $\det K$, all Leray coefficients,
the $\mathcal C=-1$ null control. Float permitted: monodromy eigenvalues, Galerkin evolution,
energy time series. Interval arithmetic: deferred to PO-13, not this pilot.

---

## G. Proof chain (10 obligations)

1. **Prop. 1 (tree no-go)** — finite linear algebra + Grönwall; Lean-able today
   (`GalerkinNoBlowup`-style). *Independent value even if the rest dies.*
2. **Prop. 2 + cycle invariant** — $\det K=\beta_a\beta_b\beta_c(1+\mathcal C)$, exact rational
   for the 512-orientation table; Lean-able as finite rational algebra.
3. **Floquet positivity with zero instantaneous eigenvalue margin** — interval enclosure of
   $\rho(M)>1$ for one exhibited orientation (PO-13 machinery, VR-L-015 comparison lemmas).
4. **Fixed-relative-width capacity** — construct an $\eta$-width cyclic carrier cloud with
   $M_{\rm eff}=\theta N^3$ and shape-averaged $|1+\mathcal C|\ge c>0$ uniformly in $N$.
   *(This is the hard step; K2 is its falsifier.)*
5. **Stage budget** — $\chi_{\rm par}(N)\,sX_P(N)-\nu N^2\ge q_*\nu N^2$ with $q_*>0$
   uniform in $N$, given $c_E\ge c_E^{\min}(\nu,\eta,\theta)$ (S3 collapse).
6. **Doubling pullback / recurrence** — the amplified child triangle is, after
   renormalisation, the parent triad one octave up; commensurability $S_{\rm dbl}=mS_{\rm pump}$,
   $m=1$.
7. **Periodic orbit of the S1 front flow** — attracting $s$-periodic $\Psi$ found by forward
   RG integration (S1 unlock: no nonconvex optimisation), with Floquet multipliers separated
   from gauge/truncation spectrum (`future_search` §9).
8. **Finite time** — $T-t=N^{-2}/2a_+$, $a_+=\nu\log2/c_p$; $\int^\infty N^{-2}ds<\infty$.
9. **Critical-norm divergence** — S4 per-octave lemma with retained fraction $e^{-\pi}$ gives
   $\|u\|_3^3\gtrsim e^{-\pi}\log N\to\infty$ and $\int\|\omega\|_\infty dt=\infty$.
   Discharges PO-11 without a separate obligation.
10. **Entry (PO-09) + non-extendability** — connect explicit trigonometric $u_0$ to the
    orbit's stable manifold; then BKM/ESS $\Rightarrow$ no extension past $T$ $\Rightarrow$
    negation of CLAY-B. **PO-09 remains strategy-less repo-wide; this candidate does not
    solve it.**

Obligations 1–3 are reachable now. Obligation 4 is the make-or-break. Obligation 10 is the
same wall every candidate in this repo faces.

---

**Permitted claim ceiling:** "a cyclic carrier network exhibits, in exact rational Fourier
algebra, positive Floquet gain with identically zero instantaneous eigenvalue margin, which
no tree can do." Nothing here is evidence of a finite-time singularity.
