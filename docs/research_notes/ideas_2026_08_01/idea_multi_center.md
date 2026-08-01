# LENS 2 — Nested Skew-Dyad Ladder (NSDL)

**Status: FORMAL ANSATZ** (core algebra at SYMBOLIC CANDIDATE: §B.1, §B.4, §B.6, §E.9 are exact finite computations; no orbit constructed, no numerics run.)

One-line thesis: the multi-center object that survives every registered no-go is **not** a pair of far-separated cores (that is provably unfixable, §B.6), but a **bounded-multiplicity skew dyad nested inside itself across scales**, with a per-generation rotation $\mathsf R_{\theta_0}$ as the phase code, lacunary Fourier bands $\rho<2$ as the exact cross-talk gate, and **scale-independent circulation** as the geometric meaning of the repo's critical normalization $E_N=c_E/N$.

---

## A. Clay target

- **Target: CLAY-D** (breakdown, $\mathbb T^3$) **with $f\equiv 0$**, hence simultaneously a negation of **CLAY-B**. Forcing: none — so F-N1/F-N2/F-N4 are structurally irrelevant, not merely evaded.
- Domain $\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3$, Fourier convention $u=\sum_k\hat u_ke^{ik\cdot x}$, $\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$, $P_k=I-k\otimes k/|k|^2$.
- Initial data: $u_0\in C^\infty(\mathbb T^3;\mathbb R^3)$, real, mean-zero, $\nabla\cdot u_0=0$, band-limited to $|k|\le \kappa\lambda_0$ (a trigonometric polynomial — legal by VR-L-012: fixed-band *datum*, infinite-band trajectory).
- Viscosity: $\nu>0$ fixed, **never** sent to $0$. Viscosity enters only through the single scale-independent dimensionless group $\mathrm{Re}_\Gamma=\Gamma/\nu$, $\Gamma=\sqrt{c_E/n}$.
- $\mathbb R^3$ transfer (TARGET-U / CLAY-C) is a *separate* obligation: the whole-space version needs the HS-5 bookkeeping and is explicitly **not** claimed here.

---

## B. Central mathematics

### B.1 Mutual axial strain of a skew dyad (exact)

A line vortex of circulation $\Gamma$, direction $t$, generates $u(x)=\frac{\Gamma}{2\pi\rho}\,\hat s$, $r=(x-y)-((x-y)\cdot t)t$, $\rho=|r|$, $\hat s=t\times\hat r$. Differentiating ($r_b=P^\perp_{bc}x_c$, $\partial_j\rho^2=2r_j$):
$$\partial_ju_i=\frac{\Gamma}{2\pi\rho^2}\big[\epsilon_{iaj}t_a-2\hat s_i\hat r_j\big].$$
The two antisymmetric contributions cancel identically, so $\nabla u=S$ is symmetric (as it must be: $\omega=0$ off the core), and
$$\boxed{\;S=-\frac{\Gamma}{2\pi\rho^2}\big(\hat r\otimes\hat s+\hat s\otimes\hat r\big)\;}\qquad\text{eigenvalues }\Big(+\tfrac{\Gamma}{2\pi\rho^2},-\tfrac{\Gamma}{2\pi\rho^2},0\Big).$$
Dimension check: $[\Gamma]/[L^2]=T^{-1}$ ✓. Trace $=0$ ✓ (incompressible). Neutral direction $=t$ ✓.

**Symmetric skew pair.** Impose the $\pi$-rotation $R:(x_1,x_2,x_3)\mapsto(-x_1,-x_2,x_3)$. Put tube 1 at $\tfrac d2e_1$ with $t_1=(sc_\varphi,ss_\varphi,c)$ ($s=\sin\theta$, $c=\cos\theta$), tube 2 at $-\tfrac d2e_1$ with $t_2=Rt_1$, equal $\Gamma$. Then $\sigma_2=\sigma_1$ automatically. With $P:=1-s^2c_\varphi^2$, $\rho^2=d^2P$, one computes
$$t_1\cdot\hat r=\frac{2sc^2c_\varphi}{\sqrt P},\qquad t_1\cdot\hat s=(t_1\times t_2)\cdot\hat r=\frac{2csz_\varphi}{\sqrt P}\Big|_{z_\varphi=s_\varphi},$$
$$\Rightarrow\quad \sigma_1=t_1\!\cdot\! S^{(2)}t_1=-\frac{\Gamma}{2\pi\rho^2}\,2(t_1\!\cdot\!\hat r)(t_1\!\cdot\!\hat s)
=\boxed{\;-\frac{2\Gamma}{\pi d^2}\cdot\frac{\sin^2\theta\,\cos^3\theta\,\sin2\varphi}{(1-\sin^2\theta\cos^2\varphi)^2}\;}$$

Sanity limits: **parallel** ($\theta=0$) $\Rightarrow\sigma=0$; **antiparallel** ($\theta=\pi/2$, $t_2=-t_1$) $\Rightarrow\sigma\propto\cos^3\theta=0$. *Straight parallel and straight antiparallel pairs exert exactly zero mutual axial strain* — this is why classical antiparallel-pair scenarios must rely on curvature/self-induction. Maximum at $\varphi_*=-\pi/4$, $\sin^2\theta_*\approx0.540$ ($\theta_*\approx47.3^\circ$):
$$\sigma_*=c_*\frac{\Gamma}{d^2},\qquad c_*=\tfrac2\pi\max_x\frac{x(1-x)^{3/2}}{(1-x/2)^2}=\tfrac2\pi(0.3162)=\mathbf{0.2013}.$$
At the optimum $t_1\cdot t_2=\cos2\theta_*=-0.08$: the tubes are **near-orthogonal**. Consequences used later: (i) mutual stretching is $O(1)$; (ii) viscous circulation annihilation $\propto-(\omega_1\cdot\omega_2)\propto-\cos2\theta_*=+0.08$ is *near-zero and of the non-cancelling sign* — the antiparallel-pair killer is structurally absent.

**Approach rate.** $\hat s\cdot e_1=0$ identically for this family ⟹ the pair does **not** approach by mutual induction and $\dot d=0$ at line-vortex order. This kills the naive collapse (see REJECTED M-2) and forces the nesting construction of §B.3.

### B.2 Critical normalization $\Leftrightarrow$ constant circulation

Let generation $j$ consist of $n_j$ blobs of size $\ell_j=\lambda_j^{-1}$, velocity $U_j$. Then $E_j=n_jU_j^2\lambda_j^{-3}$ and $\Gamma_j=U_j\ell_j$. Imposing the repo's critical law $E_j=c_E/\lambda_j$:
$$n_j\Gamma_j^2=c_E\quad\Longrightarrow\quad \Gamma_j=\sqrt{c_E/n_j}\ \ \text{— scale-independent.}$$
**The critical shell normalization $E_N=c_E/N$ is exactly Kelvin's theorem across scales.** Also $c_E=n\Gamma^2$ fixes the physical meaning of $c_E$.

**Multiplicity is forced to be $O(1)$.** $M^{\rm eff}_N=(\sum_k|\hat u_k|)^2/\|u_N\|_2^2\ \ge\ \|u_j\|_\infty^2/\|u_j\|_2^2=\dfrac{c_E\lambda_j^2/n_j}{c_E/\lambda_j}=\dfrac{\lambda_j^3}{n_j}$, and for a decoherent array this is an equality up to constants. The L-02/L-11 floor $M^{\rm eff}\gtrsim N^3$ therefore reads
$$\boxed{\,n_j=O(1)\,}$$
— **the capacity no-go itself forbids many centers per scale.** LENS 2 is thus forced into *bounded-multiplicity, nested* multi-center geometry: a dyad ($n=2$) or triad ($n=3$) at every scale, never a proliferating array.

### B.3 The generation map and the closed ladder

Let generation $j$ be a skew dyad at the fixed-point angles $(\theta_*,\varphi_*)$, separation $d_j=\lambda_j^{-1}$, circulation $\Gamma=\sqrt{c_E/2}$. It generates
$$\sigma_j=c_*\Gamma\lambda_j^2 .$$
Its own compressive eigendirection $\hat e^-_j$ compresses the *internal* transverse structure of its tubes at rate $\sigma_j$. Generation $j{+}1$ is defined as that internal structure: a sub-dyad of the *same* tubes (hence the *same* $\Gamma$ — no circulation splitting), initially at separation $d_j$, compressed to
$$d_{j+1}(t)=d_je^{-\sigma_j t}\;\Rightarrow\; d_{j+1}=d_j/\rho \text{ at }\ \Delta t_j=\frac{\ln\rho}{\sigma_j}=\frac{\ln\rho}{c_*\Gamma\lambda_j^2}.$$
Then $\lambda_{j+1}=\rho\lambda_j$, $\sigma_{j+1}=\rho^2\sigma_j$, and
$$T=\sum_j\Delta t_j=\frac{\ln\rho}{c_*\Gamma\lambda_0^2}\cdot\frac{1}{1-\rho^{-2}}<\infty,\qquad
\tau_j:=T-t_j=\frac{\ln\rho}{c_*\Gamma(1-\rho^{-2})}\lambda_j^{-2}.$$
$$\boxed{\;\lambda(t)=\Big[\tfrac{\ln\rho}{c_*\Gamma(1-\rho^{-2})}\Big]^{1/2}\tau^{-1/2}\;}\qquad \gamma=\tfrac12 .$$

**Three scale-independent ratios (this is the whole design).**
$$\frac{\sigma_j}{\nu\lambda_j^2}=\frac{c_*\Gamma}{\nu},\qquad
\frac{\ell_{\rm Burgers}}{\ell_j}=\frac{2}{\theta}\sqrt{\frac{\nu}{c_*\Gamma}},\qquad
\nu\lambda_j^2\tau_j=\frac{\nu\ln\rho}{c_*\Gamma(1-\rho^{-2})} .$$
All three are **$j$-independent**, controlled by the single group $\mathrm{Re}_\Gamma=\Gamma/\nu=\sqrt{c_E/2}/\nu$. Requiring the tube to thin below the Burgers radius ($\ell_B\le\theta\,\ell_j$, $\theta=\ell/d$) gives the explicit closure threshold
$$\boxed{\;c_E\ \ge\ \frac{32\,\nu^2}{c_*^2\theta^4}\;}\qquad(c_*=0.2013,\ \theta=\tfrac12:\ c_E\ge6.3\times10^{3}\nu^2).$$
This is an independent rederivation of S3's $c_E$-collapse ($c_E\ge2\nu^2/\chi^2$) from vortex kinematics; large $c_E$ is legal since $\sum_jE_j=c_E\lambda_0^{-1}(1-\rho^{-1})^{-1}$ is made small by starting deep ($\lambda_0$ large).

### B.4 Fourier realization: lacunary bands and the exact cross-talk gate

Carriers $k_1,k_2$ with $|k_i|=\lambda_j$ at angle $\psi$; child $k_1+k_2$ has $|k_1+k_2|=2\lambda_j\cos(\psi/2)=:\rho\lambda_j$. Bands $B_j=\{k:\,|k|\in[\lambda_j(1-\eta),\lambda_j(1+\eta)]\}$.

- **Single-carrier harmonic channel is exactly dead:** $\mathcal N_{2k}\supset-iP_{2k}(k\cdot\hat u_k)\hat u_k=0$ since $k\cdot\hat u_k=0$. Doubling *requires* two distinct carriers ⟹ **multi-center is forced by the Leray algebra itself**, not assumed.
- **Cross-generation gate.** Wake band $B_{j'}$, $j'<j$, sums with $B_j$ to radii $\le\lambda_j(1+\rho^{-1})$. Demanding this miss $B_{j+1}$:
$$\boxed{\;\rho(1-\eta)\;>\;1+\rho^{-1}\;}$$
With $\psi=30^\circ$ ($\rho=2\cos15^\circ=1.932$) and $\eta=0.15$: $1.642>1.518$ ✓. With $\psi=40^\circ$ ($\rho=1.879$), $\eta=0.20$: $1.503>1.532$ ✗ — the gate is *sharp* and selects $(\psi,\eta)$.
- **Diagonal cross-talk** ($B_j+B_{j+1}$) lands at $|k|\approx\lambda_j(1+\rho)=2.93\lambda_j$, while $B_{j+1}=[1.64,2.22]\lambda_j$ and $B_{j+2}=[3.17,4.29]\lambda_j$: it falls in a **gap**, i.e. off-ladder. This is the exact-support version of S6(b), and it replaces the repo's failed gate (measured $9\%$ signal-to-leakage). Leakage is then a *budget* item (a scale-independent loss fraction absorbed by §B.3's threshold), not a contamination.

### B.5 In-box decoherence: exact identity, and correction to the seed estimate

Clouds $\hat u_A(k)=\Phi(k-Np)Ae^{-ik\cdot x_0}$, $\hat u_B(k)=\Phi(k-Nq)Be^{-ik\cdot x_1}$, $\Delta=x_1-x_0$. Writing $k_A=Np+\alpha$, $k_B=Nq+\beta$, $\alpha+\beta=\delta$:
$$\mathcal N^{AB}_{N(p+q)+\delta}=N\,e^{i\Theta}\,P_\bullet V_0\;K_\Delta(\delta),\qquad K_\Delta(\delta)=\sum_\alpha\Phi(\alpha)\Phi(\delta-\alpha)e^{i\alpha\cdot\Delta}.$$
Generating function: $\sum_\delta K_\Delta(\delta)e^{i\delta\cdot y}=\varphi(y+\Delta)\varphi(y)$ with $\varphi=\sum_\alpha\Phi(\alpha)e^{i\alpha\cdot y}$. Hence the **exact** power-suppression factor
$$\boxed{\;S(\Delta)=\frac{\|K_\Delta\|^2_{\ell^2}}{\|K_0\|^2_{\ell^2}}=\frac{\int|\varphi(y+\Delta)|^2|\varphi(y)|^2dy}{\int|\varphi|^4dy}\;}$$
i.e. *cross-talk is exactly the physical-space envelope overlap*. For the Fejér box $\Phi(\alpha)=\prod_i(1-|\alpha_i|/W)_+$, $\varphi=\prod_iF_W$, $F_W(t)=\frac1W\frac{\sin^2(Wt/2)}{\sin^2(t/2)}$, $|F_W(t)|\le4/(Wt^2)$, $\int F_W^2\asymp W$, $\int F_W^4\asymp W^3$:
$$S(\Delta)\asymp\prod_{i:\Delta_i\neq0}(W\Delta_i)^{-4}\;\xrightarrow[\text{axis-aligned}]{}\;(W|\Delta|)^{-4}.$$
**Correction:** seed S6(a)'s $(W|\Delta|)^{-3}$ is the *random-phase mode-counting* heuristic ($W^{3}$ vs $W^{6}$); the exact overlap identity gives $-4$ (sharp box or Fejér) and $e^{-cW^2|\Delta|^2}$ for Gaussian envelopes. Use $-4$.

### B.6 THEOREM (S6(a) is unfixable) — cross-center pairing cannot beat cross-center cross-talk

Let $u_A=\mathrm{Re}[e^{iNp\cdot x}a(x-x_0)]$, $u_B$ likewise, $\hat a,\hat b$ supported in $|\alpha|\le W$. Then for **every** child projector $P_C$ (intended or not),
$$\|P_CB(u_A,u_B)\|_2\le(N|q|+W)\,\big\||a(\cdot)|\,|b(\cdot-\Delta)|\big\|_{L^2}=(N|q|+W)\,O_\Delta,$$
because $|(u_A\cdot\nabla)u_B|\le|u_A||\nabla u_B|$ pointwise and $\mathbb P$ is $L^2$-bounded. Moreover all channels are slices of the *same* sequence $K_\Delta$, so
$$\frac{\|P_IB\|_2^2}{\|P_UB\|_2^2}=\frac{\|K_\Delta\|^2_{\ell^2(I)}}{\|K_\Delta\|^2_{\ell^2(U)}}\ \xrightarrow[\ |\Delta|\ \text{large}\ ]{}\ \frac{|I|}{|U|},$$
since $\varphi(\cdot+\Delta)\varphi(\cdot)$ becomes supported on a set $\ll W^{-3}$ and its Fourier coefficients become nearly flat on the $\delta$-box. **Translation-split reduces the absolute intended channel by $S(\Delta)\to0$ while driving the intended/unintended ratio down to the pure counting ratio $|I|/|U|\le O(1)$. It strictly loses. REJECTED, and unfixable.**

**The fix actually used.** NSDL puts the two lobes at $|\Delta|=d_j=\lambda_j^{-1}$, so $W|\Delta|=\eta\lambda_j\cdot\lambda_j^{-1}=\eta<1$: the intended pairing is **fully coherent, $S=O(1)$**. Inter-center separation is *not* used to suppress anything. Suppression is instead carried by the **lacunary band geometry (B.4)** — a support (measure-zero) exclusion, not an amplitude decay. Meanwhile the low-$k$ (Reynolds-stress / Biot–Savart) channel that carries the mutual strain is *not* suppressed at all, because $\sigma_j=c_*\Gamma\lambda_j^2$ comes from the $|k|\lesssim2\eta\lambda_j$ difference channel of each lobe with itself — an $O(1)$-in-$\eta$, coherent, sign-definite quantity. **Channel separation is the new structure: cores strongly coupled at low $k$ (strain), exactly gated at high $k$ (relay).**

---

## C. Scaling table ($\tau=T-t$, $\lambda=N\sim\tau^{-1/2}$, $\Gamma=\sqrt{c_E/2}$)

| quantity | formula | exponent in $\tau$ |
|---|---|---|
| energy $E$ | $c_E\lambda_0^{-1}(1-\rho^{-1})^{-1}$ | $\tau^{0}$ (bounded ✓) |
| enstrophy | $\sum\lambda_j^2E_j\asymp c_E N$ | $\tau^{-1/2}$ |
| global $\|u\|_{L^3}$ | $\big(\mathcal R\,J\,c_E^{3/2}\big)^{1/3}$, $J=\log_\rho(N/\lambda_0)$ | $(\log\tau^{-1})^{1/3}$ **(log-divergent)** |
| $\|\omega\|_\infty$ | $\Gamma N^2$ | $\tau^{-1}$ |
| dissipation rate $\nu\|\nabla u\|_2^2$ | $\nu c_EN$ | $\tau^{-1/2}$ (so $\int_0^T<\infty$ ✓) |
| nonlinear $\|\mathbb P(u\cdot\nabla)u\|_2$ | $c_EN^{3/2}$ | $\tau^{-3/4}$ |
| pressure $\|\nabla p\|_2$ | $\le\|(u\cdot\nabla)u\|_2$, same order | $\tau^{-3/4}$ |
| physical time left | $\tau=\frac{\ln\rho}{c_*\Gamma(1-\rho^{-2})}N^{-2}$ | $\tau^{1}$ |
| Fourier bandwidth | $N$, band half-width $\eta N$ | $\tau^{-1/2}$ |
| active mode count | $M^{\rm eff}\asymp N^3/n$, $n=2$ | $\tau^{-3/2}$ |
| centers per generation | $n_j=2$ | $\tau^0$ (bounded, forced by §B.2) |
| generations alive | $J\sim\tfrac12\log_\rho\tau^{-1}$ | $\log$ |
| strain | $\sigma=c_*\Gamma N^2$ | $\tau^{-1}$ |
| Type-I indicator | $\sqrt\tau\|u\|_\infty\to\Gamma$ | $\tau^0$ **(marginal Type I)** |
| Serrin pair $2/q+3/p=1$ | $\int\|u\|_p^q\,dt\asymp\Gamma^q\!\int\tau^{-1}d\tau$ | $\log$-divergent **in every pair** |

Every critical diagnostic diverges *logarithmically* — the least-excluded possible profile, and exactly the repo's $\beta=-1$, $\sigma=\gamma$ log-critical branch, here realized geometrically.

---

## D. Closed feedback loop (each arrow a formula)

$$
\underbrace{\text{skew dyad }(\theta_*,\varphi_*),\ d_j}_{\text{geometry}}
\xrightarrow{\;\sigma_j=c_*\Gamma/d_j^2\;}
\underbrace{\text{mutual axial strain}}_{\text{low-}k\text{ channel}}
\xrightarrow{\;\dot d_{j+1}=-\sigma_jd_{j+1}\;}
\underbrace{d_{j+1}=d_j/\rho\ \text{after}\ \Delta t_j=\ln\rho/\sigma_j}_{\text{compression}}
$$
$$
\xrightarrow{\;\lambda_{j+1}=\rho\lambda_j,\ k_1+k_2\in B_{j+1}\;}
\underbrace{\text{child carriers}}_{\text{high-}k\text{ relay}}
\xrightarrow{\;\mathsf R_{\theta_0}\;}
\underbrace{\text{re-locked skew dyad}}_{\text{phase code}}
\xrightarrow{\;\Gamma_{j+1}=\Gamma_j\ (\text{Kelvin})\;}
\underbrace{\sigma_{j+1}=\rho^2\sigma_j}_{\text{closure}}
$$
Gain per turn $=\rho^2>1$; time per turn $=\ln\rho/\sigma_j\propto\rho^{-2j}$; $\sum\Delta t_j<\infty$.
Loss arrows: viscous thinning floor $\ell_B/\ell=\tfrac2\theta\sqrt{\nu/(c_*\Gamma)}$; off-ladder leakage into band gaps; cross-diffusion $\propto-\omega_1\cdot\omega_2=+0.08\Gamma^2/\ell^4$ (non-cancelling). **All losses are $j$-independent fractions**, so the loop closes iff the single inequality $c_E\ge32\nu^2/(c_*^2\theta^4)$ holds (S3 $c_E$-collapse).

The rotation $\mathsf R_{\theta_0}$ is the *anti-alignment device*: without it the child inherits the parent frame, $\omega$ aligns with the parent's stretching eigenvector, and $t\cdot St\to0$ (classical depletion). With it the child re-enters at $(\theta_*,\varphi_*)$ relative to its own partner.

---

## E. Obstruction audit — exact collision points

1. **Energy bound (F-N1).** $E$ bounded by construction; never used as a signature. No collision.
2. **Finite dissipation (F-N2).** $\nu\int_0^T\!\|\nabla u\|^2\!\propto\!\int\tau^{-1/2}d\tau<\infty$ — *satisfied*, as required. No collision.
3. **ESS $L^\infty_tL^3$.** Collision point: ESS needs $\sup_t\|u\|_3<\infty$. Here $\nu\lambda_j^2\tau_j=\nu\ln\rho/(c_*\Gamma(1-\rho^{-2}))$ is $j$-independent, so each completed band retains the fixed fraction $\mathcal R=e^{-2\nu\ln\rho/(c_*\Gamma(1-\rho^{-2}))}\in(0,1)$ until $T$; each contributes $\asymp\Gamma^3$ to $\|u\|_3^3$ (scale-free), so $\|u\|_3^3\gtrsim\mathcal R\,J(t)\,\Gamma^3\to\infty$. ESS's hypothesis fails; ESS is *used* (contrapositive) as the blow-up criterion.
4. **Fixed-finite-bandwidth no-go (F-6/F-α1/VR-L-011).** Hypothesis: trajectory confined to a *fixed* mode set. Here $N(t)=\Theta(\tau^{-1/2})\to\infty$; the datum is band-limited but the trajectory is not (VR-L-012 explicitly permits this). Evaded at the hypothesis.
5. **Galerkin global existence.** Same hypothesis, same evasion; the Grönwall bound $\|u\|\le\|u_0\|$ *holds here* and is consistent (energy is bounded).
6. **Pure-swirl $L^3$ no-go (VR-L-016 / LG-9).** Collision point: LG-9 uses rotational equivariance of $p=R_iR_j(u_iu_j)$ under axisymmetry to force $\partial_\theta p_0\equiv0$, hence $u\cdot\nabla p\equiv0$ and $F'(0)\le0$. NSDL's symmetry group is at most $C_2$ (the $\pi$-rotation of §B.1), never $SO(2)$, and $u^\theta$ is not the only component; so $\partial_\theta p\not\equiv0$ and the pressure channel $P$ is unconstrained in sign. **Verification is a pilot deliverable (F/P2-v).**
7. **One-scale self-similar (NRS1996 / Tsai1998).** Collision point: both require the *exact continuous* form $u=(2a\tau)^{-1/2}U(x/\sqrt{2a\tau})$ for all $t$ near $T$. NSDL satisfies only the **discrete, rotated, windowed** relation $u(x,t_j)\approx\rho\,\mathsf R_{\theta_0}u(\rho\mathsf R_{\theta_0}^{-1}x,\,t_{j+1})$ on the band window $[\lambda_0,\lambda_{J(t)}]$, with $\theta_0\neq0\ (\mathrm{mod}\ 2\pi)$ and $\rho=1.932\neq$ any parabolic dilation of a continuous group orbit. Additionally NRS needs $U\in L^3$; here the profile's per-octave $L^3$ mass is constant so $U\notin L^3$ (only $L^{3,\infty}$). Tsai needs finite local energy *and* exact self-similarity: local energy is finite here ($\sum_jE_j<\infty$, $\int\!\!\int|\nabla u|^2<\infty$), so the *only* escape is the failure of exact self-similarity — which is the design's rotation $\mathsf R_{\theta_0}$ and its growing window. **This is a load-bearing, checkable point, not a hand-wave.**
8. **Smooth-forcing high-frequency decay (F-N4).** $f\equiv0$. Vacuous.
9. **Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}/N^3$.** THE central collision. $\gamma=1$ here: $M^{\rm eff}\asymp mN^3$ ($m$ = fill fraction) because each generation is a *quasi-isotropic blob of size $1/N$* (§B.2), not a long tube. Insert the design's own $\tau$ (not free): $\tau=\sigma_j^{-1}\lambda_j^2=\ln\rho/(c_*\Gamma)$. Then
$$2\kappa^2\tau^2c_Em=2\kappa^2m\,c_E\frac{\ln^2\rho}{c_*^2\Gamma^2}=\frac{4\kappa^2m\ln^2\rho}{c_*^2}\quad(\text{using }c_E=2\Gamma^2),$$
**independent of $c_E$ and of $N$** — the exact $c_E$-collapse of S3. Feasibility ($\ge1/2$) needs $8\kappa^2m\ln^2\rho\ge c_*^2$: with $\rho=1.932$, $\kappa=1.15$, $c_*=0.2013$ this needs $m\ge0.0088$, i.e. a $0.9\%$ fill fraction. **Passed with two orders of magnitude of margin.** Necessary, not sufficient — stated as such.
10. **Diagonal cross-talk gate.** Passed by exact support geometry $\rho(1-\eta)>1+\rho^{-1}$ (§B.4). This replaces amplitude suppression (which §B.6 proves cannot work) by a measure-zero exclusion. The repo's failed $9\%$ gate used $\rho=2$ with box-filling clouds, for which the gate inequality reads $2(1-\eta)>1.5$ i.e. $\eta<0.25$ — never checked, and violated by the $W\!\to\!2W$ convolution broadening. Here $\eta$ is chosen *from* the inequality.
11. **CSTY2009 axisymmetric Type-I exclusion.** NSDL **is** Type I ($\sqrt\tau\|u\|_\infty\to\Gamma$). Collision point: CSTY's hypothesis is "$v$ axisymmetric strong solution". NSDL is non-axisymmetric by construction ($t_1\cdot t_2=-0.08$, symmetry group $C_2$, per-generation rotation $\theta_0$). **No registered theorem excludes non-axisymmetric Type-I blow-up.** This is the single most load-bearing evasion and the candidate's greatest exposure: if a future non-axisymmetric Type-I exclusion appears, NSDL dies immediately.
12. **KNSS2009 ancient-solution Liouville.** Applies: Type-I rescaling yields (subsequence) a nonzero bounded ancient mild solution. Ours is the DSS ancient solution $u_\infty(x,s)=\rho\mathsf R u_\infty(\rho\mathsf R^{-1}x,\rho^2s)$, $s<0$. KNSS(a) needs axisymmetric **no swirl** — fails. KNSS(b) needs axisymmetric with swirl **and** $|u|\le C/r$ — fails on both clauses ($\|u\|_\infty\asymp\Gamma\lambda$ near the axis of nothing; there is no distinguished axis). The general bounded-ancient Liouville problem is open, so KNSS gives an *obligation* (nontriviality of the limit) not a contradiction.
13. **Front-resolution threat model (TM-22, TM-04, TM-20).** Front thickness $=\ell_J=\lambda_J^{-1}$; PREREGISTERED_MIN_POINTS_PER_FRONT $=7$ ⟹ grid $\ge7\kappa\lambda_J$ per direction. A $J=2$ pilot at $\lambda_0=8$, $\rho=1.932$ gives $\lambda_2=30$, $\kappa\lambda_2=34$, grid $\ge238$ ⟹ $256^3$. Time integrator: RK4 or exponential integrator, **never Heun+central differences** (TM-20 $|G|^2=1+\alpha^4/4$). $2/3$ dealiasing (TM-03).
14. **U-X1/U-X2 (critical $L^3$ / anisotropy).** U-X1 assumes a *uniformly $L^3$-bounded single-scale* rescaling trajectory. NSDL is neither ($\|u\|_3\to\infty$, and $J(t)\to\infty$ scales are simultaneously alive). Evasion is by **wake accumulation**, not by anisotropy — note this differs from the route U-X2 anticipates, and is therefore new.

---

## F. Minimal falsification experiment ($\le$1 h)

**P1 — skew-strain fixed point (10 min, exact/interval).** Symbolically build $\sigma(\theta,\varphi)$ (§B.1) and the tangent ODE $\dot t=St-(t\!\cdot\!St)t$ reduced to $(\theta,\varphi)$ by the $C_2$ symmetry. Variables: $(\theta,\varphi)$. Compute the Jacobian eigenvalues at $(\theta_*,\varphi_*)$ in interval arithmetic (`mpmath.iv`); rationals where possible. **Success:** $\sigma_*>0$ and the fixed point has at most one unstable direction. **Kill:** $\sigma_*\le0$, or $\ge2$ unstable directions (then alignment depletion cannot be locked out and the ladder has no codim-1 shooting parameter for PO-09). No repo module needed.

**P2 — one-generation lacunary relay (40 min, float + exact gate).** Grid $192^3$ (fallback $128^3$), $\nu$ from $\mathrm{Re}_\Gamma=\Gamma/\nu\in\{20,60\}$, $\lambda_0=8$, $\psi=30^\circ$ ($\rho=1.932$), $\eta=0.15$, $\theta=1/2$, $c_E=2\Gamma^2$.
- Build the skew dyad parent: two Fejér/Gaussian-enveloped, mode-wise Leray-projected blobs at carriers $k_1,k_2$, critical energy $c_E/\lambda_0$. Reuse `leray_response_relay.leray_project`, `leray_response_relay.leray_advection` (dealiased $\mathbb P((u\cdot\nabla)u)$), `leray_response_relay.fejer_carrier_packet` (envelope only — **do not** reuse `relay_stage`/`adjoint_child_response`, that packet family is REJECTED).
- Channel tagging with **lacunary** band tags $B_{j-1},B_j,B_{j+1}$ and the gap set $G$: adapt `mesoscopic_local_fft.measure_local_fft_cloud` / `local_fft_leray_coefficients` (zero-padded, no wrap aliasing).
- Full evolution to $t=\tau\lambda_0^{-2}$, $\tau=\ln\rho/(c_*\Gamma)$, via `mesoscopic_galerkin.run_small_mesoscopic_galerkin`-style RK4$\times$16 driver.
- Measure: (i) $D_N=E_{B_{j+1}}/E_{B_j}(0)$; (ii) gap leakage $E_G/E_{B_{j+1}}$; (iii) $M^{\rm eff}$; (iv) induced low-$k$ strain $\sigma_{\rm meas}$ from $|k|\le2\eta\lambda_0$ modes vs. predicted $c_*\Gamma\lambda_0^2$ — **including sign**; (v) shape overlap $\varrho=|\langle \hat U_{j+1},\mathsf R_{\theta_0}\hat U_j\rangle|/(\|\cdot\|\|\cdot\|)$; (vi) instantaneous $L^3$ generation rate $J(u)=3\!\int\! p\,\nabla\!\cdot\!(|u|u)+V$ (VR-L-016 check that $P\not\equiv0$).
- **Success (pre-registered):** $D_N\ge0.05$; leakage ratio $\le1$; $\sigma_{\rm meas}/(\Gamma\lambda_0^2)\in[0.1,0.4]$ **and positive**; $\varrho\ge0.5$; all stable to $\le20\%$ across $\{128^3,192^3\}$ and $\lambda_0\in\{8,12\}$.
- **Kill (decisive):** $\sigma_{\rm meas}\le0$ (no compressive feedback ⟹ loop D is open ⟹ mechanism dead); or leakage $>5\times$ child energy; or $D_N$ decaying faster than $\lambda_0^{-1/4}$ between the two scales.
- **Exact arithmetic required:** band-membership / support-disjointness gate (integer lattice); the trilinear neutrality $\sum_k\Pi_k=0$ on a small $\lambda_0=2$ replica via `exact_carrier_record_verifier.verify_serialized_strict_orientation_records`; P1 eigenvalues (interval). **Float acceptable:** $D_N$, $M^{\rm eff}$, leakage, overlaps, $\sigma_{\rm meas}$.

---

## G. Proof chain (10 obligations)

1. Skew-strain fixed point $(\theta_*,\varphi_*)$ exists, $\sigma_*=c_*\Gamma/d^2>0$, codim $\le1$ unstable — interval-verified.
2. Explicit smooth divergence-free band-limited generation field $U_j$ with exact Leray sum-channel transfer coefficient $\chi>0$ into $B_{j+1}$, uniform in $j$ (the repo's #1 unproved lemma: *positive signed flux margin*).
3. Lacunary gate: $\rho(1-\eta)>1+\rho^{-1}$ ⟹ all unintended sum/difference channels lie outside $\bigcup_jB_j$; leakage fraction $\le\lambda_{\rm leak}<1$ uniformly in $j$.
4. Budget: $\chi c_E^{3/2}-\nu c_E\,\mathcal L(\text{shape})>0$ for $c_E\ge32\nu^2/(c_*^2\theta^4)$ (net positive per-generation flux; $c_E$-collapse).
5. Recurrence: the renormalized generation map $\mathfrak G=\mathcal E\circ\mathsf R_{\theta_0}^{-1}\circ\mathrm{Dil}_\rho\circ\Phi_{\Delta t_j}$ has an attracting fixed point in a Banach space of band-limited divergence-free fields (radii-polynomial $Y+Z(r)<r$; PO-04/05/07/13).
6. Zeno time: $\sum_j\Delta t_j=\ln\rho/(c_*\Gamma\lambda_0^2(1-\rho^{-2}))<\infty$ with rigorous per-generation rate bounds (PO-10).
7. Wake retention: $\nu\lambda_j^2\tau_j$ is $j$-independent ⟹ retention $\mathcal R>0$ ⟹ $\|u(t)\|_3^3\ge c\,\mathcal R\,J(t)\to\infty$ (S4 lemma, made rigorous; PO-11 in the critical norm).
8. Non-extendability: $\sup_t\|u\|_3=\infty$ and divergence in every Serrin pair ⟹ ESS/Prodi–Serrin contrapositive; corroborated by $\int_0^T\|\omega\|_\infty=\infty$ (BKM).
9. Entry (PO-09): by obligation 1 the unstable set is codim $\le1$ ⟹ reduce to one-parameter bisection (edge tracking) in the initial skew angle from explicit $u_0$; PO-02/03/08/12.
10. Assemble PO-01,02,03,12,13,14 + 1–9 ⟹ CLAY-D with $f\equiv0$ (hence $\neg$CLAY-B). $\mathbb R^3$/TARGET-U requires HS-5 separately.

---

## REJECTED sub-variants (kept, with the exact failing equation)

- **M-1 (Moffatt–Kimura single dyad).** $\dot d=-\alpha_v\Gamma/d\Rightarrow d^2=2\alpha_v\Gamma\tau$, $\sigma=a/\tau$; core $\frac{d\ell^2}{d\tau}=\frac{a\ell^2}{\tau}-4\nu\Rightarrow\ell^2\to\frac{4\nu\tau}{a-1}$, hence $\ell/d=\sqrt{2\nu/((a-1)\alpha_v\Gamma)}=\text{const}$: **exactly one-scale parabolic self-similar** ⟹ NRS/Tsai + CSTY. REJECTED.
- **M-2 (straight skew pair alone).** $\hat s\cdot e_1\equiv0\Rightarrow\dot d=0\Rightarrow\sigma=\text{const}\Rightarrow\ell\to\ell_B=\sqrt{4\nu/\sigma}$, $\omega\to\Gamma\sigma/(4\pi\nu)$ bounded. Burgers steady state, no blow-up. REJECTED (repaired by nesting, §B.3).
- **M-3 (circulation splitting, $\Gamma_j=\Gamma_02^{-j}$).** $\sigma_j/(\nu\lambda_j^2)=c_*\Gamma_0/(\nu\lambda_j)\to0$: viscosity wins at every deep scale. REJECTED.
- **M-4 (Crow crinkle array, $n_j\propto\lambda_j^2$ lobes).** $M^{\rm eff}\asymp\lambda_j^3/n_j\asymp\lambda_j\ll\lambda_j^3$, violating $M^{\rm eff}\gtrsim N^3$. REJECTED — and this is *why* §B.2 forces $n_j=O(1)$.
- **M-5 (symmetric cone triad, sum map).** Child cone: $\tan\chi'=\tfrac12\tan\chi\Rightarrow\tan\chi_j=2^{-j}\tan\chi_0\to0$; the degenerate limit is the single-carrier harmonic channel with $\mathcal N_{2k}=-iP_{2k}(k\cdot\hat u_k)\hat u_k=0$ **exactly**. REJECTED.
- **M-6 (3-carrier sum-triad similarity fixed point).** Equal moduli force $n_i\!\cdot\!n_j=c$; similarity forces $\frac{3c+1}{2+2c}=c\Rightarrow2c^2-c-1=0\Rightarrow c\in\{1,-\tfrac12\}$: $c=1$ parallel, $c=-\tfrac12$ gives $m_i=-n_{i+2}$, $\rho=1$ (no dilation). **No nontrivial 3-carrier period-1 fixed point** ⟹ $\ge4$ carriers or a period-$L$ cycle (S6(b) eq. 6.11). REJECTED as stated.
- **M-7 (S6(a) translation-split repair).** PROVED UNFIXABLE, §B.6: $\|K_\Delta\|^2_{\ell^2(I)}/\|K_\Delta\|^2_{\ell^2(U)}\to|I|/|U|$ while both collapse by $S(\Delta)\asymp(W|\Delta|)^{-4}$.
