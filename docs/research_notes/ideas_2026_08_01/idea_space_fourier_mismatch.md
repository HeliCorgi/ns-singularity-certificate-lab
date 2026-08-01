# LENS 3 — The Spectrally Fat Cell: physical width $L$ and Fourier bandwidth $N$ diverging at different rates

**STATUS: FORMAL ANSATZ.** (Lemma 2 and the corner algebra of §B.4 are SYMBOLIC CANDIDATE; sub-variants V1, V2 are REJECTED below with their exact failing equations.)

---

## A. Clay target

**Target: TARGET-U** (unforced $\mathbb R^3$ breakdown; strictly stronger than (C)); the same object transfers verbatim to **(D)** on $\mathbb T^3$ if the whole-space tail bookkeeping fails.

- **Domain:** $\mathbb R^3\times(0,T)$. Fourier/Leray bookkeeping is done with the repo's torus convention ($\hat u_k$, $P_k=I-k\otimes k/|k|^2$, $\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$) applied to the *local* band of the concentration cell; per U-X4/HS-5 this is a computational convention only and **no periodic closure is imported as a whole-space claim**.
- **Forcing:** none. (F-N4 makes smooth forcing useless for direct high-shell injection; F-N1/N2 make energy-based signatures dead. We take no forcing.)
- **Initial regularity/decay:** $u_0\in\mathcal S(\mathbb R^3)$, $\nabla\cdot u_0=0$, with **sub-Gaussian (Gevrey) decay** $|\partial^\alpha u_0|\le C_\alpha e^{-c|x|^{a}}$, $a\in(0,1)$. This satisfies Clay's condition (4) exactly. **We deliberately drop the repo's compact-support strengthening** — that is the entire price paid to dissolve Z-01/Z-02 (§B.5).
- **Viscosity:** fixed $\nu>0$, no vanishing-viscosity limit. The mechanism strengthens as $\nu\downarrow0$ but does not require it.

---

## B. Central mathematics

### B.1 State variables and the mismatch

Let $\tau=T-t$. Two **independent** scales:
$$L(\tau)=\text{physical concentration width},\qquad N(\tau)=\text{top of the active Fourier band},$$
active band $\lambda\in[L^{-1},N]$, and the **mismatch variable**
$$x(\tau):=\log\!\big(N(\tau)L(\tau)\big),\qquad\text{LENS 3 demands } x\to+\infty .$$
Wave-packet-limited objects have $x\equiv O(1)$; one-scale self-similar objects have $x\equiv$ const. Intra-band spectrum: the repo's **corrected log-critical branch** $\beta=-1$, $\sigma=\gamma$, i.e. $E_j=c_E/\lambda_j$ for $L^{-1}\le\lambda_j\le N$, with $c_E=c_E(\tau)$ now a *dynamical* quantity, not a constant.

Support measure $|\Omega|\asymp L^3$ (band-filled cell — see V2 for why not sparse). Per-octave amplitude
$$U_j=\Big(\tfrac{2E_j}{|\Omega|}\Big)^{1/2}=\Big(\tfrac{2c_E}{\lambda_j L^3}\Big)^{1/2},\qquad \text{strain of octave }j:\ S_j=\lambda_jU_j=\Big(\tfrac{2c_E}{L^3}\Big)^{1/2}\lambda_j^{1/2}. \tag{B.1}$$

### B.2 Lemma 1 (phase-space participation) — where the mesoscopic floor really bites

For $u$ with spectral support on $M$ lattice modes and spatial support of measure $|\Omega|$, $\|u\|_1\le|\Omega|^{1/2}\|u\|_2$ and $|\hat u_k|\le(2\pi)^{-3}\|u\|_1$ give
$$M^{\rm eff}=\frac{(\sum_k|\hat u_k|)^2}{\sum_k|\hat u_k|^2}\ \le\ c\,\min\big(M,\ M^2|\Omega|\big). \tag{B.2}$$
A Gevrey blob filling $B_L$ with band top $N$ saturates it: $M\asymp N^3$, $M^2|\Omega|\asymp N^6L^3=N^3e^{3x}\gg M$, so $M^{\rm eff}\asymp N^3$. The mesoscopic floor $M^{\rm eff}\gtrsim(\nu^2/c_E)N^3$ is then **not a mode-count condition at all** — it reads
$$1\ \ge\ \nu^2/c_E(\tau), \tag{B.3}$$
an *amplitude* condition. It is the growth of $c_E$, not of $M$, that must beat it.

### B.3 Lemma 2 (strain ratio) — the engine of the mismatch

From (B.1), the strain that advances the **front** and the strain that contracts the **envelope** are
$$S_N=NU_N=\Big(\tfrac{2c_E}{L^3}\Big)^{1/2}N^{1/2},\qquad S_L=L^{-1}U_L=\Big(\tfrac{2c_E}{L^3}\Big)^{1/2}L^{-1/2},$$
$$\boxed{\ \frac{S_N}{S_L}=(NL)^{1/2}=e^{x/2}\ } \tag{B.4}$$
*Dimension check:* both sides $\rm s^{-1}/s^{-1}$, RHS dimensionless ✓. *Sign check:* $S_N>S_L$ for $x>0$ ✓. Physically: a strain field of correlation length $N^{-1}$ cannot coherently compress a blob of size $L\gg N^{-1}$; the envelope only feels its own octave. Hence the two clocks
$$\frac{d\log N}{dt}\ \lesssim\ \kappa_N S_N,\qquad -\frac{d\log L}{dt}=\kappa_L S_L,\qquad \frac{dx}{dt}=S_L\big(\kappa_Ne^{x/2}-\kappa_L\big). \tag{B.5}$$
$x$ is **self-reinforcing** once $x>2\log(\kappa_L/\kappa_N)$: no fine-tuning is needed to make $NL\to\infty$. This is the positive statement of the lens.

### B.4 The constraint polytope and the unique corner

Write $L=\tau^{\mu}$, $N=\tau^{-\gamma}$, $c_E=c_0\tau^{-\varepsilon}$, $\delta=\gamma-\mu>0$.

| constraint | source | inequality |
|---|---|---|
| (P1) bounded energy $\sum_jE_j\asymp2c_EL\le E_0$ | Leray / F-N1 | $\varepsilon\le\mu$ |
| (P2) finite dissipation $\nu\!\int\!\|\nabla u\|_2^2\,dt<\infty$, rate $\asymp\nu N^2E_N$ | Leray / F-N2 | $\gamma+\varepsilon<1$ |
| (P3) Type II (escape CSTY) $\sqrt\tau\|u\|_\infty\to\infty$ | CSTY2009 | $\mu+\varepsilon/2>1/2$ |
| (P4) BKM $\int\|\omega\|_\infty dt=\infty$ | Kato–Ponce | $\gamma+\varepsilon+3\mu\ge2$ |
| (P5) envelope clock $-d\log L/dt=\kappa_LS_L$ | (B.5) | $\varepsilon+4\mu=2$ |
| (P6) front clock $d\log N/dt=\kappa_NS_N$ | (B.5) | $\gamma+\varepsilon+3\mu=2$ |
| (P7) transport bound $|\dot L|\le U_L$ | Route-12 rejection | $\kappa_L\le1$ |

**REJECTED sub-variant V1 (pure power-law self-similarity).** (P5) and (P6) are simultaneously satisfiable only if
$$(\gamma+\varepsilon+3\mu)-(\varepsilon+4\mu)=\gamma-\mu=0\ \Longrightarrow\ \boxed{\delta=0}.$$
So **no exact two-parameter power law can realise $NL\to\infty$** while both clocks are self-consistent. This is the sharp, quantitative form of the Z-01/Z-02 kinematic incompatibility, and it kills the naive lens-3 ansatz. *(Exact failing equation: $\gamma+\varepsilon+3\mu=2=\varepsilon+4\mu$.)*

**The escape.** (P6) is *not* the correct front law. The front is not free-running: it is **viscously pinned**. The front stalls where advective transfer equals viscous drain,
$$S_N=\nu N^2\ \Longleftrightarrow\ \Big(\tfrac{2c_E}{L^3}\Big)^{1/2}N^{1/2}=\nu N^2\ \Longrightarrow\ N^3=\frac{2c_E}{\nu^2L^3},$$
$$\boxed{\ NL=\Big(\frac{2c_E}{\nu^2}\Big)^{1/3}=\Big(\frac{E_0}{\nu^2L}\Big)^{1/3}\ }\tag{B.6}$$
using the **saturated** energy budget $c_E=E_0/(2L)$ from (P1). *Dimension check:* $[c_E]=[E]$; $[\nu^2]=L^4T^{-2}$; in the repo's normalization $c_E$ carries $L^3T^{-2}\cdot L^{-1}\!\cdot\!L$… concretely $2c_E/\nu^2$ has dimension $L^{-3}$ so $(2c_E/\nu^2)^{1/3}$ is dimensionless-in-$NL$ ✓. *Sign check:* $NL$ increases as $L\downarrow0$ ✓. **(B.6) is the whole idea:** the mismatch is not postulated, it is *generated* by envelope contraction feeding $c_E$. Note (B.6) differs from the Kolmogorov relation $NL\sim\mathrm{Re}_L^{3/4}$ precisely because the repo's critical shell law $E_j\propto\lambda_j^{-1}$ is steeper than K41's $\lambda_j^{-2/3}$.

Now (P6) is replaced by (B.6), i.e. $\delta=\mu/3$, and (P5) closes the system:
$$\varepsilon=\mu\ \ (\text{P1 saturated}),\qquad \varepsilon+4\mu=2\ \Rightarrow\ \mu=\varepsilon=\tfrac25,\qquad \delta=\tfrac{2}{15},\qquad \gamma=\mu+\delta=\tfrac{8}{15}.$$
Explicitly, integrating $\dot L=-\kappa_LU_L=-\kappa_L\sqrt{E_0}\,L^{-3/2}$:
$$L(\tau)=\Big(\tfrac52\kappa_L\sqrt{E_0}\,\tau\Big)^{2/5},\qquad c_E(\tau)=\tfrac{E_0}{2}L^{-1}\propto\tau^{-2/5},\qquad N(\tau)=\frac{(E_0/\nu^2)^{1/3}}{L^{4/3}}\propto\tau^{-8/15}. \tag{B.7}$$

**Verification against the remaining constraints (all strict):**
(P2) $\gamma+\varepsilon=\frac8{15}+\frac6{15}=\frac{14}{15}<1$ ✓ margin $1/15$.
(P3) $\mu+\varepsilon/2=\frac25+\frac15=\frac35>\frac12$ ✓ Type II.
(P4) $\gamma+\varepsilon+3\mu=\frac{8+6+18}{15}=\frac{32}{15}>2$ ✓ **power-law** BKM divergence, not merely logarithmic.
(P7) $\kappa_L\le1$ imposed by construction ✓.
(B.3) $\nu^2/c_E\propto\tau^{2/5}\to0$ ✓ mesoscopic floor evaded with diverging margin.

### B.5 Concrete packet construction and exact Littlewood–Paley control (Z-01/Z-02 head-on)

Z-01/Z-02 states: *compact-support divergence-free packets cannot simultaneously be exactly band-limited; LP leakage and translation cross terms are uncontrolled.* The construction is built **spectrally first**:

1. Fix a Gevrey-class radial cutoff $\chi\in G^{1/a}$, $a\in(0,1)$, supported in $[2^{-1/2},2^{1/2}]$, $\sum_j\chi(\lambda/\lambda_j)^2=1$ on $[L^{-1},N]$.
2. Set $\hat u_j(k)=s_j\,\chi(|k|/\lambda_j)\,P_k\,\hat g_j(k)$, with $\hat g_j$ a Gevrey envelope of width $L^{-1}$ in $k$ about the annulus, all octaves **co-located at the same centre $x_*$** (zero translation offsets).
3. $u=\sum_ju_j$, $s_j^2$ chosen so $\|u_j\|_2^2=2c_E/\lambda_j$.

Then: **(i)** $P_k$ is a modewise Fourier multiplier, so Leray projection **exactly preserves the annulus support** ($P_k$ is smooth and $0$-homogeneous, singular only at $k=0$, which lies in no annulus). Divergence-freeness and exact band-limitation are therefore *simultaneously exact* — the Z-01/Z-02 obstruction only forbids **compact spatial support** plus band-limitation. **(ii)** LP leakage is identically zero: $\|P_ju\|_2^2=\|u_j\|_2^2=2c_E/\lambda_j$ *exactly*, and every cross term $\langle P_ju_i,P_ju_{i'}\rangle$ vanishes by disjoint spectral support. There are no translation cross terms because $x_j\equiv x_*$. **(iii)** The price: spatial support is not compact; a Gevrey annulus cutoff gives $|u_j(x)|\lesssim e^{-c(|x-x_*|/L)^{a}}$, so tail mass outside $B_R$ is $\le Ce^{-c(R/L)^a}$ and $|\Omega|=L^3(1+o(1))$ in (B.2). This is Clay-admissible.

**REJECTED sub-variant V2 (sparse separated packet array).** Take $P$ packets of fine scale $\ell=N^{-1}$ with centres separated by $\sim L\gg\ell$. Then $\mathrm{supp}\,u_p\cap\mathrm{supp}\,u_q=\emptyset$, so $(u_p\cdot\nabla)u_q\equiv0$ pointwise and the *only* coupling is far-field pressure, of relative strength
$$\frac{p_{\rm far}(L)}{p_{\rm loc}}\asymp\Big(\frac{\ell}{L}\Big)^3=(NL)^{-3}=e^{-3x}\longrightarrow0 .$$
The sparser the array, the more the packets decouple: **the intermittent-array reading of lens 3 is dynamically inert.** *(Exact failing equation: $p_{\rm far}/p_{\rm loc}=e^{-3x}\to0$.)* Kept: the band-filled blob (V2's $\phi\equiv1$ endpoint) is the only survivor, and its lacunarity is relative to the ambient domain, not internal.

---

## C. Scaling table ($\tau=T-t$; corner $\mu=\varepsilon=2/5$, $\delta=2/15$, $\gamma=8/15$)

| quantity | formula | general exponent | corner |
|---|---|---|---|
| energy $\|u\|_2^2$ | $\asymp4c_EL$ | $\mu-\varepsilon$ | $\tau^{0}$ (bounded) |
| front-octave energy $E_N=c_E/N$ | | $\gamma-\varepsilon$ | $\tau^{+2/15}$ |
| enstrophy $\|\omega\|_2^2$ | $\asymp N^2E_N$ | $-\gamma-\varepsilon$ | $\tau^{-14/15}$ |
| dissipation rate $\nu\|\nabla u\|_2^2$ | | $-\gamma-\varepsilon$ | $\tau^{-14/15}$, $\int dt<\infty$ ✓ |
| $\|u\|_{L^3}$ | $\asymp(2c_E)^{1/2}$ | $-\varepsilon/2$ | $\tau^{-1/5}\to\infty$ ✓ |
| $\|u\|_{L^\infty}\asymp U_L$ | $\sqrt{E_0}L^{-3/2}$ | $-\tfrac32\mu$ | $\tau^{-3/5}$ |
| $\sqrt\tau\|u\|_\infty$ (Type-I test) | | $\tfrac12-\tfrac32\mu$ | $\tau^{-1/10}\to\infty$ (Type II) |
| $\|\omega\|_{L^\infty}\asymp NU_N$ | | $-\tfrac{\gamma+\varepsilon+3\mu}{2}$ | $\tau^{-16/15}$; $I_{\rm BKM}\propto\tau^{-1/15}\to\infty$ |
| $\|\mathbb P(u\cdot\nabla)u\|_2$ | $\asymp NU_N\|u_N\|_2$ | | $\tau^{-1}$ |
| $\|\nabla p\|_2=\|(I-\mathbb P)(u\cdot\nabla)u\|_2$ | | | $\tau^{-1}$ |
| physical time remaining | $\tau$ | $1$ | $\tau$ |
| front parabolic clock $s=\int N^2dt$ | | $1-2\gamma$ | $\tau^{-1/15}\to\infty$ |
| Fourier bandwidth $N$ | (B.6) | $-\gamma$ | $\tau^{-8/15}$ |
| envelope wavenumber $L^{-1}$ | | $-\mu$ | $\tau^{-2/5}$ |
| **mismatch $NL$** | $(E_0/\nu^2L)^{1/3}$ | $-\delta$ | $\tau^{-2/15}\to\infty$ |
| active mode count $M^{\rm eff}$ | $\asymp N^3$ | $-3\gamma$ | $\tau^{-8/5}$ |
| independent dof $(NL)^3$ | | $-3\delta$ | $\tau^{-2/5}$ |
| octave count $J=\log_2 NL$ | | | $\tfrac2{15}\log_2(1/\tau)$ |

**Serrin sweep.** $\|u\|_{L^p}\propto\tau^{\frac25(3/p-3/2)}$, $q=2p/(p-3)$ $\Rightarrow$ $\int\|u\|_p^qdt$ has integrand exponent $-\frac65\frac{p-2}{p-3}<-1$ for **every** $3<p\le\infty$: all Serrin pairs diverge ✓. Vorticity pairs $2/q+3/p=2$: integrand exponent $-\frac{32p-36}{30p-45}<-1$ for all $3/2<p\le\infty$ ✓.

---

## D. Mechanism as a closed feedback loop

$$
\underbrace{L\downarrow}_{\dot L=-\kappa_L\sqrt{E_0}L^{-3/2}}
\ \xrightarrow{\ c_E=E_0/2L\ }\
\underbrace{c_E\uparrow}_{\text{supercritical shell amplitude}}
\ \xrightarrow{\ N=(2c_E/\nu^2)^{1/3}L^{-1}\ }\
\underbrace{N\uparrow\ \text{faster than }L^{-1}}_{NL=(E_0/\nu^2L)^{1/3}}
$$
$$
\xrightarrow{\ S_N/S_L=(NL)^{1/2}\ }
\underbrace{\text{front strain}\gg\text{envelope strain}}_{\text{(B.4)}}
\ \xrightarrow{\ S_L=U_L/L=\sqrt{E_0}L^{-5/2}\ }\
\underbrace{|\dot L|=\kappa_L U_L\ \text{sustained}}_{\text{closes the loop}}
$$
and the two exits:
$$\|u\|_{L^3}=(2c_E)^{1/2}\uparrow\infty\ \ (\text{ESS}),\qquad \|\omega\|_\infty=NU_N\propto\tau^{-16/15},\ \int\|\omega\|_\infty dt=\infty\ \ (\text{BKM}).$$
Every arrow is one of (B.1), (B.4), (B.6), (P1), (P5). The loop's gain is $\kappa_L$; it is a genuine positive feedback because $c_E\propto L^{-1}$ while the collapse rate $\propto c_E^{1/2}L^{-1}$ — i.e. contraction accelerates superlinearly, which is exactly the Type-II ($\varepsilon>0$) signature.

---

## E. Obstruction audit (exact collision points)

1. **Energy bound (Leray/F-N1).** Collision: (P1) $\varepsilon\le\mu$. We *saturate* it, $\varepsilon=\mu=2/5$, so $\|u\|_2^2\asymp4c_EL=2E_0$ is exactly constant. We never claim energy divergence.
2. **Finite dissipation (F-N2).** Collision: rate $\propto\tau^{-(\gamma+\varepsilon)}=\tau^{-14/15}$, $\int_0\tau^{-14/15}d\tau<\infty$. Margin $1/15$; had we taken the free-running front clock (P6) the exponent would be $\ge1$ and this constraint would fire — it is the constraint that selects viscous pinning (B.6).
3. **ESS $L^3$ (U-X1).** Collision: U-X2's identity $\|u\|_3^3=A^3L_r^2L_z\|U\|_3^3$ with $A=\sqrt{E_0}L^{-3/2}$, $L_r=L_z=L$ gives $A^3L^3=E_0^{3/2}L^{-3/2}\propto\tau^{-3/5}\ne1$. We escape U-X1 **not by spatial anisotropy but by non-parabolic amplitude** ($A\ne L^{-1}$, i.e. $\varepsilon>0$). Stated precisely because the repo's default reading of U-X2 is anisotropy-based.
4. **Fixed-finite-bandwidth no-go (F-α1, VR-L-011, Galerkin).** $N(t)\propto\tau^{-8/15}\to\infty$; $\gamma=8/15\in(0,1)$ satisfies VR-L-019. Not a fixed-band trajectory.
5. **Pure-swirl $L^3$ no-go (VR-L-016).** The blob is a full Leray-projected multi-octave 3D field, not $u^\theta e_\theta$; $P\not\equiv0$. **Design constraint imported:** by LG-11, the envelope generator must be **odd in $z$** about $x_*$ on any $z$-symmetric configuration, else the pressure channel vanishes identically. Recorded as a hard constraint on step 2 of §B.5.
6. **One-scale self-similar no-go (NRS1996 / Tsai1998).** Two independent hypothesis violations: (i) NRS requires $u=(2a\tau)^{-1/2}U(x/\sqrt{2a\tau})$ with amplitude exponent $1/2$; ours is $3\mu/2=3/5\ne1/2$. (ii) In similarity variables $y=x/L$, the profile's bandwidth is $NL\propto\tau^{-2/15}\to\infty$, so **no $\tau$-independent $U$ exists at all**. $\delta>0$ *is* the negation of self-similarity, so Tsai's finite-local-energy version (whose hypotheses we otherwise satisfy — our local energy is finite) cannot be applied.
7. **Galerkin global existence.** Same as 4.
8. **Smooth-forcing high-frequency decay (F-N4).** Not applicable: no forcing.
9. **Mesoscopic $\gamma<1$ empty-child no-go $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}/N^3$.** Collision point is *not* the mode count: by Lemma 1, $M^{\rm eff}\asymp N^3$ (relative band width $\eta\equiv1$, the $\gamma=1$ boundary the repo identified as the only clean case). The bound becomes $D_N\le2\kappa^2\tau^2c_E(\tau)$ with $c_E\propto\tau^{-2/5}\to\infty$: **the exclusion is vacuous in the limit**. This is why the repo's "required $c_E\approx228$" row is not a refutation here — $c_E$ is dynamical and financed by (P1) saturation. *Honest caveat:* this shows the no-go does not bind; it is **not** a proof that the flux is positive.
10. **Diagonal cross-talk gate.** Applies to sparse carrier alphabets with $M^{\rm eff}\ll N^3$, where sumset outputs land off-chain. At band filling, the sumset of $[L^{-1},N]$ with itself is $[0,2N]$; everything below $N$ is the intended band and everything above is the intended forward flux. There is no off-chain set. The gate is **dissolved, and traded for** the open positivity-of-net-flux question ($\chi_{\rm shape}>0$ in the S1/S3 sense) — that is what §F tests.
11. **CSTY Type-I exclusion.** (2.1): $\sqrt\tau\|u\|_\infty\propto\tau^{-1/10}\to\infty$ ✓. (2.2): at $r\asymp L$ the CSTY envelope is $r^{-1+\epsilon'}\tau^{-\epsilon'/2}\propto\tau^{-2/5-\epsilon'/10}$, exceeded by our $\tau^{-3/5}$ for all $\epsilon'<2$. **Audit item (open):** $\epsilon'\ge2$ and the exterior region $r\gg L$ must be checked globally, per the repo's explicit warning that local violation is insufficient.
12. **KNSS ancient Liouville.** Swirl is retained by construction (rules out (a)). For (b), $|u|\le C/r$ at $r\asymp L$ would need $LU_L=\sqrt{E_0}L^{-1/2}\le C$, but $LU_L\propto\tau^{-1/5}\to\infty$: the $C/r$ bound is violated with exponent $1/5$. Moreover the rescaled sequence has $NL\to\infty$, so it has no bounded-bandwidth ancient limit in any structure-preserving topology.
13. **Front-resolution threat model (TM-22).** $\ge7$ points per $N^{-1}$ ⟹ grid $\ge7NL$ points across the envelope diameter. Since $NL\propto\tau^{-2/15}$, doubling $NL$ costs a factor $2^{15/2}\approx181$ in $\tau$: a $128^3$–$512^3$ sweep covers $\sim$3 doublings of the mismatch. TM-04 (spectral pile-up) is the dominant confusion risk and is why the pinning check §F(iii) is mandatory.
14. **Route-12 material-transport bound.** $|\dot L|=\kappa_LU_L\le U_L$ ⟺ $\kappa_L\le1$ ✓ (this killed the vortex-necklace scaling; here it is satisfied by construction).

---

## F. Minimal falsification experiment ($\le1$ h)

**Object.** Build the §B.5 Gevrey band-filled blob on a periodic grid with $E_j=c_E/\lambda_j$ over $[L^{-1},N]$, exactly Leray-projected modewise.

**Variables.** $L\in\{1/4,1/8\}$; $NL\in\{4,8,16,32\}$; $\nu\in\{1/50,1/200\}$; $E_0=1$; grid $128^3$ (static diagnostics) and $256^3$ for the top two $NL$; one short RK4 run of duration $0.1/S_L$.

**Measurements & success criteria (pre-registered).**
- **(i) Strain ratio (Lemma 2).** Measure $\|\omega\|_\infty/(U_L/L)$ vs $NL$. *Success:* log-log slope $=1/2\pm0.10$ over $\ge4$ points, $R^2\ge0.97$. *Kill:* slope $\le0.2$.
- **(ii) Flux sign (the load-bearing claim).** $\Pi(k)=\langle P_{>k}\mathbb P(u\cdot\nabla)u,\,P_{>k}u\rangle$; normalized $\chi=\Pi(N/2)/(c_E^{3/2}N)$. *Success:* $\chi>0$ with a scale-independent floor across the $NL$ sweep. *Kill:* $\chi\le0$, or $\chi\propto(NL)^{-p}$ with $p>0$ at $R^2\ge0.9$.
- **(iii) Viscous pinning (B.6).** Locate $N_*$ where $\Pi(k)=\nu\int_0^k2q^2E(q)dq$. *Success:* $N_*L=(2c_E/\nu^2)^{1/3}$ within $30\%$ across both $\nu$. *Kill:* $N_*L$ independent of $\nu$ (would mean the pinning is a grid artifact).
- **(iv) Dynamic mismatch exponent (the decisive one).** From the RK4 run, $d\log(N_*L)/d\log(1/L)$. **Predicted $=1/3$.** *Kill the entire lens* if this is $\le0$.

**Arithmetic.** (i), (iii), (iv) may be binary64. **(ii) must not be:** the sign of a near-cancelling triple sum is exactly the repo's flagged failure mode (TM-09). Compute $\chi$ on a reduced exactly-representable mode set in **rational arithmetic**, and enclose the full-grid value by **interval arithmetic** with outward rounding; report the enclosure, not the float.

**Reuse.** `mesoscopic_local_fft.measure_local_fft_cloud` (exact zero-padded convolution, $K=4W-3$ padding — no wrap aliasing); `leray_response_relay.leray_project` / `leray_advection` for $\mathbb P((u\cdot\nabla)u)$; `mesoscopic_galerkin.run_small_mesoscopic_galerkin` as the RK4 harness; `mesoscopic_cloud_scaling` row schema for output; `exact_carrier_record_verifier` conventions for the independent rational recheck of (ii). New code: only the Gevrey multi-octave blob builder.

---

## G. Proof chain (10 obligations)

1. **Construction.** Existence of the Gevrey band-filled divergence-free family $u_{L,N}$ with exact LP profile $\|P_ju\|_2^2=2c_E/\lambda_j$ and $\|u\|_2^2=2E_0$. *(Elementary; formalizable.)*
2. **Flux positivity.** $\chi_{\rm shape}>0$: uniform positive lower bound on $\Pi(N/2)/(c_E^{3/2}N)$ over the family. **The single hardest step**; §F(ii) screens it.
3. **Strain-ratio lemma.** Rigorous two-sided version of (B.4) with constants, including the statement that envelope compression by front-scale strain is $O(e^{-x/2})$.
4. **Pinning lemma.** Rigorous version of (B.6): $N_*L\in[c_1,c_2](2c_E/\nu^2)^{1/3}$.
5. **Clock system.** Closure of (B.5)+(B.6)+(P1) into the ODE system with $L(\tau)\asymp\tau^{2/5}$, $NL\asymp\tau^{-2/15}$, with error terms controlled by 1–4.
6. **Invariant manifold.** The clock system's solution is an attracting invariant manifold of the front flow (S1) in the $s$-clock; forward RG integration, not fixed-point optimization.
7. **Entry (PO-09).** A smooth Schwartz datum whose trajectory enters the manifold. *Strategy-less in the repo; unchanged here.*
8. **Finite time.** $T-t(\tau_0)=\tau_0<\infty$ from the $L$-clock; requires 5 with rigorous rate bounds.
9. **Norm divergence.** $\limsup_{t\uparrow T}\|u(t)\|_{L^3}=\infty$ from $\|u\|_3\asymp(2c_E)^{1/2}\propto\tau^{-1/5}$, plus $\int\|\omega\|_\infty dt=\infty$; conclude non-extendability via ESS + Kato–Ponce.
10. **TARGET-U.** Assemble 1–9 into non-extendability of the maximal $H^m$ solution, $m>5/2$, with PO-12 (not a coordinate artifact) and PO-13 (interval certificate) discharged.

**Not claimed:** any of 2, 6, 7 is closed. Obligation 2 is the repo's standing open problem ("positive nonlinear flux margin") and this lens does not solve it — it relocates it from a sparse-carrier combinatorial question to a single-blob shape question, which §F can attack in one hour.
