# The Helicity Ledger: chirality-staggered carriers and a no-go for single-sign critical cascades

**LENS 7 — helicity-based closed subsystem, general Cartesian/Fourier. Full work-out of seed S5.**

**STATUS: SYMBOLIC CANDIDATE.**
Sub-results carry their own labels:
- **T1** (helical trilinear depletion identity) — **PROOF CANDIDATE**, finite-Fourier, Lean-able, numerically confirmed to $1.8\times10^{-14}$.
- **T2** (Helicity Ledger no-go for single-sign critical cascades) — **PROOF CANDIDATE** conditional on (H1)–(H4).
- **T3** ($\gamma=\tfrac12$ pinning) — **PROOF CANDIDATE** conditional on wake-persistence being the only route to $L^3$ divergence.
- **D1** (chirality-staggered carrier alphabet killing diagonal cross-talk exactly) — **SYMBOLIC CANDIDATE**.
- **R1** (helicity-lower-bounded palinstrophy budget) — **FORMAL ANSATZ**.

---

## A. Clay target

**Primary: CLAY-B** (global regularity, $\mathbb T^3$, no forcing) — *restricted to the critical-cascade ansatz class*. T2/T3 are conditional regularity statements: they empty a sub-class of the (B)-counterexample search space.
**Secondary: CLAY-D / TARGET-U pruning.** D1 is a design rule for the breakdown lane.

- **Domain**: $\mathbb T^3$, Fourier convention $u(x)=\sum_k\hat u_ke^{ik\cdot x}$, $\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$, $P_k=I-k\otimes k/|k|^2$. **No whole-space claim is made** (HS-5 bookkeeping is not crossed).
- **Forcing**: none. Admissibility moot. F-N1/F-N2/F-N4 therefore never invoked as a resource.
- **Initial regularity/decay**: $u_0$ a trigonometric polynomial (hence $C^\omega$, mean-zero, divergence-free) — inside VR-L-012's "finite-band datum, infinite-band trajectory" region, which VR-L-011 explicitly does **not** cover.
- **Viscosity**: $\nu>0$ fixed, never sent to $0$. Heat factor $h_{\nu,\tau}(r)=(1-e^{-\nu\tau r^2})/(\nu r^2)$ throughout.

---

## B. Central mathematics

### B.1 Helical basis and the exact helicity capacity

For $k\ne0$ pick $e_1(k)\perp k$ unit, $e_2(k)=\hat k\times e_1(k)$, and set $h^\pm(k)=e_1(k)\pm ie_2(k)$. Then
$$ik\times h^\pm(k)=i|k|(e_2\mp ie_1)=\pm|k|h^\pm(k),\qquad h^\pm\!\cdot\overline{h^\pm}=2,\quad h^+\!\cdot\overline{h^-}=0 .$$
Choosing $e_1(-k)=e_1(k)$ gives $e_2(-k)=-e_2(k)$, so $h^\pm(-k)=\overline{h^\pm(k)}$ and reality reads $a^\pm_{-k}=\overline{a^\pm_k}$ with $\hat u_k=a^+_kh^+(k)+a^-_kh^-(k)$, $a^\sigma_k=\tfrac12\hat u_k\!\cdot h^{-\sigma}(k)$.

Energy and helicity ($H=\int u\cdot\omega$, $\hat\omega_k=ik\times\hat u_k$):
$$E=(2\pi)^3\!\sum_k\big(|a^+_k|^2+|a^-_k|^2\big),\qquad H=2(2\pi)^3\!\sum_k|k|\big(|a^+_k|^2-|a^-_k|^2\big).$$
Restricting to a shell supported in $|k|\in[N_-,N_+]$, $N_\pm=N(1\pm\eta)$:
$$\boxed{\;|H_N|\le 2N_+E_N\;}\qquad\text{(seed's }|H_N|\le2NE_N\text{, sharpened by the width).}$$

**Dimension check.** $[H]=L^4T^{-2}$, $[NE]=L^{-1}\cdot L^5T^{-2}=L^4T^{-2}$. ✓
**Scaling check.** Under $u_\lambda=\lambda u(\lambda x,\lambda^2t)$: $\omega_\lambda=\lambda^2\omega(\lambda x)$, so $H_\lambda=\lambda^3\lambda^{-3}H=H$ — **helicity is exactly scale-invariant (critical)**, while $E_\lambda=\lambda^{-1}E$. Hence at the repo's critical normalization $E_N=c_E/N$,
$$|H_N|\le 2(1+\eta)c_E \qquad\textbf{— scale-free capacity.}$$
Define the **helicity code** $h_j:=H_{N_j}/(2N_jE_{N_j})\in[-1,1]$ of stage $j$, so $H_{N_j}=2c_Eh_j$ when $E_{N_j}=c_E/N_j$.

### B.2 (i) Linear polarization $\Rightarrow$ exactly zero helicity

Let $u=A\cos(k\cdot x+\varphi)$, $A\in\mathbb R^3$, $A\cdot k=0$; write $A=\alpha e_1+\beta e_2$, $\alpha,\beta\in\mathbb R$. Then $\hat u_k=\tfrac12Ae^{i\varphi}$ and
$$a^+_k=\tfrac12\hat u_k\!\cdot h^-=\tfrac{e^{i\varphi}}4(\alpha-i\beta),\qquad a^-_k=\tfrac{e^{i\varphi}}4(\alpha+i\beta)\;\Rightarrow\;|a^+_k|=|a^-_k|=\tfrac{\sqrt{\alpha^2+\beta^2}}4 .$$
So $H_k=0$ **exactly, mode by mode**. More generally $\hat u_k=c_kA_k$ with $A_k$ real forces $|a^+_k|=|a^-_k|$.

> **Consequence for the repo.** Every carrier in `exact_carrier_search`, `expanded_carrier_search` and `carrier_two_stage_galerkin` uses *primitive integer real polarizations*, i.e. is linearly polarized. **The entire registered relay alphabet has identically zero helicity spectrum, mode by mode ($p=1/2$ everywhere).** Verified numerically: $H=0.0$ exactly for $(1,-1,0)\cos(x+y)$.

### B.3 (ii) The exact trilinear depletion identity (T1)

Use $(u\cdot\nabla)u=\nabla\tfrac{|u|^2}2+\omega\times u$, so $\mathcal N_k=-P_k\sum_{\ell+m=k}\hat\omega_\ell\times\hat u_m$. Symmetrising the ordered pair sum,
$$\hat\omega_\ell\times\hat u_m+\hat\omega_m\times\hat u_\ell=\sum_{s,s'}a^s_\ell a^{s'}_m\Big[s|\ell|\,h^s(\ell)\times h^{s'}(m)+s'|m|\,h^{s'}(m)\times h^{s}(\ell)\Big],$$
and since $X\times Y=-Y\times X$,
$$\boxed{\;\mathcal G^\sigma_k:=\tfrac12\mathcal N_k\!\cdot h^{-\sigma}(k)=-\tfrac14\!\!\sum_{\ell+m=k}\sum_{s,s'}\big(s|\ell|-s'|m|\big)\,a^s_\ell a^{s'}_m\,\big[h^s(\ell)\times h^{s'}(m)\big]\!\cdot h^{-\sigma}(k).\;}\tag{T1}$$
$P_k$ drops because $h^{-\sigma}(k)\perp k$. Units: $[a]^2[|\ell|]=L^2T^{-2}L^{-1}=$ velocity$^2$/length $=\partial_t$(velocity). ✓

**Vanishing rule.** The coefficient $s|\ell|-s'|m|$ is zero **iff $s=s'$ and $|\ell|=|m|$**. Therefore:
- *(Beltrami depletion, derived, not assumed)* a single-shell field with all $a^-=0$ has $\omega=Nu$, every pair has $s|\ell|-s'|m|=0$, hence $\mathcal N\equiv0$: **extremal helicity $\Rightarrow$ zero Leray transfer**.
- *(stronger)* the depletion is **pairwise**: any two modes of the same helicity sign on the same sphere contribute exactly nothing, regardless of the rest of the field.

**Numerical confirmation** (32³ pseudo-spectral, $k_A=(1,1,0)$, $k_B=(1,0,1)$, $|k|^2=2$ both):

| configuration | $\max|P((u\cdot\nabla)u)|$ |
|---|---|
| $s_A=s_B=+$, $|k_A|=|k_B|$ | $1.81\times10^{-14}$ (machine zero) |
| $s_A=+,s_B=-$, $|k_A|=|k_B|$ | $1.059$ |
| $s_A=s_B=+$, $k_B=(2,1,0)$ | $0.438$ |

and $H/(2|k|E)=1.0000000000$ for a single circular mode. All four predictions of (T1) confirmed.

### B.4 Chirality-refined Duhamel bound (sharpening L-11)

In (T1) the *single* derivative factor of the repo's Bernstein step is exactly $s|\ell|-s'|m|$. For a cloud in $[N_-,N_+]$: same-sign pairs give $\big||\ell|-|m|\big|\le2\eta N$; opposite-sign pairs give $|\ell|+|m|\le2(1+\eta)N$. With $\ell^1$ amplitudes $A^\pm=\sum_k|a^\pm_k|$ and imbalance $p=A^+/(A^++A^-)$, define
$$\frac{\kappa^{\rm hel}_N}{N}:=2\eta\big(p^2+(1-p)^2\big)+4(1+\eta)\,p(1-p).$$
Replacing $\kappa$ in the repo's exact identity gives
$$\boxed{\;D_N=\frac{E_{\rm child}(\tau N^{-2})}{E_{\rm parent}(0)}\;\le\;2\tau^2c_E\,\Big[\min\Big(\kappa,\tfrac{\kappa^{\rm hel}_N}{N}\Big)\Big]^2\frac{M^{\rm eff}_N}{N^3}.\;}\tag{B.4}$$
At $\eta=0.2$: $p=1$ (pure chirality) gives $0.4$; $p=\tfrac12$ gives $1.4$. **Chirality-pure clouds lose a factor $(1.4/0.4)^2=12.25$ in $D_N$** — the worst possible carriers. This is strictly sharper than the phase-independent bound whenever $p\notin[0.3,0.7]$ and reduces to it otherwise (take the min).

### B.5 (iii) The Helicity Ledger (T2)

Spectrally, $\dot H|_{\rm visc}=-2\nu\sum_k|k|^2H_k$; equivalently $\dot H=-2\nu\int\omega\cdot(\nabla\times\omega)$ (check: $\widehat{\nabla\times\omega}$ acts as $|k|^2$ on each helical channel, giving $\int\omega\cdot\nabla\times\omega=\sum_k|k|^2H_k$ ✓). The nonlinearity conserves $H$ exactly. With stages $N_j=2^jN_0$, $E_j(t)=g_j(t)c_E/N_j$, $H_j(t)=2c_Eh_jg_j(t)$:
$$\mathcal H(t)=\mathcal H(0)-4\nu c_E\!\int_0^t\!\sum_{j\le J(s)}N_j^2h_jg_j(s)\,ds. \tag{L}$$

**Wake retention.** The heat semigroup acts as the *same positive real scalar* $e^{-\nu|k|^2t}$ on both helical channels at each $k$, so the pointwise code $(|a^+_k|^2-|a^-_k|^2)/(|a^+_k|^2+|a^-_k|^2)$ is **exactly invariant** under the linear part; only $g_j$ decays. With $N(\tau)\asymp\tau^{-\gamma}$, stage $j$ is born at $\tau_j\asymp N_j^{-1/\gamma}$ and its residual retention is
$$g_j(T^-)\asymp\exp\!\big(-2\nu N_j^2\tau_j\big)=\exp\!\big(-2\nu N_j^{\,2-1/\gamma}\big)=\begin{cases}\to1,&\gamma<1/2\\ e^{-\nu/a}=:g_\infty>0,&\gamma=1/2\\ \to0\ \text{super-exponentially},&\gamma>1/2.\end{cases}$$

**Theorem T2 (single-sign helical critical cascade no-go).** Assume
(H1) octave stages at critical energy $E_j\ge g_\infty c_E/N_j$, $g_\infty>0$;
(H2) $H_j=2N_jE_jh_j$ with $|h_j|\le1$;
(H3) the helicity spectrum is sign-definite, $H_k\ge0$ for all $k,t$;
(H4) $h_j\ge h_*>0$ for all $j$.
Then $\dot{\mathcal H}=-2\nu\sum_k|k|^2H_k\le0$, so $\mathcal H(t)\le\mathcal H(0)$, while (H1)+(H4) force $\mathcal H(t)\ge2c_Eg_\infty h_*(J(t)+1)$. Hence
$$\boxed{\;J(t)+1\;\le\;\frac{\mathcal H(0)}{2c_Eg_\infty h_*}=\frac{e^{\nu/a}\,\mathcal H(0)}{2c_Eh_*}\quad\text{for all }t<T.\;}$$
So $N(t)\le2^{J_{\max}}N_0$ is **bounded**, the trajectory stays in a fixed finite band, and VR-L-011 / F-$\alpha$1 then forbids breakdown. ∎

Since $\mathcal H(0)\le2c_E(1+\eta)(J_0+1)$ for an initially $J_0$-octave datum, the bound reads $J_{\max}\lesssim e^{\nu/a}(J_0+1)/h_*$: **a single-sign cascade can create only $O(1/h_*)$ new octaves per initial octave.**

**Corollary (design rule).** An unbounded critical cascade must satisfy $\big|\sum_{j\le J}h_j\big|\le\mathcal H(0)/(2c_Eg_\infty)$ for every $J$, hence
$$\Big|\frac1{J+1}\sum_{j\le J}h_j\Big|\le\frac{\mathcal H(0)}{2c_Eg_\infty(J+1)}\longrightarrow0 :$$
**the helicity code must be mean-zero — alternating, or identically zero.** For $\gamma>1/2$ the weights $\kappa_j\asymp N_j^{2-1/\gamma}$ grow geometrically and the same argument tightens to $|h_J|\le C\,2^{-(2-1/\gamma)J}$.

### B.6 (T3) The exponent $\gamma$ is pinned to $1/2$

Three independent budgets:
- **BKM** $\int_0^T\|\omega\|_\infty dt=\infty$ with $\|\omega\|_\infty\asymp N^2\asymp\tau^{-2\gamma}$ $\Rightarrow\gamma\ge1/2$.
- **Finite dissipation** $\int2\nu c_EN\,dt<\infty\Rightarrow\gamma<1$.
- **Wake persistence** (§B.5): $g_\infty>0\Rightarrow\gamma\le1/2$. Since $\|u\|_3^3$ receives $O(1)$ per *retained* octave and $O(1)$ total from the live front alone, ESS's required $\|u\|_3\to\infty$ (via S4) is unreachable without a persistent wake.

$$\Rightarrow\quad\boxed{\gamma=\tfrac12\ \text{exactly}}$$
i.e. $N(t)=(2a(T-t))^{-1/2}$, matching S1/STATUS's $\dot N=kN^3$. This sits **exactly on the Type-I boundary** $\sqrt\tau\|u\|_\infty=O(1)$ — see §E.

### B.7 (D1) Chirality-staggered carriers: exact cross-talk annihilation

The registered failure: two same-scale relays with parents $A_1,A_2$ and $B_1,B_2$ (all with $|k|^2=2$) have **diagonal** pairs $A_1{+}B_2$, $A_2{+}B_1$ landing on the intended child shell with $\|B_{\rm cross}\|_2^2=2483/1890$ against $\|B_{\rm int}\|_2^2=37/315$ — intended is $8.94\%$ of cross-talk, and **all 16 real orientations fail**.

(T1) explains this: real polarizations carry *both* signs at every $k$ ($p=1/2$), so the same-shell diagonal pairs always have a nonzero opposite-sign component with coefficient $|\ell|+|m|=2\sqrt2$.

**The fix.** Use *circularly* polarized carriers $u=\mathrm{Re}\,\big(a\,h^{s}(k)e^{ik\cdot x}\big)$ (real fields; $\hat u_{\pm k}$ are pure-$s$) and 2-colour the parent set so that **cross-talk edges are monochromatic and intended edges are bichromatic**:
$$s(A_1)=s(B_2)=+,\qquad s(A_2)=s(B_1)=-.$$
Since all four have $|k|^2=2$:
- $A_1{+}B_2$: $s=s'$, $|\ell|=|m|$ $\Rightarrow$ coefficient $=0$ **exactly**;
- $A_2{+}B_1$: same $\Rightarrow$ $0$ **exactly**;
- $A_1{+}A_2$ and $B_1{+}B_2$: $s\ne s'$ $\Rightarrow$ coefficient $=|\ell|+|m|=2\sqrt2$, **maximal**.

The bipartition $\{A_1,B_2\}\mid\{A_2,B_1\}$ exists precisely because the intended-edge graph is bipartite with cross-edges inside the parts. Shell helicity: $H\propto|a_{A_1}|^2+|a_{B_2}|^2-|a_{A_2}|^2-|a_{B_1}|^2=0$ at equal amplitudes, so $h_j=0$ — **ledger-neutral (T2 satisfied) and $p=1/2$ in aggregate (B.4 optimal), while every individual mode is chirality-pure.** The three constraints are simultaneously satisfiable.

**Cloud version.** At width $\eta N$ the same-sign suppression is $2\eta N$ vs $2(1+\eta)N$, amplitude factor $\eta/(1+\eta)=1/6$ at $\eta=0.2$, power factor $1/36$. Applied to the registered numbers: cross $2483/36=69$ vs intended $222$ — **ratio $3.2$ in favour of the intended channel, a gate flip.** (Prediction; the pilot must measure whether the intended coefficient is also modified.)

### B.8 REJECTED sub-variants (kept, with the exact failing relation)

- **V1 REJECTED — "Cauchy–Schwarz forbids single-sign accumulation."** Needed $2c_E(J+1)>2\sqrt{E\Omega}=4c_E\sqrt{N_J/N_0}=4c_E2^{J/2}$. **False for every $J\ge0$** ($\log\ll\sqrt{\ }$). The seed's suggested route is vacuous; only the *dissipation* ledger (L) works.
- **V2 REJECTED — "shell chirality balance $p=1/2$ suffices for transfer."** Pure swirl $u=u^\theta(r,z)e_\theta$ has $\omega=-\partial_zu^\theta e_r+r^{-1}\partial_r(ru^\theta)e_z$, hence $u\cdot\omega\equiv0$ *pointwise*, hence $p=1/2$ in every shell — yet $J(u_0)=F'(0)\le0$ (VR-L-016). Balance is necessary-ish, never sufficient.
- **V3 REJECTED — "go maximally helical to exploit depletion of the inverse channel."** (T1)'s coefficient is symmetric in $\ell\leftrightarrow m$; forward transfer is depleted equally. (B.4): $p=1$ costs a factor $12.25$ in $D_N$ at $\eta=0.2$.
- **V4 REJECTED — "$|H(t)|\to\infty$ as the blow-up signature."** $H$ *is* scale-critical (§B.1), so this is dimensionally legitimate, but the mirror-antisymmetric class $u(Rx)=-Ru(x)$ is NS-invariant and has $H\equiv0$; a singularity there would have $H$ bounded. So $H$-divergence cannot be **necessary**. (This is also why T2 only prunes and can never close (B).)

---

## C. Scaling table ($\tau=T-t$, $\gamma=1/2$, $N=(2a\tau)^{-1/2}$, $E_j=c_E/N_j$)

| quantity | law | $\tau$-exponent |
|---|---|---|
| energy $E$ | $\le2c_E/N_0$ | $0$ (bounded) |
| enstrophy $\Omega=\sum N_j^2E_j\asymp2c_EN$ | $\asymp\tau^{-1/2}$ | $-1/2$ |
| global $\|u\|_3^3\asymp g_\infty^{3/2}(J{+}1)$ | $\asymp\tfrac12\log_2(1/\tau)$ | $0^-$ (**log-divergent**) |
| $\|\omega\|_\infty\asymp N^2$ | $\asymp\tau^{-1}$ | $-1$ |
| dissipation rate $2\nu\Omega$ | $\asymp\tau^{-1/2}$ | $-1/2$ ($\int<\infty$ ✓) |
| nonlinear $\|\mathbb PB(u,u)\|_2\asymp N^2\|u_N\|_2$ | $\asymp\tau^{-3/4}$ | $-3/4$ |
| pressure $\|\nabla p\|_2\le\|(u{\cdot}\nabla)u\|_2$ | $\asymp\tau^{-3/4}$ | $-3/4$ |
| physical time remaining $\tau$ | $\asymp N^{-2}$ | $1$ |
| Fourier bandwidth $N$ | $\asymp\tau^{-1/2}$ | $-1/2$ |
| active mode count $M_N\asymp\eta^3N^3$ | $\asymp\eta^3\tau^{-3/2}$ | $-3/2$ |
| wake mode count $\tfrac87\eta^3N^3$ | $\asymp\tau^{-3/2}$ | $-3/2$ |
| **helicity $\mathcal H=2c_Eg_\infty\sum h_j$** | **$O(1)$ by (T2)** | $0$ |
| helicity dissipation $|\dot{\mathcal H}|\asymp4\nu c_EN^2|h_J|$ | $\asymp|h_J|\tau^{-1}$ | $-1$ |
| palinstrophy $\|\nabla\omega\|_2^2\asymp2c_EN^3$ | $\asymp\tau^{-3/2}$ | $-3/2$ |
| $\mathcal C(t)=\int\|\nabla\omega\|_2^2$ | $\asymp\tau^{-1/2}$ | $\to\infty$ ✓ |
| Type-I indicator $\sqrt\tau\|u\|_\infty$ | $\asymp\tau^{1/2}\tau^{-1/2}$ | $0$ (**boundary**) |

---

## D. Closed feedback loop (every arrow a formula)

$$\textbf{S1 (stage $j$)}\ \ u_{N_j},\ E_j=\tfrac{c_E}{N_j},\ h_j=0,\ p=\tfrac12\ \text{aggregate, modes chirality-pure}$$
$$\xrightarrow[\ \text{(T1) opposite-sign pairs}\ ]{\ \mathcal G^\sigma_{2k}=-\tfrac14(|\ell|+|m|)a^+_\ell a^-_m[h^+\!\times h^-]\!\cdot h^{-\sigma},\ \ |\ell|+|m|=2N_j\ }\textbf{S2: forcing } f_{2N_j}$$
$$\xrightarrow[\ \text{Duhamel}\ ]{\ \hat v_k=N^{-2}h_{\nu,\tau}(|k|/N)\hat f_k,\quad D_N\le2\tau^2c_E\big[\min(\kappa,\kappa^{\rm hel}_N/N)\big]^2M^{\rm eff}_N/N^3\ }\textbf{S3: child at }2N_j$$
$$\xrightarrow[\ \text{cross-talk annihilation (D1)}\ ]{\ s(A_1)=s(B_2),\ |k_{A_1}|=|k_{B_2}|\ \Rightarrow\ (s|\ell|-s'|m|)=0\ \text{exactly}\ }\textbf{S4: clean child}$$
$$\xrightarrow[\ \text{renormalise}\ ]{\ E_{j+1}:=D_NE_j\stackrel{!}{=}\tfrac{c_E}{N_{j+1}},\ \ h_{j+1}=0\ \text{(colour negated)}\ }\textbf{S1 at }N_{j+1}=2N_j$$

Closure conditions, all formulas:
1. **Flux closure** $D_N\ge\tfrac12$ uniformly in $N$ (S3's $c_E$-collapse: any $\chi_{\rm shape}>0$ closes it at $c_E\ge2\nu^2/\chi_{\rm shape}^2$).
2. **Ledger closure** $\sum_{j\le J}h_j=0$ for all $J$ — satisfied *identically* by the staggered design, $h_j\equiv0$.
3. **Front clock** $\dot N=aN^3\Rightarrow N^{-2}(t)=N_0^{-2}-2at\Rightarrow T=N_0^{-2}/(2a)<\infty$.
4. **Wake retention** $g_\infty=e^{-\nu/a}>0$ (heat multiplier is helicity-blind, §B.5).
5. **Critical-norm divergence** $\|u\|_3^3\ge g_\infty^{3/2}\,c\,(J+1)=g_\infty^{3/2}c\log_2\!\big(N(t)/N_0\big)\to\infty$ (S4).

The loop is *closed*: helicity neutrality (2) is not an extra cost, it is a property of the alphabet chosen for (1); and (4) is the same fact that makes (5) work.

---

## E. Obstruction audit — exact collision points

1. **Energy bound (F-N1/N-2).** $E=\sum_jg_jc_E/N_j\le2c_E/N_0$. Bounded by construction; energy is **never** used as the blow-up signature. No collision.
2. **Finite dissipation (F-N2).** $\int_0^T2\nu c_EN\,dt=2\nu c_E\int_0^{T}(2a\tau)^{-1/2}d\tau=2\nu c_E\sqrt{2T/a}<\infty$. Collision point: $\gamma<1$; we have $\gamma=1/2$. ✓
3. **ESS $L^3$ (U-X1, VR-N-002).** ESS needs $\|u\|_{L^\infty_tL^3}=\infty$. Collision point is *exactly* the wake-retention lemma: if $g_\infty=0$, $\|u\|_3^3=O(1)$ and ESS kills the candidate. T3 shows $g_\infty>0\Leftrightarrow\gamma\le1/2$, and BKM forces $\gamma\ge1/2$ — the mechanism survives on a **single point** $\gamma=1/2$, nowhere else. This is the tightest collision in the document.
4. **Fixed-finite-bandwidth no-go (F-$\alpha$1, VR-L-011/012).** $N(t)=(2a\tau)^{-1/2}\to\infty$: not a fixed band. Datum is finite-band, trajectory is not — VR-L-012's explicitly uncovered region. Note T2 attacks candidates *by forcing them back into* this no-go.
5. **Pure-swirl $L^3$ no-go (VR-L-016).** Collision point: the dead class is exactly $\{u\cdot\omega\equiv0$ pointwise$\}$ (§B.8 V2). Our carriers have $u\cdot\omega\ne0$ pointwise ($h^\pm$ modes are individually extremal) while summing to $H_N=0$. The design is *pointwise chiral, aggregately neutral* — the complement of the dead class, and not Beltrami either (V3).
6. **One-scale self-similar (NRS / Tsai).** NRS needs $u=(2a\tau)^{-1/2}U(x/\sqrt{2a\tau})$ with $U\in L^3(\mathbb R^3)$. Collision point: our profile has $\|u\|_3^3\asymp\log(1/\tau)$, so no time-independent $U$ with $\|U\|_3<\infty$ exists — the hypothesis fails at the $L^3$ clause, not at the similarity clause. Tsai's $W^{1,2}_{\rm loc}\cap L^q$, $3<q<\infty$ version fails identically. Additionally the object is *discretely* self-similar with factor $2$ in $N$ (period $S=\log2/a$ in the RG clock $s$), and if the code alternates, factor $4$. Seregin's DSS Liouville filters must be re-checked against the period-$2S$ orbit — **flagged open**, not claimed evaded.
7. **Galerkin global existence.** Same as (4). Also T2's conclusion routes a single-sign cascade *into* this theorem, which is what makes T2 a genuine no-go rather than a heuristic.
8. **Smooth-forcing high-frequency decay (F-N4).** No forcing used. Vacuous.
9. **Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}/N^3$.** Collision point: relative width. We use the **only** surviving family, fixed-relative $W_N=\lfloor\eta N\rfloor$, $\eta\in(0,1/3)$, giving $M^{\rm eff}\asymp\eta^3N^3$ and $D_N\le2\tau^2c_E(1+2\eta)^2C\eta^3=O(1)$ in $N$ — no exponent obstruction. **Unresolved at the constant level**: the registered $\eta=0.2,N=64$ row needs $c_E\approx228$; (B.4) shows chirality-*pure* clouds would need $\approx1425$, chirality-balanced clouds are optimal. Honest status: exponent-clean, constant-open.
10. **Diagonal cross-talk gate.** The head-on collision. Registered: intended/cross $=0.0894$, 16/16 orientations fail. Collision point identified precisely: *real polarization forces $p=1/2$ at every $k$, so opposite-sign components are always present on the same-shell diagonal pairs.* Evasion: 2-colouring by helicity sign, $s(A_1)=s(B_2)$, $s(A_2)=s(B_1)$, all $|k|^2=2$ $\Rightarrow$ $(s|\ell|-s'|m|)=0$ **identically**, not approximately. Caveat: the *off*-diagonal mixed pairs $A_1{+}B_1$, $A_2{+}B_2$ remain bichromatic and must still be geometrically off-target — this must be verified, not assumed.
11. **CSTY Type-I exclusion.** CSTY2009 hypotheses require **axisymmetry**. Our object is a general Cartesian/Fourier chiral gadget with no rotational symmetry (indeed a chirality-staggered field cannot be axisymmetric-with-swirl-only, since pure swirl has $u\cdot\omega\equiv0$). CSTY does not apply. This matters, because T3 pins $\gamma=1/2$, i.e. $\sqrt\tau\|u\|_\infty=O(1)$ — **exactly Type I**. Under axisymmetry the candidate would be dead; in the general class Type-I exclusion is open. This is a real, named risk, stated as such.
12. **KNSS ancient-solution Liouville.** (a) needs axisymmetry-without-swirl; (b) needs axisymmetry-with-swirl plus $|u|\le C/r$. Neither applies (no axisymmetry). Our rescaled limit is the S1 front-flow $s$-periodic orbit $\partial_s\Psi=a(2\Psi+\xi\cdot\nabla_\xi\Psi)-\nu|\xi|^2\Psi-\mathcal Q(\Psi,\Psi)$, a *time-periodic* object in the RG clock, not a bounded ancient solution in $t$. Collision point if one later imposes axisymmetry: the whole lane dies at (5) anyway.
13. **Front-resolution threat model (TM-01/03/22).** The decisive pilot stage is **exact algebraic** (field $\mathbb Q(i,\sqrt2)$), so there is no grid, no aliasing, and no points-per-front to satisfy. The FFT confirmation stage uses `mesoscopic_local_fft`'s zero-padded linear convolution ($K=4W-3$), which is exact for the stated support — TM-03 structurally excluded, TM-22 inapplicable (no front is fitted).
14. **Kinematic incompatibility (Z-01/Z-02).** We never use compact-support packets; all objects are exact trigonometric polynomials on $\mathbb T^3$. Vacuous. Correspondingly **no whole-space claim** is made (HS-5 untouched).
15. **Rank-one moment infeasibility.** Not used; our carriers are two-mode chiral pairs with full-rank moment tensors.

---

## F. Minimal falsification experiment ($\le1$ hour)

**Stage 1 (exact, decisive, $\sim$10 min).** Extend `exact_carrier_search.py` / `expanded_carrier_search.py` with a helicity column.
- *Variables*: the 4 known parents ($|k|^2=2$), $2^4=16$ helicity assignments $\times$ the 16 registered orientations $=256$ configurations.
- *Arithmetic*: **exact, required**. For $k=(1,1,0)$ take $e_1=(1,-1,0)/\sqrt2$, $e_2=\hat k\times e_1=(0,0,-1)$ — so $h^\pm\in\mathbb Q(\sqrt2)^3+i\,\mathbb Q^3$. Clearing $\sqrt2$ puts everything in $\mathbb Z[i]$. Use `sympy` or a 20-line $\mathbb Q(i,\sqrt2)$ class. Float is *not* admissible for the "exactly zero" claim.
- *Measured*: (a) Waleffe coefficient $s|\ell|-s'|m|$ per pair; (b) $\Pi_{\rm intended}$; (c) $\|B_{\rm cross}\|_2^2$ on the target child shell; (d) per-shell $H_N$ and code $h$.
- **Success**: $\ge1$ assignment with $\Pi_{\rm intended}\ne0$ exactly, $\|B_{\rm cross}\|_2^2=0$ exactly, $H_N=0$ exactly. (Prediction: the colouring $s(A_1)=s(B_2)=+$, $s(A_2)=s(B_1)=-$ does this.)
- **Kill**: all 256 have $\Pi_{\rm intended}=0$ or $\|B_{\rm cross}\|^2\ne0$ $\Rightarrow$ (T1) mis-derived or the geometry does not bipartition $\Rightarrow$ D1 dead (T2/T3 survive independently).
- *Independent re-verification*: extend `exact_carrier_record_verifier.py` (different code path: ordered-pair convolution $\hat B(n)=i\sum_{k+\ell=n}(\hat u(k)\cdot\ell)\hat v(\ell)$). Mandatory — TM-14.

**Stage 2 (float corroboration, $\sim$5 min).** Already run here at $32^3$: same-sign/same-modulus pair gives $\max|\mathbb PB|=1.81\times10^{-14}$. Reproduce with `leray_response_relay.leray_advection`, tolerance $10^{-13}$. Float acceptable.

**Stage 3 (cloud, $\sim$30 min).** `mesoscopic_local_fft.measure_local_fft_cloud` at $N=64$, $\eta=0.2$ ($W=13$, padding $K=4W-3=49$, grid $128^3$), parent built from `leray_response_relay.helical_fejer_packet` (already exists) in a chirality-staggered configuration vs. the linear-polarization baseline. Float acceptable.
- **Success**: cross-channel power suppressed by $\ge10\times$ relative to baseline (predicted $\approx36\times$) with intended-channel power reduced by $<3\times$.
- **Kill**: intended channel drops by the same factor as cross-talk $\Rightarrow$ the depletion is not selective at cloud width $\Rightarrow$ D1 degrades to a constant-factor rule only.

**Not tested by this pilot** (state explicitly): the flux floor $D_N\ge\tfrac12$, any orbit, wake retention $g_\infty$, $\|u\|_3$ divergence, or any PDE statement. This is a *cross-talk-gate* experiment only.

---

## G. Proof chain if it works ($\le10$ obligations)

1. **T1 in Lean**: $\mathcal N_k\cdot h^{-\sigma}(k)=-\tfrac12\sum(s|\ell|-s'|m|)a^s_\ell a^{s'}_m[h^s\times h^{s'}]\cdot h^{-\sigma}$ as sorry-free finite algebra over a finite mode set, in the style of `MesoscopicDuhamelNoGo.lean`.
2. **Capacity lemma in Lean**: $|H_N|\le2N_+E_N$ from the helical Cauchy–Schwarz; hence $|h_j|\le1$.
3. **T2 in Lean**: from (2) plus $\dot{\mathcal H}=-2\nu\sum|k|^2H_k$ and $H_k\ge0$, derive $J+1\le\mathcal H(0)/(2c_Eg_\infty h_*)$. Finite nonnegative-real algebra; no PDE bridge needed.
4. **(B.4) in Lean**: chirality-refined Duhamel bound $D_N\le2\tau^2c_E[\min(\kappa,\kappa^{\rm hel}_N/N)]^2M^{\rm eff}_N/N^3$ — extends the existing `MesoscopicDuhamelNoGo` inequality.
5. **Exact-zero cross-talk certificate** in $\mathbb Q(i,\sqrt2)$ for the staggered gadget, with independent re-verification (Stage 1 above). *Certificate-verified*, not Lean.
6. **Flux floor**: exhibit $(\eta,c_E,\tau,\nu)$ with $D_N\ge\tfrac12$ uniformly in $N$ for the staggered cloud. **Currently open — the single largest gap**, shared with the whole relay lane.
7. **Two-stage recurrence**: the staggered gadget's grandchildren reproduce a scaled copy with negated colouring, i.e. $\mathfrak T^2$ has a fixed point. Test by forward RG integration (S7), then rigorously.
8. **Wake persistence**: $g_j(T^-)\ge g_\infty=e^{-\nu/a}>0$ at $\gamma=1/2$, with the helicity code exactly preserved by the heat multiplier.
9. **Flux $\to L^3$** (S4): sustained $\Pi_N\ge qc_EN$ plus (8) gives $\|u(t)\|_3^3\ge c\log N(t)\to\infty$, discharging the ESS-required divergence.
10. **Finite-Fourier $\to$ PDE bridge**: spectral-cutoff limit $T'\uparrow T$, torus Sobolev embeddings, Kato local uniqueness (F-$\delta$). **Explicitly open**; nothing above is a PDE theorem without it.

---

**Non-claims.** Nothing here proves a finite-time singularity, global regularity, or any Clay statement. T2 is an exclusion inside a stated ansatz class, not a regularity theorem. D1 is an algebraic design rule whose decisive premise (obligation 6) is untested. The numerical values in §B.3 are float observations of an exact algebraic identity, not evidence about the PDE.
