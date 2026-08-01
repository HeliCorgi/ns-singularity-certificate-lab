# LENS 10 — The Fourier Delocalization Number $\mu_N$: a critical quantity invisible to energy methods

**Status: PROOF CANDIDATE** for the dichotomy theorem (§B.1–B.5, §C) — elementary, finite-Fourier, Lean-able.
**FORMAL ANSATZ** for the blow-up feedback loop (§D). Sub-variants that died in derivation are kept in §H.

Conventions: $\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3$, normalized measure $dx/(2\pi)^3$ (so $\|u\|_2^2=\sum_k|\hat u_k|^2$, $\|u\|_p\le\|u\|_\infty$), $u=\sum_k\hat u_ke^{ik\cdot x}$, $\hat u_0=0$, $k\cdot\hat u_k=0$,
$$\partial_t\hat u_k=\mathcal N_k-\nu|k|^2\hat u_k,\qquad \mathcal N_k=-iP_k\!\!\sum_{\ell+m=k}\!(m\cdot\hat u_\ell)\hat u_m,\qquad P_k=I-\tfrac{k\otimes k}{|k|^2}.$$
Shells $S_j=\{2^j\le|k|<2^{j+1}\}$, $N=2^j$, $P_ju=\sum_{S_j}\hat u_ke^{ikx}$, $a_k=|\hat u_k|$, $W_j=\sum_{S_j}a_k$. Wiener/Fourier–Lebesgue norms $\|u\|_{X^s}=\sum_k|k|^s a_k$ ($X^{-1}$ = Lei–Lin space). Lengths are dimensionless, so $[u]=[\nu]=T^{-1}$.

---

## A. Clay target

**(B) global regularity on $\mathbb T^3$, unforced** — and, contrapositively, **(D)**. Domain $\mathbb T^3$; no forcing (all statements below remain valid with $f\in L^1_tX^0$, adding $\int\|f\|_{X^0}$ to every ledger, but forcing is *not* used). Data: $u_0\in C^\infty(\mathbb T^3)$, $\nabla\cdot u_0=0$, mean zero, no symmetry assumed (deliberately **not** axisymmetric — see §E.11–E.12). Viscosity $\nu>0$ fixed, never sent to $0$; $\nu$ appears explicitly in every constant. Whole-space transfer is *not* claimed (HS-5 discipline).

## B. Central mathematics

### B.1 The new quantity

$$\boxed{\;M_{\rm eff}(N,t)=\frac{\big(\sum_{k\in S_N}|\hat u_k|\big)^2}{\sum_{k\in S_N}|\hat u_k|^2},\qquad \mu_N(t)=\sqrt{\frac{M_{\rm eff}(N,t)}{N^3}},\qquad \mu_*(t)=\sup_{N=2^j}\mu_N(t).\;}$$
Cauchy–Schwarz gives $1\le M_{\rm eff}\le M_N^{\rm act}\le|S_N|=\tfrac{28\pi}{3}N^3(1+o(1))$, hence $0<\mu_N\le 5.42$. $\mu_N$ is **amplitude-scaling invariant** ($u\mapsto\lambda u$) and **phase-blind**; it is a pure shape number. Exact identity (definition, no loss):
$$W_N=\sqrt{M_{\rm eff}(N)}\,\|P_Nu\|_2=\mu_N N^{3/2}\|P_Nu\|_2 .\tag{B1}$$

### B.2 Wiener norm evolution and its sign structure

Write $\hat u_k=a_ke_k$, $|e_k|=1$. Since $P_ke_k=e_k$ and $\mathrm{Re}(-iz)=\mathrm{Im}\,z$,
$$\boxed{\;\dot a_k=\sum_{\ell+m=k}a_\ell a_m\,\sigma(k;\ell,m)-\nu|k|^2a_k,\qquad \sigma(k;\ell,m):=\mathrm{Im}\big[(m\cdot e_\ell)(e_m\cdot\bar e_k)\big].\;}\tag{B2}$$
$|\sigma|\le|m|$. **Sign structure:** amplitudes enter *positively*, all sign information sits in the phase triple $\sigma$; the $L^2$ cancellation $\langle u,B(u,u)\rangle=0$ (repo F-12) reads
$$G_k:=\sum_{\ell+m=k}a_\ell a_m\sigma(k;\ell,m),\qquad \sum_k a_kG_k=0,\qquad \Pi_N:=\sum_{k\in S_N}a_kG_k=\text{shell flux}.\tag{B3}$$
The Wiener production is the **unweighted** sum $\mathcal P=\sum_kG_k$. With $\bar a_N=W_N/M^{\rm act}_N$ on the active support,
$$\sum_{k\in S_N}G_k=\underbrace{\frac{\Pi_N}{\bar a_N}}_{\text{flux gain}}+\underbrace{\sum_{k\in S_N}G_k\Big(1-\frac{a_k}{\bar a_N}\Big)}_{=:R_N},\qquad \sum_{k\in S_N}\Big(1-\frac{a_k}{\bar a_N}\Big)^2=M^{\rm act}\Big(\frac{M^{\rm act}}{M_{\rm eff}}-1\Big).\tag{B4}$$
Hence the **exact per-shell Wiener ledger**
$$\boxed{\;\frac{d W_N}{dt}=\frac{M^{\rm act}_N}{W_N}\Pi_N+R_N-\nu\!\!\sum_{k\in S_N}\!|k|^2a_k,\qquad |R_N|\le\|G\|_{\ell^2(S_N)}\,M^{\rm act}\sqrt{\tfrac1{M_{\rm eff}}-\tfrac1{M^{\rm act}}}\;}\tag{B5}$$
with $\|G\|_{\ell^2}\le\|u\|_{X^0}\|\nabla u\|_2$ (Young). **Consequence (design rule):** $R_N\equiv0$ *identically* when the shell is equi-amplitude ($M_{\rm eff}=M^{\rm act}$). Energy methods see only $\sum a_kG_k$; the deficit of $\sum G_k$ from that is measured *exactly* by $1-M_{\rm eff}/M^{\rm act}$. This is the precise sense in which $M_{\rm eff}$ is what energy methods miss. Dimensions: $[\Pi_N/\bar a]=T^{-3}/T^{-1}=T^{-2}=[\dot W]$ ✓.

### B.3 Exact structure identity: $X^{-1}$ = delocalization-weighted critical Besov

On $S_j$, $2^{-j-1}<|k|^{-1}\le2^{-j}$, and $2^{-j}W_j=\mu_jb_j$ with $b_j:=2^{j/2}\|P_ju\|_2$ (critical Besov $\dot B^{1/2}_{2,1}$ coefficients). Therefore
$$\tfrac12\sum_j\mu_jb_j\;\le\;\|u\|_{X^{-1}}\;\le\;\sum_j\mu_jb_j.\tag{B6}$$
Splitting the sum at $2^J$ and Cauchy–Schwarz on each piece ($\sum_{j\le J}2^j\le2^{J+1}$, $\sum_{j>J}2^{-j}=2^{-J}$, $\sum_j2^{2j}\|P_ju\|_2^2\le\|\nabla u\|_2^2$), then optimizing over integer $J$ (loss $\le(1+\sqrt2)/2$):
$$\boxed{\;\|u\|_{X^{-1}}\le 3\,\mu_*(t)\,\big(\|u\|_{L^2}\|\nabla u\|_{L^2}\big)^{1/2}.\;}\tag{B7}$$
Both factors on the right are **energy-method quantities**; $\mu_*$ is the only new input.

### B.4 The engine: $X^{-1}$ differential inequality with constant $1$

Since $\ell\cdot\hat u_\ell=0$, $m\cdot\hat u_\ell=k\cdot\hat u_\ell$, so $|k|^{-1}|\mathcal N_k|\le\sum_{\ell+m=k}a_\ell a_m$ ($|P_kv|\le|v|$). Summing,
$$\frac{d}{dt}\|u\|_{X^{-1}}+\nu\|u\|_{X^1}\le\|u\|_{X^0}^2\le\|u\|_{X^{-1}}\|u\|_{X^{1}},\tag{B8}$$
the last step by Cauchy–Schwarz $\sum a_k=\sum(|k|^{-1/2}a_k^{1/2})(|k|^{1/2}a_k^{1/2})$. Hence $\frac{d}{dt}\|u\|_{X^{-1}}\le(\|u\|_{X^{-1}}-\nu)\|u\|_{X^1}$. If $\|u(t_0)\|_{X^{-1}}<\nu$ then $\|u\|_{X^{-1}}$ is nonincreasing thereafter and $\nu\int_{t_0}^{T}\|u\|_{X^1}\le\|u(t_0)\|_{X^{-1}}<\infty$; since $\|\omega\|_\infty\le2\|u\|_{X^1}$, BKM/Kato–Ponce continues the solution past $T$. Therefore:

> **(C1)** If the maximal smooth solution breaks down at $T_*<\infty$, then $\|u(t)\|_{X^{-1}}\ge\nu$ for **every** $t<T_*$.

(This reproves Lei–Lin smallness inline; no external theorem is imported except BKM.)

### B.5 THEOREM D (Delocalization Dichotomy)

Combining (C1) and (B7), exactly one of the following holds for smooth mean-zero divergence-free $u_0$:

**(A)** the solution is global smooth; **(B)** $T_*<\infty$ and
$$\boxed{\;\mu_*(t)\ \ge\ \frac{\nu}{3\big(\|u(t)\|_2\|\nabla u(t)\|_2\big)^{1/2}}\quad\forall t<T_*,\qquad \sup_{t<T_*}\mu_*^2\ \ge\ \frac{\sqrt2}{9}\,\frac{\nu^{5/2}T_*^{1/2}}{\|u_0\|_{2}^2}=\frac{\sqrt2}{9}\frac{\tilde T^{1/2}}{\mathrm{Re}^2}\;}$$
with $\tilde T=\nu T_*$, $\mathrm{Re}=\|u_0\|_2/\nu$. Second bound: $T\sup\mu_*^2\ge\frac{\nu^2}{9}\int(\|u\|\|\nabla u\|)^{-1}\ge\frac{\nu^2T^2}{9\int\|u\|\|\nabla u\|}$ and $\int_0^T\|u\|_2\|\nabla u\|_2\le\|u_0\|_2^2\sqrt{T/2\nu}$ by Cauchy–Schwarz + the energy identity $\int\|\nabla u\|^2\le\|u_0\|_2^2/2\nu$. Dimensions: $[\nu^{5/2}T^{1/2}]=T^{-2}=[\|u_0\|^2]$ ✓.

*Non-vacuity check (and honesty check):* since $\mu_*\le5.42$, branch (B) is refutable only if $\|u\|_2\|\nabla u\|_2<\nu^2/265$; that never holds at blow-up, so Theorem D **does not accidentally prove regularity** — as it must not.

### B.6 Seed S4 worked out: flux $\Rightarrow$ per-octave $L^3$

$\Pi_N=\langle u\otimes u,\nabla\Delta_Nu\rangle$ (integrate by parts; $\Delta_N$ = **smooth** LP block, see §H-R2). Hölder $(3/2,3)$ + Bernstein $\|\nabla\Delta_Nf\|_3\le B_1N\|\Delta_Nf\|_3$ ($B_1=\|\nabla\check\psi\|_{L^1}$, $\psi\equiv1$ on $\mathrm{supp}\chi$; torus constants by de Leeuw transference) give
$$|\Pi_N|\le B_1N\|u\|_3^2\,\|\Delta_Nu\|_3.\tag{B9}$$
If the flux is $n_0$-**scale-local** (both $u$ factors within $n_0$ octaves of $N$, fraction $\ge\theta_0$ of $\Pi_N$), then with $g_j=\|\Delta_ju\|_3$ and $g_{j'}\le g_{\max}$ over the window, $\theta_0\Pi_N\le B_1N(2n_0{+}1)^2g_{\max}^3$. A sustained critical relay budget $\Pi_N\ge q\,c_EN$ then forces
$$\boxed{\;g_{\max}\ \ge\ \varepsilon_0:=\Big(\tfrac{\theta_0qc_E}{B_1(2n_0+1)^2}\Big)^{1/3}\;}\quad\text{(scale-independent).}\tag{B10}$$
Since $\ell^2\subset\ell^3$ in the LP square function for $p=3\ge2$, $\sum_j g_j^3\le c_3\|u\|_3^3$, so with $J$ active/wake octaves $\|u\|_3^3\ge c_3^{-1}J\varepsilon_0^3$: **ESS divergence follows from the flux budget alone**, log-divergent in $N$ since $J=\log_2N$. Finally, by $\|f\|_3\le\|f\|_\infty^{1/3}\|f\|_2^{2/3}$, $\|\Delta_Nu\|_\infty\le W_N$ and (B1):
$$\|\Delta_Nu\|_3\le\mu_N^{1/3}\,N^{1/2}\|\Delta_Nu\|_2=\mu_N^{1/3}b_N\ \Longrightarrow\ \boxed{\;\mu_N\ \ge\ \varepsilon_0^3/b_N^3=\varepsilon_0^3/(2c_E)^{3/2}\;}\tag{B11}$$
on the critical wake $b_N=\sqrt{2c_E}$ — a **scale-independent delocalization floor on every flux-active octave**, obtained from the flux side, entirely independently of the repo's response-side no-go §6.

### B.7 Frozen-parent $\to$ Duhamel-with-remainder

Mild form with $t-t_0=\tau N^{-2}$: $\int_{t_0}^te^{-\nu|k|^2(t-s)}ds=N^{-2}h_{\nu,\tau}(|k|/N)$ — the repo's heat factor is *exactly* the Duhamel kernel mass. Splitting $\mathcal N_k(s)=\mathcal N_k(t_0)+\rho_k(s)$,
$$\hat u_k(t)=e^{-\nu|k|^2\Delta t}\hat u_k(t_0)+N^{-2}h_{\nu,\tau}(|k|/N)\mathcal N_k(t_0)+r_k,\quad |r_k|\le N^{-2}h_{\nu,\tau}\Delta t\sup_s|\partial_s\mathcal N_k|,$$
and $|\partial_s\mathcal N_k|\le 4N(\|u\|_{X^0}\|u\|_{X^1}+\nu\|u\|_{X^2})_{\rm loc}$ gives relative remainder $\le\tau(\|u\|_{X^1}N^{-2}+\nu)=\tau(\sqrt{2c_E}\mu+\nu)$ on the critical wake — **$O(\tau)$ uniformly in $N$**. So the repo's frozen-parent no-go is not a model: it is the true evolution up to explicit relative $O(\tau)$.

## C. Scaling table ($\tau=T-t$; front $N(\tau)=(2a\tau)^{-1/2}$, $\gamma=1/2$ (seed S1); wake $E_N=c_E/N$; $\mu_N\asymp\mu$)

| quantity | formula | exponent in $\tau$ |
|---|---|---|
| energy $\|u\|_2^2$ | $\le4c_E$ | $\tau^{0}$ (bounded ✓) |
| enstrophy $\|\nabla u\|_2^2$ | $\asymp4c_EN$ | $\tau^{-1/2}$ |
| dissipation $\nu\|\nabla u\|_2^2$ | $4\nu c_E(2a\tau)^{-1/2}$ | $\tau^{-1/2}$, $\int<\infty$ ✓ |
| global $L^3$ | $\|u\|_3^3\ge c_3^{-1}\varepsilon_0^3J$ | $(\log\tau^{-1})^{1/3}$ → $\infty$ (ESS ✓) |
| $\|\omega\|_\infty$ | $\le2\|u\|_{X^1}\asymp\sqrt{2c_E}\mu N^2$ | $\tau^{-1}$, $\int=\infty$ (BKM ✓) |
| nonlinear $\|P(u\!\cdot\!\nabla u)\|_2$ | $\le\|u\|_\infty\|\nabla u\|_2\asymp2\sqrt{c_E}N^{3/2}$ | $\tau^{-3/4}$ |
| pressure $\|\nabla p\|_2$ | $\le$ same (Riesz, const 1) | $\tau^{-3/4}$ |
| time remaining | $\tau=N^{-2}/2a$, $dN/ds=aN$ | $e^{-2as}$, $\sum<\infty$ ✓ |
| bandwidth $N$ | $(2a\tau)^{-1/2}$ | $\tau^{-1/2}$ |
| active modes $M^{\rm act}\ge M_{\rm eff}$ | $\mu^2N^3$ | $\tau^{-3/2}$ |
| **$M_{\rm eff}/N^3=\mu^2$** | $\ge\varepsilon_0^6/(2c_E)^3$ | $\tau^{0}$ |
| $\|u\|_{X^0}$ / $\|u\|_{X^{-1}}$ | $\sqrt{2c_E}\mu N$ / $\sqrt{2c_E}\mu J$ | $\tau^{-1/2}$ / $\log\tau^{-1}$ ($\ge\nu$ ✓ C1) |
| $\sqrt\tau\|u\|_\infty$ | $\asymp\sqrt{2c_E}\mu/\sqrt{2a}$ | $\tau^0$ (Type-I boundary) |

## D. Closed feedback loop (every arrow a formula)

1. **Flux $\to$ Wiener**: (B5) with equi-amplitude ($R_N=0$): $\dot W_N=M^{\rm act}\Pi_N/W_N-\nu\sum|k|^2a_k$.
2. **Wiener $\to$ amplitude**: (B1) $W_N=\mu_N\sqrt{2c_E}\,N$ at critical normalization.
3. **Budget**: gain $=q c_EN\,M^{\rm act}/W_N=q\sqrt{c_E/2}\,M^{\rm act}/\mu_N$; loss $\le4\nu N^2W_N=4\nu\sqrt{2c_E}\mu_NN^3\!/N$. With $M^{\rm act}\ge M_{\rm eff}=\mu_N^2N^3$: gain/loss $\ge qN/(8\nu)$, so the loop is self-sustaining once
$$\boxed{\;N\ \ge\ N_\nu:=8\nu/q.\;}$$
4. **Wiener $\to L^3$**: (B11) reversed, $\|\Delta_Nu\|_3\le\mu_N^{1/3}\sqrt{2c_E}$, and $\|u\|_3^3\gtrsim J\varepsilon_0^3$.
5. **$L^3\to$ front advance**: seed S2 flux closure $F(1)=\chi(2c_E)^{3/2}$, $\xi_{\max}=1+\tfrac{\chi}{\sqrt2}\sqrt{c_E}/\nu\ge2$ hands the budget to octave $2N$; $\dot N=aN^3$, $a=\Ṅ/N^3$.
6. **Back to 1** at $N\!\to\!2N$ with $\mu_{2N}=\mu_N$ (shape recurrence). Closure condition: $\chi\sqrt{c_E}\ge\sqrt2\nu$ **and** $\mu\ge\varepsilon_0^3/(2c_E)^{3/2}$ **and** $M^{\rm act}=M_{\rm eff}$ (equi-amplitude, $R_N=0$).

## E. Obstruction audit (exact collision points)

1. **Energy bound / F-N1.** No collision: $\|u\|_2$ bounded is an *input* to (B7) and to the $\tilde T^{1/2}/\mathrm{Re}^2$ floor. $\mu_N$ is amplitude-invariant, so it is not a norm and cannot be bounded by the energy inequality.
2. **Finite dissipation / F-N2, N-2.** No collision: $\int\|\nabla u\|^2\le\|u_0\|^2/2\nu$ is used to *derive* the floor. The mechanism's dissipation $\asymp\tau^{-1/2}$ is integrable ✓.
3. **ESS $L^\infty_tL^3$.** Not contradicted — §B.6 *produces* $\|u\|_3\to\infty$ (log rate). Collision point: a candidate with bounded $L^3$ and $n_0$-local sustained flux contradicts (B10)+(B11); such candidates are newly excluded.
4. **Fixed-finite-bandwidth no-go (F-$\alpha$1, VR-L-011).** Independent, no overlap: my loop requires $N(t)\to\infty$ (arrow 5–6). Explicitly, for fixed band $\|u\|_{X^{-1}}\le\sqrt{M_{\max}}\|u_0\|_2$ is bounded but not $<\nu$, so Theorem D does *not* reprove F-$\alpha$1 — I do not claim it does.
5. **Pure-swirl $L^3$ no-go (VR-L-016).** $\mu_N$ is phase-blind, hence *blind* to the pressure-channel sign. Collision point: the mechanism assumes $\Pi_N>0$, and pure swirl has $P\equiv0$ at $t=0$, so the pure-swirl class fails the hypothesis of (B10) at $t=0$. Evasion = never use pure swirl; keep LG-11 parity ($\psi_1$ odd in $z$) as a design constraint on the flux input.
6. **One-scale self-similar (NRS/Tsai).** No collision *and an admitted limitation*: $\mu_N$ is scale-invariant, so an exact Leray profile has $\mu_N\equiv$ const and satisfies Theorem D. My theorem is orthogonal to self-similarity; NRS/Tsai kill that class on other grounds and I keep them as an independent filter (the front is $s$-periodic, not one-scale stationary).
7. **Galerkin global existence.** Same as 4: the ledger (B5) is finite-mode-exact, so it is consistent with $\dot W$ being globally bounded on any *fixed* band; unboundedness enters only through $N(t)\to\infty$ in arrow 5.
8. **Smooth-forcing high-frequency decay (F-N4).** Not used (unforced target). If forcing were added: $|\dot W_N|_f\le\sum_{S_N}|\hat f_k|\le C_mN^{3-m}\to0$ for $m>3$ — the ledger is untouched at high $N$, so the mechanism neither needs nor benefits from forcing.
9. **Mesoscopic $\gamma<1$ no-go $D_N\le2\kappa^2\tau^2c_EM_{\rm eff}/N^3=2\kappa^2\tau^2c_E\mu_N^2$.** **Direct collision, same quantity, same direction — this is the strongest support.** The repo derives $\mu_N\gtrsim1$ as *necessary for the relay ansatz* (response side); (B11) derives $\mu_N\ge\varepsilon_0^3/(2c_E)^{3/2}$ as *necessary for any sustained-local-flux blow-up* (flux side). Two independent derivations of the same floor. Consequence: my construction must live at $\gamma=1$, fixed relative width $W_N=\lfloor\eta N\rfloor$, $\eta\in(0,1/3)$ — precisely the repo's only surviving regime. All $\gamma<1$ families are pre-rejected here too.
10. **Diagonal cross-talk gate.** Collision point: cross-talk children arrive with amplitudes unequal to the intended children (measured $\|B_{\rm cross}\|^2=2483/1890$ vs $\|B_{\rm int}\|^2=37/315$), which destroys equi-amplitude and hence turns on $R_N\ne0$ in (B5) with $|R_N|\le\|G\|_2M^{\rm act}\sqrt{1/M_{\rm eff}-1/M^{\rm act}}$. Evasion: seed S6(b) **scale-stagger** ($\rho=3/4$) makes the diagonal sumset geometrically miss every tagged box (exactly zero, integer-programming checkable), which simultaneously restores $M_{\rm eff}=M^{\rm act}$. S6(a) translation-split is rejected in the seed and is also rejected here (it decoheres the intended stage-2 pair by the same $(W|\Delta x|)^{-3}$).
11. **CSTY Type-I exclusion.** Applies to axisymmetric solutions on $\mathbb R^3$ only. My target is (B) on $\mathbb T^3$ with no symmetry, so it does not apply — stated as a *domain/symmetry-class* evasion, not a bound violation. Note honestly that the $\gamma=1/2$ front sits at $\sqrt\tau\|u\|_\infty\asymp$ const, i.e. exactly the Type-I boundary; any migration of this mechanism to axisymmetric $\mathbb R^3$ would collide with CSTY (2.1) head-on and is forbidden without a Type-II modification $\gamma>1/2$.
12. **KNSS ancient Liouville.** $\mu_*$ is scale-invariant, so it *passes to the rescaled ancient limit*: the limit $U$ inherits $\mu_N(U)\ge\varepsilon_0^3/(2c_E)^{3/2}$ — a new constraint on ancient solutions. KNSS kills axisym-no-swirl and $|u|\le C/r$; my limit is non-axisymmetric with $\|u\|_\infty\asymp N$ (no $C/r$ decay), so neither branch applies.
13. **Front-resolution threat model.** $\mu_N$ is a spectral diagnostic; specific threats: TM-03 aliasing inflates high-$k$ $a_k$ and hence *spuriously inflates* $M_{\rm eff}$ (dealias 2/3-rule + report $\mu$ at two cutoffs); TM-04 spectral blocking piles energy in the last shells (drop the top two octaves from $\mu_*$); TM-22 requires the front shell fully resolved, i.e. grid cutoff $K\ge4N$ so the convolution support $2S_N$ is exact; TM-21 requires $\mu_N$ logged every step, not strided (it is $O(1)$-varying).

## F. Falsification pilot ($\le$1 h)

**Variables.** $\mu_N(t),M_{\rm eff},M^{\rm act},\Pi_N,R_N,W_N$, threshold $\Theta(t)=\nu/(3(\|u\|_2\|\nabla u\|_2)^{1/2})$, ledger residual $\mathcal R=\dot W_N-M^{\rm act}\Pi_N/W_N-R_N+\nu\sum|k|^2a_k$.
**Exact vs float.** *Exact rational/interval required*: (B3) $\sum_ka_kG_k=\Pi_N$, deficit identity (B4), $R_N=0$ for equi-amplitude, and $\mu_N$ of the exact carrier alphabets. *Float permitted*: the Galerkin time series, $\Theta(t)$, and the $\mu_N$ trajectory.
**Reuse.** `mesoscopic_local_fft.local_fft_leray_coefficients` (exact zero-padded convolution → $G_k$ per mode, $K=4W{-}3$ padding already correct for TM-03); `mesoscopic_cloud_scaling.exact_sparse_leray_convolution` and `measure_mesoscopic_cloud` (parents at $\eta=0.10,0.15,0.20$, $N=16,32,48,64$); `mesoscopic_galerkin.run_small_mesoscopic_galerkin` (float time series for $\mu_*(t)$ vs $\Theta(t)$); `leray_response_relay.leray_project/leray_advection`; `exact_carrier_record_verifier` as the independent second code path for $G_k$ (mandatory, TM-14).
**Resolution.** $N\in\{16,32,48,64\}$ sparse; $64^3$ grid, cutoff $3N$, RK4$\times$16 for the Galerkin run — matches existing runtimes.
**Success.** (S1) $|\mathcal R|=0$ exactly in rationals on the $N=4$ carrier gadget; (S2) $R_N=0$ exactly for an equi-amplitude parent and $\ne0$ for the known cross-talk-contaminated parent; (S3) $\mu_*(t)\ge\Theta(t)$ at every logged time (necessary-condition consistency; violation ⟹ bug); (S4) $\mu_N$ for the $\eta$-family reported against the floor $\varepsilon_0^3/(2c_E)^{3/2}$.
**Kill.** If (S1) or (S2) fails, the ledger (B5) is wrong ⟹ kill the whole lens. If $\mu_N\to0$ as $N$ grows across $\{16,32,48,64\}$ for the surviving $\eta\in(0,1/3)$ family, then the *only* structurally uncondemned mesoscopic family also fails the delocalization floor ⟹ the entire cloud lane dies by a new, independent argument (a valuable negative result).

## G. Proof chain to a Clay statement

1. **(L)** Finite-Fourier Wiener calculus: (B2), $|P_kv|\le|v|$, $|\sigma|\le|m|$.
2. **(L)** Exact ledger (B3)–(B5) incl. deficit identity $\sum(1-a_k/\bar a)^2=M(M/M_{\rm eff}-1)$.
3. **(L)** $X^{-1}$ bilinear estimate with constant 1 and interpolation (B8).
4. **(M)** (C1): breakdown $\Rightarrow\|u\|_{X^{-1}}\ge\nu$ (needs BKM/Kato–Ponce as the one classical import).
5. **(L/M)** Split inequality (B7) and Theorem D.
6. **(M)** S4-local per-octave $L^3$ lemma (B9)–(B11); requires explicit $A_3,B_1$ via de Leeuw transference — the only non-elementary analytic import.
7. **(I)** Construct an equi-amplitude, scale-staggered ($\rho=3/4$), fixed-relative-width ($\eta<1/3$) carrier with $\Pi_N\ge qc_EN$, $R_N=0$, cross-talk geometrically empty; interval-verify $q>0$.
8. **(I)** One-octave interval stage budget closing arrow 3 ($N\ge8\nu/q$) and arrow 5 ($\chi\sqrt{c_E}\ge\sqrt2\nu$).
9. **(I)** Shape recurrence: forward RG integration of the S1 front flow to an attracting $s$-periodic orbit with $\mu$ constant; Floquet-verify.
10. **(M)** Assemble: recurrence $\Rightarrow\int_0^{T}\|\omega\|_\infty=\infty$ and $T<\infty$ (§C row 8) $\Rightarrow$ Clay (D)/(B) breakdown. *Negative branch:* if steps 7–9 show the window $\varepsilon_0^3/(2c_E)^{3/2}\le\mu_N\le(1-\theta)^{1/2}$ is empty, that is a conditional regularity theorem instead.

## H. Sub-variants REJECTED during derivation

- **R1 (seed S4 unconditional).** "$\|u\|_3^3\gtrsim\log N$ with no locality hypothesis" — **REJECTED** at $\varepsilon=qc_E/(2B_1\|u\|_3^2)$: the Hölder factor $\|u\|_3^2$ in (B9) is not $O(1)$, and self-consistency $J\varepsilon^3\le c_3\|u\|_3^3$ yields only $\|u\|_3\ge\kappa J^{1/9}$. Repaired only by the $n_0$-locality hypothesis (B10).
- **R2 (sharp shell cutoffs in $L^3$).** **REJECTED**: the ball multiplier is unbounded on $L^p(p\ne2)$ (Fefferman), so $\|P_Nu\|_3\le C\|u\|_3$ fails for sharp annuli. Fix: smooth LP blocks in §B.6; all $\ell^1$/$\ell^2$ statements (B1)–(B8) keep sharp shells, where sharp cutoffs are contractions.
- **R3 (first split).** $\|u\|_{X^{-1}}\le\sqrt2\mu_*\|\nabla u\|_2$ — **SUPERSEDED** (not false) by (B7), which is stronger whenever $\|\nabla u\|_2>\|u\|_2$, i.e. always near breakdown.
- **R4 (Re-free floor).** Attempt to remove $\mathrm{Re}^{-2}$ — **REJECTED**: the step $\int_0^T\|u\|_2\|\nabla u\|_2\,dt\le\|u_0\|_2^2\sqrt{T/2\nu}$ is saturated by Type-I profiles, so no $\mathrm{Re}$-free floor is extractable from the energy identity alone.
- **R5 (per-shell version of Theorem D).** "$M_{\rm eff}(N,t)/N^3\ge\delta$ for **all** $N$" — **REJECTED**: false, a single-mode low shell has $M_{\rm eff}=1$, $\mu_N=N^{-3/2}$. Only $\sup_N$ survives unconditionally; the per-shell floor requires the flux hypothesis (B10).
