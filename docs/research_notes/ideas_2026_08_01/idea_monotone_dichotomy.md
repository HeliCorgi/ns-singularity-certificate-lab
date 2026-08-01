# LENS 11 — The Closable Front Monotone $\Lambda$ and the Bandwidth–Dissipation Dichotomy

**Status: PROOF CANDIDATE** (regularity side).
Core objects (Identity I, Lemma K, Monotone $\Lambda$, Dichotomy D) are *proved* at the stated
hypothesis level and are exactly computable in rational arithmetic by existing repo machinery.
The Clay-(B) statement remains conditional on obligations O-1…O-6 (§G). Three sub-variants died
during derivation and are kept as **REJECTED** with their exact failing equations (§B.5).

---

## A. Clay target

* **Target: (B)** — global regularity, $\mathbb T^3$, no forcing. (Not (A): the reason is
  structural, not cosmetic — the auxiliary identity $\int_0^T N_0^2\,dt<\infty$ (§B.3) requires
  $H_0(T)>0$, which fails for whole-space self-similar profiles whose energy drains to zero. On
  $\mathbb T^3$ with fixed data this is an obligation (O-3), not a false step.)
* **Domain** $\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3$, convention $u=\sum_k\hat u_ke^{ik\cdot x}$,
  $k\in\mathbb Z^3\setminus\{0\}$ (zero mean), $k\cdot\hat u_k=0$, $\hat u_{-k}=\overline{\hat u_k}$.
* **Forcing:** none. (Forced variant in §E-8: it adds $\sum_kx_k^r\langle\hat u_k,\hat f_k\rangle$ to
  $T_r$; admissible $f\in L^1_tL^2_x$, but not used.)
* **Initial regularity:** $u_0\in C^\infty(\mathbb T^3)$, $\nabla\!\cdot u_0=0$, zero mean; maximal
  strong solution $u\in C([0,T_{\max});H^m)$, $m>5/2$.
* **Viscosity:** $\nu>0$ fixed, kept explicitly in every constant. Nothing here survives $\nu\to0$:
  the monotone's whole content is a viscous square completion.

Repo conventions used verbatim: $\partial_t\hat u_k=\mathcal N_k-\nu|k|^2\hat u_k$,
$\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m$, $P_k=I-k\otimes k/|k|^2$,
critical shell normalization $E_N=c_E/N$.

---

## B. Central mathematics

### B.1 Modal ledger (repo convention, `modal_front_actions.py`)

Paired real coefficients give, for each $k\neq0$, $e_k=\tfrac12(|c_k|^2+|s_k|^2)\ge0$,
$a_k=\tfrac12(c_k\!\cdot\!c^{\mathcal N}_k+s_k\!\cdot\!s^{\mathcal N}_k)$,
$n_k=\tfrac12(|c^{\mathcal N}_k|^2+|s^{\mathcal N}_k|^2)$, $x_k=|k|^2$, and
$$\dot e_k=2a_k-2\nu x_ke_k .$$
Moments $H_r=\sum x_k^re_k$, $T_r=\sum x_k^ra_k$, $G_r=\sum x_k^rn_k$, so
$\tfrac12\dot H_r=T_r-\nu H_{r+1}$, $H_0=\tfrac12\|u\|_2^2$, $H_1=\tfrac12\|\nabla u\|_2^2$,
$G_0=\tfrac12\|\mathbb P(u\cdot\nabla u)\|_2^2$. Radial bandwidth $N_r^2=H_{r+1}/H_r$.

**Two exact facts used repeatedly.**
(i) *Energy neutrality* $T_0=\sum_ka_k=0$ (repo F-12) $\Rightarrow \dot H_0=-2\nu H_1$.
(ii) *Pressure blindness*: the pressure enters only as $-ik\hat p_k$ and
$\langle\hat u_k,\,ik\hat p_k\rangle=i\hat p_k(k\cdot\overline{\hat u_k})=0$. Hence **every quantity
below is pressure-free** — the entire $\Lambda$-lane is immune to the pressure-channel pathologies
(LG-4/LG-9, VR-C-009) that killed the $L^3$-generation lane.

### B.2 Identity I (exact front identity with a *closable* right-hand side)

Let $p_r(k)=x_k^re_k/H_r$, $\mu=N_r^2=\mathbb E_{p_r}[x]$, $g_k=a_k/e_k$,
$V_r=\mathrm{Var}_{p_r}(x)=H_{r+2}/H_r-\mu^2$. From $\tfrac12\dot H_r=T_r-\nu H_{r+1}$,
$$\tfrac12\tfrac{d}{dt}\log N_r^2=\underbrace{\Big(\tfrac{T_{r+1}}{H_{r+1}}-\tfrac{T_r}{H_r}\Big)}_{=\,\mathrm{Cov}_{p_r}(x,g)/\mu}-\ \nu\frac{V_r}{\mu}. \tag{I.1}$$
(The covariance form uses $\sum_kx_k^r(x_k-\mu)e_k=H_{r+1}-\mu H_r=0$, which also **removes** the
mean-growth term: $\mathrm{Cov}_{p_r}(x,g)=H_r^{-1}\sum_kx_k^r(x_k-\mu)a_k$ exactly. This is the
identity already asserted and rationally verified by `modal_growth_identity`.)

Now the new step. The repo bounds $\mathrm{Cov}$ by $\sqrt{V_r\,\mathrm{Var}_{p_r}(g)}$, which is
*sharper* but **not closable**: $\mathrm{Var}_{p_r}(g)=H_r^{-1}\sum x^r(a_k-\bar g e_k)^2/e_k$ is
dominated by spectral near-nodes $e_k\to0$ and is bounded by no Sobolev norm. Instead use the
modal Cauchy–Schwarz $|a_k|\le\sqrt{e_kn_k}$ **before** dividing:
$$|\mathrm{Cov}_{p_r}(x,g)|\le\frac1{H_r}\sum_kx_k^r|x_k-\mu|\sqrt{e_kn_k}
\le\frac1{H_r}\Big(\underbrace{\sum x^r(x-\mu)^2e_k}_{=H_rV_r}\Big)^{1/2}\Big(\underbrace{\sum x^rn_k}_{=G_r}\Big)^{1/2}
=\sqrt{\frac{V_rG_r}{H_r}}. \tag{I.2}$$
Insert into (I.1) and complete the square in $\sqrt{V_r}\ge0$:
$$\tfrac12\tfrac{d}{dt}\log N_r^2\le\frac1\mu\Big[\sqrt{G_r/H_r}\,\sqrt{V_r}-\nu V_r\Big]
\le\frac1{\mu}\cdot\frac{G_r/H_r}{4\nu}=\frac{G_r}{4\nu H_{r+1}} .$$
$$\boxed{\ \frac{d}{dt}\log N_r^2\ \le\ \frac{G_r}{2\nu H_{r+1}}\ }\tag{I.3}$$
**Exact gap decomposition** (both terms $\ge0$, both exactly rational):
$$\frac{G_r}{2\nu H_{r+1}}-\frac{d}{dt}\log N_r^2=\underbrace{\frac{2}{\mu}\Big[\sqrt{V_rG_r/H_r}-\mathrm{Cov}_{p_r}(x,g)\Big]}_{\Gamma^{\rm CS}_r}
+\underbrace{\frac{2\nu}{\mu}\Big[\sqrt{V_r}-\frac{1}{2\nu}\sqrt{G_r/H_r}\Big]^2}_{\Gamma^{\rm SC}_r}. \tag{I.4}$$
**Equality (rigidity) cases.** $\Gamma^{\rm CS}_r=0$ iff $n_k=c\,(x_k-\mu)^2e_k$ for one constant
$c\ge0$ **and** $\langle\hat u_k,\mathcal N_k\rangle=\mathrm{sgn}(x_k-\mu)|\hat u_k||\mathcal N_k|$
for all $k$ (perfect modal alignment with the sign flipping exactly at the bandwidth $\mu$).
$\Gamma^{\rm SC}_r=0$ iff $\ \mathrm{Var}_{p_r}(|k|^2)=G_r/(4\nu^2H_r)$ — a single scalar tuning of
the relative spectral variance. Saturating (I.3) requires **both simultaneously**.

*Sign check:* pure heat flow has $a_k\equiv0\Rightarrow G_r=0$ and (I.1) gives
$\frac{d}{dt}\log N_r^2=-2\nu V_r/\mu\le0$: viscosity always narrows the bandwidth. ✔
*Dimension check:* $[G_r/H_{r+1}]=U^4L^{-2-2r}L^3/(U^2L^{-2-2r}L^3)=U^2$, so
$[G_r/(\nu H_{r+1})]=U^2/(L^2T^{-1})=T^{-1}$. ✔
*Scale check ($r=0$):* under $u_\lambda=\lambda u(\lambda x,\lambda^2t)$, $G_0\mapsto\lambda^3G_0$,
$H_1\mapsto\lambda H_1$, so RHS $\mapsto\lambda^2\,$RHS, and $\frac{d}{dt}\log N_0^2\mapsto\lambda^2(\cdot)$.
**(I.3) at $r=0$ is exactly scale-invariant.** ✔

### B.3 The nonlinear front wavenumber and the monotone

Define, for $D(t)=\|\nabla u\|_2^2=2H_1$,
$$K(t):=\frac{\|\mathbb P(u\cdot\nabla u)\|_{L^2}^2}{\|\nabla u\|_{L^2}^4}=\frac{G_0}{2H_1^2},\qquad [K]=L^{-1},\quad K\mapsto\lambda K .$$
$K$ is a genuine **wavenumber** (scale-covariant, rotation/translation invariant, pressure-free).
Then $G_0/(2\nu H_1)=K D/(2\nu)$ and (I.3) becomes $\frac{d}{dt}\log N_0^2\le\frac{1}{2\nu}KD$, i.e.

$$\boxed{\ \Lambda(t):=\log N_0^2(t)-\frac{1}{2\nu}\int_0^tK(s)D(s)\,ds\ \text{ is non-increasing.}\ }\tag{M}$$

Two Hölder dominations show $\int KD$ is a **strictly smaller** critical action than the classical ones:
$$KD=\frac{\|\mathbb P(u\cdot\nabla u)\|_2^2}{\|\nabla u\|_2^2}\le\|u\|_{L^\infty}^2,\qquad
KD\le\|u\|_{L^6}^2\|\nabla u\|_{L^3}^2/\|\nabla u\|_2^2\le C_S^2\|\nabla u\|_{L^3}^2 . \tag{B.3}$$
So $\int_0^T\!KD\,dt$ is dominated by both the Serrin $(p,q)=(\infty,2)$ action and the critical
vorticity $(p,q)=(3,2)$ action. Its finiteness is therefore a **weaker hypothesis**, and the
resulting criterion **stronger**, than either.

**Auxiliary exact identity (free finite action).** $T_0=0\Rightarrow\frac{d}{dt}\log H_0=-2\nu N_0^2$, hence
$$\int_0^TN_0^2\,dt=\frac1{2\nu}\log\frac{H_0(0)}{H_0(T)}<\infty\quad\text{whenever }H_0(T)>0. \tag{B.4}$$

### B.4 Lemma K (unconditional bandwidth bound on $K$)

For zero-mean $u$ spectrally supported in $0<|k|\le N$,
$\|u\|_\infty\le\sum|\hat u_k|=\sum|k|^{-1}(|k||\hat u_k|)\le(\sum_{0<|k|\le N}|k|^{-2})^{1/2}\|\nabla u\|_2$,
and the lattice sum is $\le 4\pi N+C$. With $\|\mathbb P(u\cdot\nabla u)\|_2\le\|u\|_\infty\|\nabla u\|_2$:
$$\boxed{\ K(t)\le\frac{\|u\|_\infty^2}{\|\nabla u\|_2^2}\le 4\pi N(t)+C\ }\tag{K}$$
— elementary, finite-lattice, Lean-formalizable, **no phase or coherence assumption**.

### B.5 REJECTED sub-variants (kept, with exact failing equations)

* **(a) Shell entropy $S=-\sum p_j\log p_j$, $p_j=E_j/E$.** Exact production identity (derived here):
  $$\dot S=-\frac{2}{H_0}\sum_ka_k\log p_k+2\nu\,\mathrm{Cov}_p(x,\log p).$$
  Both terms are sign-definite only under extra hypotheses (forward cascade; monotone spectrum) and
  they have **opposite** signs, so $S$ is not monotone. **Decisive kill:** under the repo's least-excluded
  log-critical shell law $E_j=\lambda_j^{-1}$, $p_j\simeq2^{-j-1}$ and
  $S=\log2\sum_j(j+1)2^{-j-1}=2\log2$ — *bounded*. Shell entropy is blind to the cascade. REJECTED as a
  monotone; retained as a diagnostic. **Salvage (kept, useful):** with $p_k=|\hat u_k|^2/\sum|\hat u|^2$,
  the repo's central no-go quantity is exactly a Rényi perplexity,
  $M^{\rm eff}=(\sum|\hat u_k|)^2/\|u\|_2^2=e^{H_{1/2}(p)}\ \ge\ e^{S}$, so
  $M^{\rm eff}\gtrsim N^3\iff H_{1/2}\ge3\log N+O(1)$: **the mesoscopic $\gamma<1$ no-go is a Rényi-$\tfrac12$
  entropy floor**, and (since $H_{1/2}\ge S$) a Shannon-entropy floor would be sufficient but is *not* necessary.
* **(b) Log-convexity ratio $\Phi=H_1/\sqrt{H_0H_2}\le1$.**
  $\frac{d}{dt}\log\Phi=2[\tfrac{T_1}{H_1}-\tfrac12(\tfrac{T_0}{H_0}+\tfrac{T_2}{H_2})]-2\nu[N_1^2-\tfrac12(N_0^2+N_2^2)]$.
  The viscous bracket is the *second difference of $N_r^2$*, which log-convexity of $r\mapsto H_r$ does
  **not** sign. Exact rational counterexample (all modes realizable: $|k|^2=1$ via $(1,0,0)$, $=2$ via
  $(1,1,0)$, $=100$ via $(6,8,0)$), masses $e=(1,1,\tfrac1{1000})$:
  $H_0=\tfrac{2001}{1000},H_1=\tfrac{31}{10},H_2=15,H_3=1009$, so
  $N_0^2=\tfrac{3100}{2001},\ N_1^2=\tfrac{150}{31},\ N_2^2=\tfrac{1009}{15}$ and
  $N_1^2-\tfrac12(N_0^2+N_2^2)=\tfrac{150}{31}-\tfrac12\big(\tfrac{3100}{2001}+\tfrac{1009}{15}\big)=-29.5692\ldots<0$:
  **viscosity alone increases $\Phi$**. Sign-indefinite ⇒ REJECTED.
* **(c) $H_3=H_0N_0^2N_1^2N_2^2$ as a monotone.** The identity is a telescope
  $H_3/H_0=\prod_{r=0}^2H_{r+1}/H_r$, i.e. tautological; it carries no dynamics by itself. REJECTED as a
  monotone. **Retained corollary** (summing (I.3) over $r=0,1,2$):
  $\log(H_3/H_0)-\frac1{2\nu}\int\sum_{r=0}^2 G_r/H_{r+1}$ is non-increasing — a legitimate higher-order
  front control, but its $r\ge1$ actions are *not* scale-critical, so it cannot feed a Clay statement.
* **(d) The repo's own $\mathrm{Var}(g)/(4\nu\mu)$ bound.** Sharper than (I.3) pointwise
  ($\mathrm{Var}_{p_r}(g)\le\mathbb E_{p_r}[g^2]\le G_r/H_r$), but REJECTED as a *closable* bound: at any
  spectral node $e_k\to0$ with $a_k\not\to0$, $g_k=a_k/e_k\to\infty$ and no Sobolev norm controls it.

---

## C. Scaling table ($\tau=T-t$; log-critical law $E_j=\lambda_j^{-1}$, $\lambda_J=N=\tau^{-\gamma}$)

| quantity | expression | exponent in $\tau$ |
|---|---|---|
| energy $H_0$ | $\sum_j\lambda_j^{-1}\le2c_E$ | $\tau^{0}$ (bounded) |
| enstrophy $H_1$ | $\asymp c_EN$ | $\tau^{-\gamma}$ |
| global $L^3$ | $\|u\|_3^3\asymp J+1$ | $\gamma\log(1/\tau)$ (log-divergent) |
| $\|\omega\|_\infty$ | $\lesssim\sum\lambda_j^{5/2}E_j^{1/2}\asymp N^2$ | $\tau^{-2\gamma}$ |
| dissipation rate $2\nu H_1$ | $\asymp\nu c_EN$ | $\nu\tau^{-\gamma}$ |
| nonlinear term $\|\mathcal N\|_2$ | $=(K)^{1/2}\|\nabla u\|_2^2\lesssim N^{1/2}\!\cdot\!N$ | $\tau^{-3\gamma/2}$ |
| pressure term $\|\nabla p\|_2$ | $=\|(I-\mathbb P)(u\!\cdot\!\nabla u)\|_2\lesssim\|u\cdot\nabla u\|_2$ | $\tau^{-3\gamma/2}$; **contributes exactly $0$ to $a_k$** |
| physical time remaining | $\tau=N^{-1/\gamma}$ | — |
| Fourier bandwidth $N$ | front top | $\tau^{-\gamma}$ |
| energy-weighted bandwidth $N_0$ | $(H_1/H_0)^{1/2}\asymp(N/2)^{1/2}$ | $\tau^{-\gamma/2}$ |
| active mode count $M$ | ball-filling $\asymp N^3$ | $\tau^{-3\gamma}$ |
| front wavenumber $K$ | **required** $\gtrsim\tau^{\gamma-1}$; **allowed** $\le4\pi N=4\pi\tau^{-\gamma}$ | $K/N\asymp\tau^{2\gamma-1}$ |

**Window derivation (new, from $\Lambda$ + Lemma K only).** $\int_0^TD\,dt<\infty$ (finite dissipation)
forces $\gamma<1$. Blow-up forces $\int KD=\infty$ (§D); with $K\le4\pi N$ this needs $\int ND\,dt=\infty$,
i.e. $\int\tau^{-2\gamma}d\tau=\infty$, i.e. $\boxed{\gamma\ge\tfrac12}$. The repo's window
$\gamma\in[\tfrac12,1)$ is thus **re-derived without BKM, without phase assumptions, without a flux sign**.

---

## D. Closed feedback loop (every arrow a formula)

$$
\underbrace{N_0^2=\tfrac{H_1}{H_0}}_{\text{bandwidth}}
\xrightarrow[\ (I.1)\ ]{\ \tfrac12\frac{d}{dt}\log N_0^2=\frac{\mathrm{Cov}_{p_0}(x,g)}{\mu}-\nu\frac{V_0}{\mu}\ }
\underbrace{V_0=\mathrm{Var}_{p_0}(|k|^2)}_{\text{spectral width}}
\xrightarrow[\ (I.2)\ ]{\ \mathrm{Cov}\le\sqrt{V_0G_0/H_0}\ }
\underbrace{G_0=\tfrac12\|\mathcal N\|_2^2}_{\text{nonlinear power}}
$$
$$
G_0\xrightarrow[\ K:=G_0/2H_1^2\ ]{}K\xrightarrow[\ (K)\ ]{K\le4\pi N}N\xrightarrow[\ \text{shell law}\ ]{H_1\asymp c_EN}H_1\xrightarrow[\ \dot H_0=-2\nu H_1\ ]{}H_0\ \curvearrowright\ N_0^2 .
$$
**The loop is a negative feedback.** Growing $N_0$ requires (I.1) positive, which requires
$\mathrm{Cov}_{p_0}(x,g)>\nu V_0$; but the same $G_0$ that supplies $\mathrm{Cov}$ is capped by Lemma K at
$K\le4\pi N$, while the viscous leg drains $H_0$ at rate $2\nu H_1$ with $\int H_1<\infty$. The only
escape is the **double-saturation corridor**: $\Gamma^{\rm CS}_0=\Gamma^{\rm SC}_0=0$ in (I.4) *and*
Bernstein saturation $\|u\|_\infty^2\asymp N\|\nabla u\|_2^2$ (full phase coherence) — and at
$\gamma=\tfrac12$ *both* must hold with no slack. For $\gamma\in(\tfrac12,1)$ the corridor has slack
$K/N\asymp\tau^{2\gamma-1}\to0$: the mechanism must run at *vanishing* front efficiency, which is the
regime the monotone is designed to test.

---

## E. Obstruction audit (exact collision points)

1. **Energy bound (F-N1)** — no collision: energy divergence is never claimed. $\dot H_0=-2\nu H_1$ is
   *used as input* to produce (B.4). Correct use of a "dead diagnostic."
2. **Finite dissipation (F-N2)** — collision point is $\int_0^TD\,dt<\infty$; this is exactly what forces
   $K\gtrsim\tau^{\gamma-1}\to\infty$ in §C. Again input, not signature.
3. **ESS $L^\infty L^3$** — not used. The regularity leg of §D closes at
   $u\in L^\infty_tH^1\subset L^\infty_tL^6$, Serrin pair $(6,\infty)$: $2/\infty+3/6=\tfrac12\le1$. No
   endpoint theorem needed, so no circularity with ESS.
4. **Fixed-finite-bandwidth no-go (F-$\alpha$1 / VR-L-011)** — collision point: on a fixed mode set $S$,
   $N_0^2\le\max_{k\in S}|k|^2$, so $\Lambda$ is trivially bounded and (M) says nothing. **Consequence for
   the pilot: fixed-band test fields are non-informative**; the pilot must use a band-*doubling* family.
   $\Lambda$ is a diagnostic on trajectories, not an ansatz, so the no-go does not apply to it.
5. **Pure-swirl $L^3$ no-go (VR-L-016/LG-9)** — no collision: that no-go kills the pressure channel $P$ of
   the $L^3$ identity. Here the pressure contributes **identically zero** to every $a_k$ (§B.1(ii)), so the
   monotone neither gains nor loses from swirl symmetry, and the "critical Reynolds $=+\infty$" pathology
   of Gate-7 cannot recur.
6. **One-scale self-similar no-go (NRS/Tsai)** — collision point is (B.4). For the exact backward
   self-similar profile $N_0^2\asymp\tau^{-1}$ and $\int N_0^2dt=\infty$; (B.4) is consistent only because
   $H_0(T)=0$ there ($\|u\|_2^2\asymp\tau^{1/2}$). This is why the target is (B) on $\mathbb T^3$ and why
   O-3 ($H_0(T)>0$) is listed rather than assumed. On $\mathbb T^3$ with fixed data, $H_0(T)=0$ would force
   $u\to0$ in $L^2$, incompatible with blow-up.
7. **Galerkin global existence** — same collision as 4; $\Lambda$ makes no finite-mode claim.
8. **Smooth-forcing high-frequency decay (F-N4)** — not invoked (unforced). If forcing were added, $T_r$
   gains $\sum x_k^r\langle\hat u_k,\hat f_k\rangle$ and (I.2) gains
   $\sqrt{V_r\Phi_r/H_r}$ with $\Phi_r=\sum x^r|\hat f_k|^2$; F-N4's super-polynomial decay makes
   $\Phi_r/G_r\to0$, so (I.3) is asymptotically forcing-blind — consistent, not exploited.
9. **Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}/N^3$** — collision point
   is the identification $M^{\rm eff}=e^{H_{1/2}(p)}$ (§B.5a). $\Lambda$ proposes no relay, so the no-go does
   not bind it; conversely the no-go's requirement is now readable as a Rényi entropy floor, which the
   pilot logs. **Non-claim:** $S\le H_{1/2}$ means a Shannon-entropy measurement can never *establish*
   $M^{\rm eff}\gtrsim N^3$.
10. **Diagonal cross-talk gate** — no collision (no carrier alphabet). The gate's *fields* are reused as
    test data in §F; their measured $\Gamma^{\rm CS}_0$ is the quantitative version of "leakage-dominated."
11. **CSTY Type-I exclusion** — collision point: §C gives saturation exactly at $\gamma=\tfrac12$, the
    Type-I boundary $\sqrt\tau\|u\|_\infty\asymp1$. CSTY excludes axisymmetric Type-I; the monotone
    independently shows $\gamma=\tfrac12$ demands *zero slack* in (I.4) **and** in Bernstein. The two
    constraints pincer the corner from different directions; the surviving open window is
    $\gamma\in(\tfrac12,1)$.
12. **KNSS ancient-solution Liouville** — no collision: $\Lambda$ is not used to build a rescaled limit.
    Noted obligation only: $N_0^2$, $K$, $D$ all rescale covariantly ($\lambda^2,\lambda,\lambda$), so
    $\int KD\,dt$ is invariant and can in principle be transported to an ancient limit — not attempted here.
13. **Front-resolution threat model (TM-22, TM-03, TM-09)** — collision points, all mitigated:
    (i) $G_0,H_2$ are the highest moments and the first to be corrupted by aliasing ⇒ symbolic tests run in
    exact `Fraction` arithmetic on finite trig fields, where $\mathbb P(u\cdot\nabla u)$ is exact;
    (ii) $V_r=H_{r+2}/H_r-\mu^2$ is a catastrophic-cancellation trap ⇒ compute as
    $\sum p_r(x-\mu)^2$, never as the difference;
    (iii) evolution runs require the $\ge7$-points-per-front rule before any exponent for $K$ is fitted.

---

## F. Minimal falsification experiment ($\le1$ h)

**New module** `src/ns_certificate_lab/spectral_front_monotone.py`, reusing
`fourier_torus.TrigVector/advection/leray` and mirroring `modal_front_actions._moments`.

* `front_gap_identity(field, order, viscosity)` → exact `Fraction` record
  $\{H_r,H_{r+1},H_{r+2},T_r,T_{r+1},G_r,V_r,\mathrm{Cov},\Gamma^{\rm CS}_r,\Gamma^{\rm SC}_r,
  \text{closable\_upper}=G_r/(4\nu H_{r+1})\}$ with assertions $\Gamma^{\rm CS}_r\ge0$,
  $\Gamma^{\rm SC}_r\ge0$, and identity (I.4) as an exact equality.
* `front_wavenumber(field)` → $K=G_0/(2H_1^2)$, and `lemma_k_margin(field)` → $4\pi N-K$ with
  $N=\max\{|k|:e_k>0\}$.
* `renyi_half_perplexity(field)` → $e^{H_{1/2}}$ and `shannon_perplexity(field)` → $e^{S}$, asserting
  $e^{S}\le e^{H_{1/2}}$.

**Test fields (all already in repo, all exact rational):** `exact_leray_relay.build_exact_relay_triad()`;
the two strict orientations from `exact_carrier_search.search_exact_carrier_gadget()`;
`carrier_two_stage_galerkin.build_partial_carrier_parent()`; the P1/P2 mesoscopic sparse parents from
`mesoscopic_cloud_scaling.build_sparse_parent()` at $\eta=0.10,0.15,0.20$, $N=16,32,48,64$.

**Variables:** $r\in\{0,1,2\}$; $\nu\in\{1/40,1/10\}$; $c_E$ irrelevant (§S3 collapse — verify:
$\Gamma^{\rm CS}_0/(G_0/2\nu H_1)$ must be $c_E$-independent up to the explicit $\nu/\sqrt{c_E}$ in
$\Gamma^{\rm SC}$).

**Arithmetic:** everything in §F is **exact rational** — no float anywhere. Only the optional evolution
check (`mesoscopic_galerkin.run_small_mesoscopic_galerkin`, $N=4$, cutoff $3N$, RK4$\times$16) is float,
and it is used solely to confirm the *sign* of $\frac{d}{dt}\log N_0^2$ vs. (I.3), never for exponents.

**Success criterion (promotes to a strengthened monotone).** The relative saturation deficit
$\mathfrak d(u):=(\Gamma^{\rm CS}_0+\Gamma^{\rm SC}_0)\big/\big(G_0/(2\nu H_1)\big)\in[0,1]$ admits a
**scale-independent positive floor** $\mathfrak d\ge\mathfrak d_0>0$ across the $N=16\ldots64$,
$\eta=0.10\ldots0.20$ grid (monotone or flat in $N$, spread $<10\%$). Then (M) upgrades to
$\frac{d}{dt}\log N_0^2\le(1-\mathfrak d_0)\frac{KD}{2\nu}$, and §C's window sharpens to
$\gamma\ge\tfrac12$ with strict inequality forced if $\mathfrak d_0>0$ survives the $N\to\infty$ audit.

**Kill conditions.** (K1) any exact field returns $\Gamma^{\rm CS}_0<0$ or $\Gamma^{\rm SC}_0<0$, or (I.4)
fails as an exact equality ⇒ the derivation is wrong, abandon. (K2) `lemma_k_margin` $<0$ on any
band-limited field ⇒ Lemma K's lattice constant is wrong, re-derive. (K3) $\mathfrak d\to0$ along the
$N$-sweep ⇒ saturation is achievable, the monotone gives no *quantitative* improvement over Serrin and
the lane is demoted to "consistency diagnostic only" (the identity survives; the sharpening dies).

**Resolution:** none needed for the exact lane (finite trig polynomials). The float cross-check needs
grid $64^3$, cutoff $\ge3N$, dealiased, per repo `mesoscopic_galerkin` settings.

---

## G. Proof chain to Clay (B)

1. **O-1 (done here).** Identity (I.1) and gap decomposition (I.4) for finite-mode divergence-free
   zero-mean fields — exact, rational, Lean-formalizable (finite sums + Cauchy–Schwarz + square completion).
2. **O-2 (done here).** Lemma K: $K\le4\pi N+C$ for band-limited fields — finite lattice sum, Lean-able.
3. **O-3.** $H_0(T_{\max})>0$ for a blowing-up solution on $\mathbb T^3$ (backward-uniqueness type). *Open,
   standard-literature target; not needed for the dichotomy, only for (B.4).*
4. **O-4.** Extend (I.1)/(I.3) from finite trig polynomials to $u\in C([0,T);H^m)$, $m>5/2$: absolute
   convergence of $H_{r+2}$, $G_r$ and legitimacy of term-by-term differentiation. *Spectral-cutoff limit
   $T'\uparrow T$ — the repo's standard open bridge; here it is only a convergence statement, no PDE
   uniqueness needed.*
5. **O-5.** Monotone (M): $\Lambda$ non-increasing on $[0,T_{\max})$.
6. **O-6.** Dichotomy **D**: if $\int_0^{T}KD\,dt<\infty$ then $\log N_0^2$ bounded, hence
   $H_1\le CH_0\le CH_0(0)$, hence $u\in L^\infty([0,T);H^1)\subset L^\infty L^6$, hence Serrin $(6,\infty)$
   ⇒ regular at $T$. **Exactly one of:** (I) $\int_0^\infty KD\,dt<\infty$ and (B) holds for that datum;
   (II) $\int_0^{T_{\max}}KD\,dt=\infty$ with $T_{\max}<\infty$.
7. **O-7.** Tail/band obligation: replace "band-limited to $N$" in Lemma K by a Besov effective bandwidth
   with a controlled tail, so that $\int KD=\infty\Rightarrow\int ND=\infty$ holds for true solutions.
8. **O-8.** Quantitative deficit: promote the pilot's $\mathfrak d_0>0$ (if it survives) from a finite
   family to a uniform lemma over divergence-free fields at critical normalization.
9. **O-9.** Close the corner: show $\gamma=\tfrac12$ requires simultaneous saturation of (I.4) and of
   Bernstein, and that the two are incompatible for divergence-free fields (this is where the repo's
   cross-talk/phase-coherence negative results become an *ingredient*, not an obstacle).
10. **O-10.** Conclude: alternative (II) is empty ⇒ **Clay (B)**. *(O-9 is the only genuinely hard step;
    O-1/O-2 are already proofs, O-4/O-5/O-6 are routine-but-unwritten, O-3/O-7 are literature transfers.)*

---

### Binding non-claims
Nothing here proves global regularity or any Clay statement. (M) is a differential inequality with a
*critical but not yet closed* right-hand side; its finiteness is a hypothesis, not a theorem. The
$\gamma\in[\tfrac12,1)$ re-derivation is conditional on the log-critical shell law being the operative
ansatz class. No numerical output of §F may be reported as evidence for or against a singularity.
