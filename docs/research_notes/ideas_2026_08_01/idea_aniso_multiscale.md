# LENS 1 — Anisotropic Multi-Rate Collapse on $\mathbb R^3$ (sheet / ribbon / tube)

**Status: FORMAL ANSATZ.** (Sub-result §C–§E "exponent polytope + anisotropy ceiling" is a
SYMBOLIC CANDIDATE: exact rational linear algebra, machine-checkable. No profile constructed.)

---

## A. Clay target

- **Statement:** CLAY-A (global regularity on $\mathbb R^3$, unforced) — attacked in the negative;
  equivalently repo TARGET-U (unforced $\mathbb R^3$ breakdown, strictly stronger than (C)).
- **Domain:** whole space $\mathbb R^3$, Cartesian, **no periodicity, no axisymmetry, no axis**.
  Track-P/torus tooling may not be imported as certification (HS-5, U-X4).
- **Forcing:** none. Deliberate: F-N1/F-N2/F-N4 make forcing useless (energy bound is free,
  high-shell injection invisible), so forcing buys nothing here.
- **Initial data:** $u_0\in C_c^\infty(\mathbb R^3;\mathbb R^3)$, $\nabla\!\cdot u_0=0$, **not**
  axisymmetric, **not** pure swirl.
- **Viscosity:** $\nu>0$ fixed, kept exactly. It is not a perturbation: it *pins one exponent*
  (§B, eq. (V)) and is the selection principle for the thin direction.

---

## B. Central mathematics

### B.1 The solenoidal anisotropic ansatz (amplitude is forced to be anisotropic)

Try first the naive ansatz $u_j(x,t)=A(t)U_j(y)$, $y_j=x_j/L_j(t)$. Then
$\nabla\!\cdot u=A\sum_j L_j^{-1}\partial_{y_j}U_j$. For a *fixed* profile this must vanish for all
$t$ while the weights $L_j^{-1}(t)$ drift independently:
$$\frac{d}{dt}\Big(\frac{L_1^{-1}}{L_3^{-1}}\Big)\neq0 \quad\text{unless } L_1\propto L_3 .$$
**REJECTED (R2, single-amplitude ansatz):** exact failing equation
$\sum_j L_j^{-1}\partial_{y_j}U_j\equiv0$ with non-proportional $L_j$ forces $\partial_{y_j}U_j\equiv0$
$\forall j$, hence (with decay) $U\equiv0$.

The unique repair is **component-wise amplitudes locked to the lengths**:
$$\boxed{\;u_j(x,t)=C(t)\,L_j(t)\,U_j(y,s),\qquad y_j=\frac{x_j}{L_j(t)},\qquad \frac{ds}{dt}=C(t)\;}\tag{1}$$
with $[C]=\text{time}^{-1}$. Then $\nabla_x\!\cdot u=C\,\nabla_y\!\cdot U$: incompressibility transfers
*exactly*, with no leftover weights. This is the only anisotropic ansatz with that property.

### B.2 Derivation of the rescaled system

Write $\alpha=\tfrac{d\log C}{ds}$, $\beta_m=\tfrac{d\log L_m}{ds}$. From (1),
$$\partial_t u_j=C L_j\Big[(\alpha+\beta_j)C\,U_j+C\,\partial_sU_j-C\sum_m\beta_m y_m\partial_{y_m}U_j\Big].$$
Advection: $\;\sum_m u_m\partial_{x_m}u_j=\sum_m CL_mU_m\,L_m^{-1}\partial_{y_m}(CL_jU_j)=C^2L_j\,(U\!\cdot\!\nabla_y)U_j$
— **the nonlinearity is exactly isotropic in $y$**, prefactor $C^2L_j$, no residual anisotropy. This is
the structural gift of (1).

Pressure: taking $\nabla_x\cdot$ of the momentum equation,
$-\Delta_x p=\sum_{j,m}\partial_{x_j}\partial_{x_m}(u_ju_m)=C^2\sum_{j,m}\partial_{y_j}\partial_{y_m}(U_jU_m)$
(again all $L$'s cancel). Set $p=C^2L_3^2\,\Pi(y,s)$ and $\delta_j:=L_3/L_j$. Then
$$-\sum_j\delta_j^2\partial_{y_j}^2\Pi=\sum_{j,m}\partial_{y_j}\partial_{y_m}(U_jU_m)=:\mathcal S(U),\qquad
-\frac{\partial_{x_j}p}{C^2L_j}=-\delta_j^2\,\partial_{y_j}\Pi. \tag{2}$$
Viscosity: $\nu\Delta_x u_j/(C^2L_j)=\mu_3\sum_m\delta_m^2\partial_{y_m}^2U_j$ with
$\mu_m:=\nu/(C L_m^2)$ (dimensionless), $\mu_m=\mu_3\delta_m^2$.

Dividing the momentum equation by $C^2L_j$ gives the **anisotropic rescaled NS system**
$$\partial_sU_j+(\alpha+\beta_j)U_j-\sum_m\beta_m y_m\partial_{y_m}U_j+(U\!\cdot\!\nabla_y)U_j
=-\delta_j^2\partial_{y_j}\Pi+\sum_m\mu_m\partial_{y_m}^2U_j,\quad \nabla_y\!\cdot U=0. \tag{R}$$

**Anisotropic Leray projector.** Eliminating $\Pi$, with $\eta$ dual to $y$,
$$(\mathfrak P^\delta V)_j=V_j-\frac{\delta_j^2\eta_j(\eta\cdot V)}{\sum_m\delta_m^2\eta_m^2}. \tag{3}$$
*Checks.* (i) $\eta\cdot\mathfrak P^\delta V=0$ and $(\mathfrak P^\delta)^2=\mathfrak P^\delta$.
(ii) It is the physical $P_k=I-k\otimes k/|k|^2$ pulled back through $k_j=\eta_j/L_j$ — verified by
substituting $v_j=CL_jV_j$ into $P_k$ and multiplying numerator/denominator by $L_3^2$; at $\delta\equiv1$
it *is* $P_k$, matching the repo convention.
(iii) $\mathfrak P^\delta$ is **not** self-adjoint on $L^2(dy)$; it is self-adjoint for
$\langle V,W\rangle_\delta=\sum_j\delta_j^{-2}\!\int V_jW_j$, and indeed
$\int|u|^2dx=C^2L_1L_2L_3L_3^2\|U\|_\delta^2$ — the weighted inner product *is* physical energy.
Energy-neutrality of advection therefore holds in $\langle\cdot,\cdot\rangle_\delta$, not in $L^2(dy)$.

**Validation.** Put $a_j\equiv\tfrac12$ (below): (R) becomes
$\tfrac12(U+y\!\cdot\!\nabla U)+(U\!\cdot\!\nabla)U+\nabla\Pi-\mu\Delta U=0$, exactly the NRS profile
equation $-\nu\Delta U+aU+a(y\cdot\nabla)U+(U\cdot\nabla)U+\nabla P=0$ with $a=\tfrac12$. ✓

### B.3 Exponents: physical time forces $C\sim1/\tau$

$T-t=\int_s^\infty ds'/C(s')<\infty$ requires $\alpha>0$; with $\alpha$ constant, $C=C_0e^{\alpha s}$ and
$$\tau:=T-t=\frac{1}{\alpha C}\;\Longrightarrow\;\boxed{C=\frac{1}{\alpha\tau}},\qquad
L_j=\ell_j\tau^{a_j},\quad a_j:=-\beta_j/\alpha>0 .$$
Normalize $\alpha=1$ (choice of $s$): $\;\alpha+\beta_j=1-a_j$, $\;-\beta_m=a_m$, $\;\tau=e^{-s}/C_0$.
Order the axes $L_1\le L_2\le L_3$, i.e. $a_1\ge a_2\ge a_3>0$ (**1** = sheet normal / thin,
**2** = ribbon width, **3** = tube axis / long). Amplitudes $A_j=CL_j\propto\tau^{a_j-1}$.

### B.4 What viscosity forces, per direction

$\mu_m=\nu/(CL_m^2)=\nu\alpha\ell_m^{-2}\,\tau^{\,1-2a_m}$.

- $a_m<\tfrac12$: $\mu_m\to0$ — direction $m$ becomes **inviscid** in the rescaled frame.
- $a_m=\tfrac12$: $\mu_m\equiv\nu\alpha/\ell_m^2$ — **viscosity-balanced**; $L_m=\ell_m\sqrt{\tau}$ is the
  viscous scale.
- $a_m>\tfrac12$: $\mu_m\to\infty$. **REJECTED (R3):** the only balance for a bounded profile is
  $\partial_{y_m}^2U_j\to0$, which with decay in $y_m$ forces $U_j\equiv0$. Failing equation:
  $\mu_m\partial_{y_m}^2U_j=O(1)$ with $\mu_m\to\infty$.

$$\textbf{(V)}\qquad a_j\le\tfrac12\quad(j=1,2,3).$$

### B.5 Multi-rate is forced, and it degenerates the pressure

If $a_1=a_2=a_3=a$ then $\delta_j$ are constants and $L_j=c_jL(t)$; the linear map $x_j\mapsto c_jx_j$
turns (1) into $u=CL\,\tilde U(x/L)$.
- $a=\tfrac12$: $CL=\ell/(\alpha\sqrt\tau)$, i.e. an **exact backward Leray profile**.
  **REJECTED (R1a):** NRS1996 ($\tilde U\in L^3$) and Tsai1998 (finite local energy) give $\tilde U\equiv0$.
  Exact collision point: $a_1=a_2=a_3=\tfrac12\Rightarrow CL^2=\text{const}$.
- $a\in[\tfrac25,\tfrac12)$: isotropic Type II. Survives §C's constraints but is *outside this lens*
  (zero anisotropy $\Rightarrow$ no mechanism in §D; degenerate face, marked OPEN-BUT-EMPTY here).

So genuine multi-rate requires $a_1>a_3$, hence $\delta_1=L_3/L_1\propto\tau^{-(a_1-a_3)}\to\infty$, and
(R) is **non-autonomous with exponentially drifting coefficients** — precisely the class
NRS/Tsai leave open.

**Proposition (anisotropy–pressure dichotomy).** Either all $a_j$ coincide (⟹ R1a/Leray or the
isotropic degenerate face), or $\delta_1\to\infty$ and (3) converges to the *unbounded* projector
$$(\mathfrak P^\infty V)_1=V_1-\frac{\eta\cdot V}{\eta_1},\qquad (\mathfrak P^\infty V)_{2,3}=V_{2,3}.$$
Consequence: **in the limit the pressure force acts only on the thin component**,
$F_1=\partial_{y_1}^{-1}\mathcal S$, $F_{2,3}=O(\delta_1^{-2})$. The limit system is
hydrostatic/Prandtl-type: (2,3)-momentum is pressureless, $U_1$ is recovered from
$\partial_1U_1=-\partial_2U_2-\partial_3U_3$, $\Pi$ is the Lagrange multiplier.

---

## C. Scaling table (clock $\tau=T-t$), general $(a_1,a_2,a_3)$

Write $\Sigma:=a_1+a_2+a_3$, $\;e:=a_1+a_2+3a_3-2$ (energy exponent),
$\;\sigma_3:=a_1+a_2+4a_3-3=e+a_3-1$ (critical-$L^3$ exponent).
Dominant component is $u_3$ ($a_3=\min a_j$); dominant gradient is $\partial_1u_3\sim C\delta_1$.

| quantity | exponent in $\tau$ | at $a^\star=(\tfrac12,\tfrac7{16},\tfrac38)$ |
|---|---|---|
| energy $E$ | $e$ | $\tau^{1/16}\to0$ |
| enstrophy $=\|\nabla u\|_2^2$ | $e-2a_1$ | $\tau^{-15/16}\to\infty$ |
| global $\|u\|_{L^3}^3$ | $\sigma_3$ | $\tau^{-9/16}\to\infty$ |
| $\|\omega\|_\infty$ | $a_3-a_1-1$ | $\tau^{-9/8}$ |
| dissipation rate $\nu\|\nabla u\|_2^2$ | $e-2a_1$ | $\tau^{-15/16}$, $\int d\tau$ finite |
| $\|(u\!\cdot\!\nabla)u\|_2$ | $a_3-2+\Sigma/2$ | $\tau^{-31/32}$ |
| $\|\nabla p\|_2$ (component 1) | $a_3-2+\Sigma/2$ | $\tau^{-31/32}$ |
| $\|\nabla p\|_2$ (components 2,3) | $+2(a_1-a_3)$ relative | suppressed by $\tau^{1/4}$ |
| physical time remaining | $\tau=e^{-s}/C_0$ | — |
| Fourier bandwidth $N=1/L_1$ | $-a_1$ ($\gamma=a_1$) | $\tau^{-1/2}$ |
| active mode count $\mathfrak M$ | $-\Sigma$ | $\tau^{-21/16}$ |
| $\sqrt\tau\|u\|_\infty$ (Type-I test) | $a_3-\tfrac12$ | $\tau^{-1/8}\to\infty$ |

**Admissible polytope.**
$$\mathcal P=\Big\{\ \tfrac12\ge a_1\ge a_2\ge a_3>0,\ \ a_1>a_3,\ \
\underbrace{a_1+a_2+3a_3\ge2}_{\text{(E) finite energy}},\ \
\underbrace{1<\Lambda:=a_2+3a_3-a_1\le\tfrac32}_{\text{(D) finite dissipation / (CKN)}}\Big\}$$
- (D) $\Lambda>1\iff e>2a_1-1\iff\int_0^{T}\!\|\nabla u\|_2^2dt<\infty$ (automatic from (V)+(E) except
  on the face $a_1=\tfrac12,e=0$).
- (CKN) $\Lambda\le\tfrac32$: the CKN $\varepsilon$-regularity necessary condition
  $r^{-1}\!\iint_{Q_r}|\nabla u|^2\ge\varepsilon_*$ at $r=\sqrt{\tau_0}$ evaluates to
  $\tau_0^{\,1/2+e-2a_1}$, so requires $e\le 2a_1-\tfrac12$. Leray ($\Lambda=\tfrac32$) saturates it. ✓
- **Automatic consequences.** $\sigma_3\le\tfrac12+\tfrac12+2-3=0$ with equality **iff** $a=(\tfrac12,\tfrac12,\tfrac12)$;
  since $a_1>a_3$, $\;\sigma_3<0$: *every genuinely multi-rate, viscously admissible collapse has
  $\|u\|_{L^3}\to\infty$ for free.* Likewise $\int\|\omega\|_\infty dt\propto\tau^{a_3-a_1}\to\infty$.
- $\Sigma$-rigidity: $a_2+3a_3\le4a_1$ and (E) give $a_1\ge\tfrac25$, strict; so $a_1\in(\tfrac25,\tfrac12]$.

**Theorem (anisotropy ceiling).** On $\mathcal P$, $\;a_1-a_3<\tfrac16$; the supremum is approached only at
$(\tfrac12,\tfrac12,\tfrac13)$, where $\Lambda=1$ and total dissipation diverges logarithmically.
*Proof.* (E) with $a_2\le a_1$ gives $a_3\ge(2-2a_1)/3$, so $a_1-a_3\le(5a_1-2)/3\le\tfrac16$. ∎
Consequence: aspect ratios can grow at most like $\tau^{-1/6}$ — **weak** anisotropy is all that
finite energy permits. This is a hard budget, not a modelling choice.

**Chosen vertex** $a^\star=(\tfrac12,\tfrac7{16},\tfrac38)$: $\Lambda=\tfrac{17}{16}\in(1,\tfrac32]$,
(E)$=\tfrac{33}{16}\ge2$, gap $a_1-a_3=\tfrac18$ (78% of the ceiling), three *distinct* rates.
Geometry: viscous sheet thickness $\ell_1\sqrt{\nu\tau}$, ribbon width $\tau^{7/16}$, tube length
$\tau^{3/8}$; $L_2/L_1=\tau^{-1/16}$, $L_3/L_2=\tau^{-1/16}$. Direction 1 is viscosity-balanced
($\mu_1$ const); directions 2,3 are asymptotically inviscid ($\mu_2\propto\tau^{1/8}$, $\mu_3\propto\tau^{1/4}$).

**All tracked Serrin pairs diverge** ($\|u\|_p\propto\tau^{(a_3-1)+\Sigma/p}$, exponent of
$\int\|u\|_p^qdt$): $(\infty,2)\!:\!-\tfrac54$, $(9,3)\!:\!-\tfrac{23}{16}$, $(6,4)\!:\!-\tfrac{13}8$,
$(4,8)\!:\!-\tfrac{19}8$, $(3,\infty)$ via $\sigma_3<0$ — all $<-1$. Vorticity endpoint $(\infty,1)\!:\!-\tfrac98$.

---

## D. The closed feedback loop

Dominant physical vorticity: $\omega_2=\partial_3u_1-\partial_1u_3=-C\delta_1\,\Omega+O(C\delta_1^{-1})$ with
$\Omega:=\partial_{y_1}U_3$ — a **vortex ribbon**: normal $e_1$, vorticity $\parallel e_2$, velocity jump
$\parallel e_3$. Differentiating (R)$_3$ in $y_1$ and using $\partial_1U_1=-(\partial_2U_2+\partial_3U_3)$
(the $-(\partial_1U_1)\Omega$ and $-\Omega\partial_3U_3$ terms partially cancel):
$$\boxed{\;\partial_s\Omega+\underbrace{(1-a_3+a_1)}_{\text{dilation drain}}\Omega
+\sum_m a_my_m\partial_{y_m}\Omega+(U\!\cdot\!\nabla_y)\Omega
=\underbrace{(\partial_{y_2}U_2)\,\Omega}_{\text{ribbon stretching}}
-\underbrace{(\partial_{y_1}U_2)(\partial_{y_2}U_3)}_{\text{tilting}}
+\mu_1\partial_{y_1}^2\Omega\;}\tag{4}$$

Loop, every arrow a formula:

1. $\Omega\;\Rightarrow\;\omega_2=-C\delta_1\Omega$, $\;\|\omega\|_\infty\propto\tau^{a_3-a_1-1}$ (drives BKM).
2. $\Omega\;\Rightarrow\;\mathcal S(U)=\sum\partial_j\partial_m(U_jU_m)\;\Rightarrow\;\Pi=-\delta_1^{-2}\partial_1^{-2}\mathcal S$
   (eq. 2, degenerate limit) $\Rightarrow$ pressure force $F_1=\partial_1^{-1}\mathcal S$ only.
3. $F_1\;\Rightarrow\;U_1$ via $\partial_1U_1=-(\partial_2U_2+\partial_3U_3)$: the thin-direction inflow.
4. $U_1\;\Rightarrow$ transports $\Omega$ into the sheet, sharpening $\partial_{y_1}$; balanced against
   $\mu_1\partial_1^2\Omega$ at $L_1=\ell_1\sqrt{\nu\tau}$ — this **is** the constraint $a_1=\tfrac12$.
5. Ribbon strain $\partial_2U_2>0$ feeds $\Omega$ at rate $\langle\partial_2U_2\rangle_{\Omega^2}$.
6. Closure: multiplying (4) by $\Omega$, integrating (advection contributes $0$ by $\nabla_y\!\cdot U=0$;
   dilation gives $+\tfrac12\Sigma\|\Omega\|_2^2$),
$$\frac{d}{ds}\log\|\Omega\|_2=\underbrace{-\Big(1-a_3+a_1-\tfrac{\Sigma}{2}\Big)}_{=-15/32\ \text{at}\ a^\star}
+\underbrace{\frac{\int(\partial_2U_2)\Omega^2-\int\Omega(\partial_1U_2)(\partial_2U_3)-\mu_1\|\partial_1\Omega\|_2^2}{\|\Omega\|_2^2}}_{=:Q(s)}.$$
$$\textbf{(K) Closure condition:}\qquad Q(s)\ \ge\ 1-a_3+a_1-\tfrac{\Sigma}{2}=\tfrac{15}{32}\ \text{ at }a^\star.$$
(K) is this lens's analogue of the repo's "positive flux margin" — but it is a **local strain
functional**, not a Fourier flux, so the Fourier relay no-gos do not bound it.

**REJECTED (R4, quasi-2D depletion).** If $\partial_{y_3}\equiv0$ the system is 2.5D: $(U_1,U_2)$ is a
2D flow and $U_3$ a passively advected scalar; 2.5D NS is globally regular, so the loop cannot close.
Failing equation: with $\partial_3\equiv0$, (4) reduces to 2D scalar-gradient stretching, and
$L_2=L_3$ makes $\delta_2=\delta_3$, killing the three-rate structure. **Design rule:** three *distinct*
rates $a_1>a_2>a_3$ and a monitored ratio $\|\partial_3U\|_2/\|\partial_1U\|_2$ bounded below.

**Honest hazard.** The limit system is Prandtl-type and is ill-posed in Sobolev without a
monotonicity/Oleinik condition. The selector is $\partial_{y_1}U_3\neq0$ across the layer (exactly the
shear-ribbon structure posited). Ill-posedness of the limit is simultaneously the growth mechanism
and the main analytic obstruction; it is *not* claimed resolved.

---

## E. Obstruction audit (collision point → evasion)

1. **Energy bound / finite dissipation (F-N1/N2, Leray).** Collision: $E\propto\tau^{e}$, $e=\tfrac1{16}\ge0$;
   $\int\|\nabla u\|_2^2dt\propto\tau^{e-2a_1+1}$, finite since $\Lambda>1$. Evasion: these are *satisfied*,
   not violated; the signature used is $\|u\|_{L^3}$ and $\int\|\omega\|_\infty dt$, never energy (O-FE clear).
2. **ESS $L^\infty_tL^3_x$ (U-X1) / anisotropic $L^3$ identity (U-X2).** Collision:
   $\|u\|_3^3=A^3L_1L_2L_3\|\tilde U\|_3^3$ with $A=CL_3$, $\tilde U_j=\delta_j^{-1}U_j$, exponent
   $\sigma_3=a_1+a_2+4a_3-3$. Isotropic-parabolic gives $\sigma_3=0$ (the U-X2 identity $\equiv1$).
   Evasion: $\sigma_3<0$ *strictly and automatically* whenever $a_1>a_3$ (§C) — this is exactly the
   "genuinely anisotropic rates" requirement U-X2 names.
3. **One-scale self-similar (NRS1996 / Tsai1998).** Collision: (R) with $\delta_j,\mu_m$ constant is the
   NRS profile equation (verified in §B.2). Evasion: $\delta_1(s)=\delta_1(0)e^{(a_1-a_3)s}\to\infty$ and
   $\mu_2,\mu_3\to0$ — (R) admits **no** $s$-independent solution. Sub-variant R1a records the exact
   failing equation for the constant-aspect-ratio case.
4. **CSTY2009 axisymmetric Type-I exclusion.** Double evasion. (i) Hypothesis fails: the ansatz is
   fully 3D Cartesian, $a_1\neq a_2$, no rotational symmetry, no axis; (2.1)/(2.2) are axisymmetric
   statements. (ii) Even if made axisymmetric, $\sqrt\tau\|u\|_\infty\propto\tau^{a_3-1/2}=\tau^{-1/8}\to\infty$
   globally (not merely locally), so bound (2.1) is violated in the exterior too. Type-I would force
   $a_3=\tfrac12$, hence $a_j\equiv\tfrac12$, hence R1a.
5. **KNSS2009 Liouville / ancient solutions.** Collision: Prop 6.1 extracts a *bounded* ancient mild
   solution from parabolic rescaling at the blowup point. Evasion: parabolic rescaling of a Type-II
   field gives $\lambda\|u\|_\infty\big|_{\lambda=\sqrt\tau}\propto\tau^{a_3-1/2}\to\infty$ — the extraction
   fails at its first step. Anisotropic rescaling *does* converge, but $(x_j\mapsto x_j/L_j)$ is not an
   NS symmetry, so the limit is (R), not an ancient NS solution. (Seregin Type-II filters: design
   warning only, per repo policy.)
6. **Pure-swirl $L^3$ no-go (VR-L-016 / LG-9).** Collision: LG-9 uses rotational equivariance of the
   Riesz pressure ⟹ $\partial_\theta p_0\equiv0$ ⟹ $P\equiv0$. Evasion: no axisymmetry; the energy
   moment $M_{jm}=\int u_ju_m$ is $\mathrm{diag}(\tau^{2a_1},\tau^{2a_2},\tau^{2a_3})\cdot(\dots)\neq c\,I$,
   so the quadrupole $p_0=(4\pi|x|^3)^{-1}(3M_{jm}n_jn_m-\mathrm{tr}M)\not\equiv0$ and the pressure
   channel is open. (Also: the rank-one moment rejection does not apply — $M$ here is full-rank.)
7. **Fixed-finite-bandwidth / Galerkin no-go (F-6, F-α1, VR-L-011).** Collision: fixed-band trajectories
   are provably global. Evasion: $N(t)=1/L_1\propto\tau^{-1/2}\to\infty$; mode count
   $\mathfrak M\propto\tau^{-21/16}\to\infty$. Repo shell exponent $\gamma=a_1=\tfrac12\in(0,1)$, consistent
   with VR-L-019 ($\gamma<1$ forced).
8. **Mesoscopic $\gamma<1$ empty-child no-go, $D_N\le2\kappa^2\tau^2c_EM^{\rm eff}_N/N^3$.** This is the
   sharpest collision and must be stated exactly. Our Fourier support is a single anisotropic box of
   $\sim(L_1L_2L_3)^{-1}$ modes, so $M^{\rm eff}/N^3\lesssim L_1^2/(L_2L_3)=\tau^{2a_1-a_2-a_3}=\tau^{3/16}\to0$:
   **if the bound applied, it would kill this candidate.** It does not apply because its hypotheses fail
   on three named points: (i) *critical shell normalization* $E_N=c_E/N$ across octaves — our spectrum is a
   single peaked anisotropic band, not a dyadic critical wake; (ii) *empty child*, $E_{\rm child}(0)=0$ —
   the band at $|k|\sim1/L_1(t)$ is continuously populated by the trajectory's own smooth tail, so the
   ratio $E_{\rm child}(\tau N^{-2})/E_{\rm parent}(0)$ is not the relevant quantity, the *logarithmic
   growth rate* of an already-warm band is; (iii) *frozen parent over one parabolic time* — here the
   parent deforms on the same $s$-scale. The bound is a one-stage cold-start Duhamel estimate; this
   mechanism is a continuous deformation, not a relay.
9. **Diagonal cross-talk gate.** Not applicable: no discrete carrier alphabet, no second-stage relay.
10. **Smooth-forcing high-frequency decay (F-N4).** Not applicable: unforced.
11. **CKN / LEI.** Not evaded — *used*: $r^{-1}\iint_{Q_r}|\nabla u|^2\propto\tau_0^{1/2+e-2a_1}\to\infty$
    satisfies $\ge\varepsilon_*$, and this becomes the upper constraint $\Lambda\le\tfrac32$ in $\mathcal P$.
    Singular set = one point, parabolic dimension 0 ✓.
12. **Front-resolution threat model (TM-22, $\ge7$ pts/front; TM-01/02/16).** The pilot runs in the
    *rescaled* frame, where the front thickness is $O(1)$ in $y_1$ by construction, so the resolution
    requirement is $s$-independent — the standard dynamic-rescaling answer. TM-16 (remesh error) is
    replaced by exponent drift, monitored via the gauge residuals of §F.
13. **Track-P/HS-5 conflation.** The pilot's periodic $y$-box is a *screening* device only; no
    whole-space claim may be derived from it (U-X4, HS-5 explicitly open).

---

## F. Minimal falsification experiment ($\le1$ hour)

**Object of test: the closure inequality (K), not blowup.**

- **Variables.** Exponents $a$ (exact rationals, default $a^\star$); seed profile $U^{(0)}$ = a
  divergence-free localized ribbon, $U_3=\Phi(y_1)G(y_2,y_3)$ with $\Phi$ an odd monotone shear layer
  ($\Phi'>0$: Oleinik selector), $U_2$ a compensating strain field with $\partial_2U_2>0$ on the core,
  $U_1$ recovered from incompressibility; $\mu_1\in\{0.02,0.05,0.1\}$; box $[-\pi,\pi]^3$ with a
  Fejér window; resolution $192\times128\times96$ (anisotropic, $\ge7$ points across the $y_1$ layer).
- **Procedure.** (a) exact-rational check that $a\in\mathcal P$ (linear algebra over $\mathbb Q$);
  (b) evaluate $Q(0)$ for the seed; (c) integrate (R) in $s$ with drifting $\delta_j(s),\mu_m(s)$ for
  $s\in[0,4]$ (RK4, dealiased), logging $Q(s)$, $\|\Omega\|_2$, shape drift
  $\rho(s)=1-|\langle U(s),U(0)\rangle|/(\|U(s)\|\|U(0)\|)$, and $\|\partial_3U\|_2/\|\partial_1U\|_2$.
- **Success criterion.** $Q(s)\ge\tfrac{15}{32}$ for all $s\in[1,4]$ with $\|\Omega(s)\|_2$
  non-decreasing and $\rho(s)$ decreasing (approach to a slow manifold), for at least one $(\mu_1$, seed$)$.
- **Kill conditions (pre-registered).** (i) $\sup_{s\ge1}Q<\tfrac{15}{32}$ across the whole scan;
  (ii) $\|\partial_3U\|_2/\|\partial_1U\|_2\to0$ (R4 depletion fires); (iii) $Q$ decreasing under
  resolution refinement (numerical origin); (iv) $\|\Omega\|_2$ monotone decreasing for all seeds.
- **Arithmetic.** *Exact rational*: membership in $\mathcal P$, all exponents in §C, and $Q(0)$ when the
  seed is a finite trigonometric polynomial. *Interval*: $Q(0)$ margin against $15/32$. *Float*: the
  $s$-evolution (screening only; never promotable past "numerical observation").
- **Repo modules to reuse.** `leray_response_relay.leray_advection` / `leray_project` (dealiased
  $\mathbb P((u\cdot\nabla)u)$ core — must be *generalized* to $\mathfrak P^\delta$ of eq. (3); at
  $\delta\equiv1$ it must reproduce the existing function bit-for-bit, which is the unit test);
  `mesoscopic_local_fft.local_fft_leray_coefficients` (zero-padded exact linear convolution, no wrap
  aliasing) for exact $Q(0)$; `modal_front_actions.modal_growth_identity` for exact rational
  $H_r,T_r,N_r^2$ bookkeeping of the $\Omega$ budget; `exact_carrier_record_verifier` as the
  independent second code path (TM-14).

---

## G. Proof chain (10 obligations)

1. **O1.** $a^\star\in\mathcal P$ in exact rational arithmetic; anisotropy-ceiling theorem formalized. *(cheap; Lean-able)*
2. **O2.** (R)+(3) is *equivalent* to NS under (1) for smooth decaying fields — round-trip identity, and
   the smoothness of the reconstructed Cartesian field (PO-01, PO-12: not a coordinate artifact).
3. **O3.** Well-posedness of the non-autonomous rescaled system on $s\in[s_0,\infty)$ under the Oleinik
   selector $\partial_{y_1}U_3>0$ (Prandtl-type; **the hard analytic core**).
4. **O4.** Existence of a bounded, localized, $s$-recurrent (slow-manifold) solution $U(\cdot,s)$ of (R)
   with $\liminf\|\Omega(s)\|_2>0$ — via closure (K) with interval-arithmetic constants (PO-04, PO-07).
5. **O5.** Uniform $y$-decay / tail bounds so that the norms of §C are the true norms of $u$ (PO-06, PO-07).
6. **O6.** $T-t(s)=\int_s^\infty ds'/C(s')<\infty$ with rigorous rate bounds (PO-10).
7. **O7.** $\limsup_{t\uparrow T}\|u(t)\|_{L^3}=\infty$ and $\int_0^T\|\omega\|_\infty dt=\infty$, hence
   non-extendability in $C([0,T);H^m)$, $m>5/2$ (PO-11).
8. **O8.** Linearized/Floquet analysis of (R) around the recurrent object: finite unstable dimension,
   gauge modes separated (PO-08).
9. **O9.** Entry: an explicit $u_0\in C_c^\infty$ whose trajectory meets the stable manifold (PO-09 —
   currently strategy-less repo-wide; the true bottleneck).
10. **O10.** Interval certificate + independent reimplementation + formalization (PO-13/14/15) ⟹
    CLAY-A resolved negatively / TARGET-U.

---

### Sub-variants killed during derivation (kept, per repo policy)

| id | variant | exact failing equation |
|---|---|---|
| **R1a** | constant aspect ratio, $a_j\equiv\tfrac12$ | $CL^2=\text{const}\Rightarrow u=(2a\tau)^{-1/2}\tilde U(x/\sqrt{2a\tau})$; NRS1996/Tsai1998 $\Rightarrow\tilde U\equiv0$ |
| **R1b** | isotropic $a_j\equiv a\in[\tfrac25,\tfrac12)$ | survives $\mathcal P$ but $\delta_j\equiv1$ ⟹ no ribbon strain in (4); degenerate face, empty of mechanism |
| **R2** | single amplitude $u_j=A(t)U_j(x_j/L_j)$ | $\sum_jL_j^{-1}\partial_{y_j}U_j\equiv0\ \forall t$ with non-proportional $L_j\Rightarrow U\equiv0$ |
| **R3** | any $a_m>\tfrac12$ | $\mu_m=\nu\alpha\ell_m^{-2}\tau^{1-2a_m}\to\infty$, balance forces $\partial_{y_m}^2U_j\equiv0$, decay $\Rightarrow U_j\equiv0$ |
| **R4** | $\partial_{y_3}\equiv0$ (2.5D sheet) | 2.5D NS globally regular; (4) degenerates to 2D passive-scalar gradient stretching |
| **R5** | max-anisotropy vertex $(\tfrac12,\tfrac12,\tfrac13)$ | $\Lambda=1$ exactly $\Rightarrow\int_0^T\|\nabla u\|_2^2dt\propto\int_0\tau^{-1}d\tau=\infty$, violates Leray |
