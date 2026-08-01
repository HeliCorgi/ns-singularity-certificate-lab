# LENS 4 — Depolarization Closure of the Pressure Hessian: the Skew-Locked Vortex Ribbon

**Status: FORMAL ANSATZ** (algebraic core at NUMERICAL-CANDIDATE grade: a one-parameter family of
exact locked states verified to residual $\le 3.4\times10^{-15}$; no PDE object constructed).

---

## A. Clay target

- **Target:** TARGET-U (unforced $\mathbb R^3$ breakdown), which implies CLAY-C. Direction: *breakdown*.
- **Domain:** $\mathbb R^3$, whole space, **general Cartesian — no axisymmetry, no swirl variable, no
  $\mathcal L_5$**. This is the load-bearing choice: every registered Type-I / Liouville exclusion in the
  constraint map (CSTY2009, KNSS2009(a),(b), Ladyzhenskaya) has *axisymmetry in its hypothesis*.
- **Forcing:** none. (So F-N1/F-N2/F-N4 are vacuous; not used as an escape hatch.)
- **Initial regularity/decay:** $u_0\in C_c^\infty(\mathbb R^3;\mathbb R^3)$, $\nabla\!\cdot u_0=0$.
- **Viscosity:** $\nu>0$ fixed, treated **non-perturbatively**. It enters as one sharp inequality
  ($p_2\le 1/2$, §C) and is *not* sent to zero.

---

## B. Central mathematics

### B.1 The exact alignment system (derived)

$A=\nabla u$, $A_{ij}=\partial_j u_i$, $S=\tfrac12(A+A^{\!\top})$, $\Omega=\tfrac12(A-A^{\!\top})$,
$\Omega_{ij}=-\tfrac12\varepsilon_{ijk}\omega_k$. From $\partial_t A+(u\cdot\nabla)A=-A^2-H+\nu\Delta A$,
$H_{ij}=\partial_i\partial_j p$, and $(\Omega^2)_{ij}=\tfrac14(\omega_i\omega_j-|\omega|^2\delta_{ij})$:
$$\frac{DS}{Dt}=-S^2-\tfrac14\big(\omega\otimes\omega-|\omega|^2 I\big)-H+\nu\Delta S. \tag{1}$$
Trace of (1): $\operatorname{tr}H=\Delta p=-\operatorname{tr}A^2=2Q$, $Q:=\tfrac14|\omega|^2-\tfrac12|S|^2$.
Dimensions: $[\lambda]=T^{-1}$, $[\lambda^2]=[\omega^2]=[H]=T^{-2}$. ✔

**Exact nonlocal operator.** $p=2\Delta^{-1}Q=R_iR_j(u_iu_j)$, and with $r=x-y$,
$$\boxed{\;H_{ij}(x)=\tfrac23\delta_{ij}Q(x)-\frac1{2\pi}\,\mathrm{P.V.}\!\!\int_{\mathbb R^3}
\frac{3\hat r_i\hat r_j-\delta_{ij}}{|r|^3}\,Q(y)\,d^3y\;} \tag{2}$$
(trace of the integrand vanishes ⇒ $\operatorname{tr}H=2Q$ ✔). Global neutrality
$\int_{\mathbb R^3}Q=\tfrac12\int\partial_i\partial_j(u_iu_j)=0$: **both signs of $Q$-charge are available**
(vorticity-dominated cores $Q>0$, strain-dominated sheaths $Q<0$).

In the $S$-eigenframe $\{e_a\}$, $\lambda_a$, $c_a=\hat\omega\cdot e_a$, and $de_a/dt=W\times e_a$ with
$W_{ab}:=e_b\!\cdot\!(W\times e_a)$:
$$\frac{d\lambda_a}{dt}=-\lambda_a^2+\tfrac14|\omega|^2(1-c_a^2)-H_{aa},\qquad
(\lambda_a-\lambda_b)W_{ab}=-\tfrac14\omega_a\omega_b-H_{ab}. \tag{3}$$
Consistency: $\sum_a\!\big[-\lambda_a^2+\tfrac14|\omega|^2(1-c_a^2)-H_{aa}\big]=-|S|^2+\tfrac12|\omega|^2-2Q=0$ ✔.

### B.2 The closure: **depolarization tensor of the $Q$-blob**

If $Q$ is uniform inside an ellipsoid and zero outside, the interior of (2) is *exactly constant*:
$$H_{ij}^{\rm self}=2\,\mathcal N_{ij}\,Q,\qquad \operatorname{tr}\mathcal N=1,\qquad
\mathcal N_a=\frac{a_1a_2a_3}{2}\!\int_0^\infty\!\frac{ds}{(a_a^2+s)\sqrt{\prod_b(a_b^2+s)}}. \tag{4}$$
(Normalization forced by $\operatorname{tr}H=2Q$.) Three exact evaluations:

| shape | $\mathcal N$ | $H^{\rm self}$ |
|---|---|---|
| sphere | $(\tfrac13,\tfrac13,\tfrac13)$ | $\tfrac23 Q\,I$ — **this is exactly restricted Euler** |
| circular needle ($a_3\!\to\!\infty$) | $(\tfrac12,\tfrac12,0)$ | $\mathrm{diag}(Q,Q,0)$ |
| elliptic cylinder $(a_1,a_2)$ | $\big(\tfrac{a_2}{a_1+a_2},\tfrac{a_1}{a_1+a_2},0\big)$ | — |
| ribbon $a_2\!\ll\! a_1\!\ll\! a_3$ | $(0,1,0)$ | $\mathrm{diag}(0,2Q,0)$ |

**So: restricted Euler is not "an approximation to NS pressure"; it is the pressure Hessian of a
*spherical* $Q$-blob.** Independent check: for a $z$-invariant field $p=2\Delta_2^{-1}Q$ gives
$H_{zz}=0$ and $H_{\perp}=\delta_{ij}Q-\tfrac1\pi\mathrm{P.V.}\!\int\frac{2\hat\rho_i\hat\rho_j-\delta_{ij}}{\rho^2}Q$,
matching the needle row; naive $z$-integration of (2) reproduces the nonlocal part exactly but misses the
needle-vs-ball exclusion shift $\mathrm{diag}(-\tfrac13,-\tfrac13,\tfrac23)Q$ — the classical
depolarization subtlety, recorded here because it is the sign-critical step.

**Companion (nonlocal) term.** An infinite $Q$-line of density $\mu=\int Q\,d^2x_\perp$, direction $\hat m$,
perpendicular offset $\rho_\perp$ with unit normal $\hat w$:
$$H^{\rm ext}=\frac{\mu}{\pi\rho_\perp^2}\,\mathcal M,\qquad
\mathcal M=I-\hat m\otimes\hat m-2\,\hat w\otimes\hat w,\qquad \operatorname{tr}\mathcal M=0. \tag{5}$$
Sign control (point-source check): $p=-\mathcal Q/2\pi|x|<0$ for $\mathcal Q>0$ — a vortex core is a
*low-pressure* core ✔ — and $H_{rr}=-\mathcal Q/\pi|x|^3<0$, i.e. $-H_{rr}>0$: **a vorticity-dominated
neighbour is a strict source of stretching along the line of centres.**

### B.3 Eigenframe locking (the alignment cone)

Require the $Q$-blob principal axes to stay in the strain eigenframe. For a materially deformed ellipsoid
$\dot G=AG+GA^{\!\top}$; off-diagonal in the frame gives $(g_a-g_b)W_{ab}=(g_b-g_a)\Omega_{ab}$, hence the
**co-frame condition** $W=\omega/2$, i.e. $W_{ab}=\tfrac12\varepsilon_{abc}\omega_c$. With
$r_a=\lambda_a/|\omega|$, $\hat h_{ab}=H_{ab}/|\omega|^2$, (3) becomes the **cone conditions**
$$\text{(C1) }2(r_1-r_2)c_3=-c_1c_2-4\hat h_{12},\quad
\text{(C2) }2(r_2-r_3)c_1=-c_2c_3-4\hat h_{23},\quad
\text{(C3) }2(r_1-r_3)c_2=c_1c_3+4\hat h_{13}. \tag{6}$$

> **REJECTED V1 (restricted Euler / any isotropic closure).** $\hat h_{ab}=0$ for $a\ne b$ reduces (6) to
> $xc_3=-c_1c_2$, $yc_1=-c_2c_3$, $zc_2=c_1c_3$ with $x+y=z$; eliminating gives the **exact failing equation**
> $$c_1^2c_3^2=-c_2^2\,(c_1^2+c_3^2),$$
> whose only real solution is $c_2=0$ and $c_1c_3=0$, forcing two equal eigenvalues and then (§B.4) $r_3=0$.
> **This is the depletion, derived:** with a spherical (or any diagonal) $Q$-blob the strain eigenframe
> *cannot* co-rotate with a tilted vorticity — it precesses, the alignment cosines oscillate, and the
> Riccati source averages out. What kills restricted Euler is not the size of $H^{\rm dev}$ but its
> **absence of off-diagonal components**.

> **REJECTED V2 (aligned needle, no companion).** $c=(0,0,1)$, $\mathcal N=(\tfrac12,\tfrac12,0)$:
> co-frame forces $r_1=r_2$, and the axial strain equation reads $2r_3^2=0\Rightarrow r_3=0$. The needle
> removes the RE sink but supplies no source: exactly marginal.

> **REJECTED V4 (planar companion).** $\hat m$ in the $e_1$–$e_3$ plane ($m_2=0$) gives
> $\hat h_{12}=\hat h_{23}=0$, so (C1),(C2) force $r_1=r_2=r_3=0$. **The companion must be skew out of the
> $(\hat\omega,\text{separation})$ plane.**

### B.4 The closed system and its solution

Take $e_3$ = tube axis, $e_2$ = ribbon normal, $e_1$ = in-ribbon transverse (separation direction);
$c=(\sin\theta,0,\cos\theta)$; ribbon limit $\mathcal N=(0,1,0)$ so
$\hat h_{11}=E_1$, $\hat h_{22}=2\hat q+E_2$, $\hat h_{33}=E_3$, with $E_a=\hat g\,\mathcal M_{aa}$,
$\hat q=\tfrac14-\tfrac12\sum r_a^2$, $\sum_a E_a=0$. Self-similar ansatz $r_a,c_a,\hat g$ constant gives
$d|\omega|/dt=\hat\sigma|\omega|^2$, $|\omega|=(\hat\sigma\tau)^{-1}$, $\hat\sigma=\sum r_ac_a^2$, and
$$\text{(★}_a)\qquad r_a\hat\sigma+r_a^2-\tfrac14(1-c_a^2)+\hat h_{aa}=0 \tag{7}$$
(the sum of the three is an identity, so two are independent). With $c_2=0$, (6) *determines* the three
off-diagonals: $\hat h_{12}=-\tfrac12(r_1-r_2)\cos\theta$, $\hat h_{13}=-\tfrac14\sin\theta\cos\theta$,
$\hat h_{23}=-\tfrac12(r_2-r_3)\sin\theta$. Inverting (5) for $(\hat g,\hat m)$: with $\rho=\hat h_{12}/\hat h_{13}$,
$$\hat h_{12}=\hat g\,m_1m_2,\quad \hat h_{13}=\hat g\,m_1m_3,\quad
\hat h_{23}=-\hat g\,m_2m_3\frac{1+m_1^2}{1-m_1^2},$$
which is solvable iff $F(m_1)=\frac{1+m_1^2}{|m_1|\sqrt{1-m_1^2}}$ attains
$\big|\hat h_{23}\sqrt{1+\rho^2}/(\hat h_{13}\rho)\big|$; since $\min F=2\sqrt2$ at $m_1^2=1/3$, the
**single-companion realizability condition** is
$$\big|r_2-r_3\big|\sqrt{4(r_1-r_2)^2+\sin^2\theta}\;\ge\;2\sqrt2\,\cos\theta\,\big|r_1-r_2\big|. \tag{8}$$
Counting: 5 equations {(★$_1$),(★$_3$),(C1),(C2),(C3)}, 6 unknowns $\{r_1,r_2,\theta,\hat g,m_1,m_2\}$ ⇒ a
**one-parameter family**, confirmed numerically (continuation in $\theta$, Newton, independent re-verification
of (7) and (6) from reconstructed $S,\omega,H$; all residuals $\le3.4\times10^{-15}$):

| $\theta$ | $r_1$ | $r_2$ | $r_3$ | $\hat\sigma$ | $\hat q$ | $\hat g$ | $E_3$ |
|---|---|---|---|---|---|---|---|
| 36° | −2.0269 | −2.0160 | +4.0429 | +1.9459 | −12.009 | +24.192 | −24.125 |
| 40° | −0.7502 | −0.7112 | +1.4614 | +0.5476 | −1.3522 | +2.9228 | −2.8328 |
| **45°** | **−0.49094** | **−0.40281** | **+0.89375** | **+0.20140** | **−0.35103** | **+0.98520** | **−0.85379** |
| 48° | −0.4229 | −0.2878 | +0.7108 | +0.0847 | −0.1334 | +0.5933 | −0.4273 |
| 50° | −0.3979 | −0.2204 | +0.6183 | +0.0220 | −0.0446 | +0.4442 | −0.2491 |

Family exists for $\theta\in(\approx35^\circ,\approx50.5^\circ)$, terminating at $\hat\sigma\to0^+$.
$\hat g>0$ ⇒ **$\mu>0$: the companion is a genuine vorticity-dominated vortex line** (not a fitted sign).
Throughout, $r_3>0$, $\hat\sigma>0$, and **$E_3<0$: the nonlocal axial pressure Hessian is a strict source.**
At $\theta=45^\circ$ the axial budget splits
$$\underbrace{r_3\hat\sigma+r_3^2}_{0.979}=\underbrace{\tfrac14\sin^2\theta}_{0.125\ (13\%)}
\;+\;\underbrace{(-E_3)}_{0.854\ (87\%)},$$
i.e. **the mechanism is pressure-driven, not vortex-stretching-driven.**
Geometry at $\theta=45^\circ$: $\hat m=(-0.9913,-0.0319,0.1280)$ — the companion vortex is aimed nearly along
the separation direction, at $\rho_\perp=0.132\,d$ (admissibility caveat: it must not intersect the core).

> **Negative scan record V6.** The second root branch of (8) ($|m_1|$ small, companion nearly perpendicular
> to $e_1$) yielded **no** admissible fixed point ($\hat\sigma>0,r_3>0$) over $\theta\in[-85^\circ,85^\circ]$,
> $|r_a|\le1$. Finite-scope negative result.

---

## C. Scaling table (clock $\tau=T-t$)

The alignment fixed point forces $\|\omega\|_\infty\sim(\hat\sigma\tau)^{-1}$. Profile scales
$\ell_a\sim\tau^{p_a}$, $s:=\sum p_a$, ribbon ordering $p_2\ge p_1>p_3\ge0$, $\|u\|_\infty\sim|\omega|\ell_2$.

| quantity | exponent | constraint | at $(p_1,p_2,p_3)=(0.44,0.45,0.41)$ |
|---|---|---|---|
| $\|\omega\|_\infty$ | $\tau^{-1}$ | forced by (7) | $\tau^{-1}$ |
| $\|u\|_\infty$ | $\tau^{p_2-1}$ | Type II iff $p_2<\tfrac12$ | $\tau^{-0.55}$ |
| energy in window | $\tau^{2p_2-2+s}$ | $\ge0$ | $\tau^{+0.20}\to0$ ✔ |
| enstrophy | $\tau^{s-2}$ | — | $\tau^{-0.70}$ |
| dissipation rate $\nu\|\nabla u\|_2^2$ | $\tau^{s-2}$ | $\int<\infty\Leftrightarrow s>1$ | $\tau^{-0.70}$, $\int\!\sim\tau^{0.30}$ ✔ |
| $\|u\|_{L^3}^3$ | $\tau^{3(p_2-1)+s}$ | ESS needs $3p_2+s<3$ | $\tau^{-0.35}\to\infty$ ✔ |
| nonlinear term $\|(u\!\cdot\!\nabla)u\|_{L^2}$ | $\tau^{p_2-2+s/2}$ | — | $\tau^{-0.90}$ |
| pressure term $\|\nabla p\|_{L^2}$ | $\tau^{p_2-2+s/2}$ | Leray-bounded, same order | $\tau^{-0.90}$ |
| $\|H^{\rm dev}\|_\infty$ | $\tau^{-2}$ | $=|\omega|^2|E|$ | $\tau^{-2}$ |
| physical time remaining | $\tau$ | trivially finite | — |
| Fourier bandwidth $N$ | $\tau^{-p_2}$ | $\gamma=p_2\in(0,1)$ | $\tau^{-0.45}$ |
| active modes $M\sim(\ell_1\ell_2\ell_3)^{-1}$ | $\tau^{-s}$ | see §E | $\tau^{-1.30}$ |

**Feasible box** (nonempty): $s>1$, $s+2p_2\ge2$, $s+3p_2<3$, $p_2\le\tfrac12$. For $p_2=0.45$:
$s\in[1.10,1.65)$. Serrin pairs all diverge: $\int\|u\|_\infty^2dt\sim\int\tau^{-1.10}$,
$\int\|u\|_9^3dt\sim\int\tau^{-1.217}$, $\int\|u\|_6^4dt\sim\int\tau^{-4/3}$ — all $=\infty$.
BKM: $\int\|\omega\|_\infty dt\sim\log(1/\tau)=\infty$ ✔.

> **REJECTED V3 (frozen material blob).** If the profile is a materially frozen blob, $p_a=-r_a/\hat\sigma$,
> and along the *entire* solved family $\min_\theta p_2^{\rm mat}=1.036>\tfrac12$ (at $\theta=36^\circ$;
> $p_2^{\rm mat}=2.00$ at $45^\circ$). The thin scale would fall far below $\sqrt{\nu\tau}$ and viscosity
> would arrest it. **The mechanism therefore requires a flux-renewed (Burgers-type) profile**, whose scales
> obey the rescaled profile PDE, not material advection. This is the largest open closure step.

> **REJECTED V5.** Any realization as a *globally exact one-scale backward self-similar* solution: NRS1996 /
> Tsai1998 give $U\equiv0$. Excluded by construction here since $p_1\ne p_2\ne p_3$.

---

## D. Feedback loop (every arrow a formula)

1. $\dfrac{d|\omega|}{dt}=\hat\sigma|\omega|^2\;\Rightarrow\;|\omega|=(\hat\sigma\tau)^{-1}$.
2. $\mu=\displaystyle\int_{\rm comp}\!\!Q\,d^2x_\perp\simeq\tfrac14|\omega|^2\ell_1\ell_2$ (companion $Q$-charge grows as $|\omega|^2$).
3. $\hat g=\dfrac{\mu}{\pi\rho_\perp^2|\omega|^2}\simeq\dfrac{\ell_1\ell_2}{4\pi\rho_\perp^2}$ — **scale-free**: constant iff the dyad collapses geometrically ($\ell_1\ell_2\propto\rho_\perp^2$).
4. $E_3=\hat g\,\mathcal M_{33}<0$ (numerically $-0.854$ at $\theta=45^\circ$) ⇒ **axial source**.
5. $\dfrac{d\lambda_3}{dt}=-\lambda_3^2+\tfrac14|\omega|^2\sin^2\theta-E_3|\omega|^2$ — Riccati with positive nonlocal forcing; fixed point $r_3\hat\sigma+r_3^2=\tfrac14\sin^2\theta-E_3$.
6. $\hat\sigma=r_1\sin^2\theta+r_3\cos^2\theta$ increases with $r_3$ ⇒ **back to 1** (loop closed).
7. **Lock:** the *same* companion supplies $\hat h_{12},\hat h_{13},\hat h_{23}$ satisfying (6), holding $\theta$ constant — without it the frame precesses and the loop opens (V1).
8. **Shape gate:** $\lambda_1>\lambda_2\Rightarrow a_2/a_1\to0\Rightarrow\mathcal N\to(0,1,0)\Rightarrow\hat h_{33}^{\rm self}=2\hat q\mathcal N_3=0$: the core's own pressure is evacuated from the axial channel, leaving step 5 governed **entirely** by the nonlocal term. This is the arrow restricted Euler cuts by setting $\mathcal N_3=\tfrac13$.

**What replaces the depletion:** (i) the *local* axial Hessian is $2\mathcal N_3Q$ with $\mathcal N_3=O(a_\perp^2/L^2)\to0$, not $\tfrac23Q$ — the low-pressure core is axially uniform and exerts no axial Hessian; (ii) the *residual* axial Hessian is purely nonlocal, $E_3=\hat g\mathcal M_{33}$, whose sign is a free geometric parameter (guaranteed realizable by $\int Q=0$); (iii) the off-diagonal nonlocal Hessian, absent in every isotropic closure, is exactly what makes eigenframe locking possible.

---

## E. Obstruction audit (collision points)

1. **Energy bound / finite dissipation (F-N1/N2, Leray).** Not used as a signature. Window energy $\tau^{2p_2-2+s}\to0$; $\int\nu\|\nabla u\|_2^2dt\sim\int\tau^{s-2}d\tau<\infty$ iff $s>1$ — an *imposed* constraint of the feasible box, not an assumption.
2. **ESS $L^\infty_tL^3_x$ (U-X1).** Requires $\|u\|_3\to\infty$: satisfied iff $3p_2+s<3$ (box constraint). Collision point: $s=3-3p_2$. U-X2 demands genuine anisotropy — satisfied since $p_1\ne p_2\ne p_3$ and the shape aspect $a_2/a_1\to0$.
3. **Fixed-finite-bandwidth no-go (F-α1, VR-L-011, Lean).** $N(t)\sim\tau^{-p_2}\to\infty$: not a fixed-band ansatz. VR-L-019 requires bandwidth exponent $<1$: $\gamma=p_2\le\tfrac12<1$ ✔.
4. **Galerkin global existence.** Same escape; the object is never confined to a fixed mode set.
5. **Pure-swirl $L^3$ no-go (VR-L-016).** Requires $\partial_\theta p_0\equiv0$, i.e. axisymmetry with $u=u^\theta e_\theta$. **Exactly the class this lens abandons.** Here the pressure channel $P=3\int p\,\nabla\!\cdot(|u|u)$ is the *driver*, not an incidentally-zero term.
6. **One-scale self-similar (NRS/Tsai).** Hypothesis: exact backward self-similar profile with $U\in L^3$ or finite local energy. Collision point: would require $p_1=p_2=p_3=\tfrac12$. Excluded by the ribbon ordering. V5 records the death of the self-similar sub-variant.
7. **Mesoscopic $\gamma<1$ empty-child no-go ($D_N\le2\kappa^2\tau^2c_EM^{\rm eff}_N/N^3$).** Transcription to this geometry: shell energy $e\sim\tau^{2p_2-2+s}$ (**not** the critical wake $c_E/N$), $N\sim\tau^{-p_2}$, $M\sim\tau^{-s}$. The L-02 necessary condition $M\gtrsim\nu^2N^2/e$ becomes $\tau^{-s}\gtrsim\tau^{2-4p_2-s}$, i.e. **exactly $p_2\le\tfrac12$** — the same inequality as viscous admissibility. The $M\gtrsim N^3$ form of the floor is *specific to critical shell normalization* $E_N=c_E/N$, which this object does not have; the general form is satisfied strictly inside the box.
8. **Diagonal cross-talk gate.** Not applicable in kind: this is a physical-space alignment mechanism with no carrier alphabet, no tagged child shell, and no frozen-parent one-step relay. The analogous leakage question — how much of the strain produced at the core is delivered to the *core* rather than the sheath — is the pilot's $\hat h$ measurement (§F), not a sumset computation.
9. **CSTY2009 Type-I exclusion.** Hypotheses: (a) **axisymmetric**, (b) $|v|\le C_*|t|^{-1/2}$ or $C_*r^{-1+\varepsilon}|t|^{-\varepsilon/2}$. This candidate fails **both**: general Cartesian, and $\sqrt\tau\|u\|_\infty\sim\tau^{p_2-1/2}\to\infty$ for $p_2<\tfrac12$ (Type II). Collision point is precisely $p_2=\tfrac12$ — the candidate must sit strictly below it, which is also what the viscous scale and item 7 require. This triple coincidence at $p_2=\tfrac12$ is the sharpest structural feature of the lane.
10. **KNSS ancient-solution Liouville.** (a) requires axisymmetric *no-swirl*; (b) requires axisymmetric with swirl **and** $|u|\le C/r$. Neither holds. General bounded ancient mild solutions are open (repo: "nontrivial ancient limits retaining swirl and violating $C/r$ remain open"). Also, Type-II growth means the standard rescaling need not produce a *bounded* ancient limit at all — this is a weakening, not a proof of evasion, and is logged as obligation 8.
11. **Front-resolution threat model.** TM-22 (≥7 points per front scale) binds on $\ell_2$; TM-17 (under-converged pressure solve) is the dominant threat since the entire mechanism is the pressure Hessian — the pilot computes $H$ spectrally and cross-checks $\operatorname{tr}H=2Q$ pointwise. TM-19 (normalization) binds on the $\mathcal N$ convention, checked against the three exact rows of §B.2.
12. **CKN / local energy.** Necessary condition satisfied: $r^{-1}\!\int_{Q_r}|\nabla u|^2\sim\tau^{s-2}\cdot\tau/\tau^{p_2}=\tau^{s-1-p_2}\to\infty$ for $s<1+p_2$ — a *fourth* constraint, cutting the box to $s\in(1,1+p_2)$, i.e. $s\in(1,1.45)$ at $p_2=0.45$. Still nonempty; intersected with $s\ge2-2p_2=1.10$ gives $s\in[1.10,1.45)$.

---

## F. Minimal falsification experiment (≤1 h)

**Question tested:** in a *true* divergence-free 3D field with the skew-dyad geometry, is the deviatoric
pressure Hessian at the core (i) axially negative and (ii) large enough to dominate the vortex-tilt term?

- **Variables:** skew angle $\chi$ of the companion (grid $10^\circ$–$80^\circ$), separation $d/\ell_1\in\{3,5,8\}$, ribbon aspect $\alpha=a_2/a_1\in\{1,\tfrac12,\tfrac14,\tfrac18\}$, tilt $\theta\in\{35^\circ,40^\circ,45^\circ,50^\circ\}$, grid $128^3$ and $192^3$.
- **Procedure:** build a dyad of Gaussian elliptic vortex tubes on $\mathbb T^3$; project with `leray_response_relay.leray_project`; compute $A=\nabla u$, $S$, $\omega$, $Q$ spectrally; compute $H$ from $\hat H_k=-2\,\frac{k\otimes k}{|k|^2}\hat Q_k$ (repo Leray convention, $P_0=I$); diagonalize $S$ at the core maximum; report $\hat h_{ab}$, $\hat q$, $E_a$, and the residuals of (C1)–(C3) and (7). Advance one step with `leray_response_relay.leray_advection` (dealiased) to get $d c_a/dt$ directly.
- **Success:** $E_3<0$ **and** $|E_3|\ge\tfrac14\sin^2\theta$ **and** $|E_3|\ge3\times|\tfrac23\hat q|$ (nonlocal beats the RE prediction) for some $(\chi,d,\alpha,\theta)$ in the scan, stable between $128^3$ and $192^3$ to $<5\%$; measured $|dc_a/dt|/(\hat\sigma|\omega|)<0.2$ (approximate locking).
- **Kill:** $E_3\ge0$ everywhere in the scan, **or** $E_3<0$ but $|E_3|<\tfrac14\sin^2\theta$ everywhere (pressure channel subdominant ⇒ the mechanism is just vortex stretching and dies with V2), **or** locking residual grows under refinement.
- **Precision:** float64 for the field diagnostics; **exact rational / interval required** for (a) the depolarization identities of §B.2 (three shape rows and $\operatorname{tr}\mathcal N=1$), (b) the algebraic fixed point of §B.4 and the realizability bound $\min F=2\sqrt2$, (c) the sign of $E_3$ at the reported family point.
- **Reuse:** `leray_response_relay.{leray_project, leray_advection, spectral_inner, gradient_l2_squared}`; `mesoscopic_local_fft.local_fft_leray_coefficients` for a wrap-free cross-check of $Q$; `mesoscopic_galerkin.run_small_mesoscopic_galerkin` pattern for the one-step evolution harness; `exact_carrier_record_verifier` pattern for an independent second code path. Solver script already written: `scratchpad/lens4_solve.py`.

---

## G. Proof chain (10 obligations)

1. **(Closure)** Prove (4) rigorously as an interior identity plus a remainder bound $\|H-2\mathcal N Q\|_{L^\infty(\text{core})}\le C\,\mathrm{osc}(Q)+C'\!\int_{\rm ext}|Q|/d^3$.
2. **(Ribbon evacuation)** Prove $\mathcal N_3\le C\,a_1^2/L^2$ for a $C^1$ tube with axial modulation $L$ — the quantitative replacement for the RE depletion.
3. **(Fixed point, rigorous)** Interval-arithmetic enclosure of one solution of {(★$_1$),(★$_3$),(C1)–(C3)} with $\hat\sigma>0$, $r_3>0$, $E_3<0$, plus (8).
4. **(Profile closure)** Replace V3: derive the rescaled anisotropic profile equations (repo `future_search.md` §4 machinery, Cartesian version) and show the exponents $(p_1,p_2,p_3,s)$ land in the box of §C ∩ item 12.
5. **(Local stability)** Linearize the alignment system at the fixed point; separate gauge/neutral directions; count unstable directions (Floquet if the family is only discretely locked).
6. **(Existence of the PDE profile)** Computer-assisted fixed point for the full rescaled profile with tail/truncation bounds (PO-04/05/06/07/13).
7. **(Entry)** Connect $u_0\in C_c^\infty$ to the stable manifold (PO-09 — no strategy exists repo-wide; this is the hardest step).
8. **(Type-II certification)** Prove $\sqrt\tau\|u\|_\infty\to\infty$ rigorously, so that no Type-I exclusion (incl. future Cartesian extensions of CSTY/Seregin) applies.
9. **(Finite time)** $T-t(\tau_0)=\int_{\tau_0}^\infty(\cdot)<\infty$ (PO-10), and $\int\nu\|\nabla u\|^2dt<\infty$ consistency.
10. **(Norm divergence)** $\limsup_{t\uparrow T}\|u(t)\|_{H^m}=\infty$ via $\|u\|_{L^3}\to\infty$ + ESS contrapositive (PO-11), then PO-12/PO-14/PO-15.

**Binding non-claim.** Nothing here is evidence of a finite-time singularity. The verified content is: an
exact geometric closure of the pressure Hessian, a derived no-go for isotropic closures under eigenframe
locking, and a numerically-verified one-parameter family of locked states of the *closed model system*.
