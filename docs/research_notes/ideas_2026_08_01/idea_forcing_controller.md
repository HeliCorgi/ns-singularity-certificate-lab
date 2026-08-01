# LENS 9 — The band-limited sea controller: legalizing Clay (D) with a two-block forced-sea / autonomous-front architecture

**STATUS: FORMAL ANSATZ.** Parts B.1–B.5 (controller admissibility, exact low-return identity, sea Grönwall, seed controllability) are **PROOF CANDIDATE** — finite-Fourier + classical, Lean-able. Part D (front recurrence) inherits the FORMAL ANSATZ status of the cloud lane and is **not** improved here.

---

## A. Clay target

**Statement (D)**: breakdown on \(\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3\), smooth forcing allowed. Conventions: \(u=\sum_k\hat u_ke^{ik\cdot x}\), \(\|u\|_2^2=\sum_k|\hat u_k|^2\), \(\partial_t\hat u_k=\mathcal N_k-\nu|k|^2\hat u_k+\hat f_k\), \(\mathcal N_k=-iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m\), \(P_k=I-k\otimes k/|k|^2\), \(P_0=I\).

**Forcing: YES**, and *strictly band-limited*: \(P_{>K_0}f\equiv0\) for one fixed \(K_0\) (eq. (10.1) of the cloud candidate). Admissibility is proved in B.2, not assumed.
**Initial data**: \(u_0\equiv0\) is admissible and is the recommended choice (see B.5). Any \(u_0\in C^\infty(\mathbb T^3)\), \(\nabla\!\cdot\!u_0=0\), also works.
**Viscosity**: \(\nu>0\) fixed, never sent to \(0\); it appears with the correct sign in every estimate below and is *used* (it is what damps the sea).

---

## B. Central mathematics

### B.1 The exact low-return identity (the load-bearing new computation)

Let \(B(u,u)=\mathbb P((u\cdot\nabla)u)\), so \(\widehat{B}(k)=-\mathcal N_k=iP_k\sum_{\ell+m=k}(m\cdot\hat u_\ell)\hat u_m\). Put \(m=k-\ell\). Incompressibility gives \(\ell\cdot\hat u_\ell=0\), hence \((k-\ell)\cdot\hat u_\ell=k\cdot\hat u_\ell\) **exactly**, and

\[
\boxed{\;\widehat{B}(k)=iP_k\sum_{\ell}(k\cdot\hat u_\ell)\,\hat u_{k-\ell}\;}
\tag{9.1}
\]

*Dimension check*: \([k][\hat u]^2\), matching \(\widehat{(u\cdot\nabla)u}\). *Sign check*: the substitution is an identity, not an estimate; no cancellation is assumed. *Special case* \(k=0\): \(\widehat B(0)=0\) exactly — the mean mode receives nothing from the nonlinearity, so \(\partial_t\hat u_0=\hat f_0\) is an exactly disturbance-free integrator.

(9.1) makes the \(|k|\) prefactor **exact** rather than Bernstein-derived. Cauchy–Schwarz on \(\ell\):

\[
|\widehat B(k)|\le|k|\;\|u\|_2^2 ,\qquad
\|P_{\le K_0}B(u,u)\|_2\le\Big(\!\!\sum_{0<|k|\le K_0}\!\!|k|^2\Big)^{1/2}\|u\|_2^2
\le\sqrt{\tfrac{8\pi}{5}}\,K_0^{5/2}\,\|u\|_2^2 .
\tag{9.2}
\]

Applied to one front octave at critical normalization \(\|u_N\|_2^2=2c_E/N\):

\[
\boxed{\;\|P_{\le K_0}B(u_N,u_N)\|_2\;\le\;\sqrt{\tfrac{8\pi}{5}}\;\frac{2c_E K_0^{5/2}}{N}\;}
\tag{9.3}
\]

**Unconditional, phase-independent, no cancellation lemma, no mode count.** Summed over the doubling ladder \(N_j=2^jN_0\): \(\sum_j\|P_{\le K_0}B(u_{N_j},u_{N_j})\|_2\le 2\sqrt{8\pi/5}\,(2c_EK_0^{5/2}/N_0)\) — a closed-form constant. Under sliver phase-incoherence the same computation with \(M_N=c_M(\eta N)^3\) modes gives the sharper \(\;\asymp c_Ec_M^{-1/2}\eta^{-3/2}(K_0/N)^{5/2}\).

### B.2 The controller and verification of (10.2)

Prescribe a **sea path** \(S^*(t)\in V_{K_0}:=\mathrm{span}\{e^{ik\cdot x}:0\le|k|\le K_0\}\cap\{\text{div-free, real}\}\), chosen as a *trigonometric polynomial in \(t\)* with base frequency \(\Omega\) and \(Q\) harmonics, multiplied by a smooth cutoff \(\chi(t)\equiv1\) on \([0,T+\delta]\), \(\chi\) Schwartz-decaying after. Define the **open-loop residual controller**

\[
\boxed{\;f(t):=\partial_tS^*(t)-\nu\Delta S^*(t)+P_{\le K_0}\mathbb P\big((S^*\!\cdot\!\nabla)S^*\big)(t)\;}
\tag{9.4}
\]

(i) *Band-limitation*: every term is \(P_{\le K_0}\) of a \(V_{K_0}\)-field except the quadratic one, which carries an explicit \(P_{\le K_0}\). So \(P_{>K_0}f\equiv0\), i.e. (10.1). ✔
(ii) *(10.2)*: on \(V_{K_0}\), \(\|g\|_{H^s}\le(1+K_0^2)^{s/2}\|g\|_2\) — **band-limitation collapses the infinite family of Clay smoothness conditions to one \(L^2\) bound per time-derivative order.** Since \(f\) is a trigonometric polynomial in \(t\) times \(\chi\),

\[
\sup_t\|\partial_t^qf(t)\|_{H^s}\le C_\chi\,(Q\Omega)^q\,(1+K_0^2)^{s/2}\,\sup_t\|f\|_2<\infty
\quad\forall q,s\ge0,
\tag{9.5}
\]

and \(\chi\) supplies Clay's \((1+|t|)^{-K}\) decay for free **because \(T<\infty\)**. ✔ This closes cloud-candidate lemma #7.

### B.3 REJECTED sub-variant V1 — exact-residual pinning (closed loop)

The natural controller *pins* the sea: choose \(f\) so that \(P_{\le K_0}u\equiv S^*\), i.e.
\(f=\partial_tS^*-\nu\Delta S^*+P_{\le K_0}B(u,u)\), reading the **true** \(u\). It satisfies (10.1) and, by (9.3), is bounded in every \(H^s\). **It fails (10.2) at \(q=1\).** Differentiating (9.1) and keeping the viscous part of \(\partial_t\hat u_\ell\) for the near-antipodal pairs \(|\ell|\simeq|k-\ell|\simeq N\):

\[
\|\partial_tf\|_2\;\gtrsim\;2\nu N(t)^2\,\|P_{\le K_0}B(u_N,u_N)\|_2
\;\asymp\;4\sqrt{\tfrac{8\pi}{5}}\,\nu c_EK_0^{5/2}\,N(t)\;\xrightarrow[t\uparrow T]{}\;\infty .
\tag{9.6}
\]

**REJECTED by (9.6).** Same failure for proportional feedback \(f=\mu(S^*-S)\) (wall at \(q=2\)) and for the first-order filtered controller \(\dot z=-\lambda z+\lambda(S^*-S),\,f=\mu z\) (wall at \(q=3\)): any controller that reads \(u(t)\) inherits \(\partial_t^qu\), which is unbounded. **Conclusion: the controller must be open-loop.** This is the single most important structural finding of this lens.

### B.4 Two-block decomposition and the interface estimate

Write \(u=S+\Phi\), \(S=P_{\le K_0}u\) (**sea**, fixed finite band), \(\Phi=P_{>K_0}u\) (**front + wake**, bandwidth \(\to\infty\)). Since \(P_{>K_0}f=0\):

\[
\partial_tS=-P_{\le K_0}B(u,u)+\nu\Delta S+f,\qquad
\partial_t\Phi=-P_{>K_0}B(u,u)+\nu\Delta\Phi .
\tag{9.7}
\]

**The front block is exactly unforced.** Set \(e:=S-S^*\) and \(d(t):=P_{\le K_0}\big[B(u,u)-B(S^*,S^*)\big]-\text{(bilinear-in-}e\text{ terms)}\); the pure-front part of \(d\) is bounded by (9.3). Subtracting (9.4) from (9.7),

\[
\partial_te=-P_{\le K_0}\mathbb P\big[(S^*\!\cdot\!\nabla)e+(e\!\cdot\!\nabla)S^*+(e\!\cdot\!\nabla)e\big]+\nu\Delta e-d .
\tag{9.8}
\]

Test with \(e\) (all terms exact):
* \(\langle P_{\le K_0}\mathbb P[(S^*\!\cdot\!\nabla)e],e\rangle=\int(S^*\!\cdot\!\nabla)e\cdot e=0\) (\(e\in V_{K_0}\), \(\mathbb Pe=e\), \(\nabla\!\cdot\!S^*=0\));
* \(\langle P_{\le K_0}\mathbb P[(e\!\cdot\!\nabla)e],e\rangle=0\) **exactly** (repo identity F-12, energy neutrality on a finite band);
* \(|\langle P_{\le K_0}\mathbb P[(e\!\cdot\!\nabla)S^*],e\rangle|\le\Lambda_*\|e\|_2^2\), \(\Lambda_*:=\|\nabla S^*\|_{L^\infty}\);
* \(\langle\nu\Delta e,e\rangle=-\nu\|\nabla e\|_2^2\le-\nu\|e\|_2^2\) (zero-mean part; \(|k|\ge1\) on \(\mathbb T^3\) — the torus spectral gap, unavailable on \(\mathbb R^3\)).

\[
\boxed{\;\frac{d}{dt}\|e\|_2\le-(\nu-\Lambda_*)\|e\|_2+\|d\|_2
\;\Longrightarrow\;
\limsup_{t\uparrow T}\|e(t)\|_2\le\frac{\sup_t\|d\|_2}{\nu-\Lambda_*}\;}
\tag{9.9}
\]

whenever \(\Lambda_*<\nu\) (design condition on \(S^*\): amplitude \(\lesssim\nu/K_0\)). By (9.3), \(\|d\|_2\lesssim c_EK_0^{5/2}/N(t)\to0\). **Interface estimate: the sea tracks its target with error \(O(N^{-1})\to0\), and the force never touches \(k>K_0\).** The front's boundary condition — the ambient strain \(\nabla S\) and low-mode phase it sees — converges to the prescribed \(\nabla S^*(t)\) as the cascade descends. This is item (iv) of the lens, closed.

### B.5 Seed preparation (PO-09) by full actuation

On \(V_{K_0}\) the control operator is the identity: the system is **fully actuated**. Given any smooth path \(t\mapsto x(t)\in V_{K_0}\) with \(x(0)=P_{\le K_0}u_0\) and \(x(t_1)=\) the desired seed, (9.4) with \(S^*=x\) drives the *truncated* sea exactly along it. Exact controllability of the low block is therefore **trivial**, and by (9.5) the required control is Clay-admissible. With \(u_0\equiv0\) and \(\Phi(0)=0\), (9.9) gives \(\|e\|\le\sup\|d\|/(\nu-\Lambda_*)\) with \(d\) small while \(\Phi\) is small: the seed is reached to relative error \(O(\|\Phi\|^2K_0^{5/2}/(\nu-\Lambda_*))\).

**Which stages the controller can manufacture outright.** A front octave at \(N_j\) with relative width \(\eta\) occupies \(|k|\in[(1-\eta)N_j,(1+\eta)N_j]\); its forward child sits at \(\le(2+2\eta)N_j\). Hence the controller *directly synthesizes* every stage with

\[
N_j\le\frac{K_0}{2+2\eta},\qquad\text{i.e. } j\le J_{\rm ctrl}:=\Big\lfloor\log_2\frac{K_0}{(2+2\eta)N_0}\Big\rfloor .
\tag{9.10}
\]

Its **full** difference band \(|k|\le2\eta N_j\) is inside the controlled zone iff \(N_j\le K_0/(2\eta)\). For \(N_0>K_0/(2+2\eta)\), \(J_{\rm ctrl}=0\): the force prepares the basin only, and every cascade stage is autonomous.

---

## C. Scaling table

Front \(N(t)=(2a(T-t))^{-1/2}\) in the default \(\gamma=1/2\) branch; \(\tau=T-t\); free knob \(a(s)\) allows \(N\asymp\tau^{-\gamma}\), \(\gamma\in[1/2,1)\). Physical core: diameter \(N^{-1}\), amplitude \(N\).

| quantity | in \(N\) | in \(\tau\) (\(\gamma=1/2\)) | general \(\gamma\) |
|---|---:|---:|---:|
| total energy | \(O(1)\) | \(\tau^0\) | \(\tau^0\) |
| front-octave energy \(c_E/N\) | \(N^{-1}\) | \(\tau^{1/2}\) | \(\tau^{\gamma}\) |
| enstrophy | \(N\) | \(\tau^{-1/2}\) | \(\tau^{-\gamma}\) |
| global \(\|u\|_3^3\) | \(\log N\) | \(\log(1/\tau)\) | \(\log(1/\tau)\) |
| \(\|\omega\|_\infty\) | \(N^2\) | \(\tau^{-1}\) | \(\tau^{-2\gamma}\) |
| dissipation rate \(\nu\|\nabla u\|_2^2\) | \(\nu c_EN\) | \(\tau^{-1/2}\) | \(\tau^{-\gamma}\) (integrable iff \(\gamma<1\)) |
| nonlinear \(\|B\|_2\) | \(2c_E\eta^{3/2}N^{3/2}\) | \(\tau^{-3/4}\) | \(\tau^{-3\gamma/2}\) |
| pressure-gradient \(\|\nabla p\|_2\) | \(\le\|B\|_2\asymp N^{3/2}\) | \(\tau^{-3/4}\) | \(\tau^{-3\gamma/2}\) |
| physical time remaining | \(N^{-2}\) | \(\tau\) | \(\tau\) |
| Fourier bandwidth | \(N\) | \(\tau^{-1/2}\) | \(\tau^{-\gamma}\) |
| active mode count | \((\eta N)^3\) | \(\tau^{-3/2}\) | \(\tau^{-3\gamma}\) |
| **forcing \(\|f\|_{H^s}\)** | \(O(1)\) | \(\tau^0\) | \(\tau^0\) |
| **\(\|\partial_t^qf\|_{H^s}\)** | \(O((Q\Omega)^q)\) | \(\tau^0\) | \(\tau^0\) |
| **absorbed low return \(\|P_{\le K_0}B\|_2\)** | \(\le\sqrt{8\pi/5}\,2c_EK_0^{5/2}N^{-1}\) | \(\tau^{1/2}\) | \(\tau^{\gamma}\) |
| **sea tracking error \(\|e\|_2\)** | \(O(N^{-1}/(\nu-\Lambda_*))\) | \(\tau^{1/2}\) | \(\tau^{\gamma}\) |
| **Type-I indicator \(\sqrt\tau\|u\|_\infty\)** | — | \(O(1)\) | \(\tau^{1/2-\gamma}\to\infty\) |

\(\int_0^TN\,dt<\infty\) (finite dissipation), \(\int_0^TN^2dt=\infty\) (BKM divergence): both hold on \(\gamma\in[1/2,1)\).

---

## D. Closed feedback loop

\[
\underbrace{f\ (\text{open-loop, }P_{>K_0}f=0)}_{(9.4)}
\;\xrightarrow[\ \text{(9.9): }\|S-S^*\|\le\|d\|/(\nu-\Lambda_*)\ ]{}\;
\underbrace{S\approx S^*(t)}_{\text{sea}}
\;\xrightarrow[\ \nabla S^*\ \text{strain + phase}\ ]{}\;
\underbrace{\Psi(\xi,s)\ \text{seed at }N_0}_{(9.10)}
\]
\[
\Psi\;\xrightarrow[\ \partial_s\Psi=a(2\Psi+\xi\!\cdot\!\nabla_\xi\Psi)-\nu|\xi|^2\Psi-\mathcal Q(\Psi,\Psi)\ (\text{S1})\ ]{}\;
\Psi(\cdot,s+S),\ S=\log2/a
\;\xrightarrow[\ \dot N=aN^3\ ]{}\;
N_{j+1}=2N_j
\]
\[
N_j\;\xrightarrow[\ (9.3)\ ]{}\;\|P_{\le K_0}B(u_{N_j},u_{N_j})\|_2\le C K_0^{5/2}c_E/N_j
\;\xrightarrow[\ \text{summable, }\sum_j\propto N_0^{-1}\ ]{}\;\|d\|\ \text{small}\;\to\;\text{(back to (9.9))}
\]
\[
\{N_j\}_{j\le J(t)}\;\xrightarrow[\ \text{S4 per-octave }\|u_{\rm band}\|_3^3\ge(q/C)c_E^{3/2}\ ]{}\;
\|u(t)\|_3^3\gtrsim J(t)\asymp\log\tfrac1\tau\;\to\;\infty .
\]

Every arrow is a displayed formula. The loop is *stabilizing on the sea* (contraction rate \(\nu-\Lambda_*>0\)) and *amplifying on the front* (gain \(a>0\)); the two are coupled only through the \(O(N^{-1})\) return, which is what makes the decomposition closeable.

---

## E. Obstruction audit — exact collision points

1. **Energy bound (F-N1).** \(\sup_t\|u\|_2\le\|u_0\|_2+\int_0^T\|f\|_2<\infty\), finite by (9.5) and \(T<\infty\). *Collision*: none — we never use energy divergence as the signature. Conceded.
2. **Finite dissipation (F-N2).** \(\nu\int_0^T\!\|\nabla u\|_2^2\asymp\nu c_E\!\int_0^T\!\tau^{-\gamma}d\tau<\infty\) for \(\gamma<1\). Conceded and consistent.
3. **ESS \(L^\infty_tL^3_x\).** *Collision point*: requires \(\limsup\|u\|_3=\infty\). Supplied by the S4 per-octave lemma + wake persistence: \(\|u\|_3^3\gtrsim\log(1/\tau)\). This is the mechanism's sole singularity signature.
4. **Fixed-finite-bandwidth no-go (F-6/F-α1, VR-L-011).** *Collision point*: hypothesis "\(u(t)\in V\) for all \(t\)". Here \(u=S+\Phi\) with \(S\in V_{K_0}\) **but** \(\Phi\not\equiv0\) and \(\operatorname{supp}\hat\Phi\) has radius \(N(t)\to\infty\). The no-go applies *verbatim to the sea alone* and we **want** it to: it certifies the sea can never blow up. The no-go says nothing about \(f\in V_{K_0}\); its hypothesis is on \(u\), and it only needs \(f\in L^1((0,T);L^2)\), which (9.5) supplies. No exclusion.
5. **Pure-swirl \(L^3\) no-go (VR-L-016).** *Collision point*: requires \(u=u^\theta(r,z)e_\theta\), whence \(\partial_\theta p_0\equiv0\) and \(P\equiv0\). Our carrier boxes sit at \(N p,Nq\) with \(p=(1,1,0)\), \(q=(1,0,1)\) — no axial symmetry, \(u\cdot\nabla p\not\equiv0\). Not in the class.
6. **One-scale self-similar no-go (NRS 1996 / Tsai 1998).** *Collision point*: NRS assumes a **time-independent** profile \(U\in L^3(\mathbb R^3)\). Our front profile is \(s\)-periodic with period \(S=\log2/a\) (discrete self-similarity), and item 3 forces \(\|U(\cdot,s)\|_3^3\gtrsim s\to\infty\) — the \(L^3\) hypothesis fails *quantitatively*, not by fiat. Tsai's finite-local-energy version assumes exact backward self-similarity (\(a\) constant, \(U\) \(s\)-independent) and no forcing; both fail here. Seregin's DSS Liouville results are preprint-status design warnings and assume unforced \(\mathbb R^3\).
7. **Galerkin global existence.** Every finite cutoff of this system exists globally; the pilot in F therefore measures the *renormalized stage gain* \(\mathcal A_j\), never "blow-up". Conceded explicitly.
8. **Smooth-forcing high-frequency decay (F-N4).** *Collision point*: \(\hat f_k\) decays faster than any polynomial, so direct high-shell injection is invisible. We impose the **strictly stronger** \(\hat f_k\equiv0\) for \(|k|>K_0\). F-N4 is **conceded, not evaded**: the entire cascade above \(K_0\) is autonomous (eq. (9.7), second line, has no \(f\)). The force's legal jobs are exactly three — (a) the finitely many stages \(j\le J_{\rm ctrl}\) of (9.10), (b) basin entry from \(u_0=0\) (B.5), (c) maintenance of the sea's phase/strain pattern via (9.9).
9. **Mesoscopic \(\gamma<1\) empty-child no-go, \(D_N\le2\kappa^2\tau^2c_EM^{\rm eff}_N/N^3\).** *Collision point*: the front's relative width exponent. We take \(W_N=\lfloor\eta N\rfloor\), \(\eta<1/3\), i.e. \(\gamma=1\) exactly, giving \(D_N\le2\kappa^2\tau^2c_E\eta^3c_M\) — a **constant**, not \(o(1)\). The no-go does not exclude us; it also gives no help. Requiring \(D_N\ge1/2\) forces \(c_E\ge[4\kappa^2\tau^2\eta^3c_M]^{-1}\); at \(\eta=3/16,\kappa=1+\eta,\tau=1,c_M=1\) this is \(c_E\gtrsim27\), legal by the \(c_E\)-collapse (S3) since total front energy \(=2c_E/N_0\) is made small by starting deep.
10. **Diagonal cross-talk gate.** *Collision point*: two same-scale relays produce diagonal parent pairs on the target child shell (\(222/2483\) intended-to-cross power). **This lens does NOT fix it.** The controller absorbs only the fraction of cross-talk landing in \(|k|\le K_0\), which by (9.3) is \(O(K_0^{5/2}/N)\) — negligible. The gate survives intact and must be cleared by scale-stagger (S6) or a wider alphabet.
11. **CSTY 2009 axisymmetric Type-I exclusion.** *Collision point*: the \(\gamma=1/2\) branch gives \(\sqrt\tau\|u\|_\infty=O(1)\) — squarely Type-I. CSTY's hypotheses are (i) axisymmetric, (ii) \(\mathbb R^3\), (iii) unforced; we violate all three. Nevertheless this is the **most dangerous** row, and the design response is structural: take \(a(s)\) increasing so \(\gamma\in(1/2,1)\), giving \(\sqrt\tau\|u\|_\infty\asymp\tau^{1/2-\gamma}\to\infty\), strict Type II, while dissipation stays integrable and BKM stays divergent.
12. **KNSS ancient-solution Liouville.** *Collision point*: KNSS(b) kills bounded ancient axisymmetric-with-swirl solutions obeying \(|u|\le C/r\). Our wake has \(|u_j|\asymp N_j\) at distance \(r\asymp N_j^{-1}\), i.e. **exactly** \(|u|\asymp C/r\) — the bound is saturated. Escape is *only* via non-axisymmetry: KNSS(b) requires axisymmetry, which our \((p,q)\) carrier geometry breaks. **Explicit warning: any axisymmetric reformulation of this candidate is killed by KNSS(b).**
13. **Front-resolution threat model (TM-22, ≥7 points per front scale).** In the RG frame \(\xi=k/N(t)\), \(ds=N^2dt\), the front profile \(\Psi(\xi,s)\) does **not** sharpen: a fixed \(\xi\)-grid resolves it uniformly in \(j\). TM-22 is evaded by construction, not by refinement. TM-20 (Heun+central-difference advection amplification) is avoided by using the repo's exact zero-padded spectral convolution.

---

## F. Minimal falsification experiment (≤1 h)

**T1 (exact, rational).** Verify identity (9.1) coefficient-by-coefficient on the exact 7-mode carrier field of `exact_leray_relay.build_exact_relay_triad()`, cross-checked against `exact_carrier_record_verifier.verify_serialized_strict_orientation_records`. *Must be exact rational arithmetic.* **Kill: any nonzero residual.** (Also check \(\widehat B(0)=0\) exactly.)

**T2 (float, the decisive measurement).** Using `mesoscopic_local_fft.measure_local_fft_cloud` (zero-padded, \(K=4W-3\)) with fixed-relative width \(\eta=3/16\), \(N\in\{16,32,48,64\}\), \(K_0\in\{2,4,8\}\), \(c_E=1\): measure \(R(N,K_0):=\|P_{\le K_0}B(u_N,u_N)\|_2\) and fit \(R\asymp N^{-p}K_0^{\,r}\).
*Predictions*: \(p=1\) unconditionally (9.3); \(p=5/2\) if the low sliver is phase-incoherent; \(r=5/2\).
**Success**: \(p\ge1\) with \(\ge4\) points, \(R^2\ge0.98\), and \(R\le\sqrt{8\pi/5}\,2c_EK_0^{5/2}/N\) at every row (the bound must hold — it is a theorem, so a violation means a code/convention error).
**Kill**: \(p<1-0.1\) at any \(K_0\) (falsifies (9.3), hence the whole lens), or \(R\) non-decreasing in \(N\).

**T3 (float, interface).** `mesoscopic_galerkin.run_small_mesoscopic_galerkin` with \(N=4\), cutoff \(3N\), \(64^3\) grid, RK4×16, plus the added open-loop force (9.4) built from a two-harmonic \(S^*\) with \(\Lambda_*=\nu/2\). Log \(\|e(t)\|_2\) and \(\|d(t)\|_2\).
**Success**: \(\|e(t)\|_2\le\|d\|_\infty/(\nu-\Lambda_*)\) at every logged step (validates (9.9)) and \(\|e\|\) non-increasing after one relaxation time \(1/(\nu-\Lambda_*)\).
**Kill**: \(\|e\|\) grows monotonically at fixed \(N\) (means \(\Lambda_*\ge\nu\); then the design condition, not the lens, is refuted — retry once at \(\Lambda_*=\nu/4\), then kill).

**T4 (float, PO-09).** From \(u_0\equiv0\), \(K_0=8\), drive the truncated \(V_{K_0}\) system along a prescribed cubic-spline path to a target seed at \(N_0=3\), \(\eta=3/16\) using (9.4); integrate the *untruncated* system at cutoff \(3K_0\).
**Success**: relative tracking error \(\|P_{\le K_0}u(t_1)-S^*(t_1)\|_2/\|S^*(t_1)\|_2\le0.1\).
**Kill**: \(\ge0.5\).

**Arithmetic split**: T1 exact rational (mandatory); T2–T4 binary64 for discovery. Promotion beyond "numerical candidate" requires interval enclosure of only two scalars — \(\Lambda_*=\|\nabla S^*\|_\infty\) and the constant \(\sqrt{8\pi/5}K_0^{5/2}\) in (9.3) — both of which are finite trigonometric-polynomial quantities, hence interval-enclosable with existing machinery (VR-C-006/007). **Reused modules**: `mesoscopic_local_fft`, `mesoscopic_galerkin`, `mesoscopic_cloud_scaling`, `leray_response_relay.leray_advection/leray_project`, `exact_leray_relay`, `exact_carrier_record_verifier`. **Not reused**: `leray_response_relay.relay_stage` (rejected Fejér-orbit route).

---

## G. Proof chain to Clay (D)

1. **(9.1)** exact low-return identity, and \(\widehat B(0)=0\). *[finite Fourier; Lean-able now]*
2. **(9.3)** unconditional bound \(\|P_{\le K_0}B(u_N,u_N)\|_2\le\sqrt{8\pi/5}\,2c_EK_0^{5/2}/N\), and its geometric summability over \(N_j=2^jN_0\). *[Cauchy–Schwarz; Lean-able]*
3. **(9.5)** the controller (9.4) is Clay-(D) admissible: \(P_{>K_0}f=0\), \(\sup_t\|\partial_t^qf\|_{H^s}<\infty\), \((1+|t|)^{-K}\) decay. *[finite-dim norm equivalence + \(T<\infty\)]*
4. **(9.9)** sea contraction: \(\Lambda_*<\nu\Rightarrow\limsup\|S-S^*\|_2\le\sup\|d\|/(\nu-\Lambda_*)=O(N^{-1})\). *[exact energy neutrality F-12 + torus spectral gap]*
5. **(9.10)+B.5** basin entry from \(u_0\equiv0\) to the seed at \(N_0\) (PO-09). *[full actuation on \(V_{K_0}\) + step 4]*
6. Existence of a fixed-relative (\(\gamma=1\), \(\eta<1/3\)) \(s\)-periodic profile \(\Psi\) of the S1 front flow with net forward efficiency \(\chi_{\rm shape}>0\). *[open — the sole remaining hard analytic step]*
7. Diagonal cross-talk of the chosen carrier graph is geometrically excluded (scale-stagger S6) or bounded below the child energy. *[open]*
8. Interval closure of the per-stage flux budget with the \(|k|\le K_0\) line item **removed** by step 2. *[CAP]*
9. \(\|u(t)\|_3^3\gtrsim\log(1/\tau)\) via S4 per-octave + wake persistence; \(T-t=\int N^{-2}ds<\infty\). *[analysis]*
10. Galerkin limit ⟹ classical solution on every \(t<T\); ESS endpoint + step 9 ⟹ non-extendability at \(T\) ⟹ **Clay (D)**.

## H. Effect on the cloud candidate's 10 unproved lemmas (§13)

| # | effect | reason |
|---|---|---|
| 7 (low-band controller satisfies (10.1)–(10.2)) | **ELIMINATED** | (9.4)+(9.5), with V1 rejected by (9.6) — the controller must be open-loop |
| 6 (future seeds form one smooth datum) | **MOSTLY ELIMINATED** | \(u_0\equiv0\) + open-loop force reaches the basin (B.5); reduced to "seeds beyond \(J_{\rm ctrl}\) come from the cascade" |
| 1 (clouds with \(M_N\!\asymp\!N^3\) satisfying (5.4)–(5.5)) | **WEAKENED** | the \(|k|\le K_0\) part of (5.4) becomes theorem (9.3) — no phase-code design constraint at low \(k\); the \(K_0<|k|\le2\eta N\) part survives untouched |
| 4 (interval budget with every output band) | **WEAKENED** | the \(|k|\le K_0\) output band leaves the budget with closed-form total \(\le2\sqrt{8\pi/5}\,2c_EK_0^{5/2}/N_0\) |
| 9 (local energy + pressure in the limit) | **HALF-FREE** | the energy half is automatic (F-N1 + step 3); the pressure half is untouched |
| 2, 3, 5, 8, 10 | **UNCHANGED** | all live in the autonomous front, which the force provably cannot reach (F-N4, conceded) |

**Net: 1 lemma eliminated, 1 mostly eliminated, 2 weakened, 1 halved, 5 untouched.**
