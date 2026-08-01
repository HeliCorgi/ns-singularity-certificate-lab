# Discovery portfolio — renormalization front flow, 14 directions, two pilots (2026-08-01, second lane)

**Label: FORMAL ANSATZ / SYMBOLIC CANDIDATE / AUDIT REQUIRED** (per-item labels below).
This note records a discovery-phase portfolio: no item is a PDE theorem, a
singularity, or a Clay result. Full derivations live in
[`ideas_2026_08_01/`](ideas_2026_08_01/); this note is the corrected synthesis.

## 1. What this session produced

1. A **carrier-frame front flow** — the continuum home of the repository's
   doubling response operator \(\mathfrak T_2\) — with an adversarial audit
   that *confirmed the algebra and killed the steady case* (§2).
2. **Fourteen worked research directions** (13 lens documents + the seed),
   each with Clay target, derived central equations, scaling table, feedback
   loop, obstruction audit, minimal falsification experiment, and proof
   chain; failed sub-variants are retained inside each document with their
   exact failing equations (§3).
3. Three independent judge reports and a **top-3 selection** (§4).
4. **Two executed pilots**: the renormalized Galerkin cascade (forward
   integration of the doubling map — negative in the scanned box, §5) and
   the exact spectral front-gap ledger (identity certificates pass; the
   fixed-relative cloud family *survives* its new delocalization kill test,
   §6).

## 2. The front flow, as corrected by adversarial audit

Ansatz \(\widehat u(k,t)=N(t)^{-2}\Psi(k/N(t),s)\), \(ds=N^2dt\), in repo
conventions. Term-by-term substitution (audit: **confirmed**, including the
Riemann factor \(N^3\), all exponents, and the sign of \(\mathcal Q\)) gives

\[
\partial_s\Psi=a(s)\,(2\Psi+\xi\cdot\nabla_\xi\Psi)-\nu|\xi|^2\Psi
-\mathcal Q(\Psi,\Psi),\qquad
a(s)=\frac{\dot N}{N^3}=-\tfrac12\frac{d}{dt}N^{-2},
\]
\[
\mathcal Q(\Psi,\Psi)(\xi)=i\,P_\xi\!\int\bigl((\xi-\lambda)\cdot\Psi(\lambda)\bigr)
\Psi(\xi-\lambda)\,d^3\lambda .
\]

- \(N^{-2}(t)=N_0^{-2}-2\int_0^ta\): any orbit with \(a\ge a_->0\) has
  finite-time bandwidth divergence with \(N\asymp(2a_+(T-t))^{-1/2}\)
  (γ = 1/2, the Type-I boundary), and \(dN/ds=aN\).
- \(\mathcal Q\) restricted to carrier cells is exactly eq. (6.6) of the
  phase-coded cloud candidate; the frozen one-step map + doubling pullback is
  one explicit-Euler step of this flow.
- **Exact scaling group (audit-confirmed):**
  \((\Psi,s,a,\nu)\mapsto(\mu\Psi,s/\mu,\mu a,\mu\nu)\). Consequently all
  quadratic channel *ratios* are \(c_E\)-independent shape functionals: the
  amplitude constant is a free lever and only **shape closure** obstructs the
  lane. (Correction: this is a high-Reynolds limit statement, not a
  reduction; the frozen-response extrapolations "\(c_E=228\) suffices" are
  invalid perturbative statements and are withdrawn.)

**The decisive audit finding (kills the steady case).** A steady profile at
constant \(a\) *is literally* a backward self-similar Leray profile; the
critical wake gives \(U\asymp|x|^{-1}\), so local energy and
\(\iint|\nabla u|^2\) are finite and **Tsai (1998) forces \(u\equiv0\) on**
\(\mathbb R^3\). The lattice \(\mathbb Z^3\) admits only the discrete
dilations \(2^{\mathbb Z}\), so the surviving object is a **discretely
self-similar (log-periodic) orbit**, and the Riemann-sum continuum limit
re-imports the no-go unless the orbit's per-period oscillation of the
critical norm stays \(O(1)\).
**Pre-registered kill gate for every future orbit search: an orbit whose
per-period critical-norm oscillation vanishes under refinement is dead by
Tsai, and convergence of any iteration to a *steady* shape is a failure
signature, not a success.**

**Corrected 1-D reduction.** The shell reduction
\(\partial_se=a(2e+\xi e')-2\nu\xi^2e-\partial_\xi F\) is exact
(audit-confirmed, including the vanishing of the dilation term on the
critical spectrum \(e_c=2c_Ec\xi^{-2}\) and the wake flux depletion
\(F'=-4\nu c_E\)). The *closure* \(F=\chi\xi e^{3/2}\) together with
\(e\equiv e_c\) is **over-determined and inconsistent on \(1<\xi<\varphi\)**
(golden ratio); the previously displayed profile and the threshold
\(\chi\sqrt{c_E}\ge\sqrt2\,\nu\) survive only as a necessary-condition
skeleton, and the correct object is the nonlinear eigenvalue problem for
\((e,a)\) jointly — open. See
[`ideas_2026_08_01/audit_front_flow_seed.md`](ideas_2026_08_01/audit_front_flow_seed.md)
for the full 12-item erratum (front energy \((2N)^{-1}\!\int e\); period
constraint \(a_+\tau_{\rm stage}=\log2\); translation-split suppression
\((W|\Delta x|)^{-3}\) in amplitude with floor \(W^{-3}\); the helicity
budget no-go REFUTED — slack \(2^{J/2}\) vs \(J\) — and retained only as the
design rule \(|h|\) bounded away from 1; multi-type cycles need
\(\prod\lambda_j=2\) per super-period with irrational per-stage ratios for
equal-ratio cycles).

**Flux→\(L^3\) per-octave lemma (PROOF CANDIDATE, repaired).** With
Littlewood–Paley square-function bookkeeping and
\((\sum a_j)^{3/2}\ge\sum a_j^{3/2}\), sustained critical transfer forces a
per-octave lower bound \(\|u_{\text{band}}\|_3^3\gtrsim q\,c_E^{3/2}\) with
\(q=2^{3/2}\chi\sqrt{c_E}\)-normalization, and heat flow (a positive real
Fourier multiplier) preserves wake phase coherence with scale-free retention
\(e^{-\nu/(2a_+)}\); hence sustained cascade transfer ⇒ logarithmic
\(L^3\) divergence. This converts proof-chain steps 8–9 of the cloud
candidate into consequences of the flux budget (step 5). Obligation: the
octave-summation lemma is stated, not yet proved.

## 3. The fourteen directions (details in `ideas_2026_08_01/`)

| # | direction | status | one-line content | fate after judging |
|---|---|---|---|---|
| 1 | aniso-multiscale | FORMAL ANSATZ | component-locked anisotropic collapse \(u_j=CL_jU_j\), anisotropic Leray projector, exponent polytope \(a_1\in(2/5,1/2]\), ceiling \(a_1-a_3<1/6\) | only real TARGET-U attack; near-fatal: \(M^{\rm eff}/N^3\propto\tau^{3/16}\to0\) vs L-11; Prandtl-ill-posed core |
| 2 | multi-center (NSDL) | FORMAL ANSATZ | exact skew-dyad mutual strain \(\sigma=-(2\Gamma/\pi d^2)\sin^2\theta\cos^3\theta\sin2\varphi/(\dots)^2\), max at \(\theta^*\!\approx\!47.3^\circ\); lacunary cross-talk gate | mechanism fatal as derived: \(\hat s\cdot e_1\equiv0\Rightarrow\dot d=0\) at line-vortex order; strain formula + gate survive as tools |
| 3 | space-fourier-mismatch | FORMAL ANSATZ | spectrally fat cell, pinning law \(NL=(E_0/\nu^2L)^{1/3}\), unique corner \(\mu=\varepsilon=2/5\) | corner arithmetic checks; fatal: pinning slaves \(N\) to \(L\), and \(c_E(\tau)\to\infty\) is the inviscid problem relabelled |
| 4 | pressure-feedback | FORMAL ANSATZ | pressure Hessian = depolarization tensor; restricted Euler ≡ spherical choice; eigenframe-locked family exists numerically | fatal to realization: companion sits inside the core; \(p_2^{\rm mat}>1/2\) ⇒ viscously arrested |
| 5 | lagrangian-degeneracy | FORMAL ANSATZ + 2 PROOF-CANDIDATE lemmas | \(\dot\sigma_i=\sigma_i(r_i\!\cdot\!Sr_i)\), Riccati \(\dot\lambda_c=\Theta\lambda_c^2\), \(\gamma=1/\Theta\) | γ-window self-contradictory (uses Cauchy formula its own L3 forbids); L1 (\(\sigma_3\to0\Rightarrow\)BKM) and L3 survive, Lean-able |
| 6 | parametric-resonance | SYMBOLIC CANDIDATE (no-go) + SPECULATION (escape) | **Prop 1: every tree/forest relay has Floquet exponent \(\le-\min\nu|k|^2\) for *all* modulations** — retro-explains every registered relay failure; cyclic holonomy \(\mathcal C\neq-1\) escape | Prop 1 is the portfolio's best new structural theorem; the exhibited gadget is fixed-cardinality (already-rejected class); escape needs \(c_E\gtrsim\nu^2N^3\) unless K2 fails |
| 7 | helicity-ledger | SYMBOLIC (design rule) | helical triad ledger; 2-colouring \(s(A_1)=s(B_2)=+\), \(s(A_2)=s(B_1)=-\) predicts **exact zero** on both diagonal cross-talk channels at \(|k|^2=2\) | T1 is Waleffe 1992 (novelty-fatal as theorem); D1 2-colouring is new, binary, and directly testable against the registered 16/16 negative |
| 8 | log-periodic-multitype | FORMAL ANSATZ | L-periodic carrier cycles alternating the two exact relays; **LP-9: odd-prime divisibility ⇒ exact sumset miss at every dyadic level**; log-4 observable signature | judges split (rigor #2, tractability #11): same-radius stage-0/1 cross-talk outside LP-9 scope; \(\eta\le1/14\) costs \(50\times\) in \(M^{\rm eff}\); exact closure needs irrational ratio — carriers can only witness G1–G3 |
| 9 | forcing-controller | FORMAL ANSATZ | exact low-band controller; open-loop necessity \(\|\partial_tf\|\gtrsim2\nu N^2\|P_{\le K_0}B\|\); decoupling theorem | clean but self-limiting: proves \(f\) cannot reach the front, so (D) reduces to (B); hygiene, not mechanism |
| 10 | critical-quantity (\(\mu_N\)) | PROOF CANDIDATE (diagnostic) | \(M^{\rm eff}\) floor re-derived from the *flux* side; interpolation \(\|u\|_{X^{-1}}\le3\mu_*(\|u\|_2\|\nabla u\|_2)^{1/2}\); frozen-parent no-go = true evolution to \(O(\tau)\) | Theorem D asymptotically vacuous (floor \(\to0\)), but the \(\mu_N\) kill test is decisive for the cloud lane — **executed, §6** |
| 11 | monotone-dichotomy (\(\Lambda\)) | **PROOF CANDIDATE** | closable front identity (I.1)–(I.4), Lemma K \(K\le S_N\), monotone \(\Lambda=\log N_0^2-\frac1{2\nu}\int KD\), dichotomy: \(\int KD<\infty\Rightarrow\) regular | **#1 across all judges; no equation-level flaw found; promoted to candidate document — executed pilot, §6** |
| 12 | two-scale-core | FORMAL ANSATZ (negative) | sea→front matching, back-reaction \(O(N^{-2})\) decoupling, LD drain law | valuable hygiene: kills "outer strain rescues the front"; saturates KNSS(b) — survival rests on non-axisymmetry |
| 13 | moving-front-ode | FORMAL ANSATZ | \(\dot N=\Phi(N,E_N,\Pi_N,\nu)\) under 3 closures; front–BKM identity \(\int\|\omega\|_\infty dt=(\theta/\chi\sqrt m)\log N\) ⇒ **BKM divergence is never numerically observable**; trade-off curve \(\chi_{\rm req}(g)\), optimum \(g^*\to2/5\) | methodologically important; accounting error at §E.8 (mixed normalizations) to repair; \(g^*=2/5\) vs γ=1/2 vs Lens-5's window is an open three-way contradiction (§7) |
| 14 | front-flow (seed) | FORMAL ANSATZ (connective tissue) | §2 above | steady case killed by Tsai; survives only as DSS-orbit frame + \(c_E\)-collapse + RG-integration method |

## 4. Selection

**Top 3 (judge consensus, all three reports):**

1. **Λ — spectral front monotone & bandwidth–dissipation dichotomy**
   (lens 11; regularity side, Clay (B) target). Only candidate whose core
   (I.1)–(I.4) + Lemma K is provable *now*, exactly rational, pressure-blind,
   Lean-able in the `MesoscopicDuhamelNoGo` style; its action
   \(\int K D\,dt\) is dominated by both the Serrin \((\infty,2)\) and the
   critical vorticity \((3,2)\) actions, so the dichotomy is strictly
   stronger than either. Promoted to
   [`../candidates/CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md`](../candidates/CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md).
2. **\(\mu_N\) — amplitude-delocalization floor** (lens 10): the second,
   flux-side derivation of the \(M^{\rm eff}\gtrsim N^3\) requirement, its
   Rényi-\(\tfrac12\) entropy reading, and the sweep test executed in §6.
3. **Parametric resonance Prop 1 + holonomy escape** (lens 6): the tree/forest
   Floquet no-go is the portfolio's strongest new *structural* theorem (it
   retro-explains every registered relay failure and is exactly checkable);
   the \(\mathcal C\neq-1\) cyclic escape is the only known way to grow a
   child with identically-zero instantaneous margin — its K2 test (does
   \(\mathcal C\to-1\) in the width-saturating limit?) is the next exact
   experiment.

**Held in reserve:** log-periodic multi-type scale-stagger (the only direct
attack on the registered diagonal cross-talk gate; sumset-miss margins are
satisfiable at \(\rho=2/3,\ \eta<1/9\) and \(\rho=3/4,\ \eta<1/12\)) and the
helicity 2-colouring D1 (binary exact test that could flip the 16/16
negative).

## 5. Pilot A — renormalized Galerkin cascade (negative in the scanned box)

`src/ns_certificate_lab/renormalized_cascade.py`,
`experiments/run_renormalized_cascade.py`,
`outputs/renormalized_cascade_v1/`. Iterates: full Galerkin evolution for
one parabolic stage \(\tau N^{-2}\) → exact lattice doubling pullback
\(w(k)=4v(2k)\) (preserves reality and exact divergence-freeness; energy
bookkeeping matches the continuum \(\lambda=2\) factor through
\(16\times\tfrac18\)) → optional sea drop → projective renormalization, with
the pre-normalization gain \(\mathcal A_j\) and lag-1/2/3 shape overlaps
recorded. A noise-floor tracker terminates once binary64 roundoff
(re-amplified by \(1/\sqrt{\mathcal A_j}\) per stage) reaches 10% of the
field.

Result over \(c_E\in\{1,30,100,228,500\}\times\{\text{front-only, keep-sea}\}\)
at \(N_0=4\), \(W=2\), \(\nu=1/40\), \(\tau=1/4\): **no positive-gain plateau
and no recurrent shape anywhere.** Nine of ten runs collapse into the float
noise floor within 6–11 stages; the strongest run (\(c_E=500\), front-only)
completes all 12 stages classified *decaying* with final gain
\(9.4\times10^{-5}\). Gains increase monotonically with \(c_E\) (final gain
\(10^{-10}\!\to\!10^{-4}\)), but nothing approaches the closure level
\(\mathcal A\ge1\).

Scope and honesty: this rejects an *attracting* positive-gain orbit of the
\(\lambda=2\) decimation map at the smallest admissible lattice realization
(54 modes) only. It does not touch the continuum \(\mathfrak T_2\) question;
the decimation pullback (audit item 11) is a lossy stand-in for the
continuum dilation; and by §2's Tsai gate even a steady plateau would not
have been a success. Next realizations that remain open: larger \(N_0/W\),
interpolating (non-decimating) pullback, multi-type \(L\ge2\) cycles, and
repulsive-orbit (shooting) searches.

## 6. Pilot B — exact spectral front-gap ledger (identities pass; cloud family survives its kill test)

`src/ns_certificate_lab/spectral_front_monotone.py`,
`experiments/run_spectral_front_ledger.py`,
`outputs/spectral_front_ledger_v1/`. All load-bearing checks are exact
`Fraction` arithmetic.

**Exact lane.** On `build_exact_relay_triad(scale=1,2)` and the four-parent
expanded-carrier field: the modal Cauchy–Schwarz certificates
\(a_k^2\le e_kn_k\), the rational gap identity (I.4-total), the closable
bound (I.3), and Lemma K \(K\le S_N=\sum_{0<|k|^2\le N^2}|k|^{-2}\) all hold
exactly (no kill condition K1/K2 fired). Measured saturation deficits
\(\mathfrak d\) at \(\nu=1/40\), \(r=0\): triad \(s=1\): **0.865**, triad
\(s=2\): **0.740** (decreasing in scale — the K3 question, whether
\(\mathfrak d\to0\) along a genuine \(N\)-sweep, is open and pre-registered);
four-parent field: \(\mathfrak d=1\) exactly (all its transfer targets empty
modes, so its own-mode ledger is trivial — a useful sanity case). Lemma-K
slack is \(\sim10^3\) on these far-from-Bernstein-saturated fields.

**Float lane (the K4 kill test for the cloud lane).** Delocalization ratio
\(\mu_N=M^{\rm eff}/N^3\) of the fixed-relative sparse parents:

| \(\eta\) | \(\mu_N\) at \(N=16\to64\) | trend |
|---|---|---|
| 0.10 | 0.0161 → 0.0213 | rising to plateau |
| 0.15 | 0.0672 → 0.0737 | rising to plateau |
| 0.20 | 0.1513 → 0.1761 | rising to plateau |
| 0.25 | 0.3847 → 0.4194 | rising to plateau |
| 0.30 | 0.3847 → 0.7035 | rising |

with \(M^{\rm eff}/M\approx0.45\)–\(0.53\) (the Fejér weighting costs a
scale-independent factor ≈ 2). **\(\mu_N\) does not vanish: the only
structurally uncondemned mesoscopic family passes the independent
delocalization floor test.** (Had \(\mu_N\to0\), the cloud lane would have
died here.)

## 7. Contradiction registry (to be resolved in the verification phase)

Three directions estimate the same viscous term at the front and disagree on
the bandwidth exponent: lens 5 selects \(\gamma=1/2\) exactly (and kills
\(\gamma>1/2\)), lens 13's trade-off optimum is \(g^*\to2/5\) (sub-Type-I,
using a different normalization branch), lens 7's T3 pins \(\gamma=1/2\)
under a sign hypothesis. At most one bookkeeping is right; the discrepancy
is order-of-magnitude, not sign, and is the cheapest theoretical dispute to
settle next.

## 8. Obligations handed to the verification phase

1. **Λ lane** (candidate doc §G): O-3 (backward uniqueness \(H_0(T)>0\) on
   \(\mathbb T^3\)), O-4 (finite-band → \(H^m\) limit), O-7 (Besov effective
   bandwidth in Lemma K), O-9 (incompatibility of double saturation) — plus
   Lean formalization of (I.1)–(I.4)+Lemma K (finite sums, Cauchy–Schwarz,
   square completion; no new axioms needed).
2. **\(\mu_N\)**: promote the sweep's positive floor to a lemma for the
   Fejér-weighted family; decide the K3 deficit question on a genuine
   \(N\)-sweep of *evolving* fields.
3. **Parametric resonance**: exact cyclic-carrier enumeration with the
   \(\mathcal C=-1\) null control; K2 (width-saturating limit of
   \(\mathcal C\)).
4. **Scale-stagger**: integer-programming sumset-miss certificate at
   \(\rho=2/3\), then a two-stage Galerkin re-run of the staggered pair.
5. **Helicity D1**: add exact helicity columns to the carrier search; test
   the 2-colouring's predicted exact zeros.
6. **Front flow**: implement an interpolating (non-decimating) pullback and
   a period-\(L\) orbit search with the Tsai oscillation gate enforced;
   repair the 1-D closure into the joint \((e,a)\) eigenvalue problem.
7. Resolve the §7 exponent contradiction.

*Binding non-claims: no singularity, no regularity theorem, no Clay
statement. Every number above is either exact rational algebra on finite
trigonometric fields or a binary64 diagnostic labelled as such.*
