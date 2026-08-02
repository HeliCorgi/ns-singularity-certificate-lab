# CANDIDATE SOLUTION — spectral front monotone and bandwidth–dissipation dichotomy

**Status: unverified solution candidate (regularity side); not a proof**

**Label: PROOF CANDIDATE for the identities and the dichotomy; the
quantitative sharpening (uniform deficit floor, old O-8) is REFUTED —
see the Verification Sprint V1 correction block below.**

> **Verification Sprint V1 (2026-08-02) correction.** Workstream C
> ([defect decomposition](../research_notes/verification_sprint_v1/lambda_O9_defect_decomposition.md))
> verified (I.1)–(I.4) and Lemma K exactly and established:
> (1) the deficit reported by `front_gap_identity` is minimised over
> viscosity at \(\nu_*=\mathrm{Cov}/(2V_r)\) with the **viscosity-free,
> scale-invariant value** \(\mathfrak d_*=1-\mathrm{Cov}^2H_r/(V_rG_r)\) —
> the pilot's 0.865/0.740 were \(\nu=1/40\) artefacts; the honest exact
> values for the relay triad are \(\mathfrak d_*^{\rm full}=3/4\),
> \(\mathfrak d_*^{\rm in}=6/11\);
> (2) **O-8 is refuted**: the exact relay family \(B=C=1\),
> \(D\to0\) has \(\mathfrak d_*^{\rm full}=15D^2/(2+18D^2)\to0\)
> (certified exact points down to \(15/524306\)), so no band-uniform
> positive deficit floor exists;
> (3) the energy-neutrality route to O-9 is **vacuous at r = 0**
> (saturation \(a_k=\gamma(x_k-\mu)e_k\) satisfies \(T_0=0\)
> automatically), and off-support leakage can be driven to \(4\times10^{-10}\);
> (4) the surviving quantitative structure is **O-8′**: near-saturation
> forces single-shell collapse, \(\mathfrak d_*\) vanishing proportionally
> to the relative spectral spread \(V_0/\mu^2\) with order-one measured
> ratios (relay family limits: full \(\to5/4\), in-support \(\to1/2\);
> smallest observed constrained ratio 0.126, from an under-converged
> optimiser, recorded as unresolved);
> (5) the \(\varepsilon\)-regularised identity for
> \(\log(N_0^2+\varepsilon)\) holds with degradation factor exactly
> \(1+\varepsilon H_0/H_1\), and \(N_0^2\ge1\) on \(\mathbb T^3\) makes
> small \(\varepsilon\) harmless unconditionally;
> (6) a ledger-convention correction: the module's \(G_r\) is the
> in-support moment; the published \(G_0=\tfrac12\|\mathbb P(u\cdot\nabla
> u)\|_2^2\) is available as `full_nonlinear_power` and the difference is
> the leakage channel.
> **Consequence: this candidate is demoted from
> "identity + dichotomy + quantitative sharpening" to
> "identity + dichotomy", with O-8′ as the only surviving quantitative
> conjecture.** The monotone (M), the dichotomy D, and Lemma K are
> untouched by the refutation.
>
> **Osgood gate (2026-08-02,
> [spectral_front_osgood_gate.md](../research_notes/verification_sprint_v1/spectral_front_osgood_gate.md)): KILL.**
> The route from this dichotomy to unconditional Clay (B) via an Osgood
> bridge \(K\le\Phi(z)+R/D\) is closed: the static maximiser satisfies
> \(K\asymp0.4\,e^{z}\) (exact-anchored family), defeating every
> Osgood-admissible \(\Phi\), and the exact identity
> \(\int KD=\int\|\partial_tu\|^2/D+\nu^2\!\int N_1^2+\nu\log(D/D_0)\)
> reduces every dynamic-depletion escape to the supercritical
> \(\dot H^1\)-bandwidth action. **Λ stands as a diagnostic criterion
> only** (the dichotomy itself remains valid and strictly stronger than
> the classical critical actions).

Ranked #1 by all three independent judge passes of the 2026-08-01 discovery
portfolio
([full derivation](../research_notes/ideas_2026_08_01/idea_monotone_dichotomy.md),
[portfolio synthesis](../research_notes/discovery_2026_08_01_front_flow_portfolio.md)).

## 1. Clay target

**(B)** — global regularity on \(\mathbb T^3\), no forcing, smooth mean-zero
divergence-free datum, fixed \(\nu>0\). (Not (A): the auxiliary identity
\(\int_0^TN_0^2\,dt=\frac1{2\nu}\log\frac{H_0(0)}{H_0(T)}\) requires
\(H_0(T)>0\), which fails for whole-space self-similar profiles.)

## 2. Complete candidate quantities

Repo conventions. Modal ledger \(e_k,a_k\) as in `modal_front_actions`, plus
the **nonlinear power** \(n_k\) (paired modal energy of
\(\mathcal N=-\mathbb P(u\cdot\nabla u)\)); moments
\(H_r=\sum x_k^re_k\), \(T_r=\sum x_k^ra_k\), \(G_r=\sum x_k^rn_k\),
\(x_k=|k|^2\); bandwidth \(N_r^2=H_{r+1}/H_r\).

**Identity I** (exact): with \(p_r(k)=x_k^re_k/H_r\), \(\mu=N_r^2\),
\(V_r=\mathrm{Var}_{p_r}(x)\),
\[
\tfrac12\tfrac{d}{dt}\log N_r^2
=\frac{\mathrm{Cov}_{p_r}(x,g)}{\mu}-\nu\frac{V_r}{\mu},
\qquad g_k=a_k/e_k .
\tag{I.1}
\]
**Closable bound**: modal Cauchy–Schwarz \(|a_k|\le\sqrt{e_kn_k}\) *before*
dividing by \(e_k\) gives \(|\mathrm{Cov}|\le\sqrt{V_rG_r/H_r}\) (I.2), and a
square completion in \(\sqrt{V_r}\):
\[
\boxed{\ \frac{d}{dt}\log N_r^2\ \le\ \frac{G_r}{2\nu H_{r+1}}\ }
\tag{I.3}
\]
with the exact nonnegative gap decomposition (I.4) whose total is rational.

**Front wavenumber and monotone**: \(K=G_0/(2H_1^2)\),
\(D=2H_1\):
\[
\boxed{\ \Lambda(t)=\log N_0^2(t)-\frac1{2\nu}\int_0^tK\,D\,ds
\ \text{ is non-increasing.}}
\tag{M}
\]
**Lemma K** (unconditional): band-limited to \(N\) ⇒
\(K\le S_N=\sum_{0<|k|^2\le N^2}|k|^{-2}\ (=4\pi N+O(1))\), no phase or
coherence assumption.

**Dichotomy D** (target theorem): either
\(\int_0^\infty K D\,dt<\infty\) and the solution is regular (via
\(N_0\) bounded ⇒ \(u\in L^\infty_tH^1\subset L^\infty_tL^6\), Serrin
\((6,\infty)\)), or \(T_{\max}<\infty\) with
\(\int_0^{T_{\max}}KD\,dt=\infty\). Since
\(KD\le\|u\|_\infty^2\) and \(KD\le C_S^2\|\nabla u\|_{L^3}^2\), the action
\(\int KD\) is **dominated by both the Serrin \((\infty,2)\) and the
critical vorticity \((3,2)\) actions**: the criterion is strictly stronger
than either classical one.

## 3. Why this can feed a Clay statement

With Lemma K, blow-up forces \(\int N D\,dt=\infty\); with finite
dissipation \(\int D<\infty\) this re-derives the bandwidth window
\(\gamma\in[\tfrac12,1)\) with **no BKM, no phase assumption, no flux
sign** — and \(\gamma=\tfrac12\) (the Type-I corner) demands *simultaneous*
saturation of the Cauchy–Schwarz gap, the square-completion gap, and
Bernstein. Obligation O-9 — proving these saturations incompatible for
divergence-free fields — is the single hard step to (B); the repo's
cross-talk and phase-coherence negative results become ingredients there.

## 4. Exact numerical support (2026-08-01 pilot)

`spectral_front_monotone.py` + `outputs/spectral_front_ledger_v1/`, all
load-bearing checks in exact `Fraction` arithmetic:

- (I.3), the rational (I.4) total, per-mode \(a_k^2\le e_kn_k\), and
  \(K\le S_N\) hold exactly on the relay triads (\(s=1,2\)) and the
  four-parent expanded-carrier field.
- Saturation deficits: the originally reported 0.865/0.740 were
  \(\nu=1/40\) artefacts (see the correction block); the viscosity-free
  exact values are \(\mathfrak d_*^{\rm full}=3/4\),
  \(\mathfrak d_*^{\rm in}=6/11\) for the relay triad, and the K3 question
  is now **answered negatively for a uniform floor** by the exact
  \(D\to0\) relay family; the open question is O-8′
  (\(\mathfrak d_*\gtrsim c\,V_0/\mu^2\)).
- Delocalization sweep: \(\mu_N=M^{\rm eff}/N^3\) has a positive
  scale-independent floor on the fixed-relative family
  (\(0.021/0.074/0.176/0.419/0.704\) at \(\eta=0.10\ldots0.30\), \(N=64\)) —
  the K4 kill did **not** fire.

## 5. Unproved lemmas / obligations (from §G of the derivation)

1. O-3: \(H_0(T_{\max})>0\) for torus blow-up (backward-uniqueness type).
2. O-4: finite-trig → \(C([0,T);H^m)\) extension of (I.1)–(I.3).
3. O-5/O-6: monotone (M) and dichotomy D on the true solution class.
4. O-7: Besov effective bandwidth replacing hard band-limitation in Lemma K.
5. ~~O-8: uniform positive deficit~~ **REFUTED** (sprint V1); replaced by
   **O-8′**: \(\mathfrak d_*\ge c\,V_0/\mu^2\) with band-uniform \(c>0\)
   (measured ratios 0.126–5/4; undecided — needs the equality-constrained
   optimiser rerun named in the sprint note).
6. O-9: incompatibility of double saturation at \(\gamma=1/2\) (hard step;
   the energy-neutrality and leakage routes are now known **vacuous** at
   \(r=0\) — the surviving route is the O-8′ single-shell-collapse
   structure: saturating fields must spectrally collapse, and a cascade
   front needs \(V_0/\mu^2\) bounded below).
7. Lean: (I.1)–(I.4) + Lemma K as finite-sum lemmas (no new axioms).

## 6. Falsifiability

- K1: any exact field with \(\Gamma^{\rm CS}<0\)/failed (I.4) — *fired never
  so far*; would abandon the lane.
- K2: negative Lemma-K margin on a band-limited field.
- K3: \(\mathfrak d\to0\) along genuine \(N\)-sweeps of evolving fields —
  demotes (M) to a consistency diagnostic (identity survives, sharpening
  dies).
- Evolution sign check: any resolved Galerkin trajectory violating (I.3)
  falsifies the derivation outright.

## 7. Next minimal experiment

Extend the ledger to *evolving* fields: log
\(\mathfrak d(t)\), \(K(t)\), \(S_{N(t)}\) along the registered small
Galerkin runs (`mesoscopic_galerkin`, `carrier_two_stage_galerkin`) and a
band-doubling family, testing K3 and the (I.3) sign on trajectories (float
lane, sign-level only, per TM-22).

## 8. Non-claims

No global regularity, no Clay statement. (M) is an inequality whose action
finiteness is a hypothesis; the window re-derivation is conditional on the
log-critical shell class. Nothing here contradicts the singularity-side
lanes: if the cloud lane ever closes, alternative (II) of the dichotomy is
the branch it must inhabit, and \(\Lambda\) then *quantifies* its required
front efficiency.
