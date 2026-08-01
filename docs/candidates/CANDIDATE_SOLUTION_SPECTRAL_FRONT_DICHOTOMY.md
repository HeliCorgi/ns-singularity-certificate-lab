# CANDIDATE SOLUTION — spectral front monotone and bandwidth–dissipation dichotomy

**Status: unverified solution candidate (regularity side); not a proof**

**Label: PROOF CANDIDATE** (core identities exactly verified; Clay bridge open)

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
- Saturation deficits \(\mathfrak d=(\text{gap})/(\text{bound})\):
  0.865 (\(s=1\)), 0.740 (\(s=2\)), 1 (four-parent, trivially).
  \(\mathfrak d\) decreasing in scale is the open K3 question.
- Delocalization sweep: \(\mu_N=M^{\rm eff}/N^3\) has a positive
  scale-independent floor on the fixed-relative family
  (\(0.021/0.074/0.176/0.419/0.704\) at \(\eta=0.10\ldots0.30\), \(N=64\)) —
  the K4 kill did **not** fire.

## 5. Unproved lemmas / obligations (from §G of the derivation)

1. O-3: \(H_0(T_{\max})>0\) for torus blow-up (backward-uniqueness type).
2. O-4: finite-trig → \(C([0,T);H^m)\) extension of (I.1)–(I.3).
3. O-5/O-6: monotone (M) and dichotomy D on the true solution class.
4. O-7: Besov effective bandwidth replacing hard band-limitation in Lemma K.
5. O-8: uniform positive deficit \(\mathfrak d_0>0\) (or its failure — K3).
6. O-9: incompatibility of double saturation at \(\gamma=1/2\) (hard step).
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
