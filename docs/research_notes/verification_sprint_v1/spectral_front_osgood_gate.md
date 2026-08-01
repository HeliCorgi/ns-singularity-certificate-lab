# Spectral-front Osgood gate — single-question verdict

# VERDICT: KILL

**The uniform Osgood bridge \(K\le\Phi(z)+R/D\) is refuted.** The static
maximiser of \(K\) at fixed \(z=\log N_0^2\) grows like \(e^{z}\) (exact
rational anchor + measured linear-in-\(N_0^2\) growth on an explicit
divergence-free family, §2), which defeats **every** Osgood-admissible
\(\Phi\) (any \(\Phi\) dominating \(c\,e^{z}\) on a sequence \(z\to\infty\)
has \(\int^\infty ds/\Phi<\infty\)) — not merely the first candidate
\(\Phi=C(1+z)\). The dynamic-depletion escape is closed by an **exact
identity** (§3): along every strong solution
\[
\int_0^{T'}KD\,dt=\int_0^{T'}\frac{\|\partial_tu\|_2^2}{D}\,dt
+\nu^2\!\int_0^{T'}N_1^2\,dt+\nu\log\frac{D(T')}{D(0)},
\]
whose only Osgood-compatible part is the logarithm; the remainder is the
\(\dot H^1\)-bandwidth action \(\int N_1^2=\int H_2/H_1\), a supercritical
quantity whose a-priori control is equivalent to the regularity being
sought. The (I.4)-defect routes (§4) fail structurally in both directions.
Short-time evolution of the counterexample family (§5) shows the
large-\(K\) configuration is not an instantaneous artefact at critical
balance (34.6% retention of \(K\) across the half parabolic window at the
scale tested), removing the measure-zero-spike rescue for \(R\); the
quantitative weight of the refutation rests on §2 and §3.
**Consequence: \(\Lambda\) is demoted to a diagnostic criterion.** The dichotomy \(\int KD\,dt<\infty\Rightarrow\)
regularity remains true and strictly stronger than the Serrin
\((\infty,2)\) and vorticity \((3,2)\) criteria, but it is not a route to
unconditional Clay (B) by Osgood-type closure. No new candidates, no new
obligations, no literature additions are introduced by this note.

Repository conventions throughout; \(D=\|\nabla u\|_2^2=2H_1\),
\(K=\|\mathbb P(u\cdot\nabla u)\|_2^2/\|\nabla u\|_2^4=G_0^{\rm full}/(2H_1^2)\),
\(N_r^2=H_{r+1}/H_r\), \(z=\log N_0^2\ge0\) on \(\mathbb T^3\).
Supporting computations: `experiments/run_osgood_gate.py`,
`outputs/verification_sprint_v1/osgood_gate/summary.json`
(binary64 + one exact `Fraction` anchor; labels per cell below).

---

## 1. Osgood closure theorem (the GO-branch content, proven and retained)

**Lemma (Osgood closure).** Let \(u\in C([0,T_{\max});H^m)\cap
C^1([0,T_{\max});H^{m-2})\), \(m>5/2\), be the maximal strong solution, and
write the bridge in product form (this avoids all \(D=0\) divisions):
\[
K(t)\,D(t)\ \le\ \Phi(z(t))\,D(t)+R(t)
\quad\text{a.e. }t\in[0,T_{\max}),
\tag{B}
\]
with \(\Phi>0\) nondecreasing, \(\int^\infty ds/\Phi(s)=\infty\),
\(R\ge0\), \(\int_0^{T_{\max}}R\,dt<\infty\). Then \(z\) is bounded on
\([0,T_{\max})\) and \(T_{\max}=\infty\).

*Proof.* If \(u_0=0\) the solution is zero; otherwise \(H_0,H_1>0\) on
\([0,T_{\max})\) (backward uniqueness of the heat-type equation is not
needed: \(H_0(t)\ge H_0(0)e^{-2\nu\int_0^tN_0^2}\) and \(N_0^2\) is
continuous on compacts, so \(H_0>0\); \(H_1=0\) would force \(u(t)=0\),
hence \(u\equiv0\) by uniqueness forward from \(t\)). Thus
\(z=\log(H_1/H_0)\) is \(C^1\) on \([0,T_{\max})\) and \(z\ge0\)
(\(|k|\ge1\)). The monotone (M) gives \(z'\le KD/(2\nu)\) pointwise;
with (B), \(z'\le\bigl(\Phi(z)D+R\bigr)/(2\nu)\). Let
\(\Omega(z)=\int_0^z ds/\Phi(s)\); \(\Omega\) is \(C^1\), increasing,
\(\Omega(\infty)=\infty\). Then
\[
\frac{d}{dt}\Omega(z)=\frac{z'}{\Phi(z)}
\le\frac{D}{2\nu}+\frac{R}{2\nu\,\Phi(z)}
\le\frac{D}{2\nu}+\frac{R}{2\nu\,\Phi(0)} ,
\]
using \(\Phi\) nondecreasing and \(z\ge0\). Both right-hand terms are
integrable on \([0,T_{\max})\): \(\int D\le E(0)/\nu\) by the energy
equality, \(\int R<\infty\) by hypothesis. Hence
\(\Omega(z(t))\le\Omega(z(0))+E(0)/(2\nu^2)+\int R/(2\nu\Phi(0))<\infty\),
so \(z\le z_*<\infty\) on \([0,T_{\max})\). Then
\(H_1(t)\le e^{z_*}H_0(t)\le e^{z_*}H_0(0)\), so
\(u\in L^\infty([0,T_{\max});H^1)\subset L^\infty L^6\), a subcritical
Serrin class (\(2/\infty+3/6=1/2<1\)), and the solution extends past any
finite \(T_{\max}\). ∎

The lemma is unconditional and elementary; **everything below is about
whether (B) can hold**. Note (B) with a *solution-dependent, unbounded*
budget \(\int R\) is vacuous (on any compact \([0,T']\), \(KD\) is
continuous, so (B) holds trivially with \(R:=(KD-\Phi(z)D)_+\); this
absorbs any finite excess). The only meaningful gate is the **uniform**
bridge: one \(\Phi\) and one budget bound
\(\int_0^{T_{\max}}R\,dt\le C_0(\nu,E(0),H_1(0))\) for all smooth periodic
solutions. That is what is killed.

## 2. Static obstruction: \(K_{\max}(z)\asymp e^{z}\)

**Family (coherent critical spectrum).** For a band \(1\le|k|\le N\) set
\[
\widehat u(k)=\frac{P_kv_0}{|k|^2},\qquad v_0=(1,2,3),\qquad
\widehat u(-k)=\overline{\widehat u(k)}\ \ (\text{all phases }1),
\]
i.e. a real, exactly divergence-free cosine field whose modulus realises
the Bernstein-extremal spectrum \(|\widehat u(k)|\asymp|k|^{-2}\) — the
unique spectral shape that saturates Lemma K
(\(\|u\|_\infty^2/\|\nabla u\|_2^2\asymp S_N\asymp4\pi N\)). Scaling
predictions: \(H_0=O(1)\) (convergent \(\sum|k|^{-4}\)),
\(H_1\asymp N\), \(N_0^2\asymp N\), \(z\asymp\log N\), and — provided the
Leray projection removes only an \(O(1)\) fraction —
\(\|\mathbb P(u\cdot\nabla u)\|_2^2\asymp N^3\), hence
\[
K\asymp N\asymp N_0^2=e^{z}.
\]

**Measured (float lane, exact anchor at \(N=4\)).**

| \(N\) | grid | \(K\) | \(N_0^2\) | \(K/(1+z)\) | \(K/N_0^2\) | Leray retention | label |
|---|---|---|---|---|---|---|---|
| 4 | 48 | 0.7884 | 3.048 | 0.3729 | 0.259 | 0.780 | float; **exact anchor agrees to 0.0e0** |
| 8 | 80 | 2.0372 | 6.090 | 0.7259 | 0.335 | 0.833 | float |
| 16 | 144 | 4.5741 | 12.183 | 1.3069 | 0.375 | 0.857 | float |
| 24 | 208 | 7.1113 | 18.269 | 1.8210 | 0.389 | 0.864 | float |
| 32 | 272 | 9.6478 | 24.352 | 2.3011 | 0.396 | 0.867 | float |

\(K/N_0^2\) approaches a constant \(\approx0.40\): **\(K\asymp0.4\,e^{z}\)**,
and \(K/(1+z)\) is unbounded (0.37 → 2.30, roughly doubling per octave).
Controls at the same measurement pipeline: single-shell fields
(\(K=1.95\times10^{-3}\) at \(N=8\), \(4.88\times10^{-4}\) at \(N=16\) —
decaying like \(N^{-2}\), as divergence-free depletion predicts), the ABC
Beltrami core with a 5% coherent perturbation (\(K=0.505\), bounded —
Beltrami kills \(\mathcal N\)), and the mesoscopic two-box parent at
\(N=16\) (\(K=4.71\times10^{-3}\), narrowband): **the broadband coherent
family is the maximiser class, and its \(K\) grows linearly in \(N_0^2\),
i.e. exponentially in \(z\), while every control stays bounded or decays.**
The exact anchor (`Fraction` arithmetic, independent `fourier_torus`
convolution path) at \(N=4\) gives \(K\) agreeing with the float pipeline
to relative difference \(0.0\) at the printed precision (\(K=0.788411\)),
after fixing a factor-2 convention slip caught by exactly this
cross-check.

**Consequence.** For every admissible Osgood \(\Phi\) there is a sequence
of divergence-free trigonometric fields with
\(K-\Phi(z)\ge\tfrac c2e^{z}\to\infty\). The static inequality
\(K\le C(1+\log N_0)\) — and every Osgood-compatible weakening of it — is
**false**. \([\)float growth law + one exact anchor; the family is the
repository's own wake-plus-coherent-core structure, not an exotic
configuration\(]\)

## 3. Dynamic depletion is closed off by an exact identity

Along any strong solution, \(\mathcal N=\partial_tu+\nu Au\)
(\(A=-\Delta\); this is the Navier–Stokes equation itself, with
\(\mathcal N=-\mathbb P(u\cdot\nabla u)\)), so
\[
\|\mathcal N\|_2^2=\|\partial_tu\|_2^2
+2\nu\langle\partial_tu,Au\rangle+\nu^2\|Au\|_2^2
=\|\partial_tu\|_2^2+\nu\frac{d}{dt}\|\nabla u\|_2^2+\nu^2\|Au\|_2^2 .
\]
Dividing by \(D=\|\nabla u\|_2^2\) and integrating (all terms continuous
on compacts of \([0,T_{\max})\)):
\[
\boxed{\;\int_0^{T'}KD\,dt
=\int_0^{T'}\frac{\|\partial_tu\|_2^2}{D}\,dt
+\nu^2\int_0^{T'}\frac{\|Au\|_2^2}{D}\,dt
+\nu\log\frac{D(T')}{D(0)}\;}
\tag{ID}
\]
with \(\|Au\|_2^2/D=H_2/H_1=N_1^2\). Every quantity is nonnegative except
the log, which is exactly the Osgood-compatible part
(\(\nu\log D=\nu z+\nu\log(2H_0)\), and \(H_0\) is bounded).

**Reading.** The bridge (B) in integrated form demands
\(\int KD\le C+C\int D(1+z)\) with \(\int D(1+z)<\infty\) whenever
\(\int D<\infty\) and \(z\) grows at most logarithmically — but (ID) says
\(\int KD\) *is*, up to the harmless log,
\(\int\|\partial_tu\|^2/D+\nu^2\int N_1^2\): the \(\dot H^1\)-bandwidth
action. Controlling \(\int N_1^2\,dt\) a priori is precisely the
half-derivative-supercritical obstruction of 3D Navier–Stokes (compare the
free identity \(\int N_0^2\,dt=\frac1{2\nu}\log(H_0(0)/H_0(T'))\), which
*is* free — one level down). None of the suggested ingredients
(\(d G_0/dt\), Leray cancellation \(T_0=0\), pressure removal — already
built into \(a_k\) —, the enstrophy equation
\(\tfrac12\dot H_1=T_1-\nu H_2\), the modal covariance gap) produces an
a-priori bound on \(\int N_1^2\): the enstrophy equation gives
\(\nu\int N_1^2=\tfrac12\log(H_1(0)/H_1(T'))+\int T_1/H_1\), and
\(\int T_1/H_1\) is the (sign-indefinite, supercritical)
vortex-stretching action. **Dynamic depletion of the Osgood excess is
therefore equivalent to the regularity problem itself; there is no
shortcut through (ID).** As a free corollary (not a new candidate):
\(\int_0^{T}\bigl(\|\partial_tu\|_2^2/D+\nu^2N_1^2\bigr)dt<\infty
\Rightarrow\) regularity, by (ID) + the dichotomy — recorded only as a
restatement of \(\int KD<\infty\).

## 4. Defect relevance test (both directions fail)

Let \(\Delta(t)=G_0/(2\nu H_1)-z'\ge0\) be the (I.4) gap.

*Route 1: \(K\le C(1+z)+C\Delta\) with \(\int\Delta D<\infty\).* Since
\(\Delta\le G_0/(2\nu H_1)=KD/(2\nu)\), the proposed majorisation is
implied by (indeed, on the counterexample family it degenerates to)
\(K\lesssim KD^2\)-type statements with non-constant \(C\); and
\(\int\Delta D\) is controlled by no identity: on broadband fields
(sprint C, O-8′) the relative deficit is order one, so
\(\Delta D\asymp KD^2/(2\nu)\), *more* supercritical than \(\int KD\).
*Route 2: \(\int KD=\infty\Rightarrow\int\Delta D=\infty\), contradicted
by an identity bounding \(\int\Delta D\).* No such identity exists in the
available ledger: \(\Delta\) contains the square-completion piece
\(\asymp(2\nu/\mu)V_r\), and \(\int VD/\mu\) involves \(H_2\) — again the
same supercritical object as (ID). Moreover sprint C's exact family shows
the complementary danger: \(\Delta\) can vanish identically fast
(single-shell collapse, deficit \(\to0\)), so \(\int KD=\infty\) does
**not** force \(\int\Delta D=\infty\) on cascades that run near
saturation. Both directions are structurally closed. \([\)no search was
run for a uniform \(\Delta\)-floor, per instruction\(]\)

## 5. Exact short-time falsification (persistence of the excess)

The static family could still be dynamically irrelevant if Navier–Stokes
evolution destroyed the large-\(K\) configuration instantaneously (then a
uniform \(R\) with small \(\int R\) could absorb measure-zero spikes).
Evolving the \(N=8\) coherent field (dealiased Fourier–Galerkin RK4,
cutoff 20, grid 96, \(\nu=1/40\)):

| run | amplitude | dt-drift | \(K(0)\) | \(\min_{[0,\,0.5N^{-2}]}K\) | \(z(0)\to z(0.5N^{-2})\) | label |
|---|---|---|---|---|---|---|
| moderate | 4 | \(1.5\times10^{-4}\) | 2.0372 | 0.7051 | 1.807 → 2.660 | float, resolved |
| strong | 16 | \(1.8\times10^{-2}\) | 2.0372 | 0.0127 | 1.807 → 5.207 | float, marginal resolution |

At moderate amplitude (nonlinear time comparable to the parabolic window —
the regime relevant for critical cascades) \(K\) retains **34.6%** of its
initial value across the whole window while the bandwidth grows
(\(z\): 1.81 → 2.66): the large-\(K\) configuration is *not* an
instantaneous transient, so no measure-zero spike absorption is available
to \(R\). At strong amplitude the coherent configuration burns out within
the window (\(K\to0.013\)) while \(z\) jumps to 5.2 — fast scrambling plus
dissipation; recorded honestly with its marginal \(1.8\times10^{-2}\)
step-halving drift.

**Scope limitation, stated plainly.** At \(N=8\) the static excess is not
yet positive (\(K/(1+z)=0.73<1\)); the crossing occurs from \(N\approx16\)
(statically, §2). A resolved evolution at \(N\ge16\) (grid \(\ge144^3\))
was not run in this sprint, so §5 establishes the *persistence mechanism*
(order-one retention over the parabolic window at critical balance), not a
dynamic realisation of the excess itself. The quantitative weight of the
KILL therefore rests on §2 (exact-anchored static family, which kills
every Osgood \(\Phi\) statically) and §3 (the exact identity closing the
dynamic-depletion escape structurally); §5 rules out the "instantaneous
transient" rescue at the one scale tested. \([\)binary64 with
step-halving; no validated integrator of this size exists in the
repository\(]\)

## 6. Decision

- **GO** — not attained: no uniform Osgood estimate exists (§2, §5).
- **CONDITIONAL** — not applicable: no restricted solution class emerged
  on which (B) holds nontrivially; the identity (ID) shows any such class
  must already control \(\int N_1^2\), i.e. be a known-regular class.
- **KILL** — **adopted.** The Osgood-level improvement of \(\Lambda\) is
  impossible: statically (\(K_{\max}(z)\asymp e^{z}\), exact anchor),
  dynamically ((ID) localises the obstruction in the supercritical
  \(\dot H^1\)-bandwidth action), and at the defect level (§4).

**Standing value of \(\Lambda\) after this gate** (no new claims): the
exact identities (I.1)–(I.4), Lemma K, the monotone (M), and the
dichotomy \(\int KD\,dt<\infty\Rightarrow\) regularity — now read through
(ID) as the equivalence \(\int KD<\infty\Leftrightarrow
\int\|\partial_tu\|^2/D+\nu^2\int N_1^2<\infty\) — remain valid as a
**diagnostic criterion**, strictly weaker as a hypothesis than the
classical critical actions it refines, and exactly computable by the
repository's rational ledger. The Clay-(B) route through \(\Lambda\) is
closed at the Osgood level; any future route must control the
supercritical action itself, which is the original problem.

*Binding non-claims: no Clay statement in either direction; §2/§5 numbers
are binary64 diagnostics anchored by one exact rational point; nothing
here proves blow-up, and the counterexample family is a family of initial
data, not a singular solution.*
