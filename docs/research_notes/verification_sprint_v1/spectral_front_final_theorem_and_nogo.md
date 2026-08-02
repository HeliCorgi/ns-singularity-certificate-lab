# The spectral-front Λ lane, frozen: final theorems and no-go

> **Superseded in detail by the refereed paper deliverables**
> ([docs/paper_lambda_dichotomy/](../../paper_lambda_dichotomy/theorem_statement.md)).
> Two adversarial referees found and the paper fixes: (1) the
> "equivalence" reading of T4 is false — only the one-sided AM–GM
> comparison \(\int KD\le2\int(\|\partial_tu\|^2/D+\nu^2N_1^2)\) holds
> (explicit decaying counterexample); (2) the positivity of \(H_0\)
> asserted in §0 requires the bootstrap argument through Lemma 7 (the
> Grönwall/forward-uniqueness sketch here is circular); (3) the
> dichotomy must be stated as {action finite ⇒ global} vs {action
> infinite}, without exhaustiveness of a blow-up/global split; (4) the
> \(S_N\) bounds are stated in the paper with crude proven constants.
> This file is retained as the freeze-directive record; where it and
> the paper differ, the paper is authoritative.

**PART I of the freeze directive. No new variants of \(\Lambda\), \(K\),
\(N_r\), Osgood, or deficit quantities are introduced; this document fixes
the existing ones as final statements with proofs, and classifies every
item as THEOREM / CONDITIONAL CRITERION / NO-GO. The designation "Clay
proof candidate" is removed everywhere (item 6); the candidate document
carries a pointer to this file as its final status.**

## 0. Classification summary

| item | statement | class |
|---|---|---|
| T1 | front identities (I.1)–(I.4) with gap decomposition and equality cases | **THEOREM** |
| T2 | monotone \(\Lambda\) and the bandwidth–dissipation dichotomy | **THEOREM** (as an implication) / **CONDITIONAL CRITERION** (as a regularity test: its hypothesis is not verifiable a priori) |
| T3 | Osgood closure (product form, \(D=0\) safe) | **THEOREM** (conditional statement, fully proven) |
| T4 | dynamic identity \(\int KD=\int\|\partial_tu\|_2^2/D+\nu^2\int N_1^2+\nu\log(D(T')/D(0))\) | **THEOREM** |
| T5 | coherent-family exact laws: \(H_0=\tfrac23c_N^2\|v_0\|^2T_N\), \(H_1=\tfrac23c_N^2\|v_0\|^2S_N\), \(N_0^2=S_N/T_N\asymp N\), \(u_N(0)=\tfrac23c_NS_Nv_0\), \(\|u_N\|_\infty\asymp N\) | **THEOREM** |
| L\* | capacity lower bound \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\), **smoothly** truncated family \(\chi(\lvert k\rvert/N)P_kv_0/\lvert k\rvert^2\) | **PROVEN** (2026-08-02), sharp exponent, uniform in \(\chi\) and \(v_0\) — `docs/paper_lambda_dichotomy/lstar/lstar_proof_main.md` Thm 7.1(3). Constants \(c_0,N_*\) non-effective |
| L\*(sharp) | the same bound for the **sharply** truncated family \(P_kv_0/\lvert k\rvert^2\) on \(1\le\lvert k\rvert\le N\) | **OPEN — and no longer used by anything.** Exact finite-\(N\) certificates and float continuation (§6) |
| T6 | impossibility of pointwise Osgood closure (\(\Phi\) must be exponential-class; \(\int ds/\Phi<\infty\)) | **NO-GO, UNCONDITIONAL** (superseded row; was "conditional on L\*") |
| — | "Clay proof candidate" designation for this lane | **REMOVED** |

Throughout: \(\mathbb T^3\), \(u=\sum_k\hat u_ke^{ik\cdot x}\), zero mean,
\(k\cdot\hat u_k=0\), \(\hat u_{-k}=\overline{\hat u_k}\);
\(\mathcal N=-\mathbb P(u\cdot\nabla u)\); paired-real modal ledger
\(e_k\ (=\) per-mode \(L^2\) contribution\()\), \(a_k=\) modal growth inner
product, \(n_k=\) per-mode \(L^2\) contribution of \(\mathcal N\);
\(H_r=\sum x^re_k\), \(T_r=\sum x^ra_k\), \(G_r^{\rm in}=\sum x^rn_k\)
(in-support), \(G_0^{\rm full}=\|\mathcal N\|_2^2\)-normalised full moment
(\(G_r^{\rm in}\le G_r^{\rm full}\)); \(x_k=|k|^2\), \(N_r^2=H_{r+1}/H_r\),
\(z=\log N_0^2\ge0\), \(D=\|\nabla u\|_2^2\),
\(K=\|\mathbb P(u\cdot\nabla u)\|_2^2/\|\nabla u\|_2^4\). For a strong
solution \(u\in C([0,T_{\max});H^m)\cap C^1([0,T_{\max});H^{m-2})\),
\(m>5/2\), every moment used below with index \(r\le m-2\) is finite and
continuously differentiable in \(t\), and term-by-term differentiation of
the (absolutely and locally uniformly convergent) lattice sums is
justified by the \(C^1H^{m-2}\)-bound; we use only \(r\in\{0,1\}\), for
which \(m>5/2\) suffices. If \(u_0\neq0\) then \(H_0(t),H_1(t)>0\) for all
\(t<T_{\max}\) (else \(u(t)=0\) and forward uniqueness forces
\(u\equiv0\)), so \(z\) is finite and \(C^1\).

## 1. T1 — front identities (THEOREM)

**Theorem 1.** Let \(u\) be a divergence-free zero-mean real
trigonometric field (or a strong solution at a fixed time), \(r\ge0\) with
\(H_r,H_{r+1}>0\), \(\nu>0\). With \(p_r(k)=x_k^re_k/H_r\),
\(\mu=N_r^2\), \(g_k=a_k/e_k\) on \(e_k>0\),
\(V_r=\sum p_r(k)(x_k-\mu)^2\), \(\mathrm{Cov}=T_{r+1}/H_r-\mu T_r/H_r\):

(I.1) \(\dfrac{d}{dt}\log N_r^2=\dfrac{2\,\mathrm{Cov}}{\mu}
-\dfrac{2\nu V_r}{\mu}\) along \(\tfrac12\dot H_r=T_r-\nu H_{r+1}\);

(I.2) \(|\mathrm{Cov}|\le\sqrt{V_rG_r^{\rm in}/H_r}\);

(I.3) \(\dfrac{d}{dt}\log N_r^2\le\dfrac{G_r^{\rm in}}{2\nu H_{r+1}}
\le\dfrac{G_r^{\rm full}}{2\nu H_{r+1}}\);

(I.4) the gap \(G_r^{\rm in}/(2\nu H_{r+1})-\frac{d}{dt}\log N_r^2
=\Gamma^{\rm CS}_r+\Gamma^{\rm SC}_r\) with
\(\Gamma^{\rm CS}_r=\frac2\mu[\sqrt{V_rG_r^{\rm in}/H_r}-\mathrm{Cov}]\ge0\),
\(\Gamma^{\rm SC}_r=\frac{2\nu}{\mu}\bigl[\sqrt{V_r}
-\frac1{2\nu}\sqrt{G_r^{\rm in}/H_r}\bigr]^2\ge0\), whose sum is rational
whenever the data are; \(\Gamma^{\rm CS}_r=0\) iff
\(a_k=\gamma(x_k-\mu)\,x_k^r e_k\)-aligned with
\(a_k^2=e_kn_k\) modally, and \(\Gamma^{\rm SC}_r=0\) iff
\(V_r=G_r^{\rm in}/(4\nu^2H_r)\).

*Proof.* (I.1): \(\frac{d}{dt}\log N_r^2
=\frac{\dot H_{r+1}}{H_{r+1}}-\frac{\dot H_r}{H_r}
=2\Bigl[\frac{T_{r+1}}{H_{r+1}}-\frac{T_r}{H_r}\Bigr]
-2\nu\Bigl[\frac{H_{r+2}}{H_{r+1}}-\frac{H_{r+1}}{H_r}\Bigr]\).
The first bracket equals \((T_{r+1}-\mu T_r)/H_{r+1}=\mathrm{Cov}/\mu\);
the second equals \((H_rH_{r+2}-H_{r+1}^2)/(H_rH_{r+1})=V_r/\mu\) since
\(H_rH_{r+2}-H_{r+1}^2=H_r^2V_r\) (direct expansion of
\(V_r=H_{r+2}/H_r-\mu^2\), which equals the centred sum by exact
bookkeeping). (I.2): modal Cauchy–Schwarz twice gives
\(|a_k|=\tfrac12|c_k\!\cdot\!c_k^{\mathcal N}+s_k\!\cdot\!s_k^{\mathcal N}|
\le\sqrt{e_kn_k}\); then
\(|\mathrm{Cov}|\le\frac1{H_r}\sum x^r|x-\mu|\sqrt{e_kn_k}
\le\sqrt{\frac{\sum x^r(x-\mu)^2e_k}{H_r}}\sqrt{\frac{\sum x^rn_k}{H_r}}
\,\sqrt{H_r}\cdot\frac1{\sqrt{H_r}}=\sqrt{V_rG_r^{\rm in}/H_r}\)
(vector Cauchy–Schwarz on \((\sqrt{x^re_k}\,|x-\mu|)\),
\((\sqrt{x^rn_k})\)). (I.3): insert (I.2) into (I.1) and complete the
square in \(\sqrt{V_r}\):
\(\sup_{V\ge0}[\sqrt{G^{\rm in}_r/H_r}\sqrt V-\nu V]
=G_r^{\rm in}/(4\nu H_r)\); divide by \(\mu/2\). The full-moment bound
follows from \(G^{\rm in}\le G^{\rm full}\) (adding nonnegative
off-support terms). (I.4): algebraic expansion; the square roots cancel
in the sum (computed and machine-verified as an exact rational identity
in `spectral_front_monotone.front_gap_identity` and the sprint-C
telescoping certificates). Equality cases: chase the two Cauchy–Schwarz
equalities and the square-completion equality. ∎

## 2. T2 — monotone and dichotomy

**Theorem 2.** For every strong solution and \(t<T_{\max}\):
\[
\Lambda(t)=\log N_0^2(t)-\frac1{2\nu}\int_0^tK\,D\,ds
\quad\text{is non-increasing},
\]
and exactly one of: (i) \(\int_0^{T_{\max}}KD\,dt<\infty\), in which case
\(z\) is bounded, \(u\in L^\infty(0,T_{\max};H^1)\subset
L^\infty L^6\) (Serrin exponent \(2/\infty+3/6=1/2<1\)) and the solution
is regular (extends past any finite \(T_{\max}\)); or (ii)
\(T_{\max}<\infty\) and \(\int_0^{T_{\max}}KD\,dt=\infty\).

*Proof.* Monotonicity: (I.3) at \(r=0\) with
\(G_0^{\rm full}/(2\nu H_1)=KD/(2\nu)\) (the full-moment form of the
convention lemma: \(K D=\|\mathcal N\|_2^2/\|\nabla u\|_2^2\), and the
pair-energy convention makes \(G_0^{\rm full}=\|\mathcal N\|_2^2\),
\(2H_1=\ldots\) — the normalisation cancels in the displayed ratio, as
calibrated exactly against the independent float pipeline). Dichotomy:
if \(\int KD<\infty\), integrate the monotone:
\(z(t)\le z(0)+\frac1{2\nu}\int_0^{T_{\max}}KD=:z_*\); then
\(H_1\le e^{z_*}H_0(0)\); Serrin \((6,\infty)\) applies; regularity and
extension follow. If \(T_{\max}<\infty\) and \(\int KD<\infty\) we just
showed extension — contradiction; hence (ii). ∎

**Classification.** As an implication, T2 is a theorem. As a regularity
test it is a **CONDITIONAL CRITERION**: the hypothesis
\(\int KD\,dt<\infty\) is weaker than the classical critical actions
(\(KD\le\|u\|_\infty^2\) and \(KD\le C_S^2\|\nabla u\|_{L^3}^2\), both by
Hölder), so the criterion is strictly stronger than the Serrin
\((\infty,2)\) and vorticity \((3,2)\) tests it refines — but T6 below
shows it cannot be upgraded to an unconditional theorem by an
Osgood-type argument.

## 3. T3 — Osgood closure (THEOREM; the conditional GO branch)

**Theorem 3.** If a strong solution satisfies, a.e. on \([0,T_{\max})\),
\(K D\le\Phi(z)D+R\) with \(\Phi>0\) nondecreasing,
\(\int^\infty ds/\Phi=\infty\), \(R\ge0\),
\(\int_0^{T_{\max}}R\,dt<\infty\), then \(z\) is bounded and
\(T_{\max}=\infty\).

*Proof.* As in the Osgood-gate note §1 (reproduced to fix conventions):
\(\Omega(z)=\int_0^zds/\Phi\) satisfies
\(\frac{d}{dt}\Omega(z)\le D/(2\nu)+R/(2\nu\Phi(0))\), both integrable
(\(\int D\le E(0)/\nu\)); hence \(\Omega(z)\) bounded, \(z\) bounded
(\(\Omega(\infty)=\infty\)), and Theorem 2(i) applies. The product form
needs no division by \(D\); \(z\ge0\) on \(\mathbb T^3\); \(H_0,H_1>0\)
as in §0. ∎

## 4. T4 — dynamic identity (THEOREM)

**Theorem 4.** For every strong solution and \(0\le T'<T_{\max}\)
(with \(u_0\neq0\)):
\[
\int_0^{T'}KD\,dt=\int_0^{T'}\frac{\|\partial_tu\|_2^2}{D}\,dt
+\nu^2\int_0^{T'}N_1^2\,dt+\nu\log\frac{D(T')}{D(0)} .
\]

*Proof.* The equation gives \(\mathcal N=\partial_tu+\nu Au\),
\(A=-\Delta\), all three terms in \(C([0,T'];L^2)\) (indeed \(H^{m-2}\)).
Expand \(\|\mathcal N\|_2^2=\|\partial_tu\|_2^2
+2\nu\langle\partial_tu,Au\rangle+\nu^2\|Au\|_2^2\), and
\(\langle\partial_tu,Au\rangle=\langle\nabla\partial_tu,\nabla u\rangle
=\tfrac12\frac{d}{dt}\|\nabla u\|_2^2\) (justified since
\(u\in C^1H^{m-2}\cap CH^m\), \(m>5/2\), so \(\nabla u\in C^1L^2\)).
Divide by \(D(t)>0\) (positivity as in §0), integrate; the middle term
integrates exactly to \(\nu\log(D(T')/D(0))\) because
\(t\mapsto D(t)\) is \(C^1\) and positive; \(\|Au\|_2^2/D=H_2/H_1=N_1^2\)
by Parseval. Boundary terms: only the displayed logarithm; no other
boundary term arises since the identity is integrated, not integrated by
parts in time. \(D=0\) cannot occur for \(u_0\neq0\); for \(u_0=0\) all
terms vanish identically. ∎

**Reading (fixed).** The only Osgood-compatible term is the logarithm;
\(\int KD<\infty\iff\int\|\partial_tu\|^2/D+\nu^2\int N_1^2<\infty\)
(given \(\log D\) bounded, which follows from either side plus energy).
The lane's criterion is exactly an a-priori bound on the
\(\dot H^1\)-bandwidth action — the half-derivative-supercritical
object. This is the precise sense in which the lane is *frozen*: any
further progress must bound \(\int N_1^2\), which is not a
\(\Lambda\)-variant question.

## 5. T5 — the coherent critical-spectrum family (THEOREM)

Fix \(v_0=(1,2,3)\) (any \(v_0\in\mathbb Z^3\setminus\{0\}\) works),
\(N\ge2\), \(c_N>0\), and let
\[
\widehat u_N(k)=c_N\,\frac{P_kv_0}{|k|^2},\qquad 1\le|k|\le N,
\qquad \widehat u_N(-k)=\overline{\widehat u_N(k)}=\widehat u_N(k).
\]

**Theorem 5.** \(u_N\) is real, zero-mean, exactly divergence-free, and
with \(S_N=\sum_{1\le|k|\le N}|k|^{-2}\),
\(T_N=\sum_{1\le|k|\le N}|k|^{-4}\):
\[
H_0=\tfrac23c_N^2\|v_0\|^2\,T_N,\qquad
H_1=\tfrac23c_N^2\|v_0\|^2\,S_N,\qquad
N_0^2=\frac{S_N}{T_N},
\]
\[
u_N(0)=\tfrac23\,c_N\,S_N\,v_0
\quad(\text{exact}),\qquad
\tfrac23S_N\,\|v_0\|\,c_N\le\|u_N\|_\infty\le S_N\,\|v_0\|\,c_N .
\]
Moreover \(4\pi(N-2)\le S_N\le 4\pi N+C_0\) (cube-comparison with the
integral \(\int_{|x|\le N}|x|^{-2}dx=4\pi N\); \(C_0\) absorbs the first
shells) and \(T_N\le T_\infty<\infty\); hence
\(N_0^2\in[c_1N,c_2N]\) with explicit \(c_1,c_2>0\), and the family
saturates the Lemma-K ratio:
\(\|u_N\|_\infty^2/\|\nabla u_N\|_2^2\ge\tfrac23S_N\gtrsim N\).

*Proof.* Divergence-freeness: \(k\cdot P_kv_0=0\) per mode. Reality:
coefficients are real and \(\pm k\)-symmetric. The band
\(B_N=\{1\le|k|\le N\}\) is invariant under coordinate sign flips and
permutations; therefore, for any radial weight \(f\),
\(\sum_{k\in B_N}k_ik_jf(|k|^2)=\delta_{ij}\tfrac13\sum_{k\in B_N}|k|^2f(|k|^2)\)
(off-diagonal terms cancel under \(k_i\mapsto-k_i\); the three diagonal
sums are equal by permutation and add to \(\sum|k|^2f\)). With
\(|P_kv_0|^2=\|v_0\|^2-(k\cdot v_0)^2/|k|^2\):
\(H_0=c_N^2\sum|P_kv_0|^2|k|^{-4}
=c_N^2[\|v_0\|^2T_N-\tfrac13\|v_0\|^2T_N]=\tfrac23c_N^2\|v_0\|^2T_N\),
and identically for \(H_1\) with \(|k|^{-2}\). The point value:
\(u_N(0)=\sum_k\widehat u_N(k)=c_N[v_0S_N-\sum_kk\,(k\cdot v_0)|k|^{-4}]
=c_N[v_0S_N-\tfrac13v_0S_N]=\tfrac23c_NS_Nv_0\), giving the sup-norm
lower bound; the upper bound is the triangle inequality
\(\|u\|_\infty\le\sum|\widehat u(k)|\le c_N\|v_0\|S_N\). The \(S_N\)
bounds: each unit cube centred at \(k\in B_N\) is contained in
\(\{|x|\le N+\sqrt3/2\}\) and contains points with
\(|x|\ge|k|-\sqrt3/2\); comparing \(\sum|k|^{-2}\) with
\(\int|x|^{-2}\) over the ball, shell by shell, yields the display (the
first two shells are estimated directly). ∎

All identities of Theorem 5 are additionally machine-verified as exact
`Fraction` equalities at \(N=4,6,8\)
(`coherent_family_certificates.md`, certificate JSON with digest).

## 6. L\* — the capacity lower bound (SUPERSEDED: now PROVEN for the smooth family)

> **Supersession notice (2026-08-02).** This section is retained as the
> freeze-directive record of what was open at the time. It is **out of
> date**. The capacity bound
> \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) is now **PROVEN**,
> with the sharp exponent \(3\), for the *smoothly* truncated family
> \(\widehat u_N(k)=\chi(|k|/N)P_kv_0/|k|^2\), uniformly in the admissible
> cutoff \(\chi\) and in \(v_0\in\mathbb R^3\setminus\{0\}\)
> (`docs/paper_lambda_dichotomy/lstar/lstar_proof_main.md`, Theorem
> 7.1(3); adversarially refereed, verdict (i)). The constants \(c_0\) and
> \(N_*\) are **non-effective**. The route below — a rigorous asymptotic
> lower bound on the explicit lattice triple sum, via a *constant*
> sweeping vector — is not the route that worked, and is in fact now a
> **proven dead end**: Theorem A of
> `lstar/lstar_direct_route_and_weakening.md` §A shows every
> constant-vector sweeping split, in any Hölder pairing and any constant
> pairing direction, is capped at \(\gtrsim N\), two powers short. The
> proof that works pairs against a *concentrated* divergence-free test
> field \(\psi_N(x)=N^{3/2}\Psi(Nx)\) at the concentration scale.
> The statement below, for the **sharply** truncated family, remains
> **OPEN** — and is used by nothing.

**Lemma L\* (statement; open for the sharp family).** There is \(c_0>0\) with
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) for all large \(N\)
(equivalently \(K(u_N)\ge c_0'N_0^2\)).

Status: **not proven**. Support: (a) exact rational certificates at
small \(N\) (\(K=0.788411\ldots\) at \(N=4\), exact = float to printed
precision; extended certificates in `coherent_family_certificates.md`);
(b) float continuation to \(N=32\) with \(K/N_0^2\) increasing
0.259 → 0.396 and Leray retention
\(\|\mathbb P(u\cdot\nabla u)\|/\|u\cdot\nabla u\|\) increasing
0.780 → 0.867; (c) the proven sweeping structure: for constant \(c\),
\(\mathbb P((c\cdot\nabla)u)=(c\cdot\nabla)u\) (the constant-coefficient
derivative commutes with each \(P_k\) and preserves divergence-freeness),
so the Leray projection cannot remove the sweeping component, and the
pairing bound
\(\|\mathbb P(u\cdot\nabla u)\|_2\ge|\langle(u\cdot\nabla)u,
(v_0\cdot\nabla)u\rangle|/\|(v_0\cdot\nabla)u\|_2\) is exactly
computable (certificates, ibid.). What is missing is a rigorous
asymptotic lower bound on that explicit lattice triple sum — a
combinatorial estimate, deliberately **not** attempted here (freeze
directive).

## 7. T6 — impossibility of pointwise Osgood closure (NO-GO, now UNCONDITIONAL)

> **Supersession notice (2026-08-02).** The hypothesis "Assume L\*" below
> has been **discharged**. Theorem 6 is now unconditional: the capacity
> bound holds for the smoothly truncated family (§6 notice), and the
> hypothesis of Theorem 6 quantifies over *all* divergence-free zero-mean
> real trigonometric fields, so one admissible family suffices. The
> authoritative statement is Theorem O of
> `docs/paper_lambda_dichotomy/theorem_statement.md`. Its constant \(c\)
> is non-effective. The honest headline is **"the no-go is
> unconditional"**, *not* "L\* is proven": the literal L\*, for the
> sharply truncated family used below, is still open.
>
> Two further differences from the paper's current proof: it uses
> \(2^{20}\)-adic rather than unit steps in \(N\) (the "asymptotically
> dense \(s_N\)" step below is replaced by an explicit two-sided bound
> \(N/44544\le N_0^2\le\tfrac{64}3N\), giving gaps in
> \([0.0985,\,27.63]\)), and it uses the smooth family's exact laws.

**Theorem 6 (as frozen; the hypothesis is now proven, see notice).**
Assume L\*. Let \(\Phi\) be nondecreasing with
\(K(u)\le\Phi(\log N_0^2(u))\) for every divergence-free zero-mean real
trigonometric field \(u\). Then there are \(c,s_0>0\) with
\(\Phi(s)\ge c\,e^{s}\) for all \(s\ge s_0\); consequently
\(\int^\infty ds/\Phi(s)<\infty\): **no Osgood-admissible \(\Phi\)
satisfies a uniform pointwise bound, and Theorem 3 can never be
activated unconditionally.**

*Proof.* By L\* and Theorem 5, the family has
\(K(u_N)\ge c_0'N_0^2(N)=c_0'e^{s_N}\) at \(s_N=\log(S_N/T_N)\).
\(\Phi(s_N)\ge c_0'e^{s_N}\). The sequence \(s_N\) is asymptotically
dense: \(s_{N+1}-s_N=\log(S_{N+1}/S_N)-\log(T_{N+1}/T_N)\to0\) (both
ratios \(\to1\); \(S_{N+1}-S_N\asymp1\ll S_N\), \(T\) convergent). For
\(s\in[s_N,s_{N+1}]\), monotonicity gives
\(\Phi(s)\ge\Phi(s_N)\ge c_0'e^{s_N}\ge c_0'e^{s-(s_{N+1}-s_N)}
\ge\tfrac{c_0'}{2}e^{s}\) for \(N\) large. Then
\(\int^\infty ds/\Phi\le\tfrac2{c_0'}\int^\infty e^{-s}ds<\infty\). ∎

Together with Theorem 4 (which reduces every dynamic-depletion escape to
the supercritical \(\dot H^1\)-bandwidth action) and the defect analysis
of the Osgood-gate note §4, this closes the Osgood route at every level
at which it was proposed. ~~The remaining logical gap of the no-go is
exactly L\*~~ — **there is no longer a logical gap**: the capacity bound
is proven for the smooth family (§6 notice), so the no-go is
unconditional. What remains open is only the *sharp*-family capacity
bound, which nothing uses.

## 8. What the lane now is

- **THEOREMS**: T1, T3, T4, T5 (+ Lemma K, the \(\varepsilon\)-identity
  with degradation factor \(1+\varepsilon H_0/H_1\), and the exact
  telescoping certificates).
- **CONDITIONAL CRITERION**: T2 — the strongest known member of its
  family of critical-action regularity tests in this repository,
  equivalent by T4 to an a-priori bound on
  \(\int\|\partial_tu\|^2/D+\nu^2\int N_1^2\).
- **NO-GO**: T6, **unconditional** (since 2026-08-02; was "conditional on
  L\*"): no Osgood upgrade of T2 exists. Its constants \(c_0,N_*\) are
  non-effective. The sharp-family capacity bound (literal L\*) is still
  open and is used by nothing.
- **Not claimed**: any Clay statement, in either direction. The former
  "PROOF CANDIDATE" and "Clay proof candidate" designations are
  withdrawn; the lane is frozen as a diagnostic engine.
