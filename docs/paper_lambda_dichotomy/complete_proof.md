# Complete proof

Conventions and the statement are in
[theorem_statement.md](theorem_statement.md). Lemmas appear in
dependency order — with the one referee-mandated exception that the
positivity Lemma 0b is placed after Lemma 7, on which it depends; each
proof is self-contained modulo the three external classical inputs
F1–F3 listed in the [audit](dependency_and_gap_audit.md) (local strong
solutions; subcritical-Serrin-class regularity; mean-zero Sobolev),
which are framework inputs, not steps invented here.

Throughout, \(u\) is the maximal strong solution with
\(u_0\in H^m_\sigma\), \(m>5/2\), \(u_0\neq0\). All modal sums run over
the **full lattice** \(k\in\mathbb Z^3\setminus\{0\}\) (no pairing
convention): \(\varepsilon_k=|\hat u_k|^2\),
\(\mathcal N=-\mathbb P(u\cdot\nabla u)\),
\(\alpha_k=\mathrm{Re}\langle\hat u_k,\hat{\mathcal N}_k\rangle\),
\(\eta_k=|\hat{\mathcal N}_k|^2\), \(x_k=|k|^2\);
\(H_r=\sum x_k^r\varepsilon_k\), \(T_r=\sum x_k^r\alpha_k\), and
\(\sum_k\eta_k=\|\mathcal N\|_2^2\). Note \(\alpha_k=0\) whenever
\(\varepsilon_k=0\).

---

**Lemma 0a (regularity of the moment functionals).**
On \([0,T_{\max})\): \(H_0,H_1,H_2\) are finite and continuous;
\(H_0,H_1\in C^1\) with
\[
\dot H_0=2\langle u,\partial_tu\rangle,\qquad
\dot H_1=2\langle Au,\partial_tu\rangle,\qquad A=-\Delta;
\]
and \(\partial_tu=\mathcal N-\nu Au\) holds in \(L^2\) at every \(t\).
On any subinterval where \(H_0>0\), also \(H_1>0\) (else
\(\hat u_k\equiv0\) and \(H_0=0\)), and \(z=\log(H_1/H_0)\) is \(C^1\)
there.

*Proof.* Finiteness/continuity: \(u\in C([0,T');H^m)\) with
\(m>5/2>2\) gives \(H_r\in C\) for \(r\le2\) (composition of the
continuous map \(t\mapsto u(t)\in H^2\) with the squared norm).
Differentiability: for \(r\in\{0,1\}\), the difference quotient of
\(H_r\) equals \(\langle A^{r}(u(t+h)+u(t)),(u(t+h)-u(t))/h\rangle_{L^2}\)
with \(A^ru\in C([0,T');L^2)\) (\(Au\in C(H^{m-2})\subset C(L^2)\))
and \((u(t+h)-u(t))/h\to\partial_tu\) in \(H^{m-2}\subset L^2\); pass
to the limit. The equation holds in \(H^{m-2}\subset L^2\) pointwise in
\(t\) by the definition of strong solution. ∎

*(Positivity of \(H_0\) on all of \([0,T_{\max})\) is Lemma 0b below,
proven after Lemma 7; Lemmas 1–4, 6, 7 are pointwise statements that
require only \(H_0,H_1>0\) at the time considered, so no circularity
arises.)*

**Lemma 1 (ledger identities; energy neutrality).**
\(\tfrac12\dot H_0=T_0-\nu H_1\) and \(\tfrac12\dot H_1=T_1-\nu H_2\),
with \(T_0=0\); consequently \(\dot H_0=-2\nu H_1\le0\) (energy
equality) and \(\int_0^{T_{\max}}D\,dt\le H_0(0)/(2\nu)\).

*Proof.* Substitute \(\partial_tu=\mathcal N-\nu Au\) into Lemma 0:
\(\dot H_r=2\langle A^ru,\mathcal N\rangle-2\nu\langle A^ru,Au\rangle
=2T_r-2\nu H_{r+1}\) for \(r=0,1\), where
\(\langle A^ru,\mathcal N\rangle=T_r\) by Parseval (definition of
\(a_k\)) and \(\langle A^ru,Au\rangle=H_{r+1}\). Energy neutrality:
\(T_0=\langle u,\mathcal N\rangle=-\langle u,\mathbb P(u\cdot\nabla
u)\rangle=-\langle u,(u\cdot\nabla)u\rangle
=-\tfrac12\int u\cdot\nabla|u|^2=0\), using \(\mathbb P u=u\), the
divergence theorem on \(\mathbb T^3\), and \(\nabla\cdot u=0\); the
integrand manipulations are classical for \(u(t)\in H^m\), \(m>5/2\)
(all products are in \(C^1\)). ∎

**Lemma 2 (identity (I.1) at \(r=0\)).** With
\(p(k)=\varepsilon_k/H_0\), \(\mu=N_0^2\),
\(\mathrm{Cov}=\sum_k(x_k-\mu)\alpha_k/H_0=T_1/H_0-\mu T_0/H_0
=T_1/H_0\), and \(V=\sum_kp(k)(x_k-\mu)^2\):
\[
\frac{d}{dt}\log N_0^2=\frac{2\,\mathrm{Cov}}{\mu}-\frac{2\nu V}{\mu}.
\]

*Proof.* \(\frac{d}{dt}\log N_0^2=\dot H_1/H_1-\dot H_0/H_0
=2\bigl[\tfrac{T_1}{H_1}-\tfrac{T_0}{H_0}\bigr]
-2\nu\bigl[\tfrac{H_2}{H_1}-\tfrac{H_1}{H_0}\bigr]\) by Lemma 1.
First bracket \(=(T_1-\mu T_0)/H_1=\mathrm{Cov}/\mu\) (with \(T_0=0\),
\(\mathrm{Cov}=T_1/H_0\)); second bracket
\(=(H_0H_2-H_1^2)/(H_0H_1)=V/\mu\), since
\(H_0H_2-H_1^2=H_0^2V\) by expanding
\(V=\sum p(x-\mu)^2=H_2/H_0-\mu^2\). ∎

**Lemma 3 (closable covariance bound (I.2)).**
\(|\mathrm{Cov}|\le\sqrt{V\,\|\mathcal N\|_2^2/H_0}\).

*Proof.* Modal Cauchy–Schwarz on \(\mathbb C^3\):
\(|\alpha_k|\le|\hat u_k||\hat{\mathcal N}_k|
=\sqrt{\varepsilon_k\eta_k}\). Hence
\(|\mathrm{Cov}|=\bigl|\sum_k(x_k-\mu)\alpha_k\bigr|/H_0
\le\frac1{H_0}\sum_k|x_k-\mu|\sqrt{\varepsilon_k}\cdot\sqrt{\eta_k}
\le\frac1{H_0}\Bigl(\sum_k(x_k-\mu)^2\varepsilon_k\Bigr)^{1/2}
\Bigl(\sum_k\eta_k\Bigr)^{1/2}
=\sqrt{\frac{V\,\|\mathcal N\|_2^2}{H_0}}\)
(vector Cauchy–Schwarz applied to
\((|x_k-\mu|\sqrt{\varepsilon_k})_k\) and \((\sqrt{\eta_k})_k\), using
\(\sum(x_k-\mu)^2\varepsilon_k=H_0V\)). ∎

**Lemma 4 (square completion; Main (a)).**
\(\frac{d}{dt}\log N_0^2\le\frac{\|\mathcal N\|_2^2}{2\nu H_1}
=\frac{K\,D}{2\nu}\), hence \(\Lambda\) is non-increasing; the defect
decomposes as
\(\frac{KD}{2\nu}-\frac{d}{dt}\log N_0^2
=\Gamma^{\rm CS}+\Gamma^{\rm SC}\ge0\) with
\(\Gamma^{\rm CS}=\frac2\mu[\sqrt{V\|\mathcal N\|^2/H_0}-\mathrm{Cov}]\),
\(\Gamma^{\rm SC}=\frac{2\nu}\mu\bigl[\sqrt V-\tfrac1{2\nu}
\sqrt{\|\mathcal N\|^2/H_0}\bigr]^2\).

*Proof.* By Lemmas 2–3,
\(\frac{d}{dt}\log N_0^2\le\frac2\mu[\sqrt{\|\mathcal N\|^2/H_0}\sqrt V
-\nu V]\le\frac2\mu\cdot\frac{\|\mathcal N\|^2/H_0}{4\nu}
=\frac{\|\mathcal N\|^2}{2\nu\mu H_0}=\frac{\|\mathcal N\|^2}{2\nu H_1}\)
(maximising the concave function \(\sqrt V\mapsto\alpha\sqrt V-\nu V\)).
\(K D=\|\mathcal N\|^2/H_1\) by definition. The decomposition is
algebra (the radicals cancel in the sum, which is therefore rational for
rational data); each part is nonnegative by Lemma 3 and by being a
square. Integrating the inequality gives the monotone. ∎

**Lemma 5 (dichotomy; Main (b)).** Exactly one of:
(b1) \(\int_0^{T_{\max}}KD\,dt<\infty\), and then \(z\) is bounded,
\(u\in L^\infty(0,T_{\max};H^1)\) and \(T_{\max}=\infty\); or
(b2) \(\int_0^{T_{\max}}KD\,dt=\infty\). In particular
\(T_{\max}<\infty\Rightarrow\int_0^{T_{\max}}KD\,dt=\infty\).
*(No claim is made that (b2) forces \(T_{\max}<\infty\); a global
solution with divergent action is not excluded by this lemma.)*

*Proof.* The two alternatives are exhaustive and exclusive by
definition. Under (b1), integrating Lemma 4:
\(z(t)\le z(0)+\frac1{2\nu}\int_0^{T_{\max}}KD=:z_*<\infty\). Then
\(H_1(t)=N_0^2(t)H_0(t)\le e^{z_*}H_0(0)\) (Lemma 1: \(H_0\)
non-increasing). So \(u\in L^\infty H^1\subset L^\infty L^6\) (mean-zero
Sobolev embedding on \(\mathbb T^3\)), a subcritical Serrin class
(\(2/q+3/p=1/2<1\) with \((p,q)=(6,\infty)\)); by F2 the solution is
regular up to any finite time and extends, so \(T_{\max}=\infty\).
The final implication is the contrapositive. ∎

**Lemma 6 (exact action representation; Main (c)).** For
\(T'<T_{\max}\):
\(\int_0^{T'}KD=\int_0^{T'}\|\partial_tu\|_2^2/D
+\nu^2\int_0^{T'}N_1^2+\nu\log(D(T')/D(0))\).

*Proof.* \(\mathcal N=\partial_tu+\nu Au\) in \(C([0,T'];L^2)\).
Expand:
\(\|\mathcal N\|_2^2=\|\partial_tu\|_2^2+2\nu\langle\partial_tu,Au\rangle
+\nu^2\|Au\|_2^2\). By Lemma 0a,
\(2\langle\partial_tu,Au\rangle=\dot H_1=\dot D\). Divide by
\(D(t)>0\) (Lemma 0b) — all three quotients are continuous — and
integrate on \([0,T']\):
\(\int\|\mathcal N\|^2/D=\int\|\partial_tu\|^2/D
+\nu\int\dot D/D+\nu^2\int\|Au\|^2/D\). The middle integral is
\(\log(D(T')/D(0))\) (\(D\in C^1\), \(D>0\));
\(\|Au\|_2^2/D=H_2/H_1=N_1^2\); and \(\|\mathcal N\|^2/D=KD\). ∎

**Corollary 6′ (one-sided comparison; the equivalence FAILS).** By
AM–GM, \(\nu\dot D/D=2\nu\langle\partial_tu,Au\rangle/D
\le\|\partial_tu\|_2^2/D+\nu^2\|Au\|_2^2/D\), so
\[
\int_0^{T'}KD\,dt\;\le\;
2\int_0^{T'}\Bigl(\frac{\|\partial_tu\|_2^2}{D}+\nu^2N_1^2\Bigr)dt :
\]
finiteness of the \(\dot H^1\)-bandwidth action implies alternative
(b1). The converse is **false**: for the explicit global solution
\(u=e^{-\nu t}(0,0,\cos x_1)\) one has \(\mathcal N\equiv0\), hence
\(\int KD=0\), while
\(\int_0^{T'}(\|\partial_tu\|^2/D+\nu^2N_1^2)=2\nu^2T'\to\infty\)
(the identity balances through \(\nu\log D\to-\infty\)). The action
representation is exact; only the one-sided comparison survives as a
criterion.

**Lemma 7 (dominations; Main (d)).**
\(KD\le\|u\|_\infty^2\) and \(KD\le C_S^2\|\nabla u\|_{L^3}^2\).

*Proof.* \(KD=\|\mathcal N\|_2^2/D\). First:
\(\|\mathcal N\|_2\le\|(u\cdot\nabla)u\|_2\le\|u\|_\infty\|\nabla u\|_2\)
(\(\mathbb P\) is an orthogonal projection on \(L^2\)); square and
divide. Second:
\(\|(u\cdot\nabla)u\|_2\le\|u\|_{L^6}\|\nabla u\|_{L^3}\le
C_S\|\nabla u\|_{L^2}\|\nabla u\|_{L^3}\) (Hölder \(6,3,2\); mean-zero
Sobolev \(H^1\hookrightarrow L^6\)); square and divide. ∎

**Lemma 0b (positivity; placed here because it uses Lemmas 1–4 and 7).**
\(H_0(t)>0\) for all \(t\in[0,T_{\max})\).

*Proof.* Let \(t_1=\inf\{t:H_0(t)=0\}\) and suppose \(t_1<T_{\max}\).
On \([0,t_1)\), \(H_0>0\) (and hence \(H_1>0\)), so Lemmas 1–4 and 7
apply there. Put \(M=\sup_{[0,t_1]}\|u\|_{H^m}<\infty\)
(\(u\in C([0,T_{\max});H^m)\), \([0,t_1]\) compact). By Lemma 7,
\(KD\le\|u\|_\infty^2\le C_mM^2\) on \([0,t_1)\); by Lemma 4,
\(z(t)\le z(0)+C_mM^2t_1/(2\nu)\) there, so \(N_0^2\) is bounded on
\([0,t_1)\) and \(\int_0^{t_1}N_0^2<\infty\). The energy identity
(Lemma 1) integrates to
\(H_0(t)\ge H_0(0)e^{-2\nu\int_0^{t}N_0^2}\) on \([0,t_1)\); letting
\(t\uparrow t_1\) and using continuity gives \(H_0(t_1)>0\), a
contradiction. ∎

**Lemma 8 (Osgood closure; Main (e)).** As stated in the theorem.

*Proof.* \(\Omega(z)=\int_0^zds/\Phi\) is \(C^1\), increasing,
\(\Omega(\infty)=\infty\). By Lemma 4 and the hypothesis,
\(z'\le(\Phi(z)D+R)/(2\nu)\), so
\(\Omega(z)'=z'/\Phi(z)\le D/(2\nu)+R/(2\nu\Phi(0))\) (using \(z\ge0\),
\(\Phi\) nondecreasing). Both terms are integrable on
\([0,T_{\max})\): \(\int D\le H_0(0)/(2\nu)\) (Lemma 1), \(\int R<\infty\).
Hence \(\Omega(z)\), and so \(z\), is bounded; Lemma 5's argument gives
(b1). ∎

**Assembly of the Main Theorem.** (a) = Lemma 4; (b) = Lemma 5;
(c) = Lemma 6; (d) = Lemma 7; (e) = Lemma 8. ∎

---

## Proposition (conditional static no-go)

**Lemma 9 (band symmetry identity).** Let
\(B_N=\{k\in\mathbb Z^3:1\le|k|\le N\}\) and \(f\ge0\) radial. Then
\(\sum_{k\in B_N}k_ik_jf(|k|^2)
=\delta_{ij}\tfrac13\sum_{k\in B_N}|k|^2f(|k|^2)\).

*Proof.* \(B_N\) is invariant under each sign flip
\(k_i\mapsto-k_i\) and under coordinate permutations. For \(i\neq j\)
the flip \(k_i\mapsto-k_i\) negates the summand, so the sum vanishes.
The three diagonal sums are equal by permutation invariance and add to
\(\sum|k|^2f\). ∎

**Lemma 10 (family laws).** For
\(\widehat u_N(k)=P_kv_0/|k|^2\) on \(B_N\) (real, since the
coefficients are real and \(\pm k\)-symmetric; divergence-free since
\(k\cdot P_kv_0=0\)):
\(H_0=\tfrac23\|v_0\|^2T_N\), \(H_1=\tfrac23\|v_0\|^2S_N\),
\(N_0^2=S_N/T_N\), \(u_N(0)=\tfrac23S_Nv_0\), and
\(\|u_N\|_\infty^2/H_1\ge\tfrac23S_N\).

*Proof.* \(|P_kv_0|^2=\|v_0\|^2-(k\cdot v_0)^2/|k|^2\); apply Lemma 9
with \(f=|k|^{-6}\) resp. \(|k|^{-4}\) to
\(\sum(k\cdot v_0)^2|k|^{-6}=\tfrac13\|v_0\|^2T_N\) resp.
\(\tfrac13\|v_0\|^2S_N\); subtract. Point value:
\(u_N(0)=\sum_k\widehat u_N(k)=v_0S_N-\sum k(k\cdot v_0)|k|^{-4}
=v_0S_N-\tfrac13v_0S_N\) (Lemma 9, \(f=|k|^{-4}\)). Then
\(\|u_N\|_\infty^2\ge|u_N(0)|^2=\tfrac49S_N^2\|v_0\|^2
=\tfrac23S_N\cdot H_1\). ∎

**Lemma 11 (lattice sums; crude explicit constants suffice).** For
\(N\ge8\):
\[
\tfrac1{250}\,N\;\le\;S_N\;\le\;432\,N,\qquad
T_N\nearrow T_\infty<\infty,
\]
and consequently, with \(s_N=\log N_0^2(u_N)=\log(S_N/T_N)\), the
dyadic gaps are bounded above and below for all large \(N\):
\[
0<c_-\le s_{2N}-s_N\le c_+<\infty .
\]

*Proof.* Upper: group \(B_N\) into dyadic shells
\(2^j\le|k|<2^{j+1}\); each shell has at most
\((2\cdot2^{j+1}+1)^3\le27\cdot8^{j+1}\) points (cube count), each
contributing \(\le4^{-j}\), so
\(S_N\le\sum_{2^j\le N}216\cdot2^j\le432N\). Lower: the all-positive
box \(\{N/(2\sqrt3)<k_i\le N/\sqrt3\}\) lies in the shell
\(N/2<|k|\le N\) and contains \(\ge(N/(2\sqrt3)-1)^3\ge(0.16N)^3\)
points for \(N\ge8\), each contributing \(\ge N^{-2}\):
\(S_N\ge0.004N\ge N/250\). \(T_\infty<\infty\): the same dyadic
grouping with weight \(16^{-j}\) gives a convergent geometric sum.
Gaps: \(S_{2N}\ge S_N+(\text{points in }N<|k|\le2N)\cdot(2N)^{-2}
\ge S_N+0.1N^3/(4N^2)=S_N+N/40\ge S_N\bigl(1+\tfrac1{40\cdot432}\bigr)\),
while \(S_{2N}/S_N\le432\cdot2N/(N/250)\); since
\(T_{2N}/T_N\to1\), the display follows for large \(N\) with
\(c_-=\tfrac12\log(1+\tfrac1{17280})\), \(c_+=\log(216000)+1\). ∎

*(Remark, no proof claimed: numerically \(S_N=4\pi N-8.7\pm0.9\) for
\(N\le200\) and \(T_\infty=16.5323\ldots\); the crude constants above
are all the Proposition uses.)*

**Proof of the Proposition.** Assume (L\*):
\(K(u_N)=\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2/H_1^2
\ge c_0N^3/(\tfrac23\|v_0\|^2S_N)^2\ge c_0'N\ge c_0''N_0^2
=c_0''e^{s_N}\) (Lemmas 10–11: \(S_N\le432N\) upward,
\(N_0^2=S_N/T_N\le432N/T_2\) downward-compatible). Take the dyadic
sequence \(N_j=2^jN_0\): \(s_{N_j}\to\infty\) with gaps
\(s_{N_{j+1}}-s_{N_j}\in[c_-,c_+]\) for large \(j\) (Lemma 11), so the
sequence \(\{s_{N_j}\}\) is increasing (eventually) and unbounded. If
\(K\le\Phi(\log N_0^2)\) uniformly, then
\(\Phi(s_{N_j})\ge c_0''e^{s_{N_j}}\); for
\(s\in[s_{N_j},s_{N_{j+1}}]\), monotonicity gives
\(\Phi(s)\ge\Phi(s_{N_j})\ge c_0''e^{s_{N_j}}
\ge c_0''e^{-c_+}\,e^{s}\). Hence \(\Phi(s)\ge ce^s\) for all large
\(s\) with \(c=c_0''e^{-c_+}\), and
\(\int^\infty ds/\Phi\le c^{-1}\int^\infty e^{-s}ds<\infty\). ∎

**Scope note (referee-mandated).** The Proposition excludes uniform
pointwise **field** inequalities \(K\le\Phi(\log N_0^2)\) (the
\(R\equiv0\) route to clause (e)). It does not by itself exclude
solution-adapted remainders \(R(t)\) with \(\int R<\infty\); that
route is constrained separately by the dynamic analysis recorded in
the repository's Osgood-gate note, which is not part of this paper's
claims.

**Certificates supporting (L\*)** (appendix of the LaTeX draft; not part
of any proof): exact rational values
\(K(u_4)=0.7884107043\ldots\) (full fraction stored),
\(K(u_6)=1.4344718\ldots\), \(K(u_8)=2.0372430\ldots\) (float pipeline
agreement \(\le1.3\times10^{-15}\)); the exact sweeping-pairing lower
bound \(\|\mathbb P(u\cdot\nabla u)\|_2\ge|\langle(u\cdot\nabla)u,
(v_0\cdot\nabla)u\rangle|/\|(v_0\cdot\nabla)u\|_2\) — valid because
\(\mathbb P((c\cdot\nabla)u)=(c\cdot\nabla)u\) for constant \(c\) (the
coefficient of \((c\cdot\nabla)u\) at \(k\) is \(i(c\cdot k)\hat u_k\perp
k\)) — captures a near-constant \(91.2\%\)–\(92.1\%\) of the norm for
\(4\le N\le32\) (measured range 91.18%–92.10%); \(K/N_0^2\) increases
\(0.259\to0.396\) over \(N=4\ldots32\).
