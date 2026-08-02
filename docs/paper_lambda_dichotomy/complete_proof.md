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

## Theorem O (unconditional static no-go)

*Formerly stated as a Proposition conditional on a capacity hypothesis
(L\*). The hypothesis has been discharged — for a smoothly truncated
member of the same coherent family, which is all the argument ever needed,
because it quantifies over **all** real zero-mean divergence-free
trigonometric fields. The capacity input is Theorem 7.1(3) of*
[lstar/lstar_proof_main.md](lstar/lstar_proof_main.md)*, cited but not
reproduced here. The sharply truncated family's own capacity bound — the
literal Hypothesis (L\*) — remains open and is no longer used; see the*
[audit](dependency_and_gap_audit.md) *§3 and §5.*

**The family.** Call \(\chi\) **admissible** if
\(\chi\in C^4([0,\infty);[0,1])\), \(\chi\equiv1\) on
\([0,\tfrac12]\) and \(\operatorname{supp}\chi\subset[0,1]\) (no
monotonicity assumed). *(The load-bearing minimum is \(C^1\), indeed
Lipschitz: only \(\chi'\) is used quantitatively, through
\(c_\chi=1+\|\chi'\|_{L^\infty}\) in the Riemann-sum rate. Orders
\(2,3,4\) enter only the far-field expansion
\(V=\pi^2W+O(|y|^{-4})\), which no proof consumes; the non-degeneracy
uses no derivative of \(\chi\) at all. \(C^4\) is stated so that the
exact-rational certificate cutoff below is admissible. No step
requires \(\chi\in C^\infty\).)* For \(v_0\in\mathbb R^3\setminus\{0\}\), an
admissible \(\chi\) and an integer \(N\ge2\), set
\[
\widehat u_N(k)=\chi\!\Big(\frac{|k|}N\Big)\frac{P_kv_0}{|k|^2}\ (k\neq0),
\qquad\widehat u_N(0)=0,
\]
a finite trigonometric polynomial banded at \(|k|\le N\). Write
\[
S_N^\chi=\sum_{k\neq0}\frac{\chi(|k|/N)^2}{|k|^2},\quad
T_N^\chi=\sum_{k\neq0}\frac{\chi(|k|/N)^2}{|k|^4},\quad
\Sigma_N^\chi=\sum_{k\neq0}\frac{\chi(|k|/N)}{|k|^2}.
\]

**Lemma 9 (band symmetry identity; general radial weight).** Let
\(f:\mathbb Z^3\setminus\{0\}\to\mathbb R\) be radial
(\(f(k)=\varphi(|k|^2)\)) with \(\sum_{k\neq0}|k|^2|f(k)|<\infty\). Then
\[
\sum_{k\neq0}k_ik_jf(k)=\delta_{ij}\,\tfrac13\sum_{k\neq0}|k|^2f(k).
\]

*Proof.* Absolute summability makes every rearrangement legitimate. The
index set \(\mathbb Z^3\setminus\{0\}\) is invariant under each sign flip
\(k_i\mapsto-k_i\) and under coordinate permutations, and \(f\) is
invariant under both. For \(i\neq j\) the flip \(k_i\mapsto-k_i\) is a
bijection of the index set that negates the summand, so that sum equals
its own negative and vanishes. The three diagonal sums are carried into
one another by permutations, hence are equal, and add to
\(\sum|k|^2f\). ∎

*(Earlier versions stated this only for \(f\ge0\) supported on the sharp
ball \(B_N=\{1\le|k|\le N\}\); that is the special case
\(f=\varphi\cdot\mathbb 1_{|k|\le N}\), \(\varphi\ge0\). No positivity is
needed, only absolute summability — which is what the smooth weight
\(\chi(|k|/N)^2|k|^{-6}\) requires, since \(\chi\) is not assumed
monotone and hence not assumed of one sign after differentiation.)*

**Lemma 10 (family laws, smooth truncation).** For the family above
(real, since the coefficients are real and \(\pm k\)-symmetric;
zero-mean by \(\widehat u_N(0)=0\); divergence-free since
\(k\cdot P_kv_0=0\)):
\[
H_0=\tfrac23\|v_0\|^2T_N^\chi,\quad
H_1=\tfrac23\|v_0\|^2S_N^\chi,\quad
N_0^2=\frac{S_N^\chi}{T_N^\chi},\quad
u_N(0)=\tfrac23\Sigma_N^\chi v_0,
\]
and \(\|u_N\|_\infty^2/H_1\ge2(\Sigma_N^\chi)^2/(3S_N^\chi)\).

*Proof.* \(|P_kv_0|^2=\|v_0\|^2-(k\cdot v_0)^2/|k|^2\); apply Lemma 9
with the radial weights \(f=\chi(|k|/N)^2|k|^{-6}\) resp.
\(\chi(|k|/N)^2|k|^{-4}\) (both finitely supported, hence absolutely
summable) to get
\(\sum\chi^2(k\cdot v_0)^2|k|^{-6}=\tfrac13\|v_0\|^2T_N^\chi\) resp.
\(\tfrac13\|v_0\|^2S_N^\chi\); subtract. Point value:
\(u_N(0)=\sum_k\widehat u_N(k)
=v_0\Sigma_N^\chi-\sum\chi\,k(k\cdot v_0)|k|^{-4}
=v_0\Sigma_N^\chi-\tfrac13v_0\Sigma_N^\chi\) (Lemma 9,
\(f=\chi(|k|/N)|k|^{-4}\)). Then
\(\|u_N\|_\infty^2\ge|u_N(0)|^2=\tfrac49(\Sigma_N^\chi)^2\|v_0\|^2
=\tfrac{2(\Sigma_N^\chi)^2}{3S_N^\chi}H_1\). ∎

*(This is Theorem 1.3 of the capacity document. The sharply truncated
family \(\widehat u_N(k)=P_kv_0/|k|^2\) on \(B_N\) is the formal case
\(\chi=\mathbb 1_{[0,1]}\), for which \(S_N^\chi=S_N\),
\(T_N^\chi=T_N\), \(\Sigma_N^\chi=S_N\), \(u_N(0)=\tfrac23S_Nv_0\), and
the Bernstein ratio is saturated two-sidedly,
\(\tfrac23S_N\le\|u_N\|_\infty^2/H_1\le S_N\). That family is no longer
used; \(\mathbb 1_{[0,1]}\) is not admissible.)*

**Lemma 11 (lattice sums; crude explicit constants suffice).** For
\(N\ge32\) and every admissible \(\chi\):
\[
\frac N{348}\;\le\;S_N^\chi\;\le\;128N,\qquad
6\;\le\;T_N^\chi\;\le\;T_\infty\le128,
\]
hence \(\dfrac N{44544}\le N_0^2(u_N)\le\dfrac{64}3N\); and consequently,
with \(s_N=\log N_0^2(u_N)\) and \(q=20\), the \(2^q\)-adic gaps are
bounded above and below:
\[
0<c_-=\log\frac{1048576}{950272}=0.0985\ldots\ \le\ s_{2^{20}N}-s_N\ \le\
c_+=\log\big(2^{20}\cdot950272\big)=27.63\ldots<\infty .
\]

*Proof.* Let \(S_N=\sum_{1\le|k|\le N}|k|^{-2}\),
\(T_\infty=\sum_{k\neq0}|k|^{-4}\). Dyadic shell count: for \(j\ge0\),
\(\#\{2^j\le|k|<2^{j+1}\}\le(2^{j+2}-1)^3<64\cdot8^j\), since such \(k\)
have \(|k_i|\le2^{j+1}-1\). Upper: with \(J=\lfloor\log_2N\rfloor\),
\(S_N\le\sum_{j\le J}64\cdot8^j4^{-j}=64(2^{J+1}-1)<128\cdot2^J\le128N\),
and \(T_\infty\le\sum_{j\ge0}64\cdot8^j16^{-j}=128\). Since
\(0\le\chi\le1\) with support in \(|k|\le N\), \(S_N^\chi\le S_N\le128N\)
and \(T_N^\chi\le T_\infty\le128\).
Lower on \(T_N^\chi\): \(\chi(1/N)=1\) for \(N\ge2\) and the six \(k\)
with \(|k|=1\) each contribute \(1\).
Lower on \(S_N^\chi\): \(\chi(|k|/N)=1\) whenever \(|k|\le N/2\), so
\(S_N^\chi\ge S_{\lfloor N/2\rfloor}\). For \(M\ge16\), the all-positive
box \(\{M/(2\sqrt3)<k_i\le M/\sqrt3\}\) lies in \(M^2/4<|k|^2\le M^2\),
so its points are counted in \(S_M\) and contribute \(\ge M^{-2}\) each;
each coordinate admits \(\ge M/(2\sqrt3)-1\ge0.28868M-0.0625M=0.226M\)
integers, giving \(\#\ge(0.226M)^3>0.01154M^3\) and
\(S_M\ge0.01154M\ge M/87\). For \(N\ge32\),
\(\lfloor N/2\rfloor\ge\max(16,N/4)\), so \(S_N^\chi\ge N/348\).
Combining, \(N_0^2=S_N^\chi/T_N^\chi\in[\tfrac{N/348}{128},\tfrac{128N}6]
=[\tfrac N{44544},\tfrac{64}3N]\).
Gaps: the two-sided bound gives, for every \(q\),
\(N_0^2(u_{2^qN})/N_0^2(u_N)\in[2^q/950272,\ 2^q\cdot950272]\) with
\(950272=44544\cdot64/3\); take \(q=20\), where
\(2^{20}=1048576>950272\), so the lower endpoint exceeds \(1\). ∎

*(Supersession note. Earlier versions proved \(N/250\le S_N\le432N\) for
\(N\ge8\) with dyadic gaps \(c_-=\tfrac12\log(1+\tfrac1{17280})\),
\(c_+=\log(216000)+1\), for the **sharply** truncated family. Those
constants are **not wrong** — they are merely superseded: \(128\) is the
tighter shell count, and the gap bound above is elementary and
\(N\)-uniform because \(T_N^\chi\ge6\) replaces the asymptotic
"\(T_{2N}/T_N\to1\)". They are recorded for continuity and used
nowhere.)*

*(Remark, no proof claimed, used nowhere: numerically
\(S_N=4\pi N-8.7\pm0.9\) for \(N\le200\) and
\(T_\infty=16.5323\ldots\).)*

**Capacity Theorem (cited, not reproduced).** *Theorem 7.1(3) of*
[lstar/lstar_proof_main.md](lstar/lstar_proof_main.md)*.* For every
\(v_0\in\mathbb R^3\setminus\{0\}\) and every admissible \(\chi\) there
exist \(c_0=c_0(\chi,v_0)>0\) and \(N_*=N_*(\chi,v_0)\) such that
\[
\big\|\mathbb P(u_N\cdot\nabla u_N)\big\|_{L^2(\mathbb T^3)}^2
\;\ge\;c_0N^3\qquad\text{for all }N\ge N_* .
\]
That is exactly and only what is used below. **Quantifiers.** The
bound holds *for every* admissible cutoff profile \(\chi\) and *every*
nonzero seed vector \(v_0\); the constants \(c_0=c_0(\chi,v_0)\) and
\(N_*=N_*(\chi,v_0)\) depend on that pair and are **not uniform** in
it (nor effective). What *is* uniform in \((\chi,v_0)\) is the
qualitative non-degeneracy \(\mathbb P(V\cdot\nabla V)\not\equiv0\).
For orientation, its proof
runs: \(\widehat u_N(k)=N^{-2}F(k/N)\) exactly, with
\(F(\xi)=\chi(|\xi|)P_\xi v_0/|\xi|^2\) (Theorem 3.2 there); an exact
lattice-Riemann identity \(u_N(y/N)=N(V+E_N)\),
\(\nabla u_N(y/N)=N^2(\nabla V+E'_N)\) with \(V=\widetilde F\) and
\(\sup_{|y|\le R}(|E_N|+|E'_N|)=O(N^{-1}\log N)\) (Lemma 4.1,
Theorem 4.2); pairing against a concentrated divergence-free test field
\(\psi_N(x)=N^{3/2}\Psi(Nx)\), whose \(L^2(\mathbb T^3)\) norm is
\(N\)-independent, giving
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge(2\pi)^{-3}
N^3|I_\Psi+\mathcal E_N|^2/\|\Psi\|_{L^2(\mathbb R^3)}^2\) with
\(I_\Psi=\int(V\cdot\nabla V)\cdot\Psi\) (Theorem 5.3, Corollary 5.4);
and finally the non-degeneracy \(\mathbb P(V\cdot\nabla V)\not\equiv0\),
proven for **all** admissible \(\chi\) and **all** \(v_0\) from the
\(\chi\)-independent leading singularity of \(F*F\) at \(\zeta=0\),
where \(\zeta\times(\tau(\zeta)\zeta)
=\tfrac{3\pi^3}8\|v_0\|^2\tfrac{\zeta_3}{|\zeta|}(\zeta\times e_3)\neq0\)
— the Fourier form of "the Oseen field is not a stationary Euler flow"
(Theorems 6.5, 6.7). The exponent \(3\) is sharp: Lemma 7 with Lemma 10
caps \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\le\|u_N\|_\infty^2H_1
\asymp N^3\).
**The constants \(c_0\) and \(N_*\) are non-effective**, because the test
field \(\Psi\) is obtained from the density of
\(C^\infty_{c,\sigma}(\mathbb R^3)\) in \(L^2_\sigma(\mathbb R^3)\); see
"Constant dependence" in [theorem_statement.md](theorem_statement.md).

**Proof of Theorem O.** Fix any admissible \(\chi\) and any
\(v_0\in\mathbb R^3\setminus\{0\}\), and use the family \(u_N\); each
\(u_N\) is a real, zero-mean, divergence-free trigonometric field
(Lemma 10), hence an admissible test object for the hypothesis on
\(\Phi\).
By Lemma 10 and Lemma 11, \(H_1=\tfrac23\|v_0\|^2S_N^\chi
\le\tfrac{256}3\|v_0\|^2N\) for \(N\ge32\), so by the Capacity Theorem,
for \(N\ge\max(N_*,32)\),
\[
K(u_N)=\frac{\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2}{H_1^2}
\ \ge\ \frac{c_0N^3}{\big(\tfrac{256}3\|v_0\|^2N\big)^2}
=\frac{9c_0}{65536\|v_0\|^4}\,N .
\]
By Lemma 11 again, \(N\ge\tfrac3{64}N_0^2(u_N)\), hence
\(K(u_N)\ge c_0''\,N_0^2(u_N)=c_0''e^{s_N}\) with
\(c_0''=\tfrac{27c_0}{4194304\|v_0\|^4}\).
Take \(N_j=2^{20j}N_\sharp\) with \(N_\sharp\ge\max(N_*,32)\): by
Lemma 11 the gaps satisfy \(s_{N_{j+1}}-s_{N_j}\in[c_-,c_+]\), so
\(\{s_{N_j}\}\) is strictly increasing and unbounded. If
\(K\le\Phi(\log N_0^2)\) uniformly, then
\(\Phi(s_{N_j})\ge K(u_{N_j})\ge c_0''e^{s_{N_j}}\); for
\(s\in[s_{N_j},s_{N_{j+1}}]\), monotonicity gives
\(\Phi(s)\ge\Phi(s_{N_j})\ge c_0''e^{s_{N_j}}\ge c_0''e^{-c_+}e^{s}\).
Hence \(\Phi(s)\ge ce^s\) for all \(s\ge s_{N_0}\) with
\(c=c_0''e^{-c_+}\), and
\(\int^\infty ds/\Phi\le c^{-1}\int^\infty e^{-s}ds<\infty\). ∎

*(No hypothesis was assumed. The argument is verbatim the one that
formerly consumed (L\*), with the sharp family replaced by the smooth one
and (L\*) replaced by the proven Capacity Theorem. \(c\) inherits the
non-effectivity of \(c_0\).)*

**Scope note (referee-mandated).** Theorem O excludes uniform
pointwise **field** inequalities \(K\le\Phi(\log N_0^2)\) (the
\(R\equiv0\) route to clause (e)). It does not by itself exclude
solution-adapted remainders \(R(t)\) with \(\int R<\infty\); that
route is constrained separately by the dynamic analysis recorded in
the repository's Osgood-gate note, which is not part of this paper's
claims.

## Certificate appendix (numerics; part of no proof)

Nothing in this section enters any proof. Both the Main Theorem and
Theorem O are now unconditional (modulo F1–F3), so — unlike in earlier
versions — these numbers no longer carry a hypothesis. They are recorded
as corroboration.

### C.1 Smooth family (the family actually used by Theorem O)

Producer: `experiments/run_smooth_family_capacity.py`; data
`outputs/lstar/smooth_family_capacity.json`
(sha256 `527b3d0b…d4b3`); tables reproduced in
[lstar/lstar_numerical_support.md](lstar/lstar_numerical_support.md).
Seed \(v_0=(1,2,3)\), \(\|v_0\|^2=14\); weight
\(\chi(r)=1-S\big(\mathrm{clamp}(\tfrac{r^2-1/4}{3/4},0,1)\big)\) with
\(S\) the degree-9 smoothstep
\(126s^5-420s^6+540s^7-315s^8+70s^9\).

This \(\chi\) is a polynomial in \(r^2\) — which is exactly what makes
the whole family **exactly rational** on the lattice — and lies in
\(C^4\setminus C^5\), with \(0\le\chi\le1\), \(\chi\equiv1\) on
\([0,\tfrac12]\), \(\operatorname{supp}\chi=[0,1]\) and
\(\|\chi'\|_\infty=5.2167\ldots\). It is therefore **admissible** as
defined above, so the exact-rational lanes certify the family the
Capacity Theorem literally covers. Degree-\((2p+1)\) smoothsteps give
rational \(\chi\in C^p\) for every \(p\), should a stronger class ever
be wanted.

Exact moment laws (`fractions.Fraction`, no rounding): the residuals of
\(H_0=\tfrac23\|v_0\|^2T_N^\chi\), \(H_1=\tfrac23\|v_0\|^2S_N^\chi\) and
\(u_N(0)=\tfrac23\Sigma_N^\chi v_0\) are the **rational \(0\)** — not
"small" — at \(N=4,6,8,10,12\), with \(u_N(0)\) exactly parallel to
\(v_0\) (components in ratio \(1:2:3\)).

Exact capacity, full \(O(|B_N|^2)\) contraction over all ordered pairs of
active band points (no truncation, no sampling):

| \(N\) | \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2/N^3\) | \(K\) | \(N_0^2\) | \(K/N_0^2\) |
|---:|---:|---:|---:|---:|
| 4  |  678.284307 | 0.61657681 | 2.316030 | 0.266221 |
| 6  | 1027.307364 | 1.13589428 | 3.455252 | 0.328744 |
| 8  | 1241.509295 | 1.67053905 | 4.585959 | 0.364273 |
| 10 | 1384.356255 | 2.20774755 | 5.718747 | 0.386054 |
| 12 | 1485.779047 | 2.74636641 | 6.852077 | 0.400808 |
| 16 | 1619.652353 | 3.82541413 | 9.119441 | 0.419479 |

Dealiased FFT continuation (grid \(4N+2\); the whole product spectrum is
represented) extends this to \(N=48\):
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2/N^3\) rises
\(1241.51\to1911.31\) and \(K/N_0^2\) rises \(0.364273\to0.457376\) over
\(N=8\to48\), both monotonically. **Exact vs. FFT** at the two shared
bands agrees to double-precision roundoff: relative difference
\(5.5\times10^{-16}\) (\(N=8\)) and \(7.0\times10^{-16}\) (\(N=16\)) in
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\); \(5.3\times10^{-16}\) and
\(8.1\times10^{-16}\) in \(K\).

Two-point Richardson limits (diagnostics, not proofs, fitted on
\(N=40,48\)): \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2/N^3\to2064.31\),
\(K/N_0^2\to0.47644\), \(H_1/N\to87.51\), \(N_0^2/N\to0.56708\). Both
target sequences increase toward positive constants, i.e. the **sharp**
exponent \(a=1\), consistent with the Capacity Theorem's \(N^3\).
The continuum non-degeneracy behind that theorem is corroborated three
independent ways (exact rational lattice modes; float 3-D quadrature of
the continuum convolution; an axisymmetric radial reduction giving
\(\mathrm{curl}(V\cdot\nabla V)=Z(r)\sin\theta\cos\theta\,e_\phi\) with
\(Z\to-12\pi^4\|v_0\|^2/r^4\)), agreeing to \(10^{-4}\) or better.

### C.2 Sharp family (corroboration only; the family is no longer used)

*These are **sharp**-family numbers. They corroborate the literal
Hypothesis (L\*), which remains **open** and which no statement in this
paper now uses. They are retained, not deleted, because they are
independently verified and were the original evidence for the
obstruction; they support nothing that is claimed here.*

Exact rational values \(K(u_4)=0.7884107043\ldots\) (full fraction
stored), \(K(u_6)=1.4344718\ldots\), \(K(u_8)=2.0372430\ldots\) (float
pipeline agreement \(\le1.3\times10^{-15}\)); the exact sweeping-pairing
lower bound \(\|\mathbb P(u\cdot\nabla u)\|_2\ge|\langle(u\cdot\nabla)u,
(v_0\cdot\nabla)u\rangle|/\|(v_0\cdot\nabla)u\|_2\) — valid because
\(\mathbb P((c\cdot\nabla)u)=(c\cdot\nabla)u\) for constant \(c\) (the
coefficient of \((c\cdot\nabla)u\) at \(k\) is \(i(c\cdot k)\hat u_k\perp
k\)) — captures a near-constant \(91.2\%\)–\(92.1\%\) of the norm for
\(4\le N\le32\) (measured range 91.18%–92.10%); \(K/N_0^2\) increases
\(0.259\to0.396\) over \(N=4\ldots32\). The smooth family's \(K/N_0^2\)
runs slightly above these at every comparable band.

A **proven negative result** about this family is recorded in
[lstar/lstar_direct_route_and_weakening.md](lstar/lstar_direct_route_and_weakening.md)
§A: no constant-vector sweeping split, in any Hölder pairing and any
constant pairing direction, can reach the sharp family's capacity bound —
every such estimate is capped at \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2
\gtrsim N\), two powers short. That is why the concentrated test field of
the Capacity Theorem, not a constant sweep, is the right instrument.
