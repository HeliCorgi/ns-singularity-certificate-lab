# Hypothesis (L\*) for the smoothly truncated coherent family — complete proof

**Status of the file.** This is the stand-alone proof document for the capacity
bound. It has been **adversarially refereed**
([lstar_referee.md](lstar_referee.md), verdict (i); no CRITICAL items; the
referee reports failing to break the route at all six points of its brief, and
confirms the pivotal identity Theorem 6.5 to 14 digits by an independent method
sharing no code, no convention and no regularisation with the derivation).
It has since been **integrated** into the paper files
([theorem_statement.md](../theorem_statement.md),
[complete_proof.md](../complete_proof.md),
[dependency_and_gap_audit.md](../dependency_and_gap_audit.md),
[paper_draft.tex](../paper_draft.tex)), which cite Theorem 7.1(3) below rather
than reproducing it.
Every statement below is labelled **PROVEN**, **NUMERICAL** or **OPEN**.
Numerical statements say whether they are float, exact-rational, or interval.

**Referee-mandated deletion.** The former Proposition 3.4 / equation (3.2)
(an Abel-regularised Poisson identity) is **FALSE** and has been **deleted**;
§3.4 now records the periodisation no-go instead. The statement was inert —
nothing cited it — so nothing else in this document changed.

**Non-effectivity.** The constants \(c_0\) and \(N_*\) of Theorem 7.1(3) are
**not effective**; see the effectivity caveat there. Every other constant in
this document is explicit.

**Headline.** Hypothesis (L\*) — the last open bridge in the Proposition of
`theorem_statement.md` — is **PROVEN** here, with the *sharp* exponent
\(N^3\), for a smoothly truncated variant of the coherent family. The
Proposition therefore becomes **unconditional**. The originally stated (L\*),
for the *sharply* truncated family, remains **OPEN**, but is no longer needed
by anything.

**One correction to the proposed route.** The reduction sketch supplied to this
task (item R3) asserts an *absolutely convergent* Poisson-summation identity
\(u_N(x)=N\sum_{m}V(N(x+2\pi m))\) with \(m\neq0\) terms \(O(1/N)\) uniformly.
That assertion is **false**: §3.4 proves that the \(m\neq0\) series **diverges**,
because \(V\) decays exactly like \(|y|^{-1}\) with a *non-oscillating,
sign-definite* leading profile. **The naive absolutely convergent periodisation
argument is therefore unavailable**; the divergent terms are moreover eventually
positive, so Abel summation in particular does not converge either (Remark 3.4).
No broader claim about summation methods is made or needed. The proof below therefore
does **not** use Poisson summation:
§4 replaces it with an exact lattice-Riemann-sum identity, which is elementary,
quantitative, and gives the \(C^1_{\rm loc}\) statement (including derivatives)
with an explicit rate. Everything downstream is unaffected.

---

## 0. Conventions, and verification of the reduction (R1)

### 0.1 Conventions (repo standard)

\(\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3\) with the **normalised** measure
\(d\mu=(2\pi)^{-3}dx\) (`fourier_torus` module docstring). For
\(u:\mathbb T^3\to\mathbb R^3\),

\[
u(x)=\sum_{k\in\mathbb Z^3}\hat u(k)e^{ik\cdot x},\qquad
\hat u(k)=\int_{\mathbb T^3}u(x)e^{-ik\cdot x}\,d\mu(x),
\]
\[
\langle f,g\rangle=\int_{\mathbb T^3}f\cdot\bar g\,d\mu=\sum_k\hat f(k)\cdot\overline{\hat g(k)},
\qquad
H_r=\sum_{k}|k|^{2r}|\hat u(k)|^2 .
\]

On \(\mathbb R^3\) we use throughout the **inverse-transform bracket**
\[
\widetilde G(y):=\int_{\mathbb R^3}G(\xi)\,e^{i\xi\cdot y}\,d\xi ,
\qquad\text{inverse: } G(\xi)=(2\pi)^{-3}\!\int_{\mathbb R^3}\widetilde G(y)e^{-i\xi\cdot y}dy,
\]
and the forward transform \(\mathcal F f(\zeta)=\int f(y)e^{-i\zeta\cdot y}dy\),
so that \(\mathcal F\widetilde G=(2\pi)^3G\) and Plancherel reads
\(\|\widetilde G\|_{L^2(\mathbb R^3)}^2=(2\pi)^3\|G\|_{L^2(\mathbb R^3)}^2\).
\(\mathbb P\) is the Leray projection: on \(\mathbb T^3\),
\(\widehat{\mathbb Pf}(k)=P_k\hat f(k)\), \(P_k=I-\frac{k\otimes k}{|k|^2}\)
(\(k\neq0\); \(k=0\) untouched, as in `fourier_torus.leray`); on \(\mathbb R^3\),
the same symbol \(P_\zeta\).

\(P_\xi\) is homogeneous of degree \(0\), smooth on \(\mathbb R^3\setminus\{0\}\),
symmetric, \(P_\xi\xi=0\), \(P_{-\xi}=P_\xi\), \(\|P_\xi\|\le1\).

### 0.2 The cutoff class

**Definition 0.1 (admissible \(\chi\)).** \(\chi\in C^4([0,\infty);[0,1])\)
with \(\chi\equiv1\) on \([0,\tfrac12]\) and \(\operatorname{supp}\chi\subset[0,1]\).
No monotonicity is assumed (it is used nowhere). Set
\(c_\chi:=1+\|\chi'\|_{L^\infty}\).

*Where the four derivatives go, and why there is no \(C^\infty\).* Only
\(\chi'\) is ever used **quantitatively**: it enters this document through the
single constant \(c_\chi\), in the Riemann-sum rate of Theorem 4.2, and
\(\chi\in C^1\) (indeed \(\chi\) merely Lipschitz) already suffices for every
statement on which Theorem 7.1 depends. Derivatives of order \(2,3,4\) are used
at exactly one place — Theorem 3.3(b), to give the far-field remainder
\(\rho=V-\pi^2W\) the decay \(O(|y|^{-4})\) — and that decay is consumed only by
the inert §3.4. Sections 1, 2, 3.2, 5 and 6 use **no** derivative of \(\chi\) at
all; in particular \(V\in C^\infty\) and the Paley–Wiener statement of
Theorem 3.2(3) are bought by the *compact support of \(F\)*, not by smoothness of
\(\chi\), and the nonvanishing Theorem 6.7 uses only \(0\le\chi\le1\) and
\(\chi\equiv1\) on \([0,\tfrac12]\). Nothing in this document requires
\(\chi\in C^\infty\).

Such \(\chi\) exist in abundance, and both cutoffs used in the repository are
admissible. The float lanes (`check_lstar.py`) use the \(C^\infty\) profile
\(\chi(t)=\Theta\!\big(\tfrac{t-1/2}{1/2}\big)\) with
\(\Theta(s)=\frac{e^{-1/(1-s)}}{e^{-1/s}+e^{-1/(1-s)}}\) on \((0,1)\).
The exact-rational lanes (`run_smooth_family_capacity.py`) use
\(\chi(r)=1-S\!\big(\tfrac{r^2-1/4}{3/4}\big)\), \(S\) the degree-9 smoothstep
\(126s^5-420s^6+540s^7-315s^8+70s^9\) clamped to \([0,1]\): a polynomial in
\(r^2\), hence \(\chi(|k|/N)\in\mathbb Q\) at every lattice point, with
\(\chi\in C^4\setminus C^5\) and \(\|\chi'\|_{L^\infty}=5.2167\ldots\), i.e.
\(c_\chi=6.2168\ldots\).

### 0.3 Verification of reduction (R1)

**Claim (R1).** If \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge cN^{2+a}\) for some
\(a>0\) along a sequence \(N\to\infty\), and \(H_1\asymp N\), \(N_0^2\asymp N\)
with two-sided constants, then no Osgood-admissible \(\Phi\) can satisfy the
uniform field inequality \(K\le\Phi(\log N_0^2)\).

**PROVEN.** \(K=\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2/H_1^2\ge
cN^{2+a}/(C_1N)^2=c'N^{a}\). With \(N_0^2\le C_2N\) we get \(N\ge N_0^2/C_2\),
so \(K\ge c''e^{as}\), \(s=\log N_0^2\). Along a subsequence with
\(s_j\to\infty\) and gaps \(s_{j+1}-s_j\le c_+\), monotonicity of \(\Phi\)
gives, for \(s\in[s_j,s_{j+1}]\),
\(\Phi(s)\ge\Phi(s_j)\ge c''e^{as_j}\ge c''e^{-ac_+}e^{as}\); hence
\(\int^\infty ds/\Phi\le\frac{e^{ac_+}}{c''}\int^\infty e^{-as}ds<\infty\). ∎

So any exponent \(2+a\), \(a>0\), suffices. **We will nevertheless obtain the
sharp \(a=1\)**, so R1 is a safety net that is not called upon. The Proposition's
proof in `complete_proof.md` is reproduced verbatim with \(a=1\) in §7.

---

## 1. The smoothly truncated family and its exact moment laws

### 1.1 Definition

**Definition 1.1.** Fix \(v_0\in\mathbb R^3\setminus\{0\}\), an admissible
\(\chi\), and an integer \(N\ge2\). Define \(u_N=u_N^{\chi,v_0}\) by
\[
\boxed{\ \widehat{u_N}(k)=\chi\!\Big(\frac{|k|}{N}\Big)\frac{P_kv_0}{|k|^2}\quad(k\neq0),
\qquad \widehat{u_N}(0)=0 .\ }
\]

Note \(\widehat{u_N}(k)=0\) for \(|k|>N\), so \(u_N\) is a **finite
trigonometric polynomial** — an admissible test field for the Proposition's
hypothesis, which quantifies over all real zero-mean divergence-free
trigonometric fields. (Real, not rational, coefficients; the Proposition's
hypothesis has no rationality requirement.)

**Lemma 1.2 (basic structure). PROVEN.** \(u_N\) is real, zero-mean, exactly
divergence-free, \(C^\infty\), and \(u_N\neq0\).

*Proof.* Zero mean: \(\widehat{u_N}(0)=0\). Real: the coefficients are real
vectors and \(\widehat{u_N}(-k)=\chi(|k|/N)P_{-k}v_0/|k|^2=\widehat{u_N}(k)
=\overline{\widehat{u_N}(k)}\), which is the reality condition
\(\hat u(-k)=\overline{\hat u(k)}\). Divergence-free: \(k\cdot P_kv_0=0\) for
every \(k\neq0\), so \(\widehat{\nabla\cdot u_N}(k)=ik\cdot\widehat{u_N}(k)=0\).
Smooth: finite sum. Nonzero: pick \(k\) with \(|k|=1\) not parallel to \(v_0\)
(possible since \(v_0\neq0\) — among \(e_1,e_2,e_3\) at least two are not
parallel to \(v_0\)); then \(\chi(1/N)=1\) (as \(1/N\le1/2\)) and
\(P_kv_0\neq0\). ∎

### 1.2 Band symmetry for a general radial weight

The existing Lemma 9 of `complete_proof.md` is stated for the sharp ball
\(B_N\) and \(f\ge0\) radial. The generalisation we need is:

**Lemma 9′ (band symmetry, general radial weight). PROVEN.** Let
\(f:\mathbb Z^3\setminus\{0\}\to\mathbb R\) be radial (i.e.\ \(f(k)=\varphi(|k|^2)\)
for some \(\varphi\)) and absolutely summable. Then
\[
\sum_{k\neq0}k_ik_j\,f(k)=\delta_{ij}\,\tfrac13\sum_{k\neq0}|k|^2f(k),
\]
provided the right side converges absolutely.

*Proof.* Absolute summability of \(|k|^2f\) makes every rearrangement legitimate.
The index set \(\mathbb Z^3\setminus\{0\}\) is invariant under each sign flip
\(\sigma_i:k_i\mapsto-k_i\) and under coordinate permutations, and \(f\) is
invariant under both (it depends on \(|k|^2\) only). For \(i\neq j\), the map
\(\sigma_i\) is a bijection of the index set that negates the summand
\(k_ik_jf(k)\); hence that sum equals its own negative, so vanishes. For \(i=j\),
the three sums \(\sum k_i^2f\) are carried into one another by coordinate
permutations, hence are equal; they add to \(\sum|k|^2f\). ∎

*(No positivity of \(f\) is needed — only absolute summability. Lemma 9 of the
paper is the special case \(\varphi=f\cdot\mathbb 1_{|k|\le N}\).)*

**Lemma 9″ (continuous band symmetry). PROVEN.** For \(g:\mathbb R^3\to\mathbb R\)
radial with \(\int|\xi|^2|g|\,d\xi<\infty\) locally as needed,
\(\int\xi_i\xi_j g(\xi)d\xi=\delta_{ij}\frac13\int|\xi|^2g\,d\xi\).
*Proof.* Identical: \(d\xi\) and \(g\) are invariant under sign flips and
permutations of coordinates. ∎

### 1.3 The exact moment laws

**Notation.**
\[
S_N^\chi:=\sum_{k\neq0}\frac{\chi(|k|/N)^2}{|k|^{2}},\qquad
T_N^\chi:=\sum_{k\neq0}\frac{\chi(|k|/N)^2}{|k|^{4}},\qquad
\Sigma_N^\chi:=\sum_{k\neq0}\frac{\chi(|k|/N)}{|k|^{2}} .
\]
All three are finite sums (support \(|k|\le N\)).

**Theorem 1.3 (family laws). PROVEN — exact identities.**
\[
H_0(u_N)=\tfrac23\|v_0\|^2\,T_N^\chi,\qquad
H_1(u_N)=\tfrac23\|v_0\|^2\,S_N^\chi,\qquad
N_0^2(u_N)=\frac{S_N^\chi}{T_N^\chi},\qquad
u_N(0)=\tfrac23\,\Sigma_N^\chi\,v_0 .
\]

*Proof.* Since \(P_k\) is an orthogonal projection,
\(|P_kv_0|^2=\|v_0\|^2-\frac{(k\cdot v_0)^2}{|k|^2}\). Hence
\[
H_0=\sum_{k\neq0}\frac{\chi(|k|/N)^2}{|k|^4}|P_kv_0|^2
=\|v_0\|^2\sum\frac{\chi^2}{|k|^4}-\sum\frac{\chi^2(k\cdot v_0)^2}{|k|^6}.
\]
Apply Lemma 9′ with the radial weight \(f(k)=\chi(|k|/N)^2|k|^{-6}\):
\(\sum \chi^2(k\cdot v_0)^2|k|^{-6}=v_{0i}v_{0j}\sum k_ik_j f
=\frac13\|v_0\|^2\sum\chi^2|k|^{-4}=\frac13\|v_0\|^2T_N^\chi\).
Therefore \(H_0=\|v_0\|^2T_N^\chi-\frac13\|v_0\|^2T_N^\chi=\frac23\|v_0\|^2T_N^\chi\).
Identically with \(f=\chi^2|k|^{-4}\) for \(H_1=\sum|k|^2|\widehat{u_N}(k)|^2\):
\(H_1=\|v_0\|^2S_N^\chi-\frac13\|v_0\|^2S_N^\chi=\frac23\|v_0\|^2S_N^\chi\).
Point value:
\(u_N(0)=\sum_k\widehat{u_N}(k)
=v_0\Sigma_N^\chi-\sum\frac{\chi\,k(k\cdot v_0)}{|k|^4}
=v_0\Sigma_N^\chi-\tfrac13v_0\Sigma_N^\chi\) by Lemma 9′ with
\(f=\chi|k|^{-4}\). ∎

So **every exact law of Lemma 10 of the paper survives verbatim** with
\(S_N\to S_N^\chi\), \(T_N\to T_N^\chi\), \(\Sigma_N^\chi\) replacing \(S_N\) in
the point value. In particular the Bernstein saturation
\(\|u_N\|_\infty^2/H_1\ge\frac{2(\Sigma_N^\chi)^2}{3S_N^\chi}\) holds; it is not
used below.

**NUMERICAL (float, `check_lstar.py` §B).** For \(v_0=(1,2,-1)\) and the
specific \(\chi\) above, the four identities of Theorem 1.3 were verified by
direct summation on a full lattice at \(N=4,8,16\) to \(10\) decimal places
(agreement in every printed digit).

---

## 2. Two-sided lattice bounds: \(N_0^2\asymp N\)

Crude explicit constants suffice; we make no attempt at sharpness.

**Lemma 2.1 (dyadic shell count). PROVEN.** For \(j\ge0\),
\(\#\{k\in\mathbb Z^3:2^j\le|k|<2^{j+1}\}\le 64\cdot8^{j}\).

*Proof.* Such \(k\) satisfy \(|k|_\infty\le|k|<2^{j+1}\), hence
\(|k_i|\le2^{j+1}-1\); the count is at most
\((2^{j+2}-1)^3<2^{3j+6}=64\cdot 8^j\). ∎

**Lemma 2.2 (upper bounds). PROVEN.** With
\(S_N:=\sum_{1\le|k|\le N}|k|^{-2}\) and \(T_\infty:=\sum_{k\neq0}|k|^{-4}\):
\[
S_N\le128N\ (N\ge1),\qquad T_\infty\le128 .
\]
*Proof.* Let \(J=\lfloor\log_2N\rfloor\). Every \(k\) with \(1\le|k|\le N\) lies
in a shell \(j\le J\), where \(|k|^{-2}\le4^{-j}\). By Lemma 2.1,
\(S_N\le\sum_{j=0}^{J}64\cdot8^j4^{-j}=64(2^{J+1}-1)<128\cdot2^{J}\le128N\).
Similarly \(T_\infty\le\sum_{j\ge0}64\cdot8^j16^{-j}=64\cdot2=128\). ∎

**Lemma 2.3 (lower bound). PROVEN.** For \(M\ge16\), \(S_M\ge M/87\).

*Proof.* Consider \(\mathcal B=\{k\in\mathbb Z^3:\ \tfrac{M}{2\sqrt3}<k_i\le\tfrac M{\sqrt3},\ i=1,2,3\}\).
For \(k\in\mathcal B\): \(|k|^2\le3(M/\sqrt3)^2=M^2\) and
\(|k|^2>3M^2/12=M^2/4\), so \(k\) is counted in \(S_M\) and contributes
\(|k|^{-2}\ge M^{-2}\). Each coordinate ranges over a half-open interval of
length \(M/(2\sqrt3)\), containing at least \(\frac{M}{2\sqrt3}-1\) integers.
For \(M\ge16\), \(\frac{M}{2\sqrt3}-1\ge0.2887M-0.0625M=0.226M\). Hence
\(\#\mathcal B\ge(0.226M)^3>0.01154M^3\) and
\(S_M\ge0.01154M^3\cdot M^{-2}=0.01154M\ge M/87\). ∎

**Proposition 2.4 (two-sided \(N_0^2\asymp N\)). PROVEN.** For \(N\ge32\),
\[
\frac{N}{348}\ \le\ S_N^\chi\ \le\ 128N,\qquad
6\ \le\ T_N^\chi\ \le\ 128,\qquad
\frac{N}{44544}\ \le\ N_0^2(u_N)\ \le\ \frac{64}{3}N .
\]

*Proof.* Upper bounds on \(S_N^\chi,T_N^\chi\): \(0\le\chi\le1\) and
\(\operatorname{supp}\chi(|\cdot|/N)\subset\{|k|\le N\}\), so
\(S_N^\chi\le S_N\le128N\) and \(T_N^\chi\le T_\infty\le128\).
Lower bound on \(T_N^\chi\): \(\chi(1/N)=1\) for \(N\ge2\) and the six
\(k\) with \(|k|=1\) each contribute \(1\), so \(T_N^\chi\ge6\).
Lower bound on \(S_N^\chi\): \(\chi(|k|/N)=1\) whenever \(|k|\le N/2\); hence
\(S_N^\chi\ge S_{\lfloor N/2\rfloor}\). For \(N\ge32\) we have
\(\lfloor N/2\rfloor\ge16\) and \(\lfloor N/2\rfloor\ge N/2-1\ge N/4\), so by
Lemma 2.3 \(S_N^\chi\ge\frac{N/4}{87}=\frac N{348}\).
Combining, \(N_0^2=S_N^\chi/T_N^\chi\in[\frac{N/348}{128},\frac{128N}{6}]\). ∎

**NUMERICAL (float).** For the specific \(\chi\) used, \(N_0^2/N\) equals
\(0.5501,\ 0.5405,\ 0.5372\) at \(N=4,8,16\) — comfortably inside
\([1/44544,\,64/3]\).

---

## 3. The profile \(V\): definition, properties, and the periodisation no-go

### 3.1 Definition of \(F\) and \(V\)

**Definition 3.1.** On \(\mathbb R^3\) set
\[
h(\xi):=\frac{P_\xi v_0}{|\xi|^2}\ \ (\xi\neq0),\qquad
F(\xi):=\chi(|\xi|)\,h(\xi)\ \ (\xi\neq0),\quad F(0):=0,
\]
\[
V(y):=\widetilde F(y)=\int_{\mathbb R^3}F(\xi)e^{i\xi\cdot y}\,d\xi .
\]

**Scaling identity. PROVEN.** \(P_\xi\) is homogeneous of degree \(0\), so
\(F(k/N)=\chi(|k|/N)P_kv_0\,N^2/|k|^2\), i.e.
\[
\boxed{\ \widehat{u_N}(k)=N^{-2}F(k/N)\quad\text{for all }k\in\mathbb Z^3.\ }
\tag{3.1}
\]
(The \(k=0\) case holds by the conventions \(\widehat{u_N}(0)=F(0)=0\).)

### 3.2 Elementary properties of \(V\)

**Theorem 3.2. PROVEN.**
1. \(|F(\xi)|\le\|v_0\|\,|\xi|^{-2}\mathbb 1_{|\xi|\le1}\); hence
   \(F\in L^1\cap L^p(\mathbb R^3)\) for every \(p<3/2\), with
   \(\|F\|_{L^1}\le4\pi\|v_0\|\).
2. \(V\) is real, even, bounded and continuous, with
   \(\|V\|_\infty\le4\pi\|v_0\|\).
3. \(V\in C^\infty(\mathbb R^3)\) and for every multi-index \(\beta\),
   \(\|\partial^\beta V\|_\infty\le\|v_0\|\int_{|\xi|\le1}|\xi|^{|\beta|-2}d\xi
   =\frac{4\pi\|v_0\|}{|\beta|+1}\); in particular
   \(\|\nabla V\|_\infty\le2\pi\|v_0\|\).
   Moreover \(V\) extends to an entire function on \(\mathbb C^3\) of
   exponential type \(\le1\) (Paley–Wiener).
   (Both statements use only that \(F\in L^1\) is supported in
   \(\{|\xi|\le1\}\); no derivative of \(\chi\) is involved, so they are
   unaffected by the regularity assumed in Definition 0.1.)
4. \(\nabla\cdot V\equiv0\).
5. \(V(0)=\frac{8\pi}{3}\Big(\int_0^\infty\chi(r)\,dr\Big)v_0\neq0\).

*Proof.* (1) \(\|P_\xi v_0\|\le\|v_0\|\) and \(0\le\chi\le1\) supported in
\([0,1]\); \(\int_{|\xi|\le1}|\xi|^{-2}d\xi=4\pi\); \(\int_{|\xi|\le1}|\xi|^{-2p}d\xi<\infty\)
iff \(2p<3\).
(2) \(F\in L^1\) gives continuity and \(\|V\|_\infty\le\|F\|_{L^1}\)
(dominated convergence). \(F\) is real and even (\(P_{-\xi}=P_\xi\)), so
\(V(y)=\int F(\xi)\cos(\xi\cdot y)d\xi\) is real and even.
(3) \(\xi^\beta F\in L^1\) for every \(\beta\) by (1) (the exponent
\(|\beta|-2>-3\)); differentiation under the integral sign is justified by
dominated convergence, giving
\(\partial^\beta V(y)=\int(i\xi)^\beta F(\xi)e^{i\xi y}d\xi\) and the stated
bound. Since \(F\in L^1\) has support in \(\{|\xi|\le1\}\), the integral
converges for complex \(y\) with \(|{\rm Im}\,y|\) arbitrary and obeys
\(|V(y)|\le\|F\|_{L^1}e^{|{\rm Im}\,y|}\): entire of exponential type \(\le1\).
(4) \(\widehat{\nabla\cdot V}\)-symbol is \(i\xi\cdot F(\xi)=\chi(|\xi|)\xi\cdot P_\xi v_0/|\xi|^2=0\);
concretely \(\nabla\cdot V(y)=\int i\,\xi\cdot F(\xi)e^{i\xi y}d\xi=0\).
(5) \(V(0)=\int F=\ \int\chi(|\xi|)\frac{v_0}{|\xi|^2}d\xi-\int\frac{\chi(|\xi|)\xi(\xi\cdot v_0)}{|\xi|^4}d\xi\).
Lemma 9″ with the radial \(g=\chi(|\xi|)|\xi|^{-4}\) gives
\(\int\xi_i\xi_j g=\frac{\delta_{ij}}3\int|\xi|^2g=\frac{\delta_{ij}}3\int\chi|\xi|^{-2}\),
so the second integral is \(\frac13v_0\int\chi|\xi|^{-2}d\xi\). Hence
\(V(0)=\frac23v_0\int_{\mathbb R^3}\frac{\chi(|\xi|)}{|\xi|^2}d\xi
=\frac23v_0\cdot4\pi\int_0^\infty\chi(r)dr\). It is nonzero because
\(\chi\ge0\) and \(\chi\equiv1\) on \([0,\frac12]\), so
\(\int_0^\infty\chi\ge\frac12\). ∎

*(Structural remark, PROVEN: \(F=\mathbb P_\xi(v_0\,\chi(|\xi|)|\xi|^{-2})\), so
\(V=\mathbb P(v_0\,g_\chi)\) with \(g_\chi\) the smoothed Newtonian potential
\(\widetilde{\chi(|\cdot|)|\cdot|^{-2}}\). Writing \(\Delta\phi=g_\chi\) with
\(\phi\) radial gives the Oseen-type normal form
\(V(y)=a(r)v_0+b(r)(v_0\cdot\hat y)\hat y\), \(a=\phi''+\phi'/r\),
\(b=\phi'/r-\phi''\). We do not need this form, but it explains §3.3.)*

### 3.3 The far field: exact \(|y|^{-1}\) profile

**Theorem 3.3 (far-field asymptotics). PROVEN.** Let
\[
W(y):=\frac{v_0}{|y|}+\frac{(v_0\cdot y)\,y}{|y|^3}\qquad(y\neq0).
\]
Then \(\widetilde h=\pi^2W\) as tempered distributions, and
\[
V(y)=\pi^2W(y)+\rho(y),\qquad
|\rho(y)|\le C_M(\chi,v_0)\,|y|^{-M}\ \ \text{for }|y|\ge1,\ \ 2\le M\le4 .
\]

*Proof.* **(a) \(\widetilde h=\pi^2W\).** The Oseen tensor
\(\mathcal O_{ij}(y)=\frac1{8\pi}\big(\frac{\delta_{ij}}{|y|}+\frac{y_iy_j}{|y|^3}\big)\)
satisfies \(\mathcal F[\mathcal O_{ij}](\xi)=\frac{\delta_{ij}}{|\xi|^2}-\frac{\xi_i\xi_j}{|\xi|^4}=h_{ij}(\xi)\)
(where \(h_i=h_{ij}v_{0j}\)). Verification from scratch: with
\(\mathcal F[|y|^{-\alpha}]=c_{3,\alpha}|\xi|^{\alpha-3}\),
\(c_{3,\alpha}=\frac{2^{3-\alpha}\pi^{3/2}\Gamma(\frac{3-\alpha}2)}{\Gamma(\frac\alpha2)}\)
(no pole for \(\alpha\in\{-1,1,2,4,6\}\)):
\(c_{3,1}=4\pi\), \(c_{3,2}=2\pi^2\), \(c_{3,-1}=-8\pi\). Hence
\(\mathcal F[|y|^{-1}\delta_{ij}/(8\pi)]=\frac{\delta_{ij}}{2|\xi|^2}\).
For the second piece use the pointwise identity
\(y_iy_j|y|^{-3}=\delta_{ij}|y|^{-1}-\partial_i\partial_j|y|\) (differentiate
\(\partial_i|y|=y_i/|y|\) once more) together with
\(\mathcal F[|y|]=c_{3,-1}|\xi|^{-4}=-8\pi|\xi|^{-4}\), so
\(\mathcal F[\partial_i\partial_j|y|]=-\xi_i\xi_j\mathcal F[|y|]=8\pi\xi_i\xi_j|\xi|^{-4}\)
and \(\mathcal F[y_iy_j|y|^{-3}]=\frac{4\pi\delta_{ij}}{|\xi|^2}-\frac{8\pi\xi_i\xi_j}{|\xi|^4}\).
Therefore
\(\mathcal F[\mathcal O_{ij}]=\frac1{8\pi}\big(\frac{4\pi\delta_{ij}}{|\xi|^2}
+\frac{4\pi\delta_{ij}}{|\xi|^2}-\frac{8\pi\xi_i\xi_j}{|\xi|^4}\big)
=\frac{\delta_{ij}}{|\xi|^2}-\frac{\xi_i\xi_j}{|\xi|^4}\), as claimed. (All
identities are between homogeneous tempered distributions; ambiguities live at
\(\xi=0\) only and are excluded by homogeneity of degree \(-2>-3\), which is
locally integrable.) Consequently
\(\widetilde h=(2\pi)^3\mathcal O v_0=\frac{(2\pi)^3}{8\pi}W=\pi^2W\).

**(b) The remainder.** \(F-h=-(1-\chi(|\cdot|))h=:-G\). \(G\) vanishes on
\(\{|\xi|<\frac12\}\), is \(C^4\) on \(\mathbb R^3\), and since \(h\) is
homogeneous of degree \(-2\) and smooth off the origin,
\(|\partial^\beta G(\xi)|\le C_\beta(\chi)\|v_0\|\,\langle\xi\rangle^{-2-|\beta|}\)
for every \(|\beta|\le4\).
For \(2\le|\beta|\le4\) this is in \(L^1(\mathbb R^3)\) (the exponent
\(2+|\beta|>3\)). From
\(\widetilde G(y)=\int G e^{i\xi y}d\xi\) (tempered) and the distributional
identity \(y^\beta\widetilde G=i^{|\beta|}\widetilde{\partial^\beta G}\), the
right-hand side is for \(2\le|\beta|\le4\) a bounded continuous function with sup
\(\le\|\partial^\beta G\|_{L^1}\). Hence \(\widetilde G\) agrees on
\(\{|y|\ge1\}\) with a function obeying \(|\widetilde G(y)|\le C_M|y|^{-M}\) for
every \(2\le M\le4\) (four derivatives are all that Definition 0.1 supplies, and
\(M=4\) is all that §3.4 consumes). Set \(\rho=-\widetilde G\); then
\(V=\widetilde F=\widetilde h-\widetilde G=\pi^2W+\rho\). ∎

So the decay of \(V\) is **exactly** \(|y|^{-1}\): the profile
\(\pi^2W\) is \(\chi\)-independent, and \(W\cdot v_0=\frac{\|v_0\|^2}{|y|}+\frac{(v_0\cdot y)^2}{|y|^3}\ge\frac{\|v_0\|^2}{|y|}>0\).

### 3.4 Why the naive periodisation argument is unavailable

**Remark 3.4 (the naive periodisation argument is unavailable — the route is
abandoned).**
The reduction sketch supplied to this task (item R3) asserts an *absolutely
convergent* Poisson-summation identity
\(u_N(x)=N\sum_{m\in\mathbb Z^3}V(N(x+2\pi m))\) with \(m\neq0\) terms
\(O(1/N)\) uniformly. That is **false**: **the naive absolutely convergent
periodisation argument is unavailable**, because the profile decays only like
\(|y|^{-1}\) with a sign-definite spherical mean, so the series of translates
diverges. (An earlier version of this document made a stronger claim about
summation methods than is needed; what is asserted below is only the divergence
of the series and the failure of Abel summation in particular.)

*The series diverges.* By Theorem 3.3 the far field of \(V\) is exactly
\(\pi^2W\) plus a rapidly decaying remainder, and
\(v_0\cdot W(y)=\frac{\|v_0\|^2}{|y|}+\frac{(v_0\cdot y)^2}{|y|^3}
\ge\frac{\|v_0\|^2}{|y|}>0\): \(V\) decays only like \(|y|^{-1}\), with a
*non-oscillating* leading profile whose spherical mean is **sign-definite in the
\(v_0\) direction**, so no cancellation is available. Applying Theorem 3.3 with
\(M=4\): for \(m\neq0\) and \(x\in[-\pi,\pi)^3\), \(|N(x+2\pi m)|\ge N\pi\ge1\),
so
\[
v_0\cdot V\big(N(x+2\pi m)\big)\ \ge\ \frac{\pi^2\|v_0\|^2}{N\,|x+2\pi m|}
\ -\ \frac{C_4\|v_0\|}{N^4|x+2\pi m|^4}.
\]
Now \(\sum_{m\neq0}|x+2\pi m|^{-4}<\infty\) (the shell \(|m|_\infty=n\) has
\(O(n^2)\) points and \(|x+2\pi m|\gtrsim n\)), whereas
\(\sum_{m\neq0}|x+2\pi m|^{-1}=\infty\) (the same shell count gives
\(\sum_nn^2\cdot n^{-1}=\infty\)). Hence
\(\sum_{m\neq0}v_0\cdot V(N(x+2\pi m))=+\infty\) for every \(x\) and every \(N\).

*And Abel regularisation does not repair it either.* The terms of the
\(m\neq0\) series are **eventually positive**, and a summation method that is
regular on series of eventually-positive terms — Abel is one — diverges together
with the series: by monotone convergence,
\(\lim_{\varepsilon\downarrow0}\sum_me^{-\varepsilon|m|}a_m=+\infty\) whenever
\(a_m\ge0\) and \(\sum a_m=\infty\). So an Abel-regularised right-hand side is
\(+\infty\) in the \(v_0\) direction while the left-hand side \(u_N(x)\) is
finite; quantitatively the constant mode blows up like \(\varepsilon^{-2}\),
because the spherical mean of \(W\) is \(\tfrac43v_0/r\neq0\). **A previous
Proposition 3.4 of this document, asserting an Abel-regularised identity
(equation (3.2)), is withdrawn as false and has been deleted.** Its offered
proof was independently unsound at two points: \(\widetilde G\notin L^1\) (it
decays like \(|y|^{-1}\)), so the periodisation's \(k\)-th Fourier coefficient
is not obtainable by Fourier inversion as written; and "Abel summation of a
distributionally convergent series with continuous sum converges pointwise" is
not a theorem.

*If a periodisation formula is wanted purely for orientation*, the correct
**renormalised** one is absolutely convergent:
\[
u_N(x)=N\Big[V(Nx)+\sum_{m\neq0}\big(V(N(x+2\pi m))-V(2\pi Nm)\big)\Big]+\text{const},
\]
since the bracketed differences are \(O(|m|^{-2})\) by Theorem 3.3. It is not
used below either.

**Consequence.** The route (R3) is **not available**: neither the naive
absolutely convergent argument nor its Abel-regularised variant works. We
therefore abandon Poisson summation entirely and use
the lattice Riemann-sum route of §4 instead, which is elementary, quantitative,
and additionally delivers the derivative statement. **Nothing anywhere in this
document depends on a periodisation identity** — this subsection is inert.

---

## 4. The \(C^1_{\rm loc}\) inner rescaling — exact identity plus explicit rate

### 4.1 The exact lattice-Riemann identity

For \(\Phi\in L^1(\mathbb R^3)\) bounded away from a lattice-measure-zero set, put
\[
\mathcal R_N[\Phi]\ :=\ N^{-3}\sum_{k\in\mathbb Z^3}\Phi(k/N)\ -\ \int_{\mathbb R^3}\Phi(\xi)\,d\xi .
\]

**Lemma 4.1 (exact inner rescaling). PROVEN.** For every \(y\in\mathbb R^3\) and
every integer \(N\ge2\),
\[
u_N(y/N)=N\big(V(y)+E_N(y)\big),\qquad
\partial_i u_{N,j}(y/N)=N^2\big(\partial_iV_j(y)+E'_{N,ij}(y)\big),
\]
where, with \(\Phi_y(\xi):=F(\xi)e^{i\xi\cdot y}\) and
\(\Phi'_{y,ij}(\xi):=i\xi_iF_j(\xi)e^{i\xi\cdot y}\) (both extended by \(0\) at
\(\xi=0\)),
\[
E_N(y)=\mathcal R_N[\Phi_y],\qquad E'_{N,ij}(y)=\mathcal R_N[\Phi'_{y,ij}].
\]

*Proof.* By (3.1),
\(u_N(y/N)=\sum_{k}N^{-2}F(k/N)e^{i(k/N)\cdot y}=N\cdot N^{-3}\sum_k\Phi_y(k/N)
=N(\int\Phi_y+\mathcal R_N[\Phi_y])=N(V(y)+E_N(y))\), using
\(V(y)=\int\Phi_y\) (Definition 3.1) and \(\Phi_y(0)=F(0)=0\).
For the gradient,
\(\partial_iu_{N,j}(x)=\sum_kik_i\widehat{u_N}_j(k)e^{ik\cdot x}\); at
\(x=y/N\),
\(=N^{-2}\!\sum_k i k_iF_j(k/N)e^{i(k/N)y}
=N^{-2}\cdot N\sum_ki(k_i/N)F_j(k/N)e^{i(k/N)y}
=N^2\cdot N^{-3}\sum_k\Phi'_{y,ij}(k/N)\), and
\(\int\Phi'_{y,ij}=\partial_iV_j(y)\) by Theorem 3.2(3). ∎

**These are identities, not approximations.** Everything now rests on bounding
a Riemann-sum defect for an explicit compactly supported integrand with a
single integrable algebraic singularity.

### 4.2 The rate

**Theorem 4.2 (\(C^1_{\rm loc}\) convergence with explicit rate). PROVEN.**
There is an absolute constant \(A\) such that for every \(R\ge1\) and every
integer \(N\ge 8R\),
\[
\sup_{|y|\le R}|E_N(y)|\ \le\ A\,\|v_0\|\,\frac{c_\chi(1+\log_2N)+R}{N},
\qquad
\sup_{|y|\le R}|E'_N(y)|\ \le\ A\,\|v_0\|\,\frac{c_\chi+R}{N}.
\]
In particular \(\sup_{|y|\le R}(|E_N|+|E'_N|)=O(N^{-1}\log N)\to0\).

*Proof.* Partition \(\mathbb R^3\) into the half-open cubes
\(Q_k=\frac kN+[0,\frac1N)^3\), \(k\in\mathbb Z^3\), of side \(1/N\) and
diameter \(\sqrt3/N\). For any \(\Phi\),
\[
\mathcal R_N[\Phi]=\sum_{k}\int_{Q_k}\big[\Phi(k/N)-\Phi(\xi)\big]\,d\xi ,
\tag{4.1}
\]
valid whenever both the sum and the integral converge absolutely — here they do,
since \(\Phi\) is supported in \(|\xi|\le1\) and \(|\Phi|\lesssim|\xi|^{-2}\).

*Bounds on \(\Phi_y\).* For \(0<|\xi|\le1\) and \(|y|\le R\):
\[
|\Phi_y(\xi)|\le\|v_0\|\,|\xi|^{-2},\qquad
|\nabla\Phi_y(\xi)|\le|\nabla F(\xi)|+R|F(\xi)|
\le 6c_\chi\|v_0\|\,|\xi|^{-3}+R\|v_0\|\,|\xi|^{-2}.
\tag{4.2}
\]
The bound \(|\nabla F|\le6c_\chi\|v_0\||\xi|^{-3}\) follows from the product rule
applied to \(F=\chi(|\xi|)\,|\xi|^{-2}P_\xi v_0\): \(|\nabla(|\xi|^{-2})|=2|\xi|^{-3}\),
\(|\nabla(P_\xi v_0)|\le2\|v_0\||\xi|^{-1}\) (\(P_\xi\) is degree-\(0\)
homogeneous and smooth on the unit sphere with \(\|\nabla P\|\le2\)), and
\(|\chi'|\le c_\chi\) with \(|\xi|\le1\).

*Core.* Let \(\mathcal K_0:=\{k:\ Q_k\cap\{|\xi|\le\sqrt3/N\}\neq\emptyset\}\);
all such \(Q_k\) lie in \(\{|\xi|\le2\sqrt3/N\}\) and \(|k|\le2\sqrt3<4\). The
core contributes to (4.1) at most
\[
\int_{|\xi|\le2\sqrt3/N}|\Phi_y|+N^{-3}\!\!\sum_{0<|k|\le3}\!|\Phi_y(k/N)|
\le \|v_0\|\Big(4\pi\tfrac{2\sqrt3}{N}+\tfrac1N\!\!\sum_{0<|k|\le3}\!|k|^{-2}\Big)
\le\frac{A_0\|v_0\|}{N},
\]
using \(\sum_{0<|k|\le3}|k|^{-2}\le S_3\le384\) (Lemma 2.2).

*Dyadic annuli.* Let \(J:=\lfloor\log_2(N/(4\sqrt3))\rfloor\) and
\(A_j:=\{2^{-j-1}<|\xi|\le2^{-j}\}\), \(0\le j\le J\). Together with the core,
these cover \(\operatorname{supp}\Phi_y\subset\{|\xi|\le1\}\) (since
\(2^{-J-1}\le\sqrt3/N\cdot\)const, absorbed in \(A_0\) by enlarging the core
constant). On \(A_j\), by (4.2),
\(\sup|\nabla\Phi_y|\le 48c_\chi\|v_0\|8^{j}+4R\|v_0\|4^{j}\), and the
\(\frac{\sqrt3}N\)-neighbourhood of \(A_j\) has volume \(\le A_1 2^{-3j}\)
because \(2^{-j}\ge 4\sqrt3/N\). Cell-wise,
\(\big|\int_{Q_k}[\Phi_y(k/N)-\Phi_y]\big|\le\frac{\sqrt3}{N}\sup_{Q_k}|\nabla\Phi_y|\,|Q_k|\)
for cells on which \(\Phi_y\in C^1\), so the \(A_j\)-contribution is
\[
\le\frac{\sqrt3}N\big(48c_\chi\|v_0\|8^j+4R\|v_0\|4^j\big)A_12^{-3j}
=\frac{A_2\|v_0\|}{N}\big(c_\chi+R\,2^{-j}\big).
\]
Summing over \(0\le j\le J\le\log_2N\):
\(\le\frac{A_2\|v_0\|}{N}\big(c_\chi(1+\log_2N)+2R\big)\).
Adding the core bound gives the first display.

*Gradient.* Repeat with \(\Phi'_y\): \(|\Phi'_y|\le\|v_0\||\xi|^{-1}\) and
\(|\nabla\Phi'_y|\le 7c_\chi\|v_0\||\xi|^{-2}+R\|v_0\||\xi|^{-1}\). The core
contributes \(\le A_0'\|v_0\|N^{-2}\); the \(A_j\)-contribution is
\(\le\frac{A_2'\|v_0\|}N(c_\chi4^{j}+R2^{j})2^{-3j}
=\frac{A_2'\|v_0\|}N(c_\chi2^{-j}+R2^{-2j})\), which sums over \(j\ge0\) to
\(\le\frac{2A_2'\|v_0\|}{N}(c_\chi+R)\) — **no logarithm**, because the
singularity of \(\Phi'_y\) is only \(|\xi|^{-1}\). ∎

*(Remark. The derivative statement genuinely needs its own argument — it is not a
consequence of the \(C^0\) statement — and it is the easier of the two: the
\(\xi\)-factor in \(\Phi'_y\) tames the origin. This is the step the proposed
route flagged as delicate.)*

---

## 5. The duality lower bound

### 5.1 The test field

**Definition 5.1.** Let \(\Psi\in C_c^\infty(\mathbb R^3;\mathbb R^3)\) with
\(\nabla\cdot\Psi=0\) and \(\operatorname{supp}\Psi\subset\{|y|\le R\}\),
\(\Psi\not\equiv0\). For \(N>R/\pi\) define on \(\mathbb T^3\)
\[
\psi_N(x):=N^{3/2}\sum_{m\in\mathbb Z^3}\Psi\big(N(x+2\pi m)\big).
\]

**Lemma 5.2. PROVEN.** For \(N>R/\pi\):
\(\psi_N\) is a well-defined \(C^\infty\) real field on \(\mathbb T^3\); it is
divergence-free; it has zero mean; on the fundamental domain
\(x\in[-\pi,\pi)^3\) exactly one term (\(m=0\)) is nonzero and
\(\operatorname{supp}\psi_N\subset\{|x|\le R/N\}\); and
\[
\|\psi_N\|_{L^2(\mathbb T^3)}=(2\pi)^{-3/2}\|\Psi\|_{L^2(\mathbb R^3)}
\quad\text{— independent of }N .
\]

*Proof.* Local finiteness: \(\Psi(N(x+2\pi m))\neq0\) forces
\(|x+2\pi m|\le R/N<\pi\); for \(x\in[-\pi,\pi)^3\) and \(m\neq0\),
\(|x+2\pi m|\ge2\pi|m|_\infty-\pi\ge\pi\), a contradiction. So the sum is
locally a single smooth term, and \(\psi_N\) is smooth and \(2\pi\mathbb Z^3\)-periodic.
Divergence-free: \(\nabla\cdot\psi_N(x)=N^{5/2}(\nabla\cdot\Psi)(N(x+2\pi m))=0\).
Zero mean: \(\int_{\mathbb T^3}\psi_Nd\mu=(2\pi)^{-3}N^{3/2}\int_{\mathbb R^3}\Psi(Ny)dy
=(2\pi)^{-3}N^{-3/2}\int_{\mathbb R^3}\Psi\), and
\(\int\Psi_i=\int\nabla\cdot(y_i\Psi)=0\) since \(\nabla\cdot\Psi=0\) and
\(y_i\Psi\in C_c^\infty\). (So a curl representation \(\Psi=\nabla\times A\) is
not needed: zero mean is automatic for compactly supported divergence-free
fields.)
Norm: \(\|\psi_N\|^2=(2\pi)^{-3}\int_{|x|\le R/N}N^3|\Psi(Nx)|^2dx
=(2\pi)^{-3}N^3N^{-3}\|\Psi\|_{L^2(\mathbb R^3)}^2\). ∎

### 5.2 The pairing

**Theorem 5.3 (duality identity with controlled error). PROVEN.** With
\(I_\Psi:=\int_{\mathbb R^3}\big(V\cdot\nabla V\big)\cdot\Psi\,dy\) and
\(R\ge1\), \(N\ge\max(8R,32)\),
\[
\big\langle \mathbb P(u_N\cdot\nabla u_N),\psi_N\big\rangle_{L^2(\mathbb T^3)}
=(2\pi)^{-3}N^{3/2}\big(I_\Psi+\mathcal E_N\big),
\]
\[
|\mathcal E_N|\ \le\ \|\Psi\|_{L^1(\mathbb R^3)}\Big[
\varepsilon_N\big(\|\nabla V\|_\infty+\varepsilon'_N\big)+\|V\|_\infty\varepsilon'_N\Big]
\ =\ O\!\Big(\frac{\log N}{N}\Big),
\]
where \(\varepsilon_N:=\sup_{|y|\le R}|E_N|\), \(\varepsilon'_N:=\sup_{|y|\le R}|E'_N|\)
are bounded by Theorem 4.2.

*Proof.* \(\mathbb P\) is an orthogonal projection on \(L^2(\mathbb T^3)\), hence
self-adjoint, and \(\mathbb P\psi_N=\psi_N\) (Lemma 5.2: zero-mean and
divergence-free). Therefore
\(\langle\mathbb P(u_N\cdot\nabla u_N),\psi_N\rangle=\langle u_N\cdot\nabla u_N,\psi_N\rangle\).
Using the normalised measure and \(\operatorname{supp}\psi_N\subset\{|x|\le R/N\}\),
then substituting \(x=y/N\) (\(dx=N^{-3}dy\)):
\[
\langle u_N\cdot\nabla u_N,\psi_N\rangle
=(2\pi)^{-3}N^{3/2}N^{-3}\int_{|y|\le R}(u_N\cdot\nabla u_N)(y/N)\cdot\Psi(y)\,dy .
\]
By Lemma 4.1,
\((u_N\cdot\nabla u_N)_j(y/N)=u_{N,i}(y/N)\,\partial_iu_{N,j}(y/N)
=N^3\big(V_i+E_{N,i}\big)\big(\partial_iV_j+E'_{N,ij}\big)\), so the prefactor
\(N^{-3}N^3\) cancels and
\[
\langle u_N\cdot\nabla u_N,\psi_N\rangle
=(2\pi)^{-3}N^{3/2}\int_{|y|\le R}\big[(V+E_N)\cdot\nabla(V)+\cdots\big]\cdot\Psi\,dy
=(2\pi)^{-3}N^{3/2}(I_\Psi+\mathcal E_N),
\]
with \(\mathcal E_N=\int_{|y|\le R}\big[E_{N,i}\partial_iV_j+V_iE'_{N,ij}
+E_{N,i}E'_{N,ij}\big]\Psi_j\,dy\), bounded as stated (all sup-norms finite by
Theorem 3.2). ∎

**Corollary 5.4 (the \(N^3\) lower bound). PROVEN, conditional on \(I_\Psi\neq0\).**
If some admissible \(\Psi\) has \(I_\Psi\neq0\), then
\[
\|\mathbb P(u_N\cdot\nabla u_N)\|_{L^2(\mathbb T^3)}^2
\ \ge\ \frac{|I_\Psi+\mathcal E_N|^2}{(2\pi)^{3}\,\|\Psi\|_{L^2(\mathbb R^3)}^2}\,N^3
\ =\ \frac{|I_\Psi|^2}{(2\pi)^{3}\|\Psi\|_{L^2(\mathbb R^3)}^2}\,N^3\big(1+o(1)\big),
\]
so \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) for all
\(N\ge N_*\), with
\[
c_0=\frac{|I_\Psi|^2}{2\,(2\pi)^{3}\,\|\Psi\|_{L^2(\mathbb R^3)}^2}
\quad\text{and }N_*\text{ determined by }|\mathcal E_N|\le(1-2^{-1/2})|I_\Psi| .
\]

*Proof.* Cauchy–Schwarz:
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2\ge|\langle\mathbb P(u_N\cdot\nabla u_N),\psi_N\rangle|/\|\psi_N\|_2\).
Insert Theorem 5.3 and Lemma 5.2:
the ratio is \((2\pi)^{-3}N^{3/2}|I_\Psi+\mathcal E_N|\big/\big((2\pi)^{-3/2}\|\Psi\|_{L^2(\mathbb R^3)}\big)\).
Square. \(\mathcal E_N\to0\) by Theorem 5.3. ∎

Everything is now reduced to one statement.

---

## 6. (V-NONDEG): \(\mathbb P_{\mathbb R^3}(V\cdot\nabla V)\not\equiv0\) — PROVEN

### 6.1 Reduction of \(I_\Psi\neq0\) to (V-NONDEG)

**Lemma 6.1. PROVEN.** \(V\cdot\nabla V=\nabla\cdot(V\otimes V)\) has the
representation
\[
(V\cdot\nabla V)(y)=\int_{\mathbb R^3}M(\zeta)e^{i\zeta\cdot y}d\zeta,\qquad
M_i(\zeta)=i\,\zeta_j\,T_{ij}(\zeta),\quad T_{ij}:=F_i*F_j ,
\]
with \(T\) real symmetric, supported in \(\{|\zeta|\le2\}\), and \(M\in L^1\cap L^2\).
Consequently \(\mathbb P(V\cdot\nabla V)\in L^2(\mathbb R^3)\) with symbol
\(P_\zeta M(\zeta)\), and
\[
\mathbb P(V\cdot\nabla V)\equiv0
\iff \zeta\times\big(T(\zeta)\zeta\big)=0\ \text{ for a.e. }\zeta .
\]

*Proof.* \(V=\widetilde F\) with \(F\in L^1\); \(\partial_jV_i=\widetilde{i\xi_jF_i}\)
with \(\xi F\in L^1\). Products of inverse transforms of \(L^1\) functions are
inverse transforms of convolutions:
\(V_j\partial_jV_i=\widetilde{\,F_j*(i\cdot_jF_i)\,}\), i.e. the symbol is
\(i\int F_j(\zeta-\eta)\eta_jF_i(\eta)d\eta\). Since
\((\zeta-\eta)\cdot F(\zeta-\eta)=0\) identically, we may replace
\(\eta_jF_j(\zeta-\eta)\) by \(\zeta_jF_j(\zeta-\eta)\), giving
\(M_i(\zeta)=i\zeta_j(F_i*F_j)(\zeta)=i\zeta_jT_{ij}(\zeta)\).
\(T\) is real (\(F\) real), symmetric (convolution is commutative), and supported
in \(|\zeta|\le2\) (sum of two supports in \(|\xi|\le1\)). Integrability: by
Theorem 3.2(1), \(F\in L^p\) for \(p<3/2\); Young's inequality with
\(p\in[\frac43,\frac32)\) gives \(T=F*F\in L^r\), \(\frac1r=\frac2p-1\le\frac12\),
so \(T\in L^r\) with \(r\ge2\); as \(T\) has compact support, \(T\in L^1\cap L^2\)
and so is \(M=i\zeta T\). Plancherel then puts
\(V\cdot\nabla V\) and \(\mathbb P(V\cdot\nabla V)\) (symbol \(P_\zeta M\), with
\(|P_\zeta M|\le|M|\)) in \(L^2(\mathbb R^3)\), and \(\mathbb P(V\cdot\nabla V)=0\)
in \(L^2\) iff \(P_\zeta M=0\) a.e. Finally, for any vector \(w\),
\(P_\zeta w=0\iff w\parallel\zeta\iff\zeta\times w=0\); here \(w=iT(\zeta)\zeta\). ∎

**Lemma 6.2 (from (V-NONDEG) to a test field). PROVEN.** If
\(\mathbb P(V\cdot\nabla V)\not\equiv0\), then there exists
\(\Psi\in C^\infty_c(\mathbb R^3;\mathbb R^3)\) with \(\nabla\cdot\Psi=0\) and
\(I_\Psi=\int(V\cdot\nabla V)\cdot\Psi\neq0\).

*Proof.* \(\mathbb P(V\cdot\nabla V)\) is a nonzero element of
\(L^2_\sigma(\mathbb R^3)\) (Lemma 6.1: in \(L^2\), symbol \(\perp\zeta\), hence
distributionally divergence-free). The space \(C^\infty_{c,\sigma}(\mathbb R^3)\)
of smooth compactly supported divergence-free fields is dense in
\(L^2_\sigma(\mathbb R^3)\) (classical; this is the standard definition of
\(L^2_\sigma\)). Hence some \(\Psi\in C^\infty_{c,\sigma}\) has
\(\langle\mathbb P(V\cdot\nabla V),\Psi\rangle\neq0\). Since
\(V\cdot\nabla V-\mathbb P(V\cdot\nabla V)=\nabla p\) for some \(p\in L^2_{\rm loc}\)
and \(\int\nabla p\cdot\Psi=-\int p\,\nabla\cdot\Psi=0\) (\(\Psi\) compactly
supported), we get \(I_\Psi=\langle\mathbb P(V\cdot\nabla V),\Psi\rangle\neq0\). ∎

*(Constructive variant: one may take \(\Psi\) with \(\widehat\Psi\) a smooth bump
supported in the explicit open cone of Theorem 6.6, mollified to be compactly
supported in \(y\); this yields an explicit, if unenlightening, \(c_0\).)*

### 6.2 The small-\(\zeta\) asymptotics of \(T\): exact and \(\chi\)-independent

Write, after rotating coordinates so that \(v_0=\|v_0\|e_3\) (all objects are
equivariant: \(F,V\) are linear in \(v_0\) and rotate with it; \(T\) is quadratic).

**Definition 6.3.** \(\tau_{ij}:=h_i*h_j\), where \(h=P_\xi v_0/|\xi|^2\)
(no cutoff).

**Lemma 6.4 (\(\tau\) is well defined and homogeneous). PROVEN.** The integral
\(\tau_{ij}(\zeta)=\int h_i(\eta)h_j(\zeta-\eta)d\eta\) converges absolutely for
every \(\zeta\neq0\), and \(\tau\) is homogeneous of degree \(-1\).

*Proof.* \(|h|\le\|v_0\||\eta|^{-2}\). Near \(\eta=0\) and \(\eta=\zeta\) the
integrand is \(O(|\eta|^{-2})\) resp. \(O(|\zeta-\eta|^{-2})\), locally
integrable in \(\mathbb R^3\); at infinity it is \(O(|\eta|^{-4})\), integrable.
Homogeneity: \(h(\lambda\eta)=\lambda^{-2}h(\eta)\), and substituting
\(\eta=\lambda\eta'\) gives \(\tau(\lambda\zeta)=\lambda^{-1}\tau(\zeta)\). ∎

**Theorem 6.5 (closed form of the obstruction at leading order). PROVEN.**
With \(v_0=\|v_0\|e_3\), for every \(\zeta\neq0\),
\[
\boxed{\ \zeta\times\big(\tau(\zeta)\zeta\big)
=\frac{3\pi^3}{8}\,\|v_0\|^2\,\frac{\zeta_3}{|\zeta|}\,(\zeta\times e_3).\ }
\]
In particular \(\zeta\times(\tau(\zeta)\zeta)\neq0\) whenever \(\zeta_3\neq0\)
and \(\zeta\not\parallel e_3\).

*Proof.* Normalise \(\|v_0\|=1\) (both sides are quadratic in \(\|v_0\|\)).
By the convolution theorem, \(\widetilde{\tau_{ij}}=\widetilde{h_i}\,\widetilde{h_j}\),
i.e. \(\tau_{ij}=(2\pi)^{-3}\mathcal F\big[\widetilde h_i\widetilde h_j\big]\).
By Theorem 3.3(a), \(\widetilde h=\pi^2W\) with
\(W(y)=\frac1r\big(e_3+\hat y_3\hat y\big)\), \(r=|y|\), \(\hat y=y/r\). Hence
\[
\widetilde h_i\widetilde h_j
=\pi^4\Big[\frac{\delta_{i3}\delta_{j3}}{r^2}
+\frac{\delta_{i3}y_3y_j+\delta_{j3}y_3y_i}{r^4}
+\frac{y_3^2y_iy_j}{r^6}\Big].
\]
Each bracket is homogeneous of degree \(-2\), locally integrable, so its Fourier
transform is an honest homogeneous distribution of degree \(-1\), given away
from \(\zeta=0\) by the analytic-continuation formula
\(\mathcal F[|y|^{-\alpha}]=c_{3,\alpha}|\zeta|^{\alpha-3}\),
\(c_{3,\alpha}=\frac{2^{3-\alpha}\pi^{3/2}\Gamma(\frac{3-\alpha}2)}{\Gamma(\frac\alpha2)}\).
With \(c_{3,2}=2\pi^2\), \(c_{3,4}=-\pi^2\), \(c_{3,6}=\frac{\pi^2}{12}\)
(no \(\Gamma\)-poles at \(\alpha=2,4,6\)) and
\(\mathcal F[y_ay_bf]=-\partial_a\partial_b\mathcal F[f]\),
\(\mathcal F[y_ay_by_cy_df]=\partial_a\partial_b\partial_c\partial_d\mathcal F[f]\):
\[
\mathcal F\Big[\frac{1}{r^2}\Big]=\frac{2\pi^2}{s},\qquad
\mathcal F\Big[\frac{y_iy_j}{r^4}\Big]=\pi^2\Big(\frac{\delta_{ij}}{s}-\frac{\zeta_i\zeta_j}{s^3}\Big),\qquad
\mathcal F\Big[\frac{y_3^2y_iy_j}{r^6}\Big]=\frac{\pi^2}{12}\,\partial_3^2\partial_i\partial_j s^3,
\]
\(s:=|\zeta|\). (Trace check on the middle one:
\(\pi^2(3/s-1/s)=2\pi^2/s=\mathcal F[r^{-2}]\) ✓.) Direct differentiation gives
\[
\partial_3^2\partial_i\partial_j s^3
=3\Big[\frac{\delta_{ij}+2\delta_{i3}\delta_{j3}}{s}
-\frac{\zeta_3^2\delta_{ij}+2\zeta_3(\delta_{i3}\zeta_j+\delta_{j3}\zeta_i)+\zeta_i\zeta_j}{s^3}
+\frac{3\zeta_i\zeta_j\zeta_3^2}{s^5}\Big].
\]
Assembling, \(\tau_{ij}=\frac{\pi^3}{8}G_{ij}\) with
\[
G_{ij}=\frac{\delta_{ij}}{4s}+\frac{9\delta_{i3}\delta_{j3}}{2s}
-\frac{\zeta_3^2\delta_{ij}}{4s^3}
-\frac{3\zeta_3(\delta_{i3}\zeta_j+\delta_{j3}\zeta_i)}{2s^3}
-\frac{\zeta_i\zeta_j}{4s^3}+\frac{3\zeta_i\zeta_j\zeta_3^2}{4s^5}.
\]
Contracting with \(\zeta_j\) and using \(\zeta_j\zeta_j=s^2\), the
\(\frac{\zeta_i}{4s}\) terms cancel and the \(\zeta_i\zeta_3^2/s^3\) coefficients
sum to \(-\frac14-\frac32+\frac34=-1\):
\[
(G\zeta)_i=\frac{3\zeta_3}{s}\delta_{i3}-\frac{\zeta_3^2}{s^3}\zeta_i,
\qquad\text{i.e.}\qquad G\zeta=\frac{3\zeta_3}{s}e_3-\frac{\zeta_3^2}{s^3}\zeta .
\]
Hence \(\zeta\times(G\zeta)=\frac{3\zeta_3}{s}(\zeta\times e_3)\), and
\(\zeta\times(\tau\zeta)=\frac{3\pi^3}{8}\frac{\zeta_3}{|\zeta|}(\zeta\times e_3)\). ∎

#### Independent verification of Theorem 6.5 in physical space

The Fourier computation uses finite-part regularisations; here is a
**completely elementary, regularisation-free confirmation**.

**Lemma 6.5′. PROVEN (elementary calculus).** For
\(W(y)=\frac{e_3}{r}+\frac{y_3y}{r^3}\) on \(\mathbb R^3\setminus\{0\}\):
\[
\nabla\cdot W=0,\qquad \omega:=\nabla\times W=\frac{2\,(e_3\times y)}{r^3},
\qquad
\nabla\times\big(W\times\omega\big)=\frac{12\,y_3\,(e_3\times y)}{r^6}\ \not\equiv0 .
\]
Since \(W\cdot\nabla W=\nabla\frac{|W|^2}2-W\times\omega\), it follows that
\(\nabla\times(W\cdot\nabla W)\neq0\), i.e. \(W\cdot\nabla W\) is **not** a
gradient, i.e. \(\mathbb P(W\cdot\nabla W)\neq0\): **the pure Oseen field is not
a stationary Euler flow.**

*Proof.* \(\nabla\times(r^{-1}e_3)=\nabla(r^{-1})\times e_3=-r^{-3}y\times e_3
=r^{-3}(e_3\times y)\).
\(\big[\nabla\times(y_3y/r^3)\big]_a=\epsilon_{abi}\partial_b(y_3y_ir^{-3})
=\epsilon_{abi}\big(\delta_{b3}y_i+y_3\delta_{bi}\big)r^{-3}-3\epsilon_{abi}y_3y_iy_br^{-5}
=\epsilon_{a3i}y_ir^{-3}=(e_3\times y)_a r^{-3}\) (the last two terms vanish by
antisymmetry). Adding, \(\omega=2(e_3\times y)/r^3\).
Then
\[
W\times\omega=\Big(\frac{e_3}{r}+\frac{y_3y}{r^3}\Big)\times\frac{2(e_3\times y)}{r^3}
=\frac{2(y_3e_3-y)}{r^4}+\frac{2y_3e_3}{r^4}-\frac{2y_3^2y}{r^6}
=\frac{4y_3e_3}{r^4}-\frac{2y}{r^4}-\frac{2y_3^2y}{r^6},
\]
using \(e_3\times(e_3\times y)=y_3e_3-y\) and \(y\times(e_3\times y)=r^2e_3-y_3y\).
Now curl each piece. \(\nabla\times(f(r)\,y)=\nabla f\times y=f'\hat y\times y=0\),
so \(-2y/r^4\) contributes nothing. For \(g\,y\) with \(g=-2y_3^2/r^6\):
\(\nabla g=-\frac{4y_3e_3}{r^6}+\frac{12y_3^2y}{r^8}\), so
\(\nabla\times(gy)=\nabla g\times y=-\frac{4y_3(e_3\times y)}{r^6}\).
For \(f\,e_3\) with \(f=4y_3/r^4\):
\(\nabla f=\frac{4e_3}{r^4}-\frac{16y_3y}{r^6}\), so
\(\nabla\times(fe_3)=\nabla f\times e_3=-\frac{16y_3(y\times e_3)}{r^6}
=\frac{16y_3(e_3\times y)}{r^6}\). Total: \(\frac{12y_3(e_3\times y)}{r^6}\). ∎

**Consistency of the two computations. PROVEN.** The curl of
\(\mathbb P(\widetilde h\cdot\nabla\widetilde h)\) has symbol
\(i\zeta\times\big(iP_\zeta\tau\zeta\big)=-\zeta\times(\tau\zeta)\) (using
\(\zeta\times P_\zeta A=\zeta\times A\)). With \(\|v_0\|=1\),
Theorem 6.5 gives \(-\frac{3\pi^3}{8}\frac{\zeta_3(\zeta\times e_3)}{|\zeta|}\).
Its inverse transform: writing the \(a\)-component as
\(\epsilon_{ab3}\zeta_b\zeta_3/|\zeta|\) and using
\(\widetilde{|\zeta|^{-1}}=c_{3,1}r^{-2}=4\pi r^{-2}\),
\(\widetilde{\zeta_b\zeta_3 G}=-\partial_b\partial_3\widetilde G\):
\[
\widetilde{\Big[\frac{\zeta_b\zeta_3}{|\zeta|}\Big]}
=-4\pi\,\partial_b\partial_3 r^{-2}
=\frac{8\pi\delta_{b3}}{r^4}-\frac{32\pi y_3y_b}{r^6},
\]
and \(\epsilon_{ab3}\delta_{b3}=0\), leaving
\(-\frac{32\pi y_3(y\times e_3)_a}{r^6}=+\frac{32\pi y_3(e_3\times y)_a}{r^6}\).
Multiplying by \(-\frac{3\pi^3}{8}\) gives \(-12\pi^4\frac{y_3(e_3\times y)}{r^6}\),
which is exactly \(-\pi^4\) times the answer of Lemma 6.5′ (recall
\(\widetilde h=\pi^2W\), so \(\widetilde h\cdot\nabla\widetilde h=\pi^4W\cdot\nabla W\)
and \(\nabla\times(\widetilde h\cdot\nabla\widetilde h)=-\pi^4\nabla\times(W\times\omega)\)).
**The two independent derivations agree exactly, constants included.**

**NUMERICAL cross-check (float, central differences \(h=10^{-4}\),
`check_lstar.py` §A).** At five random points, \(\nabla\times(W\times\omega)\)
matched \(12y_3(e_3\times y)/r^6\) with relative error \(10^{-8}\)–\(10^{-6}\)
(finite-difference truncation level); \(\nabla\cdot W=-3\times10^{-12}\) and
the \(\omega\) formula matched to \(3\times10^{-11}\).

### 6.3 Transfer to the cutoff family

**Lemma 6.6 (cutoff perturbation is bounded). PROVEN.** For all \(|\zeta|\le\frac18\)
and all \(i,j\),
\[
|T_{ij}(\zeta)-\tau_{ij}(\zeta)|\ \le\ C_1:=64\pi\|v_0\|^2 ,
\]
uniformly over admissible \(\chi\).

*Proof.* \(F=\chi(|\cdot|)h\), so
\[
T_{ij}(\zeta)-\tau_{ij}(\zeta)
=\int h_i(\eta)h_j(\zeta-\eta)\big[\chi(|\eta|)\chi(|\zeta-\eta|)-1\big]d\eta .
\]
The bracket vanishes when \(|\eta|\le\frac12\) **and** \(|\zeta-\eta|\le\frac12\).
For \(|\zeta|\le\frac18\), \(|\eta|\le\frac14\) implies \(|\zeta-\eta|\le\frac38<\frac12\),
so the bracket vanishes there; hence the integrand is supported in
\(\{|\eta|\ge\frac14\}\), where also \(|\zeta-\eta|\ge|\eta|-\frac18\ge\frac{|\eta|}2\).
With \(|h|\le\|v_0\||\cdot|^{-2}\) and \(|\chi\chi-1|\le1\),
\[
|T_{ij}-\tau_{ij}|\le\|v_0\|^2\!\!\int_{|\eta|\ge1/4}\!\!\frac{d\eta}{|\eta|^2|\zeta-\eta|^2}
\le4\|v_0\|^2\!\!\int_{|\eta|\ge1/4}\!\!\frac{d\eta}{|\eta|^4}
=4\|v_0\|^2\cdot4\pi\!\int_{1/4}^\infty\!\!\frac{dr}{r^2}=64\pi\|v_0\|^2 . \qquad\Box
\]

**Theorem 6.7 ((V-NONDEG)). PROVEN.** For every \(v_0\neq0\) and every
admissible \(\chi\),
\[
\mathbb P_{\mathbb R^3}\big(V\cdot\nabla V\big)\ \not\equiv\ 0 .
\]
Explicitly, with \(v_0=\|v_0\|e_3\), on the open set
\[
\mathcal U:=\Big\{\zeta:\ 0<|\zeta|<\tfrac{\pi^2}{2048},\
|\zeta_3|\ge\tfrac{|\zeta|}2,\ |\zeta\times e_3|\ge\tfrac{|\zeta|}2\Big\}
\quad(\text{positive measure}),
\]
one has \(\zeta\times(T(\zeta)\zeta)\neq0\), hence \(P_\zeta M(\zeta)\neq0\).

*Proof.* Let \(\zeta\in\mathcal U\). By Theorem 6.5,
\[
\big|\zeta\times(\tau(\zeta)\zeta)\big|
=\frac{3\pi^3}{8}\|v_0\|^2\frac{|\zeta_3|\,|\zeta\times e_3|}{|\zeta|}
\ \ge\ \frac{3\pi^3}{8}\|v_0\|^2\cdot\frac{|\zeta|}{4}
=\frac{3\pi^3}{32}\|v_0\|^2|\zeta| .
\]
By Lemma 6.6 (applicable since \(\frac{\pi^2}{2048}<\frac18\)),
\(\|T(\zeta)-\tau(\zeta)\|_{\rm F}\le3C_1\), so
\(\big|\zeta\times\big((T-\tau)\zeta\big)\big|\le3C_1|\zeta|^2=192\pi\|v_0\|^2|\zeta|^2\).
Therefore
\[
\big|\zeta\times(T(\zeta)\zeta)\big|
\ \ge\ \|v_0\|^2|\zeta|\Big(\frac{3\pi^3}{32}-192\pi|\zeta|\Big)\ >\ 0
\quad\text{whenever }|\zeta|<\frac{3\pi^3}{32\cdot192\pi}=\frac{\pi^2}{2048}.
\]
By Lemma 6.1 this means \(P_\zeta M\neq0\) on a set of positive measure, so
\(\mathbb P(V\cdot\nabla V)\neq0\) in \(L^2(\mathbb R^3)\). ∎

**Three remarks.**
1. **The proof is uniform in \(\chi\) and in \(v_0\):** \(\|v_0\|^2\) cancels,
   and the only properties of \(\chi\) used are \(0\le\chi\le1\) and
   \(\chi\equiv1\) on \([0,\frac12]\). No smoothness of \(\chi\) enters §6 at all
   (\(\chi'\) is used only in §4, through \(c_\chi\); the higher derivatives
   permitted by Definition 0.1 are used only in the inert Theorem 3.3(b)/§3.4).
   **This uniformity is qualitative only.** It says that the conclusion
   \(\mathbb P(V\cdot\nabla V)\not\equiv0\) holds for every \((\chi,v_0)\); it
   does **not** say that the constants \(c_0,N_*\) of Theorem 7.1(3) are
   \((\chi,v_0)\)-uniform. They are not: they are extracted only after a test
   field \(\Psi\) has been chosen for the particular profile
   \(V=V^{\chi,v_0}\) (Lemma 6.2), and they depend on that pair — and are
   non-effective (§5/§7.2).
2. **This is route (b) of the task, completed, not route (a).** No specific
   \(\chi\) and no numerically evaluated convolution integral is needed: the
   obstruction is carried entirely by the \(\chi\)-independent leading
   singularity of \(T\) at \(\zeta=0\), which is the \emph{same} for every
   admissible cutoff and is the Fourier statement of "the Oseen field is not a
   stationary Euler flow" (Lemma 6.5′). Route (a) (exact coefficient for a
   piecewise-polynomial \(\chi\)) and route (c) (interval arithmetic) are
   therefore not required.
3. Theorem 6.5 shows the obstruction is \emph{maximal} away from the two
   degenerate directions \(\zeta\parallel v_0\) and \(\zeta\perp v_0\); by the
   sign-flip symmetries of §1.2 one checks directly that
   \(T(\zeta)\zeta\parallel\zeta\) exactly on those two sets, so \(\mathcal U\)
   is essentially optimal.

---

## 7. The theorem, and what it does to the Proposition

### 7.1 Main result

**Theorem 7.1 ((L\*) for the smoothly truncated family). PROVEN.**
Let \(v_0\in\mathbb R^3\setminus\{0\}\), let \(\chi\) be admissible
(Definition 0.1), and let \(u_N=u_N^{\chi,v_0}\) be as in Definition 1.1. Then:

1. \(u_N\) is a real, zero-mean, exactly divergence-free trigonometric
   polynomial with band \(|k|\le N\), and the exact laws of Theorem 1.3 hold.
2. \(\frac{N}{44544}\le N_0^2(u_N)\le\frac{64}{3}N\) for \(N\ge32\).
3. There exist \(c_0=c_0(\chi,v_0)>0\) and \(N_*=N_*(\chi,v_0)\) with
\[
\boxed{\ \big\|\mathbb P(u_N\cdot\nabla u_N)\big\|_{L^2(\mathbb T^3)}^2\ \ge\ c_0\,N^3
\qquad\text{for all }N\ge N_* .\ }
\]
   Formally \(c_0=\dfrac{|I_\Psi|^2}{2(2\pi)^3\|\Psi\|_{L^2(\mathbb R^3)}^2}\)
   for any admissible \(\Psi\) with \(I_\Psi\neq0\), which exists by
   Theorem 6.7 + Lemma 6.2.

**Quantifiers (important).** The bound in (3) holds **for every** admissible
cutoff profile \(\chi\) and **every** nonzero seed vector \(v_0\). It is
**not** uniform in \((\chi,v_0)\): \(c_0\) and \(N_*\) depend on that pair,
through \(V=V^{\chi,v_0}\) and through the test field \(\Psi\) chosen for it.
What *is* \((\chi,v_0)\)-uniform is the qualitative non-degeneracy
\(\mathbb P(V\cdot\nabla V)\not\equiv0\) of Theorem 6.7, which is exactly what
licenses the "for every" quantifier.

**Effectivity caveat (important, and load-bearing for the paper's
"constant dependence" claim).** Part (3) is an **existence** statement:
\(c_0\) and \(N_*\) are **not effective**. Lemma 6.2 produces \(\Psi\) from the
*density* of \(C^\infty_{c,\sigma}\) in \(L^2_\sigma\), a pure existence
argument, so none of \(\Psi\), its support radius \(R\), \(I_\Psi\), \(c_0\) or
\(N_*\) is exhibited. The "constructive variant" after Lemma 6.2 is available in
principle but expensive: Theorem 6.7 certifies non-vanishing only on
\(|\zeta|<\pi^2/2048\approx4.8\times10^{-3}\), so an explicitly constructed
\(\Psi\) has Fourier support at that scale, hence spatial radius
\(R\gtrsim10^3\), hence \(N_*\ge8R\gtrsim10^4\) and a \(c_0\) many orders of
magnitude below the measured \(\approx8.4\|v_0\|^4\) of §7.5. Everything else in
this document — Theorems 1.3, 3.2, 3.3, 4.2, 5.3, 6.5, 6.7, Proposition 2.4,
Lemmas 2.1–2.3, 4.1, 5.2, 6.1, 6.6 — carries fully explicit constants; only
\(c_0\) and \(N_*\) do not.

**The exponent achieved is \(3\) — the sharp one, \(a=1\) in (R1).** No weakening
was needed.

*Proof.* (1) Lemma 1.2 + Theorem 1.3. (2) Proposition 2.4.
(3) Theorem 6.7 gives \(\mathbb P(V\cdot\nabla V)\not\equiv0\); Lemma 6.2 gives
\(\Psi\in C^\infty_{c,\sigma}\) with \(I_\Psi\neq0\); Corollary 5.4 (which rests
on Lemma 4.1, Theorem 4.2, Lemma 5.2, Theorem 5.3) converts this into the
displayed bound. ∎

### 7.2 The Proposition becomes unconditional

**Theorem 7.2 (static no-go, unconditional). PROVEN.**
Let \(\Phi:[0,\infty)\to(0,\infty)\) be nondecreasing and suppose
\(K(u)\le\Phi(\log N_0^2(u))\) for every real, zero-mean, divergence-free
trigonometric field \(u\) on \(\mathbb T^3\). Then \(\Phi(s)\ge c\,e^{s}\) for all
large \(s\), and \(\int^\infty ds/\Phi<\infty\): **no Osgood-admissible \(\Phi\)
satisfies a uniform pointwise field inequality**, so clause (e) of the Main
Theorem cannot be activated through the \(R\equiv0\) route.

*Proof.* Fix any admissible \(\chi,v_0\) and use \(u_N\). By Theorem 1.3 and
Proposition 2.4, \(H_1=\frac23\|v_0\|^2S_N^\chi\le\frac{256}{3}\|v_0\|^2N\), so
by Theorem 7.1(3), for \(N\ge N_*\),
\[
K(u_N)=\frac{\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2}{H_1^2}
\ \ge\ \frac{c_0N^3}{\big(\tfrac{256}3\|v_0\|^2N\big)^2}
=\frac{9c_0}{65536\|v_0\|^4}\,N .
\]
By Proposition 2.4, \(N\ge\frac{3}{64}N_0^2(u_N)\), hence
\(K(u_N)\ge c_0''\,N_0^2(u_N)=c_0''e^{s_N}\) with
\(c_0''=\frac{27c_0}{4194304\|v_0\|^4}\) and \(s_N=\log N_0^2(u_N)\).

*A sequence with divergent \(s\) and bounded gaps.* Proposition 2.4 gives
\(N_0^2(u_N)\in[\frac{N}{44544},\frac{64}{3}N]\), hence for any \(q\)
\[
\frac{N_0^2(u_{2^qN})}{N_0^2(u_N)}\in\Big[\frac{3\cdot2^q}{44544\cdot64},\
2^q\cdot\frac{44544\cdot64}{3}\Big]
=\Big[\frac{2^q}{950272},\ 2^q\cdot950272\Big].
\]
Choose \(q=20\) (\(2^{20}=1048576>950272\)) and \(N_j:=2^{20j}N_\sharp\) with
\(N_\sharp\ge\max(N_*,32)\). Then
\[
0<c_-:=\log\frac{1048576}{950272}=0.0985\ldots\ \le\ s_{N_{j+1}}-s_{N_j}\ \le\
\log\big(2^{20}\cdot950272\big)=27.63\ldots=:c_+ ,
\]
so \(s_{N_j}\nearrow\infty\) with gaps in \([c_-,c_+]\).

*Conclusion (verbatim the paper's argument).* \(\Phi(s_{N_j})\ge K(u_{N_j})\ge
c_0''e^{s_{N_j}}\); for \(s\in[s_{N_j},s_{N_{j+1}}]\), monotonicity gives
\(\Phi(s)\ge\Phi(s_{N_j})\ge c_0''e^{s_{N_j}}\ge c_0''e^{-c_+}e^{s}\). Hence
\(\Phi(s)\ge ce^s\) for all \(s\ge s_{N_0}\) with \(c=c_0''e^{-c_+}\), and
\(\int^\infty ds/\Phi\le c^{-1}\int^\infty e^{-s}ds<\infty\). ∎

### 7.3 Status ledger

| Statement | Status |
|---|---|
| Reduction (R1) (weaker exponent \(2+a\) suffices) | **PROVEN** (§0.3) — and not needed |
| (R2) exact moment laws for smooth truncation | **PROVEN** (Thm 1.3), verified NUMERICALLY (float) |
| Lemma 9′/9″ (general radial weight, lattice and continuum) | **PROVEN** (§1.2) |
| \(N_0^2\asymp N\), explicit crude constants | **PROVEN** (Prop 2.4) |
| \(V\) continuous, bounded, \(C^\infty\), entire of type \(\le1\), div-free, \(V(0)\neq0\) | **PROVEN** (Thm 3.2) |
| \(V=\pi^2W+O(|y|^{-4})\), \(W\) = Oseen profile (four derivatives of \(\chi\); more if \(\chi\) is smoother; used only by the inert §3.4) | **PROVEN** (Thm 3.3) |
| (R3) Poisson identity with absolutely convergent \(m\neq0\) tail, \(O(1/N)\) | **FALSE** — the naive absolutely convergent periodisation argument is unavailable: the series diverges (the profile decays only like \(|y|^{-1}\) with a sign-definite spherical mean), and since the terms are eventually positive Abel summation does not converge either (Remark 3.4). Used nowhere |
| Exact lattice-Riemann identity \(u_N(y/N)=N(V+E_N)\), \(\nabla u_N(y/N)=N^2(\nabla V+E'_N)\) | **PROVEN** (Lem 4.1) |
| \(\sup_{|y|\le R}(|E_N|+|E'_N|)=O(N^{-1}\log N)\), constants explicit | **PROVEN** (Thm 4.2) |
| Duality lower bound \(\ge c_0N^3(1+o(1))\) given \(I_\Psi\neq0\) | **PROVEN** (Cor 5.4) |
| \(\zeta\times(\tau\zeta)=\frac{3\pi^3}{8}\|v_0\|^2\frac{\zeta_3(\zeta\times e_3)}{|\zeta|}\) | **PROVEN** (Thm 6.5), cross-checked exactly in physical space (Lem 6.5′) and NUMERICALLY (float) |
| (V-NONDEG) \(\mathbb P(V\cdot\nabla V)\neq0\), all admissible \(\chi\), all \(v_0\) | **PROVEN** (Thm 6.7) |
| (L\*) for the smooth family, exponent \(3\) | **PROVEN** (Thm 7.1) |
| Proposition (static no-go) | **UNCONDITIONAL** (Thm 7.2) |
| (L\*) for the **sharply** truncated family \(\chi=\mathbb 1_{[0,1]}\) | **OPEN** — and no longer needed |

### 7.4 Why the sharp cutoff is not covered

The one and only place the smoothness of \(\chi\) is used *load-bearingly* is
Theorem 4.2, and one derivative is all it uses:
\(|\nabla F|\lesssim c_\chi|\xi|^{-3}\) fails for \(\chi=\mathbb 1_{[0,1]}\),
where \(F\)
jumps across \(|\xi|=1\). The Riemann-sum defect then carries a lattice-point
(Gauss circle) term of size \(O(N^{-1+\theta})\) *per unit* rather than
\(O(N^{-1}\log N)\), and — more seriously — \(V\) itself loses the \(|y|^{-4}\)
tail correction and acquires an oscillatory \(|y|^{-2}\) tail. None of §§5–6
would change (they are insensitive to \(\chi\)), so the sharp-cutoff case is a
purely technical gap in the inner-rescaling step, not a structural one. We flag
it as **OPEN** rather than closing it, because the Proposition does not need it.

### 7.5 Numerical corroboration (float; not used in any proof)

`check_lstar.py` §C, pseudo-spectral on an \(M^3\) grid with \(M=4N+2\)
(no aliasing: the product of two fields banded at \(N\) is banded at \(2N\)):

| \(N\) | \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\) | \(/N^3\) | \(K\) | \(K/N\) |
|---|---|---|---|---|
| 4 | 163.13 | 2.549 | 0.522 | 0.1305 |
| 8 | 2666.07 | 5.207 | 1.574 | 0.1968 |
| 12 | 10972.67 | 6.350 | 2.620 | 0.2183 |
| 16 | 28545.74 | 6.969 | 3.664 | 0.2290 |
| 24 | 105362.68 | 7.622 | 5.749 | 0.2396 |
| 32 | 260856.87 | 7.961 | 7.834 | 0.2448 |
| 40 | 522767.81 | 8.168 | 9.918 | 0.2480 |

\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2/N^3\) rises monotonically and appears to
saturate near \(8.4\pm0.2\); \(K/N\) saturates near \(0.25\). This is consistent
with the (unproven here, plausible) sharp limit
\(\lim_N N^{-3}\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2
=(2\pi)^{-3}\|\mathbb P(V\cdot\nabla V)\|_{L^2(\mathbb R^3)}^2\),
which would require controlling the pairing against \emph{all} test fields on
the torus rather than one localised \(\Psi\). **Nothing above depends on this.**

---

## 8. Summary of what changed relative to the proposed route

* **(R1)** verified; not needed (we achieve the sharp exponent).
* **(R2)** verified in full; Lemma 9 generalised to arbitrary radial weights
  (Lemma 9′) and to \(\mathbb R^3\) (Lemma 9″).
* **(R3)** the proposed absolutely convergent Poisson identity is **false**, and
  so is every regularised version of it (Remark 3.4);
  replaced by an exact lattice-Riemann identity (Lemma 4.1) with an explicit
  \(O(N^{-1}\log N)\) rate, which additionally delivers the derivative statement.
* **(R4)** carried out as proposed, with the simplification that a compactly
  supported divergence-free \(\Psi\) automatically has zero mean (no curl
  representation needed).
* **(R5)** the "only genuinely nontrivial remaining input" is **closed**, by a
  \(\chi\)-independent argument: the leading small-\(\zeta\) singularity of
  \(F*F\) is the Oseen convolution \(\tau\), for which
  \(\zeta\times(\tau\zeta)=\frac{3\pi^3}{8}\|v_0\|^2\frac{\zeta_3(\zeta\times e_3)}{|\zeta|}\neq0\)
  — verified twice, once in Fourier and once by elementary vector calculus in
  physical space (Lemma 6.5′: \(\nabla\times(W\times\omega)=12y_3(e_3\times y)/r^6\)).
