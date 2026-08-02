# Theorem statement

**A bandwidth–dissipation dichotomy for the three-dimensional periodic
Navier–Stokes equations, with an exact action representation**

This file states the single main theorem and the one clearly-separated
**unconditional** obstruction theorem (Theorem O; formerly a conditional
proposition). Proofs: [complete_proof.md](complete_proof.md); the capacity
input to Theorem O is proven in
[lstar/lstar_proof_main.md](lstar/lstar_proof_main.md).
Dependency and gap audit (external classical inputs, adversarial
re-verification, novelty comparison):
[dependency_and_gap_audit.md](dependency_and_gap_audit.md).
LaTeX version: [paper_draft.tex](paper_draft.tex).
**No Clay statement is claimed in either direction.**

## Setting

\(\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3\), normalised measure;
\(u(x,t)=\sum_{k\in\mathbb Z^3\setminus\{0\}}\hat u_k(t)e^{ik\cdot x}\),
zero mean, \(k\cdot\hat u_k=0\), \(\hat u_{-k}=\overline{\hat u_k}\);
viscosity \(\nu>0\); no forcing:
\(\partial_tu+\mathbb P(u\cdot\nabla u)=\nu\Delta u\), \(\mathbb P\) the
Leray projection. For \(u_0\in H^m_\sigma\), \(m>5/2\), \(u_0\neq0\), let
\(u\in C([0,T_{\max});H^m)\cap C^1([0,T_{\max});H^{m-2})\) be the maximal
strong solution (framework input F1 of the audit). Define, with
\(H_r(t)=\sum_k|k|^{2r}|\hat u_k(t)|^2\):

\[
D=\|\nabla u\|_2^2=H_1,\qquad
N_0^2=\frac{H_1}{H_0},\qquad z=\log N_0^2\ \ (\ge0),
\]
\[
\mathcal N=-\mathbb P(u\cdot\nabla u),\qquad
K=\frac{\|\mathcal N\|_2^2}{\|\nabla u\|_2^4},\qquad
N_1^2=\frac{H_2}{H_1}.
\]

\(K\) is a scale- and amplitude-covariant nonlinear front wavenumber
functional; \(N_0^2\) is the energy-weighted mean square wavenumber
(spectral bandwidth); all quantities are finite and, where differentiated
below, continuously differentiable on \([0,T_{\max})\) (Lemma 0).

## Main Theorem

**Theorem.** Let \(u\) be as above. Then:

**(a) (Bandwidth monotone.)**
\[
\Lambda(t)\;=\;\log N_0^2(t)-\frac1{2\nu}\int_0^tK(s)\,D(s)\,ds
\]
is non-increasing on \([0,T_{\max})\). Quantitatively,
\[
\frac{d}{dt}\log N_0^2
=\frac{2\,\mathrm{Cov}_{p}(x,g)}{N_0^2}-\frac{2\nu\,V}{N_0^2}
\;\le\;\frac{K\,D}{2\nu},
\]
where \(p(k)=|\hat u_k|^2/H_0\), \(x_k=|k|^2\),
\(\mathrm{Cov}=\sum_k(x_k-N_0^2)\,
\mathrm{Re}\langle\hat u_k,\hat{\mathcal N}_k\rangle/H_0\),
\(V=\mathrm{Var}_p(x)\), and the inequality is a modal Cauchy–Schwarz
bound followed by a square completion, with an explicit nonnegative
two-part defect.

**(b) (Dichotomy.)** Exactly one of the following holds:
1. \(\displaystyle\int_0^{T_{\max}}K\,D\,dt<\infty\); then
   \(z\) is bounded, \(u\in L^\infty(0,T_{\max};H^1)\), and
   \(T_{\max}=\infty\) (the solution is global);
2. \(\displaystyle\int_0^{T_{\max}}K\,D\,dt=\infty\).

In particular \(T_{\max}<\infty\Rightarrow\int_0^{T_{\max}}KD\,dt
=\infty\). (No claim is made that alternative 2 forces
\(T_{\max}<\infty\).)

**(c) (Exact action representation.)** For every \(T'<T_{\max}\),
\[
\int_0^{T'}K\,D\,dt
=\int_0^{T'}\frac{\|\partial_tu\|_2^2}{D}\,dt
\;+\;\nu^2\int_0^{T'}N_1^2\,dt
\;+\;\nu\,\log\frac{D(T')}{D(0)} .
\]
Moreover, by AM–GM on \(\nu\dot D/D\),
\(\int_0^{T'}KD\,dt\le2\int_0^{T'}(\|\partial_tu\|_2^2/D+\nu^2N_1^2)dt\):
**finiteness of the \(\dot H^1\)-bandwidth action implies alternative
(b1)**. The converse is false (explicit decaying counterexample
\(u=e^{-\nu t}(0,0,\cos x_1)\): \(\int KD=0\) while the bandwidth
action diverges, balanced by \(\nu\log D\to-\infty\)); only the
one-sided comparison is a criterion.

**(d) (Position among critical actions.)** Pointwise on
\([0,T_{\max})\),
\[
K\,D\;\le\;\|u\|_{L^\infty}^2,
\qquad
K\,D\;\le\;C_S^2\,\|\nabla u\|_{L^3}^2,
\]
so hypothesis (b1) is implied by the finiteness of the classical Serrin
\((\infty,2)\) action \(\int\|u\|_\infty^2dt\) and by the critical
vorticity-class action \(\int\|\nabla u\|_{L^3}^2dt\); the criterion (b)
is therefore at least as strong as both.

**(e) (Osgood closure, conditional form.)** If, in addition, a.e. on
\([0,T_{\max})\),
\(K D\le\Phi(z)\,D+R\) with \(\Phi>0\) nondecreasing,
\(\int^\infty ds/\Phi=\infty\), \(R\ge0\),
\(\int_0^{T_{\max}}R\,dt<\infty\), then case (b1) holds. *(This clause is
a proven implication; Theorem O below shows its hypothesis cannot hold
uniformly in the \(R\equiv0\) form.)*

## Theorem O (unconditional static no-go; separated from the Main Theorem)

*This statement was a conditional Proposition in earlier versions of
this paper, resting on a capacity hypothesis (L\*). It is now
**unconditional**: the capacity bound is proven for a smoothly truncated
member of the coherent family, and the no-go argument quantifies over
**all** real zero-mean divergence-free trigonometric fields, so exhibiting
the bound for that one family is all it ever needed. See*
[lstar/lstar_proof_main.md](lstar/lstar_proof_main.md) *and the*
[audit](dependency_and_gap_audit.md) *§3.*

**The smoothly truncated coherent family.** Call
\(\chi\) **admissible** if
\[
\chi\in C^\infty([0,\infty);[0,1]),\qquad
\chi\equiv1\ \text{on }[0,\tfrac12],\qquad
\operatorname{supp}\chi\subset[0,1];
\]
no monotonicity is assumed. Write \(c_\chi=1+\|\chi'\|_{L^\infty}\).
For \(v_0\in\mathbb R^3\setminus\{0\}\) (no integrality is needed), an
admissible \(\chi\), and an integer \(N\ge2\), define
\[
\widehat u_N(k)=\chi\!\Big(\frac{|k|}N\Big)\frac{P_kv_0}{|k|^2}\quad(k\neq0),
\qquad\widehat u_N(0)=0,\qquad
P_k=I-\frac{k\otimes k}{|k|^2}.
\]
Since \(\chi\) vanishes on \([1,\infty)\), \(u_N\) is a finite
trigonometric polynomial banded at \(|k|\le N\).

Then (proven, exact — Theorem 1.3 of the capacity document): \(u_N\) is
real, zero-mean and exactly divergence-free, and with
\[
S_N^\chi=\sum_{k\neq0}\frac{\chi(|k|/N)^2}{|k|^2},\quad
T_N^\chi=\sum_{k\neq0}\frac{\chi(|k|/N)^2}{|k|^4},\quad
\Sigma_N^\chi=\sum_{k\neq0}\frac{\chi(|k|/N)}{|k|^2}
\]
(all finite sums) one has
\[
H_0=\tfrac23\|v_0\|^2T_N^\chi,\qquad
H_1=\tfrac23\|v_0\|^2S_N^\chi,\qquad
N_0^2=\frac{S_N^\chi}{T_N^\chi},\qquad
u_N(0)=\tfrac23\Sigma_N^\chi v_0,
\]
with the Bernstein lower bound
\(\|u_N\|_\infty^2/\|\nabla u_N\|_2^2\ge2(\Sigma_N^\chi)^2/(3S_N^\chi)\)
(not used below), and the proven two-sided lattice bounds, for
\(N\ge32\),
\[
\frac N{348}\le S_N^\chi\le128N,\qquad
6\le T_N^\chi\le128,\qquad
\frac N{44544}\le N_0^2(u_N)\le\frac{64}3N .
\]

**Capacity Theorem (proven; no hypothesis).** For every
\(v_0\in\mathbb R^3\setminus\{0\}\) and every admissible \(\chi\) there
exist \(c_0=c_0(\chi,v_0)>0\) and \(N_*=N_*(\chi,v_0)\) with
\[
\big\|\mathbb P(u_N\cdot\nabla u_N)\big\|_2^2\;\ge\;c_0N^3
\qquad\text{for all }N\ge N_* .
\]
This is Theorem 7.1(3) of
[lstar/lstar_proof_main.md](lstar/lstar_proof_main.md); the exponent
\(3\) is the sharp one (Lemma 7 caps
\(\|\mathbb P(u\cdot\nabla u)\|_2^2\le\|u\|_\infty^2H_1\asymp N^3\)).
Its constants are **not effective** — see "Constant dependence".

**Theorem O.** If \(\Phi:[0,\infty)\to(0,\infty)\) is nondecreasing and
\(K(u)\le\Phi(\log N_0^2(u))\) for every zero-mean divergence-free real
trigonometric field \(u\) on \(\mathbb T^3\), then \(\Phi(s)\ge c\,e^{s}\)
for all large \(s\), and hence \(\int^\infty ds/\Phi<\infty\): **no
Osgood-admissible \(\Phi\) satisfies a uniform pointwise bound** —
clause (e) of the Main Theorem cannot be activated through the
\(R\equiv0\) field-inequality route. (Solution-adapted remainders
\(R(t)\) are not excluded by Theorem O; see the audit §5.)

**What remains open, and is no longer used.** The *sharply* truncated
family \(\widehat u_N(k)=P_kv_0/|k|^2\) on \(1\le|k|\le N\) satisfies the
same exact laws with \(S_N,T_N,\Sigma_N=S_N\) in place of
\(S_N^\chi,T_N^\chi,\Sigma_N^\chi\), and saturates the Bernstein ratio
two-sidedly, \(\tfrac23S_N\le\|u_N\|_\infty^2/\|\nabla u_N\|_2^2\le S_N\).
The capacity bound \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\)
for *that* family — the paper's former **Hypothesis (L\*)** — is still
**OPEN**. It is now **unused**: nothing in this paper depends on it.
The obstruction is proven without it.

## Constant dependence

Constants are explicit **except for \(c_0\) and \(N_*\) in Theorem O's
Capacity Theorem, which are non-effective**. In detail:

* (a)–(c) contain no constants beyond \(\nu\); (d) uses the mean-zero
  Sobolev constant \(C_S\) of \(H^1(\mathbb T^3)\hookrightarrow L^6\).
* Theorem O's lattice constants are explicit and proven:
  \(N/348\le S_N^\chi\le128N\), \(6\le T_N^\chi\le128\),
  \(N/44544\le N_0^2\le\tfrac{64}3N\) for \(N\ge32\), giving dyadic gaps
  \(c_-\le s_{2^{20}N}-s_N\le c_+\) with \(c_-=\log(1048576/950272)
  =0.0985\ldots\) and \(c_+=\log(2^{20}\cdot950272)=27.63\ldots\).
  (The former sharp-family bounds \(N/250\le S_N\le432N\) for \(N\ge8\)
  remain true but are superseded and unused.)
* **Non-effective:** \(c_0\) and \(N_*\). The capacity proof obtains its
  test field \(\Psi\in C^\infty_{c,\sigma}(\mathbb R^3)\) from the
  *density* of \(C^\infty_{c,\sigma}\) in \(L^2_\sigma\) — a pure
  existence argument — so neither \(\Psi\), nor its support radius, nor
  \(I_\Psi=\int(V\cdot\nabla V)\cdot\Psi\), nor \(c_0\), nor \(N_*\) is
  exhibited. Theorem O's own conclusion constant \(c\) inherits this.
  An **effective** version would require: (i) an explicitly constructed
  divergence-free \(\Psi\) with \(I_\Psi\neq0\) — available in principle,
  since the non-vanishing is certified on the explicit set
  \(\{0<|\zeta|<\pi^2/2048,\ |\zeta_3|\ge|\zeta|/2,\
  |\zeta\times e_3|\ge|\zeta|/2\}\), but at a cost: Fourier support at
  scale \(\approx4.8\times10^{-3}\) forces spatial radius
  \(R\gtrsim10^3\), hence \(N_*\ge8R\gtrsim10^4\) and a \(c_0\) many
  orders below the measured \(\approx8.4\|v_0\|^4\); and (ii) an explicit
  lower bound on \(|I_\Psi|\) for that \(\Psi\), i.e. a quantitative form
  of the non-degeneracy \(\mathbb P(V\cdot\nabla V)\not\equiv0\).
  Neither is carried out.
* Numerically (no proof claimed, used nowhere): \(S_N=4\pi N-8.7\pm0.9\)
  and \(T_\infty=16.5323\ldots\) for the sharp family.

## What is not claimed

Global regularity for all data; finite-time blow-up for some datum; any
Clay statement. Clause (b) is a conditional criterion whose hypothesis
is not verifiable a priori; Theorem O shows one specific (Osgood-type)
route to verifying it is closed. Numerical values appear only in the
certificate appendix; they enter no proof — neither the Main Theorem nor
Theorem O, both of which are now unconditional (modulo F1–F3).
The honest headline is **"the static no-go is unconditional"**, *not*
"(L\*) is proven": the paper's literal Hypothesis (L\*), for the sharply
truncated family, remains open and is recorded as such above.
