# Theorem statement

**A bandwidth–dissipation dichotomy for the three-dimensional periodic
Navier–Stokes equations, with an exact action representation**

This file states the single main theorem and the one clearly-separated
conditional proposition. Proofs: [complete_proof.md](complete_proof.md).
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
a proven implication; the following proposition shows its hypothesis
cannot hold uniformly.)*

## Proposition (conditional static no-go; separated from the Main Theorem)

Let \(v_0\in\mathbb Z^3\setminus\{0\}\) and define the coherent
critical-spectrum family
\[
\widehat u_N(k)=\frac{P_kv_0}{|k|^2},\quad 1\le|k|\le N,\qquad
P_k=I-\frac{k\otimes k}{|k|^2}.
\]
Then (proven, exact): \(u_N\) is real and divergence-free;
\(H_0=\tfrac23\|v_0\|^2T_N\), \(H_1=\tfrac23\|v_0\|^2S_N\),
\(N_0^2=S_N/T_N\asymp N\), \(u_N(0)=\tfrac23S_Nv_0\), and the
Bernstein ratio is saturated two-sidedly:
\(\tfrac23S_N\le\|u_N\|_\infty^2/\|\nabla u_N\|_2^2\le S_N\)
(lower bound from the exact point value; upper bound from
\(\|u\|_\infty\le\sum|\hat u_k|\le\sqrt{S_N}\,\|\nabla u\|_2\)),
where \(S_N=\sum_{1\le|k|\le N}|k|^{-2}\asymp N\) (crude explicit
constants proven; numerically \(4\pi N-8.7\pm0.9\)),
\(T_N=\sum|k|^{-4}\nearrow T_\infty<\infty\).

**Hypothesis (L\*)** *(open; exactly certified at \(N\le8\), measured to
\(N=32\); see the certificate appendix)*: there is \(c_0>0\) with
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) for all large \(N\).

**Proposition.** Assume (L\*). If \(\Phi\) is nondecreasing and
\(K(u)\le\Phi(\log N_0^2(u))\) for every zero-mean divergence-free real
trigonometric field \(u\), then \(\Phi(s)\ge c\,e^{s}\) for all large
\(s\), and hence \(\int^\infty ds/\Phi<\infty\): **no
Osgood-admissible \(\Phi\) satisfies a uniform pointwise bound** —
clause (e) of the Main Theorem cannot be activated through the
\(R\equiv0\) field-inequality route. (Solution-adapted remainders
\(R(t)\) are not excluded by this proposition; see the audit §5.)

## Constant dependence

All constants are explicit: (a)–(c) contain no constants beyond \(\nu\);
(d) uses the mean-zero Sobolev constant \(C_S\) of
\(H^1(\mathbb T^3)\hookrightarrow L^6\); the Proposition's constants are
\(c_0\) (from L\*) and the crude proven lattice bounds
\(N/250\le S_N\le432N\) (\(N\ge8\)) with dyadic gap bounds
\(c_-\le s_{2N}-s_N\le c_+\) (numerically \(S_N=4\pi N-8.7\pm0.9\) and
\(T_\infty=16.5323\ldots\); not used in proofs).

## What is not claimed

Global regularity for all data; finite-time blow-up for some datum; any
Clay statement. Clause (b) is a conditional criterion whose hypothesis
is not verifiable a priori; the Proposition shows one specific
(Osgood-type) route to verifying it is closed. Numerical values appear
only in the certificate appendix and support only Hypothesis (L\*),
never the Main Theorem.
