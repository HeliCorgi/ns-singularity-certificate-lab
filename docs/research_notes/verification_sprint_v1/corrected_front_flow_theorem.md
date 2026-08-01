# Corrected front-flow time/rate statements (Verification Sprint V1, workstream A)

**Label: exact elementary calculus, self-contained; no PDE content is claimed.**
This corrects and sharpens the informal rate statements of the 2026-08-01
portfolio note §2. Everything below is a statement about a positive C¹
function \(N(s)\) and the clock change \(dt = N^{-2}ds\); nothing assumes an
orbit of the front flow exists.

## A.1 Exact integral formula

Let \(s\in[s_0,S_{\max})\), \(a(s):=\frac{d}{ds}\log N(s)\),
\(A(s):=\int_{s_0}^s a\), so \(N(s)=N(s_0)e^{A(s)}\) exactly, and
\(t(s)-t_0=\int_{s_0}^s N^{-2}d\sigma\). The orbit has a finite terminal
physical time \(T:=\lim_{s\to S_{\max}}t(s)<\infty\) **iff**
\(\int^{S_{\max}}e^{-2A}<\infty\), and then for every \(s\):

\[
\boxed{\;T-t(s)\;=\;\int_s^{S_{\max}}N(\sigma)^{-2}\,d\sigma
\;=\;N(s)^{-2}\int_s^{S_{\max}}e^{-2\,(A(\sigma)-A(s))}\,d\sigma\;}
\tag{A.1}
\]

Equivalently, along the orbit in physical time,
\(\frac{d}{dt}N^{-2}=-2a\) exactly, so
\(N^{-2}(t)=N^{-2}(t_0)-2\int_{t_0}^t a\,dt'\). No case assumption is used
in (A.1).

## A.2 The three regimes, separated

**(i) One-sided bound \(a\ge a_->0\) on \([s,S_{\max})\).**
Then \(A(\sigma)-A(s)\ge a_-(\sigma-s)\), hence \(S_{\max}=\infty\),
\(N\to\infty\), \(T<\infty\), and

\[
T-t(s)\;\le\;\frac{N(s)^{-2}}{2a_-}
\qquad\Longleftrightarrow\qquad
N(t)\;\ge\;\bigl(2a_-(T-t)\bigr)^{-1/2}.
\tag{A.2}
\]

*Only* this upper bound on \(T-t\) holds; there is **no** lower bound and
**no asymptotic rate**. One-sided consequences that survive in this regime:
\(\int^T N(t)^2\,dt=\infty\) (from (A.2); the BKM-side divergence needs only
\(a\ge a_-\)). Nothing about \(\int N\,dt\) can be asserted.

**(ii) Two-sided bound \(0<a_-\le a\le a_+<\infty\) on \([s,\infty)\).**
Sandwiching the integrand of (A.1):

\[
\frac{1}{2a_+}\;\le\;N(s)^2\,(T-t(s))\;\le\;\frac{1}{2a_-},
\tag{A.3}
\]

a **Type-I window**, not a rate: \(N(t)\sqrt{T-t}\in
[(2a_+)^{-1/2},(2a_-)^{-1/2}]\) with no limit claimed. Consequences:
\(\int^T N\,dt\le\sqrt{2/a_-}\,\sqrt{T-t}<\infty\) (dissipation-side
finiteness) and \(\int^T N^2dt=\infty\) (BKM side), both with explicit
constants from (A.3).

**(iii) Convergent rate \(a(s)\to a_\infty>0\).**
For every \(\varepsilon>0\) regime (ii) applies eventually with
\(a_\pm=a_\infty\pm\varepsilon\), so (A.1) gives the genuine asymptotic

\[
\lim_{s\to\infty}N(s)^2(T-t(s))=\frac1{2a_\infty},
\qquad
N(t)=\bigl(2a_\infty(T-t)\bigr)^{-1/2}\,(1+o(1)),\ t\uparrow T.
\tag{A.4}
\]

**The coefficient asymptotic (A.4) is asserted only in this regime.** All
previous informal uses of "\(N\asymp(2a(T-t))^{-1/2}\)" outside (iii) are
replaced by (A.2) or (A.3) as appropriate.

**Exactly periodic \(a\) (the DSS orbit case).** If \(a\) is exactly
\(S\)-periodic with mean \(\bar a=\frac1S\int_0^Sa>0\), then
\(N(s+S)=e^{\bar aS}N(s)\) (a \(\lambda\)-DSS orbit with
\(\lambda=e^{\bar aS}\)) and (A.1) shows

\[
\Phi(s):=N(s)^2\,(T-t(s))\ \text{is exactly \(S\)-periodic}
\tag{A.5}
\]

(no \(o(1)\) needed: \(A(\sigma+S)-A(s+S)=A(\sigma)-A(s)\)). The
log-periodic oscillation of \(N^2(T-t)\) is therefore an **exact** signature
of the periodic orbit, bounded inside the window (A.3) with
\(a_\pm=\max/\min a\).

## A.3 What the scaling group removes, and what it does not

The audited group \((\Psi,s,a,\nu)\mapsto(\mu\Psi,s/\mu,\mu a,\mu\nu)\) is
exact on the similarity system. For the **physical** problem (fixed torus,
fixed \(\nu\)) its correct reading is: amplitude \(\times\mu\) at fixed
\(\nu\) is equivalent to amplitude 1 at viscosity \(\nu/\mu\). Hence:

**Removed as obstacles (gauge freedoms):**
1. the absolute critical-energy constant \(c_E\) (any "required \(c_E\)"
   number is gauge, not an obstruction);
2. the absolute magnitude of per-stage gain versus viscous loss — only
   channel *ratios* (shape functionals) are invariant content;
3. the absolute blow-up rate constant \(a_\infty\) (rescalable; only its
   positivity is content).

**Not removed (the real remaining obstacles):**
1. **Lattice admissibility**: \(\mathbb Z^3\) admits only integer dilations;
   a periodic orbit needs \(\lambda=e^{\bar aS}\) compatible with the
   carrier lattice (doubling \(\lambda=2^m\); multi-type cycles need the
   super-period product in \(2^{\mathbb N}\), with per-stage ratios possibly
   irrational — audit item 10);
2. **Reality/Hermitian admissibility** of the closure class
   (\(\mathcal Q_\varphi R\) only, no global phase);
3. **DSS admissibility gates** (workstream B): the steady case is killed by
   Tsai 1998; any orbit whose per-period critical-norm oscillation vanishes
   under refinement is killed by the same gate (pre-registered);
4. **Existence and transverse stability of the orbit** (pilot A: no
   attracting orbit in the scanned lattice box; shadowing needs a
   contraction margin exceeding the lattice consistency error — workstream
   E);
5. **Shape-closure sign** \(\chi_{\rm shape}>0\) at the recurrent profile —
   the gauge argument reduces closure to this sign *plus* the following
   uniformity;
6. **Inviscid-limit uniformity**: the gauge trades large \(c_E\) for small
   effective viscosity, so "closure for some large \(c_E\)" requires the
   shape problem to persist in the \(\nu_{\rm eff}\to0\) limit — a genuine
   analytic obligation, not a gauge artifact (this is the corrected form of
   audit error 4);
7. **Seed/basin entry** (PO-09) and smoothness of the assembled initial
   datum.

## A.4 Corrected-statement ledger

| old informal claim | status | corrected statement |
|---|---|---|
| \(N\asymp(2a(T-t))^{-1/2}\) whenever \(a\ge a_->0\) | **withdrawn** | (A.2): one-sided lower bound on \(N\) only |
| \(N^2(T-t)\to\) const on bounded orbits | **withdrawn** | (A.3): window only; limit needs (iii) |
| coefficient \(1/(2a_\infty)\) | kept, regime-restricted | (A.4), only for \(a\to a_\infty\) |
| log-periodic signature is asymptotic | **strengthened** | (A.5): exactly periodic, no error term |
| "\(c_E=228\) closes the budget" | **withdrawn** (audit err. 4/5) | gauge statement; real obligation is A.3-item 5+6 |
