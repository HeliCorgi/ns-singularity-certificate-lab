# Referee report A (adversarial)

Documents refereed: `theorem_statement.md`, `complete_proof.md`,
`paper_draft.tex`; dependency structure per `dependency_and_gap_audit.md`.
Numerical cross-checks against
`outputs/verification_sprint_v1/osgood_gate/{exact_family_certificates.json,summary.json}`
and against an **independent re-implementation** of the field, the
advection term, and the Leray projection written for this review.

Findings are graded CRITICAL (invalidates a stated claim) / MAJOR (gap
needing repair) / MINOR (wording, precision, bookkeeping). Correct steps
are not commented on except where a near-miss is worth recording.

---

## CRITICAL

### C1. Theorem (c), the sentence beginning "In particular ... is equivalent" is false

**Location.** `theorem_statement.md` lines 73–76:

> "In particular the regularity criterion of (b1) is equivalent (given the
> energy equality) to the a-priori finiteness of the \(\dot H^1\)-bandwidth
> action \(\int(\|\partial_tu\|_2^2/D+\nu^2N_1^2)\,dt\)."

Same defect, restated, in `paper_draft.tex` §"Position among known
criteria; limitations": *"By (c), it is exactly the a-priori control of
\(\int(\|\partial_tu\|_2^2/D+\nu^2N_1^2)dt\) ... plus a logarithm"*.

**Why it is false.** Only one implication survives. The identity (c) is
\[
\int_0^{T'}KD=A(T')+\nu\log\frac{D(T')}{D(0)},\qquad
A(T')=\int_0^{T'}\Bigl(\frac{\|\partial_tu\|_2^2}{D}+\nu^2N_1^2\Bigr)dt .
\]
The logarithm is *not* sign-definite, and for global solutions
\(D(T')\downarrow0\), so \(\nu\log(D(T')/D(0))\to-\infty\) and \(A\) may
diverge while \(\int KD\) stays finite.

**Explicit counterexample (an exact, global, smooth NS solution on
\(\mathbb T^3\)).** Take
\(u(x,t)=e^{-\nu t}\,(0,0,\cos x_1)\). Then \(\nabla\cdot u=0\),
\((u\cdot\nabla)u\equiv0\) (so \(\mathcal N\equiv0\)),
\(\Delta u=-u\), and \(\partial_tu=-\nu u=\nu\Delta u\): a genuine
solution with \(u_0\in H^m_\sigma\) for every \(m\), \(u_0\neq0\),
\(T_{\max}=\infty\).
* \(K\equiv0\), so \(\int_0^\infty KD\,dt=0<\infty\): **(b1) holds.**
* \(\partial_tu=-\nu Au\), \(H_2=H_1=H_0\), hence
  \(\|\partial_tu\|_2^2/D=\nu^2\) and \(\nu^2N_1^2=\nu^2\), so
  \(A(T')=2\nu^2T'\to\infty\): **the \(\dot H^1\)-bandwidth action
  diverges.**
* Consistency check of the identity itself:
  \(2\nu^2T'+\nu\log e^{-2\nu T'}=0=\int_0^{T'}KD\). The *identity* (c)
  is correct; only the "in particular" is wrong.

**Minimal fix.** Replace the equivalence by the (true, and sharper than
what is currently claimed) one-sided statement:

> \(\nu\log\dfrac{D(T')}{D(0)}=\displaystyle\int_0^{T'}\frac{2\nu\langle\partial_tu,Au\rangle}{D}
> \le\int_0^{T'}\Bigl(\frac{\|\partial_tu\|_2^2}{D}+\nu^2\frac{\|Au\|_2^2}{D}\Bigr)dt=A(T')\)
> by Cauchy–Schwarz and \(2ab\le a^2+b^2\); hence
> \(\int_0^{T'}KD\,dt\le 2A(T')\). **Finiteness of the \(\dot
> H^1\)-bandwidth action is therefore sufficient for (b1); it is not
> necessary** (example above).

The abstract's "which localises the criterion at the \(\dot H^1\)-bandwidth
action" survives only in this one-sided reading and should be reworded.

---

## MAJOR

### M1. Lemma 0, positivity of \(H_0\): the Grönwall/open-closed argument is circular

**Location.** `complete_proof.md` lines 43–49.

Two independent defects in the same paragraph.

1. **The forward-uniqueness clause proves nothing.** "if \(H_0(t_1)=0\)
   then \(u(t_1)=0\) and by forward uniqueness (F1) \(u\equiv0\) on
   \([t_1,\infty)\)" is a true statement that yields no contradiction with
   \(u_0\neq0\). Ruling out \(H_0(t_1)=0\) by this route would require
   *backward* uniqueness, which the audit (§2, and §"does not use ...
   backward uniqueness") explicitly declines to import.
2. **The Grönwall bound is vacuous exactly in the case to be excluded.**
   \(H_0(t)\ge H_0(0)\exp(-2\nu\int_0^tN_0^2)\) is useless if
   \(\int_0^{t_1}N_0^2=+\infty\), and \(N_0^2=H_1/H_0\to\infty\) is
   precisely what happens as \(H_0\downarrow0\). Consequently the assertion
   "the set \(\{H_0>0\}\) is open and closed in \([0,T_{\max})\) by this
   Grönwall bound" is unjustified: openness is clear from continuity;
   **closedness is not established**.

Note also that interpolation cannot repair this: \(H_1\le
H_0^{1-1/m}\|u\|_{\dot H^m}^{2/m}\) gives \(\dot H_0\ge-CH_0^{1-1/m}\),
an ODE whose solutions *do* reach \(0\) in finite time. So the gap is real,
not cosmetic.

**Minimal fix (bootstrap, no new external input).** Let
\(t_1=\inf\{t:H_0(t)=0\}\) and suppose \(t_1<T_{\max}\). On \([0,t_1)\) we
have \(H_0>0\), so Lemmas 1–4 and 7 are available there. Put
\(M=\sup_{[0,t_1]}\|u\|_{H^m}<\infty\) (\(u\in C([0,T_{\max});H^m)\),
\([0,t_1]\) compact). By Lemma 7, \(KD\le\|u\|_\infty^2\le C_mM^2\) on
\([0,t_1)\); by Lemma 4,
\(z(t)\le z(0)+C_mM^2t_1/(2\nu)\) there, i.e. \(N_0^2\) is **bounded** on
\([0,t_1)\); hence \(\int_0^{t_1}N_0^2<\infty\) and the Grönwall bound now
gives \(H_0(t_1)\ge H_0(0)e^{-2\nu\int_0^{t_1}N_0^2}>0\), a contradiction.

**Consequence for the audit.** With this fix the positivity half of L0 uses
L1–L4 and L7, so the audit's claim (§1) "Every arrow is used exactly as
drawn; no lemma uses a later one" is false and the dependency graph must be
redrawn (the standard fix is to split L0 into L0a = regularity/derivatives
and L0b = positivity, and place L0b after L7).

### M2. Theorem (b): "Exactly one of the following holds" is not proven, and the alternatives are not exhaustive as written

**Location.** `theorem_statement.md` lines 60–64; `paper_draft.tex`
Theorem (b) and abstract; `complete_proof.md` Lemma 5.

As stated, alternative 1 is "\(\int_0^{T_{\max}}KD<\infty\)" and
alternative 2 is "\(T_{\max}<\infty\) **and** \(\int KD=\infty\)". The
configuration \(T_{\max}=\infty\) **and** \(\int_0^\infty KD=\infty\)
satisfies neither, so "exactly one" is a strictly stronger claim than what
Lemma 5 proves. Lemma 5 proves only (i) \(\int KD<\infty\Rightarrow
T_{\max}=\infty\) and (ii) exclusivity. Exhaustiveness is asserted, never
argued.

(It is probably true — on \(\mathbb T^3\) Poincaré gives \(H_1\ge H_0\),
hence \(\dot H_0=-2\nu H_1\le-2\nu H_0\) and \(H_0\le H_0(0)e^{-2\nu t}\),
from which a global solution decays and \(\int KD<\infty\) — but that
argument appears nowhere and would need the \(H^1\) and \(H^2\) decay too.)

**Minimal fix (recommended).** State the dichotomy in the form that is
actually proven:

> Exactly one of: (b1) \(\int_0^{T_{\max}}KD\,dt<\infty\), and then \(z\) is
> bounded, \(u\in L^\infty(0,T_{\max};H^1)\) and \(T_{\max}=\infty\);
> (b2) \(\int_0^{T_{\max}}KD\,dt=\infty\).
> In particular \(T_{\max}<\infty\Rightarrow\int_0^{T_{\max}}KD\,dt=\infty\).

Alternatively keep the present (b2) and downgrade "Exactly one" to "At most
one, and \(T_{\max}<\infty\) implies (b2)".

### M3. Lemma 11 is not a proof; \(C_0\) is claimed "explicit" and never given

**Location.** `complete_proof.md` lines 216–227.

The displayed chain contains two literal ellipses *inside formulas*:
`\(\ge|k|^{-2}(1+\tfrac{\sqrt3}2\cdot\tfrac2{|k|}\cdot\ldots)^{-1}\)` and
`\(\int_{1\le|x|\le N\pm\sqrt3/2}|x|^{-2}dx=4\pi(N\mp\ldots)\)`. In a
document titled *Complete proof*, and for a lemma whose constants are
advertised in `theorem_statement.md` §"Constant dependence" as *explicit*,
this is a placeholder, not an argument. No value of \(C_0\) is supplied
anywhere in the three documents.

Two further problems with the sketch as given:
* \(f(x)=|x|^{-2}\) is **not convex** on \(\mathbb R^3\setminus\{0\}\)
  (Hessian eigenvalues \(6|x|^{-4}\) radially, \(-2|x|^{-4}\)
  tangentially), so "the cube average dominates the centre value" is not
  available; the direction of the cube comparison must be argued, not
  asserted.
* Shell increments are not monotone, so a naive per-shell comparison fails:
  \(S_7-S_6=11.644<4\pi=12.566\) while \(S_6-S_5=13.483>4\pi\).

**Status of the statement itself (verified numerically, \(2\le N\le200\)):
it is true, with large margin.** I computed \(S_N-4\pi N\in[-9.57,-7.90]\)
throughout, tending to \(\approx-9.0\); the claimed lower bound
\(4\pi(N-2)=4\pi N-25.13\) therefore has \(\approx16\) of slack, and the
upper bound holds with \(C_0=0\).

**Proposed minimal fix.** Replace the sketch by Abel summation against the
lattice count. With \(A_M=\#\{k:1\le|k|^2\le M\}\) and the elementary
boundary-cube bound \(|A_M-\tfrac{4\pi}3M^{3/2}|\le CM\) (explicit \(C\)),
\[
S_N=\frac{A_{N^2}}{N^2}+\sum_{M=1}^{N^2-1}\frac{A_M}{M(M+1)}
=4\pi N+O(\log N),
\]
which combined with a finite exact check for \(N\le N_0\) yields
\(4\pi(N-2)\le S_N\le 4\pi N\) with explicit constants. State \(C_0=0\).
Similarly, \(T_\infty\le16.6\) is asserted with no derivation; the true
value is \(T_\infty=16.5323\ldots\) (computed here to \(N=200\) with tail
\(4\pi/N\)), so the bound is correct but must be justified (same
comparison, with \(|x|^{-4}\)).

### M4. Overclaim in Theorem (a): "equality cases characterised"

**Location.** `theorem_statement.md` line 58.

`complete_proof.md` Lemma 4 exhibits the defect
\(\Gamma^{\rm CS}+\Gamma^{\rm SC}\) and proves each part nonnegative. It
nowhere characterises when they vanish. (For the record the
characterisation is easy and should simply be added:
\(\Gamma^{\rm SC}=0\iff\sqrt V=\tfrac1{2\nu}\sqrt{\|\mathcal N\|^2/H_0}\);
\(\Gamma^{\rm CS}=0\iff\) there is \(\lambda\ge0\) with
\(\hat{\mathcal N}_k=\lambda\,\mathrm{sgn}(x_k-\mu)\hat u_k\) for every
\(k\) with \(\varepsilon_k>0\), i.e. simultaneous saturation of the modal
and the vector Cauchy–Schwarz.) Either prove it or delete the parenthesis.

### M5. `paper_draft.tex` drifts beyond the `.md` statements in three places

1. **Abstract, lines 48–50** and **Introduction, line 68**: "conditionally
   on a single **certified** lattice-sum lower bound" / "conditionally on
   one explicitly **certified** lattice estimate". Hypothesis (L\*) is
   *open* (audit §3, §5). Only finitely many exact *values* are certified;
   an asymptotic inequality \(\ge c_0N^3\) cannot be certified at finitely
   many \(N\). Fix: "conditionally on a single unproven lattice-sum lower
   bound, for which exact finite-\(N\) evidence is given in the appendix".
2. **Introduction, line 74**: "the dichotomy is a **genuine
   strengthening** of the classical criteria". Only \(KD\le\|u\|_\infty^2\)
   and \(KD\le C_S^2\|\nabla u\|_{L^3}^2\) are proven — i.e. *at least as
   strong*, exactly as the `.md` correctly says (line 88). Strictness is
   never demonstrated (no example with \(\int KD<\infty\) and both
   classical actions infinite). Fix: "at least as strong"; or add a
   separating example.
3. **Introduction, line 71**: "**saturates** the Bernstein ratio". Lemma 10
   proves only the lower bound \(\|u_N\|_\infty^2/H_1\ge\tfrac23S_N\). Fix:
   add the one-line matching upper bound, valid for any field spectrally
   supported in \(B_N\),
   \(\|u\|_\infty\le\sum_{k\in B_N}|\hat u_k|
   \le(\sum|k|^2|\hat u_k|^2)^{1/2}(\sum_{k\in B_N}|k|^{-2})^{1/2}
   =\sqrt{H_1S_N}\), so the ratio lies in \([\tfrac23,1]\cdot S_N\); then
   "saturates up to the factor \(2/3\)" is licensed. (In fact the
   certificates show \(\|u_N\|_\infty^2/D=\tfrac23S_N\) **exactly**, since
   the sup is attained at \(x=0\); see the numerical section below. Saying
   so would be stronger *and* correct.)

### M6. Scope overclaim: the Proposition does not close clause (e)

**Location.** `theorem_statement.md` lines 95–96 ("the following
proposition shows its hypothesis cannot hold uniformly") and 121–123
("clause (e) of the Main Theorem **cannot be activated by any
field-inequality route**").

Clause (e) assumes \(KD\le\Phi(z)D+R\) **a.e. along the solution**, with
\(R\ge0\), \(\int R<\infty\). The Proposition excludes only the strictly
narrower class of **uniform, \(R=0\), all-fields** inequalities
\(K(u)\le\Phi(\log N_0^2(u))\). It does not exclude, e.g.,
\(K\le\Phi(z)+R/D\) with integrable \(R\), nor any trajectory-dependent
bound. Fix: adopt the `.tex` wording ("cannot be activated by any *uniform
pointwise field* inequality") in the `.md` as well, and say explicitly that
the \(R\)-term route and trajectory-dependent routes remain open.

### M7. The audit's "no other gaps" claim is falsified by C1, M1, M2, M3

**Location.** `dependency_and_gap_audit.md` §5: "There are no other gaps:
the Main Theorem's proof chain ... is complete at the stated regularity",
and §4.2, which repeats the defective positivity argument verbatim
("proven impossible for \(u_0\neq0\) — L0, via the Grönwall two-sidedness
of \(H_0\) and forward uniqueness"). Both must be revised; the gap table
should gain rows for M1 and M3 (M2 and C1 are statement repairs, not gaps).

### M8. Reproducibility trap: half- vs full-energy fields in the cited JSON

`outputs/.../summary.json` stores `H0`, `H1`, `H2` as **half** energies
(`run_osgood_gate.measure()` multiplies by \(0.5\)) while
`exact_family_certificates.json` stores them **full**; `K` in both files is
built from `nonlinear_sq` and `grad_sq`, which are full. A reader
recomputing \(K=\)`nonlinear_sq`\(/\)`H1`\(^2\) from `summary.json`
obtains \(4\times\) the paper's value (e.g. \(3.1536\) instead of
\(0.78841\) at \(N=4\)) — verified. Moreover the certificate's own note
*"K is convention free because measure() already divides by grad_sq^2"* is
**false as a general statement**: \(K=\|\mathcal N\|_2^2/D^2\) picks up a
factor \(2\) under a uniform half-energy convention
(\(K_{1/2}=(\tfrac12\|\mathcal N\|^2)/(\tfrac12D)^2=2K\)), and a factor
\((2\pi)^{-3}\) under Lebesgue instead of normalised measure. \(K\) *is*
invariant under amplitude scaling \(u\mapsto\lambda u\), and scales as
\(p^{-2}\) under \(k\mapsto pk\). Fix: state in the appendix that all
quoted \(K\) use full-lattice, normalised-measure quantities, and correct
the JSON note.

---

## MINOR

1. `complete_proof.md` line 39: garbled text inside the proof —
   "\((r\le1\le m-\ldots\) indeed \(Au\in C(H^{m-2})\subset C(L^2)\))".
   The substance is fine (\(m>5/2\Rightarrow Au\in C(H^{m-2})\subset
   C(L^2)\), and \((u(t+h)-u(t))/h\to\partial_tu\) in \(H^{m-2}\subset
   L^2\), so both factors converge in \(L^2\) and the pairing passes to the
   limit); rewrite the parenthesis.
2. `complete_proof.md` line 62: "by Parseval (definition of \(a_k\))" —
   should be \(\alpha_k\).
3. `theorem_statement.md` line 28: \(z=\log N_0^2\ (\ge0)\). The
   nonnegativity is Poincaré on \(\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3\)
   (\(|k|\ge1\Rightarrow H_1\ge H_0\)); it is used in Lemma 8
   (\(\Omega(z)=\int_0^z\)) and should be stated once.
4. `complete_proof.md` line 119: "which is therefore rational for rational
   data". The defect total is
   \(\frac1\mu[-2\mathrm{Cov}+2\nu V+\|\mathcal N\|^2/(2\nu H_0)]\); it
   contains \(1/(2\nu)\), so **\(\nu\) must also be rational**. Add that,
   and say "rational data" means a finite trigonometric field with rational
   coefficients.
5. `theorem_statement.md` line 63 / Lemma 5: "(the solution is global and
   smooth for \(t>0\))". Smoothing for \(t>0\) is not among F1–F3. Either
   add a parabolic-smoothing input F4 or drop "and smooth".
6. Lemma 8: "\(\Omega(z)=\int_0^zds/\Phi\) is \(C^1\)" requires \(\Phi\)
   continuous; a nondecreasing \(\Phi\) may jump. \(\Omega\) is locally
   Lipschitz and, with \(z\in C^1\), \((\Omega\circ z)'=\Omega'(z)z'\)
   a.e., which is enough — but say so, or add continuity to the hypothesis
   of (e). Also the final sentence "Lemma 5's argument gives (b1)" skips the
   one line that actually delivers *(b1)*, namely
   \(\int KD\le\Phi(z_*)\int D+\int R\le\Phi(z_*)H_0(0)/(2\nu)+\int R<\infty\).
7. Proposition proof: the interpolation over \([s_N,s_{N+1}]\) tacitly uses
   that \(s_N\) is **nondecreasing** and \(s_N\to\infty\), so the intervals
   tile \([s_{N_1},\infty)\). True (\(S_{N+1}T_N>S_NT_{N+1}\) because
   \(\Delta S\cdot T_N\ge6\Delta S>4\pi\Delta S/N\ge S_N\Delta T\) for
   \(N\ge3\)), but unstated. The factor argument itself is correct:
   \(e^{s_N}\ge e^{s}e^{-(s_{N+1}-s_N)}\ge\tfrac12e^{s}\) once
   \(s_{N+1}-s_N\le\log2\), which holds already from \(N=2\) (largest gap
   \(0.419\)).
8. `complete_proof.md` line 249 and `paper_draft.tex` Remark: "a
   near-constant \(91.2\%\)–\(92.1\%\)". The measured minimum is
   \(0.9117606\) at \(N=32\), i.e. \(91.18\%\). Write \(91.1\%\)–\(92.1\%\)
   or quote the five values.
9. `theorem_statement.md` line 100: \(v_0\in\mathbb Z^3\setminus\{0\}\) —
   integrality is never used; any \(v_0\in\mathbb R^3\setminus\{0\}\) works.
10. Notation clash: \(u_N(0)\) denotes the field evaluated at \(x=0\)
    (Lemma 10) but reads as an initial datum everywhere else. Use
    \(u_N(x{=}0)\) or \(u_N|_{x=0}\).
11. `dependency_and_gap_audit.md` §3 and `paper_draft.tex` Remark:
    \(\|(v_0\cdot\nabla)u_N\|_2^2\le\tfrac23\|v_0\|^4S_N\). Lemma 9 gives
    the sharper \(\le\tfrac13\|v_0\|^4S_N\)
    (\(\sum(v_0\cdot k)^2|k|^{-4}=\tfrac13\|v_0\|^2S_N\)); the true value at
    \(N=4\) is \(1089.27\) vs \(\tfrac13\cdot196\cdot40.698=2658.9\). The
    stated bound is valid but loose by exactly \(2\) — a factor that looks
    like a residue of the half/full convention the audit says was
    eliminated. Recommend replacing \(\tfrac23\) by \(\tfrac13\).
12. `theorem_statement.md` line 87 / Lemma 7: "critical **vorticity**-class
    action \(\int\|\nabla u\|_{L^3}^2dt\)" — \(\nabla u\) and \(\omega\)
    are comparable only via Calderón–Zygmund, which is not in F1–F3. Either
    add it or call the class "\(\int\|\nabla u\|_{L^3}^2\)" without
    "vorticity". Also \(C_S\) must be declared as the constant for the
    **normalised** measure (the paper's convention); the standard
    \(H^1(\mathbb T^3)\hookrightarrow L^6\) constants in the literature are
    for Lebesgue measure.
13. `dependency_and_gap_audit.md` §4.4 refers to "the finite-Fourier
    identities (I.1)–(I.4)"; only (I.1) and (I.2) are labelled in
    `complete_proof.md`. Dangling reference.
14. `dependency_and_gap_audit.md` §1 graph: L6 is drawn hanging off L1 but
    uses only L0; L8 uses L1 (the \(\int D\) bound) and L5 as well as L4;
    L7 uses neither L0 nor L1. With M1's fix the graph changes further.
15. `theorem_statement.md` line 36: "\(K\) is a scale- and
    amplitude-covariant ... functional". Be precise: \(K\) is **invariant**
    under \(u\mapsto\lambda u\) and satisfies \(K\mapsto p^{-2}K\),
    \(N_0^2\mapsto p^2N_0^2\) under \(u(x)\mapsto u(px)\) — worth stating,
    because it shows the obstruction family beats lattice dilation by a
    factor \(N_0^4\) and is therefore not a scaling artefact.
16. `theorem_statement.md` line 114: "(open; exactly certified at
    \(N\le8\), ...)". The exact certificates exist at \(N=4,6,8\) only
    (`preregistration.exact_bands=[4,6,8]`), not for all \(N\le8\).
17. Lemma 5's phrasing "by F2 the solution is regular on \([0,T_{\max}]\)
    and extends, so \(T_{\max}=\infty\)" writes a closed interval before
    \(T_{\max}<\infty\) has been assumed. Make the contradiction explicit.

---

## Numerical verification performed (independent of the repository code)

I rebuilt \(\hat u_N(k)=P_kv_0/|k|^2\), the advection term, and the Leray
projection from scratch and recomputed everything at two grid resolutions
per band. Convention used: full lattice \(k\in\mathbb Z^3\setminus\{0\}\),
normalised measure, \(H_r=\sum|k|^{2r}|\hat u_k|^2\),
\(K=\|\mathcal N\|_2^2/D^2\) — i.e. the paper's convention.

| quantity | paper / tex | repository JSON | this review | verdict |
|---|---|---|---|---|
| \(K(u_4)\) | \(0.7884107043\ldots\) | \(0.7884107043392768\) | \(0.7884107043392770\) (M=32), \(\ldots764\) (M=48) | correct |
| \(K(u_6)\) | \(1.4344718\ldots\) | \(1.4344718079662586\) | \(1.4344718079662584\) | correct |
| \(K(u_8)\) | \(2.0372430\ldots\) | \(2.0372430105569106\) | \(2.0372430105569124/…128\) | correct |
| exact-vs-float agreement | \(\le1.3\times10^{-15}\) | \(0.0,\,4.64\!\times\!10^{-16},\,1.31\!\times\!10^{-15}\) | — | correct |
| \(K/N_0^2\), \(N=4\to32\) | \(0.259\to0.396\) | \(0.25867\to0.39619\) | — | correct |
| \(Q/N^2\), \(N=8\to32\) | \(850.7\to1029.5\) | \(850.713\to1029.469\) | — | correct |
| Leray retention | \(0.780\to0.867\) | \(0.78018\to0.86735\) | — | correct |
| pairing captured fraction | "\(91.2\%\)–\(92.1\%\)" | \(0.92104,0.91693,0.91509,0.91265,0.91176\) | — | **low end is \(91.18\%\)** (MINOR 8) |
| band point counts | — | \(256,924,2108\) | \(256,924,2108\) | correct |

Family laws, checked **by hand in exact rational arithmetic at \(N=2\)**
(32 lattice points, \(v_0=(1,2,3)\)):
\(S_2=97/6\), \(T_2=739/72\);
\(H_0=5173/54=\tfrac23\cdot14\cdot T_2\) ✔;
\(H_1=1358/9=\tfrac23\cdot14\cdot S_2\) ✔;
\(u_2(x{=}0)=(97/9,194/9,97/3)=\tfrac23S_2v_0\) ✔;
\(N_0^2=S_2/T_2=1.5751\) ✔;
\(|u_2(0)|^2/H_1=97/9=\tfrac23S_2\) ✔ (equality, not merely \(\ge\)).
Reconfirmed numerically at \(N=4,6,8,16,32\); in particular
\(\|u_N\|_\infty^2/D=\tfrac23S_N\) **exactly** (\(27.13214\), \(44.69651\),
\(60.70872\), \(127.89161\), \(262.01400\)), so the "Bernstein-ratio
saturation" claim can be upgraded from \(\ge\) to \(=\) for this family.

Lattice sums, computed over \(|k|\le200\) (\(3.35\times10^7\) points):
\(S_N-4\pi N\in[-9.57,-7.90]\) for \(2\le N\le200\), so
\(4\pi(N-2)\le S_N\le4\pi N\) holds throughout with \(\approx16\) slack and
\(C_0=0\) is admissible; \(T_N\uparrow\), \(T_{200}=16.4695\),
\(T_\infty=16.5323\ldots\) (so "\(T_\infty\le16.6\)" is correct);
\(N_0^2/N\to0.7603\), confirming \(N_0^2\asymp N\).

Evidence bearing on (L\*): \(\|\mathcal N\|_2^2/N^3\) increases
\(1777\to2600\to2874\to3580\to3962\) over \(N=4,6,8,16,32\) (local exponent
\(3.15\) from \(N=16\) to \(32\)), so (L\*) is consistent with the data and
conservative. For orientation, \(\|\mathcal N\|_2^2/(\|u\|_\infty^2D)\)
\(=0.0335\to0.0368\) over \(N=8\to32\), i.e. the family sits a fixed factor
\(\approx27\) below the trivial \(N^3\)-order ceiling of Lemma 7 — the
right qualitative picture for (L\*).

---

## Focus-list items checked and found free of defect

Recorded for completeness only (no comment implies no finding):

* **L0 difference quotients.** \(H_r(t+h)-H_r(t)=\langle
  A^r(u(t+h)+u(t)),u(t+h)-u(t)\rangle\) by self-adjointness; \(Au\in
  C(H^{m-2})\subset C(L^2)\) for \(m>5/2\), and the quotient converges in
  \(H^{m-2}\subset L^2\); both factors converge in \(L^2\), so the pairing
  passes to the limit. Continuity of \(H_2\) needs only \(m\ge2\). Sound.
* **L1 \(T_0=0\).** \(\mathbb Pu=u\) and \(\mathbb P\) self-adjoint give
  \(\langle u,\mathbb P(u\cdot\nabla u)\rangle=\langle
  u,(u\cdot\nabla)u\rangle=\tfrac12\int u\cdot\nabla|u|^2=0\); \(m>5/2\)
  gives \(u\in C^{1,\alpha}\), so every manipulation is classical. Energy
  direction and \(\int_0^{T_{\max}}D\le H_0(0)/(2\nu)\) sound.
* **L2 algebra.** \(V=H_2/H_0-\mu^2\Rightarrow H_0H_2-H_1^2=H_0^2V\);
  first bracket \(=(T_1-\mu T_0)/H_1=\mathrm{Cov}/\mu\); second
  \(=V/\mu\). Both verified symbolically.
* **L3 summability.** \(\sum(x_k-\mu)^2\varepsilon_k=H_2-2\mu
  H_1+\mu^2H_0<\infty\) (needs \(m\ge2\)); \(\sum\eta_k=\|\mathcal
  N\|_2^2<\infty\) (needs \(u\cdot\nabla u\in L^2\), true for \(m>5/2\));
  the intermediate sum \(\sum|x_k-\mu|\sqrt{\varepsilon_k\eta_k}\)
  converges by the same Cauchy–Schwarz. Splitting is
  \((|x_k-\mu|\sqrt{\varepsilon_k})\cdot(\sqrt{\eta_k})\) — correct, no
  factor-2 slip in the full-lattice convention.
* **L4 square completion.** \(\sup_{s\ge0}(\alpha s-\nu
  s^2)=\alpha^2/(4\nu)\) at \(s=\alpha/(2\nu)\); the two-part defect sums
  to \(KD/(2\nu)-\frac{d}{dt}\log N_0^2\) with the radicals cancelling
  exactly (verified term by term). Degenerate case \(V=0\) consistent
  (\(\mathrm{Cov}=0\), \(\Gamma^{\rm CS}=0\), \(\Gamma^{\rm SC}=KD/2\nu\)).
* **L5 Serrin exponent.** \((p,q)=(6,\infty)\) gives \(2/q+3/p=1/2<1\):
  genuinely subcritical, not endpoint. \(H_0\) non-increasing is available
  from L1.
* **L6.** \(D>0\) from L0; all three quotients continuous on the compact
  \([0,T']\); \(2\langle\partial_tu,Au\rangle=\dot H_1=\dot D\);
  \(\|Au\|_2^2/D=H_2/H_1=N_1^2\); log integration legitimate
  (\(D\in C^1\), \(D>0\)). The **identity** is correct (independently
  confirmed on the exact solution of C1).
* **L7.** Hölder \(1/2=1/6+1/3\) correct; \(\mathbb P\) orthogonal on
  \(L^2\) so \(\|\mathcal N\|_2\le\|(u\cdot\nabla)u\|_2\).
* **L9.** Sign-flip and permutation invariance of \(B_N\) argument is
  correct as stated.
* **L10.** All five closed forms re-derived and verified exactly (above).
* **Convention consistency.** The full-lattice, no-\(1/2\) convention is
  used consistently in L0–L11 and matches the certificates' \(K\); see M8
  for the JSON-side caveat.

---

**REFEREE: MAJOR-REVISIONS**
