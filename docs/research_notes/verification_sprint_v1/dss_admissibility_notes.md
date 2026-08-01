# DSS literature / admissibility gate — Workstream B

**Verification sprint v1, 2026-08-02. No new ideas: verify, correct, or kill only.**

Companion file: [`dss_admissibility_matrix.csv`](dss_admissibility_matrix.csv).

**Honesty labels used throughout.**

- `[EXACT]` — closed-form scaling algebra done here by hand (no floating point anywhere in
  this document's derivations).
- `[SRC-V]` — quoted from a primary source rendering I could read as text
  (arXiv abstract pages, ar5iv HTML). Treated as verbatim or near-verbatim.
- `[SRC-S]` — secondary/paraphrase; **not** used for any load-bearing cell.
- `[REJECTED]` — an extraction attempt that produced text contradicting the paper's own
  abstract; discarded, and recorded in §8 rather than deleted.
- `[REPO]` — a number or claim taken from the repository; its own label (exact vs float)
  is carried over.

Nothing in this document is a PDE theorem, a singularity, or a Clay statement. It is an
admissibility audit: which published exclusion theorems can and cannot be applied to five
recorded candidates.

---

## 1. The three similarity classes, stated precisely

Let \(T\) be the candidate blow-up time and let the NS scaling group be
\(u\mapsto \lambda u(\lambda x,\lambda^2 t)\).

### 1.1 Class 1 — stationary (continuous) backward self-similar (Leray)

\[
u(x,t)=\frac{1}{\sqrt{2a(T-t)}}\,U\!\left(\frac{x}{\sqrt{2a(T-t)}}\right),\qquad a>0,
\]
equivalently \(u(x,t)=\lambda u(\lambda x,\lambda^2t-(\lambda^2-1)T)\) for **every**
\(\lambda>0\). The profile solves the stationary Leray system; as displayed in
Chae–Wolf (arXiv:1609.06962) `[SRC-V]`:

> "\(-\Delta U+(U\cdot\nabla)U+ay\cdot\nabla U+aU=-\nabla P,\qquad \nabla\cdot U=0\) in \(\mathbb R^3\)"

### 1.2 Class 2 — backward \(\lambda\)-DSS

Invariance only under the discrete subgroup \(\lambda^{\mathbb Z}\), \(\lambda>1\) fixed.
Chae (arXiv:1306.0305) `[SRC-V]`, blow-up time normalised to \(0\):

> "there exists \(\lambda\neq1\) such that \(\lambda v(\lambda x,\lambda^2t)=v(x,t)\) for all \((x,t)\in\mathbb R^3\times(-\infty,0)\)"

Chae–Wolf (arXiv:1610.09464) `[SRC-V]` use the same relation and call it *backward
\(\lambda\)-DSS* for \(\lambda\in(1,+\infty)\):

> "\(u(x,t)=\lambda u(\lambda x,\lambda^2 t)\)"

**CORRECTION to the sprint brief.** The brief proposed
\(u(x,t)=\lambda u(\lambda x,\lambda^2t+(\lambda^2-1)T)\). The sign is wrong. Requiring the
map to fix \((0,T)\) and preserve \(\{t<T\}\) forces \(s-T=\lambda^2(t-T)\), i.e.

\[
\boxed{\;u(x,t)=\lambda\,u\!\bigl(\lambda x,\;\lambda^2t-(\lambda^2-1)T\bigr)\;}
\qquad\text{(backward }\lambda\text{-DSS about }(0,T)).
\]
At \(T=0\) this reduces to the sources' relation. `[EXACT]`
The repository's own form in `ideas_2026_08_01/audit_front_flow_seed.md` §T.3,
"\(u(x,t)=2u(2x,\,T-4(T-t))\)", is **correct**: \(T-4(T-t)=4t-3T=\lambda^2t-(\lambda^2-1)T\)
at \(\lambda=2\). `[EXACT]`

Equivalent similarity-variable form (used below): with \(y=x/\sqrt{2a(T-t)}\),
\(s=-\log(T-t)\), the profile \(V(y,s)\) is \(S_0\)-periodic in \(s\) with
\(S_0=2\log\lambda\). For \(\lambda=2\), \(S_0=\log4\). Note the repo's front-flow clock
gives the doubling period \(S=\log2/a_+\) in the \(ds=N^2dt\) clock — a different but
consistent normalisation (\(a_+\) rescales \(s\)).

### 1.3 Class 3 — asymptotically (D)SS

No exact invariance at any finite time; the *rescaled* solution converges as \(t\uparrow T\).
Chae–Wolf ARMA 2017 `[SRC-V]` (asymptotically SS, their Thm 1.5 hypothesis):

> "\(\lim_{t\nearrow t_*}(t_*-t)^{(q-3)/(2q)}\sup_{t<\tau<t_*}\bigl\|u(\cdot,\tau)-\tfrac{1}{\sqrt{2a(t_*-\tau)}}U(\cdots)\bigr\|_{L^q}=0\)"

Chae 2013 `[SRC-V]` (locally asymptotically **discretely** self-similar; the limit object is
a *time-periodic* profile \(\bar V(y,s)=\bar V(y,s+S_0)\), \(S_0\neq0\)):

> "\(\lim_{t\uparrow0}(-t)^{(q-3)/(2q)}\sup_{t<\tau<0}\bigl\|v(\cdot,\tau)-\tfrac{1}{\sqrt{-\tau}}\bar V(\cdot/\sqrt{-\tau},-\log(-\tau))\bigr\|_{L^q(B(0,R\sqrt{-t}))}=0\)"

---

## 2. Theorem statements with sources

### T1. Nečas–Růžička–Šverák 1996

*Acta Math.* **176** (1996) 283–294, DOI [10.1007/BF02551584](https://doi.org/10.1007/BF02551584).
Hypotheses (as recorded in `docs/known_obstructions.md` §1.1 `[REPO]`, primary-source
verified there, and as restated in Chae–Wolf's introduction `[SRC-V]`):
\(U\in W^{1,2}_{\rm loc}(\mathbb R^3)\), \(\nabla\cdot U=0\), solving the Leray profile
system weakly against compactly supported divergence-free test functions, **and
\(U\in L^3(\mathbb R^3)\)**. Conclusion: \(U\equiv0\).
Chae–Wolf restate it as `[SRC-V]`: "For Nečas–Růžička–Šverák (case \(p=3\)): If
\(U\in L^3(\mathbb R^3)\), then \(U=0\)."

### T2. Tsai 1998 (+ 1999 erratum)

*Arch. Ration. Mech. Anal.* **143** (1998) 29–51,
DOI [10.1007/s002050050099](https://doi.org/10.1007/s002050050099); erratum
*ARMA* **147** (1999) 363.
Two statements matter.

(T2a) **Local-energy version.** \(u\) is a weak solution on
\(Q_1(0,T)=B_1(0)\times(T-1,T)\) of the **exact** Leray backward self-similar form, and
\[
\operatorname*{ess\,sup}_{T-1<t<T}\int_{B_1(0)}|u|^2\,dx+\int_{T-1}^{T}\!\!\int_{B_1(0)}|\nabla u|^2\,dx\,dt<\infty .
\]
Conclusion \(u\equiv0\). Chae–Wolf's restatement `[SRC-V]`:
"\(\sup_{t\in(-t_0,0)}\int_B|u(t)|^2dx+\int_{-t_0}^0\int_B|\nabla u|^2dxdt<+\infty\)".
`docs/known_obstructions.md` §1.2 `[REPO]` adds the essential rider that the ball must
**contain the similarity centre**; finiteness on a ball away from the centre is not the
hypothesis.

(T2b) **Profile-integrability version.** \(U\in W^{1,2}_{\rm loc}\cap L^q(\mathbb R^3)\):
\(U=0\) for \(3<q<\infty\); \(U\) spatially constant for \(q=\infty\). Chae–Wolf `[SRC-V]`:
"if \(U\in L^p(\mathbb R^3)\) for some \(p\in[3,+\infty]\), then \(U=0\) for \(p\in[3,\infty)\)".

### T3. Escauriaza–Seregin–Šverák 2003

*Russian Math. Surveys* **58**:2 (2003) 211–250,
[mathnet.ru/eng/rm609](https://www.mathnet.ru/eng/rm609). Abstract `[SRC-V]`:

> "It is shown that the \(L_{3,\infty}\)-solutions of the Cauchy problem for the three-dimensional Navier–Stokes equations are smooth."

**Notation warning — load-bearing.** In the Seregin-school convention
\(\|f\|_{L_{p,q}(Q)}=\bigl(\int(\int|f|^p dx)^{q/p}dt\bigr)^{1/q}\), verified verbatim from
Seregin arXiv:1201.1100 `[SRC-V]`, so \(L_{3,\infty}\) is the **mixed space-time norm**
\(\operatorname{ess\,sup}_t\|u(\cdot,t)\|_{L^3_x}\), i.e. the endpoint Prodi–Serrin pair
\((p,q)=(3,\infty)\) — **not** the Lorentz space weak-\(L^3\). Two web summarisers asserted
the Lorentz reading; both are wrong (§8). Regularity under a bounded
\(L^\infty_tL^{3,\infty}_x\) (Lorentz) norm is **open** except under a smallness condition.
This distinction is what makes every row of the matrix possible: all five candidates have
**bounded Lorentz weak-\(L^3\)** and **divergent (or infinite) \(L^3\)**.

### T4. Chae–Wolf, DSS singularity removal

"Removing discretely self-similar singularities for the 3D Navier–Stokes equations",
arXiv:[1610.09464](https://arxiv.org/abs/1610.09464), *Comm. PDE* **42** (2017).
Abstract `[SRC-V]`:

> "We study the scenario of discretely self-similar blow-up for Navier–Stokes equations. We prove that at the possible blow-up time such solutions only one point singularity. In case of the scaling parameter \(\lambda\) near \(1\) we remove the singularity."

Theorem statements `[SRC-V]` (ar5iv rendering):

> **Theorem 1.1.** "For \(3\le p<+\infty\) let \(u\in C((-\infty,0);L^p(\mathbb R^3))\cap C^\infty(Q)\) be a solution to the Navier–Stokes equations, and \(\lambda\)-DSS for some \(\lambda\in(1,+\infty)\). Then the solution \(u\) is regular on \(\bar Q\setminus\{(0,0)\}\), and satisfies the estimate \(|u(x,t)|\le C/(\sqrt{-t}+|x|)\)."

> **Theorem 1.3.** "For every \(C_+>0\) there exists \(\lambda_+>1\) depending on \(C_+\) such that if \(u\in C^\infty(Q)\) is a \(\lambda\)-DSS solution [of] the Navier–Stokes equations for \(\lambda\in(1,\lambda_+)\), which satisfies \(|u(x,t)|\le C_+/(\sqrt{-t}+|x|)\ \forall(x,t)\in Q\). Then \(u\equiv0\)."

**\(\lambda\) range: \((1,\lambda_+)\) with \(\lambda_+=\lambda_+(C_+)\) not explicit.** Per the
sprint protocol this is recorded as **UNKNOWN** for \(\lambda=2\), unresolved constant
named \(\lambda_+(C_+)\) (§7, U1).

### T5. Chae–Wolf, Liouville theorems for self-similar profiles

"On the Liouville type theorems for self-similar solutions to the Navier–Stokes equations",
arXiv:[1609.06962](https://arxiv.org/abs/1609.06962), *ARMA* **225** (2017) 549–572,
DOI [10.1007/s00205-017-1110-7](https://doi.org/10.1007/s00205-017-1110-7).
Abstract `[SRC-V]`: "We prove Liouville type theorems for the self-similar solutions to the
Navier–Stokes equations. One of our results generalizes the previous ones by
Nečas–Růžička–Šverák and Tsai."

> **Corollary 1.4.** "Let \((U,P)\in C^\infty(\mathbb R^3)^3\times C^\infty(\mathbb R^3)\) be a solution to [the Leray system]. Suppose that for some \(3/2<p<+\infty\): \(U\in L^{p,\infty}(\mathbb R^3)\)." — then \(U\) is constant. `[SRC-V]`

> **Theorem 1.5.** asymptotically self-similar blow-up with profile in \(L^{p,\infty}(\mathbb R^3)\), \(3/2<p<+\infty\), is excluded (\(U=0\), the point is not a blow-up point), under the \(L^q\) convergence condition quoted in §1.3. `[SRC-V]`

Here \(L^{p,\infty}\) **is** the Lorentz space (this paper's own convention), so \(p=3\) is
exactly weak-\(L^3\). This is the theorem that kills candidate (a) most cheaply.

### T6. Chae 2013, asymptotically discretely self-similar

"Remarks on the asymptotically discretely self-similar solutions of the Navier–Stokes and
the Euler equations", arXiv:[1306.0305](https://arxiv.org/abs/1306.0305).
Abstract `[SRC-V]`:

> "…We prove that there exists no such locally asymptotically discretely self-similar blow-up for the 3D Navier–Stokes equations if the blow-up profile is a time periodic function belonging to \(C^1(\mathbb R;L^3(\mathbb R^3)\cap C^2(\mathbb R^3))\)."

> **Theorem 1.1.** "If \(v\in C(-\infty,0;L^3(\mathbb R^3))\) a solution to (NS), which blows up at \(t=0\), then \(t=0\) is not a time for discretely self-similar blow up." `[SRC-V]`

**This is the sharpest DSS exclusion found and it carries NO \(\lambda\) restriction.** The
price is the hypothesis \(v\in C((-\infty,0);L^3(\mathbb R^3))\), i.e. a **finite \(L^3\)
norm**. See §6, CORRECTION 3.

### T7. Seregin (Barker–Seregin), local Type-I and Liouville

arXiv:[1811.00502](https://arxiv.org/abs/1811.00502), *J. Math. Fluid Mech.* (2019).
Abstract `[SRC-V]`: "We prove that suitable weak solutions of the Navier–Stokes equations
exhibit Type I singularities if and only if there exists a non-trivial mild bounded ancient
solution satisfying a Type I decay condition."

> **Theorem 1.2.** "If \(v\) is a mild ancient solution satisfying \(\sup_{k\in\mathbb N}\|v(\cdot,t_k)\|_{L^3}<\infty\) for a sequence \(t_k\downarrow-\infty\), then \(v\equiv0\)." `[SRC-V]`

And, directly on the sprint's question `[SRC-V]`:

> "many questions concerning feasible Type I scenarios, e.g., discretely self-similar blow-up, remain completely open."

### T8. Bradshaw–Tsai — DSS **existence** (relevance: no exclusion weight)

- "Forward discretely self-similar solutions of the Navier–Stokes equations II",
  arXiv:[1510.07504](https://arxiv.org/abs/1510.07504), *Ann. Henri Poincaré* (2017).
- "Discretely self-similar solutions to the Navier–Stokes equations with data in
  \(L^2_{\rm loc}\) satisfying the local energy inequality",
  arXiv:[1801.08060](https://arxiv.org/abs/1801.08060), *Anal. PDE* **12** (2019) 1943.

These construct **forward** DSS (and DSS suitable weak, local-energy-inequality) solutions
for DSS data of **arbitrarily large weak-\(L^3\) norm** `[SRC-V]`. They therefore
(i) establish that the DSS class is non-empty and robust, (ii) show that largeness in
\(L^{3,\infty}\) alone is not an obstruction, and (iii) carry **zero** exclusion weight
against backward DSS blow-up. They are *supporting* rather than *threatening* citations,
and the repo should record them as such.

### T9/T10. CSTY2009, KNSS2009

Carried over from `docs/known_obstructions.md` §2.2, §2.3 `[REPO]` (primary-source
verified there). Both require **axial symmetry** on \(\mathbb R^3\); KNSS additionally
requires a **bounded ancient** solution (and no-swirl for Thm 5.2, or \(|u|\le C/r\) for
Thm 5.3). None of the five candidates is axisymmetric, so both are inapplicable across the
board. This is a *weakness* of the candidate portfolio's protection, not a strength: it
means the axisymmetric literature simply has no purchase, in either direction.

---

## 3. Per-candidate one-line derivations `[EXACT]`

Repo conventions: \(E_N=c_E/N\); critical shell density \(e_c=2c_E\xi^{-2}\); hence
\(|\hat U(\xi)|\asymp\sqrt{c_E}\,\xi^{-2}\) and, since
\(\widehat{|x|^{-1}}\propto|\xi|^{-2}\) in \(\mathbb R^3\),
\(U(x)\asymp A|x|^{-1}\) with \(A\asymp\sqrt{c_E}\).

**(a) steady front, critical wake.**
- \(\|U\|_3^3=A^3\!\int_{|x|>1}|x|^{-3}dx=4\pi A^3\!\int_1^\infty\!\frac{dr}{r}=+\infty\) (core smooth by the UV cutoff, so only the far field diverges) ⇒ **\(L^3\) infinite, logarithmically**.
- \(|\{|U|>\alpha\}|=\frac{4\pi}{3}(A/\alpha)^3\) ⇒ \(\|U\|_{L^{3,\infty}}=(4\pi/3)^{1/3}A\) ⇒ **weak-\(L^3\) finite**.
- With \(\lambda(t)=\sqrt{2a(T-t)}\): \(\int_{B_1}|u|^2=\lambda\int_{B_{1/\lambda}}|U|^2\asymp A^2\) and \(\int_{B_1}|\nabla u|^2=\lambda^{-1}\int_{B_{1/\lambda}}|\nabla U|^2\asymp A^2\lambda^{-1}\), and \(\int_{T-1}^T\lambda^{-1}dt=\int(2a(T-t))^{-1/2}dt<\infty\) ⇒ **local energy finite**.
- \(\sqrt{T-t}\,\|u(t)\|_\infty=\|U\|_\infty/\sqrt{2a}\) ⇒ **Type-I bounded**.

**(b) \(\lambda=2\) DSS front orbit, same wake.** Same envelope, log-periodically modulated;
all four quantities identical to (a). One extra exact fact used repeatedly below:
\(\|\lambda u(\lambda\cdot)\|_{L^3}=\lambda\cdot\lambda^{-1}\|u\|_{L^3}=\|u\|_{L^3}\), so
for an exact \(\lambda\)-DSS solution \(t\mapsto\|u(\cdot,t)\|_{L^3}\) is **exactly invariant
under \(T-t\mapsto\lambda^2(T-t)\)**.

**(c) multi-type cycle, super-period ratio 2.** After one super-period the map is dilation
by \(\prod_j\lambda_j=2\); hence the object *is* backward 2-DSS and every entry equals (b).

**(d) Zeno packet relay** (\(\gamma=1/2\), \(\lambda_j=2^j\), \(E_j=\lambda_j^{-1}\),
\(\Delta t_j=a_c\lambda_j^{-2}\), \(w_j=\lambda_jR_jW(\lambda_jR_j^T(x-x_j))\)).
- \(\|\sum_{j\le J}w_j\|_3^3=(J+1)\|W\|_3^3\) (disjoint supports) ⇒ finite for each \(t<T\), \(\asymp\log\frac{1}{T-t}\).
- Lacunary weak-\(L^3\): \(|\{|u|>\alpha\}|=\sum_j|\{|w_j|>\alpha\}|\) and only \(\lambda_j\gtrsim\alpha/\|W\|_\infty\) contribute, each with measure \(\propto\lambda_j^{-3}\); the geometric sum gives \(\alpha^3|\{|u|>\alpha\}|\le(1+2^{-3}+2^{-6}+\cdots)\sup_\alpha\alpha^3|\{|W|>\alpha\}|=\tfrac87\|W\|_{3,\infty}^3\) ⇒ **weak-\(L^3\) finite and \(J\)-independent**.
- \(\sup_t\int|u|^2=\sum_j\lambda_j^{-1}\|W\|_2^2\le2\|W\|_2^2\); \(\iint|\nabla u|^2=\sum_j\lambda_j\|\nabla W\|_2^2\,a_c\lambda_j^{-2}\le2a_c\|\nabla W\|_2^2\) ⇒ **local energy finite**.
- \(\|u\|_\infty\asymp\lambda_{J}\asymp N_c(t)\asymp(T-t)^{-1/2}\) ⇒ **Type-I bounded**.

**(e) torus phase-coded cloud + wake.** Cores of diameter \(N_j^{-1}\), amplitude \(N_j\):
same lacunary computation ⇒ weak-\(L^3\) finite; \(\sum_jc_E/N_j=2c_E/N_0\) and
\(\int_0^TN\,dt<\infty\) ⇒ local energy finite; \(\|u\|_\infty\asymp N\asymp(T-t)^{-1/2}\)
⇒ Type-I bounded; \(\|u\|_3^3\gtrsim J\asymp\log\frac1{T-t}\) (candidate (8.3), conditional
on unproved (8.1)–(8.2)). **Domain \(\mathbb T^3\) ⇒ no dilation group ⇒ no similarity
class at all.**

---

## 4. Second matrix: theorems not in the CSV column set

The CSV honours the requested 19 columns exactly. The three additional theorems that turned
out to matter are tabulated here.

| candidate | Chae2013 T6 kills? | first failing hypothesis | Seregin2019 T7 Thm 1.2 kills? | first failing hypothesis | Bradshaw–Tsai T8 |
|---|---|---|---|---|---|
| (a) steady | FALSE | not DSS (it is exactly SS); and \(v\notin C_t L^3\) | FALSE | not a *bounded* ancient solution; \(\|u(t)\|_3=\infty\) | n/a (forward existence only) |
| (b) \(\lambda=2\) DSS | **FALSE — by one hypothesis only** | \(v\in C((-\infty,0);L^3(\mathbb R^3))\): here \(\|u(t)\|_3\equiv+\infty\) | FALSE | \(\sup_k\|v(t_k)\|_{L^3}<\infty\) fails; solution unbounded | supports non-emptiness of DSS class |
| (c) multi-type | FALSE | identical to (b) | FALSE | identical to (b) | same |
| (d) Zeno | FALSE | ADSS version needs profile in \(C^1(\mathbb R;L^3\cap C^2)\); stack profile is in \(L^{3,\infty}\setminus L^3\). Also \(\mathbb T^3\) + forcing. | FALSE | \(\mathbb R^3\), bounded ancient, unforced | same |
| (e) cloud | FALSE | domain \(\mathbb T^3\) (no dilation group) | FALSE | \(\mathbb R^3\), bounded ancient | same |

**Reading.** Every single FALSE in the \(\lambda\)-unrestricted column (T6) is caused by the
*same* hypothesis: the profile must lie in \(L^3\). That is the entire surviving margin of
the DSS lane.

---

## 5. Per-cell justification pointers

- `NRS_failing_hypothesis` for (b)/(c): the *first* failure is not \(U\notin L^3\) but the
  absence of a time-independent profile — a DSS orbit does not satisfy the stationary Leray
  system at any time. Recording \(U\notin L^3\) as "the" failure would be a hypothesis-order
  error, because it would suggest the candidate could be killed by making the wake
  \(L^3\)-integrable, when in fact \(L^3\)-integrability triggers T6 instead.
- `Tsai_failing_hypothesis` for (b)–(e): Tsai's local-energy theorem is stated for \(u\) of
  the exact Leray form. Since (b)/(c) *do* satisfy the finite-local-energy condition
  (verified in §3), the exact-SS-form hypothesis is the only thing standing between them and
  \(u\equiv0\). This is the repo's T.4 gate restated as a hypothesis-order statement.
- `ChaeWolf_failing_hypothesis` for (b)/(c): Thm 1.1's hypothesis
  \(u\in C((-\infty,0);L^p(\mathbb R^3))\) for some \(p\in[3,\infty)\) **is satisfied**:
  \(|x|^{-1}\in L^p(\mathbb R^3)\) iff \(p>3\), and the core is smooth. Therefore Thm 1.1
  *applies* and the candidate inherits (i) a single point singularity and (ii) the pointwise
  Type-I bound \(|u|\le C_+/(\sqrt{-t}+|x|)\). The candidate already claims both, so this is
  a consistency check that **passes**, not a kill. The kill would come from Thm 1.3, whose
  \(\lambda\) window is unresolved.
- `ESS_note` everywhere: see the notation warning in §T3. This single distinction is what
  every candidate's escape route depends on, and it was previously implicit in the repo.
- `CSTY_*`/`KNSS_*`: uniformly FALSE on axial symmetry. Retained explicitly so that no
  future note claims "the candidate survives CSTY/KNSS" as evidence — it survives them
  vacuously.

---

## 6. Corrections and errata produced by this workstream

**CORRECTION 1 (sprint brief).** Backward DSS relation sign: correct form is
\(u(x,t)=\lambda u(\lambda x,\lambda^2t-(\lambda^2-1)T)\). The repo's §T.3 form is right;
the brief's is wrong. `[EXACT]`

**CORRECTION 2 (strengthens repo audit T.2).** The audit says NRS is "evaded by a hair"
because \(|x|^{-1}\in L^{3,\infty}\setminus L^3\), and that Tsai then kills the steady case.
Both are confirmed. But the log-escape from NRS buys **nothing**: the very norm it leaves
finite, weak-\(L^3\), is the hypothesis of **Chae–Wolf ARMA 2017 Corollary 1.4**
(\(U\in L^{p,\infty}\), \(3/2<p<\infty\), \((U,P)\in C^\infty\) ⇒ \(U\) constant ⇒ \(U\equiv0\)
with decay). So the steady front is killed **twice**, and the second kill needs no local
energy hypothesis at all. Candidate (a) status: **KILLED, confirmed and strengthened.**

**CORRECTION 3 (new, structural — the most important item here).** For an exact backward
\(\lambda\)-DSS solution on \(\mathbb R^3\), \(\|u(\cdot,t)\|_{L^3}\) is invariant under
\(T-t\mapsto\lambda^2(T-t)\) `[EXACT]`. Hence it is **log-periodic**, so it is either
identically \(+\infty\) or uniformly bounded on \((-\infty,T)\) — **there is no
log-growing branch**. Consequences:

1. If it is bounded (and \(v\in C_tL^3\)), **Chae 2013 Theorem 1.1 kills it for every
   \(\lambda>1\)** — no \(\lambda\) window, no smallness. A finite-\(L^3\) DSS candidate is
   dead on arrival.
2. Therefore any \(\mathbb R^3\) \(\lambda\)-DSS candidate in this repo **must** have
   \(\|u(t)\|_{L^3}\equiv+\infty\). This is now a hard, pre-registered admissibility
   requirement, not a stylistic choice.
3. The repository's advertised signature \(\|u(t)\|_3^3\asymp\log\frac1{T-t}\to\infty\) is
   therefore **not** a property of an exactly-DSS \(\mathbb R^3\) object. It is a
   **torus/lattice artifact**: the IR cutoff \(|k|\ge N_0\) truncates the wake to \(J(t)\)
   octaves, which is exactly what breaks the dilation group. Wording in
   `CANDIDATE_SOLUTION_PHASE_CODED_LERAY_CLOUD.md` §8–§9 and in the front-flow portfolio §2
   should stop describing the \(\mathbb R^3\) limit object and the torus object with the
   same \(L^3\) law: they differ, and the difference is precisely the DSS symmetry.

**CORRECTION 4 (kills a claimed escape route).** "Multi-type log-periodic cycle with
super-period ratio 2" is **not** a distinct similarity class: after one super-period it is
backward 2-DSS, so its admissibility row is byte-identical to (b). Furthermore the DSS
literature supplies a *new* design rule: since Chae–Wolf CPDE 2017 Thm 1.3 kills
\(\lambda\in(1,\lambda_+)\), **lowering** the composite ratio toward 1 moves the candidate
*into* the known kill zone. Combined with repo audit S6.3 (equal-ratio \(L\)-cycles need the
lattice-incompatible \(\lambda=2^{1/L}\)), the conclusion is: never subdivide the composite
ratio; keep \(\prod_j\lambda_j\ge2\).

**CORRECTION 5 (notation, load-bearing).** ESS 2003's \(L_{3,\infty}\) is the mixed norm
\(L^\infty_tL^3_x\), not the Lorentz weak-\(L^3\) (§T3, verified from Seregin's own
definition of \(L_{p,q}\)). All five candidates are **uniformly bounded in Lorentz
weak-\(L^3\)** `[EXACT]`, which is legal: regularity under bounded \(L^\infty_tL^{3,\infty}_x\)
is open except under smallness. Any repo text that pairs "weak-\(L^3\)" with "ESS" must be
re-read against this. (`CANDIDATE_SOLUTION_PHASE_CODED_LERAY_CLOUD.md` §11 already lists
"Type-I, weak-\(L^3\), CKN — still open", which is correct.)

**CORRECTION 6 (upgrades the T.4 oscillation gate's justification).** The repo's gate —
*"an orbit whose per-period oscillation of the critical norm vanishes under refinement is
dead by Tsai"* — is **correct but under-cited**. The sharp statement is: a vanishing
oscillation means the rescaled solution converges to a time-independent profile, i.e. the
blow-up is *asymptotically self-similar*, and **Chae–Wolf ARMA 2017 Theorem 1.5** excludes
that for profiles in \(L^{p,\infty}\), \(3/2<p<\infty\) — which the critical wake satisfies
at \(p=3\). So the gate fires at the weak-\(L^3\) level and does not need the local-energy
route through Tsai. Recommendation: cite ARMA2017 Thm 1.5 as the primary authority for the
T.4 gate and Tsai1998 as the secondary.

**CONFIRMED (no change).** Repo audit item G8 ("DSS class is not safe; Chae–Wolf /
Bradshaw–Tsai / Seregin restrictions untested") was the right worry, and this workstream
discharges it: the class is **live but not safe**. Seregin (arXiv:1811.00502) states
verbatim that "discretely self-similar blow-up … remain[s] completely open" `[SRC-V]`, and
the only \(\lambda\)-unrestricted exclusion (Chae 2013) is evaded by exactly one hypothesis.

---

## 7. UNKNOWN cells — explicit list

| ID | UNKNOWN cell | unresolved object | what must be looked up or proved |
|---|---|---|---|
| **U1** | `ChaeWolf_kills` for (b) and (c) | **\(\lambda_+(C_+)\)** in Chae–Wolf CPDE 2017 Thm 1.3 | Read the proof of Thm 1.3 in arXiv:1610.09464 and extract the explicit dependence of \(\lambda_+\) on \(C_+\). The proof is a perturbation off \(\lambda=1\), so almost certainly \(\lambda_+\to1\) as \(C_+\to\infty\), and the repo's \(C_+\asymp\sqrt{c_E}\) with \(c_E\) large would put \(\lambda=2\) far outside the kill zone — **but this is a guess and must not be recorded as a result until the constant is read off.** Blocked here only by lack of a text-extractable copy of the proof (§8). |
| **U2** | the single margin in §4 | Can Chae 2013 Thm 1.1/1.2 (and Chae–Wolf ARMA Cor 1.4's method) be pushed from an \(L^3\) profile to an \(L^{3,\infty}\) profile? | This is the *highest-value* item in the whole workstream. If yes, candidates (a),(b),(c) all die and (d),(e) come under direct pressure. Search for post-2017 work extending DSS/ADSS Liouville theorems to Lorentz profiles; if none exists, attempt the extension (the obstruction is presumably the failure of \(L^3\) small-tail/absolute-continuity arguments in \(L^{3,\infty}\)). |
| **U3** | `ESS_note` for (d) and (e) | validity of the ESS/Seregin endpoint criterion on \(\mathbb T^3\) | `docs/known_obstructions.md` §10 already lists periodic-domain transfer as unconfirmed. Needed because the ESS cell is the *only* non-vacuous literature constraint on the two torus candidates. Locate a periodic-domain statement (Seregin's bounded-domain necessary-condition papers are the likely source). |
| **U4** | implicit in `weak_L3_finite` | the smallness constant \(\varepsilon_0\) in "\(\|u\|_{L^\infty_tL^{3,\infty}_x}\le\varepsilon_0\Rightarrow\) regular" | The candidates have \(\|u\|_{3,\infty}\asymp\sqrt{c_E}\) with \(c_E\) large `[REPO, and note audit S3.5: the specific values 228/546/1902 are invalid perturbative extrapolations]`. Almost certainly no collision, but the constant should be named and the inequality checked rather than assumed. |
| **U5** | all rows for (d),(e) | forced NS | Every theorem T1–T7 is stated for **unforced** NS. Candidates (d) and (e) are Clay (D) with a low-band force. No DSS/SS exclusion theorem for forced NS was located. This means the torus candidates are unreachable by this literature *for a second, independent reason* — and equally that no protection can be claimed. If a forced-NS version exists, it must be found. |
| **U6** | (e) `class` | is the lattice decimation \(k\mapsto2k\) ever a genuine symmetry on \(\mathbb T^3\)? | Answered NO here `[EXACT]` (the domain does not rescale), consistent with repo audit S7.1's rejection of subsampling. Recorded as resolved, but flagged because several repo notes still speak of "\(\lambda=2\) DSS on the torus", which is a category error: the torus object is at best *approximately* DSS in a moving window. |

---

## 8. Retained failures and rejected extractions (never deleted)

1. `[REJECTED]` A WebFetch summarisation of the raw PDF of arXiv:1610.09464 returned an
   "abstract" ("We study the regularity and decay properties of discretely self-similar
   solutions…") that **contradicts** the arXiv abstract page, together with an invented
   "\(\lambda_0\) is an explicit universal constant" and fabricated theorem hypotheses. It
   was discarded. Only the ar5iv HTML rendering was used. Lesson recorded: never accept a
   summariser's rendering of a compressed PDF for a hypothesis-level claim.
2. `[REJECTED]` A WebFetch summarisation of arXiv:1811.00502 defined Type I as
   \(\sup_{t<T}(T-t)\|u\|_{L^\infty}<\infty\) — wrong power; the correct normalisation is
   \(\sqrt{T-t}\,\|u\|_{L^\infty}\). Discarded; the ar5iv rendering was used instead.
3. `[REJECTED]` Two independent summarisers (a search-result summary and the mathnet fetch)
   asserted that ESS 2003's \(L_{3,\infty}\) is the **Lorentz** space. This is wrong and is
   the single most consequential error this audit had to defend against. Rebutted from
   Seregin's own definition of \(L_{p,q}\) (arXiv:1201.1100) `[SRC-V]`.
4. `[FAILED CHECK, retained]` PDF page rendering is unavailable in this environment
   (`pdftoppm` not installed), so the proof of Chae–Wolf CPDE 2017 Thm 1.3 could not be read
   and U1 could not be closed. The two PDFs are cached under the session `tool-results`
   directory. Installing poppler-utils, or obtaining the ar5iv rendering of the proof
   section, would close U1.
5. `[NEGATIVE RESULT, retained]` No theorem excluding **backward** DSS blow-up at
   \(\lambda=2\) under any of the candidates' actual hypotheses was found, after searching
   for Chae–Wolf, Bradshaw–Tsai, Bradshaw–Phelps (spatial decay of DSS, arXiv:2202.08352),
   Seregin, and 2019–2024 follow-ups. Independent confirmation from a 2018-era survey
   sentence `[SRC-V]`: "the existence of non-zero backward discretely self-similar solutions
   remains open."

---

## 9. Net effect on the portfolio

- Candidate (a) **KILLED**, twice, confirming the repo's own pre-registered gate.
- Candidates (b) and (c) are the **same object** and survive on exactly one hypothesis
  (\(L^3\) vs \(L^{3,\infty}\) profile), with one UNKNOWN (\(\lambda_+\)).
- Candidates (d) and (e) are **inadmissible to the entire SS/DSS literature** — for (e) the
  very first hypothesis (domain \(\mathbb R^3\) with a dilation group) fails. That is not
  protection; it means all exclusion pressure on them comes from ESS/Serrin/CKN and from
  repo-internal gates, and those gates must therefore be held to a higher standard.
- The repo's T.4 oscillation gate should be **retained and re-cited** (CORRECTION 6): it is
  the operational form of Chae–Wolf ARMA 2017 Thm 1.5, and it is now known to fire at the
  weak-\(L^3\) level rather than needing the local-energy route.
