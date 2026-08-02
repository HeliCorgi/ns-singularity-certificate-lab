# Hypothesis L\*: direct route for the sharp family, and a strict weakening

> **Status update (2026-08-02).** The paper's static no-go is now
> **unconditional**, via a *proven* capacity bound
> \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) with the sharp
> exponent, for the **smoothly** truncated family
> \(\widehat u_N(k)=\chi(|k|/N)P_kv_0/|k|^2\)
> ([`lstar_proof_main.md`](lstar_proof_main.md), Theorem 7.1(3)).
> Consequences for this note:
>
> * **§A (JOB A) stands, and is now more valuable.** Theorem A is a
>   **proven negative result**: for the sharp family, every
>   constant-vector sweeping split — any constant \(b\), any constant
>   pairing direction \(c\), any Hölder exponent pair — is capped at
>   \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\gtrsim N\), two powers short of
>   the target. It is recorded in the paper's audit §3a. Its scope note is
>   exactly right: it does *not* close pairing against a general
>   divergence-free test field, and the proof that succeeded takes
>   precisely that route, with the test field living at the concentration
>   scale, \(\psi_N(x)=N^{3/2}\Psi(Nx)\).
> * **§B (JOB B) is PROVEN but now MOOT.** The weakening
>   \(N^3\rightsquigarrow N^2g(N)\), \(\int^\infty dt/(t\,g(t))<\infty\),
>   is correct end to end (referee-verified), but it is not needed: the
>   full \(N^3\) is proven. Hypotheses (L\*-weak) and (L\*-min) are
>   therefore **not used** by anything.
> * **§B.3 is WITHDRAWN as an integration plan.** It prescribed keeping
>   the sharp family and replacing (L\*) by (L\*-weak). The plan actually
>   adopted switches the family instead, which makes the Proposition
>   unconditional and (L\*-weak) moot; the two plans are mutually
>   exclusive. §B.3 is retained below only as a record of the road not
>   taken — **do not apply it.**
> * **One by-product of §B survives and is used**: \(T_N\ge T_1=6\)
>   (Step 2 of §B.2), which makes the dyadic-gap bound elementary and
>   \(N\)-uniform, replacing the earlier asymptotic
>   "\(T_{2N}/T_N\to1\)". It is now part of Lemma 11 of
>   [`../complete_proof.md`](../complete_proof.md).
> * **(L\*) for the sharp family remains OPEN** — and is used by nothing.

Companion to [`../theorem_statement.md`](../theorem_statement.md) and
[`../complete_proof.md`](../complete_proof.md) (Lemmas 9–11 and Theorem O).
Two independent jobs:

- **JOB A** — attack the weak form of L\* *directly* for the sharply
  truncated family, using only the sweeping identity, the exact pairing
  bound, and exact lattice evaluation. **Result: the route is closed by a
  proven no-go.** See §A.5–A.6.
- **JOB B** — write and prove the weakest hypothesis that still drives the
  Proposition. **Result: `N^3` can be replaced by `N^2 g(N)` for any
  nondecreasing `g` with `∫^∞ dt/(t g(t)) < ∞`** — e.g.
  `N^2 (log N)^{1+ε}`. See §B.

Every claim below is labelled **PROVEN** / **NUMERICAL** / **OPEN**.
Numerics: binary64 unless labelled *exact* (`fractions.Fraction`); no
interval arithmetic anywhere in this note. Scripts are described in §D.

Throughout, \(v_0\in\mathbb Z^3\setminus\{0\}\),
\(B_N=\{k\in\mathbb Z^3:1\le|k|\le N\}\),
\[
\hat u_N(k)=\frac{P_kv_0}{|k|^2}\ (k\in B_N),\qquad
u_N(x)=\sum_{k\in B_N}\hat u_N(k)e^{ik\cdot x},
\]
full-lattice sums, \(\langle f,g\rangle=(2\pi)^{-3}\int_{\mathbb T^3}f\cdot g
=\sum_k\overline{\hat f_k}\cdot\hat g_k\),
\(H_r=\sum_k|k|^{2r}|\hat u_k|^2\),
\(S_N=\sum_{B_N}|k|^{-2}\), \(T_N=\sum_{B_N}|k|^{-4}\).
Numerical work uses \(v_0=(1,2,3)\), \(|v_0|^2=14\), \(|v_0|^4=196\).

---

## §A. JOB A — the direct route for the sharp family

### A.0 The two proven tools, restated and re-verified

**(T1) Sweeping identity.** For a constant \(c\in\mathbb R^3\) and
divergence-free zero-mean \(u\), \(\mathbb P((c\cdot\nabla)u)=(c\cdot\nabla)u\).
*Proof.* The coefficient of \((c\cdot\nabla)u\) at \(k\) is
\(i(c\cdot k)\hat u_k\), which is \(\perp k\) because \(\hat u_k\perp k\);
\(\mathbb P\) is the fibrewise projection onto \(k^\perp\). ∎ **PROVEN**
(this is the paper's own certificate remark; re-derived, correct).

**(T2) Pairing bound.** For any constant \(c\neq0\),
\[
\|\mathbb P(u\cdot\nabla u)\|_2\;\ge\;
\frac{|Q_c|}{\|(c\cdot\nabla)u\|_2},\qquad
Q_c:=\langle(u\cdot\nabla)u,\,(c\cdot\nabla)u\rangle .
\]
*Proof.* \(\langle\mathbb P(u\cdot\nabla u),(c\cdot\nabla)u\rangle
=\langle u\cdot\nabla u,\mathbb P((c\cdot\nabla)u)\rangle=Q_c\) by
self-adjointness of \(\mathbb P\) and (T1); then Cauchy–Schwarz. ∎ **PROVEN.**

**Sharpness (NUMERICAL).** With \(c=v_0\) the bound captures
\(\|\,\cdot\,\|\)-fraction \(91.2\%\to91.19\%\) over \(N=4,\dots,28\)
(binary64; reproduces the paper's stated 91.18–92.10 % range). So (T2) is
*not* the lossy step: a proof of \(Q_{v_0}\gtrsim N^2\) would give the
**full** L\*.

### A.1 Exact evaluation of \(\|(v_0\cdot\nabla)u_N\|_2^2\)

The coefficient of \((v_0\cdot\nabla)u_N\) at \(k\) is
\(i(v_0\cdot k)\hat u_N(k)\), so with
\(\cos\theta_k=(k\cdot v_0)/(|k||v_0|)\) and
\(|P_kv_0|^2=|v_0|^2\sin^2\theta_k\),
\[
\boxed{\;\|(v_0\cdot\nabla)u_N\|_2^2
=\sum_{k\in B_N}(v_0\cdot k)^2|\hat u_N(k)|^2
=|v_0|^4\sum_{k\in B_N}\frac{\cos^2\theta_k-\cos^4\theta_k}{|k|^2}\;}
\tag{A.1}
\]
**PROVEN** — this is exactly the formula announced in the task; it is a
one-line consequence of \(|P_kv_0|^2=|v_0|^2-(k\cdot v_0)^2/|k|^2\), no
symmetry lemma needed yet.

**Second moment (PROVEN, exact).** By Lemma 9 with \(f=|k|^{-4}\),
\(\sum_{B_N}k_ik_j|k|^{-4}=\delta_{ij}S_N/3\), hence
\(\sum\cos^2\theta_k/|k|^2=S_N/3\) **exactly, for every \(N\)**.

**Fourth moment — does the cubic anisotropy matter? (PROVEN, exact.)**
The hyperoctahedral group \(B_3\) (sign flips + permutations) leaves \(B_N\)
and \(|k|\) invariant but — unlike \(SO(3)\) — does **not** force the rank-4
moment tensor to be isotropic. Writing
\(M^{(4)}_{ijlm}=\sum_{B_N}k_ik_jk_lk_m|k|^{-6}\), \(B_3\)-invariance gives
\[
M^{(4)}_{ijlm}
= B_N^{\,*}\bigl(\delta_{ij}\delta_{lm}+\delta_{il}\delta_{jm}
+\delta_{im}\delta_{jl}\bigr)+(A_N-3B_N^{\,*})\,\Delta_{ijlm},
\]
\(A_N=\sum k_1^4|k|^{-6}\), \(B_N^{\,*}=\sum k_1^2k_2^2|k|^{-6}\),
\(\Delta_{ijlm}=1\) iff \(i=j=l=m\). Contracting \(i=j,\ l=m\) gives the
exact constraint
\[
3A_N+6B_N^{\,*}=S_N .
\tag{A.2}
\]
Define the **anisotropy defect** and the **quartic shape factor**
\[
\Delta_N:=\tfrac13 S_N-5B_N^{\,*}\quad(\;=0\iff A_N=3B_N^{\,*}\iff\text{isotropic}),
\qquad
\rho:=\frac{\sum_i v_{0,i}^4}{|v_0|^4}\in[\tfrac13,1].
\]
Then \(\sum\cos^4\theta_k/|k|^2=3B_N^{\,*}+(\tfrac13S_N-5B_N^{\,*})\rho\), and
substituting into (A.1):
\[
\boxed{\;\|(v_0\cdot\nabla)u_N\|_2^2
=|v_0|^4\Bigl[\tfrac{2}{15}\,S_N+\bigl(\tfrac35-\rho\bigr)\Delta_N\Bigr]\;}
\tag{A.3}
\]
**PROVEN** (pure \(B_3\) representation theory + Lemma 9 + (A.2)), and
**verified exactly** in `Fraction` arithmetic for \(N=4,6,8,12,16,24,32\)
(all assertions pass, `lstar_exact.py`).

Two structural readings of (A.3):

1. The whole anisotropy enters through the *single scalar* \(\Delta_N\),
   multiplied by \((\tfrac35-\rho)\). When the lattice is isotropic at
   fourth order the \(v_0\)-shape dependence cancels **identically** — the
   direction of \(v_0\) is invisible.
2. **NUMERICAL:** \(\Delta_N\) is **bounded**, not growing:
   \(\Delta_N=1.311,\,1.184,\,1.309,\,1.317,\,1.351,\,1.354,\,1.355\) at
   \(N=4,6,8,12,16,24,32\), while \(S_N\approx4\pi N\to\infty\). Since
   \(\Delta_N=\sum_{B_N}h(k/|k|)|k|^{-2}\) with
   \(h(\omega)=\tfrac13-5\omega_1^2\omega_2^2\) of **zero mean on \(S^2\)**,
   boundedness is exactly what standard lattice-point equidistribution in
   spherical caps (error \(O(r^2)\) per ball) predicts. *A rigorous proof of
   \(\Delta_N=O(1)\) is standard but is **not** carried out here — it is
   **not needed** below.* So: **the cubic fourth-moment anisotropy does not
   matter at leading order**; it is an \(O(1)\) additive correction to a
   quantity of size \(\Theta(N)\).

**Elementary two-sided bounds (PROVEN, no equidistribution input).** From
(A.1), \(\cos^2\theta\,\sin^2\theta\le\tfrac14\), so
\[
0\;\le\;\|(v_0\cdot\nabla)u_N\|_2^2\;\le\;\tfrac14|v_0|^4S_N .
\tag{A.4}
\]
A matching lower bound \(\ge c_-|v_0|^4S_N\) follows by restricting to
\(\{N/2<|k|\le N,\ \cos^2\theta_k\in[\tfrac14,\tfrac34]\}\) — a shell-cone
of positive solid angle \((\sqrt3-1)/2\) of the sphere, containing
\(\ge cN^3\) lattice points by the standard piecewise-smooth-boundary count
— each point contributing \(\ge(3/16)N^{-2}\). Constants left crude; only
(A.4) is used in §A.5. **PROVEN** (modulo the routine lattice count).

**Asymptotic value (NUMERICAL, and PROVEN given \(\Delta_N=O(1)\)):**
\(\|(v_0\cdot\nabla)u_N\|_2^2=\tfrac{2}{15}|v_0|^4S_N\,(1+O(1/N))\).
Measured \(\|(v_0\cdot\nabla)u_N\|^2/(|v_0|^4S_N)\): \(0.136554,\,0.135099,\,
0.134771,\,0.134265,\,0.134037,\,0.133796,\,0.133678\) at
\(N=4,6,8,12,16,24,32\) versus \(2/15=0.133\overline3\) — and (A.3) reproduces
each value **to the last printed digit** with \(\rho=1/2\) (for
\(v_0=(1,2,3)\): \(\sum v_i^4=98\), \(|v_0|^4=196\)).

*(Consequence for the paper: \(\|(v_0\cdot\nabla)u_N\|_2^2\asymp S_N\asymp N\),
so (T2) converts a bound \(Q_{v_0}\ge cN^{1+a/2+1/2}\) into
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c'N^{2+a}\). The full L\* needs
exactly \(Q_{v_0}\gtrsim N^2\).)*

### A.2 The sweeping split and its first piece

The moment matrix
\[
M_{ij}:=\sum_{k\in B_N}k_ik_j|\hat u_N(k)|^2\qquad(\text{real, symmetric, PSD}),
\qquad \operatorname{tr}M=H_1 ,
\]
satisfies \(\langle(b\cdot\nabla)u,(c\cdot\nabla)u\rangle=b^{\!\top}\!Mc\) for
constants \(b,c\). Splitting \(u_N=b+w\), \(w=u_N-b\), and using
\((u\cdot\nabla)u=(b\cdot\nabla)u+(w\cdot\nabla)u\):
\[
Q_c=b^{\!\top}\!Mc+R_{b,c},\qquad
R_{b,c}=\langle(w\cdot\nabla)u_N,(c\cdot\nabla)u_N\rangle .
\tag{A.5}
\]
With the paper's natural choice \(b=u_N(0)=\tfrac23S_Nv_0\) (Lemma 10) and
\(c=v_0\), the first piece is **exactly**
\[
b^{\!\top}\!Mv_0=\tfrac23S_N\,\|(v_0\cdot\nabla)u_N\|_2^2
\;\overset{(A.3)}{=}\;\tfrac{4}{45}|v_0|^4S_N^2
+\tfrac23S_N|v_0|^4(\tfrac35-\rho)\Delta_N\;>\;0 .
\tag{A.6}
\]
**PROVEN.** Leading constant \(\tfrac{4}{45}|v_0|^4(4\pi)^2\approx2751\)
in units of \(N^2\) — this is the "\(\sim\!2750\)" figure.

**Exact spectrum of \(M\) (PROVEN).** Expanding as in §A.1,
\[
M_{ij}=\delta_{ij}\bigl[|v_0|^2(\tfrac13S_N-B_N^{\,*})-\Delta_N v_{0,i}^2\bigr]
-2B_N^{\,*}v_{0,i}v_{0,j},
\qquad \operatorname{tr}M=\tfrac23|v_0|^2S_N=H_1\ \checkmark
\]
so in the isotropic limit \(M\) has eigenvalues
\(\tfrac15H_1\) along \(v_0\) and \(\tfrac25H_1\) twice, transverse.
Measured \(\lambda(M)/H_1\) at \(N=32\): \((0.20051,\,0.39896,\,0.40053)\)
(**NUMERICAL**, matches). And **elementarily, for every unit \(e\)**,
\[
e^{\!\top}\!Me=|v_0|^2\underbrace{\sum(e\cdot k)^2|k|^{-4}}_{=\,S_N/3\ \text{(Lemma 9)}}
-\underbrace{\sum(e\cdot k)^2(k\cdot v_0)^2|k|^{-6}}_{\ge\,0}
\;\le\;\tfrac13|v_0|^2S_N=\tfrac12H_1 ,
\tag{A.7}
\]
i.e. \(\lambda_{\max}(M)\le\tfrac12H_1\). **PROVEN, elementary.**

### A.3 Numerical status of the split (NUMERICAL, binary64)

Independently recomputed with a clean full-lattice spectral pipeline
(`lstar_direct.py`; \(n=4N+2\) grid, alias-free for all product modes).
Cross-checks passed: \(K(u_4)=0.7884\), \(K(u_8)=2.0372\) reproduce the
paper's exact certificates; the pointwise \(Q\)-density integrates to \(Q\)
to \(4\times10^{-16}\) relative.

| \(N\) | \(S_N\) | \(Q_{v_0}/N^2\) | first\(/N^2\) | remainder\(/N^2\) | \(|R|/\text{first}\) | \(\|\mathbb P(u\!\cdot\!\nabla u)\|^2/N^3\) | \(K\) | pairing frac. |
|---|---|---|---|---|---|---|---|---|
| 4  | 40.70  | 640.8  | 1847.1 | −1206.4 | 0.6531 | 1777.4 | 0.7884 | 0.9210 |
| 8  | 91.06  | 850.7  | 2281.7 | −1431.0 | 0.6272 | 2874.3 | 2.0372 | 0.9151 |
| 12 | 141.29 | 926.7  | 2432.3 | −1505.6 | 0.6190 | 3321.4 | 3.3002 | 0.9134 |
| 16 | 191.84 | 969.2  | 2517.8 | −1548.6 | 0.6151 | 3580.0 | 4.5741 | 0.9127 |
| 20 | 242.14 | 993.1  | 2564.5 | −1571.4 | 0.6128 | 3729.9 | 5.8422 | 0.9122 |
| 24 | 292.46 | 1009.4 | 2596.0 | −1586.7 | 0.6112 | 3832.8 | 7.1113 | 0.9120 |
| 28 | 342.96 | 1022.2 | 2621.6 | −1599.4 | 0.6101 | 3913.9 | 8.3852 | 0.9119 |

**The orchestrator's warning is confirmed and sharpened.** The remainder is
negative and of the *same order* as the first piece, with the ratio
**converging** to \(\approx0.607\) — not decaying. So no improvement of
constants can make a triangle-inequality split work; the failure is
asymptotic, not a small-\(N\) artifact. (\(Q_{v_0}/N^2\) is itself increasing
towards \(\approx1090\), consistent with the full L\*; Richardson
extrapolation in \(1/N\) of \(Q/(|v_0|^4S_N^2)\) gives \(0.0352\).)

Best-case audits of the remainder estimate, as multiples of the first piece
(anything \(\ge1\) makes the bound vacuous):

| \(N\) | \(L^\infty\!\times\!L^2\) with \(b=u_N(0)\) | with the *measured* Chebyshev-optimal constant | \(L^6\!\times\!L^3\) | \(L^{2.5}\) | \(L^4\) | \(L^{20}\) |
|---|---|---|---|---|---|---|
| 8  | 2.302 | 1.185 | 4.551 | 16.72 | 7.00 | 2.704 |
| 16 | 2.266 | 1.156 | 6.308 | 37.89 | 11.56 | 2.953 |
| 24 | 2.256 | 1.148 | 7.675 | 61.32 | 15.57 | 3.122 |
| 32 | 2.250 | 1.143 | 8.833 | 86.38 | 19.26 | 3.251 |

Every Hölder exponent pair \((p,q)\), \(1/p+1/q=1/2\), gives the *same power*
\(N^{3/2}\) (the estimate is scale-invariant under the concentration
\(u\sim NV(Nx)\)), so it is purely a fight over constants — and
\(p=\infty\) is the best of them, converging to \(2.250\to\sqrt5\); all
others are worse and *diverge* with \(N\).

### A.4 Why: the integrand lives exactly where the split is wrong

**NUMERICAL.** The fraction of \(Q_{v_0}\) coming from \(|x|\le r_0\):

| \(N\) | \(r_0=4/N\) | \(8/N\) | \(16/N\) |
|---|---|---|---|
| 8  | 0.496 | 0.948 | 0.997 |
| 16 | 0.461 | 0.923 | 0.987 |
| 24 | 0.469 | 0.918 | 0.982 |

\(\approx92\%\) of \(Q\) comes from a ball of radius \(8/N\), i.e. from a set
of relative volume \(\approx8/N^3\). Also **NUMERICAL, to 12 digits for all
\(N\) tested:** \(\|u_N\|_\infty=|u_N(0)|=\tfrac23S_N|v_0|\) exactly, so the
paper's Bernstein sandwich \(\tfrac23S_N\le\|u_N\|_\infty^2/H_1\le S_N\) is
saturated at its **lower** end. *(Not a triangle-inequality artifact:
\(\sum_k|\hat u_N(k)|>|\sum_k\hat u_N(k)|\) strictly, since the \(\hat u_N(k)\)
are not parallel. Free small improvement to Lemma 10, if wanted; **OPEN** as
a proof.)*

Diagnosis: the split \(u=b+w\) is a **Galilean/sweeping** approximation,
valid only when \(u\) is nearly constant on the support of the integrand.
Here \(u_N\) varies by its **full magnitude** \(\Theta(N)\) across the ball
\(|x|\lesssim1/N\) that carries all of \(Q\). The ansatz is therefore
structurally, not quantitatively, wrong.

### A.5 No-go theorem for the whole sweeping-split family

The next result makes §A.3–A.4 into a proof, and it closes the route for
**every** constant \(b\) and **every** constant pairing direction \(c\) at once.

> **Theorem A (no-go; PROVEN, elementary).**
> Fix \(v_0\) and let \(C_0:=|v_0|\sqrt{\tfrac23T_\infty}\) (an
> \(N\)-independent constant). For every pair of constant vectors
> \(b\in\mathbb R^3\), \(c\in\mathbb R^3\setminus\{0\}\), the lower bound
> obtained from (T2) + the split (A.5) + the estimate
> \(|R_{b,c}|\le\|u_N-b\|_\infty\|\nabla u_N\|_2\|(c\cdot\nabla)u_N\|_2\)
> never exceeds
> \[
> \|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\;\ge\;C_1\,S_N,
> \qquad C_1=(3+2\sqrt2)\tfrac49|v_0|^4T_\infty .
> \]
> Since \(S_N\asymp N\), this scheme can **never** deliver
> \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^{2+a}\) for any \(a>0\) —
> indeed it cannot even reach \(N^2\).

*Proof.* Write \(H_1=\tfrac23|v_0|^2S_N\) (Lemma 10). From (A.5) and (T2),
the best bound this scheme yields is
\[
\|\mathbb P(u_N\cdot\nabla u_N)\|_2\;\ge\;
\frac{b^{\!\top}\!Mc}{\sqrt{c^{\!\top}\!Mc}}
-\|u_N-b\|_\infty\sqrt{H_1}
\;=:\;\mathcal B(b,c).
\]
(Only the first term depends on \(c\).) By Cauchy–Schwarz in the
semi-inner product \(\langle a,a'\rangle_M=a^{\!\top}\!Ma\),
\(\sup_{c\neq0}b^{\!\top}\!Mc/\sqrt{c^{\!\top}\!Mc}=\sqrt{b^{\!\top}\!Mb}\),
attained at \(c=b\). Hence
\(\sup_c\mathcal B(b,c)=\sqrt{b^{\!\top}\!Mb}-\|u_N-b\|_\infty\sqrt{H_1}\).
*(Note this already shows: since \(u_N(0)\parallel v_0\), the paper's choice
\(c=v_0\) is the optimal pairing direction — no other constant direction helps.)*

Two elementary inputs:

* **(i)** \(b^{\!\top}\!Mb\le\lambda_{\max}(M)|b|^2\le\tfrac12H_1|b|^2\) by (A.7).
* **(ii)** \(\|u_N-b\|_\infty\ge|b|-C_0\). Indeed \(u_N\) is continuous and
  \((2\pi)^{-3}\!\int|u_N|^2=H_0\), so some \(x_*\) has
  \(|u_N(x_*)|\le\sqrt{H_0}=|v_0|\sqrt{\tfrac23T_N}\le C_0\) (Lemma 10 and
  \(T_N\nearrow T_\infty\)); then
  \(\|u_N-b\|_\infty\ge|b-u_N(x_*)|\ge|b|-C_0\).

Therefore
\[
\sup_c\mathcal B(b,c)\;\le\;\sqrt{H_1}\Bigl[\tfrac{1}{\sqrt2}|b|-(|b|-C_0)\Bigr]
=\sqrt{H_1}\Bigl[C_0-\bigl(1-\tfrac{1}{\sqrt2}\bigr)|b|\Bigr],
\]
which is \(>0\) only for \(|b|<(2+\sqrt2)C_0\), and on that range is
\(\le\sqrt{H_1}\,|b|/\sqrt2\le(1+\sqrt2)C_0\sqrt{H_1}\). Squaring and
inserting \(H_1=\tfrac23|v_0|^2S_N\), \(C_0^2=\tfrac23|v_0|^2T_\infty\) gives
the stated \(C_1S_N\). ∎

**Two corollaries worth recording.**

* The admissible \(b\) are *bounded*: \(|b|<(2+\sqrt2)C_0\). But
  \(|u_N(0)|=\tfrac23S_N|v_0|\to\infty\). So the paper's natural sweeping
  constant is excluded outright for all large \(N\) — the split with
  \(b=u_N(0)\) yields a **negative**, hence vacuous, lower bound. Sharp
  version of the obstruction constant, using the asymptotics of §A.1
  (**NUMERICAL** for the value of \(\lambda\), **PROVEN** as an inequality
  via (A.4)):
  \[
  \frac{\|u_N-u_N(0)\|_\infty\|\nabla u_N\|_2}
       {\tfrac23S_N\|(v_0\cdot\nabla)u_N\|_2}
  \;\longrightarrow\;\sqrt{\frac{2/3}{2/15}}=\sqrt5=2.2360679\ldots
  \]
  matching the measured column \(2.302\to2.250\) of §A.3. The elementary
  substitute of (A.4) already gives \(\ge\sqrt{8/3}=1.633>1\).
* The failure is **not** repaired by any Hölder exponent (§A.3 table), by
  any pairing direction (Cauchy–Schwarz in \(M\)), or by any constant shift
  (input (ii)).

### A.6 Verdict on JOB A

**JOB A: FAILED, and provably so — stopping here, as instructed.**

The failure is located precisely: it is *not* in (T1), (T2), or the exact
evaluation (A.1)/(A.3) — those are all proven and (T2) is 91 %-sharp. It is
in step (A.5): **a constant sweeping vector cannot resolve an \(O(1/N)\)
concentration**, and Theorem A shows every estimate in that family is capped
at \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\gtrsim N\), a full two powers below
the target.

**Scope of Theorem A (stated honestly).** It closes: pairing against
\((c\cdot\nabla)u\) for *constant* \(c\), splitting off a *constant* \(b\),
and bounding \(R_{b,c}\) by Hölder against \(\|\nabla u\|_{L^q}\). It does
**not** close: pairing against a general divergence-free test field
\(\psi\) (\(\langle\mathbb P(u\cdot\nabla u),\psi\rangle=\langle u\cdot\nabla u,\psi\rangle\)
for any divergence-free \(\psi\)), which is precisely route (R4) of the
rescaling programme. Theorem A is in fact a positive argument *for* (R4):
the correct test field must live at the concentration scale
\(\psi_N(x)=N^{3/2}\Psi(Nx)\), not be a constant-coefficient derivative of
\(u\) itself. The remaining burden for L\* therefore sits entirely in the
rescaling route, and (by §B) it is a much lighter burden than L\* as stated.

**Epilogue (2026-08-02): the prediction of this paragraph was correct.**
The rescaling route (R4) was carried out in
[`lstar_proof_main.md`](lstar_proof_main.md) — with exactly the test field
\(\psi_N(x)=N^{3/2}\Psi(Nx)\) named here — and delivered not the lighter
burden of §B but the **full** \(N^3\) with the sharp exponent, for the
smoothly truncated family. Theorem A was therefore a *productive*
negative result: it eliminated the wrong instrument and pointed at the
right one.

---

## §B. JOB B — the weakened hypothesis, and its sufficiency

### B.1 What the Proposition actually consumes

Trace the exponent through the published proof. Only **one** step uses the
number \(3\) in "\(c_0N^3\)": the conversion
\(K(u_N)\ge c_0N^3/H_1^2\ge c_0''N_0^2\). The two structural facts
\(H_1\asymp S_N\asymp N\) and \(N_0^2=S_N/T_N\asymp N\) mean that \(N^3\)
gives \(K\gtrsim N_0^2\), i.e. an *exponential* \(\Phi(s)\gtrsim e^s\). But the
Osgood integral \(\int^\infty ds/\Phi\) converges for vastly slower \(\Phi\).
The correct statement of what is needed is therefore:

> **Hypothesis (L\*-weak).** There exist \(a>0\), \(c_0>0\) and \(N_1\) such
> that \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^{2+a}\) for all \(N\ge N_1\).

and, sharper still (this is the exact threshold):

> **Hypothesis (L\*-min).** There exist \(c_0>0\), \(N_1\), and a
> nondecreasing \(g:[1,\infty)\to(0,\infty)\) with
> \(\displaystyle\int^\infty\frac{dt}{t\,g(t)}<\infty\),
> such that \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^2g(N)\) for all
> \(N\ge N_1\).

Implications: \(\text{(L\*)}\Rightarrow\text{(L\*-weak)}\) with \(a=1\)
(since \(N^3=N^2\cdot N\)), and \(\text{(L\*-weak)}\Rightarrow\text{(L\*-min)}\)
with \(g(t)=t^a\) (for which \(\int^\infty dt/t^{1+a}<\infty\)). All three are
**OPEN**; §B.2 proves that the weakest of them already suffices.

### B.2 Proposition with the exponent carried through

> **Proposition (weakened; PROVEN).** Assume (L\*-min). If \(\Phi>0\) is
> nondecreasing and \(K(u)\le\Phi(\log N_0^2(u))\) for every zero-mean
> divergence-free real trigonometric field \(u\), then there are \(c>0\),
> \(\kappa>0\) and \(s_0\) with
> \(\Phi(s)\ge c\,g(\kappa e^{s})\) for all \(s\ge s_0\), and hence
> \(\int^\infty ds/\Phi<\infty\).
> Under (L\*-weak) this specialises to \(\Phi(s)\ge c\,e^{a's}\) with
> \(a'=a\), the **same** exponent.

*Proof.* Fix \(N\ge\max(N_1,8)\) and abbreviate \(s_N=\log N_0^2(u_N)\).

**Step 1 (\(K\) from below, in \(N\)).** By Lemma 10, \(H_1(u_N)=\tfrac23|v_0|^2S_N\)
exactly; by Lemma 11, \(S_N\le432N\) for \(N\ge8\). Hence
\(H_1\le288|v_0|^2N\) and
\[
K(u_N)=\frac{\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2}{H_1^2}
\;\ge\;\frac{c_0N^2g(N)}{82944\,|v_0|^4N^2}
\;=\;c_1\,g(N),\qquad c_1=\frac{c_0}{82944\,|v_0|^4}.
\tag{B.1}
\]
*(This is the only step where the original "\(N^3\)" was consumed; \(N^2\) is
exactly the amount the two factors of \(H_1\) eat, so \(g\) is pure surplus.)*

**Step 2 (\(N\) from below, in \(N_0^2\)).** \(N_0^2(u_N)=S_N/T_N\) (Lemma 10),
with \(S_N\le432N\) (Lemma 11) and \(T_N\ge T_1=6\) (the six unit vectors
\(\pm e_i\), each contributing \(1\)). Hence \(N_0^2\le72N\), i.e.
\(N\ge\tfrac{1}{72}N_0^2=\tfrac{1}{72}e^{s_N}\). Since \(g\) is nondecreasing,
(B.1) gives
\[
K(u_N)\;\ge\;c_1\,g\!\bigl(\tfrac1{72}e^{s_N}\bigr).
\tag{B.2}
\]

**Step 3 (dyadic interpolation; the bounded-gap input).** Put
\(N_j=2^jN_2\) with \(N_2=\max(N_1,8,N_*)\), \(N_*\) the threshold in
Lemma 11. Lemma 11 gives \(0<c_-\le s_{N_{j+1}}-s_{N_j}\le c_+<\infty\) for
\(j\ge j_0\); so \((s_{N_j})_{j\ge j_0}\) is strictly increasing and (gaps
bounded **below**) unbounded. Set \(s_0=s_{N_{j_0}}\). Given \(s\ge s_0\),
choose \(j\ge j_0\) with \(s\in[s_{N_j},s_{N_{j+1}}]\). Applying the
hypothesis to \(u=u_{N_j}\) and using monotonicity of \(\Phi\) and of \(g\),
together with \(s_{N_j}\ge s-c_+\):
\[
\Phi(s)\;\ge\;\Phi(s_{N_j})\;\ge\;K(u_{N_j})
\;\overset{(B.2)}{\ge}\;c_1\,g\!\bigl(\tfrac1{72}e^{s_{N_j}}\bigr)
\;\ge\;c_1\,g\!\bigl(\kappa e^{s}\bigr),\qquad
\kappa:=\tfrac{1}{72}e^{-c_+}.
\tag{B.3}
\]

**Step 4 (convergence of the Osgood integral).** Substituting \(t=\kappa e^s\)
(so \(dt=t\,ds\)),
\[
\int_{s_0}^{\infty}\frac{ds}{\Phi(s)}
\;\le\;\frac1{c_1}\int_{s_0}^{\infty}\frac{ds}{g(\kappa e^{s})}
\;=\;\frac1{c_1}\int_{\kappa e^{s_0}}^{\infty}\frac{dt}{t\,g(t)}\;<\;\infty
\]
by the hypothesis on \(g\). On \([0,s_0]\), \(\Phi\ge\Phi(0)>0\), so that
piece is finite too. Hence \(\int^\infty ds/\Phi<\infty\).

**Specialisation to (L\*-weak).** \(g(t)=t^a\): (B.3) reads
\(\Phi(s)\ge c_1\kappa^ae^{as}=c\,e^{as}\) for \(s\ge s_0\), with
\(c=c_0(72\,e^{c_+})^{-a}/(82944|v_0|^4)\); and
\(\int_{s_0}^\infty ds/\Phi\le e^{-as_0}/(ac)<\infty\). ∎

**Threshold is sharp in \(g\).** \(\int^\infty dt/(t\,g(t))<\infty\) is not
merely sufficient for Step 4, it is what Step 4 *is* after the substitution;
so no weaker growth of \(g\) can be accommodated by this argument. Concretely:

| \(g(N)\) | \(\int^\infty\frac{dt}{tg(t)}\) | Proposition? |
|---|---|---|
| \(N^a\), any \(a>0\) | finite | **yes** |
| \((\log N)^{1+\varepsilon}\), \(\varepsilon>0\) | finite | **yes** |
| \(\log N\cdot(\log\log N)^{1+\varepsilon}\) | finite | **yes** |
| \(\log N\) | \(=\infty\) | **no** |
| \(1\) (i.e. only \(K\ge c\)) | \(=\infty\) | **no** |

So the hypothesis needed is barely more than "\(K(u_N)\to\infty\) faster than
\(\log N_0^2\) times a double-log".

### B.3 Exact text to replace (L\*) in the paper — **WITHDRAWN, do not apply**

> **Withdrawn (2026-08-02).** This subsection prescribed keeping the sharp
> family and weakening (L\*). The paper instead **switched the family** to
> the smooth truncation and **proved** the full \(N^3\) bound, which makes
> the Proposition unconditional and (L\*-weak) moot. The two plans are
> mutually exclusive. What follows is retained only as a record of the
> road not taken. The one piece that *was* adopted is the observation
> \(T_N\ge T_1=6\) noted at the end.

In [`../theorem_statement.md`](../theorem_statement.md), replace the
Hypothesis (L\*) block by:

> **Hypothesis (L\*-weak)** *(open; the sharp form \(a=1\) is exactly
> certified at \(N\le8\) and measured to \(N=32\); see the certificate
> appendix)*: there are \(a>0\), \(c_0>0\) and \(N_1\) with
> \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^{2+a}\) for all \(N\ge N_1\).
>
> *(Only the following consequence is used: \(K(u_N)\ge c_1N^a\). Any
> nondecreasing \(g\) with \(\int^\infty dt/(t g(t))<\infty\) in place of
> \(N^a\) suffices verbatim — e.g. \(N^2(\log N)^{1+\varepsilon}\); see
> `lstar/lstar_direct_route_and_weakening.md` §B. The trivial upper bound
> \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\le\|u_N\|_\infty^2H_1
> =\tfrac{8}{27}|v_0|^4S_N^3\asymp N^3\) shows \(N^3\) is the maximum
> possible, so (L\*) asserts saturation of Lemma 7 whereas (L\*-weak)
> asserts only a fixed power of surplus over the \(H_1^2\) normalisation.)*
>
> **Proposition.** Assume (L\*-weak). If \(\Phi\) is nondecreasing and
> \(K(u)\le\Phi(\log N_0^2(u))\) for every zero-mean divergence-free real
> trigonometric field \(u\), then \(\Phi(s)\ge c\,e^{as}\) for all large
> \(s\), and hence \(\int^\infty ds/\Phi<\infty\): no Osgood-admissible
> \(\Phi\) satisfies a uniform pointwise bound.

In [`../complete_proof.md`](../complete_proof.md), the *Proof of the
Proposition* is replaced by Steps 1–4 of §B.2 above. Lemmas 9, 10, 11 are
untouched and are the only inputs, exactly as before; the constant
dependence becomes \(c=c_0(72e^{c_+})^{-a}/(82944|v_0|^4)\), still explicit.

The Constant-dependence section should also gain \(T_N\ge T_1=6\) (used in
Step 2 in place of the vaguer "\(T_{2N}/T_N\to1\)" remark).

### B.4 How much headroom this buys (NUMERICAL context, not a proof)

* Ceiling (**PROVEN**, Lemma 7 + Lemma 10):
  \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\le\tfrac{8}{27}|v_0|^4S_N^3\), i.e.
  \(K(u_N)\le\tfrac23S_N\). L\* asserts this is saturated.
* Measured (**NUMERICAL**, §A.3): \(K(u_N)\approx0.316\,N\), i.e.
  \(K/N_0^2\) rising \(0.259\to\approx0.43\) over \(N=4\dots32\) —
  consistent with saturation, hence with \(a=1\).
* Needed for the Proposition (**PROVEN**, §B.2): only
  \(K(u_N)\gtrsim(\log N)(\log\log N)^{1+\varepsilon}\).

The measured growth exceeds the requirement by a factor \(\sim N/\log N\).
This is why the weakening is a *free strict improvement*: it does not weaken
the Proposition's conclusion at all, and it moves the open problem from
"prove a sharp saturation constant" to "prove any super-logarithmic growth".

---

## §C. Consolidated status ledger

| Claim | Status |
|---|---|
| (T1) sweeping identity, (T2) pairing bound | **PROVEN** |
| (A.1) \(\|(v_0\!\cdot\!\nabla)u_N\|_2^2=|v_0|^4\sum(\cos^2-\cos^4)/|k|^2\) | **PROVEN** |
| \(\sum\cos^2\theta_k/|k|^2=S_N/3\) exactly | **PROVEN** (Lemma 9) |
| (A.2) \(3A_N+6B_N^{\,*}=S_N\); (A.3) closed form with \(\Delta_N,\rho\) | **PROVEN**, exact-verified \(N\le32\) |
| (A.4) \(\|(v_0\!\cdot\!\nabla)u_N\|_2^2\le\tfrac14|v_0|^4S_N\); lower bd \(\ge c_-|v_0|^4S_N\) | **PROVEN** (lower bd modulo routine lattice count) |
| \(\Delta_N=O(1)\); hence \(\|(v_0\!\cdot\!\nabla)u_N\|^2\sim\tfrac2{15}|v_0|^4S_N\) | **NUMERICAL** (provable via standard equidistribution; not needed) |
| (A.6) first piece \(=\tfrac23S_N\|(v_0\!\cdot\!\nabla)u_N\|_2^2>0\) | **PROVEN** |
| (A.7) \(\lambda_{\max}(M)\le\tfrac12H_1\); \(\operatorname{tr}M=H_1\) | **PROVEN** |
| \(M\) spectrum \(\to(\tfrac15,\tfrac25,\tfrac25)H_1\) | **NUMERICAL** |
| remainder \(\approx-0.607\times\)first, non-decaying; \(Q_{v_0}\approx1090N^2\) | **NUMERICAL** (binary64) |
| \(\|u_N\|_\infty=|u_N(0)|\) exactly | **NUMERICAL** (12 digits) |
| **Theorem A**: sweeping split capped at \(\|\mathbb P(u\!\cdot\!\nabla u)\|^2\gtrsim N\) | **PROVEN**, elementary |
| \(c=v_0\) is the optimal constant pairing direction | **PROVEN** |
| Sharp obstruction constant \(\sqrt5\) | **NUMERICAL** (\(\ge\sqrt{8/3}>1\) is **PROVEN**) |
| **§B.2 Proposition under (L\*-weak) / (L\*-min)** | **PROVEN**, but **MOOT** — the full \(N^3\) is now proven for the smooth family |
| (L\*), (L\*-weak), (L\*-min) for the **smoothly** truncated family | **PROVEN** (`lstar_proof_main.md` Thm 7.1(3)), sharp exponent, uniform in \(\chi,v_0\); constants \(c_0,N_*\) non-effective |
| (L\*), (L\*-weak), (L\*-min) for the **sharply** truncated family | **OPEN — and used by nothing** |

## §D. Reproduction

Scripts live beside this note in `docs/paper_lambda_dichotomy/lstar/`
(binary64 + exact lanes; numpy/`fractions` only, no scipy). Run from that
directory with the repo venv, e.g.
`.venv/Scripts/python.exe lstar_direct.py 4 8 16 24`. Contents:
`lstar_direct.py` (full-lattice spectral pipeline, alias-free grid
\(n=4N+2\); table of §A.3), `lstar_exact.py` (`Fraction` verification of
(A.1)–(A.3), all assertions pass \(N\le32\)), `lstar_split.py` (remainder
audits), `lstar_holder.py` (Hölder scan and \(\operatorname{spec}M\)),
`lstar_conc.py` (\(Q\)-density concentration profile, \(\|u_N\|_\infty\)).
Cross-validation against the repository: \(K(u_4)=0.7884\), \(K(u_8)=2.0372\)
and pairing fractions \(0.9210\)–\(0.9119\) reproduce
`outputs/verification_sprint_v1/osgood_gate/` and the paper's certificate
appendix.
