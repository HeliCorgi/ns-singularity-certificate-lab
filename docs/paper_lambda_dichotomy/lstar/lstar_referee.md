# Adversarial referee report on the (L\*) dossier

Referee brief: **break the main route** (`lstar_proof_main.md`), cross-checking
against `lstar_direct_route_and_weakening.md`, `lstar_numerical_support.md`,
`../theorem_statement.md`, `../complete_proof.md`.

Every load-bearing computation below was re-derived by hand or recomputed.
Two independent verification scripts written for this report live beside it:

* [`referee_tau_check.py`](referee_tau_check.py) — direct 3-D quadrature of the
  convolution \(\tau=h*h\), testing Theorem 6.5 with **no** Fourier convention
  in the loop.
* [`referee_riemann_check.py`](referee_riemann_check.py) — direct lattice
  evaluation of \(u_N(y/N)/N\), testing Lemma 4.1 / Theorem 4.2.

**Bottom line up front: I could not break the main route.** The result is a
genuine proof, not a numerically-supported reduction. The defects found are one
false (and unused) proposition, an integration problem with the rest of the
paper, a false "all constants are explicit" claim, and cosmetics.

---

## 0. What I verified affirmatively (so the reader knows what was actually checked)

| Item | Location | Verification performed |
|---|---|---|
| Conventions, \(\mathcal F\widetilde G=(2\pi)^3G\), Plancherel | §0.1 | re-derived; self-consistent throughout |
| Lemma 9′/9″ (radial weight, no positivity) | §1.2 | re-proved; correct, and strictly more general than paper Lemma 9 |
| Theorem 1.3 exact laws | §1.3 | re-derived; matches numerics doc Table A (residual = rational 0) |
| Lemmas 2.1–2.3, Prop 2.4 | §2 | all arithmetic recomputed: shell count \((2^{j+2}-1)^3<2^{3j+6}\) ✓, \(S_N\le128N\) ✓, \(0.28868M-0.0625M=0.22618M\), \((0.226M)^3>0.01154M^3\ge M^3/87\) ✓, \(N_0^2\in[N/44544,\tfrac{64}{3}N]\) ✓ |
| \(\widehat{u_N}(k)=N^{-2}F(k/N)\) | (3.1) | correct; \(P_k=P_{k/N}\) |
| Theorem 3.2 (1)–(5) | §3.2 | \(\|F\|_1\le4\pi\|v_0\|\), \(\|\partial^\beta V\|_\infty\le4\pi\|v_0\|/(|\beta|+1)\), \(V(0)=\tfrac{8\pi}{3}(\int\chi)v_0\) — all recomputed ✓ |
| \(\widetilde h=\pi^2W\) | Thm 3.3(a) | re-derived from scratch: \(c_{3,1}=4\pi\), \(c_{3,2}=2\pi^2\), \(c_{3,-1}=-8\pi\) (\(\Gamma(-\tfrac12)=-2\sqrt\pi\)) ✓; \(y_iy_j|y|^{-3}=\delta_{ij}|y|^{-1}-\partial_i\partial_j|y|\) ✓; \(\mathcal F[\mathcal O]=h\) ✓; \(\widetilde h=(2\pi)^3\mathcal O v_0=\pi^2W\) ✓ |
| Remainder \(\rho=O(|y|^{-M})\) | Thm 3.3(b) | correct (\(|\beta|\ge2\Rightarrow\partial^\beta G\in L^1\)) |
| Lemma 4.1 exact identity | §4.1 | correct, including \(\Phi_y(0)=F(0)=0\) bookkeeping |
| Theorem 4.2 rate | §4.2 | structure correct; **numerically corroborated**: with the numerics doc's \(\chi\), \(v_0=e_3\), \(y=(1,0.5,-0.3)\), successive differences of \(u_N(y/N)/N\) at \(N=8,16,32,64\) are \(3.714\text{e-}1,\,1.857\text{e-}1,\,9.285\text{e-}2\) — halving **exactly**, i.e. \(E_N=c/N\) with no logarithm at all (Theorem 4.2 is not tight, but is valid) |
| Lemma 5.2 (\(\psi_N\) div-free, zero-mean, \(\|\psi_N\|_2=(2\pi)^{-3/2}\|\Psi\|_{L^2(\mathbb R^3)}\)) | §5.1 | all three re-proved ✓. Zero mean via \(\int\Psi_i=\int\nabla\cdot(y_i\Psi)=0\) is correct and the curl representation is indeed unnecessary |
| Theorem 5.3 / Cor 5.4 factor bookkeeping | §5.2 | every \((2\pi)^{-3}\) and every power of \(N\) recomputed: \(\|\mathbb P(u\!\cdot\!\nabla u)\|_2^2\ge(2\pi)^{-3}N^3|I_\Psi+\mathcal E_N|^2/\|\Psi\|^2_{L^2(\mathbb R^3)}\) ✓ |
| Lemma 6.1 symbol \(M_i=i\zeta_jT_{ij}\) | §6.1 | re-derived; the replacement \(\eta_jF_j(\zeta-\eta)\to\zeta_jF_j(\zeta-\eta)\) is legitimate ✓; Young's inequality gives \(T\in L^1\cap L^2\) ✓ |
| **Theorem 6.5** \(\zeta\times(\tau\zeta)=\tfrac{3\pi^3}{8}\|v_0\|^2\tfrac{\zeta_3}{|\zeta|}(\zeta\times e_3)\) | §6.2 | **confirmed to 14 digits** by `referee_tau_check.py` at \(\zeta=(0.6,0,0.8),(0.3,0.4,0.5),(1,1,2)\); rel. err. \(3.6\text{e-}14,\ 1.2\text{e-}14,\ 6.6\text{e-}15\). This is a *fully independent* test: it evaluates \(\int h_i(\eta)h_j(\zeta-\eta)d\eta\) by two-patch spherical quadrature and never uses a Fourier transform, so it closes the "finite-part regularisation" risk in §6.2 |
| Lemma 6.5′ (elementary curl computation) | §6.2 | re-derived by hand line by line: \(\nabla\cdot W=0\), \(\omega=2(e_3\times y)/r^3\), \(W\times\omega=\tfrac{4y_3e_3}{r^4}-\tfrac{2y}{r^4}-\tfrac{2y_3^2y}{r^6}\), \(\nabla\times(W\times\omega)=\tfrac{12y_3(e_3\times y)}{r^6}\) ✓ |
| Fourier↔physical consistency check | §6.2 | re-derived: \(-4\pi\partial_b\partial_3r^{-2}=8\pi\delta_{b3}r^{-4}-32\pi y_3y_br^{-6}\), \(\epsilon_{ab3}\delta_{b3}=0\), total \(-12\pi^4y_3(e_3\times y)/r^6\) ✓ — matches Lemma 6.5′ exactly. **Third, cross-document confirmation:** the numerics doc's independent axisymmetric reduction gives \(Z_\infty=-12\pi^4\|v_0\|^2/r^4\), which is exactly this field written as \(Z(r)\sin\theta\cos\theta\,e_\phi\) |
| Lemma 6.6 \(|T-\tau|\le64\pi\|v_0\|^2\) on \(|\zeta|\le\tfrac18\) | §6.3 | re-derived ✓ (\(|\eta|\le\tfrac14\Rightarrow|\zeta-\eta|\le\tfrac38<\tfrac12\); \(4\|v_0\|^2\cdot4\pi\cdot4=64\pi\|v_0\|^2\)) |
| Theorem 6.7 threshold \(\pi^2/2048\) | §6.3 | \(\tfrac{3\pi^3}{32}/(192\pi)=\tfrac{3\pi^2}{6144}=\tfrac{\pi^2}{2048}\approx0.00482<\tfrac18\) ✓ |
| Theorem 7.2 arithmetic | §7.2 | \(288^2=82944\); \(65536\cdot64=4194304\); \(44544\cdot64=2850816\), \(/3=950272\); \(2^{20}=1048576>950272\); \(c_-=\log(1.10345)=0.0985\), \(c_+=\log(2^{20}\!\cdot\!950272)=27.63\) — all ✓ |
| JOB B (§B.2 of the weakening doc) | — | \(\tfrac23\cdot432=288\), \(288^2=82944\), \(432/6=72\), \(\kappa=e^{-c_+}/72\), \(t=\kappa e^s\Rightarrow ds=dt/t\) — the substitution is exact and the threshold \(\int^\infty dt/(t\,g(t))<\infty\) is correct ✓ |
| Theorem A (JOB A no-go) | — | re-derived: \(\sup_c b^\top Mc/\sqrt{c^\top Mc}=\sqrt{b^\top Mb}\); (A.7) \(\lambda_{\max}(M)\le\tfrac12H_1\); \((1+\sqrt2)^2=3+2\sqrt2\); \(C_1=(3+2\sqrt2)\tfrac49\|v_0\|^4T_\infty\) ✓ |

The two numerical documents, which at face value disagree by a factor \(\sim\!250\)
in \(\|\mathbb P(u\!\cdot\!\nabla u)\|^2/N^3\), **reconcile**: see MINOR-1. After
rescaling, main-doc \(K/N_0^2\approx0.468\) vs numerics-doc \(0.454\) at \(N=40\)
— a 3 % difference attributable to the different \(\chi\). This is a real
cross-validation, and I count it in the proof's favour.

---

## CRITICAL

**None.** Specifically, and against my brief:

* the Poisson-summation identity is **not used** anywhere in the proof and its
  removal costs nothing (§4 replaces it by an exact identity);
* the \(C^1_{\rm loc}\) statement **is** proven, derivatives included, with an
  explicit rate — §4.2 gives the gradient its own argument (`|Φ'_y|≲|ξ|^{-1}`,
  `|∇Φ'_y|≲c_χ|ξ|^{-2}+R|ξ|^{-1}`, summing to \(O((c_\chi+R)/N)\) with no
  logarithm). It is not asserted, it is derived, and it is the easier half;
* \(\psi_N\) is genuinely divergence-free **and** zero-mean as a torus field,
  its \(L^2(\mathbb T^3)\) normalisation is right and \(N\)-independent, and on
  the fundamental domain exactly one periodisation term survives, so the
  periodisation contributes nothing. \(u_N\) is never periodised;
* **(V-NONDEG) is PROVEN, not numerically supported.** Theorem 6.7 rests on
  Theorem 6.5 (a closed-form identity, which I confirmed to 14 digits by a
  method sharing no code, no convention and no regularisation with the
  document's derivation) plus Lemma 6.6 (four lines of elementary estimation,
  re-derived). The numerics in `lstar_numerical_support.md` §F corroborate it;
  they are not load-bearing for it.

Consequently the correct conclusion is **not** "(L\*) reduced to a numerically
certain statement". It is stronger than that. See the verdict.

---

## MAJOR

### MAJOR-1. Proposition 3.4, equation (3.2): the Abel-regularised Poisson identity is **FALSE**, and contradicts the same Proposition's own divergence proof

`lstar_proof_main.md` §3.4, labelled **PROVEN**, and carried into the §7.3
ledger row "Poisson identity, Abel-regularised | **PROVEN** (3.2), unused".

The Proposition proves divergence by exhibiting sign-definiteness:
\(v_0\cdot V(N(x+2\pi m))\ge\pi^2\|v_0\|^2/(N|x+2\pi m|)-O(|m|^{-4})\), so the
terms of \(\sum_{m\neq0}v_0\cdot V(N(x+2\pi m))\) are **eventually positive**
and the series diverges to \(+\infty\). Abel summation of a series of
eventually-positive terms with divergent sum also diverges (monotone
convergence): \(\lim_{\varepsilon\downarrow0}\sum_m e^{-\varepsilon|m|}a_m=+\infty\)
whenever \(a_m\ge0\) and \(\sum a_m=\infty\). So the right-hand side of (3.2)
is \(+\infty\) in the \(v_0\) direction, while the left-hand side \(u_N(x)\) is
finite. Quantitatively the constant mode blows up like \(\varepsilon^{-2}\),
because the spherical mean of \(W\) is \(\tfrac43v_0/r\neq0\).

The offered proof is also unsound at two points: \(\widetilde G\notin L^1\)
(it decays like \(|y|^{-1}\)), so "\(k\)-th Fourier coefficient
\(=(2\pi)^{-3}\int\widetilde Ge^{-iky}dy=G(k)\) by Fourier inversion" is not
available as written; and "Abel summation of a distributionally convergent
series with continuous sum converges pointwise" is not a theorem.

*Impact:* **none on the result** — §3.4 ends with "nothing below depends on it",
and I confirmed by inspection that §§4–7 never cite (3.2).

*Minimal fix:* delete (3.2) and the ledger row, keeping only the divergence
statement. If a periodisation formula is wanted for orientation, the correct
renormalised one is absolutely convergent:
\[
u_N(x)=N\Big[V(Nx)+\sum_{m\neq0}\big(V(N(x+2\pi m))-V(2\pi Nm)\big)\Big]+\text{const},
\]
since the bracketed differences are \(O(|m|^{-2})\) by Theorem 3.3.

*Related, in `lstar_numerical_support.md` §G:* the row on Poisson summation says
the \(m\neq0\) terms' "ratio to the main term is \(O(1/N)\), which is all the
duality step used". That is not enough — there are infinitely many such terms
and their sum is \(+\infty\), not \(O(1)\). The row's verdict (**FALSE AS
STATED**) is right; its consolation clause is wrong. Delete the clause.

### MAJOR-2. The family switch is not a free "separate integration step": it orphans the paper's certificates and collides with the other author's integration plan

`../theorem_statement.md` and `../complete_proof.md` are written entirely around
the **sharply** truncated family \(\widehat u_N(k)=P_kv_0/|k|^2\) on
\(1\le|k|\le N\). Adopting Theorem 7.2 requires, at minimum:

1. **Proposition family definition** (`theorem_statement.md`, "Proposition")
   must become \(\widehat u_N(k)=\chi(|k|/N)P_kv_0/|k|^2\), with \(\chi\)
   admissible and \(v_0\in\mathbb R^3\setminus\{0\}\) (the \(v_0\in\mathbb Z^3\)
   restriction is unnecessary and should go).
2. **Hypothesis (L\*) block** must be deleted, not merely re-labelled: its
   parenthetical "*(open; exactly certified at \(N\le8\), measured to \(N=32\))*"
   points at sharp-family certificates.
3. **Certificate appendix** (`complete_proof.md` lines 300–311): \(K(u_4)=0.7884\),
   \(K(u_6)\), \(K(u_8)=2.0372\), the 91.18–92.10 % pairing fractions, and
   \(K/N_0^2:0.259\to0.396\) are all sharp-family numbers. After the switch they
   support a hypothesis that no longer appears. They must be re-labelled as
   "corroboration for the *sharp*-family (L\*), which remains open and unused",
   or replaced by the smooth-family exact table (`lstar_numerical_support.md` §B).
4. **Constant dependence section**: \(N/250\le S_N\le432N\) and the dyadic gap
   constants \(c_\pm\) come from Lemma 11 (sharp family). They must be replaced by
   Prop 2.4's \(N/348\le S_N^\chi\le128N\), \(6\le T_N^\chi\le128\)
   (\(N\ge32\)) and \(c_-=0.0985\), \(c_+=27.63\) at \(q=20\). *(No contradiction
   between the two upper bounds — 128 is simply the tighter shell count; both
   are valid.)*
5. **Lemma 9 → Lemma 9′**: the paper's Lemma 9 is stated for \(f\ge0\) on the
   sharp ball; the smooth family needs the general radial-weight version. The
   generalisation is trivial and correct (§1.2), but the paper text must change.
6. **Lemma 10 → Theorem 1.3, Lemma 11 → Prop 2.4.**

Additionally, `lstar_direct_route_and_weakening.md` §B.3 prescribes a *different*
integration: keep the sharp family, replace (L\*) by (L\*-weak), keep Lemmas
9/10/11 "untouched". **The two plans are mutually exclusive.** If the smooth
family is adopted the Proposition becomes unconditional and (L\*-weak) is moot;
§B should then be demoted to a remark ("the Proposition needs far less than the
exponent 3; recorded in case one insists on the sharp family"), which is still
worth keeping because §B.2's Step 2 (\(T_N\ge T_1=6\)) is a genuine tightening of
Lemma 11.

### MAJOR-3. "All constants are explicit" becomes false

`theorem_statement.md` §"Constant dependence" asserts all constants are explicit.
Theorem 7.1(3) writes \(c_0=|I_\Psi|^2/(2(2\pi)^3\|\Psi\|^2)\) "for any admissible
\(\Psi\) with \(I_\Psi\neq0\), which exists by Theorem 6.7 + Lemma 6.2" — but
Lemma 6.2 obtains \(\Psi\) from **density of \(C^\infty_{c,\sigma}\) in
\(L^2_\sigma\)**, a pure existence argument. Neither \(\Psi\), nor \(R\), nor
\(I_\Psi\), nor \(N_*\) is produced. The proof is therefore **non-effective**.

This is not a gap in the proof, but it *is* a false advertisement in the paper's
constant-dependence claim, and the "constructive variant" remark after Lemma 6.2
understates the cost: Theorem 6.7 certifies non-vanishing only on
\(|\zeta|<\pi^2/2048\approx4.8\times10^{-3}\), so an explicitly constructed
\(\Psi\) has Fourier support at that scale, hence spatial radius \(R\gtrsim10^3\),
hence \(N_*\ge8R\gtrsim10^4\) and a \(c_0\) many orders below the measured
\(\approx8.4\|v_0\|^4\). *Minimal fix:* state Theorem 7.1(3) as "there exists
\(c_0>0\)", and amend "Constant dependence" to "explicit except for \(c_0,N_*\)
in the Proposition, which are non-effective".

---

## MINOR

**MINOR-1 (numerical inconsistency, resolved).** The §7.5 table of
`lstar_proof_main.md` does not state \(v_0\), and is *inconsistent* with the
\(v_0=(1,2,-1)\) declared in §1.3. From the table, \(H_1=\sqrt{\|\mathbb P\|^2/K}
=229.6\) at \(N=40\), so \(S_N^\chi/N=8.6\) if and only if \(\|v_0\|=1\); with
\(\|v_0\|^2=6\) one would get \(S_N^\chi/N=1.44\), impossible (the true value is
\(4\pi\int_0^1\chi^2\approx8.4\)). So §7.5 used a **unit** \(v_0\). Once stated,
the two documents agree: rescaling by \(\|v_0\|^4=196\) gives
\(\|\mathbb P\|^2/N^3\approx1601\) vs the numerics doc's 1881 at \(N=40\)
(different \(\chi\)), and the \(\|v_0\|\)-independent quantities agree closely:
\(K=9.918\) vs \(10.309\), \(K/N_0^2\approx0.468\) vs \(0.454\). *Fix:* state
\(v_0\) in the §7.5 caption.

**MINOR-2 (over-strong hypothesis).** Definition 0.1 demands \(\chi\in C^\infty\),
but the proof uses only \(\chi\in C^1\): \(c_\chi=1+\|\chi'\|_\infty\) is the sole
smoothness input (§4.2); §6 uses only \(0\le\chi\le1\) and \(\chi\equiv1\) on
\([0,\tfrac12]\); Theorem 3.3(b), the only consumer of higher derivatives, feeds
**only** the deleted Proposition 3.4. Relaxing Definition 0.1 to \(C^1\) is free
and removes the compatibility caveat opened in `lstar_numerical_support.md` §0
(its rational \(\chi\) is \(C^4\), hence admissible after the relaxation, so the
exact-rational lanes A/B/D are evidence for the theorem as proven, not for a
weaker variant).

**MINOR-3 (constants in Theorem 4.2).** (a) \(|\nabla(P_\xi v_0)|\le2\|v_0\|/|\xi|\)
is too good: differentiating \(\xi_i\xi_j/|\xi|^2\) gives four terms, so
\(4\|v_0\|/|\xi|\) is the honest bound. (b) The sup on the \(\sqrt3/N\)-neighbourhood
of \(A_j\) is taken at radius \(\ge2^{-j-2}\), not \(2^{-j-1}\), costing another
factor 8. (c) The sentence "(since \(2^{-J-1}\le\sqrt3/N\cdot\)const, absorbed in
\(A_0\) by enlarging the core constant)" is the covering argument and should be
written out: cells not in the core have \(|\xi|>2\sqrt3/N\) and \(2^{-J-1}\ge2\sqrt3/N\),
so the annuli \(0\le j\le J\) do cover them. (d) The \(k=0\) cell deserves an
explicit word: \(\Phi_y(0):=0\) but \(\int_{Q_0}\Phi_y\neq0\); it is inside the
core and bounded there. **None of these changes the \(O(N^{-1}\log N)\) rate**,
which my numerics show is in any case not tight — the true \(E_N\) is \(\Theta(1/N)\).

**MINOR-4 (Lemma 6.2).** Should state \(\Psi\) **real**, as Definition 5.1
requires. Immediate, since \(\mathbb P(V\cdot\nabla V)\) is real and real
\(C^\infty_{c,\sigma}\) fields are dense in real \(L^2_\sigma\); but as written
the lemma delivers a complex \(\Psi\) and \(\psi_N\) would not be a real torus
field.

**MINOR-5 (rigour labels).** Two steps in §6 are "PROVEN modulo standard" in
substance though not so labelled: (a) Lemma 6.1's use of \(\widetilde A\widetilde B
=\widetilde{A*B}\) for \(F\in L^1\) is fine, but Theorem 6.5's
\(\tau=(2\pi)^{-3}\mathcal F[\widetilde h_i\widetilde h_j]\) applies it to
\(h\notin L^1\); (b) the transforms of \(r^{-4},r^{-6}\) are finite-part
continuations. Both are standard and both are *empirically closed* by
`referee_tau_check.py`, which evaluates \(\tau\) as a convergent convolution
integral and reproduces the boxed formula to \(10^{-14}\). I recommend adding
one sentence splitting \(h=h\mathbb 1_{|\xi|\le1}+h\mathbb 1_{|\xi|>1}\) (the
first in \(L^1\cap L^p\), \(p<3/2\); the second in \(L^p\cap L^\infty\), \(p>3/2\))
to make the convolution theorem literal, and citing the numerical check.

**MINOR-6 (numerics doc, §G).** \(\widehat{\psi_N}(k)=(2\pi)^{-3}N^{-3/2}
\mathcal F[\Psi](k/N)\); the \((2\pi)^{-3}\) is dropped in the row's display. The
main document does not use this route, so nothing propagates.

**MINOR-7 (`complete_proof.md` Lemma 11).** Its dyadic-gap proof leans on the
asymptotic "\(T_{2N}/T_N\to1\)". The weakening doc's observation \(T_N\ge T_1=6\)
makes the gap bound elementary and \(N\)-uniform; adopt it regardless of which
family survives. (Prop 2.4 of the main doc already does exactly this.)

**MINOR-8 (§7.5 conjecture).** The stated limit
\(N^{-3}\|\mathbb P(u_N\!\cdot\!\nabla u_N)\|^2\to(2\pi)^{-3}\|\mathbb P(V\!\cdot\!\nabla V)\|_{L^2(\mathbb R^3)}^2\)
is correctly flagged unproven; note it is the Cauchy–Schwarz-optimal version of
Corollary 5.4, so it is exactly the statement "the localised test field is
asymptotically extremal". Worth saying, since it explains why \(a=1\) is sharp
and cannot be improved.

---

## Answers to the six brief items

1. **Poisson summation.** Convention and constants are internally consistent
   (\(\mathcal F\widetilde G=(2\pi)^3G\), so the periodisation's \(k\)-th
   coefficient *would* be \(G(k)\)). The hypotheses do **not** justify it: \(F\)'s
   \(|\xi|^{-2}\) singularity forces \(V\sim|y|^{-1}\) with a sign-definite
   profile, so no summation method converges (MAJOR-1). The document's own
   verdict "R3 is false as stated" is correct; its replacement (3.2) is also
   false. **Load-bearing: no.**
2. **Decay/regularity of \(V\) and \(C^1_{\rm loc}\).** All of Theorem 3.2 and
   Theorem 3.3 check out. The derivative statement is **proven, not asserted**:
   §4.2 gives \(E'_N\) an independent argument with a strictly better (log-free)
   rate. Numerically \(E_N=\Theta(1/N)\), consistent with and better than the
   claim.
3. **Duality step.** \(\psi_N\) is div-free and mean-zero (the mean-zero proof
   via \(\int\nabla\cdot(y_i\Psi)=0\) is correct and the curl representation is
   genuinely unnecessary); \(\|\psi_N\|_2=(2\pi)^{-3/2}\|\Psi\|_{L^2(\mathbb R^3)}\)
   is right and \(N\)-independent; errors are \(O(\log N/N)\) against an \(O(1)\)
   main term; the test field's periodisation contributes exactly zero on the
   fundamental domain, and \(u_N\) is never periodised. **No defect found.**
4. **(V-NONDEG).** **PROVEN.** Theorem 6.5 verified independently to \(10^{-14}\);
   Lemma 6.5′ re-derived by hand; Lemma 6.6 elementary and correct; the
   threshold \(\pi^2/2048\) checks. The proof is uniform in \(\chi\) and \(v_0\)
   as claimed. This is *not* a numerical nonvanishing.
5. **Weakening (JOB B).** Exponent bookkeeping is correct end to end
   (\(288^2=82944\); \(432/6=72\); \(\kappa=e^{-c_+}/72\); \(t=\kappa e^s\Rightarrow
   ds=dt/t\); \(c=c_0(72e^{c_+})^{-a}/(82944\|v_0\|^4)\)). The claimed sharpness of
   the threshold \(\int^\infty dt/(tg(t))<\infty\) is sharpness *of this argument*,
   correctly qualified. Note it is now **moot** if the smooth family is adopted.
6. **Collateral damage.** See MAJOR-2 and MAJOR-3. Nothing in the **Main
   Theorem** (a)–(e) touches the family, so parts (a)–(e) are unaffected.
   Lemma 11's constants are not *wrong*, merely superseded. The exact
   sharp-family certificates become orphaned corroboration for a now-unused
   hypothesis and must be re-labelled, not deleted.

---

## Verdict

> **(i) — (L\*) is PROVEN for the smoothly truncated family, with the sharp
> exponent \(N^3\), and the paper's Proposition (static no-go) becomes
> unconditional.**

with three statements that must accompany any claim made on this basis:

* the **paper's literal Hypothesis (L\*)**, for the sharply truncated family,
  is **untouched and still open**. What has been proven is (L\*) for a
  *different, equally admissible* family — which is all the Proposition ever
  needed, because its hypothesis quantifies over *all* real zero-mean
  divergence-free trigonometric fields. The honest headline is "the Proposition
  is unconditional", not "(L\*) is proven";
* the constants \(c_0,N_*\) are **non-effective** (MAJOR-3);
* Proposition 3.4/(3.2) must be **deleted** before publication (MAJOR-1), and
  the paper-level integration (MAJOR-2) is substantive editorial work, not a
  formality.

I record explicitly that I attempted to break the route at each of the six
points in my brief and failed at all six; the single false statement I found is
inert. The strongest positive evidence is that the one identity everything
funnels through — Theorem 6.5 — survives a verification that shares no
convention, no regularisation and no code with its derivation.

*Referee, 2026-08-02.*
