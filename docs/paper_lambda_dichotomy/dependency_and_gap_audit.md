# Dependency and gap audit

For [theorem_statement.md](theorem_statement.md) and
[complete_proof.md](complete_proof.md).

## 1. Dependency graph (lemmas in proof order)

```
F1 (external) ──┐
                ├─→ L0a ─→ L1 ─→ L2 ─┐
                │         │          ├─→ L4 (Main a) ─→ L5 (Main b) ←─ F2 (external)
                │         │   L3 ────┘        │
                │         │                   └─→ L8 (Main e)
                │         └─→ L6 (Main c) + Cor 6′
                │             L7 (Main d)  [uses F3]
                │
L1,L4,L7 ──→ L0b (positivity of H_0; consumed by L5, L6, L8 and by the
              pointwise use of z on all of [0,T_max))

L9 (band symmetry, general radial weight) ─→ L10 (smooth family laws) ─┐
L11 (lattice bounds; N_0^2 ≍ N, bounded 2^20-adic gaps) ───────────────┤
                                                                       ├─→ Theorem O
CAPACITY THEOREM  ‖P(u_N·∇u_N)‖² ≥ c_0 N³  ─────────────────────────────┘  (UNCONDITIONAL
  = lstar/lstar_proof_main.md Thm 7.1(3), proven, cited not reproduced      static no-go)
  │
  ├─ Thm 1.3 (exact laws)      ← L9′ = L9 (same lemma, general weight)
  ├─ Prop 2.4 (N_0^2 ≍ N)      ← elementary shell counts
  ├─ Thm 3.2/3.3 (profile V, far field V = π²W + O(|y|^{-∞}))  ← E1
  ├─ Lem 4.1 + Thm 4.2 (exact lattice-Riemann identity, O(N^{-1} log N))
  ├─ Lem 5.2 + Thm 5.3 + Cor 5.4 (concentrated test field ψ_N, duality)
  └─ Thm 6.5 + Lem 6.6 + Thm 6.7 ((V-NONDEG): P(V·∇V) ≢ 0)   ← E1, E2, E3
```

**No conditional edge remains.** The former dashed edge `L* ──→ Proposition`
is replaced by a solid edge from a proven theorem. `L9` and `L10` were
restated (general radial weight; smooth truncation); `L11`'s constants were
replaced by the tighter, `N`-uniform ones of Prop 2.4; `L1`–`L8` and the Main
Theorem (a)–(e) are **untouched** — nothing in them ever mentioned the family.

**Referee-mandated correction (report A, M1; report B, critical 2).**
The original positivity argument for \(H_0\) was circular; it is
replaced by the bootstrap Lemma 0b, which uses L1, L4, L7 *on the
interval where positivity holds by definition of the first vanishing
time* and derives a contradiction. The graph is therefore not a
straight line: L0b sits after L7. All of L1–L4, L6, L7 are pointwise
statements needing only \(H_0,H_1>0\) at the time considered, so no
circularity remains.

## 2. External classical inputs (framework, not novelty)

| ID | statement used | where | status |
|---|---|---|---|
| F1 | local existence, uniqueness (forward), and the maximal strong solution \(u\in C([0,T_{\max});H^m)\cap C^1([0,T_{\max});H^{m-2})\), \(m>5/2\), for \(u_0\in H^m_\sigma\) on \(\mathbb T^3\) | L0 (framework; forward uniqueness for \(H_0>0\)) | classical (Kato-type); recorded as an audited external input in the repository's EXT discipline, never axiomatised |
| F2 | subcritical Serrin-class regularity: a strong solution in \(L^\infty(0,T;L^6)\) with \(2/\infty+3/6<1\) is regular on \((0,T]\) and extends | L5 | classical (Prodi–Serrin/Ladyzhenskaya family) |
| F3 | mean-zero Sobolev embedding \(H^1(\mathbb T^3)\hookrightarrow L^6\), constant \(C_S\) | L7 | classical; explicit constants available |

The Main Theorem consumes no other external result. In particular it
does not use: ESS endpoint \(L^\infty L^3\), BKM, CKN, backward
uniqueness, or any self-similar exclusion.

**Additional external inputs consumed by the Capacity Theorem** (§3),
and by nothing else in this paper:

| ID | statement used | where in `lstar/lstar_proof_main.md` | status |
|---|---|---|---|
| E1 | Fourier transform of homogeneous distributions on \(\mathbb R^3\): \(\mathcal F[\vert y\vert^{-\alpha}]=c_{3,\alpha}\vert\zeta\vert^{\alpha-3}\), \(c_{3,\alpha}=2^{3-\alpha}\pi^{3/2}\Gamma(\frac{3-\alpha}2)/\Gamma(\frac\alpha2)\), by analytic continuation / finite part, at \(\alpha=-1,1,2,4,6\) (no \(\Gamma\)-poles there) | Thm 3.3(a) (Oseen tensor), Thm 6.5 (transforms of \(r^{-2},r^{-4},r^{-6}\)) | classical (Gelfand–Shilov / Stein–Weiss). Referee: **accepted**, and re-derived by hand; flagged as "PROVEN modulo standard" in substance (MINOR-5(b)), with an independent numerical closure — see §4.8 |
| E2 | convolution theorem \(\widetilde A\widetilde B=\widetilde{A*B}\) applied to \(h\notin L^1\) (only \(h\in L^p\), \(p<3/2\), plus compact-support pieces) | Lem 6.1, Thm 6.5 | standard after splitting \(h=h\mathbb 1_{\vert\xi\vert\le1}+h\mathbb 1_{\vert\xi\vert>1}\). Referee: **accepted**; flagged "PROVEN modulo standard" (MINOR-5(a)) and **empirically closed** by an independent quadrature — see §4.8 |
| E3 | density of \(C^\infty_{c,\sigma}(\mathbb R^3)\) in \(L^2_\sigma(\mathbb R^3)\); Helmholtz–Leray decomposition on \(\mathbb R^3\) | Lem 6.2 | classical. Referee: **accepted** — but this is a pure *existence* argument and is the sole source of the non-effectivity of \(c_0,N_*\) (§5) |
| E4 | Plancherel and Young's inequality on \(\mathbb R^3\) | Lem 6.1, Cor 5.4 | classical, routine |

*(Paley–Wiener is quoted in Thm 3.2(3) for the entirety/exponential type
of \(V\); it is decorative and consumed nowhere.)*

## 3. The conditional bridge — removed

**There is no conditional bridge.** Both the Main Theorem (a)–(e) and the
static no-go (Theorem O) are unconditional modulo F1–F3 (the no-go itself
uses none of F1–F3 — it is a statement about trigonometric fields, not
about solutions).

**What replaced it.** The former **Hypothesis (L\*)** — capacity lower
bound \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) for the
*sharply* truncated coherent family — is replaced by the **Capacity
Theorem**: the same bound, with the same sharp exponent \(3\), **proven**
for the *smoothly* truncated family
\(\widehat u_N(k)=\chi(|k|/N)P_kv_0/|k|^2\), uniformly in the admissible
cutoff \(\chi\) and in \(v_0\in\mathbb R^3\setminus\{0\}\). This is
Theorem 7.1(3) of
[lstar/lstar_proof_main.md](lstar/lstar_proof_main.md).

**Why that suffices.** Theorem O's hypothesis quantifies over *all* real
zero-mean divergence-free trigonometric fields. Each \(u_N\) of the smooth
family is such a field (finite band \(|k|\le N\); real; zero mean; exactly
divergence-free). Exhibiting the capacity bound for *one* admissible
family is therefore all the argument ever required — the sharp family was
never load-bearing, only historically first.

**Honest headline (referee-mandated).** The correct statement is **"the
static no-go is unconditional"**, ***not*** "(L\*) is proven". The paper's
literal Hypothesis (L\*), for the sharply truncated family, is untouched
and still **open**; what is proven is the capacity bound for a
*different, equally admissible* family. These must not be blurred.

**New inputs the capacity proof consumes.** Beyond E1–E4 above, the proof
introduces no external results; the substantive content is its own
(Theorems 1.3, 3.2, 3.3, 4.2, 5.3, 6.5, 6.7; Prop 2.4; Lemmas 2.1–2.3,
4.1, 5.2, 6.1, 6.2, 6.6). Two steps were flagged by that document's own
author as **"PROVEN modulo standard"**; both were examined by the
adversarial referee and **both were accepted**:

1. **The symbol estimate \(V=V_\infty+O(|y|^{-p})\)** — i.e.
   \((1-\chi)h\) is a \(C^p\) symbol of order \(-2\) vanishing near the
   origin, so its transform decays like \(|y|^{-p}\) away from the
   origin. **Accepted.** In the main document this is Theorem 3.3(b),
   where it is not a gesture but a complete proof: \(|\partial^\beta G|
   \lesssim\langle\xi\rangle^{-2-|\beta|}\) is in \(L^1\) for
   \(|\beta|\ge2\), and \(y^\beta\widetilde G=i^{|\beta|}
   \widetilde{\partial^\beta G}\) is then bounded by
   \(\|\partial^\beta G\|_{L^1}\). Referee verdict: *"correct
   (\(|\beta|\ge2\Rightarrow\partial^\beta G\in L^1\))"*.
2. **Riemann integrability of \(C(\zeta)\overline{\widetilde\Psi(\zeta)}\)**
   — the Fourier-side route to the duality limit. **Accepted, and
   moreover not consumed:** the main document does *not* take that route.
   It uses instead the exact lattice-Riemann identity of Lemma 4.1 with
   the explicit defect bound of Theorem 4.2, which is quantitative rather
   than a qualitative integrability appeal, and which additionally
   delivers the derivative statement. The referee verified Lemma 4.1 as
   correct and Theorem 4.2 as structurally correct (with four constant
   corrections in MINOR-3 that *do not change* the \(O(N^{-1}\log N)\)
   rate), and recorded separately (MINOR-6) that the Fourier-side route
   of the numerics note is unused, "so nothing propagates".

The referee raised two *further* "PROVEN modulo standard" flags of its
own (E1, E2 above, its MINOR-5); it accepted both and closed them
empirically — see §4.8.

## 3a. Proven negative results about the *sharp* family (recorded, not used)

From
[lstar/lstar_direct_route_and_weakening.md](lstar/lstar_direct_route_and_weakening.md):

* **JOB A — no-go, PROVEN (Theorem A there).** For the sharply truncated
  family, *every* estimate of the form "pair against \((c\cdot\nabla)u_N\)
  for a constant \(c\), split \(u_N=b+w\) with \(b\) constant, bound the
  remainder by Hölder against \(\|\nabla u_N\|_{L^q}\)" is capped at
  \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge C_1S_N\asymp N\) — **two full
  powers short** of the target \(N^3\), and short even of \(N^2\). The
  proof is elementary: \(\sup_c b^\top Mc/\sqrt{c^\top Mc}
  =\sqrt{b^\top Mb}\), \(\lambda_{\max}(M)\le\tfrac12H_1\), and
  \(\|u_N-b\|_\infty\ge|b|-C_0\) with \(C_0=\|v_0\|\sqrt{\tfrac23T_\infty}\)
  force \(|b|<(2+\sqrt2)C_0\), while the paper's natural sweeping constant
  \(|u_N(0)|=\tfrac23S_N\|v_0\|\to\infty\) is excluded outright. Referee:
  arithmetic **re-derived and confirmed**. *Scope, stated honestly:* it
  closes constant-vector sweeping splits only; it does **not** close
  pairing against a general divergence-free test field — which is exactly
  the route the Capacity Theorem takes, so Theorem A is in fact a positive
  argument for that route.
* **JOB B — weakening, PROVEN but now MOOT.** \(N^3\) can be relaxed to
  \(N^2g(N)\) for any nondecreasing \(g\) with
  \(\int^\infty dt/(t\,g(t))<\infty\) — e.g. \(N^2(\log N)^{1+\varepsilon}\)
  — and the no-go still follows; the threshold is sharp *for that
  argument*. Referee: exponent bookkeeping **correct end to end**. It is
  **moot** because the Capacity Theorem delivers the full \(N^3\); the
  hypotheses (L\*-weak) and (L\*-min) are therefore not used and are
  demoted to a remark. One by-product of §B is retained and *is* used:
  \(T_N\ge T_1=6\), which makes the dyadic-gap bound of L11 elementary and
  \(N\)-uniform, replacing the earlier asymptotic "\(T_{2N}/T_N\to1\)".

## 4. Adversarial re-verification performed

1. **Convention traps.** The repository's paired ledger produces a
   factor-2 slip in the modal Cauchy–Schwarz (caught and eliminated: the
   paper uses full-lattice unpaired sums, where
   \(|\alpha_k|\le\sqrt{\varepsilon_k\eta_k}\) is exact); an earlier
   factor-2 slip between the exact and float \(K\)-pipelines was caught
   by the \(N=4\) exact anchor and fixed before certification (the two
   pipelines now agree to \(0.0\) at \(N=4\)).
2. **Degenerate cases.** \(u_0=0\) (excluded; trivial solution);
   \(H_0=0\) or \(H_1=0\) at some time (proven impossible for
   \(u_0\neq0\) — L0, via the Grönwall two-sidedness of \(H_0\) and
   forward uniqueness); \(D=0\) divisions (none occur: L6 divides by
   \(D>0\); L8 is stated in product form).
3. **Limit operations.** No spectral-cutoff or Galerkin limit is used
   anywhere in the Main Theorem: all identities are proven directly on
   the strong solution via \(L^2\)-pairings valid for \(m>5/2\) (L0);
   term-by-term differentiation is justified by difference quotients in
   \(C^1H^{m-2}\), not by formal exchange of sum and derivative.
4. **PDE embedding of the finite algebra.** The finite-Fourier
   identities (I.1)–(I.4) hold verbatim for the solution's (infinite)
   lattice sums because every sum involved is absolutely convergent at
   the stated regularity and the inequalities are termwise; the
   repository's exact `Fraction` telescoping certificates verify the
   algebra on finite fields independently.
5. **Counterexample cross-checks.** The coherent family was checked
   against: single-shell controls (\(K\asymp N^{-2}\), matching the
   divergence-free depletion prediction), Beltrami controls
   (\(\mathcal N\equiv0\) for exact Beltrami; perturbed value bounded),
   and the mesoscopic narrowband parent (\(K\approx5\times10^{-3}\)) —
   the family is the maximiser class among all tested classes, and its
   growth law has an exact anchor.
6. **Cross-talk with the singularity lane.** The Proposition does not
   assert that large-\(K\) fields occur along solutions; the separate
   persistence measurement (Osgood-gate note §5) is quoted only as a
   diagnostic and is **not** part of any proof.
7. **Novelty check, per statement.** (a)+(b): we are not aware of the
   bandwidth monotone \(\Lambda\) or the dichotomy in this exact form;
   both are elementary, so overlap with known Serrin-refinement
   literature is possible — the honest claim is the *packaging* (exact
   modal identities, machine-verifiable defects, and the action
   \(\int KD\) dominated by both classical critical actions).
   (c): the identity is an elementary expansion of
   \(\|\partial_tu+\nu Au\|^2\); its use to localise the criterion at
   the \(\dot H^1\)-bandwidth action appears new to us but is close in
   spirit to known \(\|\partial_tu\|\)- and \(\|Au\|\)-based criteria.
   (d): routine Hölder/Sobolev. Theorem O: the coherent
   critical-spectrum family and its exact \((2/3)\)-symmetry laws appear
   new as a *proven* Osgood obstruction; Bernstein-extremal spectra
   are classical folklore, and the profile \(V\) is the (classical) Oseen
   tensor applied to a smoothed Newtonian potential — the novelty claimed
   is the combination (exact scaling \(\widehat u_N(k)=N^{-2}F(k/N)\),
   lattice-Riemann inner rescaling, concentrated-test-field duality, and
   the \(\chi\)-uniform non-degeneracy), not any ingredient. No claim is
   made beyond this.
8. **Adversarial referee on the capacity proof**
   ([lstar/lstar_referee.md](lstar/lstar_referee.md), verdict **(i)**).
   Brief: *break the main route*. Result: **not broken at any of the six
   points in the brief.** Verified affirmatively: the Fourier conventions
   and every \((2\pi)^{-3}\) and power of \(N\); Lemma 9′/9″; the exact
   laws (Thm 1.3); every arithmetic step of §2 (shell count
   \((2^{j+2}-1)^3<2^{3j+6}\), \(S_N\le128N\), \(0.226M\) box count,
   \(N_0^2\in[N/44544,\tfrac{64}3N]\)); Thm 3.2(1)–(5); the Oseen
   constants \(c_{3,1}=4\pi\), \(c_{3,2}=2\pi^2\), \(c_{3,-1}=-8\pi\) and
   \(\widetilde h=\pi^2W\); Thm 3.3(b); Lemma 4.1 and the structure of
   Thm 4.2; Lemma 5.2 (div-free, zero mean, \(N\)-independent norm) and
   the Thm 5.3 / Cor 5.4 bookkeeping; Lemma 6.1; **Theorem 6.5 confirmed
   to 14 digits** by an independent two-patch spherical quadrature of
   \(\tau=h*h\) that shares no code, no Fourier convention and no
   regularisation with the derivation (rel. err.
   \(3.6\!\times\!10^{-14}\), \(1.2\!\times\!10^{-14}\),
   \(6.6\!\times\!10^{-15}\)) — this closes the finite-part risk flagged
   as E1/E2; Lemma 6.5′ re-derived by hand; Lemma 6.6 and the threshold
   \(\pi^2/2048\); Thm 7.2's arithmetic (\(288^2=82944\),
   \(2^{20}>950272\), \(c_-=0.0985\), \(c_+=27.63\)); and Theorem A of the
   direct-route note. An independent lattice evaluation of
   \(u_N(y/N)/N\) found the true rate \(E_N=\Theta(1/N)\), i.e. Thm 4.2 is
   valid but **not tight** (no logarithm). The two numerics documents,
   which disagree at face value by \(\sim\!250\times\), were shown to
   **reconcile** once the seed is normalised (main doc §7.5 used a unit
   \(v_0\)): \(K/N_0^2\approx0.468\) vs \(0.454\) at \(N=40\), a 3 %
   difference attributable to the different \(\chi\) — counted as a real
   cross-validation. Defects found: one **false** proposition
   (the Abel-regularised Poisson identity — see §5, now deleted; the
   referee proved that no summation method can rescue it, because the
   divergent terms are eventually positive), the paper-level integration
   burden (discharged by the present revision), the "all constants are
   explicit" claim (corrected; see §5), and cosmetics. The false
   statement was **inert**: §§4–7 never cited it.

## 5. Gaps and scope limits (complete list, post-referee)

| item | consumed by | status |
|---|---|---|
| capacity bound for the **smooth** family, \(\ge c_0N^3\) | Theorem O | **PROVEN** — Thm 7.1(3) of the capacity document; sharp exponent; uniform in \(\chi\) and \(v_0\). Referee could not break it and confirmed the pivotal identity (Thm 6.5) to 14 digits by an independent method |
| effectivity of \(c_0\) and \(N_*\) | Theorem O's constant \(c\) | **NON-EFFECTIVE**, and this is now stated in the paper. The test field \(\Psi\) comes from density of \(C^\infty_{c,\sigma}\) in \(L^2_\sigma\) (E3), a pure existence argument: neither \(\Psi\), nor its radius \(R\), nor \(I_\Psi\), nor \(c_0\), nor \(N_*\) is produced. Not a gap in the proof; it *was* a false advertisement in the old "all constants are explicit" claim, now corrected. An effective version needs an explicit \(\Psi\) plus a quantitative lower bound on \(I_\Psi\); the certified non-vanishing region \(\{|\zeta|<\pi^2/2048\}\) forces \(R\gtrsim10^3\), \(N_*\gtrsim10^4\), and a \(c_0\) far below the measured \(\approx8.4\|v_0\|^4\) |
| (L\*) for the **sharply** truncated family | — | **OPEN, and no longer used.** Untouched by the capacity proof: the smoothness of \(\chi\) is needed at exactly one point (\(|\nabla F|\lesssim|\xi|^{-3}\) in Thm 4.2, which fails for \(\chi=\mathbb 1_{[0,1]}\); \(V\) then acquires an oscillatory \(|y|^{-2}\) tail and the Riemann defect a Gauss-circle term). §§5–6 are insensitive to \(\chi\), so the gap is technical, not structural. Certified finitely (N ≤ 8 exact, N ≤ 32 float). Additionally constrained by the JOB A no-go (§3a) |
| Abel-regularised Poisson identity (former Prop 3.4 / eq. (3.2) of the capacity document) | — | **FALSE and DELETED** (referee MAJOR-1). \(V\) decays only like \(|y|^{-1}\) with a sign-definite spherical mean, so the periodisation diverges; the terms being eventually positive, *no* summation method — Abel included — rescues it, and the offered proof was independently unsound twice over. **Impact: none** — the statement was inert, cited nowhere in §§4–7; §4's exact lattice-Riemann identity does the work. Replaced by a remark recording the no-go |
| Riemann-sum rate \(O(N^{-1}\log N)\) (Thm 4.2) | capacity proof | **PROVEN but not tight**: independent lattice evaluation gives \(E_N=\Theta(1/N)\), no logarithm. Four constant corrections (referee MINOR-3) do not change the rate, and the rate is used only qualitatively (\(\mathcal E_N\to0\)) |
| exhaustiveness of a two-case blow-up/global split | — | **not claimed** (referee A M2): (b) is stated as the proven pair {action finite ⇒ global} / {action infinite}, with \(T_{\max}<\infty\Rightarrow\) action infinite; a global solution with infinite action is not excluded |
| equivalence of \(\int KD<\infty\) with the bandwidth action | — | **false and withdrawn** (referee A C1 / B critical 2): only the one-sided comparison \(\int KD\le2\int(\|\partial_tu\|^2/D+\nu^2N_1^2)\) holds (Cor 6′); explicit decaying counterexample recorded |
| equality-case characterisation in (I.4) | — | **not asserted** in this paper (removed) |
| scope of Theorem O | clause (e) | excludes only the \(R\equiv0\) field-inequality route; solution-adapted \(R\) is out of scope here |
| sharp limit \(N^{-3}\|\mathbb P(u_N\!\cdot\!\nabla u_N)\|^2\to(2\pi)^{-3}\|\mathbb P(V\!\cdot\!\nabla V)\|^2_{L^2(\mathbb R^3)}\) | — | **conjectural, flagged as such, used nowhere.** It is the Cauchy–Schwarz-optimal form of Cor 5.4 — i.e. "the localised test field is asymptotically extremal" — which is why \(a=1\) is sharp and cannot be improved |
| convention note | certificates | `summary.json` stores half-energies, `exact_family_certificates.json` full; the ratio \(K\) is convention-invariant and all quoted values were re-verified independently by referee A (one correction: captured-fraction range is 91.18%–92.10%) |

With these corrections the Main Theorem's proof chain
L0a→L1→L2→L3→L4→{L7→L0b}→L5(+F2), L6+Cor 6′, L7(+F3), L8 is complete
at the stated regularity, and Theorem O's chain
L9→L10, L11, Capacity Theorem → Theorem O is complete and
**unconditional**. Referee reports:
[referee_report_A.md](referee_report_A.md),
[referee_report_B.md](referee_report_B.md) (both MAJOR-REVISIONS; all
CRITICAL and MAJOR items are addressed by the current text — the
re-derivation referee confirmed items (i)–(iv) reproduce exactly with
no sign or constant errors), and
[lstar/lstar_referee.md](lstar/lstar_referee.md) on the capacity proof
(verdict (i); no CRITICAL; three MAJOR items, all discharged by the
present revision: the false proposition deleted, the family switch
integrated, the effectivity claim corrected).

Per the freeze directive, the former "Clay proof candidate" designation
is removed; the correct description of the result is: **an unconditional
bandwidth–dissipation dichotomy with an exact action representation and a
one-sided bandwidth-action criterion, plus an unconditional no-go showing
the dichotomy's hypothesis cannot be verified by any uniform
Osgood-pointwise field inequality.** The no-go's constants \(c_0,N_*\)
are non-effective. The literal Hypothesis (L\*), for the sharply
truncated family, remains open and is used by nothing.
