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
L9 ─→ L10 ─┐
           ├─→ Proposition (conditional no-go)  [uses L*]
L11 ───────┘
```

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

No other external result is consumed. In particular the Main Theorem
does not use: ESS endpoint \(L^\infty L^3\), BKM, CKN, backward
uniqueness, or any self-similar exclusion.

## 3. The single conditional bridge

**Hypothesis (L\*)** — capacity lower bound
\(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) for the coherent
family — is consumed **only** by the Proposition (static no-go). The
Main Theorem (a)–(e) is unconditional (modulo F1–F3). Per the freeze
directive, since this bridge is unproven the no-go is stated as a
**conditional proposition**, not folded into the theorem.

Support for L\* (certificates, no proof): exact rational
\(K(u_N)\) at \(N=4,6,8\) with float agreement \(\le1.3\times10^{-15}\);
the proven sweeping-projection identity
\(\mathbb P((c\cdot\nabla)u)=(c\cdot\nabla)u\) and the resulting exact
pairing lower bound capturing 91.2–92.1% of
\(\|\mathbb P(u\cdot\nabla u)\|_2\), essentially \(N\)-independent for
\(N\le32\); \(K/N_0^2\nearrow0.396\). Missing: a rigorous asymptotic
lower bound on the explicit lattice triple sum
\(Q(v_0)=\langle(u_N\cdot\nabla)u_N,(v_0\cdot\nabla)u_N\rangle\)
(certified positive and \(\asymp N^2\)-sized numerically at all tested
\(N\)); with the proven upper bound
\(\|(v_0\cdot\nabla)u_N\|_2^2\le\tfrac23\|v_0\|^4S_N\), the statement
\(Q(v_0)\ge cN^2\) would imply L\*. This is the complete reduction; it
is one explicit lattice-sum estimate.

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
   (d): routine Hölder/Sobolev. Proposition: the coherent
   critical-spectrum family and its exact \((2/3)\)-symmetry laws appear
   new as a *certified* Osgood obstruction; Bernstein-extremal spectra
   are classical folklore. No claim is made beyond this.

## 5. Gaps and scope limits (complete list, post-referee)

| item | consumed by | status |
|---|---|---|
| L\* (lattice triple-sum lower bound) | Proposition only | open; certified finitely (N ≤ 8 exact, N ≤ 32 float); the Main Theorem is unaffected |
| exhaustiveness of a two-case blow-up/global split | — | **not claimed** (referee A M2): (b) is stated as the proven pair {action finite ⇒ global} / {action infinite}, with \(T_{\max}<\infty\Rightarrow\) action infinite; a global solution with infinite action is not excluded |
| equivalence of \(\int KD<\infty\) with the bandwidth action | — | **false and withdrawn** (referee A C1 / B critical 2): only the one-sided comparison \(\int KD\le2\int(\|\partial_tu\|^2/D+\nu^2N_1^2)\) holds (Cor 6′); explicit decaying counterexample recorded |
| equality-case characterisation in (I.4) | — | **not asserted** in this paper (removed) |
| scope of the Proposition | clause (e) | excludes only the \(R\equiv0\) field-inequality route; solution-adapted \(R\) is out of scope here |
| convention note | certificates | `summary.json` stores half-energies, `exact_family_certificates.json` full; the ratio \(K\) is convention-invariant and all quoted values were re-verified independently by referee A (one correction: captured-fraction range is 91.18%–92.10%) |

With these corrections the Main Theorem's proof chain
L0a→L1→L2→L3→L4→{L7→L0b}→L5(+F2), L6+Cor 6′, L7(+F3), L8 is complete
at the stated regularity. Referee reports:
[referee_report_A.md](referee_report_A.md),
[referee_report_B.md](referee_report_B.md) (both MAJOR-REVISIONS; all
CRITICAL and MAJOR items are addressed by the current text — the
re-derivation referee confirmed items (i)–(iv) reproduce exactly with
no sign or constant errors). Per the freeze directive, the former
"Clay proof candidate" designation is removed; the correct description
of the result is: **an unconditional bandwidth–dissipation dichotomy
with an exact action representation and a one-sided bandwidth-action
criterion, plus a finitely-certified conditional no-go showing the
dichotomy's hypothesis cannot be verified by any uniform
Osgood-pointwise field inequality.**
