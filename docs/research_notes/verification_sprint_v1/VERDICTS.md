# Verification Sprint V1 — verdicts (2026-08-02)

**Scope: verify/correct/kill only; no new idea generation.** Workstreams:
A (orchestrator, exact calculus), B/C/D/E (independent agents), plus
orchestrator hand-verification of the load-bearing algebra of C
(\(\nu_*=\mathrm{Cov}/2V_r\), \(\mathfrak d_*=1-\mathrm{Cov}^2H_r/V_rG_r\),
\(r=0\) energy-neutrality vacuity), D (zero-diagonal lemma, edge ratio
conditions), and B (DSS \(L^3\) log-periodicity).

| candidate / claim | verdict | binding reason | deliverable |
|---|---|---|---|
| Steady continuum front profile (continuous self-similar, critical wake) | **KILLED** (twice) | Tsai 1998 (finite local energy) **and** Chae–Wolf ARMA 2017 Cor 1.4 at the weak-\(L^3\) level — the NRS log-escape buys nothing | [dss notes](dss_admissibility_notes.md) |
| \(\lambda=2\) DSS front orbit (ℝ³) | **CONDITIONAL** | survives on exactly one hypothesis: profile in \(L^{3,\infty}\setminus L^3\) (all candidates are uniformly weak-\(L^3\) bounded — exact lacunary computation); UNKNOWN: Chae–Wolf CPDE Thm 1.3 \(\lambda_+\) constant; **no log-growing \(L^3\) branch exists on ℝ³** (exact DSS log-periodicity ⇒ \(L^3\) either \(\equiv\infty\) or bounded; bounded is killed for every \(\lambda>1\)) — the torus log-growth is an IR-cutoff artefact and must be stated separately | [matrix](dss_admissibility_matrix.csv) |
| Multi-type log-periodic cycle, super-period ratio 2 | **CONDITIONAL** (= same object) | it *is* 2-DSS; identical row to the above; subdividing the composite ratio toward 1 moves it **into** the known kill zone — design rule | [dss notes](dss_admissibility_notes.md) §9 |
| Zeno packet relay / torus phase-coded cloud | **CONDITIONAL** (unprotected) | inadmissible to the entire SS/DSS literature (torus: no dilation group) — neither killed nor shielded; all exclusion pressure is ESS/Serrin/CKN + repo gates, which must be held to a higher standard | same |
| Front-flow rate statements | **CORRECTED** | one-sided/two-sided/limit regimes separated; coefficient asymptotics only for \(a\to a_\infty\); exact \(S\)-periodicity of \(N^2(T-t)\) for periodic \(a\); \(c_E\)-lever restated as high-Re limit with the inviscid-uniformity obligation | [corrected theorem](corrected_front_flow_theorem.md) |
| Λ monotone + dichotomy (identities (I.1)–(I.4), Lemma K, (M), D) | **PASS** | exact telescoping verified in `Fraction` on triads/P-families, all orders/viscosities; \(\varepsilon\)-regularised identity holds with factor \(1+\varepsilon H_0/H_1\) | [O-9 note](lambda_O9_defect_decomposition.md) |
| Λ quantitative sharpening (old O-8 uniform deficit floor) | **KILLED** | exact relay family \(\mathfrak d_*^{\rm full}=15D^2/(2+18D^2)\to0\), certified points to \(15/524306\); pilot deficits 0.865/0.740 were \(\nu\)-artefacts (true triad values \(3/4\), \(6/11\)) | same, §3 |
| Λ O-9 via energy neutrality / leakage | **KILLED** (vacuous at r=0) | saturation \(a_k=\gamma(x_k-\mu)e_k\) satisfies \(T_0=0\) automatically; leakage drivable to \(4\times10^{-10}\) | same, §5 |
| Λ O-8′ (deficit \(\ge c\,V_0/\mu^2\), single-shell collapse) | **CONDITIONAL** (new, measured not proved) | relay-family limits give ratios \(5/4\) (full), \(1/2\) (in-support); smallest observed 0.126 from an under-converged constrained optimiser — recorded as unresolved; next step named (projected-gradient rerun) | same, §6 |
| Parametric Prop 1 as stated ("forest ⇒ decay") | **KILLED as stated / PROVEN as corrected** | forest alone is neither necessary nor sufficient: per-direction modulation on a forest grows (+534 margin control); corrected hypotheses (H1) reciprocity + (H2) edge ellipticity \(\beta_{ij}\beta_{ji}<0\) + (H3) cycle balance \(\prod\rho_e=1\) ⇒ proven for **all three modulation classes at once**, Lean-able; balanced cycle (\(\mathcal C=-1\)) as dead as a tree (margin \(+4.6\times10^{-13}\)); 4.32M random + 3.60M exhaustive points, max margin \(+6.8\times10^{-10}\) = round-off | [audit](parametric_common_lyapunov_audit.md) |
| Parametric cycle escape (\(\mathcal C\neq-1\)) | **PASS** (sharpened) | \(\mathcal C=-1239/128\) exact and gauge-invariant; Prop 2 sharpened to \(\sigma_{\rm cycle}\approx|\det K|/c_1\) (15 digits); corrections: spurious \(\sqrt2\) in §B.6, ≈92 % of the staggered gain is Trotter averaging not holonomy, \(\rho(M)\ge1\) is structural (\(\det M=1\)); **(H1) remains undischarged from the NS convolution** — overall lane CONDITIONAL | same |
| Continuum-to-lattice shadowing proof route | **KILLED** (in scanned box) | no orbit exists (pilot A); one-stage map is **expanding** (\(L=14.6\)–\(15.0\) at \(c_E=1\), \(1.73\)–\(1.81\) at \(c_E=100\)); map-level consistency error is \(O(1)\) (0.65–0.84); pre-registered rule (contraction margin > consistency error) fails by over an order of magnitude; convolution-level consistency \(\sigma=1.825\) (\(R^2=0.999\), diagnostic) with a correctly-firing null control | [shadowing notes](poincare_shadowing_notes.md) |

## Net portfolio state after V1

- **Regularity side**: Λ survives as *identity + dichotomy* (demoted from
  quantitative sharpening); the single live quantitative question is O-8′.
- **Singularity side**: the only lanes not killed are (i) the torus
  phase-coded cloud (outside all SS/DSS literature, unprotected), and
  (ii) an ℝ³ \(\lambda=2\) DSS orbit **with a weak-\(L^3\)-only profile** —
  and the \(L^3\)-vs-\(L^{3,\infty}\) hypothesis is now the *entire* margin
  of that lane. The steady profile and the shadowing-based orbit route are
  dead.
- **Mechanism tools**: corrected Prop 1 ((H1)(H2)(H3)) is a proven,
  Lean-able no-go; the cycle escape survives conditionally on discharging
  (H1) from the actual NS convolution.

**UNKNOWN registry** (must be resolved before any promotion): U1
Chae–Wolf CPDE Thm 1.3 \(\lambda_+\); U2 DSS/ADSS Liouville at
\(L^{3,\infty}\) profiles; U3 ESS transfer to \(\mathbb T^3\); U4 weak-\(L^3\)
smallness constant; U5 forced-NS DSS theorems (none found); O-8′ constant;
(H1) discharge.
