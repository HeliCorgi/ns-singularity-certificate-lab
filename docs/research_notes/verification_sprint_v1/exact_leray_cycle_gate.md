# VERIFICATION SPRINT v1 — FREEZE PHASE, PART II
## Exact Leray cycle gate

> **REFEREE CORRECTIONS APPLIED (2026-08-02, per `cycle_gate_referee.md`).**
> The KILL verdict is unchanged; the referee confirmed the tree control, the
> §4/§5 statistics, single-shell rigidity, and the interval machinery reproduce
> cleanly. What changed in this revision:
>
> - **C1** — the §7.2 constant was wrong by `3^{3/2}`: corrected to
>   `C(c_E) = chi_1 sqrt(c_E)/(3 nu)` (`C(1) = 8.0826`, not `1.5555`); the
>   `c_E >= 9(1+q_0)^2 nu^2 N^3/chi_1^2` display and the zero-crossings
>   `N = 4.0/8.7/18.7/24.6` already followed from the correct constant and are
>   unchanged; the "agrees to four decimals" claim is restated against the
>   corrected constant, with the search-grid covariance caveat (M8).
> - **C2** — script/artifact mismatch repaired: `--squarings` exposed (default
>   11, i.e. `q` up to 2048, matching the archive), the script now writes
>   `float_reference` and the standalone certificate file, and a
>   `--stage certificate` mode re-runs only the certified candidate.
>   **Re-run 2026-08-02: the reproduced certificate is byte-identical to the
>   archived `cycle_gate_interval_certificate.json` (1.4 s).**
> - **C3** — the §6 "not a large-field artifact" threshold `a > 0.013975` was
>   numerically false (wrong-sign inequality); replaced by the measured
>   `a* = 0.03125`, `||U_P||_{L^2} = 0.131` (a factor ~30 below the ceiling).
> - **C4** — §3.1 restated against the archived `R1_by_kind = {73, 62, 65}`
>   (R1 does contain 4-cycles and double-triads); R2 justified on shell
>   coverage alone.
> - **C5** — the unarchived "6 535 748 geometries" figure (and its 3-way split)
>   deleted; the archived counts are quoted instead.
> - **M1** — the certified bound is now rounded **down**: `1.0163293277`.
> - **M2** — the float64 final root in the certificate is disclosed.
> - **M5** — "VIOLATED" relabelled `OUTSIDE HYPOTHESES ((H2) fails)`; §10.2
>   now says Prop 1 *does not apply*, not that its conclusion is false.
> - **M6** — the audit's registered (H1) follow-up is recorded as **open**,
>   not discharged.
> - **M8** — the "independent check" framing of the §7.2 table is withdrawn
>   (the searched schedule grid is scale-covariant by construction).
> - **M9** — the off-chain/main leakage trend is measurable from the archive
>   and is recorded: median `2.080` at `max|k|^2 <= 9` falling to `1.412` at
>   `max|k|^2 >= 25` — mildly decreasing in `N`, the **favourable** direction.
> - **M10** — step 6 is stated plainly as **never implemented** (not merely
>   unreached).
> - Minor fixes m1–m4, m8, m12–m14 (divergence-check field count, float-eig
>   labels, stage walls, file sizes, "as it must", exact-vs-float labels,
>   truncation caveat in the verdict box, Meshalkin–Sinai hedge).
> - **Open errata (would need re-runs; skipped and recorded):** M3 (exact
>   rational trace endpoints and `stage_C_refine` absent from the archive; a
>   fresh full run now stores both), M4 (the §7.1 512-orientation table has no
>   artifact; the `151/512` count is rng-state-dependent — an independent
>   re-run gives `155/512` at `nu = 1/40`), M7 (the implemented amplitude is
>   `sqrt(c_E/N)` without the stated `||C_a||` division, so absolute
>   zero-crossings in §7.2 are up to `4.5^{1/3} ~ 1.65x` optimistic).

> # VERDICT: **KILL**
>
> **Not killed at step 5 — the opposite.** Floquet growth exists and is
> **rigorously certified by interval arithmetic**: for the note's own B.5 cyclic
> carrier under its own staggered pump, in the full 4-dimensional
> divergence-free polarization space of every child mode (within the
> 12-dimensional truncation to the three cycle modes; §4 quantifies the
> discarded leakage), at `nu = 1/40` with `|a| <= 1`,
>
>     rho(M)  >=  1.0163293277   > 1   (CERTIFIED, 100-bit interval enclosure;
>                                       lower bound rounded DOWN)
>     single-parent instantaneous margins  =  -0.225, -0.225, -0.225  (< 0;
>                                       float eig — only rho(M) is interval-backed)
>
> so the note's "zero instantaneous margin" premise is **reproduced and
> certified, not refuted** — at a higher standard than the note itself reached.
> This discharges obligation 3 of the note's proof chain *affirmatively*.
> **The kill is not that nothing grows. It is that the growth means nothing.**
>
> **GROUND 1 (decisive, exact, complete) — the tree control fires.**
> A *single edge* `{k, k+p}` with one parent frozen at `|a| = 1` is a **forest**,
> on which Proposition 1 — *when its hypotheses (H1)–(H3) hold* — asserts no
> growth at any amplitude or waveform (forest topology makes (H3) vacuous). On
> the exact Leray blocks in full polarization, **3 750 of 28 480** such
> single-edge systems have `max Re lambda > 0` at `nu = 1/40`, up to
> `+3.429106`. The reason is checkable without any dynamics: hypothesis **(H2)
> fails exactly** — the `4 x 4` blocks `G_i B_ij` and `B_ji^T G_j` are *not
> proportional*, so no `rho_ij` exists at all, let alone a positive one. Prop 1
> is therefore **outside its hypotheses** and asserts nothing here: it simply
> does not apply to the exact Leray coupling. (The audit's registered top
> follow-up concerned **(H1)**, which this gate does *not* test; that
> obligation remains **open**. What is established here is narrower: (H2) fails
> for the exact Leray blocks in full polarization, so Prop 1 gives no control
> regardless of (H1).) Physically this is unsurprising: a single Fourier parent
> mode is a Kolmogorov flow, which is linearly unstable above a finite Reynolds
> number (Meshalkin–Sinai) — a heuristic attribution only; the exhibited
> counterexample is a 2-mode truncation outside that theorem's long-wave
> regime. The note's single-polarization reduction projected the unstable
> direction away.
>
> The idea's entire claim to novelty is
> *"impossible at any amplitude for every tree" -> "finite threshold for a cycle …
> that factor-infinity-to-finite jump is the whole idea."*
> **There is no such jump.** Both trees and cycles have finite thresholds; the
> cycle buys at most a constant factor.
>
> **GROUND 2 (exact, but only for the fixed-cardinality branch) — no uniformity
> in `N`.** Under critical normalization the pre-registered stage gain obeys
>
>     q_*(N)  =  chi_net / (nu N^2)  =  C(c_E) * N^{-3/2}  -  1 ,
>     C(c_E)  =  chi_1 sqrt(c_E) / (3 nu)
>
> an *identity* in this model, not a fit: the B.5 children lie on one shell, so
> `-nu D = -nu N^2 I` is scalar and `chi_net = chi_coupling - nu N^2` exactly,
> while `chi_coupling` is linear in the amplitude and in `s` for the
> scale-covariant schedule family searched (§7.2). It is strictly decreasing in
> `N` for every `c_E` and crosses zero at finite `N` — `N ~ 24.6` even at the
> repository's own measured `c_E = 228` (an upper estimate: the amplitude
> normalisation caveat of §7.2 could pull this down to `N ~ 15`; the direction
> is unaffected). Holding `q_* >= q_0 > 0` forces
> `c_E >= 9(1+q_0)^2 nu^2 N^3 / chi_1^2 -> infinity`, which is `judge_rigor`'s
> objection derived rather than estimated. The measured exponent `N^{-3/2}`
> **coincides exactly** with the repo's own `fixed_cardinality_scaling`
> rejection.
>
> **HONEST GAP — step 6 (the fattened Galerkin cloud) was NEVER IMPLEMENTED.**
> No Galerkin-cloud code exists in the script; the pre-registered ~50-minute
> compute budget was spent before it was written. Ground 2 therefore kills the
> *fixed-cardinality* branch only; the fixed-relative-width cloud
> (`W = floor(N^gamma)`, `M_eff = theta N^3`) — obligation 4 and kill condition
> K2 of the note — is **not** settled by this gate. Ground 1 does not depend on
> it. The exact measurements bearing on it are in §4: off-cycle leakage carries
> a **median `1.814x`** the in-cycle response energy over all shells
> `|k|^2 <= 50` (B.5 gadget: `2.144486...`), and — measured as a trend from the
> stage-A archive — the ratio is mildly **decreasing** in `N` (median `2.080`
> at `max|k|^2 <= 9` falling to `1.412` at `max|k|^2 >= 25`), the direction the
> GO condition asks for. Favourable to the cloud branch, but not a substitute
> for running it.
>
> **Autonomy is NOT discharged and was never tested** (flagged as required): the
> parent amplitudes are *prescribed*, from the most generous admissible class.
> Even a GO would have left autonomy as a later, untouched hurdle.

**Script:** `experiments/run_exact_leray_cycle_gate.py`
**Data:** `outputs/verification_sprint_v1/cycle_gate/`
**Read first:** `parametric_common_lyapunov_audit.md` (corrected Prop 1),
`../ideas_2026_08_01/idea_parametric_resonance.md` (the audited model).
**No new ideas, and no new variants of `Lambda`, `K`, `N_r`, Osgood or deficit
quantities were generated.** This note verifies, computes and decides only.

---

## 1. Framework (implemented exactly as specified)

**Parent field — the most generous admissible class.**

    U_P(x,t) = sum_a a_{p_a}(t) * C_a(x),
    C_a(x)   = c_a cos(p_a . x)  or  c_a sin(p_a . x),
    c_a . p_a = 0,  c_a integer primitive,
    |a_{p_a}(t)| <= 1,  piecewise constant, arbitrary switching and signs.

If even this cannot produce a gain uniform in `N`, the kill is airtight. It
cannot; see §7.

**Child band.** Integer wavevectors `k_1..k_L` closed under `k -> k +/- p_a`
along the cycle edges. **Every child mode carries its full real
divergence-free subspace, dimension 4**:

    { v1 cos(k.x), v2 cos(k.x), v1 sin(k.x), v2 sin(k.x) },
    v1 = k x e (primitive),  v2 = k x v1 (primitive),
    v1 . k = v2 . k = v1 . v2 = 0  (exact integers).

*No polarization is chosen for the children and nothing is projected away by
hand.* This is the single largest departure from the note, which fixes one real
polarization per child and thereby reduces the edges to scalars. In the honest
treatment the edge coefficients are `4 x 4` **blocks**, and the whole Lyapunov
analysis is done on the blocks (§5). The mode Gram is exactly diagonal,
`G_i = diag(|v1|^2, |v2|^2, |v1|^2, |v2|^2)/2`, so no irrational orthonormal
basis is ever needed and every projection stays in `Fraction`.

**Linearized child system.**

    xdot = A(t) x,   A(t) = -nu D + sum_a a_{p_a}(t) K_a,
    D    = blockdiag(|k_i|^2 I_4),
    K_a  from  L_a w = -P[(C_a . grad) w + (w . grad) C_a],

with **both** the advection-of-`w` and the advection-of-`U` terms, **both**
`+p` and `-p` sidebands (automatic: the real cos/sin product-to-sum identities
of `fourier_torus` emit `k+p` and `k-p` from every product), and the exact
per-mode Leray projection `P_k = I - k(x)k/|k|^2`. Convention as required:
`N_k = -i P_k sum_{l+m=k} (m.u_l) u_m`.

**Arithmetic labels.**

| quantity | label |
|---|---|
| Gram matrices, all coupling blocks `K_a[i][j]`, all cross-talk coefficients, edge ratios `rho_e`, holonomy products, divergence check, reality check, energy check | **EXACT** (`fractions.Fraction`) |
| matrix exponentials, monodromy, spectral radii, Floquet exponents, margins, scaling scan | **FLOAT** (float64; `scipy` absent, hand-rolled scaling-and-squaring Taylor `expm`) |
| the step-4 Floquet certificate | **INTERVAL** (`snapshot_certificate.Interval`, exact rational endpoints, 100-bit outward rounding) |

---

## 2. Mandatory correctness checks (all EXACT, all before any dynamics)

| check | result |
|---|---|
| **Reality** of every real-basis matrix | **895 488** block entries (in-cycle + cross-talk, all geometries, all generators); **every one an exact `Fraction`**, hence real |
| **Divergence-freeness** of every polarization | **36** reference fields (12 child basis + 3 parents x 8 generators); `max \|div\|` coefficient sum = **`0` exactly** |
| **Full-system energy conservation** (parents *dynamical*) | 18-field basis; `max` over all unordered triples of `sum_{sigma in S_3} <(f_s1.grad)f_s2, f_s3>` = **`0` exactly** |
| **div-free antisymmetry** `<(f.grad)g,h> + <(f.grad)h,g> = 0` | `max` defect over all `18^3 = 5 832` ordered triples = **`0` exactly** |

**Independent cross-check of the tensor against a published exact result.**
Re-deriving the note's own §B.1 triad
(`p=(1,1,0)`, `q=(1,0,1)`, `c=(2,1,1)`, `n=(1,-1,-1)`) through this module's
`leray_linearization` gives, in `Fraction`,

    alpha_B = -1/2,   alpha_C = -1/2,   alpha_D = +1/3,

reproducing the note's boxed values exactly, together with its detailed
conservation identity `alpha_B + alpha_C + 3 alpha_D = 0` and the energy
identity `|f_B|^2 a_B + |f_C|^2 a_C + |f_D|^2 a_D = 0` — both **exactly zero**.
The tensor implementation is validated.

---

## 3. Enumeration registry (pre-registered; the cap **binds**, reported honestly)

Pre-registered *before* the run:

| item | value |
|---|---|
| parents | canonical `p`, `1 <= \|p\|^2 <= 6` |
| children | `\|k\|^2 <= 50`, all cycle modes **distinct** |
| classes | 3-cycles (closed parent triad), 4-cycles (4 signed steps summing to zero), double-triads (`p,q,-p,-q` parallelograms, two parents alternating) |
| list **R1** | first **200** distinct geometries by increasing `max\|k\|^2` |
| per-(kind, shell) bucket cap | 24 |
| parent orientations | 8 per parent = 4 polarizations (`v1`, `v2`, `v1+v2`, `v1-v2`) x {cos, sin} |
| orientations taken to dynamics | <= 64 per geometry, ranked by the exact imbalance score |

### 3.1 The cap binds, and R1's shell coverage is too narrow — recorded, not hidden

What the archived enumeration records (`cycle_gate.json -> enumeration`): the
capped pass examined **1 094 958 / 20 973 928 / 2 343 120** raw paths
(3-cycle / 4-cycle / double-triad — *raw paths*, counted before the
distinct/valid filters, so not distinct geometries), the per-`(kind, shell)`
bucket cap of 24 overflowed on **777 378 / 11 950 156 / 1 412 724** of them,
and **3 025** distinct geometries were retained, from which R1 keeps 200. (An
earlier version of this note quoted an uncapped total of "6 535 748 distinct
geometries" with a three-way split; that figure appears in no artifact and
cannot be reproduced from them, so it is withdrawn. The caps bind hard — that
much the overflow counts establish — but the uncapped population size is
simply not archived.)

The archived R1 is *not* degenerate in cycle kind: `R1_by_kind =
{3-cycle: 73, 4-cycle: 62, double-triad: 65}` (`R1_DEGENERATE = false`) — the
per-`(kind, shell)` bucket cap of 24, listed in the registry table above, is
exactly what prevents kind-degeneracy. What *is* wrong with R1 as a sample is
its **shell coverage**: because the pre-registered ordering is "by increasing
`max|k|^2`", every member of R1 sits at `max|k|^2 <= 4`, against the
registered child bound `|k|^2 <= 50` — and R1 does not contain the note's own
B.5 gadget (`|k|^2 = 9`). That is sufficient ground for the supplement, and it
is the true one.

**Pre-registered supplement R2**, declared *before any dynamics was run*:
round-robin over the strata (cycle kind) x (`max|k|^2` shell), taking
geometries in R1 order within each stratum, 200 total, with the note's B.5
geometry forced in (`R2` covers shells `|k|^2 = 1 .. 50`). The analysed list
is `R2 union R1`, **382 distinct geometries** covering all three cycle classes
and shells `|k|^2 = 1 .. 50`. Every conclusion below is reported over this
list; the R1-only sample is retained in the JSON but is not load-bearing
anywhere.

---

## 4. Exact cross-talk: the cycle is not even the main channel

For every geometry, every parent, and every one of the 8 generators, the script
dumps **every** coupling coefficient from a cycle mode into an off-cycle mode
reachable by one parent interaction, exactly.

**Location:** `outputs/verification_sprint_v1/cycle_gate/crosstalk_tables.json`
(**34 688 exact `4 x 4` blocks**, `19.9 MB`, keyed by geometry / parent /
generator / source child / target off-cycle mode, every entry a `Fraction`
string).

The aggregate, over the 375 geometries for which both channels are nonzero:

    off-cycle (leakage) response energy  /  in-cycle response energy
       min 0.304    median 1.814    max 5.870
       B.5 gadget:  2.144486...  (exact)

**On the median geometry, and on the note's own gadget, roughly twice as much of
the exact Leray response leaves the cycle as stays in it.** Every number in
§5-§7 is computed on the truncation that discards this leakage entirely, so the
truncation is *generous* to the mechanism — and the kill survives it. It also
means the "diagonal cross-talk gate" the note claims to close exactly
(its §E.5: `X_{p_a}(t) X_{p_b}(t) = 0` pointwise under staggering) closes only
the *two-parent-simultaneous* channel; the *single*-parent leakage measured here
is untouched by staggering, because it is present whenever any one parent is on.

---

## 5. Block-Lyapunov classification (exact) — and why it is *not* good news

With `W = blockdiag(w_i G_i)`, `w_i > 0`, Lemma D0 of the audit forces, for
every edge and **every parent separately** (the parents modulate independently):

    w_i G_i B_ij + w_j B_ji^T G_j = 0,
    i.e.  Q := B_ji^T G_j  =  -(w_i/w_j) * (G_i B_ij) =: -rho_ij P .

This is the honest **block** form of (H2); `rho_ij > 0` is block ellipticity and
`prod_e rho_e = 1` is (H3).

| classification over the analysed list | count |
|---|---|
| orientations classified | **216 192** |
| class (a): common block-scalar metric exists (**dead for all modulations**) | **8** |
| class (b): no such metric (candidate for growth) | **216 184** |
| edge-block pairs with `Q` **not proportional** to `P` at all | the overwhelming majority |
| edge-block pairs proportional with `rho > 0` (block-elliptic) | a small minority; where it occurs, `rho = 1` exactly and holonomy `= 1` |

**Reading this correctly matters.** "No common metric" is *not* evidence that the
cycle works. It says the *certificate class is too narrow* — and §6 shows why:
proportionality already fails on a **single edge**, where (H3) is vacuous. The
failure is (H2) at the block level, and it has nothing to do with cycles.

---

## 6. TREE CONTROL — the decisive structural finding

A single edge `{k, k+p}` with one parent frozen at `|a| = 1` is a **forest**.
Proposition 1 / Theorem D1 assert — *under hypotheses (H1), (H2) and (H3)
jointly* — that the top exponent is `<= -nu min|k|^2` for *every* modulation
and *every* amplitude; forest topology makes (H3) vacuous, but (H1) and (H2)
are still required. Measured directly on the exact Leray blocks in full
polarization (`|p|^2 <= 6`, `|k|^2 <= 12`, all 8 generators):

| `nu` | single-edge systems | with `max Re lambda > 0` | fraction | largest rate |
|---|---|---|---|---|
| `0` | 28 480 | **5 438** | **19.1 %** | **`+3.577709`** |
| `1/40` | 28 480 | **3 750** | **13.2 %** | **`+3.429106`** |
| `1/10` | 28 480 | **2 256** | **7.9 %** | **`+3.000000`** |

The `nu = 0` row is the amplitude-independent statement: `max Re lambda(K) > 0`
means the coupling generator alone is **not elliptic**, so by linearity in the
amplitude *some* finite amplitude beats *any* fixed `nu`. On a forest, Prop 1 —
*when it applies* — forbids exactly this; here it does not apply, because (H2)
fails.

### 6.0 The counterexample in full, verified explicitly

    children   k_a = (0,1,3)   |k_a|^2 = 10
               k_b = (1,0,1)   |k_b|^2 =  2        k_a + p = k_b  (checked)
    parent     p   = (1,-1,-2) |p|^2   =  6,  c = (5,-1,3),  phase cos
               c . p = 0 exactly,   div(U_P) = 0 exactly,   ||U_P||^2 = 35/2

    coupling graph on the child set : blocks (0,1) and (1,0) only
                                    = 2 nodes, 1 edge = a FOREST
    diagonal (self-coupling) blocks : NONE  => K is zero-diagonal,
                                      exactly as Proposition 1 assumes

    (H2) block-ellipticity test :  Q = B_ji^T G_j  vs  P = G_i B_ij
                                   status = NOT PROPORTIONAL  (exact, Fraction)
                                => no scalar rho_ij exists at all,
                                   so no positive block-scalar metric exists

    max Re lambda(K)              = +3.577708764        (elliptic needs <= 0)
    max Re lambda(-nu D + K):
        nu = 0     : +3.577709   vs Prop-1 bound  -0.000000   OUTSIDE HYPOTHESES ((H2) fails)
        nu = 1/40  : +3.429106   vs Prop-1 bound  -0.050000   OUTSIDE HYPOTHESES ((H2) fails)
        nu = 1/10  : +3.000000   vs Prop-1 bound  -0.200000   OUTSIDE HYPOTHESES ((H2) fails)

    (With (H2) failing, Proposition 1 asserts nothing about this system; the
     "Prop-1 bound" column is what it WOULD assert if its hypotheses held.)

**This is not a large-field artifact.** Measured directly on the same blocks
(bisection on the amplitude; no closed-form threshold is quoted because no
valid Weyl-type inequality lower-bounds `max Re lambda` here — the
`-nu min|k|^2` bound goes the wrong way):

    a* = 0.03125   ->   max Re lambda = 0   (crossing),
    i.e. a parent of L^2 norm  ||U_P||_{L^2} = 0.131

— a factor **~30** below the `|a| <= 1` ceiling the gate grants. (An earlier
version claimed `a > 0.013975`, `||U_P|| > 0.0585`, "three orders of magnitude
below the ceiling", from an inequality with the wrong sign; at `a = 0.013975`
the measured `max Re lambda` is `-0.038` — no growth.)

So the exact Leray coupling **grows on a forest** at small amplitude, with no
cycle anywhere — a configuration where Proposition 1, *inside its hypotheses*,
would forbid growth; here (H2) fails, so Prop 1 is outside its hypotheses and
asserts nothing. This bears on residual gap #1 that the audit itself
registered as its top follow-up:

> "**(H1) from NS.** The audit *assumes* the reciprocal structure ... It does
> **not** re-derive it from `N_k = -i P_k sum (m.u_l) u_m`. ... **Registered as
> the top follow-up. If it fails, Prop 1 does not apply to the repo's relays.**"

**Result: Prop 1 gives no control here — but through (H2), not (H1), and the
registered (H1) obligation itself remains open.** This gate does not test
(H1): reciprocity (one parent amplitude driving both directions of an edge)
holds *by construction* of `A(t) = -nu D + sum_a a_{p_a} K_a` — an assertion,
not the registered re-derivation from `N_k`. What dies, measurably, is scalar
ellipticity: the `4 x 4` blocks `G_i B_ij` and `B_ji^T G_j` are simply not
proportional, so (H2) fails. The audit had moreover *already* established that
(H2) failure permits unbounded forest growth (its §2.5: "(H2) is necessary,
and its failure is not a large-amplitude effect"); the genuinely new content
here is narrower: **(H2) fails for the exact Leray blocks in full
polarization**, so Prop 1 gives no control over the repo's relays regardless
of (H1). The physical picture is classical but heuristic here: a single
Fourier parent mode is a Kolmogorov flow, which is linearly unstable above a
finite Reynolds number (Meshalkin-Sinai) — though that is a 2D long-wave
result and the exhibited counterexample (`|k_a|^2 = 10 > |p|^2 = 6`, with
`k - p` and `k + 2p` discarded by the truncation) sits outside its regime. The
note's single-polarization reduction hid the instability by projecting the
unstable direction away.

**Consequence for the idea.** The note's headline is
`"impossible at any amplitude for every tree" -> "finite threshold for a cycle"`,
and it calls that "*factor-infinity-to-finite jump ... the whole idea*". There is
no such jump. Both trees and cycles have finite thresholds. The cyclic carrier
buys, at best, a constant factor.

### 6.1 The note's own B.5 gadget is *not* an instance of the tree failure

Fairness demands the separate check, and the gadget passes it:

    B.5 gadget, children (2,2,1),(2,1,2),(1,2,2) (single shell |k|^2 = 9),
    parents |p|^2 = 2, all 512 orientations:
      max Re lambda of every single-parent generator, nu = 0:   0.000000
      (float64 eigenvalues; max over all 24 generators ~ 8.4e-16, not exact zero)
      => all 512 orientations lie in the ZERO-INSTANTANEOUS-MARGIN class.

The note's premise is **reproduced, in full polarization, not refuted**. The
gadget's kill is the scaling kill of §7, not the tree control.

---

## 7. Monodromy search, and the scaling that kills it

### 7.1 The B.5 gadget does grow at fixed `N` (512 orientations, exhaustive)

Piecewise-constant schedules, `|a| <= 1`, staggered / constant / sign-flip /
random families, interval lengths on a four-octave grid tied to `1/||K||` and to
`1/(nu|k|^2)`; `A(t)` contains `-nu D`, so `chi` below is the **net** exponent.

| `nu` | orientations | `rho(M) > 1` | best `chi` | best `rho` | best frozen single-parent margin |
|---|---|---|---|---|---|
| `0` | 512 | **480 / 512** | `+0.606196` | `1.024217` | `+0.000000` |
| `1/40` | 512 | **151 / 512** | `+0.381196` | `1.015161` | `-0.225000` |
| `1/10` | 512 | **0 / 512** | `-0.293804` | `0.988469` | `-0.900000` |

All 512 are in the zero-instantaneous-margin class at every `nu`. So at
`nu = 1/40` the mechanism delivers exactly what it advertises: **net Floquet
growth with identically zero instantaneous eigenvalue margin.** Step 5 does
**not** kill.

*(Open erratum, M4: this 512-orientation table is not dumped to any artifact —
`stage_B.top_by_chi` keeps only 40 rows, none of them geometry 0. An
independent re-run reproduces `480/512` at `nu = 0` and `0/512` at `nu = 1/10`
exactly, but gives `155/512` rather than `151/512` at `nu = 1/40`: the count
depends on the rng state, which in the archived run had been advanced by
stages A–C. The best-`chi` values, including `chi_1 = 0.606196`, reproduce.)*

### 7.1b The broad screen (382 geometries, 24 000 candidate systems)

Stage B screened **24 000** (geometry, orientation) systems — the top 64
orientations per geometry by the exact imbalance score — at three `nu`, i.e.
**72 000 monodromy searches**, in `829.8 s` (cap 20 min, **not** bound).

| `nu` | candidates | `rho(M) > 1` | in zero-margin class | growth *inside* that class | best `chi` |
|---|---|---|---|---|---|
| `0` | 24 000 | 19 118 | 15 647 | **11 281** | `+2.998115` |
| `1/40` | 24 000 | 11 334 | 20 402 | **7 770** | `+2.825720` |
| `1/10` | 24 000 | 8 112 | 21 236 | **5 567** | `+2.339524` |

Growth with zero instantaneous margin is therefore not a curiosity of the B.5
gadget — it is generic. **That is why the kill had to come from the tree control
and the scaling, not from the absence of growth.**

### 7.2 Critical normalization — the pre-registered step-6/7 question

The gate's real question is not "does it grow at `|k|^2 = 9`" but "does the gain
survive as `N -> infinity` at critical normalization". Dilating the gadget
(`k -> s k`, `p -> s p`, so `N = 3s`) and replacing `|a| <= 1` by the
amplitude `a_s = sqrt(c_E/N)` (note: *not* divided by `||C_a||` — see the
amplitude-normalisation caveat below):

Because the children sit on **one shell**, `-nu D = -nu N^2 I` is scalar and the
split is rigid — `chi_net(s) = chi_coupling(s) - nu N^2` **exactly**, no
inequality. With `chi_coupling` linear in the amplitude and linear in `s` (one
derivative in the advection; linearity in `a s` holds for the supremum over the
scale-covariant schedule family searched, by time rescaling `t -> t/(a s)`),

    chi_coupling(s) = chi_1 * s * sqrt(c_E/N) ,      N = 3s
    =>  q_*(N) := chi_net/(nu N^2) = C(c_E) * N^{-3/2} - 1 ,
        C(c_E) = chi_1 sqrt(c_E) / (3 nu) .

*(Corrected: an earlier version printed `C(c_E) = chi_1 sqrt(c_E/3)/(9 nu)`,
which is the coefficient of `s^{-3/2}`, not of `N^{-3/2}`; with `N = 3s` the
two differ by `3^{3/2}`. With `chi_1 = 0.606196` and `nu = 1/40` the corrected
`C(1) = 8.0826`.)*

**Measured** `q_*` (`nu = 1/40`; each entry an exact-block computation with a
fresh schedule search):

| `c_E` | `N=3` | `N=6` | `N=9` | `N=12` | `N=18` | monotone decreasing? | `q_* = 0` at |
|---|---|---|---|---|---|---|---|
| `1` | `+0.5555` | `-0.4500` | `-0.7006` | `-0.8056` | `-0.8942` | **yes** | `N = 4.0` |
| `10` | `+3.9189` | `+0.7391` | `-0.0534` | `-0.3851` | `-0.6653` | **yes** | `N = 8.7` |
| `100` | `+14.5550` | `+4.4995` | `+1.9936` | `+0.9444` | `+0.0584` | **yes** | `N = 18.7` |
| `228` | `+22.4875` | `+7.3041` | `+3.5202` | `+1.9359` | `+0.5981` | **yes** | `N = 24.6` |

The measured values agree with the corrected closed form
`C(c_E) N^{-3/2} - 1`, `C(c_E) = chi_1 sqrt(c_E)/(3 nu)`, to four decimal
places at every point, and the zero-crossing column above (`N = 4.0 / 8.7 /
18.7 / 24.6`) is recomputed from that corrected constant (values unchanged —
they had been derived from the correct constant all along). **This agreement
is a consistency check, not an independent one**: the searched schedule grid
is built from `deltas ~ 1/max||amp*K||`, i.e. it is scale-covariant by
construction, so `chi_coupling/(amp*s)` is forced to be constant across rows
(it is, to seven digits, in all 20). The physical content is the single-shell
rigidity of §8.2, not the fit.

*(Open erratum, M7 — amplitude normalisation.) The implemented amplitude is
`sqrt(c_E/N)` without the `||C_a||` division, and the scan maximises over all
512 orientations, so the realised parent energy is an orientation-dependent
multiple (up to `4.5x`) of `c_E/N`; the idea note's own critical amplitude
carries a further `sqrt(2/3)`. The exponent `N^{-3/2}` and the strict decrease
are scale statements and are unaffected; every **absolute** threshold in the
table — including the `q_* = 0` column and the headline `N ~ 24.6` — is
therefore an upper estimate that a corrected normalisation could pull down by
up to `4.5^{1/3} ~ 1.65` (to `N* ~ 15` at `c_E = 228`). Not re-run; the
direction of the kill does not depend on it.*

`c_E = 228` is the repository's **own measured** capacity value at `eta = 0.20`,
quoted in the note's §E.4. Even there the mechanism has **no net growth beyond
`N ~ 25`** (subject to the M7 caveat above, possibly `N ~ 15`), and `q_*` is
strictly decreasing at every step.

Because `q_* + 1` scales as `N^{-3/2}`, holding `q_* >= q_0 > 0` uniformly in `N`
requires

    c_E  >=  9 (1+q_0)^2 nu^2 N^3 / chi_1^2   ->   infinity ,

which is exactly `judge_rigor`'s `c_E >~ nu^2 N^3 / chi_par^2 -> infinity`, now
derived rather than estimated, and exactly the note's own **kill condition K2**.

### 7.3 Scope of §7.2, stated precisely

The scan **dilates the 3-mode cycle** (`k -> s k`), so it holds the mode
*cardinality* fixed at 3 and tests the **fixed-cardinality branch**. That the
measured decay is `N^{-3/2}` — numerically identical to the exponent the repo's
own `fixed_cardinality_scaling` already rejects — is a consistency check, and it
closes that branch exactly. It is **not** a test of the fixed-relative-width
cloud (`W = floor(N^gamma)`, `M_eff = theta N^3`), which is obligation 4 / K2 and
which step 6 would have measured.

### 7.4 The pre-registered step-7 GO conditions

| condition | status |
|---|---|
| `inf_N E_next/E_current > 0`, **no decreasing trend** | **FAILS** in the fixed-cardinality branch (strictly decreasing for every `c_E`; negative beyond finite `N`). **Not measured** in the fattened-cloud branch. |
| off-chain/main ratio non-increasing in `N` | **Measured from the stage-A archive: mildly DECREASING in `N`** — median `2.080` at `max\|k\|^2 <= 9` (n = 218) falling to `1.412` at `max\|k\|^2 >= 25` (n = 93). This is the direction the GO condition asks for, i.e. **favourable**. (The all-shell median is `1.814`, B.5 gadget `2.144486...`; an earlier version called the trend unmeasured and the aggregate "already adverse at the smallest shells" — both wrong: the smallest-shell median is `2.08`, and the trend was computable from the archive all along.) |
| second-stage normalized gain non-degenerate | **not measured** (step 6 never implemented) |
| output shape returns to the same cycle class | **not measured** (step 6 never implemented) |

GO requires all four to be *established*. One fails outright in the only
branch tested, one is measured and favourable, and two are unmeasured, so
**GO is not available on this evidence**; combined with GROUND 1, which is
complete and independent of step 6, the gate is **KILL**.

---

## 8. Interval certificate

The rigorous enclosure of the monodromy is built as specified: interval Taylor
scaling-and-squaring per constant piece, with the remainder applied to the
**scaled** matrix `B = A/2^s` *before* the `s` squarings,

    || e^B - T_m(B) ||_inf  <=  ||B||_inf^{m+1} e^{||B||_inf} / (m+1)!

added as a symmetric interval to every entry (entrywise magnitude is bounded by
the induced inf-norm), at 100-bit outward rounding, on the exact rational blocks
and an exactly rationalised schedule.

**One honest correction to the obvious plan.** The determinant route is *useless
here and cannot be repaired*: every `K_a` is traceless on the child band, so

    det M = exp(-nu tr(D) T)  <=  1  exactly,

and `rho(M) >= |det M|^{1/n}` can never exceed 1. (This is the same structural
fact the audit records as "`rho(M) >= 1` is structural, so only the null control
carries information", here in its viscous form.) The bound actually claimed is
the **trace-power bound**, rigorous for every `q`:

    |tr(M^q)| = |sum_i lambda_i^q| <= n rho(M)^q
      =>  rho(M) >= ( |tr(M^q)| / n )^{1/q} ,

evaluated on rigorous enclosures of `M, M^2, M^4, ...` by interval squaring.

### 8.1 The certificate SUCCEEDS — and it certifies the note's phenomenon

Candidate: the B.5 gadget, `nu = 1/40`, orientation `(2, 6, 1)`, the
**staggered** schedule (one parent per interval, `dt = 1`, period `T = 3`) — the
note's own §B.6 pump. Dimension `n = 12`, 100-bit outward rounding.

    single-parent instantaneous margins (INTERVAL-backed generators, float eig):
        [-0.225, -0.225, -0.225]        all strictly NEGATIVE
        => the candidate is inside the zero-instantaneous-margin class

    RIGOROUS enclosures of tr(M^q) and the resulting lower bounds on rho(M):

      q      |tr M^q| >=        enclosure width      rho(M) >=
      1      4.20784            7.43e-26             0.350653224
      2      0.337051           2.33e-25             0.167593435
      4      2.1468             5.85e-25             0.650358040
      8      1.96839            1.29e-24             0.797748712
      16     2.65435            3.01e-24             0.910014989
      32     3.45392            8.24e-24             0.961828982
      64     5.9642             2.97e-23             0.989135528
      128    17.7859            1.85e-22             1.003078930   <-- certified
      256    158.169            3.44e-21             1.010124166
      512    12508.6            5.69e-19             1.013665319
      1024   7.82331e+07        7.44e-15             1.015440548
      2048   3.06021e+15        6.08e-07             1.0163293277  <-- best

    CERTIFIED:  rho(M) >= 1.0163293277  >  1
    (lower bound rounded DOWN from 1.01632932774261757...; an earlier version
     printed 1.016329328, which rounds a lower bound UP and is therefore not
     itself certified).
    (float reference: rho = 1.01721889 — the rigorous bound approaches it from
     below as q grows in this instance; no monotone convergence is guaranteed,
     as the q = 2 row's near-cancellation shows.)

**One non-rigorous link, disclosed.** The final root is taken in float64:
`bound = float(lo/n) ** (1.0/q)` on the exact rational trace lower bound, with
no outward rounding. Everything upstream (the `Interval` arithmetic, `iexpm`,
the remainder, the squarings, the traces) is rigorous; this one line is not.
The induced error here is `~5e-20`, five orders below the certified margin
`1.6e-2`, so `> 1` survives comfortably — but a fully rigorous pipeline would
take the root in outward-rounded rational arithmetic, or publish the rational
`|tr M^q|/n` and let the reader take the root. Relatedly (open erratum M3),
the certificate artifact stores only `float()` images of the trace enclosures;
the exact rational endpoints are discarded at serialisation.

Total cost after the remainder fix: **1.4 s**. The relative width of the `q=2048`
trace enclosure is `2 x 10^-22`, so the certificate is nowhere near its precision
limit.

**This is a positive result, and it should be recorded as one.** It discharges
obligation 3 of the note's own proof chain —

> "*3. Floquet positivity with zero instantaneous eigenvalue margin — interval
> enclosure of `rho(M) > 1` for one exhibited orientation (PO-13 machinery).*"

— **affirmatively**, and at a strictly higher standard than the note reached
(the note had float only, and a single fixed child polarization; this is
interval arithmetic on exact rational blocks in the full 4-dimensional
polarization space). The mechanism's *phenomenon* is real and now rigorously
certified. What is dead is its *significance*, for the two reasons in the
verdict.

### 8.2 What is proved without any floating point at all

> **Exact (`Fraction`), no interval arithmetic needed:** the B.5 children lie on
> a single shell `|k|^2 = 9 s^2`, so `D = N^2 I`, the viscous term is a scalar
> multiple of the identity and commutes with everything, and every Floquet
> exponent of `A(t)` equals the corresponding exponent of the coupling part
> **minus `nu N^2` exactly**. Combined with the exact linearity of `K_a` in the
> parent amplitude and in `s` (and the scale covariance of the optimised
> schedule family, §7.2), this makes `q_*(N) = C(c_E) N^{-3/2} - 1`,
> `C(c_E) = chi_1 sqrt(c_E)/(3 nu)`, an **identity in the model, not a fit** —
> which is why §7.2 is decisive for the fixed-cardinality branch without any
> certificate at all.

---

## 9. Compute caps (pre-registered) and whether they bound

| stage | cap | actual | bound? |
|---|---|---|---|
| A — enumeration + exact blocks + classification | 15 min | `67.2 s` enumeration + `97.5 s` for 382 geometries | no |
| B — monodromy screen | 20 min | `829.8 s`, **24 000** candidate systems, 72 000 monodromy evaluations | no |
| C — refine + interval certificate | 15 min | **YES — bound** (see below; no stage-C wall time was archived) | **yes** |
| D — tree control | — | `128.4 s`, 28 480 single-edge systems x 3 `nu` | no |
| E — critical-normalization scan | — | `85.7 s` | no |

The archived stage walls sum to `67.2 + 97.5 + 829.8 + 128.4 + 85.7 = 1 208 s
~ 20 min`; no `total_wall_seconds` was archived, so the "~50-minute total"
below (which includes the abandoned first stage-C implementation) is not
verifiable from the artifacts. (A fresh run of the reconciled script does
record `total_wall_seconds` and `stage_C_refine`.)

**Two caps bound, and both are reported rather than hidden.**

1. **The geometry cap binds hard**: 200 kept in R1 out of **3 025** distinct
   geometries retained by the capped enumeration (raw paths examined: 24.4
   million; §3.1 — the previously quoted uncapped total is withdrawn), and the
   pre-registered *ordering* left R1's shell coverage at `max|k|^2 <= 4`. The
   analysis was therefore rerun on the pre-registered stratified supplement R2
   (declared before any dynamics) rather than quietly reported on R1.
2. **Stage C bound, and step 6 was never implemented.** The first
   implementation of the interval certificate carried the Taylor remainder as
   an exact `Fraction` with several thousand digits, which made `iexpm`
   pathologically slow; the run was stopped in stage C and the remainder was
   rounded outward to a short dyadic (still a rigorous upper bound), after
   which `iexpm` costs `~1 s`. By then the ~50-minute total budget was spent,
   and **the step-6 Galerkin cloud test was never implemented at all** — the
   script contains no Galerkin-cloud routine, no `W = floor(N^gamma)`, no
   `M_eff = theta N^3`; the only step-6 references are a docstring aside and
   the verdict's honest-gaps string (the archived run's hard-coded
   `step_6_cloud_test_run: false` carried no information either way). §7.3
   states exactly which GO conditions are therefore unmeasured.

---

## 10. What this gate does and does not settle

**Settled.**
1. The exact Leray tensor implementation is validated: full-system energy
   conservation and div-free antisymmetry are **exactly zero**, and the note's
   own §B.1 coefficients `(-1/2, -1/2, +1/3)` are reproduced exactly.
2. Proposition 1 **does not apply** to the exact Leray coupling in full
   polarization: hypothesis **(H2) fails exactly** at the block level, so Prop 1
   asserts nothing about these systems — and **forests grow** (§6), at parent
   `L^2` norm `0.131` (`a* = 0.03125`). The parametric audit's registered top
   follow-up concerned **(H1)** and remains **open** (§6); what is discharged
   here is the narrower fact that (H2) fails for the exact blocks, so Prop 1
   gives no control over the repo's relays regardless of (H1).
3. The note's zero-instantaneous-margin cyclic gain is **real, generic, and now
   interval-certified**: `rho(M) >= 1.0163293277` (lower bound rounded down),
   with float-evaluated instantaneous margins `~ -0.225 < 0` (§8.1; only
   `rho(M)` itself is interval-backed). Obligation 3 of the note's proof chain
   is discharged **affirmatively**, above the note's own standard.
4. In the **fixed-cardinality** branch the gain is **not uniform in `N`**:
   `q_*(N) = C(c_E) N^{-3/2} - 1`, `C(c_E) = chi_1 sqrt(c_E)/(3 nu)`, as an
   identity, strictly decreasing, negative past finite `N` (§7.2, with the M7
   amplitude-normalisation caveat on absolute thresholds). This is kill
   condition K2's fixed-cardinality half and reproduces the repo's own
   `fixed_cardinality_scaling` exponent.
5. Off-cycle leakage: median `1.814` over all shells, B.5 gadget `2.144486...`
   (float images of exact ratios; §4) — and, as a trend, mildly **decreasing**
   in `N`: median `2.080` at `max|k|^2 <= 9`, `1.412` at `max|k|^2 >= 25`
   (favourable to the cloud branch; §7.4).

**Not settled, and not claimed.**
- **Step 6 (the fattened Galerkin cloud) was never implemented** — no such
  code exists in the script; the compute budget was spent before it was
  written. The fixed-relative-width branch (`M_eff = theta N^3`, obligation
  4 / K2) is therefore **open**, and two of the four step-7 GO conditions are
  unmeasured (§7.4; the leakage-ratio trend, previously called unmeasured, is
  in fact measured and favourable). GROUND 1 of the verdict does not depend on
  any of this; GROUND 2 does.
- **Autonomy.** The parent amplitudes are prescribed. Nothing here shows a
  self-consistent flow can realise any of these schedules. Flagged in advance,
  open regardless of the verdict.
- The `4 x 4`-block truncation to the cycle band is still a truncation; §4
  quantifies what it discards but does not evolve it. In particular the tree
  control and the cycle claim are compared *within the same truncation
  convention* — which is the fair comparison, but is not a PDE statement.
- The tree-control survey covers `|p|^2 <= 6`, `|k|^2 <= 12` only.
- The geometry sample is 382 of the 3 025 distinct geometries retained by the
  capped enumeration (§3.1); the uncapped population size was not archived.
- `chi_1` in §7.2 is a float64 Floquet exponent: the `N^{-3/2}` *exponent* is
  exact (single-shell rigidity + linearity), the *prefactor* is float.
- **Nothing here is evidence for or against finite-time singularity formation.**

---

## 11. Artifacts

All paths below are relative to the repository root.

| file | size | content |
|---|---|---|
| `experiments/run_exact_leray_cycle_gate.py` | — | the whole gate; the pre-registered bounds are frozen in `REGISTRY` at the top |
| `outputs/verification_sprint_v1/cycle_gate/cycle_gate.json` | 0.08 MB | consolidated report: registry, correctness checks, enumeration, **geometry summary only** (full per-geometry tables are in `cycle_gate_stage_a.json`), monodromy screen, tree control, scaling scan, interval certificate, **verdict**. The archived copy predates the script reconciliation: it lacks `stage_C_refine` and `total_wall_seconds`, which a fresh run now records. |
| `outputs/verification_sprint_v1/cycle_gate/crosstalk_tables.json` | 19.9 MB | **the cross-talk tables** — 34 688 exact `4 x 4` blocks, every off-cycle leakage coefficient as a `Fraction` |
| `outputs/verification_sprint_v1/cycle_gate/cycle_gate_stage_a.json` | 0.8 MB | per-geometry exact edge/holonomy tables and classification |
| `outputs/verification_sprint_v1/cycle_gate/cycle_gate_stage_b.json` | 0.9 MB | monodromy screen, per-candidate single-parent margins |
| `outputs/verification_sprint_v1/cycle_gate/cycle_gate_stage_de.json` | 9 KB | tree control + critical-normalization scan |
| `outputs/verification_sprint_v1/cycle_gate/cycle_gate_interval_certificate.json` | 5 KB | the interval certificate, all `q` levels up to 2048 |

**Reproduce the certificate (cheap, ~2 s):**

    .venv/Scripts/python.exe experiments/run_exact_leray_cycle_gate.py --stage certificate

re-runs only the certified candidate (B.5 gadget, orientation `(2, 6, 1)`,
staggered `dt = 1` schedule, `nu = 1/40`, `--squarings 11`, i.e. `q` up to
2048) and rewrites `cycle_gate_interval_certificate.json`. **Verified
2026-08-02: the reproduced file is byte-identical to the archived artifact**
(`best_rho_lower_bound = 1.0163293277426175`, `certifies_rho_gt_1 = true`,
1.4 s).

**Reproduce the full pipeline (~20 min):**

    .venv/Scripts/python.exe experiments/run_exact_leray_cycle_gate.py

An earlier committed version of the script defaulted to `squarings = 6`
(`q <= 64`, best bound `0.98914 < 1`), which could not reproduce the archived
certificate and would have flipped the script's step-5 flag; the default is
now 11 and `--squarings` is exposed. A fresh full run additionally stores
`stage_C_refine`, `total_wall_seconds`, and the certificate's
`float_reference`, and writes the leakage-trend medians into the verdict
block; its `verdict`/`stage_E_verdict` prose strings carry the corrections of
this revision (rounded-down bound, "never implemented", measured trend), so
they differ from the archived copies in exactly those respects.

---

## 12. One-line summary for the sprint ledger

> The cyclic-carrier parametric mechanism produces **real, interval-certified**
> Floquet growth with negative (float-evaluated) instantaneous margins
> (`rho(M) >= 1.0163293277`) — and it does not matter, because a **single tree
> edge** does the same thing (Prop 1's hypothesis (H2) fails exactly on the real
> Leray blocks, so Prop 1 does not apply), and because the gain decays as
> `N^{-3/2}` under critical normalization. **KILL.**
