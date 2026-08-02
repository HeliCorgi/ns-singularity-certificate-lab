# Adversarial referee report — `exact_leray_cycle_gate.md`

**Document under review:** `docs/research_notes/verification_sprint_v1/exact_leray_cycle_gate.md`
**Artifacts:** `outputs/verification_sprint_v1/cycle_gate/{cycle_gate,cycle_gate_stage_a,cycle_gate_stage_b,cycle_gate_stage_de,cycle_gate_interval_certificate,crosstalk_tables}.json`
**Script:** `experiments/run_exact_leray_cycle_gate.py`
**Cross-referenced:** `parametric_common_lyapunov_audit.md`, `../ideas_2026_08_01/idea_parametric_resonance.md`
**Method:** every quoted number checked against the stored JSON; the tree-control
survey re-run in full; the §6.0 counterexample re-derived from scratch through an
independent reimplementation on `fourier_torus`; the §7.1 512-orientation sweep
re-run; the §7.2 scaling law refitted from the 20 stored rows; the interval
machinery read line by line.

Defects only. Correct steps are not commented on except where a claim about them
is wrong.

---

## CRITICAL — invalidates a statement as written

### C1. §7.2 / §8.2 / verdict GROUND 2 — the stated constant `C(c_E)` is wrong by `3^{3/2}`, and the "agrees to four decimal places" claim is false

The document prints

    q_*(N) = C(c_E) * N^{-3/2} - 1 ,   C(c_E) = chi_1 sqrt(c_E/3) / (9 nu)

and asserts "The measured values agree with the closed form `C(c_E) N^{-3/2} - 1`
to four decimal places at every point."

`C = chi_1 sqrt(c_E/3)/(9 nu)` is the coefficient of `s^{-3/2}`, **not** of
`N^{-3/2}`. With `N = 3s`, `s^{-3/2} = 3^{3/2} N^{-3/2}`, so the correct constant is

    C(c_E) = chi_1 sqrt(c_E) / (3 nu)   =   3^{3/2} x (the printed C).

Verified against all 20 rows of `stage_E_critical_normalization.rows`. With the
measured `chi_1 = 0.6061964` and `nu = 1/40`, the correct `C(1) = 8.0826`
reproduces every stored `stage_gain_q_star` to 4–5 decimals
(`N=3: +0.5555`, `6: −0.4500`, `9: −0.7006`, `12: −0.8056`, `18: −0.8942`).
The printed `C(1) = 1.5555` reproduces the value at `3N` instead of at `N`:
`1.5555·3^{-3/2} − 1 = −0.7006`, which is the stored `N = 9` entry, not `N = 3`.
The agreement claim is therefore false as written, off by a factor `5.196`.

**Internal inconsistency confirming the sign of the error.** The document's own
downstream consequence,

    c_E >= 9 (1+q_0)^2 nu^2 N^3 / chi_1^2 ,

and its own zero-crossing column (`N = 4.0 / 8.7 / 18.7 / 24.6`) both follow from
the **correct** `C`, not from the printed one. The printed `C` would give
`c_E >= 243 (1+q_0)^2 nu^2 N^3/chi_1^2` and crossings at `N = 2.3 / 5.0 / 10.8 /
14.2`. Two of the document's three formulas are mutually contradictory.

**Minimal fix.** Print `C(c_E) = chi_1 sqrt(c_E)/(3 nu)`, or state the law in `s`:
`q_* = [chi_1 sqrt(c_E/3)/(9 nu)] s^{-3/2} - 1`.

### C2. The committed script cannot produce the certificate, and run as documented it would invert the document's headline

`main()` calls `interval_certificate(Kexact, dvec, nu_f, sched, n)` — with the
default `squarings: int = 6` (line 1040). That yields `q ∈ {1,2,4,8,16,32,64}`
only. The stored certificate contains `q` up to **2048** (i.e. `squarings = 11`),
plus a `float_reference` key the script never writes, and it **lacks** the
`attempted` and `candidate` keys the script always writes.

Consequence: running `.venv/Scripts/python.exe experiments/run_exact_leray_cycle_gate.py`
as instructed in §11 stops at `q = 64`, whose best bound is `0.98914 < 1`, so
`certifies_rho_gt_1 = False`; the script's verdict logic then sets

    step_5_kill = not (any_growth and certified) = True

— the exact opposite of the document's headline "**Not killed at step 5 — the
opposite.**" The script is also inconsistent with the stored report in three
further places: `cycle_gate.json` contains `geometries_summary` where the script
writes `geometries`, contains no `stage_C_refine` key at all (the script always
writes one), carries a `verdict` block whose ten keys are entirely different from
the script's, and carries `stage_E_verdict` keys
(`q_star_strictly_decreasing_in_N_for_every_c_E`, `branch_tested`) that the script
does not emit. `REGISTRY` is the only object that matches byte-for-byte.

**Fix.** Commit the script that was actually run, or expose `--squarings` and
record it in the report; re-state §11 "Reproduce" accordingly.

### C3. §6 — the "not a large-field artifact" threshold is numerically false and rests on an inequality with the wrong sign

The document states

    a > nu min|k|^2 / max Re lambda(K) = 0.013975 ,
    i.e. a parent of L^2 norm above 0.0585 ,

and calls this "three orders of magnitude below the `|a| <= 1` ceiling".

Recomputed on the same blocks (independent reimplementation, `nu = 1/40`,
`k_a = (0,1,3)`, `k_b = (1,0,1)`, `p = (1,-1,-2)`, `c = (5,-1,3)`):

    a = 0.005000  ->  max Re lambda = -0.048413
    a = 0.010000  ->  -0.043793
    a = 0.013975  ->  -0.038197        <-- the document's claimed threshold: NO growth
    a = 0.020000  ->  -0.027037
    a = 0.050000  ->  +0.054939
    bisected crossing:  a* = 0.031250  ->  ||U_P||_{L^2} = 0.130728

The claimed threshold is wrong by a factor `2.24`. The derivation implicitly uses
`-nu D >= -nu min|k|^2 I`, which is **false**: `D >= min|k|^2 I` gives
`-nu D <= -nu min|k|^2 I`. The only Weyl-type bound available goes the other way,
`lambda_max(sym(-nuD + aK)) >= a lambda_max(sym K) - nu max|k|^2`, i.e. `max|k|^2 = 10`
and `a >= 0.0699` — and even that bounds the symmetric part, which upper-bounds
`max Re lambda` and therefore does not lower-bound it at all. No valid version of
the quoted formula exists; the number has to be measured.

"Three orders of magnitude below the ceiling" is also wrong for any of these
values: `0.0585` is 1.2 orders below 1, `0.1307` is 0.9 orders.

**Fix.** Quote the measured `a* = 0.03125`, `||U_P||_{L^2} = 0.1307`, delete the
formula, and say "a factor ~30 below the granted ceiling".

### C4. §3.1 is contradicted by its own stored output

Document: "**every member of R1 is a 3-cycle sitting at `max|k|^2 <= 4`** … it
contains **no 4-cycle and no double-triad at all**. As a sample of the hypothesis
space it is worthless", and the section heading "R1 alone is **degenerate**".

`cycle_gate.json → enumeration.selection` records:

    R1_by_kind        = {3-cycle: 73, 4-cycle: 62, double-triad: 65}
    R1_max_k_sq_range = [1, 4]
    R1_DEGENERATE     = false

R1 contains 62 four-cycles and 65 double-triads. The stated ground for
introducing the R2 supplement — R1 degeneracy in cycle *kind* — does not hold in
the archived run. (The per-`(kind, shell)` bucket cap of 24, which the document
lists in its own registry table, is exactly what prevents the degeneracy it
describes.) The narrative appears to describe the pre-bucket-cap pass while being
written as a statement about the archived R1.

**Fix.** Restate §3.1 against the archived `R1_by_kind`, and justify R2 on shell
coverage alone (`R1 max|k|^2 <= 4` vs the registered bound 50) — which is a
sufficient and true justification.

### C5. §3.1's headline enumeration figures appear in no artifact

"An uncapped pass … finds **6 535 748** distinct cycle geometries — **259 205**
3-cycles, **5 922 588** 4-cycles, **353 955** double-triads", from which the
document derives "the cap binds by a factor of `3.3 x 10^4`" and §10's "382 of
6 535 748".

No stored file contains 6535748, 259205, 5922588 or 353955. What is stored is

    enumeration.raw_paths_examined     = {3-cycle: 1 094 958, 4-cycle: 20 973 928, double-triad: 2 343 120}
    enumeration.buckets_that_overflowed= {3-cycle:   777 378, 4-cycle: 11 950 156, double-triad: 1 412 724}
    enumeration.distinct_geometries_retained = 3 025

No combination of these gives 6 535 748 or its three-way split (raw − overflow =
317 580 / 9 023 772 / 930 396). The load-bearing "worthless sample" argument and
the sampling fraction in §10 rest on an unarchived side computation. Note also
that `raw_paths_examined` counts *raw paths* (incremented before the
distinct/valid filters), so it cannot be quoted as "distinct geometries" either.

**Fix.** Archive the uncapped pass, or delete the figure and quote the binding
factor against the archived 3 025.

---

## MAJOR — gap needing repair

### M1. The certified value is rounded the wrong way

Every occurrence (verdict box, §8.1 "CERTIFIED", §10 item 3, §12) reads
`rho(M) >= 1.016329328`. The rigorous bound from the stored trace enclosure is

    (|tr M^2048| / 12)^{1/2048} = 1.016329327742617576...

and `1.016329328^2048 * 12 > |tr M^2048|_lower` (checked in 60-digit decimal), so
the printed inequality is **not** certified — it rounds a lower bound *up*. A
certified lower bound must be rounded down: `rho(M) >= 1.0163293277`.

### M2. The last step of the certificate is not rigorous

`interval_certificate` computes `bound = float(lo / n) ** (1.0 / q)` — a float64
root of an exact rational, with no outward rounding. Everything upstream
(`Interval`, `iexpm`, the remainder, the squarings, the traces) is rigorous; this
one line is not. Here the induced error is `~5e-20`, so `> 1` survives, but the
document presents the printed 9-digit figure as "CERTIFIED".
**Fix.** Take the root in outward-rounded rational arithmetic, or publish the
rational `|tr M^q|/n` and let the reader take the root.

### M3. The certificate is not independently checkable from the artifacts

`cycle_gate_interval_certificate.json` stores only `float()` images of the trace
enclosures (`trace_enclosure`, `abs_trace_lower_bound`); the exact rational
endpoints are discarded. The exact `K_a` blocks for the certified orientation are
not stored anywhere, and `stage_C_refine` — which produced the certified schedule
and the `rho = 1.01721889` float reference — is absent from *every* stored file.
`cycle_gate.json` has no `stage_C_refine` key and no `total_wall_seconds`.
Everything downstream of stage B must be taken on trust, including §9's
"Stage C bound: **YES**", for which no wall time was recorded.

### M4. §7.1 — the whole B.5 table, and with it `chi_1`, has no artifact

The 512-orientation table (`480/512`, `151/512`, `0/512`; `chi = +0.606196 /
+0.381196 / -0.293804`; `rho = 1.024217 / 1.015161 / 0.988469`) is the source of
`chi_1 = 0.606196`, which is the entire prefactor of GROUND 2. It is in no output
file: `stage_B.top_by_chi` keeps only 40 rows and contains **zero** geometry-0
entries, and no stored file contains any of those six numbers.

I re-ran it. Reproduced: `480/512` at `nu = 0` with best `chi = 0.6061964`;
`0/512` at `nu = 1/10` with best `chi = -0.2938036`. **Not** reproduced:
`155/512` at `nu = 1/40`, against the document's `151/512` — the count depends on
the rng state, which differs because in the real run `rng` had been advanced by
stages A–C. The `151` should be labelled rng-dependent, and the sweep dumped.

### M5. §6 and the verdict misstate Theorem D1 — nothing is "VIOLATED"

The document asserts, in the verdict box and in §6:

> "A *single edge* `{k, k+p}` … is a **forest**, where Proposition 1 asserts no
> growth at *any* amplitude or waveform."
> "Proposition 1 / Theorem D1 assert the top exponent is `<= -nu min|k|^2` for
> *every* modulation and *every* amplitude."

`parametric_common_lyapunov_audit.md` §2.3 and §6 state Theorem D1 / Proposition 1′
under **(H1) and (H2) and (H3) jointly**. Forest topology only makes **(H3)**
vacuous ("A forest has no cycles, so (H3) is vacuous", audit line 120); it says
nothing about (H2). With (H2) failing, Proposition 1 asserts *nothing* about the
system, so the three `VIOLATED` labels in the §6.0 box are category errors, as is
§10 item 2 "Proposition 1's *conclusion* is **false**". The document states the
correct version once ("Prop 1 simply does not apply to the exact Leray coupling")
and the incorrect version five times.
**Fix.** Replace `VIOLATED` with `OUTSIDE HYPOTHESES ((H2) fails)`, and §10.2 with
"Proposition 1 does not apply".

### M6. The registered follow-up that §6/§10.2 claims to discharge is about (H1), and was not tested

The audit's registered top follow-up, quoted verbatim in §6, is: "**(H1) from NS.**
The audit *assumes* the reciprocal structure … It does **not** re-derive it from
`N_k = -i P_k sum (m.u_l) u_m`." Nothing in the gate tests (H1). The document's
reply — "Reciprocity survives (one parent amplitude does drive both directions of
an edge)" — is true by construction of `A(t) = -nu D + sum_a a_{p_a} K_a`, but it
is an assertion, not a measurement, and it is not the re-derivation that was
registered. Claiming the follow-up is "discharged — **negatively**" through (H2)
conflates two different hypotheses; the (H1) obligation remains open.

Separately, the audit had **already** established that (H2) failure produces
unbounded forest growth ("**(H2) is necessary, and its failure is not a
large-amplitude effect**", §2.5; measured `+589.1` in the §4 table). GROUND 1's
genuine content is narrower than advertised: it is that (H2) fails *for the exact
Leray blocks in full polarization*, not that (H2) failure implies growth.

### M7. §7.2 — the implemented amplitude is not the documented critical shell law

Document: "replacing `|a| <= 1` by the critical shell law `E_p = c_E/N`, i.e.
`a_s = sqrt(c_E/N)/||C_a||`".
Code (`critical_normalization_scan`): `amp = math.sqrt(c_E / N)`, then
`Ks = [amp * blocks_to_matrix(...)]`. There is **no** `||C_a||`.

For the B.5 parents the four polarizations give `||C_a||^2 ∈ {1, 1/2, 3/2, 3/2}`
(e.g. `p = (0,1,-1)`: `pol0 = (0,1,1)`, `pol1 = (1,0,0)`, `pol2 = (1,1,1)`,
`pol3 = (1,-1,-1)`), and the scan maximises over all 512 orientations, so the
realised parent energy is `amp^2 * sum_a ||C_a||^2` — up to `4.5 amp^2`, not
`c_E/N`. Moreover the idea note's own critical parent amplitude (line 218) is
`X_P = sqrt(2 c_E / 3N)`, a further factor `sqrt(2/3)`.

Consequence: the exponent `N^{-3/2}` and the strict decrease are unaffected (both
are scale statements), but every **absolute** number in the §7.2 table — the
`q_* = 0` column `N = 4.0 / 8.7 / 18.7 / 24.6`, and the headline "no net growth
beyond `N ~ 25` even at the repository's own `c_E = 228`" — corresponds to an
unstated multiple of `c_E`. Correcting by the (orientation-dependent) factor moves
`N*` by up to `4.5^{1/3} ≈ 1.65`, i.e. `N* ≈ 15` rather than `24.6`. The
direction of the kill is unchanged; the quoted thresholds are not defensible as
printed.

### M8. §7.2 — "identity, not a fit" and "independent check" are both overstated

(i) "`chi_coupling` is exactly **linear in the amplitude** and in `s`" is false for
a fixed schedule: `rho` of a product of exponentials of non-commuting matrices is
not homogeneous in `a`. What is homogeneous of degree 1 in `a s` is the
**supremum over all piecewise-constant schedules**, by the time rescaling
`t -> t/(a s)`. That requires the searched schedule set to be closed under this
rescaling — an unstated hypothesis.

(ii) The document then claims the 20 measured rows are an independent check
("each entry an independent exact-block computation with a fresh schedule search,
*not* an evaluation of the formula … which is the check that the identity of §7.2
is real and not a coincidence of one scale"). It is not independent. The scan
builds `deltas = [f/sc for f in (0.0625,0.25,1.0,4.0)]` with
`sc = max ||amp*K||_inf ∝ a s`, and draws random schedules from the same grid, so
the searched set is **exactly** scale-covariant and the law is forced by
construction. I confirmed the degeneracy: `chi_coupling/(amp*s) = 0.6061964` to
seven digits in **all 20 rows**. The agreement measures the covariance of the
search grid, not the physics.
**Fix.** State the covariance hypothesis, and either drop the "independent check"
sentence or add a genuinely non-covariant schedule family.

### M9. §7.4 / §4 — the off-chain/main trend in `N` is claimed unmeasured, but is in the stored data and points the other way

§7.4 records "off-chain/main ratio non-increasing in `N` | **Not measured as a
trend.**", and the verdict calls the median `1.814` "evidence against the cloud
branch". `cycle_gate_stage_a.json` stores `crosstalk_to_inset_ratio_float` for all
375 geometries together with `max_k_sq`, so the trend *is* measurable from the
archived data. Computed:

    max|k|^2 <= 9  (n = 218) : median 2.080
    max|k|^2 >= 25 (n =  93) : median 1.412

i.e. the ratio is mildly **decreasing** in `N` — the direction the GO condition
asks for. Reporting the aggregate median as "already adverse" while declining a
computation the archive already supports is selective.

Also, "at the *smallest* shells" (verdict box and §7.4) is inaccurate: the median
`1.814` is over all shells `1..50`; restricted to the smallest shells it is `2.08`.

### M10. §9's account of why step 6 was not run is not accurate

"the ~50-minute total budget was spent, so **the step-6 Galerkin cloud test was
never executed**" implies implemented code that was not reached. There is no step-6
implementation anywhere in the script: no Galerkin-cloud routine, no
`W = floor(N^gamma)`, no `M_eff = theta N^3`. The only occurrences of "step 6" are
a docstring aside ("the Galerkin cloud test (step 6, only if reached)"), a comment
on stage E, and a hard-coded `"step_6_cloud_test_run": False`.
**Fix.** State that step 6 was never implemented. (The *conclusion* of the honest-gap
section — that the fixed-relative-width branch is unsettled — is correct and is
correctly reflected in `verdict.honest_gaps` of the stored report.)

---

## MINOR — wording / precision

- **m1. §2.** "Divergence-freeness of every polarization | **18** reference fields
  (12 child basis + 6 parent generators)". Stored
  `correctness_checks.divergence_free.fields_checked = 36` (12 children +
  3 parents x 8 generators). 18 is the *energy* check's basis, not the divergence
  check's.
- **m2. §6.1 / verdict box / §10.3.** "max Re lambda of every single-parent
  generator, `nu = 0`: **0.000000 EXACTLY**" is a float64 statement
  (`np.linalg.eigvals`), and the document's own arithmetic-label table puts margins
  under **FLOAT**. Recomputed maximum over all 24 B.5 single-parent generators:
  `8.40e-16`, not exactly zero. Likewise the verdict box presents
  `-0.225, -0.225, -0.225` inside a sentence headed "rigorously certified by
  interval arithmetic", and §10.3 says "interval-certified: `rho(M) >= …` **with
  all instantaneous margins `= -0.225 < 0`**". Only `rho(M)` is interval-backed;
  §8.1 discloses "float eig" but the verdict box and §10.3 do not. (The substance
  is safe — `D` is scalar on the single shell, so the margin is
  `max Re lambda(K_a) - nu N^2 ≈ -0.225` robustly — but the label is wrong.)
- **m3. §9.** Stage D "`~250 s`": the stored `wall_seconds` is **128.4**
  (`cycle_gate_stage_de.json` and `cycle_gate.json`). Stage E "fast" = 85.7 s.
  Archived totals are `67.2 + 97.5 + 829.8 + 128.4 + 85.7 = 1 208 s ≈ 20 min`; no
  `total_wall_seconds` is stored, so "the ~50-minute total budget was spent" and
  "Stage C — **YES, bound**" are both unverifiable from the artifacts.
- **m4. §11.** `cycle_gate.json` is listed at `0.9 MB`; actual size **0.081 MB**.
  (Others check out: crosstalk 19.87 MB, stage_a 0.84 MB, stage_b 0.90 MB,
  stage_de 9 KB; the certificate is 5 KB, listed as 6 KB.)
- **m5. §11.** `cycle_gate.json` is described as the "consolidated report:
  … classification, monodromy screen, …". Its per-geometry tables were replaced by
  `geometries_summary = {"count": 382, "full_per_geometry_tables_in":
  "cycle_gate_stage_a.json"}`, and `stage_C_refine` is absent (see M3). The table
  should say so.
- **m6. §7.1b.** "the **top 64** orientations per geometry": `382 x 64 = 24 448`
  but `24 000` were screened; 448 (geometry, orientation) systems were dropped
  (class-(a), or `coupling_inf_norm == 0`). "At most 64" is the accurate phrasing.
- **m7. §8.** The remainder is written `|| e^A - T_m(A) ||_inf <= ||A||^{m+1}
  e^{||A||}/(m+1)!`, but the code applies it to the **scaled** matrix
  `B = A/2^s` *before* the `s` squarings (`iexpm`, `bn = imat_inf_norm(B)`). As
  written with `A` the bound is far weaker (still valid, but not what was done).
  The implementation itself is sound: `fact = 27!` from `range(1, order+2)` with
  `order = 26`; entrywise magnitude `<= ||.||_inf` is correct;
  `_exp_upper(r <= 1/4) = 4/3 >= e^{1/4} = 1.2840`; the remainder is rounded
  outward *up*; `Interval` is exact-rational with outward rounding.
- **m8. §8.1.** "the rigorous bound converges to it from below as `q` grows, **as
  it must**". `(|tr M^q|/n)^{1/q}` need not converge, let alone monotonically —
  the document's own table shows `q=2` giving `0.1676` below `q=1`'s `0.3507`,
  because the trace can nearly cancel. Only `limsup_q` recovers `rho(M)`, and only
  when the dominant eigenvalue does not cancel in the trace. Delete "as it must".
- **m9. §5.** "class (a): common block-scalar metric exists | **8**". All 8 come
  from a *single* geometry, stage_a index 1, the trivial 3-cycle
  `{(0,0,1),(0,1,0),(1,0,0)}`. Worth stating: the "certificate class is too
  narrow" reading rests on one degenerate positive example.
- **m10. §5, last table row.** "edge-block pairs proportional with `rho > 0` … a
  small minority; **where it occurs, `rho = 1` exactly and holonomy `= 1`**". The
  archive records 25 376 `proportional` pairs (against 818 656 `not_proportional`)
  but does not split them by sign of `rho`; only 32 `rho` values survive in the
  truncated `example_edge_table` fields (all `'1'`), and exactly one geometry ever
  has a defined holonomy (value `'1'`). The claim is asserted at a generality the
  artifacts do not support.
- **m11. §6.0.** The `nu = 1/10` row is presented as part of the
  `(0,1,3)/(1,0,1)` counterexample, but the archived `nu = 1/10` worst example is a
  *different* pair, `(0,3,1)/(1,1,0)` with `p = (1,-2,-1)`. (The value coincides —
  I recomputed `+3.0000000000000013` for the quoted pair — so this is presentation
  only.)
- **m12. §4 / §10.5.** "median `1.814`, B.5 gadget `2.144486...`, **exactly**". The
  underlying coefficients are exact `Fraction`s, but the ratios stored and quoted
  are float64 (`crosstalk_to_inset_ratio_float`). "Float image of an exact ratio"
  is the accurate label.
- **m13. Verdict box.** "in the **full** 4-dimensional divergence-free
  polarization space of every child mode … `rho(M) >= …` (CERTIFIED)". The
  certificate is for the 12-dimensional truncation to the three cycle modes, which
  §4 shows discards `2.14x` more response energy than it retains. §10 states this;
  the verdict box, which is what will be quoted, does not.
- **m14. §6 / verdict box.** The physical attribution — "a single Fourier parent
  mode is a Kolmogorov flow, which is linearly unstable above a finite Reynolds
  number (Meshalkin–Sinai)" — is decorative and not established here.
  Meshalkin–Sinai is a 2D result about the *full* linearised operator in the
  long-wave regime `|k| << |p|`; the exhibited counterexample is a 2-mode
  truncation with `|k_a|^2 = 10 > |p|^2 = 6` in which `k - p` and `k + 2p` are
  discarded. Note in particular that the "forest" of §6 exists **only because of
  that truncation**: the band closed under `k -> k ± p` is a chain, not a single
  edge. §10 concedes the truncation in general but not this specific point.
- **m15. §7.2.** The shell law is applied at the *child* shell `N = 3s`, while the
  parents sit at `|p| = s sqrt(2) = N sqrt(2)/3`. Unstated modelling choice, feeding
  into M7's prefactor.
- **m16.** `select_lists` hard-codes `"B5_reference_forced_in": True` regardless of
  whether the forcing branch executed. It is a literal, not a check. (The forcing
  did happen — B.5 is stage_a geometry 0 — but the field carries no information.)
- **m17.** The script's own `REGISTRY["GO_conditions_step7"][0]` is
  "`inf_N E_next/E_current > 0` with no decreasing trend", whereas §7.2 measures
  `q_* = chi_net/(nu N^2)`. The substitution is legitimate — the idea note
  pre-registers exactly this quantity at its line 370, "`chi_par(N) s X_P(N) - nu N^2
  >= q_* nu N^2` with `q_* > 0` uniform in `N`" — but the document should say that
  the registry's phrasing and the measured quantity are the note's stage budget,
  not an energy ratio.
- **m18.** `stage_E_verdict.inf_N_stage_gain_positive` is computed as
  `all(row > 0)` over five values of `s`. It is not an infimum over `N`; the name
  invites over-reading. (It is `False` in the archive, so nothing is claimed from
  it.)

---

## Verified — checks I ran that came out clean

Listed so the referee's coverage is auditable, not as endorsement.

1. **Tree control fully re-run** (`tree_control([0, 1/40, 1/10])`): reproduces the
   archive exactly — `28 480 / 28 480 / 28 480` tested, `5 438 / 3 750 / 2 256`
   positive, largest rates `3.577708764 / 3.429106034 / 3.000000000`, same worst
   examples.
2. **Enumeration arithmetic of the tree survey**: `|shell(6)| = 40`,
   `|shell(12)| = 89`, `40 x 89 x 8 = 28 480`. No double counting — the number of
   distinct `(parent, unordered child pair)` couples is exactly 3 560, so the
   duplication factor is 1.0.
3. **§6.0 counterexample re-derived from scratch** on `fourier_torus` (independent
   reimplementation of the basis, Gram, Leray linearisation and block projection):
   `(0,1,3) + (1,-1,-2) = (1,0,1)`; `pol3 = (5,-1,3)` with `c.p = 0` and
   `||C||^2 = 35/2`; blocks present are exactly `(0,1)` and `(1,0)`, no diagonal
   block, `tr K = 0`; `max Re lambda = 3.5777087640 / 3.4291060336 / 3.0000000000`.
   **(H2) confirmed to fail exactly**: with `P = G_i B_ij`, `Q = B_ji^T G_j`,
   `Q[0][2]/P[0][2] = (25/4)/(-9/2) = -25/18` while
   `Q[0][3]/P[0][3] = (3/2)/(-13/2) = -3/13`. Non-proportional over `Fraction`,
   for both `cos` and `sin` phases.
4. **§4 statistics** recomputed from `cycle_gate_stage_a.json`: 375 geometries with
   both channels nonzero, `min 0.3041603`, `median 1.8139269`, `max 5.8697687`,
   B.5 `2.1444864`. Match.
5. **§5 counts**: `216 192 = 8 + 216 184`; edge-status totals `818 656`
   not-proportional vs `25 376` proportional. Match.
6. **§2 counts**: `895 488` reality entries; 18-field energy basis, `5 832` ordered
   triples, both defects exactly `"0"`. Match (except m1).
7. **§7.1b stage-B table**: `24 000 / 72 000 / 829.8 s`; `19 118 / 11 334 / 8 112`
   with `rho > 1`; `15 647 / 20 402 / 21 236` in the zero-margin class;
   `11 281 / 7 770 / 5 567` growing inside it; best `chi`
   `2.998115 / 2.825720 / 2.339524`. All match.
8. **Crosstalk dump**: `34 688` rows, all `4 x 4` `Fraction` strings; the write cap
   (`[:60000]`) did not bind.
9. **Interval machinery audited and sound**: `Interval` is exact-rational with
   outward `floor/ceil` rounding; `imat_mul` skips only exact zeros;
   `||e^B - T_m(B)||_inf <= ||B||^{m+1} e^{||B||}/(m+1)!` is the correct bound and
   is correctly implemented with `27!`; entrywise `<= ||.||_inf` is valid; the
   remainder is added *before* the squarings; `rho(M) >= (|tr M^q|/n)^{1/q}` is
   valid for every `q`. The only non-rigorous link is M2.
10. **Determinant remark checked**: `K_a` traceless confirmed numerically
    (`tr K = 0`), and `det M = exp(-nu tr(D) T) = exp(-8.1) = 3.0354e-4`, matching
    the stored `det_analytic`. The statement that the determinant route cannot
    certify `rho > 1` is correct.
11. **Single-shell rigidity** (the load-bearing structural claim of §8.2) is
    correct: `-nu D = -nu N^2 I` commutes with everything, so
    `chi_net = chi_coupling - nu N^2` exactly. Independently confirmed on the
    archive: for all four `c_E` at `s = 1`, `chi_net + 0.225` divided by `amp`
    equals `0.6061964` to seven digits; and across `s = 1,2,3,4,6`,
    `chi_coupling/(amp*s) = 0.6061964` in all 20 rows.
12. **Honest-gaps section matches reality** (referee item 2): step 6 was indeed not
    run — and, per M10, never implemented; autonomy is indeed untested — the parent
    amplitudes are prescribed everywhere in the script, there is no dynamical-parent
    path at all. `verdict.honest_gaps` in `cycle_gate.json` says both.
13. **§7.2 exponent** (referee item 4): the `-3/2` exponent is correct and follows
    from single-shell rigidity plus the scale covariance of the optimised exponent.
    Only the constant (C1), the amplitude normalisation (M7) and the "independent
    check" framing (M8) are defective.

---

## Assessment of the two grounds

**GROUND 1 (tree control).** The computation is reproducible and the exact (H2)
failure is real (verified item 3). What needs repair is its *framing*: Proposition 1
is not violated, it does not apply (M5); the registered follow-up it claims to
discharge is about (H1), not (H2) (M6); the audit had already established that (H2)
failure permits forest growth (M6); the small-amplitude claim is numerically wrong
(C3); and the "forest" is an artefact of a truncation that discards `k - p` (m14).
After those repairs the ground still stands, in the weaker form: *for the exact
Leray blocks in full polarization, (H2) fails on a single edge, so Proposition 1
gives no control there, and the exhibited edge grows at parent `L^2` norm `0.131`.*

**GROUND 2 (no uniformity in `N`).** The `N^{-3/2}` decay and the strict decrease
are correct and are supported by the archive. The stated identity is printed with
the wrong constant (C1), the amplitude normalisation does not implement the stated
shell law (M7), the claimed independent confirmation is forced by construction
(M8), and every absolute threshold — including the headline `N ~ 24.6` — is
therefore not defensible as printed. The *direction* of the conclusion survives all
of these.

**Certificate.** Sound in construction, wrong in the last digit (M1), non-rigorous
in the final root (M2), unreproducible from the committed script (C2) and
uncheckable from the archive (M3).

---

**REFEREE: MAJOR-REVISIONS**
