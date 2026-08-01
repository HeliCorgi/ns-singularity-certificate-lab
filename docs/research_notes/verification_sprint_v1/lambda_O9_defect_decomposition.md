# Λ obligation O-9 — defect decomposition of the closable front bound, and adversarial minimisation

**Verification Sprint V1, workstream C.** No new ideas; this note verifies,
corrects or kills claims already on the record in
[`idea_monotone_dichotomy.md`](../ideas_2026_08_01/idea_monotone_dichotomy.md)
§B.2/§F/§G and
[`CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md`](../../candidates/CANDIDATE_SOLUTION_SPECTRAL_FRONT_DICHOTOMY.md)
§5–§6.

**Arithmetic labels are mandatory and are used everywhere below.**
`EXACT` = `fractions.Fraction`; `ENCLOSED` = rigorous rational interval around
an irrational (square-root) quantity; `FLOAT` = binary64 numpy diagnostic,
never a bound. Failed checks are recorded, not deleted.

Code: [`experiments/run_lambda_o9_defect_search.py`](../../../experiments/run_lambda_o9_defect_search.py),
`ns_certificate_lab.spectral_front_monotone.front_defect_decomposition`,
tests in `tests/test_front_defect_decomposition.py` (100 passing).
Data: [`outputs/verification_sprint_v1/lambda_o9/summary.json`](../../../outputs/verification_sprint_v1/lambda_o9/summary.json).

---

## 0. Verdict up front

| claim under audit | status |
|---|---|
| (I.1) exact front identity | **PASS** (EXACT, re-verified at r=0,1,2 on 7 fields) |
| (I.2)/(I.3) chain and nonnegativity of the gap | **PASS** (EXACT + ENCLOSED four-defect telescoping) |
| §F reported deficits 𝔡 = 0.865 / 0.740 | **CORRECTED** — those are ν-artefacts; the ν-minimised deficits of the same fields are **6/137 ≈ 0.0438** (in-support G) and **15/146 ≈ 0.1027** (full G) |
| `spectral_front_monotone` computes the doc's G_r | **CORRECTED** — the shipped `_ledger` restricts G_r to supp(u); the monotone's G₀ = ½‖ℙ(u·∇u)‖₂² is the *full* sum. Both are now available and labelled |
| **O-8** — uniform positive deficit 𝔡₀ > 0 over divergence-free fields | **KILLED**, by an EXACT one-parameter family with d ↓ 0; the kill is normalisation-proof (§5.3) |
| **O-9** via "energy neutrality forbids saturation" | **KILLED at r = 0** (the critical order): T₀ = 0 is *implied* by saturation, not violated by it. Survives only for r ≥ 1 (§8) |
| **O-9** via "leakage `G₀^out` forbids saturation" | **KILLED as a band-uniform mechanism** — the free search drives `G₀^out/G₀^full` to `4.0e-10` on the sup-norm-3 band; leakage weakens as the band widens (§8.3, §6) |
| **O-9** via "saturation forces spectral degeneration" | **ALIVE, qualitative only** — every near-saturating configuration found is asymptotically single-shell (`V₀/μ² → 0` at the same order as `d_*`). The conditional constant in `d_* ≳ c·min(1,V₀/μ²)` is neither supported nor refuted by this sprint's data (§8.4); smallest observed ratio `0.126` |

**VERDICT: CONDITIONAL** (see §11).

---

## 1. Setup

Repo conventions verbatim: `∂_t û_k = N_k − ν|k|²û_k`,
`N_k = −i P_k Σ_{ℓ+m=k}(m·û_ℓ)û_m`, `P_k = I − k⊗k/|k|²`, `T³`, zero mean,
`k·û_k = 0`. Paired-real ledger (`modal_front_actions`, `spectral_front_monotone`):

```
e_k = ½(|c_k|²+|s_k|²),  a_k = ½(c_k·c^N_k + s_k·s^N_k),  n_k = ½(|c^N_k|²+|s^N_k|²),
x_k = |k|²,  H_r = Σ x_k^r e_k,  T_r = Σ x_k^r a_k,  G_r = Σ x_k^r n_k,
μ = N_r² = H_{r+1}/H_r,  p_r(k) = x_k^r e_k/H_r,  V_r = Var_{p_r}(x),  g_k = a_k/e_k.
```

Two facts used throughout: `T_0 = Σ a_k = 0` (energy neutrality, repo F-12) and
pressure blindness. In complex convention every moment is the obvious one:
`H_r = Σ_k |k|^{2r}|û_k|²`, `T_r = Σ_k |k|^{2r} Re(û_k·conj(N̂_k))`,
`G_r = Σ_k |k|^{2r}|N̂_k|²`, summed over the full lattice.

---

## 2. Task 1 — every inequality as a defect identity

Write `C := Cov_{p_r}(x,g) = H_r^{-1} Σ_k x_k^r (x_k−μ) a_k` (EXACT rational) and
introduce the intermediate quantities

```
w_k := x_k^r |x_k−μ| / H_r        (weights, EXACT ≥ 0)
A   := Σ_k w_k |a_k|              (EXACT)
B   := Σ_k w_k sqrt(e_k n_k)      (irrational; ENCLOSED)
S   := sqrt(V_r G_r / H_r)        (irrational; ENCLOSED)
```

The published chain (I.1)→(I.2)→(I.3) is exactly the chain `C ≤ A ≤ B ≤ S`
followed by a square completion. Each link is a defect identity `RHS−LHS = Δ ≥ 0`:

### Δ_sign — triangle / sign defect  (EXACT rational)

```
Δ_sign := A − C = (2/H_r) Σ_{k : (x_k−μ)a_k < 0} x_k^r |x_k−μ| |a_k|  ≥ 0.
```
It is *twice the misaligned weighted mass*: it vanishes iff every mode with
`x_k > μ` gains energy from the nonlinearity and every mode with `x_k < μ`
loses it (forward-cascade sign pattern).

### Δ_CS,modal — aggregated per-mode Cauchy–Schwarz defect  (radical-free certificate)

```
Δ_CS,modal := B − A = Σ_k w_k ( sqrt(e_k n_k) − |a_k| )
                    = Σ_k w_k · δ_k / ( sqrt(e_k n_k) + |a_k| ),
             δ_k := e_k n_k − a_k²  ≥ 0     (EXACT, per mode).
```
`δ_k ≥ 0` is the radical-free certificate; it is the paired-vector
Cauchy–Schwarz `⟨(c_k,s_k),(c^N_k,s^N_k)⟩² ≤ |·|²|·|²` and is asserted mode by
mode in EXACT arithmetic. `δ_k = 0` iff `(c^N_k, s^N_k) = λ_k (c_k, s_k)`.

### Δ_CS,vector — the vector Cauchy–Schwarz gap  (square certificate)

With `ξ_k := x_k^{r/2}|x_k−μ| sqrt(e_k)` and `η_k := x_k^{r/2} sqrt(n_k)`,
`|ξ|² = H_r V_r`, `|η|² = G_r`, `⟨ξ,η⟩ = H_r B`, so

```
Δ_CS,vector := S − B = (1/H_r)( sqrt(|ξ|²|η|²) − ⟨ξ,η⟩ ),
   |ξ|²|η|² − ⟨ξ,η⟩² = ½ Σ_{j,k} (ξ_jη_k − ξ_kη_j)²  ≥ 0   (Lagrange/Gram).
```

The Gram determinant is *not* rational (it contains `B`). Two radical-free
substitutes are recorded, and the second is what the code certifies:

1. **AM–GM linearisation.** For every rational `t>0`,
   `|x−μ| sqrt(e_k n_k) ≤ ½( t(x−μ)²e_k + n_k/t )`, hence
   `B ≤ ½( tV_r + G_r/(tH_r) )`, with the optimum at `t = sqrt(G_r/(H_rV_r))`
   returning exactly `S`. **The composite step `Δ_CS,vector + Δ_SC` is precisely
   this AM–GM at the rational value `t = 2ν`** — see (T) below.
2. **Rational combined certificate.** `A` is rational and `S² = V_rG_r/H_r` is
   rational, so
   ```
   S² − A² = V_r G_r/H_r − A²  ≥ 0      (EXACT)
   ```
   is a fully radical-free certificate for the *combined* step `A ≤ S`, i.e. for
   `Δ_CS,modal + Δ_CS,vector = (S²−A²)/(S+A)`.

### Δ_SC — square completion  (ENCLOSED)

```
Δ_SC := (2ν/μ) [ sqrt(V_r) − sqrt(G_r/H_r)/(2ν) ]²  ≥ 0.
```

### The exact telescoping (T)

Normalise the three covariance-chain defects to rate units by the factor `2/μ`:

```
    G_r/(2νH_{r+1}) − d/dt log N_r²
  = (2/μ)Δ_sign + (2/μ)Δ_CS,modal + (2/μ)Δ_CS,vector + Δ_SC .        (T)
```

*Proof.* `d/dt log N_r² = 2C/μ − 2νV_r/μ` is (I.1). Expanding
`Δ_SC = (2ν/μ)V_r − (2/μ)S + G_r/(2νμH_r)` and using `μH_r = H_{r+1}` gives
`(2/μ)(S−C) + Δ_SC = G_r/(2νH_{r+1}) − d/dt log N_r²`; the middle telescope
`S−C = (A−C)+(B−A)+(S−B)` is the chain. ∎

**Two halves of (T) are separately checkable, and both are checked.**

* *Rational half* (EXACT `Fraction` equality, asserted in code):
  ```
  gap_total = Δ̂_sign + R,   Δ̂_sign = (2/μ)(A−C),
  R := (2/μ)[ νV_r + G_r/(4νH_r) − A ]  =  Δ̂_CS,modal + Δ̂_CS,vector + Δ_SC ,
  ```
  because the two square roots `B` and `S` cancel between the last three terms
  — `R` is the AM–GM step at `t = 2ν`. Both `Δ̂_sign ≥ 0` and `R ≥ 0` are
  asserted exactly.
* *Radical half* (ENCLOSED): each of the last three defects is evaluated with
  96-bit rigorous rational `sqrt` enclosures and their sum is asserted to
  contain `R`.

**Verification.** `front_defect_decomposition` runs (T) on the relay triads
(s = 1,2,3 and child 1/32), the P1/P2/P3 families, and every rationalised
near-minimiser of §6, at r ∈ {0,1,2}, ν ∈ {1/40, 1/10, 1}, both G conventions.
All exact equalities hold; all defects are ≥ 0; the enclosure widths are
≲ 10⁻²⁸. `front_gap_identity`'s published `gap_total` is reproduced exactly in
the `in_support` convention.

---

## 3. The exact ν-reduction (this is the load-bearing simplification)

For a fixed field, the deficit
`d(u,ν) := gap_total / (G_r/(2νH_{r+1}))` is a rational function of ν alone:

```
d(u,ν) = 1 − (4H_{r+1}/(μ G_r))(νC − ν²V_r).
```

Maximising `νC − ν²V_r` gives, **EXACTLY**,

```
ν_* = C/(2V_r) ,      d_*(u) := min_{ν>0} d(u,ν) = 1 − C² H_r /( V_r G_r ).   (R)
```

`d_*` has no ν in it. Because `e_k ∝ α²`, `a_k ∝ α³`, `n_k ∝ α⁴` under `u ↦ αu`,
and `(H_r,T_r,G_r) ↦ (λ^{2r+2}, λ^{2r+4}, λ^{2r+6})·` under `u ↦ λu(λx)`,

> **`d_*` is invariant under `u ↦ αu` (with ν ↦ αν) and under `u ↦ λu(λx)`
> (ν fixed); and `d_*(u) ∈ [0,1]` with `d_* = 0` iff `C = S`, i.e. iff all three
> covariance-chain defects vanish.**

(Orientation: `N(−u) = N(u)` but `a_k(−u) = −a_k(u)`, so `C` flips sign under
`u ↦ −u`. `ν_* > 0` needs `C > 0`, so exactly one of the two orientations is
admissible; for the other, `d_* = 1` by convention because `d(u,ν) > 1` for all
`ν > 0`. This is a real trap in the numerics — see §6.)

Two consequences, both used below:

* The square-completion defect `Δ_SC` is *not* an independent obstruction: it
  can always be annihilated by choosing ν (equivalently, by rescaling the
  amplitude). The whole content of O-8/O-9 lies in `Δ_sign + Δ_CS,modal +
  Δ_CS,vector`.
* Since `d_*` is *constant along the orbits of the critical scaling group*,
  **no "critical normalisation" can change the range of `d_*`.** Any
  normalisation (fix ν = 1 and ‖u‖₂ = 1, or fix K, or fix a Reynolds number)
  is a section of that group action; every orbit meets it; the infimum is
  unchanged. This kills the escape route sketched in §G-O-8
  ("uniform lemma … at critical normalisation").

*Correction to §F/§4 of the candidate.* The reported deficits
𝔡 = 0.865 (s=1), 0.740 (s=2) were computed at ν = 1/40 and are dominated by
`Δ_SC`, i.e. they measure how far ν = 1/40 is from ν_*, not how far the field
is from saturation. For the same relay triad, `ν_* = 131/384` and
`d_* = 6/137 = 0.04380` (in-support G) / `15/146 = 0.10274` (full G) — EXACT.

---

## 4. Correction: which `G_r`?

`_ledger` iterates over the *field's* coefficient table, so the shipped
`G_r` is `G_r^{in} = Σ_{k: e_k>0} x_k^r n_k`. The published objects
`G_0 = ½‖ℙ(u·∇u)‖₂²`, `K = G_0/(2H_1²)` and Lemma K
(`K ≤ S_N`) all refer to the **full** sum `G_r^{full} = Σ_{all k} x_k^r n_k`.

Both are legitimate in (I.2) — `a_k = 0` off supp(u), so Cauchy–Schwarz may be
run on supp(u) only, and `G^{in}` gives the *sharper* true inequality — but
they are not interchangeable in the monotone: only `G^{full}` is the quantity
Lemma K bounds and the quantity in `KD ≤ ‖u‖_∞²`. The difference

```
G_r^{out} := G_r^{full} − G_r^{in}   ("leakage": nonlinear power deposited on
                                      modes the field does not occupy)
```

is a lower bound for the full-convention deficit:
`d_*^{full} ≥ G_0^{out}/G_0^{full}`. `full_nonlinear_power` and the
`convention` argument of `front_defect_decomposition` now make the choice
explicit; both are reported for every field. **Everything in §5–§9 is reported
in both conventions.**

---

## 5. Task 2 — the target rigidity lemma, and its refutation

> **Rigidity lemma (O-8), as it must be stated.** There exists `c > 0` such
> that for every divergence-free zero-mean band-limited `u` on `T³` and every
> `ν > 0`,
> ```
> Δ̂_sign + Δ̂_CS,modal + Δ̂_CS,vector + Δ_SC  ≥  c · G_0/(2νH_1),
> ```
> equivalently `d_*(u) ≥ c`, equivalently
> `sup_u C²H_0/(V_0G_0) ≤ 1−c`, equivalently
> `d/dt log N_0² ≤ (1−c)·G_0/(2νH_1)`.

### 5.1 Two-mode fields: `d_* = 1` (EXACT, complete classification)

Let `supp(u) = {±k₁, ±k₂}`.
*Collinear* (`k₂ = m k₁`): then `u = f(k₁·x)` with `k₁·f ≡ 0`, so
`(u·∇)u ≡ 0` and every moment vanishes.
*Non-collinear*: the sums `ℓ+m` with `ℓ,m ∈ supp(u)` land only on
`{0, ±2k₁, ±2k₂, ±(k₁+k₂), ±(k₁−k₂)}`, and non-collinearity excludes every one
of these from being `±k₁` or `±k₂` (`k₁±k₂ = ±k₁` or `2k₁ = ±k₂` would force
collinearity). So `N_{k} = 0` for `k ∈ supp(u)` — no self-interaction is even
needed, though `N_{2k} = −i(k·û_k)P_{2k}û_k = 0` identically by
divergence-freeness as well. Hence `a_k ≡ 0`, `C = T₁ = 0`, and `d_* = 1` — the
*maximal* deficit.
Verified EXACT on three examples in `summary.json ▸ two_mode_classification`
(`T_1 = 0`, `G_0^{in} = 0`, `G_0^{full} > 0`).

**A closing triad `k₁+k₂=k₃` is the minimal structure with a nonzero deficit
gradient.**

### 5.2 Three-mode Leray relay: closed form, and `d_* ↓ 0` (EXACT)

Take the repo's exact relay triad (`exact_leray_relay.build_exact_relay_triad`)
with general amplitudes, `p = s(1,1,0)`, `q = s(1,0,1)`, `c = p+q`, `n=(1,−1,−1)`:

```
u = B e₃ sin(p·x) + C e₂ cos(q·x) + D n cos(c·x).
```

Closed form (derived by hand from the module's exact transfer formulas, then
re-verified EXACTLY against the ledger — every `(s, D, convention)` point of
`stage_closed_form` carries a `*_matches_ledger` boolean, and
`test_relay_family_closed_form_matches_the_exact_ledger` asserts the identity
at `s ∈ {1,2,3}`, `D ∈ {1, 1/8, 1/128}`, both conventions):

```
H₀ = (B²+C²+3D²)/2,   H₁ = s²(B²+C²+9D²),   H₂ = 2s⁴(B²+C²+27D²),
T₀ = 0,               T₁ = 2s³BCD,
n_c = s²B²C²/6,  n_p = 3s²C²D²/8,  n_q = 3s²B²D²/8,  n_{p−q} = 0 (exact cancellation),
n_{p+c} = 3s²B²D²/8,  n_{q+c} = 3s²C²D²/8,  n_{2p}=n_{2q}=n_{2c}=0 (divergence-freeness).
```

With `P := B²+C²`, `H₀V₀ = 24s⁴PD²/(P+3D²)`, and therefore, **independent of
`s`** (the r = 0 bound is exactly scale-invariant):

```
d_*^{in}   = 1 − B²C²(P+3D²) / ( P[ B²C² + (9/4)D²P ] )
d_*^{full} = 1 − B²C²(P+3D²) / ( P[ B²C² + (9/2)D²P ] )
V₀/μ²      = H₀H₂/H₁² − 1 = 12 P D² / (P+9D²)² .
```

At `B = C = 1` these collapse to
`d_*^{in} = 6D²/(2+9D²)`, `d_*^{full} = 15D²/(2+18D²)`, `V₀/μ² = 24D²/(2+9D²)²`.

| D | `d_*^{in}` EXACT | `d_*^{full}` EXACT | `V₀/μ²` EXACT | `d_*^{full}/(V₀/μ²)` | `d_*^{in}/(V₀/μ²)` |
|---|---|---|---|---|---|
| 1 | 6/11 = 0.545455 | 3/4 = 0.75 | 24/121 = 0.198347 | 3.78125 | 2.75000 |
| 1/2 | 6/17 = 0.352941 | 15/26 = 0.576923 | 96/289 = 0.332180 | 1.73678 | 1.06250 |
| 1/8 | **6/137 = 0.0437956** | **15/146 = 0.1027397** | 1536/18769 = 0.0818371 | 1.25542 | 0.53516 |
| 1/32 | 6/2057 = 0.00291687 | 15/2066 = 0.00726041 | 24576/4231249 = 0.00580821 | 1.25002 | 0.50220 |
| 1/128 | 6/32777 = 1.83055e-4 | 15/32786 = 4.57512e-4 | 393216/1074331729 = 3.66010e-4 | 1.25000 | 0.50014 |
| 1/512 | 6/524297 = 1.14439e-5 | 15/524306 = 2.86092e-5 | 6291456/274887344209 = 2.28874e-5 | 1.25000 | 0.50001 |

**`d_*(u) → 0` along an EXACT one-parameter family of genuine divergence-free
zero-mean rational fields.** Asymptotically `d_*^{full} = 15D²/P + O(D⁴)`,
`d_*^{in} = 6D²/P + O(D⁴)`, `V₀/μ² = 12D²/P + O(D⁴)`, so

```
d_*^{full}/(V₀/μ²) → 5/4 = 1.25 ,     d_*^{in}/(V₀/μ²) → 1/2 .
```

Each row is a **certified EXACT upper bound on the infimum**, attained at the
exact rational viscosity `ν_* = C/(2V₀)` (e.g. `ν_* = 131/384` at `D = 1/8`).

### 5.3 Consequence

**The rigidity lemma of O-8 is FALSE**, in both G conventions, already on the
band `|k|_∞ ≤ 2`; and by §3 it stays false under every critical normalisation.
Correspondingly the promised sharpening
`d/dt log N₀² ≤ (1−𝔡₀)·KD/(2ν)` with a universal `𝔡₀ > 0` **does not exist**.

*Scope discipline.* Kill condition **K3** is worded as "`𝔡 → 0` along genuine
`N`-sweeps of **evolving** fields". That is a trajectory statement and is **not**
what is settled here. What is settled is the *field-level* obligation O-8, of
which K3 was the numerical proxy: over the class of all divergence-free
zero-mean band-limited fields the infimum is 0, as a theorem, not a trend.
Whether a Navier–Stokes trajectory ever visits such near-degenerate
configurations is untouched by this note.

---

## 6. Task 3 — adversarial search protocol and results (FLOAT, then EXACT)

**Protocol.** Fields are parameterised over all canonical `k` with
`1 ≤ |k|_∞ ≤ B`: four real parameters per mode, the components of the cosine
and sine coefficient vectors in the *integer* orthogonal basis
`t₁ = k × e_j` (first `j` with `k × e_j ≠ 0`), `t₂ = k × t₁` of `k^⊥`.
Divergence-freeness is therefore structural, and rational parameters give an
exactly divergence-free rational field. **`n_k` is never a free variable**: `N`
is recomputed from the field by an alias-free FFT convolution
(`grid = 4B+1 ≥ 2·(2B)+1`) followed by the exact Leray symbol, so the search is
over fields only, exactly as required.

Objective: maximise `J := T₁²/((H₂−H₁²/H₀)·G₀) = 1 − d_*` (eq. (R)), which is
ν-free and scale-free. **Analytic gradient** (adjoint of the bilinear form:
`grad_u ⟨N(u),ψ⟩ = −(∂_m u_j)ψ_j + (u·∇)ψ_m` for divergence-free `ψ`;
`grad_u (H₂−H₁²/H₀) = 2(x−μ)²û`), validated against central differences to
`rel 1e-4` in the test suite. Projected ascent on the unit sphere with
backtracking line search; multi-start = random Gaussian starts + structured
relay-triad seeds at four child amplitudes + the previous band's minimiser
embedded. ν ∈ {1/40, 1/10, 1} is *not* swept: eq. (R) already minimises over it
exactly (all three appear in the exact telescoping tables instead).

*Convention caveat for the float lane.* In binary64 every band mode is
generically occupied, so the float `in_band` column uses the mask
`|k|_∞ ≤ B` as a stand-in for `supp(u)`; it coincides with the exact
`in_support` moment only when the minimiser occupies the whole band. All exact
numbers below use the true support.

**Results (FLOAT, binary64 — not bounds).** See `summary.json ▸
adversarial_search`.

| band B (sup-norm cutoff) | canonical modes | parameters | `d_*` min, full G | `d_*` min, in-band G | minimiser `V₀/μ²` | minimiser leakage `G₀^out/G₀^full` | `T₀/T₁` residual |
|---|---|---|---|---|---|---|---|
| 1 | 13 | 52 | 2.7635e-03 | 1.5516e-04 | 3.7938e-03 | 2.2818e-03 | 4.3e-16 |
| 2 | 62 | 248 | 5.4230e-04 | 5.4328e-04 | 4.2931e-03 | 2.8204e-07 | 1.1e-17 |
| 3 | 171 | 684 | 1.1197e-03 | 9.9344e-04 | 5.0791e-03 | 3.9712e-10 | 3.1e-16 |

Three things in this table matter more than the minima themselves.

* **The leakage fraction at the minimiser collapses with band size**:
  `2.3e-3 → 2.8e-7 → 4.0e-10`. The optimiser spends the extra modes on
  cancelling `G₀^{out}` — i.e. **zero-leakage near-saturation is achievable**,
  and §8.3's gatekeeper is defeated in practice, not merely in the special
  relay geometry. Correspondingly `d_*^{full}` and `d_*^{in-band}` converge.
* **The spread collapses with the deficit** (`V₀/μ² ≈ 4–5e-3` at every
  minimiser): §8.4.
* `T₀ = 0` holds to `~1e-16` relative at every minimiser — energy neutrality is
  never violated by the near-saturating configurations, in agreement with §8.2.
* The `B = 3` minimum is *worse* than the `B = 2` minimum. Bands are nested, so
  this is optimiser failure, not a property of the problem, and is recorded as
  such — the `B = 3` landscape is 684-dimensional and the near-minimiser is
  near-degenerate.

**Certified EXACT points.** Float minimisers are sparsified to their strongest
modes, dyadically rationalised (common denominator `2^b`, `b ∈ {14,16,20}`) and
re-evaluated with `front_defect_decomposition`. Only those exact numbers are
upper bounds. Because the float objective `T₁²/(WG)` is sign-blind while
`d_*` needs `Cov > 0`, the admissible orientation (`u` or `−u`, which have the
same `N` and opposite `a_k`) is selected before certification; the choice is
recorded in the `orientation` field. *Recorded failure:* an earlier pass that
omitted this step returned `d_* = 1` for every rationalised minimiser — the
sign, not the rounding, was the problem. See `summary.json ▸
exact_certificates`.

All entries are EXACT rationals; the numerator/denominator digit counts are
given for the unwieldy ones (full values in `summary.json`). Every one is
attained at the exact rational `ν_*` also recorded there.

| label | exact `d_*` (full G), as a float | numerator/denominator size | `V₀/μ²` | leakage `G₀^out/G₀^full` |
|---|---|---|---|---|
| `band1_keep8_dyadic2^14` | 7.260117e-03 | 30/32 digits | 3.7783e-03 | 6.7928e-03 |
| `band1_keep13_dyadic2^16` | 2.764261e-03 | 37/40 digits | 3.7946e-03 | 2.2823e-03 |
| `band1_keep13_dyadic2^20` | **2.763470e-03** | 46/49 digits | 3.7938e-03 | 2.2818e-03 |
| `band2_keep24_dyadic2^16` | 5.917913e-04 | 40/43 digits | 4.2923e-03 | 4.0089e-05 |
| `band2_keep62_dyadic2^16` | 5.677232e-04 | 40/43 digits | 4.2925e-03 | 4.2249e-06 |
| `band2_keep62_dyadic2^20` | **5.425815e-04** | 51/55 digits | 4.2931e-03 | 3.5392e-07 |
| `band3_keep24_dyadic2^16` | 1.185587e-03 | 37/40 digits | 5.0797e-03 | 6.2972e-05 |
| `band3_keep64_dyadic2^20` | **1.118969e-03** | 59/62 digits | 5.0790e-03 | 2.6250e-07 |
| `relay_triad_D=1/8` (reference) | 15/146 = 1.027397e-01 | — | 8.1837e-02 | 6.1644e-02 |
| `relay_triad_D=1/128` (reference) | 15/32786 = 4.575124e-04 | — | 3.6601e-04 | 2.7451e-04 |

Each bold row reproduces its float minimiser to 4 significant figures, so the
rationalisation is faithful once the orientation is fixed. **These are certified
EXACT upper bounds:**
`inf d_*^{full} ≤ 5.4258e-4` on `|k|_∞ ≤ 2` and `≤ 1.1190e-3` on `|k|_∞ ≤ 3`,
each with an explicit rational field and an explicit rational `ν_*`.

The best certified exact upper bounds on `inf d_*` nevertheless remain the
relay-family values of §5.2 (down to `15/524306 = 2.86e-5`), which are also the
only ones available in closed form as a *sequence* tending to zero.

---

## 7. Task 4 — equality conditions extracted at near-minimisers

Saturation of the covariance chain at order `r` requires, simultaneously,

1. `δ_k = e_k n_k − a_k² = 0` for every occupied mode — i.e. `N_k ∥ û_k` in the
   paired-vector sense, `N_k = λ_k û_k`;
2. `n_k = γ²(x_k−μ)² e_k` for one constant `γ ≥ 0` (vector-CS equality) — in
   the **full** convention this also forces `n_k = 0` for every `k ∉ supp(u)`
   (the vector `ξ` vanishes there, so `η` must too): **zero leakage**;
3. `sgn(a_k) = sgn(x_k−μ)` (sign defect zero).

Combining (1)–(3): **`a_k = γ(x_k−μ)e_k` and `N_k = γ(x_k−μ)û_k`.** Adding the
square completion at `ν = ν_*` fixes `γ = 2ν`, so exact saturation of (I.3) at
`r = 0` is precisely the nonlinear eigen-equation

```
        ℙ((u·∇)u) = 2ν (Δ + μ) u ,      μ = H₁/H₀ = ‖∇u‖²/‖u‖² .        (E)
```

Extracted at the near-minimisers (`summary.json ▸ equality_conditions`), for
the relay triad at `D = 1/8` and `D = 1/128`:

* **Sign defect: exactly zero.** `Δ_sign = 0` at every ν and order for the
  relay family; the aligned weight fraction is 1. The two parent modes have
  `x_p − μ < 0` and `a_p < 0`; the child has `x_c − μ > 0` and `a_c > 0`. **The
  alignment sign flips exactly at μ**, as the rigidity condition demands, and
  μ sits infinitesimally above the parent shell (`x_p − μ = −12s²D²/P → 0⁻`).
* **Per-mode CS defect: exactly zero on the child, large but weightless on the
  parents.** At `D = 1/8`: `e_c n_c = a_c² = 1/256`, so `δ_child = 0` — `N_c` is
  exactly parallel to `û_c`. On each parent `e n = 3/1024`, `a² = 1/1024`, so the
  *relative* defect is `δ/(en) = 2/3` — the modal Cauchy–Schwarz is badly
  violated there. It does not matter, because the parents' weight
  `w_k ∝ |x_k − μ| = 12s²D²/P` is `O(D²)`: **saturation is achieved by making the
  badly aligned modes weightless, not by aligning them.**
* **Vector-CS ratio `n_k/((x_k−μ)²e_k)`:** on the child it is exactly
  `γ² = (BC/(12sD))²`; on the parents it is `C²P²/(192s²D²B²)`, i.e. `3γ²` at
  `B = C` — off by a *fixed factor 3*, never converging. It does not matter for
  the same reason as above: the parents contribute `72s⁴D⁴/P` to
  `|ξ|² = H₀V₀ = 24s⁴PD²/(P+3D²)` and `(3/8)s²D²P` to `G₀ ≈ s²B²C²/6`, both
  `O(D²)` *relative*, so the vector Cauchy–Schwarz is dominated by the child
  alone and is nearly saturated by a single term. **The residual deficit lives
  entirely in the parent modes and the two leakage modes `p+c`, `q+c`, all at
  relative order `D²`.**
* **Where the deficit sits.** At `ν = ν_* = 131/384`, `D = 1/8`, `r = 0`
  (EXACT / ENCLOSED): `Δ̂_sign = 0`, `Δ̂_CS,modal = 3.9158e-3`,
  `Δ̂_CS,vector = 9.0946e-3`, `Δ_SC = 0`, `gap = 1.3010e-2`,
  `G₀/(2ν_*H₁) = 1.2664e-1`, `d_* = 15/146`. In the `in_support` convention
  `Δ̂_CS,vector` drops to `1.3734e-3` and `d_* = 6/137`; the difference is
  exactly the leakage of §8.3.
* **Square completion `V₀ = G₀/(4ν²H₀)`:** holds identically at `ν = ν_*` by
  construction (eq. (R)); `ν_*` for the family is `ν_* = C/(2V₀) ∝ 1/D → ∞`,
  i.e. at fixed ν the saturating amplitude tends to zero.

**Emerging structure (this is the answer to task 4).** The near-saturating
configuration is *not* a cascade. It is

> a **two-mode core** carrying almost all the energy in a single shell, plus a
> **vanishing-amplitude tracer** at the sum wavevector carrying almost all the
> nonlinear power, with μ pinned infinitesimally above the core shell so that
> `x_k − μ` is `O(D²)` on the core and `O(1)` on the tracer.

Every saturation condition is met *because* the field degenerates, not because
it cascades.

---

## 8. Task 5 — compatibility with divergence-freeness, energy neutrality, shell concentration

### 8.1 Divergence-freeness: compatible, and it *helps*

`N_{2k} = −i(k·û_k)P_{2k}û_k = 0` identically for divergence-free fields, so the
strongest single-mode leakage channel (the self-interaction of the extreme
mode) is closed for free. Divergence-freeness is imposed structurally in the
search and exactly in every certified field. **No obstruction.**

### 8.2 Energy neutrality `T₀ = Σ a_k = 0` — the central O-9 question

Substituting the saturation profile `a_k = γ(x_k−μ_r)e_k` into `T₀`:

```
T₀ = γ Σ_k (x_k − μ_r) e_k = γ ( H₁ − μ_r H₀ ) = γ H₀ ( N₀² − N_r² ).       (N)
```

* **At `r = 0`, `μ_0 = H₁/H₀`, so `T₀ = 0` identically.** Energy neutrality is
  *automatically satisfied* by the saturating profile — it is a *consequence*
  of saturation, not a constraint on it. **Energy neutrality imposes no
  positive deficit floor at r = 0.** This is the precise answer to the O-9
  question, and it is negative.
* **At `r ≥ 1` it does obstruct.** `r ↦ N_r² = H_{r+1}/H_r` is non-decreasing
  (Cauchy–Schwarz `H_{r+1}² ≤ H_r H_{r+2}`), so (N) forces either `γ = 0` — in
  which case `n_k = γ²(x_k−μ)²e_k = 0` on supp(u), `G_r^{in} = 0`, and the
  bound is vacuous — or `N_r² = N₀²`, which forces equality in every
  Cauchy–Schwarz step, i.e. `x_k ≡ μ` on supp(u), i.e. a single-shell field, for
  which `V_r = C = 0` and both sides of (I.3) degenerate. **So for `r ≥ 1`
  energy neutrality does forbid non-degenerate exact saturation.** But `r ≥ 1`
  is not scale-critical (§C of the derivation) and cannot feed a Clay
  statement, so this survival is of no use to O-9.

Numerically, `T₀ = 0` is confirmed EXACT on every field tested and to
`|T₀|/|T₁| < 1e-12` (FLOAT) at every search minimiser.

### 8.3 Leakage — the only remaining structural obstruction

`d_*^{full} ≥ G₀^{out}/G₀^{full}`. Leakage is *not* generic-zero: for a triad
`{p,q,p+q}` the difference mode `p−q` receives `O(BC)` nonlinear power — the
same order as the child — unless it cancels exactly. In the repo's relay triad
it *does* cancel (`n_{p−q} = 0`, EXACT), which is exactly why that family can
reach `d_* → 0`; the surviving leakage `n_{p+c}, n_{q+c} = O(D²)` is of the same
order as the deficit. For a triad without the cancellation the leakage fraction
is `O(1)` and so is the deficit. **Leakage cancellation, not energy neutrality,
is the real gatekeeper of near-saturation** — and it is achievable: the free
search drives `G₀^{out}/G₀^{full}` to `2.8e-7` at `|k|_∞ ≤ 2` and `4.0e-10` at
`|k|_∞ ≤ 3` (§6), i.e. the extra modes are spent precisely on closing the
leakage channels. So leakage is a *finite-band* obstruction that grows weaker,
not stronger, as the band widens — the opposite of what O-9 needs.

### 8.4 Shell concentration — the surviving (but weakened) positive statement

**The structural fact is solid and EXACT.** Along the relay family
`V₀/μ² = 12D²/P + O(D⁴) → 0`: the deficit and the relative spectral spread
vanish together and at the same order. Saturation of (I.3) at the critical
order is achievable only by fields whose spectral measure collapses onto a
single shell; every near-minimiser found by the float search shows the same
collapse (the `minimiser_spread` column of §6 tracks `d_*` down). The natural
conditional lemma is therefore

```
        d_*(u)  ≥  c · min(1, V₀/μ²)                                       (C)
```

with, in the exact relay family, `c = 5/4` (full G) and `c = 1/2` (in-support
G) achieved in the small-`D` limit. The constrained float search (minimise
`d_*` at fixed `V₀/μ²`; `summary.json ▸ constrained_search`) gives:

| band | target `V₀/μ²` | achieved `V₀/μ²` | `d_*` full | ratio `d_*/(V₀/μ²)` |
|---|---|---|---|---|
| 2 | 0.003 | 2.9986e-03 | 3.6897e-03 | 1.2305 |
| 2 | 0.01 | 9.9854e-03 | 1.1991e-02 | 1.2009 |
| 2 | 0.03 | 2.9880e-02 | 3.3387e-02 | 1.1174 |
| 2 | 0.1 | 9.9117e-02 | 7.1035e-02 | 0.7167 |
| 2 | 0.3 | 2.9674e-01 | 5.6684e-02 | 0.1910 |
| 2 | 1 | 9.9455e-01 | 1.8222e-01 | 0.1832 |
| 2 | 3 | 2.5208e+00 | 1.0000e+00 | 0.3967 |
| 3 | 0.003 | 2.9996e-03 | 1.5763e-03 | 0.5255 |
| 3 | 0.01 | 9.9951e-03 | 4.9905e-03 | 0.4993 |
| 3 | 0.03 | 2.9956e-02 | 1.4606e-02 | 0.4876 |
| 3 | 0.1 | 9.9553e-02 | 4.1849e-02 | 0.4204 |
| 3 | 0.3 | 2.9720e-01 | 5.4284e-02 | 0.1827 |
| 3 | 1 | 9.5637e-01 | 2.6913e-01 | 0.2814 |
| 3 | 3 | 2.8872e+00 | 4.6840e-01 | 0.1622 |

**Recorded failure — this table is not trustworthy as a lower envelope.** The
*unconstrained* `B = 2` minimiser of §6 has `V₀/μ² = 4.2931e-3` and
`d_* = 5.4230e-4`, i.e. ratio `0.126`, whereas the penalised descent at the
neighbouring target `v = 0.003` returned `1.23`. A penalised descent that lands
an order of magnitude above a point the free search already found is
optimiser failure, not a measurement, so the apparent "`1.20` at `B = 2` vs
`0.49` at `B = 3`" band trend must **not** be read as evidence. The only
defensible numbers are the free minimisers themselves:

| band | `V₀/μ²` at free minimiser | `d_*` | ratio |
|---|---|---|---|
| 1 | 3.7938e-03 | 2.7635e-03 | 0.728 |
| 2 | 4.2931e-03 | 5.4230e-04 | **0.126** |
| 3 | 5.0791e-03 | 1.1197e-03 | 0.220 |

and the exact relay family (ratio `5/4` for full G, `1/2` for in-support G).
**Smallest ratio observed anywhere: `0.126`, at band 2.** It is above zero, but
the sweep is far too coarse and the `B = 3` optimisation demonstrably
under-converged, so *no positive band-uniform constant is established, and none
is refuted either.*

What *is* solid is the qualitative statement: the bound (I.3) is asymptotically
saturable **only** by spectrally collapsing fields (`V₀/μ² → 0` at every
minimiser and along the exact family), whereas a blow-up front in the
`γ ∈ [½,1)` window is required by the log-critical shell law to keep
`V₀/μ² = Θ(1)`. That is the honest residue of O-9:

> **Conjecture O-8′ (open; measured evidence weak).** There is `c > 0` with
> `d/dt log N₀² ≤ (1 − c·min(1, V₀/μ²))·G₀/(2νH₁)` for all divergence-free
> zero-mean band-limited fields. It is the only surviving quantitative
> descendant of O-8. Smallest observed value of `d_*/min(1, V₀/μ²)`: `0.126`.
> **Status: not supported and not refuted; the search protocol needs a proper
> constrained optimiser before it can decide this.**

---

## 9. Task 6 — band dependence and the demotion rate

No `B^{-β}` rate exists, because the infimum is already **exactly 0 at
`B = 2`** (the relay triad lives in `|k|_∞ ≤ 2`) and bands are nested, so

```
   inf_{|k|_∞ ≤ B} d_*  =  0    for every  B ≥ 2,   i.e.  β is undefined /
                                the correct statement is d_min(B) ≡ 0, not a rate.
```

For `B = 1` there is no exact *sequence* witness — the closing triads available
in `|k|_∞ ≤ 1` (e.g. `(1,0,0)+(0,1,0)=(1,1,0)`) all have a difference mode
`(1,−1,0)` *inside* the band, so §8.3's gatekeeper must be defeated by
cancellation among the 13 available modes rather than by the relay geometry.
The search does defeat it: the certified EXACT point
`band1_keep13_dyadic2^20` gives `d_* = 2.7635e-3` on `|k|_∞ ≤ 1`, with residual
leakage `2.3e-3`. So `inf_{B=1} d_* ≤ 2.76e-3`; **whether it is 0 is open, and
is recorded as such rather than asserted.**

The other side of the same coin is the leakage collapse of §6: at `B = 3` the
minimiser has `G₀^{out}/G₀^{full} = 4.0e-10`. **Leakage is not a band-uniform
obstruction either.**

**Recommendation.** Demote O-8 from "open obligation" to **REFUTED**, and record
the rate as `d_min(B) ≡ 0 for B ≥ 2` (no power law). The candidate's §4 line
"𝔡 decreasing in scale is the open K3 question" should be rewritten: at field
level the question is closed (𝔡 has no floor); the surviving open question is
the trajectory-level one, and it should be renumbered so it is not confused
with O-8.

---

## 10. Task 7 — the ε-regularised identity for `log(N₀² + ε)`

Let `N_{0,ε}² := (H₁ + εH₀)/H₀ = N₀² + ε` and `μ_ε := μ + ε`. Using
`Ḣ₁ = 2(T₁ − νH₂)` and — this is where `T₀ = 0` enters — `Ḣ₀ = −2νH₁`:

```
½ d/dt log(N₀²+ε) = [ (T₁ − νH₂ − ενH₁)/(H₀μ_ε) ] + νμ
                  = (1/μ_ε)[ T₁/H₀ − ν(H₂/H₀ − μ²) ]
                  = ( C − ν V₀ ) / (μ + ε) .                             (I.1_ε)
```

**The numerator is identical to (I.1); only `μ` is replaced by `μ+ε`.** Hence
the entire chain of §2 goes through verbatim, and the ε-robust analogue of
(I.3) is

```
      d/dt log(N₀² + ε)  ≤  G₀ / ( 2ν (H₁ + εH₀) )                       (I.3_ε)
```

with the four-defect telescoping (T) holding with `μ ↦ μ+ε` throughout, so

* every defect is multiplied by `μ/(μ+ε)`; and
* **the relative deficit is exactly ε-invariant: `d_ε(u,ν) = d(u,ν)`.**
  In particular §5's refutation of O-8 is unaffected by regularisation.

**Degradation factor.** `d/dt log N₀² = (1 + ε/N₀²)·d/dt log(N₀²+ε)`, and
`(I.3_ε)` RHS `= (I.3)` RHS `/(1 + εH₀/H₁)`. So

```
   d/dt log N₀²  ≤  (1 + εH₀/H₁) · G₀/(2ν(H₁+εH₀))  =  G₀/(2νH₁)   exactly.
```

The regularisation is **lossless**: the bound degrades by *exactly* the factor
`1 + εH₀/H₁ = 1 + ε/N₀²`, never more. For a uniform degradation `≤ 1+δ` one
needs `ε ≤ δ · inf_t N₀²(t)`; and on `T³` the spectral gap `|k| ≥ 1` gives
`N₀² = H₁/H₀ ≥ 1`, so

> **`ε ≤ δ` suffices unconditionally, uniformly in time and in the datum, on
> `T³`.**

The ε-monotone `Λ_ε(t) = log(N₀²(t)+ε) − (1/2ν)∫₀ᵗ G₀/(H₁+εH₀) ds` is
non-increasing whenever (M) is. **Note where the regularisation actually
matters:** on `T³` it buys nothing (no spectral gap issue); it is only needed on
`R³`, where `N₀²` may approach 0 and `log N₀²` is unbounded below — the same
place where the `H₀(T) > 0` obligation O-3 fails. This is a consistency check
on the choice of Clay target (B) rather than (A), not a new capability.

---

## 11. Recommendation

1. **O-8 is refuted.** Replace "promote 𝔡₀ > 0 to a uniform lemma" in §G by the
   theorem `inf d_* = 0` with the exact witness family of §5.2, and record K3 as
   fired. The *identity* (I.1)–(I.4) and Lemma K survive untouched; only the
   quantitative sharpening dies, exactly as kill condition K3 specified.
2. **O-9's energy-neutrality route is dead at the critical order `r = 0`**: (N)
   shows `T₀ = 0` is implied by saturation. Delete the sentence in §D
   ("the only escape is the double-saturation corridor") insofar as it suggests
   the (I.4) half is hard to saturate — it is easy; the corridor is closed, if at
   all, by the *Bernstein* half alone (`‖u‖_∞² ≍ N‖∇u‖₂²`, which the
   near-saturating two-mode-core family manifestly does **not** satisfy, since a
   bounded number of modes cannot saturate Bernstein). That is the only surviving
   O-9 ingredient and it should be restated as such.
3. **Pursue O-8′ (§8.4) only after fixing the search.** The spread-conditional
   deficit `d_* ≥ c·min(1, V₀/μ²)` is the strongest statement consistent with
   the counterexample family, it is scale-invariant, and it couples the front
   bound to broadbandness — which is what a cascade argument needs anyway. But
   this sprint's penalised descent is demonstrably unreliable (§8.4: it returns
   `1.23` where the free search already found `0.126`), so the reported `c` is
   an artefact either way. **Minimal next experiment: replace the penalty with a
   proper equality-constrained optimiser (projected gradient on the level set
   `H₀H₂/H₁² = 1+v`, whose gradient is already implemented as
   `spread_gradient`), rerun at `|k|_∞ ≤ 2,3,4`, and only then judge whether `c`
   is band-uniform.** If it is not, O-8′ dies too and the whole quantitative
   lane should be closed, leaving (I.1)–(I.4) + Lemma K + the dichotomy as the
   candidate's entire content.
4. **Correct the ledger convention** (§4) in `spectral_front_monotone`'s
   docstring and in the candidate's §2/§4, and re-report the pilot deficits at
   `ν = ν_*` rather than at `ν = 1/40` (§3).

### Binding non-claims

Nothing here proves or disproves any Clay statement. The refutation of O-8 is a
statement about a family of *fields*, not about any Navier–Stokes trajectory;
whether such near-degenerate configurations occur along a solution is not
addressed. Conjecture O-8′ is measured, not proved. No numerical output above
may be reported as evidence for or against a singularity.
