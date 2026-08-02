# Coherent critical-spectrum family — exact certificates (SHARP truncation)

> **Re-label (2026-08-02).** These are certificates for the **sharply**
> truncated family. They were the appendix to the then-open Lemma L\*.
> The paper no longer uses that family or that hypothesis: the static
> no-go is now **unconditional**, via a *proven* capacity bound
> \(\|\mathbb P(u_N\cdot\nabla u_N)\|_2^2\ge c_0N^3\) for the **smoothly**
> truncated family \(\chi(|k|/N)P_kv_0/|k|^2\)
> (`docs/paper_lambda_dichotomy/lstar/lstar_proof_main.md`, Thm 7.1(3)).
> The sharp family's capacity bound — the literal L\* — remains **OPEN**
> and is **used by nothing**. Everything below is retained, unchanged and
> still correct, as corroboration for that open statement. The smooth
> family's own exact certificates are in
> `docs/paper_lambda_dichotomy/lstar/lstar_numerical_support.md` and
> `outputs/lstar/smooth_family_capacity.json`.

Family (unchanged; `experiments/run_osgood_gate.py::coherent_field`):
\(\widehat u(k)=P_kv_0/|k|^2\) on \(1\le|k|\le N\), \(v_0=(1,2,3)\),
\(\widehat u(-k)=\widehat u(k)\).
Convention (unchanged; `exact_anchor`): pair energy \(=L^2\) mean, so
\(H_0=\|u\|_2^2\), \(H_1=\|\nabla u\|_2^2=D\),
\(K=\|\mathbb P(u\cdot\nabla u)\|_2^2/H_1^2\),
\(S_N=\sum_{1\le|k|\le N}|k|^{-2}\), \(T_N=\sum|k|^{-4}\).
Source: `experiments/run_coherent_family_certificates.py` →
`outputs/verification_sprint_v1/osgood_gate/exact_family_certificates.json`
(reproducible content sha256
`de3b7b595b4979ada59fba258f203d38316d097b95af116c0bb7cc8f2dc5e079`; the
`.sha256` sidecar covers the file including its wall-clock fields).
No new quantity is defined here.

**Sweeping-projection identity (one line).** For constant \(c\) and
divergence-free \(u\), the \(k\)-th Fourier coefficient of \((c\cdot\nabla)u\)
is \(i(c\cdot k)\widehat u(k)\perp k\) because \(k\cdot\widehat u(k)=0\); hence
\(\mathbb P((c\cdot\nabla)u)=(c\cdot\nabla)u\), so
\(Q(c)=\langle(u\cdot\nabla)u,(c\cdot\nabla)u\rangle
=\langle\mathbb P(u\cdot\nabla u),(c\cdot\nabla)u\rangle\) and Cauchy–Schwarz
gives \(\|\mathbb P(u\cdot\nabla u)\|_2\ge|Q(c)|/\|(c\cdot\nabla)u\|_2\).

---

**Table 1 — symmetry identities \(H_0=\frac23|v_0|^2T_N\),
\(H_1=\frac23|v_0|^2S_N\) (\(|v_0|^2=14\), \(\frac23|v_0|^2=\frac{28}{3}\));
all cells exact `Fraction`, full-band enumeration \(1\le|k|^2\le N^2\).**

| \(N\) | band pts | \(S_N\) | \(T_N\) | \(H_1\) | \(H_0\) | \(H_1-\frac{28}{3}S_N\) | \(H_0-\frac{28}{3}T_N\) | label |
|---|---|---|---|---|---|---|---|---|
| 4 | 256 | `4888669/120120` | `231195938461/17314577280` | `4888669/12870` | `231195938461/1855133280` | `0` | `0` | exact |
| 6 | 924 | `67892725261/1012647636` | `29642745162859844989261/2050910469392776992000` | `67892725261/108497961` | `29642745162859844989261/219740407434940392000` | `0` | `0` | exact |
| 8 | 2108 | `60280440054870170591832503/661963540931031504765600` | (JSON) | `60280440054870170591832503/70924665099753375510600` | (JSON) | `0` | `0` | exact |

**Table 1b — same rows, decimal.**

| \(N\) | \(S_N\) | \(T_N\) | \(H_0\) | \(H_1\) | \(N_0^2=H_1/H_0\) | label |
|---|---|---|---|---|---|---|
| 4 | 40.69821012 | 13.35267588 | 124.6249749 | 379.8499611 | 3.047944134 | exact→float |
| 6 | 67.04476745 | 14.45345646 | 134.8989269 | 625.7511628 | 4.638666720 | exact→float |
| 8 | 91.06308177 | 14.95291847 | 139.5605724 | 849.9220965 | 6.089987178 | exact→float |

---

**Table 2 — exact \(K\), \(Q(v_0)\), sweeping bound. Full \(|{\rm band}|^2\)
ordered-pair convolution, no truncation; every entry an exact `Fraction`
(decimals are its `float()`). \(N=6\) cost 0.51 s and \(N=8\) 3.03 s, so the
pre-registered 20-min fallback to \(N=5\) was NOT used and \(N=8\) was
certified as a bonus.**

| \(N\) | ordered pairs | sumset modes | \(K\) | \(N_0^2\) | \(K/N_0^2\) | \(\|\mathbb P(u\!\cdot\!\nabla u)\|_2\) | \(Q(v_0)\) | \(\|(v_0\!\cdot\!\nabla)u\|_2^2\) | \(|Q|/\|(v_0\!\cdot\!\nabla)u\|_2\) | captured | label |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 65 536 | 1 736 | 0.7884107043392768 | 3.047944134 | 0.2586696704 | 337.2782551 | 10252.58034507 | 1089.269007 | 310.6458870 | **0.9210373997** | exact |
| 6 | 853 776 | 6 780 | 1.4344718079662586 | 4.638666720 | 0.3092422661 | 749.4586458 | 28954.79916692 | 1775.305214 | 687.2014112 | **0.9169303938** | exact |
| 8 | 4 443 664 | 15 904 | 2.0372430105569106 | 6.089987178 | 0.3345233662 | 1213.110994 | 54445.60989576 | 2405.446015 | 1110.107559 | **0.9150914999** | exact |

**Table 2b — exact rational values (`Fraction`, reduced).**

| quantity | \(N=4\) |
|---|---|
| \(K\) | `56891190666086268129016529590244392986924943/72159333140667601589384040100247446050168000` |
| \(\|\mathbb P(u\cdot\nabla u)\|_2^2\) | `56891190666086268129016529590244392986924943/500113224073356389376340375883167200000` |
| \(Q(v_0)\) | `55131961007981920458419/5377374197755560000` |
| \(\|(v_0\cdot\nabla)u\|_2^2\) | `96751909691857/88822787625` |
| \(N_0^2\) | `704672304336/231195938461` |

\(N=6,8\) exact rationals (\(K\), \(\|\mathbb P(u\cdot\nabla u)\|_2^2\),
\(Q(v_0)\), \(\|(v_0\cdot\nabla)u\|_2^2\), \(|Q|^2/\|(v_0\cdot\nabla)u\|_2^2\),
captured\(^2\)): `exact_certificates[*]` in the JSON above.

---

**Table 3 — exact vs. float pipeline (`coherent_field` + `measure`, rebuilt
unmodified) at the certified bands. \(H_0,H_1\) compared after the convention
map \(H_0=2\cdot\)`measure["H0"]`, \(H_1=\)`measure["grad_sq"]`\(=2\cdot\)`measure["H1"]`
(\(K\) is convention-free). The exact lane evaluates \(Q(v_0)\) WITHOUT the
Leray projection and the float lane WITH it; their agreement is an independent
numerical confirmation of \(\mathbb P((v_0\cdot\nabla)u)=(v_0\cdot\nabla)u\).**

| \(N\) | grid | rel. diff \(K\) | rel. diff \(Q(v_0)\) | rel. diff \(H_0\) | rel. diff \(H_1\) | label |
|---|---|---|---|---|---|---|
| 4 | 48 | 0.0 | 1.77e-16 | 1.14e-16 | 0.0 | exact vs float |
| 6 | 48 | 4.64e-16 | 1.26e-16 | 0.0 | 1.82e-16 | exact vs float |
| 8 | 80 | 1.31e-15 | 5.35e-16 | 0.0 | 0.0 | exact vs float |

---

**Table 4 — float continuation (binary64, dealias-safe grids; \(N=4,6,8\) rows
reproduce Table 2 to the last column of Table 3, \(N=16,32\) reproduce the
published \(K=4.5741,\ 9.6478\) of `spectral_front_osgood_gate.md` §2 on a
grid-144 rerun).**

| \(N\) | grid | \(K\) | \(K/N_0^2\) | \(Q(v_0)\) | \(Q(v_0)/N^2\) | \(\|\mathbb P(u\!\cdot\!\nabla u)\|_2\) | \(\|\mathbb P(u\!\cdot\!\nabla u)\|_2/N^{3/2}\) | \(|Q|/\|(v_0\!\cdot\!\nabla)u\|_2\) | captured | label |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 48 | 0.7884107043 | 0.2586696704 | 10252.58035 | 640.786272 | 337.27826 | 42.1598 | 310.6459 | 0.9210373997 | float |
| 6 | 48 | 1.4344718080 | 0.3092422661 | 28954.79917 | 804.299977 | 749.45865 | 50.9942 | 687.2014 | 0.9169303938 | float |
| 8 | 80 | 2.0372430106 | 0.3345233662 | 54445.60990 | 850.712655 | 1213.11099 | 53.6124 | 1110.1076 | 0.9150914999 | float |
| 16 | 144 | 4.5740893285 | 0.3754340731 | 248104.76791 | 969.159250 | 3829.32674 | 59.8332 | 3494.8419 | 0.9126517995 | float |
| 32 | 144 | 9.6477871736 | 0.3961869815 | 1054176.19524 | 1029.468941 | 11393.74207 | 62.9421 | 10388.3654 | 0.9117606233 | float |

---

**Table 5 — pre-registration vs. outcome (caps declared before the run).**

| item | pre-registered | outcome | binding? |
|---|---|---|---|
| exact bands | \(N\in\{4,6,8\}\) | 4, 6, 8 certified | no |
| exact enumeration | all \(1\le|k|^2\le N^2\); all \(|{\rm band}|^2\) ordered pairs | 65 536 / 853 776 / 4 443 664 pairs, no truncation | no |
| exact wall cap | 20 min | 3.61 s | no |
| \(N=6\) fallback to \(N=5\) | if \(N=6>20\) min | not used | no |
| float wall cap | 10 min | 16.86 s | no |
| float bands | \(N\in\{8,16,32\}\) (+4, 6 for cross-check) | as declared | no |
| total wall cap | 30 min | \(\approx21\) s | no |
