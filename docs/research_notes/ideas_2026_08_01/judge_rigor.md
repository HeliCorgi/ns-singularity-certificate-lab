# JUDGE REPORT — 14 candidates, NS certificate lab (2026‑08‑01)

Axes: **Clay** = distance‑to‑Clay (5 = closest) · **Nov** = novelty vs literature *and* repo · **Clos** = closability · **NoGo** = evasion soundness · **Cmp** = computability · **IA** = interval feasibility · **Lean** · **Fals** = falsifiability · **Con** = genuine connection to singularity/regularity.

| # | candidate | Clay | Nov | Clos | NoGo | Cmp | IA | Lean | Fals | Con | Σ |
|---|---|--|--|--|--|--|--|--|--|--|--|
|11|**monotone-dichotomy**|3|4|4|5|5|4|5|5|4|**39**|
|8|**log-periodic-multitype**|3|4|3|4|5|4|4|5|3|**35**|
|13|**moving-front-ode**|3|4|3|4|4|3|4|4|4|**33**|
|7|helicity-ledger|2|2|3|3|5|4|4|5|3|31|
|1|**aniso-multiscale**|4|4|2|4|3|3|3|3|4|30|
|5|**lagrangian-degeneracy**|3|4|3|3|3|2|4|4|4|30|
|6|parametric-resonance|2|5|1|2|5|4|4|5|2|30|
|9|forcing-controller|2|3|4|3|4|4|4|4|2|30|
|10|critical-quantity|2|3|3|4|4|3|4|4|3|30|
|12|two-scale-core|3|4|2|3|4|3|3|4|4|30|
|4|pressure-feedback|3|4|1|3|3|3|2|4|2|25|
|2|multi-center|3|3|1|2|3|3|2|3|3|23|
|3|space-fourier-mismatch|4|3|1|3|3|2|1|3|3|23|
|14|front-flow (seed)|3|2|1|1|4|3|2|4|3|23|

Per‑candidate one‑liners (load‑bearing axis in bold): **11** all identities re‑derived and confirmed here — closable, pressure‑free, phase‑free. **8** LP‑5 and LP‑9 verified by hand; kills the *registered* cross‑talk gate geometrically. **13** front–BKM identity exact; capacity floor independent of BKM. **7** cheap and exact but T1 **is Waleffe (1992)** — low novelty. **1** only serious TARGET‑U attack; anisotropy ceiling exact, but Prandtl‑ill‑posed core. **5** L1/L3 correct, cheap, Lean‑able; mechanism contradicted (below). **6** Prop. 1 genuinely new; exhibited gadget does not scale. **9** admissibility proof solid; targets the weaker (D). **10** dichotomy rigorous but quantitatively vacuous. **12** LD no‑go and decoupling estimates real; saturates KNSS(b). **4** depolarization⇒restricted‑Euler is sharp; closure self‑inconsistent. **2** exact skew‑strain formula verified; generation map undefined. **3** exponent corner arithmetic checks out; mechanism narrative discarded mid‑derivation. **14** already audited; closure over‑determined.

## FATAL flaws (exact equation / no‑go collision)

- **front-flow.** `e≡e_c` and `F=χξe^{3/2}` are inconsistent: capacity `F(1)ξ⁻²` lies below requirement `F(1)(1−β(ξ−1))` for all `1<ξ<φ` at every `β≤1<2`. Plus **T.2**: the steady limit *is* a Leray profile with finite local energy ⇒ **Tsai1998 ⇒ u≡0**; T.4 makes "attracting near‑steady profile" a *failure* signature.
- **lagrangian-degeneracy.** §B.4 derives the γ‑window from `|ω|≍σ₁` — the **Cauchy formula** — while its own **L3** proves that formula inadmissible at any collapsing front (`ν∫τ^{−3γ−1}dτ=∞` for every γ>0). The γ∈[1/2,1) selection is therefore unsupported. L1/L3 survive independently and are correct.
- **parametric-resonance.** Net growth needs `χ_par sX_P>ν|k|²`; at critical normalization `X_P≍√(2c_E/3N)`, `s≍N`, giving **`c_E ≳ ν²N³/χ_par²→∞`**. The exhibited 6‑mode gadget is fixed‑cardinality — exactly the family `fixed_cardinality_scaling` already rejects at `N^{−3/2}`. All content rests on untested K2.
- **pressure-feedback.** Loop step 3 gives `ĝ≈ℓ₁ℓ₂/(4πρ_⊥²)`; the family's `ĝ≈0.985` forces `ρ_⊥≈0.28√(ℓ₁ℓ₂)` — the "external" companion sits **inside** the core, where eq. (5) is invalid. V3 additionally shows `p₂^mat>1/2` along the *entire* family, so §C's exponents are fitted, not derived.
- **multi-center.** `M_eff ≥ λ³/n` is a **lower** bound — it cannot *force* `n=O(1)` (non sequitur without the decoherence equality). Separately, a sub‑dyad of two Γ‑tubes inside each tube doubles that tube's circulation, contradicting `n_jΓ_j²=c_E` and Kelvin; the generation map is undefined.
- **space-fourier-mismatch.** (B.5)'s self‑reinforcing mismatch clock is *replaced* by (B.6) pinning, in which N is slaved to L — the advertised mechanism is not the mechanism used. And `c_E=E₀/2L→∞` ⇒ `ν_eff=ν/√c_E→0`: the corner is the **inviscid** problem relabelled.
- **critical-quantity.** Theorem D's floor `ν/(3(‖u‖₂‖∇u‖₂)^{1/2})→0` at blow‑up: branch (B) can never be contradicted. Rigorous, asymptotically vacuous. (B10)/(B11) need the unproven n₀‑locality + signed‑flux hypothesis.
- **two-scale-core.** Type‑I with rescaled `|V|≲α/|y|` **saturates KNSS(b)** exactly; survival rests on the single hypothesis "non‑axisymmetric". Its own (M2) `F(ξ₀)=0` is conceded to be an output, so forward flux is assumed.
- **aniso-multiscale.** `M_eff/N³ ∝ τ^{3/16}→0` — L‑11 *would* kill it; escape is a contested (though defensible) reading of "empty child / frozen parent". Prandtl ill‑posedness admitted.
- **forcing-controller / helicity-ledger.** No fatal error; but (9.x) leaves every hard obligation intact, and T1/T2 duplicate Waleffe and the seed audit's own S5.5 repair.
- **monotone-dichotomy, log-periodic-multitype, moving-front-ode.** No fatal flaw found. I re‑verified (I.2)→(I.3) square completion, Lemma K's `4πN` lattice sum, the γ≥1/2 window, `𝔅²=2·Id`, LP‑5 (`P_c(X+Y)=⅓(1,−1,−1)`, `P_{p−q}(Y−X)=(−1,1,1)`), LP‑9's divisibility argument, and the (CAP)→(ODE)→`g≥2/(4+a)` chain. All correct.

## Top 5

1. **monotone-dichotomy (Lens 11)** — highest rigor density; exactly rational, Lean‑able today; re‑derives γ∈[1/2,1) with no BKM, phase, or flux‑sign assumption; pressure‑blind, so immune to the LG‑9/VR‑C‑009 pathologies that killed Gate‑7/8.
2. **log-periodic-multitype (Lens 8)** — two hand‑verified exact theorems that attack the repo's *registered blocker* (16/16 cross‑talk failures) with a support‑exclusion, not an amplitude estimate; testable log‑4 signature.
3. **moving-front-ode (Lens 13)** — exact front–BKM identity plus the disciplinary consequence that BKM divergence is *never* numerically observable; kills random‑phase fronts unconditionally.
4. **aniso-multiscale (Lens 1)** — the only candidate genuinely aimed at CLAY‑A/TARGET‑U; anisotropy‑ceiling `a₁−a₃<1/6` is exact, machine‑checkable, and a real new budget.
5. **lagrangian-degeneracy (Lens 5)** — L1 (σ₃→0 ⇒ BKM divergence) and L3 (Euler‑perturbative lane permanently closed) are correct, cheap, Lean‑able, and survive their own mechanism's collapse.

## Single best 1‑hour pilot

**Lens 8's P‑LP1 + P‑LP2 (Butterfly both‑channel exact carrier search + odd‑denominator sumset‑miss enumerator).**

Why: it is the only ≤1 h experiment that hits **both** registered blockers at once — signed flux positivity (gates G2/G3) and diagonal cross‑talk (LP‑9) — entirely in exact `Fraction`/integer arithmetic, with zero new physics code. Reuse `expanded_carrier_search.canonical_waves_in_box(2)`, `.primitive_polarizations(wave,3)`, `.projected_mixed_channels`; independently re‑verify hits through `exact_carrier_record_verifier.verify_serialized_expanded_carrier_certificate` (mandatory, TM‑14). P‑LP2 is a ~30‑line integer enumeration over `P,Q≤15`.

Pre‑registered outcome: **pass** = ≥1 orientation with `P_{a+b}(X+Y)≠0`, `P_{a−b}(Y−X)≠0`, both signed fluxes >0, *and* both second‑stage channels `Π_{2a},Π_{2b}>0`, with `η_*(5/3)≥1/14` confirmed by brute force. **Kill** = zero hits at component bound 3 (Butterfly has no exact carrier realization in scope) or `η_*=0` for all admissible ρ (LP‑9 bookkeeping wrong). Free 10‑minute add‑on sharing `modal_front_actions._moments`: Lens 11's `front_gap_identity`, asserting `Γ^CS_r,Γ^SC_r≥0` and (I.4) as an exact rational equality — any negative value falsifies Lens 11 outright.