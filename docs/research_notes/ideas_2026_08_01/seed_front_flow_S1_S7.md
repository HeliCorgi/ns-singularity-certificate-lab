# SEED DERIVATIONS (orchestrator, 2026-08-01) — status: FORMAL ANSATZ unless noted

Conventions match the repo: on T³, u(x)=Σ_k û_k e^{ik·x}, NS in Fourier form
∂_t û_k = 𝒩_k − ν|k|²û_k,  𝒩_k = −i P_k Σ_{ℓ+m=k} (m·û_ℓ) û_m,  P_k = I − k⊗k/|k|²,
critical shell normalization E_N = c_E/N, heat factor h_{ν,τ}(r) = (1−e^{−ντr²})/(νr²).

## S1. Carrier-frame front flow (the continuum home of 𝔗₂) — FORMAL ANSATZ

Ansatz for the moving front: û(k,t) = N(t)^{−2} Ψ(k/N(t), s), rescaled clock ds = N(t)² dt.
Substituting (with ξ = k/N, and the lattice sum Σ_ℓ → N³∫dλ as a Riemann sum):

∂_t û = −(Ṅ/N³)(2Ψ + ξ·∇_ξΨ)·N⁰ + Ψ_s,
−ν|k|²û = −ν|ξ|²Ψ,
𝒩 → −𝒬(Ψ,Ψ)(ξ),  𝒬(Ψ,Ψ)(ξ) := i P_ξ ∫ ((ξ−λ)·Ψ(λ)) Ψ(ξ−λ) d³λ.

All terms are O(N⁰). The **front flow**:

  ∂_s Ψ = a(s)(2Ψ + ξ·∇_ξΨ) − ν|ξ|²Ψ − 𝒬(Ψ,Ψ),   a(s) := Ṅ/N³ = −(1/2) d(N^{−2})/dt.

- N^{−2}(t) = N₀^{−2} − 2∫₀ᵗ a; a(t) ≥ a₋ > 0 on the orbit ⟹ N→∞ at T ≤ N₀^{−2}/(2a₋).
  If a → a₊ constant: N ~ (2a₊(T−t))^{−1/2} (γ = 1/2, Type-I boundary; matches STATUS front ODE Ṅ = kN³).
- In s, dN/ds = aN (exponential); T − t(s) = ∫_s^∞ N^{−2} ds′ < ∞ automatically.
- The lattice ℤ³ breaks continuous dilation ⟹ the natural recurrent object is an s-periodic orbit
  with period S = log2 / a₊ (discrete doubling), i.e. exactly the projective periodic point of 𝔗₂.
  𝒬 restricted to carrier cells IS eq (6.6) of CANDIDATE_SOLUTION_PHASE_CODED_LERAY_CLOUD.md.
- KEY METHODOLOGICAL UNLOCK: an *attracting* periodic orbit of this flow is found by forward
  integration in s (renormalization-group integration) — no nonconvex fixed-point optimization.
  The frozen-parent one-step map + doubling pullback is one explicit-Euler step of this flow.

## S2. 1D shell reduction, coherence eigenvalue, front extent — FORMAL ANSATZ

Shell energy density e(ξ,s) = ξ² ∮_{S²} |Ψ(ξω,s)|² dω, front energy = N^{−1}∫e dξ.
From the flow: ∂_s e = a(2e + ξ e′) − 2νξ²e − ∂_ξ F, with radial transfer flux F (F(0)=F(∞)=0,
∫∂_ξF = 0 by exact advection energy-neutrality).
Critical wake e_c(ξ) = 2c_E ξ^{−2} (matches E_k = c_E/k at k = Nξ). **The dilation term
vanishes identically on e_c**: a(2e_c + ξe_c′) ≡ 0 — the critical spectrum is transparent
to rescaling. On the wake the steady equation gives exact linear flux depletion:
F′(ξ) = −2νξ²e_c = −4νc_E  ⟹  F(ξ) = F(1) − 4νc_E(ξ−1).
Coherent-capacity closure F(ξ) = χ ξ e^{3/2} (χ = net shape efficiency, losses subtracted;
consistent with repo capacity Π ≤ C N √M_N E_N^{3/2} after Π_phys = N·F):
F(1) = χ(2c_E)^{3/2} ⟹ **front similarity extent ξ_max = 1 + (χ/√2)·√c_E/ν**, and the
handoff condition ξ_max ≥ 2 gives the closure threshold **χ√c_E ≥ √2 ν**.
Predicted measurable profile: e(ξ) = [(F(1) − 4νc_E(ξ−1))/(χξ)]^{2/3} on 1 ≤ ξ ≤ ξ_max.

## S3. c_E-collapse: the problem reduces to shape closure only — SYMBOLIC CANDIDATE (algebraic)

Every quadratic channel (child flux, spill, difference, low return, cross-talk) scales as
c_E^{3/2}; viscous loss scales as ν c_E. All channel-to-channel RATIOS are c_E-independent
pure shape functionals. Hence:
(i) if the recurrent shape has net forward efficiency χ_shape > 0, then choosing
    c_E ≥ 2ν²/χ_shape² closes the budget (this is what the measured "c_E ≈ 228 required"
    row instantiates — large c_E is legal: total front energy Σ c_E/N_j = 2c_E/N₀ is made
    small by starting deep, N₀ large);
(ii) the ONLY scale-free obstruction is shape recurrence with χ_shape > 0.
So the decisive question for the whole lane: does the front flow (S1) admit a recurrent
(s-periodic, projective) profile with strictly positive net forward flux? Forward RG
integration answers this without optimization.

## S4. Flux→L³ per-octave lemma — PROOF CANDIDATE (finite-Fourier, Lean-able)

For band-limited shells (support in |k| ∈ [N/κ₀, κN]), Hölder + Riesz/Littlewood-Paley
boundedness + Bernstein ‖∇u_N‖₃ ≤ κN‖u_N‖₃ give
  |Π_N| = |⟨P_N B(u,u), P_N u⟩| ≤ C κ N ‖u_{band}‖₃³.
Sustained critical transfer Π_N ≥ q c_E N (needed for the relay budget) forces the
per-octave lower bound ‖u_{band N}‖₃³ ≥ (q/C) c_E^{3/2}... (constants to fix; at critical
normalization Π_N ~ N F(1) ~ χ(2c_E)^{3/2} N). Wake persistence: heat flow is a positive
real Fourier multiplier ⟹ phases are preserved; each completed octave retains a constant
energy fraction e^{−2ντ_stage-ish} until T. Hence #active-plus-wake octaves ~ log N(t) and
‖u(t)‖₃³ ≳ const · log N(t) → ∞: ESS-required divergence follows *automatically* from a
sustained flux budget. This removes proof-chain steps 8–9 of the cloud candidate as
separate obligations (they follow from step 5). New, cheap, formalizable.

## S5. Helicity ledger — SPECULATION→SYMBOLIC (cheap exact checks)

Every linearly-polarized real carrier mode (all current relay alphabets) has exactly zero
helicity. Max-helicity (Beltrami) fields have B(u,u) = ∇(|u|²/2) ⟹ zero Leray transfer:
coherent transfer requires helicity-neutrality. Per-shell helicity capacity |H_N| ≤ 2N E_N
= 2c_E: a helicity-coherent critical cascade would accumulate |H| ~ #octaves·2c_E while
|H| ≤ ‖u‖₂‖ω‖₂ ~ enstrophy^{1/2} — needs care but suggests a *design rule* (alternating or
zero helicity code) plus possibly a new no-go for single-sign helical critical cascades.
Exact rational helicity columns can be added to the carrier search at trivial cost.

## S6. Cross-talk removal by scale-stagger (multi-type cycle) — SPECULATION with exact test

Known failure: two same-scale relays have diagonal parent pairs (A₁+B₂) landing exactly on
the target child shell (leakage 2483/222 vs intended). Two candidate fixes, one bad, one testable:
(a) translation-split (different centers x₀, x₁ for the two relays): kills diagonal
    cross-talk by in-box decoherence (suppression ~ (W|Δx|)^{−3} in power) BUT also kills
    the *next-stage intended pair* (child₁ from x₀, child₂ from x₁) by the same factor —
    REJECTED unless the stage-2 pairing is redesigned.
(b) scale-stagger: run relay-2 at scale ρN, ρ ≠ 1 rational (e.g. 3/4): diagonal cross
    outputs land at (1+ρ)-type mixed radii that can be made to MISS every tagged carrier
    box (exact support/geometry exclusion, checkable by integer programming over the box
    lattice), while the intended stage-2 pair (child₁ at 2N-band, child₂ at 2ρN-band) is
    designed to land back on a scaled copy of the carrier system: an L-periodic multi-type
    cycle 𝔗_{L−1}∘…∘𝔗₀ with Πλ_j = 2 per super-period. The candidate doc's eq (6.11)
    anticipates multi-type cycles; the new content is using scale-stagger *specifically to
    turn diagonal cross-talk into a geometric (exactly zero) channel*.
Immediate exact test: enumerate small rational ρ and box widths, check sumset-miss
conditions with the existing box-tag machinery.

## S7. Pilot design: renormalized Galerkin cascade (RG iteration)

Iterate: (parent cloud at N, critical energy) → full Galerkin evolution (cutoff ≥ sumset,
dealiased) for τN^{−2} → doubling pullback (child band 2k → k, coefficient renormalization,
energy renormalized to critical; record pre-normalization gain 𝒜_j and shape overlap
ρ_j = |⟨U_{j+1},U_j⟩|/‖·‖‖·‖) → repeat. Success: 𝒜_j → 𝒜₊ > 0 (ANY positive constant, by
S3 c_E-collapse) and ρ_j → 1 for some (η, c_E, τ, ν) in a small scan. Kill: ρ_j decreasing
under refinement or 𝒜_j → 0 across the entire scan. Reuses mesoscopic_galerkin /
mesoscopic_local_fft / carrier_two_stage_galerkin infrastructure.
