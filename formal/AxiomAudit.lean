/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
import NSSingularity.ClayStatement
import NSSingularity.VelocityRecovery
import NSSingularity.FiniteTime
import NSSingularity.GalerkinNoBlowup
import NSSingularity.FiniteModeNoGo
import NSSingularity.GreenAndCascade
import NSSingularity.CertificateLayer
import NSSingularity.TimeDependentGalerkin
import NSSingularity.ControlODE
import NSSingularity.L3Generation
import NSSingularity.TrackPFourier
import NSSingularity.GaussianTransfer
import NSSingularity.TrackPChain
import NSSingularity.GalerkinPicard
import NSSingularity.KatoConstant
import NSSingularity.ChainAnalysis
import NSSingularity.MesoscopicDuhamelNoGo

/-!
# Axiom audit (P0 Lean gate)

Not imported by the library root: this file exists so that
`lake env lean AxiomAudit.lean` prints the axiom set of every core theorem.
The acceptance rule (LEAN4_VERIFICATION_POLICY.md) is that each line reports
at most `[propext, Classical.choice, Quot.sound]` — the three standard
mathlib classical axioms — and never a project-specific axiom, `sorry`
(`sorryAx`), or `admit`.

`8659 jobs` in `lake build` counts compilation jobs, not proved theorems;
the proved statements are exactly the ones listed here and in the two files.
-/

open NSSingularity

-- F-3: velocity recovery implies the divergence-free identity.
#print axioms NSSingularity.mixed_partial_comm
#print axioms NSSingularity.divergence_of_recovered_velocity_eq_zero'
#print axioms NSSingularity.divergence_of_recovered_velocity_eq_zero

-- F-2: integrable scale rate gives a finite physical blow-up time.
#print axioms NSSingularity.hasDerivAt_physicalTime
#print axioms NSSingularity.physicalTime_strictMonoOn
#print axioms NSSingularity.tendsto_physicalTime
#print axioms NSSingularity.physicalTime_lt_blowupTime
#print axioms NSSingularity.exists_finite_blowupTime
#print axioms NSSingularity.tendsto_physicalTime_atTop

-- F-6: an energy-neutral dissipative quadratic system has no finite-time
-- blow-up (the Track-F finite-mode obstruction).
#print axioms NSSingularity.norm_le_of_energy_inequality
#print axioms NSSingularity.inner_galerkin_le
#print axioms NSSingularity.galerkin_norm_le
#print axioms NSSingularity.galerkin_norm_le_of_mem
#print axioms NSSingularity.galerkin_not_tendsto_atTop

-- F-12: the trilinear cancellation in the Fourier representation.
#print axioms NSSingularity.advectionForm_eq_zero

-- F-13: norm equivalence constants on the finite Fourier space.
#print axioms NSSingularity.weighted_sq_sum_le
#print axioms NSSingularity.sq_sum_abs_le_card_mul_sum_sq
#print axioms NSSingularity.sum_abs_le_sqrt_card_mul_sqrt_sum_sq

-- F-7a: a bounded trajectory reaches the finite-time endpoint.
#print axioms NSSingularity.intervalIntegrable_of_continuousOn_bounded
#print axioms NSSingularity.exists_tendsto_nhdsWithin_of_norm_deriv_le

-- F-7b: local continuation of the autonomous Galerkin system.
#print axioms NSSingularity.contDiff_galerkinField
#print axioms NSSingularity.exists_local_galerkin_solution

-- The logical connection: a fixed-finite-mode candidate never breaks down.
#print axioms NSSingularity.not_isBreakdownCandidate_of_galerkin
#print axioms NSSingularity.galerkin_bounded_and_reaches_endpoint

-- F-14: the five-dimensional radial Green profile is harmonic off the origin.
#print axioms NSSingularity.hasDerivAt_greenProfile
#print axioms NSSingularity.hasDerivAt_greenProfileDeriv
#print axioms NSSingularity.greenProfile_radial_laplace_eq_zero

-- F-15: Newton's flux identity.
#print axioms NSSingularity.flux_newtonSlope
#print axioms NSSingularity.hasDerivAt_flux

-- F-16: the shell exponent region, with every hypothesis named.
#print axioms NSSingularity.ShellAdmissible.bandwidth_lt_one
#print axioms NSSingularity.ShellAdmissible.sigma_mem
#print axioms NSSingularity.not_shellAdmissible_of_one_le

-- Limited Clay connection, and the abandoned autonomised route to F-7c.
#print axioms NSSingularity.breakdown_time_set_empty
#print axioms NSSingularity.galerkin_solution_of_autonomised

-- F-7c, closed: local existence for a genuinely time-dependent projected force,
-- taken by the direct route through mathlib's time-dependent `IsPicardLindelof`.
#print axioms NSSingularity.galerkin_isPicardLindelof
#print axioms NSSingularity.galerkin_local_solution
#print axioms NSSingularity.galerkin_local_solution_of_continuous

-- F-17: potential error to velocity error.
#print axioms NSSingularity.velocity_radial_error_le
#print axioms NSSingularity.velocity_axial_error_le

-- F-18: the product-difference identity and the advection error.
#print axioms NSSingularity.product_difference
#print axioms NSSingularity.product_error_le
#print axioms NSSingularity.advection_error_le

-- F-19: the short-time Gronwall step.
#print axioms NSSingularity.gronwallBound_le_simple
#print axioms NSSingularity.norm_le_simple_gronwall

-- The packaged Clay restriction.
#print axioms NSSingularity.FixedBandwidthCandidate.breakdown_times_empty
#print axioms NSSingularity.FixedBandwidthCandidate.reaches_every_time

-- The control-ODE layer: the Chaplygin-Dini comparison lemma.
#print axioms NSSingularity.nonpos_of_deriv_le_mul_of_pos

-- HS-6: Gronwall with a time-dependent coefficient, in both forms.
#print axioms NSSingularity.gronwall_variable_coefficient
#print axioms NSSingularity.gronwall_variable_coefficient_integral

-- The undamped quadratic control ODE and its explicit blow-up majorant.
#print axioms NSSingularity.riccati_comparison
#print axioms NSSingularity.le_quadratic_bound

-- The rough enclosure: the rigorous replacement for the slab certificate's H2.
#print axioms NSSingularity.clampTo_mem
#print axioms NSSingularity.clampTo_eq_self
#print axioms NSSingularity.abs_clampTo_sub_clampTo_le
#print axioms NSSingularity.roughEnclosure_isPicardLindelof
#print axioms NSSingularity.exists_roughEnclosure_solution
#print axioms NSSingularity.roughEnclosure_solution_unique

-- The L^3 generation identity, algebraic core: the chain rule at exponent
-- three and the transport cancellation.
#print axioms NSSingularity.hasDerivAt_norm_cube
#print axioms NSSingularity.hasDerivAt_cube
#print axioms NSSingularity.transport_eq_one_third_deriv

-- The regularised speed used in place of |u| by the interval certificate.
#print axioms NSSingularity.eps_le_regSpeed
#print axioms NSSingularity.norm_le_regSpeed
#print axioms NSSingularity.regSpeed_le_norm_add
#print axioms NSSingularity.contDiff_regSpeed
#print axioms NSSingularity.hasDerivAt_regSpeed

-- Signs, the Kato split without division by the speed, and the scaling
-- criterion that defines Re_crit.
#print axioms NSSingularity.viscous_contribution_nonpos
#print axioms NSSingularity.kato_split_le
#print axioms NSSingularity.positive_generation_forces_pressure
#print axioms NSSingularity.generation_pos_iff_reynolds_gt

-- The pure-swirl no-go: divergence freedom and orthogonality to the gradient
-- of every axisymmetric scalar.
#print axioms NSSingularity.swirl_cartesianDiv_eq_zero
#print axioms NSSingularity.swirl_dot_grad_axisymmetric_eq_zero

-- The equality case of the no-go.
#print axioms NSSingularity.eventually_eq_of_fderiv_eq_zero_on
#print axioms NSSingularity.eq_of_locallyConstant_of_ne_zero
#print axioms NSSingularity.eq_zero_of_locallyConstant_of_tendsto_cocompact
#print axioms NSSingularity.pure_swirl_equality_case

-- Track P (periodic lane), F1: the Leray multiplier as finite algebra.
#print axioms NSSingularity.inner_leray_eq_zero
#print axioms NSSingularity.leray_eq_self_of_inner_eq_zero
#print axioms NSSingularity.leray_leray
#print axioms NSSingularity.inner_leray_left_eq_inner_leray_right
#print axioms NSSingularity.norm_leray_le

-- Track P, F2: a single Fourier mode with transverse amplitude is
-- divergence free as a classical field.
#print axioms NSSingularity.slotDivergence_cosMode
#print axioms NSSingularity.slotDivergence_sinMode

-- Track P, F3: finite trigonometric polynomials are C-infinity.
#print axioms NSSingularity.contDiff_trigPolynomial

-- Track P, F4: the fixed-band / finite-band distinction, its one-line
-- implication, its counterexample, and the restated scope of the no-go.
#print axioms NSSingularity.FixedBandTrajectory.finiteBandDatum
#print axioms NSSingularity.exists_finiteBandDatum_not_fixedBandTrajectory
#print axioms NSSingularity.FixedBandwidthCandidate.fixedBand_scope

-- Track P, F5: the homogeneous-norm ladder for weights at least one.
#print axioms NSSingularity.weighted_sum_succ_mono
#print axioms NSSingularity.weighted_sum_pow_mono

-- Track P, F6: the slab composition shell through riccati_comparison.
#print axioms NSSingularity.trackP_slab_error_le

-- G1: the Gaussian-Hermite class is closed under differentiation, with the
-- witness polynomial exhibited.
#print axioms NSSingularity.hasDerivAt_poly_mul_gaussian
#print axioms NSSingularity.gaussianDerivPoly_eval
#print axioms NSSingularity.hasDerivAt_poly_gaussian

-- G2: the pointwise J-continuity bricks of the spline-to-smooth transfer.
#print axioms NSSingularity.abs_cube_sub_cube_le
#print axioms NSSingularity.abs_norm_cube_sub_norm_cube_le
#print axioms NSSingularity.norm_smul_sub_smul_le

-- G3: the Track-P torus specialisation of the Riccati comparison.
#print axioms NSSingularity.torus_control_bound

-- Track P chain, C1-C2: the two-slab composition and the transfer triangle.
#print axioms NSSingularity.two_slab_composition
#print axioms NSSingularity.transfer_triangle

-- Track P chain, C3: the n-slab induction and its union form.
#print axioms NSSingularity.chain_composition
#print axioms NSSingularity.chain_composition_union

-- Track P chain, C4: the discrete Gronwall inequality for the datum radii.
#print axioms NSSingularity.discrete_gronwall

-- Track P chain, C5: the chained bound never exceeds the per-slab maximum.
#print axioms NSSingularity.piecewise_radius_le_max
#print axioms NSSingularity.le_foldr_max
#print axioms NSSingularity.chain_radius_le_foldr_max

-- Track P chain, C6: the Lagrange endpoint bound behind the recentring
-- transfer.
#print axioms NSSingularity.taylor_endpoint_remainder_bound

-- Galerkin Picard-Lindelof, GP1: the quadratic field is Lipschitz on a ball
-- with the explicit constant.
#print axioms NSSingularity.quadraticField_apply
#print axioms NSSingularity.quadratic_field_lipschitzOnWith

-- GP2: local existence for the quadratic Galerkin ODE, closed-interval and
-- interior-derivative forms.
#print axioms NSSingularity.quadratic_ode_local_solution
#print axioms NSSingularity.quadratic_ode_local_solution_hasDerivAt

-- GP3: uniqueness of the quadratic Galerkin ODE solution in a ball.
#print axioms NSSingularity.quadratic_ode_unique

-- Kato constant, K1-K5: the finite algebra of kato_h3_constants.md section 9.
#print axioms NSSingularity.cube_diff_bound
#print axioms NSSingularity.am_gm_split
#print axioms NSSingularity.shifted_ratio_bound
#print axioms NSSingularity.inv_pow_four_succ_le_telescope
#print axioms NSSingularity.inv_pow_tail_bound
#print axioms NSSingularity.g3_assembly_mono
#print axioms NSSingularity.g3_of_a4

-- Chain analysis, A1: the integral-inequality comparison (EXT-P2, scalar
-- half) and its Riccati instance.
#print axioms NSSingularity.integral_comparison
#print axioms NSSingularity.integral_riccati_comparison

-- Chain analysis, A2: gluing, the uniform-modulus Cauchy bridge, and the
-- endpoint extension (EXT-P3, finite halves).
#print axioms NSSingularity.glued_continuous
#print axioms NSSingularity.cauchy_map_of_uniform_modulus
#print axioms NSSingularity.exists_continuousOn_Icc_extension
#print axioms NSSingularity.exists_continuousOn_Icc_extension_of_modulus

-- Chain analysis, A3: the certificate-dependency discharge shape.
#print axioms NSSingularity.cond_to_uncond

-- Mesoscopic Duhamel no-go, MD1-MD3: only the finite nonnegative-real
-- algebra and finite coefficient Cauchy--Schwarz.  There is no PDE/Fourier
-- bridge in this module.
#print axioms NSSingularity.emptyChild_duhamel_ratio_sq_le
#print axioms NSSingularity.emptyChild_duhamel_ratio_sq_le_of_effectiveCount_le
#print axioms NSSingularity.finiteEffectiveModeCount_nonneg
#print axioms NSSingularity.finiteEffectiveModeCount_le_card
#print axioms NSSingularity.emptyChild_duhamel_ratio_sq_le_card
