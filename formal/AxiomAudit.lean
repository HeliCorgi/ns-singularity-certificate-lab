import NSSingularity.ClayStatement
import NSSingularity.VelocityRecovery
import NSSingularity.FiniteTime
import NSSingularity.GalerkinNoBlowup
import NSSingularity.FiniteModeNoGo
import NSSingularity.GreenAndCascade
import NSSingularity.CertificateLayer

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

-- Limited Clay connection and the F-7c reduction.
#print axioms NSSingularity.breakdown_time_set_empty
#print axioms NSSingularity.galerkin_solution_of_autonomised

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
