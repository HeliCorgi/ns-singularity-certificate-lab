import NSSingularity.ClayStatement
import NSSingularity.VelocityRecovery
import NSSingularity.FiniteTime
import NSSingularity.GalerkinNoBlowup

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
