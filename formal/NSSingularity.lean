/-
Copyright (c) 2026 ns-singularity-certificate-lab contributors.
Released under Apache 2.0 license as described in the file LICENSE.
-/
-- Root module for the NSSingularity formalization.
-- Stage 0 (LEAN4_VERIFICATION_POLICY.md): fix the Clay statements.
import NSSingularity.ClayStatement
import NSSingularity.VelocityRecovery
import NSSingularity.FiniteTime
import NSSingularity.GalerkinNoBlowup
import NSSingularity.FiniteModeNoGo
import NSSingularity.GreenAndCascade
import NSSingularity.CertificateLayer
import NSSingularity.TimeDependentGalerkin
import NSSingularity.L3Generation
import NSSingularity.ControlODE
import NSSingularity.TrackPFourier
import NSSingularity.GaussianTransfer
import NSSingularity.TrackPChain
import NSSingularity.GalerkinPicard
import NSSingularity.KatoConstant
import NSSingularity.ChainAnalysis
import NSSingularity.MesoscopicDuhamelNoGo
import NSSingularity.SpectralFrontIdentities
import NSSingularity.BandSymmetry
