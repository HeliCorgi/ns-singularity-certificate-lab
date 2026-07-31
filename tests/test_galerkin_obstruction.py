"""Tests for the Track-F finite-mode obstruction.

The exact-arithmetic tests are the theorem-side objects: they check integer
identities with no floating point involved.  The floating-point tests are
cross-checks of the *proved* bound, plus the fault injections that AGENTS.md
requires every defect-detection claim to carry.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from ns_certificate_lab.galerkin_obstruction import (
    apriori_norm_bound,
    build_galerkin_system,
    build_mode_set,
    derivative_amplification,
    finite_mode_no_go_certificate,
    sobolev_amplification,
    stream_apriori_bound,
    suggested_time_step,
    transverse_integer_basis,
    verify_trilinear_cancellation,
)

# Mode sets with genuine resonances (k+l+m=0 with all three in S), so the
# cancellation test has something to cancel.
RESONANT_SEEDS = {
    "planar_triad": [(1, 0, 0), (0, 1, 0), (1, 1, 0)],
    "oblique_triad": [(1, 1, 0), (1, 0, 1), (0, 1, 1), (1, -1, 0)],
    "ball_two": [
        k
        for k in itertools.product((-1, 0, 1), repeat=3)
        if 1 <= sum(c * c for c in k) <= 2
    ],
    "with_zero_mode": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)],
    "anisotropic": [(2, 0, 0), (0, 1, 0), (2, 1, 0), (1, 1, 1), (3, 1, 1)],
}


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# --------------------------------------------------------------------------- #
# Exact arithmetic: the integer geometry                                        #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "k",
    [
        (1, 0, 0),
        (0, 0, 1),
        (1, 1, 0),
        (1, -1, 2),
        (3, 5, 7),
        (2, -4, 6),
        (-1, 0, 0),
        (12, 18, 30),
    ],
)
def test_transverse_basis_is_exactly_orthogonal_and_independent(k):
    frame = transverse_integer_basis(k)
    assert len(frame) == 2
    for t in frame:
        assert _dot(k, t) == 0  # exact integer identity
    cross = (
        frame[0][1] * frame[1][2] - frame[0][2] * frame[1][1],
        frame[0][2] * frame[1][0] - frame[0][0] * frame[1][2],
        frame[0][0] * frame[1][1] - frame[0][1] * frame[1][0],
    )
    assert cross != (0, 0, 0)


def test_zero_mode_has_no_constraint():
    assert transverse_integer_basis((0, 0, 0)) == ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def test_mode_set_is_symmetric_and_dimension_counts_correctly():
    mode_set = build_mode_set([(1, 0, 0), (0, 0, 0), (2, 3, 4)])
    assert set(mode_set.modes) == {
        (1, 0, 0),
        (-1, 0, 0),
        (0, 0, 0),
        (2, 3, 4),
        (-2, -3, -4),
    }
    # two +/- pairs contribute 4 real degrees of freedom each, the zero mode 3.
    assert mode_set.dimension == 4 + 4 + 3
    assert mode_set.transverse_defect() == 0


def test_seeds_must_be_three_dimensional():
    with pytest.raises(ValueError):
        build_mode_set([(1, 0)])
    with pytest.raises(ValueError):
        build_mode_set([])


# --------------------------------------------------------------------------- #
# Exact arithmetic: Lemma 1                                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("label", sorted(RESONANT_SEEDS))
def test_trilinear_form_cancels_exactly(label):
    """``⟨u,(u·∇)u⟩ ≡ 0``: every monomial coefficient vanishes over ``ℤ[i]``."""
    certificate = verify_trilinear_cancellation(build_mode_set(RESONANT_SEEDS[label]))
    assert certificate.transverse_defect == 0
    assert certificate.surviving_monomials == 0
    assert certificate.max_abs_surviving_coefficient == 0
    assert certificate.cancels


@pytest.mark.parametrize("label", sorted(RESONANT_SEEDS))
def test_trilinear_certificate_is_not_vacuous(label):
    """Guard against a certificate that passes because nothing was expanded."""
    certificate = verify_trilinear_cancellation(build_mode_set(RESONANT_SEEDS[label]))
    assert certificate.resonant_triples > 0
    assert certificate.monomials_accumulated > 0


def test_nonresonant_mode_set_is_reported_as_having_no_triples():
    """A set with no resonances cancels trivially and must say so."""
    certificate = verify_trilinear_cancellation(build_mode_set([(1, 0, 0)]))
    assert certificate.resonant_triples == 0
    assert certificate.monomials_accumulated == 0
    assert certificate.cancels


@pytest.mark.parametrize("label", sorted(RESONANT_SEEDS))
def test_fault_injection_longitudinal_component_breaks_cancellation(label):
    """Dropping the divergence-free constraint must break the exact identity."""
    faulty = build_mode_set(RESONANT_SEEDS[label], allow_longitudinal=True)
    certificate = verify_trilinear_cancellation(faulty)
    assert certificate.transverse_defect > 0
    assert certificate.surviving_monomials > 0
    assert certificate.max_abs_surviving_coefficient > 0
    assert not certificate.cancels


def test_fault_detection_is_independent_of_the_defect_counter():
    """The monomial check alone (not the ``k·t`` audit) detects the fault."""
    faulty = build_mode_set(RESONANT_SEEDS["ball_two"], allow_longitudinal=True)
    certificate = verify_trilinear_cancellation(faulty)
    assert certificate.surviving_monomials >= 100


# --------------------------------------------------------------------------- #
# The proved constants                                                          #
# --------------------------------------------------------------------------- #


def test_apriori_bound_is_the_stated_formula():
    assert apriori_norm_bound(2.0, 3.0) == pytest.approx(5.0)
    with pytest.raises(ValueError):
        apriori_norm_bound(-1.0, 0.0)
    with pytest.raises(ValueError):
        apriori_norm_bound(0.0, -1.0)


def test_sobolev_amplification_matches_the_sharp_constant():
    mode_set = build_mode_set([(1, 0, 0), (2, 1, 0)])
    radius_squared = 5.0
    assert sobolev_amplification(mode_set, 0.0) == pytest.approx(1.0)
    assert sobolev_amplification(mode_set, 2.0) == pytest.approx(
        1.0 + 4.0 * math.pi**2 * radius_squared
    )
    with pytest.raises(ValueError):
        sobolev_amplification(mode_set, -1.0)


def test_derivative_amplification_matches_the_stated_bound():
    mode_set = build_mode_set([(1, 0, 0), (2, 1, 0)])
    assert derivative_amplification(mode_set, 0) == pytest.approx(
        math.sqrt(len(mode_set.modes))
    )
    assert derivative_amplification(mode_set, 2) == pytest.approx(
        (2.0 * math.pi * math.sqrt(5.0)) ** 2 * math.sqrt(len(mode_set.modes))
    )
    with pytest.raises(ValueError):
        derivative_amplification(mode_set, -1)


def test_derivative_amplification_dominates_a_sampled_supremum():
    """The proved ``L^∞`` constant must not be violated by an actual field."""
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"])
    system = build_galerkin_system(mode_set, 0.0)
    rng = np.random.default_rng(11)
    c = rng.standard_normal(system.dimension)
    resolution = 9
    axis = np.arange(resolution) / resolution
    points = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    ).reshape(-1, 3)
    sampled = float(np.abs(system.velocity(c, points).real).max())
    assert sampled <= derivative_amplification(mode_set, 0) * float(
        np.linalg.norm(c)
    )


# --------------------------------------------------------------------------- #
# The floating-point Galerkin system (cross-check)                              #
# --------------------------------------------------------------------------- #


def _unit_grid(resolution: int) -> np.ndarray:
    axis = np.arange(resolution) / resolution
    return np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(
        -1, 3
    )


def test_reconstructed_velocity_is_real_and_divergence_free():
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"])
    system = build_galerkin_system(mode_set, 1.0e-3)
    rng = np.random.default_rng(20260729)
    c = rng.standard_normal(system.dimension)
    field = system.velocity(c, _unit_grid(8))
    assert float(np.abs(field.imag).max()) < 1.0e-12
    coefficients = system.coefficients_to_fourier(c)
    for k in mode_set.modes:
        longitudinal = complex(coefficients[k] @ np.asarray(k, dtype=np.float64))
        assert abs(longitudinal) < 1.0e-12


def test_real_basis_is_l2_orthonormal():
    """``‖u‖_{L²(𝕋³)} = ‖c‖₂``, checked against a quadrature on a uniform grid."""
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"])
    system = build_galerkin_system(mode_set, 0.0)
    rng = np.random.default_rng(5)
    c = rng.standard_normal(system.dimension)
    resolution = 8
    field = system.velocity(c, _unit_grid(resolution)).real
    quadrature = math.sqrt(float((field**2).sum()) / resolution**3)
    assert quadrature == pytest.approx(float(np.linalg.norm(c)), rel=1.0e-12)


def test_projection_round_trip_is_the_identity_on_the_ansatz_space():
    mode_set = build_mode_set(RESONANT_SEEDS["oblique_triad"])
    system = build_galerkin_system(mode_set, 0.0)
    rng = np.random.default_rng(7)
    c = rng.standard_normal(system.dimension)
    recovered = system.fourier_to_coefficients(system.coefficients_to_fourier(c))
    assert np.allclose(recovered, c, rtol=0.0, atol=1.0e-12)


@pytest.mark.parametrize("label", sorted(RESONANT_SEEDS))
def test_energy_production_vanishes_to_roundoff(label):
    """Floating-point shadow of Lemma 1, on the same mode sets."""
    mode_set = build_mode_set(RESONANT_SEEDS[label])
    system = build_galerkin_system(mode_set, 0.0)
    rng = np.random.default_rng(3)
    for _ in range(5):
        c = rng.standard_normal(system.dimension)
        scale = float(np.linalg.norm(c)) ** 3
        assert abs(system.energy_production(c)) / scale < 1.0e-12


def test_inviscid_unforced_trajectory_conserves_energy_and_does_not_blow_up():
    """The would-be blow-up case: no viscosity, no force, large amplitude."""
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"])
    system = build_galerkin_system(mode_set, 0.0)
    rng = np.random.default_rng(20260729)
    c = rng.standard_normal(system.dimension)
    amplitude = 200.0
    c *= amplitude / float(np.linalg.norm(c))
    dt = suggested_time_step(system, amplitude)
    report = stream_apriori_bound(system, c, final_time=2000 * dt, steps=2000)
    assert not report.diverged
    assert report.bound_respected
    assert report.max_relative_energy_production < 1.0e-12
    assert report.final_norm == pytest.approx(amplitude, rel=1.0e-6)


def test_forced_trajectory_saturates_but_never_exceeds_the_proved_bound():
    """A force aligned with the state drives the trajectory onto the bound."""
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"])
    system = build_galerkin_system(mode_set, 1.0e-3)
    rng = np.random.default_rng(101)
    direction = rng.standard_normal(system.dimension)
    direction /= float(np.linalg.norm(direction))
    magnitude = 4.0
    report = stream_apriori_bound(
        system,
        0.01 * direction,
        final_time=3.0,
        steps=3000,
        force=lambda _t: magnitude * direction,
        force_l2_sup=magnitude,
    )
    assert report.bound_respected
    # Non-vacuity: the trajectory comes within 1% of the proved upper bound.
    assert report.max_bound_ratio > 0.99


def test_fault_injection_longitudinal_frame_breaks_the_energy_bound():
    """Removing the divergence-free constraint must break the proved bound."""
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"], allow_longitudinal=True)
    system = build_galerkin_system(mode_set, 0.0)
    rng = np.random.default_rng(20260729)
    c = rng.standard_normal(system.dimension)
    amplitude = 200.0
    c *= amplitude / float(np.linalg.norm(c))
    dt = suggested_time_step(system, amplitude)
    report = stream_apriori_bound(system, c, final_time=2000 * dt, steps=2000)
    assert not report.bound_respected
    # The step-size-independent detector separates fault from truncation error
    # by more than ten orders of magnitude.
    assert report.max_relative_energy_production > 1.0e-3


def test_fault_injection_viscosity_sign_breaks_the_energy_bound():
    """A negative viscosity violates ``Dissipative`` and must be caught."""
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"])
    with pytest.raises(ValueError):
        build_galerkin_system(mode_set, -1.0e-3)
    # Build the anti-dissipative operator by hand and confirm the growth.
    system = build_galerkin_system(mode_set, 1.0e-2)
    antidissipative = build_galerkin_system(mode_set, 1.0e-2)
    object.__setattr__(antidissipative, "diagonal", -system.diagonal)
    rng = np.random.default_rng(13)
    c = 0.1 * rng.standard_normal(system.dimension)
    report = stream_apriori_bound(antidissipative, c, final_time=5.0, steps=5000)
    assert not report.bound_respected


def test_stream_validates_its_inputs():
    mode_set = build_mode_set(RESONANT_SEEDS["planar_triad"])
    system = build_galerkin_system(mode_set, 0.0)
    c = np.zeros(system.dimension)
    with pytest.raises(ValueError):
        stream_apriori_bound(system, c, final_time=1.0, steps=0)
    with pytest.raises(ValueError):
        stream_apriori_bound(system, c, final_time=0.0, steps=10)
    with pytest.raises(ValueError):
        stream_apriori_bound(system, np.zeros(system.dimension + 1), final_time=1.0, steps=10)
    with pytest.raises(ValueError):
        suggested_time_step(system, -1.0)
    with pytest.raises(ValueError):
        suggested_time_step(system, 1.0, safety=0.0)


def test_suggested_time_step_shrinks_with_amplitude():
    mode_set = build_mode_set(RESONANT_SEEDS["ball_two"])
    system = build_galerkin_system(mode_set, 1.0e-3)
    assert suggested_time_step(system, 100.0) < suggested_time_step(system, 1.0)
    assert math.isinf(suggested_time_step(build_galerkin_system(mode_set, 0.0), 0.0))


# --------------------------------------------------------------------------- #
# The packaged verdict                                                          #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("label", sorted(RESONANT_SEEDS))
def test_no_go_certificate_rejects_every_finite_mode_class(label):
    certificate = finite_mode_no_go_certificate(label, RESONANT_SEEDS[label])
    assert certificate.verdict == "rejected_as_clay_cd_candidate"
    payload = certificate.as_dict()
    assert payload["label"] == label
    assert payload["dimension"] == certificate.dimension
    assert payload["trilinear"]["cancels"] is True
    assert certificate.sobolev_amplification_h3 > certificate.sobolev_amplification_h1


def test_no_go_certificate_reports_failure_under_fault_injection():
    certificate = finite_mode_no_go_certificate(
        "faulty", RESONANT_SEEDS["ball_two"], allow_longitudinal=True
    )
    assert certificate.verdict == "certificate_failed"
    assert certificate.as_dict()["verdict"] == "certificate_failed"
