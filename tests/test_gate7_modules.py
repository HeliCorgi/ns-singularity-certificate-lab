"""Tests for the Gate 7 modules: nondimensionalisation, families, ladder, slab."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.initial_data import (
    CompositeSwirlFamily,
    SwirlComponent,
    SwirlFamily,
    family_A,
    family_H,
    family_S,
)
from ns_certificate_lab.nondimensional import (
    ScalingExponents,
    Setting,
    deduplicate_settings,
    nonlinear_time_estimate,
    reynolds_number,
)
from ns_certificate_lab.picard_continuation import (
    DEPARTURE_THRESHOLDS,
    PREREGISTERED_REYNOLDS,
    PREREGISTERED_TAU,
    PicardLadder,
)
from ns_certificate_lab.slab_certificate import (
    HERMITE_BASIS_RANGES,
    SLAB_HYPOTHESES,
    build_slab_certificate,
    verify_slab_certificate,
)


def _grid(nr: int = 33, nz: int = 65) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=2.5, z_min=-3.0, z_max=3.0, periodic_z=False
    )


# --------------------------------------------------------------------------- #
# nondimensionalisation                                                        #
# --------------------------------------------------------------------------- #


def test_reynolds_number_matches_the_derivation() -> None:
    assert reynolds_number(2.0, 3.0, 0.5) == pytest.approx(2.0 * 9.0 / 0.5)


def test_reynolds_number_rejects_a_vanishing_viscosity() -> None:
    with pytest.raises(ValueError):
        reynolds_number(1.0, 1.0, 0.0)


def test_scaling_exponents_are_dimensionally_consistent() -> None:
    """Energy, the critical norm and the vorticity must be mutually consistent.

    ``E ~ A^2 L^5`` and ``|u| ~ A L`` over a volume ``~ L^3`` give
    ``E ~ (AL)^2 L^3``, and the critical norm ``(int |u|^3)^{1/3} ~ AL * L =
    A L^2``.  These are not independent choices, so a slip in one shows here.
    """
    e = ScalingExponents()
    a_v, b_v = e.velocity
    a_e, b_e = e.kinetic_energy
    assert (a_e, b_e) == (2 * a_v, 2 * b_v + 3)
    a_3, b_3 = e.l3_norm
    assert (a_3, b_3) == (a_v, b_v + 1)
    # A vorticity is a velocity over a length.
    a_w, b_w = e.vorticity
    assert (a_w, b_w) == (a_v, b_v - 1)


def test_scaling_factor_reproduces_the_exponents() -> None:
    e = ScalingExponents()
    assert e.factor("kinetic_energy", 2.0, 3.0) == pytest.approx(4.0 * 3.0**5)
    with pytest.raises(ValueError):
        e.factor("not_a_diagnostic", 1.0, 1.0)


def test_settings_with_equal_reynolds_shape_and_tau_are_duplicates() -> None:
    """Two raw settings can look different and be the same computation."""
    first = Setting(
        label="a", amplitude=1.0, length=2.0, viscosity=0.4,
        aspect_ratio=0.8, concentration=0.5, physical_time=0.1,
    )
    # Re = 1*4/0.4 = 10 and tau = 0.1.  Double the amplitude, halve the length
    # squared and quarter the time: Re = 2*2/0.4 = 10, tau = 2*0.05 = 0.1.
    second = Setting(
        label="b", amplitude=2.0, length=math.sqrt(2.0), viscosity=0.4,
        aspect_ratio=0.8, concentration=0.5, physical_time=0.05,
    )
    assert second.reynolds == pytest.approx(first.reynolds)
    assert second.dimensionless_time == pytest.approx(first.dimensionless_time)
    unique, duplicates = deduplicate_settings([first, second])
    assert len(unique) == 1
    assert duplicates == [("a", "b")]


def test_settings_differing_only_in_tau_are_not_duplicates() -> None:
    """A sweep run to a fixed physical time reaches different tau, and that
    difference is a real difference, not a duplicate."""
    base = dict(
        amplitude=1.0, length=2.0, viscosity=0.4, aspect_ratio=0.8,
        concentration=0.5,
    )
    early = Setting(label="early", physical_time=0.1, **base)
    late = Setting(label="late", physical_time=0.2, **base)
    unique, duplicates = deduplicate_settings([early, late])
    assert len(unique) == 2
    assert duplicates == []


def test_nonlinear_time_estimate_inverts_the_quadratic_law() -> None:
    # ratio = C tau^2 with C = 0.5: at tau = 0.2 the ratio is 0.02, and 10%
    # needs tau = sqrt(0.1/0.5) = 0.4472.
    tau = nonlinear_time_estimate(measured_ratio=0.02, measured_tau=0.2)
    assert tau == pytest.approx(math.sqrt(0.1 / 0.5))


def test_nonlinear_time_estimate_rejects_nonpositive_input() -> None:
    with pytest.raises(ValueError):
        nonlinear_time_estimate(measured_ratio=0.0, measured_tau=1.0)


# --------------------------------------------------------------------------- #
# the three preregistered families                                             #
# --------------------------------------------------------------------------- #


def test_family_S_reproduces_the_legacy_family_exactly() -> None:
    """Family S is the regression anchor; a drift here invalidates every
    earlier result, so the comparison is exact rather than approximate."""
    grid = _grid()
    legacy = SwirlFamily(
        amplitude=7.0, radial_support=1.2, axial_support=1.5,
        axial_concentration=0.5,
    )
    new = family_S(amplitude=7.0)
    r, z = grid.mesh()
    assert np.array_equal(new.u1(r, z), legacy.u1(r, z))


def test_moment_degeneracy_matches_the_analytic_claims() -> None:
    """S and H are even in ``u_1^2``, A is not."""
    grid = _grid()
    for builder, quadrupole_degenerate in (
        (family_S, True), (family_H, True), (family_A, False)
    ):
        report = builder(amplitude=3.0).moment_degeneracy(grid)
        assert report["monopole_degenerate"] is True
        assert report["transverse_quadrupole_degenerate"] is True
        assert report["axial_dipole_degenerate"] is False
        assert report["axial_dipole_leading"] == pytest.approx(-1.0)
        assert report["axial_quadrupole_degenerate"] is quadrupole_degenerate


def test_family_A_quadrupole_is_measurably_nonzero() -> None:
    """Non-degeneracy has to be a number, not a boolean: a quadrupole that is
    nonzero only at rounding level would not survive the elliptic solve."""
    grid = _grid()
    report = family_A(amplitude=3.0).moment_degeneracy(grid)
    assert abs(report["axial_quadrupole_leading"]) > 1.0e-3


def test_families_are_compactly_supported_and_vanish_on_the_axis() -> None:
    grid = _grid()
    r, z = grid.mesh()
    for builder in (family_S, family_A, family_H):
        family = builder(amplitude=2.0)
        values = family.u1(r, z)
        outside = np.hypot(r, z) > family.support_radius
        assert np.all(values[outside] == 0.0)
        # u^theta = r u_1 vanishes on the axis whatever u_1 does there.
        assert np.all(family.swirl_velocity(np.zeros_like(z[0]), z[0]) == 0.0)


def test_family_cartesian_velocity_is_a_pure_swirl() -> None:
    family = family_A(amplitude=2.0)
    points = np.array([[0.3, 0.4, 0.2], [-0.5, 0.1, -0.3], [0.0, 0.0, 0.1]])
    velocity = family.cartesian_velocity(points)
    # A pure swirl is orthogonal to the radial direction and has no axial part.
    assert np.allclose(velocity[:, 2], 0.0)
    radial = points[:, :2]
    assert np.allclose(np.einsum("ij,ij->i", velocity[:, :2], radial), 0.0)


def test_swirl_component_rejects_degenerate_parameters() -> None:
    with pytest.raises(ValueError):
        SwirlComponent(amplitude=1.0, radial_support=0.0, axial_support=1.0)
    with pytest.raises(ValueError):
        SwirlComponent(
            amplitude=1.0, radial_support=1.0, axial_support=1.0, radial_power=0.0
        )
    with pytest.raises(ValueError):
        CompositeSwirlFamily(name="empty", components=(), reference_length=1.0)


# --------------------------------------------------------------------------- #
# the Picard ladder                                                            #
# --------------------------------------------------------------------------- #


def test_preregistered_lists_are_ordered_and_positive() -> None:
    assert list(PREREGISTERED_TAU) == sorted(PREREGISTERED_TAU)
    assert list(PREREGISTERED_REYNOLDS) == sorted(PREREGISTERED_REYNOLDS)
    assert all(value > 0.0 for value in PREREGISTERED_TAU)
    assert DEPARTURE_THRESHOLDS["relative_distance_to_first_picard"] == 0.10


def test_ladder_rejects_periodic_grids_and_zero_viscosity() -> None:
    periodic = AxisymmetricGrid.uniform(nr=17, nz=17, r_max=1.0)
    with pytest.raises(ValueError):
        PicardLadder(periodic, viscosity=1.0e-2)
    with pytest.raises(ValueError):
        PicardLadder(_grid(17, 17), viscosity=0.0)


def test_level_zero_is_pure_diffusion_and_keeps_zero_vorticity() -> None:
    """Level 0 must not feel the nonlinearity at all: that is what makes the
    distance to level 1 a measurement of nonlinear response."""
    grid = _grid(25, 49)
    u1, _ = family_S(amplitude=4.0).initial_state(grid)
    ladder = PicardLadder(grid, viscosity=2.0e-2)
    _checkpoints, reason, state = ladder.run(
        u1, amplitude=float(np.max(np.abs(u1))), tau_checkpoints=(0.05,)
    )
    assert reason == "completed"
    assert np.max(np.abs(state[1])) == 0.0
    # Diffusion only decreases the peak of a positive-definite bump.
    assert np.max(np.abs(state[0])) <= np.max(np.abs(u1))


def test_full_solution_departs_from_the_first_picard_iterate_only_nonlinearly()\
        -> None:
    """The distance to level 1 must grow faster than linearly in ``tau``.

    If it grew linearly the ladder would be measuring a time-stepping artefact
    rather than nonlinear response.  A quadratic law gives a ratio of about
    four between successive doublings; anything below three would mean the
    signal is contaminated.
    """
    grid = _grid(33, 65)
    u1, _ = family_S(amplitude=8.0).initial_state(grid)
    amplitude = float(np.max(np.abs(u1)))
    ladder = PicardLadder(grid, viscosity=amplitude * 1.44 / 100.0)
    checkpoints, reason, _state = ladder.run(
        u1, amplitude=amplitude, tau_checkpoints=(0.1, 0.2, 0.4)
    )
    assert reason == "completed"
    distances = [c.distance_to_first_picard for c in checkpoints]
    assert distances[0] < distances[1] < distances[2]
    assert distances[1] / distances[0] > 3.0
    assert distances[2] / distances[1] > 3.0


def test_second_picard_iterate_is_closer_than_the_first() -> None:
    """A ladder that is working has the iterates converging toward the full
    solution; if level 2 were further away the recursion would be wrong."""
    grid = _grid(33, 65)
    u1, _ = family_S(amplitude=8.0).initial_state(grid)
    amplitude = float(np.max(np.abs(u1)))
    ladder = PicardLadder(grid, viscosity=amplitude * 1.44 / 100.0)
    checkpoints, _reason, _state = ladder.run(
        u1, amplitude=amplitude, tau_checkpoints=(0.2,)
    )
    final = checkpoints[-1]
    assert final.distance_to_second_picard < final.distance_to_first_picard


def test_stable_time_step_respects_the_diffusion_limit() -> None:
    grid = _grid(33, 65)
    u1, _ = family_S(amplitude=1.0).initial_state(grid)
    zero = np.zeros_like(u1)
    state = (u1, zero) * 4
    for viscosity in (1.0e-2, 1.0e-1):
        ladder = PicardLadder(grid, viscosity=viscosity)
        limit = 0.2 / (viscosity * (8.0 / grid.dr**2 + 2.0 / grid.dz**2))
        assert ladder.stable_time_step(state) <= limit * (1.0 + 1.0e-12)


def test_integrators_agree_on_a_short_interval() -> None:
    grid = _grid(33, 65)
    u1, _ = family_S(amplitude=8.0).initial_state(grid)
    amplitude = float(np.max(np.abs(u1)))
    results = {}
    for integrator in ("rk4", "ssprk3"):
        ladder = PicardLadder(grid, viscosity=amplitude * 1.44 / 100.0)
        checkpoints, _reason, _state = ladder.run(
            u1, amplitude=amplitude, tau_checkpoints=(0.2,), integrator=integrator
        )
        results[integrator] = checkpoints[-1].max_omega1
    assert results["rk4"] == pytest.approx(results["ssprk3"], rel=0.02)


def test_unknown_integrator_is_rejected() -> None:
    grid = _grid(17, 17)
    u1, _ = family_S(amplitude=1.0).initial_state(grid)
    zero = np.zeros_like(u1)
    ladder = PicardLadder(grid, viscosity=1.0e-2)
    with pytest.raises(ValueError):
        ladder.step((u1, zero) * 4, 1.0e-4, "euler")


def test_checkpoint_serialises_without_non_finite_values() -> None:
    grid = _grid(25, 49)
    u1, _ = family_S(amplitude=4.0).initial_state(grid)
    ladder = PicardLadder(grid, viscosity=2.0e-2)
    checkpoints, _reason, _state = ladder.run(
        u1, amplitude=float(np.max(np.abs(u1))), tau_checkpoints=(0.05,)
    )
    json.dumps(checkpoints[-1].as_dict(), allow_nan=False)


# --------------------------------------------------------------------------- #
# the space-time slab certificate                                              #
# --------------------------------------------------------------------------- #


def test_hermite_basis_ranges_are_the_exact_extrema() -> None:
    """These four ranges are what makes the temporal enclosure exact, so they
    are checked against the polynomials themselves rather than trusted."""
    theta = np.linspace(0.0, 1.0, 200001)
    polynomials = {
        "h00": 2 * theta**3 - 3 * theta**2 + 1,
        "h01": -2 * theta**3 + 3 * theta**2,
        "h10": theta**3 - 2 * theta**2 + theta,
        "h11": theta**3 - theta**2,
    }
    for name, values in polynomials.items():
        lower, upper = HERMITE_BASIS_RANGES[name]
        assert float(lower) <= values.min() + 1.0e-12
        assert values.max() - 1.0e-12 <= float(upper)
        # And the ranges are tight, not merely valid.
        assert values.min() == pytest.approx(float(lower), abs=1.0e-9)
        assert values.max() == pytest.approx(float(upper), abs=1.0e-9)


def test_hermite_value_bases_form_a_partition_of_unity() -> None:
    theta = np.linspace(0.0, 1.0, 1001)
    total = (2 * theta**3 - 3 * theta**2 + 1) + (-2 * theta**3 + 3 * theta**2)
    assert np.allclose(total, 1.0)


def _slab_inputs():
    grid = _grid(21, 41)
    ladder = PicardLadder(grid, viscosity=5.0e-2)
    u1, _ = family_S(amplitude=4.0).initial_state(grid)
    amplitude = float(np.max(np.abs(u1)))
    _checkpoints, reason, state = ladder.run(
        u1, amplitude=amplitude, tau_checkpoints=(0.05,)
    )
    assert reason == "completed"
    dt = ladder.stable_time_step(state)
    rhs_start = ladder.rhs(state)
    mid = ladder.step(state, 0.5 * dt, "rk4")
    rhs_mid = ladder.rhs(mid)
    end = ladder.step(state, dt, "rk4")
    rhs_end = ladder.rhs(end)
    return grid, ladder, state, end, dt, rhs_start, rhs_mid, rhs_end


def _build(grid, ladder, state, end, dt, rhs_start, rhs_mid, rhs_end, **kwargs):
    return build_slab_certificate(
        grid,
        time_start=0.5,
        time_step=dt,
        u1_start=state[6], u1_end=end[6],
        omega1_start=state[7], omega1_end=end[7],
        psi1_start=ladder.solver.solve(state[7]),
        psi1_end=ladder.solver.solve(end[7]),
        du1_start=rhs_start[6], du1_end=rhs_end[6],
        domega1_start=rhs_start[7], domega1_end=rhs_end[7],
        domega1_mid=rhs_mid[7],
        viscosity=ladder.viscosity,
        green_tail_bound=1.0e-6,
        interior_radius=1.0,
        **kwargs,
    )


def test_slab_certificate_verifies_and_encloses_cell_interiors() -> None:
    certificate = _build(*_slab_inputs())
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    assert payload["slab"]["encloses_cell_interiors"] is True
    assert payload["slab"]["cells"] > 0
    assert verify_slab_certificate(payload)["verified"]


def test_slab_certificate_records_both_hypotheses_as_unproved() -> None:
    """The point of the module is that these two are *not* theorems, so a
    payload that quietly drops or promotes one must fail the checker."""
    certificate = _build(*_slab_inputs())
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    assert set(payload["hypotheses"]) == set(SLAB_HYPOTHESES)
    assert all(
        entry["proved"] is False for entry in payload["hypotheses"].values()
    )


def test_simpson_defect_is_far_below_the_trapezoid_defect() -> None:
    """Simpson is two orders higher, so on a smooth step its defect must be
    much smaller.  Equality would mean the midpoint field is not being used."""
    from fractions import Fraction

    certificate = _build(*_slab_inputs())
    bounds = certificate.as_dict()["bounds"]
    trapezoid = Fraction(str(bounds["trapezoid_local_defect"]))
    simpson = Fraction(str(bounds["simpson_local_defect"]))
    assert 0 < simpson < trapezoid / 100


def test_slab_checker_rejects_a_narrowed_enclosure() -> None:
    certificate = _build(*_slab_inputs())
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    lower, upper = payload["enclosures"]["omega1"]
    payload["enclosures"]["omega1"] = [upper, lower]
    assert not verify_slab_certificate(payload)["verified"]


def test_slab_checker_rejects_a_forged_proof_claim() -> None:
    certificate = _build(*_slab_inputs())
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    payload["hypotheses"]["H2_hermite_remainder"]["proved"] = True
    verdict = verify_slab_certificate(payload)
    assert not verdict["verified"]
    assert any("proved" in failure for failure in verdict["failures"])


def test_slab_checker_rejects_an_inconsistent_remainder() -> None:
    """The remainder must equal ``M4 Delta^4 / 384``; a payload that inflates
    ``M4`` while leaving the remainder alone is claiming a tighter tube than
    its own hypothesis supports."""
    certificate = _build(*_slab_inputs())
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    payload["hypotheses"]["H2_hermite_remainder"]["fourth_derivative_bound"] = "1"
    assert not verify_slab_certificate(payload)["verified"]


def test_slab_checker_rejects_a_dropped_clay_disclaimer() -> None:
    certificate = _build(*_slab_inputs())
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    payload["claims"] = ["everything is fine"]
    assert not verify_slab_certificate(payload)["verified"]


def test_slab_checker_returns_a_verdict_on_garbage_rather_than_raising() -> None:
    for payload in ({}, {"certificate": "space_time_slab"}, {"slab": 3}):
        verdict = verify_slab_certificate(payload)
        assert verdict["verified"] is False
        assert verdict["failures"]


def test_slab_builder_refuses_a_periodic_grid() -> None:
    periodic = AxisymmetricGrid.uniform(nr=9, nz=9, r_max=1.0)
    zero = np.zeros((9, 9))
    with pytest.raises(ValueError):
        build_slab_certificate(
            periodic, time_start=0.0, time_step=1.0e-3,
            u1_start=zero, u1_end=zero, omega1_start=zero, omega1_end=zero,
            psi1_start=zero, psi1_end=zero, du1_start=zero, du1_end=zero,
            domega1_start=zero, domega1_end=zero,
            viscosity=1.0e-2, green_tail_bound=0.0,
        )


def test_slab_builder_refuses_a_vanishing_viscosity() -> None:
    grid = _grid(9, 9)
    zero = np.zeros((9, 9))
    with pytest.raises(ValueError):
        build_slab_certificate(
            grid, time_start=0.0, time_step=1.0e-3,
            u1_start=zero, u1_end=zero, omega1_start=zero, omega1_end=zero,
            psi1_start=zero, psi1_end=zero, du1_start=zero, du1_end=zero,
            domega1_start=zero, domega1_end=zero,
            viscosity=0.0, green_tail_bound=0.0,
        )


# --------------------------------------------------------------------------- #
# the Gate 7 experiment: reclassification, gate, promotion                     #
# --------------------------------------------------------------------------- #


def _config(tmp_path):
    from pathlib import Path

    source = Path("configs/tau_continuation_gate7.json")
    return json.loads(source.read_text(encoding="utf-8"))


def test_legacy_sweep_reclassification_reports_the_tau_actually_reached(
    tmp_path,
) -> None:
    """The Gate 6 sweep's own record must show it never left the linear regime.

    This is the finding the reclassification exists to make checkable: the
    sweep's nominal amplitudes overstated the field's peak by a factor of about
    twenty-three, so the ``tau`` it reached was two orders below one.
    """
    from pathlib import Path

    from experiments.run_tau_continuation import reclassify_legacy_sweep

    config = _config(tmp_path)
    report = reclassify_legacy_sweep(config, Path("."))
    if not report["available"]:
        pytest.skip("the Gate 6 sweep output is not present")
    assert report["points"] == 32
    assert report["max_dimensionless_time"] < 0.05
    assert report["min_dimensionless_time"] > 0.0
    # Every effective amplitude must be well below its nominal parameter.
    for row in report["detail"]:
        assert row["effective_amplitude"] < 0.2 * row["nominal_amplitude"]


def test_departure_gate_fails_a_run_still_inside_the_first_picard_iterate(
    tmp_path,
) -> None:
    from experiments.run_tau_continuation import evaluate_departure
    from ns_certificate_lab.picard_continuation import (
        ContinuationCheckpoint,
        ContinuationResult,
    )

    checkpoint = ContinuationCheckpoint(
        tau=1.0, physical_time=1.0, steps=10, accepted_time_step=0.1,
        max_u1=1.0, max_omega1=1.0, max_speed=1.0,
        distance_to_first_picard=1.0e-4,  # still linear
        distance_to_second_picard=1.0e-4,
        advection_over_source=1.0e-4, viscous_over_source=0.5,
        rhs_change_rate=0.0, kinetic_energy=1.0, enstrophy=1.0, l3_norm=1.0,
        core_l3_fraction=0.5, shell_fractions=(0.5, 0.5),
        effective_shell_count=2.0, radial_scale=1.0, axial_scale=1.0,
        shell_flux_upward=0.5, outer_band_fraction=1.0e-30,
        relative_divergence=1.0e-3, curl_defect=1.0e-3, poisson_residual=1.0e-15,
        front_points_radial=20.0, front_points_axial=20.0,
    )
    result = ContinuationResult(
        label="x", family="S", reynolds=100.0, amplitude=1.0,
        reference_length=1.0, viscosity=1.0e-2, nr=33, nz=65, r_max=2.0,
        z_max=2.0, integrator="rk4", checkpoints=(checkpoint,), reached_tau=1.0,
        stopped_reason="completed",
    )
    config = _config(tmp_path)
    gate = evaluate_departure(
        result, coarse=None, integrator_check=None, time_step_check=None,
        thresholds=config["preregistration"]["departure_gate"],
    )
    assert gate["passed"] is False
    assert "distance_to_first_picard" in gate["failed"]
    assert "advection_measurable" in gate["failed"]


def test_promotion_rejects_a_decaying_run_with_a_stated_reason(tmp_path) -> None:
    """A rejection has to say why, in numbers.  A bare ``False`` would make an
    empty result indistinguishable from a run that was never scored."""
    from experiments.run_tau_continuation import evaluate_promotion
    from ns_certificate_lab.picard_continuation import (
        ContinuationCheckpoint,
        ContinuationResult,
    )

    def make(tau: float, l3: float, shells: float, width: float):
        return ContinuationCheckpoint(
            tau=tau, physical_time=tau, steps=1, accepted_time_step=0.1,
            max_u1=1.0, max_omega1=1.0, max_speed=1.0,
            distance_to_first_picard=0.3, distance_to_second_picard=0.2,
            advection_over_source=0.4, viscous_over_source=0.3,
            rhs_change_rate=0.0, kinetic_energy=1.0, enstrophy=1.0, l3_norm=l3,
            core_l3_fraction=0.5, shell_fractions=(0.5, 0.5),
            effective_shell_count=shells, radial_scale=width, axial_scale=width,
            shell_flux_upward=0.5, outer_band_fraction=1.0e-30,
            relative_divergence=1.0e-3, curl_defect=1.0e-3,
            poisson_residual=1.0e-15, front_points_radial=20.0,
            front_points_axial=20.0,
        )

    decaying = ContinuationResult(
        label="decay", family="S", reynolds=400.0, amplitude=1.0,
        reference_length=1.0, viscosity=1.0e-2, nr=33, nz=65, r_max=2.0,
        z_max=2.0, integrator="rk4",
        checkpoints=(make(0.1, 1.0, 2.0, 1.0), make(1.0, 0.75, 2.1, 1.05)),
        reached_tau=1.0, stopped_reason="completed",
    )
    config = _config(tmp_path)
    verdict = evaluate_promotion(decaying, config["preregistration"]["promotion"])
    assert verdict["promoted"] is False
    assert len(verdict["rejection_reasons"]) == 2
    assert "0.7500" in verdict["rejection_reasons"][0]


def test_promotion_accepts_a_concentrating_run(tmp_path) -> None:
    """The criteria must be reachable: a criterion nothing can satisfy is not a
    test, and would silently convert every run into a rejection."""
    from experiments.run_tau_continuation import evaluate_promotion
    from ns_certificate_lab.picard_continuation import (
        ContinuationCheckpoint,
        ContinuationResult,
    )

    def make(tau: float, l3: float, shells: float, width: float):
        return ContinuationCheckpoint(
            tau=tau, physical_time=tau, steps=1, accepted_time_step=0.1,
            max_u1=1.0, max_omega1=1.0, max_speed=1.0,
            distance_to_first_picard=0.3, distance_to_second_picard=0.2,
            advection_over_source=0.4, viscous_over_source=0.3,
            rhs_change_rate=0.0, kinetic_energy=1.0, enstrophy=1.0, l3_norm=l3,
            core_l3_fraction=0.5, shell_fractions=(0.5, 0.5),
            effective_shell_count=shells, radial_scale=width, axial_scale=width,
            shell_flux_upward=0.5, outer_band_fraction=1.0e-30,
            relative_divergence=1.0e-3, curl_defect=1.0e-3,
            poisson_residual=1.0e-15, front_points_radial=20.0,
            front_points_axial=20.0,
        )

    concentrating = ContinuationResult(
        label="concentrate", family="S", reynolds=400.0, amplitude=1.0,
        reference_length=1.0, viscosity=1.0e-2, nr=33, nz=65, r_max=2.0,
        z_max=2.0, integrator="rk4",
        checkpoints=(make(0.1, 1.0, 2.0, 1.0), make(1.0, 1.2, 1.8, 0.8)),
        reached_tau=1.0, stopped_reason="completed",
    )
    config = _config(tmp_path)
    verdict = evaluate_promotion(
        concentrating, config["preregistration"]["promotion"]
    )
    assert verdict["promoted"] is True
    assert verdict["critical_concentration"] is True
