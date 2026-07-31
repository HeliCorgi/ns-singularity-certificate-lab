"""Tests for the small-amplitude whole-space nonlinear evolution and the toy cascade."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ns_certificate_lab.cascade_toy import CascadeModel, decompose_forcing_roles
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.whole_space_evolution import (
    INTEGRATORS,
    SwirlBump,
    WholeSpaceEvolution,
    cartesian_audit_of_state,
    initial_state,
)

BUMP = SwirlBump(center=0.0, radius=1.5, amplitude=0.05)


def _grid(nr: int, nz: int, r_max: float = 4.0, z_max: float = 4.0) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=r_max, z_min=-z_max, z_max=z_max, periodic_z=False
    )


def _run(nr=65, nz=129, r_max=4.0, z_max=4.0, steps=20, dt=2.5e-3, integrator="rk4",
         boundary_mode="monopole"):
    grid = _grid(nr, nz, r_max, z_max)
    evolution = WholeSpaceEvolution(grid, viscosity=5.0e-3, boundary_mode=boundary_mode)
    u1, omega1 = initial_state(grid, BUMP)
    return grid, evolution, evolution.run(
        u1, omega1, steps=steps, time_step=dt, integrator=integrator,
        snapshot_every=steps,
    )


def _order(coarse: float, fine: float) -> float:
    return math.log(coarse / fine) / math.log(2.0)


# --------------------------------------------------------------------------- #
# Initial data                                                                 #
# --------------------------------------------------------------------------- #


def test_initial_data_is_pure_swirl_and_compactly_supported():
    grid = _grid(65, 129)
    u1, omega1 = initial_state(grid, BUMP)
    assert np.all(omega1 == 0.0)
    r_mesh, z_mesh = grid.mesh()
    outside = np.hypot(r_mesh, z_mesh) > BUMP.radius
    assert np.all(u1[outside] == 0.0)
    assert float(np.max(u1)) == pytest.approx(BUMP.amplitude, rel=1.0e-6)
    # Boundary pinning is consistent with compact support.
    assert np.all(u1[-1, :] == 0.0) and np.all(u1[:, 0] == 0.0)


def test_swirl_bump_rejects_bad_specifications():
    with pytest.raises(ValueError):
        SwirlBump(radius=0.0)
    with pytest.raises(ValueError):
        SwirlBump(power=1)


def test_evolution_rejects_a_periodic_grid_and_bad_settings():
    periodic = AxisymmetricGrid.uniform(
        nr=17, nz=32, r_max=4.0, z_min=-4.0, z_max=4.0, periodic_z=True
    )
    with pytest.raises(ValueError):
        WholeSpaceEvolution(periodic, viscosity=1.0e-3)
    with pytest.raises(ValueError):
        WholeSpaceEvolution(_grid(33, 65), viscosity=-1.0)
    with pytest.raises(ValueError):
        WholeSpaceEvolution(_grid(33, 65), viscosity=1.0e-3, outer_band_cells=0)
    grid, evolution, _ = _run(nr=33, nz=65, steps=1)
    u1, omega1 = initial_state(grid, BUMP)
    with pytest.raises(ValueError):
        evolution.step(u1, omega1, 1.0e-3, "midpoint")
    with pytest.raises(ValueError):
        evolution.run(u1, omega1, steps=0, time_step=1.0e-3)


# --------------------------------------------------------------------------- #
# The mechanism                                                                #
# --------------------------------------------------------------------------- #


def test_vorticity_is_generated_only_by_the_audited_stretching_term():
    """From pure swirl, ``omega_1`` is produced by ``d_z(u_1^2)`` alone."""
    grid, evolution, _ = _run(steps=1)
    u1, omega1 = initial_state(grid, BUMP)
    _du, dw = evolution.rhs(u1, omega1)
    from ns_certificate_lab.operators import derivative_z

    expected = derivative_z(grid, u1 * u1)
    expected[-1, :] = 0.0
    expected[:, 0] = 0.0
    expected[:, -1] = 0.0
    assert np.allclose(dw, expected, rtol=0.0, atol=1.0e-15)


def test_energy_never_increases_and_the_run_stays_finite():
    _grid_, _evolution, summary = _run(steps=40)
    assert not summary.diverged
    assert summary.streaming_max_energy_increase <= 1.0e-12
    energies = [snapshot.kinetic_energy for snapshot in summary.snapshots]
    assert energies[0] > energies[-1] > 0.0
    assert summary.final.max_omega1 > 0.0


def test_streaming_diagnostics_are_populated():
    _grid_, _evolution, summary = _run(steps=20)
    payload = summary.as_dict()
    for key in (
        "streaming_max_relative_divergence",
        "streaming_max_curl_defect",
        "streaming_max_poisson_residual",
        "streaming_max_outer_band_fraction",
    ):
        assert key in payload
    final = summary.final
    assert final.poisson_residual_max < 1.0e-12
    assert final.outer_band_fraction < 1.0e-12
    assert 0.0 <= final.low_frequency_fraction <= 1.0
    assert final.neglected_source_fraction < 1.0e-8
    assert len(final.shell_fractions) >= 2


# --------------------------------------------------------------------------- #
# Convergence                                                                  #
# --------------------------------------------------------------------------- #


def test_joint_spatial_refinement_is_second_order_in_the_curl_defect():
    curls = []
    divergences = []
    for nr in (33, 65, 129):
        _g, _e, summary = _run(nr=nr, nz=2 * nr - 1, steps=40)
        curls.append(summary.streaming_max_curl_defect)
        divergences.append(summary.streaming_max_relative_divergence)
    assert min(_order(curls[i], curls[i + 1]) for i in range(2)) > 1.8
    assert min(_order(divergences[i], divergences[i + 1]) for i in range(2)) > 1.7


def test_radial_and_axial_refinement_separate_the_two_error_sources():
    """Refining ``dr`` alone fixes the divergence; refining ``dz`` alone fixes the curl."""
    radial = [
        _run(nr=nr, nz=129, steps=40)[2] for nr in (33, 65, 129)
    ]
    axial = [_run(nr=65, nz=nz, steps=40)[2] for nz in (65, 129, 257)]
    radial_divergence = [s.streaming_max_relative_divergence for s in radial]
    assert radial_divergence[0] > radial_divergence[1] > radial_divergence[2]
    axial_curl = [s.streaming_max_curl_defect for s in axial]
    assert axial_curl[0] > axial_curl[1] > axial_curl[2]
    # The complementary quantity saturates, which is the point of separating them.
    axial_divergence = [s.streaming_max_relative_divergence for s in axial]
    assert abs(axial_divergence[-1] - axial_divergence[0]) < 0.05 * axial_divergence[0]


def test_time_refinement_converges():
    values = []
    for factor in (1, 2, 4):
        _g, _e, summary = _run(steps=40 * factor, dt=2.5e-3 / factor)
        values.append(summary.final.max_omega1)
    assert abs(values[0] - values[-1]) > abs(values[1] - values[-1])


@pytest.mark.parametrize("integrator", list(INTEGRATORS))
def test_every_integrator_reaches_the_same_answer(integrator):
    _g, _e, reference = _run(steps=40, integrator="rk4")
    _g2, _e2, summary = _run(steps=40, integrator=integrator)
    assert not summary.diverged
    assert summary.final.max_omega1 == pytest.approx(
        reference.final.max_omega1, rel=1.0e-3
    )


def test_boundary_order_barely_changes_the_small_amplitude_answer():
    """The field never reaches the boundary, so the trace order must not matter."""
    _g, _e, monopole = _run(steps=40, boundary_mode="monopole")
    _g2, _e2, zero = _run(steps=40, boundary_mode="zero")
    assert monopole.final.max_omega1 == pytest.approx(
        zero.final.max_omega1, rel=1.0e-6
    )


@pytest.mark.parametrize("extent", [3.0, 4.0, 6.0])
def test_domain_enlargement_leaves_the_answer_unchanged(extent):
    spacing = 4.0 / 64.0
    _g, _e, summary = _run(
        nr=int(round(extent / spacing)) + 1, nz=129, r_max=extent, steps=40
    )
    _g2, _e2, reference = _run(steps=40)
    assert summary.final.max_omega1 == pytest.approx(
        reference.final.max_omega1, rel=1.0e-6
    )


# --------------------------------------------------------------------------- #
# Independent Cartesian audit                                                  #
# --------------------------------------------------------------------------- #


def test_cartesian_audit_of_a_numerical_state_is_consistent():
    grid, evolution, _summary = _run(nr=129, nz=257, steps=20)
    u1, omega1 = initial_state(grid, BUMP)
    for _ in range(20):
        u1, omega1 = evolution.step(u1, omega1, 2.5e-3, "rk4")
    report = cartesian_audit_of_state(
        grid, u1, omega1, evolution.solver.solve(omega1),
        half_width=2.0, points=41,
    )
    assert report["curl_relative"] < 0.1
    assert report["divergence_relative"] < 0.1
    with pytest.raises(ValueError):
        cartesian_audit_of_state(
            grid, u1, omega1, evolution.solver.solve(omega1),
            half_width=2.0, points=4,
        )


# --------------------------------------------------------------------------- #
# The dyadic cascade toy model                                                 #
# --------------------------------------------------------------------------- #


def _seed(shells: int) -> np.ndarray:
    out = np.zeros(shells, dtype=np.float64)
    out[0] = 1.0e-3
    return out


def test_cascade_nonlinear_term_conserves_energy_exactly():
    model = CascadeModel(shells=10, viscosity=0.0, forcing_amplitude=0.0)
    rng = np.random.default_rng(7)
    for _ in range(20):
        state = rng.standard_normal(10)
        assert model.transfer_defect(state) < 1.0e-12


def test_cascade_rejects_bad_configurations():
    with pytest.raises(ValueError):
        CascadeModel(shells=4, viscosity=1.0e-3)
    with pytest.raises(ValueError):
        CascadeModel(shells=10, viscosity=-1.0)
    with pytest.raises(ValueError):
        CascadeModel(shells=10, viscosity=1.0e-3, high_shell_start=0)
    with pytest.raises(ValueError):
        # Forcing must not reach the band whose amplitude is being measured.
        CascadeModel(
            shells=10, viscosity=1.0e-3, forced_shells=5, high_shell_start=4
        )
    model = CascadeModel(shells=10, viscosity=1.0e-3)
    with pytest.raises(ValueError):
        model.run(_seed(10), steps=0, time_step=1.0e-3)
    with pytest.raises(ValueError):
        model.run(np.zeros(3), steps=10, time_step=1.0e-3)


def test_low_mode_forcing_drives_the_high_shells_with_zero_direct_injection():
    """The central Track-F question, answered in the smallest honest model.

    The force is supported strictly in shells 0-1 and the measured band starts
    at shell 4, so ``direct_high_shell_injection`` is exactly zero.  If the
    high-shell amplitude nevertheless grows by many orders of magnitude, then
    "the force has small high-frequency coefficients" does **not** imply "the
    force cannot feed the cascade".
    """
    unforced = CascadeModel(
        shells=10, viscosity=2.0e-3, forced_shells=0, forcing_amplitude=0.0,
        high_shell_start=4,
    ).run(_seed(10), steps=100000, time_step=2.0e-4)
    forced = CascadeModel(
        shells=10, viscosity=2.0e-3, forced_shells=2, forcing_amplitude=3.0e-3,
        high_shell_start=4,
    ).run(_seed(10), steps=100000, time_step=2.0e-4)
    assert forced.direct_high_shell_injection == 0.0
    assert unforced.max_high_shell_amplitude < 1.0e-20
    assert forced.max_high_shell_amplitude > 1.0e-4
    assert forced.bandwidth_grew
    assert forced.max_active_shell > unforced.max_active_shell


def test_viscosity_competes_with_the_indirect_route():
    def run(nu: float):
        return CascadeModel(
            shells=10, viscosity=nu, forced_shells=2, forcing_amplitude=3.0e-3,
            high_shell_start=4,
        ).run(_seed(10), steps=100000, time_step=2.0e-4)

    weak = run(2.0e-3)
    strong = run(2.0e-2)
    assert strong.max_high_shell_amplitude < weak.max_high_shell_amplitude
    assert strong.max_bandwidth < weak.max_bandwidth


def test_forcing_role_decomposition_closes_exactly_one_channel():
    roles = decompose_forcing_roles()
    closed = [role for role in roles if not role.available_to_a_clay_force]
    assert len(closed) == 1
    assert closed[0].name == "direct_high_frequency_injection"
    assert len(roles) >= 5
    assert all(role.reason for role in roles)
    assert roles[0].as_dict()["available_to_a_clay_force"] is False
