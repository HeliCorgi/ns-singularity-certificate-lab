r"""P0-A / P0-B / P0-C / P1-C tests for the instrumented integrator.

Covered here:

* ``take_step`` with the cross-checking integrators (``ssprk3``, ``rk4``)
  next to ``heun``: fixed points, agreement at small ``dt``, measured time
  convergence (recorded, with the honest caveat that all methods are
  *projected* Runge--Kutta schemes through :func:`constrain_state`).
* The per-accepted-step ``step_stream`` / ``gate_summary`` machinery: every
  acceptance-critical quantity is recorded on every step, so a violation that
  occurs *between* two strided ``history`` rows is still caught.  The
  synthetic-trajectory test mandated by P0-C lives here
  (:func:`test_gate_catches_violation_between_history_rows`).
* Stage-CFL step rejection (P0-B) and the pre/predictor/post CFL records.
* The ``viscosity_sign`` fault: invisible in the Hou operating regime, it must
  be rejected by a diffusion-dominated control problem and by the
  energy-balance defect (P1-C).

Nothing here is evidence about the continuum equations; these are tests of
the discrete instrumentation.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    BINDING_CONSTRAINT_CODES,
    GATE_SUMMARY_FIELDS,
    STEP_STREAM_FIELDS,
    VALID_INTEGRATORS,
    adaptive_time_step_detail,
    constrain_state,
    directional_cfl,
    enstrophy,
    heun_step,
    hou_initial_swirl,
    integrate,
    kinetic_energy,
    swirl_energy_rates,
    take_step,
    wall_enstrophy_flux,
)


def _grid(nr: int, nz: int) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=1.0,
        z_min=0.0,
        z_max=1.0,
        periodic_z=True,
    )


def _small_state(grid: AxisymmetricGrid, amplitude: float = 100.0):
    return constrain_state(
        grid,
        hou_initial_swirl(grid, amplitude=amplitude),
        np.zeros(grid.shape),
    )


def _smooth_swirl(grid: AxisymmetricGrid, amplitude: float) -> np.ndarray:
    """A well-resolved control datum: ``A (1-r^2)^2 sin(2 pi z)``.

    Unlike the E-29 profile, whose radial FWHM is ~3 cells at nr = 17, this
    field is fully resolved on the coarse control grids, so its energy-balance
    defect is a genuine refinement quantity instead of being dominated by an
    under-resolved front.
    """

    radius, axial = grid.mesh()
    return np.asarray(
        amplitude * (1.0 - radius**2) ** 2 * np.sin(2.0 * np.pi * axial),
        dtype=np.float64,
    )


SCHEDULE = [[0.0, 5.0e-3]]


# --------------------------------------------------------------- take_step


def test_take_step_rejects_unknown_method() -> None:
    grid = _grid(9, 16)
    state = _small_state(grid)
    with pytest.raises(ValueError, match="unknown time integrator"):
        take_step(
            grid, state, method="ab2", dt=1.0e-5, time=0.0, viscosity=1.0e-3
        )
    with pytest.raises(ValueError, match="unknown time integrator"):
        integrate(
            grid,
            u1=state.u1,
            omega1=state.omega1,
            t_final=1.0e-4,
            viscosity_schedule=SCHEDULE,
            time_integrator="ab2",
        )


@pytest.mark.parametrize("method", VALID_INTEGRATORS)
def test_zero_field_is_a_fixed_point_of_every_integrator(method) -> None:
    grid = _grid(17, 32)
    state = constrain_state(grid, np.zeros(grid.shape), np.zeros(grid.shape))
    advanced, stages = take_step(
        grid, state, method=method, dt=1.0e-4, time=0.0, viscosity=5.0e-3
    )
    assert np.max(np.abs(advanced.u1)) == 0.0
    assert np.max(np.abs(advanced.omega1)) == 0.0
    expected_stage_counts = {"heun": 1, "ssprk3": 2, "rk4": 3}
    assert len(stages) == expected_stage_counts[method]


def test_heun_step_still_equals_the_take_step_heun_path() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid)
    via_wrapper = heun_step(
        grid, state, dt=2.0e-5, time=0.0, viscosity=5.0e-3
    )
    via_take, _ = take_step(
        grid, state, method="heun", dt=2.0e-5, time=0.0, viscosity=5.0e-3
    )
    assert np.array_equal(via_wrapper.u1, via_take.u1)
    assert np.array_equal(via_wrapper.omega1, via_take.omega1)


def test_integrators_agree_at_small_dt() -> None:
    """The three schemes share the spatial discretization, so their small-dt
    solutions must agree to the time-discretization error, not O(1)."""

    grid = _grid(17, 32)
    state = _small_state(grid)
    finals = {}
    for method in VALID_INTEGRATORS:
        current = state
        for step in range(8):
            current, _ = take_step(
                grid,
                current,
                method=method,
                dt=2.0e-6,
                time=step * 2.0e-6,
                viscosity=5.0e-3,
            )
        finals[method] = current
    scale = float(np.max(np.abs(finals["heun"].u1)))
    for method in ("ssprk3", "rk4"):
        difference = float(
            np.max(np.abs(finals[method].u1 - finals["heun"].u1))
        )
        assert difference / scale < 1.0e-9
    # ... but they are not bitwise identical (different schemes).
    assert not np.array_equal(finals["rk4"].u1, finals["heun"].u1)


def _run_fixed(grid, state, method: str, t_final: float, n: int):
    dt = t_final / n
    current = state
    time = 0.0
    for _ in range(n):
        current, _ = take_step(
            grid, current, method=method, dt=dt, time=time, viscosity=5.0e-3
        )
        time += dt
    return current


def test_measured_time_convergence_of_projected_schemes() -> None:
    """Richardson order against a fine-dt reference of the SAME method.

    All methods are projected through :func:`constrain_state` every stage.
    Measured on this nonlinear problem (amplitude 3000, ``17 x 32``,
    ``t = 3.2e-4``): heun 1.97/2.00, ssprk3 3.00/3.00, rk4 3.95/3.98 -- the
    projection does not degrade the observed classical orders here.  The
    assertion bounds are set below those measurements.
    """

    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=3000.0)
    t_final = 3.2e-4
    lower_bound = {"heun": 1.8, "ssprk3": 2.7, "rk4": 3.5}
    upper_bound = {"heun": 2.3, "ssprk3": 3.3, "rk4": 4.3}
    coarsest_error = {}
    for method in VALID_INTEGRATORS:
        reference = _run_fixed(grid, state, method, t_final, 256)
        errors = [
            float(
                np.max(
                    np.abs(
                        _run_fixed(grid, state, method, t_final, n).u1
                        - reference.u1
                    )
                )
            )
            for n in (8, 16, 32)
        ]
        assert errors[0] > errors[1] > errors[2] > 0.0
        for order in (
            math.log2(errors[0] / errors[1]),
            math.log2(errors[1] / errors[2]),
        ):
            assert lower_bound[method] < order < upper_bound[method]
        coarsest_error[method] = errors[0]
    # Higher-order methods are genuinely more accurate at the same dt.
    assert coarsest_error["rk4"] < coarsest_error["ssprk3"]
    assert coarsest_error["ssprk3"] < coarsest_error["heun"]


# ------------------------------------------------- step stream and gates


def test_step_stream_and_gate_summary_have_the_documented_fields() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid)
    result = integrate(
        grid,
        u1=state.u1,
        omega1=state.omega1,
        t_final=4.0e-4,
        viscosity_schedule=SCHEDULE,
        max_time_step=5.0e-5,
    )
    assert result.integrator == "heun"
    assert set(result.step_stream) == set(STEP_STREAM_FIELDS)
    lengths = {values.size for values in result.step_stream.values()}
    assert lengths == {result.step_count}
    assert set(result.gate_summary) == set(GATE_SUMMARY_FIELDS)
    assert result.gate_summary["steps_accepted"] == float(result.step_count)
    assert result.gate_summary["steps_rejected"] == 0.0
    assert result.rejected_step_count == 0
    # Binding-constraint counts add up to the accepted steps.
    total = sum(
        result.gate_summary[f"binding_count_{name}"]
        for name in BINDING_CONSTRAINT_CODES
    )
    assert total == float(result.step_count)
    # The Poisson algebraic residual of the exactly-solved system is
    # round-off, far below any discretization scale.
    assert (
        result.gate_summary["max_poisson_algebraic_residual_relative"]
        < 1.0e-12
    )


def test_streaming_maxima_match_a_stride_one_history() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid)
    kwargs = dict(
        u1=state.u1,
        omega1=state.omega1,
        t_final=4.0e-4,
        viscosity_schedule=SCHEDULE,
        max_time_step=2.0e-5,
    )
    dense = integrate(grid, diagnostic_stride=1, **kwargs)
    sparse = integrate(grid, diagnostic_stride=1000, **kwargs)
    # Identical trajectories: the stream does not depend on the stride.
    for name in STEP_STREAM_FIELDS:
        assert np.array_equal(
            dense.step_stream[name], sparse.step_stream[name]
        )
    assert dense.gate_summary == sparse.gate_summary
    # The dense history's post-state CFL maximum is exactly the streaming
    # maximum; the sparse history has almost no rows to read it from.
    dense_history_max = max(
        record["advective_cfl"] for record in dense.history[1:]
    )
    assert dense.gate_summary["max_cfl_post_state"] == pytest.approx(
        dense_history_max, rel=1.0e-12
    )
    assert len(sparse.history) < len(dense.history)


def test_gate_catches_violation_between_history_rows() -> None:
    """P0-C mandated synthetic trajectory: a forcing pulse that lives entirely
    between two strided history rows injects energy and is gone again by the
    next recorded row.  The strided history shows monotone decay; only the
    per-step stream sees the violation."""

    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=500.0)
    dt = 1.0e-5
    n_steps = 100
    stride = 50
    pulse_shape = hou_initial_swirl(grid, amplitude=1.0)
    # Inject energy over (52.5, 57.5) steps, extract it again over
    # (57.5, 62.5): the leftover is second-order small, so the strided rows
    # at steps 0/50/100 still decay monotonically while the stream sees a
    # first-order energy increase mid-pulse.
    t_on, t_mid, t_off = 52.5 * dt, 57.5 * dt, 62.5 * dt

    def pulse(t: float):
        if t_on <= t <= t_mid:
            window = math.sin(math.pi * (t - t_on) / (t_mid - t_on)) ** 2
            return 2.0e5 * window * pulse_shape
        if t_mid < t <= t_off:
            window = math.sin(math.pi * (t - t_mid) / (t_off - t_mid)) ** 2
            return -2.0e5 * window * pulse_shape
        return np.zeros(grid.shape)

    result = integrate(
        grid,
        u1=state.u1,
        omega1=state.omega1,
        t_final=n_steps * dt,
        viscosity_schedule=[[0.0, 2.0e-2]],
        fixed_time_step=dt,
        diagnostic_stride=stride,
        forcing_u1=pulse,
    )
    energies = [record["energy"] for record in result.history]
    # The recorded rows (steps 0, 50, 100) show monotone decay ...
    assert all(later < earlier for earlier, later in zip(energies, energies[1:]))
    # ... but the stream caught the pulse: some accepted step increased the
    # energy by far more than round-off.
    increase = result.gate_summary["max_energy_step_increase"]
    roundoff = 8.0 * np.finfo(np.float64).eps * energies[0]
    assert increase > 1.0e3 * roundoff
    # And the violation indeed happened between the recorded rows.
    steps = result.step_stream["step"]
    increases = result.step_stream["energy_step_increase"]
    worst_step = float(steps[int(np.argmax(increases))])
    assert 50.0 < worst_step < 100.0


def test_stage_cfl_limit_rejects_and_then_satisfies() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=2000.0)
    kwargs = dict(
        u1=state.u1,
        omega1=state.omega1,
        t_final=2.0e-3,
        viscosity_schedule=SCHEDULE,
        cfl_coefficient=0.4,
    )
    free = integrate(grid, **kwargs)
    limit = 0.9 * free.gate_summary["max_cfl_predictor_stage"]
    limited = integrate(grid, stage_cfl_limit=limit, **kwargs)
    assert limited.rejected_step_count > 0
    assert limited.gate_summary["steps_rejected"] == float(
        limited.rejected_step_count
    )
    assert limited.gate_summary["max_cfl_predictor_stage"] <= limit + 1.0e-12
    assert free.gate_summary["max_cfl_predictor_stage"] > limit


def test_stage_cfl_limit_requires_adaptive_mode() -> None:
    grid = _grid(9, 16)
    state = _small_state(grid)
    with pytest.raises(ValueError, match="adaptive"):
        integrate(
            grid,
            u1=state.u1,
            omega1=state.omega1,
            t_final=1.0e-4,
            viscosity_schedule=SCHEDULE,
            fixed_time_step=1.0e-5,
            stage_cfl_limit=0.1,
        )


def test_pre_predictor_post_cfl_are_distinct_measurements() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=2000.0)
    result = integrate(
        grid,
        u1=state.u1,
        omega1=state.omega1,
        t_final=1.0e-3,
        viscosity_schedule=SCHEDULE,
    )
    pre = result.step_stream["cfl_pre_state"]
    stage = result.step_stream["cfl_predictor_stage"]
    post = result.step_stream["cfl_post_state"]
    assert pre.size == stage.size == post.size == result.step_count
    # On a growing solution the predictor stage generally exceeds the
    # pre-state measurement somewhere; all three must be genuinely different
    # arrays, not copies.
    assert not np.array_equal(pre, stage)
    assert not np.array_equal(pre, post)
    assert float(np.max(stage)) >= float(np.max(pre)) - 1.0e-15


def test_adaptive_time_step_detail_names_the_binding_constraint() -> None:
    grid = _grid(17, 32)
    still = constrain_state(grid, np.zeros(grid.shape), np.zeros(grid.shape))
    _, binding = adaptive_time_step_detail(
        grid, still, viscosity=5.0e-3, cfl_coefficient=0.1
    )
    assert binding == "viscous"
    _, binding = adaptive_time_step_detail(
        grid,
        still,
        viscosity=5.0e-3,
        cfl_coefficient=0.1,
        max_time_step=1.0e-9,
    )
    assert binding == "max_time_step"
    # A synthetic strong meridional flow: a large omega1 drives psi1 and
    # therefore u^r, u^z through the elliptic solve, far past the viscous
    # bound of this grid.
    radius, axial = grid.mesh()
    omega1 = 1.0e6 * (1.0 - radius**2) * np.sin(2.0 * np.pi * axial)
    moving = constrain_state(grid, np.zeros(grid.shape), omega1)
    assert float(np.max(np.abs(moving.u_z))) > 1.0
    step, binding = adaptive_time_step_detail(
        grid, moving, viscosity=5.0e-3, cfl_coefficient=0.1
    )
    assert binding in ("advective_r", "advective_z")
    assert directional_cfl(grid, moving, step) == pytest.approx(0.1, rel=1e-9)


# ---------------------------------------------- energy balance and P1-C


def test_energy_balance_defect_converges_on_a_control_problem() -> None:
    """The discrete defect against ``dE/dt = -nu int |omega|^2 dV - nu (wall
    work)`` must fall under simultaneous space-time refinement (P1-C item 3).

    Measured with the smooth control datum: 8.93 -> 2.40 -> 0.61 absolute
    (relative 5.9e-2 -> 1.6e-2 -> 4.2e-3) over the three levels; the test
    runs the first two and asserts the reduction.
    """

    defects = []
    relatives = []
    for nr, nz, dt in ((17, 32, 4.0e-5), (33, 64, 1.0e-5)):
        grid = _grid(nr, nz)
        state = constrain_state(
            grid, _smooth_swirl(grid, 100.0), np.zeros(grid.shape)
        )
        result = integrate(
            grid,
            u1=state.u1,
            omega1=state.omega1,
            t_final=4.0e-4,
            viscosity_schedule=SCHEDULE,
            fixed_time_step=dt,
        )
        defects.append(result.gate_summary["max_abs_energy_balance_defect"])
        relatives.append(
            result.gate_summary["max_energy_balance_defect_relative"]
        )
    assert defects[1] < 0.5 * defects[0]
    assert relatives[1] < 0.5 * relatives[0]
    assert relatives[0] < 0.1


def test_wall_flux_term_is_the_difference_between_the_two_defects() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=500.0)
    result = integrate(
        grid,
        u1=state.u1,
        omega1=state.omega1,
        t_final=2.0e-4,
        viscosity_schedule=SCHEDULE,
        fixed_time_step=1.0e-5,
    )
    with_wall = result.step_stream["energy_balance_defect"]
    without = result.step_stream["energy_balance_defect_no_wall"]
    flux = result.step_stream["wall_enstrophy_flux"]
    nu = result.step_stream["viscosity"]
    # defect - defect_no_wall = nu * (trapezoid wall flux); check the last
    # step exactly from the stream (trapezoid of consecutive fluxes).
    reconstructed = with_wall[1:] - without[1:]
    trapezoid = 0.5 * nu[1:] * (flux[1:] + flux[:-1])
    # The reconstruction subtracts two nearly equal O(|energy rate|) sums, so
    # it only holds up to the round-off of those larger terms.
    atol = 64.0 * np.finfo(np.float64).eps * float(np.max(np.abs(without)))
    assert np.allclose(reconstructed, trapezoid, rtol=1.0e-4, atol=atol)


def test_viscosity_sign_fault_is_rejected_by_the_diffusion_control() -> None:
    """P1-C: at nu = 5e-4 with the Hou datum the viscous term is ~3e-4 of the
    term sum and a sign flip is invisible in that budget.  The
    diffusion-dominated control problem must reject it loudly instead."""

    grid = _grid(17, 32)
    state = constrain_state(
        grid, _smooth_swirl(grid, 100.0), np.zeros(grid.shape)
    )
    kwargs = dict(
        u1=state.u1,
        omega1=state.omega1,
        t_final=4.0e-4,
        viscosity_schedule=[[0.0, 2.0e-2]],
        fixed_time_step=5.0e-6,
    )
    clean = integrate(grid, **kwargs)
    corrupted = integrate(grid, fault="viscosity_sign", **kwargs)
    clean_defect = clean.gate_summary["max_energy_balance_defect_relative"]
    corrupted_defect = corrupted.gate_summary[
        "max_energy_balance_defect_relative"
    ]
    # The flipped sign makes dE/dt = +nu * dissipation, so the relative
    # defect saturates near 2 while the clean value is a refinement quantity.
    # Measured: clean 2.21e-2, corrupted 2.000 (ratio 90x).
    assert corrupted_defect > 10.0 * clean_defect
    assert corrupted_defect > 1.0
    # Energy grows monotonically under the flip and decays cleanly without.
    assert clean.gate_summary["final_energy_minus_initial"] < 0.0
    assert corrupted.gate_summary["final_energy_minus_initial"] > 0.0
    assert corrupted.gate_summary["max_energy_step_increase"] > 0.0


def test_swirl_power_decomposition_viscous_term_is_negative_semidefinite_at_rest() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=500.0)
    rates = swirl_energy_rates(grid, state, viscosity=5.0e-3)
    # The initial datum has u^r = u^z = 0 and psi1 = 0: no advection or
    # stretching power, strictly dissipative viscous power.
    assert rates["advection"] == 0.0
    assert rates["stretching"] == 0.0
    assert rates["viscous"] < 0.0


def test_wall_enstrophy_flux_matches_a_direct_quadrature() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=500.0)
    state = heun_step(grid, state, dt=1.0e-5, time=0.0, viscosity=5.0e-3)
    expected = (
        2.0
        * math.pi
        * grid.r[-1] ** 2
        * grid.dz
        * float(np.sum(state.u_z[-1] * state.omega1[-1]))
    )
    assert wall_enstrophy_flux(grid, state) == pytest.approx(expected)
    # Sanity: the initial rest state has zero flux.
    rest = _small_state(grid, amplitude=500.0)
    assert wall_enstrophy_flux(grid, rest) == 0.0


def test_energy_and_enstrophy_helpers_are_consistent_with_the_stream() -> None:
    grid = _grid(17, 32)
    state = _small_state(grid, amplitude=200.0)
    result = integrate(
        grid,
        u1=state.u1,
        omega1=state.omega1,
        t_final=1.0e-4,
        viscosity_schedule=SCHEDULE,
        fixed_time_step=1.0e-5,
    )
    assert result.step_stream["energy"][-1] == pytest.approx(
        kinetic_energy(grid, result.state), rel=1.0e-12
    )
    assert result.step_stream["enstrophy"][-1] == pytest.approx(
        enstrophy(grid, result.state), rel=1.0e-12
    )
