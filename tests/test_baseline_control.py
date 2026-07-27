from __future__ import annotations

import numpy as np

from experiments.run_baseline import (
    _blowup_guard,
    _observed_orders,
    _relative_l2,
    analytic_u1,
    evolve,
)


def test_gaussian_closed_form_satisfies_radial_heat_equation() -> None:
    viscosity = 0.07
    sigma = 1.2
    time = 0.3
    r = np.linspace(0.0, 3.0, 101)
    values = analytic_u1(
        r,
        time,
        viscosity=viscosity,
        amplitude=0.8,
        sigma=sigma,
    )
    width2 = sigma * sigma + 4.0 * viscosity * time
    exact_time_derivative = values * (
        -8.0 * viscosity / width2
        + 4.0 * viscosity * r * r / (width2 * width2)
    )
    exact_radial_operator = values * (
        -8.0 / width2 + 4.0 * r * r / (width2 * width2)
    )
    assert np.allclose(
        exact_time_derivative,
        viscosity * exact_radial_operator,
        rtol=2e-15,
        atol=2e-15,
    )


def test_independent_crank_nicolson_control_converges_near_second_order() -> None:
    resolutions = [17, 33, 65]
    errors: list[float] = []
    for resolution in resolutions:
        result = evolve(
            nr=resolution,
            r_max=5.0,
            final_time=0.2,
            viscosity=0.05,
            amplitude=1.0,
            sigma=1.0,
            dt_over_dr=0.1,
            sample_count=9,
            exact_outer_boundary=True,
        )
        exact = analytic_u1(
            result.r,
            0.2,
            viscosity=0.05,
            amplitude=1.0,
            sigma=1.0,
        )
        errors.append(
            _relative_l2(result.r, result.profiles[-1] - exact, exact)
        )
    orders = _observed_orders(errors, resolutions)
    assert all(fine < coarse for coarse, fine in zip(errors, errors[1:]))
    assert min(orders) > 1.8


def test_decaying_control_does_not_trigger_blowup_fit() -> None:
    result = evolve(
        nr=33,
        r_max=5.0,
        final_time=0.2,
        viscosity=0.05,
        amplitude=1.0,
        sigma=1.0,
        dt_over_dr=0.1,
        sample_count=9,
        exact_outer_boundary=True,
    )
    guard = _blowup_guard(result)
    assert guard["classification"] == "bounded_or_decaying_control"
    assert guard["sustained_tail_growth"] is False
    assert guard["reciprocal_blowup_fit_attempted"] is False
