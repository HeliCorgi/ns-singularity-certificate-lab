"""Acceptance tests for the production nonlinear finite-cylinder integrator.

These implement items 1-9 of ``docs/nonlinear_solver_design.md`` section 7.
Item 10 (independent solver-B residuals on stored snapshots) belongs to the
experiment test ``tests/test_hou_early_time_experiment.py``.

Manufactured derivation (all of it done by hand and transcribed below)
======================================================================

Let ``q = 2*pi`` and work on ``r in [0,1]``, ``z in [0,1)`` with period one.

    psi1*(r,z,t) = A(t) P(r) cos(q z),   A(t) = 0.15 (1 + 0.5 sin 3t)
    u1*(r,z,t)   = B(t) Q(r) sin(q z),   B(t) = 0.40 (1 + 0.3 cos 2t)

with the axis-even wall-vanishing radial profiles

    P(r) = (1-r^2)^2 (1 + 0.2 r^2) = 1 - 1.8 r^2 + 0.6 r^4 + 0.2 r^6,
    Q(r) = (1-r^2)^2 (1 + 0.1 r^2) = 1 - 1.9 r^2 + 0.8 r^4 + 0.1 r^6.

For a monomial ``r^(2k)`` the radial part of ``L5`` is

    (d_rr + (3/r) d_r) r^(2k) = [2k(2k-1) + 6k] r^(2k-2) = 4k(k+1) r^(2k-2),

so ``k = 0,1,2,3`` give the factors ``0, 8, 24, 48``.  Hence

    L5r P = -1.8*8 + 0.6*24 r^2 + 0.2*48 r^4 = -14.4 + 14.4 r^2 + 9.6 r^4,
    L5r Q = -1.9*8 + 0.8*24 r^2 + 0.1*48 r^4 = -15.2 + 19.2 r^2 + 4.8 r^4.

E-13 gives ``omega1* = -L5 psi1* = A(t) S(r) cos(q z)`` with

    S(r) = q^2 P(r) - L5r P(r)
         = (q^2 + 14.4) + (-1.8 q^2 - 14.4) r^2 + (0.6 q^2 - 9.6) r^4
           + 0.2 q^2 r^6.

Wall consistency check (E-27/E-31).  With ``P = (1-r^2)^2 (1+0.2 r^2)``,

    d_r[(1-r^2)^2] = -4r(1-r^2)          -> 0 at r=1,
    d_rr[(1-r^2)^2] = -4(1-r^2) + 8r^2   -> 8 at r=1,
    (1-r^2)^2                            -> 0 at r=1,

so ``P''(1) = 8 * (1 + 0.2) = 9.6`` and the cross terms vanish.  Directly from
the polynomial, ``P'' = -3.6 + 7.2 r^2 + 6 r^4`` gives ``P''(1) = 9.6`` as
well.  Meanwhile ``S(1) = q^2 (1 - 1.8 + 0.6 + 0.2) - 9.6 = -9.6``, so

    omega1*(1,z,t) = -9.6 A(t) cos(q z) = -psi1*_rr(1,z,t)

exactly, and ``P(1) = P'(1) = 0`` gives ``psi1*(1,z) = psi1*_r(1,z) = 0`` and
``u1*(1,z) = 0``.  The manufactured trajectory therefore satisfies every wall
condition of E-27 exactly.

Velocities (E-14):

    u^r* = -r psi1*_z = q r A P sin(q z),
    u^z* = 2 psi1* + r psi1*_r = A W(r) cos(q z),
    W(r) = 2 P + r P' = 2 - 7.2 r^2 + 3.6 r^4 + 1.6 r^6.

Forcings are the exact continuum residuals

    F_u1 = du1*/dt + u^r* u1*_r + u^z* u1*_z - 2 u1* psi1*_z - nu L5 u1*
         = B' Q sin
           + q A B [ r P Q' sin^2 + W Q cos^2 ]
           + 2 q A B P Q sin^2
           - nu B (L5r Q - q^2 Q) sin,

    F_om1 = domega1*/dt + u^r* omega1*_r + u^z* omega1*_z - d_z(u1*^2)
            - nu L5 omega1*
          = A' S cos
            + q A^2 (r P S' - W S) sin cos
            - q B^2 Q^2 sin(2 q z)
            - nu A (L5r S - q^2 S) cos,

using ``2 u1* psi1*_z = -2 q A B P Q sin^2`` and
``d_z(u1*^2) = q B^2 Q^2 sin(2 q z)``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np
import pytest

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    VALID_FAULTS,
    adaptive_time_step,
    axis_parity_defect,
    cartesian_vorticity,
    constrain_state,
    hou_initial_swirl,
    integrate,
    normalize_viscosity_schedule,
    odd_symmetry_defect,
    resume,
    save_checkpoint,
    thom_wall_vorticity,
    viscosity_at,
)
from ns_certificate_lab.operators import laplacian_5d_formal


WAVE_NUMBER = 2.0 * math.pi
MANUFACTURED_VISCOSITY = 0.05
MANUFACTURED_FINAL_TIME = 0.05
MANUFACTURED_STEP_COEFFICIENT = 0.2
MANUFACTURED_RESOLUTIONS = ((17, 32), (33, 64), (65, 128))
# Minimum observed orders.  The coarse omega1 pair is pre-asymptotic because
# the sup error sits on the E-31 wall row; the measured sequence for omega1 is
# 1.845, 1.902 and continues to 1.951 on a (129,256) grid that is too slow for
# the test suite.  u1 and psi1 are at 2.00 on every pair.
MINIMUM_OBSERVED_ORDER = 1.80
MINIMUM_FINEST_OBSERVED_ORDER = 1.85
# Preregistered rejection thresholds, set roughly five times above the clean
# sup errors measured on the (33,64) grid used by the fault tests.
FAULT_REJECTION_SUP_ERROR = {"u1": 5.0e-4, "omega1": 5.0e-2, "psi1": 4.0e-4}


def _grid(nr: int, nz: int) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr,
        nz=nz,
        r_max=1.0,
        z_min=0.0,
        z_max=1.0,
        periodic_z=True,
    )


@dataclass(frozen=True)
class Manufactured:
    """Closed-form manufactured trajectory and its exact forcings."""

    grid: AxisymmetricGrid
    viscosity: float

    def _profiles(self) -> dict[str, np.ndarray]:
        radius, axial = self.grid.mesh()
        r2 = radius * radius
        q = WAVE_NUMBER
        s0 = q * q + 14.4
        s1 = -1.8 * q * q - 14.4
        s2 = 0.6 * q * q - 9.6
        s3 = 0.2 * q * q
        return {
            "r": radius,
            "cos": np.cos(q * axial),
            "sin": np.sin(q * axial),
            "sin_double": np.sin(2.0 * q * axial),
            "P": 1.0 - 1.8 * r2 + 0.6 * r2**2 + 0.2 * r2**3,
            "W": 2.0 - 7.2 * r2 + 3.6 * r2**2 + 1.6 * r2**3,
            "S": s0 + s1 * r2 + s2 * r2**2 + s3 * r2**3,
            "Sp": 2.0 * s1 * radius
            + 4.0 * s2 * radius**3
            + 6.0 * s3 * radius**5,
            "L5rS": 8.0 * s1 + 24.0 * s2 * r2 + 48.0 * s3 * r2**2,
            "Q": 1.0 - 1.9 * r2 + 0.8 * r2**2 + 0.1 * r2**3,
            "Qp": -3.8 * radius + 3.2 * radius**3 + 0.6 * radius**5,
            "L5rQ": -15.2 + 19.2 * r2 + 4.8 * r2**2,
        }

    @staticmethod
    def amplitudes(time: float) -> tuple[float, float, float, float]:
        a = 0.15 * (1.0 + 0.5 * math.sin(3.0 * time))
        a_dot = 0.225 * math.cos(3.0 * time)
        b = 0.4 * (1.0 + 0.3 * math.cos(2.0 * time))
        b_dot = -0.24 * math.sin(2.0 * time)
        return a, a_dot, b, b_dot

    def exact(self, time: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        profiles = self._profiles()
        a, _, b, _ = self.amplitudes(time)
        return (
            b * profiles["Q"] * profiles["sin"],
            a * profiles["S"] * profiles["cos"],
            a * profiles["P"] * profiles["cos"],
        )

    def forcing_u1(self, time: float) -> np.ndarray:
        profiles = self._profiles()
        a, _, b, b_dot = self.amplitudes(time)
        q = WAVE_NUMBER
        return (
            b_dot * profiles["Q"] * profiles["sin"]
            + q
            * a
            * b
            * (
                profiles["r"] * profiles["P"] * profiles["Qp"] * profiles["sin"] ** 2
                + profiles["W"] * profiles["Q"] * profiles["cos"] ** 2
            )
            + 2.0 * q * a * b * profiles["P"] * profiles["Q"] * profiles["sin"] ** 2
            - self.viscosity
            * b
            * (profiles["L5rQ"] - q * q * profiles["Q"])
            * profiles["sin"]
        )

    def forcing_omega1(self, time: float) -> np.ndarray:
        profiles = self._profiles()
        a, a_dot, b, _ = self.amplitudes(time)
        q = WAVE_NUMBER
        return (
            a_dot * profiles["S"] * profiles["cos"]
            + q
            * a
            * a
            * (
                profiles["r"] * profiles["P"] * profiles["Sp"]
                - profiles["W"] * profiles["S"]
            )
            * profiles["sin"]
            * profiles["cos"]
            - q * b * b * profiles["Q"] ** 2 * profiles["sin_double"]
            - self.viscosity
            * a
            * (profiles["L5rS"] - q * q * profiles["S"])
            * profiles["cos"]
        )


def manufactured_errors(
    nr: int,
    nz: int,
    *,
    t_final: float = MANUFACTURED_FINAL_TIME,
    viscosity: float = MANUFACTURED_VISCOSITY,
    fault: str | None = None,
) -> dict[str, float]:
    """Run the forced manufactured problem and return sup errors at ``t_final``."""

    grid = _grid(nr, nz)
    fields = Manufactured(grid=grid, viscosity=viscosity)
    u1_initial, omega1_initial, _ = fields.exact(0.0)
    result = integrate(
        grid,
        u1=u1_initial,
        omega1=omega1_initial,
        t_final=t_final,
        viscosity_schedule=[[0.0, viscosity]],
        fixed_time_step=MANUFACTURED_STEP_COEFFICIENT * grid.dr**2,
        forcing_u1=fields.forcing_u1,
        forcing_omega1=fields.forcing_omega1,
        fault=fault,
        diagnostic_stride=10**9,
    )
    exact_u1, exact_omega1, exact_psi1 = fields.exact(t_final)
    return {
        "dr": grid.dr,
        "step_count": float(result.step_count),
        "u1": float(np.max(np.abs(result.state.u1 - exact_u1))),
        "omega1": float(np.max(np.abs(result.state.omega1 - exact_omega1))),
        "psi1": float(np.max(np.abs(result.state.psi1 - exact_psi1))),
    }


def _orders(errors: list[float], spacings: list[float]) -> list[float]:
    return [
        math.log(errors[index] / errors[index + 1])
        / math.log(spacings[index] / spacings[index + 1])
        for index in range(len(errors) - 1)
    ]


@pytest.fixture(scope="module")
def manufactured_convergence() -> list[dict[str, float]]:
    return [manufactured_errors(nr, nz) for nr, nz in MANUFACTURED_RESOLUTIONS]


@pytest.fixture(scope="module")
def unforced_hou_run():
    grid = _grid(33, 64)
    initial = hou_initial_swirl(grid, amplitude=1.0)
    result = integrate(
        grid,
        u1=initial,
        omega1=np.zeros(grid.shape),
        t_final=0.03,
        viscosity_schedule=[[0.0, 5.0e-4]],
        fixed_time_step=1.0e-4,
        diagnostic_stride=1,
    )
    return grid, initial, result


# ---------------------------------------------------------------- item 1


def test_thom_wall_vorticity_converges_at_second_order() -> None:
    """E-31 unit convergence on psi = (1-r^2)^2 phi with phi smooth and even.

    With ``phi = (1 + 0.3 r^2)(1 + 0.5 cos 2 pi z)`` the product rule gives

        d_rr[(1-r^2)^2 phi]
            = phi d_rr[(1-r^2)^2] + 2 d_r[(1-r^2)^2] d_r phi
              + (1-r^2)^2 d_rr phi,

    and at ``r=1`` the last two terms vanish because ``d_r[(1-r^2)^2] = 0`` and
    ``(1-r^2)^2 = 0`` there, while ``d_rr[(1-r^2)^2] = 8``.  Hence
    ``psi_rr(1,z) = 8 phi(1,z) = 10.4 (1 + 0.5 cos 2 pi z)`` and the exact wall
    vorticity is its negative.
    """

    errors: list[float] = []
    spacings: list[float] = []
    for nr in (17, 33, 65):
        grid = _grid(nr, 32)
        radius, axial = grid.mesh()
        phi = (1.0 + 0.3 * radius**2) * (1.0 + 0.5 * np.cos(2.0 * np.pi * axial))
        psi1 = (1.0 - radius**2) ** 2 * phi
        exact = -8.0 * 1.3 * (1.0 + 0.5 * np.cos(2.0 * np.pi * grid.z))
        computed = thom_wall_vorticity(grid, psi1)
        errors.append(float(np.max(np.abs(computed - exact))))
        spacings.append(grid.dr)
    assert errors[0] > errors[1] > errors[2]
    observed = _orders(errors, spacings)
    assert min(observed) >= 1.85, observed
    # A sign flip must not pass as a second-order discretization.
    grid = _grid(33, 32)
    radius, axial = grid.mesh()
    psi1 = (1.0 - radius**2) ** 2 * (1.0 + 0.3 * radius**2)
    flipped = thom_wall_vorticity(grid, psi1, fault="thom_sign")
    assert np.allclose(flipped, -thom_wall_vorticity(grid, psi1))


# ---------------------------------------------------------------- item 2


def test_forced_manufactured_nonlinear_converges_second_order(
    manufactured_convergence: list[dict[str, float]],
) -> None:
    rows = manufactured_convergence
    spacings = [row["dr"] for row in rows]
    for name in ("u1", "omega1", "psi1"):
        errors = [row[name] for row in rows]
        assert errors[0] > errors[1] > errors[2], (name, errors)
        observed = _orders(errors, spacings)
        assert min(observed) >= MINIMUM_OBSERVED_ORDER, (name, observed)
        assert observed[-1] >= MINIMUM_FINEST_OBSERVED_ORDER, (name, observed)
    finest = rows[-1]
    assert finest["u1"] < 1.0e-4
    assert finest["psi1"] < 1.0e-4
    assert finest["omega1"] < 1.0e-2


def test_manufactured_trajectory_satisfies_the_wall_conditions() -> None:
    """The manufactured fields obey E-27 exactly, not merely to O(dr^2)."""

    grid = _grid(65, 128)
    fields = Manufactured(grid=grid, viscosity=MANUFACTURED_VISCOSITY)
    for time in (0.0, 0.037, 0.05):
        exact_u1, exact_omega1, exact_psi1 = fields.exact(time)
        amplitude = fields.amplitudes(time)[0]
        # The expanded polynomial forms cancel to round-off at r=1 rather than
        # to a literal zero; the solver pins the trace exactly.
        assert np.max(np.abs(exact_u1[-1])) < 1.0e-15
        assert np.max(np.abs(exact_psi1[-1])) < 1.0e-15
        wall = -9.6 * amplitude * np.cos(WAVE_NUMBER * grid.z)
        assert np.allclose(exact_omega1[-1], wall, rtol=0.0, atol=1e-12)
        # The discrete E-31 reconstruction reproduces that trace to O(dr^2).
        reconstructed = thom_wall_vorticity(grid, exact_psi1)
        assert np.max(np.abs(reconstructed - wall)) < 5.0e-3 * max(
            1.0, float(np.max(np.abs(wall)))
        )


# ---------------------------------------------------------------- item 3


def test_fixed_grid_time_convergence_is_second_order() -> None:
    """Step doubling on one grid isolates the Heun temporal order."""

    grid = _grid(33, 64)
    fields = Manufactured(grid=grid, viscosity=MANUFACTURED_VISCOSITY)
    u1_initial, omega1_initial, _ = fields.exact(0.0)
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for level in range(4):
        result = integrate(
            grid,
            u1=u1_initial,
            omega1=omega1_initial,
            t_final=0.04,
            viscosity_schedule=[[0.0, MANUFACTURED_VISCOSITY]],
            fixed_time_step=1.0e-3 / 2**level,
            forcing_u1=fields.forcing_u1,
            forcing_omega1=fields.forcing_omega1,
            diagnostic_stride=10**9,
        )
        assert result.step_count == 40 * 2**level
        states.append((result.state.u1.copy(), result.state.omega1.copy()))
    for index, name in ((0, "u1"), (1, "omega1")):
        differences = [
            float(np.max(np.abs(states[k][index] - states[k + 1][index])))
            for k in range(3)
        ]
        assert differences[0] > differences[1] > differences[2] > 0.0
        observed = [
            math.log(differences[k] / differences[k + 1]) / math.log(2.0)
            for k in range(2)
        ]
        assert min(observed) >= 1.9, (name, observed)
        assert max(observed) <= 2.2, (name, observed)


# ---------------------------------------------------------------- item 4


def test_zero_state_is_an_exact_fixed_point() -> None:
    grid = _grid(17, 32)
    zeros = np.zeros(grid.shape)
    result = integrate(
        grid,
        u1=zeros,
        omega1=zeros,
        t_final=0.01,
        viscosity_schedule=[[0.0, 5.0e-4]],
        fixed_time_step=1.0e-3,
    )
    assert result.step_count == 10
    assert np.max(np.abs(result.state.u1)) == 0.0
    assert np.max(np.abs(result.state.omega1)) == 0.0
    assert np.max(np.abs(result.state.psi1)) == 0.0
    assert np.max(np.abs(result.state.u_r)) == 0.0
    assert np.max(np.abs(result.state.u_z)) == 0.0


# ---------------------------------------------------------------- item 5


def _linear_swirl_diffusion(
    grid: AxisymmetricGrid,
    initial: np.ndarray,
    *,
    viscosity: float,
    dt: float,
    steps: int,
) -> np.ndarray:
    """Independent Heun loop on ``du1/dt = nu L5 u1`` with the same wall pin."""

    field = np.array(initial, dtype=np.float64, copy=True)
    field[-1] = 0.0

    def rhs(values: np.ndarray) -> np.ndarray:
        out = viscosity * laplacian_5d_formal(grid, values)
        out[-1] = 0.0
        return out

    for _ in range(steps):
        k1 = rhs(field)
        stage = field + dt * k1
        stage[-1] = 0.0
        k2 = rhs(stage)
        field = field + 0.5 * dt * (k1 + k2)
        field[-1] = 0.0
    return field


def test_small_amplitude_limit_matches_linear_swirl_diffusion() -> None:
    grid = _grid(33, 64)
    radius, axial = grid.mesh()
    shape = (
        (1.0 - radius**2) ** 2
        * (1.0 + 0.3 * radius**2)
        * np.sin(2.0 * np.pi * axial)
    )
    viscosity = 5.0e-3
    t_final = 0.05
    dt = 5.0e-4
    steps = int(round(t_final / dt))
    deviations: dict[float, float] = {}
    for epsilon in (1.0e-2, 1.0e-3):
        initial = epsilon * shape
        result = integrate(
            grid,
            u1=initial,
            omega1=np.zeros(grid.shape),
            t_final=t_final,
            viscosity_schedule=[[0.0, viscosity]],
            fixed_time_step=dt,
            diagnostic_stride=10**9,
        )
        reference = _linear_swirl_diffusion(
            grid,
            initial,
            viscosity=viscosity,
            dt=dt,
            steps=steps,
        )
        deviation = float(np.max(np.abs(result.state.u1 - reference)))
        deviations[epsilon] = deviation
        induced = float(np.max(np.abs(result.state.omega1)))
        # omega1 is generated only by d_z(u1^2) = O(eps^2).
        assert induced <= 10.0 * epsilon**2 * t_final
        assert induced > 0.0
        relative = deviation / float(np.max(np.abs(result.state.u1)))
        assert relative <= epsilon**2
    # The deviation is cubic in the amplitude, i.e. relatively quadratic.
    ratio = deviations[1.0e-2] / deviations[1.0e-3]
    assert 900.0 <= ratio <= 1100.0, ratio


# ---------------------------------------------------------------- item 6


def test_symmetry_defects_stay_at_round_off(unforced_hou_run) -> None:
    grid, initial, result = unforced_hou_run
    amplitude = max(record["max_abs_u1"] for record in result.history)
    odd_defects = [record["odd_symmetry_defect"] for record in result.history]
    assert max(odd_defects) / amplitude < 1.0e-12
    # Measured: 3.33e-16 absolute, 1.23e-15 relative on this grid.
    assert max(odd_defects) < 1.0e-14

    # The axis parity diagnostic is a discretization-order quantity: for an
    # exactly even field the one-sided estimate evaluates to -6 b dr^3 + O(dr^5)
    # with b the quartic Taylor coefficient.  At t=0 the initial data factorizes
    # as (1-r^2)^18 f(z), so the defect is available in closed form.
    axial_factor = float(
        np.max(
            np.abs(
                np.sin(2.0 * np.pi * grid.z)
                / (1.0 + 12.5 * np.sin(np.pi * grid.z) ** 2)
            )
        )
    )
    radial = (1.0 - grid.r**2) ** 18
    predicted = abs(
        (-3.0 * radial[0] + 4.0 * radial[1] - radial[2]) / (2.0 * grid.dr)
    ) * axial_factor
    measured = [record["axis_parity_defect"] for record in result.history]
    assert measured[0] == pytest.approx(predicted, rel=1.0e-12)
    assert predicted == pytest.approx(
        6.0 * 153.0 * grid.dr**3 * axial_factor, rel=0.05
    )
    # The defect must not grow: an odd component of size e would appear as
    # O(e/dr), five orders above this floor.
    assert max(measured) <= measured[0] * (1.0 + 1.0e-12)
    assert max(measured) / amplitude < 0.05

    # The same conclusion read directly off the final state.
    for field in (result.state.u1, result.state.omega1, result.state.psi1):
        scale = max(1.0, float(np.max(np.abs(field))))
        assert odd_symmetry_defect(grid, field) / scale < 1.0e-12
        assert math.isfinite(axis_parity_defect(grid, field))


# ---------------------------------------------------------------- item 7


def test_circulation_and_energy_are_non_increasing(unforced_hou_run) -> None:
    _, _, result = unforced_hou_run
    circulation = [record["circulation_max"] for record in result.history]
    energy = [record["energy"] for record in result.history]
    assert len(circulation) == result.step_count + 1
    for index in range(len(circulation) - 1):
        # E-23: ||r^2 u1||_inf is non-increasing; 1e-12 relative per step.
        assert circulation[index + 1] <= circulation[index] * (1.0 + 1.0e-12)
        assert energy[index + 1] <= energy[index] * (1.0 + 1.0e-12)
    assert circulation[-1] <= circulation[0] * (1.0 + 1.0e-10)
    assert energy[-1] <= energy[0] * (1.0 + 1.0e-10)
    assert circulation[-1] < circulation[0]
    assert energy[-1] < energy[0]
    # E-30 qualitative target: the swirl maximum decreases at very early times.
    swirl = [record["max_abs_u1"] for record in result.history]
    assert swirl[-1] < swirl[0]


# ---------------------------------------------------------------- item 8


@pytest.fixture(scope="module")
def clean_fault_baselines() -> dict[float, dict[str, float]]:
    return {
        0.02: manufactured_errors(33, 64, t_final=0.02),
        0.10: manufactured_errors(33, 64, t_final=0.10),
    }


@pytest.mark.parametrize(
    ("fault", "diagnostic", "t_final"),
    [
        ("thom_sign", "omega1", 0.02),
        ("stretch_sign", "u1", 0.02),
        ("axis_coeff", "u1", 0.02),
        ("poisson_sign", "psi1", 0.02),
        ("wall_leak", "omega1", 0.10),
    ],
)
def test_detects_injected_faults(
    clean_fault_baselines: dict[float, dict[str, float]],
    fault: str,
    diagnostic: str,
    t_final: float,
) -> None:
    assert fault in VALID_FAULTS
    clean = clean_fault_baselines[t_final]
    broken = manufactured_errors(33, 64, t_final=t_final, fault=fault)
    assert clean[diagnostic] < FAULT_REJECTION_SUP_ERROR[diagnostic]
    assert broken[diagnostic] > FAULT_REJECTION_SUP_ERROR[diagnostic]
    assert broken[diagnostic] > 10.0 * clean[diagnostic]


def test_unknown_fault_names_are_rejected() -> None:
    grid = _grid(17, 32)
    with pytest.raises(ValueError, match="unknown injected fault"):
        constrain_state(grid, np.zeros(grid.shape), np.zeros(grid.shape), fault="nope")


# ---------------------------------------------------------------- item 9


def test_restart_reproduces_a_continuous_run(scratch_dir: Path) -> None:
    grid = _grid(17, 32)
    initial = hou_initial_swirl(grid, amplitude=1.0)
    zeros = np.zeros(grid.shape)
    schedule = [[0.0, 5.0e-4]]
    dt = 2.0e-4
    continuous = integrate(
        grid,
        u1=initial,
        omega1=zeros,
        t_final=20 * dt,
        viscosity_schedule=schedule,
        fixed_time_step=dt,
        diagnostic_stride=10**9,
    )
    partial = integrate(
        grid,
        u1=initial,
        omega1=zeros,
        t_final=8 * dt,
        viscosity_schedule=schedule,
        fixed_time_step=dt,
        diagnostic_stride=10**9,
    )
    checkpoint = scratch_dir / "restart.npz"
    save_checkpoint(
        checkpoint,
        grid=grid,
        state=partial.state,
        time=partial.time,
        viscosity=5.0e-4,
        seed=20260728,
        config={"test": "restart_fidelity"},
    )
    restarted = resume(
        checkpoint,
        t_final=20 * dt,
        viscosity_schedule=schedule,
        fixed_time_step=dt,
        diagnostic_stride=10**9,
    )
    assert restarted.step_count == continuous.step_count - partial.step_count
    assert restarted.time == pytest.approx(continuous.time, rel=1e-14)
    for name in ("u1", "omega1", "psi1"):
        reference = getattr(continuous.state, name)
        restored = getattr(restarted.state, name)
        scale = float(np.max(np.abs(reference)))
        assert float(np.max(np.abs(reference - restored))) <= 1.0e-12 * scale


# ------------------------------------------------------- solver mechanics


def test_wall_rows_are_reconstructed_not_evolved() -> None:
    grid = _grid(17, 32)
    initial = hou_initial_swirl(grid, amplitude=1.0)
    # Corrupt the wall traces: the constraint step must repair both of them.
    corrupted_u1 = initial.copy()
    corrupted_u1[-1] = 3.0
    corrupted_omega1 = np.zeros(grid.shape)
    corrupted_omega1[-1] = -7.0
    state = constrain_state(grid, corrupted_u1, corrupted_omega1)
    assert np.max(np.abs(state.u1[-1])) == 0.0
    assert np.allclose(state.omega1[-1], thom_wall_vorticity(grid, state.psi1))
    assert np.max(np.abs(state.psi1[-1])) == 0.0
    # The elliptic solve ignores the wall row of omega1, so the constrained
    # state is a fixed point of the constraint map.
    again = constrain_state(grid, state.u1, state.omega1)
    assert np.array_equal(again.psi1, state.psi1)
    assert np.array_equal(again.omega1, state.omega1)


def test_adaptive_step_and_viscosity_schedule_behave_as_designed() -> None:
    grid = _grid(17, 32)
    state = constrain_state(grid, np.zeros(grid.shape), np.zeros(grid.shape))
    viscosity = 5.0e-4
    expected = 0.1 * min(grid.dr**2, grid.dz**2) / (4.0 * viscosity)
    assert adaptive_time_step(
        grid, state, viscosity=viscosity, cfl_coefficient=0.1
    ) == pytest.approx(expected)
    assert adaptive_time_step(
        grid,
        state,
        viscosity=viscosity,
        cfl_coefficient=0.1,
        max_time_step=1.0e-7,
    ) == pytest.approx(1.0e-7)

    schedule = [[0.0, 5.0e-4], [0.002, 5.0e-3]]
    assert normalize_viscosity_schedule(schedule) == ((0.0, 5.0e-4), (0.002, 5.0e-3))
    assert viscosity_at(schedule, 0.0) == 5.0e-4
    assert viscosity_at(schedule, 0.0019) == 5.0e-4
    assert viscosity_at(schedule, 0.002) == 5.0e-3
    assert viscosity_at(schedule, 1.0) == 5.0e-3
    with pytest.raises(ValueError, match="first viscosity switch time"):
        normalize_viscosity_schedule([[0.001, 5.0e-4]])
    with pytest.raises(ValueError, match="strictly increasing"):
        normalize_viscosity_schedule([[0.0, 5.0e-4], [0.0, 5.0e-3]])

    initial = hou_initial_swirl(grid, amplitude=1.0)
    result = integrate(
        grid,
        u1=initial,
        omega1=np.zeros(grid.shape),
        t_final=0.004,
        viscosity_schedule=schedule,
        cfl_coefficient=0.1,
        max_time_step=2.0e-4,
        diagnostic_stride=1,
    )
    assert result.viscosity_switch_times == (0.002,)
    assert result.completed is True
    recorded = {record["viscosity"] for record in result.history}
    assert recorded == {5.0e-4, 5.0e-3}
    assert all(
        record["time"] <= 0.002 or record["viscosity"] == 5.0e-3
        for record in result.history[1:]
    )


def test_e29_initial_data_matches_the_audited_norms() -> None:
    """E-29b: the derived initial norms are recovered on a fine grid."""

    grid = _grid(129, 256)
    u1 = hou_initial_swirl(grid)
    state = constrain_state(grid, u1, np.zeros(grid.shape))
    assert np.max(np.abs(state.psi1)) == 0.0
    assert np.max(np.abs(state.u_r)) == 0.0
    assert np.max(np.abs(state.u_z)) == 0.0
    magnitude = np.sqrt(
        sum(
            component * component
            for component in cartesian_vorticity(grid, state.u1, state.omega1)
        )
    )
    measured_u1 = float(np.max(np.abs(u1)))
    measured_vorticity = float(np.max(magnitude))
    # Measured on this grid: 3265.5332 (1.39e-4 relative) and 7554.6095
    # (1.98e-3 relative); both are grid maxima of a smooth field.
    assert measured_u1 == pytest.approx(3265.9863, rel=0.01)
    assert measured_vorticity == pytest.approx(7569.62, rel=0.01)
    assert hou_initial_swirl(grid, amplitude_scale=0.5) == pytest.approx(0.5 * u1)


def test_max_steps_stops_an_incomplete_run() -> None:
    grid = _grid(17, 32)
    initial = hou_initial_swirl(grid, amplitude=1.0)
    result = integrate(
        grid,
        u1=initial,
        omega1=np.zeros(grid.shape),
        t_final=1.0,
        viscosity_schedule=[[0.0, 5.0e-4]],
        max_time_step=1.0e-4,
        max_steps=5,
    )
    assert result.completed is False
    assert result.step_count == 5
    assert result.time < 1.0
