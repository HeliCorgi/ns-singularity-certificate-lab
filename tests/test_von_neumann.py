r"""Verification of the frozen-coefficient von Neumann audit (P0-A).

The module under test claims three separable things, and they are checked
separately here:

1.  **Algebra.**  The stability polynomials are the ones named, and Heun's
    ``|G(ia)|^2 = 1 + a^4/4`` identity holds to roundoff.  This is checked
    against a closed form written out in this file, not against the module.
2.  **The symbol.**  ``advection_diffusion_symbol`` is the symbol of the
    centered stencil it claims.  This is checked against array propagators
    that use :func:`numpy.roll` differences and never form a symbol -- one of
    them (the two-dimensional Heun propagator below) is written out in this
    file so that the two-dimensional claim, including the frozen ``3\nu/r``
    term, is validated by an oracle the module does not supply.
3.  **The scan.**  ``scan_amplification`` really maximizes over a wavenumber
    grid that contains ``0`` and ``\pm\pi`` exactly, and its maximum is
    monotone under the refinement the module documents.

What is deliberately *not* claimed anywhere below: that a passing audit makes
the production solver stable, or that a failing audit proves it unstable.  The
model omits the variable-coefficient commutators, the axis row (E-17), the wall
rows (E-27/E-31), the elliptic coupling and the zeroth-order production terms.

Measured constants quoted in the assertions were obtained on this machine with
the pinned environment; the assertions keep margin so that they pin behavior
rather than arithmetic noise.
"""

from __future__ import annotations

import ast
import inspect
import math
from pathlib import Path

import numpy as np
import pytest

import ns_certificate_lab.von_neumann as von_neumann_module
from ns_certificate_lab.von_neumann import (
    advection_diffusion_symbol,
    audit_snapshot,
    reference_mode_amplification,
    reference_propagate,
    scan_amplification,
    stability_polynomial,
    theta_grid,
)

# The Hou-like operating point of test 7, quoted once.
HOU_DR = 1.0 / 192.0
HOU_DZ = 1.0 / 384.0
HOU_DT = 2.76e-7
HOU_MAX_U_R = 2000.0
HOU_MAX_U_Z = 7000.0
HOU_VISCOSITY = 5.0e-4


def _propagate_2d_heun(
    values: np.ndarray,
    *,
    dt: float,
    dr: float,
    dz: float,
    advection_r: float,
    advection_z: float,
    viscosity: float,
    n_steps: int,
) -> np.ndarray:
    r"""Heun on a doubly periodic 2D array with centered ``numpy.roll`` differences.

    Written out here rather than imported, so that the two-dimensional symbol
    of the module is checked against an oracle the module does not provide.
    The right-hand side is

    .. math::

       -c_r\frac{q_{i+1,j}-q_{i-1,j}}{2\Delta r}
       -c_z\frac{q_{i,j+1}-q_{i,j-1}}{2\Delta z}
       +\nu\left(\frac{q_{i+1,j}-2q_{ij}+q_{i-1,j}}{\Delta r^2}
                 +\frac{q_{i,j+1}-2q_{ij}+q_{i,j-1}}{\Delta z^2}\right),

    and the Heun update is ``q + (dt/2)(k_1 + k_2)`` with ``k_2`` evaluated at
    ``q + dt k_1``.  The frozen ``3\nu/r`` term is passed in through
    ``advection_r`` because it is algebraically just another first-derivative
    coefficient.
    """

    def rhs(field: np.ndarray) -> np.ndarray:
        r_forward = np.roll(field, -1, axis=0)
        r_backward = np.roll(field, 1, axis=0)
        z_forward = np.roll(field, -1, axis=1)
        z_backward = np.roll(field, 1, axis=1)
        return (
            -advection_r * (r_forward - r_backward) / (2.0 * dr)
            - advection_z * (z_forward - z_backward) / (2.0 * dz)
            + viscosity
            * (
                (r_forward - 2.0 * field + r_backward) / dr**2
                + (z_forward - 2.0 * field + z_backward) / dz**2
            )
        )

    current = np.array(values, dtype=np.complex128, copy=True)
    for _ in range(n_steps):
        k1 = rhs(current)
        k2 = rhs(current + dt * k1)
        current = current + 0.5 * dt * (k1 + k2)
    return current


def _mode_1d(n_points: int, wavenumber: int) -> np.ndarray:
    """Return ``e^{2 pi i k j / N}`` sampled on a periodic grid."""

    index = np.arange(n_points, dtype=np.float64)
    return np.exp(2j * np.pi * wavenumber * index / n_points)


def _mode_2d(shape: tuple[int, int], wavenumbers: tuple[int, int]) -> np.ndarray:
    """Return the separable 2D mode ``e^{i(theta_r i + theta_z j)}``."""

    n_r, n_z = shape
    k_r, k_z = wavenumbers
    row = np.arange(n_r, dtype=np.float64)[:, None]
    column = np.arange(n_z, dtype=np.float64)[None, :]
    return np.exp(2j * np.pi * (k_r * row / n_r + k_z * column / n_z))


def _project(values: np.ndarray, wavenumber: int) -> complex:
    """Return the discrete Fourier coefficient of one mode by direct projection."""

    n_points = values.size
    index = np.arange(n_points, dtype=np.float64)
    kernel = np.exp(-2j * np.pi * wavenumber * index / n_points)
    return complex(np.sum(values * kernel) / n_points)


# ---------------------------------------------------------------------------
# 1.  The Heun identity on the imaginary axis.
# ---------------------------------------------------------------------------


def test_heun_amplification_on_the_imaginary_axis_is_one_plus_quartic() -> None:
    r"""``|G_heun(ia)|^2 = 1 + a^4/4`` exactly, hence ``> 1`` for every ``a \neq 0``.

    Derivation, written out here so the test does not depend on the module for
    its oracle: ``G(ia) = 1 + ia + (ia)^2/2 = (1 - a^2/2) + ia`` and

        ``|G|^2 = (1 - a^2/2)^2 + a^2 = 1 - a^2 + a^4/4 + a^2 = 1 + a^4/4``.

    Centered differences give a purely imaginary advective symbol, so this is
    the whole reason Heun + centered + ``nu = 0`` cannot be stable.

    Measured maximum relative defect on this machine: 4.44e-16 (2 ulp).
    """

    a = np.array(
        [1.0e-3, 1.0e-2, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, np.pi],
        dtype=np.float64,
    )
    magnitude_squared = np.abs(stability_polynomial("heun", 1j * a)) ** 2
    closed_form = 1.0 + a**4 / 4.0

    assert np.max(np.abs(magnitude_squared - closed_form) / closed_form) <= 1.0e-14
    # Strictly greater than one wherever the growth is representable at all.
    assert np.all(magnitude_squared > 1.0)

    # The predictor stage is worse still: |1 + ia|^2 = 1 + a^2.
    predictor_squared = np.abs(stability_polynomial("euler", 1j * a)) ** 2
    assert np.max(np.abs(predictor_squared - (1.0 + a**2)) / (1.0 + a**2)) <= 1.0e-14

    # The polynomials are the ones named, checked coefficient by coefficient
    # against monomial sums formed here.
    z = np.array([0.3 + 0.4j, -1.1 + 0.2j, 0.0 + 1.7j, -0.5 - 0.5j])
    expectations = {
        "euler": 1.0 + z,
        "heun": 1.0 + z + z**2 / 2.0,
        "ssprk3": 1.0 + z + z**2 / 2.0 + z**3 / 6.0,
        "rk4": 1.0 + z + z**2 / 2.0 + z**3 / 6.0 + z**4 / 24.0,
    }
    for name, expected in expectations.items():
        computed = stability_polynomial(name, z)
        assert np.max(np.abs(computed - expected)) <= 1.0e-14 * np.max(
            np.abs(expected)
        )


# ---------------------------------------------------------------------------
# 2.  The mandated "detect Heun + centered growth at nu = 0" test.
# ---------------------------------------------------------------------------


def test_heun_amplifies_pure_advection_at_every_representable_cfl() -> None:
    r"""No viscosity: the Heun scan exceeds one at every nonzero CFL, and the
    array propagator reproduces the analytic factor.

    The scan maximum must equal ``sqrt(1 + CFL^4/4)`` because the maximizer
    ``\theta_z = \pm\pi/2`` (where ``|\sin| = 1``) is an exact sample of the
    grid and ``advection_r = 0`` removes every radial contribution.

    Honest floating-point limit: the growth is ``1 + CFL^4/8`` to leading
    order, so for ``CFL`` below about ``5.5e-4`` it falls under half an ulp of
    one and the strict inequality stops being representable in binary64.  The
    ladder therefore starts at ``1e-3``, where the growth is ``1.25e-13``,
    about 560 ulp.  This is a limitation of the *measurement*, not of the
    statement: the exact ``|G|^2 - 1 = a^4/4`` is positive for every ``a != 0``.

    Measured one-step propagator-versus-symbol agreement: 3.72e-16 relative,
    on a mode whose analytic one-step growth is ``|G| - 1 = 1.528e-3``.
    """

    dz = 0.05
    speed = 1.3
    for cfl in (1.0e-3, 1.0e-2, 0.05, 0.1, 0.25, 0.5, 1.0):
        dt = cfl * dz / speed
        result = scan_amplification(
            methods=("heun",),
            dt=dt,
            dr=0.02,
            dz=dz,
            advection_r=0.0,
            advection_z=speed,
            viscosity=0.0,
        )
        peak = result["methods"]["heun"]["max_amplification"]
        assert peak > 1.0
        assert result["advective_cfl_z"] == pytest.approx(cfl, rel=1.0e-14)
        assert peak == pytest.approx(math.sqrt(1.0 + cfl**4 / 4.0), rel=1.0e-14)
        assert abs(result["methods"]["heun"]["argmax_theta_z"]) == pytest.approx(
            0.5 * np.pi, rel=1.0e-14
        )
        # The predictor stage is |1 + i CFL| = sqrt(1 + CFL^2), much worse.
        assert result["methods"]["heun"][
            "predictor_max_amplification"
        ] == pytest.approx(math.sqrt(1.0 + cfl**2), rel=1.0e-14)

    # The array path must reproduce the symbol path, per step, to roundoff.
    n_points = 32
    dx = 1.0 / n_points
    wavenumber = 5
    dt = 0.4 * dx / speed
    theta = 2.0 * np.pi * wavenumber / n_points
    symbol = complex(
        advection_diffusion_symbol(
            0.0,
            theta,
            advection_r=0.0,
            advection_z=speed,
            viscosity=0.0,
            dr=dx,
            dz=dx,
        )
    )
    analytic = complex(stability_polynomial("heun", dt * symbol))
    assert abs(analytic) > 1.0

    initial = _mode_1d(n_points, wavenumber)
    stepped = reference_propagate(
        "heun",
        initial,
        dt=dt,
        dx=dx,
        advection=speed,
        viscosity=0.0,
        n_steps=1,
    )
    ratio = complex(stepped[0] / initial[0])
    assert abs(ratio - analytic) <= 1.0e-13 * abs(analytic)
    # Growth of the whole array, not just of one entry.
    growth = float(np.max(np.abs(stepped))) / float(np.max(np.abs(initial)))
    assert growth == pytest.approx(abs(analytic), rel=1.0e-13)
    assert growth > 1.0

    # The scalar-mode reference path must agree with the array path as well.
    mode_amplitude = reference_mode_amplification(
        "heun",
        wavenumber_index=wavenumber,
        n_points=n_points,
        dt=dt,
        advection=speed,
        viscosity=0.0,
        n_steps=1,
    )
    assert abs(mode_amplitude - analytic) <= 1.0e-13 * abs(analytic)


# ---------------------------------------------------------------------------
# 3.  SSPRK3 and RK4 inside their imaginary-axis intervals.
# ---------------------------------------------------------------------------


def test_ssprk3_and_rk4_are_stable_on_pure_advection_at_cfl_one() -> None:
    r"""At ``dt|c|/dz = 1`` neither three- nor four-stage method amplifies.

    The claim is *asserted from the scan*, not hardcoded: the scan is run and
    its maximum is required to be ``<= 1 + 1e-12``.  Measured: both are exactly
    ``1.0``, attained at the sawtooth ``\theta_z = -\pi`` where the centered
    advective symbol vanishes.  Heun on the same data measures
    ``1.1180339887498949 = sqrt(1 + 1/4)``, i.e. it fails where they pass.

    The interval endpoints ``\sqrt3`` and ``2\sqrt2`` are quoted in the module
    docstring; they are recomputed here from the closed forms

        ``|G_{ssprk3}(ia)|^2 = 1 - a^4/12 + a^6/36``,
        ``|G_{rk4}(ia)|^2   = 1 - a^6/72 + a^8/576``,

    which are expanded in this file rather than taken from the module.
    """

    dz = 0.05
    speed = 1.0
    dt = 1.0 * dz / speed
    result = scan_amplification(
        methods=("ssprk3", "rk4", "heun"),
        dt=dt,
        dr=0.02,
        dz=dz,
        advection_r=0.0,
        advection_z=speed,
        viscosity=0.0,
    )
    assert result["advective_cfl_z"] == pytest.approx(1.0, rel=1.0e-14)
    assert result["advective_cfl_r"] == 0.0

    for name in ("ssprk3", "rk4"):
        assert result["methods"][name]["max_amplification"] <= 1.0 + 1.0e-12
    # Non-vacuous: the method actually shipped in production fails here.
    assert result["methods"]["heun"]["max_amplification"] > 1.1

    # Independent recomputation of the two interval endpoints.
    a = np.linspace(0.0, math.sqrt(3.0), 401)
    ssprk3_squared = np.abs(stability_polynomial("ssprk3", 1j * a)) ** 2
    assert np.max(np.abs(ssprk3_squared - (1.0 - a**4 / 12.0 + a**6 / 36.0))) <= 1.0e-14
    assert np.max(ssprk3_squared) <= 1.0 + 1.0e-14
    assert abs(ssprk3_squared[-1] - 1.0) <= 1.0e-14
    assert np.abs(stability_polynomial("ssprk3", 1j * 1.05 * math.sqrt(3.0))) > 1.0

    a = np.linspace(0.0, 2.0 * math.sqrt(2.0), 401)
    rk4_squared = np.abs(stability_polynomial("rk4", 1j * a)) ** 2
    assert np.max(np.abs(rk4_squared - (1.0 - a**6 / 72.0 + a**8 / 576.0))) <= 1.0e-13
    assert np.max(rk4_squared) <= 1.0 + 1.0e-14
    assert abs(rk4_squared[-1] - 1.0) <= 1.0e-13
    assert np.abs(stability_polynomial("rk4", 1j * 1.05 * 2.0 * math.sqrt(2.0))) > 1.0


# ---------------------------------------------------------------------------
# 4.  Advection-diffusion: symbol versus a 2D propagator, and the viscosity
#     threshold at which the Heun scan starts to pass.
# ---------------------------------------------------------------------------


def test_two_dimensional_symbol_matches_an_independent_array_propagator() -> None:
    r"""``|G(dt \lambda)|`` equals the measured per-step ratio of a 2D Heun run.

    The oracle is ``_propagate_2d_heun`` in this file, which never forms a
    symbol.  The frozen radial first-derivative coefficient ``\kappa_r =
    3\nu/\Delta r`` is included, so the term that models the ``3\nu/r`` part of
    ``\mathcal L_5`` is validated too -- as an *interior* frozen coefficient.
    It remains true that the real axis row of E-17 is a boundary modification
    outside the reach of this or any von Neumann argument.

    Measured agreement on this machine: 5.68e-17 relative for a single step and
    9.72e-16 after sixteen steps.
    """

    n_r, n_z = 16, 20
    k_r, k_z = 3, 7
    dr, dz = 0.05, 0.04
    advection_r, advection_z = 0.7, 1.3
    viscosity = 1.0e-3
    dt = 0.01
    kappa_r = 3.0 * viscosity / dr

    theta_r = 2.0 * np.pi * k_r / n_r
    theta_z = 2.0 * np.pi * k_z / n_z
    symbol = complex(
        advection_diffusion_symbol(
            theta_r,
            theta_z,
            advection_r=advection_r,
            advection_z=advection_z,
            viscosity=viscosity,
            dr=dr,
            dz=dz,
            first_derivative_coefficient_r=kappa_r,
        )
    )
    # The symbol must be genuinely two-dimensional and genuinely damped.
    assert symbol.real < 0.0
    assert abs(symbol.imag) > 0.0

    initial = _mode_2d((n_r, n_z), (k_r, k_z))
    for n_steps in (1, 16):
        stepped = _propagate_2d_heun(
            initial,
            dt=dt,
            dr=dr,
            dz=dz,
            advection_r=advection_r + kappa_r,
            advection_z=advection_z,
            viscosity=viscosity,
            n_steps=n_steps,
        )
        measured = complex(stepped[0, 0] / initial[0, 0])
        analytic = complex(stability_polynomial("heun", dt * symbol)) ** n_steps
        assert abs(measured - analytic) <= 1.0e-13 * abs(analytic)
        assert abs(abs(measured) - abs(analytic)) <= 1.0e-13 * abs(analytic)

    # The additive structure the two-dimensional symbol relies on: the radial
    # and axial contributions are exactly the two one-dimensional symbols.
    radial = complex(
        advection_diffusion_symbol(
            theta_r,
            0.0,
            advection_r=advection_r,
            advection_z=advection_z,
            viscosity=viscosity,
            dr=dr,
            dz=dz,
            first_derivative_coefficient_r=kappa_r,
        )
    )
    axial = complex(
        advection_diffusion_symbol(
            0.0,
            theta_z,
            advection_r=advection_r,
            advection_z=advection_z,
            viscosity=viscosity,
            dr=dr,
            dz=dz,
            first_derivative_coefficient_r=kappa_r,
        )
    )
    assert abs((radial + axial) - symbol) <= 1.0e-15 * abs(symbol)

    # And the axial half is reproduced by the module's own 1D array path.
    n_points = n_z
    stepped_1d = reference_propagate(
        "heun",
        _mode_1d(n_points, k_z),
        dt=dt,
        dx=dz,
        advection=advection_z,
        viscosity=viscosity,
        n_steps=1,
    )
    ratio_1d = complex(stepped_1d[0] / _mode_1d(n_points, k_z)[0])
    assert abs(ratio_1d - complex(stability_polynomial("heun", dt * axial))) <= (
        1.0e-13 * abs(ratio_1d)
    )


def test_enough_viscosity_makes_the_heun_scan_pass_and_zero_viscosity_fails() -> None:
    r"""There is a viscosity at which the same ``dt`` becomes acceptable.

    The threshold is *searched for*, not assumed: a geometric ladder in ``nu``
    is scanned and the first passing value is used.  On this machine the ladder
    first passes at ``nu = 4.217e-4`` with ``dt = 0.025``, ``dr = dz = 0.05``,
    ``c_z = 1``, i.e. advective CFL ``0.5`` and viscous number ``0.016868``;
    the same ``dt`` at ``nu = 0`` measures ``1.007782218537``.

    The pass is not vacuous in the other direction either: at large enough
    ``nu`` the *parabolic* limit is violated and the scan fails again at the
    sawtooth, which the tail of this test checks.
    """

    dr = dz = 0.05
    speed = 1.0
    dt = 0.5 * dz / speed

    inviscid = scan_amplification(
        methods=("heun",),
        dt=dt,
        dr=dr,
        dz=dz,
        advection_r=0.0,
        advection_z=speed,
        viscosity=0.0,
    )
    assert inviscid["methods"]["heun"]["max_amplification"] > 1.0 + 1.0e-12
    assert inviscid["advective_cfl_z"] == pytest.approx(0.5, rel=1.0e-14)

    ladder = [1.0e-5 * (10.0 ** (index / 8.0)) for index in range(25)]
    passing: float | None = None
    for viscosity in ladder:
        result = scan_amplification(
            methods=("heun",),
            dt=dt,
            dr=dr,
            dz=dz,
            advection_r=0.0,
            advection_z=speed,
            viscosity=viscosity,
        )
        if result["methods"]["heun"]["max_amplification"] <= 1.0 + 1.0e-12:
            passing = viscosity
            break
    assert passing is not None, "no viscosity on the ladder made the Heun scan pass"
    assert 1.0e-4 < passing < 1.0e-2

    audit = audit_snapshot(
        max_abs_u_r=0.0,
        max_abs_u_z=speed,
        dr=dr,
        dz=dz,
        dt=dt,
        viscosity=passing,
    )
    assert audit["passes"] is True
    assert audit["verdict"] == "stability-verified-in-model"
    assert audit["methods"]["heun"]["passes"] is True

    failing = audit_snapshot(
        max_abs_u_r=0.0,
        max_abs_u_z=speed,
        dr=dr,
        dz=dz,
        dt=dt,
        viscosity=0.0,
    )
    assert failing["passes"] is False
    assert failing["verdict"] == "stability-unverified"

    # Too much viscosity is also rejected: the real part leaves [-2, 0] and the
    # sawtooth mode is the one that escapes.  Measured at nu = 0.03: 1.48.
    over_damped = scan_amplification(
        methods=("heun",),
        dt=dt,
        dr=dr,
        dz=dz,
        advection_r=0.0,
        advection_z=speed,
        viscosity=0.03,
    )
    assert over_damped["methods"]["heun"]["max_amplification"] > 1.0 + 1.0e-12
    assert over_damped["methods"]["heun"]["argmax_theta_z"] == -np.pi
    assert over_damped["viscous_number"] > 1.0


# ---------------------------------------------------------------------------
# 5.  Multi-mode long-time propagation against the semi-discrete solution.
# ---------------------------------------------------------------------------


def test_multi_mode_long_time_rk4_tracks_the_semi_discrete_solution() -> None:
    r"""1000 RK4 steps of a three-mode field match ``e^{\lambda_k t}`` per mode.

    The reference is the exact solution of the **semi-discrete** system -- the
    discrete symbol ``\lambda_k`` of the centered stencil, not the PDE symbol
    ``-ick - \nu k^2``.  Using the PDE symbol would fold the ``O(\Delta x^2)``
    spatial truncation error into the comparison and make it insensitive to the
    time integration, which is the thing under test here.

    Measured on this machine at CFL ``0.008``: maximum relative field error
    1.99e-11, per-mode relative amplitude error at most 4.73e-12 and per-mode
    phase error at most 1.46e-10 radians.  The residual is the
    ``O(|dt\lambda|^5)`` per-step defect of RK4 accumulated over 1000 steps, so
    it shrinks like ``dt^4`` and is not roundoff.

    The tail of the test then runs Heun at ``\nu = 0`` with the *same* ``dt``
    and confirms the amplitude grows by exactly the predicted ``|G|^n``.
    """

    n_points = 64
    dx = 1.0 / n_points
    speed = 1.0
    viscosity = 1.0e-4
    cfl = 0.008
    dt = cfl * dx / speed
    n_steps = 1000
    modes = ((1, 0.7), (5, -0.3), (11, 0.15))

    initial = np.zeros(n_points, dtype=np.complex128)
    for wavenumber, amplitude in modes:
        initial = initial + amplitude * _mode_1d(n_points, wavenumber)

    evolved = reference_propagate(
        "rk4",
        initial,
        dt=dt,
        dx=dx,
        advection=speed,
        viscosity=viscosity,
        n_steps=n_steps,
    )

    expected = np.zeros(n_points, dtype=np.complex128)
    for wavenumber, amplitude in modes:
        theta = 2.0 * np.pi * wavenumber / n_points
        symbol = complex(
            advection_diffusion_symbol(
                0.0,
                theta,
                advection_r=0.0,
                advection_z=speed,
                viscosity=viscosity,
                dr=dx,
                dz=dx,
            )
        )
        expected = expected + (
            amplitude
            * np.exp(symbol * dt * n_steps)
            * _mode_1d(n_points, wavenumber)
        )

    field_error = float(np.max(np.abs(evolved - expected))) / float(
        np.max(np.abs(expected))
    )
    assert field_error < 1.0e-8

    decay_ratios: list[float] = []
    for wavenumber, amplitude in modes:
        theta = 2.0 * np.pi * wavenumber / n_points
        symbol = complex(
            advection_diffusion_symbol(
                0.0,
                theta,
                advection_r=0.0,
                advection_z=speed,
                viscosity=viscosity,
                dr=dx,
                dz=dx,
            )
        )
        reference = amplitude * complex(np.exp(symbol * dt * n_steps))
        measured = _project(evolved, wavenumber)
        # Amplitude and phase are checked separately, as required: a pure phase
        # error would be invisible to an amplitude test and vice versa.
        assert abs(abs(measured) - abs(reference)) <= 1.0e-8 * abs(reference)
        phase_error = abs(np.angle(measured / reference))
        assert phase_error < 1.0e-8
        # Every mode decays, and every mode has been carried a long way in
        # phase, so neither half of the comparison is against a fixed point.
        assert abs(reference) < abs(amplitude)
        assert abs(symbol.imag) * dt * n_steps > 0.5
        decay_ratios.append(abs(reference) / abs(amplitude))
    # The shortest resolved mode of the three has decayed appreciably; the
    # longest one barely does, which is the physically correct ordering.
    # Measured: 0.999507 (k=1), 0.987982 (k=5), 0.947310 (k=11).
    assert decay_ratios[-1] < 0.96
    assert decay_ratios[0] > decay_ratios[1] > decay_ratios[2]

    # Heun with nu = 0 at the same dt: growth exactly |G|^n.
    wavenumber = 11
    theta = 2.0 * np.pi * wavenumber / n_points
    symbol = complex(
        advection_diffusion_symbol(
            0.0,
            theta,
            advection_r=0.0,
            advection_z=speed,
            viscosity=0.0,
            dr=dx,
            dz=dx,
        )
    )
    growth_factor = abs(complex(stability_polynomial("heun", dt * symbol)))
    predicted = growth_factor**n_steps
    single = _mode_1d(n_points, wavenumber)
    grown = reference_propagate(
        "heun",
        single,
        dt=dt,
        dx=dx,
        advection=speed,
        viscosity=0.0,
        n_steps=n_steps,
    )
    measured_growth = float(np.max(np.abs(grown)))
    # Measured: 1.0000003097342918 against a prediction of 1.0000003097342862,
    # i.e. 5.6e-15 relative.  The growth itself is 3.1e-7, about 1.4 million
    # ulp, so it is far above the noise floor of the measurement.
    assert abs(measured_growth - predicted) <= 1.0e-10 * predicted
    assert measured_growth > 1.0
    assert predicted - 1.0 > 1.0e-8

    # The scalar-mode path agrees with the array path over 1000 steps as well.
    scalar = reference_mode_amplification(
        "heun",
        wavenumber_index=wavenumber,
        n_points=n_points,
        dt=dt,
        advection=speed,
        viscosity=0.0,
        n_steps=n_steps,
    )
    assert abs(abs(scalar) - measured_growth) <= 1.0e-10 * measured_growth


# ---------------------------------------------------------------------------
# 6.  The wavenumber grid and its refinement behavior.
# ---------------------------------------------------------------------------


def test_theta_grid_contains_pi_and_the_scan_is_monotone_under_refinement() -> None:
    r"""``\pm\pi`` and ``0`` are exact samples, and refining can only raise the max.

    ``theta_grid(181)`` is a bitwise subset of ``theta_grid(721)`` (the step
    halves twice, and binary64 division by a power of two is exact), so the
    sampled maximum is monotone nondecreasing by construction rather than by
    luck.  Both facts are asserted.

    Case with an interior maximizer: pure advection in both directions, where
    the maximizer is ``\theta_r = \theta_z = -\pi/2`` -- interior, and an exact
    sample of both grids, so the two scans agree to the bit (measured relative
    difference 0.0, well inside the required 1e-6).

    Second case, reported honestly: with ``nu = 1e-3`` the maximizer moves
    *off* the coarse grid, to ``(-1.3875, -1.4399)``.  The sampled maximum then
    converges only at ``O(\Delta\theta^2)`` and the measured 181-vs-721 gap is
    5.6e-6 relative -- larger than 1e-6.  That is a property of sampling a
    smooth maximum, not a defect, and it is bounded here rather than hidden.
    """

    coarse_grid = theta_grid(181)
    fine_grid = theta_grid(721)
    for grid in (coarse_grid, fine_grid):
        assert grid[0] == -np.pi
        assert grid[-1] == np.pi
        assert np.any(grid == 0.0)
        assert np.all(np.diff(grid) > 0.0)
    assert bool(np.all(np.isin(coarse_grid, fine_grid)))

    common = dict(
        methods=("heun",),
        dt=0.02,
        dr=0.05,
        dz=0.04,
        advection_r=0.7,
        advection_z=1.3,
    )
    coarse = scan_amplification(**common, viscosity=0.0, n_theta=181)["methods"][
        "heun"
    ]
    fine = scan_amplification(**common, viscosity=0.0, n_theta=721)["methods"]["heun"]

    # Interior maximizer, not the sawtooth and not the constant mode.
    for entry in (coarse, fine):
        assert 0.0 < abs(entry["argmax_theta_r"]) < np.pi
        assert 0.0 < abs(entry["argmax_theta_z"]) < np.pi
    assert fine["max_amplification"] >= coarse["max_amplification"]
    assert fine["max_amplification"] - coarse["max_amplification"] <= (
        1.0e-6 * fine["max_amplification"]
    )
    assert fine["max_amplification"] > 1.05

    # Off-grid maximizer: monotone, and the gap is the O(dtheta^2) sampling
    # error of a smooth interior maximum.  Measured 5.6e-6 relative.
    coarse_off = scan_amplification(**common, viscosity=1.0e-3, n_theta=181)[
        "methods"
    ]["heun"]
    fine_off = scan_amplification(**common, viscosity=1.0e-3, n_theta=721)["methods"][
        "heun"
    ]
    assert fine_off["max_amplification"] >= coarse_off["max_amplification"]
    assert fine_off["max_amplification"] - coarse_off["max_amplification"] <= (
        1.0e-4 * fine_off["max_amplification"]
    )
    assert 0.0 < abs(fine_off["argmax_theta_r"]) < np.pi

    # The sawtooth is reachable as an argmax, which is the operational proof
    # that theta = pi is in the scan and not merely in the grid helper.
    sawtooth = scan_amplification(
        methods=("heun",),
        dt=0.025,
        dr=0.05,
        dz=0.05,
        advection_r=0.0,
        advection_z=1.0,
        viscosity=0.03,
    )["methods"]["heun"]
    assert sawtooth["argmax_theta_r"] == -np.pi
    assert sawtooth["argmax_theta_z"] == -np.pi

    # An even n_theta would drop theta = 0 and is rejected outright.
    with pytest.raises(ValueError, match="odd"):
        theta_grid(180)


# ---------------------------------------------------------------------------
# 7.  The Hou-like operating point.
# ---------------------------------------------------------------------------


def test_hou_like_operating_point_is_stability_unverified() -> None:
    r"""The honest measured verdict at the prescribed operating point is FAIL, twice.

    Operating point: ``dr = 1/192``, ``dz = 1/384``, ``dt = 2.76e-7``,
    ``max|u^r| = 2000``, ``max|u^z| = 7000``, ``nu = 5e-4``.  This gives
    advective CFL numbers 0.10598 (radial) and 0.74189 (axial), so the
    worst-case ``a = dt(|c_r|/dr + |c_z|/dz) = 0.8479`` and the quartic Heun
    growth ``a^4/4 = 0.1292`` overwhelms the viscous damping, whose viscous
    number is only 8.14e-05.

    Measured, both at ``theta_r = theta_z = -pi/2``:

    ==============  ==================  ==================
    viscosity       Heun ``max |G|``    predictor stage
    ==============  ==================  ==================
    ``5e-4``        1.0625773361391977  1.3110343498529449
    ``0``           1.0626380399507376  1.3110632816092440
    ==============  ==================  ==================

    So the verdict does **not** flip when the viscosity is zeroed: it is
    ``stability-unverified`` either way, and the entire contribution of
    ``nu = 5e-4`` at this step size is a reduction of 6.07e-05 in the
    amplification factor.  This comment records that as the measured verdict
    rather than asserting a flip that does not happen.

    The flip *can* happen -- it is demonstrated on a different operating point
    in ``test_enough_viscosity_makes_the_heun_scan_pass_and_zero_viscosity_fails``
    -- and it can also be reached here by shrinking ``dt``: the audit passes at
    this grid and viscosity for ``dt <= 2.4e-8``, about 11.5x below the
    prescribed step.  That number is recomputed by bisection below.

    Caveat recorded rather than smoothed over: the shipped run
    ``outputs/hou_early_time_v1`` reports ``minimum_dt = 2.7587e-07`` *together
    with* ``maximum_advective_cfl = 0.10023``, which forces
    ``max(|u^r|/dr, |u^z|/dz) = 0.10023/dt``.  The ``max|u^z| = 7000``
    prescribed here is therefore about 7.4x larger than that pairing supports;
    the operating point audited in this test is a pessimistic one, and the
    numbers above should not be read as the audit of that run.

    Nothing here says the run blew up.  It says the frozen-coefficient worst
    case is not verified at these numbers.
    """

    viscous = audit_snapshot(
        max_abs_u_r=HOU_MAX_U_R,
        max_abs_u_z=HOU_MAX_U_Z,
        dr=HOU_DR,
        dz=HOU_DZ,
        dt=HOU_DT,
        viscosity=HOU_VISCOSITY,
        tolerance=1.0e-12,
    )
    inviscid = audit_snapshot(
        max_abs_u_r=HOU_MAX_U_R,
        max_abs_u_z=HOU_MAX_U_Z,
        dr=HOU_DR,
        dz=HOU_DZ,
        dt=HOU_DT,
        viscosity=0.0,
        tolerance=1.0e-12,
    )

    assert viscous["advective_cfl_r"] == pytest.approx(0.105984, rel=1.0e-5)
    assert viscous["advective_cfl_z"] == pytest.approx(0.741888, rel=1.0e-5)
    assert viscous["viscous_number"] < 1.0e-4

    # The measured verdict.  This is a FAIL and is asserted as one.
    assert viscous["passes"] is False
    assert viscous["verdict"] == "stability-unverified"
    assert inviscid["passes"] is False
    assert inviscid["verdict"] == "stability-unverified"

    viscous_peak = viscous["methods"]["heun"]["max_amplification"]
    inviscid_peak = inviscid["methods"]["heun"]["max_amplification"]
    assert viscous_peak == pytest.approx(1.0625773361, rel=1.0e-9)
    assert inviscid_peak == pytest.approx(1.0626380400, rel=1.0e-9)
    # Viscosity helps, by an amount far too small to matter.
    assert viscous_peak < inviscid_peak
    assert 1.0e-5 < inviscid_peak - viscous_peak < 1.0e-4

    for result in (viscous, inviscid):
        entry = result["methods"]["heun"]
        assert entry["argmax_theta_r"] == pytest.approx(-0.5 * np.pi, rel=1.0e-9)
        assert entry["argmax_theta_z"] == pytest.approx(-0.5 * np.pi, rel=1.0e-9)
        # The predictor stage is far worse than the completed step.
        assert entry["predictor_max_amplification"] > 1.31
        assert entry["predictor_max_amplification"] > entry["max_amplification"]

    # The step size at which this same grid and viscosity would pass, found by
    # bisection rather than quoted.  Measured: 2.4015e-08.
    low, high = 1.0e-10, HOU_DT
    assert audit_snapshot(
        max_abs_u_r=HOU_MAX_U_R,
        max_abs_u_z=HOU_MAX_U_Z,
        dr=HOU_DR,
        dz=HOU_DZ,
        dt=low,
        viscosity=HOU_VISCOSITY,
    )["passes"]
    for _ in range(24):
        middle = 0.5 * (low + high)
        if audit_snapshot(
            max_abs_u_r=HOU_MAX_U_R,
            max_abs_u_z=HOU_MAX_U_Z,
            dr=HOU_DR,
            dz=HOU_DZ,
            dt=middle,
            viscosity=HOU_VISCOSITY,
        )["passes"]:
            low = middle
        else:
            high = middle
    assert 2.0e-8 < low < 3.0e-8
    assert HOU_DT / low > 9.0


def test_recorded_hou_cfl_pairing_is_borderline_rather_than_comfortable() -> None:
    r"""The audit at the *self-consistent* operating point of the shipped run.

    ``outputs/hou_early_time_v1/summary.json`` records, for the ``193 x 384``
    resolution, ``minimum_dt = 2.7586570871878865e-07`` together with
    ``maximum_advective_cfl = 0.10022676304450114``.  Those two numbers pin the
    dominant velocity ratio exactly:
    ``max(|u^r|/dr, |u^z|/dz) = maximum_advective_cfl / dt``.  What they do
    *not* say is which direction attains it, so all three readings are audited
    here.  The velocities are derived from the recorded pair rather than
    assumed, and no file is read: the two constants are transcribed above.

    Measured Heun worst-case amplification at ``nu = 5e-4``:

    ==========================  ==================  ======
    reading                     ``max |G|``         pass
    ==========================  ==================  ======
    radial direction binding    1.000003460085721   no
    axial direction binding     1.000000000000000   yes
    both at the recorded CFL    1.000152277605197   no
    ==========================  ==================  ======

    The point of the table is that the margin is *thin* -- at the shipped CFL
    the frozen-coefficient audit is decided in the fifth or sixth decimal, and
    which way it falls depends on which direction carries the CFL.  That is a
    reason to record the audit alongside a run, which is what this module is
    for.  It is not a claim that the run is wrong: see the module docstring for
    everything this model omits.
    """

    recorded_dt = 2.7586570871878865e-07
    recorded_cfl = 0.10022676304450114
    ratio = recorded_cfl / recorded_dt
    radial_speed = ratio * HOU_DR
    axial_speed = ratio * HOU_DZ
    assert radial_speed == pytest.approx(1892.3, rel=1.0e-3)
    assert axial_speed == pytest.approx(946.1, rel=1.0e-3)

    readings = {
        "radial": (radial_speed, 0.0),
        "axial": (0.0, axial_speed),
        "both": (radial_speed, axial_speed),
    }
    peaks = {}
    for label, (speed_r, speed_z) in readings.items():
        result = audit_snapshot(
            max_abs_u_r=speed_r,
            max_abs_u_z=speed_z,
            dr=HOU_DR,
            dz=HOU_DZ,
            dt=recorded_dt,
            viscosity=HOU_VISCOSITY,
        )
        peaks[label] = result["methods"]["heun"]["max_amplification"]
        assert result["passes"] is (
            result["methods"]["heun"]["max_amplification"] <= 1.0 + 1.0e-12
        )

    # Every reading is within 2e-4 of the acceptance threshold in either
    # direction: the audit is decided by a thin margin here, unlike the
    # prescribed operating point of the previous test where it fails by 6e-2.
    for peak in peaks.values():
        assert abs(peak - 1.0) < 2.0e-4
    assert peaks["axial"] <= 1.0 + 1.0e-12
    assert peaks["radial"] > 1.0 + 1.0e-12
    assert peaks["both"] > peaks["radial"] > peaks["axial"]


# ---------------------------------------------------------------------------
# 8.  Validation.
# ---------------------------------------------------------------------------


def test_rejects_invalid_inputs() -> None:
    """Every guard the module documents is exercised."""

    with pytest.raises(ValueError, match="unknown method"):
        stability_polynomial("rk3", 0.1j)
    with pytest.raises(ValueError, match="unknown method"):
        stability_polynomial("HEUN", 0.1j)
    with pytest.raises(TypeError, match="method must be a string"):
        stability_polynomial(2, 0.1j)
    with pytest.raises(ValueError, match="finite"):
        stability_polynomial("heun", complex(np.nan, 0.0))

    with pytest.raises(ValueError, match="unknown method"):
        scan_amplification(
            methods=("heun", "leapfrog"),
            dt=1.0e-3,
            dr=0.1,
            dz=0.1,
            advection_r=1.0,
            advection_z=1.0,
            viscosity=1.0e-3,
        )
    with pytest.raises(ValueError, match="at least one"):
        scan_amplification(
            methods=(),
            dt=1.0e-3,
            dr=0.1,
            dz=0.1,
            advection_r=1.0,
            advection_z=1.0,
            viscosity=1.0e-3,
        )
    with pytest.raises(TypeError, match="not a single string"):
        scan_amplification(
            methods="heun",
            dt=1.0e-3,
            dr=0.1,
            dz=0.1,
            advection_r=1.0,
            advection_z=1.0,
            viscosity=1.0e-3,
        )

    base = dict(
        methods=("heun",),
        dt=1.0e-3,
        dr=0.1,
        dz=0.1,
        advection_r=1.0,
        advection_z=1.0,
        viscosity=1.0e-3,
    )
    for name in ("dt", "dr", "dz"):
        for bad in (0.0, -1.0, math.inf, math.nan):
            arguments = dict(base)
            arguments[name] = bad
            with pytest.raises(ValueError, match=f"{name} must be positive"):
                scan_amplification(**arguments)
    for bad in (-1.0e-12, -1.0, math.nan):
        arguments = dict(base)
        arguments["viscosity"] = bad
        with pytest.raises(ValueError, match="viscosity must be nonnegative"):
            scan_amplification(**arguments)
    for bad in (math.inf, math.nan):
        arguments = dict(base)
        arguments["advection_z"] = bad
        with pytest.raises(ValueError, match="advection_z must be finite"):
            scan_amplification(**arguments)
    with pytest.raises(ValueError, match="odd"):
        scan_amplification(**base, n_theta=720)
    with pytest.raises(ValueError, match="at least 3"):
        scan_amplification(**base, n_theta=1)
    with pytest.raises(TypeError, match="n_theta must be an integer"):
        scan_amplification(**base, n_theta=721.0)
    with pytest.raises(TypeError, match="include_radial_first_derivative"):
        scan_amplification(**base, include_radial_first_derivative=1)

    with pytest.raises(ValueError, match="dr must be positive"):
        advection_diffusion_symbol(
            0.1,
            0.2,
            advection_r=1.0,
            advection_z=1.0,
            viscosity=1.0e-3,
            dr=0.0,
            dz=0.1,
        )
    with pytest.raises(ValueError, match="viscosity must be nonnegative"):
        advection_diffusion_symbol(
            0.1,
            0.2,
            advection_r=1.0,
            advection_z=1.0,
            viscosity=-1.0e-9,
            dr=0.1,
            dz=0.1,
        )
    with pytest.raises(ValueError, match="finite"):
        advection_diffusion_symbol(
            math.nan,
            0.2,
            advection_r=1.0,
            advection_z=1.0,
            viscosity=1.0e-3,
            dr=0.1,
            dz=0.1,
        )

    snapshot = dict(
        max_abs_u_r=1.0,
        max_abs_u_z=1.0,
        dr=0.1,
        dz=0.1,
        dt=1.0e-3,
        viscosity=1.0e-3,
    )
    for name in ("max_abs_u_r", "max_abs_u_z", "tolerance"):
        arguments = dict(snapshot)
        arguments[name] = -1.0
        with pytest.raises(ValueError, match=f"{name} must be nonnegative"):
            audit_snapshot(**arguments)
    arguments = dict(snapshot)
    arguments["viscosity"] = -1.0
    with pytest.raises(ValueError, match="viscosity must be nonnegative"):
        audit_snapshot(**arguments)

    with pytest.raises(ValueError, match="one-dimensional"):
        reference_propagate(
            "heun",
            np.zeros((4, 4)),
            dt=1.0e-3,
            dx=0.1,
            advection=1.0,
            viscosity=0.0,
            n_steps=1,
        )
    with pytest.raises(ValueError, match="at least three points"):
        reference_propagate(
            "heun",
            np.zeros(2),
            dt=1.0e-3,
            dx=0.1,
            advection=1.0,
            viscosity=0.0,
            n_steps=1,
        )
    with pytest.raises(ValueError, match="finite"):
        reference_propagate(
            "heun",
            np.array([0.0, np.nan, 1.0, 2.0]),
            dt=1.0e-3,
            dx=0.1,
            advection=1.0,
            viscosity=0.0,
            n_steps=1,
        )
    with pytest.raises(ValueError, match="n_steps must be at least 1"):
        reference_propagate(
            "heun",
            np.zeros(8),
            dt=1.0e-3,
            dx=0.1,
            advection=1.0,
            viscosity=0.0,
            n_steps=0,
        )
    with pytest.raises(ValueError, match="dx must be positive"):
        reference_propagate(
            "heun",
            np.zeros(8),
            dt=1.0e-3,
            dx=-0.1,
            advection=1.0,
            viscosity=0.0,
            n_steps=1,
        )
    with pytest.raises(ValueError, match="wavenumber_index"):
        reference_mode_amplification(
            "heun",
            wavenumber_index=9,
            n_points=16,
            dt=1.0e-3,
            advection=1.0,
            viscosity=0.0,
            n_steps=1,
        )
    with pytest.raises(ValueError, match="unknown method"):
        reference_mode_amplification(
            "midpoint",
            wavenumber_index=1,
            n_points=16,
            dt=1.0e-3,
            advection=1.0,
            viscosity=0.0,
            n_steps=1,
        )

    # A violently unstable request must raise rather than return infinities.
    with pytest.raises(FloatingPointError, match="non-finite"):
        reference_propagate(
            "euler",
            _mode_1d(16, 4),
            dt=1.0e4,
            dx=1.0e-3,
            advection=1.0e3,
            viscosity=1.0,
            n_steps=4000,
        )


# ---------------------------------------------------------------------------
# 9.  Independence of the module.
# ---------------------------------------------------------------------------


def test_module_is_numpy_only_and_imports_no_production_solver() -> None:
    """The audit must not be able to launder a solver defect into a pass.

    The module is a self-contained piece of Fourier algebra; it imports no
    production module, so a fault injected into the solver cannot change any
    number it reports, and it uses no SciPy, in line with the rest of the
    package.  The check is AST-based plus a source-text scan for transforms,
    because a transform needs no import statement of its own.
    """

    source = inspect.getsource(von_neumann_module)
    tree = ast.parse(source)
    imports: set[tuple[int, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update((0, alias.name) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add((node.level, node.module or ""))

    assert imports <= {
        (0, "__future__"),
        (0, "dataclasses"),
        (0, "math"),
        (0, "typing"),
        (0, "numpy"),
        (0, "numpy.typing"),
    }

    imported_names = {name for _, name in imports}
    forbidden = {
        "scipy",
        "grid",
        "operators",
        "pde",
        "poisson",
        "finite_cylinder_poisson",
        "realspace_poisson",
        "nonlinear_cylinder",
        "cartesian",
    }
    for name in forbidden:
        assert name not in imported_names
        assert f"ns_certificate_lab.{name}" not in imported_names

    text = Path(von_neumann_module.__file__).resolve().read_text(encoding="utf-8")
    assert "fft" not in text.lower()
    # The array reference path really is the rolled centered difference.
    assert "np.roll" in text

    # The reference propagators must not be implemented through the stability
    # polynomial, or they could not falsify it.  Checked structurally: no call
    # to ``stability_polynomial`` appears inside their function bodies.
    bodies = {
        node.name: node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef)
    }
    for name in ("_advance", "reference_propagate", "reference_mode_amplification"):
        assert name in bodies
        called = {
            child.func.id
            for child in ast.walk(bodies[name])
            if isinstance(child, ast.Call) and isinstance(child.func, ast.Name)
        }
        assert "stability_polynomial" not in called
