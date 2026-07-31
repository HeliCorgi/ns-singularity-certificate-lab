"""Tests for the gradient-form pressure pairing certificate.

The certified object is the box integral ``-3 int |u| u . grad p_h dx`` with
``p_h`` the multilinear interpolant of the computed nodal pressure; everything
about the true PDE pressure stays behind the named unproved hypothesis P1G.
Nothing in this file asserts, and nothing in the code under test claims,
anything about the Clay problem.
"""

from __future__ import annotations

import copy
from fractions import Fraction
import math

import numpy as np
import pytest

from ns_certificate_lab.gaussian_gradient_certificate import (
    GRADIENT_CERTIFICATE_KIND,
    CellPressureGradient,
    build_gradient_generation_certificate,
    compare_certificate_paths,
    pressure_cell_gradient,
    verify_gradient_generation_certificate,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.interval_gaussian import (
    GaussianFamily,
    cell_velocity_enclosure,
    gaussian_reference_field,
    seed_gaussian_family,
)
from ns_certificate_lab.l3_generation import AxisymmetricPressureSolver
from ns_certificate_lab.snapshot_certificate import Interval


def _grid(nr: int = 9, nz: int = 17) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=3.0, z_min=-3.0, z_max=3.0, periodic_z=False
    )


def _box(lower, upper) -> Interval:
    return Interval(Fraction(lower), Fraction(upper))


# --------------------------------------------------------------------------- #
# the exact interpolant gradient                                               #
# --------------------------------------------------------------------------- #


def test_cell_pressure_gradient_is_the_exact_interpolant_derivative() -> None:
    """The per-cell gradient data must reproduce the analytic derivative of
    the bilinear interpolant at sampled points, and the sub-range hulls must
    be the exact ranges (attained at endpoints), because the whole point of
    the gradient form is that no independent hull ever enters."""
    pressure = np.array([[1.0, -2.0], [0.25, 3.5]])
    dr, dz = Fraction(1, 2), Fraction(1, 4)
    cell = pressure_cell_gradient(pressure, 0, 0, dr, dz)
    # d_r at eta=0: (0.25 - 1)/dr; at eta=1: (3.5 - (-2))/dr.
    assert cell.d_r_low == Fraction(-3, 2)
    assert cell.d_r_high == Fraction(11)
    # d_z at xi=0: (-2 - 1)/dz; at xi=1: (3.5 - 0.25)/dz.
    assert cell.d_z_low == Fraction(-12)
    assert cell.d_z_high == Fraction(13)

    def bilinear_dr(xi: float, eta: float) -> float:
        return float(cell.d_r_low) * (1 - eta) + float(cell.d_r_high) * eta

    full = cell.radial_range(Fraction(0), Fraction(1))
    for eta in np.linspace(0.0, 1.0, 11):
        value = bilinear_dr(0.3, float(eta))
        assert float(full.lower) - 1e-12 <= value <= float(full.upper) + 1e-12
    # A sub-range is exactly the segment between the two endpoint values.
    sub = cell.radial_range(Fraction(1, 4), Fraction(3, 4))
    lo = float(cell.d_r_low) + float(cell.d_r_high - cell.d_r_low) * 0.25
    hi = float(cell.d_r_low) + float(cell.d_r_high - cell.d_r_low) * 0.75
    assert float(sub.lower) == pytest.approx(min(lo, hi))
    assert float(sub.upper) == pytest.approx(max(lo, hi))

    with pytest.raises(ValueError):
        cell.radial_range(Fraction(-1, 10), Fraction(1))
    with pytest.raises(ValueError):
        cell.axial_range(Fraction(1, 2), Fraction(1, 4))
    with pytest.raises(ValueError):
        pressure_cell_gradient(pressure, 0, 0, Fraction(0), dz)


# --------------------------------------------------------------------------- #
# the velocity-component cell enclosure                                        #
# --------------------------------------------------------------------------- #


def test_cell_velocity_enclosure_contains_closed_form_samples() -> None:
    """u_r, u_z, u_theta and the speed must contain a dense sample of the
    closed-form field, evaluated here from the raw exponential formula rather
    than the module's derivative closure, so a closure bug cannot poison both
    sides of the comparison."""
    family = seed_gaussian_family()
    r_box, z_box = _box(Fraction(2, 5), Fraction(1, 2)), _box(Fraction(1, 5), Fraction(3, 10))
    enclosure = cell_velocity_enclosure(family, r_box, z_box, terms=28, bits=48)
    tolerance = 1.0e-9

    def value(terms, s, z):
        total = 0.0
        for term in terms:
            poly = sum(float(c) * s**m * z**n for m, n, c in term.monomials)
            total += poly * math.exp(
                -(float(term.alpha) * s + float(term.beta) * z * z)
            )
        return total

    step = 1.0e-6
    for r in np.linspace(0.4, 0.5, 6):
        for z in np.linspace(0.2, 0.3, 6):
            s = r * r
            u_t = r * value(family.swirl, s, z)
            dz_psi = (
                value(family.stream, s, z + step)
                - value(family.stream, s, z - step)
            ) / (2 * step)
            ds_psi = (
                value(family.stream, s + step, z)
                - value(family.stream, s - step, z)
            ) / (2 * step)
            u_r = -r * dz_psi
            u_z = 2.0 * value(family.stream, s, z) + 2.0 * s * ds_psi
            speed = math.sqrt(u_r**2 + u_t**2 + u_z**2)
            fd_tol = 1.0e-6  # difference-quotient truncation dominates
            assert float(enclosure.u_theta.lower) - tolerance <= u_t
            assert u_t <= float(enclosure.u_theta.upper) + tolerance
            assert float(enclosure.u_r.lower) - fd_tol <= u_r
            assert u_r <= float(enclosure.u_r.upper) + fd_tol
            assert float(enclosure.u_z.lower) - fd_tol <= u_z
            assert u_z <= float(enclosure.u_z.upper) + fd_tol
            assert float(enclosure.speed.lower) - fd_tol <= speed
            assert speed <= float(enclosure.speed.upper) + fd_tol
    with pytest.raises(ValueError):
        cell_velocity_enclosure(family, _box(-1, 1), z_box)


# --------------------------------------------------------------------------- #
# the end-to-end certificate                                                   #
# --------------------------------------------------------------------------- #


def _payload(path: str = "plain", **kwargs):
    return build_gradient_generation_certificate(
        seed_gaussian_family(),
        _grid(),
        viscosity=5.0e-3,
        path=path,
        precision_bits=48,
        exp_terms=28,
        **kwargs,
    ).as_dict()


def test_gradient_certificate_verifies_and_reconciles() -> None:
    payload = _payload()
    report = verify_gradient_generation_certificate(payload)
    assert report["verified"] is True, report["failures"]
    assert payload["certificate"] == GRADIENT_CERTIFICATE_KIND

    rigorous = payload["rigorous"]
    assert rigorous["pressure_pairing"] == "gradient_form"
    margin = Fraction(rigorous["margin"])
    assert margin == (
        Fraction(rigorous["pressure_term_lower_bound"])
        + Fraction(rigorous["viscous_lower_bound"])
    )
    assert Fraction(rigorous["viscous_integral_upper_bound"]) == (
        Fraction(rigorous["interior_viscous_integral_upper_bound"])
        + Fraction(rigorous["exterior_viscous_tail"])
    )
    assert Fraction(rigorous["dependency_overestimation"]) == (
        Fraction(rigorous["natural_product_width"])
        - Fraction(rigorous["correlated_product_width"])
    )
    assert any("Clay" in str(claim) for claim in payload["claims"])
    # The six separated error sources are all reported as finite floats.
    diagnostics = payload["diagnostics"]
    for key in (
        "pressure_solver_residual",
        "pressure_gradient_variation_bound",
        "finite_box_viscous_tail_in_margin",
        "gaussian_velocity_tail_at_box",
        "quadrature_enclosure_width",
        "dependency_overestimation",
    ):
        assert math.isfinite(float(diagnostics[key]))


def test_gradient_enclosure_contains_a_float_midpoint_quadrature() -> None:
    """A midpoint Riemann sum of the SAME integrand (analytic |u|u paired with
    the interpolant gradient of the SAME p_h) samples the integrand inside
    each cell, so it must land inside the certified enclosure up to float
    round-off.  This is the containment test that the certificate is a
    certificate of the object it names."""
    family = seed_gaussian_family()
    grid = _grid()
    payload = build_gradient_generation_certificate(
        family, grid, viscosity=5.0e-3, precision_bits=48, exp_terms=28
    ).as_dict()
    lower = float(Fraction(payload["rigorous"]["pressure_term_lower_bound"]))
    upper = float(Fraction(payload["rigorous"]["pressure_term_upper_bound"]))

    field = gaussian_reference_field(family, grid)
    solver = AxisymmetricPressureSolver.build(grid)
    trace = solver.exterior_trace(field.u_r, field.u_theta, field.u_z)
    pressure = solver.solve(field.pressure_source(field.gradient()), trace=trace)

    # Midpoint values of |u| u at cell centres, from the analytic field.
    r_mid = 0.5 * (grid.r[:-1] + grid.r[1:])
    z_mid = 0.5 * (grid.z[:-1] + grid.z[1:])
    total = 0.0
    dr, dz = float(grid.dr), float(grid.dz)
    from ns_certificate_lab.interval_gaussian import _family_terms, _float_partials

    swirl_terms, stream_terms = _family_terms(family)
    rr, zz = np.meshgrid(r_mid, z_mid, indexing="ij")
    swirl = _float_partials(swirl_terms, rr, zz)
    stream = _float_partials(stream_terms, rr, zz)
    u_t = rr * swirl["value"]
    u_r = -rr * stream["dz"]
    u_z = 2.0 * stream["value"] + 2.0 * rr * rr * stream["ds"]
    speed = np.sqrt(u_r**2 + u_t**2 + u_z**2)
    # Interpolant gradient at cell centres: mean of the two edge quotients.
    g_r = 0.5 * (
        (pressure[1:, :-1] - pressure[:-1, :-1])
        + (pressure[1:, 1:] - pressure[:-1, 1:])
    ) / dr
    g_z = 0.5 * (
        (pressure[:-1, 1:] - pressure[:-1, :-1])
        + (pressure[1:, 1:] - pressure[1:, :-1])
    ) / dz
    total = float(
        np.sum(
            -3.0
            * 2.0
            * math.pi
            * rr
            * speed
            * (u_r * g_r + u_z * g_z)
            * dr
            * dz
        )
    )
    slack = 1.0e-12 + 1.0e-9 * abs(total)
    assert lower - slack <= total <= upper + slack


def test_taylor_and_adaptive_paths_are_never_wider_than_plain() -> None:
    plain = _payload("plain")
    taylor = _payload("taylor")
    adaptive = _payload("adaptive", max_depth=1, refine_budget=256)
    width = lambda p: Fraction(p["rigorous"]["pressure_term_width"])  # noqa: E731
    assert width(taylor) <= width(plain)
    assert width(adaptive) <= width(plain)
    assert adaptive["rigorous"]["refined_cells"] > 0
    assert adaptive["rigorous"]["leaf_evaluations"] > adaptive["rigorous"][
        "cells_enclosed"
    ]
    for payload in (plain, taylor, adaptive):
        assert verify_gradient_generation_certificate(payload)["verified"] is True
    # The adaptive sweep respects its budget: every split consumes 4 budget
    # units and nets 3 extra leaves (4 children replace 1 parent).
    spent = 256 - adaptive["rigorous"]["refine_budget_left"]
    extra = (
        adaptive["rigorous"]["leaf_evaluations"]
        - adaptive["rigorous"]["cells_enclosed"]
    )
    assert spent > 0
    assert extra == 3 * spent // 4


def test_tamper_rejections_fire() -> None:
    payload = _payload()

    forged = copy.deepcopy(payload)
    forged["rigorous"]["margin"] = str(Fraction(forged["rigorous"]["margin"]) + 1)
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["margin_is_positive"] = True
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["hypotheses"]["P1G_pressure_gradient_discretisation"]["proved"] = True
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["hypotheses"]["P1G_pressure_gradient_discretisation"][
        "statement"
    ] = "trust me"
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["claims"] = [c for c in forged["claims"] if "Clay" not in str(c)]
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["exterior_viscous_tail"] = str(
        Fraction(forged["rigorous"]["exterior_viscous_tail"]) * 2
    )
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["pressure_error_threshold"] = str(
        Fraction(forged["rigorous"]["pressure_error_threshold"]) + 1
    )
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["natural_product_width"] = str(
        Fraction(forged["rigorous"]["natural_product_width"]) + 1
    )
    assert verify_gradient_generation_certificate(forged)["verified"] is False

    assert verify_gradient_generation_certificate({"junk": True})["verified"] is False


def test_builder_rejects_bad_inputs() -> None:
    family = seed_gaussian_family()
    with pytest.raises(ValueError):
        build_gradient_generation_certificate(family, _grid(), viscosity=0.0)
    with pytest.raises(ValueError):
        build_gradient_generation_certificate(
            family, _grid(), viscosity=1e-2, path="magic"
        )
    with pytest.raises(ValueError):
        build_gradient_generation_certificate(
            family, _grid(), viscosity=1e-2, refine_budget=-1
        )
    with pytest.raises(ValueError):
        build_gradient_generation_certificate(
            family,
            AxisymmetricGrid.uniform(
                nr=9, nz=17, r_max=3.0, z_min=0.5, z_max=3.0, periodic_z=False
            ),
            viscosity=1e-2,
        )
    surrogate = GaussianFamily(
        name="surrogate",
        swirl=seed_gaussian_family().swirl,
        stream=(),
        clay_admissible=False,
    )
    with pytest.raises(ValueError):
        build_gradient_generation_certificate(surrogate, _grid(), viscosity=1e-2)


def test_compare_paths_reports_all_four() -> None:
    report = compare_certificate_paths(
        seed_gaussian_family(),
        _grid(),
        viscosity=5.0e-3,
        precision_bits=48,
        exp_terms=28,
        refine_budget=64,
        max_depth=1,
    )
    paths = report["paths"]
    assert set(paths) == {
        "divergence_hull", "gradient_plain", "gradient_adaptive", "gradient_taylor",
    }
    for name, entry in paths.items():
        assert math.isfinite(entry["margin"])
        assert entry["natural_product_width"] >= 0.0
    # The baseline margin comes from the audited divergence-form builder.
    assert paths["divergence_hull"]["payload"]["certificate"] == "l3_generation_rate"
