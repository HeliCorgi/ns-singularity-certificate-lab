"""Tests for the closed-form spectral pressure and its certificate.

The certified object is the box integral ``-3 int |u| u . grad p dx`` with
``p`` the TRUE free-space pressure of the datum, plus a bound on the exterior
remainder, so the certificate carries no pressure hypothesis at all.  Nothing
in this file asserts, and nothing in the code under test claims, anything
about the Clay problem.
"""

from __future__ import annotations

import copy
from fractions import Fraction
import math

import numpy as np
import pytest

from ns_certificate_lab.gaussian_hermite import gaussian_search_basis
from ns_certificate_lab.gaussian_spectral_pressure import (
    SPECTRAL_CERTIFICATE_KIND,
    SpectralPressure,
    _g_family,
    build_spectral_generation_certificate,
    clear_spectral_cache,
    isotropic_pressure_source,
    verify_spectral_generation_certificate,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.interval_gaussian import (
    GaussianTerm,
    GaussianFamily,
    gaussian_reference_field,
    seed_gaussian_family,
)
from ns_certificate_lab.l3_generation import AxisymmetricPressureSolver
from ns_certificate_lab.operators import derivative_r, derivative_z
from ns_certificate_lab.snapshot_certificate import Interval


#: The turn-10 search candidate: every slot nonzero, no accidental symmetry.
CANDIDATE = (-0.6197605900478841, 0.04598436893882875, 0.7202595543504078,
             -0.1228286172461205, -0.1983772283611742, -0.0845251056507508,
             0.10319690111357603, 0.15091035541730444)


def _grid(nr: int = 9, nz: int = 17) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=3.0, z_min=-3.6, z_max=3.6, periodic_z=False
    )


def _g_float(mu: float, w: float, terms: int = 600) -> float:
    """Reference g_mu by the alternating series in float (small w only)."""
    total, term = 0.0, 1.0
    for k in range(terms):
        total += term / (mu + k + 1)
        term *= -w / (k + 1)
    return total


# --------------------------------------------------------------------------- #
# the g-function enclosures                                                    #
# --------------------------------------------------------------------------- #


def test_g_family_encloses_reference_values_in_both_regimes() -> None:
    """Series regime (w < 30) and erf regime (w >= 30) must both enclose the
    float reference; the erf regime is checked against math.erf directly,
    which is an independent implementation."""
    clear_spectral_cache()
    # Small w only: the float reference series loses ~e^w * 1e-16 to
    # cancellation (the exact-rational production series does not), so beyond
    # w ~ 15 the float reference is less accurate than the enclosure.
    for w in (Fraction(0), Fraction(1, 3), Fraction(5), Fraction(12)):
        family = _g_family(w, 11, bits=53, terms=32)
        for mu2 in (-1, 1, 3, 7, 11):
            reference = _g_float(mu2 / 2.0, float(w))
            box = family[mu2]
            assert float(box.lower) - 1e-9 <= reference <= float(box.upper) + 1e-9
            assert box.upper - box.lower < Fraction(1, 10**10)
    for w in (Fraction(30), Fraction(45), Fraction(130)):
        family = _g_family(w, 11, bits=53, terms=32)
        x = math.sqrt(float(w))
        base = math.sqrt(math.pi) * math.erf(x) / x
        box = family[-1]
        assert float(box.lower) - 1e-12 <= base <= float(box.upper) + 1e-12
        # upward recursion values stay positive and below 1/(mu+1)
        for mu2 in (1, 5, 11):
            assert family[mu2].lower >= 0
            assert family[mu2].upper <= Fraction(2, mu2 + 2)


def test_g_family_is_decreasing_in_w() -> None:
    clear_spectral_cache()
    lows = [_g_family(Fraction(w), 5, bits=48, terms=28)[5] for w in (1, 2, 4)]
    assert lows[0].lower > lows[1].upper - Fraction(1, 10**9)
    assert lows[1].lower > lows[2].upper - Fraction(1, 10**9)


# --------------------------------------------------------------------------- #
# the source and the potential representation                                  #
# --------------------------------------------------------------------------- #


def test_pressure_source_matches_the_cylindrical_trace_formula() -> None:
    """sigma from the meridional closed form must equal the audited
    pressure_source of the float field at grid nodes: this pins the identity
    sigma = 2f^2 + 4sff_s + 4s^2f_s^2 - 2g^2 - 4sgg_s + 4sh_sf_z + h_z^2
    against an independently implemented formula."""
    family = gaussian_search_basis().with_amplitudes(CANDIDATE)
    source = isotropic_pressure_source(family)
    grid = _grid(11, 21)
    field = gaussian_reference_field(family, grid)
    reference = field.pressure_source(field.gradient())
    r_mesh, z_mesh = grid.mesh()
    values = np.zeros_like(reference)
    for rate, poly in source.items():
        block = np.zeros_like(reference)
        for (m, n), c in poly.items():
            block += float(c) * (r_mesh**2) ** m * z_mesh**n
        values += block * np.exp(-float(rate) * (r_mesh**2 + z_mesh**2))
    # the reference field is windowed at the outer ring/faces; compare inside
    assert np.max(np.abs(values[:-1, 1:-1] - reference[:-1, 1:-1])) < 1e-15


def test_anisotropic_generators_are_refused() -> None:
    term = GaussianTerm.from_triple({(0, 1): Fraction(1)}, Fraction(2), Fraction(3))
    family = GaussianFamily(name="aniso", swirl=(term,), stream=(term,))
    with pytest.raises(ValueError):
        isotropic_pressure_source(family)


def test_base_case_reproduces_the_erf_potential() -> None:
    """(-Delta)^{-1}[e^{-A rho^2}] = sqrt(pi) erf(sqrt(A) rho)/(4 A^{3/2} rho):
    the classical Newtonian potential of a Gaussian, evaluated through the
    g-monomial machinery on a synthetic one-term source."""
    from ns_certificate_lab.gaussian_spectral_pressure import _potential_terms

    rate = Fraction(5, 2)
    rep = _potential_terms({(0, 0): Fraction(1)}, rate)
    assert set(rep) == {(0, 0, -1)}
    assert rep[(0, 0, -1)] == Fraction(1, 4) / rate
    for rho in (0.35, 1.2, 2.4):
        w = float(rate) * rho * rho
        value = float(rep[(0, 0, -1)]) * _g_float(-0.5, w)
        reference = (
            math.sqrt(math.pi) * math.erf(math.sqrt(float(rate)) * rho)
            / (4.0 * float(rate) ** 1.5 * rho)
        )
        assert value == pytest.approx(reference, rel=1e-12)


def test_laplace_identity_holds_exactly_for_the_candidate() -> None:
    """The exact rational reduction of Delta p + sigma to zero is the
    self-verification of the whole closed form: decomposition, evolution,
    t-integration and derivative closure all stand or fall with it."""
    family = gaussian_search_basis().with_amplitudes(CANDIDATE)
    spectral = SpectralPressure.from_family(family)
    report = spectral.verify_laplace_identity()
    assert report["verified"] is True, report["failures"]


def test_laplace_identity_detects_a_corrupted_representation() -> None:
    family = seed_gaussian_family()
    spectral = SpectralPressure.from_family(family)
    rate = next(iter(spectral.p))
    corrupted = {a: dict(rep) for a, rep in spectral.p.items()}
    key = next(iter(corrupted[rate]))
    corrupted[rate][key] += Fraction(1, 7)
    from ns_certificate_lab.gaussian_spectral_pressure import (
        _rep_d_s, _rep_d_z,
    )
    lattice = {(0, 0): corrupted}
    for order in range(1, 3):
        for a in range(order + 1):
            b = order - a
            if a > 0:
                lattice[(a, b)] = {
                    rr: _rep_d_s(rep, rr) for rr, rep in lattice[(a - 1, b)].items()
                }
            else:
                lattice[(a, b)] = {
                    rr: _rep_d_z(rep, rr) for rr, rep in lattice[(a, b - 1)].items()
                }
    broken = SpectralPressure(source=spectral.source, derivatives=lattice)
    report = broken.verify_laplace_identity()
    assert report["verified"] is False


def test_spectral_gradient_matches_the_discrete_solver() -> None:
    """The spectral gradient is exact; the discrete solve converges to it, so
    at a moderate grid they must agree to discretisation accuracy.  This wires
    the closed form against the audited numerical pipeline end to end."""
    family = gaussian_search_basis().with_amplitudes(CANDIDATE)
    spectral = SpectralPressure.from_family(family)
    grid = _grid(21, 41)
    field = gaussian_reference_field(family, grid)
    solver = AxisymmetricPressureSolver.build(grid)
    trace = solver.exterior_trace(field.u_r, field.u_theta, field.u_z)
    pressure = solver.solve(field.pressure_source(field.gradient()), trace=trace)
    p_r = derivative_r(grid, pressure, even_at_axis=True)
    p_z = derivative_z(grid, pressure)
    scale = float(np.max(np.abs(p_r))) + float(np.max(np.abs(p_z)))
    for (i, j) in ((5, 14), (8, 20), (11, 26), (3, 24)):
        r0 = Fraction(float(grid.r[i]))
        z0 = Fraction(float(grid.z[j]))
        d_r, d_z = spectral.gradient_value(r0, z0, bits=48, terms=28)
        mid_r = float(d_r.lower + d_r.upper) / 2
        mid_z = float(d_z.lower + d_z.upper) / 2
        assert abs(mid_r - p_r[i, j]) < 3.0e-2 * scale
        assert abs(mid_z - p_z[i, j]) < 3.0e-2 * scale


def test_gradient_boxes_contain_thin_samples_and_taylor_is_never_wider() -> None:
    family = gaussian_search_basis().with_amplitudes(CANDIDATE)
    spectral = SpectralPressure.from_family(family)
    r_box = Interval(Fraction(3, 8), Fraction(9, 20))
    z_box = Interval(Fraction(27, 100), Fraction(9, 25))
    plain_s, plain_z = spectral.gradient_boxes(
        r_box, z_box, path="plain", bits=48, terms=28
    )
    taylor_s, taylor_z = spectral.gradient_boxes(
        r_box, z_box, path="taylor", taylor_order=4, bits=48, terms=28
    )
    assert taylor_s.upper - taylor_s.lower <= plain_s.upper - plain_s.lower
    assert taylor_z.upper - taylor_z.lower <= plain_z.upper - plain_z.lower
    for r in np.linspace(float(r_box.lower), float(r_box.upper), 5):
        for z in np.linspace(float(z_box.lower), float(z_box.upper), 5):
            s_thin, z_thin = spectral.gradient_boxes(
                Interval(Fraction(r), Fraction(r)),
                Interval(Fraction(z), Fraction(z)),
                path="plain", bits=48, terms=28,
            )
            # d_s p at the sample must lie inside both cell enclosures
            assert taylor_s.lower <= s_thin.upper
            assert s_thin.lower <= taylor_s.upper
            assert taylor_z.lower <= z_thin.upper
            assert z_thin.lower <= taylor_z.upper


def test_exterior_gradient_bound_dominates_far_samples() -> None:
    family = gaussian_search_basis().with_amplitudes(CANDIDATE)
    spectral = SpectralPressure.from_family(family)
    bound = spectral.exterior_gradient_bound(Fraction(3), bits=48)
    assert bound > 0
    for (r, z) in ((3.0, 0.5), (2.5, 2.0), (0.5, 3.2), (4.0, 1.0)):
        if r * r + z * z < 9.0:
            continue
        d_r, d_z = spectral.gradient_value(Fraction(r), Fraction(z))
        magnitude = max(abs(float(d_r.lower)), abs(float(d_r.upper))) + max(
            abs(float(d_z.lower)), abs(float(d_z.upper))
        )
        assert magnitude <= float(bound) * (1.0 + 1e-9)


# --------------------------------------------------------------------------- #
# the certificate                                                              #
# --------------------------------------------------------------------------- #


def _payload(path: str = "taylor", **kwargs):
    return build_spectral_generation_certificate(
        seed_gaussian_family(),
        _grid(),
        viscosity=5.0e-3,
        path=path,
        taylor_order=3,
        precision_bits=48,
        exp_terms=28,
        **kwargs,
    ).as_dict()


def test_spectral_certificate_verifies_and_carries_no_hypotheses() -> None:
    payload = _payload()
    report = verify_spectral_generation_certificate(payload)
    assert report["verified"] is True, report["failures"]
    assert payload["certificate"] == SPECTRAL_CERTIFICATE_KIND
    assert payload["hypotheses"] == {}
    rigorous = payload["rigorous"]
    margin = Fraction(rigorous["margin"])
    assert margin == (
        Fraction(rigorous["pressure_term_lower_bound"])
        - Fraction(rigorous["exterior_pressure_remainder"])
        + Fraction(rigorous["viscous_lower_bound"])
    )
    assert Fraction(rigorous["exterior_pressure_remainder"]) >= 0
    assert any("Clay" in str(claim) for claim in payload["claims"])
    assert any("unconditional" in str(claim) for claim in payload["claims"])
    # the float reference of the SAME box integral falls inside the enclosure
    assert payload["diagnostics"]["float_pressure_in_enclosure"] is True


def test_adaptive_path_is_never_wider_and_respects_its_budget() -> None:
    taylor = _payload("taylor")
    adaptive = _payload("adaptive", max_depth=1, refine_budget=64)
    width = lambda p: Fraction(p["rigorous"]["pressure_term_width"])  # noqa: E731
    assert width(adaptive) <= width(taylor)
    rigorous = adaptive["rigorous"]
    assert rigorous["refined_cells"] > 0
    spent = 64 - rigorous["refine_budget_left"]
    extra = rigorous["leaf_evaluations"] - rigorous["cells_enclosed"]
    assert spent > 0
    assert extra == 3 * spent // 4
    assert verify_spectral_generation_certificate(adaptive)["verified"] is True


def test_tamper_rejections_fire() -> None:
    payload = _payload()

    forged = copy.deepcopy(payload)
    forged["rigorous"]["margin"] = str(Fraction(forged["rigorous"]["margin"]) + 1)
    assert verify_spectral_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["margin_is_positive"] = True
    assert verify_spectral_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["hypotheses"] = {"P1G": {"proved": True}}
    assert verify_spectral_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["claims"] = [c for c in forged["claims"] if "Clay" not in str(c)]
    assert verify_spectral_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["exterior_pressure_remainder"] = str(Fraction(-1, 10**9))
    assert verify_spectral_generation_certificate(forged)["verified"] is False

    forged = copy.deepcopy(payload)
    forged["rigorous"]["laplace_identity_checked"] = False
    assert verify_spectral_generation_certificate(forged)["verified"] is False

    assert verify_spectral_generation_certificate({"junk": True})["verified"] is False


def test_builder_rejects_bad_inputs() -> None:
    family = seed_gaussian_family()
    with pytest.raises(ValueError):
        build_spectral_generation_certificate(family, _grid(), viscosity=0.0)
    with pytest.raises(ValueError):
        build_spectral_generation_certificate(
            family, _grid(), viscosity=1e-2, path="magic"
        )
    with pytest.raises(ValueError):
        build_spectral_generation_certificate(
            family, _grid(), viscosity=1e-2, taylor_order=9
        )
    with pytest.raises(ValueError):
        build_spectral_generation_certificate(
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
        build_spectral_generation_certificate(surrogate, _grid(), viscosity=1e-2)
