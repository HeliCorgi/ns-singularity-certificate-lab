"""Tests for the Track-P chain: composition, recurrence, checker, tampering."""

from __future__ import annotations

import copy
import json
from fractions import Fraction

import pytest

from ns_certificate_lab.fourier_torus import (
    TORUS_FAMILIES,
    taylor_coefficients,
)
from ns_certificate_lab.torus_chain import (
    ALLOWED_WORDING,
    FORBIDDEN_WORDING,
    STOP_IS_NOT_A_SINGULARITY,
    build_chain_certificate,
    curl,
    galerkin_taylor_endpoint,
    l2_inner,
    prove_chain_slab,
    recentre,
    verify_chain_certificate,
)
from ns_certificate_lab.torus_aposteriori import prove_galerkin_box

VISCOSITY = Fraction(1, 10)
CUTOFF = 2  # P1 lives in 1 <= |k|^2 <= 2; the small band keeps tests fast
STEP = Fraction(1, 2048)


@pytest.fixture(scope="module")
def chain_payload() -> dict:
    return build_chain_certificate(
        "P1",
        viscosity=VISCOSITY,
        cutoff_sq=CUTOFF,
        initial_step=STEP,
        max_slabs=2,
        allow_step_doubling=False,
    )


@pytest.fixture(scope="module")
def chain_verdict(chain_payload) -> dict:
    return verify_chain_certificate(chain_payload)


def test_chain_closes_and_verifies(chain_payload, chain_verdict):
    assert chain_payload["slab_count"] == 2
    assert chain_payload["termination"]["reason"] == "slab_budget_exhausted"
    assert chain_verdict["failures"] == []
    assert chain_verdict["verified"] and chain_verdict["proves_existence"]


def test_delta_recurrence_is_exact_and_contiguous(chain_payload):
    rows = chain_payload["delta_recurrence"]
    running = Fraction(0)
    clock = Fraction(0)
    for row in rows:
        assert Fraction(row["delta_in"]) == running
        assert Fraction(row["t_start"]) == clock
        delta_out = Fraction(row["delta_end"]) + Fraction(row["transfer"])
        assert Fraction(row["delta_out"]) == delta_out
        running = delta_out
        clock += Fraction(row["step"])
    assert Fraction(chain_payload["certified_final_time"]) == clock
    assert Fraction(chain_payload["final_delta_h4_dot"]) == running


def test_transfer_stays_at_the_remainder_floor(chain_payload):
    for entry in chain_payload["slabs"]:
        transfer = Fraction(entry["recentre"]["transfer_h4"])
        assert transfer < Fraction(1, 10**6)
    assert chain_payload["wrapping_diagnostics"]["fired"] is False


def test_chain_grows_monotonically(chain_payload):
    rows = chain_payload["delta_recurrence"]
    deltas = [Fraction(row["delta_out"]) for row in rows]
    assert all(b > a for a, b in zip(deltas, deltas[1:]))


def test_builder_output_contains_no_forbidden_wording(chain_payload):
    text = json.dumps(chain_payload)
    for phrase in FORBIDDEN_WORDING:
        assert phrase not in text
    for sentence in ALLOWED_WORDING:
        assert sentence in chain_payload["claims"]


def test_checker_rejects_tampered_delta(chain_payload):
    forged = copy.deepcopy(chain_payload)
    forged["slabs"][-1]["delta_out"] = str(
        Fraction(forged["slabs"][-1]["delta_out"]) / 2
    )
    assert not verify_chain_certificate(forged)["verified"]


def test_checker_rejects_forged_restart_point(chain_payload):
    forged = copy.deepcopy(chain_payload)
    modes = forged["slabs"][0]["next_datum_modes"]
    key = next(iter(modes))
    modes[key]["cos"][0] = str(Fraction(modes[key]["cos"][0]) + Fraction(1, 7))
    assert not verify_chain_certificate(forged)["verified"]


def test_checker_rejects_forbidden_wording(chain_payload):
    forged = copy.deepcopy(chain_payload)
    forged["claims"] = list(forged["claims"]) + [
        "this unconditionally proves the existence of the solution"
    ]
    assert not verify_chain_certificate(forged)["verified"]


def test_checker_rejects_missing_disclaimer(chain_payload):
    forged = copy.deepcopy(chain_payload)
    forged["claims"] = [
        claim for claim in forged["claims"] if claim != STOP_IS_NOT_A_SINGULARITY
    ]
    assert not verify_chain_certificate(forged)["verified"]


def test_checker_rejects_forged_external_theorem(chain_payload):
    forged = copy.deepcopy(chain_payload)
    forged["external_theorems"]["EXT-P1"]["proved"] = True
    assert not verify_chain_certificate(forged)["verified"]


def test_checker_rejects_unregistered_termination(chain_payload):
    forged = copy.deepcopy(chain_payload)
    forged["termination"] = {"reason": "made_up_reason"}
    assert not verify_chain_certificate(forged)["verified"]


def test_checker_rejects_garbage():
    assert not verify_chain_certificate({"certificate": "nonsense"})["verified"]


# --------------------------------------------------------------------------- #
# pieces                                                                       #
# --------------------------------------------------------------------------- #


def test_slab_accepts_nonzero_initial_radius():
    datum = TORUS_FAMILIES["P1"]()
    slab = prove_chain_slab(
        datum,
        delta_in=Fraction(1, 1000),
        viscosity=VISCOSITY,
        cutoff_sq=CUTOFF,
        step=STEP,
    )
    assert slab.proved
    assert slab.radius_slab > Fraction(1, 1000)
    assert slab.delta_out > slab.delta_in


def test_taylor_endpoint_tightens_with_order():
    datum = TORUS_FAMILIES["P1"]()
    box = prove_galerkin_box(
        datum, viscosity=VISCOSITY, cutoff_sq=CUTOFF, step=STEP, max_halvings=0
    )
    assert box.proved
    coarse = galerkin_taylor_endpoint(
        datum, box.box, viscosity=VISCOSITY, cutoff_sq=CUTOFF, step=STEP, order=4
    )
    fine = galerkin_taylor_endpoint(
        datum, box.box, viscosity=VISCOSITY, cutoff_sq=CUTOFF, step=STEP, order=6
    )
    # the h4dot bound floors at the sobolev_sq rounding slop, so the norm
    # diagnostic can only be asked to not grow; the coefficient widths show
    # the genuine tightening
    assert fine.remainder_h4 <= coarse.remainder_h4
    assert fine.max_width < coarse.max_width

    exact = taylor_coefficients(
        datum, viscosity=VISCOSITY, cutoff_sq=CUTOFF, order=8
    )
    reference = {}
    for i, coefficient in enumerate(exact):
        for comp, scalar in enumerate(coefficient.components):
            for k, (a, b) in scalar.terms.items():
                key = (comp, k)
                cos_v, sin_v = reference.get(key, (Fraction(0), Fraction(0)))
                reference[key] = (a * STEP**i + cos_v, b * STEP**i + sin_v)
    for comp, scalar in enumerate(fine.field.components):
        for k, (a, b) in scalar.terms.items():
            cos_ref, sin_ref = reference.get((comp, k), (Fraction(0), Fraction(0)))
            slack = Fraction(1, 10**12)
            assert a.lower - slack <= cos_ref <= a.upper + slack
            assert b.lower - slack <= sin_ref <= b.upper + slack


def test_recentre_is_exactly_divergence_free():
    datum = TORUS_FAMILIES["P1"]()
    box = prove_galerkin_box(
        datum, viscosity=VISCOSITY, cutoff_sq=CUTOFF, step=STEP, max_halvings=0
    )
    endpoint = galerkin_taylor_endpoint(
        datum, box.box, viscosity=VISCOSITY, cutoff_sq=CUTOFF, step=STEP, order=6
    )
    restart = recentre(endpoint, cutoff_sq=CUTOFF)
    assert not restart.datum.divergence().cleaned().terms
    assert restart.transfer < Fraction(1, 10**6)


def test_l2_inner_matches_norm():
    datum = TORUS_FAMILIES["P1"]()
    assert l2_inner(datum, datum) == datum.sobolev_sq(0)


def test_curl_is_divergence_free():
    datum = TORUS_FAMILIES["P3"]()
    vorticity = curl(datum)
    assert not vorticity.divergence().cleaned().terms


def test_physical_snapshot_ranges_contain_centre(chain_payload):
    for entry in chain_payload["slabs"]:
        ranges = entry["physical"]["certified_ranges_for_true_solution"]
        for name in ("h4_dot_norm", "h3_dot_norm", "energy", "enstrophy"):
            block = ranges[name]
            centre = Fraction(block["centre"])
            low, high = (Fraction(v) for v in block["certified_range"])
            assert low <= centre <= high


# --------------------------------------------------------------------------- #
# the n = 3 Kato mode                                                          #
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def kato_payload() -> dict:
    from ns_certificate_lab.kato_constant import build_kato_certificate

    payload = build_kato_certificate(lattice_cut=20)
    return payload.as_dict() if hasattr(payload, "as_dict") else payload


@pytest.fixture(scope="module")
def h3_chain_payload(kato_payload) -> dict:
    return build_chain_certificate(
        "P1",
        viscosity=VISCOSITY,
        cutoff_sq=CUTOFF,
        initial_step=Fraction(1, 256),
        max_slabs=2,
        allow_step_doubling=False,
        mode="h3_kato",
        kato_payload=kato_payload,
    )


def test_h3_chain_closes_and_verifies(h3_chain_payload):
    assert h3_chain_payload["slab_count"] == 2
    assert h3_chain_payload["mode"] == "h3_kato"
    assert h3_chain_payload["norm_order"] == 3
    verdict = verify_chain_certificate(h3_chain_payload)
    assert verdict["failures"] == []
    assert verdict["verified"] and verdict["proves_existence"]


def test_h3_linear_coefficient_beats_h4(chain_payload, h3_chain_payload):
    old = Fraction(chain_payload["slabs"][0]["constants"]["linear_coefficient"])
    new = Fraction(h3_chain_payload["slabs"][0]["constants"]["linear_coefficient"])
    assert new < old / 4  # the crude 9(K1+K2) must lose by a wide margin


def test_h3_constants_block_names_the_inequality(h3_chain_payload):
    constants = h3_chain_payload["slabs"][0]["constants"]
    linear = (
        -VISCOSITY + Fraction(constants["C_kato"]) + Fraction(constants["C_shift"])
    )
    assert Fraction(constants["linear_coefficient"]) == linear
    assert Fraction(constants["quadratic_coefficient"]) == Fraction(
        constants["G3_upper"]
    )


def test_h3_checker_rejects_forged_g3(h3_chain_payload):
    forged = copy.deepcopy(h3_chain_payload)
    forged["kato_certificate"]["g3"]["upper"] = str(
        Fraction(forged["kato_certificate"]["g3"]["upper"]) / 2
    )
    assert not verify_chain_certificate(forged)["verified"]


def test_h3_mode_requires_kato_certificate():
    with pytest.raises(ValueError):
        build_chain_certificate(
            "P1", viscosity=VISCOSITY, cutoff_sq=CUTOFF, max_slabs=1,
            mode="h3_kato", kato_payload=None,
        )


def test_h3_physical_snapshot_has_no_h4_range(h3_chain_payload):
    ranges = h3_chain_payload["slabs"][0]["physical"][
        "certified_ranges_for_true_solution"
    ]
    assert "h4_dot_norm" not in ranges
    assert "h3_dot_norm" in ranges


# --------------------------------------------------------------------------- #
# the turn-11 reissue                                                          #
# --------------------------------------------------------------------------- #


def test_reissue_verifies_and_flips_wording(chain_payload):
    from ns_certificate_lab.torus_chain import AUDITED_KIND, reissue_chain_certificate

    reissued = reissue_chain_certificate(chain_payload)
    verdict = verify_chain_certificate(reissued)
    assert verdict["failures"] == []
    assert verdict["verified"]
    assert reissued["conclusion"]["kind"] == AUDITED_KIND
    assert AUDITED_KIND in reissued["claims"]
    assert ALLOWED_WORDING[0] not in reissued["claims"]
    for name, block in reissued["external_theorems"].items():
        assert block["proved"] is True
        assert block["closure"]["lean_formalised"] is False
        assert block["closure"]["axiomatised_in_lean"] is False


def test_reissue_keeps_quantitative_content(chain_payload):
    from ns_certificate_lab.torus_chain import reissue_chain_certificate

    reissued = reissue_chain_certificate(chain_payload)
    assert reissued["delta_recurrence"] == chain_payload["delta_recurrence"]
    assert reissued["slab_count"] == chain_payload["slab_count"]
    assert (
        reissued["certified_final_time"] == chain_payload["certified_final_time"]
    )


def test_reissued_h3_names_p3_3(h3_chain_payload):
    from ns_certificate_lab.torus_chain import reissue_chain_certificate

    reissued = reissue_chain_certificate(h3_chain_payload)
    assert "P3-3" in reissued["conclusion"]["continuation_form_consumed"]
    verdict = verify_chain_certificate(reissued)
    assert verdict["failures"] == []


def test_checker_rejects_reissue_mixtures(chain_payload):
    from ns_certificate_lab.torus_chain import reissue_chain_certificate

    reissued = reissue_chain_certificate(chain_payload)

    forged = copy.deepcopy(reissued)
    forged["external_theorems"]["EXT-P1"]["proved"] = False
    assert not verify_chain_certificate(forged)["verified"]

    forged = copy.deepcopy(reissued)
    del forged["external_theorems"]["EXT-P2"]["closure"]
    assert not verify_chain_certificate(forged)["verified"]

    forged = copy.deepcopy(reissued)
    forged["external_theorems"]["EXT-P3"]["closure"]["lean_formalised"] = True
    assert not verify_chain_certificate(forged)["verified"]

    forged = copy.deepcopy(reissued)
    forged["external_theorems"]["EXT-P2"]["dini_clause"]["closed"] = True
    assert not verify_chain_certificate(forged)["verified"]

    forged = copy.deepcopy(reissued)
    forged["claims"] = list(forged["claims"]) + [ALLOWED_WORDING[0]]
    assert not verify_chain_certificate(forged)["verified"]

    forged = copy.deepcopy(reissued)
    forged["conclusion"]["kind"] = ALLOWED_WORDING[0]
    assert not verify_chain_certificate(forged)["verified"]
