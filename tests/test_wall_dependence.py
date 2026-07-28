"""Acceptance tests for the preregistered wall-dependence experiment.

The real ``T_1`` family is deliberately never integrated here.  These tests
exercise

* the E-32 initial-data family against the numbers audited in
  ``docs/equation_audit.md`` E-32 (core bit-identity, compact support, the
  E-32c deviation bound, the ``C^4`` spot check and the axis parity check),
* the configuration contract of ``experiments/run_wall_dependence.py``,
* the artifact bundle, on a grid and a time window small enough for CI,
* explicit fault injections proving each core-identity gate can fail.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from experiments import run_wall_dependence
from experiments.run_wall_dependence import (
    CORE_SNAPSHOT_FIELDS,
    EXPERIMENT_ID,
    IMPLEMENTATION_DECISION_KEYS,
    PREREGISTERED_DECISION_KEYS,
    PREREG_SECTION_6_VERBATIM,
    amplification_separations,
    argmax_displacement,
    build_members,
    build_member_grid,
    classify_wall_dependence,
    core_identity_report,
    core_reference_report,
    elliptic_nonlocal_pair,
    evaluate,
    radial_point_count,
    run,
    validate_config,
)
from ns_certificate_lab._integrity import strict_json_loads
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    DIAGNOSTIC_FIELDS,
    RELATIVE_DIAGNOSTIC_FIELDS,
    hou_initial_swirl,
    load_checkpoint,
)
from ns_certificate_lab.wall_dependence import (
    E32_BAND_FOURTH_DIFFERENCE,
    E32_CORE_RADIUS,
    E32_RHO1,
    E32_RHO2,
    E32_SUPPORT_RADIUS,
    E32_SUP_DEVIATION_BOUND,
    E32_TRANSITION_BAND,
    e29_radial_profile,
    envelope_initial_swirl,
    envelope_radial_profile,
    initial_data_acceptance,
    smooth_cutoff,
    smooth_transition,
    transition_band_fourth_difference,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "wall_dependence.json"

T1 = 0.002191729

# ``docs/equation_audit.md`` E-32 property 4: the measured sup deviation on the
# production 193-point radial grid.
E32_MEASURED_SUP_DEVIATION_NR193 = 2.578e-12


def _config() -> dict[str, Any]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="wall dependence test config",
    )
    assert isinstance(value, dict)
    return value


def _smoke_config() -> dict[str, Any]:
    """Two tiny member groups, ~20 steps each: the real T_1 window never runs."""

    config = _config()
    config["member_groups"] = [
        {
            "label": "core_dr32",
            "points_per_unit_radius": 32,
            "nz": 64,
            "wall_radii": [1.0, 1.5],
            "primary": True,
        },
        {
            "label": "core_dr16",
            "points_per_unit_radius": 16,
            "nz": 64,
            "wall_radii": [1.0, 1.5],
            "primary": False,
        },
    ]
    config["t_final"] = 2.0e-5
    config["snapshot_times"] = [0.0, 1.0e-5, 2.0e-5]
    config["max_steps"] = 200
    config["diagnostic_stride"] = 1
    # The two acceptance numbers that are resolution dependent.  Both are
    # implementation gates, not preregistered thresholds: the E-29b constants
    # are continuum maxima that a 33-point radial grid cannot resolve, and the
    # axis parity defect of a clean even field is O(dr^3).
    config["acceptance"] = dict(config["acceptance"])
    config["acceptance"]["maximum_initial_norm_relative_error"] = 0.05
    config["acceptance"]["maximum_initial_axis_parity_relative"] = 0.1
    return config


def _unit_grid(nr: int, nz: int) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=nr, nz=nz, r_max=1.0, z_min=0.0, z_max=1.0, periodic_z=True
    )


@pytest.fixture(scope="module")
def generated_run(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    config = _smoke_config()
    output = tmp_path_factory.mktemp("wall-dependence") / "evidence"
    summary = run(config, output)
    return config, output, summary


# ============================================================ E-32 properties


def test_smooth_transition_is_the_standard_flat_function() -> None:
    assert float(smooth_transition(0.0)) == 0.0
    assert float(smooth_transition(-1.0)) == 0.0
    assert float(smooth_transition(1.0)) == pytest.approx(math.exp(-1.0))
    assert float(smooth_transition(0.5)) == pytest.approx(math.exp(-2.0))
    values = np.asarray(smooth_transition([-1.0, 0.0, 0.25, 1.0]))
    assert values.shape == (4,)
    assert values[0] == 0.0 and values[1] == 0.0
    assert np.all(np.diff(values) >= 0.0)
    with pytest.raises(ValueError, match="finite arguments"):
        smooth_transition([0.0, float("nan")])


def test_smooth_cutoff_branches_are_exact_and_monotone() -> None:
    """E-32b: the plateau is exactly 1, the tail exactly 0, chi in [0,1]."""

    assert smooth_cutoff(0.0) == 1.0
    assert smooth_cutoff(E32_RHO1) == 1.0
    assert smooth_cutoff(E32_RHO1 - 1e-15) == 1.0
    assert smooth_cutoff(E32_RHO2) == 0.0
    assert smooth_cutoff(1.0) == 0.0
    assert smooth_cutoff(9.0) == 0.0
    # The construction is symmetric about the band midpoint.
    assert smooth_cutoff(0.5 * (E32_RHO1 + E32_RHO2)) == pytest.approx(0.5, abs=1e-12)

    rho = np.linspace(0.0, 1.0, 4001)
    chi = np.asarray(smooth_cutoff(rho))
    assert chi.shape == rho.shape
    assert np.all(chi >= 0.0) and np.all(chi <= 1.0)
    assert np.all(np.diff(chi) <= 0.0)
    assert np.all(chi[rho <= E32_RHO1] == 1.0)
    assert np.all(chi[rho >= E32_RHO2] == 0.0)


def test_smooth_cutoff_rejects_bad_parameters() -> None:
    with pytest.raises(ValueError, match="rho1 must be strictly below rho2"):
        smooth_cutoff(0.5, 0.9, 0.9)
    with pytest.raises(ValueError, match="rho1 must be strictly below rho2"):
        smooth_cutoff(0.5, 0.95, 0.81)
    with pytest.raises(ValueError, match="rho1 must be nonnegative"):
        smooth_cutoff(0.5, -0.1, 0.9)
    with pytest.raises(ValueError, match="rho1 must be finite"):
        smooth_cutoff(0.5, float("nan"), 0.9)
    with pytest.raises(ValueError, match="rho2 must be finite"):
        smooth_cutoff(0.5, 0.1, float("inf"))
    with pytest.raises(ValueError, match="finite arguments"):
        smooth_cutoff([0.1, float("nan")])


def test_envelope_is_bit_identical_to_e29_on_the_core() -> None:
    """E-32 property 3: multiplying by the literal 1.0 is exact."""

    for nr, nz, r_max in ((193, 384, 1.0), (289, 384, 1.5), (577, 384, 3.0)):
        grid = AxisymmetricGrid.uniform(
            nr=nr, nz=nz, r_max=r_max, z_min=0.0, z_max=1.0, periodic_z=True
        )
        envelope = envelope_initial_swirl(grid)
        plain = hou_initial_swirl(grid, amplitude=12000.0)
        core = grid.r <= E32_CORE_RADIUS
        assert np.array_equal(envelope[core], plain[core])
        # E-32 property 5: exactly zero outside the support, strictly inside
        # every wall.
        outside = grid.r >= E32_SUPPORT_RADIUS
        assert np.count_nonzero(outside) > 0
        assert np.all(envelope[outside] == 0.0)
        assert float(grid.r[-1]) == r_max


def test_envelope_core_block_is_shared_across_every_wall_radius() -> None:
    """The whole design rests on this: one core datum, four wall radii."""

    grids = [
        AxisymmetricGrid.uniform(
            nr=192 * int(2 * radius) // 2 + 1,
            nz=384,
            r_max=radius,
            z_min=0.0,
            z_max=1.0,
            periodic_z=True,
        )
        for radius in (1.0, 1.5, 2.0, 3.0)
    ]
    assert [item.nr for item in grids] == [193, 289, 385, 577]
    rows = int(np.count_nonzero(grids[0].r <= E32_CORE_RADIUS))
    reference = envelope_initial_swirl(grids[0])[:rows]
    for grid in grids[1:]:
        assert np.array_equal(grid.r[:rows], grids[0].r[:rows])
        assert np.array_equal(envelope_initial_swirl(grid)[:rows], reference)


def test_envelope_rejects_bad_parameters() -> None:
    grid = _unit_grid(65, 64)
    with pytest.raises(ValueError, match="rho1 must be strictly below rho2"):
        envelope_initial_swirl(grid, rho1=0.95, rho2=0.81)
    with pytest.raises(ValueError, match="rho1 must be strictly below rho2"):
        envelope_initial_swirl(grid, rho1=0.81, rho2=0.81)
    with pytest.raises(ValueError, match="rho2 must not exceed r_max"):
        envelope_initial_swirl(grid, rho1=0.81, rho2=1.5)
    with pytest.raises(ValueError, match="rho1 must be finite"):
        envelope_initial_swirl(grid, rho1=float("nan"))
    with pytest.raises(ValueError, match="amplitude_scale must be finite"):
        envelope_initial_swirl(grid, amplitude_scale=float("inf"))
    with pytest.raises(ValueError, match="amplitude_scale must be positive"):
        envelope_initial_swirl(grid, amplitude_scale=0.0)
    with pytest.raises(TypeError, match="AxisymmetricGrid"):
        envelope_initial_swirl(object())  # type: ignore[arg-type]


def test_e32_acceptance_reproduces_the_audited_numbers_at_nr193() -> None:
    """Every number here is quoted from ``docs/equation_audit.md`` E-32."""

    report = initial_data_acceptance(_unit_grid(193, 384))
    assert report["core_bit_identical"] is True
    assert report["exact_zero_outside_support"] is True
    assert report["cutoff_minimum"] == 0.0
    assert report["cutoff_maximum"] == 1.0
    assert report["cutoff_monotone_non_increasing"] is True

    # E-32c: the analytic bound and the measured grid value.
    assert report["sup_deviation_analytic_bound"] == pytest.approx(
        3.4008e-10, rel=1e-4
    )
    assert E32_SUP_DEVIATION_BOUND <= 3.4008e-10
    assert report["sup_deviation"] <= E32_SUP_DEVIATION_BOUND
    assert report["sup_deviation"] <= 5.0e-10
    assert report["sup_deviation"] == pytest.approx(
        E32_MEASURED_SUP_DEVIATION_NR193, rel=1e-3
    )
    assert report["sup_deviation_within_analytic_bound"] is True
    assert report["sup_deviation_within_acceptance"] is True

    # E-32 property 4: the derived E-29b norms are unchanged by the envelope.
    norms = report["derived_norms"]
    assert norms["max_abs_u1_relative_change"] == 0.0
    assert norms["max_cartesian_vorticity_relative_change"] == 0.0
    assert norms["plain_max_abs_u1"] == norms["envelope_max_abs_u1"]

    # E-32 property 7: the C^4 spot check across the transition band.
    band = report["fourth_difference"]
    assert band["identical"] is True
    assert band["bounded"] is True
    assert band["relative_difference"] == 0.0
    assert band["plain"]["maximum"] == pytest.approx(
        E32_BAND_FOURTH_DIFFERENCE[193], rel=1e-4
    )
    assert band["envelope"]["maximum"] == pytest.approx(21.739, rel=1e-4)
    assert band["plain"]["band_lower"] == E32_TRANSITION_BAND[0]

    # Preregistration (vi): the axis parity check.
    parity = report["axis_parity"]
    assert parity["relative_change"] == 0.0
    assert parity["envelope_relative"] < 1.0e-3
    assert parity["envelope_defect"] == parity["plain_defect"]


def test_e32_acceptance_reproduces_the_audited_numbers_at_nr385() -> None:
    report = initial_data_acceptance(_unit_grid(385, 768))
    assert report["core_bit_identical"] is True
    assert report["exact_zero_outside_support"] is True
    assert report["sup_deviation"] <= E32_SUP_DEVIATION_BOUND
    assert report["sup_deviation"] == pytest.approx(2.578e-12, rel=1e-3)
    band = report["fourth_difference"]
    assert band["identical"] is True
    assert band["plain"]["maximum"] == pytest.approx(
        E32_BAND_FOURTH_DIFFERENCE[385], rel=1e-4
    )
    assert band["envelope"]["maximum"] == pytest.approx(26.362, rel=1e-4)
    assert report["derived_norms"]["max_abs_u1_relative_change"] == 0.0
    assert report["axis_parity"]["relative_change"] == 0.0


def test_e32_acceptance_refuses_a_non_unit_grid() -> None:
    """E-29 is undefined past r = 1, so the comparison is refused there."""

    wide = AxisymmetricGrid.uniform(
        nr=385, nz=384, r_max=2.0, z_min=0.0, z_max=1.0, periodic_z=True
    )
    with pytest.raises(ValueError, match="use r_max = 1"):
        initial_data_acceptance(wide)


def test_transition_band_fourth_difference_contract() -> None:
    radii = np.linspace(0.0, 1.0, 193)
    plain = e29_radial_profile(radii, 0.0845842)
    envelope = envelope_radial_profile(radii, 0.0845842)
    report = transition_band_fourth_difference(radii, plain)
    assert report["maximum"] == pytest.approx(21.739, rel=1e-4)
    assert report["argmax_r"] == pytest.approx(0.8541666666666666)
    assert report["sample_count"] == 27
    assert transition_band_fourth_difference(radii, envelope)["maximum"] == (
        report["maximum"]
    )
    with pytest.raises(ValueError, match="matching 1-D samples"):
        transition_band_fourth_difference(radii, plain[:-1])
    with pytest.raises(ValueError, match="band must be nonempty"):
        transition_band_fourth_difference(radii, plain, band=(1.0, 0.85))


# ======================================================== config validation


def test_shipped_config_matches_the_preregistered_design() -> None:
    config = _config()
    validate_config(config)
    assert config["experiment_id"] == EXPERIMENT_ID
    assert config["z_period"] == 1.0
    assert config["core_radius"] == 0.9
    assert config["tail_radius"] == 1.1
    assert config["envelope"] == {"rho1": E32_RHO1, "rho2": E32_RHO2}
    assert config["amplitude_scale"] == 1.0
    assert config["viscosity_schedule"] == [[0.0, 5e-4]]
    assert config["t_final"] == pytest.approx(T1)
    # E-30 stage one only: T_1 sits below the switch time t_0.
    assert config["t_final"] < 0.00227375
    assert config["snapshot_times"] == [0.0, 0.0005, 0.001, 0.0015, T1]
    assert config["cfl_coefficient"] == 0.1
    assert config["max_time_step"] == 1e-6

    groups = config["member_groups"]
    assert [item["label"] for item in groups] == ["core_dr192", "core_dr128"]
    assert groups[0]["wall_radii"] == [1.0, 1.5, 2.0, 3.0]
    assert groups[0]["points_per_unit_radius"] == 192
    assert groups[0]["nz"] == 384
    assert groups[0]["primary"] is True
    assert groups[1]["wall_radii"] == [1.0, 2.0]
    assert groups[1]["points_per_unit_radius"] == 128
    assert groups[1]["primary"] is False

    # The preregistered nested radii give exactly the required nr values.
    members = build_members(config)
    assert [item.nr for item in members if item.group == "core_dr192"] == [
        193,
        289,
        385,
        577,
    ]
    assert [item.nr for item in members if item.group == "core_dr128"] == [129, 257]
    assert all(item.nz == 384 for item in members)

    # dr = 1/192 and dz = 1/384 exactly, and identical across the whole group.
    grids = {
        item.group: [] for item in members
    }
    for item in members:
        grids[item.group].append(build_member_grid(config, item))
    assert len(set(grid.dr for grid in grids["core_dr192"])) == 1
    assert grids["core_dr192"][0].dr == 1.0 / 192.0
    assert grids["core_dr192"][0].dz == 1.0 / 384.0
    assert len(set(grid.dr for grid in grids["core_dr128"])) == 1
    assert grids["core_dr128"][0].dr == 1.0 / 128.0

    decision = config["decision_rule"]
    assert decision["wall_dependent_amplification_ratio"] == 0.20
    assert decision["wall_dependent_argmax_relative_displacement"] == 0.10
    assert decision["wall_effect_small_amplification_ratio"] == 0.05
    assert set(decision) == set(PREREGISTERED_DECISION_KEYS) | set(
        IMPLEMENTATION_DECISION_KEYS
    )


def test_radial_point_count_refuses_a_non_integer_radius() -> None:
    assert radial_point_count(192, 1.0) == 193
    assert radial_point_count(192, 1.5) == 289
    assert radial_point_count(192, 2.0) == 385
    assert radial_point_count(192, 3.0) == 577
    assert radial_point_count(128, 1.0) == 129
    assert radial_point_count(128, 2.0) == 257
    # 2.5 IS on the 1/192 lattice (192 * 2.5 = 480), so it is accepted; the
    # rejected radii are the ones that would land r = R_wall off the lattice.
    assert radial_point_count(192, 2.5) == 481
    for radius in (1.3, 1.7, 2.4, 1.0 + 1.0e-3):
        with pytest.raises(ValueError, match="does not yield an integer"):
            radial_point_count(192, radius)
    with pytest.raises(ValueError, match="wall_radius must be positive"):
        radial_point_count(192, 0.0)


def test_detects_a_bad_radius_in_the_config() -> None:
    """A radius off the shared lattice must be refused, never rounded."""

    config = _config()
    config["member_groups"] = [dict(item) for item in config["member_groups"]]
    config["member_groups"][0] = dict(config["member_groups"][0])
    config["member_groups"][0]["wall_radii"] = [1.0, 1.5, 2.0, 2.5001]
    with pytest.raises(ValueError, match="does not yield an integer"):
        validate_config(config)


def test_detects_mislabeled_config() -> None:
    config = _config()
    config["experiment_id"] = "hou_early_time_v1"
    with pytest.raises(ValueError, match="experiment_id must be wall_dependence_v1"):
        validate_config(config)
    config["experiment_id"] = ""
    with pytest.raises(ValueError, match="experiment_id must be"):
        validate_config(config)
    config["experiment_id"] = EXPERIMENT_ID
    validate_config(config)


def test_missing_and_unknown_config_keys_are_rejected() -> None:
    for key in (
        "member_groups",
        "core_radius",
        "tail_radius",
        "envelope",
        "decision_rule",
        "acceptance",
        "snapshot_times",
        "units",
    ):
        partial = _config()
        del partial[key]
        with pytest.raises(ValueError, match="missing or unknown keys"):
            validate_config(partial)

    unknown = _config()
    unknown["r_max"] = 1.0
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(unknown)

    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(["schema_version", 1])  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"schema_version": 2}, "schema_version must be 1"),
        ({"core_radius": 0.0}, "core_radius must be positive"),
        ({"tail_radius": 0.5}, "tail_radius must exceed core_radius"),
        ({"envelope": {"rho1": 0.81}}, "exactly rho1 and rho2"),
        ({"envelope": {"rho1": 0.81, "rho2": 0.81}}, "strictly below"),
        ({"envelope": {"rho1": 0.5, "rho2": 0.9025}}, "at least core_radius"),
        ({"envelope": {"rho1": 0.81, "rho2": 1.21}}, "smallest wall radius"),
        ({"amplitude_scale": 0.0}, "amplitude_scale must be positive"),
        ({"cfl_coefficient": -1.0}, "cfl_coefficient must be positive"),
        ({"max_time_step": 0.0}, "max_time_step must be positive"),
        ({"max_steps": 0}, "max_steps must be a positive integer"),
        ({"diagnostic_stride": 1.5}, "diagnostic_stride must be a positive"),
        ({"t_final": 0.0}, "t_final must be positive"),
        ({"z_period": -1.0}, "z_period must be positive"),
        ({"seed": "abc"}, "seed must be an integer"),
        ({"units": {}}, "units must be a nonempty object"),
        ({"description": "  "}, "description must be a nonempty string"),
        ({"snapshot_times": []}, "snapshot_times must be a nonempty list"),
        ({"snapshot_times": [0.0, 0.001]}, "last snapshot time must be t_final"),
        ({"viscosity_schedule": [[0.001, 5e-4]]}, "first viscosity switch time"),
        ({"cfl_excess_tolerance": -0.01}, "cfl_excess_tolerance must be nonneg"),
    ],
)
def test_invalid_config_values_are_rejected(
    mutation: dict[str, Any],
    message: str,
) -> None:
    config = _config()
    config.update(mutation)
    with pytest.raises(ValueError, match=message):
        validate_config(config)


@pytest.mark.parametrize(
    ("groups", "message"),
    [
        ([], "member_groups must be a nonempty list"),
        (
            [
                {
                    "label": "a",
                    "points_per_unit_radius": 192,
                    "nz": 384,
                    "wall_radii": [1.0],
                    "primary": True,
                }
            ],
            "at least two wall radii",
        ),
        (
            [
                {
                    "label": "a",
                    "points_per_unit_radius": 192,
                    "nz": 384,
                    "wall_radii": [2.0, 1.0],
                    "primary": True,
                }
            ],
            "strictly increasing",
        ),
        (
            [
                {
                    "label": "a",
                    "points_per_unit_radius": 192,
                    "nz": 384,
                    "wall_radii": [1.0, 2.0],
                    "primary": False,
                }
            ],
            "exactly one member group must be marked primary",
        ),
        (
            [
                {
                    "label": "a",
                    "points_per_unit_radius": 192,
                    "nz": 384,
                    "wall_radii": [1.0, 2.0],
                    "primary": True,
                },
                {
                    "label": "a",
                    "points_per_unit_radius": 128,
                    "nz": 384,
                    "wall_radii": [1.0, 2.0],
                    "primary": False,
                },
            ],
            "labels must be unique",
        ),
        (
            [
                {
                    "label": "a",
                    "points_per_unit_radius": 192,
                    "nz": 4,
                    "wall_radii": [1.0, 2.0],
                    "primary": True,
                }
            ],
            "nz is below the AxisymmetricGrid minimum",
        ),
        (
            [
                {
                    "label": "a",
                    "points_per_unit_radius": 192,
                    "wall_radii": [1.0, 2.0],
                    "primary": True,
                }
            ],
            "member group has missing or unknown keys",
        ),
    ],
)
def test_invalid_member_groups_are_rejected(
    groups: list[dict[str, Any]],
    message: str,
) -> None:
    config = _config()
    config["member_groups"] = groups
    with pytest.raises(ValueError, match=message):
        validate_config(config)


def test_invalid_acceptance_and_decision_blocks_are_rejected() -> None:
    extra = _config()
    extra["acceptance"] = dict(extra["acceptance"])
    extra["acceptance"]["bonus"] = 1.0
    with pytest.raises(ValueError, match="acceptance object has missing or unknown"):
        validate_config(extra)

    missing = _config()
    missing["acceptance"] = dict(missing["acceptance"])
    del missing["acceptance"]["maximum_core_reference_relative_spread"]
    with pytest.raises(ValueError, match="acceptance object has missing or unknown"):
        validate_config(missing)

    decision = _config()
    decision["decision_rule"] = dict(decision["decision_rule"])
    del decision["decision_rule"]["wall_effect_small_amplification_ratio"]
    with pytest.raises(ValueError, match="decision_rule object has missing or unknown"):
        validate_config(decision)

    inverted = _config()
    inverted["decision_rule"] = dict(inverted["decision_rule"])
    inverted["decision_rule"]["wall_effect_small_amplification_ratio"] = 0.5
    with pytest.raises(ValueError, match="must be below"):
        validate_config(inverted)


# ================================================== decision-rule arithmetic


def test_amplification_separations_and_classification_follow_prereg_section_4() -> None:
    radii = [1.0, 1.5, 2.0, 3.0]

    # A separation of 0.30 on the largest pair: wall dependent.
    dependent = amplification_separations([13.0, 11.0, 10.0, 13.0], radii)
    assert dependent[-1]["separation"] == pytest.approx(3.0 / 13.0)
    verdict = classify_wall_dependence(
        dependent,
        {"relative_displacement": 0.0},
        thresholds={
            "wall_dependent_amplification_ratio": 0.2,
            "wall_dependent_argmax_relative_displacement": 0.1,
            "wall_effect_small_amplification_ratio": 0.05,
        },
        resolution_stable=True,
    )
    assert verdict["classification"] == "wall_dependent"

    # S = (0.10, 0.05, 0.02): strictly decreasing and inside 0.05 on the
    # largest pair, so wall effect small.
    small = amplification_separations([11.781, 10.71, 10.2, 10.0], radii)
    assert [item["separation"] for item in small] == pytest.approx(
        [0.10, 0.05, 0.02], rel=1e-6
    )
    verdict = classify_wall_dependence(
        small,
        {"relative_displacement": 0.0},
        thresholds={
            "wall_dependent_amplification_ratio": 0.2,
            "wall_dependent_argmax_relative_displacement": 0.1,
            "wall_effect_small_amplification_ratio": 0.05,
        },
        resolution_stable=True,
    )
    assert verdict["classification"] == "wall_effect_small"
    assert verdict["separations_strictly_decreasing"] is True
    assert verdict["largest_pair_separation"] <= 0.05

    # In between: undecided.
    middle = amplification_separations([12.0, 11.5, 11.0, 10.0], radii)
    verdict = classify_wall_dependence(
        middle,
        {"relative_displacement": 0.0},
        thresholds={
            "wall_dependent_amplification_ratio": 0.2,
            "wall_dependent_argmax_relative_displacement": 0.1,
            "wall_effect_small_amplification_ratio": 0.05,
        },
        resolution_stable=True,
    )
    assert verdict["classification"] == "undecided"


def test_a_large_argmax_shift_alone_classifies_as_wall_dependent() -> None:
    separations = amplification_separations([10.0, 10.0, 10.0, 10.0], [1.0, 1.5, 2.0, 3.0])
    displacement = argmax_displacement((0.3, 0.10), (0.0, 0.08), z_period=1.0)
    assert displacement["relative_displacement"] > 0.1
    verdict = classify_wall_dependence(
        separations,
        displacement,
        thresholds={
            "wall_dependent_amplification_ratio": 0.2,
            "wall_dependent_argmax_relative_displacement": 0.1,
            "wall_effect_small_amplification_ratio": 0.05,
        },
        resolution_stable=True,
    )
    assert verdict["classification"] == "wall_dependent"
    # Preregistration section 2: an unstable measurement holds the verdict.
    held = classify_wall_dependence(
        separations,
        displacement,
        thresholds={
            "wall_dependent_amplification_ratio": 0.2,
            "wall_dependent_argmax_relative_displacement": 0.1,
            "wall_effect_small_amplification_ratio": 0.05,
        },
        resolution_stable=False,
    )
    assert held["classification"] == "undecided"
    assert held["classification_before_qualifier"] == "wall_dependent"
    assert held["classification_with_section_2_hold"] == "undecided"

    # Section 4 attaches the consistency clause to the wall_dependent branch
    # only; section 2 states the more general hold. Both readings are recorded.
    small = amplification_separations(
        [11.781, 10.71, 10.2, 10.0], [1.0, 1.5, 2.0, 3.0]
    )
    unstable = classify_wall_dependence(
        small,
        {"relative_displacement": 0.0},
        thresholds={
            "wall_dependent_amplification_ratio": 0.2,
            "wall_dependent_argmax_relative_displacement": 0.1,
            "wall_effect_small_amplification_ratio": 0.05,
        },
        resolution_stable=False,
    )
    assert unstable["classification"] == "wall_effect_small"
    assert unstable["classification_with_section_2_hold"] == "undecided"


def test_argmax_displacement_measures_z_on_the_circle() -> None:
    report = argmax_displacement((0.0, 0.99), (0.0, 0.01), z_period=1.0)
    assert report["axial_separation"] == pytest.approx(0.02)
    assert report["displacement"] == pytest.approx(0.02)
    assert report["denominator"] == pytest.approx(0.01)


# ============================================================ smoke bundle


def test_smoke_run_produces_a_complete_accepted_bundle(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, _, summary = generated_run
    assert summary["experiment_id"] == EXPERIMENT_ID
    expected_checks = {
        "all_members_completed",
        "all_diagnostics_finite",
        "energy_non_increasing",
        "circulation_max_principle",
        "odd_symmetry_preserved",
        "advective_cfl_within_tolerance",
        "core_discretization_identical",
        "initial_core_reference_shared",
        "initial_data_family_accepted",
        "initial_norms_match_e29b",
        "elliptic_nonlocal_contribution_recorded",
        "cross_solver_elliptic_agreement_recorded",
    }
    assert set(summary["acceptance_checks"]) == expected_checks
    assert all(summary["acceptance_checks"].values())
    assert summary["accepted_as_wall_dependence_observation"] is True
    assert "not a proof" in summary["interpretation"]
    assert summary["preregistration"] == "docs/wall_dependence_prereg.md"
    assert summary["reproducibility"]["seed"] == config["seed"]
    assert summary["reproducibility"]["runtime_provenance"][
        "source_fingerprint_sha256"
    ]
    assert summary["known_gaps"]

    # Preregistration section 6, verbatim plus the English renderings.
    assert summary["limitations_preregistration_section_6_verbatim"] == list(
        PREREG_SECTION_6_VERBATIM
    )
    assert len(PREREG_SECTION_6_VERBATIM) == 4
    limitations = " ".join(summary["limitations"])
    assert "aspect ratio" in limitations
    assert "wall to infinity" in limitations
    assert "artificial" in limitations
    assert "not a reproduction" in limitations

    members = summary["members"]
    assert len(members) == 4
    for item in members:
        assert item["completed"] is True
        assert item["failure"] is None
        assert item["all_diagnostics_finite"] is True
        assert item["final_time"] == pytest.approx(config["t_final"])
        assert item["maximum_energy_growth_ratio"] <= 0.0
        assert item["maximum_wall_u1_abs"] == 0.0
        assert item["snapshot_count"] == len(config["snapshot_times"])


def test_core_identity_is_an_acceptance_check_and_holds(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    _, _, summary = generated_run
    assert summary["acceptance_checks"]["core_discretization_identical"] is True
    assert summary["acceptance_checks"]["initial_core_reference_shared"] is True
    for block in summary["analysis"]["member_groups"].values():
        identity = block["core_identity"]
        assert identity["radial_spacing_identical"] is True
        assert identity["radial_spacing_absolute_spread"] == 0.0
        assert identity["core_radial_nodes_bitwise_identical"] is True
        assert identity["axial_grid_identical"] is True
        assert identity["core_rows"] >= 2

        reference = block["core_reference"]
        assert reference["bitwise_identical"] is True
        assert reference["relative_spread"] == 0.0
        assert reference["core_equals_full_for_every_member"] is True
        assert reference["initial_core_energy_relative_spread"] == 0.0
        assert reference["reference"] > 0.0


def test_every_preregistered_core_metric_is_present_with_its_denominator(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, summary = generated_run
    with (output / "core_snapshots.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    assert list(rows[0]) == list(CORE_SNAPSHOT_FIELDS)
    # 4 members x 3 snapshots.
    assert len(rows) == 4 * len(config["snapshot_times"])

    for row in rows:
        # 1. amplification, always with the denominator it was taken against
        assert math.isfinite(float(row["core_amplification"]))
        assert float(row["core_amplification_denominator"]) > 0.0
        # 2. argmax location of u1
        assert 0.0 <= float(row["core_argmax_u1_r"]) <= config["core_radius"]
        assert 0.0 <= float(row["core_argmax_u1_z"]) < config["z_period"]
        # 3. meridional flow strength
        assert float(row["core_max_abs_psi1_z"]) >= 0.0
        # 4. core energy
        assert float(row["core_energy"]) > 0.0
        # 6. tail amplitudes, only where the preregistration asks for them
        if float(row["wall_radius"]) > 1.0:
            assert float(row["tail_max_abs_u1"]) >= 0.0
            assert float(row["tail_max_abs_omega1"]) >= 0.0
            assert int(row["tail_row_count"]) > 0
        else:
            assert row["tail_max_abs_u1"] == ""
        # the solver-B cross-check with its own denominator
        assert math.isfinite(float(row["psi_cross_solver_max_abs_difference"]))
        assert math.isfinite(float(row["psi_cross_solver_relative_denominator"]))

    # The initial snapshot is the shared datum: identical core amplification
    # for every member of a group.
    initial = {
        (row["group"], float(row["core_amplification"]))
        for row in rows
        if row["snapshot_index"] == "0"
    }
    assert len({group for group, _ in initial}) == 2
    assert len(initial) == 2

    definitions = summary["analysis"]["denominator_definitions"]
    assert "core_amplification" in definitions
    assert "elliptic_nonlocal_relative_difference" in definitions


def test_elliptic_nonlocal_contribution_is_recorded_for_every_pair(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, summary = generated_run
    entries = summary["analysis"]["elliptic_nonlocal_contribution"]
    # one pair per group, three snapshots each
    assert len(entries) == 2 * len(config["snapshot_times"])
    assert "no interpolation" in summary["analysis"]["elliptic_nonlocal_definition"]
    for item in entries:
        assert item["status"] == "measured"
        assert item["smaller_wall_radius"] < item["larger_wall_radius"]
        assert math.isfinite(item["max_abs_difference"])
        assert item["max_abs_difference"] >= 0.0
        assert math.isfinite(item["relative_denominator"])
        assert item["relative_denominator"] >= 0.0
        assert item["large_core_max_abs_psi1"] == item["relative_denominator"]
    # omega1(0) = 0, so the first snapshot's elliptic difference is exactly 0.
    first = [item for item in entries if item["snapshot_index"] == 0]
    assert first
    assert all(item["max_abs_difference"] == 0.0 for item in first)
    # Later snapshots carry a genuine nonzero nonlocal contribution.
    later = [item for item in entries if item["snapshot_index"] > 0]
    assert any(item["max_abs_difference"] > 0.0 for item in later)

    with (output / "elliptic_nonlocal.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(entries)
    assert {"relative_difference", "relative_denominator"} <= rows[0].keys()


def test_classification_is_recorded_but_never_gates(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    _, _, summary = generated_run
    classification = summary["wall_dependence_classification"]
    assert classification in {"wall_dependent", "wall_effect_small", "undecided"}
    decision = summary["wall_dependence_decision"]
    assert decision["classification"] == classification
    assert decision["classification_before_qualifier"] in decision[
        "allowed_classifications"
    ]
    assert decision["classification_with_section_2_hold"] in decision[
        "allowed_classifications"
    ]
    assert "section 2" in decision["qualifier_readings"]
    assert decision["gated"] is False
    assert "not an acceptance gate" in decision["interpretation"]
    assert decision["thresholds"]["wall_dependent_amplification_ratio"] == 0.2
    assert decision["thresholds"]["wall_effect_small_amplification_ratio"] == 0.05
    assert decision["thresholds"][
        "wall_dependent_argmax_relative_displacement"
    ] == 0.1
    assert decision["reasons"]
    assert isinstance(decision["separations"], list)
    assert len(decision["separations"]) == 1
    # The classification is absent from the acceptance gate set.
    assert "wall_dependence_classification" not in summary["acceptance_checks"]
    assert "classification" not in " ".join(summary["acceptance_checks"])

    consistency = summary["analysis"]["resolution_consistency"]
    assert consistency["tolerance_source"].startswith("not fixed by")
    assert consistency["comparisons"]
    assert isinstance(consistency["stable_under_refinement"], bool)


def test_smoke_run_artifacts_are_finite_and_checksummed(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, _ = generated_run
    with (output / "diagnostics.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    assert {"group", "wall_radius", "nr", "nz"} <= rows[0].keys()
    assert set(DIAGNOSTIC_FIELDS) <= rows[0].keys()
    assert set(RELATIVE_DIAGNOSTIC_FIELDS) <= rows[0].keys()
    # ``group`` is the only non-numeric column.
    assert all(
        math.isfinite(float(row[name]))
        for row in rows
        for name in row
        if name != "group"
    )
    assert {row["group"] for row in rows} == {"core_dr32", "core_dr16"}

    with (output / "members.csv").open(newline="", encoding="utf-8") as stream:
        member_rows = list(csv.DictReader(stream))
    assert len(member_rows) == 4
    assert {
        "wall_radius",
        "dr",
        "core_rows",
        "final_core_amplification",
        "final_core_argmax_u1_r",
        "final_core_max_abs_psi1_z",
        "final_core_energy",
        "reference_vorticity",
    } <= member_rows[0].keys()

    with np.load(output / "trajectories.npz", allow_pickle=False) as archive:
        names = set(archive.files)
        assert "r_core_dr32_R1" in names
        assert "psi1_final_core_dr32_R1p5" in names
        assert "core_amplification_core_dr32_R1" in names
        assert all(np.all(np.isfinite(archive[name])) for name in names)

    manifest = strict_json_loads(
        (output / "manifest.json").read_text(encoding="utf-8"),
        label="wall dependence manifest",
    )
    assert manifest["experiment_id"] == EXPERIMENT_ID
    assert {
        "config.snapshot.json",
        "diagnostics.csv",
        "members.csv",
        "core_snapshots.csv",
        "elliptic_nonlocal.csv",
        "summary.json",
        "trajectories.npz",
    } <= set(manifest["files"])
    assert any(name.startswith("checkpoints/") for name in manifest["files"])
    for name, record in manifest["files"].items():
        payload = (output / name).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == record["sha256"]
        assert len(payload) == record["bytes"]
    assert (output / "manifest.json.sha256").read_text(
        encoding="ascii"
    ).strip() == hashlib.sha256((output / "manifest.json").read_bytes()).hexdigest()

    snapshot = strict_json_loads(
        (output / "config.snapshot.json").read_text(encoding="utf-8"),
        label="wall dependence config snapshot",
    )
    assert snapshot == config
    assert json.loads((output / "summary.json").read_text(encoding="utf-8"))


def test_snapshot_checkpoints_reload_as_v2_candidates(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, _ = generated_run
    archives = sorted((output / "checkpoints").glob("*.npz"))
    assert len(archives) == 4 * len(config["snapshot_times"])
    wide = [path for path in archives if "R1p5" in path.name]
    artifact = load_candidate(wide[-1])
    assert set(artifact.fields) == {"u1", "omega1", "psi1"}
    assert artifact.schema.endswith("/v2")

    checkpoint = load_checkpoint(wide[-1])
    assert np.array_equal(checkpoint.state.psi1, artifact.fields["psi1"])
    assert checkpoint.metadata["experiment_id"] == EXPERIMENT_ID
    assert checkpoint.metadata["wall_radius"] == 1.5
    assert float(checkpoint.grid.r[-1]) == 1.5
    # E-27 at the outer row of the WIDER cylinder, with no solver change.
    assert float(np.max(np.abs(checkpoint.state.u1[-1]))) == 0.0
    assert float(np.max(np.abs(checkpoint.state.psi1[-1]))) == 0.0


def test_nonempty_output_directory_is_refused(scratch_dir: Path) -> None:
    config = _smoke_config()
    occupied = scratch_dir / "occupied"
    occupied.mkdir()
    marker = occupied / "preserve.txt"
    marker.write_text("keep\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(config, occupied)
    assert marker.read_text(encoding="utf-8") == "keep\n"
    assert list(occupied.iterdir()) == [marker]


# ========================================================== fault injections


def test_detects_a_member_whose_core_spacing_differs() -> None:
    """(i) A member off the shared lattice must fail the core-identity gate."""

    config = _smoke_config()
    members = build_members(config)
    good = [build_member_grid(config, item) for item in members[:2]]
    rows = int(np.count_nonzero(good[0].r <= config["core_radius"]))
    assert core_identity_report(good, members[:2], core_rows=rows)[
        "radial_spacing_identical"
    ]

    # Same nr, larger radius: dr is wrong and the core nodes no longer line up.
    broken = AxisymmetricGrid.uniform(
        nr=good[0].nr, nz=good[0].nz, r_max=1.5, z_min=0.0, z_max=1.0, periodic_z=True
    )
    report = core_identity_report(
        [good[0], broken], members[:2], core_rows=rows
    )
    assert report["radial_spacing_identical"] is False
    assert report["radial_spacing_absolute_spread"] > 0.0
    assert report["core_radial_nodes_bitwise_identical"] is False

    # The elliptic nonlocal comparison must refuse such a pair outright.
    with pytest.raises(ValueError, match="do not share the core nodes bitwise"):
        elliptic_nonlocal_pair(
            np.zeros((rows, good[0].nz)),
            small_grid=good[0],
            large_grid=broken,
            core_rows=rows,
        )


def test_detects_a_wrong_core_spacing_end_to_end(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The same fault, injected through the real pipeline, flips the gate."""

    config = _smoke_config()
    config["member_groups"] = [config["member_groups"][0]]
    original = run_wall_dependence.build_member_grid

    def corrupted(cfg: dict[str, Any], spec: Any) -> AxisymmetricGrid:
        if spec.wall_radius == 1.5:
            # One radial point too few: dr no longer matches the group.
            return AxisymmetricGrid.uniform(
                nr=spec.nr - 1,
                nz=spec.nz,
                r_max=spec.wall_radius,
                z_min=0.0,
                z_max=float(cfg["z_period"]),
                periodic_z=True,
            )
        return original(cfg, spec)

    monkeypatch.setattr(run_wall_dependence, "build_member_grid", corrupted)
    _, _, report, checks = evaluate(config)
    assert checks["core_discretization_identical"] is False
    assert checks["elliptic_nonlocal_contribution_recorded"] is False
    entries = report["elliptic_nonlocal_contribution"]
    assert entries
    assert any(item["status"] == "core_nodes_not_shared" for item in entries)


def test_detects_a_perturbed_envelope() -> None:
    """(ii) Moving rho1 to 0.5 breaks core identity and the C^4 spot check."""

    grid = _unit_grid(193, 384)
    plain = hou_initial_swirl(grid, amplitude=12000.0)
    perturbed = envelope_initial_swirl(grid, rho1=0.5, rho2=E32_RHO2)
    core = grid.r <= E32_CORE_RADIUS
    assert not np.array_equal(perturbed[core], plain[core])

    report = initial_data_acceptance(grid, rho1=0.5, rho2=E32_RHO2)
    assert report["core_bit_identical"] is False
    assert report["sup_deviation"] > E32_SUP_DEVIATION_BOUND
    assert report["sup_deviation_within_analytic_bound"] is False
    assert report["sup_deviation_within_acceptance"] is False
    assert report["fourth_difference"]["identical"] is False
    assert report["fourth_difference"]["bounded"] is False
    assert report["fourth_difference"]["relative_difference"] > 1.0e-3
    assert report["derived_norms"]["max_abs_u1_relative_change"] == 0.0

    # The config contract refuses the same perturbation before any run starts.
    config = _config()
    config["envelope"] = {"rho1": 0.5, "rho2": E32_RHO2}
    with pytest.raises(ValueError, match="at least core_radius"):
        validate_config(config)


def test_detects_a_corrupted_initial_core_reference() -> None:
    """(iii) A member whose initial core datum differs breaks the shared norm."""

    config = _smoke_config()
    members = build_members(config)[:2]
    good = {
        "max_cartesian_vorticity_core": 7317.881455790404,
        "max_cartesian_vorticity_full": 7317.881455790404,
        "core_energy": 4008.19,
    }
    clean = core_reference_report([dict(good), dict(good)], members)
    assert clean["relative_spread"] == 0.0
    assert clean["bitwise_identical"] is True

    corrupt = dict(good)
    corrupt["max_cartesian_vorticity_core"] *= 1.0 + 1.0e-9
    corrupt["max_cartesian_vorticity_full"] = corrupt["max_cartesian_vorticity_core"]
    broken = core_reference_report([dict(good), corrupt], members)
    assert broken["bitwise_identical"] is False
    assert broken["relative_spread"] > 1.0e-12


def test_detects_a_corrupted_initial_core_reference_end_to_end(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The same fault, injected through the real pipeline, flips the gate."""

    config = _smoke_config()
    config["member_groups"] = [config["member_groups"][0]]
    original = run_wall_dependence.envelope_initial_swirl

    def corrupted(grid: AxisymmetricGrid, **options: Any) -> np.ndarray:
        values = original(grid, **options)
        if float(grid.r[-1]) == 1.5:
            # A relative perturbation far above the 1e-12 acceptance bound,
            # far below anything a diagnostic plot would show.
            return values * (1.0 + 1.0e-9)
        return values

    monkeypatch.setattr(run_wall_dependence, "envelope_initial_swirl", corrupted)
    _, _, report, checks = evaluate(config)
    assert checks["initial_core_reference_shared"] is False
    reference = report["member_groups"]["core_dr32"]["core_reference"]
    assert reference["bitwise_identical"] is False
    assert reference["relative_spread"] > 1.0e-12
    # The core-identity gate on the GRIDS still passes: this fault is about the
    # datum, and the two gates are genuinely independent.
    assert checks["core_discretization_identical"] is True
