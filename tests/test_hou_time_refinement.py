"""Acceptance tests for the Hou fixed-grid time-step refinement study.

The real ``T_1`` configuration is deliberately never integrated here: these
tests exercise the configuration contract, the step-planning arithmetic, the
step-doubling analysis and the artifact bundle on a grid and a time window
small enough for continuous integration.
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

from experiments.run_hou_time_refinement import (
    ARGMAX_QUANTITIES,
    CONSTRAINT_QUANTITIES,
    EXPECTED_TIME_STEP_COUNT,
    EXPERIMENT_ID,
    SMOOTH_QUANTITIES,
    SPATIAL_AMPLIFICATION_LADDER,
    SPATIAL_LADDER_SUMMARIES,
    cfl_within_tolerance,
    evolve_time_step,
    ladder_halves_exactly,
    members_share_final_time,
    plan_ladder,
    plan_step,
    read_spatial_amplification_ladder,
    run,
    step_doubling_order,
    temporal_order_analysis,
    validate_config,
    validate_time_steps,
)
from ns_certificate_lab._integrity import strict_json_loads
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.nonlinear_cylinder import (
    DIAGNOSTIC_FIELDS,
    RELATIVE_DIAGNOSTIC_FIELDS,
    load_checkpoint,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "hou_time_refinement.json"

T1 = 0.002191729


def _config() -> dict[str, Any]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="hou time refinement test config",
    )
    assert isinstance(value, dict)
    return value


def _smoke_config() -> dict[str, Any]:
    """~10/20/40 steps on a tiny grid: the real T_1 window is never run here."""

    config = _config()
    config["resolution"] = [17, 32]
    config["t_final"] = 2.0e-5
    config["time_steps"] = [2.0e-6, 1.0e-6, 5.0e-7]
    config["max_steps"] = 200
    config["diagnostic_stride"] = 1
    return config


@pytest.fixture(scope="module")
def generated_run(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    config = _smoke_config()
    output = tmp_path_factory.mktemp("hou-time-refinement") / "evidence"
    summary = run(config, output)
    return config, output, summary


# ------------------------------------------------------- config validation


def test_shipped_config_is_valid_and_targets_the_audited_early_time() -> None:
    config = _config()
    validate_config(config)
    assert config["experiment_id"] == EXPERIMENT_ID
    assert config["r_max"] == 1.0
    assert config["z_period"] == 1.0
    assert config["amplitude_scale"] == 1.0
    assert config["resolution"] == [65, 128]
    # A single viscosity stage: T_1 lies strictly below the E-30 switch time
    # t_0 = 0.00227375, so no switch can occur inside this window.
    assert config["viscosity_schedule"] == [[0.0, 5e-4]]
    assert config["t_final"] == pytest.approx(T1)
    assert config["t_final"] < 0.00227375
    assert config["time_steps"] == [6e-7, 3e-7, 1.5e-7]
    assert config["cfl_coefficient"] == 0.1
    assert config["acceptance"]["maximum_circulation_growth_ratio"] == 1e-3
    assert config["acceptance"]["temporal_order_lower_bound"] == 1.6
    assert config["acceptance"]["temporal_order_upper_bound"] == 2.4


def test_detects_mislabeled_config() -> None:
    """A config carrying another experiment's id must never be accepted."""

    config = _config()
    config["experiment_id"] = "hou_early_time_v1"
    with pytest.raises(
        ValueError, match="experiment_id must be hou_time_refinement_v1"
    ):
        validate_config(config)

    config["experiment_id"] = ""
    with pytest.raises(ValueError, match="experiment_id must be"):
        validate_config(config)

    config["experiment_id"] = EXPERIMENT_ID
    validate_config(config)


def test_detects_non_halving_time_steps() -> None:
    """The ladder must be exactly ``dt, dt/2, dt/4``."""

    for ladder in (
        [6e-7, 3.5e-7, 1.5e-7],
        [6e-7, 1.5e-7, 7.5e-8],
        [6e-7, 3e-7, 1.0e-7],
    ):
        with pytest.raises(ValueError, match="time_steps must halve exactly"):
            validate_time_steps(ladder)
        config = _config()
        config["time_steps"] = ladder
        with pytest.raises(ValueError, match="time_steps must halve exactly"):
            validate_config(config)

    for ladder in ([6e-7, 3e-7, 3e-7], [6e-7, 6e-7, 3e-7], [3e-7, 6e-7, 1.2e-6]):
        with pytest.raises(ValueError, match="strictly decreasing"):
            validate_time_steps(ladder)

    for ladder in ([6e-7, 3e-7], [6e-7, 3e-7, 1.5e-7, 7.5e-8], [], 6e-7):
        with pytest.raises(ValueError, match="time_steps must be a list of 3"):
            validate_time_steps(ladder)

    for ladder in ([6e-7, 0.0, 1.5e-7], [6e-7, -3e-7, 1.5e-7], [6e-7, "3e-7", 1.5e-7]):
        with pytest.raises(ValueError, match="time_steps entry must be"):
            validate_time_steps(ladder)

    # The shipped ladder is accepted and returned unchanged.
    assert validate_time_steps([6e-7, 3e-7, 1.5e-7]) == (6e-7, 3e-7, 1.5e-7)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"schema_version": 2}, "schema_version must be 1"),
        ({"resolution": [65, 128, 256]}, "single \\[nr, nz\\] integer pair"),
        ({"resolution": [65]}, "single \\[nr, nz\\] integer pair"),
        ({"resolution": [65.0, 128]}, "single \\[nr, nz\\] integer pair"),
        ({"resolution": [3, 128]}, "below AxisymmetricGrid minima"),
        ({"resolution": [65, 4]}, "below AxisymmetricGrid minima"),
        ({"amplitude_scale": 0.0}, "amplitude_scale must be positive"),
        ({"cfl_coefficient": -1.0}, "cfl_coefficient must be positive"),
        ({"cfl_excess_tolerance": -0.01}, "cfl_excess_tolerance must be nonnegative"),
        ({"cfl_excess_tolerance": "loose"}, "cfl_excess_tolerance must be a finite"),
        ({"max_steps": 0}, "max_steps must be a positive integer"),
        ({"diagnostic_stride": 1.5}, "diagnostic_stride must be a positive"),
        ({"t_final": 0.0}, "t_final must be positive"),
        (
            {
                "viscosity_schedule": [[0.0, 5e-4], [0.00227375, 5e-3]],
                "t_final": 0.01,
            },
            "single viscosity stage",
        ),
        ({"viscosity_schedule": [[0.001, 5e-4]]}, "first viscosity switch time"),
        ({"viscosity_schedule": [[0.0, -1.0]]}, "viscosity values must be"),
        ({"units": {}}, "units must be a nonempty object"),
        ({"description": "  "}, "description must be a nonempty string"),
        ({"interpretation": ""}, "interpretation must be a nonempty string"),
        ({"seed": "abc"}, "seed must be an integer"),
        ({"r_max": 0.0}, "r_max must be positive"),
        ({"z_period": -1.0}, "z_period must be positive"),
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


def test_missing_and_unknown_config_keys_are_rejected() -> None:
    missing = _config()
    del missing["time_steps"]
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(missing)

    for key in ("resolution", "t_final", "acceptance", "viscosity_schedule", "units"):
        partial = _config()
        del partial[key]
        with pytest.raises(ValueError, match="missing or unknown keys"):
            validate_config(partial)

    unknown = _config()
    unknown["filter_strength"] = 0.1
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(unknown)

    not_an_object = ["schema_version", 1]
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(not_an_object)  # type: ignore[arg-type]


def test_invalid_acceptance_blocks_are_rejected() -> None:
    extra = _config()
    extra["acceptance"] = dict(extra["acceptance"])
    extra["acceptance"]["bonus"] = 1.0
    with pytest.raises(ValueError, match="acceptance object has missing or unknown"):
        validate_config(extra)

    missing = _config()
    missing["acceptance"] = dict(missing["acceptance"])
    del missing["acceptance"]["temporal_order_lower_bound"]
    with pytest.raises(ValueError, match="acceptance object has missing or unknown"):
        validate_config(missing)

    inverted = _config()
    inverted["acceptance"] = dict(inverted["acceptance"])
    inverted["acceptance"]["temporal_order_lower_bound"] = 3.0
    with pytest.raises(ValueError, match="temporal_order_lower_bound must be below"):
        validate_config(inverted)

    too_many = _config()
    too_many["acceptance"] = dict(too_many["acceptance"])
    too_many["acceptance"]["minimum_quantities_at_second_order"] = (
        len(SMOOTH_QUANTITIES) + 1
    )
    with pytest.raises(ValueError, match="exceeds the number of analysed"):
        validate_config(too_many)

    not_integer = _config()
    not_integer["acceptance"] = dict(not_integer["acceptance"])
    not_integer["acceptance"]["minimum_quantities_at_second_order"] = 2.0
    with pytest.raises(ValueError, match="must be a positive integer"):
        validate_config(not_integer)


def test_a_viscosity_switch_at_or_after_t_final_is_accepted() -> None:
    """Only a switch *inside* the window is fatal; T_1 sits below t_0."""

    config = _config()
    config["viscosity_schedule"] = [[0.0, 5e-4], [0.00227375, 5e-3]]
    validate_config(config)
    config["viscosity_schedule"] = [[0.0, 5e-4], [config["t_final"], 5e-3]]
    validate_config(config)


# ------------------------------------------------------------ step planning


def test_step_plan_lands_exactly_on_t_final() -> None:
    """Divisible windows keep the requested dt; the rest shorten it uniformly."""

    exact = plan_step(2.0e-5, 2.0e-6)
    assert exact.step_count == 10
    assert exact.requested_dt_divides_t_final is True
    assert exact.effective_dt == pytest.approx(2.0e-6, rel=1e-15)
    assert exact.step_count * exact.effective_dt == pytest.approx(2.0e-5, rel=1e-15)

    # T_1 is not an integer multiple of 6e-7: 3652.88..., so the step count is
    # rounded up and every step shrinks slightly.  Nothing is left ragged.
    clipped = plan_step(T1, 6.0e-7)
    assert clipped.step_count == 3653
    assert clipped.requested_dt_divides_t_final is False
    assert clipped.effective_dt < 6.0e-7
    assert clipped.effective_dt_relative_adjustment < 1.0e-4
    assert clipped.step_count * clipped.effective_dt == pytest.approx(T1, rel=1e-15)


def test_shipped_ladder_step_counts_double_exactly() -> None:
    plans = plan_ladder(T1, [6e-7, 3e-7, 1.5e-7])
    assert [item.step_count for item in plans] == [3653, 7306, 14612]
    assert ladder_halves_exactly(plans) is True
    for coarse, fine in zip(plans, plans[1:]):
        # Halving in binary floating point is exact, so this is an equality.
        assert fine.effective_dt == 0.5 * coarse.effective_dt

    smoke = plan_ladder(2.0e-5, [2.0e-6, 1.0e-6, 5.0e-7])
    assert [item.step_count for item in smoke] == [10, 20, 40]
    assert ladder_halves_exactly(smoke) is True
    assert all(item.requested_dt_divides_t_final for item in smoke)


# ------------------------------------------------------- step-doubling logic


def test_step_doubling_order_recovers_a_second_order_sequence() -> None:
    # Q(dt) = Q* + c dt^2 with c dt^2 = 4e-4: differences 3e-4 and 0.75e-4.
    report = step_doubling_order(
        1.0 + 4.0e-4,
        1.0 + 1.0e-4,
        1.0 + 0.25e-4,
        floor_relative=1.0e-10,
    )
    assert report["status"] == "measured"
    assert report["observed_order"] == pytest.approx(2.0, abs=1e-10)
    assert report["differences_decrease"] is True
    assert report["difference_coarse_middle"] == pytest.approx(3.0e-4)
    assert report["difference_middle_fine"] == pytest.approx(0.75e-4)


def test_step_doubling_order_reports_the_measurement_floor() -> None:
    """A temporal error below the floor is an explicit PASS, not a number."""

    report = step_doubling_order(
        1.0 + 1.0e-13,
        1.0 + 0.25e-13,
        1.0,
        floor_relative=1.0e-10,
    )
    assert report["status"] == "temporal_error_below_measurement_floor"
    assert report["observed_order"] is None
    assert report["measurement_floor"] == pytest.approx(1.0e-10)

    # Identical values are the extreme case and must not divide by zero.
    identical = step_doubling_order(2.0, 2.0, 2.0, floor_relative=0.0)
    assert identical["status"] == "temporal_error_below_measurement_floor"
    assert identical["observed_order"] is None

    non_finite = step_doubling_order(
        float("nan"), 1.0, 1.0, floor_relative=1.0e-10
    )
    assert non_finite["status"] == "non_finite"
    assert non_finite["observed_order"] is None


def test_orders_are_refused_when_a_member_stopped_early() -> None:
    """Step doubling compares members at ONE instant or reports nothing."""

    assert members_share_final_time(
        [{"final_time": 2.0e-5}, {"final_time": 2.0e-5}, {"final_time": 2.0e-5}]
    )
    assert not members_share_final_time(
        [{"final_time": 2.0e-5}, {"final_time": 1.0e-5}, {"final_time": 2.0e-5}]
    )
    assert not members_share_final_time(
        [{"final_time": 2.0e-5}, {"final_time": None}, {"final_time": 2.0e-5}]
    )

    metrics = [
        {"final_time": 2.0e-5, "final_amplification": 1.0, "dr": 0.1, "dz": 0.1},
        {"final_time": 1.0e-5, "final_amplification": 1.0, "dr": 0.1, "dz": 0.1},
        {"final_time": 2.0e-5, "final_amplification": 1.0, "dr": 0.1, "dz": 0.1},
    ]
    analysis = temporal_order_analysis(metrics, floor_relative=1.0e-10)
    assert set(analysis) == {name for name, _ in SMOOTH_QUANTITIES}
    for item in analysis.values():
        assert item["status"] == "members_disagree_on_final_time"
        assert item["observed_order"] is None


def test_cfl_within_tolerance_boundary_is_inclusive() -> None:
    assert cfl_within_tolerance(0.105, cfl_coefficient=0.1, cfl_excess_tolerance=0.05)
    assert not cfl_within_tolerance(
        0.1050001, cfl_coefficient=0.1, cfl_excess_tolerance=0.05
    )
    assert not cfl_within_tolerance(
        float("nan"), cfl_coefficient=0.1, cfl_excess_tolerance=0.05
    )
    assert not cfl_within_tolerance(
        float("inf"), cfl_coefficient=0.1, cfl_excess_tolerance=0.05
    )
    with pytest.raises(ValueError, match="cfl_excess_tolerance must be nonnegative"):
        cfl_within_tolerance(0.1, cfl_coefficient=0.1, cfl_excess_tolerance=-0.1)


# ------------------------------------------------- spatial reference ladder


def test_spatial_ladder_is_read_back_from_the_shipped_summaries() -> None:
    ladder = read_spatial_amplification_ladder(REPOSITORY_ROOT)
    assert ladder["summaries_read"] == list(SPATIAL_LADDER_SUMMARIES)
    assert ladder["summaries_missing"] == []
    assert ladder["matches_recorded_reference"] is True
    assert [item["nr"] for item in ladder["rungs"]] == [65, 129, 193, 257]
    assert all(item["read_from_summary"] for item in ladder["rungs"])
    assert ladder["final_amplifications"] == pytest.approx(
        [6.114768458177527, 12.695695243738214, 15.627954940635405, 17.258778340188396]
    )
    assert ladder["final_amplifications"] == pytest.approx(
        [item["final_amplification"] for item in SPATIAL_AMPLIFICATION_LADDER]
    )
    assert ladder["consecutive_differences"] == pytest.approx(
        [6.580926785560687, 2.932259696897191, 1.630823399552991]
    )
    assert ladder["smallest_consecutive_difference"] == pytest.approx(
        1.630823399552991
    )


# ------------------------------------------------------- end-to-end smoke run


def test_smoke_run_produces_a_complete_accepted_bundle(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, _, summary = generated_run
    assert summary["experiment_id"] == EXPERIMENT_ID
    expected_checks = {
        "all_time_steps_completed",
        "all_diagnostics_finite",
        "effective_time_steps_halve_exactly",
        "energy_non_increasing",
        "circulation_max_principle",
        "temporal_order_second_order_or_below_floor",
        "step_doubling_differences_decrease",
        "fixed_step_within_adaptive_cfl_bound",
        "cross_solver_elliptic_agreement_recorded",
    }
    assert set(summary["acceptance_checks"]) == expected_checks
    assert all(summary["acceptance_checks"].values())
    assert summary["accepted_as_temporal_refinement_study"] is True
    assert "not a proof" in summary["interpretation"]
    assert summary["reproducibility"]["seed"] == config["seed"]
    assert summary["reproducibility"]["runtime_provenance"][
        "source_fingerprint_sha256"
    ]

    limitations = " ".join(summary["limitations"])
    assert "temporal discretization error" in limitations
    assert "numerical observation" in limitations
    assert summary["known_gaps"]

    levels = summary["time_steps"]
    assert len(levels) == EXPECTED_TIME_STEP_COUNT
    assert [item["requested_dt"] for item in levels] == config["time_steps"]
    assert [item["step_count"] for item in levels] == [10, 20, 40]
    for item in levels:
        assert item["completed"] is True
        assert item["failure"] is None
        assert item["all_diagnostics_finite"] is True
        assert item["nr"] == config["resolution"][0]
        assert item["nz"] == config["resolution"][1]
        assert item["final_time"] == pytest.approx(config["t_final"])
        assert item["observed_step_count"] == item["step_count"]
        # The integrator clips only the very last step, and only by the
        # accumulated floating-point drift of t + dt.
        assert item["final_step_relative_deviation"] < 1.0e-9
        assert item["maximum_energy_growth_ratio"] <= 0.0
        assert math.isfinite(item["final_amplification"])
        assert math.isfinite(item["final_cross_solver_psi_max_abs_difference"])
        assert item["maximum_wall_u1_abs"] == 0.0

    # One spatial grid for all three members: that is what makes the mutual
    # differences purely temporal.
    grids = {(item["nr"], item["nz"], item["dr"], item["dz"]) for item in levels}
    assert len(grids) == 1


def test_temporal_orders_are_recorded_as_numbers_or_the_floor_marker(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, _, summary = generated_run
    analysis = summary["temporal_analysis"]
    orders = analysis["temporal_order"]
    assert set(orders) == {name for name, _ in SMOOTH_QUANTITIES}
    lower = config["acceptance"]["temporal_order_lower_bound"]
    upper = config["acceptance"]["temporal_order_upper_bound"]
    assert analysis["temporal_order_band"] == [lower, upper]
    assert analysis["members_share_final_time"] is True

    for name, source in SMOOTH_QUANTITIES:
        item = orders[name]
        assert item["diagnostic_key"] == source
        assert item["status"] in {
            "measured",
            "temporal_error_below_measurement_floor",
        }
        assert math.isfinite(item["difference_coarse_middle"])
        assert math.isfinite(item["difference_middle_fine"])
        assert item["measurement_floor"] == pytest.approx(
            config["acceptance"]["temporal_error_floor_relative"]
            * abs(item["fine_value"])
        )
        if item["status"] == "measured":
            assert isinstance(item["observed_order"], float)
            assert item["differences_decrease"] is True
        else:
            assert item["observed_order"] is None
            # The explicit explanation for the non-numeric outcome.
            assert min(
                item["difference_coarse_middle"], item["difference_middle_fine"]
            ) <= item["measurement_floor"]

    measured = [
        name for name, item in orders.items() if item["status"] == "measured"
    ]
    assert analysis["measured_quantity_count"] == len(measured)
    assert analysis["quantities_required_within_the_band"] == min(
        config["acceptance"]["minimum_quantities_at_second_order"], len(measured)
    )
    # Heun is second order: at least two of the four smooth quantities must
    # land in the band unless every one of them sits below the floor.
    assert (
        len(analysis["quantities_within_the_order_band"])
        >= analysis["quantities_required_within_the_band"]
    )
    assert set(analysis["quantities_within_the_order_band"]) <= set(measured)
    assert set(analysis["quantities_below_measurement_floor"]) == (
        {name for name, _ in SMOOTH_QUANTITIES} - set(measured)
    )
    assert "log2" in analysis["temporal_order_definition"]
    assert "Q(dt/4)" in analysis["measurement_floor_definition"]


def test_energy_is_non_increasing_at_every_time_step(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    _, output, summary = generated_run
    assert summary["acceptance_checks"]["energy_non_increasing"] is True
    with (output / "diagnostics.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    by_level: dict[str, list[float]] = {}
    for row in rows:
        by_level.setdefault(row["level_index"], []).append(float(row["energy"]))
    assert len(by_level) == EXPECTED_TIME_STEP_COUNT
    for series in by_level.values():
        assert len(series) > 1
        assert all(later <= earlier for earlier, later in zip(series, series[1:]))


def test_argmax_locations_and_constraint_trends_are_recorded(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    _, _, summary = generated_run
    analysis = summary["temporal_analysis"]
    argmax = analysis["argmax_locations"]
    assert set(argmax) == {name for name, _, _ in ARGMAX_QUANTITIES}
    for item in argmax.values():
        assert item["status"] == "recorded"
        assert len(item["r"]) == EXPECTED_TIME_STEP_COUNT
        assert len(item["z"]) == EXPECTED_TIME_STEP_COUNT
        assert all(0.0 <= value <= 1.0 for value in item["r"])
        assert all(0.0 <= value < 1.0 for value in item["z"])
        # Reported, never forced.  On one shared grid they are expected to be
        # identical, and agreement within one cell is the weaker claim.
        assert item["agree_within_one_grid_cell"] is True
        assert isinstance(item["identical_across_the_ladder"], bool)

    trends = analysis["constraint_metric_trends"]
    assert set(trends) == set(CONSTRAINT_QUANTITIES)
    for item in trends.values():
        assert item["status"] == "recorded"
        assert len(item["values"]) == EXPECTED_TIME_STEP_COUNT
        assert all(math.isfinite(value) for value in item["values"])
        assert item["absolute_spread"] >= 0.0
        assert item["relative_spread"] >= 0.0


def test_summary_separates_temporal_from_spatial_error(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    _, _, summary = generated_run
    block = summary["temporal_vs_spatial_separation"]
    # Promoted to the top level of the summary, not duplicated inside it.
    assert "temporal_vs_spatial_separation" not in summary["temporal_analysis"]
    assert "no spatial error at all" in block["statement"]
    assert block["fixed_resolution"] == "17x32"

    temporal = block["temporal"]
    assert temporal["quantity"] == "final_amplification"
    assert math.isfinite(
        temporal["step_doubling_difference_dt_over_2_to_dt_over_4"]
    )
    assert temporal["richardson_error_estimate_at_finest_dt"] == pytest.approx(
        temporal["step_doubling_difference_dt_over_2_to_dt_over_4"] / 3.0
    )

    spatial = block["spatial"]
    assert spatial["matches_recorded_reference"] is True
    assert spatial["final_amplifications"] == pytest.approx(
        [6.114768458177527, 12.695695243738214, 15.627954940635405, 17.258778340188396]
    )
    assert spatial["consecutive_differences"] == pytest.approx(
        [6.580926785560687, 2.932259696897191, 1.630823399552991]
    )
    assert set(spatial["summaries_read"]) == set(SPATIAL_LADDER_SUMMARIES)

    # The whole point of the study: on this grid the spatial ladder difference
    # dwarfs the temporal one.
    assert block["dominant_error_source"] == "spatial"
    assert block["spatial_to_temporal_ratio"] > 1.0
    assert block["spatial_to_temporal_ratio"] == pytest.approx(
        spatial["smallest_consecutive_difference"]
        / temporal["richardson_error_estimate_at_finest_dt"]
    )

    constraints = block["constraint_metrics_expected_spatially_dominated"]
    assert constraints["gated"] is False
    assert "independent of dt" in constraints["expectation"]
    assert set(constraints["metrics"]) == set(CONSTRAINT_QUANTITIES)
    # Spatially dominated means essentially dt independent; this is recorded
    # and asserted here rather than gated inside the experiment.
    for item in constraints["metrics"].values():
        assert item["relative_spread"] < 1.0e-3


def test_cfl_is_monitored_and_within_the_documented_tolerance(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, _, summary = generated_run
    policy = summary["cfl_policy"]
    assert policy["monitoring_only"] is True
    assert "no adaptive selection" in policy["mode"]
    assert policy["cfl_coefficient"] == config["cfl_coefficient"]
    assert policy["cfl_excess_tolerance"] == 0.05
    assert policy["cfl_excess_tolerance_source"].startswith("default")
    assert policy["accepted_effective_cfl_bound"] == pytest.approx(
        config["cfl_coefficient"] * 1.05
    )
    assert (
        policy["maximum_effective_advective_cfl"]
        <= policy["accepted_effective_cfl_bound"]
    )
    assert (
        policy["maximum_effective_viscous_cfl"]
        <= policy["accepted_effective_cfl_bound"]
    )
    assert policy["maximum_effective_advective_cfl"] == max(
        item["maximum_advective_cfl"] for item in summary["time_steps"]
    )
    assert "never raised" in policy["exceedance_policy"]
    assert summary["acceptance_checks"]["fixed_step_within_adaptive_cfl_bound"] is True


def test_smoke_run_artifacts_are_finite_and_checksummed(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, _ = generated_run
    with (output / "diagnostics.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    assert {"level_index", "requested_dt", "effective_dt"} <= rows[0].keys()
    assert set(DIAGNOSTIC_FIELDS) <= rows[0].keys()
    assert set(RELATIVE_DIAGNOSTIC_FIELDS) <= rows[0].keys()
    assert all(math.isfinite(float(value)) for row in rows for value in row.values())
    assert {row["level_index"] for row in rows} == {"0", "1", "2"}

    with (output / "time_refinement.csv").open(newline="", encoding="utf-8") as stream:
        levels = list(csv.DictReader(stream))
    assert len(levels) == EXPECTED_TIME_STEP_COUNT
    assert {
        "requested_dt",
        "effective_dt",
        "step_count",
        "final_amplification",
        "final_max_abs_u1",
        "final_energy",
        "final_argmax_u1_r",
        "final_argmax_u1_z",
        "final_divergence_residual_relative",
        "final_axis_parity_relative_u1",
        "final_cross_solver_psi_relative_difference",
        "maximum_circulation_growth_ratio",
    } <= levels[0].keys()

    with np.load(output / "trajectories.npz", allow_pickle=False) as archive:
        assert {"r", "z", "requested_dt", "effective_dt"} <= set(archive.files)
        assert any(name.startswith("psi1_final_dt") for name in archive.files)
        assert all(np.all(np.isfinite(archive[name])) for name in archive.files)

    manifest = strict_json_loads(
        (output / "manifest.json").read_text(encoding="utf-8"),
        label="hou time refinement manifest",
    )
    assert manifest["experiment_id"] == EXPERIMENT_ID
    assert {
        "config.snapshot.json",
        "diagnostics.csv",
        "summary.json",
        "time_refinement.csv",
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
        label="hou time refinement config snapshot",
    )
    assert snapshot == config
    assert json.loads((output / "summary.json").read_text(encoding="utf-8"))


def test_final_state_checkpoints_reload_as_v2_candidates(
    generated_run: tuple[dict[str, Any], Path, dict[str, Any]],
) -> None:
    config, output, summary = generated_run
    archives = sorted((output / "checkpoints").glob("*.npz"))
    assert len(archives) == EXPECTED_TIME_STEP_COUNT
    assert [item["checkpoint"] for item in summary["time_steps"]] == [
        path.name for path in archives
    ]
    artifact = load_candidate(archives[-1])
    assert set(artifact.fields) == {"u1", "omega1", "psi1"}
    assert artifact.schema.endswith("/v2")
    assert artifact.description is not None
    assert artifact.provenance is not None

    checkpoint = load_checkpoint(archives[-1])
    # The stored psi1 is the elliptic solution of the stored omega1, so
    # re-constraining the checkpoint must reproduce it bit for bit.
    assert np.array_equal(checkpoint.state.psi1, artifact.fields["psi1"])
    assert checkpoint.time == pytest.approx(config["t_final"])
    assert checkpoint.viscosity == pytest.approx(config["viscosity_schedule"][0][1])
    assert checkpoint.metadata["experiment_id"] == EXPERIMENT_ID
    assert checkpoint.metadata["step_count"] == 40
    assert np.all(np.isfinite(checkpoint.state.u1))
    assert float(np.max(np.abs(checkpoint.state.u1[-1]))) == 0.0


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


def test_max_steps_records_an_explicit_incomplete_member() -> None:
    """An exhausted step budget is recorded, not raised."""

    config = _smoke_config()
    config["max_steps"] = 3
    plan = plan_step(float(config["t_final"]), float(config["time_steps"][0]))
    result = evolve_time_step(config, index=0, plan=plan)
    assert result.completed is False
    assert result.failure is not None
    assert "before t_final" in result.failure
    assert result.observed_step_count == 3
