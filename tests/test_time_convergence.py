from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.run_baseline import analytic_u1, evolve
from experiments.run_time_convergence import (
    evaluate,
    run,
    step_doubling_diagnostics,
    temporal_orders,
)
from ns_certificate_lab._integrity import strict_json_loads


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "baseline_time_convergence.json"


def _config() -> dict[str, object]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="test time convergence config",
    )
    assert isinstance(value, dict)
    return value


def _independent_relative_grid_norm(
    numerical: np.ndarray,
    reference: np.ndarray,
) -> float:
    return float(np.linalg.norm(numerical - reference) / np.linalg.norm(reference))


def test_fixed_grid_crank_nicolson_has_temporal_convergence() -> None:
    """Check convergence without using the experiment's error/order helpers."""

    config = _config()
    nr = int(config["nr"])
    r_max = float(config["r_max"])
    dr = r_max / (nr - 1)
    final_time = float(config["final_time"])
    base_dt = float(config["base_time_step"])
    errors: list[float] = []
    grids: list[np.ndarray] = []
    final_profiles: list[np.ndarray] = []
    for divisor in (1, 2, 4):
        requested_dt = base_dt / divisor
        result = evolve(
            nr=nr,
            r_max=r_max,
            final_time=final_time,
            viscosity=float(config["viscosity"]),
            amplitude=float(config["amplitude"]),
            sigma=float(config["sigma"]),
            dt_over_dr=requested_dt / dr,
            sample_count=17,
            exact_outer_boundary=True,
        )
        exact = analytic_u1(
            result.r,
            final_time,
            viscosity=float(config["viscosity"]),
            amplitude=float(config["amplitude"]),
            sigma=float(config["sigma"]),
        )
        grids.append(result.r)
        final_profiles.append(result.profiles[-1])
        errors.append(_independent_relative_grid_norm(result.profiles[-1], exact))

    assert np.array_equal(grids[0], grids[1])
    assert np.array_equal(grids[0], grids[2])
    assert errors[0] > errors[1] > errors[2]
    independent_orders = [
        math.log(errors[index] / errors[index + 1], 2.0)
        for index in range(2)
    ]
    assert min(independent_orders) > 1.8
    step_differences = [
        np.linalg.norm(final_profiles[index] - final_profiles[index + 1])
        for index in range(2)
    ]
    assert math.log(step_differences[0] / step_differences[1], 2.0) > 1.8


def test_evaluate_records_requested_diagnostics_and_exact_dt_sequence() -> None:
    config = _config()
    runs, orders = evaluate(config)
    assert [item.actual_dt for item in runs] == pytest.approx(
        [0.5, 0.25, 0.125],
        rel=0.0,
        abs=1e-15,
    )
    assert [item.step_count for item in runs] == [2, 4, 8]
    assert min(orders) > 1.8
    reference = analytic_u1(
        runs[0].evolution.r,
        float(config["final_time"]),
        viscosity=float(config["viscosity"]),
        amplitude=float(config["amplitude"]),
        sigma=float(config["sigma"]),
    )
    differences, step_orders = step_doubling_diagnostics(runs, reference)
    assert differences[0] > differences[1]
    assert min(step_orders) > 1.8
    for item in runs:
        assert item.evolution.energy_per_z_length[-1] <= (
            item.evolution.energy_per_z_length[0] * (1.0 + 1e-12)
        )
        assert np.all(np.isfinite(item.evolution.peak_vorticity))
        assert item.auxiliary_boundary_all_time_inner_max_difference < 5e-8
        assert item.auxiliary_boundary_final_inner_max_difference < 5e-8
        assert (
            item.auxiliary_boundary_all_time_inner_max_difference
            >= item.auxiliary_boundary_final_inner_max_difference
        )


def test_time_convergence_artifacts_are_complete_and_checksummed(
    scratch_dir: Path,
) -> None:
    config = _config()
    output = scratch_dir / "time-convergence"
    summary = run(config, output)
    assert summary["accepted_as_temporal_control"] is True
    assert all(summary["acceptance_checks"].values())

    loaded_summary = strict_json_loads(
        (output / "summary.json").read_text(encoding="utf-8"),
        label="generated time convergence summary",
    )
    assert loaded_summary["temporal_convergence"]["actual_time_steps"] == [
        0.5,
        0.25,
        0.125,
    ]
    assert len(
        loaded_summary["temporal_convergence"][
            "step_doubling_relative_differences"
        ]
    ) == 2
    assert len(
        loaded_summary["temporal_convergence"]["step_doubling_observed_orders"]
    ) == 1
    assert len(loaded_summary["energy_and_vorticity"]) == 3
    boundary_summary = loaded_summary["auxiliary_boundary_radius_test"]
    assert len(
        boundary_summary["inner_differences_by_dt"]
    ) == 3
    assert "not directly" in boundary_summary["limitation"]
    assert boundary_summary["main_domain_r_max"] == 5.0
    assert boundary_summary["radii"] == [3.0, 4.0]

    with (output / "time_convergence.csv").open(
        newline="",
        encoding="utf-8",
    ) as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3
    required_columns = {
        "actual_dt",
        "relative_l2_error",
        "analytic_error_order_from_previous",
        "step_doubling_difference_to_next",
        "step_doubling_order_from_previous_difference",
        "initial_energy_per_unit_z",
        "final_energy_per_unit_z",
        "maximum_vorticity_over_time",
        "final_vorticity",
        "auxiliary_boundary_all_time_inner_max_difference",
        "auxiliary_boundary_final_inner_max_difference",
    }
    assert required_columns <= rows[0].keys()
    assert [float(row["actual_dt"]) for row in rows] == [0.5, 0.25, 0.125]

    with (output / "time_diagnostics.csv").open(
        newline="",
        encoding="utf-8",
    ) as stream:
        diagnostic_rows = list(csv.DictReader(stream))
    expected_steps = {0.5: 2, 0.25: 4, 0.125: 8}
    assert len(diagnostic_rows) == sum(value + 1 for value in expected_steps.values())
    times_by_dt: dict[float, list[float]] = {
        time_step: [] for time_step in expected_steps
    }
    for row in diagnostic_rows:
        values = {
            key: float(row[key])
            for key in ("dt", "time", "energy_per_unit_z", "peak_vorticity")
        }
        assert all(math.isfinite(value) for value in values.values())
        assert values["dt"] in times_by_dt
        times_by_dt[values["dt"]].append(values["time"])
    for time_step, step_count in expected_steps.items():
        times = times_by_dt[time_step]
        assert len(times) == step_count + 1
        assert times[0] == pytest.approx(0.0, abs=1e-15)
        assert times[-1] == pytest.approx(float(config["final_time"]), abs=1e-15)
        assert np.diff(times) == pytest.approx(
            np.full(step_count, time_step),
            rel=0.0,
            abs=2e-15,
        )

    with np.load(output / "final_profiles.npz", allow_pickle=False) as archive:
        assert set(archive.files) == {"r", "dt", "u1_final", "u1_exact"}
        assert archive["r"].shape == (int(config["nr"]),)
        assert archive["u1_final"].shape == (3, int(config["nr"]))
        assert np.array_equal(archive["dt"], np.array([0.5, 0.25, 0.125]))
        assert all(np.all(np.isfinite(archive[name])) for name in archive.files)

    manifest = strict_json_loads(
        (output / "manifest.json").read_text(encoding="utf-8"),
        label="generated time convergence manifest",
    )
    expected_manifest_files = {
        "config.snapshot.json",
        "final_profiles.npz",
        "summary.json",
        "time_convergence.csv",
        "time_diagnostics.csv",
    }
    assert set(manifest["files"]) == expected_manifest_files
    for name, record in manifest["files"].items():
        data = (output / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == record["sha256"]
        assert len(data) == record["bytes"]
    manifest_digest = hashlib.sha256((output / "manifest.json").read_bytes()).hexdigest()
    assert (
        output / "manifest.json.sha256"
    ).read_text(encoding="ascii").strip() == manifest_digest


def test_time_convergence_rejects_invalid_sequence_and_overwrite(
    scratch_dir: Path,
) -> None:
    config = _config()
    invalid = json.loads(json.dumps(config))
    invalid["time_step_divisors"] = [1, 3, 4]
    with pytest.raises(ValueError, match=r"exactly \[1, 2, 4\]"):
        evaluate(invalid)

    output = scratch_dir / "existing"
    output.mkdir()
    marker = output / "preserve.txt"
    marker.write_text("keep\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(config, output)
    assert marker.read_text(encoding="utf-8") == "keep\n"


def test_temporal_order_helper_rejects_non_refining_inputs() -> None:
    assert temporal_orders([4.0, 1.0, 0.25], [0.4, 0.2, 0.1]) == pytest.approx(
        [2.0, 2.0]
    )
    with pytest.raises(ValueError, match="strictly decreasing"):
        temporal_orders([1.0, 0.5], [0.1, 0.2])
