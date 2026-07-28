"""Tests for the pressure-free primitive Navier-Stokes residual experiment.

The smoke run genuinely integrates the E-29 datum, saves the fan of schema-v2
candidates and audits them back off disk, so every fault injection below acts
on real saved data rather than on a manufactured stand-in.  Grids are
deliberately tiny here; the shipped config uses 129x256 with a 64x64x128 audit
box and is run by the orchestrator, not by this module.

Recorded measurements that the assertions below are calibrated against (small
smoke shape ``(20,20,32)``, ``interior_margin`` 3, 33x64 source grid):

* relative curl-of-momentum residual RMS 3.23e-3, relative
  vorticity-transport residual RMS 8.79e-3;
* observed order of the centered time difference 1.99998;
* fault ratios: corrupted neighbour state 19.8x, wrong dt 7.5x.

The bounds asserted are deliberately looser than those numbers so that they
pin the order of magnitude, not the exact build.
"""

from __future__ import annotations

import ast
import copy
import dataclasses
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import experiments.run_hou_primitive_residual as experiment
from experiments.run_hou_primitive_residual import (
    EXPECTED_COORDINATE_SYSTEM,
    EXPECTED_EXPERIMENT_ID,
    STATE_OFFSETS,
    build_cartesian_grid,
    evaluate_on_grid,
    plan_states,
    reload_states,
    run,
    state_filename,
    validate_config,
)
from ns_certificate_lab import cartesian_momentum_audit as checker
from ns_certificate_lab._integrity import (
    sha256_bytes,
    sha256_file,
    strict_json_loads,
    verify_digest,
)
from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.cartesian_candidate_adapter import (
    reconstruct_loaded_candidate_on_cartesian,
)
from ns_certificate_lab.cartesian_momentum_audit import (
    adjoint_first_derivative,
    curl_momentum_residual,
    gradient_null_space,
    recover_pressure,
    region_window,
    rms,
    vorticity_transport_residual,
)
from ns_certificate_lab.cartesian_validation import (
    UniformCartesianGrid,
    cartesian_curl,
    cartesian_divergence,
    cartesian_gradient,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "hou_primitive_residual.json"

SMOKE_GRID_SPEC: dict[str, Any] = {
    "half_width": 0.7,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_periodic": True,
}
SMOKE_AUDIT_SHAPE = (20, 20, 32)
SMOKE_MARGIN = 3
SMOKE_VISCOSITY = 5.0e-4
NO_PRESSURE_RECOVERY: dict[str, Any] = {
    "enabled": False,
    "tolerance": 1.0e-9,
    "max_iterations": 100,
}

# Loose order-of-magnitude bounds; see the module docstring for the measured
# values these were derived from.
SMOKE_CURL_RESIDUAL_BOUND = 3.0e-2
SMOKE_TRANSPORT_RESIDUAL_BOUND = 6.0e-2


def _shipped_config() -> dict[str, Any]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="hou primitive residual config",
    )
    assert isinstance(value, dict)
    return value


def _smoke_config() -> dict[str, Any]:
    config = _shipped_config()
    config["resolution"] = [33, 64]
    config["time_step"] = 1.0e-6
    config["target_time"] = 8.0e-6
    config["max_steps"] = 2000
    config["diagnostic_stride"] = 1
    config["cartesian_grid"] = {
        **config["cartesian_grid"],
        "shape": list(SMOKE_AUDIT_SHAPE),
        "coarse_shape": [12, 12, 16],
        "interior_margin": SMOKE_MARGIN,
    }
    config["pressure_recovery"] = {
        "enabled": True,
        "tolerance": 1.0e-9,
        "max_iterations": 5000,
    }
    # The smoke audit box is far coarser than the shipped one, so its own
    # truncation is larger; the gates are relaxed to match, never the reverse.
    config["hard_gates"] = {
        "maximum_relative_curl_momentum_residual_rms": SMOKE_CURL_RESIDUAL_BOUND,
        "maximum_relative_vorticity_transport_residual_rms": (
            SMOKE_TRANSPORT_RESIDUAL_BOUND
        ),
        "time_difference_order_lower_bound": 1.8,
        "time_difference_order_upper_bound": 2.2,
    }
    config["record_only_thresholds"] = {
        "reference_relative_curl_momentum_residual_max": 0.3,
        "reference_relative_vorticity_transport_residual_max": 0.6,
        "reference_relative_pressure_consistency_max": 0.3,
        "reference_relative_divergence_rms": 0.05,
    }
    return config


# --------------------------------------------------------------------------
# configuration contract
# --------------------------------------------------------------------------


def test_shipped_config_is_valid() -> None:
    config = _shipped_config()
    validate_config(config)
    assert config["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert config["coordinate_system"] == EXPECTED_COORDINATE_SYSTEM
    plan = plan_states(
        target_time=float(config["target_time"]),
        time_step=float(config["time_step"]),
    )
    assert plan.steps_to_target == 500
    assert plan.step_count == 504
    assert len(set(plan.snapshot_times())) == len(STATE_OFFSETS)


@pytest.mark.parametrize(
    "key",
    (
        "schema_version",
        "experiment_id",
        "description",
        "interpretation",
        "coordinate_system",
        "normalization_description",
        "units",
        "seed",
        "r_max",
        "z_period",
        "resolution",
        "amplitude_scale",
        "viscosity_schedule",
        "time_step",
        "target_time",
        "cfl_coefficient",
        "max_steps",
        "diagnostic_stride",
        "cartesian_grid",
        "pressure_recovery",
        "hard_gates",
        "record_only_thresholds",
    ),
)
def test_config_rejects_missing_key(key: str) -> None:
    config = _shipped_config()
    del config[key]
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(config)


def test_config_rejects_unknown_key() -> None:
    config = _shipped_config()
    config["extra_setting"] = 1
    with pytest.raises(ValueError, match="missing or unknown keys"):
        validate_config(config)


def test_config_accepts_the_documented_optional_key() -> None:
    config = _shipped_config()
    config["cfl_excess_tolerance"] = 0.05
    validate_config(config)
    config["cfl_excess_tolerance"] = -1.0
    with pytest.raises(ValueError, match="nonnegative"):
        validate_config(config)


def test_config_rejects_mislabeled_experiment_id() -> None:
    config = _shipped_config()
    config["experiment_id"] = "hou_snapshot_cartesian_audit_v1"
    with pytest.raises(ValueError, match="experiment_id must be"):
        validate_config(config)


def test_config_rejects_mislabeled_coordinate_system() -> None:
    config = _shipped_config()
    config["coordinate_system"] = "axisymmetric cylindrical half-plane (r,z)"
    with pytest.raises(ValueError, match="coordinate_system"):
        validate_config(config)


def test_config_rejects_unsupported_schema_version() -> None:
    config = _shipped_config()
    config["schema_version"] = 2
    with pytest.raises(ValueError, match="schema_version"):
        validate_config(config)


@pytest.mark.parametrize("value", (0.0, -1.0e-6, "1e-6", True))
def test_config_rejects_non_positive_time_step(value: Any) -> None:
    config = _shipped_config()
    config["time_step"] = value
    with pytest.raises(ValueError, match="time_step"):
        validate_config(config)


def test_config_rejects_target_time_that_is_not_a_whole_number_of_steps() -> None:
    config = _shipped_config()
    config["target_time"] = 5.0e-4 + 3.0e-7
    with pytest.raises(ValueError, match="integer multiple of time_step"):
        validate_config(config)


def test_config_rejects_target_time_too_close_to_the_initial_time() -> None:
    config = _shipped_config()
    config["target_time"] = 4.0 * float(config["time_step"])
    with pytest.raises(ValueError, match="not reachable"):
        validate_config(config)


def test_plan_states_rejects_an_unreachable_target_directly() -> None:
    with pytest.raises(ValueError, match="not reachable"):
        plan_states(target_time=2.0e-6, time_step=1.0e-6)
    with pytest.raises(ValueError, match="integer multiple"):
        plan_states(target_time=7.5e-6, time_step=1.0e-6)


def test_config_rejects_a_step_budget_below_the_plan() -> None:
    config = _shipped_config()
    config["max_steps"] = 10
    with pytest.raises(ValueError, match="max_steps is below"):
        validate_config(config)


def test_config_rejects_a_viscosity_switch_inside_the_integration() -> None:
    config = _shipped_config()
    config["viscosity_schedule"] = [[0.0, 5.0e-4], [1.0e-4, 5.0e-3]]
    with pytest.raises(ValueError, match="single viscosity stage"):
        validate_config(config)


def test_config_rejects_a_coarse_shape_that_is_not_coarser() -> None:
    config = _shipped_config()
    config["cartesian_grid"]["coarse_shape"] = config["cartesian_grid"]["shape"]
    with pytest.raises(ValueError, match="strictly coarser"):
        validate_config(config)


def test_config_rejects_an_interior_margin_that_empties_the_region() -> None:
    config = _shipped_config()
    config["cartesian_grid"]["coarse_shape"] = [8, 8, 16]
    config["cartesian_grid"]["interior_margin"] = 4
    with pytest.raises(ValueError, match="empty norm region"):
        validate_config(config)


def test_config_rejects_a_box_whose_corners_leave_the_cylinder() -> None:
    config = _shipped_config()
    config["cartesian_grid"]["half_width"] = 0.8
    with pytest.raises(ValueError, match="leave the source cylinder"):
        validate_config(config)


def test_config_rejects_an_inverted_order_band() -> None:
    config = _shipped_config()
    config["hard_gates"]["time_difference_order_lower_bound"] = 3.0
    with pytest.raises(ValueError, match="order_lower_bound must be below"):
        validate_config(config)


def test_config_rejects_unknown_pressure_recovery_keys() -> None:
    config = _shipped_config()
    config["pressure_recovery"]["extra"] = 1
    with pytest.raises(ValueError, match="pressure_recovery has missing"):
        validate_config(config)


def test_cli_help_text_is_printable_on_a_legacy_console() -> None:
    """argparse prints the module docstring; a non-ASCII glyph crashes cp932."""

    assert experiment.__doc__ is not None
    experiment.__doc__.encode("ascii")
    with pytest.raises(SystemExit) as exit_info:
        experiment.main(["--help"])
    assert exit_info.value.code == 0


def test_evaluation_module_never_imports_a_cylindrical_operator_module() -> None:
    """The checker must be a second, independent path (TM-14)."""

    tree = ast.parse(inspect.getsource(checker))
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            modules.add(node.module or "")
            if node.level:
                modules.add(f".{node.module or ''}")
    forbidden = {
        "ns_certificate_lab.operators",
        "ns_certificate_lab.pde",
        "ns_certificate_lab.poisson",
        "ns_certificate_lab.finite_cylinder_poisson",
        "ns_certificate_lab.nonlinear_cylinder",
        "ns_certificate_lab.cartesian",
        "operators",
        "pde",
        "poisson",
        "finite_cylinder_poisson",
        "nonlinear_cylinder",
        ".operators",
        ".pde",
        ".poisson",
        ".finite_cylinder_poisson",
        ".nonlinear_cylinder",
    }
    assert not (modules & forbidden), sorted(modules & forbidden)


def test_evaluation_reads_the_saved_artifacts_not_the_in_memory_state() -> None:
    """Non-circularity: the evaluation path must start at load_candidate."""

    source = inspect.getsource(experiment)
    body = source.split("def reload_states", 1)[1].split("\ndef ", 1)[0]
    assert "load_candidate" in body
    evaluation = source.split("def evaluate_on_grid", 1)[1].split("\ndef ", 1)[0]
    for forbidden in ("integrate(", "heun_step", "ConstrainedState("):
        assert forbidden not in evaluation


# --------------------------------------------------------------------------
# smoke run: produce real states, then audit them off disk
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def smoke_bundle():
    """One real smoke run, shared by the structural and fault-injection tests.

    The experiment refuses to write outside the repository, so this mirrors the
    ``scratch_dir`` policy of ``tests/conftest.py`` at module scope instead of
    using pytest's system temporary directory.
    """

    import shutil
    import uuid

    root = (REPOSITORY_ROOT / "tests" / "_scratch").resolve()
    output_dir = (root / f"primitive_residual_{uuid.uuid4().hex}").resolve()
    if root not in output_dir.parents:
        raise RuntimeError("refusing to create a directory outside tests/_scratch")
    root.mkdir(parents=True, exist_ok=True)
    try:
        yield {"summary": run(_smoke_config(), output_dir), "output_dir": output_dir}
    finally:
        if output_dir.exists():
            shutil.rmtree(output_dir)
        if root.exists() and not any(root.iterdir()):
            root.rmdir()


def test_smoke_run_structure_and_manifest(smoke_bundle) -> None:
    summary = smoke_bundle["summary"]
    output_dir: Path = smoke_bundle["output_dir"]

    assert summary["schema_version"] == 1
    assert summary["experiment_id"] == EXPECTED_EXPERIMENT_ID
    assert set(summary["denominators"]) == set(checker.DENOMINATORS)
    assert summary["limitations"] and summary["known_gaps"]
    assert summary["pressure_recovery_policy"]["record_only"] is True
    assert "NOT A VALIDATED PRESSURE" in summary["pressure_recovery_policy"]["scope"]
    assert "still_open" in summary["closes_recorded_gap"]

    checks = summary["acceptance_checks"]
    for name in (
        "integration_completed",
        "every_requested_state_was_saved",
        "all_states_finite",
        "all_metrics_finite",
        "time_difference_order_within_band",
        "time_difference_errors_decrease_with_spacing",
        "relative_curl_momentum_residual_rms_within_gate",
        "relative_vorticity_transport_residual_rms_within_gate",
        "residuals_decrease_under_audit_grid_refinement",
        "fixed_step_within_adaptive_cfl_bound",
    ):
        assert checks[name] is True, name
    assert summary["accepted_as_pressure_free_primitive_residual_study"] is True
    assert set(summary["record_only_checks"]) == {
        "relative_curl_momentum_residual_max_below_reference",
        "relative_vorticity_transport_residual_max_below_reference",
        "relative_divergence_rms_below_reference",
        "pressure_recovery_converged",
        "relative_pressure_consistency_max_below_reference",
    }

    # every saved state is present, reloaded and digest-recorded
    assert set(summary["production"]["saved_states"]) == {
        str(offset) for offset in STATE_OFFSETS
    }
    for offset in STATE_OFFSETS:
        record = summary["reloaded_states"][str(offset)]
        target = output_dir / "checkpoints" / state_filename(offset)
        assert record["archive_sha256"] == sha256_file(target)
        assert math.isfinite(record["physical_time"])

    # every relative metric and every denominator is present and finite
    for evaluation in summary["evaluations"].values():
        for block in evaluation["regions"].values():
            for name in (
                "curl_momentum_residual",
                "vorticity_transport_residual",
                "form_difference",
                "divergence",
            ):
                report = block[name]
                for key in ("rms", "max", "term_sum_max"):
                    assert math.isfinite(report[key]), (name, key)
                assert report["term_sum_max"] > 0.0
                for key in ("relative_rms", "relative_max"):
                    assert report[key] is not None
                    assert math.isfinite(report[key]), (name, key)
                assert {"x", "y", "z", "radius", "index", "value"} <= set(
                    report["argmax"]
                )
            for value in block["scales"].values():
                assert math.isfinite(value)
            fractions = block["term_magnitudes"]["curl_momentum"][
                "fraction_of_term_rms_sum"
            ]
            assert abs(sum(fractions.values()) - 1.0) < 1.0e-9
        richardson = evaluation["time_difference_richardson"]
        assert 1.8 <= float(richardson["observed_order"]) <= 2.2
        assert richardson["differences_decrease"] is True
        assert evaluation["pressure_recovery"]["converged"] is True

    interior = summary["evaluations"]["fine"]["regions"]["interior"]
    assert 0.0 < interior["curl_momentum_residual"]["relative_rms"] < (
        SMOKE_CURL_RESIDUAL_BOUND
    )
    assert 0.0 < interior["vorticity_transport_residual"]["relative_rms"] < (
        SMOKE_TRANSPORT_RESIDUAL_BOUND
    )
    # The two forms are genuinely different discretizations, so their
    # difference is nonzero and shrinks under audit refinement.
    assert interior["form_difference"]["rms"] > 0.0
    refinement = summary["audit_grid_refinement"]
    assert refinement["spacing_ratio"] > 1.0
    for name in ("curl_momentum_residual", "vorticity_transport_residual"):
        assert refinement[name]["decreases_under_refinement"] is True

    names = sorted(path.name for path in output_dir.iterdir())
    assert names == [
        "checkpoints",
        "config.snapshot.json",
        "config.snapshot.json.sha256",
        "manifest.json",
        "manifest.json.sha256",
        "metrics.csv",
        "summary.json",
        "summary.json.sha256",
    ]
    verify_digest(output_dir / "manifest.json")
    verify_digest(output_dir / "summary.json")
    verify_digest(output_dir / "config.snapshot.json")
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    for name, entry in manifest["files"].items():
        target = output_dir / name
        assert entry["sha256"] == sha256_file(target), name
        assert entry["bytes"] == target.stat().st_size, name
    for offset in STATE_OFFSETS:
        assert f"checkpoints/{state_filename(offset)}" in manifest["files"]
    snapshot = (output_dir / "config.snapshot.json").read_bytes()
    assert summary["reproducibility"]["config_sha256"] == sha256_bytes(snapshot)


def test_run_refuses_a_nonempty_output_directory(scratch_dir: Path) -> None:
    output_dir = scratch_dir / "occupied"
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(_smoke_config(), output_dir)
    assert (output_dir / "existing.txt").read_text(encoding="utf-8") == "keep me"


def test_run_refuses_to_write_outside_the_repository(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="inside this repository"):
        run(_smoke_config(), tmp_path / "elsewhere")


# --------------------------------------------------------------------------
# fault injection on the real saved states
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def fault_baseline(smoke_bundle):
    """Reload the real smoke states once and evaluate them unperturbed."""

    checkpoint_dir = smoke_bundle["output_dir"] / "checkpoints"
    grid = build_cartesian_grid(SMOKE_GRID_SPEC, SMOKE_AUDIT_SHAPE)
    states = reload_states(checkpoint_dir, grid)
    baseline = evaluate_on_grid(
        grid,
        states,
        viscosity=SMOKE_VISCOSITY,
        interior_margin=SMOKE_MARGIN,
        pressure_recovery=NO_PRESSURE_RECOVERY,
    )
    return checkpoint_dir, grid, states, baseline


def _relative(evaluation: dict[str, Any], name: str) -> float:
    value = evaluation["regions"]["interior"][name]["relative_rms"]
    assert value is not None
    return float(value)


def _order(evaluation: dict[str, Any]) -> float | None:
    return evaluation["time_difference_richardson"]["observed_order"]


def test_unperturbed_real_states_pass_the_same_evaluation(fault_baseline) -> None:
    _, _, _, baseline = fault_baseline
    assert _relative(baseline, "curl_momentum_residual") < SMOKE_CURL_RESIDUAL_BOUND
    assert _relative(baseline, "vorticity_transport_residual") < (
        SMOKE_TRANSPORT_RESIDUAL_BOUND
    )
    order = _order(baseline)
    assert order is not None and 1.8 <= order <= 2.2


def test_detects_odd_in_r_u1_corruption_of_a_saved_neighbour_state(
    fault_baseline,
) -> None:
    r"""A tiny odd-in-``r`` ``u1`` kink in one saved state must be rejected.

    ``u1`` is even in ``r`` for every admissible state (E-16c), so
    ``u1 -> u1 + c*r*cos(2 pi z)`` is a genuine violation and not an arbitrary
    edit.  It is injected into the ``+dt`` neighbour, i.e. into a state that
    enters the evaluation only through the centered time difference, where it
    is divided by ``2 dt``.  A corruption of one part in a thousand of
    ``max|u1|`` is therefore amplified enormously, which is exactly the
    sensitivity this experiment claims for its controlled ``u_t``.
    """

    checkpoint_dir, grid, states, baseline = fault_baseline
    path = checkpoint_dir / state_filename(1)
    candidate = load_candidate(path)
    radius = candidate.grid.r[:, None]
    axial = candidate.grid.z[None, :]
    amplitude = 1.0e-3 * float(np.max(np.abs(candidate.fields["u1"])))
    corrupted = dataclasses.replace(
        candidate,
        fields={
            **candidate.fields,
            "u1": candidate.fields["u1"]
            + amplitude * radius * np.cos(2.0 * np.pi * axial),
        },
    )
    arrays = reconstruct_loaded_candidate_on_cartesian(corrupted, grid)
    perturbed = dict(states)
    perturbed[1] = dataclasses.replace(
        states[1],
        velocity=arrays.velocity,
        vorticity_e18b=arrays.vorticity_e18b,
    )
    metrics = evaluate_on_grid(
        grid,
        perturbed,
        viscosity=SMOKE_VISCOSITY,
        interior_margin=SMOKE_MARGIN,
        pressure_recovery=NO_PRESSURE_RECOVERY,
    )
    ratio = _relative(metrics, "curl_momentum_residual") / _relative(
        baseline, "curl_momentum_residual"
    )
    assert ratio > 10.0, ratio
    # The Richardson check sees it too: the dt-spaced difference no longer
    # belongs to the same expansion as the 2dt- and 4dt-spaced ones.
    order = _order(metrics)
    assert order is None or not (1.8 <= order <= 2.2), order


def test_detects_a_wrong_dt_in_the_time_difference(fault_baseline) -> None:
    """Doubling the assumed spacing halves ``u_t`` and breaks both checks."""

    _, grid, states, baseline = fault_baseline
    centre = states[0].physical_time
    step = states[1].physical_time - centre
    mislabelled = dict(states)
    mislabelled[1] = dataclasses.replace(
        states[1], physical_time=centre + 2.0 * step
    )
    mislabelled[-1] = dataclasses.replace(
        states[-1], physical_time=centre - 2.0 * step
    )
    metrics = evaluate_on_grid(
        grid,
        mislabelled,
        viscosity=SMOKE_VISCOSITY,
        interior_margin=SMOKE_MARGIN,
        pressure_recovery=NO_PRESSURE_RECOVERY,
    )
    curl_ratio = _relative(metrics, "curl_momentum_residual") / _relative(
        baseline, "curl_momentum_residual"
    )
    transport_ratio = _relative(metrics, "vorticity_transport_residual") / _relative(
        baseline, "vorticity_transport_residual"
    )
    assert curl_ratio > 3.0, curl_ratio
    assert transport_ratio > 2.0, transport_ratio
    order = _order(metrics)
    assert order is None or not (1.8 <= order <= 2.2), order


def test_viscous_sign_is_not_detectable_on_the_real_states(fault_baseline) -> None:
    """Recorded negative result, not a hidden one.

    At ``nu = 5e-4`` the viscous term is about ``1e-4`` of the time-derivative
    and advection terms, so flipping its sign moves the residual by far less
    than the audit truncation.  This asserts the measurement that justifies
    the corresponding entry in :data:`run_hou_primitive_residual.KNOWN_GAPS`,
    and the decisive version of the same injection is
    :func:`test_detects_flipped_viscous_sign_on_a_manufactured_field`.
    """

    _, grid, states, baseline = fault_baseline
    flipped = evaluate_on_grid(
        grid,
        states,
        viscosity=SMOKE_VISCOSITY,
        interior_margin=SMOKE_MARGIN,
        pressure_recovery=NO_PRESSURE_RECOVERY,
        viscous_sign=-1.0,
    )
    ratio = _relative(flipped, "curl_momentum_residual") / _relative(
        baseline, "curl_momentum_residual"
    )
    assert ratio < 1.1, ratio
    fractions = baseline["regions"]["interior"]["term_magnitudes"]["curl_momentum"][
        "fraction_of_term_rms_sum"
    ]
    assert fractions["curl_viscous"] < 1.0e-2, fractions
    gaps = " ".join(experiment.KNOWN_GAPS)
    assert "viscous term" in gaps


# --------------------------------------------------------------------------
# manufactured verification of the checker itself
# --------------------------------------------------------------------------

_MANUFACTURED_VISCOSITY = 0.03


def _manufactured_grid(count: int) -> UniformCartesianGrid:
    """Non-periodic in x and y, periodic in z: the audit box topology."""

    return UniformCartesianGrid.uniform(
        shape=(count, count, max(8, count // 2)),
        x_bounds=(-1.3, 1.1),
        y_bounds=(-0.9, 1.5),
        z_bounds=(0.0, 2.0 * math.pi),
        periodic=(False, False, True),
    )


def _taylor_green(
    grid: UniformCartesianGrid,
    time: float,
    *,
    viscosity: float = _MANUFACTURED_VISCOSITY,
) -> dict[str, np.ndarray]:
    r"""The decaying two-dimensional Taylor-Green vortex, verified in place.

    With ``F(t) = exp(-2 nu t)``,

    ``u = (F sin x cos y, -F cos x sin y, 0)``,
    ``p = (F^2/4)(cos 2x + cos 2y)``.

    ``div u = 0``; ``Laplacian u = -2 u`` so ``u_t - nu Laplacian u = 0``;
    ``(u.grad)u = ((F^2/2) sin 2x, (F^2/2) sin 2y, 0) = -grad p``.  Hence this
    solves the unforced incompressible Navier-Stokes equations exactly, and
    :func:`test_manufactured_taylor_green_really_solves_navier_stokes` checks
    every one of those claims numerically rather than trusting the algebra.
    """

    x, y, _ = grid.mesh()
    decay = math.exp(-2.0 * float(viscosity) * float(time))
    velocity = np.stack(
        (
            decay * np.sin(x) * np.cos(y),
            -decay * np.cos(x) * np.sin(y),
            np.zeros_like(x),
        ),
        axis=0,
    )
    return {
        "velocity": velocity,
        "velocity_t": -2.0 * float(viscosity) * velocity,
        "laplacian": -2.0 * velocity,
        "advection": np.stack(
            (
                0.5 * decay**2 * np.sin(2.0 * x),
                0.5 * decay**2 * np.sin(2.0 * y),
                np.zeros_like(x),
            ),
            axis=0,
        ),
        "pressure": 0.25 * decay**2 * (np.cos(2.0 * x) + np.cos(2.0 * y)),
        "pressure_gradient": np.stack(
            (
                -0.5 * decay**2 * np.sin(2.0 * x),
                -0.5 * decay**2 * np.sin(2.0 * y),
                np.zeros_like(x),
            ),
            axis=0,
        ),
    }


def _solenoidal_field(grid: UniformCartesianGrid) -> np.ndarray:
    r"""``curl`` of a smooth z-periodic vector potential, in closed form.

    Unlike the Taylor-Green vortex this field's advection is *not* a gradient,
    so both nonlinear contributions are nonzero and the two pressure-free
    forms are genuinely different.
    """

    x, y, z = grid.mesh()
    return np.stack(
        (
            0.5 * np.cos(x) * np.cos(y) * np.cos(z) - 0.7 * np.sin(x) * np.cos(z),
            -np.sin(2.0 * y) * np.sin(z) + 0.5 * np.sin(x) * np.sin(y) * np.cos(z),
            0.7 * np.cos(x) * np.sin(z) - 2.0 * np.cos(2.0 * y) * np.cos(z),
        ),
        axis=0,
    )


def _order_of(values: list[float], spacings: list[float]) -> list[float]:
    return [
        math.log(values[index] / values[index + 1])
        / math.log(spacings[index] / spacings[index + 1])
        for index in range(len(values) - 1)
    ]


def test_gradient_adjoint_is_the_exact_transpose_of_the_audit_gradient() -> None:
    """``<D u, q> == <u, D^T q>`` to machine precision on every axis."""

    generator = np.random.default_rng(20260728)
    grid = _manufactured_grid(12)
    for axis in range(3):
        left = generator.standard_normal(grid.shape)
        right = generator.standard_normal(grid.shape)
        forward = float(np.sum(cartesian_gradient(grid, left)[axis] * right))
        backward = float(
            np.sum(left * adjoint_first_derivative(grid, right, axis))
        )
        assert forward == pytest.approx(backward, rel=1.0e-12, abs=1.0e-12)


def test_gradient_null_space_is_annihilated_by_the_audit_gradient() -> None:
    grid = _manufactured_grid(12)
    basis = gradient_null_space(grid)
    # non-periodic x and y contribute constants only; periodic z with an even
    # point count adds the alternating mode.
    assert len(basis) == 2
    for vector in basis:
        assert float(np.max(np.abs(cartesian_gradient(grid, vector)))) < 1.0e-14
        assert float(np.sum(vector * vector)) == pytest.approx(1.0)


def test_manufactured_taylor_green_really_solves_navier_stokes() -> None:
    """Verify the exact-solution claim before any test relies on it."""

    grid = _manufactured_grid(64)
    fields = _taylor_green(grid, 0.3)
    window = region_window(grid, 4)
    scalar_window = window
    vector_window = (slice(None), *window)
    velocity = fields["velocity"]
    scale = float(np.max(np.abs(velocity)))

    assert rms(cartesian_divergence(grid, velocity)[scalar_window]) < 1.0e-2 * scale
    residual = checker.momentum_terms(
        grid,
        velocity=velocity,
        velocity_t=fields["velocity_t"],
        viscosity=_MANUFACTURED_VISCOSITY,
    )
    # The three analytic identities the closed form rests on.
    assert rms(
        (residual.advection - fields["advection"])[vector_window]
    ) < 1.0e-2 * rms(fields["advection"][vector_window])
    analytic_viscous = -_MANUFACTURED_VISCOSITY * fields["laplacian"]
    assert rms(
        (residual.viscous - analytic_viscous)[vector_window]
    ) < 1.0e-2 * rms(analytic_viscous[vector_window])
    analytic_gradient = cartesian_gradient(grid, fields["pressure"])
    assert rms(
        (analytic_gradient - fields["pressure_gradient"])[vector_window]
    ) < 1.0e-2 * rms(fields["pressure_gradient"][vector_window])
    # R = -grad p, the statement the whole experiment rests on.
    assert rms(
        (residual.residual + fields["pressure_gradient"])[vector_window]
    ) < 1.0e-2 * rms(fields["pressure_gradient"][vector_window])


def test_pressure_free_residuals_are_at_discretization_level_on_the_solution() -> None:
    """Both forms vanish to truncation on the exact Taylor-Green solution."""

    relative_curl: list[float] = []
    relative_transport: list[float] = []
    spacings: list[float] = []
    for count in (24, 32, 48):
        grid = _manufactured_grid(count)
        fields = _taylor_green(grid, 0.3)
        window = region_window(grid, 4)
        selector = (slice(None), *window)
        curl_form = curl_momentum_residual(
            grid,
            velocity=fields["velocity"],
            velocity_t=fields["velocity_t"],
            viscosity=_MANUFACTURED_VISCOSITY,
        )
        vorticity = cartesian_curl(grid, fields["velocity"])
        transport_form = vorticity_transport_residual(
            grid,
            velocity=fields["velocity"],
            vorticity=vorticity,
            vorticity_t=cartesian_curl(grid, fields["velocity_t"]),
            viscosity=_MANUFACTURED_VISCOSITY,
        )
        spacings.append(grid.spacings[0])
        relative_curl.append(
            rms(curl_form.residual[selector])
            / float(np.max(np.abs(curl_form.term_sum[selector])))
        )
        relative_transport.append(
            rms(transport_form.residual[selector])
            / float(np.max(np.abs(transport_form.term_sum[selector])))
        )
    assert max(relative_curl) < 1.0e-3, relative_curl
    assert max(relative_transport) < 1.0e-3, relative_transport
    assert relative_curl[0] > relative_curl[1] > relative_curl[2]
    assert relative_transport[0] > relative_transport[1] > relative_transport[2]
    assert min(_order_of(relative_curl, spacings)) > 1.3, (
        _order_of(relative_curl, spacings)
    )


def test_taylor_green_does_not_exercise_the_nonlinear_terms() -> None:
    """Documented degeneracy: both nonlinear contributions vanish there.

    ``(u.grad)u`` is an exact gradient for this field, so ``curl`` of it and
    ``(u.grad)omega - (omega.grad)u`` are both identically zero.  A test that
    only used Taylor-Green would therefore say nothing about the nonlinear
    discretizations, which is why
    :func:`test_two_pressure_free_forms_agree_at_second_order` exists.
    """

    grid = _manufactured_grid(24)
    fields = _taylor_green(grid, 0.3)
    window = region_window(grid, 4)
    selector = (slice(None), *window)
    advection = checker.cartesian_directional_derivative(
        grid, advecting=fields["velocity"], field=fields["velocity"]
    )
    vorticity = cartesian_curl(grid, fields["velocity"])
    scale = rms(advection[selector])
    assert scale > 0.1
    assert rms(cartesian_curl(grid, advection)[selector]) < 1.0e-12 * scale
    stretched = checker.cartesian_directional_derivative(
        grid, advecting=vorticity, field=fields["velocity"]
    )
    advected = checker.cartesian_directional_derivative(
        grid, advecting=fields["velocity"], field=vorticity
    )
    assert rms((advected - stretched)[selector]) < 1.0e-12 * scale
    assert "Taylor-Green" in " ".join(experiment.KNOWN_GAPS)


def test_two_pressure_free_forms_agree_at_second_order() -> None:
    """On a nontrivial solenoidal field the two forms differ, and converge."""

    differences: list[float] = []
    spacings: list[float] = []
    magnitudes: list[float] = []
    for count in (24, 32, 48):
        grid = _manufactured_grid(count)
        velocity = _solenoidal_field(grid)
        velocity_t = 0.4 * velocity
        window = region_window(grid, 4)
        selector = (slice(None), *window)
        curl_form = curl_momentum_residual(
            grid,
            velocity=velocity,
            velocity_t=velocity_t,
            viscosity=_MANUFACTURED_VISCOSITY,
        )
        vorticity = cartesian_curl(grid, velocity)
        transport_form = vorticity_transport_residual(
            grid,
            velocity=velocity,
            vorticity=vorticity,
            vorticity_t=0.4 * vorticity,
            viscosity=_MANUFACTURED_VISCOSITY,
        )
        spacings.append(grid.spacings[0])
        magnitudes.append(rms(curl_form.residual[selector]))
        differences.append(
            rms((curl_form.residual - transport_form.residual)[selector])
        )
    # They are genuinely different discretizations at every resolution ...
    assert min(differences) > 0.0
    assert differences[0] / magnitudes[0] > 1.0e-2
    # ... and their difference converges away at second order.
    orders = _order_of(differences, spacings)
    assert min(orders) > 1.7, (orders, differences)


def test_pressure_recovery_converges_at_second_order() -> None:
    """The recovered pressure gradient converges to the closed-form one."""

    errors: list[float] = []
    consistency: list[float] = []
    spacings: list[float] = []
    for count in (24, 32, 48):
        grid = _manufactured_grid(count)
        fields = _taylor_green(grid, 0.3)
        window = region_window(grid, 4)
        selector = (slice(None), *window)
        terms = checker.momentum_terms(
            grid,
            velocity=fields["velocity"],
            velocity_t=fields["velocity_t"],
            viscosity=_MANUFACTURED_VISCOSITY,
        )
        recovered = recover_pressure(
            grid,
            terms.residual,
            momentum_term_sum=terms.term_sum,
            tolerance=1.0e-12,
            max_iterations=20000,
        )
        assert recovered.converged
        spacings.append(grid.spacings[0])
        errors.append(
            rms(
                (recovered.pressure_gradient - fields["pressure_gradient"])[selector]
            )
        )
        consistency.append(rms(recovered.consistency[selector]))
    assert errors[0] > errors[1] > errors[2] > 0.0
    # Measured on this ladder: gradient orders (1.914, 1.954), consistency
    # orders (2.665, 2.627).  The bound pins second order, not the build.
    orders = _order_of(errors, spacings)
    assert min(orders) > 1.8, (orders, errors)
    assert min(_order_of(consistency, spacings)) > 1.8, (consistency, spacings)


def test_detects_flipped_viscous_sign_on_a_manufactured_field() -> None:
    """With a visible viscous term the sign flip is caught decisively."""

    grid = _manufactured_grid(32)
    viscosity = 0.5
    fields = _taylor_green(grid, 0.3, viscosity=viscosity)
    window = region_window(grid, 4)
    selector = (slice(None), *window)

    def relative(sign: float) -> float:
        form = curl_momentum_residual(
            grid,
            velocity=fields["velocity"],
            velocity_t=fields["velocity_t"],
            viscosity=viscosity,
            viscous_sign=sign,
        )
        return rms(form.residual[selector]) / float(
            np.max(np.abs(form.term_sum[selector]))
        )

    baseline = relative(1.0)
    flipped = relative(-1.0)
    assert flipped / baseline > 10.0, (baseline, flipped)


def test_recover_pressure_rejects_a_degenerate_request() -> None:
    grid = _manufactured_grid(8)
    with pytest.raises(ValueError, match="tolerance"):
        recover_pressure(grid, np.zeros((3, *grid.shape)), tolerance=0.0)
    with pytest.raises(ValueError, match="max_iterations"):
        recover_pressure(grid, np.zeros((3, *grid.shape)), max_iterations=0)
    # A residual that is exactly zero has no pressure to recover; the routine
    # must return the zero field rather than divide by a vanishing norm.
    recovered = recover_pressure(grid, np.zeros((3, *grid.shape)))
    assert float(np.max(np.abs(recovered.pressure))) == 0.0
    assert recovered.converged is True


def test_region_window_rejects_an_empty_interior() -> None:
    grid = _manufactured_grid(8)
    with pytest.raises(ValueError, match="empty norm region"):
        region_window(grid, 4)
    with pytest.raises(ValueError, match="nonnegative integer"):
        region_window(grid, -1)


def test_known_gaps_and_limitations_are_declared_not_dropped() -> None:
    text = " ".join(experiment.KNOWN_GAPS) + " " + " ".join(experiment.LIMITATIONS)
    assert "scoped" in text.lower()
    assert "bilinear" in text.lower()
    assert "shipped outputs/hou_early_time_v1" in text
    assert copy.deepcopy(experiment.KNOWN_GAPS) == experiment.KNOWN_GAPS
