from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.run_nonlinear_control import (
    evaluate,
    evolve_resolution,
    observed_orders,
    run,
)
from ns_certificate_lab._integrity import strict_json_loads


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "configs" / "nonlinear_manufactured_control.json"


def _config() -> dict[str, object]:
    value = strict_json_loads(
        CONFIG_PATH.read_text(encoding="utf-8"),
        label="nonlinear control test config",
    )
    assert isinstance(value, dict)
    return value


@pytest.fixture(scope="module")
def generated_control(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[dict[str, object], Path, dict[str, object]]:
    config = _config()
    output = tmp_path_factory.mktemp("nonlinear-control") / "evidence"
    summary = run(config, output)
    return config, output, summary


def test_forced_nonlinear_poisson_control_converges(
    generated_control: tuple[dict[str, object], Path, dict[str, object]],
) -> None:
    config, _, summary = generated_control
    assert summary["accepted_as_forced_smooth_control"] is True
    assert all(summary["acceptance_checks"].values())
    assert "Forced smooth manufactured negative control" in summary["interpretation"]
    assert "unforced general" in summary["interpretation"]
    resolutions = summary["resolutions"]
    assert [(item["nr"], item["nz"]) for item in resolutions] == list(
        zip(config["radial_resolutions"], config["axial_resolutions"])
    )
    assert all(
        item["actual_dt"]
        <= float(config["time_step_coefficient"]) * item["dr"] ** 2
        for item in resolutions
    )
    for name in ("u1", "omega1", "psi1", "velocity"):
        errors = [item["relative_errors"][name] for item in resolutions]
        assert errors[0] > errors[1] > errors[2]
        # Recompute orders without using the experiment helper.
        independent_orders = [
            math.log(errors[index] / errors[index + 1])
            / math.log(
                resolutions[index]["dr"] / resolutions[index + 1]["dr"]
            )
            for index in range(2)
        ]
        assert independent_orders == pytest.approx(summary["observed_orders"][name])
        assert min(independent_orders) > float(
            config["acceptance"]["minimum_observed_order"]
        )
    assert max(
        item["maximum_rhs_discrete_poisson_residual"] for item in resolutions
    ) < float(config["acceptance"]["maximum_discrete_poisson_residual"])
    assert summary["finest_energy_relative_error"] < float(
        config["acceptance"]["maximum_energy_relative_error"]
    )
    assert summary["finest_max_vorticity_relative_error"] < float(
        config["acceptance"]["maximum_vorticity_relative_error"]
    )


def test_poisson_solver_is_called_for_every_heun_rhs(monkeypatch: pytest.MonkeyPatch) -> None:
    import experiments.run_nonlinear_control as experiment

    config = _config()
    original = experiment.solve_streamfunction_poisson
    calls = 0

    def counted(*args: object, **kwargs: object) -> object:
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(experiment, "solve_streamfunction_poisson", counted)
    result = experiment.evolve_resolution(config, nr=9, nz=16)
    # Two RK stages per time step, plus diagnostic and final-audit solves.
    assert calls >= 2 * result.step_count + 2


@pytest.mark.parametrize(
    ("fault", "diagnostic"),
    [
        ("poisson_sign", "psi1"),
        ("velocity_z_sign", "velocity"),
        ("outer_trace", "boundary"),
    ],
)
def test_end_to_end_faults_are_rejected(
    fault: str,
    diagnostic: str,
) -> None:
    config = _config()
    nr = int(config["radial_resolutions"][-1])
    nz = int(config["axial_resolutions"][-1])
    normal = evolve_resolution(config, nr=nr, nz=nz)
    broken = evolve_resolution(config, nr=nr, nz=nz, fault=fault)
    if diagnostic == "boundary":
        error = float(
            broken.poisson_metadata["analytic_outer_trace_max_abs_error"]
        )
        assert error > float(config["acceptance"]["maximum_outer_trace_error"])
    else:
        error = broken.relative_errors[diagnostic]
        assert error > float(config["acceptance"]["maximum_finest_relative_error"])
        assert error > 10.0 * normal.relative_errors[diagnostic]


def test_artifacts_are_complete_finite_and_checksummed(
    generated_control: tuple[dict[str, object], Path, dict[str, object]],
) -> None:
    config, output, _ = generated_control
    summary = strict_json_loads(
        (output / "summary.json").read_text(encoding="utf-8"),
        label="generated nonlinear control summary",
    )
    assert summary["reproducibility"]["seed"] == config["seed"]
    assert summary["reproducibility"]["runtime_provenance"][
        "source_fingerprint_sha256"
    ]

    with (output / "convergence.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3
    required = {
        "u1_relative_l2_error",
        "omega1_relative_l2_error",
        "psi1_relative_l2_error",
        "velocity_relative_l2_error",
        "final_energy",
        "maximum_vorticity_over_samples",
        "maximum_rhs_discrete_poisson_residual",
        "maximum_rhs_independent_pde_poisson_residual",
        "analytic_psi_outer_trace_error",
    }
    assert required <= rows[0].keys()
    assert all(
        math.isfinite(float(row[key]))
        for row in rows
        for key in required
    )

    with (output / "diagnostics.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        diagnostics = list(csv.DictReader(stream))
    assert diagnostics
    assert all(
        math.isfinite(float(value))
        for row in diagnostics
        for value in row.values()
    )

    with np.load(output / "trajectories.npz", allow_pickle=False) as archive:
        assert archive.files
        assert any(name.startswith("psi1_final_") for name in archive.files)
        assert all(np.all(np.isfinite(archive[name])) for name in archive.files)

    manifest = strict_json_loads(
        (output / "manifest.json").read_text(encoding="utf-8"),
        label="generated nonlinear control manifest",
    )
    expected_files = {
        "config.snapshot.json",
        "convergence.csv",
        "diagnostics.csv",
        "summary.json",
        "trajectories.npz",
    }
    assert set(manifest["files"]) == expected_files
    for name, record in manifest["files"].items():
        payload = (output / name).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == record["sha256"]
        assert len(payload) == record["bytes"]
    assert (
        output / "manifest.json.sha256"
    ).read_text(encoding="ascii").strip() == hashlib.sha256(
        (output / "manifest.json").read_bytes()
    ).hexdigest()


def test_nonempty_output_and_invalid_config_are_rejected(
    scratch_dir: Path,
) -> None:
    config = _config()
    occupied = scratch_dir / "occupied"
    occupied.mkdir()
    marker = occupied / "preserve.txt"
    marker.write_text("keep\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run(config, occupied)
    assert marker.read_text(encoding="utf-8") == "keep\n"

    invalid = json.loads(json.dumps(config))
    invalid["radial_resolutions"] = [9, 17]
    with pytest.raises(ValueError, match="exactly three"):
        evaluate(invalid)


def test_observed_order_helper_rejects_bad_refinement() -> None:
    assert observed_orders([4.0, 1.0, 0.25], [0.4, 0.2, 0.1]) == pytest.approx(
        [2.0, 2.0]
    )
    with pytest.raises(ValueError, match="strictly decreasing"):
        observed_orders([1.0, 0.5], [0.1, 0.2])
