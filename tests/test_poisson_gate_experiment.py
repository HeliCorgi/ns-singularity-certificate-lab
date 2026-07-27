"""Tests for the finite-cylinder Poisson gate experiment runner.

The gate writes a preregistered acceptance verdict, so this module exercises
the real :func:`experiments.run_poisson_gate.run` end to end (at reduced
resolutions, to stay CPU-friendly), and pins the config validation, the
refusal to overwrite existing evidence, and the runtime provenance record.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
from typing import Any, Iterator
import uuid

import pytest

from experiments.run_poisson_gate import (
    EXPECTED_EXPERIMENT,
    _validate_config,
    run as run_gate,
)
from ns_certificate_lab.provenance import (
    PROVENANCE_SCHEMA,
    validate_runtime_provenance,
)

REPOSITORY = Path(__file__).resolve().parents[1]
SCRATCH_ROOT = (Path(__file__).resolve().parent / "_scratch").resolve()

# Reduced refinement ladder: the shipped config uses 17/33/65, which is the
# same three-grid study one octave finer.  Three grids are the documented
# minimum, so this still exercises the two-order acceptance path.
REDUCED_RESOLUTIONS = [
    {"nr": 9, "nz": 16},
    {"nr": 17, "nz": 32},
    {"nr": 33, "nz": 64},
]


def _reduced_config() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "experiment": EXPECTED_EXPERIMENT,
        "r_max": 1.7,
        "z_min": 0.0,
        "z_max": 6.283185307179586,
        "mode": 3,
        "resolutions": copy.deepcopy(REDUCED_RESOLUTIONS),
        "condition_mode_indices": [0, 1, 3],
    }


def _make_scratch_directory() -> Path:
    destination = (SCRATCH_ROOT / uuid.uuid4().hex).resolve()
    if SCRATCH_ROOT not in destination.parents:
        raise RuntimeError("refusing to create a test directory outside tests/_scratch")
    destination.mkdir(parents=True)
    return destination


def _remove_scratch_directory(destination: Path) -> None:
    if destination.exists():
        resolved = destination.resolve()
        if SCRATCH_ROOT not in resolved.parents:
            raise RuntimeError("refusing to remove a directory outside tests/_scratch")
        shutil.rmtree(resolved)
    if SCRATCH_ROOT.exists() and not any(SCRATCH_ROOT.iterdir()):
        SCRATCH_ROOT.rmdir()


@pytest.fixture(scope="module")
def module_scratch_dir() -> Iterator[Path]:
    destination = _make_scratch_directory()
    try:
        yield destination
    finally:
        _remove_scratch_directory(destination)


def write_config(directory: Path, config: dict[str, Any], *, name: str) -> Path:
    path = directory / name
    path.write_text(
        json.dumps(config, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture(scope="module")
def gate_run(module_scratch_dir: Path) -> tuple[dict[str, Any], Path]:
    """Run the real gate once at reduced resolutions."""

    config_path = write_config(
        module_scratch_dir,
        _reduced_config(),
        name="reduced_gate.json",
    )
    output_dir = module_scratch_dir / "gate-evidence"
    summary = run_gate(config_path, output_dir)
    return summary, output_dir


def test_gate_reports_full_acceptance(gate_run: tuple[dict[str, Any], Path]) -> None:
    summary, output_dir = gate_run
    acceptance = summary["acceptance"]
    assert acceptance["all_passed"] is True
    assert all(bool(value) for value in acceptance.values())

    rows = summary["rows"]
    assert [row["nr"] for row in rows] == [item["nr"] for item in REDUCED_RESOLUTIONS]
    assert rows[0]["observed_order_from_previous"] is None
    orders = [float(row["observed_order_from_previous"]) for row in rows[1:]]
    assert len(orders) == 2
    assert min(orders) >= 1.85

    for name in (
        "summary.json",
        "convergence.csv",
        "finest_solution.npz",
        "config.snapshot.json",
        "manifest.json",
        "manifest.json.sha256",
    ):
        assert (output_dir / name).is_file()

    on_disk = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
    assert on_disk["acceptance"] == acceptance
    assert on_disk["experiment"] == EXPECTED_EXPERIMENT


def test_gate_summary_records_runtime_provenance(
    gate_run: tuple[dict[str, Any], Path],
) -> None:
    _, output_dir = gate_run
    on_disk = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
    provenance = on_disk["reproducibility"]["runtime_provenance"]
    assert provenance["schema"] == PROVENANCE_SCHEMA
    assert provenance["schema"] == "ns-certificate-lab/runtime-provenance/v2"
    # Strict schema validation, so a truncated or hand-edited record fails.
    validate_runtime_provenance(provenance)
    assert len(provenance["source_fingerprint_sha256"]) == 64
    assert "outputs" not in provenance["source_fingerprint_scope"]


def test_gate_refuses_to_overwrite_a_nonempty_output_directory(
    module_scratch_dir: Path,
) -> None:
    config_path = write_config(
        module_scratch_dir,
        _reduced_config(),
        name="refusal_gate.json",
    )
    output_dir = module_scratch_dir / "existing-evidence"
    output_dir.mkdir()
    marker = output_dir / "prior-evidence.txt"
    marker.write_text("preserve\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run_gate(config_path, output_dir)

    assert marker.read_text(encoding="utf-8") == "preserve\n"
    assert set(path.name for path in output_dir.iterdir()) == {"prior-evidence.txt"}


def test_shipped_gate_config_passes_validation() -> None:
    config = json.loads(
        (REPOSITORY / "configs" / "poisson_gate.json").read_text(encoding="utf-8")
    )
    _validate_config(config)
    assert config["experiment"] == EXPECTED_EXPERIMENT
    assert len(config["resolutions"]) >= 3


@pytest.mark.parametrize(
    ("field", "wrong_value"),
    [
        ("experiment", "some_other_experiment"),
        ("experiment", "finite_cylinder_poisson_gate_v2"),
        ("schema_version", 2),
        ("schema_version", "1"),
    ],
)
def test_detects_mislabeled_config(
    module_scratch_dir: Path,
    field: str,
    wrong_value: Any,
) -> None:
    # Injected defect: the config claims to be a different experiment or a
    # different schema than the one this runner implements.  Writing an
    # acceptance verdict under the wrong label would corrupt the evidence
    # trail, so it must be refused before any output directory is created.
    config = _reduced_config()
    config[field] = wrong_value
    with pytest.raises(ValueError):
        _validate_config(config)

    label = f"mislabeled-{field}-{str(wrong_value).replace(' ', '_')}"
    config_path = write_config(module_scratch_dir, config, name=f"{label}.json")
    output_dir = module_scratch_dir / label
    with pytest.raises(ValueError):
        run_gate(config_path, output_dir)
    assert not output_dir.exists()


@pytest.mark.parametrize(
    "mutate",
    [
        pytest.param(lambda config: config.pop("mode"), id="missing-key"),
        pytest.param(
            lambda config: config.update({"unexpected_knob": 1}),
            id="unknown-key",
        ),
        pytest.param(
            lambda config: config.update({"r_max": -1.0}),
            id="nonpositive-r-max",
        ),
        pytest.param(
            lambda config: config.update({"z_max": config["z_min"]}),
            id="empty-period",
        ),
        pytest.param(
            lambda config: config.update({"mode": 0}),
            id="nonpositive-mode",
        ),
        pytest.param(
            lambda config: config.update(
                {"resolutions": [{"nr": 17, "nz": 32, "extra": 1}]}
            ),
            id="malformed-resolution",
        ),
        pytest.param(
            lambda config: config.update(
                {
                    "resolutions": [
                        {"nr": 33, "nz": 64},
                        {"nr": 17, "nz": 32},
                        {"nr": 65, "nz": 128},
                    ]
                }
            ),
            id="unsorted-resolutions",
        ),
        pytest.param(
            lambda config: config.update({"condition_mode_indices": [0, -1]}),
            id="negative-condition-index",
        ),
    ],
)
def test_detects_invalid_gate_config(mutate: Any) -> None:
    # Each case injects one config defect that would otherwise be silently
    # accepted or silently ignored by the runner.
    config = _reduced_config()
    mutate(config)
    with pytest.raises(ValueError):
        _validate_config(config)


def test_detects_too_few_resolutions() -> None:
    # Injected defect: a two-grid study.  Two grids yield only one error
    # ratio, so the reported "observed order" could not be corroborated by a
    # second refinement; the gate must refuse rather than report it.
    config = _reduced_config()
    config["resolutions"] = copy.deepcopy(REDUCED_RESOLUTIONS[:2])
    with pytest.raises(ValueError, match="at least three"):
        _validate_config(config)
