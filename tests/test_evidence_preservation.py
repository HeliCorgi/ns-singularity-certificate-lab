from __future__ import annotations

import json
from pathlib import Path

import pytest

from experiments.run_baseline import run as run_baseline
from experiments.run_manufactured import run as run_manufactured


def test_experiments_refuse_to_overwrite_existing_evidence(
    scratch_dir: Path,
) -> None:
    output = scratch_dir / "existing"
    output.mkdir()
    marker = output / "prior_result.txt"
    marker.write_text("preserve me\n", encoding="utf-8")

    baseline_config = json.loads(
        (Path(__file__).resolve().parents[1] / "configs" / "baseline.json").read_text(
            encoding="utf-8"
        )
    )
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run_baseline(baseline_config, output)
    assert marker.read_text(encoding="utf-8") == "preserve me\n"

    config_path = (
        Path(__file__).resolve().parents[1] / "configs" / "manufactured.json"
    )
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run_manufactured(config_path, output)
    assert marker.read_text(encoding="utf-8") == "preserve me\n"
