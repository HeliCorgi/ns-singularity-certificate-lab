"""Run the renormalized Galerkin cascade discovery pilot.

Forward-integrates the lattice doubling renormalization map (evolve one
parabolic stage, exact even-mode pullback ``w(k) = 4 v(2k)``, projective
energy renormalization) over a scan of critical energy constants, looking for
an attracting projective orbit with a positive gain plateau.  Binary64
discovery diagnostic only: no interval enclosure, no PDE statement, and a
decaying gain here rejects only the scanned parameter box, not the continuum
operator.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    strict_json_loads,
    write_with_digest,
)
from ns_certificate_lab.renormalized_cascade import (
    CascadeConfig,
    run_renormalized_cascade,
)

CONFIG_SCHEMA = "ns-certificate-lab/renormalized-cascade-config/v1"
OUTPUT_SCHEMA = "ns-certificate-lab/renormalized-cascade/v1"
STATUS = "BINARY64 DISCOVERY DIAGNOSTIC / NOT A PROOF"
PLATEAU_WINDOW = 4
PLATEAU_SPREAD = 0.10


def _validate_config(value: object) -> dict[str, Any]:
    required = {
        "schema",
        "scale",
        "width",
        "grid_size",
        "viscosity",
        "tau",
        "base_steps",
        "max_steps",
        "stages",
        "cfl_safety",
        "energy_constants",
        "drop_below_options",
    }
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("cascade config has missing or unknown fields")
    if value["schema"] != CONFIG_SCHEMA:
        raise ValueError("cascade config schema is invalid")
    constants = value["energy_constants"]
    if not isinstance(constants, list) or not constants:
        raise ValueError("energy_constants must be a nonempty list")
    if any(
        isinstance(item, bool)
        or not isinstance(item, (int, float))
        or not math.isfinite(float(item))
        or float(item) <= 0.0
        for item in constants
    ):
        raise ValueError("energy_constants must be finite positive numbers")
    drops = value["drop_below_options"]
    if not isinstance(drops, list) or not drops:
        raise ValueError("drop_below_options must be a nonempty list")
    for item in drops:
        if item is None:
            continue
        if isinstance(item, bool) or not isinstance(item, int) or item < 1:
            raise ValueError(
                "drop_below_options items must be null or positive integers"
            )
    return value


def _classify(gains: tuple[float, ...], terminated_reason: str) -> str:
    if terminated_reason in {"float_noise_floor_reached", "pullback_empty"}:
        return "collapsed"
    if not gains:
        return "no_stages"
    if len(gains) < PLATEAU_WINDOW + 1:
        return "too_short"
    window = gains[-PLATEAU_WINDOW:]
    if any(gain <= 0.0 for gain in window):
        return "collapsed"
    logs = [math.log(gain) for gain in window]
    spread = max(logs) - min(logs)
    slope = logs[-1] - logs[0]
    if spread <= PLATEAU_SPREAD:
        return "plateau"
    if slope < 0.0:
        return "decaying"
    return "growing"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()

    config = _validate_config(
        strict_json_loads(
            arguments.config.read_text(encoding="utf-8"),
            label="cascade config",
        )
    )
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    runs: list[dict[str, Any]] = []
    for energy_constant in config["energy_constants"]:
        for drop_below in config["drop_below_options"]:
            cascade_config = CascadeConfig(
                scale=config["scale"],
                width=config["width"],
                grid_size=config["grid_size"],
                viscosity=float(config["viscosity"]),
                energy_constant=float(energy_constant),
                tau=float(config["tau"]),
                base_steps=config["base_steps"],
                max_steps=config["max_steps"],
                stages=config["stages"],
                drop_below=drop_below,
                cfl_safety=float(config["cfl_safety"]),
            )
            result = run_renormalized_cascade(cascade_config)
            gains = result.gains
            runs.append(
                {
                    "energy_constant": float(energy_constant),
                    "drop_below": drop_below,
                    "classification": _classify(
                        gains, result.terminated_reason
                    ),
                    "final_gain": gains[-1] if gains else None,
                    "max_gain": max(gains) if gains else None,
                    "final_overlap_lag1": (
                        result.stage_records[-1].overlap_lag1
                        if result.stage_records
                        else None
                    ),
                    "final_overlap_lag2": (
                        result.stage_records[-1].overlap_lag2
                        if result.stage_records
                        else None
                    ),
                    "result": result.as_dict(),
                }
            )
            print(
                f"c_E={energy_constant} drop_below={drop_below} "
                f"-> {runs[-1]['classification']} "
                f"({result.terminated_reason}, "
                f"{result.completed_stages} stages) "
                f"final_gain={runs[-1]['final_gain']}",
                flush=True,
            )

    plateau_runs = [run for run in runs if run["classification"] == "plateau"]
    growing_runs = [run for run in runs if run["classification"] == "growing"]
    summary = {
        "schema": OUTPUT_SCHEMA,
        "status": STATUS,
        "config": config,
        "runs": runs,
        "verdict": {
            "any_plateau": bool(plateau_runs),
            "any_growing": bool(growing_runs),
            "interpretation": (
                "a positive-gain plateau or growth at any single c_E is the "
                "survival signal (by the c_E-collapse the plateau level is "
                "tunable); uniform decay rejects attracting-orbit closure in "
                "this scanned box only"
            ),
        },
    }
    write_with_digest(
        arguments.output_dir / "config.snapshot.json",
        canonical_json_bytes(config),
    )
    write_with_digest(
        arguments.output_dir / "summary.json",
        canonical_json_bytes(summary),
    )
    print(json.dumps(summary["verdict"], indent=2))


if __name__ == "__main__":
    main()
