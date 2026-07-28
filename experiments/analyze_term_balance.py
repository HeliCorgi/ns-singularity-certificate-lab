"""Post-process checkpoint series into transformed-PDE term balances.

The result is a floating-point diagnostic.  Sparse checkpoint time derivatives
are not a posteriori error bounds and cannot certify an asymptotic regime.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from ns_certificate_lab.critical_l3 import state_critical_l3
from ns_certificate_lab.nonlinear_cylinder import load_checkpoint
from ns_certificate_lab.term_balance import (
    critical_core_mask,
    derivative_sensitivity,
    snapshot_balance_report,
    time_derivatives,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _prepare_output(path: Path) -> None:
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(path)
        if any(path.iterdir()):
            raise FileExistsError(f"refusing to overwrite nonempty directory: {path}")
    else:
        path.mkdir(parents=True)


def _finite(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        raise FloatingPointError("non-finite term-balance output")
    if isinstance(value, dict):
        return {key: _finite(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(item) for item in value]
    return value


def _flatten(prefix: str, values: dict[str, Any], target: dict[str, Any]) -> None:
    for key, value in values.items():
        target[f"{prefix}_{key}"] = value


def run(checkpoints: list[Path], output_dir: Path, *, core_multiplier: float) -> dict[str, Any]:
    if not checkpoints:
        raise ValueError("at least one checkpoint is required")
    _prepare_output(output_dir)

    groups: dict[tuple[int, int], list[tuple[Path, Any]]] = {}
    for raw in checkpoints:
        path = raw.resolve()
        checkpoint = load_checkpoint(path)
        groups.setdefault((checkpoint.grid.nr, checkpoint.grid.nz), []).append(
            (path, checkpoint)
        )

    records: list[dict[str, Any]] = []
    for (nr, nz), entries in sorted(groups.items()):
        entries.sort(key=lambda item: item[1].time)
        if len(entries) < 3:
            raise ValueError(f"group {nr}x{nz} needs at least three snapshots")
        times = [item[1].time for item in entries]
        u1_series = [item[1].state.u1 for item in entries]
        omega_series = [item[1].state.omega1 for item in entries]
        du1 = time_derivatives(times, u1_series)
        domega = time_derivatives(times, omega_series)
        du1_sensitivity = derivative_sensitivity(times, u1_series, du1)
        domega_sensitivity = derivative_sensitivity(times, omega_series, domega)

        for index, (path, checkpoint) in enumerate(entries):
            critical = state_critical_l3(
                checkpoint.grid,
                u1=checkpoint.state.u1,
                u_r=checkpoint.state.u_r,
                u_z=checkpoint.state.u_z,
            )
            mask = critical_core_mask(
                checkpoint.grid,
                center_z=critical.center_z,
                radial_scale=critical.radial_scale,
                axial_scale=critical.axial_scale,
                radius_multiplier=core_multiplier,
            )
            report = snapshot_balance_report(
                checkpoint.grid,
                time=checkpoint.time,
                viscosity=checkpoint.viscosity,
                u1=checkpoint.state.u1,
                omega1=checkpoint.state.omega1,
                psi1=checkpoint.state.psi1,
                du1_dt=du1[index],
                domega1_dt=domega[index],
                du1_dt_sensitivity=du1_sensitivity[index],
                domega1_dt_sensitivity=domega_sensitivity[index],
                core_mask=mask,
            )
            raw_report = asdict(report)
            record: dict[str, Any] = {
                "checkpoint": str(path),
                "checkpoint_sha256": _sha256(path),
                "nr": nr,
                "nz": nz,
                "time": report.time,
                "viscosity": report.viscosity,
                "core_fraction": report.core_fraction,
                "critical_l3_norm": critical.l3_norm,
                "critical_radial_scale": critical.radial_scale,
                "critical_axial_scale": critical.axial_scale,
            }
            for name in ("u1_global", "omega1_global", "u1_core", "omega1_core"):
                _flatten(name, raw_report[name], record)
            records.append(_finite(record))

    records.sort(key=lambda item: (item["nr"], item["nz"], item["time"]))
    csv_path = output_dir / "term_balance.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    summary = {
        "schema_version": 1,
        "interpretation": (
            "Floating-point transformed-PDE balance diagnostic. Time derivatives "
            "come from sparse checkpoint interpolation and are not certified error "
            "bounds or proof of an asymptotic regime."
        ),
        "core_definition": (
            f"ellipse of radius {core_multiplier} in critical-density RMS coordinates"
        ),
        "classification_rule": (
            "classification requires normalized residual <= 0.25; active blocks are "
            "within a factor four of the largest of time, inviscid, diffusion"
        ),
        "records": records,
    }
    summary_path = output_dir / "term_balance.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": 1,
        "files": {
            csv_path.name: {"sha256": _sha256(csv_path), "bytes": csv_path.stat().st_size},
            summary_path.name: {
                "sha256": _sha256(summary_path),
                "bytes": summary_path.stat().st_size,
            },
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoints", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--core-multiplier", type=float, default=2.0)
    args = parser.parse_args()
    run(args.checkpoints, args.output_dir, core_multiplier=args.core_multiplier)


if __name__ == "__main__":
    main()
