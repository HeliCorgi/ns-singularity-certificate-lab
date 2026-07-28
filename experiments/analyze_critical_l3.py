"""Analyze critical-L3 diagnostics on nonlinear cylinder checkpoints.

This command is deliberately post-processing only.  Values computed on a
finite periodic cylinder are domain-truncated observations, not the global
``L^3(R^3)`` norm required by the endpoint regularity theorem.
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

from ns_certificate_lab.critical_l3 import critical_growth_series, state_critical_l3
from ns_certificate_lab.nonlinear_cylinder import load_checkpoint


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _prepare_output(path: Path) -> None:
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(f"output path is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(f"refusing to overwrite nonempty directory: {path}")
    else:
        path.mkdir(parents=True)


def _finite(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        raise FloatingPointError("non-finite value in critical-L3 output")
    if isinstance(value, dict):
        return {key: _finite(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(item) for item in value]
    return value


def run(
    checkpoints: list[Path],
    output_dir: Path,
    *,
    shell_count: int,
    outer_radial_cells: int,
) -> dict[str, Any]:
    if not checkpoints:
        raise ValueError("at least one checkpoint is required")
    _prepare_output(output_dir)

    records: list[dict[str, Any]] = []
    grouped: dict[tuple[int, int], list[tuple[float, Any]]] = {}
    for raw_path in checkpoints:
        path = raw_path.resolve()
        checkpoint = load_checkpoint(path)
        report = state_critical_l3(
            checkpoint.grid,
            u1=checkpoint.state.u1,
            u_r=checkpoint.state.u_r,
            u_z=checkpoint.state.u_z,
            shell_count=shell_count,
            outer_radial_cells=outer_radial_cells,
        )
        record = {
            "checkpoint": str(path),
            "checkpoint_sha256": _sha256(path),
            "time": checkpoint.time,
            "viscosity": checkpoint.viscosity,
            "nr": checkpoint.grid.nr,
            "nz": checkpoint.grid.nz,
            "r_max": float(checkpoint.grid.r[-1]),
            "z_period": (
                float(checkpoint.grid.z_period)
                if checkpoint.grid.periodic_z
                else None
            ),
            **asdict(report),
        }
        records.append(_finite(record))
        grouped.setdefault((checkpoint.grid.nr, checkpoint.grid.nz), []).append(
            (checkpoint.time, report)
        )

    records.sort(key=lambda item: (item["nr"], item["nz"], item["time"]))
    growth: list[dict[str, Any]] = []
    for (nr, nz), entries in sorted(grouped.items()):
        entries.sort(key=lambda item: item[0])
        # Duplicate times cannot define a growth series.  Keep every snapshot
        # in the scalar table, but skip the derivative report for that group.
        times = [item[0] for item in entries]
        if len(entries) < 2 or any(b <= a for a, b in zip(times, times[1:])):
            continue
        series = critical_growth_series(times, [item[1] for item in entries])
        growth.append({"nr": nr, "nz": nz, **asdict(series)})

    csv_path = output_dir / "critical_l3_snapshots.csv"
    scalar_fields = [
        "checkpoint",
        "checkpoint_sha256",
        "time",
        "viscosity",
        "nr",
        "nz",
        "r_max",
        "z_period",
        "l3_cubed",
        "l3_norm",
        "max_velocity",
        "center_z",
        "center_resultant_ratio",
        "radial_scale",
        "axial_scale",
        "isotropic_scale",
        "critical_scaling_product",
        "outer_radial_fraction",
        "shell_entropy",
        "effective_shell_count",
        "maximum_shell_fraction",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=scalar_fields)
        writer.writeheader()
        for record in records:
            writer.writerow({name: record[name] for name in scalar_fields})

    summary: dict[str, Any] = {
        "schema_version": 1,
        "interpretation": (
            "Floating-point critical-L3 diagnostics on represented finite "
            "domains. These are not global R3 norms, interval bounds, a "
            "singularity candidate, or a proof."
        ),
        "definitions": {
            "critical_density": "|u|^3",
            "volume_measure": "2*pi*r dr dz",
            "amplitude": "max |u| on the represented grid",
            "radial_scale": "RMS radius of |u|^3 about the symmetry axis",
            "axial_scale": "periodic RMS axial width of |u|^3",
            "critical_scaling_product": "A^3 * L_r^2 * L_z",
            "shell_geometry": (
                "axis-centered dyadic shells around the critical-density "
                "axial center"
            ),
        },
        "limitations": [
            "Periodic-z cylinder values are not L3(R3).",
            "The RMS scaling product is a shape diagnostic, not an exact "
            "self-similar coefficient unless profile stability is established.",
            "Axis-centered shells do not describe an off-axis vortex ring as "
            "a point singularity.",
            "No domain-tail, discretization, roundoff or interval enclosure "
            "is supplied by this command.",
        ],
        "snapshots": records,
        "growth_series": growth,
    }
    summary_path = output_dir / "critical_l3_summary.json"
    summary_path.write_text(
        json.dumps(_finite(summary), indent=2, sort_keys=True, allow_nan=False) + "\n",
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
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoints", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--shell-count", type=int, default=10)
    parser.add_argument("--outer-radial-cells", type=int, default=4)
    arguments = parser.parse_args(argv)
    run(
        arguments.checkpoints,
        arguments.output_dir.resolve(),
        shell_count=arguments.shell_count,
        outer_radial_cells=arguments.outer_radial_cells,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
