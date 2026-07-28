"""Fit conditional singular scaling exponents from a critical-L3 CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ns_certificate_lab.scaling_fit import scan_joint_scaling_fits, summarize_fit_scan


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--output", type=Path, default=Path("outputs/scaling_fit_summary.json"))
    parser.add_argument("--minimum-points", type=int, default=3)
    parser.add_argument("--time-padding-min", type=float, default=1.0e-6)
    parser.add_argument("--time-padding-max", type=float, default=2.0e-3)
    parser.add_argument("--time-count", type=int, default=400)
    args = parser.parse_args()

    rows = _read_rows(args.csv_path)
    groups: dict[tuple[int, int], list[dict[str, str]]] = {}
    for row in rows:
        key = (int(row["nr"]), int(row["nz"]))
        groups.setdefault(key, []).append(row)

    output: dict[str, object] = {
        "source_csv": str(args.csv_path),
        "groups": {},
        "warning": (
            "Candidate time and window are scanned. Results are diagnostics only and "
            "must not be interpreted as singularity evidence without preregistered, "
            "resolved and independently converged follow-up runs."
        ),
    }
    group_output = output["groups"]
    assert isinstance(group_output, dict)

    for key, group_rows in sorted(groups.items()):
        ordered = sorted(group_rows, key=lambda row: float(row["time"]))
        times = np.asarray([float(row["time"]) for row in ordered], dtype=np.float64)
        amplitudes = np.asarray(
            [float(row["max_velocity"]) for row in ordered], dtype=np.float64
        )
        radial = np.asarray(
            [float(row["radial_scale"]) for row in ordered], dtype=np.float64
        )
        axial = np.asarray(
            [float(row["axial_scale"]) for row in ordered], dtype=np.float64
        )
        candidate_times = np.linspace(
            float(times[-1] + args.time_padding_min),
            float(times[-1] + args.time_padding_max),
            args.time_count,
        )
        fits = scan_joint_scaling_fits(
            times,
            amplitudes,
            radial,
            axial,
            singular_times=candidate_times,
            minimum_points=args.minimum_points,
        )
        group_output[f"nr{key[0]}_nz{key[1]}"] = summarize_fit_scan(fits)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
