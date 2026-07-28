"""Scan conditional Type-II/anisotropic power-law exponent families."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ns_certificate_lab.scaling_constraints import scan_power_law_scalings


def _range(start: float, stop: float, step: float) -> np.ndarray:
    if step <= 0.0:
        raise ValueError("step must be positive")
    count = int(np.floor((stop - start) / step + 0.5))
    return start + step * np.arange(count + 1, dtype=np.float64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha-min", type=float, default=0.40)
    parser.add_argument("--alpha-max", type=float, default=0.80)
    parser.add_argument("--beta-min", type=float, default=0.20)
    parser.add_argument("--beta-max", type=float, default=0.80)
    parser.add_argument("--step", type=float, default=0.01)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/scaling_constraints"))
    args = parser.parse_args()

    alphas = _range(args.alpha_min, args.alpha_max, args.step)
    betas = _range(args.beta_min, args.beta_max, args.step)
    reports = scan_power_law_scalings(alphas, betas, betas, tolerance=args.step * 0.05)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = [report.to_dict() for report in reports]
    flat_rows: list[dict[str, object]] = []
    counts: dict[str, int] = {}
    for row in rows:
        scaling = row.pop("scaling")
        assert isinstance(scaling, dict)
        flat = {**scaling, **row}
        flat_rows.append(flat)
        family = str(flat["candidate_family"])
        counts[family] = counts.get(family, 0) + 1

    csv_path = args.output_dir / "feasible_scalings.csv"
    if flat_rows:
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0]))
            writer.writeheader()
            writer.writerows(flat_rows)
    else:
        csv_path.write_text("", encoding="utf-8")

    summary = {
        "grid": {
            "alpha_min": args.alpha_min,
            "alpha_max": args.alpha_max,
            "beta_min": args.beta_min,
            "beta_max": args.beta_max,
            "step": args.step,
            "alpha_count": int(alphas.size),
            "beta_count": int(betas.size),
            "total_grid_points": int(alphas.size * betas.size * betas.size),
        },
        "conditional_feasible_count": len(flat_rows),
        "family_counts": counts,
        "warning": (
            "These are conditional asymptotic scaling points, not PDE solutions. "
            "They assume a localized nondegenerate core and no leading cancellation."
        ),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
