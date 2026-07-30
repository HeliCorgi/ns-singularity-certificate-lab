r"""Track P — one rigorously certified slab per periodic family.

For each preregistered family and each preregistered slab length this driver
builds the full a posteriori certificate of
:mod:`ns_certificate_lab.torus_aposteriori`: the Picard enclosure of the
Galerkin trajectory, the slab-valid constants of the ``H^4`` control
inequality, the exact continuous residual, the control-ODE tube, and the
conditional conclusion that a strong periodic Navier–Stokes solution exists on
the slab within the stated ``\dot H^4`` distance of the trajectory.

Every payload is written in full and re-checked by the independent verifier
before it is accepted.  The acceptance thresholds are preregistered in the
config; a slab that fails them is recorded as failed, not dropped.

This experiment certifies **regularity near computed trajectories**.  It is not
a singularity computation and says nothing about the Clay problem.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

from ns_certificate_lab.torus_aposteriori import (
    build_torus_slab_certificate,
    verify_torus_slab_certificate,
)

EXPECTED_SCHEMA_VERSION = 1
EXPECTED_EXPERIMENT = "track_p_slab"


def _write_json(path: Path, data: object) -> None:
    text = (
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    )
    path.write_text(text, encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate(config: dict[str, Any]) -> None:
    required = {"schema_version", "experiment", "preregistration", "lattice_cut",
                "precision_bits"}
    if set(config) != required:
        raise ValueError("Track P config has missing or unknown keys")
    if config["schema_version"] != EXPECTED_SCHEMA_VERSION:
        raise ValueError("unsupported Track P config schema")
    if config["experiment"] != EXPECTED_EXPERIMENT:
        raise ValueError("experiment must equal the audited canonical value")
    if not config["preregistration"].get("recorded_before_any_run"):
        raise ValueError("the preregistration block must assert its own timing")


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    _validate(config)
    prereg = config["preregistration"]
    output_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []

    def write_json(name: str, payload: object) -> None:
        path = output_dir / name
        _write_json(path, payload)
        written.append(path)

    write_json("preregistration.json", prereg)
    write_json("config.snapshot.json", config)

    viscosity = Fraction(prereg["viscosity"])
    cutoff = int(prereg["galerkin_cutoff_sq"])
    acceptance = prereg["acceptance"]
    ratio_bar = Fraction(str(acceptance["distance_h4_sq_below_datum_h4_sq_ratio"]))

    rows: list[dict[str, Any]] = []
    for family in prereg["families"]:
        for step_text in prereg["steps"]:
            step = Fraction(step_text)
            certificate = build_torus_slab_certificate(
                family,
                viscosity=viscosity,
                cutoff_sq=cutoff,
                step=step,
                lattice_cut=int(config["lattice_cut"]),
                bits=int(config["precision_bits"]),
            )
            payload = certificate.as_dict()
            verdict = verify_torus_slab_certificate(payload)
            conclusion = payload["conclusion"]
            distance = (
                Fraction(conclusion["distance_h4_dot"])
                if conclusion.get("distance_h4_dot")
                else None
            )
            datum_sq = Fraction(payload["datum"]["h4_norm_sq"])
            ratio_ok = (
                distance is not None and distance * distance <= ratio_bar * datum_sq
            )
            accepted = bool(
                payload["galerkin_box"]["proved"]
                and conclusion.get("proved")
                and verdict["verified"]
                and ratio_ok
            )
            name = f"slab_{family}_h{step.denominator}.json"
            write_json(name, payload)
            rows.append(
                {
                    "family": family,
                    "step": str(step),
                    "box_proved": payload["galerkin_box"]["proved"],
                    "box_step": payload["galerkin_box"]["step"],
                    "conclusion_proved": bool(conclusion.get("proved")),
                    "checker_verified": bool(verdict["verified"]),
                    "distance_h4_dot": float(distance) if distance is not None else None,
                    "distance_h3": (
                        float(Fraction(conclusion["distance_h3"]))
                        if conclusion.get("distance_h3")
                        else None
                    ),
                    "datum_h4": float(datum_sq) ** 0.5,
                    "relative_distance": (
                        float(distance) / float(datum_sq) ** 0.5
                        if distance is not None
                        else None
                    ),
                    "accepted": accepted,
                    "payload": name,
                }
            )

    summary = {
        "schema_version": EXPECTED_SCHEMA_VERSION,
        "experiment": EXPECTED_EXPERIMENT,
        "preregistration": prereg,
        "rows": rows,
        "accepted": [
            f"{row['family']}@{row['step']}" for row in rows if row["accepted"]
        ],
        "limitations": [
            "Each certificate is conditional on the three named classical "
            "external theorems EXT-P1/P2/P3 recorded in its payload; they are "
            "classical, unproved here, and never axiomatised in Lean.",
            "One slab per run; chaining slabs needs Lohner-type frames and is "
            "future work.",
            "The self-contained constants are crude (a sharp Kato constant "
            "would lengthen the certified interval); this costs slab length, "
            "not validity.",
            "This certifies regularity near computed trajectories.  It is not "
            "a singularity result and says nothing about the Clay problem.",
        ],
    }
    write_json("summary.json", summary)

    manifest = {
        "schema_version": 1,
        "payloads": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in written
        },
    }
    manifest_path = output_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    (output_dir / "manifest.json.sha256").write_text(
        _sha256(manifest_path) + "\n", encoding="ascii"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.config, args.output_dir)
    print(
        json.dumps(
            {
                "accepted": summary["accepted"],
                "rows": [
                    {
                        key: row[key]
                        for key in (
                            "family", "step", "conclusion_proved",
                            "relative_distance", "accepted",
                        )
                    }
                    for row in summary["rows"]
                ],
            },
            ensure_ascii=False, indent=2, sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
