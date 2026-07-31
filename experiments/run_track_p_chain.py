"""Replay the preregistered Track-P chain runs (turn 10).

Reads ``configs/track_p_chain.json`` (preregistered before any run), builds
every chain, verifies each with the independent checker, and writes payloads,
a summary and a sha256 manifest to ``outputs/track_p_chain_v1/``.

Every chain is emitted whether it reaches its budget or stops early; a stop is
classified, never interpreted.  Nothing here bears on the Clay problem.
"""

from __future__ import annotations

import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

from ns_certificate_lab.torus_chain import (
    build_chain_certificate,
    verify_chain_certificate,
)

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "configs" / "track_p_chain.json"
OUTDIR = ROOT / "outputs" / "track_p_chain_v1"


def run_chain(name: str, family: str, viscosity: str, block: dict) -> dict:
    started = time.time()
    payload = build_chain_certificate(
        family,
        viscosity=Fraction(viscosity),
        cutoff_sq=int(block["cutoff_sq"]),
        initial_step=Fraction(block["initial_step"]),
        min_step=Fraction(block["min_step"]),
        max_slabs=int(block["max_slabs"]),
        lattice_cut=int(block["lattice_cut"]),
        bits=int(block["bits"]),
        taylor_order=int(block["taylor_order"]),
        max_taylor_order=int(block["max_taylor_order"]),
        rounding_bits=int(block["rounding_bits"]),
        relative_radius_cap=Fraction(block["relative_radius_cap"]),
        allow_step_doubling=bool(block["allow_step_doubling"]),
    )
    build_seconds = time.time() - started
    started = time.time()
    verdict = verify_chain_certificate(payload)
    verify_seconds = time.time() - started

    path = OUTDIR / f"{name}.json"
    path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    row = {
        "name": name,
        "family": family,
        "viscosity": viscosity,
        "slab_count": payload["slab_count"],
        "certified_final_time": payload["certified_final_time"],
        "horizon_over_single_slab": payload["horizon_over_single_slab"],
        "max_radius_h4_dot": payload["max_radius_h4_dot"],
        "final_delta_h4_dot": payload["final_delta_h4_dot"],
        "termination": payload["termination"]["reason"],
        "lohner_trigger_fired": payload["wrapping_diagnostics"]["fired"],
        "verified": verdict["verified"],
        "proves_existence": verdict["proves_existence"],
        "checker_failures": verdict["failures"],
        "build_seconds": round(build_seconds, 1),
        "verify_seconds": round(verify_seconds, 1),
    }
    print(
        f"{name}: slabs={row['slab_count']} T={row['certified_final_time']} "
        f"({float(Fraction(row['certified_final_time'])):.3e}) "
        f"maxR={float(Fraction(row['max_radius_h4_dot'])):.3e} "
        f"stop={row['termination']} verified={row['verified']} "
        f"[{row['build_seconds']:.0f}s+{row['verify_seconds']:.0f}s]",
        flush=True,
    )
    return row


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    OUTDIR.mkdir(parents=True, exist_ok=True)

    rows = []
    sweep = config["sweep"]
    for family in sweep["families"]:
        for viscosity in sweep["viscosities"]:
            name = f"chain_{family}_nu_{viscosity.replace('/', 'over')}"
            rows.append(run_chain(name, family, viscosity, sweep))

    long_block = config["long"]
    rows.append(
        run_chain(
            f"chain_long_{long_block['family']}_nu_"
            f"{long_block['viscosity'].replace('/', 'over')}",
            long_block["family"],
            long_block["viscosity"],
            long_block,
        )
    )

    acceptance = config["acceptance"]
    reference = Fraction(acceptance["horizon_reference_step"])
    evaluation = {
        "every_chain_at_least_slabs": all(
            row["slab_count"] >= acceptance["every_chain_at_least_slabs"]
            for row in rows
        ),
        "at_least_one_chain_horizon_over_reference": any(
            Fraction(row["certified_final_time"]) / reference
            >= Fraction(acceptance["at_least_one_chain_horizon_over_reference"])
            for row in rows
        ),
        "long_chain_target_slabs_reached": rows[-1]["slab_count"]
        >= acceptance["long_chain_target_slabs"],
        "checker_verified_every_chain": all(row["verified"] for row in rows),
        "every_termination_classified": all(
            row["termination"] is not None for row in rows
        ),
    }
    summary = {
        "config": config,
        "rows": rows,
        "acceptance_evaluation": evaluation,
        "disclaimer": (
            "conditional PDE certificates assuming EXT-P1/P2/P3; the "
            "finite-dimensional Galerkin enclosure and the scalar control ODE "
            "are verified unconditionally; nothing here bears on the Clay "
            "problem and a stopped chain is not a singularity statement"
        ),
    }
    (OUTDIR / "summary.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8"
    )

    manifest = {}
    for path in sorted(OUTDIR.glob("*.json")):
        if path.name == "manifest.json":
            continue
        manifest[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUTDIR / "manifest.json").write_text(
        json.dumps(manifest, indent=1), encoding="utf-8"
    )
    print("acceptance:", json.dumps(evaluation), flush=True)


if __name__ == "__main__":
    main()
