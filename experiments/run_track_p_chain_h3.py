"""Replay the preregistered turn-11 ``n = 3`` Kato-constant chain runs.

Builds the Kato-constant certificate once (verified), then the chains of
``configs/track_p_chain_h3.json``, verifies each with the independent
checker, and writes payloads, a comparison summary against the turn-10
``H^4`` baseline, and a sha256 manifest to ``outputs/track_p_chain_h3_v1/``.

The improvement ratio is computed from the emitted certificates; nothing here
presumes it.  Nothing here bears on the Clay problem.
"""

from __future__ import annotations

import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

from ns_certificate_lab.kato_constant import (
    build_kato_certificate,
    verify_kato_certificate,
)
from ns_certificate_lab.torus_chain import (
    build_chain_certificate,
    verify_chain_certificate,
)

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "configs" / "track_p_chain_h3.json"
BASELINE = ROOT / "outputs" / "track_p_chain_v1" / "summary.json"
OUTDIR = ROOT / "outputs" / "track_p_chain_h3_v1"


def run_chain(name: str, family: str, viscosity: str, block: dict,
              kato_payload: dict) -> dict:
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
        mode=block["mode"],
        kato_payload=kato_payload,
    )
    build_seconds = time.time() - started
    started = time.time()
    verdict = verify_chain_certificate(payload)
    verify_seconds = time.time() - started

    (OUTDIR / f"{name}.json").write_text(
        json.dumps(payload, indent=1), encoding="utf-8"
    )
    order = payload["norm_order"]
    row = {
        "name": name,
        "family": family,
        "viscosity": viscosity,
        "cutoff_sq": block["cutoff_sq"],
        "initial_step": block["initial_step"],
        "slab_count": payload["slab_count"],
        "certified_final_time": payload["certified_final_time"],
        "max_radius": payload[f"max_radius_h{order}_dot"],
        "final_delta": payload[f"final_delta_h{order}_dot"],
        "termination": payload["termination"]["reason"],
        "termination_detail": payload["termination"].get("detail"),
        "verified": verdict["verified"],
        "proves_existence": verdict["proves_existence"],
        "checker_failures": verdict["failures"],
        "build_seconds": round(build_seconds, 1),
        "verify_seconds": round(verify_seconds, 1),
        "seconds_per_slab": round(
            build_seconds / max(payload["slab_count"], 1), 1
        ),
    }
    print(
        f"{name}: slabs={row['slab_count']} "
        f"T={row['certified_final_time']} "
        f"({float(Fraction(row['certified_final_time'])):.3e}) "
        f"maxR={float(Fraction(row['max_radius'])):.3e} "
        f"stop={row['termination']} verified={row['verified']} "
        f"[{row['build_seconds']:.0f}s+{row['verify_seconds']:.0f}s]",
        flush=True,
    )
    return row


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    OUTDIR.mkdir(parents=True, exist_ok=True)

    kato_cfg = config["kato"]
    kato_payload = build_kato_certificate(
        lattice_cut=int(kato_cfg["lattice_cut"]), bits=int(kato_cfg["bits"])
    )
    if hasattr(kato_payload, "as_dict"):
        kato_payload = kato_payload.as_dict()
    kato_verdict = verify_kato_certificate(kato_payload)
    if not kato_verdict.get("verified"):
        raise SystemExit(f"Kato certificate failed: {kato_verdict['failures']}")
    (OUTDIR / "kato_certificate.json").write_text(
        json.dumps(kato_payload, indent=1), encoding="utf-8"
    )
    print(f"kato certificate: G3_upper = {kato_payload['g3']['upper']} "
          f"(~{float(Fraction(kato_payload['g3']['upper'])):.4f}), verified",
          flush=True)

    rows = []
    sweep = config["sweep"]
    for family in sweep["families"]:
        for viscosity in sweep["viscosities"]:
            name = f"h3chain_{family}_nu_{viscosity.replace('/', 'over')}"
            rows.append(run_chain(name, family, viscosity, sweep, kato_payload))

    strict = config["strict_same_step"]
    rows.append(
        run_chain("h3chain_strict_same_step_P1_nu_1over10",
                  strict["family"], strict["viscosity"], strict, kato_payload)
    )
    probe = config["cutoff_probe"]
    rows.append(
        run_chain("h3chain_cutoff6_P1_nu_1over10",
                  probe["family"], probe["viscosity"], probe, kato_payload)
    )

    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    old = {
        (row["family"], row["viscosity"]): row
        for row in baseline["rows"]
        if not row["name"].startswith("chain_long")
    }
    comparison = []
    for row in rows:
        key = (row["family"], row["viscosity"])
        if row["name"].startswith("h3chain_") and key in old and \
                row["cutoff_sq"] == 4 and row["initial_step"] == sweep["initial_step"]:
            ratio = float(
                Fraction(row["certified_final_time"])
                / Fraction(old[key]["certified_final_time"])
            )
            comparison.append(
                {
                    "family": row["family"],
                    "viscosity": row["viscosity"],
                    "old_T": old[key]["certified_final_time"],
                    "old_slabs": old[key]["slab_count"],
                    "new_T": row["certified_final_time"],
                    "new_slabs": row["slab_count"],
                    "horizon_ratio_new_over_old": round(ratio, 2),
                    "old_termination": old[key]["termination"],
                    "new_termination": row["termination"],
                }
            )

    summary = {
        "config": config,
        "kato_G3_upper": kato_payload["g3"]["upper"],
        "rows": rows,
        "comparison_vs_h4_baseline": comparison,
        "all_verified": all(row["verified"] for row in rows),
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
    print("comparison:", json.dumps(comparison), flush=True)
    print("all_verified:", summary["all_verified"], flush=True)


if __name__ == "__main__":
    main()
