"""Reissue every chain certificate with the audited EXT closure records.

Turn 11: the adversarial audit (docs/research_notes/ext_p1_p2_p3_audit.md)
closed EXT-P1*/EXT-P2-INT/EXT-P3* (+ Corollary P3-3) as audited paper
proofs.  This driver reissues the turn-10 ``H^4`` chains and the turn-11
``n = 3`` chains via :func:`torus_chain.reissue_chain_certificate`, re-runs
the full independent checker on every reissued payload (which recomputes
every link), and writes ``outputs/track_p_chain_reissued_v2/``.

The originals are kept untouched as historical conditional artifacts.
``proved: true`` means an audited paper proof of classical analysis — never
a Lean formalisation; no axiom exists or will exist for it.  Nothing here
bears on the Clay problem.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

from ns_certificate_lab.torus_chain import (
    reissue_chain_certificate,
    verify_chain_certificate,
)

ROOT = Path(__file__).resolve().parent.parent
SOURCES = [
    ROOT / "outputs" / "track_p_chain_v1",
    ROOT / "outputs" / "track_p_chain_h3_v1",
]
OUTDIR = ROOT / "outputs" / "track_p_chain_reissued_v2"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in SOURCES:
        for path in sorted(source.glob("*.json")):
            if path.name in ("summary.json", "manifest.json",
                             "kato_certificate.json"):
                continue
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("certificate") != "torus_chain_aposteriori":
                continue
            reissued = reissue_chain_certificate(payload)
            started = time.time()
            verdict = verify_chain_certificate(reissued)
            seconds = time.time() - started
            out_name = f"reissued_{path.name}"
            (OUTDIR / out_name).write_text(
                json.dumps(reissued, indent=1), encoding="utf-8"
            )
            rows.append(
                {
                    "source": f"{source.name}/{path.name}",
                    "reissued": out_name,
                    "mode": reissued.get("mode", "h4_crude"),
                    "slab_count": reissued.get("slab_count"),
                    "certified_final_time": reissued.get("certified_final_time"),
                    "verified": verdict["verified"],
                    "failures": verdict["failures"],
                    "verify_seconds": round(seconds, 1),
                }
            )
            print(
                f"{out_name}: verified={verdict['verified']} "
                f"[{seconds:.0f}s]",
                flush=True,
            )

    summary = {
        "rows": rows,
        "all_verified": all(row["verified"] for row in rows),
        "meaning": (
            "reissued under the audited closure records "
            "(EXTERNAL_THEOREMS_AUDITED); quantitative content untouched; "
            "proved:true = audited paper proof, never Lean; originals kept "
            "as historical conditional artifacts"
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
    print("all_verified:", summary["all_verified"], flush=True)


if __name__ == "__main__":
    main()
