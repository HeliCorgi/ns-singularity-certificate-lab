"""Run the exact spectral front-gap ledger pilot.

Exact lane: evaluates the closable front-gap identity (I.3)/(I.4), the
saturation deficit, and the Lemma-K lattice certificate on the repository's
exact rational fields (relay triad at two scales, the four-parent expanded
carrier field).  Float lane: sweeps the amplitude-delocalization ratio
``mu_N = M_eff/N^3`` of the fixed-relative mesoscopic parents over
``eta x N``; a scale-independent positive floor is what the mesoscopic no-go
requires of the only surviving cloud family, and ``mu_N -> 0`` would kill
that family by an independent argument.  Discovery diagnostic; not a proof.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    write_with_digest,
)
from ns_certificate_lab.exact_leray_relay import build_exact_relay_triad
from ns_certificate_lab.fourier_torus import TrigVector
from ns_certificate_lab.mesoscopic_cloud_scaling import (
    MesoscopicCloudConfig,
    build_sparse_parent,
)
from ns_certificate_lab.spectral_front_monotone import (
    front_gap_identity,
    lemma_k_certificate,
    sparse_parent_delocalization,
)

OUTPUT_SCHEMA = "ns-certificate-lab/spectral-front-ledger/v1"
STATUS = "EXACT RATIONAL LEDGER + BINARY64 SWEEP / NOT A PROOF"


def _four_parent_field() -> TrigVector:
    third = Fraction(1, 3)
    return TrigVector.from_modes(
        [
            ((1, 1, 0), (0, 0, 0), (0, 0, 1)),
            ((1, 0, 1), (0, 1, 0), (0, 0, 0)),
            ((0, 1, -1), (third, 2 * third, 2 * third), (0, 0, 0)),
            ((1, 0, -1), (2 * third, third, 2 * third), (0, 0, 0)),
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    fields = {
        "relay_triad_s1": build_exact_relay_triad(scale=1),
        "relay_triad_s2": build_exact_relay_triad(scale=2),
        "four_parent_expanded_carrier": _four_parent_field(),
    }
    exact_records: list[dict[str, Any]] = []
    for name, field in fields.items():
        entry: dict[str, Any] = {
            "field": name,
            "lemma_k": lemma_k_certificate(field),
            "orders": [],
        }
        for viscosity in (Fraction(1, 40), Fraction(1, 10)):
            for order in (0, 1, 2):
                record = front_gap_identity(
                    field, order=order, viscosity=viscosity
                )
                entry["orders"].append(record.as_dict())
        exact_records.append(entry)
        print(
            f"{name}: K={entry['lemma_k']['front_wavenumber_K']} "
            f"S_N={entry['lemma_k']['lattice_sum_S']}",
            flush=True,
        )

    sweep: list[dict[str, Any]] = []
    for eta in (0.10, 0.15, 0.20, 0.25, 0.30):
        for scale in (8, 16, 24, 32, 48, 64):
            width = max(1, math.floor(eta * scale))
            if scale <= 3 * (width - 1):
                sweep.append(
                    {
                        "eta": eta,
                        "scale": scale,
                        "width": width,
                        "status": "geometry_excluded",
                    }
                )
                continue
            config = MesoscopicCloudConfig(
                base_scale=scale, gamma=1.0, width_override=width
            )
            parent = build_sparse_parent(config, maximum_modes=400_000)
            record = sparse_parent_delocalization(parent, scale=scale)
            sweep.append(
                {
                    "eta": eta,
                    "scale": scale,
                    "width": width,
                    "status": "measured",
                    **record,
                }
            )
            print(
                f"eta={eta} N={scale} W={width} "
                f"mu_N={record['mu_delocalization']:.6f} "
                f"M_eff/M={record['effective_over_support']:.4f}",
                flush=True,
            )

    per_eta: dict[str, dict[str, float]] = {}
    for eta in (0.10, 0.15, 0.20, 0.25, 0.30):
        values = [
            row["mu_delocalization"]
            for row in sweep
            if row["status"] == "measured" and row["eta"] == eta
        ]
        if len(values) >= 3:
            per_eta[str(eta)] = {
                "count": float(len(values)),
                "first": values[0],
                "last": values[-1],
                "minimum": min(values),
                "maximum": max(values),
                "last_over_first": values[-1] / values[0],
            }

    summary = {
        "schema": OUTPUT_SCHEMA,
        "status": STATUS,
        "exact_records": exact_records,
        "delocalization_sweep": sweep,
        "per_eta_mu_trend": per_eta,
        "verdict_notes": (
            "K1/K2 kill conditions (negative gap, failed identity, negative "
            "Lemma-K margin) are enforced by assertions inside the exact "
            "lane; reaching this summary means none fired. The mu_N sweep "
            "tests the K4 question: a positive scale-independent floor keeps "
            "the fixed-relative family alive, mu_N->0 kills it."
        ),
    }
    write_with_digest(
        arguments.output_dir / "summary.json",
        canonical_json_bytes(summary),
    )
    print(json.dumps(per_eta, indent=2))


if __name__ == "__main__":
    main()
