"""Close the Track-F finite-mode ansatz class with exact certificates.

`START_NEW_SESSION_NAVIER_STOKES.md` §6 "優先候補A" proposes to *search* low-order
divergence-free Fourier ansätze for a residual force that stays smooth across a
finite-time singularity.  This experiment does not run that search: the search
space is provably empty (`docs/research_notes/track_f_finite_mode_nogo.md`,
Theorem 1 and Corollaries 1--2).  What it records instead is, for each candidate
ansatz family,

1. the **exact** verification of the trilinear cancellation `⟨u,(u·∇)u⟩ = 0`
   (integer arithmetic over `ℤ[i]`, every monomial coefficient checked), which
   is the one algebraic fact the no-go theorem rests on, and
2. the resulting rejection verdict together with the explicit constants of the
   a priori bound.

A floating-point cross-check integrates the Galerkin system to confirm that the
proved bound is respected at every step and is close to sharp, and two fault
injections confirm that both checks fire when the divergence-free hypothesis is
removed.  The cross-check is evidence about the implementation, never about the
theorem.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from ns_certificate_lab.galerkin_obstruction import (
    build_galerkin_system,
    build_mode_set,
    finite_mode_no_go_certificate,
    stream_apriori_bound,
    suggested_time_step,
)
from ns_certificate_lab.provenance import collect_runtime_provenance

EXPECTED_SCHEMA_VERSION = 1
EXPECTED_EXPERIMENT = "track_f_finite_mode_scan"
MINIMUM_FAMILIES = 3


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("config root must be an object")
    return data


def _plain_int(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    return int(value)


def _finite_number(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    converted = float(value)
    if not math.isfinite(converted):
        raise ValueError(f"{name} must be a finite number")
    return converted


def _validate_config(config: dict[str, Any]) -> None:
    """Reject unknown, missing, or mislabeled keys before any evidence is written."""
    required = {"schema_version", "experiment", "families", "cross_check"}
    if set(config) != required:
        raise ValueError("Track-F scan config has missing or unknown keys")
    if _plain_int(config["schema_version"], name="schema_version") != (
        EXPECTED_SCHEMA_VERSION
    ):
        raise ValueError("unsupported Track-F scan config schema")
    if config["experiment"] != EXPECTED_EXPERIMENT:
        raise ValueError("experiment must equal the audited canonical value")

    families = config["families"]
    if not isinstance(families, list) or len(families) < MINIMUM_FAMILIES:
        raise ValueError("families must list at least three ansatz families")
    labels: list[str] = []
    for index, family in enumerate(families):
        if not isinstance(family, dict) or "label" not in family:
            raise ValueError(f"families[{index}] must be an object with a label")
        label = family["label"]
        if not isinstance(label, str) or not label:
            raise ValueError(f"families[{index}].label must be a non-empty string")
        labels.append(label)
        keys = set(family)
        if keys == {"label", "seeds"}:
            seeds = family["seeds"]
            if not isinstance(seeds, list) or not seeds:
                raise ValueError(f"families[{index}].seeds must be a non-empty list")
            for position, seed in enumerate(seeds):
                if not isinstance(seed, list) or len(seed) != 3:
                    raise ValueError(
                        f"families[{index}].seeds[{position}] must have three entries"
                    )
                for component in seed:
                    _plain_int(component, name=f"families[{index}].seeds entry")
        elif keys == {"label", "shell_max_norm_squared"}:
            bound = _plain_int(
                family["shell_max_norm_squared"],
                name=f"families[{index}].shell_max_norm_squared",
            )
            if bound < 1:
                raise ValueError("shell_max_norm_squared must be at least 1")
        else:
            raise ValueError(
                f"families[{index}] must have exactly seeds or shell_max_norm_squared"
            )
    if len(set(labels)) != len(labels):
        raise ValueError("family labels must be unique")

    cross = config["cross_check"]
    required_cross = {
        "family",
        "seed",
        "viscosity",
        "inviscid_amplitude",
        "inviscid_steps",
        "inviscid_safety",
        "forced_magnitude",
        "forced_initial_norm",
        "forced_final_time",
        "forced_steps",
    }
    if not isinstance(cross, dict) or set(cross) != required_cross:
        raise ValueError("cross_check has missing or unknown keys")
    if cross["family"] not in labels:
        raise ValueError("cross_check.family must name one of the scanned families")
    _plain_int(cross["seed"], name="cross_check.seed")
    for key in ("inviscid_steps", "forced_steps"):
        if _plain_int(cross[key], name=f"cross_check.{key}") < 1:
            raise ValueError(f"cross_check.{key} must be positive")
    if _finite_number(cross["viscosity"], name="cross_check.viscosity") < 0.0:
        raise ValueError("cross_check.viscosity must be nonnegative")
    for key in (
        "inviscid_amplitude",
        "inviscid_safety",
        "forced_magnitude",
        "forced_initial_norm",
        "forced_final_time",
    ):
        if _finite_number(cross[key], name=f"cross_check.{key}") <= 0.0:
            raise ValueError(f"cross_check.{key} must be positive")


def _seeds_of(family: dict[str, Any]) -> list[list[int]]:
    if "seeds" in family:
        return [list(int(c) for c in seed) for seed in family["seeds"]]
    bound = int(family["shell_max_norm_squared"])
    reach = int(math.isqrt(bound))
    return [
        list(k)
        for k in itertools.product(range(-reach, reach + 1), repeat=3)
        if 1 <= sum(c * c for c in k) <= bound
    ]


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


def _cross_check(config: dict[str, Any], seeds: list[list[int]]) -> dict[str, Any]:
    cross = config["cross_check"]
    rng = np.random.default_rng(int(cross["seed"]))
    amplitude = float(cross["inviscid_amplitude"])
    safety = float(cross["inviscid_safety"])
    steps = int(cross["inviscid_steps"])

    # (a) inviscid, unforced, large amplitude: the would-be blow-up case.
    system = build_galerkin_system(build_mode_set(seeds), 0.0)
    state = rng.standard_normal(system.dimension)
    state *= amplitude / float(np.linalg.norm(state))
    dt = suggested_time_step(system, amplitude, safety=safety)
    inviscid = stream_apriori_bound(system, state, final_time=steps * dt, steps=steps)

    # (b) viscous, forced by a constant field aligned with a fixed direction:
    #     the trajectory should climb onto the proved bound without crossing it.
    viscous = build_galerkin_system(build_mode_set(seeds), float(cross["viscosity"]))
    direction = rng.standard_normal(viscous.dimension)
    direction /= float(np.linalg.norm(direction))
    magnitude = float(cross["forced_magnitude"])
    forced = stream_apriori_bound(
        viscous,
        float(cross["forced_initial_norm"]) * direction,
        final_time=float(cross["forced_final_time"]),
        steps=int(cross["forced_steps"]),
        force=lambda _t: magnitude * direction,
        force_l2_sup=magnitude,
    )

    # (c) fault injection: the same run with the divergence-free constraint
    #     removed must break the bound.
    faulty_system = build_galerkin_system(
        build_mode_set(seeds, allow_longitudinal=True), 0.0
    )
    faulty_state = rng.standard_normal(faulty_system.dimension)
    faulty_state *= amplitude / float(np.linalg.norm(faulty_state))
    faulty_dt = suggested_time_step(faulty_system, amplitude, safety=safety)
    faulted = stream_apriori_bound(
        faulty_system, faulty_state, final_time=steps * faulty_dt, steps=steps
    )

    return {
        "family": cross["family"],
        "inviscid_unforced": inviscid.as_dict(),
        "viscous_forced": forced.as_dict(),
        "fault_injection_longitudinal": faulted.as_dict(),
        "energy_production_detection_ratio": (
            faulted.max_relative_energy_production
            / max(inviscid.max_relative_energy_production, 1.0e-300)
        ),
    }


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = _load_json(config_path)
    _validate_config(config)
    if output_dir.exists():
        raise FileExistsError(
            f"refusing to overwrite existing output directory: {output_dir}"
        )
    provenance = collect_runtime_provenance()
    output_dir.mkdir(parents=True)

    rows: list[dict[str, Any]] = []
    cross_seeds: list[list[int]] | None = None
    for family in config["families"]:
        seeds = _seeds_of(family)
        certificate = finite_mode_no_go_certificate(family["label"], seeds)
        faulty = finite_mode_no_go_certificate(
            family["label"], seeds, allow_longitudinal=True
        )
        rows.append(
            {
                "label": family["label"],
                "mode_count": certificate.trilinear.mode_count,
                "dimension": certificate.dimension,
                "max_wavenumber": certificate.max_wavenumber,
                "variable_count": certificate.trilinear.variable_count,
                "resonant_triples": certificate.trilinear.resonant_triples,
                "monomials_accumulated": certificate.trilinear.monomials_accumulated,
                "surviving_monomials": certificate.trilinear.surviving_monomials,
                "transverse_defect": certificate.trilinear.transverse_defect,
                "sobolev_amplification_h1": certificate.sobolev_amplification_h1,
                "sobolev_amplification_h3": certificate.sobolev_amplification_h3,
                "gradient_amplification": certificate.gradient_amplification,
                "verdict": certificate.verdict,
                "fault_injection_surviving_monomials": (
                    faulty.trilinear.surviving_monomials
                ),
                "fault_injection_verdict": faulty.verdict,
            }
        )
        if family["label"] == config["cross_check"]["family"]:
            cross_seeds = seeds
    assert cross_seeds is not None  # guaranteed by _validate_config

    cross = _cross_check(config, cross_seeds)

    resonant_rows = [row for row in rows if int(row["resonant_triples"]) > 0]
    acceptance = {
        "every_family_cancels_exactly": all(
            int(row["surviving_monomials"]) == 0 and int(row["transverse_defect"]) == 0
            for row in rows
        ),
        "every_family_rejected": all(
            row["verdict"] == "rejected_as_clay_cd_candidate" for row in rows
        ),
        "certificate_not_vacuous": len(resonant_rows) >= MINIMUM_FAMILIES
        and all(int(row["monomials_accumulated"]) > 0 for row in resonant_rows),
        "fault_injection_breaks_every_resonant_family": all(
            int(row["fault_injection_surviving_monomials"]) > 0
            and row["fault_injection_verdict"] == "certificate_failed"
            for row in resonant_rows
        ),
        "inviscid_bound_respected": bool(cross["inviscid_unforced"]["bound_respected"]),
        "inviscid_energy_production_below_1e_12": float(
            cross["inviscid_unforced"]["max_relative_energy_production"]
        )
        < 1.0e-12,
        "forced_bound_respected": bool(cross["viscous_forced"]["bound_respected"]),
        "forced_bound_within_one_percent_of_sharp": float(
            cross["viscous_forced"]["max_bound_ratio"]
        )
        > 0.99,
        "fault_injection_breaks_the_streaming_bound": not bool(
            cross["fault_injection_longitudinal"]["bound_respected"]
        ),
        "fault_injection_detection_ratio_above_1e9": float(
            cross["energy_production_detection_ratio"]
        )
        > 1.0e9,
    }
    acceptance["all_passed"] = all(acceptance.values())

    summary = {
        "schema_version": 1,
        "experiment": EXPECTED_EXPERIMENT,
        "statement": (
            "exact rejection certificates for the Track-F finite-mode ansatz "
            "class; an exclusion theorem for a search class, not a singularity"
        ),
        "theorem": "docs/research_notes/track_f_finite_mode_nogo.md, Theorem 1",
        "formalization": (
            "formal/NSSingularity/GalerkinNoBlowup.lean (F-6): the a priori "
            "bound is proved in Lean 4; Lemma 1 and the ODE continuation are not"
        ),
        "rows": rows,
        "cross_check": cross,
        "acceptance": acceptance,
        "conclusion": (
            "No velocity field whose Fourier support stays inside a fixed finite "
            "set can be the breakdown solution of Clay (C)/(D).  The symbolic "
            "search proposed in START_NEW_SESSION_NAVIER_STOKES.md section 6 "
            "candidate A steps 1-5 is therefore provably empty and was not run."
        ),
        "limitations": [
            "the exact certificates cover only the listed finite mode sets; the "
            "general statement is the two-line proof of Lemma 1, not this scan",
            "the streaming cross-checks are binary64 without outward rounding",
            "the ODE continuation of Theorem 1(iii) is proved on paper only",
            "the Ladyzhenskaya-Prodi-Serrin consequence in section 5 of the note "
            "is a citation, neither proved here nor formalized",
            "nothing here bounds ansaetze whose Fourier bandwidth grows as t -> T",
        ],
        "reproducibility": {"runtime_provenance": provenance},
    }

    summary_path = output_dir / "summary.json"
    csv_path = output_dir / "certificates.csv"
    config_snapshot_path = output_dir / "config.snapshot.json"

    _write_json(summary_path, summary)
    _write_json(config_snapshot_path, config)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "label",
            "mode_count",
            "dimension",
            "max_wavenumber",
            "variable_count",
            "resonant_triples",
            "monomials_accumulated",
            "surviving_monomials",
            "transverse_defect",
            "sobolev_amplification_h1",
            "sobolev_amplification_h3",
            "gradient_amplification",
            "verdict",
            "fault_injection_surviving_monomials",
            "fault_injection_verdict",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})

    manifest = {
        "schema_version": 1,
        "payloads": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in (summary_path, csv_path, config_snapshot_path)
        },
    }
    manifest_path = output_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    (output_dir / "manifest.json.sha256").write_text(
        _sha256(manifest_path) + "\n", encoding="ascii"
    )

    if not acceptance["all_passed"]:
        raise RuntimeError("Track-F finite-mode scan failed acceptance checks")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.config, args.output_dir)
    print(json.dumps(summary["acceptance"], ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
