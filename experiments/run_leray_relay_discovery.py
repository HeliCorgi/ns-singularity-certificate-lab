"""Run the CPU-small exact-Leray relay discovery pilot.

The pilot has three deliberately separated layers:

1. exact rational cos/sin algebra for a true three-mode Navier--Stokes relay;
2. exact rational checks of the modal-growth covariance identity;
3. floating, dealiased Fejer-packet response-map searches at increasing scales.

The first layer proves a signed one-step transfer identity.  The third layer is
only a falsification screen for scale-uniform coherence.  No PDE orbit,
singularity, smooth-force construction, or continuation theorem is inferred.
"""

from __future__ import annotations

import argparse
import csv
import gc
import io
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    digest_sidecar,
    sha256_file,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)
from ns_certificate_lab.exact_leray_relay import (
    build_exact_relay_triad,
    exact_relay_metrics,
    fixed_cardinality_scaling,
)
from ns_certificate_lab.fourier_torus import advection, leray
from ns_certificate_lab.leray_response_relay import (
    fejer_carrier_packet,
    harmonic_carrier_mask,
    relay_stage,
)
from ns_certificate_lab.modal_front_actions import (
    h3_bandwidth_factorization,
    modal_growth_identity,
)
from ns_certificate_lab.provenance import (
    collect_runtime_provenance,
    validate_runtime_provenance,
)


CONFIG_SCHEMA = "ns-certificate-lab/leray-relay-discovery-config/v1"
OUTPUT_SCHEMA = "ns-certificate-lab/leray-relay-discovery/v1"
MANIFEST_SCHEMA = "ns-certificate-lab/leray-relay-discovery-manifest/v1"
STATUS = "SYMBOLIC CANDIDATE / AUDIT REQUIRED"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PAYLOAD_NAMES = (
    "config.snapshot.json",
    "exact_triad.json",
    "modal_actions.json",
    "response_search.csv",
    "response_ladder.csv",
    "summary.json",
)
SOURCE_PATHS = (
    "src/ns_certificate_lab/exact_leray_relay.py",
    "src/ns_certificate_lab/leray_response_relay.py",
    "src/ns_certificate_lab/modal_front_actions.py",
    "src/ns_certificate_lab/fourier_torus.py",
    "src/ns_certificate_lab/torus_chain.py",
    "src/ns_certificate_lab/_integrity.py",
    "src/ns_certificate_lab/provenance.py",
    "experiments/run_leray_relay_discovery.py",
    "tests/test_exact_leray_relay.py",
    "tests/test_leray_response_relay.py",
    "tests/test_leray_relay_experiment.py",
)


def _fraction(value: object, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise ValueError(f"{name} must be an integer or rational string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"{name} is not a valid rational") from error
    return result


def _load_config(path: Path) -> dict[str, Any]:
    value = strict_json_loads(path.read_text(encoding="utf-8"), label="relay config")
    required = {
        "schema",
        "seed",
        "viscosity",
        "exact_scale",
        "exact_parent_sine",
        "exact_parent_cosine",
        "exact_child_cosine",
        "search_grid_size",
        "recursive_grid_size",
        "carrier",
        "envelope",
        "random_polarization_samples",
        "chirp_values",
        "recursive_harmonics",
        "baseline_polarization_real",
        "baseline_polarization_imag",
    }
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("relay config has missing or unknown fields")
    if value["schema"] != CONFIG_SCHEMA:
        raise ValueError("relay config schema is invalid")
    for name in (
        "seed",
        "exact_scale",
        "search_grid_size",
        "recursive_grid_size",
        "carrier",
        "envelope",
        "random_polarization_samples",
    ):
        item = value[name]
        if isinstance(item, bool) or not isinstance(item, int) or item < 1:
            raise ValueError(f"{name} must be a positive integer")
    viscosity = _fraction(value["viscosity"], name="viscosity")
    if viscosity <= 0:
        raise ValueError("viscosity must be positive")
    exact_coefficients = {
        name: _fraction(value[name], name=name)
        for name in (
            "exact_parent_sine",
            "exact_parent_cosine",
            "exact_child_cosine",
        )
    }
    signed_product = (
        exact_coefficients["exact_parent_sine"]
        * exact_coefficients["exact_parent_cosine"]
        * exact_coefficients["exact_child_cosine"]
    )
    if signed_product <= 0:
        raise ValueError("exact triad coefficients must give positive B*C*D")
    scale = Fraction(value["exact_scale"])
    child = exact_coefficients["exact_child_cosine"]
    child_margin = (
        scale * signed_product / 2
        - 9 * viscosity * scale * scale * child * child
    )
    if child_margin <= 0:
        raise ValueError("exact triad must have a positive post-viscous child margin")
    chirps = value["chirp_values"]
    if (
        not isinstance(chirps, list)
        or not chirps
        or any(
            isinstance(item, bool)
            or not isinstance(item, (int, float))
            or not np.isfinite(float(item))
            for item in chirps
        )
    ):
        raise ValueError("chirp_values must be a nonempty finite numeric list")
    harmonics = value["recursive_harmonics"]
    if (
        not isinstance(harmonics, list)
        or not harmonics
        or any(
            isinstance(item, bool) or not isinstance(item, int) or item < 2
            for item in harmonics
        )
        or harmonics != [2 ** (index + 1) for index in range(len(harmonics))]
    ):
        raise ValueError("recursive_harmonics must be [2,4,...]")
    for name in ("baseline_polarization_real", "baseline_polarization_imag"):
        vector = value[name]
        if (
            not isinstance(vector, list)
            or len(vector) != 3
            or any(
                isinstance(item, bool)
                or not isinstance(item, (int, float))
                or not np.isfinite(float(item))
                for item in vector
            )
        ):
            raise ValueError(f"{name} must be a finite three-vector")
    baseline_norm_squared = sum(
        float(real) ** 2 + float(imag) ** 2
        for real, imag in zip(
            value["baseline_polarization_real"],
            value["baseline_polarization_imag"],
        )
    )
    if baseline_norm_squared == 0.0:
        raise ValueError("baseline polarization must be nonzero")
    carrier = value["carrier"]
    envelope = value["envelope"]
    if 2 * (carrier + envelope - 1) >= value["search_grid_size"]:
        raise ValueError("search parent support reaches the Nyquist mode")
    largest_harmonic = harmonics[-1]
    largest_edge = largest_harmonic * (carrier + envelope - 1)
    if largest_edge >= value["recursive_grid_size"] // 2:
        raise ValueError("recursive child support reaches the Nyquist mode")
    return value


def _csv_bytes(rows: list[dict[str, Any]]) -> bytes:
    if not rows:
        return b""
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _field_table(field) -> dict[str, dict[str, list[str]]]:
    output: dict[str, dict[str, list[str]]] = {}
    for wave, (cosine, sine) in field.coefficient_table().items():
        output[",".join(str(component) for component in wave)] = {
            "cosine": [str(value) for value in cosine],
            "sine": [str(value) for value in sine],
        }
    return output


def _row(stage, *, candidate: int, chirp: float, level: int) -> dict[str, Any]:
    metrics = stage.metrics
    injection = float(metrics["normalized_injection_per_scale"])
    dissipation = float(metrics["normalized_dissipation_per_scale"])
    return {
        "candidate": candidate,
        "chirp": chirp,
        "level": level,
        "parent_scale": metrics["parent_scale"],
        "parent_support_min_wavenumber": metrics[
            "parent_support_min_wavenumber"
        ],
        "parent_support_max_wavenumber": metrics[
            "parent_support_max_wavenumber"
        ],
        "parent_rms_wavenumber": metrics["parent_rms_wavenumber"],
        "parent_scale_to_rms_ratio": metrics["parent_scale_to_rms_ratio"],
        "parent_energy": metrics["parent_energy"],
        "child_energy": metrics["child_energy"],
        "normalized_injection": injection,
        "normalized_dissipation": dissipation,
        "injection_to_dissipation_shape": injection / dissipation,
        "critical_energy_constant": metrics["critical_energy_constant"],
        "child_forcing_fraction": metrics["child_forcing_fraction"],
        "off_chain_forcing_ratio": metrics["off_chain_forcing_ratio"],
        "parabolic_fill_constant": metrics["parabolic_fill_constant"],
        "best_response_alignment": metrics["best_response_alignment"],
        "forcing_resolution_floor": metrics["forcing_resolution_floor"],
        "forcing_to_resolution_floor": metrics[
            "forcing_to_resolution_floor"
        ],
        "full_populated_child_flux": metrics["full_populated_child_flux"],
        "full_nonlinear_energy_defect": metrics["full_nonlinear_energy_defect"],
        "full_populated_flux_tested": metrics["full_populated_flux_tested"],
    }


def _simple_response_orbit_passes(rows: list[dict[str, Any]]) -> bool:
    """Apply the hard-coded discovery screen, including every full flux."""

    if not rows:
        return False
    injections = [float(row["normalized_injection"]) for row in rows]
    off_ratios = [float(row["off_chain_forcing_ratio"]) for row in rows]
    full_fluxes = [row["full_populated_child_flux"] for row in rows]
    full_flux_tested = all(
        bool(row["full_populated_flux_tested"]) for row in rows
    )
    full_flux_positive = full_flux_tested and all(
        value is not None and float(value) > 0.0 for value in full_fluxes
    )
    scale_uniformity_ratio = min(injections) / max(injections)
    return (
        scale_uniformity_ratio >= 0.5
        and max(off_ratios) <= 1.0
        and full_flux_positive
    )


def _search_response(config: dict[str, Any], viscosity: float):
    rng = np.random.default_rng(int(config["seed"]))
    chirps = [float(value) for value in config["chirp_values"]]
    baseline = np.asarray(config["baseline_polarization_real"], dtype=np.float64)
    baseline = baseline + 1.0j * np.asarray(
        config["baseline_polarization_imag"], dtype=np.float64
    )
    baseline /= np.linalg.norm(baseline)
    candidates: list[tuple[complex, ...] | np.ndarray] = [baseline for _ in chirps]
    candidate_chirps = list(chirps)
    for index in range(int(config["random_polarization_samples"])):
        vector = rng.normal(size=3) + 1.0j * rng.normal(size=3)
        vector /= np.linalg.norm(vector)
        candidates.append(vector)
        candidate_chirps.append(chirps[index % len(chirps)])

    grid = int(config["search_grid_size"])
    carrier = int(config["carrier"])
    envelope = int(config["envelope"])
    mask_one = harmonic_carrier_mask(
        grid, carrier=carrier, envelope=envelope, harmonic=2
    )
    mask_two = harmonic_carrier_mask(
        grid, carrier=carrier, envelope=envelope, harmonic=4
    )
    rows: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    for index, (polarization, chirp) in enumerate(zip(candidates, candidate_chirps)):
        parent = fejer_carrier_packet(
            grid,
            carrier=carrier,
            envelope=envelope,
            polarization=polarization,
            chirp=chirp,
        )
        first = relay_stage(
            parent, mask_one, parent_scale=float(carrier), viscosity=viscosity
        )
        second = relay_stage(
            first.child,
            mask_two,
            parent_scale=float(2 * carrier),
            viscosity=viscosity,
        )
        first_row = _row(first, candidate=index, chirp=chirp, level=0)
        second_row = _row(second, candidate=index, chirp=chirp, level=1)
        rows.extend((first_row, second_row))
        score = min(
            float(first_row["injection_to_dissipation_shape"]),
            float(second_row["injection_to_dissipation_shape"]),
        )
        record = {
            "candidate": index,
            "chirp": chirp,
            "score": score,
            "polarization_real": np.asarray(polarization).real.tolist(),
            "polarization_imag": np.asarray(polarization).imag.tolist(),
            "level0": first_row,
            "level1": second_row,
        }
        if best is None or score > float(best["score"]):
            best = record
        del parent, first, second
        gc.collect()
    if best is None:
        raise AssertionError("response search produced no candidates")
    return rows, best


def _recursive_ladder(
    config: dict[str, Any], viscosity: float, best: dict[str, Any]
) -> list[dict[str, Any]]:
    grid = int(config["recursive_grid_size"])
    carrier = int(config["carrier"])
    envelope = int(config["envelope"])
    polarization = np.asarray(best["polarization_real"], dtype=np.float64)
    polarization = polarization + 1.0j * np.asarray(
        best["polarization_imag"], dtype=np.float64
    )
    parent = fejer_carrier_packet(
        grid,
        carrier=carrier,
        envelope=envelope,
        polarization=polarization,
        chirp=float(best["chirp"]),
    )
    rows: list[dict[str, Any]] = []
    for level, harmonic in enumerate(config["recursive_harmonics"]):
        mask = harmonic_carrier_mask(
            grid,
            carrier=carrier,
            envelope=envelope,
            harmonic=int(harmonic),
        )
        stage = relay_stage(
            parent,
            mask,
            parent_scale=float(carrier * 2**level),
            viscosity=viscosity,
        )
        rows.append(
            _row(
                stage,
                candidate=int(best["candidate"]),
                chirp=float(best["chirp"]),
                level=level,
            )
        )
        parent = stage.child
        del stage
        gc.collect()
    return rows


def verify_leray_relay_bundle(output_dir: Path) -> dict[str, Any]:
    """Independently verify the payload set, sidecars, and manifest."""

    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        raise ValueError("bundle directory is missing")
    data_names = (*PAYLOAD_NAMES, "manifest.json")
    expected_names = set(data_names)
    expected_names.update(digest_sidecar(Path(name)).name for name in data_names)
    observed_names = {path.name for path in output_dir.iterdir()}
    if observed_names != expected_names:
        raise ValueError("bundle file set is not exact")
    for name in data_names:
        verify_digest(output_dir / name)
    manifest = strict_json_loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8"),
        label="relay manifest",
    )
    if not isinstance(manifest, dict) or set(manifest) != {
        "schema",
        "status",
        "files",
        "seed",
        "source_config",
        "source_files",
    }:
        raise ValueError("manifest has missing or unknown fields")
    if manifest["schema"] != MANIFEST_SCHEMA or manifest["status"] != STATUS:
        raise ValueError("manifest schema or status is invalid")
    records = manifest["files"]
    if not isinstance(records, dict) or set(records) != set(PAYLOAD_NAMES):
        raise ValueError("manifest payload inventory is invalid")
    for name, record in records.items():
        if not isinstance(record, dict) or set(record) != {"bytes", "sha256"}:
            raise ValueError("manifest file record is invalid")
        path = output_dir / name
        if record["bytes"] != path.stat().st_size:
            raise ValueError(f"byte count mismatch for {name}")
        if record["sha256"] != sha256_file(path):
            raise ValueError(f"manifest checksum mismatch for {name}")
    source_files = manifest["source_files"]
    if not isinstance(source_files, dict) or set(source_files) != set(SOURCE_PATHS):
        raise ValueError("manifest source inventory is invalid")
    for name, digest in source_files.items():
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise ValueError(f"invalid source checksum for {name}")
        source_path = REPOSITORY_ROOT / name
        if not source_path.is_file() or sha256_file(source_path) != digest:
            raise ValueError(f"source checksum mismatch for {name}")
    config = strict_json_loads(
        (output_dir / "config.snapshot.json").read_text(encoding="utf-8"),
        label="relay config snapshot",
    )
    summary = strict_json_loads(
        (output_dir / "summary.json").read_text(encoding="utf-8"),
        label="relay summary",
    )
    if config.get("schema") != CONFIG_SCHEMA or summary.get("schema") != OUTPUT_SCHEMA:
        raise ValueError("bundle schema is invalid")
    if summary.get("status") != manifest["status"]:
        raise ValueError("summary status does not match manifest")
    validate_runtime_provenance(summary.get("provenance"))
    if manifest["seed"] != config["seed"]:
        raise ValueError("manifest seed does not match config snapshot")
    return {
        "verified": True,
        "schema": OUTPUT_SCHEMA,
        "payload_count": len(PAYLOAD_NAMES),
        "source_file_count": len(SOURCE_PATHS),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/leray_relay_discovery_v1.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/leray_relay_discovery_v1"),
    )
    args = parser.parse_args()
    config = _load_config(args.config)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("output directory must be new or empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    viscosity_fraction = _fraction(config["viscosity"], name="viscosity")
    viscosity = float(viscosity_fraction)
    exact_field = build_exact_relay_triad(
        scale=int(config["exact_scale"]),
        parent_sine=_fraction(config["exact_parent_sine"], name="parent_sine"),
        parent_cosine=_fraction(config["exact_parent_cosine"], name="parent_cosine"),
        child_cosine=_fraction(config["exact_child_cosine"], name="child_cosine"),
    )
    exact_metrics = exact_relay_metrics(
        viscosity=viscosity_fraction,
        scale=int(config["exact_scale"]),
        parent_sine=_fraction(config["exact_parent_sine"], name="parent_sine"),
        parent_cosine=_fraction(config["exact_parent_cosine"], name="parent_cosine"),
        child_cosine=_fraction(config["exact_child_cosine"], name="child_cosine"),
    )
    exact_nonlinear = leray(advection(exact_field, exact_field)).cleaned()
    exact_payload = {
        "schema": "ns-certificate-lab/exact-leray-triad/v1",
        "status": "EXACT FINITE ALGEBRA / NOT A PDE ORBIT",
        "claim_boundary": (
            "Only the finite coefficient table, transfer, viscosity, and "
            "off-graph identities are exact. Multiplicity and phase scaling "
            "below are analytic bounds or explicitly labelled heuristics."
        ),
        "field_coefficients": _field_table(exact_field),
        "leray_nonlinear_coefficients": _field_table(exact_nonlinear),
        "metrics": exact_metrics.as_dict(),
        "fixed_cardinality_scaling": fixed_cardinality_scaling(),
        "multiplicity_necessity": {
            "bound": "|Pi_N| <= C N sqrt(M_N) E_N^(3/2)",
            "critical_energy_law": "E_N ~ c_E/N",
            "viscous_loss": "nu N^2 E_N ~ nu c_E N",
            "necessary_mode_count": "M_N >= constant * (nu^2/c_E) * N^3",
        },
        "heuristic_phase_scaling": {
            "status": "UNVERIFIED HEURISTIC / NOT AN EXACT PAYLOAD CLAIM",
            "random_phase_flux": "N^-1/2 under independent-phase scaling",
            "coherent_target_flux": "N^1 at Bernstein-capacity saturation",
            "required_gain_over_random_phase": "N^(3/2)",
        },
    }
    modal_payload = {
        "schema": "ns-certificate-lab/modal-growth-actions/v1",
        "identities": [
            modal_growth_identity(
                exact_field, order=order, viscosity=viscosity_fraction
            ).as_dict()
            for order in range(3)
        ],
        "h3_factorization": h3_bandwidth_factorization(exact_field),
        "monotone_quantity": (
            "M_r(t)=log N_r(t)-integral sigma_r^2/(4 nu N_r^2) dt "
            "has nonpositive derivative"
        ),
    }

    search_rows, best = _search_response(config, viscosity)
    ladder_rows = _recursive_ladder(config, viscosity, best)
    injections = [float(row["normalized_injection"]) for row in ladder_rows]
    off_ratios = [float(row["off_chain_forcing_ratio"]) for row in ladder_rows]
    critical_constants = [float(row["critical_energy_constant"]) for row in ladder_rows]
    scale_uniformity_ratio = min(injections) / max(injections)
    full_flux_tested = all(
        bool(row["full_populated_flux_tested"]) for row in ladder_rows
    )
    full_flux_positive = full_flux_tested and all(
        row["full_populated_child_flux"] is not None
        and float(row["full_populated_child_flux"]) > 0.0
        for row in ladder_rows
    )
    simple_orbit_passes = _simple_response_orbit_passes(ladder_rows)
    summary = {
        "schema": OUTPUT_SCHEMA,
        "status": STATUS,
        "claim_boundary": (
            "The exact rational three-mode calculation proves a positive signed "
            "one-step transfer for the true periodic Leray nonlinearity. The "
            "finite-cardinality scale iteration is analytically rejected. The "
            "dealiased Fejer response map tests mode multiplicity and three "
            "algebraically aligned stages, but constructs no invariant phase cone, "
            "interval budget, infinite-band PDE orbit, singularity, or theorem."
        ),
        "exact_triad": exact_metrics.as_dict(),
        "modal_action_checks": {
            "orders": [0, 1, 2],
            "covariance_identities_exact": True,
            "cauchy_schwarz_gaps_nonnegative": True,
            "h3_factorization_exact": True,
        },
        "response_search": {
            "candidate_count": len(search_rows) // 2,
            "best": best,
        },
        "recursive_ladder": {
            "levels": len(ladder_rows),
            "normalized_injections": injections,
            "off_chain_forcing_ratios": off_ratios,
            "critical_energy_constants": critical_constants,
            "last_to_first_injection_ratio": injections[-1] / injections[0],
            "scale_uniformity_ratio": scale_uniformity_ratio,
            "all_full_populated_fluxes_tested": full_flux_tested,
            "all_full_populated_fluxes_positive": full_flux_positive,
            "simple_response_orbit_passes_hard_coded_screen": simple_orbit_passes,
            "verdict": (
                "SURVIVES the simple scale-uniformity/off-chain screen"
                if simple_orbit_passes
                else "REJECTED as a simple self-iterated response orbit; retain the exact triad and search phase-coded N^3 mode clouds"
            ),
        },
        "shortest_kill_condition": (
            "Reject a proposed coherent cloud if its normalized signed flux "
            "chi_N tends to zero over two scale doublings, if fewer than order "
            "N^3 modes are active under E_N~N^-1, or if low/off-chain energy "
            "reaches one half of child gain before the critical child target."
        ),
        "provenance": collect_runtime_provenance(),
    }

    files = {
        "config.snapshot.json": canonical_json_bytes(config),
        "exact_triad.json": canonical_json_bytes(exact_payload),
        "modal_actions.json": canonical_json_bytes(modal_payload),
        "response_search.csv": _csv_bytes(search_rows),
        "response_ladder.csv": _csv_bytes(ladder_rows),
        "summary.json": canonical_json_bytes(summary),
    }
    for name, payload in files.items():
        write_with_digest(args.output_dir / name, payload)
    source_files = {
        name: sha256_file(REPOSITORY_ROOT / name) for name in SOURCE_PATHS
    }
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "status": STATUS,
        "files": {
            name: {
                "bytes": (args.output_dir / name).stat().st_size,
                "sha256": sha256_file(args.output_dir / name),
            }
            for name in sorted(files)
        },
        "seed": int(config["seed"]),
        "source_config": args.config.as_posix(),
        "source_files": source_files,
    }
    write_with_digest(args.output_dir / "manifest.json", canonical_json_bytes(manifest))
    verify_leray_relay_bundle(args.output_dir)
    print(canonical_json_bytes(summary).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
