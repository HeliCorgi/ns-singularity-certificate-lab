"""Run the CPU-small Zeno packet-relay falsification pilot.

The experiment checks only three algebraic/kinematic premises:

1. a finite-low-frequency dyadic staircase can have bounded energy while its
   disjoint-packet ``L^3`` mass grows with the number of occupied scales;
2. a localized anisotropic transmitter has a far-field pressure Hessian with a
   favourable stretching direction and the required distance ``^-5`` law;
3. a scale-local net flux proportional to ``N`` leads to the formal moving-front
   law ``dN/dt = c N^3`` and hence ``N ~ (T-t)^-1/2``.

None of these checks establishes the missing Navier--Stokes cascade cell.  The
result is a FORMAL ANSATZ audit artifact, never a singularity certificate.
"""

from __future__ import annotations

import argparse
import csv
import io
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
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.zeno_packet_relay import (
    MovingFrontLaw,
    classify_finite_floor_shell_exponents,
    critical_staircase_metrics,
    gaussian_vortex_moment,
    normalized_pressure_stretch,
    pressure_hessian_quadrupole,
    zeno_scaling_exponents,
)


CONFIG_SCHEMA = "ns-certificate-lab/zeno-packet-relay-pilot-config/v2"
OUTPUT_SCHEMA = "ns-certificate-lab/zeno-packet-relay-pilot/v2"
MANIFEST_SCHEMA = "ns-certificate-lab/zeno-packet-relay-manifest/v2"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PAYLOAD_NAMES = (
    "config.snapshot.json",
    "stack_scaling.csv",
    "pressure_distance_scaling.csv",
    "front_scaling.csv",
    "summary.json",
)
SOURCE_PATHS = (
    "src/ns_certificate_lab/zeno_packet_relay.py",
    "src/ns_certificate_lab/_integrity.py",
    "src/ns_certificate_lab/provenance.py",
    "experiments/run_zeno_packet_relay_pilot.py",
    "tests/test_zeno_packet_relay.py",
    "tests/test_zeno_packet_relay_experiment.py",
)


def _load_config(path: Path) -> dict[str, Any]:
    value = strict_json_loads(path.read_text(encoding="utf-8"), label="pilot config")
    required = {
        "schema",
        "seed",
        "levels",
        "base",
        "bandwidth_exponent",
        "peak_energy_exponent",
        "spectral_slope",
        "pressure_moment_model",
        "random_orientation_samples",
        "front_initial_bandwidth",
        "front_coefficient",
        "front_remaining_fractions",
    }
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("pilot config has missing or unknown fields")
    if value["schema"] != CONFIG_SCHEMA:
        raise ValueError("pilot config schema is invalid")
    if value["pressure_moment_model"] != "gaussian_vortex_schwartz":
        raise ValueError("pressure_moment_model is invalid")
    levels = value["levels"]
    if (
        not isinstance(levels, list)
        or not levels
        or any(isinstance(item, bool) or not isinstance(item, int) or item < 0 for item in levels)
        or levels != sorted(set(levels))
    ):
        raise ValueError("levels must be a sorted unique list of nonnegative integers")
    if (
        isinstance(value["seed"], bool)
        or not isinstance(value["seed"], int)
        or value["seed"] < 0
    ):
        raise ValueError("seed must be a nonnegative integer")
    if (
        isinstance(value["random_orientation_samples"], bool)
        or not isinstance(value["random_orientation_samples"], int)
        or value["random_orientation_samples"] < 1
    ):
        raise ValueError("random_orientation_samples must be positive")
    fractions = value["front_remaining_fractions"]
    if (
        not isinstance(fractions, list)
        or not fractions
        or any(
            isinstance(item, bool)
            or not isinstance(item, (int, float))
            or not 0.0 < float(item) < 1.0
            for item in fractions
        )
    ):
        raise ValueError("front_remaining_fractions must lie strictly in (0,1)")
    return value


def _csv_bytes(rows: list[dict[str, Any]]) -> bytes:
    if not rows:
        return b""
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _orientation_search(samples: int, seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    best = -np.inf
    best_vortex_axis = np.zeros(3, dtype=np.float64)
    best_separation = np.zeros(3, dtype=np.float64)
    best_direction = np.zeros(3, dtype=np.float64)

    for _ in range(samples):
        vortex_axis = rng.normal(size=3)
        vortex_axis /= np.linalg.norm(vortex_axis)
        separation = rng.normal(size=3)
        separation /= np.linalg.norm(separation)
        moment = gaussian_vortex_moment(vortex_axis)
        hessian = pressure_hessian_quadrupole(moment, separation)
        eigenvalues, eigenvectors = np.linalg.eigh(-4.0 * np.pi * hessian)
        index = int(np.argmax(eigenvalues))
        value = float(eigenvalues[index])
        if value > best:
            best = value
            best_vortex_axis = vortex_axis.copy()
            best_separation = separation.copy()
            best_direction = eigenvectors[:, index].copy()

    return {
        "sample_count": samples,
        "best_normalized_stretch": best,
        "best_vortex_axis": best_vortex_axis.tolist(),
        "best_separation_axis": best_separation.tolist(),
        "best_stretch_axis": best_direction.tolist(),
    }


def verify_zeno_packet_relay_bundle(output_dir: Path) -> dict[str, Any]:
    """Independently verify the v2 payload set, sidecars, and manifest."""

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
        label="Zeno manifest",
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
    if manifest["schema"] != MANIFEST_SCHEMA:
        raise ValueError("manifest schema is invalid")
    if manifest["status"] != "FORMAL ANSATZ / AUDIT REQUIRED":
        raise ValueError("manifest status is invalid")
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

    config = strict_json_loads(
        (output_dir / "config.snapshot.json").read_text(encoding="utf-8"),
        label="Zeno config snapshot",
    )
    if config["schema"] != CONFIG_SCHEMA:
        raise ValueError("config snapshot schema is invalid")
    summary = strict_json_loads(
        (output_dir / "summary.json").read_text(encoding="utf-8"),
        label="Zeno summary",
    )
    if summary.get("schema") != OUTPUT_SCHEMA:
        raise ValueError("summary schema is invalid")
    if summary.get("status") != manifest["status"]:
        raise ValueError("summary status does not match manifest")
    provenance = summary.get("provenance")
    if (
        not isinstance(provenance, dict)
        or "source_fingerprint_sha256" not in provenance
    ):
        raise ValueError("summary provenance is incomplete")
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
        default=Path("configs/zeno_packet_relay_pilot_v2.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/zeno_packet_relay_pilot_v2"),
    )
    args = parser.parse_args()
    config = _load_config(args.config)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("output directory must be new or empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    base = float(config["base"])
    stack_rows: list[dict[str, Any]] = []
    for level in config["levels"]:
        metrics = critical_staircase_metrics(level, base=base)
        row = metrics.as_dict()
        row["energy_gap_to_geometric_limit"] = base / (base - 1.0) - metrics.energy
        row["critical_besov_per_shell"] = metrics.critical_besov / metrics.shell_count
        row["enstrophy_over_bandwidth"] = metrics.enstrophy / metrics.bandwidth
        stack_rows.append(row)

    # W_a(x)=(a cross x) exp(-|x|^2/2), a=e1, is an explicit
    # divergence-free Schwartz transmitter with this rank-two moment.
    vortex_axis = np.array([1.0, 0.0, 0.0])
    moment = gaussian_vortex_moment(vortex_axis)
    direction = vortex_axis
    distance_rows: list[dict[str, Any]] = []
    for distance in (2.0, 4.0, 8.0, 16.0):
        displacement = np.array([distance, 0.0, 0.0])
        hessian = pressure_hessian_quadrupole(moment, displacement)
        distance_rows.append(
            {
                "distance": distance,
                "frobenius_norm": float(np.linalg.norm(hessian)),
                "distance5_times_frobenius": float(
                    distance**5 * np.linalg.norm(hessian)
                ),
                "normalized_transverse_stretch": normalized_pressure_stretch(
                    moment, displacement, direction
                ),
                "trace": float(np.trace(hessian)),
            }
        )

    orientation = _orientation_search(
        int(config["random_orientation_samples"]), int(config["seed"])
    )
    verdict = classify_finite_floor_shell_exponents(
        float(config["bandwidth_exponent"]),
        float(config["peak_energy_exponent"]),
        float(config["spectral_slope"]),
    )

    law = MovingFrontLaw(
        initial_bandwidth=float(config["front_initial_bandwidth"]),
        coefficient=float(config["front_coefficient"]),
    )
    front_rows: list[dict[str, Any]] = []
    full_interval = law.blowup_time - law.initial_time
    for fraction in config["front_remaining_fractions"]:
        remaining = float(fraction) * full_interval
        time = law.blowup_time - remaining
        bandwidth = law.bandwidth(time)
        front_rows.append(
            {
                "time": time,
                "remaining_time": remaining,
                "bandwidth": bandwidth,
                "bandwidth_times_sqrt_remaining": bandwidth * np.sqrt(remaining),
            }
        )

    distance_invariant = [row["distance5_times_frobenius"] for row in distance_rows]
    summary = {
        "schema": OUTPUT_SCHEMA,
        "status": "FORMAL ANSATZ / AUDIT REQUIRED",
        "claim_boundary": (
            "Separate kinematic packet/shell sums, the leading quadrupole tensor "
            "of an explicit Gaussian divergence-free vortex, and a formal front "
            "ODE were checked. The pressure homogeneity check is tautological at "
            "leading order: no exact packet pressure, multipole remainder, signed "
            "Navier-Stokes flux, orbit, singularity, or theorem was constructed."
        ),
        "finite_floor_verdict": verdict.as_dict(),
        "zeno_scaling": zeno_scaling_exponents(
            float(config["bandwidth_exponent"])
        ),
        "stack_checks": {
            "max_energy": max(float(row["energy"]) for row in stack_rows),
            "critical_mass_grows_with_shell_count": all(
                float(right["critical_l3_cubed"]) > float(left["critical_l3_cubed"])
                for left, right in zip(stack_rows, stack_rows[1:])
            ),
            "besov_per_shell_max_error": max(
                abs(float(row["critical_besov_per_shell"]) - 1.0)
                for row in stack_rows
            ),
        },
        "pressure_checks": {
            "moment_model": config["pressure_moment_model"],
            "explicit_field": "W_a(x)=(a cross x) exp(-|x|^2/2)",
            "moment_eigenvalues": np.linalg.eigvalsh(moment).tolist(),
            "analytic_favourable_stretch": 12.0,
            "leading_tensor_homogeneity_relative_spread": (
                max(distance_invariant) - min(distance_invariant)
            )
            / max(distance_invariant),
            "rank_one_v1_algebraic_optimum_rejected_as_unrealizable": True,
            "exact_pressure_and_multipole_remainder_tested": False,
            "separation_remainder_versus_strength_tradeoff_tested": False,
            "orientation_search": orientation,
        },
        "front_checks": {
            "ode": "dN/dt = coefficient * N^3",
            "formal_blowup_time": law.blowup_time,
            "parabolic_invariant_relative_spread": (
                max(float(row["bandwidth_times_sqrt_remaining"]) for row in front_rows)
                - min(float(row["bandwidth_times_sqrt_remaining"]) for row in front_rows)
            )
            / max(float(row["bandwidth_times_sqrt_remaining"]) for row in front_rows),
        },
        "kill_condition": (
            "Reject the relay mechanism if a true divergence-free Navier-Stokes "
            "cascade cell cannot produce positive scale-normalized child flux after "
            "viscous loss while keeping off-chain leakage summable. This pilot does "
            "not test that decisive condition."
        ),
        "provenance": collect_runtime_provenance(),
    }

    files = {
        "config.snapshot.json": canonical_json_bytes(config),
        "stack_scaling.csv": _csv_bytes(stack_rows),
        "pressure_distance_scaling.csv": _csv_bytes(distance_rows),
        "front_scaling.csv": _csv_bytes(front_rows),
        "summary.json": canonical_json_bytes(summary),
    }
    for name, payload in files.items():
        write_with_digest(args.output_dir / name, payload)

    source_files = {
        name: sha256_file(REPOSITORY_ROOT / name) for name in SOURCE_PATHS
    }
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "status": "FORMAL ANSATZ / AUDIT REQUIRED",
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
    verify_zeno_packet_relay_bundle(args.output_dir)
    print(canonical_json_bytes(summary).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
