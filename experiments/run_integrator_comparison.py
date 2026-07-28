r"""Gate 1 (P0-A) cross-integrator comparison experiment.

Runs the E-29 early-time problem with ``heun``, ``ssprk3`` and ``rk4`` on the
*identical* spatial discretization and the *identical* fixed time-step series,
then compares amplification, peak location, ``max |u_1|`` and the final fields
pairwise.  Rationale: Heun with centered differences is von-Neumann unstable
for pure advection (``|G(i a)|^2 = 1 + a^4/4``), so amplification measured
with Heun alone must never feed a candidate decision; SSPRK3 and RK4 have
imaginary-axis stability intervals and serve as cross-checks.  If the schemes
disagree beyond the preregistered tolerances the run FAILS and that is a
Gate 1 stop condition to be reported, not repaired by tuning.

Amplification is reported under BOTH normalizations of P1-B:

* ``amplification_grid_normalized``: final over *discrete initial* maximum
  Cartesian vorticity on the same grid (resolution-dependent denominator);
* ``amplification_common_reference``: final over the continuum E-29b value
  ``24000 pi 37^{-1/2} (36/37)^{18}``, evaluated in code from that closed
  form, shared by every resolution.

Every run stores its full per-step ``step_stream`` (P0-B/P0-C): pre,
predictor-stage and post CFL numbers, viscous stability number, energy-balance
defect against ``dE/dt = -nu \int |\omega|^2 dV - nu (\text{wall work})``,
circulation, parity, divergence and Poisson algebraic residual maxima over
**all** accepted steps.

This experiment is a numerical observation.  It is not a reproduction of any
published value, not evidence of singularity formation and not a proof.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

from ns_certificate_lab._integrity import strict_json_loads, write_with_digest
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.nonlinear_cylinder import (
    DIAGNOSTIC_FIELDS,
    GATE_SUMMARY_FIELDS,
    STEP_STREAM_FIELDS,
    VALID_INTEGRATORS,
    cartesian_vorticity,
    constrain_state,
    hou_initial_swirl,
    integrate,
    save_checkpoint,
)
from ns_certificate_lab.provenance import collect_runtime_provenance

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

# E-29b derived continuum reference (closed form; see docs/equation_audit.md).
# The decimal 7569.6226982 recorded there is this expression.
E29B_MAX_CARTESIAN_VORTICITY_EXACT = (
    24000.0 * math.pi / math.sqrt(37.0) * (36.0 / 37.0) ** 18
)

LIMITATIONS: tuple[str, ...] = (
    "single 65x128 uniform grid; spatial error dominates the temporal error "
    "by several orders at this operating point, so cross-integrator "
    "agreement bounds the time-discretization risk only",
    "agreement between explicit schemes cannot prove nonlinear stability; "
    "the frozen-coefficient von Neumann audit is a separate necessary check",
    "amplification numbers remain unconverged numerical observations; the "
    "resolution series is not in its asymptotic range "
    "(ns_certificate_lab.extrapolation)",
)


def _sha256(path_or_bytes: Path | bytes) -> str:
    if isinstance(path_or_bytes, bytes):
        return hashlib.sha256(path_or_bytes).hexdigest()
    return hashlib.sha256(path_or_bytes.read_bytes()).hexdigest()


def _prepare_output(output_dir: Path) -> None:
    if output_dir.exists():
        if any(output_dir.iterdir()):
            raise FileExistsError(
                f"output directory {output_dir} is not empty; refusing to "
                "overwrite existing evidence"
            )
    else:
        output_dir.mkdir(parents=True)


def _inside_repository(path: Path) -> bool:
    try:
        path.relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def validate_config(config: dict[str, Any]) -> None:
    required = {
        "experiment_id",
        "description",
        "grid",
        "initial_amplitude",
        "amplitude_scale",
        "viscosity_schedule",
        "t_final",
        "requested_time_steps",
        "methods",
        "diagnostic_stride",
        "seed",
        "acceptance",
        "preregistration_note",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"config is missing keys: {sorted(missing)}")
    grid = config["grid"]
    if not isinstance(grid, dict) or set(grid) != {"nr", "nz"}:
        raise ValueError("config grid must be an object with nr and nz")
    methods = config["methods"]
    if (
        not isinstance(methods, list)
        or not methods
        or len(set(methods)) != len(methods)
    ):
        raise ValueError("methods must be a non-empty list without duplicates")
    for method in methods:
        if method not in VALID_INTEGRATORS:
            raise ValueError(f"unknown integrator in config: {method!r}")
    if "heun" not in methods or len(methods) < 2:
        raise ValueError(
            "the comparison requires heun plus at least one cross-checking "
            "integrator"
        )
    steps = config["requested_time_steps"]
    if not isinstance(steps, list) or len(steps) < 2:
        raise ValueError("requested_time_steps must list at least two levels")
    values = [float(item) for item in steps]
    if any(b >= a for a, b in zip(values, values[1:])):
        raise ValueError("requested_time_steps must be strictly decreasing")
    if any(not math.isfinite(v) or v <= 0.0 for v in values):
        raise ValueError("requested_time_steps must be positive and finite")
    acceptance = config["acceptance"]
    needed = {
        "pairwise_amplification_relative_tolerance",
        "pairwise_max_u1_relative_tolerance",
        "pairwise_field_linf_relative_tolerance",
        "pairwise_differences_shrink_with_dt",
        "argmax_within_one_cell",
        "max_poisson_algebraic_residual_relative",
        "max_wall_u1_abs",
    }
    if not isinstance(acceptance, dict) or set(acceptance) != needed:
        raise ValueError(
            f"acceptance must contain exactly the keys {sorted(needed)}"
        )
    if float(config["t_final"]) <= 0.0:
        raise ValueError("t_final must be positive")


def _run_single(
    config: dict[str, Any],
    *,
    method: str,
    requested_dt: float,
) -> dict[str, Any]:
    grid = AxisymmetricGrid.uniform(
        nr=int(config["grid"]["nr"]),
        nz=int(config["grid"]["nz"]),
        r_max=1.0,
        z_min=0.0,
        z_max=1.0,
        periodic_z=True,
    )
    u1 = hou_initial_swirl(
        grid,
        amplitude=float(config["initial_amplitude"]),
        amplitude_scale=float(config["amplitude_scale"]),
    )
    initial = constrain_state(grid, u1, np.zeros(grid.shape))
    initial_vorticity = float(
        np.max(
            np.sqrt(
                sum(
                    component * component
                    for component in cartesian_vorticity(
                        grid, initial.u1, initial.omega1
                    )
                )
            )
        )
    )
    result = integrate(
        grid,
        u1=initial.u1,
        omega1=initial.omega1,
        t_final=float(config["t_final"]),
        viscosity_schedule=config["viscosity_schedule"],
        fixed_time_step=float(requested_dt),
        diagnostic_stride=int(config["diagnostic_stride"]),
        time_integrator=method,
        reference_vorticity_max=initial_vorticity,
    )
    if not result.completed:
        raise RuntimeError(f"{method} run at dt={requested_dt} did not complete")
    final = result.history[-1]
    common_reference = E29B_MAX_CARTESIAN_VORTICITY_EXACT * float(
        config["amplitude_scale"]
    )
    return {
        "grid": grid,
        "result": result,
        "metrics": {
            "method": method,
            "requested_dt": float(requested_dt),
            "effective_dt": float(result.step_stream["dt"][0]),
            "steps": result.step_count,
            "initial_max_cartesian_vorticity_discrete": initial_vorticity,
            "initial_max_cartesian_vorticity_common_reference": (
                common_reference
            ),
            "final_max_cartesian_vorticity": final["max_cartesian_vorticity"],
            "amplification_grid_normalized": final["amplification"],
            "amplification_common_reference": (
                final["max_cartesian_vorticity"] / common_reference
            ),
            "final_max_abs_u1": final["max_abs_u1"],
            "argmax_u1_r": final["argmax_u1_r"],
            "argmax_u1_z": final["argmax_u1_z"],
            "argmax_vorticity_r": final["argmax_cartesian_vorticity_r"],
            "argmax_vorticity_z": final["argmax_cartesian_vorticity_z"],
            "final_energy": final["energy"],
            "gate_summary": dict(result.gate_summary),
        },
    }


def _pairwise(
    runs: dict[str, dict[str, Any]],
    methods: list[str],
) -> list[dict[str, Any]]:
    comparisons = []
    for index, first in enumerate(methods):
        for second in methods[index + 1 :]:
            a, b = runs[first], runs[second]
            metric_a, metric_b = a["metrics"], b["metrics"]
            grid = a["grid"]
            u1_a = a["result"].state.u1
            u1_b = b["result"].state.u1
            scale = max(
                float(np.max(np.abs(u1_a))), float(np.max(np.abs(u1_b)))
            )
            amp_a = metric_a["amplification_common_reference"]
            amp_b = metric_b["amplification_common_reference"]
            comparisons.append(
                {
                    "pair": f"{first}_vs_{second}",
                    "amplification_relative_difference": abs(amp_a - amp_b)
                    / max(abs(amp_b), 1.0e-300),
                    "max_u1_relative_difference": abs(
                        metric_a["final_max_abs_u1"]
                        - metric_b["final_max_abs_u1"]
                    )
                    / max(abs(metric_b["final_max_abs_u1"]), 1.0e-300),
                    "field_linf_relative_difference": float(
                        np.max(np.abs(u1_a - u1_b))
                    )
                    / max(scale, 1.0e-300),
                    "argmax_u1_cell_offset_r": abs(
                        metric_a["argmax_u1_r"] - metric_b["argmax_u1_r"]
                    )
                    / grid.dr,
                    "argmax_u1_cell_offset_z": min(
                        abs(metric_a["argmax_u1_z"] - metric_b["argmax_u1_z"]),
                        1.0
                        - abs(
                            metric_a["argmax_u1_z"] - metric_b["argmax_u1_z"]
                        ),
                    )
                    / grid.dz,
                }
            )
    return comparisons


def evaluate(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, bool]]:
    """Run the full matrix and evaluate the preregistered acceptance checks."""

    methods = [str(m) for m in config["methods"]]
    dt_levels = [float(v) for v in config["requested_time_steps"]]
    acceptance = config["acceptance"]

    levels: list[dict[str, Any]] = []
    all_runs: dict[str, dict[str, dict[str, Any]]] = {}
    for requested_dt in dt_levels:
        runs = {
            method: _run_single(
                config, method=method, requested_dt=requested_dt
            )
            for method in methods
        }
        all_runs[f"{requested_dt:.3e}"] = runs
        levels.append(
            {
                "requested_dt": requested_dt,
                "runs": {m: runs[m]["metrics"] for m in methods},
                "pairwise": _pairwise(runs, methods),
            }
        )

    amp_tol = float(acceptance["pairwise_amplification_relative_tolerance"])
    u1_tol = float(acceptance["pairwise_max_u1_relative_tolerance"])
    linf_tol = float(acceptance["pairwise_field_linf_relative_tolerance"])
    checks: dict[str, bool] = {
        "all_runs_completed": True,
        "pairwise_amplification_within_tolerance": all(
            pair["amplification_relative_difference"] <= amp_tol
            for level in levels
            for pair in level["pairwise"]
        ),
        "pairwise_max_u1_within_tolerance": all(
            pair["max_u1_relative_difference"] <= u1_tol
            for level in levels
            for pair in level["pairwise"]
        ),
        "pairwise_field_linf_within_tolerance": all(
            pair["field_linf_relative_difference"] <= linf_tol
            for level in levels
            for pair in level["pairwise"]
        ),
        "argmax_within_one_cell": all(
            pair["argmax_u1_cell_offset_r"] <= 1.0
            and pair["argmax_u1_cell_offset_z"] <= 1.0
            for level in levels
            for pair in level["pairwise"]
        )
        if bool(acceptance["argmax_within_one_cell"])
        else True,
        "poisson_algebraic_residual_bounded": all(
            level["runs"][m]["gate_summary"][
                "max_poisson_algebraic_residual_relative"
            ]
            <= float(acceptance["max_poisson_algebraic_residual_relative"])
            for level in levels
            for m in methods
        ),
        "wall_u1_bounded": all(
            level["runs"][m]["gate_summary"]["max_wall_u1_abs"]
            <= float(acceptance["max_wall_u1_abs"])
            for level in levels
            for m in methods
        ),
    }
    if bool(acceptance["pairwise_differences_shrink_with_dt"]):
        shrink = True
        for coarse, fine in zip(levels, levels[1:]):
            for pair_coarse, pair_fine in zip(
                coarse["pairwise"], fine["pairwise"]
            ):
                if (
                    pair_fine["amplification_relative_difference"]
                    > pair_coarse["amplification_relative_difference"]
                ):
                    shrink = False
        checks["pairwise_differences_shrink_with_dt"] = shrink

    payload = {
        "levels": levels,
        "runs_by_level": all_runs,
    }
    return payload, checks


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Run the comparison and write a checksummed evidence bundle."""

    validate_config(config)
    provenance = collect_runtime_provenance()
    config_bytes = (
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    _prepare_output(output_dir)
    (output_dir / "config.snapshot.json").write_bytes(config_bytes)
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir()

    payload, checks = evaluate(config)

    diagnostic_rows = []
    stream_arrays: dict[str, np.ndarray] = {}
    for dt_label, runs in payload["runs_by_level"].items():
        for method, entry in runs.items():
            result = entry["result"]
            grid = entry["grid"]
            tag = f"{method}_dt{dt_label}"
            for record in result.history:
                diagnostic_rows.append(
                    {"method": method, "requested_dt": dt_label, **record}
                )
            for name in STEP_STREAM_FIELDS:
                stream_arrays[f"{name}__{tag}"] = result.step_stream[name]
            save_checkpoint(
                checkpoint_dir / f"final_{tag}.npz",
                grid=grid,
                state=result.state,
                time=result.time,
                viscosity=float(config["viscosity_schedule"][-1][1]),
                seed=int(config["seed"]),
                config={
                    "experiment_id": config["experiment_id"],
                    "method": method,
                    "requested_dt": dt_label,
                },
                provenance=provenance,
                metadata={"integrator": method},
            )
    _write_csv(
        output_dir / "diagnostics.csv",
        ["method", "requested_dt", *DIAGNOSTIC_FIELDS],
        diagnostic_rows,
    )
    np.savez_compressed(output_dir / "step_streams.npz", **stream_arrays)

    summary: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "interpretation": (
            "Gate 1 cross-integrator comparison. Numerical observation only: "
            "agreement between explicit schemes bounds time-discretization "
            "risk at this operating point; it is not a stability proof, not "
            "a convergence claim and not evidence of singularity formation."
        ),
        "decision_rule": (
            "Amplification obtained with Heun alone must not be used for "
            "candidate decisions (P0-A item 7).  Any future candidate "
            "decision must cite this comparison or a successor."
        ),
        "common_reference_normalization": {
            "value": E29B_MAX_CARTESIAN_VORTICITY_EXACT,
            "closed_form": "24000*pi*37^(-1/2)*(36/37)^18",
            "source": "E-29b, docs/equation_audit.md",
        },
        "levels": payload["levels"],
        "gate_summary_fields": list(GATE_SUMMARY_FIELDS),
        "acceptance_checks": checks,
        "accepted_as_gate1_comparison": bool(all(checks.values())),
        "reproducibility": {
            "seed": int(config["seed"]),
            "config_sha256": _sha256(config_bytes),
            "runtime_provenance": provenance,
            "command": (
                "python experiments/run_integrator_comparison.py --config "
                "configs/integrator_comparison.json --output-dir "
                "outputs/integrator_comparison"
            ),
        },
        "limitations": list(LIMITATIONS),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    artifact_names = [
        "config.snapshot.json",
        "diagnostics.csv",
        "step_streams.npz",
        "summary.json",
    ]
    artifact_names.extend(
        sorted(
            path.relative_to(output_dir).as_posix()
            for path in checkpoint_dir.rglob("*")
            if path.is_file()
        )
    )
    manifest = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "files": {
            name: {
                "sha256": _sha256(output_dir / name),
                "bytes": (output_dir / name).stat().st_size,
            }
            for name in sorted(artifact_names)
        },
    }
    write_with_digest(
        output_dir / "manifest.json",
        (
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8"),
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "integrator_comparison.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "integrator_comparison",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="integrator comparison config",
    )
    if not isinstance(config, dict):
        parser.error("config must be a JSON object")
    summary = run(config, output_dir)
    print(json.dumps(summary["acceptance_checks"], indent=2, sort_keys=True))
    return 0 if summary["accepted_as_gate1_comparison"] else 2


if __name__ == "__main__":
    sys.exit(main())
