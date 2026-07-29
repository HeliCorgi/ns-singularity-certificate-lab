r"""Gate 7 -- leaving the Picard region and following the families to ``tau = O(1)``.

Gate 6 ended with thirty-two runs that were all inside the first Picard iterate,
diagnosed after the fact by the exact ``A^2`` scaling of ``max|omega_1|``.  This
experiment fixes the three things that caused it.

1. **The amplitude scale was wrong.**  The sweep used the *nominal* parameter of
   the datum, but the field it builds peaks at ``e^{-1} g(c)`` times that
   parameter -- about ``1/23`` for the shapes used.  So "the sweep reached
   ``tau = A t``" was measuring a ``tau`` twenty-three times larger than the one
   the equations see.  Section 1 recomputes every legacy point in the true
   dimensionless coordinates and reports the ``tau`` actually reached.
2. **The distance to the Picard iterates was inferred, not measured.**  Here the
   iterates are integrated *alongside* the full solution by
   :class:`~ns_certificate_lab.picard_continuation.PicardLadder`, with the same
   integrator and the same steps, so the distance is a directly measured field
   difference.
3. **The sweep explored ``A``, ``L`` and ``nu`` separately.**  Section 2 sweeps
   the single parameter ``Re = A L^2 / nu`` instead.

Everything preregistered -- the ``tau`` checkpoints, the ``Re`` series, the
departure thresholds and the two promotion criteria -- is in the config file and
is copied verbatim into the output before any run starts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.initial_data import (
    CompositeSwirlFamily,
    family_A,
    family_H,
    family_S,
)
from ns_certificate_lab.nondimensional import (
    ScalingExponents,
    Setting,
    deduplicate_settings,
    nonlinear_time_estimate,
    reynolds_number,
)
from ns_certificate_lab.picard_continuation import (
    ContinuationCheckpoint,
    ContinuationResult,
    PicardLadder,
)

EXPECTED_SCHEMA_VERSION = 1
EXPECTED_EXPERIMENT = "tau_continuation_gate7"

FAMILY_BUILDERS: dict[str, Callable[[float], CompositeSwirlFamily]] = {
    "S": family_S,
    "A": family_A,
    "H": family_H,
}


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
    required = {
        "schema_version", "experiment", "preregistration", "amendments",
        "reference_reynolds", "grids", "domain", "nominal_amplitude",
        "boundary_mode", "integrator", "verification_integrator", "legacy_sweep",
    }
    if set(config) != required:
        raise ValueError("Gate 7 config has missing or unknown keys")
    if config["schema_version"] != EXPECTED_SCHEMA_VERSION:
        raise ValueError("unsupported Gate 7 config schema")
    if config["experiment"] != EXPECTED_EXPERIMENT:
        raise ValueError("experiment must equal the audited canonical value")
    if not config["preregistration"].get("recorded_before_any_run"):
        raise ValueError("the preregistration block must assert its own timing")


# --------------------------------------------------------------------------- #
# Section 1 -- the legacy sweep in dimensionless coordinates                    #
# --------------------------------------------------------------------------- #


def _effective_amplitude(
    builder: Callable[..., CompositeSwirlFamily], grid: AxisymmetricGrid, nominal: float
) -> float:
    """``A = max|u_1(0)|``, the scale the equations actually see."""
    family = builder(amplitude=nominal)
    u1, _ = family.initial_state(grid)
    return float(np.max(np.abs(u1)))


def reclassify_legacy_sweep(
    config: dict[str, Any], repository_root: Path
) -> dict[str, Any]:
    """Put the thirty-two legacy points into ``(Re, aspect, c, tau)`` coordinates."""
    legacy = config["legacy_sweep"]
    path = repository_root / legacy["path"]
    if not path.exists():
        return {"available": False, "path": legacy["path"]}
    viscosity = float(legacy["viscosity"])
    physical_time = float(legacy["physical_time"])
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    # The peak of one odd component is  A * chi(0) * max_zeta chi(zeta^2) *
    # zeta / (1 + c zeta^2), independent of the supports, so the shape factor
    # can be tabulated once per concentration.
    zeta = np.linspace(-0.999, 0.999, 20001)
    chi = np.where(zeta**2 < 1.0, np.exp(-1.0 / np.maximum(1.0 - zeta**2, 1e-300)), 0.0)

    def shape_factor(concentration: float) -> float:
        profile = chi * zeta / (1.0 + concentration * zeta**2)
        return math.exp(-1.0) * float(np.max(np.abs(profile)))

    settings: list[Setting] = []
    detail: list[dict[str, Any]] = []
    for row in rows:
        nominal = float(row["amplitude"])
        radial = float(row["radial_support"])
        axial = float(row["axial_support"])
        concentration = float(row["axial_concentration"])
        amplitude = nominal * shape_factor(concentration)
        setting = Setting(
            label=row["label"],
            amplitude=amplitude,
            length=radial,
            viscosity=viscosity,
            aspect_ratio=radial / axial,
            concentration=concentration,
            physical_time=physical_time,
        )
        settings.append(setting)
        detail.append(
            {
                "label": setting.label,
                "nominal_amplitude": nominal,
                "effective_amplitude": amplitude,
                "reynolds": setting.reynolds,
                "aspect_ratio": setting.aspect_ratio,
                "concentration": concentration,
                "dimensionless_time": setting.dimensionless_time,
            }
        )
    unique, duplicates = deduplicate_settings(settings, digits=6)
    # Shape duplicates ignore tau and Re: two settings that differ only by an
    # overall length rescaling are the same *shape*, which is what a shape sweep
    # is supposed to vary.
    shape_keys = {s.shape_key() for s in settings}
    taus = [s.dimensionless_time for s in settings]
    return {
        "available": True,
        "path": legacy["path"],
        "points": len(settings),
        "distinct_computations": len(unique),
        "duplicate_pairs": [list(pair) for pair in duplicates],
        "distinct_shapes": len(shape_keys),
        "distinct_reynolds": len({round(s.reynolds, 6) for s in settings}),
        "max_dimensionless_time": max(taus),
        "min_dimensionless_time": min(taus),
        "reynolds_range": [
            min(s.reynolds for s in settings),
            max(s.reynolds for s in settings),
        ],
        "detail": detail,
    }


# --------------------------------------------------------------------------- #
# Section 2 -- the Reynolds continuation                                       #
# --------------------------------------------------------------------------- #


def run_continuation(
    *,
    family_name: str,
    reynolds: float,
    config: dict[str, Any],
    grid_name: str,
    integrator: str | None = None,
    tau_checkpoints: Sequence[float] | None = None,
    max_time_step: float | None = None,
) -> ContinuationResult:
    grids = config["grids"][grid_name]
    domain = config["domain"]
    grid = AxisymmetricGrid.uniform(
        nr=int(grids["nr"]),
        nz=int(grids["nz"]),
        r_max=float(domain["r_max"]),
        z_min=-float(domain["z_max"]),
        z_max=float(domain["z_max"]),
        periodic_z=False,
    )
    builder = FAMILY_BUILDERS[family_name]
    family = builder(amplitude=float(config["nominal_amplitude"]))
    u1, _ = family.initial_state(grid)
    amplitude = float(np.max(np.abs(u1)))
    length = family.reference_length
    viscosity = amplitude * length * length / reynolds
    ladder = PicardLadder(
        grid, viscosity=viscosity, boundary_mode=config["boundary_mode"]
    )
    checkpoints, reason, _state = ladder.run(
        u1,
        amplitude=amplitude,
        tau_checkpoints=tuple(
            tau_checkpoints
            if tau_checkpoints is not None
            else config["preregistration"]["tau_checkpoints"]
        ),
        integrator=integrator or config["integrator"],
        max_time_step=max_time_step,
    )
    return ContinuationResult(
        label=f"{family_name}_Re{reynolds:g}_{grid_name}",
        family=family_name,
        reynolds=reynolds,
        amplitude=amplitude,
        reference_length=length,
        viscosity=viscosity,
        nr=int(grids["nr"]),
        nz=int(grids["nz"]),
        r_max=float(domain["r_max"]),
        z_max=float(domain["z_max"]),
        integrator=integrator or config["integrator"],
        checkpoints=tuple(checkpoints),
        reached_tau=checkpoints[-1].tau if checkpoints else 0.0,
        stopped_reason=reason,
    )


# --------------------------------------------------------------------------- #
# Section 3 -- the departure gate and the promotion decision                    #
# --------------------------------------------------------------------------- #


def _relative(a: float, b: float) -> float:
    scale = max(abs(a), abs(b), 1.0e-300)
    return abs(a - b) / scale


def evaluate_departure(
    result: ContinuationResult,
    *,
    coarse: ContinuationResult | None,
    integrator_check: ContinuationResult | None,
    time_step_check: ContinuationResult | None,
    thresholds: dict[str, float],
) -> dict[str, Any]:
    """The preregistered first-Picard-departure gate, item by item."""
    if not result.checkpoints:
        return {"passed": False, "reason": "no checkpoint completed"}
    final = result.final

    def matched(other: ContinuationResult | None) -> ContinuationCheckpoint | None:
        if other is None:
            return None
        for checkpoint in other.checkpoints:
            if math.isclose(checkpoint.tau, final.tau, rel_tol=1e-12):
                return checkpoint
        return None

    coarse_final = matched(coarse)
    integrator_final = matched(integrator_check)
    time_step_final = matched(time_step_check)
    checks = {
        "distance_to_first_picard": {
            "value": final.distance_to_first_picard,
            "threshold": thresholds["relative_distance_to_first_picard"],
            "passed": final.distance_to_first_picard
            >= thresholds["relative_distance_to_first_picard"],
        },
        "advection_measurable": {
            "value": final.advection_over_source,
            "threshold": thresholds["measurable_advection_over_source"],
            "passed": final.advection_over_source
            >= thresholds["measurable_advection_over_source"],
        },
        "front_points_radial": {
            "value": final.front_points_radial,
            "threshold": thresholds["front_points_radial"],
            "passed": final.front_points_radial >= thresholds["front_points_radial"],
        },
        "front_points_axial": {
            "value": final.front_points_axial,
            "threshold": thresholds["front_points_axial"],
            "passed": final.front_points_axial >= thresholds["front_points_axial"],
        },
        "domain_tail_below_signal": {
            "value": final.outer_band_fraction,
            "threshold": thresholds["outer_band_fraction"],
            "passed": final.outer_band_fraction <= thresholds["outer_band_fraction"],
        },
        "divergence_free": {
            "value": final.relative_divergence,
            "threshold": thresholds["relative_divergence"],
            "passed": final.relative_divergence <= thresholds["relative_divergence"],
        },
    }
    if coarse_final is not None:
        value = _relative(
            final.distance_to_first_picard, coarse_final.distance_to_first_picard
        )
        checks["resolution_stability"] = {
            "value": value,
            "threshold": thresholds["resolution_stability"],
            "passed": value <= thresholds["resolution_stability"],
        }
    if integrator_final is not None:
        value = _relative(final.max_omega1, integrator_final.max_omega1)
        checks["integrator_agreement"] = {
            "value": value,
            "threshold": thresholds["integrator_agreement"],
            "passed": value <= thresholds["integrator_agreement"],
        }
    if time_step_final is not None:
        value = _relative(final.max_omega1, time_step_final.max_omega1)
        checks["time_step_agreement"] = {
            "value": value,
            "threshold": thresholds["time_step_agreement"],
            "passed": value <= thresholds["time_step_agreement"],
        }
    return {
        "label": result.label,
        "tau": final.tau,
        "checks": checks,
        "passed": all(entry["passed"] for entry in checks.values()),
        "failed": sorted(k for k, v in checks.items() if not v["passed"]),
    }


def evaluate_promotion(
    result: ContinuationResult, promotion: dict[str, Any]
) -> dict[str, Any]:
    """Critical-concentration or multiscale, else rejection with a reason."""
    if len(result.checkpoints) < 2:
        return {"promoted": False, "reason": "fewer than two checkpoints"}
    first, final = result.checkpoints[0], result.final
    concentration = promotion["critical_concentration"]
    multiscale = promotion["multiscale"]
    l3_growth = final.l3_norm / max(first.l3_norm, 1.0e-300)
    shell_change = final.effective_shell_count - first.effective_shell_count
    radial_ratio = final.radial_scale / max(first.radial_scale, 1.0e-300)
    axial_ratio = final.axial_scale / max(first.axial_scale, 1.0e-300)
    outer_growth = final.shell_flux_upward / max(first.shell_flux_upward, 1.0e-300)
    critical = (
        l3_growth >= concentration["l3_growth"]
        and shell_change <= concentration["shell_count_change"]
        and min(radial_ratio, axial_ratio) <= concentration["width_contraction"]
    )
    multi = (
        shell_change >= multiscale["shell_count_growth"]
        and outer_growth >= multiscale["outer_shell_fraction_growth"]
    )
    reasons = []
    if not critical:
        reasons.append(
            f"critical-concentration: L3 growth {l3_growth:.4f} "
            f"(needs >= {concentration['l3_growth']}), shell change {shell_change:+.4f} "
            f"(needs <= {concentration['shell_count_change']}), narrowest width ratio "
            f"{min(radial_ratio, axial_ratio):.4f} "
            f"(needs <= {concentration['width_contraction']})"
        )
    if not multi:
        reasons.append(
            f"multiscale: shell change {shell_change:+.4f} "
            f"(needs >= {multiscale['shell_count_growth']}), outer-shell growth "
            f"{outer_growth:.4f} (needs >= {multiscale['outer_shell_fraction_growth']})"
        )
    return {
        "promoted": bool(critical or multi),
        "critical_concentration": critical,
        "multiscale": multi,
        "l3_growth": l3_growth,
        "shell_count_change": shell_change,
        "radial_scale_ratio": radial_ratio,
        "axial_scale_ratio": axial_ratio,
        "outer_shell_growth": outer_growth,
        "rejection_reasons": reasons,
    }


# --------------------------------------------------------------------------- #
# Section 4 -- one space-time slab certificate                                  #
# --------------------------------------------------------------------------- #


def build_flagship_slab(
    *, family_name: str, reynolds: float, config: dict[str, Any], tau: float
) -> dict[str, Any]:
    """Advance to ``tau`` on the coarse grid, then certify one accepted step.

    The right-hand sides handed to the certificate are the ladder's own, taken
    at the start, the midpoint and the end of the same accepted step.  They are
    never reconstructed from stored snapshots.
    """
    grids = config["grids"]["coarse"]
    domain = config["domain"]
    grid = AxisymmetricGrid.uniform(
        nr=int(grids["nr"]), nz=int(grids["nz"]),
        r_max=float(domain["r_max"]), z_min=-float(domain["z_max"]),
        z_max=float(domain["z_max"]), periodic_z=False,
    )
    family = FAMILY_BUILDERS[family_name](amplitude=float(config["nominal_amplitude"]))
    u1, _ = family.initial_state(grid)
    amplitude = float(np.max(np.abs(u1)))
    viscosity = amplitude * family.reference_length**2 / reynolds
    ladder = PicardLadder(grid, viscosity=viscosity, boundary_mode=config["boundary_mode"])
    _checkpoints, reason, state = ladder.run(
        u1, amplitude=amplitude, tau_checkpoints=(tau,),
        integrator=config["integrator"],
    )
    if reason != "completed":
        raise RuntimeError(f"slab preparation stopped early: {reason}")

    dt = ladder.stable_time_step(state)
    rhs_start = ladder.rhs(state)
    half = ladder.step(state, 0.5 * dt, config["integrator"])
    rhs_mid = ladder.rhs(half)
    end = ladder.step(state, dt, config["integrator"])
    rhs_end = ladder.rhs(end)

    from ns_certificate_lab.free_space_recovery import recover_free_space_velocity
    from ns_certificate_lab.slab_certificate import (
        build_slab_certificate,
        verify_slab_certificate,
    )

    recovery = recover_free_space_velocity(
        ladder.solver, state[7], support_threshold=1.0e-10
    )
    certificate = build_slab_certificate(
        grid,
        time_start=tau / amplitude,
        time_step=dt,
        u1_start=state[6], u1_end=end[6],
        omega1_start=state[7], omega1_end=end[7],
        psi1_start=ladder.solver.solve(state[7]),
        psi1_end=ladder.solver.solve(end[7]),
        du1_start=rhs_start[6], du1_end=rhs_end[6],
        domega1_start=rhs_start[7], domega1_end=rhs_end[7],
        domega1_mid=rhs_mid[7],
        viscosity=viscosity,
        green_tail_bound=(
            recovery.interior_tail.gradient if recovery.tail_bound_available else 0.0
        ),
        interior_radius=1.5,
    )
    payload = json.loads(json.dumps(certificate.as_dict(), allow_nan=False))
    verdict = verify_slab_certificate(payload)

    # Tamper test: narrowing an enclosure must be rejected.
    tampered = json.loads(json.dumps(payload))
    lower, upper = tampered["enclosures"]["omega1"]
    tampered["enclosures"]["omega1"] = [upper, lower]
    narrowed = verify_slab_certificate(tampered)
    # Tamper test: claiming a hypothesis is proved must be rejected.
    forged = json.loads(json.dumps(payload))
    forged["hypotheses"]["H2_hermite_remainder"]["proved"] = True
    forged_verdict = verify_slab_certificate(forged)

    return {
        "payload": payload,
        "verified": verdict,
        "tamper_rejected_inverted_enclosure": not narrowed["verified"],
        "tamper_rejected_forged_proof": not forged_verdict["verified"],
        "family": family_name,
        "reynolds": reynolds,
        "tau": tau,
    }


# --------------------------------------------------------------------------- #
# driver                                                                       #
# --------------------------------------------------------------------------- #


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    _validate(config)
    prereg = config["preregistration"]
    output_dir.mkdir(parents=True, exist_ok=True)
    repository_root = config_path.resolve().parent.parent

    written: list[Path] = []

    def write_json(name: str, payload: object) -> None:
        path = output_dir / name
        _write_json(path, payload)
        written.append(path)

    def write_csv(name: str, rows: list[dict[str, Any]], columns: list[str]) -> None:
        path = output_dir / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in columns})
        written.append(path)

    # The preregistration is copied out before anything runs.
    write_json("preregistration.json", prereg)
    write_json("config.snapshot.json", config)

    legacy = reclassify_legacy_sweep(config, repository_root)

    # -- the Reynolds series ------------------------------------------------- #
    results: list[ContinuationResult] = []
    for family_name in prereg["families"]:
        for reynolds in prereg["reynolds_series"]:
            results.append(
                run_continuation(
                    family_name=family_name,
                    reynolds=float(reynolds),
                    config=config,
                    grid_name="medium",
                )
            )

    # -- the flagship and its verification ----------------------------------- #
    complete = [r for r in results if r.checkpoints and r.reached_tau >= 0.4]
    if not complete:
        raise RuntimeError("no run reached the preregistered tau = 0.4 checkpoint")
    flagship = max(complete, key=lambda r: r.final.distance_to_first_picard)
    verification_taus = tuple(
        tau for tau in prereg["tau_checkpoints"] if tau <= flagship.reached_tau
    )
    coarse = run_continuation(
        family_name=flagship.family, reynolds=flagship.reynolds, config=config,
        grid_name="coarse", tau_checkpoints=verification_taus,
    )
    fine = run_continuation(
        family_name=flagship.family, reynolds=flagship.reynolds, config=config,
        grid_name="fine", tau_checkpoints=verification_taus,
    )
    ultra = run_continuation(
        family_name=flagship.family, reynolds=flagship.reynolds, config=config,
        grid_name="ultra", tau_checkpoints=verification_taus,
    )
    # Both verification runs sit on the same grid as the gate, so what they
    # measure is integrator error and time-step error, not grid error.
    integrator_check = run_continuation(
        family_name=flagship.family, reynolds=flagship.reynolds, config=config,
        grid_name="ultra", integrator=config["verification_integrator"],
        tau_checkpoints=verification_taus,
    )
    halved = 0.5 * ultra.final.accepted_time_step
    time_step_check = run_continuation(
        family_name=flagship.family, reynolds=flagship.reynolds, config=config,
        grid_name="ultra", tau_checkpoints=verification_taus, max_time_step=halved,
    )
    results.extend([coarse, fine, ultra, integrator_check, time_step_check])

    # The gate is evaluated on the finest grid; resolution stability compares it
    # against the next grid down, so the comparison stays a genuine two-grid one.
    gate = evaluate_departure(
        ultra,
        coarse=fine,
        integrator_check=integrator_check,
        time_step_check=time_step_check,
        thresholds=prereg["departure_gate"],
    )
    promotions = {
        result.label: evaluate_promotion(result, prereg["promotion"])
        for result in results
    }
    promoted = sorted(label for label, v in promotions.items() if v["promoted"])

    # -- how the departure grows with tau ------------------------------------ #
    departure_law = []
    for checkpoint in flagship.checkpoints:
        if checkpoint.tau > 0.0 and checkpoint.advection_over_source > 0.0:
            departure_law.append(
                {
                    "tau": checkpoint.tau,
                    "advection_over_source": checkpoint.advection_over_source,
                    "quadratic_coefficient": checkpoint.advection_over_source
                    / checkpoint.tau**2,
                    "distance_to_first_picard": checkpoint.distance_to_first_picard,
                    "tau_for_ten_percent": nonlinear_time_estimate(
                        measured_ratio=checkpoint.advection_over_source,
                        measured_tau=checkpoint.tau,
                    ),
                }
            )

    slab = build_flagship_slab(
        family_name=flagship.family, reynolds=flagship.reynolds, config=config,
        tau=0.4,
    )
    write_json("slab_certificate.json", slab["payload"])

    exponents = ScalingExponents()
    summary = {
        "schema_version": EXPECTED_SCHEMA_VERSION,
        "experiment": EXPECTED_EXPERIMENT,
        "preregistration": prereg,
        "amendments": config["amendments"],
        "scaling_exponents": exponents.as_dict(),
        "legacy_sweep_reclassification": legacy,
        "flagship": flagship.label,
        "flagship_reached_tau": flagship.reached_tau,
        "departure_gate": gate,
        "departure_law": departure_law,
        "promotions": promotions,
        "promoted": promoted,
        "slab_certificate": {
            key: value for key, value in slab.items() if key != "payload"
        },
        "runs": [result.as_dict() for result in results],
        "limitations": [
            "This is a resolved numerical continuation, not a proof.  Nothing "
            "here establishes or refutes a Clay singularity.",
            "The departure gate measures distance to the Picard iterates of the "
            "*discrete* system on this grid; it certifies that the computation "
            "has left the linear-response regime, not that the continuous "
            "solution has.",
            "Viscosity is fixed and positive in every run; no time-dependent or "
            "two-stage viscosity protocol is used anywhere.",
        ],
    }
    write_json("summary.json", summary)

    rows = []
    for result in results:
        for checkpoint in result.checkpoints:
            row = {"label": result.label, "family": result.family,
                   "reynolds": result.reynolds, "nr": result.nr, "nz": result.nz,
                   "integrator": result.integrator}
            row.update(checkpoint.as_dict())
            row.pop("shell_fractions", None)
            rows.append(row)
    write_csv(
        "checkpoints.csv", rows,
        ["label", "family", "reynolds", "nr", "nz", "integrator", "tau",
         "physical_time", "steps", "accepted_time_step", "max_u1", "max_omega1",
         "max_speed", "distance_to_first_picard", "distance_to_second_picard",
         "advection_over_source", "viscous_over_source", "rhs_change_rate",
         "kinetic_energy", "enstrophy", "l3_norm", "core_l3_fraction",
         "effective_shell_count", "radial_scale", "axial_scale",
         "shell_flux_upward", "outer_band_fraction", "relative_divergence",
         "curl_defect", "poisson_residual", "front_points_radial",
         "front_points_axial"],
    )
    if legacy.get("available"):
        write_csv(
            "legacy_reclassified.csv", legacy["detail"],
            ["label", "nominal_amplitude", "effective_amplitude", "reynolds",
             "aspect_ratio", "concentration", "dimensionless_time"],
        )

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
                "flagship": summary["flagship"],
                "reached_tau": summary["flagship_reached_tau"],
                "departure_gate_passed": summary["departure_gate"]["passed"],
                "departure_gate_failed": summary["departure_gate"]["failed"],
                "promoted": summary["promoted"],
                "slab_verified": summary["slab_certificate"]["verified"]["verified"],
                "legacy_max_tau": summary["legacy_sweep_reclassification"].get(
                    "max_dimensionless_time"
                ),
            },
            ensure_ascii=False, sort_keys=True, indent=2,
        )
    )


if __name__ == "__main__":
    main()
