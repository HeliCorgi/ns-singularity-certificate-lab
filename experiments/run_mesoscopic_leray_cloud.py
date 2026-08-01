"""Run the mesoscopic true-Leray cloud falsification experiment.

The output deliberately separates four statements with different epistemic
status: exact finite carrier algebra, deterministic zero-padded local-FFT
measurements, a small binary64 Fourier--Galerkin time integration, and an
analytic critical-energy Duhamel obstruction.  None of them is a proof of a
Navier--Stokes singularity or a continuum PDE orbit.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from ns_certificate_lab.carrier_two_stage_galerkin import (
    run_carrier_two_stage_galerkin,
)
from ns_certificate_lab._integrity import (
    canonical_json_bytes,
    digest_sidecar,
    sha256_file,
    strict_json_loads,
    verify_digest,
    write_with_digest,
)
from ns_certificate_lab.exact_carrier_search import search_exact_carrier_gadget
from ns_certificate_lab.exact_carrier_record_verifier import (
    verify_serialized_strict_orientation_records,
)
from ns_certificate_lab.mesoscopic_cloud_scaling import (
    MesoscopicCloudConfig,
    MesoscopicCloudMetrics,
)
from ns_certificate_lab.mesoscopic_galerkin import run_small_mesoscopic_galerkin
from ns_certificate_lab.mesoscopic_local_fft import (
    measure_local_fft_mesoscopic_metrics,
)
from ns_certificate_lab.provenance import (
    collect_runtime_provenance,
    source_fingerprint,
    validate_runtime_provenance,
)


CONFIG_SCHEMA = "ns-certificate-lab/mesoscopic-leray-cloud-config/v1"
OUTPUT_SCHEMA = "ns-certificate-lab/mesoscopic-leray-cloud/v1"
FIT_SCHEMA = "ns-certificate-lab/mesoscopic-leray-scaling-fits/v1"
GALERKIN_SCHEMA = "ns-certificate-lab/mesoscopic-leray-galerkin/v1"
MANIFEST_SCHEMA = "ns-certificate-lab/mesoscopic-leray-cloud-manifest/v1"
STATUS = "NUMERICAL FALSIFICATION SCREEN / NO SINGULARITY PROOF"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MINIMUM_FIT_POINTS = 4

PAYLOAD_NAMES = (
    "config.snapshot.json",
    "exact_carrier_certificate.json",
    "mesoscopic_scaling.csv",
    "relative_width_scaling.csv",
    "scaling_fits.json",
    "galerkin_comparison.json",
    "shell_energy.svg",
    "summary.json",
)
SOURCE_PATHS = (
    "src/ns_certificate_lab/_integrity.py",
    "src/ns_certificate_lab/provenance.py",
    "src/ns_certificate_lab/fourier_torus.py",
    "src/ns_certificate_lab/torus_chain.py",
    "src/ns_certificate_lab/leray_response_relay.py",
    "src/ns_certificate_lab/exact_carrier_search.py",
    "src/ns_certificate_lab/exact_carrier_record_verifier.py",
    "src/ns_certificate_lab/carrier_two_stage_galerkin.py",
    "src/ns_certificate_lab/mesoscopic_cloud_scaling.py",
    "src/ns_certificate_lab/mesoscopic_galerkin.py",
    "src/ns_certificate_lab/mesoscopic_local_fft.py",
    "experiments/run_mesoscopic_leray_cloud.py",
    "tests/test_exact_carrier_search.py",
    "tests/test_exact_carrier_record_verifier.py",
    "tests/test_carrier_two_stage_galerkin.py",
    "tests/test_mesoscopic_cloud_scaling.py",
    "tests/test_mesoscopic_galerkin.py",
    "tests/test_mesoscopic_local_fft.py",
    "tests/test_mesoscopic_leray_experiment.py",
)

ROW_FIELDS = (
    "family",
    "family_parameter",
    "gamma",
    "relative_width_fraction",
    "base_scale",
    "width",
    "half_width",
    "realized_width_over_scale",
    "width_source",
    "support_mode_count",
    "effective_mode_count",
    "angle_sine",
    "measurement_available",
    "measurement_method",
    "local_fft_input_side_length",
    "local_fft_padded_side_length",
    "local_fft_zero_padding_alias_free",
    "local_fft_global_overlap_aware",
    "local_fft_estimated_peak_working_bytes",
    "local_fft_maximum_working_bytes",
    "targets_per_band",
    "parent_boxes_disjoint",
    "child_band_isolated",
    "channel_bands_disjoint",
    "full_output_bands_disjoint",
    "critical_energy",
    "parent_norm_squared",
    "parent_min_wavenumber",
    "parent_rms_wavenumber",
    "parent_max_wavenumber",
    "a_critical",
    "a_unit",
    "chi_support",
    "chi_effective",
    "gain_g",
    "heat_factor_h",
    "duhamel_amplitude_ratio",
    "duhamel_energy_ratio",
    "duhamel_identity_rhs",
    "duhamel_identity_relative_error",
    "duhamel_upper_bound_method",
    "support_radius_factor_kappa",
    "universal_duhamel_upper_bound",
    "effective_duhamel_upper_bound",
    "duhamel_to_universal_upper_bound",
    "duhamel_to_effective_upper_bound",
    "required_energy_constant_for_target",
    "forcing_metric_method",
    "child_forcing_fraction",
    "low_side_forcing_ratio",
    "off_chain_main_ratio",
    "off_core_main_ratio",
    "child_spill_main_ratio",
    "outside_child_full_main_ratio",
    "difference_sideband_main_ratio",
    "self_interaction_cross_interaction_ratio",
    "energy_cancellation_available",
    "energy_cancellation_method",
    "energy_cancellation_residual",
    "energy_cancellation_pairing",
    "full_nonlinear_norm",
    "normalization_relative_error",
    "divergence_relative",
    "reality_relative",
    "parent_divergence_relative",
    "parent_reality_relative",
    "nonlinear_divergence_relative",
    "nonlinear_reality_relative",
    "channel_ratios_are_orthogonal",
    "channel_norms_json",
    "channel_ratios_json",
    "box_overlap_counts_json",
    "predicted_a_unit_exponent",
    "predicted_a_critical_exponent",
    "predicted_g_exponent",
    "predicted_duhamel_exponent",
    "predicted_g_grows",
    "predicted_duhamel_decays",
    "finite_screen_pass",
    "asymptotic_screen_pass",
    "finite_screen_reasons",
    "asymptotic_screen_reasons",
    "core_fit_eligible",
    "relay_fit_eligible",
)

FIT_METRICS = (
    "support_mode_count",
    "a_unit",
    "a_critical",
    "gain_g",
    "duhamel_energy_ratio",
    "outside_child_full_main_ratio",
    "self_interaction_cross_interaction_ratio",
)


def _finite_number(value: object, *, name: str, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be numeric")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _positive_integer(value: object, *, name: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return value


def _unique_increasing_numbers(
    value: object,
    *,
    name: str,
    integer: bool,
    lower: float,
    upper: float | None,
) -> list[int] | list[float]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a nonempty list")
    output: list[int] | list[float]
    if integer:
        integers = [
            _positive_integer(item, name=f"{name} item", minimum=int(lower))
            for item in value
        ]
        output = integers
    else:
        numbers = [_finite_number(item, name=f"{name} item") for item in value]
        if any(item <= lower or (upper is not None and item >= upper) for item in numbers):
            raise ValueError(f"{name} items must lie strictly between {lower} and {upper}")
        output = numbers
    if any(right <= left for left, right in zip(output, output[1:])):
        raise ValueError(f"{name} must be strictly increasing without duplicates")
    return output


def _validate_config(value: object) -> dict[str, Any]:
    required = {
        "schema",
        "seed",
        "viscosity",
        "energy_constant",
        "tau",
        "scales",
        "gammas",
        "relative_width_fractions",
        "exact_target_limit",
        "strata_per_axis",
        "exact_energy_pair_limit",
        "local_fft_maximum_working_bytes",
        "required_duhamel_ratio",
        "small_galerkin",
        "carrier_two_stage",
    }
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("mesoscopic config has missing or unknown fields")
    if value["schema"] != CONFIG_SCHEMA:
        raise ValueError("mesoscopic config schema is invalid")
    _positive_integer(value["seed"], name="seed", minimum=0)
    if isinstance(value["viscosity"], bool) or not isinstance(
        value["viscosity"], (str, int)
    ):
        raise ValueError("viscosity must be an integer or rational string")
    try:
        viscosity = Fraction(value["viscosity"])
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("viscosity is not a valid rational") from error
    if viscosity <= 0:
        raise ValueError("viscosity must be positive")
    _finite_number(value["energy_constant"], name="energy_constant", positive=True)
    _finite_number(value["tau"], name="tau", positive=True)
    _unique_increasing_numbers(
        value["scales"], name="scales", integer=True, lower=2, upper=None
    )
    _unique_increasing_numbers(
        value["gammas"], name="gammas", integer=False, lower=0.0, upper=1.0
    )
    _unique_increasing_numbers(
        value["relative_width_fractions"],
        name="relative_width_fractions",
        integer=False,
        lower=0.0,
        upper=1.0 / 3.0,
    )
    if (
        isinstance(value["exact_target_limit"], bool)
        or not isinstance(value["exact_target_limit"], int)
        or value["exact_target_limit"] < 0
    ):
        raise ValueError("exact_target_limit must be a nonnegative integer")
    if (
        isinstance(value["exact_energy_pair_limit"], bool)
        or not isinstance(value["exact_energy_pair_limit"], int)
        or value["exact_energy_pair_limit"] < 0
    ):
        raise ValueError("exact_energy_pair_limit must be a nonnegative integer")
    _positive_integer(
        value["local_fft_maximum_working_bytes"],
        name="local_fft_maximum_working_bytes",
    )
    strata = _positive_integer(value["strata_per_axis"], name="strata_per_axis")
    if strata > 32:
        raise ValueError("strata_per_axis must not exceed 32")
    _finite_number(
        value["required_duhamel_ratio"],
        name="required_duhamel_ratio",
        positive=True,
    )
    small = value["small_galerkin"]
    small_fields = {"scale", "width", "grid_size", "tau", "steps"}
    if not isinstance(small, dict) or set(small) != small_fields:
        raise ValueError("small_galerkin has missing or unknown fields")
    _positive_integer(small["scale"], name="small_galerkin.scale", minimum=2)
    _positive_integer(small["width"], name="small_galerkin.width")
    _positive_integer(small["grid_size"], name="small_galerkin.grid_size", minimum=8)
    _finite_number(small["tau"], name="small_galerkin.tau", positive=True)
    _positive_integer(small["steps"], name="small_galerkin.steps")
    two_stage = value["carrier_two_stage"]
    two_stage_fields = {
        "scale",
        "grid_size",
        "tau",
        "time_multiples",
        "steps",
    }
    if not isinstance(two_stage, dict) or set(two_stage) != two_stage_fields:
        raise ValueError("carrier_two_stage has missing or unknown fields")
    _positive_integer(two_stage["scale"], name="carrier_two_stage.scale")
    _positive_integer(
        two_stage["grid_size"], name="carrier_two_stage.grid_size", minimum=8
    )
    _finite_number(two_stage["tau"], name="carrier_two_stage.tau", positive=True)
    _finite_number(
        two_stage["time_multiples"],
        name="carrier_two_stage.time_multiples",
        positive=True,
    )
    _positive_integer(two_stage["steps"], name="carrier_two_stage.steps")
    return value


def _load_config(path: Path) -> dict[str, Any]:
    value = strict_json_loads(
        Path(path).read_text(encoding="utf-8"), label="mesoscopic config"
    )
    return _validate_config(value)


def _compact_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _core_fit_eligible(metrics: MesoscopicCloudMetrics) -> bool:
    return bool(
        metrics.measurement_available
        and metrics.parent_boxes_disjoint
        and metrics.child_band_isolated
        and metrics.channel_bands_disjoint
    )


def _relay_fit_eligible(metrics: MesoscopicCloudMetrics) -> bool:
    return bool(_core_fit_eligible(metrics) and metrics.full_output_bands_disjoint)


def _metric_row(
    metrics: MesoscopicCloudMetrics,
    *,
    family: str,
    family_parameter: float,
    relative_width_fraction: float | None,
    energy_constant: float,
    required_duhamel_ratio: float,
    local_fft_maximum_working_bytes: int,
) -> dict[str, Any]:
    duhamel = metrics.duhamel_energy_ratio
    required_energy = (
        energy_constant * required_duhamel_ratio / duhamel
        if duhamel is not None and duhamel > 0.0
        else None
    )
    local_fft = metrics.measurement_method == (
        "exact-zero-padded-local-fft-global-combination"
    )
    padded_side = 4 * metrics.width - 3 if local_fft else None
    channel_norms = metrics.channel_norms
    if (
        channel_norms is not None
        and metrics.full_nonlinear_norm is not None
        and channel_norms.get("child", 0.0) > 0.0
    ):
        child_squared = channel_norms["child"] ** 2
        child_spill_squared = channel_norms["child_spill"] ** 2
        child_spill_main = math.sqrt(child_spill_squared / child_squared)
        outside_child_full_main = math.sqrt(
            max(
                metrics.full_nonlinear_norm**2
                - child_squared
                - child_spill_squared,
                0.0,
            )
            / child_squared
        )
    else:
        child_spill_main = None
        outside_child_full_main = None
    return {
        "family": family,
        "family_parameter": family_parameter,
        "gamma": metrics.gamma,
        "relative_width_fraction": relative_width_fraction,
        "base_scale": metrics.base_scale,
        "width": metrics.width,
        "half_width": metrics.half_width,
        "realized_width_over_scale": metrics.relative_width,
        "width_source": metrics.width_source,
        "support_mode_count": metrics.mode_count,
        "effective_mode_count": metrics.effective_mode_count,
        "angle_sine": metrics.angle_sine,
        "measurement_available": metrics.measurement_available,
        "measurement_method": metrics.measurement_method,
        "local_fft_input_side_length": 2 * metrics.width - 1 if local_fft else None,
        "local_fft_padded_side_length": padded_side,
        "local_fft_zero_padding_alias_free": local_fft,
        "local_fft_global_overlap_aware": local_fft,
        "local_fft_estimated_peak_working_bytes": (
            43 * padded_side**3 * np.dtype(np.complex128).itemsize
            if padded_side is not None
            else None
        ),
        "local_fft_maximum_working_bytes": local_fft_maximum_working_bytes,
        "targets_per_band": metrics.targets_per_band,
        "parent_boxes_disjoint": metrics.parent_boxes_disjoint,
        "child_band_isolated": metrics.child_band_isolated,
        "channel_bands_disjoint": metrics.channel_bands_disjoint,
        "full_output_bands_disjoint": metrics.full_output_bands_disjoint,
        "critical_energy": metrics.critical_energy,
        "parent_norm_squared": metrics.parent_norm_squared,
        "parent_min_wavenumber": metrics.parent_min_wavenumber,
        "parent_rms_wavenumber": metrics.parent_rms_wavenumber,
        "parent_max_wavenumber": metrics.parent_max_wavenumber,
        "a_critical": metrics.a_critical,
        "a_unit": metrics.a_unit,
        "chi_support": metrics.chi_support,
        "chi_effective": metrics.chi_effective,
        "gain_g": metrics.gain_g,
        "heat_factor_h": metrics.heat_factor_h,
        "duhamel_amplitude_ratio": (
            math.sqrt(duhamel) if duhamel is not None else None
        ),
        "duhamel_energy_ratio": duhamel,
        "duhamel_identity_rhs": metrics.duhamel_identity_rhs,
        "duhamel_identity_relative_error": (
            metrics.duhamel_identity_relative_error
        ),
        "duhamel_upper_bound_method": metrics.duhamel_upper_bound_method,
        "support_radius_factor_kappa": metrics.support_radius_factor_kappa,
        "universal_duhamel_upper_bound": (
            metrics.universal_duhamel_upper_bound
        ),
        "effective_duhamel_upper_bound": (
            metrics.effective_duhamel_upper_bound
        ),
        "duhamel_to_universal_upper_bound": (
            metrics.duhamel_to_universal_upper_bound
        ),
        "duhamel_to_effective_upper_bound": (
            metrics.duhamel_to_effective_upper_bound
        ),
        "required_energy_constant_for_target": required_energy,
        "forcing_metric_method": metrics.forcing_metric_method,
        "child_forcing_fraction": metrics.child_forcing_fraction,
        "low_side_forcing_ratio": metrics.low_side_forcing_ratio,
        "off_chain_main_ratio": metrics.off_chain_main_ratio,
        "off_core_main_ratio": metrics.off_chain_main_ratio,
        "child_spill_main_ratio": child_spill_main,
        "outside_child_full_main_ratio": outside_child_full_main,
        "difference_sideband_main_ratio": (
            metrics.difference_sideband_main_ratio
        ),
        "self_interaction_cross_interaction_ratio": (
            metrics.self_interaction_cross_interaction_ratio
        ),
        "energy_cancellation_available": metrics.energy_cancellation_available,
        "energy_cancellation_method": metrics.energy_cancellation_method,
        "energy_cancellation_residual": metrics.energy_cancellation_residual,
        "energy_cancellation_pairing": metrics.energy_cancellation_pairing,
        "full_nonlinear_norm": metrics.full_nonlinear_norm,
        "normalization_relative_error": metrics.normalization_relative_error,
        "divergence_relative": metrics.divergence_relative,
        "reality_relative": metrics.reality_relative,
        "parent_divergence_relative": metrics.parent_divergence_relative,
        "parent_reality_relative": metrics.parent_reality_relative,
        "nonlinear_divergence_relative": metrics.nonlinear_divergence_relative,
        "nonlinear_reality_relative": metrics.nonlinear_reality_relative,
        "channel_ratios_are_orthogonal": metrics.channel_ratios_are_orthogonal,
        "channel_norms_json": _compact_json(metrics.channel_norms),
        "channel_ratios_json": _compact_json(metrics.channel_ratios),
        "box_overlap_counts_json": _compact_json(metrics.box_overlap_counts),
        "predicted_a_unit_exponent": metrics.predicted_a_unit_exponent,
        "predicted_a_critical_exponent": metrics.predicted_a_critical_exponent,
        "predicted_g_exponent": metrics.predicted_g_exponent,
        "predicted_duhamel_exponent": metrics.predicted_duhamel_exponent,
        "predicted_g_grows": metrics.predicted_g_grows,
        "predicted_duhamel_decays": metrics.predicted_duhamel_decays,
        "finite_screen_pass": metrics.finite_screen_pass,
        "asymptotic_screen_pass": metrics.asymptotic_screen_pass,
        "finite_screen_reasons": "|".join(metrics.finite_screen_reasons),
        "asymptotic_screen_reasons": "|".join(
            metrics.asymptotic_screen_reasons
        ),
        "core_fit_eligible": _core_fit_eligible(metrics),
        "relay_fit_eligible": _relay_fit_eligible(metrics),
    }


def _csv_bytes(rows: list[dict[str, Any]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=ROW_FIELDS, extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _power_law_fit(scales: list[float], values: list[float]) -> dict[str, float]:
    if len(scales) < MINIMUM_FIT_POINTS or len(values) != len(scales):
        raise ValueError("power-law fit requires at least four paired points")
    if any(scale <= 0.0 for scale in scales) or any(value <= 0.0 for value in values):
        raise ValueError("power-law fit inputs must be positive")
    logarithmic_scales = np.log(np.asarray(scales, dtype=np.float64))
    logarithmic_values = np.log(np.asarray(values, dtype=np.float64))
    slope, intercept = np.polyfit(logarithmic_scales, logarithmic_values, 1)
    fitted = slope * logarithmic_scales + intercept
    residual = float(np.sum((logarithmic_values - fitted) ** 2))
    centered = float(
        np.sum((logarithmic_values - np.mean(logarithmic_values)) ** 2)
    )
    if centered <= np.finfo(float).tiny:
        r_squared = 1.0 if residual <= 256.0 * np.finfo(float).eps else 0.0
    else:
        r_squared = 1.0 - residual / centered
    return {
        "exponent": float(slope),
        "prefactor": float(math.exp(float(intercept))),
        "log_r_squared": float(r_squared),
    }


def _predicted_exponents(gamma: float) -> dict[str, float]:
    return {
        "support_mode_count": 3.0 * gamma,
        "a_unit": 1.0 + 1.5 * gamma,
        "a_critical": 1.5 * gamma,
        "gain_g": 1.5 * gamma - 1.0,
        "duhamel_energy_ratio": 3.0 * gamma - 3.0,
        "outside_child_full_main_ratio": gamma - 1.0,
    }


def _build_scaling_fits(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[tuple[str, float], list[dict[str, Any]]] = {}
    for row in rows:
        key = (str(row["family"]), float(row["family_parameter"]))
        groups.setdefault(key, []).append(row)
    records: list[dict[str, Any]] = []
    for (family, parameter), family_rows in groups.items():
        ordered = sorted(family_rows, key=lambda row: int(row["base_scale"]))
        gamma = float(ordered[0]["gamma"])

        def exact_full_metric(row: dict[str, Any]) -> bool:
            method = str(row["forcing_metric_method"])
            return method == "exact-sparse-full" or method.startswith(
                "exact-local-fft-full-global"
            )

        fits: dict[str, Any] = {}
        for metric in FIT_METRICS:
            if metric in {
                "support_mode_count",
                "a_unit",
                "a_critical",
                "gain_g",
                "duhamel_energy_ratio",
            }:
                eligibility_class = (
                    "diagnostic_core_geometry_not_relay_acceptance"
                )
                eligible = [
                    row
                    for row in ordered
                    if bool(row["core_fit_eligible"]) and row[metric] is not None
                ]
            elif metric == "outside_child_full_main_ratio":
                eligibility_class = (
                    "full_output_disjoint_or_exact_full_forcing_metric"
                )
                eligible = [
                    row
                    for row in ordered
                    if bool(row["core_fit_eligible"])
                    and row[metric] is not None
                    and (
                        bool(row["full_output_bands_disjoint"])
                        or exact_full_metric(row)
                    )
                ]
            else:
                eligibility_class = (
                    "full_output_disjoint_or_exact_interaction_decomposition"
                )
                eligible = [
                    row
                    for row in ordered
                    if bool(row["core_fit_eligible"])
                    and row[metric] is not None
                    and (
                        bool(row["full_output_bands_disjoint"])
                        or exact_full_metric(row)
                    )
                ]
            input_rows = (
                eligible[-MINIMUM_FIT_POINTS:]
                if len(eligible) >= MINIMUM_FIT_POINTS
                else []
            )
            base_record = {
                "eligibility_class": eligibility_class,
                "eligible_scales": [int(row["base_scale"]) for row in eligible],
                "input_scales": [int(row["base_scale"]) for row in input_rows],
            }
            if not input_rows:
                fits[metric] = {
                    "status": "insufficient_eligible_points",
                    **base_record,
                }
                continue
            values = [float(row[metric]) for row in input_rows]
            if any(not math.isfinite(value) or value <= 0.0 for value in values):
                fits[metric] = {
                    "status": "nonpositive_or_unavailable_metric",
                    **base_record,
                }
                continue
            fits[metric] = {
                "status": "fit_last_four_eligible_points",
                **base_record,
                **_power_law_fit(
                    [float(row["base_scale"]) for row in input_rows], values
                ),
            }
        records.append(
            {
                "family": family,
                "parameter_name": (
                    "gamma" if family == "power_width" else "relative_width_fraction"
                ),
                "parameter_value": parameter,
                "gamma": gamma,
                "analytical_classification": (
                    "REJECTED_SUBLINEAR_CRITICAL_DUHAMEL"
                    if gamma < 1.0
                    else "FIXED_RELATIVE_BOUNDARY_NOT_ANALYTICALLY_REJECTED"
                ),
                "core_fit_eligible_scales": [
                    int(row["base_scale"])
                    for row in ordered
                    if bool(row["core_fit_eligible"])
                ],
                "relay_fit_eligible_scales": [
                    int(row["base_scale"])
                    for row in ordered
                    if bool(row["relay_fit_eligible"])
                ],
                "window_policy": (
                    "each metric uses its last four appropriate eligible scales; "
                    "core fits are diagnostics and never relay acceptance"
                ),
                "predicted_exponents": _predicted_exponents(gamma),
                "fits": fits,
            }
        )
    return {
        "schema": FIT_SCHEMA,
        "status": STATUS,
        "minimum_fit_points": MINIMUM_FIT_POINTS,
        "groups": records,
    }


def _svg_bytes(history: list[dict[str, float]]) -> bytes:
    """Render the genuine two-stage carrier shell budget deterministically."""

    if len(history) < 2:
        raise ValueError("two-stage Galerkin history must contain at least two records")
    width, height = 940, 580
    left, right = 88.0, 904.0
    upper_top, upper_bottom = 64.0, 246.0
    lower_top, lower_bottom = 330.0, 516.0
    times = [float(row["scaled_time"]) for row in history]
    series = {
        "parent-energy": [float(row["parent_energy_ratio"]) for row in history],
        "total-energy": [float(row["total_energy_ratio"]) for row in history],
        "intended-child-energy": [
            float(row["first_child_energy_ratio"]) for row in history
        ],
        "cross-talk-energy": [
            float(row["cross_talk_energy_ratio"]) for row in history
        ],
        "intended-grandchild-energy": [
            float(row["grandchild_energy_ratio"]) for row in history
        ],
        "remainder-energy": [
            float(row["remainder_energy_ratio"]) for row in history
        ],
    }
    flat = times + [item for values in series.values() for item in values]
    if any(not math.isfinite(item) for item in flat) or max(times) <= 0.0:
        raise ValueError("two-stage Galerkin history contains invalid values")
    lower_max = max(
        max(series[name])
        for name in (
            "intended-child-energy",
            "cross-talk-energy",
            "intended-grandchild-energy",
            "remainder-energy",
        )
    )
    lower_max = max(lower_max * 1.08, np.finfo(float).tiny)
    total_min = min(series["total-energy"] + series["parent-energy"])
    total_floor = max(0.0, total_min - max(1.0 - total_min, 1.0e-6) * 0.08)

    def x_coordinate(time: float) -> float:
        return left + (right - left) * time / max(times)

    def upper_y(value: float) -> float:
        return upper_bottom - (upper_bottom - upper_top) * (
            (value - total_floor) / (1.0 - total_floor)
        )

    def lower_y(value: float) -> float:
        return lower_bottom - (lower_bottom - lower_top) * value / lower_max

    def points(name: str, transform) -> str:
        return " ".join(
            f"{x_coordinate(time):.6f},{transform(value):.6f}"
            for time, value in zip(times, series[name])
        )

    upper_lines = (
        ("parent-energy", "#2563eb"),
        ("total-energy", "#0f172a"),
    )
    lower_lines = (
        ("intended-child-energy", "#16a34a"),
        ("cross-talk-energy", "#dc2626"),
        ("intended-grandchild-energy", "#7c3aed"),
        ("remainder-energy", "#d97706"),
    )
    polylines = "\n".join(
        f'  <polyline id="{name}" points="{points(name, upper_y)}" '
        f'fill="none" stroke="{color}" stroke-width="2.5"/>'
        for name, color in upper_lines
    )
    polylines += "\n" + "\n".join(
        f'  <polyline id="{name}" points="{points(name, lower_y)}" '
        f'fill="none" stroke="{color}" stroke-width="2.5"/>'
        for name, color in lower_lines
    )
    svg = f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <title>Two-stage partial-carrier Fourier-Galerkin shell energies</title>
  <rect width="{width}" height="{height}" fill="#ffffff"/>
  <g stroke="#cbd5e1" stroke-width="1" fill="none">
    <rect x="{left:.1f}" y="{upper_top:.1f}" width="{right-left:.1f}" height="{upper_bottom-upper_top:.1f}"/>
    <rect x="{left:.1f}" y="{lower_top:.1f}" width="{right-left:.1f}" height="{lower_bottom-lower_top:.1f}"/>
  </g>
  <g font-family="sans-serif" fill="#0f172a">
    <text x="{left:.1f}" y="32" font-size="18">Two-stage rejected carrier — pathway-contaminated numerical diagnostic</text>
    <text x="20" y="{(upper_top+upper_bottom)/2:.1f}" font-size="13" transform="rotate(-90 20 {(upper_top+upper_bottom)/2:.1f})">parent and total / initial energy</text>
    <text x="20" y="{(lower_top+lower_bottom)/2:.1f}" font-size="13" transform="rotate(-90 20 {(lower_top+lower_bottom)/2:.1f})">other shell / initial energy</text>
    <text x="{(left+right)/2:.1f}" y="558" font-size="13">scaled time N²t</text>
    <text x="{right-270:.1f}" y="84" font-size="12" fill="#2563eb">parents</text>
    <text x="{right-180:.1f}" y="84" font-size="12" fill="#0f172a">total</text>
    <text x="{right-400:.1f}" y="350" font-size="12" fill="#16a34a">intended children</text>
    <text x="{right-275:.1f}" y="350" font-size="12" fill="#dc2626">cross-talk</text>
    <text x="{right-190:.1f}" y="350" font-size="12" fill="#7c3aed">intended grandchildren</text>
    <text x="{right-45:.1f}" y="350" font-size="12" fill="#d97706">remainder</text>
  </g>
{polylines}
</svg>
"""
    return svg.encode("utf-8")


def _read_csv(path: Path) -> list[dict[str, str]]:
    try:
        stream = io.StringIO(path.read_text(encoding="utf-8"), newline="")
    except UnicodeDecodeError as error:
        raise ValueError(f"{path.name} is not UTF-8") from error
    reader = csv.DictReader(stream)
    if tuple(reader.fieldnames or ()) != ROW_FIELDS:
        raise ValueError(f"{path.name} has an invalid header")
    rows = list(reader)
    if any(None in row or set(row) != set(ROW_FIELDS) for row in rows):
        raise ValueError(f"{path.name} has malformed rows")
    return rows


def _csv_boolean(value: str, *, name: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ValueError(f"{name} is not a CSV boolean")


def _csv_float(value: str, *, name: str, optional: bool = False) -> float | None:
    if optional and value == "":
        return None
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{name} is not numeric") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} is not finite")
    return result


def _csv_integer(value: str, *, name: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise ValueError(f"{name} is not an integer") from error
    if str(result) != value:
        raise ValueError(f"{name} is not a canonical integer")
    return result


def _typed_fit_rows(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in rows:
        output.append(
            {
                "family": row["family"],
                "family_parameter": _csv_float(
                    row["family_parameter"], name="family_parameter"
                ),
                "gamma": _csv_float(row["gamma"], name="gamma"),
                "base_scale": _csv_integer(row["base_scale"], name="base_scale"),
                "core_fit_eligible": _csv_boolean(
                    row["core_fit_eligible"], name="core_fit_eligible"
                ),
                "relay_fit_eligible": _csv_boolean(
                    row["relay_fit_eligible"], name="relay_fit_eligible"
                ),
                "full_output_bands_disjoint": _csv_boolean(
                    row["full_output_bands_disjoint"],
                    name="full_output_bands_disjoint",
                ),
                "forcing_metric_method": row["forcing_metric_method"],
                **{
                    metric: (
                        float(_csv_integer(row[metric], name=metric))
                        if metric == "support_mode_count"
                        else _csv_float(row[metric], name=metric, optional=True)
                    )
                    for metric in FIT_METRICS
                },
            }
        )
    return output


def _verify_scaling_rows(
    rows: list[dict[str, str]],
    *,
    family: str,
    config: dict[str, Any],
) -> None:
    parameters = (
        [float(value) for value in config["gammas"]]
        if family == "power_width"
        else [float(value) for value in config["relative_width_fractions"]]
    )
    expected = {(parameter, int(scale)) for parameter in parameters for scale in config["scales"]}
    observed: set[tuple[float, int]] = set()
    energy_constant = float(config["energy_constant"])
    required_ratio = float(config["required_duhamel_ratio"])
    for row in rows:
        try:
            parsed_channel_norms = strict_json_loads(
                row["channel_norms_json"], label="channel norms"
            )
            strict_json_loads(
                row["box_overlap_counts_json"], label="box overlap counts"
            )
            strict_json_loads(
                row["channel_ratios_json"], label="channel ratios"
            )
        except ValueError as error:
            raise ValueError("embedded channel JSON is invalid") from error
        if row["family"] != family:
            raise ValueError("scaling CSV contains the wrong family")
        parameter = float(_csv_float(row["family_parameter"], name="family_parameter"))
        scale = _csv_integer(row["base_scale"], name="base_scale")
        observed.add((parameter, scale))
        critical_energy = float(
            _csv_float(row["critical_energy"], name="critical_energy")
        )
        if not math.isclose(
            critical_energy,
            energy_constant / scale,
            rel_tol=2.0e-13,
            abs_tol=1.0e-15,
        ):
            raise ValueError("critical energy is inconsistent with c_E/N")
        gamma = float(_csv_float(row["gamma"], name="gamma"))
        expected_gamma = parameter if family == "power_width" else 1.0
        if not math.isclose(gamma, expected_gamma, rel_tol=0.0, abs_tol=1.0e-15):
            raise ValueError("scaling row gamma does not match its family")
        width = _csv_integer(row["width"], name="width")
        expected_width = (
            max(1, math.floor(scale**parameter))
            if family == "power_width"
            else max(1, math.floor(scale * parameter))
        )
        if width != expected_width:
            raise ValueError("scaling row width is inconsistent")
        realized = float(
            _csv_float(row["realized_width_over_scale"], name="realized width")
        )
        if not math.isclose(realized, width / scale, rel_tol=1.0e-13, abs_tol=1.0e-15):
            raise ValueError("realized relative width is inconsistent")
        measurement = _csv_boolean(
            row["measurement_available"], name="measurement_available"
        )
        resource_limit = _csv_integer(
            row["local_fft_maximum_working_bytes"],
            name="local_fft_maximum_working_bytes",
        )
        if resource_limit != int(config["local_fft_maximum_working_bytes"]):
            raise ValueError("local FFT resource cap is inconsistent")
        core_fit_eligible = _csv_boolean(
            row["core_fit_eligible"], name="core_fit_eligible"
        )
        relay_fit_eligible = _csv_boolean(
            row["relay_fit_eligible"], name="relay_fit_eligible"
        )
        core_geometry_eligible = measurement and all(
            _csv_boolean(row[name], name=name)
            for name in (
                "parent_boxes_disjoint",
                "child_band_isolated",
                "channel_bands_disjoint",
            )
        )
        full_disjoint = _csv_boolean(
            row["full_output_bands_disjoint"], name="full_output_bands_disjoint"
        )
        if (
            core_fit_eligible != core_geometry_eligible
            or relay_fit_eligible != (core_geometry_eligible and full_disjoint)
        ):
            raise ValueError("core or relay fit eligibility is inconsistent")
        cancellation_available = _csv_boolean(
            row["energy_cancellation_available"],
            name="energy_cancellation_available",
        )
        cancellation = _csv_float(
            row["energy_cancellation_residual"],
            name="energy_cancellation_residual",
            optional=True,
        )
        if cancellation_available != (cancellation is not None):
            raise ValueError("energy cancellation availability is inconsistent")
        if cancellation is not None and not 0.0 <= cancellation <= 1.0e-10:
            raise ValueError("energy cancellation residual is too large")
        if measurement:
            local_input = _csv_integer(
                row["local_fft_input_side_length"],
                name="local_fft_input_side_length",
            )
            local_padded = _csv_integer(
                row["local_fft_padded_side_length"],
                name="local_fft_padded_side_length",
            )
            local_estimate = _csv_integer(
                row["local_fft_estimated_peak_working_bytes"],
                name="local_fft_estimated_peak_working_bytes",
            )
            if (
                row["measurement_method"]
                != "exact-zero-padded-local-fft-global-combination"
                or not _csv_boolean(
                    row["local_fft_zero_padding_alias_free"],
                    name="local_fft_zero_padding_alias_free",
                )
                or not _csv_boolean(
                    row["local_fft_global_overlap_aware"],
                    name="local_fft_global_overlap_aware",
                )
                or local_input != 2 * width - 1
                or local_padded != 2 * local_input - 1
                or local_estimate
                != 43 * local_padded**3 * np.dtype(np.complex128).itemsize
                or local_estimate > resource_limit
                or not row["forcing_metric_method"].startswith(
                    "exact-local-fft-full-global"
                )
                or not cancellation_available
                or row["energy_cancellation_method"]
                != "exact-zero-padded-local-fft-full-pairing"
            ):
                raise ValueError("local FFT backend metadata is inconsistent")
            parent_norm = float(
                _csv_float(row["parent_norm_squared"], name="parent_norm_squared")
            )
            if not math.isclose(
                parent_norm,
                2.0 * energy_constant / scale,
                rel_tol=2.0e-13,
                abs_tol=1.0e-15,
            ):
                raise ValueError("parent norm is inconsistent with 2c_E/N")
            support_count = _csv_integer(
                row["support_mode_count"], name="support_mode_count"
            )
            effective_count = float(
                _csv_float(row["effective_mode_count"], name="effective_mode_count")
            )
            if not 0.0 < effective_count <= support_count * (1.0 + 2.0e-13):
                raise ValueError("effective mode count exceeds support count")
            legacy_parent_divergence = float(
                _csv_float(row["divergence_relative"], name="divergence_relative")
            )
            legacy_parent_reality = float(
                _csv_float(row["reality_relative"], name="reality_relative")
            )
            parent_divergence = float(
                _csv_float(
                    row["parent_divergence_relative"],
                    name="parent_divergence_relative",
                )
            )
            parent_reality = float(
                _csv_float(
                    row["parent_reality_relative"],
                    name="parent_reality_relative",
                )
            )
            nonlinear_divergence = float(
                _csv_float(
                    row["nonlinear_divergence_relative"],
                    name="nonlinear_divergence_relative",
                )
            )
            nonlinear_reality = float(
                _csv_float(
                    row["nonlinear_reality_relative"],
                    name="nonlinear_reality_relative",
                )
            )
            if not (
                math.isclose(
                    legacy_parent_divergence,
                    parent_divergence,
                    rel_tol=0.0,
                    abs_tol=1.0e-18,
                )
                and math.isclose(
                    legacy_parent_reality,
                    parent_reality,
                    rel_tol=0.0,
                    abs_tol=1.0e-18,
                )
            ):
                raise ValueError("legacy structural residual is not its parent alias")
            structural_defects = [
                float(
                    _csv_float(
                        row["normalization_relative_error"],
                        name="normalization_relative_error",
                    )
                ),
                parent_divergence,
                parent_reality,
                nonlinear_divergence,
                nonlinear_reality,
            ]
            full_nonlinear_norm = float(
                _csv_float(row["full_nonlinear_norm"], name="full_nonlinear_norm")
            )
            if max(structural_defects) > 1.0e-10 or full_nonlinear_norm <= 0.0:
                raise ValueError("local FFT structural residual screen failed")
            if not isinstance(parsed_channel_norms, dict) or not {
                "child",
                "child_spill",
            }.issubset(parsed_channel_norms):
                raise ValueError("measured row lacks child channel norms")
            child_norm = float(parsed_channel_norms["child"])
            spill_norm = float(parsed_channel_norms["child_spill"])
            if child_norm <= 0.0 or min(spill_norm, full_nonlinear_norm) < 0.0:
                raise ValueError("child channel norm is invalid")
            child_fraction = float(
                _csv_float(
                    row["child_forcing_fraction"], name="child_forcing_fraction"
                )
            )
            expected_child_fraction = child_norm**2 / full_nonlinear_norm**2
            if not math.isclose(
                child_fraction,
                expected_child_fraction,
                rel_tol=2.0e-11,
                abs_tol=1.0e-14,
            ):
                raise ValueError("child forcing fraction is inconsistent")
            cancellation_pairing = float(
                _csv_float(
                    row["energy_cancellation_pairing"],
                    name="energy_cancellation_pairing",
                )
            )
            expected_cancellation = abs(cancellation_pairing) / max(
                math.sqrt(parent_norm) * full_nonlinear_norm,
                np.finfo(float).tiny,
            )
            if cancellation is None or not math.isclose(
                cancellation,
                expected_cancellation,
                rel_tol=2.0e-11,
                abs_tol=1.0e-18,
            ):
                raise ValueError("energy cancellation residual is inconsistent")
            expected_spill = spill_norm / child_norm
            expected_outside = math.sqrt(
                max(
                    full_nonlinear_norm**2
                    - child_norm**2
                    - spill_norm**2,
                    0.0,
                )
                / child_norm**2
            )
            expected_off_core = math.sqrt(
                max(full_nonlinear_norm**2 - child_norm**2, 0.0)
                / child_norm**2
            )
            recorded_off = float(
                _csv_float(row["off_chain_main_ratio"], name="off_chain_main_ratio")
            )
            recorded_alias = float(
                _csv_float(row["off_core_main_ratio"], name="off_core_main_ratio")
            )
            recorded_spill = float(
                _csv_float(
                    row["child_spill_main_ratio"], name="child_spill_main_ratio"
                )
            )
            recorded_outside = float(
                _csv_float(
                    row["outside_child_full_main_ratio"],
                    name="outside_child_full_main_ratio",
                )
            )
            if not all(
                math.isclose(observed, expected, rel_tol=2.0e-11, abs_tol=1.0e-14)
                for observed, expected in (
                    (recorded_off, expected_off_core),
                    (recorded_alias, expected_off_core),
                    (recorded_spill, expected_spill),
                    (recorded_outside, expected_outside),
                )
            ) or not math.isclose(
                recorded_off**2,
                recorded_spill**2 + recorded_outside**2,
                rel_tol=2.0e-11,
                abs_tol=1.0e-14,
            ):
                raise ValueError("off-core forcing decomposition failed")
            heat = float(_csv_float(row["heat_factor_h"], name="heat_factor_h"))
            if not 0.0 < heat <= float(config["tau"]) ** 2 * (1.0 + 2.0e-13):
                raise ValueError("heat factor lies outside 0 < H <= tau^2")
            gain = float(_csv_float(row["gain_g"], name="gain_g"))
            duhamel = float(
                _csv_float(row["duhamel_energy_ratio"], name="duhamel_energy_ratio")
            )
            amplitude = float(
                _csv_float(
                    row["duhamel_amplitude_ratio"], name="duhamel_amplitude_ratio"
                )
            )
            rhs = heat * gain * gain * parent_norm
            stored_rhs = float(
                _csv_float(row["duhamel_identity_rhs"], name="duhamel_identity_rhs")
            )
            if not (
                math.isclose(duhamel, rhs, rel_tol=2.0e-11, abs_tol=1.0e-14)
                and math.isclose(stored_rhs, rhs, rel_tol=2.0e-11, abs_tol=1.0e-14)
                and math.isclose(amplitude * amplitude, duhamel, rel_tol=2.0e-11, abs_tol=1.0e-14)
            ):
                raise ValueError("frozen-parent Duhamel identity failed")
            kappa = float(
                _csv_float(
                    row["support_radius_factor_kappa"],
                    name="support_radius_factor_kappa",
                )
            )
            universal = float(
                _csv_float(
                    row["universal_duhamel_upper_bound"],
                    name="universal_duhamel_upper_bound",
                )
            )
            effective = float(
                _csv_float(
                    row["effective_duhamel_upper_bound"],
                    name="effective_duhamel_upper_bound",
                )
            )
            to_universal = float(
                _csv_float(
                    row["duhamel_to_universal_upper_bound"],
                    name="duhamel_to_universal_upper_bound",
                )
            )
            to_effective = float(
                _csv_float(
                    row["duhamel_to_effective_upper_bound"],
                    name="duhamel_to_effective_upper_bound",
                )
            )
            bound_prefactor = (
                2.0
                * kappa**2
                * float(config["tau"]) ** 2
                * energy_constant
                / scale**3
            )
            if (
                row["duhamel_upper_bound_method"]
                != "phase-independent-bernstein-heat-support-and-effective-count"
                or min(kappa, universal, effective) <= 0.0
                or duhamel > effective * (1.0 + 2.0e-11) + 1.0e-14
                or effective > universal * (1.0 + 2.0e-11) + 1.0e-14
                or not math.isclose(
                    universal,
                    bound_prefactor * support_count,
                    rel_tol=2.0e-11,
                    abs_tol=1.0e-14,
                )
                or not math.isclose(
                    effective,
                    bound_prefactor * effective_count,
                    rel_tol=2.0e-11,
                    abs_tol=1.0e-14,
                )
                or not math.isclose(
                    to_universal,
                    duhamel / universal,
                    rel_tol=2.0e-11,
                    abs_tol=1.0e-14,
                )
                or not math.isclose(
                    to_effective,
                    duhamel / effective,
                    rel_tol=2.0e-11,
                    abs_tol=1.0e-14,
                )
            ):
                raise ValueError("phase-independent Duhamel upper bound failed")
            required_energy = float(
                _csv_float(
                    row["required_energy_constant_for_target"],
                    name="required_energy_constant_for_target",
                )
            )
            expected_energy = energy_constant * required_ratio / duhamel
            if not math.isclose(
                required_energy, expected_energy, rel_tol=2.0e-11, abs_tol=1.0e-14
            ):
                raise ValueError("required energy constant is inconsistent")
        else:
            for name in (
                "effective_mode_count",
                "parent_norm_squared",
                "gain_g",
                "heat_factor_h",
                "duhamel_energy_ratio",
                "duhamel_amplitude_ratio",
                "child_forcing_fraction",
                "off_chain_main_ratio",
                "off_core_main_ratio",
                "child_spill_main_ratio",
                "outside_child_full_main_ratio",
                "energy_cancellation_residual",
                "energy_cancellation_pairing",
                "full_nonlinear_norm",
                "normalization_relative_error",
                "divergence_relative",
                "reality_relative",
                "parent_divergence_relative",
                "parent_reality_relative",
                "nonlinear_divergence_relative",
                "nonlinear_reality_relative",
                "support_radius_factor_kappa",
                "universal_duhamel_upper_bound",
                "effective_duhamel_upper_bound",
                "duhamel_to_universal_upper_bound",
                "duhamel_to_effective_upper_bound",
            ):
                if row[name] != "":
                    raise ValueError("unavailable measurement contains fabricated values")
            for name in (
                "local_fft_input_side_length",
                "local_fft_padded_side_length",
                "local_fft_estimated_peak_working_bytes",
            ):
                if row[name] != "":
                    raise ValueError("unavailable row claims a local FFT allocation")
            if (
                _csv_boolean(
                    row["local_fft_zero_padding_alias_free"],
                    name="local_fft_zero_padding_alias_free",
                )
                or _csv_boolean(
                    row["local_fft_global_overlap_aware"],
                    name="local_fft_global_overlap_aware",
                )
                or row["measurement_method"] != ""
            ):
                raise ValueError("unavailable row claims a local FFT measurement")
            if row["duhamel_upper_bound_method"] != "unavailable-parent-overlap":
                raise ValueError("unavailable Duhamel bound has a dishonest method")
        asymptotic_pass = _csv_boolean(
            row["asymptotic_screen_pass"], name="asymptotic_screen_pass"
        )
        reasons = row["asymptotic_screen_reasons"].split("|")
        if family == "power_width" and (
            gamma >= 1.0
            or asymptotic_pass
            or "generic_duhamel_ratio_decays_for_gamma_below_one" not in reasons
        ):
            raise ValueError("sublinear gamma was not analytically rejected")
    if observed != expected or len(rows) != len(expected):
        raise ValueError("scaling CSV does not contain the exact configured grid")


def _verify_one_stage_comparison(comparison: object) -> dict[str, Any]:
    if not isinstance(comparison, dict):
        raise ValueError("one-stage Galerkin comparison must be an object")
    required = {
        "status",
        "scale",
        "width",
        "grid_size",
        "galerkin_cutoff",
        "dealias_margin",
        "steps",
        "tau",
        "viscosity",
        "parent_energy",
        "initial_child_energy",
        "child_core_half_width",
        "child_full_sumset_half_width",
        "child_core_mode_count",
        "child_full_sumset_mode_count",
        "child_core_forcing_fraction",
        "frozen_child_energy_ratio",
        "full_child_energy_ratio",
        "full_to_frozen_ratio",
        "final_total_energy_ratio",
        "maximum_energy_increase",
        "final_reality_defect",
        "final_divergence_defect",
        "history",
    }
    if set(comparison) != required or comparison["status"] != (
        "SMALL BINARY64 GALERKIN CROSS-CHECK / NOT CONTINUUM"
    ):
        raise ValueError("Galerkin comparison fields are invalid")
    numeric = {
        name: float(comparison[name])
        for name in (
            "parent_energy",
            "initial_child_energy",
            "child_core_forcing_fraction",
            "frozen_child_energy_ratio",
            "full_child_energy_ratio",
            "full_to_frozen_ratio",
            "final_total_energy_ratio",
            "maximum_energy_increase",
            "final_reality_defect",
            "final_divergence_defect",
        )
    }
    if any(not math.isfinite(value) for value in numeric.values()):
        raise ValueError("Galerkin comparison contains non-finite values")
    if (
        abs(numeric["initial_child_energy"]) > 1.0e-15
        or numeric["frozen_child_energy_ratio"] <= 0.0
        or numeric["full_child_energy_ratio"] <= 0.0
        or not math.isclose(
            numeric["full_to_frozen_ratio"],
            numeric["full_child_energy_ratio"] / numeric["frozen_child_energy_ratio"],
            rel_tol=2.0e-12,
        )
        or not 0.0 < numeric["child_core_forcing_fraction"] <= 1.0
        or numeric["final_total_energy_ratio"] > 1.0 + 1.0e-10
        or numeric["maximum_energy_increase"] > 1.0e-10
        or numeric["final_reality_defect"] > 1.0e-10
        or numeric["final_divergence_defect"] > 1.0e-10
        or int(comparison["dealias_margin"]) <= 0
        or int(comparison["child_core_half_width"]) != int(comparison["width"]) - 1
        or int(comparison["child_full_sumset_half_width"])
        != 2 * int(comparison["child_core_half_width"])
        or int(comparison["child_full_sumset_mode_count"])
        < int(comparison["child_core_mode_count"])
    ):
        raise ValueError("Galerkin invariant screen failed")
    history = comparison["history"]
    if not isinstance(history, list) or len(history) != int(comparison["steps"]) + 1:
        raise ValueError("Galerkin history length is invalid")
    previous_time = -math.inf
    previous_total = math.inf
    for record in history:
        if not isinstance(record, dict) or set(record) != {
            "time",
            "child_energy_ratio",
            "total_energy_ratio",
        }:
            raise ValueError("Galerkin history record is invalid")
        time = float(record["time"])
        child = float(record["child_energy_ratio"])
        total = float(record["total_energy_ratio"])
        if (
            not all(math.isfinite(value) for value in (time, child, total))
            or time <= previous_time
            or child < 0.0
            or total > previous_total + 1.0e-10
        ):
            raise ValueError("Galerkin history invariant failed")
        previous_time, previous_total = time, total
    if not (
        math.isclose(
            float(history[-1]["child_energy_ratio"]),
            numeric["full_child_energy_ratio"],
            rel_tol=2.0e-12,
        )
        and math.isclose(
            float(history[-1]["total_energy_ratio"]),
            numeric["final_total_energy_ratio"],
            rel_tol=2.0e-12,
        )
    ):
        raise ValueError("Galerkin final history record is inconsistent")
    return comparison


def _verify_two_stage_comparison(comparison: object) -> dict[str, Any]:
    if not isinstance(comparison, dict):
        raise ValueError("two-stage Galerkin comparison must be an object")
    required = {
        "status",
        "classification",
        "interpretation",
        "scale",
        "grid_size",
        "galerkin_cutoff",
        "dealias_margin",
        "dealias_verified",
        "steps",
        "tau",
        "time_multiples",
        "final_time",
        "viscosity",
        "energy_constant",
        "initial_parent_energy",
        "initial_named_nonparent_energy",
        "mode_groups",
        "final_parent_energy_ratio",
        "final_first_child_energy_ratio",
        "final_cross_talk_energy_ratio",
        "final_grandchild_energy_ratio",
        "final_remainder_energy_ratio",
        "final_total_energy_ratio",
        "maximum_energy_increase",
        "maximum_shell_budget_residual",
        "initial_energy_identity_residual",
        "final_energy_identity_residual",
        "final_reality_defect",
        "final_divergence_defect",
        "initial_grandchild_rhs_noise_ratio",
        "grandchild_roundoff_floor_ratio",
        "grandchild_resolution_margin",
        "grandchild_binary64_resolved",
        "history",
    }
    if set(comparison) != required:
        raise ValueError("two-stage Galerkin fields are invalid")
    if (
        comparison["classification"] != "partial_rejected_cross_talk"
        or comparison["dealias_verified"] is not True
        or comparison["grandchild_binary64_resolved"] is not True
        or "pathway-contaminated" not in comparison["interpretation"]
        or "NOT A RELAY" not in comparison["status"]
    ):
        raise ValueError("two-stage Galerkin classification is invalid")
    numeric_names = (
        "initial_parent_energy",
        "initial_named_nonparent_energy",
        "final_parent_energy_ratio",
        "final_first_child_energy_ratio",
        "final_cross_talk_energy_ratio",
        "final_grandchild_energy_ratio",
        "final_remainder_energy_ratio",
        "final_total_energy_ratio",
        "maximum_energy_increase",
        "maximum_shell_budget_residual",
        "initial_energy_identity_residual",
        "final_energy_identity_residual",
        "final_reality_defect",
        "final_divergence_defect",
        "initial_grandchild_rhs_noise_ratio",
        "grandchild_roundoff_floor_ratio",
        "grandchild_resolution_margin",
    )
    numeric = {name: float(comparison[name]) for name in numeric_names}
    if any(not math.isfinite(value) for value in numeric.values()):
        raise ValueError("two-stage Galerkin comparison contains non-finite values")
    if (
        numeric["initial_parent_energy"] <= 0.0
        or abs(numeric["initial_named_nonparent_energy"]) > 1.0e-15
        or numeric["final_first_child_energy_ratio"] <= 0.0
        or numeric["final_cross_talk_energy_ratio"] <= 0.0
        or numeric["final_grandchild_energy_ratio"] <= 0.0
        or numeric["final_cross_talk_energy_ratio"]
        < 0.9 * numeric["final_first_child_energy_ratio"]
        or numeric["final_total_energy_ratio"] > 1.0 + 1.0e-10
        or numeric["maximum_energy_increase"] > 1.0e-10
        or numeric["maximum_shell_budget_residual"] > 1.0e-10
        or numeric["initial_energy_identity_residual"] > 1.0e-10
        or numeric["final_energy_identity_residual"] > 1.0e-10
        or numeric["final_reality_defect"] > 1.0e-10
        or numeric["final_divergence_defect"] > 1.0e-10
        or numeric["grandchild_resolution_margin"] < 100.0
        or int(comparison["dealias_margin"]) <= 0
    ):
        raise ValueError("two-stage Galerkin invariant screen failed")
    groups = comparison["mode_groups"]
    expected_groups = {
        "parents",
        "first_child_one",
        "first_child_two",
        "cross_talk_one",
        "cross_talk_two",
        "grandchild_sum",
        "grandchild_difference",
    }
    if not isinstance(groups, dict) or set(groups) != expected_groups:
        raise ValueError("two-stage carrier mode groups are invalid")
    history = comparison["history"]
    if not isinstance(history, list) or len(history) != int(comparison["steps"]) + 1:
        raise ValueError("two-stage Galerkin history length is invalid")
    history_fields = {
        "time",
        "scaled_time",
        "parent_energy_ratio",
        "first_child_one_energy_ratio",
        "first_child_two_energy_ratio",
        "first_child_energy_ratio",
        "cross_talk_one_energy_ratio",
        "cross_talk_two_energy_ratio",
        "cross_talk_energy_ratio",
        "grandchild_sum_energy_ratio",
        "grandchild_difference_energy_ratio",
        "grandchild_energy_ratio",
        "remainder_energy_ratio",
        "total_energy_ratio",
    }
    previous_time = -math.inf
    previous_total = math.inf
    for record in history:
        if not isinstance(record, dict) or set(record) != history_fields:
            raise ValueError("two-stage Galerkin history record is invalid")
        values = {name: float(record[name]) for name in history_fields}
        if any(not math.isfinite(value) for value in values.values()):
            raise ValueError("two-stage Galerkin history is non-finite")
        if values["time"] <= previous_time or values["total_energy_ratio"] > previous_total + 1.0e-10:
            raise ValueError("two-stage Galerkin time or energy monotonicity failed")
        if not (
            math.isclose(
                values["first_child_energy_ratio"],
                values["first_child_one_energy_ratio"]
                + values["first_child_two_energy_ratio"],
                rel_tol=2.0e-11,
                abs_tol=2.0e-14,
            )
            and math.isclose(
                values["cross_talk_energy_ratio"],
                values["cross_talk_one_energy_ratio"]
                + values["cross_talk_two_energy_ratio"],
                rel_tol=2.0e-11,
                abs_tol=2.0e-14,
            )
            and math.isclose(
                values["grandchild_energy_ratio"],
                values["grandchild_sum_energy_ratio"]
                + values["grandchild_difference_energy_ratio"],
                rel_tol=2.0e-11,
                abs_tol=2.0e-14,
            )
            and math.isclose(
                values["total_energy_ratio"],
                values["parent_energy_ratio"]
                + values["first_child_energy_ratio"]
                + values["cross_talk_energy_ratio"]
                + values["grandchild_energy_ratio"]
                + values["remainder_energy_ratio"],
                rel_tol=2.0e-11,
                abs_tol=2.0e-14,
            )
        ):
            raise ValueError("two-stage shell-energy budget failed")
        previous_time = values["time"]
        previous_total = values["total_energy_ratio"]
    final = history[-1]
    for payload_name, history_name in (
        ("final_parent_energy_ratio", "parent_energy_ratio"),
        ("final_first_child_energy_ratio", "first_child_energy_ratio"),
        ("final_cross_talk_energy_ratio", "cross_talk_energy_ratio"),
        ("final_grandchild_energy_ratio", "grandchild_energy_ratio"),
        ("final_remainder_energy_ratio", "remainder_energy_ratio"),
        ("final_total_energy_ratio", "total_energy_ratio"),
    ):
        if not math.isclose(
            float(comparison[payload_name]),
            float(final[history_name]),
            rel_tol=2.0e-12,
            abs_tol=2.0e-14,
        ):
            raise ValueError("two-stage final history record is inconsistent")
    return comparison


def _verify_galerkin(payload: object) -> tuple[dict[str, Any], dict[str, Any]]:
    if not isinstance(payload, dict) or set(payload) != {
        "schema",
        "status",
        "claim_boundary",
        "one_stage_mesoscopic",
        "two_stage_carrier",
    }:
        raise ValueError("Galerkin payload has an invalid schema")
    if payload["schema"] != GALERKIN_SCHEMA or payload["status"] != STATUS:
        raise ValueError("Galerkin payload schema or status is invalid")
    return (
        _verify_one_stage_comparison(payload["one_stage_mesoscopic"]),
        _verify_two_stage_comparison(payload["two_stage_carrier"]),
    )


def verify_mesoscopic_leray_bundle(output_dir: Path) -> dict[str, Any]:
    """Verify hashes, current sources, provenance, and mathematical screens."""

    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        raise ValueError("bundle directory is missing")
    data_names = (*PAYLOAD_NAMES, "manifest.json")
    expected_names = set(data_names)
    expected_names.update(digest_sidecar(Path(name)).name for name in data_names)
    if {path.name for path in output_dir.iterdir()} != expected_names:
        raise ValueError("bundle file set is not exact")
    for name in data_names:
        verify_digest(output_dir / name)

    manifest = strict_json_loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8"),
        label="mesoscopic manifest",
    )
    manifest_fields = {
        "schema",
        "status",
        "seed",
        "source_config",
        "source_files",
        "files",
    }
    if not isinstance(manifest, dict) or set(manifest) != manifest_fields:
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
        if record["bytes"] != path.stat().st_size or record["sha256"] != sha256_file(path):
            raise ValueError(f"manifest checksum mismatch for {name}")
    source_files = manifest["source_files"]
    if not isinstance(source_files, dict) or set(source_files) != set(SOURCE_PATHS):
        raise ValueError("manifest source inventory is invalid")
    for name, digest in source_files.items():
        source_path = REPOSITORY_ROOT / name
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or not source_path.is_file()
            or sha256_file(source_path) != digest
        ):
            raise ValueError(f"source checksum mismatch for {name}")

    config = strict_json_loads(
        (output_dir / "config.snapshot.json").read_text(encoding="utf-8"),
        label="mesoscopic config snapshot",
    )
    config = _validate_config(config)
    source_config = manifest["source_config"]
    if not isinstance(source_config, dict) or set(source_config) != {
        "requested_path",
        "snapshot_sha256",
    }:
        raise ValueError("manifest source_config is invalid")
    if (
        not isinstance(source_config["requested_path"], str)
        or not source_config["requested_path"]
    ):
        raise ValueError("manifest source config path is invalid")
    if source_config["snapshot_sha256"] != sha256_file(
        output_dir / "config.snapshot.json"
    ):
        raise ValueError("config snapshot checksum is inconsistent")
    if manifest["seed"] != config["seed"]:
        raise ValueError("manifest seed does not match config")

    carrier = strict_json_loads(
        (output_dir / "exact_carrier_certificate.json").read_text(encoding="utf-8"),
        label="carrier certificate",
    )
    recomputed_carrier = search_exact_carrier_gadget().as_dict()
    if carrier != recomputed_carrier:
        raise ValueError("exact carrier certificate does not recompute")
    if (
        carrier.get("classification") != "partial_rejected_cross_talk"
        or carrier.get("strict_search", {}).get("found") is not False
        or carrier.get("strict_search", {}).get("exhaustive") is not True
    ):
        raise ValueError("carrier certificate classification is invalid")
    independent_carrier_verification = (
        verify_serialized_strict_orientation_records(carrier)
    )
    if (
        independent_carrier_verification.get("verified") is not True
        or independent_carrier_verification.get("records_verified") != 16
        or independent_carrier_verification.get("strict_pass_count") != 0
    ):
        raise ValueError("independent strict carrier record verification failed")

    power_rows = _read_csv(output_dir / "mesoscopic_scaling.csv")
    relative_rows = _read_csv(output_dir / "relative_width_scaling.csv")
    _verify_scaling_rows(power_rows, family="power_width", config=config)
    _verify_scaling_rows(relative_rows, family="fixed_relative_width", config=config)
    fits = strict_json_loads(
        (output_dir / "scaling_fits.json").read_text(encoding="utf-8"),
        label="scaling fits",
    )
    recomputed_fits = _build_scaling_fits(
        _typed_fit_rows(power_rows + relative_rows)
    )
    if fits != recomputed_fits:
        raise ValueError("scaling fits do not use exactly the eligible last-four windows")

    galerkin = strict_json_loads(
        (output_dir / "galerkin_comparison.json").read_text(encoding="utf-8"),
        label="Galerkin comparison",
    )
    one_stage_comparison, two_stage_comparison = _verify_galerkin(galerkin)
    svg = (output_dir / "shell_energy.svg").read_bytes()
    if not (
        svg.startswith(b"<?xml version=\"1.0\"")
        and b'id="parent-energy"' in svg
        and b'id="total-energy"' in svg
        and b'id="intended-child-energy"' in svg
        and b'id="cross-talk-energy"' in svg
        and b'id="intended-grandchild-energy"' in svg
        and b'id="remainder-energy"' in svg
    ):
        raise ValueError("shell-energy SVG is invalid")

    summary = strict_json_loads(
        (output_dir / "summary.json").read_text(encoding="utf-8"),
        label="mesoscopic summary",
    )
    required_summary = {
        "schema",
        "status",
        "claim_boundary",
        "measurement_backend",
        "critical_scaling_identity",
        "table_counts",
        "sublinear_verdict",
        "shape_thickening_obstruction",
        "relative_width_boundary",
        "two_stage_comparison",
        "small_galerkin",
        "two_stage_galerkin",
        "provenance",
    }
    if not isinstance(summary, dict) or set(summary) != required_summary:
        raise ValueError("summary has missing or unknown fields")
    if summary["schema"] != OUTPUT_SCHEMA or summary["status"] != STATUS:
        raise ValueError("summary schema or status is invalid")
    if summary["measurement_backend"] != {
        "name": "exact-zero-padded-local-fft-global-combination",
        "finite_arithmetic": "complex binary64",
        "linear_convolution_padding": "L=2*W-1; K=2*L-1",
        "global_overlap_handling": "combine coefficients before norms and tags",
        "maximum_working_bytes": int(config["local_fft_maximum_working_bytes"]),
        "claim_boundary": "exact finite convolution algorithm; no continuum enclosure",
    }:
        raise ValueError("summary local FFT backend is inconsistent")
    provenance = validate_runtime_provenance(summary["provenance"])
    fingerprint, file_count, scope = source_fingerprint(REPOSITORY_ROOT)
    if (
        provenance["source_fingerprint_sha256"] != fingerprint
        or provenance["source_fingerprint_file_count"] != file_count
        or provenance["source_fingerprint_scope"] != scope
    ):
        raise ValueError("runtime provenance does not match current sources")
    if summary["table_counts"] != {
        "mesoscopic_rows": len(power_rows),
        "relative_width_rows": len(relative_rows),
        "mesoscopic_measurements_available": sum(
            _csv_boolean(row["measurement_available"], name="measurement_available")
            for row in power_rows
        ),
        "relative_width_measurements_available": sum(
            _csv_boolean(row["measurement_available"], name="measurement_available")
            for row in relative_rows
        ),
    }:
        raise ValueError("summary table counts are inconsistent")
    if summary["sublinear_verdict"].get("classification") != (
        "REJECTED_FOR_EVERY_CONFIGURED_GAMMA_BELOW_ONE"
    ) or summary["sublinear_verdict"].get("effective_count_necessity") != (
        "a scale-independent child-energy fraction requires order N^3 "
        "effective parent modes when kappa_N is uniformly bounded"
    ) or summary["sublinear_verdict"].get("configured_gammas") != [
        float(value) for value in config["gammas"]
    ]:
        raise ValueError("summary sublinear verdict is invalid")
    if summary["critical_scaling_identity"].get(
        "phase_independent_upper_bound"
    ) != (
        "D_N <= 2*kappa_N^2*tau^2*c_E*M_eff,N/N^3 "
        "<= 2*kappa_N^2*tau^2*c_E*M_N/N^3"
    ):
        raise ValueError("summary omits the phase-independent Duhamel bound")
    if summary["shape_thickening_obstruction"] != {
        "classification": "ORIGINAL_OFF_CORE_SCALING_LABEL_REJECTED",
        "exact_decomposition": (
            "off_core_main^2=child_spill_main^2+outside_child_full_main^2"
        ),
        "reason": (
            "off-core forcing contains geometric child spill, which need not "
            "decay as N^(gamma-1) for a filled Fejer box"
        ),
        "fit_policy": (
            "the gamma-1 comparison is attached only to "
            "outside_child_full_main_ratio"
        ),
    }:
        raise ValueError("summary shape-thickening obstruction is inconsistent")
    if (
        summary["relative_width_boundary"].get("width_rule")
        != "max(1,floor(rho*N))"
        or summary["relative_width_boundary"].get("fractions")
        != [float(value) for value in config["relative_width_fractions"]]
    ):
        raise ValueError("summary relative-width boundary is inconsistent")
    partial = carrier["partial_gadget"]
    strict = carrier["strict_search"]
    expected_two_stage = {
        "sequential_partial_gadget": {
            "status": "REJECTED_CROSS_TALK_AND_SIMPLE_RECURSION",
            "nonzero_next_outputs": True,
            "intended_to_full_cross_power": partial[
                "intended_fraction_of_relay_cross_power"
            ],
            "simple_binary_recursion_rejected": partial[
                "simple_binary_recursion_rejected"
            ],
        },
        "joint_strict_finite_search": {
            "status": "INFEASIBLE_IN_STATED_FINITE_ALPHABET",
            "exhaustive": strict["exhaustive"],
            "orientations_tested": strict[
                "eligible_second_relay_orientations_tested"
            ],
            "found": strict["found"],
            "joint_objective_j_n": None,
            "reason": "no admissible graph, so no optimization score is manufactured",
        },
    }
    if summary["two_stage_comparison"] != expected_two_stage:
        raise ValueError("summary two-stage comparison is inconsistent")
    if summary["small_galerkin"] != {
        "status": one_stage_comparison["status"],
        "frozen_child_energy_ratio": one_stage_comparison[
            "frozen_child_energy_ratio"
        ],
        "full_child_energy_ratio": one_stage_comparison["full_child_energy_ratio"],
        "full_to_frozen_ratio": one_stage_comparison["full_to_frozen_ratio"],
        "dealias_margin": one_stage_comparison["dealias_margin"],
        "child_core_forcing_fraction": one_stage_comparison[
            "child_core_forcing_fraction"
        ],
    }:
        raise ValueError("summary Galerkin values are inconsistent")
    if summary["two_stage_galerkin"] != {
        "status": two_stage_comparison["status"],
        "classification": two_stage_comparison["classification"],
        "interpretation": two_stage_comparison["interpretation"],
        "final_first_child_energy_ratio": two_stage_comparison[
            "final_first_child_energy_ratio"
        ],
        "final_cross_talk_energy_ratio": two_stage_comparison[
            "final_cross_talk_energy_ratio"
        ],
        "final_grandchild_energy_ratio": two_stage_comparison[
            "final_grandchild_energy_ratio"
        ],
        "grandchild_binary64_resolved": two_stage_comparison[
            "grandchild_binary64_resolved"
        ],
        "pathway_contaminated": True,
    }:
        raise ValueError("summary two-stage Galerkin values are inconsistent")
    return {
        "verified": True,
        "schema": OUTPUT_SCHEMA,
        "payload_count": len(PAYLOAD_NAMES),
        "source_file_count": len(SOURCE_PATHS),
        "power_row_count": len(power_rows),
        "relative_row_count": len(relative_rows),
        "independent_carrier_verification": independent_carrier_verification,
    }


def _measure_tables(config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    viscosity = float(Fraction(config["viscosity"]))
    common = {
        "energy_constant": float(config["energy_constant"]),
        "viscosity": viscosity,
        "tau": float(config["tau"]),
        "exact_target_limit": int(config["exact_target_limit"]),
        "strata_per_axis": int(config["strata_per_axis"]),
        "exact_energy_pair_limit": int(config["exact_energy_pair_limit"]),
        "required_duhamel_ratio": float(config["required_duhamel_ratio"]),
    }
    power_rows: list[dict[str, Any]] = []
    for gamma in config["gammas"]:
        for scale in config["scales"]:
            cloud_config = MesoscopicCloudConfig(
                base_scale=int(scale), gamma=float(gamma), **common
            )
            metrics = measure_local_fft_mesoscopic_metrics(
                cloud_config,
                maximum_working_bytes=int(
                    config["local_fft_maximum_working_bytes"]
                ),
            )
            power_rows.append(
                _metric_row(
                    metrics,
                    family="power_width",
                    family_parameter=float(gamma),
                    relative_width_fraction=None,
                    energy_constant=float(config["energy_constant"]),
                    required_duhamel_ratio=float(config["required_duhamel_ratio"]),
                    local_fft_maximum_working_bytes=int(
                        config["local_fft_maximum_working_bytes"]
                    ),
                )
            )
    relative_rows: list[dict[str, Any]] = []
    for fraction in config["relative_width_fractions"]:
        for scale in config["scales"]:
            width = max(1, math.floor(float(fraction) * int(scale)))
            cloud_config = MesoscopicCloudConfig(
                base_scale=int(scale),
                gamma=1.0,
                width_override=width,
                **common,
            )
            metrics = measure_local_fft_mesoscopic_metrics(
                cloud_config,
                maximum_working_bytes=int(
                    config["local_fft_maximum_working_bytes"]
                ),
            )
            relative_rows.append(
                _metric_row(
                    metrics,
                    family="fixed_relative_width",
                    family_parameter=float(fraction),
                    relative_width_fraction=float(fraction),
                    energy_constant=float(config["energy_constant"]),
                    required_duhamel_ratio=float(config["required_duhamel_ratio"]),
                    local_fft_maximum_working_bytes=int(
                        config["local_fft_maximum_working_bytes"]
                    ),
                )
            )
    return power_rows, relative_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/mesoscopic_leray_cloud_v1.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/mesoscopic_leray_cloud_v1"),
    )
    args = parser.parse_args()
    config = _load_config(args.config)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("output directory must be new or empty")

    carrier = search_exact_carrier_gadget().as_dict()
    power_rows, relative_rows = _measure_tables(config)
    fits = _build_scaling_fits(power_rows + relative_rows)
    small = config["small_galerkin"]
    galerkin_result = run_small_mesoscopic_galerkin(
        scale=int(small["scale"]),
        width=int(small["width"]),
        grid_size=int(small["grid_size"]),
        viscosity=float(Fraction(config["viscosity"])),
        energy_constant=float(config["energy_constant"]),
        tau=float(small["tau"]),
        steps=int(small["steps"]),
    )
    one_stage_comparison = galerkin_result.as_dict()
    two_stage_config = config["carrier_two_stage"]
    two_stage_result = run_carrier_two_stage_galerkin(
        scale=int(two_stage_config["scale"]),
        grid_size=int(two_stage_config["grid_size"]),
        viscosity=float(Fraction(config["viscosity"])),
        energy_constant=float(config["energy_constant"]),
        tau=float(two_stage_config["tau"]),
        time_multiples=float(two_stage_config["time_multiples"]),
        steps=int(two_stage_config["steps"]),
    )
    two_stage_galerkin = two_stage_result.as_dict()
    galerkin_payload = {
        "schema": GALERKIN_SCHEMA,
        "status": STATUS,
        "claim_boundary": (
            "Two binary64 finite Fourier--Galerkin RK4 diagnostics: an empty "
            "mesoscopic child versus frozen parents, and a two-stage partial "
            "carrier shell budget. No continuum, precision, attribution, or "
            "infinite-time enclosure is claimed."
        ),
        "one_stage_mesoscopic": one_stage_comparison,
        "two_stage_carrier": two_stage_galerkin,
    }
    partial = carrier["partial_gadget"]
    strict = carrier["strict_search"]
    two_stage = {
        "sequential_partial_gadget": {
            "status": "REJECTED_CROSS_TALK_AND_SIMPLE_RECURSION",
            "nonzero_next_outputs": True,
            "intended_to_full_cross_power": partial[
                "intended_fraction_of_relay_cross_power"
            ],
            "simple_binary_recursion_rejected": partial[
                "simple_binary_recursion_rejected"
            ],
        },
        "joint_strict_finite_search": {
            "status": "INFEASIBLE_IN_STATED_FINITE_ALPHABET",
            "exhaustive": strict["exhaustive"],
            "orientations_tested": strict[
                "eligible_second_relay_orientations_tested"
            ],
            "found": strict["found"],
            "joint_objective_j_n": None,
            "reason": "no admissible graph, so no optimization score is manufactured",
        },
    }
    summary = {
        "schema": OUTPUT_SCHEMA,
        "status": STATUS,
        "claim_boundary": (
            "This bundle contains exact finite algebra and reproducible numerical "
            "falsification screens. It constructs no invariant cascade, PDE orbit, "
            "singularity, or proof of any Clay statement."
        ),
        "measurement_backend": {
            "name": "exact-zero-padded-local-fft-global-combination",
            "finite_arithmetic": "complex binary64",
            "linear_convolution_padding": "L=2*W-1; K=2*L-1",
            "global_overlap_handling": "combine coefficients before norms and tags",
            "maximum_working_bytes": int(
                config["local_fft_maximum_working_bytes"]
            ),
            "claim_boundary": (
                "exact finite convolution algorithm; no continuum enclosure"
            ),
        },
        "critical_scaling_identity": {
            "energy_law": "||u_N||_2^2=2*c_E/N",
            "gain": "G_N=A_N/(N^2*||u_N||_2^2)",
            "frozen_child_energy_ratio": "D_N=(2*c_E/N)*H_N*G_N^2",
            "phase_independent_upper_bound": (
                "D_N <= 2*kappa_N^2*tau^2*c_E*M_eff,N/N^3 "
                "<= 2*kappa_N^2*tau^2*c_E*M_N/N^3"
            ),
            "necessary_effective_mode_count": (
                "D_N>=delta requires M_eff,N >= "
                "delta*N^3/(2*kappa_N^2*tau^2*c_E)"
            ),
            "optimistic_sublinear_exponent": "3*gamma-3",
            "consequence": "G_N growth alone is not a relay success condition",
        },
        "table_counts": {
            "mesoscopic_rows": len(power_rows),
            "relative_width_rows": len(relative_rows),
            "mesoscopic_measurements_available": sum(
                bool(row["measurement_available"]) for row in power_rows
            ),
            "relative_width_measurements_available": sum(
                bool(row["measurement_available"]) for row in relative_rows
            ),
        },
        "sublinear_verdict": {
            "classification": "REJECTED_FOR_EVERY_CONFIGURED_GAMMA_BELOW_ONE",
            "configured_gammas": [float(value) for value in config["gammas"]],
            "reason": (
                "The phase-independent bound D_N <= "
                "2*kappa_N^2*tau^2*c_E*M_eff,N/N^3 and "
                "M_eff,N<=M_N~N^(3*gamma) imply D_N->0 for uniformly "
                "bounded kappa_N and gamma<1; the capacity exponent "
                "3*gamma-3 is the matching heuristic scaling."
            ),
            "effective_count_necessity": (
                "a scale-independent child-energy fraction requires order N^3 "
                "effective parent modes when kappa_N is uniformly bounded"
            ),
        },
        "shape_thickening_obstruction": {
            "classification": "ORIGINAL_OFF_CORE_SCALING_LABEL_REJECTED",
            "exact_decomposition": (
                "off_core_main^2=child_spill_main^2+outside_child_full_main^2"
            ),
            "reason": (
                "off-core forcing contains geometric child spill, which need not "
                "decay as N^(gamma-1) for a filled Fejer box"
            ),
            "fit_policy": (
                "the gamma-1 comparison is attached only to "
                "outside_child_full_main_ratio"
            ),
        },
        "relative_width_boundary": {
            "classification": "NOT_ANALYTICALLY_REJECTED_BUT_NOT_ESTABLISHED",
            "width_rule": "max(1,floor(rho*N))",
            "fractions": [float(value) for value in config["relative_width_fractions"]],
            "open_requirements": [
                "scale-independent constant large enough at fixed critical energy",
                "admissible two-stage carrier closure",
                "persistence under full PDE evolution",
            ],
        },
        "two_stage_comparison": two_stage,
        "small_galerkin": {
            "status": one_stage_comparison["status"],
            "frozen_child_energy_ratio": one_stage_comparison[
                "frozen_child_energy_ratio"
            ],
            "full_child_energy_ratio": one_stage_comparison[
                "full_child_energy_ratio"
            ],
            "full_to_frozen_ratio": one_stage_comparison[
                "full_to_frozen_ratio"
            ],
            "dealias_margin": one_stage_comparison["dealias_margin"],
            "child_core_forcing_fraction": one_stage_comparison[
                "child_core_forcing_fraction"
            ],
        },
        "two_stage_galerkin": {
            "status": two_stage_galerkin["status"],
            "classification": two_stage_galerkin["classification"],
            "interpretation": two_stage_galerkin["interpretation"],
            "final_first_child_energy_ratio": two_stage_galerkin[
                "final_first_child_energy_ratio"
            ],
            "final_cross_talk_energy_ratio": two_stage_galerkin[
                "final_cross_talk_energy_ratio"
            ],
            "final_grandchild_energy_ratio": two_stage_galerkin[
                "final_grandchild_energy_ratio"
            ],
            "grandchild_binary64_resolved": two_stage_galerkin[
                "grandchild_binary64_resolved"
            ],
            "pathway_contaminated": True,
        },
        "provenance": collect_runtime_provenance(REPOSITORY_ROOT),
    }

    files = {
        "config.snapshot.json": canonical_json_bytes(config),
        "exact_carrier_certificate.json": canonical_json_bytes(carrier),
        "mesoscopic_scaling.csv": _csv_bytes(power_rows),
        "relative_width_scaling.csv": _csv_bytes(relative_rows),
        "scaling_fits.json": canonical_json_bytes(fits),
        "galerkin_comparison.json": canonical_json_bytes(galerkin_payload),
        "shell_energy.svg": _svg_bytes(list(two_stage_result.history)),
        "summary.json": canonical_json_bytes(summary),
    }
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("output directory must remain new or empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in files.items():
        write_with_digest(args.output_dir / name, payload)
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "status": STATUS,
        "seed": int(config["seed"]),
        "source_config": {
            "requested_path": args.config.as_posix(),
            "snapshot_sha256": sha256_file(
                args.output_dir / "config.snapshot.json"
            ),
        },
        "source_files": {
            name: sha256_file(REPOSITORY_ROOT / name) for name in SOURCE_PATHS
        },
        "files": {
            name: {
                "bytes": (args.output_dir / name).stat().st_size,
                "sha256": sha256_file(args.output_dir / name),
            }
            for name in sorted(files)
        },
    }
    write_with_digest(
        args.output_dir / "manifest.json", canonical_json_bytes(manifest)
    )
    verify_mesoscopic_leray_bundle(args.output_dir)
    print(canonical_json_bytes(summary).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
