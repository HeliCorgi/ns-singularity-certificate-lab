"""Power-law fitting gates for candidate Navier--Stokes singularity scalings.

The routines in this module are diagnostic.  A good log-linear fit does not
prove a singularity.  The candidate singular time and fit window must be
reported, and any inferred exponents must still pass space/time/domain and
solver refinement.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable, Sequence

import numpy as np
import numpy.typing as npt

from .scaling_constraints import (
    CandidateFamily,
    PowerLawScaling,
    ScalingConstraintReport,
    analyze_power_law_scaling,
)

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class LogPowerLawFit:
    """Fit of ``value ~ prefactor * (T-t)^signed_exponent``."""

    exponent: float
    prefactor: float
    log_rmse: float
    r_squared: float
    point_count: int
    singular_time: float
    start_time: float
    end_time: float

    def to_dict(self) -> dict[str, float | int]:
        return asdict(self)


@dataclass(frozen=True)
class JointScalingFit:
    """Joint amplitude/radial/axial fit for one time window and singular time."""

    singular_time: float
    start_index: int
    point_count: int
    amplitude_fit: LogPowerLawFit
    radial_fit: LogPowerLawFit
    axial_fit: LogPowerLawFit
    total_log_mse: float
    constraint_report: ScalingConstraintReport | None

    @property
    def alpha(self) -> float:
        return -self.amplitude_fit.exponent

    @property
    def beta_r(self) -> float:
        return self.radial_fit.exponent

    @property
    def beta_z(self) -> float:
        return self.axial_fit.exponent

    @property
    def candidate_family(self) -> CandidateFamily:
        if self.constraint_report is None:
            return CandidateFamily.NONE
        return self.constraint_report.candidate_family

    @property
    def passes_conditional_core_gates(self) -> bool:
        return bool(
            self.constraint_report is not None
            and self.constraint_report.passes_conditional_core_gates
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "singular_time": self.singular_time,
            "start_index": self.start_index,
            "point_count": self.point_count,
            "alpha": self.alpha,
            "beta_r": self.beta_r,
            "beta_z": self.beta_z,
            "total_log_mse": self.total_log_mse,
            "candidate_family": self.candidate_family.value,
            "passes_conditional_core_gates": self.passes_conditional_core_gates,
            "amplitude_fit": self.amplitude_fit.to_dict(),
            "radial_fit": self.radial_fit.to_dict(),
            "axial_fit": self.axial_fit.to_dict(),
            "constraint_report": (
                None if self.constraint_report is None else self.constraint_report.to_dict()
            ),
        }


def _validated_series(
    times: npt.ArrayLike,
    values: npt.ArrayLike,
) -> tuple[FloatArray, FloatArray]:
    time_array = np.asarray(times, dtype=np.float64)
    value_array = np.asarray(values, dtype=np.float64)
    if time_array.ndim != 1 or value_array.ndim != 1:
        raise ValueError("times and values must be one-dimensional")
    if time_array.shape != value_array.shape:
        raise ValueError("times and values must have identical shape")
    if time_array.size < 2:
        raise ValueError("at least two points are required")
    if not np.all(np.isfinite(time_array)) or not np.all(np.isfinite(value_array)):
        raise ValueError("times and values must be finite")
    if np.any(np.diff(time_array) <= 0.0):
        raise ValueError("times must be strictly increasing")
    if np.any(value_array <= 0.0):
        raise ValueError("values must be strictly positive")
    return time_array, value_array


def fit_log_power_law(
    times: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    singular_time: float,
) -> LogPowerLawFit:
    """Fit ``value = C (T-t)^m`` in logarithmic coordinates."""

    time_array, value_array = _validated_series(times, values)
    candidate_time = float(singular_time)
    if not math.isfinite(candidate_time) or candidate_time <= float(time_array[-1]):
        raise ValueError("singular_time must be finite and greater than every sample time")

    x_values = np.log(candidate_time - time_array)
    y_values = np.log(value_array)
    design = np.column_stack((np.ones_like(x_values), x_values))
    coefficients, _, _, _ = np.linalg.lstsq(design, y_values, rcond=None)
    intercept = float(coefficients[0])
    exponent = float(coefficients[1])
    prediction = design @ coefficients
    residual = y_values - prediction
    mse = float(np.mean(residual * residual))
    centered = y_values - float(np.mean(y_values))
    total = float(np.sum(centered * centered))
    residual_sum = float(np.sum(residual * residual))
    r_squared = 1.0 if total == 0.0 and residual_sum == 0.0 else (
        0.0 if total == 0.0 else 1.0 - residual_sum / total
    )
    return LogPowerLawFit(
        exponent=exponent,
        prefactor=float(math.exp(intercept)),
        log_rmse=float(math.sqrt(mse)),
        r_squared=float(r_squared),
        point_count=int(time_array.size),
        singular_time=candidate_time,
        start_time=float(time_array[0]),
        end_time=float(time_array[-1]),
    )


def fit_joint_scaling(
    times: npt.ArrayLike,
    amplitudes: npt.ArrayLike,
    radial_scales: npt.ArrayLike,
    axial_scales: npt.ArrayLike,
    *,
    singular_time: float,
    start_index: int = 0,
    balance_tolerance: float = 2.5e-2,
) -> JointScalingFit:
    """Fit amplitude and two widths over one suffix window."""

    time_array, amplitude_array = _validated_series(times, amplitudes)
    time_r, radial_array = _validated_series(times, radial_scales)
    time_z, axial_array = _validated_series(times, axial_scales)
    if not np.array_equal(time_array, time_r) or not np.array_equal(time_array, time_z):
        raise ValueError("all series must use exactly the same sample times")
    if isinstance(start_index, bool) or start_index < 0 or start_index >= time_array.size - 1:
        raise ValueError("start_index must leave at least two points")

    window = slice(start_index, None)
    amplitude_fit = fit_log_power_law(
        time_array[window], amplitude_array[window], singular_time=singular_time
    )
    radial_fit = fit_log_power_law(
        time_array[window], radial_array[window], singular_time=singular_time
    )
    axial_fit = fit_log_power_law(
        time_array[window], axial_array[window], singular_time=singular_time
    )
    total_mse = (
        amplitude_fit.log_rmse**2
        + radial_fit.log_rmse**2
        + axial_fit.log_rmse**2
    )

    alpha = -amplitude_fit.exponent
    beta_r = radial_fit.exponent
    beta_z = axial_fit.exponent
    constraints: ScalingConstraintReport | None = None
    if alpha >= 0.0 and beta_r > 0.0 and beta_z > 0.0:
        constraints = analyze_power_law_scaling(
            PowerLawScaling(alpha=alpha, beta_r=beta_r, beta_z=beta_z),
            tolerance=balance_tolerance,
        )

    return JointScalingFit(
        singular_time=float(singular_time),
        start_index=int(start_index),
        point_count=int(time_array.size - start_index),
        amplitude_fit=amplitude_fit,
        radial_fit=radial_fit,
        axial_fit=axial_fit,
        total_log_mse=float(total_mse),
        constraint_report=constraints,
    )


def scan_joint_scaling_fits(
    times: npt.ArrayLike,
    amplitudes: npt.ArrayLike,
    radial_scales: npt.ArrayLike,
    axial_scales: npt.ArrayLike,
    *,
    singular_times: Iterable[float],
    minimum_points: int = 3,
    balance_tolerance: float = 2.5e-2,
) -> list[JointScalingFit]:
    """Fit every candidate time and suffix window, sorted by joint log MSE."""

    time_array, _ = _validated_series(times, amplitudes)
    if isinstance(minimum_points, bool) or minimum_points < 3:
        raise ValueError("minimum_points must be an integer at least three")
    if minimum_points > time_array.size:
        raise ValueError("minimum_points exceeds the available samples")

    fits: list[JointScalingFit] = []
    latest_start = int(time_array.size - minimum_points)
    for candidate_time in singular_times:
        value = float(candidate_time)
        if not math.isfinite(value) or value <= float(time_array[-1]):
            continue
        for start_index in range(latest_start + 1):
            fits.append(
                fit_joint_scaling(
                    time_array,
                    amplitudes,
                    radial_scales,
                    axial_scales,
                    singular_time=value,
                    start_index=start_index,
                    balance_tolerance=balance_tolerance,
                )
            )
    return sorted(fits, key=lambda fit: fit.total_log_mse)


def summarize_fit_scan(fits: Sequence[JointScalingFit]) -> dict[str, object]:
    """Summarize best unconstrained and best gate-passing fits."""

    if not fits:
        raise ValueError("at least one fit is required")
    best = min(fits, key=lambda fit: fit.total_log_mse)
    passing = [fit for fit in fits if fit.passes_conditional_core_gates]
    family_counts: dict[str, int] = {}
    for fit in passing:
        key = fit.candidate_family.value
        family_counts[key] = family_counts.get(key, 0) + 1
    return {
        "fit_count": len(fits),
        "passing_fit_count": len(passing),
        "passing_family_counts": family_counts,
        "best_unconstrained": best.to_dict(),
        "best_passing": (
            None if not passing else min(passing, key=lambda fit: fit.total_log_mse).to_dict()
        ),
        "warning": (
            "A fit is only an exponent diagnostic. It is not a PDE solution, and scan "
            "multiplicity makes post-selected fits especially vulnerable to false positives."
        ),
    }
