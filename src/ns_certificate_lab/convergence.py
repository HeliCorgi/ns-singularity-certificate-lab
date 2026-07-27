"""Resolution-convergence diagnostics with explicit pass/fail criteria."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class ConvergenceReport:
    passed: bool
    errors: tuple[float, ...]
    spacings: tuple[float, ...]
    observed_orders: tuple[float, ...]
    minimum_required_order: float
    strictly_decreasing: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "passed": self.passed,
            "errors": list(self.errors),
            "spacings": list(self.spacings),
            "observed_orders": list(self.observed_orders),
            "minimum_required_order": self.minimum_required_order,
            "strictly_decreasing": self.strictly_decreasing,
        }


def assess_convergence(
    errors: Sequence[float],
    spacings: Sequence[float],
    *,
    minimum_order: float = 1.5,
) -> ConvergenceReport:
    """Assess monotone convergence and every adjacent observed order."""

    error_values = np.asarray(errors, dtype=np.float64)
    h_values = np.asarray(spacings, dtype=np.float64)
    if error_values.ndim != 1 or h_values.ndim != 1:
        raise ValueError("errors and spacings must be one-dimensional")
    if error_values.size < 3 or error_values.size != h_values.size:
        raise ValueError("at least three matching errors and spacings are required")
    if (
        not np.all(np.isfinite(error_values))
        or not np.all(np.isfinite(h_values))
        or np.any(error_values <= 0.0)
        or np.any(h_values <= 0.0)
    ):
        raise ValueError("errors and spacings must be positive and finite")
    if not np.isfinite(minimum_order):
        raise ValueError("minimum_order must be finite")
    if not np.all(np.diff(h_values) < 0.0):
        raise ValueError("spacings must be strictly decreasing")

    decreasing = bool(np.all(np.diff(error_values) < 0.0))
    orders = np.log(error_values[:-1] / error_values[1:]) / np.log(
        h_values[:-1] / h_values[1:]
    )
    passed = decreasing and bool(np.all(orders >= minimum_order))
    return ConvergenceReport(
        passed=passed,
        errors=tuple(float(value) for value in error_values),
        spacings=tuple(float(value) for value in h_values),
        observed_orders=tuple(float(value) for value in orders),
        minimum_required_order=float(minimum_order),
        strictly_decreasing=decreasing,
    )


def rms(field: np.ndarray) -> float:
    """Root-mean-square norm with finite-value validation."""

    values = np.asarray(field, dtype=np.float64)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("rms requires a nonempty finite array")
    return float(np.sqrt(np.mean(values * values)))


def max_abs(field: np.ndarray) -> float:
    """Maximum absolute value with finite-value validation."""

    values = np.asarray(field, dtype=np.float64)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("max_abs requires a nonempty finite array")
    return float(np.max(np.abs(values)))
