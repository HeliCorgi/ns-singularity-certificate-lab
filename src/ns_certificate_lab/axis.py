"""Checks for regularity at the cylindrical axis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid


@dataclass(frozen=True)
class AxisRegularityReport:
    """Quantitative axis checks for half-plane data."""

    passed: bool
    tolerance: float
    metrics: dict[str, float]
    violations: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "passed": self.passed,
            "tolerance": self.tolerance,
            "metrics": dict(self.metrics),
            "violations": list(self.violations),
        }


def _axis_derivative_independent(grid: AxisymmetricGrid, field: np.ndarray) -> np.ndarray:
    """One-sided derivative used only for auditing, not PDE evaluation."""

    h = grid.dr
    if grid.nr >= 5:
        return (
            -25.0 * field[0]
            + 48.0 * field[1]
            - 36.0 * field[2]
            + 16.0 * field[3]
            - 3.0 * field[4]
        ) / (12.0 * h)
    return (-11.0 * field[0] + 18.0 * field[1] - 9.0 * field[2] + 2.0 * field[3]) / (
        6.0 * h
    )


def check_axis_regularity(
    grid: AxisymmetricGrid,
    *,
    even_fields: Mapping[str, npt.ArrayLike] | None = None,
    odd_fields: Mapping[str, npt.ArrayLike] | None = None,
    tolerance: float = 1e-8,
) -> AxisRegularityReport:
    """Check necessary parity conditions at ``r=0``.

    Transformed scalars ``u1``, ``omega1`` and ``psi1`` are axis-even, as is
    ``u^z``; hence their radial derivative must vanish at the axis.  Physical
    ``u^r`` and ``u^theta = r u1`` are axis-odd and must vanish there.

    These are necessary discrete checks, not a proof that a finite array has a
    smooth even/odd extension.
    """

    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("tolerance must be nonnegative and finite")
    metrics: dict[str, float] = {}
    violations: list[str] = []

    for name, values in (even_fields or {}).items():
        field = grid.validate_field(values, name=name)
        residual = float(np.max(np.abs(_axis_derivative_independent(grid, field))))
        key = f"{name}.max_abs_dr_at_axis"
        metrics[key] = residual
        if residual > tolerance:
            violations.append(key)

    for name, values in (odd_fields or {}).items():
        field = grid.validate_field(values, name=name)
        residual = float(np.max(np.abs(field[0])))
        key = f"{name}.max_abs_value_at_axis"
        metrics[key] = residual
        if residual > tolerance:
            violations.append(key)

    return AxisRegularityReport(
        passed=not violations,
        tolerance=float(tolerance),
        metrics=metrics,
        violations=tuple(violations),
    )
