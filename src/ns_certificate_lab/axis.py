"""Checks for regularity at the cylindrical axis.

An axisymmetric scalar that is the restriction of a smooth function on ``R^3``
is even in ``r``, so *every* odd radial derivative vanishes on the axis, not
merely the first.  Until the ``derivative_order`` argument was added below,
this module tested only ``d_r f(0,z) = 0``, the ``k = 0`` member of that
family.  That omission is not cosmetic.  The norm-transfer ladder of
:mod:`ns_certificate_lab.norm_transfer` applies ``L5 = d_rr + (3/r) d_r +
d_zz`` repeatedly, and a term ``c r^3`` in an otherwise even field survives the
first-derivative test untouched -- it has zero slope at the axis -- while
producing ``15 c r`` in ``L5 f``, which is odd.  The second application of
``L5`` then differentiates an odd function as if it were even and the ladder
silently computes the wrong number.  Concretely the ``H^4`` transfer needs
``d_r^3 psi1(0,z) = 0`` and nothing in this repository checked it.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
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


def _finite_difference_weights(order: int, nodes: np.ndarray) -> np.ndarray:
    """Fornberg weights for the ``order``-th derivative at ``0`` on ``nodes``.

    The recursion is used rather than a Vandermonde solve because a Vandermonde
    system on the nodes ``0, 1, ..., N`` is badly conditioned, and the resulting
    loss of five or so digits is then multiplied by ``h^-order``.  For a third
    derivative on a fine grid that is the difference between a residual of
    ``1e-10`` and one of ``1e-6``, which is the difference between a usable
    tolerance and none at all.
    """

    count = int(nodes.size)
    weights = np.zeros((count, order + 1), dtype=np.float64)
    weights[0, 0] = 1.0
    c1 = 1.0
    c4 = float(nodes[0])
    for i in range(1, count):
        highest = min(i, order)
        c2 = 1.0
        c5 = c4
        c4 = float(nodes[i])
        for j in range(i):
            c3 = float(nodes[i] - nodes[j])
            c2 *= c3
            if j == i - 1:
                for k in range(highest, 0, -1):
                    weights[i, k] = (
                        c1 * (k * weights[i - 1, k - 1] - c5 * weights[i - 1, k]) / c2
                    )
                weights[i, 0] = -c1 * c5 * weights[i - 1, 0] / c2
            for k in range(highest, 0, -1):
                weights[j, k] = (c4 * weights[j, k] - k * weights[j, k - 1]) / c3
            weights[j, 0] = c4 * weights[j, 0] / c3
        c1 = c2
    return weights[:, order]


def _axis_derivative_of_order(
    grid: AxisymmetricGrid,
    field: np.ndarray,
    order: int,
) -> np.ndarray:
    """One-sided radial derivative of arbitrary order on the axis row.

    The stencil is one-sided on purpose.  Imposing the even reflection would
    assume exactly the parity that is under test and would return zero for
    every input, so the audit uses only interior data.

    Order one keeps the original hand-written stencil bit for bit, so recorded
    ``max_abs_dr_at_axis`` metrics in existing artifacts remain reproducible;
    higher orders use Fornberg weights on ``order + 3`` points, which is two
    orders of accuracy beyond the derivative being estimated.
    """

    if order == 1:
        return _axis_derivative_independent(grid, field)
    span = min(grid.nr, order + 3)
    if span < order + 2:
        raise ValueError(
            f"a radial derivative of order {order} at the axis needs at least "
            f"{order + 2} radial points, the grid has {grid.nr}"
        )
    nodes = np.arange(span, dtype=np.float64)
    weights = _finite_difference_weights(order, nodes) / grid.dr**order
    return np.tensordot(weights, field[:span], axes=(0, 0))


def _validated_derivative_order(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError("derivative_order must be an integer")
    order = int(value)
    if order < 0:
        raise ValueError("derivative_order must be nonnegative")
    return order


def required_axis_derivative_order(sobolev_order: int) -> int:
    """Highest odd radial derivative the ``H^n`` transfer needs to vanish.

    In the ladder of :mod:`ns_certificate_lab.norm_transfer` the stream
    generator is carried to rung ``n+1``, which applies ``L5`` exactly
    ``floor((n+1)/2)`` times, and the ``k``-th application produces an axis-even
    output only if the ``(2k-1)``-th radial derivative of its input vanishes.
    The required odd orders are therefore ``1, 3, ..., 2*floor((n+1)/2)-1``,
    that is every odd order up to the largest odd integer not exceeding ``n``.
    For ``n = 4`` that is ``1`` and ``3``; the ``3`` is the condition nothing in
    this repository used to test.  The value is never less than one, because
    even the ``L^2`` transfer differentiates ``psi1`` once.
    """

    order = _validated_derivative_order(sobolev_order)
    return max(1, 2 * math.floor((order + 1) / 2) - 1)


def check_axis_regularity(
    grid: AxisymmetricGrid,
    *,
    even_fields: Mapping[str, npt.ArrayLike] | None = None,
    odd_fields: Mapping[str, npt.ArrayLike] | None = None,
    tolerance: float = 1e-8,
    derivative_order: int = 1,
) -> AxisRegularityReport:
    """Check necessary parity conditions at ``r=0``.

    Transformed scalars ``u1``, ``omega1`` and ``psi1`` are axis-even, as is
    ``u^z``; hence their radial derivative must vanish at the axis.  Physical
    ``u^r`` and ``u^theta = r u1`` are axis-odd and must vanish there.

    ``derivative_order`` is the highest radial derivative to test.  Even fields
    are tested at the odd orders ``1, 3, ...`` not exceeding it and odd fields
    at the even orders ``0, 2, ...`` not exceeding it, those being the
    derivatives that parity forces to vanish.  The default of one reproduces
    the historical behaviour exactly, including every metric key; the extra
    orders appear as ``<name>.max_abs_dr<p>_at_axis`` only when they are asked
    for.  Use :func:`required_axis_derivative_order` to pick the order that a
    given Sobolev transfer actually needs, and see the module docstring for why
    testing only order one is not enough.

    A derivative of order ``p`` carries a factor ``h^-p``, so a single fixed
    tolerance does not mean the same thing at every order: a residual of
    ``1e-8`` in a third derivative is a far smaller defect in the field than a
    residual of ``1e-8`` in a first derivative.  The metrics are reported raw,
    in field units per length to the ``p``, and it is the caller's business to
    know which one is being read.

    These are necessary discrete checks, not a proof that a finite array has a
    smooth even/odd extension.
    """

    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("tolerance must be nonnegative and finite")
    highest = _validated_derivative_order(derivative_order)
    metrics: dict[str, float] = {}
    violations: list[str] = []

    def record(key: str, residual: float) -> None:
        metrics[key] = residual
        if residual > tolerance:
            violations.append(key)

    for name, values in (even_fields or {}).items():
        field = grid.validate_field(values, name=name)
        for order in range(1, highest + 1, 2):
            residual = float(
                np.max(np.abs(_axis_derivative_of_order(grid, field, order)))
            )
            suffix = "dr" if order == 1 else f"dr{order}"
            record(f"{name}.max_abs_{suffix}_at_axis", residual)

    for name, values in (odd_fields or {}).items():
        field = grid.validate_field(values, name=name)
        for order in range(0, highest + 1, 2):
            if order == 0:
                residual = float(np.max(np.abs(field[0])))
                key = f"{name}.max_abs_value_at_axis"
            else:
                residual = float(
                    np.max(np.abs(_axis_derivative_of_order(grid, field, order)))
                )
                key = f"{name}.max_abs_dr{order}_at_axis"
            record(key, residual)

    return AxisRegularityReport(
        passed=not violations,
        tolerance=float(tolerance),
        metrics=metrics,
        violations=tuple(violations),
    )


def check_ladder_axis_regularity(
    grid: AxisymmetricGrid,
    *,
    u1: npt.ArrayLike,
    psi1: npt.ArrayLike,
    sobolev_order: int,
    tolerance: float = 1e-8,
) -> AxisRegularityReport:
    """Check the axis conditions the ``H^n`` norm transfer actually requires.

    This is :func:`check_axis_regularity` on the two generators at the order
    returned by :func:`required_axis_derivative_order`.  It exists so that a
    caller who has chosen a Sobolev index cannot accidentally audit a weaker
    condition than the one the chosen index depends on.
    """

    return check_axis_regularity(
        grid,
        even_fields={"u1": u1, "psi1": psi1},
        tolerance=tolerance,
        derivative_order=required_axis_derivative_order(sobolev_order),
    )
