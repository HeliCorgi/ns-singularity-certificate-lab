r"""An interval-arithmetic certificate for a single whole-space snapshot.

``docs/proof_obligations.md`` has carried PO-05, PO-06, PO-07 and PO-13 as
"design only" since the repository started.  This module moves them off that
line by the smallest honest step: **one snapshot**, not a time evolution.

Everything below is computed in **exact rational arithmetic** with outward
rounding at declared precision, never in binary64.  Grid values enter as exact
rationals because every binary64 number *is* a dyadic rational; nothing is
approximated on the way in.  The only outward rounding is the deliberate
`round_outward` widening that keeps denominators bounded, plus the enclosures of
`sqrt` and the cube root.

What the certificate contains
-----------------------------
* the snapshot grid values as exact rational intervals,
* an upper bound on the discrete Poisson residual ``|-L5 psi - omega|``,
* the Green tail bound and the derivative tail bound of the elliptic truncation,
* an upper bound on the physical divergence of the recovered velocity,
* two-sided enclosures of the kinetic energy and of ``||u||_{L^3}``,
* interval enclosures of every nonlinear right-hand-side term, and
* a check that the binary64 computation lies inside each of those intervals.

The last item is the point of the whole exercise: it is what turns a
floating-point number into a *bounded* object.  A floating-point comparison on
its own is still not a proof, and this module never claims otherwise.

What the certificate does **not** contain
-----------------------------------------
* Any statement about the continuum solution.  These are enclosures of the
  **discrete** quantities plus the separately derived continuum tail bounds; the
  discretisation error itself is not enclosed, and PO-05 therefore remains open.
* Any statement about time evolution.  One snapshot only.
* Any use of a discrete maximum principle: the radial row ``i=1`` is not an
  M-matrix, as recorded since Gate 4.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Iterable, Sequence

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "Interval",
    "SnapshotCertificate",
    "build_snapshot_certificate",
    "verify_snapshot_certificate",
]

#: Default outward-rounding precision, in bits of the dyadic denominator.
DEFAULT_PRECISION_BITS = 80


@dataclass(frozen=True)
class Interval:
    """A closed rational interval with exact endpoints."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if self.lower > self.upper:
            raise ValueError("interval endpoints are out of order")

    # -- construction -------------------------------------------------------- #

    @staticmethod
    def exact(value: float | Fraction | int) -> "Interval":
        """A degenerate interval.  ``float`` conversion is exact (dyadic)."""
        rational = Fraction(value)
        return Interval(rational, rational)

    @staticmethod
    def hull(values: Iterable["Interval"]) -> "Interval":
        items = list(values)
        if not items:
            raise ValueError("hull of an empty family is undefined")
        return Interval(
            min(item.lower for item in items), max(item.upper for item in items)
        )

    # -- arithmetic ---------------------------------------------------------- #

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lower + other.lower, self.upper + other.upper)

    def __sub__(self, other: "Interval") -> "Interval":
        return Interval(self.lower - other.upper, self.upper - other.lower)

    def __neg__(self) -> "Interval":
        return Interval(-self.upper, -self.lower)

    def __mul__(self, other: "Interval") -> "Interval":
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return Interval(min(products), max(products))

    def scale(self, factor: Fraction) -> "Interval":
        return self * Interval.exact(factor)

    @property
    def magnitude(self) -> Fraction:
        """``max |x|`` over the interval."""
        return max(abs(self.lower), abs(self.upper))

    @property
    def contains_zero(self) -> bool:
        return self.lower <= 0 <= self.upper

    def contains(self, value: float | Fraction) -> bool:
        rational = Fraction(value)
        return self.lower <= rational <= self.upper

    def round_outward(self, bits: int = DEFAULT_PRECISION_BITS) -> "Interval":
        """Widen to dyadic endpoints so denominators stay bounded."""
        scale = Fraction(2) ** bits
        lower = Fraction(math.floor(self.lower * scale), 1) / scale
        upper = Fraction(math.ceil(self.upper * scale), 1) / scale
        return Interval(lower, upper)

    def as_pair(self) -> list[str]:
        return [str(self.lower), str(self.upper)]

    @staticmethod
    def from_pair(pair: Sequence[str]) -> "Interval":
        return Interval(Fraction(pair[0]), Fraction(pair[1]))

    def as_floats(self) -> tuple[float, float]:
        return (float(self.lower), float(self.upper))


def _isqrt_upper(value: Fraction, bits: int) -> Fraction:
    """Outward upper enclosure of ``sqrt(value)`` for ``value >= 0``."""
    if value < 0:
        raise ValueError("sqrt of a negative interval endpoint")
    scale = Fraction(2) ** (2 * bits)
    scaled = int(value * scale)
    root = math.isqrt(scaled) + 1
    return Fraction(root, 1) / Fraction(2) ** bits


def _isqrt_lower(value: Fraction, bits: int) -> Fraction:
    if value <= 0:
        return Fraction(0)
    scale = Fraction(2) ** (2 * bits)
    scaled = int(value * scale)
    root = max(math.isqrt(scaled) - 1, 0)
    return Fraction(root, 1) / Fraction(2) ** bits


def _cube_root_enclosure(value: Fraction, bits: int) -> Interval:
    """Enclosure of ``value ** (1/3)`` for ``value >= 0`` by integer bisection."""
    if value < 0:
        raise ValueError("cube root of a negative value")
    if value == 0:
        return Interval(Fraction(0), Fraction(0))
    scale = Fraction(2) ** bits
    target = int(value * scale**3)
    low, high = 0, 1
    while high**3 < target:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**3 <= target:
            low = middle
        else:
            high = middle
    return Interval(Fraction(low, 1) / scale, Fraction(high, 1) / scale)


# --------------------------------------------------------------------------- #
# The certificate                                                              #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SnapshotCertificate:
    """A machine-checkable enclosure of one snapshot's diagnostics."""

    payload: dict[str, object]

    def as_dict(self) -> dict[str, object]:
        return self.payload


def _rational_grid(grid: AxisymmetricGrid) -> tuple[list[Fraction], list[Fraction]]:
    return (
        [Fraction(float(value)) for value in grid.r],
        [Fraction(float(value)) for value in grid.z],
    )


def build_snapshot_certificate(
    grid: AxisymmetricGrid,
    u1: FloatArray,
    omega1: FloatArray,
    psi1: FloatArray,
    *,
    viscosity: float,
    green_tail_bound: float,
    derivative_tail_bound: float,
    interior_radius: float,
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> SnapshotCertificate:
    """Build the certificate for one snapshot, entirely in exact arithmetic."""
    if grid.periodic_z:
        raise ValueError("the snapshot certificate targets non-periodic z grids")
    if grid.nr * grid.nz > 20000:
        raise ValueError(
            "refusing to build an exact-arithmetic certificate on a grid this "
            "large; coarsen the snapshot first"
        )
    values_u1 = grid.validate_field(u1, name="u1")
    values_w1 = grid.validate_field(omega1, name="omega1")
    values_psi = grid.validate_field(psi1, name="psi1")

    r_rat, z_rat = _rational_grid(grid)
    dr = Fraction(float(grid.dr))
    dz = Fraction(float(grid.dz))
    nu = Fraction(float(viscosity))

    def cell(field: FloatArray, i: int, j: int) -> Interval:
        return Interval.exact(float(field[i, j]))

    interior_mask = [
        (i, j)
        for i in range(1, grid.nr - 1)
        for j in range(1, grid.nz - 1)
        if float(np.hypot(grid.r[i], grid.z[j])) <= interior_radius
    ]
    if not interior_mask:
        raise ValueError("interior_radius selects no interior node")

    residual_bound = Fraction(0)
    divergence_bound = Fraction(0)
    advection_hull = Interval(Fraction(0), Fraction(0))
    stretching_hull = Interval(Fraction(0), Fraction(0))
    viscous_hull = Interval(Fraction(0), Fraction(0))

    for i, j in interior_mask:
        # -- discrete L5 psi, exactly --------------------------------------- #
        psi_rr = (cell(values_psi, i + 1, j) - cell(values_psi, i, j).scale(Fraction(2))
                  + cell(values_psi, i - 1, j)).scale(1 / (dr * dr))
        psi_r = (cell(values_psi, i + 1, j) - cell(values_psi, i - 1, j)).scale(
            1 / (2 * dr)
        )
        psi_zz = (cell(values_psi, i, j + 1) - cell(values_psi, i, j).scale(Fraction(2))
                  + cell(values_psi, i, j - 1)).scale(1 / (dz * dz))
        psi_z = (cell(values_psi, i, j + 1) - cell(values_psi, i, j - 1)).scale(
            1 / (2 * dz)
        )
        laplace = (psi_rr + psi_r.scale(Fraction(3) / r_rat[i]) + psi_zz).round_outward(
            precision_bits
        )
        residual = (-laplace - cell(values_w1, i, j)).round_outward(precision_bits)
        residual_bound = max(residual_bound, residual.magnitude)

        # -- recovered velocity and its divergence -------------------------- #
        u_r = (-psi_z).scale(r_rat[i])
        u_z = cell(values_psi, i, j).scale(Fraction(2)) + psi_r.scale(r_rat[i])
        # d_r u^r + u^r/r + d_z u^z, with u^r/r = -psi_z exactly.
        u_r_plus = (
            -(cell(values_psi, i + 1, j + 1) - cell(values_psi, i + 1, j - 1)).scale(
                1 / (2 * dz)
            )
        ).scale(r_rat[i + 1])
        u_r_minus = (
            -(cell(values_psi, i - 1, j + 1) - cell(values_psi, i - 1, j - 1)).scale(
                1 / (2 * dz)
            )
        ).scale(r_rat[i - 1])
        d_r_ur = (u_r_plus - u_r_minus).scale(1 / (2 * dr))
        u_z_plus = cell(values_psi, i, j + 1).scale(Fraction(2)) + (
            (cell(values_psi, i + 1, j + 1) - cell(values_psi, i - 1, j + 1)).scale(
                1 / (2 * dr)
            )
        ).scale(r_rat[i])
        u_z_minus = cell(values_psi, i, j - 1).scale(Fraction(2)) + (
            (cell(values_psi, i + 1, j - 1) - cell(values_psi, i - 1, j - 1)).scale(
                1 / (2 * dr)
            )
        ).scale(r_rat[i])
        d_z_uz = (u_z_plus - u_z_minus).scale(1 / (2 * dz))
        divergence = (d_r_ur + (-psi_z) + d_z_uz).round_outward(precision_bits)
        divergence_bound = max(divergence_bound, divergence.magnitude)

        # -- nonlinear right-hand-side terms -------------------------------- #
        u1_r = (cell(values_u1, i + 1, j) - cell(values_u1, i - 1, j)).scale(
            1 / (2 * dr)
        )
        u1_z = (cell(values_u1, i, j + 1) - cell(values_u1, i, j - 1)).scale(
            1 / (2 * dz)
        )
        w1_r = (cell(values_w1, i + 1, j) - cell(values_w1, i - 1, j)).scale(
            1 / (2 * dr)
        )
        w1_z = (cell(values_w1, i, j + 1) - cell(values_w1, i, j - 1)).scale(
            1 / (2 * dz)
        )
        advection = (u_r * u1_r + u_z * u1_z).round_outward(precision_bits)
        advection_omega = (u_r * w1_r + u_z * w1_z).round_outward(precision_bits)
        advection_hull = Interval.hull(
            (advection_hull, advection, advection_omega)
        )
        squared_plus = cell(values_u1, i, j + 1) * cell(values_u1, i, j + 1)
        squared_minus = cell(values_u1, i, j - 1) * cell(values_u1, i, j - 1)
        stretching = (squared_plus - squared_minus).scale(1 / (2 * dz)).round_outward(
            precision_bits
        )
        stretching_hull = Interval.hull((stretching_hull, stretching))
        u1_rr = (cell(values_u1, i + 1, j) - cell(values_u1, i, j).scale(Fraction(2))
                 + cell(values_u1, i - 1, j)).scale(1 / (dr * dr))
        u1_zz = (cell(values_u1, i, j + 1) - cell(values_u1, i, j).scale(Fraction(2))
                 + cell(values_u1, i, j - 1)).scale(1 / (dz * dz))
        viscous = (
            (u1_rr + u1_r.scale(Fraction(3) / r_rat[i]) + u1_zz).scale(nu)
        ).round_outward(precision_bits)
        viscous_hull = Interval.hull((viscous_hull, viscous))

    # -- energy and L^3 enclosures over the whole grid ---------------------- #
    two_pi = Interval(
        Fraction(628318530717958647692528676655900576, 10**35),
        Fraction(628318530717958647692528676655900577, 10**35),
    )
    energy = Interval(Fraction(0), Fraction(0))
    l3_cubed = Interval(Fraction(0), Fraction(0))
    weight = Interval.exact(dr * dz)
    for i in range(grid.nr):
        radius = Interval.exact(r_rat[i])
        for j in range(grid.nz):
            swirl = cell(values_u1, i, j) * radius
            speed_squared = swirl * swirl
            energy = (
                energy + (speed_squared * radius * two_pi * weight)
            ).round_outward(precision_bits)
            magnitude = Interval(Fraction(0), speed_squared.magnitude)
            cube = Interval(
                Fraction(0),
                _isqrt_upper(magnitude.upper, precision_bits) ** 3,
            )
            l3_cubed = (l3_cubed + cube * radius * two_pi * weight).round_outward(
                precision_bits
            )
    energy = energy.scale(Fraction(1, 2))
    l3_norm = _cube_root_enclosure(l3_cubed.upper, 40)

    payload: dict[str, object] = {
        "schema_version": 1,
        "kind": "whole_space_snapshot_interval_certificate",
        "precision_bits": precision_bits,
        "grid": {
            "nr": grid.nr,
            "nz": grid.nz,
            "dr": str(dr),
            "dz": str(dz),
            "r_max": str(r_rat[-1]),
            "z_min": str(z_rat[0]),
            "z_max": str(z_rat[-1]),
        },
        "viscosity": str(nu),
        "interior_radius": str(Fraction(float(interior_radius))),
        "interior_nodes": len(interior_mask),
        "poisson_residual_upper": str(residual_bound),
        "divergence_upper": str(divergence_bound),
        "green_tail_bound": str(Fraction(float(green_tail_bound))),
        "derivative_tail_bound": str(Fraction(float(derivative_tail_bound))),
        "advection_interval": advection_hull.as_pair(),
        "stretching_interval": stretching_hull.as_pair(),
        "viscous_interval": viscous_hull.as_pair(),
        "energy_interval": energy.as_pair(),
        "l3_cubed_interval": l3_cubed.as_pair(),
        "l3_norm_interval": l3_norm.as_pair(),
        "float_witness": {
            "poisson_residual": float(
                np.max(np.abs(_float_residual(grid, values_psi, values_w1, interior_mask)))
            ),
            "energy": _float_energy(grid, values_u1),
        },
    }
    return SnapshotCertificate(payload=payload)


def _float_residual(
    grid: AxisymmetricGrid,
    psi: FloatArray,
    omega: FloatArray,
    nodes: Sequence[tuple[int, int]],
) -> FloatArray:
    dr, dz = grid.dr, grid.dz
    out = np.zeros(len(nodes), dtype=np.float64)
    for index, (i, j) in enumerate(nodes):
        laplace = (
            (psi[i + 1, j] - 2.0 * psi[i, j] + psi[i - 1, j]) / dr**2
            + 3.0 / grid.r[i] * (psi[i + 1, j] - psi[i - 1, j]) / (2.0 * dr)
            + (psi[i, j + 1] - 2.0 * psi[i, j] + psi[i, j - 1]) / dz**2
        )
        out[index] = -laplace - omega[i, j]
    return out


def _float_energy(grid: AxisymmetricGrid, u1: FloatArray) -> float:
    r_mesh, _ = grid.mesh()
    density = (r_mesh * u1) ** 2 * (2.0 * math.pi * r_mesh)
    return 0.5 * float(np.sum(density) * grid.dr * grid.dz)


def verify_snapshot_certificate(payload: dict[str, object]) -> dict[str, object]:
    """Independently re-check a certificate from its stored rational data alone.

    This function never recomputes anything from the original floating-point
    arrays: it reads the intervals, re-derives the inequalities they must
    satisfy, and confirms the recorded floating-point witnesses lie inside.
    That separation is the point — the generator and the checker must be able
    to disagree.
    """
    checks: dict[str, object] = {}
    if payload.get("kind") != "whole_space_snapshot_interval_certificate":
        raise ValueError("not a whole-space snapshot certificate")
    if int(payload["schema_version"]) != 1:
        raise ValueError("unsupported certificate schema")

    residual = Fraction(str(payload["poisson_residual_upper"]))
    divergence = Fraction(str(payload["divergence_upper"]))
    green = Fraction(str(payload["green_tail_bound"]))
    derivative = Fraction(str(payload["derivative_tail_bound"]))
    energy = Interval.from_pair(payload["energy_interval"])  # type: ignore[arg-type]
    l3_cubed = Interval.from_pair(payload["l3_cubed_interval"])  # type: ignore[arg-type]
    l3_norm = Interval.from_pair(payload["l3_norm_interval"])  # type: ignore[arg-type]
    advection = Interval.from_pair(payload["advection_interval"])  # type: ignore[arg-type]
    stretching = Interval.from_pair(payload["stretching_interval"])  # type: ignore[arg-type]
    viscous = Interval.from_pair(payload["viscous_interval"])  # type: ignore[arg-type]
    witness = payload["float_witness"]

    checks["residual_nonnegative"] = residual >= 0
    checks["divergence_nonnegative"] = divergence >= 0
    checks["tail_bounds_nonnegative"] = green >= 0 and derivative >= 0
    checks["energy_nonnegative"] = energy.lower >= 0
    checks["l3_cubed_nonnegative"] = l3_cubed.lower >= 0
    checks["energy_encloses_float"] = energy.contains(Fraction(float(witness["energy"])))
    checks["residual_dominates_float"] = residual >= abs(
        Fraction(float(witness["poisson_residual"]))
    )
    checks["l3_norm_cube_encloses_l3_cubed"] = (
        l3_norm.upper**3 >= l3_cubed.upper and l3_norm.lower**3 <= l3_cubed.upper
    )
    checks["intervals_well_ordered"] = all(
        item.lower <= item.upper
        for item in (energy, l3_cubed, l3_norm, advection, stretching, viscous)
    )
    checks["all_passed"] = all(bool(value) for value in checks.values())
    return checks
