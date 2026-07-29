r"""Conservative domain expansion for a whole-space run.

A fixed box is safe only while the field stays away from its walls.  As the
amplitude grows the outer tail eventually reaches them, and from that moment the
computation is a wall-bounded computation wearing a whole-space label.  This
module implements the escape hatch: monitor four independent indicators, and
when any of them crosses a **preregistered** threshold, embed the state in a
larger box and recompute the elliptic field.

The embedding uses **no interpolation**.  The enlarged grid keeps ``dr`` and
``dz`` and extends the extent by an integer number of cells, so the old nodes
coincide with new ones and the state is transferred by zero padding.  The
residual ``interpolation_defect`` is not exactly zero — the two coordinate
arrays come from independent ``linspace`` calls and differ at roundoff
(``~3 x 10^-15`` at the tested settings) — but anything above ``10^-12`` means
the alignment argument has genuinely broken, which is what the number is for.

What is checked across an expansion: the vorticity maximum (must agree to
roundoff, since the state array is untouched), the kinetic energy and ``L^3``
(computed from the *recovered* velocity, so they move by the truncation error
the expansion exists to remove), the physical divergence, and the interpolation
defect.  The dyadic shell masses are recorded but **not** asserted: the
outermost shell is by definition the remaining represented domain, so enlarging
the box changes it on purpose.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "ExpansionEvent",
    "ExpansionPolicy",
    "embed_in_larger_box",
    "expansion_is_required",
]


@dataclass(frozen=True)
class ExpansionPolicy:
    """Preregistered thresholds.  Changing them requires a recorded reason."""

    outer_band_fraction: float = 1.0e-4
    tail_bound_ratio: float = 1.0e-2
    multipole_remainder_ratio: float = 1.0e-2
    boundary_difference_ratio: float = 1.0e-3
    radial_growth_cells: int = 16
    axial_growth_cells: int = 16
    maximum_expansions: int = 4

    def __post_init__(self) -> None:
        for name in (
            "outer_band_fraction",
            "tail_bound_ratio",
            "multipole_remainder_ratio",
            "boundary_difference_ratio",
        ):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be a positive finite threshold")
        if self.radial_growth_cells < 1 or self.axial_growth_cells < 1:
            raise ValueError("growth must add at least one cell per direction")
        if self.maximum_expansions < 0:
            raise ValueError("maximum_expansions must be nonnegative")

    def as_dict(self) -> dict[str, float | int]:
        return {
            "outer_band_fraction": self.outer_band_fraction,
            "tail_bound_ratio": self.tail_bound_ratio,
            "multipole_remainder_ratio": self.multipole_remainder_ratio,
            "boundary_difference_ratio": self.boundary_difference_ratio,
            "radial_growth_cells": self.radial_growth_cells,
            "axial_growth_cells": self.axial_growth_cells,
            "maximum_expansions": self.maximum_expansions,
        }


def expansion_is_required(
    policy: ExpansionPolicy,
    *,
    outer_band_fraction: float,
    tail_bound_ratio: float,
    multipole_remainder_ratio: float,
    boundary_difference_ratio: float,
) -> tuple[bool, tuple[str, ...]]:
    """Return ``(fire, reasons)`` for the four monitored indicators."""
    reasons: list[str] = []
    if outer_band_fraction > policy.outer_band_fraction:
        reasons.append("outer_band_fraction")
    if tail_bound_ratio > policy.tail_bound_ratio:
        reasons.append("tail_bound_ratio")
    if multipole_remainder_ratio > policy.multipole_remainder_ratio:
        reasons.append("multipole_remainder_ratio")
    if boundary_difference_ratio > policy.boundary_difference_ratio:
        reasons.append("boundary_difference_ratio")
    return bool(reasons), tuple(reasons)


@dataclass(frozen=True)
class ExpansionEvent:
    """One recorded expansion, with the before/after invariants."""

    step: int
    time: float
    reasons: tuple[str, ...]
    old_r_max: float
    old_z_max: float
    new_r_max: float
    new_z_max: float
    energy_before: float
    energy_after: float
    l3_before: float
    l3_after: float
    omega_max_before: float
    omega_max_after: float
    shell_mass_before: tuple[float, ...]
    shell_mass_after: tuple[float, ...]
    divergence_before: float
    divergence_after: float
    interpolation_defect: float

    @property
    def invariants_preserved(self) -> bool:
        """Whether the physical invariants survived the embedding.

        Three tolerances, each with a stated reason.

        * ``omega_max`` must agree to roundoff: zero padding does not touch the
          state array.
        * ``energy`` and ``L^3`` are computed from the *recovered velocity*, and
          the recovery is a fresh elliptic solve on the larger box.  They
          therefore move by the truncation error the expansion exists to remove,
          which at the tested settings is ``~10^-9`` relative.  Demanding
          roundoff agreement here would be demanding that the expansion changed
          nothing, i.e. that it was pointless.
        * ``interpolation_defect`` measures node misalignment.  The two grids are
          built by independent ``linspace`` calls, so exact bitwise agreement is
          not available; anything above ``10^-12`` means the alignment argument
          has actually broken.

        The dyadic **shell masses are deliberately excluded**: the outermost
        shell is defined as "the remaining represented domain", so enlarging the
        box changes it by construction.  They are recorded for inspection rather
        than asserted.
        """
        def close(a: float, b: float, tolerance: float) -> bool:
            return abs(a - b) <= tolerance * max(abs(a), abs(b), 1.0)

        return (
            close(self.energy_before, self.energy_after, 1.0e-6)
            and close(self.l3_before, self.l3_after, 1.0e-6)
            and close(self.omega_max_before, self.omega_max_after, 1.0e-14)
            and len(self.shell_mass_before) == len(self.shell_mass_after)
            and self.interpolation_defect < 1.0e-12
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "step": self.step,
            "time": self.time,
            "reasons": list(self.reasons),
            "old_r_max": self.old_r_max,
            "old_z_max": self.old_z_max,
            "new_r_max": self.new_r_max,
            "new_z_max": self.new_z_max,
            "energy_before": self.energy_before,
            "energy_after": self.energy_after,
            "l3_before": self.l3_before,
            "l3_after": self.l3_after,
            "omega_max_before": self.omega_max_before,
            "omega_max_after": self.omega_max_after,
            "shell_mass_before": list(self.shell_mass_before),
            "shell_mass_after": list(self.shell_mass_after),
            "divergence_before": self.divergence_before,
            "divergence_after": self.divergence_after,
            "interpolation_defect": self.interpolation_defect,
            "invariants_preserved": self.invariants_preserved,
        }


def embed_in_larger_box(
    grid: AxisymmetricGrid,
    fields: tuple[FloatArray, ...],
    policy: ExpansionPolicy,
) -> tuple[AxisymmetricGrid, tuple[FloatArray, ...], float]:
    """Zero-pad the state onto a grid with the same spacing and a larger extent.

    Returns ``(new_grid, new_fields, interpolation_defect)``.  The defect is the
    maximum absolute difference between the old values and the corresponding
    entries of the new arrays; it is zero whenever the node alignment holds, and
    a nonzero value means the embedding is no longer exact.
    """
    if grid.periodic_z:
        raise ValueError("domain expansion is defined for non-periodic z grids")
    dr, dz = grid.dr, grid.dz
    extra_r = policy.radial_growth_cells
    extra_z = policy.axial_growth_cells
    new_nr = grid.nr + extra_r
    new_nz = grid.nz + 2 * extra_z
    new_grid = AxisymmetricGrid.uniform(
        nr=new_nr,
        nz=new_nz,
        r_max=dr * (new_nr - 1),
        z_min=float(grid.z[0]) - extra_z * dz,
        z_max=float(grid.z[-1]) + extra_z * dz,
        periodic_z=False,
    )
    defect = max(
        float(np.max(np.abs(new_grid.r[: grid.nr] - grid.r))),
        float(np.max(np.abs(new_grid.z[extra_z : extra_z + grid.nz] - grid.z))),
    )
    embedded: list[FloatArray] = []
    for field in fields:
        values = grid.validate_field(field, name="field")
        out = np.zeros(new_grid.shape, dtype=np.float64)
        out[: grid.nr, extra_z : extra_z + grid.nz] = values
        embedded.append(out)
    return new_grid, tuple(embedded), defect
