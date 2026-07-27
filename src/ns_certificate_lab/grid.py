"""Grid definitions and array validation for the half-plane ``r >= 0``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class AxisymmetricGrid:
    """Uniform tensor-product grid with the cylindrical axis included.

    Parameters
    ----------
    r:
        One-dimensional radial coordinates.  The first point must be exactly
        zero and coordinates must be uniformly spaced and increasing.
    z:
        One-dimensional axial coordinates, uniformly spaced and increasing.
        For a periodic grid the right endpoint is not included.
    periodic_z:
        Whether axial derivatives wrap around.
    z_period:
        Period used by the periodic grid.  It is inferred as ``n_z * dz`` when
        omitted.
    """

    r: FloatArray
    z: FloatArray
    periodic_z: bool = True
    z_period: float | None = None

    def __post_init__(self) -> None:
        r = np.asarray(self.r, dtype=np.float64)
        z = np.asarray(self.z, dtype=np.float64)
        if r.ndim != 1 or z.ndim != 1:
            raise ValueError("r and z coordinates must be one-dimensional")
        if r.size < 4:
            raise ValueError("at least four radial points are required")
        if z.size < 5:
            raise ValueError("at least five axial points are required")
        if not np.all(np.isfinite(r)) or not np.all(np.isfinite(z)):
            raise ValueError("grid coordinates must be finite")
        if r[0] != 0.0:
            raise ValueError("the cylindrical axis r=0 must be included")
        if np.any(np.diff(r) <= 0.0) or np.any(np.diff(z) <= 0.0):
            raise ValueError("grid coordinates must be strictly increasing")
        dr_values = np.diff(r)
        dz_values = np.diff(z)
        if not np.allclose(dr_values, dr_values[0], rtol=1e-12, atol=1e-14):
            raise ValueError("radial grid must be uniform")
        if not np.allclose(dz_values, dz_values[0], rtol=1e-12, atol=1e-14):
            raise ValueError("axial grid must be uniform")

        inferred_period = float(z.size * dz_values[0])
        period = inferred_period if self.z_period is None else float(self.z_period)
        if self.periodic_z:
            if not np.isfinite(period) or period <= 0.0:
                raise ValueError("z_period must be positive and finite")
            if not np.isclose(period, inferred_period, rtol=1e-11, atol=1e-13):
                raise ValueError(
                    "periodic z grid must omit its endpoint: z_period must equal n_z * dz"
                )
        elif self.z_period is not None:
            raise ValueError("z_period is only meaningful when periodic_z=True")

        r.setflags(write=False)
        z.setflags(write=False)
        object.__setattr__(self, "r", r)
        object.__setattr__(self, "z", z)
        object.__setattr__(self, "z_period", period if self.periodic_z else None)

    @classmethod
    def uniform(
        cls,
        *,
        nr: int,
        nz: int,
        r_max: float,
        z_min: float = 0.0,
        z_max: float = 2.0 * np.pi,
        periodic_z: bool = True,
    ) -> "AxisymmetricGrid":
        """Construct a uniform grid with ``r=0`` and ``r=r_max`` included."""

        if nr < 4 or nz < 5:
            raise ValueError("nr >= 4 and nz >= 5 are required")
        if not np.isfinite(r_max) or r_max <= 0.0:
            raise ValueError("r_max must be positive and finite")
        if not (np.isfinite(z_min) and np.isfinite(z_max) and z_max > z_min):
            raise ValueError("z_min and z_max must be finite with z_max > z_min")
        r = np.linspace(0.0, r_max, nr, dtype=np.float64)
        z = np.linspace(
            z_min,
            z_max,
            nz,
            endpoint=not periodic_z,
            dtype=np.float64,
        )
        return cls(
            r=r,
            z=z,
            periodic_z=periodic_z,
            z_period=(z_max - z_min) if periodic_z else None,
        )

    @property
    def nr(self) -> int:
        return int(self.r.size)

    @property
    def nz(self) -> int:
        return int(self.z.size)

    @property
    def dr(self) -> float:
        return float(self.r[1] - self.r[0])

    @property
    def dz(self) -> float:
        return float(self.z[1] - self.z[0])

    @property
    def shape(self) -> tuple[int, int]:
        return (self.nr, self.nz)

    def mesh(self) -> tuple[FloatArray, FloatArray]:
        """Return ``R, Z`` arrays with shape ``(n_r, n_z)``."""

        return np.meshgrid(self.r, self.z, indexing="ij")

    def validate_field(
        self,
        field: npt.ArrayLike,
        *,
        name: str = "field",
        require_finite: bool = True,
    ) -> FloatArray:
        array = np.asarray(field, dtype=np.float64)
        if array.shape != self.shape:
            raise ValueError(f"{name} must have shape {self.shape}, got {array.shape}")
        if require_finite and not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must contain only finite values")
        return array
