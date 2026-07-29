r"""An explicit whole-space initial-data family for Clay-admissible candidates.

One family is fixed here and used for every whole-space candidate run, so that
"which initial data" is never a free parameter of a reported result.

.. math::

   u_1(0,r,z) = A\,
       \chi\!\left(\frac{r^2}{R_0^2}\right)
       \chi\!\left(\frac{z^2}{Z_0^2}\right)
       \frac{z/Z_0}{1 + c\,(z/Z_0)^2},
   \qquad \omega_1(0,r,z) = 0 ,

with the standard ``C^\infty`` bump

.. math::

   \chi(s) = \begin{cases}\exp\!\bigl(-\tfrac{1}{1-s}\bigr), & 0 \le s < 1,\\
             0, & s \ge 1.\end{cases}

**The radial factor is a function of ``r^2``.**  That is the whole point: the
physical swirl is ``u^\theta = r u_1``, whose Cartesian form is
``u_1\,(-y, x, 0)``, and this is ``C^\infty`` at the origin exactly when ``u_1``
is a smooth function of ``(r^2, z)``.  Writing ``\chi(r/R_0)`` instead would give
a field that is only Lipschitz across the axis.

Properties proved on paper and checked in ``tests/test_initial_data.py``:

* **Cartesian smoothness** — ``u = u_1(r^2,z)\,(-y,x,0)`` with ``u_1`` a
  composition of ``C^\infty`` functions of ``x^2+y^2`` and ``z``.
* **Divergence free, exactly** — a pure axisymmetric swirl has
  ``\nabla\cdot u = r^{-1}\partial_\theta u^\theta = 0`` identically, not merely
  to discretisation order.
* **Compact support** — contained in ``\{r \le R_0\} \cap \{|z| \le Z_0\}``.
* **Finite energy and finite ``L^3``** — both integrals are over a compact set
  with a bounded integrand.
* **Axis regularity** — ``u^\theta = r u_1`` vanishes on the axis together with
  the whole Cartesian vector, and ``\partial_r u_1|_{r=0} = 0`` because ``u_1``
  is even in ``r``.

The viscosity is a **fixed positive constant** for every production run.  Neither
a time-dependent viscosity nor the two-stage protocol of the Hou setup is used
for a Clay candidate: those are reproduction devices, not admissible data.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "SwirlFamily",
    "smooth_bump",
    "smooth_bump_derivative",
]


def smooth_bump(s: npt.ArrayLike) -> FloatArray:
    r"""``\chi(s) = \exp(-1/(1-s))`` for ``0 \le s < 1``, else ``0``.

    Evaluated on the *squared* argument by every caller, so the composition is a
    ``C^\infty`` function of the Cartesian coordinates.  All derivatives vanish
    as ``s \to 1^-``, so the extension by zero is ``C^\infty``.
    """
    values = np.asarray(s, dtype=np.float64)
    out = np.zeros(np.shape(values), dtype=np.float64)
    inside = (values < 1.0) & (values >= 0.0)
    if np.any(inside):
        out[inside] = np.exp(-1.0 / (1.0 - values[inside]))
    return out


def smooth_bump_derivative(s: npt.ArrayLike) -> FloatArray:
    r"""``\chi'(s) = -\chi(s)/(1-s)^2``."""
    values = np.asarray(s, dtype=np.float64)
    out = np.zeros(np.shape(values), dtype=np.float64)
    inside = (values < 1.0) & (values >= 0.0)
    if np.any(inside):
        gap = 1.0 - values[inside]
        out[inside] = -np.exp(-1.0 / gap) / gap**2
    return out


@dataclass(frozen=True)
class SwirlFamily:
    """The fixed candidate family, parameterised by amplitude and shape."""

    amplitude: float
    radial_support: float
    axial_support: float
    axial_concentration: float = 0.0

    def __post_init__(self) -> None:
        for name in ("amplitude", "radial_support", "axial_support", "axial_concentration"):
            value = getattr(self, name)
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if self.radial_support <= 0.0 or self.axial_support <= 0.0:
            raise ValueError("radial_support and axial_support must be positive")
        if self.axial_concentration < 0.0:
            raise ValueError("axial_concentration must be nonnegative")

    @property
    def aspect_ratio(self) -> float:
        return self.radial_support / self.axial_support

    @property
    def support_radius(self) -> float:
        """Five-dimensional radius of a ball containing the support."""
        return math.hypot(self.radial_support, self.axial_support)

    def label(self) -> str:
        return (
            f"A{self.amplitude:g}_R{self.radial_support:g}"
            f"_Z{self.axial_support:g}_c{self.axial_concentration:g}"
        )

    def as_dict(self) -> dict[str, float]:
        return {
            "amplitude": self.amplitude,
            "radial_support": self.radial_support,
            "axial_support": self.axial_support,
            "axial_concentration": self.axial_concentration,
            "aspect_ratio": self.aspect_ratio,
            "support_radius": self.support_radius,
        }

    # -- evaluation --------------------------------------------------------- #

    def u1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        """``u_1`` at the requested points; smooth in ``r^2`` by construction."""
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        radial = smooth_bump((r_array / self.radial_support) ** 2)
        scaled = z_array / self.axial_support
        axial = smooth_bump(scaled**2)
        profile = scaled / (1.0 + self.axial_concentration * scaled**2)
        return self.amplitude * radial * axial * profile

    def swirl_velocity(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        """The physical azimuthal velocity ``u^\\theta = r u_1``."""
        return np.asarray(r, dtype=np.float64) * self.u1(r, z)

    def initial_state(self, grid: AxisymmetricGrid) -> tuple[FloatArray, FloatArray]:
        """``(u_1, \\omega_1)`` on the grid, with ``\\omega_1 = 0`` exactly."""
        r_mesh, z_mesh = grid.mesh()
        u1 = self.u1(r_mesh, z_mesh)
        # The datum is compactly supported; pin the boundary rows so the run
        # starts exactly inside the Dirichlet frame the time-stepper imposes.
        u1[-1, :] = 0.0
        u1[:, 0] = 0.0
        u1[:, -1] = 0.0
        return u1, np.zeros_like(u1)

    def cartesian_velocity(self, points: npt.ArrayLike) -> FloatArray:
        """``u`` at Cartesian points of shape ``(n,3)``: ``u_1 (-y, x, 0)``."""
        coordinates = np.asarray(points, dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("points must have shape (n, 3)")
        x, y, z = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]
        swirl = self.u1(np.hypot(x, y), z)
        return np.stack((-y * swirl, x * swirl, np.zeros_like(swirl)), axis=1)

    # -- exact scalar invariants -------------------------------------------- #

    def kinetic_energy(self, grid: AxisymmetricGrid) -> float:
        """``\\tfrac12\\int |u|^2\\,dV`` on the grid, by the trapezoid rule."""
        r_mesh, _z_mesh = grid.mesh()
        u1, _ = self.initial_state(grid)
        density = (r_mesh * u1) ** 2 * (2.0 * math.pi * r_mesh)
        return 0.5 * float(np.sum(density) * grid.dr * grid.dz)

    def l3_norm(self, grid: AxisymmetricGrid) -> float:
        """``\\|u\\|_{L^3}`` on the grid."""
        r_mesh, _z_mesh = grid.mesh()
        u1, _ = self.initial_state(grid)
        density = np.abs(r_mesh * u1) ** 3 * (2.0 * math.pi * r_mesh)
        return float(np.sum(density) * grid.dr * grid.dz) ** (1.0 / 3.0)
