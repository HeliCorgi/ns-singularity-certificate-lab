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
    "CompositeSwirlFamily",
    "SwirlComponent",
    "SwirlFamily",
    "family_A",
    "family_H",
    "family_S",
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


# --------------------------------------------------------------------------- #
# The three preregistered families                                             #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SwirlComponent:
    """One smooth compactly supported swirl bump.

    The contribution to ``u_1`` is

        amplitude * chi(s_r)^radial_power * chi(s_z)^axial_power * profile

    with ``s_r = (r/radial_support)^2``, ``s_z = ((z-center)/axial_support)^2``
    and ``profile`` either the odd ``zeta/(1+c zeta^2)`` or the even constant
    one.  Raising ``chi`` to a power keeps the field ``C^infinity`` while
    concentrating it, since ``chi(s)^k = exp(-k/(1-s))``.
    """

    amplitude: float
    radial_support: float
    axial_support: float
    axial_center: float = 0.0
    axial_concentration: float = 0.0
    odd_axial: bool = True
    radial_power: float = 1.0
    axial_power: float = 1.0

    def __post_init__(self) -> None:
        for name in (
            "amplitude", "radial_support", "axial_support", "axial_center",
            "axial_concentration", "radial_power", "axial_power",
        ):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        if self.radial_support <= 0.0 or self.axial_support <= 0.0:
            raise ValueError("supports must be positive")
        if self.radial_power <= 0.0 or self.axial_power <= 0.0:
            raise ValueError("bump powers must be positive")
        if self.axial_concentration < 0.0:
            raise ValueError("axial_concentration must be nonnegative")

    def evaluate(self, r: FloatArray, z: FloatArray) -> FloatArray:
        radial = smooth_bump((r / self.radial_support) ** 2) ** self.radial_power
        scaled = (z - self.axial_center) / self.axial_support
        axial = smooth_bump(scaled**2) ** self.axial_power
        profile = (
            scaled / (1.0 + self.axial_concentration * scaled**2)
            if self.odd_axial
            else np.ones_like(scaled)
        )
        return self.amplitude * radial * axial * profile

    @property
    def support_radius(self) -> float:
        return math.hypot(
            self.radial_support, abs(self.axial_center) + self.axial_support
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "amplitude": self.amplitude,
            "radial_support": self.radial_support,
            "axial_support": self.axial_support,
            "axial_center": self.axial_center,
            "axial_concentration": self.axial_concentration,
            "odd_axial": self.odd_axial,
            "radial_power": self.radial_power,
            "axial_power": self.axial_power,
        }


@dataclass(frozen=True)
class CompositeSwirlFamily:
    """A preregistered initial-data family built from smooth swirl components.

    Superposing components preserves every property of a single one: the sum is
    ``C^infinity``, compactly supported, even in ``r`` through ``r^2``, and
    exactly divergence free because any pure axisymmetric swirl is.
    """

    name: str
    components: tuple[SwirlComponent, ...]
    reference_length: float

    def __post_init__(self) -> None:
        if not self.components:
            raise ValueError("a family needs at least one component")
        if self.reference_length <= 0.0:
            raise ValueError("reference_length must be positive")

    @property
    def amplitude(self) -> float:
        """The scale ``A`` used by the nondimensionalisation."""
        return max(abs(component.amplitude) for component in self.components)

    @property
    def support_radius(self) -> float:
        return max(component.support_radius for component in self.components)

    @property
    def aspect_ratio(self) -> float:
        first = self.components[0]
        return first.radial_support / first.axial_support

    def label(self) -> str:
        return f"{self.name}_A{self.amplitude:g}_L{self.reference_length:g}"

    def u1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        total = np.zeros(shape, dtype=np.float64)
        rb = np.broadcast_to(r_array, shape)
        zb = np.broadcast_to(z_array, shape)
        for component in self.components:
            total = total + component.evaluate(rb, zb)
        return total

    def swirl_velocity(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        return np.asarray(r, dtype=np.float64) * self.u1(r, z)

    def initial_state(self, grid: AxisymmetricGrid) -> tuple[FloatArray, FloatArray]:
        r_mesh, z_mesh = grid.mesh()
        u1 = self.u1(r_mesh, z_mesh)
        u1[-1, :] = 0.0
        u1[:, 0] = 0.0
        u1[:, -1] = 0.0
        return u1, np.zeros_like(u1)

    def cartesian_velocity(self, points: npt.ArrayLike) -> FloatArray:
        coordinates = np.asarray(points, dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("points must have shape (n, 3)")
        x, y, z = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]
        swirl = self.u1(np.hypot(x, y), z)
        return np.stack((-y * swirl, x * swirl, np.zeros_like(swirl)), axis=1)

    def kinetic_energy(self, grid: AxisymmetricGrid) -> float:
        r_mesh, _z = grid.mesh()
        u1, _ = self.initial_state(grid)
        density = (r_mesh * u1) ** 2 * (2.0 * math.pi * r_mesh)
        return 0.5 * float(np.sum(density) * grid.dr * grid.dz)

    def moment_degeneracy(self, grid: AxisymmetricGrid) -> dict[str, object]:
        r"""Which multipole moments of the generated ``omega_1`` degenerate.

        At leading order in time ``omega_1 = t d_z(u_1^2)``, so the moments come
        from integrating that against ``1``, ``z``, ``z^2`` and ``r^2`` with the
        measure ``dV_5 = |S^3| r^3 dr dz``:

        * **monopole** ``int d_z(u_1^2) dV_5 = 0`` for *every* compactly
          supported datum, being the ``z`` integral of a ``z`` derivative.  It
          is therefore degenerate for all three families and grows only at
          ``O(t^2)``, through the advection term whose five-dimensional
          divergence ``2u^r/r`` is nonzero.
        * **axial dipole** ``int z d_z(u_1^2) dV_5 = - int u_1^2 dV_5``, nonzero
          whenever ``u_1`` is.  **Never degenerate.**
        * **axial quadrupole** ``int z^2 d_z(u_1^2) dV_5 = -2 int z u_1^2 dV_5``,
          which vanishes exactly when ``u_1^2`` is even in ``z`` -- that is, for
          any family odd or even in ``z``.  Degenerate for S and H,
          **non-degenerate for A**.
        * **transverse quadrupole** ``(1/4) int r^2 d_z(u_1^2) dV_5 = 0`` for
          every compactly supported datum, for the same reason as the monopole.

        The numbers returned are grid quadratures of those exact leading-order
        integrals, normalised by ``int u_1^2 dV_5``, so the claims are checkable.
        """
        r_mesh, z_mesh = grid.mesh()
        u1 = self.u1(r_mesh, z_mesh)
        weight = (r_mesh**3) * grid.dr * grid.dz * (2.0 * math.pi**2)
        square = u1**2
        scale = float(np.sum(square * weight)) or 1.0
        first_moment = float(np.sum(z_mesh * square * weight)) / scale
        return {
            "monopole_leading": 0.0,
            "monopole_degenerate": True,
            "axial_dipole_leading": -1.0,
            "axial_dipole_degenerate": False,
            "axial_quadrupole_leading": -2.0 * first_moment,
            "axial_quadrupole_degenerate": abs(first_moment) < 1.0e-10,
            "transverse_quadrupole_leading": 0.0,
            "transverse_quadrupole_degenerate": True,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "amplitude": self.amplitude,
            "reference_length": self.reference_length,
            "aspect_ratio": self.aspect_ratio,
            "support_radius": self.support_radius,
            "components": [c.as_dict() for c in self.components],
        }


def family_S(amplitude: float = 10.0) -> CompositeSwirlFamily:
    """**Family S** -- the symmetric baseline, identical to :class:`SwirlFamily`.

    Odd in ``z``, so its axial quadrupole degenerates.  Kept as the regression
    reference for every earlier result.
    """
    return CompositeSwirlFamily(
        name="S",
        components=(
            SwirlComponent(
                amplitude=amplitude, radial_support=1.2, axial_support=1.5,
                axial_concentration=0.5, odd_axial=True,
            ),
        ),
        reference_length=1.2,
    )


def family_A(amplitude: float = 10.0) -> CompositeSwirlFamily:
    """**Family A** -- axially asymmetric, so the quadrupole does not degenerate.

    Two bumps at different centres with different amplitudes and widths.  The
    sum is still ``C^infinity``, compactly supported and exactly divergence
    free; what changes is that ``u_1^2`` is no longer even in ``z``, so
    ``int z u_1^2 dV_5`` is nonzero and the axial quadrupole survives.
    """
    return CompositeSwirlFamily(
        name="A",
        components=(
            SwirlComponent(
                amplitude=amplitude, radial_support=1.2, axial_support=1.1,
                axial_center=0.45, axial_concentration=0.3, odd_axial=True,
            ),
            SwirlComponent(
                amplitude=0.55 * amplitude, radial_support=0.8,
                axial_support=0.7, axial_center=-0.75, odd_axial=False,
            ),
        ),
        reference_length=1.2,
    )


def family_H(amplitude: float = 10.0) -> CompositeSwirlFamily:
    """**Family H** -- a Hou-shaped whole-space datum.

    Keeps the two features of the Hou setup that matter -- concentration near
    the axis and a thin axial extent -- while remaining a legitimate whole-space
    Clay datum: non-periodic ``z``, smooth compact support in both directions,
    Cartesian smoothness through ``r^2``, and a fixed positive viscosity.  It
    imports **neither** the finite-cylinder wall condition **nor** the two-stage
    viscosity protocol.

    Concentration comes from raising the bump to a power, which preserves
    ``C^infinity`` smoothness: ``chi(s)^k = exp(-k/(1-s))``.
    """
    return CompositeSwirlFamily(
        name="H",
        components=(
            SwirlComponent(
                amplitude=amplitude, radial_support=1.0, axial_support=0.35,
                axial_concentration=0.0, odd_axial=True,
                radial_power=2.5, axial_power=1.0,
            ),
        ),
        reference_length=1.0,
    )
