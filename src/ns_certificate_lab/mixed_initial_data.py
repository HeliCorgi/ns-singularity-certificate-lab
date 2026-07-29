r"""Mixed axisymmetric initial data: swirl **and** meridional flow.

Every whole-space family the repository has run so far is a pure swirl, and
:mod:`ns_certificate_lab.l3_generation` shows why that was a dead end: the
pressure contribution to ``\frac{d}{dt}\|u\|_{L^3}^3`` vanishes identically for a
purely azimuthal field, and the advective contribution vanishes for *every*
divergence-free field, so the initial rate is forced non-positive.  Leaving that
regime requires meridional flow in the datum itself.

Construction
------------
Two scalar generators, both smooth compactly supported functions of ``(r^2, z)``:

.. math::

   u^\theta = r\,u_1(r^2,z),\qquad
   u^r = -r\,\partial_z\psi_1(r^2,z),\qquad
   u^z = 2\psi_1 + r\,\partial_r\psi_1 .

The meridional part is divergence free for **any** ``\psi_1``:

.. math::

   \partial_r u^r + \frac{u^r}{r} + \partial_z u^z
   = -r\partial_{rz}\psi_1 - \partial_z\psi_1
     + 2\partial_z\psi_1 + \partial_z\psi_1 + r\partial_{rz}\psi_1 = 0 ,

and the swirl part is divergence free because it has no ``\theta`` dependence.
So the sum is exactly divergence free, with no projection step and no residual.

Why ``(r^2, z)`` and not ``(r, z)``
-----------------------------------
The Cartesian field is

.. math::

   u = u_1(r^2,z)\,(-y,x,0)
     + \bigl(-\partial_z\psi_1(r^2,z)\bigr)(x,y,0)
     + \bigl(2\psi_1 + r\partial_r\psi_1\bigr)(0,0,1) .

The first two brackets are smooth functions of ``(x^2+y^2, z)`` multiplying
polynomials, hence ``C^\infty``.  The third needs care: ``r\partial_r\psi_1``.
Writing ``\psi_1 = P(s,z)`` with ``s = r^2`` gives ``\partial_r\psi_1 =
2r\,\partial_sP``, so ``r\partial_r\psi_1 = 2s\,\partial_sP`` — a smooth function
of ``(s,z)``, hence ``C^\infty`` in Cartesian coordinates.  Had ``\psi_1`` been a
function of ``r`` rather than ``r^2``, ``u^z`` would only be Lipschitz across the
axis.  This is the same requirement that already governs ``u_1``, and it is why
:class:`MeridionalComponent` is parameterised by ``s = (r/R)^2`` throughout.

The three preregistered families
--------------------------------
* **M1** — a weak meridional flow added to the existing swirl of family S.  The
  meridional amplitude is a free knob so the pure-swirl limit is recoverable
  exactly, which makes M1 the regression bridge to every Gate 7 result.
* **M2** — the compression/stretching pattern aligned with the swirl peak.  The
  meridional flow is placed so its axial compression sits where ``|u^\theta|`` is
  largest, which is where the pressure term has the best chance of turning
  positive.
* **M3** — axially asymmetric with a non-degenerate dipole *and* quadrupole,
  built by giving swirl and meridional parts different axial centres.

The axial-parity selection rule
-------------------------------
A datum whose generators have definite parity in ``z`` obeys a rule that decides
in advance whether the pressure term of
:func:`~ns_certificate_lab.l3_generation.l3_generation_rate` can be nonzero at
all.  Let ``\sigma_u`` and ``\sigma_\psi`` be the parities of ``u_1`` and
``\psi_1``.  Then

* ``u^\theta`` has parity ``\sigma_u``, ``u^r`` has parity ``-\sigma_\psi``,
  ``u^z`` has parity ``\sigma_\psi``;
* ``|u|^2`` is a sum of squares, hence **always even**, so ``|u|`` is even;
* every term of ``\partial_iu_j\partial_ju_i`` is either a square or one of the
  two cross products ``(\partial_ru^\theta)(u^\theta/r)`` and
  ``(\partial_ru^z)(\partial_zu^r)``, each a product of two factors of equal
  parity, so the pressure source is **always even** and therefore so is ``p``;
* the pressure integrand ``-3|u|(u^r\partial_rp + u^z\partial_zp)`` then has
  parity ``-\sigma_\psi``.

Integrated over a ``z``-symmetric domain it therefore **vanishes identically
unless** ``\sigma_\psi = -1``: the stream generator must be **odd** in ``z``.
The swirl parity does not enter.  See
:func:`~ns_certificate_lab.l3_generation.parity_selection_rule`.

An odd ``\psi_1`` is also the physically right choice.  It gives ``u^r`` even and
``u^z`` odd, which is the axisymmetric strain flow that converges radially at
every height and ejects along both ``\pm z`` — the configuration that stretches
vortex lines.  An even ``\psi_1`` gives inflow on one side and outflow on the
other, and the two halves cancel.  All three families below therefore use odd
stream generators, or break the parity outright.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .initial_data import smooth_bump, smooth_bump_derivative

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "GeneratorComponent",
    "MixedFamily",
    "family_M1",
    "family_M2",
    "family_M3",
    "MIXED_FAMILIES",
]


@dataclass(frozen=True)
class GeneratorComponent:
    r"""One smooth compactly supported generator bump.

    .. math::

       g(r,z) = a\,\chi(s_r)^{p}\,\chi(s_z)^{q}\,\Pi(\zeta),
       \quad s_r = \frac{r^2}{R^2},\ \ s_z = \zeta^2,\ \
       \zeta = \frac{z - z_0}{Z} ,

    with ``\Pi`` either the odd ``\zeta/(1+c\zeta^2)`` or the even ``1``.
    Raising ``\chi`` to a power keeps the result ``C^\infty`` because
    ``\chi(s)^k = \exp(-k/(1-s))``, and concentrates it.

    Everything is a function of ``r^2``, never of ``r``.
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
        s_r = (r / self.radial_support) ** 2
        zeta = (z - self.axial_center) / self.axial_support
        radial = smooth_bump(s_r) ** self.radial_power
        axial = smooth_bump(zeta**2) ** self.axial_power
        profile = (
            zeta / (1.0 + self.axial_concentration * zeta**2)
            if self.odd_axial
            else np.ones_like(zeta)
        )
        return self.amplitude * radial * axial * profile

    def d_ds(self, r: FloatArray, z: FloatArray) -> FloatArray:
        r"""``\partial g/\partial s`` with ``s = r^2``.

        Needed because ``u^z`` contains ``r\partial_r\psi_1 = 2s\,\partial_s
        \psi_1``, and evaluating that through ``\partial_s`` rather than
        ``\partial_r`` is what makes the axis value exact rather than a
        one-sided difference.
        """
        radius_sq = self.radial_support**2
        s_r = (r * r) / radius_sq
        chi = smooth_bump(s_r)
        # d/ds [chi(s/R^2)^p] = p chi^{p-1} chi'(s/R^2) / R^2
        radial_derivative = np.where(
            chi > 0.0,
            self.radial_power
            * chi ** (self.radial_power - 1.0)
            * smooth_bump_derivative(s_r)
            / radius_sq,
            0.0,
        )
        zeta = (z - self.axial_center) / self.axial_support
        axial = smooth_bump(zeta**2) ** self.axial_power
        profile = (
            zeta / (1.0 + self.axial_concentration * zeta**2)
            if self.odd_axial
            else np.ones_like(zeta)
        )
        return self.amplitude * radial_derivative * axial * profile

    # -- analytic partial derivatives ---------------------------------------- #
    #
    # Every component factors as ``a X(s) Z(z)`` with ``s = r^2``, so every mixed
    # partial factors too: ``d_s^m d_z^n g = a X^{(m)}(s) Z^{(n)}(z)``.  That is
    # why the generators are built as a product and not as a general function of
    # two variables -- it makes the exact velocity gradient a few lines of
    # algebra instead of an automatic-differentiation dependency, and it is what
    # lets the constructed field be divergence free to machine precision rather
    # than to the truncation error of a difference operator.

    def _radial_factors(self, r: FloatArray) -> tuple[FloatArray, FloatArray, FloatArray]:
        r"""``(X, X_s, X_{ss})`` with ``X(s) = \chi(s/R^2)^p``, ``s = r^2``."""
        radius_sq = self.radial_support**2
        sigma = (r * r) / radius_sq
        chi = smooth_bump(sigma)
        inside = chi > 0.0
        gap = np.where(inside, 1.0 - sigma, 1.0)
        chi_1 = np.where(inside, -chi / gap**2, 0.0)
        chi_2 = np.where(inside, chi * (1.0 / gap**4 - 2.0 / gap**3), 0.0)
        power = self.radial_power
        base = np.where(inside, chi, 1.0)
        value = np.where(inside, chi**power, 0.0)
        first = np.where(inside, power * base ** (power - 1.0) * chi_1, 0.0)
        second = np.where(
            inside,
            power * (power - 1.0) * base ** (power - 2.0) * chi_1**2
            + power * base ** (power - 1.0) * chi_2,
            0.0,
        )
        return value, first / radius_sq, second / radius_sq**2

    def _axial_factors(self, z: FloatArray) -> tuple[FloatArray, FloatArray, FloatArray]:
        r"""``(Z, Z_z, Z_{zz})`` with ``Z(z) = \chi(\zeta^2)^q\,\Pi(\zeta)``."""
        width = self.axial_support
        zeta = (z - self.axial_center) / width
        argument = zeta**2
        chi = smooth_bump(argument)
        inside = chi > 0.0
        gap = np.where(inside, 1.0 - argument, 1.0)
        chi_1 = np.where(inside, -chi / gap**2, 0.0)
        chi_2 = np.where(inside, chi * (1.0 / gap**4 - 2.0 / gap**3), 0.0)
        power = self.axial_power
        base = np.where(inside, chi, 1.0)
        bump = np.where(inside, chi**power, 0.0)
        # d/dzeta chi(zeta^2)^q  and its second derivative.
        d_chi = chi_1 * 2.0 * zeta
        d2_chi = chi_2 * 4.0 * argument + chi_1 * 2.0
        bump_1 = np.where(inside, power * base ** (power - 1.0) * d_chi, 0.0)
        bump_2 = np.where(
            inside,
            power * (power - 1.0) * base ** (power - 2.0) * d_chi**2
            + power * base ** (power - 1.0) * d2_chi,
            0.0,
        )
        if self.odd_axial:
            c = self.axial_concentration
            denominator = 1.0 + c * argument
            profile = zeta / denominator
            profile_1 = (1.0 - c * argument) / denominator**2
            profile_2 = -2.0 * c * zeta * (3.0 - c * argument) / denominator**3
        else:
            profile = np.ones_like(zeta)
            profile_1 = np.zeros_like(zeta)
            profile_2 = np.zeros_like(zeta)
        value = bump * profile
        first = bump_1 * profile + bump * profile_1
        second = bump_2 * profile + 2.0 * bump_1 * profile_1 + bump * profile_2
        return value, first / width, second / width**2

    def partials(self, r: FloatArray, z: FloatArray) -> dict[str, FloatArray]:
        r"""``value``, ``ds``, ``dss``, ``dz``, ``dzz``, ``dsz`` -- all exact."""
        x, x_s, x_ss = self._radial_factors(r)
        y, y_z, y_zz = self._axial_factors(z)
        a = self.amplitude
        return {
            "value": a * x * y,
            "ds": a * x_s * y,
            "dss": a * x_ss * y,
            "dz": a * x * y_z,
            "dzz": a * x * y_zz,
            "dsz": a * x_s * y_z,
        }

    def d_dz(self, r: FloatArray, z: FloatArray) -> FloatArray:
        r"""``\partial g/\partial z``, analytically.

        ``u^r = -r\partial_z\psi_1`` is the only velocity component that needs
        an axial derivative of a generator.  Taking it analytically rather than
        by finite differences means the constructed field is divergence free as
        an *exact algebraic identity*, so the divergence the diagnostics report
        is the finite-difference operator's error and nothing else.
        """
        radial = smooth_bump((r / self.radial_support) ** 2) ** self.radial_power
        width = self.axial_support
        zeta = (z - self.axial_center) / width
        chi = smooth_bump(zeta**2)
        axial = chi**self.axial_power
        d_axial = np.where(
            chi > 0.0,
            self.axial_power
            * chi ** (self.axial_power - 1.0)
            * smooth_bump_derivative(zeta**2)
            * 2.0
            * zeta
            / width,
            0.0,
        )
        if self.odd_axial:
            c = self.axial_concentration
            denominator = 1.0 + c * zeta**2
            profile = zeta / denominator
            d_profile = (1.0 - c * zeta**2) / (denominator**2 * width)
        else:
            profile = np.ones_like(zeta)
            d_profile = np.zeros_like(zeta)
        return self.amplitude * radial * (d_axial * profile + axial * d_profile)

    @property
    def support_radius(self) -> float:
        return math.hypot(
            self.radial_support, abs(self.axial_center) + self.axial_support
        )

    @property
    def axial_parity(self) -> str:
        """``'odd'``, ``'even'`` or ``'none'`` -- the parity in ``z``."""
        if self.axial_center != 0.0:
            return "none"
        return "odd" if self.odd_axial else "even"

    def as_dict(self) -> dict[str, object]:
        payload = {name: getattr(self, name) for name in self.__dataclass_fields__}
        payload["axial_parity"] = self.axial_parity
        return payload


@dataclass(frozen=True)
class MixedFamily:
    r"""A preregistered mixed family: a swirl generator and a stream generator."""

    name: str
    swirl: tuple[GeneratorComponent, ...]
    stream: tuple[GeneratorComponent, ...]
    reference_length: float

    def __post_init__(self) -> None:
        if not self.swirl and not self.stream:
            raise ValueError("a family needs at least one generator component")
        if self.reference_length <= 0.0:
            raise ValueError("reference_length must be positive")

    # -- generators ---------------------------------------------------------- #

    def u1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        return self._sum(self.swirl, r, z, "evaluate")

    def psi1(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        return self._sum(self.stream, r, z, "evaluate")

    def dpsi1_ds(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        return self._sum(self.stream, r, z, "d_ds")

    def dpsi1_dz(self, r: npt.ArrayLike, z: npt.ArrayLike) -> FloatArray:
        return self._sum(self.stream, r, z, "d_dz")

    def partials(
        self, components, r: npt.ArrayLike, z: npt.ArrayLike
    ) -> dict[str, FloatArray]:
        """Summed exact partials of a generator group."""
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        rb = np.broadcast_to(r_array, shape)
        zb = np.broadcast_to(z_array, shape)
        keys = ("value", "ds", "dss", "dz", "dzz", "dsz")
        total = {key: np.zeros(shape, dtype=np.float64) for key in keys}
        for component in components:
            piece = component.partials(rb, zb)
            for key in keys:
                total[key] = total[key] + piece[key]
        return total

    def exact_gradient(self, grid: AxisymmetricGrid) -> dict[str, FloatArray]:
        r"""The nine physical components of ``\nabla u``, analytically.

        With ``s = r^2`` so that ``\partial_r = 2r\,\partial_s``:

        .. math::

           \partial_r u^\theta = u_1 + 2s\,\partial_su_1, \quad
           \partial_z u^\theta = r\,\partial_zu_1, \quad
           u^\theta/r = u_1, \\
           \partial_r u^r = -(\partial_z\psi_1 + 2s\,\partial_s\partial_z\psi_1),
           \quad \partial_z u^r = -r\,\partial_{zz}\psi_1, \quad
           u^r/r = -\partial_z\psi_1, \\
           \partial_r u^z = 2r\,(4\partial_s\psi_1 + 2s\,\partial_{ss}\psi_1),
           \quad \partial_z u^z
             = 2\partial_z\psi_1 + 2s\,\partial_s\partial_z\psi_1 .

        Adding the ``\theta`` row gives ``\nabla\cdot u = \partial_ru^r +
        u^r/r + \partial_zu^z = -2\partial_z\psi_1 - 2s\partial_s\partial_z\psi_1
        + 2\partial_z\psi_1 + 2s\partial_s\partial_z\psi_1 = 0`` exactly, at
        every point including the axis.
        """
        r_mesh, z_mesh = grid.mesh()
        swirl = self.partials(self.swirl, r_mesh, z_mesh)
        stream = self.partials(self.stream, r_mesh, z_mesh)
        s = r_mesh**2
        return {
            "rr": -(stream["dz"] + 2.0 * s * stream["dsz"]),
            "rt": swirl["value"] + 2.0 * s * swirl["ds"],
            "rz": 2.0 * r_mesh * (4.0 * stream["ds"] + 2.0 * s * stream["dss"]),
            "tr": -swirl["value"],
            "tt": -stream["dz"],
            "tz": np.zeros_like(r_mesh),
            "zr": -r_mesh * stream["dzz"],
            "zt": r_mesh * swirl["dz"],
            "zz": 2.0 * stream["dz"] + 2.0 * s * stream["dsz"],
        }

    def _sum(self, components, r, z, method: str) -> FloatArray:
        r_array = np.asarray(r, dtype=np.float64)
        z_array = np.asarray(z, dtype=np.float64)
        shape = np.broadcast(r_array, z_array).shape
        total = np.zeros(shape, dtype=np.float64)
        rb = np.broadcast_to(r_array, shape)
        zb = np.broadcast_to(z_array, shape)
        for component in components:
            total = total + getattr(component, method)(rb, zb)
        return total

    # -- the field ----------------------------------------------------------- #

    def field(self, grid: AxisymmetricGrid):
        r"""The :class:`~ns_certificate_lab.l3_generation.MixedField`.

        Both meridional components come from **analytic** derivatives of the
        generators: ``u^r = -r\partial_z\psi_1`` and
        ``u^z = 2\psi_1 + 2s\,\partial_s\psi_1``.  The field is then exactly
        divergence free as an algebraic identity, so any divergence the
        diagnostics report is the finite-difference operator's truncation error
        and not a defect of the datum.  Using ``\partial_s`` rather than
        ``\partial_r`` also makes the axis value exact instead of a one-sided
        difference.
        """
        from .l3_generation import MixedField

        r_mesh, z_mesh = grid.mesh()
        swirl = self.partials(self.swirl, r_mesh, z_mesh)
        stream = self.partials(self.stream, r_mesh, z_mesh)
        u_theta = r_mesh * swirl["value"]
        u_r = -r_mesh * stream["dz"]
        u_z = 2.0 * stream["value"] + 2.0 * (r_mesh**2) * stream["ds"]
        return MixedField(
            grid=grid, u_r=u_r, u_theta=u_theta, u_z=u_z,
            analytic_gradient=self.exact_gradient(grid),
        )

    def initial_state(self, grid: AxisymmetricGrid) -> tuple[FloatArray, FloatArray]:
        r"""``(u_1, \omega_1)`` for the evolution, with ``\omega_1 = -\mathcal L_5\psi_1``."""
        from .operators import laplacian_5d_formal

        r_mesh, z_mesh = grid.mesh()
        u1 = self.u1(r_mesh, z_mesh)
        psi1 = self.psi1(r_mesh, z_mesh)
        omega1 = -laplacian_5d_formal(grid, psi1)
        for field in (u1, omega1):
            field[-1, :] = 0.0
            field[:, 0] = 0.0
            field[:, -1] = 0.0
        return u1, omega1

    def cartesian_velocity(self, points: npt.ArrayLike) -> FloatArray:
        """Cartesian velocity at arbitrary points, from the analytic generators."""
        coordinates = np.asarray(points, dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("points must have shape (n, 3)")
        x, y, z = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]
        radius = np.hypot(x, y)
        swirl = self.u1(radius, z)
        dpsi_dz = self.dpsi1_dz(radius, z)
        stream = self.psi1(radius, z)
        ds = self.dpsi1_ds(radius, z)
        return np.stack(
            (
                -y * swirl - x * dpsi_dz,
                x * swirl - y * dpsi_dz,
                2.0 * stream + 2.0 * radius**2 * ds,
            ),
            axis=1,
        )

    # -- bookkeeping --------------------------------------------------------- #

    @property
    def support_radius(self) -> float:
        return max(
            component.support_radius
            for component in (*self.swirl, *self.stream)
        )

    @property
    def meridional_amplitude(self) -> float:
        return max((abs(c.amplitude) for c in self.stream), default=0.0)

    @property
    def swirl_amplitude(self) -> float:
        return max((abs(c.amplitude) for c in self.swirl), default=0.0)

    def scaled(self, *, swirl: float = 1.0, stream: float = 1.0) -> "MixedFamily":
        """A copy with the two generator amplitudes rescaled independently."""
        return replace(
            self,
            swirl=tuple(
                replace(c, amplitude=c.amplitude * swirl) for c in self.swirl
            ),
            stream=tuple(
                replace(c, amplitude=c.amplitude * stream) for c in self.stream
            ),
        )

    def with_amplitudes(self, values: npt.ArrayLike) -> "MixedFamily":
        """A copy with every generator amplitude replaced, in order.

        The optimiser's design variables are exactly these amplitudes, so the
        search space is a finite set of coefficients over a fixed
        compactly-supported basis and every point of it is automatically
        ``C^\\infty``, divergence free and compactly supported.
        """
        array = np.asarray(values, dtype=np.float64).ravel()
        total = len(self.swirl) + len(self.stream)
        if array.size != total:
            raise ValueError(f"expected {total} amplitudes, received {array.size}")
        swirl = tuple(
            replace(c, amplitude=float(a))
            for c, a in zip(self.swirl, array[: len(self.swirl)])
        )
        stream = tuple(
            replace(c, amplitude=float(a))
            for c, a in zip(self.stream, array[len(self.swirl) :])
        )
        return replace(self, swirl=swirl, stream=stream)

    @property
    def amplitudes(self) -> FloatArray:
        return np.array(
            [c.amplitude for c in (*self.swirl, *self.stream)], dtype=np.float64
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "reference_length": self.reference_length,
            "support_radius": self.support_radius,
            "swirl_amplitude": self.swirl_amplitude,
            "meridional_amplitude": self.meridional_amplitude,
            "swirl": [c.as_dict() for c in self.swirl],
            "stream": [c.as_dict() for c in self.stream],
        }


# --------------------------------------------------------------------------- #
# the three preregistered families                                             #
# --------------------------------------------------------------------------- #


def family_M1(
    swirl_amplitude: float = 10.0, meridional_amplitude: float = 1.0
) -> MixedFamily:
    r"""**M1** — a weak meridional flow added to the existing swirl.

    The swirl generator is byte-identical to family S, so setting
    ``meridional_amplitude = 0`` recovers the Gate 7 baseline exactly and the
    whole of M1 is a one-parameter deformation away from it.  That makes M1 the
    family in which "what does meridional flow buy" is a controlled question
    rather than a change of everything at once.
    """
    return MixedFamily(
        name="M1",
        swirl=(
            GeneratorComponent(
                amplitude=swirl_amplitude, radial_support=1.2, axial_support=1.5,
                axial_concentration=0.5, odd_axial=True,
            ),
        ),
        stream=(
            GeneratorComponent(
                amplitude=meridional_amplitude, radial_support=1.2,
                axial_support=1.5, odd_axial=True,
            ),
        ),
        reference_length=1.2,
    )


def family_M2(
    swirl_amplitude: float = 10.0, meridional_amplitude: float = 4.0
) -> MixedFamily:
    r"""**M2** — compression aligned with the swirl peak.

    An **odd** stream generator gives ``u^r`` even and ``u^z`` odd: radial
    convergence at every height, ejection along both ``\pm z``.  That is the
    axisymmetric strain flow, and it is also the only parity for which the
    pressure term survives integration (see the module docstring).  The axial
    strain is then extremal on the mid-plane, where the swirl generator of
    family S has its steepest gradient: a strong swirl sitting in a strong
    strain is where the pressure integrand ``-|u|\,u\cdot\nabla p`` has the best
    chance of beating viscosity.

    Two stream bumps of opposite sign and different widths concentrate the
    strain instead of spreading it over the whole support.
    """
    return MixedFamily(
        name="M2",
        swirl=(
            GeneratorComponent(
                amplitude=swirl_amplitude, radial_support=1.2, axial_support=1.5,
                axial_concentration=0.5, odd_axial=True,
            ),
        ),
        stream=(
            GeneratorComponent(
                amplitude=meridional_amplitude, radial_support=0.9,
                axial_support=0.8, odd_axial=True, radial_power=1.0,
            ),
            GeneratorComponent(
                amplitude=-0.4 * meridional_amplitude, radial_support=1.2,
                axial_support=1.4, odd_axial=True,
            ),
        ),
        reference_length=1.2,
    )


def family_M3(
    swirl_amplitude: float = 10.0, meridional_amplitude: float = 4.0
) -> MixedFamily:
    r"""**M3** — axially asymmetric, with dipole *and* quadrupole non-degenerate.

    The axial quadrupole of the generated ``\omega_1`` degenerates whenever
    ``u_1^2`` is even in ``z``; the dipole never does.  Offsetting the swirl
    components from each other breaks the first degeneracy, and offsetting the
    stream generator from the swirl breaks the alignment symmetry as well, so
    the meridional flow is stronger on one side of the mid-plane than the other.
    """
    return MixedFamily(
        name="M3",
        swirl=(
            GeneratorComponent(
                amplitude=swirl_amplitude, radial_support=1.2, axial_support=1.1,
                axial_center=0.45, axial_concentration=0.3, odd_axial=True,
            ),
            GeneratorComponent(
                amplitude=0.55 * swirl_amplitude, radial_support=0.8,
                axial_support=0.7, axial_center=-0.75, odd_axial=False,
            ),
        ),
        stream=(
            GeneratorComponent(
                amplitude=meridional_amplitude, radial_support=1.0,
                axial_support=0.9, axial_center=0.3, odd_axial=True,
            ),
            GeneratorComponent(
                amplitude=-0.6 * meridional_amplitude, radial_support=0.7,
                axial_support=0.6, axial_center=-0.5, odd_axial=True,
            ),
        ),
        reference_length=1.2,
    )


#: The preregistered family builders, by name.
MIXED_FAMILIES = {"M1": family_M1, "M2": family_M2, "M3": family_M3}
