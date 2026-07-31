r"""Nondimensionalisation of the axisymmetric swirl system.

Raw runs are parameterised by ``(A, L, nu, t)``.  That is four numbers for a
problem that has two.  This module carries out the reduction once, so that no
later sweep can waste effort on runs that are the same run in disguise.

Substituting

.. math::

   r = L\rho,\quad z = L\zeta,\quad \tau = At,\quad
   u_1 = A\,U,\quad \omega_1 = (A/L)\,W,\quad \psi_1 = AL\,\Psi

into the audited system (E-11--E-14) gives, term by term:

* elliptic:  ``-L_5 psi_1 = omega_1`` becomes ``-\mathcal L_5 \Psi = W`` with
  ``\mathcal L_5 = \partial_{\rho\rho} + 3\rho^{-1}\partial_\rho +
  \partial_{\zeta\zeta}``, because ``L_5 = L^{-2}\mathcal L_5`` and both sides
  carry the factor ``A/L``;
* recovery:  ``u^r = AL\,U^\rho`` with ``U^\rho = -\rho\,\partial_\zeta\Psi``,
  ``u^z = AL\,U^\zeta`` with ``U^\zeta = 2\Psi + \rho\,\partial_\rho\Psi``;
* swirl:  every term of the ``u_1`` equation carries ``A^2`` except the viscous
  one, which carries ``\nu A/L^2``;
* vorticity:  every term of the ``\omega_1`` equation carries ``A^2/L`` except
  the viscous one, which carries ``\nu A/L^3``.

Dividing by ``A^2`` and ``A^2/L``:

.. math::

   \partial_\tau U + U^\rho\partial_\rho U + U^\zeta\partial_\zeta U
     = 2U\,\partial_\zeta\Psi + \mathrm{Re}^{-1}\mathcal L_5 U, \\
   \partial_\tau W + U^\rho\partial_\rho W + U^\zeta\partial_\zeta W
     = \partial_\zeta(U^2) + \mathrm{Re}^{-1}\mathcal L_5 W,

with the **single** parameter

.. math::   \mathrm{Re} = \frac{A L^2}{\nu}.

Everything else is dimensionless shape.  Two raw settings with equal ``Re`` and
equal dimensionless shape are the same computation, and
:func:`deduplicate_settings` finds them.

A warning that this module exists to make unavoidable.  ``\tau = At`` is the
natural time of the *swirl* variable, but pure swirl advects nothing: the
meridional velocity that does the advecting is itself generated, at order
``\tau``.  The advection-to-source ratio therefore behaves like ``C\tau^2`` with
a small ``C``, and ``\tau = O(1)`` does **not** by itself imply nonlinearity.
:func:`nonlinear_time_estimate` turns a measured ``C`` into the ``\tau`` a run
actually has to reach; ignoring it is how the previous sweep spent thirty-two
runs inside the first Picard iterate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

__all__ = [
    "ScalingExponents",
    "Setting",
    "deduplicate_settings",
    "nonlinear_time_estimate",
    "reynolds_number",
]


def reynolds_number(amplitude: float, length: float, viscosity: float) -> float:
    """``Re = A L^2 / nu``, the only parameter the system retains."""
    if viscosity <= 0.0:
        raise ValueError("a Clay candidate requires a fixed positive viscosity")
    if amplitude <= 0.0 or length <= 0.0:
        raise ValueError("amplitude and length must be positive")
    return amplitude * length * length / viscosity


@dataclass(frozen=True)
class ScalingExponents:
    r"""How each physical diagnostic scales with ``(A, L)`` at fixed ``Re``.

    Each entry is ``(a, b)`` meaning the quantity equals ``A^a L^b`` times a
    dimensionless number that depends only on ``Re``, the shape and ``\tau``.
    """

    kinetic_energy: tuple[int, int] = (2, 5)
    l3_norm: tuple[int, int] = (1, 2)
    vorticity: tuple[int, int] = (1, 0)
    front_width: tuple[int, int] = (0, 1)
    swirl_variable: tuple[int, int] = (1, 0)
    vorticity_variable: tuple[int, int] = (1, -1)
    stream_variable: tuple[int, int] = (1, 1)
    velocity: tuple[int, int] = (1, 1)
    time: tuple[int, int] = (-1, 0)

    def factor(self, name: str, amplitude: float, length: float) -> float:
        """The dimensional factor ``A^a L^b`` for the named diagnostic."""
        if not hasattr(self, name):
            raise ValueError(f"unknown diagnostic {name!r}")
        a, b = getattr(self, name)
        return amplitude**a * length**b

    def as_dict(self) -> dict[str, list[int]]:
        return {
            name: list(getattr(self, name))
            for name in (
                "kinetic_energy", "l3_norm", "vorticity", "front_width",
                "swirl_variable", "vorticity_variable", "stream_variable",
                "velocity", "time",
            )
        }


@dataclass(frozen=True)
class Setting:
    """One raw run setting together with its dimensionless coordinates."""

    label: str
    amplitude: float
    length: float
    viscosity: float
    aspect_ratio: float
    concentration: float
    physical_time: float

    @property
    def reynolds(self) -> float:
        return reynolds_number(self.amplitude, self.length, self.viscosity)

    @property
    def dimensionless_time(self) -> float:
        """``tau = A t``."""
        return self.amplitude * self.physical_time

    def shape_key(self, digits: int = 6) -> tuple[float, float]:
        return (round(self.aspect_ratio, digits), round(self.concentration, digits))

    def dimensionless_key(self, digits: int = 6) -> tuple[float, float, float, float]:
        """Two settings sharing this key are the *same* computation."""
        return (
            round(self.reynolds, digits),
            *self.shape_key(digits),
            round(self.dimensionless_time, digits),
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "amplitude": self.amplitude,
            "length": self.length,
            "viscosity": self.viscosity,
            "aspect_ratio": self.aspect_ratio,
            "concentration": self.concentration,
            "physical_time": self.physical_time,
            "reynolds": self.reynolds,
            "dimensionless_time": self.dimensionless_time,
        }


def deduplicate_settings(
    settings: Sequence[Setting], *, digits: int = 6
) -> tuple[list[Setting], list[tuple[str, str]]]:
    """Split a sweep into distinct computations and duplicate pairs.

    Returns ``(unique, duplicates)`` where ``duplicates`` lists
    ``(kept_label, dropped_label)`` pairs.  Two settings are duplicates when
    their ``Re``, their dimensionless shape **and** their dimensionless time all
    agree — the last condition matters, because a sweep run to a fixed physical
    time reaches different ``tau`` at different amplitudes.
    """
    seen: dict[tuple[float, float, float, float], Setting] = {}
    unique: list[Setting] = []
    duplicates: list[tuple[str, str]] = []
    for setting in settings:
        key = setting.dimensionless_key(digits)
        if key in seen:
            duplicates.append((seen[key].label, setting.label))
        else:
            seen[key] = setting
            unique.append(setting)
    return unique, duplicates


def nonlinear_time_estimate(
    *, measured_ratio: float, measured_tau: float, target_ratio: float = 0.1
) -> float:
    r"""The ``tau`` at which advection reaches ``target_ratio`` of the source.

    The advection-to-source ratio of a pure-swirl datum grows like ``C tau^2``:
    the vorticity grows linearly in ``tau`` from the stretching source, the
    stream function inherits that growth through the elliptic solve, and the
    advecting meridional velocity is bilinear in the two.  Measuring the ratio
    once fixes ``C`` and hence the answer

    .. math::  \tau_\ast = \tau_{\mathrm{meas}}
                 \sqrt{\text{target} / \text{measured}} .

    This is an *empirical* extrapolation of a quadratic law, not a theorem.  It
    is here because the alternative — assuming ``tau = O(1)`` suffices — is what
    kept the previous sweep inside the first Picard iterate.
    """
    if measured_ratio <= 0.0 or measured_tau <= 0.0 or target_ratio <= 0.0:
        raise ValueError("measured ratio, tau and target must all be positive")
    return measured_tau * math.sqrt(target_ratio / measured_ratio)
