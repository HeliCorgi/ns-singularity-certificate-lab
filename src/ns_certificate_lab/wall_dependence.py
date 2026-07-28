r"""E-32 :math:`C^\infty` compact-support initial-data family for wall dependence.

``docs/wall_dependence_prereg.md`` fixes the family used to probe whether the
Hou growth mechanism depends on the wall at ``r = 1``.  The E-29 datum may not
be extended past ``r = 1``: ``(1-r^2)^18`` *grows* there and is a different
function, and a plain zero extension is only ``C^17`` at ``r = 1``.  E-32
instead multiplies E-29 by a standard smooth partition-of-unity cutoff of
``rho = r^2``,

.. math::

   \tilde u_1(0,r,z) = u_1^{E29}(0,r,z)\,\chi_c(r^2),
   \qquad \tilde\omega_1(0,r,z) = 0,
   \tag{E-32a}

.. math::

   \chi_c(\rho) =
   \begin{cases}
   1, & \rho \le \rho_1,\\
   \dfrac{\theta\!\left(\frac{\rho_2-\rho}{\rho_2-\rho_1}\right)}
         {\theta\!\left(\frac{\rho_2-\rho}{\rho_2-\rho_1}\right)
          + \theta\!\left(\frac{\rho-\rho_1}{\rho_2-\rho_1}\right)},
   & \rho_1 < \rho < \rho_2,\\
   0, & \rho \ge \rho_2,
   \end{cases}
   \qquad
   \theta(s) = \begin{cases} e^{-1/s}, & s > 0\\ 0, & s \le 0\end{cases}
   \tag{E-32b}

with the audited defaults ``rho1 = 0.81`` (``r = 0.9``) and
``rho2 = 0.9025`` (``r = 0.95``).

Why ``rho = r^2``
-----------------

Defining the cutoff as a function of ``r^2`` makes it automatically even in
``r``, so the E-16 axis regularity of the E-29 datum survives unchanged.

What this module guarantees
---------------------------

* **core bit-identity.**  For ``r <= 0.9`` the cutoff is the literal constant
  ``1.0`` taken from a branch, and multiplication by ``1.0`` is exact in IEEE
  arithmetic, so :func:`envelope_initial_swirl` returns *bit-identical* values
  to :func:`~ns_certificate_lab.nonlinear_cylinder.hou_initial_swirl` on the
  core.  Every wall radius therefore starts from the identical core datum.
* **compact support.**  For ``r >= 0.95`` the cutoff is the literal ``0.0``,
  so the datum vanishes exactly, strictly inside every wall
  ``R_wall >= 1``.
* **E-32c deviation bound.**  ``0 <= chi_c <= 1`` gives
  ``sup |u~1 - u1^E29| <= 12000 (1-rho1)^18 max_z |g(z)| = 3.4008e-10`` on
  ``r <= 1``, a relative perturbation of ``1.0413e-13`` against the field
  maximum ``3265.9863``.

Nothing in this module is a proof or evidence of singularity formation.  It is
initial data for a numerical observation.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid, FloatArray
from .nonlinear_cylinder import (
    cartesian_vorticity,
    constrain_state,
    hou_initial_swirl,
    relative_axis_parity,
)

# ------------------------------------------------------------- E-32 constants

#: E-29 amplitude; part of the audited equation, not a tunable input.
E29_AMPLITUDE = 12000.0

#: E-32b default inner plateau radius squared (``r = 0.9``).
E32_RHO1 = 0.81
#: E-32b default outer support radius squared (``r = 0.95``).
E32_RHO2 = 0.9025
#: ``sqrt(E32_RHO1)``: the datum is bit-identical to E-29 inside this radius.
E32_CORE_RADIUS = 0.9
#: ``sqrt(E32_RHO2)``: the datum is exactly zero outside this radius.
E32_SUPPORT_RADIUS = 0.95

#: ``max_z |sin(2 pi z) / (1 + 12.5 sin^2(pi z))|`` (E-32 property 4).
E32_MAX_ABS_AXIAL_FACTOR = 0.272165526975908
#: the ``z`` at which that maximum is attained.
E32_Z_OF_MAX_ABS_AXIAL_FACTOR = 0.0845842
#: ``12000 * max|g|``; identical to the E-29b swirl norm.
E32_AMPLITUDE_TIMES_MAX_ABS_AXIAL_FACTOR = 3265.986323710896
#: ``(1 - rho1)^18``.
E32_TRANSITION_DECAY = 1.0412735029791071e-13
#: E-32c: ``12000 (1-rho1)^18 max|g|``, the sup deviation bound on ``r <= 1``.
E32_SUP_DEVIATION_BOUND = 3.400785019972301e-10
#: E-32c relative form, against the field maximum ``3265.9863``.
E32_SUP_DEVIATION_RELATIVE_BOUND = 1.0412735029791071e-13
#: The preregistered acceptance bound of ``docs/wall_dependence_prereg.md`` (ii).
E32_SUP_DEVIATION_ACCEPTANCE = 5.0e-10

#: E-32 property 7 spot check: the radial band the fourth difference is taken on.
E32_TRANSITION_BAND = (0.85, 1.0)
#: Audited fourth-derivative estimates in that band at the ``z`` of ``max|g|``.
E32_BAND_FOURTH_DIFFERENCE = {193: 21.739, 385: 26.362}
#: Preregistration section 2 check (v) asks for *boundedness* of the high-order
#: difference across the transition band, not for exact equality.  At the
#: production core spacings ``dr = 1/192`` and ``dr = 1/128`` the band maximum
#: sits at a stencil whose whole support has ``chi_c == 1``, so the measured
#: relative change is exactly zero (E-32 property 7).  On a coarse grid the
#: band maximum's stencil does reach into the transition, where the E-29
#: amplitude is already ``1e-13`` relative; the change stays far below this
#: tolerance (``3.8e-5`` at ``dr = 1/32``).
E32_FOURTH_DIFFERENCE_RELATIVE_TOLERANCE = 1.0e-3


# ---------------------------------------------------------------- primitives


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.floating)):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _cutoff_parameters(rho1: object, rho2: object) -> tuple[float, float]:
    lower = _finite(rho1, name="rho1")
    upper = _finite(rho2, name="rho2")
    if lower < 0.0:
        raise ValueError("rho1 must be nonnegative: rho = r^2")
    if lower >= upper:
        raise ValueError("rho1 must be strictly below rho2")
    return lower, upper


def smooth_transition(s: npt.ArrayLike) -> FloatArray:
    r"""Return the standard flat function ``theta(s) = e^{-1/s}`` for ``s > 0``.

    ``theta(s) = 0`` for ``s <= 0``.  Every derivative vanishes at the origin,
    which is what makes the partition of unity built from it ``C^infinity``.
    """

    values = np.asarray(s, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError("smooth_transition requires finite arguments")
    out = np.zeros(values.shape, dtype=np.float64)
    positive = values > 0.0
    out[positive] = np.exp(-1.0 / values[positive])
    return out


def smooth_cutoff(
    rho: npt.ArrayLike,
    rho1: float = E32_RHO1,
    rho2: float = E32_RHO2,
) -> FloatArray | float:
    r"""Return the E-32b cutoff ``chi_c(rho)``.

    The three branches are evaluated exactly as written in E-32b: the plateau
    returns the literal ``1.0`` and the tail the literal ``0.0``, so a field
    multiplied by this cutoff is bit-identical to the unmultiplied field on
    ``rho <= rho1`` and exactly zero on ``rho >= rho2``.

    A scalar argument returns a Python ``float``; an array argument returns an
    array of the same shape.
    """

    lower, upper = _cutoff_parameters(rho1, rho2)
    values = np.asarray(rho, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError("smooth_cutoff requires finite arguments")
    out = np.zeros(values.shape, dtype=np.float64)
    plateau = values <= lower
    out[plateau] = 1.0
    band = (values > lower) & (values < upper)
    if np.any(band):
        width = upper - lower
        inner = smooth_transition((upper - values[band]) / width)
        outer = smooth_transition((values[band] - lower) / width)
        out[band] = inner / (inner + outer)
    if out.ndim == 0:
        return float(out)
    return out


# ------------------------------------------------------------- initial datum


def envelope_initial_swirl(
    grid: AxisymmetricGrid,
    *,
    amplitude_scale: float = 1.0,
    rho1: float = E32_RHO1,
    rho2: float = E32_RHO2,
) -> FloatArray:
    r"""Return the E-32a initial swirl ``u~1(0,r,z) = u1^{E29} chi_c(r^2)``.

    The companion condition is ``omega_1(0) = 0``, exactly as in E-29, so
    ``psi_1(0) = 0`` and ``u^r(0) = u^z(0) = 0``.

    ``rho2`` may not exceed ``r_max^2``: the cutoff has to close *inside* the
    domain, otherwise the datum would not be compactly supported strictly
    inside the wall and property 5 of E-32 would fail.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    lower, upper = _cutoff_parameters(rho1, rho2)
    scale = _finite(amplitude_scale, name="amplitude_scale")
    if scale <= 0.0:
        raise ValueError("amplitude_scale must be positive")
    r_max = float(grid.r[-1])
    if upper > r_max * r_max:
        raise ValueError(
            "rho2 must not exceed r_max^2: the E-32 cutoff must close strictly "
            "inside the domain"
        )
    base = hou_initial_swirl(
        grid,
        amplitude=E29_AMPLITUDE,
        amplitude_scale=scale,
    )
    radius, _ = grid.mesh()
    cutoff = smooth_cutoff(radius * radius, lower, upper)
    return np.asarray(base * cutoff, dtype=np.float64)


def core_index_count(grid: AxisymmetricGrid, *, core_radius: float = E32_CORE_RADIUS) -> int:
    """Return the number of leading radial rows inside ``r <= core_radius``."""

    radius = _finite(core_radius, name="core_radius")
    if radius <= 0.0:
        raise ValueError("core_radius must be positive")
    count = int(np.count_nonzero(grid.r <= radius))
    if count < 2:
        raise ValueError("the core must contain at least two radial rows")
    return count


def core_mask(
    grid: AxisymmetricGrid,
    *,
    core_radius: float = E32_CORE_RADIUS,
) -> npt.NDArray[np.bool_]:
    """Return the boolean radial mask of ``Omega_c = {r <= core_radius}``."""

    return np.asarray(grid.r <= _finite(core_radius, name="core_radius"))


# ------------------------------------------------- analytic radial profiles
#
# The E-32 property-7 spot check is stated at the *continuum* z of ``max |g|``
# rather than at a grid node, so these helpers evaluate the analytic radial
# profile at an arbitrary z on an arbitrary radial sample.


def axial_factor(z: npt.ArrayLike) -> FloatArray | float:
    """Return the E-29 axial factor ``g(z) = sin(2 pi z)/(1+12.5 sin^2(pi z))``."""

    values = np.asarray(z, dtype=np.float64)
    out = np.sin(2.0 * np.pi * values) / (
        1.0 + 12.5 * np.sin(np.pi * values) ** 2
    )
    if out.ndim == 0:
        return float(out)
    return np.asarray(out, dtype=np.float64)


def e29_radial_profile(
    radii: npt.ArrayLike,
    z: float,
    *,
    amplitude_scale: float = 1.0,
) -> FloatArray:
    """Return the plain E-29 radial profile ``u1(0, r, z)`` at one fixed ``z``."""

    r = np.asarray(radii, dtype=np.float64)
    scale = E29_AMPLITUDE * _finite(amplitude_scale, name="amplitude_scale")
    return np.asarray(
        scale * (1.0 - r * r) ** 18 * float(axial_factor(float(z))),
        dtype=np.float64,
    )


def envelope_radial_profile(
    radii: npt.ArrayLike,
    z: float,
    *,
    amplitude_scale: float = 1.0,
    rho1: float = E32_RHO1,
    rho2: float = E32_RHO2,
) -> FloatArray:
    """Return the E-32a radial profile at one fixed ``z``."""

    r = np.asarray(radii, dtype=np.float64)
    lower, upper = _cutoff_parameters(rho1, rho2)
    base = e29_radial_profile(r, z, amplitude_scale=amplitude_scale)
    return np.asarray(base * smooth_cutoff(r * r, lower, upper), dtype=np.float64)


def fourth_difference(values: npt.ArrayLike) -> FloatArray:
    """Return the undivided fourth central difference of a 1-D sample.

    ``D4 f_i = f_{i+2} - 4 f_{i+1} + 6 f_i - 4 f_{i-1} + f_{i-2}``; the result
    is indexed by the stencil centres, i.e. it has four fewer entries than the
    input.
    """

    f = np.asarray(values, dtype=np.float64)
    if f.ndim != 1 or f.size < 5:
        raise ValueError("fourth_difference needs a 1-D sample of at least five points")
    return np.asarray(
        f[4:] - 4.0 * f[3:-1] + 6.0 * f[2:-2] - 4.0 * f[1:-3] + f[:-4],
        dtype=np.float64,
    )


def transition_band_fourth_difference(
    radii: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    band: tuple[float, float] = E32_TRANSITION_BAND,
) -> dict[str, float]:
    """Return the largest fourth-derivative estimate inside the transition band.

    The reported ``maximum`` is ``max |D4 f| / dr^4`` over the stencil centres
    strictly inside ``band``.  E-32 property 7 asserts that this number is the
    *same* for the plain E-29 datum and for the E-32 envelope, because the
    largest value sits at the inner edge of the band where ``chi_c`` is exactly
    one.
    """

    r = np.asarray(radii, dtype=np.float64)
    f = np.asarray(values, dtype=np.float64)
    if r.ndim != 1 or f.shape != r.shape:
        raise ValueError("radii and values must be matching 1-D samples")
    lower, upper = (_finite(band[0], name="band lower"), _finite(band[1], name="band upper"))
    if lower >= upper:
        raise ValueError("the transition band must be nonempty")
    spacing = float(r[1] - r[0])
    difference = fourth_difference(f)
    centres = r[2:-2]
    inside = (centres > lower) & (centres < upper)
    if not np.any(inside):
        raise ValueError("the transition band contains no stencil centre")
    scaled = np.abs(difference[inside]) / spacing**4
    index = int(np.argmax(scaled))
    return {
        "maximum": float(scaled[index]),
        "argmax_r": float(centres[inside][index]),
        "band_lower": lower,
        "band_upper": upper,
        "dr": spacing,
        "sample_count": int(np.count_nonzero(inside)),
    }


# ------------------------------------------ preregistered acceptance measures


def initial_data_acceptance(
    grid: AxisymmetricGrid,
    *,
    amplitude_scale: float = 1.0,
    rho1: float = E32_RHO1,
    rho2: float = E32_RHO2,
    core_radius: float = E32_CORE_RADIUS,
    band: tuple[float, float] = E32_TRANSITION_BAND,
) -> dict[str, Any]:
    """Measure the six preregistered E-32 initial-data checks (prereg section 2).

    ``grid`` must be a *unit* cylinder ``r_max = 1``: the comparison partner is
    the plain E-29 datum, and E-29 is only defined on ``r <= 1``.  Extending
    ``(1-r^2)^18`` past the unit radius is precisely the error E-32 exists to
    avoid, so measuring against it there would be meaningless.

    Returned keys mirror the preregistration:

    ``core_bit_identical``
        (i) ``u~1 == u1^{E29}`` bitwise on ``r <= core_radius``.
    ``sup_deviation``
        (ii) ``sup |u~1 - u1^{E29}|`` over the whole unit cylinder, with the
        E-32c analytic bound and the preregistered ``5e-10`` acceptance.
    ``exact_zero_outside_support``
        (iii) ``u~1 == 0`` exactly on ``r >= sqrt(rho2)``.
    ``derived_norms``
        (iv) the E-29b derived norms are unchanged by the envelope: the
        envelope and the plain datum have the same ``max |u1|`` and the same
        Cartesian ``||omega(0)||_inf`` on this grid, because both maxima sit at
        ``r = 0`` and ``r = 1/sqrt(37)``, deep inside the plateau.  The
        separate grid-convergence comparison against the E-29b constants is
        reported alongside as ``e29b_*_relative_error``.
    ``fourth_difference``
        (v) the discrete ``C^4`` spot check across the transition band.
    ``axis_parity``
        (vi) the E-16c axis parity check on the enveloped datum.
    """

    if not isinstance(grid, AxisymmetricGrid):
        raise TypeError("grid must be an AxisymmetricGrid")
    if float(grid.r[-1]) != 1.0:
        raise ValueError(
            "the E-32 initial-data acceptance measures compare against E-29, "
            "which is defined only on the unit cylinder: use r_max = 1"
        )
    lower, upper = _cutoff_parameters(rho1, rho2)
    plain = hou_initial_swirl(
        grid,
        amplitude=E29_AMPLITUDE,
        amplitude_scale=amplitude_scale,
    )
    envelope = envelope_initial_swirl(
        grid,
        amplitude_scale=amplitude_scale,
        rho1=lower,
        rho2=upper,
    )
    core = core_mask(grid, core_radius=core_radius)
    support_radius = math.sqrt(upper)
    outside = grid.r >= support_radius

    deviation = float(np.max(np.abs(envelope - plain)))
    field_maximum = float(np.max(np.abs(plain)))

    plain_state = constrain_state(grid, plain, np.zeros(grid.shape))
    envelope_state = constrain_state(grid, envelope, np.zeros(grid.shape))

    def vorticity_norm(state: Any) -> float:
        return float(
            np.max(
                np.sqrt(
                    sum(
                        component * component
                        for component in cartesian_vorticity(
                            grid, state.u1, state.omega1
                        )
                    )
                )
            )
        )

    plain_u1_max = float(np.max(np.abs(plain)))
    envelope_u1_max = float(np.max(np.abs(envelope)))
    plain_vorticity = vorticity_norm(plain_state)
    envelope_vorticity = vorticity_norm(envelope_state)

    z_star = E32_Z_OF_MAX_ABS_AXIAL_FACTOR
    plain_profile = e29_radial_profile(grid.r, z_star, amplitude_scale=amplitude_scale)
    envelope_profile = envelope_radial_profile(
        grid.r,
        z_star,
        amplitude_scale=amplitude_scale,
        rho1=lower,
        rho2=upper,
    )
    plain_band = transition_band_fourth_difference(grid.r, plain_profile, band=band)
    envelope_band = transition_band_fourth_difference(
        grid.r, envelope_profile, band=band
    )

    parity = relative_axis_parity(grid, envelope)
    parity_plain = relative_axis_parity(grid, plain)

    cutoff_samples = smooth_cutoff(grid.r * grid.r, lower, upper)
    cutoff_array = np.asarray(cutoff_samples, dtype=np.float64)

    return {
        "grid": {"nr": grid.nr, "nz": grid.nz, "dr": grid.dr, "dz": grid.dz},
        "rho1": lower,
        "rho2": upper,
        "core_radius": float(core_radius),
        "support_radius": support_radius,
        # (i)
        "core_bit_identical": bool(
            np.array_equal(envelope[core], plain[core])
        ),
        "core_row_count": int(np.count_nonzero(core)),
        # (ii)
        "sup_deviation": deviation,
        "sup_deviation_relative": deviation / field_maximum,
        "sup_deviation_analytic_bound": E32_SUP_DEVIATION_BOUND,
        "sup_deviation_acceptance": E32_SUP_DEVIATION_ACCEPTANCE,
        "sup_deviation_within_analytic_bound": bool(
            deviation <= E32_SUP_DEVIATION_BOUND
        ),
        "sup_deviation_within_acceptance": bool(
            deviation <= E32_SUP_DEVIATION_ACCEPTANCE
        ),
        "sup_deviation_comparison_domain": (
            "r <= 1, the only radial range on which E-29 is defined; "
            "extending (1-r^2)^18 beyond the unit radius is the error E-32 "
            "exists to avoid"
        ),
        # (iii)
        "exact_zero_outside_support": bool(np.all(envelope[outside] == 0.0)),
        "outside_support_row_count": int(np.count_nonzero(outside)),
        # cutoff shape
        "cutoff_minimum": float(np.min(cutoff_array)),
        "cutoff_maximum": float(np.max(cutoff_array)),
        "cutoff_monotone_non_increasing": bool(
            np.all(np.diff(cutoff_array) <= 0.0)
        ),
        # (iv)
        "derived_norms": {
            "plain_max_abs_u1": plain_u1_max,
            "envelope_max_abs_u1": envelope_u1_max,
            "max_abs_u1_relative_change": abs(envelope_u1_max - plain_u1_max)
            / plain_u1_max,
            "plain_max_cartesian_vorticity": plain_vorticity,
            "envelope_max_cartesian_vorticity": envelope_vorticity,
            "max_cartesian_vorticity_relative_change": abs(
                envelope_vorticity - plain_vorticity
            )
            / plain_vorticity,
            "e29b_max_abs_u1": E32_AMPLITUDE_TIMES_MAX_ABS_AXIAL_FACTOR,
            "e29b_max_abs_u1_relative_error": abs(
                envelope_u1_max - E32_AMPLITUDE_TIMES_MAX_ABS_AXIAL_FACTOR
            )
            / E32_AMPLITUDE_TIMES_MAX_ABS_AXIAL_FACTOR,
            "interpretation": (
                "the relative CHANGE columns compare the envelope with the "
                "plain E-29 datum on this same grid and are the preregistered "
                "1e-12 check; the e29b_* column is the separate grid "
                "convergence of the derived E-29b constant, which a uniform "
                "grid of this size cannot resolve to 1e-12"
            ),
        },
        # (v)
        "fourth_difference": {
            "z": z_star,
            "z_definition": "the continuum argmax of |g(z)| recorded by E-32",
            "plain": plain_band,
            "envelope": envelope_band,
            "identical": bool(plain_band["maximum"] == envelope_band["maximum"]),
            "relative_difference": abs(
                envelope_band["maximum"] - plain_band["maximum"]
            )
            / max(plain_band["maximum"], 1.0e-300),
            "relative_tolerance": E32_FOURTH_DIFFERENCE_RELATIVE_TOLERANCE,
            "bounded": bool(
                math.isfinite(envelope_band["maximum"])
                and abs(envelope_band["maximum"] - plain_band["maximum"])
                <= E32_FOURTH_DIFFERENCE_RELATIVE_TOLERANCE
                * max(plain_band["maximum"], 1.0e-300)
            ),
            "check": (
                "preregistration (v) asks for boundedness of the high-order "
                "difference across the transition band; exact identity is the "
                "stronger observation E-32 property 7 records at dr = 1/192 "
                "and dr = 1/128 and is reported separately as identical"
            ),
            "finite": bool(math.isfinite(envelope_band["maximum"])),
        },
        # (vi)
        "axis_parity": {
            "envelope_defect": parity["defect"],
            "envelope_denominator": parity["denominator"],
            "envelope_relative": parity["relative"],
            "plain_defect": parity_plain["defect"],
            "plain_relative": parity_plain["relative"],
            "relative_change": abs(parity["relative"] - parity_plain["relative"]),
        },
    }
