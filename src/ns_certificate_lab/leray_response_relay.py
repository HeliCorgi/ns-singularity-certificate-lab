r"""Adjoint-seeded response maps for a true periodic Leray interaction.

This module studies a deliberately narrow discovery question.  Given a
divergence-free parent Fourier packet ``p`` and a disjoint child band ``C``,
define

.. math::

   g_C(p)=-P_C\mathbb P((p\cdot\nabla)p),\qquad
   c=\sqrt{2E_C}\,g_C/\|g_C\|_2.

Then ``<c, g_C> > 0`` is an algebraic identity: the child is the Riesz best
response to the *actual* Navier--Stokes/Leray quadratic term.  This does not
establish a cascade.  Viscosity, the other output bands, interactions with an
already populated child, persistence over an interval, and iteration of the
response map are separate tests and are all reported explicitly.

Fourier coefficients use the normalized torus convention
``u(x)=sum_k u_hat[k] exp(i k.x)``.  Hence the spatial mean of ``|u|^2`` is
``sum_k |u_hat[k]|^2``.  Products are evaluated on a grid only after an exact
support check guarantees that no Fourier aliasing is possible.

This discovery lane is restricted to real, mean-zero, divergence-free fields.
Reality is represented by Hermitian symmetry, ``u_hat[-k]=conj(u_hat[k])``.
The zero Fourier coefficient is therefore rejected at the response-map
boundary, and :func:`leray_project` sets the zero coefficient to zero after
checking that it is only roundoff.  This is a mean-zero convention, not the
definition of the Leray projector on a torus with arbitrary constant fields.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


ComplexArray = npt.NDArray[np.complex128]
BoolArray = npt.NDArray[np.bool_]
FloatArray = npt.NDArray[np.float64]


_STRUCTURE_RELATIVE_TOLERANCE = 1.0e-11
_SUPPORT_RELATIVE_TOLERANCE = 1024.0 * np.finfo(np.float64).eps
_FORCING_ROUNDOFF_FACTOR = 4096.0

__all__ = [
    "RelayStage",
    "adjoint_child_response",
    "fejer_carrier_packet",
    "gradient_l2_squared",
    "harmonic_carrier_mask",
    "helical_fejer_packet",
    "leray_advection",
    "leray_project",
    "mean_energy",
    "relay_stage",
    "spectral_inner",
]


def _grid_size(field: ComplexArray) -> int:
    value = np.asarray(field)
    if (
        value.ndim != 4
        or value.shape[0] != 3
        or value.shape[1] != value.shape[2]
        or value.shape[1] != value.shape[3]
        or value.shape[1] < 8
    ):
        raise ValueError("a Fourier field must have shape (3,n,n,n), n >= 8")
    if not np.all(np.isfinite(value)):
        raise ValueError("Fourier coefficients must be finite")
    return int(value.shape[1])


def _frequencies(grid_size: int) -> FloatArray:
    if isinstance(grid_size, bool) or not isinstance(grid_size, int) or grid_size < 8:
        raise ValueError("grid_size must be an integer at least 8")
    return np.asarray(np.fft.fftfreq(grid_size, d=1.0 / grid_size), dtype=np.float64)


def _frequency_mesh(grid_size: int) -> tuple[FloatArray, FloatArray, FloatArray]:
    frequency = _frequencies(grid_size)
    return tuple(np.meshgrid(frequency, frequency, frequency, indexing="ij"))  # type: ignore[return-value]


def spectral_inner(left: ComplexArray, right: ComplexArray) -> float:
    """Return the real normalized-torus ``L2`` inner product."""

    n = _grid_size(left)
    if _grid_size(right) != n:
        raise ValueError("Fourier fields must use the same grid")
    value = float(np.vdot(left, right).real)
    if not math.isfinite(value):
        raise ValueError("the spectral inner product overflowed")
    return value


def mean_energy(field: ComplexArray) -> float:
    """Return ``(1/2) mean_x |u|^2`` by Parseval."""

    return 0.5 * spectral_inner(field, field)


def gradient_l2_squared(field: ComplexArray) -> float:
    """Return ``mean_x |grad u|^2`` by Parseval."""

    n = _grid_size(field)
    kx, ky, kz = _frequency_mesh(n)
    weight = kx * kx + ky * ky + kz * kz
    value = float(np.sum(weight[None, ...] * np.abs(field) ** 2).real)
    if not math.isfinite(value):
        raise ValueError("the gradient norm overflowed")
    return value


def _support_limits(field: ComplexArray) -> tuple[int, int, int]:
    n = _grid_size(field)
    frequency = np.rint(_frequencies(n)).astype(np.int64)
    occupied = np.any(np.abs(field) > 0.0, axis=0)
    if not np.any(occupied):
        return (0, 0, 0)
    limits = []
    for axis in range(3):
        shape = [1, 1, 1]
        shape[axis] = n
        values = np.broadcast_to(np.abs(frequency).reshape(shape), occupied.shape)
        limits.append(int(np.max(values[occupied])))
    return tuple(limits)  # type: ignore[return-value]


def _require_dealiased(left: ComplexArray, right: ComplexArray) -> None:
    n = _grid_size(left)
    if _grid_size(right) != n:
        raise ValueError("Fourier fields must use the same grid")
    left_limits = _support_limits(left)
    right_limits = _support_limits(right)
    nyquist = n // 2
    for axis, (left_limit, right_limit) in enumerate(zip(left_limits, right_limits)):
        if left_limit + right_limit >= nyquist:
            raise ValueError(
                "Fourier product is not dealiased on axis "
                f"{axis}: {left_limit}+{right_limit} >= {nyquist}"
            )


def leray_project(field: ComplexArray) -> ComplexArray:
    """Apply ``I-k k^T/|k|^2`` on the mean-zero Fourier subspace.

    A resolved nonzero mean is rejected.  The zero coefficient is then set to
    zero so that roundoff from FFT products cannot leak out of the mean-zero
    research lane.  Constant vector fields require a different convention.
    """

    n = _grid_size(field)
    value = np.asarray(field, dtype=np.complex128)
    field_norm = float(np.linalg.norm(value))
    mean_norm = float(np.linalg.norm(value[:, 0, 0, 0]))
    if mean_norm > _STRUCTURE_RELATIVE_TOLERANCE * max(
        field_norm, np.finfo(np.float64).tiny
    ):
        raise ValueError("Leray response fields must be mean zero")
    kx, ky, kz = _frequency_mesh(n)
    wave = np.stack((kx, ky, kz), axis=0)
    norm_squared = kx * kx + ky * ky + kz * kz
    dot = np.sum(wave * value, axis=0)
    denominator = np.where(norm_squared > 0.0, norm_squared, 1.0)
    projected = np.asarray(
        value - wave * (dot / denominator)[None, ...], dtype=np.complex128
    )
    projected[:, 0, 0, 0] = 0.0
    return projected


def leray_advection(advecting: ComplexArray, advected: ComplexArray) -> ComplexArray:
    r"""Return ``P((advecting . grad) advected)`` without Fourier aliasing."""

    n = _grid_size(advecting)
    if _grid_size(advected) != n:
        raise ValueError("Fourier fields must use the same grid")
    left = np.asarray(advecting, dtype=np.complex128)
    right = np.asarray(advected, dtype=np.complex128)
    _require_dealiased(left, right)
    axes = (1, 2, 3)
    scale = float(n**3)
    frequencies = _frequency_mesh(n)
    physical_advecting = np.fft.ifftn(left, axes=axes) * scale
    output = np.empty(right.shape, dtype=np.complex128)
    for component in range(3):
        accumulated = np.zeros((n, n, n), dtype=np.complex128)
        for direction in range(3):
            derivative_hat = (
                1.0j * frequencies[direction] * right[component]
            )
            derivative = np.fft.ifftn(derivative_hat, axes=(0, 1, 2)) * scale
            accumulated += physical_advecting[direction] * derivative
        output[component] = np.fft.fftn(accumulated, axes=(0, 1, 2)) / scale
    return leray_project(output)


def _negative_frequency_reflection(field: ComplexArray) -> ComplexArray:
    n = _grid_size(field)
    indices = (-np.arange(n, dtype=np.int64)) % n
    reflected = np.take(field, indices, axis=1)
    reflected = np.take(reflected, indices, axis=2)
    return np.take(reflected, indices, axis=3)


def _significant_support(field: ComplexArray) -> BoolArray:
    value = np.asarray(field, dtype=np.complex128)
    _grid_size(value)
    amplitude = np.max(np.abs(value), axis=0)
    peak = float(np.max(amplitude))
    if not peak > 0.0:
        return np.zeros(amplitude.shape, dtype=np.bool_)
    return np.asarray(
        amplitude > _SUPPORT_RELATIVE_TOLERANCE * peak, dtype=np.bool_
    )


def _quadratic_support(support: BoolArray) -> BoolArray:
    """Return the cyclic sumset of a scalar Fourier support mask."""

    mask = np.asarray(support, dtype=np.bool_)
    if mask.ndim != 3 or not (mask.shape[0] == mask.shape[1] == mask.shape[2]):
        raise ValueError("support must be a cubic scalar Fourier mask")
    coordinates = np.argwhere(mask)
    reachable = np.zeros(mask.shape, dtype=np.bool_)
    if coordinates.size == 0:
        return reachable
    # Sparse response packets are much cheaper to sum directly.  The FFT path
    # prevents a quadratic blow-up if this diagnostic is used on a dense cloud.
    if coordinates.shape[0] <= 4096:
        n = mask.shape[0]
        for left in coordinates:
            sums = (coordinates + left[None, :]) % n
            reachable[sums[:, 0], sums[:, 1], sums[:, 2]] = True
        return reachable
    counts = np.fft.ifftn(np.fft.fftn(mask.astype(np.float64)) ** 2).real
    return np.asarray(counts > 0.5, dtype=np.bool_)


def _support_wavenumber_range(field: ComplexArray) -> tuple[float, float, float]:
    n = _grid_size(field)
    support = _significant_support(field)
    if not np.any(support):
        raise ValueError("Fourier field has no resolved nonzero support")
    kx, ky, kz = _frequency_mesh(n)
    radius = np.sqrt(kx * kx + ky * ky + kz * kz)
    occupied_radius = radius[support]
    if np.any(occupied_radius == 0.0):
        nonzero = occupied_radius[occupied_radius > 0.0]
    else:
        nonzero = occupied_radius
    if nonzero.size == 0:
        raise ValueError("Fourier field has only a mean mode")
    l2_squared = spectral_inner(field, field)
    rms = math.sqrt(gradient_l2_squared(field) / l2_squared)
    return float(np.min(nonzero)), float(np.max(nonzero)), rms


def _require_physical_parent(parent: ComplexArray) -> ComplexArray:
    value = np.asarray(parent, dtype=np.complex128)
    _grid_size(value)
    field_norm = float(np.linalg.norm(value))
    if not field_norm > 0.0:
        raise ValueError("parent must be a nonzero Fourier field")
    reflected = _negative_frequency_reflection(value)
    reality_defect = float(np.linalg.norm(value - np.conjugate(reflected)))
    if reality_defect > _STRUCTURE_RELATIVE_TOLERANCE * field_norm:
        raise ValueError("parent Fourier coefficients are not Hermitian symmetric")
    mean_norm = float(np.linalg.norm(value[:, 0, 0, 0]))
    if mean_norm > _STRUCTURE_RELATIVE_TOLERANCE * field_norm:
        raise ValueError("parent must be mean zero")
    value = value.copy()
    value[:, 0, 0, 0] = 0.0
    gradient_norm = math.sqrt(gradient_l2_squared(value))
    divergence = _divergence_l2(value)
    if divergence > _STRUCTURE_RELATIVE_TOLERANCE * max(
        gradient_norm, np.finfo(np.float64).tiny
    ):
        raise ValueError("parent must be divergence free")
    return value


def _validated_child_mask(parent: ComplexArray, child_mask: BoolArray) -> BoolArray:
    n = _grid_size(parent)
    mask = np.asarray(child_mask, dtype=np.bool_)
    if mask.shape != (n, n, n) or not np.any(mask):
        raise ValueError("child_mask must be a nonempty scalar Fourier mask")
    indices = (-np.arange(n, dtype=np.int64)) % n
    reflected = np.take(mask, indices, axis=0)
    reflected = np.take(reflected, indices, axis=1)
    reflected = np.take(reflected, indices, axis=2)
    if not np.array_equal(mask, reflected):
        raise ValueError("child_mask must be symmetric under k -> -k")
    support = _significant_support(parent)
    if np.any(mask & support):
        raise ValueError("child_mask must be disjoint from the parent support")
    reachable = _quadratic_support(support)
    active = np.asarray(mask & reachable, dtype=np.bool_)
    if not np.any(active):
        raise ValueError("child_mask does not intersect the quadratic parent support")
    return active


def _forcing_resolution_floor(parent: ComplexArray) -> float:
    n = _grid_size(parent)
    _, maximum, _ = _support_wavenumber_range(parent)
    parent_l2_squared = spectral_inner(parent, parent)
    return (
        _FORCING_ROUNDOFF_FACTOR
        * np.finfo(np.float64).eps
        * max(1.0, math.log2(float(n)))
        * maximum
        * parent_l2_squared
    )


def _energy_rescaled(field: ComplexArray, target_energy: float) -> ComplexArray:
    if not math.isfinite(target_energy) or target_energy <= 0.0:
        raise ValueError("target_energy must be finite and positive")
    energy = mean_energy(field)
    if not energy > 0.0:
        raise ValueError("cannot normalize a zero Fourier field")
    return np.asarray(field * math.sqrt(target_energy / energy), dtype=np.complex128)


def fejer_carrier_packet(
    grid_size: int,
    *,
    carrier: int,
    envelope: int,
    polarization: npt.ArrayLike,
    chirp: float = 0.0,
    energy_constant: float = 1.0,
) -> ComplexArray:
    r"""Construct a divergence-free, Fourier-localized Fejer carrier packet.

    A vector potential has coefficients

    ``A_hat(q +/- K e1) = (1/2) prod_i(1-|q_i|/m) (z or conj(z))``

    for ``|q_i| < m``.  The velocity is ``curl A``.  Its energy is normalized
    to ``energy_constant / carrier``, the critical packet energy law.
    ``polarization`` is a complex three-vector and supplies both cosine and
    sine phases while conjugate symmetry keeps the physical field real.  A
    nonzero ``chirp`` multiplies the positive carrier coefficient by
    ``exp(i chirp |q|^2)``; this changes aperiodic autocorrelation and is not a
    translation.
    """

    n = int(grid_size)
    _frequencies(n)
    if isinstance(carrier, bool) or not isinstance(carrier, int) or carrier < 2:
        raise ValueError("carrier must be an integer at least 2")
    if isinstance(envelope, bool) or not isinstance(envelope, int) or envelope < 1:
        raise ValueError("envelope must be a positive integer")
    if carrier <= 2 * (envelope - 1):
        raise ValueError("the positive and negative carrier boxes must be disjoint")
    if 2 * (carrier + envelope - 1) >= n:
        raise ValueError("the parent carrier does not fit below the Nyquist mode")
    if not math.isfinite(energy_constant) or energy_constant <= 0.0:
        raise ValueError("energy_constant must be finite and positive")
    if not math.isfinite(chirp):
        raise ValueError("chirp must be finite")
    z = np.asarray(polarization, dtype=np.complex128)
    if z.shape != (3,) or not np.all(np.isfinite(z)) or np.linalg.norm(z) == 0.0:
        raise ValueError("polarization must be a finite nonzero complex 3-vector")
    z = z / np.linalg.norm(z)

    potential = np.zeros((3, n, n, n), dtype=np.complex128)
    support = range(-(envelope - 1), envelope)
    for qx in support:
        wx = 1.0 - abs(qx) / envelope
        for qy in support:
            wy = 1.0 - abs(qy) / envelope
            for qz in support:
                wz = 1.0 - abs(qz) / envelope
                phase = np.exp(1.0j * chirp * float(qx * qx + qy * qy + qz * qz))
                coefficient = 0.5 * wx * wy * wz
                positive = ((carrier + qx) % n, qy % n, qz % n)
                negative = ((-carrier + qx) % n, qy % n, qz % n)
                potential[(slice(None), *positive)] += coefficient * phase * z
                potential[(slice(None), *negative)] += coefficient * np.conjugate(
                    phase * z
                )

    kx, ky, kz = _frequency_mesh(n)
    wave = np.stack((kx, ky, kz), axis=0)
    velocity = 1.0j * np.cross(
        np.moveaxis(wave, 0, -1), np.moveaxis(potential, 0, -1)
    )
    velocity = np.moveaxis(velocity, -1, 0)
    return _energy_rescaled(velocity, energy_constant / carrier)


def helical_fejer_packet(
    grid_size: int,
    *,
    carrier: int,
    envelope: int,
    helicity: int = 1,
    chirp: float = 0.0,
    energy_constant: float = 1.0,
) -> ComplexArray:
    r"""Construct a same-helicity Fejer packet with an optional quadratic chirp.

    For every mode ``k=K e1+q`` in the positive carrier box, the coefficient
    is a Fejer weight times the normalized helical vector
    ``h_s(k)=(e_perp+i s khat cross e_perp)/sqrt(2)``.  Its negative mode is
    set by conjugate symmetry.  Equal-radius, equal-helicity modes would form
    a Beltrami field and have zero projected self-interaction; the envelope's
    radial spread deliberately breaks that cancellation.  ``chirp`` inserts
    ``exp(i chirp |q|^2)`` and is not a mere translation when nonzero.
    """

    n = int(grid_size)
    _frequencies(n)
    if isinstance(carrier, bool) or not isinstance(carrier, int) or carrier < 2:
        raise ValueError("carrier must be an integer at least 2")
    if isinstance(envelope, bool) or not isinstance(envelope, int) or envelope < 1:
        raise ValueError("envelope must be a positive integer")
    if carrier <= 2 * (envelope - 1):
        raise ValueError("the positive and negative carrier boxes must be disjoint")
    if 2 * (carrier + envelope - 1) >= n:
        raise ValueError("the parent carrier does not fit below the Nyquist mode")
    if helicity not in (-1, 1):
        raise ValueError("helicity must be -1 or +1")
    if not math.isfinite(chirp):
        raise ValueError("chirp must be finite")
    if not math.isfinite(energy_constant) or energy_constant <= 0.0:
        raise ValueError("energy_constant must be finite and positive")

    velocity = np.zeros((3, n, n, n), dtype=np.complex128)
    support = range(-(envelope - 1), envelope)
    reference = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    for qx in support:
        wx = 1.0 - abs(qx) / envelope
        for qy in support:
            wy = 1.0 - abs(qy) / envelope
            for qz in support:
                wz = 1.0 - abs(qz) / envelope
                wave = np.array([carrier + qx, qy, qz], dtype=np.float64)
                wave_norm = float(np.linalg.norm(wave))
                unit = wave / wave_norm
                perpendicular = reference - float(np.dot(reference, unit)) * unit
                perpendicular /= np.linalg.norm(perpendicular)
                second = np.cross(unit, perpendicular)
                helical = (
                    perpendicular + 1.0j * float(helicity) * second
                ) / math.sqrt(2.0)
                phase = np.exp(1.0j * chirp * float(qx * qx + qy * qy + qz * qz))
                coefficient = wx * wy * wz * phase
                positive = ((carrier + qx) % n, qy % n, qz % n)
                negative = ((-carrier - qx) % n, (-qy) % n, (-qz) % n)
                velocity[(slice(None), *positive)] += coefficient * helical
                velocity[(slice(None), *negative)] += np.conjugate(
                    coefficient * helical
                )
    return _energy_rescaled(velocity, energy_constant / carrier)


def harmonic_carrier_mask(
    grid_size: int,
    *,
    carrier: int,
    envelope: int,
    harmonic: int,
) -> BoolArray:
    """Select the two boxes around ``+/- harmonic * carrier * e1``."""

    n = int(grid_size)
    kx, ky, kz = _frequency_mesh(n)
    if isinstance(harmonic, bool) or not isinstance(harmonic, int) or harmonic < 1:
        raise ValueError("harmonic must be a positive integer")
    width = harmonic * (envelope - 1)
    centre = harmonic * carrier
    if centre + width >= n // 2:
        raise ValueError("harmonic carrier box reaches the Nyquist mode")
    return np.asarray(
        (np.abs(np.abs(kx) - centre) <= width)
        & (np.abs(ky) <= width)
        & (np.abs(kz) <= width),
        dtype=np.bool_,
    )


def adjoint_child_response(
    parent: ComplexArray,
    child_mask: BoolArray,
    *,
    target_energy: float,
) -> tuple[ComplexArray, ComplexArray, ComplexArray]:
    r"""Return the normalized best-response child, forcing, and ``B(p,p)``."""

    value = _require_physical_parent(parent)
    _require_dealiased(value, value)
    mask = _validated_child_mask(value, child_mask)
    nonlinear = leray_advection(value, value)
    forcing = np.asarray(-nonlinear * mask[None, ...], dtype=np.complex128)
    forcing_norm = math.sqrt(spectral_inner(forcing, forcing))
    resolution_floor = _forcing_resolution_floor(value)
    if forcing_norm <= resolution_floor:
        raise ValueError(
            "child forcing is below the resolved Fourier roundoff floor"
        )
    child = _energy_rescaled(forcing, target_energy)
    return child, forcing, nonlinear


@dataclass(frozen=True)
class RelayStage:
    """One algebraically positive parent-to-child response stage."""

    parent: ComplexArray
    child: ComplexArray
    parent_nonlinear: ComplexArray
    metrics: dict[str, float | int | bool | None]


def relay_stage(
    parent: ComplexArray,
    child_mask: BoolArray,
    *,
    parent_scale: float,
    viscosity: float,
) -> RelayStage:
    """Build and diagnose one adjoint-seeded Leray relay stage."""

    if (
        isinstance(parent_scale, bool)
        or not math.isfinite(parent_scale)
        or parent_scale <= 0.0
    ):
        raise ValueError("parent_scale must be finite and positive")
    if (
        isinstance(viscosity, bool)
        or not math.isfinite(viscosity)
        or viscosity <= 0.0
    ):
        raise ValueError("viscosity must be finite and positive")
    parent = _require_physical_parent(parent)
    minimum_scale, maximum_scale, rms_scale = _support_wavenumber_range(parent)
    scale_tolerance = _STRUCTURE_RELATIVE_TOLERANCE * max(1.0, maximum_scale)
    if not (
        minimum_scale - scale_tolerance
        <= parent_scale
        <= maximum_scale + scale_tolerance
    ):
        raise ValueError("parent_scale must lie inside the resolved parent support")
    parent_energy = mean_energy(parent)
    child, forcing, parent_nonlinear = adjoint_child_response(
        parent, child_mask, target_energy=0.5 * parent_energy
    )
    child_energy = mean_energy(child)
    injection = spectral_inner(child, forcing)
    forcing_sq = spectral_inner(forcing, forcing)
    forcing_norm = math.sqrt(forcing_sq)
    forcing_resolution_floor = _forcing_resolution_floor(parent)
    nonlinear_sq = spectral_inner(parent_nonlinear, parent_nonlinear)
    active_child_mask = _significant_support(child)
    off_chain = np.asarray(
        parent_nonlinear * (~active_child_mask)[None, ...],
        dtype=np.complex128,
    )
    off_chain_sq = spectral_inner(off_chain, off_chain)
    child_dissipation_shape = gradient_l2_squared(child)
    energy_constant = parent_scale * parent_energy
    injection_shape = injection / energy_constant**1.5
    dissipation_shape = child_dissipation_shape / energy_constant
    critical_energy_constant = (
        viscosity * dissipation_shape / injection_shape
    ) ** 2

    parent_defect = abs(spectral_inner(parent, parent_nonlinear))
    child_divergence = _divergence_l2(child)
    parent_divergence = _divergence_l2(parent)
    alignment = injection / math.sqrt(
        max(spectral_inner(child, child) * forcing_sq, np.finfo(float).tiny)
    )
    fill_time = child_energy / injection

    full_flux: float | None = None
    full_parent_flux: float | None = None
    full_energy_defect: float | None = None
    full_critical_energy_constant: float | None = None
    try:
        total = parent + child
        total_nonlinear = leray_advection(total, total)
        full_flux = -spectral_inner(child, total_nonlinear)
        full_parent_flux = -spectral_inner(parent, total_nonlinear)
        full_energy_defect = abs(spectral_inner(total, total_nonlinear))
        if full_flux > 0.0:
            full_shape = full_flux / energy_constant**1.5
            full_critical_energy_constant = (
                viscosity * dissipation_shape / full_shape
            ) ** 2
    except ValueError as error:
        if "not dealiased" not in str(error):
            raise

    metrics: dict[str, float | int | bool | None] = {
        "grid_size": _grid_size(parent),
        "parent_scale": float(parent_scale),
        "parent_support_min_wavenumber": minimum_scale,
        "parent_support_max_wavenumber": maximum_scale,
        "parent_rms_wavenumber": rms_scale,
        "parent_scale_to_rms_ratio": parent_scale / rms_scale,
        "parent_energy": parent_energy,
        "child_energy": child_energy,
        "energy_constant": energy_constant,
        "parent_parent_injection": injection,
        "parent_parent_forcing_l2": forcing_norm,
        "forcing_resolution_floor": forcing_resolution_floor,
        "forcing_to_resolution_floor": forcing_norm / forcing_resolution_floor,
        "normalized_injection_per_scale": injection / parent_scale,
        "child_gradient_l2_squared": child_dissipation_shape,
        "normalized_dissipation_per_scale": child_dissipation_shape / parent_scale,
        "critical_energy_constant": critical_energy_constant,
        "child_forcing_fraction": forcing_sq / nonlinear_sq,
        "off_chain_forcing_ratio": math.sqrt(off_chain_sq / forcing_sq),
        "parabolic_fill_constant": fill_time * parent_scale**2,
        "best_response_alignment": alignment,
        "parent_energy_cancellation_defect": parent_defect,
        "parent_divergence_l2": parent_divergence,
        "child_divergence_l2": child_divergence,
        "full_populated_child_flux": full_flux,
        "full_populated_parent_flux": full_parent_flux,
        "full_nonlinear_energy_defect": full_energy_defect,
        "full_critical_energy_constant": full_critical_energy_constant,
        "full_populated_flux_tested": full_flux is not None,
    }
    return RelayStage(
        parent=np.asarray(parent, dtype=np.complex128),
        child=child,
        parent_nonlinear=parent_nonlinear,
        metrics=metrics,
    )


def _divergence_l2(field: ComplexArray) -> float:
    n = _grid_size(field)
    wave = np.stack(_frequency_mesh(n), axis=0)
    divergence = 1.0j * np.sum(wave * field, axis=0)
    return float(np.sqrt(np.sum(np.abs(divergence) ** 2)))
