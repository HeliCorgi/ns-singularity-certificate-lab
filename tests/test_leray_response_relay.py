import numpy as np
import pytest

from ns_certificate_lab.leray_response_relay import (
    adjoint_child_response,
    fejer_carrier_packet,
    harmonic_carrier_mask,
    helical_fejer_packet,
    leray_advection,
    leray_project,
    mean_energy,
    relay_stage,
    spectral_inner,
)


POLARIZATION = np.array(
    [
        0.63377081 - 0.71460668j,
        -0.12052672 - 0.13277247j,
        -0.23055982 - 0.04856754j,
    ],
    dtype=np.complex128,
)


def _frequency_index(wave: tuple[int, int, int], grid_size: int) -> tuple[int, int, int]:
    return tuple(component % grid_size for component in wave)


def _exact_relay_fourier(grid_size: int, *, include_child: bool = True) -> np.ndarray:
    field = np.zeros((3, grid_size, grid_size, grid_size), dtype=np.complex128)
    parent_sine = (1, 1, 0)
    parent_cosine = (1, 0, 1)
    child = (2, 1, 1)
    normal = np.array([1.0, -1.0, -1.0])
    sine_direction = np.array([0.0, 0.0, 1.0])
    cosine_direction = np.array([0.0, 1.0, 0.0])
    field[(slice(None), *_frequency_index(parent_sine, grid_size))] = (
        -0.5j * sine_direction
    )
    field[(slice(None), *_frequency_index(tuple(-x for x in parent_sine), grid_size))] = (
        0.5j * sine_direction
    )
    for orientation in (parent_cosine, tuple(-x for x in parent_cosine)):
        field[(slice(None), *_frequency_index(orientation, grid_size))] = (
            0.5 * cosine_direction
        )
    if include_child:
        for orientation in (child, tuple(-x for x in child)):
            field[(slice(None), *_frequency_index(orientation, grid_size))] = (
                normal / 16.0
            )
    return field


def _direct_sparse_leray_advection(field: np.ndarray) -> np.ndarray:
    """Independent coefficient convolution for a small dealiased support."""

    grid_size = field.shape[1]
    frequency = np.rint(
        np.fft.fftfreq(grid_size, d=1.0 / grid_size)
    ).astype(np.int64)
    occupied = np.argwhere(np.any(np.abs(field) > 0.0, axis=0))
    output = np.zeros_like(field, dtype=np.complex128)
    for left_index in occupied:
        left_wave = frequency[left_index]
        left = field[(slice(None), *left_index)]
        for right_index in occupied:
            right_wave = frequency[right_index]
            right = field[(slice(None), *right_index)]
            target_wave = left_wave + right_wave
            coefficient = 1.0j * np.dot(left, right_wave) * right
            norm_squared = int(np.dot(target_wave, target_wave))
            if norm_squared == 0:
                projected = np.zeros(3, dtype=np.complex128)
            else:
                projected = coefficient - target_wave * (
                    np.dot(target_wave, coefficient) / norm_squared
                )
            target = tuple(int(component) % grid_size for component in target_wave)
            output[(slice(None), *target)] += projected
    return output


def _divergence_l2(field: np.ndarray) -> float:
    grid_size = field.shape[1]
    frequency = np.fft.fftfreq(grid_size, d=1.0 / grid_size)
    wave = np.stack(
        np.meshgrid(frequency, frequency, frequency, indexing="ij"), axis=0
    )
    return float(np.linalg.norm(1.0j * np.sum(wave * field, axis=0)))


def test_fejer_best_response_has_true_positive_leray_flux() -> None:
    parent = fejer_carrier_packet(
        48, carrier=4, envelope=2, polarization=POLARIZATION
    )
    mask = harmonic_carrier_mask(48, carrier=4, envelope=2, harmonic=2)
    stage = relay_stage(parent, mask, parent_scale=4.0, viscosity=1.0)
    metrics = stage.metrics
    physical_parent = np.fft.ifftn(parent, axes=(1, 2, 3)) * 48**3
    assert mean_energy(parent) == pytest.approx(0.25, rel=1.0e-13)
    assert np.max(np.abs(physical_parent.imag)) < 1.0e-11
    assert metrics["parent_parent_injection"] > 0.0
    assert metrics["best_response_alignment"] == pytest.approx(1.0, abs=1.0e-12)
    assert metrics["parent_divergence_l2"] < 1.0e-12
    assert metrics["child_divergence_l2"] < 1.0e-12
    assert metrics["parent_energy_cancellation_defect"] < 1.0e-12
    assert metrics["full_populated_flux_tested"]
    assert metrics["full_populated_child_flux"] > 0.0
    assert metrics["full_nonlinear_energy_defect"] < 1.0e-12
    assert metrics["forcing_to_resolution_floor"] > 1.0e8
    assert (
        metrics["parent_support_min_wavenumber"]
        <= metrics["parent_scale"]
        <= metrics["parent_support_max_wavenumber"]
    )


def test_quadratic_chirp_changes_low_high_autocorrelation_balance() -> None:
    mask = harmonic_carrier_mask(48, carrier=4, envelope=2, harmonic=2)
    plain = relay_stage(
        fejer_carrier_packet(
            48, carrier=4, envelope=2, polarization=POLARIZATION, chirp=0.0
        ),
        mask,
        parent_scale=4.0,
        viscosity=1.0,
    )
    chirped = relay_stage(
        fejer_carrier_packet(
            48, carrier=4, envelope=2, polarization=POLARIZATION, chirp=1.5
        ),
        mask,
        parent_scale=4.0,
        viscosity=1.0,
    )
    assert (
        chirped.metrics["off_chain_forcing_ratio"]
        < plain.metrics["off_chain_forcing_ratio"]
    )


def test_same_helicity_packet_is_real_and_divergence_free() -> None:
    packet = helical_fejer_packet(
        32, carrier=4, envelope=2, helicity=1, chirp=0.75
    )
    physical = np.fft.ifftn(packet, axes=(1, 2, 3)) * 32**3
    assert np.max(np.abs(physical.imag)) < 1.0e-11
    assert _divergence_l2(packet) < 1.0e-12
    frequency = np.fft.fftfreq(32, d=1.0 / 32)
    wave = np.stack(
        np.meshgrid(frequency, frequency, frequency, indexing="ij"), axis=0
    )
    curl = 1.0j * np.moveaxis(
        np.cross(np.moveaxis(wave, 0, -1), np.moveaxis(packet, 0, -1)),
        -1,
        0,
    )
    expected = np.sqrt(np.sum(wave * wave, axis=0))[None, ...] * packet
    assert np.max(np.abs(curl - expected)) < 1.0e-12
    nonlinear = leray_advection(packet, packet)
    assert np.all(np.isfinite(nonlinear))


def test_fft_advection_matches_independent_sparse_convolution_and_sign() -> None:
    field = _exact_relay_fourier(32)
    computed = leray_advection(field, field)
    independent = _direct_sparse_leray_advection(field)
    assert np.max(np.abs(computed - independent)) < 2.0e-15

    child_wave = _frequency_index((2, 1, 1), 32)
    assert computed[(slice(None), *child_wave)] == pytest.approx(
        np.array([-1.0, 1.0, 1.0]) / 6.0, abs=2.0e-15
    )
    child = np.zeros_like(field)
    for orientation in ((2, 1, 1), (-2, -1, -1)):
        index = _frequency_index(orientation, 32)
        child[(slice(None), *index)] = field[(slice(None), *index)]
    assert -spectral_inner(child, computed) == pytest.approx(1.0 / 16.0)


def test_real_fourier_dtype_is_promoted_before_complex_convolution() -> None:
    grid_size = 16
    field = np.zeros((3, grid_size, grid_size, grid_size), dtype=np.float64)
    for orientation in (1, -1):
        field[1, orientation % grid_size, 0, 0] = 0.5
        field[0, 0, orientation % grid_size, orientation % grid_size] = 0.5
    promoted = leray_advection(field, field)
    reference = leray_advection(
        field.astype(np.complex128), field.astype(np.complex128)
    )
    assert np.linalg.norm(promoted) == pytest.approx(1.0 / np.sqrt(6.0))
    assert promoted == pytest.approx(reference, abs=2.0e-15)


def test_nonphysical_parent_and_child_masks_are_rejected() -> None:
    parent = fejer_carrier_packet(
        48, carrier=4, envelope=2, polarization=POLARIZATION
    )
    mask = harmonic_carrier_mask(48, carrier=4, envelope=2, harmonic=2)
    frequency = np.fft.fftfreq(48, d=1.0 / 48)
    kx = np.meshgrid(frequency, frequency, frequency, indexing="ij")[0]
    with pytest.raises(ValueError, match="symmetric"):
        adjoint_child_response(
            parent, mask & (kx > 0.0), target_energy=0.125
        )

    compressible = parent.copy()
    compressible[0, 4, 0, 0] += 0.1
    compressible[0, -4 % 48, 0, 0] += 0.1
    with pytest.raises(ValueError, match="divergence free"):
        relay_stage(
            compressible, mask, parent_scale=4.0, viscosity=1.0
        )

    overlapping_parent = fejer_carrier_packet(
        32, carrier=3, envelope=2, polarization=POLARIZATION
    )
    overlapping_mask = harmonic_carrier_mask(
        32, carrier=3, envelope=2, harmonic=2
    )
    with pytest.raises(ValueError, match="disjoint"):
        relay_stage(
            overlapping_parent,
            overlapping_mask,
            parent_scale=3.0,
            viscosity=1.0,
        )


def test_unreachable_and_roundoff_only_response_bands_are_rejected() -> None:
    parent = fejer_carrier_packet(
        64, carrier=4, envelope=2, polarization=POLARIZATION
    )
    unreachable = np.zeros((64, 64, 64), dtype=np.bool_)
    unreachable[15, 0, 0] = True
    unreachable[-15 % 64, 0, 0] = True
    with pytest.raises(ValueError, match="does not intersect"):
        adjoint_child_response(parent, unreachable, target_energy=0.125)

    cancelling_parent = _exact_relay_fourier(32, include_child=False)
    cancelling_band = np.zeros((32, 32, 32), dtype=np.bool_)
    difference = (0, 1, -1)
    cancelling_band[_frequency_index(difference, 32)] = True
    cancelling_band[_frequency_index(tuple(-x for x in difference), 32)] = True
    with pytest.raises(ValueError, match="roundoff floor"):
        adjoint_child_response(
            cancelling_parent, cancelling_band, target_energy=0.25
        )


def test_mean_mode_and_inconsistent_parent_scale_are_rejected() -> None:
    mean = np.zeros((3, 8, 8, 8), dtype=np.complex128)
    mean[:, 0, 0, 0] = (1.0, 2.0, 3.0)
    with pytest.raises(ValueError, match="mean zero"):
        leray_project(mean)

    parent = fejer_carrier_packet(
        48, carrier=4, envelope=2, polarization=POLARIZATION
    )
    mask = harmonic_carrier_mask(48, carrier=4, envelope=2, harmonic=2)
    with pytest.raises(ValueError, match="resolved parent support"):
        relay_stage(parent, mask, parent_scale=40.0, viscosity=1.0)


def test_aliasing_and_invalid_masks_are_rejected() -> None:
    packet = fejer_carrier_packet(
        32, carrier=4, envelope=2, polarization=POLARIZATION
    )
    mask = harmonic_carrier_mask(32, carrier=4, envelope=2, harmonic=2)
    stage = relay_stage(
        packet, mask, parent_scale=4.0, viscosity=1.0
    )
    with pytest.raises(ValueError, match="not dealiased"):
        leray_advection(stage.child, stage.child)
    with pytest.raises(ValueError):
        relay_stage(
            packet,
            np.zeros((32, 32, 32), dtype=np.bool_),
            parent_scale=6.0,
            viscosity=1.0,
        )
