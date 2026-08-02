"""Exact rational certificates for the coherent critical-spectrum family.

Supports Lemma L* of the final theorem document.  No new quantities are
introduced: everything below is ``H_0``, ``H_1``, ``K``, and the already
defined sweeping pairing ``Q(c)=<(u.grad)u,(c.grad)u>`` evaluated on the
family already fixed by ``experiments/run_osgood_gate.py``
(``coherent_field`` / ``exact_anchor``),

    u(x) = sum_{1<=|k|<=N} a_k e^{i k.x},   a_k = P_k v_0 / |k|^2,
    v_0 = (1,2,3),   P_k = I - k k^T/|k|^2,   a_{-k} = a_k  (real cosine field).

Conventions (identical to ``exact_anchor``): the pair energy of a mode is its
``L^2`` *mean*, so with ``a_k`` stored at every full-band lattice point

    H_0 = sum_{full band} |a_k|^2 = ||u||_2^2,
    H_1 = sum_{full band} |k|^2 |a_k|^2 = ||grad u||_2^2 = D,
    K   = ||P(u.grad u)||_2^2 / ||grad u||_2^4 = ||N||_2^2 / H_1^2.

Exact lane.  All rationals are carried as Python integers over one common
denominator ``M = lcm{|k|^2}^2`` (so ``a_k = A_k / M`` with ``A_k`` an integer
3-vector), and converted to ``fractions.Fraction`` only for the final
contractions.  The rational value is therefore *exact*, not interval.

Pre-registered enumeration bounds and compute caps (fixed before the run, see
``preregistration`` in the emitted JSON):

  * exact symmetry identities  H_0 = (2/3)|v_0|^2 T_N, H_1 = (2/3)|v_0|^2 S_N
    at N in {4, 6, 8}; enumeration = every lattice point with 1<=|k|^2<=N^2,
    no truncation, no sampling;
  * exact K, Q(v_0), sweeping norm at N in {4, 6, 8}; pair enumeration =
    |band|^2 ordered pairs (N=4: 65_536, N=6: 853_776, N=8: 4_443_664),
    no truncation;
  * exact-lane wall cap 20 min, float-lane wall cap 10 min, total cap 30 min;
  * float continuation at N in {8, 16, 32} on dealias-safe grids
    (48, 48, 80, 144, 144); no other N is examined.

Float lane reuses ``coherent_field`` and ``measure`` from
``experiments.run_osgood_gate`` unmodified.  Every emitted number is labelled
``exact`` or ``float``.  Diagnostic computation; no PDE theorem is asserted
here.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from ns_certificate_lab._integrity import (  # noqa: E402
    canonical_json_bytes,
    sha256_bytes,
    write_with_digest,
)
from ns_certificate_lab.leray_response_relay import leray_advection  # noqa: E402

from experiments.run_osgood_gate import (  # noqa: E402
    coherent_field,
    exact_anchor,
    measure,
)

V0 = (1, 2, 3)
V0_NORM_SQ = sum(c * c for c in V0)

EXACT_BANDS = (4, 6, 8)
FLOAT_ROWS = ((4, 48), (6, 48), (8, 80), (16, 144), (32, 144))
EXACT_WALL_CAP_SECONDS = 20 * 60
FLOAT_WALL_CAP_SECONDS = 10 * 60


# --------------------------------------------------------------------------
# exact lane
# --------------------------------------------------------------------------


def full_band(band: int) -> list[tuple[int, int, int]]:
    """Every lattice point with ``1 <= |k|^2 <= band^2`` (both signs)."""

    waves = []
    for kx in range(-band, band + 1):
        for ky in range(-band, band + 1):
            for kz in range(-band, band + 1):
                norm_sq = kx * kx + ky * ky + kz * kz
                if norm_sq == 0 or norm_sq > band * band:
                    continue
                waves.append((kx, ky, kz))
    return waves


def scaled_coefficients(
    waves: list[tuple[int, int, int]]
) -> tuple[int, dict[tuple[int, int, int], tuple[int, int, int]]]:
    """Return ``(M, A)`` with ``a_k = P_k v_0 / |k|^2 = A_k / M`` exactly.

    ``M = L^2`` where ``L = lcm{|k|^2}``, so ``M / |k|^4 = (L/|k|^2)^2`` is an
    integer and ``A_k = (v_0 |k|^2 - k (k.v_0)) (L/|k|^2)^2`` is an integer
    vector.  Pure integer arithmetic downstream; no rounding anywhere.
    """

    lcm = 1
    for wave in waves:
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        lcm = lcm * norm_sq // math.gcd(lcm, norm_sq)
    scale = lcm * lcm
    table: dict[tuple[int, int, int], tuple[int, int, int]] = {}
    for wave in waves:
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        dot = sum(c * v for c, v in zip(wave, V0))
        factor = (lcm // norm_sq) ** 2
        vector = tuple((V0[i] * norm_sq - wave[i] * dot) * factor for i in range(3))
        if sum(c * v for c, v in zip(wave, vector)) != 0:
            raise AssertionError("exact coefficient is not divergence free")
        table[wave] = vector  # type: ignore[assignment]
    return scale, table


def symmetry_certificate(band: int) -> dict[str, object]:
    """Exact check of ``H_0=(2/3)|v_0|^2 T_N`` and ``H_1=(2/3)|v_0|^2 S_N``."""

    waves = full_band(band)
    scale, table = scaled_coefficients(waves)
    h0 = Fraction(0)
    h1 = Fraction(0)
    s_n = Fraction(0)
    t_n = Fraction(0)
    scale_sq = Fraction(scale) ** 2
    for wave in waves:
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        vector = table[wave]
        modulus_sq = Fraction(sum(c * c for c in vector)) / scale_sq
        h0 += modulus_sq
        h1 += norm_sq * modulus_sq
        s_n += Fraction(1, norm_sq)
        t_n += Fraction(1, norm_sq * norm_sq)
    prefactor = Fraction(2, 3) * V0_NORM_SQ
    h0_rhs = prefactor * t_n
    h1_rhs = prefactor * s_n
    return {
        "band": band,
        "band_points": len(waves),
        "label": "exact",
        "S_N": str(s_n),
        "S_N_float": float(s_n),
        "T_N": str(t_n),
        "T_N_float": float(t_n),
        "H0": str(h0),
        "H0_float": float(h0),
        "H1": str(h1),
        "H1_float": float(h1),
        "H0_rhs_two_thirds_v0sq_TN": str(h0_rhs),
        "H1_rhs_two_thirds_v0sq_SN": str(h1_rhs),
        "H0_identity_holds_exactly": bool(h0 == h0_rhs),
        "H1_identity_holds_exactly": bool(h1 == h1_rhs),
        "H0_residual": str(h0 - h0_rhs),
        "H1_residual": str(h1 - h1_rhs),
        "N0_sq": str(h1 / h0),
        "N0_sq_float": float(h1 / h0),
    }


def exact_certificate(band: int) -> dict[str, object]:
    """Exact ``K``, ``Q(v_0)``, sweeping norm and pairing bound at one band."""

    started = time.time()
    waves = full_band(band)
    scale, table = scaled_coefficients(waves)
    items = list(table.items())

    # b_k = sum_{l+m=k} (a_l . m) a_m  carried as B_k = b_k * M^2 (integer).
    accumulator: dict[tuple[int, int, int], list[int]] = {}
    for left, left_vector in items:
        c0, c1, c2 = left_vector
        for right, right_vector in items:
            weight = c0 * right[0] + c1 * right[1] + c2 * right[2]
            if weight == 0:
                continue
            key = (
                left[0] + right[0],
                left[1] + right[1],
                left[2] + right[2],
            )
            current = accumulator.get(key)
            if current is None:
                accumulator[key] = [
                    weight * right_vector[0],
                    weight * right_vector[1],
                    weight * right_vector[2],
                ]
            else:
                current[0] += weight * right_vector[0]
                current[1] += weight * right_vector[1]
                current[2] += weight * right_vector[2]
    convolution_seconds = time.time() - started

    scale_pow4 = Fraction(scale) ** 4
    scale_pow2 = Fraction(scale) ** 2
    scale_pow3 = Fraction(scale) ** 3

    # ||P(u.grad u)||_2^2 = sum_k |P_k b_k|^2, and N_hat(k) = i P_k b_k.
    nonlinear_sq = Fraction(0)
    for wave, vector in accumulator.items():
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        if norm_sq == 0:
            if any(component != 0 for component in vector):
                raise AssertionError("mean mode of u.grad u is nonzero")
            continue
        dot = wave[0] * vector[0] + wave[1] * vector[1] + wave[2] * vector[2]
        modulus_sq = vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
        nonlinear_sq += Fraction(norm_sq * modulus_sq - dot * dot, norm_sq)
    nonlinear_sq /= scale_pow4

    # H_0, H_1 and Q(v_0) = sum_k (v_0.k) (b_k . a_k); sweeping norm
    # ||(v_0.grad)u||_2^2 = sum_k (v_0.k)^2 |a_k|^2.
    h0 = Fraction(0)
    h1 = Fraction(0)
    pairing = Fraction(0)
    sweep_sq = Fraction(0)
    for wave, vector in items:
        norm_sq = wave[0] ** 2 + wave[1] ** 2 + wave[2] ** 2
        modulus_sq = Fraction(sum(c * c for c in vector)) / scale_pow2
        h0 += modulus_sq
        h1 += norm_sq * modulus_sq
        seed_dot = sum(c * v for c, v in zip(wave, V0))
        sweep_sq += seed_dot * seed_dot * modulus_sq
        response = accumulator.get(wave)
        if response is None:
            continue
        pairing += Fraction(
            seed_dot * sum(p * q for p, q in zip(response, vector))
        ) / scale_pow3

    k_exact = nonlinear_sq / h1**2
    # Cauchy-Schwarz with P((v_0.grad)u) = (v_0.grad)u (see note):
    #   ||P(u.grad u)||_2 >= |Q(v_0)| / ||(v_0.grad)u||_2.
    bound_sq = pairing * pairing / sweep_sq
    captured_sq = bound_sq / nonlinear_sq
    return {
        "band": band,
        "band_points": len(waves),
        "ordered_pairs_enumerated": len(waves) ** 2,
        "sumset_modes": len(accumulator),
        "label": "exact",
        "common_denominator_bits": int(scale.bit_length()),
        "convolution_seconds": convolution_seconds,
        "total_seconds": time.time() - started,
        "H0": str(h0),
        "H1": str(h1),
        "H1_float": float(h1),
        "nonlinear_sq": str(nonlinear_sq),
        "nonlinear_sq_float": float(nonlinear_sq),
        "nonlinear_norm_float": math.sqrt(float(nonlinear_sq)),
        "K": str(k_exact),
        "K_float": float(k_exact),
        "N0_sq": str(h1 / h0),
        "N0_sq_float": float(h1 / h0),
        "K_over_N0_sq_float": float(k_exact / (h1 / h0)),
        "Q_v0": str(pairing),
        "Q_v0_float": float(pairing),
        "Q_v0_over_band_sq_float": float(pairing) / float(band * band),
        "sweep_norm_sq": str(sweep_sq),
        "sweep_norm_sq_float": float(sweep_sq),
        "pairing_bound_sq": str(bound_sq),
        "pairing_bound_float": math.sqrt(float(bound_sq)),
        "pairing_captured_fraction_sq": str(captured_sq),
        "pairing_captured_fraction_float": math.sqrt(float(captured_sq)),
    }


# --------------------------------------------------------------------------
# float lane (reuses run_osgood_gate)
# --------------------------------------------------------------------------


def _integer_frequency_mesh(grid: int):
    frequency = np.rint(np.fft.fftfreq(grid, d=1.0 / grid)).astype(np.float64)
    return np.meshgrid(frequency, frequency, frequency, indexing="ij")


def float_row(band: int, grid: int) -> dict[str, object]:
    field = coherent_field(grid, band)
    metrics = measure(field)
    kx, ky, kz = _integer_frequency_mesh(grid)
    seed_dot = V0[0] * kx + V0[1] * ky + V0[2] * kz
    sweeping = 1.0j * seed_dot[None, ...] * field
    nonlinear = leray_advection(field, field)
    pairing = float(np.vdot(nonlinear, sweeping).real)
    sweep_sq = float(np.vdot(sweeping, sweeping).real)
    nonlinear_norm = math.sqrt(metrics["nonlinear_sq"])
    bound = abs(pairing) / math.sqrt(sweep_sq)
    return {
        "band": band,
        "grid": grid,
        "label": "float",
        # measure() stores the *half* energies; the pair-energy (L2 mean)
        # convention used by exact_anchor and by this script is twice those.
        "H0": 2.0 * metrics["H0"],
        "H1": metrics["grad_sq"],
        "measure_half_H0": metrics["H0"],
        "measure_half_H1": metrics["H1"],
        "K": metrics["K"],
        "N0_sq": metrics["N0_sq"],
        "K_over_N0_sq": metrics["K"] / metrics["N0_sq"],
        "nonlinear_sq": metrics["nonlinear_sq"],
        "nonlinear_norm": nonlinear_norm,
        "Q_v0": pairing,
        "Q_v0_over_band_sq": pairing / float(band * band),
        "sweep_norm_sq": sweep_sq,
        "pairing_bound": bound,
        "pairing_captured_fraction": bound / nonlinear_norm,
    }


# --------------------------------------------------------------------------


def _without_timings(value):
    """Drop wall-clock fields so the mathematical content hashes stably."""

    if isinstance(value, dict):
        return {
            key: _without_timings(item)
            for key, item in value.items()
            if not key.endswith("_seconds") and not key.endswith("cap_bound")
        }
    if isinstance(value, list):
        return [_without_timings(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=_REPO_ROOT
        / "outputs"
        / "verification_sprint_v1"
        / "osgood_gate"
        / "exact_family_certificates.json",
    )
    args = parser.parse_args()

    preregistration = {
        "exact_bands": list(EXACT_BANDS),
        "exact_enumeration": (
            "every lattice point with 1<=|k|^2<=N^2; all |band|^2 ordered "
            "pairs in the convolution; no truncation, no sampling"
        ),
        "float_rows": [list(row) for row in FLOAT_ROWS],
        "exact_wall_cap_seconds": EXACT_WALL_CAP_SECONDS,
        "float_wall_cap_seconds": FLOAT_WALL_CAP_SECONDS,
        "total_wall_cap_seconds": 30 * 60,
        "fallback_rule": (
            "if the N=6 exact convolution exceeds 20 minutes, certify N=5 "
            "instead and say so"
        ),
    }

    exact_start = time.time()
    symmetry = [symmetry_certificate(band) for band in EXACT_BANDS]
    for row in symmetry:
        print(
            f"symmetry N={row['band']}: H0={row['H0']} "
            f"identity={row['H0_identity_holds_exactly']} "
            f"H1={row['H1']} identity={row['H1_identity_holds_exactly']}",
            flush=True,
        )
    certificates = []
    for band in EXACT_BANDS:
        row = exact_certificate(band)
        certificates.append(row)
        print(
            f"exact N={band}: K={row['K_float']:.10f} "
            f"Q(v0)={row['Q_v0_float']:.6f} "
            f"captured={row['pairing_captured_fraction_float']:.6f} "
            f"({row['total_seconds']:.2f}s)",
            flush=True,
        )
    exact_seconds = time.time() - exact_start

    trigvector_anchor = exact_anchor(4)

    float_start = time.time()
    float_rows = []
    for band, grid in FLOAT_ROWS:
        row = float_row(band, grid)
        float_rows.append(row)
        print(
            f"float N={band} grid={grid}: K={row['K']:.10f} "
            f"K/N0^2={row['K_over_N0_sq']:.6f} "
            f"Q/N^2={row['Q_v0_over_band_sq']:.6f} "
            f"captured={row['pairing_captured_fraction']:.6f}",
            flush=True,
        )
    float_seconds = time.time() - float_start

    comparisons = []
    for row, sym in zip(certificates, symmetry):
        match = next(r for r in float_rows if r["band"] == row["band"])
        comparisons.append(
            {
                "band": row["band"],
                "grid": match["grid"],
                "K_exact_float": row["K_float"],
                "K_float_pipeline": match["K"],
                "K_relative_difference": abs(row["K_float"] - match["K"])
                / row["K_float"],
                "Q_exact_float": row["Q_v0_float"],
                "Q_float_pipeline": match["Q_v0"],
                "Q_relative_difference": abs(row["Q_v0_float"] - match["Q_v0"])
                / abs(row["Q_v0_float"]),
                "H0_relative_difference": abs(
                    sym["H0_float"] - match["H0"]
                )
                / sym["H0_float"],
                "H1_relative_difference": abs(
                    sym["H1_float"] - match["H1"]
                )
                / sym["H1_float"],
                "captured_exact": row["pairing_captured_fraction_float"],
                "captured_float": match["pairing_captured_fraction"],
            }
        )
        print(
            f"cross-check N={row['band']}: K rel diff "
            f"{comparisons[-1]['K_relative_difference']:.3e}, Q rel diff "
            f"{comparisons[-1]['Q_relative_difference']:.3e}",
            flush=True,
        )

    summary = {
        "schema": "ns-certificate-lab/coherent-family-certificates/v1",
        "status": "EXACT RATIONAL CERTIFICATES + FLOAT CONTINUATION",
        "family": (
            "u_hat(k) = P_k v0/|k|^2 on 1<=|k|<=N, v0=(1,2,3), all phases 1; "
            "identical to experiments/run_osgood_gate.py coherent_field"
        ),
        "convention": (
            "pair energy = L2 mean; H0=||u||_2^2, H1=||grad u||_2^2=D, "
            "K=||P(u.grad u)||_2^2/H1^2 (same as exact_anchor)"
        ),
        "preregistration": preregistration,
        "caps_binding": {
            "exact_seconds": exact_seconds,
            "float_seconds": float_seconds,
            "exact_cap_bound": exact_seconds > EXACT_WALL_CAP_SECONDS,
            "float_cap_bound": float_seconds > FLOAT_WALL_CAP_SECONDS,
            "fallback_to_N5_used": False,
        },
        "symmetry_identities": symmetry,
        "exact_certificates": certificates,
        "trigvector_anchor_band4": trigvector_anchor,
        "float_continuation": float_rows,
        "exact_vs_float": comparisons,
        "sweeping_projection_note": (
            "For constant c and divergence-free u the k-th Fourier "
            "coefficient of (c.grad)u is i(c.k) u_hat(k), which is orthogonal "
            "to k because k.u_hat(k)=0; hence P((c.grad)u)=(c.grad)u and "
            "Q(c)=<(u.grad)u,(c.grad)u>=<P(u.grad u),(c.grad)u>, so "
            "Cauchy-Schwarz gives ||P(u.grad u)||_2 >= |Q(c)|/||(c.grad)u||_2."
        ),
        "projection_identity_check": (
            "the exact lane evaluates Q(v0) WITHOUT the Leray projection "
            "(sum_k (v0.k)(b_k.a_k)) while the float lane evaluates it WITH "
            "the projection (vdot(leray_advection(u,u), i(v0.k)u_hat)); their "
            "agreement to ~1e-16 in exact_vs_float.Q_relative_difference is an "
            "independent numerical confirmation of P((v0.grad)u)=(v0.grad)u"
        ),
        "convention_note": (
            "run_osgood_gate.measure() returns half energies in its 'H0' and "
            "'H1' keys; the L2-mean pair-energy convention used here and by "
            "exact_anchor gives H0 = 2*measure['H0'] and H1 = "
            "measure['grad_sq'] = 2*measure['H1'] = D.  K is convention free "
            "because measure() already divides by grad_sq^2."
        ),
    }
    summary["content_digest_sha256"] = sha256_bytes(
        canonical_json_bytes(_without_timings(summary))
    )
    summary["content_digest_note"] = (
        "sha256 of this document with every *_seconds and *cap_bound field "
        "removed; the file-level .sha256 sidecar varies between runs only "
        "through those wall-clock fields"
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_with_digest(args.output, canonical_json_bytes(summary))
    print("content digest", summary["content_digest_sha256"])
    print("written", args.output)


if __name__ == "__main__":
    main()
