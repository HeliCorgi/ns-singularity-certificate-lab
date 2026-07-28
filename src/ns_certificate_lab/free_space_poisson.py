"""Free-radial Green solver for the axisymmetric five-dimensional Poisson operator.

The regularized axisymmetric stream-function equation is

    -L5 psi = source,
    L5 = d_rr + (3/r) d_r + d_zz.

For a Fourier mode exp(i k z), the radial operator is the radial Laplacian in
four dimensions.  Its whole-space Green kernel is explicit:

    G_0(r, rho) = 1 / (2 max(r, rho)^2),

and, for k > 0,

    G_k(r, rho) = I_1(k r_<) K_1(k r_>) / (r_< r_>).

The apparent r_< = 0 singularity is removable because I_1(x) / x -> 1/2.
The source is integrated with the four-dimensional radial measure rho^3 d rho.

The z direction remains periodic after zero padding.  Increasing ``pad_factor``
therefore provides a controlled period-image sensitivity study; it is not, by
itself, a rigorous free-space-z certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class FreeSpacePoissonDiagnostics:
    """Metadata needed to audit a free-radial/periodized-z solve."""

    radial_points: int
    original_z_points: int
    padded_z_points: int
    radial_max: float
    dz: float
    padded_period: float
    pad_factor: int
    quadrature: str
    bessel_backend: str
    source_linf: float
    solution_linf: float


def _i1e_cephes(x: np.ndarray) -> np.ndarray:
    """Approximate exp(-x) I_1(x) for x >= 0 without overflow.

    Coefficients are the classical Cephes/Numerical-Recipes piecewise
    approximants.  The function is deterministic and has roughly single-
    precision relative accuracy, sufficient for the independent prototype
    gate.  A rigorous solver must eventually replace this approximation with
    interval-certified special functions.
    """

    x = np.asarray(x, dtype=float)
    if np.any(x < 0.0):
        raise ValueError("_i1e_cephes expects x >= 0")

    out = np.empty_like(x)
    small = x < 3.75
    if np.any(small):
        xs = x[small]
        y = (xs / 3.75) ** 2
        poly = 0.5 + y * (
            0.87890594
            + y
            * (
                0.51498869
                + y
                * (
                    0.15084934
                    + y * (0.02658733 + y * (0.00301532 + y * 0.00032411))
                )
            )
        )
        out[small] = xs * poly * np.exp(-xs)

    large = ~small
    if np.any(large):
        xl = x[large]
        y = 3.75 / xl
        poly = 0.39894228 + y * (
            -0.03988024
            + y
            * (
                -0.00362018
                + y
                * (
                    0.00163801
                    + y
                    * (
                        -0.01031555
                        + y
                        * (
                            0.02282967
                            + y * (-0.02895312 + y * (0.01787654 - y * 0.00420059))
                        )
                    )
                )
            )
        )
        out[large] = poly / np.sqrt(xl)

    return out


def _k1e_cephes(x: np.ndarray) -> np.ndarray:
    """Approximate exp(x) K_1(x) for x > 0 without underflow."""

    x = np.asarray(x, dtype=float)
    if np.any(x <= 0.0):
        raise ValueError("_k1e_cephes expects x > 0")

    out = np.empty_like(x)
    small = x <= 2.0
    if np.any(small):
        xs = x[small]
        y = 0.25 * xs * xs
        i1 = _i1e_cephes(xs) * np.exp(xs)
        poly = 1.0 + y * (
            0.15443144
            + y
            * (
                -0.67278579
                + y
                * (
                    -0.18156897
                    + y
                    * (
                        -0.01919402
                        + y * (-0.00110404 + y * (-0.00004686))
                    )
                )
            )
        )
        k1 = np.log(0.5 * xs) * i1 + poly / xs
        out[small] = np.exp(xs) * k1

    large = ~small
    if np.any(large):
        xl = x[large]
        y = 2.0 / xl
        poly = 1.25331414 + y * (
            0.23498619
            + y
            * (
                -0.03655620
                + y
                * (
                    0.01504268
                    + y * (-0.00780353 + y * (0.00325614 + y * (-0.00068245)))
                )
            )
        )
        out[large] = poly / np.sqrt(xl)

    return out


def _bessel_product_i1k1(
    k: float,
    r_less: np.ndarray,
    r_greater: np.ndarray,
    *,
    backend: str,
) -> np.ndarray:
    """Return I1(k*r_less) K1(k*r_greater)/(r_less*r_greater)."""

    if not k > 0.0:
        raise ValueError("k must be positive")
    a = k * np.asarray(r_less, dtype=float)
    b = k * np.asarray(r_greater, dtype=float)
    if np.any(a > b + 1e-14):
        raise ValueError("r_less must not exceed r_greater")

    positive_b = b > 0.0
    scaled_i = np.zeros_like(a, dtype=float)
    scaled_k = np.zeros_like(b, dtype=float)
    if backend == "scipy":
        try:
            from scipy.special import ive, kve  # type: ignore
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise RuntimeError("SciPy backend requested but scipy is unavailable") from exc
        scaled_i[...] = ive(1, a)
        scaled_k[positive_b] = kve(1, b[positive_b])
    elif backend == "cephes":
        scaled_i[...] = _i1e_cephes(a)
        scaled_k[positive_b] = _k1e_cephes(b[positive_b])
    else:
        raise ValueError("backend must be 'cephes' or 'scipy'")

    exponential = np.exp(a - b)
    out = np.empty(np.broadcast(a, b).shape, dtype=float)
    rl = np.broadcast_to(np.asarray(r_less, dtype=float), out.shape)
    rg = np.broadcast_to(np.asarray(r_greater, dtype=float), out.shape)
    positive = rl > 0.0
    out[positive] = (
        scaled_i[positive]
        * scaled_k[positive]
        * exponential[positive]
        / (rl[positive] * rg[positive])
    )

    # I1(k r)/r -> k/2 as r -> 0.
    zero = ~positive
    if np.any(zero):
        bz = b[zero]
        if np.any(bz <= 0.0):
            # The only such entry is r = rho = 0.  Its quadrature weight rho^3
            # is exactly zero, so setting the kernel to zero is harmless and
            # avoids an undefined point value.
            valid = bz > 0.0
            values = np.zeros_like(bz)
            if np.any(valid):
                if backend == "scipy":
                    from scipy.special import kve  # type: ignore

                    kscaled = kve(1, bz[valid])
                else:
                    kscaled = _k1e_cephes(bz[valid])
                values[valid] = 0.5 * k * kscaled * np.exp(-bz[valid]) / rg[zero][valid]
            out[zero] = values
        else:
            out[zero] = 0.5 * k * scaled_k[zero] * np.exp(-b[zero]) / rg[zero]

    return out


def radial_green_matrix(
    r: np.ndarray,
    k: float,
    *,
    bessel_backend: str = "cephes",
) -> np.ndarray:
    """Build the whole-space radial Green matrix for one Fourier wavenumber."""

    r = np.asarray(r, dtype=float)
    if r.ndim != 1 or r.size < 2:
        raise ValueError("r must be a one-dimensional grid with at least two points")
    if r[0] != 0.0 or np.any(np.diff(r) <= 0.0):
        raise ValueError("r must start at zero and be strictly increasing")
    if k < 0.0 or not math.isfinite(k):
        raise ValueError("k must be finite and nonnegative")

    ri = r[:, None]
    rho = r[None, :]
    r_less = np.minimum(ri, rho)
    r_greater = np.maximum(ri, rho)

    if k == 0.0:
        out = np.zeros((r.size, r.size), dtype=float)
        positive = r_greater > 0.0
        out[positive] = 0.5 / (r_greater[positive] ** 2)
        return out

    return _bessel_product_i1k1(
        k,
        r_less,
        r_greater,
        backend=bessel_backend,
    )


def _radial_quadrature_weights(r: np.ndarray) -> tuple[np.ndarray, str]:
    """Return composite-trapezoid weights on a uniform radial grid.

    The modal Green kernel has a derivative jump at ``rho = r``.  A global
    Simpson rule straddles that kink and is not uniformly higher order; the
    trapezoid rule is the conservative independent choice for this prototype.
    """

    if r.ndim != 1 or r.size < 2:
        raise ValueError("r must contain at least two points")
    dr = np.diff(r)
    if not np.allclose(dr, dr[0], rtol=1e-12, atol=1e-14):
        raise ValueError("the prototype currently requires a uniform radial grid")
    weights = np.full(r.size, float(dr[0]), dtype=float)
    weights[[0, -1]] *= 0.5
    return weights, "uniform composite trapezoid in rho with rho^3 measure"


def solve_l5_free_radial_periodized_z(
    source: np.ndarray,
    r: np.ndarray,
    dz: float,
    *,
    pad_factor: int = 1,
    bessel_backend: str = "cephes",
) -> tuple[np.ndarray, FreeSpacePoissonDiagnostics]:
    """Solve ``-L5 psi = source`` with free radial decay and padded-periodic z.

    Parameters
    ----------
    source:
        Real array with shape ``(Nr, Nz)``.  The z grid is uniform and excludes
        the repeated periodic endpoint.
    r:
        Uniform radial node grid starting at zero.
    dz:
        z grid spacing.
    pad_factor:
        Integer zero-padding factor in z.  ``1`` solves the periodic-z problem.
        Larger factors separate periodic images while retaining the same local
        spacing.
    bessel_backend:
        ``"cephes"`` is dependency-free.  ``"scipy"`` uses scaled special
        functions when SciPy is available.
    """

    source = np.asarray(source, dtype=float)
    r = np.asarray(r, dtype=float)
    if source.ndim != 2 or source.shape[0] != r.size:
        raise ValueError("source must have shape (len(r), Nz)")
    if not np.all(np.isfinite(source)):
        raise ValueError("source must contain only finite values")
    if not (math.isfinite(dz) and dz > 0.0):
        raise ValueError("dz must be finite and positive")
    if not isinstance(pad_factor, int) or pad_factor < 1:
        raise ValueError("pad_factor must be a positive integer")
    if source.shape[1] < 4:
        raise ValueError("at least four z points are required")

    base_weights, quadrature_name = _radial_quadrature_weights(r)
    radial_weights = base_weights * (r**3)
    nr, nz = source.shape
    nz_padded = pad_factor * nz
    start = (nz_padded - nz) // 2
    padded = np.zeros((nr, nz_padded), dtype=float)
    padded[:, start : start + nz] = source

    source_hat = np.fft.rfft(padded, axis=1)
    solution_hat = np.empty_like(source_hat, dtype=complex)
    wavenumbers = 2.0 * math.pi * np.fft.rfftfreq(nz_padded, d=dz)

    weighted_source = source_hat * radial_weights[:, None]
    for mode, k in enumerate(wavenumbers):
        green = radial_green_matrix(r, float(k), bessel_backend=bessel_backend)
        solution_hat[:, mode] = green @ weighted_source[:, mode]

    padded_solution = np.fft.irfft(solution_hat, n=nz_padded, axis=1)
    solution = padded_solution[:, start : start + nz].copy()

    diagnostics = FreeSpacePoissonDiagnostics(
        radial_points=nr,
        original_z_points=nz,
        padded_z_points=nz_padded,
        radial_max=float(r[-1]),
        dz=float(dz),
        padded_period=float(nz_padded * dz),
        pad_factor=pad_factor,
        quadrature=quadrature_name,
        bessel_backend=bessel_backend,
        source_linf=float(np.max(np.abs(source))),
        solution_linf=float(np.max(np.abs(solution))),
    )
    return solution, diagnostics
