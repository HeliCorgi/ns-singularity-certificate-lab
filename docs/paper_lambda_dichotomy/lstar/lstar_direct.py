"""JOB A numerics: sweeping-split diagnostics for the sharply truncated family.

Convention (full lattice, matches complete_proof.md):
    u(x) = sum_{k in Z^3\{0}} uhat_k e^{i k.x},   uhat_k = P_k v0 / |k|^2, 1<=|k|<=N
    <f,g> = (2pi)^{-3} int f.g = sum_k conj(fhat_k).ghat_k
    H_r   = sum_k |k|^{2r} |uhat_k|^2
All float64 unless labelled exact.
"""
import sys
import numpy as np

V0 = np.array([1.0, 2.0, 3.0])


def freqs(n):
    return np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(np.int64)


def build(n, N):
    f = freqs(n).astype(np.float64)
    KX, KY, KZ = np.meshgrid(f, f, f, indexing="ij")
    K2 = KX * KX + KY * KY + KZ * KZ
    band = (K2 >= 1.0) & (K2 <= N * N + 1e-9)
    kdotv = KX * V0[0] + KY * V0[1] + KZ * V0[2]
    uhat = np.zeros((3, n, n, n), dtype=np.complex128)
    inv = np.zeros_like(K2)
    inv[band] = 1.0 / K2[band]
    for c, (KC, vc) in enumerate(zip((KX, KY, KZ), V0)):
        uhat[c] = np.where(band, (vc - KC * kdotv * inv) * inv, 0.0)
    return uhat, (KX, KY, KZ), K2, band


def to_phys(a, n):
    return np.fft.ifftn(a, axes=(-3, -2, -1)) * float(n ** 3)


def to_spec(a, n):
    return np.fft.fftn(a, axes=(-3, -2, -1)) / float(n ** 3)


def advect(uhat, Kg, n):
    """Fourier coefficients of (u.grad)u."""
    up = to_phys(uhat, n)
    out = np.zeros_like(uhat)
    for comp in range(3):
        acc = np.zeros((n, n, n), dtype=np.complex128)
        for d in range(3):
            acc += up[d] * to_phys(1j * Kg[d] * uhat[comp], n)
        out[comp] = to_spec(acc, n)
    return out


def leray(what, Kg, K2):
    kd = np.zeros_like(what[0])
    for d in range(3):
        kd += Kg[d] * what[d]
    safe = np.where(K2 > 0, K2, 1.0)
    out = np.empty_like(what)
    for d in range(3):
        out[d] = np.where(K2 > 0, what[d] - Kg[d] * kd / safe, 0.0)
    return out


def ip(a, b):
    return float(np.real(np.sum(np.conj(a) * b)))


def run(N, full=True):
    n = (4 * N + 2) if full else (3 * N + 2)
    uhat, Kg, K2, band = build(n, N)
    KX, KY, KZ = Kg
    dens = np.sum(np.abs(uhat) ** 2, axis=0)
    H0 = float(np.sum(dens))
    H1 = float(np.sum(K2 * dens))
    S_N = float(np.sum(np.where(band, 1.0 / np.where(K2 > 0, K2, 1.0), 0.0)))
    T_N = float(np.sum(np.where(band, 1.0 / np.where(K2 > 0, K2, 1.0) ** 2, 0.0)))
    kdotv = KX * V0[0] + KY * V0[1] + KZ * V0[2]
    Du = 1j * kdotv[None, ...] * uhat                      # (v0.grad)u
    Du2 = float(np.sum(kdotv ** 2 * dens))                 # ||(v0.grad)u||^2
    w = advect(uhat, Kg, n)                                # (u.grad)u
    Q = ip(w, Du)
    first = (2.0 / 3.0) * S_N * Du2                        # <(u(0).grad)u, Du>
    rem = Q - first
    out = dict(N=N, n=n, S_N=S_N, T_N=T_N, H0=H0, H1=H1, N0sq=H1 / H0,
               Du2=Du2, Q=Q, first=first, rem=rem)
    if full:
        Pw = leray(w, Kg, K2)
        nsq = ip(Pw, Pw)
        out["PN2"] = nsq
        out["K"] = nsq / H1 ** 2
        out["pair_lb"] = Q * Q / Du2
        out["pair_frac"] = np.sqrt(Q * Q / Du2 / nsq)
        out["raw2"] = ip(w, w)
    # sup norm and point value
    up = to_phys(uhat, n)
    out["sup_u"] = float(np.max(np.sqrt(np.sum(np.abs(up) ** 2, axis=0))))
    out["u0_norm"] = (2.0 / 3.0) * S_N * float(np.linalg.norm(V0))
    return out


if __name__ == "__main__":
    Ns = [int(x) for x in sys.argv[1:]] or [4, 6, 8, 12, 16, 20, 24]
    print(f"{'N':>4} {'S_N':>9} {'Du2/S_N':>10} {'Q/N^2':>10} {'first/N^2':>10} "
          f"{'rem/N^2':>11} {'|rem|/first':>11} {'PN2/N^3':>9} {'K':>8} "
          f"{'pairfrac':>8}")
    for N in Ns:
        r = run(N)
        v4 = float(V0 @ V0) ** 2
        print(f"{r['N']:>4} {r['S_N']:>9.3f} {r['Du2']/r['S_N']/v4:>10.6f} "
              f"{r['Q']/N**2:>10.2f} {r['first']/N**2:>10.2f} "
              f"{r['rem']/N**2:>11.2f} {abs(r['rem'])/r['first']:>11.4f} "
              f"{r['PN2']/N**3:>9.3f} {r['K']:>8.4f} {r['pair_frac']:>8.4f}")
