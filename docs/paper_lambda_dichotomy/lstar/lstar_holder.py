"""Scan the whole Holder family for the remainder bound; also M's spectrum."""
import numpy as np
from lstar_direct import build, to_phys, advect, ip, V0


def scan(N):
    n = 3 * N + 2
    uhat, Kg, K2, band = build(n, N)
    KX, KY, KZ = Kg
    dens = np.sum(np.abs(uhat) ** 2, axis=0)
    H1 = float(np.sum(K2 * dens))
    S_N = float(np.sum(np.where(band, 1.0 / np.where(K2 > 0, K2, 1.0), 0.0)))
    kdotv = KX * V0[0] + KY * V0[1] + KZ * V0[2]
    Du2 = float(np.sum(kdotv ** 2 * dens))
    M = np.array([[float(np.sum(a * b * dens)) for b in (KX, KY, KZ)]
                  for a in (KX, KY, KZ)])
    ev = np.linalg.eigvalsh(M)
    up = to_phys(uhat, n).real
    b = (2.0 / 3.0) * S_N * V0
    ut = np.linalg.norm(up - b[:, None, None, None], axis=0)
    gr = np.stack([to_phys(1j * Kg[d] * uhat[c], n).real
                   for c in range(3) for d in range(3)])
    gm = np.sqrt(np.sum(gr ** 2, axis=0))
    first = (2.0 / 3.0) * S_N * Du2
    out = []
    for p in (2.5, 3.0, 4.0, 6.0, 10.0, 20.0, np.inf):
        if np.isinf(p):
            np_, nq = float(ut.max()), np.sqrt(H1)
        else:
            q = 1.0 / (0.5 - 1.0 / p)
            np_ = float(np.mean(ut ** p)) ** (1 / p)
            nq = float(np.mean(gm ** q)) ** (1 / q)
        out.append(np_ * nq * np.sqrt(Du2) / first)
    print(f"N={N:3d} lam(M)/H1={ev/H1}  |  bound/first by p="
          + "  ".join(f"{p if not np.isinf(p) else 'inf'}:{v:.3f}"
                      for p, v in zip((2.5, 3, 4, 6, 10, 20, np.inf), out)))


for N in (8, 16, 24, 32):
    scan(N)
print("\nsqrt(5) =", np.sqrt(5), "  sqrt(8/3) =", np.sqrt(8 / 3))
