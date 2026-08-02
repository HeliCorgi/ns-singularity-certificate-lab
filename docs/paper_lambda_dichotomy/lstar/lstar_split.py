"""How badly does the sweeping split fail? Measure every triangle-type bound."""
import numpy as np
from lstar_direct import build, to_phys, advect, ip, V0, freqs


def report(N):
    n = 3 * N + 2
    uhat, Kg, K2, band = build(n, N)
    KX, KY, KZ = Kg
    dens = np.sum(np.abs(uhat) ** 2, axis=0)
    H1 = float(np.sum(K2 * dens))
    S_N = float(np.sum(np.where(band, 1.0 / np.where(K2 > 0, K2, 1.0), 0.0)))
    kdotv = KX * V0[0] + KY * V0[1] + KZ * V0[2]
    Du = 1j * kdotv[None, ...] * uhat
    Du2 = float(np.sum(kdotv ** 2 * dens))
    w = advect(uhat, Kg, n)
    Q = ip(w, Du)
    c = (2.0 / 3.0) * S_N * V0                     # u_N(0)
    first = (2.0 / 3.0) * S_N * Du2
    rem = Q - first
    up = to_phys(uhat, n).real                     # field is real
    ut = up - c[:, None, None, None]               # u - u(0)
    sup_u = float(np.max(np.linalg.norm(up, axis=0)))
    sup_ut = float(np.max(np.linalg.norm(ut, axis=0)))
    # best possible constant shift (Chebyshev centre, over the sampled grid)
    lo = up.reshape(3, -1).min(axis=1); hi = up.reshape(3, -1).max(axis=1)
    cheb = 0.5 * (lo + hi)
    sup_cheb = float(np.max(np.linalg.norm(up - cheb[:, None, None, None], axis=0)))
    # L^6 / L^3 route
    gradu = np.stack([to_phys(1j * Kg[d] * uhat[cc], n).real
                      for cc in range(3) for d in range(3)])
    g2 = np.sum(gradu ** 2, axis=0)
    L3grad = float(np.mean(g2 ** 1.5)) ** (1 / 3)
    ut6 = float(np.mean(np.sum(ut ** 2, axis=0) ** 3)) ** (1 / 6)
    v4 = float(V0 @ V0) ** 2
    print(f"N={N:3d} S={S_N:8.3f} | first={first/(v4*S_N**2):.5f}*v4*S^2 "
          f"Q={Q/(v4*S_N**2):+.5f} rem={rem/(v4*S_N**2):+.5f} | "
          f"sup|u|/(|v0|S)={sup_u/(np.linalg.norm(V0)*S_N):.4f} "
          f"sup|u-u(0)|/(|v0|S)={sup_ut/(np.linalg.norm(V0)*S_N):.4f} "
          f"cheb={sup_cheb/(np.linalg.norm(V0)*S_N):.4f}")
    b_inf = sup_ut * np.sqrt(H1) * np.sqrt(Du2)
    b_cheb_pref = sup_cheb * np.sqrt(H1) * np.sqrt(Du2)
    b_63 = ut6 * L3grad * np.sqrt(Du2)
    print(f"      naive |rem| bounds / first :  Linf={b_inf/first:7.3f}   "
          f"Cheb-best={b_cheb_pref/first:7.3f}   L6L3L2={b_63/first:7.3f}   "
          f"|  TRUE |rem|/first = {abs(rem)/first:.4f}")


for N in (4, 8, 12, 16, 24, 32):
    report(N)
