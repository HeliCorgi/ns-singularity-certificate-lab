"""Where does Q live, and is c=v0 the optimal sweeping direction?"""
import numpy as np
from lstar_direct import build, to_phys, advect, ip, V0


def conc(N):
    n = 3 * N + 2
    uhat, Kg, K2, band = build(n, N)
    KX, KY, KZ = Kg
    S_N = float(np.sum(np.where(band, 1.0 / np.where(K2 > 0, K2, 1.0), 0.0)))
    kdotv = KX * V0[0] + KY * V0[1] + KZ * V0[2]
    w = advect(uhat, Kg, n)
    Du = 1j * kdotv[None, ...] * uhat
    Q = ip(w, Du)
    wp = to_phys(w, n).real
    Dup = to_phys(Du, n).real
    integ = np.sum(wp * Dup, axis=0)                # pointwise Q-density
    # M_ij = sum k_i k_j |uhat_k|^2  (checks c=v0 optimality claim)
    dens = np.sum(np.abs(uhat) ** 2, axis=0)
    M = np.array([[float(np.sum(a * b * dens)) for b in (KX, KY, KZ)]
                  for a in (KX, KY, KZ)])
    # sup |u| vs exact point value (2/3) S_N |v0|
    up = to_phys(uhat, n).real
    sup = float(np.max(np.linalg.norm(up, axis=0)))
    exact_pt = (2.0 / 3.0) * S_N * float(np.linalg.norm(V0))
    # radial profile of the Q density about x=0
    g = np.arange(n)
    x = np.where(g <= n // 2, g, g - n) * (2 * np.pi / n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X * X + Y * Y + Z * Z)
    tot = float(np.mean(integ))
    print(f"\nN={N}: Q={Q:.6e}  mean(integrand)={tot:.6e} (agree "
          f"{abs(tot-Q)/abs(Q):.2e})   sup|u|/[(2/3)S|v0|]={sup/exact_pt:.12f}")
    print(f"   M/H1 eigenstructure: M v0 || v0 ? cos="
          f"{float(M@V0 @ V0)/np.linalg.norm(M@V0)/np.linalg.norm(V0):.6f}")
    for mult in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        r0 = mult / N
        frac = float(np.mean(np.where(R <= r0, integ, 0.0))) / tot
        print(f"   frac of Q from |x| <= {mult:4.1f}/N : {frac:8.4f}")


for N in (8, 16, 24):
    conc(N)
