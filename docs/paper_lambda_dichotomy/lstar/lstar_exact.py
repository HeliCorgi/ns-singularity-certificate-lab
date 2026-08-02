"""Exact (Fraction) verification of the lattice identities used in JOB A."""
from fractions import Fraction as F
import itertools, sys

V0 = (1, 2, 3)


def lattice(N):
    for k in itertools.product(range(-N, N + 1), repeat=3):
        n2 = k[0] ** 2 + k[1] ** 2 + k[2] ** 2
        if 1 <= n2 <= N * N:
            yield k, n2


def sums(N):
    S = T = F(0); B = A = F(0); Du2 = F(0); Sig4 = F(0)
    v2 = sum(v * v for v in V0)
    for k, n2 in lattice(N):
        S += F(1, n2); T += F(1, n2 * n2)
        A += F(k[0] ** 4, n2 ** 3)
        B += F(k[0] ** 2 * k[1] ** 2, n2 ** 3)
        kv = sum(a * b for a, b in zip(k, V0))
        # |P_k v0|^2 = v2 - kv^2/n2 ; |uhat|^2 = |P_k v0|^2/n2^2
        pk2 = F(v2) - F(kv * kv, n2)
        Du2 += F(kv * kv) * pk2 / F(n2 * n2)
        Sig4 += F(kv ** 4, n2 ** 3)
    return S, T, A, B, Du2, Sig4


for N in [int(x) for x in (sys.argv[1:] or ["4", "6", "8", "10", "12"])]:
    S, T, A, B, Du2, Sig4 = sums(N)
    v2 = F(sum(v * v for v in V0)); v4 = v2 * v2
    rho = F(sum(v ** 4 for v in V0), int(v2) ** 2)
    # claimed: Sig4_normalised = sum cos^4/|k|^2 = 3B + (S/3 - 5B) rho
    claim_cos4 = 3 * B + (S / 3 - 5 * B) * rho
    assert Sig4 / v4 == claim_cos4, (N, Sig4 / v4, claim_cos4)
    # claimed: ||Du||^2 = v4 * (S/3 - sum cos^4/|k|^2)
    claim_Du2 = v4 * (S / 3 - claim_cos4)
    assert Du2 == claim_Du2, (N, Du2, claim_Du2)
    # contraction identity 3A + 6B = S
    assert 3 * A + 6 * B == S
    aniso = S / 3 - 5 * B                       # 0 iff isotropic 4th moment
    print(f"N={N:3d}  S_N={float(S):9.4f}  B/S_N={float(B/S):.8f} (iso 1/15="
          f"{1/15:.8f})  (S/3-5B)={float(aniso):+9.5f}  "
          f"(S/3-5B)/S={float(aniso/S):+.6f}  ||Du||^2/(v4 S)={float(Du2/(v4*S)):.8f}"
          f"  [2/15={2/15:.8f}]")
