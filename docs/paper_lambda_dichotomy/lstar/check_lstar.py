"""Numerical corroboration for the L* smooth-family proof (float pipeline).

Checks:
 (A) physical-space:  curl(W x omega) == 12 y_3 (e_3 x y)/r^6   for W = e3/r + y3 y/r^3
 (B) exact moment laws H_0, H_1, u_N(0) for the smoothly truncated family
 (C) trend of ||P(u_N . grad u_N)||^2 / N^3   (pseudo-spectral, dealiased exactly)
"""
import numpy as np

# ---------------------------------------------------------------- (A)
def W(y):
    r = np.linalg.norm(y)
    e3 = np.array([0.0, 0.0, 1.0])
    return e3 / r + y[2] * y / r**3

def curl(f, y, h=1e-5):
    J = np.zeros((3, 3))
    for i in range(3):
        dy = np.zeros(3); dy[i] = h
        J[:, i] = (f(y + dy) - f(y - dy)) / (2 * h)      # J[j,i] = d_i f_j
    return np.array([J[2,1]-J[1,2], J[0,2]-J[2,0], J[1,0]-J[0,1]])

def omega(y):
    return curl(W, y)

def Wxomega(y):
    return np.cross(W(y), omega(y))

print("=== (A) curl(W x omega) vs 12 y3 (e3 x y)/r^6 ===")
rng = np.random.default_rng(0)
for _ in range(5):
    y = rng.normal(size=3) * 1.3
    r = np.linalg.norm(y)
    lhs = curl(Wxomega, y, h=1e-4)
    rhs = 12 * y[2] * np.cross(np.array([0.,0.,1.]), y) / r**6
    print(f"  y={np.round(y,3)}  lhs={np.round(lhs,6)}  rhs={np.round(rhs,6)}  "
          f"relerr={np.linalg.norm(lhs-rhs)/max(np.linalg.norm(rhs),1e-300):.2e}")

# also: divergence-free and omega formula
def div(f, y, h=1e-5):
    s = 0.0
    for i in range(3):
        dy = np.zeros(3); dy[i] = h
        s += (f(y+dy)-f(y-dy))[i] / (2*h)
    return s
y = np.array([0.7, -1.1, 0.4])
print("  div W =", div(W, y),
      " omega err =", np.linalg.norm(omega(y) - 2*np.cross([0,0,1], y)/np.linalg.norm(y)**3))

# ---------------------------------------------------------------- (B)/(C)
def chi(t):
    """smooth, ==1 on [0,1/2], supp in [0,1], 0<=chi<=1, nonincreasing."""
    t = np.asarray(t, dtype=float)
    s = np.clip((t - 0.5) / 0.5, 0.0, 1.0)
    # C^inf transition
    out = np.ones_like(s)
    m = (s > 0) & (s < 1)
    a = np.exp(-1.0 / s[m]); b = np.exp(-1.0 / (1.0 - s[m]))
    out[m] = b / (a + b)
    out[s >= 1] = 0.0
    return out

def family(N, v0=np.array([1.0, 0.0, 0.0]), M=None):
    """Return grid-space u_N and its spectral data on an M^3 grid."""
    if M is None:
        M = 4 * N + 2
    k1 = np.fft.fftfreq(M, d=1.0 / M).astype(int)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    K2 = (KX**2 + KY**2 + KZ**2).astype(float)
    Kn = np.sqrt(K2)
    w = np.zeros_like(K2)
    nz = K2 > 0
    w[nz] = chi(Kn[nz] / N) / K2[nz]
    Kv = KX * v0[0] + KY * v0[1] + KZ * v0[2]
    uh = np.empty((3,) + K2.shape)
    Ks = [KX, KY, KZ]
    for i in range(3):
        comp = np.zeros_like(K2)
        comp[nz] = w[nz] * (v0[i] - Ks[i][nz] * Kv[nz] / K2[nz])
        uh[i] = comp
    return uh, (KX, KY, KZ), K2, M

print("\n=== (B) exact moment laws (float) ===")
v0 = np.array([1.0, 2.0, -1.0])
for N in (4, 8, 16):
    uh, (KX, KY, KZ), K2, M = family(N, v0)
    nz = K2 > 0
    Kn = np.sqrt(K2)
    w2 = np.zeros_like(K2); w2[nz] = chi(Kn[nz]/N)**2
    S = np.sum(w2[nz] / K2[nz])            # sum chi^2/|k|^2
    T = np.sum(w2[nz] / K2[nz]**2)         # sum chi^2/|k|^4
    Sc = np.sum(np.where(nz, np.nan_to_num(chi(np.sqrt(np.where(nz,K2,1)))/1), 0))  # unused
    chi1 = np.zeros_like(K2); chi1[nz] = chi(Kn[nz]/N)
    S1 = np.sum(chi1[nz] / K2[nz])         # sum chi/|k|^2
    H0 = np.sum(uh**2)
    H1 = np.sum(K2 * uh**2)
    u0 = uh.sum(axis=(1, 2, 3))
    nv = v0 @ v0
    print(f"  N={N:3d}  H0 vs (2/3)|v0|^2 T : {H0:.10f} {2/3*nv*T:.10f}"
          f"   H1 vs (2/3)|v0|^2 S : {H1:.10f} {2/3*nv*S:.10f}")
    print(f"          u_N(0) vs (2/3) S1 v0 : {np.round(u0,8)}  {np.round(2/3*S1*v0,8)}")
    print(f"          N_0^2 = S/T = {S/T:.4f}   ratio to N = {S/T/N:.4f}")

print("\n=== (C) ||P(u.grad u)||_2^2 / N^3  (normalised measure) ===")
prev = None
for N in (4, 8, 12, 16, 24, 32, 40):
    uh, (KX, KY, KZ), K2, M = family(N, np.array([1.0, 0.0, 0.0]))
    u = np.array([np.fft.ifftn(uh[i]).real * M**3 for i in range(3)])
    # grad u in physical space
    adv = np.zeros_like(u)
    Ks = [KX, KY, KZ]
    for j in range(3):
        for i in range(3):
            duj_di = np.fft.ifftn(1j * Ks[i] * uh[j]).real * M**3
            adv[j] += u[i] * duj_di
    advh = np.array([np.fft.fftn(adv[j]) / M**3 for j in range(3)])
    # Leray project
    nz = K2 > 0
    dot = np.zeros_like(K2, dtype=complex)
    for i in range(3):
        dot += Ks[i] * advh[i]
    proj = advh.copy()
    for i in range(3):
        proj[i][nz] -= Ks[i][nz] * dot[nz] / K2[nz]
    proj[:, ~nz] = 0.0
    nrm2 = float(np.sum(np.abs(proj) ** 2))
    H1 = float(np.sum(K2 * uh**2))
    Kfun = nrm2 / H1**2
    print(f"  N={N:3d}  ||P(u.grad u)||^2={nrm2:12.5f}   /N^3 = {nrm2/N**3:.6f}"
          f"   K={Kfun:.5f}   K/N={Kfun/N:.5f}")
