import numpy as np

def hvec(eta):
    # eta: (...,3); h = (e3 - eta*eta3/|eta|^2)/|eta|^2
    r2 = np.sum(eta*eta,axis=-1)
    out = np.zeros_like(eta)
    out[...,2] = 1.0
    out = out - eta*(eta[...,2]/r2)[...,None]
    return out/r2[...,None]

def gl(n):
    x,w = np.polynomial.legendre.leggauss(n)
    return x,w

def tau(zeta, nr=180, nt=90, np_=90):
    zeta=np.asarray(zeta,float)
    # nodes
    xr,wr = gl(nr); xt,wt = gl(nt); xp,wp = gl(np_)
    t = 0.5*(xr+1); wt_r = 0.5*wr           # t in (0,1)
    mu = xt; wmu = wt                        # cos polar in (-1,1)
    ph = np.pi*(xp+1); wph = np.pi*wp        # azimuth (0,2pi)
    T = np.zeros((3,3))
    # ---- Part A: spherical about 0, radius r = t/(1-t)
    r = t/(1-t); jr = 1.0/(1-t)**2
    R,M,P = np.meshgrid(r,mu,ph,indexing='ij')
    WR,WM,WP = np.meshgrid(wt_r*jr,wmu,wph,indexing='ij')
    s = np.sqrt(1-M**2)
    eta = np.stack([R*s*np.cos(P), R*s*np.sin(P), R*M],axis=-1)
    d0 = np.linalg.norm(eta,axis=-1); dz = np.linalg.norm(eta-zeta,axis=-1)
    w0 = dz**4/(d0**4+dz**4)
    ha = hvec(eta); hb = hvec(zeta-eta)
    wgt = WR*WM*WP*(R**2)*w0
    T += np.einsum('ni,nj,n->ij', ha.reshape(-1,3), hb.reshape(-1,3), wgt.reshape(-1))
    # ---- Part B: spherical about zeta, radius rho = t/(1-t)
    eta = zeta + np.stack([R*s*np.cos(P), R*s*np.sin(P), R*M],axis=-1)
    d0 = np.linalg.norm(eta,axis=-1); dz = np.linalg.norm(eta-zeta,axis=-1)
    wz = d0**4/(d0**4+dz**4)
    ha = hvec(eta); hb = hvec(zeta-eta)
    wgt = WR*WM*WP*(R**2)*wz
    T += np.einsum('ni,nj,n->ij', ha.reshape(-1,3), hb.reshape(-1,3), wgt.reshape(-1))
    return T

for z in [np.array([0.6,0.0,0.8]), np.array([0.3,0.4,0.5]), np.array([1.0,1.0,2.0])]:
    for n in [(120,60,60),(200,100,100)]:
        T = tau(z,*n)
        cross = np.cross(z, T@z)
        pred = (3*np.pi**3/8)*(z[2]/np.linalg.norm(z))*np.cross(z,[0,0,1])
        print(z, n, "num",np.round(cross,5), "pred",np.round(pred,5),
              "relerr", np.linalg.norm(cross-pred)/np.linalg.norm(pred))
    print("  sym resid", np.abs(T-T.T).max())
