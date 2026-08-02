import numpy as np
def chi(r):
    s=np.clip((r*r-0.25)/0.75,0,1)
    S=126*s**5-420*s**6+540*s**7-315*s**8+70*s**9
    return 1-S
v0=np.array([0.,0.,1.])
def uN_over_N(N,y):
    a=np.arange(-N,N+1)
    K=np.stack(np.meshgrid(a,a,a,indexing='ij'),axis=-1).reshape(-1,3).astype(float)
    r2=np.sum(K*K,axis=1); m=(r2>0)&(r2<=N*N); K=K[m]; r2=r2[m]
    c=chi(np.sqrt(r2)/N)
    Pv = v0[None,:] - K*(K@v0/r2)[:,None]
    coef=(c/r2)[:,None]*Pv
    ph=np.exp(1j*(K@y)/N)
    return (coef*ph[:,None]).sum(0).real/N
y=np.array([1.0,0.5,-0.3])
prev=None
for N in [8,16,32,64]:
    val=uN_over_N(N,y)
    print(N, np.round(val,7), "" if prev is None else "diff %.3e"%np.linalg.norm(val-prev))
    prev=val
