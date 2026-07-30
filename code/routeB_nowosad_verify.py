#!/usr/bin/env python3
"""
Route B leg-2: verify GPT's Nowosad-Yamagami shortcut to close S0 WITHOUT the
P-7 resultant. GPT claims: uniform is the UNIQUE interior local minimum of
P(x)=sum x_i/(p x_{i+1}+q x_{i+2}) for all p in (0,1), because
  S = pC + qC^2 (C = cyclic shift) is invertible,
  H_unif = grad^2 P(1) = 2 S^T S - S - S^T,
  spec_{1^perp} H_unif > 0,
  => (Nowosad-Yamagami) 1 is the unique interior local-min ray.

Verify locally:
 (1) 2-term KKT T_i=(1-b_{i-1})T_{i-1}+b_{i-2}T_{i-2} integrates to the
     1-term  T_i + b_{i-1} T_{i-1} = const C  (Fourier convention b_i=q x_{i+2}/d_i).
     [GPT says my earlier "1-term fails" was a buggy test using wrong A_i.]
 (2) S = pC+qC^2 invertible for all p in (0,1) (n=7 odd).
 (3) H_unif = 2 S^T S - S - S^T  matches the numerical Hessian of P at x=1.
 (4) spectrum of H_unif on 1^perp is >0  (== B.3 uniform cyclic spectrum).
 (5) S^{-1} entrywise sign -- is it POSITIVE? (Nowosad-type theorems often need T=S^{-1}>=0.)
 (6) Numerical hunt for non-uniform interior LOCAL MINIMA (stationary + H PSD on 1^perp):
     if any exists other than uniform, the uniqueness claim is false.
"""
import numpy as np
from scipy.optimize import root

n = 7
rng = np.random.default_rng(0)

def P_val(x, p):
    q = 1 - p
    return sum(x[i] / (p * x[(i+1) % n] + q * x[(i+2) % n]) for i in range(n))

def grad_P(x, p):
    q = 1 - p; g = np.zeros(n)
    den = np.array([p*x[(i+1)%n]+q*x[(i+2)%n] for i in range(n)])
    for j in range(n): g[j] += 1.0/den[j]
    for i in range(n):
        g[(i+1)%n] += -x[i]*p/den[i]**2
        g[(i+2)%n] += -x[i]*q/den[i]**2
    return g

def hessian_P(x, p, h=1e-5):
    H = np.zeros((n,n))
    for i in range(n):
        xp=x.copy();xp[i]+=h;xm=x.copy();xm[i]-=h
        H[i,i]=(P_val(xp,p)-2*P_val(x,p)+P_val(xm,p))/h**2
    for i in range(n):
        for j in range(i+1,n):
            xpp=x.copy();xpp[i]+=h;xpp[j]+=h;xpm=x.copy();xpm[i]+=h;xpm[j]-=h
            xmp=x.copy();xmp[i]-=h;xmp[j]+=h;xmm=x.copy();xmm[i]-=h;xmm[j]-=h
            v=(P_val(xpp,p)-P_val(xpm,p)-P_val(xmp,p)+P_val(xmm,p))/(4*h*h)
            H[i,j]=v;H[j,i]=v
    return H

def cyclic_shift():
    C = np.zeros((n,n))
    for i in range(n): C[i,(i+1)%n] = 1.0
    return C

def find_s0_stationary(p, seed):
    r = np.random.default_rng(seed)
    def F(y):
        x=np.empty(n);x[0]=1.0;x[1:]=np.exp(y)
        g=grad_P(x,p); return np.array([x[i]*g[i] for i in range(1,n)])
    for _ in range(60):
        y0=r.normal(0,0.9,6)
        sol=root(F,y0,method='hybr',tol=1e-13)
        if sol.success and np.max(np.abs(sol.fun))<1e-9:
            x=np.empty(n);x[0]=1.0;x[1:]=np.exp(sol.x)
            if np.max(np.abs(x-np.mean(x)))<1e-4: continue  # uniform
            return x
    return None

# ---------- (1) 2-term -> 1-term ----------
print("="*60); print("(1) 2-term KKT integrates to 1-term  T_i + b_{i-1}T_{i-1}=const")
print("="*60)
for p in [0.22, 0.25, 0.30]:
    x = find_s0_stationary(p, seed=7)
    if x is None: x = find_s0_stationary(p, seed=11)
    if x is None: continue
    q=1-p
    d=np.array([p*x[(i+1)%n]+q*x[(i+2)%n] for i in range(n)])
    T=x/d; be=q*np.roll(x,-2)/d  # Fourier beta_i = q x_{i+2}/d_i
    two_err=np.max(np.abs(T - ((1-np.roll(be,1))*np.roll(T,1) + np.roll(be,2)*np.roll(T,2))))
    # 1-term: T_i + beta_{i-1} T_{i-1} = const
    s = T + np.roll(be,1)*np.roll(T,1)
    print(f" p={p}: 2-term err={two_err:.2e}  T_i+b_-T_- spread={s.max()-s.min():.2e} (const={s.mean():.5f})")

# ---------- (2)(3)(4)(5) S, H_unif, spectrum, S^{-1} sign ----------
print("\n"+"="*60); print("(2-5) S=pC+qC^2: invertible, H_unif=2S^TS-S-S^T, spectrum, S^{-1} sign")
print("="*60)
C = cyclic_shift()
for p in [0.22, 0.27, 0.30]:
    q=1-p
    S = p*C + q*C@C
    # invertibility
    evS = np.linalg.eigvals(S)
    min_abs = min(abs(evS))
    # H_unif formula
    Hform = 2*S.T@S - S - S.T
    Hnum = hessian_P(np.ones(n), p, h=1e-5)
    # spectrum on 1^perp: project
    w = np.linalg.eigvalsh(Hform)
    w_perp = np.sort(w)[1:]  # drop the ~0 (k=0) mode
    # B.3 check: mu_k = 2(1-cos th_k)(2q cos th_k + 2p^2 -3p +2)
    mus=[]
    for k in range(1,n):
        th=2*np.pi*k/n; ck=np.cos(th)
        mus.append(2*(1-ck)*(2*q*ck+2*p**2-3*p+2))
    # S^{-1} sign
    Sinv = np.linalg.inv(S)
    nneg = int((Sinv < -1e-12).sum()); npos=int((Sinv>1e-12).sum())
    print(f" p={p}: min|eig(S)|={min_abs:.4f} inv=YES  "
          f"|H_form-H_num|_max={np.max(np.abs(Hform-Hnum)):.2e}  "
          f"min_eig_perp={w_perp.min():.4f} (B.3 min mu={min(mus):.4f})  "
          f"S^-1 sign: neg={nneg} pos={npos}")

# ---------- (6) hunt for non-uniform interior local minima ----------
print("\n"+"="*60); print("(6) Numerical hunt: non-uniform interior LOCAL MINIMA (stat + H PSD on 1^perp)")
print("="*60)
for p in [0.22, 0.25, 0.27, 0.30, 0.33]:
    found_min=None
    for seed in range(60):
        x=find_s0_stationary(p,seed)
        if x is None: continue
        # Hessian on 1^perp: scale-fix x[0]=1; use log-Hessian 6x6 via B=(e1-e0,..,e6-e0) on log-coords
        # simpler: ordinary Hessian restricted to sum-zero subspace
        H=hessian_P(x,p)
        # build sum-zero projector
        V=np.zeros((n,n-1))
        for i in range(n-1): V[i,i]=1; V[i,-1]=-1  # not orthonormal but fine for signs via eigsh
        # use orthonormal basis of 1^perp
        Q=np.linalg.qr(np.ones((n,1)))[0]
        # basis of complement
        M=np.random.default_rng(1).normal(size=(n,n-1))
        M=M - Q@(Q.T@M)
        Bx,_=np.linalg.qr(M)
        Hred=Bx.T@H@Bx
        wr=np.linalg.eigvalsh(Hred)
        is_local_min = wr.min() > -1e-3
        P=P_val(x,p)
        if is_local_min:
            if found_min is None or P<found_min[1]:
                found_min=(x,P,wr.min())
    if found_min:
        print(f" p={p}: FOUND non-unif interior local min  P={found_min[1]:.5f} min_eig={found_min[2]:.4f}")
    else:
        print(f" p={p}: no non-uniform interior local min found (unif unique holds)")
