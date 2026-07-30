#!/usr/bin/env python3
"""
S1 (1-zero {0}) Morse-index analysis across the band.
Trace the S1 positive stationary branch (support {1..6}, x1=1, vars b,c,d,e,f,
grad P = 0 via Euler lambda=0) across p in (a7,b7). Rescale to simplex sum=1.
Compute the reduced Hessian: full 7x7 numerical Hessian of P projected onto the
simplex-face tangent V = {v: v0=0, sum_{1..6} v_i = 0} (dim 5). Eigenvalues ->
Morse index (#negative). If >=1 negative everywhere in band -> S1 stationary
points are saddles (not local mins).
Also do S0 (full support) for comparison.
"""
import numpy as np
from scipy.optimize import root
n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s

def hessian_full(x,p):
    H=np.zeros((n,n)); h=1e-5
    for i in range(n):
        for j in range(n):
            xp=x.copy();xm=x.copy();xpp=x.copy();xmm=x.copy()
            xp[i]+=h;xp[j]+=h;xm[i]-=h;xm[j]-=h
            xpp[i]+=h;xpp[j]-=h;xmm[i]-=h;xmm[j]+=h
            H[i,j]=(Pval(xp,p)-Pval(xm,p)-Pval(xpp,p)+Pval(xmm,p))/(4*h*h)
    return H

# tangent basis for V={v0=0, sum_{1..6}v=0}: use e_i - e_6 for i=1..5 (indices 1..5), v0=0
def tangent_basis(free_indices):
    # free_indices = support indices (v0=0 fixed if zero). constraint sum over support=0.
    # basis vectors on the support subspace summing to 0
    k=len(free_indices)
    B=np.zeros((n,k-1))
    for j in range(k-1):
        B[free_indices[j],j]=1
        B[free_indices[-1],j]=-1
    return B

def reduced_eigs(x, p, free_indices):
    H=hessian_full(x,p)
    B=tangent_basis(free_indices)
    return np.linalg.eigvalsh(B.T@H@B)

# ---- S1: zero {0}, support [1,2,3,4,5,6], x1=1, vars x2..x6 ----
def solve_S1(p, init):
    sup=[1,2,3,4,5,6]; free=[2,3,4,5,6]
    def grad(v):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(free)):
            xp=x.copy(); xp[free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14})
    if r.success and max(abs(r.fun))<1e-7 and all(r.x>1e-8):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x, r.x
    return None,None

print("=== S1 branch across band (a7,b7) ===")
a7,b7=0.214273520909841,0.328627677916592
prev=None
for pp in np.linspace(a7+0.005, b7-0.005, 12):
    init = prev if prev is not None else [0.27,0.68,1.55,0.066,1.30]
    x,v = solve_S1(pp, init)
    if x is None:
        # try GPT init
        x,v = solve_S1(pp, [0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029])
    if x is None:
        print(f"  p={pp:.4f} no convergence"); continue
    prev=v
    xs=x/x.sum()  # rescale to simplex
    Pv=Pval(xs,pp)
    eigs=reduced_eigs(xs,pp,[1,2,3,4,5,6])
    neg=int((eigs<-1e-4).sum())
    print(f"  p={pp:.4f} P={Pv:.6f} (P-7={Pv-7:+.5f}) reduced-Hess eigs={np.round(eigs,3)} neg={neg} -> {'SADDLE' if neg>0 else 'min?'}")

print("\n=== S0 branch (full support) across band ===")
def solve_S0(p, init):
    free=[1,2,3,4,5,6]
    def grad(v):
        x=np.zeros(n); x[0]=1
        for j,idx in enumerate(free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(free)):
            xp=x.copy(); xp[free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14})
    if r.success and max(abs(r.fun))<1e-7 and all(r.x>1e-8):
        x=np.zeros(n); x[0]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x, r.x
    return None,None
prev0=None
for pp in [0.22,0.25,0.27,0.30,0.329]:
    init=prev0 if prev0 is not None else [0.2899598915706492,0.7447488449365604,0.8875842694414185,0.3764381331194338,1.1665179227192974,0.1591166088535238]
    x,v=solve_S0(pp,init)
    if x is None: print(f"  p={pp:.4f} no convergence"); continue
    prev0=v
    xs=x/x.sum(); Pv=Pval(xs,pp)
    eigs=reduced_eigs(xs,pp,[0,1,2,3,4,5,6])
    neg=int((eigs<-1e-4).sum())
    print(f"  p={pp:.4f} P={Pv:.6f} (P-7={Pv-7:+.5f}) reduced-Hess eigs={np.round(eigs,3)} neg={neg} -> {'SADDLE' if neg>0 else 'min?'}")
