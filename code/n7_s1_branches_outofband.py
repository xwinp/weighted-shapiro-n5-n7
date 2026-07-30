#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reconcile the S1 branch structure out-of-band.  At several p (incl. p=0.5 in
the H_C range (0.389,0.537)) enumerate ALL positive S1 stationary points by heavy
multi-start, and for each report P, reduced-Hessian eigenvalues, Morse index, det,
and beta-coords (u,v,w,z).  Decides whether EVERY out-of-band S1 stationary point
is non-PSD (Morse>=1), and whether H_B / H_C branches are both saddles.
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
def hessian_full(x,p,h=1e-5):
    H=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            xp=x.copy();xm=x.copy();xpp=x.copy();xmm=x.copy()
            xp[i]+=h;xp[j]+=h;xm[i]-=h;xm[j]-=h;xpp[i]+=h;xpp[j]+=h;xmm[i]-=h;xmm[j]+=h
            H[i,j]=(Pval(xp,p)-Pval(xm,p)-Pval(xpp,p)+Pval(xmm,p))/(4*h*h)
    return H
def tangent_basis(free_indices):
    k=len(free_indices); B=np.zeros((n,k-1))
    for j in range(k-1):
        B[free_indices[j],j]=1.0
        B[free_indices[-1],j]=-1.0
    return B
def reduced_eigs(x,p,free_indices,h=1e-5):
    H=hessian_full(x,p,h); B=tangent_basis(free_indices)
    M=B.T@H@B
    return np.linalg.eigvalsh(M), M

sup=[1,2,3,4,5,6]; free=[2,3,4,5,6]
def solve_S1(p, init):
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
        return x
    return None

def beta_coords(x,p):
    b,c,d,e,f=x[2],x[3],x[4],x[5],x[6]
    rho=[b,c/b,d/c,e/d,f/e]; q=1-p
    bet=[q*rh/(p+q*rh) for rh in rho[1:]]
    return bet[0],bet[1],bet[2],bet[3]

for pp in [0.10,0.34,0.40,0.50,0.55,0.70,0.90]:
    print("\n=== p=%.2f ==="%pp)
    found=[]
    rng=np.random.default_rng(0)
    inits=[rng.uniform(0.05,2.5,size=5) for _ in range(200)]
    inits.append([0.3,0.7,1.5,0.07,1.3]); inits.append([1,1,1,1,1])
    inits.append([0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029])
    for init in inits:
        x=solve_S1(pp,init)
        if x is None: continue
        xs=x/x.sum()
        Pv=Pval(xs,pp)
        if any(abs(Pv-f["P"])<1e-5 for f in found): continue
        eigs,M=reduced_eigs(xs,pp,sup,h=1e-5)
        neg=int((eigs<-1e-3).sum())
        detM=np.linalg.det(M)
        u,v,w,z=beta_coords(xs,pp)
        found.append({"P":Pv})
        print("  P=%.6f P-7=%+.6f Morse=%d det=%+.4e eig=%s"%(Pv,Pv-7,neg,detM,np.array2string(eigs,precision=3)))
        print("     u,v,w,z=(%.5f,%.5f,%.5f,%.5f)  z-in-inband(0.9488,0.9911)? %s"%(u,v,w,z, 0.9488<z<0.9911))
print("DONE")
