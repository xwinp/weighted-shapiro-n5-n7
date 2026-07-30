#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Enumerate ALL positive S1 (one-zero {0}, support {1..6}) stationary points
of n=7 P across the FULL p in (0,1), recording P and Morse index (reduced
Hessian eigenvalues on the support tangent).  Purpose: map the out-of-band
structure that §4.3 must cover (GPT review: prove p notin (a7,b7) => inf_S1 P>=7).

Multi-start root-find of grad P=0 (Euler lambda=0) at each p; dedup; record
(P-7, Morse=#neg).  Output the p-range and P-range of every distinct branch.
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
            xp[i]+=h;xp[j]+=h;xm[i]-=h;xm[j]-=h;xpp[i]+=h;xpp[j]-=h;xmm[i]-=h;xmm[j]+=h
            H[i,j]=(Pval(xp,p)-Pval(xm,p)-Pval(xpp,p)+Pval(xmm,p))/(4*h*h)
    return H
def tangent_basis(free_indices):
    k=len(free_indices); B=np.zeros((n,k-1))
    for j in range(k-1):
        B[free_indices[j],j]=1; B[free_indices[-1],j]=-1
    return B
def reduced_eigs(x,p,free_indices):
    H=hessian_full(x,p); B=tangent_basis(free_indices)
    return np.linalg.eigvalsh(B.T@H@B)

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

a7,b7=0.214273520909841,0.328627677916592
# multi-start inits spanning plausible stationary geometries
inits=[]
rng=np.random.default_rng(0)
for _ in range(60):
    inits.append(rng.uniform(0.05,2.5,size=5))
# plus a few structured inits near known in-band point and near uniform
inits.append([0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029])
inits+=[ [0.3,0.7,1.5,0.07,1.3],[0.5,0.5,1.0,0.5,0.5],[1,1,1,1,1],[0.2,0.5,2.0,0.05,1.5],[0.8,1.2,0.8,0.8,1.2] ]

rows=[]
ps=np.linspace(0.02,0.98,49)
for pp in ps:
    found=[]
    for init in inits:
        x=solve_S1(pp,init)
        if x is None: continue
        xs=x/x.sum()
        Pv=Pval(xs,pp)
        # dedup
        if any(abs(Pv-f[0])<1e-5 for f in found): continue
        eigs=reduced_eigs(xs,pp,sup)
        neg=int((eigs<-1e-3).sum())
        found.append((Pv,neg,xs))
    for (Pv,neg,xs) in found:
        rows.append((pp,Pv,Pv-7,neg))

# print branch summary: group by Morse and by p-region
print("p-range scan: %d stationary (p,P,P-7,Morse) entries"%len(rows))
print("in-band (a7,b7):")
for r in rows:
    if a7<r[0]<b7: print("  p=%.3f P=%.5f P-7=%+.5f Morse=%d"%(r[0],r[1],r[2],r[3]))
print("out-of-band p<a7 or p>b7:")
for r in rows:
    if not(a7<r[0]<b7): print("  p=%.3f P=%.5f P-7=%+.5f Morse=%d"%(r[0],r[1],r[2],r[3]))

# what is min P out-of-band?
ob=[r for r in rows if not(a7<r[0]<b7)]
if ob:
    mp=min(ob,key=lambda r:r[1])
    print("\nOUT-OF-BAND min P = %.6f at p=%.3f (Morse=%d)  -> %s"%(
        mp[1],mp[0],mp[3], "P>=7 OK" if mp[1]>=7-1e-6 else "P<7 GAP!!"))
print("DONE")
