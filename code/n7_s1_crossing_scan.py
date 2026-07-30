#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Locate precisely where the S1 interior stationary value P_{S1}^stat crosses 7,
relative to the S2 band edges a7,b7.  Decisive for the out-of-band theorem:
if P_{S1}^stat crosses 7 INSIDE [a7,b7], then out-of-band P>=7 (clean).
If it crosses OUTSIDE the band, we must check the Morse index there.

S1 = one-zero face {0}, support {1..6}, x1=1 fixed (degree-0 normalization).
Solve grad_{b,c,d,e,f} P = 0; record P and reduced-Hessian Morse index.
"""
import numpy as np
from scipy.optimize import root, brentq

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
# good seed near in-band branch
seed=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]

def Pstat_and_morse(p):
    x=solve_S1(p,seed)
    if x is None:
        # try a few random inits
        rng=np.random.default_rng(int(p*100000))
        for _ in range(20):
            x=solve_S1(p,rng.uniform(0.05,2.5,size=5))
            if x is not None: break
    if x is None: return None
    xs=x/x.sum()
    Pv=Pval(xs,p)
    eigs=reduced_eigs(xs,p,sup)
    neg=int((eigs<-1e-3).sum())
    return Pv,neg

print("Fine scan near a7=%.6f and b7=%.6f"%(a7,b7))
print("\n--- near a7 (crossing low side) ---")
for p in np.linspace(0.15,0.24,37):
    r=Pstat_and_morse(p)
    if r: print("  p=%.4f  P=%.6f  P-7=%+.6f  Morse=%d  %s"%(p,r[0],r[0]-7,r[1],"IN" if a7<p<b7 else "out"))
print("\n--- near b7 (crossing high side) ---")
for p in np.linspace(0.30,0.40,41):
    r=Pstat_and_morse(p)
    if r: print("  p=%.4f  P=%.6f  P-7=%+.6f  Morse=%d  %s"%(p,r[0],r[0]-7,r[1],"IN" if a7<p<b7 else "out"))

# Refine the crossing near b7 by bisection on (P-7) sign, using continuation seed
print("\n--- refine P=7 crossing near b7 ---")
def Pm7(p):
    r=Pstat_and_morse(p)
    return r[0]-7 if r else None
# scan to bracket
ps=np.linspace(0.30,0.40,101); prev=None; brackets=[]
for p in ps:
    v=Pm7(p)
    if v is not None and prev is not None and prev[1] is not None:
        if (prev[1]<0)!=(v<0):
            brackets.append((prev[0],p))
    prev=(p,v)
print("  sign-change brackets:",brackets)
for (lo,hi) in brackets:
    try:
        pc=brentq(lambda p:Pm7(p), lo, hi, xtol=1e-10)
        print("  P=7 crossing at p=%.8f   b7=%.8f   diff=%+.8f  (b7-p=%+.8f)"%(pc,b7,pc-b7,b7-pc))
    except Exception as e:
        print("  brentq failed:",e)
print("DONE")
