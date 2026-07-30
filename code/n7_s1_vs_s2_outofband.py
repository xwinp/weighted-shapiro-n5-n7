#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test option (C): is inf S2 <= P_{S1}^stat for p>b7 (out-of-band high)?
If so, the S1 global min out-of-band high-p is controlled by the S2 boundary
(>=7), sidestepping the P_{S1}^stat>7 crossing resultant.

S2 = zero-set {0,2}, support {1,3,4,5,6}.  Compare S2 stationary value vs S1
stationary value across out-of-band p.  Also S3 for completeness.
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

def solve_face(p, zero_set, init):
    """zero_set: list of indices forced to 0. x normalized x1=1 if 1 not in zero_set
       else pick first nonzero support index =1."""
    support=[i for i in range(n) if i not in zero_set]
    # fix scale: set the first support var =1
    fixed=support[0]
    freev=support[1:]
    def grad(v):
        x=np.zeros(n); x[fixed]=1.0
        for j,idx in enumerate(freev): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(freev)):
            xp=x.copy(); xp[freev[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad,init,method='hybr',options={'xtol':1e-14})
    if r.success and max(abs(r.fun))<1e-7 and all(r.x>1e-8):
        x=np.zeros(n); x[fixed]=1.0
        for j,idx in enumerate(freev): x[idx]=r.x[j]
        return x
    return None

a7,b7=0.214273520909841,0.328627677916592
# S1 init (4 free after x1=1): b,c,d,e,f
s1init=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
print("p     S1_stat   S2_stat   S3_stat   S1-S2   S1-S3   region")
for pp in [0.10,0.20,0.25,0.30,0.34,0.40,0.45,0.50,0.60,0.70,0.80,0.90]:
    # S1: zero {0}, support {1,2,3,4,5,6}, fixed=1, free=[2,3,4,5,6] (5 vars)
    x1=solve_face(pp,[0],s1init)
    P1=Pval(x1/x1.sum(),pp) if x1 is not None else float('nan')
    if x1 is not None:
        xs=x1/x1.sum(); s1init=[xs[2],xs[3],xs[4],xs[5],xs[6]]
    # S2: zero {0,2}, support {1,3,4,5,6}, fixed=1, free=[3,4,5,6] (4 vars)
    best2=None
    rng=np.random.default_rng(int(pp*999))
    for init2 in [[0.7,1.5,0.07,1.3]]+[list(rng.uniform(0.05,2.5,size=4)) for _ in range(40)]:
        x2=solve_face(pp,[0,2],init2)
        if x2 is None: continue
        Pv=Pval(x2/x2.sum(),pp)
        if best2 is None or Pv<best2: best2=Pv
    # S3: zero {0,3}, support {1,2,4,5,6}, fixed=1, free=[2,4,5,6] (4 vars)
    best3=None
    for init3 in [[0.7,1.5,0.07,1.3]]+[list(rng.uniform(0.05,2.5,size=4)) for _ in range(40)]:
        x3=solve_face(pp,[0,3],init3)
        if x3 is None: continue
        Pv=Pval(x3/x3.sum(),pp)
        if best3 is None or Pv<best3: best3=Pv
    reg="IN" if a7<pp<b7 else "out"
    print("%.3f  %8.5f  %8.5f  %8.5f  %+7.4f  %+7.4f  %s"%(pp,P1,best2 if best2 else 0,best3 if best3 else 0,P1-(best2 or 0),P1-(best3 or 0),reg))
print("DONE")
