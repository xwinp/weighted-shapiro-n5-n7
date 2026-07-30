#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CORRECT reduced Hessian for the S1 stationary branch, degree-0 gauge
(x0=0, x1=1 fixed; variables x2..x6).  Uses the PROPER second-derivative formula:
  diag:  [P(x+h ei) - 2P(x) + P(x-h ei)] / h^2
  off:   [P(++) - P(+-) - P(-+) + P(--)] / (4 h^2)
Reports P, Morse (#neg eigenvalues), trace, det, eigenvalues across p in (0,1).
Also prints the diagonal entries (which are individually >=0 since P is convex
in each x_i) so trace = sum diag + 0 (off-diag don't enter trace) >= 0 always --
a sanity check that the Hessian is now correct.
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

def Hred(x,p,h=1e-5):
    """5x5 gauge Hessian, vars indices [2,3,4,5,6]. Correct formula."""
    vi=[2,3,4,5,6]; k=5
    def Pp(d):
        xx=x.copy()
        for idx,sgn in d: xx[idx]+=sgn*h
        return Pval(xx,p)
    H=np.zeros((k,k))
    for a in range(k):
        for b in range(k):
            ia=vi[a]; ib=vi[b]
            if a==b:
                H[a,b]=(Pp([(ia,1)])-2*Pval(x,p)+Pp([(ia,-1)]))/(h*h)
            else:
                H[a,b]=(Pp([(ia,1),(ib,1)])-Pp([(ia,1),(ib,-1)])
                        -Pp([(ia,-1),(ib,1)])+Pp([(ia,-1),(ib,-1)]))/(4*h*h)
    return H

free=[2,3,4,5,6]
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
seed=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
print("p      P        P-7       Morse  trace      det        diag-sum  region")
alltrace_pos=True
for pp in np.linspace(0.12,0.78,34):
    x=solve_S1(pp,seed)
    if x is None:
        rng=np.random.default_rng(int(pp*1e6))
        for _ in range(40):
            x=solve_S1(pp,rng.uniform(0.05,2.5,size=5))
            if x is not None: break
    if x is None:
        print("%5.3f  (no conv)"%pp); continue
    xs=x/x.sum()
    Pv=Pval(xs,pp)
    H=Hred(xs,pp,h=1e-5)
    eigs=np.linalg.eigvalsh(H)
    neg=int((eigs<-1e-2).sum())
    tr=np.trace(H); det=np.linalg.det(H); dsum=np.sum(np.diag(H))
    reg="IN" if a7<pp<b7 else "out"
    if tr<0: alltrace_pos=False
    print("%5.3f  %8.5f  %+8.5f  %d   %+10.3f  %+11.4e  %+10.3f  %s"%(pp,Pv,Pv-7,neg,tr,det,dsum,reg))
    seed=[xs[2],xs[3],xs[4],xs[5],xs[6]]
print("\ntrace>=0 everywhere (as expected, sum of convex diagonals)? ", alltrace_pos)
print("DONE")
