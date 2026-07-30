#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reconstruct H_C-branch S1 stationary points from (u,v,w,z) = (beta1..beta4),
evaluate P and CORRECT reduced-Hessian Morse index.  H_C admissible lifts live at
p~p0 (branch 1, w~0.22) and p in [0.514,0.537] (branch 2, w~0.05) per routeA_verify_F.

Given (u,v,w,z,p): rho2=(p/q)u/(1-u), rho3=(p/q)v/(1-v), rho4=(p/q)w/(1-w),
rho5=(p/q)z/(1-z); rho1 (=x2) is the free height fixed by the closure/KKT. We seed
rho1=1 and refine with solve_S1 (which solves the full 5-var KKT), then verify the
refined point still has beta=(u,v,w,z) and lies on H_C.
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
def reduced_eigs(x,p,h):
    free=[2,3,4,5,6]; k=5
    def Pv(v):
        xx=np.zeros(n); xx[1]=1.0
        for j,idx in enumerate(free): xx[idx]=v[j]
        return Pval(xx,p)
    v0=np.array([x[idx] for idx in free]); H=np.zeros((k,k))
    for i in range(k):
        ei=np.zeros(k); ei[i]=h
        H[i,i]=(Pv(v0+ei)-2*Pv(v0)+Pv(v0-ei))/h**2
    for i in range(k):
        for j in range(i+1,k):
            ei=np.zeros(k); ej=np.zeros(k); ei[i]=h; ej[j]=h
            H[i,j]=H[j,i]=(Pv(v0+ei+ej)-Pv(v0+ei-ej)-Pv(v0-ei+ej)+Pv(v0-ei-ej))/(4*h**2)
    return np.linalg.eigvalsh(H)
def solve_S1(p, init):
    free=[2,3,4,5,6]
    def grad(v):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(free)):
            xp=x.copy(); xp[free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14,'maxfev':40000})
    if r.success and max(abs(r.fun))<1e-7 and all(r.x>1e-8):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x
    return None
def beta_of(x,p):
    b,c,d,e,f=x[2],x[3],x[4],x[5],x[6]; q=1-p
    rho=[b,c/b,d/c,e/d,f/e]
    return [q*rh/(p+q*rh) for rh in rho[1:]]  # beta1..4 = u,v,w,z

# H_C branch points from routeA_verify_F: (z, w, v, u, p)
# branch1 (p~p0):  (0.85,0.33640,0.62826,0.86857,0.38982),(0.90,0.22413,0.69450,0.72280,0.38853),(0.95,0.11214,0.75198,0.57853,0.39000)
# branch2 (p~0.52): (0.85,0.07495,0.60588,0.42797,0.53680),(0.90,0.05188,0.63532,0.42552,0.52495),(0.95,0.02686,0.66406,0.41920,0.51395)
hc_pts=[
 (0.85,0.33640,0.62826,0.86857,0.38982),(0.90,0.22413,0.69450,0.72280,0.38853),(0.95,0.11214,0.75198,0.57853,0.39000),
 (0.85,0.07495,0.60588,0.42797,0.53680),(0.90,0.05188,0.63532,0.42552,0.52495),(0.95,0.02686,0.66406,0.41920,0.51395),
]
print("H_C branch reconstruction: P and Morse index")
print("z      w       p       P        P-7     Morse  betas(u,v,w,z) match?")
for (zv,wv,vv,uu,pp) in hc_pts:
    q=1-pp
    # rho2..rho5 from beta=u,v,w,z
    rho2=(pp/q)*uu/(1-uu); rho3=(pp/q)*vv/(1-vv); rho4=(pp/q)*wv/(1-wv); rho5=(pp/q)*zv/(1-zv)
    # seed rho1=1 -> x=(0,1,1,rho2,rho2*rho3,rho2*rho3*rho4,rho2*rho3*rho4*rho5)
    x2=1.0; x3=rho2; x4=rho2*rho3; x5=rho2*rho3*rho4; x6=rho2*rho3*rho4*rho5
    init=[x2,x3,x4,x5,x6]
    # try a few rho1 scales
    x=None
    for s in [1.0,0.3,3.0,0.1,10.0]:
        x=solve_S1(pp,[s,s*rho2,s*rho2*rho3,s*rho2*rho3*rho4,s*rho2*rho3*rho4*rho5])
        if x is not None:
            xs=x/x.sum(); bt=beta_of(xs,pp)
            # check it matches (uu,vv,wv,zv) -> H_C point
            if max(abs(bt[0]-uu),abs(bt[1]-vv),abs(bt[2]-wv),abs(bt[3]-zv))<1e-3:
                break
            x=None
    if x is None:
        print("%.2f %.5f %.5f  -- did not converge to H_C point"%(zv,wv,pp)); continue
    xs=x/x.sum(); Pv=Pval(xs,pp)
    eigs=reduced_eigs(xs,pp,1e-5)
    neg=int((eigs<-1e-3).sum())
    bt=beta_of(xs,pp)
    match=max(abs(bt[0]-uu),abs(bt[1]-vv),abs(bt[2]-wv),abs(bt[3]-zv))
    print("%.2f %.5f %.5f  %.6f %+.6f  %s  (u,v,w,z)=(%.4f,%.4f,%.4f,%.4f) match=%.1e eigs=%s"%(
        zv,wv,pp,Pv,Pv-7,"MIN" if neg==0 else "SADDLE%d"%neg,bt[0],bt[1],bt[2],bt[3],match,np.array2string(eigs,precision=2)))
print("DONE")
