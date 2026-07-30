#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Land on H_C stationary points exactly by 1D root-find of the free height rho1.

Given (u,v,w,z)=beta1..4 and p: rho2..rho5 are fixed via rho_{i+1}=(p/q)beta_i/(1-beta_i).
With x1=1 gauge, x=(0,1,rho1, rho1*rho2, rho1*rho2*rho3, rho1*rho2*rho3*rho4,
rho1*rho2*rho3*rho4*rho5). The single remaining KKT equation g1=rho1*dP/drho1=0
determines rho1.  Solve 1D, then verify all 5 KKT residuals ~0, evaluate P and
CORRECT reduced-Hessian Morse index.
"""
import numpy as np
from scipy.optimize import brentq, root

n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s
def grad_all(x,p):
    free=[2,3,4,5,6]; g=[]
    h=1e-7; f0=Pval(x,p)
    for idx in free:
        xp=x.copy(); xp[idx]+=h; g.append((Pval(xp,p)-f0)/h)
    return g
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
            H[i,j]=H[j,i]=(Pv(v0+ei+ej)-Pv(v0+ei-ej)-Pv(v0-ei+ej)-Pv(v0-ei-ej)*0+Pv(v0-ei-ej)*0-0) if False else (Pv(v0+ei+ej)-Pv(v0+ei-ej)-Pv(v0-ei+ej)+Pv(v0-ei-ej))/(4*h**2)
    return np.linalg.eigvalsh(H)

def x_from_rho1(rho1, rhos, p):
    x=np.zeros(n); x[1]=1.0
    cum=1.0
    for i in range(5):
        cum*=rhos[i] if i==0 else rhos[i]
    # build: x2=rho1, x3=rho1*rho2, ...
    x[2]=rho1
    c=rho1
    for i in range(1,5):
        c=c*rhos[i]; x[2+i]=c
    return x

hc_pts=[
 (0.85,0.33640,0.62826,0.86857,0.38982),(0.90,0.22413,0.69450,0.72280,0.38853),(0.95,0.11214,0.75198,0.57853,0.39000),
 (0.85,0.07495,0.60588,0.42797,0.53680),(0.90,0.05188,0.63532,0.42552,0.52495),(0.95,0.02686,0.66406,0.41920,0.51395),
]
print("H_C exact reconstruction (1D rho1 root-find)")
print("z      w       p       rho1    P        P-7     Morse  KKTmax  eigs")
for (zv,wv,vv,uu,pp) in hc_pts:
    q=1-pp
    rhos=[None, (pp/q)*uu/(1-uu), (pp/q)*vv/(1-vv), (pp/q)*wv/(1-wv), (pp/q)*zv/(1-zv)]
    def g1(r1):
        x=np.zeros(n); x[1]=1.0; x[2]=r1; c=r1
        for i in range(1,5): c=c*rhos[i]; x[2+i]=c
        h=1e-7; return (Pval(x+np.eye(1,n,2)[0]*h,pp)-Pval(x,pp))/h  # dP/dx2 ~ g1 (rho1=x2)
    # scan rho1 for sign change of g1
    r1_lo=1e-3; r1_hi=50.0
    # find bracket
    xs=np.linspace(1e-3,30,2000)
    gs=np.array([g1(r) for r in xs])
    found=False
    for i in range(len(xs)-1):
        if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and gs[i]*gs[i+1]<0:
            try:
                r1=brentq(g1,xs[i],xs[i+1],xtol=1e-14)
                found=True; break
            except: pass
    if not found:
        print("%.2f %.5f %.5f  no rho1 bracket"%(zv,wv,pp)); continue
    x=np.zeros(n); x[1]=1.0; x[2]=r1; c=r1
    for i in range(1,5): c=c*rhos[i]; x[2+i]=c
    xs_=x/x.sum()
    Pv=Pval(xs_,pp)
    gr=grad_all(x,pp)  # un-normalized x; grad scale-equivariant? P degree0 so dP/dx scale-equiv
    kkt=max(abs(g) for g in gr)
    eigs=reduced_eigs(xs_,pp,1e-5)
    neg=int((eigs<-1e-3).sum())
    print("%.2f %.5f %.5f %.5f  %.6f %+.6f  %s  %.1e  %s"%(
        zv,wv,pp,r1,Pv,Pv-7,"MIN" if neg==0 else "SADDLE%d"%neg,kkt,np.array2string(eigs,precision=2)))
print("DONE")
