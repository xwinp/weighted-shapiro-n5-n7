#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Task B decisive: track {0}-face interior saddle branch by continuation,
compute Hessian det at p=0.4 (must be <0), find P=9 crossing p_P < b9."""
import numpy as np
from scipy.optimize import minimize, root
import sys, io, warnings
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

def Px(x,p):
    q=1-p
    return sum(x[i]/(p*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9))
def grad_P(x,p):
    q=1-p; g=np.zeros(9)
    for i in range(9):
        Dj=p*x[(i+1)%9]+q*x[(i+2)%9]; g[i]+=1.0/Dj
        for l in range(9):
            Dl=p*x[(l+1)%9]+q*x[(l+2)%9]
            if (l+1)%9==i: g[i]-=p*x[l]/Dl**2
            if (l+2)%9==i: g[i]-=q*x[l]/Dl**2
    return g

def refine(x0,p):
    """Newton-ish refine on ∇P=0 via root on log-coords (8 free, x0=0)."""
    def F(v):
        x=np.array([0.0]+list(np.exp(v))); g=grad_P(x,p)
        return g[1:]
    v0=np.log(np.clip(x0[1:],1e-12,None))
    r=root(F,v0,method='hybr',options={'xtol':1e-13,'maxfev':20000})
    x=np.array([0.0]+list(np.exp(r.x)))
    return x,r

def find_p04(seed=3,nstarts=600):
    def F(v):
        x=np.array([0.0]+list(np.exp(v))); g=grad_P(x,0.4); return np.sum(g[1:]**2)
    rng=np.random.default_rng(seed); best=None
    for _ in range(nstarts):
        v0=np.log(rng.dirichlet(np.ones(8))+1e-9)
        r=minimize(F,v0,method='Nelder-Mead',options={'maxiter':8000,'xatol':1e-12,'fatol':1e-16})
        if r.fun<1e-13:
            x=np.array([0.0]+list(np.exp(r.x))); P=Px(x,0.4)
            if 7<P<9.5 and (best is None or abs(P-8.83)<abs(best[0]-8.83)):
                best=(P,x)
    return best

def hess7(x,p):
    def P7(v):
        xx=np.array([0.0,1.0]+list(v)); return Px(xx,p)
    v0=x[2:9]; h=1e-4; H=np.zeros((7,7))
    for i in range(7):
        for j in range(7):
            ei=np.zeros(7);ej=np.zeros(7);ei[i]=h;ej[j]=h
            if i==j: H[i,j]=(P7(v0+ei)-2*P7(v0)+P7(v0-ei))/(h*h)
            else: H[i,j]=(P7(v0+ei+ej)-P7(v0+ei-ej)-P7(v0-ei+ej)+P7(v0-ei-ej))/(4*h*h)
    return H

b9=0.43388588203369836
print("="*70); print("TASK B decisive: {0}-face interior saddle"); print("="*70)
b=find_p04()
print(f"\np=0.4 saddle: P={b[0]:.9f}  (GPT 8.829694260862)")
x,_=refine(b[1],0.4); P=Px(x,0.4)
q=0.6; T=np.array([x[i]/(0.4*x[(i+1)%9]+0.6*x[(i+2)%9]) for i in range(9)])
seq=np.round(T[1:]/T[4],4)
print(f"  refined P={P:.9f}  T1..T8/T4={seq}  (1,A,B,C,C,B,A,1 with C=1)")
H=hess7(x,0.4); ev=np.sort(np.linalg.eigvalsh(H))
print(f"  7x7 red Hessian eig: {np.round(ev,5)}")
print(f"  det(H_red)={np.linalg.det(H):.4e}  Morse_idx={int((ev<-1e-3).sum())}  (GPT: saddle det<0)")

# continuation: warm-start from p=0.4 saddle, step p up
print("\nbranch continuation (warm-start):")
xw=x.copy(); pts=[]
for p in np.linspace(0.40,0.445,19):
    xr,rr=refine(xw,p)
    if rr.success and np.max(np.abs(rr.fun))<1e-7:
        P=Px(xr,p); xw=xr.copy(); pts.append((p,P))
        print(f"  p={p:.4f}  P_saddle={P:.7f}  {'<9' if P<9 else '>=9'}")
    else:
        print(f"  p={p:.4f}  refine failed")
# crossing
for i in range(len(pts)-1):
    if (pts[i][1]-9)*(pts[i+1][1]-9)<0:
        from scipy.optimize import brentq
        def gp(p):
            xr,rr=refine(xw if p<pts[i+1][0] else pts[i][1] and xw,p)
            # robust: warm from nearest stored
            return P-9
        # simpler linear
        t=(9-pts[i][1])/(pts[i+1][1]-pts[i][1])
        pP=pts[i][0]+t*(pts[i+1][0]-pts[i][0])
        print(f"\n  >> P=9 crossing p_P ≈ {pP:.7f}  (GPT 0.4318363763)  < b9={b9:.7f}? {pP<b9}")
        break
