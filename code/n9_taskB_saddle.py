#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Task B: find the INTERIOR saddle of {0}-face via ∇P=0 (reflection-symmetric),
confirm det H_red<0 when P<9, and P=9 crossing p_P < b9."""
import numpy as np
from scipy.optimize import root, brentq
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

b9 = 0.43388588203369836

def Px(x, p):
    q=1-p
    return sum(x[i]/(p*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9))

def grad_P(x, p):
    q=1-p
    g=np.zeros(9)
    for i in range(9):
        j=(i+1)%9; k=(i+2)%9
        Dj=p*x[j]+q*x[k]
        g[i]+=1.0/Dj
        # d/dx_m of x_l/(p x_{l+1}+q x_{l+2})
        for l in range(9):
            Dl=p*x[(l+1)%9]+q*x[(l+2)%9]
            if (l+1)%9==i: g[i]-=p*x[l]/Dl**2
            if (l+2)%9==i: g[i]-=q*x[l]/Dl**2
    return g

def saddle(p, nstarts=200, seed=1):
    """Find reflection-symmetric interior stationary x=(0,a,b,c,d,d,c,b,a)."""
    # free params (b/a,c/a,d/a); solve grad on a,b,c,d with a=1 gauge.
    rng=np.random.default_rng(seed)
    cands=[]
    for _ in range(nstarts):
        v0=rng.uniform(0.05,3.0,3)
        def F(v):
            b,c,d=v
            x=np.array([0.0,1.0,b,c,d,d,c,b,1.0])
            g=grad_P(x,p)
            # independent residual eqns: g1-g8, g2-g7, g3-g6 (stationarity up to gauge) + use g4
            # use g1,g2,g3,g4 projected: actually gauge fix a=1, stationarity ⟺ g on free vars {1,2,3,4}
            return [g[1],g[2],g[3]]   # 3 eqns, 3 unknowns (scale-invariance handles gauge)
        res=root(F,v0,method='hybr')
        if res.success and np.max(np.abs(res.fun))<1e-8:
            b,c,d=res.x
            x=np.array([0.0,1.0,b,c,d,d,c,b,1.0])
            if np.all(x>0) and b<50:
                cands.append((Px(x,p),x))
    if not cands: return None
    # pick the one with P closest to but distinct from boundary; return all unique
    cands.sort(key=lambda t:t[0])
    return cands

def hess_red4(x,p):
    """Hessian of P wrt free vars (b,c,d) with a=1 gauge: 3x3, plus check."""
    def P3(v):
        b,c,d=v
        xx=np.array([0.0,1.0,b,c,d,d,c,b,1.0])
        return Px(xx,p)
    v0=np.array([x[2],x[3],x[4]]); h=1e-5
    H=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            ei=np.zeros(3);ej=np.zeros(3);ei[i]=h;ej[j]=h
            if i==j: H[i,j]=(P3(v0+ei)-2*P3(v0)+P3(v0-ei))/(h*h)
            else: H[i,j]=(P3(v0+ei+ej)-P3(v0+ei-ej)-P3(v0-ei+ej)+P3(v0-ei-ej))/(4*h*h)
    return H

print("="*70); print("TASK B: interior saddle of {0}-face (∇P=0, reflection-symmetric)"); print("="*70)
# p=0.4
p=0.4
cands=saddle(p,nstarts=300,seed=7)
print(f"\np=0.4: distinct interior stationary found: {len(cands) if cands else 0}")
if cands:
    for Pm,x in cands[:5]:
        q=1-p
        T=np.array([x[i]/(p*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9)])
        seq=np.round(T[1:]/T[1],5)  # normalize by T1
        print(f"  P={Pm:.9f}  x=(0,{x[1]:.4f},{x[2]:.4f},{x[3]:.4f},{x[4]:.4f},{x[5]:.4f},{x[6]:.4f},{x[7]:.4f},{x[8]:.4f})")
        print(f"    T1..T8/T1 = {seq}  (compare 1,A,B,C,C,B,A,1)")
        H=hess_red4(x,p); ev=np.linalg.eigvalsh(H)
        print(f"    3x3 red Hessian eig: {np.round(np.sort(ev),5)}  det={np.linalg.det(H):.4e}  Morse_idx={int((ev<-1e-3).sum())}")

# scan p for P=9 crossing of the interior saddle (pick the saddle branch near P~8.8)
print("\nscan interior-saddle P vs p (pick branch), find P=9 crossing p_P:")
pts=[]
for p in np.linspace(0.10,0.435,30):
    cc=saddle(p,nstarts=120,seed=11)
    if not cc: continue
    # pick the stationary point with P in (7,9.5) closest to 9-region
    Pm,x=min(cc,key=lambda t:abs(t[0]-8.8))
    pts.append((p,Pm))
    print(f"  p={p:.4f}  P_saddle={Pm:.6f}")
# find crossing
ps=[t[0] for t in pts]; Ps=[t[1] for t in pts]
for i in range(len(ps)-1):
    if (Ps[i]-9)*(Ps[i+1]-9)<0:
        pP=brentq(lambda p: (lambda cc: min(cc,key=lambda t:abs(t[0]-8.8))[0] if cc else 9.0)(saddle(p,nstarts=120,seed=11))-9, ps[i],ps[i+1])
        print(f"  >> P=9 crossing at p_P = {pP:.10f}  (GPT 0.4318363763)  < b9? {pP<b9}")
        break
