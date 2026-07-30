#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Enumerate ALL positive S1 stationary points at each p (heavy multi-start) and
report P, CORRECT reduced-Hessian eigenvalues, Morse index, det, beta (u,v,w,z).

Uses the CORRECT second-derivative formula (the bug in n7_s1_branches_outofband.py
zeroed all off-diagonals: xpp was a copy of xp). Correct:
  diag:  [P(x+h e_i)-2P(x)+P(x-h e_i)]/h^2
  off:   [P(++)-P(+-)-P(-+)+P(--)]/(4 h^2)
Gauge: x0=0 fixed, x1=1 fixed (scale), free=[2,3,4,5,6], reduced Hessian on the
5 free vars (tangent to the gauge-fixing x1=const; the scale direction is quotiented
by degree-0 homogeneity so 5x5 is correct).
Three step sizes h=1e-4,1e-5,1e-6 for stability check.
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
    """5x5 reduced Hessian on free vars [2,3,4,5,6], x1=1 gauge."""
    free=[2,3,4,5,6]; k=5
    def Pv(v):
        xx=np.zeros(n); xx[1]=1.0
        for j,idx in enumerate(free): xx[idx]=v[j]
        return Pval(xx,p)
    v0=np.array([x[idx] for idx in free])
    H=np.zeros((k,k))
    for i in range(k):
        ei=np.zeros(k); ei[i]=h
        H[i,i]=(Pv(v0+ei)-2*Pv(v0)+Pv(v0-ei))/h**2
    for i in range(k):
        for j in range(i+1,k):
            ei=np.zeros(k); ej=np.zeros(k); ei[i]=h; ej[j]=h
            cpp=Pv(v0+ei+ej); cpm=Pv(v0+ei-ej); cmp=Pv(v0-ei+ej); cmm=Pv(v0-ei-ej)
            H[i,j]=H[j,i]=(cpp-cpm-cmp+cmm)/(4*h**2)
    return np.linalg.eigvalsh(H), H

def solve_S1(p, init):
    free=[2,3,4,5,6]
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
    # beta_i = q*rho_{i+1}/(p+q*rho_{i+1}), i=1..4 -> rho2..rho5
    bet=[q*rh/(p+q*rh) for rh in rho[1:]]
    return bet[0],bet[1],bet[2],bet[3]   # u,v,w,z  (beta1,beta2,beta3,beta4)

a7,b7=0.214273520909841,0.328627677916592
p0=0.388528131361
print("Enumerating all positive S1 stationary points (CORRECT Hessian, 3 step sizes)")
print("p0=%.6f  band=(%.6f,%.6f)"%(p0,a7,b7))
for pp in [0.10,0.20,0.27,0.34,0.389,0.40,0.45,0.50,0.52,0.55,0.60,0.70,0.80,0.90]:
    found=[]
    rng=np.random.default_rng(0)
    inits=[rng.uniform(0.05,2.5,size=5) for _ in range(400)]
    inits.append([0.3,0.7,1.5,0.07,1.3]); inits.append([1,1,1,1,1])
    inits.append([0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029])
    for init in inits:
        x=solve_S1(pp,init)
        if x is None: continue
        xs=x/x.sum()
        Pv=Pval(xs,pp)
        if any(abs(Pv-f["P"])<1e-5 for f in found): continue
        eigs5={}
        for h in [1e-4,1e-5,1e-6]:
            eigs5[h],_=reduced_eigs(xs,pp,h)
        e1=eigs5[1e-5]
        neg=int((e1<-1e-3).sum())
        detM=float(np.linalg.det(np.diag(e1)))  # eigs product
        u,v,w,z=beta_coords(xs,pp)
        # H_B branch has a2=a5, a3=a4 symmetry -> x2/x1==x5/x6? check via beta
        # symmetric branch: u~?  just tag by z range
        hb = (z>0.5)  # crude
        found.append({"P":Pv})
        stab="MIN" if neg==0 else "SADDLE%d"%neg
        print("p=%.3f P=%.6f P-7=%+.6f %s det=%+.3e z=%.4f u,v,w,z=(%.4f,%.4f,%.4f,%.4f) eigs=%s"%(
            pp,Pv,Pv-7,stab,detM,z,u,v,w,z,np.array2string(e1,precision=2)))
        # check consistency across step sizes
        n4=int((eigs5[1e-4]<-1e-3).sum()); n6=int((eigs5[1e-6]<-1e-3).sum())
        if not (n4==neg==n6):
            print("   WARN step-size Morse mismatch: h=1e-4->%d, 1e-5->%d, 1e-6->%d"%(n4,neg,n6))
print("DONE")
