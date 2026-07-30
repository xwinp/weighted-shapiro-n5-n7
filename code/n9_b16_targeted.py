#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Targeted (B.16) hunt: seed near the known palindromic saddle with ANTI-symmetric
perturbations (breaking the palindrome) + dense random, to specifically surface any
non-palindromic stationary point. If none with P<9 & PSD, strengthens (B.16)."""
import numpy as np
from scipy.optimize import minimize, root
import sys, warnings
warnings.filterwarnings('ignore')

def Px(x,p):
    q=1-p; return sum(x[i]/(p*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9))
def grad_P(x,p):
    q=1-p; g=np.zeros(9)
    for i in range(9):
        Dj=p*x[(i+1)%9]+q*x[(i+2)%9]; g[i]+=1.0/Dj
        for l in range(9):
            Dl=p*x[(l+1)%9]+q*x[(l+2)%9]
            if (l+1)%9==i: g[i]-=p*x[l]/Dl**2
            if (l+2)%9==i: g[i]-=q*x[l]/Dl**2
    return g
def polish(x,p):
    def F(v):
        xx=np.array([0.0,1.0]+list(v)); g=grad_P(xx,p); return g[2:]
    sol=root(F, x[2:9], method='hybr', options={'xtol':1e-14,'maxfev':8000})
    if sol.success:
        xx=np.array([0.0,1.0]+list(sol.x)); gn=np.sum(grad_P(xx,p)[2:]**2)
        if gn<1e-20: return xx,gn
    return None,None
def hess7(x,p):
    def P7(v):
        xx=np.array([0.0,1.0]+list(v)); return Px(xx,p)
    v0=x[2:9]; h=1e-5; H=np.zeros((7,7))
    for i in range(7):
        for j in range(i,7):
            ei=np.zeros(7);ej=np.zeros(7);ei[i]=h;ej[j]=h
            if i==j: H[i,j]=(P7(v0+ei)-2*P7(v0)+P7(v0-ei))/(h*h)
            else:
                v=(P7(v0+ei+ej)-P7(v0+ei-ej)-P7(v0-ei+ej)+P7(v0-ei-ej))/(4*h*h); H[i,j]=H[j,i]=v
    return H

def find_pal_saddle(p, seed=1):
    def F(v):
        x=np.array([0.0]+list(np.exp(v))); g=grad_P(x,p); return np.sum(g[1:]**2)
    rng=np.random.default_rng(seed)
    for _ in range(300):
        v0=np.log(rng.dirichlet(np.ones(8))+1e-9)
        r=minimize(F,v0,method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-11,'fatol':1e-15})
        if r.fun<1e-12:
            x=np.array([0.0]+list(np.exp(r.x))); xp,gn=polish(x,p)
            if xp is not None and 8<Px(xp,p)<9: return xp
    return None

print("="*70); print("Targeted (B.16) non-palindromic hunt (antisymmetric + dense starts)"); print("="*70, flush=True)
for p in [0.40, 0.43]:
    xs=find_pal_saddle(p)
    if xs is None: print(f"p={p}: no palindromic saddle found"); continue
    # normalize x[1]=1
    xs=xs/xs[1]
    print(f"\np={p}: palindromic saddle P={Px(xs,p):.5f}", flush=True)
    rng=np.random.default_rng(42); found=[]; ntry=0
    # antisymmetric perturbation starts: x* * exp(eta), eta antisymmetric (eta[i]=-eta[9-i])
    for _ in range(300):
        eta=np.zeros(9); e=rng.normal(0,0.15,4); eta[1:5]=e; eta[5:9]=-e[::-1]; eta[0]=0
        x0=xs*np.exp(eta); x0=x0/x0[1]
        xp,gn=polish(x0,p)
        if xp is None: continue
        ntry+=1
        sc=xp[1:].max()
        if (xp[1:]/sc).min()<1e-3: continue
        P=Px(xp,p)
        if not (5<P<12): continue
        q=1-p; T=np.array([xp[i]/(p*xp[(i+1)%9]+q*xp[(i+2)%9]) for i in range(9)]); Tn=T[1:]/T[1]
        pal=np.max(np.abs(Tn-Tn[::-1]))<0.02
        key=round(P,3)
        if not any(abs(key-f[0])<0.01 for f in found):
            found.append((key,P,xp,pal,gn))
    # also 300 dense random
    for _ in range(300):
        v0=np.log(rng.dirichlet(np.ones(8))+1e-9)
        x=np.array([0.0]+list(np.exp(v0)))
        def F(v):
            xx=np.array([0.0]+list(np.exp(v))); g=grad_P(xx,p); return np.sum(g[1:]**2)
        r=minimize(F,v0,method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-10,'fatol':1e-14})
        if r.fun<1e-10:
            x=np.array([0.0]+list(np.exp(r.x))); xp,gn=polish(x,p)
            if xp is None: continue
            ntry+=1; sc=xp[1:].max()
            if (xp[1:]/sc).min()<1e-3: continue
            P=Px(xp,p)
            if not (5<P<12): continue
            q=1-p; T=np.array([xp[i]/(p*xp[(i+1)%9]+q*xp[(i+2)%9]) for i in range(9)]); Tn=T[1:]/T[1]
            pal=np.max(np.abs(Tn-Tn[::-1]))<0.02
            key=round(P,3)
            if not any(abs(key-f[0])<0.01 for f in found):
                found.append((key,P,xp,pal,gn))
    print(f"  polished stationary points found: {len(found)} (from {ntry} convergences)", flush=True)
    viol=0
    for key,P,xp,pal,gn in sorted(found):
        H=hess7(xp,p); ev=np.sort(np.linalg.eigvalsh(H)); psd=ev.min()>-1e-2; morse=int((ev<-1e-2).sum())
        flag=''
        if P<9 and psd: flag=' *** B.16 VIOLATION ***'; viol+=1
        elif P<9 and not psd: flag=' (P<9 saddle)'
        else: flag=' (P>=9)'
        print(f"    P={P:.5f} pal={pal} morse={morse} eig_min={ev[0]:.3f}{flag}", flush=True)
    if viol==0: print(f"  -> p={p}: NO non-palindromic P<9 local minimum. Supports (B.16).", flush=True)
print("\nDone.", flush=True)
