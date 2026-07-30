#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""(B.16) numerical emptiness evidence: on the {0}-face (L=9), search ALL
stationary points (incl. non-palindromic) via full 8-var grad P=0; check whether
any has P<9 AND H_red PSD (local min, not saddle). If none across the band,
strong evidence for (B.16) emptiness (non-palindromic degenerate minima absent)."""
import numpy as np
from scipy.optimize import minimize, root
import sys, warnings
warnings.filterwarnings('ignore')

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

def polish(x,p):
    # refine to true stationary point via root on grad (gauge-fix x[1]=1)
    def F(v):
        xx=np.array([0.0,1.0]+list(v)); g=grad_P(xx,p); return g[2:]
    sol=root(F, x[2:9], method='hybr', options={'xtol':1e-14,'maxfev':5000})
    if sol.success:
        xx=np.array([0.0,1.0]+list(sol.x))
        gn=np.sum(grad_P(xx,p)[2:]**2)
        if gn<1e-18: return xx,gn
    return None,None

def find_stationary(p, nstarts=400, seed=1):
    def F(v):
        x=np.array([0.0]+list(np.exp(v))); g=grad_P(x,p); return np.sum(g[1:]**2)
    rng=np.random.default_rng(seed); found=[]; cand=0
    for _ in range(nstarts):
        v0=np.log(rng.dirichlet(np.ones(8))+1e-9)
        r=minimize(F,v0,method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-10,'fatol':1e-14})
        if r.fun<1e-10:
            x=np.array([0.0]+list(np.exp(r.x)))
            sc=x[1:].max()
            if (x[1:]/sc).min() < 1e-3:  # drifted to sub-face boundary
                continue
            xp,gn=polish(x,p)
            if xp is None: continue
            sc=xp[1:].max()
            if (xp[1:]/sc).min() < 1e-3: continue
            P=Px(xp,p)
            if not (5<P<12): continue
            key=round(P,3)
            if not any(abs(key-f[0])<0.01 for f in found):
                found.append((key,P,xp,gn)); cand+=1
    return found

def hess7(x,p):
    def P7(v):
        xx=np.array([0.0,1.0]+list(v)); return Px(xx,p)
    v0=x[2:9]; h=1e-5; H=np.zeros((7,7))
    for i in range(7):
        for j in range(i,7):
            ei=np.zeros(7);ej=np.zeros(7);ei[i]=h;ej[j]=h
            if i==j: H[i,j]=(P7(v0+ei)-2*P7(v0)+P7(v0-ei))/(h*h)
            else:
                v=(P7(v0+ei+ej)-P7(v0+ei-ej)-P7(v0-ei+ej)+P7(v0-ei-ej))/(4*h*h)
                H[i,j]=H[j,i]=v
    return H

print("="*72); print("(B.16) search: all {0}-face stationary points, P<9 & H_red PSD?"); print("="*72)
print("  (palindromic branch GPT-certified saddle det<0; looking for NON-palindromic P<9 local minima)\n", flush=True)

ps=[0.20,0.30,0.35,0.40,0.43]
for p in ps:
    found=find_stationary(p,nstarts=400,seed=hash(str(p))%2**31)
    print(f"--- p={p}: {len(found)} distinct interior stationary points ---", flush=True)
    for key,P,x,gn in sorted(found):
        H=hess7(x,p); ev=np.sort(np.linalg.eigvalsh(H))
        q=1-p; T=np.array([x[i]/(p*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9)])
        Tn=T[1:]/T[1]
        pal = np.max(np.abs(Tn - Tn[::-1])) < 0.02
        psd = ev.min() > -1e-2
        morse = int((ev<-1e-2).sum())
        flag=''
        if P<9 and psd: flag=' *** P<9 & PSD (B.16 VIOLATION!) ***'
        elif P<9 and not psd: flag=' (P<9 saddle)'
        elif P>=9: flag=' (P>=9, irrelevant to m9<9)'
        print(f"   P={P:.5f} ||g||^2={gn:.1e} pal={pal} morse={morse} eig_min={ev[0]:.3f}{flag}", flush=True)
print("\nDone. If no 'B.16 VIOLATION' lines, numerical evidence supports (B.16) emptiness.", flush=True)
