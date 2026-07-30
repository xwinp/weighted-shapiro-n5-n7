#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract (p, rho, A, B, C, lambda, P) along the {0}-face palindromic saddle
branch, to identify the correct P-formula (B.7) and lambda(A,B,C,rho).
Term values T_1..T_8 = (lam, lam*A, lam*B, lam*C, lam*C, lam*B, lam*A, lam)."""
import numpy as np
from scipy.optimize import minimize, root
import sys, warnings, json
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
        if gn<1e-20: return xx
    return None
def find_saddle(p, seed=1):
    def F(v):
        x=np.array([0.0]+list(np.exp(v))); g=grad_P(x,p); return np.sum(g[1:]**2)
    rng=np.random.default_rng(seed)
    best=None
    for _ in range(400):
        v0=np.log(rng.dirichlet(np.ones(8))+1e-9)
        r=minimize(F,v0,method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-11,'fatol':1e-15})
        if r.fun<1e-12:
            x=np.array([0.0]+list(np.exp(r.x))); xp=polish(x,p)
            if xp is not None and 8<Px(xp,p)<9.05:
                return xp
    return None

rows=[]
for p in [0.38,0.40,0.41,0.42,0.425,0.43,0.4318363763]:
    xs=find_saddle(p, seed=hash(str(p))%2**31)
    if xs is None:
        print(f"p={p}: no saddle found"); continue
    xs=xs/xs[1]
    q=1-p; rho=p/q
    T=np.array([xs[i]/(p*xs[(i+1)%9]+q*xs[(i+2)%9]) for i in range(9)])
    lam=T[1]; A=T[2]/T[1]; B=T[3]/T[1]; C=T[4]/T[1]
    P=Px(xs,p)
    # check palindrome
    Tn=T[1:]/T[1]; pal=np.max(np.abs(Tn-Tn[::-1]))
    # check closure (B.6)
    rho9=(2*C-1)*(A+B-1)**4*(B+C-1)**2/(A**2*(1-A)**4*(1-B)**2*(1-C)**2)
    row=dict(p=p,rho=rho,A=A,B=B,C=C,lam=lam,P=P,pal=pal,rho9=rho9,rho9_9th=rho9**(1/9))
    rows.append(row)
    print(f"p={p:.10f} rho={rho:.6f} A={A:.10f} B={B:.10f} C={C:.10f} lam={lam:.10f} P={P:.10f} pal={pal:.1e} rho9^(1/9)={rho9**(1/9):.6f} (rho={rho:.6f})", flush=True)

# Try to identify lambda. Candidate forms and residuals:
print("\n=== lambda formula identification ===", flush=True)
for r in rows:
    A,B,C,lam,rho=r['A'],r['B'],r['C'],r['lam'],r['rho']
    # GPT (B.7) lambda-equiv: rho(1+rho)A(1-A)/(A+B-1)
    g1=rho*(1+rho)*A*(1-A)/(A+B-1)
    print(f"p={r['p']:.4f}: lam={lam:.8f}  rho(1+rho)A(1-A)/(A+B-1)={g1:.8f}  ratio={lam/g1:.6f}", flush=True)
json.dump(rows, open('paper/_branch_data.json','w'), indent=2)
print("\nsaved paper/_branch_data.json", flush=True)
