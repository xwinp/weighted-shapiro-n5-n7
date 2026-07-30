#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Follow-up: (1) L=5/L=7 cross at p*=0.3924290; (2) L=9 interior saddle value."""
import numpy as np
from scipy.optimize import minimize
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
n=9
def P_val(x,p):
    q=1-p; s=0.0
    for i in range(n):
        d=p*x[(i+1)%n]+q*x[(i+2)%n]
        if d<=0: return 1e18
        s+=x[i]/d
    return s
def face_min(zeros,p,nstarts=600,seed=0):
    supp=[i for i in range(n) if i not in zeros]; k=len(supp)
    rng=np.random.default_rng(seed); best=np.inf
    for _ in range(nstarts):
        y=rng.dirichlet(np.ones(k)); x=np.zeros(n); x[supp]=y
        v=P_val(x,p)
        if v<best: best=v
    def f(y):
        x=np.zeros(n); x[supp]=np.abs(y); s=x.sum()
        if s<1e-15: return 1e18
        x=x/s; return P_val(x,p)
    res=minimize(f, np.abs(np.random.default_rng(99).dirichlet(np.ones(k))), method='Nelder-Mead',
                 options={'maxiter':8000,'xatol':1e-12,'fatol':1e-13})
    return min(best, f(res.x))

# (1) cross at p*=0.3924290
pstar=0.3924290
m5=face_min({0,2,4},pstar,seed=3); m7=face_min({0,2},pstar,seed=5)
print(f"p*={pstar}: M_9,5={m5:.8f}  M_9,7={m7:.8f}  (GPT cross value 8.78637382)")
# scan finer to locate true cross
ps=np.linspace(0.30,0.42,25)
prev=None
for p in ps:
    a=face_min({0,2,4},p,seed=3); b=face_min({0,2},p,seed=5)
    diff=a-b
    print(f"  p={p:.5f}: M95={a:.6f} M97={b:.6f} diff={diff:+.6f} {'<-- cross' if prev is not None and prev*diff<0 else ''}")
    prev=diff

# (2) L=9 interior saddle: minimize on {0} face from STRICTLY interior starts, keep min-component large
print("\nL=9 one-zero face {0}: find interior stationary (saddle) vs boundary global min")
for p in [0.4, 0.43]:
    supp=[i for i in range(n) if i!=0]; k=8
    # boundary global min (= L=7):
    mb=face_min({0,2},p,seed=5)
    # search interior stationary: many Nelder-Mead from interior Dirichlet, report those with min-comp>0.02
    cands=[]
    rng=np.random.default_rng(7)
    for _ in range(2000):
        y=rng.dirichlet(np.ones(k)*2.0)  # more uniform -> stays interior
        x=np.zeros(n); x[supp]=y
        def f(z):
            xx=np.zeros(n); xx[supp]=np.abs(z); s=xx.sum()
            if s<1e-15: return 1e18
            xx=xx/s; return P_val(xx,p)
        res=minimize(f, y, method='Nelder-Mead', options={'maxiter':6000,'xatol':1e-12,'fatol':1e-13})
        z=np.abs(res.x); z=z/z.sum(); x=np.zeros(n); x[supp]=z
        if x.min()>0.02:  # strictly interior
            cands.append((P_val(x,p), x.min(), res.fun))
    cands.sort()
    print(f" p={p}: boundary(L7) global min={mb:.6f}")
    if cands:
        print(f"   best interior stationary (min-comp>0.02): P={cands[0][0]:.6f}  mincomp={cands[0][1]:.4f}  (GPT saddle 8.829694@0.4, 8.990733@0.43)")
    else:
        print("   no strictly-interior stationary found (all drift to boundary)")
