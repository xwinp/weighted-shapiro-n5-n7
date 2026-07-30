#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Task B numeric ground-truth: minimize P on {0} face (L=9) directly.
Confirm: (a) interior stationary is a saddle w/ det H_red<0 when P<9;
         (b) P=9 crossing p_P < b9=0.4338858820; (c) at p=0.4, P~8.829694."""
import numpy as np
from scipy.optimize import minimize
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

def P_full(x, p):
    q = 1-p
    s = 0.0
    for i in range(9):
        d = p*x[(i+1)%9] + q*x[(i+2)%9]
        if d <= 1e-300: return 1e18
        s += x[i]/d
    return s

def face0_min(p, nstarts=60, seed=1):
    """Minimize P on {0} face: x0=0, x1..x8>0. Parametrize log x1..x8, fix gauge sum=1."""
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(nstarts):
        x0 = rng.dirichlet(np.ones(8))
        def F(v):
            xv = np.exp(v - v.mean())       # positive, gauge-invariant-ish
            xx = np.concatenate(([0.0], xv))
            return P_full(xx, p)
        v0 = np.log(x0 + 1e-9)
        res = minimize(F, v0, method='Nelder-Mead',
                       options={'maxiter':20000,'xatol':1e-12,'fatol':1e-12})
        if best is None or res.fun < best[0]:
            xv = np.exp(res.x - res.x.mean())
            xx = np.concatenate(([0.0], xv))
            best = (res.fun, xx)
    return best

def hessian_red(x, p):
    """Reduced Hessian of P on {0} face wrt x2..x8 (x0=0,x1=1 gauge)."""
    def P7(v):
        xx = np.array([0.0,1.0]+list(v))
        return P_full(xx, p)
    v0 = x[2:9]
    h = 1e-4
    H = np.zeros((7,7))
    for i in range(7):
        for j in range(7):
            ei=np.zeros(7); ej=np.zeros(7); ei[i]=h; ej[j]=h
            if i==j:
                H[i,j]=(P7(v0+ei)-2*P7(v0)+P7(v0-ei))/(h*h)
            else:
                H[i,j]=(P7(v0+ei+ej)-P7(v0+ei-ej)-P7(v0-ei+ej)+P7(v0-ei-ej))/(4*h*h)
    return H

b9 = 0.43388588203369836
print("="*70); print("TASK B numeric: L=9 ({0}-face) interior stationary"); print("="*70)

# (c) p=0.4
p=0.4
Pmin, xmin = face0_min(p, nstarts=80, seed=3)
print(f"\n[c] p=0.4:  P_min = {Pmin:.9f}  (GPT 8.829694260862)  <9? {Pmin<9}")
T = np.array([xmin[i]/(p*xmin[(i+1)%9]+(1-p)*xmin[(i+2)%9]) for i in range(9)])
print(f"   term values T_i (i=0..8): {np.round(T,5)}")
print(f"   T0={T[0]:.2e} (should be ~0);  T1..T8 = {np.round(T[1:],5)}")
# normalize by T[8] to compare with (1,A,B,C,C,B,A,1)
seq = T[1:]/T[8]
print(f"   T1..T8 / T8 = {np.round(seq,5)}  (compare 1,A,B,C,C,B,A,1 → first=1, last=1)")
H = hessian_red(xmin, p)
ev = np.linalg.eigvalsh(H)
print(f"   reduced Hessian eig: {np.round(np.sort(ev),4)}")
print(f"   det(H_red)={np.linalg.det(H):.4e}  Morse_idx={int((ev<-1e-3).sum())}  (GPT: saddle det<0)")

# (a)(b) scan p, find interior min crossing P=9
print("\n[a][b] scan: P_min({0}-face) vs p, find P=9 crossing p_P")
ps = np.linspace(0.05, 0.60, 24)
cross=[]
prev=None
for p in ps:
    Pm,_ = face0_min(p, nstarts=40, seed=5)
    flag = '<9 LOW' if Pm<9-1e-4 else '>=9'
    print(f"   p={p:.4f}  P_min={Pm:.6f}  {flag}")
    if prev is not None and (prev-9)*(Pm-9)<0:
        cross.append((p-0.0125, p))  # bracket
    prev=Pm
print(f"\n   crossings near: {cross}")
print(f"   GPT p_P=0.4318363763 < b9={b9:.10f}? {0.4318363763<b9}")
