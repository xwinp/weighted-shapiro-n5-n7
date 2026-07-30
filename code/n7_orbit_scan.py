#!/usr/bin/env python3
"""
n=7 global proof — per-orbit numerical scan.
For each C_7 orbit of the zero-set (support), minimize P_{7,p,1-p} over the
support (zeros fixed at 0, nonzeros free on simplex) via SLSQP, scanning p.
Report min P over p in (0,1) for each orbit. Only the main orbit should dip <7.
Orbit representatives (zero index sets, cyclic-normalized):
  O0 no-zero      : {}
  O1 1-zero       : {0}
  O2 2-zero dist2 : {0,5}   <- main (band (a7,b7))
  O3 2-zero dist3 : {0,4}
  O4 3-zero 2,2,3 : {0,2,4}
"""
import numpy as np
from scipy.optimize import minimize

n = 7
orbits = {
    "O0 nozero":   set(),
    "O1 1zero":    {0},
    "O2 2z-d2 MAIN": {0,5},
    "O3 2z-d3":    {0,4},
    "O4 3z-223":   {0,2,4},
}

def Pval(x, p):
    q = 1-p
    s = 0.0
    for i in range(n):
        den = p*x[(i+1)%n] + q*x[(i+2)%n]
        if den <= 1e-15: return 1e6
        s += x[i]/den
    return s

def min_over_orbit(zeros, p, nstarts=8):
    free = [i for i in range(n) if i not in zeros]
    k = len(free)
    best = 1e9
    rng = np.random.RandomState(12345)
    for s in range(nstarts):
        # init: uniform on free coords, tiny on zeros
        y = rng.rand(k) + 0.1
        y = y/ y.sum()
        x0 = np.zeros(n)
        for j,idx in enumerate(free): x0[idx] = y[j]
        # constraint sum(x)=1, bounds: zeros fixed 0, free >=1e-9
        bounds = [(1e-9,1.0) if i in free else (0.0,0.0) for i in range(n)]
        cons = [{'type':'eq','fun':lambda x: x.sum()-1.0}]
        try:
            r = minimize(lambda x: Pval(x,p), x0, method='SLSQP',
                         bounds=bounds, constraints=cons,
                         options={'maxiter':400,'ftol':1e-12})
            if r.fun < best: best = r.fun
        except Exception:
            pass
    return best

ps = np.linspace(0.05, 0.49, 23)
print(f"{'orbit':16s} " + " ".join(f"p={pp:.3f}" for pp in ps[::2]))
for name, z in orbits.items():
    row = [min_over_orbit(z, pp) for pp in ps]
    mn = min(row)
    tag = "  <-- <7!" if mn < 7-1e-4 else ""
    print(f"{name:16s} minP={mn:.5f} at p={ps[int(np.argmin(row))]:.3f}{tag}")
    # print compact row at every other p
    print("   " + " ".join(f"{v:.3f}" for v in row[::2]))
