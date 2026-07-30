#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Does the non-palindromic stationary CURVE {G=S=0} strict-interior exist?

For fixed sigma, G(C,D,sigma)=S(C,D,sigma)=0 is 2 eqs in (C,D) -> isolated pts.
Scan sigma over (0, sigma_max) with multi-start; collect strict-interior solutions
(where also a2,a6,a7,h_i in range). If the curve is EMPTY, B.16 closes trivially
(no non-palindromic stationary points at all); if non-empty, we must show Theta!=0 on it.
"""
import time, random
from pathlib import Path
import sympy as sp
import mpmath as mp
mp.mp.dps = 12

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
X, Y, sigma = sp.symbols('X Y sigma')
C, D, sig = sp.symbols('C D sigma')

def load_small(name):
    s = sp.symbols("s")
    text = (HERE/name).read_text(encoding="utf-8").strip()
    return sp.Poly(sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma),
                   X, Y, sigma, domain=sp.ZZ)
G = load_small("nonpal_G_clean.txt"); S = load_small("nonpal_S_clean.txt")
Gcd = sp.expand(G.as_expr().subs({X: C+D, Y: C*D}))
Scd = sp.expand(S.as_expr().subs({X: C+D, Y: C*D}))
fG = sp.lambdify((C,D,sig), Gcd, 'mpmath')
fS = sp.lambdify((C,D,sig), Scd, 'mpmath')

def center_lift(c, d, s):
    gap = c + d - 1
    a3 = 1 - c + s*c**2*(1-c)/gap
    a6 = 1 - d + s*d**2*(1-d)/gap
    a2 = 1 - a3 + s*a3**2*(1-a3)/(a3 + c - 1)
    a7 = 1 - a6 + s*a6**2*(1-a6)/(a6 + d - 1)
    return [1,a2,a3,c,d,a6,a7,1]

def strict_ok(a):
    if not all(0 < float(ai) < 1 for ai in a): return False
    if not all(float(a[i]+a[i+1]) > 1 for i in range(7)): return False
    return True

random.seed(11)
sols = []   # (c,d,s, res)
t0 = time.time()
sigmas = [0.03,0.08,0.15,0.25,0.4,0.6,0.9,1.3,1.8,2.5,3.5]
tries = 0
for s0 in sigmas:
    for _ in range(25):
        tries += 1
        c0 = random.uniform(0.5, 0.98); d0 = random.uniform(0.02, 0.48)
        if c0+d0 <= 1.02 or c0 <= d0: continue
        try:
            sol = mp.findroot(lambda c,d:(fG(c,d,s0), fS(c,d,s0)),
                              (mp.mpf(c0), mp.mpf(d0)), tol=1e-25, maxsteps=50)
            c1,d1 = float(sol[0]), float(sol[1])
            if not (0 < d1 < c1 < 1 and c1+d1 > 1.0): continue
            res = max(abs(float(fG(c1,d1,s0))), abs(float(fS(c1,d1,s0))))
            if res > 1e-8: continue
            a = center_lift(c1,d1,s0)
            if not strict_ok(a): continue
            sols.append((c1,d1,s0,res))
        except Exception:
            pass
print(f"tries={tries}  strict-interior {{G=S=0}} solutions: {len(sols)}", flush=True)
# dedup
uniq = []
for c,d,s,r in sorted(sols, key=lambda x:(x[2],x[0],x[1])):
    if any(abs(c-uc)<1e-4 and abs(d-ud)<1e-4 and abs(s-us)<1e-4 for uc,ud,us,_ in uniq): continue
    uniq.append((c,d,s,r))
for c,d,s,r in uniq[:40]:
    print(f"  C={c:.6f} D={d:.6f} sig={s:.4f} res={r:.1e}", flush=True)
print(f"unique={len(uniq)}  t={time.time()-t0:.0f}s DONE", flush=True)
