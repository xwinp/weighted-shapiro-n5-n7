#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lean Theta verification (no expansion): symmetry, {G=S=Theta=0} emptiness."""
import time, random
from pathlib import Path
import sympy as sp
import mpmath as mp
mp.mp.dps=30

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
X, Y, sigma = sp.symbols('X Y sigma')
C, D, sig = sp.symbols('C D sigma')

def load_small(name):
    s = sp.symbols("s")
    text = (HERE/name).read_text(encoding="utf-8").strip()
    expr = sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma)
    return sp.Poly(expr, X, Y, sigma, domain=sp.ZZ)

G = load_small("nonpal_G_clean.txt")
S = load_small("nonpal_S_clean.txt")
Gcd = sp.expand(G.as_expr().subs({X: C+D, Y: C*D}))
Scd = sp.expand(S.as_expr().subs({X: C+D, Y: C*D}))

def center_lift(C, D, sigma):
    gap = C + D - 1
    a3 = 1 - C + sigma*C**2*(1-C)/gap
    a6 = 1 - D + sigma*D**2*(1-D)/gap
    a2 = 1 - a3 + sigma*a3**2*(1-a3)/(a3 + C - 1)
    a7 = 1 - a6 + sigma*a6**2*(1-a6)/(a6 + D - 1)
    return a2, a3, a6, a7
def full_center_lift(C, D, sigma):
    a2,a3,a6,a7 = center_lift(C, D, sigma)
    return (sp.Integer(1),a2,a3,C,D,a6,a7,sp.Integer(1))
def compact_theta_from_terms(a, sigma):
    (a1,a2,a3,a4,a5,a6,a7,a8)=a
    h2=a2+a3-1; h3=a3+a4-1; h4=a4+a5-1; h5=a5+a6-1; h6=a6+a7-1
    W23=((a6*a7/h6)*(sigma+h2*h6/(sigma*a2*a3*a6*a7)
        -h2*h3*h4*h5*h6/(sigma**4*a2*a3**2*a4**2*a5**2*a6**2*a7))
        *(1+a2*a3*h6/(h2*a6*a7)))
    W24=(h2*h3/(sigma**2*a2*a3**2*a4)*(h4/(sigma*a4*a5)-1)
        *(1+a2*a3**2*a4*h5*h6/(h2*h3*a5*a6**2*a7)))
    W34=((a5*a6/h5)*(sigma+h3*h5/(sigma*a3*a4*a5*a6)
        -h3*h4*h5/(sigma**2*a3*a4**2*a5**2*a6))
        *(1+a3*a4*h5/(h3*a5*a6)))
    return W23,W24,W34, W23*W24+W23*W34+W24*W34

print("building Theta...", flush=True)
t0=time.time()
W23,W24,W34,Theta = compact_theta_from_terms(full_center_lift(C,D,sig), sig)
print(f"  built {time.time()-t0:.1f}s", flush=True)

# lambdify for fast numeric eval
fG = sp.lambdify((C,D,sig), Gcd, 'mpmath')
fS = sp.lambdify((C,D,sig), Scd, 'mpmath')
fT = sp.lambdify((C,D,sig), Theta, 'mpmath')

# 1. symmetry of Theta (numeric): Theta(C,D,s) vs Theta(D,C,s)
random.seed(1); sym_ok=True; nchk=0
for _ in range(8):
    c=random.uniform(0.55,0.95); d=random.uniform(0.15,0.45); s=random.uniform(0.2,1.5)
    if c+d<=1.05 or c<=d: continue
    nchk+=1
    v1=float(fT(c,d,s)); v2=float(fT(d,c,s))
    if abs(v1-v2)>1e-6*max(1.0,abs(v1)):
        sym_ok=False; print(f"  ASYM ({c:.3f},{d:.3f},{s:.3f}): {v1:.6g} vs {v2:.6g}")
print(f"Theta symmetric in C,D? {sym_ok} (checked {nchk} pts)", flush=True)

# 2. relation to Delta=0: on the G=S=0 curve, sample and check Theta sign / non-vanishing
#    Find points on G=S=0 curve (vary sigma, solve G=S=0 for C,D), evaluate Theta.
print("sampling G=S=0 curve, checking Theta...", flush=True)
random.seed(2)
curve_pts=[]
tries=0
while len(curve_pts)<20 and tries<400:
    tries+=1
    c0=random.uniform(0.5,0.98); d0=random.uniform(0.05,0.45); s0=random.uniform(0.1,1.5)
    if c0+d0<=1.02 or c0<=d0: continue
    try:
        sol=mp.findroot(lambda c,d: (fG(c,d,s0), fS(c,d,s0)), (c0,d0), tol=1e-20, maxsteps=50)
        c1,d1=float(sol[0]),float(sol[1])
        if 0<d1<c1<1 and c1+d1>1.0:
            # verify lift positivity
            a2,a3,a6,a7=[float(x.subs({C:c1,D:d1,sig:s0})) for x in center_lift(C,D,sig)]
            aa=[1,a2,a3,c1,d1,a6,a7,1]
            if all(0<ai<1 for ai in aa) and all(aa[i]+aa[i+1]>1 for i in range(7)):
                tv=float(fT(c1,d1,s0))
                curve_pts.append((c1,d1,s0,tv))
    except Exception:
        pass
print(f"  got {len(curve_pts)} strict-interior G=S=0 curve points", flush=True)
if curve_pts:
    tvs=[p[3] for p in curve_pts]
    print(f"  Theta on curve: min={min(tvs):.6g} max={max(tvs):.6g} (all nonzero? {all(abs(t)>1e-8 for t in tvs)})", flush=True)
    nz=sum(1 for t in tvs if abs(t)>1e-6)
    print(f"  |Theta|>1e-6 on {nz}/{len(curve_pts)} pts", flush=True)

# 3. multi-start solve {G=S=Theta=0}
print("multi-start solve {G=S=Theta=0}...", flush=True)
random.seed(3)
roots=[]; tries=0
while tries<200 and len(roots)<5:
    tries+=1
    c0=random.uniform(0.5,0.98); d0=random.uniform(0.05,0.45); s0=random.uniform(0.1,1.5)
    if c0+d0<=1.02 or c0<=d0: continue
    try:
        sol=mp.findroot(lambda c,d,s: (fG(c,d,s),fS(c,d,s),fT(c,d,s)),
                        (c0,d0,s0), tol=1e-25, maxsteps=60)
        c1,d1,s1=float(sol[0]),float(sol[1]),float(sol[2])
        if 0<d1<c1<1 and c1+d1>1.0 and s1>0:
            a2,a3,a6,a7=[float(x.subs({C:c1,D:d1,sig:s1})) for x in center_lift(C,D,sig)]
            aa=[1,a2,a3,c1,d1,a6,a7,1]
            if all(0<ai<1 for ai in aa) and all(aa[i]+aa[i+1]>1 for i in range(7)):
                res=max(abs(float(fG(c1,d1,s1))),abs(float(fS(c1,d1,s1))),abs(float(fT(c1,d1,s1))))
                roots.append((c1,d1,s1,res))
    except Exception:
        pass
print(f"  {tries} starts -> {len(roots)} strict-interior roots of {{G=S=Theta=0}}", flush=True)
for r in roots[:10]: print(f"    C={r[0]:.6f} D={r[1]:.6f} sig={r[2]:.6f} res={r[3]:.2e}", flush=True)
print("DONE", flush=True)
