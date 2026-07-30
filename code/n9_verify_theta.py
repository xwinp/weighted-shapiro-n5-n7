#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify Theta locally:
  1. Build Theta from GPT's code (center_lift + compact_theta).
  2. Check C<->D symmetry of Theta_num (GPT claims numerator symmetric in C,D).
  3. Numerical multi-start: {G=S=Theta=0} in strict interior -> expect empty.
  4. Check Theta=0 and Delta=0 coincide on the G=S=0 curve (sample curve points).
  5. Report Theta_num degree/size for assessing the certificate path.
"""
import re, time
from pathlib import Path
import sympy as sp
from sympy import Rational as R

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
Xs, Ys, ss = sp.symbols('X Y s')
X, Y, sigma = sp.symbols('X Y sigma')
C, D, sig = sp.symbols('C D sigma')

def load_small(name):
    s = sp.symbols("s")
    text = (HERE/name).read_text(encoding="utf-8").strip()
    expr = sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma)
    return sp.Poly(expr, X, Y, sigma, domain=sp.ZZ)

G = load_small("nonpal_G_clean.txt")
S = load_small("nonpal_S_clean.txt")
# substitute X->C+D, Y->C*D to get polys in (C,D,sigma)
Gcd = sp.Poly(sp.expand(G.as_expr().subs({X: C+D, Y: C*D})), C, D, sigma, domain=sp.ZZ)
Scd = sp.Poly(sp.expand(S.as_expr().subs({X: C+D, Y: C*D})), C, D, sigma, domain=sp.ZZ)
print("Gcd deg/terms:", Gcd.total_degree(), len(Gcd.terms()), flush=True)
print("Scd deg/terms:", Scd.total_degree(), len(Scd.terms()), flush=True)

# ---- Theta construction (verbatim from GPT reply #4) ----
def center_lift(C, D, sigma):
    gap = C + D - 1
    a3 = 1 - C + sigma*C**2*(1-C)/gap
    a6 = 1 - D + sigma*D**2*(1-D)/gap
    a2 = 1 - a3 + sigma*a3**2*(1-a3)/(a3 + C - 1)
    a7 = 1 - a6 + sigma*a6**2*(1-a6)/(a6 + D - 1)
    return a2, a3, a6, a7

def full_center_lift(C, D, sigma):
    a2, a3, a6, a7 = center_lift(C, D, sigma)
    return (sp.Integer(1), a2, a3, C, D, a6, a7, sp.Integer(1))

def compact_theta_from_terms(a, sigma):
    (a1,a2,a3,a4,a5,a6,a7,a8) = a
    h2=a2+a3-1; h3=a3+a4-1; h4=a4+a5-1; h5=a5+a6-1; h6=a6+a7-1
    W23 = ((a6*a7/h6) * (sigma + h2*h6/(sigma*a2*a3*a6*a7)
          - h2*h3*h4*h5*h6/(sigma**4*a2*a3**2*a4**2*a5**2*a6**2*a7))
          * (1 + a2*a3*h6/(h2*a6*a7)))
    W24 = (h2*h3/(sigma**2*a2*a3**2*a4) * (h4/(sigma*a4*a5) - 1)
          * (1 + a2*a3**2*a4*h5*h6/(h2*h3*a5*a6**2*a7)))
    W34 = ((a5*a6/h5) * (sigma + h3*h5/(sigma*a3*a4*a5*a6)
          - h3*h4*h5/(sigma**2*a3*a4**2*a5**2*a6))
          * (1 + a3*a4*h5/(h3*a5*a6)))
    Theta = W23*W24 + W23*W34 + W24*W34
    return W23, W24, W34, Theta

print("building Theta...", flush=True)
t0=time.time()
W23,W24,W34,Theta = compact_theta_from_terms(full_center_lift(C,D,sigma), sigma)
print(f"  built in {time.time()-t0:.1f}s", flush=True)

# numerical Theta eval
def num_theta(c,d,s):
    return float(Theta.subs({C:c,D:d,sigma:s}))

# sanity: pick a strict-interior point, check Theta finite & well-defined
c0,d0,s0 = R(6,10), R(4,10), R(1,2)  # C=0.6,D=0.4 (C+D=1.0... need >1)
c0,d0,s0 = R(7,10), R(4,10), R(1,2)  # C=0.7,D=0.4, sum=1.1>1, C>D
try:
    tv = num_theta(c0,d0,s0)
    print(f"Theta(0.7,0.4,0.5) = {tv}", flush=True)
except Exception as e:
    print("Theta eval error:", e, flush=True)

# symmetry: Theta_num(C,D,s) vs Theta_num(D,C,s)
print("checking symmetry of Theta_num (this may take a while)...", flush=True)
t0=time.time()
Theta_t = sp.together(Theta)
Theta_num, Theta_den = Theta_t.as_numer_denom()
Theta_num = sp.expand(Theta_num)
Theta_den = sp.expand(Theta_den)
print(f"  together done in {time.time()-t0:.1f}s", flush=True)
# degree/size
Pnum = sp.Poly(Theta_num, C, D, sigma)
Pden = sp.Poly(Theta_den, C, D, sigma)
print(f"Theta_num: tot_deg={Pnum.total_degree()} terms={len(Pnum.terms())} degC={Pnum.degree(C)} degD={Pnum.degree(D)} degS={Pnum.degree(sigma)}", flush=True)
print(f"Theta_den: tot_deg={Pden.total_degree()} terms={len(Pden.terms())}", flush=True)
# symmetry check via numerical substitution at random points
import random
random.seed(1)
sym_ok=True
for _ in range(5):
    c=random.uniform(0.55,0.95); d=random.uniform(0.1,0.45); s=random.uniform(0.1,1.5)
    if c+d<=1.05: continue
    v1=float(Theta_num.subs({C:c,D:d,sigma:s}))
    v2=float(Theta_num.subs({C:d,D:c,sigma:s}))
    if abs(v1-v2)>1e-6*max(1,abs(v1)):
        sym_ok=False; print(f"  ASYM at ({c},{d},{s}): {v1} vs {v2}", flush=True)
print(f"Theta_num symmetric in C,D? {sym_ok}", flush=True)
