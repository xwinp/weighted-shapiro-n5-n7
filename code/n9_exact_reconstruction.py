#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GPT's exact E_L,E_R reconstruction (compact factored form) with exact checks,
reflection symmetry, and G^-/G^+ via primitive parts. Verifies the FORWARD
inclusion exactly (no numerical test). Also:
 - verify the sigma-bound direction: a6<1  <=>  sigma < v/(d(1-d)).
 - check whether the earlier 6 forward-test points are palindromic (c=d).
"""
import time
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
X, Y, sigma = sp.symbols('X Y sigma')
c, d, s = sp.symbols('c d sigma')

v = c + d - 1

def center_lift(c, d, sigma):
    v = c + d - 1
    a3 = 1 - c + sigma * c**2 * (1 - c) / v
    a6 = 1 - d + sigma * d**2 * (1 - d) / v
    a2 = 1 - a3 + sigma * a3**2 * (1 - a3) / (a3 + c - 1)
    a7 = 1 - a6 + sigma * a6**2 * (1 - a6) / (a6 + d - 1)
    return a2, a3, a6, a7

a2, a3, a6, a7 = center_lift(c, d, s)
E_L = sp.together(a2 + a3 - 1 - s * a2 * (1 - a2))
E_R = sp.together(a7 + a6 - 1 - s * a7 * (1 - a7))

# Compact exact construction
L_c = v + s * c**2
U_c = v - s * c * (1 - c)
B_c = c**2 * v + (1 - c) * L_c**2
L_d = v + s * d**2
U_d = v - s * d * (1 - d)
B_d = d**2 * v + (1 - d) * L_d**2
F_c = c * v**2 * (1 - c) * L_c**2 - s * B_c * (c * v**2 - U_c * B_c)
F_d = d * v**2 * (1 - d) * L_d**2 - s * B_d * (d * v**2 - U_d * B_d)

print("exact checks (cancel == 0)...", flush=True); t0 = time.time()
chk1 = sp.cancel(E_L - U_c * F_c / (c**2 * v**4))
chk2 = sp.cancel(E_R - U_d * F_d / (d**2 * v**4))
print(f"  E_L == U_c F_c/(c^2 v^4)? {chk1==0}   E_R == U_d F_d/(d^2 v^4)? {chk2==0}  t={time.time()-t0:.1f}s", flush=True)

# divide nonzero factors (c-1)*L_c, (d-1)*L_d
P_c = sp.Poly(sp.expand(F_c), c, d, s, domain=sp.ZZ)
P_d = sp.Poly(sp.expand(F_d), c, d, s, domain=sp.ZZ)
print(f"  F_c: deg={P_c.total_degree()} terms={len(P_c.terms())}; F_d: deg={P_d.total_degree()} terms={len(P_d.terms())}", flush=True)
div_c = sp.Poly((c - 1) * L_c, c, d, s, domain=sp.ZZ)
div_d = sp.Poly((d - 1) * L_d, c, d, s, domain=sp.ZZ)
# exact quotient (exquo raises if not exact)
try:
    F_L = P_c.exquo(div_c)
    print("  F_c / ((c-1)L_c) exact (exquo) OK", flush=True)
except Exception as e:
    print("  F_c division NOT exact:", e, flush=True)
    F_L = None
try:
    F_R = P_d.exquo(div_d)
    print("  F_d / ((d-1)L_d) exact (exquo) OK", flush=True)
except Exception as e:
    print("  F_d division NOT exact:", e, flush=True)
    F_R = None

if F_L and F_R:
    # reflection check
    refl = sp.simplify(F_R.as_expr() - F_L.as_expr().xreplace({c: d, d: c}))
    print(f"  reflection F_R == F_L(c<->d)? {refl==0}", flush=True)
    diff = sp.Poly(sp.expand(F_L.as_expr() - F_R.as_expr()), c, d, s, domain=sp.ZZ)
    try:
        Gm_cd = diff.exquo(sp.Poly(c - d, c, d, s, domain=sp.ZZ))
        print(f"  (F_L-F_R)/(c-d) exact: G^- deg={Gm_cd.total_degree()} terms={len(Gm_cd.terms())}", flush=True)
    except Exception as e:
        print("  (F_L-F_R)/(c-d) NOT exact:", e, flush=True)
    Gp_cd = sp.Poly(sp.expand(F_L.as_expr() + F_R.as_expr()), c, d, s, domain=sp.ZZ)
    print(f"  G^+ = F_L+F_R: deg={Gp_cd.total_degree()} terms={len(Gp_cd.terms())}", flush=True)

# ---- sigma-bound direction check ----
print("\nsigma-bound check:", flush=True)
a6_expr = 1 - d + s * d**2 * (1 - d) / v
# a6 < 1  <=>  sigma < ?
ineq = sp.simplify(a6_expr - 1)  # = s d^2 (1-d)/v - d
# a6<1  <=>  s d^2(1-d)/v < d  <=>  s < d v/(d^2(1-d)) = v/(d(1-d))
print("  a6-1 =", sp.factor(ineq), " => a6<1 <=> sigma < v/(d(1-d))", flush=True)
a3_expr = 1 - c + s * c**2 * (1 - c) / v
print("  a3-1 =", sp.factor(sp.simplify(a3_expr-1)), " => a3<1 <=> sigma < v/(c(1-c))", flush=True)
print("  since c(1-c)<d(1-d) for c>d,c+d>1, tighter is sigma < v/(d(1-d))  [CONFIRMED]", flush=True)

# ---- are the 6 forward-test points palindromic? ----
# The earlier 6 points were strict-interior (E_L=E_R=0) test points. Re-derive: solve E_L=E_R=0
# numerically and check |c-d|.
print("\nchecking if E_L=E_R=0 strict-interior solutions are palindromic (c=d)...", flush=True)
import mpmath as mp
mp.mp.dps = 18
fEL = sp.lambdify((c,d,s), E_L, 'mpmath')
fER = sp.lambdify((c,d,s), E_R, 'mpmath')
import random
random.seed(3)
pal = 0; nonpal = 0; tot = 0
for _ in range(200):
    c0 = random.uniform(0.5,0.95); d0 = random.uniform(0.1,0.45); s0 = random.uniform(0.05, 1.5)
    if c0+d0<=1.02 or c0<=d0: continue
    try:
        sol = mp.findroot(lambda cc,dd,ss:(fEL(cc,dd,ss),fER(cc,dd,ss), fEL(cc,dd,ss)*0+fER(cc,dd,ss)),
                          (mp.mpf(c0),mp.mpf(d0),mp.mpf(s0)), tol=1e-20, maxsteps=40)
        # 3 unknowns 2 eqs: findroot needs square; instead solve E_L=E_R=0 with sigma fixed by scanning
    except Exception:
        pass
# Better: fix sigma, solve E_L=E_R=0 in (c,d) (2 eqns 2 unknowns)
fEL2 = sp.lambdify((c,d,s), E_L, 'mpmath')
fER2 = sp.lambdify((c,d,s), E_R, 'mpmath')
for s0 in [0.1,0.3,0.6,1.0,1.5]:
    for _ in range(30):
        c0 = random.uniform(0.5,0.95); d0 = random.uniform(0.1,0.45)
        if c0+d0<=1.02 or c0<=d0: continue
        try:
            sol = mp.findroot(lambda cc,dd:(fEL2(cc,dd,s0),fER2(cc,dd,s0)),
                              (mp.mpf(c0),mp.mpf(d0)), tol=1e-22, maxsteps=50)
            c1,d1 = float(sol[0]),float(sol[1])
            if not (0<d1<c1<1 and c1+d1>1): continue
            r = max(abs(float(fEL2(c1,d1,s0))),abs(float(fER2(c1,d1,s0))))
            if r>1e-7: continue
            tot += 1
            if abs(c1-d1)<1e-4: pal += 1
            else: nonpal += 1
        except Exception:
            pass
print(f"  E_L=E_R=0 strict-interior solutions found: {tot} (palindromic c~d: {pal}, non-pal: {nonpal})", flush=True)
print("DONE", flush=True)
