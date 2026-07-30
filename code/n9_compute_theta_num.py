#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compute Theta_num(C,D,sigma) = numerator of together(Theta), report size,
and attempt symmetrization to (X,Y,sigma). One-time heavy computation."""
import time, sys
from pathlib import Path
import sympy as sp

C, D, sig = sp.symbols('C D sigma')
X, Y, sigma = sp.symbols('X Y sigma')

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
    return W23*W24+W23*W34+W24*W34

print("building Theta...", flush=True)
t0=time.time()
Theta = compact_theta_from_terms(full_center_lift(C,D,sig), sig)
print(f"  built {time.time()-t0:.1f}s", flush=True)

print("together...", flush=True); t0=time.time()
Theta_t = sp.together(Theta)
print(f"  together {time.time()-t0:.1f}s", flush=True)
Theta_num, Theta_den = Theta_t.as_numer_denom()
print("expand num...", flush=True); t0=time.time()
Theta_num = sp.expand(Theta_num)
Theta_den = sp.expand(Theta_den)
print(f"  expand {time.time()-t0:.1f}s", flush=True)

Pn = sp.Poly(Theta_num, C, D, sig)
Pd = sp.Poly(Theta_den, C, D, sig)
print(f"Theta_num: tot_deg={Pn.total_degree()} terms={len(Pn.terms())} degC={Pn.degree(C)} degD={Pn.degree(D)} degS={Pn.degree(sig)}", flush=True)
print(f"Theta_den: tot_deg={Pd.total_degree()} terms={len(Pd.terms())} degS={Pd.degree(sig)}", flush=True)
sys.stdout.flush()

# save numerator as term list (C,D,sigma exponents + coeff)
out = Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'/'theta_num_CD.txt'
with open(out,'w') as f:
    for (dc,dd,ds),c in Pn.terms():
        f.write(f"{c}*C**{dc}*D**{dd}*s**{ds}\n")
print(f"saved theta_num_CD.txt ({len(Pn.terms())} terms)", flush=True)

# symmetrize: check symmetric, then convert C,D -> X=C+D, Y=CD
print("checking symmetry & symmetrizing...", flush=True); t0=time.time()
# symmetric if Pn(C,D) == Pn(D,C); test numerically
import random
random.seed(3); sym_ok=True
for _ in range(6):
    c=random.uniform(0.5,0.95); d=random.uniform(0.1,0.45); s=random.uniform(0.05,1.0)
    if c+d<=1.05: continue
    v1=float(Pn.eval({C:c,D:d,sig:s}))
    v2=float(Pn.eval({C:d,D:c,sig:s}))
    if abs(v1-v2)>1e-6*max(1,abs(v1)): sym_ok=False; break
print(f"  Theta_num symmetric? {sym_ok} ({time.time()-t0:.1f}s)", flush=True)

if sym_ok:
    print("converting to (X,Y,sigma) via C+D=X, CD=Y ...", flush=True); t0=time.time()
    # use sympy symmetric poly rewrite: substitute D=X-C, then Y appears as C*(X-C)
    # Better: use the fact symmetric poly in C,D = poly in e1=X,e2=Y.
    # Convert each monomial C^a D^b: if a==b -> Y^a; else use power-sum reduction.
    # Use sympy's poly in e1,e2 via .rewrite with syms.
    try:
        Pxy = sp.Poly(Pn.as_expr(), C, D, sig)
        # symmetric reduction: introduce e1,e2
        e1,e2 = sp.symbols('e1 e2')
        # replace using groebner / direct: use sp.poly over C,D and reduce by e1-C-D, e2-CD
        from sympy import groebner
        # Build target: express Pn in terms of e1=C+D, e2=C*D
        expr = Pn.as_expr()
        # Use the substitution method: D = e1 - C, CD=e2 => C(e1-C)=e2 => C^2 - e1*C + e2=0
        # Reduce expr (in C) modulo (C^2 - e1*C + e2), treating e1,e2 as symbols, result independent of C
        C2red = sp.Poly(C**2 - e1*C + e2, C)
        # reduce expr as poly in C with coeffs in (D... no D=e1-C). Substitute D->e1-C first.
        expr2 = sp.expand(expr.subs(D, e1 - C))
        Pexpr = sp.Poly(expr2, C)
        rem = sp.rem(Pexpr, C2red, C)  # remainder mod C^2-e1*C+e2
        # rem should be degree<2 in C; for symmetric expr, rem is constant in C (or linear that cancels)
        # Express: if rem has C term, it's not symmetric (shouldn't happen). Take the C^0 part.
        rem_expr = rem.as_expr()
        # collect
        coeff_C = rem_expr.coeff(C,1)
        const = rem_expr.coeff(C,0)
        if coeff_C!=0:
            print(f"  WARNING: residual C term {coeff_C} (not fully symmetric?)", flush=True)
        Theta_num_sym = sp.expand(const.subs({e1:X, e2:Y}))
        Psym = sp.Poly(Theta_num_sym, X, Y, sigma)
        print(f"  Theta_num_sym(X,Y,sigma): tot_deg={Psym.total_degree()} terms={len(Psym.terms())} degX={Psym.degree(X)} degY={Psym.degree(Y)} degS={Psym.degree(sigma)} ({time.time()-t0:.1f}s)", flush=True)
        out2 = Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'/'theta_num_XY.txt'
        with open(out2,'w') as f:
            for (dx,dy,ds),c in Psym.terms():
                f.write(f"{c}*X**{dx}*Y**{dy}*s**{ds}\n")
        print(f"  saved theta_num_XY.txt ({len(Psym.terms())} terms)", flush=True)
    except Exception as e:
        print("  symmetrize failed:", repr(e), flush=True)
print("ALL DONE", flush=True)
