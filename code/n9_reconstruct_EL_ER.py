#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reconstruct E_L, E_R (pre-pp boundary equations) from the center_lift recurrence,
verify G^- = pp((E_L-E_R)/(C-D)) and G^+ = pp(E_L+E_R) match the loaded polys,
and report their structure for the saturated-ideal identity.

Recurrence: h_{i-1} h_i = sigma a_i^2 (1-a_i), h_i = a_i+a_{i+1}-1, i=2..7, a1=a8=1.
  i=2: h1 h2 = sigma a2^2(1-a2), h1=a1+a2-1=a2  =>  h2 = sigma a2(1-a2)
       => E_L := (a2+a3-1) - sigma*a2*(1-a2) = 0
  i=7: h6 h7 = sigma a7^2(1-a7), h7=a7+a8-1=a7  =>  h6 = sigma a7(1-a7)
       => E_R := (a6+a7-1) - sigma*a7*(1-a7) = 0
"""
import time
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
X, Y, sigma = sp.symbols('X Y sigma')
C, D, sig = sp.symbols('C D sigma')

def load_small(name):
    s = sp.symbols("s")
    text = (HERE/name).read_text(encoding="utf-8").strip()
    return sp.Poly(sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma),
                   X, Y, sigma, domain=sp.ZZ)
Gload = load_small("nonpal_G_clean.txt")   # G^- in (X,Y,sigma)
Sload = load_small("nonpal_S_clean.txt")   # G^+ in (X,Y,sigma)

def center_lift(c, d, s):
    gap = c + d - 1
    a3 = 1 - c + s*c**2*(1-c)/gap
    a6 = 1 - d + s*d**2*(1-d)/gap
    a2 = 1 - a3 + s*a3**2*(1-a3)/(a3 + c - 1)
    a7 = 1 - a6 + s*a6**2*(1-a6)/(a6 + d - 1)
    return a2, a3, a6, a7

a2, a3, a6, a7 = center_lift(C, D, sig)
E_L = (a2 + a3 - 1) - sig*a2*(1-a2)
E_R = (a6 + a7 - 1) - sig*a7*(1-a7)

t0 = time.time()
print("together(E_L - E_R)...", flush=True)
ELmER = sp.together(E_L - E_R)
num_m, den_m = sp.fraction(ELmER)
num_m = sp.expand(num_m); den_m = sp.expand(den_m)
print(f"  done {time.time()-t0:.1f}s. num terms={len(sp.Poly(num_m,C,D,sig).terms()) if num_m!=0 else 0}", flush=True)

# check (C-D) divides num_m
q, r = sp.div(sp.Poly(num_m, C, D, sig), sp.Poly(C-D, C, D, sig), domain='ZZ')
print(f"  (E_L-E_R) num divisible by (C-D)? remainder zero = {r.is_zero}", flush=True)

print("together(E_L + E_R)...", flush=True)
ELp = sp.together(E_L + E_R)
num_p, den_p = sp.fraction(ELp)
num_p = sp.expand(num_p); den_p = sp.expand(den_p)
print(f"  done {time.time()-t0:.1f}s.", flush=True)

# Primitive parts: G^- = pp(num_m/(C-D) / den_m) etc. But den_m, den_p are the common lift denominators.
# The loaded G^- , G^+ are INTEGER polys in (X,Y,sigma). To compare, express everything in (X,Y,sigma)
# via C+D=X, C*D=Y, and clear denominators, take primitive integer part.
print("converting to (X,Y,sigma) and comparing to loaded G^-, G^+ ...", flush=True)
# Build G^- candidate: (num_m/(C-D)) with denominator den_m; the poly part is num_m/(C-D).
# Loaded G^- is an integer poly in X,Y,sigma. The relation: G^- = pp( (num_m/(C-D)) expressed in X,Y,sigma
#   with the denominator factored out and content removed, after saturating X-1).
# Simpler check: evaluate both at random (C,D,sigma) with C!=D, compare up to a (C,D,sigma)-dependent
# rational factor? No—loaded G^- is the primitive INTEGER poly. Let's instead verify the IDEAL relation
# numerically: at strict-interior points, {E_L=E_R=0} iff {G^-=G^+=0} (already known).
# Here we verify the algebraic CONSTRUCTION: substitute X=C+D,Y=C*D into loaded G^-, multiply by
# den_m*(C-D), and check it equals a constant times num_m. I.e. loaded_Gm(C,D,sig)*(C-D)*den_m == const * num_m.
Gm_CD = sp.expand(Gload.as_expr().subs({X: C+D, Y: C*D}))
Gp_CD = sp.expand(Sload.as_expr().subs({X: C+D, Y: C*D}))
# Check: Gm_CD * (C-D) * den_m  vs  num_m  (up to constant factor)
lhs = sp.expand(Gm_CD * (C-D) * den_m)
# find constant ratio by comparing one monomial coeff
P_lhs = sp.Poly(lhs, C, D, sig); P_num = sp.Poly(num_m, C, D, sig)
# ratio via a leading term
t_lhs = next(iter(P_lhs.terms())); t_num = next(iter(P_num.terms()))
if t_lhs[0] != t_num[0]:
    # find a common term
    common = set(P_lhs.terms()) & set(P_num.terms())
    if common:
        tt = next(iter(common))
        ratio = P_lhs.coeff_monomial(tt[0]) / P_num.coeff_monomial(tt[0])
    else:
        ratio = None
else:
    ratio = P_lhs.coeff_monomial(t_lhs[0]) / P_num.coeff_monomial(t_lhs[0])
print(f"  ratio Gm*(C-D)*den_m / num_m = {ratio}", flush=True)
if ratio is not None:
    diff = sp.expand(lhs - ratio*num_m)
    print(f"  G^- construction verified (lhs == ratio*num_m)? {diff==0}", flush=True)

# same for G^+
lhs2 = sp.expand(Gp_CD * den_p)
P_lhs2 = sp.Poly(lhs2, C, D, sig); P_num2 = sp.Poly(num_p, C, D, sig)
common2 = set(P_lhs2.terms()) & set(P_num2.terms())
ratio2 = None
if common2:
    tt = next(iter(common2)); ratio2 = P_lhs2.coeff_monomial(tt[0]) / P_num2.coeff_monomial(tt[0])
    diff2 = sp.expand(lhs2 - ratio2*num_p)
    print(f"  ratio G^+*den_p / num_p = {ratio2}; verified? {diff2==0}", flush=True)
else:
    print("  G^+: no common term to compare", flush=True)

print(f"\nE_L, E_R reconstruction {'VERIFIED' if (ratio is not None and (sp.expand(lhs-ratio*num_m)==0)) else 'MISMATCH'}  t={time.time()-t0:.1f}s", flush=True)
# Save num E_L, num E_R (as polynomials in C,D,sigma) for the saturation identity
def save_poly(expr, name):
    P = sp.Poly(sp.expand(expr), C, D, sig)
    out = HERE/name
    with open(out,'w') as f:
        f.write(str(P.as_expr())+"\n")
    print(f"  saved {name}: deg={P.total_degree()} terms={len(P.terms())}", flush=True)
# num E_L, num E_R individually
nL, dL = sp.fraction(sp.together(E_L)); nR, dR = sp.fraction(sp.together(E_R))
save_poly(sp.expand(nL), "nonpal_EL_num.txt")
save_poly(sp.expand(nR), "nonpal_ER_num.txt")
print("DONE", flush=True)
