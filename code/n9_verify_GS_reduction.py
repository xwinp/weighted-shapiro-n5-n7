#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify GPT's G/S reduction against the DEFINITION (4.1).

Definition (4.1): for normalized term values a=(1,a2,a3,C,D,a6,a7,1) with a1=a8=1,
  (a_{i-1}+a_i-1)(a_i+a_{i+1}-1) = sigma * a_i^2 * (1-a_i),  i=2..7
center_lift recovers a3,a6,a2,a7 from C=a4,D=a5,sigma using i=3,4,5,6 eqns.
The remaining boundary eqns are:
  E_L (i=2): a2+a3-1 - sigma*a2^2*(1-a2) = 0   [since (a1+a2-1)=a2, (a2+a3-1) is LHS factor]
    wait: i=2: (a1+a2-1)(a2+a3-1) = sigma a2^2(1-a2); a1=1 so (a1+a2-1)=a2.
    => a2*(a2+a3-1) = sigma*a2^2*(1-a2) => (a2+a3-1) = sigma*a2*(1-a2)  [if a2!=0]
    => E_L := a2+a3-1 - sigma*a2*(1-a2) = 0
  i=7: (a6+a7-1)(a7+a8-1)=sigma a7^2(1-a7); a8=1 so (a7+a8-1)=a7
    => (a6+a7-1)*a7 = sigma*a7^2*(1-a7) => E_R := a6+a7-1 - sigma*a7*(1-a7) = 0

GPT: G = primitive( (E_L - E_R)/(C-D) ) with (X-1) factor saturated.
      S = primitive of the other combination (E_L+E_R or E_L alone after saturation).
We verify (E_L-E_R)/(C-D) as a symmetric poly in (X,Y,sigma) equals G * (X-1)^k."""
import re, sympy as sp
from pathlib import Path

HERE=Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'
X,Y,s=sp.symbols('X Y s')
C,D,sigma=sp.symbols('C D sigma')

def parse_poly(name, symX, symY, syms):
    txt=(HERE/name).read_text().strip()
    terms=re.findall(r'[+-][^+-]+|^[^+-]+', txt)
    monos={}
    for t in terms:
        t=t.strip()
        if not t: continue
        sign=1
        if t.startswith('-'): sign=-1; t=t[1:].strip()
        elif t.startswith('+'): t=t[1:].strip()
        xp=yp=sp_=0
        for fac in re.finditer(r'(X|Y|s)(?:\*\*(\d+))?', t):
            base=fac.group(1); exp=int(fac.group(2)) if fac.group(2) else 1
            if base=='X': xp=exp
            elif base=='Y': yp=exp
            else: sp_=exp
        lead=re.match(r'(\d+)', t)
        coeff=sign*int(lead.group(1)) if lead else sign
        key=(xp,yp,sp_); monos[key]=monos.get(key,0)+coeff
    return sp.Poly({k:v for k,v in monos.items()}, symX,symY,syms, domain=sp.ZZ)

Gpoly=parse_poly('nonpal_G_clean.txt',X,Y,s)
Spoly=parse_poly('nonpal_S_clean.txt',X,Y,s)
print(f"G: deg={Gpoly.total_degree()} terms={len(Gpoly.terms())}", flush=True)
print(f"S: deg={Spoly.total_degree()} terms={len(Spoly.terms())}", flush=True)

# center_lift in C,D,sigma (from GPT system.py)
def center_lift_expr(C,D,sigma):
    den=C+D-1
    a3=1-C+sigma*C**2*(1-C)/den
    a6=1-D+sigma*D**2*(1-D)/den
    a2=1-a3+sigma*a3**2*(1-a3)/(a3+C-1)
    a7=1-a6+sigma*a6**2*(1-a6)/(a6+D-1)
    return a2,a3,a6,a7

a2,a3,a6,a7=center_lift_expr(C,D,sigma)
# E_L, E_R from the i=2, i=7 boundary equations (definition form)
E_L = a2+a3-1 - sigma*a2*(1-a2)
E_R = a6+a7-1 - sigma*a7*(1-a7)
E_L=sp.together(E_L); E_R=sp.together(E_R)
# numerator over common denominator
E_Ln=sp.factor(sp.numer(sp.together(E_L)))
E_Rn=sp.factor(sp.numer(sp.together(E_R)))
print("E_L num factored:", str(E_Ln)[:200], flush=True)
print("E_R num factored:", str(E_Rn)[:200], flush=True)

# (E_L - E_R) should vanish when C=D (palindrome), so divisible by (C-D).
diff = sp.together(E_L - E_R)
diffn = sp.numer(diff)
diffd = sp.denom(diff)
# divide by (C-D)
q,r = sp.div(sp.expand(diffn), C-D, C, D)
print("diff/(C-D) remainder (should be 0):", sp.simplify(r), flush=True)
quot = sp.expand(q)  # symmetric in C,D? check
# check symmetry: swap C<->D, should be invariant (since (E_L-E_R) is antisymmetric, /(C-D) symmetric)
quot_swap = sp.expand(quot.subs({C:D,D:C}))
print("quot symmetric under C<->D? (diff):", sp.simplify(quot-quot_swap)==0, flush=True)

# Convert quot (symmetric in C,D) to (X=C+D, Y=CD). Use resultant/substitution.
# Express quot in terms of X,Y via C+D=X, CD=Y: replace using D=X-C then it's a poly in C with
# coefficients in X; symmetrize. Easier: use sympy's symmetrize.
sym = sp.symmetrize(quot, (C,D), formal=True)[0]
# sym is in terms of elementary symmetric polys e1=C+D, e2=CD
e1,e2=sp.symbols('e1 e2')
symXY = sp.expand(sym.subs({e1:X, e2:Y}))
print("quot in (X,Y,sigma) degree check; sample terms:", str(symXY)[:200], flush=True)

# Now compare symXY to Gpoly*(X-1)^k. Find k by evaluating at random points.
import random
random.seed(1)
Gexpr=Gpoly.as_expr(); Sexpr=Spoly.as_expr()
rng_pts=[]
for _ in range(8):
    Cv=random.uniform(0.3,0.9); Dv=random.uniform(0.3,0.9)
    if abs(Cv-Dv)<0.05 or Cv+Dv<=1: continue
    sv=random.uniform(0.1,2.0)
    Xv=Cv+Dv; Yv=Cv*Dv
    g=float(Gexpr.subs({X:Xv,Y:Yv,s:sv}))
    qv=float(symXY.subs({X:Xv,Y:Yv,s:sv}))
    if abs(g)<1e-6: continue
    ratio=qv/g
    xm1=Xv-1
    rng_pts.append((Xv,Yv,sv,ratio,xm1))
    print(f"  C={Cv:.3f} D={Dv:.3f} sig={sv:.3f}: quot/G={ratio:.6g}, (X-1)={xm1:.4f}, ratio/(X-1)^k?", flush=True)
# determine k: ratio_i / ratio_j = ((X_i-1)/(X_j-1))^k
if len(rng_pts)>=2:
    X1,_,_,r1,x1=rng_pts[0]; X2,_,_,r2,x2=rng_pts[1]
    if x1!=x2 and x2!=1:
        import math
        k=math.log(r1/r2)/math.log(x1/x2)
        print(f"  estimated k = {k:.4f}", flush=True)
