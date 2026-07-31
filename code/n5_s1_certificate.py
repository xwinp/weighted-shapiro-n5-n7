#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Re-derive the n=5 S1 (one-zero {0}, support {1,2,3,4}) stationary branch from
FIRST PRINCIPLES and certify P>5 on it, with NO GPT-supplied formula.

Face: x0=0, x4=1 (homogeneity), vars a=x1,b=x2,c=x3, p+q=1.
  P = a/(p b+q c) + b/(p c+q) + c/p + 1/(q a)
KKT (Euler lambda=0): dP/da=dP/db=dP/dc=0.  Hand-derived reduction:
  eq1: p b+q c = q a^2
  eq2: p c+q   = q^2 a^3 / p        =>  c = (q^2 a^3 - p q)/p^2,  b = q(a^2-c)/p
  eq3 (stationary curve):  q^3 a^6 - p^3 a^2 - p^2 q = 0
Dehomogenize p:q = r:1 (q=1,p=r):  curve  a^6 - r^3 a^2 - r^2 = 0.
Bound is n/(p+q); homog-deg-0 objective  M = (p+q)P - 5 = (r+1)P - 5.
P=5  <=>  M=0.  Certificate: Res_a(curve, num(M)) = prefactor * Q(r);
Q all-positive coeffs => Q>0 for r>0 => no crossing => sign constant
(checked at r=1: P-5>0) => P>5 on the whole branch.
"""
import sympy as sp
from sympy import Poly, resultant, factor, expand, together, Rational

a, r = sp.symbols('a r', positive=True)
p, q = r, 1                      # dehomogenized p:q = r:1
# c, b from eq1,eq2 (q=1, p=r)
c = (q**2 * a**3 - p*q)/p**2          # = (a^3 - r)/r^2
b = q*(a**2 - c)/p                    # = (a^2 - c)/r
# P on the branch (four terms, x4=1)
P = a/(p*b+q*c) + b/(p*c+q) + c/p + 1/(q*a)
P = sp.simplify(P)
M = sp.simplify((p+q)*P - 5)          # homog-deg-0: (r+1)P - 5
print("P (branch) =", P)
print("M = (r+1)P-5 =", M)

# numerator of M (rational in a,r)
Mnum = sp.together(M).as_numer_denom()[0]
Mnum = sp.expand(Mnum)
curve = a**6 - r**3*a**2 - r**2
print("\ncurve =", curve)
print("Mnum  =", Mnum)

# Reduce Mnum modulo the curve (degree in a < 6) to drop the curve factor cleanly
Pcurve = Poly(curve, a, domain=sp.QQ[r])
Pnum   = Poly(Mnum,  a, domain=sp.QQ[r])
rem = sp.rem(Pnum, Pcurve, a)        # remainder of Mnum mod curve
print("Mnum mod curve (remainder, deg_a<6) =", sp.expand(rem.as_expr()))

# Resultant eliminating a
Res = resultant(Poly(curve, a, domain=sp.QQ[r]), Poly(Mnum, a, domain=sp.QQ[r]), a)
Res = sp.expand(Res)
print("\nRes_a(curve, Mnum) expanded.  total deg(r) =", Poly(Res, r, domain=sp.QQ).degree())
print("factored:", factor(Res))

# Pull out the primitive part in r (drop pure power-of-r prefactors and content)
Pr = Poly(Res, r, domain=sp.QQ)
# remove content and r-powers
content = Pr.content()
prim = Poly(sp.quo(Res, content, r), r, domain=sp.QQ)
# factor out the largest r^k
k = 0
coeffs = prim.all_coeffs()  # high->low; r^k prefactor shows as k low-degree zeros
while coeffs and coeffs[-1] == 0:
    coeffs = coeffs[:-1]; k += 1
prim = Poly(coeffs, r, domain=sp.QQ) if coeffs else prim
print("\nprimitive Q(r): deg", prim.degree(), " content_k=", k)
print("Q coeffs (high->low deg):", prim.all_coeffs())
print("ALL POSITIVE?", all(c > 0 for c in prim.all_coeffs()))
print("ALL SAME SIGN (positive)?", all(c > 0 for c in prim.all_coeffs()))

# sign check at r=1 (p=q): need a solving curve, compute P
import mpmath as mp
mp.mp.dps = 40
def curve_pos_a(rv):
    # a^6 - r^3 a^2 - r^2 = 0, let u=a^2: u^3 - r^3 u - r^2 = 0, unique positive u
    f = lambda u: u**3 - rv**3*u - rv**2
    # bisection on (0, hi): f(0)=-r^2<0, f(hi)>0 -> unique real positive root (no complex)
    lo, hi = mp.mpf(0), mp.mpf(max(1.0, rv**1.5 + rv + 1.0))
    while f(hi) <= 0:
        hi *= 2
    for _ in range(150):
        mid = (lo + hi) / 2
        if f(mid) < 0: lo = mid
        else: hi = mid
    return mp.sqrt((lo + hi) / 2)
a1 = curve_pos_a(1.0)
c1 = (a1**3 - 1)/1**2
b1 = (a1**2 - c1)/1
P1 = a1/(1*b1+1*c1) + b1/(1*c1+1) + c1/1 + 1/(1*a1)
M1 = (1.0+1.0)*P1 - 5   # M=(r+1)P-5 at r=1 (p=q=1/2): the homogeneity-degree-0 target
print("\nr=1 (p=q=1/2): a=%s  P=%s  M=(r+1)P-5=%s  (>0? %s)" % (mp.nstr(a1,15), mp.nstr(P1,15), mp.nstr(M1,8), bool(M1>0)))

# cross-check: numerical M=(r+1)P(r,1)-5 at several r, must be >0 (sign-constant, M(1)>0)
# ALSO verify the stationary point is genuine: dP/da=dP/db=dP/dc ~ 0 on the ORIGINAL face.
print("\nsanity M=(r+1)P-5 across r  +  KKT residual check:")
for rv in [0.05,0.1,0.3,0.5,1.0,2.0,5.0,10.0,50.0,200.0]:
    av = curve_pos_a(rv)
    cv = (av**3 - rv)/rv**2
    bv = (av**2 - cv)/rv
    Pv = av/(rv*bv+cv) + bv/(rv*cv+1) + cv/rv + 1/(av)   # P(r,1), dehomogenized p:q=r:1
    Mv = (rv+1)*Pv - 5                                    # = P_normalized - 5
    # KKT residuals on original face P(a,b,c; p=rv/(1+rv), q=1/(1+rv)), x4=1
    pp = rv/(1+rv); qq = 1/(1+rv)
    dPa = 1/(pp*bv+qq*cv) - 1/(qq*av**2)
    dPb = -av*pp/(pp*bv+qq*cv)**2 + 1/(pp*cv+qq)
    dPc = -av*qq/(pp*bv+qq*cv)**2 - bv*pp/(pp*cv+qq)**2 + 1/pp
    print("  r=%7.2f  M=Pnorm-5=%+9.3e  KKT(|dPa|,|dPb|,|dPc|)=%.1e,%.1e,%.1e" % (
        rv, float(Mv), abs(float(dPa)), abs(float(dPb)), abs(float(dPc))))
print("DONE")
