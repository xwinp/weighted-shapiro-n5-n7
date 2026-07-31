#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify GPT's clean monotonicity certificate for rho(z) on the H_B positive-w branch.

Claim (GPT): on w=w_+(z)>0 (positive root of H_B=zw^2+(1-z^2)w+z^2-z=0),
  K(z)=w z^5/((1-z)D^3) = rho(z)^7,  D=1-z^2+w z^2,
  K'/K = -2 L(w,z) / (w z D H_w),   H_w = 2 z w + 1 - z^2 > 0,
  L(w,z) = (2z^5-3z^4-z^3-2z^2+2) w  -2z^5+3z^4-z^3+5z^2-5z,
  Res_w(H_B, L) = z(z-1) Q_7(z),  Q_7=8z^7-24z^6+20z^5-9z^4+30z^3-15z^2-6.
Q_7 has unique (0,1)-root z7; L=0 lift there is w_L(z7)<0, but w_+(z7)>0,
so L has no zero on the positive branch; sample z=1/2 gives L<0  =>  K'>0, rho'>0.

This script VERIFIES every link:
  (A) re-derive K'/K symbolically and check it equals -2L/(w z D H_w) with GPT's L;
  (B) Res_w(H_B, L) == z(z-1) Q_7  (exact coefficient match);
  (C) w_L(z7) < 0 rigorously (mpmath.iv), w_+(z7) > 0;
  (D) L(w_+(1/2), 1/2) < 0 rigorously.
"""
import sys
import sympy as sp
import mpmath as mp

z, w = sp.symbols('z w', positive=True)

HB = z*w**2 + (1-z**2)*w + z**2 - z
disc = z**4 - 4*z**3 + 2*z**2 + 1          # H_B discriminant in w
wp = (-(1-z**2) + sp.sqrt(disc))/(2*z)     # positive root
D = 1 - z**2 + wp*z**2
K = wp*z**5/((1-z)*D**3)                   # = rho^7

# ---- (A) re-derive K'/K and compare to GPT's L form ----
logK = sp.log(wp) + 5*sp.log(z) - sp.log(1-z) - 3*sp.log(D)
dlogK = sp.diff(logK, z)
# Put over common denominator; express with sqrt(disc)->S to separate the
# rational-numerator structure. We instead directly test GPT's identity:
Lgpt = (2*z**5 - 3*z**4 - z**3 - 2*z**2 + 2)*w + (-2*z**5 + 3*z**4 - z**3 + 5*z**2 - 5*z)
Hw = 2*z*w + 1 - z**2
rhs = -2*Lgpt/(w*z*D*Hw)
# Compare dlogK (with w=wp) to rhs (with w=wp). Both contain sqrt(disc) via wp,D.
diff = sp.simplify(sp.together(dlogK - rhs.subs(w, wp)))
A_ok = (diff == 0)
print("(A) K'/K - (GPT rhs) simplified =", diff, "  => identity holds:", A_ok)

# ---- (B) Res_w(H_B, L) == z(z-1) Q_7 ----
# Use H_B as polynomial in w (deg 2), L linear in w. Resultant wrt w.
Res = sp.resultant(sp.Poly(HB, w), sp.Poly(Lgpt, w), w)
Res = sp.expand(Res)
Q7 = 8*z**7 - 24*z**6 + 20*z**5 - 9*z**4 + 30*z**3 - 15*z**2 - 6
target = sp.expand(z*(z-1)*Q7)
print("(B) Res_w(H_B,L) =", Res)
print("    z(z-1)Q_7    =", target)
print("    coeff match  :", sp.expand(Res - target) == 0)
B_ok = (sp.expand(Res - target) == 0)
# also confirm Q_7 is THE det-certificate Q_7 (same coeffs as rigorous_certs)
Q7_check = 8*z**7-24*z**6+20*z**5-9*z**4+30*z**3-15*z**2-6
print("    Q_7 == det-cert Q_7 :", sp.expand(Q7-Q7_check)==0)
# Sturm: unique root of Q_7 in (0,1)
Q7p = sp.Poly(Q7, z)
print("    Q_7 real roots in (0,1):", [sp.N(r,12) for r in sp.real_roots(Q7p) if 0<r<1])

# ---- (C) w_L(z7) < 0 rigorously, w_+(z7) > 0 ----
mp.mp.ivprec = 80
z7 = sp.nsolve(Q7, z, 0.88, prec=60)
z7f = float(z7)
print("\n(C) z7 = %.12f" % z7f)
# isolate z7 in tight rational interval, verify Q7 brackets
z7lo = mp.mpf(int(z7f*10**14)-200)/mp.mpf(10)**14
z7hi = mp.mpf(int(z7f*10**14)+200)/mp.mpf(10)**14
ziv = mp.iv.mpf([z7lo, z7hi])
def Q7v(t): return 8*t**7-24*t**6+20*t**5-9*t**4+30*t**3-15*t**2-6
qlo, qhi = Q7v(z7lo), Q7v(z7hi)
print("    Q7(z7lo)=%s  Q7(z7hi)=%s  bracket root?"%(mp.nstr(qlo,4),mp.nstr(qhi,4)),
      (qlo<0<qhi) or (qhi<0<qlo))
# w_L(z) = -b/a  with a=(2z^5-3z^4-z^3-2z^2+2), b=(-2z^5+3z^4-z^3+5z^2-5z)
a_iv = 2*ziv**5-3*ziv**4-ziv**3-2*ziv**2+2
b_iv = -2*ziv**5+3*ziv**4-ziv**3+5*ziv**2-5*ziv
print("    a(z7) in [%s,%s]"%(mp.nstr(a_iv.a,6),mp.nstr(a_iv.b,6)))
print("    b(z7) in [%s,%s]"%(mp.nstr(b_iv.a,6),mp.nstr(b_iv.b,6)))
# w_L = -b/a . a,b sign-definite => use direct rigorous interval division.
wL = None
if (a_iv.a>0 or a_iv.b<0) and not (a_iv.a<=0<=a_iv.b):
    wL = (-b_iv)/a_iv     # rigorous interval division (a excludes 0)
    print("    w_L(z7) = -b/a in [%s, %s]  -> %s"%(
        mp.nstr(wL.a,8), mp.nstr(wL.b,8),
        "<0 (NEGATIVE, non-admissible)" if wL.b<0 else (">0?!" if wL.a>0 else "AMBIGUOUS")))
else:
    print("    a sign ambiguous; need tighter isolation")
# positive H_B root w_+(z7)
def iv_sqrt(a):
    return mp.iv.mpf([mp.sqrt(mp.mpf(a.a)), mp.sqrt(mp.mpf(a.b))])
disc_iv = ziv**4-4*ziv**3+2*ziv**2+1
sq = iv_sqrt(disc_iv)
wp_iv = (-(1-ziv**2)+sq)/(2*ziv)
print("    w_+(z7) in [%s, %s]  (>0? %s)"%(mp.nstr(wp_iv.a,8),mp.nstr(wp_iv.b,8), wp_iv.a>0))
C_ok = (wL is not None) and (wL.b < 0) and (wp_iv.a > 0)
print("    w_L(z7) != w_+(z7):", C_ok)

# ---- (D) L(w_+(1/2), 1/2) < 0 rigorously ----
print("\n(D) sample z=1/2")
z0 = mp.mpf(1)/2
d0 = z0**4-4*z0**3+2*z0**2+1
# exact: w_+ = (-3+sqrt(17))/4 ; L = (35 - 5 sqrt 17)/... check sign
import sympy as sp2
zS=sp2.Rational(1,2)
discS = zS**4-4*zS**3+2*zS**2+1   # = 17/16
wpS = (-(1-zS**2)+sp2.sqrt(discS))/(2*zS)   # = (-3+sqrt17)/4
LS = (2*zS**5-3*zS**4-zS**3-2*zS**2+2)*wpS + (-2*zS**5+3*zS**4-zS**3+5*zS**2-5*zS)
LS = sp2.simplify(LS)
D_ok = bool(sp2.N(LS, 20) < 0)
print("    L(w_+(1/2),1/2) exact =", LS, "  numerical =", sp2.N(LS,20))
print("    <0 ?", D_ok)
ok = A_ok and B_ok and C_ok and D_ok
print("\nCERTIFICATE: K'/K identity=%s; Res_w(H_B,L)=z(z-1)Q_7=%s; w_L(z7)<0<w_+(z7)=%s; "
      "L(w_+(1/2),1/2)<0=%s  => rho(z) strictly increasing: %s" % (A_ok, B_ok, C_ok, D_ok, ok))
assert ok, "rho-monotonicity L-certificate failed"
print("DONE-RHO-MONO-LCERT")
sys.exit(0 if ok else 1)
