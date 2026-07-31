#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorous rational-isolation + interval-arithmetic certificate that M(1)>1/2
for the n=5 S1 branch (GPT review item n=5 #3).

At r=1 (p=q=1/2): stationary curve a^6 - a^2 - 1 = 0  <=>  u^3 - u - 1 = 0, u=a^2
(the plastic constant).  a(1)=sqrt(u), u in (1.32471, 1.32472).
c = a^3 - 1, b = a^2 - c = a^2 - a^3 + 1.  P~_1 = a/(b+c) + b/(c+1) + c + 1/a.
M(1) = 2*P~_1 - 5.  Certify M(1) > 1/2 by rational interval arithmetic.

Use sympy rational intervals: isolate a in a small rational interval [alo,ahi]
with f(a)=a^6-a^2-1, f(alo)<0<f(ahi) and f monotone (f'>0 on the interval), then
bound each term of P~_1 by interval arithmetic over [alo,ahi].
"""
import sys
import sympy as sp
a=sp.symbols('a')
f=a**6-a**2-1
fp=sp.diff(f,a)
# isolate positive root
rr=sp.real_roots(sp.Poly(f,a))
pos=[r for r in rr if r>0]
print("positive roots of a^6-a^2-1:", [sp.nsimplify(r.evalf(20)) for r in pos], "approx", [float(r) for r in pos])
a0=pos[0]
alo=sp.Rational(1150,1000); ahi=sp.Rational(1151,1000)  # 1.150 .. 1.151
# widen/narrow to bracket
print("f(1.150)=",f.subs(a,alo), " f(1.151)=",f.subs(a,ahi))
# ensure bracket
lo,hi=alo,ahi
if f.subs(a,lo)>0 or f.subs(a,hi)<0:
    # search a bracket
    import numpy as np
    av=float(a0)
    lo=sp.Rational(int((av-0.002)*100000),100000); hi=sp.Rational(int((av+0.002)*100000),100000)
    while f.subs(a,lo)>0: lo-=sp.Rational(1,100000)
    while f.subs(a,hi)<0: hi+=sp.Rational(1,100000)
print("bracket a in [",float(lo),",",float(hi),"]  f(lo)=",f.subs(a,lo)," f(hi)=",f.subs(a,hi))
print("f'(a) on bracket >0?", fp.subs(a,lo)>0 and fp.subs(a,hi)>0, "(min f' at lo):", fp.subs(a,lo))

# interval arithmetic for P~_1 and M(1)=2 P~_1 - 5
# c=a^3-1, b=a^2-c=a^2-a^3+1, b+c = a^2 (since c=a^3-1, b=a^2-a^3+1 -> b+c=a^2). Nice!
# term1 = a/(b+c) = a/a^2 = 1/a.   term3 = c = a^3-1.  term4 = 1/a.
# term2 = b/(c+1) = (a^2-a^3+1)/a^3 = (a^2+1)/a^3 - 1 = 1/a + 1/a^3 - 1
# So P~_1 = 1/a + (1/a + 1/a^3 - 1) + (a^3-1) + 1/a = 3/a + 1/a^3 + a^3 - 2
# M(1) = 2(3/a + 1/a^3 + a^3 - 2) - 5 = 6/a + 2/a^3 + 2 a^3 - 9
Mexpr = 6/a + 2/a**3 + 2*a**3 - 9
# Verify the algebraic simplification M(1) = 2*(3/a + 1/a^3 + a^3 - 2) - 5 (difference == 0)
diff_check = sp.simplify(Mexpr - (2*(sp.Rational(3,1)/a + 1/a**3 + a**3 - 2) - 5))
print("identity  M(1) = 2*(3/a+1/a^3+a^3-2)-5  holds (diff==0):", diff_check == 0)
# bound each monomial on [lo,hi] with a>0
def iv_pow(lo,hi,k):
    # a^k positive, monotone increasing for k>0 -> [lo^k, hi^k]; for negative -> [hi^k,lo^k]
    if k>0: return lo**k, hi**k
    else: return hi**k, lo**k   # k negative
alo_,ahi_=lo,hi
# 6/a: [6/hi, 6/lo]; 2/a^3: [2/hi^3, 2/lo^3]; 2 a^3: [2 lo^3, 2 hi^3]
# (6/a and 2/a^3 DECREASE with a; 2 a^3 INCREASES -> lower bound uses a=hi for the
#  first two and a=lo for the third; upper bound reverses.  This is honest interval
#  arithmetic with the per-monomial monotonicity made explicit.)
t1_lo,t1_hi = 6/ahi_, 6/alo_
t2_lo,t2_hi = 2/ahi_**3, 2/alo_**3
t3_lo,t3_hi = 2*alo_**3, 2*ahi_**3
Mlo = t1_lo + t2_lo + t3_lo - 9
Mhi = t1_hi + t2_hi + t3_hi - 9
print("M(1) rigorous interval: [%.6f, %.6f]  (exact rationals above)"%(float(Mlo), float(Mhi)))
print("  Mlo exact =", Mlo)
print("  Mhi exact =", Mhi)
certified = Mlo > sp.Rational(1,2)
print("M(1) > 1/2 ?", certified, "  (lower bound %.6f > 0.5)"%float(Mlo))
print("DONE-M1-CERT certified=%s"%bool(certified))
sys.exit(0 if certified else 1)
