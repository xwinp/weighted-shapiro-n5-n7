#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorously verify {p1,p2} = S1-crossings of P=7 lie inside (a7,b7) = S2 failure band,
using exact Sturm isolation with tight rational brackets.

  a7,b7 : two (0,1)-roots of degree-15 F(p).
  p1,p2 : p=1/(1+rho) for the two positive roots rho1,rho2 of Phi_35(rho).
          Phi_35 from n7_s1_crossing_resultant.py (coeffs match GPT's).
We isolate each root in a tight rational interval by exact Sturm sign changes,
then check the four inequalities  a7 < p2 < b7  and  a7 < p1 < b7  rigorously
(interval endpoints are rational, comparisons exact).
"""
import sympy as sp

p = sp.symbols('p')
rho = sp.symbols('rho')

# degree-15 F(p)
F = (5764801*p**15 - 47765494*p**14 + 190003135*p**13 - 486209703*p**12
     + 901678743*p**11 - 1287828143*p**10 + 1464952167*p**9 - 1351039522*p**8
     + 1017028633*p**7 - 624621984*p**6 + 310300032*p**5 - 122238368*p**4
     + 36836352*p**3 - 7952896*p**2 + 1073408*p - 65536)

# Phi_35(rho) coefficients (rho^35 .. rho^0), from n7_s1_crossing_resultant.py / GPT
gpt_coeffs = [262144,1211392,19453952,57874432,444107776,547314432,2185738240,-9171843072,
-34693234688,-166554596992,-327985403648,-543801267968,88546441088,2306869918304,
7848268705920,15843550970240,25279161341952,34334035751596,43526037225048,51848980402968,
56982453166940,55494399759599,46846168825232,33774085926224,20701541228760,10772897763040,
4815489002744,1952186204080,823732168256,410624553696,221073893824,107322284272,41692980224,
11883053056,2193551360,191102976]
Phi = sum(int(c)*rho**(35-i) for i,c in enumerate(gpt_coeffs))

def sturm_count(poly, x, a, b):
    """number of distinct real roots in (a,b) via sympy (exact rational polys)."""
    return sp.count_roots(poly, a, b)

def isolate(poly, x, lo, hi, prec=12):
    """Bisect [lo,hi] using exact Sturm to isolate a single root; return rational (a,b) with a<root<b."""
    # ensure exactly one root in (lo,hi)
    assert sturm_count(poly, x, lo, hi) == 1, (lo, hi, sturm_count(poly, x, lo, hi))
    for _ in range(prec + 8):
        mid = (lo + hi)/2
        nlo = sturm_count(poly, x, lo, mid)
        if nlo == 1:
            hi = mid
        else:
            lo = mid
    return lo, hi

print("Isolating a7, b7 (roots of F in (0,1))...")
# F has exactly two (0,1)-roots. Split at 1/4.
print("  F roots in (0,1/4):", sturm_count(F, p, 0, sp.Rational(1,4)))
print("  F roots in (1/4,1/3):", sturm_count(F, p, sp.Rational(1,4), sp.Rational(1,3)))
print("  F roots in (1/3,1):", sturm_count(F, p, sp.Rational(1,3), 1))
a7_lo, a7_hi = isolate(F, p, sp.Rational(1,5), sp.Rational(1,4), 14)
b7_lo, b7_hi = isolate(F, p, sp.Rational(1,4), sp.Rational(1,3), 14)
print("  a7 in (%.12f, %.12f)" % (float(a7_lo), float(a7_hi)))
print("  b7 in (%.12f, %.12f)" % (float(b7_lo), float(b7_hi)))

print("\nIsolating rho1, rho2 (positive roots of Phi_35)...")
print("  Phi roots in (2,5/2):", sturm_count(Phi, rho, 2, sp.Rational(5,2)))
print("  Phi roots in (3,7/2):", sturm_count(Phi, rho, 3, sp.Rational(7,2)))
r1_lo, r1_hi = isolate(Phi, rho, 2, sp.Rational(5,2), 14)
r2_lo, r2_hi = isolate(Phi, rho, 3, sp.Rational(7,2), 14)
print("  rho1 in (%.12f, %.12f)" % (float(r1_lo), float(r1_hi)))
print("  rho2 in (%.12f, %.12f)" % (float(r2_lo), float(r2_hi)))

# p = 1/(1+rho). rho increasing => p decreasing. So rho1 in (2,5/2) gives p1 = 1/(1+rho1) in
#   (1/(1+5/2), 1/(1+2)) = (2/7, 1/3) = (0.2857, 0.3333).
# rho2 in (3,7/2) gives p2 in (1/(1+7/2), 1/(1+3)) = (2/9, 1/4) = (0.2222, 0.25).
# Tighten: p1 = 1/(1+rho1) with rho1 in (r1_lo,r1_hi) => p1 in (1/(1+r1_hi), 1/(1+r1_lo)).
p1_lo = 1/(1+r1_hi); p1_hi = 1/(1+r1_lo)
p2_lo = 1/(1+r2_hi); p2_hi = 1/(1+r2_lo)
print("  p1 = 1/(1+rho1) in (%.12f, %.12f)" % (float(p1_lo), float(p1_hi)))
print("  p2 = 1/(1+rho2) in (%.12f, %.12f)" % (float(p2_lo), float(p2_hi)))

print("\nRigorous inclusion check {p1,p2} subset (a7,b7):")
c1 = a7_hi < p2_lo        # a7 < p2
c2 = p2_hi < b7_lo        # p2 < b7
c3 = a7_hi < p1_lo        # a7 < p1
c4 = p1_hi < b7_lo        # p1 < b7
print("  a7 < p2 :  a7_hi=%.10f < p2_lo=%.10f  -> %s" % (float(a7_hi), float(p2_lo), c1))
print("  p2 < b7 :  p2_hi=%.10f < b7_lo=%.10f  -> %s" % (float(p2_hi), float(b7_lo), c2))
print("  a7 < p1 :  a7_hi=%.10f < p1_lo=%.10f  -> %s" % (float(a7_hi), float(p1_lo), c3))
print("  p1 < b7 :  p1_hi=%.10f < b7_lo=%.10f  -> %s" % (float(p1_hi), float(b7_lo), c4))
print("\nALL FOUR INEQUALITIES HOLD (rigorous):", c1 and c2 and c3 and c4)
print("  => {p1,p2} subset (a7,b7) :", c1 and c2 and c3 and c4)

# also confirm both crossings have p < 1/3 < p0  (p0 in (3/8,2/5))
print("\n  p1 < 1/3 :", p1_hi < sp.Rational(1,3), "  p2 < 1/3 :", p2_hi < sp.Rational(1,3))
print("  => both crossings < 1/3 < 3/8 < p0  :", p1_hi < sp.Rational(1,3) and p2_hi < sp.Rational(1,3))
print("DONE")
