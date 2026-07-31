#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorously assert the four decimal endpoints cited in the paper.

  a_7 = 0.21427352091...   (root of degree-15 F(p), in (1/5, 1/4))
  b_7 = 0.32862767792...   (root of degree-15 F(p), in (1/4, 1/3))
  p_1 = 0.300187927...     ( = 1/(1+rho_1), rho_1 root of Phi_35 in (2, 5/2) )
  p_2 = 0.235880343...     ( = 1/(1+rho_2), rho_2 root of Phi_35 in (3, 7/2) )

Each cited decimal is confirmed by exact Sturm isolation to a rational interval
(lo, hi) of width < 10^-(k+2) (k = cited digits) that ROUNDS to the cited value:
both endpoints lie within 0.5 ulp = 0.5*10^-k of the cited decimal, so the true
root -- which lies in (lo, hi) -- rounds to the cited k-digit decimal.  This is
an exact-rational certificate (Sturm via sympy.count_roots); no floats are used
in the comparison.  exit 0 iff all four decimals are confirmed.
"""
import sys
from fractions import Fraction
import sympy as sp

p = sp.symbols('p')
rho = sp.symbols('rho')

# degree-15 F(p) (from n7_resultant.py: Res_t(R,B) = p^15 (p-1)^6 F)
F = (5764801*p**15 - 47765494*p**14 + 190003135*p**13 - 486209703*p**12
     + 901678743*p**11 - 1287828143*p**10 + 1464952167*p**9 - 1351039522*p**8
     + 1017028633*p**7 - 624621984*p**6 + 310300032*p**5 - 122238368*p**4
     + 36836352*p**3 - 7952896*p**2 + 1073408*p - 65536)

# Phi_35(rho) coefficients (rho^35 .. rho^0), from n7_s1_crossing_resultant.py
gpt_coeffs = [262144,1211392,19453952,57874432,444107776,547314432,2185738240,-9171843072,
-34693234688,-166554596992,-327985403648,-543801267968,88546441088,2306869918304,
7848268705920,15843550970240,25279161341952,34334035751596,43526037225048,51848980402968,
56982453166940,55494399759599,46846168825232,33774085926224,20701541228760,10772897763040,
4815489002744,1952186204080,823732168256,410624553696,221073893824,107322284272,41692980224,
11883053056,2193551360,191102976]
Phi = sum(int(c)*rho**(35-i) for i, c in enumerate(gpt_coeffs))

def sturm_count(poly, x, a, b):
    return sp.count_roots(poly, a, b)

def isolate(poly, x, lo, hi, target_width):
    """Bisect [lo,hi] (exact Sturm) until a single root lies in (lo,hi) with
    hi-lo < target_width.  Returns rational (lo, hi)."""
    assert sturm_count(poly, x, lo, hi) == 1, "not a single-root interval"
    while hi - lo > target_width:
        mid = (lo + hi) / 2
        if sturm_count(poly, x, lo, mid) == 1:
            hi = mid
        else:
            lo = mid
    return lo, hi

def confirms_decimal(lo, hi, cited_frac, k):
    """True iff the root in (lo,hi) rounds to the cited k-digit decimal:
    both endpoints within 0.5*10^-k of cited_frac (exact rational)."""
    half_ulp = Fraction(1, 2) * Fraction(10)**(-k)
    return (abs(lo - cited_frac) < half_ulp) and (abs(hi - cited_frac) < half_ulp)

print("=== Endpoint decimal-confirmation certificate ===\n")

# ---- a_7, b_7 : roots of F in (0,1) ----
print("F(p) real (0,1)-roots:", sturm_count(F, p, 0, 1), "(expect 2)")
# split at 1/4: a_7 in (1/5,1/4), b_7 in (1/4,1/3)
assert sturm_count(F, p, sp.Rational(1,5), sp.Rational(1,4)) == 1
assert sturm_count(F, p, sp.Rational(1,4), sp.Rational(1,3)) == 1
a7_lo, a7_hi = isolate(F, p, sp.Rational(1,5), sp.Rational(1,4), Fraction(1, 10**15))
b7_lo, b7_hi = isolate(F, p, sp.Rational(1,4), sp.Rational(1,3), Fraction(1, 10**15))

# cited decimals (paper abstract / Prop 4.2): a7=0.21427352091 (11), b7=0.32862767792 (11)
a7_cited = Fraction(21427352091, 10**11)   # 0.21427352091
b7_cited = Fraction(32862767792, 10**11)   # 0.32862767792
a7_ok = confirms_decimal(a7_lo, a7_hi, a7_cited, 11)
b7_ok = confirms_decimal(b7_lo, b7_hi, b7_cited, 11)
print("  a_7 in (%.14f, %.14f)  width=%.2e" % (float(a7_lo), float(a7_hi), float(a7_hi-a7_lo)))
print("  a_7 = 0.21427352091 confirmed (both endpoints round to it, 11 digits): %s" % a7_ok)
print("  b_7 in (%.14f, %.14f)  width=%.2e" % (float(b7_lo), float(b7_hi), float(b7_hi-b7_lo)))
print("  b_7 = 0.32862767792 confirmed (both endpoints round to it, 11 digits): %s" % b7_ok)

# ---- rho_1, rho_2 : positive roots of Phi_35; p = 1/(1+rho) ----
print("\nPhi_35 positive roots in (2,5/2):", sturm_count(Phi, rho, 2, sp.Rational(5,2)),
      "  in (3,7/2):", sturm_count(Phi, rho, 3, sp.Rational(7,2)))
assert sturm_count(Phi, rho, 2, sp.Rational(5,2)) == 1
assert sturm_count(Phi, rho, 3, sp.Rational(7,2)) == 1
r1_lo, r1_hi = isolate(Phi, rho, 2, sp.Rational(5,2), Fraction(1, 10**15))
r2_lo, r2_hi = isolate(Phi, rho, 3, sp.Rational(7,2), Fraction(1, 10**15))
# p = 1/(1+rho), rho increasing => p decreasing => p in (1/(1+rho_hi), 1/(1+rho_lo))
p1_lo, p1_hi = 1/(1+r1_hi), 1/(1+r1_lo)
p2_lo, p2_hi = 1/(1+r2_hi), 1/(1+r2_lo)

# cited decimals: p_1=0.300187927 (9), p_2=0.235880343 (9)
p1_cited = Fraction(300187927, 10**9)   # 0.300187927
p2_cited = Fraction(235880343, 10**9)   # 0.235880343
p1_ok = confirms_decimal(p1_lo, p1_hi, p1_cited, 9)
p2_ok = confirms_decimal(p2_lo, p2_hi, p2_cited, 9)
print("  p_1 in (%.13f, %.13f)  width=%.2e" % (float(p1_lo), float(p1_hi), float(p1_hi-p1_lo)))
print("  p_1 = 0.300187927 confirmed (9 digits): %s" % p1_ok)
print("  p_2 in (%.13f, %.13f)  width=%.2e" % (float(p2_lo), float(p2_hi), float(p2_hi-p2_lo)))
print("  p_2 = 0.235880343 confirmed (9 digits): %s" % p2_ok)

# ---- inclusion {p1,p2} subset (a7,b7) (rigorous rational) ----
incl = (a7_hi < p2_lo and p2_hi < b7_lo) and (a7_hi < p1_lo and p1_hi < b7_lo)
print("\n  {p1,p2} subset (a7,b7) (rigorous): %s" % incl)

ok = a7_ok and b7_ok and p1_ok and p2_ok and incl
print("\nALL ENDPOINT DECIMALS CONFIRMED (a7,b7 to 11 digits; p1,p2 to 9 digits): %s" % ok)
print("DONE-ENDPOINTS")
sys.exit(0 if ok else 1)
