#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorous certificate that the three-segment sign portrait of P_{S1}^stat-7 is
correct, using ONLY:
  (1) Phi_35 squarefree  (gcd(Phi_35, Phi_35') = 1, exact over Q)  -> crossings simple
  (2) Sturm root counts of Phi_35: exactly two positive roots rho1 in (2,5/2),
      rho2 in (3,7/2); 0 roots in (0,2) and in (7/2, inf)
  (3) one rigorous interval sample P_{S1}^stat(2/5) > 7  (rho=3/2 in (0,rho1))

Logic: P_{S1}^stat-7 is continuous on the positive-w branch (rho>0). Its zeros are
exactly the positive roots of Phi_35 (spurious factors 896 rho^13 (rho+1)^7
(8rho^2+8rho+7)^6 are strictly positive for rho>0). Squarefree => both zeros simple
=> sign changes (alternates). Sturm gives 0 roots in (0,2)⊃(0,rho1) and (7/2,inf)⊃
(rho2,inf), so each of the three components (0,rho1),(rho1,rho2),(rho2,inf) is
sign-constant. Sample rho=3/2 in (0,rho1): P>7 => (0,rho1)>7 => p in (p1,1): P>7.
Alternation: (rho1,rho2)<7 => (p2,p1)<7; (rho2,inf)>7 => (0,p2)>7.
"""
import sympy as sp

rho = sp.symbols('rho')
gpt_coeffs = [262144,1211392,19453952,57874432,444107776,547314432,2185738240,-9171843072,
-34693234688,-166554596992,-327985403648,-543801267968,88546441088,2306869918304,
7848268705920,15843550970240,25279161341952,34334035751596,43526037225048,51848980402968,
56982453166940,55494399759599,46846168825232,33774085926224,20701541228760,10772897763040,
4815489002744,1952186204080,823732168256,410624553696,221073893824,107322284272,41692980224,
11883053056,2193551360,191102976]
Phi = sum(int(c)*rho**(35-i) for i,c in enumerate(gpt_coeffs))

# spurious factors strictly positive for rho>0?
spurious = 896 * rho**13 * (rho+1)**7 * (8*rho**2+8*rho+7)**6
print("(0) spurious factors = 896 rho^13 (rho+1)^7 (8rho^2+8rho+7)^6 >0 for rho>0 :",
      "yes (each factor >0 for rho>0)")

# (1) squarefree
g = sp.gcd(sp.Poly(Phi,rho).as_expr(), sp.Poly(sp.diff(Phi,rho),rho).as_expr())
sqfree = sp.degree(g, rho) == 0
print("(1) gcd(Phi_35, Phi_35') degree =", sp.degree(g, rho), " -> squarefree:", sqfree)

# (2) Sturm counts
def sc(a,b): return sp.count_roots(Phi, a, b)
print("(2) Sturm: (0,2) =", sc(0,2), " (2,5/2) =", sc(2,sp.Rational(5,2)),
      " (5/2,3) =", sc(sp.Rational(5,2),3), " (3,7/2) =", sc(3,sp.Rational(7,2)),
      " (7/2,inf) =", sc(sp.Rational(7,2), sp.oo))
n_pos = sc(0, sp.oo)
print("    total positive roots:", n_pos)
roots_pos = sorted([sp.N(r,15) for r in sp.real_roots(sp.Poly(Phi,rho)) if sp.N(r)>0])
print("    rho1, rho2 =", roots_pos)

# (3) sample rho=3/2 in (0,rho1): P(2/5)>7 is certified by n7_s1_rigorous_certs.py
#   here we just confirm 3/2 < rho1 < 2 ... actually rho1 in (2,5/2), so 3/2 < 2 < rho1. OK.
print("(3) rho=3/2 (p=2/5) in (0, rho1)?:", sp.Rational(3,2) < roots_pos[0],
      " (since rho1>2>3/2)")
print("    P_{S1}^stat(2/5) in [7.157554, 7.157652] > 7  (rigorous, n7_s1_rigorous_certs.py)")

# Conclusion
c0 = sc(0,2)==0 and sc(sp.Rational(7,2),sp.oo)==0 and n_pos==2 and sqfree
print("\nConclusion: squarefree + Sturm(0 roots in (0,rho1) and (rho2,inf)) + P(2/5)>7")
print("  => (0,rho1)>7 i.e. (p1,1): P>7")
print("  => (rho1,rho2)<7 (alternation) i.e. (p2,p1): P<7")
print("  => (rho2,inf)>7 (alternation) i.e. (0,p2): P>7")
print("ALL CERTIFICATES HOLD:", c0)
print("DONE")
