#!/usr/bin/env python3
"""
Independently verify s2's upgrade of the Tuan-Thuong counterexample:

Claim (s2): keeping the witness vector x=(0,3,4,1,5,0,4) and setting q=1-p
(so the conjectured lower bound stays n/(p+q)=7), the weighted Shapiro sum

    P(p) = sum_i x_i / (p*x_{i+1} + (1-p)*x_{i+2})    (indices mod 7)

satisfies

    P(p) - 7 = -N(p) / [60*p*(p-1)*(3p+1)*(4p-5)],

    N(p) = 5040 p^4 - 7944 p^3 + 2231 p^2 + 113 p - 60,

with two internal real roots p_- ~ 0.219621593897538, p_+ ~ 0.306299360631756,
and P(p) < 7 on (p_-, p_+). Intersecting with the n=5 proven-valid interval
[(5-sqrt5)/10, (5+sqrt5)/10] gives [ (5-sqrt5)/10, p_+ ) where n=5 holds but
n=7 fails.

This script checks all of that with exact sympy rational arithmetic + root
isolation, and also re-verifies the strictly-positive variant
x=(1,300,400,100,500,1,400).
"""
from sympy import symbols, Rational, sqrt, simplify, Poly, real_roots, nsimplify
from fractions import Fraction as F

x = (0, 3, 4, 1, 5, 0, 4)
N = 7
p = symbols('p')

# --- 1. Build P(p) exactly as a sympy rational function of p ---
terms = []
for i in range(N):
    num = x[i]
    den = p * x[(i + 1) % N] + (1 - p) * x[(i + 2) % N]
    terms.append(num / den)
P = sum(terms)
P_minus_7 = simplify(P - 7)
print("P(p) - 7 (simplified) =")
print(P_minus_7)
print()

# --- 2. Factor / compare to s2's claimed closed form ---
numP = simplify(P_minus_7.as_numer_denom()[0])
denP = simplify(P_minus_7.as_numer_denom()[1])
print("Numerator   =", numP)
print("Denominator =", denP)
print()

claimed_num = -5040*p**4 + 7944*p**3 - 2231*p**2 - 113*p + 60
claimed_den = 60*p*(p-1)*(3*p+1)*(4*p-5)
print("s2 claimed num   =", claimed_num)
print("s2 claimed den   =", simplify(claimed_den))
print()
print("num matches s2:", simplify(numP - claimed_num) == 0)
print("den matches s2:", simplify(denP - claimed_den) == 0)
print()

# --- 3. Real roots of the numerator ---
# s2 wrote P(p)-7 = -N(p)/D(p) with N(p)=5040p^4-...; numerator of P-7 is -N(p).
num_poly = Poly(numP, p)
print("Numerator polynomial coeffs (high->low deg):", num_poly.all_coeffs())
roots = real_roots(num_poly)
print("real roots of numerator:", [float(r) for r in roots])
print()

# --- 4. Sign of P(p)-7 on the interval between the two interior roots ---
r_sorted = sorted(roots)
# interior positive roots
pos_roots = [r for r in r_sorted if float(r) > 0 and float(r) < 1]
print("positive interior roots:", [float(r) for r in pos_roots])
if len(pos_roots) >= 2:
    pa, pb = pos_roots[0], pos_roots[1]
    mid = (pa + pb) / 2
    print(f"midpoint p*={float(mid)} -> P-7 = {float(P_minus_7.subs(p, mid))}  (negative means cex)")
    print(f"s2 p_- = {float(pa)}   (claimed 0.219621593897538)")
    print(f"s2 p_+ = {float(pb)}   (claimed 0.306299360631756)")
print()

# --- 5. Intersection with n=5 proven interval ---
lo5 = (5 - sqrt(5)) / 10
hi5 = (5 + sqrt(5)) / 10
print(f"n=5 valid interval = [{float(lo5)}, {float(hi5)}]")
if len(pos_roots) >= 2:
    pb = pos_roots[1]
    inter_lo = lo5
    inter_hi = pb
    print(f"intersection [n=5 holds] x [n=7 fails] = [{float(inter_lo)}, {float(inter_hi)})")
    # check a test point inside the intersection, e.g. p=0.29
    for ptest in [Rational(29,100), lo5, Rational(3,10)]:
        val = P_minus_7.subs(p, ptest)
        print(f"  p={float(ptest):.6f}: P-7 = {float(val):.6e}  (<0 => cex)   in n=5 interval? "
              f"{bool(lo5 <= ptest <= hi5)}")
print()

# --- 6. Re-verify the strictly-positive variant x=(1,300,400,100,500,1,400) ---
xp = (1, 300, 400, 100, 500, 1, 400)
s = F(0)
for i in range(N):
    d = F(3,10)*xp[(i+1) % N] + F(7,10)*xp[(i+2) % N]
    s += F(xp[i]) / d
print(f"strictly-positive variant x={xp}")
print(f"  P = {s} = {float(s):.10f}   <7? {s < 7}   (s2 claimed 527034044794263379/75376890537112230)")
print(f"  s2 fraction equals? {s == F(527034044794263379, 75376890537112230)}")
