#!/usr/bin/env python3
"""Analyze the degree-15 factor F(p): irreducibility, real roots in (0,1), confirm a_7,b_7."""
import sympy as sp
import mpmath as mp
mp.mp.dps = 60

p = sp.symbols('p')
F = (5764801*p**15 - 47765494*p**14 + 190003135*p**13 - 486209703*p**12 + 901678743*p**11
     - 1287828143*p**10 + 1464952167*p**9 - 1351039522*p**8 + 1017028633*p**7
     - 624621984*p**6 + 310300032*p**5 - 122238368*p**4 + 36836352*p**3
     - 7952896*p**2 + 1073408*p - 65536)

print("=== factor F over Q ===")
ff = sp.factor(F)
print(f"factored: {ff}")
print(f"irreducible over Q? {sp.Poly(F,p).is_irreducible}")

print("\n=== real roots of F in (0,1) ===")
Fp = sp.Poly(F, p)
# use mpmath to find all real roots
Ff = mp.mpf
coeffs = [Fp.nth(i) for i in range(Fp.degree(), -1, -1)]  # high->low
import mpmath
roots = mpmath.polyroots(coeffs, maxsteps=500, extraprec=20)
real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-40 and 0 < r.real < 1])
print(f"real roots in (0,1): {len(real_roots)}")
a7 = mp.mpf('0.21427352090984096558774231909001838797745454906426773562141322120851430')
b7 = mp.mpf('0.32862767791659196669734085647711086500878730337285459186534900232409311442317')
for i, r in enumerate(real_roots):
    tag = ""
    if abs(r-a7) < 1e-20: tag = "  <-- a_7"
    if abs(r-b7) < 1e-20: tag = "  <-- b_7"
    print(f"  root {i}: {mp.nstr(r, 40)}{tag}")

# verify F(a_7), F(b_7) ~ 0 at high precision
print(f"\nF(a_7) = {mp.nstr(sum(c*a7**(14-i) for i,c in enumerate(coeffs)),6)}")
print(f"F(b_7) = {mp.nstr(sum(c*b7**(14-i) for i,c in enumerate(coeffs)),6)}")

# how many real roots total, and in (0,1/2)?
print(f"\ntotal real roots: {len([r for r in roots if abs(r.imag)<1e-40])}")
print(f"real roots in (0,1/2): {len([r.real for r in roots if abs(r.imag)<1e-40 and 0<r.real<0.5])}")
