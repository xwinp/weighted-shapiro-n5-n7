#!/usr/bin/env python3
"""Analyze the degree-15 factor F(p): irreducibility, real roots in (0,1), confirm a_7,b_7.

Exact-symbolic throughout (no mpmath.polyroots): root counting and isolation use
SymPy's exact rational Sturm sequences (Poly.count_roots / Poly.intervals), and
irreducibility is certified both over Q and as a single degree-15 factor mod 23.
"""
import sympy as sp

p = sp.symbols('p')
F = (5764801*p**15 - 47765494*p**14 + 190003135*p**13 - 486209703*p**12 + 901678743*p**11
     - 1287828143*p**10 + 1464952167*p**9 - 1351039522*p**8 + 1017028633*p**7
     - 624621984*p**6 + 310300032*p**5 - 122238368*p**4 + 36836352*p**3
     - 7952896*p**2 + 1073408*p - 65536)
Fp = sp.Poly(F, p, domain=sp.QQ)

print("=== factor F over Q ===")
print(f"factored: {sp.factor(F)}")
print(f"irreducible over Q? {Fp.is_irreducible}")

print("\n=== finite-field certificate: F mod 23 ===")
fl = sp.factor_list(F, modulus=23)
print(f"factor_list(mod 23): {fl}")
deg23 = [sp.degree(f, p) for (f, k) in fl[1]]
print(f"single irreducible factor mod 23, degree {deg23[0]}? "
      f"{len(deg23)==1 and deg23[0]==15 and fl[1][0][1]==1}")

print("\n=== exact real roots of F in (0,1) (rational Sturm) ===")
n01 = Fp.count_roots(sp.Integer(0), sp.Integer(1))
print(f"count_roots(0,1) = {n01}  (==2: {n01==2})")
ivs = sorted([item[0] for item in Fp.intervals(inf=sp.Integer(0), sup=sp.Integer(1))], key=lambda x: x[0])
print(f"isolating intervals: {[(str(lo), str(hi)) for (lo, hi) in ivs]}")
_real = [(r, r.evalf(45)) for r in Fp.all_roots() if r.is_real]
roots01 = sorted([(r, v) for (r, v) in _real if 0 < v < 1], key=lambda x: x[1])
print(f"a_7 ~ {sp.N(roots01[0][1], 40)}  in ({ivs[0][0]}, {ivs[0][1]})")
print(f"b_7 ~ {sp.N(roots01[1][1], 40)}  in ({ivs[1][0]}, {ivs[1][1]})")

print(f"\ntotal real roots: {Fp.count_roots()}  (in (0,1/2): {Fp.count_roots(sp.Integer(0), sp.Rational(1,2))})")
print("DONE")
