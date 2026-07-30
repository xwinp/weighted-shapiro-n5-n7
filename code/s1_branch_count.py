#!/usr/bin/env python3
"""
S1 branch-completeness check at fixed rational p.
S1={0}, x=(0,1,b,c,d,e,f). KKT (Euler lambda=0): support gradients = 0.
Clear denominators -> 5 polynomial equations in b,c,d,e,f.
Compute Groebner basis (lex) -> univariate in one var; count positive real roots.
If exactly 1 positive real stationary point at each sampled p -> no missed branch.
"""
import sympy as sp
from sympy import groebner, Poly, symbols, Rational, solve, nsolve
import mpmath as mp
mp.mp.dps = 30

b,c,d,e,f = symbols('b c d e f', positive=True)
p = Rational(1,4); q = 1-p
x = {0:0, 1:sp.Integer(1), 2:b, 3:c, 4:d, 5:e, 6:f}
P = sum(x[i]/(p*x[(i+1)%7] + q*x[(i+2)%7]) for i in range(7))
grads = [sp.diff(P, v) for v in [b,c,d,e,f]]
# clear denominators
nums = []
for g in grads:
    n = sp.together(g).as_numer_denom()[0]
    nums.append(sp.expand(n))
print("KKT numerators (5 eqs):")
for i,nm in enumerate(nums): print(f"  eq{i} deg_vars:", nm.free_symbols, "total_deg:", sp.total_degree(nm, *[b,c,d,e,f]))

print("\nComputing Groebner basis (lex b>c>d>e>f)... this may take a while", flush=True)
import sys
try:
    G = groebner(nums, b, c, d, e, f, order='lex', domain='QQ')
    print("Groebner basis computed. #polys:", len(G.polys), flush=True)
    # find univariate (in f) element
    uni = [g for g in G.polys if len(g.free_symbols)==1]
    print("univariate elements:", [sp.total_degree(g, *list(g.free_symbols)) for g in uni])
    for g in uni:
        v = list(g.free_symbols)[0]
        poly = Poly(g, v)
        rr = sp.nroots(poly)
        pos = [complex(r).real for r in rr if abs(complex(r).imag)<1e-9 and complex(r).real>1e-12]
        print(f"  univariate in {v}, deg {poly.degree()}: real positive roots = {len(pos)}: {[f'{r:.6f}' for r in pos]}", flush=True)
except Exception as ex:
    import traceback; traceback.print_exc()
    print("Groebner failed/slow:", ex, flush=True)
    # fallback: numerical polynomial system, count solutions via mpmath
    print("fallback: nsolve from known init + random multistart")
