#!/usr/bin/env python3
"""Timeboxed Groebner elimination: reduced n=7 KKT (grad P=0 on support {1,2,3,4,6})
   + sum=1 + P=7, eliminate a,b,c,d,e -> univariate in p."""
import sympy as sp
import signal, sys

a,b,c,d,e,p = sp.symbols('a b c d e p', positive=True)
q = 1 - p
P = a/(p*b+q*c) + b/(p*c+q*d) + c/(p*d) + d/(q*e) + e/(q*a)
grads = [sp.diff(P,v) for v in (a,b,c,d,e)]

# clear denominators -> polynomial numerators
def numify(expr):
    return sp.together(expr).as_numer_denom()[0]

polys = [sp.expand(numify(g)) for g in grads]
polys.append(sp.expand(a+b+c+d+e-1))
polys.append(sp.expand(numify(P-7)))

print("Polynomial system (cleared denominators):")
for i,pp in enumerate(polys):
    print(f"  f{i} = {pp}")
print(f"\nVariables: a,b,c,d,e  (eliminate) ; p (keep)")
sys.stdout.flush()

class TO(Exception): pass
def handler(s,f): raise TO()
signal.signal(signal.SIGALRM, handler)
signal.alarm(240)  # 4-minute hard budget

try:
    G = sp.groebner(polys, [a,b,c,d,e,p], order='lex')
    print("\nGroebner basis computed. Univariate-in-p elements:")
    for g in G:
        if g.free_symbols <= {p}:
            print(f"  [p-only] degree {sp.degree(g,p)}: {g}")
    # also try factor
    for g in G:
        if g.free_symbols <= {p} and sp.degree(g,p) >= 1:
            f = sp.factor(g)
            print(f"  factored: {f}")
except TO:
    print("\n>>> Groebner timed out at 240s (system too heavy for direct elimination).")
except Exception as ex:
    print(f"\n>>> Groebner error: {ex}")
finally:
    signal.alarm(0)
