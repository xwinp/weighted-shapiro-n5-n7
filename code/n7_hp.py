#!/usr/bin/env python3
"""High-precision n=7 boundary: re-solve KKT at dps 120, then PSLQ deg 4..12 with verification."""
import sympy as sp
import mpmath as mp
mp.mp.dps = 120

a, b, c, d, e, lam, p = sp.symbols('a b c d e lam p', positive=True)
q = 1 - p
P = (a/(p*b+q*c) + b/(p*c+q*d) + c/(p*d) + d/(q*e) + e/(q*a))
Pa, Pb, Pc, Pd, Pe = [sp.diff(P, v) for v in (a, b, c, d, e)]
eqs = [Pa-lam, Pb-lam, Pc-lam, Pd-lam, Pe-lam, a+b+c+d+e-1, P-7]
vars7 = [a, b, c, d, e, lam, p]

def P_expr(av,bv,cv,dv,ev,pv):
    qv=1-pv
    return av/(pv*bv+qv*cv)+bv/(pv*cv+qv*dv)+cv/(pv*dv)+dv/(qv*ev)+ev/(qv*av)

left_init  = [0.17,0.23,0.06,0.31,0.23,6.0,0.22]
right_init = [0.20,0.20,0.08,0.28,0.24,6.0,0.33]

roots = {}
for name, ini in [("a_7", left_init), ("b_7", right_init)]:
    sol = sp.nsolve(eqs, vars7, ini, prec=130, tol=mp.mpf('1e-115'), maxsteps=300)
    pv = mp.mpf(sp.N(sol[6], 120))
    av,bv,cv,dv,ev = [mp.mpf(sp.N(x,120)) for x in sol[:5]]
    Pv = P_expr(av,bv,cv,dv,ev,pv)
    print(f"{name} = {mp.nstr(pv, 70)}")
    print(f"  P check = {mp.nstr(Pv, 30)}")
    roots[name] = pv

a7, b7 = roots["a_7"], roots["b_7"]

print("\n=== PSLQ with verification (dps=120) ===")
for name, t in [("a_7", a7), ("b_7", b7)]:
    print(f"\n{name}:")
    for deg in [4, 5, 6, 7, 8, 9, 10, 12]:
        vec = [t**i for i in range(deg+1)]
        rel = mp.pslq(vec, maxcoeff=10**5, maxsteps=10000)
        if rel:
            resid = sum(rel[i]*t**i for i in range(deg+1))
            print(f"  deg {deg}: coeffs={rel}")
            print(f"        residual={mp.nstr(resid, 8)}")
            # verify on the OTHER root
            other = b7 if name=="a_7" else a7
            resid2 = sum(rel[i]*other**i for i in range(deg+1))
            print(f"        resid at other root={mp.nstr(resid2, 8)}")
            if abs(resid) < mp.mpf('1e-30'):
                print(f"  *** candidate minimal poly deg {deg} ***")
                break
