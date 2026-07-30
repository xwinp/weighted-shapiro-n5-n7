#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent computation of the p0 minimal polynomial (Q_p) by direct
elimination, bypassing GPT's (inverted) P-formula entirely.

Palindromic {0}-face stationary ansatz: x=(0,1,x2,x3,x4,x4,x3,x2,1), gauge x1=1.
Stationarity (reflection symmetry): g2=g3=g4=0 (3 eqns). P=9 crossing (1 eqn).
Eliminate x4,x3,x2 via chained resultants -> univariate T(p) whose roots
include p0. Sturm-count roots in (a9,b9). If exactly the expected crossings,
p0 is rigorous INDEPENDENTLY of GPT's Q97.

Heavy: runs in background. Degrees are high (GPT's Q97 is degree 97 in C)."""
import sympy as sp, sys, io, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

p = sp.symbols('p')
x2,x3,x4 = sp.symbols('x2 x3 x4')
q = 1-p
# palindromic x: indices 0..8
# x0=0, x1=1, x2, x3, x4, x5=x4, x6=x3, x7=x2, x8=1
X = [sp.Integer(0), sp.Integer(1), x2, x3, x4, x4, x3, x2, sp.Integer(1)]
n=9
def Di(i):
    return p*X[(i+1)%n] + q*X[(i+2)%n]
# P = sum_{i=0}^{8} X[i]/Di(i); i=0 term = 0/Di(0)=0
P = sum(X[i]/Di(i) for i in range(n))
# grad: g_k = dP/dX_k = 1/D_k - sum_{i: (i+1)%9=k} p*X[i]/Di(i)^2 - sum_{i:(i+2)%9=k} q*X[i]/Di(i)^2
# But X depends on x2,x3,x4 with chain rule. Easier: differentiate P wrt x2,x3,x4 directly
# (P as function of x2,x3,x4,p). Stationarity = dP/dx2=dP/dx3=dP/dx4=0.
g2 = sp.diff(P, x2)
g3 = sp.diff(P, x3)
g4 = sp.diff(P, x4)
P9 = P - 9

print("clearing denominators...", flush=True)
t0=time.time()
# clear denominators per equation
g2n = sp.together(g2); g2n = sp.fraction(g2n)[0]
g3n = sp.together(g3); g3n = sp.fraction(g3n)[0]
g4n = sp.together(g4); g4n = sp.fraction(g4n)[0]
P9n = sp.together(P9); P9n = sp.fraction(P9n)[0]
g2n=sp.expand(g2n); g3n=sp.expand(g3n); g4n=sp.expand(g4n); P9n=sp.expand(P9n)
print(f"  cleared in {time.time()-t0:.0f}s", flush=True)
for name,poly in [('g2',g2n),('g3',g3n),('g4',g4n),('P9',P9n)]:
    pp=sp.Poly(poly, x2,x3,x4,p)
    print(f"  {name}: total_deg~{poly.count(x2)+poly.count(x3)+poly.count(x4)} terms={len(pp.as_dict())}", flush=True)

print("\neliminating x4 (resultant pairs)...", flush=True)
t0=time.time()
R1 = sp.resultant(g2n, g3n, x4)  # in x2,x3,p
print(f"  R1=Res_x4(g2,g3) done {time.time()-t0:.0f}s deg_check", flush=True)
R2 = sp.resultant(g2n, g4n, x4)
print(f"  R2=Res_x4(g2,g4) done", flush=True)
R3 = sp.resultant(g2n, P9n, x4)
print(f"  R3=Res_x4(g2,P9) done", flush=True)
R1=sp.expand(R1); R2=sp.expand(R2); R3=sp.expand(R3)

print("\neliminating x3...", flush=True)
t0=time.time()
S1 = sp.resultant(R1, R2, x3)  # in x2,p
print(f"  S1=Res_x3(R1,R2) done {time.time()-t0:.0f}s", flush=True)
S2 = sp.resultant(R1, R3, x3)
print(f"  S2=Res_x3(R1,R3) done", flush=True)
S1=sp.expand(S1); S2=sp.expand(S2)

print("\neliminating x2...", flush=True)
t0=time.time()
T = sp.resultant(S1, S2, x2)  # univariate in p
T = sp.expand(T)
print(f"  T=Res_x2(S1,S2) done {time.time()-t0:.0f}s", flush=True)
Tpoly = sp.Poly(T, p)
print(f"\nT(p) degree = {Tpoly.degree()}", flush=True)
# factor out extraneous p-powers, (p-1)-powers
Tpoly = Tpoly
# primitive part
Tprim = sp.primitive(Tpoly.as_expr())[1]
Tprim = sp.Poly(Tprim, p)
print(f"T primitive degree = {Tprim.degree()}", flush=True)
# count roots in (a9,b9)
a9=sp.Rational(708264058366074,10**16) if False else sp.nsimplify(0.0708264058)
b9=sp.nsimplify(0.43388588203)
try:
    n_band = sp.Poly(Tprim,p).count_roots(sp.Rational(7,100), sp.Rational(44,100))
    print(f"roots in (0.07,0.44): {n_band}", flush=True)
except Exception as e:
    print(f"Sturm failed: {e}", flush=True)
# save coefficients
coeffs = Tprim.all_coeffs()
open('paper/_Qp_independent_coeffs.txt','w').write(str([int(c) for c in coeffs]))
print(f"\nsaved {len(coeffs)} coeffs to paper/_Qp_independent_coeffs.txt", flush=True)
print(f"degree {Tprim.degree()}", flush=True)
