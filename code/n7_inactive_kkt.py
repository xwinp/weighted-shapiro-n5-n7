#!/usr/bin/env python3
"""
n=7 inactive KKT verification (main support {1,2,3,4,6}, x0=x5=0).
Recompute dP/∂x0, dP/∂x5 from scratch (do NOT trust GPT's formula).
Substitute the reduction (a=1, e=t, d=t^2, c=p^2 t^2/[q(q-p t^4)],
 b=pt/q - p^3 t^2/[q^2(q-p t^4)]) and check D0,D5 >= 0 along R(p,t)=0 on (a_7,b_7).
Numerical check + resultant/Sturm setup for the proof.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 40

# ---- full P with x0=x5=0, vars a=x1,b=x2,c=x3,d=x4,e=x6 ----
a,b,c,d,e,p,t = sp.symbols('a b c d e p t', positive=True)
q = 1 - p
# x = (0, a, b, c, d, 0, e)
x = {0:0, 1:a, 2:b, 3:c, 4:d, 5:0, 6:e}
def xi(i): return x[i]
Pfull = sum(xi(i)/(p*xi((i+1)%7)+q*xi((i+2)%7)) for i in range(7))

# inactive derivatives wrt x0 and x5 (treat x0,x5 as separate symbols momentarily)
x0,x5 = sp.symbols('x0 x5', nonnegative=True)
xx = {0:x0, 1:a, 2:b, 3:c, 4:d, 5:x5, 6:e}
Pxx = sum(xx[i]/(p*xx[(i+1)%7]+q*xx[(i+2)%7]) for i in range(7))
dP0 = sp.diff(Pxx, x0).subs([(x0,0),(x5,0)])
dP5 = sp.diff(Pxx, x5).subs([(x0,0),(x5,0)])
dP0 = sp.simplify(dP0)
dP5 = sp.simplify(dP5)
print("dP/∂x0 (x0=x5=0) =", dP0)
print("dP/∂x5 (x0=x5=0) =", dP5)

# reduction substitutions
q1 = 1 - p
c_expr = p**2*t**2/(q1*(q1 - p*t**4))
b_expr = p*t/q1 - p**3*t**2/(q1**2*(q1 - p*t**4))
d_expr = t**2
e_expr = t
a_expr = 1
subs = [(a,1),(b,b_expr),(c,c_expr),(d,d_expr),(e,e_expr)]
D0 = sp.simplify(dP0.subs(subs))
D5 = sp.simplify(dP5.subs(subs))
print("\nD0(p,t) =", D0)
print("D5(p,t) =", D5)

R = q1**3 - p**3*t**5 - p**2*q1*t**8

# ---- numerical check along the curve on (a_7, b_7) ----
# use mpmath: for given p, solve R=0 for t near ~1.36, eval D0,D5
a7 = mp.mpf('0.21427352090984096558774231909001838797745455')
b7 = mp.mpf('0.32862767791659196669734085647711086500878730')
def Dt(pp):
    pv=mp.mpf(pp); qv=1-pv
    # solve R=0 in t (use mpmath findroot near previous)
    f=lambda tv: qv**3 - pv**3*tv**5 - pv**2*qv*tv**8
    tv=mp.findroot(f, mp.mpf('1.36'))
    return tv
def evalD(D, pp, tv):
    pv=mp.mpf(pp); qv=1-pv
    return mp.mpf(sp.N(D.subs([(p,sp.Float(pv,40)),(t,sp.Float(tv,40))]),40))
print("\n=== numerical D0,D5 along curve on (a_7,b_7) ===")
import numpy as np
prev_t = mp.mpf('1.3663')
for pp in np.linspace(0.215, 0.328, 12):
    try:
        tv = mp.findroot(lambda tv:(1-pp)**3 - pp**3*tv**5 - pp**2*(1-pp)*tv**8, prev_t)
        prev_t = tv
        v0 = evalD(D0, pp, tv); v5 = evalD(D5, pp, tv)
        print(f"  p={pp:.4f} t={float(tv):.6f}  D0={float(v0):+.4e}  D5={float(v5):+.4e}")
    except Exception as ex:
        print(f"  p={pp:.4f} FAILED {ex}")

# ---- resultant for proof: eliminate t between R and D0 (and D5) ----
print("\n=== resultant Res_t(R, numerator(D0)) for Sturm proof ===")
D0n = sp.together(D0).as_numer_denom()[0]
D5n = sp.together(D5).as_numer_denom()[0]
G0 = sp.resultant(R, sp.expand(D0n), t)
G5 = sp.resultant(R, sp.expand(D5n), t)
G0 = sp.factor(G0); G5 = sp.factor(G5)
print("G0(p) =", G0)
print("G5(p) =", G5)
