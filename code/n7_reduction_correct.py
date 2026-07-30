#!/usr/bin/env python3
"""Re-derive n=7 reduction with HAND-derived b,c (GPT inverted them). No sp.solve."""
import sympy as sp
import mpmath as mp
mp.mp.dps = 50

p,t,b,c,d = sp.symbols('p t b c d', positive=True)
q = 1 - p
# CORRECT (hand-derived, verified): a=1, e=t, d=t^2
d_expr = t**2
c_expr = q*(q - p*t**4)/(p**2*t**2)
b_expr = q/(p*t) - q**2*(q - p*t**4)/(p**3*t**2)

# P on support with a=1
P = 1/(p*b+q*c) + b/(p*c+q*d) + c/(p*d) + d/(q*t) + t/q
# gradients wrt b,c,d (d symbolic before subs)
Pd = P.subs([(b,b_expr),(c,c_expr),(d,d_expr)])
# curve R: from dP/dd=0 (use full P with d symbolic)
Pb = sp.diff(P, b).subs([(b,b_expr),(c,c_expr),(d,d_expr)])
Pc = sp.diff(P, c).subs([(b,b_expr),(c,c_expr),(d,d_expr)])
Pd_full = sp.diff(P, d).subs([(b,b_expr),(c,c_expr),(d,d_expr)])
# each should vanish identically on the curve; the curve eq is the numerator of (Pb or Pc or Pd) after simplification
def numify(e): return sp.expand(sp.together(e).as_numer_denom()[0])
Rb, Rc, Rd = numify(Pb), numify(Pc), numify(Pd_full)
print("R from dP/db (num):", sp.factor(Rb))
print("R from dP/dc (num):", sp.factor(Rc))
print("R from dP/dd (num):", sp.factor(Rd))

# boundary P=7
B = numify(sp.simplify(Pd) - 7)
print("B (P=7, num):", sp.factor(B))

# verify at a_7
a7=mp.mpf('0.21427352090984096558774231909001838797745454906426773562')
tv=mp.mpf('1.36634936829077')  # = e/a from nsolve
def fv(expr,pp,tt): return mp.mpf(sp.N(expr.subs([(p,sp.Float(pp,45)),(t,sp.Float(tt,45))]),45))
print(f"\nverify a_7: b={mp.nstr(fv(b_expr,a7,tv),10)} c={mp.nstr(fv(c_expr,a7,tv),10)} d={mp.nstr(fv(d_expr,a7,tv),10)}")
print(f"  Rb={mp.nstr(fv(Rb,a7,tv),4)} Rc={mp.nstr(fv(Rc,a7,tv),4)} Rd={mp.nstr(fv(Rd,a7,tv),4)} B={mp.nstr(fv(B,a7,tv),4)}")

# Use the CLEAN irreducible curve factor Rc = q^3 - p^3 t^5 - p^2 q t^8.
# (Rd = dP/dd numerator has spurious extra factors (p-1)(p t^4+p-1); using it
#  lets findroot jump to the pt^4=q branch where c->0. Rc is the true stationary curve.)
R = Rc
print("\nUsing R = Rc (clean curve q^3-p^3 t^5-p^2 q t^8).")
# resultant
E = sp.factor(sp.resultant(R, B, t))
print("E(p) = Res_t(R,B) =", E)

# inactive KKT D0, D5 with correct b,c,d
x0,x5 = sp.symbols('x0 x5', nonnegative=True)
xx={0:x0,1:sp.Integer(1),2:b_expr,3:c_expr,4:d_expr,5:x5,6:t}
Pxx = sum(xx[i]/(p*xx[(i+1)%7]+q*xx[(i+2)%7]) for i in range(7))
D0 = sp.simplify(sp.diff(Pxx,x0).subs([(x0,0),(x5,0)]))
D5 = sp.simplify(sp.diff(Pxx,x5).subs([(x0,0),(x5,0)]))
print("\nD0 =", sp.simplify(D0))
print("D5 =", sp.simplify(D5))
import numpy as np
prev=mp.mpf('1.366')
print("\n=== D0,D5 + b,c positivity along R=0 on (a7,b7) ===")
for pp in np.linspace(0.216,0.327,12):
    try:
        tt=mp.findroot(lambda tt: fv(R,pp,tt), prev); prev=tt
        v0=fv(D0,pp,tt); v5=fv(D5,pp,tt); bb=fv(b_expr,pp,tt); cc=fv(c_expr,pp,tt)
        print(f"  p={pp:.4f} t={float(tt):.5f} b={float(bb):+.3f} c={float(cc):+.3f} D0={float(v0):+.3e} D5={float(v5):+.3e}")
    except Exception as ex:
        print(f"  p={pp:.4f} FAILED {ex}")
