#!/usr/bin/env python3
"""Verify GPT's CORRECTED S3={0,3} closed form.
x=(0,1,b,0,c,d,e), r=(q/p)^{1/5} (INVERTED from GPT's earlier wrong version),
b=r^2, c=r^{-1}, e=r^3, d=r-r^8=r(1-r^7). KKT: c=(p/q)b^2, e=(p/q)b^4,
d=(p/q)b^3-b^4, p^2 b^5 = q^2. d>0 <=> r<1 <=> p>1/2.
P_S3^stat = (1+r^5)(5-r^7)/r^2, claimed always >7 for 0<r<1.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 40

p, q = sp.symbols('p q', positive=True)
b, c, d, e = sp.symbols('b c d e', positive=True)
# S3 = {0,3}: x0=0,x1=1,x2=b,x3=0,x4=c,x5=d,x6=e
x = {0:0, 1:sp.Integer(1), 2:b, 3:0, 4:c, 5:d, 6:e}
P = sum(x[i]/(p*x[(i+1)%7] + q*x[(i+2)%7]) for i in range(7))
P = sp.simplify(P)
print("P_S3 =", P)

# gradients wrt support vars b,c,d,e (Euler lambda=0 -> grad=0 on support)
grads = {v: sp.simplify(sp.diff(P, v)) for v in [b,c,d,e]}
for v,g in grads.items(): print(f"  dP/d{v} =", g)

# substitute closed form with r, p=1/(1+r^5), q=r^5/(1+r^5)
r = sp.symbols('r', positive=True)
subs = [(b, r**2), (c, 1/r), (e, r**3), (d, r - r**8),
        (p, 1/(1+r**5)), (q, r**5/(1+r**5))]
print("\n-- substitute closed form into gradients (must be 0) --")
for v,g in grads.items():
    gv = sp.simplify(g.subs(subs))
    print(f"  dP/d{v} = {gv}   (zero: {gv==0})")

# P value
Pv = sp.simplify(P.subs(subs))
print("\nP_S3^stat(r) =", sp.factor(Pv))
target = (1+r**5)*(5-r**7)/r**2
print("matches (1+r^5)(5-r^7)/r^2 ?", sp.simplify(Pv - target)==0)

# check >7 for 0<r<1: f(r)-7 = N(r)/r^2, N(r)=(1+r^5)(5-r^7)-7r^2
fr = mp.mpf
def f(rr): return (1+rr**5)*(5-rr**7)/rr**2
def N(rr): return (1+rr**5)*(5-rr**7) - 7*rr**2
import numpy as np
rs = np.linspace(0.001, 0.999, 20000)
vals = [f(fr(rr)) for rr in rs]
imin = int(np.argmin(vals))
print(f"\nscan min of P_S3^stat on (0,1): at r={rs[imin]:.5f}, P={float(vals[imin]):.8f}  (>7: {vals[imin]>7})")
# P is decreasing on (0,1) (scan min at largest r), min = 8 at boundary r=1 (p=1/2).
# Confirm N(r)>0 on (0,1) and monotonicity of P via sign of f'.
rsm = sp.symbols('rsm', positive=True)
Nsym = (1+rsm**5)*(5-rsm**7) - 7*rsm**2
fp = sp.diff((1+rsm**5)*(5-rsm**7)/rsm**2, rsm)   # f'(r)
print("at r=1 (p=1/2): P =", float(f(fr(1))), " N(1) =", float(N(fr(1))))
Nvals = [N(fr(rr)) for rr in rs]
print(f"min of N(r) on (0,1) scan: {float(min(Nvals)):.6e}  (>0: {min(Nvals)>0})")
fpvals = [float(fp.subs(rsm, float(rr))) for rr in rs]
print(f"f'(r) sign on (0,1): min={min(fpvals):.4e} max={max(fpvals):.4e}  (all<0 -> P decreasing -> min at r=1=8)")
# Rigorous: Sturm count of N' roots in (0,1) to certify no interior min below boundary
Np = sp.Poly(sp.diff(Nsym, rsm), rsm)
print("N'(r) roots in (0,1):", Np.count_roots(0, 1), " (0 => N monotonic-ish, min on boundary)")
print("N(0) =", float(Nsym.subs(rsm,0)), " N(1) =", float(Nsym.subs(rsm,1)), " => N>0 on (0,1)")

# confirm d>0 <=> p>1/2, and band (p in (a7,b7) ~0.21-0.33, all <1/2) -> r>1 -> d<0 -> no positive stat
a7,b7 = mp.mpf('0.214273520909841'), mp.mpf('0.328627677916592')
for pp in [a7, mp.mpf('0.25'), b7, mp.mpf('0.5'), mp.mpf('0.6'), mp.mpf('0.7')]:
    qq = 1-pp; rr = (qq/pp)**mp.mpf('0.2')
    dd = rr - rr**8
    print(f"  p={float(pp):.4f} r={float(rr):.4f} d=r-r^8={float(dd):+.5f} {'(positive stat exists)' if dd>0 else '(no positive stat)'}")
