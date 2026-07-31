#!/usr/bin/env python3
"""
n=7 rigorous inactive-KKT proof via Sturm.
Main support {1,2,3,4,6}, x0=x5=0. Reduction (CORRECT b,c):
  a=1, e=t, d=t^2, c=q(q-p t^4)/(p^2 t^2), b=q/(p t)-q^2(q-p t^4)/(p^3 t^2),
  R(p,t)=q^3-p^3 t^5-p^2 q t^8=0 (stationary curve), B(p,t): P=7.
D0 = dP/dx0|x0=x5=0, D5 = dP/dx5|x0=x5=0  (inactive KKT; need >=0 on (a7,b7)).
Eliminate t: G0(p)=Res_t(R, num(D0)), G5(p)=Res_t(R, num(D5)).
Sturm-count roots of G0,G5 in (a7,b7); check sign at band center.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 40

p, t = sp.symbols('p t', positive=True)
q = 1 - p
d_expr = t**2
c_expr = q*(q - p*t**4)/(p**2*t**2)
b_expr = q/(p*t) - q**2*(q - p*t**4)/(p**3*t**2)

# full P with x0,x5 symbolic
x0, x5 = sp.symbols('x0 x5', nonnegative=True)
xx = {0:x0, 1:sp.Integer(1), 2:b_expr, 3:c_expr, 4:d_expr, 5:x5, 6:t}
Pxx = sum(xx[i]/(p*xx[(i+1)%7]+q*xx[(i+2)%7]) for i in range(7))
D0 = sp.simplify(sp.diff(Pxx, x0).subs([(x0,0),(x5,0)]))
D5 = sp.simplify(sp.diff(Pxx, x5).subs([(x0,0),(x5,0)]))

R = q**3 - p**3*t**5 - p**2*q*t**8
D0n = sp.together(D0).as_numer_denom()[0]
D5n = sp.together(D5).as_numer_denom()[0]
D0n = sp.expand(D0n)
D5n = sp.expand(D5n)
print("deg_t R =", sp.degree(R, t), "deg_t D0n =", sp.degree(D0n, t), "deg_t D5n =", sp.degree(D5n, t))

print("Computing G0 = Res_t(R, D0n) ...")
G0 = sp.factor(sp.resultant(R, D0n, t))
print("G0 =", G0)
print("Computing G5 = Res_t(R, D5n) ...")
G5 = sp.factor(sp.resultant(R, D5n, t))
print("G5 =", G5)

# band endpoints: exact rational Sturm isolation of F's two (0,1)-roots (no polyroots)
F = (5764801*p**15 - 47765494*p**14 + 190003135*p**13 - 486209703*p**12
     + 901678743*p**11 - 1287828143*p**10 + 1464952167*p**9 - 1351039522*p**8
     + 1017028633*p**7 - 624621984*p**6 + 310300032*p**5 - 122238368*p**4
     + 36836352*p**3 - 7952896*p**2 + 1073408*p - 65536)
Fp = sp.Poly(F, p, domain=sp.QQ)
n_F_01 = Fp.count_roots(sp.Integer(0), sp.Integer(1))
print(f"\nF: exact count_roots(0,1) = {n_F_01} (expect 2); irreducible over Q? {Fp.is_irreducible}")
ivs = sorted([item[0] for item in Fp.intervals(inf=sp.Integer(0), sup=sp.Integer(1))], key=lambda x: x[0])
# high-precision value of each isolated (0,1)-root via CRootOf.evalf (exact root, not polyroots)
_real = [(r, r.evalf(45)) for r in Fp.all_roots() if r.is_real]
roots01 = sorted([(r, v) for (r, v) in _real if 0 < v < 1], key=lambda x: x[1])
a7 = mp.mpf(roots01[0][1])
b7 = mp.mpf(roots01[1][1])
print(f"a7 = {mp.nstr(a7,25)}  in ({ivs[0][0]}, {ivs[0][1]})")
print(f"b7 = {mp.nstr(b7,25)}  in ({ivs[1][0]}, {ivs[1][1]})")

# strip the (irrelevant) p^k(p-1)^l powers from G0,G5 and Sturm on the remaining factor
def strip_ppow(g):
    g = sp.factor(g)
    # divide out all p and (p-1) factors
    while True:
        r = sp.rem(sp.Poly(g,p), sp.Poly(p,p), p)
        if r.is_zero:
            g = sp.quo(sp.Poly(g,p), sp.Poly(p,p), p).as_expr()
        else:
            break
    while True:
        r = sp.rem(sp.Poly(g,p), sp.Poly(p-1,p), p)
        if r.is_zero:
            g = sp.quo(sp.Poly(g,p), sp.Poly(p-1,p), p).as_expr()
        else:
            break
    return sp.factor(g)

G0c = strip_ppow(G0)
G5c = strip_ppow(G5)
print("\nG0 (stripped) =", G0c)
print("G5 (stripped) =", G5c)
print("deg G0c =", sp.degree(G0c,p), "deg G5c =", sp.degree(G5c,p))

def sturm_count_exact(poly_expr, lo, hi):
    """Exact rational Sturm root count in (lo, hi) via Poly.count_roots (no float truncation)."""
    return sp.Poly(poly_expr, p, domain=sp.QQ).count_roots(sp.Rational(lo), sp.Rational(hi))

# (a7,b7) subset (1/5,1/3); 0 roots in the larger rational interval => 0 in the band
n0 = sturm_count_exact(G0c, sp.Rational(1,5), sp.Rational(1,3))
n5 = sturm_count_exact(G5c, sp.Rational(1,5), sp.Rational(1,3))
print(f"\nSturm roots of G0 in (1/5,1/3): {n0}  (==0: {n0==0})")
print(f"Sturm roots of G5 in (1/5,1/3): {n5}  (==0: {n5==0})")
assert n0 == 0 and n5 == 0, "inactive-KKT Sturm count must be 0"

# sign at band center
pc = (a7+b7)/2
def sgn(g, val):
    num = mp.mpf(sp.N(g.subs(p, sp.Float(val,35)), 35))
    den_dummy = 1
    return num
# D0,D5 sign directly along curve at center
tv = mp.findroot(lambda tt:(1-pc)**3 - pc**3*tt**5 - pc**2*(1-pc)*tt**8, mp.mpf('1.27'))
def fv(expr,pp,tt):
    return mp.mpf(sp.N(expr.subs([(p,sp.Float(pp,35)),(t,sp.Float(tt,35))]),35))
print(f"\nband center p={mp.nstr(pc,15)}: D0={mp.nstr(fv(D0,pc,tv),8)} D5={mp.nstr(fv(D5,pc,tv),8)} (both must be >0)")
print(f"G0c(center) sign = {mp.sign(mp.mpf(sp.N(G0c.subs(p,sp.Float(pc,35)),35)))}")
print(f"G5c(center) sign = {mp.sign(mp.mpf(sp.N(G5c.subs(p,sp.Float(pc,35)),35)))}")
