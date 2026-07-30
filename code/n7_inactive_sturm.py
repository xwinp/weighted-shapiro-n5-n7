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

# band endpoints (roots of F in (0,1))
F = (5764801*p**15 - 47765494*p**14 + 190003135*p**13 - 486209703*p**12
     + 901678743*p**11 - 1287828143*p**10 + 1464952167*p**9 - 1351039522*p**8
     + 1017028633*p**7 - 624621984*p**6 + 310300032*p**5 - 122238368*p**4
     + 36836352*p**3 - 7952896*p**2 + 1073408*p - 65536)
Froots = [mp.mpf(r) for r in mp.polyroots([int(c) for c in sp.Poly(F,p).all_coeffs()], extraprec=60) if abs(mp.im(r))<1e-30]
Froots_real = sorted([mp.re(r) for r in Froots if 0 < mp.re(r) < 1])
a7, b7 = Froots_real[0], Froots_real[1]
print(f"\na7={mp.nstr(a7,25)}\nb7={mp.nstr(b7,25)}")

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

def sturm_count(poly_expr, lo, hi):
    Pp = sp.Poly(poly_expr, p)
    seq = sp.sturm(Pp)
    # convert each sturm term to a mpmath-evaluable coeff list in p
    def coeffs_of(expr):
        Poly = sp.Poly(expr, p)
        d = Poly.degree()
        cl = [mp.mpf(0)]*(d+1)
        for pw, co in Poly.as_dict().items():
            cl[pw[0]] = mp.mpf(int(co))
        return cl
    seq_c = [coeffs_of(tt) for tt in seq]
    def evalc(cl, val):
        v = mp.mpf(0); xv = mp.mpf(1)
        for c in cl:
            v += c*xv; xv *= val
        return v
    def sign_chg_at(val):
        s = [mp.sign(evalc(cl, val)) for cl in seq_c]
        return sum(1 for i in range(len(s)-1) if s[i]*s[i+1] < 0)
    return sign_chg_at(lo) - sign_chg_at(hi)

n0 = sturm_count(G0c, a7, b7)
n5 = sturm_count(G5c, a7, b7)
print(f"\nSturm roots of G0 in (a7,b7): {n0}")
print(f"Sturm roots of G5 in (a7,b7): {n5}")

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
