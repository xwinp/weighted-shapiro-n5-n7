#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorous three-sample certificate for the H_B sign portrait (GPT "way B").

For each of the three rho-intervals cut by the Phi_35 crossings rho1, rho2 we take
a rational sample rho (rho=3/2, 3, 4  <=>  p=2/5, 1/4, 1/5) and prove the sign of
P_{S1}^stat-7 by rigorous interval arithmetic:

  rho=3/2 in (0,rho1)   -> expect P>7   (segment (p1,1))
  rho=3   in (rho1,rho2) -> expect P<7   (segment (p2,p1))
  rho=4   in (rho2,inf)  -> expect P>7   (segment (0,p2))

The resultant Res_z(A,B)=896 rho^13 (rho+1)^7 (8rho^2+8rho+7)^6 Phi_35 proves that
P_{S1}^stat=7 has NO solution except at rho1,rho2; hence on each OPEN interval the
sign is constant, and one rigorous sample per interval fixes it. This does NOT rely
on alternation / simple-root transversality.

Method (all rigorous, mpmath.iv with NATIVE iv.sqrt which is outward-directed):
  1. Krawczyk isolate the unique z solving the closure  g(z;rho)=0  on (0,1),
       g = rho^7 (1-z) D^3 - w z^5,   w=(z^2-1+sqrt(Delta))/(2z),
       Delta=z^4-4z^3+2z^2+1,  D=1-z^2+w z^2.
  2. Bound  P = 2 rho (1+rho)/z * [3-2z-z^2+w z(1+z)]  over the certified z-box.
  3. Report the rigorous P interval and its sign vs 7.
"""
import mpmath as mp
import sympy as sp

mp.mp.ivprec = 120
IV = mp.iv
mp.mp.prec = 120

zs = sp.symbols('zs')
rhoS = sp.symbols('rhoS', positive=True)
Delta_sym = zs**4 - 4*zs**3 + 2*zs**2 + 1
wexpr = (-(1-zs**2) + sp.sqrt(Delta_sym)) / (2*zs)
Dexpr = 1 - zs**2 + wexpr*zs**2
gexpr = rhoS**7 * (1-zs) * Dexpr**3 - wexpr * zs**5
Pexpr = 2*rhoS*(1+rhoS)/zs * (3 - 2*zs - zs**2 + wexpr*zs*(1+zs))
gprime = sp.diff(gexpr, zs)

# Hand-written iv evaluators (native IV.sqrt -> outward, rigorous)
def iv_w(z):
    # w = (z^2 - 1 + sqrt(Delta)) / (2z)   [since -(1-z^2)=z^2-1]
    Del = z**4 - 4*z**3 + 2*z**2 + 1
    s = IV.sqrt(Del)
    return (z**2 - 1 + s) / (2*z)
def iv_D(z):
    return 1 - z**2 + iv_w(z)*z**2
def iv_g(z, rho):
    D = iv_D(z); w = iv_w(z)
    return rho**7 * (1-z) * D**3 - w * z**5
def iv_P(z, rho):
    w = iv_w(z)
    return 2*rho*(1+rho)/z * (3 - 2*z - z**2 + w*z*(1+z))

# g'(z) via sympy-derivative lambdified to a plain python function on floats,
# then evaluated on iv by re-implementing with iv ops is messy; instead get the
# derivative expression and lambdify with an iv-aware namespace.
_iv_ns = {'sqrt': IV.sqrt, 'mpf': IV.mpf}
gpfunc = sp.lambdify((zs, rhoS), gprime, modules=[_iv_ns, 'mpmath'])
gfunc = sp.lambdify((zs, rhoS), gexpr, modules=[_iv_ns, 'mpmath'])
Pfunc = sp.lambdify((zs, rhoS), Pexpr, modules=[_iv_ns, 'mpmath'])

def iv_mid(b):
    return (b.a + b.b) / 2

def krawczyk_isolate(rho_val, z_lo, z_hi, iters=80):
    """Krawczyk on g(z)=0. Returns (tight iv box, unique_flag).
    unique_flag is True if at some step K(W) subset int(W) (proves a unique root
    in the box at that step, hence in all refinements contained in it)."""
    rho = IV.mpf([rho_val, rho_val])
    W = IV.mpf([z_lo, z_hi])
    unique = False
    for _ in range(iters):
        mid = IV.mpf([iv_mid(W), iv_mid(W)])
        gm = iv_g(mid, rho)
        gpW = gpfunc(W, rho)
        if not (gpW.a > 0 or gpW.b < 0):
            return W, unique
        K = mid - gm / gpW
        # Krawczyk uniqueness: K subset interior(W)
        if K.a > W.a and K.b < W.b:
            unique = True
        lo = max(K.a, W.a); hi = min(K.b, W.b)
        if lo > hi:
            return W, unique
        Wn = IV.mpf([lo, hi])
        if Wn.a == W.a and Wn.b == W.b:
            # no further contraction via Krawczyk; tighten by bisection sign change
            c = iv_mid(W)
            gl = iv_g(IV.mpf([W.a, W.a]), rho)
            if (gl.a <= 0 <= gl.b):
                Wn = IV.mpf([W.a, c])
            else:
                Wn = IV.mpf([c, W.b])
        W = Wn
        if W.b - W.a < mp.mpf('1e-40'):
            break
    return W, unique

# approximate z to seed the Krawczyk bracket (real bisection, robust)
def wnum(z):
    D=z**4-4*z**3+2*z**2+1
    return (-(1-z**2)+mp.sqrt(D))/(2*z)
def Dnum(z,rho):
    w=wnum(z); return 1-z**2+w*z**2
def gnum(z,rho):
    return rho**7*(1-z)*Dnum(z,rho)**3 - wnum(z)*z**5
def seed_bisect(rho, lo, hi):
    flo=gnum(lo,rho); fhi=gnum(hi,rho)
    for _ in range(200):
        mid=(lo+hi)/2; fm=gnum(mid,rho)
        if flo*fm<=0: hi=mid; fhi=fm
        else: lo=mid; flo=fm
    return (lo+hi)/2

samples = [
    (mp.mpf(3)/2, mp.mpf('0.5'), mp.mpf('0.95'), "rho=3/2 (p=2/5)", ">7"),
    (mp.mpf(3),   mp.mpf('0.9'),  mp.mpf('0.999'), "rho=3   (p=1/4)", "<7"),
    (mp.mpf(4),   mp.mpf('0.95'), mp.mpf('0.9999'),"rho=4   (p=1/5)", ">7"),
]

print("Rigorous three-sample H_B sign certificate (GPT way B)")
print("="*64)
all_ok = True
for rho_val, zlo, zhi, label, expect in samples:
    # numerical seed by real bisection
    znum = seed_bisect(rho_val, zlo, zhi)
    # bracket tightly around seed
    w = mp.mpf('1e-6')
    Z, unique = krawczyk_isolate(rho_val, float(znum - w), float(znum + w))
    if not unique:
        # fall back: wider bracket
        Z, unique = krawczyk_isolate(rho_val, float(zlo), float(zhi))
    rho_iv = IV.mpf([rho_val, rho_val])
    # bound P
    P = iv_P(Z, rho_iv)
    sign = ">7" if P.a > 7 else ("<7" if P.b < 7 else "AMBIGUOUS")
    ok = (sign == expect) and unique
    all_ok = all_ok and ok
    print("%s: z in [%.12f, %.12f]  unique=%s"%(label, float(Z.a), float(Z.b), unique))
    print("    P in [%.10f, %.10f]  -> %s  (expect %s)  %s"%(
        float(P.a), float(P.b), sign, expect, "OK" if ok else "FAIL"))

print("="*64)
print("ALL THREE SAMPLES RIGOROUS (native iv.sqrt, Krawczyk-unique z):", all_ok)
print("  => (0,rho1)>7  (p1,1): P>7")
print("  => (rho1,rho2)<7  (p2,p1): P<7")
print("  => (rho2,inf)>7  (0,p2): P>7")
print("DONE")
