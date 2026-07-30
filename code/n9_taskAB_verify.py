#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verify GPT's n=9 Task A+B closure (paper/gpt_reply_n9_taskAB.txt).
Task A: L=7=(2,7) endpoint via Ψ_45 (degree-45) resultant + Sturm.
Task B: L=9 one-zero face, det H_red<0 when P<9; transpose-dual + double-min lemma.
Every polynomial re-derived/sturm-checked locally. GPT inversion-error pattern known.
"""
import numpy as np
import sympy as sp
from scipy.optimize import brentq
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

r = sp.symbols('r', positive=True)
z = sp.symbols('z', real=True)
C = sp.symbols('C', real=True)
B = sp.symbols('B', real=True)

def sturm_count(poly, var, lo, hi):
    return sp.Poly(poly, var).count_roots(lo, hi)

def pos_roots(poly, var):
    return sorted([float(rt) for rt in sp.nroots(poly) if abs(sp.im(rt))<1e-6 and sp.re(rt)>0])

print("="*72); print("TASK A: L=7=(2,7) endpoint — Ψ_45 resultant + Sturm"); print("="*72)

# --- Ψ_45 coefficients (degree 45 .. 0), from GPT (1.12) ---
psi45_c = [
1876133896875,33770410143750,276152284968750,1377167753278125,4744913498227503,
12209430389074209,25009654827225285,44203164300791226,78325378780461003,
168594440802523688,436884674861909052,1131480825950468961,2608520200049103096,
5210667018497202633,9082196476378635483,13975080411142541487,19073279249731452069,
23035584713813137782,24497422734518138141,23019959193342628671,19421778899959715124,
15236926325884110258,11283752514161064672,7567481348377939224,3672293077437648912,
39265012871584416,-2638304982140416128,-3694526160829148560,-3440201891407263072,
-2463061427097282720,-1475515623489448704,-739773827916903936,-324724071701567232,
-121165277323960320,-40426609796061696,-11339762983610112,-2907049009348608,
-600362233108992,-118526735204352,-16711740696576,-2547862401024,-181825007616,
-23174381568,969117696,-22929408,25600000]
assert len(psi45_c)==46, len(psi45_c)
Psi45 = sum(psi45_c[i]*r**(45-i) for i in range(46))
print(f" Ψ_45 degree = {sp.degree(Psi45)}")
n01 = sturm_count(Psi45, r, 0, 1)
n1inf = sturm_count(Psi45, r, 1, 10**9)
print(f" Sturm: #roots in (0,1)={n01}  (GPT 2);  in (1,∞)={n1inf}  (GPT 0)")
pr = pos_roots(Psi45, r)
print(f" positive roots: {pr}")
if len(pr)>=2:
    rm, rp = pr[0], pr[1]
    pm = rm/(1+rm); pp = rp/(1+rp)
    print(f" r_-={rm:.16f} (GPT 0.07650147965082752)  r_+={rp:.16f} (GPT 0.76642830175721873)")
    print(f" p_-=r/(1+r)={pm:.16f} (GPT 0.07106490896384223)  p_+={pp:.16f} (GPT 0.43388588203369836)")
    print(f" I_9,7=({pm:.10f},{pp:.10f})  (numerical face_min earlier: (0.0710649,0.4338859))")

# --- other factors of Res_z(F_H,F_K) positive on (0,∞)? ---
f1 = 18*r**2+25*r+25; f2 = 9*r**2+2*r+2; f3 = 81*r**4+108*r**3+112*r**2+8*r+4
print(f" factor 18r²+25r+25: disc={25**2-4*18*25}>0? roots={pos_roots(f1,r)} (need >0 on (0,∞))")
print(f" factor 9r²+2r+2: disc={4-72}<0 → >0 everywhere ✓")
print(f" factor 81r⁴+108r³+112r²+8r+4: min on (0,∞) via eval:", float(sp.minimum(f3, r, sp.Interval(0, sp.oo))))

# --- J_7 branch monotonicity ---
J7 = 60*z**7-184*z**6+174*z**5-72*z**4+178*z**3-96*z**2-35
print(f" J_7(z) roots in (0,1): {sturm_count(J7,z,0,1)}  (GPT: unique root, neg w-lift)")
# verify: at the (0,1) root, w from H_B is negative
J7r = [float(rt) for rt in sp.nroots(J7) if abs(sp.im(rt))<1e-6 and 0<sp.re(rt)<1]
for zr in J7r:
    # H_B = z w^2 + (1-z^2) w + (z^2-z); solve for w
    a=zr; b=1-zr**2; c=zr**2-zr
    disc=b**2-4*a*c
    ws=[(-b+disc**0.5)/(2*a),(-b-disc**0.5)/(2*a)]
    print(f"   J_7 root z={zr:.6f}: H_B w-roots={ws} (expect the relevant one negative)")

# --- F_H, F_K: verify Ψ_45 roots give common z-root and P_9,7=9 ---
print("\n -- F_H/F_K at Ψ_45 roots: common z? P_9,7=9? --")
FH = (-81*r**4*z**2+27*r**3*z**4-36*r**3*z**3-99*r**3*z**2+108*r**3*z+33*r**2*z**4
      -74*r**2*z**3-101*r**2*z**2+177*r**2*z-35*r**2+12*r*z**4-76*r*z**3-4*r*z**2
      +138*r*z-70*r+6*z**4-38*z**3-2*z**2+69*z-35)
FK = (9*r**11*z**7+6*r**10*z**8+6*r**10*z**7-7*r**10*z**6+6*r**9*z**8+6*r**9*z**7-7*r**9*z**6
      +9*r**2*z**13-54*r**2*z**12+81*r**2*z**11+63*r**2*z**10-225*r**2*z**9+54*r**2*z**8
      +216*r**2*z**7-135*r**2*z**6-90*r**2*z**5+90*r**2*z**4-9*r**2*z
      +2*r*z**13-22*r*z**12+70*r*z**11-47*r*z**10-127*r*z**9+198*r*z**8+39*r*z**7
      -218*r*z**6+65*r*z**5+90*r*z**4-52*r*z**3+6*r*z**2-9*r*z+5*r
      +2*z**13-22*z**12+70*z**11-47*z**10-127*z**9+198*z**8+39*z**7-218*z**6+65*z**5
      +90*z**4-52*z**3+6*z**2-9*z+5)
for rv in pr[:2]:
    fh = sp.Poly(FH.subs(r, sp.Rational(rv).limit(sp.oo,0) if False else rv), z) if False else sp.Poly(sp.nsimplify(FH.subs(r,rv), rational=False), z)
    # just use numeric
    fhn = sp.Poly(FH.subs(r, rv), z, domain='RR')
    fkn = sp.Poly(FK.subs(r, rv), z, domain='RR')
    rg = sp.gcd(fhn, fkn)
    common = sp.nroots(rg) if rg.degree()>0 else []
    common = [float(x) for x in common if abs(sp.im(x))<1e-6 and 0<sp.re(x)<1]
    print(f"  r={rv:.6f}: common z-roots in (0,1) = {common}")
    for zr in common:
        # recover w from H_B (positive root), then u,v, then P
        a_=zr; b_=1-zr**2; c_=zr**2-zr; disc=b_**2-4*a_*c_
        wp=(-b_+disc**0.5)/(2*a_) if disc>=0 else None
        wm=(-b_-disc**0.5)/(2*a_)
        wv = wp if wp and wp>0 else wm
        if wv and wv>0:
            v=zr*(1-wv); u=(2*zr**2-1)/(zr*(zr**2+zr-1-zr*wv))
            a3=1; a4=1-u; a5=1-v*a4; a6=1-wv*a5; a7=1-zr*a6
            P = rv**2/(1+rv)*u/(1-v)*(3+a4+a5+a6+a7)
            print(f"     z={zr:.6f} w={wv:.6f} u={u:.6f} v={v:.6f}  P_9,7={P:.6f} (expect 9)")

# --- L=7 fixed-branch P at p=0.4 (inside, should be <9) ---
print("\n -- L=7 fixed branch P at p=0.4 (expect <9) --")
p=0.4; rv=(1-p)/p
# find z s.t. closure r^9 = u^7(1-w)(1-z)/((1-v)^6 a5 a6 a7) holds with H_B
# scan z in (0,1), compute w,v,u,a's, compute rhs, compare to r^9
zs=np.linspace(0.01,0.99,4000); best=None
for zr in zs:
    a_=zr; b_=1-zr**2; c_=zr**2-zr; disc=b_**2-4*a_*c_
    if disc<0: continue
    wp=(-b_+disc**0.5)/(2*a_); wm=(-b_-disc**0.5)/(2*a_)
    wv=wp if wp>0 else wm
    if wv<=0: continue
    v=zr*(1-wv)
    if v<=0 or v>=1: continue
    den=zr*(zr**2+zr-1-zr*wv)
    if den==0: continue
    u=(2*zr**2-1)/den
    if u<=0 or u>=1: continue
    a3=1;a4=1-u;a5=1-v*a4;a6=1-wv*a5;a7=1-zr*a6
    if min(a4,a5,a6,a7)<=0: continue
    rhs=u**7*(1-wv)*(1-zr)/((1-v)**6*a5*a6*a7)
    P=rv**2/(1+rv)*u/(1-v)*(3+a4+a5+a6+a7)
    if abs(rhs-rv**9)<0.05 and (best is None or abs(rhs-rv**9)<best[1]):
        best=(zr,abs(rhs-rv**9),P)
if best:
    print(f"  p=0.4: z={best[0]:.5f} closure_err={best[1]:.4f} P_9,7={best[2]:.6f} (GPT <9, table 8.823816)")
else:
    print("  p=0.4: no fixed-branch stationary found in scan (check)")

# --- Task A resultant divisibility: Res_z(F_H,F_K) divisible by Ψ_45? ---
print("\n -- Res_z(F_H,F_K) divisibility by Ψ_45 (heavy) --")
try:
    res = sp.resultant(FH, FK, z)
    q_, mod = sp.div(sp.Poly(res, r), sp.Poly(Psi45, r), r)
    print(f"  Res_z degree={sp.degree(res,r)}; divisible by Ψ_45? remainder is zero: {mod.is_zero}")
    # also check the listed prefactor
    pref = 5*r**9*(r+1)**3*(18*r**2+25*r+25)**2*(9*r**2+2*r+2)**5*(81*r**4+108*r**3+112*r**2+8*r+4)
    q2, mod2 = sp.div(sp.Poly(res,r), sp.Poly(pref*Psi45, r), r)
    print(f"  Res_z = pref·Ψ_45·(quotient deg {q2.degree() if q2 else '?'}); exact factorization? remainder zero: {mod2.is_zero}")
except Exception as e:
    print(f"  resultant computation failed/slow: {e}")
