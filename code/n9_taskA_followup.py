#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Task A follow-up: confirm Ψ_45 roots give F_H=F_K=0 with P_9,7=9; and L=7 P at p=0.4."""
import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar, brentq
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

r=sp.symbols('r',positive=True); z=sp.symbols('z')
psi45_c=[1876133896875,33770410143750,276152284968750,1377167753278125,4744913498227503,
12209430389074209,25009654827225285,44203164300791226,78325378780461003,168594440802523688,
436884674861909052,1131480825950468961,2608520200049103096,5210667018497202633,9082196476378635483,
13975080411142541487,19073279249731452069,23035584713813137782,24497422734518138141,
23019959193342628671,19421778899959715124,15236926325884110258,11283752514161064672,
7567481348377939224,3672293077437648912,39265012871584416,-2638304982140416128,-3694526160829148560,
-3440201891407263072,-2463061427097282720,-1475515623489448704,-739773827916903936,
-324724071701567232,-121165277323960320,-40426609796061696,-11339762983610112,-2907049009348608,
-600362233108992,-118526735204352,-16711740696576,-2547862401024,-181825007616,-23174381568,
969117696,-22929408,25600000]
Psi45=sum(psi45_c[i]*r**(45-i) for i in range(46))
FH=(-81*r**4*z**2+27*r**3*z**4-36*r**3*z**3-99*r**3*z**2+108*r**3*z+33*r**2*z**4-74*r**2*z**3
    -101*r**2*z**2+177*r**2*z-35*r**2+12*r*z**4-76*r*z**3-4*r*z**2+138*r*z-70*r+6*z**4-38*z**3
    -2*z**2+69*z-35)
FK=(9*r**11*z**7+6*r**10*z**8+6*r**10*z**7-7*r**10*z**6+6*r**9*z**8+6*r**9*z**7-7*r**9*z**6
    +9*r**2*z**13-54*r**2*z**12+81*r**2*z**11+63*r**2*z**10-225*r**2*z**9+54*r**2*z**8+216*r**2*z**7
    -135*r**2*z**6-90*r**2*z**5+90*r**2*z**4-9*r**2*z+2*r*z**13-22*r*z**12+70*r*z**11-47*r*z**10
    -127*r*z**9+198*r*z**8+39*r*z**7-218*r*z**6+65*r*z**5+90*r*z**4-52*r*z**3+6*r*z**2-9*r*z+5*r
    +2*z**13-22*z**12+70*z**11-47*z**10-127*z**9+198*z**8+39*z**7-218*z**6+65*z**5+90*z**4
    -52*z**3+6*z**2-9*z+5)
rroots=sorted([float(x) for x in sp.nroots(Psi45) if abs(sp.im(x))<1e-6 and sp.re(x)>0])
print("Ψ_45 positive roots:",rroots)

def HB_w(zr):
    a=zr; b=1-zr**2; c=zr**2-zr; disc=b*b-4*a*c
    if disc<0: return None
    wp=(-b+disc**0.5)/(2*a); wm=(-b-disc**0.5)/(2*a)
    return wp if wp>0 else wm

def recover(zr, rv):
    w=HB_w(zr)
    if w is None or w<=0: return None
    v=zr*(1-w)
    if v<=0 or v>=1: return None
    den=zr*(zr**2+zr-1-zr*w)
    if abs(den)<1e-15: return None
    u=(2*zr**2-1)/den
    if u<=0 or u>=1: return None
    a3=1;a4=1-u;a5=1-v*a4;a6=1-w*a5;a7=1-zr*a6
    if min(a4,a5,a6,a7)<=0: return None
    P=rv**2/(1+rv)*u/(1-v)*(3+a4+a5+a6+a7)
    # closure rhs
    rhs=u**7*(1-w)*(1-zr)/((1-v)**6*a5*a6*a7)
    return dict(z=zr,w=w,u=u,v=v,P=P,rhs=rhs,a=(a4,a5,a6,a7))

print("\n-- at each Ψ_45 root, find z with F_H=F_K=0, check P=9 --")
for rv in rroots:
    fh=lambda zz: float(FH.subs([(r,rv),(z,zz)]))
    fk=lambda zz: float(FK.subs([(r,rv),(z,zz)]))
    # grid search for common root (both small)
    zs=np.linspace(1e-4,0.9999,8000); best=None
    for zz in zs:
        v=abs(fh(zz))+abs(fk(zz))
        if best is None or v<best[1]: best=(zz,v)
    # refine
    zz0=best[0]
    res=minimize_scalar(lambda zz:(abs(fh(zz))+abs(fk(zz)))**2 if 0<zz<1 else 1e18, bracket=(max(1e-6,zz0-0.02),zz0,min(0.999,zz0+0.02)), method='brent')
    zr=res.x
    fhv=fh(zr); fkv=fk(zr)
    print(f"  r={rv:.8f}: z*={zr:.8f}  |F_H|={abs(fhv):.3e}  |F_K|={abs(fkv):.3e}")
    rc=recover(zr, rv)
    if rc:
        print(f"     w={rc['w']:.6f} u={rc['u']:.6f} v={rc['v']:.6f}  P_9,7={rc['P']:.6f} (expect 9)  closure rhs/r^9={rc['rhs']/rv**9:.6f}")
    else:
        print(f"     recover failed (branch/positivity)")

print("\n-- L=7 fixed branch at p=0.4 (r=1.5), find z via closure r^9=rhs, expect P<9 (8.823816) --")
rv=1.5; target=rv**9
zs=np.linspace(1e-3,0.999,20000); best=None
for zz in zs:
    rc=recover(zz,rv)
    if rc is None: continue
    err=abs(rc['rhs']-target)
    if best is None or err<best[1]: best=(zz,err,rc)
zz,err,rc=best
print(f"  best z={zz:.6f} closure_err={err:.4f}  P_9,7={rc['P']:.6f}  (GPT table 8.823816)")
# refine via brentq on (rhs-target)
def g(zz):
    rc=recover(zz,rv)
    return rc['rhs']-target if rc else None
# find sign changes
prev=None; crossings=[]
for zz in zs:
    rc=recover(zz,rv)
    val=rc['rhs']-target if rc else None
    if val is None: prev=None; continue
    if prev is not None and prev*val<0: crossings.append(zz)
    prev=val
print(f"  closure crossings: {len(crossings)}")
for cr in crossings[:3]:
    try:
        zr=brentq(lambda zz: g(zz) if g(zz) is not None else (1 if zz>0.5 else -1), cr-0.001, cr+0.001)
        rc=recover(zr,rv)
        print(f"    z={zr:.8f} P_9,7={rc['P']:.9f}  (expect 8.823816)")
    except Exception as e:
        print(f"    refine failed: {e}")
