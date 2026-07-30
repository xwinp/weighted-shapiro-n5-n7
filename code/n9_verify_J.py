#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify J = (1/const) * det(3x3) from definition (5.2), numerically.

(5.2): J_crit = det | G_X   S_X   D*N_X - N*D_X |
                     | G_Y   S_Y   D*N_Y - N*D_Y |
                     | G_s   S_s   D*N_s - N*D_s |
where N=rho9_num, D=rho9_den, s=sigma.  GPT: remove common factor sigma^3 -> deg 53.
So J_file should be proportional to det * sigma^(-3)  (i.e. det = const * sigma^3 * J_file)
up to a rational constant.  Find the power of sigma and the constant by sampling."""
import re, sympy as sp, numpy as np
from pathlib import Path
HERE=Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'
Xs,Ys,ss=sp.symbols('X Y s')
def parse_poly(name):
    txt=(HERE/name).read_text().strip()
    terms=re.findall(r'[+-][^+-]+|^[^+-]+', txt)
    monos={}
    for t in terms:
        t=t.strip()
        if not t: continue
        sign=1
        if t.startswith('-'): sign=-1; t=t[1:].strip()
        elif t.startswith('+'): t=t[1:].strip()
        xp=yp=sp_=0
        for fac in re.finditer(r'(X|Y|s)(?:\*\*(\d+))?', t):
            base=fac.group(1); exp=int(fac.group(2)) if fac.group(2) else 1
            if base=='X': xp=exp
            elif base=='Y': yp=exp
            else: sp_=exp
        lead=re.match(r'(\d+)', t)
        coeff=sign*int(lead.group(1)) if lead else sign
        key=(xp,yp,sp_); monos[key]=monos.get(key,0)+coeff
    return sp.Poly({k:v for k,v in monos.items()},Xs,Ys,ss,domain=sp.ZZ)
G=parse_poly('nonpal_G_clean.txt'); S=parse_poly('nonpal_S_clean.txt')
J=parse_poly('nonpal_J_clean.txt'); N=parse_poly('nonpal_rho9_num.txt'); D=parse_poly('nonpal_rho9_den.txt')

# derivatives
Gx,Gy,Gs=G.diff(Xs),G.diff(Ys),G.diff(ss)
Sx,Sy,Ss=S.diff(Xs),S.diff(Ys),S.diff(ss)
Nx,Ny,Ns=N.diff(Xs),N.diff(Ys),N.diff(ss)
Dx,Dy,Ds=D.diff(Xs),D.diff(Ys),D.diff(ss)
# third column entries: D*N_x - N*D_x
Cx=D*Nx-N*Dx; Cy=D*Ny-N*Dy; Cz=D*Ns-N*Ds
# determinant
det = Gx*(Sy*Cz-Cy*Ss) - Gy*(Sx*Cz-Cx*Ss) + Gs*(Sx*Cy-Cx*Sy)
det=sp.expand(det.as_expr()) if hasattr(det,'as_expr') else sp.expand(det)
Jexpr=J.as_expr()
print(f"det total degree ~ checking; J deg={J.total_degree()}", flush=True)

# direct numerical evaluation of a Poly via its term dict
def ev(P,Xv,Yv,sv):
    t=0.0
    for (a,b,c),coef in P.terms():
        t+=float(coef)*(Xv**a)*(Yv**b)*(sv**c)
    return t
# build det as a Poly too
detP=sp.Poly(det,Xs,Ys,ss,domain=sp.ZZ)
print(f"det total_degree={detP.total_degree()} terms={len(detP.terms())}  (J deg={J.total_degree()})", flush=True)

rng=np.random.default_rng(11)
print("\nSampling det / (sigma^k * J) at random interior points:", flush=True)
samples=[]
for _ in range(12):
    Xv=rng.uniform(1.1,1.9)
    ylo=Xv-1+1e-3; yhi=Xv**2/4-1e-3
    if yhi<=ylo: continue
    Yv=rng.uniform(ylo,yhi)
    sv=rng.uniform(0.2,3.0)
    try:
        dv=ev(detP,Xv,Yv,sv); jv=ev(J,Xv,Yv,sv)
    except: continue
    if abs(jv)<1e-6 or abs(dv)<1e-6: continue
    r=dv/jv
    samples.append((Xv,Yv,sv,r))
    print(f"  X={Xv:.4f} Y={Yv:.4f} s={sv:.4f}: det/J = {r:.6g}", flush=True)

# det/J should = const * sigma^k.  Find k from ratio of pairs.
import math
if len(samples)>=2:
    X1,Y1,s1,r1=samples[0]; X2,Y2,s2,r2=samples[1]
    # r = const * s^k  => r1/r2 = (s1/s2)^k
    if r1>0 and r2>0 and s1>0 and s2>0:
        k=math.log(abs(r1/r2))/math.log(s1/s2)
        print(f"\n  estimated power k = {k:.5f}  (GPT claims factor sigma^3, so det/J ~ sigma^3)", flush=True)
        kr=round(k)
        const=r1/(s1**kr)
        good=all(abs(const*(sv**kr)-r)<1e-6*max(1,abs(r)) for _,_,sv,r in samples)
        print(f"  -> k={kr}, const={const:.6g}: {'VERIFIED det = const*sigma^'+str(kr)+'*J' if good else 'MISMATCH'}", flush=True)
        print(f"  J (criticality numerator) {'MATCHES definition (5.2)' if good else 'does NOT match'}", flush=True)
