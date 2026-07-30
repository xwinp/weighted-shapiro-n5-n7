#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test whether {G=S=0} characterizes the non-palindromic stationary curve {E_L=E_R=0, C!=D}.

Method:
 1. Find points on E_L=E_R=0 by picking sigma, solving 2 eqns in (C,D) numerically.
    At each solution check G(X,Y,sigma)=0 and S=0  (forward: stationary => G=S=0).
 2. Find points on G=S=0 by picking sigma, solving G=S=0 in (X,Y) [C,D from X,Y],
    then check E_L=E_R=0  (backward: G=S => stationary).
If both hold, {G,S} correctly characterizes the curve and can be trusted for the
empty-set certificate."""
import re, sympy as sp, numpy as np
from scipy.optimize import root
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
Gpoly=parse_poly('nonpal_G_clean.txt'); Spoly=parse_poly('nonpal_S_clean.txt')
Gf=sp.lambdify((Xs,Ys,ss),Gpoly.as_expr(),'numpy')
Sf=sp.lambdify((Xs,Ys,ss),Spoly.as_expr(),'numpy')

C,D,sigma=sp.symbols('C D sigma')
def center_lift(C,D,sigma):
    den=C+D-1
    a3=1-C+sigma*C**2*(1-C)/den
    a6=1-D+sigma*D**2*(1-D)/den
    a2=1-a3+sigma*a3**2*(1-a3)/(a3+C-1)
    a7=1-a6+sigma*a6**2*(1-a6)/(a6+D-1)
    return a2,a3,a6,a7
a2,a3,a6,a7=center_lift(C,D,sigma)
E_L = a2+a3-1 - sigma*a2*(1-a2)
E_R = a6+a7-1 - sigma*a7*(1-a7)
ELf=sp.lambdify((C,D,sigma),sp.sympify(E_L),'numpy')
ERf=sp.lambdify((C,D,sigma),sp.sympify(E_R),'numpy')

print("=== Forward: E_L=E_R=0  =>  G=S=0 ? ===", flush=True)
rng=np.random.default_rng(3)
fwd_ok=0; fwd_bad=0
for sig in [0.3,0.6,1.0,1.5,2.0,3.0]:
    # solve E_L=E_R=0 in (C,D) with C,D in (0,1), C+D>1, C!=D
    found=False
    for _ in range(200):
        c0=rng.uniform(0.2,0.98); d0=rng.uniform(0.2,0.98)
        if abs(c0-d0)<0.05 or c0+d0<=1.02: continue
        def F(v):
            c,d=v
            try:
                return [float(ELf(c,d,sig)), float(ERf(c,d,sig))]
            except: return [1e10,1e10]
        try:
            sol=root(F,[c0,d0],method='hybr',options={'xtol':1e-14,'maxfev':5000})
        except: continue
        if sol.success and np.linalg.norm(F(sol.x))<1e-10:
            c,d=sol.x
            if 0<c<1 and 0<d<1 and c+d>1 and abs(c-d)>1e-4:
                Xv=c+d; Yv=c*d
                gv=float(Gf(Xv,Yv,sig)); sv=float(Sf(Xv,Yv,sig))
                # scale: G,S are large; compare relative to local scale
                gok = abs(gv)<1e-6; sok=abs(sv)<1e-6
                tag = "OK" if (gok and sok) else "MISMATCH"
                print(f"  sig={sig}: C={c:.5f} D={d:.5f} X={Xv:.5f} Y={Yv:.5f} |G|={abs(gv):.3e} |S|={abs(sv):.3e} {tag}", flush=True)
                if gok and sok: fwd_ok+=1
                else: fwd_bad+=1
                found=True
                break
    if not found: print(f"  sig={sig}: no E_L=E_R=0 solution found", flush=True)
print(f"Forward: {fwd_ok} OK, {fwd_bad} MISMATCH", flush=True)

print("\n=== Backward: G=S=0  =>  E_L=E_R=0 ? ===", flush=True)
bwd_ok=0; bwd_bad=0
for sig in [0.3,0.6,1.0,1.5,2.0,3.0]:
    found=False
    for _ in range(300):
        X0=rng.uniform(1.05,1.95)
        # Y in (X-1, X^2/4)
        ylo=max(1e-3,X0-1); yhi=X0**2/4*0.98
        if yhi<=ylo: continue
        Y0=rng.uniform(ylo,yhi)
        def F(v):
            Xv,Yv=v
            try:
                return [float(Gf(Xv,Yv,sig)), float(Sf(Xv,Yv,sig))]
            except: return [1e10,1e10]
        sol=root(F,[X0,Y0],method='hybr',options={'xtol':1e-14,'maxfev':8000})
        if sol.success and np.linalg.norm(F(sol.x))<1e-9:
            Xv,Yv=sol.x
            if 1<Xv<2 and Xv-1<Yv<Xv**2/4 and (Xv**2-4*Yv)>1e-6:
                disc=Xv**2-4*Yv
                c=(Xv+np.sqrt(disc))/2; d=(Xv-np.sqrt(disc))/2  # C>D
                if 0<c<1 and 0<d<1:
                    el=float(ELf(c,d,sig)); er=float(ERf(c,d,sig))
                    eok = abs(el)<1e-6 and abs(er)<1e-6
                    tag="OK" if eok else "MISMATCH"
                    print(f"  sig={sig}: X={Xv:.5f} Y={Yv:.5f} C={c:.5f} D={d:.5f} |E_L|={abs(el):.3e} |E_R|={abs(er):.3e} {tag}", flush=True)
                    if eok: bwd_ok+=1
                    else: bwd_bad+=1
                    found=True
                    break
    if not found: print(f"  sig={sig}: no G=S=0 solution found", flush=True)
print(f"Backward: {bwd_ok} OK, {bwd_bad} MISMATCH", flush=True)
