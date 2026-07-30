#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Numerical check: does {G=S=det=0} have interior solutions?
G,S verified to contain the non-palindromic stationary curve.
det = det[grad G, grad S, D^2 grad(rho^9)] = rho-criticality condition (definition).
If {G=S=det=0, interior, a_i in (0,1)} is numerically EMPTY -> B.16 likely empty.
Use multi-start scipy root on (X,Y,sigma) with interior domain check."""
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
G=parse_poly('nonpal_G_clean.txt'); S=parse_poly('nonpal_S_clean.txt')
N=parse_poly('nonpal_rho9_num.txt'); D=parse_poly('nonpal_rho9_den.txt')
Gx,Gy,Gs=G.diff(Xs),G.diff(Ys),G.diff(ss)
Sx,Sy,Ss=S.diff(Xs),S.diff(Ys),S.diff(ss)
Nx,Ny,Ns=N.diff(Xs),N.diff(Ys),N.diff(ss)
Dx,Dy,Ds=D.diff(Xs),D.diff(Ys),D.diff(ss)
Cx=D*Nx-N*Dx; Cy=D*Ny-N*Dy; Cz=D*Ns-N*Ds
det = Gx*(Sy*Cz-Cy*Ss) - Gy*(Sx*Cz-Cx*Ss) + Gs*(Sx*Cy-Cx*Sy)
detP=sp.Poly(sp.expand(det.as_expr()),Xs,Ys,ss,domain=sp.ZZ)
print(f"G deg{G.total_degree()} S deg{S.total_degree()} det deg{detP.total_degree()}", flush=True)

def ev(P,x,y,s):
    t=0.0
    for (a,b,c),coef in P.terms():
        t+=float(coef)*(x**a)*(y**b)*(s**c)
    return t
def F(v):
    x,y,s=v
    return [ev(G,x,y,s),ev(S,x,y,s),ev(detP,x,y,s)]

def center_lift_num(C,D,sig):
    den=C+D-1
    a3=1-C+sig*C*C*(1-C)/den
    a6=1-D+sig*D*D*(1-D)/den
    a2=1-a3+sig*a3*a3*(1-a3)/(a3+C-1)
    a7=1-a6+sig*a6*a6*(1-a6)/(a6+D-1)
    return [1,a2,a3,C,D,a6,a7,1]

rng=np.random.default_rng(2024)
sols=[]
ntry=0
for _ in range(60):
    X0=rng.uniform(1.05,1.95)
    ylo=X0-1+1e-3; yhi=X0**2/4-1e-3
    if yhi<=ylo: continue
    Y0=rng.uniform(ylo,yhi); s0=rng.uniform(0.1,8.0)
    ntry+=1
    sol=root(F,[X0,Y0,s0],method='hybr',options={'xtol':1e-15,'maxfev':20000})
    if sol.success and np.linalg.norm(F(sol.x))<1e-8:
        Xv,Yv,sv=sol.x
        if 1<Xv<2 and Xv-1<Yv<Xv**2/4-1e-8 and sv>0:
            disc=Xv**2-4*Yv
            if disc>1e-8:
                c=(Xv+np.sqrt(disc))/2; d=(Xv-np.sqrt(disc))/2
                if 0<c<1 and 0<d<1 and c+d>1:
                    try:
                        a=center_lift_num(c,d,sv)
                    except: continue
                    if all(0<ai<1 for ai in a) and all(a[i]+a[i+1]>1 for i in range(7)):
                        sols.append((Xv,Yv,sv,c,d))
                        print(f"  INTERIOR SOL: X={Xv:.6f} Y={Yv:.6f} sig={sv:.6f} C={c:.5f} D={d:.5f}", flush=True)
print(f"\nTried {ntry} starts. Interior solutions found: {len(sols)}", flush=True)
if not sols:
    print("==> {G=S=det=0} NUMERICALLY EMPTY in interior -> B.16 strongly supported empty.", flush=True)
