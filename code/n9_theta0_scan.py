#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Decisive B.16 test: scan {G^-=G^+=Theta=0} in (v,xi,tau) (3 eqns, 3 unknowns).
The real B.16 target is emptiness of this set in Omega_np (non-palindromic local
min => Theta=0 AND J=0; so empty {G^-=G^+=Theta=0} => no non-palindromic local min).
The saddle found by the G^-=G^+=0 scan had Theta=-0.0023 != 0, so it does NOT
appear here. If this scan finds zero strict-interior non-palindromic solutions,
B.16 stands (and the 'stronger target' is simply false -- saddles exist).
"""
import random, mpmath as mp, sympy as sp
mp.mp.dps = 25
c,d,s=sp.symbols('c d sigma'); v,xi,tau=sp.symbols('v xi tau')
vv=c+d-1
Lc=vv+s*c**2; Uc=vv-s*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+s*d**2; Ud=vv-s*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - s*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - s*Bd*(d*vv**2-Ud*Bd)
Pc=sp.Poly(sp.expand(Fc),c,d,s,domain=sp.ZZ); Pd=sp.Poly(sp.expand(Fd),c,d,s,domain=sp.ZZ)
FL=Pc.exquo(sp.Poly((c-1)*Lc,c,d,s,domain=sp.ZZ)); FR=Pd.exquo(sp.Poly((d-1)*Ld,c,d,s,domain=sp.ZZ))
Gm=sp.Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,s,domain=sp.ZZ).exquo(sp.Poly(c-d,c,d,s,domain=sp.ZZ))
Gp=sp.Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,s,domain=sp.ZZ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,s:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,s:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])

# Theta in (c,d,sigma) then substitute
def theta_expr():
    a1,a2,a3,a4,a5,a6,a7,a8=sp.symbols('a1 a2 a3 a4 a5 a6 a7 a8')
    # rebuild via center_lift symbolic in c,d,s
    g=c+d-1
    A3=1-c+s*c**2*(1-c)/g; A6=1-d+s*d**2*(1-d)/g
    A2=1-A3+s*A3**2*(1-A3)/(A3+c-1); A7=1-A6+s*A6**2*(1-A6)/(A6+d-1)
    A1=sp.Integer(1); A8=sp.Integer(1); A4=c; A5=d
    h2=A2+A3-1; h3=A3+A4-1; h4=A4+A5-1; h5=A5+A6-1; h6=A6+A7-1
    W23=((A6*A7/h6)*(s+h2*h6/(s*A2*A3*A6*A7)-h2*h3*h4*h5*h6/(s**4*A2*A3**2*A4**2*A5**2*A6**2*A7))*(1+A2*A3*h6/(h2*A6*A7)))
    W24=(h2*h3/(s**2*A2*A3**2*A4)*(h4/(s*A4*A5)-1)*(1+A2*A3**2*A4*h5*h6/(h2*h3*A5*A6**2*A7)))
    W34=((A5*A6/h5)*(s+h3*h5/(s*A3*A4*A5*A6)-h3*h4*h5/(s**2*A3*A4**2*A5**2*A6))*(1+A3*A4*h5/(h3*A5*A6)))
    return W23*W24+W23*W34+W24*W34
print("building Theta in (v,xi,tau)...",flush=True)
Thetasub=sp.together(theta_expr().subs({c:cc,d:dd,s:sig}))
Thn=sp.expand(sp.fraction(Thetasub)[0])
fGm=sp.lambdify((v,xi,tau),Gmt,'mpmath'); fGp=sp.lambdify((v,xi,tau),Gpt,'mpmath')
fTh=sp.lambdify((v,xi,tau),Thn,'mpmath')
print("  lambdified.",flush=True)

def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]

random.seed(7)
n=0; valid=0; th0=0
import itertools
# 3-eqns 3-unknowns: solve (Gm,Gp,Theta) from many starts
for t0 in [0.2,0.5,0.8,0.95]:
    for _ in range(120):
        v0=random.uniform(0.05,0.95); x0=random.uniform(0.02,0.98)
        try:
            sol=mp.findroot(lambda vv,xx,tt:(fGm(vv,xx,tt),fGp(vv,xx,tt),fTh(vv,xx,tt)),
                            (mp.mpf(v0),mp.mpf(x0),mp.mpf(t0)),tol=1e-25,maxsteps=80)
            v1,x1,t1=float(sol[0]),float(sol[1]),float(sol[2])
            if not(0<v1<1 and 0<x1<1 and 0<t1<1): continue
            r=max(abs(float(fGm(v1,x1,t1))),abs(float(fGp(v1,x1,t1))),abs(float(fTh(v1,x1,t1))))
            if r>1e-7: continue
            a=fullvec(v1,x1,t1)
            ok=all(0<float(a[i])<1 for i in range(1,7)) and all(float(a[i]+a[i+1])>1 for i in range(7))
            n+=1
            if ok:
                valid+=1
                pal=abs(float((1+v1+(1-v1)*x1)/2)-float((1+v1-(1-v1)*x1)/2))<1e-4
                tag="PAL" if pal else "NONPAL"
                print(f"  FOUND tau={t1:.4f} v={v1:.4f} xi={x1:.4f} {tag} a={[round(float(z),4) for z in a]}",flush=True)
                if not pal: th0+=1
        except Exception:
            pass
    print(f"  [tau-seed {t0}: cumulative {n} sols, {valid} valid, {th0} nonpal]",flush=True)
print(f"\nTOTAL: {n} G^-=G^+=Theta=0 sols, {valid} strict-interior, {th0} non-palindromic",flush=True)
print("DONE",flush=True)
