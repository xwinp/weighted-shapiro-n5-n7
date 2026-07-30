#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check whether a real positive-v non-palindromic curve branch approaches the
corner (xi->1, v->0) at tau near 1/2.  The refined blow-up suggests NO real
branch at tau=1/2 (longitudinal forces v<0). Numerically solve {G^-=G^+=0} in
(v,xi) for tau in {0.45,0.48,0.49,0.495,0.499,0.4999,0.5,0.5001,0.501,0.505,0.51,0.52,0.55},
small-v & xi->1 starts; report smallest-v solution per tau.
"""
import random, mpmath as mp, sympy as sp
mp.mp.dps=18
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
fGm=sp.lambdify((v,xi,tau),Gmt,'mpmath'); fGp=sp.lambdify((v,xi,tau),Gpt,'mpmath')
def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]
random.seed(31)
taus=[0.45,0.48,0.49,0.495,0.499,0.4999,0.5,0.5001,0.501,0.505,0.51,0.52,0.55]
for t0 in taus:
    found=[]
    for _ in range(60):
        v0=random.uniform(1e-5,0.4); x0=random.uniform(0.6,0.999)
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,mp.mpf(t0)),fGp(vv,xx,mp.mpf(t0))),(mp.mpf(v0),mp.mpf(x0)),tol=1e-18,maxsteps=50)
            v1,x1=float(sol[0]),float(sol[1])
            if not(1e-7<v1<1 and 0<x1<1): continue
            r=max(abs(float(fGm(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))),abs(float(fGp(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))))
            if r>1e-9: continue
            a=fullvec(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0))
            if not all(0<float(a[i])<1 for i in range(1,7)): continue
            if not all(float(a[i]+a[i+1])>1 for i in range(7)): continue
            pal=abs(float((cc.subs({v:mp.mpf(v1),xi:mp.mpf(x1)})-dd.subs({v:mp.mpf(v1),xi:mp.mpf(x1)}))))<1e-3
            if pal: continue
            found.append((v1,x1))
        except Exception: pass
    if found:
        found.sort(key=lambda q:q[0])
        b=found[0]
        # also report the largest-xi (most corner-like) solution
        bx=max(found,key=lambda q:q[1])
        print(f"tau={t0:.4f}: n={len(found)} | smallest v: v={b[0]:.3e} xi={b[1]:.5f} eta={1-b[1]:.2e} | largest xi: v={bx[0]:.3e} xi={bx[1]:.5f}",flush=True)
    else:
        print(f"tau={t0:.4f}: no nonpal-interior solution",flush=True)
print("DONE",flush=True)
