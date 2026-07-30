#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stepwise resultant elimination for S1 KKT  (Issue 7, step 2: after u-elim)."""
import sympy as sp
r1,u,v,w,z,rho = sp.symbols('r1 u v w z rho', positive=True)
A=(1-u)*(1+rho)
B=(1+rho)*u*v*w*z/(rho**5*(1-u)*(1-v)*(1-w)*(1-z))
S=(1-v)*(1-u)/u+(1-w)*(1-v)/v+(1-z)*(1-w)/w+(1-z)/z
C=rho*(1+rho)*S
dA={s:sp.together(sp.diff(A,s)) for s in (u,v,w,z)}
dB={s:sp.together(sp.diff(B,s)) for s in (u,v,w,z)}
dC={s:sp.together(sp.diff(C,s)) for s in (u,v,w,z)}
def L(j):
    expr=B*r1*dC[j]+B*dA[j]+A*dB[j]
    return sp.expand(sp.together(expr).as_numer_denom()[0])
Lu,Lv,Lw,Lz=L(u),L(v),L(w),L(z)
def ab(Lx):
    c1=sp.expand(Lx.coeff(r1,1)); c0=sp.expand(Lx-c1*r1); return c1,c0
au,bu=ab(Lu); av,bv=ab(Lv); aw,bw=ab(Lw); az,bz=ab(Lz)
def elim(ai,bi,aj,bj): return sp.expand(sp.together(ai*bj-aj*bi).as_numer_denom()[0])
E_uw=elim(au,bu,aw,bw); E_uz=elim(au,bu,az,bz); E_vw=elim(av,bv,aw,bw)
E_vz=elim(av,bv,az,bz); E_wz=elim(aw,bw,az,bz); E_g1=sp.expand(sp.together(B*bu**2-A*au**2).as_numer_denom()[0])

# E3 = u v^2 + u w - u - v^2 + v = 0  ->  u = v(1-v)/((1-w)-v^2)   [usol]
usol = sp.together(v*(1-v)/((1-w)-v**2))
def subs_u(E): return sp.expand(sp.together(E.subs(u, usol)).as_numer_denom()[0])
print("substituting u = v(1-v)/((1-w)-v^2)  [from E3=0]", flush=True)
for nm,E in [('E_uw',E_uw),('E_uz',E_uz),('E_vw',E_vw),('E_vz',E_vz),('E_wz',E_wz),('E_g1',E_g1)]:
    Es=subs_u(E); print("\n%s -> factored:"%nm, flush=True); print(sp.factor(Es), flush=True)
print("\nDONE-STEP2", flush=True)
