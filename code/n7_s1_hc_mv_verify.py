#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the mean-value P bound is a VALID lower bound (MV_lo <= P_true) and
that direct (definitely-rigorous) eval also gives >7 where MV does."""
import mpmath as mp, numpy as np, sympy as sp
mp.mp.ivprec=110; IV=mp.iv; mp.mp.prec=80
exec(open('code/n7_s1_hc_rigorous_cert.py').read().split('with open')[0])

def Ptrue(u,v,w,z,rho):
    A=(1-u)*(1+rho); B=(1+rho)*u*v*w*z/(rho**5*(1-u)*(1-v)*(1-w)*(1-z))
    C=rho*(1+rho)*((1-v)*(1-u)/u+(1-w)*(1-v)/v+(1-z)*(1-w)/w+(1-z)/z)
    return C+2*np.sqrt(A*B)

print("s        c      v      branch   P_true   MV_lo    direct_lo  valid(MV<=Pt)")
for sc in [0.198,0.30,0.40,0.11,0.05,0.01,1e-3,1e-4]:
    lifts=admissible_lifts(sc)
    for (c0,v0,u0) in lifts:
        zv=1-sc; w0=c0*sc
        a5v=1-zv+zv*w0-zv*v0*w0+zv*u0*v0*w0
        K=u0*v0*w0*zv**3*a5v**2/((1-v0)*(1-w0)*(1-zv)**3)
        rho=K**(1/7); Pt=Ptrue(u0,v0,w0,zv,rho)
        hw=mp.mpf(2)*mp.mpf(10)**(-5)
        S=IV.mpf([mp.mpf(sc)-hw, mp.mpf(sc)+hw])
        C=krawczyk_c(S,c0); V=krawczyk_v(C,S,v0) if C is not None else None
        if V is None: continue
        mv=P_box(V,C,S); dr=P_direct(iv(V),iv(C),iv(S))
        br="real" if 2.0<c0<2.5 else "spur"
        valid = (mv is None) or (mv <= Pt+1e-9)
        print("%.4e %.4f %.4f %s  %.5f  %s  %s  %s"%(
            sc,c0,v0,br,Pt,('%.5f'%mv) if mv else 'None',('%.5f'%dr) if dr else 'None',valid))
print("DONE-MVVERIFY")
