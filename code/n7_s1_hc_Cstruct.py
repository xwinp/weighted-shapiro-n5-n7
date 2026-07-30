#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Probe structure of P = C + 2*sqrt(AB) at admissible H_C lifts.
Is C>7 always (trivializing the cert)? How small is 2sqrt(AB)?"""
import numpy as np
from scipy.optimize import brentq
import sympy as sp
n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s
def P_closed(u,v,w,z,rho):
    A=(1-u)*(1+rho)
    B=(1+rho)*u*v*w*z/(rho**5*(1-u)*(1-v)*(1-w)*(1-z))
    C=rho*(1+rho)*((1-v)*(1-u)/u + (1-w)*(1-v)/v + (1-z)*(1-w)/w + (1-z)/z)
    return C + 2*np.sqrt(A*B), A, B, C
zS,wS,uS,vS=sp.symbols('z w u v')
HC=zS*wS**3+wS**2*zS**3-wS**2*zS+wS*zS**4-3*wS*zS**3+2*wS*zS**2+wS*zS-wS-zS**4+3*zS**3-3*zS**2+zS
a3=1-vS+uS*vS; a5=1-zS+zS*wS-zS*vS*wS+zS*uS*vS*wS
E3=a3*vS-uS*(1-wS); E2=uS*(1-zS)-zS*a5*(1-vS)
usol=vS*(1-vS)/((1-wS)-vS**2)
E2u=sp.together(E2.subs(uS,usol)); E2u_num=sp.expand(E2u.as_numer_denom()[0])
print("z      p       tag      C        2sqrt(AB)  P      C-7")
for zv in [0.30,0.40,0.50,0.60,0.65,0.70,0.78,0.80,0.82,0.85,0.88,0.90,0.92,0.95,0.97,0.99]:
    HCz=sp.Poly(HC.subs(zS,zv),wS)
    wroots=[float(sp.re(r)) for r in sp.nroots(HCz,n=20) if abs(sp.im(r))<1e-9]
    for wv in wroots:
        if not (0<wv<1): continue
        Ev=sp.Poly(E2u_num.subs({zS:zv,wS:wv}), vS)
        vroots=[float(sp.re(r)) for r in sp.nroots(Ev,n=20) if abs(sp.im(r))<1e-9]
        for vv in vroots:
            if not (0<vv<1): continue
            denom=(1-wv)-vv**2
            if abs(denom)<1e-9: continue
            uu=vv*(1-vv)/denom
            if not (0<uu<1): continue
            a5v=1-zv+zv*wv-zv*vv*wv+zv*uu*vv*wv
            if a5v<=0: continue
            K=uu*vv*wv*(zv**3)*a5v**2/((1-vv)*(1-wv)*(1-zv)**3)
            if K<=0: continue
            rho=K**(1/7); pp=1/(1+rho)
            Pc,A,B,C=P_closed(uu,vv,wv,zv,rho)
            print("%.2f  %.5f %s  %.5f  %.5f  %.5f  %+.4f"%(
                zv,pp, "  ",C,2*np.sqrt(A*B),Pc,C-7))
print("DONE-CSTRUCT")
