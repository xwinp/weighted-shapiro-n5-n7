#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the clean algebraic form  P = C + 2*sqrt(A*B)  on the H_C branch
(no rho1 variable needed). Derivation:
  beta=(u,v,w,z)=(b1,b2,b3,b4);  rho_{i+1}=rho*beta_i/(1-beta_i) (i=1..4) -> rho2..rho5.
  p=1/(1+rho), q=rho/(1+rho).
  P = sum_{i=1}^4 1/(rho_i (p+q rho_{i+1})) + 1/(p rho5) + prod_rho/q.
  rho1-terms: term1=(1-b1)/(p rho1)=A/rho1,  prod_rho/q=rho1*(rho2 rho3 rho4 rho5)/q = rho1*B.
  => P = A/rho1 + B*rho1 + C,  g1=rho1 dP/drho1=-A/rho1+B rho1=0 => rho1=sqrt(A/B) (>0).
  => P = C + 2 sqrt(A B).
  A=(1-u)/p=(1-u)(1+rho).
  B=(rho2 rho3 rho4 rho5)/q = rho^3(1+rho) u v w z / ((1-u)(1-v)(1-w)(1-z)).
  => A B = rho^3 (1+rho)^2  u v w z / ((1-v)(1-w)(1-z))   [ (1-u) cancels ].
  C=(1+rho)/rho * [ (1-v)(1-u)/u + (1-w)(1-v)/v + (1-z)(1-w)/w + (1-z)/z ].
Crossing P=7  <=>  (C-7)^2 = 4 A B  (polynomial, no sqrt, no rho1).

This script NUMERICALLY checks P=C+2sqrt(AB) against direct Pval at:
  (i) H_B points (symmetric branch; AB should be a perfect square there -> sqrt-free),
  (ii) the 8 H_C admissible points from the definitive trace.
Then reports degrees of the crossing poly (C-7)^2-4AB in (w,z,rho) after substituting
the rational u,v lifts on H_C.
"""
import numpy as np
from scipy.optimize import root, brentq
import sympy as sp

n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s

# ---- direct P from beta+rho via the closed form  P = C + 2 sqrt(A B) ----
# rho_{i+1} = beta_i / (rho (1-beta_i))   (since p/q = 1/rho)
def P_closed(u,v,w,z,rho):
    A=(1-u)*(1+rho)                                          # = (1-b1)/p
    B=(1+rho)*u*v*w*z/(rho**5*(1-u)*(1-v)*(1-w)*(1-z))       # = (rho2 rho3 rho4 rho5)/q
    C=rho*(1+rho)*((1-v)*(1-u)/u + (1-w)*(1-v)/v + (1-z)*(1-w)/w + (1-z)/z)
    return C + 2*np.sqrt(A*B), A, B, C

# ---- solve S1 stationary at given p ----
sup_free=[2,3,4,5,6]
def solve_S1(p, init):
    def grad(v):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(sup_free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(sup_free)):
            xp=x.copy(); xp[sup_free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14,'maxfev':40000})
    if r.success and max(abs(r.fun))<1e-7 and all(ri>1e-8 for ri in r.x):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(sup_free): x[idx]=r.x[j]
        return x/x.sum()
    return None

print("=== (i) H_B points: P_closed vs Pval ===")
init=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
prev=init
for pp in [0.20,0.27,0.329,0.40,0.60,0.90]:
    x=solve_S1(pp, prev)
    if x is None:
        for ii in [[0.3,0.7,1.5,0.07,1.3],[1,1,1,1,1]]:
            x=solve_S1(pp,ii)
            if x is not None: break
    if x is None: print("  p=%.3f NO CONVERGE"%pp); continue
    prev=[x[2],x[3],x[4],x[5],x[6]]
    b,c,d,e,f=x[2],x[3],x[4],x[5],x[6]
    rl=[b,c/b,d/c,e/d,f/e]; q=1-pp; rho=q/pp
    u=q*rl[1]/(pp+q*rl[1]); v=q*rl[2]/(pp+q*rl[2]); w=q*rl[3]/(pp+q*rl[3]); z=q*rl[4]/(pp+q*rl[4])
    # note: beta1=u corresponds to rho2; but our (u,v,w,z) here are beta1..beta4
    Pc,A,B,C=P_closed(u,v,w,z,rho)
    print("  p=%.3f  Pval=%.7f  P_closed=%.7f  diff=%+.2e  sqrt(AB)=%+.5f (H_B: near-integer? AB=%+.5f)"%(
        pp, Pval(x,pp), Pc, Pc-Pval(x,pp), np.sqrt(A*B), A*B))

print("\n=== (ii) H_C points: P_closed vs Pval (reconstruct via g1) ===")
# replicate H_C trace
zS,wS,uS,vS=sp.symbols('z w u v')
HC=zS*wS**3+wS**2*zS**3-wS**2*zS+wS*zS**4-3*wS*zS**3+2*wS*zS**2+wS*zS-wS-zS**4+3*zS**3-3*zS**2+zS
a3=1-vS+uS*vS; a5=1-zS+zS*wS-zS*vS*wS+zS*uS*vS*wS
E3=a3*vS-uS*(1-wS); E2=uS*(1-zS)-zS*a5*(1-vS)
usol=vS*(1-vS)/((1-wS)-vS**2)
E2u=sp.together(E2.subs(uS,usol)); E2u_num=sp.expand(E2u.as_numer_denom()[0])
print("  H_C trace:")
hc_pts=[]
for zv in [0.80,0.82,0.85,0.88,0.90,0.92,0.95,0.97,0.99]:
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
            pp=1/(1+K**(1/7)); rho=K**(1/7)
            # reconstruct x via g1=0 (rho1) to compute direct Pval
            r2=uu/(rho*(1-uu)); r3=vv/(rho*(1-vv)); r4=wv/(rho*(1-wv)); r5=zv/(rho*(1-zv))
            def g1(r1):
                x=np.zeros(n); x[1]=1.0; x[2]=r1
                x[3]=r1*r2; x[4]=r1*r2*r3; x[5]=r1*r2*r3*r4; x[6]=r1*r2*r3*r4*r5
                h=1e-7; return (Pval(x+np.eye(1,n,2)[0]*h,pp)-Pval(x,pp))/h
            xs=np.linspace(1e-3,30,3000); gs=np.array([g1(r) for r in xs]); br=None
            for i in range(len(xs)-1):
                if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and gs[i]*gs[i+1]<0:
                    try: br=brentq(g1,xs[i],xs[i+1],xtol=1e-13); break
                    except: continue
            if br is None: continue
            r1=br
            x=np.zeros(n); x[1]=1.0; x[2]=r1
            x[3]=r1*r2; x[4]=r1*r2*r3; x[5]=r1*r2*r3*r4; x[6]=r1*r2*r3*r4*r5
            xs_=x/x.sum(); Pd=Pval(xs_,pp)
            # closed form uses beta=(u,v,w,z)=(beta1..beta4) BUT note: in trace, u,v,w,z are beta1,beta2,beta3,beta4?
            # verify mapping: r2=rho*u/(1-u) => u=beta1. yes.
            Pc,A,B,C=P_closed(uu,vv,wv,zv,rho)
            print("  z=%.2f p=%.5f Pval=%.6f P_closed=%.6f diff=%+.2e"%(zv,pp,Pd,Pc,Pc-Pd))
            hc_pts.append((zv,pp,Pd,Pc))
print("\n  H_C points found: %d, all P_closed match Pval? %s"%(
    len(hc_pts), all(abs(p[2]-p[3])<1e-4 for p in hc_pts)))
print("DONE-CHECK")
