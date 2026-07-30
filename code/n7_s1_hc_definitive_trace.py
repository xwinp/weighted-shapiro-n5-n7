#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Definitive H_C-branch trace. For z in a grid over (zeta,1):
  (1) solve H_C(w,z)=0 (cubic in w) for real roots;
  (2) for each w in (0,1), solve E2,E3 (the two remaining p-free KKT reductions)
      for (u,v) in (0,1)^2 -- these come from the beta-recurrence;
  (3) closure K=(q/p)^7 = u v w z^3 a5^2 / ((1-v)(1-w)(1-z)^3)  -> p=1/(1+K^{1/7});
  (4) reconstruct x (rho1 from 1D g1=0) and check FULL 5-var KKT residual;
  (5) if KKT residual small -> real stationary point: report P, Morse (correct Hessian).

H_C(w,z) = z w^3 + z^3 w^2 - z w^2 + z^4 w - 3z^3 w + 2z^2 w + z w - w - z^4 + 3z^3 - 3z^2 + z
E3: a3 v - u(1-w) = 0,  a3 = 1 - v + u v
E2: u(1-z) - z a5 (1-v) = 0,  a5 = 1 - z + z w - z v w + z u v w
Solve E3 for u = a3 v/(1-w) = (1-v+uv)v/(1-w)  -> u(1-w)=v-v^2+uv^2 -> u((1-w)-v^2)=v-v^2
   u = v(1-v)/((1-w)-v^2).  Sub into E2.
"""
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
def grad_all(x,p):
    free=[2,3,4,5,6]; g=[]; h=1e-7; f0=Pval(x,p)
    for idx in free:
        xp=x.copy(); xp[idx]+=h; g.append((Pval(xp,p)-f0)/h)
    return g
def reduced_eigs(x,p,h):
    free=[2,3,4,5,6]; k=5
    def Pv(v):
        xx=np.zeros(n); xx[1]=1.0
        for j,idx in enumerate(free): xx[idx]=v[j]
        return Pval(xx,p)
    v0=np.array([x[idx] for idx in free]); H=np.zeros((k,k))
    for i in range(k):
        ei=np.zeros(k); ei[i]=h; H[i,i]=(Pv(v0+ei)-2*Pv(v0)+Pv(v0-ei))/h**2
    for i in range(k):
        for j in range(i+1,k):
            ei=np.zeros(k); ej=np.zeros(k); ei[i]=h; ej[j]=h
            H[i,j]=H[j,i]=(Pv(v0+ei+ej)-Pv(v0+ei-ej)-Pv(v0-ei+ej)+Pv(v0-ei-ej))/(4*h**2)
    return np.linalg.eigvalsh(H)

z,w,u,v=sp.symbols('z w u v')
HC=z*w**3+w**2*z**3-w**2*z+w*z**4-3*w*z**3+2*w*z**2+w*z-w-z**4+3*z**3-3*z**2+z
a3=1-v+u*v
a5=1-z+z*w-z*v*w+z*u*v*w
E3=a3*v-u*(1-w)
E2=u*(1-z)-z*a5*(1-v)
# u = v(1-v)/((1-w)-v^2)  from E3=0
usol = v*(1-v)/((1-w)-v**2)
E2u = sp.simplify(E2.subs(u,usol))
E2u_num = sp.together(E2u).as_numer_denom()[0]
E2u_num = sp.expand(E2u_num)
print("E2(u(v)) numerator in v (given z,w):", "deg_v =", sp.Poly(E2u_num,v).degree() if v in E2u_num.free_symbols else 0)

real_pts=[]
zs_grid=[0.80,0.82,0.85,0.88,0.90,0.92,0.95,0.97,0.99]
print("\nz      w-root   u        v        p        P        P-7     Morse  KKTmax")
for zv in zs_grid:
    HCz=sp.Poly(HC.subs(z,zv),w)
    wroots=[float(sp.re(r)) for r in sp.nroots(HCz,n=20) if abs(sp.im(r))<1e-9]
    for wv in wroots:
        if not (0<wv<1): continue
        # E2u_num as poly in v
        Ev=sp.Poly(E2u_num.subs({z:zv,w:wv}), v)
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
            pp=1/(1+K**(1/7))
            # reconstruct x
            q=1-pp
            rhos=[None,(pp/q)*uu/(1-uu),(pp/q)*vv/(1-vv),(pp/q)*wv/(1-wv),(pp/q)*zv/(1-zv)]
            def g1(r1):
                x=np.zeros(n); x[1]=1.0; x[2]=r1; c=r1
                for i in range(1,5): c=c*rhos[i]; x[2+i]=c
                h=1e-7; return (Pval(x+np.eye(1,n,2)[0]*h,pp)-Pval(x,pp))/h
            xs=np.linspace(1e-3,30,3000); gs=np.array([g1(r) for r in xs])
            br=None
            for i in range(len(xs)-1):
                if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and gs[i]*gs[i+1]<0:
                    try: br=brentq(g1,xs[i],xs[i+1],xtol=1e-13); break
                    except: continue
            if br is None: continue
            r1=br
            x=np.zeros(n); x[1]=1.0; x[2]=r1; c=r1
            for i in range(1,5): c=c*rhos[i]; x[2+i]=c
            xs_=x/x.sum(); Pv=Pval(xs_,pp)
            kkt=max(abs(g) for g in grad_all(x,pp))
            if kkt>0.05: continue  # not a real KKT point
            eigs=reduced_eigs(xs_,pp,1e-5); neg=int((eigs<-1e-3).sum())
            real_pts.append((zv,wv,uu,vv,pp,Pv,neg,kkt))
            print("%.2f  %.5f  %.5f  %.5f  %.6f  %.6f %+.6f  %s  %.1e"%(
                zv,wv,uu,vv,pp,Pv,Pv-7,"MIN" if neg==0 else "SADDLE%d"%neg,kkt))
print("\nTotal real H_C stationary points found: %d"%len(real_pts))
# summary: any with P<7?
bad=[p for p in real_pts if p[5]<7-1e-6]
print("H_C points with P<7:", len(bad))
mins=[p for p in real_pts if p[6]==0]
print("H_C local minima: %d, all P>7? %s"%(len(mins), all(p[5]>7 for p in mins)))
print("min P among H_C mins:", min((p[5] for p in mins),default=None))
print("DONE")
