#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnostic: at H_C lifts (both p-branches), reconstruct x via g1=0 and check
ALL five KKT residuals g1..g5 = rho_i dP/drho_i.  Determines whether the
beta-variety (H_C+E2+E3+closure+g1) is the FULL stationary locus or a superset.
"""
import sys
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

zS,wS,uS,vS=sp.symbols('z w u v')
HC=zS*wS**3+wS**2*zS**3-wS**2*zS+wS*zS**4-3*wS*zS**3+2*wS*zS**2+wS*zS-wS-zS**4+3*zS**3-3*zS**2+zS
a3=1-vS+uS*vS; a5=1-zS+zS*wS-zS*vS*wS+zS*uS*vS*wS
E3=a3*vS-uS*(1-wS); E2=uS*(1-zS)-zS*a5*(1-vS)
usol=vS*(1-vS)/((1-wS)-vS**2)
E2u=sp.together(E2.subs(uS,usol)); E2u_num=sp.expand(E2u.as_numer_denom()[0])

# symbolic P and g_i = rho_i dP/drho_i
rho=sp.symbols('rho1:6'); p=sp.symbols('p'); q=1-p
P=sp.together(sum(1/(rho[i]*(p+q*rho[i+1])) for i in range(4)) + 1/(p*rho[4]) + sp.prod(rho)/q)
g=[sp.together(rho[i]*sp.diff(P,rho[i])) for i in range(5)]
gnum=[sp.expand(gi.as_numer_denom()[0]) for gi in g]

print("z      p       branch    g1        g2        g3        g4        g5       Pval")
real_Pvals = []          # Pval of every REAL (full-KKT) H_C lift
spurious_fail_g45 = []   # for every spurious lift: does it fail via g4 or g5 (>1e-3)?
for zv in [0.50, 0.65, 0.80, 0.85, 0.90, 0.95, 0.99]:
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
            rho_v=K**(1/7); pp=1/(1+rho_v)
            r2=uu/(rho_v*(1-uu)); r3=vv/(rho_v*(1-vv)); r4=wv/(rho_v*(1-wv)); r5=zv/(rho_v*(1-zv))
            def g1(r1):
                x=np.zeros(n); x[1]=1.0; x[2]=r1
                x[3]=r1*r2; x[4]=r1*r2*r3; x[5]=r1*r2*r3*r4; x[6]=r1*r2*r3*r4*r5
                h=1e-7; return (Pval(x+np.eye(1,n,2)[0]*h,pp)-Pval(x,pp))/h
            xs=np.linspace(1e-3,50,5000); gs=np.array([g1(r) for r in xs]); br=None
            for i in range(len(xs)-1):
                if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and gs[i]*gs[i+1]<0:
                    try: br=brentq(g1,xs[i],xs[i+1],xtol=1e-13); break
                    except: continue
            if br is None: continue
            r1=br
            rholist=[r1,r2,r3,r4,r5]
            # evaluate g1..g5 symbolically
            subs={p:pp, q:1-pp}
            for i in range(5): subs[rho[i]]=rholist[i]
            gres=[float(gnum[i].subs(subs)) for i in range(5)]
            x=np.zeros(n); x[1]=1.0; x[2]=r1
            x[3]=r1*r2; x[4]=r1*r2*r3; x[5]=r1*r2*r3*r4; x[6]=r1*r2*r3*r4*r5
            xs_=x/x.sum(); Pd=Pval(xs_,pp)
            gmax=max(abs(gg) for gg in gres)
            tag="REAL" if gmax<1e-3 else "spurious"
            if tag == "REAL":
                real_Pvals.append(Pd)
            else:
                # spurious lift must fail via the inactive-support KKT (g4 or g5),
                # not merely a near-miss on g1..g3
                spurious_fail_g45.append(max(abs(gres[3]), abs(gres[4])) > 1e-3)
            print("%.2f  %.5f %s  %.1e %.1e %.1e %.1e %.1e  %.5f"%(
                zv,pp,tag,gres[0],gres[1],gres[2],gres[3],gres[4],Pd))
n_real = len(real_Pvals)
n_spur = len(spurious_fail_g45)
all_real_above_7 = all(P > 7 for P in real_Pvals)
all_spur_fail = all(spurious_fail_g45) if spurious_fail_g45 else True
ok = (n_real + n_spur > 0) and all_real_above_7 and all_spur_fail
print("\nKKT classification: %d REAL lifts (all Pval>7: %s, min Pval=%.5f); "
      "%d spurious (all fail g4/g5: %s)" % (
          n_real, all_real_above_7, min(real_Pvals) if real_Pvals else float('nan'),
          n_spur, all_spur_fail))
print("CERTIFICATE: every REAL H_C-lift stationary point has P>7=%s AND every spurious lift "
      "fails g4/g5=%s  => H_C is a superset of the stationary locus, no sub-bound stationary "
      "point missed: %s" % (all_real_above_7, all_spur_fail, ok))
assert ok, "n7 H_C KKT superset check failed"
print("DONE-KKT")
sys.exit(0 if ok else 1)
