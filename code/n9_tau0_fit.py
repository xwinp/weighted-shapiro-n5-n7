#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SAFE numerical Puiseux slope for tau->0 corner branch (v2).
Corner branch has eta=1-xi ~ 4 tau^2 (tiny), v ~ tau^ell.  So we must select the
solution with SMALLEST eta (closest to corner), not smallest v.  Starts biased
hard toward xi->1.  mpmath only, bounded memory.
"""
import random, mpmath as mp, sympy as sp, math
mp.mp.dps=16
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
random.seed(11)
# corner branch eta~4tau^2: at tau=0.01 eta~4e-4, tau=0.1 eta~4e-2.  stay in findable range.
taus=[0.01,0.015,0.02,0.03,0.05,0.07,0.1,0.13,0.16,0.2]
rows=[]
for t0 in taus:
    found=[]
    for _ in range(60):
        # bias hard to corner: xi in [0.999, 0.9999999], v tiny
        x0=1-random.uniform(1e-6,5e-3); v0=random.uniform(1e-7,2e-2)
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,mp.mpf(t0)),fGp(vv,xx,mp.mpf(t0))),(mp.mpf(v0),mp.mpf(x0)),tol=1e-18,maxsteps=45)
            v1,x1=float(sol[0]),float(sol[1])
            if not(1e-9<v1<1 and 0.9<x1<1): continue
            r=max(abs(float(fGm(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))),abs(float(fGp(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))))
            if r>1e-11: continue
            a=fullvec(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0))
            if not all(0<float(a[i])<1 for i in range(1,7)): continue
            if not all(float(a[i]+a[i+1])>1 for i in range(7)): continue
            pal=abs(float((cc.subs({v:mp.mpf(v1),xi:mp.mpf(x1)})-dd.subs({v:mp.mpf(v1),xi:mp.mpf(x1)}))))<1e-3
            if pal: continue
            found.append((v1,x1,1-x1))
        except Exception: pass
    if found:
        # corner branch = smallest eta
        found.sort(key=lambda q:q[2])
        b=found[0]
        rows.append((t0,b[0],b[1],b[2]))
        print(f"tau={t0:.4f}: n={len(found)} | smallest-eta: v={b[0]:.4e} xi={b[1]:.7f} eta={b[2]:.4e} eta/4tau^2={b[2]/(4*t0*t0):.3f}",flush=True)
    else:
        print(f"tau={t0:.4f}: none",flush=True)
print("---",flush=True)
if len(rows)>=3:
    pts=[(math.log(r[0]),math.log(r[1])) for r in rows]
    n=len(pts); sx=sum(p[0] for p in pts); sy=sum(p[1] for p in pts)
    sxx=sum(p[0]**2 for p in pts); sxy=sum(p[0]*p[1] for p in pts)
    ell=(n*sxy-sx*sy)/(n*sxx-sx*sx)
    print("fitted ell (v~tau^ell) over %d pts: %.4f"%(n,ell),flush=True)
    peta=[(math.log(r[0]),math.log(r[3])) for r in rows if r[3]>0]
    n=len(peta); sx=sum(p[0] for p in peta); sy=sum(p[1] for p in peta)
    sxx=sum(p[0]**2 for p in peta); sxy=sum(p[0]*p[1] for p in peta)
    m=(n*sxy-sx*sy)/(n*sxx-sx*sx)
    print("fitted m  (eta~tau^m)  over %d pts: %.4f  (expect ~2)"%(n,m),flush=True)
print("DONE",flush=True)
