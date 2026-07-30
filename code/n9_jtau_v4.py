#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""J_tau v4 (lean): classify strict-interior vs boundary, sign split.
Question: does J_tau change sign among STRICT-interior (v>=BDY) curve points?
  - if single sign in strict interior -> Case 1 likely (turning pts only at v=0 boundary).
  - if mixed -> Case 2 (genuine turning pts in strict interior).
Lower dps (15) + Jtau via low-order finite difference (fast) instead of mp.diff.
Then targeted 3-eq findroot only between opposite-sign strict pairs.
"""
import random, mpmath as mp, sympy as sp
mp.mp.dps = 18
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
print("lambdified.",flush=True)

# fast Jtau via central difference (order ~1e-6 abs, enough for sign away from 0)
def Jtau_fd(vv,xx,tt,h=1e-5):
    gmv=(fGm(vv+h,xx,tt)-fGm(vv-h,xx,tt))/(2*h)
    gmx=(fGm(vv,xx+h,tt)-fGm(vv,xx-h,tt))/(2*h)
    gpv=(fGp(vv+h,xx,tt)-fGp(vv-h,xx,tt))/(2*h)
    gpx=(fGp(vv,xx+h,tt)-fGp(vv,xx-h,tt))/(2*h)
    return gmv*gpx-gmx*gpv

def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]

def adm(v1,x1,t1):
    a=fullvec(v1,x1,t1)
    if not all(0<float(a[i])<1 for i in range(1,7)): return False,None
    if not all(float(a[i]+a[i+1])>1 for i in range(7)): return False,None
    pal=abs(float((1+v1+(1-v1)*x1)/2)-float((1+v1-(1-v1)*x1)/2))<1e-3
    return True,pal

random.seed(321)
pts=[]
for t0 in [0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.9,0.95,0.97]:
    for _ in range(40):
        v0=random.uniform(0.02,0.95); x0=random.uniform(0.02,0.98)
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,t0),fGp(vv,xx,t0)),
                            (mp.mpf(v0),mp.mpf(x0)),tol=1e-20,maxsteps=50)
            v1,x1=float(sol[0]),float(sol[1])
            if not(0<v1<1 and 0<x1<1): continue
            r=max(abs(float(fGm(v1,x1,t0))),abs(float(fGp(v1,x1,t0))))
            if r>1e-8: continue
            ok,pal=adm(v1,x1,t0)
            if not ok or pal: continue
            Jv=float(Jtau_fd(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))
            pts.append((t0,v1,x1,Jv))
        except Exception: pass

BDY=1e-3
strict=[p for p in pts if p[1]>=BDY]; bdy=[p for p in pts if p[1]<BDY]
print(f"\n{len(pts)} nonpal-interior curve pts: {len(strict)} strict (v>={BDY}), {len(bdy)} boundary (v<{BDY})",flush=True)
if strict:
    ss=[p[3] for p in strict]
    sp_=sum(1 for x in ss if x>0); sn=sum(1 for x in ss if x<0)
    print(f"STRICT-INTERIOR J_tau: pos={sp_} neg={sn} min|J|={min(abs(x) for x in ss):.4e} max|J|={max(abs(x) for x in ss):.4e}",flush=True)
    pmin=min(strict,key=lambda q:abs(q[3]))
    print(f"  strict min|J| at tau={pmin[0]} v={pmin[1]:.5f} xi={pmin[2]:.5f} J={pmin[3]:+.4e}",flush=True)
    # show a few strict points sorted by tau
    for q in sorted(strict,key=lambda q:(q[0],q[1]))[:25]:
        print(f"  tau={q[0]:.3f} v={q[1]:.4f} xi={q[2]:.4f} J={q[3]:+.4e}",flush=True)
if bdy:
    bs=[p[3] for p in bdy]
    print(f"BOUNDARY J_tau: pos={sum(1 for x in bs if x>0)} neg={sum(1 for x in bs if x<0)} min|J|={min(abs(x) for x in bs):.4e}",flush=True)

# targeted 3-eq: only seed between opposite-sign strict pairs that are close in (tau,v,xi)
def Jtau_mp(vv,xx,tt):
    gmv=mp.diff(lambda z:fGm(z,xx,tt),vv)
    gmx=mp.diff(lambda z:fGm(vv,z,tt),xx)
    gpv=mp.diff(lambda z:fGp(z,xx,tt),vv)
    gpx=mp.diff(lambda z:fGp(vv,z,tt),xx)
    return gmv*gpx-gmx*gpv

strict_sorted=sorted(strict,key=lambda q:(q[0],q[1]))
turns=[]
for i in range(len(strict_sorted)-1):
    a=strict_sorted[i]; b=strict_sorted[i+1]
    if a[3]*b[3]<0 and abs(a[0]-b[0])<0.12 and abs(a[1]-b[1])<0.2 and abs(a[2]-b[2])<0.2:
        # sign change between close strict points -> search for Jt=0 between
        vt=(a[1]+b[1])/2; xt=(a[2]+b[2])/2; tt=(a[0]+b[0])/2
        try:
            sol=mp.findroot(lambda vv,xx,ttt:(fGm(vv,xx,ttt),fGp(vv,xx,ttt),Jtau_mp(vv,xx,ttt)),
                            (mp.mpf(vt),mp.mpf(xt),mp.mpf(tt)),tol=1e-18,maxsteps=40)
            v1,x1,t1=float(sol[0]),float(sol[1]),float(sol[2])
            if 0<v1<1 and 0<x1<1 and 0<t1<1:
                r=max(abs(float(fGm(v1,x1,t1))),abs(float(fGp(v1,x1,t1))),abs(float(Jtau_mp(mp.mpf(v1),mp.mpf(x1),mp.mpf(t1)))))
                if r<1e-6:
                    ok,pal=adm(v1,x1,t1)
                    if ok and v1>=BDY:
                        turns.append((t1,v1,x1,pal))
                        print(f"  TURN(strict): tau={t1:.5f} v={v1:.5f} xi={x1:.5f} {'PAL' if pal else 'NONPAL'}",flush=True)
        except Exception: pass

print(f"\n{len(turns)} strict-interior turning points (Jt=0, non-boundary).",flush=True)
print("DONE",flush=True)
