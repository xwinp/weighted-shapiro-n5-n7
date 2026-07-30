#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""J_tau global-parameter check (v2): eval J_tau via mpmath high-precision diff
(symbolic lambdify hit RecursionError). On {G^-=G^+=0} curve points, J_tau =
Gm_v*Gp_xi - Gm_xi*Gp_v. If nonzero (bounded from 0) on the Omega_np curve,
tau is a global parameter -> branches are analytic graphs (v(tau),xi(tau)).
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
cc=(1+v+(1-v)*xi)/2; dd=(1+v+(1-v)*xi)/2  # placeholder
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,s:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,s:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])
fGm=sp.lambdify((v,xi,tau),Gmt,'mpmath'); fGp=sp.lambdify((v,xi,tau),Gpt,'mpmath')
print("lambdified G^-,G^+.",flush=True)

def Jtau(vv,xx,tt):
    gmv=mp.diff(lambda z:fGm(z,xx,tt),vv)
    gmx=mp.diff(lambda z:fGm(vv,z,tt),xx)
    gpv=mp.diff(lambda z:fGp(z,xx,tt),vv)
    gpx=mp.diff(lambda z:fGp(vv,z,tt),xx)
    return gmv*gpx-gmx*gpv

def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]

random.seed(321)
Jvals=[]; near=[]
for t0 in [0.25,0.5,0.75,0.9,0.95]:
    for _ in range(30):
        v0=random.uniform(0.05,0.95); x0=random.uniform(0.02,0.98)
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,t0),fGp(vv,xx,t0)),
                            (mp.mpf(v0),mp.mpf(x0)),tol=1e-22,maxsteps=50)
            v1,x1=float(sol[0]),float(sol[1])
            if not(0<v1<1 and 0<x1<1): continue
            r=max(abs(float(fGm(v1,x1,t0))),abs(float(fGp(v1,x1,t0))))
            if r>1e-8: continue
            a=fullvec(v1,x1,t0)
            if not all(0<float(a[i])<1 for i in range(1,7)): continue
            if not all(float(a[i]+a[i+1])>1 for i in range(7)): continue
            pal=abs(float((1+v1+(1-v1)*x1)/2)-float((1+v1-(1-v1)*x1)/2))<1e-3
            if pal: continue
            Jv=float(Jtau(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))
            Jvals.append((t0,v1,x1,Jv))
            if abs(Jv)<1e-2: near.append((t0,v1,x1,Jv))
        except Exception: pass
    print(f"  tau={t0}: {sum(1 for p in Jvals if p[0]==t0)} pts",flush=True)

print(f"\n{len(Jvals)} nonpal-interior curve pts.",flush=True)
if Jvals:
    js=[p[3] for p in Jvals]
    pos=sum(1 for x in js if x>0); neg=sum(1 for x in js if x<0)
    print(f"J_tau signs: pos={pos} neg={neg}  min|J|={min(abs(x) for x in js):.4e}  max|J|={max(abs(x) for x in js):.4e}",flush=True)
    print(f"near-zero |J|<1e-2: {len(near)}",flush=True)
    for z in sorted(near,key=lambda q:abs(q[3]))[:8]:
        print(f"  tau={z[0]} v={z[1]:.4f} xi={z[2]:.4f} J={z[3]:+.4e}",flush=True)
print("DONE",flush=True)
