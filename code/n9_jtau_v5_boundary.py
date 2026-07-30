#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""J_tau v5: (a) cross-check J_tau>0 with high-prec mp.diff at strict-interior pts;
(b) find boundary (v->0) non-palindromic curve points (seed small v0), eval Theta.
Theta routine reused verbatim from n9_theta_on_curve.py (verified, 111 pts Theta<0).
"""
import random, mpmath as mp, sympy as sp
mp.mp.dps = 20
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

def Jtau_mp(vv,xx,tt):
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

def Th_direct(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    a4=c1; a5=d1
    h2=a2+a3-1; h3=a3+a4-1; h4=a4+a5-1; h5=a5+a6-1; h6=a6+a7-1
    try:
        W23=((a6*a7/h6)*(s1+h2*h6/(s1*a2*a3*a6*a7)-h2*h3*h4*h5*h6/(s1**4*a2*a3**2*a4**2*a5**2*a6**2*a7))*(1+a2*a3*h6/(h2*a6*a7)))
        W24=(h2*h3/(s1**2*a2*a3**2*a4)*(h4/(s1*a4*a5)-1)*(1+a2*a3**2*a4*h5*h6/(h2*h3*a5*a6**2*a7)))
        W34=((a5*a6/h5)*(s1+h3*h5/(s1*a3*a4*a5*a6)-h3*h4*h5/(s1**2*a3*a4**2*a5**2*a6))*(1+a3*a4*h5/(h3*a5*a6)))
        return W23*W24+W23*W34+W24*W34
    except Exception:
        return mp.mpf('nan')

def adm(v1,x1,t1):
    a=fullvec(v1,x1,t1)
    if not all(0<float(a[i])<1 for i in range(1,7)): return False,None,None
    if not all(float(a[i]+a[i+1])>1 for i in range(7)): return False,None,None
    pal=abs(float((1+v1+(1-v1)*x1)/2)-float((1+v1-(1-v1)*x1)/2))<1e-3
    return True,pal,a

# (a) strict-interior: cross-check J_tau>0 with mp.diff at a handful
random.seed(7)
strict_chk=[]
for t0 in [0.5,0.85,0.95,0.97]:
    for _ in range(60):
        v0=random.uniform(0.1,0.9); x0=random.uniform(0.05,0.95)
        if len([q for q in strict_chk if q[0]==t0])>=2: break
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,t0),fGp(vv,xx,t0)),
                            (mp.mpf(v0),mp.mpf(x0)),tol=1e-22,maxsteps=50)
            v1,x1=float(sol[0]),float(sol[1])
            if not(0<v1<1 and 0<x1<1): continue
            r=max(abs(float(fGm(v1,x1,t0))),abs(float(fGp(v1,x1,t0))))
            if r>1e-10: continue
            ok,pal,a=adm(v1,x1,t0)
            if not ok or pal or v1<0.05: continue
            Jv=float(Jtau_mp(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))
            strict_chk.append((t0,v1,x1,Jv))
        except Exception: pass
print("\n(a) strict-interior J_tau (mp.diff, dps25):",flush=True)
for q in strict_chk:
    print(f"  tau={q[0]:.3f} v={q[1]:.4f} xi={q[2]:.4f} J={q[3]:+.5e}",flush=True)
if strict_chk:
    print(f"  all positive: {all(q[3]>0 for q in strict_chk)}",flush=True)

# (b) boundary: seed small v0 to find v->0 non-palindromic curve points, eval Theta
random.seed(55)
bpts=[]
for t0 in [0.3,0.5,0.7,0.85,0.95]:
    for _ in range(80):
        v0=random.uniform(1e-6,0.02); x0=random.uniform(0.02,0.98)
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,t0),fGp(vv,xx,t0)),
                            (mp.mpf(v0),mp.mpf(x0)),tol=1e-20,maxsteps=60)
            v1,x1=float(sol[0]),float(sol[1])
            if not(0<v1<1 and 0<x1<1): continue
            r=max(abs(float(fGm(v1,x1,t0))),abs(float(fGp(v1,x1,t0))))
            if r>1e-8: continue
            ok,pal,a=adm(v1,x1,t0)
            if not ok or pal: continue
            Th=float(Th_direct(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))
            bpts.append((t0,v1,x1,Th))
        except Exception: pass
    print(f"  tau={t0}: {sum(1 for q in bpts if q[0]==t0)} boundary nonpal pts so far",flush=True)

print(f"\n(b) {len(bpts)} boundary (v small) non-palindromic curve pts.",flush=True)
if bpts:
    Ths=[q[3] for q in bpts if q[3]==q[3]]
    pos=sum(1 for t in Ths if t>0); neg=sum(1 for t in Ths if t<0); nz=sum(1 for t in Ths if abs(t)<1e-6)
    print(f"  Theta at boundary: pos={pos} neg={neg} |Theta|<1e-6={nz} min|Theta|={min(abs(t) for t in Ths):.3e}",flush=True)
    for q in sorted(bpts,key=lambda q:q[1])[:15]:
        print(f"  tau={q[0]:.3f} v={q[1]:.6f} xi={q[2]:.4f} Theta={q[3]:+.4e}",flush=True)
print("DONE",flush=True)
