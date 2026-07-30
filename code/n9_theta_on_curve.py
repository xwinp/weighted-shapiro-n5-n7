#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B.16 decisive test via the 1D curve: find {G^-=G^+=0} points (fix tau, solve
v,xi -- FAST 2-eq), keep strict-interior non-palindromic, evaluate Theta.
A non-palindromic local min needs Theta=0 (necessary). So if Theta has NO sign
change along the Omega-portion of the curve, no Theta=0 point -> B.16 evidence.
Also: directly hunt Theta=0 by tracking sign of Theta across tau-slices.
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

def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]

random.seed(123)
pts=[]   # (tau, v, xi, Theta, sum_a)
for t0 in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]:
    for _ in range(50):
        v0=random.uniform(0.05,0.95); x0=random.uniform(0.02,0.98)
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,t0),fGp(vv,xx,t0)),
                            (mp.mpf(v0),mp.mpf(x0)),tol=1e-20,maxsteps=50)
            v1,x1=float(sol[0]),float(sol[1])
            if not(0<v1<1 and 0<x1<1): continue
            r=max(abs(float(fGm(v1,x1,t0))),abs(float(fGp(v1,x1,t0))))
            if r>1e-8: continue
            a=fullvec(v1,x1,t0)
            if not all(0<float(a[i])<1 for i in range(1,7)): continue
            if not all(float(a[i]+a[i+1])>1 for i in range(7)): continue
            pal=abs(float((1+v1+(1-v1)*x1)/2)-float((1+v1-(1-v1)*x1)/2))<1e-3
            if pal: continue
            Th=float(Th_direct(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))
            pts.append((t0,v1,x1,Th,float(sum(a))))
        except Exception:
            pass
    print(f"  tau={t0}: {sum(1 for p in pts if p[0]==t0)} nonpal-interior curve pts so far",flush=True)

print(f"\n{len(pts)} strict-interior non-palindromic G^-=G^+=0 curve points found.",flush=True)
# Theta sign analysis
import statistics
Ths=[p[3] for p in pts if p[3]==p[3]]  # drop nan
if Ths:
    pos=sum(1 for t in Ths if t>0); neg=sum(1 for t in Ths if t<0); near=sum(1 for t in Ths if abs(t)<1e-6)
    print(f"Theta signs: pos={pos} neg={neg} |Theta|<1e-6={near}  min|Theta|={min(abs(t) for t in Ths):.3e}",flush=True)
print("sample (tau,v,xi,Theta,sum_a):")
for p in sorted(pts,key=lambda q:abs(q[3]))[:12]:
    print(f"  tau={p[0]:.3f} v={p[1]:.4f} xi={p[2]:.4f} Theta={p[3]:+.5e} sum_a={p[4]:.4f}",flush=True)
print("DONE",flush=True)
