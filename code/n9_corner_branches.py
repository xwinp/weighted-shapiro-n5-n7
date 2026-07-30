#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Trace corner branches by fixing tau (near 0 and near 1/2), solving (v,xi).
Find small-v solutions with xi->1, record Theta. Confirms tau->{0,1/2} & Theta<0.
"""
import random, mpmath as mp, sympy as sp
mp.mp.dps=15
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
    W23=((a6*a7/h6)*(s1+h2*h6/(s1*a2*a3*a6*a7)-h2*h3*h4*h5*h6/(s1**4*a2*a3**2*a4**2*a5**2*a6**2*a7))*(1+a2*a3*h6/(h2*a6*a7)))
    W24=(h2*h3/(s1**2*a2*a3**2*a4)*(h4/(s1*a4*a5)-1)*(1+a2*a3**2*a4*h5*h6/(h2*h3*a5*a6**2*a7)))
    W34=((a5*a6/h5)*(s1+h3*h5/(s1*a3*a4*a5*a6)-h3*h4*h5/(s1**2*a3*a4**2*a5**2*a6))*(1+a3*a4*h5/(h3*a5*a6)))
    return W23*W24+W23*W34+W24*W34
def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]

random.seed(7)
# tau grid: near 0 and near 1/2
taus=[0.001,0.003,0.005,0.01,0.02,0.03,0.05,0.08,0.1,0.15,0.2,
      0.40,0.45,0.48,0.49,0.495,0.499,0.4999,0.5001,0.501,0.505,0.51,0.52,0.55,0.6]
rows=[]
for t0 in taus:
    found=[]
    for _ in range(35):
        v0=random.uniform(1e-5,0.5); x0=random.uniform(0.5,0.999)
        try:
            sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,mp.mpf(t0)),fGp(vv,xx,mp.mpf(t0))),(mp.mpf(v0),mp.mpf(x0)),tol=1e-18,maxsteps=40)
            v1,x1=float(sol[0]),float(sol[1])
            if not(1e-7<v1<1 and 0<x1<1): continue
            r=max(abs(float(fGm(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))),abs(float(fGp(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))))
            if r>1e-9: continue
            a=fullvec(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0))
            if not all(0<float(a[i])<1 for i in range(1,7)): continue
            if not all(float(a[i]+a[i+1])>1 for i in range(7)): continue
            pal=abs(float((cc.subs({v:mp.mpf(v1),xi:mp.mpf(x1)})-dd.subs({v:mp.mpf(v1),xi:mp.mpf(x1)}))))<1e-3
            if pal: continue
            Th=float(Th_direct(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))
            found.append((v1,x1,Th))
        except Exception: pass
    # keep smallest-v solution
    if found:
        found.sort(key=lambda q:q[0])
        b=found[0]
        rows.append((t0,b[0],b[1],1-b[1],b[2],len(found)))
        print(f"tau={t0:.4f}: smallest v={b[0]:.3e} xi={b[1]:.5f} eta={1-b[1]:.2e} Theta={b[2]:+.4e} (nfound={len(found)})",flush=True)
    else:
        print(f"tau={t0:.4f}: no small-v nonpal solution",flush=True)
print("\n--- summary: eta vs tau (does eta->0 only at tau->0,1/2?) ---",flush=True)
for r in rows:
    print(f"  tau={r[0]:.4f} eta={r[3]:.3e} v={r[1]:.3e} Theta={r[4]:+.3e}",flush=True)
# Theta sign
Ths=[r[4] for r in rows if r[4]==r[4]]
if Ths:
    print(f"\nTheta: pos={sum(1 for t in Ths if t>0)} neg={sum(1 for t in Ths if t<0)} min|Theta|={min(abs(t) for t in Ths):.3e}",flush=True)
print("DONE",flush=True)
