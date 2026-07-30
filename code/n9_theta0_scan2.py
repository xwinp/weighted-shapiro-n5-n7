#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Decisive B.16 test (v2, fast): {G^-=G^+=Theta=0} in (v,xi,tau).
Avoids lambdifying the giant Theta rational: evaluate Theta directly via the
mpmath compact formula at each (c,d,sigma). 3 eqns / 3 unknowns findroot.
If zero strict-interior non-palindromic solutions -> B.16 stands (the saddle
found earlier has Theta=-0.0023 != 0, so it won't satisfy Theta=0).
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
    a1=mp.mpf(1); a8=mp.mpf(1); a4=c1; a5=d1
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
n=0; valid=0; nonpal=0
for t0 in [0.2,0.5,0.8,0.95]:
    cnt=0
    for _ in range(150):
        v0=random.uniform(0.05,0.95); x0=random.uniform(0.02,0.98)
        try:
            sol=mp.findroot(lambda vv,xx,tt:(fGm(vv,xx,tt),fGp(vv,xx,tt),Th_direct(vv,xx,tt)),
                            (mp.mpf(v0),mp.mpf(x0),mp.mpf(t0)),tol=1e-18,maxsteps=35)
            v1,x1,t1=float(sol[0]),float(sol[1]),float(sol[2])
            if not(0<v1<1 and 0<x1<1 and 0<t1<1): continue
            r=max(abs(float(fGm(v1,x1,t1))),abs(float(fGp(v1,x1,t1))),abs(float(Th_direct(v1,x1,t1))))
            if r>1e-6: continue
            a=fullvec(v1,x1,t1)
            ok=all(0<float(a[i])<1 for i in range(1,7)) and all(float(a[i]+a[i+1])>1 for i in range(7))
            n+=1
            if ok:
                valid+=1
                pal=abs(float((1+v1+(1-v1)*x1)/2)-float((1+v1-(1-v1)*x1)/2))<1e-4
                if not pal: nonpal+=1
                print(f"  FOUND tau={t1:.5f} v={v1:.5f} xi={x1:.5f} {'PAL' if pal else 'NONPAL'} Th={float(Th_direct(v1,x1,t1)):.2e} a={[round(float(z),4) for z in a]}",flush=True)
        except Exception:
            pass
        cnt+=1
        if cnt%25==0: print(f"    ...tau={t0} {cnt}/150",flush=True)
    print(f"  [seed tau={t0}: cum {n} sols, {valid} valid, {nonpal} nonpal]",flush=True)
print(f"\nTOTAL Gm=Gp=Th=0: {n} sols, {valid} strict-interior, {nonpal} non-palindromic",flush=True)
print("DONE",flush=True)
