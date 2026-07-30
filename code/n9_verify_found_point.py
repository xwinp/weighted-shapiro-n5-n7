#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the suspected strict-interior non-palindromic G^-=G^+=0 point found
by the tau-scan (tau~0.95, v~0.1621, xi~0.3563). Re-solve to high precision,
check ALL admissibility, compute Theta and J (criticality). If Theta!=0 or
J!=0, the point is a SADDLE -> B.16 safe, but the "stronger target"
{G^-=G^+=0}cap Omega_np = empty is FALSE (revert to {G^-=G^+=Theta=0}).
"""
import mpmath as mp, sympy as sp
mp.mp.dps = 40

c,d,s = sp.symbols('c d sigma')
v,xi,tau = sp.symbols('v xi tau')
vv = c+d-1
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

# re-solve near (v=0.1621, xi=0.3563, tau=0.95)
t0=mp.mpf('0.95')
sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,t0),fGp(vv,xx,t0)),
                (mp.mpf('0.16'),mp.mpf('0.36')),tol=1e-35,maxsteps=200)
v1,x1=sol[0],sol[1]
res=max(abs(fGm(v1,x1,t0)),abs(fGp(v1,x1,t0)))
print(f"solved: v={v1} xi={x1} tau={t0}  resid={res}", flush=True)
c1=(1+v1+(1-v1)*x1)/2; d1=(1+v1-(1-v1)*x1)/2; s1=t0*v1/(d1*(1-d1)); g=c1+d1-1
a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
a=[1,a2,a3,c1,d1,a6,a7,1]
print(f"c={c1} d={d1} sigma={s1}", flush=True)
print(f"a={[float(z) for z in a]}", flush=True)
print(f"sum a_i (P?) = {sum(a)}", flush=True)
# admissibility
adm=[]
for i in range(1,7):
    if not(0<a[i]<1): adm.append(f"a{i+1} not in (0,1): {float(a[i])}")
for i in range(7):
    h=a[i]+a[i+1]-1
    if h<=0: adm.append(f"h{i+1}=a{i+1}+a{i+2}-1<=0: {float(h)}")
print(f"admissibility fails: {adm if adm else 'NONE -> strictly interior'}", flush=True)
print(f"non-palindromic (c!=d)? {abs(c1-d1)>1e-6}  |c-d|={float(abs(c1-d1))}", flush=True)

# ---- Theta at this point ----
def theta_num_eval(a, sigma):
    (a1,a2,a3,a4,a5,a6,a7,a8)=a
    h2=a2+a3-1; h3=a3+a4-1; h4=a4+a5-1; h5=a5+a6-1; h6=a6+a7-1
    W23=((a6*a7/h6)*(sigma+h2*h6/(sigma*a2*a3*a6*a7)
        -h2*h3*h4*h5*h6/(sigma**4*a2*a3**2*a4**2*a5**2*a6**2*a7))
        *(1+a2*a3*h6/(h2*a6*a7)))
    W24=(h2*h3/(sigma**2*a2*a3**2*a4)*(h4/(sigma*a4*a5)-1)
        *(1+a2*a3**2*a4*h5*h6/(h2*h3*a5*a6**2*a7)))
    W34=((a5*a6/h5)*(sigma+h3*h5/(sigma*a3*a4*a5*a6)
        -h3*h4*h5/(sigma**2*a3*a4**2*a5**2*a6))
        *(1+a3*a4*h5/(h3*a5*a6)))
    return W23*W24+W23*W34+W24*W34
Th = theta_num_eval(a, s1)
print(f"\nTheta = {Th}", flush=True)
print(f"|Theta| = {abs(Th)}", flush=True)
print(f"Theta == 0? {abs(Th) < 1e-8}", flush=True)

# ---- J at this point (criticality) ----
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'
Xs,Ys,Ss=sp.symbols('X Y sigma')
Jload = sp.sympify((HERE/'nonpal_J_clean.txt').read_text(encoding='utf-8').strip(),
                   locals={'X':Xs,'Y':Ys,'s':Ss}).subs(Ss,s)
# evaluate at X=c+d, Y=c*d
fJ=sp.lambdify((Xs,Ys),Jload,'mpmath')
Jv=fJ(c1+d1, c1*d1)
print(f"J(X,Y) = {Jv}   |J|={abs(Jv)}   J==0? {abs(Jv)<1e-6}", flush=True)
print("\nSUMMARY:", flush=True)
print(f"  G^-=G^+=0 strict-interior non-palindromic point: FOUND (resid {float(res):.2e})", flush=True)
print(f"  Theta==0? {abs(Th)<1e-8}    J==0? {abs(Jv)<1e-6}", flush=True)
if abs(Th)>1e-8 or abs(Jv)>1e-6:
    print("  => point is a SADDLE (not a non-palindromic local min): B.16 UNVIOLATED.", flush=True)
    print("  => but {G^-=G^+=0}cap Omega_np is NONEMPTY: the 'stronger target' is FALSE.", flush=True)
    print("  => must revert B.16 target to {G^-=G^+=Theta=0}cap Omega = empty (orig, 198-start empty).", flush=True)
print("DONE", flush=True)
