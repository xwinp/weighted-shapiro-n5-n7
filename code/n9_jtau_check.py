#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Certificate step 1: is tau a global parameter on {G^-=G^+=0} cap Omega_np?
J_tau = det d(v,xi)/d(G^-,G^+) = [[Gm_v,Gm_xi],[Gp_v,Gp_xi]] det, of the cleared
polys. At a curve zero, J_tau(cleared)=den^2 * J_tau(rational), so they vanish
together. If J_tau != 0 everywhere on the curve in Omega_np, every connected
component is an analytic graph (v(tau),xi(tau)) -> huge simplification.
Approach: (a) eval J_tau on the 111-style curve points; (b) 3-eq scan
{G^-=G^+=J_tau=0} for J_tau=0 points.
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
print("diff J_tau...",flush=True)
Jt=sp.expand(sp.diff(Gmt,v)*sp.diff(Gpt,xi)-sp.diff(Gmt,xi)*sp.diff(Gpt,v))
fGm=sp.lambdify((v,xi,tau),Gmt,'mpmath'); fGp=sp.lambdify((v,xi,tau),Gpt,'mpmath')
fJt=sp.lambdify((v,xi,tau),Jt,'mpmath')
print("  lambdified.",flush=True)

def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]

# (a) eval J_tau on curve points found by 2-eq solve
random.seed(321)
Jvals=[]; zero_hits=[]
for t0 in [0.25,0.5,0.75,0.9,0.95]:
    for _ in range(40):
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
            Jv=float(fJt(v1,x1,t1 if False else t0))
            Jvals.append(Jv)
            if abs(Jv)<1e-3: zero_hits.append((t0,v1,x1,Jv))
        except Exception: pass
print(f"(a) {len(Jvals)} nonpal-interior curve pts.",flush=True)
if Jvals:
    pos=sum(1 for x in Jvals if x>0); neg=sum(1 for x in Jvals if x<0)
    print(f"  J_tau signs: pos={pos} neg={neg} min|J|={min(abs(x) for x in Jvals):.3e}",flush=True)
    print(f"  J_tau near-zero (|J|<1e-3): {len(zero_hits)}",flush=True)
    for z in zero_hits[:6]: print(f"    tau={z[0]} v={z[1]:.4f} xi={z[2]:.4f} J={z[3]:.3e}",flush=True)

# (b) 3-eq scan {G^-=G^+=J_tau=0} -- is there ANY point where tau fails as param?
random.seed(99)
n3=0; v3=0
for t0 in [0.3,0.6,0.9]:
    for _ in range(60):
        v0=random.uniform(0.05,0.95); x0=random.uniform(0.02,0.98)
        try:
            sol=mp.findroot(lambda vv,xx,tt:(fGm(vv,xx,tt),fGp(vv,xx,tt),fJt(vv,xx,tt)),
                            (mp.mpf(v0),mp.mpf(x0),mp.mpf(t0)),tol=1e-18,maxsteps=35)
            v1,x1,t1=float(sol[0]),float(sol[1]),float(sol[2])
            if not(0<v1<1 and 0<x1<1 and 0<t1<1): continue
            r=max(abs(float(fGm(v1,x1,t1))),abs(float(fGp(v1,x1,t1))),abs(float(fJt(v1,x1,t1))))
            if r>1e-6: continue
            a=fullvec(v1,x1,t1)
            ok=all(0<float(a[i])<1 for i in range(1,7)) and all(float(a[i]+a[i+1])>1 for i in range(7))
            n3+=1
            if ok:
                pal=abs(float((1+v1+(1-v1)*x1)/2)-float((1+v1-(1-v1)*x1)/2))<1e-3
                if not pal: v3+=1
                print(f"  Jt=0 HIT: tau={t1:.4f} v={v1:.4f} xi={x1:.4f} {'PAL' if pal else 'NONPAL'}",flush=True)
        except Exception: pass
    print(f"  [3-eq seed {t0}: cum {n3} sols, {v3} nonpal-interior]",flush=True)
print(f"\n(b) {{Gm=Gp=Jt=0}}: {n3} sols, {v3} strict-interior non-palindromic",flush=True)
print("DONE",flush=True)
