#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GPT reply#9 step 1-2: exact v-adic orders of F(=Gmt), H(=Gpt); v=0 blow-up and
boundary limiting system. Also report structure of F(0,xi,tau), H(0,xi,tau).
J_tau, Theta v-adic orders deferred (need symbolic exprs; do separately).
"""
import sympy as sp
from sympy import Poly, resultant, symbols, factor, gcd, ZZ
v,xi,tau=symbols('v xi tau')
c,d,s=symbols('c d sigma')
vv=c+d-1
Lc=vv+s*c**2; Uc=vv-s*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+s*d**2; Ud=vv-s*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - s*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - s*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(sp.expand(Fc),c,d,s,domain=ZZ); Pd=Poly(sp.expand(Fd),c,d,s,domain=ZZ)
FL=Pc.exquo(Poly((c-1)*Lc,c,d,s,domain=ZZ)); FR=Pd.exquo(Poly((d-1)*Ld,c,d,s,domain=ZZ))
Gm=Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,s,domain=ZZ).exquo(Poly(c-d,c,d,s,domain=ZZ))
Gp=Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,s,domain=ZZ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,s:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,s:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])
F=Poly(Gmt,v,xi,tau,domain=ZZ); H=Poly(Gpt,v,xi,tau,domain=ZZ)
print("F: terms",len(F.as_dict()),"deg_v",F.degree(v),"deg_xi",F.degree(xi),"deg_tau",F.degree(tau),flush=True)
print("H: terms",len(H.as_dict()),"deg_v",H.degree(v),"deg_xi",H.degree(xi),"deg_tau",H.degree(tau),flush=True)

def vadic(P):
    k=0
    cur=P
    while True:
        r=cur.eval({v:0})
        if r!=0: break
        cur=cur.exquo(Poly(v,v,xi,tau,domain=ZZ)); k+=1
        if k>300: break
    return k,cur

kF,F0=vadic(F); kH,H0=vadic(H)
print("\nv-adic: nu_v(F)=",kF," nu_v(H)=",kH,flush=True)
print("F0=v^-kF*F: deg_v",F0.degree(v),"deg_xi",F0.degree(xi),"deg_tau",F0.degree(tau),"terms",len(F0.as_dict()),flush=True)
print("H0=v^-kH*H: deg_v",H0.degree(v),"deg_xi",H0.degree(xi),"deg_tau",H0.degree(tau),"terms",len(H0.as_dict()),flush=True)

# boundary limiting system at v=0
Fb=F0.eval({v:0}); Hb=H0.eval({v:0})
PFb=Poly(Fb,xi,tau,domain=ZZ); PHb=Poly(Hb,xi,tau,domain=ZZ)
print("\nBoundary limiting system F0(0,xi,tau), H0(0,xi,tau):",flush=True)
print("  F0(0): terms",len(PFb.as_dict()),"deg_xi",PFb.degree(xi),"deg_tau",PFb.degree(tau),flush=True)
print("  H0(0): terms",len(PHb.as_dict()),"deg_xi",PHb.degree(xi),"deg_tau",PHb.degree(tau),flush=True)
g=gcd(PFb,PHb)
print("  gcd(F0(0),H0(0)) terms:",len(g.as_dict()) if g else 0,flush=True)

# eliminate xi on boundary: Res_xi(F0(0), H0(0)) -> poly in tau
print("\ncomputing Res_xi(F0(0),H0(0)) -> poly in tau ...",flush=True)
try:
    Rb=resultant(PFb.as_expr(),PHb.as_expr(),xi)
    PRb=Poly(Rb,tau,domain=ZZ)
    print("  Res_xi(F0(0),H0(0)): degree",PRb.degree(tau),"terms",len(PRb.as_dict()),flush=True)
    # tau-adic valuation
    k=0; cur=PRb
    while True:
        if cur.eval({tau:0})!=0: break
        cur=cur.exquo(Poly(tau,tau,domain=ZZ)); k+=1
        if k>200: break
    print("  tau-adic val",k,flush=True)
    PRr=cur
    print("  reduced degree",PRr.degree(tau),"terms",len(PRr.as_dict()),flush=True)
    # count real roots in (0,1)
    n01=PRr.count_roots(0,1)
    print("  real roots in (0,1):",n01,flush=True)
    open("paper/_gpt_artifacts/boundary_Res_xi_F0H0_tau.txt","w").write(str(PRb.as_expr()))
    print("  saved.",flush=True)
except Exception as e:
    print("  FAILED:",repr(e),flush=True)

# Also: palindromic boundary xi=0 -- what is F(xi=0),H(xi=0)?
print("\nxi=0 (palindromic diagonal) check:",flush=True)
Fi=Poly(F.eval({xi:0}),v,tau,domain=ZZ); Hi=Poly(H.eval({xi:0}),v,tau,domain=ZZ)
print("  F(xi=0) terms",len(Fi.as_dict()),"deg_v",Fi.degree(v),"deg_tau",Fi.degree(tau),flush=True)
g2=gcd(Fi,Hi)
print("  gcd(F(xi=0),H(xi=0)) terms:",len(g2.as_dict()) if g2 else 0,flush=True)
print("DONE",flush=True)
