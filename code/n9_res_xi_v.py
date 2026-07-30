#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Certificate step 2 (branch enumeration prep): compute projection resultants.
Res_xi(G^-,G^+) -> poly in (v,tau); Res_v(G^-,G^+) -> poly in (xi,tau).
Report degree, term count, adic valuations, factor structure (square-free part).
This chooses the main projection (answers Q B locally) and seeds Sturm isolation.
G^-,G^+ are the cleared-denominator regularized polys Gmt,Gpt in (v,xi,tau).
"""
import sympy as sp
from sympy import Poly, resultant, symbols
v,xi,tau=symbols('v xi tau')
c,d,s=symbols('c d sigma')
vv=c+d-1
Lc=vv+s*c**2; Uc=vv-s*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+s*d**2; Ud=vv-s*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - s*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - s*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(sp.expand(Fc),c,d,s,domain=sp.ZZ); Pd=Poly(sp.expand(Fd),c,d,s,domain=sp.ZZ)
FL=Pc.exquo(Poly((c-1)*Lc,c,d,s,domain=sp.ZZ)); FR=Pd.exquo(Poly((d-1)*Ld,c,d,s,domain=sp.ZZ))
Gm=Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,s,domain=sp.ZZ).exquo(Poly(c-d,c,d,s,domain=sp.ZZ))
Gp=Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,s,domain=sp.ZZ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,s:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,s:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])
Pmt=Poly(Gmt,v,xi,tau,domain=sp.ZZ); Ppt=Poly(Gpt,v,xi,tau,domain=sp.ZZ)
print("G^-: total_deg",Pmt.total_degree(),"terms",len(Pmt.as_dict()),
      "deg_v",Pmt.degree(v),"deg_xi",Pmt.degree(xi),"deg_tau",Pmt.degree(tau),flush=True)
print("G^+: total_deg",Ppt.total_degree(),"terms",len(Ppt.as_dict()),
      "deg_v",Ppt.degree(v),"deg_xi",Ppt.degree(xi),"deg_tau",Ppt.degree(tau),flush=True)

import sys
def adic(p,var):
    # valuation in var (largest k with var^k dividing p)
    P=Poly(p,var,domain=sp.ZZ)
    return P.eval(0)  # not val; compute below
def val(P,var):
    k=0
    while True:
        r=P.eval({var:0}) if False else P.eval(0)
        if r!=0: break
        P=P.exquo(Poly(var,var,domain=sp.ZZ)); k+=1
        if k>200: break
    return k

print("\ncomputing Res_xi(G^-,G^+) ...",flush=True)
try:
    Rxi=resultant(Gmt,Gpt,xi)
    Rx=Poly(Rxi,v,tau,domain=sp.ZZ)
    print("Res_xi: total_deg",Rx.total_degree(),"terms",len(Rx.as_dict()),
          "deg_v",Rx.degree(v),"deg_tau",Rx.degree(tau),flush=True)
    vx=val(Rx,v); vtau=val(Rx,tau)
    print("  v-adic val",vx," tau-adic val",vtau,flush=True)
    Rxr=Rx.exquo(Poly(v**vx*tau**vtau,v,tau,domain=sp.ZZ)) if (vx or vtau) else Rx
    print("  reduced: deg_v",Rxr.degree(v),"deg_tau",Rxr.degree(tau),"terms",len(Poly(Rxr,v,tau,domain=sp.ZZ).as_dict()),flush=True)
    open("paper/_gpt_artifacts/Res_xi_vtau.txt","w").write(str(Rx.as_expr()))
    print("  saved Res_xi.",flush=True)
except Exception as e:
    print("Res_xi FAILED:",repr(e),flush=True)

print("\ncomputing Res_v(G^-,G^+) ...",flush=True)
try:
    Rv=resultant(Gmt,Gpt,v)
    Rvp=Poly(Rv,xi,tau,domain=sp.ZZ)
    print("Res_v: total_deg",Rvp.total_degree(),"terms",len(Rvp.as_dict()),
          "deg_xi",Rvp.degree(xi),"deg_tau",Rvp.degree(tau),flush=True)
    vxi=val(Rvp,xi); vtau2=val(Rvp,tau)
    print("  xi-adic val",vxi," tau-adic val",vtau2,flush=True)
    open("paper/_gpt_artifacts/Res_v_xitau.txt","w").write(str(Rvp.as_expr()))
    print("  saved Res_v.",flush=True)
except Exception as e:
    print("Res_v FAILED:",repr(e),flush=True)
print("DONE",flush=True)
