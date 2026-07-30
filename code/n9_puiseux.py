#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Newton-Puiseux blow-up at the corner (v,xi)=(0,1).  GPT reply#10 procedure.
Robust version: v-adic order via term-scan; substitutions via expr.subs.
eta = 1-xi.  m = nu_v(F0-H0).  A=(F0+H0)/2, B=(F0-H0)/(2 v^m).
Newton polygon (v,eta); verify generic edge; v=s^5, eta=s*z; edge system at s=0.
"""
import sympy as sp
from sympy import Poly, symbols, ZZ, factor, real_roots
v,xi,tau,eta,z,s=symbols('v xi tau eta z s')
c,d,sg=symbols('c d sigma')   # sg = sigma
vv=c+d-1
Lc=vv+sg*c**2; Uc=vv-sg*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+sg*d**2; Ud=vv-sg*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - sg*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - sg*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(sp.expand(Fc),c,d,sg,domain=ZZ); Pd=Poly(sp.expand(Fd),c,d,sg,domain=ZZ)
FL=Pc.exquo(Poly((c-1)*Lc,c,d,sg,domain=ZZ)); FR=Pd.exquo(Poly((d-1)*Ld,c,d,sg,domain=ZZ))
Gm=Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,sg,domain=ZZ).exquo(Poly(c-d,c,d,sg,domain=ZZ))
Gp=Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,sg,domain=ZZ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,sg:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,sg:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])
F=Poly(Gmt,v,xi,tau,domain=ZZ); H=Poly(Gpt,v,xi,tau,domain=ZZ)

def v_order(P,var):
    """min exponent of var among all terms; robust multivariate."""
    m=None
    for monom,coef in P.terms():
        if coef==0: continue
        idx=P.gens.index(var)
        e=monom[idx]
        if m is None or e<m: m=e
    return m if m is not None else 0
def clear_vadic(P,var,k):
    # divide P by var^k by dropping terms with exponent<k and shifting
    idx=P.gens.index(var)
    d={}
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if e>=k:
            nm=list(monom); nm[idx]=e-k; d[tuple(nm)]=d.get(tuple(nm),0)+coef
    return Poly(d,*P.gens,domain=ZZ)

nuF=v_order(F,v); nuH=v_order(H,v)
print("nu_v(F)=",nuF," nu_v(H)=",nuH,flush=True)
F0=clear_vadic(F,v,nuF); H0=clear_vadic(H,v,nuH)
# sanity: F0(0) nonzero in v
print("v_order(F0)=",v_order(F0,v)," v_order(H0)=",v_order(H0,v),flush=True)
F0at0=Poly(sp.expand(F0.as_expr().subs({v:0})),xi,tau,domain=ZZ)
H0at0=Poly(sp.expand(H0.as_expr().subs({v:0})),xi,tau,domain=ZZ)
print("F0(0,xi,tau) factored:",factor(F0at0.as_expr()),flush=True)
print("H0(0,xi,tau) factored:",factor(H0at0.as_expr()),flush=True)
DH=Poly(sp.expand(F0.as_expr()-H0.as_expr()),v,xi,tau,domain=ZZ)
m=v_order(DH,v)
print("m = nu_v(F0-H0) =",m,flush=True)

# A=(F0+H0)/2, B=(F0-H0)/(2 v^m)
Aexpr=sp.expand((F0.as_expr()+H0.as_expr())/2)
Bexpr=sp.expand((F0.as_expr()-H0.as_expr())/(2*v**m))
A=Poly(Aexpr,v,xi,tau,domain=ZZ); B=Poly(Bexpr,v,xi,tau,domain=ZZ)
# to (v,eta,tau): xi=1-eta
Ae=Poly(sp.expand(A.as_expr().subs({xi:1-eta})),v,eta,tau,domain=ZZ)
Be=Poly(sp.expand(B.as_expr().subs({xi:1-eta})),v,eta,tau,domain=ZZ)
A0e=Poly(sp.expand(Ae.as_expr().subs({v:0})),eta,tau,domain=ZZ)
B0e=Poly(sp.expand(Be.as_expr().subs({v:0})),eta,tau,domain=ZZ)
print("\nA(0,eta,tau) factored:",factor(A0e.as_expr()),flush=True)
print("B(0,eta,tau) factored:",factor(B0e.as_expr()),flush=True)

# Newton polygon support (deg_v,deg_eta)
def support(P):
    iv_=P.gens.index(v); ie_=P.gens.index(eta)
    pts=set()
    for monom,coef in P.terms():
        if coef==0: continue
        pts.add((monom[iv_],monom[ie_]))
    return pts
sA=support(Ae); sB=support(Be)
# lower hull: for each deg_v, min deg_eta
def lower_frontier(pts):
    byv={}
    for (iv,ie) in pts: byv[iv]=min(byv.get(iv,ie),ie)
    items=sorted(byv.items())
    return items
print("\nA lower frontier (deg_v -> min deg_eta):",flush=True)
print("  ",lower_frontier(sA),flush=True)
print("B lower frontier (deg_v -> min deg_eta):",flush=True)
print("  ",lower_frontier(sB),flush=True)

# blow-up v=s^5, eta=s*z
As=Poly(sp.expand(Ae.as_expr().subs({v:s**5,eta:s*z})),s,z,tau,domain=ZZ)
Bs=Poly(sp.expand(Be.as_expr().subs({v:s**5,eta:s*z})),s,z,tau,domain=ZZ)
kA=v_order(As,s); kB=v_order(Bs,s)
print("\nlowest s-power: A s^%d  B s^%d"%(kA,kB),flush=True)
Ared=clear_vadic(As,s,kA); Bred=clear_vadic(Bs,s,kB)
A0=Poly(sp.expand(Ared.as_expr().subs({s:0})),z,tau,domain=ZZ)
B0=Poly(sp.expand(Bred.as_expr().subs({s:0})),z,tau,domain=ZZ)
print("A0(z,tau): deg_z",A0.degree(z),"deg_tau",A0.degree(tau),"terms",len(A0.as_dict()),flush=True)
print("  A0 =",A0.as_expr(),flush=True)
print("B0(z,tau): deg_z",B0.degree(z),"deg_tau",B0.degree(tau),"terms",len(B0.as_dict()),flush=True)
print("  B0 =",B0.as_expr(),flush=True)
print("\nB0 factored:",factor(B0.as_expr()),flush=True)
print("A0 factored:",factor(A0.as_expr()),flush=True)
try:
    Rz=sp.resultant(A0.as_expr(),B0.as_expr(),z)
    PR=Poly(sp.expand(Rz),tau,domain=ZZ)
    print("\nRes_z(A0,B0): deg_tau",PR.degree(tau),"terms",len(PR.as_dict()),flush=True)
    print("  factored:",factor(PR.as_expr()),flush=True)
    rr=[float(r) for r in real_roots(PR.as_expr(),tau) if 0<float(r)<1]
    print("  real roots in (0,1):",rr,flush=True)
except Exception as e:
    print("Res_z failed:",repr(e),flush=True)
print("DONE",flush=True)
