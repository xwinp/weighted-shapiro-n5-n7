#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Second blow-up: leading edge is rank-1 degenerate.
F0(0)=8(xi^2-1)^5, H0(0)=16(xi^2-1)^5  -> proportional (ratio 1:2), NOT equal.
Common leading part killed by combo U = 2F0 - H0  (= A+3B), which vanishes at v=0.
p = nu_v(U); U_p = leading v-part of U.  U_p(eta,tau)=0 is the true corner-approach
transverse relation (replaces the WRONG F1-H1=0 from the equal-leading assumption).
Longitudinal direction W=B gives v ~ eta^5.
"""
import sympy as sp
from sympy import Poly, symbols, ZZ, factor, real_roots
v,xi,tau,eta=symbols('v xi tau eta')
c,d,sg=symbols('c d sigma')
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
    m=None; idx=P.gens.index(var)
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if m is None or e<m: m=e
    return m if m is not None else 0
def clear_vadic(P,var,k):
    idx=P.gens.index(var); d={}
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if e>=k:
            nm=list(monom); nm[idx]=e-k; d[tuple(nm)]=d.get(tuple(nm),0)+coef
    return Poly(d,*P.gens,domain=ZZ)

nuF=v_order(F,v); nuH=v_order(H,v)
F0=clear_vadic(F,v,nuF); H0=clear_vadic(H,v,nuH)
# U = 2F0 - H0  (kills common leading part)
U=Poly(sp.expand(2*F0.as_expr()-H0.as_expr()),v,xi,tau,domain=ZZ)
p=v_order(U,v)
print("nu_v(F)=",nuF,"nu_v(H)=",nuH,flush=True)
print("p = nu_v(2F0-H0) =",p,flush=True)
Up=clear_vadic(U,v,p)
# Up at v=0 (should be nonzero in v now)
print("v_order(Up)=",v_order(Up,v),flush=True)
Up0=Poly(sp.expand(Up.as_expr().subs({v:0})),xi,tau,domain=ZZ)
print("\nUp(0,xi,tau) = leading transverse part:",flush=True)
print("  deg_xi",Up0.degree(xi),"deg_tau",Up0.degree(tau),"terms",len(Up0.as_dict()),flush=True)
print("  factored:",factor(Up0.as_expr()),flush=True)
# substitute xi=1-eta, look near eta=0
Up0e=Poly(sp.expand(Up0.as_expr().subs({xi:1-eta})),eta,tau,domain=ZZ)
print("\nUp(0,1-eta,tau) factored:",factor(Up0e.as_expr()),flush=True)
# eta-order at eta=0
eo=v_order(Up0e,eta)
print("eta-order of Up at eta=0:",eo,flush=True)
# value at xi=1 (eta=0) as poly in tau  -> vanishes (xi-1)^4 factor; extract R(1,tau)
Rexpr=factor(Up0.as_expr())
print("\nUp(0,xi,tau) factored:",Rexpr,flush=True)
# pull R(xi,tau) = Up0 / (16 (xi-1)^4 (xi+1)^4)
Up0p=Up0
for fac,k in [((xi-1),4),((xi+1),4)]:
    for _ in range(k):
        Up0p=Up0p.exquo(Poly(fac,xi,tau,domain=ZZ))
Up0p=Up0p.exquo(Poly(16,xi,tau,domain=ZZ))
R=Poly(Up0p.as_expr(),xi,tau,domain=ZZ)
print("R(xi,tau) = Up0/[16(xi-1)^4(xi+1)^4]:",flush=True)
print("  R =",R.as_expr(),flush=True)
R1=Poly(sp.expand(R.as_expr().subs({xi:1})),tau,domain=ZZ)
print("  R(1,tau) =",R1.as_expr(),"  factored:",factor(R1.as_expr()),flush=True)
print("  R(1,tau) real roots in (0,1):",[float(r) for r in real_roots(R1.as_expr(),tau) if 0<float(r)<1],flush=True)
# R(1-eta,tau) expanded, to study branches tau->1/2 and tau->0
Re=Poly(sp.expand(R.as_expr().subs({xi:1-eta})),eta,tau,domain=ZZ)
print("\nR(1-eta,tau) factored:",factor(Re.as_expr()),flush=True)
# eta-order of R at eta=0
print("  eta-order of R at 0:",v_order(Re,eta),flush=True)
# coeff eta^0 of R(1-eta,tau) = R(1,tau); coeff eta^1
def coeff_eta(P,ie):
    # P in (eta,tau); collect coeff of eta^ie as poly in tau
    idx=P.gens.index(eta); tidx=P.gens.index(tau); d={}
    for monom,coef in P.terms():
        if coef==0: continue
        if monom[idx]==ie:
            d[monom[tidx]]=d.get(monom[tidx],0)+coef
    return Poly(d,tau,domain=ZZ) if d else Poly(0,tau,domain=ZZ)
R_e0=coeff_eta(Re,0); R_e1=coeff_eta(Re,1)
print("  R(1-eta,tau) eta^0 coeff (=R(1,tau)):",factor(R_e0.as_expr()),flush=True)
print("  R(1-eta,tau) eta^1 coeff:",factor(R_e1.as_expr()),flush=True)

# Also the longitudinal: B=(F0-H0)/2 leading v^0 part in eta, and its v^1 coeff B1(0,tau)
# v ~ eta^5 relation needs B1(eta,tau) nonzero at eta=0.
Bexpr=sp.expand((F0.as_expr()-H0.as_expr())/2)
B=Poly(Bexpr,v,xi,tau,domain=ZZ)
# B = B0(xi,tau) + v B1 + ... ; B1 = coeff v^1
B1=Poly(sp.expand(clear_vadic(B,v,1).as_expr().subs({v:0}) if v_order(B,v)>=1 else 0),xi,tau,domain=ZZ) if v_order(B,v)>=1 else Poly(0,xi,tau,domain=ZZ)
# simpler: B1 = d(B)/dv at v=0
B1expr=sp.expand(sp.diff(B.as_expr(),v).subs({v:0}))
B1=Poly(B1expr,xi,tau,domain=ZZ)
B1e=Poly(sp.expand(B1.as_expr().subs({xi:1-eta})),eta,tau,domain=ZZ)
B1_eta0=Poly(sp.expand(B1.as_expr().subs({xi:1})),tau,domain=ZZ)
print("\nB1(eta,tau)=coeff v^1 of B, at eta=0 (xi=1):",flush=True)
print("  B1(0,tau) =",B1_eta0.as_expr() if B1_eta0!=0 else 0,flush=True)
# eta-order of B1 at eta=0 (xi=1-eta)
B1e=Poly(sp.expand(B1.as_expr().subs({xi:1-eta})),eta,tau,domain=ZZ)
print("  eta-order of B1 at 0:",v_order(B1e,eta),flush=True)
# leading eta-term of B1 as poly in tau
q=v_order(B1e,eta)
if q is not None and q<10:
    B1lead=coeff_eta(B1e,q)
    print("  B1 leading eta^%d coeff (in tau):"%q,factor(B1lead.as_expr()),flush=True)
    print("    value at tau=1/2:",float(B1lead.eval(sp.Rational(1,2))),flush=True)
    print("    value at tau=0:",float(B1lead.eval(0)),flush=True)
print("DONE",flush=True)
