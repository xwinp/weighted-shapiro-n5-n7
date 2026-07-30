#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GPT reply#11 prescribed corner analysis via E=F0, T=(2F0-H0)/v.
E(0,eta,tau) = -8 eta^5 (2-eta)^5            (5th-order longitudinal, tau-free)
T(0,eta,tau) = 16 eta^4 (2-eta)^4 R(1-eta,tau)  (4th-order transverse)
Branch tau->1/2: mu=tau-1/2=eta*w, w0=-1/4 (R_xi/R_tau at (1,1/2)); subst into E,
  Newton polygon (v,eta) -> kappa (v ~ eta^kappa).
Branch tau->0:   eta=tau^2*z, z0=4; subst into E, Newton polygon (v,tau) -> ell.
Also: list common lower Newton faces of E,T in (v,eta) for the full-face check.
"""
import sympy as sp
from sympy import Poly, symbols, factor, Rational
v,xi,tau,eta,w,z=symbols('v xi tau eta w z')
c,d,sg=symbols('c d sigma')
vv=c+d-1
Lc=vv+sg*c**2; Uc=vv-sg*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+sg*d**2; Ud=vv-sg*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - sg*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - sg*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(sp.expand(Fc),c,d,sg,domain=sp.QQ); Pd=Poly(sp.expand(Fd),c,d,sg,domain=sp.QQ)
FL=Pc.exquo(Poly((c-1)*Lc,c,d,sg,domain=sp.QQ)); FR=Pd.exquo(Poly((d-1)*Ld,c,d,sg,domain=sp.QQ))
Gm=Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,sg,domain=sp.QQ).exquo(Poly(c-d,c,d,sg,domain=sp.QQ))
Gp=Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,sg,domain=sp.QQ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,sg:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,sg:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])
F=Poly(Gmt,v,xi,tau,domain=sp.QQ); H=Poly(Gpt,v,xi,tau,domain=sp.QQ)
def v_order(P,var):
    m=None; idx=P.gens.index(var)
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if m is None or e<m: m=e
    return m if m is not None else 0
def shift_v(P,var,k):
    idx=P.gens.index(var); d={}
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if e>=k:
            nm=list(monom); nm[idx]=e-k; d[tuple(nm)]=d.get(tuple(nm),0)+coef
    return Poly(d,*P.gens,domain=sp.QQ)
nuF=v_order(F,v); nuH=v_order(H,v)
F0=shift_v(F,v,nuF); H0=shift_v(H,v,nuH)
# E=F0, T=(2F0-H0)/v  (nu_v(2F0-H0)=1)
E=Poly(sp.expand(F0.as_expr()),v,xi,tau,domain=sp.QQ)
U=Poly(sp.expand(2*F0.as_expr()-H0.as_expr()),v,xi,tau,domain=sp.QQ)
T=shift_v(U,v,1)   # divide by v
# to (v,eta,tau)
Ee=Poly(sp.expand(E.as_expr().subs({xi:1-eta})),v,eta,tau,domain=sp.QQ)
Te=Poly(sp.expand(T.as_expr().subs({xi:1-eta})),v,eta,tau,domain=sp.QQ)
print("nu_v(E)=",v_order(Ee,v)," nu_v(T)=",v_order(Te,v),flush=True)
print("E(0,eta,tau) =",factor(Ee.eval({v:0}).as_expr()),flush=True)
print("T(0,eta,tau) =",factor(Te.eval({v:0}).as_expr()),flush=True)

# ---- common Newton faces of E,T in (v,eta): support (deg_v,deg_eta) ----
def lower_frontier(P):
    iv=P.gens.index(v); ie=P.gens.index(eta); byv={}
    for monom,coef in P.terms():
        if coef==0: continue
        ivv=monom[iv]; iee=monom[ie]
        byv[ivv]=min(byv.get(ivv,iee),iee)
    return sorted(byv.items())
print("\nE lower frontier (deg_v->min deg_eta):",lower_frontier(Ee),flush=True)
print("T lower frontier (deg_v->min deg_eta):",lower_frontier(Te),flush=True)

# ---- branch tau=1/2: subst tau=1/2 + eta*w, w0=-1/4 ----
print("\n===== branch tau->1/2 =====",flush=True)
# verify w0 = R_xi/R_tau at (1,1/2)
Rexpr=4*tau**3*xi**2+12*tau**3-8*tau**2*xi**2+4*tau*xi**2-4*tau-xi**2+1
Rxi=sp.diff(Rexpr,xi); Rtau=sp.diff(Rexpr,tau)
w0=Rxi.subs({xi:1,tau:Rational(1,2)})/Rtau.subs({xi:1,tau:Rational(1,2)})
print("w0 = R_xi/R_tau @ (1,1/2) =",w0,flush=True)
# subst tau = 1/2 + eta*w0  into E (initial transverse solution)
Ehalf=Poly(sp.expand(Ee.as_expr().subs({tau:Rational(1,2)+eta*w0})),v,eta,domain=sp.QQ)
print("E(tau=1/2-eta/4) lower frontier:",lower_frontier(Ehalf),flush=True)
print("  E(0,eta) factored:",factor(Ehalf.eval({v:0}).as_expr()),flush=True)
# Newton polygon: correct lower convex hull + edge initial forms.
def lower_convex_hull(pts):
    pts=sorted(set(pts))
    if len(pts)<=2: return pts
    hull=[]
    for p in pts:
        while len(hull)>=2:
            a=hull[-2]; b=hull[-1]
            cross=(b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
            if cross<=0: hull.pop()
            else: break
        hull.append(p)
    return hull
def newton_analysis(P,var2):
    iv=P.gens.index(v); i2=P.gens.index(var2); byv={}
    for monom,coef in P.terms():
        if coef==0: continue
        byv[monom[iv]]=min(byv.get(monom[iv],monom[i2]),monom[i2])
    hull=lower_convex_hull(list(byv.items()))
    edges=[]
    for k in range(len(hull)-1):
        a=hull[k]; b=hull[k+1]
        sl=(b[1]-a[1])/(b[0]-a[0]) if b[0]!=a[0] else None
        di=b[0]-a[0]; dj=b[1]-a[1]; const=dj*a[0]-di*a[1]
        terms={}
        for monom,coef in P.terms():
            if coef==0: continue
            if dj*monom[iv]-di*monom[i2]==const:
                terms[(monom[iv],monom[i2])]=terms.get((monom[iv],monom[i2]),0)+coef
        edges.append((sl,a,b,terms))
    return hull,edges
hull,edges=newton_analysis(Ehalf,eta)
print("  E(tau=1/2-eta/4) lower hull:",hull,flush=True)
for sl,a,b,terms in edges:
    kappa=-1/sl if sl not in (None,0) else None
    print("    edge",a,"->",b,"slope",sl,"kappa=",kappa,flush=True)
    print("      initial-form terms:",{k:sp.Rational(t) for k,t in terms.items()},flush=True)

# ---- branch tau=0: subst eta=tau^2*z, z0=4 ----
print("\n===== branch tau->0 =====",flush=True)
z0=4
Ezero=Poly(sp.expand(Ee.as_expr().subs({eta:4*tau**2})),v,tau,domain=sp.QQ)
print("  E(0,tau) factored:",factor(Ezero.eval({v:0}).as_expr()),flush=True)
hull2,edges2=newton_analysis(Ezero,tau)
print("  E(eta=4tau^2) lower hull (v,tau):",hull2,flush=True)
for sl,a,b,terms in edges2:
    ell=-1/sl if sl not in (None,0) else None
    print("    edge",a,"->",b,"slope",sl,"ell=",ell,flush=True)
    print("      initial-form terms:",{k:sp.Rational(t) for k,t in terms.items()},flush=True)
print("DONE",flush=True)
