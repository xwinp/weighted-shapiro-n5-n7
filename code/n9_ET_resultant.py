#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full transverse elimination per GPT reply#11: substitute the branch blow-up
variable, then eliminate it via resultant between E and T, giving a single
bivariate polynomial whose Newton polygon yields the true Puiseux slope
(handles the leading-coefficient cancellation that w0/z0 alone cannot).

Branch tau->1/2: tau = 1/2 + eta*w.  E,T -> polys in (v,eta,w).  Res_w(E,T) in (v,eta).
Branch tau->0:   eta = tau^2 * z.     E,T -> polys in (v,tau,z).  Res_z(E,T) in (v,tau).
Then Newton polygon of the resultant -> v ~ eta^kappa (resp. v ~ tau^ell), and the
edge initial form must have a real POSITIVE root for an admissible branch.
"""
import sympy as sp
from sympy import Poly, symbols, factor, Rational, real_roots
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
E=Poly(sp.expand(F0.as_expr()),v,xi,tau,domain=sp.QQ)
U=Poly(sp.expand(2*F0.as_expr()-H0.as_expr()),v,xi,tau,domain=sp.QQ)
T=shift_v(U,v,1)

def lower_convex_hull(pts):
    pts=sorted(set(pts))
    if len(pts)<=2: return pts
    hull=[]
    for p in pts:
        while len(hull)>=2:
            a=hull[-2]; b=hull[-1]
            if (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])<=0: hull.pop()
            else: break
        hull.append(p)
    return hull
def newton(P,var2):
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

def report(P,var2,label):
    print("\n--- %s ---"%label,flush=True)
    print("  resultant terms:",len(P.as_dict())," deg_v",P.degree(v)," deg",var2,P.degree(var2),flush=True)
    hull,edges=newton(P,var2)
    print("  lower hull:",hull,flush=True)
    for sl,a,b,terms in edges:
        kappa=-sl if sl not in (None,0) else None
        print("    edge",a,"->",b,"slope",sl," kappa(v~var^kappa)=",-sl if sl else sl,flush=True)
        # initial form in V where v = V * var^kappa.  Along edge, v^i var^j with j+sl*i=const.
        # substitute v=V*var^(-sl): exponent of var = j + (-sl)*i = const. coeff * V^i.
        if sl and sl<0:
            polyV=sum(coef*v**i for (i,j),coef in terms.items())
            polyV=sp.expand(polyV)
            print("      initial form in V:",sp.factor(polyV),flush=True)
            rts=[float(r) for r in real_roots(polyV,v) if r>0]
            print("      positive real V roots:",rts,flush=True)

# ---- branch tau=1/2: tau=1/2+eta*w, eliminate w ----
print("Building branch tau->1/2: tau=1/2+eta*w",flush=True)
Eh=Poly(sp.expand(E.as_expr().subs({xi:1-eta,tau:Rational(1,2)+eta*w})),v,eta,w,domain=sp.QQ)
Th=Poly(sp.expand(T.as_expr().subs({xi:1-eta,tau:Rational(1,2)+eta*w})),v,eta,w,domain=sp.QQ)
print("  deg_w E=",Eh.degree(w)," deg_w T=",Th.degree(w),flush=True)
try:
    R=sp.resultant(Eh.as_expr(),Th.as_expr(),w)
    RP=Poly(sp.expand(R),v,eta,domain=sp.QQ)
    report(RP,eta,"branch tau->1/2: Res_w(E,T) in (v,eta)")
except Exception as e:
    print("  Res_w failed:",repr(e),flush=True)

# ---- branch tau=0: eta=tau^2*z, eliminate z ----
print("\nBuilding branch tau->0: eta=tau^2*z",flush=True)
Ez=Poly(sp.expand(E.as_expr().subs({xi:1-tau**2*z})),v,tau,z,domain=sp.QQ)   # xi=1-eta=1-tau^2 z
Tz=Poly(sp.expand(T.as_expr().subs({xi:1-tau**2*z})),v,tau,z,domain=sp.QQ)
print("  deg_z E=",Ez.degree(z)," deg_z T=",Tz.degree(z),flush=True)
try:
    Rz=sp.resultant(Ez.as_expr(),Tz.as_expr(),z)
    RPz=Poly(sp.expand(Rz),v,tau,domain=sp.QQ)
    report(RPz,tau,"branch tau->0: Res_z(E,T) in (v,tau)")
except Exception as e:
    print("  Res_z failed:",repr(e),flush=True)
print("DONE",flush=True)
