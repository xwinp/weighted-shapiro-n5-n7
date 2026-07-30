#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BOUNDED analytic Puiseux slope for tau->0 corner branch (v2, BUG-FIXED).
Previous version summed all (V,z)-monomial coeffs into one Rational, destroying
the V,z structure and falsely reporting constant leading forms.  Fixed: lead is
kept as a proper polynomial in (V,z).

Parametrize  v = V*s^ell, tau = s, xi = 1 - s^2 * z  (eta = s^2 z, z->4).
For each trial integer ell, lowest s-power of E and T gives leading system
{E_lead(V,z)=0, T_lead(V,z)=0}; an admissible branch needs a real solution
V>0, z>0 (z near 4).  Pure polynomial expansion -- bounded memory, no resultant.
"""
import sympy as sp
from sympy import Poly, symbols, expand, factor
import mpmath as mp
mp.mp.dps=30
v,xi,tau,z,s,V=symbols('v xi tau z s V')
c,d,sg=symbols('c d sigma')
vv=c+d-1
Lc=vv+sg*c**2; Uc=vv-sg*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+sg*d**2; Ud=vv-sg*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - sg*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - sg*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(expand(Fc),c,d,sg,domain=sp.QQ); Pd=Poly(expand(Fd),c,d,sg,domain=sp.QQ)
FL=Pc.exquo(Poly((c-1)*Lc,c,d,sg,domain=sp.QQ)); FR=Pd.exquo(Poly((d-1)*Ld,c,d,sg,domain=sp.QQ))
Gm=Poly(expand(FL.as_expr()-FR.as_expr()),c,d,sg,domain=sp.QQ).exquo(Poly(c-d,c,d,sg,domain=sp.QQ))
Gp=Poly(expand(FL.as_expr()+FR.as_expr()),c,d,sg,domain=sp.QQ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gmt=expand(sp.fraction(sp.together(Gm.as_expr().subs({c:cc,d:dd,sg:sig})))[0])
Gpt=expand(sp.fraction(sp.together(Gp.as_expr().subs({c:cc,d:dd,sg:sig})))[0])
F=Poly(Gmt,v,xi,tau,domain=sp.QQ); H=Poly(Gpt,v,xi,tau,domain=sp.QQ)
def v_order(P):
    m=None
    for monom,coef in P.terms():
        if coef==0: continue
        if m is None or monom[0]<m: m=monom[0]
    return m
def shift_v(P,k):
    d={}
    for monom,coef in P.terms():
        if coef==0: continue
        if monom[0]>=k:
            nm=list(monom); nm[0]=monom[0]-k; d[tuple(nm)]=d.get(tuple(nm),0)+coef
    return Poly(d,*P.gens,domain=sp.QQ)
F0=shift_v(F,v_order(F)); H0=shift_v(H,v_order(H))
E=Poly(expand(F0.as_expr()),v,xi,tau,domain=sp.QQ)
U=Poly(expand(2*F0.as_expr()-H0.as_expr()),v,xi,tau,domain=sp.QQ); T=shift_v(U,1)

def lowest_s_system(P, ell):
    """Return (kmin, Poly(lead,V,z)) with lead = coeff of lowest s-power, V,z preserved."""
    expr=P.as_expr().subs({v:V*s**ell, tau:s, xi:1-s**2*z})
    Pe=Poly(expand(expr),s,V,z,domain=sp.QQ)
    is_=Pe.gens.index(s); iV=Pe.gens.index(V); iz=Pe.gens.index(z)
    bypow={}
    for monom,coef in Pe.terms():
        if coef==0: continue
        e=monom[is_]
        bypow.setdefault(e,[]).append((monom[iV],monom[iz],coef))
    kmin=min(bypow)
    lead=expand(sum(coef*V**b*z**cz for b,cz,coef in bypow[kmin]))
    return kmin,Poly(lead,V,z,domain=sp.QQ)

def real_pos_roots_2d(Elead,Tlead):
    """Find (V>0,z>0) with Elead=Tlead=0.  Scan z, root Elead in V, check Tlead."""
    Eexpr=Elead.as_expr(); Texpr=Tlead.as_expr()
    Ef=sp.lambdify((V,z),Eexpr,'mpmath')
    Tf=sp.lambdify((V,z),Texpr,'mpmath')
    hits=[]
    for zv in [0.2,0.5,1.0,1.5,2.0,3.0,4.0,5.0,6.0,8.0,10.0,15.0,20.0]:
        # build univariate poly in V by substituting z=zv (rational) into the expression
        Esub=sp.Poly(sp.expand(Eexpr.subs({z:sp.Rational(zv).limit_denominator(1000)})),V,domain=sp.QQ)
        coeffs=Esub.all_coeffs()
        if len(coeffs)<=1: continue
        try:
            rr=mp.polyroots([mp.mpf(float(cc)) for cc in coeffs],maxsteps=120,extraprec=30)
        except Exception:
            continue
        for root in rr:
            rv=complex(root)
            if abs(rv.imag)<1e-6*max(1,abs(rv.real)) and rv.real>1e-9:
                Vv=rv.real; tv=complex(Tf(Vv,zv))
                scale=abs(complex(Ef(Vv,zv)))
                if abs(tv)<1e-6*max(1,scale):
                    hits.append((zv,Vv,float(tv.real)))
    return hits

for ell in range(1,8):
    ke,Elead=lowest_s_system(E,ell)
    kt,Tlead=lowest_s_system(T,ell)
    print("===== ell=%d : E lowest s^%d , T lowest s^%d ====="%(ell,ke,kt),flush=True)
    print("  E_lead deg_V=%d deg_z=%d  T_lead deg_V=%d deg_z=%d"%(
        Elead.degree(V),Elead.degree(z),Tlead.degree(V),Tlead.degree(z)),flush=True)
    print("  E_lead =",factor(Elead.as_expr()),flush=True)
    print("  T_lead =",factor(Tlead.as_expr()),flush=True)
    hits=real_pos_roots_2d(Elead,Tlead)
    print("  common (z>0,V>0) hits:",hits[:6],flush=True)
print("DONE",flush=True)
