#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Next-order blow-up at the corner v=0.
F0 = P(xi) + v*F1(xi,tau) + v^2*...,  H0 = P(xi) + v*H1(xi,tau) + ...
(P=8(xi^2-1)^5, since F0(0)=H0(0)=P coincide.)
On the curve (small v>0): subtract => v*(F1-H1)+O(v^2)=0 => limiting branch F1(xi,tau)=H1(xi,tau).
Compute D=F1-H1 = (dF0/dv - dH0/dv)|_{v=0}, analyze structure, roots near xi=1.
"""
import sympy as sp
from sympy import Poly, symbols, ZZ, real_roots, factor, gcd
v,xi,tau=symbols('v xi tau'); c,d,s=symbols('c d sigma')
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
F0=F.exquo(Poly(v**3,v,xi,tau,domain=ZZ)); H0=H.exquo(Poly(v**3,v,xi,tau,domain=ZZ))
# F1 = coeff of v^1 in F0  (= dF0/dv at v=0)
dF0=Poly(sp.diff(F0.as_expr(),v),v,xi,tau,domain=ZZ).eval({v:0})
dH0=Poly(sp.diff(H0.as_expr(),v),v,xi,tau,domain=ZZ).eval({v:0})
F1=Poly(dF0,xi,tau,domain=ZZ); H1=Poly(dH0,xi,tau,domain=ZZ)
print("F1: deg_xi",F1.degree(xi),"deg_tau",F1.degree(tau),"terms",len(F1.as_dict()),flush=True)
print("H1: deg_xi",H1.degree(xi),"deg_tau",H1.degree(tau),"terms",len(H1.as_dict()),flush=True)
D=Poly(sp.expand(F1.as_expr()-H1.as_expr()),xi,tau,domain=ZZ)
print("\nD=F1-H1: deg_xi",D.degree(xi),"deg_tau",D.degree(tau),"terms",len(D.as_dict()),flush=True)
print("D factored:",sp.factor(D.as_expr()),flush=True)

# D(xi,tau)=0 is the corner-approach relation. Analyze as poly in xi for given tau,
# or eliminate tau.
# tau-degree of D:
print("D deg_tau",D.degree(tau),flush=True)
# If D depends on tau, solve D=0 for tau as function of xi, or Res.
# Try: for several xi near 1, what tau solve D=0?
import mpmath as mp
mp.mp.dps=20
fD=sp.lambdify((xi,tau),D.as_expr(),'mpmath')
print("\nD(xi,tau) sampled near xi=1:",flush=True)
for x0 in [0.5,0.8,0.9,0.95,0.98,0.99,0.999]:
    # solve D(x0,tau)=0 for tau in (0,1)
    found=[]
    for t0 in [0.1,0.3,0.5,0.7,0.9]:
        try:
            sol=mp.findroot(lambda tt:fD(x0,tt),mp.mpf(t0),tol=1e-20,maxsteps=50)
            tf=float(sol)
            if 0<tf<1 and abs(float(fD(x0,tf)))<1e-8:
                if not any(abs(tf-f)<1e-4 for f in found): found.append(tf)
        except Exception: pass
    print(f"  xi={x0}: D=0 at tau={[round(f,4) for f in sorted(found)]}",flush=True)

# Also check D at xi=1 exactly (the corner)
print("\nD(1,tau) =",sp.simplify(D.eval({xi:1}).as_expr()) if D.eval({xi:1})!=0 else 0,flush=True)
print("D(xi,0) =",sp.factor(D.eval({tau:0}).as_expr()),flush=True)
print("DONE",flush=True)
