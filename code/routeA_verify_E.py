import sympy as sp
from sympy import Poly
z=sp.Symbol('z')
Q5=2*z**5+2*z**3-2*z**2-1
Q7=8*z**7-24*z**6+20*z**5-9*z**4+30*z**3-15*z**2-6
print("=== Exact Sturm root counts (sympy Poly.count_roots, no floats) ===")
for name,Q in [("Q5",Q5),("Q7",Q7)]:
    P=Poly(Q,z)
    n_91=P.count_roots(sp.Rational(9,10),1)
    n_full=P.count_roots(0,1)
    print(f"  {name}: roots in (9/10,1) = {n_91};  roots in (0,1) = {n_full}")
print(f"\n  Q5(9/10) = {sp.Rational(Q5.subs(z,sp.Rational(9,10)))} = {float(Q5.subs(z,sp.Rational(9,10))):.6f}")
print(f"  Q7(9/10) = {sp.Rational(Q7.subs(z,sp.Rational(9,10)))} = {float(Q7.subs(z,sp.Rational(9,10))):.6f}")
# verify z-interval of band: compute z at p near a7 and b7 on S1 branch
import numpy as np
b,c,d,e,f,p=sp.symbols('b c d e f p',positive=True);q=1-p
P=1/(p*b+q*c)+b/(p*c+q*d)+c/(p*d+q*e)+d/(p*e+q*f)+e/(p*f)+f/q
gs=[sp.diff(P,v) for v in [b,c,d,e,f]]
def z_at(pp,init):
    gsub=[g.subs(p,sp.Float(pp)) for g in gs]
    sol=sp.nsolve(gsub,[b,c,d,e,f],init,prec=40,tol=1e-35,maxsteps=200)
    bb,cc,dd,ee,ff=[float(s) for s in sol]
    rho5=ff/ee;qv=1-pp
    return qv*rho5/(pp+qv*rho5),sol
a7,b7=0.214273520909841,0.328627677916592
prev=[0.27,0.68,1.5,0.07,1.3]
print("\n=== z(p) on S1 branch near band endpoints ===")
for pp in [a7+0.0005, 0.27, b7-0.0005]:
    zv,sol=z_at(pp,prev);prev=[float(s) for s in sol]
    print(f"  p={pp:.5f}: z={zv:.10f}")
# H_C exclusion: check K on H_C positive branch. Solve H_C=0 for w given z in (zeta,1), compute K, check >1/128
HC=w**3*z+w**2*z**3-w**2*z+w*z**4-3*w*z**3+2*w*z**2+w*z-w-z**4+3*z**3-3*z**2+z  # need w sym
w=sp.Symbol('w')
HC=z*w**3+w**2*z**3-w**2*z+w*z**4-3*w*z**3+2*w*z**2+w*z-w-z**4+3*z**3-3*z**2+z
print("\n=== H_C branch K check ===")
# zeta = root of z^3+2z^2-z-1 in (0,1)
zeta=sp.nsolve(z**3+2*z**2-z-1,z,0.8,prec=30)
print(f"  zeta={float(zeta):.10f}")
# sample z=0.9 in (zeta,1), find positive w root of HC, then u,v via E3,E2, compute K
for zv in [0.9,0.95]:
    HCz=HC.subs(z,zv)
    wrots=[complex(r) for r in sp.nroots(HCz)]
    for wr in wrots:
        if abs(wr.imag)<1e-8 and 0<wr.real<1:
            wr=wr.real
            # v=z(1-w) (does H_C branch also have v=z(1-w)? NO - that was H_B-specific)
            # For H_C, solve E3,E2 for u,v given w,z
            u,v=sp.symbols('u v')
            a3g=1-v+u*v;a5g=1-z+z*w-z*v*w+z*u*v*w
            E2=u*(1-z)-z*a5g*(1-v);E3=a3g*v-u*(1-w)
            sol=sp.nsolve([E2,E3],[u,v],[0.5,0.5],prec=30,maxsteps=200)
            uv,vv=float(sol[0]),float(sol[1])
            a5v=1-z+z*w-z*vv*w+z*uv*vv*w
            K=uv*vv*w*(z**3)*a5v**2/((1-vv)*(1-w)*(1-z)**3)
            pval=sp.nsolve((sp.sqrt if False else 1)*(1),z,0.5)
            print(f"  z={zv}: w={wr:.6f} u={uv:.6f} v={vv:.6f} K={K:.6e}  >1/128={1/128:.6e}? {K>1/128}")
            break
