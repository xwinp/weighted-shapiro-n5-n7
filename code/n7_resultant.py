#!/usr/bin/env python3
"""
Verify GPT's n=7 hand-reduction and compute the EXACT resultant E(p)=Res_t(R,B).
R(p,t) = q^3 - p^3 t^5 - p^2 q t^8        (stationary curve, q=1-p)
B(p,t) = 5 p^2 t + 2 p q t^4 - 2 q^2 - 7 q p^2   (boundary P=7)
a_7,b_7 = p-coords of real solutions {R=0,B=0}.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 60

p, t = sp.symbols('p t')
q = 1 - p
R = q**3 - p**3*t**5 - p**2*q*t**8
B = 5*p**2*t + 2*p*q*t**4 - 2*q**2 - 7*q*p**2

# ---- 1. verify reduction against the high-precision nsolve boundary solution ----
# from phase3: at a_7, (a,b,c,d,e)=(x1,x2,x3,x4,x6) summing to 1:
a7_val = mp.mpf('0.21427352090984096558774231909001838797745454906426773562141322120851430')
b7_val = mp.mpf('0.32862767791659196669734085647711086500878730337285459186534900232409311442317')
# minimizer vars at a_7 (from phase3 nsolve, ~10 digits): re-derive at high precision via nsolve
import sympy as sp2
aa,bb,cc,dd,ee,lam,pv = sp2.symbols('aa bb cc dd ee lam pv', positive=True)
qv = 1 - pv
Pv = aa/(pv*bb+qv*cc)+bb/(pv*cc+qv*dd)+cc/(pv*dd)+dd/(qv*ee)+ee/(qv*aa)
gs=[sp2.diff(Pv,v) for v in (aa,bb,cc,dd,ee)]
eqs=[g-lam for g in gs]+[aa+bb+cc+dd+ee-1, Pv-7]
sol=sp2.nsolve(eqs,[aa,bb,cc,dd,ee,lam,pv],[0.17,0.23,0.06,0.31,0.23,6.0,0.22],prec=60,tol=mp.mpf('1e-50'),maxsteps=300)
a_n,b_n,c_n,d_n,e_n = [mp.mpf(sp2.N(x,55)) for x in sol[:5]]
p_n = mp.mpf(sp2.N(sol[6],55))
print(f"nsolve a_7 branch: p={mp.nstr(p_n,40)}")
print(f"  a,b,c,d,e = {mp.nstr(a_n,12)},{mp.nstr(b_n,12)},{mp.nstr(c_n,12)},{mp.nstr(d_n,12)},{mp.nstr(e_n,12)}")
# normalize a=1, t=e/a
t_val = e_n / a_n
d_pred = t_val**2
print(f"  t=e/a={mp.nstr(t_val,15)}  d={mp.nstr(d_n/a_n,15)}  t^2={mp.nstr(d_pred,15)}  (should match)")
# evaluate R,B at (p_n, t_val)
def Rf(pv,tv): qv=1-pv; return qv**3 - pv**3*tv**5 - pv**2*qv*tv**8
def Bf(pv,tv): qv=1-pv; return 5*pv**2*tv + 2*pv*qv*tv**4 - 2*qv**2 - 7*qv*pv**2
print(f"  R(p,t) = {mp.nstr(Rf(p_n,t_val),8)}  (should be ~0)")
print(f"  B(p,t) = {mp.nstr(Bf(p_n,t_val),8)}  (should be ~0)")

# ---- 2. exact resultant ----
print("\n=== computing E(p) = Res_t(R, B) ===")
E = sp.resultant(R, B, t)
E = sp.expand(E)
print(f"resultant degree in p: {sp.degree(E, p)}")
Ef = sp.factor(E)
print(f"factored resultant:\n{Ef}")

# collect p-only factors
facs = [f for f in sp.Mul.make_args(Ef) if f.free_symbols <= {p} or (hasattr(f,'is_Pow') and f.base.free_symbols<={p})]
# simpler: get factors as (poly,exp)
pf = sp.Poly(E, p)
print(f"\nE(p) as polynomial, degree {pf.degree()}")
# square-free factorization
sqf = sp.sqf_list(E, p)
print(f"square-free factorization: {[(sp.degree(f,p),e) for f,e in sqf[1]]}")

# find which factor has roots a_7, b_7
print("\n=== identifying factor containing a_7,b_7 ===")
for f, e in sqf[1]:
    fpoly = sp.Poly(f, p)
    deg = fpoly.degree()
    if deg == 0: continue
    # evaluate at a_7, b_7 numerically
    fa = complex(fpoly.eval(sp.Float(a7_val, 50)))
    fb = complex(fpoly.eval(sp.Float(b7_val, 50)))
    print(f"  deg {deg} factor: |f(a_7)|={abs(fa):.3e}  |f(b_7)|={abs(fb):.3e}")
    if abs(fa) < 1e-20 and abs(fb) < 1e-20:
        print(f"  *** FOUND: degree-{deg} factor contains a_7 and b_7 ***")
        print(f"  factor: {f}")
        with open("code/n7_minpoly_factor.txt","w") as fh:
            fh.write(str(f)+"\n")
