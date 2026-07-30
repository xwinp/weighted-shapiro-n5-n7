#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Explore the H_C branch structure for the rigorous P>7 certificate.

Goal: understand the algebraic structure so we can pick the safest rigorous route.
On H_C (asymmetric), the stationary point is determined by z (1 param):
  - H_C(w,z)=0  (cubic in w)
  - u = U_C(w,z), v = V_C(w,z)   (rational, from E2,E3 p-free beta-recurrence)
  - closure:  rho^7 = K_C(w,z)    (rational in w,z)
  - g1(rho1; w,z,rho)=0           (KKT for the height ratio rho1=b; polynomial)
  - P = P(w,z,rho,rho1)           (rational)
So P_C is algebraic in z (rho1 is an algebraic function of z).

This script:
 (1) build u(v), v(w,z) rational lifts and the closure K_C;
 (2) build g1 and P symbolically, report degrees in rho1;
 (3) eliminate rho1 from {g1=0, P-7=0} -> relation R(w,z,rho); report degrees;
 (4) numerically on the 8 H_C points: confirm rho1-root gives P, and that the
     univariate-in-z crossing poly (after full elimination) has no admissible root.
Aborts any Sylvester > 14 to avoid the memory-explosion lesson.
"""
import sympy as sp
import numpy as np
from scipy.optimize import brentq

n=7
z,w,u,v, rho, rho1, p = sp.symbols('z w u v rho rho1 p', positive=True)
q = 1 - p
# ---- H_C branch curve (cubic in w) ----
HC = z*w**3 + w**2*z**3 - w**2*z + w*z**4 - 3*w*z**3 + 2*w*z**2 + w*z - w \
     - z**4 + 3*z**3 - 3*z**2 + z
# ---- p-free beta-recurrence lifts (from n7_s1_hc_definitive_trace) ----
a3 = 1 - v + u*v
a5 = 1 - z + z*w - z*v*w + z*u*v*w
E3 = a3*v - u*(1-w)        # => u = v(1-v)/((1-w)-v^2)
E2 = u*(1-z) - z*a5*(1-v)
usol = v*(1-v)/((1-w)-v**2)
E2u = sp.together(E2.subs(u, usol))
E2u_num = sp.expand(E2u.as_numer_denom()[0])   # poly in v (given z,w)
print("(1) E2(u(v)) numerator: deg_v =", sp.Poly(E2u_num, v).degree(),
      " deg_z =", sp.Poly(E2u_num, z).degree(), " deg_w =", sp.Poly(E2u_num, w).degree())

# closure K_C = u v w z^3 a5^2 / ((1-v)(1-w)(1-z)^3) = rho^7
Kexpr = usol*v*w*z**3*a5.subs(u,usol)**2 / ((1-v)*(1-w)*(1-z)**3)
Knum = sp.expand(sp.together(Kexpr).as_numer_denom()[0])
Kden = sp.expand(sp.together(Kexpr).as_numer_denom()[1])
print("    closure num: total deg =", sp.total_degree(Knum, w, z, v),
      " den total deg =", sp.total_degree(Kden, w, z, v))

# ---- P in (rho1, rho, w,z,u,v) with rho2..rho5 from beta ----
# beta_i = q*rho_{i+1}/(p+q rho_{i+1}) = w? Need mapping. From trace: w=beta3, z=beta4,
#   u=beta1?, v=beta2?  rho_{i+1}=(p/q)*beta_i/(1-beta_i)=(rho)*beta_i/(1-beta_i) since p/q=1/rho.
# beta = (u,v,w,z) = (beta1,beta2,beta3,beta4) -> rho2=rho*u/(1-u), rho3=rho*v/(1-v),
#   rho4=rho*w/(1-w), rho5=rho*z/(1-z).
r2 = rho*u/(1-u); r3 = rho*v/(1-v); r4 = rho*w/(1-w); r5 = rho*z/(1-z)
# P = sum_{i=1}^4 1/(rho_i (p+q rho_{i+1})) + 1/(p rho5) + prod_rho/q
# with rho=(rho1, r2, r3, r4, r5).  p = 1/(1+rho), q = rho/(1+rho).
pp = 1/(1+rho); qq = rho/(1+rho)
prod_rho = rho1*r2*r3*r4*r5
# term i: 1/(rho_i (p+q rho_{i+1}));  p+q*rho_{i+1} = p + qq*r_{i+1}
def term(ri, rnext):
    return 1/(ri*(pp + qq*rnext))
Pexpr = term(rho1, r2) + term(r2, r3) + term(r3, r4) + term(r4, r5) + 1/(pp*r5) + prod_rho/qq
Pexpr = sp.together(Pexpr)
# g1 = rho1 * dP/drho1
g1 = sp.together(rho1 * sp.diff(Pexpr, rho1))
g1_num = sp.expand(g1.as_numer_denom()[0])
g1_den = sp.expand(g1.as_numer_denom()[1])
print("(2) g1 = rho1 dP/drho1: num deg_rho1 =", sp.Poly(g1_num, rho1).degree(),
      " den deg_rho1 =", sp.Poly(g1_den, rho1).degree())
# P-7 numerator
Pm7 = sp.together(Pexpr - 7)
Pm7_num = sp.expand(Pm7.as_numer_denom()[0])
Pm7_den = sp.expand(Pm7.as_numer_denom()[1])
print("    P-7: num deg_rho1 =", sp.Poly(Pm7_num, rho1).degree(),
      " den deg_rho1 =", sp.Poly(Pm7_den, rho1).degree())

# substitute usol into g1_num, Pm7_num (u rational in v,w,z). Keep v symbolic for now;
# but v is itself algebraic in (w,z) via E2u_num=0. We'll handle v later.
print("\n(3) degree check after substituting u=usol (rational in v)...")
g1_sub = sp.together(g1_num.subs(u, usol))
Pm7_sub = sp.together(Pm7_num.subs(u, usol))
g1s_num = sp.expand(g1_sub.as_numer_denom()[0])
Pm7s_num = sp.expand(Pm7_sub.as_numer_denom()[0])
print("    g1 (u subst) num: deg_rho1=%d deg_v=%d deg_w=%d deg_z=%d deg_rho=%d"%(
    sp.Poly(g1s_num, rho1).degree(), sp.Poly(g1s_num, v).degree(),
    sp.Poly(g1s_num, w).degree(), sp.Poly(g1s_num, z).degree(), sp.Poly(g1s_num, rho).degree()))
print("    P-7 (u subst) num: deg_rho1=%d deg_v=%d deg_w=%d deg_z=%d deg_rho=%d"%(
    sp.Poly(Pm7s_num, rho1).degree(), sp.Poly(Pm7s_num, v).degree(),
    sp.Poly(Pm7s_num, w).degree(), sp.Poly(Pm7s_num, z).degree(), sp.Poly(Pm7s_num, rho).degree()))
# Eliminate rho1: Res_{rho1}(g1s_num, Pm7s_num).  Size = deg_rho1(g1)+deg_rho1(Pm7).
dg1 = sp.Poly(g1s_num, rho1).degree(); dP = sp.Poly(Pm7s_num, rho1).degree()
print("    Sylvester size for Res_rho1 = %d + %d = %d"%(dg1, dP, dg1+dP))
if dg1+dP <= 14:
    print("    => SAFE. Computing Res_rho1...")
    R1 = sp.resultant(sp.Poly(g1s_num, rho1), sp.Poly(Pm7s_num, rho1), rho1)
    R1 = sp.expand(R1)
    print("    Res_rho1 computed. deg_v=%d deg_w=%d deg_z=%d deg_rho=%d"%(
        sp.Poly(R1, v).degree(), sp.Poly(R1, w).degree(),
        sp.Poly(R1, z).degree(), sp.Poly(R1, rho).degree()))
    sp.save('code/_hc_R1_rho1elim.pickle', R1)
    print("    saved R1 to code/_hc_R1_rho1elim.pickle")
else:
    print("    => TOO LARGE, abort (memory risk).")
print("DONE-EXPLORE")
