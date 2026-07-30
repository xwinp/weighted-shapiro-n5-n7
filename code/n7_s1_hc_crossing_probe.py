#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the H_C crossing polynomial  (C-7)^2 - 4 A B = 0  with the verified
closed form  P = C + 2 sqrt(A B), and probe elimination degrees.

Variables: v (algebraic via E2u_num=0, cubic in v), w (via H_C=0, cubic in w),
z (parameter), rho.  u = usol = v(1-v)/((1-w)-v^2)  (rational).
  A=(1-u)(1+rho)
  B=(1+rho) u v w z / (rho^5 (1-u)(1-v)(1-w)(1-z))
  C=rho(1+rho)[ (1-v)(1-u)/u + (1-w)(1-v)/v + (1-z)(1-w)/w + (1-z)/z ]
  closure: rho^7 = u v w z^3 a5^2 / ((1-v)(1-w)(1-z)^3),  a5=1-z+zw-zvw+zuwvw
Crossing poly (no sqrt, no rho1):  (C-7)^2 - 4 A B  (clear denoms).
"""
import sympy as sp
v,w,z,rho = sp.symbols('v w z rho', positive=True)
u = v*(1-v)/((1-w)-v**2)
a5 = 1 - z + z*w - z*v*w + z*u*v*w
HC = z*w**3 + w**2*z**3 - w**2*z + w*z**4 - 3*w*z**3 + 2*w*z**2 + w*z - w - z**4 + 3*z**3 - 3*z**2 + z
# E2u_num (v-polynomial): from E2 with u=usol
E2 = u*(1-z) - z*a5*(1-v)
E2t = sp.together(E2)
E2u_num = sp.expand(E2t.as_numer_denom()[0])
print("E2u_num: deg_v=%d deg_w=%d deg_z=%d"%(sp.Poly(E2u_num,v).degree(),sp.Poly(E2u_num,w).degree(),sp.Poly(E2u_num,z).degree()))

A = (1-u)*(1+rho)
B = (1+rho)*u*v*w*z/(rho**5*(1-u)*(1-v)*(1-w)*(1-z))
C = rho*(1+rho)*((1-v)*(1-u)/u + (1-w)*(1-v)/v + (1-z)*(1-w)/w + (1-z)/z)

# crossing = (C-7)^2 - 4 A B ; put over common denom, keep numerator
cross = sp.together((C-7)**2 - 4*A*B)
cross_num = sp.expand(cross.as_numer_denom()[0])
cross_den = sp.expand(cross.as_numer_denom()[1])
print("cross_num: deg_v=%d deg_w=%d deg_z=%d deg_rho=%d"%(
    sp.Poly(cross_num,v).degree(),sp.Poly(cross_num,w).degree(),
    sp.Poly(cross_num,z).degree(),sp.Poly(cross_num,rho).degree()))
# closure num
Kexpr = u*v*w*z**3*a5**2/((1-v)*(1-w)*(1-z)**3)
Kt = sp.together(Kexpr)
Knum = sp.expand(Kt.as_numer_denom()[0])
Kden = sp.expand(Kt.as_numer_denom()[1])
closure_num = sp.expand(rho**7*Kden - Knum)   # rho^7 * den - num = 0
print("closure_num: deg_v=%d deg_w=%d deg_z=%d deg_rho=%d"%(
    sp.Poly(closure_num,v).degree(),sp.Poly(closure_num,w).degree(),
    sp.Poly(closure_num,z).degree(),sp.Poly(closure_num,rho).degree()))
print("HC: deg_w=%d deg_z=%d"%(sp.Poly(HC,w).degree(),sp.Poly(HC,z).degree()))

# Save the key polynomials for the elimination script
import pickle
with open('code/_hc_polys.pickle','wb') as f:
    pickle.dump({'E2u_num':E2u_num,'cross_num':cross_num,'closure_num':closure_num,'HC':HC}, f)
print("saved polys. Now probe elimination sizes:")
# Eliminate v first: Res_v(E2u_num, cross_num), Res_v(E2u_num, closure_num)
dv_c = sp.Poly(cross_num,v).degree(); dv_cl = sp.Poly(closure_num,v).degree(); dv_E = sp.Poly(E2u_num,v).degree()
print("  Res_v(E2u,cross): Sylvester %d+%d=%d"%(dv_E,dv_c,dv_E+dv_c))
print("  Res_v(E2u,closure): Sylvester %d+%d=%d"%(dv_E,dv_cl,dv_E+dv_cl))
if dv_E+dv_c<=16 and dv_E+dv_cl<=16:
    print("  computing Res_v(E2u, closure_num)...")
    Rv_cl = sp.resultant(sp.Poly(E2u_num,v), sp.Poly(closure_num,v), v)
    Rv_cl = sp.expand(Rv_cl)
    print("  Rv_closure: deg_w=%d deg_z=%d deg_rho=%d"%(
        sp.Poly(Rv_cl,w).degree(), sp.Poly(Rv_cl,z).degree(), sp.Poly(Rv_cl,rho).degree()))
    with open('code/_hc_Rv_closure.pickle','wb') as f: pickle.dump(Rv_cl,f)
    print("  computing Res_v(E2u, cross_num)...")
    Rv_cr = sp.resultant(sp.Poly(E2u_num,v), sp.Poly(cross_num,v), v)
    Rv_cr = sp.expand(Rv_cr)
    print("  Rv_cross: deg_w=%d deg_z=%d deg_rho=%d"%(
        sp.Poly(Rv_cr,w).degree(), sp.Poly(Rv_cr,z).degree(), sp.Poly(Rv_cr,rho).degree()))
    with open('code/_hc_Rv_cross.pickle','wb') as f: pickle.dump(Rv_cr,f)
    drho_cl = sp.Poly(Rv_cl,rho).degree(); drho_cr = sp.Poly(Rv_cr,rho).degree()
    print("  Res_rho(Rv_closure,Rv_cross): Sylvester %d+%d=%d"%(drho_cl,drho_cr,drho_cl+drho_cr))
print("DONE-PROBE")
