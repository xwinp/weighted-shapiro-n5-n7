#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stepwise resultant elimination for S1 KKT  (Issue 7, step 3: v-elim -> H_B*H_C)."""
import sympy as sp
v,w,z,rho = sp.symbols('v w z rho', positive=True)
# F1 = E2_red (the v-relation after u-elim), F2,F3 the other two v-relations (from step2)
F1 = v**2*w + v*w**3 - v*w**2 + v*w*z - v*w - v*z + v - w**3 + 2*w**2 - w
F2 = v**2*w*z + v*w*z**3 - v*w*z**2 - v*w*z + v*w - v*z**3 + v*z**2 + w**2*z - w*z
F3 = -v*w**3*z + v*w**2*z + v*w*z**3 - 2*v*w*z**2 + v*w - v*z**3 + 2*v*z**2 - v*z + w**3*z - w**2*z

# Cross-check F1 == E2u_num from n7_s1_hc_kkt_check.py (up to a prefactor):
uS,vS,wS,zS=sp.symbols('u v w z')
a3=1-vS+uS*vS; a5=1-zS+zS*wS-zS*vS*wS+zS*uS*vS*wS
E2=uS*(1-zS)-zS*a5*(1-vS)
usol=sp.together(vS*(1-vS)/((1-wS)-vS**2))
E2u_num=sp.expand(sp.together(E2.subs(uS,usol)).as_numer_denom()[0])
print("F1 == E2u_num (E2_red):", sp.expand(F1 - E2u_num.subs({vS:v,wS:w,zS:z}))==0, flush=True)
# also check F2, F3 are multiples/consequences of F1
print("F2 in (v,w,z), degree in v:", sp.Poly(F2,v).degree(), flush=True)
print("F3 in (v,w,z), degree in v:", sp.Poly(F3,v).degree(), flush=True)

# Eliminate v:  Res_v(F1, F2) and Res_v(F1, F3)  -> (w,z) relation = H_B * H_C (times prefactors)
print("\nComputing Res_v(F1,F2)...", flush=True)
R12 = sp.resultant(F1, F2, v)
print("  factored:", flush=True); print("  ", sp.factor(R12), flush=True)
print("\nComputing Res_v(F1,F3)...", flush=True)
R13 = sp.resultant(F1, F3, v)
print("  factored:", flush=True); print("  ", sp.factor(R13), flush=True)

# H_B and H_C (paper definitions)
H_B = z*w**2 + (1-z**2)*w + z**2 - z
H_C = z*w**3 + w**2*z**3 - w**2*z + w*z**4 - 3*w*z**3 + 2*w*z**2 + w*z - w - z**4 + 3*z**3 - 3*z**2 + z
print("\nH_B*H_C factored:", sp.factor(sp.expand(H_B*H_C)), flush=True)
# Check H_B*H_C divides the v-resultant (after dropping prefactors)
q12,_=sp.div(sp.expand(R12), sp.expand(H_B*H_C), w, z)
print("\nDONE-STEP3", flush=True)
