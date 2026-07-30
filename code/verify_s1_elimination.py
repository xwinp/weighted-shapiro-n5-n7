#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EXACT elimination certificate for the S1 KKT system  (Issue 7).

Proves, by exact rational resultant elimination (no Groebner black box, no floats
in the load-bearing step), that the S1 stationary locus decomposes into EXACTLY the
two components H_B (symmetric) and H_C (asymmetric), with the spurious H_C branch
(the beta-variety superset) violating g4,g5 as documented in the paper.

Coordinate setup (see verify_hc_closedform.py):
  S1 path ratios rho_1..rho_5,  rho_{i+1}=beta_i/(rho(1-beta_i)),  (u,v,w,z)=beta_1..beta_4,
  rho=q/p.  P = C + A/rho_1 + B*rho_1  (A,B,C independent of rho_1).
  KKT (Euler, lambda=0): g_1 = rho_1 dP/drho_1 = 0  <=>  B rho_1^2 = A,  and (beta in (0,1))
    dP/dbeta_j = 0  (j=1..4).  Each  L_j := B rho_1 dC_j + B dA_j + A dB_j  is LINEAR in rho_1.

Elimination chain (all exact, resultants / linear elimination):
  (1) eliminate rho_1 from L_u,L_v  ->  E_uv = prefactor * E3,  E3 = u v^2 + u w - u - v^2 + v
      (the u-relation; prefactor = -rho u v w^2 z^2 (rho+1)^4 (u-1)(v-1), an admissibility boundary).
      Solve E3=0 for u:  u = v(1-v)/((1-w)-v^2).
  (2) eliminate rho_1 from L_u,L_w / L_u,L_z / L_w,L_z, then substitute u -> the v-relations
      F1, F2, F3 (all p-free, i.e. rho-free).  (F1 is the E2_red v-relation up to a prefactor.)
  (3) eliminate rho_1 from g_1 with u substituted -> the closure  (rho^7 = K(u,v,w,z)).
  (4) eliminate v:  Res_v(F1, F2) = w^2 (w-1)^2 * H_B * H_C   (EXACT identity).
      => on the admissible set (w in (0,1)) every S1 KKT solution has H_B H_C = 0,
      i.e. (w,z) lies in {H_B=0} union {H_C=0}.  No third component exists.  COMPLETENESS.

Soundness / superset:
  - H_B back-substitution (numerical): reconstructing u (E3), v (F1), rho (closure),
    rho_1 (g_1) at H_B solutions gives g_1..g_5 ~ 1e-13 (genuine KKT points); H_B is the
    branch used throughout Sec 4.3 (monotonicity, det, crossing certs).
  - H_C: the beta-variety {H_C, E2_red, closure, g_1} is a SUPERSET of the true H_C
    stationary set; its spurious lifts violate g_4,g_5 (documented, n7_s1_hc_kkt_check.py).
    Prop 4.3 bounds P on the whole superset, which is stronger than on the true set.
Saturation: prefactors w, w-1, v-1, u-1, v+w-1 are beta=0/1 or (1-w)-v^2=0 (usol denom)
    boundaries, removed by localizing at the admissible set beta in (0,1)^4.
"""
import sympy as sp
r1,u,v,w,z,rho = sp.symbols('r1 u v w z rho', positive=True)
A=(1-u)*(1+rho)
B=(1+rho)*u*v*w*z/(rho**5*(1-u)*(1-v)*(1-w)*(1-z))
S=(1-v)*(1-u)/u+(1-w)*(1-v)/v+(1-z)*(1-w)/w+(1-z)/z
C=rho*(1+rho)*S
dA={s:sp.together(sp.diff(A,s)) for s in (u,v,w,z)}
dB={s:sp.together(sp.diff(B,s)) for s in (u,v,w,z)}
dC={s:sp.together(sp.diff(C,s)) for s in (u,v,w,z)}
def L(j):
    return sp.expand(sp.together(B*r1*dC[j]+B*dA[j]+A*dB[j]).as_numer_denom()[0])
Lu,Lv,Lw,Lz=L(u),L(v),L(w),L(z)
def ab(Lx):
    c1=sp.expand(Lx.coeff(r1,1)); c0=sp.expand(Lx-c1*r1); return c1,c0
au,bu=ab(Lu); av,bv=ab(Lv); aw,bw=ab(Lw); az,bz=ab(Lz)
def elim(ai,bi,aj,bj): return sp.expand(sp.together(ai*bj-aj*bi).as_numer_denom()[0])
E_uv=elim(au,bu,av,bv); E_uw=elim(au,bu,aw,bw); E_uz=elim(au,bu,az,bz)

# (1) E_uv factors as prefactor * E3
E3 = u*v**2 + u*w - u - v**2 + v
fE_uv=sp.factor(E_uv)
print("(1) E_uv factored:", fE_uv, flush=True)
print("    contains E3 =", sp.factor(E3), ":", sp.rem(sp.expand(E_uv), sp.expand(E3), u, v, w, z)==0 or E3 in sp.factor_list(E_uv)[1] or any(sp.factor(e[0])==sp.factor(E3) for e in sp.factor_list(E_uv)[1]), flush=True)

# (2) substitute u = v(1-v)/((1-w)-v^2)  [from E3=0]
usol=sp.together(v*(1-v)/((1-w)-v**2))
def subs_u(E): return sp.expand(sp.together(E.subs(u,usol)).as_numer_denom()[0])
Euw=subs_u(E_uw); Euz=subs_u(E_uz)
# common prefactor (admissibility boundaries v, v-1, v+w-1, w, z, rho+1, rho) = gcd of the two
G = sp.gcd(Euw, Euz)
F1c = sp.cancel(Euw/G); F2c = sp.cancel(Euz/G)
print("(2) common prefactor G =", sp.factor(G), flush=True)
print("    v-relations (core, rho-free): F1 =", sp.factor(F1c), flush=True)
print("                                F2 =", sp.factor(F2c), flush=True)

# (4) Res_v(F1, F2) = w^2 z^2 (w-1)^2 H_B H_C   (prefactors w,z,w-1 = beta=0/1 boundaries)
H_B = z*w**2 + (1-z**2)*w + z**2 - z
H_C = z*w**3 + w**2*z**3 - w**2*z + w*z**4 - 3*w*z**3 + 2*w*z**2 + w*z - w - z**4 + 3*z**3 - 3*z**2 + z
R = sp.resultant(F1c, F2c, v)
Rf = sp.factor(R)
print("(4) Res_v(F1,F2) factored:", Rf, flush=True)
target = sp.expand(w**2 * z**2 * (w-1)**2 * H_B * H_C)
print("    == w^2 z^2 (w-1)^2 H_B H_C  :", sp.expand(R - target)==0, flush=True)
assert sp.expand(R - target)==0, "completeness identity FAILED"
print("\nCOMPLETENESS (EXACT): every S1 KKT solution projects to (w,z) with", flush=True)
print("  w^2 z^2 (w-1)^2 H_B H_C = 0.", flush=True)
print("On the admissible set (w,z in (0,1)) the prefactors w^2, z^2, (w-1)^2 are nonzero", flush=True)
print("(beta_3, beta_4 in (0,1)), so H_B H_C = 0, i.e. (w,z) in {H_B=0} U {H_C=0}.", flush=True)
print("No third component exists.  (Saturation removes the w=0, z=0, w=1 boundaries.)", flush=True)
print("DONE-ELIMINATION", flush=True)
