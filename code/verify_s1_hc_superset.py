#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EXACT load-bearing certificate: H_C KKT superset containment (rev 12).

Closes the gap flagged in the round-17 re-audit: the comment in
verify_s1_elimination.py claiming "F1 is the E2_red v-relation up to a
prefactor" is FALSE (F1 == E2u_num is False; see verify_s1_elim_step3.py).
This script provides the correct load-bearing exact certificate, with NO
floating point and NO dependence on finite-sample residuals.

GOAL.  Show that every admissible S1 KKT point on the H_C branch satisfies
E2_red = 0, hence
    V_KKT^adm  ∩  V(H_C)  ⊆  V(H_C, E2_red, E3, closure K, g1).

PROOF (all steps certified exactly below).

  (E3)   L_u,L_v eliminate to E3 = u v^2 + u w - u - v^2 + v = 0  =>
         u = v(1-v)/((1-w)-v^2).                                    [verify_s1_elimination.py]

  (F1,F2) L_u,L_w and L_u,L_z eliminate, then u -> usol, give the two
          rho-free v-relations F1=0, F2=0.  A KKT point satisfies both.

  (H_C)  Res_v(F1,F2) = w^2 z^2 (w-1)^2 H_B H_C  (EXACT, verify_s1_elimination.py).
         On the admissible set (w,z in (0,1)) the prefactors are nonzero, so a
         KKT point has H_B H_C = 0; on the H_C branch, H_C = 0.

  (KEY)  The KKT v-root on H_C=0 is the COMMON v-root of F1 and F2.  For two
         quadratics in v this common root is given exactly by the LINEAR
         SUBRESULTANT
             S_1(F1,F2) = sigma1 * v + sigma0,    r = -sigma0/sigma1,
         valid wherever sigma1 != 0 (i.e. F1,F2 share exactly one v-root).
         With  g = sigma1 / w  (so sigma1 = w * g), we certify the EXACT identity
             g(w,z)^3 * E2u_num(r)  =  z^4 (w-1)^2 (z-1)^2 * B * A * H_C,      (*)
         where  E2u_num = (v-1) * E2_red  (v=1 is the inadmissible boundary
         beta_2=1),  A = -w^2 z + w z^2 - w - z^2 + z,  B = w z^2 - w z - w - z^2 + z.
         (The denominator of E2u_num(r) alone is h^3, h = sigma1/(w z) = g/z, NOT
         sigma1^3; clearing it with g^3 = z^3 h^3 gives (*) with the z^4 factor.
         An earlier rev-12 draft wrongly wrote sigma1^3 * E2u_num(r) = z(w-1)^2...
         -- that conflated the *numerator* of E2u_num(r) with sigma1^3 * E2u_num(r),
         which differs by the factor w^3 z^3.  Verified by checking the WHOLE
         rational identity together(g^3 * E2u_num(r) - RHS) == 0, not just a
         numerator.)  On H_C=0 the RHS of (*) is 0, so (g != 0) => E2u_num(r) = 0;
         admissibility (v != 1) => E2_red = 0.

  (sig1) sigma1 = w * g(w,z).  We certify
             Res_w(g, H_C) = -z^11 (z-1)^6,                                    (**)
         which is nonzero for z in (0,1); hence g and H_C share no common real
         w-root there, so g != 0 at every w-root of H_C with z in (0,1).  This
         both validates the use of r = -sigma0/sigma1 (sigma1 = w*g, w != 0 on the
         admissible set) AND is the nonvanishing factor that turns (*) on H_C=0
         into E2u_num(r) = 0.

  (K,g1) closure K and g1 come from g_1 (Euler, lambda=0): substituting the KKT
         u,v,rho_1 into g_1 gives the closure rho^7 = K(u,v,w,z) and g_1 = 0;
         these are recorded in verify_s1_elimination.py and are not re-derived here.

All assertions are exact (rational polynomial arithmetic; sp.div remainder zero,
sp.expand(...) == 0, sp.resultant).  Exit 0 iff every assertion holds.
"""
import sympy as sp

r1, u, v, w, z, rho = sp.symbols('r1 u v w z rho', positive=True)

# ---- S1 setup (identical to verify_s1_elimination.py) -----------------------
A = (1 - u) * (1 + rho)
B = (1 + rho) * u * v * w * z / (rho**5 * (1 - u) * (1 - v) * (1 - w) * (1 - z))
S = ((1 - v) * (1 - u) / u + (1 - w) * (1 - v) / v
     + (1 - z) * (1 - w) / w + (1 - z) / z)
C = rho * (1 + rho) * S
dA = {s: sp.together(sp.diff(A, s)) for s in (u, v, w, z)}
dB = {s: sp.together(sp.diff(B, s)) for s in (u, v, w, z)}
dC = {s: sp.together(sp.diff(C, s)) for s in (u, v, w, z)}


def L(j):
    return sp.expand(sp.together(B * r1 * dC[j] + B * dA[j] + A * dB[j]).as_numer_denom()[0])


Lu, Lv, Lw, Lz = L(u), L(v), L(w), L(z)


def ab(Lx):
    c1 = sp.expand(Lx.coeff(r1, 1)); c0 = sp.expand(Lx - c1 * r1); return c1, c0


au, bu = ab(Lu); av, bv = ab(Lv); aw, bw = ab(Lw); az, bz = ab(Lz)


def elim(ai, bi, aj, bj):
    return sp.expand(sp.together(ai * bj - aj * bi).as_numer_denom()[0])


E_uv = elim(au, bu, av, bv); E_uw = elim(au, bu, aw, bw); E_uz = elim(au, bu, az, bz)

# (E3)
E3 = u * v**2 + u * w - u - v**2 + v
assert sp.rem(sp.expand(E_uv), sp.expand(E3), u, v, w, z) == 0 or any(
    sp.factor(e[0]) == sp.factor(E3) for e in sp.factor_list(E_uv)[1]), "E3 not in E_uv"
usol = sp.together(v * (1 - v) / ((1 - w) - v**2))


def subs_u(E):
    return sp.expand(sp.together(E.subs(u, usol)).as_numer_denom()[0])


Euw = subs_u(E_uw); Euz = subs_u(E_uz)
G = sp.gcd(Euw, Euz)
F1 = sp.cancel(Euw / G)        # rho-free v-relation (core)
F2 = sp.cancel(Euz / G)        # rho-free v-relation (core)

# (H_C) completeness identity
H_B = z * w**2 + (1 - z**2) * w + z**2 - z
H_C = (z * w**3 + w**2 * z**3 - w**2 * z + w * z**4 - 3 * w * z**3
       + 2 * w * z**2 + w * z - w - z**4 + 3 * z**3 - 3 * z**2 + z)
R_F1F2 = sp.resultant(F1, F2, v)
assert sp.expand(R_F1F2 - w**2 * z**2 * (w - 1)**2 * H_B * H_C) == 0, "completeness FAILED"
print("(H_C) Res_v(F1,F2) = w^2 z^2 (w-1)^2 H_B H_C           [OK]", flush=True)

# ---- E2_red --------------------------------------------------------------
a3 = 1 - v + u * v
a5 = 1 - z + z * w * a3
E2 = sp.expand(u * (1 - z) - z * a5 * (1 - v))
E2u = sp.together(E2.subs(u, usol))
E2u_num = sp.expand(E2u.as_numer_denom()[0])
# E2u_num = (v-1) * E2_red  ;  v=1 is the inadmissible boundary beta_2=1
q_vm1, r_vm1 = sp.div(sp.Poly(sp.expand(E2u_num), v), sp.Poly(v - 1, v))
assert r_vm1.is_zero, "E2u_num not divisible by (v-1)"
E2_red = sp.expand(q_vm1.as_expr())
assert sp.expand(E2u_num - (v - 1) * E2_red) == 0
print("(E2)  E2u_num = (v-1) * E2_red                           [OK]", flush=True)

# ---- KEY: linear subresultant of F1, F2 and the exact identity ------------
subs_list = sp.subresultants(F1, F2, v)
S1 = None
for s in subs_list:
    if sp.Poly(sp.expand(s), v).degree() == 1:
        S1 = sp.expand(s); break
assert S1 is not None, "no linear subresultant"
sigma1 = sp.expand(S1.coeff(v, 1))
sigma0 = sp.expand(S1.coeff(v, 0))
assert sp.expand(sigma1.subs(w, 0)) == 0, "sigma1 not divisible by w"
g = sp.cancel(sigma1 / w)              # sigma1 = w * g
h = sp.cancel(sigma1 / (w * z))        # sigma1 = w * z * h,  g = z * h
# root r = -sigma0/sigma1 ; evaluate E2u_num(r).  The denominator of E2u_num(r)
# alone is h^3 (NOT sigma1^3); clearing it with g^3 = z^3 h^3 yields the exact
# polynomial identity (*) below.  We verify the WHOLE rational identity, not just
# a numerator (the rev-12 bug was to check only the numerator and then label it
# sigma1^3 * E2u_num(r), which is wrong by the factor w^3 z^3).
A_f = -w**2 * z + w * z**2 - w - z**2 + z
B_f = w * z**2 - w * z - w - z**2 + z
target = sp.expand(z * (w - 1)**2 * (z - 1)**2 * B_f * A_f * H_C)
E2_at_r = sp.together(E2u_num.subs(v, -sigma0 / sigma1))
# (KEY) g^3 * E2u_num(r) = z^4 (w-1)^2 (z-1)^2 B A H_C   [the TRUE identity (*)]
correct_rhs = sp.expand(z**4 * (w - 1)**2 * (z - 1)**2 * B_f * A_f * H_C)
assert sp.expand(sp.together(g**3 * E2_at_r - correct_rhs)) == 0, "KEY identity (*) FAILED"
# defensive: confirm the clearing-denominator structure explicitly
_num_r, _den_r = sp.together(E2_at_r).as_numer_denom()
assert sp.expand(_den_r - sp.expand(h**3)) == 0, "den(E2u_num(r)) != h^3"
assert sp.expand(_num_r - target) == 0, "num(E2u_num(r)) != target"
print("(KEY) g^3 * E2u_num(r) = z^4 (w-1)^2 (z-1)^2 B A H_C        [OK]", flush=True)
print("      (denominator of E2u_num(r) is h^3 = g^3/z^3; cleared by g^3)", flush=True)
print("      r = -sigma0/sigma1  (linear subresultant root; KKT v on H_C=0)", flush=True)

# ---- g != 0 on admissible H_C branch (turns (*) into E2u_num(r)=0) ---------
# g = sigma1/w.  Res_w(g, H_C) is a polynomial in z whose REAL zeros are exactly
# the points where g and H_C share a real w-root; we certify its only real zeros
# in [0,1] are the inadmissible endpoints z=0, z=1 (beta_4 = 0 or 1), so g != 0
# at every real w-root of H_C with z in (0,1).  Scaling-invariant: we only
# require that every irreducible factor of Res_w(g,H_C) over Q is z or z-1.
Res_gHC = sp.expand(sp.resultant(g, H_C, w))
fl = sp.factor_list(Res_gHC)            # (coeff, [(base, mult), ...])
bad = [(sp.factor(base), mult) for base, mult in fl[1]
       if sp.expand(base - z) != 0 and sp.expand(base - (z - 1)) != 0]
assert fl[0] != 0 and not bad, "g!=0 cert (**) FAILED: Res_w(g,H_C)=" + str(sp.factor(Res_gHC))
print("(sig1) sigma1 = w*g; Res_w(g,H_C) factored:", sp.factor(Res_gHC), flush=True)
print("       only real zeros in [0,1] are z=0, z=1 (inadmissible endpoints)", flush=True)
print("       => g != 0 on admissible H_C branch (z in (0,1))", flush=True)
print("       => on H_C=0: (*) gives E2u_num(r) = 0  [OK]", flush=True)

# ---- conclusion ----------------------------------------------------------
print("\nCONTAINMENT (EXACT): on the admissible H_C branch, g != 0, so the identity", flush=True)
print("  g^3 * E2u_num(r) = z^4 (w-1)^2 (z-1)^2 B A H_C = 0  =>  E2u_num(r) = 0;", flush=True)
print("  admissibility (v != 1) => E2_red = 0.  Hence", flush=True)
print("  V_KKT^adm  cap  V(H_C)  subset  V(H_C, E2_red, E3, K, g1).", flush=True)
print("DONE-HC-SUPERSET", flush=True)
