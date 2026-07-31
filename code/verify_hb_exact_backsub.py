#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Forward-construction RECORD on the H_B branch  (NON-LOAD-BEARING).

This script is NOT a load-bearing verifier.  The load-bearing forward
containment is `verify_s1_elimination.py`, which proves the exact resultant
identity

      Res_v(F_1, F_2)  =  w^2 z^2 (w-1)^2  H_B(w,z) H_C(w,z)        (E)

as an exact polynomial identity (`expand(R - target) == 0`), so every genuine
S_1 KKT point projects into {H_B = 0} u {H_C = 0} (containment, not equality).

What this script reconfirms (exact, `sp.div` remainder zero):
  (B)  H_B  exactly divides  Res_v(F_1, F_2)   (univariate division in w,
      coefficients in Q[z]; remainder == 0), and the cofactor is the
      polynomial  w^2 z^2 (w-1)^2 H_C.  This is the forward-containment
      identity (E), rechecked here by genuine polynomial division.

What this script does NOT claim (and previously claimed in error):
  (A)  It does NOT claim that the w- or z-KKT numerator L_w, L_z is a
       polynomial multiple of F_1, F_2.  In fact F_1 does NOT divide
       L_w|_{u=usol} as a v-polynomial (`sp.div` remainder != 0): L_w
       depends on rho, while F_1 is rho-free, so the relation between
       L_w = 0 and F_1 = 0 at a KKT point is mediated by the resultant (E)
       and the r_1-consistency, NOT by a direct divisibility L_w = prefactor*F_1.
       The earlier `cancel(Lw_n/F1)` check was tautological (q := Lw_n/F_1
       => Lw_n - q*F_1 == 0 by definition) and masked this; it is removed.
  (rev) It makes NO reverse-soundness claim (H_B = 0 => a genuine KKT lift
       exists).  A vanishing specialised resultant can also reflect a
       projective root at infinity, a complex common root, or a real root
       outside (0,1).  The theorem analyses the whole positive H_B-candidate
       graph, which is stronger than the true KKT subset.

No floating-point residual is invoked.  See verify_s1_elimination.py for the
load-bearing forward-containment certificate.
"""
import sympy as sp

v, w, z = sp.symbols('v w z')

# branch relations (match verify_s1_elimination.py exactly)
F1 = v**2 * w + v * w**3 - v * w**2 + v * w * z - v * w - v * z + v - w**3 + 2 * w**2 - w
F2 = v**2 * w * z + v * w * z**3 - v * w * z**2 - v * w * z + v * w - v * z**3 + v * z**2 + w**2 * z - w * z
H_B = z * w**2 + (1 - z**2) * w + z**2 - z

print("=== (B) forward-containment resultant identity (load-bearing: verify_s1_elimination.py) ===", flush=True)
R = sp.resultant(F1, F2, v)
# genuine polynomial division: H_B | Res_v(F1,F2) as w-polynomials (coeffs in Q[z])
q_HB, r_HB = sp.div(sp.Poly(sp.expand(R), w), sp.Poly(sp.expand(H_B), w))
HB_divides = r_HB.is_zero
print("  Res_v(F1,F2) / H_B  (sp.div, univariate in w) remainder == 0:", HB_divides, flush=True)

# (B) reconfirmation by GENUINE polynomial division (this replaces the tautological
# cancel(L_w/F_1) check).  H_B exactly divides R = Res_v(F_1,F_2) as a w-polynomial
# (coeffs in Q[z]); the quotient Q = R/H_B is the polynomial w^2 (w-1)^2 H_C.
#
# Scaling note (why no z^2 here): verify_hb uses the *unscaled* F_1 (no -z
# prefactor), whereas the load-bearing verify_s1_elimination.py uses F_1 scaled
# by -z, so its resultant carries an extra z^2 and reads w^2 z^2 (w-1)^2 H_B H_C.
# Both express the SAME containment H_B | Res_v(F_1,F_2); they differ only by the
# z^2 scaling inherited from F_1.  The cofactor check below uses verify_hb's own
# (unscaled) resultant, so its cofactor is w^2 (w-1)^2 H_C.
H_C = (w**3*z + w**2*z**3 - w**2*z + w*z**4 - 3*w*z**3 + 2*w*z**2
       + w*z - w - z**4 + 3*z**3 - 3*z**2 + z)
q_HB_expr = sp.expand(q_HB.as_expr())          # exact quotient Q = R / H_B (in w)
cofactor_identity = sp.expand(q_HB_expr - w**2 * (w - 1)**2 * H_C) == 0
print("  H_B | Res_v(F1,F2)  (sp.div, univariate in w, rem == 0):", HB_divides, flush=True)
print("  Res_v(F1,F2) / H_B  ==  w^2 (w-1)^2 H_C  (exact quotient, expand == 0):",
      cofactor_identity, flush=True)
print("  [canonical w^2 z^2 (w-1)^2 H_B H_C: load-bearing verify_s1_elimination.py]",
      flush=True)

print(flush=True)
print("=== (A) NOT claimed: L_w, L_z are NOT polynomial multiples of F_1, F_2 ===", flush=True)
# demonstrate (for honesty) that F_1 does NOT divide L_w|_{u=usol}: build L_w and check.
r1, u, rho = sp.symbols('r1 u rho')
A = (1 - u) * (1 + rho)
B = (1 + rho) * u * v * w * z / (rho**5 * (1 - u) * (1 - v) * (1 - w) * (1 - z))
S = (1 - v) * (1 - u) / u + (1 - w) * (1 - v) / v + (1 - z) * (1 - w) / w + (1 - z) / z
C = rho * (1 + rho) * S
dA = {s: sp.together(sp.diff(A, s)) for s in (u, v, w, z)}
dB = {s: sp.together(sp.diff(B, s)) for s in (u, v, w, z)}
dC = {s: sp.together(sp.diff(C, s)) for s in (u, v, w, z)}
def Lnum(j):
    expr = B * r1 * dC[j] + B * dA[j] + A * dB[j]
    return sp.expand(sp.together(expr).as_numer_denom()[0])
Lw = Lnum(w)
E3 = u * v**2 + u * w - u - v**2 + v
usol = sp.together(v * (1 - v) / ((1 - w) - v**2))
Lw_n = sp.expand(sp.together(Lw.subs(u, usol)).as_numer_denom()[0])
_, r_Lw_F1 = sp.div(sp.Poly(Lw_n, v), sp.Poly(F1, v))
print("  F_1 divides L_w|_{u=usol} (sp.div rem==0):", r_Lw_F1.is_zero,
      "  <- FALSE; L_w is NOT a polynomial multiple of F_1.", flush=True)
print("  (The KKT<->F_1=0 relation is mediated by the resultant (E), not by L_w = prefactor*F_1.)", flush=True)

print(flush=True)
print("=== (C) closure + g1 fix rho, r1 > 0 at a genuine KKT point on H_B (qualitative) ===", flush=True)
print("  closure: rho^7 = w z^5 / ((1-z) D^3), D=1-z^2+w z^2; on H_B (w,z in (0,1)) RHS > 0", flush=True)
print("  g1: r1^2 = A/B; A=(1-u)(1+rho)>0, B>0 on admissible set => r1 = +sqrt(A/B) > 0", flush=True)
print("  (positivity certified by interval arithmetic in the cover checker)", flush=True)

print(flush=True)
ok = bool(HB_divides and cofactor_identity)
if ok:
    print("FORWARD-CONTAINMENT IDENTITY RECONFIRMED (non-load-bearing record):", flush=True)
    print("  (B) H_B | Res_v(F1,F2) (sp.div rem 0); cofactor R/H_B = w^2 (w-1)^2 H_C (exact).", flush=True)
    print("  Canonical w^2 z^2 (w-1)^2 H_B H_C (z^2 from F_1 scaling): verify_s1_elimination.py.", flush=True)
    print("  (A) NOT claimed: L_w is not a polynomial multiple of F_1 (no tautological cancel check).", flush=True)
    print("  No reverse soundness (H_B => KKT lift) is claimed; no float residual is invoked.", flush=True)
    print("DONE-HB-FORWARD-RECORD", flush=True)
else:
    print("FORWARD-CONTAINMENT IDENTITY FAILED — check above.", flush=True)
    raise SystemExit(1)
