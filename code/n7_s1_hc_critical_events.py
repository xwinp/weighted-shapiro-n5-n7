#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Enumerate the critical-z events of the admissible H_C lift set, to partition
(0,1) into z-intervals on each of which the admissible (w,v) lifts form fixed
continuous arcs.  All event resultants are SMALL (sizes <= ~8), so memory-safe.

Lift set: z in (0,1); w in (0,1) root of H_C(w,z)=0 (cubic in w);
v in (0,1) root of E2u_num(v,w,z)=0 (cubic in v); u=usol=v(1-v)/((1-w)-v^2) in (0,1);
a5>0; K_C>0.  Events (z-values where the admissible lift set changes):
  E1 w-double : Res_w(H_C, dH_C/dw)
  E2 v-double : Res_w(H_C, Res_v(E2u, dE2u/dv))
  E3 w=1      : H_C(1,z)
  E4 v=0      : Res_w(H_C, E2u_num(0,w,z))
  E5 v=1      : Res_w(H_C, E2u_num(1,w,z))
  E6 u=1      : Res_w(H_C, Res_v(E2u, num(usol-1)))    [u=0 subsumed by v=0,1]
  E7 a5=0     : Res_w(H_C, Res_v(E2u, num(a5)))         [a5 with u=usol]
  E8 u-denom=0: Res_w(H_C, Res_v(E2u, (1-w)-v^2))
  E9 K_C=0    : Res_w(H_C, Res_v(E2u, num(Kexpr)))
Collect all (0,1)-roots -> sorted critical z list.  (z=0,1 added as boundaries.)

The computation is deterministic (exact symbolic resultants + real_roots) and
is exposed via ``compute_critical_z`` so downstream scripts
(cover_dump / completeness / rigorous_cert) recompute it directly instead of
reading an uncommitted pickle — there is no implicit temp-file dependency.
"""
import sympy as sp

v,w,z = sp.symbols('v w z', positive=True)
u = v*(1-v)/((1-w)-v**2)
a5 = 1 - z + z*w - z*v*w + z*u*v*w
HC = z*w**3 + w**2*z**3 - w**2*z + w*z**4 - 3*w*z**3 + 2*w*z**2 + w*z - w - z**4 + 3*z**3 - 3*z**2 + z
E2 = u*(1-z) - z*a5*(1-v)
E2u_num = sp.expand(sp.together(E2).as_numer_denom()[0])
Kexpr = u*v*w*z**3*a5**2/((1-v)*(1-w)*(1-z)**3)
Knum = sp.expand(sp.together(Kexpr).as_numer_denom()[0])

def roots_in_01(poly, sym):
    poly = sp.Poly(sp.expand(poly), sym)
    if poly.degree()==0: return []
    return [sp.N(r,14) for r in sp.real_roots(poly) if 0<r<1]

def resw_z(poly_in_v_w_z, label, verbose=True):
    """Eliminate v then w; return (0,1)-roots in z."""
    pv = sp.Poly(sp.expand(poly_in_v_w_z), v)
    # resultant wrt v with E2u_num
    Rv = sp.expand(sp.resultant(sp.Poly(E2u_num,v), pv, v))
    Rw = sp.expand(sp.resultant(sp.Poly(HC,w), sp.Poly(Rv,w), w))
    rts = roots_in_01(Rw, z)
    if verbose:
        print("  %-12s deg_z(R)=%d  (0,1)-roots: %s"%(label, sp.Poly(Rw,z).degree() if Rw!=0 else 0,
            ["%.6f"%float(r) for r in rts]))
    return rts

def compute_critical_z(verbose=True):
    """Return the sorted list of critical z-values in (0,1) at which the
    admissible H_C lift set changes (E1..E9 event roots). Deterministic; no
    external file dependency. ``verbose`` controls the per-event diagnostic
    prints (callers that only need the values pass verbose=False)."""
    all_roots=set()
    # E1: w-double
    R1 = sp.expand(sp.resultant(sp.Poly(HC,w), sp.Poly(sp.diff(HC,w),w), w))
    r1 = roots_in_01(R1,z)
    if verbose: print("E1 w-double: %s"%["%.6f"%float(r) for r in r1])
    all_roots|=set(r1)
    # E2: v-double
    dE2u = sp.diff(E2u_num,v)
    Rv2 = sp.expand(sp.resultant(sp.Poly(E2u_num,v), sp.Poly(dE2u,v), v))
    R2 = sp.expand(sp.resultant(sp.Poly(HC,w), sp.Poly(Rv2,w), w))
    r2 = roots_in_01(R2,z)
    if verbose: print("E2 v-double: %s"%["%.6f"%float(r) for r in r2])
    all_roots|=set(r2)
    # E3: w=1
    R3 = sp.expand(HC.subs(w,1))
    r3 = roots_in_01(R3,z)
    if verbose: print("E3 w=1: %s"%["%.6f"%float(r) for r in r3])
    all_roots|=set(r3)
    # E4: v=0
    r4 = resw_z(E2u_num.subs(v,0), "E4 v=0", verbose); all_roots|=set(r4)
    # E5: v=1
    r5 = resw_z(E2u_num.subs(v,1), "E5 v=1", verbose); all_roots|=set(r5)
    # E6: u=1  -> usol-1=0 -> v(1-v)=((1-w)-v^2) -> v - v^2 -1 +w +v^2=0 -> v+w-1=0. So u=1 <=> v+w=1.
    #     Substitute v=1-w into E2u_num, eliminate w with HC.
    E6cond = sp.expand(E2u_num.subs(v, 1-w))
    R6 = sp.expand(sp.resultant(sp.Poly(HC,w), sp.Poly(E6cond,w), w))
    r6 = roots_in_01(R6,z)
    if verbose: print("E6 u=1(v+w=1): %s"%["%.6f"%float(r) for r in r6])
    all_roots|=set(r6)
    # E7: a5=0 (numerator, with u=usol). a5_num = numerator of a5 (u=usol).
    a5_num = sp.expand(sp.together(a5).as_numer_denom()[0])
    r7 = resw_z(a5_num, "E7 a5=0", verbose); all_roots|=set(r7)
    # E8: u-denom (1-w)-v^2=0
    r8 = resw_z((1-w)-v**2, "E8 u-denom", verbose); all_roots|=set(r8)
    # E9: K_C=0
    r9 = resw_z(Knum, "E9 K=0", verbose); all_roots|=set(r9)

    crit = sorted(float(r) for r in all_roots)
    if verbose:
        print("\nAll critical z in (0,1), sorted (%d):"%len(crit))
        print(["%.6f"%c for c in crit])
    return crit

if __name__ == "__main__":
    crit = compute_critical_z(verbose=True)
    print("DONE-EVENTS")
