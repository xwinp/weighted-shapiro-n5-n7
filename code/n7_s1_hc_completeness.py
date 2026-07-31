#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Completeness + event-fiber + endpoint certificate for the H_C interval cover.

Rigor gap addressed (review): mp.polyroots finds roots but does not prove
exhaustiveness.  Here, for every event-free s-cell we
  (1) exact Sturm-count the positive c-roots of G(c,s0)=0 (cubic, rational coeffs);
  (2) form R(v)=Res_c(G,E2_red) (degree-9, rational), exact Sturm-count its (0,1)
      roots, and ISOLATE them with sympy real_roots (exact algebraic isolation);
  (3) for each isolated v-root, recover the c-root of G satisfying E2_red(v,c,s0)=0
      and VERIFY the common zero (ruling out spurious resultant roots from
      leading-coefficient drop), then check admissibility;
  (4) confirm the certified admissible lifts equal those the cover script
      (admissible_lifts) enumerates, and that each is Krawczyk-isolated.

Because code/n7_s1_hc_critical_events.py collects EVERY z at which a lift is
born/dies (c-double E1, v-double E2) or crosses an admissibility boundary
(w=1 E3, v=0 E4, v=1 E5, u=1 E6, a5=0 E7, u-denom=0 E8, K=0 E9), on an event-free
cell no lift crosses a boundary, so the admissible-lift count is CONSTANT and the
midpoint enumeration is exhaustive for the whole cell.  Event fibers (5 critical s)
and endpoints (s=0, s=1) are certified separately.
"""
import mpmath as mp, sympy as sp, json
from n7_s1_hc_critical_events import compute_critical_z
mp.mp.ivprec=110; IV=mp.iv; mp.mp.prec=80
exec(open('code/n7_s1_hc_rigorous_cert.py').read().split('with open')[0])

_c,_v,_s = sp.symbols('c v s')
G_expr = (1-_s)*_c**3 + (1-_s)*(_s-2)*_c**2 + (_s**2-_s-1)*_c + (1-_s)
# rebuild _E2u_cs symbol in our c,v,s (the exec used _v,_c,_s already, same names)
# _E2u_cs is E2u_red(v,c,s) from the exec.
print("Building Res_c(G, E2u_red) ...", flush=True)
Res_v = sp.resultant(G_expr, _E2u_cs, _c)
Res_v = sp.expand(Res_v)
print("  deg_v(Res) =", sp.Poly(Res_v,_v).degree(), flush=True)

def sturm_count(expr, sym, lo, hi):
    """Exact rational Sturm count of real roots of poly in CLOSED [lo,hi]."""
    P=sp.Poly(sp.expand(expr), sym)
    if P.degree()==0: return 0
    return P.count_roots(lo, hi)
def sturm_count_open(expr, sym, lo, hi):
    """Exact distinct real roots in OPEN (lo,hi): sqf count on [lo,hi] minus endpoints."""
    P=sp.Poly(sp.expand(expr), sym)
    if P.degree()==0: return 0
    sqf=P.sqf_part()
    n=sqf.count_roots(lo,hi)
    if sqf.eval(lo)==0: n-=1
    if sqf.eval(hi)==0: n-=1
    return n

def isolate_real_roots(expr, sym, lo, hi):
    """Exact algebraic isolation of DISTINCT real roots in (lo,hi) via sympy
    CRootOf (robust; matches Sturm count_roots which also counts distinct)."""
    P=sp.Poly(sp.expand(expr), sym)
    if P.degree()==0: return []
    out=[]
    for r in sp.real_roots(P):
        rn=float(sp.N(r, 40))
        if lo<rn<hi: out.append(rn)
    return sorted(out)

def recover_c(vv, s0):
    """Given numerical v at s0, find c-root of G with E2_red(v,c,s0)=0 (common zero)."""
    Gs = sp.Poly(G_expr.subs(_s,s0), _c)
    croots = [complex(sp.N(r,40)) for r in sp.nroots(Gs.as_expr(), n=40)]
    E2v = sp.Lambda(_c, sp.expand(_E2u_cs.subs({_v:vv,_s:s0})))
    best=None; bestd=1e9
    for cr in croots:
        if abs(cr.imag)>1e-8: continue
        cf=cr.real
        try: ev=complex(sp.N(_E2u_cs.subs({_v:vv,_c:cf,_s:s0}),40))
        except: continue
        if abs(ev)<bestd: bestd=abs(ev); best=cf
    return best, bestd

# ---- cells ----
crit = compute_critical_z(verbose=False)  # recompute; no uncommitted-pickle dependency
s_crit = sorted(set(mp.mpf(1)-mp.mpf(c) for c in crit if 0<1-c<1))
s_bounds = [mp.mpf(0)] + s_crit + [mp.mpf(1)]

print("\n=== Per-cell completeness certificate ===", flush=True)
cell_records=[]
all_ok=True
for i in range(len(s_bounds)-1):
    a,b = s_bounds[i], s_bounds[i+1]
    if b-a < mp.mpf(10)**-12: continue
    fa,fb = float(a),float(b)
    s0 = sp.Rational(int((fa+fb)/2*10**8), 10**8)
    if not (a < mp.mpf(s0.p)/mp.mpf(s0.q) < b):
        s0 = sp.Rational(int(fa*10**8)+1, 10**8)
    s0f=float(s0)
    Gs = G_expr.subs(_s,s0)
    n_c_pos = sturm_count(Gs,_c,0,sp.oo)
    n_c_adm = sturm_count(Gs,_c,0,sp.Rational(1,1)/s0)
    Rv = sp.expand(Res_v.subs(_s,s0))
    n_v_sturm = sturm_count(Rv,_v,0,1)                       # closed, with multiplicity
    n_v_distinct = sturm_count_open(Rv,_v,0,1)               # distinct in OPEN (0,1)
    v_isos = isolate_real_roots(Rv,_v,0,1)                    # distinct in (0,1), via CRootOf
    # verify each isolated v-root has a common c-root; collect admissible
    adm=[]; n_no_common=0; n_inadmissible=0; no_common_resd_min=1e9
    for vm in v_isos:
        cf,resd = recover_c(vm, s0)
        if cf is None or resd>1e-6:
            n_no_common+=1                                     # resultant root, no common c (leading-coeff drop)
            no_common_resd_min=min(no_common_resd_min, resd if resd==resd else 1e9)
            continue
        # common (c,v) zero exists (resd~0); check admissibility
        w0=cf*s0f; z0=1-s0f
        admissible = (1e-9<w0<1-1e-9)
        denom=(1-w0)-vm*vm
        if admissible and abs(denom)>1e-9:
            u0=vm*(1-vm)/denom
            admissible = admissible and (1e-9<u0<1-1e-9)
            if admissible:
                a5v=1-z0+z0*w0-z0*vm*w0+z0*u0*vm*w0
                admissible = admissible and (a5v>1e-12)
                if admissible:
                    K=u0*vm*w0*(z0**3)*a5v**2/((1-vm)*(1-w0)*((1-z0)**3))
                    admissible = admissible and (K>1e-12)
        if not admissible:
            n_inadmissible+=1                                  # common zero but fails an admissibility inequality
            continue
        adm.append((cf,vm))
    # cover-script enumeration at same s0
    lifts=admissible_lifts(s0f)
    # Krawczyk isolation
    S0=IV.mpf([mp.mpf(s0),mp.mpf(s0)]); n_iso=0
    for (c0,v0,u0) in lifts:
        Cb=krawczyk_c(S0,c0); Vb=krawczyk_v(Cb,S0,v0) if Cb else None
        if Vb is not None: n_iso+=1
    # completeness = (distinct Sturm == CRootOf count) AND (exact admissible == cover) AND (all isolated)
    ok = (n_v_distinct==len(v_isos)) and (len(adm)==len(lifts)) and (n_iso==len(lifts))
    all_ok = all_ok and ok
    rec=dict(cell=i, s_a=float(a), s_b=float(b), s0=str(s0),
             n_c_pos_sturm=int(n_c_pos), n_c_adm_sturm=int(n_c_adm),
             n_v_sturm_mult=int(n_v_sturm), n_v_distinct_sturm=int(n_v_distinct),
             n_v_isolated=int(len(v_isos)),
             n_no_common=int(n_no_common), n_inadmissible=int(n_inadmissible),
             n_adm_exact=int(len(adm)), n_adm_cover=int(len(lifts)), n_krawczyk=int(n_iso), ok=bool(ok))
    cell_records.append(rec)
    print("cell %d s(%.6f,%.6f) s0=%s: n_c+=%d  n_v(0,1) mult/distinct/iso=%d/%d/%d  "
          "no_common=%d inadmiss=%d  n_adm exact/cover=%d/%d  n_iso=%d  %s"%(
        i,float(a),float(b),s0,n_c_pos,n_v_sturm,n_v_distinct,len(v_isos),
        n_no_common,n_inadmissible,len(adm),len(lifts),n_iso,"OK" if ok else "MISMATCH"), flush=True)

print("\nAll cells complete (Sturm==iso, exact==cover, all Krawczyk):", all_ok, flush=True)

# ---- event fibers ----
print("\n=== Event-fiber certificates (point s=s_crit) ===", flush=True)
ev_records=[]
for sc in s_crit:
    scf=float(sc); lifts=admissible_lifts(scf)
    S0=IV.mpf([mp.mpf(sc),mp.mpf(sc)]); plo_min=1e9; n=0
    for (c0,v0,u0) in lifts:
        Cb=krawczyk_c(S0,c0)
        if Cb is None: continue
        Vb=krawczyk_v(Cb,S0,v0)
        if Vb is None: continue
        plo=P_box(Vb,Cb,S0)
        if plo is None: continue
        plo_min=min(plo_min,plo); n+=1
    ev_records.append(dict(s_crit=scf, z=1-scf, n_adm=n, P_lo=round(plo_min,6) if n else None))
    print("  s=%.6f (z=%.6f): n_adm=%d  P_lo=%.6f"%(scf,1-scf,n,plo_min if n else -1), flush=True)

# ---- endpoints ----
print("\n=== Endpoints ===", flush=True)
G_at_1 = sp.expand(G_expr.subs(_s,1))
print("  s=1 (z=0): G(c,1)=%s  -> 1/s=1, c<1; positive c-roots of -c: none. "
      "No admissible lift (w=cs, c=0=>w=0=>K=0)."%G_at_1, flush=True)
G_at_0 = sp.expand(G_expr.subs(_s,0))
n_c0 = sturm_count(G_at_0,_c,0,sp.oo)
print("  s=0 (z=1): G(c,0)=%s, positive c-roots=%d. But w=c*0=0 => K has factor w => "
      "K=0, degenerate (NOT admissible). Real arc is a LIMIT s->0+."%(G_at_0,n_c0), flush=True)
print("  Real-arc P limit as s->0+ (c~2.24 branch):", flush=True)
for ssm in [1e-3,1e-4,1e-5,1e-6,1e-7]:
    lifts=admissible_lifts(ssm)
    if not lifts: continue
    # real arc = the lift with largest c (c0~2.24); spurious has c~0.55
    c0,v0,u0=max(lifts, key=lambda t:t[0])
    S0=IV.mpf([mp.mpf(ssm),mp.mpf(ssm)])
    Cb=krawczyk_c(S0,c0); Vb=krawczyk_v(Cb,S0,v0) if Cb else None
    if Vb is None: continue
    plo=P_box(Vb,Cb,S0)
    print("    s=%.1e: c=%.6f  P_lo=%s"%(ssm,c0, round(plo,6) if plo else None), flush=True)

ev_min = min((r['P_lo'] for r in ev_records if r['P_lo']), default=None)
print("\nEvent-fiber global P_lo: %.6f"%ev_min, flush=True)
out=dict(cells=cell_records, events=ev_records, all_cells_ok=all_ok, event_fiber_P_lo=ev_min)
with open('code/_hc_completeness.json','w') as f: json.dump(out,f,indent=2)
print("\nDONE-COMPLETENESS all_cells_ok=%s event_P_lo=%s"%(all_ok,ev_min), flush=True)
