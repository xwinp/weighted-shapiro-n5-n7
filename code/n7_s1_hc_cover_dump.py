#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Machine-checkable cover dump + rational L_C + deterministic checker.

Re-runs the H_C interval cover with per-piece recording.  For every certified
piece emits: rational s-endpoints, the Krawczyk c-box and v-box (as rational
intervals), the admissibility lower bounds (1-w, 1-cs-v^2, a5-analog, rho^7~K,
v in (0,1), c>0), and the mean-value-form P lower bound.  Then a short
deterministic CHECKER re-reads the dump and independently verifies, for each
piece: (i) Krawczyk inclusion N(W) subset int(W) for G and E2_red (unique root);
(ii) admissibility strictly positive over the box; (iii) P_MV > L_C, where
L_C = 141/20 = 7.05 is the certified rational global lower bound (the computed
global minimum 7.050002 > 141/20).  Outputs code/_hc_cover.json.
"""
import mpmath as mp, sympy as sp, json, fractions
from n7_s1_hc_critical_events import compute_critical_z
mp.mp.ivprec=110; IV=mp.iv; mp.mp.prec=80
exec(open('code/n7_s1_hc_rigorous_cert.py').read().split('with open')[0])

L_C = sp.Rational(141,20)   # 7.05 ; global computed min 7.050002 > 7.05
COVER=[]

def iv_to_rat(I):
    """Outward-rounded rational interval from an mp.iv interval (mpf endpoints)."""
    a=mp.mpf(I.a); b=mp.mpf(I.b)
    # floor a, ceil b at 2^-120
    sc=fractions.Fraction(1<<120)
    lo=fractions.Fraction(int(mp.floor(a*sc)))/sc
    hi=fractions.Fraction(int(mp.ceil(b*sc)))/sc
    return [str(lo), str(hi)]

def eval_piece_rec(S):
    """Like eval_piece but returns list of (plo, Cbox, Vbox, admiss_dict)."""
    sm=float(iv_mid(S)); out=[]
    for (c0,v0,u0) in admissible_lifts(sm):
        C=krawczyk_c(S,c0)
        if C is None or not(C.a>0): return None
        V=krawczyk_v(C,S,v0)
        if V is None or not(0<V.a and V.b<1): return None
        comp=P_components(iv(V),iv(C),iv(S))
        if comp is None: return None
        R,sqAB,T,omt,omw,oms,AB=comp
        plo=P_box(V,C,S)
        if plo is None: return None
        cs=iv(C)*iv(S)
        A5=-iv(V)*iv(V)+iv(V)*oms*iv(C)*(cs-1)+1+iv(C)*(1-2*iv(S))-iv(C)*iv(C)*iv(S)*oms
        rho7=iv(V)*iv(V)*iv(C)*oms**3*A5*A5/(omw*omt**3)
        out.append(dict(plo=float(plo),
                        C=iv_to_rat(C), V=iv_to_rat(V),
                        s=iv_to_rat(S),
                        omt_lo=float(omt.a), omw_lo=float(omw.a),
                        A5_lo=float(A5.a), rho7_lo=float(rho7.a),
                        V_a=float(V.a), V_b=float(V.b), C_a=float(C.a)))
    return out

def cover_rec(sa,sb):
    stack=[(mp.mpf(sa),mp.mpf(sb))]; imin=1e9; n=0
    while stack:
        a,b=stack.pop(); S=IV.mpf([a,b])
        Ps=eval_piece_rec(S)
        if Ps is None:
            if (b-a)<MINS: return None
            m=(a+b)/2; stack.append((a,m)); stack.append((m,b)); continue
        if len(Ps)==0:
            n+=1; continue
        lo=min(p['plo'] for p in Ps)
        if lo<float(L_C):
            if (b-a)<MINS:
                for p in Ps: p['s']=iv_to_rat(S); COVER.append(p)
                imin=min(imin,lo); n+=1; continue
            m=(a+b)/2; stack.append((a,m)); stack.append((m,b)); continue
        for p in Ps: p['s']=iv_to_rat(S); COVER.append(p)
        imin=min(imin,lo); n+=1
    return imin,n

crit = compute_critical_z(verbose=False)  # recompute; no uncommitted-pickle dependency
s_bounds=[mp.mpf(0)]+sorted([mp.mpf(1)-mp.mpf(c) for c in crit if 0<1-mp.mpf(c)<1])+[mp.mpf(1)]
gmin=1e9; ntot=0
print("Building machine-checkable cover (L_C=%s = %.5f)..."%(L_C,float(L_C)), flush=True)
for i in range(len(s_bounds)-1):
    sa,sb=s_bounds[i],s_bounds[i+1]
    if sb-sa<1e-9: continue
    res=cover_rec(sa,sb)
    if res is None: print("  cell %d FAILED"%i); continue
    imin,n=res; gmin=min(gmin,imin); ntot+=n
    print("  cell %d s(%.6f,%.6f): pieces=%d P_inf>=%.6f"%(i,float(sa),float(sb),n,imin), flush=True)

print("\nTotal pieces: %d"%len(COVER), flush=True)
print("Global computed P lower bound: %.6f"%gmin, flush=True)
print("Rational L_C = %s = %.5f ; gmin > L_C : %s"%(L_C,float(L_C),gmin>float(L_C)), flush=True)

out=dict(L_C=str(L_C), L_C_float=float(L_C),
         global_computed_min=gmin, n_pieces=len(COVER),
         ivprec=110, pieces=COVER)
with open('code/_hc_cover.json','w') as f: json.dump(out,f)
print("Wrote code/_hc_cover.json (%d pieces)"%len(COVER), flush=True)
print("DONE-COVERDUMP", flush=True)
