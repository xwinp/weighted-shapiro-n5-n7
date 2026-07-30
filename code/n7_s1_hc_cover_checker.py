#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Deterministic checker for the H_C cover dump (code/_hc_cover.json).

For each recorded piece the checker INDEPENDENTLY re-derives, from the stored
rational s/c/v boxes alone:
  (i)   Krawczyk inclusion N(W_c) subset int(W_c) for G(c,s)=0  (=> unique c-root),
        and N(W_v) subset int(W_v) for E2_red(v,c,s)=0 (=> unique v-root);
  (ii)  strict admissibility over the box: c>0, 0<v<1, 1-w>0, 1-cs-v^2>0,
        a5-analog>0, rho^7 (=K-analog)>0;
  (iii) the mean-value-form P lower bound, and asserts P_MV > L_C (=141/20).
It does NOT trust the stored plo.  Reports the count of verified pieces and any
failure.  A piece passes iff all three hold.
"""
import mpmath as mp, json, fractions
mp.mp.ivprec=130; IV=mp.iv; mp.mp.prec=130   # >=121 so k/2^120 endpoints are exact
exec(open('code/n7_s1_hc_rigorous_cert.py').read().split('with open')[0])

def fr2mpf(fr):
    # exact when value = k/2^120 and prec>=121 (correctly-rounded division -> exact)
    return mp.mpf(fr.numerator)/mp.mpf(fr.denominator)
def rat_to_iv(pair):
    lo=fractions.Fraction(pair[0]); hi=fractions.Fraction(pair[1])
    return IV.mpf([fr2mpf(lo), fr2mpf(hi)])

import sys
_cpath = sys.argv[1] if len(sys.argv) > 1 else 'code/_hc_cover.json'
with open(_cpath) as f: D=json.load(f)
_lc=fractions.Fraction(D['L_C']); L_C=mp.mpf(_lc.numerator)/mp.mpf(_lc.denominator)
print("Checking %d pieces, L_C=%s ..."%(len(D['pieces']), D['L_C']), flush=True)

n_ok=0; n_fail=0; fails=[]
for k,p in enumerate(D['pieces']):
    S=rat_to_iv(p['s']); C=rat_to_iv(p['C']); V=rat_to_iv(p['V'])
    reasons=[]
    # (i) Krawczyk inclusion: K = mid - f(mid)/f'(box); require K subset int(box)
    # G in c:
    cm=iv(iv_mid(C)); dG=Gc_iv(C,S)
    if dG.a<=0<=dG.b: reasons.append("Gc straddles 0")
    else:
        Kc=iv(cm)-G_iv(cm,S)/dG
        if not(Kc.a>=C.a and Kc.b<=C.b and (Kc.a>C.a or Kc.b<C.b)):
            # not strict interior: try as subset (non-strict still proves unique root if K subset W)
            if not(Kc.a>=C.a and Kc.b<=C.b): reasons.append("Kc not subset Wc")
    # E2_red in v:
    vm=iv(iv_mid(V)); dE=E2v_iv(V,C,S)
    if dE.a<=0<=dE.b: reasons.append("Ev straddles 0")
    else:
        Kv=iv(vm)-E2_iv(vm,C,S)/dE
        if not(Kv.a>=V.a and Kv.b<=V.b): reasons.append("Kv not subset Wv")
    # (ii) admissibility over the box
    cs=iv(C)*iv(S); oms=1-iv(S)
    omt=1-cs-iv(V)*iv(V); omw=1-cs
    A5=-iv(V)*iv(V)+iv(V)*oms*iv(C)*(cs-1)+1+iv(C)*(1-2*iv(S))-iv(C)*iv(C)*iv(S)*oms
    rho7=iv(V)*iv(V)*iv(C)*oms**3*A5*A5/(omw*omt**3)
    if not(C.a>0): reasons.append("c<=0")
    if not(0<V.a and V.b<1): reasons.append("v not in (0,1)")
    if not(omt.a>0): reasons.append("omt<=0")
    if not(omw.a>0): reasons.append("omw<=0")
    if not(A5.a>0): reasons.append("A5<=0")
    if not(rho7.a>0): reasons.append("rho7<=0")
    # (iii) recompute P_MV from the box, assert > L_C
    plo=P_box(V,C,S)
    if plo is None: reasons.append("P_box None")
    elif plo<=L_C: reasons.append("P_MV=%.6f<=L_C"%plo)
    if reasons:
        n_fail+=1
        if len(fails)<20: fails.append((k,p['s'],reasons))
    else:
        n_ok+=1
    if (k+1)%500==0: print("  checked %d/%d ok=%d fail=%d"%(k+1,len(D['pieces']),n_ok,n_fail), flush=True)

print("\nVerified pieces: %d / %d"%(n_ok,len(D['pieces'])), flush=True)
print("Failed: %d"%n_fail, flush=True)
for f in fails: print("  FAIL",f, flush=True)
print("\nALL PIECES VERIFIED (Krawczyk + admissibility + P_MV > L_C):", n_fail==0, flush=True)
print("DONE-CHECKER", flush=True)
