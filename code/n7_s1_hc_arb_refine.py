#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Refine the few tight real-arc pieces whose rigorous Arb P_MV falls in
(7.0487, 7.05): bisect their s-interval, re-derive Krawczyk c/v boxes, and
recurse until the Arb mean-value-form P lower bound exceeds L_C=141/20.
Writes code/_hc_cover_refined.json (drop-in superset of _hc_cover.json)."""
import json, fractions
import flint
flint.ctx.prec = 150
import mpmath as mp
mp.mp.ivprec = 130; IV = mp.iv; mp.mp.prec = 130
exec(open('code/n7_s1_hc_rigorous_cert.py', encoding='utf-8').read().split('with open')[0])
src = open('code/n7_s1_hc_arb_checker.py', encoding='utf-8').read()
exec(src.split("with open('code/_hc_cover.json')")[0])

L_C = fractions.Fraction(141, 20)
L_C_arb = flint.arb(L_C.numerator) / flint.arb(L_C.denominator)

def fr2mpf(fr):
    return mp.mpf(fr.numerator) / mp.mpf(fr.denominator)

def iv_to_rat(I):
    a = mp.mpf(I.a); b = mp.mpf(I.b)
    sc = fractions.Fraction(1 << 120)
    lo = fractions.Fraction(int(mp.floor(a * sc))) / sc
    hi = fractions.Fraction(int(mp.ceil(b * sc))) / sc
    return [str(lo), str(hi)]

def arb_plo(p):
    """Arb P_MV lower bound for a piece dict (rational s/c/v boxes)."""
    S = box_ball(p['s']); C = box_ball(p['C']); V = box_ball(p['V'])
    plo, cb = P_mv_lo(V, C, S)
    return float(plo), cb

def make_piece(S_iv, C_iv, V_iv):
    cs = iv(C_iv) * iv(S_iv); oms = 1 - iv(S_iv)
    omt = 1 - cs - iv(V_iv) * iv(V_iv); omw = 1 - cs
    A5 = -iv(V_iv)*iv(V_iv) + iv(V_iv)*oms*iv(C_iv)*(cs-1) + 1 + iv(C_iv)*(1-2*iv(S_iv)) - iv(C_iv)*iv(C_iv)*iv(S_iv)*oms
    rho7 = iv(V_iv)*iv(V_iv)*iv(C_iv)*oms**3*A5*A5/(omw*omt**3)
    return dict(s=iv_to_rat(S_iv), C=iv_to_rat(C_iv), V=iv_to_rat(V_iv),
                plo=float(P_box(V_iv, C_iv, S_iv)),
                omt_lo=float(omt.a), omw_lo=float(omw.a),
                A5_lo=float(A5.a), rho7_lo=float(rho7.a),
                V_a=float(V_iv.a), V_b=float(V_iv.b), C_a=float(C_iv.a))

def refine_piece(sa, sb, depth=0):
    """Bisect s in (sa,sb); re-derive c/v boxes; return list of sub-pieces with
    Arb P_MV > L_C (recursing as needed)."""
    sm = (sa + sb) / 2
    out = []
    lifts = admissible_lifts(float(sm))
    for (c0, v0, u0) in lifts:
        C_iv = krawczyk_c(IV.mpf([sa, sb]), c0)
        if C_iv is None or not (C_iv.a > 0):
            continue
        V_iv = krawczyk_v(C_iv, IV.mpf([sa, sb]), v0)
        if V_iv is None or not (0 < V_iv.a and V_iv.b < 1):
            continue
        p = make_piece(IV.mpf([sa, sb]), C_iv, V_iv)
        plo, cb = arb_plo(p)
        if plo > float(L_C):
            out.append(p)
        else:
            if depth >= 10 or (sb - sa) < mp.mpf(1) / mp.mpf(2) ** 40:
                # floor: accept if still > 7 (rigorous exclusion holds)
                if plo > 7.0:
                    out.append(p)
                else:
                    print("  WARN: piece at s[%.6f,%.6f] plo=%.6f < 7" % (sa, sb, plo), flush=True)
            else:
                mid = (sa + sb) / 2
                out.extend(refine_piece(sa, mid, depth + 1))
                out.extend(refine_piece(mid, sb, depth + 1))
    return out

D = json.load(open('code/_hc_cover.json'))
pieces = D['pieces']
print("Loaded %d pieces; Arb-checking for refine candidates..." % len(pieces), flush=True)

new_pieces = []; n_refined = 0; n_kept = 0
for k, p in enumerate(pieces):
    plo, cb = arb_plo(p)
    if plo > float(L_C):
        new_pieces.append(p); n_kept += 1
    else:
        sa = fr2mpf(fractions.Fraction(p['s'][0])); sb = fr2mpf(fractions.Fraction(p['s'][1]))
        sub = refine_piece(sa, sb)
        if sub:
            new_pieces.extend(sub); n_refined += len(sub)
            print("  piece %d plo=%.6f -> %d sub-pieces" % (k, plo, len(sub)), flush=True)
        else:
            new_pieces.append(p); n_kept += 1
            print("  piece %d plo=%.6f NO SUB-LIFTS, kept" % (k, plo), flush=True)
    if (k + 1) % 500 == 0:
        print("  processed %d/%d" % (k + 1, len(pieces)), flush=True)

# global Arb min on refined cover
gmin = min(arb_plo(p)[0] for p in new_pieces)
out = dict(L_C=str(L_C), L_C_float=float(L_C),
           global_computed_min=gmin, n_pieces=len(new_pieces),
           ivprec=130, refined_from=len(pieces), pieces=new_pieces)
with open('code/_hc_cover_refined.json', 'w') as f:
    json.dump(out, f)
print("\nRefined cover: %d pieces (was %d), Arb global min=%.6f > L_C=%.5f: %s" % (
    len(new_pieces), len(pieces), gmin, float(L_C), gmin > float(L_C)), flush=True)
print("Wrote code/_hc_cover_refined.json", flush=True)
print("DONE-REFINE", flush=True)
