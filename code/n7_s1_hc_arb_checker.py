#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Arb (python-flint) independent re-verification of the H_C cover.

Re-reads code/_hc_cover.json.  Per piece, with rigorous Arb BALL arithmetic
(rigorous for +,-,*,/, sqrt, 7th-root — unlike mpmath.iv's heuristic non-basic
ops), independently re-derives:
  (i)   iterative Krawczyk on G(c,s)=0 and E2_red(v,c,s)=0: contract and test
        K(W) subset int(W)  (=> unique root in the box);
  (ii)  strict admissibility (lower endpoint > 0);
  (iii) the MEAN-VALUE-FORM P lower bound (the form used by the certificate,
        not loose direct eval): P(box) subset P(mid)+sum partial_i*(box_i-mid),
        with partials and P(mid) evaluated by Arb; assert P_lo > L_C=141/20.
"""
import json, fractions
import flint
flint.ctx.prec = 150
exec(open('code/n7_s1_hc_rigorous_cert.py').read().split('with open')[0])

def ar(x):
    if isinstance(x, str):
        fr = fractions.Fraction(x)
        return flint.arb(fr.numerator) / flint.arb(fr.denominator)
    return flint.arb(x)
def box_ball(pair):
    return ar(pair[0]).union(ar(pair[1]))
def midpoint(b):
    return (ar(b.lower()) + ar(b.upper())) / 2     # near-point ball
def intersect(a, b):
    lo = a.lower(); hi = a.upper()
    glo = b.lower(); ghi = b.upper()
    nlo = lo if lo > glo else glo
    nhi = hi if hi < ghi else ghi
    if nlo > nhi:
        return None
    return ar(nlo).union(ar(nhi))
def contains_int(W, K):
    # K subset interior(W)  <=>  K.lo>W.lo and K.hi<W.hi
    return K.lower() > W.lower() and K.upper() < W.upper()

def arb_eval(expr, subs):
    if expr.is_Number:
        return flint.arb(int(expr)) if expr.is_Integer else flint.arb(str(expr))
    if expr.is_Symbol:
        return subs[expr]
    if expr.is_Pow:
        b = arb_eval(expr.base, subs); e = int(expr.exp)
        r = b
        for _ in range(abs(e) - 1):
            r = r * b
        return r if e > 0 else (flint.arb(1) / r)
    if expr.is_Mul:
        r = arb_eval(expr.args[0], subs)
        for a in expr.args[1:]:
            r = r * arb_eval(a, subs)
        return r
    if expr.is_Add:
        r = arb_eval(expr.args[0], subs)
        for a in expr.args[1:]:
            r = r + arb_eval(a, subs)
        return r
    raise ValueError("unhandled %r" % type(expr))

G_sym = (1-_s)*_c**3 + (1-_s)*(_s-2)*_c**2 + (_s**2-_s-1)*_c + (1-_s)
Gc_sym = sp.diff(G_sym, _c)
# lo/hi interval evaluators (tight; ball arithmetic over-widens (1-s) and Gc,
# breaking Krawczyk contraction on the spurious branch).
def G_lh(c, s):
    C = pf(c); S = pf(s); one = pc(1); two = pc(2)
    oms = psub1(s)                                   # 1-s >0
    s2sms1 = sneg(padd(one, pmul(S, oms)))           # s^2-s-1 = -(1+s(1-s))  (factored, tight)
    c2 = ppow(C, 2); c3 = ppow(C, 3)
    s_m_2 = ssub(S, two)                             # s-2 (signed)
    inner = sadd(sadd(c3, smul(s_m_2, c2)), one)     # c^3+(s-2)c^2+1
    return sadd(smul(oms, inner), smul(s2sms1, C))
def Gc_lh(c, s):
    C = pf(c); S = pf(s); one = pc(1); two = pc(2)
    oms = psub1(s)
    s2sms1 = sneg(padd(one, pmul(S, oms)))
    c2 = ppow(C, 2)
    s_m_2 = ssub(S, two)
    inner_g = sadd(smul(pc(3), c2), smul(smul(two, s_m_2), C))  # 3c^2+2(s-2)c
    return sadd(smul(oms, inner_g), s2sms1)
def G_ball(c, s):
    return to_ball(G_lh(c, s))
def Gc_ball(c, s):
    return to_ball(Gc_lh(c, s))
def E2_ball(v, c, s):
    r = arb_eval(_e[0], {_c: c, _s: s})
    for k in range(1, len(_e)):
        r = r + arb_eval(_e[k], {_c: c, _s: s}) * (v ** k)
    return r
def E2v_ball(v, c, s):
    r = arb_eval(_e[1], {_c: c, _s: s})
    for k in range(2, len(_e)):
        r = r + flint.arb(k) * arb_eval(_e[k], {_c: c, _s: s}) * (v ** (k - 1))
    return r

# ---- lo/hi interval arithmetic over Arb (tight, like mpmath.iv; rigorous for
#      all ops incl. sqrt / nth-root, unlike mpmath.iv's heuristic non-basic ops).
# An interval is a (lo, hi) pair of point arb balls.  Positive operands use
# monotone endpoint ops (exact range, full precision); signed use 4-product.
def pf(b):           # positive interval from a ball
    return (ar(b.lower()), ar(b.upper()))
def pc(x):           # point constant
    a = ar(x); return (a, a)
def padd(a, b):      # positive + positive
    return (a[0] + b[0], a[1] + b[1])
def pmul(a, b):      # positive * positive  -> [a0*b0, a1*b1]
    return (a[0] * b[0], a[1] * b[1])
def pdiv(a, b):      # positive / positive  -> [a0/b1, a1/b0]
    return (a[0] / b[1], a[1] / b[0])
def ppow(a, n):      # positive ** n (n>=1) -> [a0^n, a1^n]
    return (a[0] ** n, a[1] ** n)
def psqrt(a):        # positive sqrt -> [sqrt(a0), sqrt(a1)]
    return (a[0].sqrt(), a[1].sqrt())
def proot7(a):       # positive 7th root -> [a0^(1/7), a1^(1/7)]
    return (a[0].root(7), a[1].root(7))
def psub1(b):        # 1 - positive ball b  (b<1) -> [1-b.hi, 1-b.lo]
    return (flint.arb(1) - ar(b.upper()), flint.arb(1) - ar(b.lower()))
# signed ops
def sadd(a, b):
    return (a[0] + b[0], a[1] + b[1])
def ssub(a, b):
    return (a[0] - b[1], a[1] - b[0])
def sneg(a):
    return (-a[1], -a[0])
def smul(a, b):
    ps = [a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]]
    lo = min(p.lower() for p in ps); hi = max(p.upper() for p in ps)
    return (ar(lo), ar(hi))
def sdiv(a, b):   # b not containing 0
    ps = [a[0]/b[0], a[0]/b[1], a[1]/b[0], a[1]/b[1]]
    lo = min(p.lower() for p in ps); hi = max(p.upper() for p in ps)
    return (ar(lo), ar(hi))
def to_ball(a):
    return a[0].union(a[1])

def components(v, c, s):
    # v, c, s are BALLS; build (lo,hi) positive intervals
    V = pf(v); C = pf(c); S = pf(s)
    cs = pmul(C, S)
    oms = psub1(s)
    omt = ssub(pc(1), sadd(cs, pmul(V, V)))      # 1 - cs - v^2  (>0 admissible)
    omw = ssub(pc(1), cs)                          # 1 - cs
    # A5 = -v^2 + v*oms*c*(cs-1) + 1 + c*(1-2s) - c^2*s*oms   (signed)
    v2 = pmul(V, V); c2 = pmul(C, C)
    t1 = sneg(v2)
    t2 = smul(smul(smul(V, oms), C), ssub(cs, pc(1)))
    t3 = pc(1)
    t4 = smul(C, ssub(pc(1), padd(padd(S, S), pc(0))))   # c*(1-2s) ; 1-2s = 1-(s+s)
    t5 = sneg(smul(smul(c2, S), oms))
    A5 = sadd(sadd(sadd(sadd(t1, t2), t3), t4), t5)
    # rho7 = v^2 * c * oms^3 * A5^2 / (omw * omt^3)   (positive)
    A5sq = pmul(A5, A5)   # A5>0 admissible
    oms3 = ppow(oms, 3); omt3 = ppow(omt, 3)
    num = pmul(pmul(pmul(v2, C), oms3), A5sq)
    den = pmul(omw, omt3)
    rho7 = pdiv(num, den)
    R = proot7(rho7)
    # AB = (1+R)^2 * v^2 * c * oms / (R^5 * omt * omw)   (positive)
    one = pc(1)
    oneR = padd(one, R)
    oneR2 = ppow(oneR, 2)
    R5 = ppow(R, 5)
    num2 = pmul(pmul(pmul(oneR2, v2), C), oms)
    den2 = pmul(pmul(R5, omt), omw)
    AB = pdiv(num2, den2)
    sqAB = psqrt(AB)
    # T = 1/c + (s^2+2s-2)/oms + (2+cs*v-2cs)/v   (signed)
    s2 = pmul(S, S)
    s2p2sm2 = ssub(sadd(s2, padd(S, S)), pc(2))         # s^2+2s-2
    term_T1 = pdiv(pc(1), C)
    term_T2 = sdiv(s2p2sm2, oms)
    csv = pmul(cs, V)
    two_cs = padd(cs, cs)
    # 2 + c*s*v - 2*c*s = 2 + csv - 2*cs
    num3 = sadd(pc(2), ssub(csv, two_cs))
    term_T3 = sdiv(num3, V)
    T = sadd(sadd(term_T1, term_T2), term_T3)
    Cc = smul(smul(R, oneR), T)
    P = sadd(Cc, sadd(sqAB, sqAB))                       # C + 2*sqAB
    return dict(cs=cs, oms=oms, omt=omt, omw=omw, A5=A5, rho7=rho7, R=R,
                AB=AB, sqAB=sqAB, T=T, C=Cc, P=P)

def P_mv_lo(V, C, S):
    """Mean-value-form lower bound on P over the box (Arb lo/hi intervals)."""
    cb = components(V, C, S)
    vm, cm, sm = midpoint(V), midpoint(C), midpoint(S)
    mb = components(vm, cm, sm)
    R, sqAB, T, omt, omw, oms = cb['R'], cb['sqAB'], cb['T'], cb['omt'], cb['omw'], cb['oms']
    cs = cb['cs']
    v, c, s = V, C, S   # balls; use as positive intervals via pf in ops below
    # partials over the box (same factored form as certificate), as (lo,hi)
    two = pc(2)
    cs_m1 = ssub(cs, pc(1))
    Vv = pf(v); Cc = pf(c); Ss = pf(s)
    v2 = pmul(Vv, Vv); c2 = pmul(Cc, Cc); oms2 = ppow(oms, 2)
    RP1 = padd(pc(1), R)
    RP1R = pmul(R, RP1)
    dPdv = sadd(smul(RP1R, smul(two, sdiv(cs_m1, v2))),
                smul(sqAB, sadd(pdiv(two, Vv), pdiv(pmul(two, Vv), omt))))
    dPdc = sadd(smul(RP1R, sadd(sneg(pdiv(pc(1), c2)), sdiv(smul(Ss, ssub(Vv, two)), Vv))),
                smul(sqAB, sadd(sadd(pdiv(pc(1), Cc), pdiv(Ss, omt)), pdiv(Ss, omw))))
    dPds = sadd(smul(RP1R, sadd(sdiv(smul(Ss, ssub(two, Ss)), oms2), sdiv(smul(Cc, ssub(Vv, two)), Vv))),
                smul(sqAB, sadd(sadd(sneg(pdiv(pc(1), oms)), pdiv(Cc, omt)), pdiv(Cc, omw))))
    dPdr = sadd(smul(sadd(pc(1), padd(R, R)), T),
                smul(sqAB, ssub(pdiv(two, RP1), pdiv(pc(5), R))))
    Pm = mb['P']
    acc = sadd(sadd(sadd(
        smul(dPdv, ssub(Vv, pf(vm))),
        smul(dPdc, ssub(Cc, pf(cm)))),
        smul(dPds, ssub(Ss, pf(sm)))),
        smul(dPdr, ssub(R, mb['R'])))
    # P(box) subset Pm + acc  =>  lower bound = (Pm.lo + acc.lo)
    lo = (Pm[0] + acc[0]).lower()
    return lo, cb

def krawczyk(W, f, fp, S_or_C, which):
    """Iterative Krawczyk on f(var, ...)=0 in 1D.  Returns True if unique root."""
    Wc = W
    for _ in range(60):
        d = fp(Wc) if which == 'c' else fp(Wc, S_or_C[0], S_or_C[1])
        if d.contains(flint.arb(0)):
            return False
        m = midpoint(Wc)
        fm = f(m, S_or_C) if which == 'v' else f(m, S_or_C)
        K = m - fm / d
        if contains_int(Wc, K):
            return True
        Wn = intersect(K, Wc)
        if Wn is None:
            return False
        Wc = Wn
    return False

import sys
_cover_path = sys.argv[1] if len(sys.argv) > 1 else 'code/_hc_cover.json'
with open(_cover_path) as f:
    D = json.load(f)
L_C = fractions.Fraction(D['L_C'])
L_C_arb = flint.arb(L_C.numerator) / flint.arb(L_C.denominator)
print("Arb re-verification: %d pieces, L_C=%s, prec=%d" % (
    len(D['pieces']), D['L_C'], flint.ctx.prec), flush=True)

n_ok = 0; n_fail = 0; fails = []; arb_min = None; arb_min_piece = -1
from collections import Counter
rcnt = Counter()
for k, p in enumerate(D['pieces']):
    S = box_ball(p['s']); C = box_ball(p['C']); V = box_ball(p['V'])
    reasons = []
    # (i) Krawczyk uniqueness (iterative)
    if not krawczyk(C, lambda c, S: G_ball(c, S), lambda c: Gc_ball(c, S), S, 'c'):
        reasons.append("Kc-unique-fail")
    if not krawczyk(V, lambda v, CS: E2_ball(v, CS[0], CS[1]),
                    lambda v, C, S: E2v_ball(v, C, S), (C, S), 'v'):
        reasons.append("Kv-unique-fail")
    # (ii) admissibility
    plo, cb = P_mv_lo(V, C, S)
    if arb_min is None or plo < arb_min:
        arb_min = plo; arb_min_piece = k
    if not (C.lower() > 0): reasons.append("c<=0")
    if not (V.lower() > 0 and V.upper() < 1): reasons.append("v not in (0,1)")
    if not (cb['omt'][0] > 0): reasons.append("omt<=0")
    if not (cb['omw'][0] > 0): reasons.append("omw<=0")
    if not (cb['A5'][0] > 0): reasons.append("A5<=0")
    if not (cb['rho7'][0] > 0): reasons.append("rho7<=0")
    # (iii) MV-form P > L_C
    if not (plo > L_C_arb):
        reasons.append("P_MV<=L_C")
    for r in reasons:
        rcnt[r] += 1
    if reasons:
        n_fail += 1
        if len(fails) < 20: fails.append((k, reasons))
    else:
        n_ok += 1
    if (k + 1) % 500 == 0:
        print("  arb-checked %d/%d ok=%d fail=%d" % (k + 1, len(D['pieces']), n_ok, n_fail), flush=True)

print("\nArb-verified pieces: %d / %d" % (n_ok, len(D['pieces'])), flush=True)
print("Failed: %d" % n_fail, flush=True)
print("ARB GLOBAL MIN P_MV = %.6f at piece %d" % (float(arb_min), arb_min_piece), flush=True)
print("ARB GLOBAL MIN > 7   : %s" % (float(arb_min) > 7.0), flush=True)
print("REASON TALLY:", dict(rcnt), flush=True)
for f in fails:
    print("  FAIL", f, flush=True)
print("\nALL PIECES Arb-VERIFIED (Krawczyk-unique + admissibility + MV P > L_C):", n_fail == 0, flush=True)
print("DONE-ARBCHECKER", flush=True)
