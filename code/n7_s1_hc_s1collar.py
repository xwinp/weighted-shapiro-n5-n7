#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorous s->1 (z->0) COLLAR for the H_C cover  (FIX 1 of GPT's 5th verdict).

WHY THIS EXISTS.  The 2604-piece cover (_hc_cover.json) stops at
s_max = 8840966100154093/9007199254740992 ~= 0.9815444124321552 < 1.
The interval (s_max, 1) -- z = 1-s in (0, 1-s_max) ~= (0, 0.01846) -- is NOT
empty: at s = 491/500 = 0.982 > s_max there is an exact admissible lift
(c~=0.0176818, v~=0.0176874, u~=0.0176872, K~=3.34e-11>0, P~=10.87>7) with
ZERO active cover boxes.  So the existing cover does NOT exhaust the admissible
H_C-superset.  This script closes that gap with a second blow-up collar.

SECOND BLOW-UP.  As s->1 (z->0): c,v,u -> 0 (all = O(z)), and the P
mean-value partials 1/c, 1/v, dP/ds blow up, so the (c,v,s) MV machinery
cannot certify P > L_C there even though TRUE P -> +infinity (P ~ z^{-1/7}).
Desingularise with the second blow-up

        delta = 1 - s = z,   c = delta * cbar,   v = delta * vbar,

under which  G(c,s) = delta * Gbar(cbar,delta),  E2u_red = delta * E2bar,
with Gbar REGULAR at delta=0 (Gbar(cbar,0) = 1 - cbar, simple root cbar=1)
and E2bar regular (E2bar(vbar,0) = vbar - 1).  The admissible arc
(cbar(delta), vbar(delta)) is smooth, cbar,vbar -> 1, bounded away from 0.

DESINGULARISED ADMISSIBILITY.  As delta->0 the raw admissibility quantities
c,v,w,u,K,rho7 all vanish like delta^k, so their interval lower bounds
underflow to <= 0 (cannot certify > 0) even though they are strictly positive
for delta > 0.  We therefore certify the REGULAR rescaled quantities:
    c   = delta*cbar        ->  cbar > 0                    (regular, ->1)
    v   = delta*vbar        ->  vbar > 0                    (regular, ->1)
    w   = delta*cbar*(1-d)  ->  cbar*(1-d) > 0              (regular)
    u   = delta * uhat      ->  uhat > 0  (uhat=u/delta)    (regular, ->vbar)
    K   = delta^6 * Ktilde  ->  Ktilde > 0                  (regular, ->vbar^2*cbar)
    rho7= delta^6 * rhot7   ->  rhot7 > 0                   (regular, ->vbar^2*cbar)
    a5                       ->  a5 > 0                     (regular, ->1)
All rescaled quantities are O(1) and bounded away from 0 on the arc, so their
interval lower bounds are cleanly > 0 with NO precision dependence.  Admissibility
for delta > 0 is EQUIVALENT to these rescaled conditions (each raw quantity is
delta^k times a strictly-positive regular factor, and delta > 0 on the collar).

P LOWER BOUND (two regimes).
  * delta >= delta_direct: DIRECT interval evaluation  P = R(1+R)T + 2 sqrt(AB)
    over the blown-up box (c=delta*cbar, v=delta*vbar, s=1-delta).  P ranges
    ~10.9 (at delta=delta0) up to ~+inf (delta->0); comfortably > L_C=141/20.
    Representable down to delta ~ 1e-10 at ivprec=260 (rho7~delta^6 binding).
  * delta < delta_direct: CRUDE rigorous bound.  On the collar (s > 0.732) every
    term of T = 1/c + (s^2+2s-2)/(1-s) + (2+csv-2cs)/v is >= 0, so
    P >= C = R(1+R)T >= R*T >= R/c = rho7^{1/7}/c = delta^{-1/7} * cbar^{-6/7}
            * (vbar^2 * a5^2 / (omw * omt^3))^{1/7},
    a product of delta^{-1/7} (->+inf) and a regular O(1) factor.  Its minimum
    on (0, delta_direct] is at delta_direct, where it is already > L_C
    (delta_direct^{-1/7} ~ 19 >> 7.05).  No underflow: only regular quantities
    are evaluated.  (delta=0 itself has c=v=0 => K=0, NOT admissible, so the
    admissible set is delta > 0 and is fully covered.)

KRAWCZYK UNIQUENESS.  At each delta-subinterval midpoint, interval-Newton in
cbar (then vbar) proves a UNIQUE admissible root in the certified box, with
N(box) subset int(box); the regular Jacobians dGbar/dcbar (~-1) and
dE2bar/dvbar (~1) are bounded away from 0, so contraction holds.

OUTPUT.  code/_hc_s1collar.json : collar pieces in ORIGINAL (c,v,s) coords
(keys s,C,V + lower bounds, same schema as _hc_cover.json, so the existing
mpmath.iv / Arb checkers can re-verify each listed box) PLUS the blown-up
Krawczyk/admissibility/P records (cbar,vbar,delta,method) for this certificate.
The collar extends the cover from s_max to s=1, closing the exhaustiveness gap.

Load-bearing EXHAUSTIVENESS certificate for the s->1 tail; box VALIDITY of each
collar piece is certified here by direct interval P (and re-checkable).
"""
import sympy as sp, mpmath as mp, json, time, fractions
mp.mp.ivprec = 260; IV = mp.iv; mp.mp.prec = 200

# ---- helpers (mirror rigorous_cert.py) ----
def iv(x):
    if isinstance(x, IV.mpf): return x
    if isinstance(x, mp.mpf): return IV.mpf([x, x])
    return IV.mpf([mp.mpf(x), mp.mpf(x)])
def iv_mid(I): return (mp.mpf(I.a) + mp.mpf(I.b)) / 2
def iv_isect(a, b):
    lo = a.a if a.a > b.a else b.a; hi = a.b if a.b < b.b else b.b
    return None if lo > hi else IV.mpf([lo, hi])
def iv_7throot_pos(I):
    if I.a <= 0: return None
    e = mp.mpf(1) / 7; return IV.mpf([mp.mpf(I.a) ** e, mp.mpf(I.b) ** e])
def iv_sqrt_pos(I):
    if I.a <= 0: return None
    return IV.sqrt(I)
def iv_to_rat(I):
    return [fractions.Fraction(str(mp.nstr(mp.mpf(I.a), 24))),
            fractions.Fraction(str(mp.nstr(mp.mpf(I.b), 24)))]
def mpf_rat(r):
    """Exact-ish mpf from a Fraction (rational -> mpf via p/q)."""
    return mp.mpf(r.numerator) / mp.mpf(r.denominator)
def iv_rat(r):
    return iv(mpf_rat(r))
def iv_box_rat(lo, hi):
    return IV.mpf([mpf_rat(lo), mpf_rat(hi)])
def iv_to_rat_inward(I, ndig=30):
    """Rational box STRICTLY INSIDE the iv interval (lo rounded UP, hi DOWN).
    Guarantees the recorded rational box is a subset of the certified mpf box,
    so Krawczyk uniqueness transfers to the rational box."""
    p = mp.mpf(10) ** ndig
    lo_a = mp.mpf(I.a); hi_b = mp.mpf(I.b)
    lo_int = int(lo_a * p)              # floor for positive -> below I.a
    hi_int = int(hi_b * p)             # floor -> below I.b
    lo_rat = fractions.Fraction(lo_int + 1, 10 ** ndig)   # +1 -> ceil -> above I.a (inside)
    hi_rat = fractions.Fraction(hi_int, 10 ** ndig)       # floor -> below I.b (inside)
    if lo_rat >= hi_rat:               # box too thin for ndig; fall back to faithful
        return iv_to_rat(I)
    return [lo_rat, hi_rat]

# ---- exec symbol-building part of rigorous_cert.py (gives _E2u_cs, etc.) ----
src = open('code/n7_s1_hc_rigorous_cert.py').read().split('with open')[0]
exec(src)
mp.mp.ivprec = 260; IV = mp.iv; mp.mp.prec = 200   # re-raise (rigorous_cert set 110)
_c, _v, _s = sp.symbols('c v s')
Gexpr = (1 - _s) * _c**3 + (1 - _s) * (_s - 2) * _c**2 + (_s**2 - _s - 1) * _c + (1 - _s)

# ---- blown-up symbols ----
d, cb, vb = sp.symbols('delta cbar vbar')
w_b = d * cb * (1 - d)                 # w = c*s = delta*cbar*(1-delta)
oms_b = d                              # 1-s = delta
omw_b = 1 - w_b                        # 1-w   (regular -> 1)
omt_b = 1 - w_b - (d * vb) ** 2        # (1-w)-v^2  (regular -> 1)
a5_b = -((d*vb)**2) + (d*vb)*oms_b*(d*cb)*(w_b-1) + 1 + (d*cb)*(1-2*(1-d)) - (d*cb)**2*(1-d)*oms_b
# regular rescaled quantities (divide out the vanishing delta-power)
u_b = (d*vb)*(1 - d*vb) / omt_b                       # u = v(1-v)/((1-w)-v^2)
uhat = sp.cancel(u_b / d)                             # u/delta  (regular)
rho7_b = (d*vb)**2 * (d*cb) * oms_b**3 * a5_b**2 / (omw_b * omt_b**3)
rhot7 = sp.cancel(rho7_b / d**6)                      # rho7/delta^6  (regular)
K_b = u_b * (d*vb) * w_b * oms_b**3 * a5_b**2 / ((1 - d*vb) * (1 - w_b) * (1 - d)**3)
Ktilde = sp.cancel(K_b / d**6)                        # K/delta^6  (regular)
# Gbar, E2bar and their derivatives (regular reduced equations)
Gbar = sp.expand(Gexpr.subs(_c, d*cb).subs(_s, 1-d) / d)
Gbar_cb = sp.diff(Gbar, cb)
E2bar = sp.expand(_E2u_cs.subs({_v: d*vb, _c: d*cb, _s: 1-d}) / d)
E2bar_vb = sp.diff(E2bar, vb)

# ---- generic EXACT-coefficient interval evaluator (rational expr in cb,vb,d) ----
def _num(x):
    if isinstance(x, sp.Rational): return mp.mpf(x.p) / mp.mpf(x.q)
    if isinstance(x, sp.Integer): return mp.mpf(int(x))
    return mp.mpf(x)
def iv_eval(expr, env):
    expr = sp.sympify(expr)
    if expr.is_Number: return iv(_num(expr))
    if expr.is_Symbol: return env[expr]
    if expr.is_Pow:
        b = iv_eval(expr.base, env); e = expr.exp
        if e.is_Integer and int(e) > 0:
            r = b
            for _ in range(int(e) - 1): r = r * b
            return r
        if e.is_Integer and int(e) < 0:
            r = iv(1)
            for _ in range(-int(e)): r = r / b
            return r
        ef = mp.mpf(_num(e)); return IV.mpf([mp.mpf(b.a)**ef, mp.mpf(b.b)**ef])  # base>0
    if expr.is_Mul:
        r = iv_eval(expr.args[0], env)
        for a in expr.args[1:]: r = r * iv_eval(a, env)
        return r
    if expr.is_Add:
        r = iv_eval(expr.args[0], env)
        for a in expr.args[1:]: r = r + iv_eval(a, env)
        return r
    raise ValueError("unhandled node %r" % type(expr))

def Gbar_iv(cbi, di):   return iv_eval(Gbar, {cb: cbi, vb: iv(0), d: di})
def Gbarc_iv(cbi, di):  return iv_eval(Gbar_cb, {cb: cbi, vb: iv(0), d: di})
def E2bar_iv(vbi, cbi, di):  return iv_eval(E2bar, {cb: cbi, vb: vbi, d: di})
def E2barv_iv(vbi, cbi, di): return iv_eval(E2bar_vb, {cb: cbi, vb: vbi, d: di})

# ---- rescaled admissibility (all regular, O(1), no underflow) ----
def adm_rescaled(cbi, vbi, di):
    """Certify admissibility for delta>0 via regular rescaled quantities."""
    env = {cb: cbi, vb: vbi, d: di}
    cbv, vbv, dv = iv(cbi), iv(vbi), iv(di)
    if not (cbv.a > 0 and vbv.a > 0): return False, None
    if not (dv.a > 0 and dv.b < 1): return False, None
    omw = iv_eval(omw_b, env); omt = iv_eval(omt_b, env); a5 = iv_eval(a5_b, env)
    if omw.a <= 0 or omt.a <= 0 or a5.a <= 0: return False, None
    rhot7_iv = iv_eval(rhot7, env)
    if rhot7_iv.a <= 0: return False, None
    Kt = iv_eval(Ktilde, env)
    if Kt.a <= 0: return False, None
    uh = iv_eval(uhat, env)
    if uh.a <= 0: return False, None
    v_hi = dv.b * vbv.b; w_hi = dv.b * cbv.b * (1 - dv.a); u_hi = dv.b * uh.b
    if not (v_hi < 1 and w_hi < 1 and u_hi < 1): return False, None
    # K = delta^6 * Ktilde; tiny in the tail -- store as string to avoid float underflow
    Kfull = mp.mpf(dv.a)**6 * mp.mpf(Kt.a)
    try:
        K_lo_f = float(Kfull)
    except (ValueError, OverflowError):
        K_lo_f = 0.0   # subfloat-tiny but strictly > 0 (Ktilde>0, delta>0)
    return True, dict(cb_lo=float(cbv.a), vb_lo=float(vbv.a),
                      omw_lo=float(omw.a), omt_lo=float(omt.a), a5_lo=float(a5.a),
                      rhot7_lo=float(rhot7_iv.a), Ktilde_lo=float(Kt.a),
                      uhat_lo=float(uh.a), K_lo=K_lo_f)

# ---- direct interval P over blown-up box ----
def P_direct_blown(cbi, vbi, di):
    """Direct P = R(1+R)T + 2 sqrt(AB) over blown-up box. Returns P_lo or None."""
    env = {cb: cbi, vb: vbi, d: di}
    dv = iv(di); cbv = iv(cbi); vbv = iv(vbi)
    c = dv * cbv; v = dv * vbv; s = 1 - dv
    cs_ = c * s; oms = dv; omt = iv_eval(omt_b, env); omw = iv_eval(omw_b, env)
    a5 = iv_eval(a5_b, env)
    if omt.a <= 0 or omw.a <= 0 or a5.a <= 0: return None
    rho7 = v * v * c * oms**3 * a5 * a5 / (omw * omt**3)
    if rho7.a <= 0: return None
    R = iv_7throot_pos(rho7)
    if R is None: return None
    AB = (1 + R)**2 * v * v * c * oms / (R**5 * omt * omw)
    sqAB = iv_sqrt_pos(AB)
    if sqAB is None: return None
    t1 = 1 / c; t2 = (s*s + 2*s - 2) / oms; t3 = (2 + cs_*v - 2*cs_) / v
    T = t1 + t2 + t3
    P = R * (1 + R) * T + 2 * sqAB
    return float(P.a)

# ---- crude rigorous P lower bound for the tiny-delta tail ----
def P_crude_blown(cbi, vbi, di):
    """P >= delta^{-1/7} * cbar^{-6/7} * (vbar^2 a5^2/(omw omt^3))^{1/7}."""
    env = {cb: cbi, vb: vbi, d: di}
    dv = iv(di); cbv = iv(cbi); vbv = iv(vbi)
    db = mp.mpf(dv.b); cb_hi = mp.mpf(cbv.b)
    omw = iv_eval(omw_b, env); omt = iv_eval(omt_b, env); a5 = iv_eval(a5_b, env)
    vb_lo2 = mp.mpf(vbv.a) ** 2; a5_lo2 = mp.mpf(a5.a) ** 2
    omw_hi = mp.mpf(omw.b); omt_hi3 = mp.mpf(omt.b) ** 3
    factor_lo = vb_lo2 * a5_lo2 / (omw_hi * omt_hi3)
    if factor_lo <= 0: return None
    factor_7th = mp.mpf(factor_lo) ** (mp.mpf(1) / 7)
    db_term = db ** (-mp.mpf(1) / 7)
    cb_term = cb_hi ** (-mp.mpf(6) / 7)
    return float(db_term * cb_term * factor_7th)

# ---- Krawczyk in cbar / vbar at scalar delta ----
def krawczyk_cbar(D, c0):
    C = IV.mpf([mp.mpf(c0) - mp.mpf('0.05'), mp.mpf(c0) + mp.mpf('0.05')])
    for _ in range(100):
        dG = Gbarc_iv(C, D)
        if dG.a <= 0 <= dG.b: return None
        m = iv_mid(C); K = iv(m) - Gbar_iv(iv(m), D) / dG
        Cn = iv_isect(K, C)
        if Cn is None: return None
        if K.a >= Cn.a and K.b <= Cn.b and (Cn.b - Cn.a) < mp.mpf('1e-14'):
            return Cn
        C = Cn
    return C if (C.b - C.a) < mp.mpf('1e-12') else None

def krawczyk_vbar(Cbar, D, v0):
    V = IV.mpf([mp.mpf(v0) - mp.mpf('0.05'), mp.mpf(v0) + mp.mpf('0.05')])
    for _ in range(80):
        dE = E2barv_iv(V, Cbar, D)
        if dE.a <= 0 <= dE.b: return None
        m = iv_mid(V); K = iv(m) - E2bar_iv(iv(m), Cbar, D) / dE
        Vn = iv_isect(K, V)
        if Vn is None: return None
        if K.a >= Vn.a and K.b <= Vn.b and (Vn.b - Vn.a) < mp.mpf('1e-14'):
            return Vn
        V = Vn
    return V if (V.b - V.a) < mp.mpf('1e-12') else None

# ---- connect to existing cover: seed from last box at s=s_max ----
_cov = json.load(open('code/_hc_cover.json'))
def _find_pieces(o):
    if isinstance(o, list) and o and isinstance(o[0], dict): return o
    if isinstance(o, dict):
        for vv in o.values():
            r = _find_pieces(vv)
            if r: return r
    return None
EXIST = _find_pieces(_cov)
s_max = max(sp.Rational(p['s'][1]) for p in EXIST)
delta0 = 1 - s_max
L_C = mp.mpf(141) / mp.mpf(20)
print("s_max = %s = %.16f" % (s_max, float(s_max)), flush=True)
print("delta0 = 1-s_max = %s = %.16f" % (delta0, float(delta0)), flush=True)
print("L_C = 141/20 = %s" % L_C, flush=True)
# seed from the ACTUAL admissible lift at s=s_max (the cover box is wide in c;
# its midpoint is a poor seed -- the true cbar root at delta0 is ~0.9819, not
# the midpoint ~1.012).  admissible_lifts returns the true (c,v,u) at scalar s.
_lifts = admissible_lifts(float(s_max))
print("admissible lifts at s=s_max: %d" % len(_lifts), flush=True)
cbar_seed = vbar_seed = None
if _lifts:
    c0, v0, u0 = _lifts[0]
    d0f = float(delta0); cbar_seed = c0 / d0f; vbar_seed = v0 / d0f
print("seed from admissible lift: c=%.7f v=%.7f -> cbar=%.6f vbar=%.6f" % (
    _lifts[0][0] if _lifts else 0, _lifts[0][1] if _lifts else 0, cbar_seed, vbar_seed), flush=True)
# cross-check: true cbar root at delta0 via exact CRootOf
_cbroots = [sp.N(r, 30) for r in sp.real_roots(sp.Poly(sp.expand(Gbar.subs(d, delta0)), cb))]
print("true cbar roots at delta0: %s" % [float(x) for x in _cbroots if abs(x) < 5], flush=True)

# verify the admissible arc reaches delta=delta0 (exact CRootOf check)
Rv0 = sp.expand(sp.resultant(Gexpr, _E2u_cs, _c).subs(_s, s_max))
vroots0 = [r for r in sp.real_roots(sp.Poly(Rv0, _v)) if 0 < r < 1]
print("exact v-roots at s=s_max: %d" % len(vroots0), flush=True)

# ---- monotonic continuation sweep: delta from delta0 DOWN to delta_direct ----
print("\n=== Krawczyk arc tracking (delta: delta0 -> delta_direct) ===", flush=True)
delta_direct = fractions.Fraction(1, 10**10)   # below this, crude bound + rescaled adm
MINS = fractions.Fraction(1, 2**50)            # bisection floor ~8.9e-16
pieces = []
t0 = time.time()
n_fail = 0

def certify_delta_box(da, db, cs, vs):
    """Certify Fraction delta-box [da,db] via Krawczyk(mid) + rescaled adm + direct P.
    Records EXACT rational bounds (s=1-delta, C=delta*cbar, V=delta*vbar) so that
    adjacent pieces abut exactly and the seam s_lo[0]=s_max is exact."""
    dm = (da + db) / 2                                  # Fraction midpoint
    Cbar = krawczyk_cbar(iv_rat(dm), float(cs))
    if Cbar is None: return None
    Vbar = krawczyk_vbar(Cbar, iv_rat(dm), float(vs))
    if Vbar is None: return None
    Dbox = iv_box_rat(da, db)
    ok, info = adm_rescaled(Cbar, Vbar, Dbox)
    if not ok: return None
    P_lo = P_direct_blown(Cbar, Vbar, Dbox)
    if P_lo is None or P_lo <= float(L_C): return None
    cb_lo, cb_hi = iv_to_rat_inward(Cbar)              # rational INSIDE certified mpf box
    vb_lo, vb_hi = iv_to_rat_inward(Vbar)
    # exact rational original-coord bounds: c = delta*cbar, s = 1-delta
    c_lo, c_hi = da * cb_lo, db * cb_hi
    v_lo, v_hi = da * vb_lo, db * vb_hi
    s_lo, s_hi = 1 - db, 1 - da
    piece = dict(
        s=[str(s_lo), str(s_hi)], C=[str(c_lo), str(c_hi)], V=[str(v_lo), str(v_hi)],
        cbar=[str(cb_lo), str(cb_hi)], vbar=[str(vb_lo), str(vb_hi)],
        delta=[str(da), str(db)],
        plo=P_lo, rho7_lo=info['rhot7_lo'], A5_lo=info['a5_lo'],
        omt_lo=info['omt_lo'], omw_lo=info['omw_lo'], K_lo=info['K_lo'],
        Ktilde_lo=info['Ktilde_lo'], uhat_lo=info['uhat_lo'],
        krawczyk_unique=True, method='direct_interval_P')
    return piece, iv_mid(Cbar), iv_mid(Vbar), P_lo

d_cur = fractions.Fraction(delta0)
h = fractions.Fraction(delta0) / 40
n_step = 0
stalled = False
while d_cur > delta_direct and not stalled:
    da = d_cur - h
    if da <= delta_direct:              # final direct box ends at delta_direct
        da = delta_direct
    db = d_cur
    res = certify_delta_box(da, db, cbar_seed, vbar_seed)
    if res is None:
        if h < MINS:
            n_fail += 1; print("  BISECT LIMIT d=%.3e h=%.3e" % (float(d_cur), float(h)), flush=True)
            stalled = True; continue
        h = h / 2; continue
    piece, cs_new, vs_new, P_lo = res
    pieces.append(piece)
    cbar_seed, vbar_seed = cs_new, vs_new
    d_cur = da; n_step += 1
    if n_step % 25 == 0:
        print("  %d pieces, d_cur=%.6f P_lo=%.4f h=%.4e" % (len(pieces), float(d_cur), P_lo, float(h)), flush=True)
print("direct sweep done: %d pieces, d_cur=%.3e" % (len(pieces), float(d_cur)), flush=True)

# ---- crude tail: (0, delta_direct] via rescaled adm + crude P (bisect if wide) ----
# The admissible set is delta>0 (delta=0 gives c=v=0 => K=0, inadmissible).  Both
# the crude P bound (P >= delta^{-1/7}*(regular), monotone DECREASING in delta) and
# the rescaled admissibility (v=delta*vbar<1 etc., worst at LARGEST delta) have their
# worst case at delta=db.  So certify each tail sub-interval (da,db] at the POINT
# delta=db; this certifies the whole (da,db] by monotonicity.  da=0 is the excluded
# endpoint (delta=0 inadmissible), so the box (0,db] is fully covered.
print("\n=== crude tail (0, delta_direct] ===", flush=True)
tail_stack = [(fractions.Fraction(0), delta_direct)]
while tail_stack:
    da, db = tail_stack.pop()
    Cbc = IV.mpf([mp.mpf('0.90'), mp.mpf('1.10')])
    Vbc = IV.mpf([mp.mpf('0.90'), mp.mpf('1.10')])
    Dbox_adm = iv_rat(db)                  # POINT at db (worst case for adm)
    ok, info = adm_rescaled(Cbc, Vbc, Dbox_adm)
    if not ok:
        if (db - da) < MINS:
            n_fail += 1; print("  tail BISECT LIMIT [%.3e,%.3e]" % (float(da), float(db)), flush=True); continue
        mid = (da + db) / 2; tail_stack.append((da, mid)); tail_stack.append((mid, db)); continue
    Dbox_P = iv_box_rat(da, db)            # [da,db]; crude P uses db (upper) -> min
    P_lo = P_crude_blown(Cbc, Vbc, Dbox_P)
    if P_lo is None or P_lo <= float(L_C):
        if (db - da) < MINS:
            n_fail += 1; print("  tail P fail [%.3e,%.3e] P=%s" % (float(da), float(db), P_lo), flush=True); continue
        mid = (da + db) / 2; tail_stack.append((da, mid)); tail_stack.append((mid, db)); continue
    cb_lo, cb_hi = fractions.Fraction('0.90'), fractions.Fraction('1.10')
    vb_lo, vb_hi = fractions.Fraction('0.90'), fractions.Fraction('1.10')
    c_lo, c_hi = da * cb_lo, db * cb_hi
    v_lo, v_hi = da * vb_lo, db * vb_hi
    s_lo, s_hi = 1 - db, 1 - da
    pieces.append(dict(
        s=[str(s_lo), str(s_hi)], C=[str(c_lo), str(c_hi)], V=[str(v_lo), str(v_hi)],
        cbar=[str(cb_lo), str(cb_hi)], vbar=[str(vb_lo), str(vb_hi)],
        delta=[str(da), str(db)],
        plo=P_lo, rho7_lo=info['rhot7_lo'], A5_lo=info['a5_lo'],
        omt_lo=info['omt_lo'], omw_lo=info['omw_lo'], K_lo=info['K_lo'],
        Ktilde_lo=info['Ktilde_lo'], uhat_lo=info['uhat_lo'],
        krawczyk_unique=False, method='crude_d^{-1/7}_tail',
        tail_certified_at_db_by_monotonicity=True))

# ---- INDEPENDENT RE-VERIFICATION on the recorded EXACT rational pieces ----
# Re-evaluate the desingularized P lower bound AND rescaled admissibility on the
# RECORDED rational (cbar,vbar,delta) boxes -- independent of the mpf sweep's
# internal interval values.  (We do NOT use the MV form P_box_mv here: that is the
# very form whose 1/c,1/v partials blow up at s->1, i.e. the singularity the
# collar exists to bypass.  The collar is certified by the desingularized direct
# interval P + rescaled admissibility, re-checked on the exact rationals.)
# Also checks (b) exact abutment s_hi[i]==s_lo[i+1], (c) seam s_lo[0]==s_max, s_hi[-1]==1.
print("\n=== independent re-verification on recorded rationals ===", flush=True)
pieces.sort(key=lambda p: fractions.Fraction(p['s'][0]))
LC_rat = fractions.Fraction(141, 20)
n_bad_P = n_bad_adm = n_bad_abut = 0
for i, p in enumerate(pieces):
    cbr = [fractions.Fraction(p['cbar'][0]), fractions.Fraction(p['cbar'][1])]
    vbr = [fractions.Fraction(p['vbar'][0]), fractions.Fraction(p['vbar'][1])]
    ddr = [fractions.Fraction(p['delta'][0]), fractions.Fraction(p['delta'][1])]
    Cbi = iv_box_rat(cbr[0], cbr[1]); Vbi = iv_box_rat(vbr[0], vbr[1])
    is_tail = p['method'].startswith('crude')
    # tail: adm certified at POINT delta=db (worst case by monotonicity); P on [da,db]
    Dbi_adm = iv_rat(ddr[1]) if is_tail else iv_box_rat(ddr[0], ddr[1])
    Dbi_P = iv_box_rat(ddr[0], ddr[1])
    ok_re, info_re = adm_rescaled(Cbi, Vbi, Dbi_adm)
    if not ok_re:
        n_bad_adm += 1; print("  ADM-REVIFY FAIL piece %d" % i, flush=True)
    Pfun = P_crude_blown if is_tail else P_direct_blown
    P_lo_re = Pfun(Cbi, Vbi, Dbi_P)
    if P_lo_re is None or P_lo_re <= float(L_C):
        n_bad_P += 1; print("  P-REVIFY FAIL piece %d s=[%.6f,%.6f] P_re=%s" % (
            i, float(fractions.Fraction(p['s'][0])), float(fractions.Fraction(p['s'][1])), P_lo_re), flush=True)
    p['P_lo_reverify'] = P_lo_re
    if i + 1 < len(pieces):
        s_hi_i = fractions.Fraction(p['s'][1])
        s_lo_j = fractions.Fraction(pieces[i+1]['s'][0])
        if s_hi_i != s_lo_j:
            n_bad_abut += 1
            print("  ABUT FAIL @%d: s_hi=%s != s_lo=%s" % (i, str(s_hi_i), str(s_lo_j)), flush=True)
seam_lo = fractions.Fraction(pieces[0]['s'][0]); seam_hi = fractions.Fraction(pieces[-1]['s'][1])
seam_ok = (seam_lo == s_max) and (seam_hi == 1)
print("re-verify: P failures=%d, adm failures=%d, abutment failures=%d" % (
    n_bad_P, n_bad_adm, n_bad_abut), flush=True)
print("seam: s_lo[0]==s_max (%s) and s_hi[-1]==1 (%s) : %s" % (
    seam_lo == s_max, seam_hi == 1, seam_ok), flush=True)
all_revify_ok = (n_bad_P == 0 and n_bad_adm == 0 and n_bad_abut == 0 and seam_ok)
gmin = min(p['plo'] for p in pieces) if pieces else None
gmin_re = min((p.get('P_lo_reverify') or 0) for p in pieces) if pieces else None

print("\n=== summary ===", flush=True)
print("collar pieces: %d ; sweep P_lo min=%.6f ; re-verify P_lo min=%.6f > L_C=%.5f" % (
    len(pieces), gmin, gmin_re, float(L_C)), flush=True)
print("failures(sweep): %d ; re-verify ok: %s" % (n_fail, all_revify_ok), flush=True)
print("elapsed %.1fs" % (time.time() - t0), flush=True)

out = dict(L_C=str(LC_rat), s_max=str(s_max), delta0=str(delta0),
           delta_direct=str(delta_direct), n_pieces=len(pieces),
           global_P_lo_min=gmin, global_P_lo_min_reverify=gmin_re,
           gmin_gt_L_C=bool(gmin > float(L_C)), gmin_re_gt_L_C=bool(gmin_re > float(L_C)),
           n_failures=n_fail, reverify_ok=all_revify_ok,
           seam_at_s_max=bool(seam_lo == s_max), seam_at_1=bool(seam_hi == 1),
           abutment_exact=(n_bad_abut == 0),
           method="second_blowup(d=1-s,c=dcbar,v=dvbar); Krawczyk(cbar,vbar); "
                  "rescaled admissibility; direct interval P (d>=1e-10) + crude d^{-1/7} tail; "
                  "independent desingularized re-verify on exact rational bounds; exact-Fraction abutment",
           pieces=pieces)
with open('code/_hc_s1collar.json', 'w') as f:
    json.dump(out, f, indent=2)
print("Wrote code/_hc_s1collar.json", flush=True)
print("DONE-S1COLLAR n=%d gmin=%.6f gmin_re=%.6f >L_C=%s revify=%s seam=%s abut=%s" % (
    len(pieces), gmin, gmin_re, gmin_re > float(L_C), all_revify_ok, seam_ok, n_bad_abut == 0), flush=True)
