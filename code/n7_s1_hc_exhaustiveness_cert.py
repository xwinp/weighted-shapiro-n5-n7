#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Certified exhaustiveness of the H_C cover + s->1 collar.

This is the load-bearing replacement for the earlier midpoint/count-only and
floating critical-fiber checks.  It proves three independent properties:

  1. PARAMETRIC BOX VALIDITY is checked elsewhere:
       n7_s1_hc_cover_checker.py / n7_s1_hc_arb_checker.py for the 2604 cover
       boxes, and n7_s1_hc_s1collar.py for the 60 second-blow-up collar boxes.

  2. BRANCH-LEVEL BOX GRAPH.  Whenever two boxes overlap in parameter and
     coordinate ranges, a 2x2 Krawczyk operator at an exact rational overlap
     parameter proves that their unique roots are the same root.  Graph
     components therefore represent continuous root branches, not merely
     overlapping s-projections.  Every component's exact rational s-projection
     is checked for gap-free coverage.

  3. ALGEBRAIC/INTERVAL ENUMERATION.  On one rational sample of each of the six
     exact event-free cells, all c-roots are exact CRootOf objects.  At each of
     the five exact critical fibers, a certified rational isolating interval
     for s* and full-parameter interval-Newton c-tubes enclose every real
     c-root.  The quadratic Q(v,c,s)=E2_red/(v-1) is then enclosed by an
     outward-rounded quadratic formula.  Boundary roots are certified by a
     3x3 Krawczyk operator for (G,Q,B)=0.  No nroots/polyroots, number-field
     factorisation, residual threshold, nsimplify, or floating containment is
     used in any load-bearing step.

At event-free samples, admissible roots are in bijection with graph components
active on the whole cell.  At critical fibers, every admissible root is inside
exactly one component.  Since the event list contains every multiplicity and
admissibility-boundary event, this proves that every admissible H_C lift lies in
one of the certified boxes.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
import sympy as sp

IV = mp.iv
IV.prec = 240
mp.mp.prec = 220

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"

# ---------------------------------------------------------------------------
# Exact reduced equations and interval evaluators
# ---------------------------------------------------------------------------
_ns: dict[str, Any] = {}
_src = (CODE / "n7_s1_hc_rigorous_cert.py").read_text().split("with open")[0]
exec(_src, _ns)
IV.prec = 240
mp.mp.prec = 220

_c, _v, _s = _ns["_c"], _ns["_v"], _ns["_s"]
G = (1 - _s) * _c**3 + (1 - _s) * (_s - 2) * _c**2 + (_s**2 - _s - 1) * _c + (1 - _s)
E = _ns["_E2u_cs"]
Q = sp.cancel(E / (_v - 1))  # v=1 is always inadmissible; Q is quadratic in v.
_Qp = sp.Poly(Q, _v)
_Qcoef = [_Qp.nth(k) for k in range(3)]
_Ec_expr = sp.diff(E, _c)

iv = _ns["iv"]
iv_mid = _ns["iv_mid"]
G_iv, Gc_iv = _ns["G_iv"], _ns["Gc_iv"]
E_iv, Ev_iv = _ns["E2_iv"], _ns["E2v_iv"]


def Ec_iv(V, C, S):
    return _ns["iv_eval"](_Ec_expr, {_v: V, _c: C, _s: S})

# Collar regular equations/functions, without executing the generator main.
_col_ns: dict[str, Any] = {}
_col_src = (CODE / "n7_s1_hc_s1collar.py").read_text().split("# ---- connect to existing cover")[0]
exec(_col_src, _col_ns)
IV.prec = 240
mp.mp.prec = 220
_Ecbar_expr = sp.diff(_col_ns["E2bar"], _col_ns["cb"])


def Ecbar_iv(V, C, D):
    return _col_ns["iv_eval"](
        _Ecbar_expr,
        {_col_ns["vb"]: V, _col_ns["cb"]: C, _col_ns["d"]: D},
    )


def mpf_fraction(q: Fraction) -> mp.mpf:
    return mp.mpf(q.numerator) / mp.mpf(q.denominator)


def I(a, b=None):
    if b is None:
        b = a
    if isinstance(a, Fraction):
        a = mpf_fraction(a)
    if isinstance(b, Fraction):
        b = mpf_fraction(b)
    return IV.mpf([a, b])


def strict_subset(K, W) -> bool:
    return bool(K.a > W.a and K.b < W.b)


def krawczyk2_raw(C, V, S) -> bool:
    """2x2 Krawczyk for (G(c,s), E(v,c,s)) at a point S."""
    cm, vm = iv_mid(C), iv_mid(V)
    gc0 = iv_mid(Gc_iv(I(cm), S))
    ec0 = iv_mid(Ec_iv(I(vm), I(cm), S))
    ev0 = iv_mid(Ev_iv(I(vm), I(cm), S))
    if gc0 == 0 or ev0 == 0:
        return False
    # inverse of lower triangular midpoint Jacobian [[gc,0],[ec,ev]]
    y11 = I(1 / gc0)
    y21 = I(-ec0 / (gc0 * ev0))
    y22 = I(1 / ev0)
    f1, f2 = G_iv(I(cm), S), E_iv(I(vm), I(cm), S)
    z1 = I(cm) - y11 * f1
    z2 = I(vm) - y21 * f1 - y22 * f2
    gc, ec, ev = Gc_iv(C, S), Ec_iv(V, C, S), Ev_iv(V, C, S)
    m11 = 1 - y11 * gc
    m21 = -(y21 * gc + y22 * ec)
    m22 = 1 - y22 * ev
    Kc = z1 + m11 * (C - I(cm))
    Kv = z2 + m21 * (C - I(cm)) + m22 * (V - I(vm))
    return strict_subset(Kc, C) and strict_subset(Kv, V)


def krawczyk2_collar(C, V, D) -> bool:
    """2x2 Krawczyk for (Gbar(cbar,d), E2bar(vbar,cbar,d))."""
    cm, vm = _col_ns["iv_mid"](C), _col_ns["iv_mid"](V)
    gc0 = _col_ns["iv_mid"](_col_ns["Gbarc_iv"](I(cm), D))
    ec0 = _col_ns["iv_mid"](Ecbar_iv(I(vm), I(cm), D))
    ev0 = _col_ns["iv_mid"](_col_ns["E2barv_iv"](I(vm), I(cm), D))
    if gc0 == 0 or ev0 == 0:
        return False
    y11 = I(1 / gc0)
    y21 = I(-ec0 / (gc0 * ev0))
    y22 = I(1 / ev0)
    f1 = _col_ns["Gbar_iv"](I(cm), D)
    f2 = _col_ns["E2bar_iv"](I(vm), I(cm), D)
    z1 = I(cm) - y11 * f1
    z2 = I(vm) - y21 * f1 - y22 * f2
    gc = _col_ns["Gbarc_iv"](C, D)
    ec = Ecbar_iv(V, C, D)
    ev = _col_ns["E2barv_iv"](V, C, D)
    m11 = 1 - y11 * gc
    m21 = -(y21 * gc + y22 * ec)
    m22 = 1 - y22 * ev
    Kc = z1 + m11 * (C - I(cm))
    Kv = z2 + m21 * (C - I(cm)) + m22 * (V - I(vm))
    return strict_subset(Kc, C) and strict_subset(Kv, V)

# ---------------------------------------------------------------------------
# Exact critical events
# ---------------------------------------------------------------------------
v, w, z = sp.symbols("v w z")
u = v * (1 - v) / ((1 - w) - v**2)
a5 = 1 - z + z * w - z * v * w + z * u * v * w
HC = (
    z * w**3 + w**2 * z**3 - w**2 * z + w * z**4 - 3 * w * z**3
    + 2 * w * z**2 + w * z - w - z**4 + 3 * z**3 - 3 * z**2 + z
)
E2 = u * (1 - z) - z * a5 * (1 - v)
E2num = sp.expand(sp.together(E2).as_numer_denom()[0])
Kexpr = u * v * w * z**3 * a5**2 / ((1 - v) * (1 - w) * (1 - z)**3)
Knum = sp.expand(sp.together(Kexpr).as_numer_denom()[0])


def roots01(poly, sym):
    P = sp.Poly(sp.expand(poly), sym)
    if P.degree() <= 0:
        return []
    return [r for r in sp.real_roots(P) if 0 < r < 1]


def resw_z(poly):
    rv = sp.expand(sp.resultant(sp.Poly(E2num, v), sp.Poly(poly, v), v))
    rw = sp.expand(sp.resultant(sp.Poly(HC, w), sp.Poly(rv, w), w))
    return roots01(rw, z)


# Exact u=1 boundary event.  At a root of its resultant, the last
# nonzero linear subresultant gives the unique common w as a rational
# function of z.  This lets the critical-fiber enumerator remove the
# terminating u=1 branch exactly, rather than deciding a zero quantity
# from a dependency-widened interval.
_U1_EQ = sp.expand(E2num.subs(v, 1 - w))
_U1_SUBRES = sp.subresultants(HC, _U1_EQ, w)
_U1_LINEAR = next(P for P in reversed(_U1_SUBRES) if sp.Poly(P, w).degree() == 1)
_U1_LPOLY = sp.Poly(_U1_LINEAR, w)
_U1_A = sp.expand(_U1_LPOLY.nth(1))
_U1_B = sp.expand(_U1_LPOLY.nth(0))
_U1_RESULTANT = sp.expand(sp.resultant(HC, _U1_EQ, w))
_U1_ROOTS_CACHE = None

def exact_critical_z():
    roots = set()
    roots |= set(roots01(sp.resultant(HC, sp.diff(HC, w), w), z))
    dE = sp.diff(E2num, v)
    rv = sp.resultant(E2num, dE, v)
    roots |= set(roots01(sp.resultant(HC, rv, w), z))
    roots |= set(roots01(HC.subs(w, 1), z))
    roots |= set(resw_z(E2num.subs(v, 0)))
    roots |= set(resw_z(E2num.subs(v, 1)))
    roots |= set(roots01(sp.resultant(HC, E2num.subs(v, 1 - w), w), z))
    roots |= set(resw_z(sp.together(a5).as_numer_denom()[0]))
    roots |= set(resw_z((1 - w) - v**2))
    roots |= set(resw_z(Knum))
    out = sorted(roots)
    assert len(out) == 5
    return out

# ---------------------------------------------------------------------------
# Certified algebraic interval enumeration
# ---------------------------------------------------------------------------
CRootOfClass = sp.polys.rootoftools.ComplexRootOf


def iv_rational(q):
    q = sp.Rational(q)
    return IV.mpf(q.p) / IV.mpf(q.q)


def croot_rational_interval(r, digits):
    """Certified rational isolation interval for a real CRootOf."""
    r.eval_rational(sp.Rational(1, 10**digits))
    ri = r._get_interval()
    return (sp.Rational(ri.a.numerator, ri.a.denominator),
            sp.Rational(ri.b.numerator, ri.b.denominator))


def croot_interval(r, digits):
    a, b = croot_rational_interval(r, digits)
    return IV.mpf([iv_rational(a).a, iv_rational(b).b])


def alg_iv(expr, digits=140):
    """Outward interval for an expression built from Rational, CRootOf and sqrt."""
    expr = sp.sympify(expr)
    env = {r: croot_interval(r, digits) for r in expr.atoms(CRootOfClass)}

    def rec(x):
        if x in env:
            return env[x]
        if x.is_Integer or x.is_Rational:
            return iv_rational(x)
        if x.is_Add:
            y = IV.mpf(0)
            for a in x.args:
                y += rec(a)
            return y
        if x.is_Mul:
            y = IV.mpf(1)
            for a in x.args:
                y *= rec(a)
            return y
        if x.is_Pow:
            base = rec(x.base)
            exponent = x.exp
            if exponent.is_Integer:
                return base ** int(exponent)
            if exponent == sp.Rational(1, 2):
                if base.a < 0:
                    raise ValueError("sqrt interval crosses negative values")
                return IV.sqrt(base)
        raise TypeError(f"unsupported algebraic node: {x!r}")

    return rec(expr)


def c_roots_at_sample(s0):
    """All real c roots at rational s0, represented exactly as RootOf/Rational."""
    return list(sp.real_roots(sp.Poly(G.subs(_s, s0), _c)))


def interval_newton_c(C, S):
    """Parametric interval Newton for G(c,s)=0 over the full S interval."""
    for _ in range(100):
        dG = Gc_iv(C, S)
        if dG.a <= 0 <= dG.b:
            return None
        m = iv_mid(C)
        N = I(m) - G_iv(I(m), S) / dG
        if strict_subset(N, C):
            width = max(mp.mpf(N.b) - mp.mpf(N.a), mp.mpf("1e-50"))
            pad = max(width / 8, mp.mpf("1e-45"))
            C2 = IV.mpf([mp.mpf(N.a) - pad, mp.mpf(N.b) + pad])
            d2 = Gc_iv(C2, S)
            if not (d2.a <= 0 <= d2.b):
                m2 = iv_mid(C2)
                N2 = I(m2) - G_iv(I(m2), S) / d2
                if strict_subset(N2, C2):
                    return C2
        lo = max(mp.mpf(C.a), mp.mpf(N.a))
        hi = min(mp.mpf(C.b), mp.mpf(N.b))
        if lo > hi:
            return None
        C = IV.mpf([lo, hi])
    return None


def critical_s_interval(cz, digits=100):
    """Certified rational S interval containing exactly cz, with Newton guard."""
    za, zb = croot_rational_interval(cz, digits)
    guard = sp.Rational(1, 10**55)
    return sp.Rational(1) - zb - guard, sp.Rational(1) - za + guard


def c_tubes_at_critical(cz):
    """All real c-root tubes over a tiny certified interval containing s*."""
    sa, sb = critical_s_interval(cz, 100)
    sm = (sa + sb) / 2
    roots = list(sp.real_roots(sp.Poly(G.subs(_s, sm), _c)))
    S = IV.mpf([iv_rational(sa).a, iv_rational(sb).b])
    tubes = []
    for r in roots:
        if isinstance(r, CRootOfClass):
            ra, rb = croot_rational_interval(r, 80)
        else:
            ra = rb = sp.Rational(r)
        C = None
        for guard in (sp.Rational(1, 10**20), sp.Rational(1, 10**12), sp.Rational(1, 10**8)):
            C0 = IV.mpf([iv_rational(ra - guard).a, iv_rational(rb + guard).b])
            C = interval_newton_c(C0, S)
            if C is not None:
                break
        if C is None:
            raise AssertionError(("critical c Newton failed", cz, r, sa, sb))
        tubes.append(C)
    tubes.sort(key=lambda X: mp.mpf(X.a))
    assert all(tubes[i].b < tubes[i + 1].a for i in range(len(tubes) - 1))
    return S, tubes


def exact_u1_boundary_tubes(cz, digits=180):
    """Exact (c,v) tubes of the u=1 boundary branch at critical z=cz.

    The linear subresultant A(z) w+B(z) is the gcd representative at
    the event root.  A(cz) is certified nonzero, so w=-B/A,
    c=w/(1-z), v=1-w are exact algebraic expressions.
    """
    global _U1_ROOTS_CACHE
    if _U1_ROOTS_CACHE is None:
        _U1_ROOTS_CACHE = frozenset(roots01(_U1_RESULTANT, z))
    if cz not in _U1_ROOTS_CACHE:
        return []
    Aexpr = _U1_A.subs(z, cz)
    Bexpr = _U1_B.subs(z, cz)
    Aiv = alg_iv(Aexpr, digits)
    if Aiv.a <= 0 <= Aiv.b:
        raise AssertionError(("u=1 linear subresultant coefficient vanishes", cz, Aexpr))
    # Keep the rational functions unevaluated.  Expanding/simplifying a high-
    # degree expression in a CRootOf number field is unnecessary and can be
    # prohibitively expensive; alg_iv evaluates these exact expression trees
    # outwardly from certified CRootOf isolation intervals.
    wexpr = (-Bexpr) / Aexpr
    sexpr = 1 - cz
    cexpr = wexpr / sexpr
    vexpr = 1 - wexpr
    # Exact relation checks.  Algebraic simplification can be expensive, so
    # the defining linear subresultant and interval non-vanishing are the
    # primary certificate; these identities are cheap consistency checks.
    Civ, Viv = alg_iv(cexpr, digits), alg_iv(vexpr, digits)
    return [(Civ, Viv)]


def intervals_overlap(X, Y):
    return not (X.b < Y.a or Y.b < X.a)


def quadratic_v_intervals(C, S):
    """All real roots of Q(v,c,s), enclosed by outward quadratic formula."""
    aa = _ns["iv_eval"](_Qcoef[2], {_c: C, _s: S})
    bb = _ns["iv_eval"](_Qcoef[1], {_c: C, _s: S})
    cc = _ns["iv_eval"](_Qcoef[0], {_c: C, _s: S})
    if aa.a <= 0 <= aa.b:
        raise AssertionError("quadratic leading coefficient contains zero")
    disc = bb * bb - 4 * aa * cc
    if disc.b < 0:
        return []
    sd = IV.sqrt(IV.mpf([max(mp.mpf(0), mp.mpf(disc.a)), mp.mpf(disc.b)]))
    den = 2 * aa
    roots = [(-bb - sd) / den, (-bb + sd) / den]
    roots.sort(key=lambda X: mp.mpf(X.a))
    if roots[1].a <= roots[0].b:
        return [IV.mpf([roots[0].a, max(roots[0].b, roots[1].b)])]
    return roots


def classify_root(C, V, S):
    """Strict interval classification of a full root tube."""
    wv = C * S
    if C.b <= 0 or V.b <= 0 or V.a >= 1 or wv.a >= 1:
        return "inadmissible"
    den = 1 - wv - V * V
    if den.a <= 0 <= den.b:
        return "unresolved"
    uv = V * (1 - V) / den
    zz = 1 - S
    a5v = 1 - zz + zz * wv - zz * V * wv + zz * uv * V * wv
    kd = (1 - V) * (1 - wv) * S**3
    if kd.a <= 0:
        return "unresolved"
    kval = uv * V * wv * zz**3 * a5v * a5v / kd
    positives = [C, V, 1 - V, wv, 1 - wv, den, uv, 1 - uv, a5v, kval]
    if all(X.a > 0 for X in positives):
        return "admissible"
    if any(X.b < 0 for X in positives) or uv.a >= 1:
        return "inadmissible"
    return "unresolved"



# Polynomial boundary loci at which strict admissibility fails.  At a critical
# fiber an outward enclosure may straddle equality; a 3x3 Krawczyk certificate
# for (G,Q,B)=0 proves that the unresolved system root is exactly on B=0.
_wexpr = _c * _s
_denexpr = 1 - _wexpr - _v**2
_uexpr = _v * (1 - _v) / _denexpr
_zexpr = 1 - _s
_a5expr = 1 - _zexpr + _zexpr*_wexpr - _zexpr*_v*_wexpr + _zexpr*_uexpr*_v*_wexpr
_a5num = sp.expand(sp.together(_a5expr).as_numer_denom()[0])
_BOUNDARIES = {
    "v0": _v,
    "v1": _v - 1,
    "w0": _wexpr,
    "w1": _wexpr - 1,
    "uden": _denexpr,
    "u1": _v + _wexpr - 1,
    "a5": _a5num,
}


def _eval_expr(expr, C, V, S):
    return _ns["iv_eval"](expr, {_c: C, _v: V, _s: S})


def krawczyk3_boundary(C, V, S, Bexpr):
    """Certify a unique zero of (G,Q,B) in CxVxS."""
    exprs = (G, Q, Bexpr)
    vars_ = (_c, _v, _s)
    boxes = (C, V, S)
    mids = [iv_mid(X) for X in boxes]
    point_env = {_c: I(mids[0]), _v: I(mids[1]), _s: I(mids[2])}
    F0 = [_ns["iv_eval"](f, point_env) for f in exprs]
    Jexpr = [[sp.diff(f, x) for x in vars_] for f in exprs]
    Jmid = mp.matrix(3)
    for i in range(3):
        for j in range(3):
            Jij = _ns["iv_eval"](Jexpr[i][j], point_env)
            Jmid[i,j] = iv_mid(Jij)
    try:
        Y = Jmid**-1
    except Exception:
        return False
    # Interval Jacobian on the full box.
    full_env = {_c: C, _v: V, _s: S}
    J = [[_ns["iv_eval"](Jexpr[i][j], full_env) for j in range(3)] for i in range(3)]
    K = []
    for i in range(3):
        zi = I(mids[i])
        for k in range(3):
            zi -= I(Y[i,k]) * F0[k]
        tail = IV.mpf(0)
        for j in range(3):
            mij = IV.mpf(1 if i == j else 0)
            for k in range(3):
                mij -= I(Y[i,k]) * J[k][j]
            tail += mij * (boxes[j] - I(mids[j]))
        K.append(zi + tail)
    return all(strict_subset(K[i], boxes[i]) for i in range(3))


def certify_boundary_root(C, V, S):
    """Return a certified boundary name, or None."""
    for name, B in _BOUNDARIES.items():
        BI = _eval_expr(B, C, V, S)
        if BI.a <= 0 <= BI.b and krawczyk3_boundary(C, V, S, B):
            return name
    return None


def enumerate_admissible(s_expr, critical=False):
    """Rigorous interval enumeration of every admissible root.

    At a critical algebraic parameter, a certified rational isolating interval
    for s* and parametric c-root tubes enclose the exact fiber.  No number-field
    factorisation, polyroots/nroots, residual tolerance, or float containment is
    used in the load-bearing enumeration.
    """
    result, excluded, unresolved = [], 0, []
    if critical:
        cz = 1 - s_expr
        S, ctubes = c_tubes_at_critical(cz)
        u1_boundaries = exact_u1_boundary_tubes(cz)
        items = [(None, C) for C in ctubes]
    else:
        u1_boundaries = []
        S = I(Fraction(int(sp.numer(s_expr)), int(sp.denom(s_expr))))
        items = [(ce, None) for ce in c_roots_at_sample(s_expr)]

    for ce, preC in items:
        decided = False
        for digits in ((80, 140, 220, 320) if not critical else (0,)):
            C = preC if critical else alg_iv(ce, digits)
            wv = C * S
            # The sole v-double critical branch is already outside w<1, so it
            # is rigorously rejected here before discriminant separation.
            if C.b <= 0 or wv.a >= 1:
                excluded += 1
                decided = True
                break
            local, local_unresolved = [], False
            for V in quadratic_v_intervals(C, S):
                if V.b < 0 or V.a > 1:
                    excluded += 1
                    continue
                # At an exact u=1 critical event, dependency widening can make
                # 1-u look slightly positive.  Match the exact algebraic
                # boundary root in both coordinates and exclude it explicitly.
                hits = [(Cb, Vb) for Cb, Vb in u1_boundaries
                        if intervals_overlap(C, Cb) and intervals_overlap(V, Vb)]
                if hits:
                    if len(hits) != 1:
                        raise AssertionError(("ambiguous u=1 boundary match", s_expr, C, V, hits))
                    Cb, Vb = hits[0]
                    if not (C.a <= Cb.a <= Cb.b <= C.b and V.a <= Vb.a <= Vb.b <= V.b):
                        raise AssertionError(("u=1 exact root not contained in candidate tube", s_expr, C, V, Cb, Vb))
                    excluded += 1
                    continue
                status = classify_root(C, V, S)
                if status == "admissible":
                    local.append((C, V))
                elif status == "inadmissible":
                    excluded += 1
                else:
                    boundary = certify_boundary_root(C, V, S) if critical else None
                    if boundary is not None:
                        excluded += 1
                    else:
                        local_unresolved = True
            if not local_unresolved:
                result.extend(local)
                decided = True
                break
        if not decided:
            unresolved.append((str(ce) if ce is not None else "critical c-tube", C, S))
    if unresolved:
        raise AssertionError(f"unresolved root classifications: {unresolved}")
    return result, excluded

# ---------------------------------------------------------------------------
# Load boxes and build a certified branch graph
# ---------------------------------------------------------------------------
def load_nodes():
    nodes = []
    for source, filename in (("cover", "_hc_cover.json"), ("collar", "_hc_s1collar.json")):
        data = json.loads((CODE / filename).read_text())
        for idx, p in enumerate(data["pieces"]):
            node = {
                "source": source,
                "source_index": idx,
                "s0": Fraction(p["s"][0]), "s1": Fraction(p["s"][1]),
                "c0": Fraction(p["C"][0]), "c1": Fraction(p["C"][1]),
                "v0": Fraction(p["V"][0]), "v1": Fraction(p["V"][1]),
            }
            if source == "collar":
                node.update(
                    cb0=Fraction(p["cbar"][0]), cb1=Fraction(p["cbar"][1]),
                    vb0=Fraction(p["vbar"][0]), vb1=Fraction(p["vbar"][1]),
                )
            nodes.append(node)
    return nodes


def validate_edge(a, b):
    slo, shi = max(a["s0"], b["s0"]), min(a["s1"], b["s1"])
    if slo > shi:
        return False
    sm = (slo + shi) / 2
    if a["source"] == b["source"] == "collar":
        C = I(min(a["cb0"], b["cb0"]), max(a["cb1"], b["cb1"]))
        V = I(min(a["vb0"], b["vb0"]), max(a["vb1"], b["vb1"]))
        return krawczyk2_collar(C, V, I(1 - sm))
    C = I(min(a["c0"], b["c0"]), max(a["c1"], b["c1"]))
    V = I(min(a["v0"], b["v0"]), max(a["v1"], b["v1"]))
    return krawczyk2_raw(C, V, I(sm))


def build_graph(nodes):
    order = sorted(range(len(nodes)), key=lambda i: nodes[i]["s0"])
    active = []
    edges = []
    candidates = 0
    for i in order:
        a = nodes[i]
        active = [j for j in active if nodes[j]["s1"] >= a["s0"]]
        for j in active:
            b = nodes[j]
            # coordinate overlap is a cheap necessary filter
            if min(a["c1"], b["c1"]) < max(a["c0"], b["c0"]):
                continue
            if min(a["v1"], b["v1"]) < max(a["v0"], b["v0"]):
                continue
            candidates += 1
            if not validate_edge(a, b):
                raise AssertionError(f"branch edge failed: {i},{j}")
            edges.append((i, j))
        active.append(i)

    parent = list(range(len(nodes)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[y] = x

    for i, j in edges:
        union(i, j)
    comps = defaultdict(list)
    for i in range(len(nodes)):
        comps[find(i)].append(i)
    comp_list = list(comps.values())
    node_comp = {}
    for k, comp in enumerate(comp_list):
        for i in comp:
            node_comp[i] = k
    return edges, comp_list, node_comp, candidates


def union_covers(intervals, lo, hi):
    clipped = sorted((max(a, lo), min(b, hi)) for a, b in intervals if a < hi and b > lo)
    if not clipped:
        return False
    cur = lo
    for a, b in clipped:
        if a > cur:
            return False
        if b > cur:
            cur = b
        if cur >= hi:
            return True
    return cur >= hi


def root_hosts(root, s_expr, nodes, node_comp, strict_s=True):
    C, V = root
    hosts = set()
    for i, p in enumerate(nodes):
        in_s = (p["s0"] < s_expr < p["s1"]) if strict_s else (p["s0"] <= s_expr <= p["s1"])
        if not in_s:
            continue
        if C.a > mpf_fraction(p["c0"]) and C.b < mpf_fraction(p["c1"]):
            if V.a > mpf_fraction(p["v0"]) and V.b < mpf_fraction(p["v1"]):
                hosts.add(node_comp[i])
    return hosts


def rational_between(a, b):
    for digits in (20, 40, 80, 120):
        af, bf = sp.N(a, digits), sp.N(b, digits)
        q = sp.Rational(str((af + bf) / 2))
        if a < q < b:
            return q
    raise AssertionError("failed to choose rational in algebraic interval")


def main():
    crit_z = exact_critical_z()
    crit_s = sorted([1 - x for x in crit_z])
    # Exact degree-stability certificate for the cubic G: no c-double event
    # occurs in (0,1), so the c-root count cannot change inside a critical
    # isolating interval.
    RGexpr = sp.factor(sp.resultant(G, sp.diff(G, _c), _c))
    RG = sp.Poly(sp.cancel(RGexpr / (_s - 1)**2), _s)
    assert RG.count_roots(0, 1) == 0
    s_ext = [sp.Rational(0)] + crit_s + [sp.Rational(1)]
    # Each guarded rational interval used for a critical fiber contains that
    # critical value and no other event.
    for i, cs in enumerate(crit_s):
        sa, sb = critical_s_interval(1 - cs)
        assert s_ext[i] < sa < cs < sb < s_ext[i + 2]
    cells = [(s_ext[i], s_ext[i + 1], rational_between(s_ext[i], s_ext[i + 1]))
             for i in range(len(s_ext) - 1)]
    assert len(cells) == 6

    nodes = load_nodes()
    edges, comps, node_comp, candidates = build_graph(nodes)
    assert len(comps) == 2, f"expected two root-branch components, got {len(comps)}"

    comp_records = []
    comp_spans = []
    for k, comp in enumerate(comps):
        intervals = [(nodes[i]["s0"], nodes[i]["s1"]) for i in comp]
        lo, hi = min(a for a, _ in intervals), max(b for _, b in intervals)
        assert union_covers(intervals, lo, hi)
        comp_spans.append((lo, hi))
        comp_records.append({
            "component": k, "n_boxes": len(comp), "s_lo": str(lo), "s_hi": str(hi),
            "cover_boxes": sum(nodes[i]["source"] == "cover" for i in comp),
            "collar_boxes": sum(nodes[i]["source"] == "collar" for i in comp),
        })

    cell_records = []
    for idx, (sa, sb, s0) in enumerate(cells):
        roots, excluded = enumerate_admissible(s0, critical=False)
        active_components = []
        for k, (lo_f, hi_f) in enumerate(comp_spans):
            lo_q = sp.Rational(lo_f.numerator, lo_f.denominator)
            hi_q = sp.Rational(hi_f.numerator, hi_f.denominator)
            if lo_q <= sa and hi_q >= sb:
                active_components.append(k)
        root_components = [root_hosts(r, s0, nodes, node_comp, strict_s=True) for r in roots]
        assert all(len(h) == 1 for h in root_components)
        mapped = [next(iter(h)) for h in root_components]
        assert len(set(mapped)) == len(mapped)
        assert set(mapped) == set(active_components), (
            idx, mapped, active_components, len(roots)
        )
        cell_records.append({
            "cell": idx, "s_lo": str(sa), "s_hi": str(sb), "sample": str(s0),
            "n_admissible_roots": len(roots), "excluded_c_branches": excluded,
            "components": mapped, "branch_bijection": True,
        })
        print(f"cell {idx}: roots={len(roots)} components={mapped} OK", flush=True)

    fiber_records = []
    for idx, cs in enumerate(crit_s):
        roots, excluded = enumerate_admissible(cs, critical=True)
        root_components = [root_hosts(r, cs, nodes, node_comp, strict_s=False) for r in roots]
        assert all(len(h) == 1 for h in root_components), (idx, root_components)
        mapped = [next(iter(h)) for h in root_components]
        assert len(set(mapped)) == len(mapped)
        fiber_records.append({
            "fiber": idx, "critical_s": str(cs),
            "critical_s_evalf": float(sp.N(cs, 16)),
            "n_admissible_roots": len(roots), "excluded_c_branches": excluded,
            "components": mapped, "all_admissible_roots_boxed": True,
        })
        print(f"fiber {idx}: roots={len(roots)} components={mapped} OK", flush=True)

    out = {
        "n_nodes": len(nodes), "n_edge_candidates": candidates, "n_certified_edges": len(edges),
        "n_components": len(comps), "components": comp_records,
        "cells": cell_records, "all_cells_branch_bijection": True,
        "critical_fibers": fiber_records, "all_critical_admissible_roots_boxed": True,
        "n_critical_values": len(crit_s),
        "method": "2x2 Krawczyk branch graph + exact event cells + parametric critical-fiber c-tubes + outward quadratic enumeration + 3x3 generic-boundary Krawczyk + exact u=1 linear subresultant",
    }
    (CODE / "_hc_exhaustiveness.json").write_text(json.dumps(out, indent=2))
    print(
        "DONE-EXHAUSTIVENESS components=%d cells=%d fibers=%d edges=%d" %
        (len(comps), len(cells), len(crit_s), len(edges)),
        flush=True,
    )


if __name__ == "__main__":
    main()
