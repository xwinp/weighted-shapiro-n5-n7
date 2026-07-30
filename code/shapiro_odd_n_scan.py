#!/usr/bin/env python3
"""
Extend the counterexample search to all odd n=7,9,11,13,15,17,19,21 at p=0.3,q=0.7.

Two prongs:
 (A) Test a generalized alternating-zero pattern family explicitly.
 (B) Small-box exhaustive search (B=1..3) for each n to catch any witness.
All hits are exact-verified with fractions.Fraction, then the p-interval
(over q=1-p, same vector) is derived by solving P(p)=n symbolically.
"""
import numpy as np
from fractions import Fraction as F
from sympy import symbols, simplify, Poly, real_roots, Rational as sR, sqrt
from math import gcd

P, Q = F(3, 10), F(7, 10)
n_sym = None


def exact_P(xv, p=F(3, 10), q=F(7, 10)):
    n = len(xv)
    s = F(0)
    for i in range(n):
        d = p * xv[(i + 1) % n] + q * xv[(i + 2) % n]
        if d == 0:
            if xv[i] == 0:
                continue
            return None
        s += F(xv[i]) / d
    return s


def valid_domain(xv):
    n = len(xv)
    for i in range(n):
        d = 3 * xv[(i + 1) % n] + 7 * xv[(i + 2) % n]
        if d <= 0:
            return False
    return True


def box_search(n, Bmax):
    for B in range(1, Bmax + 1):
        k = B + 1
        rng = [range(k)] * n
        # iterate via numpy grid only if small enough
        M = k ** n
        if M > 6_000_000:
            print(f"  [box] n={n} B={B} too large ({M}), skipping", flush=True)
            continue
        grid = np.stack(np.meshgrid(*[np.arange(k)] * n, indexing='ij')).reshape(n, -1).T
        x = grid.astype(np.float64)
        x1 = np.roll(x, -1, axis=1); x2 = np.roll(x, -2, axis=1)
        den = 0.3 * x1 + 0.7 * x2
        invalid = (den <= 0).any(axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            term = np.where(den > 0, x / den, 0.0)
        Psum = term.sum(axis=1); Psum[invalid] = np.inf
        cand = np.where(Psum < n - 1e-9)[0]
        if len(cand) == 0:
            continue
        best = None
        for idx in cand:
            xv = tuple(int(v) for v in grid[idx])
            if not valid_domain(xv):
                continue
            s = exact_P(xv)
            if s is not None and s < n:
                if best is None or xv < best[0]:
                    best = (xv, s)
        if best:
            return B, best[0], best[1]
    return None, None, None


def p_interval(xv):
    """For witness xv (q=1-p), solve P(p)=n symbolically; return (p_lo,p_hi) interval where P<n."""
    n = len(xv)
    p = symbols('p')
    expr = 0
    for i in range(n):
        num = xv[i]
        den = p * xv[(i + 1) % n] + (1 - p) * xv[(i + 2) % n]
        expr += num / den
    diff = simplify(expr - n)
    numP = simplify(diff.as_numer_denom()[0])
    roots = real_roots(Poly(numP, p))
    pos = sorted(float(r) for r in roots if 0 < float(r) < 1)
    # find an interval (a,b) where P-n<0 containing p=0.3
    # sample sign between consecutive roots
    cex_intervals = []
    test_pts = [0.0] + [(pos[i] + pos[i + 1]) / 2 for i in range(len(pos) - 1)] + [1.0]
    bounds = [0.0] + pos + [1.0]
    for i in range(len(bounds) - 1):
        mid = (bounds[i] + bounds[i + 1]) / 2
        val = float(diff.subs(p, mid))
        if val < -1e-12:
            cex_intervals.append((bounds[i], bounds[i + 1]))
    return cex_intervals, pos


def pattern_witness(n):
    """Generalized alternating pattern: even positions 0, odd positions 1, last = c.
       Try a few c values."""
    cands = []
    for c in range(0, 6):
        xv = [0] * n
        for j in range(n):
            if j % 2 == 1:
                xv[j] = 1
        xv[-1] = c
        if valid_domain(xv) and sum(xv) > 0:
            cands.append(tuple(xv))
    return cands


if __name__ == "__main__":
    print("=== Odd-n scan at p=0.3, q=0.7 ===\n", flush=True)
    for n in [7, 9, 11, 13, 15, 17, 19, 21]:
        print(f"--- n={n} ---", flush=True)
        # prong A: pattern
        hit = None
        for xv in pattern_witness(n):
            s = exact_P(xv)
            if s is not None and s < n:
                hit = (xv, s)
                break
        if hit:
            xv, s = hit
            print(f"  pattern witness: x={xv}", flush=True)
            print(f"  P={s}={float(s):.6f} < {n}  (deficit {n-s})", flush=True)
        else:
            print(f"  pattern: no cex among tried c", flush=True)
        # prong B: box
        B, bxv, bs = box_search(n, 3)
        if bxv:
            print(f"  box B={B}: x={bxv}", flush=True)
            print(f"  P={bs}={float(bs):.6f} < {n}", flush=True)
            # use the smaller-deficit witness for interval (prefer box if found)
            chosen = bxv if (hit is None or bs < hit[1] - 1e-12) else hit[0]
            chosen_s = bs if chosen is bxv else hit[1]
        elif hit:
            chosen = hit[0]; chosen_s = hit[1]
        else:
            print(f"  NO counterexample found for n={n}\n", flush=True)
            continue
        # derive p-interval
        cex_int, pos = p_interval(chosen)
        print(f"  p-interval (q=1-p) where this vector is cex: {cex_int}", flush=True)
        print(f"  roots: {pos}\n", flush=True)
