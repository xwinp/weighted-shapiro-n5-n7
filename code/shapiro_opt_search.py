#!/usr/bin/env python3
"""
Multi-start numerical optimization for the Tuan-Thuong weighted Shapiro sum

    P_{n,p,q}(x) = sum_i x_i / (p*x_{i+1} + q*x_{i+2}),   q = 1 - p

Homogeneous of degree 0, so we minimize over the simplex sum_i x_i = 1, x_i>=0.
Claimed lower bound is n/(p+q) = n.

Track 1 (n=7 critical curve): for a grid of p, find min_x P; the TRUE failure
region is {p : min_x P < 7}. Our single-witness interval (0.21962, 0.30630) is
only a lower bound on this region.

Track 2 (existence at p=0.3 for n=9,11,13): if min_x P < n, a counterexample
exists (then we rationalize the argmin and exact-verify).
"""
import numpy as np
from scipy.optimize import minimize
from fractions import Fraction as F
from math import gcd
import multiprocessing as mp


def P_float(x, p, n):
    q = 1.0 - p
    # Valid domain requires every denominator > 0. Exclude degenerate 0/0.
    for i in range(n):
        d = p * x[(i + 1) % n] + q * x[(i + 2) % n]
        if d <= 1e-15:
            return 1e18
    s = 0.0
    for i in range(n):
        d = p * x[(i + 1) % n] + q * x[(i + 2) % n]
        s += x[i] / d
    return s


def min_P(n, p, n_starts=400, seed=0):
    rng = np.random.default_rng(seed)
    cons = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1.0},)
    bounds = [(0.0, 1.0)] * n
    best = 1e18
    bestx = None
    for s in range(n_starts):
        x0 = rng.random(n)
        x0 /= x0.sum()
        try:
            res = minimize(P_float, x0, args=(p, n), method='SLSQP',
                           bounds=bounds, constraints=cons,
                           options={'maxiter': 300, 'ftol': 1e-12})
            if res.fun < best:
                best = res.fun
                bestx = res.x
        except Exception:
            pass
    return best, bestx


def rationalize(x, denom=1000):
    """Round floats to small rationals; return primitive int vector (scale-free)."""
    xs = np.round(np.array(x) * denom).astype(int)
    # try a few denominators to get a clean small vector
    best = xs
    for d in (10, 20, 50, 100, 200, 500, 1000, 2000):
        v = np.round(np.array(x) * d).astype(int)
        nz = v[v != 0]
        if len(nz) == 0:
            continue
        g = nz[0]
        for u in nz[1:]:
            g = gcd(g, u)
        if g > 1:
            v = v // g
        if v.max() < best.max() and v.max() > 0:
            best = v
    return tuple(int(v) for v in best)


def exact_P(xv, p, q, n):
    s = F(0)
    for i in range(n):
        d = p * xv[(i + 1) % n] + q * xv[(i + 2) % n]
        if d == 0:
            if xv[i] == 0:
                continue
            return None
        s += F(xv[i]) / d
    return s


# ---------- Track 1: n=7 critical curve ----------
def n7_point(p):
    val, x = min_P(7, p, n_starts=600, seed=12345 + int(p * 10000))
    return p, val, x


def track1_n7():
    print("=== Track 1: n=7 critical curve (min_x P vs p) ===", flush=True)
    ps = np.arange(0.02, 0.98, 0.01)
    with mp.Pool(processes=min(16, mp.cpu_count())) as pool:
        res = pool.map(n7_point, ps)
    res.sort()
    print(f"\n{'p':>8} {'minP':>12} {'<7?':>5}", flush=True)
    crossings = []
    prev_fail = None
    for p, val, x in res:
        fail = val < 7 - 1e-7
        print(f"{p:8.4f} {val:12.6f} {'YES' if fail else 'no':>5}", flush=True)
        if prev_fail is not None and fail != prev_fail:
            crossings.append(p)
        prev_fail = fail
    print(f"\ncrossings near p ~ {crossings}", flush=True)
    # refine: for the apparent failure interval, find best witness and rationalize
    fail_pts = [(p, val, x) for p, val, x in res if val < 7 - 1e-7]
    if fail_pts:
        pmin = min(p for p, _, _ in fail_pts)
        pmax = max(p for p, _, _ in fail_pts)
        print(f"apparent numerical failure interval: [{pmin:.4f}, {pmax:.4f}]", flush=True)
        print(f"(single-witness rigorous interval was [0.2196, 0.3063])", flush=True)
        # rationalize the strongest witness
        best = min(fail_pts, key=lambda t: t[1])
        p, val, x = best
        xv = rationalize(x)
        ex = exact_P(xv, F(3, 10), F(7, 10), 7)
        print(f"strongest witness at p={p:.4f}: rationalized x={xv}, exact P(0.3)={ex}", flush=True)


# ---------- Track 2: existence at p=0.3 for n=9,11,13 ----------
def track2_existence():
    print("\n=== Track 2: existence at p=0.3 (q=0.7) ===", flush=True)
    for n in (9, 11, 13):
        val, x = min_P(n, 0.3, n_starts=1500, seed=99999 + n)
        fail = val < n - 1e-7
        print(f"n={n}: minP ~ {val:.6f}  bound={n}  cex? {'YES' if fail else 'no'}", flush=True)
        if fail:
            xv = rationalize(x, denom=2000)
            ex = exact_P(xv, F(3, 10), F(7, 10), n)
            print(f"   rationalized x={xv}", flush=True)
            print(f"   exact P = {ex} ({float(ex) if ex else 'undef'})  < {n}? "
                  f"{ex is not None and ex < n}", flush=True)


if __name__ == "__main__":
    track1_n7()
    track2_existence()
