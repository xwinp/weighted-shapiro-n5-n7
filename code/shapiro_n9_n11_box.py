#!/usr/bin/env python3
"""
Exhaustive small-integer-box search for counterexamples to the Tuan-Thuong
weighted Shapiro inequality

    P_{n,p,q}(x) = sum_i x_i / (p*x_{i+1} + q*x_{i+2})  >=  n/(p+q)

at p = 3/10, q = 7/10  (so the claimed lower bound is n).

Homogeneous of degree 0 => if ANY integer vector is a counterexample, its
primitive reduction is too and lives in a smaller box, so we do NOT need a
primitivity filter for existence: just find any x with P < n.

Vectorized over the whole box with numpy float64 prescreen; survivors are
exact-verified with fractions.Fraction.
"""
import sys, time
import numpy as np
from fractions import Fraction as F

P = F(3, 10)
Q = F(7, 10)
PF, QF = 3 / 10, 7 / 10


def box_search(n, Bmax):
    """Exhaustive x_i in 0..B for B=1..Bmax. Return first cex per B (lex)."""
    for B in range(1, Bmax + 1):
        t0 = time.time()
        k = B + 1
        M = k ** n
        print(f"[n={n}] B={B}  grid={k}^{n} = {M:,} vectors ...", flush=True)
        # build grid (M, n) int8
        ranges = [np.arange(k, dtype=np.int16)] * n
        grid = np.stack(np.meshgrid(*ranges, indexing='ij')).reshape(n, -1).T
        x = grid.astype(np.float64)
        x1 = np.roll(x, -1, axis=1)
        x2 = np.roll(x, -2, axis=1)
        den = PF * x1 + QF * x2
        # Valid domain: every denominator must be strictly positive (no two
        # consecutive cyclic zeros). Exclude any vector with a zero denom.
        invalid = (den <= 0).any(axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            term = np.where(den > 0, x / den, 0.0)
        Psum = term.sum(axis=1)
        Psum[invalid] = np.inf
        # prescreen: strictly below n - tiny epsilon
        cand = np.where(Psum < n - 1e-9)[0]
        dt = time.time() - t0
        print(f"[n={n}] B={B}  prescreen survivors={len(cand)}  ({dt:.1f}s)", flush=True)
        if len(cand) == 0:
            continue
        # exact-verify survivors, collect lex-first
        best = None
        for idx in cand:
            xv = grid[idx].tolist()
            s = F(0)
            ok = True
            for i in range(n):
                d = P * xv[(i + 1) % n] + Q * xv[(i + 2) % n]
                if d == 0:
                    if xv[i] == 0:
                        continue
                    ok = False
                    break
                s += F(xv[i]) / d
            if ok and s < n:
                if best is None or tuple(xv) < best[0]:
                    best = (tuple(xv), s)
        if best is not None:
            xv, s = best
            print(f"[n={n}] B={B}  COUNTEREXAMPLE FOUND", flush=True)
            print(f"        x = {xv}", flush=True)
            print(f"        P = {s} = {float(s):.10f}  <  {n}   (deficit {n - s})", flush=True)
            return n, B, xv, s
    print(f"[n={n}] no counterexample up to B={Bmax}", flush=True)
    return n, None, None, None


if __name__ == "__main__":
    n9_Bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    n11_Bmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    results = []
    print("=== n=9 ===", flush=True)
    results.append(box_search(9, n9_Bmax))
    print("=== n=11 ===", flush=True)
    results.append(box_search(11, n11_Bmax))
    print("\n=== SUMMARY ===", flush=True)
    for n, B, xv, s in results:
        if xv is None:
            print(f"n={n}: no cex (Bmax searched)")
        else:
            print(f"n={n}: cex at B={B}, x={xv}, P={s} ({float(s):.6f}) < {n}")
