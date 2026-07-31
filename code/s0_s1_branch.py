#!/usr/bin/env python3
"""
# !! DIAGONAL BUG (2026-07-29): the 4-point mixed-difference H_ij below degenerates
# to f'/h noise for i=j (diagonal). Diagonal entries are WRONG. Use code/sym_hess3.py
# (symbolic) or code/det_scan2.py (three-point diagonal) instead. Off-diagonal entries are OK.
# Corrected result: S1 Morse=1, S0-nonunif Morse=1 (not 3/4); cert target is det(H)<0.
Dense continuation of S0 (full support) and S1 ({0}) positive-stationary branches
across the failure band (a7,b7), starting from the GPT p=1/4 inits.
Goal: confirm (numerically) that EVERY S0 non-uniform stationary point and EVERY
S1 stationary point in the band is a SADDLE (reduced Hessian has a negative
eigenvalue). This completes the Morse-index picture for both remaining gaps.

Strategy: tiny p-steps; each solve uses the previous solution as init; on failure,
halve the step and try again (up to a few retries); also try random perturbations.
"""
import numpy as np
from scipy.optimize import root
n = 7

def Pval(x, p):
    q = 1 - p; s = 0.0
    for i in range(n):
        den = p*x[(i+1)%n] + q*x[(i+2)%n]
        if abs(den) < 1e-15: return 1e6
        s += x[i]/den
    return s

def hessian_full(x, p):
    H = np.zeros((n, n)); h = 1e-5
    for i in range(n):
        for j in range(n):
            xp=x.copy(); xm=x.copy(); xpp=x.copy(); xmm=x.copy()
            xp[i]+=h; xp[j]+=h; xm[i]-=h; xm[j]-=h
            xpp[i]+=h; xpp[j]-=h; xmm[i]-=h; xmm[j]+=h
            H[i,j] = (Pval(xp,p)-Pval(xm,p)-Pval(xpp,p)+Pval(xmm,p))/(4*h*h)
    return H

def tangent_basis(free_indices):
    k = len(free_indices); B = np.zeros((n, k-1))
    for j in range(k-1):
        B[free_indices[j], j] = 1
        B[free_indices[-1], j] = -1
    return B

def reduced_eigs(x, p, free_indices):
    H = hessian_full(x, p); B = tangent_basis(free_indices)
    return np.linalg.eigvalsh(B.T @ H @ B)

def make_solver(free, anchor_idx, anchor_val):
    """grad=0 on coords `free` (relative to a normalized anchor)."""
    def grad(v):
        x = np.zeros(n); x[anchor_idx] = anchor_val
        for j, idx in enumerate(free): x[idx] = v[j]
        h = 1e-7; f0 = Pval(x, p_cur[0]); g = []
        for j in range(len(free)):
            xp = x.copy(); xp[free[j]] += h
            g.append((Pval(xp, p_cur[0]) - f0)/h)
        return g
    return grad

p_cur = [0.25]  # mutable so grad closure sees current p

def solve_branch(free, anchor_idx, p0, v0, p_list, label):
    """Continue from p0 with init v0 across p_list (sorted toward p0)."""
    grad = make_solver(free, anchor_idx, 1.0)
    results = []
    v = list(v0)
    # build path from p0 outward
    for pp in p_list:
        p_cur[0] = pp
        ok = False
        tries = [v]
        # add a few perturbations
        for _ in range(6):
            tries.append([vi*(1+0.05*(np.random.random()-0.5)*2) for vi in v])
        for init in tries:
            r = root(grad, init, method='hybr', options={'xtol':1e-14, 'maxfev':20000})
            if r.success and max(abs(r.fun))<1e-7 and all(ri>1e-7 for ri in r.x):
                v = list(r.x); ok = True; break
        if not ok:
            results.append((pp, None)); continue
        x = np.zeros(n); x[anchor_idx] = 1.0
        for j, idx in enumerate(free): x[idx] = v[j]
        xs = x/x.sum(); Pv = Pval(xs, pp)
        eigs = reduced_eigs(xs, pp, list(range(n)) if anchor_idx!=0 or len(free)==6
                            else [i for i in range(n) if i!=0])
        # proper free_indices for tangent: the support indices
        if anchor_idx == 0 and len(free) == 6:  # S0, support all
            fi = list(range(n))
        elif anchor_idx == 1 and len(free) == 5:  # S1, support {1..6}
            fi = [1,2,3,4,5,6]
        else:
            fi = free + [anchor_idx]
        eigs = reduced_eigs(xs, pp, fi)
        neg = int((eigs < -1e-4).sum())
        results.append((pp, (Pv, eigs, neg, xs)))
        print(f"  [{label}] p={pp:.5f} P={Pv:.6f} (P-7={Pv-7:+.5f}) neg_eig={neg} -> {'SADDLE' if neg>0 else 'MIN?'}")
    return results

if __name__ == '__main__':
    np.random.seed(1)
    a7, b7 = 0.214273520909841, 0.328627677916592
    # ---- S1 ----
    print("=== S1 branch ({0}) dense continuation from p=1/4 ===")
    s1_free = [2,3,4,5,6]; s1_anchor = 1
    s1_v0 = [0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
    # upward then downward from 0.25
    up = list(np.linspace(0.25, b7-0.003, 18))
    down = list(np.linspace(0.25, a7+0.003, 18))[::-1]
    solve_branch(s1_free, s1_anchor, 0.25, s1_v0, up[1:], "S1↑")
    solve_branch(s1_free, s1_anchor, 0.25, s1_v0, down[1:], "S1↓")

    # ---- S0 ----
    print("\n=== S0 branch (full support) dense continuation from p=1/4 ===")
    s0_free = [1,2,3,4,5,6]; s0_anchor = 0
    s0_v0 = [0.2899598915706492,0.7447488449365604,0.8875842694414185,0.3764381331194338,1.1665179227192974,0.1591166088535238]
    up0 = list(np.linspace(0.25, b7-0.003, 18))
    down0 = list(np.linspace(0.25, a7+0.003, 18))[::-1]
    solve_branch(s0_free, s0_anchor, 0.25, s0_v0, up0[1:], "S0↑")
    solve_branch(s0_free, s0_anchor, 0.25, s0_v0, down0[1:], "S0↓")
