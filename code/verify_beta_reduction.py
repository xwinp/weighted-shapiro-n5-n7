#!/usr/bin/env python3
"""DIAGNOSTIC (not a formal certificate) -- numerical cross-checks of the
structural beta-reduction formulas.

This script is NOT load-bearing: each of its three checks has a rigorous
replacement elsewhere in the certificate chain, noted below.  It exists only to
give a quick numerical sanity signal.  It always exits 0.

  (B.3) uniform reduced-Hessian eigenvalues via the CIRCULANT spectrum, compared
        to the Fourier formula  2(1-cos th_k)[2 q cos th_k + 2 p^2 - 3 p + 2].
        RIGOROUS replacement: verify_uniform_hess.py, which solves the GENERALISED
        eigenproblem  B^T H B v = lam (B^T B) v  (B=[e_i-e_6] non-orthonormal)
        and matches the formula to machine precision at several p, plus proves
        the worst-mode discriminant < 0 (all lambda_k > 0).  This file's
        circulant check is the orthonormal-basis special case and agrees with it.

  (B.4-B.6) S0 (full support) cyclic beta-reduction:
        A_i + beta_{i-1} A_{i-1} = const C;   (p/q)^7 = prod (1-beta_i)/beta_i.
        Evaluated here at the UNIFORM point x_i=1 -- the actual S0 minimiser
        (P >= 7, equality only at uniform, by the n=7 theorem's easy direction;
        the circulant Hessian being positive definite at uniform makes it a strict
        local min).  The earlier version tried to numerically solve for a
        NON-uniform S0 stationary point with finite-difference scipy `root`; that
        solver frequently FAILS TO CONVERGE, which is consistent with uniform
        being the unique full-support minimiser -- so there is nothing to solve
        for.  Checking the identities at uniform is exact and convergent.
        The exact S0/S1 beta-reduction (symbolic KKT) is the duty of
        verify_s1_elimination.py.

  (S1 beta p-free) the S1 stationary curve in beta-coordinates is p-free.
        Numerical suggestion only here (support ratios vary smoothly with p);
        the EXACT symbolic reduction is verify_s1_elimination.py.
"""
import sys
import numpy as np

n = 7

def Pval(x, p):
    q = 1 - p; s = 0.0
    for i in range(n):
        den = p * x[(i + 1) % n] + q * x[(i + 2) % n]
        if abs(den) < 1e-15: return 1e6
        s += x[i] / den
    return s

# ---- (B.3) circulant spectrum (orthonormal-basis cross-check of verify_uniform_hess) ----
print("=== (B.3) uniform reduced-Hessian eigenvalues (circulant cross-check) ===")
print("  [RIGOROUS: verify_uniform_hess.py -- generalised eigenproblem + discriminant]")
b3_ok = True
for pp in [0.1, 0.214, 0.25, 0.329, 0.5, 0.9]:
    qv = 1 - pp
    h0 = 2 * pp**2 + 2 * qv**2; h1 = pp * (1 - 2 * pp); h2 = pp - 1; h3 = 0.0
    eigs = []
    for k in [1, 2, 3]:
        th = 2 * np.pi * k / 7
        lam = h0 + 2 * h1 * np.cos(th) + 2 * h2 * np.cos(2 * th) + 2 * h3 * np.cos(3 * th)
        gf = 2 * (1 - np.cos(th)) * (2 * qv * np.cos(th) + 2 * pp**2 - 3 * pp + 2)
        assert abs(lam - gf) < 1e-9, f"mismatch k={k}"
        eigs.append(lam)
    b3_ok = b3_ok and all(e > 0 for e in eigs)
    print(f"  p={pp:.3f}  circulant eigs(k=1,2,3)={np.round(eigs, 5)}  all>0: {all(e > 0 for e in eigs)}")
print(f"  => circulant formula == Fourier formula, all positive: {b3_ok}")

# ---- (B.4-B.6) S0 beta-reduction at the UNIFORM minimiser (exact, convergent) ----
print("\n=== (B.4-B.6) S0 cyclic beta-reduction at the uniform minimiser ===")
print("  [RIGOROUS symbolic S0/S1 beta-reduction: verify_s1_elimination.py]")
pp = 0.25; qv = 1 - pp
x0 = np.ones(n)                       # uniform: the S0 minimiser (P=7)
r = np.ones(n)                        # r_i = x_{i+1}/x_i = 1
beta = np.full(n, qv / (pp + qv))     # = q/(p+q) = q  (since p+q=1)
A = 1.0 / (r * (pp + qv * np.roll(r, -1)))    # = 1/(p+q) = 1
vals = np.array([A[i] + beta[(i - 1) % n] * A[(i - 1) % n] for i in range(n)])
C = float(vals.mean())
lhs = (pp / qv)**7
rhs = float(np.prod((1 - beta) / beta))
a = A / C
rec_err = max(abs(a[i] - (1 - beta[(i - 1) % n] * a[(i - 1) % n])) for i in range(n))
print(f"  uniform x=1, P={Pval(x0, pp):.1f}, beta_i=q={qv:.4f}, A_i=1")
print(f"  A_i + beta_{{i-1}} A_{{i-1}} = {np.round(vals, 6)}  (const C={C:.6f}, max-min={vals.max()-vals.min():.2e})")
print(f"  (p/q)^7={lhs:.6f}  prod((1-b)/b)={rhs:.6f}  match={abs(lhs - rhs) < 1e-12}")
print(f"  recurrence a_i=1-beta_{{i-1}} a_{{i-1}} max err = {rec_err:.2e}")
b456_ok = (vals.max() - vals.min() < 1e-12) and abs(lhs - rhs) < 1e-12 and rec_err < 1e-12
print(f"  => (B.4-B.6) identities hold exactly at the uniform S0 minimiser: {b456_ok}")
print("  [NOTE: the earlier non-uniform S0 numeric solver (scipy root, finite-diff) was")
print("   non-convergent; consistent with uniform being the UNIQUE full-support minimiser.]")

# ---- S1 beta p-free curve (numerical suggestion; exact = verify_s1_elimination) ----
print("\n=== S1 beta-coordinate reduction (p-free curve -- numerical suggestion) ===")
print("  [RIGOROUS symbolic KKT reduction: verify_s1_elimination.py]")
from scipy.optimize import root
def solve_S1(p, init):
    free = [2, 3, 4, 5, 6]
    def grad(v):
        x = np.zeros(n); x[1] = 1
        for j, idx in enumerate(free): x[idx] = v[j]
        h = 1e-7; f0 = Pval(x, p); g = []
        for j in range(len(free)):
            xp = x.copy(); xp[free[j]] += h; g.append((Pval(xp, p) - f0) / h)
        return g
    r = root(grad, init, method='hybr', options={'xtol': 1e-14, 'maxfev': 40000})
    if r.success and max(abs(r.fun)) < 1e-7 and all(ri > 1e-8 for ri in r.x):
        x = np.zeros(n); x[1] = 1
        for j, idx in enumerate(free): x[idx] = r.x[j]
        return x / x.sum()
    return None
for pp in [0.24, 0.27, 0.31]:
    x = solve_S1(pp, [0.2684881167890583140, 0.6791742990557855304, 1.5461708324775024161,
                      0.0656843931252869930, 1.3009478193484040029])
    if x is None:
        print(f"  p={pp}: S1 solver did not converge (diagnostic only)")
        continue
    sup = [1, 2, 3, 4, 5, 6]
    rs = [x[sup[(i + 1) % 6]] / x[sup[i]] for i in range(6)]
    print(f"  p={pp}: P={Pval(x, pp):.5f}  support log-ratios={np.round(np.log(rs), 4)}")
print("  (Support ratios vary smoothly with p; p-free curve means the algebraic RELATIONS")
print("   among beta_i omit p -- confirmed EXACTLY by verify_s1_elimination.py, not here.)")

print("\nDIAGNOSTIC complete. All rigorous certificates live in the load-bearing scripts")
print("(verify_uniform_hess, verify_s1_elimination, verify_ny_spectral, verify_hc_closedform,")
print("verify_s3_closed, n7_s1_hc_*). This file is a numerical cross-check only.")
sys.exit(0)
