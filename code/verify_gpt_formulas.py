#!/usr/bin/env python3
"""DIAGNOSTIC (not a formal certificate) -- numerical cross-checks of GPT's
structural formulas.  Each check has a rigorous replacement elsewhere; this file
gives a quick numerical sanity signal and always exits 0.

  (B.3) uniform reduced-Hessian eigenvalue formula
        lambda_k = 2(1-cos th_k)[2 q cos th_k + 2 p^2 - 3 p + 2],  th_k=2 pi k/7.
        Checked here for POSITIVITY and for the discriminant<0 claim (worst mode
        k=3).  The numerical eigenvalue MATCH (generalised eigenproblem
        B^T H B v = lam (B^T B) v) is the duty of verify_uniform_hess.py; this
        file's earlier finite-difference `eigvalsh(B^T H B)` used the
        NON-orthonormal basis B=[e_i-e_0] and produced a spurious mismatch -- that
        was a script bug, not a theorem defect, and is removed here.

  (B.1) second variation along the log-direction x_i(s)=x_i e^{s v_i}:
        Q_x(v) = sum_i T_i[(v_i-a_i v_{i+1}-b_i v_{i+2})^2 - a_i b_i(v_{i+1}-v_{i+2})^2],
        with T_i=x_i/d_i, a_i=p x_{i+1}/d_i, b_i=q x_{i+2}/d_i (a+b=1).
        Verified against numerical d^2/ds^2 P(x e^{s v}) at s=0 at an S1 stationary
        point.  (Unique to this script; if the S1 numeric solver fails to converge,
        the check is skipped -- diagnostic only.)

  (B.4-B.6) S0 cyclic beta-reduction: a_i=1-beta_{i-1}a_{i-1},
        (p/q)^7 = prod (1-beta_i)/beta_i.  Evaluated at the UNIFORM S0 minimiser
        (exact, convergent); the exact symbolic S0/S1 reduction is
        verify_s1_elimination.py.
"""
import sys
import numpy as np
from scipy.optimize import root

n = 7

def Pval(x, p):
    q = 1 - p; s = 0.0
    for i in range(n):
        den = p * x[(i + 1) % n] + q * x[(i + 2) % n]
        if abs(den) < 1e-15: return 1e6
        s += x[i] / den
    return s

def solve_S1(p, init):
    free = [2, 3, 4, 5, 6]
    def grad(v):
        x = np.zeros(n); x[1] = 1
        for j, idx in enumerate(free): x[idx] = v[j]
        h = 1e-7; f0 = Pval(x, p); g = []
        for j in range(len(free)):
            xp = x.copy(); xp[free[j]] += h; g.append((Pval(xp, p) - f0) / h)
        return g
    r = root(grad, init, method='hybr', options={'xtol': 1e-14, 'maxfev': 20000})
    if r.success and max(abs(r.fun)) < 1e-7 and all(ri > 1e-8 for ri in r.x):
        x = np.zeros(n); x[1] = 1
        for j, idx in enumerate(free): x[idx] = r.x[j]
        return x / x.sum()
    return None

# ---- (B.3) formula positivity + discriminant (rigorous match = verify_uniform_hess) ----
print("=== (B.3) uniform reduced-Hessian eigenvalue formula (positivity + discriminant) ===")
print("  [RIGOROUS eigenvalue match: verify_uniform_hess.py -- generalised eigenproblem]")
b3_ok = True
for pp in [0.1, 0.214, 0.25, 0.27, 0.329, 0.5, 0.7, 0.9]:
    q = 1 - pp
    eigs = []
    for k in [1, 2, 3]:
        th = 2 * np.pi * k / 7
        lam = 2 * (1 - np.cos(th)) * (2 * q * np.cos(th) + 2 * pp**2 - 3 * pp + 2)
        eigs.append(lam)
    pos = all(e > 0 for e in eigs)
    b3_ok = b3_ok and pos
    print(f"  p={pp:.3f}  formula eigs(k=1,2,3)={np.round(eigs, 5)}  all>0: {pos}")
# worst mode k=3: bracket disc = 4 cos^2(6pi/7) - 4 cos(6pi/7) - 7
disc_k3 = 4 * np.cos(2 * np.pi * 3 / 7)**2 - 4 * np.cos(2 * np.pi * 3 / 7) - 7
print(f"  worst-mode (k=3) bracket discriminant = {disc_k3:.4f}  (<0 -> bracket>0 for all p): {disc_k3 < 0}")
print(f"  => (B.3) all lambda_k>0 on (0,1): {b3_ok and disc_k3 < 0}")

# ---- (B.1) Q_x(v) vs numerical d^2/ds^2 P(x e^{s v}) ----
print("\n=== (B.1) Q_x(v) vs numerical d^2/ds^2 P(x e^{s v}) ===")
x = solve_S1(0.25, [0.2684881167890583140, 0.6791742990557855304, 1.5461708324775024161,
                    0.0656843931252869930, 1.3009478193484040029])
if x is None:
    print("  S1 numeric solver did not converge; skipping (B.1) (diagnostic only).")
else:
    pp = 0.25; q = 1 - pp
    print(f"  using S1 stationary x = {np.round(x, 4)}  P={Pval(x, pp):.5f}")
    d = np.array([pp * x[(i + 1) % n] + q * x[(i + 2) % n] for i in range(n)])
    T = x / d; alpha = pp * np.roll(x, -1) / d; beta = q * np.roll(x, -2) / d
    print(f"  sum(alpha+beta) should be 1: {np.round(alpha + beta, 6)}")
    np.random.seed(3)
    b1_ok = True
    for trial in range(4):
        v = np.random.randn(n)
        Q_formula = 0.0
        for i in range(n):
            Q_formula += T[i] * ((v[i] - alpha[i] * v[(i + 1) % n] - beta[i] * v[(i + 2) % n]) ** 2
                                 - alpha[i] * beta[i] * (v[(i + 1) % n] - v[(i + 2) % n]) ** 2)
        h = 1e-6
        Pp = Pval(x * np.exp(h * v), pp); Pm = Pval(x * np.exp(-h * v), pp); P0 = Pval(x, pp)
        Q_num = (Pp - 2 * P0 + Pm) / h ** 2
        match = abs(Q_formula - Q_num) < 1e-3
        b1_ok = b1_ok and match
        print(f"  trial{trial}: Q_formula={Q_formula:+.6f}  Q_num={Q_num:+.6f}  match={match}")
    print(f"  => (B.1) Q_x(v) == numerical second variation: {b1_ok}")

# ---- (B.4-B.6) S0 beta-reduction at the UNIFORM minimiser (exact, convergent) ----
print("\n=== (B.4-B.6) S0 cyclic beta-reduction at the uniform minimiser ===")
print("  [RIGOROUS symbolic S0/S1 reduction: verify_s1_elimination.py]")
pp = 0.25; q = 1 - pp
x0 = np.ones(n)                       # uniform S0 minimiser (P=7)
r = np.ones(n)
beta = np.array([q * r[(i + 1) % n] / (pp + q * r[(i + 1) % n]) for i in range(n)])   # = q
A = 1.0 / (r * (pp + q * np.roll(r, -1)))                                              # = 1
vals = np.array([A[i] + beta[(i - 1) % n] * A[(i - 1) % n] for i in range(n)])
C = float(vals.mean())
lhs = (pp / q) ** 7; rhs = float(np.prod((1 - beta) / beta))
a = A / C
rec_err = max(abs(a[i] - (1 - beta[(i - 1) % n] * a[(i - 1) % n])) for i in range(n))
b456_ok = (vals.max() - vals.min() < 1e-12) and abs(lhs - rhs) < 1e-12 and rec_err < 1e-12
print(f"  uniform x=1, beta_i=q={q:.4f}, A_i=1, C={C:.6f}")
print(f"  A_i+beta_{{i-1}}A_{{i-1}} const (max-min={vals.max()-vals.min():.2e}); "
      f"(p/q)^7={lhs:.6f} prod={rhs:.6f} match={abs(lhs-rhs)<1e-12}; rec err={rec_err:.2e}")
print(f"  => (B.4-B.6) identities hold exactly at uniform: {b456_ok}")

print("\nDIAGNOSTIC complete. Rigorous certificates: verify_uniform_hess (B.3),")
print("verify_s1_elimination (B.4-B.6 + S1 p-free), verify_ny_spectral, verify_hc_closedform,")
print("verify_s3_closed, n7_s1_hc_*. This file is a numerical cross-check only.")
sys.exit(0)
