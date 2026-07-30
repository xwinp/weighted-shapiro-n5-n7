#!/usr/bin/env python3
"""
Route B leg-1 local prep: test the "average curvature = const*(P-7)" identity
for the n=7 weighted Shapiro cyclic functional, INDEPENDENTLY of GPT.

P_{n,p,q}(x) = sum_{i=1}^n x_i / (p*x_{i+1} + q*x_{i+2}),  indices mod n, q=1-p.
Normalize p+q=1, target bound n=7.

Question (Nowosad / averaged quadratic form):
  Does there exist a probability vector of directions {v^(m)} (e.g. DFT/cosine basis)
  and weights c_m >= 0 such that
      Sum_m c_m * v^(m)^T H_x v^(m)  =  C * (P(x) - 7)            (*)
  for ALL x>0, with C>0 constant?  If yes, then P<7 => some direction has negative
  curvature => H_x not PSD => x not a local min -- the Route B certificate, no
  stationary-point enumeration needed.

Strategy:
  (a) Build the FULL 7x7 Hessian of P at x (symbolic-free, via autodiff / finite
      differences with the CORRECT three-point diagonal formula -- avoid the fd-bug).
  (b) Compute the DFT basis v^(m)_k = cos(2*pi*m*k/7) and sin (real basis of R^7,
      excluding the constant mode v^(0)=1 which is the scaling direction, degenerate
      because P is degree-0 -> curvature 0 along 1).
  (c) For random x>0 and a few p in the band, compute:
        A(x) := Sum_{m=1..6} v^(m)^T H_x v^(m) / (v^(m)^T v^(m))
      and compare A(x) / (P(x)-7) across x. If the ratio is a p-only (or universal)
      constant, identity (*) holds with the uniform weight c_m = 1/6.
  (d) Also test the "Laplacian-trace" version: trace(H_x) restricted to the
      sum-zero subspace = Sum of all eigenvalues except the scaling one. For a
      circulant-ish Hessian this equals the DFT average. Check trace_relation.

If A(x)/(P-7) is constant -> clean certificate found locally; report C(p).
If not constant -> uniform DFT weights don't work; we'll need GPT's tuned c_m.
"""
import numpy as np

n = 7

def P_val(x, p):
    q = 1.0 - p
    s = 0.0
    for i in range(n):
        s += x[i] / (p * x[(i+1) % n] + q * x[(i+2) % n])
    return s

def grad_P(x, p):
    q = 1.0 - p
    g = np.zeros(n)
    # dP/dx_j: term i contributes when j is the numerator (i=j): 1/(p x_{j+1}+q x_{j+2})
    #          and when j appears in a denominator: -x_i * (p*delta_{j,i+1}+q*delta_{j,i+2}) / (den)^2
    den = np.array([p * x[(i+1) % n] + q * x[(i+2) % n] for i in range(n)])
    for j in range(n):
        # numerator derivative
        g[j] += 1.0 / den[j]
    for i in range(n):
        # denominator derivatives: d(den_i)/dx_{i+1} = p, d/dx_{i+2} = q
        g[(i+1) % n] += -x[i] * p / den[i]**2
        g[(i+2) % n] += -x[i] * q / den[i]**2
    return g

def hessian_P(x, p, h=1e-5):
    """Full 7x7 Hessian. Off-diagonal: 4-point mixed diff. Diagonal: 3-point.
    (Avoids the fd-bug: never use the mixed formula on i==j.)"""
    q = 1.0 - p
    H = np.zeros((n, n))
    # diagonal: 3-point  [f(+h)-2f0+f(-h)]/h^2
    for i in range(n):
        xp = x.copy(); xp[i] += h
        xm = x.copy(); xm[i] -= h
        H[i, i] = (P_val(xp, p) - 2 * P_val(x, p) + P_val(xm, p)) / (h * h)
    # off-diagonal: 4-point mixed  [f(++)-f(+-)-f(-+)+f(--)]/(4h^2)
    for i in range(n):
        for j in range(i+1, n):
            xpp = x.copy(); xpp[i] += h; xpp[j] += h
            xpm = x.copy(); xpm[i] += h; xpm[j] -= h
            xmp = x.copy(); xmp[i] -= h; xmp[j] += h
            xmm = x.copy(); xmm[i] -= h; xmm[j] -= h
            val = (P_val(xpp, p) - P_val(xpm, p) - P_val(xmp, p) + P_val(xmm, p)) / (4 * h * h)
            H[i, j] = val
            H[j, i] = val
    return H

def dft_basis():
    """Real orthonormal basis of the sum-zero subspace of R^7 (6 vectors),
    plus the constant mode (scaling). v_m,k = sqrt(2/7) cos(2 pi m k/7), m=1,2,3;
    sin for m=1,2,3; (m=4,5,6 are redundant with cos/sin by symmetry for n=7 prime)."""
    V = []
    for m in range(1, 4):
        c = np.sqrt(2.0 / n) * np.cos(2 * np.pi * m * np.arange(n) / n)
        s = np.sqrt(2.0 / n) * np.sin(2 * np.pi * m * np.arange(n) / n)
        V.append(c)
        V.append(s)
    # constant mode (scaling direction, P degree-0 -> should have ~0 curvature)
    ones = np.ones(n) / np.sqrt(n)
    return np.array(V), ones  # V: 6 x 7

def main():
    V, ones = dft_basis()
    print("=== Test identity  Sum_{m=1..6} v^T H v / (v^T v)  ?=  C * (P-7) ===")
    print(f"(using orthonormal DFT basis, v^T v = 1)\n")
    rng = np.random.default_rng(0)
    for p in [0.22, 0.25, 0.27, 0.30, 0.33]:
        print(f"--- p = {p} ---")
        ratios = []
        for trial in range(6):
            if trial == 0:
                x = np.ones(n)  # uniform point: P=7 exactly, curvature baseline
            else:
                x = np.exp(rng.normal(0, 0.5, n))
            P = P_val(x, p)
            H = hessian_P(x, p)
            # average curvature over the 6 DFT (sum-zero) directions
            A = sum(float(v @ H @ v) for v in V)  # weights uniform 1/6 -> factor 1/6
            A_avg = A / 6.0
            # also curvature along constant mode (should be ~0)
            curv_const = float(ones @ H @ ones)
            if abs(P - 7) > 1e-9:
                ratio = A_avg / (P - 7)
            else:
                ratio = float('nan')
            ratios.append(ratio)
            print(f"  x{trial}: P-7={P-7:+.5f}  A_avg={A_avg:+.6f}  "
                  f"ratio={ratio:+.6f}  curv_const={curv_const:+.2e}")
        # consistency: do all ratios agree (-> constant C)?
        valid = [r for r in ratios if not np.isnan(r)]
        if valid:
            print(f"  -> C range [{min(valid):+.6f}, {max(valid):+.6f}]  "
                  f"spread={max(valid)-min(valid):.2e}")
        print()

    print("=== Also: trace(H) on sum-zero subspace (= sum of 6 non-const eigenvalues) ===")
    p = 0.27
    rng = np.random.default_rng(1)
    for trial in range(5):
        x = np.exp(rng.normal(0, 0.5, n))
        P = P_val(x, p)
        H = hessian_P(x, p)
        w = np.linalg.eigvalsh(H)
        # smallest eigenvalue is the scaling direction (degree-0 -> ~0); sort
        w_sorted = np.sort(w)
        trace_sumzero = sum(w_sorted[1:])  # drop the ~0 scaling eigenvalue
        print(f"  P-7={P-7:+.5f}  trace_sumzero={trace_sumzero:+.6f}  "
              f"ratio={trace_sumzero/(P-7) if abs(P-7)>1e-9 else float('nan'):+.6f}  "
              f"eig0={w_sorted[0]:+.2e}")

if __name__ == "__main__":
    main()
