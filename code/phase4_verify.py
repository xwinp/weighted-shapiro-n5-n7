#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 4 (general odd n) verification of GPT's claims. LOCAL CHECKS ONLY.
GPT inversion-error pattern is known; EVERY explicit formula re-derived here.

Claims to verify:
 (1) n=5: m_5(p)=5 for all p (H_5=(0,1)). Tuan-Thuong interval only sufficient.
     - S2 face (2 zeros): P=a/(qb)+b/(pc)+c/(qa) >= 3/(p^{1/3} q^{2/3}), min over p at p=1/3 = 5.6696.
     - S1 face (1 zero): positive coeffs Q_11(r) => branch never crosses 5; M_{5,5}(r=1)=5.574.
     - S0: Nowosad-Yamagami (uniform Hessian PSD for n=5).
 (2) L=3 closed form: M_{n,3} = (m+1)/(p q^m)^{1/(m+1)}, n=2m+1.
     Failure condition (re-derived): p q^m > ((m+1)/(2m+1))^{m+1}.  [GPT (3.3) INVERTED this]
     Check n=7: never fails (M_{7,3}=4/(p^{1/4}q^{3/4})=16/3^{3/4}>7). n=9: fails.
 (3) Uniform Hessian spectral instability: mu_k=2(1-cos th_k)[2q cos th_k+2p^2-3p+2].
     Dangerous mode k=(n-1)/2, c_n=-cos(pi/n). Delta_n=4c_n^2-4c_n-7.
     n=5,7: Delta<0 (stable all p). n>=9: Delta>0, uniform saddle in (p_n-,p_n+).
     p_n+/- = [(3+2 c_n) +/- sqrt(Delta_n)]/4.   [GPT (7.1) rendering garbled]
 (4) Left-endpoint table (L=3 roots) a_n for n=7..15; asymptotic coeff: sqrt(e) NOT e.
 (5) p_infty = 0.5250803166... root of (1-p)(2-p)=exp(-p/(2-p)).  b_n -> p_infty > 1/2.
 (6) Numerical m_9(p): uniform saddle + L=3 face => m_9<9 for p in failure band; b_9~0.434.
"""
import numpy as np
from scipy.optimize import minimize, brentq
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

def P_val(x, n, p):
    q = 1 - p
    s = 0.0
    for i in range(n):
        d = p*x[(i+1) % n] + q*x[(i+2) % n]
        s += x[i] / d
    return s

def m_n_numerical(n, p, nstarts=200, seed=0):
    """Numerical inf of P over simplex x_i>=0 (interior + boundary via many starts)."""
    rng = np.random.default_rng(seed)
    best = np.inf
    for _ in range(nstarts):
        # random support: drop a random independent-set of zeros is hard; just random positive + boundary
        x = rng.dirichlet(np.ones(n))
        # also try sparse: set some entries ~0
        if _ % 3 == 0:
            k = rng.integers(0, n//2)
            idx = rng.choice(n, k, replace=False)
            x[idx] = 1e-9
            x = x / x.sum()
        try:
            res = minimize(lambda xx: P_val(np.abs(xx), n, p), x, method='Nelder-Mead',
                           options={'maxiter': 4000, 'xatol': 1e-10, 'fatol': 1e-11})
            v = P_val(np.abs(res.x), n, p)
            if v < best:
                best = v
        except Exception:
            pass
    return best

# ---------- (1) n=5 ----------
print("="*70); print("(1) n=5: m_5(p)=5 ?  (H_5=(0,1))"); print("="*70)
# S2 face AM-GM: 3/(p^{1/3} q^{2/3}), min over p
ps = np.linspace(0.001, 0.999, 2000)
s2_n5 = 3.0/(ps**(1/3)*(1-ps)**(2/3))
print(f" S2(n=5) AM-GM min = {s2_n5.min():.6f} at p={ps[s2_n5.argmin()]:.4f} (expect 5.6696 at 1/3)")
# S1 face Q_11(r) all-positive?
Q11 = [1024,3644,6360,11580,15895,28772,23892,21120,9360,7680,2368,1728]
print(f" S1(n=5) Q_11 all coeffs positive: {all(c>0 for c in Q11)}  (deg {len(Q11)-1})")
# M_{5,5}(r=1): need y solving (y^2-1)^3 = r^5 y^2 with r=1
from numpy.polynomial import polynomial as Pp
def y_solve_n5(r):
    f = lambda y: (y**2-1)**3 - r**5 * y**2
    return brentq(f, 1.0, 10.0)
def M55(r):
    y = y_solve_n5(r)
    return 2*r**2/(1+r) * (y-1) * (8 + 1.0/y**4)
print(f" M_{{5,5}}(r=1) = {M55(1.0):.6f}  (expect 5.574158)")
# full numerical m_5 across p
for p in [0.1, 0.2, 0.25, 0.3, 0.333, 0.4, 0.5, 0.6, 0.7, 0.8]:
    m = m_n_numerical(5, p, nstarts=300, seed=1)
    print(f" m_5(p={p:.3f}) ~ {m:.5f}  ({'OK>=5' if m>=5-1e-3 else 'FAIL<5'})")

# ---------- (2) L=3 closed form & failure condition ----------
print("\n"+"="*70); print("(2) L=3 closed form M_{n,3}=(m+1)/(p q^m)^{1/(m+1)}, failure cond"); print("="*70)
for n in [5,7,9,11,13,15]:
    m = (n-1)//2
    # M_{n,3} min over p: minimize (m+1)/(p q^m)^{1/(m+1)}  ==  maximize p q^m
    # max p q^m at p = 1/(m+1)
    pstar = 1.0/(m+1)
    maxpq = pstar*(1-pstar)**m
    Mmin = (m+1)/maxpq**(1.0/(m+1))
    # failure cond (re-derived): p q^m > ((m+1)/(2m+1))^{m+1}
    thr = ((m+1)/(2*m+1))**(m+1)
    fails = maxpq > thr
    print(f" n={n} m={m}: M_{{n,3}}^min={Mmin:.5f} {'<' if Mmin<n else '>='}{n}  "
          f"maxpq={maxpq:.6f} thr={thr:.6f} -> {'FAILS' if fails else 'no-fail'}  "
          f"[GPT (3.3) thr_inv={1/thr:.2e} WRONG]")
# specifically n=7 M_{7,3} = 16/3^{3/4}
print(f" M_{{7,3}}^min = {4/( (1/4)**(1/4)*(3/4)**(3/4) ):.6f} = 16/3^(3/4) = {16/3**0.75:.6f}")

# ---------- (3) uniform Hessian instability ----------
print("\n"+"="*70); print("(3) Uniform Hessian: Delta_n, p_n+/-  (formula [(3+2c)+/-sqrt(Delta)]/4)"); print("="*70)
for n in [5,7,9,11,13,15,25]:
    cn = -np.cos(np.pi/n)
    Delta = 4*cn**2 - 4*cn - 7
    if Delta > 0:
        pm = ((3+2*cn) - np.sqrt(Delta))/4
        pp = ((3+2*cn) + np.sqrt(Delta))/4
        print(f" n={n}: c_n={cn:.6f} Delta={Delta:+.4f}>0  p_n-={pm:.6f} p_n+={pp:.6f}  uniform SADDLE")
    else:
        print(f" n={n}: c_n={cn:.6f} Delta={Delta:+.4f}<0  uniform stable all p")
# verify uniform Hessian actually has negative eigenvalue at n=9, p=0.3
print("  check n=9,p=0.3 uniform Hessian eigenvalues (1^perp):")
n=9; p=0.3; q=1-p
mus=[]
for k in range(1,n):
    th=2*np.pi*k/n
    mus.append(2*(1-np.cos(th))*(2*q*np.cos(th)+2*p**2-3*p+2))
print(f"    min mu_k = {min(mus):.4f}  ({'NEGATIVE -> saddle' if min(mus)<0 else 'positive'})")

# ---------- (4) left-endpoint table (L=3 small root) ----------
print("\n"+"="*70); print("(4) Left endpoints a_n = small root of p q^m = ((m+1)/(2m+1))^{m+1}"); print("="*70)
gpt_table = {7:0.2142735209, 9:0.0708264, 11:0.0307780, 13:0.0143085, 15:0.00686968}
print(f" {'n':>3} {'a_n(computed)':>16} {'a_n(GPT)':>14} {'sqrt(e)*2^-(n+1)/2':>22} {'e*2^-(n+1)/2':>16}")
for n in [7,9,11,13,15]:
    m=(n-1)//2
    thr = ((m+1)/(2*m+1))**(m+1)
    # small root of p(1-p)^m = thr
    f = lambda p: p*(1-p)**m - thr
    # for n=7 L=3 never fails (f<0 everywhere), so no root -- a_7 is L=5 not L=3
    if n==7:
        print(f" {n:>3} {'(L=5, not L=3)':>16} {gpt_table[n]:>14} {np.exp(0.5)*2**(-(n+1)/2):>22.6f} {np.exp(1)*2**(-(n+1)/2):>16.6f}")
        continue
    a = brentq(f, 1e-12, 1.0)
    sq = np.exp(0.5)*2**(-(n+1)/2)
    e1 = np.exp(1)*2**(-(n+1)/2)
    print(f" {n:>3} {a:>16.7f} {gpt_table[n]:>14.7f} {sq:>22.6f} {e1:>16.6f}")
# ratio test for coefficient
print("  ratio a_n / 2^-(n+1)/2 (should -> sqrt(e)=1.6487, NOT e=2.718):")
for n in [9,11,13,15]:
    m=(n-1)//2; thr=((m+1)/(2*m+1))**(m+1)
    a = brentq(lambda p: p*(1-p)**m - thr, 1e-12, 1.0)
    print(f"    n={n}: a_n/2^-(n+1)/2 = {a/2**(-(n+1)/2):.4f}")

# ---------- (5) p_infty ----------
print("\n"+"="*70); print("(5) p_infty root of (1-p)(2-p)=exp(-p/(2-p))"); print("="*70)
f = lambda p: (1-p)*(2-p) - np.exp(-p/(2-p))
pinf = brentq(f, 0.5, 1.0)
print(f" p_infty = {pinf:.16f}  (GPT 0.5250803166496057)")
print(f" check: (1-p)(2-p)={(1-pinf)*(2-pinf):.6f}  exp(-p/(2-p))={np.exp(-pinf/(2-pinf)):.6f}")
print(f" p_infty > 1/2 ? {pinf>0.5}  -> b_n -> p_infty > 1/2 (my 1/2 guess WRONG)")

# ---------- (6) numerical m_9 ----------
print("\n"+"="*70); print("(6) Numerical m_9(p) — should be <9 in failure band"); print("="*70)
for p in [0.10, 0.1453, 0.20, 0.30, 0.4150, 0.45, 0.5]:
    m = m_n_numerical(9, p, nstarts=400, seed=2)
    print(f" m_9(p={p:.4f}) ~ {m:.5f}  ({'<9 FAIL' if m<9-1e-3 else '>=9'})")
