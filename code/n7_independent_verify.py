#!/usr/bin/env python3
"""
INDEPENDENT verification of the H_7 classification (fresh code, no reuse).
Claims to verify:
  (C1) a7,b7 are the 2 roots in (0,1) of the degree-15 F (irreducible).
  (C2) Global SLSQP minimum m_7(p) < 7  <=>  p in (a7,b7).
  (C3) Inside band, global minimizer support == {1,2,3,4,6} (2z-d2, zeros {0,5}).
       Outside band, minimizer == uniform (P=7).
  (C4) P_S4(p)=4/(p^{1/4}q^{3/4}) >= 16/3^{3/4} > 7  (so S4 never the global min).
  (C5) On the S2 stationary curve, P<7 exactly on (a7,b7); endpoints P=7.
"""
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
mp.mp.dps = 50

n = 7
def Pval(x, p):
    q = 1 - p
    s = 0.0
    for i in range(n):
        den = p*x[(i+1)%n] + q*x[(i+2)%n]
        if den <= 1e-15: return 1e6
        s += x[i]/den
    return s

def global_min(p, nstarts=40):
    best = 1e9; bx = None
    rng = np.random.RandomState(98765)
    # include uniform start + many random
    starts = [np.full(n, 1.0/n)]
    for _ in range(nstarts):
        y = rng.rand(n)+1e-3; y/=y.sum(); starts.append(y)
    for x0 in starts:
        cons = [{'type':'eq','fun':lambda x:x.sum()-1.0}]
        bounds = [(0,1)]*n
        r = minimize(lambda x:Pval(x,p), x0, method='SLSQP', bounds=bounds,
                     constraints=cons, options={'maxiter':600,'ftol':1e-14})
        if r.fun < best: best = r.fun; bx = r.x
    return best, bx

# C1: F, its real roots in (0,1)
Fcoeffs = [5764801,-47765494,190003135,-486209703,901678743,-1287828143,1464952167,
 -1351039522,1017028633,-624621984,310300032,-122238368,36836352,-7952896,1073408,-65536]
Froots = [mp.re(r) for r in mp.polyroots(Fcoeffs, extraprec=80) if abs(mp.im(r))<1e-30 and 0<mp.re(r)<1]
Froots.sort()
a7, b7 = Froots[0], Froots[1]
print(f"(C1) F real roots in (0,1): {[mp.nstr(r,15) for r in Froots]}")
print(f"     a7={mp.nstr(a7,20)}  b7={mp.nstr(b7,20)}")

# C4: S4 closed form
def P_S4(p): return 4/(p**0.25*(1-p)**0.75)
print(f"\n(C4) P_S4 min = 16/3^(3/4) = {mp.nstr(mp.mpf(16)/mp.power(3,mp.mpf(3)/4),15)}  (>7: {16**4 > 7**4*3**3})")

# C2 + C3: global SLSQP scan
print("\n(C2/C3) global SLSQP scan vs predicted band (a7,b7):")
print(f"  predicted band = ({mp.nstr(a7,8)}, {mp.nstr(b7,8)})")
print(f"  {'p':>7s} {'m_7(p) SLSQP':>13s} {'<7?':>5s} {'in band?':>9s} {'match':>6s} {'min support zeros':>20s}")
allmatch = True
for pp in np.linspace(0.05, 0.49, 19):
    mp_, bx = global_min(pp)
    in_band = (a7 < mp.mpf(pp) < b7)
    pred_fail = in_band
    actual_fail = mp_ < 7 - 1e-6
    match = (pred_fail == actual_fail)
    allmatch &= match
    zeros = sorted([i for i in range(n) if bx[i] < 1e-4])
    # canonicalize zeros to a dihedral rep to compare with {0,5}/{0,2}
    print(f"  {pp:7.3f} {mp_:13.6f} {'YES' if actual_fail else 'no':>5s} {'yes' if in_band else 'no':>9s} {'OK' if match else 'XX':>6s} zeros={zeros}")
print(f"\n  ALL p match predicted band: {allmatch}")

# C5: on S2 curve, P at endpoints == 7, <7 inside
# curve R=q^3-p^3 t^5 - p^2 q t^8=0; for given p solve t, compute P via reduction
def curve_P(pp):
    pv=mp.mpf(pp); qv=1-pv
    f=lambda t: qv**3 - pv**3*t**5 - pv**2*qv*t**8
    tv=mp.findroot(f, mp.mpf('1.3'))
    b = qv/(pv*tv) - qv**2*(qv-pv*tv**4)/(pv**3*tv**2)
    c = qv*(qv-pv*tv**4)/(pv**2*tv**2)
    d = tv**2
    P = 1/(pv*b+qv*c) + b/(pv*c+qv*d) + c/(pv*d) + d/(qv*tv) + tv/qv
    return P
print("\n(C5) S2 curve value P at endpoints (should be 7) and midband (should be <7):")
for pp in [float(a7)+1e-4, 0.27, float(b7)-1e-4]:
    print(f"  p={pp:.5f}  P_curve={float(curve_P(pp)):.8f}")
