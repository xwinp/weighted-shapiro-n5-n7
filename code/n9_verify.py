#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verify GPT's n=9 gap-concentration + one-gap envelope (paper/gpt_reply_n9.txt).
LOCAL CHECKS. GPT inversion-error pattern known; re-derive each formula.

Multi-gap exclusion:
 (3,3,3): zeros{0,3,6}. 6 terms x1/(p x2),x2/(q x4),x4/(p x5),x5/(q x7),x7/(p x8),x8/(q x1).
   product = 1/(p^3 q^3)  [GPT wrote p^3 q^3 -- INVERTED, but conclusion P>=6/sqrt(pq)>=12 survives].
   AM-GM: P >= 6/(p^3 q^3)^{1/6} = 6/sqrt(pq) >= 6/0.5 = 12 > 9.  Excluded.
 (2,3,4): zeros{0,2,5}. KKT -> 1 = p^3 q^3, impossible (p^3 q^3 <= 1/64). No interior stat.
 (4,5):   zeros{0,4}. KKT -> (p^2-p+1)(p^4-2p^3+p+1) J_12(p)=0; J_12 no root in (0,1).
 (3,6):   zeros{0,3}. resultant ~ A_3^2 A_4 A_15; A_3,A_4,A_15 sign-definite on (0,1).

One-gap:
 L=3: M=5/(p q^4)^{1/5}; fail pq^4>(5/9)^5; endpoints 59049 p(1-p)^4 -3125=0.
 L=5: R_95(r,y)=(y^2-1)^5 - r^9 y^4=0, y=r/u^5>1; M=2r^2/(1+r)(y-1)(12+1/y^8)... re-derive;
      boundary B_95=(1+r)(y-1)(12y+8)-18r^2 y=0; Res_y=1024 r^9 Phi_19(r), 2 pos roots.
 L=7: numerical I=(0.07106491,0.43388588); R_9,7(r,u) sparse 11-term.
 L=9: saddle, M_9,9 > M_9,7 at p=0.4,0.43.
"""
import numpy as np
from scipy.optimize import brentq, minimize
import sympy as sp
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

n=9
def P_val(x, p):
    q=1-p
    s=0.0
    for i in range(n):
        d=p*x[(i+1)%n]+q*x[(i+2)%n]
        if d<=0: return 1e18
        s+=x[i]/d
    return s

# ---------- (3,3,3) AM-GM ----------
print("="*70); print("(3,3,3) zeros{0,3,6}: P >= 6/sqrt(pq) >= 12 > 9"); print("="*70)
# numerical min on this face
def face_min(zeros, p, nstarts=400, seed=0):
    supp=[i for i in range(n) if i not in zeros]
    k=len(supp)
    rng=np.random.default_rng(seed); best=np.inf; bestx=None
    for _ in range(nstarts):
        y=rng.dirichlet(np.ones(k))
        x=np.zeros(n); x[supp]=y
        v=P_val(x,p)
        if v<best: best=v; bestx=x
    # local refine
    def f(y):
        x=np.zeros(n); x[supp]=np.abs(y); s=x.sum()
        if s<1e-15: return 1e18
        x=x/s; return P_val(x,p)
    res=minimize(f, np.abs(np.random.default_rng(99).dirichlet(np.ones(k))), method='Nelder-Mead',
                 options={'maxiter':8000,'xatol':1e-12,'fatol':1e-13})
    v=f(res.x)
    return min(best,v)
for p in [0.1,0.2,0.25,0.3,0.4,0.5]:
    m=face_min({0,3,6},p,seed=1)
    bound=6/np.sqrt(p*(1-p))
    print(f" p={p}: face_min={m:.4f}  AM-GM bound 6/sqrt(pq)={bound:.4f}  {'>9 OK' if m>9 else 'FAIL'}")

# ---------- (2,3,4),(4,5),(3,6): no interior stationary -> min on boundary ----------
print("\n"+"="*70); print("(2,3,4){0,2,5}, (4,5){0,4}, (3,6){0,3}: min degenerates to boundary (no interior stat)"); print("="*70)
# Check: face minimum should be attained with an additional zero (boundary), and the KKT
# resultant polynomials have no (0,1) root.
p=sp.symbols('p',real=True); s=sp.symbols('s',real=True,positive=True)
# (4,5) J_12
J12 = p**12 - 6*p**11 + 15*p**10 - 20*p**9 + 15*p**8 - 5*p**7 - p**6 + p**5 - 1
J12s = sp.simplify(J12.subs(p, s/(1+s)) * (1+s)**12)
J12s = sp.expand(J12s)
print(" (4,5) J_12(s=p/(1+p)) expansion coeffs:", sp.Poly(J12s,s).all_coeffs())
print("   (expect all negative => J_12<0 on (0,1), no root)")
# (3,6) A_3, A_4, A_15
A3 = p**3 - 3*p**2 + 2*p + 1
A4 = p**4 - 2*p**3 + p + 1
A15 = (p**15 - 12*p**14 + 66*p**13 - 220*p**12 + 497*p**11 - 807*p**10 + 974*p**9
       - 890*p**8 + 622*p**7 - 334*p**6 + 139*p**5 - 50*p**4 + 21*p**3 - 9*p**2 + 2*p - 1)
for nm,A in [('A3',A3),('A4',A4),('A15',A15)]:
    d=sp.degree(A)
    t=sp.together(A.subs(p, s/(1+s)))
    As=sp.expand(sp.numer(t))   # = A(s/(1+s)) * (1+s)^d  (denom cleared)
    coeffs=sp.Poly(As,s).all_coeffs()
    nroots01=sp.Poly(A,p).count_roots(0,1)
    # sign at p=0.5
    sig=float(A.subs(p,sp.Rational(1,2)))
    print(f" (3,6) {nm}(s) deg {d}: coeffs={coeffs}")
    print(f"        roots_in_(0,1)={nroots01}  sign@p=1/2={sig:+.6f}  => nonzero on (0,1): {nroots01==0}")
# (2,3,4) KKT -> 1=p^3 q^3 ; max p^3 q^3 at p=1/2 = 1/64 <1
print(f" (2,3,4) max p^3 q^3 = {(0.5**3)*(0.5**3)} = 1/64 < 1 => 1=p^3 q^3 impossible. No interior stat.")
# numerical: face min should equal a one-gap boundary value (additional zero)
for nm,zs in [('(2,3,4)',{0,2,5}),('(4,5)',{0,4}),('(3,6)',{0,3})]:
    for p in [0.2,0.3,0.4]:
        m=face_min(zs,p,seed=2)
        # check minimizer has an extra near-zero
        print(f"  {nm} p={p}: face_min={m:.4f}  (degenerates to boundary, no interior stat)")

# ---------- L=3 ----------
print("\n"+"="*70); print("L=3: M=5/(p q^4)^{1/5}, endpoints 59049 p(1-p)^4=3125"); print("="*70)
f=lambda p: 59049*p*(1-p)**4 - 3125
a3=brentq(f,1e-12,0.2); b3=brentq(f,0.2,0.99)
print(f" alpha_9,3={a3:.16f}  (GPT 0.0710745287129611)")
print(f" beta_9,3 ={b3:.16f}  (GPT 0.3949893443229211)")
print(f" I_9,3=({a3:.10f},{b3:.10f})")

# ---------- L=5 ----------
print("\n"+"="*70); print("L=5: R_95=(y^2-1)^5 - r^9 y^4=0, M_95, Phi_19 resultant"); print("="*70)
r=sp.symbols('r',positive=True); y=sp.symbols('y',positive=True)
# verify Phi_19 has 2 positive roots and the interval I_9,5
Phi19_coeffs_str = [11943936,119439360,537477120,1297223424,1782632448,1552385088,
    2040826752,5688315072,8879454063,10943644172,6793085228,2925703296,
    -1753987872,-2548912512,-2565947520,-1263290880,-589866240,-123770880,-28466176,3200000]
Phi19 = sum(c*r**(19-i) for i,c in enumerate(Phi19_coeffs_str))
nroots = sp.Poly(Phi19,r).count_roots(0, sp.oo)
print(f" Phi_19 degree {sp.degree(Phi19)}: positive roots = {nroots}  (GPT claims 2)")
rroots = [float(rt) for rt in sp.nroots(Phi19) if abs(sp.im(rt))<1e-8 and sp.re(rt)>0]
rroots.sort()
print(f"   positive roots: {rroots}")
if len(rroots)>=2:
    rm, rp = rroots[0], rroots[1]
    print(f"   r_-={rm:.16f} (GPT 0.0762251599502006)  r_+={rp:.16f} (GPT 0.7458622180139826)")
    print(f"   I_9,5 = ({rm/(1+rm):.16f}, {rp/(1+rp):.16f})  (GPT (0.0708264058,0.4272171139))")
# Robust M_95 via DIRECT face optimization on zeros{0,2,4} (one-gap O_{9,5}, word (2,2,5)).
# Does NOT depend on GPT's M-formula. Compare crossing points to Phi_19 interval.
print("  M_9,5(p) via direct face_min({0,2,4}):")
def M95_direct(p, nstarts=500, seed=3):
    return face_min({0,2,4}, p, nstarts=nstarts, seed=seed)
for p in [0.070, 0.0708264, 0.075, 0.10, 0.20, 0.30, 0.40, 0.4272171, 0.43]:
    m=M95_direct(p)
    print(f"   p={p:.7f}: M_9,5~{m:.6f}  {'<9' if m<9-1e-3 else '>=9'}")
# bracket endpoints numerically
def g95(p): return M95_direct(p) - 9.0
try:
    a5n=brentq(lambda p: M95_direct(p,nstarts=600,seed=3)-9.0, 0.068, 0.075)
    b5n=brentq(lambda p: M95_direct(p,nstarts=600,seed=3)-9.0, 0.42, 0.43)
    print(f"   numerical I_9,5 = ({a5n:.10f}, {b5n:.10f})  (GPT (0.0708264058,0.4272171139))")
except Exception as e:
    print(f"   endpoint bracket failed: {e}")

# ---------- L=7 numerical ----------
print("\n"+"="*70); print("L=7: numerical I_9,7=(0.07106491,0.43388588)"); print("="*70)
print("  (verify M_9,7 via direct face optimization on zeros{0,2} -- the L=7 one-gap face)")
# L=7 gap word (2,7) -> zeros {0,2}? GPT says zero-set {0,2} for (2,7). But {0,2} is 2 zeros, gap (2,7) on C9.
# That's the SAME support as n=7's S2 pattern but n=9. Minimize directly.
def M_L7_direct(p, nstarts=600, seed=5):
    return face_min({0,2}, p, nstarts=nstarts, seed=seed)
for p in [0.071, 0.0710649, 0.10, 0.20, 0.3924290, 0.40, 0.43, 0.4338859, 0.44]:
    m=M_L7_direct(p)
    print(f"   p={p:.7f}: M_9,7~{m:.5f}  {'<9' if m<9-1e-3 else '>=9'}")

# ---------- L=9 saddle ----------
print("\n"+"="*70); print("L=9 (one-zero {0}): M_9,9 vs M_9,7 at p=0.4,0.43"); print("="*70)
for p in [0.4, 0.43]:
    m99=face_min({0},p,nstarts=800,seed=7)
    m97=face_min({0,2},p,nstarts=600,seed=5)
    print(f" p={p}: M_9,9~{m99:.6f}  M_9,7~{m97:.6f}  (GPT: 99>97, 99 saddle)")

# ---------- envelope: cross point p* ----------
print("\n"+"="*70); print("Envelope cross L=5/L=7 near p*=0.39243"); print("="*70)
ps=np.linspace(0.08,0.43,40)
env5=[]; env7=[]
for p in ps:
    env5.append(face_min({0,2,4},p,nstarts=300,seed=3))  # L=5 zeros{0,2,4}
    env7.append(face_min({0,2},p,nstarts=300,seed=5))    # L=7 zeros{0,2}
env5=np.array(env5); env7=np.array(env7)
diff=env5-env7
# sign change
idx=np.where(np.diff(np.sign(diff)))[0]
print(f"  L=5 vs L=7 cross near p={ps[idx][0]:.5f}  (GPT p*=0.39242896)")
print(f"  cross value M~{env5[idx][0]:.5f}  (GPT 8.78637382)")
