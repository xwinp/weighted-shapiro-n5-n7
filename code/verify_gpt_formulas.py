#!/usr/bin/env python3
"""
Verify GPT's key formulas for the formal certificate:
 (B.3) uniform reduced-Hessian eigenvalues:
   lambda_k = 2(1-cos theta_k)[2 q cos theta_k + 2 p^2 - 3 p + 2],  k=1,2,3
   theta_k = 2 pi k /7.  Claim: >0 for all p in (0,1) (uniform always strict local min).
 (B.1) second variation along log-direction x_i(s)=x_i e^{s v_i}:
   Q_x(v) = sum_i T_i [ (v_i - alpha_i v_{i+1} - beta_i v_{i+2})^2
                         - alpha_i beta_i (v_{i+1}-v_{i+2})^2 ]
   where d_i = p x_{i+1}+q x_{i+2}, T_i = x_i/d_i,
         alpha_i = p x_{i+1}/d_i, beta_i = q x_{i+2}/d_i (alpha+beta=1).
   Verify against numerical d^2/ds^2 P(x e^{s v}) at s=0.
 (B.4-B.6) S0 cyclic beta-reduction: a_i = 1 - beta_{i-1} a_{i-1},
   (p/q)^7 = prod (1-beta_i)/beta_i. Verify at a known S0 stationary point.
"""
import numpy as np
import mpmath as mp
from scipy.optimize import root
mp.mp.dps = 30
n = 7

def Pval(x, p):
    q = 1-p; s = 0.0
    for i in range(n):
        den = p*x[(i+1)%n] + q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s += x[i]/den
    return s

# ---- (B.3) uniform eigenvalues ----
print("=== (B.3) uniform reduced-Hessian eigenvalues ===")
pp_list = [0.1, 0.214, 0.25, 0.27, 0.329, 0.5, 0.7, 0.9]
for pp in pp_list:
    q = 1-pp
    x = np.ones(n)
    # numerical Hessian at uniform
    H = np.zeros((n,n)); h=1e-5
    for i in range(n):
        for j in range(n):
            xp=x.copy(); xm=x.copy(); xpp=x.copy(); xmm=x.copy()
            xp[i]+=h; xp[j]+=h; xm[i]-=h; xm[j]-=h
            xpp[i]+=h; xpp[j]-=h; xmm[i]-=h; xmm[j]+=h
            H[i,j]=(Pval(xp,pp)-Pval(xm,pp)-Pval(xpp,pp)+Pval(xmm,pp))/(4*h*h)
    # restrict to sum=0 subspace (Fourier modes), eigenvalues
    # tangent basis: e_i - e_0 for i=1..6
    B = np.zeros((n,6))
    for i in range(6): B[i+1,i]=1; B[0,i]=-1
    eigs_num = np.linalg.eigvalsh(B.T@H@B)
    # GPT formula (B.3) for k=1,2,3 (each twice, since real sin/cos pair)
    eigs_formula = []
    for k in [1,2,3]:
        th = 2*np.pi*k/7
        lam = 2*(1-np.cos(th))*(2*q*np.cos(th) + 2*pp**2 - 3*pp + 2)
        eigs_formula.append(lam)
    eigs_formula_sorted = sorted(np.repeat(eigs_formula,2))
    print(f"  p={pp:.3f}  num eigs={np.round(sorted(eigs_num),4)}")
    print(f"           formula  ={np.round(eigs_formula_sorted,4)}  all>0: {all(e>0 for e in eigs_formula)}")
# check discriminant of bracket 2q cos + 2p^2-3p+2 <0 claim: disc = 4cos^2-4cos-7
print("  bracket disc 4cos^2-4cos-7 for k=1,2,3:",
      [round(4*np.cos(2*np.pi*k/7)**2 - 4*np.cos(2*np.pi*k/7) - 7, 4) for k in [1,2,3]], "(all<0 -> bracket>0 always)")

# ---- (B.1) second variation Q_x(v) ----
print("\n=== (B.1) Q_x(v) vs numerical d^2/ds^2 P(x e^{s v}) ===")
# use the S1 stationary point at p=1/4 as a nontrivial x
def solve_S1(p, init):
    free=[2,3,4,5,6]
    def grad(v):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(free)):
            xp=x.copy(); xp[free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14,'maxfev':20000})
    if r.success and max(abs(r.fun))<1e-7 and all(ri>1e-8 for ri in r.x):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x/x.sum()
    return None
x = solve_S1(0.25, [0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029])
pp=0.25; q=1-pp
print("  using S1 stationary x =", np.round(x,4))
# T_i, alpha_i, beta_i
d = np.array([pp*x[(i+1)%n]+q*x[(i+2)%n] for i in range(n)])
T = x/d; alpha = pp*np.roll(x,-1)/d; beta = q*np.roll(x,-2)/d
print("  sum(alpha+beta) should be 1:", np.round(alpha+beta,6))
np.random.seed(3)
for trial in range(4):
    v = np.random.randn(n)
    # GPT formula
    Q_formula = 0.0
    for i in range(n):
        Q_formula += T[i]*((v[i]-alpha[i]*v[(i+1)%n]-beta[i]*v[(i+2)%n])**2 - alpha[i]*beta[i]*(v[(i+1)%n]-v[(i+2)%n])**2)
    # numerical d^2/ds^2 P(x e^{s v}) at s=0
    h=1e-6
    Pp = Pval(x*np.exp(h*v), pp); Pm = Pval(x*np.exp(-h*v), pp); P0 = Pval(x,pp)
    Q_num = (Pp - 2*P0 + Pm)/h**2
    print(f"  trial{trial}: Q_formula={Q_formula:+.6f}  Q_num={Q_num:+.6f}  match={abs(Q_formula-Q_num)<1e-3}")

# ---- (B.4-B.6) S0 beta-reduction at a known S0 stationary point ----
print("\n=== (B.4-B.6) S0 cyclic beta-reduction ===")
def solve_S0(p, init):
    free=[1,2,3,4,5,6]
    def grad(v):
        x=np.zeros(n); x[0]=1
        for j,idx in enumerate(free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(free)):
            xp=x.copy(); xp[free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14,'maxfev':20000})
    if r.success and max(abs(r.fun))<1e-7 and all(ri>1e-8 for ri in r.x):
        x=np.zeros(n); x[0]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x/x.sum()
    return None
x0 = solve_S0(0.25, [0.2899598915706492,0.7447488449365604,0.8875842694414185,0.3764381331194338,1.1665179227192974,0.1591166088535238])
print("  S0 non-uniform stationary x =", np.round(x0,4))
r = np.roll(x0,-1)/x0  # r_i = x_{i+1}/x_i
beta = q*r/ (pp + q*r)   # beta_i = q r_{i+1}/(p+q r_{i+1})? GPT: beta_i = q r_{i+1}/(p+q r_{i+1})
# GPT defines r_i=x_{i+1}/x_i, A_i=1/(r_i(p+q r_{i+1})), beta_i = q r_{i+1}/(p+q r_{i+1})
# so beta_i uses r_{i+1}: beta_i = q*r[(i+1)%n]/(p+q*r[(i+1)%n])
beta = np.array([q*r[(i+1)%n]/(pp+q*r[(i+1)%n]) for i in range(n)])
A = 1/(r*(pp+q*np.roll(r,-1)))
# recurrence a_i = 1 - beta_{i-1} a_{i-1}, with a_i = A_i/C. Check A_i + beta_{i-1} A_{i-1} = const
vals = np.array([A[i] + beta[(i-1)%n]*A[(i-1)%n] for i in range(n)])
print("  A_i + beta_{i-1} A_{i-1} (should be const C):", np.round(vals,6))
print("  const-ness (max-min):", float(vals.max()-vals.min()))
# closure (p/q)^7 = prod (1-beta_i)/beta_i
lhs = (pp/q)**7; rhs = np.prod((1-beta)/beta)
print(f"  (p/q)^7 = {lhs:.6f}   prod(1-beta_i/beta_i) = {float(rhs):.6f}  match={abs(lhs-rhs)<1e-4}")
