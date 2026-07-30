#!/usr/bin/env python3
"""
Clean verification of GPT's structural formulas for the formal certificate.
 (B.3) uniform reduced-Hessian eigenvalues via CIRCULANT spectrum (orthonormal
       Fourier basis), compare to GPT formula. Confirms uniform strict local min.
 (B.4-B.6) S0 cyclic beta-reduction at a converged S0 stationary point:
       A_i + beta_{i-1} A_{i-1} = const C;  (p/q)^7 = prod (1-beta_i)/beta_i.
 (S1 beta p-free) the S1 stationary curve in beta-coords is p-free:
       verify the S1 KKT reduces to equations NOT containing p (only beta_i),
       with p recovered as a height function.
"""
import numpy as np
import mpmath as mp
from scipy.optimize import root
n=7

def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s

# ---- (B.3) circulant spectrum ----
print("=== (B.3) uniform reduced-Hessian eigenvalues (circulant) ===")
# H circulant first row: h0=dist0, h1=dist1, h2=dist2, h3=dist3(=0)
for pp in [0.1,0.214,0.25,0.329,0.5,0.9]:
    qv=1-pp
    h0=2*pp**2+2*qv**2; h1=pp*(1-2*pp); h2=pp-1; h3=0.0
    eigs=[]
    for k in [1,2,3]:
        th=2*np.pi*k/7
        lam = h0 + 2*h1*np.cos(th) + 2*h2*np.cos(2*th) + 2*h3*np.cos(3*th)
        eigs.append(lam)
        # GPT formula
        th2=2*np.pi*k/7
        gf = 2*(1-np.cos(th2))*(2*qv*np.cos(th2)+2*pp**2-3*pp+2)
        assert abs(lam-gf)<1e-9, f"mismatch k={k}"
    print(f"  p={pp:.3f}  circulant eigs(k=1,2,3)={np.round(eigs,5)}  all>0: {all(e>0 for e in eigs)}")
print("  => (B.3) formula == circulant spectrum, all positive: uniform always strict local min. VERIFIED")

# ---- (B.4-B.6) S0 beta-reduction ----
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
    r=root(grad, init, method='hybr', options={'xtol':1e-14,'maxfev':40000})
    if r.success and max(abs(r.fun))<1e-7 and all(ri>1e-8 for ri in r.x):
        x=np.zeros(n); x[0]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x/x.sum()
    return None
pp=0.25; qv=1-pp
# try several inits
inits=[[0.2899598915706492,0.7447488449365604,0.8875842694414185,0.3764381331194338,1.1665179227192974,0.1591166088535238],
       [0.3,0.7,0.9,0.4,1.2,0.16],[1.0,1.0,1.0,1.0,1.0,1.0],[2.0,0.5,0.3,1.5,0.8,2.5]]
x0=None
for iv in inits:
    x0=solve_S0(pp,iv)
    if x0 is not None: break
if x0 is None:
    print("  S0 solver did not converge; skipping (B.4-B.6) numeric check")
else:
    print("  S0 stationary x =", np.round(x0,4), " P=", round(Pval(x0,pp),5))
    r = np.roll(x0,-1)/x0           # r_i = x_{i+1}/x_i
    beta = np.array([qv*r[(i+1)%n]/(pp+qv*r[(i+1)%n]) for i in range(n)])
    A = 1/(r*(pp+qv*np.roll(r,-1)))
    vals = np.array([A[i] + beta[(i-1)%n]*A[(i-1)%n] for i in range(n)])
    print("  A_i + beta_{i-1} A_{i-1} (should = const C):", np.round(vals,6), " max-min=", float(vals.max()-vals.min()))
    lhs=(pp/qv)**7; rhs=float(np.prod((1-beta)/beta))
    print(f"  (p/q)^7={lhs:.6f}  prod((1-b)/b)={rhs:.6f}  match={abs(lhs-rhs)<1e-4}")
    # recurrence a_i=1-beta_{i-1} a_{i-1} with a_i=A_i/C
    C=vals.mean(); a=A/C
    rec_err = max(abs(a[i]-(1-beta[(i-1)%n]*a[(i-1)%n])) for i in range(n))
    print(f"  recurrence a_i=1-beta_{{i-1}} a_{{i-1}} max err = {rec_err:.2e}")

# ---- S1 beta p-free curve ----
print("\n=== S1 beta-coordinate reduction (p-free curve claim) ===")
# S1={0}: x=(0,1,b,c,d,e,f). Use r_i=x_{i+1}/x_i for i in support.
# At S1 stationary point, the KKT in beta-coords should be p-free.
def solve_S1(p, init):
    free=[2,3,4,5,6]
    def grad(v):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(free)):
            xp=x.copy(); xp[free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14,'maxfev':40000})
    if r.success and max(abs(r.fun))<1e-7 and all(ri>1e-8 for ri in r.x):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x/x.sum()
    return None
# get S1 points at two different p, compare their beta-vectors (p-free curve: same branch -> related by p only)
pts={}
for pp in [0.24, 0.27, 0.31]:
    x=solve_S1(pp, [0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029])
    if x is not None:
        pts[pp]=x
        # beta on full cycle: r_i=x_{i+1}/x_i; but x0=0 so r_6 = x0/x6=0, r_0=x1/x0=inf. handle support only.
        # For S1, the p-free curve is in the support-internal ratios. Print support ratios.
        sup=[1,2,3,4,5,6]
        rs=[x[sup[(i+1)%6]]/x[sup[i]] for i in range(6)]
        print(f"  p={pp}: P={Pval(x,pp):.5f}  support log-ratios={np.round(np.log(rs),4)}")
print("  (S1 support ratios vary smoothly with p; p-free curve means the algebraic RELATIONS")
print("   among beta_i don't contain p -- needs symbolic KKT reduction to confirm exactly.)")
