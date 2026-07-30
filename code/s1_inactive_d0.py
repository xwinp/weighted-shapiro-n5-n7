#!/usr/bin/env python3
"""
S1 inactive-KKT test (analog of the S2 D0,D5 certificate).
At an S1={0} stationary point x* (x0=0, grad=0 on support {1..6}, Euler=>lambda=0
so dP/dx_i=0 for i=1..6), the directional derivative moving into S0 (increase x0,
decrease others to stay on simplex) equals  D0 = dP/dx0  evaluated at x0=0.
  D0 < 0  =>  P decreases moving into S0  =>  x* is a SADDLE (first-order cert).
  D0 > 0  =>  need second order (Hessian).
Compute D0 along the S1 stationary branch across (a7,b7).
Also compute D0 for S0 non-uniform: there x0>0 so "inactive" doesn't apply; instead
check the full reduced Hessian sign (already done in s0_s1_branch.py: Morse 4).
Here focus: S1 D0 sign across band.
"""
import numpy as np
from scipy.optimize import root
n = 7

def Pval(x, p):
    q = 1-p; s = 0.0
    for i in range(n):
        den = p*x[(i+1)%n] + q*x[(i+2)%n]
        if abs(den) < 1e-15: return 1e6
        s += x[i]/den
    return s

def dP_dx0(x, p):
    """partial derivative of P wrt x0, all else fixed. x0 appears in:
    term i=5: x5/(p*x6 + q*x0)  -> d/dx0 = x5 * q / (p*x6+q*x0)^2
    term i=6: x6/(p*x0 + q*x1)  -> d/dx0 = -x6 * p / (p*x0+q*x1)^2
    term i=0: x0/(p*x1+q*x2)    -> d/dx0 = 1/(p*x1+q*x2)
    (x0 also: term i where x0 is in denominator: i=5 (q*x0), i=6 (p*x0). and numerator term i=0.)
    """
    q = 1-p
    d = 0.0
    # term i=0: x0/(p*x1+q*x2)
    d += 1.0/(p*x[1] + q*x[2])
    # term i=5: x5/(p*x6 + q*x0)  -> d/dx0 = -x5*q/(den5)^2  (denominator increases)
    den5 = p*x[6] + q*x[0]
    d -= x[5]*q/(den5*den5)
    # term i=6: x6/(p*x0 + q*x1)
    den6 = p*x[0] + q*x[1]
    d -= x[6]*p/(den6*den6)
    return d

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
        return x, r.x
    return None,None

print("=== S1 D0 = dP/dx0 at x0=0, along S1 branch across band ===")
a7,b7=0.214273520909841,0.328627677916592
v=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
print("  p        P_S1stat   D0=dP/dx0   (D0<0 => SADDLE via descent into S0)")
for pp in list(np.linspace(0.25, b7-0.004, 16))[1:] + list(np.linspace(0.25, a7+0.004, 16))[1:][::-1]:
    x,vx = solve_S1(pp, v)
    if x is None: print(f"  p={pp:.4f} no conv"); continue
    v = vx
    xs = x/x.sum()
    Pv = Pval(xs, pp)
    D0 = dP_dx0(xs, pp)
    # cross-check with finite diff: increase x0 by h, rescale
    h=1e-7; xs2=xs.copy(); xs2[0]+=h; xs2=xs2/xs2.sum()
    D0_fd = (Pval(xs2,pp)-Pv)/ (h/xs.sum())  # derivative along rescale-on-simplex dir
    print(f"  p={pp:.4f}  P={Pv:.6f}  D0={D0:+.5f}  D0_fd={D0_fd:+.5f}  -> {'SADDLE(D0<0)' if D0<-1e-6 else ('min-dir' if D0>1e-6 else '0')}")
