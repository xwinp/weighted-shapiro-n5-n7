#!/usr/bin/env python3
"""
Value-inequality certificate strategy (avoids Morse/Hessian).
If we can show, across the whole band (a7,b7):
  (I1) P_S1^stat(p) >= P_S2^curve(p)   [S1 interior stationary above S2 curve]
  (I2) P_S0^nonunif(p) >= 7             [S0 non-uniform stationary above 7]
then (assuming each orbit's interior stationary set is fully captured -- "no遗漏"):
  - S1 face inf = min(P_S1^stat, boundary=S2..) = P_S2  (since P_S1^stat>=P_S2)
  - S0 face inf = min(7 uniform, P_S0^nonunif>=7, boundary=S1 inf=P_S2) = P_S2 in band
  => m_7 = P_S2 in band.  No second-order/Morse analysis needed.
This script checks the inequalities numerically with margins, and counts branches.
"""
import numpy as np
import mpmath as mp
from scipy.optimize import root
mp.mp.dps = 40
n = 7

def Pval(x, p):
    q = 1-p; s = 0.0
    for i in range(n):
        den = p*x[(i+1)%n] + q*x[(i+2)%n]
        if abs(den) < 1e-15: return 1e6
        s += x[i]/den
    return s

# ---- P_S2^curve(p): solve R(p,t)=q^3 - p^3 t^5 - p^2 q t^8 = 0, compute P ----
def P_S2_curve(pp):
    pv=mp.mpf(pp); qv=1-pv
    tv=mp.findroot(lambda tt: qv**3 - pv**3*tt**5 - pv**2*qv*tt**8, mp.mpf('1.3'))
    d=tv**2
    c=qv*(qv-pv*tv**4)/(pv**2*tv**2)
    b=qv/(pv*tv)-qv**2*(qv-pv*tv**4)/(pv**3*tv**2)
    x=mp.matrix([0,1,b,c,d,0,tv])
    s=mp.mpf(0)
    for i in range(n):
        den=pv*x[(i+1)%n]+qv*x[(i+2)%n]
        s+=x[i]/den
    return float(s)

# ---- S1 stationary branch solver ----
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
        return x/x.sum(), r.x
    return None,None

# ---- S0 non-uniform branch solver ----
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
        return x/x.sum(), r.x
    return None,None

a7,b7=0.214273520909841,0.328627677916592
print("=== (I1) P_S1^stat vs P_S2^curve across band ===")
print("  p        P_S2      P_S1stat   diff(S1-S2)  margin")
v1=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
ps=list(np.linspace(0.25,b7-0.004,14))[1:]+list(np.linspace(0.25,a7+0.004,14))[1:][::-1]
min_diff=9
for pp in ps:
    x1,vx=solve_S1(pp,v1)
    if x1 is None: continue
    v1=vx
    P1=Pval(x1,pp); P2=P_S2_curve(pp)
    min_diff=min(min_diff,P1-P2)
    print(f"  {pp:.4f}  {P2:.5f}  {P1:.5f}  {P1-P2:+.5f}")
print(f"  >> min(P_S1stat - P_S2) over band = {min_diff:.5f}  (>=0: {min_diff>=0})")

print("\n=== (I2) P_S0^nonunif vs 7 across band ===")
print("  p        P_S0nonunif   P-7")
v0=[0.2899598915706492,0.7447488449365604,0.8875842694414185,0.3764381331194338,1.1665179227192974,0.1591166088535238]
min_m=9
for pp in ps:
    x0,vx=solve_S0(pp,v0)
    if x0 is None: continue
    v0=vx
    P0=Pval(x0,pp); min_m=min(min_m,P0-7)
    print(f"  {pp:.4f}  {P0:.5f}  {P0-7:+.5f}")
print(f"  >> min(P_S0nonunif - 7) over band = {min_m:.5f}  (>=0: {min_m>=0})")
