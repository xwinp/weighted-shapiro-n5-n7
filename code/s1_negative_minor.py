#!/usr/bin/env python3
"""
# !! DIAGONAL BUG (2026-07-29): the 4-point mixed-difference H_ij below degenerates
# to f'/h noise for i=j (diagonal). Diagonal entries are WRONG. Use code/sym_hess3.py
# (symbolic) or code/det_scan2.py (three-point diagonal) instead. Off-diagonal entries are OK.
# Corrected result: S1 Morse=1, S0-nonunif Morse=1 (not 3/4); cert target is det(H)<0.
Route A Phase 1 -- numerical discovery step.
At the S1={0} stationary branch (x=(0,1,b,c,d,e,f), x1=1, 5 vars, degree-0 so
this 5x5 Hessian IS the scaling-quotient reduced Hessian), compute the 5x5 Hessian
and find a PRINCIPAL MINOR that is strictly negative across the whole band.
That minor is the target for the symbolic resultant+Sturm certificate
(prove C_S1 ∩ band ∩ {H_S1 ⪰ 0} = ∅ via one negative minor).
"""
import numpy as np
from scipy.optimize import root
n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s

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
        return r.x  # (b,c,d,e,f) with x1=1
    return None

def hess5(v, p):
    # 5x5 Hessian of P in vars (b,c,d,e,f) = x[2..6], x0=0,x1=1
    H=np.zeros((5,5)); h=1e-5
    def f(vv):
        x=np.zeros(n); x[1]=1
        for j in range(5): x[2+j]=vv[j]
        return Pval(x,p)
    for i in range(5):
        for j in range(5):
            vp=v.copy(); vm=v.copy(); vpp=v.copy(); vmm=v.copy()
            vp[i]+=h; vp[j]+=h; vm[i]-=h; vm[j]-=h
            vpp[i]+=h; vpp[j]-=h; vmm[i]-=h; vmm[j]+=h
            H[i,j]=(f(vp)-f(vm)-f(vpp)+f(vmm))/(4*h*h)
    return H

from itertools import combinations
a7,b7=0.214273520909841,0.328627677916592
v=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
# collect Hessian minors across band
ps=list(np.linspace(0.25,b7-0.004,12))[1:]+list(np.linspace(0.25,a7+0.004,12))[1:][::-1]
diag_signs=np.zeros(5)  # track min of each diagonal entry
minor2_negcount={}
samples=[]
for pp in ps:
    sol=solve_S1(pp,v)
    if sol is None: continue
    v=sol
    H=hess5(sol,pp)
    samples.append((pp,H))
    # diagonal entries
    for i in range(5):
        diag_signs[i]=min(diag_signs[i], H[i,i]) if samples!=[(pp,H)] else H[i,i]
# recompute diag min properly
diag_min=np.full(5, np.inf); diag_max=np.full(5, -np.inf)
for pp,H in samples:
    for i in range(5):
        diag_min[i]=min(diag_min[i], H[i,i]); diag_max[i]=max(diag_max[i], H[i,i])
print("=== diagonal entries H_ii (vars b,c,d,e,f) over band: min / max ===")
for i in range(5):
    print(f"  H[{i},{i}]  min={diag_min[i]:+.5f}  max={diag_max[i]:+.5f}  {'NEG THROUGHOUT -> certificate target!' if diag_max[i]<0 else ('sign-changing' if diag_min[i]<0 else 'positive')}")

# 2x2 principal minors: det of [[H_ii,H_ij],[H_ji,H_jj]]
print("\n=== 2x2 principal minors (i<j), min over band ===")
minor2_min={}
for i,j in combinations(range(5),2):
    minor2_min[(i,j)]=np.inf
for pp,H in samples:
    for i,j in combinations(range(5),2):
        m=H[i,i]*H[j,j]-H[i,j]**2
        minor2_min[(i,j)]=min(minor2_min[(i,j)], m)
cand=[(k,m) for k,m in minor2_min.items() if m<0]
cand.sort(key=lambda x:x[1])
for (i,j),m in cand[:8]:
    print(f"  minor[{i},{j}] min = {m:+.6f}  NEG candidate")
if not cand:
    print("  no 2x2 minor negative; check 3x3")

# eigenvalues at band center for reference
pp_mid=0.27; sol=solve_S1(pp_mid, v)
if sol is not None:
    H=hess5(sol,pp_mid)
    print(f"\n  at p={pp_mid}: Hessian eigs={np.round(np.linalg.eigvalsh(H),5)}  (neg count = Morse index)")
    print(f"  diagonal = {np.round(np.diag(H),5)}")
