#!/usr/bin/env python3
"""
n=7 COMPREHENSIVE orbit scan. Enumerate all 2^7 subsets as zero-sets, group by
dihedral orbit (C_7 + reflection), pick a representative, SLSQP-minimize P over
the support (zeros=0, free>=1e-6) at band center p=0.27. Detect whether the
minimizer is a GENUINE interior point of that support (all free coords > 1e-3)
or DEGENERATES (some free coord hits ~0 -> belongs to smaller support).
Report each orbit's min P and degeneracy status. Goal: confirm only the 2z-d2
orbit ({0,5}) has interior min P<7.
"""
import numpy as np
from scipy.optimize import minimize
from itertools import combinations

n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if den<=1e-15: return 1e6
        s+=x[i]/den
    return s

def orbit_key(zeros):
    z=frozenset(zeros)
    reps=[]
    for rot in range(n):
        r=frozenset((i+rot)%n for i in z)
        reps.append(tuple(sorted(r)))
        reps.append(tuple(sorted((n-i)%n for i in r)))  # reflection
    return min(reps)

def min_orbit(zeros,p,nstarts=16,lo=1e-6):
    free=[i for i in range(n) if i not in zeros]
    if not free: return None,None
    k=len(free); best=1e9; bx=None
    rng=np.random.RandomState(0)
    for s in range(nstarts):
        y=rng.rand(k)+0.05; y/=y.sum()
        x0=np.zeros(n)
        for j,idx in enumerate(free): x0[idx]=y[j]
        bounds=[(lo,1.0) if i in free else (0.0,0.0) for i in range(n)]
        cons=[{'type':'eq','fun':lambda x:x.sum()-1.0}]
        try:
            r=minimize(lambda x:Pval(x,p),x0,method='SLSQP',bounds=bounds,constraints=cons,options={'maxiter':500,'ftol':1e-13})
            if r.fun<best: best=r.fun; bx=r.x
        except Exception: pass
    return best,bx

# enumerate all orbits
seen={}
for k in range(0,5):  # 0..4 zeros (5+ zeros -> <=2 nonzeros, P huge, skip)
    for z in combinations(range(n),k):
        key=orbit_key(z)
        if key not in seen:
            seen[key]=set(z)

p=0.27
print(f"p={p}. {'orbit(zeros)':24s} {'#nz':>4s} {'minP':>9s} {'genuine?':>9s} {'extraZeros':>22s}")
results=[]
for key,zset in sorted(seen.items(), key=lambda kv:(len(kv[0]),kv[0])):
    nz=n-len(zset)
    if nz<=1: continue
    v,x=min_orbit(zset,p)
    if v is None: continue
    free=[i for i in range(n) if i not in zset]
    # genuine interior = all free coords clearly >0
    extra=sorted([i for i in free if x[i]<1e-3])
    genuine = (len(extra)==0)
    tag = "INTERIOR" if genuine else "degen->"
    zstr=','.join(map(str,sorted(zset))) or '-'
    print(f"  z={{{zstr:20s}}} {nz:4d} {v:9.5f} {tag:>9s} extra={extra}")
    results.append((zset,v,genuine,extra))

print("\nInterior-stationary orbits (genuine, not degenerate):")
for zset,v,g,ex in results:
    if g:
        zstr=','.join(map(str,sorted(zset))) or '-'
        flag = "  *** P<7 ***" if v<7-1e-4 else ""
        print(f"  z={{{zstr:20s}}} minP={v:.6f}{flag}")
