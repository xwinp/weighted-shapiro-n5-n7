#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extend the S1 det<0 certificate to the FULL p in (0,1).
Route A: on the H_B branch, det(H_S1)=0  <=>  Q5(z)=0 or Q7(z)=0  (z=last beta ratio).
In-band z in (0.9488,0.9911) subset (9/10,1) where Q5,Q7 have 0 roots.
Here: (a) all real roots of Q5,Q7; (b) z-range of the full S1 branch over p in (0,1);
(c) check Q5,Q7 have no root in the full z-range => det != 0 on whole branch =>
Morse constant on whole branch => det<0 everywhere (sample neg) => S1 saddle for all p.
"""
import numpy as np
from scipy.optimize import root
import sympy as sp
n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s
sup=[1,2,3,4,5,6]; free=[2,3,4,5,6]
def solve_S1(p, init):
    def grad(v):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(free)):
            xp=x.copy(); xp[free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14})
    if r.success and max(abs(r.fun))<1e-7 and all(r.x>1e-8):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(free): x[idx]=r.x[j]
        return x
    return None

# (a) all real roots of Q5,Q7
z=sp.symbols('z')
Q5=2*z**5+2*z**3-2*z**2-1
Q7=8*z**7-24*z**6+20*z**5-9*z**4+30*z**3-15*z**2-6
r5=[float(r) for r in sp.real_roots(sp.Poly(Q5,z))]
r7=[float(r) for r in sp.real_roots(sp.Poly(Q7,z))]
print("Q5 real roots:", r5)
print("Q7 real roots:", r7)
# how many in (0,1)?
print("Q5 roots in (0,1):", [r for r in r5 if 0<r<1])
print("Q7 roots in (0,1):", [r for r in r7 if 0<r<1])

# (b) trace full branch, compute z = q*rho4/(p+q*rho4), rho=[b, c/b, d/c, e/d, f/e]
# x = (0,1,b,c,d,e,f)  -> indices: x1=1,x2=b,x3=c,x4=d,x5=e,x6=f
inits=[]
rng=np.random.default_rng(0)
for _ in range(40): inits.append(rng.uniform(0.05,2.5,size=5))
inits.append([0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029])
inits+=[[0.3,0.7,1.5,0.07,1.3],[1,1,1,1,1],[0.5,0.5,1.0,0.5,0.5]]
zs=[]
ps=np.linspace(0.02,0.98,97)
prev=None
for pp in ps:
    candidates=[]
    init_list = [prev] if prev is not None else []
    init_list += inits[:20]
    for init in init_list:
        x=solve_S1(pp,init)
        if x is None: continue
        xs=x/x.sum()
        Pv=Pval(xs,pp)
        if any(abs(Pv-c[0])<1e-4 for c in candidates): continue
        candidates.append((Pv,xs))
    if candidates:
        # pick the one closest to 7 (the main branch; avoids spurious far ones)
        candidates.sort(key=lambda c:abs(c[0]-7))
        Pv,xs=candidates[0]
        prev=[xs[2],xs[3],xs[4],xs[5],xs[6]]  # b,c,d,e,f as next seed (un-normalized-ish)
        b=xs[2];c=xs[3];d=xs[4];e=xs[5];f=xs[6]
        rho=[b,c/b,d/c,e/d,f/e]
        q=1-pp
        z0=q*rho[4]/(pp+q*rho[4])
        zs.append((pp,z0,Pv))
zs_arr=np.array([t[1] for t in zs])
print("\nbranch z-range over p in (0.02,0.98): z in (%.6f, %.6f)"%(zs_arr.min(),zs_arr.max()))
print("z at p=0.02:%.5f  p=0.214(a7):%.5f  p=0.27:%.5f  p=0.329(b7):%.5f  p=0.98:%.5f"%(
    zs[0][1], [t for t in zs if abs(t[0]-0.22)<0.03][0][1] if any(abs(t[0]-0.22)<0.03 for t in zs) else -1,
    [t for t in zs if abs(t[0]-0.27)<0.02][0][1] if any(abs(t[0]-0.27)<0.02 for t in zs) else -1,
    [t for t in zs if abs(t[0]-0.33)<0.03][0][1] if any(abs(t[0]-0.33)<0.03 for t in zs) else -1,
    zs[-1][1]))
zlo,zhi=zs_arr.min(),zs_arr.max()
# (c) Q5,Q7 roots in (zlo,zhi)?
in5=[r for r in r5 if zlo-1e-9<r<zhi+1e-9]
in7=[r for r in r7 if zlo-1e-9<r<zhi+1e-9]
print("\nQ5 roots in branch z-range (%.5f,%.5f): %s  -> %s"%(zlo,zhi,in5,"NONE => det!=0" if not in5 else "GAP"))
print("Q7 roots in branch z-range (%.5f,%.5f): %s  -> %s"%(zlo,zhi,in7,"NONE => det!=0" if not in7 else "GAP"))
# sample det sign at an out-of-band point (p=0.5)
Pcoeff=[-8,100,-528,1564,-2972,4112,-4806,5120,-4696,3540,-2298,1302,-576,188,-52,10] # z^16..z^1
Qcoeff=[6,-48,194,-568,1332,-2536,4056,-5644,6786,-6968,6270,-5016,3340,-1648,536,-100,8] # z^0..z^16
def num_red(zv,wv):
    Pv=sum(Pcoeff[i]*zv**(16-i) for i in range(16)) # z^16 down to z^1
    Qv=sum(Qcoeff[i]*zv**i for i in range(17))
    return Pv*wv+Qv
# at p=0.5 out-of-band
cand=[t for t in zs if abs(t[0]-0.5)<0.03]
if cand:
    pp0,z0,Pv0=cand[0]
    # need w too; recompute w = q*rho3/(p+q*rho3)
    # re-solve to get x
    x=solve_S1(0.5,[0.4,0.6,1.2,0.3,1.0])
    if x is not None:
        xs=x/x.sum(); b,c,d,e,f=xs[2],xs[3],xs[4],xs[5],xs[6]
        rho=[b,c/b,d/c,e/d,f/e]; q=0.5
        w0=q*rho[3]/(0.5+q*rho[3]); z0=q*rho[4]/(0.5+q*rho[4])
        db=num_red(z0,w0)/(1-z0**2+z0**2*w0)**5
        print("out-of-band p=0.5: z=%.5f w=%.5f  D_B=num_red/den=%.4e  -> %s"%(z0,w0,db,"<0 (saddle)" if db<0 else ">0"))
print("DONE")
