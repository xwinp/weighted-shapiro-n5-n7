#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent rigor for p0 (P=9 crossing on {0}-face palindromic branch):
solve the 4-eqn system [stationarity (3, by reflection symmetry) + P=9] in mpmath
high precision, then PSLQ to recover the minimal polynomial of p0 (and C0).
Bypasses GPT's Q97 entirely. If PSLQ finds a modest-degree minimal polynomial,
Sturm verifies the unique root in the band -> p0 rigorous, independent of GPT."""
import numpy as np
from scipy.optimize import minimize, root
import mpmath as mp
import sys, warnings
warnings.filterwarnings('ignore')
mp.mp.prec = 200  # ~60 decimal digits

def Px_np(x,p):
    q=1-p; return sum(x[i]/(p*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9))
def grad_P_np(x,p):
    q=1-p; g=np.zeros(9)
    for i in range(9):
        Dj=p*x[(i+1)%9]+q*x[(i+2)%9]; g[i]+=1.0/Dj
        for l in range(9):
            Dl=p*x[(l+1)%9]+q*x[(l+2)%9]
            if (l+1)%9==i: g[i]-=p*x[l]/Dl**2
            if (l+2)%9==i: g[i]-=q*x[l]/Dl**2
    return g

# --- get initial guess (double precision) ---
def get_guess(p):
    def F(v):
        x=np.array([0.0]+list(np.exp(v))); g=grad_P_np(x,p); return np.sum(g[1:]**2)
    rng=np.random.default_rng(1)
    for _ in range(500):
        v0=np.log(rng.dirichlet(np.ones(8))+1e-9)
        r=minimize(F,v0,method='Nelder-Mead',options={'maxiter':5000,'xatol':1e-12,'fatol':1e-16})
        if r.fun<1e-14:
            x=np.array([0.0]+list(np.exp(r.x)))
            # polish
            def G(v):
                xx=np.array([0.0,1.0]+list(v)); g=grad_P_np(xx,p); return g[2:]
            sol=root(G, x[2:9], method='hybr', options={'xtol':1e-14,'maxfev':8000})
            if sol.success:
                xx=np.array([0.0,1.0]+list(sol.x))
                if np.sum(grad_P_np(xx,p)[2:]**2)<1e-20 and 8<Px_np(xx,p)<9.05:
                    return xx
    return None

print("finding initial guess near p0...", flush=True)
# find saddle near crossing; bisect p where P=9
def saddleP(p):
    x=get_guess(p)
    return (Px_np(x,p) if x is not None else None), x

# bracket crossing
lo,hi=0.40,0.4318363763
x_lo=get_guess(lo); x_hi=get_guess(hi)
print(f"  P({lo})={Px_np(x_lo,lo):.6f}, P({hi})={Px_np(x_hi,hi):.6f}", flush=True)

# --- high-precision solve: 4 eqns (stat x2,x3,x4 =0  [reflection-symmetric] + P=9), 4 unknowns (x2,x3,x4,p) ---
# x palindromic: x0=0, x1=1(gauge), x2,x3,x4,x5=x4,x6=x3,x7=x2,x8=1
def build(x2,x3,x4,p):
    x=[mp.mpf(0),mp.mpf(1),x2,x3,x4,x4,x3,x2,mp.mpf(1)]
    return x
def P_mp(x,p):
    q=1-p; return sum(x[i]/(p*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9))
def grad_mp(x,p):
    q=1-p; g=[mp.mpf(0)]*9
    for i in range(9):
        Dj=p*x[(i+1)%9]+q*x[(i+2)%9]; g[i]+=1/Dj
        for l in range(9):
            Dl=p*x[(l+1)%9]+q*x[(l+2)%9]
            if (l+1)%9==i: g[i]-=p*x[l]/Dl**2
            if (l+2)%9==i: g[i]-=q*x[l]/Dl**2
    return g

def system(x2,x3,x4,p):
    x=build(x2,x3,x4,p)
    g=grad_mp(x,p)
    # stationarity for x2,x3,x4 (gauge x1 fixed; symmetry handles x5..x8)
    P=P_mp(x,p)
    return [g[2], g[3], g[4], P-9]

# initial guess from double-precision saddle at hi (P~9)
xg=x_hi
x2g,x3g,x4g,pg=xg[2],xg[3],xg[4],0.4318363763
print(f"initial guess: x2={x2g:.6f} x3={x3g:.6f} x4={x4g:.6f} p={pg:.10f}", flush=True)

# refine in double precision first: 4-eqn scipy root
def system_np(v):
    x2,x3,x4,p=v
    x=np.array([0.0,1.0,x2,x3,x4,x4,x3,x2,1.0])
    g=grad_P_np(x,p); P=Px_np(x,p)
    return [g[2],g[3],g[4],P-9]
print("scipy refine of 4-eqn system...", flush=True)
sr=root(system_np, [x2g,x3g,x4g,pg], method='hybr', options={'xtol':1e-13,'maxfev':20000})
print(f"  scipy residual={np.linalg.norm(system_np(sr.x)):.3e}, sol={sr.x}", flush=True)
x2g,x3g,x4g,pg=[float(v) for v in sr.x]

print("high-precision solve...", flush=True)
try:
    # gradual precision: start at low dps with loose tol, then refine
    mp.mp.prec=120
    sol=mp.findroot(system, [mp.mpf(x2g),mp.mpf(x3g),mp.mpf(x4g),mp.mpf(pg)], tol=mp.mpf(10)**(-30), maxsteps=400)
    x2s,x3s,x4s,ps=[mp.mpf(s) for s in sol]
    mp.mp.prec=220
    sol=mp.findroot(system, [x2s,x3s,x4s,ps], tol=mp.mpf(10)**(-60), maxsteps=400)
    x2s,x3s,x4s,ps=[mp.mpf(s) for s in sol]
    x2s,x3s,x4s,ps=sol
    x=build(x2s,x3s,x4s,ps)
    Pval=P_mp(x,ps)
    print(f"  SOLVED: p0={mp.nstr(ps,50)}", flush=True)
    print(f"          P={mp.nstr(Pval,30)} (expect 9)", flush=True)
    # C0 = T4/T1
    q=1-ps; T=[x[i]/(ps*x[(i+1)%9]+q*x[(i+2)%9]) for i in range(9)]
    C0=T[3]/T[1]; A0=T[2]/T[1]; B0=T[4]/T[1]  # wait T4=T5=lam*C; T1=lam
    C0v=T[4]/T[1]
    print(f"  A0={mp.nstr(A0,30)} B0={mp.nstr(B0,30)} C0={mp.nstr(C0v,30)}", flush=True)
    # save high-precision p0, C0
    open('paper/_p0_highprec.txt','w').write(f"p0={mp.nstr(ps,60)}\nC0={mp.nstr(C0v,60)}\nA0={mp.nstr(A0,40)}\nB0={mp.nstr(B0,40)}\n")
    # --- PSLQ for minimal polynomial of p0 ---
    print("\n=== PSLQ for minimal polynomial of p0 ===", flush=True)
    for deg in range(6, 56):
        powers=[ps**k for k in range(deg+1)]
        rel=mp.pslq(powers)
        if rel is not None:
            # verify
            val=sum(c*ps**k for k,c in enumerate(rel))
            if abs(val)<mp.mpf(10)**(-40):
                print(f"  degree {deg}: coeffs={rel}", flush=True)
                print(f"  residual={mp.nstr(val,5)}", flush=True)
                # leading coeff should be nonzero; print poly
                poly_str=" + ".join(f"({c})*p^{k}" for k,c in enumerate(rel) if c!=0)
                print(f"  poly: {poly_str}", flush=True)
                break
    else:
        print("  no minimal poly found up to degree 55 (p0 may be high-degree)", flush=True)
    # --- PSLQ for C0 ---
    print("\n=== PSLQ for minimal polynomial of C0 ===", flush=True)
    for deg in range(6, 56):
        powers=[C0v**k for k in range(deg+1)]
        rel=mp.pslq(powers)
        if rel is not None:
            val=sum(c*C0v**k for k,c in enumerate(rel))
            if abs(val)<mp.mpf(10)**(-40):
                print(f"  degree {deg}: coeffs={rel}", flush=True)
                print(f"  poly: {' + '.join(f'({c})*C^{k}' for k,c in enumerate(rel) if c!=0)}", flush=True)
                break
    else:
        print("  no minimal poly found up to degree 55", flush=True)
except Exception as e:
    print(f"  high-precision solve failed: {e}", flush=True)
    import traceback; traceback.print_exc()
