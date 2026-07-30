#!/usr/bin/env python3
"""
n=7 structural verification: trace the support-{1,2,3,4,6} critical curve x*(p)
across p in (0,1) by continuation. Confirm:
  (1) grad P = 0 on support (lambda=0 by Euler) at every p on the curve,
  (2) P_curve(p) < 7 exactly on (a_7, b_7) and >=7 outside,
  (3) P_curve matches the SLSQP global min inside the band.
Square system for fixed p: grad_i = lambda (i=a..e) [5] + sum=1 [1] = 6 eqs, 6 unk (a,b,c,d,e,lam).
"""
import sympy as sp
import numpy as np
from scipy.optimize import minimize

a,b,c,d,e,lam,p = sp.symbols('a b c d e lam p', positive=True)
q=1-p
P=(a/(p*b+q*c)+b/(p*c+q*d)+c/(p*d)+d/(q*e)+e/(q*a))
gs=[sp.diff(P,v) for v in (a,b,c,d,e)]
# fixed-p system: substitute p=p0
def solve_at(p0, init):
    p0=sp.Rational(p0).limit_denominator(10**12) if False else sp.Float(p0)
    eqs=[g.subs(p,p0)-lam for g in gs]+[a+b+c+d+e-1]
    vs=[a,b,c,d,e,lam]
    sol=sp.nsolve(eqs, vs, init, prec=40, tol=sp.Float('1e-30'), maxsteps=200)
    return [float(x) for x in sol]

def P_num(xv, pv):
    a,b,c,d,e=xv; q=1-pv
    return a/(pv*b+q*c)+b/(pv*c+q*d)+c/(pv*d)+d/(q*e)+e/(q*a)

# SLSQP global min for comparison
def P_float(x,n,pv):
    q=1-pv; s=0.0
    for i in range(n):
        dd=pv*x[(i+1)%n]+q*x[(i+2)%n]
        if dd<=1e-12: return 1e18
        s+=x[i]/dd
    return s
def min_true(n,pv,ns=200,seed=1):
    rng=np.random.default_rng(seed)
    cons=({'type':'eq','fun':lambda x:np.sum(x)-1},)
    bd=[(0,1)]*n; best=1e18; bx=None
    seeds=[np.full(n,1/n),np.array([1 if i%2==0 else .5 for i in range(n)],float)]
    k=0
    while k<ns:
        x0=seeds[k].copy() if k<len(seeds) else rng.random(n)
        x0=np.array(x0,float); x0/=x0.sum()
        try:
            r=minimize(P_float,x0,args=(n,pv),method='SLSQP',bounds=bd,constraints=cons,options={'maxiter':1000,'ftol':1e-14})
            if r.fun<best: best=r.fun; bx=r.x
        except: pass
        k+=1
    return best,bx

# continuation from mid-band
ps_grid=np.linspace(0.05,0.50,46)
init=[0.188,0.21,0.072,0.295,0.235,0.0]  # from phase2 mid, lam=0
curve=[]
print(f"{'p':>7} {'P_curve':>12} {'lam':>10} {'min_5var?':>9} {'SLSQP':>10} {'match':>6}")
prev=init
for pv in ps_grid:
    try:
        sol=solve_at(float(pv), prev)
        av,bv,cv,dv,ev,lamv=sol
        Pv=P_num([av,bv,cv,dv,ev], float(pv))
        # verify all 5 positive (support intact)
        pos = all(v>1e-6 for v in [av,bv,cv,dv,ev])
        curve.append((pv,Pv,lamv,pos))
        # SLSQP global
        ms,_=min_true(7,float(pv),ns=120,seed=2)
        match = abs(Pv-ms)<1e-3 if Pv<7 else (ms>=7-1e-3)
        print(f"{pv:7.3f} {Pv:12.5f} {lamv:10.2e} {str(pos):>9} {ms:10.5f} {str(match):>6}")
        prev=sol
    except Exception as ex:
        print(f"{pv:7.3f}  FAILED: {ex}")
        curve.append((pv,None,None,None))
