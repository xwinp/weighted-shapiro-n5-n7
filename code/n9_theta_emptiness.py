#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fast numerical emptiness check for {G=S=Theta=0} strict interior + P<9."""
import time, random
from pathlib import Path
import sympy as sp
import mpmath as mp
mp.mp.dps=15

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
X, Y, sigma = sp.symbols('X Y sigma')
C, D, sig = sp.symbols('C D sigma')

def load_small(name):
    s = sp.symbols("s")
    text = (HERE/name).read_text(encoding="utf-8").strip()
    return sp.Poly(sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma),
                   X, Y, sigma, domain=sp.ZZ)
G = load_small("nonpal_G_clean.txt"); S = load_small("nonpal_S_clean.txt")
Gcd = sp.expand(G.as_expr().subs({X: C+D, Y: C*D}))
Scd = sp.expand(S.as_expr().subs({X: C+D, Y: C*D}))
# load N,D for P<9 (rho^9=N/D) and stationary value
N = load_small("nonpal_rho9_num.txt"); Den = load_small("nonpal_rho9_den.txt")
Ncd = sp.expand(N.as_expr().subs({X: C+D, Y: C*D}))
Dcd = sp.expand(Den.as_expr().subs({X: C+D, Y: C*D}))

def center_lift(C, D, sigma):
    gap = C + D - 1
    a3 = 1 - C + sigma*C**2*(1-C)/gap
    a6 = 1 - D + sigma*D**2*(1-D)/gap
    a2 = 1 - a3 + sigma*a3**2*(1-a3)/(a3 + C - 1)
    a7 = 1 - a6 + sigma*a6**2*(1-a6)/(a6 + D - 1)
    return a2, a3, a6, a7
def full_center_lift(C, D, sigma):
    a2,a3,a6,a7 = center_lift(C, D, sigma)
    return (sp.Integer(1),a2,a3,C,D,a6,a7,sp.Integer(1))
def compact_theta_from_terms(a, sigma):
    (a1,a2,a3,a4,a5,a6,a7,a8)=a
    h2=a2+a3-1; h3=a3+a4-1; h4=a4+a5-1; h5=a5+a6-1; h6=a6+a7-1
    W23=((a6*a7/h6)*(sigma+h2*h6/(sigma*a2*a3*a6*a7)
        -h2*h3*h4*h5*h6/(sigma**4*a2*a3**2*a4**2*a5**2*a6**2*a7))
        *(1+a2*a3*h6/(h2*a6*a7)))
    W24=(h2*h3/(sigma**2*a2*a3**2*a4)*(h4/(sigma*a4*a5)-1)
        *(1+a2*a3**2*a4*h5*h6/(h2*h3*a5*a6**2*a7)))
    W34=((a5*a6/h5)*(sigma+h3*h5/(sigma*a3*a4*a5*a6)
        -h3*h4*h5/(sigma**2*a3*a4**2*a5**2*a6))
        *(1+a3*a4*h5/(h3*a5*a6)))
    return W23*W24+W23*W34+W24*W34
Theta = compact_theta_from_terms(full_center_lift(C,D,sig), sig)
fG=sp.lambdify((C,D,sig),Gcd,'mpmath'); fS=sp.lambdify((C,D,sig),Scd,'mpmath')
fT=sp.lambdify((C,D,sig),Theta,'mpmath')
fN=sp.lambdify((C,D,sig),Ncd,'mpmath'); fDd=sp.lambdify((C,D,sig),Dcd,'mpmath')

def lift_vec(c,d,s):
    a2,a3,a6,a7=[float(x.subs({C:c,D:d,sig:s})) for x in center_lift(C,D,sig)]
    return [1,a2,a3,c,d,a6,a7,1]

random.seed(7)
roots=[]; tries=0; t0=time.time()
# seed grid: vary sigma widely, C in (0.5,0.98), D in (0.05,0.45)
sigmas=[0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.3,1.7,2.2,3.0]
while tries<400 and time.time()-t0<240:
    tries+=1
    s0=random.choice(sigmas)*random.uniform(0.7,1.3)
    c0=random.uniform(0.5,0.98); d0=random.uniform(0.05,0.45)
    if c0+d0<=1.02 or c0<=d0: continue
    try:
        sol=mp.findroot(lambda c,d,s:(fG(c,d,s),fS(c,d,s),fT(c,d,s)),
                        (mp.mpf(c0),mp.mpf(d0),mp.mpf(s0)),tol=1e-20,maxsteps=40)
        c1,d1,s1=float(sol[0]),float(sol[1]),float(sol[2])
        if not (0<d1<c1<1 and c1+d1>1.0 and s1>0): continue
        aa=lift_vec(c1,d1,s1)
        if not (all(0<ai<1 for ai in aa) and all(aa[i]+aa[i+1]>1 for i in range(7))): continue
        res=max(abs(float(fG(c1,d1,s1))),abs(float(fS(c1,d1,s1))),abs(float(fT(c1,d1,s1))))
        if res>1e-6: continue
        # P<9 check: rho^9=N/D, rho=(N/D)^(1/9), P=rho*(1+rho)*sigma*sum(a)
        rho9=float(fN(c1,d1,s1))/float(fDd(c1,d1,s1))
        if rho9<=0: roots.append((c1,d1,s1,res,None,"rho9<=0")); continue
        rho=rho9**(1/9)
        Pval=rho*(1+rho)*s1*sum(aa)
        roots.append((c1,d1,s1,res,Pval,"P<9" if Pval<9 else "P>=9"))
    except Exception:
        pass
print(f"{tries} starts -> {len(roots)} strict-interior {{G=S=Theta=0}} roots",flush=True)
for r in roots[:20]:
    print(f"  C={r[0]:.6f} D={r[1]:.6f} sig={r[2]:.6f} res={r[3]:.1e} P={r[4]} {r[5]}",flush=True)
print("DONE",flush=True)
