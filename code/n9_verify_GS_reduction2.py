#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Numerically verify G = primitive((E_L-E_R)/(C-D)) up to (X-1)^k, and identify S.
Evaluate the FULL rational (E_L-E_R)/(C-D) at random (C,D,sigma), compare to G(X,Y,sigma)."""
import re, sympy as sp, random, math
from pathlib import Path
HERE=Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'
X,Y,s=sp.symbols('X Y s')
def parse_poly(name):
    txt=(HERE/name).read_text().strip()
    terms=re.findall(r'[+-][^+-]+|^[^+-]+', txt)
    monos={}
    for t in terms:
        t=t.strip()
        if not t: continue
        sign=1
        if t.startswith('-'): sign=-1; t=t[1:].strip()
        elif t.startswith('+'): t=t[1:].strip()
        xp=yp=sp_=0
        for fac in re.finditer(r'(X|Y|s)(?:\*\*(\d+))?', t):
            base=fac.group(1); exp=int(fac.group(2)) if fac.group(2) else 1
            if base=='X': xp=exp
            elif base=='Y': yp=exp
            else: sp_=exp
        lead=re.match(r'(\d+)', t)
        coeff=sign*int(lead.group(1)) if lead else sign
        key=(xp,yp,sp_); monos[key]=monos.get(key,0)+coeff
    return sp.Poly({k:v for k,v in monos.items()},X,Y,s,domain=sp.ZZ).as_expr()
Gexpr=parse_poly('nonpal_G_clean.txt')
Sexpr=parse_poly('nonpal_S_clean.txt')

C,D,sigma=sp.symbols('C D sigma')
def center_lift(C,D,sigma):
    den=C+D-1
    a3=1-C+sigma*C**2*(1-C)/den
    a6=1-D+sigma*D**2*(1-D)/den
    a2=1-a3+sigma*a3**2*(1-a3)/(a3+C-1)
    a7=1-a6+sigma*a6**2*(1-a6)/(a6+D-1)
    return a2,a3,a6,a7
a2,a3,a6,a7=center_lift(C,D,sigma)
E_L = a2+a3-1 - sigma*a2*(1-a2)
E_R = a6+a7-1 - sigma*a7*(1-a7)
diff = sp.together(E_L - E_R)         # full rational
summ = sp.together(E_L + E_R)
EL   = sp.together(E_L)

random.seed(7)
print("Testing (E_L-E_R)/(C-D) vs G*(X-1)^k:", flush=True)
ratios=[]
for _ in range(10):
    while True:
        Cv=random.uniform(0.35,0.95); Dv=random.uniform(0.35,0.95)
        if abs(Cv-Dv)>0.05 and Cv+Dv>1.05 and Cv+Dv<1.95: break
    sv=random.uniform(0.2,2.0)
    subs={C:Cv,D:Dv,sigma:sv}
    Xv=Cv+Dv; Yv=Cv*Dv
    diffv=float(diff.subs(subs));   # = (E_L-E_R) full rational
    CDv=Cv-Dv
    quotv=diffv/CDv                  # (E_L-E_R)/(C-D)
    gv=float(Gexpr.subs({X:Xv,Y:Yv,s:sv}))
    if abs(gv)<1e-6 or abs(quotv)<1e-9: continue
    r=quotv/gv
    ratios.append((Xv,Yv,sv,r,Xv-1))
    print(f"  C={Cv:.3f} D={Dv:.3f} sig={sv:.3f} X-1={Xv-1:.4f}: (E_L-E_R)/(C-D)/G = {r:.6g}", flush=True)

# determine k from ratios: r_i = const * (X_i-1)^k  => r_i/r_j = ((X_i-1)/(X_j-1))^k
if len(ratios)>=2:
    const_k=None
    ok=True
    for i in range(1,len(ratios)):
        X1,_,_,r1,x1=ratios[0]; X2,_,_,r2,x2=ratios[i]
        if x1<=0 or x2<=0 or r1<=0 or r2<=0: ok=False; break
        k=math.log(r1/r2)/math.log(x1/x2)
        const_k = k if const_k is None else const_k
        print(f"  k from pair(0,{i}) = {k:.5f}", flush=True)
    # check constant of proportionality = r/(X-1)^k
    if ratios:
        X1,_,_,r1,x1=ratios[0]
        k_round=round(const_k) if const_k else 0
        const=r1/(x1**k_round)
        print(f"  -> k={k_round}, const={const:.6g}; checking all: ", end="", flush=True)
        good=True
        for Xv,Yv,sv,r,xm1 in ratios:
            pred=const*(xm1**k_round)
            if abs(pred-r)>1e-6*max(1,abs(r)): good=False
        print("ALL MATCH" if good else "MISMATCH", flush=True)
        print(f"  VERDICT: (E_L-E_R)/(C-D) = const * (X-1)^{k_round} * G  ==> G reduction {'VERIFIED' if good else 'FAILED'}", flush=True)

# Now identify S: test (E_L+E_R) and E_L against S*(X-1)^m
print("\nTesting (E_L+E_R) and E_L vs S:", flush=True)
for label, expr in [('E_L+E_R', summ), ('E_L', EL)]:
    rs=[]
    for _ in range(8):
        while True:
            Cv=random.uniform(0.35,0.95); Dv=random.uniform(0.35,0.95)
            if abs(Cv-Dv)>0.05 and Cv+Dv>1.05 and Cv+Dv<1.95: break
        sv=random.uniform(0.2,2.0)
        subs={C:Cv,D:Dv,sigma:sv}; Xv=Cv+Dv; Yv=Cv*Dv
        ev=float(expr.subs(subs)); sv_s=float(Sexpr.subs({X:Xv,Y:Yv,s:sv}))
        if abs(sv_s)<1e-6 or abs(ev)<1e-9: continue
        rs.append((Xv,Yv,sv,ev/sv_s,Xv-1))
    if len(rs)>=2:
        # find m
        X1,_,_,r1,x1=rs[0]; X2,_,_,r2,x2=rs[1]
        if x1>0 and x2>0 and r1>0 and r2>0:
            m=round(math.log(r1/r2)/math.log(x1/x2))
            const=r1/(x1**m); good=all(abs(const*(xm**m)-r)<1e-6*max(1,abs(r)) for _,_,_,r,xm in rs)
            print(f"  {label}/S = const*(X-1)^{m}: {'VERIFIED' if good else 'no'}", flush=True)
