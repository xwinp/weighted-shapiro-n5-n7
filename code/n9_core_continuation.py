#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core tube scaffolding (GPT#12): dense numerical continuation of the
non-palindromic stationary curve {G^-=G^+=0} in the CORE v>=3/20=0.15.
Track each branch over tau in (0.15-range .. 0.97) by continuation (previous
solution as next seed).  Record (tau, v, xi, Theta, J_tau) per branch.
Output branch structure: count, tau-range, min|Theta|, min J_tau.
mpmath only -- bounded memory, safe.
"""
import random, mpmath as mp, sympy as sp, math, json
mp.mp.dps=22
c,d,s=sp.symbols('c d sigma'); v,xi,tau=sp.symbols('v xi tau')
vv=c+d-1
Lc=vv+s*c**2; Uc=vv-s*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+s*d**2; Ud=vv-s*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - s*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - s*Bd*(d*vv**2-Ud*Bd)
Pc=sp.Poly(sp.expand(Fc),c,d,s,domain=sp.ZZ); Pd=sp.Poly(sp.expand(Fd),c,d,s,domain=sp.ZZ)
FL=Pc.exquo(sp.Poly((c-1)*Lc,c,d,s,domain=sp.ZZ)); FR=Pd.exquo(sp.Poly((d-1)*Ld,c,d,s,domain=sp.ZZ))
Gm=sp.Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,s,domain=sp.ZZ).exquo(sp.Poly(c-d,c,d,s,domain=sp.ZZ))
Gp=sp.Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,s,domain=sp.ZZ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,s:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,s:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])
fGm=sp.lambdify((v,xi,tau),Gmt,'mpmath'); fGp=sp.lambdify((v,xi,tau),Gpt,'mpmath')
# partials for J_tau = d(v,xi)/d(G-,G+) = (Gm_v*Gp_xi - Gm_xi*Gp_v)
fGm_v=sp.lambdify((v,xi,tau),sp.diff(Gmt,v),'mpmath'); fGm_xi=sp.lambdify((v,xi,tau),sp.diff(Gmt,xi),'mpmath')
fGp_v=sp.lambdify((v,xi,tau),sp.diff(Gpt,v),'mpmath'); fGp_xi=sp.lambdify((v,xi,tau),sp.diff(Gpt,xi),'mpmath')
def Jtau(vv,xx,tt):
    return fGm_v(vv,xx,tt)*fGp_xi(vv,xx,tt)-fGm_xi(vv,xx,tt)*fGp_v(vv,xx,tt)
def Th_direct(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    a4=c1; a5=d1
    h2=a2+a3-1; h3=a3+a4-1; h4=a4+a5-1; h5=a5+a6-1; h6=a6+a7-1
    W23=((a6*a7/h6)*(s1+h2*h6/(s1*a2*a3*a6*a7)-h2*h3*h4*h5*h6/(s1**4*a2*a3**2*a4**2*a5**2*a6**2*a7))*(1+a2*a3*h6/(h2*a6*a7)))
    W24=(h2*h3/(s1**2*a2*a3**2*a4)*(h4/(s1*a4*a5)-1)*(1+a2*a3**2*a4*h5*h6/(h2*h3*a5*a6**2*a7)))
    W34=((a5*a6/h5)*(s1+h3*h5/(s1*a3*a4*a5*a6)-h3*h4*h5/(s1**2*a3*a4**2*a5**2*a6))*(1+a3*a4*h5/(h3*a5*a6)))
    return W23*W24+W23*W34+W24*W34
def fullvec(vv,xx,tt):
    c1=(1+vv+(1-vv)*xx)/2; d1=(1+vv-(1-vv)*xx)/2; s1=tt*vv/(d1*(1-d1)); g=c1+d1-1
    a3=1-c1+s1*c1**2*(1-c1)/g; a6=1-d1+s1*d1**2*(1-d1)/g
    a2=1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7=1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]
def admiss(vv,xx,tt):
    if not(0<vv<1 and 0<xx<1): return False
    a=fullvec(vv,xx,tt)
    if not all(0<float(a[i])<1 for i in range(1,7)): return False
    if not all(float(a[i]+a[i+1])>1 for i in range(7)): return False
    if abs(float(cc.subs({v:mp.mpf(vv),xi:mp.mpf(xx)})-dd.subs({v:mp.mpf(vv),xi:mp.mpf(xx)})))<1e-3: return False  # non-pal
    return True
def solve(v0,x0,t0):
    try:
        sol=mp.findroot(lambda vv,xx:(fGm(vv,xx,mp.mpf(t0)),fGp(vv,xx,mp.mpf(t0))),(mp.mpf(v0),mp.mpf(x0)),tol=1e-22,maxsteps=60)
        v1,x1=float(sol[0]),float(sol[1])
        r=max(abs(float(fGm(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))),abs(float(fGp(mp.mpf(v1),mp.mpf(x1),mp.mpf(t0)))))
        if r>1e-14: return None
        if not admiss(v1,x1,t0): return None
        if v1<0.15: return None  # core only
        return (v1,x1)
    except Exception:
        return None

# 1) find all core branch points at tau=0.5 by multi-start
random.seed(2024)
seeds=[]
for _ in range(400):
    v0=random.uniform(0.16,0.95); x0=random.uniform(0.02,0.98)
    r=solve(v0,x0,0.5)
    if r:
        # dedup
        if not any(abs(r[0]-s[0])<1e-6 and abs(r[1]-s[1])<1e-6 for s in seeds):
            seeds.append(r)
print("tau=0.5 core branch points:",len(seeds),[(round(s[0],4),round(s[1],4)) for s in seeds],flush=True)

# 2) continue each branch forward (tau up) and backward (tau down)
dtau=0.01
def continue_branch(seed, direction):
    pts=[(0.5,seed[0],seed[1])]
    t=0.5
    v_cur,x_cur=seed
    while True:
        t=t+direction*dtau
        if t<=0.15 or t>=0.98: break
        r=solve(v_cur,x_cur,t)
        if r is None:
            # try a tiny perturbation / smaller step
            r=solve(v_cur,x_cur,t)
            if r is None: break
        v_cur,x_cur=r
        pts.append((t,v_cur,x_cur))
    return pts

branches=[]
for sd in seeds:
    fwd=continue_branch(sd, +1)
    bwd=continue_branch(sd, -1)
    full=list(reversed(bwd))+fwd[1:]   # bwd came back from 0.5; reversed gives increasing tau
    branches.append(full)

# 3) compute Theta, J_tau along each branch; summarize
print("\n=== branch summary ===",flush=True)
allpts=[]
for i,br in enumerate(branches):
    recs=[]
    for (t,vv,xx) in br:
        th=float(Th_direct(mp.mpf(vv),mp.mpf(xx),mp.mpf(t)))
        jt=float(Jtau(mp.mpf(vv),mp.mpf(xx),mp.mpf(t)))
        recs.append((t,vv,xx,th,jt))
    if not recs: continue
    taus=[r[0] for r in recs]; vs=[r[1] for r in recs]
    ths=[r[3] for r in recs if r[3]==r[3]]; jts=[r[4] for r in recs if r[4]==r[4]
    allpts.append(recs)
    print("branch %d: n=%d tau=[%.3f,%.3f] v=[%.4f,%.4f] min|Theta|=%.4e Theta_max=%.4e min_Jtau=%.4e max_Jtau=%.4e"%(
        i,len(recs),min(taus),max(taus),min(vs),max(vs),
        min(abs(x) for x in ths),max(ths),min(jts),max(jts)),flush=True)
    # Theta sign check
    sign_changes=0
    for j in range(1,len(recs)):
        if recs[j-1][3]*recs[j][3]<0: sign_changes+=1
    print("   Theta sign changes: %d  (all Theta<0? %s)"%(sign_changes, all(r[3]<0 for r in recs if r[3]==r[3])),flush=True)

# save all points for tube building
with open('code/_core_branches.json','w') as f:
    json.dump([[(round(t,5),round(vv,6),round(xx,6),round(th,8),round(jt,8)) for (t,vv,xx,th,jt) in br] for br in allpts],f)
print("\nsaved code/_core_branches.json (%d branches, %d total pts)"%(len(allpts),sum(len(b) for b in allpts)),flush=True)
print("DONE",flush=True)
