#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B.16 empty-set certificate via interval B&B in (C,D,sigma) with factored Theta.

Variables: 0<D<C<1, C+D>1, sigma>0.  sigma bounded above by (6.1).
Equations: G_cd=0, S_cd=0 (polynomials in C,D,sigma), Theta=0 (factored rational,
  all denominators positive in strict interior -> clean interval division).

Pruning per box (exclude if provably empty):
  - domain: D>=C or D<=0 or C>=1 or C+D<=1
  - sigma bound (6.1): sigma >= min{C(1-C)/gap, D(1-D)/gap} everywhere
  - lift positivity: a3,a6,a2,a7 in (0,1); h_i=a_i+a_{i+1}-1>0; gap>0; a3+C-1>0; a6+D-1>0
  - G_cd, S_cd, Theta intervals each must contain 0
Subdivide largest axis. Survivors below MINW -> mpmath.iv re-check.
If all boxes excluded -> {G=S=Theta=0} strict interior certified EMPTY.
"""
import re, time
from pathlib import Path
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
X, Y, sigma = sp.symbols('X Y sigma')
C, D, sig = sp.symbols('C D sigma')

def load_small(name):
    s = sp.symbols("s")
    text = (HERE/name).read_text(encoding="utf-8").strip()
    return sp.Poly(sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma),
                   X, Y, sigma, domain=sp.ZZ)
G = load_small("nonpal_G_clean.txt"); S = load_small("nonpal_S_clean.txt")
Gcd = sp.Poly(sp.expand(G.as_expr().subs({X: C+D, Y: C*D})), C, D, sigma, domain=sp.ZZ)
Scd = sp.Poly(sp.expand(S.as_expr().subs({X: C+D, Y: C*D})), C, D, sigma, domain=sp.ZZ)
def tl(P): return [(tuple(map(int,k)),int(v)) for k,v in P.terms()]
Gterms = tl(Gcd); Sterms = tl(Scd)
print(f"Gcd deg={Gcd.total_degree()} terms={len(Gterms)}; Scd deg={Scd.total_degree()} terms={len(Sterms)}", flush=True)

# ---- interval ops (vectorized over numpy arrays of box bounds) ----
REL=1e-10
def iv_pow(lo,hi,n):  # positive vars
    return lo**n, hi**n
def iv_mul(alo,ahi,blo,bhi):
    t1=alo*blo; t2=alo*bhi; t3=ahi*blo; t4=ahi*bhi
    return np.minimum(np.minimum(t1,t2),np.minimum(t3,t4)), np.maximum(np.maximum(t1,t2),np.maximum(t3,t4))
def iv_div_pos(alo,ahi,blo,bhi):  # blo>0
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(bhi>0, alo/bhi, -np.inf), np.where(blo>0, ahi/blo, np.inf)
def iv_add(alo,ahi,blo,bhi): return alo+blo, ahi+bhi
def iv_sub(alo,ahi,blo,bhi): return alo-bhi, ahi-blo
def iv_recip_pos(blo,bhi):  # 1/[blo,bhi], blo>0
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(bhi>0,1.0/bhi,-np.inf), np.where(blo>0,1.0/blo,np.inf)
def inflate(lo,hi):
    w=hi-lo; return lo-w*REL-1e-12, hi+w*REL+1e-12

# poly interval eval (positive vars C,D,sigma)
def iv_poly(terms, Clo,Chi,Dlo,Dhi,Slo,Shi):
    lo=np.zeros_like(Clo); hi=np.zeros_like(Clo)
    for (dc,dd,ds),c in terms:
        clo=1.0 if dc==0 else Clo**dc; chi=1.0 if dc==0 else Chi**dc
        dlo=1.0 if dd==0 else Dlo**dd; dhi=1.0 if dd==0 else Dhi**dd
        slo=1.0 if ds==0 else Slo**ds; shi=1.0 if ds==0 else Shi**ds
        mlo=clo*dlo*slo; mhi=chi*dhi*shi
        if c>=0: lo=lo+c*mlo; hi=hi+c*mhi
        else: lo=lo+c*mhi; hi=hi+c*mlo
    return inflate(lo,hi)

# ---- Theta factored interval eval ----
# a3 = 1-C + sigma*C^2*(1-C)/gap,  gap=C+D-1
# a6 = 1-D + sigma*D^2*(1-D)/gap
# a2 = 1-a3 + sigma*a3^2*(1-a3)/(a3+C-1)
# a7 = 1-a6 + sigma*a6^2*(1-a6)/(a6+D-1)
def iv_theta(Clo,Chi,Dlo,Dhi,Slo,Shi):
    One_lo=np.ones_like(Clo); One_hi=np.ones_like(Clo)
    # C, D, sigma intervals
    # gap = C+D-1 >0
    glo,ghi = iv_sub(*iv_add(Clo,Chi,Dlo,Dhi), One_lo,One_hi)
    # C^2, (1-C), D^2, (1-D)
    C2lo,C2hi=iv_pow(Clo,Chi,2); D2lo,D2hi=iv_pow(Dlo,Dhi,2)
    C1m_lo,C1m_hi=iv_sub(One_lo,One_hi,Clo,Chi)   # 1-C
    D1m_lo,D1m_hi=iv_sub(One_lo,One_hi,Dlo,Dhi)   # 1-D
    # sigma*C^2*(1-C)
    sC2lo,sC2hi=iv_mul(Slo,Shi,C2lo,C2hi)
    sC2C1lo,sC2C1hi=iv_mul(sC2lo,sC2hi,C1m_lo,C1m_hi)
    num3lo,num3hi=iv_div_pos(sC2C1lo,sC2C1hi,glo,ghi)  # /gap
    a3lo,a3hi=iv_add(*iv_sub(One_lo,One_hi,Clo,Chi), num3lo,num3hi)  # 1-C + ...
    # a6
    sD2lo,sD2hi=iv_mul(Slo,Shi,D2lo,D2hi)
    sD2D1lo,sD2D1hi=iv_mul(sD2lo,sD2hi,D1m_lo,D1m_hi)
    num6lo,num6hi=iv_div_pos(sD2D1lo,sD2D1hi,glo,ghi)
    a6lo,a6hi=iv_add(*iv_sub(One_lo,One_hi,Dlo,Dhi), num6lo,num6hi)
    # a2 = 1-a3 + sigma*a3^2*(1-a3)/(a3+C-1)
    a3Cm1lo,a3Cm1hi=iv_sub(*iv_add(a3lo,a3hi,Clo,Chi), One_lo,One_hi)  # a3+C-1 >0
    a32lo,a32hi=iv_pow(a3lo,a3hi,2)
    a31mlo,a31mhi=iv_sub(One_lo,One_hi,a3lo,a3hi)  # 1-a3
    sa32lo,sa32hi=iv_mul(Slo,Shi,a32lo,a32hi)
    sa32a1lo,sa32a1hi=iv_mul(sa32lo,sa32hi,a31mlo,a31mhi)
    q2lo,q2hi=iv_div_pos(sa32a1lo,sa32a1hi,a3Cm1lo,a3Cm1hi)
    a2lo,a2hi=iv_add(*iv_sub(One_lo,One_hi,a3lo,a3hi), q2lo,q2hi)
    # a7 = 1-a6 + sigma*a6^2*(1-a6)/(a6+D-1)
    a6Dm1lo,a6Dm1hi=iv_sub(*iv_add(a6lo,a6hi,Dlo,Dhi), One_lo,One_hi)
    a62lo,a62hi=iv_pow(a6lo,a6hi,2)
    a61mlo,a61mhi=iv_sub(One_lo,One_hi,a6lo,a6hi)
    sa62lo,sa62hi=iv_mul(Slo,Shi,a62lo,a62hi)
    sa62a1lo,sa62a1hi=iv_mul(sa62lo,sa62hi,a61mlo,a61mhi)
    q7lo,q7hi=iv_div_pos(sa62a1lo,sa62a1hi,a6Dm1lo,a6Dm1hi)
    a7lo,a7hi=iv_add(*iv_sub(One_lo,One_hi,a6lo,a6hi), q7lo,q7hi)
    # a4=C, a5=D, a1=a8=1
    a4lo,a4hi=Clo,Chi; a5lo,a5hi=Dlo,Dhi
    a1lo,a1hi=One_lo,One_hi
    # h_i = a_i + a_{i+1} - 1
    h2lo,h2hi=iv_sub(*iv_add(a2lo,a2hi,a3lo,a3hi),One_lo,One_hi)
    h3lo,h3hi=iv_sub(*iv_add(a3lo,a3hi,a4lo,a4hi),One_lo,One_hi)
    h4lo,h4hi=iv_sub(*iv_add(a4lo,a4hi,a5lo,a5hi),One_lo,One_hi)
    h5lo,h5hi=iv_sub(*iv_add(a5lo,a5hi,a6lo,a6hi),One_lo,One_hi)
    h6lo,h6hi=iv_sub(*iv_add(a6lo,a6hi,a7lo,a7hi),One_lo,One_hi)
    # W23 = (a6*a7/h6) * B1 * B2
    a6a7lo,a6a7hi=iv_mul(a6lo,a6hi,a7lo,a7hi)
    a6a7h6lo,a6a7hi_=iv_div_pos(a6a7lo,a6a7hi,h6lo,h6hi)
    # B1 = sigma + h2*h6/(sigma*a2*a3*a6*a7) - h2*h3*h4*h5*h6/(sigma^4*a2*a3^2*a4^2*a5^2*a6^2*a7)
    h2h6lo,h2h6hi=iv_mul(h2lo,h2hi,h6lo,h6hi)
    s_a2a3a6a7lo,s_a2a3a6a7hi=iv_mul(*iv_mul(Slo,Shi,a2lo,a2hi),*iv_mul(a3lo,a3hi,a6a7lo,a6a7hi))
    t1lo,t1hi=iv_div_pos(h2h6lo,h2h6hi,s_a2a3a6a7lo,s_a2a3a6a7hi)
    h2h3h4h5h6lo,h2h3h4h5h6hi=iv_mul(*iv_mul(h2lo,h2hi,h3lo,h3hi),*iv_mul(*iv_mul(h4lo,h4hi,h5lo,h5hi),h6lo,h6hi))
    a3sqlo,a3sqhi=iv_pow(a3lo,a3hi,2); a4sqlo,a4sqhi=iv_pow(a4lo,a4hi,2); a5sqlo,a5sqhi=iv_pow(a5lo,a5hi,2); a6sqlo,a6sqhi=iv_pow(a6lo,a6hi,2)
    s4lo,s4hi=iv_pow(Slo,Shi,4)
    den_B1_2lo,den_B1_2hi=iv_mul(*iv_mul(*iv_mul(s4lo,s4hi,a2lo,a2hi),*iv_mul(a3sqlo,a3sqhi,a4sqlo,a4sqhi)),*iv_mul(*iv_mul(a5sqlo,a5sqhi,a6sqlo,a6sqhi),a7lo,a7hi))
    t2lo,t2hi=iv_div_pos(h2h3h4h5h6lo,h2h3h4h5h6hi,den_B1_2lo,den_B1_2hi)
    B1lo,B1hi=iv_sub(*iv_sub(*iv_add(Slo,Shi,t1lo,t1hi),t2lo,t2hi), np.zeros_like(Clo),np.zeros_like(Clo))  # sigma+t1-t2 (sub 0 for the -t2 already done)
    B1lo,B1hi=iv_sub(*iv_add(Slo,Shi,t1lo,t1hi), t2lo,t2hi)
    # B2 = 1 + a2*a3*h6/(h2*a6*a7)
    a2a3h6lo,a2a3h6hi=iv_mul(*iv_mul(a2lo,a2hi,a3lo,a3hi),h6lo,h6hi)
    h2a6a7lo,h2a6a7hi=iv_mul(h2lo,h2hi,a6a7lo,a6a7hi)
    t3lo,t3hi=iv_div_pos(a2a3h6lo,a2a3h6hi,h2a6a7lo,h2a6a7hi)
    B2lo,B2hi=iv_add(One_lo,One_hi,t3lo,t3hi)
    W23lo,W23hi=iv_mul(*iv_mul(a6a7h6lo,a6a7hi_,B1lo,B1hi),B2lo,B2hi)
    # W24 = (h2*h3/(sigma^2*a2*a3^2*a4)) * B3 * B4
    h2h3lo,h2h3hi=iv_mul(h2lo,h2hi,h3lo,h3hi)
    s2lo,s2hi=iv_pow(Slo,Shi,2)
    den24lo,den24hi=iv_mul(*iv_mul(*iv_mul(s2lo,s2hi,a2lo,a2hi),*iv_mul(a3sqlo,a3sqhi,a4lo,a4hi)),One_lo,One_hi)
    fac24lo,fac24hi=iv_div_pos(h2h3lo,h2h3hi,den24lo,den24hi)
    # B3 = h4/(sigma*a4*a5) - 1
    sa4a5lo,sa4a5hi=iv_mul(*iv_mul(Slo,Shi,a4lo,a4hi),a5lo,a5hi)
    t4lo,t4hi=iv_div_pos(h4lo,h4hi,sa4a5lo,sa4a5hi)
    B3lo,B3hi=iv_sub(t4lo,t4hi,One_lo,One_hi)
    # B4 = 1 + a2*a3^2*a4*h5*h6/(h2*h3*a5*a6^2*a7)
    numB4lo,numB4hi=iv_mul(*iv_mul(*iv_mul(a2lo,a2hi,a3sqlo,a3sqhi),*iv_mul(a4lo,a4hi,h5lo,h5hi)),h6lo,h6hi)
    denB4lo,denB4hi=iv_mul(*iv_mul(*iv_mul(h2lo,h2hi,h3lo,h3hi),*iv_mul(a5lo,a5hi,a6sqlo,a6sqhi)),a7lo,a7hi)
    t5lo,t5hi=iv_div_pos(numB4lo,numB4hi,denB4lo,denB4hi)
    B4lo,B4hi=iv_add(One_lo,One_hi,t5lo,t5hi)
    W24lo,W24hi=iv_mul(*iv_mul(fac24lo,fac24hi,B3lo,B3hi),B4lo,B4hi)
    # W34 = (a5*a6/h5) * B5 * B6
    a5a6lo,a5a6hi=iv_mul(a5lo,a5hi,a6lo,a6hi)
    a5a6h5lo,a5a6h5hi=iv_div_pos(a5a6lo,a5a6hi,h5lo,h5hi)
    # B5 = sigma + h3*h5/(sigma*a3*a4*a5*a6) - h3*h4*h5/(sigma^2*a3*a4^2*a5^2*a6)
    h3h5lo,h3h5hi=iv_mul(h3lo,h3hi,h5lo,h5hi)
    den5_1lo,den5_1hi=iv_mul(*iv_mul(*iv_mul(Slo,Shi,a3lo,a3hi),*iv_mul(a4lo,a4hi,a5lo,a5hi)),a6lo,a6hi)
    u1lo,u1hi=iv_div_pos(h3h5lo,h3h5hi,den5_1lo,den5_1hi)
    h3h4h5lo,h3h4h5hi=iv_mul(*iv_mul(h3lo,h3hi,h4lo,h4hi),h5lo,h5hi)
    den5_2lo,den5_2hi=iv_mul(*iv_mul(*iv_mul(s2lo,s2hi,a3lo,a3hi),*iv_mul(a4sqlo,a4sqhi,a5sqlo,a5sqhi)),a6lo,a6hi)
    u2lo,u2hi=iv_div_pos(h3h4h5lo,h3h4h5hi,den5_2lo,den5_2hi)
    B5lo,B5hi=iv_sub(*iv_add(Slo,Shi,u1lo,u1hi),u2lo,u2hi)
    # B6 = 1 + a3*a4*h5/(h3*a5*a6)
    numB6lo,numB6hi=iv_mul(*iv_mul(a3lo,a3hi,a4lo,a4hi),h5lo,h5hi)
    denB6lo,denB6hi=iv_mul(*iv_mul(h3lo,h3hi,a5lo,a5hi),a6lo,a6hi)
    u3lo,u3hi=iv_div_pos(numB6lo,numB6hi,denB6lo,denB6hi)
    B6lo,B6hi=iv_add(One_lo,One_hi,u3lo,u3hi)
    W34lo,W34hi=iv_mul(*iv_mul(a5a6h5lo,a5a6h5hi,B5lo,B5hi),B6lo,B6hi)
    # Theta = W23*W24 + W23*W34 + W24*W34
    p1lo,p1hi=iv_mul(W23lo,W23hi,W24lo,W24hi)
    p2lo,p2hi=iv_mul(W23lo,W23hi,W34lo,W34hi)
    p3lo,p3hi=iv_mul(W24lo,W24hi,W34lo,W34hi)
    Tlo,Thi=iv_add(*iv_add(p1lo,p1hi,p2lo,p2hi),p3lo,p3hi)
    Tlo,Thi=inflate(Tlo,Thi)
    return Tlo,Thi,(a2lo,a2hi,a3lo,a3hi,a6lo,a6hi,a7lo,a7hi,
                              h2lo,h2hi,h3lo,h3hi,h4lo,h4hi,h5lo,h5hi,h6lo,h6hi,
                              glo,ghi,a3Cm1lo,a3Cm1hi,a6Dm1lo,a6Dm1hi)

def cheap_excl(Clo,Chi,Dlo,Dhi,Slo,Shi):
    excl=np.zeros(len(Clo),dtype=bool)
    One=np.ones_like(Clo)
    # domain: D<C, 0<D, C<1, C+D>1  (exclude only when provably violated everywhere)
    excl |= (Dlo>=Chi)|(Dhi<=0)|(Clo>=1)|(Chi+Dhi<=1)
    if excl.all(): return excl
    # gap=C+D-1>0 (already implied by C+D>1)
    # sigma bound (6.1): sigma < min{C(1-C)/gap, D(1-D)/gap}, gap=C+D-1>0
    # ONLY valid where gap>0 everywhere (glo>0); else subdivide (don't apply)
    glo=Clo+Dlo-1; ghi=Chi+Dhi-1
    valid = glo>0  # gap definitely >0
    if np.any(valid):
        with np.errstate(divide='ignore',invalid='ignore'):
          C1m_lo=np.where(valid,1-Chi,0); C1m_hi=np.where(valid,1-Clo,0)
          D1m_lo=np.where(valid,1-Dhi,0); D1m_hi=np.where(valid,1-Dlo,0)
          C2lo=np.where(valid,Clo**2,0); C2hi=np.where(valid,Chi**2,0)
          D2lo=np.where(valid,Dlo**2,0); D2hi=np.where(valid,Dhi**2,0)
          # gap>0 so C(1-C)/gap: lower=(1-Chi)*Clo^2/ghi, upper=(1-Clo)*Chi^2/glo
          bndC_lo=np.where(valid,C1m_lo*C2lo/ghi, np.inf)
          bndC_hi=np.where(valid,C1m_hi*C2hi/glo, -np.inf)
          bndD_lo=np.where(valid,D1m_lo*D2lo/ghi, np.inf)
          bndD_hi=np.where(valid,D1m_hi*D2hi/glo, -np.inf)
          bnd_hi=np.minimum(bndC_hi,bndD_hi)  # upper bound on sigma_max
          excl |= valid & (Slo>=bnd_hi)
    if excl.all(): return excl
    return excl

def prune_batch(Clo,Chi,Dlo,Dhi,Slo,Shi):
    excl=cheap_excl(Clo,Chi,Dlo,Dhi,Slo,Shi)
    idx=np.where(~excl)[0]
    if len(idx)==0: return excl
    def sub(a,i=idx): return a[i]
    # G, S contain 0
    g_lo,g_hi=iv_poly(Gterms,sub(Clo),sub(Chi),sub(Dlo),sub(Dhi),sub(Slo),sub(Shi))
    excl[idx] |= (g_lo>0)|(g_hi<0)
    idx=idx[(g_lo<=0)&(g_hi>=0)]
    if len(idx)==0: return excl
    def sub2(a,i=idx): return a[i]
    s_lo,s_hi=iv_poly(Sterms,sub2(Clo),sub2(Chi),sub2(Dlo),sub2(Dhi),sub2(Slo),sub2(Shi))
    excl[idx] |= (s_lo>0)|(s_hi<0)
    idx=idx[(s_lo<=0)&(s_hi>=0)]
    if len(idx)==0: return excl
    # Theta + lift positivity on remaining
    def sub3(a,i=idx): return a[i]
    T_lo,T_hi,aux=iv_theta(sub3(Clo),sub3(Chi),sub3(Dlo),sub3(Dhi),sub3(Slo),sub3(Shi))
    excl[idx] |= (T_lo>0)|(T_hi<0)
    # lift positivity: a_i in (0,1), h_i>0, a3+C-1>0, a6+D-1>0
    # exclude only when PROVABLY violated (definitely outside); straddling -> subdivide
    (a2lo,a2hi,a3lo,a3hi,a6lo,a6hi,a7lo,a7hi,h2lo,h2hi,h3lo,h3hi,h4lo,h4hi,h5lo,h5hi,h6lo,h6hi,glo,ghi,a3Cm1lo,a3Cm1hi,a6Dm1lo,a6Dm1hi)=aux
    lex = (a2hi<=0)|(a2lo>=1)|(a3hi<=0)|(a3lo>=1)|(a6hi<=0)|(a6lo>=1)|(a7hi<=0)|(a7lo>=1)
    lex |= (h2hi<=0)|(h3hi<=0)|(h4hi<=0)|(h5hi<=0)|(h6hi<=0)
    lex |= (a3Cm1hi<=0)|(a6Dm1hi<=0)
    excl[idx] |= lex
    return excl

# ---- B&B driver ----
def run():
  # initial box: C in (0.5,1), D in (0,0.5), sigma in (0, SMAX)
  # refine: C+D>1 so C>0.5 when D small. Use C(0.5,1) D(0,0.5) sigma(0,4)
  SMAX=4.0
  queue=[(0.5,1.0, 0.0,0.5, 0.0,SMAX)]
  survivors=[]; npruned=0; nproc=0; t0=time.time()
  MAXTIME=600; MINW=1e-4
  print(f"start B&B (C,D,sigma), SMAX={SMAX} minw={MINW}", flush=True)
  while queue and time.time()-t0<MAXTIME:
    batch=queue[:8192]; queue=queue[8192:]
    Clo=np.array([b[0] for b in batch]);Chi=np.array([b[1] for b in batch])
    Dlo=np.array([b[2] for b in batch]);Dhi=np.array([b[3] for b in batch])
    Slo=np.array([b[4] for b in batch]);Shi=np.array([b[5] for b in batch])
    excl=prune_batch(Clo,Chi,Dlo,Dhi,Slo,Shi)
    nproc+=len(batch)
    for i,b in enumerate(batch):
        if excl[i]: npruned+=1; continue
        w=max(b[1]-b[0],b[3]-b[2],b[5]-b[4])
        if w<2*MINW: survivors.append(b)
        else:
            ax=b[1]-b[0];ay=b[3]-b[2];az=b[5]-b[4]
            if ax>=ay and ax>=az:
                m=(b[0]+b[1])/2; queue.append((b[0],m,b[2],b[3],b[4],b[5]));queue.append((m,b[1],b[2],b[3],b[4],b[5]))
            elif ay>=az:
                m=(b[2]+b[3])/2; queue.append((b[0],b[1],b[2],m,b[4],b[5]));queue.append((b[0],b[1],m,b[3],b[4],b[5]))
            else:
                m=(b[4]+b[5])/2; queue.append((b[0],b[1],b[2],b[3],b[4],m));queue.append((b[0],b[1],b[2],b[3],m,b[5]))
    if nproc%100000<8192:
        print(f"  proc={nproc} pruned={npruned} surv={len(survivors)} Q={len(queue)} t={time.time()-t0:.0f}s", flush=True)
  print(f"\nDONE proc={nproc} pruned={npruned} surv={len(survivors)} Q={len(queue)} t={time.time()-t0:.0f}s", flush=True)
  if not survivors and not queue:
    print("==> ALL BOXES EXCLUDED: {G=S=Theta=0} strict interior certified EMPTY.", flush=True)
  elif survivors:
    print(f"==> {len(survivors)} survivors below minw {MINW}:", flush=True)
    for b in survivors[:25]:
        print(f"   C[{b[0]:.6f},{b[1]:.6f}] D[{b[2]:.6f},{b[3]:.6f}] sig[{b[4]:.6f},{b[5]:.6f}]", flush=True)

if __name__=='__main__':
    run()

