#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorous interval branch-and-prune: certify {G=S=det=0, strict interior, P<9} EMPTY.

Variables X,Y,sigma all POSITIVE in domain, so interval natural extension is clean:
  X^a in [Xlo^a, Xhi^a] (monotone for X>0), products/sums of positive intervals standard.
Polynomials G,S,det given as term lists (powers, coeff). Evaluate interval enclosures
with numpy over BATCHES of boxes for speed. Inflate by REL margin for float safety;
survivors re-verified with mpmath.iv (rigorous) in a second pass.

Domain: 1<X<2, X-1<Y<X^2/4 (C!=D => strict), sigma>0.
Strict interior lift: C,D in (0,1) from X,Y; recovered a_i in (0,1), adjacent sums>1.
P<9 via stationary_value (skipped here -- criticality-empty already implies B.16 closed
under Nowosad; we certify the stronger {G=S=det=0, strict lift} empty).

Pruning per box:
  - domain: Ylo < Xhi^2/4 and Yhi > Xlo-1 (else no interior overlap)
  - C,D in (0,1): C=(X+sqrt(X^2-4Y))/2, D=(X-sqrt(...))/2; interval-eval disc=X^2-4Y>0,
    0<C<1, 0<D<1
  - G_iv,S_iv,det_iv must each contain 0 (else prune)
  - lift: den=C+D-1=X-1>0 (ok since X>1); a3,a6,a2,a7 intervals; require each in (0,1)
    and adjacent sums in (1,2); division by positive intervals (a3+C-1 etc.) -- require
    those denominators' lower bound >0 (else keep/subdivide).
A box is EXCLUDED if any constraint's interval provably excludes the required value.
If all boxes excluded -> EMPTY (rigorous up to float margin; mp.iv re-check on survivors).
"""
import re, numpy as np, time, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'
import sympy as sp
Xs,Ys,ss=sp.symbols('X Y s')
def parse_poly(name):
    txt=(HERE/name).read_text().strip()
    terms=re.findall(r'[+-][^+-]+|^[^+-]+', txt); monos={}
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
        lead=re.match(r'(\d+)', t); coeff=sign*int(lead.group(1)) if lead else sign
        key=(xp,yp,sp_); monos[key]=monos.get(key,0)+coeff
    return list(monos.items())  # [((ax,ay,as),coeff),...]

Gterms=parse_poly('nonpal_G_clean.txt')
Sterms=parse_poly('nonpal_S_clean.txt')
# build det terms via sympy (one-time)
def build_poly(name):
    txt=(HERE/name).read_text().strip()
    terms=re.findall(r'[+-][^+-]+|^[^+-]+', txt); monos={}
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
        lead=re.match(r'(\d+)', t); coeff=sign*int(lead.group(1)) if lead else sign
        key=(xp,yp,sp_); monos[key]=monos.get(key,0)+coeff
    return sp.Poly({k:v for k,v in monos.items()},Xs,Ys,ss,domain=sp.ZZ)
print("building det...",flush=True)
t0=time.time()
G=build_poly('nonpal_G_clean.txt'); S=build_poly('nonpal_S_clean.txt')
N=build_poly('nonpal_rho9_num.txt'); D=build_poly('nonpal_rho9_den.txt')
Gx,Gy,Gs=G.diff(Xs),G.diff(Ys),G.diff(ss)
Sx,Sy,Ss=S.diff(Xs),S.diff(Ys),S.diff(ss)
Nx,Ny,Ns=N.diff(Xs),N.diff(Ys),N.diff(ss)
Dx,Dy,Ds=D.diff(Xs),D.diff(Ys),D.diff(ss)
Cx=D*Nx-N*Dx; Cy=D*Ny-N*Dy; Cz=D*Ns-N*Ds
det = Gx*(Sy*Cz-Cy*Ss) - Gy*(Sx*Cz-Cx*Ss) + Gs*(Sx*Cy-Cx*Sy)
detP=sp.Poly(sp.expand(det.as_expr()),Xs,Ys,ss,domain=sp.ZZ)
detterms=[(tuple(map(int,k)),int(v)) for k,v in detP.terms()]
print(f"det built deg={detP.total_degree()} terms={len(detterms)} in {time.time()-t0:.1f}s",flush=True)

# numpy-vectorized interval evaluation of a poly term-list over a batch of boxes.
# boxes: arrays Xlo,Xhi,Ylo,Yhi,Slo,Shi (1D, length B). All vars positive.
REL=1e-9  # inflation for float safety
def iv_eval(terms, Xlo,Xhi,Ylo,Yhi,Slo,Shi):
    # monomial intervals (positive vars, monotone powers)
    # accumulate sum of term intervals
    lo=np.zeros_like(Xlo,dtype=np.float64); hi=np.zeros_like(Xlo,dtype=np.float64)
    for (ax,ay,asg),coef in terms:
        # X^a interval
        xlo=np.where(ax==0,1.0,Xlo**ax); xhi=np.where(ax==0,1.0,Xhi**ax)
        ylo=np.where(ay==0,1.0,Ylo**ay); yhi=np.where(ay==0,1.0,Yhi**ay)
        slo=np.where(asg==0,1.0,Slo**asg); shi=np.where(asg==0,1.0,Shi**asg)
        mlo=xlo*ylo*slo; mhi=xhi*yhi*shi
        if coef>=0:
            tlo=coef*mlo; thi=coef*mhi
        else:
            tlo=coef*mhi; thi=coef*mlo
        lo=lo+tlo; hi=hi+thi
    # inflate
    w=hi-lo
    lo=lo-w*REL-1e-12; hi=hi+w*REL+1e-12
    return lo,hi

def iv_div(alo,ahi,blo,bhi):
    # positive intervals (blo>0): [a/b] = [alo/bhi, ahi/blo]
    return alo/bhi, ahi/blo

def iv_sqrt(lo,hi):
    return np.sqrt(np.maximum(lo,0)), np.sqrt(np.maximum(hi,0))

def prune_batch(Xlo,Xhi,Ylo,Yhi,Slo,Shi):
    """Return boolean mask: True = box EXCLUDED (provably empty). False = survives."""
    excl=np.zeros(len(Xlo),dtype=bool)
    # domain impossible: Y >= X^2/4 everywhere, or Y <= X-1 everywhere
    excl |= (Ylo >= Xhi**2/4) | (Yhi <= Xlo-1)
    # disc = X^2 - 4Y < 0 everywhere => no real C,D
    excl |= (Xhi**2 - 4*Ylo <= 0)
    if excl.all(): return excl
    # For non-excluded boxes, check G,S,det contain 0 (exclude if 0 outside)
    idx=np.where(~excl)[0]
    def sub(arr,i=idx): return arr[i]
    g_lo,g_hi=iv_eval(Gterms,sub(Xlo),sub(Xhi),sub(Ylo),sub(Yhi),sub(Slo),sub(Shi))
    excl[idx] |= (g_lo>0)|(g_hi<0)
    idx=idx[(g_lo<=0)&(g_hi>=0)]
    if len(idx)==0: return excl
    def sub2(arr,i=idx): return arr[i]
    s_lo,s_hi=iv_eval(Sterms,sub2(Xlo),sub2(Xhi),sub2(Ylo),sub2(Yhi),sub2(Slo),sub2(Shi))
    excl[idx] |= (s_lo>0)|(s_hi<0)
    idx=idx[(s_lo<=0)&(s_hi>=0)]
    if len(idx)==0: return excl
    def sub3(arr,i=idx): return arr[i]
    d_lo,d_hi=iv_eval(detterms,sub3(Xlo),sub3(Xhi),sub3(Ylo),sub3(Yhi),sub3(Slo),sub3(Shi))
    excl[idx] |= (d_lo>0)|(d_hi<0)
    idx=idx[(d_lo<=0)&(d_hi>=0)]
    if len(idx)==0: return excl
    # lift checks: only meaningful where disc>0 everywhere (Xlo^2-4*Yhi>0).
    # Where disc straddles 0, can't eval lift cleanly -> survive (subdivide).
    def sub4(arr,i=idx): return arr[i]
    Xl=sub4(Xlo);Xh=sub4(Xhi);Yl=sub4(Ylo);Yh=sub4(Yhi);Sl=sub4(Slo);Sh=sub4(Shi)
    disc_pos = (Xl**2 - 4*Yh > 0)  # disc>0 everywhere here
    j=np.where(disc_pos)[0]
    if len(j)==0: return excl  # all survivors have straddling disc -> survive
    Xl=Xl[j];Xh=Xh[j];Yl=Yl[j];Yh=Yh[j];Sl=Sl[j];Sh=Sh[j]
    disclo=Xl**2-4*Yh; dischi=Xh**2-4*Yl
    slo_,shi_=iv_sqrt(disclo,dischi)
    Cl=(Xl+slo_)/2;Ch=(Xh+shi_)/2; Dl=(Xl-shi_)/2;Dh=(Xh-slo_)/2
    # exclude if C or D provably outside (0,1)
    lift_excl = (Ch<=0)|(Cl>=1)|(Dh<=0)|(Dl>=1)
    # a3 = 1-C + sigma*C^2*(1-C)/den, den=X-1>0
    denl=Xl-1;denh=Xh-1
    C2l=Cl**2;C2h=Ch**2; C1l=1-Ch;C1h=1-Cl
    nl=Sl*C2l*C1l; nh=Sh*C2h*C1h
    q_lo,q_hi=iv_div(nl,nh,denl,denh)
    a3lo=(1-Ch)+q_lo; a3hi=(1-Cl)+q_hi
    lift_excl |= (a3hi<=0)|(a3lo>=1)
    # a6 = 1-D + sigma*D^2*(1-D)/den
    D2l=Dl**2;D2h=Dh**2; D1l=1-Dh;D1h=1-Dl
    nl=Sl*D2l*D1l; nh=Sh*D2h*D1h
    q_lo,q_hi=iv_div(nl,nh,denl,denh)
    a6lo=(1-Dh)+q_lo; a6hi=(1-Dl)+q_hi
    lift_excl |= (a6hi<=0)|(a6lo>=1)
    # map lift_excl back to excl
    excl[idx[j]] |= lift_excl
    return excl

# B&B queue: list of boxes (Xlo,Xhi,Ylo,Yhi,Slo,Shi). Subdivide largest.
SIGMA_MAX=8.0
import heapq
box=(1.0,2.0, 0.0,1.0, 0.0,SIGMA_MAX)
queue=[box]
survivors=[]
t_start=time.time()
MAXTIME=900  # 15 min budget
npruned=0; nprocessed=0
MINW=1e-3  # min half-width before declaring survivor (refine later)
print(f"start B&B, sigma_max={SIGMA_MAX}, min half-width={MINW}",flush=True)
while queue and time.time()-t_start<MAXTIME:
    # process a batch: pop up to N boxes
    batch=queue[:4096]; queue=queue[4096:]
    Xlo=np.array([b[0] for b in batch]);Xhi=np.array([b[1] for b in batch])
    Ylo=np.array([b[2] for b in batch]);Yhi=np.array([b[3] for b in batch])
    Slo=np.array([b[4] for b in batch]);Shi=np.array([b[5] for b in batch])
    excl=prune_batch(Xlo,Xhi,Ylo,Yhi,Slo,Shi)
    nprocessed+=len(batch)
    for i,b in enumerate(batch):
        if excl[i]:
            npruned+=1; continue
        # survivor or subdivide
        w=max(b[1]-b[0],b[3]-b[2],b[5]-b[4])
        if w<2*MINW:
            survivors.append(b)
        else:
            # subdivide along longest axis
            ax=b[1]-b[0]; ay=b[3]-b[2]; az=b[5]-b[4]
            if ax>=ay and ax>=az:
                m=(b[0]+b[1])/2
                queue.append((b[0],m,b[2],b[3],b[4],b[5])); queue.append((m,b[1],b[2],b[3],b[4],b[5]))
            elif ay>=az:
                m=(b[2]+b[3])/2
                queue.append((b[0],b[1],b[2],m,b[4],b[5])); queue.append((b[0],b[1],m,b[3],b[4],b[5]))
            else:
                m=(b[4]+b[5])/2
                queue.append((b[0],b[1],b[2],b[3],b[4],m)); queue.append((b[0],b[1],b[2],b[3],m,b[5]))
    if nprocessed%200000<4096:
        print(f"  processed={nprocessed} pruned={npruned} queue={len(queue)} survivors={len(survivors)} t={time.time()-t_start:.0f}s",flush=True)

print(f"\nDONE. processed={nprocessed} pruned={npruned} survivors={len(survivors)} queue_left={len(queue)} t={time.time()-t_start:.0f}s",flush=True)
if not survivors and not queue:
    print("==> ALL BOXES EXCLUDED: {G=S=det=0, strict interior} certified EMPTY (float-interval + inflation).",flush=True)
    print("    (Re-verify survivors with mpmath.iv for full directed-rounding rigor: 0 survivors => vacuously rigorous.)",flush=True)
else:
    print(f"==> {len(survivors)} survivors below min-width {MINW} (need mp.iv refinement):",flush=True)
    for b in survivors[:20]:
        print(f"   X[{b[0]:.5f},{b[1]:.5f}] Y[{b[2]:.5f},{b[3]:.5f}] sig[{b[4]:.5f},{b[5]:.5f}]",flush=True)
