#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Krawczyk interval-Newton branch-and-prune: certify {G=S=det=0, strict interior} EMPTY.

Krawczyk step for F=(G,S,det), J=Jacobian, midpoint m, box B:
  K = m - J(m)^{-1} F(m) + (I - J(m)^{-1} J(B)) (B - m)
  K∩B=∅  => no root in B (EXCLUDE)
  K⊂B     => unique root in B (SURVIVOR -> check lift)
Excludes root-free boxes in one step regardless of polynomial degree (uses midpoint).
Domain/lift range-pruning (cheap) applied first.
"""
import re, numpy as np, time, sys
from pathlib import Path
import sympy as sp
HERE=Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'
Xs,Ys,ss=sp.symbols('X Y s')
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
print("building polys...",flush=True)
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
# term lists: F components and Jacobian rows
def tl(P): return [(tuple(map(int,k)),int(v)) for k,v in P.terms()]
Fterms=[tl(G),tl(S),tl(detP)]                 # F = (G,S,det)
Jterms=[[tl(Gx),tl(Gy),tl(Gs)],[tl(Sx),tl(Sy),tl(Ss)],
        [tl(detP.diff(Xs)),tl(detP.diff(Ys)),tl(detP.diff(ss))]]
print(f"built in {time.time()-t0:.1f}s. det deg={detP.total_degree()}",flush=True)

# point evaluation (float)
def peval(terms,x,y,s):
    t=0.0
    for (ax,ay,asg),c in terms: t+=c*(x**ax)*(y**ay)*(s**asg)
    return t
# interval evaluation (positive vars), numpy batched
REL=1e-10
def iveval(terms,Xlo,Xhi,Ylo,Yhi,Slo,Shi):
    lo=np.zeros_like(Xlo); hi=np.zeros_like(Xlo)
    for (ax,ay,asg),c in terms:
        xlo=1.0 if ax==0 else Xlo**ax; xhi=1.0 if ax==0 else Xhi**ax
        ylo=1.0 if ay==0 else Ylo**ay; yhi=1.0 if ay==0 else Yhi**ay
        slo=1.0 if asg==0 else Slo**asg; shi=1.0 if asg==0 else Shi**asg
        mlo=xlo*ylo*slo; mhi=xhi*yhi*shi
        if c>=0: lo=lo+c*mlo; hi=hi+c*mhi
        else: lo=lo+c*mhi; hi=hi+c*mlo
    w=hi-lo; lo=lo-w*REL-1e-12; hi=hi+w*REL+1e-12
    return lo,hi

def iv_div(alo,ahi,blo,bhi):  # blo>0
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(bhi>0,alo/bhi,-np.inf), np.where(blo>0,ahi/blo,np.inf)
def iv_sqrt(lo,hi):
    return np.sqrt(np.maximum(lo,0)), np.sqrt(np.maximum(hi,0))

# cheap domain/lift exclusion (range-based, exclude only provably empty)
def cheap_excl(Xlo,Xhi,Ylo,Yhi,Slo,Shi):
    excl=np.zeros(len(Xlo),dtype=bool)
    excl |= (Ylo>=Xhi**2/4)|(Yhi<=Xlo-1)         # domain
    excl |= (Xhi**2-4*Ylo<=0)                     # disc<0 everywhere
    # lift only where disc>0 everywhere
    idx=np.where(~excl)[0]
    if len(idx)==0: return excl
    Xl=Xlo[idx];Xh=Xhi[idx];Yl=Ylo[idx];Yh=Yhi[idx]
    dp=(Xl**2-4*Yh>0)
    j=np.where(dp)[0]
    if len(j)==0: return excl
    Xl=Xl[j];Xh=Xh[j];Yl=Yl[j];Yh=Yh[j];Sl=Slo[idx][j];Sh=Shi[idx][j]
    disclo=Xl**2-4*Yh; dischi=Xh**2-4*Yl
    slo,shi=iv_sqrt(disclo,dischi)
    Cl=(Xl+slo)/2;Ch=(Xh+shi)/2; Dl=(Xl-shi)/2;Dh=(Xh-slo)/2
    lex=(Ch<=0)|(Cl>=1)|(Dh<=0)|(Dl>=1)
    denl=Xl-1;denh=Xh-1
    # a3
    C2l=Cl**2;C2h=Ch**2; C1l=1-Ch;C1h=1-Cl
    nl=Sl*C2l*C1l; nh=Sh*C2h*C1h
    qlo,qhi=iv_div(nl,nh,denl,denh)
    a3lo=(1-Ch)+qlo; a3hi=(1-Cl)+qhi
    lex|=(a3hi<=0)|(a3lo>=1)
    # a6
    D2l=Dl**2;D2h=Dh**2; D1l=1-Dh;D1h=1-Dl
    nl=Sl*D2l*D1l; nh=Sh*D2h*D1h
    qlo,qhi=iv_div(nl,nh,denl,denh)
    a6lo=(1-Dh)+qlo; a6hi=(1-Dl)+qhi
    lex|=(a6hi<=0)|(a6lo>=1)
    excl[idx[j]]|=lex
    return excl

def krawczyk_batch(boxes):
    """boxes: list of (Xlo,Xhi,Ylo,Yhi,Slo,Shi). Return list of (status, box) where
    status in {'excl','root','sub'}.  Batched numpy."""
    n=len(boxes)
    Xlo=np.array([b[0] for b in boxes]);Xhi=np.array([b[1] for b in boxes])
    Ylo=np.array([b[2] for b in boxes]);Yhi=np.array([b[3] for b in boxes])
    Slo=np.array([b[4] for b in boxes]);Shi=np.array([b[5] for b in boxes])
    # cheap exclude
    excl=cheap_excl(Xlo,Xhi,Ylo,Yhi,Slo,Shi)
    statuses=[]
    # midpoints
    Xm=(Xlo+Xhi)/2; Ym=(Ylo+Yhi)/2; Sm=(Slo+Shi)/2
    # F(m) shape (3,n)
    Fm=np.array([[peval(ft,x,y,s) for x,y,s in zip(Xm,Ym,Sm)] for ft in Fterms])
    # J(m) shape (3,3,n) float
    Jm=np.array([[[peval(jt,x,y,s) for x,y,s in zip(Xm,Ym,Sm)] for jt in row] for row in Jterms])
    # J(B) interval: for each (i,k) component, (lo,hi) over boxes -> shape (3,3,n)
    JBlo=np.zeros((3,3,n)); JBhi=np.zeros((3,3,n))
    for i in range(3):
        for k in range(3):
            lo,hi=iveval(Jterms[i][k],Xlo,Xhi,Ylo,Yhi,Slo,Shi)
            JBlo[i,k]=lo; JBhi[i,k]=hi
    # invert J(m) per box
    JmT=np.transpose(Jm,(2,0,1))  # (n,3,3)
    out=[]
    for idx in range(n):
        if excl[idx]:
            out.append(('excl',boxes[idx])); continue
        M=JmT[idx]
        try:
            Minv=np.linalg.inv(M)
        except np.linalg.LinAlgError:
            out.append(('sub',boxes[idx])); continue
        if not np.all(np.isfinite(Minv)):
            out.append(('sub',boxes[idx])); continue
        fm=Fm[:,idx]  # (3,)
        # center c = m - Minv @ fm
        c=np.array([Xm[idx],Ym[idx],Sm[idx]])-Minv@fm
        # (I - Minv @ J(B)) @ (B - m), interval
        # Blo interval matrix (3,3) for this box
        JBl=JBlo[:,:,idx]; JBh=JBhi[:,:,idx]
        # Minv @ J(B): interval = [Minv@JBl, Minv@JBh] is wrong (Minv has signs);
        # compute interval product properly: for each entry sum_k Minv[i,k]*JB[k,j]
        # Minv float, JB interval -> product interval [min,max] over sign of Minv[i,k]
        Plo=np.zeros((3,3)); Phi=np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                vals_lo=0.0; vals_hi=0.0
                for k in range(3):
                    mik=Minv[i,k]
                    if mik>=0:
                        vals_lo+=mik*JBl[k,j]; vals_hi+=mik*JBh[k,j]
                    else:
                        vals_lo+=mik*JBh[k,j]; vals_hi+=mik*JBl[k,j]
                Plo[i,j]=vals_lo; Phi[i,j]=vals_hi
        # I - P
        Elo=1.0-Phi; Ehi=1.0-Plo  # identity minus interval [Plo,Phi]
        # (B-m) interval vector
        dlo=np.array([Xlo[idx]-Xm[idx],Ylo[idx]-Ym[idx],Slo[idx]-Sm[idx]])
        dhi=np.array([Xhi[idx]-Xm[idx],Yhi[idx]-Ym[idx],Shi[idx]-Sm[idx]])
        # E @ d interval
        Klo=np.zeros(3); Khi=np.zeros(3)
        for i in range(3):
            vlo=0.0; vhi=0.0
            for k in range(3):
                eik_lo=Elo[i,k]; eik_hi=Ehi[i,k]
                # eik * dk interval
                terms=[eik_lo*dlo[k],eik_lo*dhi[k],eik_hi*dlo[k],eik_hi*dhi[k]]
                vlo+=min(terms); vhi+=max(terms)
            Klo[i]=vlo; Khi[i]=vhi
        # K = c + [Klo,Khi]
        K_lo=c+Klo; K_hi=c+Khi
        # B interval
        Blo=np.array([Xlo[idx],Ylo[idx],Slo[idx]]); Bhi=np.array([Xhi[idx],Yhi[idx],Shi[idx]])
        # K∩B
        inter_lo=np.maximum(K_lo,Blo); inter_hi=np.minimum(K_hi,Bhi)
        if np.any(inter_lo>inter_hi):
            out.append(('excl',boxes[idx])); continue
        # K⊂B ?
        inside=np.all(K_lo>=Blo-1e-15) and np.all(K_hi<=Bhi+1e-15)
        if inside:
            out.append(('root',boxes[idx]))
        else:
            out.append(('sub',boxes[idx]))
    return out

SIGMA_MAX=8.0
MINW=1e-5
MAXTIME=900
queue=[(1.0,2.0,0.0,1.0,0.0,SIGMA_MAX)]
survivors=[]; roots=[]
npruned=0; nproc=0; t0=time.time()
print(f"start Krawczyk B&B sigma_max={SIGMA_MAX} minw={MINW}",flush=True)
while queue and time.time()-t0<MAXTIME:
    batch=queue[:2048]; queue=queue[2048:]
    res=krawczyk_batch(batch)
    nproc+=len(batch)
    for st,b in res:
        if st=='excl': npruned+=1
        elif st=='root': roots.append(b)
        else:
            w=max(b[1]-b[0],b[3]-b[2],b[5]-b[4])
            if w<2*MINW:
                survivors.append(b)
            else:
                ax=b[1]-b[0];ay=b[3]-b[2];az=b[5]-b[4]
                if ax>=ay and ax>=az:
                    m=(b[0]+b[1])/2; queue.append((b[0],m,b[2],b[3],b[4],b[5]));queue.append((m,b[1],b[2],b[3],b[4],b[5]))
                elif ay>=az:
                    m=(b[2]+b[3])/2; queue.append((b[0],b[1],b[2],m,b[4],b[5]));queue.append((b[0],b[1],m,b[3],b[4],b[5]))
                else:
                    m=(b[4]+b[5])/2; queue.append((b[0],b[1],b[2],b[3],b[4],m));queue.append((b[0],b[1],b[2],b[3],m,b[5]))
    if nproc%50000<2048:
        print(f"  proc={nproc} pruned={npruned} roots={len(roots)} surv={len(survivors)} Q={len(queue)} t={time.time()-t0:.0f}s",flush=True)
print(f"\nDONE proc={nproc} pruned={npruned} roots={len(roots)} surv={len(survivors)} Q={len(queue)} t={time.time()-t0:.0f}s",flush=True)
if roots:
    print(f"==> {len(roots)} Krawczyk-confirmed roots (check lift/P<9):",flush=True)
    for b in roots[:30]: print(f"   X[{b[0]:.6f},{b[1]:.6f}] Y[{b[2]:.6f},{b[3]:.6f}] sig[{b[4]:.6f},{b[5]:.6f}]",flush=True)
if not roots and not survivors and not queue:
    print("==> ALL BOXES EXCLUDED (Krawczyk): {G=S=det=0, strict interior} certified EMPTY.",flush=True)
elif survivors:
    print(f"==> {len(survivors)} unresolved survivors (need finer MINW):",flush=True)
    for b in survivors[:20]: print(f"   X[{b[0]:.6f},{b[1]:.6f}] Y[{b[2]:.6f},{b[3]:.6f}] sig[{b[4]:.6f},{b[5]:.6f}]",flush=True)
