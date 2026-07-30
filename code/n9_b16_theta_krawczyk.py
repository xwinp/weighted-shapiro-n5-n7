#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B.16 empty-set certificate: 3D Krawczyk interval-Newton in (C,D,sigma) with Theta.

Why Krawczyk converges here (unlike range-pruning): near the G=S=0 curve, G(m),S(m)
are small but Theta(m) != 0 (numerically verified), so F(m)=(G,S,Theta) is large in
the 3rd component -> center m - J(m)^{-1}F(m) escapes the box -> K cap B = empty ->
one-step exclusion. Range-pruning failed because it used the loose Theta *interval*.

F=(G_cd, S_cd, Theta).  Gradients: G_cd,S_cd via Poly.diff (polynomial); Theta via
sympy.diff (rational).  All evaluated through lambdify+'mpmath', which transparently
supports both mp.mpf (point, for F(m),J(m)) and iv.mpf (interval, for J(B)).

numpy cheap-prune (domain + sigma-bound (6.1) + lift positivity) eliminates
non-strict-interior boxes fast before the expensive per-box Krawczyk.
"""
import time
from pathlib import Path
import numpy as np
import sympy as sp
import mpmath as mp
from mpmath import iv
iv.dps=18; mp.mp.dps=30

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
print("diff Theta...", flush=True); t0=time.time()
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
print(f"  Theta built {time.time()-t0:.1f}s; diff...", flush=True)
dTheta = [sp.diff(Theta, v) for v in (C,D,sig)]
dGcd = [sp.diff(Gcd, v) for v in (C,D,sig)]
dScd = [sp.diff(Scd, v) for v in (C,D,sig)]
print(f"  diffs done {time.time()-t0:.1f}s; lambdify...", flush=True)
fG = sp.lambdify((C,D,sig), Gcd, 'mpmath')
fS = sp.lambdify((C,D,sig), Scd, 'mpmath')
fT = sp.lambdify((C,D,sig), Theta, 'mpmath')
fJ = [[sp.lambdify((C,D,sig), g, 'mpmath') for g in dGcd],
      [sp.lambdify((C,D,sig), g, 'mpmath') for g in dScd],
      [sp.lambdify((C,D,sig), g, 'mpmath') for g in dTheta]]
print(f"  lambdify done {time.time()-t0:.1f}s", flush=True)

# ---- numpy cheap prune: domain + sigma-bound + a3/a6 lift (only where gap>0 everywhere) ----
def _mul(al,ah,bl,bh):
    t1=al*bl;t2=al*bh;t3=ah*bl;t4=ah*bh
    return np.minimum(np.minimum(t1,t2),np.minimum(t3,t4)), np.maximum(np.maximum(t1,t2),np.maximum(t3,t4))
def cheap_excl(Clo,Chi,Dlo,Dhi,Slo,Shi):
    excl=np.zeros(len(Clo),dtype=bool)
    excl |= (Dlo>=Chi)|(Dhi<=0)|(Clo>=1)|(Chi+Dhi<=1)|(Shi<=0)
    if excl.all(): return excl
    glo=Clo+Dlo-1; ghi=Chi+Dhi-1
    # sigma-bound (6.1) and a3/a6 lift only valid where gap>0 everywhere (glo>0)
    valid = glo>0
    if np.any(valid):
        with np.errstate(divide='ignore',invalid='ignore'):
            bndC=(1-Clo)*Chi**2/glo   # upper bound C(1-C)/gap (gap>0 -> /glo for upper)
            bndD=(1-Dlo)*Dhi**2/glo
            bnd_hi=np.where(valid, np.minimum(bndC,bndD), np.inf)
            excl |= valid & (Slo>=bnd_hi)
        if excl.all(): return excl
        # a3,a6 interval (gap>0 so /gh,/gl valid positive)
        idx=np.where(valid & ~excl)[0]
        if len(idx)>0:
            Cl=Clo[idx];Ch=Chi[idx];Dl=Dlo[idx];Dh=Dhi[idx];Sl=Slo[idx];Sh=Shi[idx]
            gl=glo[idx];gh=ghi[idx]
            C1l=1-Ch; C1h=1-Cl; C2l=Cl**2; C2h=Ch**2
            sC2l,sC2h=_mul(Sl,Sh,C2l,C2h); sC2C1l,sC2C1h=_mul(sC2l,sC2h,C1l,C1h)
            a3l=(1-Ch)+sC2C1l/gh; a3h=(1-Cl)+sC2C1h/gl
            D1l=1-Dh; D1h=1-Dl; D2l=Dl**2; D2h=Dh**2
            sD2l,sD2h=_mul(Sl,Sh,D2l,D2h); sD2D1l,sD2D1h=_mul(sD2l,sD2h,D1l,D1h)
            a6l=(1-Dh)+sD2D1l/gh; a6h=(1-Dl)+sD2D1h/gl
            # exclude if a3 or a6 provably outside (0,1); straddling -> subdivide
            lex=(a3h<=0)|(a3l>=1)|(a6h<=0)|(a6l>=1)
            excl[idx]|=lex
    return excl

# ---- Krawczyk per-box (mpmath.iv) ----
def iv_eval(f, c, d, s):
    return f(iv.mpf([float(c.a),float(c.b)]), iv.mpf([float(d.a),float(d.b)]), iv.mpf([float(s.a),float(s.b)]))
def pt_eval(f, c, d, s):
    return f(mp.mpf(c), mp.mpf(d), mp.mpf(s))

def krawczyk(box):
    Clo,Chi,Dlo,Dhi,Slo,Shi = box
    # cheap excl already done before calling; here do Krawczyk
    cm=(Clo+Chi)/2; dm=(Dlo+Dhi)/2; sm=(Slo+Shi)/2
    # point F(m), J(m)
    try:
        Fm = [pt_eval(f, cm,dm,sm) for f in (fG,fS,fT)]
        Jm = [[pt_eval(fJ[i][j], cm,dm,sm) for j in range(3)] for i in range(3)]
    except Exception:
        return 'sub'
    # invert J(m)
    try:
        JM = mp.matrix([[Jm[i][j] for j in range(3)] for i in range(3)])
        Jinv = JM**(-1)
    except Exception:
        return 'sub'
    # if det ~0 -> subdivide
    detJ = mp.det(JM)
    if abs(detJ) < mp.mpf('1e-12'):
        return 'sub'
    if not all(mp.isfinite(x) for x in Fm):
        return 'sub'
    # interval J(B)
    cB=iv.mpf([Clo,Chi]); dB=iv.mpf([Dlo,Dhi]); sB=iv.mpf([Slo,Shi])
    try:
        JB = [[fJ[i][j](cB,dB,sB) for j in range(3)] for i in range(3)]
    except Exception:
        return 'sub'
    # check any inf in JB -> subdivide
    def is_inf(ivv):
        return (ivv.a==mp.ninf or ivv.b==mp.inf) if hasattr(ivv,'a') else False
    if any(is_inf(JB[i][j]) for i in range(3) for j in range(3)):
        return 'sub'
    # center = m - Jinv @ Fm
    center=[cm,dm,sm]
    for i in range(3):
        s=mp.mpf('0')
        for k in range(3):
            s += Jinv[i,k]*Fm[k]
        center[i]=center[i]-s
    # (I - Jinv @ JB) @ (B - m), interval
    Blo=[Clo,Dlo,Slo]; Bhi=[Chi,Dhi,Shi]; Bm=[cm,dm,sm]
    Klo=[mp.mpf(0)]*3; Khi=[mp.mpf(0)]*3
    for i in range(3):
        # accumulate interval sum over k of (delta_ik - Jinv[i,k]*JB[k,j])*(B_j - m_j)
        lo=mp.mpf(0); hi=mp.mpf(0)
        for j in range(3):
            # coefficient interval for (i,j): delta_ij - sum_k Jinv[i,k]*JB[k,j]
            clo=mp.mpf(0); chi=mp.mpf(0)
            for k in range(3):
                # Jinv[i,k] (mpf, point) * JB[k,j] (interval)
                jv=Jinv[i,k]
                ja=float(JB[k][j].a); jb=float(JB[k][j].b)
                if jv>=0:
                    plo=jv*ja; phi=jv*jb
                else:
                    plo=jv*jb; phi=jv*ja
                clo+=plo; chi+=phi
            if i==j:
                clo+=1; chi+=1
            # (B_j - m_j) interval
            dlo=Blo[j]-Bm[j]; dhi=Bhi[j]-Bm[j]
            # product interval [clo,chi]*[dlo,dhi]
            ps=[clo*dlo,clo*dhi,chi*dlo,chi*dhi]
            lo+=min(ps); hi+=max(ps)
        Klo[i]=center[i]+lo; Khi[i]=center[i]+hi
    # K cap B
    inter_lo=[max(Klo[i],Blo[i]) for i in range(3)]
    inter_hi=[min(Khi[i],Bhi[i]) for i in range(3)]
    if any(inter_lo[i]>inter_hi[i] for i in range(3)):
        return 'excl'
    # K subset B ?
    if all(Klo[i]>=Blo[i] and Khi[i]<=Bhi[i] for i in range(3)):
        return 'root'   # unique root in box (should not happen; check lift/P<9)
    return 'sub'

# ---- B&B ----
def run():
  SMAX=4.0
  queue=[(0.5,1.0, 0.0,0.5, 0.0,SMAX)]
  survivors=[]; roots=[]; npruned=0; nexcl=0; nproc=0; t0=time.time()
  MAXTIME=600; MINW=1e-5
  print(f"start Krawczyk B&B, SMAX={SMAX} minw={MINW}", flush=True)
  while queue and time.time()-t0<MAXTIME:
    # pop a batch, numpy cheap-excl
    batch=queue[:1024]; queue=queue[1024:]
    Clo=np.array([b[0] for b in batch]);Chi=np.array([b[1] for b in batch])
    Dlo=np.array([b[2] for b in batch]);Dhi=np.array([b[3] for b in batch])
    Slo=np.array([b[4] for b in batch]);Shi=np.array([b[5] for b in batch])
    excl=cheap_excl(Clo,Chi,Dlo,Dhi,Slo,Shi)
    nproc+=len(batch)
    for i,b in enumerate(batch):
        if excl[i]: npruned+=1; continue
        st=krawczyk(b)
        if st=='excl': nexcl+=1
        elif st=='root': roots.append(b)
        else:
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
    if nproc%2048<1024:
        print(f"  proc={nproc} cheap_pruned={npruned} kraw_excl={nexcl} roots={len(roots)} surv={len(survivors)} Q={len(queue)} t={time.time()-t0:.0f}s", flush=True)
  print(f"\nDONE proc={nproc} cheap_pruned={npruned} kraw_excl={nexcl} roots={len(roots)} surv={len(survivors)} Q={len(queue)} t={time.time()-t0:.0f}s", flush=True)
  if not survivors and not queue and not roots:
      print("==> CERTIFIED EMPTY: {G=S=Theta=0} strict interior.", flush=True)
  if roots:
      print(f"==> {len(roots)} Krawczyk-confirmed roots (check lift/P<9):", flush=True)
      for b in roots[:20]: print(f"   C[{b[0]},{b[1]}] D[{b[2]},{b[3]}] sig[{b[4]},{b[5]}]", flush=True)
  if survivors:
      print(f"==> {len(survivors)} survivors:", flush=True)
      for b in survivors[:20]: print(f"   C[{b[0]:.6f},{b[1]:.6f}] D[{b[2]:.6f},{b[3]:.6f}] sig[{b[4]:.6f},{b[5]:.6f}]", flush=True)

if __name__=='__main__':
    run()
