#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Variable regularization of {G^-=G^+=0} per GPT: v=c+d-1, u=c-d=(1-v)*xi,
sigma = tau*d*(1-d)/v, mapping strict interior to (v,xi,tau) in (0,1)^3.

Build Gtilde^-, Gtilde^+ (cleared-denominator polynomials in v,xi,tau), report
degree/terms, and verify interval eval is finite (no straddling denominators) on
a coarse grid -- the obstruction that stalled the (C,D,sigma) Krawczyk.
"""
import time
from pathlib import Path
import sympy as sp
import mpmath as mp
from mpmath import iv
iv.dps = 15; mp.mp.dps = 20

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
X, Y, sigma = sp.symbols('X Y sigma')
v, xi, tau = sp.symbols('v xi tau')

def load_small(name):
    s = sp.symbols("s")
    text = (HERE/name).read_text(encoding="utf-8").strip()
    return sp.Poly(sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma),
                   X, Y, sigma, domain=sp.ZZ)
Gm = load_small("nonpal_G_clean.txt")
Gp = load_small("nonpal_S_clean.txt")

# regularization
c = (1 + v + (1-v)*xi)/2
d = (1 + v - (1-v)*xi)/2
sig = tau*d*(1-d)/v
subs = {X: c+d, Y: c*d, sigma: sig}
print("substituting regularization into G^-, G^+ ...", flush=True); t0=time.time()
Gm_r = sp.together(Gm.as_expr().subs(subs))
Gp_r = sp.together(Gp.as_expr().subs(subs))
nm, dm = sp.fraction(Gm_r); np_, dp = sp.fraction(Gp_r)
Gmtilde = sp.expand(nm); Gptilde = sp.expand(np_)
Pm = sp.Poly(Gmtilde, v, xi, tau); Pp = sp.Poly(Gptilde, v, xi, tau)
print(f"  Gtilde^-: deg={Pm.total_degree()} terms={len(Pm.terms())} "
      f"degv={Pm.degree(v)} degxi={Pm.degree(xi)} degtau={Pm.degree(tau)}", flush=True)
print(f"  Gtilde^+: deg={Pp.total_degree()} terms={len(Pp.terms())} "
      f"degv={Pp.degree(v)} degxi={Pp.degree(xi)} degtau={Pp.degree(tau)}", flush=True)
print(f"  den_m factors sample: {sp.factor(dm) if dm!=1 else '1'}", flush=True)
print(f"  build t={time.time()-t0:.1f}s", flush=True)

# lambdify for interval eval
fGm = sp.lambdify((v,xi,tau), Gmtilde, 'mpmath')
fGp = sp.lambdify((v,xi,tau), Gptilde, 'mpmath')

# test: interval eval on a coarse 3D grid of (0,1)^3 sub-boxes; count inf/nan
def ivbox(lo,hi):
    return iv.mpf([float(lo),float(hi)])
ninf = mp.ninf; pinf = mp.inf
nboxes = 0; nfin = 0; nz_m = 0; nz_p = 0; nz_both = 0
import itertools
N = 8  # 8x8x8 = 512 boxes
for i,j,k in itertools.product(range(N),repeat=3):
    vlo,vhi = i/N,(i+1)/N
    xlo,xhi = j/N,(j+1)/N
    tlo,thi = k/N,(k+1)/N
    vb=ivbox(vlo,vhi); xb=ivbox(xlo,xhi); tb=ivbox(tlo,thi)
    nboxes += 1
    try:
        gm = fGm(vb,xb,tb); gp = fGp(vb,xb,tb)
        gmfin = not (gm.a==ninf or gm.b==pinf)
        gpfin = not (gp.a==ninf or gp.b==pinf)
        if gmfin: nfin+=1
        if gmfin and (gm.a<=0<=gm.b): nz_m+=1
        if gpfin and (gp.a<=0<=gp.b): nz_p+=1
        if gmfin and gpfin and (gm.a<=0<=gm.b) and (gp.a<=0<=gp.b): nz_both+=1
    except Exception as e:
        pass
print(f"\nN={N}^3={nboxes} boxes: finite_eval={nfin}  G^-∋0:{nz_m}  G^+∋0:{nz_p}  both∋0:{nz_both}", flush=True)
print("(both∋0 are candidate boxes needing subdivision; if denominators controlled, no inf)", flush=True)
print("DONE", flush=True)
