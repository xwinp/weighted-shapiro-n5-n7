#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""For fixed tau in (0,1), solve G^-=G^+=0 in (v,xi), then print the full lift
vector (a1..a8) and which admissibility condition fails. Goal: find a UNIFORM
lift-failure mode on {G^-=G^+=0} cap {tau in (0,1)} -> cheap B.16 certificate.
"""
import random, mpmath as mp
import sympy as sp

mp.mp.dps = 18
c, d, s = sp.symbols('c d sigma')
v, xi, tau = sp.symbols('v xi tau')
vv = c + d - 1
Lc = vv + s*c**2; Uc = vv - s*c*(1-c); Bc = c**2*vv + (1-c)*Lc**2
Ld = vv + s*d**2; Ud = vv - s*d*(1-d); Bd = d**2*vv + (1-d)*Ld**2
Fc = c*vv**2*(1-c)*Lc**2 - s*Bc*(c*vv**2 - Uc*Bc)
Fd = d*vv**2*(1-d)*Ld**2 - s*Bd*(d*vv**2 - Ud*Bd)
Pc = sp.Poly(sp.expand(Fc), c, d, s, domain=sp.ZZ)
Pd = sp.Poly(sp.expand(Fd), c, d, s, domain=sp.ZZ)
FL = Pc.exquo(sp.Poly((c-1)*Lc, c, d, s, domain=sp.ZZ))
FR = Pd.exquo(sp.Poly((d-1)*Ld, c, d, s, domain=sp.ZZ))
Gm = sp.Poly(sp.expand(FL.as_expr()-FR.as_expr()), c, d, s, domain=sp.ZZ).exquo(sp.Poly(c-d, c, d, s, domain=sp.ZZ))
Gp = sp.Poly(sp.expand(FL.as_expr()+FR.as_expr()), c, d, s, domain=sp.ZZ)
cc = (1 + v + (1-v)*xi)/2; dd = (1 + v - (1-v)*xi)/2; sig = tau*v/(dd*(1-dd))
Gm_r = sp.together(Gm.as_expr().subs({c:cc, d:dd, s:sig})); Gp_r = sp.together(Gp.as_expr().subs({c:cc, d:dd, s:sig}))
Gmt = sp.expand(sp.fraction(Gm_r)[0]); Gpt = sp.expand(sp.fraction(Gp_r)[0])
fGm = sp.lambdify((v,xi,tau), Gmt, 'mpmath'); fGp = sp.lambdify((v,xi,tau), Gpt, 'mpmath')

def fullvec(vv, xx, tt):
    c1 = (1+vv+(1-vv)*xx)/2; d1 = (1+vv-(1-vv)*xx)/2; s1 = tt*vv/(d1*(1-d1)); g = c1+d1-1
    a3 = 1-c1+s1*c1**2*(1-c1)/g; a6 = 1-d1+s1*d1**2*(1-d1)/g
    a2 = 1-a3+s1*a3**2*(1-a3)/(a3+c1-1); a7 = 1-a6+s1*a6**2*(1-a6)/(a6+d1-1)
    return [1,a2,a3,c1,d1,a6,a7,1]

random.seed(11)
seen = {}
for t0 in [0.15, 0.3, 0.5, 0.7, 0.85, 0.95]:
    nfound = 0
    for _ in range(60):
        v0 = random.uniform(0.05,0.95); x0 = random.uniform(0.02,0.98)
        try:
            sol = mp.findroot(lambda vv,xx:(fGm(vv,xx,t0),fGp(vv,xx,t0)),
                              (mp.mpf(v0),mp.mpf(x0)), tol=1e-22, maxsteps=60)
            v1, x1 = float(sol[0]), float(sol[1])
            if not(0<v1<1 and 0<x1<1): continue
            r = max(abs(float(fGm(v1,x1,t0))),abs(float(fGp(v1,x1,t0))))
            if r>1e-7: continue
            # dedup
            key = (round(v1,3), round(x1,3))
            if key in seen: continue
            seen[key]=1
            a = fullvec(v1,x1,t0)
            fails = []
            for i in range(8):
                if not (0 < float(a[i]) < 1): fails.append(f'a{i+1}={float(a[i]):.4f}!in(0,1)')
            for i in range(7):
                if float(a[i]+a[i+1]) <= 1: fails.append(f'h{i+1}=a{i+1}+a{i+2}-1<=0')
            nfound += 1
            tag = "VALID" if not fails else "INVALID"
            print(f'tau={t0} v={v1:.4f} xi={x1:.4f} {tag} a={[round(float(z),4) for z in a]}', flush=True)
            if fails: print(f'      FAILS: {fails[:4]}', flush=True)
        except Exception:
            pass
    print(f'  [tau={t0}: {nfound} distinct solutions]', flush=True)
print('DONE', flush=True)
