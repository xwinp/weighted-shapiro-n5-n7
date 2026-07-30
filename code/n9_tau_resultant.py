#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compute R(v,xi) = Res_tau(G^-, G^+) for the regularized compact system, and
inspect its v-adic valuation. If R is v^k * (curve with no (0,1)^2 pt), the
certificate is: {G^-=G^+=0, tau in (0,1)} projects to v=0 boundary -> empty in
strict interior. Also test: does G^-=G^+=0 force v=0? (substitute v=0, check
common tau-roots exist; substitute small v>0, check none).
"""
import time, sympy as sp

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
print(f"Gtilde^-: tau-deg={sp.Poly(Gmt,tau).degree()}  Gtilde^+: tau-deg={sp.Poly(Gpt,tau).degree()}", flush=True)

# remove the common content v^3(v-1) (the gcd found earlier) to reduce
# Actually keep raw; compute resultant in tau over ZZ[v,xi].
Pm = sp.Poly(Gmt, tau, v, xi, domain=sp.ZZ)
Pp = sp.Poly(Gpt, tau, v, xi, domain=sp.ZZ)
print("computing Res_tau(G^-,G^+) ...", flush=True); t0=time.time()
R = sp.resultant(Pm, Pp, tau)
R = sp.expand(R)
PR = sp.Poly(R, v, xi, domain=sp.ZZ)
print(f"  done t={time.time()-t0:.1f}s. R: total_deg={PR.total_degree()} terms={len(PR.terms())} deg_v={PR.degree(v)} deg_xi={PR.degree(xi)}", flush=True)

# v-adic valuation
val = 0; cur = R
while True:
    q, r = sp.div(sp.Poly(cur, v, xi, domain=sp.ZZ), sp.Poly(v, v, xi, domain=sp.ZZ), domain='ZZ')
    if r.is_zero:
        val += 1; cur = q.as_expr()
    else:
        break
print(f"  v-adic valuation of R = {val}", flush=True)
Rred = sp.Poly(cur, v, xi, domain=sp.ZZ)
print(f"  R/v^{val}: total_deg={Rred.total_degree()} terms={len(Rred.terms())} deg_v={Rred.degree(v)} deg_xi={Rred.degree(xi)}", flush=True)
# is Rred free of v? (curve entirely on v=0)
print(f"  Rred is independent of v? {Rred.degree(v)==0}", flush=True)
if Rred.degree(v)==0:
    print(f"  => R = v^{val} * f(xi); curve projects to v=0 OR f(xi)=0.", flush=True)
    print(f"  f(xi) = {sp.factor(Rred.as_expr())}", flush=True)
# save
from pathlib import Path
out = Path(__file__).resolve().parent.parent/'paper'/'_gpt_artifacts'
(out/'nonpal_Rtau_vxi.txt').write_text(str(R)+"\n", encoding='utf-8')
print("  saved nonpal_Rtau_vxi.txt", flush=True)
print("DONE", flush=True)
