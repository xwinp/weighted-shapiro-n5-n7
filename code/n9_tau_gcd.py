#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test the clean B.16 certificate: in regularized coords (v,xi,tau) with the
CORRECT sigma = tau*v/(d(1-d)), does {G^-=G^+=0} force tau=1 (i.e. a6=1)?

If gcd_tau(G^-, G^+) over Q(v,xi) is linear in tau with root tau=1, then the
ONLY common tau is tau=1, so {G^-=G^+=0} has no point with tau in (0,1) ->
{G^-=G^+=0} cap Omega_np = empty -> B.16 closed.
"""
import time
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'
c, d, s = sp.symbols('c d sigma')
v, xi, tau = sp.symbols('v xi tau')

# compact G^-, G^+ in (c,d,sigma)
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
print(f"G^- deg={Gm.total_degree()} sigma-deg={Gm.degree(s)}; G^+ deg={Gp.total_degree()} sigma-deg={Gp.degree(s)}", flush=True)

# regularization with CORRECT sigma = tau*v/(d(1-d))
cc = (1 + v + (1-v)*xi)/2
dd = (1 + v - (1-v)*xi)/2
sig = tau*v/(dd*(1-dd))
print("substituting regularization (correct sigma=tau v/(d(1-d)))...", flush=True); t0=time.time()
Gm_r = sp.together(Gm.as_expr().subs({c:cc, d:dd, s:sig}))
Gp_r = sp.together(Gp.as_expr().subs({c:cc, d:dd, s:sig}))
nm, dm = sp.fraction(Gm_r); np_, dp = sp.fraction(Gp_r)
Gmtilde = sp.expand(nm); Gptilde = sp.expand(np_)
Pm = sp.Poly(Gmtilde, v, xi, tau); Pp = sp.Poly(Gptilde, v, xi, tau)
print(f"  Gtilde^-: deg={Pm.total_degree()} terms={len(Pm.terms())} tau-deg={Pm.degree(tau)}", flush=True)
print(f"  Gtilde^+: deg={Pp.total_degree()} terms={len(Pp.terms())} tau-deg={Pp.degree(tau)}", flush=True)
print(f"  den_m: {sp.factor(dm)}", flush=True)
print(f"  den_p: {sp.factor(dp)}", flush=True)
print(f"  build t={time.time()-t0:.1f}s", flush=True)

# gcd in tau over Q(v,xi): treat as polys in tau with coeff in Q(v,xi)
print("computing gcd_tau(G^-,G^+) over Q(v,xi)...", flush=True); t1=time.time()
Pm_t = sp.Poly(Gmtilde, tau, v, xi, domain='QQ')
Pp_t = sp.Poly(Gptilde, tau, v, xi, domain='QQ')
g = sp.gcd(Pm_t, Pp_t)
print(f"  gcd deg in tau = {g.degree(tau) if g!=0 else 'zero'}  t={time.time()-t1:.1f}s", flush=True)
print(f"  gcd (as expr, head): {str(g.as_expr())[:300]}", flush=True)
# if linear in tau, root = -B/A; check if root == 1
if g != 0 and g.degree(tau) == 1:
    # g = A(v,xi)*tau + B(v,xi)
    A = g.nth(1); B = g.nth(0)
    print(f"  linear: A={sp.factor(A)}", flush=True)
    print(f"           B={sp.factor(B)}", flush=True)
    root = sp.simplify(-B/A)
    print(f"  root tau = -B/A = {sp.factor(root)}", flush=True)
    print(f"  root == 1? {sp.simplify(root-1)==0}", flush=True)
print("DONE", flush=True)
