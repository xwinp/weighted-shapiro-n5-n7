#!/usr/bin/env python3
"""
Derive the S1 beta-reduction locally (to cross-check GPT and enable resultant+Sturm).
S1: x=(0,1,b,c,d,e,f), path ratios rho_1=b, rho_2=c/b, rho_3=d/c, rho_4=e/d, rho_5=f/e.
P = sum_{i=1}^4 1/(rho_i*(p+q*rho_{i+1}))  +  1/(p*rho_5)  +  (prod rho_i)/q,   q=1-p.
KKT in log-ratio coords y_i=log rho_i:  g_i := rho_i * dP/drho_i = 0  (i=1..5).
beta_i = q*rho_{i+1}/(p+q*rho_{i+1}) for i=1..4  (interior edges).
Goal: eliminate p (and rho) to get p-free relations among beta (the curve E2,E3,E4),
and the H_S1(beta) diagonal entries (GPT A.4): H_jj = 1 + a_j + a_{j-1}*beta_{j-1}*(2*beta_{j-1}-1).
Strategy: compute g_i, substitute rho_{i+1} = (p/q)*beta_i/(1-beta_i), simplify, and
look for p-free combinations. Also compute the log-Hessian diagonal numerically at the
known S1 branch to confirm which diagonal is negative (certificate target in beta-frame).
"""
import sympy as sp
import numpy as np
from scipy.optimize import root
n=7
rho = sp.symbols('rho1:6', positive=True)  # rho1..rho5
p = sp.symbols('p', positive=True); q = 1-p
prod_rho = sp.prod(rho)
# P in rho-coords
P = sum(1/(rho[i]*(p+q*rho[i+1])) for i in range(4)) + 1/(p*rho[4]) + prod_rho/q
P = sp.together(P)
# KKT g_i = rho_i * dP/drho_i
g = [sp.simplify(sp.together(rho[i]*sp.diff(P, rho[i]))) for i in range(5)]
print("=== KKT g_i = rho_i dP/drho_i (numerators) ===")
gnum = []
for i in range(5):
    num = g[i].as_numer_denom()[0]
    gnum.append(sp.expand(num))
    print(f"  g_{i+1} num: total_deg_in_rho={sp.total_degree(num,*rho)}, has_p={'yes' if p in num.free_symbols else 'no'}")

# beta_i = q*rho_{i+1}/(p+q*rho_{i+1}), i=1..4  -> rho_{i+1} = (p/q)*beta_i/(1-beta_i)
beta = sp.symbols('beta1:5', positive=True)  # beta1..beta4
subs_rho = {rho[i+1]: (p/q)*beta[i]/(1-beta[i]) for i in range(4)}  # rho2..rho5 in terms of beta1..beta4
# but rho1 is free (the "height"/scale). Since P is degree-0, rho1 should drop out of KKT.
# Substitute rho2..rho5
g_sub = [sp.simplify(sp.together(gnum[i].subs(subs_rho))) for i in range(5)]
print("\n=== g_i after substituting rho_{i+1}=(p/q)*beta_i/(1-beta_i) ===")
for i in range(5):
    print(f"  g_{i+1}: has_rho1={'yes' if rho[0] in g_sub[i].free_symbols else 'no(degree0 OK)'}, has_p={'yes' if p in g_sub[i].free_symbols else 'NO (p-free!)'}")

# Factor out powers of p,q,rho1 from each g_i to find the p-free core
print("\n=== p-free core of each g_i (factor out p,q,rho1 powers) ===")
for i in range(5):
    f = sp.factor(g_sub[i])
    print(f"  g_{i+1} factored: {f}")
