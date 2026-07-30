#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact (symbolic) proof that the H_C stationary value is P = C + 2 sqrt(A B).

Removes the old 'verified to 1e-14' gap.  The identity is elementary algebra:

S1 path ratios  rho_1, rho_2, ..., rho_5  with
    x = (0, 1, rho_1, rho_1 rho_2, ..., rho_1 rho_2 rho_3 rho_4 rho_5),
    rho_{i+1} = beta_i / (rho (1 - beta_i))   (i = 1..4),   p = 1/(1+rho), q = rho/(1+rho),
where  rho := q/p  and  beta_i in (0,1).  Write (u,v,w,z) = (beta_1,beta_2,beta_3,beta_4).

P in rho-coords is
    P = sum_{i=1}^4 1/(rho_i (p + q rho_{i+1}))  +  1/(p rho_5)  +  (prod rho_i)/q.

Claim:  P = C + A/rho_1 + B rho_1                                   (rational identity)
with
    A = (1-u)(1+rho),
    B = (1+rho) u v w z / (rho^5 (1-u)(1-v)(1-w)(1-z)),
    C = rho (1+rho) S,
    S = (1-v)(1-u)/u + (1-w)(1-v)/v + (1-z)(1-w)/w + (1-z)/z.

Then  g_1 := rho_1 dP/drho_1 = rho_1 (B - A/rho_1^2) = B rho_1 - A/rho_1,
so the first KKT equation g_1 = 0  <=>  rho_1^2 = A/B  (rho_1 > 0).
With rho_1 = sqrt(A/B):
    A/rho_1 + B rho_1 = A sqrt(B/A) + B sqrt(A/B) = 2 sqrt(A B),
hence  P = C + 2 sqrt(A B)   EXACTLY.

This script verifies the rational identity P - (C + A/rho_1 + B rho_1) == 0 and the
factorisation g_1 = rho_1 (B - A/rho_1^2) by exact SymPy expansion (no floats).
"""
import sympy as sp

# --- variables ---
r1 = sp.symbols('r1', positive=True)          # rho_1  (free height; >0 as a ratio of positive x's)
u, v, w, z = sp.symbols('u v w z', positive=True)   # beta_1..beta_4 in (0,1)
rho = sp.symbols('rho', positive=True)         # q/p (>0)
p = 1/(1+rho); q = rho/(1+rho)

# path ratios rho_2..rho_5 in terms of beta and rho
r2 = u/(rho*(1-u))
r3 = v/(rho*(1-v))
r4 = w/(rho*(1-w))
r5 = z/(rho*(1-z))
R = [r1, r2, r3, r4, r5]

# P in rho-coords (same definition as s1_beta_reduction.py / kkt_check)
P = sp.together(sum(1/(R[i]*(p + q*R[i+1])) for i in range(4))
                + 1/(p*R[4]) + sp.prod(R)/q)

# claimed closed-form pieces
A = (1-u)*(1+rho)
B = (1+rho)*u*v*w*z / (rho**5*(1-u)*(1-v)*(1-w)*(1-z))
S = (1-v)*(1-u)/u + (1-w)*(1-v)/v + (1-z)*(1-w)/w + (1-z)/z
C = rho*(1+rho)*S
P_closed = sp.together(C + A/r1 + B*r1)

# 1) rational identity  P == C + A/r1 + B r1
diff = sp.together(P - P_closed)
num, den = sp.fraction(diff)
num = sp.expand(num)
print("P == C + A/r1 + B*r1  (rational identity):",
      "TRUE (numerator identically zero)" if num == 0 else "FALSE")
assert num == 0, "closed-form rational identity FAILED"

# 2) first KKT equation  g_1 = r1 * dP/dr1  factorises as  r1*(B - A/r1^2)
g1 = sp.together(r1 * sp.diff(P, r1))
g1_claim = sp.together(r1*(B - A/r1**2))
d2 = sp.together(g1 - g1_claim)
n2, _ = sp.fraction(d2); n2 = sp.expand(n2)
print("g_1 == r1*(B - A/r1^2):",
      "TRUE (numerator identically zero)" if n2 == 0 else "FALSE")
assert n2 == 0, "g1 factorisation FAILED"

# 3) therefore g_1 = 0  <=>  r1^2 = A/B  (r1>0), and the value becomes C + 2 sqrt(AB).
#    The rho_1-dependent part is  A/r1 + B*r1 = (A + B*r1^2)/r1.  Verify by squaring
#    (a rational check, no radicals): set t = r1^2 = A/B, then
#        (A/r1 + B*r1)^2 = (A + B*t)^2 / t  ==  4 A B   at t = A/B,
#    and A/r1 + B*r1 > 0 on the admissible set (A,B,r1>0) => = +2 sqrt(AB).
t = sp.symbols('t', positive=True)            # t = r1^2
part_sq = sp.together(((A + B*t)/sp.sqrt(t))**2)   # (A/r1 + B*r1)^2 with r1=sqrt(t)
part_sq_sub = sp.together(part_sq.subs(t, A/B))    # impose g_1=0: t = A/B
n3, _ = sp.fraction(part_sq_sub - 4*A*B); n3 = sp.expand(n3)
print("at r1^2=A/B:  (A/r1 + B*r1)^2 == 4 A B:",
      "TRUE (numerator identically zero)" if n3 == 0 else "FALSE")
assert n3 == 0, "squared identity FAILED"
print("A/r1 + B*r1 > 0 on admissible set (A,B,r1>0)  =>  A/r1 + B*r1 = +2 sqrt(AB).")

# 4) sanity: on the admissible set A,B,radicand>0 so sqrt is the positive real root
#    (admissibility already certifies u,v,w,z in (0,1), rho>0, denominators>0).
print("\nA = (1-u)(1+rho)                       >0  on admissible set  (1-u>0, 1+rho>0)")
print("B = (1+rho)uvwz/(rho^5 (1-u)(1-v)(1-w)(1-z)) >0  on admissible set (all factors >0)")
print("=> rho_1 = +sqrt(A/B) is the unique positive solution of g_1=0.")
print("\nCONCLUSION:  P = C + 2 sqrt(A B)  is an EXACT algebraic identity on the H_C")
print("stationary set (no numerical verification involved).")
print("DONE-HC-CLOSEDFORM")
