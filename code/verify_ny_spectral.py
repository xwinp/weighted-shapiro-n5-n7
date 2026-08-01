#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact verification of the N-Y spectral bridge cited in Lemma 2.3 (NY-2).

Checks (all exact SymPy, no floats):
  (A) The Hessian eigenvalue identity
        mu_k = 2|lambda_k|^2 - 2 Re(lambda_k) = 2(1-cos th_k)[2q cos th_k + 2p^2-3p+2]
      where lambda_k = p e^{i th} + q e^{2 i th}, q=1-p, c=cos th.
  (B) Disk equivalence:  mu_k > 0  <=>  |lambda_k - 1/2|^2 > 1/4
        (i.e. lambda_k strictly outside the closed disk |lambda-1/2|<=1/2).
      mu_k = 2(|lambda|^2 - Re lambda) = 2(|lambda-1/2|^2 - 1/4), so sign(mu_k)=sign(|lambda-1/2|^2-1/4).
  (C) Delta_n = 4 c_n^2 - 4 c_n - 7 < 0 for n=5,7  (c_n = -cos(pi/n)), so the dangerous-mode
      bracket 2q c + 2p^2-3p+2 is a quadratic in p with negative discriminant and positive
      leading coefficient, hence >0 for all p in (0,1)  =>  every mu_k > 0  =>  N-Y applies.
"""
import sympy as sp
p, c = sp.symbols('p c', real=True)   # c = cos(theta_k)
q = 1 - p

# lambda_k = p e^{i th} + q e^{2 i th};  Re(lambda)=p c + q (2c^2-1);  |lambda|^2 = p^2+q^2+2 p q c
Re_lam = p*c + q*(2*c**2 - 1)
abs2   = p**2 + q**2 + 2*p*q*c

mu_NY = sp.expand(2*abs2 - 2*Re_lam)                      # 2|lambda|^2 - 2 Re(lambda)
mu_paper = sp.expand(2*(1 - c)*(2*q*c + 2*p**2 - 3*p + 2))  # paper's formula (line 61)
print("(A) mu_NY == mu_paper:", sp.expand(mu_NY - mu_paper) == 0, flush=True)
assert sp.expand(mu_NY - mu_paper) == 0

# (B) disk: |lambda-1/2|^2 - 1/4 = |lambda|^2 - Re(lambda) = mu_NY/2
disk_excess = sp.expand(abs2 - Re_lam)
print("(B) |lambda-1/2|^2 - 1/4 == mu/2:", sp.expand(2*disk_excess - mu_NY) == 0, flush=True)
print("    => sign(mu) = sign(|lambda-1/2|^2 - 1/4):  mu>0 <=> lambda outside closed disk", flush=True)
assert sp.expand(2*disk_excess - mu_NY) == 0

# (C) Delta_n for n=5,7 at the dangerous mode k=(n-1)/2  (theta = pi - pi/n, cos th = -cos(pi/n) = c_n)
for n in (5, 7):
    cn = -sp.cos(sp.pi/n)
    Delta = sp.nsimplify(4*cn**2 - 4*cn - 7)
    # bracket at dangerous mode: 2q c_n + 2p^2 - 3p + 2, quadratic in p, lead coeff 2>0
    bracket = sp.expand(2*(1-p)*cn + 2*p**2 - 3*p + 2)
    disc = sp.discriminant(bracket, p)
    print("(C) n=%d: c_n=%s  Delta_n=%s<0  bracket_disc=%s==Delta_n: %s" % (
        n, sp.simplify(cn), Delta, disc, sp.simplify(disc - Delta) == 0), flush=True)
    assert Delta < 0
    # bracket > 0 for all real p (lead 2>0, disc<0) => in particular p in (0,1); 1-c>0 since k!=0
    # mu_k = 2(1-c)*bracket > 0  (1-c = 1-cos th > 0 for k!=0)
    print("    bracket lead coeff=%d>0, disc<0 => bracket>0 for all p; 1-cos th>0 => mu_k>0" %
          sp.expand(bracket).coeff(p, 2), flush=True)
print("DONE-NY-SPECTRAL", flush=True)
