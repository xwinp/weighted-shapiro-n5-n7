#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check monotonicity of rho(z) on the H_B positive-w branch, and set up a
rigorous cert.  rho(z)^7 = F(z) = w_+(z) z^5 / ((1-z) D(z)^3),
  w_+(z) = (-(1-z^2)+sqrt(disc))/(2z),  disc = z^4-4z^3+2z^2+1,  D=1-z^2+w z^2.
Branch is a graph over p  <=>  F'(z) has no zero on (0,1)  <=>  rho strictly monotone.
Numerically scan F'(z); if always one sign, build the polynomial G(z)=A^2-B^2*disc
where F'=A+B*sqrt(disc), Sturm-count its roots in (0,1), and sign-check.
"""
import numpy as np
import sympy as sp

z=sp.symbols('z',positive=True)
disc = z**4 - 4*z**3 + 2*z**2 + 1
w = (-(1-z**2)+sp.sqrt(disc))/(2*z)
D = 1 - z**2 + w*z**2
F = w*z**5/((1-z)*D**3)   # = rho^7
logF = sp.log(w) + 5*sp.log(z) - sp.log(1-z) - 3*sp.log(D)
dlogF = sp.diff(logF, z)
dlogF_s = sp.simplify(dlogF)
# Numerical scan
f = sp.lambdify(z, dlogF_s, 'numpy')
zs = np.linspace(0.01,0.99,500)
vals = f(zs)
vals = np.where(np.isfinite(vals), vals, np.nan)
print("d(log F)/dz numerical scan on (0.01,0.99):")
print("  min=%.6f  max=%.6f  all>0? %s  all<0? %s"%(np.nanmin(vals),np.nanmax(vals),np.nanmin(vals)>0,np.nanmax(vals)<0))
# also F itself monotone?
Ff = sp.lambdify(z, sp.simplify(F), 'numpy')
Fvals=Ff(zs)
print("  F(z): F(0.01)=%.4f F(0.5)=%.4f F(0.99)=%.4f (rho=F^(1/7))"%(Ff(0.01),Ff(0.5),Ff(0.99)))

# Build F' = A + B sqrt(disc).  Compute F' symbolically, substitute sqrt(disc)->S.
Fp = sp.diff(F, z)
Fp_together = sp.together(sp.expand(Fp))
num = sp.expand(sp.together(Fp).as_numer_denom()[0])
den = sp.expand(sp.together(Fp).as_numer_denom()[1])
S = sp.symbols('S')
numS = sp.expand(num.subs(sp.sqrt(disc), S))
numS = sp.collect(numS, S)
A = sp.expand(numS.coeff(S,0))
B = sp.expand(numS.coeff(S,1))
print("\nF' numerator = A + B*sqrt(disc); deg_z A=%d deg_z B=%d"%(sp.Poly(A,z).degree(), sp.Poly(B,z).degree() if B!=0 else 0))
# F'=0 (with den!=0) => A + B sqrt(disc)=0.  If B=0 need A=0.  General: A^2 - B^2*disc=0.
G = sp.expand(A**2 - B**2*disc)
GPoly = sp.Poly(G, z)
print("G = A^2 - B^2*disc, deg_z =", GPoly.degree())
rts = [r for r in sp.real_roots(GPoly) if 0<r<1]
print("G roots in (0,1):", [sp.N(r,12) for r in rts])
# For each root, check it's a REAL F'=0 root (sign of A and -B*sqrt(disc) match) and disc>0
for r in rts:
    Av=float(A.subs(z,r)); Bv=float(B.subs(z,r)); dv=float(disc.subs(z,r))
    if dv<0: print("  r=%.6f: disc<0, not on branch"%r); continue
    sq=dv**0.5
    fpv = Av + Bv*sq
    print("  r=%.6f: A=%+.6e B=%+.6e sqrt(disc)=%.6f  F'num=%+.6e  (real zero? %s)"%(r,Av,Bv,sq,fpv,abs(fpv)<1e-6))
print("DONE")
