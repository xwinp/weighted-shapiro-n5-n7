#!/usr/bin/env python3
"""
Verify (B.3): the uniform-point Hessian of P restricted to sum x_i = 0 is
diagonalised by the Fourier modes of C_7, with eigenvalue
    lambda_k = 2(1-cos theta_k)[2 q cos theta_k + 2 p^2 - 3 p + 2],
    theta_k = 2 pi k/7,  k = 1,2,3   (each of multiplicity 2).

IMPORTANT: the constraint basis B = [e_i - e_6]_{i=0..5} is NOT orthonormal
(B^T B != I), so the ORDINARY eigenvalues of B^T H B do NOT equal the spectrum
of the restricted operator. The correct computation is the GENERALISED
eigenproblem
    B^T H B v = lambda (B^T B) v                         (Rayleigh-Ritz, mass M = B^T B)
equivalently orthonormalising B first. The earlier script compared
eigvalsh(B^T H B) to the Fourier formula and (correctly) found a mismatch --
that mismatch was the bug, not a defect in the theorem: with the generalised
eigenproblem the two agree to machine precision.

All lambda_k > 0 on (0,1) (the bracket 2q cos theta_k + 2p^2 - 3p + 2 has
negative discriminant for the worst mode k=3, cos(6pi/7)) => the uniform point
is a strict local minimum on sum=0 (Nowosad-Yamagami interior-disk, symbol-dual).
"""
import sys
import sympy as sp
import numpy as np
from scipy.linalg import eigh as gen_eigh

n = 7
p = sp.symbols('p', positive=True); q = 1 - p
xs = sp.symbols('x0:7', positive=True)
P = sum(xs[i]/(p*xs[(i+1)%n] + q*xs[(i+2)%n]) for i in range(n))

subsU = {xs[i]: sp.Integer(1) for i in range(n)}
H = sp.Matrix(n, n, lambda i, j: sp.diff(P, xs[i], xs[j]).subs(subsU))
# B = [e_i - e_6, i=0..5]: a (non-orthonormal) basis of sum=0
B = sp.zeros(n, 6)
for j in range(6):
    B[j, j] = 1; B[6, j] = -1
Hred = sp.simplify(B.T * H * B)      # 6x6  = B^T H B
BtB  = sp.simplify(B.T * B)          # mass matrix M = B^T B  (6x6)

def formula_eig(pp):
    qv = 1 - pp; form = []
    for k in [1, 2, 3]:
        th = 2 * np.pi * k / 7
        form.append(2 * (1 - np.cos(th)) * (2 * qv * np.cos(th) + 2 * pp**2 - 3 * pp + 2))
    return sorted(np.repeat(form, 2))

print("Uniform reduced Hessian on sum=0: generalised eigenproblem B^T H B v = lam (B^T B) v")
ok = True
for pp in [0.25, 0.1, 0.5, 0.35, 0.75]:
    Hn = np.array(Hred.subs(p, sp.Rational(pp).limit_denominator(100)).evalf(), dtype=float)
    Mn = np.array(BtB.subs(p, sp.Rational(pp).limit_denominator(100)).evalf(), dtype=float)
    ev = sorted(gen_eigh(Hn, Mn)[0])            # generalised eigenvalues
    form = formula_eig(pp)
    match = np.allclose(ev, form, atol=1e-7)
    ok = ok and match
    print(f"  p={pp}: gen-eig ={np.round(ev, 5)}")
    print(f"         formula  ={np.round(form, 5)}  match={match}")

# strict local minimum: all lambda_k > 0  (check the worst mode k=3 symbolically)
k3 = 3; th3 = 2 * sp.pi * k3 / 7
bracket = 2 * q * sp.cos(th3) + 2 * p**2 - 3 * p + 2          # most negative cos -> smallest lambda
disc = sp.discriminant(bracket, p)                            # should be < 0 => no real roots, always > 0
lam_pos = bool(sp.simplify(disc) < 0) and bool(sp.limit(bracket, p, sp.Rational(1, 2)) > 0)
print(f"\nWorst-mode bracket (k=3) discriminant = {disc}  (<0 => bracket always >0): {bool(sp.simplify(disc)<0)}")
print("All reduced-Hessian eigenvalues > 0 (uniform strict local min on sum=0):", lam_pos)
print("SPECTRUM MATCH:", ok)
sys.exit(0 if (ok and lam_pos) else 1)
