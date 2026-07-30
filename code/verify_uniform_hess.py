#!/usr/bin/env python3
"""
Carefully verify (B.3): uniform reduced-Hessian eigenvalues.
Compute the Hessian of P at uniform SYMBOLICALLY, restrict to sum=0, get eigenvalues,
compare to GPT formula  lambda_k = 2(1-cos theta_k)[2 q cos theta_k + 2 p^2 - 3 p + 2].
Also compare to the log-Hessian (x_i -> e^{y_i}, Hessian in y at y=0).
Discrepancy earlier: numerical x-Hessian on sum=0 gave ~0. Investigate.
"""
import sympy as sp
import numpy as np
n=7
p = sp.symbols('p', positive=True); q = 1-p
xs = sp.symbols('x0:7', positive=True)
P = sum(xs[i]/(p*xs[(i+1)%n] + q*xs[(i+2)%n]) for i in range(n))

# Hessian at uniform x_i=1
subsU = {xs[i]: sp.Integer(1) for i in range(n)}
H = sp.Matrix(n, n, lambda i,j: sp.diff(P, xs[i], xs[j]).subs(subsU))
print("Uniform Hessian H (symbolic, in p):")
sp.pprint(H)
# restrict to sum=0: project  H_red = B^T H B, B = [e_i - e_6, i=0..5]  (sum=0 basis)
B = sp.Matrix(n, 6, lambda i,j: (1 if i==j else 0) - (1 if i==6 else 0))
# careful: B[:,j] has 1 at row j, -1 at row 6
B = sp.zeros(n,6)
for j in range(6):
    B[j,j] = 1; B[6,j] = -1
Hred = sp.simplify(B.T * H * B)
print("\nHred = B^T H B (6x6) on sum=0 subspace:")
# eigenvalues symbolically may be heavy; evaluate numerically at p=0.25
for pp in [0.25, 0.1, 0.5]:
    Hred_n = np.array(Hred.subs(p, sp.Rational(pp).limit_denominator(100)).evalf(), dtype=float)
    ev = np.linalg.eigvalsh(Hred_n)
    # GPT formula
    qv=1-pp; form=[]
    for k in [1,2,3]:
        th=2*np.pi*k/7
        form.append(2*(1-np.cos(th))*(2*qv*np.cos(th)+2*pp**2-3*pp+2))
    form=sorted(np.repeat(form,2))
    print(f"  p={pp}: Hred eig={np.round(sorted(ev),5)}")
    print(f"         formula ={np.round(form,5)}")

# Also: is the x-Hessian at uniform actually degenerate? Check trace / a diagonal entry
print("\nH[0,0] at uniform:", sp.simplify(H[0,0]))
print("H[0,1] at uniform:", sp.simplify(H[0,1]))
print("H[1,2] at uniform:", sp.simplify(H[1,2]))
