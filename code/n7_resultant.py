#!/usr/bin/env python3
"""
Verify GPT's n=7 hand-reduction and compute the EXACT resultant E(p)=Res_t(R,B).
R(p,t) = q^3 - p^3 t^5 - p^2 q t^8        (stationary curve, q=1-p)
B(p,t) = 5 p^2 t + 2 p q t^4 - 2 q^2 - 7 q p^2   (boundary P=7)
a_7,b_7 = p-coords of real solutions {R=0,B=0}.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 60

p, t = sp.symbols('p t')
q = 1 - p
R = q**3 - p**3*t**5 - p**2*q*t**8
B = 5*p**2*t + 2*p*q*t**4 - 2*q**2 - 7*q*p**2

# ---- 1. verify reduction against the high-precision nsolve boundary solution ----
# from phase3: at a_7, (a,b,c,d,e)=(x1,x2,x3,x4,x6) summing to 1:
a7_val = mp.mpf('0.21427352090984096558774231909001838797745454906426773562141322120851430')
b7_val = mp.mpf('0.32862767791659196669734085647711086500878730337285459186534900232409311442317')
# minimizer vars at a_7 (from phase3 nsolve, ~10 digits): re-derive at high precision via nsolve
import sympy as sp2
aa,bb,cc,dd,ee,lam,pv = sp2.symbols('aa bb cc dd ee lam pv', positive=True)
qv = 1 - pv
Pv = aa/(pv*bb+qv*cc)+bb/(pv*cc+qv*dd)+cc/(pv*dd)+dd/(qv*ee)+ee/(qv*aa)
gs=[sp2.diff(Pv,v) for v in (aa,bb,cc,dd,ee)]
eqs=[g-lam for g in gs]+[aa+bb+cc+dd+ee-1, Pv-7]
sol=sp2.nsolve(eqs,[aa,bb,cc,dd,ee,lam,pv],[0.17,0.23,0.06,0.31,0.23,6.0,0.22],prec=60,tol=mp.mpf('1e-50'),maxsteps=300)
a_n,b_n,c_n,d_n,e_n = [mp.mpf(sp2.N(x,55)) for x in sol[:5]]
p_n = mp.mpf(sp2.N(sol[6],55))
print(f"nsolve a_7 branch: p={mp.nstr(p_n,40)}")
print(f"  a,b,c,d,e = {mp.nstr(a_n,12)},{mp.nstr(b_n,12)},{mp.nstr(c_n,12)},{mp.nstr(d_n,12)},{mp.nstr(e_n,12)}")
# normalize a=1, t=e/a
t_val = e_n / a_n
d_pred = t_val**2
print(f"  t=e/a={mp.nstr(t_val,15)}  d={mp.nstr(d_n/a_n,15)}  t^2={mp.nstr(d_pred,15)}  (should match)")
# evaluate R,B at (p_n, t_val)
def Rf(pv,tv): qv=1-pv; return qv**3 - pv**3*tv**5 - pv**2*qv*tv**8
def Bf(pv,tv): qv=1-pv; return 5*pv**2*tv + 2*pv*qv*tv**4 - 2*qv**2 - 7*qv*pv**2
print(f"  R(p,t) = {mp.nstr(Rf(p_n,t_val),8)}  (should be ~0)")
print(f"  B(p,t) = {mp.nstr(Bf(p_n,t_val),8)}  (should be ~0)")

# ---- 2. exact resultant, via exact quotient/remainder ----
# E(p) = Res_t(R, B).  Claim (paper): E = p^15 * (p-1)^6 * F(p), F degree 15,
# and the two real (0,1)-roots of F are exactly a_7, b_7.
# Certified here by EXACT polynomial division (quotient + zero remainder), not by
# numerical factor evaluation.
import sys
print("\n=== E(p) = Res_t(R, B) via exact quotient/remainder ===")
E = sp.expand(sp.resultant(R, B, t))
deg_E = sp.degree(E, p)
print(f"deg_p E = {deg_E}")

# exact division: E / p^15, remainder must be 0
q1, r1 = sp.div(E, p**15, p)
assert sp.expand(r1) == 0, "p^15 does not divide E exactly"
print(f"E / p^15 : exact (remainder 0); deg_p quotient = {sp.degree(q1, p)}")

# exact division: q1 / (p-1)^6, remainder must be 0
q2, r2 = sp.div(q1, (p - 1)**6, p)
assert sp.expand(r2) == 0, "(p-1)^6 does not divide E/p^15 exactly"
F = sp.expand(q2)
deg_F = sp.degree(F, p)
print(f"(E/p^15) / (p-1)^6 : exact (remainder 0); deg_p F = {deg_F}")
assert deg_F == 15, "F is not degree 15"
print(f"\nF(p) = {F}")

# F squarefree (15 distinct roots): gcd(F, F') = 1
Fp = sp.diff(F, p)
g = sp.gcd(F, Fp)
sqfree = (sp.expand(g) == 1)
print(f"\ngcd(F, F') = 1 (F squarefree, 15 distinct roots): {sqfree}")
assert sqfree

# F irreducible over Q via finite-field certificate: F mod 23 is irreducible.
# Leading coeff 7^8 mod 23 = 12 != 0, so the degree is preserved mod 23; by Gauss
# irreducibility over F_23 lifts to irreducibility over Q.
lc_mod23 = sp.Mod(sp.LC(F, p), 23)
Fp23 = sp.Poly(F, p, modulus=23)
irred_mod23 = Fp23.is_irreducible
print(f"leading coeff 7^8 mod 23 = {lc_mod23} (!= 0, degree preserved)")
print(f"F mod 23 irreducible: {irred_mod23}  => F irreducible over Q (Gauss)")
assert irred_mod23, "F mod 23 is not irreducible"

# exact Sturm count of F's real roots in (0,1): must be 2 (= a_7, b_7)
FP = sp.Poly(F, p)
n_01 = FP.count_roots(0, 1)
# exclude endpoints: F(0) = -65536 != 0, F(1) = ?
F_at_0 = FP.eval(0); F_at_1 = FP.eval(1)
n_open_01 = n_01 - (1 if F_at_0 == 0 else 0) - (1 if F_at_1 == 0 else 0)
print(f"F(0) = {F_at_0}  F(1) = {F_at_1}  (neither 0 -> no endpoint roots)")
print(f"Sturm count of F real roots in (0,1): {n_open_01}  (== 2 => exactly a_7, b_7)")
assert n_open_01 == 2, "F does not have exactly 2 real roots in (0,1)"

# exact algebraic isolation of the two (0,1)-roots via CRootOf; confirm they
# match the high-precision nsolve a_7, b_7 (sanity: isolating intervals contain them)
roots_01 = sorted(float(sp.N(r, 40)) for r in sp.real_roots(FP) if 0 < r < 1)
match_a = abs(roots_01[0] - float(a7_val)) < 1e-20
match_b = abs(roots_01[1] - float(b7_val)) < 1e-20
print(f"  isolated (0,1)-roots: {roots_01[0]:.20f}  {roots_01[1]:.20f}")
print(f"  nsolve a_7         : {float(a7_val):.20f}")
print(f"  nsolve b_7         : {float(b7_val):.20f}")
print(f"  match a_7: {match_a}   match b_7: {match_b}")
assert match_a and match_b

with open("code/n7_minpoly_factor.txt", "w") as fh:
    fh.write(str(F) + "\n")
print("\nALL EXACT CHECKS PASS: E = p^15 (p-1)^6 F (exact quotient), deg F=15, "
      "F squarefree, F mod 23 irreducible (=> irreducible over Q), exactly 2 real (0,1)-roots = a_7,b_7.")
print("DONE-RESULTANT")
sys.exit(0)
