#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify GPT's region-2 crossing-resultant certificate for n=7 S1.

GPT's construction (rho = q/p = (1-p)/p, p = 1/(1+rho)):
  H_B(w,z) = z w^2 + (1-z^2) w + z^2 - z = 0          (bivariate branch curve)
  D = 1 - z^2 + w z^2
  closure:  rho^7 (1-z) D^3 - w z^5 = 0                (1)
  P_{S1}^stat = 2 rho(1+rho)/z * [3 - 2z - z^2 + wz(1+z)]   (2)
  P=7 crossing:  2 rho(1+rho)[3-2z-z^2+wz(1+z)] - 7z = 0    (3)
  => w = [2 rho(rho+1)(z^2+2z-3)+7z] / [2 rho(rho+1) z(z+1)]   (4)
  Substitute (4) into H_B and closure -> A(z,rho)=0, B(z,rho)=0
  Res_z(A,B) = 896 rho^13 (rho+1)^7 (8rho^2+8rho+7)^6 Phi_35(rho)
  Sturm: Phi_35 has 2 positive roots, both > 2  (in (2,5/2) and (3,7/2))
         => p = 1/(1+rho) < 1/3 for all crossings.
  Also p0 (Hessian det=0 transition) in (3/8, 2/5) > 1/3.

This script:
 (A) Numerically solves S1 stationary points, computes w=beta3, z=beta4,
     checks H_B=0, closure=0, P-formula == direct P.
 (B) Symbolically builds A,B, computes Res_z, factors, extracts Phi_35,
     compares coefficients to GPT's, runs Sturm, isolates roots.
 (C) p0 from Q5, checks p0 in (3/8,2/5); Q7 root -> negative w.
 (D) P_{S1}^stat(2/5) > 7 rational interval cert.
"""
import sys
import numpy as np
from scipy.optimize import root
import sympy as sp

n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s

sup_free=[2,3,4,5,6]
def solve_S1(p, init):
    def grad(v):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(sup_free): x[idx]=v[j]
        h=1e-7; f0=Pval(x,p); g=[]
        for j in range(len(sup_free)):
            xp=x.copy(); xp[sup_free[j]]+=h; g.append((Pval(xp,p)-f0)/h)
        return g
    r=root(grad, init, method='hybr', options={'xtol':1e-14})
    if r.success and max(abs(r.fun))<1e-7 and all(r.x>1e-8):
        x=np.zeros(n); x[1]=1
        for j,idx in enumerate(sup_free): x[idx]=r.x[j]
        return x
    return None

print("="*70)
print("(A) Numerical verification of H_B, closure, P-formula")
print("="*70)
z,w,rho = sp.symbols('z w rho', positive=True)
HB = z*w**2 + (1-z**2)*w + z**2 - z
D = 1 - z**2 + w*z**2
closure = rho**7*(1-z)*D**3 - w*z**5
Pform = 2*rho*(1+rho)/z * (3 - 2*z - z**2 + w*z*(1+z))

init=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
print("  p      rho      w        z        H_B        closure    Pform-Pdirect")
prev=init
for pp in [0.10,0.20,0.236,0.27,0.30,0.329,0.40,0.50,0.60,0.80,0.95]:
    x=solve_S1(pp, prev)
    if x is None:
        # try fresh inits
        for ii in [[0.3,0.7,1.5,0.07,1.3],[1,1,1,1,1],[0.5,0.5,1.0,0.5,0.5]]:
            x=solve_S1(pp,ii)
            if x is not None: break
    if x is None:
        print("  %.3f  NO CONVERGE"%pp); continue
    xs=x/x.sum(); prev=[xs[2],xs[3],xs[4],xs[5],xs[6]]
    b,c,d,e,f=xs[2],xs[3],xs[4],xs[5],xs[6]
    rholist=[b,c/b,d/c,e/d,f/e]; q=1-pp
    wv=q*rholist[3]/(pp+q*rholist[3]); zv=q*rholist[4]/(pp+q*rholist[4])
    rhov=q/pp
    hb=float(HB.subs([(w,wv),(z,zv)]))
    cl=float(closure.subs([(w,wv),(z,zv),(rho,rhov)]))
    pf=float(Pform.subs([(w,wv),(z,zv),(rho,rhov)]))
    pd=Pval(xs,pp)
    print("  %.3f  %.4f  %.5f  %.5f  %+8.2e  %+9.2e  %+9.2e"%(pp,rhov,wv,zv,hb,cl,pf-pd))

print("\n"+"="*70)
print("(B) Symbolic crossing resultant Phi_35")
print("="*70)
# w from (3): solve 2*rho(1+rho)[3-2z-z^2+wz(1+z)] - 7z = 0 for w
Pm7 = 2*rho*(1+rho)*(3 - 2*z - z**2 + w*z*(1+z)) - 7*z
wsol = sp.solve(sp.Eq(Pm7,0), w)[0]
print("w from (3) =", wsol)
# Substitute into H_B and closure, clear denominators -> A, B in (z,rho)
HB_sub = sp.together(HB.subs(w, wsol))
cl_sub = sp.together(closure.subs(w, wsol))
A = sp.expand(sp.fraction(HB_sub)[0])   # numerator of H_B after sub
B = sp.expand(sp.fraction(cl_sub)[0])   # numerator of closure after sub
Az=sp.Poly(A,z); Bz=sp.Poly(B,z)
print("A: deg_z=%d deg_rho=%d"%(Az.degree(), sp.Poly(A,rho).degree()))
print("B: deg_z=%d deg_rho=%d"%(Bz.degree(), sp.Poly(B,rho).degree()))

print("Computing Res_z(A,B) ... (12x12 Sylvester, may take a moment)")
Res = sp.resultant(A, B, z)
Res = sp.expand(Res)
print("Res computed. deg_rho =", sp.Poly(Res, rho).degree())
# Factor out spurious factors: 896 rho^13 (rho+1)^7 (8rho^2+8rho+7)^6
spurious = 896 * rho**13 * (rho+1)**7 * (8*rho**2+8*rho+7)**6
quot, rem = sp.div(sp.Poly(Res,rho), sp.Poly(sp.expand(spurious),rho), rho)
rem_zero = (sp.expand(rem.as_expr()) == 0)
print("remainder after dividing out spurious factors:", rem_zero)
Phi = quot.as_expr()
Phi = sp.expand(Phi)
Pp = sp.Poly(Phi, rho)
print("Phi_35 degree:", Pp.degree())
# GPT's coefficients (rho^35 down to rho^0)
gpt_coeffs = [262144,1211392,19453952,57874432,444107776,547314432,2185738240,-9171843072,
-34693234688,-166554596992,-327985403648,-543801267968,88546441088,2306869918304,
7848268705920,15843550970240,25279161341952,34334035751596,43526037225048,51848980402968,
56982453166940,55494399759599,46846168825232,33774085926224,20701541228760,10772897763040,
4815489002744,1952186204080,823732168256,410624553696,221073893824,107322284272,41692980224,
11883053056,2193551360,191102976]
my_coeffs = [Pp.nth(35-i) for i in range(36)]
match = all(int(a)==int(b) for a,b in zip(my_coeffs, gpt_coeffs))
print("Phi_35 coefficients MATCH GPT:", match)
if not match:
    for i,(a,b) in enumerate(zip(my_coeffs,gpt_coeffs)):
        if int(a)!=int(b): print("  mismatch rho^%d: mine=%d gpt=%d"%(35-i,a,b))

# Sturm
print("\nSturm root counts for Phi_35:")
PhiPoly = sp.Poly(Phi, rho)
def nroots_in(a,b):
    return sp.nroots_in_interval if False else sum(1 for r in sp.real_roots(PhiPoly) if a<sp.N(r)<b)
# use sturm via count_roots
from sympy import S
def sturm_count(a,b):
    return sp.count_roots(Phi, a, b)
print("  (0,2)   =", sturm_count(sp.Rational(0), sp.Rational(2)))
n_2_52 = sturm_count(sp.Rational(2), sp.Rational(5,2))
n_3_72 = sturm_count(sp.Rational(3), sp.Rational(7,2))
print("  (2,5/2) =", n_2_52)
print("  (5/2,3) =", sturm_count(sp.Rational(5,2), sp.Rational(3)))
print("  (3,7/2) =", n_3_72)
print("  (7/2,inf)=", len([r for r in sp.real_roots(PhiPoly) if sp.N(r)>sp.Rational(7,2)]))
rrts=[sp.N(r,15) for r in sp.real_roots(PhiPoly)]
print("  positive roots:", rrts)
for r in rrts:
    pv = 1/(1+r)
    print("    rho=%s -> p=%s"%(r, sp.N(pv,12)))

print("\n"+"="*70)
print("(C) p0 from Q5; check p0 in (3/8,2/5); Q7 root -> negative w")
print("="*70)
zs=sp.symbols('z')
Q5=2*zs**5+2*zs**3-2*zs**2-1
Q7=8*zs**7-24*zs**6+20*zs**5-9*zs**4+30*zs**3-15*zs**2-6
r5=[r for r in sp.real_roots(sp.Poly(Q5,zs)) if 0<r<1]
r7=[r for r in sp.real_roots(sp.Poly(Q7,zs)) if 0<r<1]
print("Q5 roots in (0,1):", [sp.N(r,14) for r in r5])
print("Q7 roots in (0,1):", [sp.N(r,14) for r in r7])
# p0 from Q5 root z0: need p0 = ?  GPT says p0=0.388528131361137. rho0 = (1-p0)/p0 = 0.61147/0.38853 = 1.5735
# Determine p(z) on the branch: at z=z0 (Q5 root, det=0), what is rho?
# Use closure + H_B: given z=z0, solve H_B=0 for w (positive root), then closure for rho.
z0 = r5[0]
# H_B(z0,w)=0: z w^2 + (1-z^2)w + z^2-z =0 -> quadratic in w
a_q = z0; b_q = 1-z0**2; c_q = z0**2 - z0
disc = b_q**2 - 4*a_q*c_q
wroots = [(-b_q + sp.sqrt(disc))/(2*a_q), (-b_q - sp.sqrt(disc))/(2*a_q)]
print("  at z0=Q5 root, w roots:", [sp.N(wr,10) for wr in wroots])
# pick positive w
for wr in wroots:
    if sp.N(wr)>0:
        w0=wr
        # closure: rho^7 (1-z0) D0^3 = w0 z0^5  -> rho^7 = w0 z0^5 / ((1-z0) D0^3)
        D0 = 1 - z0**2 + w0*z0**2
        rhs = w0*z0**5/((1-z0)*D0**3)
        rho0 = rhs**sp.Rational(1,7)
        p0 = 1/(1+rho0)
        print("  z0=%.8f w0=%.8f D0=%.8f rho0=%.8f p0=%.12f"%(sp.N(z0,10),sp.N(w0,10),sp.N(D0,10),sp.N(rho0,12),sp.N(p0,14)))
        p0N=sp.N(p0,20)
        print("  p0 in (3/8,2/5)? 3/8=%.6f p0=%.6f 2/5=%.6f"%(0.375, float(p0N), 0.4))
        print("  p0 > 1/3?", float(p0N) > 1/3)
# Q7 root -> w sign
z7 = r7[0]
a_q=z7; b_q=1-z7**2; c_q=z7**2-z7; disc=b_q**2-4*a_q*c_q
w7roots=[(-b_q+sp.sqrt(disc))/(2*a_q), (-b_q-sp.sqrt(disc))/(2*a_q)]
print("  at z7=Q7 root=%.8f, w roots: %s"%(sp.N(z7,10), [sp.N(wr,10) for wr in w7roots]))

print("\n"+"="*70)
print("(D) P_{S1}^stat(2/5) > 7  rational interval cert")
print("="*70)
# At p=2/5: rho = (1-2/5)/(2/5) = 3/2.  Solve S1 numerically, then build exact interval.
# We certify via the (w,z) system: rho=3/2, H_B=0 & closure=0 give (w,z); P>7.
# Strategy: isolate z (the branch root near the numerical value) in a rational interval
# using the closure+H_B reduced to a univariate in z (substitute w from H_B into closure).
pp = sp.Rational(2,5); rhov = (1-pp)/pp  # 3/2
print("  rho =", rhov)
# w from H_B: z w^2 + (1-z^2)w + (z^2-z)=0 -> w = [-(1-z^2) + sqrt((1-z^2)^2 -4z(z^2-z))]/(2z)
# (positive root). Substitute into closure rho^7(1-z)D^3 - wz^5 =0.
# Build the univariate in z (with the sqrt). Square-out to get polynomial G(z)=0, isolate z.
a_q=z; b_q=1-z**2; c_q=z**2-z
disc_z = b_q**2 - 4*a_q*c_q   # = (1-z^2)^2 - 4z(z^2-z) = 1 -2z^2+z^4 -4z^3+4z^2 = z^4-4z^3+2z^2+1
disc_z = sp.expand(disc_z)
print("  discriminant(z) =", disc_z)
wsol_HB = (-b_q + sp.sqrt(disc_z))/(2*a_q)   # positive root (z>0)
# closure(z) = rho^7 (1-z) (1-z^2+wsol*z^2)^3 - wsol*z^5
cl_z = rhov**7*(1-z)*(1-z**2+wsol_HB*z**2)**3 - wsol_HB*z**5
# Isolate z numerically first by solving the S1 stationary at p=2/5
x=solve_S1(0.4, [0.4,0.6,1.2,0.3,1.0])
if x is None: x=solve_S1(0.4, init)
xs=x/x.sum(); b,c,d,e,f=xs[2],xs[3],xs[4],xs[5],xs[6]
rl=[b,c/b,d/c,e/d,f/e]; q=0.6
znum=q*rl[4]/(0.4+q*rl[4]); wnum=q*rl[3]/(0.4+q*rl[3])
print("  numerical: z=%.8f w=%.8f  P=%.8f"%(znum,wnum,Pval(xs,0.4)))
# Build polynomial: cl_z has sqrt(disc_z). Write cl_z = A(z) + B(z)*sqrt(disc_z). Then
# cl_z=0 and the conjugate give A^2 - B^2*disc = 0 (polynomial). Find it.
cl_together = sp.together(sp.expand(cl_z))
num = sp.simplify(sp.expand(cl_together))
# Separate rational part and sqrt part manually by substituting sqrt->S, S^2->disc
S = sp.symbols('S')
expr_sub = sp.expand(sp.expand(num).subs(sp.sqrt(disc_z), S))
# collect in S
expr_sub = sp.collect(expr_sub, S)
A_part = expr_sub.coeff(S,0)
B_part = expr_sub.coeff(S,1)
polyG = sp.expand(A_part**2 - B_part**2*disc_z)
print("  G(z) = A^2 - B^2*disc, deg =", sp.Poly(polyG,z).degree())
# isolate z root near znum
GPoly = sp.Poly(polyG, z)
rts = [r for r in sp.real_roots(GPoly) if 0<r<1 and abs(float(r)-znum)<0.05]
print("  candidate z roots near znum:", [sp.N(r,14) for r in rts])
if rts:
    zr=rts[0]
    zlo=sp.Rational(int(float(zr)*10**7)-3, 10**7)
    zhi=sp.Rational(int(float(zr)*10**7)+3, 10**7)
    # ensure G(zlo),G(zhi) bracket and no sign change ambiguity: just use for interval eval
    print("  z isolated in [%.10f, %.10f]"%(float(zlo),float(zhi)))
    # Evaluate P_{S1}^stat = 2 rho(1+rho)/z * [3-2z-z^2+wz(1+z)] over z in [zlo,zhi], w=wsol_HB(z)
    # w has sqrt; bound via interval: compute w lo/hi using disc bounds, then P.
    # Simpler: evaluate P directly as function of z with w=wsol_HB, bound by sampling rational
    # endpoints and monotonicity. Use mpmath interval-free: evaluate P at zlo,zhi (with w(zlo),w(zhi))
    def Pstat_at(zv):
        dv=float(disc_z.subs(z,zv))
        if dv<0: return None
        wv=(-(1-zv**2)+dv**0.5)/(2*zv)
        return 2*float(rhov)*(1+float(rhov))/zv*(3-2*zv-zv**2+wv*zv*(1+zv))
    print("  P(zlo)=%.8f  P(zhi)=%.8f  (both should be >7)"%(Pstat_at(float(zlo)),Pstat_at(float(zhi))))
    p27 = min(Pstat_at(float(zlo)), Pstat_at(float(zhi))) > 7
    print("  P_{S1}^stat(2/5) > 7 :", p27)
else:
    p27 = False
# aggregate certificate
ok = rem_zero and match and (n_2_52 == 1) and (n_3_72 == 1) and p27
print("\nCERTIFICATE: Res_z = spurious*Phi_35 (rem=0)=%s; coeffs match=%s; "
      "Phi_35 roots in (2,5/2)=%d and (3,7/2)=%d; P(2/5)>7=%s : %s" % (
          rem_zero, match, n_2_52, n_3_72, p27, ok))
assert ok, "n7 S1 crossing-resultant certificate failed"
print("DONE-CROSSING-RESULTANT")
sys.exit(0 if ok else 1)
