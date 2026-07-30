#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify GPT n=9 Task B (L=9 one-zero face, det H_red<0 when P<9)."""
import numpy as np
import sympy as sp
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

C, B = sp.symbols('C B', real=True)
def sc(poly, var, lo, hi): return sp.Poly(poly,var).count_roots(lo,hi)
def roots_in(poly, var, lo, hi):
    return sorted([float(x) for x in sp.nroots(poly) if abs(sp.im(x))<1e-6 and lo<sp.re(x)<hi])

print("="*72); print("TASK B: L=9 one-zero face — branch uniqueness + P=9/Hessian ordering"); print("="*72)

# C_0 = root of C^3-2C^2-C+1 in (1/2,1)
fC0 = C**3-2*C**2-C+1
print(f" C³-2C²-C+1 roots in (1/2,1): {sc(fC0,C,sp.Rational(1,2),1)}  value={roots_in(fC0,C,0.5,1)}  (GPT C_0=0.554958132087371)")

# Δ_14
D14 = (2806*C**14-11761*C**13+23719*C**12-30454*C**11+25754*C**10-11507*C**9-3109*C**8
       +10573*C**7-10461*C**6+6764*C**5-3192*C**4+1120*C**3-284*C**2+48*C-4)
print(f" Δ_14 roots in (1/2,1): {sc(D14,C,sp.Rational(1,2),1)}  (GPT 0)")

# M_22
M22 = (681858*C**22-5358069*C**21+19862874*C**20-45061200*C**19+66655008*C**18-58294348*C**17
       +3492304*C**16+82970486*C**15-158508822*C**14+181624473*C**13-141797994*C**12
       +65248658*C**11+6629860*C**10-45956452*C**9+50826852*C**8-36600888*C**7+19649040*C**6
       -8123536*C**5+2581588*C**4-614956*C**3+104016*C**2-11184*C+576)
print(f" M_22 roots in (1/2,1): {sc(M22,C,sp.Rational(1,2),1)}  (GPT 0)  → ρ,p monotone")

# D_17
D17 = (19616*C**17-82144*C**16+104080*C**15+111792*C**14-635984*C**13+1211840*C**12
       -1415716*C**11+1057038*C**10-375247*C**9-197120*C**8+418348*C**7-353452*C**6
       +196768*C**5-78254*C**4+22344*C**3-4408*C**2+544*C-32)
nd17=sc(D17,C,sp.Rational(1,2),1); rd17=roots_in(D17,C,0.5,1)
print(f" D_17 roots in (1/2,1): {nd17}  (GPT 1)  C_D={rd17}  (GPT 0.588063026780902)")

# H(B,C) curve (2.6)
HBC = (4*B**5*C**2-4*B**5*C+B**5+4*B**4*C**3-12*B**4*C**2+9*B**4*C-2*B**4
       -4*B**3*C**4+2*B**3*C**3+6*B**3*C**2-5*B**3*C+B**3-4*B**2*C**5+12*B**2*C**4
       -11*B**2*C**3+3*B**2*C**2+2*B*C**5-5*B*C**4+4*B*C**3-B*C**2
       +C**7-3*C**6+3*C**5-C**4)

def fixed_branch_at_C(cv):
    """Given C, find legit B root of H(B,C)=0 in (0,1), then A from (2.5), rho from (2.8)."""
    Hp = HBC.subs(C, cv)
    Brts = [float(x) for x in sp.nroots(Hp) if abs(sp.im(x))<1e-6 and 0<sp.re(x)<1]
    cands=[]
    for bv in Brts:
        # GPT (2.5) A=C^2(C-1)/[(B-1)(...)] is INVERTED (gives A>1).
        # Correct A = reciprocal, verified as common root of (2.3)&(2.4). A in (0,1).
        den=(bv-1)*(2*bv**2*cv-bv**2-cv**3+cv**2)
        if abs(den)<1e-12: continue
        Av=den/(cv**2*(cv-1))   # reciprocal of GPT's (2.5)
        if not (0<Av<1): continue
        if not (Av+bv>1 and bv+cv>1): continue
        # rho^9 = A^2(1-A)^4(1-B)^2(1-C)^2(2C-1)/[(A+B-1)^4(B+C-1)^2]
        num=Av**2*(1-Av)**4*(1-bv)**2*(1-cv)**2*(2*cv-1)
        dd=(Av+bv-1)**4*(bv+cv-1)**2
        if dd<=0 or num<=0: continue
        rho9=num/dd; rho=rho9**(1/9)
        # P via (2.9) and via 2(1+A+B+C) consistency
        P29=rho**2/(2*(1+rho))*Av*(1-Av)/(Av+bv-1)*(1+Av+bv+cv)
        Psum=2*(1+Av+bv+cv)
        cands.append((cv,bv,Av,rho,P29,Psum))
    return cands

print("\n -- trace L=9 fixed branch, find C_P (P=9 crossing) --")
C0=roots_in(fC0,C,0.5,1)[0]
Cs=np.linspace(C0+0.002, 0.999, 300)
branch=[]
for cv in Cs:
    cs=fixed_branch_at_C(cv)
    if cs:
        branch.append(cs[0])
print(f"  branch points found: {len(branch)} (C from {branch[0][0]:.4f} to {branch[-1][0]:.4f})")
# check P29 vs Psum consistency
maxdiff=max(abs(b[4]-b[5]) for b in branch)
print(f"  max|P(2.9) - 2(1+A+B+C)| = {maxdiff:.3e}  (consistency of GPT's P formula)")
# find C_P where P=9
Pvals=[b[4] for b in branch]
Cp=None
for i in range(len(branch)-1):
    if (Pvals[i]-9)*(Pvals[i+1]-9)<0:
        # linear interp
        t=(9-Pvals[i])/(Pvals[i+1]-Pvals[i])
        Cp=branch[i][0]+t*(branch[i+1][0]-branch[i][0])
        rho_p=branch[i][3]+t*(branch[i+1][3]-branch[i][3])
        pP=1/(1+rho_p)
        break
if Cp is not None:
    print(f"  C_P (P=9 crossing) = {Cp:.16f}  (GPT 0.5769552413176371)")
    print(f"  p_P = {pP:.16f}  (GPT 0.4318363763332235)")
    print(f"  p_P < b_9=0.4338858820? {pP < 0.4338858820}")
    print(f"  C_P < C_D=0.5880630268? {Cp < 0.5880630268}")
else:
    print("  C_P not found in scan range!")

# p=0.4 stationary: rho=1.5, find C with rho=1.5
print("\n -- L=9 fixed branch at p=0.4 (rho=1.5), expect P=8.829694, det<0 --")
target_rho=1.5
rhos=[b[3] for b in branch]
s40=None
for i in range(len(branch)-1):
    if (rhos[i]-target_rho)*(rhos[i+1]-target_rho)<0:
        t=(target_rho-rhos[i])/(rhos[i+1]-rhos[i])
        s40=tuple(branch[i][k]+t*(branch[i+1][k]-branch[i][k]) for k in range(6))
        break
if s40:
    cv,bv,Av,rho,P29,Psum=s40
    print(f"  C={cv:.8f} B={bv:.8f} A={Av:.8f} rho={rho:.8f}")
    print(f"  P_9,9^stat = {P29:.9f}  (GPT 8.829694260862)  <9? {P29<9}")
    # reconstruct x from term values T=(1,A,B,C,C,B,A,1) for x_1..x_8, solve x_i=T_i(p x_{i+1}+q x_{i+2})
    p=0.4; q=0.6; T=[1,Av,bv,cv,cv,bv,Av,1]
    # indices 1..8 (0-based 0..7), x0=0, x9=x0=0, x10=x1
    M=np.zeros((8,8))
    for i in range(8):
        M[i,i]+=1
        j1=(i+1)  # x_{i+2} in 1-based = index i+1 in 0-based for x_1.. ; careful
    # x_i (1-based i=1..8) = T_{i} (p x_{i+1} + q x_{i+2}); 0-based idx = i-1
    M=np.zeros((8,8))
    for i in range(1,9):
        ii=i-1
        M[ii,ii]+=1
        # x_{i+1}: 1-based i+1; if i+1==9 -> x_9=x_0=0 (skip); else 0-based (i+1)-1=i
        if i+1<=8:
            M[ii, (i+1)-1] += -T[ii]*p
        # x_{i+2}: 1-based i+2; if i+2==9 -> x_0=0; if i+2==10 -> x_1 (0-based 0)
        if i+2<=8:
            M[ii, (i+2)-1] += -T[ii]*q
        elif i+2==10:
            M[ii, 0] += -T[ii]*q
    # null space
    U,S,Vt=np.linalg.svd(M)
    x=Vt[-1]
    x=np.abs(x); x=x/x[0]  # normalize x_1=1
    xn=np.concatenate(([0],x))  # x_0..x_8
    def Pval(xx,p):
        q=1-p; s=0.0
        for i in range(9):
            d=p*xx[(i+1)%9]+q*xx[(i+2)%9]
            if d<=0: return 1e18
            s+=xx[i]/d
        return s
    Pdirect=Pval(xn,p)
    print(f"  reconstructed x P (direct) = {Pdirect:.9f}  (should match P_9,9^stat)")
    # numerical Hessian on {0} face (7 free vars x_1..x_8, fix x_1=1 gauge -> 7 vars? use x_2..x_8 free, x_1=1)
    # reduced Hessian wrt x_2..x_8 (7 vars), x_0=0,x_1=1
    def P7(v):
        xx=np.array([0.0,1.0]+list(v))
        return Pval(xx,p)
    v0=x[1:]  # x_2..x_8
    h=1e-4
    H=np.zeros((7,7))
    for i in range(7):
        for j in range(7):
            ei=np.zeros(7); ej=np.zeros(7); ei[i]=h; ej[j]=h
            if i==j:
                H[i,j]=(P7(v0+ei)-2*P7(v0)+P7(v0-ei))/(h*h)
            else:
                H[i,j]=(P7(v0+ei+ej)-P7(v0+ei-ej)-P7(v0-ei+ej)+P7(v0-ei-ej))/(4*h*h)
    ev=np.linalg.eigvalsh(H)
    print(f"  reduced Hessian (7x7) eigenvalues: {np.sort(ev)}")
    print(f"  det(H_red) = {np.linalg.det(H):.6e}  (GPT: <0 at p=0.4)")
    print(f"  neg eigenvalue count (Morse index): {int((ev<-1e-4).sum())}  (GPT saddle)")

# double-min lemma identity: X_T(y)=X_{T^T}(y^{-1})
print("\n -- double-min lemma identity X_T(y)=X_{T^T}(y^{-1}) (formal) --")
import sympy as sp
y=sp.IndexedBase('y'); n_=sp.symbols('n',positive=True,integer=True)
# verify numerically on a random T (3x3) that sum y_i/(Ty)_i = sum (1/y_i)/((T^T)(1/y))_i
Tm=np.array([[2,1,0],[0,2,1],[1,0,2]],dtype=float)
yy=np.array([1.3,0.7,2.1])
lhs=sum(yy[i]/(Tm@yy)[i] for i in range(3))
yinv=1/yy
rhs=sum(yinv[i]/(Tm.T@yinv)[i] for i in range(3))
print(f"  random 3x3 T: X_T(y)={lhs:.10f}, X_{{T^T}}(y^-1)={rhs:.10f}, match? {abs(lhs-rhs)<1e-10}")
# stationary dual x†=(Sx)^2/x : verify gradient condition implies it
print("  (x†=(Sx)²/x dual stationary — formal consequence of ∇f_S=0; GPT uses Yamagami Thm1, already verified in n=7 closure)")
