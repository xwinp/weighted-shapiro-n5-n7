#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify NEW content from GPT's COMPLETE n=9 Task A+B reply.
- Task A: R9,7 (A.1), E9,7; lifts (r±,u±,b±,c±,d±,e±) → R9,7=0, E9,7≠0, P=9; sample r=2/3.
- Task B: (B.4) inversion; (B.6) ρ⁹; H(B,C) (B.5); D-,D+ (B.9,B.10);
  Res_B(H,D-)=C¹²(C-1)¹²(2C-1)⁹Q17, Res_B(H,D+)=...Q22; D-<0/D+>0 at p=0.4;
  ρ(C_D)→p_D>p0; (B.2)(B.3) at (A0,B0,C0)."""
import numpy as np, sympy as sp, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

r,u = sp.symbols('r u', positive=True)
B,C,A = sp.symbols('B C A', real=True)

# ---- Task A: R9,7, E9,7 ----
E97 = r**5*u**11 + r**4*u**16 - 2*r**4*u**7 - 2*r**3*u**12 + r**3*u**3 + r**2*u**8 - 2*r*u**4 + 1
R97 = (r**8*u**6 + r**7*u**11 + 2*r**5*u**3 + 6*r**4*u**8 + 2*r**3*u**13 + 4*r**3*u**4
       - r**2*u**18 + r**2 - 4*r*u**14 + r*u**5 - 2*u**10)
print("="*70); print("TASK A: R9,7 / E9,7 at the two positive lifts"); print("="*70)
lifts = {
 'r-': (0.07650147965082751680, 0.59807863514696202007, 0.02684655935602285087, 0.21187776230867383272, 1.65628022730008034404, 0.00122799090542040579),
 'r+': (0.76642830175721873100, 1.03108996087378196415, 0.47762861602421434072, 0.73013160168322879863, 0.91514389630561341642, 0.24822594458142154662),
}
# KKT reconstruction: c=u^3-rb; d=r/u^4-rc; D2=rd+e=(u^4+r^3 b)/u^8; e=D2-rd; D3=re+1
# M9,7 = (1+r)[ u^2 + (rb+c)/u^2 + (rc+d)/b + (rd+e)/c + (re+1)/d + r/e ]
def M97(rv,uv,bv,cv,dv,ev):
    return (1+rv)*(uv**2 + (rv*bv+cv)/uv**2 + (rv*cv+dv)/bv + (rv*dv+ev)/cv + (rv*ev+1)/dv + rv/ev)
for key,(rv,uv,bv,cv,dv,ev) in lifts.items():
    Rval = float(R97.subs([(r,sp.Rational(rv).limit(sp.oo,0) if False else rv),(u,uv)]))
    Eval = float(E97.subs([(r,rv),(u,uv)]))
    # check KKT recurrence consistency
    c_kkt = uv**3 - rv*bv
    d_kkt = rv/uv**4 - rv*cv
    D2 = (uv**4 + rv**3*bv)/uv**8
    e_kkt = D2 - rv*dv
    M = M97(rv,uv,bv,cv,dv,ev)
    print(f" {key}: R9,7={Rval:.3e} (≈0?), E9,7={Eval:.6f} (≠0?), M9,7={M:.8f} (expect 9)")
    print(f"      KKT chk: c_rec={c_kkt:.6f} vs c={cv:.6f}; d_rec={d_kkt:.6f} vs d={dv:.6f}; e_rec={e_kkt:.6f} vs e={ev:.6f}")

# sample r=2/3: need a positive lift u on R9,7; find numerically
print("\n sample r=2/3: find u with R9,7=0, compute M9,7 (expect 8.8237..8.8239)")
rv=sp.Rational(2,3)
Rr = R97.subs(r, rv)
uroots = [float(x) for x in sp.nroots(Rr) if abs(sp.im(x))<1e-6 and sp.re(x)>0]
for uv in uroots:
    # recover b,c,d,e from KKT? need b. b is free param — actually need full lift.
    pass
print(f"  R9,7(2/3,u) positive roots (u): {len(uroots)} -> {[round(x,5) for x in uroots]}")

# ---- Task B ----
print("\n"+"="*70); print("TASK B: (B.4) inversion, (B.6) ρ⁹, H/D-/D+, resultants"); print("="*70)
HBC = (4*B**5*C**2-4*B**5*C+B**5+4*B**4*C**3-12*B**4*C**2+9*B**4*C-2*B**4
       -4*B**3*C**4+2*B**3*C**3+6*B**3*C**2-5*B**3*C+B**3-4*B**2*C**5+12*B**2*C**4
       -11*B**2*C**3+3*B**2*C**2+2*B*C**5-5*B*C**4+4*B*C**3-B*C**2
       +C**7-3*C**6+3*C**5-C**4)
Dm = (24*B**5*C**3-40*B**5*C**2+22*B**5*C-4*B**5+24*B**4*C**4-84*B**4*C**3+98*B**4*C**2
      -47*B**4*C+8*B**4-12*B**3*C**5+10*B**3*C**4+30*B**3*C**3-48*B**3*C**2+24*B**3*C-4*B**3
      -8*B**2*C**6+28*B**2*C**5-36*B**2*C**4+20*B**2*C**3-4*B**2*C**2-2*C**8+8*C**7-12*C**6+8*C**5-2*C**4)
Dp = (48*B**5*C**4-88*B**5*C**3+68*B**5*C**2-26*B**5*C+4*B**5+36*B**4*C**5-140*B**4*C**4
      +205*B**4*C**3-150*B**4*C**2+55*B**4*C-8*B**4-24*B**3*C**6+48*B**3*C**5-10*B**3*C**4
      -52*B**3*C**3+60*B**3*C**2-26*B**3*C+4*B**3-12*B**2*C**7+46*B**2*C**6-72*B**2*C**5
      +58*B**2*C**4-24*B**2*C**3+4*B**2*C**2-3*C**9+16*C**8-34*C**7+36*C**6-19*C**5+4*C**4)

def A_correct(bv,cv):  # correct (reciprocal of B.4), from (B.3)
    return 1-bv + bv**2*(1-bv)*(2*cv-1)/(cv**2*(1-cv))
def A_gpt(bv,cv):      # GPT's (B.4) as written
    return 1-bv + cv**2*(1-cv)/(bv**2*(1-bv)*(2*cv-1))

# (B.2),(B.3) check at (A0,B0,C0) and (A_D,B_D,C_D)
for name,(Av,Bv,Cv) in [('P=9 cross',(0.12685437921712335662,0.93513808432609586522,0.57695524131763708970)),
                        ('Hess deg ' ,(0.17954511354954542726,0.91142492877540238805,0.58806302678090147401))]:
    e2 = Bv**2*(1-Bv) - Av*(1-Av)*(Bv+Cv-1)
    e3 = Cv**2*(1-Cv)*(Av+Bv-1) - Bv**2*(1-Bv)*(2*Cv-1)
    Arec = A_correct(Bv,Cv)
    print(f" {name}: (B.2)res={e2:.3e} (B.3)res={e3:.3e}  A_given={Av:.6f} A_correct(recip)={Arec:.6f} A_gpt(B.4)={A_gpt(Bv,Cv):.4f}")

# (B.6) ρ⁹ at p=0.4 saddle: A=0.017256,B=0.9904,C=0.5574 → ρ=2/3
Av,Bv,Cv = 0.017256,0.9904,0.5574
rho9 = (2*Cv-1)*(Av+Bv-1)**4*(Bv+Cv-1)**2/(Av**2*(1-Av)**4*(1-Bv)**2*(1-Cv)**2)
print(f"\n (B.6) ρ⁹ at p=0.4 saddle = {rho9:.5e} → ρ={rho9**(1/9):.5f} (expect 2/3={2/3:.5f})")

# D-, D+ at p=0.4 saddle (B,C)
Dm_v = float(Dm.subs([(B,Bv),(C,Cv)])); Dp_v = float(Dp.subs([(B,Bv),(C,Cv)]))
print(f" p=0.4: D-={Dm_v:.4e} (<0?), D+={Dp_v:.4e} (>0?)  → det H ∝ D-·D+ sign={Dm_v*Dp_v:+.1e}")

# ρ(C_D) via (B.6) → p_D
CD,BD,AD = 0.58806302678090147401,0.91142492877540238805,0.17954511354954542726
rho9_D = (2*CD-1)*(AD+BD-1)**4*(BD+CD-1)**2/(AD**2*(1-AD)**4*(1-BD)**2*(1-CD)**2)
rho_D = rho9_D**(1/9); p_D = rho_D/(1+rho_D)
print(f" C_D→ρ={rho_D:.6f}→p_D={p_D:.10f} (GPT 0.4483521955); p0=0.4318363763 < p_D? {0.4318363763<p_D}")

# Resultants Res_B(H,D-), Res_B(H,D+) → Q17, Q22
print("\n computing Res_B(H,D-) and Res_B(H,D+) (heavy)...")
import signal
def handler(signum,frame): raise TimeoutError()
try:
    Rm = sp.resultant(HBC, Dm, B)
    Rp = sp.resultant(HBC, Dp, B)
    Rm = sp.expand(Rm); Rp = sp.expand(Rp)
    # factor out C^12(C-1)^12(2C-1)^9
    def strip(poly):
        p = sp.Poly(poly, C)
        # divide by C^12
        for _ in range(20):
            if p.eval(0)==0: p = sp.quo(p, sp.Poly(C,C), C)
            else: break
        # divide by (C-1)^12 and (2C-1)^9 by repeated substitution check
        return p
    Q17_gpt = sp.Poly(19616*C**17-82144*C**16+104080*C**15+111792*C**14-635984*C**13+1211840*C**12-1415716*C**11+1057038*C**10-375247*C**9-197120*C**8+418348*C**7-353452*C**6+196768*C**5-78254*C**4+22344*C**3-4408*C**2+544*C-32, C)
    Q22_gpt = sp.Poly(681858*C**22-5358069*C**21+19862874*C**20-45061200*C**19+66655008*C**18-58294348*C**17+3492304*C**16+82970486*C**15-158508822*C**14+181624473*C**13-141797994*C**12+65248658*C**11+6629860*C**10-45956452*C**9+50826852*C**8-36600888*C**7+19649040*C**6-8123536*C**5+2581588*C**4-614956*C**3+104016*C**2-11184*C+576, C)
    # check divisibility: Res / [C^12(C-1)^12(2C-1)^9 Q17] remainder zero?
    pref = C**12*(C-1)**12*(2*C-1)**9
    qm, rm = sp.div(sp.Poly(Rm,C), sp.Poly(pref*Q17_gpt.as_expr(), C), C)
    qp, rp = sp.div(sp.Poly(Rp,C), sp.Poly(pref*Q22_gpt.as_expr(), C), C)
    print(f" Res_B(H,D-) deg={sp.degree(Rm,C)}; = pref·Q17·(quot deg {qm.degree()})? remainder zero: {rm.is_zero}")
    print(f" Res_B(H,D+) deg={sp.degree(Rp,C)}; = pref·Q22·(quot deg {qp.degree()})? remainder zero: {rp.is_zero}")
    # Sturm
    print(f" Q17 roots in (1/2,1): {Q17_gpt.count_roots(sp.Rational(1,2),1)} (GPT 1)")
    print(f" Q22 roots in (1/2,1): {Q22_gpt.count_roots(sp.Rational(1,2),1)} (GPT 0)")
except Exception as e:
    print(f" resultant computation issue: {e}")
