#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-branch refined blow-up at the two exceptional tau = 1/2 and tau = 0.
At these, the leading B1 term -512 tau^2 (2 tau-1) eta^3 vanishes, so the
v~eta balance from the generic collar is degenerate. Resolve exactly.

U = 2F0 - H0  (transverse), vanishes at v=0, nu_v=1:
   U = v*Up + v^2*Up2 + v^3*Up3 + ...
   Up  = 16 eta^4 (eta-2)^4 R(1-eta,tau)            [known]
   Up2 = coeff v^2 of U.
Transverse equation (U/v = 0):  Up(eta,tau) + v*Up2(eta,tau) + ... = 0.
Longitudinal B = (F0-H0)/2 = B0 + v*B1 + v^2*B2 + ...
   B0 = -4 eta^5 (eta-2)^5
   B1 = coeff v^1 (eta-order 3 generically, leading -512 tau^2(2tau-1) eta^3)
At tau=1/2 and tau=0, expand Up, Up2, B0, B1 in eta and solve the 2x2 leading
system in (v, eta) for the Puiseux slope v ~ eta^alpha.
"""
import sympy as sp
from sympy import Poly, symbols, ZZ, factor, Rational
v,xi,tau,eta=symbols('v xi tau eta')
c,d,sg=symbols('c d sigma')
vv=c+d-1
Lc=vv+sg*c**2; Uc=vv-sg*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+sg*d**2; Ud=vv-sg*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - sg*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - sg*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(sp.expand(Fc),c,d,sg,domain=sp.QQ); Pd=Poly(sp.expand(Fd),c,d,sg,domain=sp.QQ)
FL=Pc.exquo(Poly((c-1)*Lc,c,d,sg,domain=sp.QQ)); FR=Pd.exquo(Poly((d-1)*Ld,c,d,sg,domain=sp.QQ))
Gm=Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,sg,domain=sp.QQ).exquo(Poly(c-d,c,d,sg,domain=sp.QQ))
Gp=Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,sg,domain=sp.QQ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gm_r=sp.together(Gm.as_expr().subs({c:cc,d:dd,sg:sig})); Gp_r=sp.together(Gp.as_expr().subs({c:cc,d:dd,sg:sig}))
Gmt=sp.expand(sp.fraction(Gm_r)[0]); Gpt=sp.expand(sp.fraction(Gp_r)[0])
F=Poly(Gmt,v,xi,tau,domain=sp.QQ); H=Poly(Gpt,v,xi,tau,domain=sp.QQ)

def v_order(P,var):
    m=None; idx=P.gens.index(var)
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if m is None or e<m: m=e
    return m if m is not None else 0
def vcoeff(P,var,k):
    """return Poly (in remaining gens) = coeff of var^k, dropping var from gens."""
    idx=P.gens.index(var)
    d={}
    for monom,coef in P.terms():
        if coef==0: continue
        if monom[idx]==k:
            nm=[monom[i] for i in range(len(monom)) if i!=idx]
            d[tuple(nm)]=d.get(tuple(nm),0)+coef
    gens=[g for i,g in enumerate(P.gens) if i!=idx]
    return Poly(d,*gens,domain=sp.QQ) if d else Poly(0,*gens,domain=sp.QQ)

nuF=v_order(F,v); nuH=v_order(H,v)
# Rebuild F0,H0 as polys in (v,xi,tau) by clearing v^nu (shift).
def shift_v(P,var,k):
    idx=P.gens.index(var); d={}
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if e>=k:
            nm=list(monom); nm[idx]=e-k; d[tuple(nm)]=d.get(tuple(nm),0)+coef
    return Poly(d,*P.gens,domain=sp.QQ)
F0=shift_v(F,v,nuF); H0=shift_v(H,v,nuH)
U=Poly(sp.expand(2*F0.as_expr()-H0.as_expr()),v,xi,tau,domain=sp.QQ)
B=Poly(sp.expand((F0.as_expr()-H0.as_expr())/2),v,xi,tau,domain=sp.QQ)
p=v_order(U,v)
Up =vcoeff(U,v,p)    # poly in (xi,tau) -- but vcoeff returns without v gen
# vcoeff returns Poly in (xi,tau). Good. But careful: U has gens (v,xi,tau); vcoeff drops v -> (xi,tau).
Up2=vcoeff(U,v,p+1)  # coeff of v^{p+1}
B0p=vcoeff(B,v,0)
B1p=vcoeff(B,v,1)
B2p=vcoeff(B,v,2)
print("nu_v(F)=",nuF,"nu_v(H)=",nuH," p=nu_v(U)=",p,flush=True)

# substitute xi=1-eta
def to_eta(Pxi_tau):
    return Poly(sp.expand(Pxi_tau.as_expr().subs({xi:1-eta})),eta,tau,domain=sp.QQ)
Up_e =to_eta(Up);  Up2_e=to_eta(Up2)
B0_e =to_eta(B0p); B1_e =to_eta(B1p); B2_e=to_eta(B2p)

def eta_terms(P,emax=8):
    """return dict eta_pow -> Poly in tau (coeff), for eta pow 0..emax."""
    idx=P.gens.index(eta); tidx=P.gens.index(tau); d={}
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if e<=emax:
            d.setdefault(e,{})
            d[e][monom[tidx]]=d[e].get(monom[tidx],0)+coef
    out={}
    for e in d:
        out[e]=Poly(d[e],tau,domain=sp.QQ)
    return out

for tval,lab in [(Rational(1,2),'tau=1/2'),(0,'tau=0')]:
    print("\n================ %s ================"%lab,flush=True)
    Up_t =Up_e.eval({tau:tval});  Up2_t =Up2_e.eval({tau:tval})
    B0_t =B0_e.eval({tau:tval});  B1_t =B1_e.eval({tau:tval}); B2_t=B2_e.eval({tau:tval})
    Up_t=Poly(sp.expand(Up_t.as_expr() if hasattr(Up_t,'as_expr') else Up_t),eta,domain=sp.QQ)
    Up2_t=Poly(sp.expand(Up2_t.as_expr() if hasattr(Up2_t,'as_expr') else Up2_t),eta,domain=sp.QQ)
    B0_t=Poly(sp.expand(B0_t.as_expr() if hasattr(B0_t,'as_expr') else B0_t),eta,domain=sp.QQ)
    B1_t=Poly(sp.expand(B1_t.as_expr() if hasattr(B1_t,'as_expr') else B1_t),eta,domain=sp.QQ)
    B2_t=Poly(sp.expand(B2_t.as_expr() if hasattr(B2_t,'as_expr') else B2_t),eta,domain=sp.QQ)
    # eta-orders
    print("  ord_eta: Up=%d Up2=%d  B0=%d B1=%d B2=%d"%(
        v_order(Up_t,eta),v_order(Up2_t,eta),v_order(B0_t,eta),v_order(B1_t,eta),v_order(B2_t,eta)),flush=True)
    # leading eta coeffs (constants)
    def lead(P):
        k=v_order(P,eta); return k, P.eval({eta:0}) if False else sp.LC(P,eta)
    for nm,Pp in [("Up",Up_t),("Up2",Up2_t),("B0",B0_t),("B1",B1_t),("B2",B2_t)]:
        k=v_order(Pp,eta)
        # coeff of eta^k
        ck=Poly(sp.expand(Pp.as_expr()),eta,domain=sp.QQ).nth(k) if Pp!=0 else 0
        print("    %s: eta-order %d, leading coeff eta^%d = %s"%(nm,k,k,ck),flush=True)
    # Transverse (U/v=0): Up + v*Up2 + ... =0.  If Up ~ a eta^r, Up2(0)=c2:
    #   a eta^r + v*c2 ~ 0  => v ~ -(a/c2) eta^r   (if c2!=0)
    rA=v_order(Up_t,eta); aA=Poly(sp.expand(Up_t.as_expr()),eta,domain=sp.QQ).nth(rA)
    c2 =Up2_t.eval({eta:0}) if Up2_t!=0 else 0
    print("  Transverse: Up~%s*eta^%d ; Up2(0)=%s"%(aA,rA,c2),flush=True)
    if c2!=0:
        print("    => v ~ -(Up_lead/Up2(0)) eta^%d = %s * eta^%d"%(rA, -sp.Rational(aA,c2) if c2!=0 else '?', rA),flush=True)
    # Longitudinal B=0: B0 + v*B1 + ... =0. B0~b0 eta^s0, B1~b1 eta^s1.
    s0=v_order(B0_t,eta); b0=Poly(sp.expand(B0_t.as_expr()),eta,domain=sp.QQ).nth(s0)
    s1=v_order(B1_t,eta); b1=Poly(sp.expand(B1_t.as_expr()),eta,domain=sp.QQ).nth(s1) if B1_t!=0 else 0
    print("  Longitudinal: B0~%s*eta^%d ; B1~%s*eta^%d"%(b0,s0,b1,s1),flush=True)
    if s1 is not None and b1!=0:
        # v ~ -(b0/b1) eta^(s0-s1)
        print("    => v ~ -(B0_lead/B1_lead) eta^%d = %s * eta^%d"%(s0-s1,-sp.Rational(b0,b1),s0-s1),flush=True)
print("DONE",flush=True)
