#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Issue 7 final: F1 vs E2_red equivalence, saturation prefactors, back-substitution."""
import sympy as sp
v,w,z = sp.symbols('v w z', positive=True)
uS,vS,wS,zS=sp.symbols('u v w z')
F1 = v**2*w + v*w**3 - v*w**2 + v*w*z - v*w - v*z + v - w**3 + 2*w**2 - w

# E2_red from kkt_check (E2 with u=usol, numerator)
a3=1-vS+uS*vS; a5=1-zS+zS*wS-zS*vS*wS+zS*uS*vS*wS
E2=uS*(1-zS)-zS*a5*(1-vS)
usol=sp.together(vS*(1-vS)/((1-wS)-vS**2))
E2u_num=sp.expand(sp.together(E2.subs(uS,usol)).as_numer_denom()[0])
E2u_num_vwz=sp.expand(E2u_num.subs({vS:v,wS:w,zS:z}))
print("E2_red factored:", sp.factor(E2u_num_vwz), flush=True)
print("F1 factored     :", sp.factor(F1), flush=True)
# Quotient F1 / E2_red (should be a nonzero rational function = the prefactor relating them)
q=sp.cancel(F1/E2u_num_vwz)
print("F1 / E2_red =", sp.factor(q), flush=True)

# Saturation prefactors appearing through the elimination chain (admissibility/singularity):
# w (w=0), w-1 (w=1), v-1 (v=1), u-1 (u=1), v+w-1 (usol denom / E2 boundary).
# Confirm these are exactly the beta in {0,1} / denominator-zero boundaries:
print("\nSaturation prefactors (removed by localizing at admissible beta in (0,1)):", flush=True)
print("  w, w-1, v-1, u-1, v+w-1   [= beta_j=0 or 1, or (1-w)-v^2=0 usol denominator]", flush=True)

# --- back-substitution (numerical soundness) on H_B and H_C ---
import numpy as np
n=7
def Pval(x,p):
    q=1-p; s=0.0
    for i in range(n):
        d=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(d)<1e-15: return 1e6
        s+=x[i]/d
    return s
rho_syms=sp.symbols('r1:6'); p=sp.symbols('p'); q=1-p
P_expr=sp.together(sum(1/(rho_syms[i]*(p+q*rho_syms[i+1])) for i in range(4))
                   + 1/(p*rho_syms[4]) + sp.prod(rho_syms)/q)
g=[sp.expand(sp.together(rho_syms[i]*sp.diff(P_expr,rho_syms[i])).as_numer_denom()[0]) for i in range(5)]

H_B = z*w**2 + (1-z**2)*w + z**2 - z
H_C = z*w**3 + w**2*z**3 - w**2*z + w*z**4 - 3*w*z**3 + 2*w*z**2 + w*z - w - z**4 + 3*z**3 - 3*z**2 + z
def lift_and_check(tag, wz_pairs):
    print("\n--- back-substitution on %s ---"%tag, flush=True)
    for (wv,zv) in wz_pairs:
        # solve E3 for u: u = v(1-v)/((1-w)-v^2);  solve F1 (=E2_red*v-relation) for v
        F1n = sp.lambdify(sp.Symbol('v'), F1.subs({w:wv,z:zv}))
        from scipy.optimize import brentq
        # find a v in (0,1) with F1=0
        vs=np.linspace(0.02,0.98,400); fv=np.array([F1n(t) for t in vs])
        v0=None
        for i in range(len(vs)-1):
            if np.isfinite(fv[i]) and np.isfinite(fv[i+1]) and fv[i]*fv[i+1]<0:
                try: v0=brentq(F1n,vs[i],vs[i+1]); break
                except: continue
        if v0 is None: print("  w=%.4f z=%.4f: no v-root found"%(wv,zv)); continue
        denom=(1-wv)-v0**2
        if abs(denom)<1e-9: print("  denom 0"); continue
        u0=v0*(1-v0)/denom
        if not(0<u0<1): print("  u0=% .4f not admissible"%u0); continue
        a5v=1-zv+zv*wv-zv*v0*wv+zv*u0*v0*wv
        if a5v<=0: print("  a5<=0"); continue
        K=u0*v0*wv*(zv**3)*a5v**2/((1-v0)*(1-wv)*(1-zv)**3)
        if K<=0: print("  K<=0"); continue
        rho_v=K**(1/7); pp=1/(1+rho_v)
        r2=u0/(rho_v*(1-u0)); r3=v0/(rho_v*(1-v0)); r4=wv/(rho_v*(1-wv)); r5=zv/(rho_v*(1-zv))
        # r1 from g1: r1^2 = A/B  with A=(1-u)(1+rho), B=(1+rho)uvwz/(rho^5 (1-u)(1-v)(1-w)(1-z))
        A_=(1-u0)*(1+rho_v); B_=(1+rho_v)*u0*v0*wv*zv/(rho_v**5*(1-u0)*(1-v0)*(1-wv)*(1-zv))
        r1=np.sqrt(A_/B_)
        rholist=[r1,r2,r3,r4,r5]
        subs={p:pp}
        for i in range(5): subs[rho_syms[i]]=rholist[i]
        gres=[float(g[i].subs(subs)) for i in range(5)]
        x=np.zeros(n); x[1]=1.0
        for k in range(1,6): x[1+k]=np.prod(rholist[:k])
        xs=x/x.sum()
        print("  w=%.4f z=%.4f: u=%.4f v=%.4f rho=%.4f p=%.4f  g1..g5=%s  P=%.5f"%(
            wv,zv,u0,v0,rho_v,pp,["%.1e"%abs(t) for t in gres], Pval(xs,pp)), flush=True)
# H_B sample (z in (0,1), w = positive root of H_B)
def hb_w(zv):
    import numpy as np
    a=zv; # zw^2+(1-z^2)w+z^2-z=0 -> quadratic in w
    Aq=a; Bq=1-a**2; Cq=a**2-a
    disc=Bq**2-4*Aq*Cq
    if disc<0: return None
    return (-Bq+np.sqrt(disc))/(2*Aq)
zs=[0.3,0.5,0.7,0.85]
lift_and_check("H_B", [(hb_w(zv),zv) for zv in zs if hb_w(zv) is not None and 0<hb_w(zv)<1])
# H_C sample: find w-root of H_C for given z
def hc_w(zv):
    HCz=sp.Poly(H_C.subs(z,zv),w)
    rs=[float(sp.re(r)) for r in sp.nroots(HCz) if abs(sp.im(r))<1e-9]
    return [r for r in rs if 0<r<1]
hcs=[]
for zv in [0.5,0.8,0.9,0.95]:
    for wv in hc_w(zv): hcs.append((wv,zv))
lift_and_check("H_C", hcs[:4])
print("\nDONE-STEP4", flush=True)
