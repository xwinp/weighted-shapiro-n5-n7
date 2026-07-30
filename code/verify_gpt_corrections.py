#!/usr/bin/env python3
"""Independently verify GPT's corrections to the n=7 classification."""
import sympy as sp
import mpmath as mp
import numpy as np
from scipy.optimize import root
mp.mp.dps = 50

p, t = sp.symbols('p t', positive=True)
q = 1 - p
R = q**3 - p**3*t**5 - p**2*q*t**8
Bsimp = 5*p**2*t + 2*p*q*t**4 - 2*q**2 - 7*q*p**2

# ---------- (A) Resultant object check ----------
print("=== (A) resultant object check ===")
# B_simplified (deg 4 in t): Res_t(R, Bsimp)
RB = sp.factor(sp.resultant(R, Bsimp, t))
print("Res_t(R, B_simplified) =", RB)
# N = num(P_curve - 7) before R-simplification. Build P_curve from reduction.
d_e = t**2
c_e = q*(q - p*t**4)/(p**2*t**2)
b_e = q/(p*t) - q**2*(q - p*t**4)/(p**3*t**2)
# P on S2 with a=1,e=t
P = 1/(p*b_e+q*c_e) + b_e/(p*c_e+q*d_e) + c_e/(p*d_e) + d_e/(q*t) + t/q
N = sp.together(sp.simplify(P - 7)).as_numer_denom()[0]
N = sp.expand(N)
RN = sp.factor(sp.resultant(R, N, t))
print("Res_t(R, N=num(P_curve-7)) =", RN)
# also verify the identity P-7 = B/(q p^2) + R*(...): on R=0, P-7 = B/(q p^2)?
# check numerically at a7
a7=mp.mpf('0.21427352090984096558774231909001838797745454906427')
tv=mp.findroot(lambda tt:(1-a7)**3 - a7**3*tt**5 - a7**2*(1-a7)*tt**8, mp.mpf('1.36'))
def fv(e,pp,tt): return mp.mpf(sp.N(e.subs([(p,sp.Float(pp,45)),(t,sp.Float(tt,45))]),45))
print("at a7: P-7 =", mp.nstr(fv(P,a7,tv)-7,4), " B/(q p^2) =", mp.nstr(fv(Bsimp,a7,tv)/( (1-a7)*a7**2 ),4), " R =", mp.nstr(fv(R,a7,tv),4))

# ---------- (B) S3 = {0,3} closed form ----------
print("\n=== (B) S3={0,3} closed form (GPT: r=(p/q)^{1/5}, b=r^2,c=r^-1,e=r^3,d=r-r^8) ===")
def P_S3(pp):
    pv=mp.mpf(pp); qv=1-pv; r=(pv/qv)**mp.mpf('0.2')
    b=r**2; c=r**(-1); e=r**3; d=r-r**8
    P = 1/(pv*b) + b/(qv*c) + c/(pv*d+qv*e) + d/(pv*e) + e/qv
    return P, r, d
print("  p    r        d>0?    P_S3      P-7")
for pp in [0.1,0.214,0.25,0.329,0.4,0.49]:
    Pv,r,d = P_S3(pp)
    print(f"  {pp:.3f}  {float(r):.4f}  {d>0}  {float(Pv):.6f}  {float(Pv-7):+.6f}")
# verify it's a stationary point: grad P on S3 support {1,2,4,5,6}, vars b,c,d,e (x1=1)
bs,cs,ds,es = sp.symbols('b c d e', positive=True)
x={0:0,1:sp.Integer(1),2:bs,3:0,4:cs,5:ds,6:es}
PS3=sum(x[i]/(p*x[(i+1)%7]+q*x[(i+2)%7]) for i in range(7))
grads=[sp.simplify(sp.diff(PS3,v)) for v in [bs,cs,ds,es]]
# substitute GPT closed form symbolically with r
rs=sp.symbols('r',positive=True)
subs_sym=[(bs,rs**2),(cs,rs**(-1)),(ds,rs-rs**8),(es,rs**3),(p, rs**5/(1+rs**5)),(q,1/(1+rs**5))]
# p = r^5 q, q=1-p -> p = r^5/(1+r^5)
for i,g in enumerate(grads):
    gsub=sp.simplify(g.subs(subs_sym))
    print(f"  grad[{i}] after subs (should be 0):", sp.simplify(gsub))

# ---------- (C) S1 = {0} saddle at p=1/4 ----------
print("\n=== (C) S1={0} stationary point at p=1/4 (GPT: P=6.98998<7, saddle) ===")
def P_full(x,p):
    q=1-p; n=7; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if abs(den)<1e-15: return 1e6
        s+=x[i]/den
    return s
pp=0.25
# support {1,2,3,4,5,6}, x0=0, x1=1, vars b,c,d,e,f
from numpy import array
def grad_S1(v, p):
    # v=[b,c,d,e,f]; x=[0,1,b,c,d,e,f]
    h=1e-7; g=[]
    x0=array([0,1,*v],dtype=float)
    f0=P_full(x0,p)
    for i in range(5):
        xp=x0.copy(); xp[2+i]+=h; g.append((P_full(xp,p)-f0)/h)
    return g
# use GPT's init
init=[0.2684881167890583140,0.6791742990557855304,1.5461708324775024161,0.0656843931252869930,1.3009478193484040029]
r=root(lambda v: grad_S1(v,pp), init, method='hybr', options={'xtol':1e-13})
xsol=[0,1,*r.x]; Pval_s1=P_full(array(xsol),pp)
print(f"  S1 stationary at p=1/4: P={Pval_s1:.10f}  (<7: {Pval_s1<7})")
print(f"  b,c,d,e,f = {r.x}")
# Hessian (5x5) of P on S1 face (fix x1=1)
h=1e-5
H=np.zeros((5,5)); x0=array([0,1,*r.x],dtype=float)
for i in range(5):
    for j in range(5):
        xp=x0.copy(); xm=x0.copy(); xpp=x0.copy(); xmm=x0.copy()
        xp[2+i]+=h; xp[2+j]+=h; xm[2+i]-=h; xm[2+j]-=h
        xpp[2+i]+=h; xpp[2+j]-=h; xmm[2+i]-=h; xmm[2+j]+=h
        H[i,j]=(P_full(xp,pp)-P_full(xm,pp)-P_full(xpp,pp)+P_full(xmm,pp))/(4*h*h)
eig=np.linalg.eigvalsh(H)
print(f"  Hessian eigenvalues: {np.round(eig,4)}  (any negative -> saddle: {(eig<0).any()})")

# ---------- (D) S0 non-uniform stationary at p=1/4 ----------
print("\n=== (D) S0 non-uniform stationary at p=1/4 (GPT: P=7.0514, saddle) ===")
# full support, x0=1, vars b,c,d,e,f,g
init0=[0.2899598915706492,0.7447488449365604,0.8875842694414185,0.3764381331194338,1.1665179227192974,0.1591166088535238]
def grad_S0(v,p):
    h=1e-7; g=[]; x0=array([1,*v],dtype=float); f0=P_full(x0,p)
    for i in range(6):
        xp=x0.copy(); xp[1+i]+=h; g.append((P_full(xp,p)-f0)/h)
    return g
r0=root(lambda v: grad_S0(v,pp), init0, method='hybr', options={'xtol':1e-13})
xsol0=[1,*r0.x]; Pv0=P_full(array(xsol0),pp)
print(f"  S0 non-uniform stationary at p=1/4: P={Pv0:.10f}  (uniform P=7)")
h=1e-5
H0=np.zeros((6,6)); x0=array([1,*r0.x],dtype=float)
for i in range(6):
    for j in range(6):
        xp=x0.copy(); xm=x0.copy(); xpp=x0.copy(); xmm=x0.copy()
        xp[1+i]+=h; xp[1+j]+=h; xm[1+i]-=h; xm[1+j]-=h
        xpp[1+i]+=h; xpp[1+j]-=h; xmm[1+i]-=h; xmm[1+j]+=h
        H0[i,j]=(P_full(xp,pp)-P_full(xm,pp)-P_full(xpp,pp)+P_full(xmm,pp))/(4*h*h)
eig0=np.linalg.eigvalsh(H0)
print(f"  Hessian eigenvalues: {np.round(eig0,4)}  (any negative -> saddle: {(eig0<0).any()})")
