#!/usr/bin/env python3
"""
n=7 orbit O3 (2-zero distance 3, zeros {0,4}, support {1,2,3,5,6}) KKT study.
x=(0,a,b,c,0,d,e), a=1, e=t:
  P = 1/(pb+qc) + b/(pc) + c/(qd) + d/(pt) + t/q
KKT (Euler lambda=0): dP/db=dP/dc=dP/dd=0. Solve numerically across p,
confirm interior stationary (grad~0), inactive KKT for zeros {0,4} >=0,
and P>7 everywhere (no failure band). Compare to O2.
"""
import numpy as np
from scipy.optimize import minimize, root
import mpmath as mp

def Pvec(x, p):
    q=1-p; n=7; s=0.0
    for i in range(n):
        den=p*x[(i+1)%n]+q*x[(i+2)%n]
        if den<=1e-15: return 1e6
        s+=x[i]/den
    return s

# zeros {0,4}; free coords 1,2,3,5,6 = (a,b,c,d,e)
def solve_kkt(p, init=None):
    # vars y=[b,c,d,t] with a=1,e=t;  KKT dP/d{b,c,d}=0 (3 eqs) + P-coupling? use 3 eqs
    q=1-p
    def Pbcd(b,c,d,t):
        return 1/(p*b+q*c) + b/(p*c) + c/(q*d) + d/(p*t) + t/q
    def grad(b,c,d,t):
        # numerical gradient wrt b,c,d
        h=1e-7
        f0=Pbcd(b,c,d,t)
        gb=(Pbcd(b+h,c,d,t)-f0)/h
        gc=(Pbcd(b,c+h,d,t)-f0)/h
        gd=(Pbcd(b,c,d+h,t)-f0)/h
        return [gb,gc,gd]
    # also need dP/dt? t=e is a free var too (e not normalized away). Actually a=1 fixed, e=t free -> dP/dt=0 too (4th eq). But Euler makes one dep. Use 4 eqs, 4 unknowns b,c,d,t.
    def grad4(b,c,d,t):
        h=1e-7; f0=Pbcd(b,c,d,t)
        return [(Pbcd(b+h,c,d,t)-f0)/h,(Pbcd(b,c+h,d,t)-f0)/h,(Pbcd(b,c,d+h,t)-f0)/h,(Pbcd(b,c,d,t+h)-f0)/h]
    if init is None:
        init=[1.0,0.4,1.8,1.3]
    try:
        r=root(lambda v: grad4(*v), init, method='hybr', options={'xtol':1e-12})
        if r.success and all(abs(g)<1e-6 for g in r.fun):
            return r.x, Pbcd(*r.x)
    except Exception:
        pass
    return None, None

ps=np.linspace(0.05,0.49,45)
print("p      b      c      d      t      P       P-7")
prev=None
best=(9,None)
for pp in ps:
    v,Pv=solve_kkt(pp, init=prev)
    if v is None:
        # try fresh
        v,Pv=solve_kkt(pp)
    if v is not None:
        prev=v
        b,c,d,t=v
        print(f"{pp:.3f}  {b:.4f} {c:.4f} {d:.4f} {t:.4f}  {Pv:.5f}  {Pv-7:+.5f}")
        if Pv<best[0]: best=(Pv,pp)
    else:
        print(f"{pp:.3f}  no convergence")
print(f"\nO3 min P = {best[0]:.6f} at p={best[1]}")

# inactive KKT for O3: zeros at {0,4}. Build full P with x0,x4 symbolic, compute dP/dx0, dP/dx4 at 0.
import sympy as sp
ps_,t_ = sp.symbols('p t', positive=True)
q_ = 1-ps_
b_,c_,d_ = sp.symbols('b c d', positive=True)
x0s,x4s = sp.symbols('x0 x4', nonnegative=True)
xx={0:x0s,1:sp.Integer(1),2:b_,3:c_,4:x4s,5:d_,6:t_}
def xi(i): return xx[i]
Pxx=sum(xi(i)/(ps_*xi((i+1)%7)+q_*xi((i+2)%7)) for i in range(7))
# substitute a KKT solution at p=0.25 (near min)
pp=0.25
v,Pv=solve_kkt(pp)
b_v,c_v,d_v,t_v=v
subs=[(b_,sp.Float(b_v,20)),(c_,sp.Float(c_v,20)),(d_,sp.Float(d_v,20)),(t_,sp.Float(t_v,20)),(ps_,sp.Float(pp,20))]
D0=sp.diff(Pxx,x0s).subs([(x0s,0),(x4s,0)])
D4=sp.diff(Pxx,x4s).subs([(x0s,0),(x4s,0)])
print(f"\nO3 inactive KKT at p={pp}: dP/dx0={float(sp.N(D0.subs(subs))):.4f}  dP/dx4={float(sp.N(D4.subs(subs))):.4f}  (need >=0 for zeros to stay)")
