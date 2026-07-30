#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rigorous rational-interval certificates for the n=7 S2 face (GPT review items
4.2 #2, #3): (a) three-sign samples of P_curve-7 at rational p=1/5,1/4,1/3
  (in the three regions p<a7, a7<p<b7, p>b7, using a7 in (1/5,1/4), b7 in (1/4,1/3));
(b) inactive derivatives D0,D5 > 0 at p=1/4.

For each rational p: isolate the unique positive root t of R(p,t)=q^3-p^3 t^5-p^2 q t^8
(q=1-p) in a rational interval [tlo,thi] (R(tlo)>0>R(thi), R monotone decreasing),
then evaluate B=5 p^2 t+2 p q t^4-2 q^2-7 q p^2  (sign of P_curve-7 = sign of B,
since P_curve-7 = B/(q p^2) + R*(...)) and D0,D5 via rational interval arithmetic.
"""
import sympy as sp
t,p=sp.symbols('t p'); q=1-p
# KKT recovery on S2: x=(0,1,b,c,d,0,e=t), d=t^2
d=t**2
c=q*(q-p*t**4)/(p**2*t**2)
b=q/(p*t)-q**2*(q-p*t**4)/(p**3*t**2)
# P on S2 (x1=1): terms i=0..6
x0,x1,x2,x3,x4,x5,x6=0,1,b,c,d,0,t
xs=[x0,x1,x2,x3,x4,x5,x6]
n=7
P=sum(xs[i]/(p*xs[(i+1)%n]+q*xs[(i+2)%n]) for i in range(n))
P=sp.together(P)
# B = num of (P-7) after substituting R=0; paper says P_curve-7 = B/(q p^2) (mod R)
B=5*p**2*t+2*p*q*t**4-2*q**2-7*q*p**2
# inactive derivatives D0 = dP/dx0, D5 = dP/dx5  at the S2 point (x0=x5=0)
# compute symbolically as functions of (p,t) with b,c,d above
Px=[sp.diff(P,sp.Symbol('xx%d'%i)) for i in range(0)] # placeholder; recompute directly
# Recompute P with symbolic x to get derivatives
xb=sp.symbols('x0:7')
Ps=sum(xb[i]/(p*xb[(i+1)%7]+q*xb[(i+2)%7]) for i in range(7))
dP=[sp.diff(Ps,xb[i]) for i in range(7)]
subs={xb[0]:0,xb[1]:1,xb[2]:b,xb[3]:c,xb[4]:d,xb[5]:0,xb[6]:t}
D0=sp.together(dP[0].subs(subs))
D5=sp.together(dP[5].subs(subs))

def isol_t(pval, denom=10_000_000, half=200):
    Rv=q**3-p**3*t**5-p**2*q*t**8
    Rv=Rv.subs(p,pval)
    rts=sp.real_roots(sp.Poly(Rv,t))
    pos=[r for r in rts if r>0]
    if not pos: return None
    tv=float(pos[0])
    lo=sp.Rational(int((tv)*denom)-half, denom) if tv*denom>half else sp.Rational(1,denom)
    hi=sp.Rational(int((tv)*denom)+half, denom)
    while Rv.subs(t,lo)<=0: lo-=sp.Rational(1,denom)
    while Rv.subs(t,hi)>=0: hi+=sp.Rational(1,denom)
    return lo,hi,tv

def iv_eval(expr,pval,lo,hi):
    """Rational interval eval: p=pval, t in [lo,hi]. Bound num/den (polys in t)
    monomial-wise: k=0 term is the constant ck itself; k>0 term ck*[lo^k,hi^k]
    (reversed if ck<0).  Assumes den>0 on [lo,hi]."""
    e=sp.together(expr).subs(p,pval)
    num,den=sp.fraction(e)
    num=sp.Poly(num,t); denv=sp.Poly(den,t)
    def bound_poly(Po):
        if Po.degree()<0:
            v=Po.as_expr(); return v,v
        terms=[]
        for k in range(Po.degree()+1):
            ck=Po.nth(k)
            if k==0:
                terms.append((ck,ck))
            else:
                lo_k=lo**k; hi_k=hi**k
                if ck>=0: terms.append((ck*lo_k,ck*hi_k))
                else: terms.append((ck*hi_k,ck*lo_k))
        return sum(a for a,_ in terms), sum(b for _,b in terms)
    nlo,nhi=bound_poly(num); dlo,dhi=bound_poly(denv)
    if dlo<=0 or dhi<=0: return None  # denominator not sign-definite positive
    return nlo/dhi, nhi/dlo   # den>0: [num_lo/den_hi, num_hi/den_lo]

for pval,region in [(sp.Rational(1,5),"p=1/5<a7 (expect P>7, B>0)"),
                    (sp.Rational(1,4),"p=1/4 in (a7,b7) (expect P<7, B<0)"),
                    (sp.Rational(1,3),"p=1/3>b7 (expect P>7, B>0)")]:
    res=isol_t(pval)
    if not res: print(region,"no root"); continue
    lo,hi,tv=res
    Biv=iv_eval(B,pval,lo,hi)
    sgn = "+" if (Biv and Biv[0]>0) else ("-" if (Biv and Biv[1]<0) else "?")
    print("%-42s t~%.5f in[%.5f,%.5f]  B in[%s,%s] -> sign %s"%(region,tv,float(lo),float(hi),
          sp.nsimplify(Biv[0]),sp.nsimplify(Biv[1]),sgn) if Biv else region+" eval fail")

# D0,D5 at p=1/4
pval=sp.Rational(1,4)
lo,hi,tv=isol_t(pval)
for name,Dexpr in [("D0",D0),("D5",D5)]:
    iv=iv_eval(Dexpr,pval,lo,hi)
    if iv:
        print("p=1/4 %s in [%s, %s]  >0? %s  (float [%.5f,%.5f])"%(name,sp.nsimplify(iv[0]),sp.nsimplify(iv[1]),iv[0]>0,float(iv[0]),float(iv[1])))
    else:
        print(name,"eval fail")
print("DONE")
