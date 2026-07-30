#!/usr/bin/env python3
"""Push algebraic identification: higher dps, larger coeffs, wider degree range."""
import mpmath as mp
mp.mp.dps = 200

# re-solve at dps 200 for max reliability
import sympy as sp
a,b,c,d,e,lam,p = sp.symbols('a b c d e lam p', positive=True)
q=1-p
P=(a/(p*b+q*c)+b/(p*c+q*d)+c/(p*d)+d/(q*e)+e/(q*a))
gs=[sp.diff(P,v) for v in (a,b,c,d,e)]
eqs=[g-lam for g in gs]+[a+b+c+d+e-1, P-7]
vs=[a,b,c,d,e,lam,p]
def Pe(av,bv,cv,dv,ev,pv):
    qv=1-pv; return av/(pv*bv+qv*cv)+bv/(pv*cv+qv*dv)+cv/(pv*dv)+dv/(qv*ev)+ev/(qv*av)

roots={}
for name,ini in [("a_7",[0.17,0.23,0.06,0.31,0.23,6.0,0.22]),("b_7",[0.20,0.20,0.08,0.28,0.24,6.0,0.33])]:
    sol=sp.nsolve(eqs,vs,ini,prec=210,tol=mp.mpf('1e-200'),maxsteps=400)
    pv=mp.mpf(sp.N(sol[6],200))
    roots[name]=pv
    print(f"{name} = {mp.nstr(pv,80)}")

a7,b7=roots["a_7"],roots["b_7"]
print(f"\nsum  = {mp.nstr(a7+b7,60)}")
print(f"prod = {mp.nstr(a7*b7,60)}")

print("\n=== PSLQ wide sweep (dps=200, maxcoeff=1e8) ===")
for name,t in [("a_7",a7),("b_7",b7)]:
    print(f"\n{name}:")
    hit=None
    for deg in [6,8,10,12,14,16,18,20,24]:
        vec=[t**i for i in range(deg+1)]
        rel=mp.pslq(vec,maxcoeff=10**8,maxsteps=20000)
        if rel:
            resid=sum(rel[i]*t**i for i in range(deg+1))
            other=b7 if name=="a_7" else a7
            resid2=sum(rel[i]*other**i for i in range(deg+1))
            ok = abs(resid)<mp.mpf('1e-40') and abs(resid2)<mp.mpf('1e-20')
            print(f"  deg {deg}: resid={mp.nstr(resid,4)} other_resid={mp.nstr(resid2,4)} {'*** COMMON' if ok else ''} coeffs={rel[:6]}...")
            if ok and hit is None: hit=(deg,rel)
    if hit:
        deg,rel=hit
        print(f"  >>> minimal poly candidate deg {deg}: {rel}")
