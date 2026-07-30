import sympy as sp
import numpy as np
b,c,d,e,f,p=sp.symbols('b c d e f p',positive=True)
q=1-p
P = 1/(p*b+q*c) + b/(p*c+q*d) + c/(p*d+q*e) + d/(p*e+q*f) + e/(p*f) + f/q
vars5=[b,c,d,e,f]
grads=[sp.diff(P,v) for v in vars5]
Hsym=sp.Matrix(5,5,lambda i,j: sp.diff(P,vars5[i],vars5[j]))
for pp in [0.24,0.27,0.30,0.32]:
    gsub=[g.subs(p,sp.Rational(pp).limit_denom(1000) if False else sp.Float(pp)) for g in grads]
    Hsub=Hsym.subs(p,sp.Float(pp))
    init=[0.27,0.68,1.5,0.07,1.3]
    try:
        sol=sp.nsolve(gsub,vars5,init,prec=50,tol=1e-45,maxsteps=300)
    except Exception as ex:
        print(f"p={pp} nsolve fail: {ex}");continue
    subs={b:sol[0],c:sol[1],d:sol[2],e:sol[3],f:sol[4]}
    gv=[float(g.evalf(subs=subs)) for g in gsub]
    Hn=np.array(Hsub.evalf(subs=subs),dtype=float)
    ev=np.linalg.eigvalsh(Hn);det=float(np.linalg.det(Hn));Pv=float(P.subs(p,sp.Float(pp)).evalf(subs=subs))
    print(f"p={pp}: sol={[round(float(s),6) for s in sol]} grad_max={max(abs(x) for x in gv):.1e} P={Pv:.6f} det={det:+.4f} Morse={int(np.sum(ev<0))} eigs={np.round(ev,4)}")
