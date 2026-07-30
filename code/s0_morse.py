import sympy as sp
import numpy as np
# S0 full interior, x0=1 anchor, vars x1..x6
xs=sp.symbols('x1:7',positive=True)  # x1..x6
p=sp.symbols('p',positive=True);q=1-p
xfull=[1]+list(xs)  # x0=1, x1..x6
n=7
P=sum(xfull[i]/(p*xfull[(i+1)%n]+q*xfull[(i+2)%n]) for i in range(n))
P=sp.together(P)
grads=[sp.diff(P,v) for v in xs]
Hsym=sp.Matrix(6,6,lambda i,j: sp.diff(P,xs[i],xs[j]))
# tangent basis for reduced (quotient radial in x0..x6): use sum-zero on x1..x6 with x0 fixed
# reduced Hessian = restrict to {sum x1..x6 = const} transverse; use B: 6x5, cols e_j - e_5
def tangent6():
    B=np.zeros((6,5))
    for j in range(5): B[j,j]=1;B[5,j]=-1
    return B
B=tangent6()
for pp in [0.25,0.27,0.30,0.32]:
    gsub=[g.subs(p,sp.Float(pp)) for g in grads]
    Hsub=Hsym.subs(p,sp.Float(pp))
    # init: GPT p=1/4 nonunif [x0=1,0.28996,0.74475,0.88758,0.37644,1.16652,0.15912]
    init=[0.28996,0.74475,0.88758,0.37644,1.16652,0.15912]
    try:
        sol=sp.nsolve(gsub,list(xs),init,prec=40,tol=1e-35,maxsteps=300)
    except Exception as ex:
        print(f"p={pp} nsolve fail: {ex}");continue
    subs={xs[i]:sol[i] for i in range(6)}
    gv=[float(g.evalf(subs=subs)) for g in gsub]
    Hn=np.array(Hsub.evalf(subs=subs),dtype=float)
    Hred=B.T@Hn@B
    ev=np.linalg.eigvalsh(Hred)
    Pv=float(P.subs(p,sp.Float(pp)).evalf(subs=subs))
    print(f"p={pp}: grad_max={max(abs(x) for x in gv):.1e} P={Pv:.6f} Morse={int(np.sum(ev<0))} eigs={np.round(ev,4)}")
    print(f"   sol={[round(float(s),5) for s in sol]}")
