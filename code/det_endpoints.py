import sympy as sp
import numpy as np
b,c,d,e,f,p=sp.symbols('b c d e f p',positive=True)
q=1-p
P = 1/(p*b+q*c) + b/(p*c+q*d) + c/(p*d+q*e) + d/(p*e+q*f) + e/(p*f) + f/q
vars5=[b,c,d,e,f]
grads=[sp.diff(P,v) for v in vars5]
Hsym=sp.Matrix(5,5,lambda i,j: sp.diff(P,vars5[i],vars5[j]))
a7,b7=0.214273520909841,0.328627677916592
def hess_num(sol,pp):
    def Pn(vv):
        qv=1-pp;x=[0,1,vv[0],vv[1],vv[2],vv[3],vv[4]];s=0.0
        for i in range(7):
            den=pp*x[(i+1)%7]+qv*x[(i+2)%7];s+=x[i]/den
        return s
    H=np.zeros((5,5));h=1e-4;v=np.array([float(s) for s in sol])
    for i in range(5):
        ei=np.zeros(5);ei[i]=1
        H[i,i]=(Pn(v+h*ei)-2*Pn(v)+Pn(v-h*ei))/h**2
        for j in range(i+1,5):
            ej=np.zeros(5);ej[j]=1
            H[i,j]=H[j,i]=(Pn(v+h*ei+h*ej)-Pn(v+h*ei-h*ej)-Pn(v-h*ei+h*ej)+Pn(v-h*ei-h*ej))/(4*h*h)
    return H
# scan very close to endpoints, continuation from inside
prev=[0.27,0.68,1.5,0.07,1.3]
# go toward a7 from 0.25
ps_a=list(np.linspace(0.25,a7+0.001,10))
ps_b=list(np.linspace(0.30,b7-0.001,10))
print("toward a7:")
for pp in ps_a:
    gsub=[g.subs(p,sp.Float(pp)) for g in grads]
    try: sol=sp.nsolve(gsub,vars5,prev,prec=35,tol=1e-30,maxsteps=200)
    except: continue
    prev=[float(s) for s in sol]
    H=hess_num(sol,pp);dv=float(np.linalg.det(H))
    print(f"  p={pp:.5f} det={dv:+.4f} min(x)={min(float(s) for s in sol):.4f}")
print("toward b7:")
prev=[0.30,0.69,1.3,0.10,1.21]
for pp in ps_b:
    gsub=[g.subs(p,sp.Float(pp)) for g in grads]
    try: sol=sp.nsolve(gsub,vars5,prev,prec=35,tol=1e-30,maxsteps=200)
    except: continue
    prev=[float(s) for s in sol]
    H=hess_num(sol,pp);dv=float(np.linalg.det(H))
    print(f"  p={pp:.5f} det={dv:+.4f} min(x)={min(float(s) for s in sol):.4f}")
