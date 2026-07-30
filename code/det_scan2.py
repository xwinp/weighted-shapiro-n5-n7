import sympy as sp
import numpy as np
b,c,d,e,f,p=sp.symbols('b c d e f p',positive=True)
q=1-p
P = 1/(p*b+q*c) + b/(p*c+q*d) + c/(p*d+q*e) + d/(p*e+q*f) + e/(p*f) + f/q
vars5=[b,c,d,e,f]
grads=[sp.diff(P,v) for v in vars5]
a7,b7=0.214273520909841,0.328627677916592
def hess_num(sol,pp):
    def Pn(vv):
        qv=1-pp;x=[0,1,vv[0],vv[1],vv[2],vv[3],vv[4]]
        s=0.0
        for i in range(7):
            den=pp*x[(i+1)%7]+qv*x[(i+2)%7];s+=x[i]/den
        return s
    H=np.zeros((5,5));h=1e-4;v=np.array([float(s) for s in sol])
    for i in range(5):
        ei=np.zeros(5);ei[i]=1
        H[i,i]=(Pn(v+h*ei)-2*Pn(v)+Pn(v-h*ei))/h**2
        for j in range(i+1,5):
            ej=np.zeros(5);ej[j]=1
            Hij=(Pn(v+h*ei+h*ej)-Pn(v+h*ei-h*ej)-Pn(v-h*ei+h*ej)+Pn(v-h*ei-h*ej))/(4*h*h)
            H[i,j]=H[j,i]=Hij
    return H
pp_list=np.linspace(a7+0.003,b7-0.003,16)
prev=[0.27,0.68,1.5,0.07,1.3]
dets=[]
print(f"{'p':>7} {'det(H)':>12} {'P_S1':>10} {'Morse':>6}")
for pp in pp_list:
    gsub=[g.subs(p,sp.Float(pp)) for g in grads]
    try: sol=sp.nsolve(gsub,vars5,prev,prec=35,tol=1e-30,maxsteps=150)
    except: 
        try: sol=sp.nsolve(gsub,vars5,[0.27,0.68,1.5,0.07,1.3],prec=35,tol=1e-30,maxsteps=150)
        except Exception as ex: print(f"{pp:.4f} fail {ex}");continue
    prev=[float(s) for s in sol]
    H=hess_num(sol,pp);ev=np.linalg.eigvalsh(H);dv=float(np.linalg.det(H))
    Pv=float(P.subs(p,sp.Float(pp)).subs({b:sol[0],c:sol[1],d:sol[2],e:sol[3],f:sol[4]}))
    dets.append(dv)
    print(f"{pp:7.4f} {dv:12.4f} {Pv:10.5f} {int(np.sum(ev<0)):6d}")
print(f"\ndet range: [{min(dets):.4f}, {max(dets):.4f}]  -> {'NEG THROUGHOUT (cert OK)' if max(dets)<0 else 'NOT all neg'}")
