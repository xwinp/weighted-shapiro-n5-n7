import sympy as sp
import numpy as np
b,c,d,e,f,p=sp.symbols('b c d e f p',positive=True)
q=1-p
P = 1/(p*b+q*c) + b/(p*c+q*d) + c/(p*d+q*e) + d/(p*e+q*f) + e/(p*f) + f/q
vars5=[b,c,d,e,f]
grads=[sp.diff(P,v) for v in vars5]
a7,b7=0.214273520909841,0.328627677916592
# scan p, track S1 stationary P-value. Use continuation: feed prev sol as init.
pp_list=np.linspace(0.40,0.15,26)
init=[0.34,0.70,1.2,0.12,1.17]
print(f"{'p':>7} {'P_S1':>10} {'P-7':>10}  region")
prev=None
for pp in pp_list:
    gsub=[g.subs(p,sp.Float(pp)) for g in grads]
    use_init = prev if prev is not None else init
    try:
        sol=sp.nsolve(gsub,vars5,use_init,prec=40,tol=1e-35,maxsteps=200)
    except Exception:
        # retry with base init
        try:
            sol=sp.nsolve(gsub,vars5,init,prec=40,tol=1e-35,maxsteps=200)
        except Exception as ex:
            print(f"{pp:7.4f}   NSOLVE FAIL"); continue
    subs={b:sol[0],c:sol[1],d:sol[2],e:sol[3],f:sol[4]}
    Pv=float(P.subs(p,sp.Float(pp)).evalf(subs=subs))
    pos = all(float(s)>0 for s in sol)
    prev=[float(s) for s in sol]
    region = 'BAND' if a7<pp<b7 else ('hold' if (pp<=a7 or pp>=b7) else '?')
    print(f"{pp:7.4f} {Pv:10.5f} {Pv-7:+10.5f}  {region}  pos={pos}")
