import sympy as sp
import numpy as np
# I2: P_S0^nonunif - 7  (S0, x0=1, vars x1..x6)
xs=sp.symbols('x1:7',positive=True);p=sp.symbols('p',positive=True);q=1-p
xf=[1]+list(xs);n=7
P0=sum(xf[i]/(p*xf[(i+1)%n]+q*xf[(i+2)%n]) for i in range(n))
g0=[sp.diff(P0,v) for v in xs]
# I1: P_S1^stat - P_S2^curve. S1: x0=0,x1=1,vars b,c,d,e,f
b,c,d,e,f=sp.symbols('b c d e f',positive=True)
P1=1/(p*b+q*c)+b/(p*c+q*d)+c/(p*d+q*e)+d/(p*e+q*f)+e/(p*f)+f/q
g1=[sp.diff(P1,v) for v in [b,c,d,e,f]]
# S2 curve: x=(0,a,1,b,c,0,d) a=1... use prior R(p,t)=0 param. Easier: solve S2 KKT numerically.
# S2: x=(0,1,B,C,D,0,E) vars B,C,D,E (x0=0,x1=1,x5=0). 
B,C,D,E=sp.symbols('B C D E',positive=True)
# P_S2 = x0/(..)=0 + x1/(p B+q C)+ B/(p C+q D)+ C/(p D+q*0)+ D/(p*0+q E)+ x5/(..)=0 + E/(q*0+p*1)=E/p
# wait indices: x=(0,1,B,C,D,0,E). terms:
# i0:0/(p*1+q*B)=0 ; i1:1/(p*B+q*C); i2:B/(p*C+q*D); i3:C/(p*D+q*0)=C/(pD); i4:D/(p*0+q*E)=D/(qE);
# i5:0/(p*E+q*0)=0 ; i6:E/(p*0+q*1)=E/q
P2=1/(p*B+q*C)+B/(p*C+q*D)+C/(p*D)+D/(q*E)+E/q
g2=[sp.diff(P2,v) for v in [B,C,D,E]]
a7,b7=0.214273520909841,0.328627677916592
s0_prev=[0.28996,0.74475,0.88758,0.37644,1.16652,0.15912]
s1_prev=[0.27,0.68,1.5,0.07,1.3]
s2_prev=[0.3,0.5,0.1,1.2]
print(f"{'p':>6} {'P_S0-7':>10} {'P_S1-P_S2':>12} {'P_S2-7':>10}")
for pp in [0.24,0.26,0.28,0.30,0.32]:
    # S0
    try:
        s0=sp.nsolve([g.subs(p,sp.Float(pp)) for g in g0],list(xs),s0_prev,prec=35,tol=1e-30,maxsteps=200)
        s0_prev=[float(s) for s in s0]
        P0v=float(P0.subs(p,sp.Float(pp)).subs({xs[i]:s0[i] for i in range(6)}))
        I2=P0v-7
    except:I2=float('nan')
    # S1
    try:
        s1=sp.nsolve([g.subs(p,sp.Float(pp)) for g in g1],[b,c,d,e,f],s1_prev,prec=35,tol=1e-30,maxsteps=200)
        s1_prev=[float(s) for s in s1]
        P1v=float(P1.subs(p,sp.Float(pp)).subs({b:s1[0],c:s1[1],d:s1[2],e:s1[3],f:s1[4]}))
    except:P1v=float('nan')
    # S2
    try:
        s2=sp.nsolve([g.subs(p,sp.Float(pp)) for g in g2],[B,C,D,E],s2_prev,prec=35,tol=1e-30,maxsteps=200)
        s2_prev=[float(s) for s in s2]
        P2v=float(P2.subs(p,sp.Float(pp)).subs({B:s2[0],C:s2[1],D:s2[2],E:s2[3]}))
    except:P2v=float('nan')
    I1=P1v-P2v
    print(f"{pp:6.3f} {I2:10.5f} {I1:12.5f} {P2v-7:10.5f}")
