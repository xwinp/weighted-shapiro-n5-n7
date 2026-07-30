import sympy as sp
import numpy as np
# get high-precision S1 stationary point at p=1/4 and p=0.27
b,c,d,e,f,p=sp.symbols('b c d e f p',positive=True)
q=1-p
P=1/(p*b+q*c)+b/(p*c+q*d)+c/(p*d+q*e)+d/(p*e+q*f)+e/(p*f)+f/q
vs=[b,c,d,e,f]
gs=[sp.diff(P,v) for v in vs]
Hs=sp.Matrix(5,5,lambda i,j: sp.diff(P,vs[i],vs[j]))
def get_sol(pp,init):
    gsub=[g.subs(p,sp.Float(pp)) for g in gs]
    sol=sp.nsolve(gsub,vs,init,prec=50,tol=1e-45,maxsteps=300)
    return [sp.Float(sol[i],50) for i in range(5)]
sol27=get_sol(0.27,[0.27,0.68,1.5,0.07,1.3])
sol25=get_sol(0.25,[0.28996,0.74475,0.88758,0.37644,1.16652,0.15912][:5])
for label,sol,pp in [("p=0.27",sol27,0.27),("p=0.25",sol25,0.25)]:
    bb,cc,dd,ee,ff=sol
    rho=[bb, cc/bb, dd/cc, ee/dd, ff/ee]
    ppf=sp.Float(pp);qv=1-ppf
    beta=[qv*rho[i+1]/(ppf+qv*rho[i+1]) for i in range(4)]
    u,v,w,z=[float(x) for x in beta]
    a1=1.0;a2=1-u;a3=1-v+u*v;a4=1-w+v*w-u*v*w;a5=1-z+z*w-z*v*w+z*u*v*w
    E2=u*(1-z)-z*a5*(1-v)
    E3=a3*v-u*(1-w)
    E4=a4*w-z*a5*(1-w)
    K=num=u*v*w*(z**3)*a5**2/((1-v)*(1-w)*(1-z)**3)
    Ktarget=(qv/ppf)**7
    print(f"\n{label}: u={u:.6f} v={v:.6f} w={w:.6f} z={z:.6f}")
    print(f"  a5={a5:.6f}")
    print(f"  E2={E2:.3e} E3={E3:.3e} E4={E4:.3e}  (should ~0)")
    print(f"  K={float(K):.6e}  (q/p)^7={float(Ktarget):.6e}  ratio={float(K/Ktarget):.10f}")
    # H_B parametrization check
    HB = z*w**2+(1-z**2)*w+z**2-z
    print(f"  H_B(w,z)={HB:.3e} (should ~0)")
    wp = (2*z)/(z**2-1) + sp.sqrt((1-z)*(1+z+3*z**2-z**3))/( ... )  # placeholder
