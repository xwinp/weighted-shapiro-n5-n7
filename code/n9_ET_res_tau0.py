# tau->0 branch only, longer timeout. eta=tau^2*z, Res_z(E,T) in (v,tau).
import sympy as sp
from sympy import Poly, symbols, factor, real_roots
v,xi,tau,eta,z=symbols('v xi tau eta z'); c,d,sg=symbols('c d sigma')
vv=c+d-1
Lc=vv+sg*c**2; Uc=vv-sg*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+sg*d**2; Ud=vv-sg*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - sg*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - sg*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(sp.expand(Fc),c,d,sg,domain=sp.QQ); Pd=Poly(sp.expand(Fd),c,d,sg,domain=sp.QQ)
FL=Pc.exquo(Poly((c-1)*Lc,c,d,sg,domain=sp.QQ)); FR=Pd.exquo(Poly((d-1)*Ld,c,d,sg,domain=sp.QQ))
Gm=Poly(sp.expand(FL.as_expr()-FR.as_expr()),c,d,sg,domain=sp.QQ).exquo(Poly(c-d,c,d,sg,domain=sp.QQ))
Gp=Poly(sp.expand(FL.as_expr()+FR.as_expr()),c,d,sg,domain=sp.QQ)
cc=(1+v+(1-v)*xi)/2; dd=(1+v-(1-v)*xi)/2; sig=tau*v/(dd*(1-dd))
Gmt=sp.expand(sp.fraction(sp.together(Gm.as_expr().subs({c:cc,d:dd,sg:sig})))[0])
Gpt=sp.expand(sp.fraction(sp.together(Gp.as_expr().subs({c:cc,d:dd,sg:sig})))[0])
F=Poly(Gmt,v,xi,tau,domain=sp.QQ); H=Poly(Gpt,v,xi,tau,domain=sp.QQ)
def v_order(P,var):
    m=None; idx=P.gens.index(var)
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if m is None or e<m: m=e
    return m if m is not None else 0
def shift_v(P,var,k):
    idx=P.gens.index(var); d={}
    for monom,coef in P.terms():
        if coef==0: continue
        e=monom[idx]
        if e>=k:
            nm=list(monom); nm[idx]=e-k; d[tuple(nm)]=d.get(tuple(nm),0)+coef
    return Poly(d,*P.gens,domain=sp.QQ)
F0=shift_v(F,v,v_order(F,v)); H0=shift_v(H,v,v_order(H,v))
E=Poly(sp.expand(F0.as_expr()),v,xi,tau,domain=sp.QQ)
U=Poly(sp.expand(2*F0.as_expr()-H0.as_expr()),v,xi,tau,domain=sp.QQ); T=shift_v(U,v,1)
# eta = tau^2 z  =>  xi = 1 - eta = 1 - tau^2 z
Ez=Poly(sp.expand(E.as_expr().subs({xi:1-tau**2*z})),v,tau,z,domain=sp.QQ)
Tz=Poly(sp.expand(T.as_expr().subs({xi:1-tau**2*z})),v,tau,z,domain=sp.QQ)
print("deg_z E=",Ez.degree(z)," deg_z T=",Tz.degree(z),flush=True)
# Use Poly.resultant (matrix-free, often faster) w.r.t. z
Rz=Poly(Ez.as_expr(),v,tau,z,domain=sp.QQ).resultant(Poly(Tz.as_expr(),v,tau,z,domain=sp.QQ),z)
RPz=Poly(sp.expand(Rz.as_expr()),v,tau,domain=sp.QQ)
print("resultant terms",len(RPz.as_dict()),"deg_v",RPz.degree(v),"deg_tau",RPz.degree(tau),flush=True)
# newton
iv=RPz.gens.index(v); it=RPz.gens.index(tau); byv={}
for monom,coef in RPz.terms():
    if coef==0: continue
    byv[monom[iv]]=min(byv.get(monom[iv],monom[it]),monom[it])
pts=sorted(byv.items())
hull=[]
for p in pts:
    while len(hull)>=2:
        a=hull[-2]; b=hull[-1]
        if (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])<=0: hull.pop()
        else: break
    hull.append(p)
print("lower hull (deg_v,deg_tau):",hull,flush=True)
for k in range(len(hull)-1):
    a=hull[k]; b=hull[k+1]
    sl=(b[1]-a[1])/(b[0]-a[0]) if b[0]!=a[0] else None
    di=b[0]-a[0]; dj=b[1]-a[1]; const=dj*a[0]-di*a[1]
    terms={}
    for monom,coef in RPz.terms():
        if coef==0: continue
        if dj*monom[iv]-di*monom[it]==const:
            terms[(monom[iv],monom[it])]=terms.get((monom[iv],monom[it]),0)+coef
    polyV=sp.expand(sum(coef*v**i for (i,j),coef in terms.items()))
    rts=[float(r) for r in real_roots(polyV,v) if r>0]
    print("edge",a,"->",b,"slope",sl,"ell(v~tau^ell)=",-sl if sl else sl,flush=True)
    print("  initial form:",sp.factor(polyV)," positive V roots:",rts,flush=True)
print("DONE",flush=True)
