# -*- coding: utf-8 -*-
"""tau->0 branch per GPT#11 exact prescription:
   eta = tau^2 * z,  solve T=0 for z = 4 + ... (series in v,tau),
   substitute into E, compute Newton polygon of E_red(v,tau) -> v ~ tau^ell.
   E = F0, T = (2F0-H0)/v.  Branch z->4 is simple root of leading transverse eq.
"""
import sympy as sp
from sympy import Poly, symbols, factor, Rational, real_roots
v,xi,tau,eta,z=symbols('v xi tau eta z')
c,d,sg=symbols('c d sigma')
vv=c+d-1
Lc=vv+sg*c**2; Uc=vv-sg*c*(1-c); Bc=c**2*vv+(1-c)*Lc**2
Ld=vv+sg*d**2; Ud=vv-sg*d*(1-d); Bd=d**2*vv+(1-d)*Ld**2
Fc=c*vv**2*(1-c)*Lc**2 - sg*Bc*(c*vv**2-Uc*Bc)
Fd=d*vv**2*(1-d)*Ld**2 - sg*Bd*(d*vv**2-Ud*Bd)
Pc=Poly(sp.expand(Fc),c,d,sg,domain=sp.QQ); Pd=sp.Poly(sp.expand(Fd),c,d,sg,domain=sp.QQ)
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
# eta = tau^2 z, xi = 1 - eta = 1 - tau^2 z
Ez=Poly(sp.expand(E.as_expr().subs({xi:1-tau**2*z})),v,tau,z,domain=sp.QQ)
Tz=Poly(sp.expand(T.as_expr().subs({xi:1-tau**2*z})),v,tau,z,domain=sp.QQ)
print("deg_z E=",Ez.degree(z),"deg_z T=",Tz.degree(z),flush=True)
# verify z=4 is root of leading transverse (v=0,tau=0): Tz(0,0,4)
T00=Poly(sp.expand(Tz.as_expr().subs({v:0,tau:0})),z,domain=sp.QQ)
print("T(0,0,z) =",factor(T00.as_expr()),flush=True)
print("  z=4 root?",T00.eval({z:4})," ; simple?",sp.gcd(T00,Poly(sp.diff(T00.as_expr(),z),z,domain=sp.QQ)).degree(),flush=True)

# ---- series z = 4 + delta, delta = sum a_{ij} v^i tau^j ----
# Truncation grid: v up to VMAX, tau up to TMAX.
VMAX=6; TMAX=14
# Build unknown coeffs a[i][j]
from sympy import symbols as syms
av=[[Rational(0)]*(TMAX+1) for _ in range(VMAX+1)]
adict={}
for i in range(VMAX+1):
    for j in range(TMAX+1):
        if i==0 and j==0: continue   # delta has no constant (z->4)
        s=syms('a%d_%d'%(i,j))
        av[i][j]=s; adict[(i,j)]=s
delta=sum(av[i][j]*v**i*tau**j for i in range(VMAX+1) for j in range(TMAX+1) if (i,j)!=(0,0))
zser=4+delta
# substitute z=zser into Tz, expand, truncate to (v^i tau^j) with i<=VMAX,j<=TMAX
Tsub=sp.expand(Tz.as_expr().subs({z:zser}))
# collect coeffs of v^i tau^j (treat a's as symbols). Use Poly in v,tau over QQ[a's].
# Simpler: iteratively substitute. Use series truncation by truncating monomials.
def trunc(expr,V,T):
    P=Poly(sp.expand(expr),v,tau,domain=sp.QQ)
    d={}
    for monom,coef in P.terms():
        if monom[0]<=V and monom[1]<=T:
            d[monom]=d.get(monom,0)+coef
    return Poly(d,v,tau,domain=sp.QQ)
Ttr=trunc(Tsub,VMAX,TMAX)
print("truncated T terms:",len(Ttr.as_dict()),flush=True)
# Now Ttr must ==0 identically (as series). Solve linearly for a's.
# Collect equations: coeff of each v^i tau^j (i<=VMAX,j<=TMAX) is linear in a's? 
# T is linear in z? No, T is degree 10 in z, so substitution gives nonlinear in a's.
# But near z=4, Tz is LINEAR in delta to leading order (simple root). The nonlinear
# terms are higher order. Use Newton/Hensel iteration: solve order-by-order.
# Strategy: build z series incrementally by total degree in (v,tau).
print("Solving z-series by Hensel/Newton order-by-order...",flush=True)
# Represent zser as a plain sympy expr in v,tau with rational coeffs, built up.
zcur=Rational(4)
# We solve Tz(z)=0. At each total degree D, add correction c*v^i*tau^j and determine c
# from the linear term: T(zcur)+ c*v^i*tau^j * dT/dz(zcur) = 0 (to that order).
for D in range(1,VMAX+TMAX+1):
    for i in range(VMAX+1):
        for j in range(TMAX+1):
            if i+j!=D: continue
            if i>VMAX or j>TMAX: continue
            # residual at current zcur truncated to order D
            res=trunc(sp.expand(Tz.as_expr().subs({z:zcur})),i,j)  # poly up to v^i tau^j
            # coefficient of v^i tau^j in residual:
            coef_res=res.nth(i,j) if (i<=res.degree(v) or True) else 0
            # need the (i,j) coeff; nth works
            try:
                coef_res=Poly(sp.expand(Tz.as_expr().subs({z:zcur})),v,tau,domain=sp.QQ).nth(i,j)
            except Exception:
                coef_res=Rational(0)
            if coef_res==0: 
                continue
            # dT/dz at zcur, coeff of v^0 tau^0 (leading) -> but need coeff at (i,j) offset 0
            dTdz=sp.diff(Tz.as_expr(),z).subs({z:zcur})
            try:
                lin=Poly(sp.expand(dTdz),v,tau,domain=sp.QQ).nth(0,0)
            except Exception:
                lin=Rational(0)
            if lin==0:
                print("  WARN linear coeff 0 at (%d,%d)"%(i,j),flush=True)
                continue
            c=Rational(-coef_res,lin)
            zcur=sp.expand(zcur+c*v**i*tau**j)
print("z series (truncated) =",zcur,flush=True)
# ---- substitute z(v,tau) into E, truncate, Newton polygon ----
Ered=trunc(sp.expand(Ez.as_expr().subs({z:zcur})),VMAX,TMAX)
print("E_red terms:",len(Ered.as_dict()),flush=True)
print("E_red(0,0) =",Ered.nth(0,0),flush=True)
# Newton polygon lower hull in (deg_v, deg_tau)
iv=0; it=1; byv={}
for monom,coef in Ered.terms():
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
    for monom,coef in Ered.terms():
        if coef==0: continue
        if dj*monom[iv]-di*monom[it]==const:
            terms[(monom[iv],monom[it])]=terms.get((monom[iv],monom[it]),0)+coef
    polyV=sp.expand(sum(coef*v**i for (i,j),coef in terms.items()))
    rts=[float(r) for r in real_roots(polyV,v) if r>0]
    print("edge",a,"->",b,"slope",sl," ell(v~tau^ell)=",-sl if sl else sl,flush=True)
    print("  initial form:",sp.factor(polyV)," positive V roots:",rts,flush=True)
print("DONE",flush=True)
