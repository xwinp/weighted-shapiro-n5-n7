import sympy as sp
w,z,u,v=sp.symbols('w z u v')
HC=z*w**3+w**2*z**3-w**2*z+w*z**4-3*w*z**3+2*w*z**2+w*z-w-z**4+3*z**3-3*z**2+z
a7,b7=0.214273520909841,0.328627677916592
for zv in [0.85,0.9,0.95]:
    HCz=sp.Poly(HC.subs(z,zv),w)
    print(f"\nz={zv}: H_C real roots in (0,1):")
    for wr in sp.nroots(HCz,n=15):
        if abs(sp.im(wr))<1e-8:
            wrf=float(sp.re(wr))
            tag="(0,1)" if 0<wrf<1 else ""
            print(f"   w={wrf:.6f} {tag}")
    # for each w in (0,1), solve E2,E3 over a grid
    HCz_e=HC.subs(z,zv)
    wroots=[float(sp.re(r)) for r in sp.nroots(HCz_e,n=15) if abs(sp.im(r))<1e-8 and 0<float(sp.re(r))<1]
    for wrf in wroots:
        a3g=1-v+u*v;a5g=1-z+z*w-z*v*w+z*u*v*w
        E2=u*(1-z)-z*a5g*(1-v);E3=a3g*v-u*(1-w)
        E2z=E2.subs({z:zv,w:wrf});E3z=E3.subs({z:zv,w:wrf})
        # solve E3 for u: u=(v-v^2)/(1-w-v^2) ; substitute into E2
        usol=sp.solve(E3z,u)[0]
        E2u=sp.simplify(E2z.subs(u,usol))
        vroots=[float(sp.re(r)) for r in sp.nroots(sp.Poly(sp.together(E2u).as_numer_denom()[0],v),n=15) if abs(sp.im(r))<1e-8 and 0<float(sp.re(r))<1]
        for vv in vroots:
            uu=float(usol.subs(v,vv))
            if 0<uu<1:
                a5v=1-zv+zv*wrf-zv*vv*wrf+zv*uu*vv*wrf
                K=uu*vv*wrf*(zv**3)*a5v**2/((1-vv)*(1-wrf)*(1-zv)**3)
                pval=1/(1+K**(1/7))
                print(f"   w={wrf:.5f} v={vv:.5f} u={uu:.5f} K={K:.4e} p={pval:.5f} {'<<IN BAND' if a7<pval<b7 else ''}")
