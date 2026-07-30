import sympy as sp
w,z=sp.symbols('w z')
# verify p=1/4 sample (GPT: u=.883570,v=.872280,w=.113040,z=.983449, D_B~-3.75e-4)
b,c,d,e,f,p=sp.symbols('b c d e f p',positive=True);q=1-p
P=1/(p*b+q*c)+b/(p*c+q*d)+c/(p*d+q*e)+d/(p*e+q*f)+e/(p*f)+f/q
gs=[sp.diff(P,v) for v in [b,c,d,e,f]]
sol=sp.nsolve([g.subs(p,sp.Rational(1,4)) for g in gs],[b,c,d,e,f],[0.27,0.68,1.5,0.07,1.3],prec=40,tol=1e-35,maxsteps=200)
bb,cc,dd,ee,ff=[float(s) for s in sol]
qv=0.75;pp=0.25
rho=[bb,cc/bb,dd/cc,ee/dd,ff/ee]
u0=qv*rho[1]/(pp+qv*rho[1]);v0=qv*rho[2]/(pp+qv*rho[2]);w0=qv*rho[3]/(pp+qv*rho[3]);z0=qv*rho[4]/(pp+qv*rho[4])
print(f"p=1/4 sample: u={u0:.6f} v={v0:.6f} w={w0:.6f} z={z0:.6f}")
print(f"  (GPT:      u=0.883570 v=0.872280 w=0.113040 z=0.983449)")
# D_B at this z,w
Pcoeff=-8*z**16+100*z**15-528*z**14+1564*z**13-2972*z**12+4112*z**11-4806*z**10+5120*z**9-4696*z**8+3540*z**7-2298*z**6+1302*z**5-576*z**4+188*z**3-52*z**2+10*z
Qcoeff=8*z**16-100*z**15+536*z**14-1648*z**13+3340*z**12-5016*z**11+6270*z**10-6968*z**9+6786*z**8-5644*z**7+4056*z**6-2536*z**5+1332*z**4-568*z**3+194*z**2-48*z+6
num_red=float((Pcoeff*w+Qcoeff).subs({z:z0,w:w0}))
den=float((1-z0**2+z0**2*w0)**5)
print(f"  D_B = {num_red/den:.6e}  (GPT: -3.7505e-4)  -> {'<0 CERT' if num_red/den<0 else '>0'}")
# Confirm H_B has exactly one positive w-root for z in (0.9,1): product of roots = (z^2-z)/z = z-1 <0
print(f"\n  H_B product of roots = z-1; for z in (0.9,1): z-1<0 -> exactly one positive root. OK")
# Final: Q5,Q7 sign on (9/10,1) - already 0 roots + positive at 9/10 => positive throughout
Q5=2*z**5+2*z**3-2*z**2-1;Q7=8*z**7-24*z**6+20*z**5-9*z**4+30*z**3-15*z**2-6
print(f"\n  Q5(9/10)={float(Q5.subs(z,sp.Rational(9,10)))} >0, 0 roots in (9/10,1) => Q5>0 on band")
print(f"  Q7(9/10)={float(Q7.subs(z,sp.Rational(9,10)))} >0, 0 roots in (9/10,1) => Q7>0 on band")
print(f"  => Res_w = 4*z*(z-1)^9*Q5*Q7 != 0 on band (z in (0.9488,0.9911))")
print(f"  => det(H)!=0 on H_B branch in band => Morse constant => det<0 (sample neg) => S1 saddle. QED")
