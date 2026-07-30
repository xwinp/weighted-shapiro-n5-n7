import sympy as sp
w,z=sp.symbols('w z')
P=-8*z**16+100*z**15-528*z**14+1564*z**13-2972*z**12+4112*z**11-4806*z**10+5120*z**9-4696*z**8+3540*z**7-2298*z**6+1302*z**5-576*z**4+188*z**3-52*z**2+10*z
Q=8*z**16-100*z**15+536*z**14-1648*z**13+3340*z**12-5016*z**11+6270*z**10-6968*z**9+6786*z**8-5644*z**7+4056*z**6-2536*z**5+1332*z**4-568*z**3+194*z**2-48*z+6
num_red=P*w+Q
HB=z*w**2+(1-z**2)*w+z**2-z
Res=sp.resultant(sp.expand(HB),sp.expand(num_red),w)
Res=sp.factor(Res)
print("Res_w(H_B, num_red) factored:")
print(Res)
# also check sign of num_red at z=0.977,w=0.1296
zv=0.977242;wv=0.129580
nv=float(num_red.subs({z:zv,w:wv}))
dv=float((1-zv**2+zv**2*wv)**5)
print(f"\nat z={zv},w={wv}: num_red={nv:.6e}, denom^5={dv:.6e}, D_B={nv/dv:.6e} -> {'<0 CERT' if nv/dv<0 else '>0'}")
