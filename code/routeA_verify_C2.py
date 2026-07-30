import sympy as sp
u,v,w,z=sp.symbols('u v w z')
a1=sp.Integer(1);a2=1-u;a3=1-v+u*v;a4=1-w+v*w-u*v*w;a5=1-z+z*w-z*v*w+z*u*v*w
bet=[u,v,w,z]
H=sp.ones(5,5)
for i in range(4):
    bi=bet[i];ai=[a1,a2,a3,a4][i]
    H[i,i]+=ai;H[i,i+1]+=ai*bi;H[i+1,i]+=ai*bi;H[i+1,i+1]+=ai*bi*(2*bi-1)
H[4,4]+=a5
D=sp.expand(sp.det(H))
HB=z*w**2+(1-z**2)*w+z**2-z
# parametrization
v_expr=z*(1-w)
u_expr=z*(1-z+z*w)/(1-z**2+z**2*w)
Dsub=sp.together(D.subs({v:v_expr,u:u_expr}))
num,den=sp.fraction(Dsub)
num=sp.expand(num);den=sp.expand(den)
print("denominator of D_B:",sp.factor(den))
# reduce numerator mod H_B (in w)
num_red=sp.rem(num,HB,w)
num_red=sp.expand(num_red)
print("num(D_B) mod H_B, nterms=",len(num_red.as_ordered_terms()))
print("num_red =",num_red)
# Try to see structure: collect in w. Since H_B is deg2 in w, num_red is linear in w: A(z)+w*B(z)
num_colw=sp.collect(num_red,w)
print("\ncollected in w:",num_colw)
