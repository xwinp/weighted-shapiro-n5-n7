import sympy as sp
u,v,w,z=sp.symbols('u v w z')
a1=sp.Integer(1);a2=1-u;a3=1-v+u*v;a4=1-w+v*w-u*v*w;a5=1-z+z*w-z*v*w+z*u*v*w
# Normalized Hessian H = 11^T + sum_{i=1}^4 a_i B_i + a5 e5e5^T
# B_i at (i,i+1): [[1,beta_i],[beta_i,beta_i(2beta_i-1)]], beta=(u,v,w,z) for i=1..4
bet=[u,v,w,z]
H=sp.ones(5,5)
for i in range(4):
    bi=bet[i];ai=[a1,a2,a3,a4][i]
    H[i,i]+=ai*1
    H[i,i+1]+=ai*bi
    H[i+1,i]+=ai*bi
    H[i+1,i+1]+=ai*bi*(2*bi-1)
H[4,4]+=a5
print("H[0,0]=",sp.simplify(H[0,0]),"H[1,1]=",sp.simplify(H[1,1]))
D=sp.det(H)
D=sp.expand(D)
print("det(H) computed, total degree check, nterms(approx)=",len(D.as_ordered_terms()))
# factor 2(1-u)^2 ?
F=sp.factor(D)
print("det(H) factored:",F)
