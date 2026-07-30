import sympy as sp
u,v,w,z=sp.symbols('u v w z')
a1=1
a2=1-u
a3=1-v+u*v
a4=1-w+v*w-u*v*w
a5=1-z+z*w-z*v*w+z*u*v*w
E2=u*(1-z)-z*a5*(1-v)
E3=a3*v-u*(1-w)
E4=a4*w-z*a5*(1-w)
# Branch: v=z(1-w) (verified numerically). Substitute, then E3 gives u.
v_expr=z*(1-w)
E3v=sp.simplify(E3.subs(v,v_expr))
# E3=0: solve for u
u_expr=sp.solve(E3v,u)[0]
u_expr=sp.simplify(u_expr)
print("v =",v_expr)
print("u (from E3, v=z(1-w)) =",u_expr)
# Substitute v,u into E2,E4, reduce mod H_B
HB=z*w**2+(1-z**2)*w+z**2-z
E2red=sp.simplify(sp.rem(sp.numer(sp.together(E2.subs({v:v_expr,u:u_expr}))),HB,w))
E4red=sp.simplify(sp.rem(sp.numer(sp.together(E4.subs({v:v_expr,u:u_expr}))),HB,w))
E3red=sp.simplify(sp.rem(sp.numer(sp.together(E3.subs({v:v_expr,u:u_expr}))),HB,w))
print("E2 mod H_B =",E2red," (should be 0)")
print("E3 mod H_B =",E3red," (should be 0)")
print("E4 mod H_B =",E4red," (should be 0)")
# Now verify branch decomposition: eliminate u from E2,E3,E4, then resultant in v
# Solve E3 for u (general), substitute into E2,E4
u_gen=sp.solve(E3,u)[0]
E2u=sp.simplify(E2.subs(u,u_gen))
E4u=sp.simplify(E4.subs(u,u_gen))
# resultant in v of numerators
n2=sp.numer(sp.together(E2u));n4=sp.numer(sp.together(E4u))
print("\ncomputing Res_v(E2,E4) (after u eliminated via E3)...")
Res=sp.resultant(sp.expand(n2),sp.expand(n4),v)
Res=sp.factor(Res)
print("Res_v =",Res)
