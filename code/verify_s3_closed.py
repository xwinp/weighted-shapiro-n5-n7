#!/usr/bin/env python3
"""Verify S3={0,3} closed form AND the EXACT positivity certificate used in the
paper (Sec 4.1).

The closed form: with r=(q/p)^{1/5}, the KKT solution of S3={0,3} is
    x=(0,1,b,0,c,d,e),  b=r^2, c=r^{-1}, e=r^3, d=r(1-r^7),
    p=1/(1+r^5), q=r^5/(1+r^5),  and d>0 <=> r<1 <=> p>1/2.
The stationary value is  P_S3^stat(r) = (1+r^5)(5-r^7)/r^2.

EXACT certificate (replaces the earlier non-rigorous N(0)/N(1) endpoint scan --
"one interior critical point and two positive endpoints" does NOT imply positivity
throughout; the scan was only suggestive).  The paper instead proves the STRICTER
bound P_S3^stat > 8 > 7 via the identity

    (1+r^5)(5-r^7) - 8 r^2 = (1-r) H(r),

    H(r) = r^11 + r^10 + r^9 + r^8 + r^7 + 2 r^6 + 2 r^5
             - 3 r^4 - 3 r^3 - 3 r^2 + 5 r + 5,

and the exact positivity chain on 0<r<1:
    H(r) - (5+5r-9r^2) = r^2 * B(r),
        B(r) = r^9+r^8+r^7+r^6+r^5+2r^4+2r^3-3r^2-3r+6,
    B(r) = 3(2-r^2-r) + (r^9+r^8+r^7+r^6+r^5+2r^4+2r^3),
    2-r^2-r = (1-r)(r+2) > 0  and  the residual sum >= 0  =>  B(r) > 0,
    5+5r-9r^2 = 1 + (1-r)(9r+4) > 1                                     =>  H(r) > 0,
    (1-r)/r^2 > 0                                                       =>  P_S3^stat - 8 > 0.

Every step is verified to be the ZERO polynomial (sympy expand/simplify == 0),
so the certificate is exact, not numerical.  exit 0 iff all identities hold and
the strict-positivity factors are positive on (0,1).
"""
import sys
import sympy as sp

p, q = sp.symbols('p q', positive=True)
b, c, d, e = sp.symbols('b c d e', positive=True)
r = sp.symbols('r', positive=True)

# ---- 1. KKT: build P on S3={0,3}, verify the closed form zeroes the gradient ----
# S3 = {0,3}: x0=0,x1=1,x2=b,x3=0,x4=c,x5=d,x6=e
x = {0: 0, 1: sp.Integer(1), 2: b, 3: 0, 4: c, 5: d, 6: e}
P = sp.simplify(sum(x[i] / (p * x[(i + 1) % 7] + q * x[(i + 2) % 7]) for i in range(7)))
grads = {v: sp.simplify(sp.diff(P, v)) for v in [b, c, d, e]}

subs = [(b, r**2), (c, 1 / r), (e, r**3), (d, r - r**8),
        (p, 1 / (1 + r**5)), (q, r**5 / (1 + r**5))]
kkt_ok = True
print("-- KKT: gradients on the support vanish at the closed form --")
for v, g in grads.items():
    gv = sp.simplify(g.subs(subs))
    is_zero = sp.simplify(gv) == 0
    kkt_ok = kkt_ok and is_zero
    print(f"  dP/d{v} = {gv}   (zero: {is_zero})")

Pv = sp.simplify(P.subs(subs))
target = (1 + r**5) * (5 - r**7) / r**2
form_ok = sp.simplify(Pv - target) == 0
print(f"\nP_S3^stat(r) = {sp.factor(Pv)}")
print(f"matches (1+r^5)(5-r^7)/r^2 ? {form_ok}")

# ---- 2. EXACT identity:  (1+r^5)(5-r^7) - 8 r^2 == (1-r) H(r) ----
H = (r**11 + r**10 + r**9 + r**8 + r**7 + 2 * r**6 + 2 * r**5
     - 3 * r**4 - 3 * r**3 - 3 * r**2 + 5 * r + 5)
lhs = sp.expand((1 + r**5) * (5 - r**7) - 8 * r**2)
rhs = sp.expand((1 - r) * H)
id_ok = sp.simplify(lhs - rhs) == 0
print(f"\n-- EXACT identity (1+r^5)(5-r^7)-8r^2 = (1-r)H(r) : {id_ok}")

# ---- 3. EXACT positivity chain (each link == 0 polynomial) ----
# 3a. H(r) - (5+5r-9r^2) == r^2 * B(r)
B = (r**9 + r**8 + r**7 + r**6 + r**5 + 2 * r**4 + 2 * r**3
     - 3 * r**2 - 3 * r + 6)
link_a = sp.simplify(sp.expand(H - (5 + 5 * r - 9 * r**2)) - sp.expand(r**2 * B)) == 0
print(f"  H - (5+5r-9r^2) = r^2 * B(r)              : {link_a}")

# 3b. B(r) == 3(2-r^2-r) + residual,  residual = r^9+...+2r^3 (>=0 on (0,1))
residual = sp.expand(r**9 + r**8 + r**7 + r**6 + r**5 + 2 * r**4 + 2 * r**3)
link_b = sp.simplify(sp.expand(B - (3 * (2 - r**2 - r) + residual))) == 0
print(f"  B = 3(2-r^2-r) + (r^9+...+2r^3)           : {link_b}")
# residual >= 0 for r>0 (every term has positive coeff, positive power)
res_terms_pos = all(coef > 0 for coef in sp.Poly(residual, r).coeffs())
print(f"  residual terms all positive coeff         : {res_terms_pos}")

# 3c. 2-r^2-r == (1-r)(r+2)  > 0 on (0,1)
link_c = sp.simplify(sp.expand(2 - r**2 - r) - sp.expand((1 - r) * (r + 2))) == 0
print(f"  2-r^2-r = (1-r)(r+2)                      : {link_c}")

# 3d. 5+5r-9r^2 == 1 + (1-r)(9r+4) > 1 on (0,1)
link_d = sp.simplify(sp.expand(5 + 5 * r - 9 * r**2) - sp.expand(1 + (1 - r) * (9 * r + 4))) == 0
print(f"  5+5r-9r^2 = 1 + (1-r)(9r+4)               : {link_d}")

# 3e. (1-r)/r^2 > 0 on (0,1)
print(f"  (1-r)/r^2 > 0 on (0,1)                    : True  (1-r>0, r^2>0)")

chain_ok = id_ok and link_a and link_b and link_c and link_d and res_terms_pos
print(f"\nEXACT chain holds => H(r)>0 => P_S3^stat - 8 = (1-r)H(r)/r^2 > 0 : {chain_ok}")
print(f"=> P_S3^stat(r) > 8 > 7  for all 0<r<1 (i.e. p>1/2)              : {chain_ok}")

# ---- 4. The failure band has NO positive S3 stationary point (p<=1/2) ----
# d = r(1-r^7); d>0 <=> r<1 <=> q<p <=> p>1/2.  Band (a7,b7) subset (0,1/3) subset (0,1/2).
a7, b7 = sp.Rational('0.214273520909841').limit_denominator(10**12), \
         sp.Rational('0.328627677916592').limit_denominator(10**12)
band_below_half = b7 < sp.Rational(1, 2)
print(f"\n-- failure band (a7,b7) lies below p=1/2: b7={float(b7):.6f} < 0.5 : {band_below_half}")
print("   => for p in (a7,b7): p<1/2 => r>1 => d=r(1-r^7)<0 => no positive S3 stat pt;")
print("      inf_S3 P attained on dS3 subset S4 where P>7 (Sec 4.1).  S3 never fails.")

# ---- 5. UNIQUENESS of the positive KKT point (b=r^2 is the only positive root) ----
# Sequential linear-fractional elimination (each gradient is affine in one variable):
#   gb=0 (affine in 1/c) -> c = b^2/r^5
#   gc=0 & gd=0          -> e = b^4/r^5,  d = b^3(1-b r^5)/r^5
#   ge=0 under these     -> -(b^5-r^10)(r^5+1) = 0, i.e. (b-r^2)*quartic*(r^5+1)=0.
# Recompute P and gradients with sp.together (form needed for the exact substitution).
# Substitute p=1/(1+r^5), q=r^5/(1+r^5) so the gradients are functions of (b,c,d,e,r) only.
P_t = sp.together(sum(x[i] / (p * x[(i + 1) % 7] + q * x[(i + 2) % 7]) for i in range(7)))
P_t = sp.together(P_t.subs({p: 1 / (1 + r**5), q: r**5 / (1 + r**5)}))
g_t = {v: sp.together(sp.diff(P_t, v)) for v in [b, c, d, e]}
def num_of(expr):
    return sp.expand(sp.together(expr).as_numer_denom()[0])
gb_n, gc_n, gd_n, ge_n = [num_of(g_t[v]) for v in [b, c, d, e]]
# gb = 1/(q c) - 1/(p b^2) (independent of d,e), strictly decreasing in c
# => gb=0 has a UNIQUE positive solution c = p b^2 / q = b^2/r^5.
p_e = 1 / (1 + r**5); q_e = r**5 / (1 + r**5)
gb_form = sp.simplify(g_t[b] - (1 / (q_e * c) - 1 / (p_e * b**2))) == 0
gb_decr = sp.simplify(sp.diff(g_t[b], c) + 1 / (q_e * c**2)) == 0   # dgb/dc = -1/(q c^2) < 0
c_expr = b**2 / r**5
e_expr = b**4 / r**5
d_expr = b**3 * (1 - b * r**5) / r**5
gb_zero = num_of(g_t[b].subs(c, c_expr)) == 0
# with c = b^2/r^5, the gradients gc, gd factor EXPLICITLY (no sp.solve on the system):
#   gc_ce = -e (r+1) (cyclot) (d r^5 + e r^10 - b^3)        -> eq1: d r^5 + e r^10 = b^3
#   gd_ce =  b (r+1) (cyclot) (r^5 (d + e r^5)^2 - b^2 e)   -> eq2: r^5 (d + e r^5)^2 = b^2 e
# (prefactors e, b, r+1, cyclot=(r^5+1)/(r+1) are all > 0 on the positive domain), so the KKT
# conditions gc=0 & gd=0 are EQUIVALENT to {eq1=0, eq2=0}.  eq1 is linear in d and eq2 linear in
# e (given d), so the two pin (e,d) UNIQUELY.
gc_ce = num_of(g_t[c].subs(c, c_expr)); gd_ce = num_of(g_t[d].subs(c, c_expr))
cyclot = r**4 - r**3 + r**2 - r + 1          # = (r^5+1)/(r+1) > 0 for r>0
eq1 = sp.expand(d*r**5 + e*r**10 - b**3)
eq2 = sp.expand(r**5*(d + e*r**5)**2 - b**2*e)
gc_factor = sp.expand(gc_ce - (-e*(r + 1)*cyclot*eq1)) == 0
gd_factor = sp.expand(gd_ce - ( b*(r + 1)*cyclot*eq2)) == 0
ed_sol = sp.solve([eq1, eq2], [e, d], dict=True)
ed_unique = (len(ed_sol) == 1 and sp.simplify(ed_sol[0][e] - e_expr) == 0
             and sp.simplify(ed_sol[0][d] - d_expr) == 0)
gc_zero = num_of(g_t[c].subs([(c, c_expr), (e, e_expr), (d, d_expr)])) == 0
gd_zero = num_of(g_t[d].subs([(c, c_expr), (e, e_expr), (d, d_expr)])) == 0
ge_num = num_of(g_t[e].subs([(c, c_expr), (e, e_expr), (d, d_expr)]))
quartic = b**4 + b**3 * r**2 + b**2 * r**4 + b * r**6 + r**8
target = sp.expand(-(-b + r**2) * (r + 1) * quartic * cyclot)
uniq_factor = sp.expand(ge_num - target) == 0
quartic_pos = all(co > 0 for co in sp.Poly(quartic, b, r).coeffs())
cyclot_pos = sp.expand((r + 1) * cyclot - (r**5 + 1)) == 0  # (r+1)(r^4-...+1)=r^5+1>0
# b^5 - r^10 == (b-r^2)*quartic  (so b=r^2 is the unique positive b^5=r^10 root)
b5_factor = sp.expand((b - r**2) * quartic - (b**5 - r**10)) == 0
print("\n-- UNIQUENESS of the positive S3 KKT point --")
print(f"  gb = 1/(qc)-1/(pb^2), strictly decreasing in c (unique c) : {gb_form and gb_decr}")
print(f"  gb=0 -> c=b^2/r^5                                          : {gb_zero}")
print(f"  gc factors: gc = -e(r+1)(cyclot)*(d r^5+e r^10-b^3)       : {gc_factor}")
print(f"  gd factors: gd =  b(r+1)(cyclot)*(r^5(d+e r^5)^2-b^2 e)   : {gd_factor}")
print(f"  eq1,eq2 -> UNIQUE (e,d)=(b^4/r^5, b^3(1-br^5)/r^5)        : {ed_unique}")
print(f"  ge=0 -> -(b-r^2)(r+1)(quartic)(r^4-r^3+r^2-r+1) [exact]   : {uniq_factor}")
print(f"  b^5-r^10 = (b-r^2)*quartic                                 : {b5_factor}")
print(f"  quartic coeffs all positive (=> >0 for b,r>0)             : {quartic_pos}")
print(f"  (r+1)(r^4-r^3+r^2-r+1) = r^5+1 > 0                        : {cyclot_pos}")
print("  => b=r^2 is the UNIQUE positive root; c=r^-1, e=r^3, d=r(1-r^7) unique.")

uniq_ok = (gb_form and gb_decr and gb_zero and gc_factor and gd_factor and ed_unique
           and gc_zero and gd_zero and uniq_factor
           and quartic_pos and cyclot_pos and b5_factor)

ok = kkt_ok and form_ok and chain_ok and band_below_half and uniq_ok
print(f"\nALL CHECKS PASS (KKT + closed form + exact H(r) chain + band + UNIQUENESS): {ok}")
sys.exit(0 if ok else 1)
