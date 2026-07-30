#!/usr/bin/env python3
"""
Exhaustive exact verification of the GPT-generated paper
'Explicit Odd-Dimensional Counterexamples to a Weighted Shapiro Cyclic Inequality'.
All rational arithmetic via fractions.Fraction. sympy for interval algebra.
"""
from fractions import Fraction as F
import sympy as sp

p = F(3, 10)
q = F(7, 10)

def P(x, p=p, q=q):
    n = len(x)
    s = F(0)
    for i in range(n):
        d = p * x[(i + 1) % n] + q * x[(i + 2) % n]
        assert d > 0, f"nonpositive denominator at i={i}: {d}"
        s += F(x[i]) / d
    return s

ok = True
def chk(name, cond, val=""):
    global ok
    ok = ok and cond
    print(f"[{'OK ' if cond else 'FAIL'}] {name}  {val}")

# ---------- §2.1 Proposition 2.1: n=7 ----------
x7 = (0, 3, 4, 1, 5, 0, 4)
P7 = P(x7)
chk("n=7 P == 1859/266", P7 == F(1859, 266), f"P={P7}")
chk("n=7 deficit == 3/266", 7 - P7 == F(3, 266), f"deficit={7-P7}")
# seven denominators
dens7 = [p*x7[(i+1)%7]+q*x7[(i+2)%7] for i in range(7)]
chk("n=7 denominators match paper",
    dens7 == [F(37,10),F(19,10),F(19,5),F(3,2),F(14,5),F(6,5),F(21,10)],
    str([str(d) for d in dens7]))

# ---------- §2.1 Proposition 2.2: n=9 ----------
x9 = (0,1,0,1,0,1,0,1,2)
P9 = P(x9)
chk("n=9 P == 185/21", P9 == F(185, 21), f"P={P9}")
chk("n=9 deficit == 4/21", 9 - P9 == F(4, 21), f"deficit={9-P9}")

# ---------- §2.2 Theorem 2.3: closed-form family ----------
def family(m):
    # n=2m+1, x=(0,1,0,1,...,0,1,1)
    n = 2*m+1
    x = [0]*n
    for j in range(1, m+1):
        x[2*j-1] = 1
    x[2*m] = 1
    return tuple(x)

for m in range(5, 30):
    x = family(m)
    Pfam = P(x)
    chk(f"family m={m} P==(30m+70)/21", Pfam == F(30*m+70, 21), f"P={Pfam}")
    chk(f"family m={m} deficit==(12m-49)/21>0",
        (2*m+1) - Pfam == F(12*m-49, 21) and F(12*m-49,21) > 0,
        f"deficit={(2*m+1)-Pfam}")

# ---------- §3 table n=4..13 ----------
table = {
    4: ((1,2,1,2), F(860,221)),
    6: ((1,2,1,2,1,2), F(1290,221)),
    7: (x7, F(1859,266)),
    8: ((1,2,1,2,1,2,1,2), F(1720,221)),
    9: (x9, F(185,21)),
    10: ((1,2)*5, F(2150,221)),
    11: (family(5), F(220,21)),
    12: ((1,2)*6, F(2580,221)),
    13: (family(6), F(250,21)),
}
for n,(x,exp) in table.items():
    Pv = P(x)
    chk(f"table n={n} P=={exp}", Pv == exp, f"P={Pv}")

# ---------- §3 even-n witness: P=215n/221, deficit 6n/221 ----------
for n in range(4, 30, 2):
    x = (1,2)*(n//2)
    Pv = P(x)
    chk(f"even n={n} P==215n/221", Pv == F(215*n,221), f"P={Pv}")
    chk(f"even n={n} deficit==6n/221", n-Pv == F(6*n,221), f"def={n-Pv}")

# ---------- §4 general-p analysis ----------
m, psp = sp.symbols('m p', real=True)
# P(p) = m/(1-p) + 1/p
Pgen = m/(1-psp) + 1/psp
# P<2m+1  <=> (2m+1)p^2-(m+2)p+1<0
ineq = sp.simplify(Pgen - (2*m+1))
# bring to common form
F_expr = sp.simplify(sp.together(ineq) * psp * (1-psp))  # numerator over p(1-p)
chk("§4 (2m+1)p^2-(m+2)p+1 is the numerator of P-(2m+1)",
    sp.expand(F_expr) == sp.expand((2*m+1)*psp**2 - (m+2)*psp + 1),
    f"got {sp.expand(F_expr)}")

# discriminant = m(m-4)
disc = sp.expand((m+2)**2 - 4*(2*m+1))
chk("§4 discriminant == m(m-4)", disc == sp.expand(m*(m-4)), f"disc={disc}")

# p_\pm(m)
p_minus = (m+2 - sp.sqrt(m*(m-4)))/(4*m+2)
p_plus  = (m+2 + sp.sqrt(m*(m-4)))/(4*m+2)
# verify sample intervals
import sympy
samples = {11:(0.216542,0.419821),13:(0.174458,0.440927),15:(0.147247,0.452753),
           17:(0.127740,0.460496),21:(0.101287,0.470142)}
for n,(lo,hi) in samples.items():
    mm = (n-1)//2
    pl = float(p_minus.subs(m,mm)); ph = float(p_plus.subs(m,mm))
    chk(f"§4 interval n={n} low~{lo}", abs(pl-lo)<2e-6, f"p_-={pl}")
    chk(f"§4 interval n={n} high~{hi}", abs(ph-hi)<2e-6, f"p_+={ph}")

# F_{m+1}-F_m = p(2p-1)
Fm = (2*m+1)*psp**2 - (m+2)*psp + 1
diff = sp.expand(Fm.subs(m,m+1) - Fm)
chk("§4 F_{m+1}-F_m == p(2p-1)", diff == sp.expand(psp*(2*psp-1)), f"diff={diff}")

# F_m = m*p*(2p-1) + (1-p)^2  (rewritten form)
Fm_alt = m*psp*(2*psp-1) + (1-psp)**2
chk("§4 F_m == m*p*(2p-1)+(1-p)^2", sp.expand(Fm) == sp.expand(Fm_alt))

# Corollary 4.1 threshold: F_m<0  <=>  m > (1-p)^2/(p(1-2p))  for 0<p<1/2
threshold = (1-psp)**2 / (psp*(1-2*psp))
# verify: F_m<0  <=>  m>threshold
# F_m = -m*p(1-2p) + (1-p)^2 ; <0 <=> m p(1-2p) > (1-p)^2 <=> m > thr (since p(1-2p)>0)
import random
rng = random.Random(1)
for _ in range(200):
    pp = F(rng.randint(1,4999), 10000)  # 0<p<1/2
    mm = rng.randint(5, 500)
    Fval = (2*mm+1)*pp**2 - (mm+2)*pp + 1
    thr = (1-pp)**2 / (pp*(1-2*pp))
    cond_F = Fval < 0
    cond_thr = mm > thr
    if cond_F != cond_thr:
        chk("Cor4.1 threshold equivalence", False, f"p={pp},m={mm},F={Fval},thr={thr}")
        break
else:
    chk("Cor4.1 threshold equivalence (200 random)", True)

# F_m(3/10) = (49-12m)/100
for mm in range(5, 30):
    Fval = (2*mm+1)*F(3,10)**2 - (mm+2)*F(3,10) + 1
    chk(f"§4 F_m(0.3)==(49-12m)/100 m={mm}", Fval == F(49-12*mm, 100), f"F={Fval}")

# ---------- §5 strictly positive n=7 variant ----------
x7pos = (1,300,400,100,500,1,400)
Ppos = P(x7pos)
chk("§5 strictly-positive n=7 P == 527034044794263379/75376890537112230",
    Ppos == F(527034044794263379, 75376890537112230), f"P={Ppos} ~{float(Ppos):.9f}")
chk("§5 strictly-positive n=7 P<7", Ppos < 7)

print()
print("ALL OK" if ok else "SOME FAILURES")
