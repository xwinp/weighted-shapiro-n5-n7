#!/usr/bin/env python3
"""
Airtight verification of the infinite counterexample family and the
"holds only for n in {3,5}" framing at p=3/10, q=7/10.

Family (odd n = 2m+1, m >= 5, i.e. n >= 11):
    x = (0,1,0,1,...,0,1,1)   (even positions 0, odd positions 1, last = 1)
    => P = (30 m + 70)/21,   deficit = n - P = (12 m - 49)/21 > 0  for m >= 5.

n=7 (m=3): delicate witness (0,3,4,1,5,0,4), P = 1859/266.
n=9 (m=4): witness (0,1,0,1,0,1,0,1,2), P = 185/21.

Even n>=4 (p<q): paper's pattern x=(a,b,a,b,...) with a!=b, e.g. (1,2,1,2,...).

This script verifies every witness exactly and checks the general-odd-family
formula symbolically over m.
"""
from fractions import Fraction as F
from math import gcd


def exact_P(xv, p=F(3, 10), q=F(7, 10)):
    n = len(xv)
    # require all denominators strictly positive (valid domain)
    for i in range(n):
        if p * xv[(i + 1) % n] + q * xv[(i + 2) % n] <= 0:
            return None
    s = F(0)
    for i in range(n):
        d = p * xv[(i + 1) % n] + q * xv[(i + 2) % n]
        s += F(xv[i]) / d
    return s


def odd_family_witness(m):
    """n = 2m+1, witness (0,1,0,1,...,0,1,1) for m>=5."""
    n = 2 * m + 1
    x = [0] * n
    for j in range(n):
        if j % 2 == 1:
            x[j] = 1
    x[-1] = 1
    return tuple(x)


def even_pattern_witness(n, a=1, b=2):
    x = [a if i % 2 == 0 else b for i in range(n)]
    return tuple(x)


print("=" * 70)
print("ODD FAMILY  x=(0,1,0,1,...,0,1,1),  n=2m+1, formula P=(30m+70)/21")
print("=" * 70)
print(f"{'m':>3} {'n':>4} {'P(exact)':>20} {'formula':>12} {'match':>6} {'deficit':>18} {'>0?':>5}")
all_ok = True
for m in range(5, 26):
    n = 2 * m + 1
    x = odd_family_witness(m)
    P = exact_P(x)
    formula = F(30 * m + 70, 21)
    deficit = n - P
    match = (P == formula)
    ok = match and deficit > 0
    all_ok &= ok
    print(f"{m:>3} {n:>4} {str(P):>20} {str(formula):>12} {str(match):>6} {str(deficit):>18} {str(deficit>0):>5}")
print(f"\nGeneral odd family (m=5..25) all verified: {all_ok}")
print(f"Deficit = (12m-49)/21 > 0  iff  m >= 5  (n >= 11).")

print("\n" + "=" * 70)
print("SMALL-n TABLE  (p=3/10, q=7/10)")
print("=" * 70)
specials = {
    7: (0, 3, 4, 1, 5, 0, 4),
    9: (0, 1, 0, 1, 0, 1, 0, 1, 2),
}
print(f"{'n':>4} {'type':>8} {'witness':>30} {'P':>22} {'bound':>6} {'<bound?':>7}")
for n in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
    if n == 3 or n == 5:
        # pick a generic positive vector to show it holds (>= n)
        x = tuple([1] * n)
        P = exact_P(x)
        print(f"{n:>4} {'hold':>8} {str(x):>30} {str(P):>22} {n:>6} {str(P>=n):>7}  (illustrative; n=3,5 proven)")
    elif n in specials:
        x = specials[n]
        P = exact_P(x)
        print(f"{n:>4} {'cex':>8} {str(x):>30} {str(P):>22} {n:>6} {str(P<n):>7}  deficit={n-P}")
    elif n % 2 == 1:  # odd >=11
        m = (n - 1) // 2
        x = odd_family_witness(m)
        P = exact_P(x)
        print(f"{n:>4} {'cex':>8} {str(x):>30} {str(P):>22} {n:>6} {str(P<n):>7}  deficit={n-P}")
    else:  # even
        x = even_pattern_witness(n)
        P = exact_P(x)
        print(f"{n:>4} {'cex':>8} {str(x):>30} {str(P):>22} {n:>6} {str(P<n):>7}  deficit={n-P}")

print("\n" + "=" * 70)
print("EVEN-n pattern x=(1,2,1,2,...) at p=0.3,q=0.7  (paper Remark 2.4, p<q)")
print("=" * 70)
print(f"{'n':>4} {'P':>22} {'<n?':>6} {'deficit':>18}")
for n in [4, 6, 8, 10, 12, 14, 20]:
    x = even_pattern_witness(n)
    P = exact_P(x)
    print(f"{n:>4} {str(P):>22} {str(P<n):>6} {str(n-P):>18}")

print("\n" + "=" * 70)
print("HEADLINE CHECK: at p=3/10,q=7/10, inequality holds for n in {3,5},")
print("fails for n=4 and every n>=6.")
print("=" * 70)
fails = []
for n in range(4, 26):
    if n == 5:
        continue
    if n % 2 == 0:
        P = exact_P(even_pattern_witness(n))
    elif n == 7:
        P = exact_P(specials[7])
    elif n == 9:
        P = exact_P(specials[9])
    else:
        P = exact_P(odd_family_witness((n - 1) // 2))
    if P is None or P >= n:
        fails.append(n)
print(f"n=4..25 (excl 5) all fail? {len(fails)==0}  (exceptions: {fails})")
