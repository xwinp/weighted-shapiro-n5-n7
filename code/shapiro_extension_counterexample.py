#!/usr/bin/env python3
"""
VERIFIED counterexample to the Tuan-Thuong (2009) weighted Shapiro cyclic inequality.

Source:
  Nguyen Minh Tuan & Le Quy Thuong,
  "On an Extension of Shapiro's Cyclic Inequality",
  Journal of Inequalities and Applications (2009), DOI 10.1155/2009/491576.

Inequality instance (open per the paper's n>=5 odd classification):
  For x_1,...,x_7 >= 0 (cyclic), with p = 3/10, q = 7/10 (so p+q = 1),
  the claim  P = sum_{i=1}^{7} 10*x_i / (3*x_{i+1} + 7*x_{i+2})  >=  n/(p+q) = 7
  is FALSE.

Minimal (smallest-box, lexicographically-first primitive integer) counterexample:
  x = (0, 3, 4, 1, 5, 0, 4)
  P = 1859/266 = 6.988721804511278...  <  7      (deficit = 3/266)

Independently verified by exhaustive search: 0 counterexamples for box bound B=1..4,
7 counterexamples at B=5; (0,3,4,1,5,0,4) is the lexicographically first.
"""
from fractions import Fraction as F
from math import gcd
from itertools import product

N = 7


def P_value(x):
    """Exact rational value of the weighted Shapiro sum (None if undefined)."""
    s = F(0)
    for i in range(N):
        d = 3 * x[(i + 1) % N] + 7 * x[(i + 2) % N]
        if d == 0:
            if x[i] == 0:
                continue
            return None
        s += F(10 * x[i], d)
    return s


def is_primitive(x):
    nz = [v for v in x if v != 0]
    if not nz:
        return False
    g = nz[0]
    for v in nz[1:]:
        g = gcd(g, v)
    return g == 1


def main():
    witness = (0, 3, 4, 1, 5, 0, 4)
    P = P_value(witness)
    print("Witness  x =", witness)
    print("P =", P, "=", float(P))
    print("7 - P =", 7 - P, "  ->  counterexample:", P < 7)
    print()
    print("Per-term breakdown:")
    for i in range(N):
        d = 3 * witness[(i + 1) % N] + 7 * witness[(i + 2) % N]
        print(f"  i={i+1}: 10*{witness[i]}/(3*{witness[(i+1)%N]}+7*{witness[(i+2)%N]}) "
              f"= {F(10*witness[i], d) if d else 0}")
    print()
    print("Exhaustive primitive-integer search (homogeneous -> scale-invariant):")
    for B in range(1, 7):
        cex = []
        for x in product(range(B + 1), repeat=N):
            if not is_primitive(x):
                continue
            P = P_value(x)
            if P is None or P >= 7:
                continue
            cex.append(x)
        print(f"  B={B}: counterexamples = {len(cex)}")
        if cex:
            print(f"        lexicographically first = {cex[0]}, P = {P_value(cex[0])}")
            break


if __name__ == "__main__":
    main()
