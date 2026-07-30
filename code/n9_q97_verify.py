#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q97 (P=9 crossing on {0}-face palindromic branch) verification harness.
GPT claims: Q97(C) degree 97; under Mobius C=(t+2)/(t+1), the transformed
poly (t+2)^97 Q97((t+2)/(t+1)) [cleared to integer coeffs in t] has exactly ONE
sign change in its coefficient sequence -> exactly 1 positive root -> p0 rigorous.

USAGE: fill in Q97_COEFFS (list, degree-97 to degree-0) when GPT provides them.
Then run. Verifies: (1) Mobius transform, (2) sign-change count = 1,
(3) Sturm root count in (0, +inf) and in the t-interval mapping to C in (1/2,1).
Mobius C=(t+2)/(t+1): C=1/2 -> t=0; C->1^- -> t->+inf. So C in (1/2,1) <=> t in (0,+inf).
"""
import sympy as sp, sys
C, t = sp.symbols('C t')

# ---- FILL IN when GPT provides Q97 coefficients (degree 97 .. 0) ----
Q97_COEFFS = []  # placeholder; e.g. [a97, a96, ..., a1, a0]

def build_q97():
    if not Q97_COEFFS:
        print("Q97_COEFFS empty - awaiting GPT coefficient table. Nothing to verify yet.")
        return None
    assert len(Q97_COEFFS)==98, f"expect 98 coeffs (deg 97), got {len(Q97_COEFFS)}"
    return sum(int(a)*C**(97-i) for i,a in enumerate(Q97_COEFFS))

def main():
    Q = build_q97()
    if Q is None: return
    print(f"Q97 degree: {sp.degree(Q, C)}")
    # Mobius: C=(t+2)/(t+1). Transformed = (t+1)^97 * Q((t+2)/(t+1))  [clears denominator]
    # GPT writes (t+2)^97 Q((t+2)/(t+1)) - try both; the one giving integer coeffs is correct.
    for label, T in [("(t+1)^97 Q((t+2)/(t+1))", (t+1)**97 * Q.subs(C,(t+2)/(t+1))),
                     ("(t+2)^97 Q((t+2)/(t+1))", (t+2)**97 * Q.subs(C,(t+2)/(t+1)))]:
        T = sp.expand(T)
        p = sp.Poly(T, t)
        coeffs = p.all_coeffs()
        # check integrality
        allint = all(c.is_integer for c in coeffs)
        # sign changes
        sc = 0
        nz = [c for c in coeffs if c != 0]
        for i in range(len(nz)-1):
            if (nz[i]>0) != (nz[i+1]>0): sc += 1
        print(f"\n{label}: integer_coeffs={allint}, #sign_changes={sc} (GPT claims 1)")
        if allint:
            P = sp.Poly(T, t)
            npos = P.count_roots(0, sp.oo)
            print(f"  Sturm: positive roots (t>0) = {npos}  [C in (1/2,1) <=> t>0]")
            # also full real root count
            print(f"  Sturm: total real roots = {P.count_roots(-sp.oo, sp.oo)}")
            if npos==1:
                rt = [float(x) for x in sp.nroots(T) if abs(sp.im(x))<1e-6 and sp.re(x)>0]
                if rt:
                    tt=rt[0]; Cval=(tt+2)/(tt+1)
                    # C->rho via (B.6) needs A(B,C); skip, just report C root
                    print(f"  unique positive t-root = {tt:.6f} -> C = {float(Cval):.10f}")
                    print(f"  (compare C0 at P=9 crossing = 0.5769552413)")

if __name__ == "__main__":
    main()
