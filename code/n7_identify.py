#!/usr/bin/env python3
"""Identify a_7, b_7 algebraically via PSLQ."""
import mpmath as mp
mp.mp.dps = 60

a7 = mp.mpf('0.214273520909840965587742319090018387977454549')
b7 = mp.mpf('0.328627677916591966697340856477110865008787303')

S = a7 + b7
R = a7 * b7
print(f"a_7 = {mp.nstr(a7, 45)}")
print(f"b_7 = {mp.nstr(b7, 45)}")
print(f"S=a+b = {mp.nstr(S, 45)}")
print(f"R=ab  = {mp.nstr(R, 45)}")
print(f"a,b roots of t^2 - S t + R; disc = S^2-4R = {mp.nstr(S*S-4*R, 20)}")
print(f"sqrt(disc) = {mp.nstr(mp.sqrt(S*S-4*R), 30)}")

print("\n-- PSLQ: is S rational? (relation among 1, S) --")
for den_bound in [1]:
    rel = mp.pslq([mp.mpf(1), S], maxcoeff=10**8, maxsteps=5000)
    print(f"  pslq([1,S]) = {rel}  => S = {-rel[0]/rel[1] if rel else None}")

print("\n-- PSLQ: is R rational? --")
rel = mp.pslq([mp.mpf(1), R], maxcoeff=10**8, maxsteps=5000)
print(f"  pslq([1,R]) = {rel}  => R = {-rel[0]/rel[1] if rel else None}")

print("\n-- PSLQ: are S,R both rational => quadratic --")
# If S=p/q, R=r/s, polynomial q*s*t^2 - ... let's just find integer relation on [1, t, t^2] for each
for name, t in [("a_7", a7), ("b_7", b7)]:
    rel = mp.pslq([mp.mpf(1), t, t**2], maxcoeff=10**7, maxsteps=5000)
    print(f"  {name}: pslq([1,t,t^2]) = {rel}")
    if rel:
        c0, c1, c2 = rel
        print(f"    => {c2} t^2 + {c1} t + {c0} = 0")

print("\n-- higher degree attempt on a_7 (in case not quadratic) --")
for deg in [3, 4, 6, 8]:
    vec = [a7**i for i in range(deg + 1)]
    rel = mp.pslq(vec, maxcoeff=10**6, maxsteps=5000)
    if rel:
        print(f"  deg {deg}: {rel}")

# Check if a_7,b_7 satisfy the SAME quadratic
print("\n-- common quadratic? evaluate candidate from S,R --")
# Build quadratic with integer coeffs from S,R via pslq on [1, S, R]
rel = mp.pslq([mp.mpf(1), S, R], maxcoeff=10**6, maxsteps=5000)
print(f"  pslq([1, S, R]) = {rel}  (relation c0 + c1*S + c2*R = 0)")
