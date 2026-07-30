#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GPT's Δ = sigma^3 * J verification script (item 3 of reply #4), run locally.
HERE pointed at paper/_gpt_artifacts. Definitive resolution of the deg-52-vs-56
discrepancy: if assert passes, GPT is right and my earlier det had a parsing bug."""
from pathlib import Path
import re
import sympy as sp

HERE = Path(__file__).resolve().parent.parent / 'paper' / '_gpt_artifacts'

X, Y, sigma = sp.symbols("X Y sigma")


def load_small_clean_poly(filename):
    s = sp.symbols("s")
    text = (HERE / filename).read_text(encoding="utf-8").strip()
    expr = sp.sympify(text, locals={"X": X, "Y": Y, "s": s}).subs(s, sigma)
    return sp.Poly(expr, X, Y, sigma, domain=sp.ZZ)


def load_large_J_sparse(filename):
    text = (HERE / filename).read_text(encoding="utf-8").replace(" ", "").replace("\n", "")
    if not text.startswith(("+", "-")):
        text = "+" + text
    terms = {}
    for sign, term in re.findall(r"([+-])([^+-]+)", text):
        coeff = -1 if sign == "-" else 1
        exponents = []
        remainder = term
        for name in ("X", "Y", "s"):
            match = re.search(rf"(?<![A-Za-z]){name}(?:\*\*(\d+))?", remainder)
            if match:
                exponent = int(match.group(1) or 1)
                exponents.append(exponent)
                remainder = remainder[:match.start()] + remainder[match.end():]
            else:
                exponents.append(0)
        remainder = remainder.replace("*", "")
        if remainder:
            coeff *= int(remainder)
        monomial = tuple(exponents)
        terms[monomial] = terms.get(monomial, 0) + coeff
    return sp.Poly.from_dict(terms, (X, Y, sigma), domain=sp.ZZ)


G = load_small_clean_poly("nonpal_G_clean.txt")
S = load_small_clean_poly("nonpal_S_clean.txt")
N = load_small_clean_poly("nonpal_rho9_num.txt")
Den = load_small_clean_poly("nonpal_rho9_den.txt")
J = load_large_J_sparse("nonpal_J_clean.txt")

print("G: tot/terms/degX/degY/degS =", G.total_degree(), len(G.terms()), G.degree(X), G.degree(Y), G.degree(sigma))
print("S: tot/terms/degX/degY/degS =", S.total_degree(), len(S.terms()), S.degree(X), S.degree(Y), S.degree(sigma))
print("N: tot/terms/degS =", N.total_degree(), len(N.terms()), N.degree(sigma))
print("D: tot/terms/degS =", Den.total_degree(), len(Den.terms()), Den.degree(sigma))
print("J: tot/terms/degS =", J.total_degree(), len(J.terms()), J.degree(sigma))

variables = (X, Y, sigma)
grad_G = [G.diff(v) for v in variables]
grad_S = [S.diff(v) for v in variables]
grad_K_num = [Den * N.diff(v) - N * Den.diff(v) for v in variables]

Delta = (
    grad_G[0] * (grad_S[1] * grad_K_num[2] - grad_S[2] * grad_K_num[1])
    - grad_G[1] * (grad_S[0] * grad_K_num[2] - grad_S[2] * grad_K_num[0])
    + grad_G[2] * (grad_S[0] * grad_K_num[1] - grad_S[1] * grad_K_num[0])
)
DeltaP = sp.Poly(sp.expand(Delta.as_expr()), X, Y, sigma, domain=sp.ZZ)
sigma3 = sp.Poly(sigma**3, X, Y, sigma, domain=sp.ZZ)
expected = sigma3 * J

print("Delta: tot/terms/degS =", DeltaP.total_degree(), len(DeltaP.terms()), DeltaP.degree(sigma))
print("expected sigma^3*J: tot/terms/degS =", expected.total_degree(), len(expected.terms()), expected.degree(sigma))
print("Delta == sigma^3*J ?", DeltaP == expected)
