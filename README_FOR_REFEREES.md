# For referees — quick verification guide

This supplement supports the manuscript
**"The Tuan–Thuong weighted Shapiro cyclic inequality for n = 5 and n = 7"**.
It lists the shortest path to re-verify the computer-assisted components and
where to focus review. The complete environment, per-file listings, the
`Φ_35` coefficient vector, and Sturm counts are in `SUPPLEMENT.md`; the
authoritative script role map is in `code/README.md`.

## Version

Submission tag **`v1.0-submission`** (exact commit hash in `COMMIT.txt`; also
`git rev-parse v1.0-submission`). Every file in this archive is listed in
`sha256_manifest.txt` (run `sha256sum -c sha256_manifest.txt` to confirm
integrity).

## Reproduce in one command

```bash
pip install -r code/requirements_lock.txt     # Python 3.13, SymPy 1.14, mpmath 1.3, python-flint 0.9
python code/run_all_certificates.py           # 22 certificates, asserts exit 0
```

Expected final line: `ALL PASS (every certificate exited 0).`
Captured output: `code/_archival_run.txt` (22 × `[exit 0]`).

The full chain runs in a few minutes on a laptop; the slowest script
(`n7_s1_hc_arb_checker.py`, Arb/python-flint) takes ~1–2 min.

## What is certified (22 scripts, two tiers)

**Tier 1 — load-bearing core (8).** The `H_C` cover and the exact elimination:

| Script | Certificate |
|---|---|
| `verify_ny_spectral.py` | Nowosad–Yamagami spectral bridge (Lemma 2.2): `μ_k` identity, disk equivalence, `Δ_5, Δ_7 < 0` ⇒ all `μ_k > 0` |
| `verify_hc_closedform.py` | exact identity `P = C + 2√(AB)` (numerator ≡ 0) |
| `verify_s1_elimination.py` | exact forward containment `Res_v(F_1,F_2) = w²z²(w−1)² H_B H_C` (no third component) |
| `verify_s1_hc_superset.py` | **exact `H_C` containment** (load-bearing): `g(w,z)³·E_2^u(r) = z⁴(w−1)²(z−1)²·B·A·H_C`, `Res_w(g,H_C)=−z¹¹(z−1)⁶` ⇒ `g≠0` ⇒ `H_C=0 ⇒ E_2^red=0` |
| `n7_s1_hc_s1collar.py` | 60 full-parameter interval-Newton boxes tiling `[s_max,1]`; `g_min = 7.344054 > L_C` |
| `n7_s1_hc_exhaustiveness_cert.py` | 2×2-Krawczyk branch graph: `components=2, cells=6, fibers=5, edges=2823` |
| `n7_s1_hc_cover_checker.py` | mpmath re-check: `2604/2604` boxes, Krawczyk/admissibility/`P_MV > 141/20` |
| `n7_s1_hc_arb_checker.py` | **independent** Arb/python-flint check: `2604/2604`, `min 7.050005 > 141/20` |

**Tier 2 — formal theorem certificates (14).** n=5 `S_1`; n=7 `S_2` band
(resultant `F` deg-15 irreducible over Q via mod-23, exactly two `(0,1)`-roots
`a_7, b_7`), `S_3` (`P_S_3^stat > 8 > 7` via the exact `(1−r)H(r)` chain), `S_0`
uniform Hessian, and the `S_1` crossing/band/sign/monotonicity pieces. See
`code/README.md` for the per-script certificate table.

## Where to focus review

1. **`H_C` containment (load-bearing exact).** `verify_s1_hc_superset.py`
   asserts the **whole rational identity**
   `together(g³·E_2^u(r) − z⁴(w−1)²(z−1)²·B·A·H_C) == 0` (not just a numerator)
   plus the defensive `den == h³` / `num == target` checks; the nonvanishing
   factor is `g` (`Res_w(g,H_C) = −z¹¹(z−1)⁶`, nonzero on the admissible
   `z ∈ (0,1)` branch). This is the `S_1` forward-containment certificate.
2. **`H_C` cover exhaustiveness.** `n7_s1_hc_exhaustiveness_cert.py` builds the
   branch graph with 2×2 Krawczyk at rational overlap parameters (proving same
   root across overlaps) and union-finds `components=2`. Each box is
   Krawczyk-unique; the cover is checked by two independent backends.
3. **`S_3` uniqueness.** `verify_s3_closed.py` §5: the `gc, gd` gradients factor
   explicitly to `eq1: dr⁵+er¹⁰=b³`, `eq2: r⁵(d+er⁵)²=b²e` (linear in `d, e`),
   pinning the unique positive KKT point; `ge=0` then factors to
   `−(b−r²)(r+1)(quartic)(r⁴−r³+r²−r+1)=0`, all factors positive for `b,r>0`
   except `b−r²`.

## Machine records

`code/_hc_cover.json` (the 2604-piece cover), `code/_hc_s1collar.json` (the
60-box collar), `code/_hc_exhaustiveness.json` (the branch graph). These are
committed canonical LF; the Tier 1 checkers verify them read-only by default.

## Notes

- No load-bearing step uses floating point. Exact certificates use SymPy
  rational algebra (`sp.expand(...) == 0`, `sp.div` remainder zero,
  `sp.resultant`); interval certificates use outward-rounded `mpmath.iv` / Arb.
- Record-*generators* (`n7_s1_hc_critical_events.py`, `n7_s1_hc_cover_dump.py`,
  `n7_s1_hc_arb_refine.py`) build the committed records; the Tier 1 checkers
  verify their output, so the generators are not in `run_all_certificates.py`.
- Exploratory / companion-paper (`n=9`) / scratch scripts are excluded from
  this archive (they are not load-bearing for the n=5/n=7 theorems); the full
  repository history is available at the submission tag.
