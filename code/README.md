# Code supplement — role map

Scripts are split into **formal certificates** (load-bearing for the theorems;
each is exact-rational or rigorous-interval and must pass) and **exploration /
numerical** scripts (used to discover the structure; *not* load-bearing — a
timeout or display quirk in one of these does not affect any theorem).

The paper cites only **formal certificate** scripts. Frozen re-run logs and
SHA-256 hashes live in `code/_frozen_logs/`.

## Formal certificates (load-bearing)

| Theorem piece | Script | Certificate |
|---|---|---|
| n=5, S₁ | `n5_s1_m1_interval.py` | M(1)∈[0.5662,0.5800]>1/2 (rational interval) |
| n=5, S₁ | `n5_s1_certificate.py` | Res_a(curve,Mnum)=r¹⁷·Q₁₁, Q₁₁ all-positive coeffs *(the closing sanity loop is display-only and not part of the certificate)* |
| n=7, S₂ band | `n7_resultant.py` | F (degree 15), degree-15 factor contains a₇,b₇ |
| n=7, S₂ signs | `n7_s2_rational_signs.py` | three-sign rational-interval samples; D0,D5>0 at p=1/4 |
| n=7, inactive KKT | `n7_inactive_sturm.py` | G (degree 21), 0 roots in (a₇,b₇) |
| n=7, S₁ crossing | `n7_s1_crossing_resultant.py` | Res_z=896ρ¹³(ρ+1)⁷(8ρ²+8ρ+7)⁶Φ₃₅; Sturm counts |
| n=7, S₁ band | `_verify_crossings_in_band.py` | exact Sturm isolation of a₇,b₇,p₁,p₂ ⇒ {p₁,p₂}⊂(a₇,b₇) |
| n=7, S₁ sign portrait | `_verify_s1_three_samples.py` | three rigorous interval samples (Krawczyk-unique z, native outward √) |
| n=7, S₁ corroboration | `_verify_phi35_squarefree.py` | gcd(Φ₃₅,Φ₃₅′)=1 (corroboration only) |
| n=7, S₁ det/p₀ | `n7_s1_rigorous_certs.py` | P(2/5)>7; det signs; w_det(z₇)<0; p₀∈(3/8,2/5) |
| n=7, S₁ monotonicity | `n7_s1_rho_monotone_Lcert.py` | Res_w(H_B,L)=z(z−1)Q₇; ρ(z) strictly increasing |
| n=7, S₃ | `verify_s3_closed.py` | P_S₃^stat=(1+r⁵)(5−r⁷)/r²; min 8.009975>7 |
| n=7, S₀ uniform | `verify_uniform_hess.py` | 6×6 Hred eigenvalues = closed-form spectrum |
| n=7, H_C completeness | `n7_s1_hc_completeness.py` | per-cell Sturm==iso, all_cells_ok=True |
| n=7, H_C cover (Arb, **formal backend**) | `n7_s1_hc_arb_checker.py` | 2604/2604, global min 7.050005>141/20 |
| n=7, H_C cover (mpmath, re-check) | `n7_s1_hc_cover_checker.py` | 2604/2604 independent re-verification |
| n=7, H_C cover build | `n7_s1_hc_arb_refine.py` / `n7_s1_hc_cover_dump.py` | canonical 2604-piece cover → `_hc_cover.json` |
| n=7, H_C superset | `n7_s1_hc_kkt_check.py` | β-variety is a superset (spurious branch fails g₄,g₅) |
| n=7, H_C events | `n7_s1_hc_critical_events.py` | five critical z-events partitioning (0,1) |
| β-reduction | `s1_beta_reduction.py` / `verify_beta_reduction.py` | β-coordinate KKT reduction |

## Exploration / numerical (NOT load-bearing)

`n7_roots.py` (numerical root display — slow, may time out; the *exact* a₇,b₇
isolation lives in `_verify_crossings_in_band.py`), `n7_s1_correct_morse_trace.py`,
`n7_s1_all_branches_correct_hess.py`, `n7_s1_hc_definitive_trace.py`,
`n7_s1_hc_Cstruct.py`, `n7_s1_hc_mv_verify.py`, `n7_s1_rho_monotone_check.py`,
`n7_s1_hc_branch_reconstruct.py`, `n7_s1_hc_crossing_probe.py`,
`n7_s1_hc_reconstruct2.py`, `n7_s1_hc_struct_explore.py`, `n7_s1_hc_Pformula_check.py`,
`n7_s1_crossing_scan.py`, `n7_s1_fullrange.py`, `n7_s1_zrange.py`,
`n7_s1_branches_outofband.py`, `n7_s1_vs_s2_outofband.py`, and all `n7_*` /
`classify_*` / `phase*` / `s0_*` / `s1_*` scripts not listed above.

## Reproducibility

```
python code/_frozen_logs/run_n5_s2.*   # see code/_frozen_logs/ for the exact commands
```

Software: Python 3.13.9, SymPy 1.14, mpmath 1.3, python-flint 0.9 (Arb),
NumPy 2.4.4, SciPy 1.18. SHA-256 of every formal script and of `paper/n57_paper.md`
are recorded in `code/_frozen_logs/sha256.txt`.
