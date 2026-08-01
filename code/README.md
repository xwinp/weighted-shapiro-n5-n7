# Code supplement — role map

Scripts are classified into **three tiers** (authoritative list in `ARCHIVAL.md` → "Script tiering"):

- **Tier 1 — load-bearing core (8):** the computational core whose independent re-verification is the focus of the archival manifest. Must pass.
- **Tier 2 — formal theorem certificates:** exact-rational or rigorous-interval scripts certifying theorem pieces outside the Tier 1 core (n=5, n=7 S₂/S₃/S₀, and the n=7 S₁ crossing/band/sign/monotonicity pieces). Must also pass.
- **Tier 3 — diagnostic / corroboration / exploratory:** numerical cross-checks, corroboration-only scripts, superseded cross-checks, traces, scans, and companion-paper / scratch code. *Not* load-bearing — a timeout or display quirk in Tier 3 affects no theorem.

`code/run_all_certificates.py` runs every Tier 1 + Tier 2 *certificate* script and asserts a zero exit code for each. The captured verification log is `code/_archival_run.txt`; repository-wide hashes are in `sha256_manifest.txt`.

## Tier 1 — load-bearing core (8)

| Script | Certificate |
|---|---|
| `verify_ny_spectral.py` | exact N–Y spectral bridge (Lemma 2.2): μ_k identity, disk equivalence, Δ₅,Δ₇<0 ⇒ all μ_k>0 |
| `verify_hc_closedform.py` | exact identity P = C + 2√(AB) (numerator ≡0; squared check) |
| `verify_s1_elimination.py` | exact forward containment Res_v(F1,F2) = w²z²(w−1)²H_B H_C (no third component) |
| `verify_s1_hc_superset.py` | exact H_C containment: g³·E₂^u(r)=z⁴(w−1)²(z−1)²·B·A·H_C, Res_w(g,H_C)=−z¹¹(z−1)⁶ ⇒ H_C=0 ⇒ E₂^red=0 |
| `n7_s1_hc_s1collar.py` | 60 full-parameter interval-Newton boxes tiling [s_max,1], regular δ=0 tail |
| `n7_s1_hc_exhaustiveness_cert.py` | 2×2-Krawczyk branch graph (2 components), six-cell bijection, certified critical fibers |
| `n7_s1_hc_cover_checker.py` | mpmath re-check: 2604/2604 boxes, Krawczyk/admissibility/P_mv>141/20 |
| `n7_s1_hc_arb_checker.py` | independent Arb/python-flint check: 2604/2604, min 7.050005>141/20 |

## Tier 2 — formal theorem certificates (must also pass)

| Theorem piece | Script | Certificate |
|---|---|---|
| n=5, S₁ | `n5_s1_certificate.py` | Res_a(curve,Mnum)=r¹⁷·Q₁₁, Q₁₁ all-positive coeffs *(closing sanity loop display-only)* |
| n=5, S₁ | `n5_s1_m1_interval.py` | M(1)∈[0.5662,0.5821]>1/2 (rational interval) |
| n=7, S₂ band | `n7_resultant.py` | exact quotient Res_t(R,B)=p¹⁵(p−1)⁶F; F deg-15 squarefree; F mod 23 irreducible ⟹ irreducible over Q; exactly two (0,1)-roots = a₇,b₇ |
| n=7, endpoints | `verify_endpoints.py` | exact Sturm isolation: a₇,b₇ (11 dig), p₁,p₂ (9 dig) within 0.5 ulp; {p₁,p₂}⊂(a₇,b₇) |
| n=7, S₂ signs | `n7_s2_rational_signs.py` | three-sign rational-interval samples; D0,D5>0 at p=1/4 |
| n=7, inactive KKT | `n7_inactive_sturm.py` | G (degree 21), 0 roots in (a₇,b₇) |
| n=7, S₃ | `verify_s3_closed.py` | exact (1−r)H(r) identity chain ⇒ P_S₃^stat > 8 > 7 |
| n=7, S₀ uniform | `verify_uniform_hess.py` | generalised eigenproblem = Fourier formula; discriminant<0 ⇒ strict local min |
| n=7, S₁ crossing | `n7_s1_crossing_resultant.py` | Res_z=896ρ¹³(ρ+1)⁷(8ρ²+8ρ+7)⁶Φ₃₅; Sturm counts |
| n=7, S₁ band | `_verify_crossings_in_band.py` | exact Sturm isolation of a₇,b₇,p₁,p₂ ⇒ {p₁,p₂}⊂(a₇,b₇) |
| n=7, S₁ sign portrait | `_verify_s1_three_samples.py` | three rigorous interval samples (Krawczyk-unique z, native outward √) |
| n=7, S₁ det/p₀ | `n7_s1_rigorous_certs.py` | P(2/5)>7; det signs; w_det(z₇)<0; p₀∈(3/8,2/5) |
| n=7, S₁ monotonicity | `n7_s1_rho_monotone_Lcert.py` | Res_w(H_B,L)=z(z−1)Q₇; ρ(z) strictly increasing |
| n=7, H_C superset | `n7_s1_hc_kkt_check.py` | β-variety is a superset (spurious branch fails g₄,g₅) |

Tier 2 also includes the H_C **build / events** scripts, which are formal but *produce* records rather than *verify* them (so they are not in `run_all_certificates.py` — their output is verified by the Tier 1 checkers): `n7_s1_hc_critical_events.py` (5 critical z-events), `n7_s1_hc_cover_dump.py` + `n7_s1_hc_arb_refine.py` (build canonical 2604-piece `_hc_cover.json`).

## Tier 3 — diagnostic / corroboration / exploratory (NOT load-bearing)

- **Corroboration only:** `_verify_phi35_squarefree.py` (gcd(Φ₃₅,Φ₃₅′)=1; sign portrait certified by `_verify_s1_three_samples.py`, not by squarefree alternation).
- **Superseded numerical cross-check:** `n7_s1_hc_completeness.py` (superseded by `n7_s1_hc_exhaustiveness_cert.py`).
- **Non-load-bearing diagnostic / numerical cross-check:** `verify_hb_exact_backsub.py`, `verify_beta_reduction.py`, `verify_gpt_formulas.py`, `s1_beta_reduction.py` (rigorous replacements: `verify_uniform_hess.py` / `verify_s1_elimination.py`).
- **Traces / scans / exploratory:** `n7_s1_correct_morse_trace.py`, `n7_s1_all_branches_correct_hess.py`, `n7_s1_hc_definitive_trace.py`, `n7_s1_hc_Cstruct.py`, `n7_s1_hc_mv_verify.py`, `n7_s1_hc_Pformula_check.py`, `n7_s1_rho_monotone_check.py`, `n7_roots.py` (numerical root display — slow; exact a₇,b₇ isolation lives in `_verify_crossings_in_band.py`), `n7_s1_crossing_scan.py`, `n7_s1_fullrange.py`, `n7_s1_zrange.py`, `n7_s1_branches_outofband.py`, `n7_s1_vs_s2_outofband.py`, `n7_s1_hc_branch_reconstruct.py`, `n7_s1_hc_crossing_probe.py`, `n7_s1_hc_reconstruct2.py`, `n7_s1_hc_struct_explore.py`, and all `n9_*` / `routeA_*` / `routeB_*` / `classify_*` / `phase*` / `s0_*` / `s1_*` / `_scratch_*` / `shapiro_*` / `a295196_*` / `agrawal_*` scripts not listed in Tier 1 or 2.

## Reproducibility

Run the Tier 1 + Tier 2 certificate chain from the repository root:

```bash
python code/run_all_certificates.py      # runs every Tier 1 + Tier 2 certificate, asserts exit 0
```

or individually (Tier 1 core):

```bash
python code/verify_ny_spectral.py
python code/verify_hc_closedform.py
python code/verify_s1_elimination.py
python code/verify_s1_hc_superset.py
python code/n7_s1_hc_s1collar.py
python code/n7_s1_hc_exhaustiveness_cert.py
python code/n7_s1_hc_cover_checker.py
python code/n7_s1_hc_arb_checker.py
```

Exact package versions are in `code/requirements_lock.txt`; captured output is in `code/_archival_run.txt`.
