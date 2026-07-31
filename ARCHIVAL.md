# Archival manifest - Tuan-Thuong weighted Shapiro cyclic inequality, n=5 and n=7

**Author:** Weipeng Xue (薛炜鹏), Sun Yat-sen University, xuewp5@mail2.sysu.edu.cn  
**Paper source:** `paper/n57_paper.md`  
**Typeset paper:** `paper/n57_paper.pdf` (21 pages)

This repository contains the source paper and every load-bearing computer-assisted certificate.

## Load-bearing verification chain

Run from the repository root after installing `code/requirements_lock.txt`:

```bash
python code/verify_ny_spectral.py
python code/verify_hc_closedform.py
python code/verify_s1_elimination.py
python code/n7_s1_hc_s1collar.py
python code/n7_s1_hc_exhaustiveness_cert.py
python code/n7_s1_hc_cover_checker.py
python code/n7_s1_hc_arb_checker.py
```

`verify_hb_exact_backsub.py` is retained only as a non-load-bearing diagnostic.

### Division of responsibility

- `verify_ny_spectral.py`: exact Nowosad-Yamagami spectral identities.
- `verify_hc_closedform.py`: exact identity `P=C+2 sqrt(AB)`.
- `verify_s1_elimination.py`: exact forward containment
  `Res_v(F1,F2)=w^2 z^2 (w-1)^2 H_B H_C`.
- `n7_s1_hc_cover_checker.py` and `n7_s1_hc_arb_checker.py`: independent validity checks for all 2604 original cover boxes.
- `n7_s1_hc_s1collar.py`: 60 full-parameter interval-Newton boxes covering `[s_max,1]`, including the regular closed tail at `delta=0`; exact seam and abutment; rechecked lower bound above `141/20`.
- `n7_s1_hc_exhaustiveness_cert.py`: global exhaustiveness. It builds a branch graph whose edges are certified by 2x2 Krawczyk tests, proves two gap-free root components, enumerates all admissible roots at exact rational samples in the six event-free cells, and certifies all five critical fibers using exact `CRootOf` event isolation, full-parameter `c`-root tubes, outward quadratic `v`-enclosures, 3x3 Krawczyk certificates for generic boundary roots, and an exact linear-subresultant certificate for the terminating `u=1` branch. It uses no `nroots`, `polyroots`, `nsimplify`, residual threshold, or floating containment in a load-bearing step.

Machine records:

- `code/_hc_cover.json` - 2604 original boxes.
- `code/_hc_s1collar.json` - 60 second-blow-up collar boxes.
- `code/_hc_exhaustiveness.json` - branch graph, six cell bijections, and five critical-fiber records.
- `code/_archival_run.txt` - captured verification output.

## Script tiering (3 tiers)

To resolve the "7 vs 38" appearance (the paper highlights 7 load-bearing scripts while `code/` contains many more), every script is classified into one of three tiers. **Tier 1** is the load-bearing core above; **Tier 2** are the remaining formal theorem certificates (also exact-rational or rigorous-interval, also must-pass, but certifying theorem pieces outside the H_C cover core); **Tier 3** are diagnostics / corroboration / exploration (never load-bearing — a timeout or display quirk in Tier 3 affects no theorem). `code/run_all_certificates.py` runs every Tier 1 + Tier 2 *certificate* script and asserts a zero exit code for each.

### Tier 1 — load-bearing core (7; the chain above)

`verify_ny_spectral.py`, `verify_hc_closedform.py`, `verify_s1_elimination.py`, `n7_s1_hc_s1collar.py`, `n7_s1_hc_exhaustiveness_cert.py`, `n7_s1_hc_cover_checker.py`, `n7_s1_hc_arb_checker.py`.

### Tier 2 — formal theorem certificates (must also pass)

- **n=5 S₁:** `n5_s1_certificate.py` (Res_a = r¹⁷Q₁₁, Q₁₁ all-positive), `n5_s1_m1_interval.py` (M(1)∈[0.5662,0.5821]>½, rational interval).
- **n=7 S₂ band:** `n7_resultant.py` (exact quotient Res_t(R,B)=p¹⁵(p−1)⁶F via exact polynomial division; F deg-15 squarefree; F mod 23 irreducible ⟹ irreducible over Q; exactly two (0,1)-roots = a₇,b₇), `verify_endpoints.py` (a₇,b₇ to 11 digits, p₁,p₂ to 9 digits — each within 0.5 ulp via exact Sturm; {p₁,p₂}⊂(a₇,b₇)), `n7_s2_rational_signs.py` (D₀,D₅>0 at p=1/4), `n7_inactive_sturm.py` (G deg-21, 0 roots in (a₇,b₇)).
- **n=7 S₃:** `verify_s3_closed.py` (exact (1−r)H(r) identity chain ⇒ P_S₃^stat > 8 > 7).
- **n=7 S₀ uniform:** `verify_uniform_hess.py` (generalised eigenproblem B^T H B v = λ(B^T B)v matches Fourier formula; worst-mode discriminant <0 ⇒ uniform strict local min).
- **n=7 S₁ crossing / band / sign / monotonicity:** `n7_s1_crossing_resultant.py` (Res_z = 896ρ¹³(ρ+1)⁷(8ρ²+8ρ+7)⁶Φ₃₅, Sturm), `_verify_crossings_in_band.py` ({p₁,p₂}⊂(a₇,b₇)), `_verify_s1_three_samples.py` (three rigorous interval samples), `n7_s1_rho_monotone_Lcert.py` (ρ(z) strictly increasing), `n7_s1_rigorous_certs.py` (P(2/5)>7; det signs; p₀∈(3/8,2/5)).
- **n=7 H_C build / events / superset (rigorous, but *produce* records rather than *verify* them — run by the Tier 1 checkers, not by `run_all_certificates.py`):** `n7_s1_hc_critical_events.py` (5 critical z-events), `n7_s1_hc_cover_dump.py` + `n7_s1_hc_arb_refine.py` (build the canonical 2604-piece `_hc_cover.json`), `n7_s1_hc_kkt_check.py` (β-variety is a superset). These are formal but are *generators*; the load-bearing verification of their output is done by `n7_s1_hc_cover_checker.py` / `n7_s1_hc_arb_checker.py` (Tier 1).

### Tier 3 — diagnostic / corroboration / exploratory (NOT load-bearing)

- **Corroboration only:** `_verify_phi35_squarefree.py` (gcd(Φ₃₅,Φ₃₅′)=1; the sign portrait is certified by `_verify_s1_three_samples.py`, not by squarefree alternation).
- **Superseded numerical cross-check:** `n7_s1_hc_completeness.py` (retained only as a non-load-bearing midpoint/count cross-check; superseded by `n7_s1_hc_exhaustiveness_cert.py`).
- **Non-load-bearing diagnostic / forward-construction record:** `verify_hb_exact_backsub.py`, `verify_beta_reduction.py`, `verify_gpt_formulas.py`, `s1_beta_reduction.py` (numerical cross-checks of structural formulas whose rigorous replacements are `verify_uniform_hess.py` / `verify_s1_elimination.py`).
- **Traces / scans / exploratory (companion n=9 paper or scratch):** `n7_s1_correct_morse_trace.py`, `n7_s1_all_branches_correct_hess.py`, `n7_s1_hc_definitive_trace.py`, `n7_s1_hc_Cstruct.py`, `n7_s1_hc_Pformula_check.py`, `n7_s1_hc_mv_verify.py`, `n7_s1_rho_monotone_check.py`, `n7_roots.py`, `n7_s1_crossing_scan.py`, `n7_s1_fullrange.py`, `n7_s1_zrange.py`, `n7_s1_branches_outofband.py`, `n7_s1_vs_s2_outofband.py`, `n7_s1_hc_branch_reconstruct.py`, `n7_s1_hc_crossing_probe.py`, `n7_s1_hc_reconstruct2.py`, `n7_s1_hc_struct_explore.py`, and all `n9_*` / `routeA_*` / `routeB_*` / `classify_*` / `phase*` / `s0_*` / `s1_*` / `_scratch_*` / `shapiro_*` / `a295196_*` / `agrawal_*` scripts not listed in Tier 1 or 2.

## Integrity

`sha256_manifest.txt` uses standard `HASH  path` format. On a fresh extraction:

```bash
sha256sum -c sha256_manifest.txt
```

The ZIP includes `.git`; `git rev-parse HEAD` identifies the snapshot. The SHA-256 manifest remains authoritative even if `.git` is stripped.

**Archived commit (rev 6, prior snapshot):** `36b2691eb7e4e0870800e65a02698926e2401da8`. **Rev 7** supersedes rev 6; **rev 8 (current `HEAD`)** supersedes rev 7. **No mathematics changed**: the certificate logic of all seven load-bearing scripts is identical to rev 4/5/6 (the 3 repaired exact-symbolic diagnostic scripts and the rev-5 editorial corrections are retained). Rev 7 is a packaging/reproducibility cleanup addressing three artifact-integrity findings:
1. **PDF formula integrity.** Several wide display equations previously overflowed the text margin (truncated) or were shrunk by `\resizebox` to an unreadable single line. They are now hand-broken with amsmath `aligned` (multi-line): the four flagged polynomials $Q_{11}(r)$ (p.7), $F(p)$ degree-15 (p.9), the $Q_7$ determinant-resultant line (p.12), $G(p)$ degree-21, plus the $H(r)$ identity and the $\Phi_{35}$ coefficient vector $(c_0,\ldots,c_{35})$ (the last retains its deliberate `\scriptsize`). `code/md2tex.py` gains a `has_multiline` exemption: a display equation using `aligned`/`split`/`multline`/`align` is rendered verbatim — no `\resizebox` (which would shrink it to an unreadable single line and defeat the explicit line breaks) — and the author's `\small`/`\scriptsize` is preserved as a deliberate size choice for genuinely wide content; only the `\resizebox` fallback (wide single-line equations) strips size commands. The $H_n$ definition spacing was tightened (`\qquad`→`\quad`), the author block is set on two lines (name / affiliation+email), and `\emergencystretch` is raised to `4em`. Rebuilt PDF: 21 pages, **0 Overfull hboxes** (xelatex log).
2. **Cross-platform machine records.** The two record-generating load-bearing scripts (`n7_s1_hc_s1collar.py`, `n7_s1_hc_exhaustiveness_cert.py`) were the only scripts that wrote indented JSON, and on Windows they wrote CRLF — so a Linux re-run rewrote the files as LF, dirtying the tree and failing `sha256sum -c`. They are now **read-only by default** (a clean-checkout re-run computes and reports the certificate invariants but does not overwrite the committed record, which is independently verified by `n7_s1_hc_cover_checker.py` / `n7_s1_hc_arb_checker.py`); `--write-record` opts in and writes **canonical LF** (`newline='\n'`). The committed `_hc_s1collar.json` / `_hc_exhaustiveness.json` are converted to LF and are **byte-identically reproduced** by `--write-record` (verified by SHA-256 match on the author machine); a read-only re-run leaves them unchanged.
3. **Repository compile-cleanliness.** `code/n9_core_continuation.py:113` had a missing `]` (`compileall` failed); fixed. `python -m compileall code/` is now clean. (n=9 code is exploratory, not load-bearing for the n=5/n=7 theorems.)

The certificate invariants on the author machine (python-flint 0.9.0) are identical to rev 4/5/6: collar `gmin=7.344054>L_C`, exhaustiveness `components=2 cells=6 fibers=5 edges=2823`, Arb `2604/2604` pieces verified. The current snapshot is identified by `git rev-parse HEAD`; the SHA-256 manifest is authoritative even if `.git` is stripped. Earlier snapshots (rev 5 `41227cc`, rev 4 `495b41c`) remain in git history.

**Archived commit (rev 7, prior snapshot):** `ec7643280cb28af911290dae3d4f5858d3f9a7c0`. **Rev 8 (current `HEAD`)** supersedes rev 7. **No mathematics changed**: all theorem logic and certificate invariants are identical to rev 4–7 (collar `gmin=7.344054>L_C`, exhaustiveness `components=2 cells=6 fibers=5 edges=2823`, Arb `2604/2604`, `min 7.050005>141/20`). Rev 8 is a script-hardening + reproducibility cleanup responding to GPT round-14 artifact findings:
1. **Formal-assert sweep.** Every Tier 1 + Tier 2 *certificate* script now terminates with `assert <invariants>` + `sys.exit(0 if ok else 1)`, so pass/fail is machine-decisive; `code/run_all_certificates.py` asserts exit 0 for each (now 7 Tier 1 + 14 Tier 2 = 21 certificates, including the new `verify_endpoints.py`).
2. **`n7_resultant.py` exact quotient.** `Res_t(R,B)=p^15(p-1)^6 F` is certified by *exact polynomial division* (quotient + zero remainder), not numerical factor evaluation; `F` squarefree (`gcd(F,F')=1`); `F mod 23` irreducible (finite-field cert ⟹ irreducible over Q, Gauss); Sturm count 2 in `(0,1)`. The paper's mod-23 claim is thus re-certified alongside the quotient.
3. **`verify_endpoints.py` (new, Tier 2).** Exact Sturm isolation confirms the four cited decimals to their stated precision (`a_7,b_7` to 11 digits, `p_1,p_2` to 9 digits, each within `0.5` ulp) and `{p_1,p_2}⊂(a_7,b_7)`.
4. **Cross-platform stdout.** `run_all_certificates.py` sets `PYTHONUTF8=1` for child processes so non-ASCII verdict lines survive on CJK (cp936/gbk) Windows consoles; one non-ASCII superset glyph in `n7_inactive_sturm.py` is replaced with ASCII.
5. **Pickle dependency removed** (rev-7 working-tree change, now committed). The three H_C generators (`n7_s1_hc_cover_dump.py`, `n7_s1_hc_rigorous_cert.py`, `n7_s1_hc_completeness.py`) no longer `pickle.load` the uncommitted `_hc_critical_z.pickle`; they import `compute_critical_z` from `n7_s1_hc_critical_events.py` and recompute (5 critical z-values). The committed `_hc_cover.json` / `_hc_s1collar.json` / `_hc_exhaustiveness.json` are unchanged.
6. **PDF rebuilt.** `paper/n57_paper.pdf` regenerated via `md2tex.py` + XeLaTeX (xeCJK/SimSun canonical): 21 pages, 0 Overfull hboxes (xelatex log). `paper/n57_paper.md` §4.2, `code/README.md`, and the Tier-2 list above are updated for (2) and (3).

The captured certificate log is `code/_archival_run.txt` (rev 8: 21 certificates, full per-script output with `[exit 0]` markers). The current snapshot is identified by `git rev-parse HEAD`; the SHA-256 manifest is authoritative even if `.git` is stripped.

## Build the paper

```bash
python code/md2tex.py
cd paper
xelatex -interaction=nonstopmode -halt-on-error n57_paper.tex
xelatex -interaction=nonstopmode -halt-on-error n57_paper.tex
```

The generated `.tex` is regenerable and is not load-bearing.

**Canonical-PDF policy.** The committed `paper/n57_paper.pdf` (built with XeLaTeX + `xeCJK`, SimSun on the author's Windows machine) is the **canonical artifact**: its SHA-256 is recorded in `sha256_manifest.txt` and it is what the SHA-256 manifest and fresh-extract verification check against. The build uses a `\IfFontExistsTF` CJK-font fallback (Noto Serif CJK SC → SimSun), so a source rebuild on a machine with a *different* CJK font available (e.g. Linux with Noto but not SimSun) selects a different CJK font and produces **different pagination and TOC page numbers** (though identical mathematical content, 21 pages, and 0 Overfull hboxes). This is a declared property of the cross-platform source build, not a defect: the canonical PDF is byte-stable and is the artifact of record; a cross-platform source rebuild is guaranteed to reproduce the **content** and the **no-overflow** property, but is **not** guaranteed to reproduce pagination or byte-identity. To obtain the exact canonical PDF, build on the author's locked environment (Windows + MiKTeX + SimSun) or use the committed PDF directly.

## Dependencies

Pinned in `code/requirements_lock.txt`: Python 3.13.9, SymPy 1.14.0, mpmath 1.3.0, NumPy 2.4.4, SciPy 1.18.0, and python-flint 0.9.0.

## Deposition

The Zenodo DOI will be recorded here upon deposition. The DOI is a reproducibility convenience, not a proof ingredient.
