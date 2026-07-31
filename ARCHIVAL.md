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

## Integrity

`sha256_manifest.txt` uses standard `HASH  path` format. On a fresh extraction:

```bash
sha256sum -c sha256_manifest.txt
```

The ZIP includes `.git`; `git rev-parse HEAD` identifies the snapshot. The SHA-256 manifest remains authoritative even if `.git` is stripped.

**Archived commit (rev 5, prior snapshot):** `41227cc2d9a970dc3f742a33b383742128f0c2e6`. **Rev 6 (current `HEAD`)** supersedes rev 5: the seven load-bearing certificate scripts are identical to rev 4/5 (unchanged); three diagnostic scripts (`n5_s1_certificate.py`, `n7_inactive_sturm.py`, `n7_roots.py`) are repaired to be exact-symbolic (no `mpmath.polyroots`, no float-truncated Sturm — `Poly.count_roots` over `QQ` with rational endpoints, `Poly.intervals` isolation, `factor_list(..., modulus=23)` finite-field cert); `md2tex.py` strips the markdown blockquote marker so no literal `>` appears in the typeset paper; and the rev-5 editorial corrections (intro root ordering, Lemma 2.3 `inf∅=+∞`, S₃ identity `r⁵` term) are retained. The current snapshot is identified by `git rev-parse HEAD`; the SHA-256 manifest is authoritative even if `.git` is stripped.

## Build the paper

```bash
python code/md2tex.py
cd paper
xelatex -interaction=nonstopmode -halt-on-error n57_paper.tex
xelatex -interaction=nonstopmode -halt-on-error n57_paper.tex
```

The generated `.tex` is regenerable and is not load-bearing.

## Dependencies

Pinned in `code/requirements_lock.txt`: Python 3.13.9, SymPy 1.14.0, mpmath 1.3.0, NumPy 2.4.4, SciPy 1.18.0, and python-flint 0.9.0.

## Deposition

The Zenodo DOI will be recorded here upon deposition. The DOI is a reproducibility convenience, not a proof ingredient.
