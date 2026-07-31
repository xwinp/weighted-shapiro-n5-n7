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

**Archived commit (rev 6, prior snapshot):** `36b2691eb7e4e0870800e65a02698926e2401da8`. **Rev 7 (current `HEAD`)** supersedes rev 6. **No mathematics changed**: the certificate logic of all seven load-bearing scripts is identical to rev 4/5/6 (the 3 repaired exact-symbolic diagnostic scripts and the rev-5 editorial corrections are retained). Rev 7 is a packaging/reproducibility cleanup addressing three artifact-integrity findings:
1. **PDF formula integrity.** Several wide display equations previously overflowed the text margin (truncated) or were shrunk by `\resizebox` to an unreadable single line. They are now hand-broken with amsmath `aligned` (multi-line): the four flagged polynomials $Q_{11}(r)$ (p.7), $F(p)$ degree-15 (p.9), the $Q_7$ determinant-resultant line (p.12), $G(p)$ degree-21, plus the $H(r)$ identity and the $\Phi_{35}$ coefficient vector $(c_0,\ldots,c_{35})$ (the last retains its deliberate `\scriptsize`). `code/md2tex.py` gains a `has_multiline` exemption: a display equation using `aligned`/`split`/`multline`/`align` is rendered verbatim — no `\resizebox` (which would shrink it to an unreadable single line and defeat the explicit line breaks) — and the author's `\small`/`\scriptsize` is preserved as a deliberate size choice for genuinely wide content; only the `\resizebox` fallback (wide single-line equations) strips size commands. The $H_n$ definition spacing was tightened (`\qquad`→`\quad`), the author block is set on two lines (name / affiliation+email), and `\emergencystretch` is raised to `4em`. Rebuilt PDF: 21 pages, **0 Overfull hboxes** (xelatex log).
2. **Cross-platform machine records.** The two record-generating load-bearing scripts (`n7_s1_hc_s1collar.py`, `n7_s1_hc_exhaustiveness_cert.py`) were the only scripts that wrote indented JSON, and on Windows they wrote CRLF — so a Linux re-run rewrote the files as LF, dirtying the tree and failing `sha256sum -c`. They are now **read-only by default** (a clean-checkout re-run computes and reports the certificate invariants but does not overwrite the committed record, which is independently verified by `n7_s1_hc_cover_checker.py` / `n7_s1_hc_arb_checker.py`); `--write-record` opts in and writes **canonical LF** (`newline='\n'`). The committed `_hc_s1collar.json` / `_hc_exhaustiveness.json` are converted to LF and are **byte-identically reproduced** by `--write-record` (verified by SHA-256 match on the author machine); a read-only re-run leaves them unchanged.
3. **Repository compile-cleanliness.** `code/n9_core_continuation.py:113` had a missing `]` (`compileall` failed); fixed. `python -m compileall code/` is now clean. (n=9 code is exploratory, not load-bearing for the n=5/n=7 theorems.)

The certificate invariants on the author machine (python-flint 0.9.0) are identical to rev 4/5/6: collar `gmin=7.344054>L_C`, exhaustiveness `components=2 cells=6 fibers=5 edges=2823`, Arb `2604/2604` pieces verified. The current snapshot is identified by `git rev-parse HEAD`; the SHA-256 manifest is authoritative even if `.git` is stripped. Earlier snapshots (rev 5 `41227cc`, rev 4 `495b41c`) remain in git history.

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
