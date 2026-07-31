# Archival manifest — Tuan–Thuong weighted Shapiro cyclic inequality, $n=5$ and $n=7$

**Author:** Weipeng Xue (薛炜鹏), Sun Yat-sen University, xuewp5@mail2.sysu.edu.cn
**Paper:** `paper/n57_paper.md` (Theorem 1: $n=5$, $H_5=(0,1)$; Theorem 2: $n=7$, $H_7=(0,a_7]\cup[b_7,1)$)
**Purpose:** This document makes every load-bearing artifact of the proof machine-checkable and independently reproducible. It is the companion to §5 of the paper.

---

## 1. What is archived

| Group | Path | Count | Role |
|---|---|---|---|
| Paper (source of truth) | `paper/n57_paper.md` | 1 | Full proof, $n=5$ and $n=7$ |
| Paper (typeset deliverable) | `paper/n57_paper.pdf` | 1 | 20-page PDF built from the `.md` via `code/md2tex.py` + xelatex |
| Exact symbolic certificates | `code/verify_*.py` | 14 | N–Y spectral bridge, $P=C+2\sqrt{AB}$, $S_1$ elimination completeness, $H_B$ forward record, $\beta$-reduction, uniform Hessian, $S_3$ closed form, GPT-formula cross-checks |
| $n=5$ certificates | `code/n5_*.py` | 2 | stationary curve, $M(1)$ interval |
| $n=7$ certificates | `code/n7_*.py` | 42 | $S_2$ resultant/Sturm/rational-signs, $S_1$ crossing/monotonicity/determinant/three-sample, $H_C$ critical-events/cover/Krawczyk/Arb/exact-completeness, inactive-KKT, etc. |
| Machine-checkable cover | `code/_hc_cover.json` | 1 | 2604-piece $H_C$ validated cover (rational $s/c/v$ boxes, admissibility lower bounds, $L_C=141/20$) |
| $s\to1$ collar (desingularised) | `code/_hc_s1collar.json` | 1 | 60-piece second-blow-up collar tiling $[s_{\max},1]$, $\min P_{\rm lo}^{\rm re}=7.351180>L_C$ |
| Exact exhaustiveness cert (output) | `code/_hc_completeness_exact.json` | 1 | per-cell exact Sturm/`CRootOf`/box-chain $s$-tiling/critical-fiber $(c,v)$ match records (`all_cells_ok=True`, `s_tiling_ok=True`, `all_event_fibers_cv_matched=True`, `collar_integrated=True`) |
| Verifier/checker run output | `code/_archival_run.txt` | 1 | captured stdout of all 8 pass scripts (this run) |
| Dependency lock | `code/requirements_lock.txt` | 1 | exact package versions |
| Code README | `code/README.md` | 1 | per-script index |
| Integrity manifest | `sha256_manifest.txt` | 1 | SHA-256 of every archived file (standard `HASH  path` form; this file excluded) |
| This document | `ARCHIVAL.md` | 1 | this file |

**The eight load-bearing reproducible scripts** (the ones captured in `_archival_run.txt`, §3): six exact SymPy — `verify_ny_spectral.py`, `verify_hc_closedform.py`, `verify_s1_elimination.py`, `verify_hb_exact_backsub.py` (non-load-bearing forward record), `n7_s1_hc_completeness_exact.py` (exact-algebraic exhaustiveness: box-chain $s$-tiling + critical-fiber $(c,v)$ match), `n7_s1_hc_s1collar.py` ($s\to1$ second-blow-up collar) — and two independent interval backends — `n7_s1_hc_cover_checker.py` (mpmath.iv), `n7_s1_hc_arb_checker.py` (Arb/python-flint).

**Excluded (regenerable or out of scope):** `__pycache__`/`*.pyc`; `*.pickle` (regenerable; see the dependency note below); browser-act/ChatGPT driver scripts (`_paste*.js`, `_send*.js`, `_poll*.sh`, …) and session logs/screenshots; exploration `*_out.txt`/`*_fb.txt`/`_gpt_msg*.txt`/`_gpt_verdict.txt` logs (superseded by `_archival_run.txt`; `_gpt_verdict.txt` is deleted from the formal artifact); `_hc_cover_refined.json` (duplicate of `_hc_cover.json`); and the unrelated number-theory-falsification scripts (`a295196_*`, `agrewal_*`, `classify_*`, `classical_shapiro_check.py`). The intermediate `paper/n57_paper.tex` is regenerable from `n57_paper.md` via `code/md2tex.py` and is not load-bearing.

**Dependency note (`_hc_critical_z.pickle`).** The superseded float-based `n7_s1_hc_completeness.py` reads `code/_hc_critical_z.pickle` (regenerable by `n7_s1_hc_critical_events.py`). The **load-bearing** `n7_s1_hc_completeness_exact.py` is self-contained: it recomputes the critical $z$-events itself as exact `CRootOf` algebraic numbers and does **not** read the pickle. The pickle is therefore excluded from the manifest as a regenerable intermediate.

---

## 2. Integrity (SHA-256)

`sha256_manifest.txt` lists every archived file in the **standard `HASH  path` form** (two-space separator, no byte suffix), so `sha256sum -c sha256_manifest.txt` passes verbatim on a fresh extract. It covers 170 files (it does not list itself, by the no-self-reference convention). Key entries:

```
0ccfceb0579d00ce10f55e9388323b1d3986c704f02ef1d83c0fae99f71fa60d  .gitignore
8dcfcf439cae6380a5f08a399a52b53bd37a1a306479ad96bd96ba920326b574  paper/n57_paper.md
a4e6640883a20c069457d62f36c0a63382dcd8654499f1f57fb4344afb9ccf42  paper/n57_paper.pdf
d84ba16037baa621c5da01bbc54c418dda08c9324707795efb26b236089f4541  code/_hc_cover.json
e6abf7a3563728bef34b434e0904aed2d39f146ce4cb70b1693a988cfcca4c68  code/_hc_s1collar.json
bebbeb9a95eaa23a06ec2fae748e1c62da53de47c5d93b850cbc05465c29f46a  code/_hc_completeness_exact.json
23561ab920f4cebd4b35814eacdf8c28a56933820f20398387034bd309deba39  code/n7_s1_hc_completeness_exact.py
3ea75c3d0acf39793f8901927e02e05222ac5473f4acf0554aebe77b84ff827d  code/n7_s1_hc_s1collar.py
ed1e20b2d0a943b38ab21964453b08ffb3ddc6a02f1331a4cc16339f4b3adbe5  code/verify_hb_exact_backsub.py
76a0a5e53a71ef43536542e2a74bdd3f988db2bff482ba16cb68807a1c14cbe0  code/_archival_run.txt
001bfb29e9611bc2c2eac7fd9db318bbc7aad3bdcbdf42f075d75bfdd0de6e80  code/requirements_lock.txt
4ea6bb42524b5686e732309e2c03b6722618da3f5ceeaeeb21b1154baf38f271  code/README.md
```

Verify on any checkout: `sha256sum -c sha256_manifest.txt` → `170 OK, 0 FAILED`.

---

## 3. Reproduce (8 scripts, all pass)

Environment: Python 3.13.9, Windows 11. `pip install -r code/requirements_lock.txt`. Run from the **repository root**:

```bash
python code/verify_ny_spectral.py            # EXACT: Lemma 2.2 N-Y spectral bridge (mu_k identity, disk equiv, Delta_n<0)
python code/verify_hc_closedform.py          # EXACT: P = C + 2sqrt(AB)  (Issue 6)
python code/verify_s1_elimination.py         # EXACT: Res_v(F1,F2) = w^2 z^2 (w-1)^2 H_B H_C  (Issue 7; load-bearing forward containment)
python code/verify_hb_exact_backsub.py       # EXACT (non-load-bearing record): H_B | Res_v(F1,F2) (sp.div rem 0); F_1 does NOT divide L_w (honest)
python code/n7_s1_hc_completeness_exact.py   # EXACT: algebraic EXHAUSTIVENESS cert (CRootOf events, box-chain s-tiling of [0,1], bijective lift<->box, critical-fiber (c,v) match, 5 event fibers)
python code/n7_s1_hc_s1collar.py             # EXACT: s->1 second blow-up (delta=1-s, c=delta*cbar, v=delta*vbar); 60-piece collar tiling [s_max,1]; min P_lo^re = 7.351180 > L_C  (s->1 tail VALIDITY)
python code/n7_s1_hc_cover_checker.py        # INTERVAL (mpmath.iv): 2604/2604 pieces, P_MV > L_C = 141/20  (box VALIDITY)
python code/n7_s1_hc_arb_checker.py          # INTERVAL (independent Arb/python-flint): 2604/2604, global min 7.050005023 > 7  (box VALIDITY)
```

Expected: each prints `DONE-*` and (for the checkers) `ALL PIECES ... VERIFIED ... : True`; the exact-completeness script prints `DONE-COMPLETENESS-EXACT all_cells_ok=True s_tiling_ok=True event_fibers_cv_matched=True collar_integrated=True`; the collar script prints `DONE-S1COLLAR n=60 gmin=7.351180 gmin_re=7.351180 >L_C=True revify=True seam=True abut=True`. The captured output of one full run is in `code/_archival_run.txt`.

The six `EXACT` scripts use only SymPy rational algebra / exact `CRootOf` algebraic numbers (no floating point in any load-bearing step); the two interval scripts use outward-rounded `mpmath.iv` and an independent hand-written Arb interval layer respectively, and agree on 2604/2604 pieces. The two interval checkers certify box **validity** on $(0,s_{\max}]$; `n7_s1_hc_s1collar.py` certifies the $s\to1$ tail $(s_{\max},1]$ by the second blow-up; `n7_s1_hc_completeness_exact.py` certifies the orthogonal **exhaustiveness** property — every admissible lift on $(0,s_{\max}]$ is enumerated by a cover box (bijective lift$\leftrightarrow$box match per event-free cell), the cover$\cup$collar union tiles $[0,1]$ with no gaps (box-chain $s$-tiling, `s_tiling_ok=True`), and at each of the 5 critical $s^\ast$ every system root $(c,v)$ is matched to a straddling box with every admissible root hosted and every unhosted root inadmissible (`all_event_fibers_cv_matched=True`). Together they close both halves of the $H_C$ cover argument on all of $(0,1)$.

---

## 4. Version control

This is a self-contained local git repository (independent of any parent repo). The proof artifacts are committed in a single snapshot whose SHA-256 manifest is §2.

- **Authoritative integrity check:** `sha256sum -c sha256_manifest.txt` (§2) — self-contained, works on **any** extract, including the bare zip (no git required). This is the primary integrity check.
- **Current commit:** `git rev-parse HEAD`. The proof-artifact zip **includes `.git`**, so `git rev-parse HEAD` works on a fresh extract and returns the artifact commit hash. (If `.git` is stripped, fall back to the SHA-256 manifest, which is self-contained.)
- **Zenodo DOI:** `<to be minted on deposition>` (a reproducibility convenience, not a proof ingredient; recorded here upon deposition)

---

## 5. Dependency lock

`code/requirements_lock.txt` pins: `sympy==1.14.0`, `mpmath==1.3.0`, `numpy==2.4.4`, `scipy==1.18.0`, `python-flint==0.9.0` (Python 3.13.9). The EXACT certificates are version-robust across SymPy ≥ 1.12; the interval certificates are pinned for bit-reproducibility.
