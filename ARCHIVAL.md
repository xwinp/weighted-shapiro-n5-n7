# Archival manifest — Tuan–Thuong weighted Shapiro cyclic inequality, $n=5$ and $n=7$

**Author:** Weipeng Xue (薛炜鹏), Sun Yat-sen University, xuewp5@mail2.sysu.edu.cn
**Paper:** `paper/n57_paper.md` (Theorem 1: $n=5$, $H_5=(0,1)$; Theorem 2: $n=7$, $H_7=(0,a_7]\cup[b_7,1)$)
**Purpose:** This document makes every load-bearing artifact of the proof machine-checkable and independently reproducible. It is the companion to §5 of the paper.

---

## 1. What is archived

| Group | Path | Count | Role |
|---|---|---|---|
| Paper (source of truth) | `paper/n57_paper.md` | 1 | Full proof, $n=5$ and $n=7$ |
| Exact symbolic certificates | `code/verify_*.py` | 8 | N–Y spectral bridge, $P=C+2\sqrt{AB}$, $S_1$ elimination completeness, $\beta$-reduction, uniform Hessian, $S_3$ closed form, GPT-formula cross-checks |
| $n=5$ certificates | `code/n5_*.py` | 2 | stationary curve, $M(1)$ interval |
| $n=7$ certificates | `code/n7_*.py` | 145 | $S_2$ resultant/Sturm/rational-signs, $S_1$ crossing/monotonicity/determinant/three-sample, $H_C$ critical-events/cover/Krawczyk/Arb, inactive-KKT, etc. |
| Machine-checkable cover | `code/_hc_cover.json` | 1 | 2604-piece $H_C$ validated cover (rational $s/c/v$ boxes, admissibility lower bounds, $L_C=141/20$) |
| Completeness summary | `code/_hc_completeness.json` | 1 | per-cell Sturm counts (`all_cells_ok=True`) |
| Verifier/checker run output | `code/_archival_run.txt` | 1 | captured stdout of all 5 pass scripts (this run) |
| Dependency lock | `code/requirements_lock.txt` | 1 | exact package versions |
| Code README | `code/README.md` | 1 | per-script index |
| Integrity manifest | `sha256_manifest.txt` | 1 | SHA-256 of every archived file (this file excluded) |
| This document | `ARCHIVAL.md` | 1 | this file |

**Excluded (regenerable or out of scope):** `__pycache__`, `*.pickle`, browser-act/ChatGPT driver scripts (`_paste*.js`, `_send*.js`, …) and session logs/screenshots, exploration `*_out.txt`/`*_fb.txt` logs (superseded by `_archival_run.txt`), `_hc_cover_refined.json` (duplicate of `_hc_cover.json`), and the unrelated number-theory-falsification scripts (`a295196_*`, `agrawal_*`, `classify_*`, `classical_shapiro_check.py`). The `.tex`/`.pdf` in `paper/` are excluded as stale; regenerate from `n57_paper.md` before final submission.

---

## 2. Integrity (SHA-256)

`sha256_manifest.txt` lists `sha256  path  (bytes)` for all 162 archived files. Key entries:

```
53a518fe099591a43a7da2ecc31e130ad61a6c0eebb9833f3a51ed8192d054bd  paper/n57_paper.md            (63243 bytes)
d84ba16037baa621c5da01bbc54c418dda08c9324707795efb26b236089f4541  code/_hc_cover.json           (1237848 bytes)
cc504b280272e8b68d2a8fa79a29e1a64f1b3a6e7625b20962bf607e08b4146e  code/_hc_completeness.json    (3230 bytes)
001bfb29e9611bc2c2eac7fd9db318bbc7aad3bdcbdf42f075d75bfdd0de6e80  code/requirements_lock.txt    (709 bytes)
d719d65ad3b1f4eeb489f9ae89857f8fc977cda05de90869c7f646458a9f7276  code/_archival_run.txt        (2394 bytes)
4ea6bb42524b5686e732309e2c03b6722618da3f5ceeaeeb21b1154baf38f271  code/README.md                (3811 bytes)
09c944f621f26df4f0bed5323e4370241a437df5cb56d119e540d2a784225913  .gitignore                    (1257 bytes)
```

Verify on any checkout: `sha256sum -c sha256_manifest.txt` (the manifest does not include itself or `ARCHIVAL.md`, by the standard no-self-reference convention).

---

## 3. Reproduce (5 scripts, all pass)

Environment: Python 3.13.9, Windows 11. `pip install -r code/requirements_lock.txt`. Run from the **repository root**:

```bash
python code/verify_ny_spectral.py            # EXACT: Lemma 2.2 N-Y spectral bridge (μ_k identity, disk equiv, Δ_n<0)
python code/verify_hc_closedform.py          # EXACT: P = C + 2√(AB)  (Issue 6)
python code/verify_s1_elimination.py         # EXACT: Res_v(F1,F2) = w²z²(w-1)² H_B H_C  (Issue 7)
python code/n7_s1_hc_cover_checker.py        # INTERVAL (mpmath.iv): 2604/2604 pieces, P_MV > L_C = 141/20
python code/n7_s1_hc_arb_checker.py          # INTERVAL (independent Arb/python-flint): 2604/2604, global min 7.050005 > 7
```

Expected: each prints `DONE-*` and (for the checkers) `ALL PIECES ... VERIFIED ... : True`. The captured output of one full run is in `code/_archival_run.txt`.

The three `EXACT` scripts use only SymPy rational algebra (no floating point in any load-bearing step); the two interval scripts use outward-rounded `mpmath.iv` and an independent hand-written Arb interval layer respectively, and agree on 2604/2604 pieces.

---

## 4. Version control

This is a self-contained local git repository (independent of any parent repo). The proof artifacts are committed in a single snapshot whose SHA-256 manifest is §2.

- **Proof-artifact commit:** `<filled by the doc-pointer commit — see git log>` (contains every file listed in §1 and hashed in `sha256_manifest.txt`).
- **Authoritative hash on any checkout:** `git rev-parse HEAD` (the doc-pointer commit on top of the proof-artifact commit only updates this file's hash field; the proof artifacts are unchanged, so `sha256_manifest.txt` is identical).
- **Zenodo DOI:** `<to be minted on deposition>`

---

## 5. Dependency lock

`code/requirements_lock.txt` pins: `sympy==1.14.0`, `mpmath==1.3.0`, `numpy==2.4.4`, `scipy==1.18.0`, `python-flint==0.9.0` (Python 3.13.9). The EXACT certificates are version-robust across SymPy ≥ 1.12; the interval certificates are pinned for bit-reproducibility.
