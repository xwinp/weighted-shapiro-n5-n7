# Weighted Shapiro cyclic inequality (n = 5, 7) — proof artifact

This repository is the computer-assisted proof artifact for the paper
**"The Tuan–Thuong weighted Shapiro cyclic inequality for n = 5 and n = 7"**
by Weipeng Xue (Sun Yat-sen University).

It proves that the Tuan–Thuong weighted Shapiro cyclic inequality holds for
`n = 5` (the whole band `H_5 = (0,1)`) and for `n = 7`
(`H_7 = (0, a_7] ∪ [b_7, 1)` with `a_7 ≈ 0.21427352091`, `b_7 ≈ 0.32862767792`),
i.e. `P_{n,p,q}(x) ≥ n` throughout the holding region, with equality only at
the uniform point `x = 1`.

## What is here

- `paper/n57_paper.md` — the manuscript (source of truth).
- `paper/n57_paper.pdf` — the typeset manuscript (canonical artifact).
- `code/` — the verification scripts and machine records.
- `SUPPLEMENT.md` — reader-facing software environment, per-file listings,
  the `Φ_35` coefficient vector, and reproducibility commands.
- `code/README.md` — the authoritative three-tier script role map.
- `code/requirements_lock.txt` — pinned dependencies.
- `sha256_manifest.txt` — SHA-256 of every tracked file.
- `ARCHIVAL.md` — revision history, tier classification, build/deposition policy.

## Verify the theorems (one command)

```bash
pip install -r code/requirements_lock.txt
python code/run_all_certificates.py
```

This runs **22 certificates** (8 Tier 1 load-bearing + 14 Tier 2 formal theorem
certificates) and asserts a zero exit code for each. Expected output ends with
`ALL PASS (every certificate exited 0).` Key invariants it certifies:

| invariant | value |
|---|---|
| `S_1` collar `g_min` | `7.344054 > L_C` |
| `H_C` exhaustiveness | `components=2, cells=6, fibers=5, edges=2823` |
| `H_C` Arb cover | `2604/2604` boxes, `min 7.050005 > 141/20` |

The captured run log is `code/_archival_run.txt`; the machine-readable cover
records are `code/_hc_cover.json`, `code/_hc_s1collar.json`,
`code/_hc_exhaustiveness.json`.

## Version of record

The submission version is the git tag **`v1.0-submission`**:

```bash
git rev-parse v1.0-submission      # the exact commit hash
git checkout v1.0-submission
```

The complete archive (including the packed `.git` history) is deposited on
Zenodo (DOI recorded in `SUPPLEMENT.md` / `ARCHIVAL.md` upon deposition); the
archived artifact is the authoritative version associated with the manuscript.

## Arithmetic backends

Two independent, outward-directed backends are used: **exact rational**
(SymPy) and **interval arithmetic** (`mpmath.iv` and Arb via `python-flint`).
The `H_C` cover is verified by both backends independently. No load-bearing
step uses floating point.

## AI assistance

Some computer-assisted verification scripts were developed with the assistance
of a large language model (ChatGPT) and were adopted only after the author's
independent audit and exact recomputation; all mathematical content, theorems,
and proofs are the author's own. See `paper/n57_paper.md` §5.
