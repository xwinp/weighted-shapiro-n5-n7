#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run every Tier 1 + Tier 2 *certificate* script and assert a zero exit code.

Tier 1 (load-bearing core, 8) and Tier 2 (formal theorem certificates) are the
scripts whose pass/fail determines the theorems; Tier 3 (diagnostic /
corroboration / exploratory) is intentionally excluded.  See ARCHIVAL.md ->
"Script tiering" for the authoritative 3-tier list.

Record-*generators* (Tier 2 build scripts that WRITE committed records and would
overwrite the canonical cover) are also excluded -- their output is verified by
the Tier 1 checkers (n7_s1_hc_cover_checker.py / n7_s1_hc_arb_checker.py):
  n7_s1_hc_cover_dump.py, n7_s1_hc_arb_refine.py  (build _hc_cover.json)

Usage (from the repository root):
    python code/run_all_certificates.py            # Tier 1 + Tier 2 (default)
    python code/run_all_certificates.py --tier 1   # Tier 1 only
    python code/run_all_certificates.py --tier 2   # Tier 2 only

Exits 0 iff every requested certificate script exits 0.
"""
import sys
import os
import subprocess

# (tier, script) -- certificate scripts only (no record-generators, no Tier 3)
CERTIFICATES = [
    # ---- Tier 1: load-bearing core (7) ----
    (1, "verify_ny_spectral.py"),
    (1, "verify_hc_closedform.py"),
    (1, "verify_s1_elimination.py"),
    (1, "verify_s1_hc_superset.py"),
    (1, "n7_s1_hc_s1collar.py"),
    (1, "n7_s1_hc_exhaustiveness_cert.py"),
    (1, "n7_s1_hc_cover_checker.py"),
    (1, "n7_s1_hc_arb_checker.py"),
    # ---- Tier 2: formal theorem certificates ----
    (2, "n5_s1_certificate.py"),
    (2, "n5_s1_m1_interval.py"),
    (2, "n7_resultant.py"),
    (2, "verify_endpoints.py"),
    (2, "n7_s2_rational_signs.py"),
    (2, "n7_inactive_sturm.py"),
    (2, "verify_s3_closed.py"),
    (2, "verify_uniform_hess.py"),
    (2, "n7_s1_crossing_resultant.py"),
    (2, "_verify_crossings_in_band.py"),
    (2, "_verify_s1_three_samples.py"),
    (2, "n7_s1_rho_monotone_Lcert.py"),
    (2, "n7_s1_rigorous_certs.py"),
    (2, "n7_s1_hc_kkt_check.py"),
]

def main():
    tier_filter = None
    if "--tier" in sys.argv:
        i = sys.argv.index("--tier")
        if i + 1 < len(sys.argv):
            tier_filter = int(sys.argv[i + 1])

    scripts = [(t, s) for (t, s) in CERTIFICATES if tier_filter is None or t == tier_filter]
    print("run_all_certificates.py -- running %d certificate script(s) (tier_filter=%s)\n" % (
        len(scripts), tier_filter if tier_filter else "all"))

    failures = []
    for tier, script in scripts:
        label = "T%d %s" % (tier, script)
        print("=== %s ===" % label, flush=True)
        try:
            # Force UTF-8 stdio in the child so non-ASCII output (math symbols in
            # verdict lines) survives on CJK Windows consoles (cp936/gbk).
            env = dict(os.environ)
            env["PYTHONUTF8"] = "1"
            r = subprocess.run([sys.executable, "code/" + script],
                               capture_output=True, text=True,
                               encoding="utf-8", errors="replace",
                               env=env, timeout=900)
        except subprocess.TimeoutExpired:
            print("  TIMEOUT (>900s) -- FAIL")
            failures.append((tier, script, "timeout"))
            continue
        # echo the tail of stdout so the certificate's verdict line is visible
        tail = "\n".join(r.stdout.strip().splitlines()[-3:])
        print("  exit=%d" % r.returncode)
        if tail:
            print("  " + tail.replace("\n", "\n  "))
        if r.returncode != 0:
            err_tail = "\n".join(r.stderr.strip().splitlines()[-5:])
            if err_tail:
                print("  [stderr] " + err_tail.replace("\n", "\n  "))
            failures.append((tier, script, "exit=%d" % r.returncode))
        print()

    print("=" * 60)
    n_t1 = sum(1 for t, s in scripts if t == 1)
    n_t2 = sum(1 for t, s in scripts if t == 2)
    print("Ran %d Tier 1 + %d Tier 2 certificate scripts." % (n_t1, n_t2))
    if not failures:
        print("ALL PASS (every certificate exited 0).")
        sys.exit(0)
    print("FAILURES (%d):" % len(failures))
    for tier, script, why in failures:
        print("  T%d %s : %s" % (tier, script, why))
    sys.exit(1)

if __name__ == "__main__":
    main()
