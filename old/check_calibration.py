#!/usr/bin/env python
"""
check_calibration.py - Is the estimator actually running, and what does it say?

Answers three questions without running the pipeline:

  1. Does the installed calibrate.py contain the per-cluster estimator at all?
  2. Do the installed files match the shipped checksums?
  3. Given a real fragments.tsv, what threshold would the current code choose,
     and how many fragments would qualify?

Question 3 is the useful one. It reads the debug table from a previous run and
reports what the installed code would do with it, so the answer does not depend
on snakemake having decided to rerun anything.

Usage:
    python check_calibration.py --eukdetect_dir /path/to/EukDetect
    python check_calibration.py --eukdetect_dir /path/to/EukDetect \\
        --fragments /path/to/ANIC000643_fragments.tsv
"""

import argparse
import hashlib
import importlib.util
import os
import sys
from collections import defaultdict

EXPECTED = {
    "eukdetect/util/calibrate.py": "estimate_anchor_threshold",
    "eukdetect/util/realign.py": "substitutions",
    "eukdetect/util/eukfrac_calc.py": "cluster_threshold",
}


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eukdetect_dir", required=True,
                   help="Repo root, the directory containing eukdetect/")
    p.add_argument("--fragments",
                   help="A *_fragments.tsv from a realignment debug directory")
    p.add_argument("--manifest", help="manifest.md5 to compare against")
    p.add_argument("--config",
                   help="The config .yml the run used. The code defaults mean "
                        "nothing if the config overrides them, which is the "
                        "single most common reason output does not change.")
    args = p.parse_args(argv)

    root = args.eukdetect_dir
    fail = False

    print("== 1. is the estimator present in the installed code ==")
    for rel, marker in EXPECTED.items():
        path = os.path.join(root, rel)
        if not os.path.exists(path):
            print("  MISSING  %s" % rel)
            fail = True
            continue
        text = open(path).read()
        n = text.count(marker)
        ok = n > 0
        fail = fail or not ok
        print("  %-8s %-34s '%s' x%d" % ("OK" if ok else "STALE", rel, marker, n))
    print()

    if args.manifest and os.path.exists(args.manifest):
        print("== 2. checksums ==")
        want = {}
        for line in open(args.manifest):
            parts = line.split()
            if len(parts) == 2:
                want[parts[1]] = parts[0]
        for rel in EXPECTED:
            path = os.path.join(root, rel)
            if not os.path.exists(path):
                continue
            got = md5(path)
            exp = want.get(rel)
            if exp is None:
                print("  ?        %s not in manifest" % rel)
            elif got == exp:
                print("  OK       %s" % rel)
            else:
                print("  DIFFERS  %s\n             installed %s\n             expected  %s"
                      % (rel, got, exp))
                fail = True
        print()

    print("== 3. what the installed code would do ==")
    cal = os.path.join(root, "eukdetect", "util", "calibrate.py")
    if not os.path.exists(cal):
        print("  cannot load calibrate.py")
        return 1
    spec = importlib.util.spec_from_file_location("calibrate", cal)
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)

    if not hasattr(C, "estimate_anchor_threshold"):
        print("  estimate_anchor_threshold is ABSENT: this is the old code.")
        print("  Every fragment below the fixed threshold falls through to the")
        print("  spillover comparison, which is the behaviour you are seeing.")
        return 1

    floor = C.DEFAULT_ANCHOR_IDENTITY
    quant = C.DEFAULT_ANCHOR_QUANTILE
    slack = C.DEFAULT_ANCHOR_SLACK
    source = "code default"

    if args.config and os.path.exists(args.config):
        try:
            import yaml
            cfg = yaml.safe_load(open(args.config)) or {}
        except Exception as e:
            cfg = {}
            print("  could not read config (%s)" % e)
        if "anchor_identity" in cfg:
            floor = float(cfg["anchor_identity"])
            source = "config: %s" % args.config
        quant = float(cfg.get("anchor_quantile", quant))
        slack = float(cfg.get("anchor_slack", slack))

    print("  estimate_anchor_threshold present")
    print("  floor (anchor_identity)         = %.2f   [%s]" % (floor, source))
    print("  quantile                        = %.2f" % quant)
    print("  slack                           = %.2f" % slack)
    if floor >= 0.99:
        print()
        print("  PROBLEM: the floor is %.2f. anchor_identity is a FLOOR now, not" % floor)
        print("           the threshold. At this value the estimate can never")
        print("           lower the bar, so the run behaves exactly as it did")
        print("           before the estimator existed.")
        print("           Set anchor_identity: 0.97 or remove it from the config.")
        fail = True
    print()

    if not args.fragments:
        print("  (pass --fragments to see what it would do on your data)")
        return 1 if fail else 0

    if not os.path.exists(args.fragments):
        print("  fragments file not found: %s" % args.fragments)
        return 1

    # Rebuild the per-fragment structure the estimator expects.
    frag = {}
    by_tax = defaultdict(int)
    with open(args.fragments) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_q, i_tax = header.index("Fragment"), header.index("Taxid")
            i_id, i_best = header.index("Identity"), header.index("Is_best")
            i_win = header.index("Is_winner")
        except ValueError as e:
            print("  unexpected columns in fragments file: %s" % e)
            return 1
        for line in f:
            r = line.rstrip("\n").split("\t")
            if len(r) <= max(i_q, i_tax, i_id, i_best, i_win):
                continue
            q = r[i_q]
            rec = frag.setdefault(q, {"winners": [], "identity": {},
                                      "scores": {}, "best_marker": {}})
            try:
                ident = float(r[i_id]) / 100.0
            except ValueError:
                continue
            rec["identity"][r[i_tax]] = ident
            rec["best_marker"][r[i_tax]] = "m"
            if r[i_win] == "yes":
                rec["winners"].append(r[i_tax])
            if r[i_best] == "yes":
                by_tax[r[i_tax]] += 1
    for rec in frag.values():
        rec["winners"] = sorted(set(rec["winners"]))

    # The floor that the RUN would use, which is not necessarily the module
    # default: a config key overrides it, and an old value silently disables
    # the estimate.
    th, level, n = C.estimate_anchor_threshold(
        frag, quantile=quant, slack=slack, floor=floor)
    print("  fragments read                  %d" % len(frag))
    print("  estimated within-species level  %.2f%%" % (level * 100 if level else 0))
    print("  anchor threshold it would use   %.2f%%" % (th * 100))
    if level is not None and th <= floor:
        print("                                  (raised to the floor)")
    print()

    qual = defaultdict(int)
    old = defaultdict(int)
    for rec in frag.values():
        if len(rec["winners"]) != 1:
            continue
        t = rec["winners"][0]
        i = rec["identity"].get(t)
        if i is None:
            continue
        if i >= th:
            qual[t] += 1
        if i >= 0.99:
            old[t] += 1
    print("  %-12s %10s %14s %14s" % ("taxid", "fragments", "qualify now", "qualified at .99"))
    for t in sorted(by_tax, key=lambda x: -by_tax[x]):
        print("  %-12s %10d %14d %14d" % (t, by_tax[t], qual[t], old[t]))
    print()
    print("  'qualified at .99' should match the Unique_reads in your last")
    print("  filtered_hits table. 'qualify now' is what the current code gives.")
    print("  If your table still shows the .99 numbers, taxonomize did not rerun:")
    print("      snakemake ... --forcerun taxonomize")
    print("  or with the wrapper, add --force")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
