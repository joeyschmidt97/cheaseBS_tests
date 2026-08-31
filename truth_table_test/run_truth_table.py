#!/usr/bin/env python3
"""Score completed cheaseBS runs against the profile-scaling truth table.

Point it at run directories that already exist -- it never launches CHEASE.

    # one point
    python run_truth_table.py --run /path/to/Te_ped_scale_1.300

    # one campaign: every point under it, plus the cross-case ordering checks
    python run_truth_table.py --campaign /path/to/runs/129015_20260824_20-28-57

    # a whole tree of campaigns: reshape_convergence/runs/ holds one directory
    # per shot per timestamp, and each of those holds the four box-edge points
    python run_truth_table.py --campaign /path/to/reshape_convergence/runs --summary-only

    # score the documented NSTX replay point instead of the final iteration
    python run_truth_table.py --campaign <dir> --iteration 00

--campaign is recursive. It walks the tree, treats any directory carrying an
EQDSK*.OUT / iteration_NN / cheasebs_run_config.json as one run point, and does
not descend into a point's own artifact subdirectories. Points are grouped by
their parent directory, so a tree of four shots gets one ordering block per
campaign rather than one meaningless block spanning all of them.

--iteration only chooses which equilibrium the control and response blocks are
scored on. The convergence block always reads the whole path: residuals against
the run's own tol_ip_rel / tol_bs / tol_q / tol_a, ringing versus saturation,
and the iteration-00-to-final drift.

Exit status is 0 when nothing FAILed, 1 otherwise, so a wrapper can tell "ran"
from "verified". FLAG never fails the run: it marks the rows that need a human.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import checks  # noqa: E402
import truth_table as tt  # noqa: E402

# Subdirectories of a run point, or of a campaign, that are never themselves a
# run point. Pruning them keeps the walk cheap on a tree carrying hundreds of
# per-iteration artifact directories.
PRUNE = {"artifacts", "cheasebs_baseline", "reference_profiles", "__pycache__",
         ".git", ".ipynb_checkpoints", "seed_profiles", "_reconstruction",
         "chease_mapping_artifacts", "chease_stage2_artifacts",
         "bootstrap_diagnostic_artifacts"}

# Campaign records that carry the axis and scale of each point -- which is what
# tells the checker which way each response is supposed to move.
CAMPAIGN_JSONS = ("reshape_convergence.json", "scale_symmetry_results.json")


def _is_run_point(path):
    return bool(glob.glob(os.path.join(path, "EQDSK*.OUT"))
                or glob.glob(os.path.join(path, "iteration_*"))
                or os.path.isfile(os.path.join(path, "cheasebs_run_config.json")))


def _campaign_meta(campaign_dir):
    """{abs point dir: {axis, scale, error}} from a campaign's own JSON record.

    The recorded savedir paths are absolute and usually from the machine the
    campaign ran on, so the same-named local subdirectory wins when it exists.
    """
    meta = {}
    for name in CAMPAIGN_JSONS:
        jpath = os.path.join(campaign_dir, name)
        if not os.path.isfile(jpath):
            continue
        try:
            rec = json.load(open(jpath))
        except Exception as exc:                                # noqa: BLE001
            print("warning: %s unreadable: %s" % (jpath, exc))
            continue
        for r in rec.get("rows", []):
            saved = str(r.get("savedir") or "")
            if not saved:
                continue
            local = os.path.join(campaign_dir, os.path.basename(saved))
            path = local if os.path.isdir(local) else saved
            meta[os.path.abspath(path)] = {"axis": r.get("axis"),
                                           "scale": r.get("scale"),
                                           "error": r.get("error"),
                                           "campaign": campaign_dir}
    return meta


def discover_cases(root):
    """Every run point under root, at any depth, with its axis/scale when known.

    Handles all three shapes without being told which one it got: a single run
    point, a campaign of points, and a tree of campaigns.
    """
    root = os.path.abspath(root)
    if _is_run_point(root):
        m = _campaign_meta(os.path.dirname(root)).get(root, {})
        return [{"path": root, "axis": m.get("axis"), "scale": m.get("scale"),
                 "campaign": os.path.dirname(root)}]

    cases, meta, seen = [], {}, set()
    for dirpath, dirnames, _ in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in PRUNE and not d.startswith("iteration_")]
        if any(os.path.isfile(os.path.join(dirpath, n)) for n in CAMPAIGN_JSONS):
            meta.update(_campaign_meta(dirpath))
        keep = []
        for d in sorted(dirnames):
            path = os.path.abspath(os.path.join(dirpath, d))
            if _is_run_point(path):
                cases.append({"path": path, "campaign": dirpath})
                seen.add(path)
            else:
                keep.append(d)                # not a point: keep walking into it
        dirnames[:] = keep

    for c in cases:
        m = meta.get(c["path"], {})
        c["axis"], c["scale"] = m.get("axis"), m.get("scale")

    # Points a campaign recorded that are not on this machine. Reported rather
    # than skipped: a scan with a missing edge is not a scan.
    for path, m in sorted(meta.items()):
        if path not in seen and not os.path.isdir(path):
            cases.append({"path": path, "axis": m.get("axis"),
                          "scale": m.get("scale"),
                          "campaign": m.get("campaign", root), "missing": True})
    return cases


def _tol_overrides(args):
    out = {}
    if args.tol_ip_rel is not None:
        out["ip_rel"] = args.tol_ip_rel
    if args.dead is not None:
        out["dead"] = args.dead
    return out


def score(run_dir, args, axis=None, scale=None):
    """Locate, build and return (files, rows) for one run directory."""
    files = tt.locate(run_dir, source_eqdsk=args.source_eqdsk,
                      iteration=args.iteration)
    if axis is not None:
        files["case"]["axis"] = axis
    if scale is not None:
        files["case"]["scale"] = scale
    return files, checks.build_rows(files, tol=_tol_overrides(args))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Truth table for completed cheaseBS profile-scaling runs.")
    ap.add_argument("--run", action="append", default=[],
                    help="a single completed run directory; repeatable")
    ap.add_argument("--campaign", action="append", default=[],
                    help="a campaign directory, or a tree of them; walked "
                         "recursively; repeatable")
    ap.add_argument("--iteration", default="final",
                    help="which equilibrium the control/response blocks score: "
                         "'final' (default) or an index such as 00. The "
                         "convergence block always reads the whole path")
    ap.add_argument("--source-eqdsk", default=None,
                    help="source EQDSK to compare against, when the copy in the "
                         "run directory is missing")
    ap.add_argument("--tol-ip-rel", type=float, default=None,
                    help="override the Ip tolerance (default: the run's own "
                         "tol_ip_rel, else 0.02)")
    ap.add_argument("--dead", type=float, default=None,
                    help="relative change below which a response counts as "
                         "'never moved' (default 1e-3)")
    ap.add_argument("--out", default=None,
                    help="write every row as JSON here (default: "
                         "truth_table.json at the first target)")
    ap.add_argument("--quiet-notes", action="store_true",
                    help="table only, no per-row explanation")
    ap.add_argument("--summary-only", action="store_true",
                    help="one line per run plus the cross-case blocks; for a "
                         "tree of campaigns, where full tables run to pages")
    args = ap.parse_args(argv)

    if not args.run and not args.campaign:
        ap.error("nothing to check: pass --run and/or --campaign")

    targets = []
    for p in args.run:
        path = os.path.abspath(p)
        targets.append({"path": path, "axis": None, "scale": None,
                        "campaign": os.path.dirname(path)})
    for c in args.campaign:
        found = discover_cases(c)
        if not found:
            print("no run directories found under %s" % c)
        targets.extend(found)

    scored = []
    totals = {tt.PASS: 0, tt.FLAG: 0, tt.FAIL: 0, tt.SKIP: 0}

    def bump(counts):
        for k, v in counts.items():
            totals[k] = totals.get(k, 0) + v

    for t in targets:
        path = t["path"]
        label = os.path.basename(path.rstrip("/\\")) or path
        if not os.path.isdir(path):
            print("\n=== %s ===" % label)
            print("  directory missing: %s" % path)
            print("  a campaign record lists this point but it is not on this "
                  "machine; copy the run directory here, or run the checker there")
            totals[tt.FAIL] += 1
            scored.append({"label": label, "path": path, "campaign": t["campaign"],
                           "axis": t.get("axis"), "scale": t.get("scale"),
                           "rows": [], "missing": True})
            continue

        files, rows = score(path, args, t.get("axis"), t.get("scale"))
        counts = checks.tally(rows)
        bump(counts)
        line = "  ".join("%s %d" % (k, counts.get(k, 0))
                         for k in (tt.PASS, tt.FLAG, tt.FAIL, tt.SKIP))
        if args.summary_only:
            failed = [r["check"] for r in rows if r["verdict"] == tt.FAIL]
            print("%-34s axis=%-14s scale=%-5s %s%s"
                  % (label[:34], files["case"]["axis"], files["case"]["scale"],
                     line, ("   FAIL: " + ", ".join(failed)) if failed else ""))
        else:
            print("\n=== %s ===" % label)
            print("  axis=%s scale=%s  scored iteration=%s"
                  % (files["case"]["axis"], files["case"]["scale"],
                     files["iteration_index"]))
            for n in files["notes"]:
                print("  note: %s" % n)
            print(checks.render(rows, notes=not args.quiet_notes))
            print("  " + line)
        scored.append({"label": label, "path": path, "campaign": t["campaign"],
                       "axis": files["case"]["axis"],
                       "scale": files["case"]["scale"],
                       "iteration": files["iteration_index"], "rows": rows})

    # Cross-case checks per campaign: a tree of four shots must not be scored as
    # one scan.
    cross, groups = {}, {}
    for s in scored:
        if s["rows"]:
            groups.setdefault(s["campaign"], []).append(s)
    for campaign, group in sorted(groups.items()):
        if len(group) < 2:
            continue
        rows = checks.campaign_rows(group)
        if not rows:
            continue
        cross[campaign] = rows
        print("\n=== across cases: %s ==="
              % (os.path.basename(campaign.rstrip("/\\")) or campaign))
        print(checks.render(rows, notes=not args.quiet_notes))
        counts = checks.tally(rows)
        bump(counts)
        print("  " + "  ".join("%s %d" % (k, counts.get(k, 0))
                               for k in (tt.PASS, tt.FLAG, tt.FAIL, tt.SKIP)))

    print("\n=== totals over %d run(s) ===" % len(scored))
    print("  " + "  ".join("%s %d" % (k, totals.get(k, 0))
                           for k in (tt.PASS, tt.FLAG, tt.FAIL, tt.SKIP)))
    if totals.get(tt.FAIL):
        print("  FAIL rows present: read those before taking any physics from "
              "this scan. Controls first, then responses.")
    elif totals.get(tt.FLAG):
        print("  no failures; FLAG rows need a human (edge-q offset, iteration "
              "cap, ringing residual, sub-grid-noise response).")

    out = args.out
    if out is None:
        base = os.path.abspath(args.campaign[0] if args.campaign
                               else targets[0]["path"])
        root = base if os.path.isdir(base) else os.path.dirname(base)
        out = os.path.join(root, "truth_table.json")
    payload = {"scored_iteration": args.iteration, "totals": totals,
               "cases": scored, "across_cases": cross}
    try:
        with open(out, "w") as fh:
            json.dump(payload, fh, indent=1, default=str)
        print("  wrote %s" % out)
    except OSError as exc:
        print("  could not write %s: %s" % (out, exc))

    return 1 if totals.get(tt.FAIL) else 0


if __name__ == "__main__":
    sys.exit(main())
