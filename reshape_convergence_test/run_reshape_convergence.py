"""Reshape limits and cheaseBS convergence, per discharge, per variable.

Straight discharge object -> transform -> cheaseBS. No sparse grid, no campaign,
no sampler: those add refinement logic between the reshape and the solve, and
the question here is only about the solve.

WHAT THIS ANSWERS

1. How many iterations does cheaseBS need, as a function of how hard the profile
   was reshaped? The 2026-08-18 five-point 132543 run showed four points
   converging in 2 iterations and one (ne x1.3) needing 19 — a 10x spread with
   no explanation. If iteration count rises smoothly with |scale - 1| the solver
   is behaving; if one variable or one direction is erratic, it is not.

2. Where is the reshape limit? Scales run wider than the campaign box on
   purpose. The limit shows up as the scale at which cheaseBS stops converging
   inside max_iter, the acceptance gate starts rejecting, or the solve raises.
   Nobody has established that boundary, so the campaign bounds are currently
   set from literature alone with no check that the reconstruction can follow.

3. Is cheaseBS set up correctly at all? A monotone, symmetric, converging
   response across four discharges is the assurance. Anything else is a finding.

WHAT IT DOES NOT ANSWER

Nothing about GENE, growth rates, or whether an equilibrium is a good input for
a linear run. It is a solver-behaviour test.

Usage (NERSC):
    python run_reshape_convergence.py
    python run_reshape_convergence.py --shots 132543 --vars ne
    python run_reshape_convergence.py --scales 0.5,0.7,0.9,1.1,1.3,1.5
    python run_reshape_convergence.py --allow-missing-chease   # harness check only
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime

import numpy as np

from TPED.projects.discharge_tools.src.discharge_data import DischargeData
from TPED.projects.discharge_tools.src.discharge_physics import DischargePhysics
from TPED.projects.discharge_tools.src.cheasebs_runner import (
    CheasebsAcceptance, read_acceptance)
from TPED.projects.discharge_tools.src.transforms.mtanh_transforms import fit_mtanh

DEFAULT_ROOT = "/global/homes/j/joeschm/data/ST_research/NSTXU_discharges"

# Wider than the campaign box (0.7-1.3) on purpose: the reshape limit cannot be
# found from inside the range already assumed to work. Symmetric about 1.0 so
# each +delta has a -delta partner.
DEFAULT_SCALES = (0.5, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5)

# analysis_radii are where the gate scores q -- the radii GENE will run at.
# 129038 has no scan of record, so its entry is a placeholder and its gate
# verdict should not be read as meaningful.
DISCHARGES = {
    129015: {"radii": (0.85,), "placeholder_radii": False, "extra": {}},
    129038: {"radii": (0.85,), "placeholder_radii": True,
             "extra": {"pfile_name": "p129038.00400"}},
    132543: {"radii": (0.736, 0.825), "placeholder_radii": False, "extra": {}},
    132588: {"radii": (0.736, 0.825), "placeholder_radii": False, "extra": {}},
}

# Same settings the campaign uses, so the numbers transfer directly.
CHEASEBS = {
    "max_iter": 25,
    "tol_bs": 1e-4,
    "tol_q": 1e-4,
    "tol_ip_rel": 0.02,
    "bootstrap_mix": 0.1,
    "istar_mix": 0.05,
    "plot_errors": True,
}

PED = (0.90, 0.98)
MTANH_FIT_KWARGS = dict(pedestal_weight=8.0)


def ped_top(x, y):
    m = (x >= PED[0]) & (x <= PED[1])
    return float(np.max(y[m])) if m.any() else float("nan")


def load(shot, root):
    """Build the DischargeData, honouring per-shot explicit files.

    129038's directory holds five pfiles, so auto-discovery refuses to guess --
    correctly, since picking one silently would put an unannounced rotation
    variant under the whole test.
    """
    d = os.path.join(root, str(shot))
    kwargs = {"input_dir": d}
    name = DISCHARGES[shot]["extra"].get("pfile_name")
    if name:
        kwargs["pfile"] = os.path.join(d, name)
    return DischargeData(**kwargs)


def run_one(phys0, fit, var, scale, savedir, radii, allow_missing):
    """One reshape, one cheaseBS solve. Returns the row for the results table."""
    os.makedirs(savedir, exist_ok=True)
    row = {"var": var, "scale": scale, "savedir": savedir}

    phys = phys0.apply_mtanh_ped(var, fit=fit, scale_height=scale,
                                 enforce_quasineutrality=True, qz=6.0)
    x = np.asarray(phys.rhot.values)
    row["ped_top"] = ped_top(x, np.asarray(getattr(phys, var).values))

    t0 = time.time()
    try:
        phys.output_gfile(savedir=savedir, run_cheasebs=True,
                          cheasebs_acceptance=CheasebsAcceptance.production(
                              analysis_radii=radii),
                          cheasebs_strict=False, comment=f"{var}_{scale:.3f}",
                          **CHEASEBS)
    except Exception as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["wall_s"] = time.time() - t0
        if not allow_missing:
            # A solve that raises IS a reshape-limit datapoint, so it is
            # recorded and the sweep continues rather than aborting.
            print(f"      RAISED: {row['error']}")
        return row
    row["wall_s"] = time.time() - t0

    rec = read_acceptance(savedir) or {}
    summary_path = os.path.join(savedir, "convergence_summary.json")
    summary = json.load(open(summary_path)) if os.path.exists(summary_path) else {}
    row.update({
        "accepted": rec.get("accepted"),
        "reasons": rec.get("reasons"),
        "ip_error_rel": rec.get("ip_error_rel"),
        "q_errors_rel": rec.get("q_errors_rel"),
        "q_edge_error_rel": rec.get("q_edge_error_rel"),
        "iterations": summary.get("iterations"),
        "converged": summary.get("final_converged"),
        "bootstrap_change_rel": summary.get("final_bootstrap_change_rel"),
        "q_change_rel": summary.get("final_q_change_rel"),
        "final_ip_a": summary.get("final_ip_a"),
        "target_ip_a": summary.get("target_ip_a"),
    })
    return row


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--outroot", default=os.path.join(os.path.dirname(__file__), "runs"))
    ap.add_argument("--shots", default=",".join(str(s) for s in DISCHARGES))
    ap.add_argument("--vars", default="ne,Te")
    ap.add_argument("--scales", default=",".join(str(s) for s in DEFAULT_SCALES))
    ap.add_argument("--allow-missing-chease", action="store_true",
                    help="keep going when cheaseBS cannot run, to check the "
                         "harness itself. Produces no convergence data.")
    args = ap.parse_args()

    shots = [int(s) for s in args.shots.split(",")]
    variables = [v.strip() for v in args.vars.split(",")]
    scales = [float(s) for s in args.scales.split(",")]
    stamp = datetime.now().strftime("%Y%m%d_%H-%M-%S")
    workdir = os.path.abspath(os.path.join(args.outroot, stamp))
    os.makedirs(workdir, exist_ok=True)

    print(f"[test] workdir {workdir}")
    print(f"[test] shots   {shots}")
    print(f"[test] vars    {variables}")
    print(f"[test] scales  {scales}")
    print(f"[test] cheaseBS {CHEASEBS}")

    results = {}
    for shot in shots:
        print("\n" + "=" * 74 + f"\n  {shot}\n" + "=" * 74)
        try:
            phys0 = DischargePhysics(load(shot, args.root))
        except Exception as exc:
            print(f"  setup FAILED: {type(exc).__name__}: {exc}")
            results[shot] = {"error": f"{type(exc).__name__}: {exc}"}
            continue

        x0 = np.asarray(phys0.rhot.values)
        radii = DISCHARGES[shot]["radii"]
        rows, fits = [], {}
        for var in variables:
            fit, rec = fit_mtanh(phys0.ds, var, **MTANH_FIT_KWARGS)
            fits[var] = {"rms_relative": rec["rms_relative"],
                         "fit_params": rec["fit_params"],
                         "ped_top": ped_top(x0, np.asarray(getattr(phys0, var).values))}
            print(f"  {var}: fit rms_relative {rec['rms_relative']:.4f}, "
                  f"pedestal top {fits[var]['ped_top']:.4g}")
            for s in scales:
                tag = f"{shot}_{var}_{s:.3f}"
                print(f"    {tag}")
                row = run_one(phys0, fit, var, s,
                              os.path.join(workdir, str(shot), f"{var}_{s:.3f}"),
                              radii, args.allow_missing_chease)
                row["ped_top_ratio"] = row["ped_top"] / fits[var]["ped_top"]
                rows.append(row)
                if "iterations" in row:
                    print(f"      iters {row['iterations']}  "
                          f"converged {row['converged']}  "
                          f"accepted {row['accepted']}  "
                          f"Ip err {row['ip_error_rel']:.4%}  "
                          f"{row['wall_s']:.0f}s")
        results[shot] = {"radii": list(radii),
                         "placeholder_radii": DISCHARGES[shot]["placeholder_radii"],
                         "fits": fits, "rows": rows}

    out = os.path.join(workdir, "reshape_convergence.json")
    with open(out, "w") as f:
        json.dump({"root": args.root, "scales": scales, "variables": variables,
                   "cheasebs": CHEASEBS, "results": results}, f,
                  indent=1, default=str)

    # ---- tables -------------------------------------------------------
    print("\n" + "=" * 78)
    print("iterations vs reshape  (blank = solve raised)")
    hdr = f"{'shot':>7} {'var':<4} " + " ".join(f"{s:>7.2f}" for s in scales)
    print(hdr); print("-" * len(hdr))
    for shot in shots:
        r = results.get(shot, {})
        for var in variables:
            cells = []
            for s in scales:
                row = next((x for x in r.get("rows", [])
                            if x["var"] == var and x["scale"] == s), None)
                cells.append(f"{row.get('iterations', ''):>7}" if row and
                             row.get("iterations") is not None else f"{'--':>7}")
            print(f"{shot:>7} {var:<4} " + " ".join(cells))
    print("  A smooth rise with |scale-1| is the solver behaving. An erratic row, "
          "or one direction much worse than the other, is not.")

    print("\nacceptance vs reshape  (. accepted, R rejected, X raised)")
    print(hdr); print("-" * len(hdr))
    for shot in shots:
        r = results.get(shot, {})
        for var in variables:
            cells = []
            for s in scales:
                row = next((x for x in r.get("rows", [])
                            if x["var"] == var and x["scale"] == s), None)
                if row is None or "error" in row:
                    cells.append(f"{'X':>7}")
                else:
                    cells.append(f"{'.' if row.get('accepted') else 'R':>7}")
            print(f"{shot:>7} {var:<4} " + " ".join(cells))
    print("  The reshape limit is where a row turns to R or X. Campaign bounds "
          "should sit inside it, with margin.")

    print("\nIp error vs reshape")
    print(hdr); print("-" * len(hdr))
    for shot in shots:
        r = results.get(shot, {})
        for var in variables:
            cells = []
            for s in scales:
                row = next((x for x in r.get("rows", [])
                            if x["var"] == var and x["scale"] == s), None)
                v = row.get("ip_error_rel") if row else None
                cells.append(f"{v * 100:>6.2f}%" if isinstance(v, (int, float))
                             else f"{'--':>7}")
            print(f"{shot:>7} {var:<4} " + " ".join(cells))
    print("  Under qspec=on CHEASE imposes q and Ip is an output, so this is the "
          "solver's response, not a target it is chasing. Compare +delta against "
          "-delta: a one-sided response is the open asymmetry.")

    print(f"\n[test] results  {out}")
    print(f"[test] per-run iteration_errors.png under {workdir}/<shot>/<var>_<scale>/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
